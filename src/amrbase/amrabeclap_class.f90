!> @file amrabeclap_class.f90
!> @brief Variable coefficient linear solver for AMR: (α·A - β·∇·(B∇))ϕ = f
!> @details Wraps AMReX's amrex_abeclaplacian operator with amrex_multigrid solver.
!>          A is cell-centered, B is face-centered. Default: α=0, β=1 (pure diffusion).
!>          Note: AMReX uses 2nd-order central differences for the Laplacian internally.
module amrabeclap_class
   use precision,          only: WP
   use amrgrid_class,      only: amrgrid
   use amrdata_class,      only: amrdata
   use amrex_amr_module,   only: amrex_multifab, amrex_geometry, amrex_boxarray, &
      amrex_distromap
   use amrex_linear_solver_module
   implicit none
   private

   ! Public constants (same as amrmg for consistency)
   integer, parameter, public :: abeclap_bc_interior   = 0
   integer, parameter, public :: abeclap_bc_dirichlet  = 101  ! AMREX_LO_DIRICHLET
   integer, parameter, public :: abeclap_bc_neumann    = 102  ! AMREX_LO_NEUMANN
   integer, parameter, public :: abeclap_bc_periodic   = 200  ! AMREX_LO_PERIODIC

   ! Bottom solver types
   integer, parameter, public :: abeclap_bottom_bicgstab = 1
   integer, parameter, public :: abeclap_bottom_cg       = 2
   integer, parameter, public :: abeclap_bottom_hypre    = 3
   integer, parameter, public :: abeclap_bottom_smoother = 4

   !> Variable coefficient AMR linear solver
   !> Solves: (α·A - β·∇·(B∇))ϕ = f
   type, public :: amrabeclap
      ! Solver configuration (public for direct access)
      real(WP) :: alpha = 0.0_WP    !< Scalar coefficient for A term
      real(WP) :: beta  = 1.0_WP    !< Scalar coefficient for ∇·(B∇) term
      real(WP) :: tol_rel = 1.0e-10_WP
      real(WP) :: tol_abs = 0.0_WP
      integer  :: verbose = 0
      integer  :: max_iter = 200
      integer  :: bottom_solver = abeclap_bottom_bicgstab
      integer  :: maxorder = 2      !< BC/interpolation order

      ! Output from solve
      real(WP) :: res   = 0.0_WP    !< Final residual after solve
      integer  :: niter = 0         !< Number of iterations

      ! Private internals
      type(amrgrid), pointer, private :: amr => null()
      integer, private :: lo_bc(3), hi_bc(3)
      logical, private :: initialized = .false.
   contains
      procedure :: initialize
      procedure :: solve
      procedure :: finalize
   end type amrabeclap

contains

   !> Initialize solver with grid and boundary conditions
   subroutine initialize(this, amr, lo_bc, hi_bc)
      class(amrabeclap), intent(inout) :: this
      type(amrgrid), target, intent(in) :: amr
      integer, intent(in) :: lo_bc(3), hi_bc(3)

      this%amr => amr
      this%lo_bc = lo_bc
      this%hi_bc = hi_bc
      this%initialized = .true.

   end subroutine initialize

   !> Solve (α·A - β·∇·(B∇))ϕ = rhs
   !> @param phi Solution field (in: initial guess with BC in ghosts, out: solution)
   !> @param rhs Right-hand side field
   !> @param acoef Optional cell-centered A coefficient (default: uniform 1)
   !> @param bcoef_x,bcoef_y,bcoef_z Optional face-centered B coefficients (default: uniform 1)
   subroutine solve(this, phi, rhs, acoef, bcoef_x, bcoef_y, bcoef_z)
      use amrex_abeclaplacian_module
      use amrex_multigrid_module
      use messager, only: die
      class(amrabeclap), intent(inout) :: this
      type(amrdata), intent(inout) :: phi
      type(amrdata), intent(in) :: rhs
      type(amrdata), intent(in), optional :: acoef
      type(amrdata), intent(in), optional :: bcoef_x, bcoef_y, bcoef_z

      type(amrex_abeclaplacian) :: linop
      type(amrex_multigrid) :: multigrid
      type(amrex_geometry), allocatable :: geom(:)
      type(amrex_boxarray), allocatable :: ba(:)
      type(amrex_distromap), allocatable :: dm(:)
      type(amrex_multifab), allocatable :: sol(:), rhsmf(:)
      type(amrex_multifab) :: bcoef(3)
      integer :: lev

      if (.not. this%initialized) call die('amrabeclap not initialized')

      ! Build AMReX wrapper arrays from amrgrid
      allocate(geom(0:this%amr%maxlvl))
      allocate(ba(0:this%amr%maxlvl))
      allocate(dm(0:this%amr%maxlvl))
      allocate(sol(0:this%amr%maxlvl))
      allocate(rhsmf(0:this%amr%maxlvl))

      do lev = 0, this%amr%maxlvl
         geom(lev) = this%amr%geom(lev)
         ba(lev) = this%amr%get_boxarray(lev)
         dm(lev) = this%amr%get_distromap(lev)
         sol(lev) = phi%mf(lev)
         rhsmf(lev) = rhs%mf(lev)
      end do

      ! Build the ABecLaplacian operator
      call amrex_abeclaplacian_build(linop, geom, ba, dm, &
         metric_term=.false., agglomeration=.true., consolidation=.true., &
         max_coarsening_level=30)

      ! Set domain boundary conditions
      call linop%set_domain_bc(this%lo_bc, this%hi_bc)

      ! Set maxorder for BC interpolation
      call linop%set_maxorder(this%maxorder)

      ! Set scalar coefficients
      call linop%set_scalars(this%alpha, this%beta)

      ! Set spatially-varying A coefficient if provided
      if (present(acoef)) then
         do lev = 0, this%amr%maxlvl
            call linop%set_acoeffs(lev, acoef%mf(lev))
         end do
      end if

      ! Set spatially-varying B coefficients if provided
      if (present(bcoef_x) .and. present(bcoef_y) .and. present(bcoef_z)) then
         do lev = 0, this%amr%maxlvl
            bcoef(1) = bcoef_x%mf(lev)
            bcoef(2) = bcoef_y%mf(lev)
            bcoef(3) = bcoef_z%mf(lev)
            call linop%set_bcoeffs(lev, bcoef)
         end do
      end if

      ! Set level BCs from solution ghost cells
      do lev = 0, this%amr%maxlvl
         call linop%set_level_bc(lev, sol(lev))
      end do

      ! Build multigrid solver
      call amrex_multigrid_build(multigrid, linop)
      call multigrid%set_verbose(this%verbose)
      call multigrid%set_max_iter(this%max_iter)
      call multigrid%set_bottom_solver(this%bottom_solver)

      ! Solve
      this%res = multigrid%solve(sol, rhsmf, this%tol_rel, this%tol_abs)

      ! Clean up
      call amrex_multigrid_destroy(multigrid)
      call amrex_abeclaplacian_destroy(linop)

      this%niter = 0  ! TODO: Get actual count from MLMG
   end subroutine solve

   !> Finalize and release resources
   subroutine finalize(this)
      class(amrabeclap), intent(inout) :: this
      this%amr => null()
      this%initialized = .false.
   end subroutine finalize

end module amrabeclap_class

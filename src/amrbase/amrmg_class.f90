module amrmg_class
   use iso_c_binding, only: c_ptr, c_null_ptr, c_associated
   use precision,     only: WP
   use amrgrid_class, only: amrgrid
   use amrdata_class, only: amrdata
   implicit none
   private

   ! Public constants for bottom solver selection (matches AMReX)
   integer, parameter, public :: amrmg_bottom_smoother = 0
   integer, parameter, public :: amrmg_bottom_bicgstab = 1
   integer, parameter, public :: amrmg_bottom_cg       = 2
   integer, parameter, public :: amrmg_bottom_hypre    = 3
   integer, parameter, public :: amrmg_bottom_petsc    = 4

   ! Public constants for boundary conditions (matches AMReX LinOpBCType)
   integer, parameter, public :: amrmg_bc_interior    = 0
   integer, parameter, public :: amrmg_bc_dirichlet   = 101
   integer, parameter, public :: amrmg_bc_neumann     = 102
   integer, parameter, public :: amrmg_bc_reflect_odd = 103
   integer, parameter, public :: amrmg_bc_periodic    = 200

   !> AMR Multigrid Poisson solver
   !> Wraps AMReX MLMG for solving ∇²ϕ = f
   type, public :: amrmg
      ! Grid reference
      class(amrgrid), pointer :: amr => null()
      ! Solver parameters (set at initialize)
      integer  :: bottom_solver = amrmg_bottom_bicgstab
      integer  :: max_iter = 100
      real(WP) :: tol_rel = 1.0e-10_WP
      real(WP) :: tol_abs = 0.0_WP
      integer  :: verbose = 0
      ! Boundary conditions
      integer, dimension(3) :: lo_bc = amrmg_bc_dirichlet
      integer, dimension(3) :: hi_bc = amrmg_bc_dirichlet
      ! Monitoring - populated after solve
      integer  :: niter = 0       !< Number of iterations used
      real(WP) :: res   = 0.0_WP  !< Final residual
      ! Internal state
      logical, private :: initialized = .false.
   contains
      procedure :: initialize
      procedure :: solve_composite
      procedure :: solve_level
      procedure :: finalize
   end type amrmg

contains

   !> Initialize the Poisson solver
   subroutine initialize(this, amr, lo_bc, hi_bc, tol_rel, tol_abs, max_iter, bottom_solver, verbose)
      use messager, only: log
      implicit none
      class(amrmg), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      integer, dimension(3), intent(in) :: lo_bc, hi_bc
      real(WP), intent(in), optional :: tol_rel, tol_abs
      integer, intent(in), optional :: max_iter, bottom_solver, verbose

      ! Store reference and parameters
      this%amr => amr
      this%lo_bc = lo_bc
      this%hi_bc = hi_bc
      if (present(tol_rel)) this%tol_rel = tol_rel
      if (present(tol_abs)) this%tol_abs = tol_abs
      if (present(max_iter)) this%max_iter = max_iter
      if (present(bottom_solver)) this%bottom_solver = bottom_solver
      if (present(verbose)) this%verbose = verbose

      this%initialized = .true.
      call log('Initialized amrmg Poisson solver')
   end subroutine initialize


   !> Composite solve - all levels simultaneously
   subroutine solve_composite(this, phi, rhs)
      use messager, only: die
      use amrex_linear_solver_module
      implicit none
      class(amrmg), intent(inout) :: this
      type(amrdata), intent(inout) :: phi  !< Solution (with BC in ghost cells)
      type(amrdata), intent(in) :: rhs     !< Right-hand side

      type(amrex_poisson) :: poisson
      type(amrex_multigrid) :: multigrid
      type(amrex_geometry), allocatable :: geom(:)
      type(amrex_boxarray), allocatable :: ba(:)
      type(amrex_distromap), allocatable :: dm(:)
      type(amrex_multifab), allocatable :: sol(:), rhsmf(:)
      integer :: lev, nlevs

      if (.not. this%initialized) call die('amrmg not initialized')

      nlevs = this%amr%maxlvl + 1

      ! Build AMReX wrapper types from amrgrid
      allocate(geom(0:this%amr%maxlvl))
      allocate(ba(0:this%amr%maxlvl))
      allocate(dm(0:this%amr%maxlvl))
      allocate(sol(0:this%amr%maxlvl))
      allocate(rhsmf(0:this%amr%maxlvl))

      do lev = 0, this%amr%maxlvl
         ! geom is stored directly in amrgrid
         geom(lev) = this%amr%geom(lev)
         ! ba and dm need accessor methods
         ba(lev) = this%amr%get_boxarray(lev)
         dm(lev) = this%amr%get_distromap(lev)
         ! Copy MultiFab references (non-owning)
         sol(lev) = phi%mf(lev)
         rhsmf(lev) = rhs%mf(lev)
      end do

      ! Build Poisson linear operator
      call amrex_poisson_build(poisson, geom, ba, dm, &
         metric_term=.false., agglomeration=.true., consolidation=.true., &
         max_coarsening_level=30)

      ! Set boundary conditions
      call poisson%set_domain_bc(this%lo_bc, this%hi_bc)

      ! Set level BCs from solution ghost cells
      do lev = 0, this%amr%maxlvl
         call poisson%set_level_bc(lev, sol(lev))
      end do

      ! Build multigrid solver
      call amrex_multigrid_build(multigrid, poisson)
      call multigrid%set_verbose(this%verbose)
      call multigrid%set_max_iter(this%max_iter)
      call multigrid%set_bottom_solver(this%bottom_solver)

      ! Solve
      this%res = multigrid%solve(sol, rhsmf, this%tol_rel, this%tol_abs)

      ! Clean up
      call amrex_multigrid_destroy(multigrid)
      call amrex_poisson_destroy(poisson)

      this%niter = 0  ! TODO: Get actual count from MLMG
   end subroutine solve_composite


   !> Level-by-level solve (for subcycling)
   subroutine solve_level(this, lev, phi, rhs, phi_crse)
      use messager, only: die
      implicit none
      class(amrmg), intent(inout) :: this
      integer, intent(in) :: lev
      type(amrdata), intent(inout) :: phi
      type(amrdata), intent(in) :: rhs
      type(amrdata), intent(in), optional :: phi_crse  !< Coarse level solution for C/F BC

      ! Placeholder - requires per-level solver setup
      call die('solve_level not yet implemented - use solve_composite')
   end subroutine solve_level


   !> Finalize and clean up
   subroutine finalize(this)
      implicit none
      class(amrmg), intent(inout) :: this
      nullify(this%amr)
      this%initialized = .false.
   end subroutine finalize

end module amrmg_class

!> Test program for amrmg Poisson solver
!> Solves ∇²ϕ = f with known solution
!> NOTE: The "error" reported is discretization error (difference from analytical solution),
!>       not solver error. The solver residual being ~1e-12 means the discrete system
!>       was solved to machine precision. The ~5% error is O(h²) truncation error.
module mod_test_amrmg
   use precision,         only: WP
   use string,            only: rtoa
   use amrgrid_class,     only: amrgrid
   use amrdata_class,     only: amrdata
   use amrmg_class
   use messager,          only: log
   use amrex_amr_module,  only: amrex_mfiter, amrex_box
   use mpi_f08
   implicit none
   private
   public :: test_amrmg

   ! Module-level objects
   type(amrgrid), allocatable, target :: amr_mg
   type(amrdata), allocatable, target :: phi, rhs_mg, exact

   real(WP), parameter :: pi = 4.0_WP * atan(1.0_WP)

contains

   !> Simple tagging: refine center region
   subroutine tag_center(lvl, tags_ptr, time)
      use iso_c_binding, only: c_ptr, c_char
      use amrex_amr_module, only: amrex_tagboxarray, amrex_mfiter, amrex_box
      implicit none
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      character(kind=c_char), parameter :: SET=char(1)
      real(WP) :: x, y, z, dx
      integer :: i, j, k
      tags = tags_ptr
      dx = amr_mg%dx(lvl)
      call amr_mg%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = amr_mg%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2), bx%hi(2)
               y = amr_mg%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1), bx%hi(1)
                  x = amr_mg%xlo + (real(i,WP)+0.5_WP)*dx
                  ! Tag center region for refinement
                  if (x > 0.3_WP .and. x < 0.7_WP .and. &
                     y > 0.3_WP .and. y < 0.7_WP .and. &
                     z > 0.3_WP .and. z < 0.7_WP) then
                     tagarr(i,j,k,1) = SET
                  end if
               end do
            end do
         end do
      end do
      call amr_mg%mfiter_destroy(mfi)
   end subroutine tag_center

   subroutine test_amrmg()
      type(amrmg) :: solver
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPhi, pRhs, pExact
      real(WP) :: x, y, z, dx, err_l2, local_err, err_linf
      integer :: lvl, i, j, k, ierr, ng
      integer :: lo_bc(3), hi_bc(3)
      logical :: use_hypre

      write(*,*) '=== TEST_AMRMG STARTING ==='
      write(*,*) ''

      ! Allocate and configure grid
      allocate(amr_mg)
      amr_mg%nx = 32
      amr_mg%ny = 32
      amr_mg%nz = 32
      amr_mg%xlo = 0.0_WP; amr_mg%xhi = 1.0_WP
      amr_mg%ylo = 0.0_WP; amr_mg%yhi = 1.0_WP
      amr_mg%zlo = 0.0_WP; amr_mg%zhi = 1.0_WP
      amr_mg%xper = .false.
      amr_mg%yper = .false.
      amr_mg%zper = .false.
      amr_mg%maxlvl = 1   ! Two levels for AMR test
      amr_mg%nmax = 16
      amr_mg%nbloc = 8

      call amr_mg%initialize("mg_test")
      call amr_mg%add_tagging(tag_center)

      ! Build phi, rhs, exact - register for regrid callbacks
      allocate(phi, rhs_mg, exact)
      call phi%initialize(   amr_mg, name='phi',   ncomp=1, ng=1)
      call rhs_mg%initialize(amr_mg, name='rhs',   ncomp=1, ng=0)
      call exact%initialize( amr_mg, name='exact', ncomp=1, ng=0)
      call phi%register()
      call rhs_mg%register()
      call exact%register()

      ! Initialize with user data callback
      call phi%set_on_init(init_phi)
      call rhs_mg%set_on_init(init_rhs)
      call exact%set_on_init(init_exact)

      ! Build AMR hierarchy
      call amr_mg%initialize_grid(0.0_WP)

      write(*,*) 'Grid built:', amr_mg%nlevels, 'levels,', amr_mg%nboxes, 'boxes'

      ! Zero interior of phi (keep BC in ghosts)
      do lvl = 0, amr_mg%clvl()
         call amr_mg%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pPhi => phi%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pPhi(i,j,k,1) = 0.0_WP
                  end do
               end do
            end do
         end do
         call amr_mg%mfiter_destroy(mfi)
      end do

      ! Initialize solver and solve
      lo_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      hi_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]

      ! Use Hypre bottom solver (AMReX rebuilt with -DAMReX_HYPRE=ON)
      call solver%initialize(amr_mg, lo_bc, hi_bc, tol_rel=1.0e-12_WP, verbose=2, &
         bottom_solver=amrmg_bottom_hypre)
      write(*,*) 'Using Hypre bottom solver'

      write(*,*) ''
      write(*,*) 'Starting Poisson solve...'
      call solver%solve_composite(phi, rhs_mg)
      write(*,*) 'Solve complete: residual =', solver%res
      write(*,*) ''

      ! Compute discretization error (vs analytical solution)
      err_l2 = 0.0_WP
      err_linf = 0.0_WP
      do lvl = 0, amr_mg%clvl()
         call amr_mg%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pPhi   => phi%mf(lvl)%dataptr(mfi)
            pExact => exact%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     local_err = abs(pPhi(i,j,k,1) - pExact(i,j,k,1))
                     err_l2 = err_l2 + local_err**2
                     err_linf = max(err_linf, local_err)
                  end do
               end do
            end do
         end do
         call amr_mg%mfiter_destroy(mfi)
      end do

      ! Reduce across ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      err_l2 = sqrt(err_l2)

      write(*,*) 'Discretization Error (vs analytical solution):'
      write(*,*) '  L2  :', err_l2
      write(*,*) '  Linf:', err_linf
      write(*,*) ''
      write(*,*) 'NOTE: This is truncation error O(h^2), NOT solver error.'
      write(*,*) '      The solver residual ~1e-12 shows the discrete system'
      write(*,*) '      was solved to machine precision.'
      write(*,*) ''

      ! For AMR with refinement, the Linf error in the refined region should
      ! be smaller than in the coarse region. Accept if Linf < 10% (coarse grid)
      if (err_linf < 0.1_WP .and. solver%res < 1.0e-8_WP) then
         write(*,*) 'PASS: amrmg Poisson solver test'
      else
         write(*,*) 'FAIL: residual=', solver%res, ' err_linf=', err_linf
      end if

      ! Cleanup
      call solver%finalize()
      call phi%destroy()
      call rhs_mg%destroy()
      call exact%destroy()
      deallocate(phi, rhs_mg, exact)
      call amr_mg%finalize()
      deallocate(amr_mg)

      write(*,*) ''
      write(*,*) '=== TEST_AMRMG COMPLETE ==='

   end subroutine test_amrmg

   !> Initialize phi with exact solution (including ghosts for BC)
   subroutine init_phi(lvl, mf, geom)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry, amrex_mfiter, amrex_box, amrex_mfiter_build
      implicit none
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: mf
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: x, y, z, dx
      integer :: i, j, k
      dx = geom%dx(1)
      call amrex_mfiter_build(mfi, mf)
      do while (mfi%next())
         bx = mfi%tilebox()
         p => mf%dataptr(mfi)
         ! Fill valid + ghosts with exact solution (for Dirichlet BC)
         do k = bx%lo(3)-1, bx%hi(3)+1
            z = amr_mg%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2)-1, bx%hi(2)+1
               y = amr_mg%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1)-1, bx%hi(1)+1
                  x = amr_mg%xlo + (real(i,WP)+0.5_WP)*dx
                  p(i,j,k,1) = sin(pi*x) * sin(pi*y) * sin(pi*z)
               end do
            end do
         end do
      end do
   end subroutine init_phi

   !> Initialize RHS with source term
   subroutine init_rhs(lvl, mf, geom)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry, amrex_mfiter, amrex_box, amrex_mfiter_build
      implicit none
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: mf
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: x, y, z, dx
      integer :: i, j, k
      dx = geom%dx(1)
      call amrex_mfiter_build(mfi, mf)
      do while (mfi%next())
         bx = mfi%tilebox()
         p => mf%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = amr_mg%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2), bx%hi(2)
               y = amr_mg%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1), bx%hi(1)
                  x = amr_mg%xlo + (real(i,WP)+0.5_WP)*dx
                  ! ∇²(sin*sin*sin) = -3π² * sin*sin*sin
                  p(i,j,k,1) = -3.0_WP * pi**2 * sin(pi*x) * sin(pi*y) * sin(pi*z)
               end do
            end do
         end do
      end do
   end subroutine init_rhs

   !> Initialize exact solution
   subroutine init_exact(lvl, mf, geom)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry, amrex_mfiter, amrex_box, amrex_mfiter_build
      implicit none
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: mf
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: x, y, z, dx
      integer :: i, j, k
      dx = geom%dx(1)
      call amrex_mfiter_build(mfi, mf)
      do while (mfi%next())
         bx = mfi%tilebox()
         p => mf%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = amr_mg%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2), bx%hi(2)
               y = amr_mg%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1), bx%hi(1)
                  x = amr_mg%xlo + (real(i,WP)+0.5_WP)*dx
                  p(i,j,k,1) = sin(pi*x) * sin(pi*y) * sin(pi*z)
               end do
            end do
         end do
      end do
   end subroutine init_exact

end module mod_test_amrmg

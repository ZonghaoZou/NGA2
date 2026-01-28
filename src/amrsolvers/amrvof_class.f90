!> AMR VOF Solver class
!> Provides Volume-of-Fluid advection for two-phase flow with AMReX
!> IRL-free implementation using native cutting geometry
module amrvof_class
   use iso_c_binding,    only: c_ptr, c_null_ptr, c_loc, c_f_pointer, c_char
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrvof_geometry,  only: VFlo, VFhi
   use surfmesh_class,   only: surfmesh
   use amrex_amr_module, only: amrex_multifab, amrex_mfiter, amrex_box, &
   &                           amrex_boxarray, amrex_distromap, amrex_geometry
   implicit none
   private

   ! Expose type and dispatchers
   public :: amrvof, VFlo, VFhi
   public :: amrvof_on_init, amrvof_on_coarse, amrvof_on_remake
   public :: amrvof_on_clear, amrvof_tagging, amrvof_postregrid

   ! PLIC boundary condition types
   integer, parameter, public :: BC_LIQ     = 1  !< All liquid in ghost
   integer, parameter, public :: BC_GAS     = 2  !< All gas in ghost
   integer, parameter, public :: BC_REFLECT = 3  !< Symmetry (mirror across boundary)
   integer, parameter, public :: BC_USER    = 4  !< User-defined callback

   !> AMR VOF solver type
   type, extends(amrsolver) :: amrvof
      ! User-configurable callbacks
      procedure(vof_init_iface), pointer, nopass :: user_init => null()
      procedure(vof_tagging_iface), pointer, nopass :: user_tagging => null()
      procedure(vof_bc_iface), pointer, nopass :: user_vof_bc => null()

      ! PLIC boundary conditions (per face, only used if direction is non-periodic)
      integer :: bc_xlo = BC_REFLECT
      integer :: bc_xhi = BC_REFLECT
      integer :: bc_ylo = BC_REFLECT
      integer :: bc_yhi = BC_REFLECT
      integer :: bc_zlo = BC_REFLECT
      integer :: bc_zhi = BC_REFLECT

      ! VOF data (solver owns these - 4 MultiFabs as per plan)
      type(amrdata) :: VF           !< Volume fraction (cell-centered)
      type(amrdata) :: Cliq         !< Liquid barycenter (3 components)
      type(amrdata) :: Cgas         !< Gas barycenter (3 components)
      type(amrdata) :: PLIC         !< PLIC plane (4 components: nx, ny, nz, d)

      ! Old data for time stepping
      type(amrdata) :: VFold
      type(amrdata) :: Cliqold
      type(amrdata) :: Cgasold
      type(amrdata) :: PLICold

      ! Monitoring quantities
      real(WP) :: VFmin = 0.0_WP    !< Minimum VF
      real(WP) :: VFmax = 0.0_WP    !< Maximum VF
      real(WP) :: VFint = 0.0_WP    !< Integral of VF (liquid volume)
      
      ! Surface mesh for visualization (finest level polygons)
      type(surfmesh) :: smesh

   contains
      procedure :: initialize
      procedure :: finalize
      ! Override internal type-bound callbacks from amrsolver
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      ! Deferred from amrsolver base class
      procedure :: get_info
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
      ! VOF-specific procedures
      procedure :: build_plic             !< Reconstruct PLIC from VF and barycenters
      procedure :: advance_vof            !< Advect VF using staggered or collocated velocity
      procedure :: fill_moments_lvl       !< Fill VF/Cliq/Cgas ghosts at level (sync + BC)
      procedure :: sync_moments_lvl       !< Sync VF/Cliq/Cgas ghosts at level + fix periodic barycenters
      procedure :: sync_moments           !< Sync VF/Cliq/Cgas ghosts all levels + fix periodic barycenters
      procedure :: sync_plic_lvl          !< Sync PLIC ghosts at level + fix periodic plane distance
      procedure :: sync_plic              !< Sync PLIC ghosts all levels + fix periodic plane distance
      procedure :: fill_plic_lvl          !< Fill PLIC ghosts at level (sync + BC)
      procedure :: average_down_moments   !< Average down VF/Cliq/Cgas to coarse levels
      procedure :: reset_moments          !< Recompute VF/barycenters from PLIC
      procedure :: get_cfl                !< Compute advective CFL at finest level
      procedure :: print => amrvof_print  !< Print solver info
   end type amrvof

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine vof_init_iface(solver, lvl, time, ba, dm)
         import :: amrvof, WP, amrex_boxarray, amrex_distromap
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine vof_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine vof_tagging_iface(solver, lvl, tags, time)
         import :: amrvof, c_ptr, WP
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags
         real(WP), intent(in) :: time
      end subroutine vof_tagging_iface
   end interface

   !> Abstract interface for user-defined VOF boundary condition
   !> User must set VF, Cliq, Cgas, and PLIC consistently in ghost cells
   abstract interface
      subroutine vof_bc_iface(solver, lvl, time, dir, side, pVF, pCliq, pCgas, pPLIC, ilo, ihi, jlo, jhi, klo, khi)
         import :: amrvof, WP
         class(amrvof), intent(inout) :: solver
         integer, intent(in) :: lvl       !< Level
         real(WP), intent(in) :: time     !< Current simulation time
         integer, intent(in) :: dir       !< Direction (1=x, 2=y, 3=z)
         integer, intent(in) :: side      !< Side (-1=lo, +1=hi)
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pVF, pCliq, pCgas, pPLIC
         integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi  !< Ghost region bounds
      end subroutine vof_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrvof type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrvof_on_init(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(lvl, time, ba, dm)
      if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
   end subroutine amrvof_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrvof_on_coarse(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_coarse(lvl, time, ba, dm)
   end subroutine amrvof_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrvof_on_remake(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_remake(lvl, time, ba, dm)
   end subroutine amrvof_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrvof_on_clear(ctx, lvl)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_clear(lvl)
   end subroutine amrvof_on_clear

   !> Dispatch tagging: tag cells near interface (3x3x3 stencil check)
   subroutine amrvof_tagging(ctx, lvl, tags, time)
      use amrex_amr_module, only: amrex_tagboxarray, amrex_mfiter, amrex_box
      use amrgrid_class, only: SETtag
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrvof), pointer :: this
      type(amrex_tagboxarray) :: tba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      logical :: near_interface
      integer :: i, j, k, ii, jj, kk
      real(WP) :: vf_center, vf_neighbor
      
      call c_f_pointer(ctx, this)
      
      ! Tag cells at interface using 3x3x3 stencil
      tba = tags
      call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tba%dataPtr(mfi)
         pVF => this%VF%mf(lvl)%dataptr(mfi)
         
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  vf_center = pVF(i,j,k,1)
                  near_interface = .false.
                  
                  ! Check if cell itself is interfacial
                  if (vf_center .gt. VFlo .and. vf_center .lt. VFhi) then
                     near_interface = .true.
                  else
                     ! Pure cell - check 3x3x3 neighborhood for different phase
                     stencil: do kk = k-1, k+1
                        do jj = j-1, j+1
                           do ii = i-1, i+1
                              vf_neighbor = pVF(ii,jj,kk,1)
                              ! If center is pure gas and neighbor has liquid, or vice versa
                              if (vf_center .lt. VFlo .and. vf_neighbor .gt. VFlo) then
                                 near_interface = .true.
                                 exit stencil
                              end if
                              if (vf_center .gt. VFhi .and. vf_neighbor .lt. VFhi) then
                                 near_interface = .true.
                                 exit stencil
                              end if
                           end do
                        end do
                     end do stencil
                  end if
                  
                  if (near_interface) tagarr(i,j,k,1) = SETtag
               end do
            end do
         end do
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Call user tagging if provided
      if (associated(this%user_tagging)) call this%user_tagging(this, lvl, tags, time)
   end subroutine amrvof_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrvof_postregrid(ctx, lbase, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrvof), pointer :: this
      call c_f_pointer(ctx, this)
      call this%post_regrid(lbase, time)
   end subroutine amrvof_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the VOF solver
   subroutine initialize(this, amr, name)
      class(amrvof), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'UNNAMED_VOF'
      end if

      ! Store amrgrid pointer
      this%amr => amr

      ! Initialize VOF data (cell-centered)
      call this%VF%initialize(amr, name='VF', ncomp=1, ng=1)
      call this%Cliq%initialize(amr, name='Cliq', ncomp=3, ng=1)
      call this%Cgas%initialize(amr, name='Cgas', ncomp=3, ng=1)
      call this%PLIC%initialize(amr, name='PLIC', ncomp=4, ng=2)

      ! Initialize old data
      call this%VFold%initialize(amr, name='VFold', ncomp=1, ng=1)
      call this%Cliqold%initialize(amr, name='Cliqold', ncomp=3, ng=1)
      call this%Cgasold%initialize(amr, name='Cgasold', ncomp=3, ng=1)
      call this%PLICold%initialize(amr, name='PLICold', ncomp=4, ng=2)
      
      ! Initialize surface mesh for visualization
      this%smesh%name = trim(this%name)//'_plic'

      ! Set parent pointers for callback context access
      this%VF%parent => this
      this%Cliq%parent => this
      this%Cgas%parent => this
      this%PLIC%parent => this
      this%VFold%parent => this
      this%Cliqold%parent => this
      this%Cgasold%parent => this
      this%PLICold%parent => this

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      ! NOTE: Standard amrdata interpolation is fine IF the tagger ensures the
      !       interface never escapes the finest level between regrids. Tag with
      !       sufficient buffer: interface cells + CFL*dt*regrid_interval.
      select type (this)
       type is (amrvof)
         call this%amr%add_on_init   (amrvof_on_init,    c_loc(this))
         call this%amr%add_on_coarse (amrvof_on_coarse,  c_loc(this))
         call this%amr%add_on_remake (amrvof_on_remake,  c_loc(this))
         call this%amr%add_on_clear  (amrvof_on_clear,   c_loc(this))
         call this%amr%add_tagging   (amrvof_tagging,    c_loc(this))
         call this%amr%add_postregrid(amrvof_postregrid, c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts
      call this%VF%reset_level(lvl, ba, dm)
      call this%Cliq%reset_level(lvl, ba, dm)
      call this%Cgas%reset_level(lvl, ba, dm)
      call this%PLIC%reset_level(lvl, ba, dm)
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Set to zero
      call this%VF%setval(val=0.0_WP, lvl=lvl)
      call this%Cliq%setval(val=0.0_WP, lvl=lvl)
      call this%Cgas%setval(val=0.0_WP, lvl=lvl)
      call this%PLIC%setval(val=0.0_WP, lvl=lvl)
      call this%VFold%setval(val=0.0_WP, lvl=lvl)
      call this%Cliqold%setval(val=0.0_WP, lvl=lvl)
      call this%Cgasold%setval(val=0.0_WP, lvl=lvl)
      call this%PLICold%setval(val=0.0_WP, lvl=lvl)
      ! Call user init to set VF
      if (associated(this%user_init)) call this%user_init(this, lvl, time, ba, dm)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse
   subroutine on_coarse(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Allocate and interpolate moments from coarse
      call this%VF%on_coarse(this%VF, lvl, time, ba, dm)
      call this%Cliq%on_coarse(this%Cliq, lvl, time, ba, dm)
      call this%Cgas%on_coarse(this%Cgas, lvl, time, ba, dm)
      ! PLIC: just allocate, then set to trivial planes (new cells are pure, away from interface)
      call this%PLIC%reset_level(lvl, ba, dm)
      call set_trivial_plic()
      ! Old data just needs geometry
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Fill moment and PLIC ghosts
      call this%fill_moments_lvl(lvl, time)
      call this%fill_plic_lvl(lvl, time)
   contains
      !> Set PLIC to trivial planes based on VF
      subroutine set_trivial_plic()
         use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         integer :: i, j, k
         call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     pPLIC(i,j,k,1:3) = 0.0_WP
                     pPLIC(i,j,k,4) = sign(1.0e10_WP, pVF(i,j,k,1) - 0.5_WP)
                  end do
               end do
            end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine set_trivial_plic
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid
   subroutine on_remake(this, lvl, time, ba, dm)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Remake all data (copies existing + fills from coarse for new areas)
      call this%VF%on_remake(this%VF, lvl, time, ba, dm)
      call this%Cliq%on_remake(this%Cliq, lvl, time, ba, dm)
      call this%Cgas%on_remake(this%Cgas, lvl, time, ba, dm)
      call this%PLIC%on_remake(this%PLIC, lvl, time, ba, dm)
      ! Fix PLIC for pure cells (new cells from coarse have interpolated PLIC which is wrong)
      call fix_pure_plic()
      ! Old data just needs geometry
      call this%VFold%reset_level(lvl, ba, dm)
      call this%Cliqold%reset_level(lvl, ba, dm)
      call this%Cgasold%reset_level(lvl, ba, dm)
      call this%PLICold%reset_level(lvl, ba, dm)
      ! Fill moment and PLIC ghosts
      call this%fill_moments_lvl(lvl, time)
      call this%fill_plic_lvl(lvl, time)
   contains
      !> Set PLIC to trivial for pure cells
      subroutine fix_pure_plic()
         use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         integer :: i, j, k
         real(WP) :: vf
         call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     vf = pVF(i,j,k,1)
                     if (vf.lt.VFlo .or. vf.gt.VFhi) then
                        pPLIC(i,j,k,1:3) = 0.0_WP
                        pPLIC(i,j,k,4) = sign(1.0e10_WP, vf - 0.5_WP)
                     end if
                  end do
               end do
            end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine fix_pure_plic
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this, lvl)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%VF%clear_level(lvl)
      call this%Cliq%clear_level(lvl)
      call this%Cgas%clear_level(lvl)
      call this%PLIC%clear_level(lvl)
      call this%VFold%clear_level(lvl)
      call this%Cliqold%clear_level(lvl)
      call this%Cgasold%clear_level(lvl)
      call this%PLICold%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this, lbase, time)
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      call this%average_down_moments(lbase)
   end subroutine post_regrid

   !> Average down VF/Cliq/Cgas from finest to lbase, then sync ghost cells
   subroutine average_down_moments(this, lbase)
      use amrex_interface, only: amrmfab_average_down_cell
      class(amrvof), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb = 0; if (present(lbase)) lb = lbase
      ! Average valid cells from fine to coarse
      do lvl = this%amr%clvl()-1, lb, -1
         call amrmfab_average_down_cell(fmf=this%VF%mf(lvl+1)  , cmf=this%VF%mf(lvl)  , rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_cell(fmf=this%Cliq%mf(lvl+1), cmf=this%Cliq%mf(lvl), rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_cell(fmf=this%Cgas%mf(lvl+1), cmf=this%Cgas%mf(lvl), rr=this%amr%rref(lvl), cgeom=this%amr%geom(lvl))
      end do
      ! Sync ghost cells on all levels + fix periodic barycenters
      call this%sync_moments()
   end subroutine average_down_moments

   !> Sync VF/Cliq/Cgas ghosts on all levels + fix periodic barycenters
   subroutine sync_moments(this)
      class(amrvof), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_moments_lvl(lvl)
      end do
   end subroutine sync_moments

   !> Sync VF/Cliq/Cgas ghosts at level + fix periodic barycenters
   subroutine sync_moments_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box, amrex_geometry
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pCliq, pCgas
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL
      
      ! Sync ghosts (periodic + MPI exchange)
      call this%VF%sync_lvl(lvl)
      call this%Cliq%sync_lvl(lvl)
      call this%Cgas%sync_lvl(lvl)
      
      ! Fix barycenter positions in periodic ghost cells
      geom = this%amr%geom(lvl)
      dlo = geom%domain%lo
      dhi = geom%domain%hi
      xL = this%amr%xhi - this%amr%xlo
      yL = this%amr%yhi - this%amr%ylo
      zL = this%amr%zhi - this%amr%zlo
      
      call amrex_mfiter_build(mfi, this%Cliq%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         ilo = lbound(pCliq,1); ihi = ubound(pCliq,1)
         jlo = lbound(pCliq,2); jhi = ubound(pCliq,2)
         klo = lbound(pCliq,3); khi = ubound(pCliq,3)
         ! X-periodic
         if (this%amr%xper) then
            if (ilo .lt. dlo(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = ilo, dlo(1)-1
                  pCliq(ig,jg,kg,1) = pCliq(ig,jg,kg,1) - xL
                  pCgas(ig,jg,kg,1) = pCgas(ig,jg,kg,1) - xL
               end do; end do; end do
            end if
            if (ihi .gt. dhi(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = dhi(1)+1, ihi
                  pCliq(ig,jg,kg,1) = pCliq(ig,jg,kg,1) + xL
                  pCgas(ig,jg,kg,1) = pCgas(ig,jg,kg,1) + xL
               end do; end do; end do
            end if
         end if
         ! Y-periodic
         if (this%amr%yper) then
            if (jlo .lt. dlo(2)) then
               do kg = klo, khi; do jg = jlo, dlo(2)-1; do ig = ilo, ihi
                  pCliq(ig,jg,kg,2) = pCliq(ig,jg,kg,2) - yL
                  pCgas(ig,jg,kg,2) = pCgas(ig,jg,kg,2) - yL
               end do; end do; end do
            end if
            if (jhi .gt. dhi(2)) then
               do kg = klo, khi; do jg = dhi(2)+1, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,2) = pCliq(ig,jg,kg,2) + yL
                  pCgas(ig,jg,kg,2) = pCgas(ig,jg,kg,2) + yL
               end do; end do; end do
            end if
         end if
         ! Z-periodic
         if (this%amr%zper) then
            if (klo .lt. dlo(3)) then
               do kg = klo, dlo(3)-1; do jg = jlo, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,3) = pCliq(ig,jg,kg,3) - zL
                  pCgas(ig,jg,kg,3) = pCgas(ig,jg,kg,3) - zL
               end do; end do; end do
            end if
            if (khi .gt. dhi(3)) then
               do kg = dhi(3)+1, khi; do jg = jlo, jhi; do ig = ilo, ihi
                  pCliq(ig,jg,kg,3) = pCliq(ig,jg,kg,3) + zL
                  pCgas(ig,jg,kg,3) = pCgas(ig,jg,kg,3) + zL
               end do; end do; end do
            end if
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine sync_moments_lvl

   !> Sync PLIC ghosts at level + fix periodic plane distance
   !> The plane distance d must be corrected in periodic ghost cells:
   !> d ← d ± n·L where L is domain length and n is normal component
   subroutine sync_plic_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL
      
      ! Sync ghosts (periodic + MPI exchange)
      call this%PLIC%sync_lvl(lvl)
      
      ! Get geometry and domain bounds
      geom = this%amr%geom(lvl)
      dlo = geom%domain%lo
      dhi = geom%domain%hi
      
      xL = this%amr%xhi - this%amr%xlo
      yL = this%amr%yhi - this%amr%ylo
      zL = this%amr%zhi - this%amr%zlo
      
      ! Fix periodic plane distance
      call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pP => this%PLIC%mf(lvl)%dataptr(mfi)
         ilo = lbound(pP,1); ihi = ubound(pP,1)
         jlo = lbound(pP,2); jhi = ubound(pP,2)
         klo = lbound(pP,3); khi = ubound(pP,3)
         
         ! X-periodic
         if (this%amr%xper) then
            if (ilo .lt. dlo(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = ilo, dlo(1)-1
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,1)*xL
               end do; end do; end do
            end if
            if (ihi .gt. dhi(1)) then
               do kg = klo, khi; do jg = jlo, jhi; do ig = dhi(1)+1, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,1)*xL
               end do; end do; end do
            end if
         end if
         
         ! Y-periodic
         if (this%amr%yper) then
            if (jlo .lt. dlo(2)) then
               do kg = klo, khi; do jg = jlo, dlo(2)-1; do ig = ilo, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,2)*yL
               end do; end do; end do
            end if
            if (jhi .gt. dhi(2)) then
               do kg = klo, khi; do jg = dhi(2)+1, jhi; do ig = ilo, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,2)*yL
               end do; end do; end do
            end if
         end if
         
         ! Z-periodic
         if (this%amr%zper) then
            if (klo .lt. dlo(3)) then
               do kg = klo, dlo(3)-1; do jg = jlo, jhi; do ig = ilo, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) - pP(ig,jg,kg,3)*zL
               end do; end do; end do
            end if
            if (khi .gt. dhi(3)) then
               do kg = dhi(3)+1, khi; do jg = jlo, jhi; do ig = ilo, ihi
                  pP(ig,jg,kg,4) = pP(ig,jg,kg,4) + pP(ig,jg,kg,3)*zL
               end do; end do; end do
            end if
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine sync_plic_lvl

   !> Sync PLIC ghosts on all levels + fix periodic plane distance
   subroutine sync_plic(this)
      class(amrvof), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_plic_lvl(lvl)
      end do
   end subroutine sync_plic

   !> Fill PLIC ghosts at a level (sync + physical BC)
   !> Handles BC_LIQ (trivial d=+∞), BC_GAS (trivial d=-∞), 
   !> BC_REFLECT (mirror + flip normal), BC_USER (user callback)
   subroutine fill_plic_lvl(this, lvl, time)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      
      ! Sync PLIC ghosts + fix periodic plane distance
      call this%sync_plic_lvl(lvl)
      
      ! Apply physical BC for PLIC
      call apply_plic_bc()
      
   contains
      
      !> Apply physical BC to PLIC at domain boundaries
      subroutine apply_plic_bc()
         type(amrex_mfiter) :: mfi
         type(amrex_geometry) :: geom
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
         integer :: ilo, ihi, jlo, jhi, klo, khi
         integer :: dlo(3), dhi(3)
         real(WP) :: dx, dy, dz
         
         geom = this%amr%geom(lvl)
         dlo = geom%domain%lo
         dhi = geom%domain%hi
         dx = this%amr%dx(lvl)
         dy = this%amr%dy(lvl)
         dz = this%amr%dz(lvl)
         
         call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
         do while (mfi%next())
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
            pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            ilo = lbound(pPLIC,1); ihi = ubound(pPLIC,1)
            jlo = lbound(pPLIC,2); jhi = ubound(pPLIC,2)
            klo = lbound(pPLIC,3); khi = ubound(pPLIC,3)
            
            ! X-boundaries
            if (.not.this%amr%xper) then
               if (ilo.lt.dlo(1)) call apply_plic_bc_face(pVF, pCliq, pCgas, pPLIC, 1, -1, this%bc_xlo, &
               &   ilo, dlo(1)-1, jlo, jhi, klo, khi, dlo(1), this%amr%xlo, dx, dy, dz)
               if (ihi.gt.dhi(1)) call apply_plic_bc_face(pVF, pCliq, pCgas, pPLIC, 1, +1, this%bc_xhi, &
               &   dhi(1)+1, ihi, jlo, jhi, klo, khi, dhi(1), this%amr%xhi, dx, dy, dz)
            end if
            
            ! Y-boundaries
            if (.not.this%amr%yper) then
               if (jlo.lt.dlo(2)) call apply_plic_bc_face(pVF, pCliq, pCgas, pPLIC, 2, -1, this%bc_ylo, &
               &   ilo, ihi, jlo, dlo(2)-1, klo, khi, dlo(2), this%amr%ylo, dx, dy, dz)
               if (jhi.gt.dhi(2)) call apply_plic_bc_face(pVF, pCliq, pCgas, pPLIC, 2, +1, this%bc_yhi, &
               &   ilo, ihi, dhi(2)+1, jhi, klo, khi, dhi(2), this%amr%yhi, dx, dy, dz)
            end if
            
            ! Z-boundaries
            if (.not.this%amr%zper) then
               if (klo.lt.dlo(3)) call apply_plic_bc_face(pVF, pCliq, pCgas, pPLIC, 3, -1, this%bc_zlo, &
               &   ilo, ihi, jlo, jhi, klo, dlo(3)-1, dlo(3), this%amr%zlo, dx, dy, dz)
               if (khi.gt.dhi(3)) call apply_plic_bc_face(pVF, pCliq, pCgas, pPLIC, 3, +1, this%bc_zhi, &
               &   ilo, ihi, jlo, jhi, dhi(3)+1, khi, dhi(3), this%amr%zhi, dx, dy, dz)
            end if
            
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine apply_plic_bc
      
      !> Apply BC to PLIC on a single face
      subroutine apply_plic_bc_face(pVF, pCliq, pCgas, pPLIC, dir, side, bc_type, &
      &   i1, i2, j1, j2, k1, k2, bnd, x_bnd, dx, dy, dz)
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pVF, pCliq, pCgas, pPLIC
         integer, intent(in) :: dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: x_bnd, dx, dy, dz
         integer :: ig, jg, kg, isrc, jsrc, ksrc
         
         select case (bc_type)
         
          case (BC_LIQ)
            ! Trivial PLIC: full liquid
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               pPLIC(ig,jg,kg,1:3) = 0.0_WP
               pPLIC(ig,jg,kg,4) = 1.0e10_WP
            end do; end do; end do
            
          case (BC_GAS)
            ! Trivial PLIC: full gas
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               pPLIC(ig,jg,kg,1:3) = 0.0_WP
               pPLIC(ig,jg,kg,4) = -1.0e10_WP
            end do; end do; end do
            
          case (BC_REFLECT)
            ! Mirror PLIC from interior + flip normal component
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               isrc = ig; jsrc = jg; ksrc = kg
               if (dir.eq.1) isrc = 2*bnd - ig - side
               if (dir.eq.2) jsrc = 2*bnd - jg - side
               if (dir.eq.3) ksrc = 2*bnd - kg - side
               ! Copy plane
               pPLIC(ig,jg,kg,1:4) = pPLIC(isrc,jsrc,ksrc,1:4)
               ! Flip normal component
               pPLIC(ig,jg,kg,dir) = -pPLIC(ig,jg,kg,dir)
               ! Correct plane distance
               pPLIC(ig,jg,kg,4) = pPLIC(ig,jg,kg,4) - 2.0_WP*pPLIC(isrc,jsrc,ksrc,dir)*x_bnd
            end do; end do; end do
            
          case (BC_USER)
            ! User callback sets PLIC (user_vof_bc receives all arrays)
            if (associated(this%user_vof_bc)) then
               call this%user_vof_bc(this, lvl, time, dir, side, pVF, pCliq, pCgas, pPLIC, i1, i2, j1, j2, k1, k2)
            end if
            
          case default
            ! Do nothing
            
         end select
         
      end subroutine apply_plic_bc_face
      
   end subroutine fill_plic_lvl

   subroutine fill_moments_lvl(this, lvl, time)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry
      implicit none
      class(amrvof), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      
      ! Sync VF/Cliq/Cgas ghosts + fix periodic barycenters
      call this%sync_moments_lvl(lvl)
      
      ! Apply physical BC for moments
      call apply_moments_bc()
      
   contains
      
      !> Apply physical BC to VF/Cliq/Cgas at domain boundaries
      subroutine apply_moments_bc()
         type(amrex_mfiter) :: mfi
         type(amrex_geometry) :: geom
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
         integer :: ilo, ihi, jlo, jhi, klo, khi
         integer :: dlo(3), dhi(3)
         real(WP) :: dx, dy, dz
         
         geom = this%amr%geom(lvl)
         dlo = geom%domain%lo
         dhi = geom%domain%hi
         dx = this%amr%dx(lvl)
         dy = this%amr%dy(lvl)
         dz = this%amr%dz(lvl)
         
         call amrex_mfiter_build(mfi, this%VF%mf(lvl), tiling=.false.)
         do while (mfi%next())
            pVF => this%VF%mf(lvl)%dataptr(mfi)
            pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
            pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            ilo = lbound(pVF,1); ihi = ubound(pVF,1)
            jlo = lbound(pVF,2); jhi = ubound(pVF,2)
            klo = lbound(pVF,3); khi = ubound(pVF,3)
            
            ! X-boundaries
            if (.not.this%amr%xper) then
               if (ilo.lt.dlo(1)) call apply_moments_bc_face(pVF, pCliq, pCgas, pPLIC, 1, -1, this%bc_xlo, &
               &   ilo, dlo(1)-1, jlo, jhi, klo, khi, dlo(1), dx, dy, dz)
               if (ihi.gt.dhi(1)) call apply_moments_bc_face(pVF, pCliq, pCgas, pPLIC, 1, +1, this%bc_xhi, &
               &   dhi(1)+1, ihi, jlo, jhi, klo, khi, dhi(1), dx, dy, dz)
            end if
            
            ! Y-boundaries
            if (.not.this%amr%yper) then
               if (jlo.lt.dlo(2)) call apply_moments_bc_face(pVF, pCliq, pCgas, pPLIC, 2, -1, this%bc_ylo, &
               &   ilo, ihi, jlo, dlo(2)-1, klo, khi, dlo(2), dx, dy, dz)
               if (jhi.gt.dhi(2)) call apply_moments_bc_face(pVF, pCliq, pCgas, pPLIC, 2, +1, this%bc_yhi, &
               &   ilo, ihi, dhi(2)+1, jhi, klo, khi, dhi(2), dx, dy, dz)
            end if
            
            ! Z-boundaries
            if (.not.this%amr%zper) then
               if (klo.lt.dlo(3)) call apply_moments_bc_face(pVF, pCliq, pCgas, pPLIC, 3, -1, this%bc_zlo, &
               &   ilo, ihi, jlo, jhi, klo, dlo(3)-1, dlo(3), dx, dy, dz)
               if (khi.gt.dhi(3)) call apply_moments_bc_face(pVF, pCliq, pCgas, pPLIC, 3, +1, this%bc_zhi, &
               &   ilo, ihi, jlo, jhi, dhi(3)+1, khi, dhi(3), dx, dy, dz)
            end if
            
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine apply_moments_bc
      
      !> Apply BC to VF/Cliq/Cgas on a single face (moments only)
      !> For BC_USER, user callback receives all arrays including PLIC for convenience
      subroutine apply_moments_bc_face(pVF, pCliq, pCgas, pPLIC, dir, side, bc_type, &
      &   i1, i2, j1, j2, k1, k2, bnd, dx, dy, dz)
         real(WP), dimension(:,:,:,:), contiguous, intent(inout) :: pVF, pCliq, pCgas, pPLIC
         integer, intent(in) :: dir, side, bc_type, i1, i2, j1, j2, k1, k2, bnd
         real(WP), intent(in) :: dx, dy, dz
         integer :: ig, jg, kg, isrc, jsrc, ksrc
         real(WP), dimension(3) :: center
         
         select case (bc_type)
         
          case (BC_LIQ)
            ! Full liquid: VF=1, barycenters at cell center
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               center = [this%amr%xlo + (real(ig,WP)+0.5_WP)*dx, &
               &         this%amr%ylo + (real(jg,WP)+0.5_WP)*dy, &
               &         this%amr%zlo + (real(kg,WP)+0.5_WP)*dz]
               pVF(ig,jg,kg,1) = 1.0_WP
               pCliq(ig,jg,kg,1:3) = center
               pCgas(ig,jg,kg,1:3) = center
            end do; end do; end do
            
          case (BC_GAS)
            ! Full gas: VF=0, barycenters at cell center
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               center = [this%amr%xlo + (real(ig,WP)+0.5_WP)*dx, &
               &         this%amr%ylo + (real(jg,WP)+0.5_WP)*dy, &
               &         this%amr%zlo + (real(kg,WP)+0.5_WP)*dz]
               pVF(ig,jg,kg,1) = 0.0_WP
               pCliq(ig,jg,kg,1:3) = center
               pCgas(ig,jg,kg,1:3) = center
            end do; end do; end do
            
          case (BC_REFLECT)
            ! Mirror VF/Cliq/Cgas from interior
            do kg = k1, k2; do jg = j1, j2; do ig = i1, i2
               isrc = ig; jsrc = jg; ksrc = kg
               if (dir.eq.1) isrc = 2*bnd - ig - side
               if (dir.eq.2) jsrc = 2*bnd - jg - side
               if (dir.eq.3) ksrc = 2*bnd - kg - side
               pVF(ig,jg,kg,1) = pVF(isrc,jsrc,ksrc,1)
               pCliq(ig,jg,kg,1:3) = pCliq(isrc,jsrc,ksrc,1:3)
               pCgas(ig,jg,kg,1:3) = pCgas(isrc,jsrc,ksrc,1:3)
            end do; end do; end do
            
          case (BC_USER)
            ! User callback sets moments (and optionally PLIC)
            if (associated(this%user_vof_bc)) then
               call this%user_vof_bc(this, lvl, time, dir, side, pVF, pCliq, pCgas, pPLIC, i1, i2, j1, j2, k1, k2)
            end if
            
          case default
            ! Do nothing
            
         end select
         
      end subroutine apply_moments_bc_face
      
   end subroutine fill_moments_lvl

   !> Finalize the VOF solver
   subroutine finalize(this)
      class(amrvof), intent(inout) :: this
      call this%VF%finalize()
      call this%Cliq%finalize()
      call this%Cgas%finalize()
      call this%PLIC%finalize()
      call this%VFold%finalize()
      call this%Cliqold%finalize()
      call this%Cgasold%finalize()
      call this%PLICold%finalize()
      call this%smesh%reset()
      nullify(this%amr)
      nullify(this%user_init)
      nullify(this%user_tagging)
   end subroutine finalize

   ! ============================================================================
   ! VOF-SPECIFIC METHODS (STUBS)
   ! ============================================================================

   !> Build PLIC reconstruction from VF and barycenters using PLICnet
   !> Also extracts PLIC polygons and accumulates them to smesh
   subroutine build_plic(this, time)
      use plicnet, only: get_normal, reflect_moments
      use mathtools, only: normalize
      use amrvof_geometry, only: get_plane_dist, cut_hex_polygon
      class(amrvof), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl
      real(WP) :: dx, dy, dz
      
      ! Only build at finest level
      lvl = this%amr%clvl()
      
      ! Get cell size at this level
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      
      ! ========== Pass 1: Compute PLIC planes ==========
      plic_reconstruction: block
         integer :: i, j, k, ii, jj, kk, direction, direction2
         real(WP), dimension(0:188) :: moments
         real(WP), dimension(3) :: normal, center, lo, hi
         real(WP) :: m000, m100, m010, m001, temp, vf_cell
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
         logical :: flip
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
         
            ! Get pointers (with ghost cells for stencil access)
            pVF   => this%VF%mf(lvl)%dataptr(mfi)
            pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
            pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         
            ! Loop over cells in this box
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                  
                     vf_cell = pVF(i,j,k,1)
                  
                     ! Handle full cells: set trivial plane
                     if (vf_cell.lt.VFlo .or. vf_cell.gt.VFhi) then
                        pPLIC(i,j,k,1) = 0.0_WP  ! nx
                        pPLIC(i,j,k,2) = 0.0_WP  ! ny
                        pPLIC(i,j,k,3) = 0.0_WP  ! nz
                        pPLIC(i,j,k,4) = sign(1.0e10_WP, vf_cell - 0.5_WP)  ! d
                        cycle
                     end if
                  
                     ! Liquid-gas symmetry
                     flip = .false.
                     if (vf_cell.ge.0.5_WP) flip = .true.
                  
                     ! Initialize geometric moments
                     m000 = 0.0_WP; m100 = 0.0_WP; m010 = 0.0_WP; m001 = 0.0_WP
                  
                     ! Construct neighborhood of volume moments (3x3x3 stencil)
                     if (flip) then
                        do kk = k-1, k+1
                           do jj = j-1, j+1
                              do ii = i-1, i+1
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))   = 1.0_WP - pVF(ii,jj,kk,1)
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1) = (pCgas(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) / dx
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2) = (pCgas(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) / dy
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3) = (pCgas(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) / dz
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4) = (pCliq(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) / dx
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5) = (pCliq(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) / dy
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6) = (pCliq(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) / dz
                                 m000 = m000 + moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m100 = m100 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m010 = m010 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m001 = m001 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              end do
                           end do
                        end do
                     else
                        do kk = k-1, k+1
                           do jj = j-1, j+1
                              do ii = i-1, i+1
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))   = pVF(ii,jj,kk,1)
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1) = (pCliq(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) / dx
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2) = (pCliq(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) / dy
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3) = (pCliq(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) / dz
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4) = (pCgas(ii,jj,kk,1) - (this%amr%xlo + (real(ii,WP)+0.5_WP)*dx)) / dx
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5) = (pCgas(ii,jj,kk,2) - (this%amr%ylo + (real(jj,WP)+0.5_WP)*dy)) / dy
                                 moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6) = (pCgas(ii,jj,kk,3) - (this%amr%zlo + (real(kk,WP)+0.5_WP)*dz)) / dz
                                 m000 = m000 + moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m100 = m100 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m010 = m010 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                                 m001 = m001 + (moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k)) * moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                              end do
                           end do
                        end do
                     end if
                  
                     ! Geometric center of neighborhood
                     if (m000.gt.tiny(1.0_WP)) then
                        center = [m100, m010, m001] / m000
                     else
                        center = 0.0_WP
                     end if
                  
                     ! Apply symmetry (48 symmetries via reflect_moments)
                     call reflect_moments(moments, center, direction, direction2)
                  
                     ! Get normal from neural network
                     call get_normal(moments, normal)
                     normal = normalize(normal)
                  
                     ! Undo direction2 rotation (axis permutation)
                     if (direction2.eq.1) then
                        temp = normal(1); normal(1) = normal(2); normal(2) = temp
                     else if (direction2.eq.2) then
                        temp = normal(2); normal(2) = normal(3); normal(3) = temp
                     else if (direction2.eq.3) then
                        temp = normal(1); normal(1) = normal(3); normal(3) = temp
                     else if (direction2.eq.4) then
                        temp = normal(2); normal(2) = normal(3); normal(3) = temp
                        temp = normal(1); normal(1) = normal(2); normal(2) = temp
                     else if (direction2.eq.5) then
                        temp = normal(1); normal(1) = normal(3); normal(3) = temp
                        temp = normal(1); normal(1) = normal(2); normal(2) = temp
                     end if
                  
                     ! Undo direction reflection (octant)
                     if (direction.eq.1) then
                        normal(1) = -normal(1)
                     else if (direction.eq.2) then
                        normal(2) = -normal(2)
                     else if (direction.eq.3) then
                        normal(3) = -normal(3)
                     else if (direction.eq.4) then
                        normal(1) = -normal(1); normal(2) = -normal(2)
                     else if (direction.eq.5) then
                        normal(1) = -normal(1); normal(3) = -normal(3)
                     else if (direction.eq.6) then
                        normal(2) = -normal(2); normal(3) = -normal(3)
                     else if (direction.eq.7) then
                        normal(1) = -normal(1); normal(2) = -normal(2); normal(3) = -normal(3)
                     end if
                  
                     ! Undo liquid-gas flip
                     if (.not.flip) normal = -normal
                  
                     ! Renormalize
                     normal = normalize(normal)
                  
                     ! Cell bounds
                     lo = [this%amr%xlo + real(i,WP)*dx, this%amr%ylo + real(j,WP)*dy, this%amr%zlo + real(k,WP)*dz]
                     hi = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  
                     ! Store PLIC plane: (nx, ny, nz, d)
                     pPLIC(i,j,k,1) = normal(1)
                     pPLIC(i,j,k,2) = normal(2)
                     pPLIC(i,j,k,3) = normal(3)
                     pPLIC(i,j,k,4) = get_plane_dist(normal, lo, hi, vf_cell)
                  
                  end do
               end do
            end do
         
         end do
         call this%amr%mfiter_destroy(mfi)
      end block plic_reconstruction
      
      ! Fill PLIC ghosts (sync + periodic correction + physical BC)
      call this%fill_plic_lvl(lvl, time)
      
      ! ========== Pass 2: Per-FAB polygon extraction ==========
      call this%smesh%reset()
      
      polygon_extraction: block
         integer :: i, j, k, ii, jj, kk
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
         real(WP), dimension(3) :: lo, hi
         real(WP), dimension(4) :: plane
         real(WP), dimension(3,8) :: hex
         real(WP) :: vf_cell
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx, gbx
         
         ! Per-FAB polygon storage (allocatable, indexed by cell)
         real(WP), dimension(:,:,:,:,:), allocatable :: polygon_local  ! (3, 6, ilo:ihi, jlo:jhi, klo:khi)
         integer, dimension(:,:,:), allocatable :: poly_nv_local       ! (ilo:ihi, jlo:jhi, klo:khi)
         real(WP), dimension(3,6) :: poly_verts
         integer :: poly_nv
         integer :: glo(3), ghi(3)
         
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            gbx = mfi%growntilebox(2)  ! Grown by 2 for 5x5x5 stencil
            pVF   => this%VF%mf(lvl)%dataptr(mfi)
            pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
            
            ! Get bounds for per-FAB allocation
            glo = gbx%lo
            ghi = gbx%hi
            
            ! ----- Step A: Allocate per-FAB polygon storage -----
            allocate(polygon_local(1:3, 1:6, glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)))
            allocate(poly_nv_local(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3)))
            polygon_local = 0.0_WP
            poly_nv_local = 0
            
            ! ----- Step B: Extract polygons (grown box including ghosts) -----
            do kk = glo(3), ghi(3)
               do jj = glo(2), ghi(2)
                  do ii = glo(1), ghi(1)
                     vf_cell = pVF(ii,jj,kk,1)
                     if (vf_cell.lt.VFlo .or. vf_cell.gt.VFhi) cycle
                     
                     ! Build hex and plane
                     lo = [this%amr%xlo + real(ii,WP)*dx, this%amr%ylo + real(jj,WP)*dy, this%amr%zlo + real(kk,WP)*dz]
                     hi = [this%amr%xlo + real(ii+1,WP)*dx, this%amr%ylo + real(jj+1,WP)*dy, this%amr%zlo + real(kk+1,WP)*dz]
                     plane = [pPLIC(ii,jj,kk,1), pPLIC(ii,jj,kk,2), pPLIC(ii,jj,kk,3), pPLIC(ii,jj,kk,4)]
                     hex(:,1) = [hi(1), lo(2), lo(3)]
                     hex(:,2) = [hi(1), hi(2), lo(3)]
                     hex(:,3) = [hi(1), hi(2), hi(3)]
                     hex(:,4) = [hi(1), lo(2), hi(3)]
                     hex(:,5) = [lo(1), lo(2), lo(3)]
                     hex(:,6) = [lo(1), hi(2), lo(3)]
                     hex(:,7) = [lo(1), hi(2), hi(3)]
                     hex(:,8) = [lo(1), lo(2), hi(3)]
                     
                     call cut_hex_polygon(hex, plane, poly_nv, poly_verts)
                     
                     ! Store in per-FAB array
                     poly_nv_local(ii,jj,kk) = poly_nv
                     if (poly_nv.ge.3) then
                        polygon_local(:,1:poly_nv,ii,jj,kk) = poly_verts(:,1:poly_nv)
                     end if
                  end do
               end do
            end do
            
            ! ----- Step C: Compute curvature (valid cells, stencil access) -----
            ! TODO: curvature = f(polygon_local stencil around i,j,k)
            ! For now, skip curvature computation
            
            ! ----- Step D: Append to smesh (valid cells only) -----
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     poly_nv = poly_nv_local(i,j,k)
                     if (poly_nv.ge.3) then
                        call this%smesh%add_polygon(polygon_local(:,1:poly_nv,i,j,k), poly_nv)
                     end if
                  end do
               end do
            end do
            
            ! ----- Step E: Deallocate per-FAB storage -----
            deallocate(polygon_local, poly_nv_local)
            
         end do
         call this%amr%mfiter_destroy(mfi)
      end block polygon_extraction
      
   end subroutine build_plic
   
   !> Reset VF and barycenters from PLIC plane to ensure consistency
   !> Computes in valid + ghost cells from PLIC (which is already filled)
   !> Then averages down to coarse levels
   subroutine reset_moments(this)
      use amrvof_geometry, only: cut_hex_vol, VFlo, VFhi
      class(amrvof), intent(inout) :: this
      integer :: lvl, i, j, k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane
      real(WP) :: vol_liq, vol_gas, cell_vol, dx, dy, dz
      real(WP), dimension(3) :: bary_liq, bary_gas, cell_center
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: ilo, ihi, jlo, jhi, klo, khi
      
      ! Only work at finest level
      lvl = this%amr%clvl()
      
      dx = this%amr%dx(lvl)
      dy = this%amr%dy(lvl)
      dz = this%amr%dz(lvl)
      cell_vol = dx * dy * dz
      
      ! Use tiling=.false. to iterate over full FABs including ghosts
      call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
      do while (mfi%next())
         bx = mfi%fabbox()  ! Full box including ghosts
         
         pVF   => this%VF%mf(lvl)%dataptr(mfi)
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         
         ! Get array bounds (same as bx but from array)
         ilo = lbound(pVF,1); ihi = ubound(pVF,1)
         jlo = lbound(pVF,2); jhi = ubound(pVF,2)
         klo = lbound(pVF,3); khi = ubound(pVF,3)
         
         do k = klo, khi
            do j = jlo, jhi
               do i = ilo, ihi
                  
                  ! Cell center
                  cell_center = [this%amr%xlo + (real(i,WP)+0.5_WP)*dx, &
                  &              this%amr%ylo + (real(j,WP)+0.5_WP)*dy, &
                  &              this%amr%zlo + (real(k,WP)+0.5_WP)*dz]
                  
                  ! Build hex cell (8 vertices, standard ordering)
                  hex(:,1) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                  hex(:,2) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                  hex(:,3) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                  hex(:,4) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k  ,WP)*dz]
                  hex(:,5) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  hex(:,6) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j  ,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  hex(:,7) = [this%amr%xlo + real(i+1,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  hex(:,8) = [this%amr%xlo + real(i  ,WP)*dx, this%amr%ylo + real(j+1,WP)*dy, this%amr%zlo + real(k+1,WP)*dz]
                  
                  ! Get plane from PLIC (already filled including ghosts)
                  plane(1:3) = pPLIC(i,j,k,1:3)
                  plane(4)   = pPLIC(i,j,k,4)
                  
                  ! Skip cutting for full cells (trivial PLIC with large distance)
                  if (abs(plane(4)) .ge. 1.0e9_WP) then
                     if (plane(4) .gt. 0.0_WP) then
                        pVF(i,j,k,1) = 1.0_WP  ! All liquid
                     else
                        pVF(i,j,k,1) = 0.0_WP  ! All gas
                     end if
                     pCliq(i,j,k,1:3) = cell_center
                     pCgas(i,j,k,1:3) = cell_center
                     cycle
                  end if
                  
                  ! Cut hex by plane
                  call cut_hex_vol(hex, plane, vol_liq, vol_gas, bary_liq, bary_gas)
                  
                  ! Update VF
                  pVF(i,j,k,1) = vol_liq / cell_vol
                  
                  ! Update barycenters
                  pCliq(i,j,k,1:3) = bary_liq
                  pCgas(i,j,k,1:3) = bary_gas
                  
                  ! Clean up edge cases
                  if (pVF(i,j,k,1).lt.VFlo) then
                     pVF(i,j,k,1) = 0.0_WP
                     pCliq(i,j,k,1:3) = cell_center
                     pCgas(i,j,k,1:3) = cell_center
                  end if
                  if (pVF(i,j,k,1).gt.VFhi) then
                     pVF(i,j,k,1) = 1.0_WP
                     pCliq(i,j,k,1:3) = cell_center
                     pCgas(i,j,k,1:3) = cell_center
                  end if
                  
               end do
            end do
         end do
         
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Average down to coarse levels (uses ghost-capable average_down)
      call this%average_down_moments()
      
   end subroutine reset_moments

   !> Advect VF using velocity field (staggered or collocated, auto-detected from nodality)
   !> User must provide MultiFabs at finest level with >= 2 ghost cells filled
   subroutine advance_vof(this, U, V, W, dt, time)
      use amrvof_geometry, only: cut_tet_vol
      use amrex_amr_module, only: amrex_multifab
      implicit none
      class(amrvof), intent(inout) :: this
      type(amrex_multifab), intent(in) :: U, V, W
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      logical :: is_staggered
      integer :: lvl, i, j, k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold, pCliqold, pCgasold, pPLICold
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      real(WP), dimension(:,:,:,:), allocatable :: Fx, Fy, Fz  ! Volume fluxes (Lvol, Gvol, Lbar(3), Gbar(3))
      real(WP) :: Lvol_old, Gvol_old, Lvol_new, Gvol_new, Lvol_flux, Gvol_flux
      real(WP), dimension(3) :: Lbar_old, Gbar_old, Lbar_new, Gbar_new, Lbar_flux, Gbar_flux
      real(WP) :: dx, dy, dz, vol, face_vel
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: ilo, ihi, jlo, jhi, klo, khi
      ! Shared variables for internal functions
      real(WP), dimension(3,9) :: face     ! 4 current + 4 projected + 1 center
      real(WP) :: dxi, dyi, dzi
      integer :: ilo_, ihi_, jlo_, jhi_, klo_, khi_  ! Velocity bounds

      ! Check velocity centering and ghost cell requirements
      check_velocity: block
         use messager, only: die
         logical, dimension(3) :: nodal_U, nodal_V, nodal_W
         nodal_U = U%nodal_type()
         nodal_V = V%nodal_type()
         nodal_W = W%nodal_type()
         if (all(nodal_U .eqv. [.true. , .false., .false.]) .and. & 
         &   all(nodal_V .eqv. [.false., .true. , .false.]) .and. & 
         &   all(nodal_W .eqv. [.false., .false., .true. ])) then
            is_staggered = .true.
         else if (.not.any(nodal_U) .and. .not.any(nodal_V) .and. .not.any(nodal_W)) then
            is_staggered = .false.
         else
            call die('[advance_vof] velocity must be either staggered (face-centered) or collocated (cell-centered)')
         end if
         if (U%nghost() .lt. 2 .or. V%nghost() .lt. 2 .or. W%nghost() .lt. 2) then
            call die('[advance_vof] velocity requires >= 2 ghost cells')
         end if
      end block check_velocity
      
      ! Only advect at finest level
      lvl = this%amr%clvl()
      
      ! Get cell size
      dx = this%amr%dx(lvl); dxi=1.0_WP/dx
      dy = this%amr%dy(lvl); dyi=1.0_WP/dy
      dz = this%amr%dz(lvl); dzi=1.0_WP/dz
      vol = dx * dy * dz
      
      ! Iterate over boxes
      call this%amr%mfiter_build(lvl, mfi, tiling=.false.)
      do while (mfi%next())
         bx = mfi%tilebox()
         
         ! Get pointers
         pVF => this%VF%mf(lvl)%dataptr(mfi)
         pCliq => this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC => this%PLIC%mf(lvl)%dataptr(mfi)
         pVFold => this%VFold%mf(lvl)%dataptr(mfi)
         pCliqold => this%Cliqold%mf(lvl)%dataptr(mfi)
         pCgasold => this%Cgasold%mf(lvl)%dataptr(mfi)
         pPLICold => this%PLICold%mf(lvl)%dataptr(mfi)
         pU => U%dataptr(mfi)
         pV => V%dataptr(mfi)
         pW => W%dataptr(mfi)
         
         ilo = lbound(pVF,1); ihi = ubound(pVF,1)
         jlo = lbound(pVF,2); jhi = ubound(pVF,2)
         klo = lbound(pVF,3); khi = ubound(pVF,3)
         
         ! Allocate flux arrays (8 components: Lvol, Gvol, Lbar(3), Gbar(3))
         allocate(Fx(8,ilo:ihi+1,jlo:jhi,klo:khi)); Fx = 0.0_WP
         allocate(Fy(8,ilo:ihi,jlo:jhi+1,klo:khi)); Fy = 0.0_WP
         allocate(Fz(8,ilo:ihi,jlo:jhi,klo:khi+1)); Fz = 0.0_WP
         
         ! Compute fluxes at each face
         do k = bx%lo(3), bx%hi(3)+1
            do j = bx%lo(2), bx%hi(2)+1
               do i = bx%lo(1), bx%hi(1)+1
                  
                  ! X-flux at face (i,j,k)
                  if (i.le.bx%hi(1)+1 .and. j.le.bx%hi(2) .and. k.le.bx%hi(3)) then
                     if (is_staggered) then
                        face_vel = pU(i,j,k,1)
                     else
                        face_vel = 0.5_WP*(pU(i-1,j,k,1) + pU(i,j,k,1))
                     end if
                     call compute_flux(i, j, k, 1, face_vel, Fx(:,i,j,k))
                  end if
                  
                  ! Y-flux at face (i,j,k)
                  if (i.le.bx%hi(1) .and. j.le.bx%hi(2)+1 .and. k.le.bx%hi(3)) then
                     if (is_staggered) then
                        face_vel = pV(i,j,k,1)
                     else
                        face_vel = 0.5_WP*(pV(i,j-1,k,1) + pV(i,j,k,1))
                     end if
                     call compute_flux(i, j, k, 2, face_vel, Fy(:,i,j,k))
                  end if
                  
                  ! Z-flux at face (i,j,k)
                  if (i.le.bx%hi(1) .and. j.le.bx%hi(2) .and. k.le.bx%hi(3)+1) then
                     if (is_staggered) then
                        face_vel = pW(i,j,k,1)
                     else
                        face_vel = 0.5_WP*(pW(i,j,k-1,1) + pW(i,j,k,1))
                     end if
                     call compute_flux(i, j, k, 3, face_vel, Fz(:,i,j,k))
                  end if
                  
               end do
            end do
         end do
         
         ! Update volume moments
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  ! Old volumes
                  Lvol_old = pVFold(i,j,k,1) * vol
                  Gvol_old = (1.0_WP - pVFold(i,j,k,1)) * vol
                  Lbar_old = pCliqold(i,j,k,1:3)
                  Gbar_old = pCgasold(i,j,k,1:3)
                  
                  ! Net flux (outflow positive)
                  Lvol_flux = Fx(1,i+1,j,k) - Fx(1,i,j,k) + Fy(1,i,j+1,k) - Fy(1,i,j,k) + Fz(1,i,j,k+1) - Fz(1,i,j,k)
                  Gvol_flux = Fx(2,i+1,j,k) - Fx(2,i,j,k) + Fy(2,i,j+1,k) - Fy(2,i,j,k) + Fz(2,i,j,k+1) - Fz(2,i,j,k)
                  Lbar_flux = Fx(3:5,i+1,j,k) - Fx(3:5,i,j,k) + Fy(3:5,i,j+1,k) - Fy(3:5,i,j,k) + Fz(3:5,i,j,k+1) - Fz(3:5,i,j,k)
                  Gbar_flux = Fx(6:8,i+1,j,k) - Fx(6:8,i,j,k) + Fy(6:8,i,j+1,k) - Fy(6:8,i,j,k) + Fz(6:8,i,j,k+1) - Fz(6:8,i,j,k)
                  
                  ! New volumes
                  Lvol_new = Lvol_old - Lvol_flux
                  Gvol_new = Gvol_old - Gvol_flux
                  
                  ! New VF
                  pVF(i,j,k,1) = Lvol_new / (Lvol_new + Gvol_new)
                  
                  ! Clip and update barycenters
                  if (pVF(i,j,k,1) .lt. VFlo) then
                     pVF(i,j,k,1) = 0.0_WP
                     pCliq(i,j,k,1:3) = [this%amr%xlo + (real(i,WP)+0.5_WP)*dx, &
                     &                   this%amr%ylo + (real(j,WP)+0.5_WP)*dy, &
                     &                   this%amr%zlo + (real(k,WP)+0.5_WP)*dz]
                     pCgas(i,j,k,1:3) = pCliq(i,j,k,1:3)
                  else if (pVF(i,j,k,1) .gt. VFhi) then
                     pVF(i,j,k,1) = 1.0_WP
                     pCliq(i,j,k,1:3) = [this%amr%xlo + (real(i,WP)+0.5_WP)*dx, &
                     &                   this%amr%ylo + (real(j,WP)+0.5_WP)*dy, &
                     &                   this%amr%zlo + (real(k,WP)+0.5_WP)*dz]
                     pCgas(i,j,k,1:3) = pCliq(i,j,k,1:3)
                  else
                     ! Update barycenters from moment conservation
                     Lbar_new = (Lbar_old * Lvol_old - Lbar_flux) / Lvol_new
                     Gbar_new = (Gbar_old * Gvol_old - Gbar_flux) / Gvol_new
                     pCliq(i,j,k,1:3) = Lbar_new
                     pCgas(i,j,k,1:3) = Gbar_new
                  end if
               end do
            end do
         end do
         
         deallocate(Fx, Fy, Fz)
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Sync and apply BC
      call this%fill_moments_lvl(lvl, time)
      
   contains
      
      !> Compute flux through a face using full semi-Lagrangian RK2 with recursive tet cutting
      subroutine compute_flux(fi, fj, fk, dir, fv, flux)
         use amrvof_geometry, only: tet_map, cut_side, cut_v1, cut_v2, cut_vtet, cut_ntets, cut_nvert
         integer, intent(in) :: fi, fj, fk, dir
         real(WP), intent(in) :: fv
         real(WP), dimension(8), intent(out) :: flux
         real(WP), dimension(3,4) :: tetra
         integer, dimension(3,4) :: myijk
         real(WP) :: vol_expected, my_tet_sign
         integer :: src_i, src_j, src_k, ntet, nn
         real(WP), dimension(8) :: tet_flux
         
         flux = 0.0_WP
         
         ! Precompute inverse cell sizes and velocity array bounds
         dxi = 1.0_WP / dx; dyi = 1.0_WP / dy; dzi = 1.0_WP / dz
         ilo_ = lbound(pU,1); ihi_ = ubound(pU,1) - 1
         jlo_ = lbound(pU,2); jhi_ = ubound(pU,2) - 1
         klo_ = lbound(pU,3); khi_ = ubound(pU,3) - 1
         
         ! Determine upwind cell
         if (dir .eq. 1) then
            if (fv .ge. 0.0_WP) then
               src_i = fi-1; src_j = fj; src_k = fk
            else
               src_i = fi; src_j = fj; src_k = fk
            end if
         else if (dir .eq. 2) then
            if (fv .ge. 0.0_WP) then
               src_i = fi; src_j = fj-1; src_k = fk
            else
               src_i = fi; src_j = fj; src_k = fk
            end if
         else
            if (fv .ge. 0.0_WP) then
               src_i = fi; src_j = fj; src_k = fk-1
            else
               src_i = fi; src_j = fj; src_k = fk
            end if
         end if
         
         ! Skip pure cells - fast path
         if (pVFold(src_i,src_j,src_k,1) .lt. VFlo .or. pVFold(src_i,src_j,src_k,1) .gt. VFhi) then
            vol_expected = abs(fv) * dt
            if (dir .eq. 1) vol_expected = vol_expected * dy * dz
            if (dir .eq. 2) vol_expected = vol_expected * dx * dz
            if (dir .eq. 3) vol_expected = vol_expected * dx * dy
            if (pVFold(src_i,src_j,src_k,1) .lt. VFlo) then
               flux(2) = sign(vol_expected, fv)  ! Gas volume
            else
               flux(1) = sign(vol_expected, fv)  ! Liquid volume
            end if
            return
         end if
         
         ! Build face vertices (1-4 at current position, 5-8 projected back)
         if (dir .eq. 1) then
            face(:,1) = [this%amr%xlo + real(fi,WP)*dx, this%amr%ylo + real(fj+1,WP)*dy, this%amr%zlo + real(fk  ,WP)*dz]
            face(:,2) = [this%amr%xlo + real(fi,WP)*dx, this%amr%ylo + real(fj+1,WP)*dy, this%amr%zlo + real(fk+1,WP)*dz]
            face(:,3) = [this%amr%xlo + real(fi,WP)*dx, this%amr%ylo + real(fj  ,WP)*dy, this%amr%zlo + real(fk+1,WP)*dz]
            face(:,4) = [this%amr%xlo + real(fi,WP)*dx, this%amr%ylo + real(fj  ,WP)*dy, this%amr%zlo + real(fk  ,WP)*dz]
            vol_expected = -fv * dt * dy * dz
         else if (dir .eq. 2) then
            face(:,1) = [this%amr%xlo + real(fi+1,WP)*dx, this%amr%ylo + real(fj,WP)*dy, this%amr%zlo + real(fk  ,WP)*dz]
            face(:,2) = [this%amr%xlo + real(fi  ,WP)*dx, this%amr%ylo + real(fj,WP)*dy, this%amr%zlo + real(fk  ,WP)*dz]
            face(:,3) = [this%amr%xlo + real(fi  ,WP)*dx, this%amr%ylo + real(fj,WP)*dy, this%amr%zlo + real(fk+1,WP)*dz]
            face(:,4) = [this%amr%xlo + real(fi+1,WP)*dx, this%amr%ylo + real(fj,WP)*dy, this%amr%zlo + real(fk+1,WP)*dz]
            vol_expected = -fv * dt * dx * dz
         else
            face(:,1) = [this%amr%xlo + real(fi  ,WP)*dx, this%amr%ylo + real(fj+1,WP)*dy, this%amr%zlo + real(fk,WP)*dz]
            face(:,2) = [this%amr%xlo + real(fi  ,WP)*dx, this%amr%ylo + real(fj  ,WP)*dy, this%amr%zlo + real(fk,WP)*dz]
            face(:,3) = [this%amr%xlo + real(fi+1,WP)*dx, this%amr%ylo + real(fj  ,WP)*dy, this%amr%zlo + real(fk,WP)*dz]
            face(:,4) = [this%amr%xlo + real(fi+1,WP)*dx, this%amr%ylo + real(fj+1,WP)*dy, this%amr%zlo + real(fk,WP)*dz]
            vol_expected = -fv * dt * dx * dy
         end if
         
         ! Project vertices back in time with RK2
         face(:,5) = project(face(:,1), -dt)
         face(:,6) = project(face(:,2), -dt)
         face(:,7) = project(face(:,3), -dt)
         face(:,8) = project(face(:,4), -dt)
         
         ! Center point (will be adjusted by volume_correct)
         face(:,9) = 0.25_WP * (face(:,5) + face(:,6) + face(:,7) + face(:,8))
         
         ! Apply volume correction for exact flux
         if (dir .eq. 1) then
            call volume_correct_x(vol_expected)
         else if (dir .eq. 2) then
            call volume_correct_y(vol_expected)
         else
            call volume_correct_z(vol_expected)
         end if
         
         ! Process each of the 8 tets with recursive cutting
         do ntet = 1, 8
            ! Build tet from face vertices
            do nn = 1, 4
               tetra(:,nn) = face(:,tet_map(nn,ntet))
               if (dir .eq. 1) then
                  myijk(1,nn) = src_i
                  myijk(2,nn) = min(max(floor((tetra(2,nn) - this%amr%ylo) * dyi), jlo_), jhi_)
                  myijk(3,nn) = min(max(floor((tetra(3,nn) - this%amr%zlo) * dzi), klo_), khi_)
               else if (dir .eq. 2) then
                  myijk(1,nn) = min(max(floor((tetra(1,nn) - this%amr%xlo) * dxi), ilo_), ihi_)
                  myijk(2,nn) = src_j
                  myijk(3,nn) = min(max(floor((tetra(3,nn) - this%amr%zlo) * dzi), klo_), khi_)
               else
                  myijk(1,nn) = min(max(floor((tetra(1,nn) - this%amr%xlo) * dxi), ilo_), ihi_)
                  myijk(2,nn) = min(max(floor((tetra(2,nn) - this%amr%ylo) * dyi), jlo_), jhi_)
                  myijk(3,nn) = src_k
               end if
            end do
            
            ! Compute tet sign
            my_tet_sign = tet_sign_func(tetra)
            
            ! Recursively cut tet by grid and PLIC, accumulate flux
            tet_flux = tet2flux(tetra, myijk)
            
            ! Apply flux sign based on direction
            if (dir .eq. 1) then
               flux = flux + my_tet_sign * tet_flux
            else if (dir .eq. 2) then
               flux = flux - my_tet_sign * tet_flux
            else
               flux = flux + my_tet_sign * tet_flux
            end if
         end do
         
      end subroutine compute_flux
      
      !> Compute tet sign (positive = right-handed, negative = left-handed)
      pure function tet_sign_func(v) result(s)
         real(WP), dimension(3,4), intent(in) :: v
         real(WP) :: s
         real(WP), dimension(3) :: a, b, c
         a = v(:,1) - v(:,4)
         b = v(:,2) - v(:,4)
         c = v(:,3) - v(:,4)
         s = sign(1.0_WP, -(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP)
      end function tet_sign_func
      
      !> Recursive function that cuts a tet by grid planes to compute fluxes
      recursive function tet2flux(mytet, myind) result(myflux)
         use amrvof_geometry, only: cut_side, cut_v1, cut_v2, cut_vtet, cut_ntets, cut_nvert
         real(WP), dimension(3,4), intent(in) :: mytet
         integer, dimension(3,4), intent(in) :: myind
         real(WP), dimension(8) :: myflux
         integer :: dir, cut_ind, icase, n1, n2, v1, v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         integer, dimension(3,8,2) :: vert_ind
         real(WP) :: denom, mu, my_vol
         real(WP), dimension(3,4) :: newtet
         integer, dimension(3,4) :: newind
         real(WP), dimension(3) :: a, b, c
         real(WP) :: xcut, ycut, zcut
         
         myflux = 0.0_WP
         
         ! Determine if tet spans multiple cells and needs cutting
         if (maxval(myind(1,:)) - minval(myind(1,:)) .gt. 0) then
            ! Cut by x planes
            dir = 1
            cut_ind = maxval(myind(1,:))
            xcut = this%amr%xlo + real(cut_ind,WP) * dx
            dd(:) = mytet(1,:) - xcut
         else if (maxval(myind(2,:)) - minval(myind(2,:)) .gt. 0) then
            ! Cut by y planes
            dir = 2
            cut_ind = maxval(myind(2,:))
            ycut = this%amr%ylo + real(cut_ind,WP) * dy
            dd(:) = mytet(2,:) - ycut
         else if (maxval(myind(3,:)) - minval(myind(3,:)) .gt. 0) then
            ! Cut by z planes
            dir = 3
            cut_ind = maxval(myind(3,:))
            zcut = this%amr%zlo + real(cut_ind,WP) * dz
            dd(:) = mytet(3,:) - zcut
         else
            ! All vertices in same cell - cut by PLIC and return
            myflux = tet2flux_plic(mytet, myind(1,1), myind(2,1), myind(3,1))
            return
         end if
         
         ! Find cut case (1-indexed: 1-16)
         icase = 1 + int(0.5_WP + sign(0.5_WP, dd(1))) &
               + 2 * int(0.5_WP + sign(0.5_WP, dd(2))) &
               + 4 * int(0.5_WP + sign(0.5_WP, dd(3))) &
               + 8 * int(0.5_WP + sign(0.5_WP, dd(4)))
         
         ! Copy vertices and indices
         do n1 = 1, 4
            vert(:, n1) = mytet(:, n1)
            vert_ind(:, n1, 1) = myind(:, n1)
            vert_ind(:, n1, 2) = myind(:, n1)
            ! Enforce boundedness at cut plane
            vert_ind(dir, n1, 1) = min(vert_ind(dir, n1, 1), cut_ind - 1)
            vert_ind(dir, n1, 2) = max(vert_ind(dir, n1, 1), cut_ind)
         end do
         
         ! Create interpolated vertices on cut plane
         do n1 = 1, cut_nvert(icase)
            v1 = cut_v1(n1, icase); v2 = cut_v2(n1, icase)
            denom = dd(v2) - dd(v1)
            if (abs(denom) .ge. tiny(1.0_WP)) then
               mu = max(0.0_WP, min(1.0_WP, -dd(v1) / denom))
            else
               mu = 0.0_WP
            end if
            vert(:, 4 + n1) = (1.0_WP - mu) * vert(:, v1) + mu * vert(:, v2)
            ! Compute index for interpolated vertex
            vert_ind(1, 4+n1, 1) = min(max(floor((vert(1, 4+n1) - this%amr%xlo) * dxi), ilo_), ihi_)
            vert_ind(2, 4+n1, 1) = min(max(floor((vert(2, 4+n1) - this%amr%ylo) * dyi), jlo_), jhi_)
            vert_ind(3, 4+n1, 1) = min(max(floor((vert(3, 4+n1) - this%amr%zlo) * dzi), klo_), khi_)
            ! Enforce boundedness
            vert_ind(:, 4+n1, 1) = max(vert_ind(:, 4+n1, 1), min(vert_ind(:, v1, 1), vert_ind(:, v2, 1)))
            vert_ind(:, 4+n1, 1) = min(vert_ind(:, 4+n1, 1), max(vert_ind(:, v1, 1), vert_ind(:, v2, 1)))
            ! Set +/- indices in cut direction
            vert_ind(:, 4+n1, 2) = vert_ind(:, 4+n1, 1)
            vert_ind(dir, 4+n1, 1) = cut_ind - 1
            vert_ind(dir, 4+n1, 2) = cut_ind
         end do
         
         ! Create and process sub-tets
         do n1 = 1, cut_ntets(icase)
            do n2 = 1, 4
               newtet(:, n2) = vert(:, cut_vtet(n2, n1, icase))
               newind(:, n2) = vert_ind(:, cut_vtet(n2, n1, icase), cut_side(n1, icase))
            end do
            ! Check for zero-volume tet
            a = newtet(:,1) - newtet(:,4)
            b = newtet(:,2) - newtet(:,4)
            c = newtet(:,3) - newtet(:,4)
            my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            if (my_vol .lt. 1.0e-20_WP) cycle
            ! Recursively process sub-tet
            myflux = myflux + tet2flux(newtet, newind)
         end do
         
      end function tet2flux
      
      !> Cut tet by PLIC and compute flux (base case of recursion)
      function tet2flux_plic(mytet, i0, j0, k0) result(myflux)
         use amrvof_geometry, only: cut_v1, cut_v2, cut_vtet, cut_ntets, cut_nvert, cut_nntet
         real(WP), dimension(3,4), intent(in) :: mytet
         integer, intent(in) :: i0, j0, k0
         real(WP), dimension(8) :: myflux
         integer :: icase, n1, v1, v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         real(WP), dimension(3) :: a, b, c, bary, normal
         real(WP) :: denom, mu, my_vol, dist
         
         myflux = 0.0_WP
         
         ! Get PLIC from this cell
         normal = pPLICold(i0, j0, k0, 1:3)
         dist = pPLICold(i0, j0, k0, 4)
         
         ! Compute signed distance to plane for each vertex
         dd(1) = normal(1)*mytet(1,1) + normal(2)*mytet(2,1) + normal(3)*mytet(3,1) - dist
         dd(2) = normal(1)*mytet(1,2) + normal(2)*mytet(2,2) + normal(3)*mytet(3,2) - dist
         dd(3) = normal(1)*mytet(1,3) + normal(2)*mytet(2,3) + normal(3)*mytet(3,3) - dist
         dd(4) = normal(1)*mytet(1,4) + normal(2)*mytet(2,4) + normal(3)*mytet(3,4) - dist
         
         ! Find cut case
         icase = 1 + int(0.5_WP + sign(0.5_WP, dd(1))) &
               + 2 * int(0.5_WP + sign(0.5_WP, dd(2))) &
               + 4 * int(0.5_WP + sign(0.5_WP, dd(3))) &
               + 8 * int(0.5_WP + sign(0.5_WP, dd(4)))
         
         ! Copy vertices
         vert(:, 1:4) = mytet(:, 1:4)
         
         ! Create interpolated vertices on cut plane
         do n1 = 1, cut_nvert(icase)
            v1 = cut_v1(n1, icase); v2 = cut_v2(n1, icase)
            denom = dd(v2) - dd(v1)
            if (abs(denom) .ge. tiny(1.0_WP)) then
               mu = max(0.0_WP, min(1.0_WP, -dd(v1) / denom))
            else
               mu = 0.0_WP
            end if
            vert(:, 4 + n1) = (1.0_WP - mu) * vert(:, v1) + mu * vert(:, v2)
         end do
         
         ! Gas tets: from 1 to cut_nntet-1
         do n1 = 1, cut_nntet(icase) - 1
            a = vert(:, cut_vtet(1, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            b = vert(:, cut_vtet(2, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            c = vert(:, cut_vtet(3, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            bary = 0.25_WP * (vert(:, cut_vtet(1, n1, icase)) + vert(:, cut_vtet(2, n1, icase)) &
            &               + vert(:, cut_vtet(3, n1, icase)) + vert(:, cut_vtet(4, n1, icase)))
            myflux(2) = myflux(2) + my_vol
            myflux(6:8) = myflux(6:8) + my_vol * bary
         end do
         
         ! Liquid tets: from cut_ntets down to cut_nntet
         do n1 = cut_ntets(icase), cut_nntet(icase), -1
            a = vert(:, cut_vtet(1, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            b = vert(:, cut_vtet(2, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            c = vert(:, cut_vtet(3, n1, icase)) - vert(:, cut_vtet(4, n1, icase))
            my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            bary = 0.25_WP * (vert(:, cut_vtet(1, n1, icase)) + vert(:, cut_vtet(2, n1, icase)) &
            &               + vert(:, cut_vtet(3, n1, icase)) + vert(:, cut_vtet(4, n1, icase)))
            myflux(1) = myflux(1) + my_vol
            myflux(3:5) = myflux(3:5) + my_vol * bary
         end do
         
      end function tet2flux_plic
      
      !> RK2 vertex projection back in time
      function project(p1,mydt) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3)             :: p2
         real(WP),               intent(in) :: mydt
         p2=p1+mydt*interp_velocity(        p1    )
         p2=p1+mydt*interp_velocity(0.5_WP*(p1+p2))
      end function project
      
      !> Trilinear interpolation of velocity (handles staggered or collocated)
      function interp_velocity(pos) result(vel)
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer :: ip, jp, kp
         real(WP) :: wx1, wy1, wz1, wx2, wy2, wz2
         integer :: ipu, jpv, kpw
         real(WP) :: wxu1, wyv1, wzw1, wxu2, wyv2, wzw2
         
         ! Cell-centered indices and weights (always needed)
         ip = max(ilo_, min(ihi_, floor((pos(1) - this%amr%xlo) * dxi)))
         jp = max(jlo_, min(jhi_, floor((pos(2) - this%amr%ylo) * dyi)))
         kp = max(klo_, min(khi_, floor((pos(3) - this%amr%zlo) * dzi)))
         wx1 = max(0.0_WP, min(1.0_WP, (pos(1) - (this%amr%xlo + (real(ip,WP)+0.5_WP)*dx)) * dxi + 0.5_WP))
         wy1 = max(0.0_WP, min(1.0_WP, (pos(2) - (this%amr%ylo + (real(jp,WP)+0.5_WP)*dy)) * dyi + 0.5_WP))
         wz1 = max(0.0_WP, min(1.0_WP, (pos(3) - (this%amr%zlo + (real(kp,WP)+0.5_WP)*dz)) * dzi + 0.5_WP))
         wx2 = 1.0_WP - wx1; wy2 = 1.0_WP - wy1; wz2 = 1.0_WP - wz1
         
         if (is_staggered) then
            ! Face-centered indices and weights for each component
            ipu = max(ilo_, min(ihi_, floor((pos(1) - this%amr%xlo) * dxi)))
            jpv = max(jlo_, min(jhi_, floor((pos(2) - this%amr%ylo) * dyi)))
            kpw = max(klo_, min(khi_, floor((pos(3) - this%amr%zlo) * dzi)))
            wxu1 = max(0.0_WP, min(1.0_WP, (pos(1) - (this%amr%xlo + real(ipu,WP)*dx)) * dxi))
            wyv1 = max(0.0_WP, min(1.0_WP, (pos(2) - (this%amr%ylo + real(jpv,WP)*dy)) * dyi))
            wzw1 = max(0.0_WP, min(1.0_WP, (pos(3) - (this%amr%zlo + real(kpw,WP)*dz)) * dzi))
            wxu2 = 1.0_WP - wxu1; wyv2 = 1.0_WP - wyv1; wzw2 = 1.0_WP - wzw1
            ! U at x-faces
            vel(1) = wz1*(wy1*(wxu1*pU(ipu+1,jp+1,kp+1,1)+wxu2*pU(ipu,jp+1,kp+1,1)) + &
            &             wy2*(wxu1*pU(ipu+1,jp  ,kp+1,1)+wxu2*pU(ipu,jp  ,kp+1,1))) + &
            &        wz2*(wy1*(wxu1*pU(ipu+1,jp+1,kp  ,1)+wxu2*pU(ipu,jp+1,kp  ,1)) + &
            &             wy2*(wxu1*pU(ipu+1,jp  ,kp  ,1)+wxu2*pU(ipu,jp  ,kp  ,1)))
            ! V at y-faces
            vel(2) = wz1*(wyv1*(wx1*pV(ip+1,jpv+1,kp+1,1)+wx2*pV(ip,jpv+1,kp+1,1)) + &
            &             wyv2*(wx1*pV(ip+1,jpv  ,kp+1,1)+wx2*pV(ip,jpv  ,kp+1,1))) + &
            &        wz2*(wyv1*(wx1*pV(ip+1,jpv+1,kp  ,1)+wx2*pV(ip,jpv+1,kp  ,1)) + &
            &             wyv2*(wx1*pV(ip+1,jpv  ,kp  ,1)+wx2*pV(ip,jpv  ,kp  ,1)))
            ! W at z-faces
            vel(3) = wzw1*(wy1*(wx1*pW(ip+1,jp+1,kpw+1,1)+wx2*pW(ip,jp+1,kpw+1,1)) + &
            &              wy2*(wx1*pW(ip+1,jp  ,kpw+1,1)+wx2*pW(ip,jp  ,kpw+1,1))) + &
            &        wzw2*(wy1*(wx1*pW(ip+1,jp+1,kpw  ,1)+wx2*pW(ip,jp+1,kpw  ,1)) + &
            &              wy2*(wx1*pW(ip+1,jp  ,kpw  ,1)+wx2*pW(ip,jp  ,kpw  ,1)))
         else
            ! All cell-centered
            vel(1) = wz1*(wy1*(wx1*pU(ip+1,jp+1,kp+1,1)+wx2*pU(ip,jp+1,kp+1,1)) + &
            &             wy2*(wx1*pU(ip+1,jp  ,kp+1,1)+wx2*pU(ip,jp  ,kp+1,1))) + &
            &        wz2*(wy1*(wx1*pU(ip+1,jp+1,kp  ,1)+wx2*pU(ip,jp+1,kp  ,1)) + &
            &             wy2*(wx1*pU(ip+1,jp  ,kp  ,1)+wx2*pU(ip,jp  ,kp  ,1)))
            vel(2) = wz1*(wy1*(wx1*pV(ip+1,jp+1,kp+1,1)+wx2*pV(ip,jp+1,kp+1,1)) + &
            &             wy2*(wx1*pV(ip+1,jp  ,kp+1,1)+wx2*pV(ip,jp  ,kp+1,1))) + &
            &        wz2*(wy1*(wx1*pV(ip+1,jp+1,kp  ,1)+wx2*pV(ip,jp+1,kp  ,1)) + &
            &             wy2*(wx1*pV(ip+1,jp  ,kp  ,1)+wx2*pV(ip,jp  ,kp  ,1)))
            vel(3) = wz1*(wy1*(wx1*pW(ip+1,jp+1,kp+1,1)+wx2*pW(ip,jp+1,kp+1,1)) + &
            &             wy2*(wx1*pW(ip+1,jp  ,kp+1,1)+wx2*pW(ip,jp  ,kp+1,1))) + &
            &        wz2*(wy1*(wx1*pW(ip+1,jp+1,kp  ,1)+wx2*pW(ip,jp+1,kp  ,1)) + &
            &             wy2*(wx1*pW(ip+1,jp  ,kp  ,1)+wx2*pW(ip,jp  ,kp  ,1)))
         end if
      end function interp_velocity

   end subroutine advance_vof


   ! ============================================================================
   ! DEFERRED METHODS
   ! ============================================================================

   !> Get solver information
   subroutine get_info(this)
      use parallel, only: MPI_REAL_WP
      use mpi_f08
      class(amrvof), intent(inout) :: this
      integer :: lvl, ierr

      ! Initialize
      this%VFmin = huge(1.0_WP)
      this%VFmax = -huge(1.0_WP)
      this%VFint = 0.0_WP

      ! Loop over levels
       do lvl = 0, this%amr%clvl()
         this%VFmin = min(this%VFmin, this%VF%get_min(lvl=lvl))
         this%VFmax = max(this%VFmax, this%VF%get_max(lvl=lvl))
      end do

      ! Compute volume integral at level 0
      this%VFint = this%VF%get_sum(lvl=0) * this%amr%dx(0) * this%amr%dy(0) * this%amr%dz(0)

      ! Reduce across MPI ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%VFint, 1, MPI_REAL_WP, MPI_SUM, this%amr%comm, ierr)
   end subroutine get_info

   !> Register checkpoint
   subroutine register_checkpoint(this, io)
      use amrio_class, only: amrio
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%VF, 'VF')
      call io%add_data(this%Cliq, 'Cliq')
      call io%add_data(this%Cgas, 'Cgas')
      call io%add_data(this%PLIC, 'PLIC')
   end subroutine register_checkpoint

   !> Restore checkpoint
   subroutine restore_checkpoint(this, io, dirname)
      use amrio_class, only: amrio
      class(amrvof), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      call io%read_data(dirname, this%VF, 'VF')
      call io%read_data(dirname, this%Cliq, 'Cliq')
      call io%read_data(dirname, this%Cgas, 'Cgas')
      call io%read_data(dirname, this%PLIC, 'PLIC')
   end subroutine restore_checkpoint

   !> Compute advective CFL at finest level
   !> Takes external staggered velocity MultiFabs and returns max CFL
   subroutine get_cfl(this, U, V, W, dt, cfl)
      class(amrvof), intent(inout) :: this
      type(amrex_multifab), intent(in) :: U, V, W  !< Staggered velocity at clvl
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      real(WP) :: Umax, Vmax, Wmax, CFLx, CFLy, CFLz
      integer :: lvl
      ! Get finest level metrics
      lvl = this%amr%clvl()
      ! Get max velocity norms
      Umax = U%norm0()
      Vmax = V%norm0()
      Wmax = W%norm0()
      ! Compute directional CFLs
      CFLx = dt * Umax / this%amr%dx(lvl)
      CFLy = dt * Vmax / this%amr%dy(lvl)
      CFLz = dt * Wmax / this%amr%dz(lvl)
      ! Return max CFL
      cfl = max(CFLx, CFLy, CFLz)
   end subroutine get_cfl

   !> Print solver info to screen
   subroutine amrvof_print(this)
      use messager, only: log
      use string, only: str_long
      class(amrvof), intent(in) :: this
      character(len=str_long) :: message
      call log("VOF solver: "//trim(this%name))
      write(message,'("  VF range: [",ES12.5,", ",ES12.5,"]")') VFlo, VFhi
      call log(trim(message))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrvof_print

end module amrvof_class

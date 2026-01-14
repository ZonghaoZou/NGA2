!> Data class wrapping AMReX MultiFab
!> Designed to be managed by amrgrid Registry
module amrdata_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrex_amr_module, only: amrex_multifab,amrex_boxarray,amrex_distromap,&
   &                           amrex_multifab_build,amrex_multifab_destroy,amrex_geometry,&
   &                           amrex_interp_cell_cons
   implicit none
   private

   public :: amrdata
   public :: amrdata_on_init,amrdata_on_coarse,amrdata_on_remake,amrdata_on_clear

   !> Abstract interface for on_init callback
   !> User provides: subroutine my_init(lvl, mf, geom)
   abstract interface
      subroutine on_init_iface(lvl, mf, geom)
         import :: amrex_multifab, amrex_geometry
         integer, intent(in) :: lvl
         type(amrex_multifab), intent(inout) :: mf
         type(amrex_geometry), intent(in) :: geom
      end subroutine on_init_iface
   end interface

   !> Abstract interface for on_coarse callback (new fine level from coarse interpolation)
   abstract interface
      subroutine on_coarse_iface(lvl, mf, time)
         import :: amrex_multifab, WP
         integer, intent(in) :: lvl
         type(amrex_multifab), intent(inout) :: mf
         real(WP), intent(in) :: time
      end subroutine on_coarse_iface
   end interface

   !> Abstract interface for on_remake callback (regrid existing level)
   abstract interface
      subroutine on_remake_iface(lvl, mf_old, mf_new, time)
         import :: amrex_multifab, WP
         integer, intent(in) :: lvl
         type(amrex_multifab), intent(in) :: mf_old
         type(amrex_multifab), intent(inout) :: mf_new
         real(WP), intent(in) :: time
      end subroutine on_remake_iface
   end interface

   !> Abstract interface for on_clear callback
   abstract interface
      subroutine on_clear_iface(lvl, mf)
         import :: amrex_multifab
         integer, intent(in) :: lvl
         type(amrex_multifab), intent(inout) :: mf
      end subroutine on_clear_iface
   end interface

   !> Abstract interface for fillbc callback (physical BCs)
   abstract interface
      subroutine fillbc_iface(lvl, mf, geom, time, scomp, ncomp)
         import :: amrex_multifab, amrex_geometry, WP
         integer, intent(in) :: lvl
         type(amrex_multifab), intent(inout) :: mf
         type(amrex_geometry), intent(in) :: geom
         real(WP), intent(in) :: time
         integer, intent(in) :: scomp, ncomp
      end subroutine fillbc_iface
   end interface

   !> Data object wrapping a MultiFab hierarchy (one MultiFab per level)
   type :: amrdata
      ! The underlying AMReX objects (array of MultiFabs, one per level)
      type(amrex_multifab), dimension(:), allocatable :: mf
      ! Pointer to parent amrgrid (set by initialize)
      class(amrgrid), pointer :: amr => null()
      ! Metadata
      character(len=str_medium) :: name='UNNAMED_AMRDATA'
      integer :: ncomp=1                                   !< Number of components
      integer :: ng=0                                      !< Number of ghost cells
      integer :: interp=amrex_interp_cell_cons             !< Interpolation method
      integer, dimension(:,:), allocatable :: lo_bc,hi_bc  !< Boundary conditions: lo_bc(3,ncomp), hi_bc(3,ncomp)
      ! User-overridable callback pointers
      procedure(on_init_iface), pointer, nopass :: user_on_init => null()
      procedure(on_coarse_iface), pointer, nopass :: user_on_coarse => null()
      procedure(on_remake_iface), pointer, nopass :: user_on_remake => null()
      procedure(on_clear_iface), pointer, nopass :: user_on_clear => null()
      procedure(fillbc_iface), pointer, nopass :: user_fillbc => null()
   contains
      ! Core methods
      procedure :: initialize      !< Initialize amrdata with amrgrid and parameters
      procedure :: define
      procedure :: destroy
      procedure :: clear_level
      procedure :: fill_ghosts
      procedure :: fillbc          !< Default physical BC callback (uses amrex_filcc)
      procedure :: get_data_ptr
      procedure :: on_regrid
      procedure :: register        !< Register regrid callbacks with amrgrid
      ! Fill operations
      procedure :: fill            !< Fill ghost cells and coarse-fine boundary data
      procedure :: fill_mfab       !< Fill into a target MultiFab
      procedure :: fill_from_coarse !< Interpolate from coarse level only
      procedure :: average_down    !< Average from finest to coarsest level
      procedure :: average_downto  !< Average level lvl+1 down to level lvl
      ! Callback setters
      procedure :: set_on_init
      procedure :: set_on_coarse
      procedure :: set_on_remake
      procedure :: set_on_clear
      procedure :: set_fillbc
      ! Callback methods (handle user override internally)
      procedure :: on_init      !< Create level from scratch
      procedure :: on_coarse    !< Create level from coarse interpolation
      procedure :: on_remake    !< Remake existing level (regrid)
      procedure :: on_clear     !< Clear level data
   end type amrdata

   ! Public dispatch for BC callbacks
   public :: dispatch_fillbc

contains

   !> Define the data on a specific grid level (BoxArray + DistroMap)
   !> ncomp and ng should already be set via amrgrid%register
   subroutine define(this,lvl,ba,dm)
      use iso_c_binding, only: c_associated
      use messager, only: die
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Verify allocation (should be handled by amrgrid%register)
      if (.not.allocated(this%mf)) call die('[amrdata define] mf not allocated - use amrgrid%register')
      if (lvl.lt.0.or.lvl.gt.ubound(this%mf,1)) call die('[amrdata define] lvl out of bounds')
      ! Destroy if already allocated at this level
      if (c_associated(this%mf(lvl)%p)) call amrex_multifab_destroy(this%mf(lvl))
      ! Build the MultiFab
      call amrex_multifab_build(this%mf(lvl),ba,dm,this%ncomp,this%ng)
   end subroutine define

   !> Destroy data at a specific level
   subroutine clear_level(this,lvl)
      use iso_c_binding, only: c_associated
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      if (allocated(this%mf)) then
         if (lvl.ge.0.and.lvl.le.ubound(this%mf,1)) then
            if (c_associated(this%mf(lvl)%p)) call amrex_multifab_destroy(this%mf(lvl))
         end if
      end if
   end subroutine clear_level

   !> Destroy the underlying MultiFabs and BC arrays
   subroutine destroy(this)
      use iso_c_binding, only: c_associated
      class(amrdata), intent(inout) :: this
      integer :: i
      if (allocated(this%mf)) then
         do i=0,ubound(this%mf,1)
            if (c_associated(this%mf(i)%p)) call amrex_multifab_destroy(this%mf(i))
         end do
         deallocate(this%mf)
      end if
      if (allocated(this%lo_bc)) deallocate(this%lo_bc)
      if (allocated(this%hi_bc)) deallocate(this%hi_bc)
   end subroutine destroy

   !> Fill ghost cells (boundary conditions) for a specific level
   subroutine fill_ghosts(this,lvl,geom)
      use iso_c_binding, only: c_associated
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_geometry), intent(in) :: geom
      if (allocated(this%mf).and.lvl.ge.0.and.lvl.le.ubound(this%mf,1)) then
         if (c_associated(this%mf(lvl)%p)) call this%mf(lvl)%fill_boundary(geom)
      end if
   end subroutine fill_ghosts

   !> Get access to raw data pointer for a specific block and level
   subroutine get_data_ptr(this,lvl,mfi,ptr)
      use amrex_amr_module, only: amrex_mfiter
      use iso_c_binding, only: c_double
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      real(c_double), dimension(:,:,:,:), pointer, intent(out) :: ptr
      ptr=>this%mf(lvl)%dataPtr(mfi)
   end subroutine get_data_ptr

   !> Hook called by Registry after regridding
   subroutine on_regrid(this,lvl,time)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      ! Default: do nothing (override for fillpatch, etc.)
   end subroutine on_regrid

   !> Default physical BC callback - applies amrex_filcc using lo_bc/hi_bc
   !> This is called by amrgrid%fill() via the C++ FillPatch wrapper
   !> If user_fillbc is set, it is called instead
   subroutine fillbc(this,mf_ptr,geom_ptr,time,scomp,ncomp)
      use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build
      use iso_c_binding, only: c_ptr,c_int,c_f_pointer
      class(amrdata), intent(inout) :: this
      type(c_ptr), value :: mf_ptr
      type(c_ptr), value :: geom_ptr
      real(WP), value :: time
      integer, value :: scomp,ncomp
      type(amrex_geometry) :: geom
      type(amrex_multifab) :: mf
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer, dimension(4) :: plo,phi
      integer :: lvl
      ! Convert C pointers to Fortran types
      geom=geom_ptr
      mf=mf_ptr
      ! Determine level from geometry (lvl=0 for first, etc.)
      lvl = -1  ! We don't have level info in this callback
      ! If user override is set, call it and return
      if (associated(this%user_fillbc)) then
         call this%user_fillbc(lvl, mf, geom, time, scomp, ncomp)
         return
      end if
      ! Default: Loop over boxes and apply filcc
      call amrex_mfiter_build(mfi,mf)
      do while(mfi%next())
         p=>mf%dataptr(mfi)
         ! Check if part of box is outside the domain
         if (.not.geom%domain%contains(p)) then
            plo=lbound(p); phi=ubound(p)
            call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,&
            &   geom%get_physical_location(plo),this%lo_bc,this%hi_bc)
         end if
      end do
   end subroutine fillbc

   !> Set user-provided callback for on_init
   subroutine set_on_init(this, func)
      class(amrdata), intent(inout) :: this
      procedure(on_init_iface) :: func
      this%user_on_init => func
   end subroutine set_on_init

   !> Set user-provided callback for on_coarse
   subroutine set_on_coarse(this, func)
      class(amrdata), intent(inout) :: this
      procedure(on_coarse_iface) :: func
      this%user_on_coarse => func
   end subroutine set_on_coarse

   !> Set user-provided callback for on_remake
   subroutine set_on_remake(this, func)
      class(amrdata), intent(inout) :: this
      procedure(on_remake_iface) :: func
      this%user_on_remake => func
   end subroutine set_on_remake

   !> Set user-provided callback for on_clear
   subroutine set_on_clear(this, func)
      class(amrdata), intent(inout) :: this
      procedure(on_clear_iface) :: func
      this%user_on_clear => func
   end subroutine set_on_clear

   !> Set user-provided callback for fillbc (physical BCs)
   subroutine set_fillbc(this, func)
      class(amrdata), intent(inout) :: this
      procedure(fillbc_iface) :: func
      this%user_fillbc => func
   end subroutine set_fillbc

   !> Default on_init callback: define MultiFab, then call user override if set
   subroutine on_init(this, lvl, ba, dm, geom)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_geometry), intent(in) :: geom
      ! Build the MultiFab
      call this%define(lvl, ba, dm)
      ! Call user override if set
      if (associated(this%user_on_init)) then
         call this%user_on_init(lvl, this%mf(lvl), geom)
      end if
   end subroutine on_init


   !> Default on_coarse callback: define MultiFab, then fill from coarse (or user override)
   subroutine on_coarse(this, lvl, ba, dm, time)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      real(WP), intent(in) :: time
      ! Build the MultiFab
      call this%define(lvl, ba, dm)
      ! Call user override if set, otherwise fill from coarse
      if (associated(this%user_on_coarse)) then
         call this%user_on_coarse(lvl, this%mf(lvl), time)
      else
         call this%fill_from_coarse(lvl, time)
      end if
   end subroutine on_coarse


   !> Default on_remake callback: build new MultiFab, fill from old+coarse, copy back
   subroutine on_remake(this, lvl, ba, dm, time)
      use amrex_amr_module, only: amrex_multifab_build, amrex_multifab_destroy
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      real(WP), intent(in) :: time
      type(amrex_multifab) :: mf_new
      ! Build temp MultiFab with new geometry
      call amrex_multifab_build(mf_new, ba, dm, this%ncomp, 0)
      ! Fill temp from old data (handles C-F interface)
      call this%fill_mfab(mf_new, lvl, time)
      ! Destroy old
      call this%clear_level(lvl)
      ! Build new with proper ghost count
      call amrex_multifab_build(this%mf(lvl), ba, dm, this%ncomp, this%ng)
      ! Copy from temp
      call this%mf(lvl)%copy(mf_new, 1, 1, this%ncomp, 0)
      ! Destroy temp
      call amrex_multifab_destroy(mf_new)
      ! Call user override if set
      if (associated(this%user_on_remake)) then
         call this%user_on_remake(lvl, mf_new, this%mf(lvl), time)
      end if
   end subroutine on_remake


   !> Default on_clear callback: clear level data
   subroutine on_clear(this, lvl)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Call user override if set (before clearing)
      if (associated(this%user_on_clear)) then
         call this%user_on_clear(lvl, this%mf(lvl))
      end if
      ! Always clear the level
      call this%clear_level(lvl)
   end subroutine on_clear

   !> Initialize amrdata with amrgrid reference and parameters
   !> This sets the amr pointer, allocates mf array, and configures BCs
   !> Required before calling fill, fill_from_coarse, fill_mfab
   subroutine initialize(this, amr, name, ncomp, ng)
      use amrex_amr_module, only: amrex_bc_int_dir
      implicit none
      class(amrdata), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in) :: name
      integer, intent(in) :: ncomp
      integer, intent(in), optional :: ng
      ! Store pointer to amrgrid
      this%amr => amr
      ! Set metadata
      this%name = name
      this%ncomp = ncomp
      if (present(ng)) this%ng = ng
      ! Allocate mf array
      if (.not.allocated(this%mf)) allocate(this%mf(0:amr%maxlvl))
      ! Allocate and set default BCs (periodic if grid is periodic)
      if (.not.allocated(this%lo_bc)) allocate(this%lo_bc(3,ncomp))
      if (.not.allocated(this%hi_bc)) allocate(this%hi_bc(3,ncomp))
      this%lo_bc = 0  ! Default: reflect_even
      this%hi_bc = 0
      if (amr%xper) this%lo_bc(1,:) = amrex_bc_int_dir
      if (amr%yper) this%lo_bc(2,:) = amrex_bc_int_dir
      if (amr%zper) this%lo_bc(3,:) = amrex_bc_int_dir
      this%hi_bc = this%lo_bc
   end subroutine initialize

   !> Register regrid callbacks with amrgrid (optional, for auto-regrid handling)
   !> Requires initialize was called first
   subroutine register(this)
      use iso_c_binding, only: c_loc
      use messager, only: die
      implicit none
      class(amrdata), target, intent(inout) :: this
      ! Verify amr is set
      if (.not.associated(this%amr)) call die('[amrdata%register] Must call initialize first')
      ! Register all 4 callbacks with amrgrid
      ! Note: Need select type since c_loc requires non-polymorphic type
      select type (this)
       type is (amrdata)
         call this%amr%add_on_init(amrdata_on_init, c_loc(this))
         call this%amr%add_on_coarse(amrdata_on_coarse, c_loc(this))
         call this%amr%add_on_remake(amrdata_on_remake, c_loc(this))
         call this%amr%add_on_clear(amrdata_on_clear, c_loc(this))
      end select
   end subroutine register


   !> Dispatch callback for on_init (matches level_callback interface)
   subroutine amrdata_on_init(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(lvl, ba, dm, this%amr%geom(lvl))
   end subroutine amrdata_on_init


   !> Dispatch callback for on_coarse (matches level_callback interface)
   subroutine amrdata_on_coarse(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      ! Define MultiFab first
      call this%define(lvl, ba, dm)
      ! Call user override if set, otherwise use default fill_from_coarse
      if (associated(this%user_on_coarse)) then
         call this%user_on_coarse(lvl, this%mf(lvl), time)
      else
         call this%fill_from_coarse(lvl, time)
      end if
   end subroutine amrdata_on_coarse


   !> Dispatch callback for on_remake (matches level_callback interface)
   subroutine amrdata_on_remake(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr,c_f_pointer,c_associated
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      type(amrex_multifab) :: newmf
      call c_f_pointer(ctx, this)
      ! Build temp MultiFab with new geometry
      call amrex_multifab_build(newmf, ba, dm, this%ncomp, this%ng)
      ! Call user override if set, otherwise use default fill_mfab
      if (associated(this%user_on_remake)) then
         call this%user_on_remake(lvl, this%mf(lvl), newmf, time)
      else
         call this%fill_mfab(newmf, lvl, time)
      end if
      ! Destroy old MultiFab
      if (c_associated(this%mf(lvl)%p)) call amrex_multifab_destroy(this%mf(lvl))
      ! Assign new as current
      this%mf(lvl) = newmf
   end subroutine amrdata_on_remake


   !> Dispatch callback for on_clear (matches clear_callback interface)
   subroutine amrdata_on_clear(ctx, lvl)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      ! Call user override if set, otherwise use default clear_level
      if (associated(this%user_on_clear)) then
         call this%user_on_clear(lvl, this%mf(lvl))
      else
         call this%clear_level(lvl)
      end if
   end subroutine amrdata_on_clear


   !> Fill ghost cells and coarse-fine boundary data
   !> For lvl=0: fills ghost cells using boundary conditions
   !> For lvl>0: interpolates from coarser level and fills ghosts
   subroutine fill(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two
      implicit none
      class(amrdata), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      real(WP) :: t_old, t_new
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      t_old = time - 1.0e200_WP
      t_new = time
      ! Get context pointer (this amrdata)
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(dispatch_fillbc)
      ! Call appropriate FillPatch (scomp/dcomp use 1-indexed Fortran convention)
      if (lvl .eq. 0) then
         call amrmfab_fillpatch_single(this%mf(0)%p, t_old, this%mf(0)%p, &
         &   t_new, this%mf(0)%p, this%amr%geom(0)%p, data_ctx, bc_dispatch_ptr, &
         &   time, 1, 1, this%ncomp)
      else
         call amrmfab_fillpatch_two(this%mf(lvl)%p, t_old, this%mf(lvl-1)%p, &
         &   t_new, this%mf(lvl-1)%p, this%amr%geom(lvl-1)%p, &
         &   t_old, this%mf(lvl)%p, t_new, this%mf(lvl)%p, this%amr%geom(lvl)%p, &
         &   data_ctx, bc_dispatch_ptr, time, 1, 1, this%ncomp, &
         &   this%amr%rref(lvl-1), this%interp, this%lo_bc, this%hi_bc, this%ncomp)
      end if
   end subroutine fill


   !> Fill into a target MultiFab from this amrdata (for regrid callbacks)
   subroutine fill_mfab(this, dest, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_amr_module, only: amrex_multifab
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two
      implicit none
      class(amrdata), target, intent(inout) :: this
      type(amrex_multifab), intent(inout) :: dest
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      real(WP) :: t_old, t_new
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      t_old = time - 1.0e200_WP
      t_new = time
      ! Get context pointer
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(dispatch_fillbc)
      ! Call appropriate FillPatch
      if (lvl .eq. 0) then
         call amrmfab_fillpatch_single(dest%p, t_old, this%mf(0)%p, &
         &   t_new, this%mf(0)%p, this%amr%geom(0)%p, data_ctx, bc_dispatch_ptr, &
         &   time, 1, 1, this%ncomp)
      else
         call amrmfab_fillpatch_two(dest%p, t_old, this%mf(lvl-1)%p, &
         &   t_new, this%mf(lvl-1)%p, this%amr%geom(lvl-1)%p, &
         &   t_old, this%mf(lvl)%p, t_new, this%mf(lvl)%p, this%amr%geom(lvl)%p, &
         &   data_ctx, bc_dispatch_ptr, time, 1, 1, this%ncomp, &
         &   this%amr%rref(lvl-1), this%interp, this%lo_bc, this%hi_bc, this%ncomp)
      end if
   end subroutine fill_mfab


   !> Fill fine level from coarse only (for creating new fine levels)
   subroutine fill_from_coarse(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillcoarsepatch
      implicit none
      class(amrdata), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      ! Only valid for fine levels
      if (lvl .lt. 1) return
      ! Get context pointer
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(dispatch_fillbc)
      ! Call C++ wrapper
      call amrmfab_fillcoarsepatch(this%mf(lvl)%p, time, this%mf(lvl-1)%p, &
      &   this%amr%geom(lvl-1)%p, this%amr%geom(lvl)%p, data_ctx, bc_dispatch_ptr, &
      &   1, 1, this%ncomp, this%amr%rref(lvl-1), this%interp, this%lo_bc, this%hi_bc, this%ncomp)
   end subroutine fill_from_coarse


   !> Average down from finest level to coarsest (ensures level consistency)
   subroutine average_down(this)
      use amrex_amr_module, only: amrex_average_down
      implicit none
      class(amrdata), intent(inout) :: this
      integer :: lvl
      if (.not.associated(this%amr)) return
      ! Loop from finest to coarsest
      do lvl = this%amr%clvl()-1, 0, -1
         call amrex_average_down(this%mf(lvl+1), this%mf(lvl), &
         &   this%amr%geom(lvl+1), this%amr%geom(lvl), 1, this%ncomp, this%amr%rref(lvl))
      end do
   end subroutine average_down


   !> Average level lvl+1 down to level lvl
   subroutine average_downto(this, lvl)
      use amrex_amr_module, only: amrex_average_down
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      if (.not.associated(this%amr)) return
      if (lvl.lt.0.or.lvl.ge.this%amr%clvl()) return
      call amrex_average_down(this%mf(lvl+1), this%mf(lvl), &
      &   this%amr%geom(lvl+1), this%amr%geom(lvl), 1, this%ncomp, this%amr%rref(lvl))
   end subroutine average_downto

   !> Dispatch FillPatch BC callback to amrdata%fillbc
   subroutine dispatch_fillbc(data_ctx, mf_ptr, geom_ptr, time, scomp, ncomp) bind(c)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_int
      implicit none
      type(c_ptr), value, intent(in) :: data_ctx
      type(c_ptr), value, intent(in) :: mf_ptr
      type(c_ptr), value, intent(in) :: geom_ptr
      real(c_double), value, intent(in) :: time
      integer(c_int), value, intent(in) :: scomp
      integer(c_int), value, intent(in) :: ncomp
      type(amrdata), pointer :: data
      call c_f_pointer(data_ctx, data)
      call data%fillbc(mf_ptr, geom_ptr, real(time, WP), int(scomp), int(ncomp))
   end subroutine dispatch_fillbc


end module amrdata_class

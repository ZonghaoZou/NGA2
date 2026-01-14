!> AMR config object is defined based on amrcore AMREX object
!> Amrconfig differs quite a bit from other configs
module amrgrid_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_char,c_int,c_funptr,c_funloc,c_loc
   use precision,        only: WP
   use string,           only: str_medium
   use amrex_amr_module, only: amrex_geometry,amrex_boxarray,amrex_distromap
   use mpi_f08,          only: MPI_Comm
   implicit none
   private


   ! Expose type/constructor/methods
   public :: amrgrid

   !> Abstract interface for user-provided tagging callback
   abstract interface
      subroutine tagging_callback(lvl,tags,time)
         use iso_c_binding, only: c_ptr
         use precision, only: WP
         implicit none
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags  !< amrex_tagboxarray C pointer
         real(WP), intent(in) :: time
      end subroutine tagging_callback
   end interface

   !> Abstract interface for post-regrid callback (with context)
   abstract interface
      subroutine postregrid_callback(ctx)
         use iso_c_binding, only: c_ptr
         implicit none
         type(c_ptr), intent(in) :: ctx
      end subroutine postregrid_callback
   end interface

   !> Abstract interface for level callbacks (init, coarse, remake)
   abstract interface
      subroutine level_callback(ctx,lvl,time,ba,dm)
         use iso_c_binding, only: c_ptr
         use precision, only: WP
         use amrex_amr_module, only: amrex_boxarray,amrex_distromap
         implicit none
         type(c_ptr), intent(in) :: ctx
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine level_callback
   end interface

   !> Abstract interface for clear level callback
   abstract interface
      subroutine clear_callback(ctx,lvl)
         use iso_c_binding, only: c_ptr
         implicit none
         type(c_ptr), intent(in) :: ctx
         integer, intent(in) :: lvl
      end subroutine clear_callback
   end interface

   ! Wrappers for callback lists
   type :: tagger_wrapper
      procedure(tagging_callback), pointer, nopass :: f=>null()
   end type tagger_wrapper
   type :: postregrid_wrapper
      procedure(postregrid_callback), pointer, nopass :: f=>null()
      type(c_ptr) :: ctx=c_null_ptr
   end type postregrid_wrapper
   type :: level_cb_wrapper
      procedure(level_callback), pointer, nopass :: f=>null()
      type(c_ptr) :: ctx=c_null_ptr
   end type level_cb_wrapper
   type :: clear_cb_wrapper
      procedure(clear_callback), pointer, nopass :: f=>null()
      type(c_ptr) :: ctx=c_null_ptr
   end type clear_cb_wrapper

   !> Amrgrid object definition based on AMReX's amrcore
   type :: amrgrid
      ! Name of amrgrid
      character(len=str_medium) :: name='UNNAMED_AMRGRID'
      ! Pointer to AMReX's amrcore object
      type(c_ptr) :: amrcore=c_null_ptr
      ! Self-pointer for C interop (set early before `this` becomes polymorphic)
      type(c_ptr), private :: self_ptr=c_null_ptr
      ! Coordinate system
      integer :: coordsys=0
      ! Domain periodicity
      logical :: xper,yper,zper
      ! Domain extent
      real(WP) :: xlo,ylo,zlo
      real(WP) :: xhi,yhi,zhi
      ! Level 0 mesh dimensions
      integer :: nx,ny,nz
      ! Number of refinement levels
      integer :: maxlvl
      ! Max grid size
      integer :: nmax=32
      ! Blocking factor
      integer :: nbloc=8
      ! Refinement ratio
      integer, dimension(:), allocatable :: rref
      ! Geometry object at each level
      type(amrex_geometry), dimension(:), allocatable :: geom
      ! Shortcut for domain volume
      real(WP) :: vol
      ! Shortcut to cell size per level
      real(WP), dimension(:), allocatable :: dx,dy,dz
      real(WP) :: min_meshsize
      ! Parallel info
      type(MPI_Comm) :: comm            !< Communicator for our group
      integer        :: nproc           !< Number of processors
      integer        :: rank            !< Processor grid rank
      logical        :: amRoot          !< Am I the root?
      ! Monitoring info - call get_info() to update
      integer  :: nlevels=-1            !< Current total number of levels
      integer  :: nboxes =-1            !< Current total number of boxes
      real(WP) :: ncells =-1.0_WP       !< Current total number of cells (real!)
      real(WP) :: compression=-1.0_WP   !< Current compression ratio (ncells/total cells with uniform mesh)
      ! Level callback lists (solvers register their handlers)
      type(level_cb_wrapper), dimension(:), allocatable :: on_init
      type(level_cb_wrapper), dimension(:), allocatable :: on_coarse
      type(level_cb_wrapper), dimension(:), allocatable :: on_remake
      type(clear_cb_wrapper), dimension(:), allocatable :: on_clear
      ! Tagging and post-regrid callback lists
      type(tagger_wrapper), dimension(:), allocatable :: taggers
      type(postregrid_wrapper), dimension(:), allocatable :: postregrid_funcs
   contains
      procedure :: initialize                !< Initialization of amrgrid object
      procedure :: finalize                  !< Finalization of amrgrid object
      procedure :: initialize_grid           !< Initialize data on armgrid according to registered function
      procedure :: regrid                    !< Perform regriding operation on level baselvl
      procedure :: get_info                  !< Calculate various information on our amrgrid object
      procedure :: print                     !< Print out grid info
      ! Level callbacks - solvers register their handlers
      procedure :: add_on_init               !< Add init level callback
      procedure :: add_on_coarse             !< Add coarse level callback
      procedure :: add_on_remake             !< Add remake level callback
      procedure :: add_on_clear              !< Add clear level callback
      procedure :: clear_level_callbacks     !< Clear all level callbacks

      procedure :: get_boxarray              !< Obtain box array at a given level
      procedure :: get_distromap             !< Obtain distromap at a given level

      procedure :: mfiter_build              !< Build mfiter at a given level
      procedure :: mfiter_destroy            !< Destroy mfiter
      procedure :: clvl                      !< Return current finest level
      procedure :: mfab_build                !< Build multifab at a given level
      procedure :: mfab_destroy              !< Destroy multifab
      procedure :: add_tagging               !< Add a tagging callback
      procedure :: add_postregrid            !< Add a post-regrid callback
      procedure :: clear_tagging             !< Clear all tagging callbacks
      procedure :: clear_postregrid          !< Clear all post-regrid callbacks
      procedure :: init_from_scratch         !< Build all levels by looping regrid
   end type amrgrid

   ! Instance counter for automated AMReX lifecycle management
   integer :: amr_instance_count=0

contains


   !> Initialization of an amrgrid object
   subroutine initialize(this,name)
      use messager, only: die
      implicit none
      class(amrgrid), target, intent(inout) :: this
      character(len=*), intent(in), optional :: name
      ! First of all, confirm that nga2 and amrex reals are compatible
      check_real: block
         use messager,         only: die
         use precision,        only: WP
         use amrex_amr_module, only: amrex_real
         if (amrex_real.ne.WP) call die('Incompatible real type between nga2 and amrex!')
      end block check_real
      ! Check if AMReX has been initialized - if not, initialize it here
      check_if_init: block
         use amrex_amr_module, only: amrex_initialized,amrex_init
         use parallel, only: comm
         ! Increment instance counter
         amr_instance_count=amr_instance_count+1
         ! Initialize AMReX if it's the first instance or not yet initialized
         if (.not.amrex_initialized()) call amrex_init(comm=comm%MPI_VAL,arg_parmparse=.false.)
      end block check_if_init
      ! Set parameters for amrcore and geometry
      set_params: block
         use amrex_amr_module, only: amrex_parmparse,amrex_parmparse_build,amrex_parmparse,amrex_parmparse_destroy
         type(amrex_parmparse) :: pp
         integer, dimension(3) :: per
         call amrex_parmparse_build(pp,'amr')
         call pp%addarr('n_cell'         ,[this%nx,this%ny,this%nz])
         if (this%maxlvl.lt.0) call die('[amrgrid initialize] maxlvl must be >= 0')
         call pp%add   ('max_level'      ,this%maxlvl)
         call pp%add   ('blocking_factor',this%nbloc)
         call pp%add   ('max_grid_size'  ,this%nmax)
         if (.not.allocated(this%rref)) this%rref=[2]
         call pp%addarr('ref_ratio'      ,this%rref)
         call amrex_parmparse_destroy(pp)
         call amrex_parmparse_build(pp,'geometry')
         call pp%add   ('coord_sys'      ,this%coordsys)
         per=0
         if (this%xper) per(1)=1
         if (this%yper) per(2)=1
         if (this%zper) per(3)=1
         call pp%addarr('is_periodic',per)
         call pp%addarr('prob_lo'    ,[this%xlo,this%ylo,this%zlo])
         call pp%addarr('prob_hi'    ,[this%xhi,this%yhi,this%zhi])
         call amrex_parmparse_destroy(pp)
      end block set_params
      ! Create an amrcore object using our C++ interface
      create_amrcore_obj: block
         use amrex_interface, only: amrcore_create,amrcore_set_owner, &
         &   amrcore_set_on_init_dispatch,amrcore_set_on_coarse_dispatch, &
         &   amrcore_set_on_remake_dispatch,amrcore_set_on_clear_dispatch, &
         &   amrcore_set_on_tag_dispatch,amrcore_set_on_postregrid_dispatch
         this%amrcore=amrcore_create()
         ! Set owner pointer - use select type to get c_loc of concrete type
         select type (this)
          type is (amrgrid)
            this%self_ptr=c_loc(this)
         end select
         call amrcore_set_owner(this%amrcore,this%self_ptr)
         call amrcore_set_on_init_dispatch(this%amrcore,c_funloc(dispatch_mak_lvl_init))
         call amrcore_set_on_coarse_dispatch(this%amrcore,c_funloc(dispatch_mak_lvl_crse))
         call amrcore_set_on_remake_dispatch(this%amrcore,c_funloc(dispatch_mak_lvl_remk))
         call amrcore_set_on_clear_dispatch(this%amrcore,c_funloc(dispatch_clr_lvl))
         call amrcore_set_on_tag_dispatch(this%amrcore,c_funloc(dispatch_err_est))
         call amrcore_set_on_postregrid_dispatch(this%amrcore,c_funloc(dispatch_postregrid))
         if (present(name)) this%name=trim(adjustl(name))
      end block create_amrcore_obj
      ! Get back geometry objects
      store_geometries: block
         use amrex_amr_module, only: amrex_geometry_init_data
         use amrex_interface, only: amrcore_get_geometry
         integer :: n
         allocate(this%geom(0:this%maxlvl))
         do n=0,this%maxlvl
            call amrcore_get_geometry(this%geom(n)%p,n,this%amrcore)
            call amrex_geometry_init_data(this%geom(n))
         end do
      end block store_geometries
      ! Store effective refinement ratio
      store_ref_ratio: block
         use amrex_interface, only: amrcore_get_ref_ratio
         if (allocated(this%rref)) deallocate(this%rref); allocate(this%rref(0:this%maxlvl-1))
         call amrcore_get_ref_ratio(this%rref,this%amrcore)
      end block store_ref_ratio
      ! Store parallel info
      store_parallel_info: block
         use parallel, only: comm,nproc,rank,amRoot
         this%comm=comm
         this%nproc=nproc
         this%rank=rank
         this%amRoot=amRoot
      end block store_parallel_info
      ! Compute shortcuts
      compute_shortcuts: block
         integer :: lvl
         ! Mesh size at each level
         allocate(this%dx(0:this%maxlvl),this%dy(0:this%maxlvl),this%dz(0:this%maxlvl))
         do lvl=0,this%maxlvl
            this%dx(lvl)=this%geom(lvl)%dx(1)
            this%dy(lvl)=this%geom(lvl)%dx(2)
            this%dz(lvl)=this%geom(lvl)%dx(3)
         end do
         ! Total domain volume
         this%vol=(this%xhi-this%xlo)*(this%yhi-this%ylo)*(this%zhi-this%zlo)
         ! Smallest mesh size
         this%min_meshsize=min(this%dx(this%maxlvl),this%dy(this%maxlvl),this%dz(this%maxlvl))
      end block compute_shortcuts
   end subroutine initialize


   !> Finalization of amrgrid object
   impure elemental subroutine finalize(this)
      use mpi_f08,       only: MPI_COMM_NULL
      use iso_c_binding, only: c_associated
      use amrex_interface, only: amrcore_destroy
      implicit none
      class(amrgrid), intent(inout) :: this
      ! Delete amrcore object if it exists
      if (c_associated(this%amrcore)) then
         call amrcore_destroy(this%amrcore)
         this%amrcore=c_null_ptr
      end if
      ! Deallocate allocatable arrays
      if (allocated(this%rref)) deallocate(this%rref)
      if (allocated(this%geom)) deallocate(this%geom)
      if (allocated(this%dx))   deallocate(this%dx)
      if (allocated(this%dy))   deallocate(this%dy)
      if (allocated(this%dz))   deallocate(this%dz)
      ! Deallocate callback lists
      if (allocated(this%on_init))   deallocate(this%on_init)
      if (allocated(this%on_coarse)) deallocate(this%on_coarse)
      if (allocated(this%on_remake)) deallocate(this%on_remake)
      if (allocated(this%on_clear))  deallocate(this%on_clear)
      if (allocated(this%taggers)) deallocate(this%taggers)
      if (allocated(this%postregrid_funcs)) deallocate(this%postregrid_funcs)
      ! Do not free comm as it was passed to us
      this%comm=MPI_COMM_NULL
      ! Handle automated AMReX finalization
      auto_finalize: block
         use amrex_amr_module, only: amrex_initialized,amrex_finalize
         amr_instance_count=amr_instance_count-1
         if (amr_instance_count.eq.0.and.amrex_initialized()) call amrex_finalize()
      end block auto_finalize
   end subroutine finalize


   !> Add a level init callback (called for MakeNewLevelFromScratch)
   subroutine add_on_init(this,callback,ctx)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrgrid), intent(inout) :: this
      procedure(level_callback) :: callback
      type(c_ptr), intent(in) :: ctx
      type(level_cb_wrapper), dimension(:), allocatable :: tmp
      integer :: n
      if (allocated(this%on_init)) then
         n=size(this%on_init)
         allocate(tmp(n+1))
         tmp(1:n)=this%on_init
         tmp(n+1)%f=>callback
         tmp(n+1)%ctx=ctx
         call move_alloc(tmp,this%on_init)
      else
         allocate(this%on_init(1))
         this%on_init(1)%f=>callback
         this%on_init(1)%ctx=ctx
      end if
   end subroutine add_on_init

   !> Add a coarse level callback (called for MakeNewLevelFromCoarse)
   subroutine add_on_coarse(this,callback,ctx)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrgrid), intent(inout) :: this
      procedure(level_callback) :: callback
      type(c_ptr), intent(in) :: ctx
      type(level_cb_wrapper), dimension(:), allocatable :: tmp
      integer :: n
      if (allocated(this%on_coarse)) then
         n=size(this%on_coarse)
         allocate(tmp(n+1))
         tmp(1:n)=this%on_coarse
         tmp(n+1)%f=>callback
         tmp(n+1)%ctx=ctx
         call move_alloc(tmp,this%on_coarse)
      else
         allocate(this%on_coarse(1))
         this%on_coarse(1)%f=>callback
         this%on_coarse(1)%ctx=ctx
      end if
   end subroutine add_on_coarse

   !> Add a remake level callback (called for RemakeLevel)
   subroutine add_on_remake(this,callback,ctx)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrgrid), intent(inout) :: this
      procedure(level_callback) :: callback
      type(c_ptr), intent(in) :: ctx
      type(level_cb_wrapper), dimension(:), allocatable :: tmp
      integer :: n
      if (allocated(this%on_remake)) then
         n=size(this%on_remake)
         allocate(tmp(n+1))
         tmp(1:n)=this%on_remake
         tmp(n+1)%f=>callback
         tmp(n+1)%ctx=ctx
         call move_alloc(tmp,this%on_remake)
      else
         allocate(this%on_remake(1))
         this%on_remake(1)%f=>callback
         this%on_remake(1)%ctx=ctx
      end if
   end subroutine add_on_remake

   !> Add a clear level callback (called for ClearLevel)
   subroutine add_on_clear(this,callback,ctx)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrgrid), intent(inout) :: this
      procedure(clear_callback) :: callback
      type(c_ptr), intent(in) :: ctx
      type(clear_cb_wrapper), dimension(:), allocatable :: tmp
      integer :: n
      if (allocated(this%on_clear)) then
         n=size(this%on_clear)
         allocate(tmp(n+1))
         tmp(1:n)=this%on_clear
         tmp(n+1)%f=>callback
         tmp(n+1)%ctx=ctx
         call move_alloc(tmp,this%on_clear)
      else
         allocate(this%on_clear(1))
         this%on_clear(1)%f=>callback
         this%on_clear(1)%ctx=ctx
      end if
   end subroutine add_on_clear

   !> Clear all level callbacks
   subroutine clear_level_callbacks(this)
      implicit none
      class(amrgrid), intent(inout) :: this
      if (allocated(this%on_init))   deallocate(this%on_init)
      if (allocated(this%on_coarse)) deallocate(this%on_coarse)
      if (allocated(this%on_remake)) deallocate(this%on_remake)
      if (allocated(this%on_clear))  deallocate(this%on_clear)
   end subroutine clear_level_callbacks

   !> Initialize grid on an amrgrid object
   subroutine initialize_grid(this,time)
      use amrex_interface, only: amrcore_init_from_scratch
      implicit none
      class(amrgrid), target, intent(inout) :: this
      real(WP), intent(in) :: time
      ! Generate grid and allocate data
      call amrcore_init_from_scratch(this%amrcore,time)
      ! Generate info about grid
      call this%get_info()
   end subroutine initialize_grid


   !> Perform regriding operation on baselvl
   subroutine regrid(this,baselvl,time)
      use amrex_interface, only: amrcore_regrid
      implicit none
      class(amrgrid), target, intent(inout) :: this
      integer,  intent(in) :: baselvl
      real(WP), intent(in) :: time
      ! Regenerate grid and resize/transfer data
      call amrcore_regrid(this%amrcore,baselvl,time)
      ! Generate info about grid
      call this%get_info()
   end subroutine regrid

   !> Get info on amrgrid object
   subroutine get_info(this)
      use amrex_amr_module, only: amrex_boxarray,amrex_box
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_boxarray) :: ba
      type(amrex_box)      :: bx
      integer :: n,m,nb
      integer, dimension(3) :: lo,hi
      ! Update number of levels
      this%nlevels=this%clvl()+1
      ! Update number of boxes and cells
      this%nboxes=0
      this%ncells=0.0_WP
      do n=0,this%clvl()
         ! Get boxarray
         ba=this%get_boxarray(lvl=n)
         ! Increment boxes
         nb=int(ba%nboxes())
         this%nboxes=this%nboxes+nb
         ! Traverse boxes and count cells
         do m=1,nb
            ! Get box
            bx=ba%get_box(m-1)
            ! Increment cells
            this%ncells=this%ncells+real((bx%hi(1)-bx%lo(1)+1),WP)*&
            &                       real((bx%hi(2)-bx%lo(2)+1),WP)*&
            &                       real((bx%hi(3)-bx%lo(3)+1),WP)
         end do
      end do
      ! Update compression ratio
      this%compression=this%geom(this%clvl())%dx(1)*&
      &                this%geom(this%clvl())%dx(2)*&
      &                this%geom(this%clvl())%dx(3)*&
      &                this%ncells/((this%xhi-this%xlo)*(this%yhi-this%ylo)*(this%zhi-this%zlo))
   end subroutine get_info


   !> Print amrgrid object
   subroutine print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use amrex_amr_module, only: amrex_boxarray,amrex_box
      use parallel, only: amRoot
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_boxarray) :: ba
      type(amrex_box)      :: bx
      integer :: n,m,nb
      integer, dimension(3) :: lo,hi
      if (this%amRoot) then
         write(output_unit,'("AMR Cartesian grid [",a,"]")') trim(this%name)
         write(output_unit,'(" > amr level = ",i2)') this%clvl()
         write(output_unit,'(" > max level = ",i2)') this%maxlvl
         if (allocated(this%rref)) write(output_unit,'(" > ref ratio = ",100(" ",i0))') this%rref
         write(output_unit,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%xlo,this%xhi,this%ylo,this%yhi,this%zlo,this%zhi
         write(output_unit,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper
         ! Loop over levels
         do n=0,this%clvl()
            ! Get boxarray at that level
            ba=this%get_boxarray(lvl=n)
            write(output_unit,'(" >>> Level ",i2,"/",i2," with [dx =",es12.5,",dy =",es12.5,",dz =",es12.5,"]")') n,this%clvl(),this%geom(n)%dx(1),this%geom(n)%dx(2),this%geom(n)%dx(3)
            ! Loop over boxes
            nb=int(ba%nboxes())
            do m=1,nb
               bx=ba%get_box(m-1)
               write(output_unit,'(" >>> Box ",i0,"/",i0," from [",i3,",",i3,",",i3,"] to [",i3,",",i3,",",i3,"]")') m,nb,bx%lo,bx%hi
            end do
         end do
      end if
   end subroutine print


   !> Obtain box array at a level
   function get_boxarray(this,lvl) result(ba)
      use amrex_amr_module, only: amrex_boxarray
      use amrex_interface, only: amrcore_get_boxarray
      implicit none
      class(amrgrid), intent(inout) :: this
      integer, intent(in)  :: lvl
      type(amrex_boxarray) :: ba
      call amrcore_get_boxarray(ba%p,lvl,this%amrcore)
   end function get_boxarray


   !> Obtain distromap at a level
   function get_distromap(this,lvl) result(dm)
      use amrex_amr_module, only: amrex_distromap
      use amrex_interface, only: amrcore_get_distromap
      implicit none
      class(amrgrid), intent(inout) :: this
      integer, intent(in)  :: lvl
      type(amrex_distromap) :: dm
      call amrcore_get_distromap(dm%p,lvl,this%amrcore)
   end function get_distromap


   !> Build mfiter at a level
   subroutine mfiter_build(this,lvl,mfi)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_mfiter_build
      implicit none
      class(amrgrid), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(out) :: mfi
      type(amrex_boxarray) :: ba
      type(amrex_distromap) :: dm
      ba=this%get_boxarray (lvl)
      dm=this%get_distromap(lvl)
      call amrex_mfiter_build(mfi,ba,dm)
   end subroutine mfiter_build


   !> Destroy mfiter
   subroutine mfiter_destroy(this,mfi)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_destroy
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_mfiter), intent(inout) :: mfi
      call amrex_mfiter_destroy(mfi)
   end subroutine mfiter_destroy

   ! --- STATIC DISPATCHERS (bind(c)) ---
   ! These receive calls from C++ with the owner pointer
   subroutine dispatch_mak_lvl_init(owner,lvl,time,ba,dm) bind(c)
      use iso_c_binding, only: c_ptr,c_int,c_double,c_f_pointer
      implicit none
      type(c_ptr), value, intent(in) :: owner
      integer(c_int), value, intent(in) :: lvl
      real(c_double), value, intent(in) :: time
      type(c_ptr), value, intent(in) :: ba,dm
      type(amrex_boxarray) :: ba_obj
      type(amrex_distromap) :: dm_obj
      type(amrgrid), pointer :: this_grid
      integer :: i
      call c_f_pointer(owner,this_grid)
      ba_obj=ba
      dm_obj=dm
      if (allocated(this_grid%on_init)) then
         do i=1,size(this_grid%on_init)
            call this_grid%on_init(i)%f(this_grid%on_init(i)%ctx,int(lvl),real(time,WP),ba_obj,dm_obj)
         end do
      end if
   end subroutine dispatch_mak_lvl_init

   subroutine dispatch_mak_lvl_crse(owner,lvl,time,ba,dm) bind(c)
      use iso_c_binding, only: c_ptr,c_int,c_double,c_f_pointer
      implicit none
      type(c_ptr), value, intent(in) :: owner
      integer(c_int), value, intent(in) :: lvl
      real(c_double), value, intent(in) :: time
      type(c_ptr), value, intent(in) :: ba,dm
      type(amrex_boxarray) :: ba_obj
      type(amrex_distromap) :: dm_obj
      type(amrgrid), pointer :: this_grid
      integer :: i
      call c_f_pointer(owner,this_grid)
      ba_obj=ba
      dm_obj=dm
      if (allocated(this_grid%on_coarse)) then
         do i=1,size(this_grid%on_coarse)
            call this_grid%on_coarse(i)%f(this_grid%on_coarse(i)%ctx,int(lvl),real(time,WP),ba_obj,dm_obj)
         end do
      end if
   end subroutine dispatch_mak_lvl_crse

   subroutine dispatch_mak_lvl_remk(owner,lvl,time,ba,dm) bind(c)
      use iso_c_binding, only: c_ptr,c_int,c_double,c_f_pointer
      implicit none
      type(c_ptr), value, intent(in) :: owner
      integer(c_int), value, intent(in) :: lvl
      real(c_double), value, intent(in) :: time
      type(c_ptr), value, intent(in) :: ba,dm
      type(amrex_boxarray) :: ba_obj
      type(amrex_distromap) :: dm_obj
      type(amrgrid), pointer :: this_grid
      integer :: i
      call c_f_pointer(owner,this_grid)
      ba_obj=ba
      dm_obj=dm
      if (allocated(this_grid%on_remake)) then
         do i=1,size(this_grid%on_remake)
            call this_grid%on_remake(i)%f(this_grid%on_remake(i)%ctx,int(lvl),real(time,WP),ba_obj,dm_obj)
         end do
      end if
   end subroutine dispatch_mak_lvl_remk

   subroutine dispatch_clr_lvl(owner,lvl) bind(c)
      use iso_c_binding, only: c_ptr,c_int,c_f_pointer
      implicit none
      type(c_ptr), value, intent(in) :: owner
      integer(c_int), value, intent(in) :: lvl
      type(amrgrid), pointer :: this_grid
      integer :: i
      call c_f_pointer(owner,this_grid)
      if (allocated(this_grid%on_clear)) then
         do i=1,size(this_grid%on_clear)
            call this_grid%on_clear(i)%f(this_grid%on_clear(i)%ctx,int(lvl))
         end do
      end if
   end subroutine dispatch_clr_lvl



   !> Error estimation callback - calls all registered taggers
   subroutine dispatch_err_est(owner, lvl, tags, time, tagval, clrval) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int, c_double, c_f_pointer
      implicit none
      type(c_ptr), value, intent(in) :: owner
      integer(c_int), value, intent(in) :: lvl
      type(c_ptr), value, intent(in) :: tags
      real(c_double), value, intent(in) :: time
      character(kind=c_char), value, intent(in) :: tagval, clrval
      type(amrgrid), pointer :: this_grid
      integer :: i
      call c_f_pointer(owner,this_grid)
      ! Call all registered tagging callbacks
      if (allocated(this_grid%taggers)) then
         do i=1,size(this_grid%taggers)
            call this_grid%taggers(i)%f(int(lvl), tags, real(time, WP))
         end do
      end if
   end subroutine dispatch_err_est


   !> Post-regrid callback - calls all registered post-regrid callbacks
   subroutine dispatch_postregrid(owner) bind(c)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), value, intent(in) :: owner
      type(amrgrid), pointer :: this_grid
      integer :: i
      call c_f_pointer(owner,this_grid)
      ! Call all registered post-regrid callbacks with their context
      if (allocated(this_grid%postregrid_funcs)) then
         do i=1,size(this_grid%postregrid_funcs)
            call this_grid%postregrid_funcs(i)%f(this_grid%postregrid_funcs(i)%ctx)
         end do
      end if
   end subroutine dispatch_postregrid


   !> Add a tagging callback to the list
   subroutine add_tagging(this,tagger)
      implicit none
      class(amrgrid), intent(inout) :: this
      procedure(tagging_callback) :: tagger
      type(tagger_wrapper), dimension(:), allocatable :: tmp
      integer :: n
      if (allocated(this%taggers)) then
         n=size(this%taggers)
         allocate(tmp(n+1))
         tmp(1:n)=this%taggers
         tmp(n+1)%f=>tagger
         call move_alloc(tmp,this%taggers)
      else
         allocate(this%taggers(1))
         this%taggers(1)%f=>tagger
      end if
   end subroutine add_tagging


   !> Clear all tagging callbacks
   subroutine clear_tagging(this)
      implicit none
      class(amrgrid), intent(inout) :: this
      if (allocated(this%taggers)) deallocate(this%taggers)
   end subroutine clear_tagging


   !> Add a post-regrid callback to the list
   subroutine add_postregrid(this,callback,ctx)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrgrid), intent(inout) :: this
      procedure(postregrid_callback) :: callback
      type(c_ptr), intent(in) :: ctx
      type(postregrid_wrapper), dimension(:), allocatable :: tmp
      integer :: n
      if (allocated(this%postregrid_funcs)) then
         n=size(this%postregrid_funcs)
         allocate(tmp(n+1))
         tmp(1:n)=this%postregrid_funcs
         tmp(n+1)%f=>callback
         tmp(n+1)%ctx=ctx
         call move_alloc(tmp,this%postregrid_funcs)
      else
         allocate(this%postregrid_funcs(1))
         this%postregrid_funcs(1)%f=>callback
         this%postregrid_funcs(1)%ctx=ctx
      end if
   end subroutine add_postregrid


   !> Clear all post-regrid callbacks
   subroutine clear_postregrid(this)
      implicit none
      class(amrgrid), intent(inout) :: this
      if (allocated(this%postregrid_funcs)) deallocate(this%postregrid_funcs)
   end subroutine clear_postregrid

   !> Initialize grid hierarchy from scratch using AMReX's native init
   !> This calls MakeNewLevelFromScratch for ALL levels (triggers on_init callbacks)
   !> If do_postregrid=.true. (default), fires post-regrid callbacks after init
   subroutine init_from_scratch(this, time, do_postregrid)
      use iso_c_binding, only: c_loc
      use amrex_interface, only: amrcore_init_from_scratch
      implicit none
      class(amrgrid), target, intent(inout) :: this
      real(WP), intent(in) :: time
      logical, intent(in), optional :: do_postregrid
      logical :: fire_postregrid
      ! Default: fire post-regrid callbacks
      fire_postregrid = .true.
      if (present(do_postregrid)) fire_postregrid = do_postregrid
      ! Use AMReX's native init_from_scratch which properly builds all levels
      call amrcore_init_from_scratch(this%amrcore, time)
      ! Optionally fire post-regrid callbacks (for average_down, etc.)
      if (fire_postregrid) then
         select type (this)
          type is (amrgrid)
            call dispatch_postregrid(c_loc(this))
         end select
      end if
   end subroutine init_from_scratch


   !> Return current finest level
   function clvl(this) result(cl)
      use amrex_amr_module, only: amrex_boxarray
      use amrex_interface, only: amrcore_finest_level
      implicit none
      class(amrgrid), intent(inout) :: this
      integer :: cl
      cl=amrcore_finest_level(this%amrcore)
   end function clvl


   !> Build multifab at a level
   subroutine mfab_build(this,lvl,mfab,ncomp,nover,atface)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_multifab,amrex_multifab_build
      implicit none
      class(amrgrid), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: mfab
      integer, intent(in) :: ncomp, nover
      logical, intent(in), optional :: atface(3)
      type(amrex_boxarray) :: ba
      type(amrex_distromap) :: dm
      ba=this%get_boxarray(lvl)
      dm=this%get_distromap(lvl)
      call amrex_multifab_build(mfab,ba,dm,ncomp,nover,atface)
   end subroutine mfab_build

   !> Destroy multifab
   subroutine mfab_destroy(this,mfab)
      use amrex_amr_module, only: amrex_multifab,amrex_multifab_destroy
      implicit none
      class(amrgrid), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mfab
      call amrex_multifab_destroy(mfab)
   end subroutine mfab_destroy


   !> Finalization of amrex
   subroutine finalize_amrex()
      use amrex_amr_module, only: amrex_finalize
      implicit none
      call amrex_finalize()
   end subroutine finalize_amrex


end module amrgrid_class

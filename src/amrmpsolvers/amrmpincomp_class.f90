!> AMR Incompressible Multiphase solver class
!> Grown from amrmpcomp_class by removing compressibility and energy equations
module amrmpincomp_class
   use iso_c_binding,    only: c_ptr
   use precision,        only: WP
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use surfmesh_class,   only: surfmesh
   use amrex_amr_module, only: amrex_boxarray,amrex_distromap
   implicit none
   private

   ! VF clipping thresholds
   real(WP), parameter, public :: VFlo=1.0e-12_WP   !< Below this VF, liquid vanishes
   real(WP), parameter, public :: VFhi=1.0_WP-VFlo  !< Above this VF, gas vanishes
   real(WP), parameter, public :: vol_eps=1.0e-8_WP !< Volume epsilon for division by zero

   ! PLIC boundary condition types
   integer, parameter, public :: BC_LIQ    =1       !< All liquid in ghost
   integer, parameter, public :: BC_GAS    =2       !< All gas in ghost
   integer, parameter, public :: BC_REFLECT=3       !< Symmetry (mirror across boundary)
   integer, parameter, public :: BC_USER   =4       !< User-defined callback

   ! Expose type
   public :: amrmpincomp

   !> AMR Incompressible Multiphase solver type
   type, extends(amrsolver) :: amrmpincomp
      ! User-configurable callbacks
      procedure(incomp_init_iface), pointer, nopass :: user_init=>null()
      procedure(incomp_tagging_iface), pointer, nopass :: user_tagging=>null()
      procedure(incomp_bc_iface), pointer, nopass :: user_bc=>null()
      procedure(vof_bc_iface), pointer, nopass :: user_vof_bc=>null()

      ! PLIC boundary conditions (per face, only used if direction is non-periodic)
      integer :: vof_lo_bc(3)=BC_REFLECT
      integer :: vof_hi_bc(3)=BC_REFLECT

      ! Constant phasic densities
      real(WP) :: rhoL=1.0_WP         !< Liquid density
      real(WP) :: rhoG=1.0_WP         !< Gas density

      ! Collocated velocity
      type(amrdata) :: U,V,W          !< Current velocity
      type(amrdata) :: Uold,Vold,Wold !< Old velocity

      ! Pressure
      type(amrdata) :: P

      ! Volume fraction and barycenters
      type(amrdata) :: VF,VFold            !< Liquid volume fraction
      type(amrdata) :: Cliq,Cliqold        !< Liquid barycenter (3 components)
      type(amrdata) :: Cgas,Cgasold        !< Gas barycenter (3 components)

      ! PLIC interface (4 components: nx, ny, nz, d)
      type(amrdata) :: PLIC,PLICold

      ! Physical properties
      type(amrdata) :: visc                !< Dynamic viscosity

      ! CFL numbers
      real(WP) :: CFLc_x=0.0_WP,CFLc_y=0.0_WP,CFLc_z=0.0_WP  !< Convective
      real(WP) :: CFLv_x=0.0_WP,CFLv_y=0.0_WP,CFLv_z=0.0_WP  !< Viscous

      ! Monitoring quantities
      real(WP) :: Umax=0.0_WP,Vmax=0.0_WP,Wmax=0.0_WP
      real(WP) :: Pmax=0.0_WP
      real(WP) :: VFint=0.0_WP,VFmin=0.0_WP,VFmax=0.0_WP
      real(WP) :: rhoKint=0.0_WP

      ! Number of overlap cells (2 for WENO3 stencil)
      integer :: nover=2

      ! Cost for mixed cells in load balancing (1.0 = same as pure, higher = more expensive)
      real(WP) :: SLcost=100.0_WP

      ! Tagging parameter for VOF
      integer :: regrid_buffer=10  !< Number of cells to buffer around interface for tagging

      ! Surface mesh for visualization
      type(surfmesh) :: smesh

      ! Load distribution diagnostics
      real(WP) :: ncells_max=0.0_WP, ncells_min=0.0_WP  !< Finest-level cells per rank (max/min)
      real(WP) :: nmixed_max=0.0_WP, nmixed_min=0.0_WP  !< Mixed cells per rank (max/min)
      ! Per-rank timing
      real(WP) :: wt_plic=0.0_WP       !< Full build_plic
      real(WP) :: wt_plicnet=0.0_WP    !< PLICnet reconstruction loop
      real(WP) :: wt_polygon=0.0_WP    !< Polygon extraction loop
      real(WP) :: wt_visc=0.0_WP       !< Viscosity models
      ! Reduced timing
      real(WP) :: wtmax_plic   =0.0_WP, wtmin_plic   =0.0_WP
      real(WP) :: wtmax_plicnet=0.0_WP, wtmin_plicnet=0.0_WP
      real(WP) :: wtmax_polygon=0.0_WP, wtmin_polygon=0.0_WP
      real(WP) :: wtmax_visc   =0.0_WP, wtmin_visc   =0.0_WP

   contains
      procedure :: initialize
      procedure :: finalize
      ! Lifecycle callbacks
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      procedure :: tagging
      procedure :: get_cost
      ! VOF sync/fill/average utilities
      procedure :: fill_moments_lvl
      procedure :: sync_moments_lvl
      procedure :: sync_moments
      procedure :: fill_plic_lvl
      procedure :: sync_plic_lvl
      procedure :: sync_plic
      procedure :: vof_average_down
      ! Physics
      procedure :: get_cfl
      procedure :: build_plic
      procedure :: build_polygons
      procedure :: reset_moments
      procedure :: add_vreman
      ! Print
      procedure :: get_info
      procedure :: print=>amrmpincomp_print
      ! Checkpoint I/O (deferred from amrsolver)
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrmpincomp

   !> Abstract interface for user init callback
   abstract interface
      subroutine incomp_init_iface(solver,lvl,time,ba,dm)
         import :: amrmpincomp,WP,amrex_boxarray,amrex_distromap
         class(amrmpincomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine incomp_init_iface
   end interface

   !> Abstract interface for tagging callback
   abstract interface
      subroutine incomp_tagging_iface(solver,lvl,tags,time)
         import :: amrmpincomp,c_ptr,WP
         class(amrmpincomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         type(c_ptr), intent(in) :: tags
         real(WP), intent(in) :: time
      end subroutine incomp_tagging_iface
   end interface

   !> Abstract interface for user BC callback
   abstract interface
      subroutine incomp_bc_iface(solver,pU,bc_bx,face,time)
         use amrex_amr_module, only: amrex_box
         import :: amrmpincomp,WP
         class(amrmpincomp), intent(inout) :: solver
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU
         type(amrex_box), intent(in) :: bc_bx
         integer, intent(in) :: face
         real(WP), intent(in) :: time
      end subroutine incomp_bc_iface
   end interface

   !> Abstract interface for VOF boundary condition callback
   abstract interface
      subroutine vof_bc_iface(solver,bx,pVF,pCliq,pCgas,pPLIC,face,time,what)
         use amrex_amr_module, only: amrex_box
         import :: amrmpincomp,WP
         class(amrmpincomp), intent(inout) :: solver
         type(amrex_box), intent(in) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCliq,pCgas,pPLIC
         integer, intent(in) :: face      !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         real(WP), intent(in) :: time
         integer, intent(in) :: what      !< 1=PLIC, 2=moments
      end subroutine vof_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrmpincomp_on_init(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_init)) call this%user_init(this,lvl,time,ba,dm)
   end subroutine amrmpincomp_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrmpincomp_on_coarse(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrmpincomp_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrmpincomp_on_remake(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrmpincomp_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrmpincomp_on_clear(ctx,lvl)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrmpincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrmpincomp_on_clear

   !> Dispatch tagging: calls type-bound method then user callback
   subroutine amrmpincomp_tagging(ctx,lvl,tags,time)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrmpincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%tagging(lvl,tags,time)
      if (associated(this%user_tagging)) call this%user_tagging(this,lvl,tags,time)
   end subroutine amrmpincomp_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrmpincomp_postregrid(ctx,lbase,time)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrmpincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrmpincomp_postregrid

   !> Dispatch cost: calls type-bound method
   subroutine amrmpincomp_get_cost(ctx,lvl,nboxes,costs,ba)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrmpincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%get_cost(lvl,nboxes,costs,ba)
   end subroutine amrmpincomp_get_cost

   !> Internal fillbc for velocity - calls default_fillbc first, then user_bc for ext_dir faces
   subroutine velocity_fillbc(this,mf,scomp,ncomp,time,geom)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,&
      &                           amrex_box,amrex_geometry,amrex_multifab,amrex_bc_ext_dir
      use amrdata_class, only: default_fillbc
      implicit none
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp, ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bc_bx
      class(amrmpincomp), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer :: dlo(3),dhi(3)
      
      ! First apply default BC handling (foextrap,hoextrap,reflect,etc.)
      call default_fillbc(this,mf,scomp,ncomp,time,geom)
      
      ! Access parent solver
      select type (s=>this%parent)
       class is (amrmpincomp)
         solver=>s
      end select
      
      ! Check if user callback exists
      if (.not.associated(solver%user_bc)) return

      ! Get domain bounds
      dlo=geom%domain%lo
      dhi=geom%domain%hi

      ! Loop over FABs and apply user_bc for ext_dir faces
      call amrex_mfiter_build(mfi,mf,tiling=.false.)
      do while (mfi%next())
         p=>mf%dataptr(mfi)
         ilo=lbound(p,1); ihi=ubound(p,1)
         jlo=lbound(p,2); jhi=ubound(p,2)
         klo=lbound(p,3); khi=ubound(p,3)

         ! X-LOW (face=1)
         if (this%lo_bc(1,1).eq.amrex_bc_ext_dir .and. ilo.lt.dlo(1)) then
            bc_bx=amrex_box([ilo,jlo,klo],[dlo(1)-1,jhi,khi])
            call solver%user_bc(solver,p,bc_bx,1,time)
         end if
         ! X-HIGH (face=2)
         if (this%hi_bc(1,1).eq.amrex_bc_ext_dir .and. ihi.gt.dhi(1)) then
            bc_bx=amrex_box([dhi(1)+1,jlo,klo],[ihi,jhi,khi])
            call solver%user_bc(solver,p,bc_bx,2,time)
         end if
         ! Y-LOW (face=3)
         if (this%lo_bc(2,1).eq.amrex_bc_ext_dir .and. jlo.lt.dlo(2)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,dlo(2)-1,khi])
            call solver%user_bc(solver,p,bc_bx,3,time)
         end if
         ! Y-HIGH (face=4)
         if (this%hi_bc(2,1).eq.amrex_bc_ext_dir .and. jhi.gt.dhi(2)) then
            bc_bx=amrex_box([ilo,dhi(2)+1,klo],[ihi,jhi,khi])
            call solver%user_bc(solver,p,bc_bx,4,time)
         end if
         ! Z-LOW (face=5)
         if (this%lo_bc(3,1).eq.amrex_bc_ext_dir .and. klo.lt.dlo(3)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,jhi,dlo(3)-1])
            call solver%user_bc(solver,p,bc_bx,5,time)
         end if
         ! Z-HIGH (face=6)
         if (this%hi_bc(3,1).eq.amrex_bc_ext_dir .and. khi.gt.dhi(3)) then
            bc_bx=amrex_box([ilo,jlo,dhi(3)+1],[ihi,jhi,khi])
            call solver%user_bc(solver,p,bc_bx,6,time)
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine velocity_fillbc

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the incompressible multiphase solver
   subroutine initialize(this,amr,name)
      use iso_c_binding, only: c_loc
      use amrex_amr_module, only: amrex_bc_foextrap
      implicit none
      class(amrmpincomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name=trim(name)
      else
         this%name='UNNAMED_MPINCOMP'
      end if
      this%amr=>amr

      ! Initialize collocated velocity (cell-centered)
      call this%U%initialize(amr,name='U',ncomp=1,ng=this%nover); this%U%parent=>this
      call this%V%initialize(amr,name='V',ncomp=1,ng=this%nover); this%V%parent=>this
      call this%W%initialize(amr,name='W',ncomp=1,ng=this%nover); this%W%parent=>this
      call this%Uold%initialize(amr,name='Uold',ncomp=1,ng=this%nover); this%Uold%parent=>this
      call this%Vold%initialize(amr,name='Vold',ncomp=1,ng=this%nover); this%Vold%parent=>this
      call this%Wold%initialize(amr,name='Wold',ncomp=1,ng=this%nover); this%Wold%parent=>this

      ! Initialize pressure with Neumann BCs
      call this%P%initialize(amr,name='P',ncomp=1,ng=1); this%P%parent=>this
      if (.not.amr%xper) then; this%P%lo_bc(1,1)=amrex_bc_foextrap; this%P%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.amr%yper) then; this%P%lo_bc(2,1)=amrex_bc_foextrap; this%P%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.amr%zper) then; this%P%lo_bc(3,1)=amrex_bc_foextrap; this%P%hi_bc(3,1)=amrex_bc_foextrap; end if

      ! Initialize volume fraction and barycenters
      call this%VF%initialize     (amr,name='VF'     ,ncomp=1,ng=this%nover); this%VF%parent     =>this
      call this%VFold%initialize  (amr,name='VFold'  ,ncomp=1,ng=this%nover); this%VFold%parent  =>this
      call this%Cliq%initialize   (amr,name='Cliq'   ,ncomp=3,ng=this%nover); this%Cliq%parent   =>this
      call this%Cliqold%initialize(amr,name='Cliqold',ncomp=3,ng=this%nover); this%Cliqold%parent=>this
      call this%Cgas%initialize   (amr,name='Cgas'   ,ncomp=3,ng=this%nover); this%Cgas%parent   =>this
      call this%Cgasold%initialize(amr,name='Cgasold',ncomp=3,ng=this%nover); this%Cgasold%parent=>this

      ! Initialize PLIC (4 components: nx, ny, nz, d)
      call this%PLIC%initialize   (amr,name='PLIC'   ,ncomp=4,ng=this%nover); this%PLIC%parent   =>this
      call this%PLICold%initialize(amr,name='PLICold',ncomp=4,ng=this%nover); this%PLICold%parent=>this

      ! Initialize dynamic viscosity with Neumann BCs
      call this%visc%initialize(amr,name='visc',ncomp=1,ng=this%nover); this%visc%parent=>this
      if (.not.amr%xper) then; this%visc%lo_bc(1,1)=amrex_bc_foextrap; this%visc%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.amr%yper) then; this%visc%lo_bc(2,1)=amrex_bc_foextrap; this%visc%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.amr%zper) then; this%visc%lo_bc(3,1)=amrex_bc_foextrap; this%visc%hi_bc(3,1)=amrex_bc_foextrap; end if

      ! Set velocity fillbc callbacks to internal handler
      this%U%fillbc=>velocity_fillbc
      this%V%fillbc=>velocity_fillbc
      this%W%fillbc=>velocity_fillbc

      ! Register callbacks with amrgrid
      select type (this)
       type is (amrmpincomp)
         call this%amr%add_on_init   (amrmpincomp_on_init,   c_loc(this))
         call this%amr%add_on_coarse (amrmpincomp_on_coarse, c_loc(this))
         call this%amr%add_on_remake (amrmpincomp_on_remake, c_loc(this))
         call this%amr%add_on_clear  (amrmpincomp_on_clear,  c_loc(this))
         call this%amr%add_tagging   (amrmpincomp_tagging,   c_loc(this))
         call this%amr%add_postregrid(amrmpincomp_postregrid,c_loc(this))
         call this%amr%set_get_cost  (amrmpincomp_get_cost,  c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the incompressible multiphase solver
   subroutine finalize(this)
      implicit none
      class(amrmpincomp), intent(inout) :: this
      ! Velocity
      call this%U%finalize(); call this%V%finalize(); call this%W%finalize()
      call this%Uold%finalize(); call this%Vold%finalize(); call this%Wold%finalize()
      ! Pressure
      call this%P%finalize()
      ! VOF moments
      call this%VF%finalize(); call this%VFold%finalize()
      call this%Cliq%finalize(); call this%Cliqold%finalize()
      call this%Cgas%finalize(); call this%Cgasold%finalize()
      ! PLIC
      call this%PLIC%finalize(); call this%PLICold%finalize()
      ! Physical properties
      call this%visc%finalize()
      ! Nullify pointers
      nullify(this%amr); nullify(this%user_init); nullify(this%user_tagging)
      nullify(this%user_bc); nullify(this%user_vof_bc)
   end subroutine finalize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      implicit none
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset velocity
      call this%U%reset_level(lvl,ba,dm); call this%V%reset_level(lvl,ba,dm); call this%W%reset_level(lvl,ba,dm)
      call this%Uold%reset_level(lvl,ba,dm); call this%Vold%reset_level(lvl,ba,dm); call this%Wold%reset_level(lvl,ba,dm)
      ! Reset pressure
      call this%P%reset_level(lvl,ba,dm)
      ! Reset VOF moments
      call this%VF%reset_level(lvl,ba,dm); call this%VFold%reset_level(lvl,ba,dm)
      call this%Cliq%reset_level(lvl,ba,dm); call this%Cliqold%reset_level(lvl,ba,dm)
      call this%Cgas%reset_level(lvl,ba,dm); call this%Cgasold%reset_level(lvl,ba,dm)
      ! Reset PLIC
      call this%PLIC%reset_level(lvl,ba,dm); call this%PLICold%reset_level(lvl,ba,dm)
      ! Reset physical properties
      call this%visc%reset_level(lvl,ba,dm)
      ! Zero out everything
      call this%U%setval(val=0.0_WP,lvl=lvl); call this%V%setval(val=0.0_WP,lvl=lvl); call this%W%setval(val=0.0_WP,lvl=lvl)
      call this%Uold%setval(val=0.0_WP,lvl=lvl); call this%Vold%setval(val=0.0_WP,lvl=lvl); call this%Wold%setval(val=0.0_WP,lvl=lvl)
      call this%P%setval(val=0.0_WP,lvl=lvl)
      call this%VF%setval(val=0.0_WP,lvl=lvl); call this%VFold%setval(val=0.0_WP,lvl=lvl)
      call this%Cliq%setval(val=0.0_WP,lvl=lvl); call this%Cliqold%setval(val=0.0_WP,lvl=lvl)
      call this%Cgas%setval(val=0.0_WP,lvl=lvl); call this%Cgasold%setval(val=0.0_WP,lvl=lvl)
      call this%PLIC%setval(val=0.0_WP,lvl=lvl); call this%PLICold%setval(val=0.0_WP,lvl=lvl)
      call this%visc%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using conservative interpolation
   subroutine on_coarse(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      
      ! --- Velocity: interpolate from coarse ---
      call this%U%on_coarse(this%U,lvl,time,ba,dm)
      call this%V%on_coarse(this%V,lvl,time,ba,dm)
      call this%W%on_coarse(this%W,lvl,time,ba,dm)
      call this%Uold%reset_level(lvl,ba,dm)
      call this%Vold%reset_level(lvl,ba,dm)
      call this%Wold%reset_level(lvl,ba,dm)
      
      ! --- Pressure: interpolate from coarse ---
      call this%P%on_coarse(this%P,lvl,time,ba,dm)
      
      ! --- Physical properties (just reset, will be recomputed) ---
      call this%visc%reset_level(lvl,ba,dm)
      
      ! --- VOF moments: interpolate from coarse ---
      call this%VF%on_coarse(this%VF,lvl,time,ba,dm)
      call this%Cliq%on_coarse(this%Cliq,lvl,time,ba,dm)
      call this%Cgas%on_coarse(this%Cgas,lvl,time,ba,dm)
      call this%VFold%reset_level(lvl,ba,dm)
      call this%Cliqold%reset_level(lvl,ba,dm)
      call this%Cgasold%reset_level(lvl,ba,dm)
      
      ! --- PLIC: allocate and set trivial planes for pure cells ---
      call this%PLIC%reset_level(lvl,ba,dm)
      call this%PLICold%reset_level(lvl,ba,dm)
      call set_trivial_plic()
      
      ! --- Fill VOF ghosts (sync + periodic corrections + physical BC) ---
      call this%fill_moments_lvl(lvl,time)
      call this%fill_plic_lvl(lvl,time)
      
   contains
      !> Set PLIC to trivial planes based on VF
      subroutine set_trivial_plic()
         use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPLIC
         integer :: i,j,k
         call amrex_mfiter_build(mfi,this%VF%mf(lvl),tiling=.false.)
         do while (mfi%next())
            bx=mfi%tilebox()
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pPLIC(i,j,k,1:3)=0.0_WP
               pPLIC(i,j,k,4)=sign(1.0e10_WP,pVF(i,j,k,1)-0.5_WP)
            end do; end do; end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine set_trivial_plic
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using conservative interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      
      ! --- Velocity: copy existing + fill new areas from coarse ---
      call this%U%on_remake(this%U,lvl,time,ba,dm)
      call this%V%on_remake(this%V,lvl,time,ba,dm)
      call this%W%on_remake(this%W,lvl,time,ba,dm)
      call this%Uold%reset_level(lvl,ba,dm)
      call this%Vold%reset_level(lvl,ba,dm)
      call this%Wold%reset_level(lvl,ba,dm)
      
      ! --- Pressure: copy existing + fill new areas ---
      call this%P%on_remake(this%P,lvl,time,ba,dm)
      
      ! --- Physical properties (just reset, will be recomputed) ---
      call this%visc%reset_level(lvl,ba,dm)
      
      ! --- VOF moments: copy existing + fill new areas from coarse ---
      call this%VF%on_remake(this%VF,lvl,time,ba,dm)
      call this%Cliq%on_remake(this%Cliq,lvl,time,ba,dm)
      call this%Cgas%on_remake(this%Cgas,lvl,time,ba,dm)
      call this%VFold%reset_level(lvl,ba,dm)
      call this%Cliqold%reset_level(lvl,ba,dm)
      call this%Cgasold%reset_level(lvl,ba,dm)
      
      ! --- PLIC: copy existing + fill new areas, then fix pure cells ---
      call this%PLIC%on_remake(this%PLIC,lvl,time,ba,dm)
      call this%PLICold%reset_level(lvl,ba,dm)
      call fix_pure_plic()
      
      ! --- Fill VOF ghosts (sync + periodic corrections + physical BC) ---
      call this%fill_moments_lvl(lvl,time)
      call this%fill_plic_lvl(lvl,time)
      
   contains
      !> Fix PLIC for pure cells (new cells from coarse have wrong interpolated PLIC)
      subroutine fix_pure_plic()
         use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPLIC
         integer :: i,j,k
         real(WP) :: vf
         call amrex_mfiter_build(mfi,this%VF%mf(lvl),tiling=.false.)
         do while (mfi%next())
            bx=mfi%tilebox()
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               vf=pVF(i,j,k,1)
               if (vf.lt.VFlo.or.vf.gt.VFhi) then
                  pPLIC(i,j,k,1:3)=0.0_WP
                  pPLIC(i,j,k,4)=sign(1.0e10_WP,vf-0.5_WP)
               end if
            end do; end do; end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine fix_pure_plic
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%clear_level(lvl); call this%V%clear_level(lvl); call this%W%clear_level(lvl)
      call this%Uold%clear_level(lvl); call this%Vold%clear_level(lvl); call this%Wold%clear_level(lvl)
      call this%P%clear_level(lvl)
      call this%VF%clear_level(lvl); call this%VFold%clear_level(lvl)
      call this%Cliq%clear_level(lvl); call this%Cliqold%clear_level(lvl)
      call this%Cgas%clear_level(lvl); call this%Cgasold%clear_level(lvl)
      call this%PLIC%clear_level(lvl); call this%PLICold%clear_level(lvl)
      call this%visc%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      integer :: lvl
      ! Average down VOF (VF/Cliq/Cgas + trivial PLIC cleanup at coarse levels)
      call this%vof_average_down(lbase)
      ! Average down velocity for C/F consistency
      do lvl=this%amr%clvl()-1,lbase,-1
         call this%U%average_downto(lvl)
         call this%V%average_downto(lvl)
         call this%W%average_downto(lvl)
         call this%P%average_downto(lvl)
      end do
      ! Fill velocity ghosts
      call this%U%fill(time)
      call this%V%fill(time)
      call this%W%fill(time)
   end subroutine post_regrid

   !> Tag cells near interface with regrid_buffer layer growth
   subroutine tagging(this,lvl,tags,time)
      use amrex_amr_module, only: amrex_tagboxarray,amrex_mfiter,amrex_box,amrex_multifab
      use amrgrid_class, only: SETtag
      use iso_c_binding, only: c_char
      implicit none
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tba
      type(amrex_multifab) :: band
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pBand
      integer :: i,j,k,dir,n,layer
      integer, dimension(3) :: ind
      integer :: eff_buffer
      
      tba=tags
      
      ! Build band MultiFab with 1 ghost cell
      call this%amr%mfab_build(lvl=lvl,mfab=band,ncomp=1,nover=1)
      call band%setval(0.0_WP)
      
      ! Effective buffer (shrink at coarser levels)
      eff_buffer=max(1,this%regrid_buffer/(2**(this%amr%clvl()-lvl)))
      
      ! Pass 1: Mark interface cells (band=1)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         pVF  =>this%VF%mf(lvl)%dataptr(mfi)
         pBand=>band%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Flag all obvious mixture cells
            if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
               pBand(i,j,k,1)=1.0_WP
               cycle
            end if
            ! Check for implicit interfaces (pure cell adjacent to opposite phase)
            do dir=1,3; do n=-1,+1,2
               ind=[i,j,k]; ind(dir)=ind(dir)+n
               if (pVF(i,j,k,1).lt.VFlo.and.pVF(ind(1),ind(2),ind(3),1).gt.VFhi.or.&
               &   pVF(i,j,k,1).gt.VFhi.and.pVF(ind(1),ind(2),ind(3),1).lt.VFlo) then
                  pBand(i,j,k,1)=1.0_WP
                  cycle
               end if
            end do; end do
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      call band%fill_boundary(this%amr%geom(lvl))
      
      ! Pass 2: Grow band by effective buffer layers
      do layer=2,eff_buffer
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pBand=>band%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pBand(i,j,k,1).eq.0.0_WP.and.any(pBand(i-1:i+1,j-1:j+1,k-1:k+1,1).eq.real(layer-1,WP))) pBand(i,j,k,1)=real(layer,WP)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         call band%fill_boundary(this%amr%geom(lvl))
      end do
      
      ! Pass 3: Set tags from band
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tba%dataPtr(mfi)
         pBand =>band%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            if (pBand(i,j,k,1).gt.0.0_WP) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Cleanup
      call this%amr%mfab_destroy(band)
   end subroutine tagging

   !> Estimate per-box costs for load balancing
   !> Cost based on number of mixed cells (SL is expensive) vs pure cells
   subroutine get_cost(this,lvl,nboxes,costs,ba)
      use iso_c_binding, only: c_associated
      use amrex_amr_module, only: amrex_boxarray,amrex_box,amrex_mfiter,&
      &                           amrex_mfiter_build,amrex_mfiter_destroy,amrex_intersection
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      implicit none
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: old_bx,new_bx,isect
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer :: n,i,j,k,ierr
      ! Guard: if VF data doesn't exist yet, return uniform costs
      if (.not.allocated(this%VF%mf)) then; costs=1.0_WP; return; end if
      if (.not.c_associated(this%VF%mf(lvl)%p)) then; costs=1.0_WP; return; end if
      ! Coarser levels are never mixed
      if (lvl.lt.this%amr%clvl()) then
         do n=1,nboxes
            new_bx=ba%get_box(n-1)
            costs(n)=real(new_bx%numpts(),WP)
         end do
         return
      end if
      ! At finest level, count mixed cells per new box from local old data
      costs=0.0_WP
      call amrex_mfiter_build(mfi,this%VF%mf(lvl),tiling=.false.)
      do while (mfi%next())
         old_bx=mfi%tilebox()
         pVF=>this%VF%mf(lvl)%dataptr(mfi)
         do n=1,nboxes
            new_bx=ba%get_box(n-1)
            if (.not.old_bx%intersects(new_bx)) cycle
            isect=amrex_intersection(old_bx,new_bx)
            do k=isect%lo(3),isect%hi(3); do j=isect%lo(2),isect%hi(2); do i=isect%lo(1),isect%hi(1)
               if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
                  costs(n)=costs(n)+this%SLcost
               else
                  costs(n)=costs(n)+1.0_WP
               end if
            end do; end do; end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE,costs,nboxes,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
   end subroutine get_cost

   ! ============================================================================
   ! VOF SYNC / FILL / AVERAGE_DOWN METHODS
   ! ============================================================================

   !> Sync VF/Cliq/Cgas ghosts on all levels + fix periodic barycenters
   subroutine sync_moments(this)
      class(amrmpincomp), intent(inout) :: this
      integer :: lvl
      do lvl=0,this%amr%clvl()
         call this%sync_moments_lvl(lvl)
      end do
   end subroutine sync_moments

   !> Sync VF/Cliq/Cgas ghosts at level + fix periodic barycenters
   subroutine sync_moments_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box, amrex_geometry
      class(amrmpincomp), intent(inout) :: this
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
      dlo = geom%domain%lo; dhi = geom%domain%hi
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
         if (this%amr%xper) then
            if (ilo.lt.dlo(1)) then; do kg=klo,khi; do jg=jlo,jhi; do ig=ilo,dlo(1)-1
               pCliq(ig,jg,kg,1)=pCliq(ig,jg,kg,1)-xL; pCgas(ig,jg,kg,1)=pCgas(ig,jg,kg,1)-xL
            end do; end do; end do; end if
            if (ihi.gt.dhi(1)) then; do kg=klo,khi; do jg=jlo,jhi; do ig=dhi(1)+1,ihi
               pCliq(ig,jg,kg,1)=pCliq(ig,jg,kg,1)+xL; pCgas(ig,jg,kg,1)=pCgas(ig,jg,kg,1)+xL
            end do; end do; end do; end if
         end if
         if (this%amr%yper) then
            if (jlo.lt.dlo(2)) then; do kg=klo,khi; do jg=jlo,dlo(2)-1; do ig=ilo,ihi
               pCliq(ig,jg,kg,2)=pCliq(ig,jg,kg,2)-yL; pCgas(ig,jg,kg,2)=pCgas(ig,jg,kg,2)-yL
            end do; end do; end do; end if
            if (jhi.gt.dhi(2)) then; do kg=klo,khi; do jg=dhi(2)+1,jhi; do ig=ilo,ihi
               pCliq(ig,jg,kg,2)=pCliq(ig,jg,kg,2)+yL; pCgas(ig,jg,kg,2)=pCgas(ig,jg,kg,2)+yL
            end do; end do; end do; end if
         end if
         if (this%amr%zper) then
            if (klo.lt.dlo(3)) then; do kg=klo,dlo(3)-1; do jg=jlo,jhi; do ig=ilo,ihi
               pCliq(ig,jg,kg,3)=pCliq(ig,jg,kg,3)-zL; pCgas(ig,jg,kg,3)=pCgas(ig,jg,kg,3)-zL
            end do; end do; end do; end if
            if (khi.gt.dhi(3)) then; do kg=dhi(3)+1,khi; do jg=jlo,jhi; do ig=ilo,ihi
               pCliq(ig,jg,kg,3)=pCliq(ig,jg,kg,3)+zL; pCgas(ig,jg,kg,3)=pCgas(ig,jg,kg,3)+zL
            end do; end do; end do; end if
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine sync_moments_lvl

   !> Sync PLIC ghosts on all levels + fix periodic plane distance
   subroutine sync_plic(this)
      class(amrmpincomp), intent(inout) :: this
      integer :: lvl
      do lvl = 0, this%amr%clvl()
         call this%sync_plic_lvl(lvl)
      end do
   end subroutine sync_plic

   !> Sync PLIC ghosts at level + fix periodic plane distance
   subroutine sync_plic_lvl(this, lvl)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry
      implicit none
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL
      call this%PLIC%sync_lvl(lvl)
      geom=this%amr%geom(lvl); dlo=geom%domain%lo; dhi=geom%domain%hi
      xL=this%amr%xhi-this%amr%xlo; yL=this%amr%yhi-this%amr%ylo; zL=this%amr%zhi-this%amr%zlo
      call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pP => this%PLIC%mf(lvl)%dataptr(mfi)
         ilo=lbound(pP,1); ihi=ubound(pP,1); jlo=lbound(pP,2); jhi=ubound(pP,2); klo=lbound(pP,3); khi=ubound(pP,3)
         if (this%amr%xper) then
            if (ilo.lt.dlo(1)) then; do kg=klo,khi; do jg=jlo,jhi; do ig=ilo,dlo(1)-1
               pP(ig,jg,kg,4)=pP(ig,jg,kg,4)-pP(ig,jg,kg,1)*xL
            end do; end do; end do; end if
            if (ihi.gt.dhi(1)) then; do kg=klo,khi; do jg=jlo,jhi; do ig=dhi(1)+1,ihi
               pP(ig,jg,kg,4)=pP(ig,jg,kg,4)+pP(ig,jg,kg,1)*xL
            end do; end do; end do; end if
         end if
         if (this%amr%yper) then
            if (jlo.lt.dlo(2)) then; do kg=klo,khi; do jg=jlo,dlo(2)-1; do ig=ilo,ihi
               pP(ig,jg,kg,4)=pP(ig,jg,kg,4)-pP(ig,jg,kg,2)*yL
            end do; end do; end do; end if
            if (jhi.gt.dhi(2)) then; do kg=klo,khi; do jg=dhi(2)+1,jhi; do ig=ilo,ihi
               pP(ig,jg,kg,4)=pP(ig,jg,kg,4)+pP(ig,jg,kg,2)*yL
            end do; end do; end do; end if
         end if
         if (this%amr%zper) then
            if (klo.lt.dlo(3)) then; do kg=klo,dlo(3)-1; do jg=jlo,jhi; do ig=ilo,ihi
               pP(ig,jg,kg,4)=pP(ig,jg,kg,4)-pP(ig,jg,kg,3)*zL
            end do; end do; end do; end if
            if (khi.gt.dhi(3)) then; do kg=dhi(3)+1,khi; do jg=jlo,jhi; do ig=ilo,ihi
               pP(ig,jg,kg,4)=pP(ig,jg,kg,4)+pP(ig,jg,kg,3)*zL
            end do; end do; end do; end if
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine sync_plic_lvl

   !> Fill PLIC ghosts at a level (fill + physical BC)
   subroutine fill_plic_lvl(this, lvl, time)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_geometry, amrex_box
      implicit none
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL, dx, dy, dz
      call this%PLIC%fill_lvl(lvl, time)
      geom=this%amr%geom(lvl); dlo=geom%domain%lo; dhi=geom%domain%hi
      xL=this%amr%xhi-this%amr%xlo; yL=this%amr%yhi-this%amr%ylo; zL=this%amr%zhi-this%amr%zlo
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      call amrex_mfiter_build(mfi, this%PLIC%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pVF  =>this%VF%mf(lvl)%dataptr(mfi); pCliq=>this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas=>this%Cgas%mf(lvl)%dataptr(mfi); pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
         ilo=lbound(pPLIC,1); ihi=ubound(pPLIC,1); jlo=lbound(pPLIC,2); jhi=ubound(pPLIC,2)
         klo=lbound(pPLIC,3); khi=ubound(pPLIC,3)
         ! Fix periodic plane distance
         if (this%amr%xper) then
            if (ilo.lt.dlo(1)) then; do kg=klo,khi; do jg=jlo,jhi; do ig=ilo,dlo(1)-1
               pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)-pPLIC(ig,jg,kg,1)*xL
            end do; end do; end do; end if
            if (ihi.gt.dhi(1)) then; do kg=klo,khi; do jg=jlo,jhi; do ig=dhi(1)+1,ihi
               pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)+pPLIC(ig,jg,kg,1)*xL
            end do; end do; end do; end if
         end if
         if (this%amr%yper) then
            if (jlo.lt.dlo(2)) then; do kg=klo,khi; do jg=jlo,dlo(2)-1; do ig=ilo,ihi
               pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)-pPLIC(ig,jg,kg,2)*yL
            end do; end do; end do; end if
            if (jhi.gt.dhi(2)) then; do kg=klo,khi; do jg=dhi(2)+1,jhi; do ig=ilo,ihi
               pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)+pPLIC(ig,jg,kg,2)*yL
            end do; end do; end do; end if
         end if
         if (this%amr%zper) then
            if (klo.lt.dlo(3)) then; do kg=klo,dlo(3)-1; do jg=jlo,jhi; do ig=ilo,ihi
               pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)-pPLIC(ig,jg,kg,3)*zL
            end do; end do; end do; end if
            if (khi.gt.dhi(3)) then; do kg=dhi(3)+1,khi; do jg=jlo,jhi; do ig=ilo,ihi
               pPLIC(ig,jg,kg,4)=pPLIC(ig,jg,kg,4)+pPLIC(ig,jg,kg,3)*zL
            end do; end do; end do; end if
         end if
         ! Physical BC for PLIC
         if (.not.this%amr%xper) then
            if (ilo.lt.dlo(1)) call apply_bc_face(1,-1,this%vof_lo_bc(1),ilo,dlo(1)-1,jlo,jhi,klo,khi,dlo(1),this%amr%xlo)
            if (ihi.gt.dhi(1)) call apply_bc_face(1,+1,this%vof_hi_bc(1),dhi(1)+1,ihi,jlo,jhi,klo,khi,dhi(1),this%amr%xhi)
         end if
         if (.not.this%amr%yper) then
            if (jlo.lt.dlo(2)) call apply_bc_face(2,-1,this%vof_lo_bc(2),ilo,ihi,jlo,dlo(2)-1,klo,khi,dlo(2),this%amr%ylo)
            if (jhi.gt.dhi(2)) call apply_bc_face(2,+1,this%vof_hi_bc(2),ilo,ihi,dhi(2)+1,jhi,klo,khi,dhi(2),this%amr%yhi)
         end if
         if (.not.this%amr%zper) then
            if (klo.lt.dlo(3)) call apply_bc_face(3,-1,this%vof_lo_bc(3),ilo,ihi,jlo,jhi,klo,dlo(3)-1,dlo(3),this%amr%zlo)
            if (khi.gt.dhi(3)) call apply_bc_face(3,+1,this%vof_hi_bc(3),ilo,ihi,jlo,jhi,dhi(3)+1,khi,dhi(3),this%amr%zhi)
         end if
         ! Enforce strict VF/PLIC consistency
         do kg=klo,khi; do jg=jlo,jhi; do ig=ilo,ihi
            if (pVF(ig,jg,kg,1).lt.VFlo.or.pVF(ig,jg,kg,1).gt.VFhi) &
            & pPLIC(ig,jg,kg,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0e10_WP,pVF(ig,jg,kg,1)-0.5_WP)]
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   contains
      subroutine apply_bc_face(dir,side,bc_type,i1,i2,j1,j2,k1,k2,bnd,x_bnd)
         integer, intent(in) :: dir,side,bc_type,i1,i2,j1,j2,k1,k2,bnd
         real(WP), intent(in) :: x_bnd
         integer :: ig2,jg2,kg2,isrc,jsrc,ksrc,face
         type(amrex_box) :: bc_bx
         select case (bc_type)
          case (BC_LIQ)
            do kg2=k1,k2; do jg2=j1,j2; do ig2=i1,i2
               pPLIC(ig2,jg2,kg2,1:3)=0.0_WP; pPLIC(ig2,jg2,kg2,4)=1.0e10_WP
            end do; end do; end do
          case (BC_GAS)
            do kg2=k1,k2; do jg2=j1,j2; do ig2=i1,i2
               pPLIC(ig2,jg2,kg2,1:3)=0.0_WP; pPLIC(ig2,jg2,kg2,4)=-1.0e10_WP
            end do; end do; end do
          case (BC_REFLECT)
            do kg2=k1,k2; do jg2=j1,j2; do ig2=i1,i2
               isrc=ig2; jsrc=jg2; ksrc=kg2
               if (dir.eq.1) isrc=2*bnd-ig2+side
               if (dir.eq.2) jsrc=2*bnd-jg2+side
               if (dir.eq.3) ksrc=2*bnd-kg2+side
               pPLIC(ig2,jg2,kg2,1:4)=pPLIC(isrc,jsrc,ksrc,1:4)
               pPLIC(ig2,jg2,kg2,dir)=-pPLIC(ig2,jg2,kg2,dir)
               pPLIC(ig2,jg2,kg2,4)=pPLIC(ig2,jg2,kg2,4)-2.0_WP*pPLIC(isrc,jsrc,ksrc,dir)*x_bnd
            end do; end do; end do
          case (BC_USER)
            if (associated(this%user_vof_bc)) then
               bc_bx=amrex_box([i1,j1,k1],[i2,j2,k2])
               face=2*dir-1+(1+side)/2
               call this%user_vof_bc(this,bc_bx,pVF,pCliq,pCgas,pPLIC,face,time,1)
            end if
         end select
      end subroutine apply_bc_face
   end subroutine fill_plic_lvl

   !> Fill moments ghosts at a level (fill + physical BC)
   subroutine fill_moments_lvl(this, lvl, time)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box, amrex_geometry
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_mfiter) :: mfi
      type(amrex_geometry) :: geom
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas, pPLIC
      integer :: ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi
      integer :: dlo(3), dhi(3)
      real(WP) :: xL, yL, zL, dx, dy, dz
      call this%VF%fill_lvl(lvl,time); call this%Cliq%fill_lvl(lvl,time); call this%Cgas%fill_lvl(lvl,time)
      geom=this%amr%geom(lvl); dlo=geom%domain%lo; dhi=geom%domain%hi
      xL=this%amr%xhi-this%amr%xlo; yL=this%amr%yhi-this%amr%ylo; zL=this%amr%zhi-this%amr%zlo
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      call amrex_mfiter_build(mfi, this%VF%mf(lvl), tiling=.false.)
      do while (mfi%next())
         pVF  =>this%VF%mf(lvl)%dataptr(mfi); pCliq=>this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas=>this%Cgas%mf(lvl)%dataptr(mfi); pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
         ilo=lbound(pVF,1); ihi=ubound(pVF,1); jlo=lbound(pVF,2); jhi=ubound(pVF,2)
         klo=lbound(pVF,3); khi=ubound(pVF,3)
         ! Fix periodic barycenters
         if (this%amr%xper) then
            if (ilo.lt.dlo(1)) then; do kg=klo,khi; do jg=jlo,jhi; do ig=ilo,dlo(1)-1
               pCliq(ig,jg,kg,1)=pCliq(ig,jg,kg,1)-xL; pCgas(ig,jg,kg,1)=pCgas(ig,jg,kg,1)-xL
            end do; end do; end do; end if
            if (ihi.gt.dhi(1)) then; do kg=klo,khi; do jg=jlo,jhi; do ig=dhi(1)+1,ihi
               pCliq(ig,jg,kg,1)=pCliq(ig,jg,kg,1)+xL; pCgas(ig,jg,kg,1)=pCgas(ig,jg,kg,1)+xL
            end do; end do; end do; end if
         end if
         if (this%amr%yper) then
            if (jlo.lt.dlo(2)) then; do kg=klo,khi; do jg=jlo,dlo(2)-1; do ig=ilo,ihi
               pCliq(ig,jg,kg,2)=pCliq(ig,jg,kg,2)-yL; pCgas(ig,jg,kg,2)=pCgas(ig,jg,kg,2)-yL
            end do; end do; end do; end if
            if (jhi.gt.dhi(2)) then; do kg=klo,khi; do jg=dhi(2)+1,jhi; do ig=ilo,ihi
               pCliq(ig,jg,kg,2)=pCliq(ig,jg,kg,2)+yL; pCgas(ig,jg,kg,2)=pCgas(ig,jg,kg,2)+yL
            end do; end do; end do; end if
         end if
         if (this%amr%zper) then
            if (klo.lt.dlo(3)) then; do kg=klo,dlo(3)-1; do jg=jlo,jhi; do ig=ilo,ihi
               pCliq(ig,jg,kg,3)=pCliq(ig,jg,kg,3)-zL; pCgas(ig,jg,kg,3)=pCgas(ig,jg,kg,3)-zL
            end do; end do; end do; end if
            if (khi.gt.dhi(3)) then; do kg=dhi(3)+1,khi; do jg=jlo,jhi; do ig=ilo,ihi
               pCliq(ig,jg,kg,3)=pCliq(ig,jg,kg,3)+zL; pCgas(ig,jg,kg,3)=pCgas(ig,jg,kg,3)+zL
            end do; end do; end do; end if
         end if
         ! Physical BC for moments
         if (.not.this%amr%xper) then
            if (ilo.lt.dlo(1)) call apply_bc_face(1,-1,this%vof_lo_bc(1),ilo,dlo(1)-1,jlo,jhi,klo,khi,dlo(1),this%amr%xlo)
            if (ihi.gt.dhi(1)) call apply_bc_face(1,+1,this%vof_hi_bc(1),dhi(1)+1,ihi,jlo,jhi,klo,khi,dhi(1),this%amr%xhi)
         end if
         if (.not.this%amr%yper) then
            if (jlo.lt.dlo(2)) call apply_bc_face(2,-1,this%vof_lo_bc(2),ilo,ihi,jlo,dlo(2)-1,klo,khi,dlo(2),this%amr%ylo)
            if (jhi.gt.dhi(2)) call apply_bc_face(2,+1,this%vof_hi_bc(2),ilo,ihi,dhi(2)+1,jhi,klo,khi,dhi(2),this%amr%yhi)
         end if
         if (.not.this%amr%zper) then
            if (klo.lt.dlo(3)) call apply_bc_face(3,-1,this%vof_lo_bc(3),ilo,ihi,jlo,jhi,klo,dlo(3)-1,dlo(3),this%amr%zlo)
            if (khi.gt.dhi(3)) call apply_bc_face(3,+1,this%vof_hi_bc(3),ilo,ihi,jlo,jhi,dhi(3)+1,khi,dhi(3),this%amr%zhi)
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   contains
      subroutine apply_bc_face(dir,side,bc_type,i1,i2,j1,j2,k1,k2,bnd,x_bnd)
         integer, intent(in) :: dir,side,bc_type,i1,i2,j1,j2,k1,k2,bnd
         real(WP), intent(in) :: x_bnd
         integer :: ig2,jg2,kg2,isrc,jsrc,ksrc,face
         real(WP), dimension(3) :: center
         type(amrex_box) :: bc_bx
         select case (bc_type)
          case (BC_LIQ)
            do kg2=k1,k2; do jg2=j1,j2; do ig2=i1,i2
               center=[this%amr%xlo+(real(ig2,WP)+0.5_WP)*dx,this%amr%ylo+(real(jg2,WP)+0.5_WP)*dy,this%amr%zlo+(real(kg2,WP)+0.5_WP)*dz]
               pVF(ig2,jg2,kg2,1)=1.0_WP; pCliq(ig2,jg2,kg2,1:3)=center; pCgas(ig2,jg2,kg2,1:3)=center
            end do; end do; end do
          case (BC_GAS)
            do kg2=k1,k2; do jg2=j1,j2; do ig2=i1,i2
               center=[this%amr%xlo+(real(ig2,WP)+0.5_WP)*dx,this%amr%ylo+(real(jg2,WP)+0.5_WP)*dy,this%amr%zlo+(real(kg2,WP)+0.5_WP)*dz]
               pVF(ig2,jg2,kg2,1)=0.0_WP; pCliq(ig2,jg2,kg2,1:3)=center; pCgas(ig2,jg2,kg2,1:3)=center
            end do; end do; end do
          case (BC_REFLECT)
            do kg2=k1,k2; do jg2=j1,j2; do ig2=i1,i2
               isrc=ig2; jsrc=jg2; ksrc=kg2
               if (dir.eq.1) isrc=2*bnd-ig2+side
               if (dir.eq.2) jsrc=2*bnd-jg2+side
               if (dir.eq.3) ksrc=2*bnd-kg2+side
               pVF(ig2,jg2,kg2,1)=pVF(isrc,jsrc,ksrc,1)
               pCliq(ig2,jg2,kg2,1:3)=pCliq(isrc,jsrc,ksrc,1:3); pCgas(ig2,jg2,kg2,1:3)=pCgas(isrc,jsrc,ksrc,1:3)
               pCliq(ig2,jg2,kg2,dir)=2.0_WP*x_bnd-pCliq(isrc,jsrc,ksrc,dir)
               pCgas(ig2,jg2,kg2,dir)=2.0_WP*x_bnd-pCgas(isrc,jsrc,ksrc,dir)
            end do; end do; end do
          case (BC_USER)
            if (associated(this%user_vof_bc)) then
               bc_bx=amrex_box([i1,j1,k1],[i2,j2,k2])
               face=2*dir-1+(1+side)/2
               call this%user_vof_bc(this,bc_bx,pVF,pCliq,pCgas,pPLIC,face,time,2)
            end if
         end select
      end subroutine apply_bc_face
   end subroutine fill_moments_lvl

   !> Average down VF/Cliq/Cgas from finest to lbase, then sync ghost cells
   subroutine vof_average_down(this,lbase)
      use amrex_interface, only: amrmfab_average_down_cell
      class(amrmpincomp), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl, lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=this%amr%clvl()-1,lb,-1
         call amrmfab_average_down_cell(fmf=this%VF%mf(lvl+1)  ,cmf=this%VF%mf(lvl)  ,rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_cell(fmf=this%Cliq%mf(lvl+1),cmf=this%Cliq%mf(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_cell(fmf=this%Cgas%mf(lvl+1),cmf=this%Cgas%mf(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
      end do
      call this%sync_moments()
      do lvl=this%amr%clvl()-1,lb,-1
         call set_trivial_plic()
         call this%sync_plic_lvl(lvl)
      end do
   contains
      subroutine set_trivial_plic()
         use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPLIC
         integer :: i,j,k
         call amrex_mfiter_build(mfi,this%PLIC%mf(lvl),tiling=.false.)
         do while (mfi%next())
            pVF  =>this%VF%mf(lvl)%dataptr(mfi); pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0e10_WP,pVF(i,j,k,1)-0.5_WP)]
            end do; end do; end do
         end do
         call amrex_mfiter_destroy(mfi)
      end subroutine set_trivial_plic
   end subroutine vof_average_down

   ! ============================================================================
   ! PLIC RECONSTRUCTION AND POLYGON METHODS
   ! ============================================================================

   !> Build PLIC interface at finest level using PLICnet
   subroutine build_plic(this,time)
      use plicnet, only: get_normal,reflect_moments
      use mathtools, only: normalize
      use amrvof_geometry, only: get_plane_dist
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08, only: MPI_Wtime
      class(amrmpincomp), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      real(WP) :: t0,t1

      ! Start full routine timer
      t0=MPI_Wtime()

      ! Only build at finest level
      lvl=this%amr%clvl()

      ! Get cell size at this level
      dx=this%amr%dx(lvl); dxi=1.0_WP/dx
      dy=this%amr%dy(lvl); dyi=1.0_WP/dy
      dz=this%amr%dz(lvl); dzi=1.0_WP/dz

      ! Compute PLIC planes
      t1=MPI_Wtime()
      plic_reconstruction: block
         integer :: i,j,k,ii,jj,kk,direction,direction2
         real(WP), dimension(0:188) :: moments
         real(WP), dimension(3) :: normal,center,lo,hi
         real(WP) :: m000,m100,m010,m001,temp
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCliq,pCgas,pPLIC
         logical :: flip
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx

         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()

            ! Get pointers to data
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pCliq=>this%Cliq%mf(lvl)%dataptr(mfi)
            pCgas=>this%Cgas%mf(lvl)%dataptr(mfi)
            pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)

            ! Loop over cells in this box
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)

               ! Handle full cells: set trivial plane
               if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi) then
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0e10_WP,pVF(i,j,k,1)-0.5_WP)]
                  cycle
               end if

               ! Handle partial cells: compute moments
               flip=.false.; if (pVF(i,j,k,1).ge.0.5_WP) flip=.true.
               m000=0.0_WP; m100=0.0_WP; m010=0.0_WP; m001=0.0_WP
               if (flip) then
                  do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+0)=1.0_WP-pVF(ii,jj,kk,1)
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(pCgas(ii,jj,kk,1)-(this%amr%xlo+(real(ii,WP)+0.5_WP)*dx))*dxi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(pCgas(ii,jj,kk,2)-(this%amr%ylo+(real(jj,WP)+0.5_WP)*dy))*dyi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(pCgas(ii,jj,kk,3)-(this%amr%zlo+(real(kk,WP)+0.5_WP)*dz))*dzi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(pCliq(ii,jj,kk,1)-(this%amr%xlo+(real(ii,WP)+0.5_WP)*dx))*dxi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(pCliq(ii,jj,kk,2)-(this%amr%ylo+(real(jj,WP)+0.5_WP)*dy))*dyi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(pCliq(ii,jj,kk,3)-(this%amr%zlo+(real(kk,WP)+0.5_WP)*dz))*dzi
                     m000=m000+moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                     m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                     m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                     m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                  end do; end do; end do
               else
                  do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+0)=pVF(ii,jj,kk,1)
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(pCliq(ii,jj,kk,1)-(this%amr%xlo+(real(ii,WP)+0.5_WP)*dx))*dxi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(pCliq(ii,jj,kk,2)-(this%amr%ylo+(real(jj,WP)+0.5_WP)*dy))*dyi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(pCliq(ii,jj,kk,3)-(this%amr%zlo+(real(kk,WP)+0.5_WP)*dz))*dzi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(pCgas(ii,jj,kk,1)-(this%amr%xlo+(real(ii,WP)+0.5_WP)*dx))*dxi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(pCgas(ii,jj,kk,2)-(this%amr%ylo+(real(jj,WP)+0.5_WP)*dy))*dyi
                     moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(pCgas(ii,jj,kk,3)-(this%amr%zlo+(real(kk,WP)+0.5_WP)*dz))*dzi
                     m000=m000+moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                     m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                     m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                     m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))
                  end do; end do; end do
               end if
               if (m000.gt.tiny(1.0_WP)) then; center=[m100,m010,m001]/m000; else; center=0.0_WP; end if
               call reflect_moments(moments,center,direction,direction2)
               call get_normal(moments,normal); normal=normalize(normal)
               ! Undo direction2 rotation
               if      (direction2.eq.1) then; temp=normal(1); normal(1)=normal(2); normal(2)=temp
               else if (direction2.eq.2) then; temp=normal(2); normal(2)=normal(3); normal(3)=temp
               else if (direction2.eq.3) then; temp=normal(1); normal(1)=normal(3); normal(3)=temp
               else if (direction2.eq.4) then; temp=normal(2); normal(2)=normal(3); normal(3)=temp; temp=normal(1); normal(1)=normal(2); normal(2)=temp
               else if (direction2.eq.5) then; temp=normal(1); normal(1)=normal(3); normal(3)=temp; temp=normal(1); normal(1)=normal(2); normal(2)=temp
               end if
               ! Undo direction reflection
               if      (direction.eq.1) then; normal(1)=-normal(1)
               else if (direction.eq.2) then; normal(2)=-normal(2)
               else if (direction.eq.3) then; normal(3)=-normal(3)
               else if (direction.eq.4) then; normal(1)=-normal(1); normal(2)=-normal(2)
               else if (direction.eq.5) then; normal(1)=-normal(1); normal(3)=-normal(3)
               else if (direction.eq.6) then; normal(2)=-normal(2); normal(3)=-normal(3)
               else if (direction.eq.7) then; normal(1)=-normal(1); normal(2)=-normal(2); normal(3)=-normal(3)
               end if
               if (.not.flip) normal=-normal

               ! General final PLIC data
               normal=normalize(normal)
               lo=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
               hi=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
               pPLIC(i,j,k,:)=[normal(1),normal(2),normal(3),get_plane_dist(normal,lo,hi,pVF(i,j,k,1))]

            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end block plic_reconstruction
      this%wt_plicnet=this%wt_plicnet+(MPI_Wtime()-t1)

      ! Fill PLIC ghosts (sync + periodic correction + physical BC)
      call this%fill_plic_lvl(lvl, time)

      ! Build polygons from PLIC
      call this%build_polygons()

      ! Stop full routine timer
      this%wt_plic=this%wt_plic+(MPI_Wtime()-t0)

   end subroutine build_plic

   !> Build polygons from PLIC planes
   subroutine build_polygons(this)
      use amrvof_geometry, only: cut_hex_polygon
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08, only: MPI_Wtime
      class(amrmpincomp), intent(inout) :: this
      integer :: lvl,i,j,k
      real(WP) :: dx,dy,dz,t0
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLIC
      real(WP), dimension(3) :: lo,hi
      real(WP), dimension(4) :: plane
      real(WP), dimension(3,8) :: hex
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx,gbx
      real(WP), dimension(:,:,:,:,:), allocatable :: polygon_local
      integer, dimension(:,:,:), allocatable :: poly_nv_local
      real(WP), dimension(3,6) :: poly_verts
      integer :: poly_nv

      ! Start timer
      t0=MPI_Wtime()

      ! Get level and cell size
      lvl=this%amr%clvl()
      dx=this%amr%dx(lvl)
      dy=this%amr%dy(lvl)
      dz=this%amr%dz(lvl)

      ! Reset polygon storage
      call this%smesh%reset()

      ! Compute new polygons
      call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
      do while (mfi%next())
         ! Get local and grown boxes
         bx=mfi%tilebox()
         gbx=mfi%growntilebox(2)

         ! Get pointer to PLIC data
         pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)

         ! Allocate per-grownbox polygon storage
         allocate(polygon_local(1:3,1:6,gbx%lo(1):gbx%hi(1),gbx%lo(2):gbx%hi(2),gbx%lo(3):gbx%hi(3))); polygon_local=0.0_WP
         allocate(poly_nv_local        (gbx%lo(1):gbx%hi(1),gbx%lo(2):gbx%hi(2),gbx%lo(3):gbx%hi(3))); poly_nv_local=0

         ! Extract polygons on grown box including ghosts
         do k=gbx%lo(3),gbx%hi(3); do j=gbx%lo(2),gbx%hi(2); do i=gbx%lo(1),gbx%hi(1)
            ! Skip cells with no interface
            if (abs(pPLIC(i,j,k,4)).gt.1.0e+9_WP) cycle
            ! Build hex and plane
            lo=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hi=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            plane=[pPLIC(i,j,k,1),pPLIC(i,j,k,2),pPLIC(i,j,k,3),pPLIC(i,j,k,4)]
            hex(:,1)=[hi(1),lo(2),lo(3)]; hex(:,2)=[hi(1),hi(2),lo(3)]; hex(:,3)=[hi(1),hi(2),hi(3)]; hex(:,4)=[hi(1),lo(2),hi(3)]
            hex(:,5)=[lo(1),lo(2),lo(3)]; hex(:,6)=[lo(1),hi(2),lo(3)]; hex(:,7)=[lo(1),hi(2),hi(3)]; hex(:,8)=[lo(1),lo(2),hi(3)]
            call cut_hex_polygon(hex,plane,poly_nv,poly_verts)
            ! Store polygon
            poly_nv_local(i,j,k)=poly_nv
            if (poly_nv.ge.3) polygon_local(:,1:poly_nv,i,j,k)=poly_verts(:,1:poly_nv)
         end do; end do; end do

         ! Compute curvature here

         ! Append polygons to smesh
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            poly_nv=poly_nv_local(i,j,k)
            if (poly_nv.ge.3) call this%smesh%add_polygon(polygon_local(:,1:poly_nv,i,j,k),poly_nv)
         end do; end do; end do

         ! Deallocate local arrays
         deallocate(polygon_local,poly_nv_local)

      end do
      call this%amr%mfiter_destroy(mfi)

      ! Stop timer
      this%wt_polygon=this%wt_polygon+(MPI_Wtime()-t0)
      
   end subroutine build_polygons

   !> Reset VF and barycenters from PLIC plane to ensure consistency
   !> Computes in valid + ghost cells from PLIC (which is already filled)
   subroutine reset_moments(this)
      use amrvof_geometry, only: cut_hex_vol
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      class(amrmpincomp), intent(inout) :: this
      integer :: lvl,i,j,k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCliq,pCgas,pPLIC
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane
      real(WP) :: vol_liq,vol_gas,cell_vol,dx,dy,dz
      real(WP), dimension(3) :: bary_liq,bary_gas,cell_center
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      
      ! Only work at finest level
      lvl=this%amr%clvl()
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      cell_vol=dx*dy*dz
      
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pVF  =>this%VF%mf(lvl)%dataptr(mfi)
         pCliq=>this%Cliq%mf(lvl)%dataptr(mfi)
         pCgas=>this%Cgas%mf(lvl)%dataptr(mfi)
         pPLIC=>this%PLIC%mf(lvl)%dataptr(mfi)
         ! Loop over grown tiles
         bx=mfi%growntilebox(this%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Cell center
            cell_center=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx, &
            &            this%amr%ylo+(real(j,WP)+0.5_WP)*dy, &
            &            this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
            ! Build hex cell (8 vertices)
            hex(:,1)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,2)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,3)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,4)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hex(:,5)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,6)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,7)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            hex(:,8)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            ! Get plane from PLIC
            plane=pPLIC(i,j,k,:)
            ! Skip cutting for full cells (trivial PLIC with large distance)
            if (abs(plane(4)).ge.1.0e9_WP) then
               if (plane(4).gt.0.0_WP) then
                  pVF(i,j,k,1)=1.0_WP
               else
                  pVF(i,j,k,1)=0.0_WP
               end if
               pCliq(i,j,k,1:3)=cell_center
               pCgas(i,j,k,1:3)=cell_center
               cycle
            end if
            ! Cut hex by plane
            call cut_hex_vol(hex,plane,vol_liq,vol_gas,bary_liq,bary_gas)
            ! Update VF and barycenters
            pVF(i,j,k,1)=vol_liq/cell_vol
            pCliq(i,j,k,1:3)=bary_liq
            pCgas(i,j,k,1:3)=bary_gas
            ! Clean up edge cases
            if (pVF(i,j,k,1).lt.VFlo) then
               pVF(i,j,k,1)=0.0_WP
               pCliq(i,j,k,1:3)=cell_center
               pCgas(i,j,k,1:3)=cell_center
            end if
            if (pVF(i,j,k,1).gt.VFhi) then
               pVF(i,j,k,1)=1.0_WP
               pCliq(i,j,k,1:3)=cell_center
               pCgas(i,j,k,1:3)=cell_center
            end if
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      
      ! Average down VF/Cliq/Cgas to coarse levels + clean up PLIC
      call this%vof_average_down()
      
   end subroutine reset_moments

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Vreman SGS viscosity model adapted for incompressible multiphase
   subroutine add_vreman(this,dt,Cs)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpincomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      ! Local variables
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: visc_t,scratch
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc_t,pScratch,pU,pV,pW,pVisc,pVF
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_visc,Cmodel,Aij,Bij,t0
      real(WP), dimension(1:3,1:3) :: gradU,betaij
      integer :: lvl,i,j,k,si,sj,sk,n
      ! Parameters
      real(WP), parameter :: max_cfl=0.5_WP
      integer, parameter :: nfilter=2
      real(WP), dimension(-1:+1), parameter :: filter=[1.0_WP/6.0_WP,2.0_WP/3.0_WP,1.0_WP/6.0_WP]
      
      ! Start timer
      t0=MPI_Wtime()

      ! Model constant: c=2.5*Cs**2 (Vreman uses c=0.07 which corresponds to Cs=0.17)
      if (present(Cs)) then; Cmodel=2.5_WP*Cs**2; else; Cmodel=2.5_WP*0.17_WP**2; end if
      
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         
         ! Grid spacings
         dx=this%amr%dx(lvl); dxi=1.0_WP/dx
         dy=this%amr%dy(lvl); dyi=1.0_WP/dy
         dz=this%amr%dz(lvl); dzi=1.0_WP/dz
         
         ! Max visc from CFL
         max_visc=max_cfl*this%amr%min_meshsize(lvl)**2/(4.0_WP*dt)
         
         ! Build temp multifab for eddy viscosity at this level (nover ghost cells)
         call this%amr%mfab_build(lvl=lvl,mfab=visc_t,ncomp=1,nover=this%nover); call visc_t%setval(0.0_WP)
         
         ! Phase 1: Compute kinematic eddy viscosity
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            ! Get data pointers
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pVisc_t=>visc_t%dataptr(mfi)
            ! Loop over interior tiles
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Compute cell-centered velocity gradient tensor
               gradU(1,1)=0.5_WP*dxi*(pU(i+1,j,k,1)-pU(i-1,j,k,1))
               gradU(2,1)=0.5_WP*dyi*(pU(i,j+1,k,1)-pU(i,j-1,k,1))
               gradU(3,1)=0.5_WP*dzi*(pU(i,j,k+1,1)-pU(i,j,k-1,1))
               gradU(1,2)=0.5_WP*dxi*(pV(i+1,j,k,1)-pV(i-1,j,k,1))
               gradU(2,2)=0.5_WP*dyi*(pV(i,j+1,k,1)-pV(i,j-1,k,1))
               gradU(3,2)=0.5_WP*dzi*(pV(i,j,k+1,1)-pV(i,j,k-1,1))
               gradU(1,3)=0.5_WP*dxi*(pW(i+1,j,k,1)-pW(i-1,j,k,1))
               gradU(2,3)=0.5_WP*dyi*(pW(i,j+1,k,1)-pW(i,j-1,k,1))
               gradU(3,3)=0.5_WP*dzi*(pW(i,j,k+1,1)-pW(i,j,k-1,1))
               ! Compute A=gradU_ij*gradU_ij invariant
               Aij=sum(gradU**2)
               ! Compute beta_ij=dx_m^2*gradU_mi*gradU_mj
               do sj=1,3; do si=1,3
                  betaij(si,sj)=dx**2*gradU(1,si)*gradU(1,sj)+dy**2*gradU(2,si)*gradU(2,sj)+dz**2*gradU(3,si)*gradU(3,sj)
               end do; end do
               ! Compute B invariant
               Bij=betaij(1,1)*betaij(2,2)-betaij(1,2)**2+betaij(1,1)*betaij(3,3)-betaij(1,3)**2+betaij(2,2)*betaij(3,3)-betaij(2,3)**2
               ! Assemble eddy viscosity
               if (Bij.gt.0.0_WP) then
                  pVisc_t(i,j,k,1)=Cmodel*sqrt(Bij/Aij)
               end if
               ! Clip to CFL limit
               pVisc_t(i,j,k,1)=min(pVisc_t(i,j,k,1),max_visc)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         
         ! Phase 2: Filter visc_t
         call this%amr%mfab_build(lvl=lvl,mfab=scratch,ncomp=1,nover=1); call scratch%setval(0.0_WP)
         do n=1,nfilter
            call scratch%setval(0.0_WP)
            call scratch%copy(srcmf=visc_t,srccomp=1,dstcomp=1,nc=1,ng=0)
            call scratch%fill_boundary(this%amr%geom(lvl))
            call this%amr%mfab_foextrap(lvl=lvl,mfab=scratch)
            call this%amr%mfiter_build(lvl,mfi)
            do while(mfi%next())
               pScratch=>scratch%dataptr(mfi)
               pVisc_t=>visc_t%dataptr(mfi)
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pVisc_t(i,j,k,1)=0.0_WP
                  do sk=-1,+1; do sj=-1,+1; do si=-1,+1
                     pVisc_t(i,j,k,1)=pVisc_t(i,j,k,1)+filter(si)*filter(sj)*filter(sk)*pScratch(i+si,j+sj,k+sk,1)
                  end do; end do; end do
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
         call amrex_multifab_destroy(scratch)
         
         ! Fill visc_t ghosts after filtering
         call visc_t%fill_boundary(this%amr%geom(lvl))
         call this%amr%mfab_foextrap(lvl=lvl,mfab=visc_t)

         ! Phase 3: Convert to dynamic viscosity via harmonic averaging and add to this%visc
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pVisc_t=>visc_t%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pVisc(i,j,k,1)=pVisc(i,j,k,1)+pVisc_t(i,j,k,1)/(pVF(i,j,k,1)/this%rhoL+(1.0_WP-pVF(i,j,k,1))/this%rhoG)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)

         ! Destroy temp multifab
         call amrex_multifab_destroy(visc_t)
         
      end do
      
      ! End timer
      this%wt_visc=this%wt_visc+(MPI_Wtime()-t0)

   end subroutine add_vreman

   !> Calculate CFL numbers (convective + viscous, no acoustic)
   subroutine get_cfl(this,dt,cfl)
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      class(amrmpincomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      integer :: lvl,i,j,k,ierr
      real(WP) :: Umax,Vmax,Wmax,viscmax,rho
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc,pVF
      ! Reset CFLs
      this%CFLc_x=0.0_WP; this%CFLc_y=0.0_WP; this%CFLc_z=0.0_WP
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      ! Compute CFL at each level (finest level determines dt)
      do lvl=0,this%amr%clvl()
         ! Max velocity
         Umax=this%U%norm0(lvl=lvl)
         Vmax=this%V%norm0(lvl=lvl)
         Wmax=this%W%norm0(lvl=lvl)
         ! Max kinematic viscosity
         get_viscmax: block
            viscmax=0.0_WP
            call this%amr%mfiter_build(lvl,mfi)
            do while(mfi%next())
               ! Get data pointers
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               ! Loop over interior tiles
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho=pVF(i,j,k,1)*this%rhoL+(1.0_WP-pVF(i,j,k,1))*this%rhoG
                  viscmax=max(viscmax,pVisc(i,j,k,1)/rho)
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            call MPI_ALLREDUCE(MPI_IN_PLACE,viscmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         end block get_viscmax
         ! Convective
         if (this%amr%nx.gt.1) this%CFLc_x=max(this%CFLc_x,Umax*dt/this%amr%dx(lvl))
         if (this%amr%ny.gt.1) this%CFLc_y=max(this%CFLc_y,Vmax*dt/this%amr%dy(lvl))
         if (this%amr%nz.gt.1) this%CFLc_z=max(this%CFLc_z,Wmax*dt/this%amr%dz(lvl))
         ! Viscous
         if (this%amr%nx.gt.1) this%CFLv_x=max(this%CFLv_x,4.0_WP*viscmax*dt/this%amr%dx(lvl)**2)
         if (this%amr%ny.gt.1) this%CFLv_y=max(this%CFLv_y,4.0_WP*viscmax*dt/this%amr%dy(lvl)**2)
         if (this%amr%nz.gt.1) this%CFLv_z=max(this%CFLv_z,4.0_WP*viscmax*dt/this%amr%dz(lvl)**2)
      end do
      ! Return max CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
   end subroutine get_cfl

   !> Calculate monitoring info
   subroutine get_info(this)
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_MAX,MPI_MIN
      implicit none
      class(amrmpincomp), intent(inout) :: this
      integer :: lvl,ierr
      real(WP) :: dV

      ! Velocity/pressure/VF extrema
      extrema: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pU,pV,pW
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         integer :: i,j,k
         real(WP) :: rho
         ! Initialize extrema
         this%Umax=0.0_WP; this%Vmax=0.0_WP; this%Wmax=0.0_WP
         this%Pmax=0.0_WP
         this%VFmin=huge(1.0_WP); this%VFmax=-huge(1.0_WP); this%VFint=0.0_WP
         this%rhoKint=0.0_WP
         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Velocity norm 0
            this%Umax=max(this%Umax,this%U%norm0(lvl=lvl))
            this%Vmax=max(this%Vmax,this%V%norm0(lvl=lvl))
            this%Wmax=max(this%Wmax,this%W%norm0(lvl=lvl))
            ! Pressure norm 0
            this%Pmax=max(this%Pmax,this%P%norm0(lvl=lvl))
            ! Extrema of volume fraction
            this%VFmin=min(this%VFmin,this%VF%get_min(lvl=lvl))
            this%VFmax=max(this%VFmax,this%VF%get_max(lvl=lvl))
         end do
      end block extrema

      ! Conserved integrals at base level
      dV=this%amr%cell_vol(0)
      this%VFint=this%VF%get_sum(lvl=0)*dV

      ! Kinetic energy integral: 0.5 * rho * (U^2 + V^2 + W^2) * dV
      ! Uses composite integration with fine masking to avoid double-counting
      get_rhoKint: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pU,pV,pW
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         integer :: i,j,k
         real(WP) :: rho
         do lvl=0,this%amr%clvl()
            ! Get cell volume
            dV=this%amr%cell_vol(lvl)
            ! Build fine mask for this level (if not finest)
            if (lvl.lt.this%amr%clvl()) then
               call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
            end if
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               bx=mfi%tilebox()
               ! Get pointers to data
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               pU =>this%U%mf(lvl)%dataptr(mfi)
               pV =>this%V%mf(lvl)%dataptr(mfi)
               pW =>this%W%mf(lvl)%dataptr(mfi)
               ! Get pointer to fine mask
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over cells
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
                  ! Accumulate kinetic energy (rho = VF*rhoL + (1-VF)*rhoG)
                  rho=pVF(i,j,k,1)*this%rhoL+(1.0_WP-pVF(i,j,k,1))*this%rhoG
                  this%rhoKint=this%rhoKint+0.5_WP*rho*(pU(i,j,k,1)**2+pV(i,j,k,1)**2+pW(i,j,k,1)**2)*dV
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         ! Reduce across MPI ranks
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block get_rhoKint

      ! Load distribution diagnostics at finest level
      load_distribution: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: ncells_local,nmixed_local
         integer :: i,j,k
         ncells_local=0.0_WP;nmixed_local=0.0_WP
         lvl=this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ncells_local=ncells_local+1.0_WP
               if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) nmixed_local=nmixed_local+1.0_WP
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         this%ncells_max=ncells_local; this%ncells_min=ncells_local
         this%nmixed_max=nmixed_local; this%nmixed_min=nmixed_local
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%ncells_max,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%ncells_min,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nmixed_max,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nmixed_min,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      end block load_distribution

      ! Reduce per-rank timing to min/max across ranks
      call MPI_ALLREDUCE(this%wt_plic,   this%wtmax_plic,   1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_plic,   this%wtmin_plic,   1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_plicnet,this%wtmax_plicnet,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_plicnet,this%wtmin_plicnet,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_polygon,this%wtmax_polygon,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_polygon,this%wtmin_polygon,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_visc,   this%wtmax_visc,   1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_visc,   this%wtmin_visc,   1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      ! Reset per-rank timing accumulators for next interval
      this%wt_plic=0.0_WP; this%wt_plicnet=0.0_WP; this%wt_polygon=0.0_WP; this%wt_visc=0.0_WP

   end subroutine get_info

   !> Print solver info to screen
   subroutine amrmpincomp_print(this)
      use messager, only: log
      implicit none
      class(amrmpincomp), intent(in) :: this
      call log("Incompressible Multiphase solver: "//trim(this%name))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrmpincomp_print

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrmpincomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%U,'U')
      call io%add_data(this%V,'V')
      call io%add_data(this%W,'W')
      call io%add_data(this%P,'P')
      call io%add_data(this%VF,'VF')
      call io%add_data(this%Cliq,'Cliq')
      call io%add_data(this%Cgas,'Cgas')
      call io%add_data(this%PLIC,'PLIC')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrmpincomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      integer :: lvl
      call io%read_data(dirname,this%U,'U')
      call io%read_data(dirname,this%V,'V')
      call io%read_data(dirname,this%W,'W')
      call io%read_data(dirname,this%P,'P')
      call io%read_data(dirname,this%VF,'VF')
      call io%read_data(dirname,this%Cliq,'Cliq')
      call io%read_data(dirname,this%Cgas,'Cgas')
      call io%read_data(dirname,this%PLIC,'PLIC')
      ! Fill ghost cells
      call this%U%fill(time)
      call this%V%fill(time)
      call this%W%fill(time)
      do lvl=0,this%amr%clvl()
         call this%fill_moments_lvl(lvl,time)
         call this%fill_plic_lvl(lvl,time)
      end do
      ! Rebuild polygons from restored PLIC
      call this%build_polygons()
   end subroutine restore_checkpoint

end module amrmpincomp_class



module amrscalar_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrflux_class,    only: amrflux
   use amrsolver_class,  only: amrsolver
   use amrex_amr_module, only: amrex_multifab,amrex_boxarray,amrex_distromap
   implicit none
   private

   ! Expose type/constructor/methods
   public :: amrscalar

   !> Constant density scalar solver object definition
   type, extends(amrsolver) :: amrscalar

      ! This is our amrgrid
      class(amrgrid), pointer :: amr=>null()

      ! Scalar variable definition
      integer :: nscalar
      character(len=str_medium), dimension(:), allocatable :: SCname

      ! Scalar data containers
      type(amrdata) :: SC
      type(amrdata) :: SCold

      ! Flux register for interlevel conservation
      type(amrflux) :: flux

      ! Overlap size
      integer :: nover=2

      ! Monitoring quantities
      real(WP), dimension(:), allocatable :: SCmax,SCmin,SCint

   contains
      ! Basic procedures
      procedure :: initialize
      procedure :: finalize
      procedure :: get_info
      ! Physics procedures
      procedure :: get_dSCdt
      procedure :: copy2old
      procedure :: reflux_avg_lvl
      procedure :: reflux_avg
      ! Fill with solver context (for coupled BCs)
      procedure :: fill
      ! Implement amrsolver abstract methods
      procedure :: fillbc => amrscalar_fillbc
      procedure :: on_regrid => amrscalar_on_regrid
      procedure :: register_checkpoint => amrscalar_register_checkpoint
      procedure :: restore_checkpoint => amrscalar_restore_checkpoint
   end type amrscalar

contains


   !> Initialization for amrscalar solver
   subroutine initialize(this,amr,nscalar,name)
      use iso_c_binding,    only: c_loc
      use messager,         only: die
      implicit none
      class(amrscalar), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      integer, intent(in) :: nscalar
      character(len=*), optional :: name

      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))

      ! Point to amrgrid object
      this%amr=>amr

      ! Set the number of scalars
      this%nscalar=nscalar
      if (this%nscalar.le.0) call die('[amrscalar initialize] At least 1 scalar is required')

      ! Initialize scalar names
      allocate(this%SCname(1:this%nscalar))
      this%SCname=''

      ! Initialize info storage
      allocate(this%SCmin(1:this%nscalar)); this%SCmin=+huge(1.0_WP)
      allocate(this%SCmax(1:this%nscalar)); this%SCmax=-huge(1.0_WP)
      allocate(this%SCint(1:this%nscalar)); this%SCint= 0.0_WP

      ! Initialize SC and SCold using amrdata%initialize
      call this%SC%initialize(amr, name='SC', ncomp=this%nscalar, ng=0)
      call this%SCold%initialize(amr, name='SCold', ncomp=this%nscalar, ng=0)

      ! Initialize flux register
      call this%flux%initialize(amr%maxlvl, name='SC_flux', ncomp=this%nscalar)

      ! Register callbacks with amrgrid (pass self as context)
      ! Note: Need select type since c_loc requires non-polymorphic type
      select type (this)
       type is (amrscalar)
         call this%amr%add_on_init(on_init_level,c_loc(this))
         call this%amr%add_on_coarse(on_coarse_level,c_loc(this))
         call this%amr%add_on_remake(on_remake_level,c_loc(this))
         call this%amr%add_on_clear(on_clear_level,c_loc(this))
         call this%amr%add_postregrid(on_postregrid,c_loc(this))
      end select

   end subroutine initialize


   !> Callback: create level from scratch
   subroutine on_init_level(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx,this)
      ! Initialize SC and SCold (calls define + user_on_init if set)
      call this%SC%on_init(lvl,ba,dm,this%amr%geom(lvl))
      call this%SCold%on_init(lvl,ba,dm,this%amr%geom(lvl))
      ! Build flux register for fine levels
      if (lvl.ge.1) call this%flux%build_level(lvl,ba,dm,this%amr%rref(lvl-1))
   end subroutine on_init_level


   !> Callback: create level from coarse
   subroutine on_coarse_level(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx,this)
      ! Initialize SC (define + fill_from_coarse or user override)
      call this%SC%on_coarse(lvl, ba, dm, time)
      ! SCold only needs to be defined (no fill needed)
      call this%SCold%define(lvl, ba, dm)
      ! Build flux register for fine levels
      if (lvl.ge.1) call this%flux%build_level(lvl,ba,dm,this%amr%rref(lvl-1))
   end subroutine on_coarse_level



   !> Callback: remake level (after regrid)
   subroutine on_remake_level(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx,this)
      ! Remake SC (handles fill from old data internally)
      call this%SC%on_remake(lvl, ba, dm, time)
      ! SCold just needs new geometry (define clears if needed)
      call this%SCold%define(lvl, ba, dm)
      ! Rebuild flux register
      if (lvl.ge.1) then
         call this%flux%destroy_level(lvl)
         call this%flux%build_level(lvl, ba, dm, this%amr%rref(lvl-1))
      end if
   end subroutine on_remake_level


   !> Callback: clear level
   subroutine on_clear_level(ctx,lvl)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx,this)
      call this%SC%on_clear(lvl)
      call this%SCold%on_clear(lvl)
      if (lvl.ge.1) call this%flux%destroy_level(lvl)
   end subroutine on_clear_level


   !> Callback: post-regrid (called once after all levels are processed)
   subroutine on_postregrid(ctx)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx,this)
      ! Average down to ensure level consistency after regrid
      call this%SC%average_down()
   end subroutine on_postregrid


   !> Finalization for amrscalar solver
   impure elemental subroutine finalize(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      ! Destroy containers
      call this%SC%destroy()
      call this%SCold%destroy()
      call this%flux%destroy()
      ! Deallocate arrays
      if (allocated(this%SCmax)) deallocate(this%SCmax)
      if (allocated(this%SCmin)) deallocate(this%SCmin)
      if (allocated(this%SCint)) deallocate(this%SCint)
      if (allocated(this%SCname)) deallocate(this%SCname)
      nullify(this%amr)
   end subroutine finalize


   !> Calculate various information on our amrscalar object
   subroutine get_info(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl,nsc
      this%SCmin=+huge(1.0_WP)
      this%SCmax=-huge(1.0_WP)
      this%SCint= 0.0_WP
      do nsc=1,this%nscalar
         do lvl=0,this%amr%clvl()
            this%SCmin(nsc)=min(this%SCmin(nsc),this%SC%mf(lvl)%min(comp=nsc))
            this%SCmax(nsc)=max(this%SCmax(nsc),this%SC%mf(lvl)%max(comp=nsc))
         end do
         this%SCint(nsc)=this%SC%mf(0)%sum(comp=nsc)*(this%amr%dx(0)*this%amr%dy(0)*this%amr%dz(0))/this%amr%vol
      end do
   end subroutine get_info


   !> Calculate dSC/dt at level (lvl)
   subroutine get_dSCdt(this,lvl,dSCdt,SC,U,V,W)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: dSCdt
      type(amrex_multifab), intent(in) :: SC,U,V,W
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab), dimension(3) :: flx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: rhs,pSC,FX,FY,FZ,pU,pV,pW
      integer :: i,j,k,nsc

      ! Prepare flux multifabs
      call this%amr%mfab_build(lvl=lvl,mfab=flx(1),ncomp=this%nscalar,nover=0,atface=[.true. ,.false.,.false.])
      call this%amr%mfab_build(lvl=lvl,mfab=flx(2),ncomp=this%nscalar,nover=0,atface=[.false.,.true. ,.false.])
      call this%amr%mfab_build(lvl=lvl,mfab=flx(3),ncomp=this%nscalar,nover=0,atface=[.false.,.false.,.true. ])

      ! Loop over boxes
      call this%amr%mfiter_build(lvl,mfi)
      do while(mfi%next())
         bx=mfi%tilebox()
         pSC=>SC%dataptr(mfi)
         rhs=>dSCdt%dataptr(mfi)
         pU=>U%dataptr(mfi)
         pV=>V%dataptr(mfi)
         pW=>W%dataptr(mfi)
         FX=>flx(1)%dataptr(mfi)
         FY=>flx(2)%dataptr(mfi)
         FZ=>flx(3)%dataptr(mfi)

         do nsc=1,this%nscalar
            ! X fluxes
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)+1
                     FX(i,j,k,nsc)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i-2,j,k,nsc)+5.0_WP/6.0_WP*pSC(i-1,j,k,nsc)+2.0_WP/6.0_WP*pSC(i  ,j,k,nsc)) &
                     &             -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i-1,j,k,nsc)+5.0_WP/6.0_WP*pSC(i  ,j,k,nsc)-1.0_WP/6.0_WP*pSC(i+1,j,k,nsc))
                  end do
               end do
            end do
            ! Y fluxes
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)+1
                  do i=bx%lo(1),bx%hi(1)
                     FY(i,j,k,nsc)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i,j-2,k,nsc)+5.0_WP/6.0_WP*pSC(i,j-1,k,nsc)+2.0_WP/6.0_WP*pSC(i,j  ,k,nsc)) &
                     &             -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i,j-1,k,nsc)+5.0_WP/6.0_WP*pSC(i,j  ,k,nsc)-1.0_WP/6.0_WP*pSC(i,j+1,k,nsc))
                  end do
               end do
            end do
            ! Z fluxes
            do k=bx%lo(3),bx%hi(3)+1
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)
                     FZ(i,j,k,nsc)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i,j,k-2,nsc)+5.0_WP/6.0_WP*pSC(i,j,k-1,nsc)+2.0_WP/6.0_WP*pSC(i,j,k  ,nsc)) &
                     &             -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i,j,k-1,nsc)+5.0_WP/6.0_WP*pSC(i,j,k  ,nsc)-1.0_WP/6.0_WP*pSC(i,j,k+1,nsc))
                  end do
               end do
            end do
            ! RHS
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)
                     rhs(i,j,k,nsc)=(FX(i+1,j,k,nsc)-FX(i,j,k,nsc))/this%amr%dx(lvl)+&
                     &              (FY(i,j+1,k,nsc)-FY(i,j,k,nsc))/this%amr%dy(lvl)+&
                     &              (FZ(i,j,k+1,nsc)-FZ(i,j,k,nsc))/this%amr%dz(lvl)
                  end do
               end do
            end do
         end do
         ! Scale fluxes for refluxing
         FX=FX*this%amr%dy(lvl)*this%amr%dz(lvl)
         FY=FY*this%amr%dz(lvl)*this%amr%dx(lvl)
         FZ=FZ*this%amr%dx(lvl)*this%amr%dy(lvl)
      end do
      call this%amr%mfiter_destroy(mfi)

      ! Handle refluxing
      if (lvl.gt.0)               call this%flux%fr(lvl  )%fineadd (flx,-1.0_WP)
      if (lvl.lt.this%amr%clvl()) call this%flux%fr(lvl+1)%crseinit(flx,+1.0_WP)

      ! Cleanup
      call this%amr%mfab_destroy(flx(1))
      call this%amr%mfab_destroy(flx(2))
      call this%amr%mfab_destroy(flx(3))

   end subroutine get_dSCdt


   !> Copy SC in SCold
   subroutine copy2old(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl
      do lvl=0,this%amr%clvl()
         call this%SCold%mf(lvl)%copy(srcmf=this%SC%mf(lvl),srccomp=1,dstcomp=1,nc=this%nscalar,ng=0)
      end do
   end subroutine copy2old


   !> Perform refluxing and averaging at a given level
   subroutine reflux_avg_lvl(this,lvl,dt)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer, intent(in) :: lvl
      interface
         subroutine amrex_fi_fluxregister_reflux(fr,mf,scale,geom) bind(c)
            import :: c_ptr,WP
            type(c_ptr), value :: fr,mf,geom
            real(WP), value :: scale
         end subroutine amrex_fi_fluxregister_reflux
      end interface
      call amrex_fi_fluxregister_reflux(this%flux%fr(lvl+1)%p,this%SC%mf(lvl)%p,dt,this%amr%geom(lvl)%p)
      call this%SC%average_downto(lvl)
   end subroutine reflux_avg_lvl


   !> Perform successive refluxing and averaging at all levels
   subroutine reflux_avg(this,dt)
      implicit none
      class(amrscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: lvl
      do lvl=this%amr%clvl()-1,0,-1
         call this%reflux_avg_lvl(lvl,dt)
      end do
   end subroutine reflux_avg


   !> Implement amrsolver%fillbc - apply physical BCs to scalar data
   subroutine amrscalar_fillbc(this,mf_ptr,geom_ptr,time,scomp,ncomp)
      use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrscalar), intent(inout) :: this
      type(c_ptr), value :: mf_ptr
      type(c_ptr), value :: geom_ptr
      real(WP), value :: time
      integer, value :: scomp,ncomp
      type(amrex_geometry) :: geom
      type(amrex_multifab) :: mf
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer, dimension(4) :: plo,phi
      ! Skip if fully periodic
      if (all([this%amr%xper,this%amr%yper,this%amr%zper])) return
      ! Convert pointers
      geom=geom_ptr; mf=mf_ptr
      ! Loop over boxes and apply filcc
      call amrex_mfiter_build(mfi,mf)
      do while(mfi%next())
         p=>mf%dataptr(mfi)
         if (.not.geom%domain%contains(p)) then
            plo=lbound(p); phi=ubound(p)
            call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,&
            &   geom%get_physical_location(plo),this%SC%lo_bc,this%SC%hi_bc)
         end if
      end do
   end subroutine amrscalar_fillbc


   !> Implement amrsolver%on_regrid - called after regrid for each level
   subroutine amrscalar_on_regrid(this,lvl,time)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, value :: lvl
      real(WP), value :: time
      ! Default: nothing special needed - rebuild is done in callbacks
   end subroutine amrscalar_on_regrid


   !> Fill ghost cells for an amrdata field using solver context
   !> This allows the BC callback to access all solver fields (for coupled BCs)
   subroutine fill(this, data, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two
      implicit none
      class(amrscalar), target, intent(inout) :: this
      class(amrdata), intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      real(WP) :: t_old, t_new
      type(c_ptr) :: solver_ctx
      type(c_funptr) :: bc_dispatch_ptr
      t_old = time - 1.0e200_WP
      t_new = time
      ! Get solver context (use select type for c_loc compatibility)
      select type (this)
       type is (amrscalar)
         solver_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(dispatch_amrscalar_fillbc)
      ! Call appropriate FillPatch
      if (lvl .eq. 0) then
         call amrmfab_fillpatch_single(data%mf(0)%p, t_old, data%mf(0)%p, &
         &   t_new, data%mf(0)%p, data%amr%geom(0)%p, solver_ctx, bc_dispatch_ptr, &
         &   time, 1, 1, data%ncomp)
      else
         call amrmfab_fillpatch_two(data%mf(lvl)%p, t_old, data%mf(lvl-1)%p, &
         &   t_new, data%mf(lvl-1)%p, data%amr%geom(lvl-1)%p, &
         &   t_old, data%mf(lvl)%p, t_new, data%mf(lvl)%p, data%amr%geom(lvl)%p, &
         &   solver_ctx, bc_dispatch_ptr, time, 1, 1, data%ncomp, &
         &   data%amr%rref(lvl-1), data%interp, data%lo_bc, data%hi_bc, data%ncomp)
      end if
   end subroutine fill


   !> BC dispatch for solver-level fill - recovers amrscalar and calls fillbc
   subroutine dispatch_amrscalar_fillbc(solver_ctx, mf_ptr, geom_ptr, time, scomp, ncomp) bind(c)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_int
      implicit none
      type(c_ptr), value, intent(in) :: solver_ctx
      type(c_ptr), value, intent(in) :: mf_ptr
      type(c_ptr), value, intent(in) :: geom_ptr
      real(c_double), value, intent(in) :: time
      integer(c_int), value, intent(in) :: scomp
      integer(c_int), value, intent(in) :: ncomp
      type(amrscalar), pointer :: sc
      call c_f_pointer(solver_ctx, sc)
      call sc%fillbc(mf_ptr, geom_ptr, real(time, WP), int(scomp), int(ncomp))
   end subroutine dispatch_amrscalar_fillbc


   !> Register fields for checkpoint with amrio
   subroutine amrscalar_register_checkpoint(this, io)
      use amrio_class, only: amrio
      implicit none
      class(amrscalar), intent(inout) :: this
      class(amrio), intent(inout) :: io
      ! Register SC and SCold containers (each contains all nscalar components)
      call io%add_data(this%SC, 'SC')
      call io%add_data(this%SCold, 'SCold')
   end subroutine amrscalar_register_checkpoint


   !> Restore fields from checkpoint
   subroutine amrscalar_restore_checkpoint(this, io, dirname)
      use amrio_class, only: amrio
      implicit none
      class(amrscalar), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      ! Restore SC and SCold containers
      call io%read_data(dirname, this%SC, 'SC')
      call io%read_data(dirname, this%SCold, 'SCold')
   end subroutine amrscalar_restore_checkpoint


end module amrscalar_class

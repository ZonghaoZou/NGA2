!> AMR Falling Cylinder - Free-falling 3D cylinder in quiescent gas
!> Periodic in X/Z, Dirichlet inflow at Y-, Foextrap outflow at Y+
!> PID-controlled moving domain to keep cylinder centered
module simulation
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrmpinc_class,    only: amrmpinc
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   use amrio_class,       only: amrio
   use string,            only: str_medium
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! Time integration
   type(timetracker) :: time

   ! Solver data
   type(amrmpinc), target :: fs
   type(amrdata) :: resUVW,Umag

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Re_tag=huge(1.0_WP)

   ! Monitoring
   type(monitor) :: mfile,cflfile,gridfile,dropfile

   ! Restart data
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   ! Physical parameters
   real(WP) :: viscL_mol,viscG_mol
   real(WP), dimension(3) :: gravity

   ! Moving domain (PID controller)
   logical :: moving_domain
   real(WP) :: Ycent,Vfall                  !< Droplet centroid Y and fall velocity
   real(WP) :: Yin,Vin,Vin_old              !< Inflow velocity for PID
   real(WP) :: Vfall_ref,Ycent_ref          !< PID reference values
   real(WP) :: PID_G,PID_ti                 !< PID gains

   ! Cylinder parameters
   real(WP) :: cyl_radius
   real(WP), dimension(3) :: cyl_center

contains

   !> Levelset function for 3D cylinder (axis along Z)
   function cylinder_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=cyl_radius-sqrt((xyz(1)-cyl_center(1))**2+(xyz(2)-cyl_center(2))**2)
   end function cylinder_levelset

   !> Compute viscosity
   subroutine get_viscosity()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVisc
      real(WP), parameter :: myeps=1.0e-15_WP
      ! Loop over levels
      do lvl=0,amr%clvl()
         ! Loop over domain
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Use harmonic averaging
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(viscL_mol,myeps)+(1.0_WP-pVF(i,j,k,1))/max(viscG_mol,myeps))
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosity

   !> Add gravity + frame acceleration source to velocity RHS (resUVW)
   !> Since get_dUVWdt already divides by density, this adds acceleration directly
   subroutine add_gravity_src(grav)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      real(WP), dimension(3), intent(in) :: grav
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRes
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pRes=>resUVW%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pRes(i,j,k,1)=pRes(i,j,k,1)+grav(1)
               pRes(i,j,k,2)=pRes(i,j,k,2)+grav(2)
               pRes(i,j,k,3)=pRes(i,j,k,3)+grav(3)
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine add_gravity_src

   !> Compute droplet centroid (Y) and fall velocity using AMR VF and UVW
   subroutine compute_droplet_stats()
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface, only: amrmask_make_fine
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      integer :: lvl,i,j,k,ierr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pUVW
      integer, dimension(:,:,:,:), contiguous, pointer :: pMask
      type(amrex_imultifab) :: mask
      real(WP) :: myYcent,myVfall,myvol,drop_vol
      real(WP) :: dy,cell_vol
      myYcent=0.0_WP; myVfall=0.0_WP; myvol=0.0_WP
      do lvl=0,amr%clvl()
         dy=amr%dy(lvl); cell_vol=amr%cell_vol(lvl)
         ! Build fine mask to avoid double-counting cells covered by finer level
         if (lvl.lt.amr%clvl()) then
            call amrex_imultifab_build(mask,amr%ba(lvl),amr%dm(lvl),1,0)
            call amrmask_make_fine(mask,amr%ba(lvl+1),[amr%rrefx(lvl),amr%rrefy(lvl),amr%rrefz(lvl)],0,1)
         end if
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pUVW=>fs%UVW%mf(lvl)%dataptr(mfi)
            if (lvl.lt.amr%clvl()) pMask=>mask%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Skip cells covered by finer level
               if (lvl.lt.amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
               myYcent=myYcent+(amr%ylo+(real(j,WP)+0.5_WP)*dy)*pVF(i,j,k,1)*cell_vol
               myVfall=myVfall+pUVW(i,j,k,2)*pVF(i,j,k,1)*cell_vol
               myvol=myvol+pVF(i,j,k,1)*cell_vol
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
         if (lvl.lt.amr%clvl()) call amrex_imultifab_destroy(mask)
      end do
      call MPI_ALLREDUCE(myYcent,Ycent   ,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      call MPI_ALLREDUCE(myVfall,Vfall   ,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      call MPI_ALLREDUCE(myvol  ,drop_vol,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      if (drop_vol.gt.0.0_WP) then
         Ycent=Ycent/drop_vol
         Vfall=Vfall/drop_vol
      end if
   end subroutine compute_droplet_stats

   !> Tagger for this case based on velocity gradient magnitude
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrmpinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVW
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,gradU_mag,Re_cell
      real(WP), dimension(3,3) :: gradU
      integer :: i,j,k
      tags=tags_ptr
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         tagarr=>tags%dataPtr(mfi)
         pUVW=>solver%UVW%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! No refinement near the outflow (y+)
            if (solver%amr%ylo+(real(j,WP)+0.5_WP)*dy.gt.solver%amr%yhi-5.0_WP*dy*2**solver%amr%maxlvl/real(2**lvl,WP)) cycle
            ! Velocity gradient tensor
            gradU(1,1)=0.5_WP*dxi*(pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))
            gradU(2,1)=0.5_WP*dyi*(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
            gradU(3,1)=0.5_WP*dzi*(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))
            gradU(1,2)=0.5_WP*dxi*(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))
            gradU(2,2)=0.5_WP*dyi*(pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))
            gradU(3,2)=0.5_WP*dzi*(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
            gradU(1,3)=0.5_WP*dxi*(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
            gradU(2,3)=0.5_WP*dyi*(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))
            gradU(3,3)=0.5_WP*dzi*(pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
            gradU_mag=sqrt(sum(gradU**2))
            Re_cell=solver%rhoG*gradU_mag*solver%amr%min_meshsize(lvl)**2/viscG_mol
            if (Re_cell.gt.Re_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Dirichlet BC: PID-controlled inflow at Y-, zero tangential
   subroutine dirichlet_velocity(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrmpinc), intent(in) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), pointer, intent(inout) :: p
      integer :: i,j,k,ic
      if (size(p,4).eq.1) then; ic=1
      else; ic=merge(1,merge(2,3,comp.eq.'V'),comp.eq.'U')
      end if
      select case (face)
       case (3)  ! Inflow at Y- (face=3 for lo_bc(2))
         select case (comp)
          case ('V')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,ic)=Vin
            end do; end do; end do
          case ('U','W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,ic)=0.0_WP
            end do; end do; end do
         end select
      end select
   end subroutine dirichlet_velocity

   !> User-provided initialization for drop
   subroutine cylinder_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpinc_class, only: VFlo
      class(amrmpinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pU,pUVW
      real(WP), dimension(3) :: BL,BG  ! Dummy barycenters
      real(WP) :: dx,dy,dz,VF
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         pUVW=>solver%UVW%mf(lvl)%dataptr(mfi)
         pU=>solver%U%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute VF and barycenters from levelset
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=cylinder_levelset,time=time,level=nref,VFlo=VFlo,VF=VF,BL=BL,BG=BG)
            ! Store volume fraction
            pVF(i,j,k,1)=VF
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine cylinder_init

   !> Initialization hook
   subroutine simulation_init()
      use param, only: param_read
      implicit none
      
      ! Create amrgrid
      create_amrgrid: block
         real(WP) :: Lx,Ly,Lz
         amr%name='amrcylinder'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         call param_read('Lz',Lz)
         amr%xlo=-0.5_WP*Lx; amr%xhi=+0.5_WP*Lx
         amr%ylo=-0.5_WP*Ly; amr%yhi=+0.5_WP*Ly
         amr%zlo=-0.5_WP*Lz; amr%zhi=+0.5_WP*Lz
         call param_read('Max level',amr%maxlvl)
         ! Periodic in X, non-periodic in Y
         amr%xper=.true.; amr%yper=.false.; amr%zper=.true.
         amr%nbloc=4
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart/saves
      handle_restart: block
         integer :: restart_step
         call io%initialize(amr=amr,nfiles=1)
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time integration
      initialize_time: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
      end block initialize_time
      
      ! Read physical parameters
      read_physics: block
         call param_read('Gravity',gravity)
         call param_read('Liquid dynamic viscosity',viscL_mol)
         call param_read('Gas dynamic viscosity',viscG_mol)
         call param_read('Moving domain',moving_domain,default=.false.)
      end block read_physics

      ! Set up cylinder parameters
      setup_cylinder: block
         call param_read('Droplet diameter',cyl_radius); cyl_radius=cyl_radius*0.5_WP
         cyl_center=[0.0_WP,0.0_WP,0.0_WP]  ! Domain center
      end block setup_cylinder

      ! Create flow solver
      create_flow_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrdata_class,    only: amrex_interp_face_linear
         use amrmpinc_class,   only: BC_GAS
         use amrmg_class,      only: amrmg_outer_pcg_mlmg
         call fs%initialize(amr,name='fallingcyl')
         fs%user_mpinc_init=>cylinder_init
         ! Set densities
         call param_read('Liquid density',fs%rhoL)
         call param_read('Gas density',fs%rhoG)
         ! Set surface tension
         call param_read('Surface tension coefficient',fs%sigma)
         ! Set pressure convergence
         fs%psolver%outer_solver=amrmg_outer_pcg_mlmg
         fs%psolver%tol_rel=1.0e-5_WP
         ! Set boundary conditions: Y- = inflow (Dirichlet), Y+ = outflow (foextrap)
         fs%lo_bc(2)=BC_GAS
         ! Y-direction BCs for cell-centered velocity (UVW)
         fs%UVW%lo_bc(2,:)=amrex_bc_ext_dir
         fs%UVW%hi_bc(2,:)=amrex_bc_foextrap
         ! Y-direction BCs for face velocities
         fs%U%lo_bc(2,1)=amrex_bc_ext_dir
         fs%V%lo_bc(2,1)=amrex_bc_ext_dir
         fs%W%lo_bc(2,1)=amrex_bc_ext_dir
         fs%U%hi_bc(2,1)=amrex_bc_foextrap
         fs%V%hi_bc(2,1)=amrex_bc_foextrap
         fs%W%hi_bc(2,1)=amrex_bc_foextrap
         ! Set user BC callback
         fs%user_mpinc_bc=>dirichlet_velocity
      end block create_flow_solver

      ! Create workspace array
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call resUVW%initialize(amr,name='resUVW',ncomp=3,ng=0,interp=amrex_interp_none); call resUVW%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=amrex_interp_none); call Umag%register()
      end block create_workspace

      ! Initialize PID controller
      init_pid: block
         Vin=0.0_WP; Vin_old=0.0_WP
         Ycent_ref=cyl_center(2)
         Vfall_ref=0.0_WP
         PID_G=0.25_WP
         PID_ti=time%dtmax
      end block init_pid

      ! Initialize regridding
      init_regridding: block
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         fs%user_mpinc_tagging=>my_tagger
         call param_read('Tagging Reynolds',Re_tag)
         if (restarted) then
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            call amr%init_from_scratch(time=time%t)
            call fs%build_plic(time%t)
            call fs%build_subVF()
            call fs%interp_vel_to_face()
            call fs%average_down_velocity()
            call fs%fill_velocity(time=time%t)
         end if
         ! Set viscosity (molecular only, no SGS)
         call get_viscosity()
         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         call fs%register_checkpoint(io)
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize visualization
      create_visualization: block
         call viz%initialize(amr,'amrfallingcyl',use_hdf5=.false.)
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_surfmesh(fs%smesh,'plic')
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_visualization

      ! Create monitor
      create_monitor: block
         call fs%get_info()
         call fs%get_cfl(time%dt,time%cfl)
         call compute_droplet_stats()
         ! Simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(fs%CFL,'CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%psolver%res,'Pressure residual')
         call mfile%add_column(fs%psolver%niter,'Pressure iterations')
         call mfile%add_column(fs%divmax,'Divergence')
         call mfile%write()
         ! CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLst,'CFLst')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
         ! Droplet monitor
         dropfile=monitor(amRoot=amr%amRoot,name='droplet')
         call dropfile%add_column(time%n,'Timestep')
         call dropfile%add_column(time%t,'Time')
         call dropfile%add_column(Ycent,'Y centroid')
         call dropfile%add_column(Vfall,'Fall velocity')
         call dropfile%add_column(Vin,'Inflow velocity')
         call dropfile%write()
      end block create_monitor

   end subroutine simulation_init

   !> Run the simulation
   subroutine simulation_run()
      implicit none
      real(WP), dimension(3) :: grav_eff

      ! Time integration loop
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! PID controller: compute inflow velocity for moving domain
         ! Uses Ycent/Vfall from end of previous timestep (or init)
         if (moving_domain) then
            Vin_old=Vin
            Vin=PID_G*((Vfall_ref-Vfall)+(Ycent_ref-Ycent)/PID_ti)
         end if

         ! Store old interface and velocities
         call fs%store_old()

         ! Compute effective gravity (physical + frame acceleration)
         grav_eff=gravity
         if (moving_domain) grav_eff(2)=grav_eff(2)+(Vin-Vin_old)/time%dt

         ! Sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity: U^{mid} = 0.5*(U + Uold)
            call fs%UVW%lincomb(a=0.5_WP,src1=fs%UVWold,b=0.5_WP,src2=fs%UVW)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Increment velocity with advection+viscous terms
            call fs%get_dUVWdt(dUVWdt=resUVW,dt=time%dt,time=time%t)

            ! Add gravity + frame acceleration source term
            call add_gravity_src(grav_eff)

            ! Update velocity
            call fs%UVW%lincomb(a=1.0_WP,src1=fs%UVWold,b=time%dt,src2=resUVW)
            call fs%UVW%average_down(); call fs%UVW%fill(time%t)

            ! Rebuild PLIC and sub-cell VF
            call fs%build_plic(time%t)
            call fs%build_subVF()

            ! Interpolate velocity to the faces
            call fs%interp_vel_to_face()

            ! Increment both velocities with current pressure term
            call fs%correct_both_velocities(scale=time%dt,phi=fs%P)

            ! Add surface tension to both velocities
            call fs%add_surface_tension(scale=time%dt)

            ! Average down and fill ghosts
            call fs%UVW%average_down(); call fs%UVW%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Correct outflow for mass conservation
            call fs%correct_outflow()

            ! Prepare and solve pressure Poisson
            call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
            call fs%prepare_psolver()
            call fs%psolver%solve(rhs=fs%div)

            ! Correct both velocities with pressure increment
            call fs%correct_both_velocities(scale=time%dt)

            ! Add pressure increment
            call fs%P%add(src=fs%psolver%sol)

            ! Average down and fill ghosts
            call fs%UVW%average_down(); call fs%UVW%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Update viscosity (molecular only)
         call get_viscosity()

         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)

         ! Compute droplet stats for monitoring
         call compute_droplet_stats()

         ! Monitor output
         call fs%get_info()
         call mfile%write()
         call cflfile%write()
         call dropfile%write()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/cyl_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if
         
      end do

   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      implicit none
      call time%finalize()
      call amr%finalize()
      call regrid_evt%finalize()
      call fs%finalize()
      call resUVW%finalize()
      call Umag%finalize()
      call viz%finalize()
      call viz_evt%finalize()
      call save_evt%finalize()
      call io%finalize()
      call mfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
      call dropfile%finalize()
   end subroutine simulation_final

end module simulation

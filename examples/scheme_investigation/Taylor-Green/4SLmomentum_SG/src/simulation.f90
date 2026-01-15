!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpcons_class,      only: tpcons
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   public :: simulation_init,simulation_run,simulation_final
   
   !> Flow solver objects
   type(hypre_str),   public :: ps     !< Structured Hypre linear solver for pressure
   type(ddadi),       public :: vs     !< DDADI solver for velocity
   type(tpcons),      public :: fs     !< Two-phase flow solver
   type(vfs),         public :: vf     !< Volume fraction solver
   type(timetracker), public :: time   !< Time info
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh     !< Surface mesh for interface
   type(ensight)  :: ens_out   !< Ensight output for flow variables
   type(event)    :: ens_evt   !< Event trigger for Ensight output
   
   !> Monitoring files
   type(monitor) :: mfile      !< General simulation monitoring
   type(monitor) :: cflfile    !< CFL monitoring
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities (based on Umid)
   real(WP), dimension(:,:,:), allocatable :: psi
   
   real(WP) :: KE,VFInt,rhoInt,rhoUInt,rhoVInt,rhoWInt
   real(WP) :: KE_l,KE_g,rhoU_l,rhoV_l,rhoW_l,rhoU_g,rhoV_g,rhoW_g
   character(len=20) :: filename='output.csv'

contains
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      real(WP), dimension(3) :: center
      real(WP) :: radius,depth
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(psi (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))         
      end block allocate_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: plicnet,VFhi,VFlo,flux,lvira
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=plicnet,transport_method=flux,name='VOF')
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_square,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use mathtools, only: pi
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Two-phase NS')
         ! Read in parameters
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver
      
      
      ! Generate initial conditions for velocity
      initialize_velocity: block
         use mathtools,       only: Pi,twoPi
         use vfs_class,       only: VFlo
         use random,          only: random_normal
         integer :: i,j,k
         ! Initialize density
         resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_faceRHO(rho=resU)
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_ 
                  psi(i,j,k)=sin(cfg%x(i))*sin(cfg%y(j))*cos(cfg%z(k))
               end do
            end do
         end do
         call cfg%sync(psi)

         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_ 
                  fs%U(i,j,k)= sum(fs%divp_y(:,i,j,k)*psi(i,j:j+1,k))
                  fs%V(i,j,k)=-sum(fs%divp_x(:,i,j,k)*psi(i:i+1,j,k))
               end do
            end do
         end do
         fs%W=0.0_WP
         call cfg%sync(fs%U)
         call cfg%sync(fs%V)
         ! Calculate cell-centered velocities and divergence
         call fs%get_div()
         if (cfg%amRoot) then
            open(unit=10, file=filename, status="replace", action="write")
            write(10, '(A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24)') &
            & 'time', 'timestep', 'KE', 'KE_l', 'KE_g', 'rhoU', 'rhoV', 'rhoW', 'rhoU_l', 'rhoV_l', 'rhoW_l','rhoU_g', 'rhoV_g', 'rhoW_g','VF','rho'
            close(unit=10)
         end if
         call get_KE()
      end block initialize_velocity
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='FallingDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
   contains
      
      !> Function that defines a level set function for a falling drop problem
      function levelset_square(xyz,t) result(G)
         use mathtools, only: pi
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=      1.0_WP-xyz(1)
         G=min(G,1.0_WP+xyz(1))
      end function levelset_square
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpcons_class, only: arithmetic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         vf%bandold=vf%band
         
         ! Remember old velocities and sRHOs
         fs%Uold=fs%U; fs%RHOXold=fs%RHOX
         fs%Vold=fs%V; fs%RHOYold=fs%RHOY
         fs%Wold=fs%W; fs%RHOZold=fs%RHOZ
         
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho_l=fs%rho_l,rho_g=fs%rho_g)
         ! Update sqrt(face density) and momentum vector
         resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_faceRHO(rho=resU)
         fs%rhoU=fs%rho_l*vf%UFl(1,:,:,:)+fs%rho_g*vf%UFg(1,:,:,:)
         fs%rhoV=fs%rho_l*vf%UFl(2,:,:,:)+fs%rho_g*vf%UFg(2,:,:,:)
         fs%rhoW=fs%rho_l*vf%UFl(3,:,:,:)+fs%rho_g*vf%UFg(3,:,:,:)

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW,vf,1)
            ! Assemble explicit residual
            resU=-2.0_WP*fs%RHOX*fs%U+(fs%RHOX+fs%RHOXold)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%RHOY*fs%V+(fs%RHOY+fs%RHOYold)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%RHOZ*fs%W+(fs%RHOZ+fs%RHOZold)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW,vf,1)
            
            ! Compute predictor U
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%RHOX
            fs%V=fs%V-time%dt*resV/fs%RHOY
            fs%W=fs%W-time%dt*resW/fs%RHOZ
         
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do

         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)  
         call fs%get_div()
         
         call get_KE()
         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   subroutine get_KE
      implicit none
      integer :: i,j,k,ierr
      real(WP), dimension(:,:,:), allocatable:: Utmp !U direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: Vtmp !V direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: Wtmp !W direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: FX   !U direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: FY   !V direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: FZ   !W direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: Stmp !Summed properties for integration
      allocate(Utmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Utmp=0.0_WP
      allocate(Vtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Vtmp=0.0_WP
      allocate(Wtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Wtmp=0.0_WP
      allocate(Stmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Stmp=0.0_WP
      allocate(FX(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FX=0.0_WP
      allocate(FY(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FY=0.0_WP
      allocate(FZ(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FZ=0.0_WP
      calculate_KE: block
         Utmp=fs%RHOX*fs%U*fs%U/2.0_WP
         Vtmp=fs%RHOY*fs%V*fs%V/2.0_WP
         Wtmp=fs%RHOZ*fs%W*fs%W/2.0_WP
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))+&
                  &           sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))+&
                  &           sum(fs%itpw_z(:,i,j,k)*Wtmp(i,j,k:k+1))
               end do
            end do
         end do
         call cfg%integrate(Stmp,integral=KE)
      end block calculate_KE

      Utmp=0.0_WP;Vtmp=0.0_WP;Wtmp=0.0_WP
      calculate_momentumconservation: block
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*(fs%RHOX(i:i+1,j,k)*fs%U(i:i+1,j,k)))
                  Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*(fs%RHOY(i,j:j+1,k)*fs%V(i,j:j+1,k)))
                  Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*(fs%RHOZ(i,j,k:k+1)*fs%W(i,j,k:k+1)))
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoUInt)
         call fs%cfg%integrate(Vtmp,integral=rhoVInt)
         call fs%cfg%integrate(Wtmp,integral=rhoWInt)
      end block calculate_momentumconservation 

      ! Get staggered liquid face density
      FX=0.0_WP; FY=0.0_WP; FZ=0.0_WP; Stmp=0.0_WP; resU=fs%rho_l*vf%VF
      do k=fs%cfg%kmino_  ,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_  ,fs%cfg%jmaxo_
            do i=fs%cfg%imino_+1,fs%cfg%imaxo_
               FX(i,j,k)=sum(fs%itpr_x(:,i,j,k)*resU(i-1:i,j,k))
            end do
         end do
      end do
      do k=fs%cfg%kmino_  ,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_+1,fs%cfg%jmaxo_
            do i=fs%cfg%imino_  ,fs%cfg%imaxo_
               FY(i,j,k)=sum(fs%itpr_y(:,i,j,k)*resU(i,j-1:j,k))
            end do
         end do
      end do
      do k=fs%cfg%kmino_+1,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_  ,fs%cfg%jmaxo_
            do i=fs%cfg%imino_  ,fs%cfg%imaxo_
               FZ(i,j,k)=sum(fs%itpr_z(:,i,j,k)*resU(i,j,k-1:k))
            end do
         end do
      end do
      ! Handle non-periodic borders
      if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.1) FX(fs%cfg%imino,:,:)=resU(fs%cfg%imino,:,:)
      if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.1) FY(:,fs%cfg%jmino,:)=resU(:,fs%cfg%jmino,:)
      if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.1) FZ(:,:,fs%cfg%kmino)=resU(:,:,fs%cfg%kmino)
      ! Synchronize boundaries
      call fs%cfg%sync(FX)
      call fs%cfg%sync(FY)
      call fs%cfg%sync(FZ)

      ! Get KE_l
      calculate_l_KE: block
         Utmp=FX*fs%U*fs%U/2.0_WP
         Vtmp=FY*fs%V*fs%V/2.0_WP
         Wtmp=FZ*fs%W*fs%W/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))+&
                  &           sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))+&
                  &           sum(fs%itpw_z(:,i,j,k)*Wtmp(i,j,k:k+1))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_l)
      end block calculate_l_KE


      Utmp=0.0_WP;Vtmp=0.0_WP;Wtmp=0.0_WP
      calculate_l_momentumconservation: block
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*(FX(i:i+1,j,k)*fs%U(i:i+1,j,k)))
                  Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*(FY(i,j:j+1,k)*fs%V(i,j:j+1,k)))
                  Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*(FZ(i,j,k:k+1)*fs%W(i,j,k:k+1)))
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoU_l)
         call fs%cfg%integrate(Vtmp,integral=rhoV_l)
         call fs%cfg%integrate(Wtmp,integral=rhoW_l)
      end block calculate_l_momentumconservation 


      ! Get staggered gas face density
      FX=0.0_WP; FY=0.0_WP; FZ=0.0_WP; Stmp=0.0_WP; resU=fs%rho_g*(1.0_WP-vf%VF)
      do k=fs%cfg%kmino_  ,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_  ,fs%cfg%jmaxo_
            do i=fs%cfg%imino_+1,fs%cfg%imaxo_
               FX(i,j,k)=sum(fs%itpr_x(:,i,j,k)*resU(i-1:i,j,k))
            end do
         end do
      end do
      do k=fs%cfg%kmino_  ,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_+1,fs%cfg%jmaxo_
            do i=fs%cfg%imino_  ,fs%cfg%imaxo_
               FY(i,j,k)=sum(fs%itpr_y(:,i,j,k)*resU(i,j-1:j,k))
            end do
         end do
      end do
      do k=fs%cfg%kmino_+1,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_  ,fs%cfg%jmaxo_
            do i=fs%cfg%imino_  ,fs%cfg%imaxo_
               FZ(i,j,k)=sum(fs%itpr_z(:,i,j,k)*resU(i,j,k-1:k))
            end do
         end do
      end do
      ! Handle non-periodic borders
      if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.1) FX(fs%cfg%imino,:,:)=resU(fs%cfg%imino,:,:)
      if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.1) FY(:,fs%cfg%jmino,:)=resU(:,fs%cfg%jmino,:)
      if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.1) FZ(:,:,fs%cfg%kmino)=resU(:,:,fs%cfg%kmino)
      ! Synchronize boundaries
      call fs%cfg%sync(FX)
      call fs%cfg%sync(FY)
      call fs%cfg%sync(FZ)

      calculate_g_KE: block
         Utmp=FX*fs%U*fs%U/2.0_WP
         Vtmp=FY*fs%V*fs%V/2.0_WP
         Wtmp=FZ*fs%W*fs%W/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))+&
                  &           sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))+&
                  &           sum(fs%itpw_z(:,i,j,k)*Wtmp(i,j,k:k+1))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_g)
      end block calculate_g_KE


      Utmp=0.0_WP;Vtmp=0.0_WP;Wtmp=0.0_WP
      calculate_g_momentumconservation: block
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*(FX(i:i+1,j,k)*fs%U(i:i+1,j,k)))
                  Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*(FY(i,j:j+1,k)*fs%V(i,j:j+1,k)))
                  Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*(FZ(i,j,k:k+1)*fs%W(i,j,k:k+1)))
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoU_g)
         call fs%cfg%integrate(Vtmp,integral=rhoV_g)
         call fs%cfg%integrate(Wtmp,integral=rhoW_g)
      end block calculate_g_momentumconservation 

      call fs%cfg%integrate(vf%VF,integral=VFInt)
      call fs%cfg%integrate(fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF),integral=rhoInt)


       if (cfg%amRoot) then
          ! Open file dynamically with append mode
          open(unit=10, file=filename, status="unknown", position="append", action="write")
          ! Write data with 16-digit precision in CSV format
          write(10, '(F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16)')&
          time%t,time%dt,KE,KE_l,KE_g,rhoUInt,rhoVInt,rhoWInt,rhoU_l,rhoV_l,rhoW_l,rhoU_g,rhoV_g,rhoW_g,VFInt,rhoInt
          close(10)
       end if

   end subroutine get_KE
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
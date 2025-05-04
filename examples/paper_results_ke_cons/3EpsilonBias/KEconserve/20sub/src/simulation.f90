!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
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
   type(tpns),        public :: fs     !< Two-phase flow solver
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
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities (based on Umid)
   real(WP), dimension(:,:,:,:), allocatable :: vel               !< Other cell-centered velocity (based on U)
   real(WP), dimension(:,:,:), allocatable :: psi
   
   real(WP) :: dKEdt,KE,VFInt,rhoInt,rhoUInt,rhoVInt,rhoWInt
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
         allocate(vel (1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(psi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))         
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
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux,name='VOF')
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
         real(WP) :: epsilon
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Two-phase NS')
         call param_read('epsilon',epsilon)
         fs%theta=0.5_WP+epsilon

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
         ! vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)!,implicit_solver=vs)

         if (cfg%amRoot) then
            open(unit=10, file='output.csv', status="replace", action="write")
            write(10, '(A10,A10,A10,A10,A10,A10,A10,A10,A10)')'time','timestep','dKEdt', 'KE', 'VFint', 'rhoint','rhoUint','rhoVint','rhoWint'
            close(unit=10)
         end if

      end block create_flow_solver
      
      
      ! Generate initial conditions for velocity
      initialize_velocity: block
         use mathtools,       only: Pi,twoPi
         use vfs_class,       only: VFlo
         use random,          only: random_normal
         integer :: i,j,k,seed_size
         integer, allocatable, dimension(:)  :: seed
         real(WP) :: randx,randy,randamp,randtheta
         ! Initialize density
         resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_density(rho=resU)
         call random_seed(size=seed_size)
         allocate(seed(seed_size))
         seed = 13798
         call random_seed(put=seed)
         deallocate(seed)

         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_ 
                  randx=random_normal(m=0.0_WP,sd=1.0_WP)
                  randy=random_normal(m=0.0_WP,sd=1.0_WP)
                  randamp=random_normal(m=0.0_WP,sd=1.0_WP)
                  randtheta=random_normal(m=0.0_WP,sd=twoPi**2)
                  psi(i,j,k)=randamp*sin(randx*cfg%x(i)+randy*cfg%y(j)+randtheta)
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
         call cfg%sync(fs%U)
         call cfg%sync(fs%V)
         fs%W=0.0_WP
         fs%Umid=fs%U;fs%Vmid=fs%V;fs%Wmid=fs%W

         ! Calculate cell-centered velocities and divergence
         call fs%interp_velmid(Ui,Vi,Wi)
         call fs%interp_vel(vel(1,:,:,:),vel(2,:,:,:),vel(3,:,:,:))
         call fs%get_div()
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
         call ens_out%add_vector('othervel',vel(1,:,:,:),vel(2,:,:,:),vel(3,:,:,:))
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
         ! Add Y direction
         G=0.5_WP*pi+xyz(2) 
         G=min(G,0.5_WP*pi-xyz(2))
         ! Add X direction
         G=min(G,0.5_WP*pi-xyz(1))
         G=min(G,0.5_WP*pi+xyz(1))
      end function levelset_square
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocities and sRHOs
         fs%Uold=fs%U; fs%sRHOxold=fs%sRHOx
         fs%Vold=fs%V; fs%sRHOyold=fs%sRHOy
         fs%Wold=fs%W; fs%sRHOzold=fs%sRHOz
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! VOF equation ====================================================
            ! Advance VOF equation
            vf%VF=vf%VFold
            if (time%it.eq.time%itmax) then   
               call vf%advance(dt=time%dt,U=fs%Umid,V=fs%Vmid,W=fs%Wmid)
            else
               call vf%advance_tmp(dt=time%dt,U=fs%Umid,V=fs%Vmid,W=fs%Wmid)
            end if
            
            ! Update sqrt(face density) and momentum vector
            resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_density(rho=resU)
            fs%rhoU=fs%rho_l*vf%UFl(1,:,:,:)+fs%rho_g*vf%UFg(1,:,:,:)
            fs%rhoV=fs%rho_l*vf%UFl(2,:,:,:)+fs%rho_g*vf%UFg(2,:,:,:)
            fs%rhoW=fs%rho_l*vf%UFl(3,:,:,:)+fs%rho_g*vf%UFg(3,:,:,:)
            
            ! Prepare new staggered viscosity (at n+1)
            call fs%get_viscosity(vf=vf,strat=arithmetic_visc)

            ! Momentum equation ===============================================
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-(fs%U*fs%sRHOX**2-fs%Uold*fs%sRHOXold**2)+time%dt*resU
            resV=-(fs%V*fs%sRHOY**2-fs%Vold*fs%sRHOYold**2)+time%dt*resV
            resW=-(fs%W*fs%sRHOZ**2-fs%Wold*fs%sRHOZold**2)+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Compute predictor U
            fs%U=fs%U+resU
            fs%V=fs%V+resV
            fs%W=fs%W+resW
            
            ! Sync and apply boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Poisson equation ================================================
            ! Compute Umid from U and Uold
            call fs%get_Umid()
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct pressure and Umid
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            ! fs%P=fs%P+fs%psolv%sol
            fs%P=fs%psolv%sol
            fs%Umid=fs%Umid-time%dt*resU/((fs%sRHOX+fs%sRHOXold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOX)
            fs%Vmid=fs%Vmid-time%dt*resV/((fs%sRHOY+fs%sRHOYold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOY)
            fs%Wmid=fs%Wmid-time%dt*resW/((fs%sRHOZ+fs%sRHOZold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOZ)
            
            ! Regenerate U from Umid and Uold
            call fs%get_U()
            
            ! Increment sub-iteration counter =================================
            time%it=time%it+1
            
         end do

         call get_KE()
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_velmid(Ui,Vi,Wi)
         call fs%interp_vel(vel(1,:,:,:),vel(2,:,:,:),vel(3,:,:,:))
         call fs%get_div()
         
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
      use mpi_f08,    only: MPI_ALLREDUCE,MPI_MAX
      use parallel,   only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,ierr
      real(WP), dimension(:,:,:), allocatable:: Utmp !U direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: Vtmp !V direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: Wtmp !W direction temp field for calculating properties 
      real(WP), dimension(:,:,:), allocatable:: Stmp !Summed properties for integration
      character(len=20) :: filename
      real(WP) :: Umax,myUmax
      allocate(Utmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Stmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      if (time%t.gt.0.0_WP) then
         calculate_dKEdt: block
            Utmp=((fs%sRHOX**2)*fs%U**2-(fs%sRHOXold**2)*fs%Uold**2)/(2.0_WP*time%dt)
            Vtmp=((fs%sRHOY**2)*fs%V**2-(fs%sRHOYold**2)*fs%Vold**2)/(2.0_WP*time%dt)
            Wtmp=((fs%sRHOZ**2)*fs%W**2-(fs%sRHOZold**2)*fs%Wold**2)/(2.0_WP*time%dt)
            do k=cfg%kmin_,cfg%kmax_
               do j=cfg%jmin_,cfg%jmax_
                  do i=cfg%imin_,cfg%imax_
                     Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))+&
                     &           sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))+&
                     &           sum(fs%itpw_z(:,i,j,k)*Wtmp(i,j,k:k+1))
                  end do
               end do
            end do
            call cfg%integrate_without_VF(Stmp,integral=dKEdt)
         end block calculate_dKEdt


         calculate_KE: block
            Utmp=(fs%sRHOX**2)*fs%U*fs%U/2.0_WP
            Vtmp=(fs%sRHOY**2)*fs%V*fs%V/2.0_WP
            Wtmp=(fs%sRHOZ**2)*fs%W*fs%W/2.0_WP
            do k=cfg%kmin_,cfg%kmax_
               do j=cfg%jmin_,cfg%jmax_
                  do i=cfg%imin_,cfg%imax_
                     Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))+&
                     &           sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))+&
                     &           sum(fs%itpw_z(:,i,j,k)*Wtmp(i,j,k:k+1))
                  end do
               end do
            end do
            call cfg%integrate_without_VF(Stmp,integral=KE)
         end block calculate_KE

         call cfg%integrate_without_VF((vf%VF-vf%VFold)/time%dt,integral=VFInt)
         
         call cfg%integrate_without_VF((vf%VF*fs%rho_l+(1.0_WP-vf%VF)*fs%rho_g -(vf%VFold*fs%rho_l+(1.0_WP-vf%VFold)*fs%rho_g))/time%dt,integral=rhoInt)

         calculate_momentumconservation: block
            stmp=0.0_WP
            do k=cfg%kmin_,cfg%kmax_
               do j=cfg%jmin_,cfg%jmax_
                  do i=cfg%imin_,cfg%imax_
                     Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*((fs%sRHOX(i:i+1,j,k)**2)*fs%U(i:i+1,j,k)-(fs%sRHOXold(i:i+1,j,k)**2)*fs%Uold(i:i+1,j,k)))/time%dt
                     Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*((fs%sRHOY(i,j:j+1,k)**2)*fs%V(i,j:j+1,k)-(fs%sRHOYold(i,j:j+1,k)**2)*fs%Vold(i,j:j+1,k)))/time%dt
                     Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*((fs%sRHOZ(i,j,k:k+1)**2)*fs%W(i,j,k:k+1)-(fs%sRHOZold(i,j,k:k+1)**2)*fs%Wold(i,j,k:k+1)))/time%dt
                  end do
               end do
            end do
            call cfg%integrate_without_VF(Utmp,integral=rhoUInt)
            call cfg%integrate_without_VF(Vtmp,integral=rhoVInt)
            call cfg%integrate_without_VF(Wtmp,integral=rhoWInt)
         end block calculate_momentumconservation 

         filename='output.csv'
         if (cfg%amRoot) then
            ! Open file dynamically with append mode
            open(unit=10, file=filename, status="unknown", position="append", action="write")
            ! Write data with 16-digit precision in CSV format
            write(10, '(F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16)') time%t,time%dt,dKEdt,KE,VFInt,rhoInt,rhoUInt,rhoVInt,rhoWInt
            close(10)

         end if
      else 
         call cfg%integrate_without_VF(vf%VF,integral=VFInt)
         call cfg%integrate_without_VF(vf%VF*fs%rho_l+(1.0_WP-vf%VF)*fs%rho_g,integral=rhoInt)
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*((fs%sRHOX(i:i+1,j,k)**2)*fs%U(i:i+1,j,k)))
                  Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*((fs%sRHOY(i,j:j+1,k)**2)*fs%V(i,j:j+1,k)))
                  Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*((fs%sRHOZ(i,j,k:k+1)**2)*fs%W(i,j,k:k+1)))
               end do
            end do
         end do
         call cfg%integrate_without_VF(Utmp,integral=rhoUInt)
         call cfg%integrate_without_VF(Vtmp,integral=rhoVInt)
         call cfg%integrate_without_VF(Wtmp,integral=rhoWInt)
         myUmax=0.0_WP
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  myUmax=max(myUmax,sqrt(((fs%U(i,j,k)+fs%U(i+1,j,k))/2.0_WP)**2+((fs%V(i,j,k)+fs%V(i,j+1,k))/2.0_WP)**2+((fs%W(i,j,k)+fs%W(i,j,k+1))/2.0_WP)**2))
               end do
            end do
         end do
         
         ! Get the parallel max
         call MPI_ALLREDUCE(myUmax,Umax,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)

         if (cfg%amRoot) print *, VFInt,rhoInt,rhoUInt,rhoVInt,rhoWInt,Umax
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
!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use iterator_class,    only: iterator
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpcons_class,      only: tpcons
   use vfs_class,         only: vfs
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   public :: simulation_init,simulation_run,simulation_final,output_info
   
   !> Flow solver objects
   type(hypre_str),   public :: ps     !< Structured Hypre linear solver for pressure
   type(ddadi),       public :: vs     !< DDADI solver for velocity
   type(tpcons),      public :: fs     !< Two-phase conservative flow solver
   type(vfs),         public :: vf     !< Volume fraction solver
   type(timetracker), public :: time   !< Time info
   
   !> SGS modeling
   logical        :: use_sgs   !< Is an LES model used?
   type(sgsmodel) :: sgs       !< SGS model for eddy viscosity
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh     !< Surface mesh for interface
   type(ensight)  :: ens_out   !< Ensight output for flow variables
   type(event)    :: ens_evt   !< Event trigger for Ensight output
   
   !> Monitoring files
   type(monitor) :: mfile      !< General simulation monitoring
   type(monitor) :: cflfile    !< CFL monitoring
   
   !> Iterator for VOF removal
   type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
   integer        :: nlayer=4           !< Size of buffer layer for VOF removal
   real(WP)       :: vof_removed        !< Integral of VOF removed
   
   integer :: type_method
   logical :: falling,Implicit,STflag
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
   real(WP), dimension(:,:,:), allocatable :: indicator      !< Residuals
   
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
         allocate(indicator(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));indicator=0.0_WP
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
      
      if (cfg%amRoot) print *, "Here-2"
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: plicnet,VFhi,VFlo,flux,flux_storage
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=plicnet,transport_method=flux,name='VOF')
         ! vf%cons_correct=.false.
         ! Initialize to a droplet and a pool
         radius=0.5_WP
         call param_read('Droplet center',center)
         call param_read('Pool depth',depth)
         call param_read('Falling Drop Activated',falling)
         call param_read('Surface tension Activated',STflag)
         call param_read('Implicit',Implicit)
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_falling_drop,0.0_WP,amr_ref_lvl)
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
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         vof_removal_layer=iterator(cfg,'VOF removal',vof_removal_layer_locator)
         vof_removed=0.0_WP
      end block create_iterator
      
      ! Create a two-phase flow solver without bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         real(WP) :: Re,Fr,We,r,m
         call param_read('Method type',type_method)
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Two-phase NS')
         ! Read in adimensional parameters
         call param_read('Froude number',Fr)     ! Fr=U^2/(g*D)=1/g
         call param_read('Weber number',We)      ! We=rho_l*U^2*D/sigma=1/sigma
         call param_read('Reynolds number',Re)   ! Re=rho_l*U*D/mu_l=1/mu_l
         call param_read('Density ratio',r)      ! r=rho_l/rho_g=1/rho_g
         call param_read('Viscosity ratio',m)    ! m=mu_l/mu_g
         ! Assign fluid properties to each phase
         if (falling) then
            fs%gravity=[0.0_WP,-Fr**(-1.0_WP),0.0_WP]
            fs%sigma=We**(-1.0_WP)
            fs%rho_l=1.0_WP
            fs%rho_g=fs%rho_l/r
            fs%visc_l=Re**(-1.0_WP)
            fs%visc_g=fs%visc_l/m
         else 
            fs%gravity=[0.0_WP,0.0_WP,0.0_WP]
            fs%sigma=We**(-1.0_WP)
            fs%rho_l=1.0_WP
            fs%rho_g=fs%rho_l/r
            ! fs%rho_l=r
            ! fs%rho_g=1.0_WP
            fs%visc_l=0.0_WP
            fs%visc_g=0.0_WP
         end if
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=24
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         if (Implicit) then
            vs=ddadi(cfg=cfg,name='Velocity',nst=7)
            call fs%setup(pressure_solver=ps,implicit_solver=vs)
         else
            call fs%setup(pressure_solver=ps)
         end if
      end block create_flow_solver
      
      ! Generate initial conditions for velocity
      initialize_velocity: block
         use vfs_class, only: VFlo
         integer :: i,j,k
         ! Initialize density
         fs%rho=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF)
         if (falling) then
            ! Initialize velocity
            do k=cfg%kmin_,cfg%kmax_
               do j=cfg%jmin_,cfg%jmax_
                  do i=cfg%imin_,cfg%imax_
                     if (cfg%ym(j).gt.depth.and.vf%VF(i,j,k).gt.VFlo) fs%V(i,j,k)=-1.0_WP
                  end do
               end do
            end do
            
            ! Apply all other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%construction_correction(time%dt)
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            ! Corrector step for face velocities
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%Uf=fs%Uf-resU/fs%rho_x
            fs%Vf=fs%Vf-resV/fs%rho_y
            fs%Wf=fs%Wf-resW/fs%rho_z

            ! Corrector step for center velocities
            call fs%get_pgrad_wide(fs%psolv%sol,resU,resV,resW)
            fs%U=fs%U-resU/fs%rho
            fs%V=fs%V-resV/fs%rho
            fs%W=fs%W-resW/fs%rho
         else
            fs%V=-0.0_WP
            fs%U=-1.0_WP
            call fs%construction_correction(time%dt)
         end if
         ! Calculate cell-centered velocities and divergence
         call fs%get_div()
      end block initialize_velocity
      
      ! ! Create an LES model
      ! create_sgs: block
      !    call param_read('Use SGS model',use_sgs)
      !    if (use_sgs) then
      !       allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      !       sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      !    end if
      ! end block create_sgs
      
      
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
         call ens_out%add_vector('velocity',fs%U,fs%V,fs%W)
         call ens_out%add_vector('facevelocity',fs%Uf,fs%Vf,fs%Wf)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('indicator',indicator)
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
         call mfile%add_column(vof_removed,'VOF removed')
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
      function levelset_falling_drop(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         ! Create the droplet
         G=radius-sqrt(sum((xyz-center)**2))
         ! Add the pool
         if (falling) G=max(G,depth-xyz(2))
      end function levelset_falling_drop
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.ge.pg%jmax-nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc
      implicit none
      integer :: i,j,k
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF; fs%rhoold=fs%rho
         ! Remember old velocity and face densities
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! VOF advection
         call vf%advance(dt=time%dt,U=fs%Uf,V=fs%Vf,W=fs%Wf)
         ! Update face density and momentum vector
         fs%rho=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF)
         fs%rhoU=fs%rho_l*vf%UFl(1,:,:,:)+fs%rho_g*vf%UFg(1,:,:,:)
         fs%rhoV=fs%rho_l*vf%UFl(2,:,:,:)+fs%rho_g*vf%UFg(2,:,:,:)
         fs%rhoW=fs%rho_l*vf%UFl(3,:,:,:)+fs%rho_g*vf%UFg(3,:,:,:)
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
         ! ! Turbulence modeling
         ! if (use_sgs) then
         !    sgs_modeling: block
         !       use sgsmodel_class, only: vreman
         !       integer :: i,j,k
         !       resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF)
         !       call fs%get_gradU(gradU)
         !       call sgs%get_visc(type=vreman,dt=time%dt,rho=resU,gradu=gradU)
         !       do k=fs%cfg%kmino_+1,fs%cfg%kmaxo_
         !          do j=fs%cfg%jmino_+1,fs%cfg%jmaxo_
         !             do i=fs%cfg%imino_+1,fs%cfg%imaxo_
         !                fs%visc(i,j,k)   =fs%visc(i,j,k)   +sgs%visc(i,j,k)
         !                fs%visc_xy(i,j,k)=fs%visc_xy(i,j,k)+sum(fs%itp_xy(:,:,i,j,k)*sgs%visc(i-1:i,j-1:j,k))
         !                fs%visc_yz(i,j,k)=fs%visc_yz(i,j,k)+sum(fs%itp_yz(:,:,i,j,k)*sgs%visc(i,j-1:j,k-1:k))
         !                fs%visc_zx(i,j,k)=fs%visc_zx(i,j,k)+sum(fs%itp_xz(:,:,i,j,k)*sgs%visc(i-1:i,j,k-1:k))
         !             end do
         !          end do
         !       end do
         !    end block sgs_modeling
         ! end if
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
          
            ! Build mid-time velocity =========================================
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho*fs%U+(fs%rho+fs%rhoold)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho*fs%V+(fs%rho+fs%rhoold)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho*fs%W+(fs%rho+fs%rhoold)*fs%Wold+time%dt*resW
            
            ! call output_info()
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Compute predictor U
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Sync and apply boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%construction_correction(time%dt)
            call fs%get_div()
            if (STflag) call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Corrector step for face velocities
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%Uf=fs%Uf-time%dt*resU/fs%rho_x
            fs%Vf=fs%Vf-time%dt*resV/fs%rho_y
            fs%Wf=fs%Wf-time%dt*resW/fs%rho_z

            ! Corrector step for center velocities
            call fs%get_pgrad_wide(fs%psolv%sol,resU,resV,resW)
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            
            ! call output_info()
            ! Increment sub-iteration counter =================================
            time%it=time%it+1
         end do
         ! Recompute interpolated velocity and divergence
         ! call fs%interp_vel(Ui,Vi,Wi)
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
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW)
      if (use_sgs) deallocate(gradU)
      
   end subroutine simulation_final


   subroutine output_info
      implicit none
      integer :: i,j,k
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               ! if (abs(fs%U(i,j,k)).gt.10894_WP) then
               !    print  *, "This is the U velocity location", i,j,k,fs%U(i,j,k)
               ! end if


               ! if (abs(fs%U(i,j,k)).gt.3.066) then
               !    print  *, "This is the U velocity location", i,j,k,fs%U(i,j,k)
               ! end if
               if ((i.eq.74).and.(j.eq.68).and.(k.eq.1)) then
                  indicator(i,j,k)=1.0_WP
                  print*, fs%rhoU(i,j,k), fs%rho_l*vf%UFl(1,i,j,k)+fs%rho_g*vf%UFg(1,i,j,k)
                  ! print *, fs%rhoU(i,j,k),fs%rhoU(i+1,j,k),fs%rho_l*vf%UFl(1,i:i+1,j,k)+fs%rho_g*vf%UFg(1,i:i+1,j,k)
                  print *, sum(fs%itp_x(:,i,j,k)*fs%U(i-1:i,j,k)),sum(fs%itp_x(:,i+1,j,k)*fs%U(i:i+1,j,k))
                  ! print *, resU(i-2:i+2,j,k)
                  ! print *, fs%U(i-2:i+2,j,k)
                  ! print *, fs%rhoU(i,j,k)*sum(fs%itp_x(:,i,j,k)*fs%U(i-1:i,j,k)),fs%rhoU(i+1,j,k)*sum(fs%itp_x(:,i+1,j,k)*fs%U(i:i+1,j,k))
                  ! print *, resU(i,j,k),fs%rho(i,j,k),fs%U(i,j,k)
               !    ! fs%indicator(i,j,k)=1
               !    print  *, "This is the U velocity location", i,j,k,fs%U(i,j,k)
               !    print *, "U info", vf%MFX(1,i-1:i+1,j,k), fs%U(i-1:i+1,j,k)
               !    print *, "U info", vf%MFY(1,i-1:i,j,k),vf%MFY(1,i-1:i,j+1,k)
               !    print *, 'Band', vf%band(i-2:i+1,j,k)

               !    print  *, "This is the V velocity location", i,j,k,fs%V(i,j,k)
               !    print *, "V info", vf%MFY(2,i,j-1:j+1,k)
               !    print *, "V info", vf%MFX(2,i,j-1:j,k),vf%MFX(2,i+1,j-1:j,k)
               end if

               ! if ((i.eq.54).and.(j.eq.58).and.(k.eq.1)) then
               !    print *, "U info", vf%MFX(1,i-1:i+1,j,k)
               !    print *, "U info", vf%MFX(1,i-1:i,j:j+1,k)
               ! end if

            end do 
         end do 
      end do
   end subroutine output_info
   
   
end module simulation
!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use iterator_class,    only: iterator
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
   type(tpcons),      public :: fs     !< Two-phase conservative flow solver
   type(vfs),         public :: vf     !< Volume fraction solver
   type(timetracker), public :: time   !< Time info
   
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
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
   real(WP) :: V_drop

   
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
      
      
      ! ! Create an iterator for removing VOF at edges
      ! create_iterator: block
      !    vof_removal_layer=iterator(cfg,'VOF removal',vof_removal_layer_locator)
      !    vof_removed=0.0_WP
      ! end block create_iterator
      
      
      ! Create a two-phase flow solver without bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         real(WP) :: Re,Fr,We,r,m
         call param_read('Method type',type_method)
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Two-phase NS')
         ! Read in physical parameters
         call param_read('Gravity',fs%gravity)
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         call param_read('Surface tension coefficient',fs%sigma)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         if (Implicit) then
            vs=ddadi(cfg=cfg,name='Velocity',nst=7)
            call fs%setup(pressure_solver=ps,implicit_solver=vs)
         else
            call fs%setup(pressure_solver=ps)
         end if
      end block create_flow_solver
      
      
      ! Initialize velocity to zero
      initialize_velocity: block
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         fs%Uf=0.0_WP; fs%Vf=0.0_WP; fs%Wf=0.0_WP
         
         ! Initialize density
         fs%rho=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_faceRHO(vf=vf,rho=fs%rho)
         
         call fs%get_div()
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
         call ens_out%add_vector('velocity',fs%U,fs%V,fs%W)
         call ens_out%add_vector('velface',fs%Uf,fs%Vf,fs%Wf)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('band',vf%band)
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
         call mfile%add_column(V_drop,'Droplet Velocity')
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
      end function levelset_falling_drop
      
      
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
         vf%VFold=vf%VF;fs%rhoold=fs%rho
         vf%bandold=vf%band
         ! Remember old velocity and face densities
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         call vf%advance(dt=time%dt,U=fs%Uf,V=fs%Vf,W=fs%Wf,Uc=fs%U,Vc=fs%V,Wc=fs%W,rho_l=fs%rho_l,rho_g=fs%rho_g)
         ! Update face density and momentum vector
         fs%rho=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_faceRHO(vf=vf,rho=fs%rho)
         fs%rhoU=fs%rho_l*vf%UFl(1,:,:,:)+fs%rho_g*vf%UFg(1,:,:,:)
         fs%rhoV=fs%rho_l*vf%UFl(2,:,:,:)+fs%rho_g*vf%UFg(2,:,:,:)
         fs%rhoW=fs%rho_l*vf%UFl(3,:,:,:)+fs%rho_g*vf%UFg(3,:,:,:)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            ! Build mid-time velocity =========================================
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(vf,resU,resV,resW)

            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)

            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho*fs%U+(fs%rho+fs%rhoold)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho*fs%V+(fs%rho+fs%rhoold)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho*fs%W+(fs%rho+fs%rhoold)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(vf,time%dt,resU,resV,resW)

            ! Compute predictor U
            fs%U=2.0_WP*fs%U-fs%Uold+resU; call cfg%sync(fs%U)
            fs%V=2.0_WP*fs%V-fs%Vold+resV; call cfg%sync(fs%V)
            fs%W=2.0_WP*fs%W-fs%Wold+resW; call cfg%sync(fs%W)

            ! ! Update viscosity explictly
            ! call fs%viscosity_gravity_explict(vf,time%dt)

            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%update_faceU(vf,fs%U,fs%V,fs%W,fs%Uf,fs%Vf,fs%Wf)
            call fs%update_pgrad_all(vf,time%dt)
            ! call fs%update_faceU_correction(vf,time%dt)
            call fs%apply_bcond(time%dt,'face')
            call fs%correct_mfr()
            call fs%get_div()
            if (STflag) call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            ! Corrector step
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%Uf=fs%Uf-time%dt*resU/fs%RHOX
            fs%Vf=fs%Vf-time%dt*resV/fs%RHOY
            fs%Wf=fs%Wf-time%dt*resW/fs%RHOZ

            call fs%get_cell_pgrad(vf,fs%psolv%sol,resU,resV,resW,.true.)
            fs%U=fs%U-time%dt*resU/fs%rho; call cfg%sync(fs%U)
            fs%V=fs%V-time%dt*resV/fs%rho; call cfg%sync(fs%V)
            fs%W=fs%W-time%dt*resW/fs%rho; call cfg%sync(fs%W)
            call fs%apply_bcond(time%dt,'cell')
            ! Increment sub-iteration counter =================================
            time%it=time%it+1
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         
         call get_stats()
         
         call mfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   subroutine get_stats()
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools,  only: Pi
      implicit none
      integer :: i,j,k,ierr
      real(WP) ::  vol_drop
      V_drop=0.0_WP; vol_drop=0.0_WP
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               V_drop=V_drop+vf%VF(i,j,k)*cfg%vol(i,j,k)*fs%V(i,j,k)
               vol_drop=vol_drop+vf%VF(i,j,k)*cfg%vol(i,j,k)
            end do
         end do
      end do
      
      call MPI_ALLREDUCE(MPI_IN_PLACE, V_drop  , 1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, vol_drop, 1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      
      if (vol_drop.gt.0.0_WP) then
         V_drop = V_drop / vol_drop
      else
         V_drop=0.0_WP
      end if
      
   end subroutine get_stats
   
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
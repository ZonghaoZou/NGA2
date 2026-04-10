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
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi

   !> Problem definition
   real(WP), dimension(3) :: center
   real(WP) :: radius
   real(WP) :: U0  ! Translation velocity component

   !> Global Statistics for Monitor
   real(WP), public :: TKE          ! Total Kinetic Energy
   real(WP), public :: vol_liq      ! Total Liquid Volume
   real(WP), public :: xc, yc, zc   ! Instantaneous Center of Mass
   real(WP), public :: lmax=0.0_WP  ! Semi-axis max extent
   real(WP), public :: lymax=0.0_WP ! Semi-axis y extent
   real(WP), public :: vel_perturb_max=0.0_WP  ! Max velocity perturbation |u-u0|
   real(WP), public :: vel_perturb_rms=0.0_WP  ! RMS velocity perturbation
   real(WP), public :: vel_perturb_mean=0.0_WP ! Mean velocity perturbation
   real(WP), public :: Ca_max=0.0_WP  ! Max capillary number
   real(WP), public :: Ca_rms=0.0_WP  ! RMS capillary number
   real(WP), public :: Ca_mean=0.0_WP ! Mean capillary number
   
contains
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
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
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: plicnet,lvira,VFhi,VFlo,remap,flux
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux,name='VOF')
         call param_read('Droplet center',center)
         call param_read('Droplet diameter',radius); radius=radius/2.0_WP
         call param_read('Translation velocity',U0)

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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
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
      
      
      ! Create a two-phase flow solver without bconds (all periodic)
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
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
         
         ! Initialize with uniform velocity field u0 = [U0, U0, U0]
         fs%U=U0; fs%V=U0; fs%W=U0
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

      end block create_and_initialize_flow_solver
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='TranslatingDrop')
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
         call mfile%add_column(vel_perturb_max,'Max vel perturbation')
         call mfile%add_column(vel_perturb_rms,'RMS vel perturbation')
         call mfile%add_column(vel_perturb_mean,'Mean vel perturbation')
         call mfile%add_column(Ca_max,'Ca_max')
         call mfile%add_column(Ca_rms,'Ca_rms')
         call mfile%add_column(Ca_mean,'Ca_mean')
         call mfile%add_column(TKE, 'Kinetic Energy')
         call mfile%add_column(lmax,'Semi-major Axis')
         call mfile%add_column(lymax,'Semi-minor Axis')
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
      !> Level set for a perfect sphere (no perturbation)
      function levelset_sphere(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G = radius - sqrt(sum((xyz-center)**2))
      end function levelset_sphere
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - Default NGA2 staggered scheme
   subroutine simulation_run
      use tpns_class, only: static_contact,harmonic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         call get_stats()

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
      
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
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
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! ! Output to ensight
         ! if (ens_evt%occurs()) then
         !    call vf%update_surfmesh(smesh)
         !    call ens_out%write_data(time%t)
         ! end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()

         ! Write final Ca to CSV at the end
         if (time%done()) call get_Ca()
         
      end do
      
   end subroutine simulation_run
   
   !> Calculate velocity perturbation and capillary number diagnostics
   !> Following the same pattern as the spurious current get_Ca method
   subroutine get_stats
      use irl_fortran_interface
      use mpi_f08,   only: MPI_ALLREDUCE, MPI_MAX, MPI_SUM, MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use vfs_class, only: VFlo, VFhi
      implicit none
      
      real(WP), dimension(3) :: ploc
      real(WP) :: cell_vol, rho_cell, v_sq
      real(WP) :: mom_x, mom_y, mom_z
      real(WP) :: U_perturb_sq, U_perturb_max_sq, U_perturb_sum_sq, U_perturb_sum
      real(WP) :: ncell
      integer :: i,j,k,nplane,ierr

      ! ========================================================================
      ! PASS 1: Integrals (TKE, Volume, CoM, Velocity perturbation)
      ! ========================================================================
      TKE     = 0.0_WP
      vol_liq = 0.0_WP
      mom_x   = 0.0_WP
      mom_y   = 0.0_WP
      mom_z   = 0.0_WP
      U_perturb_max_sq = -huge(1.0_WP)
      U_perturb_sum_sq = 0.0_WP
      U_perturb_sum    = 0.0_WP
      ncell            = 0.0_WP
      
      do k = fs%cfg%kmin_, fs%cfg%kmax_
         do j = fs%cfg%jmin_, fs%cfg%jmax_
            do i = fs%cfg%imin_, fs%cfg%imax_
               
               cell_vol = vf%cfg%vol(i,j,k)
               ncell = ncell + 1.0_WP
               
               ! --- Liquid Volume & Center of Mass Accumulation ---
               if (vf%VF(i,j,k) > epsilon(1.0_WP)) then
                  vol_liq = vol_liq + vf%VF(i,j,k) * cell_vol
                  mom_x = mom_x + (vf%VF(i,j,k) * cell_vol) * fs%cfg%xm(i)
                  mom_y = mom_y + (vf%VF(i,j,k) * cell_vol) * fs%cfg%ym(j)
                  mom_z = mom_z + (vf%VF(i,j,k) * cell_vol) * fs%cfg%zm(k)
               end if

               ! --- Total Kinetic Energy ---
               rho_cell = fs%rho_l * vf%VF(i,j,k) + fs%rho_g * (1.0_WP - vf%VF(i,j,k))
               v_sq = Ui(i,j,k)**2 + Vi(i,j,k)**2 + Wi(i,j,k)**2
               TKE = TKE + 0.5_WP * rho_cell * v_sq * cell_vol
               
               ! --- Velocity perturbation: |u - u0|^2 ---
               ! Following the spurious current pattern (accumulate squared, take sqrt after)
               U_perturb_sq = (Ui(i,j,k) - U0)**2 + (Vi(i,j,k) - U0)**2 + (Wi(i,j,k) - U0)**2
               U_perturb_max_sq = max(U_perturb_max_sq, U_perturb_sq)
               U_perturb_sum_sq = U_perturb_sum_sq + U_perturb_sq
               U_perturb_sum    = U_perturb_sum + sqrt(U_perturb_sq)
               
            end do
         end do
      end do

      ! Reduce Integrals
      call MPI_ALLREDUCE(MPI_IN_PLACE, TKE,     1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, vol_liq, 1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, mom_x,   1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, mom_y,   1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, mom_z,   1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, U_perturb_max_sq, 1, MPI_REAL_WP, MPI_MAX, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, U_perturb_sum_sq, 1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, U_perturb_sum,    1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, ncell,            1, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)

      ! Compute velocity perturbation metrics (following spurious current reference)
      vel_perturb_max  = sqrt(U_perturb_max_sq)
      vel_perturb_rms  = sqrt(U_perturb_sum_sq / ncell)
      vel_perturb_mean = U_perturb_sum / ncell

      ! Capillary numbers: Ca = |u - u0| * mu_l / sigma  (using visc_l per reference)
      Ca_max  = vel_perturb_max  * fs%visc_l / fs%sigma
      Ca_rms  = vel_perturb_rms  * fs%visc_l / fs%sigma
      Ca_mean = vel_perturb_mean * fs%visc_l / fs%sigma

      ! Center of Mass
      if (vol_liq > epsilon(1.0_WP)) then
         xc = mom_x / vol_liq
         yc = mom_y / vol_liq
         zc = mom_z / vol_liq
      else
         xc = center(1)
         yc = center(2)
         zc = center(3)
      end if

      ! ========================================================================
      ! PASS 2: Interface Extents (Relative to Instantaneous CoM)
      ! ========================================================================
      lmax  = -huge(1.0_WP)
      lymax = -huge(1.0_WP)
      
      do k = vf%cfg%kmin_, vf%cfg%kmax_
         do j = vf%cfg%jmin_, vf%cfg%jmax_
            do i = vf%cfg%imin_, vf%cfg%imax_
               
               if(vf%VF(i,j,k).gt.VFlo .and. vf%VF(i,j,k).lt.VFhi) then
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        ploc = calculateCentroid(vf%interface_polygon(nplane,i,j,k))
                        if (abs(ploc(1) - xc) .gt. lmax) then
                           lmax = abs(ploc(1) - xc)
                        end if
                        if (abs(ploc(2) - yc) .gt. lymax) then
                           lymax = abs(ploc(2) - yc)
                        end if
                     end if
                  end do
               end if
            end do
         end do
      end do

      ! Reduce Extents
      call MPI_ALLREDUCE(MPI_IN_PLACE, lmax,  1, MPI_REAL_WP, MPI_MAX, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, lymax, 1, MPI_REAL_WP, MPI_MAX, vf%cfg%comm, ierr)

   end subroutine get_stats

   !> Write final Ca numbers to a CSV file (following spurious current pattern)
   subroutine get_Ca
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use param, only: param_read
      implicit none
      integer :: i,j,k
      real(WP) :: U_max,U_sum,U_mean_sum,Utmp,ncell
      integer :: ierr,meshsize
      character(len=20) :: filename
      U_max=-huge(1.0_WP);ncell=0.0_WP;U_sum=0.0_WP;U_mean_sum=0.0_WP
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               ncell=ncell+1.0_WP
               ! Get sum of velocity perturbation squared for Ca_rms
               Utmp=(Ui(i,j,k)-U0)**2+(Vi(i,j,k)-U0)**2+(Wi(i,j,k)-U0)**2
               U_sum=U_sum+Utmp
               U_max=max(U_max,Utmp)
               U_mean_sum=U_mean_sum+sqrt(Utmp)
            end do
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,U_sum,1,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,ncell,1,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,U_max,1,MPI_REAL_WP,MPI_MAX,vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,U_mean_sum,1,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
      U_sum=U_sum/ncell
      Ca_rms =sqrt(U_sum)*fs%visc_l/fs%sigma
      Ca_max =sqrt(U_max)*fs%visc_l/fs%sigma
      Ca_mean=U_mean_sum/ncell*fs%visc_l/fs%sigma

      if (cfg%amRoot) then
         call param_read("nx",meshsize); write(filename, '(I0, ".csv")') meshsize
         open(unit=10, file=filename, status="unknown", position="append", action="write")
         write(10, '(F24.16,A,F24.16,A,F24.16)') Ca_rms, ',', Ca_max, ',', Ca_mean
         close(10)
      end if
   end subroutine get_Ca

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation

!> Various definitions and tools for running an NGA2 simulation
!> Dam break — Scheme 3: SL Momentum (Collocated Grid)
!>
!> Water column (VF=1) in bottom-left corner: [0,2a] x [0,2a]
!> Air (VF=0) elsewhere. Gravity g=9.81 downward
!> Benchmark: Martin & Moyce (1952)
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
   public :: simulation_init,simulation_run,simulation_final,get_KE
   
   !> Flow solver objects
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs
   type(tpcons),      public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Monitoring files
   type(monitor) :: mfile,cflfile
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   
   !> Problem definition: Dam break
   real(WP) :: dam_width   ! Width of water column from wall = semi-base a
   real(WP) :: dam_height  ! Height of water column = n^2 * a
   real(WP) :: a_half      ! Semi-base a (= dam_width, for non-dimensionalization)
   
   !> Stats for monitoring
   real(WP) :: front_x    ! X-position of surge front
   real(WP) :: column_h   ! Height of remaining water column at left wall
   real(WP) :: Z_nd       ! Non-dimensional front position z/a (starts at 1)
   real(WP) :: T_nd       ! Non-dimensional time n*t*sqrt(g/a)
   real(WP) :: KE,VFInt,rhoInt,rhoUInt,rhoVInt,rhoWInt
   real(WP) :: KE_l,KE_g,rhoU_l,rhoV_l,rhoW_l,rhoU_g,rhoV_g,rhoW_g
   real(WP) :: EN_l,EN_g,PE_l,PE_g
   real(WP) :: KE_x,KE_y,KE_xl,KE_xg,KE_yl,KE_yg
   character(len=20) :: filename='output.csv'
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
         use vfs_class, only: lvira,VFhi,VFlo,flux_storage,flux,plicnet
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux,name='VOF')
         ! Read dam break parameters
         call param_read('Dam width',dam_width)
         call param_read('Dam height',dam_height)
         a_half = dam_width  ! semi-base a = column width from wall
         ! Initialize the interface via level set
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
                  ! Call adaptive refinement code to get volume and barycenters
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_dambreak,0.0_WP,amr_ref_lvl)
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
         ! Reset moments
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver with wall BCs
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         use tpcons_class,    only: dirichlet
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Two-phase NS')
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
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP;fs%V=0.0_WP;fs%W=0.0_WP
         fs%Uf=0.0_WP;fs%Vf=0.0_WP;fs%Wf=0.0_WP
         ! Calculate divergence
         call fs%get_div()
         fs%rho=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_faceRHO(vf=vf,rho=fs%rho)
         if (cfg%amRoot) then
            open(unit=10, file=filename, status="replace", action="write")
            write(10, '(A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24)') &
            & 'time', 'timestep', 'KE', 'KE_l', 'KE_g', 'rhoU', 'rhoV', 'rhoW', 'rhoU_l', 'rhoV_l', 'rhoW_l','rhoU_g', 'rhoV_g', 'rhoW_g','PE_l','PE_g','EN_l','EN_g'
            close(unit=10)
         end if
         call get_KE()
      end block create_flow_solver
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='DamBreak')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%U,fs%V,fs%W)
         call ens_out%add_vector('velface',fs%Uf,fs%Vf,fs%Wf)
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
         call get_front()
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
         call mfile%add_column(front_x,'Front X')
         call mfile%add_column(column_h,'Column H')
         call mfile%add_column(Z_nd,'Z nondim')
         call mfile%add_column(T_nd,'T nondim')
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
      
      !> Level set for dam break
      !> Rectangular water column: [0, dam_width] x [0, dam_height]
      function levelset_dambreak(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G = min(dam_width - xyz(1), dam_height - xyz(2))
      end function levelset_dambreak
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation — SL Momentum collocated scheme
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
          
            ! Build mid-time velocity
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

            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%update_faceU(vf,fs%U,fs%V,fs%W,fs%Uf,fs%Vf,fs%Wf)
            call fs%update_pgrad_all(vf,time%dt)
            call fs%apply_bcond(time%dt,'face')
            call fs%correct_mfr()

            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
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
         
         ! Recompute divergence
         call fs%get_div()

         call get_front()
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
   
   
   !> Track surge front position and water column height
   !> Uses IRL interface polygon centroids for accurate tracking
   subroutine get_front()
      use irl_fortran_interface
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      use vfs_class, only: VFlo,VFhi
      implicit none
      real(WP), dimension(3) :: ploc
      integer :: i,j,k,nplane,ierr
      
      front_x  = 0.0_WP
      column_h = 0.0_WP
      
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               if (vf%VF(i,j,k).gt.VFlo .and. vf%VF(i,j,k).lt.VFhi) then
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        ploc = calculateCentroid(vf%interface_polygon(nplane,i,j,k))
                        if (ploc(1) .gt. front_x) front_x = ploc(1)
                        if (ploc(2) .gt. column_h) column_h = ploc(2)
                     end if
                  end do
               end if
            end do
         end do
      end do
      
      call MPI_ALLREDUCE(MPI_IN_PLACE, front_x,  1, MPI_REAL_WP, MPI_MAX, vf%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, column_h, 1, MPI_REAL_WP, MPI_MAX, vf%cfg%comm, ierr)
      
      ! Z = z/a (starts at 1),  T = n*t*sqrt(g/a) where n^2 = H/a
      Z_nd = front_x / a_half
      T_nd = time%t * sqrt((dam_height / a_half) * abs(fs%gravity(2)) / a_half)
      
   end subroutine get_front
   
   
   subroutine get_KE
      implicit none
      integer :: i,j,k,ierr
      real(WP), dimension(:,:,:), allocatable:: Utmp
      real(WP), dimension(:,:,:), allocatable:: Vtmp
      real(WP), dimension(:,:,:), allocatable:: Wtmp
      real(WP), dimension(:,:,:), allocatable:: FX
      real(WP), dimension(:,:,:), allocatable:: FY
      real(WP), dimension(:,:,:), allocatable:: FZ
      real(WP), dimension(:,:,:), allocatable:: Stmp
      real(WP), dimension(:,:,:), allocatable:: Stmp2
      allocate(Utmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Utmp=0.0_WP
      allocate(Vtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Vtmp=0.0_WP
      allocate(Wtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Wtmp=0.0_WP
      allocate(Stmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Stmp=0.0_WP
      allocate(Stmp2(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Stmp2=0.0_WP
      allocate(FX(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FX=0.0_WP
      allocate(FY(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FY=0.0_WP
      allocate(FZ(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FZ=0.0_WP
      calculate_KE: block
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  Stmp(i,j,k)=fs%rho(i,j,k)*(fs%U(i,j,k)**2+fs%V(i,j,k)**2+fs%W(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call cfg%integrate(Stmp,integral=KE)
      end block calculate_KE

      calculate_KEl: block
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  Stmp(i,j,k)=fs%rho(i,j,k)*(fs%U(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call cfg%integrate(Stmp,integral=KE_x)
      end block calculate_KEl

      calculate_KEg: block
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  Stmp(i,j,k)=fs%rho(i,j,k)*(fs%V(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call cfg%integrate(Stmp,integral=KE_y)
      end block calculate_KEg
      
      calculate_momentumconservation: block
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=fs%rho(i,j,k)*fs%U(i,j,k)
                  Vtmp(i,j,k)=fs%rho(i,j,k)*fs%V(i,j,k)
                  Wtmp(i,j,k)=fs%rho(i,j,k)*fs%W(i,j,k)
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoUInt)
         call fs%cfg%integrate(Vtmp,integral=rhoVInt)
         call fs%cfg%integrate(Wtmp,integral=rhoWInt)
      end block calculate_momentumconservation 

      ! Get KE_l
      calculate_l_KE: block
         Stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=fs%rho_l*vf%VF(i,j,k)*(fs%U(i,j,k)**2+fs%V(i,j,k)**2+fs%W(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_l)
      end block calculate_l_KE


      calculate_l_KEx: block
         Stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=fs%rho_l*vf%VF(i,j,k)*(fs%U(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_xl)
      end block calculate_l_KEx


      calculate_l_KEy: block
         Stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=fs%rho_l*vf%VF(i,j,k)*(fs%V(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_yl)
      end block calculate_l_KEy


      calculate_l_momentumconservation: block
         Utmp=0.0_WP;Vtmp=0.0_WP;Wtmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=fs%rho_l*vf%VF(i,j,k)*fs%U(i,j,k)
                  Vtmp(i,j,k)=fs%rho_l*vf%VF(i,j,k)*fs%V(i,j,k)
                  Wtmp(i,j,k)=fs%rho_l*vf%VF(i,j,k)*fs%W(i,j,k)
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoU_l)
         call fs%cfg%integrate(Vtmp,integral=rhoV_l)
         call fs%cfg%integrate(Wtmp,integral=rhoW_l)
      end block calculate_l_momentumconservation 

      calculate_g_KE: block
         Stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=fs%rho_g*(1.0_WP-vf%VF(i,j,k))*(fs%U(i,j,k)**2+fs%V(i,j,k)**2+fs%W(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_g)
      end block calculate_g_KE

      calculate_g_KEx: block
         Stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=fs%rho_g*(1.0_WP-vf%VF(i,j,k))*(fs%U(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_xg)
      end block calculate_g_KEx

      calculate_g_KEy: block
         Stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=fs%rho_g*(1.0_WP-vf%VF(i,j,k))*(fs%V(i,j,k)**2)/2.0_WP
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_yg)
      end block calculate_g_KEy



      calculate_g_momentumconservation: block
      Utmp=0.0_WP;Vtmp=0.0_WP;Wtmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=fs%rho_g*(1.0_WP-vf%VF(i,j,k))*fs%U(i,j,k)
                  Vtmp(i,j,k)=fs%rho_g*(1.0_WP-vf%VF(i,j,k))*fs%V(i,j,k)
                  Wtmp(i,j,k)=fs%rho_g*(1.0_WP-vf%VF(i,j,k))*fs%W(i,j,k)
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoU_g)
         call fs%cfg%integrate(Vtmp,integral=rhoV_g)
         call fs%cfg%integrate(Wtmp,integral=rhoW_g)
      end block calculate_g_momentumconservation 

      call fs%cfg%integrate(vf%VF,integral=VFInt)
      call fs%cfg%integrate(fs%rho,integral=rhoInt)


      calculate_PE: block
      Utmp=0.0_WP;Vtmp=0.0_WP
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               Utmp(i,j,k)=vf%VF(i,j,k)*fs%rho_l*cfg%ym(j)*abs(fs%gravity(2))
               Vtmp(i,j,k)=(1.0_WP-vf%VF(i,j,k))*fs%rho_g*cfg%ym(j)*abs(fs%gravity(2))
            end do
         end do
      end do
      call fs%cfg%integrate(Utmp,integral=PE_l)
      call fs%cfg%integrate(Vtmp,integral=PE_g)
   end block calculate_PE

   calculate_enstrophy: block
         real(WP) :: wxc,wyc,wzc
         Utmp=0.0_WP;Vtmp=0.0_WP;Wtmp=0.0_WP
         do k=cfg%kmin_,cfg%kmax_+1
            do j=cfg%jmin_,cfg%jmax_+1
               do i=cfg%imin_,cfg%imax_+1
                  Utmp(i,j,k)=sum(fs%grdw_y(:,i,j,k)*fs%Wf(i,j-1:j,k))-sum(fs%grdv_z(:,i,j,k)*fs%Vf(i,j,k-1:k))
                  Vtmp(i,j,k)=sum(fs%grdu_z(:,i,j,k)*fs%Uf(i,j,k-1:k))-sum(fs%grdw_x(:,i,j,k)*fs%Wf(i-1:i,j,k))
                  Wtmp(i,j,k)=sum(fs%grdv_x(:,i,j,k)*fs%Vf(i-1:i,j,k))-sum(fs%grdu_y(:,i,j,k)*fs%Uf(i,j-1:j,k))
               end do
            end do
         end do
         call cfg%sync(Utmp)
         call cfg%sync(Vtmp)
         call cfg%sync(Wtmp)
         Stmp=0.0_WP;Stmp2=0.0_WP
         do k = cfg%kmin_, cfg%kmax_
            do j = cfg%jmin_, cfg%jmax_
               do i = cfg%imin_, cfg%imax_
                  wxc = 0.25_wp * ( Utmp(i,   j  ,   k  )   &
                                 +   Utmp(i,   j+1,   k  ) &
                                 +   Utmp(i,   j  ,   k+1) &
                                 +   Utmp(i,   j+1,   k+1) )
                  wyc = 0.25_wp * ( Vtmp(i,   j  ,   k  )   &
                                 +   Vtmp(i+1, j  ,   k  ) &
                                 +   Vtmp(i,   j  ,   k+1) &
                                 +   Vtmp(i+1, j  ,   k+1) )
                  wzc = 0.25_wp * ( Wtmp(i,   j  ,   k  )   &
                                 +   Wtmp(i+1, j  ,   k  ) &
                                 +   Wtmp(i,   j+1, k  )   &
                                 +   Wtmp(i+1, j+1, k  ) )
                  Stmp(i,j,k) = 0.5_wp * (wxc*wxc + wyc*wyc + wzc*wzc)*vf%VF(i,j,k)
                  Stmp2(i,j,k) = 0.5_wp * (wxc*wxc + wyc*wyc + wzc*wzc)*(1.0_WP-vf%VF(i,j,k))
               end do
            end do
         end do
         call cfg%integrate(Stmp,integral=EN_l)
         call cfg%integrate(Stmp2,integral=EN_g)
      end block calculate_enstrophy


      if (cfg%amRoot) then
         open(unit=10, file=filename, status="unknown", position="append", action="write")
         write(10, '(F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16)')&
         time%t,time%dt,KE,KE_l,KE_g,rhoUInt,rhoVInt,rhoWInt,rhoU_l,rhoV_l,rhoW_l,rhoU_g,rhoV_g,rhoW_g,PE_l,PE_g,EN_l,EN_g,KE_x,KE_xl,KE_xg,KE_y,KE_yl,KE_yg
         close(10)
      end if

   end subroutine get_KE
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      deallocate(resU,resV,resW)
   end subroutine simulation_final
   
   
end module simulation

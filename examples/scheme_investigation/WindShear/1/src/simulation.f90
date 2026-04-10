!> Various definitions and tools for running an NGA2 simulation
!> Rayleigh-Taylor instability — Scheme 1: Default NGA2 (Staggered Grid)
!>
!> Heavy fluid (rho_l=1.225) on top, light fluid (rho_g=0.1694) below
!> Interface at y=0 with cosine perturbation A0*cos(2*pi*x)
!> No surface tension, gravity g=9.81 downward
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
   
   !> Solvers
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
   
   public :: simulation_init,simulation_run,simulation_final,get_KE
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   
   real(WP), dimension(:), allocatable :: all_time, all_amp,all_amp_height
   real(WP) :: amp, grate, amp0, amp_height
   integer :: nx
   !> Problem definition: KH instability
   real(WP) :: Ll,Lg,Ug,Ul,dg,dl,alpha,Lx
   ! real(WP) :: dg,dl,Lx,Ug,Ll,Lg,Ul
   real(WP) :: eps_kh  ! perturbation amplitude = factor * Ug
   
   
   real(WP) :: KE,VFInt,rhoInt,rhoUInt,rhoVInt,rhoWInt
   real(WP) :: KE_l,KE_g,rhoU_l,rhoV_l,rhoW_l,rhoU_g,rhoV_g,rhoW_g
   real(WP) :: EN_l,EN_g,PE_l,PE_g
   real(WP) :: KE_x,KE_y,KE_xl,KE_xg,KE_yl,KE_yg
   character(len=20) :: filename='output.csv'
   
   !> Traveling wave parameters
   logical :: run_traveling_wave
   real(WP) :: tw_H, tw_omega, tw_z0a, tw_ustar, tw_K


contains   
   subroutine postproc_data()
     use irl_fortran_interface
     use mathtools, only: Pi
     use string,    only: str_medium
     use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX
     use parallel,  only: MPI_REAL_WP
     implicit none
     
     integer :: ierr, i, j, k, my_size
     real(WP), dimension(:), allocatable :: temp
     real(WP), dimension(1:nx) :: local_map, global_map
     
     ! Variables for least-squares fitting
     real(WP) :: k_wave, x_coord, yi, Si, Ci, dx_loc
     real(WP) :: M11, M12, M22, R1, R2, det, C1, C2
     real(WP) :: b_amp, mean_height
     real(WP) :: maxheight,minheight

     ! Calculate new amplitude (store old for growth rate)
     grate=amp

     ! 1. Compute liquid column height at each x-column: integral of VF*dy
     local_map = 0.0_WP
     do k = vf%cfg%kmin_, vf%cfg%kmax_
        do i = vf%cfg%imin_, vf%cfg%imax_
           do j = vf%cfg%jmin_, vf%cfg%jmax_
               if (vf%cfg%z(k).le.0.0_WP.and.vf%cfg%z(k+1).gt.0.0_WP) then
                  local_map(i) = local_map(i) + vf%VF(i,j,k)*vf%cfg%dy(j)
               end if
           end do
        end do
     end do
     call MPI_ALLREDUCE(local_map, global_map, nx, MPI_REAL_WP, MPI_SUM, vf%cfg%comm, ierr)

     ! 2. Compute mean height (should be ~Ly/2 for unperturbed interface at y=0)
     mean_height = sum(global_map) / real(nx, WP)

     maxheight=maxval(global_map)
     minheight=minval(global_map)
     ! 3. Least-squares fit: η(x) = C1*sin(αx) + C2*cos(αx)
     !    Compute x-coordinates analytically (valid on ALL processors)
     k_wave = alpha  ! = 2π/Lx
     dx_loc = Lx / real(nx, WP)

     M11 = 0.0_WP; M12 = 0.0_WP; M22 = 0.0_WP
     R1  = 0.0_WP; R2  = 0.0_WP

     do i = 1, nx
        ! Analytical cell-center coordinate for uniform grid centered at origin
        x_coord = (real(i, WP) - 0.5_WP) * dx_loc - 0.5_WP * Lx

        ! Perturbation from mean unperturbed height
        yi = global_map(i) - mean_height

        Si = sin(k_wave * x_coord)
        Ci = cos(k_wave * x_coord)

        M11 = M11 + Si * Si
        M12 = M12 + Si * Ci
        M22 = M22 + Ci * Ci
        R1  = R1  + yi * Si
        R2  = R2  + yi * Ci
     end do

     ! Solve 2x2 normal equation: [M11 M12; M12 M22] * [C1; C2] = [R1; R2]
     det = M11 * M22 - M12 * M12

     if (abs(det) > 1.0e-30_WP) then
        C1 = (R1 * M22 - R2 * M12) / det
        C2 = (R2 * M11 - R1 * M12) / det
        b_amp = sqrt(C1**2 + C2**2)
     else
        b_amp = 0.0_WP
     end if

     ! Set the new amplitude
   !   if ((maxheight-minheight)/2.0_WP .gt. b_amp) then
   !       amp=(maxheight-minheight)/2.0_WP
   !   else
   !       amp = b_amp
   !   end if
     amp = b_amp
     amp_height=(maxheight-minheight)/2.0_WP

     ! Estimate growth rate
     if (time%t .gt. 0.0_WP) then
        grate = (amp - grate) / time%dt
     else
        grate = 0.0_WP
     end if
     ! Store time and amplitude series
     if (.not.allocated(all_time)) then
        my_size=0
     else
        my_size=size(all_time,dim=1)
     end if
     allocate(temp(my_size+1)); temp(1:my_size)=all_time; temp(my_size+1)=time%t; call MOVE_ALLOC(temp,all_time)
     allocate(temp(my_size+1)); temp(1:my_size)=all_amp ; temp(my_size+1)=amp   ; call MOVE_ALLOC(temp,all_amp )
     allocate(temp(my_size+1)); temp(1:my_size)=all_amp_height ; temp(my_size+1)=amp_height   ; call MOVE_ALLOC(temp,all_amp_height )
   end subroutine postproc_data
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      theoretical_KH: block
         use orr_sommerfeld_kh, only: setup_kh
         use mathtools, only: twoPi
         integer :: N_in,i
         real(WP) :: rho_g_in,rho_l_in,mu_g_in,mu_l_in,sigma_in
         real(WP) :: dwn, alpha_nd
         complex(WP) :: c_eval
         logical :: run_sweep
         
         real(WP), dimension(:), allocatable:: alpha_delta_g, nd_growthrate
         
         call param_read('Ug', Ug)
         call param_read('Delta g', dg)
         call param_read('Delta l', dl)
         call param_read('Gas density', rho_g_in)
         call param_read('Density ratio', rho_l_in); rho_l_in=rho_g_in/rho_l_in
         call param_read('Gas dynamic viscosity', mu_g_in)
         call param_read('Viscosity ratio', mu_l_in); mu_l_in=mu_g_in/mu_l_in
         call param_read('Surface tension coefficient', sigma_in)
         call param_read('Lg', Lg)
         call param_read('Ll', Ll)
         
         Ll=Ll*dl; Lg=Lg*dg
         Ul = Ug * mu_g_in * dl / (dg * mu_l_in)
         N_in=128
         
         ! Get simulation wavenumber from input
         call param_read('Wavenumber', alpha_nd)  ! non-dimensional α·δ_g
         alpha = alpha_nd / dg                     ! dimensional α [1/m]
         Lx = twoPi / alpha                        ! wavelength [m]
         nx = cfg%nx                                ! actual grid size from geometry
         
         ! Check if user wants theoretical sweep
         call param_read('Run theoretical KH', run_sweep, default=.false.)
         
         if (run_sweep .and. cfg%amRoot) then
            allocate(alpha_delta_g(1:200)); alpha_delta_g=0.0_WP
            allocate(nd_growthrate(1:200)); nd_growthrate=0.0_WP
            
            open(unit=2, file='sol.csv', status="replace", action="write")
            write(2, '(A24,A24)') 'alpha_delta_g', 'nd_growthrate'
            
            dwn=4.0_WP/200
            do i = 1,200
               call setup_kh(N_in,dwn*i/dg,Ug,Ul,dg,dl,rho_g_in,&
               rho_l_in,mu_g_in,mu_l_in,sigma_in,Lg,Ll,c_eval)
               
               alpha_delta_g(i) = dwn*i
               nd_growthrate(i) = (dwn*i/dg) * aimag(c_eval) * dg / Ug
               
               write(2, '(ES15.6, 1X, ES15.6)') alpha_delta_g(i), nd_growthrate(i)
            end do
            
            close(2)
            print *, "Done theoretical KH calculation"
            deallocate(alpha_delta_g, nd_growthrate)
         end if

         ! Solve OS for the SIMULATION wavenumber
         call param_read('Perturbation factor', eps_kh, default=1.0e-3_WP)
         eps_kh = eps_kh * Ug
         
         call setup_kh(N_in, alpha, Ug, Ul, dg, dl, rho_g_in, &
              rho_l_in, mu_g_in, mu_l_in, sigma_in, Lg, Ll, c_eval)
         
         if (cfg%amRoot) then
            print '(A,ES12.4)', ' OS eigenvalue c_r = ', real(c_eval)
            print '(A,ES12.4)', ' OS eigenvalue c_i = ', aimag(c_eval)
            print '(A,ES12.4)', ' ND growth rate    = ', alpha * aimag(c_eval) * dg / Ug
            print '(A,ES12.4)', ' Perturbation eps  = ', eps_kh
         end if
         
      end block theoretical_KH

      read_traveling_wave: block
         call param_read('Run traveling wave', run_traveling_wave, default=.false.)
         if (run_traveling_wave) then
            call param_read('Traveling wave H', tw_H, default=0.01_WP)
            call param_read('Traveling wave omega', tw_omega, default=20.39_WP)
            call param_read('Traveling wave z0a', tw_z0a, default=0.003_WP)
            call param_read('Traveling wave ustar', tw_ustar, default=0.3_WP)
            call param_read('Traveling wave K', tw_K, default=0.4_WP)
            if (cfg%amRoot) then
               print '(A)', ' --- Traveling Wave Initialization Enabled --- '
               print '(A,ES12.4)', ' H     = ', tw_H
               print '(A,ES12.4)', ' omega = ', tw_omega
               print '(A,ES12.4)', ' z0a   = ', tw_z0a
               print '(A,ES12.4)', ' ustar = ', tw_ustar
               print '(A,ES12.4)', ' K     = ', tw_K
            end if
         end if
      end block read_traveling_wave



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
         use vfs_class, only: lvira,VFhi,VFlo,remap,r2p,flux
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux,name='VOF')
         ! amp0 for normalization (used in ODRPACK fit)
         ! Initialize the interface via level set (using OS eigenvector)
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_KH,0.0_WP,amr_ref_lvl)
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
      
      
      ! Create a two-phase flow solver with wall BCs
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: dirichlet,slip
         ! (Uint, di removed — no longer needed with OS init)
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Gas dynamic viscosity',fs%visc_g)
         call param_read('Viscosity ratio',fs%visc_l); fs%visc_l=fs%visc_g/fs%visc_l
         ! Assign constant density to each phase
         call param_read('Gas density',fs%rho_g)
         call param_read('Density ratio',fs%rho_l);  fs%rho_l=fs%rho_g/fs%rho_l
         ! Read in surface tension coefficient 
         call param_read('Surface tension coefficient',fs%sigma)
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Add slip conditions top and bottom
         call fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         call fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         
         ! Initialize velocity using OS eigenvectors
         init_velocity: block
            use orr_sommerfeld_kh, only: eval_kh
            real(WP) :: u_tmp, v_tmp, G_tmp
            integer :: i,j,k
            if (run_traveling_wave) then
               call init_traveling_wave_velocity()
            else
               ! U vel (face x): eval at (x(i), ym(j))
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_+1
                        call eval_kh(fs%cfg%x(i), fs%cfg%ym(j), eps_kh, u_tmp, v_tmp, G_tmp)
                        fs%U(i,j,k) = u_tmp
                     end do
                  end do
               end do
               ! V vel (face y): eval at (xm(i), y(j))
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_+1
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        call eval_kh(fs%cfg%xm(i), fs%cfg%y(j), eps_kh, u_tmp, v_tmp, G_tmp)
                        fs%V(i,j,k) = v_tmp
                     end do
                  end do
               end do
               ! W = 0 (2D problem)
            end if
         end block init_velocity
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         ! Make it solenoidal
         call fs%get_olddensity(vf)
         fs%rho_U=fs%rho_Uold
         fs%rho_V=fs%rho_Vold
         fs%rho_W=fs%rho_Wold
         call fs%update_laplacian()
         call fs%get_div()
         fs%psolv%rhs=-fs%cfg%vol*fs%div
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         fs%U=fs%U-resU/fs%rho_U
         fs%V=fs%V-resV/fs%rho_V
         fs%W=fs%W-resW/fs%rho_W
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         ! Print out mixing layer definition
         ! if (fs%cfg%amRoot) then
         !    write(message,'("[Initial conditions] => Gas thickness =",es12.5)')   dg; call log(message)
         !    write(message,'("[Initial conditions] => Gas  velocity =",es12.5)')   Ug; call log(message)
         !    write(message,'("[Initial conditions] => Liq thickness =",es12.5)')   dl; call log(message)
         !    write(message,'("[Initial conditions] => Liq  velocity =",es12.5)')   Ul; call log(message)
         !    write(message,'("[Initial conditions] => Deficit param =",es12.5)')   di; call log(message)
         !    write(message,'("[Initial conditions] => Interface vel =",es12.5)') Uint; call log(message)
         ! end if

         ! if (cfg%amRoot) then
         !    open(unit=10, file=filename, status="replace", action="write")
         !    write(10, '(A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24)') &
         !    & 'time', 'timestep', 'KE', 'KE_l', 'KE_g', 'rhoU', 'rhoV', 'rhoW', 'rhoU_l', 'rhoV_l', 'rhoW_l','rhoU_g', 'rhoV_g', 'rhoW_g','VF','rho'
         !    close(unit=10)
         ! end if
         ! call get_KE()
      end block create_and_initialize_flow_solver
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='MPKH')
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
         call mfile%add_column(amp,'amplitude')
         call mfile%add_column(amp_height,'amplitude_height')
         
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
      
      function levelset_KH(xyz,t) result(G)
         use orr_sommerfeld_kh, only: eval_kh
         use mathtools, only : twoPi
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         real(WP) :: u_tmp, v_tmp
         if (run_traveling_wave) then
            G = (tw_H / 2.0_WP) * cos(alpha * xyz(1)) - xyz(2)
         else
            call eval_kh(xyz(1), xyz(2), eps_kh, u_tmp, v_tmp, G)
         end if
      end function levelset_KH
      
      !> Function that localizes the top (y+) of the domain
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      
      !> Function that localizes the bottom (y-) of the domain
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      subroutine init_traveling_wave_velocity()
         implicit none
         integer :: i,j,k
         real(WP) :: GammaH, abs_gradH
         real(WP), dimension(:,:,:), allocatable :: H_val, gradH_x, gradH_y
         real(WP), dimension(:,:,:), allocatable :: psi
         real(WP) :: y_dist, u_tmp

         allocate(H_val(fs%cfg%imino_:fs%cfg%imaxo_, fs%cfg%jmino_:fs%cfg%jmaxo_, fs%cfg%kmino_:fs%cfg%kmaxo_)); H_val = 0.0_WP
         allocate(gradH_x(fs%cfg%imino_:fs%cfg%imaxo_, fs%cfg%jmino_:fs%cfg%jmaxo_, fs%cfg%kmino_:fs%cfg%kmaxo_)); gradH_x = 0.0_WP
         allocate(gradH_y(fs%cfg%imino_:fs%cfg%imaxo_, fs%cfg%jmino_:fs%cfg%jmaxo_, fs%cfg%kmino_:fs%cfg%kmaxo_)); gradH_y = 0.0_WP
         allocate(psi(fs%cfg%imino_:fs%cfg%imaxo_, fs%cfg%jmino_:fs%cfg%jmaxo_, fs%cfg%kmino_:fs%cfg%kmaxo_)); psi = 0.0_WP

         ! Construct H based on VF.
         do k = fs%cfg%kmin_, fs%cfg%kmax_
            do j = fs%cfg%jmin_, fs%cfg%jmax_
               do i = fs%cfg%imin_, fs%cfg%imax_
                  if (vf%VF(i,j,k) > 0.0_WP) then
                     H_val(i,j,k) = 1.0_WP
                  else
                     H_val(i,j,k) = 0.0_WP
                  end if
               end do
            end do
         end do
         call fs%cfg%sync(H_val)

         ! Calculate |grad H| at cell centers for RHS
         do k = fs%cfg%kmin_, fs%cfg%kmax_
            do j = fs%cfg%jmin_, fs%cfg%jmax_
               do i = fs%cfg%imin_, fs%cfg%imax_
                  gradH_x(i,j,k) = (max(H_val(i+1,j,k), H_val(i,j,k)) - max(H_val(i,j,k), H_val(i-1,j,k))) * fs%cfg%dxi(i)
                  gradH_y(i,j,k) = (max(H_val(i,j+1,k), H_val(i,j,k)) - max(H_val(i,j,k), H_val(i,j-1,k))) * fs%cfg%dyi(j)
               end do
            end do
         end do

         ! Temporarily set density to 1
         fs%rho_U = 1.0_WP
         fs%rho_V = 1.0_WP
         fs%rho_W = 1.0_WP
         call fs%update_laplacian()

         if (.not. fs%cfg%yper) then
            if (fs%cfg%jproc == 1) then
               do k = fs%cfg%kmin_, fs%cfg%kmax_
                  do i = fs%cfg%imin_, fs%cfg%imax_
                     fs%psolv%opr(1,i,fs%cfg%jmin_,k) = fs%psolv%opr(1,i,fs%cfg%jmin_,k) + 2.0_WP * fs%cfg%vol(i,fs%cfg%jmin_,k) / (fs%cfg%dy(fs%cfg%jmin_) * fs%cfg%dy(fs%cfg%jmin_))
                  end do
               end do
            end if
            if (fs%cfg%jproc == fs%cfg%npy) then
               do k = fs%cfg%kmin_, fs%cfg%kmax_
                  do i = fs%cfg%imin_, fs%cfg%imax_
                     fs%psolv%opr(1,i,fs%cfg%jmax_,k) = fs%psolv%opr(1,i,fs%cfg%jmax_,k) + 2.0_WP * fs%cfg%vol(i,fs%cfg%jmax_,k) / (fs%cfg%dy(fs%cfg%jmax_) * fs%cfg%dy(fs%cfg%jmax_))
                  end do
               end do
            end if
         end if

         ! Re-setup psolv with updated operator matrices
         call fs%psolv%setup()

         ! Set RHS
         do k = fs%cfg%kmin_, fs%cfg%kmax_
            do j = fs%cfg%jmin_, fs%cfg%jmax_
               do i = fs%cfg%imin_, fs%cfg%imax_
                  abs_gradH = sqrt(gradH_x(i,j,k)**2 + gradH_y(i,j,k)**2)
                  GammaH = tw_H * tw_omega * cos(alpha * fs%cfg%xm(i))
                  fs%psolv%rhs(i,j,k) = fs%cfg%vol(i,j,k) * GammaH * abs_gradH
                  fs%psolv%sol(i,j,k) = 0.0_WP
               end do
            end do
         end do

         ! Solve
         call fs%psolv%solve()
         psi = fs%psolv%sol
         call fs%cfg%sync(psi)

         ! Initialize Liquid velocity using psi, Gas velocity using log profile
         do k = fs%cfg%kmin_, fs%cfg%kmax_
            do j = fs%cfg%jmin_, fs%cfg%jmax_
               do i = fs%cfg%imin_, fs%cfg%imax_+1
                  ! U liquid
                  fs%U(i,j,k) = ( psi(i,j+1,k) + psi(i-1,j+1,k) - psi(i,j-1,k) - psi(i-1,j-1,k) ) / (2.0_WP * (fs%cfg%ym(j+1)-fs%cfg%ym(j-1)))
                  ! Gas velocity
                  y_dist = fs%cfg%ym(j) - (tw_H / 2.0_WP) * cos(alpha * fs%cfg%x(i))
                  if (y_dist > tw_z0a) then
                     u_tmp = (tw_ustar / tw_K) * log(y_dist / tw_z0a)
                  else
                     u_tmp = 0.0_WP
                  end if
                  ! Blend
                  if (max(H_val(i,j,k), H_val(i-1,j,k)) > 0.5_WP) then
                     fs%U(i,j,k) = fs%U(i,j,k) ! Liquid
                  else
                     fs%U(i,j,k) = u_tmp ! Gas
                  end if
               end do
            end do
         end do

         do k = fs%cfg%kmin_, fs%cfg%kmax_
            do j = fs%cfg%jmin_, fs%cfg%jmax_+1
               do i = fs%cfg%imin_, fs%cfg%imax_
                  ! V liquid = -d(psi)/dx
                  fs%V(i,j,k) = -( psi(i+1,j,k) + psi(i+1,j-1,k) - psi(i-1,j,k) - psi(i-1,j-1,k) ) / (2.0_WP * (fs%cfg%xm(i+1)-fs%cfg%xm(i-1)))
                  ! Gas velocity
                  u_tmp = 0.0_WP
                  ! Blend
                  if (max(H_val(i,j,k), H_val(i,j-1,k)) > 0.5_WP) then
                     fs%V(i,j,k) = fs%V(i,j,k) ! Liquid
                  else
                     fs%V(i,j,k) = u_tmp ! Gas
                  end if
               end do
            end do
         end do

         fs%W = 0.0_WP

         deallocate(H_val, gradH_x, gradH_y, psi)

      end subroutine init_traveling_wave_velocity

   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation — Default staggered grid scheme
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
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
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
         
         ! ! Compute stats after
         ! call get_KE()
         
         call postproc_data()
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
      linear_staiblity: block
         use param, only: param_read
         integer :: count, start_idx,N,i
         real(WP) :: sum_x, sum_y, sum_xy, sum_xx
         real(WP) :: x, y, slope,lc,tau, threshold, tskip, skipratio
         integer :: file_unit = 42
         logical :: file_exists
         ! 1. Determine safe start point
         ! Since we used the exact eigenfunction, the transient is tiny.
         ! We simply skip the first 10 data points to let the Poisson solver settle.
         call param_read('Fitting threshold',threshold)
         call param_read('Skip Ratio', skipratio)
         ! tskip=dg/Ug+Lx/Ug*skipratio
         tskip=dg/Ug+Lx/Ug*skipratio
         tskip = max(tskip, 20.0_WP * time%dtmax)
         ! call param_read('Skip Points', start_idx)
         N = size(all_time, dim=1)
         ! if (N < start_idx) start_idx = 2 ! Fallback for very short simulations
         start_idx = 1
         do i = 1, N
             if (all_time(i) >= tskip) then
                 start_idx = i
                 exit
             end if
         end do
         start_idx=max(start_idx,200)
         count = 0
         sum_x = 0.0_WP; sum_y = 0.0_WP; sum_xy = 0.0_WP; sum_xx = 0.0_WP
         
         ! 2. Filter data and accumulate OLS sums
         do i = start_idx, N
             ! Stop accumulating the moment we hit the non-linear vortex roll-up
             if (all_amp(i) > threshold * dg .or. all_amp_height(i) > threshold * dg ) exit
             
             ! Only use strictly positive amplitudes for the logarithm
             if (all_amp(i) > 0.0_WP) then
                 x = all_time(i)
                 y = log(all_amp(i))
                 
                 sum_x  = sum_x + x
                 sum_y  = sum_y + y
                 sum_xy = sum_xy + x * y
                 sum_xx = sum_xx + x * x
                 count  = count + 1
             end if
         end do
         
         ! 3. Calculate the slope (Growth Rate) analytically
         if (count > 1) then
             slope = (real(count, WP) * sum_xy - sum_x * sum_y) / &
                     (real(count, WP) * sum_xx - sum_x * sum_x)
         else
             slope = 0.0_WP
         end if
         
         ! 4. Output the result to CSV
         lc = dg
         tau = dg / Ug
         
         if (fs%cfg%amRoot) then
             ! Print a diagnostic to the terminal so you can verify the fit quality
            !  print '(A,I6,A,ES12.4)', ' -> OLS Fit applied to ', count, ' points. Raw Slope = ', slope
             print *, start_idx, tskip
             inquire(file='sweep_results.csv', exist=file_exists)
             open(unit=file_unit, file='sweep_results.csv', status='unknown', position='append', action='write')
             if (.not. file_exists) then
                write(file_unit, '(A12,A12)')  'Wavenumber', 'Growth_Rate'
             end if
             ! Write normalized wavenumber and normalized growth rate
             write(file_unit, '(es12.5, es12.5)') alpha*lc, slope*tau
             close(file_unit)
         end if
      end block linear_staiblity
      ! ! Post-process growth rate using ODRPACK
      ! odr_fit: block
      !    use, intrinsic :: iso_fortran_env, only: output_unit
      !    use mathtools, only: twoPi
      !    use messager,  only: log
      !    use string,    only: str_long
      !    character(len=str_long) :: message
      !    integer :: i
      !    ! ODRPACK variables - explicit model based on exponential of time
      !    integer                       :: N                      !> Number of observations (number of polygons)
      !    integer , parameter           :: M=1                    !> Number of elements per explanatory variables (1 time)
      !    integer , parameter           :: NP=2                   !> Number of parameters in our model (2 for a normalized exponential in time with time shift)
      !    integer , parameter           :: NQ=1                   !> Number of response per observation (only 1, the normalized amplitude)
      !    real(WP), dimension(NP)       :: BETA=0.0_WP            !> Array of model parameter values (the growth rate and time shift)
      !    real(WP), dimension(:,:)  , allocatable :: YY           !> Value of response variable (of size LDYYxNQ)
      !    integer                       :: LDYY                   !> Leading dimension of YY (equals N since an explicit model is used)
      !    real(WP), dimension(:,:)  , allocatable :: XX           !> Value of explanatory variable (of size LDXXxM)
      !    integer                       :: LDXX                   !> Leading dimension of XX (equals N)
      !    real(WP), dimension(:,:,:), allocatable :: WE           !> Weighting of response data (of size LDWExLD2WExNQ)
      !    integer                       :: LDWE                   !> Leading dimension of WE (equals N since an explicit model is used)
      !    integer                       :: LD2WE                  !> Second dimension of WE (equals NQ)
      !    real(WP), dimension(:,:,:), allocatable :: WD           !> Weighting of explanatory data (of size LDWDxLD2WDxM)
      !    integer                       :: LDWD                   !> Leading dimension of WD (equals N)
      !    integer                       :: LD2WD                  !> Second dimension of WD (equals 1)
      !    integer , dimension(NP)       :: IFIXB=-1               !> Whether any model parameters has to be kept constant
      !    integer , parameter           :: LDIFX=1                !> Leading dimension of IFIXX (equals 1)
      !    integer , dimension(LDIFX,M)  :: IFIXX=-1               !> Whether any explanatory variable data is to be treated as "fixed"
      !    integer                       :: JOB=00030              !> 5-digit parameter flag that controls execution (this invokes analytical Jacobian with explicit model)
      !    integer                       :: NDIGIT=1               !> Number of reliable digits in our model - let ODRPACK figure it out on its own
      !    real(WP)                      :: TAUFAC=0.0_WP          !> To control size of first step (ignored here)
      !    real(WP)                      :: SSTOL=-1.0_WP          !> Relative cvg of sum of squares: this sets it to 1e-8             ********* Need to change to sth else
      !    real(WP)                      :: PARTOL=-1.0_WP         !> Relative cvg for model parameters: this sets it to 1e-11         ********* Need to change to sth else
      !    integer                       :: MAXIT=-1               !> Maximum number of iterations                                     ********* Need to change to sth else
      !    integer                       :: IPRINT=0               !> 4-digit parameter flag for controlling printing (default is -1)
      !    integer                       :: LUNERR=10              !> Logical unit for error reporting (6 by default)
      !    integer                       :: LUNRPT=10              !> Logical unit for reporting
      !    real(WP), dimension(NP)       :: STPB=0.0_WP            !> Relative step sizes for Jacobian for model parameters (here, default)
      !    integer , parameter           :: LDSTPD=1               !> Leading dimension of STPD, either 1 or N (here, 1)
      !    real(WP), dimension(LDSTPD,1) :: STPD=0.0_WP            !> Relative step sizes for Jacobian for input errors (here, default)
      !    real(WP), dimension(NP)       :: SCLB=1.0_WP            !> Scaling for the model parameters (here, not default but set to 1.0 to avoid rescaling 0 coefficients)
      !    real(WP), dimension(:,:)  , allocatable :: SCLD         !> Scaling for the input errors (here, not default but set to 1.0 to avoid rescaling 0 coefficients)
      !    integer                       :: LDSCLD                 !> Leading dimension of SCLD, either 1 or N (here, N)
      !    integer                       :: LWORK                  !> Size of WORK array
      !    real(WP), dimension(:)    , allocatable :: WORK         !> WORK array
      !    integer , parameter           :: LiWORK=20+NP+NQ*(NP+M) !> Size of IWORK array
      !    integer , dimension(LiWORK)   :: iWORK                  !> iWORK array
      !    integer                       :: INFO                   !> Why the calculations stopped
      !    real(WP)                      :: lc, tau                ! <--- ADD THIS LINE
      !    integer :: file_unit
      !    logical :: file_exists
      !    file_unit = 42
      !    ! Copy over data and sizes
      !    N=size(all_time,dim=1)
      !    amp0=all_amp(1)
      !    LDYY=N; allocate(YY(LDYY,NQ)); YY(:,1)=all_amp/amp0
      !    LDXX=N; allocate(XX(LDXX,M )); XX(:,1)=all_time
      !    LDWE=N; LD2WE=NQ; allocate(WE(LDWE,LD2WE,NQ)); WE=1.0_WP
      !    LDWD=N; LD2WD=1 ; allocate(WD(LDWD,LD2WD,M )); WD=1.0_WP
      !    LDSCLD=N; allocate(SCLD(LDSCLD,M)); SCLD=1.0_WP
      !    LWORK=18+11*NP+NP**2+M+M**2+4*N*NQ+6*N*M+2*N*NQ*NP+2*N*NQ*M+NQ**2+5*NQ+NQ*(NP+M)+(LDWE*LD2WE)*NQ; allocate(WORK(LWORK))
      !    ! Call ODRPACK to find time shift
      !    call DODRC(exponential_model,N,M,NP,NQ,BETA,YY,LDYY,XX,LDXX,WE,LDWE,LD2WE,WD,LDWD,LD2WD,IFIXB,IFIXX,LDIFX,JOB,NDIGIT,TAUFAC,&
      !    &          SSTOL,PARTOL,MAXIT,IPRINT,LUNERR,LUNRPT,STPB,STPD,LDSTPD,SCLB,SCLD,LDSCLD,WORK,LWORK,iWORK,LiWORK,INFO)
      !    ! Adjust weights to eliminate the early non-exponential part
      !    print *, "Early non-exponetial time", 2.0_WP*BETA(2)
      !    do i=1,size(all_time,dim=1)
      !       if (all_time(i).le.2.0_WP*BETA(2)) then
      !          WE(i,1,1)=0.0_WP
      !          WD(i,1,1)=0.0_WP
      !       else if (all_amp(i) > 0.3_WP * dg) then
      !          WE(i,1,1) = 0.0_WP
      !          WD(i,1,1) = 0.0_WP   
      !       ! 3. Keep everything in between
      !       else
      !          WE(i,1,1) = 1.0_WP
      !          WD(i,1,1) = 1.0_WP
      !       end if
      !    end do
         
      !    ! Call ODRPACK again to find growth rate
      !    call DODRC(exponential_model,N,M,NP,NQ,BETA,YY,LDYY,XX,LDXX,WE,LDWE,LD2WE,WD,LDWD,LD2WD,IFIXB,IFIXX,LDIFX,JOB,NDIGIT,TAUFAC,&
      !    &          SSTOL,PARTOL,MAXIT,IPRINT,LUNERR,LUNRPT,STPB,STPD,LDSTPD,SCLB,SCLD,LDSCLD,WORK,LWORK,iWORK,LiWORK,INFO)
      !    lc=dg
      !    tau=dg/Ug
      !    ! Get back growth rate
      !    if (fs%cfg%amRoot) then
      !       inquire(file='sweep_results.csv', exist=file_exists)
      !       open(unit=file_unit, file='sweep_results.csv', status='unknown', position='append', action='write')
      !       if (.not. file_exists) then
      !          write(file_unit, '(A12,A12)')  'Wavenumber', 'Growth_Rate'
      !       end if
      !       write(file_unit, '(es12.5, es12.5)') alpha*lc, BETA(1)*tau
      !       close(file_unit)
      !       ! write(output_unit,'(es12.5,x,es12.5,x,es12.5,x,es12.5)') lc,tau,twoPi/fs%cfg%xL*lc,BETA(1)*tau
      !       ! write(message    ,'("Reference time scale   = ",es12.5)') tau               ; call log(message)
      !       ! write(message    ,'("Cut-off length scale   = ",es12.5)') lc                ; call log(message)
      !       ! write(message    ,'("Normalized growth rate = ",es12.5)') BETA(1)*tau       ; call log(message)
      !       ! write(message    ,'("Normalized wave number = ",es12.5)') twoPi/fs%cfg%xL*lc; call log(message)
      !    end if
      ! end block odr_fit
      
      ! odr_hybrid_fit: block
      !    use, intrinsic :: iso_fortran_env, only: output_unit
      !    use mathtools, only: twoPi
      !    use string,    only: str_long
      !    character(len=str_long) :: message
      !    integer :: i, start_idx, end_step, step_T, num_points
      !    ! ODRPACK variables
      !    integer                       :: N, M=1, NP=2, NQ=1
      !    real(WP), dimension(2)        :: BETA=0.0_WP
      !    real(WP), dimension(:,:)  , allocatable :: YY, XX
      !    integer                       :: LDYY, LDXX
      !    real(WP), dimension(:,:,:), allocatable :: WE, WD
      !    integer                       :: LDWE, LD2WE, LDWD, LD2WD
      !    integer , dimension(2)        :: IFIXB=-1
      !    integer , parameter           :: LDIFX=1
      !    integer , dimension(1,1)      :: IFIXX=-1
      !    integer                       :: JOB=00030, NDIGIT=1
      !    real(WP)                      :: TAUFAC=0.0_WP, SSTOL=-1.0_WP, PARTOL=-1.0_WP
      !    integer                       :: MAXIT=-1, IPRINT=0, LUNERR=10, LUNRPT=10
      !    real(WP), dimension(2)        :: STPB=0.0_WP
      !    integer , parameter           :: LDSTPD=1
      !    real(WP), dimension(1,1)      :: STPD=0.0_WP
      !    real(WP), dimension(2)        :: SCLB=1.0_WP
      !    real(WP), dimension(:,:)  , allocatable :: SCLD
      !    integer                       :: LDSCLD, LWORK, INFO
      !    real(WP), dimension(:)    , allocatable :: WORK
      !    integer , parameter           :: LiWORK=20+2+1*(2+1) ! Size is 25 for NP=2, M=1, NQ=1
      !    integer , dimension(25)       :: iWORK
         
      !    ! Literature R^2 variables
      !    real(WP) :: mean_t, mean_logb, dt, dlogb, num, den_t, den_logb, R2, R2_prev
      !    real(WP), allocatable :: log_beta(:)
         
      !    ! Output variables
      !    real(WP) :: lc, tau
      !    integer :: file_unit = 42
      !    logical :: file_exists

      !    N = size(all_time, dim=1)
         
      !    ! Safety fallback for initial amplitude (prevents divide-by-zero or log(0) crash)
      !    amp0 = all_amp(1)
      !    if (amp0 <= 0.0_WP) then
      !        ! Find the first non-zero amplitude to use as normalization
      !        do i = 1, N
      !            if (all_amp(i) > 0.0_WP) then
      !                amp0 = all_amp(i)
      !                exit
      !            end if
      !        end do
      !    end if

      !    ! Allocate ODR arrays
      !    LDYY=N; allocate(YY(LDYY,NQ)); YY(:,1)=all_amp/amp0
      !    LDXX=N; allocate(XX(LDXX,M )); XX(:,1)=all_time
      !    LDWE=N; LD2WE=NQ; allocate(WE(LDWE,LD2WE,NQ)); WE=1.0_WP
      !    LDWD=N; LD2WD=1 ; allocate(WD(LDWD,LD2WD,M )); WD=1.0_WP
      !    LDSCLD=N; allocate(SCLD(LDSCLD,M)); SCLD=1.0_WP
         
      !    LWORK=18+11*NP+NP**2+M+M**2+4*N*NQ+6*N*M+2*N*NQ*NP+2*N*NQ*M+NQ**2+5*NQ+NQ*(NP+M)+(LDWE*LD2WE)*NQ
      !    allocate(WORK(LWORK))

      !    ! =================================================================
      !    ! STEP 1: First ODRPACK Pass (Find initial time shift BETA(2))
      !    ! =================================================================
      !    call DODRC(exponential_model,N,M,NP,NQ,BETA,YY,LDYY,XX,LDXX,WE,LDWE,LD2WE,WD,LDWD,LD2WD,IFIXB,IFIXX,LDIFX,JOB,NDIGIT,TAUFAC,&
      !    &          SSTOL,PARTOL,MAXIT,IPRINT,LUNERR,LUNRPT,STPB,STPD,LDSTPD,SCLB,SCLD,LDSCLD,WORK,LWORK,iWORK,LiWORK,INFO)

      !    ! =================================================================
      !    ! STEP 2: Find the start index (Skipping the transient)
      !    ! =================================================================
      !    start_idx = 1
      !    do i = 1, N
      !       if (all_time(i) > 2.0_WP * BETA(2)) then
      !          start_idx = i
      !          exit
      !       end if
      !    end do

      !    ! =================================================================
      !    ! STEP 3: R^2 Logic (Find the end of the linear regime)
      !    ! =================================================================
      !    allocate(log_beta(N))
      !    do i = 1, N
      !       if (all_amp(i) > 0.0_WP) then
      !          log_beta(i) = log(all_amp(i))
      !       ! else
      !       !    log_beta(i) = -99.0_WP ! Safe arbitrary low value for log(0) cases
      !       end if
      !    end do

      !    R2_prev = -1.0_WP
      !    end_step = N

      !    do step_T = start_idx + 2, N
      !       num_points = (step_T - start_idx) + 1
            
      !       mean_t    = sum(all_time(start_idx:step_T)) / real(num_points, kind=WP)
      !       mean_logb = sum(log_beta(start_idx:step_T)) / real(num_points, kind=WP)

      !       num      = 0.0_WP
      !       den_t    = 0.0_WP
      !       den_logb = 0.0_WP

      !       do i = start_idx, step_T
      !          dt = all_time(i) - mean_t
      !          dlogb = log_beta(i) - mean_logb
      !          num = num + (dt * dlogb)
      !          den_t = den_t + (dt**2)
      !          den_logb = den_logb + (dlogb**2)
      !       end do

      !       if (den_t * den_logb > 0.0_WP) then
      !          R2 = (num**2) / (den_t * den_logb)
      !       else
      !          R2 = 0.0_WP
      !       end if

      !       ! Trigger the cutoff if R^2 drops
      !       if (num_points > 3 .and. R2 < R2_prev) then
      !          end_step = step_T - 1
      !          exit
      !       end if

      !       R2_prev = R2
      !    end do

      !    ! =================================================================
      !    ! STEP 4: Second ODRPACK Pass (Isolate the valid regime)
      !    ! =================================================================
      !    ! Turn off the weighting for points outside our newly found [start_idx, end_step] bounds
      !    do i = 1, N
      !       if (i < start_idx .or. i > end_step) then
      !          WE(i,1,1) = 0.0_WP
      !          WD(i,1,1) = 0.0_WP
      !       else
      !          WE(i,1,1) = 1.0_WP
      !          WD(i,1,1) = 1.0_WP
      !       end if
      !    end do

      !    call DODRC(exponential_model,N,M,NP,NQ,BETA,YY,LDYY,XX,LDXX,WE,LDWE,LD2WE,WD,LDWD,LD2WD,IFIXB,IFIXX,LDIFX,JOB,NDIGIT,TAUFAC,&
      !    &          SSTOL,PARTOL,MAXIT,IPRINT,LUNERR,LUNRPT,STPB,STPD,LDSTPD,SCLB,SCLD,LDSCLD,WORK,LWORK,iWORK,LiWORK,INFO)

      !    ! =================================================================
      !    ! STEP 5: Data Output
      !    ! =================================================================
      !    lc = dg
      !    tau = dg / Ug
         
      !    if (fs%cfg%amRoot) then
      !       inquire(file='sweep_results.csv', exist=file_exists)
      !       open(unit=file_unit, file='sweep_results.csv', status='unknown', position='append', action='write')
            
      !       if (.not. file_exists) then
      !          write(file_unit, '(A12,A12)')  'Wavenumber', 'Growth_Rate'
      !       end if
            
      !       write(file_unit, '(es12.5, es12.5)') alpha*lc, BETA(1)*tau
      !       close(file_unit)
      !    end if

      !    ! Clean up memory
      !    deallocate(log_beta, YY, XX, WE, WD, SCLD, WORK)
         
      ! end block odr_hybrid_fit
   end subroutine simulation_run
   
   
   !> Definition of our exponential function of time model
   subroutine exponential_model(N,M,NP,NQ,LDN,LDM,LDNP,BETA,XPLUSD,IFIXB,IFIXX,LDFIX,IDEVAL,F,FJACB,FJACD,ISTOP)
      implicit none
      ! Input parameters
      integer , intent(in) :: IDEVAL,LDFIX,LDM,LDN,LDNP,M,N,NP,NQ
      integer , dimension(NP)     , intent(in) :: IFIXB
      integer , dimension(LDFIX,M), intent(in) :: IFIXX
      real(WP), dimension(NP)     , intent(in) :: BETA
      real(WP), dimension(LDN,M)  , intent(in) :: XPLUSD
      ! Output parameters
      real(WP), dimension(LDN,NQ) :: F
      real(WP), dimension(LDN,LDNP,NQ) :: FJACB
      real(WP), dimension(LDN,LDM ,NQ) :: FJACD
      integer :: ISTOP,i
      ! Check stopping condition - all values are acceptable
      ISTOP=0
      ! Compute model value
      if (mod(IDEVAL,10).ge.1) then
         do i=1,N
            F(i,1)=exp(BETA(1)*(XPLUSD(i,1)-BETA(2)))
         end do
      end if
      ! Compute model derivatives with respect to BETA
      if (mod(IDEVAL/10,10).GE.1) then
         do i=1,N
            FJACB(i,1,1)=(XPLUSD(i,1)-BETA(2))*exp(BETA(1)*(XPLUSD(i,1)-BETA(2)))
            FJACB(i,2,1)=            -BETA(1) *exp(BETA(1)*(XPLUSD(i,1)-BETA(2)))
         end do
      end if
      ! Compute model derivatives with respect to input
      if (mod(IDEVAL/100,10).GE.1) then
         do i=1,N
            FJACD(i,1,1)=BETA(1)*exp(BETA(1)*(XPLUSD(i,1)-BETA(2)))
         end do
      end if
   end subroutine exponential_model

   subroutine get_KE()
      implicit none
      real(WP), dimension(:,:,:), allocatable:: Utmp 
      real(WP), dimension(:,:,:), allocatable:: Vtmp 
      real(WP), dimension(:,:,:), allocatable:: Wtmp 
      real(WP), dimension(:,:,:), allocatable:: Stmp,Stmp2 
      integer :: i,j,k,ii,jj,kk
      allocate(Utmp(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(Vtmp(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(Wtmp(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(Stmp(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(Stmp2(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))

      ! Get KE 
      calculate_KE: block
         Utmp=fs%rho_U*fs%U*fs%U/2.0_WP
         Vtmp=fs%rho_V*fs%V*fs%V/2.0_WP
         Wtmp=fs%rho_W*fs%W*fs%W/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))+&
                  &           sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))+&
                  &           sum(fs%itpw_z(:,i,j,k)*Wtmp(i,j,k:k+1))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE)
      end block calculate_KE

      ! Get KE_l
      calculate_KE_l: block
         Utmp=fs%rho_Ul*fs%U*fs%U/2.0_WP
         Vtmp=fs%rho_Vl*fs%V*fs%V/2.0_WP
         Wtmp=fs%rho_Wl*fs%W*fs%W/2.0_WP
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
      end block calculate_KE_l

      ! Get KE_g
      calculate_KE_g: block
         Utmp=fs%rho_Ug*fs%U*fs%U/2.0_WP
         Vtmp=fs%rho_Vg*fs%V*fs%V/2.0_WP
         Wtmp=fs%rho_Wg*fs%W*fs%W/2.0_WP
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
      end block calculate_KE_g

      ! Get KE 
      calculate_KEx: block
         Utmp=fs%rho_U*fs%U*fs%U/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_x)
      end block calculate_KEx

      ! Get KE_l
      calculate_KE_xl: block
         Utmp=fs%rho_Ul*fs%U*fs%U/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_xl)
      end block calculate_KE_xl

      ! Get KE_g
      calculate_KE_xg: block
         Utmp=fs%rho_Ug*fs%U*fs%U/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_xg)
      end block calculate_KE_xg


      ! Get KE 
      calculate_KEy: block
         Vtmp=fs%rho_V*fs%V*fs%V/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_y)
      end block calculate_KEy

      ! Get KE_l
      calculate_KE_ly: block
         Vtmp=fs%rho_Vl*fs%V*fs%V/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_yl)
      end block calculate_KE_ly

      ! Get KE_g
      calculate_KE_gy: block
         Vtmp=fs%rho_Vg*fs%V*fs%V/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_yg)
      end block calculate_KE_gy

      
       ! Get total momentum
      calculate_rhoU: block
         stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*(fs%rho_U(i:i+1,j,k)*fs%U(i:i+1,j,k)))
                  Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*(fs%rho_V(i,j:j+1,k)*fs%V(i,j:j+1,k)))
                  Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*(fs%rho_W(i,j,k:k+1)*fs%W(i,j,k:k+1)))
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoUInt)
         call fs%cfg%integrate(Vtmp,integral=rhoVInt)
         call fs%cfg%integrate(Wtmp,integral=rhoWInt)
      end block calculate_rhoU 

      calculate_rhoU_l: block
         stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*(fs%rho_Ul(i:i+1,j,k)*fs%U(i:i+1,j,k)))
                  Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*(fs%rho_Vl(i,j:j+1,k)*fs%V(i,j:j+1,k)))
                  Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*(fs%rho_Wl(i,j,k:k+1)*fs%W(i,j,k:k+1)))
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoU_l)
         call fs%cfg%integrate(Vtmp,integral=rhoV_l)
         call fs%cfg%integrate(Wtmp,integral=rhoW_l)
      end block calculate_rhoU_l 

      calculate_rhoU_g: block
         stmp=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*(fs%rho_Ug(i:i+1,j,k)*fs%U(i:i+1,j,k)))
                  Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*(fs%rho_Vg(i,j:j+1,k)*fs%V(i,j:j+1,k)))
                  Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*(fs%rho_Wg(i,j,k:k+1)*fs%W(i,j,k:k+1)))
               end do
            end do
         end do
         call fs%cfg%integrate(Utmp,integral=rhoU_g)
         call fs%cfg%integrate(Vtmp,integral=rhoV_g)
         call fs%cfg%integrate(Wtmp,integral=rhoW_g)
      end block calculate_rhoU_g

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
                     Utmp(i,j,k)=sum(fs%grdw_y(:,i,j,k)*fs%W(i,j-1:j,k))-sum(fs%grdv_z(:,i,j,k)*fs%V(i,j,k-1:k))
                     Vtmp(i,j,k)=sum(fs%grdu_z(:,i,j,k)*fs%U(i,j,k-1:k))-sum(fs%grdw_x(:,i,j,k)*fs%W(i-1:i,j,k))
                     Wtmp(i,j,k)=sum(fs%grdv_x(:,i,j,k)*fs%V(i-1:i,j,k))-sum(fs%grdu_y(:,i,j,k)*fs%U(i,j-1:j,k))
                  end do
               end do
            end do
            call cfg%sync(Utmp)
            call cfg%sync(Vtmp)
            call cfg%sync(Wtmp)
            Stmp=0.0_WP
            do k = cfg%kmin_, cfg%kmax_
               do j = cfg%jmin_, cfg%jmax_
                  do i = cfg%imin_, cfg%imax_
            
                     ! omega_x at cell center: average of 4 edges
                     wxc = 0.25_wp * ( Utmp(i,   j  ,   k  )   &
                                    +   Utmp(i,   j+1,   k  ) &
                                    +   Utmp(i,   j  ,   k+1) &
                                    +   Utmp(i,   j+1,   k+1) )

                     ! omega_y at cell center: average of 4 edges
                     wyc = 0.25_wp * ( Vtmp(i,   j  ,   k  )   &
                                    +   Vtmp(i+1, j  ,   k  ) &
                                    +   Vtmp(i,   j  ,   k+1) &
                                    +   Vtmp(i+1, j  ,   k+1) )
            
                     ! omega_z at cell center: average of 4 edges
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
         ! Open file dynamically with append mode
         open(unit=10, file=filename, status="unknown", position="append", action="write")
         ! Write data with 16-digit precision in CSV format
         write(10, '(F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16,F24.16)')&
         time%t,time%dt,KE,KE_l,KE_g,rhoUInt,rhoVInt,rhoWInt,rhoU_l,rhoV_l,rhoW_l,rhoU_g,rhoV_g,rhoW_g,PE_l,PE_g,EN_l,EN_g,KE_x,KE_xl,KE_xg,KE_y,KE_yl,KE_yg
         close(10)
      end if

      deallocate(Utmp,Vtmp,Wtmp,stmp)
   end subroutine get_KE
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   

   
   
end module simulation

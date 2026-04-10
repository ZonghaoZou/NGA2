!> Various definitions and tools for running an NGA2 simulation
!> Rayleigh-Taylor instability — Scheme 2: KE Conservative (Staggered Grid)
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
   real(WP), dimension(:,:,:,:), allocatable :: vel
   
   real(WP), dimension(:), allocatable :: all_time, all_amp,all_amp_height
   real(WP) :: amp, grate, amp0, amp_height
   integer :: nx
   !> Problem definition: KH instability
   real(WP) :: Ll,Lg,Ug,Ul,dg,dl,alpha,Lx
   real(WP) :: eps_kh 


   real(WP) :: KE,VFInt,rhoInt,rhoUInt,rhoVInt,rhoWInt
   real(WP) :: KE_l,KE_g,rhoU_l,rhoV_l,rhoW_l,rhoU_g,rhoV_g,rhoW_g
   real(WP) :: EN_l,EN_g,PE_l,PE_g
   real(WP) :: KE_x,KE_y,KE_xl,KE_xg,KE_yl,KE_yg
   character(len=20) :: filename='output.csv'
   
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

      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(vel (1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Two-phase NS')
         fs%theta=0.5_WP+0.01_WP
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
         resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_density(rho=resU)

         init_velocity: block
            use orr_sommerfeld_kh, only: eval_kh
            real(WP) :: u_tmp, v_tmp, G_tmp
            integer :: i,j,k
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
         end block init_velocity
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%apply_bcond(time%t,time%dt)
         ! Poisson equation ================================================
         call fs%get_Umid()
         call fs%update_laplacian()
         call fs%correct_mfr()
         call fs%get_div()
         fs%psolv%rhs=-fs%cfg%vol*fs%div
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         fs%Umid=fs%Umid-resU/fs%sRHOX**2
         fs%Vmid=fs%Vmid-resV/fs%sRHOY**2
         fs%Wmid=fs%Wmid-resW/fs%sRHOZ**2
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         ! call fs%interp_vel(vel(1,:,:,:),vel(2,:,:,:),vel(3,:,:,:))
         call fs%get_div()
         ! if (cfg%amRoot) then
         !    open(unit=10, file=filename, status="replace", action="write")
         !    write(10, '(A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24)') &
         !    & 'time', 'timestep', 'KE', 'KE_l', 'KE_g', 'rhoU', 'rhoV', 'rhoW', 'rhoU_l', 'rhoV_l', 'rhoW_l','rhoU_g', 'rhoV_g', 'rhoW_g','PE_l','PE_g','EN_l','EN_g'
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
         call cflfile%add_column(fs%CFLv_z,'Convective zCFL')
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
         call eval_kh(xyz(1), xyz(2), eps_kh, u_tmp, v_tmp, G)
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
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation — KE Conservative staggered scheme
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
            
            ! Apply boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Poisson equation ================================================
            call fs%get_Umid()
            
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
            fs%P=fs%P+fs%psolv%sol
            fs%Umid=fs%Umid-time%dt*resU/((fs%sRHOX+fs%sRHOXold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOX)
            fs%Vmid=fs%Vmid-time%dt*resV/((fs%sRHOY+fs%sRHOYold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOY)
            fs%Wmid=fs%Wmid-time%dt*resW/((fs%sRHOZ+fs%sRHOZold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOZ)
            
            ! Regenerate U from Umid and Uold
            call fs%get_U()
            
            ! Increment sub-iteration counter =================================
            time%it=time%it+1
         end do
         
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         ! call fs%interp_vel(vel(1,:,:,:),vel(2,:,:,:),vel(3,:,:,:))
         call fs%get_div()
         
          ! Compute stats
         ! call get_stats()
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
      real(WP), dimension(:,:,:), allocatable:: Stmp,Stmp2 
      allocate(Utmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Utmp=0.0_WP
      allocate(Vtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Vtmp=0.0_WP
      allocate(Wtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Wtmp=0.0_WP
      allocate(Stmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Stmp=0.0_WP
      allocate(FX(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FX=0.0_WP
      allocate(FY(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FY=0.0_WP
      allocate(FZ(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));FZ=0.0_WP
      allocate(Stmp2(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
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
         call cfg%integrate(Stmp,integral=KE)
      end block calculate_KE

      calculate_KEx: block
         Utmp=(fs%sRHOX**2)*fs%U*fs%U/2.0_WP
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))
               end do
            end do
         end do
         call cfg%integrate(Stmp,integral=KE_x)
      end block calculate_KEx


      calculate_KEy: block
         Vtmp=(fs%sRHOY**2)*fs%V*fs%V/2.0_WP
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))
               end do
            end do
         end do
         call cfg%integrate(Stmp,integral=KE_y)
      end block calculate_KEy


      Utmp=0.0_WP;Vtmp=0.0_WP;Wtmp=0.0_WP
      calculate_momentumconservation: block
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Utmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*(fs%sRHOX(i:i+1,j,k)**2*fs%U(i:i+1,j,k)))
                  Vtmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*(fs%sRHOY(i,j:j+1,k)**2*fs%V(i,j:j+1,k)))
                  Wtmp(i,j,k)=sum(fs%itpw_z(:,i,j,k)*(fs%sRHOZ(i,j,k:k+1)**2*fs%W(i,j,k:k+1)))
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


      ! Get KE_l
      calculate_l_KEx: block
         Utmp=FX*fs%U*fs%U/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_xl)
      end block calculate_l_KEx


      ! Get KE_l
      calculate_l_KEy: block
         Vtmp=FY*fs%V*fs%V/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_yl)
      end block calculate_l_KEy


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

      calculate_g_KEx: block
         Utmp=FX*fs%U*fs%U/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Utmp(i:i+1,j,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_xg)
      end block calculate_g_KEx


      calculate_g_KEy: block
         Vtmp=FY*fs%V*fs%V/2.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Stmp(i,j,k)=sum(fs%itpv_y(:,i,j,k)*Vtmp(i,j:j+1,k))
               end do
            end do
         end do
         call fs%cfg%integrate(Stmp,integral=KE_yg)
      end block calculate_g_KEy

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

   end subroutine get_KE

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      deallocate(resU,resV,resW,Ui,Vi,Wi,vel)
   end subroutine simulation_final
   
   
end module simulation

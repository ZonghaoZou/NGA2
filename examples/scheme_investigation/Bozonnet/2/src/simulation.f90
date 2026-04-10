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
   use iterator_class,    only: iterator
   use pardata_class,     only: pardata
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
   
   type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
   real(WP) :: vof_removed              !< Integral of VOF removed
   integer  :: nlayer=10                !< Size of buffer layer for VOF removal

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   !> Provide a pardata object for restarts
   logical       :: restarted  !< Is the simulation restarted?
   type(pardata) :: df         !< Pardata object for restart I/O
   type(event)   :: save_evt   !< Event to trigger restart I/O

   public :: simulation_init,simulation_run,simulation_final,get_KE
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: vel
   
   real(WP) :: Hg,Hl,Ug,Ul,dg,dl,deficit,Uint,a_vel,Ucf,Usm
   ! For data processing
   real(WP) :: tauflash,t_last,t_start

   real(WP) :: KE,VFInt,rhoInt,rhoUInt,rhoVInt,rhoWInt
   real(WP) :: KE_l,KE_g,rhoU_l,rhoV_l,rhoW_l,rhoU_g,rhoV_g,rhoW_g
   real(WP) :: EN_l,EN_g,PE_l,PE_g
   real(WP) :: KE_x,KE_y,KE_xl,KE_xg,KE_yl,KE_yg
   character(len=20) :: filename='output.csv'
   
contains
   
   subroutine processdata
      use mpi_f08, only: MPI_ALLREDUCE, MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i, j, k, ierr
      real(WP), dimension(:), allocatable :: h_local, h_global
      if (time%t.gt.tauflash.and.time%t - t_last>=1.0e-4) then
         if (t_start==0.0_WP) t_start=time%t
         allocate(h_local (1:cfg%nx));h_local =0.0_WP
         allocate(h_global(1:cfg%nx));h_global=0.0_WP
         do k = vf%cfg%kmin_, vf%cfg%kmax_
            do i = vf%cfg%imin_, vf%cfg%imax_
               do j = vf%cfg%jmin_, vf%cfg%jmax_
                     if (vf%cfg%z(k).le.0.0_WP.and.vf%cfg%z(k+1).gt.0.0_WP) then
                        h_local(i)=h_local(i)+vf%VF(i,j,k)*vf%cfg%dy(j)
                     end if
               end do
            end do
         end do
         call MPI_ALLREDUCE(h_local,h_global,cfg%nx,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
         if (cfg%amRoot) then
            open(unit=20, file='interface_height.csv', position='append', status='unknown')
            write(20, *) time%t, h_global(1:cfg%nx)
            close(20)
         end if
         t_last=time%t
         if (time%t-t_start>=1.2_WP*40_WP/33.0_WP) time%tmax=time%t
         deallocate(h_local,h_global)
      end if   
   end subroutine processdata


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
         ! Read in liquid height for levelset
         call param_read('Liquid height', Hl)
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_AWML,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
                  if (i.ge.cfg%imin+2) then
                     vf%VF(i,j,k)=0.0_WP
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
         ! Now apply Neumann condition on interface at inlet to have proper round injection
         neumann_irl: block
            use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
            &                                setNumberOfPlanes,setPlane,matchVolumeFraction
            real(WP), dimension(1:4) :: plane
            type(RectCub_type) :: cell
            call new(cell)
            if (vf%cfg%iproc.eq.1) then
               do k=vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i=vf%cfg%imino,vf%cfg%imin-1
                        ! Extract plane data and copy in overlap
                        plane=getPlane(vf%liquid_gas_interface(vf%cfg%imin,j,k),0)
                        call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                        &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
                        plane(4)=dot_product(plane(1:3),[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                        call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
                     end do
                  end do
               end do
            end if
         end block neumann_irl
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
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: dirichlet,slip,clipped_neumann
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Two-phase NS')
         fs%theta=0.5_WP+0.01_WP
         ! Assign constant viscosity to each phase
         call param_read('Gas dynamic viscosity',fs%visc_g)
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         ! Assign constant density to each phase
         call param_read('Gas density',fs%rho_g)
         call param_read('Liquid density',fs%rho_l)
         ! Read in surface tension coefficient 
         call param_read('Surface tension coefficient',fs%sigma)
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Slip at the top, no-slip wall is implied via cfg%VF 
         call fs%add_bcond(name='bc_yp'  ,type=slip           ,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         ! Dirichlet inflow
         call fs%add_bcond(name='inlow'  ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Clipped-Neumann outflow
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
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


      initialize_velocity: block
         use tpns_class, only: bcond
         use random,   only: random_uniform
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k,ierr
         real(WP) :: yval,Ud
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_density(rho=resU)
         call param_read('Gas height', Hg); call param_read('Liquid vorticity thickness', dl); call param_read('Deficit parameter', deficit)
         call param_read('Velocity disturbance', a_vel);  call param_read('Liquid velocity',Ul); call param_read('Gas velocity',Ug)
         Ucf=0.1_WP*Ug; Usm=(Ucf+Ug)/2.0_WP; dg=6.0_WP*Hg/sqrt(fs%rho_g*Ug*Hg/fs%visc_g)
         Uint=dl*deficit*(Ug*fs%visc_g/dg + Ul*fs%visc_l/dl)/(fs%visc_l+fs%visc_g)
         Ud=(sqrt(fs%rho_l)*Ul+sqrt(fs%rho_g)*Ug)/(sqrt(fs%rho_l)+sqrt(fs%rho_g))
         tauflash=3.0_WP*dg*240.0_WP/Ud; t_last=0.0_WP
         call fs%get_bcond('inlow',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            yval=cfg%ym(j)
            if (yval.le.Hl) then
               fs%U(i,j,k)=(Ul*erf((Hl-yval)/dl) + Uint*(1.0_WP + erf((yval-Hl)/(dl*deficit))))*erf(yval/dl)
            else if (yval.gt.Hl.and.yval.le.Hl+Hg) then
               fs%U(i,j,k)=(Ug*erf((yval-Hl)/dg) + Uint*(1.0_WP - erf((yval-Hl)/(dl*deficit))))* &
               &           (erf((Hg+Hl-yval)/dg))+Usm*(1.0_WP +erf((yval-(Hg+Hl))/dg))
            else
               fs%U(i,j,k)=Usm*(1-erf((yval-(Hl+Hg))/dg)) + Ucf*erf((yval-(Hl+Hg))/dg)
            end if
            if (abs(yval - Hl).le.0.5_WP*dg) fs%U(i,j,k) = fs%U(i,j,k)+random_uniform(lo=-0.5_WP*a_vel,hi=+0.5_WP*a_vel)
            if (abs(cfg%y(j) - Hl).le.0.5_WP*dg) fs%V(i,j,k) = random_uniform(lo=-0.5_WP*a_vel,hi=+0.5_WP*a_vel)       
         end do
         do k=fs%cfg%kmin_,fs%cfg%kmax_+1
            do j=fs%cfg%jmin_,fs%cfg%jmax_+1
               do i=fs%cfg%imin_,fs%cfg%imax_+1
                  if (fs%mask(i,j-1,k).eq.1.and.fs%mask(i-1,j-1,k).eq.1) fs%U(i,j-1,k)=0.0_WP
                  if (fs%mask(i,j-1,k).eq.1.and.fs%mask(i,j-1,k-1).eq.1) fs%W(i,j-1,k)=0.0_WP
               end do 
            end do
         end do 
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
         call fs%get_div()
         ! if (cfg%amRoot) then
         !    open(unit=10, file=filename, status="replace", action="write")
         !    write(10, '(A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24,A24)') &
         !    & 'time', 'timestep', 'KE', 'KE_l', 'KE_g', 'rhoU', 'rhoV', 'rhoW', 'rhoU_l', 'rhoV_l', 'rhoW_l','rhoU_g', 'rhoV_g', 'rhoW_g','PE_l','PE_g','EN_l','EN_g'
         !    close(unit=10)
         ! end if
         ! call get_KE()
      end block initialize_velocity
      
      handle_restart: block
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         use tpns_class,            only: bcond
         use random,                only: random_uniform
         character(len=str_medium) :: filename_tmp
         integer, dimension(3) :: iopartition
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         integer :: i,j,k,n
         type(bcond), pointer :: mybc
         real(WP) :: yval
         ! Create event for saving restart files
         save_evt=event(time,'Restart output')
         call param_read('Restart output period',save_evt%tper)
         ! Read in the I/O partition
         call param_read('I/O partition',iopartition)
         ! Check if we are restarting
         call param_read('Restart from',filename_tmp,default='')
         restarted=.false.; if (len_trim(filename_tmp).gt.0) restarted=.true.
         ! Perform pardata initialization
         if (restarted) then
            ! Read in the file
            call df%initialize(pg=cfg,iopartition=iopartition,fdata=trim(filename_tmp))
            ! Read in the planes directly and set the IRL interface
            allocate(P11(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); call df%pull(name='P11',var=P11)
            allocate(P12(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); call df%pull(name='P12',var=P12)
            allocate(P13(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); call df%pull(name='P13',var=P13)
            allocate(P14(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); call df%pull(name='P14',var=P14)
            do k=vf%cfg%kmin_,vf%cfg%kmax_
               do j=vf%cfg%jmin_,vf%cfg%jmax_
                  do i=vf%cfg%imin_,vf%cfg%imax_
                     call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                     call setPlane(vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                  end do
               end do
            end do
            call vf%sync_interface()
            deallocate(P11,P12,P13,P14)
            ! Reset moments
            call vf%reset_volume_moments()
            ! Update the band
            call vf%update_band()
            ! Create discontinuous polygon mesh from IRL interface
            call vf%polygonalize_interface()
            ! Calculate distance from polygons
            call vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call vf%subcell_vol()
            ! Calculate curvature
            call vf%get_curvature()
            ! Recalculate density
            resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_density(rho=resU)
            ! Now read in the velocity solver data
            call df%pull(name='U',var=fs%U)
            call df%pull(name='V',var=fs%V)
            call df%pull(name='W',var=fs%W)
            call df%pull(name='Umid',var=fs%Umid)
            call df%pull(name='Vmid',var=fs%Vmid)
            call df%pull(name='Wmid',var=fs%Wmid)
            call df%pull(name='P',var=fs%P)
            call df%pull(name='Pjx',var=fs%Pjx)
            call df%pull(name='Pjy',var=fs%Pjy)
            call df%pull(name='Pjz',var=fs%Pjz)
            ! Re-apply inflow velocity profile
            call fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%n_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               yval=cfg%ym(j)
               if (yval.le.Hl) then
                  fs%U(i,j,k)=(Ul*erf((Hl-yval)/dl) + Uint*(1.0_WP + erf((yval-Hl)/(dl*deficit))))*erf(yval/dl)
               else if (yval.gt.Hl.and.yval.le.Hl+Hg) then
                  fs%U(i,j,k)=(Ug*erf((yval-Hl)/dg) + Uint*(1.0_WP - erf((yval-Hl)/(dl*deficit))))* &
                  &           (erf((Hg+Hl-yval)/dg))+Usm*(1.0_WP +erf((yval-(Hg+Hl))/dg))
               else
                  fs%U(i,j,k)=Usm*(1-erf((yval-(Hl+Hg))/dg)) + Ucf*erf((yval-(Hl+Hg))/dg)
               end if
               if (abs(yval - Hl).le.0.5_WP*dg) fs%U(i,j,k) = fs%U(i,j,k)+random_uniform(lo=-0.5_WP*a_vel,hi=+0.5_WP*a_vel)
               if (abs(cfg%y(j) - Hl).le.0.5_WP*dg) fs%V(i,j,k) = random_uniform(lo=-0.5_WP*a_vel,hi=+0.5_WP*a_vel)       
            end do
            do k=fs%cfg%kmin_,fs%cfg%kmax_+1
               do j=fs%cfg%jmin_,fs%cfg%jmax_+1
                  do i=fs%cfg%imin_,fs%cfg%imax_+1
                     if (fs%mask(i,j-1,k).eq.1.and.fs%mask(i-1,j-1,k).eq.1) fs%U(i,j-1,k)=0.0_WP
                     if (fs%mask(i,j-1,k).eq.1.and.fs%mask(i,j-1,k-1).eq.1) fs%W(i,j-1,k)=0.0_WP
                  end do 
               end do
            end do 
            call fs%cfg%sync(fs%U)
            call fs%cfg%sync(fs%V)
            ! Apply all other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            ! Adjust MFR for global mass balance
            call fs%correct_mfr()
            ! Compute cell-centered velocity
            call fs%interp_vel(Ui,Vi,Wi)
            ! Compute divergence
            call fs%get_div()
            ! Also update time
            call df%pull(name='t' ,val=time%t )
            call df%pull(name='dt',val=time%dt)
            time%told=time%t-time%dt
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call df%initialize(pg=cfg,iopartition=iopartition,filename=trim(cfg%name),nval=2,nvar=14)
            df%valname=['t ','dt']
            df%varname=['U   ','V   ','W   ','Umid','Vmid','Wmid','P   ','Pjx ','Pjy ','Pjz ','P11 ','P12 ','P13 ','P14 ']
         end if
      end block handle_restart
      

      create_iterator: block
         vof_removal_layer=iterator(cfg,'VOF removal',vof_removal_layer_locator)
         vof_removed=0.0_WP
      end block create_iterator


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
      
      function levelset_AWML(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G = -(xyz(2)-Hl)
      end function levelset_AWML
      
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
      
      !> Function that localizes the left (x-) of the domain
      function xm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin) isIn=.true.
      end function xm_locator

      function xp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function xp_locator
      
   end subroutine simulation_init
   
      !> Function that localizes region of VOF removal
   function vof_removal_layer_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.ge.pg%imax-nlayer) isIn=.true.
   end function vof_removal_layer_locator
   
   !> Perform an NGA2 simulation — KE Conservative staggered scheme
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc,bcond
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
         fs%Uold=fs%U; fs%sRHOXold=fs%sRHOX
         fs%Vold=fs%V; fs%sRHOYold=fs%sRHOY
         fs%Wold=fs%W; fs%sRHOZold=fs%sRHOZ
         

         update_temporal_variation: block
            use random,   only: random_uniform
            integer :: n,i,j,k
            real(WP) :: yval
            type(bcond), pointer :: mybc
            call fs%get_bcond('inlow',mybc)
            do n=1,mybc%itr%n_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               yval=cfg%ym(j)
               if (yval.le.Hl) then
                  fs%U(i,j,k)=(Ul*erf((Hl-yval)/dl) + Uint*(1.0_WP + erf((yval-Hl)/(dl*deficit))))*erf(yval/dl)
               else if (yval.gt.Hl.and.yval.le.Hl+Hg) then
                  fs%U(i,j,k)=(Ug*erf((yval-Hl)/dg) + Uint*(1.0_WP - erf((yval-Hl)/(dl*deficit))))* &
                  &           (erf((Hg+Hl-yval)/dg))+Usm*(1.0_WP +erf((yval-(Hg+Hl))/dg))
               else
                  fs%U(i,j,k)=Usm*(1-erf((yval-(Hl+Hg))/dg)) + Ucf*erf((yval-(Hl+Hg))/dg)
               end if
               if (abs(yval - Hl).le.0.5_WP*dg) fs%U(i,j,k) = fs%U(i,j,k)+random_uniform(lo=-0.5_WP*a_vel,hi=+0.5_WP*a_vel)
               if (abs(cfg%y(j) - Hl).le.0.5_WP*dg) fs%V(i,j,k) = random_uniform(lo=-0.5_WP*a_vel,hi=+0.5_WP*a_vel)       
            end do
            call fs%cfg%sync(fs%U)
            call fs%cfg%sync(fs%V)
         end block update_temporal_variation


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
         call fs%get_div()
         
          ! Compute stats
         ! call get_stats()
         ! call get_KE()
         ! Remove VOF at edge of domain
         remove_vof: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
            use parallel, only: MPI_REAL_WP
            integer :: n,i,j,k,ierr
            ! Adjust VOF field
            vof_removed=0.0_WP
            do n=1,vof_removal_layer%no_
               i=vof_removal_layer%map(1,n)
               j=vof_removal_layer%map(2,n)
               k=vof_removal_layer%map(3,n)
               if (n.le.vof_removal_layer%n_) vof_removed=vof_removed+cfg%vol(i,j,k)*vf%VF(i,j,k)
               vf%VF(i,j,k)=0.0_WP
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,vof_removed,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
            call vf%clean_irl_and_band()
            ! Also adjust density
            resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_density(rho=resU)
         end block remove_vof

         call processdata()

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
         
         if (save_evt%occurs()) then
            save_restart: block
               use irl_fortran_interface
               use string, only: str_medium
               character(len=str_medium) :: timestamp
               real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
               integer :: i,j,k
               real(WP), dimension(4) :: plane
               ! Handle IRL data
               allocate(P11(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(P12(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(P13(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(P14(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               do k=vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i=vf%cfg%imino_,vf%cfg%imaxo_
                        plane=getPlane(vf%liquid_gas_interface(i,j,k),0)
                        P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                     end do
                  end do
               end do
               ! Prefix for files
               write(timestamp,'(es12.5)') time%t
               ! Populate df and write it
               call df%push(name='t'   ,val=time%t )
               call df%push(name='dt'  ,val=time%dt)
               call df%push(name='U'   ,var=fs%U   )
               call df%push(name='V'   ,var=fs%V   )
               call df%push(name='W'   ,var=fs%W   )
               call df%push(name='Umid',var=fs%Umid)
               call df%push(name='Vmid',var=fs%Vmid)
               call df%push(name='Wmid',var=fs%Wmid)
               call df%push(name='P'   ,var=fs%P   )
               call df%push(name='Pjx' ,var=fs%Pjx )
               call df%push(name='Pjy' ,var=fs%Pjy )
               call df%push(name='Pjz' ,var=fs%Pjz )
               call df%push(name='P11' ,var=P11         )
               call df%push(name='P12' ,var=P12         )
               call df%push(name='P13' ,var=P13         )
               call df%push(name='P14' ,var=P14         )
               call df%write(fdata='restart/AWML_'//trim(adjustl(timestamp)))
               ! Deallocate
               deallocate(P11,P12,P13,P14)
            end block save_restart
         end if
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

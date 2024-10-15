!> Definition for a ligament atomization class
module ligament_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use iterator_class,    only: iterator
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use partmesh_class,    only: partmesh
   use hypre_str_class,   only: hypre_str
   !use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use lpt_class,         only: lpt
   use breakup_class,     only: breakup
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use pardata_class,     only: pardata
   use monitor_class,     only: monitor
   use string,            only: str_medium
   implicit none
   private
   
   public :: ligament
   
   !> Ligament object
   type :: ligament
      
      !> Provide a pardata and an event tracker for saving restarts
      type(event)    :: save_evt
      type(pardata)  :: df
      character(len=str_medium) :: lpt_file
      logical :: restarted

      !> Input file for the simulation
      type(inputfile) :: input

      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf    !< Volume fraction solver
      type(tpns)        :: fs    !< Two-phase flow solver
      type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
      !type(ddadi)       :: vs    !< DDADI solver for velocity
      type(timetracker) :: time  !< Time info
      
      !> Break-up modeling
      type(lpt)         :: lp    !< Lagrangian particle solver
      type(breakup)     :: bu    !< SGS break-up model

      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(partmesh) :: pmesh    !< For spray conversion
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
      
      
   contains
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
   end type ligament
   
   
   !> Hardcode size of buffer layer for VOF removal
   integer, parameter :: nlayer=5
   

contains
   
   
   !> Function that defines a level set function for a droplet
   function levelset_droplet(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
   end function levelset_droplet


   !> Function that defines a level set function for a ligament
   function levelset_ligament(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2)
   end function levelset_ligament
   
   
   !> Initialization of ligament simulation
   subroutine init(this)
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Setup an input file
      read_input: block
         use parallel, only: amRoot
         this%input=inputfile(amRoot=amRoot,filename='ligament.input')
      end block read_input
      
      ! Create the ligament mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xlig
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x(nx+1));  call this%input%read('X ligament',xlig)
                  call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y(ny+1))
                  call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xlig
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='Ligament')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         call this%input%read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      ! Handle restart/saves here
      restart_and_save: block
         use filesys, only: makedir,isdir
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         ! this%save_evt%tper = 5.0_WP
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call this%input%read('I/O partition',iopartition)
         ! Perform pardata initialization
         if (this%restarted) then
            ! We are restarting, read the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_'//trim(adjustl(timestamp)))
            ! this%lpt_file='restart/datalpt_'//trim(adjustl(timestamp))
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=15)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24']
         end if
      end block restart_and_save

      update_timetracker: block
         if (this%restarted) then
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
         end if
      end block update_timetracker

      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
       ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: VFlo,VFhi,elvira,r2p,flux,remap,r2plig,r2pnetlig
         use irl_fortran_interface
         use mms_geom,  only: cube_refine_vol
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         real(WP) :: vol,area,nx
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=r2plig,transport_method=flux,name='VOF')
         ! Initialize to a ligament
         if (this%restarted) then
            ! Read in the planes directly and set the IRL interface
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P21',var=P21)
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P22',var=P22)
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P23',var=P23)
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P24',var=P24)
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! Check if the second plane is meaningful
                     if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
                     else
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                     end if
                     ! For this restart, I want to attempt to remove all gas inclusions next to the walls
                     !rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                     !if (this%vf%cfg%xm(i).lt.-0.0015_WP.and.rad.le.0.002_WP.and.abs(this%cfg%Gib(i,j,k)).lt.2.0_WP*this%cfg%min_meshsize) then
                     !   call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                     !   call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],1.0_WP)
                     !else if (this%vf%cfg%xm(i).lt.0.0_WP.and.rad.le.0.00143_WP.and.abs(this%cfg%Gib(i,j,k)).lt.2.0_WP*this%cfg%min_meshsize) then
                     !   call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                     !   call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],1.0_WP)
                     !end if
                     ! For this restart, I want to remove all liquid outside the nozzle
                     !if (this%vf%cfg%xm(i).gt.0.0_WP) then
                     !   call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                     !   call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],-1.0_WP)
                     !end if
                     !rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                     !if (this%vf%cfg%xm(i).gt.-0.0015_WP.and.rad.gt.0.00143_WP) then
                     !   call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                     !   call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],-1.0_WP)
                     !end if
                  end do
               end do
            end do
            call this%vf%sync_interface()
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Ensure that boundaries are correct
            ! if (this%vf%cfg%iproc.eq.1)               this%vf%VF(this%vf%cfg%imino:this%vf%cfg%imin-1,:,:)=0.0_WP
            ! if (this%vf%cfg%iproc.eq.this%vf%cfg%npx) this%vf%VF(this%vf%cfg%imax+1:this%vf%cfg%imaxo,:,:)=0.0_WP
            ! if (this%vf%cfg%jproc.eq.1)               this%vf%VF(:,this%vf%cfg%jmino:this%vf%cfg%jmin-1,:)=0.0_WP
            ! if (this%vf%cfg%jproc.eq.this%vf%cfg%npy) this%vf%VF(:,this%vf%cfg%jmax+1:this%vf%cfg%jmaxo,:)=0.0_WP
            ! if (this%vf%cfg%kproc.eq.1)               this%vf%VF(:,:,this%vf%cfg%kmino:this%vf%cfg%kmin-1)=0.0_WP
            ! if (this%vf%cfg%kproc.eq.this%vf%cfg%npz) this%vf%VF(:,:,this%vf%cfg%kmax+1:this%vf%cfg%kmaxo)=0.0_WP
            ! do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            !    do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
            !       do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
            !          rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
            !          if (i.lt.this%vf%cfg%imin.and.rad.le.0.002_WP) this%vf%VF(i,j,k)=1.0_WP
            !       end do
            !    end do
            ! end do
            ! Update the band
            call this%vf%update_band()
            ! Set interface planes at the boundaries
            call this%vf%set_full_bcond()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
         else
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! Set cube vertices
                     n=0
                     do sk=0,1
                        do sj=0,1
                           do si=0,1
                              n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
                           end do
                        end do
                     end do
                     ! Call adaptive refinement code to get volume and barycenters recursively
                     vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_droplet,0.0_WP,amr_ref_lvl)
                     this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
                     if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                        this%vf%Lbary(:,i,j,k)=v_cent
                        this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
                     else
                        this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                        this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     end if
                  end do
               end do
            end do
            ! Update the band
            call this%vf%update_band()
            ! Perform interface reconstruction from VOF field
            call this%vf%build_interface()
            ! Set interface planes at the boundaries
            call this%vf%set_full_bcond()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Reset moments to guarantee compatibility with interface reconstruction
            call this%vf%reset_volume_moments()
         end if
      end block create_and_initialize_vof
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator

      
      ! Create a multiphase flow solver with bconds
      create_flow_solver: block
         use tpns_class,      only: dirichlet,clipped_neumann,bcond
         use hypre_str_class, only: pcg_pfmg2
         type(bcond), pointer :: mybc
         integer :: n,i,j,k      
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set fluid properties
         this%fs%rho_g=1.0_WP; call this%input%read('Density ratio',this%fs%rho_l)
         call this%input%read('Reynolds number',this%fs%visc_g); this%fs%visc_g=1.0_WP/this%fs%visc_g
         call this%input%read('Viscosity ratio',this%fs%visc_l); this%fs%visc_l=this%fs%visc_g*this%fs%visc_l
         call this%input%read('Weber number',this%fs%sigma); this%fs%sigma=1.0_WP/this%fs%sigma
         ! Define inflow boundary condition on the left
         call this%fs%add_bcond(name='inflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Define outflow boundary condition on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Configure implicit velocity solver
         !this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)!,implicit_solver=this%vs)
      end block create_flow_solver
      
       ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         ! Zero velocity except if restarting
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         if (this%restarted) then
            ! Read data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
            ! Apply boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
         else
            call this%fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               this%fs%U(i,j,k)=1.0_WP
            end do
         end if
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
      end block initialize_velocity

      ! Create lpt solver
      create_lpt_solver: block
         ! Create the solver
         this%lp=lpt(cfg=this%cfg,name='spray')
         ! Get particle density from the flow solver
         this%lp%rho=this%fs%rho_l
         ! Turn off drag
         this%lp%drag_model='none'         
         ! Initialize with zero particles
         call this%lp%resize(0)
         ! Get initial particle volume fraction
         call this%lp%update_VF()               
         ! Handle restarts
         ! if (this%restarted) call this%lp%read(filename=trim(this%lpt_file))
         call this%lp%get_max()
      end block create_lpt_solver
      
      
      ! Create breakup model
      create_breakup: block
        call this%bu%initialize(vf=this%vf,fs=this%fs,lp=this%lp)
      end block create_breakup
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         this%smesh=surfmesh(nvar=10,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='curv'
         this%smesh%varname(3)='edge_sensor'
         this%smesh%varname(4)='thin_sensor'
         this%smesh%varname(5)='thickness'
         this%smesh%varname(6)='struct_type'
         this%smesh%varname(7)='id_ccl_film'
         ! this%smesh%varname(8)='id_ccl_ligament'
         this%smesh%varname(8)='id_ccl'
         this%smesh%varname(9)='lig_ind'
         this%smesh%varname(10)='struct_type_ori'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh(this%smesh)
         ! Also populate nplane variable
         this%smesh%var(1,:)=1.0_WP
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                        this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                        this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                        this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                        this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                        this%smesh%var(6,np)=real(this%bu%struct_type(i,j,k),WP)
                        this%smesh%var(7,np)=real(this%bu%ccl_film%id(i,j,k),WP)
                        ! this%smesh%var(8,np)=real(this%vf%ccl%id(i,j,k),WP)
                        this%smesh%var(8,np)=real(this%bu%ccl%id(i,j,k),WP)
                        this%smesh%var(9,np)=real(this%vf%lig_ind(i,j,k),WP)
                        this%smesh%var(10,np)=real(this%vf%struct_type(i,j,k),WP) 
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh


       ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         ! Include an extra variable for droplet diameter
         this%pmesh=partmesh(nvar=2,nvec=1,name='lpt')
         this%pmesh%varname(1)='diameter'
         this%pmesh%varname(2)='id'
         this%pmesh%vecname(1)='velocity'
         ! Transfer particles to pmesh
         call this%lp%update_partmesh(this%pmesh)
         ! Also populate diameter variable
         do i=1,this%lp%np_
            this%pmesh%var(1,i)=this%lp%p(i)%d
            this%pmesh%var(2,i)=this%lp%p(i)%id
            this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
         end do
      end block create_pmesh

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='ligament')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_surface('plic',this%smesh)
         call this%ens_out%add_particle('spray',this%pmesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_atom')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFmax,'VOF maximum')
         call this%mfile%add_column(this%vf%VFmin,'VOF minimum')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%bu%vol_convert,'VOL converted')
         call this%mfile%add_column(this%vf%flotsam_error,'Flotsam error')
         call this%mfile%add_column(this%vf%thinstruct_error,'Film error')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_atom')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc,harmonic_visc
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Advance our spray
      this%resU=this%fs%rho_g; this%resV=this%fs%visc_g
      call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)

      ! Remember old VOF
      this%vf%VFold=this%vf%VF

      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
         
      ! VOF solver step
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Form implicit residuals
         !call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU/this%fs%rho_U
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV/this%fs%rho_V
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW/this%fs%rho_W
         
         ! Apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         !call this%fs%update_laplacian(pinpoint=[this%fs%cfg%imin,this%fs%cfg%jmin,this%fs%cfg%kmin])
         call this%fs%correct_mfr()
         call this%fs%get_div()
         !call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         ! call this%fs%add_surface_tension_jump_thin(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         call this%fs%add_surface_tension_jump_twoVF(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         !if (this%cfg%amRoot) this%fs%psolv%rhs(this%cfg%imin,this%cfg%jmin,this%cfg%kmin)=0.0_WP
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         ! Apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      
      ! Breakup
      call this%bu%attempt_breakup()

      ! Remove VOF at edge of domain
      remove_vof: block
         integer :: n
         do n=1,this%vof_removal_layer%no_
            this%vf%VF(this%vof_removal_layer%map(1,n),this%vof_removal_layer%map(2,n),this%vof_removal_layer%map(3,n))=0.0_WP
         end do
      end block remove_vof
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update surfmesh object
         update_smesh: block
            use irl_fortran_interface
            integer :: i,j,k,nplane,np
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh(this%smesh)
            ! Also populate nplane variable
            this%smesh%var(1,:)=1.0_WP
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                           this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                           this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                           this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                           this%smesh%var(6,np)=real(this%bu%struct_type(i,j,k),WP)
                           this%smesh%var(7,np)=real(this%bu%ccl_film%id(i,j,k),WP)
                           ! this%smesh%var(8,np)=real(this%vf%ccl%id(i,j,k),WP)
                           this%smesh%var(8,np)=real(this%bu%ccl%id(i,j,k),WP)
                           this%smesh%var(9,np)=real(this%vf%lig_ind(i,j,k),WP)
                           this%smesh%var(10,np)=real(this%vf%struct_type(i,j,k),WP) 
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Update partmesh object
         update_pmesh: block
            integer :: i
            ! Transfer particles to pmesh
            call this%lp%update_partmesh(this%pmesh)
            ! Also populate diameter variable
            do i=1,this%lp%np_
               this%pmesh%var(1,i)=this%lp%p(i)%d
               this%pmesh%var(2,i)=this%lp%p(i)%id
               this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
            end do
         end block update_pmesh         
         ! Perform ensight output
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            character(len=str_medium) :: timestamp
            real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
            real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
            integer :: i,j,k
            real(WP), dimension(4) :: plane
            ! Handle IRL data
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! First plane
                     plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                     P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                     ! Second plane
                     plane=0.0_WP
                     if (getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.2) plane=getPlane(this%vf%liquid_gas_interface(i,j,k),1)
                     P21(i,j,k)=plane(1); P22(i,j,k)=plane(2); P23(i,j,k)=plane(3); P24(i,j,k)=plane(4)
                  end do
               end do
            end do
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t'  ,val=this%time%t )
            call this%df%push(name='dt' ,val=this%time%dt)
            call this%df%push(name='U'  ,var=this%fs%U   )
            call this%df%push(name='V'  ,var=this%fs%V   )
            call this%df%push(name='W'  ,var=this%fs%W   )
            call this%df%push(name='P'  ,var=this%fs%P   )
            call this%df%push(name='Pjx',var=this%fs%Pjx )
            call this%df%push(name='Pjy',var=this%fs%Pjy )
            call this%df%push(name='Pjz',var=this%fs%Pjz )
            call this%df%push(name='P11',var=P11         )
            call this%df%push(name='P12',var=P12         )
            call this%df%push(name='P13',var=P13         )
            call this%df%push(name='P14',var=P14         )
            call this%df%push(name='P21',var=P21         )
            call this%df%push(name='P22',var=P22         )
            call this%df%push(name='P23',var=P23         )
            call this%df%push(name='P24',var=P24         )
            call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
            call this%lp%write(filename='restart/datalpt_'//trim(adjustl(timestamp)))
            ! Deallocate
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
         end block save_restart
      end if
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      
   end subroutine final
   
   
   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   
   
   !> Function that localizes region of VOF removal
   function vof_removal_layer_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.ge.pg%imax-nlayer) isIn=.true.
   end function vof_removal_layer_locator
   
   
end module ligament_class
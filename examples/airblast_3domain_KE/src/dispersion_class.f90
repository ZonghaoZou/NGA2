!> Definition for a dispersion class
module dispersion_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use pardata_class,     only: pardata
   use monitor_class,     only: monitor
   use lpt_class,         only: lpt
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use partmesh_class,    only: partmesh
   implicit none
   private
  !  type(lpt), public :: lp         !< Lagrangian particle tracking    
   public :: dispersion
   
   !> Nozzle object
   type :: dispersion
      
      !> Provide a datafile and an event tracker for saving restarts
      type(event)    :: save_evt
      type(pardata)  :: df
      logical :: restarted
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB
      type(ibconfig) :: cfg
      
      !> Surface mesh for IB
      type(surfmesh) :: plymesh
      
      !> Flow solver
      type(incomp)      :: fs    !< Incompressible flow solver
      type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
      type(sgsmodel)    :: sgs   !< SGS model for eddy viscosity
      type(timetracker) :: time  !< Time info

      type(lpt)      :: lp 
      type(partmesh) :: pmesh      !< Particle mesh for lpt

      !> Ensight postprocessing
      type(ensight) :: ens_out  !< Ensight output for flow variables
      type(event)   :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> IB velocity and mass source
      real(WP), dimension(:,:,:), allocatable :: Uib,Vib,Wib
      real(WP), dimension(:,:,:), allocatable :: U2on3,V2on3,W2on3
      
      !> Fluid definition
      real(WP) :: visc
      
   contains
      procedure, private :: geometry_init          !< Initialize geometry for dispersion
      procedure, private :: simulation_init        !< Initialize simulation for dispersion
      procedure :: init                            !< Initialize dispersion simulation
      procedure :: step                            !< Advance dispersion simulation by one time step
      procedure :: final                           !< Finalize dispersion simulation
      procedure :: record_droplet
   end type dispersion
   

   !> Hardcode inlet positions used in locator functions at x=-0.01
  ! real(WP), parameter, public :: dl=0.0025_WP   ! Liquid pipe diameter ~(inner+outer)/2
  real(WP), parameter, public :: dl=0.003_WP   ! Liquid outer pipe diameter 
  ! real(WP), parameter, public :: dg=0.0100_WP   ! Gas pipe diameter ~(inner+outer)/2
  ! 0.0206 doesn't seem right, it seems to be over estimating
  real(WP), parameter, public :: dg=0.0206_WP   ! Gas pipe diameter ~(inner+outer)/2
  real(WP), parameter, public :: rl=0.0010_WP   ! Liquid pipe inner radius
  real(WP), parameter, public :: rlo=0.0015_WP   ! Liquid pipe outer radius
  real(WP), parameter, public :: rgi=0.005_WP   ! Liquid pipe outer radius
  real(WP) :: rho_l,visc_l
   
contains
   
 !> Initialization of dispersion simulation
    subroutine record_droplet(this)
     use parallel,  only: MPI_REAL_WP
     use string,   only: str_medium
     use messager, only: die
     use mpi_f08
     implicit none
     class(dispersion), intent(inout) :: this
     character(len=str_medium) :: filename
     real(WP), dimension(:,:), allocatable :: pinfo,pinfo_
     integer, dimension(:), allocatable:: plist,dispels
     real(WP) :: xloc_90,xloc_100,xloc_150,xloc_200,input_xloc
     integer:: n,count_90,count_100,count_150,count_200,totalcount,input_count,i
     integer:: rank,count,ierr,iunit
     xloc_90=90e-3_WP;xloc_100=100e-3_WP;xloc_150=150e-3_WP;xloc_200=200e-3_WP
     count_90=0;count_100=0;count_150=0;count_200=0 
     allocate(plist(0:this%cfg%nproc-1))
     ! For each particle on each processor, count how many have passed the different x locations
     do n =1, this%lp%np_
       if (this%lp%p(n)%pos(1).lt.xloc_90 .and. this%lp%p(n)%pos(1)+this%time%dt*this%lp%p(n)%vel(1).ge.xloc_90) count_90=count_90+1
       if (this%lp%p(n)%pos(1).lt.xloc_100 .and. this%lp%p(n)%pos(1)+this%time%dt*this%lp%p(n)%vel(1).ge.xloc_100) count_100=count_100+1
       if (this%lp%p(n)%pos(1).lt.xloc_150 .and. this%lp%p(n)%pos(1)+this%time%dt*this%lp%p(n)%vel(1).ge.xloc_150) count_150=count_150+1
       if (this%lp%p(n)%pos(1).lt.xloc_200 .and. this%lp%p(n)%pos(1)+this%time%dt*this%lp%p(n)%vel(1).ge.xloc_200) count_200=count_200+1
     end do
     
     input_count=count_90; input_xloc=xloc_90; call output()
     input_count=count_100; input_xloc=xloc_100; call output()
     input_count=count_150; input_xloc=xloc_150; call output()
     input_count=count_200; input_xloc=xloc_200; call output()

     contains
     
     subroutine output()
        implicit none 
        ! Lets first deal with xloc = 30e-3
        call MPI_AllGATHER(input_count,1,MPI_INTEGER,plist,1,MPI_INTEGER,this%cfg%comm,ierr)
        totalcount=sum(plist)
        if (totalcount .gt. 0) then
           allocate(pinfo_(1:8,1:input_count))
           allocate(pinfo(1:8,1:totalcount))
           allocate(dispels(0:this%cfg%nproc-1))
           input_count=0
           do n =1, this%lp%np_
              if (this%lp%p(n)%pos(1).lt.input_xloc .and. this%lp%p(n)%pos(1)+this%time%dt*this%lp%p(n)%vel(1).ge.input_xloc) then
                 input_count=input_count+1
                 pinfo_(1,  input_count)=this%lp%p(n)%d
                 pinfo_(2:4,input_count)=this%lp%p(n)%vel
                 pinfo_(5:7,input_count)=this%lp%p(n)%pos
                 pinfo_(8,  input_count)=this%lp%p(n)%id
              end if
           end do
           ! Calculate dispels
           count = 0
           do rank=0,this%cfg%nproc-1
              dispels(rank) = count
              count = count + plist(rank)
           end do
           ! Communicate to root
           do i = 1,8
              call MPI_GATHERV(pinfo_(i,:),input_count,MPI_REAL_WP,pinfo(i,:),plist,dispels,MPI_REAL_WP,0,this%cfg%comm)
           end do
           !!! Write to droplet list !!!
           if (this%cfg%amRoot)  then
            !   filename='spray-disper/x=30e-3'
              write(filename, '("spray-disper/x=",ES10.3)') input_xloc
              open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
              if (ierr.ne.0) call die('[Dipersion stat analysis] Could not open file: '//trim(filename))
              do i = 1,totalcount
              write(iunit,'(f24.16,1x,f24.16,1x,f24.16,1x,f24.16,1x,f24.16,f24.16,1x,f24.16,1x,f24.16,1x,I2)')this%time%t,pinfo(1,i),pinfo(2,i),pinfo(3,i)&
              &,pinfo(4,i),pinfo(5,i),pinfo(6,i),pinfo(7,i),INT(pinfo(8,i))
              end do
              close(iunit)
           end if
           deallocate(pinfo,pinfo_,dispels)
        end if 
     end subroutine
     
   end subroutine record_droplet
   
   !> Initialization of dispersion simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(dispersion), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_dispersion')
      
      ! Initialize the geometry
      call this%geometry_init()
      
      ! Initialize the simulation
      call this%simulation_init()
      
   end subroutine init
   
   
   !> Initialize geometry
   subroutine geometry_init(this)
      use sgrid_class, only: sgrid
      implicit none
      class(dispersion) :: this
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xshift
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x(nx+1)); call this%input%read('X shift',xshift)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='dispersion')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call this%input%read('Partition',partition)

         ! Create partitioned grid
         this%cfg=ibconfig(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      
      ! Read in the PLY geometry
      read_ply: block
         use string,   only: str_medium
         use parallel, only: MPI_REAL_WP
         use mpi_f08
         character(len=str_medium) :: plyfile
         integer :: ierr,size_conn
         
         ! Read in ply filename
         call this%input%read('PLY filename',plyfile)
         
         ! Root creates surface mesh from ply, other process create empty surface mesh
         if (this%cfg%amRoot) then
            this%plymesh=surfmesh(plyfile=plyfile,nvar=0,name='ply')
         else
            this%plymesh=surfmesh(nvar=0,name='ply')
         end if
         
         ! Go through parallel broadcast of surface mesh
         call MPI_BCAST(this%plymesh%nVert   ,1                 ,MPI_INTEGER,0,this%cfg%comm,ierr)
         call MPI_BCAST(this%plymesh%nPoly   ,1                 ,MPI_INTEGER,0,this%cfg%comm,ierr)
         if (.not.this%cfg%amRoot) call this%plymesh%set_size(nvert=this%plymesh%nVert,npoly=this%plymesh%nPoly)
         call MPI_BCAST(this%plymesh%xVert   ,this%plymesh%nVert,MPI_REAL_WP,0,this%cfg%comm,ierr)
         call MPI_BCAST(this%plymesh%yVert   ,this%plymesh%nVert,MPI_REAL_WP,0,this%cfg%comm,ierr)
         call MPI_BCAST(this%plymesh%zVert   ,this%plymesh%nVert,MPI_REAL_WP,0,this%cfg%comm,ierr)
         call MPI_BCAST(this%plymesh%polySize,this%plymesh%nPoly,MPI_INTEGER,0,this%cfg%comm,ierr)
         if (this%cfg%amRoot) size_conn=size(this%plymesh%polyConn)
         call MPI_BCAST(size_conn            ,1                 ,MPI_INTEGER,0,this%cfg%comm,ierr)
         if (.not.this%cfg%amRoot) allocate(this%plymesh%polyConn(size_conn))
         call MPI_BCAST(this%plymesh%polyConn,size_conn         ,MPI_INTEGER,0,this%cfg%comm,ierr)
         
      end block read_ply
      
      
      ! Create IB walls for this config
      create_walls: block
         use ibconfig_class, only: sharp
         use messager,       only: die
         use mathtools,      only: cross_product,normalize
         use irl_fortran_interface
         real(WP), parameter :: safe_coeff=3.0_WP
         integer :: i,j,k,np,nv,iv,ip
         real(WP) :: mydist
         real(WP), dimension(3) :: pos,nearest_pt,mynearest,mynorm
         type(Poly_type), dimension(:), allocatable :: poly
         real(WP), dimension(:,:), allocatable :: bary,vert,nvec
         
         ! Preprocess surface mesh data using IRL
         allocate(poly(1:this%plymesh%nPoly))
         allocate(vert(1:3,1:maxval(this%plymesh%polySize)))
         allocate(bary(1:3,1:this%plymesh%nPoly))
         allocate(nvec(1:3,1:this%plymesh%nPoly))
         do np=1,this%plymesh%nPoly
            ! Allocate polygon
            call new(poly(np))
            ! Fill it up
            do nv=1,this%plymesh%polySize(np)
               iv=sum(this%plymesh%polySize(1:np-1))+nv
               vert(:,nv)=[this%plymesh%xVert(this%plymesh%polyConn(iv)),this%plymesh%yVert(this%plymesh%polyConn(iv)),this%plymesh%zVert(this%plymesh%polyConn(iv))]
            end do
            call construct(poly(np),this%plymesh%polySize(np),vert(1:3,1:this%plymesh%polySize(np)))
            mynorm=normalize(cross_product(vert(:,2)-vert(:,1),vert(:,3)-vert(:,2)))
            mydist=dot_product(mynorm,vert(:,1))
            call setPlaneOfExistence(poly(np),[mynorm(1),mynorm(2),mynorm(3),mydist])
            ! Also store its barycenter
            bary(:,np)=calculateCentroid(poly(np))
            nvec(:,np)=calculateNormal  (poly(np))
         end do
         deallocate(vert)
         
         ! Create IB distance field using IRL
         this%cfg%Gib=huge(1.0_WP)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Store cell center position
                  pos=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  ! Traverse all polygons
                  do np=1,this%plymesh%nPoly
                     ! Calculate distance to centroid
                     nearest_pt=pos-bary(:,np)
                     mydist=dot_product(nearest_pt,nearest_pt)
                     ! If close enough, compute exact distance to the polygon instead
                     if (mydist.lt.(safe_coeff*this%cfg%min_meshsize)**2) then
                        nearest_pt=calculateNearestPtOnSurface(poly(np),pos)
                        nearest_pt=pos-nearest_pt
                        mydist=dot_product(nearest_pt,nearest_pt)
                     end if
                     ! Remember closest distance
                     if (mydist.lt.this%cfg%Gib(i,j,k)) then
                        this%cfg%Gib(i,j,k)=mydist
                        mynearest=nearest_pt
                        ip=np
                     end if
                  end do
                  ! Take the square root
                  this%cfg%Gib(i,j,k)=sqrt(this%cfg%Gib(i,j,k))
                  ! Find the sign
                  if (dot_product(mynearest,nvec(:,ip)).gt.0.0_WP) this%cfg%Gib(i,j,k)=-this%cfg%Gib(i,j,k)
               end do
            end do
         end do
         deallocate(bary,nvec,poly)
         
         ! Get normal vector
         call this%cfg%calculate_normal()
         
         ! Get VF field
         call this%cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
        
      end block create_walls
      
      
   end subroutine geometry_init
   
   
   !> Initialize simulation
   subroutine simulation_init(this)
      implicit none
      class(dispersion), intent(inout) :: this
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Uib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))  
         allocate(this%U2on3  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%V2on3  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%W2on3  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%x_cell(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%x_cell=0.0_WP
         allocate(this%y_cell(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%y_cell=0.0_WP
         allocate(this%z_cell(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%z_cell=0.0_WP
         allocate(this%x_face(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%x_face=0.0_WP
         allocate(this%y_face(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%y_face=0.0_WP
         allocate(this%z_face(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%z_face=0.0_WP
      end block allocate_work_arrays
      
      ! Handle restart/saves here
      restart_and_save: block
         use string, only: str_medium
         use filesys,  only: makedir,isdir
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call this%input%read('I/O partition',iopartition)
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_dispersion_'//trim(timestamp))
         else
            ! Prepare a new directory for storing files for restart
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=6)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','LM ','MM ']
         end if
      end block restart_and_save
      
      
      ! Revisit timetracker to adjust time and time step values if this is a restart
      update_timetracker: block
         if (this%restarted) then
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
         end if
      end block update_timetracker

      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         use incomp_class, only: clipped_neumann,dirichlet
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         call this%input%read('Gas density',this%fs%rho)
         call this%input%read('Gas dynamic viscosity',this%visc); this%fs%visc=this%visc
         call this%input%read('Liquid density',rho_l)
         call this%input%read('Liquid dynamic viscosity',visc_l)
         ! Define gas and liquid inlet boundary conditions
        call this%fs%add_bcond(name='gas_inlet',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=gas_inlet)
        ! call this%fs%add_bcond(name='liq_inlet',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=liq_inlet)
        ! Outflow on the right
        call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=right_boundary)
        ! Slip on the sides
        call this%fs%add_bcond(name='bc_yp',type=clipped_neumann,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
        call this%fs%add_bcond(name='bc_ym',type=clipped_neumann,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
        call this%fs%add_bcond(name='bc_zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
        call this%fs%add_bcond(name='bc_zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
        this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg,nst=7)
        this%ps%maxlevel=16
        call this%input%read('Pressure iteration',this%ps%maxit)
        call this%input%read('Pressure tolerance',this%ps%rcvg)
        ! Configure implicit velocity solver
        ! this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
        ! Setup the solver
        call this%fs%setup(pressure_solver=this%ps)!,implicit_solver=this%vs)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
     initialize_velocity: block
        use mpi_f08,    only: MPI_ALLREDUCE,MPI_SUM
        use parallel,   only: MPI_REAL_WP
        use incomp_class, only: bcond
        type(bcond), pointer :: mybc
        integer  :: n,i,j,k,ierr
        real(WP) :: Ugas,myAgas,Agas,Uliq,myAliq,Aliq
        real(WP) :: Qgas,Qliq,myU
        real(WP), parameter :: SLPM2SI=1.66667E-5_WP
        ! Zero initial field
        this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
        if (this%restarted) then
           ! Read data
           call this%df%pull(name='U',var=this%fs%U)
           call this%df%pull(name='V',var=this%fs%V)
           call this%df%pull(name='W',var=this%fs%W)
           call this%df%pull(name='P',var=this%fs%P)  !< Reset pressure upon restart because I've noticed IB is causing drift...
        end if
        ! Read in gas flow rate and convert to SI
        call this%input%read('Gas flow rate (SLPM)',Qgas)
        Qgas=Qgas*SLPM2SI
        ! Calculate gas flow area - no overlap here!
        myAgas=0.0_WP
        call this%fs%get_bcond('gas_inlet',mybc)
        do n=1,mybc%itr%n_
           i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
           myAgas=myAgas+this%cfg%dy(j)*this%cfg%dz(k)*sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
        end do
        call MPI_ALLREDUCE(myAgas,Agas,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
        ! Calculate bulk gas velocity
        Ugas=Qgas/Agas
        ! Apply Dirichlet at gas inlet
        call this%fs%get_bcond('gas_inlet',mybc)
        do n=1,mybc%itr%no_
           i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
           this%fs%U(i,j,k)=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*Ugas
        end do
        ! Apply all other boundary conditions
        call this%fs%apply_bcond(this%time%t,this%time%dt)
        ! Compute MFR through all boundary conditions
        call this%fs%get_mfr()
        ! Adjust MFR for global mass balance
        call this%fs%correct_mfr()
        ! Compute cell-centered velocity
        call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
        ! Compute divergence
        call this%fs%get_div()
     end block initialize_velocity
      
      
     initialize_lpt: block
        use string, only: str_medium
        use filesys,  only: makedir,isdir
        use messager, only: die
        character(len=str_medium) :: timestamp,filename
        logical :: partfile_exists
        integer :: ierr,iunit
        real(WP) :: input_xloc
        this%lp=lpt(cfg=this%cfg,name='spray_dispersion')
        this%lp%rho=rho_l
        call this%lp%resize(0)
        ! this%lp%filter_width=3.5_WP*this%cfg%min_meshsize
        if (this%restarted) then
           call this%input%read('Restart from',timestamp,default='')
           inquire(file='restart/part_dispersion_'//trim(timestamp),exist=partfile_exists)
           ! If so, read it
           if (partfile_exists) call this%lp%read(filename='restart/part_dispersion_'//trim(timestamp))
        end if

        if (this%lp%cfg%amroot) then
           if (.not.isdir('spray-disper')) call makedir('spray-disper')
           input_xloc = 90e-3_WP; write(filename, '("spray-disper/x=",ES10.3)') input_xloc
           open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',access='stream',iostat=ierr)
           if (ierr.ne.0) call die('[Dipersion stat analysis] Could not open file: '//trim(filename))
           close(iunit)

           input_xloc = 100e-3_WP; write(filename, '("spray-disper/x=",ES10.3)') input_xloc
           open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',access='stream',iostat=ierr)
           if (ierr.ne.0) call die('[Dipersion stat analysis] Could not open file: '//trim(filename))
           close(iunit)

           input_xloc = 150e-3_WP; write(filename, '("spray-disper/x=",ES10.3)') input_xloc
           open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',access='stream',iostat=ierr)
           if (ierr.ne.0) call die('[Dipersion stat analysis] Could not open file: '//trim(filename))
           close(iunit)     
           
           input_xloc = 200e-3_WP; write(filename, '("spray-disper/x=",ES10.3)') input_xloc
           open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',access='stream',iostat=ierr)
           if (ierr.ne.0) call die('[Dipersion stat analysis] Could not open file: '//trim(filename))
           close(iunit)     
        end if
     end block initialize_lpt

     
     create_pmesh: block
        integer :: i
        this%pmesh=partmesh(nvar=2,nvec=1,name='lpt')
        this%pmesh%varname(1)='radius'
        this%pmesh%varname(2)='id'
        this%pmesh%vecname(1)='velocity'
        call this%lp%update_partmesh(this%pmesh)
        do i=1,this%lp%np_
           this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
           this%pmesh%var(2,i)=this%lp%p(i)%id
           this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
        end do
     end block create_pmesh

      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
         if (this%restarted) then
            call this%df%pull(name='LM',var=this%sgs%LM)
            call this%df%pull(name='MM',var=this%sgs%MM)
         end if
      end block create_sgs
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='dispersion')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_particle('part',this%pmesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_dispersion')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_dispersion')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   

   !> Take one time step
   subroutine step(this,cfga2d)
      implicit none
      class(dispersion), intent(inout) :: this
      type(ibconfig), intent(inout) :: cfga2d
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
     
      call this%record_droplet() 
      
      this%resU=this%fs%rho
      this%resV=this%fs%visc
      call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)

      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Turbulence modeling
      sgs_modeling: block
         use sgsmodel_class, only: vreman
         this%resU=this%fs%rho
         call this%fs%get_gradu(this%gradU)
         call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
         this%fs%visc=this%visc+this%sgs%visc
      end block sgs_modeling
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*(this%fs%rho*this%fs%U-this%fs%rho*this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%rho*this%fs%V-this%fs%rho*this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%rho*this%fs%W-this%fs%rho*this%fs%Wold)+this%time%dt*this%resW
         

         nudge: block
           integer :: i,j,k
           real(WP) :: xcoord,ycoord,zcoord
           do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
              do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                 do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                    if (this%fs%umask(i,j,k).eq.0) then
                       xcoord=max((cfga2d%x(cfga2d%imax+1)-    this%fs%cfg%x (i) )/(cfga2d%x(cfga2d%imax+1)),0.0_WP)
                       ycoord=max((        0.5_WP*cfga2d%yL-abs(this%fs%cfg%ym(j)))/(0.5_WP*cfga2d%yL        ),0.0_WP)
                       zcoord=max((        0.5_WP*cfga2d%zL-abs(this%fs%cfg%zm(k)))/(0.5_WP*cfga2d%zL        ),0.0_WP)
                       this%resU(i,j,k)=this%resU(i,j,k)+(this%U2on3(i,j,k)-this%fs%U(i,j,k))*(xcoord*ycoord*zcoord)**2
                    end if
                    if (this%fs%vmask(i,j,k).eq.0) then
                       xcoord=max((cfga2d%x(cfga2d%imax+1)-    this%fs%cfg%xm(i) )/(cfga2d%x(cfga2d%imax+1)),0.0_WP)
                       ycoord=max((        0.5_WP*cfga2d%yL-abs(this%fs%cfg%y (j)))/(0.5_WP*cfga2d%yL        ),0.0_WP)
                       zcoord=max((        0.5_WP*cfga2d%zL-abs(this%fs%cfg%zm(k)))/(0.5_WP*cfga2d%zL        ),0.0_WP)
                       this%resV(i,j,k)=this%resV(i,j,k)+(this%V2on3(i,j,k)-this%fs%V(i,j,k))*(xcoord*ycoord*zcoord)**2
                    end if
                    if (this%fs%wmask(i,j,k).eq.0) then
                       xcoord=max((cfga2d%x(cfga2d%imax+1)-    this%fs%cfg%xm(i) )/(cfga2d%x(cfga2d%imax+1)),0.0_WP)
                       ycoord=max((        0.5_WP*cfga2d%yL-abs(this%fs%cfg%ym(j)))/(0.5_WP*cfga2d%yL        ),0.0_WP)
                       zcoord=max((        0.5_WP*cfga2d%zL-abs(this%fs%cfg%z (k)))/(0.5_WP*cfga2d%zL        ),0.0_WP)
                       this%resW(i,j,k)=this%resW(i,j,k)+(this%W2on3(i,j,k)-this%fs%W(i,j,k))*(xcoord*ycoord*zcoord)**2
                    end if
                 end do
              end do
           end do
        end block nudge
         ! Form implicit residuals
        !  call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU/this%fs%rho
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV/this%fs%rho
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW/this%fs%rho
         
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Solve Poisson equation
         call this%fs%correct_mfr()
         call this%fs%get_div()  !< a volume source term to div
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div*this%fs%rho/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()  !< a volume source term to div
      
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then 
           update_pmesh: block
              integer :: i
              call this%lp%update_partmesh(this%pmesh)
              do i=1,this%lp%np_
                 this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
                 this%pmesh%var(2,i)=this%lp%p(i)%id
                 this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
              end do
           end block update_pmesh 
           call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t' ,val=this%time%t )
            call this%df%push(name='dt',val=this%time%dt)
            call this%df%push(name='U' ,var=this%fs%U   )
            call this%df%push(name='V' ,var=this%fs%V   )
            call this%df%push(name='W' ,var=this%fs%W   )
            call this%df%push(name='P' ,var=this%fs%P   )
            call this%df%push(name='LM', var=this%sgs%LM)
            call this%df%push(name='MM', var=this%sgs%MM)
            call this%df%write(fdata='restart/data_dispersion_'//trim(adjustl(timestamp)))
            call this%lp%write(filename='restart/part_dispersion_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
      
   end subroutine step
   

   !> Finalize dispersion simulation
   subroutine final(this)
      implicit none
      class(dispersion), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      deallocate(this%Uib,this%Vib,this%Wib)
      
   end subroutine final
   
   
   !> Function that localizes the right domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imax+1) isIn=.true.
  end function right_boundary
  

  !> Function that localizes liquid stream at -x
  function liq_inlet(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     real(WP) :: rad
     isIn=.false.
     rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
     if (rad.lt.0.5_WP*dl.and.i.eq.pg%imin) isIn=.true.
  end function liq_inlet
  
  
  !> Function that localizes gas stream at -x
  function gas_inlet(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     real(WP) :: rad
     isIn=.false.
     rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
     if (rad.ge.0.5_WP*dl.and.rad.lt.0.5_WP*dg.and.i.eq.pg%imin) isIn=.true.
  end function gas_inlet
   
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
  
  
  !> Function that localizes the top (z+) of the domain
  function zp_locator(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (k.eq.pg%kmax+1) isIn=.true.
  end function zp_locator
  
  
  !> Function that localizes the bottom (z-) of the domain
  function zm_locator(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (k.eq.pg%kmin) isIn=.true.
  end function zm_locator

   
end module dispersion_class
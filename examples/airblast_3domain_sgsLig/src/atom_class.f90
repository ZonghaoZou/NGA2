!> Definition for an atomization class
module atom_class
   use precision,         only: WP
   use string,            only: str_medium
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use iterator_class,    only: iterator
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use pardata_class,     only: pardata
   use monitor_class,     only: monitor
   use timer_class,       only: timer
   use lpt_class,         only: lpt
   use breakup_class,     only: breakup
   implicit none
   private
   
   public :: atom
   
   !> Atom object
   type :: atom
      
      !> Provide a datafile and an event tracker for saving restarts
      type(event)    :: save_evt
      type(pardata) :: df
      logical :: restarted

      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB
      type(ibconfig) :: cfg
      
      !> Surface mesh for IB
      type(surfmesh) :: plymesh
      
      !> Flow solver
      type(vfs)         :: vf    !< Volume fraction solver
      type(tpns)        :: fs    !< Two-phase flow solver
      type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
      type(breakup)     :: bu    !< Breakup model
      type(ddadi)       :: vs    !< DDADI solver for velocity
      type(sgsmodel)    :: sgs   !< SGS model for eddy viscosity
      type(timetracker) :: time  !< Time info

      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      ! Stats recording
      type(event)    :: stat_evt
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Timing info
      type(monitor) :: timefile !< Timing monitoring
      type(timer)   :: tstep    !< Timer for step
      type(timer)   :: tvel     !< Timer for velocity
      type(timer)   :: tpres    !< Timer for pressure
      type(timer)   :: tvof     !< Timer for VOF
      type(timer)   :: tdrop   !< Timer for transfer
      type(timer)   :: tfilm   !< Timer for transfer
      type(timer)   :: tlig   !< Timer for transfer

      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities

      !> Record stats
      real(WP), dimension(:), allocatable :: EPLmean,EPL2mean,Umean1em3,U2mean1em3
      real(WP) :: timeint,stat_begin

      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
      integer  :: nlayer=4                 !< Size of buffer layer for VOF removal

      logical :: restart_from_sgslig,using_sgslig
   contains
      procedure, private :: geometry_init          !< Initialize geometry for nozzle
      procedure, private :: simulation_init        !< Initialize simulation for nozzle
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
      procedure :: record_stats
   end type atom

   
   !> Hardcode inlet positions used in locator functions at x=-0.01
   real(WP), parameter, public :: dl=0.003_WP   ! Liquid outer pipe diameter 
   real(WP), parameter, public :: dg=0.0206_WP   ! Gas pipe diameter ~(inner+outer)/2
   real(WP), parameter, public :: rl=0.0010_WP   ! Liquid pipe inner radius
   real(WP), parameter, public :: rlo=0.0015_WP   ! Liquid pipe outer radius
   
contains
   !> Initialization of dispersion simulation
   subroutine record_stats(this)
      use parallel,  only: MPI_REAL_WP
      use string,   only: str_medium
      use messager, only: die
      use mpi_f08
      implicit none
      class(atom), intent(inout) :: this
      character(len=str_medium) :: filename
      integer:: ierr,t_output,i,j,k
      real(WP) :: xloc
      real(WP), dimension(:), allocatable :: VFint,Uint
      if (this%time%t.lt.this%stat_begin) return
      ! Array extends to the entire domain
      allocate(VFint(this%cfg%imin:this%cfg%imax));VFint=0.0_WP
      allocate(Uint(this%cfg%jmin:this%cfg%jmax));Uint=0.0_WP
      ! Accumulate the data
      xloc=1.0e-3
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%cfg%ym(j).ge.0.0_WP.and.this%cfg%ym(j-1).lt.0.0_WP) then
                  VFint(i)=VFint(i)+this%vf%VF(i,j,k)*this%cfg%dz(k)
               end if
               if (this%cfg%zm(k).ge.0.0_WP.and.this%cfg%zm(k-1).lt.0.0_WP)then
                  if (this%cfg%xm(i).ge.xloc.and.this%cfg%xm(i-1).lt.xloc) then
                     Uint(j)=this%fs%U(i,j,k)
                  end if 
               end if
            end do
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFint,size(VFint),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Uint,size(Uint),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      this%EPLmean =this%EPLmean +this%time%dt*VFint
      this%EPL2mean=this%EPL2mean+this%time%dt*VFint**2
      this%Umean1em3 =this%Umean1em3 +this%time%dt*Uint
      this%U2mean1em3=this%U2mean1em3+this%time%dt*Uint**2
      this%timeint=this%timeint+this%time%dt
      deallocate(VFint,Uint)
      ! If it is time to record stats
      if (this%stat_evt%occurs().and.this%timeint.gt.0.0_WP) then
         if (this%cfg%amRoot) then
            t_output=int(this%time%t*1.0e4)

            write(filename, '(A,I0,A)') 'stats/EPL_', t_output, '.csv'
            open(unit=10, file=filename, status="replace", action="write")
            do i=this%cfg%imin,this%cfg%imax
               write(10, '(F20.12, ",", F20.12, ",", F20.12, ",", F20.12)') this%cfg%xm(i),this%EPLmean(i)/this%timeint,&
               & sqrt(max(0.0_WP,this%EPL2mean(i)/this%timeint-(this%EPLmean(i)/this%timeint)**2)),this%EPL2mean(i)/this%timeint
            end do
            close(unit=10)

            write(filename, '(A,I0,A)') 'stats/U1em3_', t_output, '.csv'
            open(unit=10, file=filename, status="replace", action="write")
            do j=this%cfg%jmin,this%cfg%jmax
               write(10, '(F20.12, ",", F20.12, ",", F20.12, ",", F20.12)') this%cfg%ym(j),this%Umean1em3(j)/this%timeint,&
               & sqrt(max(0.0_WP,this%U2mean1em3(j)/this%timeint-(this%Umean1em3(j)/this%timeint)**2)),this%U2mean1em3(j)/this%timeint
            end do
            close(unit=10)
         end if

      end if
   end subroutine record_stats
   
   !> Initialization of atom simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(atom), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_atom')
      
      ! Initialize the geometry
      call this%geometry_init()
      
      ! Initialize the simulation
      call this%simulation_init()
      
   end subroutine init


   !> Initialize geometry
   subroutine geometry_init(this)
      use sgrid_class, only: sgrid
      implicit none
      class(atom) :: this
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,ns_yz,ns_x
         real(WP) :: Lx,Ly,Lz,xshift,sratio_yz,sratio_x
         real(WP), dimension(:), allocatable :: x_uni,y_uni,z_uni
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x_uni(nx+1)); call this%input%read('X shift',xshift)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y_uni(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z_uni(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x_uni(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y_uni(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z_uni(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! Add stretching
         call this%input%read('Stretched cells in yz',ns_yz,default=0)
         if (ns_yz.gt.0) call this%input%read('Stretch ratio in yz',sratio_yz)
         call this%input%read('Stretched cells in x' ,ns_x ,default=0)
         if (ns_x .gt.0) call this%input%read('Stretch ratio in x' ,sratio_x )
         allocate(x(nx+1+1*ns_x )); x(      1:      1+nx)=x_uni
         allocate(y(ny+1+2*ns_yz)); y(ns_yz+1:ns_yz+1+ny)=y_uni
         allocate(z(nz+1+2*ns_yz)); z(ns_yz+1:ns_yz+1+nz)=z_uni
         do i=nx+2,nx+1+ns_x
            x(i)=x(i-1)+sratio_x*(x(i-1)-x(i-2))
         end do
         do j=ns_yz,1,-1
            y(j)=y(j+1)+sratio_yz*(y(j+1)-y(j+2))
         end do
         do j=ns_yz+2+ny,ny+1+2*ns_yz
            y(j)=y(j-1)+sratio_yz*(y(j-1)-y(j-2))
         end do
         do k=ns_yz,1,-1
            z(k)=z(k+1)+sratio_yz*(z(k+1)-z(k+2))
         end do
         do k=ns_yz+2+nz,nz+1+2*ns_yz
            z(k)=z(k-1)+sratio_yz*(z(k-1)-z(k-2))
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='atom')
         
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
         if (.not.this%cfg%amRoot) call this%plymesh%set_size(nvert=this%plymesh%nVert,npoly=this%plymesh%nPoly,nbeziertri=0)
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
      class(atom), intent(inout) :: this
      
      
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
         allocate(this%Umean1em3   (this%cfg%jmin:this%cfg%jmax));this%Umean1em3=0.0_WP
         allocate(this%U2mean1em3  (this%cfg%jmin:this%cfg%jmax));this%U2mean1em3=0.0_WP
         allocate(this%EPLmean (this%cfg%imin:this%cfg%imax));this%EPLmean=0.0_WP
         allocate(this%EPL2mean(this%cfg%imin:this%cfg%imax));this%EPL2mean=0.0_WP
         this%timeint=0.0_WP
      end block allocate_work_arrays
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: elvira,r2p,remap,r2pnet,r2p_cylinder
         integer :: i,j,k
         real(WP) :: xloc,rad
         ! Create a VOF solver with r2pnet
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=r2p_cylinder,transport_method=remap,name='VOF')
         this%vf%thin_thld_min=0.0_WP
         this%vf%flotsam_thld=0.0_WP
         ! Initialize to flat interface in liquid needle
         xloc=0.0_WP !< Interface initially at x=0
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                  if (this%vf%cfg%xm(i).lt.xloc.and.rad.le.0.5_WP*dl) then
                     this%vf%VF(i,j,k)=1.0_WP
                  else
                     this%vf%VF(i,j,k)=0.0_WP
                  end if
                  this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
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
      end block create_and_initialize_vof
           
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator

      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg,pcg_pfmg2
         use tpns_class,      only: dirichlet,clipped_neumann,slip
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set the flow properties
         call this%input%read('Liquid dynamic viscosity',this%fs%visc_l)
         call this%input%read('Gas dynamic viscosity'   ,this%fs%visc_g)
         call this%input%read('Liquid density',this%fs%rho_l)
         call this%input%read('Gas density'   ,this%fs%rho_g)
         call this%input%read('Surface tension coefficient',this%fs%sigma)
         ! Define gas and liquid inlet boundary conditions
         call this%fs%add_bcond(name='gas_inlet',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=gas_inlet)
         call this%fs%add_bcond(name='liq_inlet',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=liq_inlet)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=right_boundary)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Configure implicit velocity solver
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
      end block create_flow_solver
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mpi_f08,    only: MPI_ALLREDUCE,MPI_SUM
         use parallel,   only: MPI_REAL_WP
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k,ierr
         real(WP) :: Ugas,myAgas,Agas,Uliq,myAliq,Aliq
         real(WP) :: Qgas,Qliq,myU
         real(WP), parameter :: SLPM2SI=1.66667E-5_WP
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
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
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=+sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*Ugas
         end do
         ! Read in liquid flow rate and convert to SI
         call this%input%read('Liquid flow rate (SLPM)',Qliq)
         Qliq=Qliq*SLPM2SI
         ! Calculate liquid flow area - no overlap here!
         myAliq=0.0_WP
         call this%fs%get_bcond('liq_inlet',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAliq=myAliq+this%cfg%dy(j)*this%cfg%dz(k)*sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
         end do
         call MPI_ALLREDUCE(myAliq,Aliq,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         ! Calculate bulk axial velocity
         Uliq=Qliq/Aliq
         ! Apply Dirichlet at liquid injector port
         call this%fs%get_bcond('liq_inlet',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myU=2.0_WP*Uliq*(1.0_WP-min((this%fs%cfg%ym(j)**2+this%fs%cfg%zm(k)**2)/rl**2,1.0_WP))
            this%fs%U(i,j,k)=+sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*myU
         end do
         ! Sync the boundary velocities
         call this%cfg%sync(this%fs%U)
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

      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs
      
      ! Handle restart/saves here
      handle_restart: block
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use irl_fortran_interface
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         real(WP), dimension(:,:,:), allocatable :: Cr,Co1,Co2,Co3,Cn11,Cn12,Cn13,RT!,Cf
         real(WP), dimension(3) :: datum
         real(WP), dimension(3) :: ref1,ref2,ref3
         logical :: partfile_exists
         real(WP) :: n1, n2, n3
         real(WP) :: mag
         real(WP), dimension(3) :: v1
         real(WP) :: rad
         integer :: i,j,k,n
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call this%input%read('I/O partition',iopartition)

         call this%input%read('Restart from SGS lig',this%restart_from_sgslig,default=.false.)
         call this%input%read('Using SGS lig',this%using_sgslig,default=.false.)
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_atom_'//trim(timestamp))
            ! Read in the planes directly and set the IRL interface
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P21',var=P21)
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P22',var=P22)
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P23',var=P23)
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P24',var=P24)
            if (this%restart_from_sgslig) then
               allocate(Cr(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='Cr',var=Cr)
               allocate(Co1(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='Co1',var=Co1)
               allocate(Co2(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='Co2',var=Co2)
               allocate(Co3(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='Co3',var=Co3)
               allocate(Cn11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='C11',var=Cn11)
               allocate(Cn12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='C12',var=Cn12)
               allocate(Cn13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='C13',var=Cn13)
               allocate(RT(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='RT',var=RT)            
               !allocate(Cf(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='Cf',var=Cf)  
               do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
                  do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                     do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                        if (i.lt.this%vf%cfg%imin.or.i.gt.this%vf%cfg%imax.or.j.lt.this%vf%cfg%jmin.or.j.gt.this%vf%cfg%jmax.or.k.lt.this%vf%cfg%kmin.or.k.gt.this%vf%cfg%kmax) then
                           P11(i,j,k) = 0.0_WP; P12(i,j,k) = 0.0_WP; P13(i,j,k) = 0.0_WP; P14(i,j,k) = 0.0_WP
                           P21(i,j,k) = 0.0_WP; P22(i,j,k) = 0.0_WP; P23(i,j,k) = 0.0_WP; P24(i,j,k) = 0.0_WP
                           Cr(i,j,k)  = 0.0_WP; !Cf(i,j,k)  = 0.0_WP
                           Co1(i,j,k) = 0.0_WP; Co2(i,j,k) = 0.0_WP; Co3(i,j,k) = 0.0_WP
                           Cn11(i,j,k)= 0.0_WP; Cn12(i,j,k)= 0.0_WP; Cn13(i,j,k)= 0.0_WP; RT(i,j,k)  = 0.0_WP
                        end if
                        ! Check if the second plane is meaningful
                        if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                           call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                           call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                           call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
                        else if (Cr(i,j,k).gt.0.0_WP) then
                           call setAlignedCylinder(this%vf%liquid_gas_interface(i,j,k),1.0_WP,Cr(i,j,k))
                           datum(1) = Co1(i,j,k); datum(2) = Co2(i,j,k); datum(3) = Co3(i,j,k)
                           call setDatum(this%vf%liquid_gas_interface(i,j,k),datum)
                           ref1(1) = Cn11(i,j,k); ref1(2) = Cn12(i,j,k); ref1(3) = Cn13(i,j,k)
                           n1 = 0.0_WP; n2 = 0.0_WP; n3 = 0.0_WP; v1 = 0.0_WP
                           if (abs(ref1(1)) .ge. abs(ref1(2)) .and. abs(ref1(1)) .ge. abs(ref1(3))) then
                              n2 = ref1(1) / sqrt(ref1(2)**2 + ref1(1)**2)
                              n1 = (-n2 * ref1(2)) / ref1(1)
                              v1(1) = n1
                              v1(2) = n2
                              v1(3) = 0.0_WP
                           else if (abs(ref1(2)) .ge. abs(ref1(1)) .and. abs(ref1(2)) .ge. abs(ref1(3))) then
                              n1 = ref1(2) / sqrt(ref1(2)**2 + ref1(1)**2)
                              n2 = (-n1 * ref1(1)) / ref1(2)
                              v1(1) = n1
                              v1(2) = n2
                              v1(3) = 0.0_WP
                           else if (abs(ref1(3)) .ge. abs(ref1(1)) .and. abs(ref1(3)) .ge. abs(ref1(2))) then
                              n2 = ref1(3) / sqrt(ref1(2)**2 + ref1(3)**2)
                              n3 = (-n2 * ref1(2)) / ref1(3)
                              v1(1) = 0.0_WP
                              v1(2) = n2
                              v1(3) = n3
                           else
                              v1(1) = 0.0_WP; v1(2) = 0.0_WP; v1(3) = 0.0_WP
                           end if
                           ref3(1) = ref1(2) * v1(3) - ref1(3) * v1(2)
                           ref3(2) = ref1(3) * v1(1) - ref1(1) * v1(3)
                           ref3(3) = ref1(1) * v1(2) - ref1(2) * v1(1)
                           mag = sqrt(ref3(1)**2 + ref3(2)**2 + ref3(3)**2)
                           if (mag .ge. 1.0e-12) then
                              ref3 = ref3 / mag
                           end if
                           ref2(1) = ref3(2) * ref1(3) - ref3(3) * ref1(2)
                           ref2(2) = ref3(3) * ref1(1) - ref3(1) * ref1(3)
                           ref2(3) = ref3(1) * ref1(2) - ref3(2) * ref1(1)
                           mag = sqrt(ref2(1)**2 + ref2(2)**2 + ref2(3)**2)
                           if (mag .ge. 1.0e-12) then
                              ref2 = ref2 / mag
                           end if
                           call setReferenceFrame(this%vf%liquid_gas_interface(i,j,k),ref1,ref2,ref3)
                        else
                           call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                        end if
                        this%vf%recon_type(i,j,k)=RT(i,j,k)
                     end do
                  end do
               end do
            else
               do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
                  do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                     do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                        ! Check if the second plane is meaningful
                        if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                           call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                           call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                           call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
                        else
                           call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                        end if
                     end do
                  end do
               end do
            end if
            call this%vf%sync_interface()
            if (this%restart_from_sgslig) then
               deallocate(P11,P12,P13,P14,P21,P22,P23,P24,Cr,Co1,Co2,Co3,Cn11,Cn12,Cn13,RT)!,Cf)
            else
               deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            end if
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Update the band
            call this%vf%update_band()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! ! Calculate distance from polygons
            if (.not.this%restart_from_sgslig) then
               call this%vf%distance_from_polygon()
            end if
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Now read in the velocity solver data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
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
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            ! Update SGS viscousity
            call this%df%pull(name='LM',var=this%sgs%LM)
            call this%df%pull(name='MM',var=this%sgs%MM)
            if (this%using_sgslig.and. .not. this%restart_from_sgslig) then
               if (allocated(this%df%valname)) deallocate(this%df%valname)
               if (allocated(this%df%val    )) deallocate(this%df%val    )
               if (allocated(this%df%varname)) deallocate(this%df%varname)
               if (allocated(this%df%var    )) deallocate(this%df%var    )
               call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=25)
               this%df%valname=['t ','dt']
               this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24','LM ','MM ','Cr ','Co1','Co2','Co3','C11','C12','C13','RT ']
            end if
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            if (this%using_sgslig) then
               call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=25)
               this%df%valname=['t ','dt']
               this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24','LM ','MM ','Cr ','Co1','Co2','Co3','C11','C12','C13','RT ']
            else
               ! Prepare pardata object for saving restart files
               call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=17)
               this%df%valname=['t ','dt']
               this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24','LM ','MM ']
            end if
         end if
      end block handle_restart

      initialize_bu: block
         ! Initialize the breakup model
         call this%bu%init(this%input,this%cfg,this%vf,this%fs)
      end block initialize_bu

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
	  	   integer :: i,j,k,nplane,np
         this%smesh=surfmesh(nvar=5,name='plic')
		   this%smesh%varname(1)='recon_type'
		   this%smesh%varname(2)='thickness'
         this%smesh%varname(3)='curv'
         this%smesh%varname(4)='ccl_film'
         this%smesh%varname(5)='ccl_lig'
         call this%vf%update_surfmesh(this%smesh)
	      this%smesh%var(1,:)=0.0_WP
	      this%smesh%var(2,:)=0.0_WP    
         this%smesh%var(3,:)=0.0_WP
         this%smesh%var(4,:)=0.0_WP
         this%smesh%var(5,:)=0.0_WP    
         call add_surfgrid_variable(this,this%smesh,1,real(this%vf%recon_type,WP))  
         call add_surfgrid_variable(this,this%smesh,2,this%vf%thickness) 
         call add_surfgrid_variable(this,this%smesh,3,this%vf%curv) 
         call add_surfgrid_variable(this,this%smesh,4,real(this%bu%ccl_film%id,WP))  
         call add_surfgrid_variable(this,this%smesh,5,real(this%bu%ccl_lig%id,WP))  
      end block create_smesh

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='atom')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_surface('plic',this%smesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      create_stat: block
         use filesys,  only: makedir,isdir
         if (this%cfg%amRoot) then
            if (.not.isdir('stats')) call makedir('stats')
         end if
         this%stat_evt=event(time=this%time,name='stat output')
         call this%input%read('Stat output period',this%stat_evt%tper)
         call this%input%read('Stat begin',this%stat_begin)
      end block create_stat

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
         call this%mfile%add_column(this%vof_removed,'VOF removed')
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
      
      ! Create a timing monitor
      create_timing: block
         ! Create timers
         this%tstep=timer(comm=this%cfg%comm,name='Timestep')
         this%tvof =timer(comm=this%cfg%comm,name='VOFsolve')
         this%tvel =timer(comm=this%cfg%comm,name='Velocity')
         this%tpres=timer(comm=this%cfg%comm,name='Pressure')
         this%tdrop=timer(comm=this%cfg%comm,name='TranDrop')
         this%tfilm=timer(comm=this%cfg%comm,name='TranFilm')
         this%tlig =timer(comm=this%cfg%comm,name='TranLig')
         ! Create corresponding monitor file
         this%timefile=monitor(this%fs%cfg%amRoot,'timing')
         call this%timefile%add_column(this%time%n,'Timestep number')
         call this%timefile%add_column(this%time%t,'Time')
         call this%timefile%add_column(this%tstep%time ,trim(this%tstep%name))
         call this%timefile%add_column(this%tvof%time  ,trim(this%tvof%name))
         call this%timefile%add_column(this%tvel%time  ,trim(this%tvel%name))
         call this%timefile%add_column(this%tpres%time ,trim(this%tpres%name))
         call this%timefile%add_column(this%tdrop%time ,trim(this%tdrop%name))
         call this%timefile%add_column(this%tfilm%time ,trim(this%tfilm%name))
         call this%timefile%add_column(this%tlig%time  ,trim(this%tlig%name))
      end block create_timing

      contains
         
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


      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.ge.pg%imax-this%nlayer.or.&
         &   j.le.pg%jmin+this%nlayer.or.&
         &   j.ge.pg%jmax-this%nlayer.or.&
         &   k.le.pg%kmin+this%nlayer.or.&
         &   k.ge.pg%kmax-this%nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
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
      
   end subroutine simulation_init
   

   !> Take one time step
   subroutine step(this,lp)
      implicit none
      class(atom), intent(inout) :: this
      class(lpt), intent(inout) :: lp 
      
      call this%tstep%reset()
      call this%tvof%reset()
      call this%tvel%reset()
      call this%tpres%reset()
      call this%tdrop%reset()
      call this%tfilm%reset()
      call this%tlig%reset()
      call this%tstep%start()

      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF

      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
         
      ! VOF solver step
      call this%tvof%start()
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%tvof%stop()
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf)

      ! Turbulence modeling
      sgs_modeling: block
         use sgsmodel_class, only: vreman
         integer :: i,j,k
         this%resU=this%fs%rho_g
         call this%fs%get_gradu(this%gradU)
         call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
         do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_
            do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_
               do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
                  this%fs%visc(i,j,k)   =this%fs%visc(i,j,k)   +this%sgs%visc(i,j,k)
                  this%fs%visc_xy(i,j,k)=this%fs%visc_xy(i,j,k)+sum(this%fs%itp_xy(:,:,i,j,k)*this%sgs%visc(i-1:i,j-1:j,k))
                  this%fs%visc_yz(i,j,k)=this%fs%visc_yz(i,j,k)+sum(this%fs%itp_yz(:,:,i,j,k)*this%sgs%visc(i,j-1:j,k-1:k))
                  this%fs%visc_zx(i,j,k)=this%fs%visc_zx(i,j,k)+sum(this%fs%itp_xz(:,:,i,j,k)*this%sgs%visc(i-1:i,j,k-1:k))
               end do
            end do
         end do
      end block sgs_modeling
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         call this%tvel%start()
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
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply IB forcing to enforce BC at the pipe walls
         ibforcing: block
            integer :: i,j,k
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     if (this%fs%umask(i,j,k).eq.0) this%fs%U(i,j,k)=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*this%fs%U(i,j,k)
                     if (this%fs%vmask(i,j,k).eq.0) this%fs%V(i,j,k)=sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))*this%fs%V(i,j,k)
                     if (this%fs%wmask(i,j,k).eq.0) this%fs%W(i,j,k)=sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))*this%fs%W(i,j,k)
                  end do
               end do
            end do
            call this%fs%cfg%sync(this%fs%U)
            call this%fs%cfg%sync(this%fs%V)
            call this%fs%cfg%sync(this%fs%W)
         end block ibforcing
        
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         call this%tvel%stop()
         call this%tpres%start()

         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         call this%fs%add_surface_tension_jump_twoVF(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         call this%tpres%stop()
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! if (this%cfg%amRoot) print *, "before transfer"
      ! attempt transfter
      attempt_transfer : block
         call this%tfilm%start()
         if (this%bu%use_film_transfer) call this%bu%transfer_films(lp,this%time%t)
         call this%tfilm%stop()
         call this%tlig%start()
         if (this%bu%use_lig_transfer) call this%bu%transfer_ligs(lp,this%time%t)
         call this%tlig%stop()
         call this%tdrop%start()
         if (this%bu%use_drop_transfer) call this%bu%transfer_drops(lp,this%time%t)
         call this%tdrop%stop()
      end block attempt_transfer

      ! Recording stats
      call this%record_stats()

      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         real(WP) :: my_vof_removed
         my_vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            my_vof_removed=my_vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(my_vof_removed,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      end block remove_vof
      
      call this%tstep%stop()

      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update surface mesh
         update_smesh: block
            use irl_fortran_interface
            integer :: i,j,k,nplane,np
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh(this%smesh)
            ! Also populate nplane variable
            this%smesh%var(1,:)=0.0_WP
            this%smesh%var(2,:)=0.0_WP
            this%smesh%var(3,:)=0.0_WP
            this%smesh%var(4,:)=0.0_WP
            this%smesh%var(5,:)=0.0_WP
            call add_surfgrid_variable(this,this%smesh,1,real(this%vf%recon_type,WP))  
            call add_surfgrid_variable(this,this%smesh,2,this%vf%thickness) 
            call add_surfgrid_variable(this,this%smesh,3,this%vf%curv) 
            call add_surfgrid_variable(this,this%smesh,4,real(this%bu%ccl_film%id,WP))  
            call add_surfgrid_variable(this,this%smesh,5,real(this%bu%ccl_lig%id,WP))  
         end block update_smesh

         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      call this%timefile%write()
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
            real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
            real(WP), dimension(:,:,:), allocatable :: Cr,Co1,Co2,Co3,Cn11,Cn12,Cn13,RT!,Cf
            integer :: i,j,k
            real(WP), dimension(4) :: plane
            real(WP), dimension(3) :: aligned_cyl
            real(WP), dimension(3) :: datum
            real(WP), dimension(9) :: ref
            ! Handle IRL data
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));P11=0.0_WP
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));P12=0.0_WP
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));P13=0.0_WP
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));P14=0.0_WP
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));P21=0.0_WP
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));P22=0.0_WP
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));P23=0.0_WP
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));P24=0.0_WP
            if (this%using_sgslig) then
               allocate(Cr(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Cr=0.0_WP
               allocate(Co1(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Co1=0.0_WP
               allocate(Co2(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Co2=0.0_WP
               allocate(Co3(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Co3=0.0_WP
               allocate(Cn11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Cn11=0.0_WP
               allocate(Cn12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Cn12=0.0_WP
               allocate(Cn13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Cn13=0.0_WP
               allocate(RT(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); RT=0.0_WP
               ! allocate(Cf(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Cf=0.0_WP
               do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
                  do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                     do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                        RT(i,j,k)=2!this%vf%det%recon_type(i,j,k)
                        if (this%vf%recon_type(i,j,k).eq.1) then
                           aligned_cyl=getAlignedCylinder(this%vf%liquid_gas_interface(i,j,k))
                           datum=getDatum(this%vf%liquid_gas_interface(i,j,k))
                           ref=getReferenceFrame(this%vf%liquid_gas_interface(i,j,k))
                           Cr(i,j,k)=aligned_cyl(1)
                           ! Cf(i,j,k)=aligned_cyl(3)
                           Co1(i,j,k)=datum(1)
                           Co2(i,j,k)=datum(2)
                           Co3(i,j,k)=datum(3)
                           Cn11(i,j,k)=ref(1)
                           Cn12(i,j,k)=ref(2)
                           Cn13(i,j,k)=ref(3)
                        else
                           ! First plane
                           plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                           P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                           ! Second plane
                           plane=0.0_WP
                           if (getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.2) plane=getPlane(this%vf%liquid_gas_interface(i,j,k),1)
                           P21(i,j,k)=plane(1); P22(i,j,k)=plane(2); P23(i,j,k)=plane(3); P24(i,j,k)=plane(4)
                        end if
                     end do
                  end do
               end do
            else
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
            end if
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
            call this%df%push(name='LM', var=this%sgs%LM )
            call this%df%push(name='MM', var=this%sgs%MM )
            if (this%using_sgslig) then
               call this%df%push(name='Cr',var=Cr           )
               call this%df%push(name='Co1',var=Co1         )
               call this%df%push(name='Co2',var=Co2         )
               call this%df%push(name='Co3',var=Co3         )
               call this%df%push(name='C11',var=Cn11         )
               call this%df%push(name='C12',var=Cn12         )
               call this%df%push(name='C13',var=Cn13         )
               call this%df%push(name='RT',var=RT         )
            end if
            call this%df%write(fdata='restart/data_atom_'//trim(adjustl(timestamp)))
            ! Deallocate
            if (this%using_sgslig) then
               deallocate(P11,P12,P13,P14,P21,P22,P23,P24,Cr,Co1,Co2,Co3,Cn11,Cn12,Cn13,RT)
            else
               deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            end if
         end block save_restart
      end if
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(atom), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      
   end subroutine final
   
   !> Make a surface scalar variable
	subroutine add_surfgrid_variable(this,smesh,var_index,A)
		use irl_fortran_interface
		use vfs_class,only: VFhi,VFlo
		implicit none
      class(atom), intent(inout) :: this
		class(surfmesh), intent(inout) :: smesh
		integer, intent(in) :: var_index
		real(WP), dimension(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_), intent(in) :: A 
		integer :: i,j,k,n,shape,np,nplane,nbt
		
		! Fill out arrays
		if ((smesh%nPoly+smesh%nBezierTri).gt.0) then
		   np=0; nbt = 0
		   ! Start with quadratic surfaces
		   do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
			  do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
				 do i=this%vf%cfg%imin_,this%vf%cfg%imax_
					shape=getNumberOfTriangles(this%vf%interface_mixed_surface(i,j,k))
					if (shape.gt.0) then
					   do n=1,shape
						  nbt=nbt+1
						  smesh%var(var_index,nbt)=A(i,j,k)
					   end do
					end if
				 end do
			  end do
		   end do
		   ! Then do planes
		   do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
			  do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
				 do i=this%vf%cfg%imin_,this%vf%cfg%imax_
					do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
					   shape=getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k))
					   if (shape.gt.0) then
						  ! Increment polygon counter
						  np=np+1
						  ! Set nplane variable
						  smesh%var(var_index,nbt+np)=A(i,j,k)
					   end if
					end do
				 end do
			  end do
		   end do
		else
		   smesh%var(var_index,1)=1
		end if      
		
	 end subroutine add_surfgrid_variable
   

end module atom_class
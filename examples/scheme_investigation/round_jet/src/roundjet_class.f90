!> Definition for a round jet class
module roundjet_class
   use precision,         only: WP
   use config_class,      only: config
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: roundjet
   
   !> roundjet object
   type :: roundjet
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(incomp)      :: fs     !< Two-phase flow solver
      type(hypre_str)   :: ps     !< Structured Hypre linear solver for pressure
      type(timetracker) :: time   !< Time info
      
      !> Implicit solver
      logical     :: use_implicit !< Is an implicit solver used?
      type(ddadi) :: vs           !< DDADI solver for velocity
      
      !> SGS modeling
      logical        :: use_sgs   !< Is an LES model used?
      type(sgsmodel) :: sgs       !< SGS model for eddy viscosity
      
      !> Ensight postprocessing
      type(ensight)  :: ens_out   !< Ensight output for flow variables
      type(event)    :: ens_evt   !< Event trigger for Ensight output
      
      !> Simulation monitoring files
      type(monitor) :: mfile      !< General simulation monitoring
      type(monitor) :: cflfile    !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals

      !> Fluid definition
      real(WP) :: visc
      
   contains
      procedure :: init                            !< Initialize round jet simulation
      procedure :: step                            !< Advance round jet simulation by one time step
      procedure :: final                           !< Finalize round jet simulation
   end type roundjet
   
contains
   
   !> Initialization of roundjet simulation
   subroutine init(this)
      use param, only: param_read
      implicit none
      class(roundjet), intent(inout) :: this
      
      
      ! Initialize the config
      initialize_config: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         integer :: i,j,k,nx,ny,nz,ns_yz,ns_x
         real(WP) :: Lx,Ly,Lz,sratio_yz,sratio_x
         real(WP), dimension(:), allocatable :: x_uni,y_uni,z_uni
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         
         ! Read in grid definition
         call param_read('[Jet] Lx',Lx); call param_read('[Jet] nx',nx); allocate(x_uni(nx+1))
         call param_read('[Jet] Ly',Ly); call param_read('[Jet] ny',ny); allocate(y_uni(ny+1))
         call param_read('[Jet] Lz',Lz); call param_read('[Jet] nz',nz); allocate(z_uni(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x_uni(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y_uni(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z_uni(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x_uni,y=y_uni,z=z_uni,xper=.false.,yper=.false.,zper=.false.,name='jet')
         
         ! Read in partition
         call param_read('[Jet] Partition',partition)
         
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
         
         ! No walls in the atomization domain
         this%cfg%VF=1.0_WP

      end block initialize_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot,name='Jet')
         call param_read('[Jet] Max timestep size',this%time%dtmax)
         call param_read('[Jet] Max cfl number',this%time%cflmax)
         call param_read('[Jet] Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resU         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use incomp_class,      only: dirichlet,clipped_neumann,slip
         ! Create flow solver
         call this%fs%initialize(cfg=this%cfg,name='Incompressible NS')
         ! Set fluid properties
         call param_read('Density',this%fs%rho)
         call param_read('Dynamic viscosity',this%visc); this%fs%visc=this%visc

         ! Inflow on the left
         call this%fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp'  ,type=slip           ,face='y',dir=+1,canCorrect=.true. ,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym'  ,type=slip           ,face='y',dir=-1,canCorrect=.true. ,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp'  ,type=slip           ,face='z',dir=+1,canCorrect=.true. ,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm'  ,type=slip           ,face='z',dir=-1,canCorrect=.true. ,locator=zm_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call param_read('[Jet] Pressure iteration',this%ps%maxit)
         call param_read('[Jet] Pressure tolerance',this%ps%rcvg)
         ! Check if we want to use an implicit solver
         call param_read('[Jet] Use implicit solver',this%use_implicit)
         if (this%use_implicit) then
            ! Configure implicit velocity solver
            this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
            ! Setup the solver
            call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         else
            ! Setup the solver
            call this%fs%setup(pressure_solver=this%ps)
         end if
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         this%fs%Uf=0.0_WP; this%fs%Vf=0.0_WP; this%fs%Wf=0.0_WP
         ! Apply convective velocity
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            if (sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2).le.0.5_WP) then
               this%fs%Uf(i,j,k) =1.0_WP
               this%fs%U(i-1,j,k)=1.0_WP
               this%fs%U(i  ,j,k)=1.0_WP
            end if 
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%dt,'face')
         call this%fs%apply_bcond(this%time%dt,'cell')
         ! Compute MFR through all boundary conditions
         call this%fs%get_mfr()
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute divergence
         call this%fs%get_div()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('[Jet] Use SGS model',this%use_sgs)
         if (this%use_sgs) this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='jet')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('[Jet] Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%fs%U,this%fs%V,this%fs%W)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'jet')
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
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
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

   contains
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      implicit none
      class(roundjet), intent(inout) :: this

      
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      applyinflow : block
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            if (sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2).le.0.5_WP) then
               this%fs%Uf(i,j,k) =1.0_WP
               this%fs%U(i-1,j,k)=1.0_WP
               this%fs%U(i  ,j,k)=1.0_WP
            end if 
         end do
      end block applyinflow
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity()
      ! Turbulence modeling
      sgs_modeling: block
         use sgsmodel_class, only: vreman
         integer :: i,j,k
         this%resU=this%fs%rho
         call this%fs%get_gradu(this%gradU)
         call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
         this%fs%visc=this%visc+this%sgs%visc
         do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
                  this%fs%visc  (i,j,k)=this%fs%visc  (i,j,k)+this%sgs%visc(i,j,k)
                  this%fs%visc_x(i,j,k)=this%fs%visc_x(i,j,k)+sum(this%fs%itpr_x(:,i,j,k)*this%sgs%visc(i-1:i,j,k))
                  this%fs%visc_y(i,j,k)=this%fs%visc_y(i,j,k)+sum(this%fs%itpr_y(:,i,j,k)*this%sgs%visc(i,j-1:j,k))
                  this%fs%visc_z(i,j,k)=this%fs%visc_z(i,j,k)+sum(this%fs%itpr_z(:,i,j,k)*this%sgs%visc(i,j,k-1:k))
         end do; end do; end do
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
         
         ! ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Compute predictor U
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU; call this%cfg%sync(this%fs%U)
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV; call this%cfg%sync(this%fs%V)
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW; call this%cfg%sync(this%fs%W)
         
         
         call this%fs%update_faceU(this%fs%U,this%fs%V,this%fs%W,this%fs%Uf,this%fs%Vf,this%fs%Wf)
         call this%fs%update_pgrad_all(this%time%dt)
         call this%fs%apply_bcond(this%time%dt,'face')
         
         ! Solve Poisson equation
         call this%fs%correct_mfr()
         call this%fs%get_div() 
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div*this%fs%rho/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct face velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%Uf=this%fs%Uf-this%time%dt*this%resU/this%fs%rho
         this%fs%Vf=this%fs%Vf-this%time%dt*this%resV/this%fs%rho
         this%fs%Wf=this%fs%Wf-this%time%dt*this%resW/this%fs%rho
            
         ! Correct center velocity
         call this%fs%get_cell_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW,.true.)
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho; call this%cfg%sync(this%fs%U)
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho; call this%cfg%sync(this%fs%V)
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho; call this%cfg%sync(this%fs%W)
         call this%fs%apply_bcond(this%time%dt,'cell')

         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute divergence for monitoring
      call this%fs%get_div()

      ! Output to ensight
      if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)

      ! Perform and output monitoring
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
   contains
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(roundjet), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%gradU)
      
   end subroutine final
   
   
   !> Function that localizes the right (x+) of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   

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
   
   
   !> Function that localizes the front (z+) of the domain
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator
   
   
   !> Function that localizes the back (z-) of the domain
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator
   
   
end module roundjet_class
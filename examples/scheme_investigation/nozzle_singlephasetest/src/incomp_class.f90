!> Incompressible flow solver class:
!> Provides support for various BC, RHS calculation,
!> implicit solver, and pressure solution
!> Assumes constant viscosity and density.
module incomp_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use linsol_class,   only: linsol
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: incomp,bcond
   
   ! List of known available bcond types for this solver
   integer, parameter, public :: wall=1              !< Dirichlet at zero condition
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   integer, parameter, public :: convective=4        !< Convective outflow condition
   integer, parameter, public :: clipped_neumann=5   !< Clipped Neumann condition (outflow only)
   integer, parameter, public :: slip=6              !< Free-slip condition

   !> Boundary conditions for the incompressible solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      type(iterator) :: itr                               !< This is the iterator for the bcond - this identifies the (i,j,k)
      character(len=1) :: face                            !< Bcond face (x/y/z)
      integer :: dir                                      !< Bcond direction (+1,-1,0 for interior)
      real(WP) :: rdir                                    !< Bcond direction (real variable)
      logical :: canCorrect                               !< Can this bcond be corrected for global conservation?

      integer :: dir_rhs,dir_lhs,fdir_rhs,fdir_lhs  !< Keys for shift array
      character(len=2) :: celldir                   !< Helper string (e.g. 'xp')
   end type bcond
   integer, dimension(3,10), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1,-2,0,0,0,-2,0,0,0,-2,0,0,0],shape(shift))
   !> Incompressible solver object definition
   type :: incomp
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_INCOMP'  !< Solver name (default=UNNAMED_INCOMP)
      
      ! Constant property fluid
      real(WP) :: rho                                     !< This is our constant fluid density
      ! Viscosity fields
      real(WP), dimension(:,:,:), allocatable :: visc     !< Viscosity field on P-cell
      real(WP), dimension(:,:,:), allocatable :: visc_x   !< Viscosity field on U-cell
      real(WP), dimension(:,:,:), allocatable :: visc_y   !< Viscosity field on V-cell
      real(WP), dimension(:,:,:), allocatable :: visc_z   !< Viscosity field on W-cell
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      real(WP), dimension(:), allocatable :: mfr          !< MFR through each bcond
      real(WP), dimension(:), allocatable :: area         !< Area for each bcond
      real(WP) :: correctable_area                        !< Area of bcond that can be corrected
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! Flow variables
      real(WP), dimension(:,:,:), allocatable :: U        !< U velocity array
      real(WP), dimension(:,:,:), allocatable :: V        !< V velocity array
      real(WP), dimension(:,:,:), allocatable :: W        !< W velocity array
      real(WP), dimension(:,:,:), allocatable :: Uf       !< Uf velocity array
      real(WP), dimension(:,:,:), allocatable :: Vf       !< Vf velocity array
      real(WP), dimension(:,:,:), allocatable :: Wf       !< Wf velocity array
      real(WP), dimension(:,:,:), allocatable :: P        !< Pressure array
      
      ! Old flow variables
      real(WP), dimension(:,:,:), allocatable :: Uold     !< Uold velocity array
      real(WP), dimension(:,:,:), allocatable :: Vold     !< Vold velocity array
      real(WP), dimension(:,:,:), allocatable :: Wold     !< Wold velocity array
      
      ! Flow divergence
      real(WP), dimension(:,:,:), allocatable :: div      !< Divergence array
      
      ! Pressure solver
      class(linsol), pointer :: psolv                     !< Iterative linear solver object for the pressure Poisson equation
      
      ! Implicit velocity solver
      class(linsol), pointer :: implicit                  !< Iterative linear solver object for an implicit prediction of the NS residual
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: itpr_x,itpr_y,itpr_z   !< Interpolation for density
      real(WP), dimension(:,:,:,:), allocatable :: divp_x,divp_y,divp_z   !< Divergence for P-cell
      real(WP), dimension(:,:,:,:), allocatable :: divu_x,divv_y,divw_z   !< Divergence for W-cell
      real(WP), dimension(:,:,:,:), allocatable :: grdu_x,grdu_y,grdu_z   !< Velocity gradient for U
      real(WP), dimension(:,:,:,:), allocatable :: grdv_x,grdv_y,grdv_z   !< Velocity gradient for V
      real(WP), dimension(:,:,:,:), allocatable :: grdw_x,grdw_y,grdw_z   !< Velocity gradient for W
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable ::  mask                     !< Integer array used for modifying P metrics
      integer, dimension(:,:,:), allocatable :: umask                     !< Integer array used for modifying U metrics
      integer, dimension(:,:,:), allocatable :: vmask                     !< Integer array used for modifying V metrics
      integer, dimension(:,:,:), allocatable :: wmask                     !< Integer array used for modifying W metrics
      
      ! CFL numbers
      real(WP) :: CFLc_x,CFLc_y,CFLc_z                                    !< Convective CFL numbers
      real(WP) :: CFLv_x,CFLv_y,CFLv_z                                    !< Viscous CFL numbers
      
      ! Monitoring quantities
      real(WP) :: Umax,Vmax,Wmax,Pmax,divmax                              !< Maximum velocity, pressure, divergence
      
   contains
      procedure :: print=>incomp_print                    !< Output solver to the screen
      procedure :: initialize                             !< Initialize the flow solver
      procedure :: setup                                  !< Finish configuring the flow solver
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: get_dmomdt                             !< Calculate dmom/dt
      procedure :: get_div                                !< Calculate velocity divergence
      procedure :: get_pgrad                              !< Calculate pressure gradient
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_max                                !< Calculate maximum field values
      procedure :: get_strainrate                         !< Calculate deviatoric part of strain rate tensor
      procedure :: get_gradu                              !< Calculate velocity gradient tensor
      procedure :: get_vorticity                          !< Calculate vorticity tensor
      procedure :: get_mfr                                !< Calculate outgoing MFR through each bcond
      procedure :: correct_mfr                            !< Correct for mfr mismatch to ensure global conservation
      procedure :: shift_p                                !< Shift pressure to have zero average
      procedure :: solve_implicit                         !< Solve for the velocity residuals implicitly

      ! New Subroutines
      procedure :: apply_bcond_facepressure
      procedure :: get_viscosity
      procedure :: update_faceU
      procedure :: update_pgrad_all
      procedure :: get_cell_pgrad
   end type incomp
   
   
contains
   
   
   !> Default constructor for incompressible flow solver
   subroutine initialize(this,cfg,name)
      implicit none
      class(incomp), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to pgrid object
      this%cfg=>cfg
      
      ! Nullify bcond list
      this%nbc=0
      this%first_bc=>NULL()
      
      ! Allocate flow variables
      allocate(this%U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
      allocate(this%V(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%V=0.0_WP
      allocate(this%W(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%W=0.0_WP
      allocate(this%Uf(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Uf=0.0_WP
      allocate(this%Vf(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Vf=0.0_WP
      allocate(this%Wf(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Wf=0.0_WP
      allocate(this%P(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%P=0.0_WP
      
      ! Allocate flow divergence
      allocate(this%div(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%div=0.0_WP
      
      ! Allocate fluid viscosity
      allocate(this%visc  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc  =0.0_WP
      allocate(this%visc_x(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_x=0.0_WP
      allocate(this%visc_y(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_y=0.0_WP
      allocate(this%visc_z(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_z=0.0_WP
      
      ! Allocate old flow variables
      allocate(this%Uold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Uold=0.0_WP
      allocate(this%Vold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Vold=0.0_WP
      allocate(this%Wold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Wold=0.0_WP
      
      ! Prepare default metrics
      call this%init_metrics()
      
      ! Prepare P-cell masks
      allocate(this%mask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%mask=0
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.           1) this%mask(:this%cfg%imin-1,:,:)=2
         if (this%cfg%iproc.eq.this%cfg%npx) this%mask(this%cfg%imax+1:,:,:)=2
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.           1) this%mask(:,:this%cfg%jmin-1,:)=2
         if (this%cfg%jproc.eq.this%cfg%npy) this%mask(:,this%cfg%jmax+1:,:)=2
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.           1) this%mask(:,:,:this%cfg%kmin-1)=2
         if (this%cfg%kproc.eq.this%cfg%npz) this%mask(:,:,this%cfg%kmax+1:)=2
      end if
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%cfg%VF(i,j,k).eq.0.0_WP) this%mask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%mask)
      
      ! Prepare face mask for U
      allocate(this%umask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%umask=0
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.           1) this%umask(this%cfg%imin  ,:,:)=2
         if (this%cfg%iproc.eq.this%cfg%npx) this%umask(this%cfg%imax+1,:,:)=2
      end if
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               if (minval(this%cfg%VF(i-1:i,j,k)).eq.0.0_WP) this%umask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%umask)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%umask(this%cfg%imino,:,:)=this%umask(this%cfg%imino+1,:,:)
      
      ! Prepare face mask for V
      allocate(this%vmask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%vmask=0
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.           1) this%vmask(:,this%cfg%jmin  ,:)=2
         if (this%cfg%jproc.eq.this%cfg%npy) this%vmask(:,this%cfg%jmax+1,:)=2
      end if
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (minval(this%cfg%VF(i,j-1:j,k)).eq.0.0_WP) this%vmask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%vmask)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%vmask(:,this%cfg%jmino,:)=this%vmask(:,this%cfg%jmino+1,:)
      
      ! Prepare face mask for W
      allocate(this%wmask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%wmask=0
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.           1) this%wmask(:,:,this%cfg%kmin  )=2
         if (this%cfg%kproc.eq.this%cfg%npz) this%wmask(:,:,this%cfg%kmax+1)=2
      end if
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (minval(this%cfg%VF(i,j,k-1:k)).eq.0.0_WP) this%wmask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%wmask)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%wmask(:,:,this%cfg%kmino)=this%wmask(:,:,this%cfg%kmino+1)
      
   end subroutine initialize
      
   
   !> Metric initialization with no awareness of walls nor bcond
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k,st1,st2
      real(WP), dimension(-1:0) :: itpx,itpy,itpz
      
      ! Allocate finite difference density interpolation coefficients to cell faces
      allocate(this%itpr_x(-1:0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< X-face-centered
      allocate(this%itpr_y(-1:0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Y-face-centered
      allocate(this%itpr_z(-1:0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Z-face-centered
      ! Create density interpolation coefficients to cell face in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%itpr_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
            end do
         end do
      end do
      ! Create density interpolation coefficients to cell face in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpr_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
            end do
         end do
      end do
      ! Create density interpolation coefficients to cell face in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpr_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do

      ! Allocate finite volume divergence operators
      allocate(this%divp_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divp_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divp_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm] or tangent to cell face
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%divp_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divp_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divp_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
            end do
         end do
      end do

      allocate(this%divu_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (x)
      allocate(this%divv_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (y)
      allocate(this%divw_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Face-centered (z)
      ! Create divergence operator perpendicular to cell face [x ,ym,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%divu_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      ! Create divergence operator perpendicular to cell face [xm,y ,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%divv_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      ! Create divergence operator perpendicular to cell face [xm,ym,z ]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%divw_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      ! Allocate finite difference velocity gradient operators
      allocate(this%grdu_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdv_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdw_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%grdw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%grdu_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%grdw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (yz)
      allocate(this%grdu_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%grdv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (yz)
      ! Create gradient coefficients to cell center [xm,ym,zm]
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%grdu_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of U from [x ,ym,zm]
               this%grdv_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y of V from [xm,y ,zm]
               this%grdw_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%grdv_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of V from [xm,y ,zm]
               this%grdw_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%grdu_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of U from [x ,ym,zm]
               this%grdv_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of V from [xm,y ,zm]
            end do
         end do
      end do
      
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%grdu_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of U from [x ,ym,zm]
               this%grdv_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of V from [xm,y ,zm]
            end do
         end do
      end do
   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls
   subroutine adjust_metrics(this)
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k,st1,st2
      real(WP) :: delta,mysum
      
      ! Sync up u/v/wmasks
      call this%cfg%sync(this%umask)
      call this%cfg%sync(this%vmask)
      call this%cfg%sync(this%wmask)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%umask(this%cfg%imino,:,:)=this%umask(this%cfg%imino+1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%vmask(:,this%cfg%jmino,:)=this%vmask(:,this%cfg%jmino+1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%wmask(:,:,this%cfg%kmino)=this%wmask(:,:,this%cfg%kmino+1)
      
      ! I am assuming here that we do not really need to zero out wall cells
      ! as they could be used for Dirichlet (then the density needs to be available! could be problematic if we do not have an explicit BC for scalars, e.g. for a Couette flow)
      ! or outflow condition (then the density needs to be available but it should be directly calculated)
      ! or used for a real no-slip wall (then density is always multiplied by zero)
      ! Adjust density interpolation coefficients to cell faces in the presence of walls (only walls!)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1
              ! Linear interpolation in x
              if (this%mask(i,j,k).eq.0.and.this%mask(i-1,j,k).gt.0) this%itpr_x(:,i,j,k)=[0.0_WP,1.0_WP]
              if (this%mask(i,j,k).gt.0.and.this%mask(i-1,j,k).eq.0) this%itpr_x(:,i,j,k)=[1.0_WP,0.0_WP]
              ! Linear interpolation in y               
              if (this%mask(i,j,k).eq.0.and.this%mask(i,j-1,k).gt.0) this%itpr_y(:,i,j,k)=[0.0_WP,1.0_WP]
              if (this%mask(i,j,k).gt.0.and.this%mask(i,j-1,k).eq.0) this%itpr_y(:,i,j,k)=[1.0_WP,0.0_WP]
              ! Linear interpolation in z
              if (this%mask(i,j,k).eq.0.and.this%mask(i,j,k-1).gt.0) this%itpr_z(:,i,j,k)=[0.0_WP,1.0_WP]
              if (this%mask(i,j,k).gt.0.and.this%mask(i,j,k-1).eq.0) this%itpr_z(:,i,j,k)=[1.0_WP,0.0_WP]
           end do
        end do
      end do
      
      ! Loop over the domain and adjust divergence for P cell
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).gt.0) then
                  this%divp_x(:,i,j,k)=0.0_WP
                  this%divp_y(:,i,j,k)=0.0_WP
                  this%divp_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to U metrics
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               if (this%umask(i,j,k).gt.0) then
                  this%divu_x(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to V metrics
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (this%vmask(i,j,k).gt.0) then
                  this%divv_y(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to W metrics
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (this%wmask(i,j,k).gt.0) then
                  this%divw_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      

      ! Adjust gradient coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               ! FD gradient in x of V from [xm,y ,zm]
               if (maxval(this%vmask(i-1:i,j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%vmask(i  ,j,k).eq.0) delta=delta+(this%cfg%xm(i)-this%cfg%x (i  ))
                  if (this%vmask(i-1,j,k).eq.0) delta=delta+(this%cfg%x (i)-this%cfg%xm(i-1))
                  if (delta.gt.0.0_WP) then
                     this%grdv_x(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdv_x(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in x of W from [xm,ym,z ]
               if (maxval(this%wmask(i-1:i,j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%wmask(i  ,j,k).eq.0) delta=delta+(this%cfg%xm(i)-this%cfg%x (i  ))
                  if (this%wmask(i-1,j,k).eq.0) delta=delta+(this%cfg%x (i)-this%cfg%xm(i-1))
                  if (delta.gt.0.0_WP) then
                     this%grdw_x(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdw_x(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! FD gradient in y of U from [x ,ym,zm]
               if (maxval(this%umask(i,j-1:j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%umask(i,j  ,k).eq.0) delta=delta+(this%cfg%ym(j)-this%cfg%y (j  ))
                  if (this%umask(i,j-1,k).eq.0) delta=delta+(this%cfg%y (j)-this%cfg%ym(j-1))
                  if (delta.gt.0.0_WP) then
                     this%grdu_y(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdu_y(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in y of W from [xm,ym,z ]
               if (maxval(this%wmask(i,j-1:j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%wmask(i,j  ,k).eq.0) delta=delta+(this%cfg%ym(j)-this%cfg%y (j  ))
                  if (this%wmask(i,j-1,k).eq.0) delta=delta+(this%cfg%y (j)-this%cfg%ym(j-1))
                  if (delta.gt.0.0_WP) then
                     this%grdw_y(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdw_y(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! FD gradient in z of U from [x ,ym,zm]
               if (maxval(this%umask(i,j,k-1:k)).gt.0) then
                  delta=0.0_WP
                  if (this%umask(i,j,k  ).eq.0) delta=delta+(this%cfg%zm(k)-this%cfg%z (k  ))
                  if (this%umask(i,j,k-1).eq.0) delta=delta+(this%cfg%z (k)-this%cfg%zm(k-1))
                  if (delta.gt.0.0_WP) then
                     this%grdu_z(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdu_z(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in z of V from [xm,y ,zm]
               if (maxval(this%vmask(i,j,k-1:k)).gt.0) then
                  delta=0.0_WP
                  if (this%vmask(i,j,k  ).eq.0) delta=delta+(this%cfg%zm(k)-this%cfg%z (k  ))
                  if (this%vmask(i,j,k-1).eq.0) delta=delta+(this%cfg%z (k)-this%cfg%zm(k-1))
                  if (delta.gt.0.0_WP) then
                     this%grdv_z(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdv_z(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do

      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%divp_x=0.0_WP
         this%divu_x=0.0_WP
         this%grdu_x=0.0_WP
         this%grdv_x=0.0_WP
         this%grdw_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%divp_y=0.0_WP
         this%divv_y=0.0_WP
         this%grdu_y=0.0_WP
         this%grdv_y=0.0_WP
         this%grdw_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%divp_z=0.0_WP
         this%divw_z=0.0_WP
         this%grdu_z=0.0_WP
         this%grdv_z=0.0_WP
         this%grdw_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the flow solver now that bconds have been defined
   subroutine setup(this,pressure_solver,implicit_solver)
      implicit none
      class(incomp), intent(inout) :: this
      class(linsol), target, intent(in) :: pressure_solver                      !< A pressure solver is required
      class(linsol), target, intent(in), optional :: implicit_solver            !< An implicit solver can be provided
      integer :: i,j,k
      
      ! Adjust metrics based on bcflag array
      call this%adjust_metrics()
      
      ! Point to pressure solver linsol object
      this%psolv=>pressure_solver
      
      ! Set 7-pt stencil map for the pressure solver
      this%psolv%stc(1,:)=[ 0, 0, 0]
      this%psolv%stc(2,:)=[+1, 0, 0]
      this%psolv%stc(3,:)=[-1, 0, 0]
      this%psolv%stc(4,:)=[ 0,+1, 0]
      this%psolv%stc(5,:)=[ 0,-1, 0]
      this%psolv%stc(6,:)=[ 0, 0,+1]
      this%psolv%stc(7,:)=[ 0, 0,-1]
      
      ! Setup the scaled Laplacian operator from incomp metrics: lap(*)=-vol*div(grad(*))
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Set Laplacian
               this%psolv%opr(1,i,j,k)=this%divp_x(1,i,j,k)*this%divu_x(-1,i+1,j,k)+&
               &                       this%divp_x(0,i,j,k)*this%divu_x( 0,i  ,j,k)+&
               &                       this%divp_y(1,i,j,k)*this%divv_y(-1,i,j+1,k)+&
               &                       this%divp_y(0,i,j,k)*this%divv_y( 0,i,j  ,k)+&
               &                       this%divp_z(1,i,j,k)*this%divw_z(-1,i,j,k+1)+&
               &                       this%divp_z(0,i,j,k)*this%divw_z( 0,i,j,k  )
               this%psolv%opr(2,i,j,k)=this%divp_x(1,i,j,k)*this%divu_x( 0,i+1,j,k)
               this%psolv%opr(3,i,j,k)=this%divp_x(0,i,j,k)*this%divu_x(-1,i  ,j,k)
               this%psolv%opr(4,i,j,k)=this%divp_y(1,i,j,k)*this%divv_y( 0,i,j+1,k)
               this%psolv%opr(5,i,j,k)=this%divp_y(0,i,j,k)*this%divv_y(-1,i,j  ,k)
               this%psolv%opr(6,i,j,k)=this%divp_z(1,i,j,k)*this%divw_z( 0,i,j,k+1)
               this%psolv%opr(7,i,j,k)=this%divp_z(0,i,j,k)*this%divw_z(-1,i,j,k  )
               ! Scale it by the cell volume
               this%psolv%opr(:,i,j,k)=-this%psolv%opr(:,i,j,k)*this%cfg%vol(i,j,k)
            end do
         end do
      end do
      
      ! Initialize the pressure Poisson solver
      call this%psolv%init()
      call this%psolv%setup()
      
      ! Prepare implicit solver if it had been provided
      if (present(implicit_solver)) then
         
         ! Point to implicit solver linsol object
         this%implicit=>implicit_solver
         
         ! Set 7-pt stencil map for the velocity solver
         this%implicit%stc(1,:)=[ 0, 0, 0]
         this%implicit%stc(2,:)=[+1, 0, 0]
         this%implicit%stc(3,:)=[-1, 0, 0]
         this%implicit%stc(4,:)=[ 0,+1, 0]
         this%implicit%stc(5,:)=[ 0,-1, 0]
         this%implicit%stc(6,:)=[ 0, 0,+1]
         this%implicit%stc(7,:)=[ 0, 0,-1]
         
         ! Set the diagonal to 1 to make sure all cells participate in solver
         this%implicit%opr(1,:,:,:)=1.0_WP
         
         ! Initialize the implicit velocity solver
         call this%implicit%init()
         
      else
         
         ! Point to implicit solver linsol object
         this%implicit=>NULL()
         
      end if
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,face,dir,canCorrect)
      use string,         only: lowercase
      use messager,       only: die
      use iterator_class, only: locator_ftype
      implicit none
      class(incomp), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: type
      procedure(locator_ftype) :: locator
      character(len=1), intent(in) :: face
      integer, intent(in) :: dir
      logical, intent(in) :: canCorrect
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type

      ! 1. Store Direction and Correction capabilities (incomp specific)
      select case (dir) ! Outward-oriented
      case (+1); new_bc%dir = +1
      case (-1); new_bc%dir = -1
      case ( 0); new_bc%dir =  0
      case default; call die('[incomp add_bcond] Unknown bcond dir - expecting -1, +1, or 0')
      end select
      new_bc%rdir = real(new_bc%dir, WP)
      new_bc%canCorrect = canCorrect

      ! 2. Set up Shift Keys (MAST Logic for Collocated Updates)
      select case (lowercase(face))
      case ('x')
         new_bc%face = 'x'
         ! Logic: if dir=+1 (Right), LHS=10(Self), RHS=1(Right Neighbor). Ghost copies Neighbor.
         new_bc%fdir_rhs = 1 + max(0, -new_bc%dir)
         new_bc%dir_lhs  = (1 - max(0, -new_bc%dir))*10 + max(0, -new_bc%dir)*1
         new_bc%dir_rhs  = (1 - max(0, -new_bc%dir))*1  + max(0, -new_bc%dir)*10
      case ('y')
         new_bc%face = 'y'
         new_bc%fdir_rhs = 3 + max(0, -new_bc%dir)
         new_bc%dir_lhs  = (1 - max(0, -new_bc%dir))*10 + max(0, -new_bc%dir)*3
         new_bc%dir_rhs  = (1 - max(0, -new_bc%dir))*3  + max(0, -new_bc%dir)*10
      case ('z')
         new_bc%face = 'z'
         new_bc%fdir_rhs = 5 + max(0, -new_bc%dir)
         new_bc%dir_lhs  = (1 - max(0, -new_bc%dir))*10 + max(0, -new_bc%dir)*5
         new_bc%dir_rhs  = (1 - max(0, -new_bc%dir))*5  + max(0, -new_bc%dir)*10
      case default
         call die('[incomp add_bcond] Unknown bcond face - expecting x, y, or z')
      end select
      
      ! Default LHS face shift (usually 10 = [0,0,0], meaning no shift)
      new_bc%fdir_lhs = 10 

      ! 3. Create Iterator
      new_bc%itr = iterator(pg=this%cfg, name=new_bc%name, locator=locator, face=new_bc%face)
      
      ! 4. Insert into Linked List
      new_bc%next => this%first_bc
      this%first_bc => new_bc
      
      ! Increment counter
      this%nbc = this%nbc + 1
      
      ! 5. Apply Masks immediately (if Dirichlet)
      select case (new_bc%type)
      case (dirichlet) ! Added 'wall' here as it implies fixed value
         
         do n = 1, new_bc%itr%n_
            i = new_bc%itr%map(1,n)
            j = new_bc%itr%map(2,n)
            k = new_bc%itr%map(3,n)
            
            ! Mask Face Velocity (Scalar/Staggered Face)
            select case (new_bc%face)
            case ('x')
               this%umask(i,j,k) = 2
               ! Mask Cell Center (Pressure/Collocated Velocity)
               ! Logic: Shift index to target the Ghost Cell if iterator points to an interior face
               this%mask(i + min(0, new_bc%dir), j, k) = 2
            case ('y')
               this%vmask(i,j,k) = 2
               this%mask(i, j + min(0, new_bc%dir), k) = 2
            case ('z')
               this%wmask(i,j,k) = 2
               this%mask(i, j, k + min(0, new_bc%dir)) = 2
            end select
         end do
         
      case (neumann)
         ! No masking needed for Neumann
      case (clipped_neumann)
         ! No masking needed
      case (slip)
         ! do n = 1, new_bc%itr%n_
         !    i = new_bc%itr%map(1,n); j = new_bc%itr%map(2,n); k = new_bc%itr%map(3,n)
            
         !    select case (new_bc%face)
         !    case ('x'); this%umask(i,j,k) = 0
         !    case ('y'); this%vmask(i,j,k) = 0
         !    case ('z'); this%wmask(i,j,k) = 0
         !    end select
         !    ! Note: We do NOT touch this%mask here.
         ! end do
         ! Slip usually requires flow, so we generally don't mask the pressure solver
      case (convective)
      case default
         call die('[incomp add_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(incomp), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[incomp get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,dt,scope)
      use messager, only: die
      use string,   only: lowercase
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      character(len=*), optional, intent(in) :: scope
      integer :: i,j,k,n,stag
      type(bcond), pointer :: my_bc
      logical :: do_cell, do_face
      character(len=10) :: mode
      real(WP) :: flip_u, flip_v, flip_w

      ! 1. Determine Scope
      if (present(scope)) then
         mode = lowercase(scope)
      else
         mode = 'all'
      end if

      select case (trim(mode))
      case ('cell')
         do_cell = .true.; do_face = .false.
      case ('face')
         do_cell = .false.; do_face = .true.
      case ('all')
         do_cell = .true.; do_face = .true.
      case default
         call die('[incomp apply_bcond] Unknown scope - use cell, face, or all')
      end select
      ! 2. Traverse Boundary List
      my_bc => this%first_bc
      do while (associated(my_bc))
         ! Only execute if this process owns the boundary cells
         if (my_bc%itr%amIn) then
            select case (my_bc%type)
            ! Dirichlet/Wall: Values are fixed. 
            ! We assume they are set during init and masked out, so we skip them here.
            case (dirichlet)
               ! Do nothing. Masks set in add_bcond protect these values.
            ! Neumann/Slip/Outflow: We must actively update the boundary values.
            case (neumann, clipped_neumann, slip, convective)
               ! =========================================================
               ! A. UPDATE CELL-CENTERED VELOCITIES (Ghost Cells)
               ! =========================================================
               if (do_cell) then
                  ! 1. Determine Reflection Factors based on BC Type and Face
                  ! Default: Copy everything (Neumann/Outflow behavior)
                  flip_u = 1.0_WP
                  flip_v = 1.0_WP
                  flip_w = 1.0_WP

                  ! Special handling for Free Slip:
                  ! Reflect (flip sign) the Normal component.
                  ! Copy (keep sign) the Tangential component.
                  if (my_bc%type .eq. slip) then
                     select case (my_bc%face)
                     case ('x'); flip_u = -1.0_WP
                     case ('y'); flip_v = -1.0_WP
                     case ('z'); flip_w = -1.0_WP
                     end select
                  end if
         
                  ! 2. Apply Update
                  do n = 1, my_bc%itr%n_
                     i = my_bc%itr%map(1,n); j = my_bc%itr%map(2,n); k = my_bc%itr%map(3,n)
                     
                     ! Update U (Normal or Tangential depending on face)
                     this%U(i-shift(1,my_bc%dir_lhs), j-shift(2,my_bc%dir_lhs), k-shift(3,my_bc%dir_lhs)) &
                        = flip_u * this%U(i-shift(1,my_bc%dir_rhs), j-shift(2,my_bc%dir_rhs), k-shift(3,my_bc%dir_rhs))
         
                     ! Update V
                     this%V(i-shift(1,my_bc%dir_lhs), j-shift(2,my_bc%dir_lhs), k-shift(3,my_bc%dir_lhs)) &
                        = flip_v * this%V(i-shift(1,my_bc%dir_rhs), j-shift(2,my_bc%dir_rhs), k-shift(3,my_bc%dir_rhs))
         
                     ! Update W
                     this%W(i-shift(1,my_bc%dir_lhs), j-shift(2,my_bc%dir_lhs), k-shift(3,my_bc%dir_lhs)) &
                        = flip_w * this%W(i-shift(1,my_bc%dir_rhs), j-shift(2,my_bc%dir_rhs), k-shift(3,my_bc%dir_rhs))
                     
                     ! Clipping for cell-centered velocities (Normal component only for simplicity/robustness)
                     if (my_bc%type == clipped_neumann) then
                        select case (my_bc%face)
                        case ('x')
                           if (this%U(i-shift(1,my_bc%dir_lhs), j, k) * my_bc%rdir .lt. 0.0_WP) &
                              this%U(i-shift(1,my_bc%dir_lhs), j, k) = -this%U(i-shift(1,my_bc%dir_rhs), j-shift(2,my_bc%dir_rhs), k-shift(3,my_bc%dir_rhs))
                        case ('y')
                           if (this%V(i, j-shift(2,my_bc%dir_lhs), k) * my_bc%rdir .lt. 0.0_WP) &
                              this%V(i, j-shift(2,my_bc%dir_lhs), k) = -this%V(i-shift(1,my_bc%dir_rhs), j-shift(2,my_bc%dir_rhs), k-shift(3,my_bc%dir_rhs))
                        case ('z')
                           if (this%W(i, j, k-shift(3,my_bc%dir_lhs)) * my_bc%rdir .lt. 0.0_WP) &
                              this%W(i, j, k-shift(3,my_bc%dir_lhs)) = -this%W(i-shift(1,my_bc%dir_rhs), j-shift(2,my_bc%dir_rhs), k-shift(3,my_bc%dir_rhs))
                        end select
                     end if
                  end do
               end if

               ! =========================================================
               ! B. UPDATE FACE-CENTERED VELOCITIES (Boundary Fluxes)
               ! =========================================================
               if (do_face) then
                  stag = min(my_bc%dir, 0)
                  
                  select case (my_bc%face)
                  case ('x')
                     do n = 1, my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        
                        ! 1. Normal Component (Uf)
                        if (my_bc%type == slip) then
                           this%Uf(i,j,k) = 0.0_WP  ! No penetration
                        else
                           this%Uf(i,j,k) = this%Uf(i-my_bc%dir, j, k) ! Zero gradient
                        end if
                        
                        ! 2. Tangential Components (Vf, Wf)
                        this%Vf(i+stag, j:j+1, k) = this%Vf(i-my_bc%dir+stag, j:j+1, k)
                        this%Wf(i+stag, j, k:k+1) = this%Wf(i-my_bc%dir+stag, j, k:k+1)
                        
                        ! 3. Clipping (Outflow only)
                        if (my_bc%type == clipped_neumann) then
                           if (this%Uf(i,j,k) * my_bc%rdir .lt. 0.0_WP) this%Uf(i,j,k) = 0.0_WP
                        end if
                     end do
                     
                  case ('y')
                     do n = 1, my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        
                        ! 1. Normal Component (Vf)
                        if (my_bc%type == slip) then
                           this%Vf(i,j,k) = 0.0_WP
                        else
                           this%Vf(i,j,k) = this%Vf(i, j-my_bc%dir, k)
                        end if
                        
                        ! 2. Tangential Components (Uf, Wf)
                        this%Uf(i:i+1, j+stag, k) = this%Uf(i:i+1, j-my_bc%dir+stag, k)
                        this%Wf(i, j+stag, k:k+1) = this%Wf(i, j-my_bc%dir+stag, k:k+1)
                        
                        ! 3. Clipping
                        if (my_bc%type == clipped_neumann) then
                           if (this%Vf(i,j,k) * my_bc%rdir .lt. 0.0_WP) this%Vf(i,j,k) = 0.0_WP
                        end if
                     end do
                     
                  case ('z')
                     do n = 1, my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        
                        ! 1. Normal Component (Wf)
                        if (my_bc%type == slip) then
                           this%Wf(i,j,k) = 0.0_WP
                        else
                           this%Wf(i,j,k) = this%Wf(i, j, k-my_bc%dir)
                        end if
                        
                        ! 2. Tangential Components (Uf, Vf)
                        this%Uf(i:i+1, j, k+stag) = this%Uf(i:i+1, j, k-my_bc%dir+stag)
                        this%Vf(i, j:j+1, k+stag) = this%Vf(i, j:j+1, k-my_bc%dir+stag)
                        
                        ! 3. Clipping
                        if (my_bc%type == clipped_neumann) then
                           if (this%Wf(i,j,k) * my_bc%rdir .lt. 0.0_WP) this%Wf(i,j,k) = 0.0_WP
                        end if
                     end do
                  end select
               end if ! End do_face logic

            case default
               call die('[incomp apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
      ! 3. Sync Fields (Only what was updated)
      if (do_cell) then
         call this%cfg%sync(this%U)
         call this%cfg%sync(this%V)
         call this%cfg%sync(this%W)
      end if
      
      if (do_face) then
         call this%cfg%sync(this%Uf)
         call this%cfg%sync(this%Vf)
         call this%cfg%sync(this%Wf)
      end if
      
   end subroutine apply_bcond
   
   
   !> Enforce boundary conditions on face pressure with Hydrostatic Correction
   subroutine apply_bcond_facepressure(this,P,Px,Py,Pz)
      use messager, only: die
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: P
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Px,Py,Pz
      integer :: i,j,k,n
      type(bcond), pointer :: my_bc
      real(WP) :: hydro_inc, dist

      ! Traverse bcond list
      my_bc=>this%first_bc
      do while (associated(my_bc))
         ! Only processes inside the bcond work here
         if (my_bc%itr%amIn) then
            ! Select appropriate action based on the bcond type
            select case (my_bc%type)
            case (dirichlet)               !< Apply Dirichlet conditions (nothing to do here)
            case (neumann)                 !< no treatment of face pressure with neumann
            
            ! =======================================================
            ! FREE-SLIP: Enforce Hydrostatic Balance (P_face = P_cell + rho*g*h)
            ! =======================================================
            case (slip)                    
            !   do n=1,my_bc%itr%n_
            !     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                
            !     select case (my_bc%face)
            !     case ('x')
            !       ! Distance from Cell Center to Face is 0.5 * dx * direction (+1/-1)
            !       dist = 0.5_WP * real(my_bc%dir, WP) * this%cfg%dx(i)
            !       hydro_inc = this%rho * this%gravity(1) * dist
                  
            !       ! Set Face Pressure based on Cell Pressure + Hydrostatic Jump
            !       Px(i-shift(1,my_bc%fdir_lhs),j,k) = P(i,j,k) + hydro_inc
                  
            !     case ('y')
            !       dist = 0.5_WP * real(my_bc%dir, WP) * this%cfg%dy(j)
            !       hydro_inc = this%rho * this%gravity(2) * dist
                  
            !       Py(i,j-shift(2,my_bc%fdir_lhs),k) = P(i,j,k) + hydro_inc
                  
            !     case ('z')
            !       dist = 0.5_WP * real(my_bc%dir, WP) * this%cfg%dz(k)
            !       hydro_inc = this%rho * this%gravity(3) * dist
                  
            !       Pz(i,j,k-shift(3,my_bc%fdir_lhs)) = P(i,j,k) + hydro_inc
                  
            !     end select
            !   end do

            case (clipped_neumann)         !< Apply clipped Neumann condition
              do n=1,my_bc%itr%n_
                i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                select case (my_bc%face)
                case ('x')
                  Px   (i-shift(1,my_bc%fdir_lhs),j,k) &
                    =Px(i-shift(1,my_bc%fdir_rhs),j,k)
                case ('y')
                  Py   (i,j-shift(2,my_bc%fdir_lhs),k) &
                    =Py(i,j-shift(2,my_bc%fdir_rhs),k)
                case ('z')
                  Pz   (i,j,k-shift(3,my_bc%fdir_lhs)) &
                    =Pz(i,j,k-shift(3,my_bc%fdir_rhs))
                end select
              end do
            !case (convective)   ! Not implemented yet!
            !case (adiabatic)    ! Not implemented yet!
            case default
               call die('[mast apply_bcond_facepressure] Unknown bcond type')
            end select
         end if
         ! Move on to the next bcond
         my_bc=>my_bc%next
      end do
      
      ! Sync fields
      call this%cfg%sync(Px)
      call this%cfg%sync(Py)
      call this%cfg%sync(Pz)
      
   end subroutine apply_bcond_facepressure
   
   
   !> Calculate the explicit momentum time derivative based on U/V/W/P
   subroutine get_dmomdt(this,drhoUdt,drhoVdt,drhoWdt)
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      
      ! Zero out drhoUVW/dt arrays
      drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Flux of rhoU
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FX(i,j,k)=-this%rho*this%Uf(i,j,k)*sum(this%itpr_x(:,i,j,k)*this%U(i-1:i,j,k))
               FY(i,j,k)=-this%rho*this%Vf(i,j,k)*sum(this%itpr_y(:,i,j,k)*this%U(i,j-1:j,k)) 
               FZ(i,j,k)=-this%rho*this%Wf(i,j,k)*sum(this%itpr_z(:,i,j,k)*this%U(i,j,k-1:k))
               ! X Viscous Flux
               FX(i,j,k)=FX(i,j,k)+this%visc_x(i,j,k)*(sum(this%divu_x(:,i,j,k)*this%U(i-1:i,j,k))+sum(this%divu_x(:,i,j,k)*this%U(i-1:i,j,k)))
               ! Y Viscous Flux
               FY(i,j,k)=FY(i,j,k)+this%visc_y(i,j,k)*(sum(this%divv_y(:,i,j,k)*this%U(i,j-1:j,k))+sum(this%itpr_y(:,i,j,k)* &
               &[0.5_WP*(sum(this%divu_x(:,i,j-1,k)*this%V(i-1:i,j-1,k))+sum(this%divu_x(:,i+1,j-1,k)*this%V(i:i+1,j-1,k))), & 
               & 0.5_WP*(sum(this%divu_x(:,i,j  ,k)*this%V(i-1:i,j  ,k))+sum(this%divu_x(:,i+1,j  ,k)*this%V(i:i+1,j  ,k)))]))
               ! Z Viscous Flux
               FZ(i,j,k)=FZ(i,j,k)+this%visc_z(i,j,k)*(sum(this%divw_z(:,i,j,k)*this%U(i,j,k-1:k))+sum(this%itpr_z(:,i,j,k)* &
               &[0.5_WP*(sum(this%divu_x(:,i,j,k-1)*this%W(i-1:i,j,k-1))+sum(this%divu_x(:,i+1,j,k-1)*this%W(i:i+1,j,k-1))), & 
               & 0.5_WP*(sum(this%divu_x(:,i,j,k  )*this%W(i-1:i,j,k  ))+sum(this%divu_x(:,i+1,j,k  )*this%W(i:i+1,j,k  )))]))
   
            end do
         end do
      end do
      ! Time derivative of rhoU
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoUdt(i,j,k)=sum(this%divp_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &              sum(this%divp_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &              sum(this%divp_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoUdt)
      
      ! Flux of rhoV
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FX(i,j,k)=-this%rho*this%Uf(i,j,k)*sum(this%itpr_x(:,i,j,k)*this%V(i-1:i,j,k))
               FY(i,j,k)=-this%rho*this%Vf(i,j,k)*sum(this%itpr_y(:,i,j,k)*this%V(i,j-1:j,k))
               FZ(i,j,k)=-this%rho*this%Wf(i,j,k)*sum(this%itpr_z(:,i,j,k)*this%V(i,j,k-1:k))
               ! X Viscous Flux
               FX(i,j,k)=FX(i,j,k)+this%visc_x(i,j,k)*(sum(this%divu_x(:,i,j,k)*this%V(i-1:i,j,k))+sum(this%itpr_x(:,i,j,k)* &
               &[0.5_WP*(sum(this%divv_y(:,i-1,j,k)*this%U(i-1,j-1:j,k))+sum(this%divv_y(:,i-1,j+1,k)*this%U(i-1,j:j+1,k))), & 
               & 0.5_WP*(sum(this%divv_y(:,i  ,j,k)*this%U(i  ,j-1:j,k))+sum(this%divv_y(:,i  ,j+1,k)*this%U(i  ,j:j+1,k)))]))
               ! Y Viscous Flux
               FY(i,j,k)=FY(i,j,k)+this%visc_y(i,j,k)*(sum(this%divv_y(:,i,j,k)*this%V(i,j-1:j,k))+sum(this%divv_y(:,i,j,k)*this%V(i,j-1:j,k)))
               ! Z Viscous Flux
               FZ(i,j,k)=FZ(i,j,k)+this%visc_z(i,j,k)*(sum(this%divw_z(:,i,j,k)*this%V(i,j,k-1:k))+sum(this%itpr_z(:,i,j,k)* &
               &[0.5_WP*(sum(this%divv_y(:,i,j,k-1)*this%W(i,j-1:j,k-1))+sum(this%divv_y(:,i,j+1,k-1)*this%W(i,j:j+1,k-1))), & 
               & 0.5_WP*(sum(this%divv_y(:,i,j,k  )*this%W(i,j-1:j,k  ))+sum(this%divv_y(:,i,j+1,k  )*this%W(i,j:j+1,k  )))]))
            end do
         end do
      end do
      ! Time derivative of rhoV
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoVdt(i,j,k)=sum(this%divp_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &              sum(this%divp_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &              sum(this%divp_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoVdt)
      
      ! Flux of rhoW
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FX(i,j,k)=-this%rho*this%Uf(i,j,k)*sum(this%itpr_x(:,i,j,k)*this%W(i-1:i,j,k))
               FY(i,j,k)=-this%rho*this%Vf(i,j,k)*sum(this%itpr_y(:,i,j,k)*this%W(i,j-1:j,k))
               FZ(i,j,k)=-this%rho*this%Wf(i,j,k)*sum(this%itpr_z(:,i,j,k)*this%W(i,j,k-1:k))
               ! X Viscous Flux
               FX(i,j,k)=FX(i,j,k)+this%visc_x(i,j,k)*(sum(this%divu_x(:,i,j,k)*this%W(i-1:i,j,k))+sum(this%itpr_x(:,i,j,k)* &
               &[0.5_WP*(sum(this%divw_z(:,i-1,j,k)*this%U(i-1,j,k-1:k))+sum(this%divw_z(:,i-1,j,k+1)*this%U(i-1,j,k:k+1))), & 
               & 0.5_WP*(sum(this%divw_z(:,i  ,j,k)*this%U(i  ,j,k-1:k))+sum(this%divw_z(:,i  ,j,k+1)*this%U(i  ,j,k:k+1)))]))
               ! Y Viscous Flux
               FY(i,j,k)=FY(i,j,k)+this%visc_y(i,j,k)*(sum(this%divv_y(:,i,j,k)*this%W(i,j-1:j,k))+sum(this%itpr_y(:,i,j,k)* &
               &[0.5_WP*(sum(this%divw_z(:,i,j-1,k)*this%V(i,j-1,k-1:k))+sum(this%divw_z(:,i,j-1,k+1)*this%V(i,j-1,k:k+1))), & 
               & 0.5_WP*(sum(this%divw_z(:,i,j  ,k)*this%V(i,j  ,k-1:k))+sum(this%divw_z(:,i,j  ,k+1)*this%V(i,j  ,k:k+1)))]))
               ! Z Viscous Flux
               FZ(i,j,k)=FZ(i,j,k)+this%visc_z(i,j,k)*(sum(this%divw_z(:,i,j,k)*this%W(i,j,k-1:k))+sum(this%divw_z(:,i,j,k)*this%W(i,j,k-1:k)))
            end do
         end do
      end do
      ! Time derivative of rhoW
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoWdt(i,j,k)=sum(this%divp_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &              sum(this%divp_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &              sum(this%divp_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoWdt)
      
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      
   end subroutine get_dmomdt
   
   
   !> Calculate the velocity divergence based on U/V/W
   subroutine get_div(this,src)
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), optional :: src !< Mass source term
      integer :: i,j,k
      ! Calculate divergence of velocity
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%div(i,j,k)=sum(this%divp_x(:,i,j,k)*this%Uf(i:i+1,j,k))+&
               &               sum(this%divp_y(:,i,j,k)*this%Vf(i,j:j+1,k))+&
               &               sum(this%divp_z(:,i,j,k)*this%Wf(i,j,k:k+1))
            end do
         end do
      end do
      ! If present, account for mass source
      if (present(src)) then
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%div(i,j,k)=this%div(i,j,k)-src(i,j,k)
               end do
            end do
         end do
      end if
      ! Sync it
      call this%cfg%sync(this%div)
   end subroutine get_div
   
   
   !> Calculate the pressure gradient based on P
   subroutine get_pgrad(this,P,Pgradx,Pgrady,Pgradz)
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: P      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgradx !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgrady !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgradz !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      Pgradx=0.0_WP; Pgrady=0.0_WP; Pgradz=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               Pgradx(i,j,k)=sum(this%divu_x(:,i,j,k)*P(i-1:i,j,k))
               Pgrady(i,j,k)=sum(this%divv_y(:,i,j,k)*P(i,j-1:j,k))
               Pgradz(i,j,k)=sum(this%divw_z(:,i,j,k)*P(i,j,k-1:k))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(Pgradx)
      call this%cfg%sync(Pgrady)
      call this%cfg%sync(Pgradz)
   end subroutine get_pgrad
   
   
   !> Get face viscosity for the purpose of sgs models
   subroutine get_viscosity(this)
      use messager,  only: die
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%visc_x(i,j,k)=sum(this%itpr_x(:,i,j,k)*this%visc(i-1:i,j,k))
               this%visc_y(i,j,k)=sum(this%itpr_y(:,i,j,k)*this%visc(i,j-1:j,k))
               this%visc_z(i,j,k)=sum(this%itpr_z(:,i,j,k)*this%visc(i,j,k-1:k))
            end do
         end do
      end do
      ! Synchronize boundaries - not really needed...
      call this%cfg%sync(this%visc)
      call this%cfg%sync(this%visc_x)
      call this%cfg%sync(this%visc_y)
      call this%cfg%sync(this%visc_z)
   end subroutine get_viscosity

   ! !> Compute face velocities from U/V/W field based on Michael's code
   subroutine update_faceU(this,U,V,W,Uface,Vface,Wface)
      use vfs_class, only : vfs
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Uface
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Vface
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Wface
      integer :: i,j,k
      real(WP) :: vol_r,vol_l,rho_r,rho_l
      logical  :: flag
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (this%umask(i,j,k).ne.2) Uface(i,j,k)=sum(this%itpr_x(:,i,j,k)*U(i-1:i,j,k))
               if (this%vmask(i,j,k).ne.2) Vface(i,j,k)=sum(this%itpr_y(:,i,j,k)*V(i,j-1:j,k))
               if (this%wmask(i,j,k).ne.2) Wface(i,j,k)=sum(this%itpr_z(:,i,j,k)*W(i,j,k-1:k))
            end do
         end do
      end do
      call this%cfg%sync(Uface)
      call this%cfg%sync(Vface)
      call this%cfg%sync(Wface)
   end subroutine update_faceU

   subroutine update_pgrad_all(this,dt)
      use vfs_class, only: vfs
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(:,:,:), allocatable :: PgradX,PgradY,PgradZ

      ! Allocate flux arrays
      allocate(PgradX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));PgradX=0.0_WP
      allocate(PgradY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));PgradY=0.0_WP
      allocate(PgradZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));PgradZ=0.0_WP

      call this%get_cell_pgrad(this%P,PgradX,PgradY,PgradZ,.false.)
      this%U=this%U-dt*PgradX/this%rho
      this%V=this%V-dt*PgradY/this%rho
      this%W=this%W-dt*PgradZ/this%rho

      call this%get_pgrad(this%P,PgradX,PgradY,PgradZ)
      
      this%Uf=this%Uf-dt*PgradX/this%rho
      this%Vf=this%Vf-dt*PgradY/this%rho
      this%Wf=this%Wf-dt*PgradZ/this%rho
      
   end subroutine update_pgrad_all

   subroutine get_cell_pgrad(this,P,Pgradx,Pgrady,Pgradz,is_correction)
      use vfs_class, only: vfs
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: P !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Pgradx !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Pgrady !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Pgradz !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      logical, intent(in) :: is_correction
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: PX,PY,PZ
      ! Allocate flux arrays
      allocate(PX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));PX=0.0_WP
      allocate(PY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));PY=0.0_WP
      allocate(PZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));PZ=0.0_WP

      ! Get density weighted face pressure
      update_faceP: block
         real(WP), dimension(0:1) :: rho_f
         real(WP) :: vol_l,vol_r
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Update face pressure and density in X
                  PX(i,j,k)=sum(this%itpr_x(:,i,j,k)*P(i-1:i,j,k))
                  ! Update face pressure and density in Y
                  PY(i,j,k)=sum(this%itpr_y(:,i,j,k)*P(i,j-1:j,k))
                  ! Update face pressure and density in Z
                  PZ(i,j,k)=sum(this%itpr_z(:,i,j,k)*P(i,j,k-1:k))
               end do
            end do
         end do
      end block update_faceP
      
      
      ! ! Apply hydrostatic pressure correction at solid walls
      ! apply_hydrostatic_pressure: block
      !    real(WP) :: hydro_inc, dist
         ! if (is_correction) then
      !       do k = this%cfg%kmin_, this%cfg%kmax_
      !          do j = this%cfg%jmin_, this%cfg%jmax_
      !             do i = this%cfg%imin_, this%cfg%imax_
      !                ! Check X-
      !                if (this%mask(i,j,k).eq.0.and.this%mask(i-1,j,k).eq.1) then
      !                   dist = -0.5_WP * this%cfg%dx(i)
      !                   hydro_inc = this%rho * this%gravity(1) * dist
      !                   Px(i,j,k) = P(i,j,k) + hydro_inc
      !                end if
      !                ! Check X+
      !                if (this%mask(i,j,k).eq.0.and.this%mask(i+1,j,k).eq.1) then
      !                   dist = +0.5_WP * this%cfg%dx(i)
      !                   hydro_inc = this%rho * this%gravity(1) * dist
      !                   Px(i+1,j,k) = P(i,j,k) + hydro_inc
      !                end if
      !                ! Check Y-
      !                if (this%mask(i,j,k).eq.0.and.this%mask(i,j-1,k).eq.1) then
      !                   dist = -0.5_WP * this%cfg%dy(j)
      !                   hydro_inc = this%rho * this%gravity(2) * dist
      !                   Py(i,j,k) = P(i,j,k) + hydro_inc
      !                end if
      !                ! Check Y+
      !                if (this%mask(i,j,k).eq.0.and.this%mask(i,j+1,k).eq.1) then
      !                   dist = +0.5_WP * this%cfg%dy(j)
      !                   hydro_inc = this%rho * this%gravity(2) * dist
      !                   Py(i,j+1,k) = P(i,j,k) + hydro_inc
      !                end if
      !                ! Check Z-
      !                if (this%mask(i,j,k).eq.0.and.this%mask(i,j,k-1).eq.1) then
      !                   dist = -0.5_WP * this%cfg%dz(k)
      !                   hydro_inc = this%rho * this%gravity(3) * dist
      !                   Pz(i,j,k) = P(i,j,k) + hydro_inc
      !                end if
      !                ! Check Z+
      !                if (this%mask(i,j,k).eq.0.and.this%mask(i,j,k+1).eq.1) then
      !                   dist = +0.5_WP * this%cfg%dz(k)
      !                   hydro_inc = this%rho * this%gravity(3) * dist
      !                   Pz(i,j,k+1) = P(i,j,k) + hydro_inc
      !                end if
      !             end do
      !          end do
      !       end do
      !       call this%cfg%sync(Px)
      !       call this%cfg%sync(Py)
      !       call this%cfg%sync(Pz)
      !       ! Apply boundary conditions on boundary faces (free-slip walls and clipped neuman bcs)
            call this%apply_bcond_facepressure(P,Px,PY,PZ)
         ! end if
      ! end block apply_hydrostatic_pressure

      

      ! Compute the pressure gradient at cell centers
      get_pgrad_cellcenter: block
         Pgradx=0.0_WP; Pgrady=0.0_WP; Pgradz=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%mask(i,j,k).gt.0) cycle
                  Pgradx(i,j,k)=sum(this%divp_x(:,i,j,k)*PX(i:i+1,j,k))
                  Pgrady(i,j,k)=sum(this%divp_y(:,i,j,k)*PY(i,j:j+1,k))
                  Pgradz(i,j,k)=sum(this%divp_z(:,i,j,k)*PZ(i,j,k:k+1))
               end do
            end do
         end do
         call this%cfg%sync(Pgradx)
         call this%cfg%sync(Pgrady)
         call this%cfg%sync(Pgradz)
      end block get_pgrad_cellcenter

   end subroutine get_cell_pgrad

   
   !> Calculate the deviatoric part of the strain rate tensor from U/V/W
   !> 1: du/dx-div/3
   !> 2: dv/dy-div/3
   !> 3: dw/dz-div/3
   !> 4: (du/dy+dv/dx)/2
   !> 5: (dv/dz+dw/dy)/2
   !> 6: (dw/dx+du/dz)/2
   subroutine get_strainrate(this,SR)
      use messager, only: die
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: SR  !< Needs to be (1:6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      real(WP) :: div
      integer :: i,j,k
      
      ! Check SR's first dimension
	   if (size(SR,dim=1).ne.6) call die('[incomp get_strainrate] SR should be of size (1:6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')
      
      ! Compute dudx, dvdy, and dwdz first
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               SR(1,i,j,k)=sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))
               SR(2,i,j,k)=sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))
               SR(3,i,j,k)=sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))
               div=sum(SR(1:3,i,j,k))/3.0_WP
               SR(1,i,j,k)=SR(1,i,j,k)-div
               SR(2,i,j,k)=SR(2,i,j,k)-div
               SR(3,i,j,k)=SR(3,i,j,k)-div
            end do
         end do
      end do
      
      ! Allocate velocity gradient components
	   allocate(dudy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dudz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Calculate components of the velocity gradient at their natural locations with an extra cell for interpolation
	   do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               dudy(i,j,k)=sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))
               dudz(i,j,k)=sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k))
               dvdx(i,j,k)=sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))
               dvdz(i,j,k)=sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k))
               dwdx(i,j,k)=sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))
               dwdy(i,j,k)=sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))
            end do
         end do
      end do
      
      ! Interpolate off-diagonal components of the velocity gradient to the cell center and store strain rate
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               SR(4,i,j,k)=0.125_WP*(sum(dudy(i:i+1,j:j+1,k    ))+sum(dvdx(i:i+1,j:j+1,k    )))
               SR(5,i,j,k)=0.125_WP*(sum(dvdz(i    ,j:j+1,k:k+1))+sum(dwdy(i    ,j:j+1,k:k+1)))
               SR(6,i,j,k)=0.125_WP*(sum(dwdx(i:i+1,j    ,k:k+1))+sum(dudz(i:i+1,j    ,k:k+1)))
            end do
         end do
      end do
      
      ! Apply a Neumann condition in non-periodic directions
	   if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            SR(:,this%cfg%imin-1,:,:)=SR(:,this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) SR(:,this%cfg%imax+1,:,:)=SR(:,this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            SR(:,:,this%cfg%jmin-1,:)=SR(:,:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) SR(:,:,this%cfg%jmax+1,:)=SR(:,:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            SR(:,:,:,this%cfg%kmin-1)=SR(:,:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) SR(:,:,:,this%cfg%kmax+1)=SR(:,:,:,this%cfg%kmax)
      end if
      
      ! Ensure zero in walls
	   do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) SR(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Sync it
	   call this%cfg%sync(SR)
      
      ! Deallocate velocity gradient storage
	   deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
      
   end subroutine get_strainrate

   
   !> Calculate the velocity gradient tensor from U/V/W
   !> Note that gradu(i,j)=duj/dxi
   subroutine get_gradu(this,gradu)
      use messager, only: die
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: gradu  !< Needs to be (1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      
      ! Check gradu's first two dimensions
	   if (size(gradu,dim=1).ne.3.or.size(gradu,dim=2).ne.3) call die('[incomp get_gradu] gradu should be of size (1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')
      
      ! Compute dudx, dvdy, and dwdz first
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               gradu(1,1,i,j,k)=sum(this%grdu_x(:,i,j,k)*this%Uf(i:i+1,j,k))
               gradu(2,2,i,j,k)=sum(this%grdv_y(:,i,j,k)*this%Vf(i,j:j+1,k))
               gradu(3,3,i,j,k)=sum(this%grdw_z(:,i,j,k)*this%Wf(i,j,k:k+1))
            end do
         end do
      end do
      
      ! Allocate velocity gradient components
	   allocate(dudy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dudz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Calculate components of the velocity gradient at their natural locations with an extra cell for interpolation
	   do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               dudy(i,j,k)=sum(this%grdu_y(:,i,j,k)*this%Uf(i,j-1:j,k))
               dudz(i,j,k)=sum(this%grdu_z(:,i,j,k)*this%Uf(i,j,k-1:k))
               dvdx(i,j,k)=sum(this%grdv_x(:,i,j,k)*this%Vf(i-1:i,j,k))
               dvdz(i,j,k)=sum(this%grdv_z(:,i,j,k)*this%Vf(i,j,k-1:k))
               dwdx(i,j,k)=sum(this%grdw_x(:,i,j,k)*this%Wf(i-1:i,j,k))
               dwdy(i,j,k)=sum(this%grdw_y(:,i,j,k)*this%Wf(i,j-1:j,k))
            end do
         end do
      end do
      
      ! Interpolate off-diagonal components of the velocity gradient to the cell center
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               gradu(2,1,i,j,k)=0.25_WP*sum(dudy(i:i+1,j:j+1,k))
               gradu(3,1,i,j,k)=0.25_WP*sum(dudz(i:i+1,j,k:k+1))
               gradu(1,2,i,j,k)=0.25_WP*sum(dvdx(i:i+1,j:j+1,k))
               gradu(3,2,i,j,k)=0.25_WP*sum(dvdz(i,j:j+1,k:k+1))
               gradu(1,3,i,j,k)=0.25_WP*sum(dwdx(i:i+1,j,k:k+1))
               gradu(2,3,i,j,k)=0.25_WP*sum(dwdy(i,j:j+1,k:k+1))
            end do
         end do
      end do
      
      ! Apply a Neumann condition in non-periodic directions
	   if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            gradu(:,:,this%cfg%imin-1,:,:)=gradu(:,:,this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) gradu(:,:,this%cfg%imax+1,:,:)=gradu(:,:,this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            gradu(:,:,:,this%cfg%jmin-1,:)=gradu(:,:,:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) gradu(:,:,:,this%cfg%jmax+1,:)=gradu(:,:,:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            gradu(:,:,:,:,this%cfg%kmin-1)=gradu(:,:,:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) gradu(:,:,:,:,this%cfg%kmax+1)=gradu(:,:,:,:,this%cfg%kmax)
      end if
      
      ! Ensure zero in walls
	   do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) gradu(:,:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Sync it
	   call this%cfg%sync(gradu)
      
      ! Deallocate velocity gradient storage
	   deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
      
   end subroutine get_gradu
   
   
   !> Calculate vorticity vector
   subroutine get_vorticity(this,vort)
      use messager, only: die
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: vort  !< Needs to be (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      
      ! Check vort's first two dimensions
      if (size(vort,dim=1).ne.3) call die('[incomp get_vorticity] vort should be of size (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')

      ! Allocate velocity gradient components
      allocate(dudy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dudz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Calculate components of the velocity gradient at their natural locations with an extra cell for interpolation
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               dudy(i,j,k)=sum(this%grdu_y(:,i,j,k)*this%Uf(i,j-1:j,k))
               dudz(i,j,k)=sum(this%grdu_z(:,i,j,k)*this%Uf(i,j,k-1:k))
               dvdx(i,j,k)=sum(this%grdv_x(:,i,j,k)*this%Vf(i-1:i,j,k))
               dvdz(i,j,k)=sum(this%grdv_z(:,i,j,k)*this%Vf(i,j,k-1:k))
               dwdx(i,j,k)=sum(this%grdw_x(:,i,j,k)*this%Wf(i-1:i,j,k))
               dwdy(i,j,k)=sum(this%grdw_y(:,i,j,k)*this%Wf(i,j-1:j,k))
            end do
         end do
      end do

      ! Interpolate off-diagonal components of the velocity gradient to the cell center
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               vort(1,i,j,k)=0.25_WP*(sum(dwdy(i,j:j+1,k:k+1))-sum(dvdz(i,j:j+1,k:k+1)))
               vort(2,i,j,k)=0.25_WP*(sum(dudz(i:i+1,j,k:k+1))-sum(dwdx(i:i+1,j,k:k+1)))
               vort(3,i,j,k)=0.25_WP*(sum(dvdx(i:i+1,j:j+1,k))-sum(dudy(i:i+1,j:j+1,k)))
            end do
         end do
      end do
      
      ! Apply a Neumann condition in non-periodic directions
	   if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            vort(:,this%cfg%imin-1,:,:)=vort(:,this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) vort(:,this%cfg%imax+1,:,:)=vort(:,this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            vort(:,:,this%cfg%jmin-1,:)=vort(:,:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) vort(:,:,this%cfg%jmax+1,:)=vort(:,:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            vort(:,:,:,this%cfg%kmin-1)=vort(:,:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) vort(:,:,:,this%cfg%kmax+1)=vort(:,:,:,this%cfg%kmax)
      end if
      
      ! Ensure zero in walls
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) vort(:,i,j,k)=0.0_WP
            end do
         end do
      end do

      ! Sync it
      call this%cfg%sync(vort)

      ! Deallocate velocity gradient storage
      deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)

   end subroutine get_vorticity
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cflc,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cflc
      real(WP), optional :: cfl
      integer :: i,j,k,ierr
      real(WP) :: my_CFLc_x,my_CFLc_y,my_CFLc_z,my_CFLv_x,my_CFLv_y,my_CFLv_z
      
      ! Set the CFLs to zero
      my_CFLc_x=0.0_WP; my_CFLc_y=0.0_WP; my_CFLc_z=0.0_WP
      my_CFLv_x=0.0_WP; my_CFLv_y=0.0_WP; my_CFLv_z=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_CFLc_x=max(my_CFLc_x,abs(this%U(i,j,k))*this%cfg%dxmi(i))
               my_CFLc_y=max(my_CFLc_y,abs(this%V(i,j,k))*this%cfg%dymi(j))
               my_CFLc_z=max(my_CFLc_z,abs(this%W(i,j,k))*this%cfg%dzmi(k))
               my_CFLv_x=max(my_CFLv_x,4.0_WP*this%visc(i,j,k)*this%cfg%dxi(i)**2/this%rho)
               my_CFLv_y=max(my_CFLv_y,4.0_WP*this%visc(i,j,k)*this%cfg%dyi(j)**2/this%rho)
               my_CFLv_z=max(my_CFLv_z,4.0_WP*this%visc(i,j,k)*this%cfg%dzi(k)**2/this%rho)
            end do
         end do
      end do
      my_CFLc_x=my_CFLc_x*dt; my_CFLc_y=my_CFLc_y*dt; my_CFLc_z=my_CFLc_z*dt
      my_CFLv_x=my_CFLv_x*dt; my_CFLv_y=my_CFLv_y*dt; my_CFLv_z=my_CFLv_z*dt
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFLc_x,this%CFLc_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLc_y,this%CFLc_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLc_z,this%CFLc_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_x,this%CFLv_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_y,this%CFLv_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_z,this%CFLv_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
      ! Return the maximum convective CFL
      cflc=max(this%CFLc_x,this%CFLc_y,this%CFLc_z)
      
      ! If asked for, also return the maximum overall CFL
      if (present(CFL)) cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
      
   end subroutine get_cfl
   
   
   !> Calculate the max of our fields
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k,ierr
      real(WP) :: my_Umax,my_Vmax,my_Wmax,my_Pmax,my_divmax
      
      ! Set all to zero
      my_Umax=0.0_WP; my_Vmax=0.0_WP; my_Wmax=0.0_WP; my_Pmax=0.0_WP; my_divmax=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_Umax  =max(my_Umax  ,abs(this%U(i,j,k)  ))
               my_Vmax  =max(my_Vmax  ,abs(this%V(i,j,k)  ))
               my_Wmax  =max(my_Wmax  ,abs(this%W(i,j,k)  ))
               if (this%cfg%VF(i,j,k).gt.0.0_WP) my_Pmax  =max(my_Pmax  ,abs(this%P(i,j,k)  ))
               if (this%cfg%VF(i,j,k).gt.0.0_WP) my_divmax=max(my_divmax,abs(this%div(i,j,k)))
            end do
         end do
      end do
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_Umax  ,this%Umax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Vmax  ,this%Vmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Wmax  ,this%Wmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Pmax  ,this%Pmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_divmax,this%divmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_max
   
   
   !> Compute MFR through all bcs
   subroutine get_mfr(this)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k,n,ibc,ierr
      type(bcond), pointer :: my_bc
      real(WP), dimension(:), allocatable :: my_mfr,my_area
      real(WP), dimension(:), allocatable :: canCorrect
      
      ! Ensure this%mfr is of proper size
      if (.not.allocated(this%mfr)) then
         allocate(this%mfr(this%nbc))
      else
         if (size(this%mfr).ne.this%nbc) then
            deallocate(this%mfr); allocate(this%mfr(this%nbc))
         end if
      end if
      
      ! Ensure this%area is of proper size
      if (.not.allocated(this%area)) then
         allocate(this%area(this%nbc))
      else
         if (size(this%area).ne.this%nbc) then
            deallocate(this%area); allocate(this%area(this%nbc))
         end if
      end if
      
      ! Allocate temp array for communication
      allocate(my_mfr(this%nbc))
      allocate(my_area(this%nbc))
      allocate(canCorrect(this%nbc))
      
      ! Traverse bcond list and integrate local outgoing MFR
      my_bc=>this%first_bc; ibc=1
      do while (associated(my_bc))
         
         ! Set zero local MFR and area
         my_mfr(ibc)=0.0_WP
         my_area(ibc)=0.0_WP
         if (my_bc%canCorrect) then
            canCorrect(ibc)=1.0_WP
         else
            canCorrect(ibc)=0.0_WP
         end if
         
         ! Only processes inside the bcond have a non-zero MFR
         if (my_bc%itr%amIn) then
            
            ! Implement based on bcond face and dir, loop over interior only
            select case (my_bc%face)
            case ('x')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  my_mfr(ibc)=my_mfr(ibc)+my_bc%rdir*this%Uf(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
                  my_area(ibc)=my_area(ibc)+this%cfg%dy(j)*this%cfg%dz(k)
               end do
            case ('y')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  my_mfr(ibc)=my_mfr(ibc)+my_bc%rdir*this%Vf(i,j,k)*this%cfg%dz(k)*this%cfg%dx(i)
                  my_area(ibc)=my_area(ibc)+this%cfg%dz(k)*this%cfg%dx(i)
               end do
            case ('z')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  my_mfr(ibc)=my_mfr(ibc)+my_bc%rdir*this%Wf(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j)
                  my_area(ibc)=my_area(ibc)+this%cfg%dx(i)*this%cfg%dy(j)
               end do
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next; ibc=ibc+1
         
      end do
      
      ! Sum up all values
      call MPI_ALLREDUCE(my_mfr ,this%mfr ,this%nbc,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_area,this%area,this%nbc,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      
      ! Compute the correctable area
      this%correctable_area=sum(this%area*canCorrect)
      
      ! Deallocate temp array
      deallocate(my_mfr,my_area,canCorrect)
      
   end subroutine get_mfr
   
   
   !> Correct MFR through correctable bconds
   subroutine correct_mfr(this,src)
      use mpi_f08, only: MPI_SUM
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), optional :: src !< Mass source term
      real(WP) :: mfr_error,vel_correction,int
      integer :: i,j,k,n,ii,jj,kk
      type(bcond), pointer :: my_bc
      
      ! Evaluate MFR mismatch and velocity correction
      call this%get_mfr()
      mfr_error=sum(this%mfr)
      if (present(src)) then
         ! Also account for provided source term
         call this%cfg%integrate_without_VF(src,int)
         mfr_error=mfr_error-int
      end if
      if (abs(mfr_error).lt.10.0_WP*epsilon(1.0_WP).or.abs(this%correctable_area).lt.10.0_WP*epsilon(1.0_WP)) return
      vel_correction=-mfr_error/(this%correctable_area)
      
      ! Traverse bcond list and correct bcond MFR
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside correctable bcond need to work
         if (my_bc%itr%amIn.and.my_bc%canCorrect) then
            
            ! Implement based on bcond direction, loop over all cell
            select case (my_bc%face)
            case ('x')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  ! ii = i - max(0, my_bc%dir)
                  ! this%U(ii,j,k)=this%U(ii,j,k)+my_bc%rdir*vel_correction
                  this%Uf(i,j,k)=this%Uf(i,j,k)+my_bc%rdir*vel_correction
               end do
            case ('y')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  ! jj = j - max(0, my_bc%dir)
                  ! this%V(i,jj,k)=this%V(i,jj,k)+my_bc%rdir*vel_correction
                  this%Vf(i,j,k)=this%Vf(i,j,k)+my_bc%rdir*vel_correction
               end do
            case ('z')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  ! kk = k - max(0, my_bc%dir)
                  ! this%W(i,j,kk)=this%W(i,j,kk)+my_bc%rdir*vel_correction
                  this%Wf(i,j,k)=this%Wf(i,j,k)+my_bc%rdir*vel_correction
               end do
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
      ! Sync full fields
      call this%cfg%sync(this%U);call this%cfg%sync(this%Uf)
      call this%cfg%sync(this%V);call this%cfg%sync(this%Vf)
      call this%cfg%sync(this%W);call this%cfg%sync(this%Wf)
      
   end subroutine correct_mfr
   
   
   !> Shift pressure to ensure zero average
   subroutine shift_p(this,pressure)
      implicit none
      class(incomp), intent(in) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: pressure !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: pressure_tot
      integer :: i,j,k
      
      ! Compute volume-averaged pressure
      call this%cfg%integrate(A=pressure,integral=pressure_tot); pressure_tot=pressure_tot/this%cfg%fluid_vol
      
      ! Shift the pressure
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%cfg%VF(i,j,k).gt.0.0_WP) pressure(i,j,k)=pressure(i,j,k)-pressure_tot
            end do
         end do
      end do
      call this%cfg%sync(pressure)
      
   end subroutine shift_p
   
   
   !> Solve for implicit velocity residual
   subroutine solve_implicit(this,dt,resU,resV,resW)
      implicit none
      class(incomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP) :: rhoUp,rhoUm,rhoVp,rhoVm,rhoWp,rhoWm
      
      ! If no implicit solver available, just divide by density and return
      if (.not.associated(this%implicit)) then
         resU=resU/this%rho
         resV=resV/this%rho
         resW=resW/this%rho
         call this%cfg%sync(resU)
         call this%cfg%sync(resV)
         call this%cfg%sync(resW)
         return
      end if
      ! X direction
      this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=this%rho*this%Uf(i+1,j,k); rhoVp=this%rho*this%Vf(i,j+1,k); rhoWp=this%rho*this%Wf(i,j,k+1)
               rhoUm=this%rho*this%Uf(i  ,j,k); rhoVm=this%rho*this%Vf(i,j  ,k); rhoWm=this%rho*this%Wf(i,j,k  )
               ! +X face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_x(+1,i,j,k)*this%itpr_x(-1,i+1,j,k)*rhoUp*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+dt*this%divp_x(+1,i,j,k)*this%itpr_x( 0,i+1,j,k)*rhoUp*0.5_WP
               ! -X face   
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_x( 0,i,j,k)*this%itpr_x( 0,i  ,j,k)*rhoUm*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+dt*this%divp_x( 0,i,j,k)*this%itpr_x(-1,i  ,j,k)*rhoUm*0.5_WP
               ! +Y face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_y(+1,i,j,k)*this%itpr_y(-1,i,j+1,k)*rhoVp*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+dt*this%divp_y(+1,i,j,k)*this%itpr_y( 0,i,j+1,k)*rhoVp*0.5_WP
               ! -Y face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_y( 0,i,j,k)*this%itpr_y( 0,i,j  ,k)*rhoVm*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+dt*this%divp_y( 0,i,j,k)*this%itpr_y(-1,i,j  ,k)*rhoVm*0.5_WP
               ! +Z face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_z(+1,i,j,k)*this%itpr_z(-1,i,j,k+1)*rhoWp*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+dt*this%divp_z(+1,i,j,k)*this%itpr_z( 0,i,j,k+1)*rhoWp*0.5_WP
               ! -Z face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_z( 0,i,j,k)*this%itpr_z( 0,i,j,k  )*rhoWm*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+dt*this%divp_z( 0,i,j,k)*this%itpr_z(-1,i,j,k  )*rhoWm*0.5_WP
            end do
         end do
      end do
      ! Viscosity
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divp_x(1,i,j,k)*2.0_WP*this%visc_x(i+1,j,k)*this%divu_x(-1,i+1,j,k)+&
               &                                                                this%divp_x(0,i,j,k)*2.0_WP*this%visc_x(i  ,j,k)*this%divu_x( 0,i  ,j,k)+&
               &                                                                this%divp_y(1,i,j,k)*       this%visc_y(i,j+1,k)*this%divv_y(-1,i,j+1,k)+&
               &                                                                this%divp_y(0,i,j,k)*       this%visc_y(i,j  ,k)*this%divv_y( 0,i,j  ,k)+&
               &                                                                this%divp_z(1,i,j,k)*       this%visc_z(i,j,k+1)*this%divw_z(-1,i,j,k+1)+&
               &                                                                this%divp_z(0,i,j,k)*       this%visc_z(i,j,k  )*this%divw_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divp_x(1,i,j,k)*2.0_WP*this%visc_x(i+1,j,k)*this%divu_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divp_x(0,i,j,k)*2.0_WP*this%visc_x(i  ,j,k)*this%divu_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divp_y(1,i,j,k)*       this%visc_y(i,j+1,k)*this%divv_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divp_y(0,i,j,k)*       this%visc_y(i,j  ,k)*this%divv_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divp_z(1,i,j,k)*       this%visc_z(i,j,k+1)*this%divw_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divp_z(0,i,j,k)*       this%visc_z(i,j,k  )*this%divw_z(-1,i,j,k  ))
            end do
         end do
      end do
      ! Solve implicit U problem
      call this%implicit%setup()
      this%implicit%rhs=resU
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resU=this%implicit%sol

      ! Y direction
      this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=this%rho*this%Uf(i+1,j,k); rhoVp=this%rho*this%Vf(i,j+1,k); rhoWp=this%rho*this%Wf(i,j,k+1)
               rhoUm=this%rho*this%Uf(i  ,j,k); rhoVm=this%rho*this%Vf(i,j  ,k); rhoWm=this%rho*this%Wf(i,j,k  )
               ! +X face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_x(+1,i,j,k)*this%itpr_x(-1,i+1,j,k)*rhoUp*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+dt*this%divp_x(+1,i,j,k)*this%itpr_x( 0,i+1,j,k)*rhoUp*0.5_WP
               ! -X face   
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_x( 0,i,j,k)*this%itpr_x( 0,i  ,j,k)*rhoUm*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+dt*this%divp_x( 0,i,j,k)*this%itpr_x(-1,i  ,j,k)*rhoUm*0.5_WP
               ! +Y face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_y(+1,i,j,k)*this%itpr_y(-1,i,j+1,k)*rhoVp*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+dt*this%divp_y(+1,i,j,k)*this%itpr_y( 0,i,j+1,k)*rhoVp*0.5_WP
               ! -Y face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_y( 0,i,j,k)*this%itpr_y( 0,i,j  ,k)*rhoVm*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+dt*this%divp_y( 0,i,j,k)*this%itpr_y(-1,i,j  ,k)*rhoVm*0.5_WP
               ! +Z face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_z(+1,i,j,k)*this%itpr_z(-1,i,j,k+1)*rhoWp*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+dt*this%divp_z(+1,i,j,k)*this%itpr_z( 0,i,j,k+1)*rhoWp*0.5_WP
               ! -Z face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_z( 0,i,j,k)*this%itpr_z( 0,i,j,k  )*rhoWm*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+dt*this%divp_z( 0,i,j,k)*this%itpr_z(-1,i,j,k  )*rhoWm*0.5_WP
            end do
         end do
      end do
      ! Viscosity
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divp_x(1,i,j,k)*       this%visc_x(i+1,j,k)*this%divu_x(-1,i+1,j,k)+&
               &                                                                this%divp_x(0,i,j,k)*       this%visc_x(i  ,j,k)*this%divu_x( 0,i  ,j,k)+&
               &                                                                this%divp_y(1,i,j,k)*2.0_WP*this%visc_y(i,j+1,k)*this%divv_y(-1,i,j+1,k)+&
               &                                                                this%divp_y(0,i,j,k)*2.0_WP*this%visc_y(i,j  ,k)*this%divv_y( 0,i,j  ,k)+&
               &                                                                this%divp_z(1,i,j,k)*       this%visc_z(i,j,k+1)*this%divw_z(-1,i,j,k+1)+&
               &                                                                this%divp_z(0,i,j,k)*       this%visc_z(i,j,k  )*this%divw_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divp_x(1,i,j,k)*       this%visc_x(i+1,j,k)*this%divu_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divp_x(0,i,j,k)*       this%visc_x(i  ,j,k)*this%divu_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divp_y(1,i,j,k)*2.0_WP*this%visc_y(i,j+1,k)*this%divv_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divp_y(0,i,j,k)*2.0_WP*this%visc_y(i,j  ,k)*this%divv_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divp_z(1,i,j,k)*       this%visc_z(i,j,k+1)*this%divw_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divp_z(0,i,j,k)*       this%visc_z(i,j,k  )*this%divw_z(-1,i,j,k  ))
            end do
         end do
      end do

      ! Solve implicit V problem
      call this%implicit%setup()
      this%implicit%rhs=resV
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resV=this%implicit%sol
      ! Y direction
      this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=this%rho*this%Uf(i+1,j,k); rhoVp=this%rho*this%Vf(i,j+1,k); rhoWp=this%rho*this%Wf(i,j,k+1)
               rhoUm=this%rho*this%Uf(i  ,j,k); rhoVm=this%rho*this%Vf(i,j  ,k); rhoWm=this%rho*this%Wf(i,j,k  )
               ! +X face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_x(+1,i,j,k)*this%itpr_x(-1,i+1,j,k)*rhoUp*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+dt*this%divp_x(+1,i,j,k)*this%itpr_x( 0,i+1,j,k)*rhoUp*0.5_WP
               ! -X face   
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_x( 0,i,j,k)*this%itpr_x( 0,i  ,j,k)*rhoUm*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+dt*this%divp_x( 0,i,j,k)*this%itpr_x(-1,i  ,j,k)*rhoUm*0.5_WP
               ! +Y face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_y(+1,i,j,k)*this%itpr_y(-1,i,j+1,k)*rhoVp*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+dt*this%divp_y(+1,i,j,k)*this%itpr_y( 0,i,j+1,k)*rhoVp*0.5_WP
               ! -Y face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_y( 0,i,j,k)*this%itpr_y( 0,i,j  ,k)*rhoVm*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+dt*this%divp_y( 0,i,j,k)*this%itpr_y(-1,i,j  ,k)*rhoVm*0.5_WP
               ! +Z face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_z(+1,i,j,k)*this%itpr_z(-1,i,j,k+1)*rhoWp*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+dt*this%divp_z(+1,i,j,k)*this%itpr_z( 0,i,j,k+1)*rhoWp*0.5_WP
               ! -Z face
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*this%divp_z( 0,i,j,k)*this%itpr_z( 0,i,j,k  )*rhoWm*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+dt*this%divp_z( 0,i,j,k)*this%itpr_z(-1,i,j,k  )*rhoWm*0.5_WP
            end do
         end do
      end do
      ! Viscosity
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divp_x(1,i,j,k)*       this%visc_x(i+1,j,k)*this%divu_x(-1,i+1,j,k)+&
               &                                                                this%divp_x(0,i,j,k)*       this%visc_x(i  ,j,k)*this%divu_x( 0,i  ,j,k)+&
               &                                                                this%divp_y(1,i,j,k)*       this%visc_y(i,j+1,k)*this%divv_y(-1,i,j+1,k)+&
               &                                                                this%divp_y(0,i,j,k)*       this%visc_y(i,j  ,k)*this%divv_y( 0,i,j  ,k)+&
               &                                                                this%divp_z(1,i,j,k)*2.0_WP*this%visc_z(i,j,k+1)*this%divw_z(-1,i,j,k+1)+&
               &                                                                this%divp_z(0,i,j,k)*2.0_WP*this%visc_z(i,j,k  )*this%divw_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divp_x(1,i,j,k)*       this%visc_x(i+1,j,k)*this%divu_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divp_x(0,i,j,k)*       this%visc_x(i  ,j,k)*this%divu_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divp_y(1,i,j,k)*       this%visc_y(i,j+1,k)*this%divv_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divp_y(0,i,j,k)*       this%visc_y(i,j  ,k)*this%divv_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divp_z(1,i,j,k)*2.0_WP*this%visc_z(i,j,k+1)*this%divw_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divp_z(0,i,j,k)*2.0_WP*this%visc_z(i,j,k  )*this%divw_z(-1,i,j,k  ))
            end do
         end do
      end do

      ! Solve implicit W problem
      call this%implicit%setup()
      this%implicit%rhs=resW
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resW=this%implicit%sol
      
   end subroutine solve_implicit
   
   
   !> Print out info for incompressible flow solver
   subroutine incomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(incomp), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Incompressible solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >   density = ",es12.5)') this%rho
      end if
      
   end subroutine incomp_print
   
   
end module incomp_class

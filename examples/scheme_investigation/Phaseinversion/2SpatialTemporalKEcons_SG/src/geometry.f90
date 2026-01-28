!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,n
         real(WP) :: L
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('H',L); call param_read('nx',n)
         allocate(x(n+1)); allocate(y(n+1)); allocate(z(n+1))
         
         ! Create simple rectilinear grid
         do i=1,n+1
            x(i)=real(i-1,WP)/real(n,WP)*L
         end do
         do j=1,n+1
            y(j)=real(j-1,WP)/real(n,WP)*L
         end do
         do k=1,n+1
            z(k)=real(k-1,WP)/real(n,WP)*L
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='Phaseinversion')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Create masks for this config
      create_walls: block
         ! Put walls all around
         cfg%VF=1.0_WP
         ! cfg%VF=0.0_WP
         ! cfg%VF(cfg%imin_:cfg%imax_,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_)=1.0_WP
         ! call cfg%sync(cfg%VF)
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry

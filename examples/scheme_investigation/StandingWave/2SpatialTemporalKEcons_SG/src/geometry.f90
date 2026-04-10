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
      use param,       only: param_read,param_exists
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use mathtools, only : twoPi
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         
         Lx = twoPi; Ly = twoPi; Lz = twoPi;
         ! Read in grid definition: 2D problem with different orientations
         call param_read('nx',nx); allocate(x(nx+1))
         if (param_exists('ny')) then
           call param_read('ny',ny); allocate(y(ny+1))
         else
           if (.not.param_exists('nz')) then
             ny = nx
           else
             ny = 1
           end if
           allocate(y(ny+1))
         end if
         if (param_exists('nz')) then
           call param_read('nz',nz); allocate(z(nz+1))
         else
           nz = 1; allocate(z(nz+1))
         end if
         ! Reduce length of unused dimension
         if (nx.eq.1) then
           Lx = Ly/ny
         elseif (ny.eq.1) then
           Ly = Lx/nx
         elseif (nz.eq.1) then
           Lz = Lx/nx
         end if
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='Randomsolenoidal')
         
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
         use mathtools, only: twoPi
         integer :: i,j,k
         cfg%VF=1.0_WP
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  if (cfg%ym(j).lt.0.0_WP) then
                     cfg%VF(i,j,k)=0.0_WP
                  end if
               end do
            end do
         end do
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
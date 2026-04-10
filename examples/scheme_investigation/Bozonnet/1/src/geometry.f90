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
   !> Rayleigh-Taylor instability: domain [−0.5,0.5] × [−2,2]
   !> x-periodic, y-walled, z-periodic (pseudo-2D, nz=1)
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         use mathtools, only: twoPi
         integer :: i,j,k,nx,ny,nz,n_res
         real(WP) :: Lx,Ly,Lz
         real(WP) :: dg,Hg,rho_g,Ug,mu_g,dx
         real(WP), dimension(:), allocatable :: x,y,z
         
         call param_read('Gas height', Hg);  Hg=min(0.01,Hg)
         call param_read('Gas velocity', Ug)
         call param_read('Gas density', rho_g) 
         call param_read('Gas dynamic viscosity', mu_g)
         dg = 6.0_WP*Hg/sqrt(rho_g*Ug*Hg/mu_g)
         call param_read('n', n_res)
         Lx=240*dg; Ly=Lx*1.3
         nx=240*n_res; ny=13*24*n_res
         allocate(x(nx+1)); allocate(y(ny+1))
         nz = 1; Lz = Lx / real(nx,WP); allocate(z(nz+1))

         ! Create simple rectilinear grid with origin at (0,0,0)
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! x-periodic, y-walled (slip), z-periodic (pseudo-2D)
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.true.,name='AWML')
         
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
      ! Bottom wall is no slip
      create_walls: block
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

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
         integer :: i,j,k,nx,ny,nz,ny_min,nx_ideal,ny_ideal
         real(WP) :: Lx,Ly,Lz,dg,dl,alpha_dim,alpha_nd,dx,dx_target,Ly_min
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read inputs
         ! call param_read('nx',nx)
         call param_read('Delta g',dg)
         call param_read('Delta l',dl)
         call param_read('Wavenumber',alpha_nd)  ! non-dimensional: α·δ_g
         
         ! Compute dimensional wavenumber and wavelength
         alpha_dim = alpha_nd / dg
         Lx = twoPi / alpha_dim

         ! Establish the target dx based on the boundary layer requirement
         dx_target = dg / 16.0_WP
         nx_ideal = Lx / dx_target

         ! Enforce a minimum of 128 cells across Lx.
         ! Note: Switched from 'nint' to 'ceiling' to guarantee dx <= dx_target
         nx = max(128, ceiling(nx_ideal / 16.0_WP) * 16)
         dx = Lx / real(nx, WP)

         ! Compute ny, maintaining a uniform grid (dy = dx)
         Ly_min = 12.0_WP * dg
         ny_ideal = Ly_min / dx
         ny = max(16, ceiling(ny_ideal / 16.0_WP) * 16)
         Ly = real(ny, WP) * dx



         ! ! Compute ny: ensure Ly >= 12·δ_g, round up to next even number
         ! dx = Lx / real(nx, WP)
         ! ny_min = ceiling(12.0_WP * dg / dx)
         ! if (mod(ny_min, 2) /= 0) ny_min = ny_min + 1
         ! ny = ny_min
         ! Ly = real(ny, WP) * dx

         ! Ly = Lx; ny=nx

         ! z-direction: pseudo-2D (1 cell)
         nz = 1
         Lz = dx
         
         ! ! Print grid info
         ! print '(A)',        '=== Grid Setup ==='
         ! print '(A,I6)',     '  nx       = ', nx
         ! print '(A,I6)',     '  ny       = ', ny
         ! print '(A,ES12.4)', '  Lx       = ', Lx
         ! print '(A,ES12.4)', '  Ly       = ', Ly
         ! print '(A,ES12.4)', '  dx = dy  = ', dx
         ! print '(A,ES12.4)', '  Ly/dg    = ', Ly/dg
         ! print '(A,ES12.4)', '  alpha_nd = ', alpha_nd
         ! print '(A)',        '==================='
         
         ! Build grid arrays centered at origin
         allocate(x(nx+1), y(ny+1), z(nz+1))
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! x-periodic, y-walled (slip), z-periodic (pseudo-2D)
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='MPKH')
         
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
      
      
      ! Create masks for this config - walls top and bottom
      create_walls: block
         ! cfg%VF=0.0_WP; cfg%VF(cfg%imin_:cfg%imax_,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_)=1.0_WP; call cfg%sync(cfg%VF)
         cfg%VF=1.0_WP
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry

!> Test VOF advection with time loop
module mod_test_amrvof
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrvof_class,      only: amrvof
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   use string,            only: itoa, rtoa
   use amrex_amr_module,  only: amrex_multifab, amrex_multifab_destroy
   implicit none
   private
   public :: test_amrvof

   ! Grid
   type(amrgrid), target :: amr

   ! VOF solver
   type(timetracker) :: time
   type(amrvof), target :: vof

   ! Velocity MultiFabs (finest level only, staggered)
   type(amrex_multifab) :: U, V, W
   integer :: vel_ng = 2  ! Ghost cells for velocity

   ! Sphere parameters
   real(WP) :: sphere_xc, sphere_yc, sphere_zc, sphere_radius

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regrid parameters (disabled by default)
   type(event) :: regrid_evt

   ! Monitoring
   type(monitor) :: mfile

contains

   !> Initialize VF field with sphere using levelset-based moments
   subroutine sphere_init(solver, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_distromap, &
      &                           amrex_mfiter_build, amrex_mfiter_destroy
      use mms_geom,         only: initialize_volume_moments
      use amrvof_geometry,  only: VFlo
      class(amrvof), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pCliq, pCgas
      real(WP) :: dx, dy, dz
      real(WP), dimension(3) :: BL, BG
      integer :: i, j, k
      
      dx = solver%amr%dx(lvl)
      dy = solver%amr%dy(lvl)
      dz = solver%amr%dz(lvl)
      
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi, ba, dm, tiling=.false.)
      do while (mfi%next())
         bx = mfi%growntilebox(solver%VF%ng)  ! Include ghost cells
         pVF => solver%VF%mf(lvl)%dataptr(mfi)
         pCliq => solver%Cliq%mf(lvl)%dataptr(mfi)
         pCgas => solver%Cgas%mf(lvl)%dataptr(mfi)
         
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  ! Compute VF and barycenters from levelset with 4 levels of refinement
                  call initialize_volume_moments(lo=[solver%amr%xlo + real(i  ,WP)*dx, solver%amr%ylo + real(j  ,WP)*dy, solver%amr%zlo + real(k  ,WP)*dz], &
                  &                              hi=[solver%amr%xlo + real(i+1,WP)*dx, solver%amr%ylo + real(j+1,WP)*dy, solver%amr%zlo + real(k+1,WP)*dz], &
                  &                              levelset=sphere_levelset, time=time, level=4, VFlo=VFlo, VF=pVF(i,j,k,1), BL=BL, BG=BG)
                  ! Store barycenters
                  pCliq(i,j,k,1:3) = BL
                  pCgas(i,j,k,1:3) = BG
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine sphere_init
   
   !> Sphere levelset function with periodic distance
   function sphere_levelset(xyz, t) result(G)
      implicit none
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP), dimension(3) :: d, L
      L = [amr%xhi - amr%xlo, amr%yhi - amr%ylo, amr%zhi - amr%zlo]
      d = xyz - [sphere_xc, sphere_yc, sphere_zc]
      d = d - L * nint(d / L)  ! Nearest image
      G = sphere_radius - sqrt(sum(d**2))
   end function sphere_levelset


   !> Main test routine
   subroutine test_amrvof()
      implicit none

      call log("=== VOF Advection Test ===")

      ! Create amrgrid
      create_amrgrid: block
         use param, only: param_read
         amr%name = 'vof_advect'
         call param_read('Base nx', amr%nx, default=32)
         call param_read('Base ny', amr%ny, default=32)
         call param_read('Base nz', amr%nz, default=32)
         amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
         amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
         amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
         amr%xper = .true.; amr%yper = .true.; amr%zper = .true.
         call param_read('Max level', amr%maxlvl, default=0)
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         use param, only: param_read
         time = timetracker(amRoot=amr%amRoot)
         call param_read('Max time', time%tmax, default=1.0_WP)
         call param_read('Max dt', time%dtmax, default=0.01_WP)
         call param_read('Max CFL', time%cflmax, default=0.5_WP)
         time%dt = time%dtmax
         time%itmax = 1  ! No sub-iterations
      end block initialize_timetracker

      ! Setup sphere parameters
      setup_sphere: block
         use param, only: param_read
         call param_read('Sphere radius', sphere_radius, default=0.25_WP)
         sphere_xc = 0.5_WP * (amr%xlo + amr%xhi)
         sphere_yc = 0.5_WP * (amr%ylo + amr%yhi)
         sphere_zc = 0.5_WP * (amr%zlo + amr%zhi)
         call log("  Sphere center: ("//trim(rtoa(sphere_xc))//", "// &
         &        trim(rtoa(sphere_yc))//", "//trim(rtoa(sphere_zc))//")")
         call log("  Sphere radius: "//trim(rtoa(sphere_radius)))
      end block setup_sphere

      ! Create VOF solver
      create_vof_solver: block
         vof%user_init => sphere_init
         call vof%initialize(amr, name='sphere_vof')
      end block create_vof_solver

      ! Initialize grid
      initialize: block
         call amr%init_from_scratch(time=time%t)
      end block initialize

      ! Build initial PLIC
      build_initial_plic: block
         call vof%build_plic(time%t)
         call log("  Initial PLIC constructed")
      end block build_initial_plic

      ! Create velocity MultiFabs at finest level (staggered)
      create_velocity: block
         integer :: lvl
         
         lvl = amr%clvl()
         
         ! Build staggered MultiFabs: U (x-face), V (y-face), W (z-face)
         call amr%mfab_build(lvl, U, ncomp=1, nover=vel_ng, atface=[.true., .false., .false.])
         call amr%mfab_build(lvl, V, ncomp=1, nover=vel_ng, atface=[.false., .true., .false.])
         call amr%mfab_build(lvl, W, ncomp=1, nover=vel_ng, atface=[.false., .false., .true.])
         
         ! Initialize: uniform translation U=1, V=0, W=0
         call U%setval(1.0_WP)
         call V%setval(0.0_WP)
         call W%setval(0.0_WP)
         
         call log("  Velocity field initialized: U=1, V=W=0")
      end block create_velocity

      ! Create visualization
      create_visualization: block
         use param, only: param_read
         call viz%initialize(amr, 'vof_advect')
         call viz%add_scalar(vof%VF, 1, 'VF')
         call viz%add_surfmesh(vof%smesh, 'plic')
         ! Create visualization output event
         viz_evt = event(time=time, name='Visualization output')
         call param_read('Output period', viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_visualization

      ! Regridding setup (disabled by default)
      regrid_setup: block
         use param, only: param_read
         regrid_evt = event(time=time, name='Regrid')
         call param_read('Regrid nsteps', regrid_evt%nper, default=0)
      end block regrid_setup

      ! Create monitor
      create_monitor: block
         call vof%get_info()
         mfile = monitor(amRoot=amr%amRoot, name='vof_advect')
         call mfile%add_column(time%n, 'Timestep')
         call mfile%add_column(time%t, 'Time')
         call mfile%add_column(time%dt, 'dt')
         call mfile%add_column(vof%VFint, 'VFint')
         call mfile%add_column(vof%VFmin, 'VFmin')
         call mfile%add_column(vof%VFmax, 'VFmax')
         call mfile%write()
      end block create_monitor

      ! Time integration loop
      time_loop: do while (.not.time%done())

         ! Compute CFL and update dt based on CFL constraint
         call vof%get_cfl(U=U, V=V, W=W, dt=time%dt, cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Store old state
         call vof%VFold%copy(src=vof%VF)
         call vof%Cliqold%copy(src=vof%Cliq)
         call vof%Cgasold%copy(src=vof%Cgas)
         call vof%PLICold%copy(src=vof%PLIC)

         ! Advect VOF
         call vof%advance_vof(U=U, V=V, W=W, dt=time%dt, time=time%t)

         ! Rebuild PLIC and reset moments for consistency
         call vof%build_plic(time%t)
         call vof%reset_moments()

         ! Regrid if event triggers (disabled by default)
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0, time=time%t)
         end if

         ! Monitor output
         call vof%get_info()
         call mfile%write()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)
         
      end do time_loop

      ! Final summary
      final_summary: block
         call log("  Final VFint: "//trim(rtoa(vof%VFint)))
         call log("  Final VFmin: "//trim(rtoa(vof%VFmin)))
         call log("  Final VFmax: "//trim(rtoa(vof%VFmax)))
      end block final_summary

      ! Cleanup
      cleanup: block
         call amrex_multifab_destroy(U)
         call amrex_multifab_destroy(V)
         call amrex_multifab_destroy(W)
         call viz%finalize()
         call vof%finalize()
         call amr%finalize()
         call time%finalize()
         call mfile%finalize()
         call viz_evt%finalize()
         call regrid_evt%finalize()
      end block cleanup

      call log("=== VOF Advection Test Complete ===")

   end subroutine test_amrvof

end module mod_test_amrvof

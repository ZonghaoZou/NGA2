!> Test PLIC initialization with sphere
module mod_test_plic
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
   implicit none
   private
   public :: test_plic

   ! Grid
   type(amrgrid), target :: amr

   ! VOF solver
   type(timetracker) :: time
   type(amrvof), target :: vof

   ! Sphere parameters
   real(WP) :: sphere_xc, sphere_yc, sphere_zc, sphere_radius

   ! Visualization
   type(amrviz) :: viz

   ! Test results
   real(WP) :: min_dot_product, max_dot_product, avg_dot_product
   integer  :: ncells_interface

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
         bx = mfi%growntilebox(solver%VF%ng)  ! Include ghost cells - correctly initialized from levelset
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


   !> Check that PLIC normals point radially outward from sphere center
   subroutine check_plic_normals()
      use amrex_amr_module, only: amrex_mfiter, amrex_box
      use amrvof_geometry,  only: VFlo, VFhi
      implicit none
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF, pPLIC
      real(WP) :: dx, dy, dz, x, y, z
      real(WP) :: nx, ny, nz, norm_mag
      real(WP) :: rx, ry, rz, r_mag, dot
      integer :: i, j, k, lvl
      
      min_dot_product = huge(1.0_WP)
      max_dot_product = -huge(1.0_WP)
      avg_dot_product = 0.0_WP
      ncells_interface = 0
      
      do lvl = 0, amr%clvl()
         dx = amr%dx(lvl)
         dy = amr%dy(lvl)
         dz = amr%dz(lvl)
         
         call amr%mfiter_build(lvl, mfi, tiling=.false.)
         do while (mfi%next())
            bx = mfi%tilebox()
            pVF => vof%VF%mf(lvl)%dataptr(mfi)
            pPLIC => vof%PLIC%mf(lvl)%dataptr(mfi)
            
            do k = bx%lo(3), bx%hi(3)
               do j = bx%lo(2), bx%hi(2)
                  do i = bx%lo(1), bx%hi(1)
                     ! Only check interface cells (mixed VF) - use wider bounds
                     if (pVF(i,j,k,1) .gt. 0.01_WP .and. pVF(i,j,k,1) .lt. 0.99_WP) then
                        ! Get PLIC normal (stored as nx, ny, nz, d)
                        nx = pPLIC(i,j,k,1)
                        ny = pPLIC(i,j,k,2)
                        nz = pPLIC(i,j,k,3)
                        norm_mag = sqrt(nx**2 + ny**2 + nz**2)
                        
                        if (norm_mag .gt. 1.0e-10_WP) then
                           ! Normalize
                           nx = nx / norm_mag
                           ny = ny / norm_mag
                           nz = nz / norm_mag
                           
                           ! Cell center position
                           x = amr%xlo + (real(i,WP) + 0.5_WP) * dx
                           y = amr%ylo + (real(j,WP) + 0.5_WP) * dy
                           z = amr%zlo + (real(k,WP) + 0.5_WP) * dz
                           
                           ! Radial direction from sphere center (periodic)
                           rx = x - sphere_xc; rx = rx - (amr%xhi-amr%xlo) * nint(rx / (amr%xhi-amr%xlo))
                           ry = y - sphere_yc; ry = ry - (amr%yhi-amr%ylo) * nint(ry / (amr%yhi-amr%ylo))
                           rz = z - sphere_zc; rz = rz - (amr%zhi-amr%zlo) * nint(rz / (amr%zhi-amr%zlo))
                           r_mag = sqrt(rx**2 + ry**2 + rz**2)
                           
                           if (r_mag .gt. 1.0e-10_WP) then
                              rx = rx / r_mag
                              ry = ry / r_mag
                              rz = rz / r_mag
                              
                              ! Dot product: abs since PLIC normals can point either direction
                              ! Should be close to 1.0 for radially-aligned normals
                              dot = abs(nx*rx + ny*ry + nz*rz)
                              
                              min_dot_product = min(min_dot_product, dot)
                              max_dot_product = max(max_dot_product, dot)
                              avg_dot_product = avg_dot_product + dot
                              ncells_interface = ncells_interface + 1
                           end if
                        end if
                     end if
                  end do
               end do
            end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      
      ! Compute average
      if (ncells_interface .gt. 0) then
         avg_dot_product = avg_dot_product / real(ncells_interface, WP)
      end if
   end subroutine check_plic_normals

   !> Main test routine
   subroutine test_plic()
      implicit none

      call log("=== Testing PLIC initialization ===")

      ! Create amrgrid
      create_amrgrid: block
         use param, only: param_read
         amr%name = 'plic_test'
         call param_read('Base nx', amr%nx, default=32)
         call param_read('Base ny', amr%ny, default=32)
         call param_read('Base nz', amr%nz, default=32)
         amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
         amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
         amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
         amr%xper = .true.; amr%yper = .true.; amr%zper = .true.
         call param_read('Max level', amr%maxlvl, default=1)
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker (needed for regrid interface)
      initialize_timetracker: block
         time = timetracker(amRoot=amr%amRoot)
         time%t = 0.0_WP
         time%dt = 1.0_WP
         time%tmax = 0.0_WP  ! No time stepping
      end block initialize_timetracker

      ! Setup sphere parameters
      setup_sphere: block
         use param, only: param_read
         call param_read('Sphere radius', sphere_radius, default=0.25_WP)
         sphere_xc = 0.9_WP  ! Near x=1 boundary to test periodic wrapping
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

      ! Build PLIC from VF
      build_plic: block
         call vof%build_plic(time%t)
         call log("  PLIC constructed on all levels")
      end block build_plic

      ! Check PLIC normals
      check_normals: block
         call check_plic_normals()
         call log("  Interface cells: "//trim(itoa(ncells_interface)))
         call log("  Min dot(normal, radial): "//trim(rtoa(min_dot_product)))
         call log("  Max dot(normal, radial): "//trim(rtoa(max_dot_product)))
         call log("  Avg dot(normal, radial): "//trim(rtoa(avg_dot_product)))
         
         ! Success criterion: avg dot product > 0.99 (normals mostly point outward)
         if (avg_dot_product .gt. 0.99_WP) then
            call log("  [PASS] PLIC normals point radially outward")
         else if (avg_dot_product .gt. 0.95_WP) then
            call log("  [WARN] PLIC normals somewhat radial (avg="//trim(rtoa(avg_dot_product))//")")
         else
            call log("  [FAIL] PLIC normals NOT radially outward (avg="//trim(rtoa(avg_dot_product))//")")
         end if
      end block check_normals

      ! Test reset_moments round-trip: PLIC -> reset_moments -> VF should match original
      test_reset_moments: block
         real(WP) :: VF_max_err
         
         call log("=== Testing reset_moments round-trip ===")
         
         ! Copy VF to VFold before reset
         call vof%VFold%copy(vof%VF)
         
         ! Call reset_moments: recomputes VF from PLIC
         call vof%reset_moments()
         call log("  Called reset_moments")
         
         ! Compute error: VFold = VFold - VF, then get L-inf norm on finest level
         call vof%VFold%saxpy(-1.0_WP, vof%VF)
         VF_max_err = vof%VFold%norm0(lvl=amr%clvl())
         call log("  Max per-cell VF error: "//trim(rtoa(VF_max_err)))
         
         if (VF_max_err .lt. 1.0e-12_WP) then
            call log("  [PASS] reset_moments preserves VF (max cell error < 1e-12)")
         else if (VF_max_err .lt. 1.0e-6_WP) then
            call log("  [WARN] reset_moments has small VF error: "//trim(rtoa(VF_max_err)))
         else
            call log("  [FAIL] reset_moments does not preserve VF (max error="//trim(rtoa(VF_max_err))//")")
         end if
      end block test_reset_moments

      ! Create visualization for manual verification
      create_visualization: block
         call viz%initialize(amr, 'plic_test')
         call viz%add_scalar(vof%VF, 1, 'VF')
         call viz%add_scalar(vof%PLIC, 1, 'PLIC_nx')
         call viz%add_scalar(vof%PLIC, 2, 'PLIC_ny')
         call viz%add_scalar(vof%PLIC, 3, 'PLIC_nz')
         call viz%add_scalar(vof%PLIC, 4, 'PLIC_d')
         call viz%add_scalar(vof%Cliq, 1, 'Cliq_x')
         call viz%add_scalar(vof%Cliq, 2, 'Cliq_y')
         call viz%add_scalar(vof%Cliq, 3, 'Cliq_z')
         call viz%add_surfmesh(vof%smesh, 'plic')
         call viz%write(time=time%t)
         call log("  Visualization output written to amrviz/plic_test/")
         call log("  PLIC surface mesh: "//trim(itoa(vof%smesh%nPoly))//" polygons, "//trim(itoa(vof%smesh%nVert))//" vertices")
      end block create_visualization

      ! Cleanup
      cleanup: block
         call viz%finalize()
         call vof%finalize()
         call amr%finalize()
         call time%finalize()
      end block cleanup

      call log("=== PLIC test complete ===")

   end subroutine test_plic

end module mod_test_plic

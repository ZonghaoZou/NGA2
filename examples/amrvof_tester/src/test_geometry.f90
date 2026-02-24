!> Test amrvof geometry functions - cut_tet_vol, get_plane_dist
module mod_test_geometry
   use precision, only: WP
   use messager, only: log, die
   use string, only: rtoa
   implicit none
   private
   public :: test_geometry

contains

   !> Main test routine for geometry functions
   subroutine test_geometry()
      use amrvof_geometry, only: cut_tet_vol, get_plane_dist
      real(WP) :: tol, max_err
      
      call log("=== Testing amrvof_geometry ===")
      
      tol = 1.0e-12_WP
      max_err = 0.0_WP
      
      ! Test 1: cut_tet_vol with plane completely above tet (full tet below plane)
      test_full_tet: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_liq, vol_gas, vol_exact, err
         real(WP), dimension(3) :: bary_liq, bary_gas
         
         ! Unit tet at origin: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Plane above: n=(0,0,1), d=2 (z=2, all vertices have z<2 so all liquid)
         plane = [0.0_WP, 0.0_WP, 1.0_WP, 2.0_WP]
         
         call cut_tet_vol(tet, plane, vol_liq, vol_gas, bary_liq, bary_gas)
         vol_exact = 1.0_WP / 6.0_WP  ! Volume of unit tet
         err = abs(vol_liq - vol_exact)
         max_err = max(max_err, err)
         
         call log("  Test 1 (full tet): vol_liq = "//rtoa(vol_liq)//", err = "//rtoa(err))
         if (err .gt. tol) call die("Test 1 failed: cut_tet_vol full tet")
      end block test_full_tet
      
      ! Test 2: cut_tet_vol with plane completely below tet (zero liquid volume)
      test_empty_tet: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_liq, vol_gas, err
         real(WP), dimension(3) :: bary_liq, bary_gas
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Plane below: n=(0,0,1), d=-1 (z=-1, all vertices have z>-1 so all gas)
         plane = [0.0_WP, 0.0_WP, 1.0_WP, -1.0_WP]
         
         call cut_tet_vol(tet, plane, vol_liq, vol_gas, bary_liq, bary_gas)
         err = abs(vol_liq)
         max_err = max(max_err, err)
         
         call log("  Test 2 (empty tet): vol_liq = "//rtoa(vol_liq)//", err = "//rtoa(err))
         if (err .gt. tol) call die("Test 2 failed: cut_tet_vol empty tet")
      end block test_empty_tet
      
      ! Test 3: cut_tet_vol conservation (vol_liq + vol_gas = total)
      test_conservation: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_liq, vol_gas, vol_total, vol_exact, err
         real(WP), dimension(3) :: bary_liq, bary_gas
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Plane at z=0.3
         plane = [0.0_WP, 0.0_WP, 1.0_WP, 0.3_WP]
         
         call cut_tet_vol(tet, plane, vol_liq, vol_gas, bary_liq, bary_gas)
         vol_total = vol_liq + vol_gas
         vol_exact = 1.0_WP / 6.0_WP
         err = abs(vol_total - vol_exact)
         max_err = max(max_err, err)
         
         call log("  Test 3 (conservation): vol_liq+vol_gas = "//rtoa(vol_total)//", err = "//rtoa(err))
         if (err .gt. tol) call die("Test 3 failed: volume conservation")
      end block test_conservation
      
      ! Test 4: get_plane_dist with VF = 0.5 on unit cube
      test_plane_dist: block
         real(WP), dimension(3) :: normal, lo, hi
         real(WP) :: VF, dist, dist_expected, err
         
         ! Unit cube from (0,0,0) to (1,1,1), plane perpendicular to x-axis
         lo = [0.0_WP, 0.0_WP, 0.0_WP]
         hi = [1.0_WP, 1.0_WP, 1.0_WP]
         normal = [1.0_WP, 0.0_WP, 0.0_WP]
         VF = 0.5_WP
         
         dist = get_plane_dist(normal, lo, hi, VF)
         dist_expected = 0.5_WP
         err = abs(dist - dist_expected)
         max_err = max(max_err, err)
         
         call log("  Test 4 (plane_dist VF=0.5): dist = "//rtoa(dist)//" (expected "//rtoa(dist_expected)//")")
         if (err .gt. 1.0e-6_WP) call die("Test 4 failed: get_plane_dist")
      end block test_plane_dist
      
      ! Test 5: get_plane_dist with VF = 0 and VF = 1
      test_extreme_vf: block
         real(WP), dimension(3) :: normal, lo, hi
         real(WP) :: dist0, dist1
         
         lo = [0.0_WP, 0.0_WP, 0.0_WP]
         hi = [1.0_WP, 1.0_WP, 1.0_WP]
         normal = [1.0_WP, 0.0_WP, 0.0_WP]
         
         dist0 = get_plane_dist(normal, lo, hi, 0.0_WP)
         dist1 = get_plane_dist(normal, lo, hi, 1.0_WP)
         
         call log("  Test 5 (extreme VF): dist0 = "//rtoa(dist0)//", dist1 = "//rtoa(dist1))
         ! dist0 should be <= 0 (plane at or before origin)
         ! dist1 should be >= 1 (plane at or after far end)
         if (dist0 .gt. 1.0e-10_WP) call die("Test 5 failed: VF=0 should give dist<=0")
         if (dist1 .lt. 1.0_WP - 1.0e-10_WP) call die("Test 5 failed: VF=1 should give dist>=1")
      end block test_extreme_vf
      
      call log("=== All geometry tests passed! Max error: "//rtoa(max_err)//" ===")
      
      ! === Polyhedron clipping tests ===
      call test_poly_clipping()
      
   end subroutine test_geometry
   
   !> Test polyhedron clipping routines (VOFTools-style in-place clipping)
   subroutine test_poly_clipping()
      use amrvof_geometry, only: tet_vol, convex_poly, tet_to_poly, &
         clip_poly_by_plane, poly_vol_centroid, poly_vol
      real(WP) :: tol
      
      call log("=== Testing poly clipping ===")
      tol = 1.0e-12_WP
      
      ! Test P1: tet_to_poly volume matches tet_vol
      test_tet_to_poly: block
         real(WP), dimension(3,4) :: tet
         type(convex_poly) :: p
         real(WP) :: vol_tet, vol_poly, err
         real(WP), dimension(3) :: centroid
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         vol_tet = tet_vol(tet)
         call tet_to_poly(tet, p)
         call poly_vol_centroid(p, vol_poly, centroid)
         
         err = abs(vol_poly - vol_tet)
         call log("  P1 (tet_to_poly): vol_poly="//rtoa(vol_poly)//" vol_tet="//rtoa(vol_tet)//" err="//rtoa(err))
         if (err.gt.tol) call die("P1 failed: tet_to_poly volume mismatch")
         
         ! Also check centroid = (0.25, 0.25, 0.25)
         err = maxval(abs(centroid - [0.25_WP, 0.25_WP, 0.25_WP]))
         call log("  P1 centroid: ("//rtoa(centroid(1))//","//rtoa(centroid(2))//","//rtoa(centroid(3))//") err="//rtoa(err))
         if (err.gt.tol) call die("P1 failed: tet centroid mismatch")
      end block test_tet_to_poly
      
      ! Test P2: clip tet by x=0.5, volumes sum to original
      ! clip_poly_by_plane(poly, normal, c, icontp, icontn)
      ! Plane: n·x + c = 0. Keep n·x + c > 0 side.
      ! For "keep x < 0.5": normal=(-1,0,0), c=0.5 => -x+0.5>0 => x<0.5
      ! For "keep x > 0.5": normal=(1,0,0), c=-0.5 => x-0.5>0 => x>0.5
      test_clip_x: block
         real(WP), dimension(3,4) :: tet
         type(convex_poly) :: p_lo, p_hi
         real(WP) :: vol_tet, vol_lo, vol_hi, err
         real(WP), dimension(3) :: c_lo, c_hi
         integer :: icontp, icontn, ierr
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         vol_tet = tet_vol(tet)
         
         ! Get x < 0.5 piece
         call tet_to_poly(tet, p_lo)
         call clip_poly_by_plane(p_lo, [1.0_WP, 0.0_WP, 0.0_WP], 0.5_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p_lo, vol_lo, c_lo)
         
         ! Get x > 0.5 piece
         call tet_to_poly(tet, p_hi)
         call clip_poly_by_plane(p_hi, [-1.0_WP, 0.0_WP, 0.0_WP], -0.5_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p_hi, vol_hi, c_hi)
         
         err = abs(vol_lo + vol_hi - vol_tet)
         call log("  P2 (clip x=0.5): V_lo="//rtoa(vol_lo)//" V_hi="//rtoa(vol_hi)//" err="//rtoa(err))
         if (err.gt.tol) call die("P2 failed: volume conservation after x-clip")
      end block test_clip_x
      
      ! Test P3: clip tet by y=0.3, volumes sum to original
      test_clip_y: block
         real(WP), dimension(3,4) :: tet
         type(convex_poly) :: p_lo, p_hi
         real(WP) :: vol_tet, vol_lo, vol_hi, err
         real(WP), dimension(3) :: c_lo, c_hi
         integer :: icontp, icontn, ierr
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         vol_tet = tet_vol(tet)
         
         call tet_to_poly(tet, p_lo)
         call clip_poly_by_plane(p_lo, [0.0_WP, 1.0_WP, 0.0_WP], 0.3_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p_lo, vol_lo, c_lo)
         
         call tet_to_poly(tet, p_hi)
         call clip_poly_by_plane(p_hi, [0.0_WP, -1.0_WP, 0.0_WP], -0.3_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p_hi, vol_hi, c_hi)
         
         err = abs(vol_lo + vol_hi - vol_tet)
         call log("  P3 (clip y=0.3): V_lo="//rtoa(vol_lo)//" V_hi="//rtoa(vol_hi)//" err="//rtoa(err))
         if (err.gt.tol) call die("P3 failed: volume conservation after y-clip")
      end block test_clip_y
      
      ! Test P4: trivial clip (plane past tet — all on one side)
      test_clip_trivial: block
         real(WP), dimension(3,4) :: tet
         type(convex_poly) :: p
         real(WP) :: vol_tet, vol_clip, err
         real(WP), dimension(3) :: centroid
         integer :: icontp, icontn, ierr
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         vol_tet = tet_vol(tet)
         
         ! Keep x < 2.0 side (all vertices have x < 2, so full tet is kept)
         call tet_to_poly(tet, p)
         call clip_poly_by_plane(p, [1.0_WP, 0.0_WP, 0.0_WP], 2.0_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p, vol_clip, centroid)
         
         err = abs(vol_clip - vol_tet)
         call log("  P4 (trivial): vol_clip="//rtoa(vol_clip)//" icontp="//rtoa(real(icontp,WP))//" err="//rtoa(err))
         if (err.gt.tol) call die("P4 failed: trivial clip")
      end block test_clip_trivial
      
      ! Test P5: sequential clips — clip by x<0.3, then clip remainder by y<0.4
      test_sequential_clips: block
         real(WP), dimension(3,4) :: tet
         type(convex_poly) :: p1, p2, p3
         real(WP) :: vol_tet, v1, v2, v3, err
         real(WP), dimension(3) :: c1, c2, c3
         integer :: icontp, icontn, ierr
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         vol_tet = tet_vol(tet)
         
         ! Piece 1: x < 0.3
         call tet_to_poly(tet, p1)
         call clip_poly_by_plane(p1, [1.0_WP, 0.0_WP, 0.0_WP], 0.3_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p1, v1, c1)
         
         ! Piece 2: x >= 0.3 AND y < 0.4
         call tet_to_poly(tet, p2)
         call clip_poly_by_plane(p2, [-1.0_WP, 0.0_WP, 0.0_WP], -0.3_WP, icontp, icontn, ierr)
         call clip_poly_by_plane(p2, [0.0_WP, 1.0_WP, 0.0_WP], 0.4_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p2, v2, c2)
         
         ! Piece 3: x >= 0.3 AND y >= 0.4
         call tet_to_poly(tet, p3)
         call clip_poly_by_plane(p3, [-1.0_WP, 0.0_WP, 0.0_WP], -0.3_WP, icontp, icontn, ierr)
         call clip_poly_by_plane(p3, [0.0_WP, -1.0_WP, 0.0_WP], -0.4_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p3, v3, c3)
         
         err = abs(v1 + v2 + v3 - vol_tet)
         call log("  P5 (sequential): pieces="//rtoa(v1)//","//rtoa(v2)//","//rtoa(v3)//" err="//rtoa(err))
         if (err.gt.tol) call die("P5 failed: sequential clip volume conservation")
      end block test_sequential_clips
      
      
      call log("=== All poly clipping tests passed! ===")
      
      ! === Poly vs Tet decomposition comparison tests ===
      call test_poly_vs_tet()
      
   end subroutine test_poly_clipping
   
   
   !> Compare poly clipping path against cut_tet_vol (tet decomposition) path
   !! to verify identical geometric results
   subroutine test_poly_vs_tet()
      use amrvof_geometry, only: tet_vol, convex_poly, tet_to_poly, &
         clip_poly_by_plane, poly_vol_centroid, poly_vol, cut_tet_vol
      real(WP) :: tol
      
      call log("=== Testing poly vs tet decomposition ===")
      tol = 1.0e-12_WP
      
      ! T1: PLIC-only — diagonal plane through unit tet
      ! cut_tet_vol: plane=[nx,ny,nz,dist], liquid = n·x < dist
      ! clip_poly_by_plane: keeps n·x < dist — direct match!
      test_plic_only: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_liq_tet, vol_gas_tet
         real(WP), dimension(3) :: bary_liq_tet, bary_gas_tet
         type(convex_poly) :: p_liq, p_gas
         real(WP) :: vol_liq_poly, vol_gas_poly
         real(WP), dimension(3) :: bary_liq_poly, bary_gas_poly
         integer :: icontp, icontn, ierr
         real(WP) :: err_vliq, err_vgas, err_bliq, err_bgas
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Diagonal plane: n=(1,1,1)/sqrt(3), dist=0.3
         plane(1:3) = [1.0_WP, 1.0_WP, 1.0_WP] / sqrt(3.0_WP)
         plane(4) = 0.3_WP
         
         ! Tet decomposition path
         call cut_tet_vol(tet, plane, vol_liq_tet, vol_gas_tet, bary_liq_tet, bary_gas_tet)
         
         ! Poly clipping path — liquid side (n·x < dist) — direct call!
         call tet_to_poly(tet, p_liq)
         call clip_poly_by_plane(p_liq, plane(1:3), plane(4), icontp, icontn, ierr)
         call poly_vol_centroid(p_liq, vol_liq_poly, bary_liq_poly)
         vol_liq_poly = abs(vol_liq_poly)
         
         ! Poly clipping path — gas side (n·x > dist)
         call tet_to_poly(tet, p_gas)
         call clip_poly_by_plane(p_gas, -plane(1:3), -plane(4), icontp, icontn, ierr)
         call poly_vol_centroid(p_gas, vol_gas_poly, bary_gas_poly)
         vol_gas_poly = abs(vol_gas_poly)
         
         err_vliq = abs(vol_liq_poly - vol_liq_tet)
         err_vgas = abs(vol_gas_poly - vol_gas_tet)
         err_bliq = maxval(abs(bary_liq_poly - bary_liq_tet))
         err_bgas = maxval(abs(bary_gas_poly - bary_gas_tet))
         
         call log("  T1 (PLIC only): Vliq err="//rtoa(err_vliq)//" Vgas err="//rtoa(err_vgas))
         call log("  T1 centroids:   Bliq err="//rtoa(err_bliq)//" Bgas err="//rtoa(err_bgas))
         if (err_vliq.gt.tol) call die("T1 failed: liquid volume mismatch")
         if (err_vgas.gt.tol) call die("T1 failed: gas volume mismatch")
         if (err_bliq.gt.tol) call die("T1 failed: liquid centroid mismatch")
         if (err_bgas.gt.tol) call die("T1 failed: gas centroid mismatch")
      end block test_plic_only
      
      ! T2: Grid-plane split — tet spanning x=[0,1], cut at x=0.4
      ! cut_tet_vol with plane n=(1,0,0), dist=0.4:
      !   liquid = x < 0.4 (lo), gas = x > 0.4 (hi)
      test_grid_plane: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_lo_tet, vol_hi_tet
         real(WP), dimension(3) :: bary_lo_tet, bary_hi_tet
         type(convex_poly) :: p_lo, p_hi
         real(WP) :: vol_lo_poly, vol_hi_poly
         real(WP), dimension(3) :: bary_lo_poly, bary_hi_poly
         integer :: icontp, icontn, ierr
         real(WP) :: err_vlo, err_vhi, err_blo, err_bhi
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Grid plane at x=0.4: n=(1,0,0), dist=0.4
         plane = [1.0_WP, 0.0_WP, 0.0_WP, 0.4_WP]
         
         ! Tet decomposition: liq = x<0.4 (lo), gas = x>0.4 (hi)
         call cut_tet_vol(tet, plane, vol_lo_tet, vol_hi_tet, bary_lo_tet, bary_hi_tet)
         
         ! Poly: x < 0.4 piece
         call tet_to_poly(tet, p_lo)
         call clip_poly_by_plane(p_lo, [1.0_WP, 0.0_WP, 0.0_WP], 0.4_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p_lo, vol_lo_poly, bary_lo_poly)
         vol_lo_poly = abs(vol_lo_poly)
         
         ! Poly: x > 0.4 piece
         call tet_to_poly(tet, p_hi)
         call clip_poly_by_plane(p_hi, [-1.0_WP, 0.0_WP, 0.0_WP], -0.4_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p_hi, vol_hi_poly, bary_hi_poly)
         vol_hi_poly = abs(vol_hi_poly)
         
         err_vlo = abs(vol_lo_poly - vol_lo_tet)
         err_vhi = abs(vol_hi_poly - vol_hi_tet)
         err_blo = maxval(abs(bary_lo_poly - bary_lo_tet))
         err_bhi = maxval(abs(bary_hi_poly - bary_hi_tet))
         
         call log("  T2 (grid plane): Vlo err="//rtoa(err_vlo)//" Vhi err="//rtoa(err_vhi))
         call log("  T2 centroids:    Blo err="//rtoa(err_blo)//" Bhi err="//rtoa(err_bhi))
         if (err_vlo.gt.tol) call die("T2 failed: lo volume mismatch")
         if (err_vhi.gt.tol) call die("T2 failed: hi volume mismatch")
         if (err_blo.gt.tol) call die("T2 failed: lo centroid mismatch")
         if (err_bhi.gt.tol) call die("T2 failed: hi centroid mismatch")
      end block test_grid_plane
      
      ! T3: Grid-plane + PLIC — tet cut by x=0.4, then each piece cut by diagonal PLIC
      test_grid_plus_plic: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plic_plane
         type(convex_poly) :: p_lo, p_hi, p_lo_liq, p_lo_gas, p_hi_liq, p_hi_gas
         real(WP) :: vol_lo_liq, vol_lo_gas, vol_hi_liq, vol_hi_gas
         real(WP), dimension(3) :: bary_lo_liq, bary_lo_gas, bary_hi_liq, bary_hi_gas
         real(WP) :: vol_tot_liq, vol_tot_gas, vol_tet, err
         integer :: icontp, icontn, ierr
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         vol_tet = abs(tet_vol(tet))
         
         plic_plane(1:3) = [1.0_WP, 1.0_WP, 1.0_WP] / sqrt(3.0_WP)
         plic_plane(4) = 0.3_WP
         
         ! Step 1: Cut by x=0.4 (keep lo=x<0.4, hi=x>0.4)
         call tet_to_poly(tet, p_lo)
         call clip_poly_by_plane(p_lo, [1.0_WP, 0.0_WP, 0.0_WP], 0.4_WP, icontp, icontn, ierr)
         
         call tet_to_poly(tet, p_hi)
         call clip_poly_by_plane(p_hi, [-1.0_WP, 0.0_WP, 0.0_WP], -0.4_WP, icontp, icontn, ierr)
         
         ! Step 2: Cut each piece by PLIC (liquid = n·x < dist — direct call!)
         ! lo-piece liquid
         p_lo_liq = p_lo
         call clip_poly_by_plane(p_lo_liq, plic_plane(1:3), plic_plane(4), icontp, icontn, ierr)
         call poly_vol_centroid(p_lo_liq, vol_lo_liq, bary_lo_liq)
         vol_lo_liq = abs(vol_lo_liq)
         
         ! lo-piece gas (n·x > dist)
         p_lo_gas = p_lo
         call clip_poly_by_plane(p_lo_gas, -plic_plane(1:3), -plic_plane(4), icontp, icontn, ierr)
         call poly_vol_centroid(p_lo_gas, vol_lo_gas, bary_lo_gas)
         vol_lo_gas = abs(vol_lo_gas)
         
         ! hi-piece liquid (direct call)
         p_hi_liq = p_hi
         call clip_poly_by_plane(p_hi_liq, plic_plane(1:3), plic_plane(4), icontp, icontn, ierr)
         call poly_vol_centroid(p_hi_liq, vol_hi_liq, bary_hi_liq)
         vol_hi_liq = abs(vol_hi_liq)
         
         ! hi-piece gas (n·x > dist)
         p_hi_gas = p_hi
         call clip_poly_by_plane(p_hi_gas, -plic_plane(1:3), -plic_plane(4), icontp, icontn, ierr)
         call poly_vol_centroid(p_hi_gas, vol_hi_gas, bary_hi_gas)
         vol_hi_gas = abs(vol_hi_gas)
         
         ! Total liq and gas should sum to vol_tet
         vol_tot_liq = vol_lo_liq + vol_hi_liq
         vol_tot_gas = vol_lo_gas + vol_hi_gas
         err = abs(vol_tot_liq + vol_tot_gas - vol_tet)
         
         call log("  T3 (grid+PLIC): Vliq="//rtoa(vol_tot_liq)//" Vgas="//rtoa(vol_tot_gas)//" conservation err="//rtoa(err))
         if (err.gt.tol) call die("T3 failed: total volume conservation")
         
         ! Also compare total liq/gas against direct PLIC cut of full tet
         compare_with_direct: block
            real(WP) :: vol_liq_direct, vol_gas_direct
            real(WP), dimension(3) :: bary_liq_direct, bary_gas_direct
            real(WP) :: err_liq, err_gas
            
            call cut_tet_vol(tet, plic_plane, vol_liq_direct, vol_gas_direct, bary_liq_direct, bary_gas_direct)
            err_liq = abs(vol_tot_liq - vol_liq_direct)
            err_gas = abs(vol_tot_gas - vol_gas_direct)
            
            call log("  T3 vs direct:   Vliq err="//rtoa(err_liq)//" Vgas err="//rtoa(err_gas))
            if (err_liq.gt.tol) call die("T3 failed: liq volume vs direct PLIC mismatch")
            if (err_gas.gt.tol) call die("T3 failed: gas volume vs direct PLIC mismatch")
         end block compare_with_direct
      end block test_grid_plus_plic
      
      
      ! T4: Sign preservation — verify poly_vol sign matches tet_vol sign
      !     through tet_to_poly AND through clipping
      test_sign_preservation: block
         real(WP), dimension(3,4) :: tet_pos, tet_neg
         type(convex_poly) :: p_pos, p_neg, p_pos_clip, p_neg_clip
         real(WP) :: vol_pos_tet, vol_neg_tet
         real(WP) :: vol_pos_poly, vol_neg_poly
         real(WP) :: vol_pos_clip, vol_neg_clip
         real(WP), dimension(3) :: bary
         integer :: icontp, icontn, ierr
         
         ! Positive-orientation tet
         tet_pos(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet_pos(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet_pos(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet_pos(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         vol_pos_tet = tet_vol(tet_pos)
         
         ! Negative-orientation tet (swap v3 and v4)
         tet_neg(:,1) = tet_pos(:,1)
         tet_neg(:,2) = tet_pos(:,2)
         tet_neg(:,3) = tet_pos(:,4)  ! swapped
         tet_neg(:,4) = tet_pos(:,3)  ! swapped
         vol_neg_tet = tet_vol(tet_neg)
         
         ! Verify tet_vol signs
         if (vol_pos_tet.le.0.0_WP) call die("T4 failed: pos tet should have positive tet_vol")
         if (vol_neg_tet.ge.0.0_WP) call die("T4 failed: neg tet should have negative tet_vol")
         
         ! Check poly_vol matches tet_vol sign (unclipped)
         call tet_to_poly(tet_pos, p_pos)
         vol_pos_poly = poly_vol(p_pos)
         call tet_to_poly(tet_neg, p_neg)
         vol_neg_poly = poly_vol(p_neg)
         
         if (vol_pos_poly.le.0.0_WP) call die("T4 failed: pos poly_vol should be positive")
         if (vol_neg_poly.ge.0.0_WP) call die("T4 failed: neg poly_vol should be negative")
         if (abs(vol_pos_poly - vol_pos_tet).gt.tol) call die("T4 failed: pos poly_vol != tet_vol")
         if (abs(vol_neg_poly - vol_neg_tet).gt.tol) call die("T4 failed: neg poly_vol != tet_vol")
         
         ! Now clip both by a PLIC plane and check sign is preserved
         p_pos_clip = p_pos
         call clip_poly_by_plane(p_pos_clip, [1.0_WP, 1.0_WP, 1.0_WP]/sqrt(3.0_WP), 0.3_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p_pos_clip, vol_pos_clip, bary)
         
         p_neg_clip = p_neg
         call clip_poly_by_plane(p_neg_clip, [1.0_WP, 1.0_WP, 1.0_WP]/sqrt(3.0_WP), 0.3_WP, icontp, icontn, ierr)
         call poly_vol_centroid(p_neg_clip, vol_neg_clip, bary)
         
         call log("  T4 (sign):      pos_vol="//rtoa(vol_pos_clip)//" neg_vol="//rtoa(vol_neg_clip))
         if (vol_pos_clip.le.0.0_WP) call die("T4 failed: clipped pos poly should have positive vol")
         if (vol_neg_clip.ge.0.0_WP) call die("T4 failed: clipped neg poly should have negative vol")
         if (abs(abs(vol_pos_clip) - abs(vol_neg_clip)).gt.tol) call die("T4 failed: magnitudes should match")
      end block test_sign_preservation
      
      call log("=== All poly vs tet tests passed! ===")
      
   end subroutine test_poly_vs_tet

end module mod_test_geometry

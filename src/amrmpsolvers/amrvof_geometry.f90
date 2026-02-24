!> AMR VOF Geometry module
!> Provides cutting tables and volume computation for native VOF geometry
module amrvof_geometry
   use precision, only: WP
   implicit none
   private

   ! Expose tables and routines
   public :: tet_map, cut_side, cut_v1, cut_v2, cut_vtet
   public :: cut_ntets, cut_nvert, cut_nntet
   public :: get_plane_dist
   public :: tet_vol, tet_sign, cut_tet_vol
   public :: flux_polyhedron_vol, cut_hex_vol
   public :: correct_flux_poly
   public :: cut_hex_polygon, hex_poly_nvert, get_hex_poly_nvert
   public :: convex_poly, tet_to_poly, clip_poly_by_plane
   public :: poly_vol_centroid, poly_vol

   ! Cutting tables from mpcomp_class_noirl
   ! tet_map: maps a hex cell (8 vertices + center) to 8 tetrahedra
   integer, dimension(4,8), parameter :: tet_map = reshape([ &
      7, 4, 3, 6, 6, 3, 2, 4, 6, 2, 1, 4, 7, 8, 4, 6, &
      6, 5, 8, 4, 6, 5, 4, 1, 5, 6, 8, 9, 6, 7, 8, 9], shape(tet_map))

   ! cut_side: which side of plane each tet vertex is on (1=below, 2=above)
   integer, dimension(6,16), parameter :: cut_side = reshape([ &
      1,-1,-1,-1,-1,-1, 2, 1, 1, 1,-1,-1, 2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, &
      2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, &
      2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, &
      2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, 2, 2, 2, 1,-1,-1, 2,-1,-1,-1,-1,-1], shape(cut_side))

   ! cut_v1, cut_v2: edge endpoints for cut vertices
   integer, dimension(4,16), parameter :: cut_v1 = reshape([ &
      -1,-1,-1,-1, 1, 1, 1,-1, 2, 2, 2,-1, 1, 2, 1, 2, 3, 3, 3,-1, 1, 3, 1, 3, &
      2, 3, 2, 3, 4, 4, 4,-1, 4, 4, 4,-1, 1, 4, 1, 4, 2, 4, 2, 4, 3, 3, 3,-1, &
      3, 4, 3, 4, 2, 2, 2,-1, 1, 1, 1,-1,-1,-1,-1,-1], shape(cut_v1))

   integer, dimension(4,16), parameter :: cut_v2 = reshape([ &
      -1,-1,-1,-1, 2, 3, 4,-1, 3, 4, 1,-1, 4, 4, 3, 3, 4, 1, 2,-1, 4, 4, 2, 2, &
      4, 4, 1, 1, 1, 2, 3,-1, 1, 2, 3,-1, 3, 3, 2, 2, 3, 3, 1, 1, 4, 1, 2,-1, &
      2, 2, 1, 1, 3, 4, 1,-1, 2, 3, 4,-1,-1,-1,-1,-1], shape(cut_v2))

   ! cut_vtet: sub-tetrahedra for each cutting configuration
   integer, dimension(4,6,16), parameter :: cut_vtet = reshape([ &
      1, 2, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 7, 6, 1, 6, 2, 3, 4, 4, 2, 5, 6, 5, 6, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1, &
      7, 5, 6, 2, 1, 3, 4, 6, 1, 5, 3, 6, 5, 7, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 8, 6, 2, 5, 7, 8, 1, 5, 1, 8, 2, 5, 6, 8, 4, 5, 8, 7, 3, 5, 8, 3, 4, &
      6, 5, 7, 3, 2, 1, 4, 6, 6, 5, 4, 2, 6, 7, 5, 2,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 6, 8, 3, 5, 8, 7, 1, 5, 8, 1, 3, 5, 8, 6, 4, 5, 7, 8, 2, 5, 8, 4, 2, &
      8, 6, 5, 3, 5, 7, 8, 2, 8, 5, 2, 3, 8, 5, 6, 4, 5, 8, 7, 1, 5, 8, 1, 4, &
      1, 2, 3, 7, 1, 2, 7, 6, 5, 7, 6, 1, 5, 6, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 6, 7, 4, 1, 2, 3, 6, 5, 1, 3, 6, 5, 7, 6, 3,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 8, 6, 4, 5, 7, 8, 1, 5, 8, 4, 1, 5, 6, 8, 3, 5, 8, 7, 2, 5, 8, 2, 3, &
      8, 5, 6, 4, 5, 8, 7, 2, 8, 2, 5, 4, 8, 6, 5, 3, 5, 7, 8, 1, 5, 8, 3, 1, &
      1, 4, 2, 7, 4, 1, 6, 7, 6, 7, 5, 4, 6, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1,-1, &
      8, 6, 5, 4, 5, 7, 8, 3, 8, 4, 5, 3, 8, 5, 6, 2, 5, 8, 7, 1, 5, 8, 1, 2, &
      3, 4, 1, 7, 7, 6, 3, 4, 7, 6, 5, 3, 7, 5, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1, &
      7, 4, 2, 3, 2, 3, 6, 7, 5, 6, 7, 2, 5, 7, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1, &
      1, 2, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1], shape(cut_vtet))

   ! Number of resulting tetrahedra for each config
   integer, dimension(16), parameter :: cut_ntets = [1,4,4,6,4,6,6,4,4,6,6,4,6,4,4,1]
   ! Number of cut vertices for each config
   integer, dimension(16), parameter :: cut_nvert = [0,3,3,4,3,4,4,3,3,4,4,3,4,3,3,0]
   ! Number of resulting tetrahedra on one side
   integer, dimension(16), parameter :: cut_nntet = [1,2,2,4,2,4,4,4,2,4,4,4,4,4,4,2]

   !> ============================================================================
   !> Hex polygon cutting tables (from IRL lookup_tables.h)
   !> For cutting a plane by a hexahedron to get the intersection polygon
   !> Vertex ordering (x,y,z):
   !>   0:(+,-,-), 1:(+,+,-), 2:(+,+,+), 3:(+,-,+)
   !>   4:(-,-,-), 5:(-,+,-), 6:(-,+,+), 7:(-,-,+)
   !> ============================================================================
   
   ! Number of polygon vertices for each of 256 cases
   ! -1 indicates invalid/ambiguous case
   integer, dimension(0:255), parameter :: hex_poly_nvert = [ &
      0,  3,  3,  4,  3, -1,  4,  5,  3,  4, -1,  5,  4,  5,  5,  4, &
      3,  4, -1,  5, -1, -1, -1, -1, -1,  5, -1,  6, -1, -1, -1,  5, &
      3, -1,  4,  5, -1, -1,  5,  6, -1, -1, -1, -1, -1, -1, -1,  5, &
      4,  5,  5,  4, -1, -1, -1,  5, -1, -1, -1,  5, -1, -1, -1,  4, &
      3, -1, -1, -1,  4, -1,  5, -1, -1, -1, -1, -1,  5, -1,  6,  5, &
     -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      4, -1,  5, -1,  5, -1,  4,  5, -1, -1, -1, -1, -1, -1,  5,  4, &
      5, -1,  6,  5, -1, -1,  5,  4, -1, -1, -1, -1, -1, -1, -1,  3, &
      3, -1, -1, -1, -1, -1, -1, -1,  4,  5, -1, -1,  5,  6, -1,  5, &
      4,  5, -1, -1, -1, -1, -1, -1,  5,  4, -1,  5, -1,  5, -1,  4, &
     -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      5,  6, -1,  5, -1, -1, -1, -1, -1,  5, -1,  4, -1, -1, -1,  3, &
      4, -1, -1, -1,  5, -1, -1, -1,  5, -1, -1, -1,  4,  5,  5,  4, &
      5, -1, -1, -1, -1, -1, -1, -1,  6,  5, -1, -1,  5,  4, -1,  3, &
      5, -1, -1, -1,  6, -1,  5, -1, -1, -1, -1, -1,  5, -1,  4,  3, &
      4,  5,  5,  4,  5, -1,  4,  3,  5,  4, -1,  3,  4,  3,  3,  0]

   ! Edge pairs for polygon vertices: hex_poly_edges(1:2, vertex, case)
   ! Each pair (v1,v2) defines an edge; interpolate to find intersection point
   ! Stored as hex_poly_v1 and hex_poly_v2 for the two endpoints
   integer, dimension(6,0:255), parameter :: hex_poly_v1 = reshape([ &
      -1,-1,-1,-1,-1,-1,  0, 3, 0,-1,-1,-1,  0, 1, 1,-1,-1,-1,  1, 3, 0, 1,-1,-1, &
       1, 2, 2,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 1, 2, 2,-1,-1,  2, 3, 0, 1, 2,-1, &
       2, 3, 3,-1,-1,-1,  0, 2, 3, 0,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 0, 1,-1, &
       1, 2, 3, 3,-1,-1,  0, 1, 2, 3, 0,-1,  0, 1, 2, 3, 3,-1,  0, 1, 2, 3,-1,-1, &
       4, 0, 7,-1,-1,-1,  0, 3, 7, 4,-1,-1, -1,-1,-1,-1,-1,-1,  1, 3, 7, 4, 1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 2, 3, 7, 4,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 7, 4, 1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  4, 1, 2, 3, 7,-1, &
       4, 5, 1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 5, 1,-1,-1,  1, 3, 0, 4, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 5, 2, 2,-1,  2, 3, 0, 4, 5, 2, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  4, 5, 2, 3, 0,-1, &
       5, 1, 0, 7,-1,-1,  0, 3, 7, 5, 1,-1,  0, 0, 7, 5, 1,-1,  1, 3, 7, 5,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 5, 2,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 7, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 2, 3, 7,-1,-1, &
       5, 6, 2,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 6, 2,-1,-1, -1,-1,-1,-1,-1,-1,  0, 1, 5, 6, 2,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 6, 3, 3,-1, -1,-1,-1,-1,-1,-1,  0, 1, 5, 6, 3, 3,  5, 6, 3, 0, 1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       4, 6, 2, 1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 2, 1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 6, 2,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 2,-1,-1,  2, 3, 0, 4, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 3, 3,-1,  4, 6, 3, 0,-1,-1, &
       6, 2, 1, 0, 7,-1, -1,-1,-1,-1,-1,-1,  0, 0, 7, 6, 2, 1,  1, 3, 7, 6, 2,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 0, 7, 6, 2,-1,  2, 3, 7, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 3, 7,-1,-1,-1, &
       6, 7, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 3,-1,-1,  0, 2, 6, 7, 0,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 2, 6, 7, 3,-1,  0, 1, 2, 6, 7, 0, -1,-1,-1,-1,-1,-1,  6, 7, 0, 1, 2,-1, &
       4, 0, 3, 6,-1,-1,  0, 3, 3, 6, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 0, 3,-1,  0, 2, 6, 4,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 6, 4, 1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 1, 2, 6, 4,-1, -1,-1,-1,-1,-1,-1,  4, 1, 2, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       5, 1, 0, 3, 6,-1,  0, 3, 3, 6, 5, 1, -1,-1,-1,-1,-1,-1,  1, 3, 3, 6, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 2, 6, 5, 1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 6, 5,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 2, 6,-1,-1,-1, &
       5, 7, 3, 2,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 7, 3, 2,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 2, 5, 7, 3,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 7, 3,-1,-1,  0, 1, 5, 7, 0,-1,  0, 1, 5, 7, 3,-1,  5, 7, 0, 1,-1,-1, &
       4, 0, 3, 2, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 2, 5, 4, 0, 3,  0, 2, 2, 5, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 4, 0, 3,-1,  0, 1, 5, 4,-1,-1, -1,-1,-1,-1,-1,-1,  4, 1, 5,-1,-1,-1, &
       4, 7, 3, 2, 1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 7, 3, 2, -1,-1,-1,-1,-1,-1,  0, 4, 7, 3, 2,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1,  0, 4, 7, 3,-1,-1,  4, 7, 0,-1,-1,-1, &
       0, 3, 2, 1,-1,-1,  0, 3, 3, 2, 1,-1,  0, 0, 3, 2, 1,-1,  1, 3, 3, 2,-1,-1, &
       1, 1, 0, 3, 2,-1, -1,-1,-1,-1,-1,-1,  0, 0, 3, 2,-1,-1,  2, 3, 3,-1,-1,-1, &
       2, 2, 1, 0, 3,-1,  0, 2, 2, 1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 2,-1,-1,-1, &
       1, 1, 0, 3,-1,-1,  0, 1, 1,-1,-1,-1,  0, 0, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1], &
      shape=[6,256])

   integer, dimension(6,0:255), parameter :: hex_poly_v2 = reshape([ &
      -1,-1,-1,-1,-1,-1,  1, 0, 4,-1,-1,-1,  1, 5, 2,-1,-1,-1,  2, 0, 4, 5,-1,-1, &
       2, 6, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 3,-1,-1,  3, 0, 4, 5, 6,-1, &
       3, 7, 0,-1,-1,-1,  1, 3, 7, 4,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 5,-1, &
       2, 6, 7, 0,-1,-1,  1, 2, 6, 7, 4,-1,  1, 5, 6, 7, 0,-1,  4, 5, 6, 7,-1,-1, &
       5, 4, 4,-1,-1,-1,  1, 0, 4, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 0, 4, 5, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  1, 3, 7, 4, 5,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 5, 5, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6, 7, 4,-1, &
       5, 6, 5,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 2,-1,-1,  2, 0, 4, 5, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 6, 3,-1,  3, 0, 4, 5, 6, 6, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 6, 6, 7, 4,-1, &
       6, 5, 4, 4,-1,-1,  1, 0, 4, 6, 5,-1,  1, 4, 4, 6, 2,-1,  2, 0, 4, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  3, 0, 4, 6, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 6, 7, 4,-1,-1, &
       6, 7, 6,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 3,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 7, 3,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 7, 0,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 7, 7, 0,  6, 7, 7, 4, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       5, 7, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 6, 2,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 7, 3,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 3,-1,-1,  3, 0, 4, 5, 7,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 7, 0,-1,  5, 7, 7, 4,-1,-1, &
       7, 6, 5, 4, 4,-1, -1,-1,-1,-1,-1,-1,  1, 4, 4, 7, 6, 2,  2, 0, 4, 7, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 4, 4, 7, 3,-1,  3, 0, 4, 7,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  7, 7, 4,-1,-1,-1, &
       7, 4, 7,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 7, 4, 0,-1,-1,  1, 3, 7, 4, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 4, 0,-1,  1, 2, 6, 7, 4, 4, -1,-1,-1,-1,-1,-1,  7, 4, 4, 5, 6,-1, &
       5, 4, 7, 7,-1,-1,  1, 0, 7, 7, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 7, 5, 4, 0,-1,  1, 3, 7, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 5, 5,-1, &
      -1,-1,-1,-1,-1,-1,  1, 2, 6, 7, 5,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6, 7,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       6, 5, 4, 7, 7,-1,  1, 0, 7, 7, 6, 5, -1,-1,-1,-1,-1,-1,  2, 0, 7, 7, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  1, 3, 7, 6, 5,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 6, 7,-1,-1,-1, &
       6, 4, 7, 6,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 6, 6, 4, 0,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 0,-1,-1,  1, 2, 6, 4, 4,-1,  1, 5, 6, 4, 0,-1,  6, 4, 4, 5,-1,-1, &
       5, 4, 7, 6, 6,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 6, 6, 5, 4, 0,  1, 3, 6, 6, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 5, 4, 0,-1,  1, 2, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6,-1,-1,-1, &
       5, 4, 7, 6, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 4, 7, 3, -1,-1,-1,-1,-1,-1,  1, 5, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 4, 0,-1, -1,-1,-1,-1,-1,-1,  1, 5, 4, 0,-1,-1,  5, 4, 4,-1,-1,-1, &
       4, 7, 6, 5,-1,-1,  1, 0, 7, 6, 5,-1,  1, 4, 7, 6, 2,-1,  2, 0, 7, 6,-1,-1, &
       2, 5, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1,  1, 4, 7, 3,-1,-1,  3, 0, 7,-1,-1,-1, &
       3, 6, 5, 4, 0,-1,  1, 3, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 6,-1,-1,-1, &
       2, 5, 4, 0,-1,-1,  1, 2, 5,-1,-1,-1,  1, 4, 0,-1,-1,-1, -1,-1,-1,-1,-1,-1], &
      shape=[6,256])

   ! =========================================================================
   ! Convex polyhedron type for polyhedron clipping approach
   ! =========================================================================
   integer, parameter, public :: POLY_MAXV  = 40  ! max vertex indices
   integer, parameter, public :: POLY_MAXF  = 24  ! max faces
   integer, parameter, public :: POLY_MAXFV = 8   ! max vertices per face
   
   type, public :: convex_poly
      integer  :: nv  = 0                         ! number of active vertices
      integer  :: nf  = 0                         ! number of faces
      integer  :: ntp = 0                         ! last vertex index (can exceed nv)
      real(WP) :: v(3,POLY_MAXV)                  ! vertex coordinates
      integer  :: fnv(POLY_MAXF)                  ! vertices per face
      integer  :: fv(POLY_MAXFV, POLY_MAXF)       ! face->vertex indices
   end type convex_poly

contains

   !> Compute plane distance d such that plane (n,d) cuts cell to give target VF
   !> Cell is axis-aligned box defined by corners lo and hi
   !> Normal n must be unit length
   !> Returns d in global coordinates: n.x = d
   !> Reference: Scardovelli & Zaleski, JCP 164, 228-247 (2000)
   pure function get_plane_dist(normal, lo, hi, VF) result(d)
      implicit none
      real(WP), dimension(3), intent(in) :: normal  !< Unit normal
      real(WP), dimension(3), intent(in) :: lo      !< Lower cell corner
      real(WP), dimension(3), intent(in) :: hi      !< Upper cell corner
      real(WP), intent(in) :: VF                    !< Target volume fraction
      real(WP) :: d
      real(WP) :: m1, m2, m3, mm1, mm2, mm3, tmp
      real(WP) :: norm, factor, alpha, VOFo
      real(WP) :: m12, V1, V2, V3
      real(WP) :: a0, a1, a2, p0, q0, theta
      real(WP), parameter :: eps = 1.0e-15_WP
      real(WP), dimension(3) :: cellsize, cellctr
      
      ! Compute cell size and centroid from corners
      cellsize = hi - lo
      cellctr  = 0.5_WP * (lo + hi)
      
      ! Form scaled normal components
      mm1 = normal(1) * cellsize(1)
      mm2 = normal(2) * cellsize(2)
      mm3 = normal(3) * cellsize(3)
      
      ! Normalize to sum of absolute values = 1
      norm = abs(mm1) + abs(mm2) + abs(mm3)
      if (norm.lt.eps) then
         d = sign(1.0e10_WP, VF - 0.5_WP)
         return
      end if
      mm1 = mm1 / norm
      mm2 = mm2 / norm
      mm3 = mm3 / norm

      
      ! Track offset for negative components
      factor = 0.0_WP
      if (mm1.lt.0.0_WP) factor = factor + mm1
      if (mm2.lt.0.0_WP) factor = factor + mm2
      if (mm3.lt.0.0_WP) factor = factor + mm3
      
      ! Handle VF > 0.5 by symmetry
      VOFo = VF
      if (VF.gt.0.5_WP) VOFo = 1.0_WP - VF
      
      ! === Core S&Z algorithm ===
      ! Take absolute value and sort so that m1 <= m2 <= m3
      m1 = abs(mm1); m2 = abs(mm2); m3 = abs(mm3)
      if (m2.lt.m1) then; tmp = m2; m2 = m1; m1 = tmp; end if
      if (m3.lt.m2) then; tmp = m3; m3 = m2; m2 = tmp; end if
      if (m2.lt.m1) then; tmp = m2; m2 = m1; m1 = tmp; end if
      
      ! Form volume thresholds
      m12 = m1 + m2
      V1 = m1**2 / max(6.0_WP * m2 * m3, eps)
      V2 = V1 + 0.5_WP * (m2 - m1) / max(m3, eps)
      if (m12.le.m3) then
         V3 = 0.5_WP * m12 / max(m3, eps)
      else
         V3 = (m3**2 * (3.0_WP*m12 - m3) + m1**2 * (m1 - 3.0_WP*m3) + &
              m2**2 * (m2 - 3.0_WP*m3)) / max(6.0_WP * m1 * m2 * m3, eps)
      end if
      
      ! Calculate alpha based on regime
      if (VOFo.ge.0.0_WP.and.VOFo.lt.V1) then
         alpha = (6.0_WP * m1 * m2 * m3 * VOFo)**(1.0_WP/3.0_WP)
      else if (VOFo.lt.V2) then
         alpha = 0.5_WP * (m1 + sqrt(m1**2 + 8.0_WP * m2 * m3 * (VOFo - V1)))
      else if (VOFo.lt.V3) then
         a0 = -(m1**3 + m2**3 - 6.0_WP * m1 * m2 * m3 * VOFo)
         a1 = 3.0_WP * (m1**2 + m2**2)
         a2 = -3.0_WP * m12
         p0 = -(a1 / 3.0_WP - a2**2 / 9.0_WP)
         q0 = (a1 * a2 - 3.0_WP * a0) / 6.0_WP - a2**3 / 27.0_WP
         theta = acos(max(-1.0_WP, min(1.0_WP, q0 / sqrt(p0**3)))) / 3.0_WP
         alpha = sqrt(p0) * (sqrt(3.0_WP) * sin(theta) - cos(theta)) - a2 / 3.0_WP
      else
         if (m12.le.m3) then
            alpha = m3 * VOFo + 0.5_WP * m12
         else
            a0 = -0.5_WP * (m1**3 + m2**3 + m3**3 - 6.0_WP * m1 * m2 * m3 * VOFo)
            a1 = 1.5_WP * (m1**2 + m2**2 + m3**2)
            a2 = -1.5_WP
            p0 = -(a1 / 3.0_WP - a2**2 / 9.0_WP)
            q0 = (a1 * a2 - 3.0_WP * a0) / 6.0_WP - a2**3 / 27.0_WP
            theta = acos(max(-1.0_WP, min(1.0_WP, q0 / sqrt(p0**3)))) / 3.0_WP
            alpha = sqrt(p0) * (sqrt(3.0_WP) * sin(theta) - cos(theta)) - a2 / 3.0_WP
         end if
      end if
      ! === End core S&Z algorithm ===
      
      ! Undo VF > 0.5 symmetry
      if (VF.gt.0.5_WP) alpha = 1.0_WP - alpha
      
      ! Adjust for negative normal components
      alpha = alpha + factor
      
      ! Transform from cell-corner origin to cell-centroid origin
      alpha = alpha - 0.5_WP * (mm1 + mm2 + mm3)
      
      ! Scale back to physical coordinates and add centroid offset
      d = alpha * norm + dot_product(normal, cellctr)
      
   end function get_plane_dist
   
   !> Compute volume of tetrahedron given 4 vertices
   !> v(:,1:4) are the vertex coordinates
   pure function tet_vol(v) result(vol)
      implicit none
      real(WP), dimension(3,4), intent(in) :: v
      real(WP) :: vol
      real(WP), dimension(3) :: a,b,c
      a=v(:,1)-v(:,4)
      b=v(:,2)-v(:,4)
      c=v(:,3)-v(:,4)
      vol=(-a(1)*(b(2)*c(3)-c(2)*b(3))&
      &    +a(2)*(b(1)*c(3)-c(1)*b(3))&
      &    -a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
   end function tet_vol

   !> Function that calculates the sign of a tet
   function tet_sign(vert) result(s)
      implicit none
      real(WP) :: s
      real(WP), dimension(3,4), intent(in) :: vert
      real(WP), dimension(3) :: a,b,c
      a=vert(:,1)-vert(:,4); b=vert(:,2)-vert(:,4); c=vert(:,3)-vert(:,4)
      s=sign(1.0_WP,-(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2))))
   end function tet_sign
   
   !> Cut a tetrahedron by a plane and return liquid/gas volumes and barycenters
   !> Input:  v(:,1:4) = 4 tet vertices, plane(1:4) = [nx,ny,nz,d] where n.x=d defines plane
   !> Output: vol_liq, vol_gas = phase volumes
   !>         bary_liq(3), bary_gas(3) = volume-weighted barycenters
   pure subroutine cut_tet_vol(v, plane, vol_liq, vol_gas, bary_liq, bary_gas)
      implicit none
      real(WP), dimension(3,4), intent(in)  :: v
      real(WP), dimension(4),   intent(in)  :: plane
      real(WP),                 intent(out) :: vol_liq, vol_gas
      real(WP), dimension(3),   intent(out) :: bary_liq, bary_gas
      real(WP), dimension(4)   :: d
      real(WP), dimension(3,8) :: vert
      real(WP), dimension(3)   :: a, b, c, bary
      real(WP) :: denom, mu, my_vol
      integer  :: icase, n1, v1, v2
      
      ! Compute signed distance of each vertex to plane
      d(1) = plane(1)*v(1,1) + plane(2)*v(2,1) + plane(3)*v(3,1) - plane(4)
      d(2) = plane(1)*v(1,2) + plane(2)*v(2,2) + plane(3)*v(3,2) - plane(4)
      d(3) = plane(1)*v(1,3) + plane(2)*v(2,3) + plane(3)*v(3,3) - plane(4)
      d(4) = plane(1)*v(1,4) + plane(2)*v(2,4) + plane(3)*v(3,4) - plane(4)
      
      ! Determine cutting case (1-16 based on which vertices are above plane)
      icase = 1 + int(0.5_WP+sign(0.5_WP,d(1))) + 2*int(0.5_WP+sign(0.5_WP,d(2))) &
                + 4*int(0.5_WP+sign(0.5_WP,d(3))) + 8*int(0.5_WP+sign(0.5_WP,d(4)))
      
      ! Copy original vertices
      vert(:,1:4) = v(:,1:4)
      
      ! Create interpolated vertices on cut plane
      do n1 = 1, cut_nvert(icase)
         v1 = cut_v1(n1,icase); v2 = cut_v2(n1,icase)
         denom = d(v2) - d(v1)
         if (abs(denom).ge.tiny(1.0_WP)) then
            mu = max(0.0_WP, min(1.0_WP, -d(v1)/denom))
         else
            mu = 0.0_WP
         end if
         vert(:,4+n1) = (1.0_WP-mu)*vert(:,v1) + mu*vert(:,v2)
      end do
      
      ! Initialize outputs
      vol_liq = 0.0_WP; vol_gas = 0.0_WP
      bary_liq = 0.0_WP; bary_gas = 0.0_WP
      
      ! Gas tets: from 1 to cut_nntet-1
      do n1 = 1, cut_nntet(icase)-1
         a = vert(:,cut_vtet(1,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         b = vert(:,cut_vtet(2,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         c = vert(:,cut_vtet(3,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
         bary = 0.25_WP * (vert(:,cut_vtet(1,n1,icase)) + vert(:,cut_vtet(2,n1,icase)) &
                         + vert(:,cut_vtet(3,n1,icase)) + vert(:,cut_vtet(4,n1,icase)))
         vol_gas = vol_gas + my_vol
         bary_gas = bary_gas + my_vol * bary
      end do
      
      ! Liquid tets: from cut_ntets down to cut_nntet
      do n1 = cut_ntets(icase), cut_nntet(icase), -1
         a = vert(:,cut_vtet(1,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         b = vert(:,cut_vtet(2,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         c = vert(:,cut_vtet(3,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
         bary = 0.25_WP * (vert(:,cut_vtet(1,n1,icase)) + vert(:,cut_vtet(2,n1,icase)) &
                         + vert(:,cut_vtet(3,n1,icase)) + vert(:,cut_vtet(4,n1,icase)))
         vol_liq = vol_liq + my_vol
         bary_liq = bary_liq + my_vol * bary
      end do
      
      ! Normalize barycenters
      if (vol_liq.gt.tiny(1.0_WP)) bary_liq = bary_liq / vol_liq
      if (vol_gas.gt.tiny(1.0_WP)) bary_gas = bary_gas / vol_gas
      
   end subroutine cut_tet_vol
   
   !> Compute signed volume of flux polyhedron
   !> face(:,1:8) = 8 vertices (4 at time t, 4 back-projected)
   !> Decomposes into 6 tets using tet_map and sums signed volumes
   pure function flux_polyhedron_vol(face) result(vol)
      implicit none
      real(WP), dimension(3,8), intent(in) :: face
      real(WP) :: vol
      real(WP), dimension(3) :: a, b, c
      integer :: ntet
      vol = 0.0_WP
      do ntet = 1, 6
         a = face(:,tet_map(1,ntet)) - face(:,tet_map(4,ntet))
         b = face(:,tet_map(2,ntet)) - face(:,tet_map(4,ntet))
         c = face(:,tet_map(3,ntet)) - face(:,tet_map(4,ntet))
         vol = vol + (a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
      end do
   end function flux_polyhedron_vol

   !> Adjust flux polyhedron to enforce target volume
   !> Uses IRL's direction-independent approach: moves vertex 9 along back-face normal
   subroutine correct_flux_poly(poly,target_volume)
      implicit none
      real(WP), dimension(3,9), intent(inout) :: poly
      real(WP), intent(in) :: target_volume
      real(WP) :: starting_volume,needed_change,mag,adjustment
      real(WP), dimension(3) :: cross_sum,dir
      real(WP), dimension(3) :: e1,e2,c1,c2
      integer :: n
      ! Compute starting volume
      starting_volume=0.0_WP
      do n=1,8; starting_volume=starting_volume+tet_vol([poly(:,tet_map(1,n)),poly(:,tet_map(2,n)),poly(:,tet_map(3,n)),poly(:,tet_map(4,n))]); end do
      needed_change=target_volume-starting_volume
      ! Compute volume gradient for tets 7 and 8 (the only tets using vertex 9)
      ! For tet (A,B,C,D), gradient w.r.t. D is -(1/6)*((B-A)×(C-A))
      ! Tet 7 = (5,6,8,9): gradient = -((v6-v5)×(v8-v5)). Take edges from vertex 5
      e1=poly(:,6)-poly(:,5); e2=poly(:,8)-poly(:,5)
      c1=[e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)]
      ! Tet 8 = (6,7,8,9): gradient = -((v7-v6)×(v8-v6)). Take edges from vertex 6
      e1=poly(:,7)-poly(:,6); e2=poly(:,8)-poly(:,6)
      c2=[e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)]
      ! Total gradient (negative sign absorbed into adjustment formula)
      cross_sum=c1+c2
      ! Compute adjustment along normal direction
      mag=sqrt(cross_sum(1)**2+cross_sum(2)**2+cross_sum(3)**2)
      adjustment=6.0_WP*needed_change/max(mag,tiny(1.0_WP))
      dir=cross_sum/max(mag,tiny(1.0_WP))
      ! Move vertex 9
      poly(:,9)=poly(:,9)+adjustment*dir
   end subroutine correct_flux_poly

   !> Cut a hex cell by a plane and compute liquid/gas volumes and barycenters
   !> hex(:,1:8) = 8 vertices of hex cell (standard ordering)
   !> plane(1:3) = normal, plane(4) = distance (n·x = d)
   !> Returns vol_liq/vol_gas and bary_liq/bary_gas
   subroutine cut_hex_vol(hex, plane, vol_liq, vol_gas, bary_liq, bary_gas)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex
      real(WP), dimension(4), intent(in) :: plane
      real(WP), intent(out) :: vol_liq, vol_gas
      real(WP), dimension(3), intent(out) :: bary_liq, bary_gas
      ! Local
      real(WP), dimension(3,4) :: tet
      real(WP), dimension(3) :: center
      real(WP) :: vl, vg
      real(WP), dimension(3) :: bl, bg
      integer :: ntet, v
      ! Hex->tet decomposition: 6 tets using center point
      ! hex_tet_map: indices into hex vertices for 6 tets (each shares center)
      integer, dimension(3,6), parameter :: hex_tet_map = reshape([ &
         1, 2, 3,  1, 3, 4,  1, 2, 6,  1, 6, 5, &
         5, 6, 7,  5, 7, 8], shape(hex_tet_map))
      integer, dimension(3,6), parameter :: hex_tet_top = reshape([ &
         5, 6, 7,  5, 7, 8,  4, 3, 7,  4, 7, 8, &
         1, 2, 3,  1, 3, 4], shape(hex_tet_top))
      ! Actually use standard 5-tet decomposition of hex
      integer, dimension(4,5), parameter :: tet5 = reshape([ &
         1, 2, 4, 5,  2, 3, 4, 7,  2, 5, 6, 7, &
         4, 5, 7, 8,  2, 4, 5, 7], shape(tet5))
      
      vol_liq = 0.0_WP; vol_gas = 0.0_WP
      bary_liq = 0.0_WP; bary_gas = 0.0_WP
      
      ! Decompose hex into 5 tetrahedra
      do ntet = 1, 5
         ! Build tetrahedron
         do v = 1, 4
            tet(:,v) = hex(:,tet5(v,ntet))
         end do
         ! Cut tet by plane
         call cut_tet_vol(tet, plane, vl, vg, bl, bg)
         ! Accumulate volumes and volume-weighted barycenters
         bary_liq = bary_liq + vl * bl
         bary_gas = bary_gas + vg * bg
         vol_liq = vol_liq + vl
         vol_gas = vol_gas + vg
      end do
      
      ! Normalize barycenters
      if (vol_liq.gt.tiny(1.0_WP)) bary_liq = bary_liq / vol_liq
      if (vol_gas.gt.tiny(1.0_WP)) bary_gas = bary_gas / vol_gas
      
   end subroutine cut_hex_vol

   !> =========================================================================
   !> Extract PLIC polygon from hex cell using IRL's algorithm
   !> Given 8 hex vertices and a plane (normal, distance), compute the 
   !> intersection polygon vertices.
   !> 
   !> Hex vertex ordering (IRL convention):
   !>   0:(+,-,-), 1:(+,+,-), 2:(+,+,+), 3:(+,-,+)
   !>   4:(-,-,-), 5:(-,+,-), 6:(-,+,+), 7:(-,-,+)
   !>
   !> For axis-aligned cell [xlo,ylo,zlo] to [xhi,yhi,zhi]:
   !>   hex(:,1) = [xhi, ylo, zlo]  ! vertex 0
   !>   hex(:,2) = [xhi, yhi, zlo]  ! vertex 1
   !>   hex(:,3) = [xhi, yhi, zhi]  ! vertex 2
   !>   hex(:,4) = [xhi, ylo, zhi]  ! vertex 3
   !>   hex(:,5) = [xlo, ylo, zlo]  ! vertex 4
   !>   hex(:,6) = [xlo, yhi, zlo]  ! vertex 5
   !>   hex(:,7) = [xlo, yhi, zhi]  ! vertex 6
   !>   hex(:,8) = [xlo, ylo, zhi]  ! vertex 7
   !>
   !> Returns:
   !>   nvert: number of polygon vertices (0-6)
   !>   poly(:,1:nvert): polygon vertices in cyclic order
   !> =========================================================================
   pure subroutine cut_hex_polygon(hex, plane, nvert, poly)
      implicit none
      real(WP), dimension(3,8), intent(in)  :: hex      !< 8 hex vertices
      real(WP), dimension(4),   intent(in)  :: plane    !< (nx, ny, nz, d)
      integer,                  intent(out) :: nvert    !< Number of polygon vertices
      real(WP), dimension(3,6), intent(out) :: poly     !< Polygon vertices (up to 6)
      
      ! Local variables
      real(WP), dimension(8) :: dist       ! Signed distance to plane for each vertex
      integer :: icase                     ! Cutting case (0-255)
      integer :: n, v1, v2
      real(WP) :: mu, denom
      
      ! Initialize output
      nvert = 0
      poly = 0.0_WP
      
      ! Compute signed distance from plane for each vertex
      ! Plane equation: n.x - d = 0, so distance = n.x - d
      do n = 1, 8
         dist(n) = plane(1)*hex(1,n) + plane(2)*hex(2,n) + plane(3)*hex(3,n) - plane(4)
      end do
      
      ! Determine cutting case from sign pattern (0-indexed vertices)
      ! Case bit i is set if vertex i is above plane (dist > 0)
      icase = 0
      if (dist(1).gt.0.0_WP) icase = icase + 1    ! vertex 0
      if (dist(2).gt.0.0_WP) icase = icase + 2    ! vertex 1
      if (dist(3).gt.0.0_WP) icase = icase + 4    ! vertex 2
      if (dist(4).gt.0.0_WP) icase = icase + 8    ! vertex 3
      if (dist(5).gt.0.0_WP) icase = icase + 16   ! vertex 4
      if (dist(6).gt.0.0_WP) icase = icase + 32   ! vertex 5
      if (dist(7).gt.0.0_WP) icase = icase + 64   ! vertex 6
      if (dist(8).gt.0.0_WP) icase = icase + 128  ! vertex 7
      
      ! Get number of polygon vertices for this case
      nvert = hex_poly_nvert(icase)
      
      ! Handle no intersection or ambiguous cases
      if (nvert.le.0) then
         nvert = 0
         return
      end if
      
      ! Compute intersection points along cut edges
      do n = 1, nvert
         ! Get edge endpoints (0-indexed in tables, convert to 1-indexed)
         v1 = hex_poly_v1(n, icase) + 1
         v2 = hex_poly_v2(n, icase) + 1
         
         ! Linear interpolation: find point where dist = 0
         ! mu = -dist(v1) / (dist(v2) - dist(v1))
         denom = dist(v2) - dist(v1)
         if (abs(denom).ge.tiny(1.0_WP)) then
            mu = max(0.0_WP, min(1.0_WP, -dist(v1)/denom))
         else
            mu = 0.5_WP
         end if
         
         ! Interpolated point
         poly(:,n) = (1.0_WP - mu)*hex(:,v1) + mu*hex(:,v2)
      end do
      
   end subroutine cut_hex_polygon
   
   !> =========================================================================
   !> Cheap vertex count for hex-plane intersection (no interpolation)
   !> Returns the number of polygon vertices from cutting case lookup only.
   !> This is much cheaper than cut_hex_polygon since it skips interpolation.
   !> =========================================================================
   pure function get_hex_poly_nvert(hex, plane) result(nvert)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex      !< 8 hex vertices
      real(WP), dimension(4),   intent(in) :: plane    !< (nx, ny, nz, d)
      integer :: nvert
      
      real(WP), dimension(8) :: dist
      integer :: icase, n
      
      ! Compute signed distance from plane for each vertex
      do n = 1, 8
         dist(n) = plane(1)*hex(1,n) + plane(2)*hex(2,n) + plane(3)*hex(3,n) - plane(4)
      end do
      
      ! Determine cutting case from sign pattern
      icase = 0
      if (dist(1).gt.0.0_WP) icase = icase + 1
      if (dist(2).gt.0.0_WP) icase = icase + 2
      if (dist(3).gt.0.0_WP) icase = icase + 4
      if (dist(4).gt.0.0_WP) icase = icase + 8
      if (dist(5).gt.0.0_WP) icase = icase + 16
      if (dist(6).gt.0.0_WP) icase = icase + 32
      if (dist(7).gt.0.0_WP) icase = icase + 64
      if (dist(8).gt.0.0_WP) icase = icase + 128
      
      ! Get vertex count from lookup table
      nvert = hex_poly_nvert(icase)
      if (nvert.lt.0) nvert = 0
      
   end function get_hex_poly_nvert


   ! =========================================================================
   ! Convert a tetrahedron to a convex_poly
   ! Face winding is CCW when viewed from outside for positive-volume tets
   ! =========================================================================
   pure subroutine tet_to_poly(tet, poly)
      implicit none
      real(WP), dimension(3,4), intent(in)  :: tet
      type(convex_poly),        intent(out) :: poly
      poly%nv  = 4
      poly%nf  = 4
      poly%ntp = 4
      poly%v(:,1:4) = tet(:,1:4)
      poly%fnv(1:4) = 3
      poly%fv(1:3,1) = [1,3,2]  ! opposite vertex 4
      poly%fv(1:3,2) = [1,2,4]  ! opposite vertex 3
      poly%fv(1:3,3) = [2,3,4]  ! opposite vertex 1
      poly%fv(1:3,4) = [1,4,3]  ! opposite vertex 2
   end subroutine tet_to_poly

   ! =========================================================================
   ! Clip convex polyhedron by plane: n·x = dist (keep n·x < dist side)
   ! Following VOFTools INTE3D + NEWPOL3D (Lopez & Hernandez)
   !
   ! Convention: plane is n·x = dist. Vertices with n·x < dist are KEPT.
   ! This directly matches PLIC: call with normal and dist to get liquid.
   !
   ! The polyhedron is modified IN-PLACE to contain only the kept side.
   ! Returns icontp (kept vertices) and icontn (removed vertices) counts.
   ! =========================================================================
   pure subroutine clip_poly_by_plane(poly, normal, dist, icontp, icontn, ierr)
      implicit none
      type(convex_poly),      intent(inout) :: poly
      real(WP), dimension(3), intent(in)    :: normal
      real(WP),               intent(in)    :: dist
      integer,                intent(out)   :: icontp, icontn
      integer,                intent(out)   :: ierr  ! 0=ok, 1=MAXV, 2=MAXFV, 3=MAXF
      ! Shared work arrays (needed across multiple blocks)
      integer  :: ia(POLY_MAXV)              ! vertex classification: 1=keep, 0=cut
      real(WP) :: phiv(POLY_MAXV)            ! signed distance to plane
      integer  :: iscut(POLY_MAXF)           ! face is cut flag
      integer  :: nedge(POLY_MAXF)           ! number of cut edges per face
      integer  :: ipia0(POLY_MAXV)           ! parent vertex 0 of new vertex (IA=0 side)
      integer  :: ipia1(POLY_MAXV)           ! parent vertex 1 of new vertex (IA=1 side)
      integer  :: ipv1(POLY_MAXFV,POLY_MAXF) ! work face vertex lists
      integer  :: nipv1(POLY_MAXF)           ! work face vertex counts
      integer  :: nts0, nts00                ! original face count
      integer  :: nipnew                     ! last vertex index after clipping
      
      ierr = 0
      
      ! ---------------------------------------------------------------
      ! Block 1: Classify vertices by signed distance to plane
      ! ---------------------------------------------------------------
      classify_vertices: block
         integer :: is, iv, ip
         ia = -1
         icontp = 0; icontn = 0
         do is = 1, poly%nf
            do iv = 1, poly%fnv(is)
               ip = poly%fv(iv, is)
               if (ia(ip).eq.(-1)) then
                  phiv(ip) = dist - (normal(1)*poly%v(1,ip) + normal(2)*poly%v(2,ip) + normal(3)*poly%v(3,ip))
                  if (phiv(ip).gt.0.0_WP) then
                     ia(ip) = 1
                     icontp = icontp + 1
                  else
                     ia(ip) = 0
                     icontn = icontn + 1
                  end if
               end if
            end do
         end do
      end block classify_vertices
      
      ! Trivial cases: all vertices on one side
      if (icontp.eq.0 .or. icontn.eq.0) return
      
      nts0  = poly%nf
      nts00 = nts0
      
      ! ---------------------------------------------------------------
      ! Block 2: Identify which faces are cut by the plane
      ! ---------------------------------------------------------------
      identify_cut_faces: block
         integer :: is, iv, iv1, ip, ip1
         do is = 1, nts0
            nedge(is) = 0
            if (poly%fnv(is).gt.0) then
               iscut(is) = 0
               do iv = 1, poly%fnv(is)
                  ip = poly%fv(iv, is)
                  iv1 = iv + 1; if (iv.eq.poly%fnv(is)) iv1 = 1
                  ip1 = poly%fv(iv1, is)
                  if (ia(ip).ne.ia(ip1)) then
                     iscut(is) = 1
                     nedge(is) = nedge(is) + 1
                  end if
               end do
               ! Mark uncut faces entirely on the removed side
               if (iscut(is).eq.0 .and. ia(poly%fv(1,is)).eq.0) &
                  poly%fnv(is) = -poly%fnv(is)
            end if
         end do
      end block identify_cut_faces
      
      ! ---------------------------------------------------------------
      ! Block 3: Clip cut faces and create intersection vertices
      ! ---------------------------------------------------------------
      clip_faces: block
         integer :: ise(POLY_MAXF,POLY_MAXFV)    ! intersection vertex per cut edge
         integer :: ivise(POLY_MAXF,POLY_MAXV)   ! position of intersection vertex in face
         integer :: ipise(POLY_MAXV,2)            ! face adjacency of intersection vertex
         integer :: is, iv, iv1, ip, ip1, ip0i, ip1i, itype
         integer :: niv, nint, is1, ie, ipnew, ip0n, ip1n
         logical :: found
         
         nipnew = poly%ntp
         do is = 1, nts0
            if (iscut(is).eq.1) then
               niv = 0; nint = 0
               do iv = 1, poly%fnv(is)
                  ip = poly%fv(iv, is)
                  iv1 = iv + 1; if (iv1.gt.poly%fnv(is)) iv1 = 1
                  ip1 = poly%fv(iv1, is)
                  ! Keep vertices on the kept side
                  if (ia(ip).eq.1) then
                     niv = niv + 1
                     if (niv.gt.POLY_MAXFV) then; ierr = 2; return; end if
                     ipv1(niv, is) = poly%fv(iv, is)
                  end if
                  ! Edge crosses plane — find or create intersection vertex
                  if (ia(ip).ne.ia(ip1)) then
                     nint = nint + 1
                     niv = niv + 1
                     if (niv.gt.POLY_MAXFV) then; ierr = 2; return; end if
                     if (ia(ip).eq.1) then
                        ip1i = ip; ip0i = ip1; itype = 2
                     else
                        ip1i = ip1; ip0i = ip; itype = 1
                     end if
                     ! Search previous faces for existing intersection on same edge
                     found = .false.
                     edge_search: do is1 = 1, is-1
                        do ie = 1, nedge(is1)
                           ipnew = ise(is1, ie)
                           ip0n = ipia0(ipnew); ip1n = ipia1(ipnew)
                           if (ip0n.eq.ip0i .and. ip1n.eq.ip1i) then
                              ise(is, nint) = ipnew
                              ipv1(niv, is) = ipnew
                              ivise(is, ipnew) = niv
                              ipise(ipnew, itype) = is
                              found = .true.
                              exit edge_search
                           end if
                        end do
                     end do edge_search
                     ! Create new intersection vertex
                     if (.not.found) then
                        nipnew = nipnew + 1
                        if (nipnew.gt.POLY_MAXV) then; ierr = 1; return; end if
                        ipia0(nipnew) = ip0i
                        ipia1(nipnew) = ip1i
                        ipv1(niv, is) = nipnew
                        ise(is, nint) = nipnew
                        ivise(is, nipnew) = niv
                        ipise(nipnew, itype) = is
                     end if
                  end if
               end do
               nipv1(is) = niv
            end if
         end do
         
         ! ---------------------------------------------------------------
         ! Block 4: Construct cap faces by chasing face adjacency
         ! ---------------------------------------------------------------
         build_cap_faces: block
            integer :: ipmark(POLY_MAXV)
            integer :: nivnew, isnew, ivnew, ivnewt, ipini
            
            nivnew = nipnew - poly%ntp
            isnew = nts0
            do ip = poly%ntp+1, nipnew
               ipmark(ip) = 0
            end do
            ivnewt = 0
            ipnew = poly%ntp + 1
            
            cap_loop: do while (ivnewt.lt.nivnew)
               ivnew = 1; ivnewt = ivnewt + 1
               isnew = isnew + 1
               if (isnew.gt.POLY_MAXF) then; ierr = 3; return; end if
               ipini = ipnew
               poly%fv(ivnew, isnew) = ipnew
               ipmark(ipnew) = 1
               chase: do while (ivnewt.lt.nivnew)
                  is = ipise(ipnew, 1)
                  iv = ivise(is, ipnew)
                  iv1 = iv - 1; if (iv1.eq.0) iv1 = nipv1(is)
                  ipnew = ipv1(iv1, is)
                  if (ipnew.eq.ipini) exit chase
                  ivnew = ivnew + 1; ivnewt = ivnewt + 1
                  poly%fv(ivnew, isnew) = ipnew
                  ipmark(ipnew) = 1
               end do chase
               poly%fnv(isnew) = ivnew
               do ipnew = poly%ntp+2, nipnew
                  if (ipmark(ipnew).eq.0) cycle cap_loop
               end do
               exit cap_loop
            end do cap_loop
            poly%fnv(isnew) = ivnew
            
            nts0 = isnew  ! update for finalize block
         end block build_cap_faces
         
      end block clip_faces
      
      ! ---------------------------------------------------------------
      ! Block 5: Finalize — update face lists, compute new vertex
      !           positions, set cap normals, update counts
      ! ---------------------------------------------------------------
      finalize: block
         integer :: is, iv, ip, niv, ip0i, ip1i
         real(WP) :: mu
         
         ! Copy clipped face vertex lists back to poly
         niv = nts0 - nts00  ! start with number of new (intersection) vertices
         do is = 1, nts00
            if (poly%fnv(is).gt.0) then
               if (iscut(is).eq.1) then
                  poly%fnv(is) = nipv1(is)
                  do iv = 1, nipv1(is)
                     poly%fv(iv, is) = ipv1(iv, is)
                     if (ia(ipv1(iv, is)).eq.1) then
                        niv = niv + 1
                        ia(ipv1(iv, is)) = -1
                     end if
                  end do
               else
                  if (iscut(is).eq.0 .and. poly%fnv(is).lt.0) poly%fnv(is) = 0
                  do iv = 1, poly%fnv(is)
                     if (ia(poly%fv(iv, is)).eq.1) then
                        niv = niv + 1
                        ia(poly%fv(iv, is)) = -1
                     end if
                  end do
               end if
            end if
         end do
         
         ! Reset IA flags
         do ip = 1, poly%ntp
            if (ia(ip).eq.(-1)) ia(ip) = 1
         end do
         
         ! Suppress degenerate cap faces
         do is = nts00+1, nts0
            if (poly%fnv(is).lt.3) poly%fnv(is) = 0
         end do
         
         ! Compute new vertex positions via linear interpolation
         do is = nts00+1, nts0
            do iv = 1, poly%fnv(is)
               ip = poly%fv(iv, is)
               ip0i = ipia0(ip); ip1i = ipia1(ip)
               mu = -phiv(ip0i) / (phiv(ip1i) - phiv(ip0i))
               poly%v(:,ip) = poly%v(:,ip0i) + mu * (poly%v(:,ip1i) - poly%v(:,ip0i))
            end do
         end do
         
         ! Update counts
         poly%nv  = niv
         poly%ntp = nipnew
         poly%nf  = nts0
      end block finalize
      
   end subroutine clip_poly_by_plane
   
   ! =========================================================================
   ! Compute volume and centroid of a convex polyhedron
   ! Uses divergence theorem via fan decomposition: each face is triangulated
   ! from its first vertex, each triangle + a reference point forms a signed
   ! tet whose volume and centroid are accumulated.
   ! =========================================================================
   pure subroutine poly_vol_centroid(poly, vol, centroid)
      implicit none
      type(convex_poly),      intent(in)  :: poly
      real(WP),               intent(out) :: vol
      real(WP), dimension(3), intent(out) :: centroid
      ! Local
      real(WP), dimension(3) :: ref, a, b, c, tc
      real(WP) :: tvol
      integer  :: f, e, v1i, v2i, v3i
      
      ref = poly%v(:,poly%fv(1,1))  ! any vertex as reference
      vol = 0.0_WP
      centroid = 0.0_WP
      
      do f = 1, poly%nf
         if (poly%fnv(f).le.0) cycle
         v1i = poly%fv(1, f)
         do e = 2, poly%fnv(f) - 1
            v2i = poly%fv(e, f)
            v3i = poly%fv(e+1, f)
            a = poly%v(:,v1i) - ref
            b = poly%v(:,v2i) - ref
            c = poly%v(:,v3i) - ref
            tvol = (a(1)*(b(2)*c(3)-c(2)*b(3)) &
            &      -a(2)*(b(1)*c(3)-c(1)*b(3)) &
            &      +a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
            tc = 0.25_WP * (ref + poly%v(:,v1i) + poly%v(:,v2i) + poly%v(:,v3i))
            vol = vol + tvol
            centroid = centroid + tvol * tc
         end do
      end do
      
      ! Normalize centroid
      if (abs(vol).gt.tiny(1.0_WP)) then
         centroid = centroid / vol
      end if
      
   end subroutine poly_vol_centroid

   ! =========================================================================
   ! Compute volume of a convex polyhedron (volume only, no centroid)
   ! =========================================================================
   pure function poly_vol(poly) result(vol)
      implicit none
      type(convex_poly), intent(in) :: poly
      real(WP) :: vol
      real(WP), dimension(3) :: ref, a, b, c
      integer :: f, e, v1i, v2i, v3i
      ref = poly%v(:,poly%fv(1,1))
      vol = 0.0_WP
      do f = 1, poly%nf
         if (poly%fnv(f).le.0) cycle
         v1i = poly%fv(1, f)
         do e = 2, poly%fnv(f) - 1
            v2i = poly%fv(e, f)
            v3i = poly%fv(e+1, f)
            a = poly%v(:,v1i) - ref
            b = poly%v(:,v2i) - ref
            c = poly%v(:,v3i) - ref
            vol = vol + (a(1)*(b(2)*c(3)-c(2)*b(3)) &
            &           -a(2)*(b(1)*c(3)-c(1)*b(3)) &
            &           +a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
         end do
      end do
   end function poly_vol

end module amrvof_geometry
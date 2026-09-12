! Grid definition for FESOM's YAC component.
!
! Builds the unstructured YAC grid from the FESOM2 mesh: one polygonal cell
! per node (element centres + edge midpoints, plus the node itself on the
! coast) and the node coordinates as the point set on which every field is
! defined. Also sets the YAC start datetime from the FESOM clock.
!
! Shared by every coupling interface FESOM registers on this component
! (atm_coupling_interface, ice_coupling_interface). The grid name is an
! argument so that the atmosphere-facing and ice-facing exchanges can move
! to different grids later without an API change.
module yac_grid_utils
#if defined (__cpl_yac)

  use yac
  use o_PARAM, only: WP, PI
  use g_clock, only: yearnew, month, day_in_month, timenew

  implicit none
  private

  public :: compute_midpoint, compute_center, cpl_yac_define_unstr_generic

contains

  subroutine compute_midpoint(geo_a, geo_b, mid)
    real(kind=WP), intent(in)  :: geo_a(2), geo_b(2)
    real(kind=WP), intent(out) :: mid(2)
    real(kind=WP) :: cos_lat_a, cos_lat_b, x_mid, y_mid, z_mid

    cos_lat_a = cos(geo_a(2))
    cos_lat_b = cos(geo_b(2))
    x_mid = 0.5_WP * (cos_lat_a * cos(geo_a(1)) + cos_lat_b * cos(geo_b(1)))
    y_mid = 0.5_WP * (cos_lat_a * sin(geo_a(1)) + cos_lat_b * sin(geo_b(1)))
    z_mid = 0.5_WP * (sin(geo_a(2)) + sin(geo_b(2)))

    mid(1) = atan2(y_mid, x_mid)
    mid(2) = PI/2 - acos(z_mid / sqrt(x_mid*x_mid + y_mid*y_mid + z_mid*z_mid))
  end subroutine compute_midpoint

  subroutine compute_center(geo_a, geo_b, geo_c, mid)
    real(kind=WP), intent(in)  :: geo_a(2), geo_b(2), geo_c(2)
    real(kind=WP), intent(out) :: mid(2)
    real(kind=WP) :: cos_lat_a, cos_lat_b, cos_lat_c, x_mid, y_mid, z_mid

    cos_lat_a = cos(geo_a(2))
    cos_lat_b = cos(geo_b(2))
    cos_lat_c = cos(geo_c(2))
    x_mid = (cos_lat_a * cos(geo_a(1)) + cos_lat_b * cos(geo_b(1)) + cos_lat_c * cos(geo_c(1))) / 3._WP
    y_mid = (cos_lat_a * sin(geo_a(1)) + cos_lat_b * sin(geo_b(1)) + cos_lat_c * sin(geo_c(1))) / 3._WP
    z_mid = (sin(geo_a(2)) + sin(geo_b(2)) + sin(geo_c(2))) / 3._WP

    mid(1) = atan2(y_mid, x_mid)
    mid(2) = PI/2 - acos(z_mid / sqrt(x_mid*x_mid + y_mid*y_mid + z_mid*z_mid))
  end subroutine compute_center

  ! Define the YAC datetime and the FESOM-mesh-derived unstructured grid.
  ! Callers add their own field definitions on the returned points_id.
  !
  ! Assumption: every interface on this component shares the FESOM2
  ! triangular mesh. If the sea ice ever moves to another mesh topology,
  ! this routine must be specialised per grid name.
  subroutine cpl_yac_define_unstr_generic(partit, mesh, grid_name, grid_id, points_id)
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_rotate_grid
    implicit none

    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    character(len=*), intent(in)  :: grid_name
    integer,          intent(out) :: grid_id, points_id

    real(kind=WP), allocatable :: x_vertices(:), y_vertices(:)
    real(kind=WP) :: mid(2)
    integer, allocatable :: nbr_vertices_per_cell(:), cell_to_vertex(:)
    integer :: i, j, k, nbr_vertices, nbr_boundary_nodes, nbr_connections, vtx_idx, c2v_idx
    integer :: curr_elem, curr_edge
    logical, allocatable :: node_is_boundary(:)
    character(LEN=24) :: startdatetime

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    WRITE(startdatetime, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,".",I3.3, "Z")') &
         yearnew, month, day_in_month, &
         INT(timenew)/3600, MODULO(INT(timenew), 3600)/60, MODULO(INT(timenew), 60), &
         INT(MODULO(timenew, 1.0_WP)*1000)

    CALL yac_fdef_datetime(startdatetime)

    ! find boundary nodes
    ALLOCATE(node_is_boundary(myDim_nod2D))
    node_is_boundary = .FALSE.
    DO i=1,myDim_edge2D
       IF (myList_edge2D(i) > edge2D_in) THEN
          IF (edges(1,i) <= myDim_nod2D) &
               node_is_boundary(edges(1,i)) = .TRUE.
          IF (edges(2,i) <= myDim_nod2D) &
               node_is_boundary(edges(2,i)) = .TRUE.
       END IF
    END DO
    nbr_boundary_nodes = COUNT(node_is_boundary)

    nbr_vertices = myDim_elem2D + myDim_edge2D + nbr_boundary_nodes
    nbr_connections = 2*SUM(nod_in_elem2D_num(1:myDim_nod2D)) + 2*nbr_boundary_nodes

    ALLOCATE(x_vertices(nbr_vertices))
    ALLOCATE(y_vertices(nbr_vertices))
    ALLOCATE(nbr_vertices_per_cell(myDim_nod2D))
    ALLOCATE(cell_to_vertex(nbr_connections))

    ! compute vertices
    ! 1 element centers
    DO i=1,myDim_elem2D
       CALL compute_center(geo_coord_nod2D(1:2,elem2D_nodes(1, i)), &
            geo_coord_nod2D(1:2,elem2D_nodes(2, i)), &
            geo_coord_nod2D(1:2,elem2D_nodes(3, i)), mid)
       x_vertices(i) = mid(1)
       y_vertices(i) = mid(2)
    END DO
    ! 2 edges midpoints
    DO i=1,myDim_edge2D
       CALL compute_midpoint(geo_coord_nod2D(1:2,edges(1,i)), geo_coord_nod2D(1:2,edges(2,i)), mid)
       x_vertices(myDim_elem2D + i) = mid(1)
       y_vertices(myDim_elem2D + i) = mid(2)
    END DO
    ! 3 boundary nodes -> are added on the fly when the cells are computed

    ! compute cells
    vtx_idx = myDim_elem2D + myDim_edge2D
    c2v_idx = 0
    DO i=1,myDim_nod2D
       nbr_vertices_per_cell(i) = 2*nod_in_elem2D_num(i)
       curr_elem = nod_in_elem2D(1,i)
       curr_edge = elem_edges(1, curr_elem)
       IF (ALL(edges(1:2,curr_edge) /= i)) curr_edge = elem_edges(2, curr_elem)
       IF (node_is_boundary(i)) THEN
          nbr_vertices_per_cell(i) = nbr_vertices_per_cell(i) + 2
          ! we're starting with the boundary node
          vtx_idx = vtx_idx+1
          x_vertices(vtx_idx) = geo_coord_nod2D(1,i)
          y_vertices(vtx_idx) = geo_coord_nod2D(2,i)
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = vtx_idx
          elem_loop: DO j=1,nod_in_elem2D_num(i)
             DO k=1,3
                IF (myList_edge2D(elem_edges(k,nod_in_elem2D(j,i))) > edge2D_in .AND. &
                     ANY(edges(1:2,elem_edges(k,nod_in_elem2D(j,i))) == i)) THEN
                   curr_elem = nod_in_elem2D(j,i)
                   curr_edge = elem_edges(k, curr_elem)
                   EXIT elem_loop
                END IF
             END DO
          END DO elem_loop
       END IF

       DO j=1,nod_in_elem2D_num(i)
          IF (curr_elem == 0) THEN
             WRITE (0,*) "Error: cells with two or more coast lines are currently not supported by the yac setup code"
             STOP
          ENDIF
          ! add the midpoint of curr_edge
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = myDim_elem2D + curr_edge
          ! add the center of curr_elem
          c2v_idx = c2v_idx + 1
          IF (curr_elem < 1) THEN
             WRITE (0,*) "elem 0 detected in cell ", i, " with ", nbr_vertices_per_cell(i), " vertices - j is ", j, " is_boundary ", node_is_boundary(i)
          END IF
          cell_to_vertex(c2v_idx) = curr_elem
          ! find next edge
          edge_loop: DO k=1,3
             IF (elem_edges(k, curr_elem) /= curr_edge .AND. &
                  ANY(edges(1:2,elem_edges(k, curr_elem)) == i)) THEN
                curr_edge = elem_edges(k, curr_elem)
                EXIT edge_loop
             END IF
          END DO edge_loop
          IF (edge_tri(1, curr_edge) /= curr_elem) THEN
             curr_elem = edge_tri(1, curr_edge)
          ELSE
             curr_elem = edge_tri(2, curr_edge)
          END IF
       END DO

       ! in case of a boundary volume, we add the others edges midpoint
       IF (node_is_boundary(i)) THEN
          ! add the midpoint of curr_edge
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = myDim_elem2D + curr_edge
       END IF
    END DO

    CALL yac_fdef_grid( &
         trim(grid_name), &
         nbr_vertices, &
         myDim_nod2D, &
         nbr_connections, &
         nbr_vertices_per_cell, &
         x_vertices, &
         y_vertices, &
         cell_to_vertex, &
         grid_id)

    CALL yac_fset_global_index( &
         partit%myList_nod2D - 1, &
         YAC_LOCATION_CELL, &
         grid_id)

    CALL yac_fdef_points(grid_id, &
         myDim_nod2D, &
         YAC_LOCATION_CELL, &
         geo_coord_nod2D(1,1:myDim_nod2D), &
         geo_coord_nod2D(2,1:myDim_nod2D), &
         points_id)

  end subroutine cpl_yac_define_unstr_generic

#endif
end module yac_grid_utils

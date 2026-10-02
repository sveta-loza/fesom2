module cpl_yac_driver
#if defined (__cpl_yac)

  USE yac
  USE o_PARAM
  USE g_clock
  USE cpl_config, only : cpl_comp_name, cpl_grid_name, cpl_config_file
  USE cpl_fields, only : nsend, nrecv, cpl_send, cpl_recv, &
                         cpl_send_collection_size, cpl_recv_collection_size, &
                         cpl_send_partner, cpl_send_steps, cpl_recv_steps, &
                         CPL_PARTNER_ATMOSPHERE

  implicit none
  save

  integer :: comp_id, grid_id, points_id
  integer, allocatable, dimension(:) :: send_field_id, recv_field_id
  real(kind=WP), dimension(:,:),   allocatable   :: a2o_fcorr_stat  !flux correction statistics for the output
  integer                    :: source_root, target_root   !this root/source in MPI_COMM_WORLD
  logical                    :: commRank       ! true for ranks doing YAC communication

  ! One-shot check, on rank 0, that coupling_period in coupling.yaml equals the
  ! exchange interval of the ocean <-> sea-ice fields (dt * cpl_stride).
  ! Nothing else enforces it: with a longer period YAC couples only every Nth
  ! call, and the component silently reuses the fields of the last coupling
  ! (zero-order hold), which drove surface run-aways on high-resolution meshes.
  integer, private :: n_put_ice = 0, n_put_ice_coupled = 0
  logical, private :: period_checked = .false.
  integer, private :: driver_mype = -1

  public nsend, nrecv
  public cpl_send, cpl_recv
  public cpl_send_collection_size, cpl_recv_collection_size
  public source_root, target_root, commRank
  public a2o_fcorr_stat


contains

  subroutine cpl_yac_init( localCommunicator )
    use ifs_mpmd_world, only: yac_world_comm, ifs_mpmd_have_yac_world
    implicit none

    integer, intent(out) :: localCommunicator

#ifdef VERBOSE
    print *, '================================================='
    print *, 'cpl_yac_init : coupler initialization for YAC'
    print *, '*************************************************'
#endif /* VERBOSE */

    ! IFS + FESIM: MPI_COMM_WORLD also holds the IFS IO servers and the FESOM
    ! multio servers, which never call YAC; both executables initialise YAC on
    ! the communicator ifs_mpmd_world joined them into instead.
    if (ifs_mpmd_have_yac_world()) then
       CALL yac_finit_comm(yac_world_comm)
    else
       CALL yac_finit()
    end if
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fread_config_yaml (trim(cpl_config_file))
    CALL yac_fdef_comp(trim(cpl_comp_name), comp_id)

    CALL yac_fget_comp_comm(comp_id, localCommunicator)

  end subroutine cpl_yac_init

  subroutine compute_midpoint(geo_a, geo_b, mid)
    implicit none

    REAL(kind=WP), intent(in) :: geo_a(2), geo_b(2)
    REAL(kind=WP), intent(out) :: mid(2)
    REAL(kind=WP) :: cos_lat_a, cos_lat_b, x_mid, y_mid, z_mid

    cos_lat_a = cos(geo_a(2))
    cos_lat_b = cos(geo_b(2))
    x_mid = 0.5*(cos_lat_a * cos(geo_a(1)) + cos_lat_b * cos(geo_b(1)))
    y_mid = 0.5*(cos_lat_a * sin(geo_a(1)) + cos_lat_b * sin(geo_b(1)))
    z_mid = 0.5*(sin(geo_a(2)) + sin(geo_b(2)))

    mid(1) = atan2(y_mid , x_mid)
    mid(2) = PI/2 - acos(z_mid/sqrt(x_mid*x_mid + y_mid*y_mid + z_mid*z_mid))
  end subroutine compute_midpoint

  subroutine compute_center(geo_a, geo_b, geo_c, mid)
    implicit none

    REAL(kind=WP), intent(in) :: geo_a(2), geo_b(2), geo_c(2)
    REAL(kind=WP), intent(out) :: mid(2)
    REAL(kind=WP) :: cos_lat_a, cos_lat_b, cos_lat_c, x_mid, y_mid, z_mid

    cos_lat_a = cos(geo_a(2))
    cos_lat_b = cos(geo_b(2))
    cos_lat_c= cos(geo_c(2))
    x_mid = (cos_lat_a * cos(geo_a(1)) + cos_lat_b * cos(geo_b(1)) + cos_lat_c * cos(geo_c(1))) / 3.
    y_mid = (cos_lat_a * sin(geo_a(1)) + cos_lat_b * sin(geo_b(1)) + cos_lat_c * sin(geo_c(1)) ) / 3.
    z_mid = (sin(geo_a(2)) + sin(geo_b(2)) + sin(geo_c(2))) /3.

    mid(1) = atan2(y_mid , x_mid)
    mid(2) = PI/2 - acos(z_mid/sqrt(x_mid*x_mid + y_mid*y_mid + z_mid*z_mid))
  end subroutine compute_center

  !> nsteps: the time steps of this run. Without an atmosphere among the YAC
  !> components (ICON-A) there is no YAC time master, and the components
  !> define the end of the run as well: start + nsteps*dt, from each one's own
  !> clock and namelist.config, so coupling.yaml carries no dates (a date
  !> pinned there breaks every restart). With ICON-A, ICON sets the end.
  subroutine cpl_yac_define_unstr(partit, mesh, nsteps)
    use mod_mesh
    USE MOD_PARTIT
    use g_rotate_grid
    use g_config, only: dt
    use cpl_config, only: is_coupled_to_icon_a

    implicit none

    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP), allocatable :: x_vertices(:), y_vertices(:)
    real(kind=WP) :: mid(2)
    integer, allocatable :: nbr_vertices_per_cell(:), cell_to_vertex(:)
    integer :: i, j, k, nbr_vertices, nbr_boundary_nodes
    integer :: nbr_connections, vtx_idx, c2v_idx
    integer :: curr_elem, curr_edge
    logical, allocatable :: node_is_boundary(:)
    integer, allocatable :: nbnd_edges(:)        ! boundary edges per node: 0 interior, 2 one coastal arc, 2k k arcs
    integer, allocatable :: star_vtx(:)
    real(kind=WP), allocatable :: star_ang(:)
    integer :: n_star, ed, l, itmp
    real(kind=WP) :: dx, dy, atmp
    integer, intent(in) :: nsteps
    character(len=8)           :: dt_str
    character(LEN=24) :: startdatetime, enddatetime

    driver_mype = partit%mype

    WRITE(startdatetime, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,".",I3.3, "Z")') &
         yearnew, month, day_in_month, &
         INT(timenew)/3600, MODULO(INT(timenew), 3600)/60, MODULO(INT(timenew), 60), &
         INT(MODULO(timenew, 1.0_WP)*1000)

    if (is_coupled_to_icon_a) then
       CALL yac_fdef_datetime(startdatetime)
    else
       call run_end_datetime(nsteps, enddatetime)
       CALL yac_fdef_datetime(startdatetime, enddatetime)
       if (partit%mype == 0) write(*,*) '     YAC run period: ', startdatetime, ' -> ', enddatetime
    end if

    ! find boundary nodes and count the boundary edges per node: 2 = one coastal
    ! arc (the usual coastal node), 2k = k separate arcs (a node touching land in
    ! k places, e.g. at a one-element-wide strait; tORCA025 has 7 such nodes,
    ! CORE2 and DART none)
    ALLOCATE(nbnd_edges(partit%myDim_nod2D), node_is_boundary(partit%myDim_nod2D))
    nbnd_edges = 0
    DO i=1,partit%myDim_edge2D
       IF (partit%myList_edge2D(i) > mesh%edge2D_in) THEN
          IF (mesh%edges(1,i) <= partit%myDim_nod2D) &
               nbnd_edges(mesh%edges(1,i)) = nbnd_edges(mesh%edges(1,i)) + 1
          IF (mesh%edges(2,i) <= partit%myDim_nod2D) &
               nbnd_edges(mesh%edges(2,i)) = nbnd_edges(mesh%edges(2,i)) + 1
       END IF
    END DO
    ! single-arc coastal nodes get the node itself as an extra vertex and are
    ! walked around element by element below; multi-arc nodes become a star
    ! polygon without a node vertex
    node_is_boundary = (nbnd_edges == 2)
    nbr_boundary_nodes = COUNT(node_is_boundary)

    nbr_vertices = partit%myDim_elem2D + partit%myDim_edge2D + nbr_boundary_nodes
    nbr_connections = 2*SUM(mesh%nod_in_elem2D_num(1:partit%myDim_nod2D)) + 2*nbr_boundary_nodes
    DO i=1,partit%myDim_nod2D
       IF (nbnd_edges(i) > 2) nbr_connections = nbr_connections + nbnd_edges(i)/2
    END DO
    IF (COUNT(nbnd_edges > 2) > 0) WRITE(*,*) 'YAC grid: rank ', partit%mype, ' builds ', &
         COUNT(nbnd_edges > 2), ' star cell(s) for nodes touching the coast in several places'
    ALLOCATE(star_vtx(2*MAXVAL(mesh%nod_in_elem2D_num(1:partit%myDim_nod2D)) + 16), &
             star_ang(2*MAXVAL(mesh%nod_in_elem2D_num(1:partit%myDim_nod2D)) + 16))

    ALLOCATE(x_vertices(nbr_vertices))
    ALLOCATE(y_vertices(nbr_vertices))
    ALLOCATE(nbr_vertices_per_cell(partit%myDim_nod2D))
    ALLOCATE(cell_to_vertex(nbr_connections))

    ! compute vertices
    ! 1 element centers
    DO i=1,partit%myDim_elem2D
       CALL compute_center(mesh%geo_coord_nod2D(1:2,mesh%elem2D_nodes(1, i)), &
            mesh%geo_coord_nod2D(1:2,mesh%elem2D_nodes(2, i)), &
            mesh%geo_coord_nod2D(1:2,mesh%elem2D_nodes(3, i)), mid)
       x_vertices(i) = mid(1)
       y_vertices(i) = mid(2)
    END DO
    ! 2 edges midpoints
    DO i=1,partit%myDim_edge2D
       CALL compute_midpoint(mesh%geo_coord_nod2D(1:2,mesh%edges(1,i)), mesh%geo_coord_nod2D(1:2,mesh%edges(2,i)), mid)
       x_vertices(partit%myDim_elem2D + i) = mid(1)
       y_vertices(partit%myDim_elem2D + i) = mid(2)
    END DO
    ! 3 boundary nodes -> are added on the fly when the cells are computed

    ! compute cells
    vtx_idx = partit%myDim_elem2D + partit%myDim_edge2D
    c2v_idx = 0
    DO i=1,partit%myDim_nod2D
       nbr_vertices_per_cell(i) = 2*mesh%nod_in_elem2D_num(i)
       IF (nbnd_edges(i) > 2) THEN
          ! star polygon: all element centres and edge midpoints around the
          ! node, ordered by azimuth. The element fan of such a node is not
          ! connected, so the edge walk below cannot close the polygon; this one
          ! spans the land wedges between the coastal arcs (a small over-estimate
          ! of the cell area at these few nodes; the field points stay the node
          ! coordinates, so nearest-neighbour exchanges are unaffected).
          n_star = 0
          DO j=1,mesh%nod_in_elem2D_num(i)
             n_star = n_star + 1
             star_vtx(n_star) = mesh%nod_in_elem2D(j,i)
             DO k=1,3
                ed = mesh%elem_edges(k, mesh%nod_in_elem2D(j,i))
                IF (ALL(mesh%edges(1:2,ed) /= i)) CYCLE
                IF (ANY(star_vtx(1:n_star) == partit%myDim_elem2D + ed)) CYCLE
                n_star = n_star + 1
                star_vtx(n_star) = partit%myDim_elem2D + ed
             END DO
          END DO
          DO l=1,n_star
             dx = x_vertices(star_vtx(l)) - mesh%geo_coord_nod2D(1,i)
             IF (dx >  PI) dx = dx - 2.0_WP*PI
             IF (dx < -PI) dx = dx + 2.0_WP*PI
             dx = dx*COS(mesh%geo_coord_nod2D(2,i))
             dy = y_vertices(star_vtx(l)) - mesh%geo_coord_nod2D(2,i)
             star_ang(l) = ATAN2(dy, dx)
          END DO
          DO l=2,n_star                      ! insertion sort by azimuth
             atmp = star_ang(l); itmp = star_vtx(l); k = l-1
             DO WHILE (k >= 1)
                IF (star_ang(k) <= atmp) EXIT
                star_ang(k+1) = star_ang(k); star_vtx(k+1) = star_vtx(k); k = k-1
             END DO
             star_ang(k+1) = atmp; star_vtx(k+1) = itmp
          END DO
          nbr_vertices_per_cell(i) = n_star
          cell_to_vertex(c2v_idx+1:c2v_idx+n_star) = star_vtx(1:n_star)
          c2v_idx = c2v_idx + n_star
          CYCLE
       END IF
       curr_elem = mesh%nod_in_elem2D(1,i)
       curr_edge = mesh%elem_edges(1, curr_elem)
       IF (ALL(mesh%edges(1:2,curr_edge) /= i)) curr_edge = mesh%elem_edges(2, curr_elem)
       IF (node_is_boundary(i)) THEN
          nbr_vertices_per_cell(i) = nbr_vertices_per_cell(i) + 2
          ! we're starting with the boundary node
          vtx_idx = vtx_idx+1
          x_vertices(vtx_idx) = mesh%geo_coord_nod2D(1,i)
          y_vertices(vtx_idx) = mesh%geo_coord_nod2D(2,i)
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = vtx_idx
          elem_loop: DO j=1,mesh%nod_in_elem2D_num(i)
             DO k=1,3
                IF (partit%myList_edge2D(mesh%elem_edges(k,mesh%nod_in_elem2D(j,i))) > mesh%edge2D_in .AND. &
                     ANY(mesh%edges(1:2,mesh%elem_edges(k,mesh%nod_in_elem2D(j,i))) == i)) THEN
                   curr_elem = mesh%nod_in_elem2D(j,i)
                   curr_edge = mesh%elem_edges(k, curr_elem)
                   EXIT elem_loop
                END IF
             END DO
          END DO elem_loop
       END IF

       DO j=1,mesh%nod_in_elem2D_num(i)
          IF (curr_elem == 0) THEN
             WRITE (0,*) "Error: cells with two or more coast lines are currently not supported by the yac setup code"
             STOP
          ENDIF
          ! add the midpoint of curr_edge
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = partit%myDim_elem2D + curr_edge
          ! add the center of curr_elem
          c2v_idx = c2v_idx + 1
          IF (curr_elem < 1) THEN
             WRITE (0,*) "elem 0 detected in cell ", i, " with ", nbr_vertices_per_cell(i), " vertices - j is ", j, " is_boundary ", node_is_boundary(i)
          END IF
          cell_to_vertex(c2v_idx) = curr_elem
          ! find next edge
          edge_loop: DO k=1,3
             IF (mesh%elem_edges(k, curr_elem) /= curr_edge .AND. &
                  ANY(mesh%edges(1:2,mesh%elem_edges(k, curr_elem)) == i)) THEN
                curr_edge = mesh%elem_edges(k, curr_elem)
                EXIT edge_loop
             END IF
          END DO edge_loop
          IF (mesh%edge_tri(1, curr_edge) /= curr_elem) THEN
             curr_elem = mesh%edge_tri(1, curr_edge)
          ELSE
             curr_elem = mesh%edge_tri(2, curr_edge)
          END IF
       END DO

       ! in case of a boundary volume, we add the others edges midpoint
       IF (node_is_boundary(i)) THEN
          ! add the midpoint of curr_edge
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = partit%myDim_elem2D + curr_edge
       END IF
    END DO

    DEALLOCATE(star_vtx, star_ang, nbnd_edges)

    CALL yac_fdef_grid( &
         trim(cpl_grid_name), &
         nbr_vertices, &   ! nbr_vertices
         partit%myDim_nod2D, &  ! nbr_cells
         nbr_connections, & ! nbr_connections
         nbr_vertices_per_cell, & ! nbr_vertices_per_cell
         x_vertices, & ! x_vertices
         y_vertices, & ! y_vertices
         cell_to_vertex, &
         grid_id)

    CALL yac_fset_global_index( &
         partit%myList_nod2D - 1,   &
         YAC_LOCATION_CELL,     &
         grid_id )

    CALL yac_fdef_points(grid_id, &
         partit%myDim_nod2D, &
         YAC_LOCATION_CELL, &
         mesh%geo_coord_nod2D(1,1:partit%myDim_nod2D), &
         mesh%geo_coord_nod2D(2,1:partit%myDim_nod2D), &
         points_id)

    ! The field names and collection sizes were registered by
    ! define_cpl_fields (cpl_fields). They must match coupling.yaml.
    ALLOCATE(send_field_id(nsend), recv_field_id(nrecv))

    ! Each field is registered with its own exchange interval (cpl_fields).
    DO i=1,nsend
       WRITE(dt_str, '(I8.8)') INT(dt*cpl_send_steps(i)*1000)
       CALL yac_fdef_field(cpl_send(i), comp_id, [points_id], 1, &
            cpl_send_collection_size(i), &
            dt_str, YAC_TIME_UNIT_MILLISECOND, send_field_id(i))
    END DO

    DO i=1,nrecv
       WRITE(dt_str, '(I8.8)') INT(dt*cpl_recv_steps(i)*1000)
       CALL yac_fdef_field(cpl_recv(i), comp_id, [points_id], 1, cpl_recv_collection_size(i), &
            dt_str, YAC_TIME_UNIT_MILLISECOND, recv_field_id(i))
    END DO

    CALL yac_fenddef()

  end subroutine cpl_yac_define_unstr


  subroutine cpl_yac_send(ind, data_array, action)
    implicit none

    integer, intent( IN )          :: ind       ! variable Id
    logical, intent( OUT )         :: action    !
    real(kind=WP),  intent(IN)     :: data_array(:, :)
    integer :: info, ierr

    call yac_fput(send_field_id(ind), SIZE(data_array, 1), SIZE(data_array, 2), &
         data_array, info, ierr)
    action = info .eq. YAC_ACTION_COUPLING
    if (info .eq. YAC_ACTION_OUT_OF_BOUND) call out_of_bound_abort('put', cpl_send(ind))
    if (cpl_send_partner(ind) /= CPL_PARTNER_ATMOSPHERE) then
       n_put_ice = n_put_ice + 1
       if (action) n_put_ice_coupled = n_put_ice_coupled + 1
       call check_coupling_period()
    end if
  end subroutine cpl_yac_send

  subroutine cpl_yac_recv(ind, data_array, action)
    implicit none

    integer, intent( IN )  :: ind       ! variable Id
    logical, intent( OUT ) :: action    ! 
    real(kind=WP), intent( INOUT )    :: data_array(:,:)
    integer :: info, ierr

    call yac_fget(recv_field_id(ind), SIZE(data_array, 1), SIZE(data_array, 2),&
         data_array, info, ierr)
    action = info .eq. YAC_ACTION_COUPLING
    if (info .eq. YAC_ACTION_OUT_OF_BOUND) call out_of_bound_abort('get', cpl_recv(ind))
  end subroutine cpl_yac_recv

  !> End of the run: the clock walked forward by nsteps*dt seconds, as YAC
  !> datetime string.
  subroutine run_end_datetime(nsteps, enddatetime)
    use g_config, only: dt
    integer,           intent(in)  :: nsteps
    character(len=24), intent(out) :: enddatetime
    integer       :: eyear, eday, emonth, eday_in_month, fleap
    real(kind=WP) :: esec

    eyear = yearnew
    eday  = daynew
    esec  = timenew + real(nsteps, WP)*dt
    call check_fleapyr(eyear, fleap)
    do while (esec >= 86400.0_WP)
       esec = esec - 86400.0_WP
       eday = eday + 1
       if (eday > 365 + fleap) then
          eday  = 1
          eyear = eyear + 1
          call check_fleapyr(eyear, fleap)
       end if
    end do
    eday_in_month = eday
    do emonth = 1, 12
       if (eday_in_month <= num_day_in_month(fleap, emonth)) exit
       eday_in_month = eday_in_month - num_day_in_month(fleap, emonth)
    end do
    WRITE(enddatetime, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,".",I3.3, "Z")') &
         eyear, emonth, eday_in_month, &
         INT(esec)/3600, MODULO(INT(esec), 3600)/60, MODULO(INT(esec), 60), &
         INT(MODULO(esec, 1.0_WP)*1000)
  end subroutine run_end_datetime

  !> YAC answers a put/get outside the run period with YAC_ACTION_OUT_OF_BOUND
  !> and only warns; the component then runs on with stale or zero partner
  !> fields. A run that has left its coupling window is not valid: stop.
  subroutine out_of_bound_abort(what, name)
    use mpi
    character(len=*), intent(in) :: what, name
    integer :: ierr
    write(*,'(a)') 'FATAL: YAC '//what//' of '//trim(name)//' returned YAC_ACTION_OUT_OF_BOUND:'
    write(*,'(a)') '       the model clock is outside the run period known to YAC.'
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end subroutine out_of_bound_abort

  !> Fires once, some 40 sea-ice puts into the run, on rank 0 only.
  subroutine check_coupling_period()
    if (period_checked .or. n_put_ice < 40) return
    period_checked = .true.
    if (driver_mype /= 0) return
    if (n_put_ice_coupled >= n_put_ice - 4) return   ! ~1:1, allowing for src/tgt_lag
    write(*,*) 'WARNING: coupling_period in coupling.yaml does not match dt*cpl_stride:'
    write(*,*) '  of ', n_put_ice, ' ocean/sea-ice puts only ', n_put_ice_coupled, ' coupled.'
    write(*,*) '  The components then reuse the fields of the last coupling in between'
    write(*,*) '  (zero-order hold), a known cause of surface run-aways. Set'
    write(*,*) '  coupling_period to dt*cpl_stride, or average over the period with'
    write(*,*) '  "time_reduction: average". Each component reads the coupling.yaml'
    write(*,*) '  of its own run directory.'
  end subroutine check_coupling_period

  subroutine cpl_yac_finalize ()
    implicit none
#ifdef VERBOSE
    print *, '================================================='
    print *, 'cpl_yac_finalize : coupler finalization for YAC'
    print *, '*************************************************'
#endif /* VERBOSE */
    CALL yac_ffinalize()
  end subroutine cpl_yac_finalize

#endif
end module cpl_yac_driver

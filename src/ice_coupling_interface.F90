! Sea-ice coupling interface (FESOM ocean side, YAC).
!
! The seam through which FESOM exchanges fields with FESIM, the sea ice
! running as its own YAC component (is_coupled_to_fesim). The sibling module
! on the FESIM side is ocean_coupling_interface.
!
! The native ocean state sent to the ice is the same for every atmosphere.
! What the ocean forwards from the atmosphere follows the runtime partner
! selection in namelist.cpl (cpl_config), so the send layout is decided in
! ice_cpl_define rather than by the preprocessor:
!
!   is_coupled_to_icon_a: the ocean receives pre-computed atmospheric FLUXES
!     from ICON over YAC and forwards them. Sends:
!       native:    sst_feom_to_ice, ocean_to_ice_bundle, ocean_to_ice_uv
!       forwarded: taux_to_ice, tauy_to_ice, surface_fresh_water_flux_to_ice,
!                  total_heat_flux_to_ice, atmosphere_sea_ice_bundle_to_ice
!
!   is_coupled_to_ifs: the atmosphere reaches the ocean through the IFS
!     interface, which deposits fluxes already in FESOM units (precip/evap in
!     m/s, stresses in Pa and rotated to the FESOM grid, ice-fraction
!     weighted). They are forwarded 1:1 as one 10-component bundle; FESIM does
!     not re-convert or re-rotate. Components:
!       1 stress_atmice_x [Pa]   2 stress_atmice_y [Pa]
!       3 oce_heat_flux [W/m2]   4 ice_heat_flux [W/m2]   5 shortwave [W/m2]
!       6 prec_rain [m/s]        7 prec_snow [m/s]         8 evap_no_ifrac [m/s]
!       9 sublimation [m/s]     10 enthalpyoffuse [W/m2]
!     stress_atmoce is not forwarded: it drives ocean momentum, which lives
!     here. The sea-ice bundle received from FESIM then carries ice surface
!     temperature and albedo as well, for the relay back to IFS.
!
!   no atmosphere (FESOM + FESIM alone): the ocean reads its forcing files
!     and forwards the raw atmospheric STATE as one 9-component bundle; FESIM
!     computes the bulk fluxes itself. Components (forcing-file units):
!       1 u_wind [m/s]   2 v_wind [m/s]   3 t_air [K]    4 shum [kg/kg]
!       5 shortwave down [W/m2]           6 longwave down [W/m2]
!       7 prec_rain [kg/m2/s]             8 prec_snow [kg/m2/s]   9 mslp [Pa]
!
!   recvs from ice (every partner): sea_ice_bundle, ice_to_ocean_stress,
!     ice_to_ocean_flux
!
! river_runoff is not forwarded: it is ocean-only. Forwarded fields travel in
! their source units; FESIM converts on receive.
!
! `sst_feom_to_ice` duplicates the atmosphere-bound `sst_feom` so that each
! interface owns its YAC fields end to end.
!
! ice_cpl_define takes the field time step in seconds. FESIM registers its
! own time step on its side; the coupling_period in coupling.yaml has to
! equal dt*cpl_stride on both sides (cpl_check_period warns when it does
! not).
module ice_coupling_interface
#if defined (__cpl_yac)

  use yac
  use o_PARAM, only: WP

  implicit none
  private

  ! --- public field-index parameters ------------------------------------
  ! Send slots (1-based). Native ocean state first, then the forwarded
  ! atmospheric fields, whose slots depend on the atmosphere partner. The
  ! three layouts reuse slot 4 onwards, so a slot constant is only
  ! meaningful under the partner it belongs to.
  integer, parameter, public :: ICE_SEND_SST_FEOM             = 1
  integer, parameter, public :: ICE_SEND_OCEAN_TO_ICE_BUNDLE  = 2
  integer, parameter, public :: ICE_SEND_OCEAN_TO_ICE_UV      = 3
  ! is_coupled_to_icon_a: forwarded ICON fluxes
  integer, parameter, public :: ICE_SEND_TAUX                 = 4
  integer, parameter, public :: ICE_SEND_TAUY                 = 5
  integer, parameter, public :: ICE_SEND_FRESH_WATER          = 6
  integer, parameter, public :: ICE_SEND_HEAT_FLUX            = 7
  integer, parameter, public :: ICE_SEND_ATM_SEA_ICE_BUNDLE   = 8
  ! is_coupled_to_ifs: forwarded IFS fluxes
  integer, parameter, public :: ICE_SEND_ATM_ICE_FLUX         = 4
  ! no atmosphere: raw atmospheric state from the forcing files
  integer, parameter, public :: ICE_SEND_ATM_STATE            = 4
  integer, parameter, public :: ICE_NSEND_MAX                 = 8

  ! Recv slots (1-based), the same for every partner: ice state, the
  ! ice->ocean momentum drag, and FESIM's net heat + freshwater flux.
  integer, parameter, public :: ICE_RECV_SEA_ICE_BUNDLE       = 1
  integer, parameter, public :: ICE_RECV_ICE_STRESS           = 2
  integer, parameter, public :: ICE_RECV_ICE_FLUX             = 3
  integer, parameter, public :: ICE_NRECV                     = 3

  ! --- runtime layout, set by ice_cpl_define from the partner selection ---
  integer,           public, protected, save :: ICE_NSEND = 0
  integer,           public, protected, save :: ice_send_collection_size(ICE_NSEND_MAX) = 0
  character(len=32), public, protected, save :: ice_send_names(ICE_NSEND_MAX) = ''
  ! sea_ice_bundle: 3 (m_ice, m_snow, a_ice), or 5 with ice_temp + ice_alb
  ! under IFS; ice_to_ocean_stress: 2; ice_to_ocean_flux: 2.
  integer,           public, protected, save :: ice_recv_collection_size(ICE_NRECV) = [3, 2, 2]
  character(len=32), parameter, public :: ice_recv_names(ICE_NRECV) = [character(len=32) :: &
       'sea_ice_bundle', &
       'ice_to_ocean_stress', &
       'ice_to_ocean_flux' ]

  ! --- private module state ---------------------------------------------
  integer, save :: ice_send_field_id(ICE_NSEND_MAX) = -1
  integer, save :: ice_recv_field_id(ICE_NRECV)     = -1
  integer, save :: ice_points_id_local              = -1

  ! --- coupling-cost instrumentation ------------------------------------
  ! yac_fput/yac_fget are the only points at which the pair exchanges, so
  ! timing them here attributes the coupler cost without touching callers.
  ! The remainder of the caller's bracket (cpl_call_tic/toc) is then halo
  ! exchanges, unit conversion and copies.
  real(kind=WP), public, save :: cpl_time_put = 0.0_WP   ! s, this rank, cumulative
  real(kind=WP), public, save :: cpl_time_get = 0.0_WP
  integer,       public, save :: cpl_n_put    = 0        ! calls (not exchanges:
  integer,       public, save :: cpl_n_get    = 0        !  YAC no-ops off-period)
  integer,       public, save :: cpl_n_put_act = 0       ! calls that actually coupled
  integer,       public, save :: cpl_n_get_act = 0
  real(kind=WP), public, save :: cpl_time_call = 0.0_WP  ! whole exchange routine
  ! One-shot check that coupling_period equals dt*cpl_stride. The Fortran
  ! side registers its fields with dt*cpl_stride, but coupling_period lives
  ! in coupling.yaml and nothing enforces agreement. A larger yaml period
  ! makes YAC couple only every Nth call; the other calls return no-action
  ! and the component silently reuses stale fields. That zero-order hold
  ! produced surface runaways within ~100 steps on high-resolution meshes.
  integer, save :: cpl_mype    = -1
  logical, save :: cpl_checked = .false.

  public :: ice_cpl_init, ice_cpl_define, ice_cpl_send, ice_cpl_recv, ice_cpl_finalize
  public :: cpl_timers_report, cpl_call_tic, cpl_call_toc

contains

  subroutine ice_cpl_init(localCommunicator)
    use yac_component_runtime, only: yac_runtime_init
    integer, intent(out) :: localCommunicator
    call yac_runtime_init(localCommunicator)
  end subroutine ice_cpl_init

  ! Pick the send/recv layout from the atmosphere partner (namelist.cpl).
  subroutine ice_cpl_set_layout()
    use cpl_config, only: is_coupled_to_icon_a, is_coupled_to_ifs

    ice_send_names            = ''
    ice_send_collection_size  = 0

    ice_send_names(1:3) = [character(len=32) :: &
         'sst_feom_to_ice', &
         'ocean_to_ice_bundle', &
         'ocean_to_ice_uv' ]
    ice_send_collection_size(1:3) = [1, 2, 2]

    if (is_coupled_to_icon_a) then
       ICE_NSEND = 8
       ice_send_names(4:8) = [character(len=32) :: &
            'taux_to_ice', &
            'tauy_to_ice', &
            'surface_fresh_water_flux_to_ice', &
            'total_heat_flux_to_ice', &
            'atmosphere_sea_ice_bundle_to_ice' ]
       ! mirror atm_recv_collection_size (taux, tauy, fresh_water, heat_flux,
       ! atm_sea_ice_bundle)
       ice_send_collection_size(4:8) = [2, 2, 3, 4, 2]
    else if (is_coupled_to_ifs) then
       ICE_NSEND = 4
       ice_send_names(4)           = 'atm_ice_flux_to_ice'
       ice_send_collection_size(4) = 10
    else
       ICE_NSEND = 4
       ice_send_names(4)           = 'atm_state_to_ice'
       ice_send_collection_size(4) = 9
    end if

    ice_recv_collection_size = [3, 2, 2]
    if (is_coupled_to_ifs) ice_recv_collection_size(ICE_RECV_SEA_ICE_BUNDLE) = 5
  end subroutine ice_cpl_set_layout

  ! dt_seconds: the time step the fields are registered with (in seconds;
  ! fractional values are kept, YAC gets milliseconds).
  subroutine ice_cpl_define(partit, mesh, dt_seconds)
    use MOD_MESH,              only: t_mesh
    use MOD_PARTIT,            only: t_partit
    use cpl_config,            only: cpl_grid_name
    use yac_component_runtime, only: yac_runtime_ensure_grid, yac_runtime_comp_id
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP),  intent(in)            :: dt_seconds

    character(len=8) :: dt_str
    integer          :: grid_id, i

    cpl_mype = partit%mype

    call ice_cpl_set_layout()

    call yac_runtime_ensure_grid(trim(cpl_grid_name), partit, mesh, grid_id, ice_points_id_local)

    write(dt_str, '(I8.8)') INT(dt_seconds*1000)

    do i = 1, ICE_NSEND
       call yac_fdef_field(ice_send_names(i), yac_runtime_comp_id(), &
            [ice_points_id_local], 1, ice_send_collection_size(i), &
            dt_str, YAC_TIME_UNIT_MILLISECOND, ice_send_field_id(i))
    end do

    do i = 1, ICE_NRECV
       call yac_fdef_field(ice_recv_names(i), yac_runtime_comp_id(), &
            [ice_points_id_local], 1, ice_recv_collection_size(i), &
            dt_str, YAC_TIME_UNIT_MILLISECOND, ice_recv_field_id(i))
    end do
  end subroutine ice_cpl_define

  subroutine ice_cpl_send(ind, data_array, action)
    use mpi, only: MPI_Wtime
    integer,       intent(in)  :: ind
    real(kind=WP), intent(in)  :: data_array(:,:)
    logical,       intent(out) :: action
    integer :: info, ierr
    real(kind=WP) :: t_cpl0
    t_cpl0 = MPI_Wtime()
    call yac_fput(ice_send_field_id(ind), size(data_array, 1), size(data_array, 2), &
         data_array, info, ierr)
    cpl_time_put = cpl_time_put + (MPI_Wtime() - t_cpl0)
    cpl_n_put = cpl_n_put + 1
    action = info == YAC_ACTION_COUPLING
    if (action) cpl_n_put_act = cpl_n_put_act + 1
    call cpl_check_period()
  end subroutine ice_cpl_send

  subroutine ice_cpl_recv(ind, data_array, action)
    use mpi, only: MPI_Wtime
    integer,       intent(in)    :: ind
    real(kind=WP), intent(inout) :: data_array(:,:)
    logical,       intent(out)   :: action
    integer :: info, ierr
    real(kind=WP) :: t_cpl0
    t_cpl0 = MPI_Wtime()
    call yac_fget(ice_recv_field_id(ind), size(data_array, 1), size(data_array, 2), &
         data_array, info, ierr)
    cpl_time_get = cpl_time_get + (MPI_Wtime() - t_cpl0)
    cpl_n_get = cpl_n_get + 1
    action = info == YAC_ACTION_COUPLING
    if (action) cpl_n_get_act = cpl_n_get_act + 1
  end subroutine ice_cpl_recv

  subroutine ice_cpl_finalize()
    use yac_component_runtime, only: yac_runtime_finalize
    call yac_runtime_finalize()
  end subroutine ice_cpl_finalize

  ! Bracket the caller's whole exchange routine (tic at entry, toc at exit).
  subroutine cpl_call_tic(t0)
    use mpi, only: MPI_Wtime
    real(kind=WP), intent(out) :: t0
    t0 = MPI_Wtime()
  end subroutine cpl_call_tic

  subroutine cpl_call_toc(t0)
    use mpi, only: MPI_Wtime
    real(kind=WP), intent(in) :: t0
    cpl_time_call = cpl_time_call + (MPI_Wtime() - t0)
  end subroutine cpl_call_toc

  ! Collective report of the coupler cost, printed next to the model's own
  ! per-task runtime block: mean/min/max over ranks.
  subroutine cpl_timers_report(comm, mype, npes, label)
    use mpi
    integer,          intent(in) :: comm, mype, npes
    character(len=*), intent(in) :: label
    ! 1 fput, 2 fget, 3 whole exchange routine, 4 the routine minus YAC.
    ! (4) is formed per rank before reducing: reducing (3) and (1)+(2)
    ! separately and subtracting would not give the min/max of the difference.
    real(kind=WP) :: v(4), vsum(4), vmin(4), vmax(4)
    integer       :: c(4), csum(4), ierr

    v(1) = cpl_time_put
    v(2) = cpl_time_get
    v(3) = cpl_time_call
    v(4) = cpl_time_call - cpl_time_put - cpl_time_get
    c = [cpl_n_put, cpl_n_get, cpl_n_put_act, cpl_n_get_act]
    call MPI_Allreduce(v, vsum, 4, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    call MPI_Allreduce(v, vmin, 4, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    call MPI_Allreduce(v, vmax, 4, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    call MPI_Allreduce(c, csum, 4, MPI_INTEGER,          MPI_SUM, comm, ierr)
    if (mype /= 0) return
    vsum = vsum / real(npes, WP)
    print '(a)',        '___COUPLER COST ('//label//') per task [seconds]__mean___________min___________max_'
    print '(a,3f14.4)', '   yac_fput                  :', vsum(1), vmin(1), vmax(1)
    print '(a,3f14.4)', '   yac_fget                  :', vsum(2), vmin(2), vmax(2)
    print '(a,3f14.4)', '   yac_fput + yac_fget       :', vsum(1)+vsum(2), vmin(1)+vmin(2), vmax(1)+vmax(2)
    if (cpl_time_call > 0.0_WP) then
       print '(a,3f14.4)', '   exchange routine total    :', vsum(3), vmin(3), vmax(3)
       print '(a,3f14.4)', '   ... of which NOT yac      :', vsum(4), vmin(4), vmax(4)
    end if
    print '(a,2i12)',   '   fput calls / of which coupling :', csum(1)/npes, csum(3)/npes
    print '(a,2i12)',   '   fget calls / of which coupling :', csum(2)/npes, csum(4)/npes
  end subroutine cpl_timers_report

  ! Fires once, ~10 coupling calls in, on rank 0 only.
  subroutine cpl_check_period()
    if (cpl_checked .or. cpl_n_put < 40) return
    cpl_checked = .true.
    if (cpl_mype /= 0) return
    if (cpl_n_put_act >= cpl_n_put - 4) return   ! ~1:1 allowing for src/tgt_lag
    write(*,*) '**********************************************************************'
    write(*,*) '*  WARNING: coupling_period does NOT match dt * cpl_stride            *'
    write(*,*) '**********************************************************************'
    write(*,*) '  yac_fput calls so far      :', cpl_n_put
    write(*,*) '  ... of which really coupled:', cpl_n_put_act
    write(*,*) '  Expected these to be ~equal. They are not, so YAC is coupling only'
    write(*,*) '  every Nth call and this component is reusing STALE fields in'
    write(*,*) '  between (zero-order hold). That is a known cause of surface'
    write(*,*) '  runaways at high resolution.'
    write(*,*) '  FIX: set coupling_period in coupling.yaml equal to dt*cpl_stride,'
    write(*,*) '  or, for genuine asynchronous coupling, add "time_reduction:'
    write(*,*) '  average" to the couple so the source averages over the period'
    write(*,*) '  instead of sending an instantaneous snapshot.'
    write(*,*) '  NB each component reads the coupling.yaml in ITS OWN directory.'
    write(*,*) '**********************************************************************'
  end subroutine cpl_check_period

#endif
end module ice_coupling_interface

! Ocean coupling interface (FESIM sea-ice side, YAC).
!
! The only coupling seam of the sea-ice component: FESIM talks to the ocean
! and never to the atmosphere directly. The sibling module on the FESOM side
! is ice_coupling_interface. Both register on the grid built by
! yac_grid_utils; component, grid and configuration-file names come from
! &coupling_yac in namelist.cpl (cpl_config).
!
! What the ocean forwards from the atmosphere follows the runtime partner
! selection in namelist.cpl -- the ocean's atmosphere, read by this
! executable as well -- so the receive layout is decided in ocn_cpl_define
! rather than by the preprocessor. It must match the send layout of the
! ocean's ice_coupling_interface:
!
!   is_coupled_to_icon_a: pre-computed ICON FLUXES, forwarded by the ocean.
!     recvs: sst_feom_to_ice, taux_to_ice, tauy_to_ice,
!            surface_fresh_water_flux_to_ice, total_heat_flux_to_ice,
!            atmosphere_sea_ice_bundle_to_ice, ocean_to_ice_bundle,
!            ocean_to_ice_uv
!   is_coupled_to_ifs: the IFS fluxes as one 10-component bundle, already in
!     FESOM units and rotated (no conversion on receive).
!     recvs: sst_feom_to_ice, atm_ice_flux_to_ice, ocean_to_ice_bundle,
!            ocean_to_ice_uv
!     The sea-ice bundle sent up then also carries ice surface temperature
!     and albedo (5 components instead of 3) for the relay to IFS.
!   no atmosphere: the raw atmospheric STATE from the ocean's forcing files
!     as one 9-component bundle; FESIM computes the bulk fluxes itself.
!     recvs: sst_feom_to_ice, atm_state_to_ice, ocean_to_ice_bundle,
!            ocean_to_ice_uv
!
!   sends (every partner): sea_ice_bundle, ice_to_ocean_stress,
!     ice_to_ocean_flux
!
! ocn_cpl_define takes the field time step in seconds; the coupling_period
! in coupling.yaml has to equal dt*cpl_stride on both sides
! (cpl_check_period warns when it does not).
module ocean_coupling_interface
#if defined (__cpl_yac)

  use yac
  use o_PARAM, only: WP

  implicit none
  private

  ! --- public field-index parameters ------------------------------------
  ! Send slots (1-based), the same for every partner.
  integer, parameter, public :: OCN_SEND_SEA_ICE_BUNDLE       = 1
  integer, parameter, public :: OCN_SEND_ICE_STRESS           = 2
  integer, parameter, public :: OCN_SEND_ICE_FLUX             = 3
  integer, parameter, public :: OCN_NSEND                     = 3

  ! Recv slots (1-based). SST first, then the forwarded atmospheric fields,
  ! whose slots depend on the atmosphere partner, then the native ocean
  ! state. The three layouts reuse slot 2 onwards, so a slot constant is
  ! only meaningful under the partner it belongs to.
  integer, parameter, public :: OCN_RECV_SST_FEOM             = 1
  ! is_coupled_to_icon_a
  integer, parameter, public :: OCN_RECV_TAUX                 = 2
  integer, parameter, public :: OCN_RECV_TAUY                 = 3
  integer, parameter, public :: OCN_RECV_FRESH_WATER          = 4
  integer, parameter, public :: OCN_RECV_HEAT_FLUX            = 5
  integer, parameter, public :: OCN_RECV_ATM_SEA_ICE_BUNDLE   = 6
  ! is_coupled_to_ifs
  integer, parameter, public :: OCN_RECV_ATM_ICE_FLUX         = 2
  ! no atmosphere
  integer, parameter, public :: OCN_RECV_ATM_STATE            = 2
  ! native ocean state: slot depends on the partner, see ocn_cpl_set_layout
  integer,           public, protected, save :: OCN_RECV_OCEAN_TO_ICE_BUNDLE = 0
  integer,           public, protected, save :: OCN_RECV_OCEAN_TO_ICE_UV     = 0
  integer, parameter, public :: OCN_NRECV_MAX                 = 8

  ! --- runtime layout, set by ocn_cpl_define from the partner selection ---
  integer,           public, protected, save :: OCN_NRECV = 0
  integer,           public, protected, save :: ocn_recv_collection_size(OCN_NRECV_MAX) = 0
  character(len=32), public, protected, save :: ocn_recv_names(OCN_NRECV_MAX) = ''
  ! sea_ice_bundle: 3 (m_ice, m_snow, a_ice), or 5 with ice_temp + ice_alb
  ! under IFS; ice_to_ocean_stress: 2; ice_to_ocean_flux: 2.
  integer,           public, protected, save :: ocn_send_collection_size(OCN_NSEND) = [3, 2, 2]
  character(len=32), parameter, public :: ocn_send_names(OCN_NSEND) = [character(len=32) :: &
       'sea_ice_bundle', &
       'ice_to_ocean_stress', &
       'ice_to_ocean_flux' ]

  ! Kept for the flux-correction routines that still name these symbols
  ! (force_flux_consv, net_rec_from_atm). Never set under YAC.
  real(kind=WP), allocatable, public :: a2o_fcorr_stat(:,:)
  integer,                   public :: source_root  = 0
  integer,                   public :: target_root  = 0
  logical,                   public :: commRank     = .false.

  ! --- private module state ---------------------------------------------
  integer, save :: ocn_comp_id          = -1
  integer, save :: ocn_local_comm       = -1
  integer, save :: ocn_grid_id          = -1
  integer, save :: ocn_points_id        = -1
  integer, save :: ocn_send_field_id(OCN_NSEND)     = -1
  integer, save :: ocn_recv_field_id(OCN_NRECV_MAX) = -1
  logical, save :: ocn_inited           = .false.

  ! --- coupling-cost instrumentation ------------------------------------
  ! yac_fput/yac_fget are the only points at which the pair exchanges, so
  ! timing them here attributes the coupler cost without touching callers.
  real(kind=WP), public, save :: cpl_time_put = 0.0_WP   ! s, this rank, cumulative
  real(kind=WP), public, save :: cpl_time_get = 0.0_WP
  integer,       public, save :: cpl_n_put    = 0        ! calls (not exchanges:
  integer,       public, save :: cpl_n_get    = 0        !  YAC no-ops off-period)
  integer,       public, save :: cpl_n_put_act = 0       ! calls that actually coupled
  integer,       public, save :: cpl_n_get_act = 0
  ! One-shot check that coupling_period equals dt*cpl_stride: a larger yaml
  ! period makes YAC couple only every Nth call; the other calls return
  ! no-action and the component silently reuses stale fields. That
  ! zero-order hold produced surface runaways within ~100 steps on
  ! high-resolution meshes.
  integer, save :: cpl_mype    = -1
  logical, save :: cpl_checked = .false.

  public :: ocn_cpl_init, ocn_cpl_define, ocn_cpl_send, ocn_cpl_recv, ocn_cpl_finalize
  public :: cpl_timers_report

contains

  subroutine ocn_cpl_init(localCommunicator)
    use cpl_config, only: cpl_comp_name, cpl_config_file
    integer, intent(out) :: localCommunicator
    if (.not. ocn_inited) then
#ifdef VERBOSE
       print *, '================================================='
       print *, 'ocn_cpl_init : coupler initialization for YAC'
       print *, '*************************************************'
#endif
       call yac_finit()
       call yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
       call yac_fread_config_yaml(trim(cpl_config_file))
       call yac_fdef_comp(trim(cpl_comp_name), ocn_comp_id)
       call yac_fget_comp_comm(ocn_comp_id, ocn_local_comm)
       ocn_inited = .true.
    end if
    localCommunicator = ocn_local_comm
  end subroutine ocn_cpl_init

  ! Pick the recv/send layout from the ocean's atmosphere (namelist.cpl).
  subroutine ocn_cpl_set_layout()
    use cpl_config, only: is_coupled_to_icon_a, is_coupled_to_ifs

    ocn_recv_names           = ''
    ocn_recv_collection_size = 0

    ocn_recv_names(1)           = 'sst_feom_to_ice'
    ocn_recv_collection_size(1) = 1

    if (is_coupled_to_icon_a) then
       ocn_recv_names(2:6) = [character(len=32) :: &
            'taux_to_ice', &
            'tauy_to_ice', &
            'surface_fresh_water_flux_to_ice', &
            'total_heat_flux_to_ice', &
            'atmosphere_sea_ice_bundle_to_ice' ]
       ocn_recv_collection_size(2:6) = [2, 2, 3, 4, 2]
       OCN_RECV_OCEAN_TO_ICE_BUNDLE = 7
       OCN_RECV_OCEAN_TO_ICE_UV     = 8
    else if (is_coupled_to_ifs) then
       ocn_recv_names(2)           = 'atm_ice_flux_to_ice'
       ocn_recv_collection_size(2) = 10
       OCN_RECV_OCEAN_TO_ICE_BUNDLE = 3
       OCN_RECV_OCEAN_TO_ICE_UV     = 4
    else
       ocn_recv_names(2)           = 'atm_state_to_ice'
       ocn_recv_collection_size(2) = 9
       OCN_RECV_OCEAN_TO_ICE_BUNDLE = 3
       OCN_RECV_OCEAN_TO_ICE_UV     = 4
    end if
    ocn_recv_names(OCN_RECV_OCEAN_TO_ICE_BUNDLE)           = 'ocean_to_ice_bundle'
    ocn_recv_collection_size(OCN_RECV_OCEAN_TO_ICE_BUNDLE) = 2
    ocn_recv_names(OCN_RECV_OCEAN_TO_ICE_UV)               = 'ocean_to_ice_uv'
    ocn_recv_collection_size(OCN_RECV_OCEAN_TO_ICE_UV)     = 2
    OCN_NRECV = OCN_RECV_OCEAN_TO_ICE_UV

    ocn_send_collection_size = [3, 2, 2]
    if (is_coupled_to_ifs) ocn_send_collection_size(OCN_SEND_SEA_ICE_BUNDLE) = 5
  end subroutine ocn_cpl_set_layout

  ! dt_seconds: the time step the fields are registered with (in seconds;
  ! fractional values are kept, YAC gets milliseconds).
  subroutine ocn_cpl_define(partit, mesh, dt_seconds)
    use MOD_MESH,       only: t_mesh
    use MOD_PARTIT,     only: t_partit
    use cpl_config,     only: cpl_grid_name
    use yac_grid_utils, only: cpl_yac_define_unstr_generic
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP),  intent(in)            :: dt_seconds

    character(len=8) :: dt_str
    integer          :: i

    cpl_mype = partit%mype

    call ocn_cpl_set_layout()

    call cpl_yac_define_unstr_generic(partit, mesh, trim(cpl_grid_name), ocn_grid_id, ocn_points_id)

    write(dt_str, '(I8.8)') INT(dt_seconds*1000)

    do i = 1, OCN_NSEND
       call yac_fdef_field(ocn_send_names(i), ocn_comp_id, &
            [ocn_points_id], 1, ocn_send_collection_size(i), &
            dt_str, YAC_TIME_UNIT_MILLISECOND, ocn_send_field_id(i))
    end do

    do i = 1, OCN_NRECV
       call yac_fdef_field(ocn_recv_names(i), ocn_comp_id, &
            [ocn_points_id], 1, ocn_recv_collection_size(i), &
            dt_str, YAC_TIME_UNIT_MILLISECOND, ocn_recv_field_id(i))
    end do

    call yac_fenddef()
  end subroutine ocn_cpl_define

  subroutine ocn_cpl_send(ind, data_array, action)
    use mpi, only: MPI_Wtime
    integer,       intent(in)  :: ind
    real(kind=WP), intent(in)  :: data_array(:,:)
    logical,       intent(out) :: action
    integer :: info, ierr
    real(kind=WP) :: t_cpl0
    t_cpl0 = MPI_Wtime()
    call yac_fput(ocn_send_field_id(ind), size(data_array, 1), size(data_array, 2), &
         data_array, info, ierr)
    cpl_time_put = cpl_time_put + (MPI_Wtime() - t_cpl0)
    cpl_n_put = cpl_n_put + 1
    action = info == YAC_ACTION_COUPLING
    if (action) cpl_n_put_act = cpl_n_put_act + 1
    call cpl_check_period()
  end subroutine ocn_cpl_send

  subroutine ocn_cpl_recv(ind, data_array, action)
    use mpi, only: MPI_Wtime
    integer,       intent(in)    :: ind
    real(kind=WP), intent(inout) :: data_array(:,:)
    logical,       intent(out)   :: action
    integer :: info, ierr
    real(kind=WP) :: t_cpl0
    t_cpl0 = MPI_Wtime()
    call yac_fget(ocn_recv_field_id(ind), size(data_array, 1), size(data_array, 2), &
         data_array, info, ierr)
    cpl_time_get = cpl_time_get + (MPI_Wtime() - t_cpl0)
    cpl_n_get = cpl_n_get + 1
    action = info == YAC_ACTION_COUPLING
    if (action) cpl_n_get_act = cpl_n_get_act + 1
  end subroutine ocn_cpl_recv

  subroutine ocn_cpl_finalize()
    if (ocn_inited) then
#ifdef VERBOSE
       print *, '================================================='
       print *, 'ocn_cpl_finalize : coupler finalization for YAC'
       print *, '*************************************************'
#endif
       call yac_ffinalize()
       ocn_inited = .false.
    end if
  end subroutine ocn_cpl_finalize

  ! Collective report of the coupler cost, printed next to the model's own
  ! per-task runtime block: mean/min/max over ranks.
  subroutine cpl_timers_report(comm, mype, npes, label)
    use mpi
    integer,          intent(in) :: comm, mype, npes
    character(len=*), intent(in) :: label
    real(kind=WP) :: v(2), vsum(2), vmin(2), vmax(2)
    integer       :: c(4), csum(4), ierr

    v = [cpl_time_put, cpl_time_get]
    c = [cpl_n_put, cpl_n_get, cpl_n_put_act, cpl_n_get_act]
    call MPI_Allreduce(v, vsum, 2, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    call MPI_Allreduce(v, vmin, 2, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    call MPI_Allreduce(v, vmax, 2, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    call MPI_Allreduce(c, csum, 4, MPI_INTEGER,          MPI_SUM, comm, ierr)
    if (mype /= 0) return
    vsum = vsum / real(npes, WP)
    print '(a)',          '___COUPLER COST ('//label//') per task [seconds]_mean______min______max_'
    print '(a,3f14.4)',   '   yac_fput                  :', vsum(1), vmin(1), vmax(1)
    print '(a,3f14.4)',   '   yac_fget                  :', vsum(2), vmin(2), vmax(2)
    print '(a,3f14.4)',   '   yac_fput + yac_fget       :', vsum(1)+vsum(2), &
                                                            vmin(1)+vmin(2), vmax(1)+vmax(2)
    print '(a,2i12)',     '   fput calls / of which coupling :', csum(1)/npes, csum(3)/npes
    print '(a,2i12)',     '   fget calls / of which coupling :', csum(2)/npes, csum(4)/npes
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
end module ocean_coupling_interface

! Ocean coupling interface (FESIM sea-ice side).
!
! This module is the sole coupling seam on the FESIM side: sea-ice talks
! only to the ocean, never directly to the atmosphere. The backend is
! always YAC. The sibling module on the FESOM side is ice_coupling_interface.
!
! Field inventory (ocean side from the ice's perspective). ALL fields arrive
! from the ocean (FESOM) — FESIM never talks to the atmosphere directly. What
! the ocean forwards from the atmosphere depends on __yac_atm:
!   sends to ocean (both): sea_ice_bundle
!   __yac_atm (ICON, Stage 3): recv pre-computed atm FLUXES — sst_feom_to_ice,
!     taux_to_ice, tauy_to_ice, surface_fresh_water_flux_to_ice,
!     total_heat_flux_to_ice, atmosphere_sea_ice_bundle_to_ice,
!     ocean_to_ice_bundle, ocean_to_ice_uv
!   !__yac_atm (standalone, Stage 4 config #2): recv raw atm STATE — sst_feom_to_ice,
!     atm_state_to_ice (9-comp), ocean_to_ice_bundle, ocean_to_ice_uv. FESIM runs
!     its own bulk (gen_bulk_formulae) on the received state.
!
! As of Stage 3 there is no atm2ice YAC channel; the ocean is the only source.
!
! This means FESIM is now compatible with any atm-coupler that talks to
! FESOM (YAC for ICON, __ifsinterface for IFS-FESOM, ...) — FESIM sees
! only the ocean.
!
! Time-stepping: ocn_cpl_define takes the FESIM sea-ice dt in seconds. The
! interface makes NO assumption that FESIM's dt matches FESOM's dt or the
! YAC coupling_period (set in coupling.yaml).
!
! Since FESIM has only one coupling interface, this module owns the YAC
! component init/finalize directly (no shared yac_component_runtime is
! needed here, unlike the FESOM side).
module ocean_coupling_interface
#if defined(__yac)

  use yac
  use o_PARAM, only: WP

  implicit none
  private

  character(len=*), parameter, public :: OCN_COMP_NAME = "fesim"
  character(len=*), parameter, public :: OCN_GRID_NAME = "fesim_grid"

  ! --- public field-index parameters ------------------------------------
  ! Send slots (1-based). Sent in ALL configs (the ice->ocean momentum stress is
  ! ice-model-derived regardless of atm source).
  integer, parameter, public :: OCN_SEND_SEA_ICE_BUNDLE       = 1
  integer, parameter, public :: OCN_SEND_ICE_STRESS           = 2
  ! FESIM's net heat + freshwater flux to the ocean, sent in ALL configs. FESIM
  ! runs the same ice thermodynamics the monolithic FESOM2 ran (fed the atm fluxes
  ! from ICON/IFS/forcing), so ice%flx_h/flx_fw reproduce the standard-FESOM ocean
  ! surface flux; the ocean applies them via the unchanged oce_fluxes.
  integer, parameter, public :: OCN_SEND_ICE_FLUX             = 3
  integer, parameter, public :: OCN_NSEND                     = 3

  ! Recv slots (1-based). What the ocean forwards from the atmosphere depends on
  ! __yac_atm (must match the FESOM-side ice_coupling_interface send layout).
  integer, parameter, public :: OCN_RECV_SST_FEOM             = 1
#if defined(__yac_atm)
  ! ICON: receive pre-computed atm FLUXES (Stage 3 layout).
  integer, parameter, public :: OCN_RECV_TAUX                 = 2
  integer, parameter, public :: OCN_RECV_TAUY                 = 3
  integer, parameter, public :: OCN_RECV_FRESH_WATER          = 4
  integer, parameter, public :: OCN_RECV_HEAT_FLUX            = 5
  integer, parameter, public :: OCN_RECV_ATM_SEA_ICE_BUNDLE   = 6
  integer, parameter, public :: OCN_RECV_OCEAN_TO_ICE_BUNDLE  = 7
  integer, parameter, public :: OCN_RECV_OCEAN_TO_ICE_UV      = 8
  integer, parameter, public :: OCN_NRECV                     = 8
#elif defined(__ifs_fwd)
  ! IFS (Stage 4 config #1): receive the IFS atm fluxes the ocean forwards as a
  ! single 10-component bundle, already in FESIM internal units + rotated (NO
  ! conversion on receive). Must match the FESOM ice_coupling_interface send.
  ! Components: 1 stress_atmice_x 2 stress_atmice_y 3 oce_heat_flux
  !   4 ice_heat_flux 5 shortwave 6 prec_rain 7 prec_snow 8 evap_no_ifrac
  !   9 sublimation 10 enthalpyoffuse   (stress_atmoce is NOT forwarded)
  integer, parameter, public :: OCN_RECV_ATM_ICE_FLUX         = 2
  integer, parameter, public :: OCN_RECV_OCEAN_TO_ICE_BUNDLE  = 3
  integer, parameter, public :: OCN_RECV_OCEAN_TO_ICE_UV      = 4
  integer, parameter, public :: OCN_NRECV                     = 4
#else
  ! Standalone (Stage 4 config #2): receive the raw atm STATE bundle and run our
  ! own bulk. Components (FESIM converts units on receive — see gen_forcing_couple):
  !   1 u_wind [m/s]   2 v_wind [m/s]   3 t_air [K]    4 shum [kg/kg]
  !   5 shortwave [W/m2] 6 longwave [W/m2] 7 prec_rain [kg/m2/s] 8 prec_snow [kg/m2/s] 9 mslp [Pa]
  integer, parameter, public :: OCN_RECV_ATM_STATE            = 2
  integer, parameter, public :: OCN_RECV_OCEAN_TO_ICE_BUNDLE  = 3
  integer, parameter, public :: OCN_RECV_OCEAN_TO_ICE_UV      = 4
  integer, parameter, public :: OCN_NRECV                     = 4
#endif

  ! Collection sizes per field (indexed by slot constants above).
  ! IFS config #1 also relays ice surface temperature + albedo up to the ocean
  ! (onward to IFS), so the sea-ice bundle carries 5 components there.
! sea_ice_bundle (3 or 5 for IFS) + ice_to_ocean_stress (2: stress_iceoce_x/y)
! + ice_to_ocean_flux (2: net_heat_flux, fresh_wa_flux) in standalone only.
#if defined(__ifs_fwd)
  integer, parameter, public :: ocn_send_collection_size(OCN_NSEND) = [5, 2, 2]
#else
  integer, parameter, public :: ocn_send_collection_size(OCN_NSEND) = [3, 2, 2]
#endif
#if defined(__yac_atm)
  integer, parameter, public :: ocn_recv_collection_size(OCN_NRECV) = [1, 2, 2, 3, 4, 2, 2, 2]
#elif defined(__ifs_fwd)
  integer, parameter, public :: ocn_recv_collection_size(OCN_NRECV) = [1, 10, 2, 2]
#else
  integer, parameter, public :: ocn_recv_collection_size(OCN_NRECV) = [1, 9, 2, 2]
#endif

  ! YAC field names per slot (32-char strings, padded). Public so downstream
  ! diagnostic/flux-correction code can print field names.
  character(len=32), parameter, public :: ocn_send_names(OCN_NSEND) = [character(len=32) :: &
       'sea_ice_bundle', &
       'ice_to_ocean_stress', &
       'ice_to_ocean_flux' ]
#if defined(__yac_atm)
  character(len=32), parameter, public :: ocn_recv_names(OCN_NRECV) = [character(len=32) :: &
       'sst_feom_to_ice', &
       'taux_to_ice', &
       'tauy_to_ice', &
       'surface_fresh_water_flux_to_ice', &
       'total_heat_flux_to_ice', &
       'atmosphere_sea_ice_bundle_to_ice', &
       'ocean_to_ice_bundle', &
       'ocean_to_ice_uv' ]
#elif defined(__ifs_fwd)
  character(len=32), parameter, public :: ocn_recv_names(OCN_NRECV) = [character(len=32) :: &
       'sst_feom_to_ice', &
       'atm_ice_flux_to_ice', &
       'ocean_to_ice_bundle', &
       'ocean_to_ice_uv' ]
#else
  character(len=32), parameter, public :: ocn_recv_names(OCN_NRECV) = [character(len=32) :: &
       'sst_feom_to_ice', &
       'atm_state_to_ice', &
       'ocean_to_ice_bundle', &
       'ocean_to_ice_uv' ]
#endif

  ! Legacy state preserved for downstream flux-correction routines that
  ! still expect these symbols (carried forward from cpl_yac_driver_fesim
  ! verbatim). Under the YAC backend these are declared but never set —
  ! same as in the legacy driver.
  real(kind=WP), allocatable, public :: a2o_fcorr_stat(:,:)
  integer,                   public :: source_root  = 0
  integer,                   public :: target_root  = 0
  logical,                   public :: commRank     = .false.

  ! --- private module state ---------------------------------------------
  integer, save :: ocn_comp_id          = -1
  integer, save :: ocn_local_comm       = -1
  integer, save :: ocn_grid_id          = -1
  integer, save :: ocn_points_id        = -1
  integer, save :: ocn_send_field_id(OCN_NSEND) = -1
  integer, save :: ocn_recv_field_id(OCN_NRECV) = -1
  logical, save :: ocn_inited           = .false.


  ! --- coupling-cost instrumentation (2026-09-11, §10.2) ------------------
  ! yac_fput/yac_fget are the only points at which the pair actually exchanges,
  ! so timing them here attributes the coupler cost without touching callers.
  ! Whatever else the caller's exchange routine does -- halo exchanges, unit
  ! conversion, copies -- then shows up as the remainder of the caller's own
  ! timer, which is what makes the split interpretable. Before this, the ocean
  ! side timed the exchange nowhere at all and FESIM lumped it in with its wait.
  real(kind=WP), public, save :: cpl_time_put = 0.0_WP   ! s, this rank, cumulative
  real(kind=WP), public, save :: cpl_time_get = 0.0_WP
  integer,       public, save :: cpl_n_put    = 0        ! calls (not exchanges:
  integer,       public, save :: cpl_n_get    = 0        !  YAC no-ops off-period)
  ! One-shot consistency check: is coupling_period actually equal to
  ! dt*cpl_stride? The Fortran side registers its fields with dt*cpl_stride, but
  ! coupling_period lives in coupling.yaml and NOTHING enforces agreement. If the
  ! yaml period is larger, YAC couples only every Nth call; the rest return
  ! no-action and the component silently reuses STALE fields. That zero-order
  ! hold blew up both high-resolution meshes for weeks and was misdiagnosed as a
  ! physics bug -- see the 2026-09-11 A/B test (PT5M ran 864 steps clean, PT30M
  ! against dt=300 died at step 101 with a 100x surface heat flux).
  integer, save :: cpl_mype    = -1
  logical, save :: cpl_checked = .false.

  integer,       public, save :: cpl_n_put_act = 0       ! calls that actually coupled
  integer,       public, save :: cpl_n_get_act = 0

  ! --- public API -------------------------------------------------------
  public :: ocn_cpl_init, ocn_cpl_define, ocn_cpl_send, ocn_cpl_recv, ocn_cpl_finalize
  public :: cpl_timers_report

contains

  subroutine ocn_cpl_init(localCommunicator)
    integer, intent(out) :: localCommunicator
    if (.not. ocn_inited) then
#ifdef VERBOSE
       print *, '================================================='
       print *, 'ocn_cpl_init : coupler initialization for YAC'
       print *, '*************************************************'
#endif
       call yac_finit()
       call yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
       call yac_fread_config_yaml("coupling.yaml")
       call yac_fdef_comp(OCN_COMP_NAME, ocn_comp_id)
       call yac_fget_comp_comm(ocn_comp_id, ocn_local_comm)
       ocn_inited = .true.
    end if
    localCommunicator = ocn_local_comm
  end subroutine ocn_cpl_init

  subroutine ocn_cpl_define(partit, mesh, dt_seconds)
    use MOD_MESH,       only: t_mesh
    use MOD_PARTIT,     only: t_partit
    use yac_grid_utils, only: cpl_yac_define_unstr_generic
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    integer,        intent(in)            :: dt_seconds

    character(len=4) :: dt_str
    integer          :: ierr, i

    call cpl_yac_define_unstr_generic(partit, mesh, OCN_GRID_NAME, ocn_grid_id, ocn_points_id)

    cpl_mype = partit%mype

    write(dt_str, '(I4.4)') dt_seconds

    do i = 1, OCN_NSEND
       call yac_fdef_field(ocn_send_names(i), ocn_comp_id, &
            [ocn_points_id], 1, ocn_send_collection_size(i), &
            dt_str, YAC_TIME_UNIT_SECOND, ocn_send_field_id(i))
    end do

    do i = 1, OCN_NRECV
       call yac_fdef_field(ocn_recv_names(i), ocn_comp_id, &
            [ocn_points_id], 1, ocn_recv_collection_size(i), &
            dt_str, YAC_TIME_UNIT_SECOND, ocn_recv_field_id(i))
    end do

    call yac_fenddef(ierr)
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
  ! per-task runtime block. mean/min/max over ranks, as in that block.
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
    write(*,*) '  runaways at high resolution -- it killed fArc/DARS at ~step 100.'
    write(*,*) '  FIX: set coupling_period in coupling.yaml equal to dt*cpl_stride,'
    write(*,*) '  or, for genuine asynchronous coupling, add "time_reduction:'
    write(*,*) '  average" to the couple so the source averages over the period'
    write(*,*) '  instead of sending an instantaneous snapshot.'
    write(*,*) '  NB each component reads the coupling.yaml in ITS OWN directory.'
    write(*,*) '**********************************************************************'
  end subroutine cpl_check_period


#endif
end module ocean_coupling_interface

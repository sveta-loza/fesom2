! Sea-ice coupling interface (FESOM ocean side).
!
! This module is the seam through which FESOM exchanges fields with FESIM.
! The backend is always YAC. The sibling module on the FESIM side is
! ocean_coupling_interface.
!
! Field inventory (ice side from the ocean's perspective). The native ocean
! state is the same in all configs; what the ocean forwards from the atmosphere
! depends on how the atmosphere reaches the ocean (compile-time __yac_atm):
!
!   __yac_atm (ICON, Stage 3): the ocean receives pre-computed atm FLUXES from
!     ICON via YAC and re-forwards them. sends to ice:
!       - native: sst_feom_to_ice, ocean_to_ice_bundle, ocean_to_ice_uv
!       - forwarded fluxes: taux_to_ice, tauy_to_ice,
!         surface_fresh_water_flux_to_ice, total_heat_flux_to_ice,
!         atmosphere_sea_ice_bundle_to_ice
!
!   !__yac_atm (standalone FESOM-FESIM, Stage 4 config #2): the ocean reads atm
!     forcing from files and forwards the raw atm STATE; FESIM runs its own bulk.
!     sends to ice:
!       - native: sst_feom_to_ice, ocean_to_ice_bundle, ocean_to_ice_uv
!       - forwarded state: atm_state_to_ice (9-component bundle)
!
!   recvs from ice (both): sea_ice_bundle
!
! river_runoff is intentionally NOT forwarded — it is ocean-only, no FESIM
! consumer. Unit convention: forwarded fields travel in raw/source units; FESIM
! converts on receive (see ocean_coupling_interface on the FESIM side).
!
! NB: `sst_feom_to_ice` is the *duplicate* of the atm-bound `sst_feom`.
! Per Stage 2a design decision, the ocean code calls both
! atm_cpl_send(ATM_SEND_SST_FEOM, ...) and
! ice_cpl_send(ICE_SEND_SST_FEOM, ...) with the same data each timestep,
! so each interface owns its end-to-end YAC field. The yaml topology
! needs the oce2ice block to read `src: sst_feom_to_ice` instead of
! `sst_feom` — coupling.yaml update lands in 2c alongside caller wiring.
!
! Time-stepping: ice_cpl_define takes the FESOM ocean dt in seconds. The
! ice-side fields here are *emitted by FESOM*, so they carry FESOM's dt as
! their yac_fdef_field time step. FESIM's dt is registered on the sibling
! ocean_coupling_interface in the FESIM tree.
!
! YAC component init and grid definition are shared with
! atm_coupling_interface via yac_component_runtime.
!
! NB: callers must invoke yac_runtime_enddef() once after both
! atm_cpl_define and ice_cpl_define have run.
module ice_coupling_interface
#if defined(__yac)

  use yac
  use o_PARAM, only: WP

  implicit none
  private

  ! --- public field-index parameters ------------------------------------
  ! Send slots (1-based). Native ocean state first (same in all configs), then
  ! the forwarded atm fields whose layout depends on __yac_atm.
  integer, parameter, public :: ICE_SEND_SST_FEOM             = 1
  integer, parameter, public :: ICE_SEND_OCEAN_TO_ICE_BUNDLE  = 2
  integer, parameter, public :: ICE_SEND_OCEAN_TO_ICE_UV      = 3
#if defined(__yac_atm)
  ! ICON: forward pre-computed atm fluxes (Stage 3).
  integer, parameter, public :: ICE_SEND_TAUX                 = 4
  integer, parameter, public :: ICE_SEND_TAUY                 = 5
  integer, parameter, public :: ICE_SEND_FRESH_WATER          = 6
  integer, parameter, public :: ICE_SEND_HEAT_FLUX            = 7
  integer, parameter, public :: ICE_SEND_ATM_SEA_ICE_BUNDLE   = 8
  integer, parameter, public :: ICE_NSEND                     = 8
#elif defined(__ifs_fwd)
  ! IFS (Stage 4 config #1): the atmosphere reaches the ocean via the IFS
  ! interface, which deposits *already FESOM-unit* fluxes (precip/evap in m/s,
  ! stresses in Pa and already rotated to the FESOM grid, ice-fraction weighted).
  ! We forward them 1:1 as a single 10-component bundle; FESIM does NOT re-convert
  ! or re-rotate. stress_atmoce is intentionally NOT forwarded — it drives ocean
  ! momentum, which lives in this ocean component, not in FESIM. Components (wire
  ! units == FESOM internal units):
  !   1 stress_atmice_x [Pa]   2 stress_atmice_y [Pa]
  !   3 oce_heat_flux [W/m2]   4 ice_heat_flux [W/m2]   5 shortwave [W/m2]
  !   6 prec_rain [m/s]        7 prec_snow [m/s]         8 evap_no_ifrac [m/s]
  !   9 sublimation [m/s]     10 enthalpyoffuse [W/m2]
  integer, parameter, public :: ICE_SEND_ATM_ICE_FLUX         = 4
  integer, parameter, public :: ICE_NSEND                     = 4
#else
  ! Standalone: forward raw atm state as one 9-component bundle (FESIM bulks it).
  ! Components (wire units; FESIM converts on receive):
  !   1 u_wind [m/s]   2 v_wind [m/s]   3 t_air [K]    4 shum [kg/kg]
  !   5 shortwave [W/m2] down           6 longwave [W/m2] down
  !   7 prec_rain [kg/m2/s]             8 prec_snow [kg/m2/s]   9 mslp [Pa]
  integer, parameter, public :: ICE_SEND_ATM_STATE            = 4
  integer, parameter, public :: ICE_NSEND                     = 4
#endif

  ! Recv slots (1-based). ice_to_ocean_stress (the ice->ocean momentum drag) is
  ! received in ALL configs — it is ice-model-derived regardless of atm source.
  integer, parameter, public :: ICE_RECV_SEA_ICE_BUNDLE       = 1
  integer, parameter, public :: ICE_RECV_ICE_STRESS           = 2
  ! FESIM's net heat + freshwater flux, received in ALL configs (FESIM's ice thermo
  ! produces the standard-FESOM ocean surface flux regardless of atm source).
  integer, parameter, public :: ICE_RECV_ICE_FLUX             = 3
  integer, parameter, public :: ICE_NRECV                     = 3

  ! Collection sizes per field.
#if defined(__yac_atm)
  ! Forwarded flux sizes mirror atm_recv_collection_size (taux=2, tauy=2,
  ! fresh_water=3, heat_flux=4, atm_sea_ice_bundle=2).
  integer, parameter, public :: ice_send_collection_size(ICE_NSEND) = [1, 2, 2, 2, 2, 3, 4, 2]
#elif defined(__ifs_fwd)
  integer, parameter, public :: ice_send_collection_size(ICE_NSEND) = [1, 2, 2, 10]
#else
  integer, parameter, public :: ice_send_collection_size(ICE_NSEND) = [1, 2, 2, 9]
#endif
  ! IFS config #1 also relays the ice surface temperature + albedo up to IFS, so
  ! the sea-ice bundle carries 5 components there (m_ice, m_snow, a_ice, ice_temp,
  ! ice_alb) instead of 3.
  ! sea_ice_bundle (3, or 5 for IFS) + ice_to_ocean_stress (2: stress_iceoce_x/y)
  ! + ice_to_ocean_flux (2: net_heat_flux, fresh_wa_flux) in standalone only.
#if defined(__ifs_fwd)
  integer, parameter, public :: ice_recv_collection_size(ICE_NRECV) = [5, 2, 2]
#else
  integer, parameter, public :: ice_recv_collection_size(ICE_NRECV) = [3, 2, 2]
#endif

  ! YAC field names per slot. The "_to_ice" suffix distinguishes forwarded
  ! fields from the atm-bound originals (each YAC component+grid needs unique
  ! field names). The FESIM-side ocean_coupling_interface receives these names.
#if defined(__yac_atm)
  character(len=32), parameter, public :: ice_send_names(ICE_NSEND) = [character(len=32) :: &
       'sst_feom_to_ice', &
       'ocean_to_ice_bundle', &
       'ocean_to_ice_uv', &
       'taux_to_ice', &
       'tauy_to_ice', &
       'surface_fresh_water_flux_to_ice', &
       'total_heat_flux_to_ice', &
       'atmosphere_sea_ice_bundle_to_ice' ]
#elif defined(__ifs_fwd)
  character(len=32), parameter, public :: ice_send_names(ICE_NSEND) = [character(len=32) :: &
       'sst_feom_to_ice', &
       'ocean_to_ice_bundle', &
       'ocean_to_ice_uv', &
       'atm_ice_flux_to_ice' ]
#else
  character(len=32), parameter, public :: ice_send_names(ICE_NSEND) = [character(len=32) :: &
       'sst_feom_to_ice', &
       'ocean_to_ice_bundle', &
       'ocean_to_ice_uv', &
       'atm_state_to_ice' ]
#endif
  character(len=32), parameter, public :: ice_recv_names(ICE_NRECV) = [character(len=32) :: &
       'sea_ice_bundle', &
       'ice_to_ocean_stress', &
       'ice_to_ocean_flux' ]

  ! Grid this interface registers its fields on. Today identical to
  ! ATM_GRID_NAME; could diverge in the future if atm and ice need
  ! different meshes.
  character(len=*), parameter, public :: ICE_GRID_NAME = "fesom_grid"

  ! --- private module state ---------------------------------------------
  integer, save :: ice_send_field_id(ICE_NSEND) = -1
  integer, save :: ice_recv_field_id(ICE_NRECV) = -1
  integer, save :: ice_points_id_local          = -1


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
  integer,       public, save :: cpl_n_put_act = 0       ! calls that actually coupled
  integer,       public, save :: cpl_n_get_act = 0
  ! Whole-routine bracket for exchange_oce_ice_yac. The ocean prints no timer
  ! covering the exchange (rtime_ice comes from ice_timestep, which the decoupled
  ! ocean never calls), so without this the halo exchanges and copies around the
  ! YAC calls can only be got at by subtracting from the loop wall.
  real(kind=WP), public, save :: cpl_time_call = 0.0_WP

  ! --- public API -------------------------------------------------------
  public :: ice_cpl_init, ice_cpl_define, ice_cpl_send, ice_cpl_recv, ice_cpl_finalize
  public :: cpl_timers_report, cpl_call_tic, cpl_call_toc

contains

  subroutine ice_cpl_init(localCommunicator)
    use yac_component_runtime, only: yac_runtime_init
    integer, intent(out) :: localCommunicator
    call yac_runtime_init(localCommunicator)
  end subroutine ice_cpl_init

  subroutine ice_cpl_define(partit, mesh, dt_seconds)
    use MOD_MESH,              only: t_mesh
    use MOD_PARTIT,            only: t_partit
    use yac_component_runtime, only: yac_runtime_ensure_grid, yac_runtime_comp_id
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    integer,        intent(in)            :: dt_seconds

    character(len=4) :: dt_str
    integer          :: grid_id, i

    call yac_runtime_ensure_grid(ICE_GRID_NAME, partit, mesh, grid_id, ice_points_id_local)

    write(dt_str, '(I4.4)') dt_seconds

    do i = 1, ICE_NSEND
       call yac_fdef_field(ice_send_names(i), yac_runtime_comp_id(), &
            [ice_points_id_local], 1, ice_send_collection_size(i), &
            dt_str, YAC_TIME_UNIT_SECOND, ice_send_field_id(i))
    end do

    do i = 1, ICE_NRECV
       call yac_fdef_field(ice_recv_names(i), yac_runtime_comp_id(), &
            [ice_points_id_local], 1, ice_recv_collection_size(i), &
            dt_str, YAC_TIME_UNIT_SECOND, ice_recv_field_id(i))
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


  ! Collective report of the coupler cost, printed next to the model's own
  ! per-task runtime block. mean/min/max over ranks, as in that block.
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

  subroutine cpl_timers_report(comm, mype, npes, label)
    use mpi
    integer,          intent(in) :: comm, mype, npes
    character(len=*), intent(in) :: label
    ! 1 fput, 2 fget, 3 whole exchange routine, 4 the routine minus YAC (halo
    ! exchanges + copies + unit work). (4) is formed per rank BEFORE reducing --
    ! reducing (3) and (1)+(2) separately and subtracting would not give the
    ! min/max of the difference.
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

#endif
end module ice_coupling_interface

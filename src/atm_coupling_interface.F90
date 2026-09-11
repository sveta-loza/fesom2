! Atmosphere coupling interface (FESOM ocean side).
!
! This module is the seam through which FESOM exchanges fields with the
! atmosphere. Today the only backend is YAC (for ICON). Stage 4 will add an
! __ifsinterface backend so IFS-FESOM can coexist with the same interface.
!
! The interface is intentionally thin: it presents the same shape as the
! legacy cpl_yac_driver (init/define/send/recv/finalize plus field metadata)
! but its public symbols are now organised by *which counterpart* the data
! is for. The sibling ice_coupling_interface owns ice-bound traffic.
!
! Field inventory (atm side, YAC backend):
!   sends to atm: sst_feom, ocean_sea_ice_bundle
!   recvs from atm: taux, tauy, surface_fresh_water_flux, total_heat_flux,
!                   atmosphere_sea_ice_bundle, river_runoff
!
! Time-stepping: atm_cpl_define takes the FESOM model time step in seconds
! as an explicit argument. This interface makes NO assumption that the
! atm, ocean and sea-ice components share a common dt, nor that the YAC
! coupling_period (set in coupling.yaml) matches any component's dt.
!
! YAC component init (yac_finit etc.) and the FESOM grid definition live
! in yac_component_runtime. atm_cpl_init / atm_cpl_define / atm_cpl_finalize
! delegate to the runtime, which is idempotent.
!
! NB: callers must invoke yac_component_runtime%yac_runtime_enddef() once
! after both atm_cpl_define and ice_cpl_define have run.
module atm_coupling_interface
#if defined(__yac)

  use yac
  use o_PARAM, only: WP

  implicit none
  private

  ! --- public field-index parameters ------------------------------------
  ! Send slots (1-based, used as the `ind` argument to atm_cpl_send).
  integer, parameter, public :: ATM_SEND_SST_FEOM             = 1
  integer, parameter, public :: ATM_SEND_OCEAN_SEA_ICE_BUNDLE = 2
  integer, parameter, public :: ATM_NSEND                     = 2

  ! Recv slots (1-based, used as the `ind` argument to atm_cpl_recv).
  integer, parameter, public :: ATM_RECV_TAUX                 = 1
  integer, parameter, public :: ATM_RECV_TAUY                 = 2
  integer, parameter, public :: ATM_RECV_FRESH_WATER          = 3
  integer, parameter, public :: ATM_RECV_HEAT_FLUX            = 4
  integer, parameter, public :: ATM_RECV_ATM_SEA_ICE_BUNDLE   = 5
  integer, parameter, public :: ATM_RECV_RIVER_RUNOFF         = 6
  integer, parameter, public :: ATM_NRECV                     = 6

  ! Collection sizes per field (indexed by slot constants above). Public so
  ! callers can size their exchange buffers without round-tripping through
  ! the module. The values match the yac_fdef_field collection_size arg
  ! used in atm_cpl_define.
  integer, parameter, public :: atm_send_collection_size(ATM_NSEND) = [1, 3]
  integer, parameter, public :: atm_recv_collection_size(ATM_NRECV) = [2, 2, 3, 4, 2, 1]

  ! YAC field names per slot. Public so downstream flux-correction code
  ! can print field names.
  character(len=32), parameter, public :: atm_send_names(ATM_NSEND) = [character(len=32) :: &
       'sst_feom', &
       'ocean_sea_ice_bundle' ]
  character(len=32), parameter, public :: atm_recv_names(ATM_NRECV) = [character(len=32) :: &
       'taux', &
       'tauy', &
       'surface_fresh_water_flux', &
       'total_heat_flux', &
       'atmosphere_sea_ice_bundle', &
       'river_runoff' ]

  ! Grid this interface registers its fields on. Today shared with
  ! ice_coupling_interface; could diverge in the future without an API
  ! change. See project memory "Future mesh split between atm- and
  ! ice-facing exchanges".
  character(len=*), parameter, public :: ATM_GRID_NAME = "fesom_grid"

  ! --- public module state ----------------------------------------------
  ! Flux correction statistics for output. Atm-side only (river runoff /
  ! freshwater bookkeeping). Allocated by caller (e.g. gen_forcing_couple).
  real(kind=WP), allocatable, public :: a2o_fcorr_stat(:,:)

  ! Legacy state preserved for downstream flux-correction routines that
  ! still expect these symbols (carried forward from cpl_yac_driver
  ! verbatim). Under the YAC backend these are declared but never set —
  ! same as in the legacy driver.
  integer, public :: source_root = 0
  integer, public :: target_root = 0
  logical, public :: commRank    = .false.

  ! --- private module state ---------------------------------------------
  ! YAC field handles, populated by atm_cpl_define and consumed by
  ! atm_cpl_send / atm_cpl_recv.
  integer, save :: atm_send_field_id(ATM_NSEND) = -1
  integer, save :: atm_recv_field_id(ATM_NRECV) = -1
  integer, save :: atm_points_id_local          = -1

  ! --- public API -------------------------------------------------------
  public :: atm_cpl_init, atm_cpl_define, atm_cpl_send, atm_cpl_recv, atm_cpl_finalize

contains

  subroutine atm_cpl_init(localCommunicator)
    use yac_component_runtime, only: yac_runtime_init
    integer, intent(out) :: localCommunicator
    call yac_runtime_init(localCommunicator)
  end subroutine atm_cpl_init

  subroutine atm_cpl_define(partit, mesh, dt_seconds)
    use MOD_MESH,              only: t_mesh
    use MOD_PARTIT,            only: t_partit
    use yac_component_runtime, only: yac_runtime_ensure_grid, yac_runtime_comp_id
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    integer,        intent(in)            :: dt_seconds

    character(len=4) :: dt_str
    integer          :: grid_id, i

    call yac_runtime_ensure_grid(ATM_GRID_NAME, partit, mesh, grid_id, atm_points_id_local)

    write(dt_str, '(I4.4)') dt_seconds

    do i = 1, ATM_NSEND
       call yac_fdef_field(atm_send_names(i), yac_runtime_comp_id(), &
            [atm_points_id_local], 1, atm_send_collection_size(i), &
            dt_str, YAC_TIME_UNIT_SECOND, atm_send_field_id(i))
    end do

    do i = 1, ATM_NRECV
       call yac_fdef_field(atm_recv_names(i), yac_runtime_comp_id(), &
            [atm_points_id_local], 1, atm_recv_collection_size(i), &
            dt_str, YAC_TIME_UNIT_SECOND, atm_recv_field_id(i))
    end do
  end subroutine atm_cpl_define

  subroutine atm_cpl_send(ind, data_array, action)
    integer,       intent(in)  :: ind
    real(kind=WP), intent(in)  :: data_array(:,:)
    logical,       intent(out) :: action
    integer :: info, ierr
    call yac_fput(atm_send_field_id(ind), size(data_array, 1), size(data_array, 2), &
         data_array, info, ierr)
    action = info == YAC_ACTION_COUPLING
  end subroutine atm_cpl_send

  subroutine atm_cpl_recv(ind, data_array, action)
    integer,       intent(in)    :: ind
    real(kind=WP), intent(inout) :: data_array(:,:)
    logical,       intent(out)   :: action
    integer :: info, ierr
    call yac_fget(atm_recv_field_id(ind), size(data_array, 1), size(data_array, 2), &
         data_array, info, ierr)
    action = info == YAC_ACTION_COUPLING
  end subroutine atm_cpl_recv

  subroutine atm_cpl_finalize()
    use yac_component_runtime, only: yac_runtime_finalize
    call yac_runtime_finalize()
  end subroutine atm_cpl_finalize

#endif
end module atm_coupling_interface

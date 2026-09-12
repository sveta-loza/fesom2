! Atmosphere coupling interface (FESOM ocean side, YAC).
!
! The seam through which FESOM exchanges fields with a YAC-coupled
! atmosphere (ICON-A, is_coupled_to_icon_a). The sibling
! ice_coupling_interface owns the sea-ice-bound traffic (FESIM); both
! register on the same YAC component through yac_component_runtime.
!
! Field inventory:
!   sends to atm:   sst_feom, ocean_sea_ice_bundle
!   recvs from atm: taux, tauy, surface_fresh_water_flux, total_heat_flux,
!                   atmosphere_sea_ice_bundle, river_runoff
!
! atm_cpl_define takes the field time step in seconds as an argument and
! assumes nothing about the other components' time steps or about the
! coupling_period in coupling.yaml.
!
! NB: fesom_module calls yac_runtime_enddef() once after every interface
! has registered its fields.
module atm_coupling_interface
#if defined (__cpl_yac)

  use yac
  use o_PARAM, only: WP

  implicit none
  private

  ! --- public field-index parameters ------------------------------------
  ! Send slots (1-based, the `ind` argument of atm_cpl_send).
  integer, parameter, public :: ATM_SEND_SST_FEOM             = 1
  integer, parameter, public :: ATM_SEND_OCEAN_SEA_ICE_BUNDLE = 2
  integer, parameter, public :: ATM_NSEND                     = 2

  ! Recv slots (1-based, the `ind` argument of atm_cpl_recv).
  integer, parameter, public :: ATM_RECV_TAUX                 = 1
  integer, parameter, public :: ATM_RECV_TAUY                 = 2
  integer, parameter, public :: ATM_RECV_FRESH_WATER          = 3
  integer, parameter, public :: ATM_RECV_HEAT_FLUX            = 4
  integer, parameter, public :: ATM_RECV_ATM_SEA_ICE_BUNDLE   = 5
  integer, parameter, public :: ATM_RECV_RIVER_RUNOFF         = 6
  integer, parameter, public :: ATM_NRECV                     = 6

  ! Collection sizes per field, indexed by the slot constants. Public so
  ! callers can size their exchange buffers. They are the collection_size
  ! arguments of yac_fdef_field in atm_cpl_define.
  integer, parameter, public :: atm_send_collection_size(ATM_NSEND) = [1, 3]
  integer, parameter, public :: atm_recv_collection_size(ATM_NRECV) = [2, 2, 3, 4, 2, 1]

  ! YAC field names per slot. Must match coupling.yaml.
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

  ! --- public module state ----------------------------------------------
  ! Flux-correction statistics for output (allocated by the caller).
  real(kind=WP), allocatable, public :: a2o_fcorr_stat(:,:)

  ! Kept for the flux-correction routines that still name these symbols
  ! (force_flux_consv, net_rec_from_atm). Never set under YAC, as before.
  integer, public :: source_root = 0
  integer, public :: target_root = 0
  logical, public :: commRank    = .false.

  ! --- private module state ---------------------------------------------
  integer, save :: atm_send_field_id(ATM_NSEND) = -1
  integer, save :: atm_recv_field_id(ATM_NRECV) = -1
  integer, save :: atm_points_id_local          = -1

  public :: atm_cpl_init, atm_cpl_define, atm_cpl_send, atm_cpl_recv, atm_cpl_finalize

contains

  subroutine atm_cpl_init(localCommunicator)
    use yac_component_runtime, only: yac_runtime_init
    integer, intent(out) :: localCommunicator
    call yac_runtime_init(localCommunicator)
  end subroutine atm_cpl_init

  ! dt_seconds: the time step the fields are registered with (in seconds;
  ! fractional values are kept, YAC gets milliseconds).
  subroutine atm_cpl_define(partit, mesh, dt_seconds, nsteps)
    use MOD_MESH,              only: t_mesh
    use MOD_PARTIT,            only: t_partit
    use cpl_config,            only: cpl_grid_name
    use yac_component_runtime, only: yac_runtime_ensure_grid, yac_runtime_comp_id
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP),  intent(in)            :: dt_seconds
    integer,        intent(in)            :: nsteps       ! steps of this run (YAC end datetime)

    character(len=8) :: dt_str
    integer          :: grid_id, i

    call yac_runtime_ensure_grid(trim(cpl_grid_name), partit, mesh, grid_id, atm_points_id_local, nsteps)

    write(dt_str, '(I8.8)') INT(dt_seconds*1000)

    do i = 1, ATM_NSEND
       call yac_fdef_field(atm_send_names(i), yac_runtime_comp_id(), &
            [atm_points_id_local], 1, atm_send_collection_size(i), &
            dt_str, YAC_TIME_UNIT_MILLISECOND, atm_send_field_id(i))
    end do

    do i = 1, ATM_NRECV
       call yac_fdef_field(atm_recv_names(i), yac_runtime_comp_id(), &
            [atm_points_id_local], 1, atm_recv_collection_size(i), &
            dt_str, YAC_TIME_UNIT_MILLISECOND, atm_recv_field_id(i))
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

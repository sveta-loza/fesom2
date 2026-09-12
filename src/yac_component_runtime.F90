! Shared YAC component runtime for the FESOM ocean.
!
! FESOM registers fields on ONE YAC component through two interfaces,
! atm_coupling_interface (atmosphere-bound traffic) and
! ice_coupling_interface (sea-ice-bound traffic, FESIM). The YAC calls that
! must happen exactly once per component -- yac_finit, yac_fdef_calendar,
! yac_fread_config_yaml, yac_fdef_comp, yac_fget_comp_comm, yac_fenddef,
! yac_ffinalize -- and the grid definition the interfaces share live here.
!
! Every entry point is idempotent: the first caller does the work, later
! callers see the cached state, so the interfaces need no coordination.
!
! Caller pattern (fesom_module):
!     call atm_cpl_init(localCommunicator)       ! -> yac_runtime_init
!     call atm_cpl_define(partit, mesh, dt)      ! atm fields (if an atmosphere is coupled)
!     call ice_cpl_define(partit, mesh, dt)      ! ice fields (if FESIM is coupled)
!     call yac_runtime_enddef()                  ! closes the definition phase
!     ... time loop ...
!     call atm_cpl_finalize()                    ! -> yac_runtime_finalize
!
! Component, grid and configuration-file names come from &coupling_yac in
! namelist.cpl (cpl_config). yac_runtime_ensure_grid takes the grid name as
! an argument and caches one grid; when the atmosphere-facing and ice-facing
! exchanges need different grids, callers pass different names and the
! cache becomes an array.
!
! This module is dt-agnostic: each interface sets its own field time step in
! yac_fdef_field.
module yac_component_runtime
#if defined (__cpl_yac)

  use yac
  use yac_grid_utils, only: cpl_yac_define_unstr_generic
  use cpl_config,     only: cpl_comp_name, cpl_config_file

  implicit none
  private

  ! Cached state, exposed through the getters below.
  integer, save :: comp_id_         = -1
  integer, save :: local_comm_      = -1
  integer, save :: grid_id_         = -1
  integer, save :: points_id_       = -1
  character(len=64), save :: cached_grid_name_ = ""
  logical, save :: runtime_inited_  = .false.
  logical, save :: grid_defined_    = .false.
  logical, save :: enddef_done_     = .false.

  public :: yac_runtime_init, yac_runtime_ensure_grid
  public :: yac_runtime_enddef, yac_runtime_finalize
  public :: yac_runtime_comp_id, yac_runtime_local_comm
  public :: yac_runtime_grid_id, yac_runtime_points_id
  public :: yac_runtime_is_inited, yac_runtime_grid_is_defined

contains

  ! First call: yac_finit, calendar, config file, component definition and
  ! the component communicator. Later calls return the cached communicator.
  subroutine yac_runtime_init(localCommunicator)
    integer, intent(out) :: localCommunicator

    if (.not. runtime_inited_) then
#ifdef VERBOSE
       print *, '================================================='
       print *, 'yac_runtime_init : coupler initialization for YAC'
       print *, '*************************************************'
#endif
       call yac_finit()
       call yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
       call yac_fread_config_yaml(trim(cpl_config_file))
       call yac_fdef_comp(trim(cpl_comp_name), comp_id_)
       call yac_fget_comp_comm(comp_id_, local_comm_)
       runtime_inited_ = .true.
    end if
    localCommunicator = local_comm_
  end subroutine yac_runtime_init

  ! First call per grid name: defines the grid (and the YAC datetime).
  ! Later calls with the same name return the cached ids; a different name
  ! is a fatal error until multi-grid support is needed.
  subroutine yac_runtime_ensure_grid(grid_name, partit, mesh, grid_id, points_id, nsteps)
    use MOD_MESH,   only: t_mesh
    use MOD_PARTIT, only: t_partit
    character(len=*), intent(in)            :: grid_name
    type(t_mesh),     intent(in),    target :: mesh
    type(t_partit),   intent(inout), target :: partit
    integer,          intent(out)           :: grid_id, points_id
    integer,          intent(in)            :: nsteps

    if (.not. grid_defined_) then
       call cpl_yac_define_unstr_generic(partit, mesh, grid_name, grid_id_, points_id_, nsteps)
       cached_grid_name_ = grid_name
       grid_defined_     = .true.
    else if (trim(cached_grid_name_) /= trim(grid_name)) then
       write(0,*) "yac_runtime_ensure_grid: multi-grid not yet supported. ", &
            "Cached='", trim(cached_grid_name_), "' requested='", trim(grid_name), "'"
       stop 1
    end if
    grid_id   = grid_id_
    points_id = points_id_
  end subroutine yac_runtime_ensure_grid

  ! First call: yac_fenddef. Must run once after every interface has
  ! registered its fields.
  subroutine yac_runtime_enddef()
    if (.not. enddef_done_) then
       call yac_fenddef()
       enddef_done_ = .true.
    end if
  end subroutine yac_runtime_enddef

  ! First call: yac_ffinalize. Later calls no-op.
  subroutine yac_runtime_finalize()
    if (runtime_inited_) then
#ifdef VERBOSE
       print *, '================================================='
       print *, 'yac_runtime_finalize : coupler finalization for YAC'
       print *, '*************************************************'
#endif
       call yac_ffinalize()
       runtime_inited_ = .false.
       grid_defined_   = .false.
       enddef_done_    = .false.
    end if
  end subroutine yac_runtime_finalize

  pure integer function yac_runtime_comp_id()
    yac_runtime_comp_id = comp_id_
  end function yac_runtime_comp_id

  pure integer function yac_runtime_local_comm()
    yac_runtime_local_comm = local_comm_
  end function yac_runtime_local_comm

  pure integer function yac_runtime_grid_id()
    yac_runtime_grid_id = grid_id_
  end function yac_runtime_grid_id

  pure integer function yac_runtime_points_id()
    yac_runtime_points_id = points_id_
  end function yac_runtime_points_id

  pure logical function yac_runtime_is_inited()
    yac_runtime_is_inited = runtime_inited_
  end function yac_runtime_is_inited

  pure logical function yac_runtime_grid_is_defined()
    yac_runtime_grid_is_defined = grid_defined_
  end function yac_runtime_grid_is_defined

#endif
end module yac_component_runtime

! Shared YAC component runtime for the FESOM ocean side.
!
! The FESOM tree has TWO coupling interfaces (atm_coupling_interface +
! ice_coupling_interface) that both register fields on the same YAC
! component "fesom2". The YAC primitives that must be called exactly once
! per component (yac_finit, yac_fdef_calendar, yac_fread_config_yaml,
! yac_fdef_comp, yac_fget_comp_comm, yac_ffinalize) and the grid
! definition that today is shared between them live here.
!
! All public entry points are idempotent: the first caller does the work,
! subsequent callers see cached state. This lets each interface call
! `yac_runtime_init` / `yac_runtime_ensure_grid` / `yac_runtime_enddef` /
! `yac_runtime_finalize` independently without coordination, in any order.
!
! Caller pattern (typically from fesom_main / fesom_module):
!     call atm_cpl_init(localCommunicator)    ! delegates to yac_runtime_init
!     call atm_cpl_define(partit, mesh, dt)   ! registers atm fields
!     call ice_cpl_define(partit, mesh, dt)   ! registers ice fields
!     call yac_runtime_enddef()               ! closes YAC field-define phase
!     ... timeloop ...
!     call atm_cpl_finalize()                 ! delegates to yac_runtime_finalize
!
! Future-mesh hook (kept design-only for now — see
! yac_grid_utils.F90:cpl_yac_define_unstr_generic for the matching
! mesh-topology caveat): `yac_runtime_ensure_grid` takes the grid name
! as an argument and caches per name. Today every caller passes
! "fesom_grid". When atm-facing and ice-facing channels need different
! grids, they pass different names and this module manages the cache.
! The current implementation supports only one cached grid; multi-grid
! support is a one-array bump kept for the future stage.
!
! Time-stepping: this module is dt-agnostic. Per-field dt is set by each
! interface's _define routine when calling yac_fdef_field. See project
! memory "Component timesteps are independent".
module yac_component_runtime
#if defined(__yac)

  use yac
  use yac_grid_utils, only: cpl_yac_define_unstr_generic

  implicit none
  private

  character(len=*), parameter, public :: yac_comp_name = "fesom2"

  ! Cached state, exposed as read-only getters below.
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

  ! Idempotent: first call runs yac_finit, yac_fdef_calendar,
  ! yac_fread_config_yaml("coupling.yaml"), yac_fdef_comp(yac_comp_name)
  ! and yac_fget_comp_comm. Subsequent calls return the cached local
  ! communicator without re-entering YAC.
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
       call yac_fread_config_yaml("coupling.yaml")
       call yac_fdef_comp(yac_comp_name, comp_id_)
       call yac_fget_comp_comm(comp_id_, local_comm_)
       runtime_inited_ = .true.
    end if
    localCommunicator = local_comm_
  end subroutine yac_runtime_init

  ! Idempotent per grid_name: first call defines the grid via
  ! cpl_yac_define_unstr_generic (which also sets yac_fdef_datetime).
  ! Subsequent calls with the same grid_name return the cached
  ! grid_id/points_id. A different grid_name on a subsequent call is
  ! currently a fatal error — see the multi-grid future hook in the
  ! module header.
  subroutine yac_runtime_ensure_grid(grid_name, partit, mesh, grid_id, points_id)
    use MOD_MESH,   only: t_mesh
    use MOD_PARTIT, only: t_partit
    character(len=*), intent(in)            :: grid_name
    type(t_mesh),     intent(in),    target :: mesh
    type(t_partit),   intent(inout), target :: partit
    integer,          intent(out)           :: grid_id, points_id

    if (.not. grid_defined_) then
       call cpl_yac_define_unstr_generic(partit, mesh, grid_name, grid_id_, points_id_)
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

  ! Idempotent: first call runs yac_fenddef. Subsequent calls no-op.
  ! Must be called once after all interfaces have registered their fields.
  subroutine yac_runtime_enddef()
    integer :: ierr
    if (.not. enddef_done_) then
       call yac_fenddef(ierr)
       enddef_done_ = .true.
    end if
  end subroutine yac_runtime_enddef

  ! Idempotent: first call runs yac_ffinalize. Subsequent calls no-op.
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

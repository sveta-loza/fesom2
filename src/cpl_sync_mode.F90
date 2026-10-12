! Ocean <-> sea-ice exchange mode (identical module in the ocean and the sea-ice tree).
!
! The split's purpose is compute performance: ocean and sea ice run CONCURRENTLY, each
! side sends the fields of its previous step before it receives (lag 1). Decision
! 2026-10-12: the concurrent (asynchronous) exchange is the DEFAULT.
!
! Known property of the lag-1 exchange with an interactive atmosphere whose fluxes are
! piecewise constant over 3 h (IFS, config #1): the ocean applies the sea-ice-returned
! surface fluxes one coupling step late, which deepens the summer mixed layer (IFS
! tco79/CORE2 one-month test: Southern Ocean SST -0.11 K vs the monolithic model). It
! vanishes with the synchronous exchange ('all'), which reproduces the monolithic
! sequence exactly but puts the sea ice on the critical path. With prescribed forcing
! (config #2) the concurrent mode is equivalent to the monolithic model.
!
! FESOM_CPL_SYNC selects the mode (same value in both executables):
!   'async' - (default) concurrent exchange, lag 1
!   'all'   - sequential exchange every coupling step: ocean sends, sea ice
!             receives/computes/sends, ocean receives -> the monolithic sequence
!   'first' - sequential on the first step only (no cold-start step of zero fluxes),
!             concurrent afterwards
module cpl_sync_mode
  implicit none
  private
  public :: cpl_sync_step, cpl_sync_mode_name
  logical, save :: initialized = .false.
  character(len=16), save :: mode = 'async'
contains
  subroutine cpl_sync_init(mype)
    integer, intent(in) :: mype
    character(len=16) :: env
    integer :: istat
    initialized = .true.
    call get_environment_variable('FESOM_CPL_SYNC', env, status=istat)
    if (istat == 0 .and. len_trim(env) > 0) then
       mode = trim(adjustl(env))
    else
       mode = 'async'
    end if
    if (mode /= 'all' .and. mode /= 'first' .and. mode /= 'async') then
       if (mype == 0) print *, 'TerraDT: FESOM_CPL_SYNC=', trim(mode), ' not understood -> using the default ''async'''
       mode = 'async'
    end if
    if (mype == 0) print *, 'TerraDT: ocean/sea-ice exchange mode FESOM_CPL_SYNC=''', trim(mode), ''' (async = concurrent, default)'
  end subroutine cpl_sync_init

  logical function cpl_sync_step(istep, mype)
    integer, intent(in) :: istep, mype
    if (.not. initialized) call cpl_sync_init(mype)
    cpl_sync_step = (mode == 'all') .or. (mode == 'first' .and. istep == 1)
  end function cpl_sync_step

  function cpl_sync_mode_name() result(s)
    character(len=16) :: s
    s = mode
  end function cpl_sync_mode_name
end module cpl_sync_mode

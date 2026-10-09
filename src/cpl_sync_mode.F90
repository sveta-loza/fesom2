! Ocean <-> sea-ice exchange mode (identical module in the ocean and the sea-ice tree).
!
! Asynchronous exchange (each side sends the fields of its previous step before it
! receives, lag 1, both components run concurrently) was found to bias the summer
! mixed layer: the ocean applies the sea-ice-returned surface fluxes one coupling step
! late, which deepens the mixed layer irreversibly (IFS tco79/CORE2 one-month test:
! Southern Ocean SST -0.12 K vs the monolithic model, vanishes with the synchronous
! exchange). Decision 2026-10-09: synchronous is the DEFAULT.
!
! FESOM_CPL_SYNC selects the mode (same value in both executables):
!   'all'   - (default) sequential exchange every coupling step: ocean sends, sea ice
!             receives/computes/sends, ocean receives -> the monolithic sequence
!   'first' - sequential on the first step only (no cold-start step of zero fluxes),
!             asynchronous afterwards
!   'async' - the former asynchronous behaviour (lag 1, concurrent)
module cpl_sync_mode
  implicit none
  private
  public :: cpl_sync_step, cpl_sync_mode_name
  logical, save :: initialized = .false.
  character(len=16), save :: mode = 'all'
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
       mode = 'all'
    end if
    if (mode /= 'all' .and. mode /= 'first' .and. mode /= 'async') then
       if (mype == 0) print *, 'TerraDT: FESOM_CPL_SYNC=', trim(mode), ' not understood -> using the default ''all'''
       mode = 'all'
    end if
    if (mype == 0) print *, 'TerraDT: ocean/sea-ice exchange mode FESOM_CPL_SYNC=''', trim(mode), ''' (all = synchronous, default)'
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

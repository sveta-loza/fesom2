! TerraDT diagnostic switch (2026-10-08), identical in the ocean and the sea-ice tree.
!
! The ocean/sea-ice YAC exchange is asynchronous: each side first sends the fields
! of its previous step and then receives, so the partner's data are always one
! coupling step old (lag 1) and the two components run concurrently. A cold start
! therefore gives the ocean one step of zero surface fluxes before the sea ice has
! computed anything. FESOM_CPL_SYNC selects a sequential (synchronous) exchange:
!   'all'   - every coupling step: ocean sends, sea ice receives/computes/sends,
!             ocean receives (monolithic sequence; no lag, no concurrency)
!   'first' - only on the first step (removes the cold-start step of missing
!             fluxes, keeps the lag afterwards)
!   unset   - unchanged asynchronous behaviour
module cpl_sync_mode
  implicit none
  private
  public :: cpl_sync_step, cpl_sync_mode_name
  logical, save :: initialized = .false.
  character(len=16), save :: mode = ''
contains
  subroutine cpl_sync_init(mype)
    integer, intent(in) :: mype
    character(len=16) :: env
    integer :: istat
    initialized = .true.
    call get_environment_variable('FESOM_CPL_SYNC', env, status=istat)
    if (istat == 0) then
       mode = trim(adjustl(env))
    else
       mode = ''
    end if
    if (mode /= 'all' .and. mode /= 'first' .and. mode /= '') then
       if (mype == 0) print *, 'TerraDT: FESOM_CPL_SYNC=', trim(mode), ' not understood -> ignored (async)'
       mode = ''
    end if
    if (mype == 0) print *, 'TerraDT: ocean/sea-ice exchange mode FESOM_CPL_SYNC=''', trim(mode), ''' (empty = async, lag 1)'
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

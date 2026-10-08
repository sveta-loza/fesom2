! TerraDT diagnostic switch (2026-10-08), ocean side of the ocean/sea-ice split.
!
! With the sea ice as a separate YAC component (src_lag = tgt_lag = 1) the ocean
! receives the thermodynamic surface fluxes (heat, fresh water, ice-ocean stress,
! ice fraction) one coupling step after the atmospheric state they were computed
! from, while it applies the atmosphere->ocean momentum flux of the current step
! immediately. Surface buoyancy and momentum forcing are therefore out of phase by
! one coupling step. FESOM_CPL_LAG_ATM_STRESS=1 delays stress_atmoce_x/y by one
! step around the surface-flux stage so that everything the ocean applies refers
! to the same instant (the one the sea ice saw). Off by default: no change.
module cpl_atm_stress_lag
  use o_PARAM, only: WP
  implicit none
  private
  public :: atm_stress_lag_begin, atm_stress_lag_end
  logical, save :: initialized = .false., enabled = .false., have_prev = .false.
  real(kind=WP), allocatable, save :: cur_x(:), cur_y(:), prev_x(:), prev_y(:)
contains
  subroutine atm_stress_lag_begin(mype)
    use o_ARRAYS, only: stress_atmoce_x, stress_atmoce_y
    integer, intent(in) :: mype
    character(len=32) :: env
    integer :: istat
    if (.not. initialized) then
       initialized = .true.
       call get_environment_variable('FESOM_CPL_LAG_ATM_STRESS', env, status=istat)
       enabled = (istat == 0 .and. len_trim(env) > 0 .and. env(1:1) /= '0')
       if (enabled) then
          allocate(cur_x(size(stress_atmoce_x)), cur_y(size(stress_atmoce_y)), &
                   prev_x(size(stress_atmoce_x)), prev_y(size(stress_atmoce_y)))
       end if
       if (mype == 0) print *, 'TerraDT: atm->ocean stress lagged by one coupling step (FESOM_CPL_LAG_ATM_STRESS) = ', enabled
    end if
    if (.not. enabled) return
    cur_x = stress_atmoce_x
    cur_y = stress_atmoce_y
    if (have_prev) then
       stress_atmoce_x = prev_x
       stress_atmoce_y = prev_y
    end if
  end subroutine atm_stress_lag_begin

  subroutine atm_stress_lag_end()
    use o_ARRAYS, only: stress_atmoce_x, stress_atmoce_y
    if (.not. enabled) return
    stress_atmoce_x = cur_x
    stress_atmoce_y = cur_y
    prev_x = cur_x
    prev_y = cur_y
    have_prev = .true.
  end subroutine atm_stress_lag_end
end module cpl_atm_stress_lag

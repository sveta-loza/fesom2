module ice_initial_state_interface
    interface
        subroutine ice_initial_state(ice, tracers, partit, mesh)
        USE MOD_ICE
        USE MOD_TRACER
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine ice_initial_state
    end interface
end module ice_initial_state_interface

module ice_setup_interface
    interface
        subroutine ice_setup(dynamics, ice, tracers, partit, mesh)
        USE MOD_ICE
        USE MOD_TRACER
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        use MOD_DYN
        type(t_ice)   , intent(inout), target :: ice
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_tracer), intent(inout)   , target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine ice_setup
    end interface
end module ice_setup_interface

module ice_timestep_interface
    interface
        subroutine ice_timestep(istep, ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        integer       , intent(in)            :: istep
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine ice_timestep
    end interface
end module ice_timestep_interface

module tracer_init_interface
    interface
        subroutine tracer_init(tracers, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in)  ,  target :: mesh
        end subroutine tracer_init
    end interface
end module tracer_init_interface

module dynamics_init_interface
    interface
        subroutine dynamics_init(dynamics, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        use MOD_DYN
        type(t_dyn)   , intent(inout), target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine dynamics_init
    end interface
end module dynamics_init_interface

!
!_______________________________________________________________________________
! ice initialization + array allocation + time stepping
subroutine ice_setup(dynamics, ice, tracers, partit, mesh)
    USE MOD_ICE
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use MOD_DYN
    use o_param
    use g_CONFIG
    use ice_initial_state_interface
    use ice_init_ale_interface
    use ice_fct_interfaces
    implicit none
    type(t_ice)   , intent(inout),   target :: ice
    type(t_tracer), intent(inout),   target :: tracers
    type(t_dyn), intent(inout), target :: dynamics
    type(t_mesh)  , intent(inout),   target :: mesh
    type(t_partit), intent(inout),   target :: partit
    !___________________________________________________________________________
    integer                               :: i, n

    !___________________________________________________________________________
    ! initialize arrays for ALE
    if (partit%mype==0) then
       write(*,*) '____________________________________________________________'
       write(*,*) ' --> sparse SSH stiff matrix'
       write(*,*)
    end if

!sl    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_ale'//achar(27)//'[0m'
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_ale'//achar(27)//'[0m'
!sl       call init_ale(dynamics, partit, mesh)
    call init_ale_ice(dynamics, partit, mesh)
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call init_stiff_mat_ale'//achar(27)//'[0m'
    call init_stiff_mat_ale(partit, mesh) !!PS test  
    

    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call muscl_init'//achar(27)//'[0m'    
!sl    call oce_adv_tra_fct_init(tracers%work, partit, mesh)
    call muscl_adv_init(tracers%work, partit, mesh) !!PS test    

    !___________________________________________________________________________
    ! initialise ice derived type
!sl    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_init'//achar(27)//'[0m'
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_init'//achar(27)//'[0m'
    call ice_init(ice, partit, mesh)

    !___________________________________________________________________________
    ! DO not change
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> calc ice_dt'//achar(27)//'[0m'
    ice%ice_dt   = real(ice%ice_ave_steps,WP)*dt
    ! ice_dt=dt
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> calc ice_Tevp_inv'//achar(27)//'[0m'
    ice%Tevp_inv = 3.0_WP/ice%ice_dt
    ! This is combination it always enters
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> calc ice_Clim_evp'//achar(27)//'[0m'
    ice%Clim_evp = ice%Clim_evp*(ice%evp_rheol_steps/ice%ice_dt)**2/ice%Tevp_inv

    !___________________________________________________________________________
!sl    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_init'//achar(27)//'[0m'
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_init'//achar(27)//'[0m'
    call ice_mass_matrix_fill(ice, partit, mesh)

    !___________________________________________________________________________
    ! Initialization routine, user input is required
    !call ice_init_fields_test
!sl    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_initial_state'//achar(27)//'[0m'
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_initial_state'//achar(27)//'[0m'
    call ice_initial_state(ice, tracers, partit, mesh)   ! Use it unless running test example

    if(partit%mype==0) write(*,*) 'Ice is initialized'
end subroutine ice_setup
!
!
!_______________________________________________________________________________
! Sea ice model step
subroutine ice_timestep(step, ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_CONFIG
    use ice_EVPdynamics_interface
    use ice_maEVPdynamics_interface
    use ice_fct_interfaces
    use ice_thermodynamics_interfaces
    use cavity_interfaces
#if defined (__icepack)
    use icedrv_main,   only: step_icepack
#endif
#if defined (FESIM_PROFILING)
    use fesim_profiler
#endif
    implicit none
    integer       , intent(in)            :: step
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: i
    REAL(kind=WP)                         :: t0,t1, t2, t3
#if defined (__icepack)
    real(kind=WP)                         :: time_evp, time_advec, time_therm
#endif
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
    !LA 2023-03-08
    real(kind=WP), dimension(:), pointer  :: u_ice_ib, v_ice_ib
#if defined (__oifs) || defined (__ifsinterface)
    real(kind=WP), dimension(:), pointer  :: a_ice, ice_temp
    !LA 2023-03-08
    real(kind=WP), dimension(:), pointer  :: a_ice_ib
#endif
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

!---------------------------------------------
      flag_debug =.true. !sl
  if (ib_async_mode == 0) then
      u_ice    => ice%uice(:)
      v_ice    => ice%vice(:)
  else
!$omp parallel sections num_threads(2)
! kh 19.02.21 support "first touch" idea
!$omp section
      u_ice    => ice%uice(:)
      v_ice    => ice%vice(:)
      u_ice    = 0._WP
      v_ice    = 0._WP
!$omp section
      if (use_icebergs) then
        if (allocated(ice%uice_ib)) then
          u_ice_ib => ice%uice_ib(:)
          u_ice_ib = 0._WP
        end if
        if (allocated(ice%vice_ib)) then
          v_ice_ib => ice%vice_ib(:)
          v_ice_ib = 0._WP
        end if
      end if
!$omp end parallel sections
  end if
!---------------------------------------------
#if defined (__oifs) || defined (__ifsinterface)
    a_ice    => ice%data(1)%values(:)    
    ice_temp => ice%data(4)%values(:)
#endif
    !___________________________________________________________________________
    t0=MPI_Wtime()
#if defined (FESIM_PROFILING)
    call fesim_profiler_start("ice_dynamics")
#endif
#if defined (__icepack)
    call step_icepack(ice, mesh, time_evp, time_advec, time_therm) ! EVP, advection and thermodynamic parts
#else

    !$ACC UPDATE DEVICE (ice%work%fct_massmatrix) &
    !$ACC DEVICE (ice%delta_min, ice%Tevp_inv, ice%cd_oce_ice) &
    !$ACC DEVICE (ice%work%fct_tmax, ice%work%fct_tmin) &
    !$ACC DEVICE (ice%work%fct_fluxes, ice%work%fct_plus, ice%work%fct_minus) &
    !$ACC DEVICE (ice%work%eps11, ice%work%eps12, ice%work%eps22) &
    !$ACC DEVICE (ice%work%sigma11, ice%work%sigma12, ice%work%sigma22) &
    !$ACC DEVICE (ice%work%ice_strength, ice%stress_atmice_x, ice%stress_atmice_y) &
    !$ACC DEVICE (ice%thermo%rhosno, ice%thermo%rhoice, ice%thermo%inv_rhowat) &
    !$ACC DEVICE (ice%srfoce_ssh, ice%pstar, ice%c_pressure) &
    !$ACC DEVICE (ice%work%inv_areamass, ice%work%inv_mass, ice%uice_rhs, ice%vice_rhs) &
    !$ACC DEVICE (ice%uice, ice%vice, ice%srfoce_u, ice%srfoce_v, ice%uice_old, ice%vice_old) &
    !$ACC DEVICE (ice%data(1)%values, ice%data(2)%values, ice%data(3)%values) &
    !$ACC DEVICE (ice%data(1)%valuesl, ice%data(2)%valuesl, ice%data(3)%valuesl) &
    !$ACC DEVICE (ice%data(1)%dvalues, ice%data(2)%dvalues, ice%data(3)%dvalues) &
    !$ACC DEVICE (ice%data(1)%values_rhs, ice%data(2)%values_rhs, ice%data(3)%values_rhs) &
    !$ACC DEVICE (ice%data(1)%values_div_rhs, ice%data(2)%values_div_rhs, ice%data(3)%values_div_rhs)
#if defined (__oifs) || defined (__ifsinterface)
    !$ACC UPDATE DEVICE (ice%data(4)%values, ice%data(4)%valuesl, ice%data(4)%dvalues, ice%data(4)%values_rhs, ice%data(4)%values_div_rhs)
#endif
    !___________________________________________________________________________
    ! ===== Dynamics

    SELECT CASE (ice%whichEVP)
    CASE (0)
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics...'//achar(27)//'[0m'
#if defined(_CRAYFTN)
	!dir$ noinline
#endif
        call EVPdynamics  (ice, partit, mesh)
    CASE (1)
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics_m...'//achar(27)//'[0m'
        call EVPdynamics_m(ice, partit, mesh)
    CASE (2)
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics_a...'//achar(27)//'[0m'
        call EVPdynamics_a(ice, partit, mesh)
    CASE DEFAULT
        if (mype==0) write(*,*) 'a non existing EVP scheme specified!'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    END SELECT

    if (use_cavity) call cavity_ice_clean_vel(ice, partit, mesh)
    t1=MPI_Wtime()
#if defined (FESIM_PROFILING)
    call fesim_profiler_end("ice_dynamics")
    call fesim_profiler_start("ice_advection")
#endif

    !___________________________________________________________________________
    ! ===== Advection part
    ! old FCT routines
    ! call ice_TG_rhs
    ! call ice_fct_solve
    ! call cut_off
    ! new FCT routines from Sergey Danilov 08.05.2018
#if defined (__oifs) || defined (__ifsinterface)
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else
!$ACC parallel loop present(ice_temp, a_ice)
#endif
    do i=1,myDim_nod2D+eDim_nod2D
        ice_temp(i) = ice_temp(i)*a_ice(i)
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
!$ACC END parallel loop
#endif
#endif /* (__oifs) */
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_TG_rhs_div...'//achar(27)//'[0m'
    call ice_TG_rhs    (ice, partit, mesh)

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_solve...'//achar(27)//'[0m'
    call ice_fct_solve     (ice, partit, mesh)

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_update_for_div...'//achar(27)//'[0m'
!   call ice_update_for_div(ice, partit, mesh)

    !$ACC UPDATE HOST (ice%work%fct_massmatrix) &
    !$ACC HOST (ice%delta_min, ice%Tevp_inv, ice%cd_oce_ice) &
    !$ACC HOST (ice%work%fct_tmax, ice%work%fct_tmin) &
    !$ACC HOST (ice%work%fct_fluxes, ice%work%fct_plus, ice%work%fct_minus) &
    !$ACC HOST (ice%work%eps11, ice%work%eps12, ice%work%eps22) &
    !$ACC HOST (ice%work%sigma11, ice%work%sigma12, ice%work%sigma22) &
    !$ACC HOST (ice%work%ice_strength, ice%stress_atmice_x, ice%stress_atmice_y) &
    !$ACC HOST (ice%thermo%rhosno, ice%thermo%rhoice, ice%thermo%inv_rhowat) &
    !$ACC HOST (ice%srfoce_ssh, ice%pstar, ice%c_pressure) &
    !$ACC HOST (ice%work%inv_areamass, ice%work%inv_mass, ice%uice_rhs, ice%vice_rhs) &
    !$ACC HOST (ice%uice, ice%vice, ice%srfoce_u, ice%srfoce_v, ice%uice_old, ice%vice_old) &
    !$ACC HOST (ice%data(1)%values, ice%data(2)%values, ice%data(3)%values) &
    !$ACC HOST (ice%data(1)%valuesl, ice%data(2)%valuesl, ice%data(3)%valuesl) &
    !$ACC HOST (ice%data(1)%dvalues, ice%data(2)%dvalues, ice%data(3)%dvalues) &
    !$ACC HOST (ice%data(1)%values_rhs, ice%data(2)%values_rhs, ice%data(3)%values_rhs) &
    !$ACC HOST (ice%data(1)%values_div_rhs, ice%data(2)%values_div_rhs, ice%data(3)%values_div_rhs)
#if defined (__oifs) || defined (__ifsinterface)
    !$ACC UPDATE HOST (ice%data(4)%values, ice%data(4)%valuesl, ice%data(4)%dvalues, ice%data(4)%values_rhs, ice%data(4)%values_div_rhs)
#endif

#if defined (__oifs) || defined (__ifsinterface)
!$OMP PARALLEL DO
    do i=1,myDim_nod2D+eDim_nod2D
        if (a_ice(i)>0.0_WP) ice_temp(i) = ice_temp(i)/max(a_ice(i), 1.e-6_WP)
    end do
!$OMP END PARALLEL DO
#endif /* (__oifs) */

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call cut_off...'//achar(27)//'[0m'
    call cut_off(ice, partit, mesh)

    if (use_cavity) call cavity_ice_clean_ma(ice, partit, mesh)
    t2=MPI_Wtime()
#if defined (FESIM_PROFILING)
    call fesim_profiler_end("ice_advection")
    call fesim_profiler_start("ice_thermodynamics")
#endif

    !___________________________________________________________________________
    ! ===== Thermodynamic part
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call thermodynamics...'//achar(27)//'[0m'
    call thermodynamics(ice, partit, mesh)
#endif /* (__icepack) */


    !___________________________________________________________________________
!$OMP PARALLEL DO
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> after thermodynamics...'//achar(27)//'[0m'
    do i=1,myDim_nod2D+eDim_nod2D
        ice%h_ice(i) =ice%data(2)%values(i)/max(ice%data(1)%values(i), 1.e-3)
        ice%h_snow(i)=ice%data(3)%values(i)/max(ice%data(1)%values(i), 1.e-3)
        if ( ( U_ice(i)/=0.0_WP .and. mesh%ulevels_nod2d(i)>1) .or. (V_ice(i)/=0.0_WP .and. mesh%ulevels_nod2d(i)>1) ) then
            write(*,*) " --> found cavity velocity /= 0.0_WP , ", mype
            write(*,*) " ulevels_nod2d(n) = ", mesh%ulevels_nod2d(i)
            write(*,*) " U_ice(n) = ", U_ice(i)
            write(*,*) " V_ice(n) = ", V_ice(i)
            write(*,*)
        end if
    end do
!$OMP END PARALLEL DO
    t3=MPI_Wtime()
#if defined (FESIM_PROFILING)
    call fesim_profiler_end("ice_thermodynamics")
#endif
    rtime_ice = rtime_ice + (t3-t0)
    rtime_tot = rtime_tot + (t3-t0)
    if(mod(step,logfile_outfreq)==0 .and. mype==0) then
        write(*,*) '___ICE STEP EXECUTION TIMES____________________________'
#if defined (__icepack)
        write(*,"(A, ES10.3)") '	Ice Dyn.        :', time_evp
                write(*,"(A, ES10.3)") '        Ice Advect.     :', time_advec
                write(*,"(A, ES10.3)") '        Ice Thermodyn.  :', time_therm
#else
        write(*,"(A, ES10.3)") '	Ice Dyn.        :', t1-t0
        write(*,"(A, ES10.3)") '	Ice Advect.     :', t2-t1
        write(*,"(A, ES10.3)") '	Ice Thermodyn.  :', t3-t2
#endif /* (__icepack) */
        write(*,*) '   _______________________________'
        write(*,"(A, ES10.3)") '	Ice TOTAL       :', t3-t0
        write(*,*)
     endif
     flag_debug =.false.   !sl
end subroutine ice_timestep
!
!
!_______________________________________________________________________________
! sets inital values or reads restart file for ice model
!sl it uses SST originally as tracers
!sl consider using SST
subroutine ice_initial_state(ice, tracers, partit, mesh)
    USE MOD_ICE
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use o_arrays
    use g_CONFIG
    USE g_read_other_NetCDF, only: read_other_NetCDF
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_tracer), intent(in)   , target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: i
    character(MAX_PATH)                   :: filename
    real(kind=WP), external               :: TFrez  ! Sea water freeze temperature.
!============== namelistatmdata variables ================
   integer, save                                :: nm_ic_unit     = 107 ! unit to open namelist file
   integer                                      :: iost                 !I/O status
   integer, parameter                           :: ic_max=10
   logical                                      :: ic_cyclic=.true.
   integer,             save                    :: n_ic2d
   integer,             save, dimension(ic_max) :: idlist
   character(MAX_PATH), save, dimension(ic_max) :: filelist
   logical                                      :: ini_ice_from_file=.false., file_exist=.false.
   character(50),       save, dimension(ic_max) :: varlist
   integer                                      :: current_tracer
   namelist / tracer_init2d / n_ic2d, idlist, filelist, varlist, ini_ice_from_file

    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
    !LA 2023-03-07
    real(kind=WP), dimension(:), pointer  :: a_ice_ib, m_ice_ib
    real(kind=WP), dimension(:), pointer  :: u_ice_ib, v_ice_ib
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

! LA: 2023-01-31 add asynchronous icebergs
!---------------------------------------------
m_snow       => ice%data(3)%values(:)    
m_snow=0._WP
if (.not.use_icebergs) then
    u_ice        => ice%uice(:)
    v_ice        => ice%vice(:)
    a_ice        => ice%data(1)%values(:)
    m_ice        => ice%data(2)%values(:)
    m_snow       => ice%data(3)%values(:)
    !___________________________________________________________________________
    m_ice =0._WP
    a_ice =0._WP
    u_ice =0._WP
    v_ice =0._WP
else
  if (ib_async_mode == 0) then
    u_ice        => ice%uice(:)
    v_ice        => ice%vice(:)
    a_ice        => ice%data(1)%values(:)
    m_ice        => ice%data(2)%values(:)
    m_ice        = 0._WP
    a_ice        = 0._WP
    u_ice        = 0._WP
    v_ice        = 0._WP
    if (allocated(ice%uice_ib)) then
        u_ice_ib     => ice%uice_ib(:)
        u_ice_ib     = 0._WP
    end if
    if (allocated(ice%vice_ib)) then
        v_ice_ib     => ice%vice_ib(:)
        v_ice_ib     = 0._WP
    end if
    a_ice_ib     => ice%data(size(ice%data)-1)%values(:)
    a_ice_ib     = 0._WP
    m_ice_ib     => ice%data(size(ice%data))%values(:)
    m_ice_ib     = 0._WP
  else
! kh 19.02.21 support "first touch" idea
!$omp parallel sections num_threads(2)
!$omp section
      !allocate(m_ice(n_size), a_ice(n_size))
      !do i = 1, n_size
      !    m_ice(i) = 0._WP
      !    a_ice(i) = 0._WP
      !end do
      u_ice        => ice%uice(:)
      v_ice        => ice%vice(:)
      a_ice        => ice%data(1)%values(:)
      m_ice        => ice%data(2)%values(:)
      !allocate(m_ice(n_size), a_ice(n_size))
      !allocate(m_ice_ib(n_size), a_ice_ib(n_size))
      m_ice        = 0._WP
      a_ice        = 0._WP
      u_ice        = 0._WP
      v_ice        = 0._WP
!$omp section
      !allocate(m_ice_ib(n_size), a_ice_ib(n_size))
      !do i = 1, n_size
      !    m_ice_ib(i) = 0._WP
      !    a_ice_ib(i) = 0._WP
      !end do
      u_ice_ib     => ice%uice_ib(:)
      v_ice_ib     => ice%vice_ib(:)
      a_ice_ib     => ice%data(size(ice%data)-1)%values(:)
      m_ice_ib     => ice%data(size(ice%data))%values(:)
      !allocate(m_ice(n_size), a_ice(n_size))
      !allocate(m_ice_ib(n_size), a_ice_ib(n_size))
      u_ice_ib     = 0._WP
      v_ice_ib     = 0._WP
      m_ice_ib     = 0._WP
      a_ice_ib     = 0._WP
!$omp end parallel sections
  end if
end if
! LA: 2023-01-31 add asynchronous icebergs
!---------------------------------------------


    !___________________________________________________________________________
    ! OPEN and read namelist for I/O
    open( unit=nm_ic_unit, file='namelist.tra', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
        if (mype==0) WRITE(*,*) '     file   : ', 'namelist.tra',' open ok'
        else
        if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist.tra',' ; iostat=',iost
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    end if
    read(nm_ic_unit, nml=tracer_init2d,   iostat=iost)
    close(nm_ic_unit)
    
    !
    !
    !___________________________________________________________________________
    ! switch for making sea-ice initialisation from regular gridded files and 
    ! do interpolation to fesom grid or to initialise them with a constant value
    if (.not. ini_ice_from_file) then
        if(mype==0) write(*,*) 'initialize the sea ice: cold start'
        !___________________________________________________________________________
        do i=1,myDim_nod2D+eDim_nod2D
            !_______________________________________________________________________
            ! if cavity, no sea ice, no initial state
            if (ulevels_nod2d(i)>1) cycle 
            
            !_______________________________________________________________________
            !sl consider SST that is passed from ocean
            !sl how would it relate to tracers
            if (tracers%data(1)%values(1,i)< 0.0_WP) then
                if (geo_coord_nod2D(2,i)>0._WP) then
                    m_ice(i) = 1.0_WP
                    m_snow(i)= 0.1_WP
                else
                    m_ice(i) = 2.0_WP
                    m_snow(i)= 0.5_WP
                end if
                a_ice(i) = 0.9_WP
                u_ice(i) = 0.0_WP
                v_ice(i) = 0.0_WP
            endif
        enddo
        
    else ! --> if (.not. ini_ice_from_file) then
        if (mype==0) write(*,*) 'initialize the sea ice: from file'
        do i=1, n_ic2d
            do current_tracer=1, ice%num_itracers
                if (ice%data(current_tracer)%ID==idlist(i)) then
                    
                    !___________________________________________________________
                    ! check if regular gridded sea-ice initialisation file exists
                    ! if not throw error message
                    file_exist=.False.
                    inquire(file=trim(trim(ClimateDataPath)//trim(filelist(i))), exist=file_exist) 
                    if (file_exist) then   
                        if (mype==0) then
                            write(*,*) ' --> reading 2D variable: ', trim(varlist(i)), ' into 2D tracer ID=', current_tracer
                            write(*,*) '     from file ',trim(ClimateDataPath)//trim(filelist(i))
                        end if 
                        ! read 2d sea ice variable from file and interpoalte it to 
                        ! fesom grid 
                        call read_other_NetCDF(trim(ClimateDataPath)//trim(filelist(i)), varlist(i),  1, ice%data(current_tracer)%values(:), .false., .true., partit, mesh)
                    
                    else    
                        if (mype==0) then
                            write(*,*) '____________________________________________________________________'
                            write(*,*) ' ERROR: sea-ice initialisation file not found! '
                            write(*,*) '        ', trim(ClimateDataPath)//trim(filelist(i))
                            write(*,*) '        --> check your namelist.config (ClimateDataPath=...) and namelist.tra'
                            write(*,*) '            (&tracer_init2d'
                            write(*,*) '            ....'
                            write(*,*) '            filelist= ...'
                            write(*,*) '            ....'
                            write(*,*) '            /)'
                            write(*,*) '____________________________________________________________________'
                        end if
                        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
                        
                    end if ! --> if (file_exist) then       
                end if ! --> IF (ice%data(current_tracer)%ID==idlist(i)) then
            end do ! --> DO current_tracer=1, ice%num_itracers
        end do ! --> DO i=1, n_ic2d
    end if ! --> if (.not. ini_ice_from_file) then
end subroutine ice_initial_state
!_______________________________________________________________________________
SUBROUTINE arrays_init(num_tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE o_ARRAYS
    USE o_PARAM
    use g_comm_auto
    use g_config
    use g_forcing_arrays
    use o_mixing_kpp_mod ! KPP
    USE g_forcing_param, only: use_virt_salt
    use diagnostics,     only: ldiag_dMOC, ldiag_DVD
#if defined(__recom)
    use recom_glovar
    use recom_config
    use recom_ciso
#endif

    IMPLICIT NONE
    integer,        intent(in)            :: num_tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    integer                               :: elem_size, node_size
    integer                               :: n, nt
    !___________________________________________________________________________
    ! define dynamics namelist parameter
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
! #include "associate_mesh_ass.h"
nl              => mesh%nl

    !___________________________________________________________________________
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D

    ! ================
    ! Velocities
    ! ================
    !allocate(stress_diag(2, elem_size))!delete me
    !!PS allocate(Visc(nl-1, elem_size))
    allocate(t_star(node_size))
    allocate(qsr_c(node_size))
    ! ================
    ! elevation and its rhs
    ! ================

    ! ================
    ! Monin-Obukhov
    ! ================
    if (use_ice .and. use_momix) allocate(mo(nl,node_size),mixlength(node_size))
    if (use_ice .and. use_momix) mixlength=0.
    ! ================
    ! Vertical velocity and pressure
    ! ================
    allocate( hpressure(nl,node_size))
    allocate(bvfreq(nl,node_size),mixlay_dep(node_size),bv_ref(node_size))
    ! ================
    ! Ocean forcing arrays
    ! ================
    allocate(Tclim(nl-1,node_size), Sclim(nl-1, node_size))
    !---
    ! LA: add iceberg tracers 2023-02-08
    allocate(Tclim_ib(nl-1,node_size), Sclim_ib(nl-1, node_size))
    !---
    allocate(stress_surf(2,myDim_elem2D))    !!! Attention, it is shorter !!!
    allocate(stress_node_surf(2,node_size))
    allocate(stress_atmoce_x(node_size), stress_atmoce_y(node_size))
    allocate(relax2clim(node_size))
    allocate(heat_flux(node_size), Tsurf(node_size))
    allocate(water_flux(node_size), Ssurf(node_size))
    allocate(relax_salt(node_size))
    allocate(virtual_salt(node_size))

    allocate(heat_flux_in(node_size))
    allocate(real_salt_flux(node_size)) !PS

    ! =================
    ! =================
    ! Arrays used to organize surface forcing
    ! =================
    allocate(Tsurf_t(node_size,2), Ssurf_t(node_size,2))
    allocate(tau_x_t(node_size,2), tau_y_t(node_size,2))

    ! ================
    ! REcoM forcing arrays
    ! ================
#if defined(__recom)
    allocate(dtr_bf    ( nl-1, node_size ))
    allocate(str_bf    ( nl-1, node_size ))
    allocate(vert_sink ( nl-1, node_size ))
    allocate(Alk_surf  (       node_size ))
#endif
    ! =================
    ! Visc and Diff coefs
    ! =================

    allocate(Av(nl,elem_size), Kv(nl,node_size))

    Av=0.0_WP
    Kv=0.0_WP
    if (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
    allocate(Kv_double(nl,node_size, num_tracers))
    Kv_double=0.0_WP
    !!PS call oce_mixing_kpp_init ! Setup constants, allocate arrays and construct look up table
    end if
    ! tracer gradients & RHS  
    allocate(tr_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
    allocate(tr_z(nl,myDim_nod2D+eDim_nod2D))

    ! neutral slope etc. to be used in Redi formulation
    allocate(neutral_slope(3, nl-1, node_size))
    allocate(slope_tapered(3, nl-1, node_size))
    allocate(Ki(nl-1, node_size))

    do n=1, node_size
        !  Ki(n)=K_hor*area(1,n)/scale_area
        Ki(:,n)=K_hor*(mesh%mesh_resolution(n)/100000.0_WP)**2
    end do
    call exchange_nod(Ki, partit)

    neutral_slope=0.0_WP
    slope_tapered=0.0_WP

    allocate(MLD1(node_size), MLD2(node_size), MLD3(node_size))
    allocate(MLD1_ind(node_size), MLD2_ind(node_size), MLD3_ind(node_size))
    if (use_global_tides) then
    allocate(ssh_gp(node_size))
    ssh_gp=0.
    end if
    ! xy gradient of a neutral surface
    allocate(sigma_xy(2, nl-1, node_size))
    sigma_xy=0.0_WP
    ! alpha and beta in the EoS
    allocate(sw_beta(nl-1, node_size), sw_alpha(nl-1, node_size))
    allocate(dens_flux(node_size))
    sw_beta  =0.0_WP
    sw_alpha =0.0_WP
    dens_flux=0.0_WP

    if (Fer_GM) then
    allocate(fer_c(node_size),fer_scal(node_size), fer_gamma(2, nl, node_size), fer_K(nl, node_size))
    fer_gamma=0.0_WP
    fer_K=500._WP
    fer_c=1._WP
    fer_scal = 0.0_WP
    end if

    if (SPP) then
    allocate(ice_rejected_salt(node_size))
    ice_rejected_salt=0._WP
    end if

    ! =================
    ! Initialize with zeros 
    ! =================
    hpressure=0.0_WP
!
    heat_flux=0.0_WP
    heat_flux_in=0.0_WP
    Tsurf=0.0_WP

    water_flux=0.0_WP
    relax_salt=0.0_WP
    virtual_salt=0.0_WP

    Ssurf=0.0_WP

    real_salt_flux=0.0_WP

    stress_surf      =0.0_WP
    stress_node_surf =0.0_WP
    stress_atmoce_x  =0.0_WP
    stress_atmoce_y  =0.0_WP

    bvfreq=0.0_WP
    mixlay_dep=0.0_WP
    bv_ref=0.0_WP

    MLD1   =0.0_WP
    MLD2   =0.0_WP
    MLD1_ind=0.0_WP
    MLD2_ind=0.0_WP

    relax2clim=0.0_WP

    Tsurf_t=0.0_WP
    Ssurf_t=0.0_WP
    tau_x_t=0.0_WP
    tau_y_t=0.0_WP

! ================
! RECOM forcing arrays
! ================
#if defined(__recom)
    dtr_bf              = 0.0_WP
    str_bf              = 0.0_WP
    vert_sink           = 0.0_WP
    Alk_surf            = 0.0_WP
#endif

    ! init field for pressure force 
    allocate(density_ref(nl-1,node_size))
    density_ref = density_0
    allocate(density_m_rho0(nl-1, node_size))
    allocate(density_m_rho0_slev(nl-1, node_size)) !!PS
    if (ldiag_dMOC) then
       allocate(density_dMOC       (nl-1, node_size))
    end if
    allocate(pgf_x(nl-1, elem_size),pgf_y(nl-1, elem_size))
    density_m_rho0=0.0_WP
    density_m_rho0_slev=0.0_WP !!PS
    if (ldiag_dMOC) then
       density_dMOC       =0.0_WP
    end if
    pgf_x = 0.0_WP
    pgf_y = 0.0_WP

!!PS     ! init dummy arrays
!!PS     allocate(dum_2d_n(node_size), dum_3d_n(nl-1,node_size))
!!PS     allocate(dum_2d_e(elem_size), dum_3d_e(nl-1,elem_size)) 
!!PS     dum_2d_n = 0.0_WP
!!PS     dum_3d_n = 0.0_WP
!!PS     dum_2d_e = 0.0_WP
!!PS     dum_3d_e = 0.0_WP

    !---wiso-code
    if (lwiso) then
      allocate(tr_arr_ice(node_size,3))  ! add sea ice tracers
      allocate(wiso_flux_oce(node_size,3))
      allocate(wiso_flux_ice(node_size,3))

      ! initialize sea ice isotopes with 0. permill
      ! absolute tracer values are increased by factor 1000. for numerical reasons
      ! (see also routine oce_fluxes in ice_oce_coupling.F90)
      do nt = 1,3
         tr_arr_ice(:,nt)=wiso_smow(nt) * 1000.0_WP
      end do

      ! initialize atmospheric fluxes over open ocean and sea ice
      wiso_flux_oce=0.0_WP
      wiso_flux_ice=0.0_WP
    end if
    !---wiso-code-end
                                                           
END SUBROUTINE arrays_init
!_______________________________________________________________________________
SUBROUTINE dynamics_init(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE o_param
    IMPLICIT NONE
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_dyn)   , intent(inout), target :: dynamics
    !___________________________________________________________________________
    integer        :: elem_size, node_size
    integer, save  :: nm_unit  = 105       ! unit to open namelist file, skip 100-102 for cray
    integer        :: iost
    !___________________________________________________________________________
    ! define dynamics namelist parameter
    integer        :: opt_visc
    real(kind=WP)  :: visc_gamma0, visc_gamma1, visc_gamma2
    real(kind=WP)  :: visc_easybsreturn
    logical        :: use_ivertvisc=.true.
    logical        :: uke_scaling=.true.
    real(kind=WP)  :: uke_scaling_factor=1._WP
    logical        :: uke_advection=.false.
    real(kind=WP)  :: rosb_dis=1._WP
    integer        :: smooth_back=2
    integer        :: smooth_dis=2
    integer        :: smooth_back_tend=4
    real(kind=WP)  :: K_back=600._WP
    real(kind=WP)  :: c_back=0.1_8
    integer        :: momadv_opt
    logical        :: use_freeslip =.false.
    logical        :: use_wsplit   =.false.
    logical        :: ldiag_KE     =.false.
    integer        :: AB_order     = 2
    logical        :: check_opt_visc=.true.
    real(kind=WP)  :: wsplit_maxcfl
    logical        :: use_ssh_se_subcycl=.false.
    integer        :: se_BTsteps
    real(kind=WP)  :: se_BTtheta
    logical        :: se_visc, se_bottdrag, se_bdrag_si
    real(kind=WP)  :: se_visc_gamma0, se_visc_gamma1, se_visc_gamma2

    namelist /dynamics_visc   / opt_visc, check_opt_visc, visc_gamma0, visc_gamma1, visc_gamma2,  &
                                use_ivertvisc, visc_easybsreturn, &
                                uke_scaling, uke_scaling_factor, uke_advection, &
                                rosb_dis, smooth_back, smooth_dis, smooth_back_tend, K_back, c_back

    namelist /dynamics_general/ momadv_opt, use_freeslip, use_wsplit, wsplit_maxcfl, &
                                ldiag_KE, AB_order,                                  &
                                use_ssh_se_subcycl, se_BTsteps, se_BTtheta,          &
                                se_bottdrag, se_bdrag_si, se_visc, se_visc_gamma0,   &
                                se_visc_gamma1, se_visc_gamma2

    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
! #include "associate_mesh_ass.h"
nl => mesh%nl

    !___________________________________________________________________________
    ! open and read namelist for I/O
    open(unit=nm_unit, file='namelist.dyn', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
        if (mype==0) write(*,*) '     file   : ', 'namelist.dyn',' open ok'
    else
        if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ', 'namelist.dyn',' ; iostat=',iost
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    end if
    read(nm_unit, nml=dynamics_visc,    iostat=iost)
    read(nm_unit, nml=dynamics_general, iostat=iost)
    close(nm_unit)

    !___________________________________________________________________________
    ! set parameters in derived type
    dynamics%opt_visc            = opt_visc
    dynamics%check_opt_visc      = check_opt_visc
    dynamics%visc_gamma0         = visc_gamma0
    dynamics%visc_gamma1         = visc_gamma1
    dynamics%visc_gamma2         = visc_gamma2
    dynamics%visc_easybsreturn   = visc_easybsreturn
    dynamics%uke_scaling         = uke_scaling
    dynamics%uke_scaling_factor  = uke_scaling_factor
    dynamics%uke_advection       = uke_advection
    dynamics%rosb_dis            = rosb_dis
    dynamics%smooth_back         = smooth_back
    dynamics%smooth_dis          = smooth_dis
    dynamics%smooth_back_tend    = smooth_back_tend
    dynamics%K_back              = K_back
    dynamics%c_back              = c_back
    dynamics%use_ivertvisc       = use_ivertvisc
    dynamics%momadv_opt          = momadv_opt
    dynamics%use_freeslip        = use_freeslip
    dynamics%use_wsplit          = use_wsplit
    dynamics%wsplit_maxcfl       = wsplit_maxcfl
    dynamics%ldiag_KE            = ldiag_KE
    dynamics%AB_order            = AB_order
    dynamics%use_ssh_se_subcycl  = use_ssh_se_subcycl
    if (dynamics%use_ssh_se_subcycl) then
        dynamics%se_BTsteps      = se_BTsteps
        dynamics%se_BTtheta      = se_BTtheta
        dynamics%se_bottdrag     = se_bottdrag
        dynamics%se_bdrag_si     = se_bdrag_si
        dynamics%se_visc         = se_visc
        dynamics%se_visc_gamma0  = se_visc_gamma0
        dynamics%se_visc_gamma1  = se_visc_gamma1
        dynamics%se_visc_gamma2  = se_visc_gamma2
        if (mype==0) then
            write(*,*) " ___Split-Explicit barotropic subcycling_________", dynamics%se_BTsteps
            write(*,*) "     se_BTsteps     = ", dynamics%se_BTsteps
            write(*,*) "     se_BTtheta     = ", dynamics%se_BTtheta
            write(*,*) "     se_bottdrag    = ", dynamics%se_bottdrag
            write(*,*) "     se_bdrag_si    = ", dynamics%se_bdrag_si
            write(*,*) "     se_visc        = ", dynamics%se_visc
            write(*,*) "     se_visc_gamma0 = ", dynamics%se_visc_gamma0
            write(*,*) "     se_visc_gamma1 = ", dynamics%se_visc_gamma1
            write(*,*) "     se_visc_gamma2 = ", dynamics%se_visc_gamma2
        end if
    end if

    !___________________________________________________________________________
    ! define local vertice & elem array size
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D

    !___________________________________________________________________________
    ! allocate/initialise horizontal velocity arrays in derived type
    allocate(dynamics%uv(        2, nl-1, elem_size))
    allocate(dynamics%uv_rhs(    2, nl-1, elem_size))
    allocate(dynamics%uv_rhsAB(  dynamics%AB_order-1, 2, nl-1, elem_size))
    allocate(dynamics%uvnode(    2, nl-1, node_size))
    dynamics%uv              = 0.0_WP
    dynamics%uv_rhs          = 0.0_WP
    dynamics%uv_rhsAB        = 0.0_WP
    dynamics%uvnode          = 0.0_WP
    if (Fer_GM) then
        allocate(dynamics%fer_uv(2, nl-1, elem_size))
        dynamics%fer_uv      = 0.0_WP
    end if

    !___________________________________________________________________________
    ! allocate/initialise vertical velocity arrays in derived type
    allocate(dynamics%w(              nl, node_size))
    if (dynamics%ldiag_ke) then
       allocate(dynamics%w_old(       nl, node_size))
    end if
    allocate(dynamics%w_e(            nl, node_size))
    allocate(dynamics%w_i(            nl, node_size))
    allocate(dynamics%cfl_z(          nl, node_size))
    dynamics%w               = 0.0_WP
    dynamics%w_e             = 0.0_WP
    dynamics%w_i             = 0.0_WP
    dynamics%cfl_z           = 0.0_WP
    if (Fer_GM) then
        allocate(dynamics%fer_w(      nl, node_size))
        dynamics%fer_w       = 0.0_WP
    end if

    !___________________________________________________________________________
    ! allocate/initialise ssh arrays in derived type
    allocate(dynamics%eta_n(      node_size))
    dynamics%eta_n           = 0.0_WP
    if (dynamics%use_ssh_se_subcycl) then
        allocate(dynamics%se_uvh(       2, nl-1, elem_size))
        allocate(dynamics%se_uvBT_rhs(  2,       elem_size))
        allocate(dynamics%se_uvBT_4AB(  4,       elem_size))
        allocate(dynamics%se_uvBT(      2,       elem_size))
        allocate(dynamics%se_uvBT_theta(2,       elem_size))
        allocate(dynamics%se_uvBT_mean( 2,       elem_size))
        allocate(dynamics%se_uvBT_12(   2,       elem_size))
        dynamics%se_uvh         = 0.0_WP
        dynamics%se_uvBT_rhs    = 0.0_WP
        dynamics%se_uvBT_4AB    = 0.0_WP
        dynamics%se_uvBT        = 0.0_WP
        dynamics%se_uvBT_theta  = 0.0_WP
        dynamics%se_uvBT_mean   = 0.0_WP
        dynamics%se_uvBT_12     = 0.0_WP
        if (dynamics%se_visc) then
            allocate(dynamics%se_uvBT_stab_hvisc(2, elem_size))
            dynamics%se_uvBT_stab_hvisc = 0.0_WP
        end if
        if (dynamics%se_bottdrag) then
            allocate(dynamics%se_uvBT_stab_bdrag(elem_size))
            dynamics%se_uvBT_stab_bdrag = 0.0_WP
        end if
    else
        allocate(dynamics%d_eta(      node_size))
        allocate(dynamics%ssh_rhs(    node_size))
        dynamics%d_eta          = 0.0_WP
        dynamics%ssh_rhs        = 0.0_WP
        !!PS     allocate(dynamics%ssh_rhs_old(node_size))
        !!PS     dynamics%ssh_rhs_old= 0.0_WP   
    end if

    !___________________________________________________________________________
    ! inititalise working arrays
    allocate(dynamics%work%uvnode_rhs(2, nl-1, node_size))
    allocate(dynamics%work%u_c(nl-1, elem_size))
    allocate(dynamics%work%v_c(nl-1, elem_size))
    dynamics%work%uvnode_rhs = 0.0_WP
    dynamics%work%u_c = 0.0_WP
    dynamics%work%v_c = 0.0_WP
    if (dynamics%opt_visc==5) then
        allocate(dynamics%work%u_b(nl-1, elem_size))
        allocate(dynamics%work%v_b(nl-1, elem_size))
        dynamics%work%u_b = 0.0_WP
        dynamics%work%v_b = 0.0_WP
    end if

    if (dynamics%ldiag_ke) then
       allocate(dynamics%ke_adv    (2, nl-1, elem_size))
       allocate(dynamics%ke_cor    (2, nl-1, elem_size))
       allocate(dynamics%ke_pre    (2, nl-1, elem_size))
       allocate(dynamics%ke_hvis   (2, nl-1, elem_size))
       allocate(dynamics%ke_vvis   (2, nl-1, elem_size))
       allocate(dynamics%ke_umean  (2, nl-1, elem_size))
       allocate(dynamics%ke_u2mean (2, nl-1, elem_size))
       allocate(dynamics%ke_du2    (2, nl-1, elem_size))
       allocate(dynamics%ke_adv_AB (dynamics%AB_order-1, 2, nl-1, elem_size))
       allocate(dynamics%ke_cor_AB (dynamics%AB_order-1, 2, nl-1, elem_size))
       allocate(dynamics%ke_rhs_bak(2, nl-1, elem_size))
       allocate(dynamics%ke_wrho   (nl-1, node_size))
       allocate(dynamics%ke_dW     (nl-1, node_size))
       allocate(dynamics%ke_Pfull  (nl-1, node_size))
       allocate(dynamics%ke_wind   (2, elem_size))
       allocate(dynamics%ke_drag   (2, elem_size))

       allocate(dynamics%ke_pre_xVEL (2, nl-1, elem_size))
       allocate(dynamics%ke_adv_xVEL (2, nl-1, elem_size))
       allocate(dynamics%ke_cor_xVEL (2, nl-1, elem_size))
       allocate(dynamics%ke_hvis_xVEL(2, nl-1, elem_size))
       allocate(dynamics%ke_vvis_xVEL(2, nl-1, elem_size))
       allocate(dynamics%ke_wind_xVEL(2, elem_size))
       allocate(dynamics%ke_drag_xVEL(2, elem_size))
       allocate(dynamics%ke_J(node_size),  dynamics%ke_D(node_size),   dynamics%ke_G(node_size),  &
                dynamics%ke_D2(node_size), dynamics%ke_n0(node_size),  dynamics%ke_JD(node_size), &
                dynamics%ke_GD(node_size), dynamics%ke_swA(node_size), dynamics%ke_swB(node_size))

       dynamics%ke_adv      =0.0_WP
       dynamics%ke_cor      =0.0_WP
       dynamics%ke_pre      =0.0_WP
       dynamics%ke_hvis     =0.0_WP
       dynamics%ke_vvis     =0.0_WP
       dynamics%ke_du2      =0.0_WP
       dynamics%ke_umean    =0.0_WP
       dynamics%ke_u2mean   =0.0_WP
       dynamics%ke_adv_AB   =0.0_WP
       dynamics%ke_cor_AB   =0.0_WP
       dynamics%ke_rhs_bak  =0.0_WP
       dynamics%ke_wrho     =0.0_WP
       dynamics%ke_wind     =0.0_WP
       dynamics%ke_drag     =0.0_WP
       dynamics%ke_pre_xVEL =0.0_WP
       dynamics%ke_adv_xVEL =0.0_WP
       dynamics%ke_cor_xVEL =0.0_WP
       dynamics%ke_hvis_xVEL=0.0_WP
       dynamics%ke_vvis_xVEL=0.0_WP
       dynamics%ke_wind_xVEL=0.0_WP
       dynamics%ke_drag_xVEL=0.0_WP
       dynamics%ke_dW       =0.0_WP
       dynamics%ke_Pfull    =0.0_WP
       dynamics%ke_J        =0.0_WP
       dynamics%ke_D        =0.0_WP
       dynamics%ke_G        =0.0_WP
       dynamics%ke_D2       =0.0_WP
       dynamics%ke_n0       =0.0_WP
       dynamics%ke_JD       =0.0_WP
       dynamics%ke_GD       =0.0_WP
       dynamics%ke_swA      =0.0_WP
       dynamics%ke_swB      =0.0_WP
    end if
END SUBROUTINE dynamics_init
!
!
!_______________________________________________________________________________
SUBROUTINE tracer_init(tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_TRACER
    USE DIAGNOSTICS, only: ldiag_DVD
    USE g_ic3d
    use g_forcing_param, only: use_age_tracer !---age-code
    use g_config, only : lwiso, use_transit   ! add lwiso switch and switch for transient tracers
    use mod_transit, only : index_transit_r14c, index_transit_r39ar, index_transit_f11, index_transit_f12, index_transit_sf6, l_r14c, l_r39ar, l_f11, l_f12, l_sf6
    IMPLICIT NONE
    type(t_tracer), intent(inout), target               :: tracers
    type(t_partit), intent(inout), target               :: partit
    type(t_mesh),   intent(in) ,   target               :: mesh
    type(nml_tracer_list_type),    target, allocatable  :: nml_tracer_list(:)
    !___________________________________________________________________________
    integer        :: elem_size, node_size
    integer, save  :: nm_unit  = 104       ! unit to open namelist file, skip 100-102 for cray
    integer        :: iost
    integer        :: n
    !___________________________________________________________________________
    ! define tracer namelist parameter
    integer        :: num_tracers
    logical        :: i_vert_diff, smooth_bh_tra
    real(kind=WP)  :: gamma0_tra, gamma1_tra, gamma2_tra
    integer        :: AB_order = 2
    namelist /tracer_listsize/ num_tracers
    namelist /tracer_list    / nml_tracer_list
    namelist /tracer_general / smooth_bh_tra, gamma0_tra, gamma1_tra, gamma2_tra, i_vert_diff, AB_order
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
! #include "associate_mesh_ass.h"
nl => mesh%nl

    !___________________________________________________________________________
    ! OPEN and read namelist for I/O
    open( unit=nm_unit, file='namelist.tra', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
        if (mype==0) WRITE(*,*) '     file   : ', 'namelist.tra',' open ok'
        else
        if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist.tra',' ; iostat=',iost
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    end if

    READ(nm_unit, nml=tracer_listsize, iostat=iost)
    allocate(nml_tracer_list(num_tracers))
    READ(nm_unit, nml=tracer_list,     iostat=iost)
    read(nm_unit, nml=tracer_init3d,   iostat=iost)
    READ(nm_unit, nml=tracer_general,  iostat=iost)
    close(nm_unit)

    do n=1, num_tracers
    if (nml_tracer_list(n)%id==-1) then
        if (mype==0) write(*,*) 'number of tracers will be changed from ', num_tracers, ' to ', n-1, '!'
        num_tracers=n-1
        EXIT
    end if
    end do

    !---wiso-code
    !=====================
    ! set necessary water isotope variables
    !=====================
    IF (lwiso) THEN
      ! always assume 3 water isotope tracers in the order H218O, HD16O, H216O
      ! tracers simulated in the model
      nml_tracer_list(num_tracers+1) = nml_tracer_list(1) ! use the same scheme as temperature
      nml_tracer_list(num_tracers+2) = nml_tracer_list(1)
      nml_tracer_list(num_tracers+3) = nml_tracer_list(1)

      nml_tracer_list(num_tracers+1)%id = 101
      nml_tracer_list(num_tracers+2)%id = 102
      nml_tracer_list(num_tracers+3)%id = 103

      index_wiso_tracers(1) = num_tracers+1
      index_wiso_tracers(2) = num_tracers+2
      index_wiso_tracers(3) = num_tracers+3

      num_tracers = num_tracers + 3

      ! tracers initialised from file
      idlist((n_ic3d+1):(n_ic3d+3)) = (/101, 102, 103/)
      filelist((n_ic3d+1):(n_ic3d+3)) = (/'wiso.nc', 'wiso.nc', 'wiso.nc'/)
      varlist((n_ic3d+1):(n_ic3d+3))  = (/'h2o18', 'hDo16', 'h2o16'/)

      n_ic3d = n_ic3d + 3

      if (mype==0) write(*,*) '3 water isotope tracers will be used in FESOM'
    END IF
    !---wiso-code-end

    !---age-code-begin
    if (use_age_tracer) then
      ! add age tracer in the model
      nml_tracer_list(num_tracers+1) = nml_tracer_list(1)
      nml_tracer_list(num_tracers+1)%id = 100
      index_age_tracer = num_tracers+1
      num_tracers = num_tracers + 1

      if (mype==0) write(*,*) '1 water age tracer will be used in FESOM'
    endif
    !---age-code-end

    ! Transient tracers
    if (use_transit) then
      ! add transient tracers to the model
      if (l_sf6) then
        nml_tracer_list(num_tracers+1) = nml_tracer_list(1)
        nml_tracer_list(num_tracers+1)%id = 6
        index_transit_sf6 = num_tracers+1
        num_tracers = num_tracers + 1
      endif

      if (l_f11) then
        nml_tracer_list(num_tracers+1) = nml_tracer_list(1)
        nml_tracer_list(num_tracers+1)%id = 11
        index_transit_f11 = num_tracers+1
        num_tracers = num_tracers + 1
      endif

      if (l_f12) then
        nml_tracer_list(num_tracers+1) = nml_tracer_list(1)
        nml_tracer_list(num_tracers+1)%id = 12
        index_transit_f12 = num_tracers+1
        num_tracers = num_tracers + 1
      endif
      if (l_r14c) then
        nml_tracer_list(num_tracers+1) = nml_tracer_list(1)
        nml_tracer_list(num_tracers+1)%id = 14
        index_transit_r14c = num_tracers+1
        num_tracers = num_tracers + 1
      endif

      if (l_r39ar) then
        nml_tracer_list(num_tracers+1) = nml_tracer_list(1)
        nml_tracer_list(num_tracers+1)%id = 39
        index_transit_r39ar = num_tracers+1
        num_tracers = num_tracers + 1
      endif

      ! tracers initialised from file
      idlist((n_ic3d+1):(n_ic3d+1)) = (/14/)
      filelist((n_ic3d+1):(n_ic3d+1)) = (/'R14C.nc'/)
      varlist((n_ic3d+1):(n_ic3d+1))  = (/'R14C'/)

      if (mype==0) write(*,*) 'XXX Transient tracers will be used in FESOM'
    endif
    ! 'use_transit' end


    if (mype==0) write(*,*) 'total number of tracers is: ', num_tracers

    !___________________________________________________________________________
    ! define local vertice & elem array size + number of tracers
    elem_size=myDim_elem2D + eDim_elem2D
    node_size=myDim_nod2D  + eDim_nod2D
    tracers%num_tracers=num_tracers

    !___________________________________________________________________________
    ! allocate/initialise horizontal velocity arrays in derived type
    ! Temperature (index=1), Salinity (index=2), etc.
    allocate(tracers%data(num_tracers))
    do n=1, tracers%num_tracers
        allocate(tracers%data(n)%values   (                             nl-1, node_size))
        allocate(tracers%data(n)%valuesAB (                             nl-1, node_size))
        tracers%data(n)%AB_order      = AB_order
        allocate(tracers%data(n)%valuesold(tracers%data(n)%AB_order-1,  nl-1, node_size))
        tracers%data(n)%ID            = nml_tracer_list(n)%id
        tracers%data(n)%tra_adv_hor   = TRIM(nml_tracer_list(n)%adv_hor)
        tracers%data(n)%tra_adv_ver   = TRIM(nml_tracer_list(n)%adv_ver)
        tracers%data(n)%tra_adv_lim   = TRIM(nml_tracer_list(n)%adv_lim)
        tracers%data(n)%tra_adv_ph    = nml_tracer_list(n)%adv_ph
        tracers%data(n)%tra_adv_pv    = nml_tracer_list(n)%adv_pv
        tracers%data(n)%smooth_bh_tra = smooth_bh_tra
        tracers%data(n)%gamma0_tra    = gamma0_tra
        tracers%data(n)%gamma1_tra    = gamma1_tra
        tracers%data(n)%gamma2_tra    = gamma2_tra
        tracers%data(n)%values        = 0.
        tracers%data(n)%valuesAB      = 0.
        tracers%data(n)%valuesold     = 0.
        tracers%data(n)%i_vert_diff   = i_vert_diff
    end do
    allocate(tracers%work%del_ttf(nl-1,node_size))
    allocate(tracers%work%del_ttf_advhoriz(nl-1,node_size),tracers%work%del_ttf_advvert(nl-1,node_size))
    tracers%work%del_ttf          = 0.0_WP
    tracers%work%del_ttf_advhoriz = 0.0_WP
    tracers%work%del_ttf_advvert  = 0.0_WP
    if (ldiag_DVD) then
        allocate(tracers%work%dvd_trflx_hor(nl-1, myDim_edge2D, 2))
        allocate(tracers%work%dvd_trflx_ver(nl  , myDim_nod2D , 2))
        tracers%work%dvd_trflx_hor = 0.0_WP
        tracers%work%dvd_trflx_ver = 0.0_WP
    end if
END SUBROUTINE tracer_init


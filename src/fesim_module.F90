! synopsis: save any derived types we initialize
!           so they can be reused after fesom_init/fesim_init
!module fesim_main_storage_module
module fesom_main_storage_module
  use iceberg_step
  USE MOD_MESH
  USE MOD_ICE
  USE MOD_TRACER
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_DYN
  USE o_ARRAYS
  USE o_PARAM
  use g_clock
  use g_config
  use g_comm_auto
  use g_forcing_arrays
  use io_RESTART
  use io_MEANDATA
  use io_mesh_info
  use diagnostics
  use mo_tidal
  use tracer_init_interface
!sl  use ocean_setup_interface
  use ice_setup_interface
  use oce_fluxes_interface
  use update_atm_forcing_interface
!sl  use before_oce_step_interface
  use oce_timestep_ale_interface
  use read_mesh_interface
  use fesom_version_info_module
! use fesim_version_info_module  
  use command_line_options_module
  use, intrinsic :: iso_fortran_env, only : real32
  use g_forcing_param, only: use_landice_water, use_age_tracer
  use landice_water_init_interface
  use age_tracer_init_interface
  use iceberg_params
  use iceberg_step
  use iceberg_ocean_coupling
  use Toy_Channel_Soufflet, only: compute_zonal_mean
  ! Define icepack module

#if defined (__icepack)
  use icedrv_main,          only: set_icepack, init_icepack, alloc_icepack
#endif

#if defined (__oasis)
  use cpl_driver
#endif
#if defined (__yac)
  use ocean_coupling_interface
#endif


!sl ! Transient tracers
!sl use mod_transit, only: year_ce, r14c_nh, r14c_tz, r14c_sh, r14c_ti, xCO2_ti, xf11_nh, xf11_sh, xf12_nh, xf12_sh, xsf6_nh, xsf6_sh, ti_transit, anthro_transit

  implicit none
    
  type :: fesom_main_storage_type

    integer           :: n, from_nstep, offset, row, i, provided, id
    integer           :: which_readr ! read which restart files (0=netcdf, 1=core dump,2=dtype)
    integer           :: total_nsteps
    integer, pointer  :: mype, npes, MPIerr, MPI_COMM_FESOM, MPI_COMM_WORLD, MPI_COMM_FESOM_IB
!sl    integer, pointer  :: mype, npes, MPIerr, MPI_COMM_FESIM, MPI_COMM_WORLD, MPI_COMM_FESIM_IB
    real(kind=WP)     :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t0_ice, t1_ice, t0_frc, t1_frc
    real(kind=WP)     :: rtime_fullice,    rtime_write_restart, rtime_write_means, rtime_compute_diag, rtime_read_forcing
    real(kind=real32) :: rtime_setup_mesh, rtime_setup_ocean, rtime_setup_forcing 
    real(kind=real32) :: rtime_setup_ice,  rtime_setup_other, rtime_setup_restart
    real(kind=real32) :: runtime_alltimesteps

    type(t_mesh)   mesh
    type(t_tracer) tracers
    type(t_dyn)    dynamics
    type(t_partit) partit
    type(t_ice)    ice


    character(LEN=256)               :: dump_dir, dump_filename
    logical                          :: L_EXISTS
    type(t_mesh)   mesh_copy
    type(t_tracer) tracers_copy
    type(t_dyn)    dynamics_copy
    type(t_ice)    ice_copy

    character(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version_txt
    integer mpi_version_len
    logical fesim_did_mpi_init
    
  end type fesom_main_storage_type
  type(fesom_main_storage_type), save, target :: f

end module fesom_main_storage_module


! synopsis: main FESIM program split into 3 parts
!           this way FESIM can e.g. be used as a library with an external time loop driver
!           used with IFS-FESOM
module fesim_module
  use fesom_main_storage_module      
#if defined  __ifsinterface
  use, intrinsic :: ieee_exceptions
#endif
  ! Enhanced profiler integration
!sl#if defined (FESOM_PROFILING)
!sl  use fesom_profiler
!sl#endif
#if defined (FESIM_PROFILING)
  use fesim_profiler
#endif
  implicit none
  public fesim_init, fesim_runloop, fesim_finalize
  private

contains
 
  subroutine fesim_init(fesim_total_nsteps)
      use fesom_main_storage_module
#if defined(__MULTIO)
      use iom
#endif
      integer, intent(out) :: fesim_total_nsteps
      ! EO parameters
      logical mpi_is_initialized
      integer              :: tr_num
#if !defined  __ifsinterface
      if(command_argument_count() > 0) then
        call command_line_options%parse()
        stop
      end if
#endif

      mpi_is_initialized = .false.
      f%fesim_did_mpi_init = .false.

#ifndef __oifs
        !ECHAM6-FESOM2 coupling: cpl_oasis3mct_init is called here in order to avoid circular dependencies between modules (cpl_driver and g_PARSUP)
        !OIFS-FESOM2 coupling: does not require MPI_INIT here as this is done by OASIS
        call MPI_Initialized(mpi_is_initialized, f%i)
        if(.not. mpi_is_initialized) then
            ! TODO: do not initialize MPI here if it has been initialized already, e.g. via IFS when fesom is called as library (__ifsinterface is defined)
            call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, f%provided, f%i)
            f%fesim_did_mpi_init = .true.
            !f%fesom_did_mpi_init = .true.
        end if
#endif

#if defined (__oasis)
        call cpl_oasis3mct_init(f%partit,f%partit%MPI_COMM_FESOM)
#elif defined (__yac)
        call ocn_cpl_init(f%partit%MPI_COMM_FESOM)
#endif

        f%t1 = MPI_Wtime()

        ! Initialize enhanced profiler
#if defined (FESIM_PROFILING)
        call fesim_profiler_init(.true.)
        call fesim_profiler_start("fesim_init_total")
#endif

#if defined (FESIM_PROFILING)
        call fesim_profiler_start("par_init")
#endif
        call par_init(f%partit)
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("par_init")
#endif

        f%mype          =>f%partit%mype
        f%MPIerr        =>f%partit%MPIerr
        f%MPI_COMM_FESOM=>f%partit%MPI_COMM_FESOM
        f%MPI_COMM_FESOM_IB=>f%partit%MPI_COMM_FESOM_IB
!sl        f%MPI_COMM_FESIM=>f%partit%MPI_COMM_FESIM
!sl        f%MPI_COMM_FESIM_IB=>f%partit%MPI_COMM_FESIM_IB        
        f%MPI_COMM_WORLD=>f%partit%MPI_COMM_WORLD

        f%npes          =>f%partit%npes


        if(f%mype==0) then
            write(*,*)
!sl TOBE further edited            
            print *,"FESOM2 git SHA: "//fesom_git_sha()
            call MPI_Get_library_version(f%mpi_version_txt, f%mpi_version_len, f%MPIERR)
            print *,"MPI library version: "//trim(f%mpi_version_txt)
            print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
            print *, achar(27)//'[7;32m'//' --> FESIM BUILDS UP MODEL CONFIGURATION                    '//achar(27)//'[0m'
        end if
        !=====================
        ! Read configuration data,  
        ! load the mesh and fill in 
        ! auxiliary mesh arrays
        !=====================
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call setup_model'//achar(27)//'[0m'
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("setup_model")
#endif
        call setup_model(f%partit)  ! Read Namelists, always before clock_init
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("setup_model")
#endif
        
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call clock_init'//achar(27)//'[0m'
        call clock_init(f%partit)   ! read the clock file
        
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call get_run_steps'//achar(27)//'[0m'
        call get_run_steps(fesim_total_nsteps, f%partit)
        f%total_nsteps=fesim_total_nsteps
#if defined (FESIM_PROFILING)
        call fesim_profiler_set_timesteps(fesim_total_nsteps)
        ! Set timestep size in seconds for SYPD calculation: 86400 seconds/day / steps_per_day
        call fesim_profiler_set_timestep_size(86400.0d0 / real(step_per_day, kind=8))
#endif
        
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call mesh_setup'//achar(27)//'[0m'
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("mesh_setup")
#endif
        call mesh_setup(f%partit, f%mesh)
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("mesh_setup")
#endif

        if (f%mype==0) write(*,*) 'FESIM mesh_setup... complete'

        !=====================
        ! Allocate field variables 
        ! and additional arrays needed for 
        ! fancy advection etc.  
        !=====================
!sl#if defined (__oasis)
!sl        !---wiso-code
!sl        IF (lwiso) THEN
!sl          nsend = nsend + 6       ! add number of water isotope tracers to coupling parameter nsend, nrecv
!sl          nrecv = nrecv + 6
!sl        END IF
!sl        !---wiso-code-end
!sl#if !defined (__oifs)
!sl        IF (use_icebergs) THEN
!sl          nrecv = nrecv + 2
!sl        END IF
!sl#endif
!sl#endif

        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call check_mesh_consistency'//achar(27)//'[0m'
        call check_mesh_consistency(f%partit, f%mesh)
        if (f%mype==0) f%t2=MPI_Wtime()

!sl        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call dynamics_init'//achar(27)//'[0m'
           if (f%mype==0)  print *, achar(27)//'[34m'//' --> call dynamics_init'//achar(27)//'[0m'
           call dynamics_init(f%dynamics, f%partit, f%mesh)
        
!sl        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call tracer_init'//achar(27)//'[0m'
           if (f%mype==0)  print *, achar(27)//'[34m'//' --> call tracer_init'//achar(27)//'[0m'
           call tracer_init(f%tracers, f%partit, f%mesh)                ! allocate array of ocean tracers (derived type "t_tracer")

!sl        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call arrays_init'//achar(27)//'[0m'
           if (f%mype==0)  print *, achar(27)//'[34m'//' --> call arrays_init'//achar(27)//'[0m'
           call arrays_init(f%tracers%num_tracers, f%partit, f%mesh)    ! allocate other arrays (to be refactured same as tracers in the future)
        
!sl        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ocean_setup'//achar(27)//'[0m'
!sl#if defined (FESOM_PROFILING)
!sl        call fesom_profiler_start("ocean_setup")
!sl        call fesom_profiler_start("dynamics_init")
!sl#endif
!sl        call ocean_setup(f%dynamics, f%tracers, f%partit, f%mesh)
!sl#if defined (FESOM_PROFILING)
!sl        call fesom_profiler_end("dynamics_init")
!sl        call fesom_profiler_end("ocean_setup")
!sl#endif

!sl        ! global tides
!sl        if (use_global_tides) then
!sl           call foreph_ini(yearnew, month, f%partit)
!sl        end if

        !_______________________________________________________________________
        ! NB the cold-start sea-ice initial state needs an ocean SST, which FESIM
        ! has none of at this point. It is now read directly from the T/S
        ! climatology inside ice_initial_state (sst_ini_file in &tracer_init2d),
        ! which avoids the ocean's 3-D vertical machinery. An earlier attempt to
        ! call do_ic3d here instead failed: it interpolates onto mesh%Z_3d_n,
        ! which only init_ale allocates, and FESIM calls init_ale_ice.
        ! [ocean-leftover audit 2026-09-11]

        call forcing_setup(f%partit, f%mesh)

        if (f%mype==0) f%t4=MPI_Wtime()
        if (use_ice) then        
!sl            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_setup'//achar(27)//'[0m'
               if (f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_setup'//achar(27)//'[0m'
               call ice_setup(f%dynamics,f%ice, f%tracers, f%partit, f%mesh)
!sl            call ice_init_toyocean_dummy(f%ice, f%partit, f%mesh)
!sl            f%ice%ice_steps_since_upd = f%ice%ice_ave_steps-1
!sl            f%ice%ice_update=.true.
            if (f%mype==0) write(*,*) 'EVP scheme option=', f%ice%whichEVP
        else 
            ! create a dummy ice derived type with only a_ice, m_ice, m_snow and 
            ! uvice since oce_timesteps still needs in moment
            ! ice as an input for mo_convect(ice, partit, mesh), call 
            ! compute_vel_rhs(ice, dynamics, partit, mesh),  
            ! call write_step_info(...) and call check_blowup(...)
            call ice_init_toyocean_dummy(f%ice, f%partit, f%mesh)            
        endif 
        
        if (f%mype==0) f%t5=MPI_Wtime()

!sl        call compute_diagnostics(0,f%ice, f%partit, f%mesh) ! allocate arrays for diagnostic
        
!sl#if defined (__oasis)

!sl        call cpl_oasis3mct_define_unstr(f%partit, f%mesh)

!sl        if(f%mype==0)  write(*,*) 'FESOM ---->     cpl_oasis3mct_define_unstr nsend, nrecv:',nsend, nrecv
!sl#endif
    

#if defined (__yac)
        call ocn_cpl_define(f%partit, f%mesh, INT(dt)*cpl_stride)   ! field dt = coupling period (every cpl_stride model steps)
        if(f%mype==0)  write(*,*) 'FESIM ---->     ocn_cpl_define OCN_NSEND, OCN_NRECV:', OCN_NSEND, OCN_NRECV
#endif

#if defined (__icepack)
        !=====================
        ! Setup icepack
        !=====================
        if (f%mype==0) write(*,*) 'Icepack: reading namelists from namelist.icepack'
        call set_icepack(f%ice, f%partit)
        call alloc_icepack
        call init_icepack(f%ice, f%tracers%data(1), f%mesh)
        if (f%mype==0) write(*,*) 'Icepack: setup complete'
#endif
        call clock_newyear                        ! check if it is a new year
        if (f%mype==0) f%t6=MPI_Wtime()
        !___READ INITIAL CONDITIONS IF THIS IS A RESTART RUN________________________
        if (r_restart) then
!sl            call read_initial_conditions(f%which_readr, f%ice, f%partit, f%mesh)
        call read_initial_conditions(f%which_readr, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
        end if
        if (f%mype==0) f%t7=MPI_Wtime()
        
        ! store grid information into netcdf file
        if (.not. r_restart) call write_mesh_info(f%partit, f%mesh)

        if (f%mype==0) then
           f%t8=MPI_Wtime()
    
           f%rtime_setup_mesh    = real( f%t2 - f%t1              ,real32)
           f%rtime_setup_ice     = real( f%t5 - f%t4              ,real32)
           f%rtime_setup_restart = real( f%t7 - f%t6              ,real32)
           f%rtime_setup_other   = real((f%t8 - f%t7) + (f%t6 - f%t5) ,real32)


           write(*,*) '=========================================='
           write(*,*) 'MODEL SETUP took on mype=0 [seconds]      '
           write(*,*) 'runtime setup total      ',real(f%t8-f%t1,real32)      
           write(*,*) ' > runtime setup mesh    ',f%rtime_setup_mesh   
           write(*,*) ' > runtime setup ice     ',f%rtime_setup_ice    
           write(*,*) ' > runtime setup restart ',f%rtime_setup_restart
!sl           write(*,*) ' > runtime setup other   ',f%rtime_setup_other
!sl#if defined (__recom)
!sl           write(*,*) ' > runtime setup recom   ',f%rtime_setup_recom
!sl#endif
            write(*,*) '============================================' 
        endif

#if defined(__MULTIO)
          call iom_send_fesom_domains(f%partit, f%mesh)
#endif

    ! Initialize timers
    f%rtime_fullice       = 0._WP
    f%rtime_write_restart = 0._WP
    f%rtime_write_means   = 0._WP
    f%rtime_compute_diag  = 0._WP
    f%rtime_read_forcing  = 0._WP

    f%from_nstep = 1

    ! End initialization profiling
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("fesim_init_total")
#endif

!sl#if defined (__oifs) || defined (__ifsinterface)
    !$ACC ENTER DATA CREATE (f%ice%data(4)%values, f%ice%data(4)%valuesl, f%ice%data(4)%dvalues, f%ice%data(4)%values_rhs, f%ice%data(4)%values_div_rhs)
!sl#endif
    !$ACC ENTER DATA COPYIN (f%dynamics)
    !$ACC ENTER DATA CREATE (f%dynamics%w, f%dynamics%w_e, f%dynamics%uv)
    !$ACC ENTER DATA CREATE (f%tracers%work%del_ttf)
    !$ACC ENTER DATA CREATE (f%tracers%data, f%tracers%work) 
!sl    do tr_num=1, f%tracers%num_tracers
    !$ACC ENTER DATA CREATE (f%tracers%data(tr_num)%values, f%tracers%data(tr_num)%valuesAB)
    !$ACC ENTER DATA CREATE (f%tracers%data(tr_num)%valuesold)
    !$ACC ENTER DATA CREATE (f%tracers%data(tr_num)%tra_adv_ph, f%tracers%data(tr_num)%tra_adv_pv)
!sl    end do
    !$ACC ENTER DATA CREATE (f%tracers%work%fct_ttf_min, f%tracers%work%fct_ttf_max, f%tracers%work%fct_plus, f%tracers%work%fct_minus)
    !$ACC ENTER DATA CREATE (f%tracers%work%adv_flux_hor, f%tracers%work%adv_flux_ver, f%tracers%work%fct_LO)
    !$ACC ENTER DATA CREATE (f%tracers%work%del_ttf_advvert, f%tracers%work%del_ttf_advhoriz, f%tracers%work%edge_up_dn_grad)
    !$ACC ENTER DATA CREATE (tr_xy, tr_z, relax2clim, Sclim, Tclim)
  end subroutine fesim_init


  subroutine fesim_runloop(current_nsteps)
    use fesom_main_storage_module
!   use openacc_lib
    integer, intent(in) :: current_nsteps 
    ! EO parameters
    integer n, nstart, ntotal, tr_num

    !=====================
    ! Time stepping
    !=====================

    ! --------------
    ! LA icebergs: 2023-05-17 
    f%MPI_COMM_FESOM_IB = f%MPI_COMM_FESOM
!sl    f%MPI_COMM_FESIM_IB = f%MPI_COMM_FESIM
    if (f%mype==0) then
!        write (*,*) 'ib_async_mode, initial omp_num_threads ', ib_async_mode, omp_get_num_threads()
        write (*,*) 'current_nsteps, steps_per_ib_step, icb_outfreq :', current_nsteps, steps_per_ib_step, icb_outfreq
    end if
    ! --------------

    if (f%mype==0) write(*,*) 'FESIM start iteration before the barrier...'
    ! Cross-component sync before starting the loop timer: both ocean and ice
    ! resume from the SAME instant, so the ocean's longer init (mesh+forcing) no
    ! longer shows up as ice first-step yac_fget wait (the ~21 s loop-timer skew).
    !
    ! ONLY valid when MPI_COMM_WORLD is exactly ocean+ice, i.e. the standalone
    ! config #2. MPI_Barrier is collective over the WHOLE communicator, so under
    ! config #0 (ICON in the same MPMD world) or config #1 (IFS likewise) the
    ! ocean and ice ranks would block here forever waiting for atmosphere ranks
    ! that never call it. Config #0/#1 therefore sync on their own component
    ! communicator only.                     [gated 2026-09-09, ICON three-way prep]
#if !defined(__yac_atm) && !defined(__ifs_fwd)
    call MPI_Barrier(MPI_COMM_WORLD, f%MPIERR)      ! [scalability shared-clock fix 2026-07]
#endif
    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)
    if (f%mype==0) then
       write(*,*) 'FESIM start iteration after the barrier...'
       f%t0 = MPI_Wtime()
    endif
    if(f%mype==0) then
        write(*,*)
        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;32m'//' --> FESIM STARTS TIME LOOP                                 '//achar(27)//'[0m'
    end if
    
    ! Start main time loop profiling
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("fesim_runloop_total")
#endif
    !___MODEL TIME STEPPING LOOP________________________________________________
    nstart=f%from_nstep
    ntotal=f%from_nstep-1+current_nsteps

    do n=nstart, ntotal
        mstep = n
        if (mod(n,logfile_outfreq)==0 .and. f%mype==0) then
            write(*,*) 'FESIM ======================================================='
!             write(*,*) 'FESIM step:',n,' day:', n*dt/24./3600.,
            write(*,*) 'FESIM step:',n,' day:', daynew,' year:',yearnew 
            write(*,*)
        end if
!sl?? should it be:
!sl #if defined (__oifs) || defined (__oasis) || defined (__yac)
#if defined (__oifs) || defined (__oasis)
            seconds_til_now=INT(dt)*(n-1)
#endif
        call clock      
        ! --------------
        !___model sea-ice step__________________________________________________
        f%t1 = MPI_Wtime()
!SL        if(use_ice) then
            ! No ocean2ice call here: in the decoupled design the ocean surface
            ! state (SST, SSS, SSH, surface u/v) arrives over YAC in
            ! update_atm_forcing_yac below, and the ice->ocean fluxes go back the
            ! same way. ocean2ice/oce_fluxes were deleted 2026-09-11 -- see the
            ! header of ice_oce_coupling_fesim.F90.
            !___compute update of atmospheric forcing____________________________
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call update_atm_forcing(n)'//achar(27)//'[0m'
            f%t0_frc = MPI_Wtime()
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("update_atm_forcing")
#endif
#if defined (__yac)
           if (f%mype==0)  print *, achar(27)//'[34m'//' --> call update_atm_forcing(n)'//achar(27)//'[0m'
            ! call-gating: only exchange on coupling steps (cpl_stride). Off-steps
            ! reuse the cached ocean/atm fields (zero-order hold). Ocean side gates
            ! identically (same n, same cpl_stride) so yac_fput/fget stay paired.
            if (mod(n-1, cpl_stride) == 0) &
            call update_atm_forcing_yac(n, f%ice, f%tracers, f%dynamics, f%partit, f%mesh)
           if (f%mype==0)  print *, achar(27)//'[34m'//' --> after update_atm_forcing(n)'//achar(27)//'[0m'
#endif 
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("update_ice_oce_atm_forcing")
#endif
            f%t1_frc = MPI_Wtime()
            !___compute ice step________________________________________________
            if (f%ice%ice_steps_since_upd>=f%ice%ice_ave_steps-1) then
                f%ice%ice_update=.true.
                f%ice%ice_steps_since_upd = 0
            else
                f%ice%ice_update=.false.
                f%ice%ice_steps_since_upd=f%ice%ice_steps_since_upd+1
            endif
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_timestep(n)'//achar(27)//'[0m'
            if (f%ice%ice_update) then
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("ice_timestep")
#endif
        if (f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_timestep(n)'//achar(27)//'[0m'
                call ice_timestep(n, f%ice, f%partit, f%mesh)
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("ice_timestep")
#endif
            endif
            !SLThe following will be calculated from the ocean side
            !___compute fluxes to the ocean: heat, freshwater, momentum_________
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call oce_fluxes_mom...'//achar(27)//'[0m'
            call oce_fluxes_mom(f%ice, f%dynamics, f%partit, f%mesh) ! momentum only: fills ice%stress_iceoce_x/y for the ice->ocean YAC send
            ! (no oce_fluxes: the heat/freshwater flux to the ocean is sent as
            !  `ice_to_ocean_flux` from gen_forcing_couple, not applied locally)
!sl        end if  !sl??
        f%t2 = MPI_Wtime()

        
        !___model ice step____________________________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_timestep_ale'//achar(27)//'[0m'
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("ice_timestep_ale")
#endif
!sl!!        call ice_timestep_ale(n, f%ice, f%partit, f%mesh)
        call ice_timestep_ale(n, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("ice_timestep_ale")
#endif

        f%t3 = MPI_Wtime()
        !___compute energy diagnostics..._______________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call compute_diagnostics(1)'//achar(27)//'[0m'
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("compute_diagnostics_fesim")
#endif
        call compute_diagnostics(1, f%dynamics, f%tracers, f%ice, f%partit, f%mesh)
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("compute_diagnostics_fesim")
#endif

        f%t4 = MPI_Wtime()
        !___prepare output______________________________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call output (n)'//achar(27)//'[0m'
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("output")
#endif
        call output(n, f%ice,f%dynamics,f%tracers,f%partit, f%mesh)
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("output")
#endif

        f%t5 = MPI_Wtime()
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("restart")
#endif
!sl        call write_initial_conditions(n, nstart, f%total_nsteps, f%which_readr, f%ice, f%partit, f%mesh)
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("restart")
#endif
        f%t6 = MPI_Wtime()
        
        f%rtime_fullice       = f%rtime_fullice       + f%t3 - f%t1   ! fold ice_timestep_ale (t2->t3) into ice compute [scalability instr. 2026-07]
        f%rtime_compute_diag  = f%rtime_compute_diag  + f%t4 - f%t3
        f%rtime_write_means   = f%rtime_write_means   + f%t5 - f%t4
        f%rtime_write_restart = f%rtime_write_restart + f%t6 - f%t5
        f%rtime_read_forcing  = f%rtime_read_forcing  + f%t1_frc - f%t0_frc

    end do
!call cray_acc_set_debug_global_level(3)    
    f%from_nstep = f%from_nstep+current_nsteps
!call cray_acc_set_debug_global_level(0)    
!   write(0,*) 'f%from_nstep after the loop:', f%from_nstep    
    
    ! End main time loop profiling
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("fesim_runloop_total")
#endif
  end subroutine fesim_runloop


  subroutine fesim_finalize()
    use fesom_main_storage_module
#if defined(__MULTIO)
    use iom
    use mpp_io
#endif
    ! EO parameters
    real(kind=real32) :: mean_rtime(15), max_rtime(15), min_rtime(15)
    integer           :: tr_num
    
    ! Start finalization profiling
#if defined (FESIM_PROFILING)
        call fesim_profiler_start("fesim_finalize_total")
#endif

    call finalize_output()
    call finalize_restart()

    !___FINISH MODEL RUN________________________________________________________

    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)
    !$ACC EXIT DATA DELETE (f%ice%delta_min, f%ice%Tevp_inv, f%ice%cd_oce_ice)
    !$ACC EXIT DATA DELETE (f%ice%work%fct_tmax, f%ice%work%fct_tmin)
    !$ACC EXIT DATA DELETE (f%ice%work%fct_fluxes, f%ice%work%fct_plus, f%ice%work%fct_minus)
    !$ACC EXIT DATA DELETE (f%ice%work%eps11, f%ice%work%eps12, f%ice%work%eps22)
    !$ACC EXIT DATA DELETE (f%ice%work%sigma11, f%ice%work%sigma12, f%ice%work%sigma22)
    !$ACC EXIT DATA DELETE (f%ice%work%ice_strength, f%ice%stress_atmice_x, f%ice%stress_atmice_y)
    !$ACC EXIT DATA DELETE (f%ice%thermo%rhosno, f%ice%thermo%rhoice, f%ice%thermo%inv_rhowat)
    !$ACC EXIT DATA DELETE (f%ice%srfoce_ssh, f%ice%pstar, f%ice%c_pressure)
    !$ACC EXIT DATA DELETE (f%ice%work%inv_areamass, f%ice%work%inv_mass, f%ice%uice_rhs, f%ice%vice_rhs)
    !$ACC EXIT DATA DELETE (f%ice%uice, f%ice%vice, f%ice%srfoce_u, f%ice%srfoce_v, f%ice%uice_old, f%ice%vice_old)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%values, f%ice%data(2)%values, f%ice%data(3)%values)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%valuesl, f%ice%data(2)%valuesl, f%ice%data(3)%valuesl)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%dvalues, f%ice%data(2)%dvalues, f%ice%data(3)%dvalues)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%values_rhs, f%ice%data(2)%values_rhs, f%ice%data(3)%values_rhs)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%values_div_rhs, f%ice%data(2)%values_div_rhs, f%ice%data(3)%values_div_rhs)
#if defined (__oifs) || defined (__ifsinterface)
    !$ACC EXIT DATA DELETE (f%ice%data(4)%values, f%ice%data(4)%valuesl, f%ice%data(4)%dvalues, f%ice%data(4)%values_rhs, f%ice%data(4)%values_div_rhs)
#endif
    !$ACC EXIT DATA DELETE (f%ice%data, f%ice%work, f%ice%work%fct_massmatrix)
    !$ACC EXIT DATA DELETE (f%ice)
    !$ACC EXIT DATA DELETE (f%tracers%work%fct_ttf_min, f%tracers%work%fct_ttf_max, f%tracers%work%fct_plus, f%tracers%work%fct_minus)
    !$ACC EXIT DATA DELETE (f%tracers%work%adv_flux_hor, f%tracers%work%adv_flux_ver, f%tracers%work%fct_LO)
    !$ACC EXIT DATA DELETE (f%tracers%work%del_ttf_advvert, f%tracers%work%del_ttf_advhoriz, f%tracers%work%edge_up_dn_grad)
    !$ACC EXIT DATA DELETE (f%tracers%work%del_ttf)
    !$ACC EXIT DATA DELETE (tr_xy, tr_z, relax2clim, Sclim, Tclim)
    !$ACC EXIT DATA DELETE (f%tracers%data, f%tracers%work)
    !$ACC EXIT DATA DELETE (f%dynamics%w, f%dynamics%w_e, f%dynamics%uv)
    !$ACC EXIT DATA DELETE (f%dynamics, f%tracers)

    !delete mesh and partit data.
    !$ACC EXIT DATA DELETE (f%mesh%coriolis_node, f%mesh%nn_num, f%mesh%nn_pos) 
    !$ACC EXIT DATA DELETE (f%mesh%ssh_stiff, f%mesh%ssh_stiff%rowptr) 
    !$ACC EXIT DATA DELETE (f%mesh%gradient_sca, f%mesh%metric_factor, f%mesh%elem_area, f%mesh%area, f%mesh%edge2D_in) 
    !$ACC EXIT DATA DELETE (f%mesh%elem2D_nodes, f%mesh%ulevels, f%mesh%ulevels_nod2d, f%mesh%edges, f%mesh%edge_tri) 
    !$ACC EXIT DATA DELETE (f%mesh%helem, f%mesh%elem_cos, f%mesh%edge_cross_dxdy, f%mesh%elem2d_nodes, f%mesh%nl) 
    !$ACC EXIT DATA DELETE (f%mesh%nlevels_nod2D, f%mesh%nod_in_elem2D, f%mesh%nod_in_elem2D_num) 
    !$ACC EXIT DATA DELETE (f%mesh%edge_dxdy, f%mesh%nlevels, f%mesh%hnode, f%mesh%hnode_new, f%mesh%ulevels_nod2D_max) 
    !$ACC EXIT DATA DELETE (f%mesh%zbar_3d_n, f%mesh%z_3d_n, f%mesh%areasvol, f%mesh%nlevels_nod2D_min) 
    !$ACC EXIT DATA DELETE (f%partit%eDim_nod2D, f%partit%myDim_edge2D) 
    !$ACC EXIT DATA DELETE (f%partit%myDim_elem2D, f%partit%myDim_nod2D, f%partit%myList_edge2D) 
    !$ACC EXIT DATA DELETE (f%mesh, f%partit, f)
    if (f%mype==0) then
       f%t1 = MPI_Wtime()
       f%runtime_alltimesteps = real(f%t1-f%t0,real32)
       write(*,*) 'FESIM Run is finished, updating clock'
    endif

!    mean_rtime(1)  = rtime_oce         
!    mean_rtime(2)  = rtime_oce_mixpres 
!    mean_rtime(3)  = rtime_oce_dyn     
!    mean_rtime(4)  = rtime_oce_dynssh  
!    mean_rtime(5)  = rtime_oce_solvessh
!    mean_rtime(6)  = rtime_oce_GMRedi  
!    mean_rtime(7)  = rtime_oce_solvetra
    mean_rtime(8)  = rtime_ice         
!    mean_rtime(9)  = rtime_tot  
    mean_rtime(10) = f%rtime_fullice - f%rtime_read_forcing 
    mean_rtime(11) = f%rtime_compute_diag
    mean_rtime(12) = f%rtime_write_means
    mean_rtime(13) = f%rtime_write_restart
    mean_rtime(14) = f%rtime_read_forcing
    ! FESIM scalability instrumentation (2026-07): fill the indices the summary
    ! actually prints. (1) = pure ice compute (EVP dyn + thermo + ale advection,
    ! excluding the yac recv wait); (14) = f%rtime_read_forcing already brackets
    ! update_atm_forcing_yac, i.e. the yac_fget block = time ice waits on the
    ! ocean; (9) = FESIM per-task loop total (compute + coupling wait + diag/io).
    mean_rtime(1)  = f%rtime_fullice - f%rtime_read_forcing
    mean_rtime(9)  = f%rtime_fullice + f%rtime_compute_diag &
                   + f%rtime_write_means + f%rtime_write_restart
    max_rtime(1:14) = mean_rtime(1:14)
    min_rtime(1:14) = mean_rtime(1:14)

    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 14, MPI_REAL, MPI_SUM, f%MPI_COMM_FESOM, f%MPIerr)
!sl    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 14, MPI_REAL, MPI_SUM, f%MPI_COMM_FESIM, f%MPIerr)
    mean_rtime(1:14) = mean_rtime(1:14) / real(f%npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  14, MPI_REAL, MPI_MAX, f%MPI_COMM_FESOM, f%MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  14, MPI_REAL, MPI_MIN, f%MPI_COMM_FESOM, f%MPIerr)
!sl    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  14, MPI_REAL, MPI_MAX, f%MPI_COMM_FESIM, f%MPIerr)
!sl    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  14, MPI_REAL, MPI_MIN, f%MPI_COMM_FESIM, f%MPIerr)

#if defined (__yac)
    ! Coupler cost, attributed (2026-09-11, ANALYSIS.md §10.2). Collective, so it
    ! must run before par_ex below finalizes MPI. "runtime yac recv/wait" above is
    ! NOT pure wait -- it also holds these calls plus 14 halo exchanges and the
    ! unit conversions; this block separates out the YAC part.
    call cpl_timers_report(f%MPI_COMM_FESOM, f%mype, f%npes, 'fesim')
#endif
    
!sl#if defined (__oifs) 
!sl    ! OpenIFS coupled version has to call oasis_terminate through par_ex
!sl    call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype)
!sl#endif

#if defined(__MULTIO) && !defined(__ifsinterface) && !defined(__oasis)
   call mpp_stop
#endif
    ! Generate enhanced profiler report BEFORE MPI finalization
!sl#if defined (FESOM_PROFILING)
!sl        call fesom_profiler_end("fesom_finalize_total")
!sl        call fesom_profiler_report(f%MPI_COMM_FESOM, f%mype)
!sl#endif
#if defined (FESIM_PROFILING)
        call fesim_profiler_end("fesim_finalize_total")
        call fesim_profiler_report(f%MPI_COMM_FESOM, f%mype)
        ! Note: Do NOT call fesim_profiler_finalize here as it would duplicate the report
#endif
    
    if(f%fesim_did_mpi_init) call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype) ! finalize MPI before FESOM prints its stats block, otherwise there is sometimes output from other processes from an earlier time in the programm AFTER the starts block (with parastationMPI)
!sl    if(f%fesim_did_mpi_init) call par_ex(f%partit%MPI_COMM_FESIM, f%partit%mype) ! finalize MPI before FESOM prints its stats block, otherwise there is sometimes output from other processes from an earlier time in the programm AFTER the starts block (with parastationMPI)
    if (f%mype==0) then
        41 format (a35,a10,2a15) !Format for table heading
        42 format (a30,3f15.4)   !Format for table content

        print 41, '___MODEL RUNTIME per task [seconds]','_____mean_','___________min_', '___________max_'
        print 42, '  runtime ice compute       :',    mean_rtime(1),     min_rtime(1),      max_rtime(1)
        print 42, '  runtime yac recv/wait     :',    mean_rtime(14),    min_rtime(14),     max_rtime(14)
        print 42, '  runtime diag              :',    mean_rtime(11),    min_rtime(11),     max_rtime(11)
        print 42, '  runtime output            :',    mean_rtime(12),    min_rtime(12),     max_rtime(12)
        print 42, '  runtime restart           :',    mean_rtime(13),    min_rtime(13),     max_rtime(13)
        print 42, '  runtime total (fesim)     :',    mean_rtime(9),     min_rtime(9),      max_rtime(9)
!sl#if defined (__recom)
!sl        print 42, '  runtime recom:              ',    mean_rtime(15),    min_rtime(15),     max_rtime(15)
!sl#endif

        43 format (a33,i15)        !Format Ncores
        44 format (a33,i15)        !Format OMP threads
        45 format (a33,f15.4,a4)   !Format runtime

        write(*,*)
        write(*,*) '======================================================'
        write(*,*) '================ BENCHMARK RUNTIME ==================='
        print 43, '    Number of cores :            ',f%npes
#if defined(_OPENMP)
        print 44, '    Max OpenMP threads :         ',OMP_GET_MAX_THREADS()
#endif
        print 45, '    Runtime for all timesteps :  ',f%runtime_alltimesteps,' sec'
        write(*,*) '======================================================'
        write(*,*)
    end if    
!   call clock_finish  
    
    ! Enhanced profiler is already finalized above before MPI finalization
  end subroutine fesim_finalize

end module fesim_module

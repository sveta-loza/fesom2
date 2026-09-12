! synopsis: FESIM -- the sea-ice component of the FESOM/FESIM split.
!
! Built from the FESOM source tree with this driver in place of fesom_module.
! Same three-part structure (init / runloop / finalize). The ocean surface
! state and the atmospheric forcing arrive from the ocean over YAC
! (ocean_coupling_interface, update_atm_forcing_yac); the ice state, the
! ice-ocean drag and the ice-ocean heat/freshwater flux go back the same way.
! What the ocean forwards from its atmosphere is chosen at run time in
! namelist.cpl (cpl_config).
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
  use dynamics_init_interface
  use init_ale_interface
  use ice_setup_interface
  use ice_timestep_interface
  use oce_fluxes_interface
  use update_atm_forcing_interface
  use read_mesh_interface
  use fesom_version_info_module
  use command_line_options_module
  use, intrinsic :: iso_fortran_env, only : real32
  use g_forcing_param, only: use_landice_water, use_age_tracer
  use iceberg_params
  use iceberg_ocean_coupling

#if defined (__icepack)
  use icedrv_main,          only: set_icepack, init_icepack, alloc_icepack
#endif
#if defined (__cpl_yac)
  use cpl_config, only: read_cpl_namelist, check_cpl_config, cpl_component_is_sea_ice, &
                        cpl_has_atmosphere
  use ocean_coupling_interface
#endif

  implicit none

  type :: fesom_main_storage_type

    integer           :: n, from_nstep, offset, row, i, provided, id
    integer           :: which_readr ! read which restart files (0=netcdf, 1=core dump,2=dtype)
    integer           :: total_nsteps
    integer, pointer  :: mype, npes, MPIerr, MPI_COMM_FESOM, MPI_COMM_WORLD, MPI_COMM_FESOM_IB
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

    character(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version_txt
    integer mpi_version_len
    logical fesim_did_mpi_init

  end type fesom_main_storage_type
  type(fesom_main_storage_type), save, target :: f

end module fesom_main_storage_module


module fesim_module
  use fesom_main_storage_module
#if defined (FESOM_PROFILING)
  use fesom_profiler
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

      if(command_argument_count() > 0) then
        call command_line_options%parse()
        stop
      end if

      mpi_is_initialized = .false.
      f%fesim_did_mpi_init = .false.

#if defined (__cpl_yac)
      ! Before the coupler is initialised: yac_fdef_comp needs the component
      ! name from &coupling_yac, and the partner flags describe what the
      ! ocean forwards to this component.
      cpl_component_is_sea_ice = .true.
      call read_cpl_namelist()
#endif

      call MPI_Initialized(mpi_is_initialized, f%i)
      if(.not. mpi_is_initialized) then
          call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, f%provided, f%i)
          f%fesim_did_mpi_init = .true.
      end if

#if defined (__cpl_yac)
      call ocn_cpl_init(f%partit%MPI_COMM_FESOM)
#endif

      f%t1 = MPI_Wtime()

#if defined (FESOM_PROFILING)
      call fesom_profiler_init(.true.)
      call fesom_profiler_start("fesom_init_total")
      call fesom_profiler_start("par_init")
#endif
      call par_init(f%partit)
#if defined (FESOM_PROFILING)
      call fesom_profiler_end("par_init")
#endif

      f%mype              => f%partit%mype
      f%MPIerr            => f%partit%MPIerr
      f%MPI_COMM_FESOM    => f%partit%MPI_COMM_FESOM
      f%MPI_COMM_FESOM_IB => f%partit%MPI_COMM_FESOM_IB
      f%MPI_COMM_WORLD    => f%partit%MPI_COMM_WORLD
      f%npes              => f%partit%npes

      if(f%mype==0) then
          write(*,*)
          print *,"FESIM (FESOM2 source) git SHA: "//fesom_git_sha()
          call MPI_Get_library_version(f%mpi_version_txt, f%mpi_version_len, f%MPIERR)
          print *,"MPI library version: "//trim(f%mpi_version_txt)
          print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
          print *, achar(27)//'[7;32m'//' --> FESIM BUILDS UP MODEL CONFIGURATION                    '//achar(27)//'[0m'
      end if
      !=====================
      ! Read configuration data, load the mesh and fill in auxiliary mesh arrays
      !=====================
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call setup_model'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
      call fesom_profiler_start("setup_model")
#endif
      call setup_model(f%partit)  ! Read Namelists, always before clock_init
#if defined (__cpl_yac)
      call check_cpl_config(.false., f%partit%MPI_COMM_FESOM, f%mype)
#endif
      ! The sea-ice component restarts from the portable (netCDF) ice restart
      ! only. The raw and binary restarts dump the ocean's dynamics and tracer
      ! containers as well, which hold nothing here, so they are switched off
      ! whatever namelist.config says.
      if (raw_restart_length_unit /= 'off' .or. bin_restart_length_unit /= 'off') then
         if (f%mype==0) write(*,*) 'FESIM: raw and binary restarts are not used by the sea-ice', &
                                   ' component; only the netCDF ice restart is written'
         raw_restart_length_unit = 'off'
         bin_restart_length_unit = 'off'
      end if
#if defined (FESOM_PROFILING)
      call fesom_profiler_end("setup_model")
#endif

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call clock_init'//achar(27)//'[0m'
      call clock_init(f%partit)   ! read the clock file

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call get_run_steps'//achar(27)//'[0m'
      call get_run_steps(fesim_total_nsteps, f%partit)
      f%total_nsteps=fesim_total_nsteps
#if defined (FESOM_PROFILING)
      call fesom_profiler_set_timesteps(fesim_total_nsteps)
      call fesom_profiler_set_timestep_size(86400.0d0 / real(step_per_day, kind=8))
#endif

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call mesh_setup'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
      call fesom_profiler_start("mesh_setup")
#endif
      call mesh_setup(f%partit, f%mesh)
#if defined (FESOM_PROFILING)
      call fesom_profiler_end("mesh_setup")
#endif
      if (f%mype==0) write(*,*) 'FESIM mesh_setup... complete'

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call check_mesh_consistency'//achar(27)//'[0m'
      call check_mesh_consistency(f%partit, f%mesh)
      if (f%mype==0) f%t2=MPI_Wtime()

      !=====================
      ! The parts of ocean_setup the sea ice needs: the dynamics/tracer/array
      ! containers the ice and the I/O address, the ALE layer arrays, the SSH
      ! stiffness matrix (its neighbourhood lists size the FCT mass matrix)
      ! and the MUSCL neighbourhood used by the ice advection. No ocean
      ! initial state is read: the ocean surface arrives over YAC.
      !=====================
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call dynamics_init'//achar(27)//'[0m'
      call dynamics_init(f%dynamics, f%partit, f%mesh)
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call tracer_init'//achar(27)//'[0m'
      call tracer_init(f%tracers, f%partit, f%mesh)
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call arrays_init'//achar(27)//'[0m'
      call arrays_init(f%tracers%num_tracers, f%partit, f%mesh)
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call init_ale'//achar(27)//'[0m'
      call init_ale(f%dynamics, f%partit, f%mesh)
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call init_stiff_mat_ale'//achar(27)//'[0m'
      call init_stiff_mat_ale(f%partit, f%mesh)
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call muscl_adv_init'//achar(27)//'[0m'
      call muscl_adv_init(f%tracers%work, f%partit, f%mesh)

      call forcing_setup(f%partit, f%mesh)

      if (f%mype==0) f%t4=MPI_Wtime()
      if (use_ice) then
          if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_setup'//achar(27)//'[0m'
          call ice_setup(f%ice, f%tracers, f%partit, f%mesh)
          if (f%mype==0) write(*,*) 'EVP scheme option=', f%ice%whichEVP
      else
          call ice_init_toyocean_dummy(f%ice, f%partit, f%mesh)
      endif
      if (f%mype==0) f%t5=MPI_Wtime()

      call compute_diagnostics(0, f%dynamics, f%tracers, f%ice, f%partit, f%mesh) ! allocate arrays for diagnostic

#if defined (__cpl_yac)
      ! The ocean<->ice fields are exchanged every cpl_stride steps.
      call ocn_cpl_define(f%partit, f%mesh, dt*cpl_stride, f%total_nsteps)
      if(f%mype==0)  write(*,*) 'FESIM ---->     YAC fields defined, nsend/nrecv:', OCN_NSEND, OCN_NRECV
#endif

#if defined (__icepack)
      if (f%mype==0) write(*,*) 'Icepack: reading namelists from namelist.icepack'
      call set_icepack(f%ice, f%partit)
      call alloc_icepack
      call init_icepack(f%ice, f%tracers%data(1), f%mesh)
      if (f%mype==0) write(*,*) 'Icepack: setup complete'
#endif
      call clock_newyear                        ! check if it is a new year
      if (f%mype==0) f%t6=MPI_Wtime()
      !___READ INITIAL CONDITIONS IF THIS IS A RESTART RUN______________________
      if (r_restart) then
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

#if defined (FESOM_PROFILING)
      call fesom_profiler_end("fesom_init_total")
#endif
  end subroutine fesim_init


  subroutine fesim_runloop(current_nsteps)
    use fesom_main_storage_module
    integer, intent(in) :: current_nsteps
    ! EO parameters
    integer n, nstart, ntotal

    f%MPI_COMM_FESOM_IB = f%MPI_COMM_FESOM
    if (f%mype==0) then
        write (*,*) 'current_nsteps, steps_per_ib_step, icb_outfreq :', current_nsteps, steps_per_ib_step, icb_outfreq
    end if

    if (f%mype==0) write(*,*) 'FESIM start iteration before the barrier...'
#if defined (__cpl_yac)
    ! Cross-component sync before the loop timer starts, so that ocean and
    ! ice resume from the same instant and the ocean's longer init does not
    ! show up as ice-side wait in the first step. MPI_Barrier is collective
    ! over the whole communicator, so this is only possible when
    ! MPI_COMM_WORLD is exactly ocean + ice: with an atmosphere in the same
    ! MPMD world its ranks never call it and the barrier would hang.
    if (.not. cpl_has_atmosphere()) call MPI_Barrier(MPI_COMM_WORLD, f%MPIERR)
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

#if defined (FESOM_PROFILING)
    call fesom_profiler_start("fesom_runloop_total")
#endif
    !___MODEL TIME STEPPING LOOP________________________________________________
    nstart=f%from_nstep
    ntotal=f%from_nstep-1+current_nsteps

    do n=nstart, ntotal
        mstep = n
        if (mod(n,logfile_outfreq)==0 .and. f%mype==0) then
            write(*,*) 'FESIM ======================================================='
            write(*,*) 'FESIM step:',n,' day:', daynew,' year:',yearnew
            write(*,*)
        end if
        call clock

        !___exchange with the ocean: send the ice state, receive the ocean
        !   surface state and the atmospheric forcing (every cpl_stride steps;
        !   the ocean gates identically so that yac_fput/yac_fget stay paired)
        f%t1 = MPI_Wtime()
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call update_atm_forcing_yac(n)'//achar(27)//'[0m'
        f%t0_frc = MPI_Wtime()
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("update_atm_forcing")
#endif
#if defined (__cpl_yac)
        if (mod(n-1, cpl_stride) == 0) &
        call update_atm_forcing_yac(n, f%ice, f%tracers, f%dynamics, f%partit, f%mesh)
#endif
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("update_atm_forcing")
#endif
        f%t1_frc = MPI_Wtime()

        !___sea-ice step________________________________________________________
        if (f%ice%ice_steps_since_upd>=f%ice%ice_ave_steps-1) then
            f%ice%ice_update=.true.
            f%ice%ice_steps_since_upd = 0
        else
            f%ice%ice_update=.false.
            f%ice%ice_steps_since_upd=f%ice%ice_steps_since_upd+1
        endif
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_timestep(n)'//achar(27)//'[0m'
        if (f%ice%ice_update) then
#if defined (FESOM_PROFILING)
            call fesom_profiler_start("ice_timestep")
#endif
            call ice_timestep(n, f%ice, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
            call fesom_profiler_end("ice_timestep")
#endif
        endif

        !___ice-ocean momentum drag, for the send to the ocean__________________
        ! (the heat and freshwater flux is sent from the thermodynamics arrays;
        !  nothing is applied to a local ocean here)
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call oce_fluxes_mom'//achar(27)//'[0m'
        call oce_fluxes_mom(f%ice, f%dynamics, f%partit, f%mesh)
        f%t2 = MPI_Wtime()

        ! collective blow-up status (the ocean step's check, kept for the ice)
        call status_check(f%partit)
        f%t3 = MPI_Wtime()

        !___diagnostics_________________________________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call compute_diagnostics(1)'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("compute_diagnostics")
#endif
        call compute_diagnostics(1, f%dynamics, f%tracers, f%ice, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("compute_diagnostics")
#endif
        f%t4 = MPI_Wtime()

        !___output______________________________________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call output (n)'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("output")
#endif
        call output(n, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("output")
#endif
        f%t5 = MPI_Wtime()
        !___restart (netCDF ice restart + clock; no ocean group is registered)__
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call write_initial_conditions(n,...)'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("restart")
#endif
        call write_initial_conditions(n, nstart, f%total_nsteps, f%which_readr, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("restart")
#endif
        f%t6 = MPI_Wtime()

        f%rtime_fullice       = f%rtime_fullice       + f%t3 - f%t1
        f%rtime_compute_diag  = f%rtime_compute_diag  + f%t4 - f%t3
        f%rtime_write_means   = f%rtime_write_means   + f%t5 - f%t4
        f%rtime_write_restart = f%rtime_write_restart + f%t6 - f%t5
        f%rtime_read_forcing  = f%rtime_read_forcing  + f%t1_frc - f%t0_frc
    end do
    f%from_nstep = f%from_nstep+current_nsteps

#if defined (FESOM_PROFILING)
    call fesom_profiler_end("fesom_runloop_total")
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

#if defined (FESOM_PROFILING)
    call fesom_profiler_start("fesom_finalize_total")
#endif

    call finalize_output()
    call finalize_restart()

    !___FINISH MODEL RUN________________________________________________________
    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)
    if (f%mype==0) then
       f%t1 = MPI_Wtime()
       f%runtime_alltimesteps = real(f%t1-f%t0,real32)
       write(*,*) 'FESIM Run is finished, updating clock'
    endif

    mean_rtime = 0.0_real32
    ! (1) pure ice compute (EVP dyn + thermo + advection, excluding the yac
    ! receive wait); (14) = the exchange with the ocean, i.e. the time the ice
    ! waits on the ocean; (9) = per-task loop total.
    mean_rtime(8)  = rtime_ice
    mean_rtime(10) = f%rtime_fullice - f%rtime_read_forcing
    mean_rtime(11) = f%rtime_compute_diag
    mean_rtime(12) = f%rtime_write_means
    mean_rtime(13) = f%rtime_write_restart
    mean_rtime(14) = f%rtime_read_forcing
    mean_rtime(1)  = f%rtime_fullice - f%rtime_read_forcing
    mean_rtime(9)  = f%rtime_fullice + f%rtime_compute_diag &
                   + f%rtime_write_means + f%rtime_write_restart
    max_rtime(1:14) = mean_rtime(1:14)
    min_rtime(1:14) = mean_rtime(1:14)

    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 14, MPI_REAL, MPI_SUM, f%MPI_COMM_FESOM, f%MPIerr)
    mean_rtime(1:14) = mean_rtime(1:14) / real(f%npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  14, MPI_REAL, MPI_MAX, f%MPI_COMM_FESOM, f%MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  14, MPI_REAL, MPI_MIN, f%MPI_COMM_FESOM, f%MPIerr)

#if defined (__cpl_yac)
    ! Coupler cost, collective: must run before par_ex finalizes MPI. The
    ! "yac recv/wait" line below is not pure wait -- it also holds the halo
    ! exchanges and unit conversions; this block separates out the YAC part.
    call cpl_timers_report(f%MPI_COMM_FESOM, f%mype, f%npes, 'fesim')
#endif

#if defined(__MULTIO) && !defined (__cpl_direct) && !defined (__cpl_oasis)
   call mpp_stop
#endif
#if defined (FESOM_PROFILING)
    call fesom_profiler_end("fesom_finalize_total")
    call fesom_profiler_report(f%MPI_COMM_FESOM, f%mype)
#endif

    if(f%fesim_did_mpi_init) call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype) ! finalize MPI before the stats block

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
  end subroutine fesim_finalize

end module fesim_module

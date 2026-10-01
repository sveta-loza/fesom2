!=====================================================================
! ifs_mpmd_world -- MPI_COMM_WORLD choreography for the MPMD layout
!   "IFS executable (ifsMASTER, FESOM linked in as a library)" + "FESIM executable"
! (Stage 4 config #1, __ifs_fwd). Identical copy in the FESOM and FESIM trees.
!
! Background. In the IFS-FESOM single executable FESOM runs ON the IFS compute
! tasks. The only MPI_COMM_WORLD-collective operations that executable performs
! are, in this order (verified against DE_CY48R1 ClimateDT, 2026-10-01):
!   (1) mpp_io_init (ifs_interface/mpp_io.F90, FESOM code; first thing done by
!       arpifs/programs/master.F90 through ININEMOIO when NEMOIOSERVER=yes):
!       MPI_Comm_split(MPI_COMM_WORLD) separating the compute tasks from the
!       FESOM multio-server tasks. The compute part becomes IFS's MPLUSERCOMM
!       and IFS/MPL run on that sub-communicator from then on.
!   (2) master.F90: MPI_Barrier(MPI_COMM_WORLD) iff the environment variable
!       MPI_EPOCH is non-blank (raps always exports it).
!   (3) master.F90: one MPI_Comm_split(MPI_COMM_WORLD) (compute vs IFS IO servers).
! multio works on sub-communicators. A separate FESIM executable that shares
! MPI_COMM_WORLD with ifsMASTER (srun heterogeneous block / --multi-prog) must
! therefore mirror exactly these three calls, in this order, before anything
! else: that is fesim_mpmd_join_ifs_world.
!
! The YAC world (FESOM compute tasks + FESIM tasks) can NOT be built by yet
! another split of MPI_COMM_WORLD: the only point where all ranks are still in
! lockstep is (1), and there the IFS IO servers are not yet separated from the
! compute tasks (that happens in (3)), so they would end up inside the YAC
! communicator and never enter YAC (found the hard way, 2026-10-01). Instead the
! two groups are joined once both know their own communicator, with
! MPI_Intercomm_create + MPI_Intercomm_merge: these are collective over the two
! local communicators only, plus leader-to-leader messages on MPI_COMM_WORLD.
! The IFS IO servers and the multio servers are not involved.
!
! Rank layout in MPI_COMM_WORLD, fixed by the launcher:
!   [ IFS compute + IFS IO servers ][ FESOM multio servers ][ FESIM ]
! FESIM occupies the LAST ranks. The IFS executable learns their number from
! the environment variable FESIM_NTASKS (default 0 = no FESIM: every code path
! below is then skipped and the behaviour is identical to upstream FESOM).
!=====================================================================
module ifs_mpmd_world
  use mpi
  implicit none
  private

  ! split (1) uses colours 1 (compute) and 3 (server); split (3) uses 1 and 2.
  integer, parameter :: COLOR_FESIM_IN_IOSPLIT = 7
  integer, parameter :: COLOR_FESIM_IN_MASTER  = 9
  integer, parameter :: TAG_LEADER_RANK = 7731
  integer, parameter :: TAG_INTERCOMM   = 7732

  ! Communicator spanning FESOM compute tasks + FESIM tasks (FESOM ranks first);
  ! MPI_COMM_NULL on ranks that are not part of it and whenever no FESIM is present.
  integer, save, public :: yac_world_comm = MPI_COMM_NULL

  public :: ifs_mpmd_fesim_ntasks, ifs_mpmd_have_yac_world
  public :: ifs_mpmd_connect_yac_world, fesim_mpmd_join_ifs_world

contains

  ! Number of FESIM tasks appended to MPI_COMM_WORLD (env FESIM_NTASKS, default 0).
  integer function ifs_mpmd_fesim_ntasks()
    character(len=32) :: buf
    integer :: l, stat, n
    ifs_mpmd_fesim_ntasks = 0
    call get_environment_variable('FESIM_NTASKS', buf, l, stat)
    if (stat /= 0 .or. l <= 0) return
    read(buf(1:l), *, iostat=stat) n
    if (stat == 0 .and. n > 0) ifs_mpmd_fesim_ntasks = n
  end function ifs_mpmd_fesim_ntasks

  logical function ifs_mpmd_have_yac_world()
    ifs_mpmd_have_yac_world = (yac_world_comm /= MPI_COMM_NULL)
  end function ifs_mpmd_have_yac_world

  ! Join the two groups into yac_world_comm. Collective over local_comm on each
  ! side: the FESOM side calls it with its compute communicator (IFS's icomm,
  ! after split (3)), the FESIM side with its own communicator. The FESOM
  ! leader tells the FESIM leader (world rank world_size - FESIM_NTASKS) its
  ! world rank, then both groups create and merge an inter-communicator.
  subroutine ifs_mpmd_connect_yac_world(local_comm, is_fesim)
    integer, intent(in) :: local_comm
    logical, intent(in) :: is_fesim
    integer :: ierr, world_rank, world_size, local_rank, remote_leader, inter, ntask_fesim
    integer :: status(MPI_STATUS_SIZE)

    call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)
    call MPI_Comm_rank(local_comm, local_rank, ierr)

    if (is_fesim) then
       if (local_rank == 0) call MPI_Recv(remote_leader, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                                          TAG_LEADER_RANK, MPI_COMM_WORLD, status, ierr)
       call MPI_Bcast(remote_leader, 1, MPI_INTEGER, 0, local_comm, ierr)
    else
       ntask_fesim = ifs_mpmd_fesim_ntasks()
       if (ntask_fesim <= 0) then
          yac_world_comm = MPI_COMM_NULL
          return
       end if
       remote_leader = world_size - ntask_fesim
       if (local_rank == 0) call MPI_Send(world_rank, 1, MPI_INTEGER, remote_leader, &
                                          TAG_LEADER_RANK, MPI_COMM_WORLD, ierr)
    end if

    call MPI_Intercomm_create(local_comm, 0, MPI_COMM_WORLD, remote_leader, TAG_INTERCOMM, inter, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(*,*) 'ifs_mpmd_connect_yac_world: MPI_Intercomm_create failed, ierr=', ierr
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    ! FESOM ranks come first in the merged communicator, FESIM ranks last.
    call MPI_Intercomm_merge(inter, is_fesim, yac_world_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(*,*) 'ifs_mpmd_connect_yac_world: MPI_Intercomm_merge failed, ierr=', ierr
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    call MPI_Comm_free(inter, ierr)
    if (local_rank == 0) then
       call MPI_Comm_size(yac_world_comm, world_size, ierr)
       write(*,*) 'ifs_mpmd_connect_yac_world: YAC world communicator built, size', world_size, &
                  ' (fesim side: ', is_fesim, ')'
    end if
  end subroutine ifs_mpmd_connect_yac_world

  ! FESIM-executable side: mirror the three world-collective calls (1)-(3) of
  ! ifsMASTER, then connect to the FESOM compute tasks. Call right after
  ! MPI_Init and before any YAC call; afterwards yac_world_comm is the
  ! communicator to hand to yac_finit_comm.
  subroutine fesim_mpmd_join_ifs_world()
    integer :: ierr, key, dummy, fesim_comm
    character(len=64) :: epoch
    call MPI_Comm_rank(MPI_COMM_WORLD, key, ierr)
    ! (1) mpp_io_init's compute/server split -- FESIM is neither
    call MPI_Comm_split(MPI_COMM_WORLD, COLOR_FESIM_IN_IOSPLIT, key, dummy, ierr)
    call MPI_Comm_free(dummy, ierr)
    ! (2) master.F90's MPI-startup-cost barrier, taken iff MPI_EPOCH is non-blank
    epoch = ' '
    call get_environment_variable('MPI_EPOCH', epoch)
    if (epoch /= ' ') call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! (3) master.F90's compute / IO-server split -- FESIM is neither; the result
    !     is FESIM's own communicator.
    call MPI_Comm_split(MPI_COMM_WORLD, COLOR_FESIM_IN_MASTER, key, fesim_comm, ierr)
    ! join the FESOM compute tasks
    call ifs_mpmd_connect_yac_world(fesim_comm, is_fesim = .true.)
    call MPI_Comm_free(fesim_comm, ierr)
    if (key == 0 .or. yac_world_comm == MPI_COMM_NULL) then
       write(*,*) 'fesim_mpmd_join_ifs_world: joined the IFS MPI world, yac_world_comm valid: ', &
            yac_world_comm /= MPI_COMM_NULL
    end if
  end subroutine fesim_mpmd_join_ifs_world

end module ifs_mpmd_world

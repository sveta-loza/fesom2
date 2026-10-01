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
!   (2) [added by this module] MPI_Comm_split(MPI_COMM_WORLD) that builds the
!       YAC world = FESOM compute tasks + FESIM tasks, still inside mpp_io_init,
!       i.e. while every rank of MPI_COMM_WORLD is still in lockstep.
!   (3) master.F90: MPI_Barrier(MPI_COMM_WORLD) iff the environment variable
!       MPI_EPOCH is non-blank (raps always exports it).
!   (4) master.F90: one MPI_Comm_split(MPI_COMM_WORLD) (compute vs IFS IO servers).
! multio works on sub-communicators. A separate FESIM executable that shares
! MPI_COMM_WORLD with ifsMASTER (srun heterogeneous block / --multi-prog) must
! therefore mirror exactly these four calls, in this order, before anything
! else: that is fesim_mpmd_join_ifs_world.
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

  ! split (1) uses colours 1 (compute) and 3 (server); split (4) uses 1 and 2.
  integer, parameter :: COLOR_FESIM_IN_IOSPLIT = 7
  integer, parameter :: COLOR_YAC_WORLD        = 20
  integer, parameter :: COLOR_FESIM_IN_MASTER  = 9

  ! Communicator spanning FESOM compute tasks + FESIM tasks; MPI_COMM_NULL
  ! on ranks that are not part of it and whenever no FESIM is present.
  integer, save, public :: yac_world_comm = MPI_COMM_NULL

  public :: ifs_mpmd_fesim_ntasks, ifs_mpmd_have_yac_world
  public :: ifs_mpmd_split_yac_world, fesim_mpmd_join_ifs_world

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

  ! Split (2), IFS-executable side. Must be called by EVERY rank of the IFS
  ! executable (compute and multio-server tasks alike); only the compute
  ! tasks pass is_member=.true. and receive the YAC world communicator.
  subroutine ifs_mpmd_split_yac_world(is_member)
    logical, intent(in) :: is_member
    integer :: color, key, ierr
    color = MPI_UNDEFINED
    if (is_member) color = COLOR_YAC_WORLD
    call MPI_Comm_rank(MPI_COMM_WORLD, key, ierr)
    call MPI_Comm_split(MPI_COMM_WORLD, color, key, yac_world_comm, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(*,*) 'ifs_mpmd_split_yac_world: MPI_Comm_split failed, ierr=', ierr
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    if (.not. is_member) yac_world_comm = MPI_COMM_NULL
  end subroutine ifs_mpmd_split_yac_world

  ! FESIM-executable side: mirror the four world-collective calls (1)-(4) of
  ! ifsMASTER. Call right after MPI_Init and before any YAC call; afterwards
  ! yac_world_comm is the communicator to hand to yac_finit_comm.
  subroutine fesim_mpmd_join_ifs_world()
    integer :: ierr, key, dummy
    character(len=64) :: epoch
    call MPI_Comm_rank(MPI_COMM_WORLD, key, ierr)
    ! (1) mpp_io_init's compute/server split -- FESIM is neither
    call MPI_Comm_split(MPI_COMM_WORLD, COLOR_FESIM_IN_IOSPLIT, key, dummy, ierr)
    call MPI_Comm_free(dummy, ierr)
    ! (2) the YAC world
    call MPI_Comm_split(MPI_COMM_WORLD, COLOR_YAC_WORLD, key, yac_world_comm, ierr)
    ! (3) master.F90's MPI-startup-cost barrier, taken iff MPI_EPOCH is non-blank
    epoch = ' '
    call get_environment_variable('MPI_EPOCH', epoch)
    if (epoch /= ' ') call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! (4) master.F90's compute / IO-server split -- FESIM is neither
    call MPI_Comm_split(MPI_COMM_WORLD, COLOR_FESIM_IN_MASTER, key, dummy, ierr)
    call MPI_Comm_free(dummy, ierr)
    if (key == 0 .or. yac_world_comm == MPI_COMM_NULL) then
       write(*,*) 'fesim_mpmd_join_ifs_world: joined the IFS MPI world, yac_world_comm valid: ', &
            yac_world_comm /= MPI_COMM_NULL
    end if
  end subroutine fesim_mpmd_join_ifs_world

end module ifs_mpmd_world

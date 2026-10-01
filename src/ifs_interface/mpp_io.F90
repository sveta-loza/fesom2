!=====================================================
! Ocean output intialisation.
!
! -Original code for NEMOv40 by Kristian Mogensen, ECMWF.
! -Adapted to FESOM2 by Razvan Aguridan, ECMWF, 2023.
!-----------------------------------------------------

MODULE mpp_io
    USE mpi
#if defined(__MULTIO)        
    USE iom, only : iom_enable_multio, iom_initialize, iom_init_server, iom_finalize
#endif
#if defined(__yac)
    ! TerraDT config #1: a separate FESIM executable shares MPI_COMM_WORLD with
    ! ifsMASTER (see ifs_mpmd_world.F90 for the full choreography).
    USE ifs_mpmd_world, ONLY : ifs_mpmd_fesim_ntasks, ifs_mpmd_split_yac_world
#endif
    IMPLICIT NONE
    PRIVATE

    PUBLIC &
        & mpp_io_init, &
        & mpp_io_init_2, &
        & mpp_stop

    INTEGER :: ntask_multio  = 0
    INTEGER :: ntask_xios    = 0
    INTEGER :: ntask_fesim   = 0   ! TerraDT: FESIM tasks appended to MPI_COMM_WORLD (env FESIM_NTASKS)
    LOGICAL, PUBLIC :: lioserver, lmultioserver, lmultiproc
    INTEGER :: ntask_notio
    INTEGER, SAVE :: mppallrank, mppallsize, mppiorank, mppiosize
    INTEGER, SAVE :: mppmultiorank, mppmultiosize
    INTEGER, SAVE :: mppcomprank, mppcompsize
    INTEGER, SAVE :: pcommworld, pcommworldmultio
    
    CONTAINS
    
    SUBROUTINE mpp_io_init( iicomm,  lio, irequired, iprovided, lmpi1 ) 

        INTEGER, INTENT(INOUT) :: iicomm
        LOGICAL, INTENT(INOUT) :: lio
        INTEGER, INTENT(INOUT) :: irequired, iprovided
        LOGICAL, INTENT(IN) :: lmpi1

        INTEGER :: icode, ierr, icolor
        LOGICAL :: mpi_called
        CHARACTER(len=128) :: cdlogfile
        INTEGER :: ji
        NAMELIST/namio/ntask_multio,ntask_xios
        
        CALL mpi_initialized( mpi_called, icode )
        IF ( icode /= MPI_SUCCESS ) THEN
            WRITE(*,*)' mpp_io_init: Error in routine mpi_initialized'
            CALL mpi_abort( mpi_comm_world, icode, ierr )
        ENDIF

        IF(  mpi_called ) THEN
            WRITE(*,*)' mpi_io_init assumes that it is initialising MPI'
            CALL mpi_abort( mpi_comm_world, 1, ierr )
        ENDIF

        IF (lmpi1) THEN
            CALL mpi_init( icode )
        ELSE
#ifdef MPI1
            WRITE(0,*)'mpp_io_init:'
            WRITE(0,*)'MPI1 defined but lmpi1 is false'
            CALL abort
#else
            CALL mpi_init_thread(irequired,iprovided,icode)
#endif
        ENDIF

        IF ( icode /= MPI_SUCCESS ) THEN
            WRITE(*,*)' mpp_io_init: Error in routine mpi_init'
            CALL mpi_abort( mpi_comm_world, icode, ierr )
        ENDIF

        CALL mpi_comm_rank( mpi_comm_world, mppallrank, ierr )
        CALL mpi_comm_size( mpi_comm_world, mppallsize, ierr )

        OPEN(10,file='namio.in')
        READ(10,namio)
        WRITE(*,namio)
        CLOSE(10)

#if defined(__MULTIO)
        IF (ntask_multio /= 0) THEN
            CALL iom_enable_multio()
        ENDIF
#endif

#if defined(__yac)
        ! TerraDT config #1: the LAST ntask_fesim ranks of MPI_COMM_WORLD belong to
        ! the FESIM executable. They take part in the split below with their own
        ! colour (see ifs_mpmd_world) and must not be counted as compute tasks.
        ntask_fesim = ifs_mpmd_fesim_ntasks()
        IF (mppallrank == 0 .AND. ntask_fesim > 0) &
            WRITE(*,*) ' mpp_io_init: FESIM_NTASKS =', ntask_fesim, &
                       ' -> the last', ntask_fesim, ' ranks of MPI_COMM_WORLD are the FESIM executable'
#endif

        IF ( ntask_xios + ntask_multio + ntask_fesim == 0 ) THEN
            iicomm = mpi_comm_world
            lio=.FALSE.
            RETURN
        ENDIF

        ntask_notio = mppallsize - ntask_xios - ntask_multio - ntask_fesim
        IF ((mppallrank+1)<=ntask_notio) THEN
            icolor=1
            lioserver=.FALSE.
            lmultioserver=.FALSE.
        ELSE
            icolor=3
            lioserver=.TRUE.
            lmultioserver=.TRUE.
        ENDIF
        lio=lioserver

        CALL mpi_comm_split( mpi_comm_world, icolor, 0, iicomm, icode )
        IF ( icode /= MPI_SUCCESS ) THEN
            WRITE(*,*)' mpp_io_init: Error in routine mpi_comm_split'
            CALL mpi_abort( mpi_comm_world, icode, ierr )
        ENDIF
#if defined(__yac)
        ! TerraDT config #1: build the YAC world (compute tasks + FESIM) while all
        ! ranks of MPI_COMM_WORLD are still in lockstep; FESIM mirrors this split.
        IF (ntask_fesim > 0) CALL ifs_mpmd_split_yac_world(is_member = .NOT. lioserver)
#endif
        IF (lioserver) THEN
            CALL mpi_comm_rank( iicomm, mppiorank, ierr )
            CALL mpi_comm_size( iicomm, mppiosize, ierr )
            WRITE(cdlogfile,'(A,I4.4,A)')'nemo_io_server.',mppiorank,'.log'
        ELSE
            mppiorank=0
            mppiosize=0
        ENDIF
        lio=lioserver

    END SUBROUTINE mpp_io_init

    SUBROUTINE mpp_io_init_2( iicomm )
        
        INTEGER, INTENT(INOUT) :: iicomm

        INTEGER :: icode, ierr, icolor, iicommx, iicommm, iicommo
        INTEGER :: ji,inum
        LOGICAL :: lcompp

        ! Construct multio server communicator

        IF (lmultioserver.OR..NOT.lioserver) THEN
            icolor=12
        ELSE
            icolor=13
        ENDIF

        CALL mpi_comm_split( iicomm, icolor, 0, pcommworldmultio, icode )
        IF ( icode /= MPI_SUCCESS ) THEN
            WRITE(*,*)' mpp_io_init2: Error in routine mpi_comm_split'
            CALL mpi_abort( mpi_comm_world, icode, ierr )
        ENDIF

        CALL mpi_comm_rank( pcommworldmultio, mppmultiorank, ierr )
        CALL mpi_comm_size( pcommworldmultio, mppmultiosize, ierr )

        ! Construct compute communicator

        IF (.NOT.lioserver) THEN
            icolor=14
            lcompp=.TRUE.
        ELSE
            icolor=15
            lcompp=.FALSE.
        ENDIF

        CALL mpi_comm_split( iicomm, icolor, 0, iicommo, icode )
        IF ( icode /= MPI_SUCCESS ) THEN
            WRITE(*,*)' mpp_io_init2: Error in routine mpi_comm_split'
            CALL mpi_abort( mpi_comm_world, icode, ierr )
        ENDIF
        
        CALL mpi_comm_rank( iicommo, mppcomprank, ierr )
        CALL mpi_comm_size( iicommo, mppcompsize, ierr )

#if defined(__MULTIO)
        IF (.NOT.lioserver) THEN
            CALL iom_initialize( "for_xios_mpi_id", return_comm=iicommm, global_comm = pcommworldmultio )    ! nemo local communicator given by xios
        ELSE
            ! For io-server tasks start an run the right server
            CALL iom_init_server( server_comm = pcommworldmultio )
        ENDIF
#endif

        ! Return to the model with iicomm being compute only tasks
        iicomm = iicommo

    END SUBROUTINE mpp_io_init_2

    SUBROUTINE mpp_stop
        INTEGER :: ierr
        
#if defined(__MULTIO)
        IF (.NOT.lioserver) THEN
            call iom_finalize()
        ENDIF
#endif

        CALL mpi_finalize( ierr )
    END SUBROUTINE mpp_stop

END MODULE mpp_io

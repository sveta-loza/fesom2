! Cold start of the sea ice in the sea-ice component (FESIM).
!
! The ocean with built-in sea ice seeds the initial ice from its initial SST
! (ice_initial_state reads tracers%data(1)). FESIM has no ocean state at
! initialisation: the SST arrives from the ocean once the run starts. The cold
! start is therefore completed from the first SST received, which is the
! ocean's initial condition -- the same information the built-in ice uses.
!
! ice_initial_state (ice_setup_step.F90) marks the cold start pending instead
! of seeding from the empty tracers; the exchange with the ocean
! (fesim_ocean_exchange.F90) calls seed_ice_from_sst after the first SST
! receive. Restarts and file-initialised ice are untouched.
module fesim_ice_init
  use o_PARAM, only: WP
  implicit none
  private

  logical, public, save :: ice_cold_start_pending = .false.

  public :: seed_ice_from_sst

contains

  ! Same rule as ice_initial_state: where the sea surface is below 0 degC,
  ! 1 m of ice under 0.1 m of snow in the north, 2 m under 0.5 m in the south,
  ! at 90 % concentration.
  subroutine seed_ice_from_sst(ice, sst, partit, mesh)
    use MOD_ICE
    use MOD_PARTIT
    use MOD_MESH
    use mpi
    type(t_ice),    intent(inout), target :: ice
    real(kind=WP),  intent(in)            :: sst(:)        ! degC, node-sized
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
    integer :: i, nseeded, ierr

    a_ice  => ice%data(1)%values(:)
    m_ice  => ice%data(2)%values(:)
    m_snow => ice%data(3)%values(:)

    nseeded = 0
    do i = 1, partit%myDim_nod2D+partit%eDim_nod2D
       if (mesh%ulevels_nod2d(i) > 1) cycle          ! cavity: no sea ice
       if (sst(i) < 0.0_WP) then
          if (mesh%geo_coord_nod2D(2,i) > 0._WP) then
             m_ice(i)  = 1.0_WP
             m_snow(i) = 0.1_WP
          else
             m_ice(i)  = 2.0_WP
             m_snow(i) = 0.5_WP
          end if
          a_ice(i)    = 0.9_WP
          ice%uice(i) = 0.0_WP
          ice%vice(i) = 0.0_WP
          if (i <= partit%myDim_nod2D) nseeded = nseeded + 1
       end if
    end do
    call MPI_Allreduce(MPI_IN_PLACE, nseeded, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype == 0) write(*,*) 'FESIM cold start: sea ice seeded from the first ocean SST at', &
                                     nseeded, 'nodes'
    ice_cold_start_pending = .false.
  end subroutine seed_ice_from_sst

end module fesim_ice_init

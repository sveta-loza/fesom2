module cpl_fields
  !============================================================================
  ! The fields FESOM exchanges with its external components, by name.
  !
  ! define_cpl_fields registers the field set of the configured partner once,
  ! in the order the coupler announces and exchanges them. The send and
  ! receive loops dispatch on the name instead of the position, so a field
  ! keeps its meaning whatever optional fields are registered before it, and
  ! nsend/nrecv follow from what was registered.
  !
  ! Modelled on CPLNG2_ADD_FLD in OpenIFS. The names are the ones the
  ! namcouple (OASIS) or coupling.yaml (YAC) refers to.
  !============================================================================
  implicit none
  save
  private

  public :: define_cpl_fields
  public :: nsend, nrecv, cpl_send, cpl_recv
  public :: cpl_send_collection_size, cpl_recv_collection_size
  public :: cpl_send_partner, cpl_recv_partner
  public :: cpl_send_steps, cpl_recv_steps
  public :: CPL_PARTNER_ATMOSPHERE, CPL_PARTNER_SEA_ICE, CPL_PARTNER_OCEAN

  integer                                      :: nsend = 0
  integer                                      :: nrecv = 0
  character(len=32), allocatable, dimension(:) :: cpl_send
  character(len=32), allocatable, dimension(:) :: cpl_recv
  ! Number of 2D fields bundled under one name. Always 1 for OASIS; YAC
  ! exchanges some fields as collections.
  integer,           allocatable, dimension(:) :: cpl_send_collection_size
  integer,           allocatable, dimension(:) :: cpl_recv_collection_size
  ! The component a field goes to or comes from. The atmosphere and the sea
  ! ice (FESIM) are exchanged at different points of the time step, so each
  ! exchange routine handles the fields of its own partner only.
  integer, parameter :: CPL_PARTNER_ATMOSPHERE = 1
  integer, parameter :: CPL_PARTNER_SEA_ICE    = 2   ! seen from the ocean
  integer, parameter :: CPL_PARTNER_OCEAN      = 3   ! seen from the sea ice
  integer,           allocatable, dimension(:) :: cpl_send_partner
  integer,           allocatable, dimension(:) :: cpl_recv_partner
  ! Model time steps between two exchanges of the field (YAC registers the
  ! field with this multiple of dt). 1 for the atmosphere; cpl_stride for the
  ! ocean <-> sea-ice exchange.
  integer,           allocatable, dimension(:) :: cpl_send_steps
  integer,           allocatable, dimension(:) :: cpl_recv_steps

contains

  !____________________________________________________________________________
  !> Register the field set of the configured partner. Called once, after the
  !> namelists are read and before the coupler definition phase.
  subroutine define_cpl_fields()
    use cpl_config, only: is_coupled_to_echam, is_coupled_to_oifs, &
                          is_coupled_to_icon_a, is_coupled_to_ifs, &
                          is_coupled_to_fesim, is_coupled_to_fesom
    use g_config,   only: lwiso, use_icebergs, use_atm_ice_tskin, use_atm_mslp, &
                          cpl_stride

    nsend = 0
    nrecv = 0
    allocate(cpl_send(0), cpl_recv(0))
    allocate(cpl_send_collection_size(0), cpl_recv_collection_size(0))
    allocate(cpl_send_partner(0), cpl_recv_partner(0))
    allocate(cpl_send_steps(0), cpl_recv_steps(0))

    ! The sea ice talks to the ocean only; its atmosphere flags name what the
    ! ocean forwards (see cpl_config).
    if (is_coupled_to_fesom) then
       call define_sea_ice_side_fields(is_coupled_to_icon_a, is_coupled_to_ifs, cpl_stride)
       return
    end if

    if (is_coupled_to_oifs) then
       ! OpenIFS (AWI-CM3, AWI-ESM3) over OASIS3-MCT 5
       call add_send('sst_feom')   ! sea surface temperature [K]
       call add_send('sie_feom')   ! sea ice concentration [0-1]
       call add_send('snt_feom')   ! snow thickness [m]
       call add_send('ist_feom')   ! ice surface temperature * concentration [K]
       call add_send('sia_feom')   ! sea ice albedo [0-1]
       call add_send('u_feom')     ! eastward surface velocity [m/s]
       call add_send('v_feom')     ! northward surface velocity [m/s]
       call add_send('sit_feom')   ! effective sea ice thickness [m]
#if defined (__recom)
       call add_send('FCO2_feom')  ! CO2 flux [kgCO2 m-2 s-1]
#endif
       call add_common_atm_recv()
       call add_recv('calv_oce')   ! solid runoff (calving)
       call add_recv('u10w_oce')   ! zonal 10 m wind
       call add_recv('v10w_oce')   ! meridional 10 m wind
#if defined (__recom)
       call add_recv('XCO2_oce')   ! atmospheric CO2 mass mixing ratio
#endif
       if (use_atm_ice_tskin) call add_recv('tsk_ico')  ! ice-tile skin temperature
       if (use_atm_mslp)      call add_recv('mslp_oce') ! mean sea-level pressure

    else if (is_coupled_to_echam) then
       ! ECHAM6 (AWI-CM2, AWI-ESM2) over OASIS3-MCT 2.8
       call add_send('sst_feom')   ! sea surface temperature [degC]
       call add_send('sit_feom')   ! effective sea ice thickness [m]
       call add_send('sie_feom')   ! sea ice concentration [0-1]
       call add_send('snt_feom')   ! snow thickness [m]
       if (lwiso) then
          call add_send('o18w_oce') ! delta 18O of surface water
          call add_send('hdow_oce') ! delta D of surface water
          call add_send('o16w_oce') ! delta H216O of surface water
          call add_send('o18i_oce') ! delta 18O of sea ice
          call add_send('hdoi_oce') ! delta D of sea ice
          call add_send('o16i_oce') ! delta H216O of sea ice
       end if
       call add_common_atm_recv()
       if (lwiso) then
          call add_recv('w1_oce')  ! H218O precipitation over water
          call add_recv('w2_oce')  ! HDO precipitation over water
          call add_recv('w3_oce')  ! H216O precipitation over water
          call add_recv('i1_oce')  ! H218O snowfall over sea ice
          call add_recv('i2_oce')  ! HDO snowfall over sea ice
          call add_recv('i3_oce')  ! H216O snowfall over sea ice
       end if
       if (use_icebergs) then
          call add_recv('u10w_oce') ! zonal 10 m wind
          call add_recv('v10w_oce') ! meridional 10 m wind
       end if

    else if (is_coupled_to_icon_a) then
       ! ICON-A over YAC
       call add_send('sst_feom')                     ! sea surface temperature [K]
       call add_send('ocean_sea_ice_bundle', 3)      ! ice thickness, snow thickness, concentration
       call add_recv('taux', 2)                      ! zonal stress over water, ice
       call add_recv('tauy', 2)                      ! meridional stress over water, ice
       call add_recv('surface_fresh_water_flux', 3)  ! rain, snow, evaporation
       call add_recv('total_heat_flux', 4)           ! shortwave, longwave, sensible, latent
       call add_recv('atmosphere_sea_ice_bundle', 2) ! heat flux into the ice
       call add_recv('river_runoff')                 ! river runoff
    end if

    if (is_coupled_to_fesim) &
       call define_ocean_side_ice_fields(is_coupled_to_icon_a, is_coupled_to_ifs, cpl_stride)
  end subroutine define_cpl_fields

  !____________________________________________________________________________
  !> The ocean <-> sea-ice (FESIM) fields, seen from the ocean. Exchanged every
  !> `steps` ocean steps.
  !>
  !> The native ocean state sent to the sea ice is the same for every
  !> atmosphere; what the ocean forwards from the atmosphere depends on it:
  !>  - ICON-A: the atmospheric fluxes the ocean received, as received.
  !>  - IFS: the fluxes the IFS interface deposited in the ocean's arrays,
  !>    already in FESOM units and rotated, as one 10-component bundle. The
  !>    sea ice then also returns its surface temperature and albedo, which the
  !>    ocean relays to IFS.
  !>  - none: the raw atmospheric state from the ocean's forcing files, as one
  !>    9-component bundle; the sea ice computes its fluxes itself.
  !> Forwarded fields travel in their source units; the sea ice converts.
  subroutine define_ocean_side_ice_fields(icon_a, ifs, steps)
    logical, intent(in) :: icon_a, ifs
    integer, intent(in) :: steps
    integer, parameter  :: p = CPL_PARTNER_SEA_ICE

    call add_send('sst_feom_to_ice',     1, p, steps)  ! SST [degC]
    call add_send('ocean_to_ice_bundle', 2, p, steps)  ! SSS [psu], SSH [m]
    call add_send('ocean_to_ice_uv',     2, p, steps)  ! surface velocity [m/s]
    call add_send('runoff_to_ice',       1, p, steps)  ! river runoff [m/s]
    if (icon_a) then
       call add_send('taux_to_ice',                      2, p, steps)
       call add_send('tauy_to_ice',                      2, p, steps)
       call add_send('surface_fresh_water_flux_to_ice',  3, p, steps)
       call add_send('total_heat_flux_to_ice',           4, p, steps)
       call add_send('atmosphere_sea_ice_bundle_to_ice', 2, p, steps)
    else if (ifs) then
       call add_send('atm_ice_flux_to_ice', 10, p, steps)
    else
       call add_send('atm_state_to_ice',     9, p, steps)
    end if

    call add_ice_to_ocean_fields(.false., ifs, steps, p)
  end subroutine define_ocean_side_ice_fields

  !____________________________________________________________________________
  !> The same fields seen from the sea ice: what the ocean sends is received
  !> and vice versa. Every field of the sea ice is exchanged with the ocean.
  subroutine define_sea_ice_side_fields(icon_a, ifs, steps)
    logical, intent(in) :: icon_a, ifs
    integer, intent(in) :: steps
    integer, parameter  :: p = CPL_PARTNER_OCEAN

    call add_ice_to_ocean_fields(.true., ifs, steps, p)

    call add_recv('sst_feom_to_ice',     1, p, steps)
    call add_recv('ocean_to_ice_bundle', 2, p, steps)
    call add_recv('ocean_to_ice_uv',     2, p, steps)
    call add_recv('runoff_to_ice',       1, p, steps)
    if (icon_a) then
       call add_recv('taux_to_ice',                      2, p, steps)
       call add_recv('tauy_to_ice',                      2, p, steps)
       call add_recv('surface_fresh_water_flux_to_ice',  3, p, steps)
       call add_recv('total_heat_flux_to_ice',           4, p, steps)
       call add_recv('atmosphere_sea_ice_bundle_to_ice', 2, p, steps)
    else if (ifs) then
       call add_recv('atm_ice_flux_to_ice', 10, p, steps)
    else
       call add_recv('atm_state_to_ice',     9, p, steps)
    end if
  end subroutine define_sea_ice_side_fields

  !____________________________________________________________________________
  !> The fields the sea ice returns to the ocean: sent by the sea ice
  !> (as_send), received by the ocean.
  !>  - sea_ice_bundle: m_ice, m_snow, a_ice, [ice_temp, ice_alb under IFS,
  !>    for the relay to IFS], uice, vice (the ocean's under-ice mixing needs
  !>    the ice drift)
  !>  - ice_to_ocean_stress: ice-ocean drag x, y [Pa]
  !>  - ice_to_ocean_flux: net heat [W/m2] and freshwater [m/s] flux
  !>  - ice_to_ocean_thermo: evaporation, ice_sublimation, thdgr, thdgrsn,
  !>    a_ice before the thermodynamics, brine salt flux -- the terms of the
  !>    ice thermodynamics the ocean's freshwater balancing integrates
  subroutine add_ice_to_ocean_fields(as_send, ifs, steps, p)
    logical, intent(in) :: as_send, ifs
    integer, intent(in) :: steps, p
    character(len=32)   :: names(4)
    integer             :: sizes(4), i

    names = [character(len=32) :: 'sea_ice_bundle', 'ice_to_ocean_stress', &
                                  'ice_to_ocean_flux', 'ice_to_ocean_thermo']
    sizes = [5, 2, 2, 6]
    if (ifs) sizes(1) = 7
    do i = 1, 4
       if (as_send) then
          call add_send(trim(names(i)), sizes(i), p, steps)
       else
          call add_recv(trim(names(i)), sizes(i), p, steps)
       end if
    end do
  end subroutine add_ice_to_ocean_fields

  !____________________________________________________________________________
  !> The atmosphere-to-ocean fields ECHAM6 and OpenIFS both send, in the order
  !> both namcouples use.
  subroutine add_common_atm_recv()
    call add_recv('taux_oce')   ! zonal stress over water
    call add_recv('tauy_oce')   ! meridional stress over water
    call add_recv('taux_ico')   ! zonal stress over ice
    call add_recv('tauy_ico')   ! meridional stress over ice
    call add_recv('prec_oce')   ! total precipitation
    call add_recv('snow_oce')   ! snowfall
    call add_recv('evap_oce')   ! evaporation over water
    call add_recv('subl_oce')   ! sublimation over ice
    call add_recv('heat_oce')   ! net heat flux over water
    call add_recv('heat_ico')   ! net heat flux over ice
    call add_recv('heat_swo')   ! shortwave over water
    call add_recv('hydr_oce')   ! runoff
  end subroutine add_common_atm_recv

  !____________________________________________________________________________
  !> partner defaults to the atmosphere and steps to 1 (every time step).
  subroutine add_send(name, collection_size, partner, steps)
    character(len=*), intent(in)           :: name
    integer,          intent(in), optional :: collection_size, partner, steps

    call append(cpl_send, cpl_send_collection_size, nsend, name, collection_size)
    call append_tag(cpl_send_partner, CPL_PARTNER_ATMOSPHERE, partner)
    call append_tag(cpl_send_steps,   1,                      steps)
  end subroutine add_send

  !____________________________________________________________________________
  subroutine add_recv(name, collection_size, partner, steps)
    character(len=*), intent(in)           :: name
    integer,          intent(in), optional :: collection_size, partner, steps

    call append(cpl_recv, cpl_recv_collection_size, nrecv, name, collection_size)
    call append_tag(cpl_recv_partner, CPL_PARTNER_ATMOSPHERE, partner)
    call append_tag(cpl_recv_steps,   1,                      steps)
  end subroutine add_recv

  !____________________________________________________________________________
  subroutine append(names, sizes, n, name, collection_size)
    character(len=32), allocatable, intent(inout) :: names(:)
    integer,           allocatable, intent(inout) :: sizes(:)
    integer,                        intent(inout) :: n
    character(len=*),               intent(in)    :: name
    integer,          optional,     intent(in)    :: collection_size

    integer :: size_of_this

    size_of_this = 1
    if (present(collection_size)) size_of_this = collection_size

    names = [character(len=32) :: names, name]
    sizes = [sizes, size_of_this]
    n     = n + 1
  end subroutine append

  !____________________________________________________________________________
  subroutine append_tag(tags, default_value, value)
    integer, allocatable, intent(inout) :: tags(:)
    integer,              intent(in)    :: default_value
    integer, optional,    intent(in)    :: value

    if (present(value)) then
       tags = [tags, value]
    else
       tags = [tags, default_value]
    end if
  end subroutine append_tag

end module cpl_fields

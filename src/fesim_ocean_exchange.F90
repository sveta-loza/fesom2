! The sea-ice component (FESIM): its exchange with the ocean over YAC.
!
! The counterpart of exchange_oce_ice_yac (gen_forcing_couple.F90) on the
! ocean side. FESIM sends its ice state, the ice-ocean drag, its heat and
! freshwater flux and the thermodynamic terms of the ocean's freshwater
! balance; it receives the ocean's surface state (SST, SSS, SSH, velocity),
! the river runoff and what the ocean forwards from its atmosphere, which
! follows the atmosphere flag set together with is_coupled_to_fesom in
! namelist.cpl:
!   is_coupled_to_icon_a  ICON fluxes, as ICON sent them (geographic frame,
!                         forcing units), converted here as the ocean does
!   is_coupled_to_ifs     IFS fluxes, already in FESOM units and rotated
!   neither               the raw atmospheric state from the ocean's forcing
!                         files; the bulk formulae below turn it into wind
!                         stresses, the bulk ice thermodynamics into fluxes
! The fields are the ones define_cpl_fields registered (cpl_fields.F90). On
! a cold start the sea ice is seeded from the first SST received.
module fesim_ocean_exchange_module
  use o_PARAM
  use MOD_MESH
  use MOD_PARTIT
  use MOD_TRACER
  use MOD_ICE
  use MOD_DYN
  use o_arrays
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_comm_auto
  use g_rotate_grid
  use g_sbf,              only: l_mslp, l_xstre, l_ystre
  use gen_bulk,           only: ncar_ocean_fluxes_mode, cal_wind_drag_coeff
  use par_support_module, only: par_abort
#if defined (__cpl_yac)
  use cpl_yac_driver,     only: nsend, nrecv, cpl_send, cpl_recv, &
                                cpl_send_collection_size, cpl_recv_collection_size, &
                                cpl_yac_send, cpl_yac_recv
  use cpl_fields,         only: cpl_send_partner, cpl_recv_partner, CPL_PARTNER_OCEAN
  use cpl_config,         only: is_coupled_to_icon_a, is_coupled_to_ifs, is_coupled_to_atmosphere
  use fesim_ice_init,     only: ice_cold_start_pending, seed_ice_from_sst
#endif
  implicit none
  private
  public :: exchange_ice_oce_yac

contains

#if defined (__cpl_yac)
subroutine exchange_ice_oce_yac(istep, ice, tracers, dynamics, partit, mesh)
  integer,        intent(in)            :: istep
  type(t_ice)   , intent(inout), target :: ice
  type(t_tracer), intent(in),    target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  type(t_dyn)   , intent(in),    target :: dynamics
  !_____________________________________________________________________________
  integer                  :: i, n, nod2
  real(kind=WP)            :: aux, dux, dvy
  real(kind=WP), dimension(:,:), allocatable , save  :: exchange
  logical                                          :: action
  logical                                          :: do_rotate_oce_wind=.false.
  logical                                          :: do_rotate_ice_wind=.false.
  real(kind=WP), dimension(:), pointer  :: u_ice, v_ice, u_w, v_w
  real(kind=WP), dimension(:), pointer  :: stress_atmice_x, stress_atmice_y
  real(kind=WP), dimension(:), pointer  :: oce_heat_flux, ice_heat_flux
  real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
  real(kind=WP), dimension(:), pointer  :: t_oce, s_oce, eta_n
  real(kind=WP), dimension(:), pointer  :: ice_temp, ice_alb, enthalpyoffuse
  real(kind=WP)              , pointer  :: rhoair

  nod2             =  partit%myDim_nod2D
  u_ice            => ice%uice(:)
  v_ice            => ice%vice(:)
  u_w              => ice%srfoce_u(:)
  v_w              => ice%srfoce_v(:)
  stress_atmice_x  => ice%stress_atmice_x(:)
  stress_atmice_y  => ice%stress_atmice_y(:)
  a_ice            => ice%data(1)%values(:)
  m_ice            => ice%data(2)%values(:)
  m_snow           => ice%data(3)%values(:)
  oce_heat_flux    => ice%atmcoupl%oce_flx_h(:)
  ice_heat_flux    => ice%atmcoupl%ice_flx_h(:)
  rhoair           => ice%thermo%rhoair
  t_oce            => ice%srfoce_temp(:)
  s_oce            => ice%srfoce_salt(:)
  eta_n            => ice%srfoce_ssh(:)
  ice_alb          => ice%atmcoupl%ice_alb(:)
  enthalpyoffuse   => ice%atmcoupl%enthalpyoffuse(:)
  ! the ice surface temperature is a tracer when FESOM provides it (IFS);
  ! it is relayed to the ocean then
  if (ice%ist_itracer_idx > 0) ice_temp => ice%data(ice%ist_itracer_idx)%values(:)

  ! NB: the wind-stress arrays are not zeroed here: between coupling instants
  ! they keep the last received values (zero-order hold) instead of leaving
  ! the sea ice without wind stress. They are zeroed once at setup.
  if (.NOT. ALLOCATED(exchange)) then
     ALLOCATE(exchange(nod2, MAX(MAXVAL(cpl_send_collection_size), MAXVAL(cpl_recv_collection_size))))
     exchange = 0
  end if

  ! ---- send to the ocean -----------------------------------------------------
  do i=1,nsend
     if (cpl_send_partner(i) /= CPL_PARTNER_OCEAN) cycle
     exchange = 0.
     select case (trim(cpl_send(i)))
     case ('sea_ice_bundle')
        exchange(:,1) = m_ice(1:nod2)              ! ice thickness [m]
        exchange(:,2) = m_snow(1:nod2)             ! snow thickness [m]
        exchange(:,3) = a_ice(1:nod2)              ! ice concentration [0-1]
        if (is_coupled_to_ifs) then
           exchange(:,4) = ice_temp(1:nod2)        ! ice surface temperature [K], for IFS
           exchange(:,5) = ice_alb(1:nod2)         ! ice albedo [0-1], for IFS
           exchange(:,6) = ice%uice(1:nod2)        ! ice drift [m/s], for the ocean's
           exchange(:,7) = ice%vice(1:nod2)        !   under-ice mixing
        else
           exchange(:,4) = ice%uice(1:nod2)
           exchange(:,5) = ice%vice(1:nod2)
        end if
     case ('ice_to_ocean_stress')
        ! ice-ocean drag, formed by oce_fluxes_mom from the ice velocity and the
        ! received ocean velocity; in the rotated model frame
        exchange(:,1) = ice%stress_iceoce_x(1:nod2) ! [Pa]
        exchange(:,2) = ice%stress_iceoce_y(1:nod2) ! [Pa]
     case ('ice_to_ocean_flux')
        ! net heat and freshwater flux of the ice thermodynamics; the ocean
        ! applies them in oce_fluxes (heat_flux=-flx_h, water_flux=-flx_fw)
        exchange(:,1) = ice%flx_h(1:nod2)   ! [W/m2]
        exchange(:,2) = ice%flx_fw(1:nod2)  ! [m/s]
     case ('ice_to_ocean_thermo')
        ! the terms of the ice thermodynamics the ocean's global freshwater
        ! balancing (oce_fluxes) integrates
        exchange(:,1) = evaporation(1:nod2)            ! [m/s], negative up
        exchange(:,2) = ice_sublimation(1:nod2)        ! [m/s]
        exchange(:,3) = ice%thermo%thdgr(1:nod2)       ! ice growth rate  [m/s]
        exchange(:,4) = ice%thermo%thdgrsn(1:nod2)     ! snow growth rate [m/s]
        exchange(:,5) = ice%data(1)%values_old(1:nod2) ! a_ice before the thermodynamics
        exchange(:,6) = real_salt_flux(1:nod2)         ! brine salt flux [m/s psu] (0 with use_virt_salt)
     case default
        write(*,*) 'exchange_ice_oce_yac: no handler for the sent field ', trim(cpl_send(i))
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
     end select
     call cpl_yac_send(i, exchange(:,1:cpl_send_collection_size(i)), action)
     if (flag_debug .and. partit%mype==0) write(*,*) 'FESIM SEND: ', trim(cpl_send(i)), &
                                                     ' max val:', maxval(exchange), ' . ACTION? ', action
  enddo

  ! ---- receive from the ocean -------------------------------------------------
  do i=1,nrecv
     if (cpl_recv_partner(i) /= CPL_PARTNER_OCEAN) cycle
     exchange = 0.0
     call cpl_yac_recv(i, exchange(:,1:cpl_recv_collection_size(i)), action)
     if (flag_debug .and. partit%mype==0) write(*,*) 'FESIM RECV: ', trim(cpl_recv(i)), &
                                                     ' max val:', maxval(exchange), ' . ACTION? ', action
     if (.not. action) cycle
     select case (trim(cpl_recv(i)))
     case ('sst_feom_to_ice')
        t_oce(1:nod2) = exchange(:,1)   ! SST [degC], the ocean's native unit
        call exchange_nod(t_oce, partit)
        if (ice_cold_start_pending) call seed_ice_from_sst(ice, t_oce, partit, mesh)
     case ('runoff_to_ice')
        ! river runoff [m/s]; the ice thermodynamics adds it to the freshwater
        ! flux, which is how it reaches the ocean with the built-in sea ice
        runoff(1:nod2) = exchange(:,1)
        call exchange_nod(runoff, partit)
     case ('ocean_to_ice_bundle')
        s_oce(1:nod2) = exchange(:,1)   ! SSS [psu]
        eta_n(1:nod2) = exchange(:,2)   ! SSH [m]
        call exchange_nod(s_oce, partit)
        call exchange_nod(eta_n, partit)
     case ('ocean_to_ice_uv')
        ! in the rotated model frame already: both components share the mesh
        u_w(1:nod2) = exchange(:,1)
        v_w(1:nod2) = exchange(:,2)
        call exchange_nod(u_w, partit)
        call exchange_nod(v_w, partit)
     ! ICON fluxes forwarded by the ocean: converted as the ocean converts them
     case ('taux_to_ice')
        stress_atmoce_x(1:nod2) = exchange(:,1)
        call exchange_nod(stress_atmoce_x, partit)
        stress_atmice_x(1:nod2) = exchange(:,2)
        call exchange_nod(stress_atmice_x, partit)
        do_rotate_oce_wind=.true.
        do_rotate_ice_wind=.true.
     case ('tauy_to_ice')
        stress_atmoce_y(1:nod2) = exchange(:,1)
        call exchange_nod(stress_atmoce_y, partit)
        stress_atmice_y(1:nod2) = exchange(:,2)
        call exchange_nod(stress_atmice_y, partit)
        do_rotate_oce_wind=.true.
        do_rotate_ice_wind=.true.
     case ('surface_fresh_water_flux_to_ice')
        prec_rain(1:nod2)     = exchange(:,1)/1000   ! kg m-2 s-1 -> m/s
        call exchange_nod(prec_rain, partit)
        prec_snow(1:nod2)     = exchange(:,2)/1000
        call exchange_nod(prec_snow, partit)
        evap_no_ifrac(1:nod2) = exchange(:,3)/1000
        call exchange_nod(evap_no_ifrac, partit)
     case ('total_heat_flux_to_ice')
        oce_heat_flux(1:nod2) = exchange(:,2) + exchange(:,3) + exchange(:,4)
        call exchange_nod(oce_heat_flux, partit)
        shortwave(1:nod2)     = exchange(:,1)
        call exchange_nod(shortwave, partit)
     case ('atmosphere_sea_ice_bundle_to_ice')
        ice_heat_flux(1:nod2) = exchange(:,2) + exchange(:,1)
        call exchange_nod(ice_heat_flux, partit)
     ! IFS fluxes forwarded by the ocean, in FESOM units and rotated: 1:1
     case ('atm_ice_flux_to_ice')
        stress_atmice_x(1:nod2) = exchange(:,1)   ! [Pa]
        stress_atmice_y(1:nod2) = exchange(:,2)   ! [Pa]
        oce_heat_flux(1:nod2)   = exchange(:,3)   ! [W/m2]
        ice_heat_flux(1:nod2)   = exchange(:,4)   ! [W/m2]
        shortwave(1:nod2)       = exchange(:,5)   ! [W/m2]
        prec_rain(1:nod2)       = exchange(:,6)   ! [m/s]
        prec_snow(1:nod2)       = exchange(:,7)   ! [m/s]
        evap_no_ifrac(1:nod2)   = exchange(:,8)   ! [m/s]
        sublimation(1:nod2)     = exchange(:,9)   ! [m/s]
        enthalpyoffuse(1:nod2)  = exchange(:,10)  ! [W/m2]
        call exchange_nod(stress_atmice_x, partit)
        call exchange_nod(stress_atmice_y, partit)
        call exchange_nod(oce_heat_flux, partit)
        call exchange_nod(ice_heat_flux, partit)
        call exchange_nod(shortwave, partit)
        call exchange_nod(prec_rain, partit)
        call exchange_nod(prec_snow, partit)
        call exchange_nod(evap_no_ifrac, partit)
        call exchange_nod(sublimation, partit)
        call exchange_nod(enthalpyoffuse, partit)
     ! the raw atmospheric state from the ocean's forcing files, in the files'
     ! units: converted as update_atm_forcing converts them
     case ('atm_state_to_ice')
        u_wind(1:nod2)    = exchange(:,1)               ! [m/s]
        v_wind(1:nod2)    = exchange(:,2)               ! [m/s]
        Tair(1:nod2)      = exchange(:,3) - 273.15_WP   ! K -> degC
        shum(1:nod2)      = exchange(:,4)               ! [kg/kg]
        shortwave(1:nod2) = exchange(:,5)               ! [W/m2] downward
        longwave(1:nod2)  = exchange(:,6)               ! [W/m2] downward
        prec_rain(1:nod2) = exchange(:,7)/1000._WP      ! kg/m2/s -> m/s
        prec_snow(1:nod2) = exchange(:,8)/1000._WP      ! kg/m2/s -> m/s
        if (l_mslp) press_air(1:nod2) = exchange(:,9)   ! [Pa]
        call exchange_nod(u_wind, partit)
        call exchange_nod(v_wind, partit)
        call exchange_nod(Tair, partit)
        call exchange_nod(shum, partit)
        call exchange_nod(shortwave, partit)
        call exchange_nod(longwave, partit)
        call exchange_nod(prec_rain, partit)
        call exchange_nod(prec_snow, partit)
        if (l_mslp) call exchange_nod(press_air, partit)
     case default
        write(*,*) 'exchange_ice_oce_yac: no handler for the received field ', trim(cpl_recv(i))
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
     end select
  end do

  if ((do_rotate_oce_wind .AND. do_rotate_ice_wind) .AND. rotated_grid) then
     ! ICON sends the stresses in the geographic frame: one g2r into the
     ! rotated model frame, as the ocean does for its own copy
     do n=1, partit%myDim_nod2D+partit%eDim_nod2D
        call vector_g2r(stress_atmoce_x(n), stress_atmoce_y(n), mesh%coord_nod2D(1, n), mesh%coord_nod2D(2, n), 0)
        call vector_g2r(stress_atmice_x(n), stress_atmice_y(n), mesh%coord_nod2D(1, n), mesh%coord_nod2D(2, n), 0)
     end do
     do_rotate_oce_wind=.false.
     do_rotate_ice_wind=.false.
  end if

  if (.not. is_coupled_to_atmosphere()) then
     ! No atmosphere: the received state becomes drag and heat-exchange
     ! coefficients and wind stress here, as update_atm_forcing does for
     ! forcing files; the heat and freshwater fluxes follow in the bulk ice
     ! thermodynamics.
     if (use_cavity) then
!$OMP PARALLEL DO
        do i=1,partit%myDim_nod2d+partit%eDim_nod2d
           if (mesh%ulevels_nod2d(i)>1) then
              u_wind(i)   = 0.0_WP
              v_wind(i)   = 0.0_WP
              shum(i)     = 0.0_WP
              shortwave(i)= 0.0_WP
              longwave(i) = 0.0_WP
              Tair(i)     = 0.0_WP
              prec_rain(i)= 0.0_WP
              prec_snow(i)= 0.0_WP
              if (l_mslp) press_air(i)= 0.0_WP
           end if
        end do
!$OMP END PARALLEL DO
     endif
     if (AOMIP_drag_coeff) call cal_wind_drag_coeff(partit)
     if (ncar_bulk_formulae) then
        call ncar_ocean_fluxes_mode(ice, partit, mesh)
     elseif (AOMIP_drag_coeff) then
        cd_atm_oce_arr = cd_atm_ice_arr
     end if
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dux, dvy, aux)
     do i=1,partit%myDim_nod2d+partit%eDim_nod2d
        if (mesh%ulevels_nod2d(i)>1) then
           stress_atmoce_x(i)=0.0_WP
           stress_atmoce_y(i)=0.0_WP
           stress_atmice_x(i)=0.0_WP
           stress_atmice_y(i)=0.0_WP
           cycle
        end if
        dux=u_wind(i)-(1.0_WP-Swind)*u_w(i)
        dvy=v_wind(i)-(1.0_WP-Swind)*v_w(i)
        aux=sqrt(dux**2+dvy**2)*rhoair
        if (.not. l_xstre) stress_atmoce_x(i) = Cd_atm_oce_arr(i)*aux*dux
        if (.not. l_ystre) stress_atmoce_y(i) = Cd_atm_oce_arr(i)*aux*dvy
        dux=u_wind(i)-u_ice(i)
        dvy=v_wind(i)-v_ice(i)
        aux=sqrt(dux**2+dvy**2)*rhoair
        stress_atmice_x(i) = Cd_atm_ice_arr(i)*aux*dux
        stress_atmice_y(i) = Cd_atm_ice_arr(i)*aux*dvy
     end do
!$OMP END PARALLEL DO
  end if

end subroutine exchange_ice_oce_yac
#else
subroutine exchange_ice_oce_yac()
end subroutine exchange_ice_oce_yac
#endif /* __cpl_yac */

end module fesim_ocean_exchange_module

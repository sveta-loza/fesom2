module gen_forcing_couple_module
    USE o_PARAM
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_TRACER
    USE MOD_ICE
    USE MOD_DYN
    USE o_arrays
    USE g_forcing_param
    USE g_forcing_arrays
    USE g_clock
    USE g_config
    USE g_comm_auto
    USE g_rotate_grid
    USE g_sbf, only: sbc_do, atmdata, i_totfl, i_xwind, i_ywind, i_xstre, i_ystre, i_humi, &
            i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_cloud, i_snow, l_xwind, l_ywind, l_xstre, &
            l_ystre, l_humi, l_qsr, l_qlw, l_tair, l_prec, l_mslp, l_cloud, l_snow
    USE cpl_yac_driver
    USE gen_bulk
    USE g_support, only: integrate_nod
    USE par_support_module, only: par_abort
    USE cpl_config, only: is_coupled_to_echam, is_coupled_to_oifs
#if defined (__cpl_yac)
    USE cpl_config, only: is_coupled_to_fesim, is_coupled_to_ifs
    USE cpl_fields, only: cpl_send_partner, cpl_recv_partner, &
                          CPL_PARTNER_ATMOSPHERE, CPL_PARTNER_SEA_ICE
#endif
#if defined (__recom)
  use g_sbf, only: sbc_do_recom
#endif
#if defined (__cpl_oasis)
  use cpl_driver
#endif
#if defined (__recom)
  use REcoM_GloVar, only: x_co2atm, GloCO2flux_seaicemask
#endif
#if defined (__cpl_oasis)
  use cpl_driver,	 only : nrecv, cpl_recv, a2o_fcorr_stat
#elif defined (__cpl_yac)
  use cpl_yac_driver,	 only : nrecv, cpl_recv, a2o_fcorr_stat
#endif
#if defined (__cpl_oasis)
  use cpl_driver
#elif defined (__cpl_yac)
  use cpl_yac_driver
#endif
#if defined(__recom) && defined(__usetp)
  use g_config, only: num_fesom_groups
#endif

    implicit none

    private
    public :: update_atm_forcing
#if defined (__cpl_yac)
    public :: update_atm_forcing_yac, exchange_atm_ice_yac, exchange_oce_ice_yac

    ! The atmospheric fluxes as last received from ICON, kept for forwarding
    ! to the sea ice when it runs in FESIM (collection sizes as registered).
    real(kind=WP), allocatable, save :: fwd_taux(:,:), fwd_tauy(:,:)
    real(kind=WP), allocatable, save :: fwd_fresh_water(:,:), fwd_heat_flux(:,:)
    real(kind=WP), allocatable, save :: fwd_atm_sea_ice_bundle(:,:)
#endif
#if defined (__cpl_coupler)
    public :: force_flux_consv, compute_residual, integrate_2D, &
              net_rec_from_atm
#endif

contains

! Routines for updating ocean surface forcing fields
!-------------------------------------------------------------------------
#if defined (__cpl_yac)
!> The atmosphere is ICON over YAC, the sea ice inside the ocean: send the
!> ocean and ice state to the atmosphere, then receive its fluxes.
subroutine update_atm_forcing_yac(istep, ice, tracers, dynamics, partit, mesh)

  implicit none
  integer,        intent(in)            :: istep
  type(t_ice)   , intent(inout), target :: ice
  type(t_tracer), intent(in),    target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  type(t_dyn)   , intent(in), target :: dynamics

  stress_atmoce_x = 0.
  stress_atmoce_y = 0.
  ice%stress_atmice_x = 0.
  ice%stress_atmice_y = 0.
  call send_to_atm_yac(ice, tracers, partit)
  call recv_from_atm_yac(istep, ice, partit, mesh)

end subroutine update_atm_forcing_yac

!-------------------------------------------------------------------------
!> The atmosphere is ICON over YAC and the sea ice runs in FESIM. Called every
!> step: receive the atmosphere's fluxes, exchange with the sea ice every
!> cpl_stride steps (forwarding the fluxes just received), then send the
!> ocean and the received ice state to the atmosphere.
!>
!> The wind stresses are not zeroed before the receive, unlike in
!> update_atm_forcing_yac: a step on which YAC does not couple then keeps the
!> stresses of the last coupling (zero-order hold) instead of applying none.
subroutine exchange_atm_ice_yac(istep, ice, tracers, dynamics, partit, mesh)

  implicit none
  integer,        intent(in)            :: istep
  type(t_ice)   , intent(inout), target :: ice
  type(t_tracer), intent(in),    target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  type(t_dyn)   , intent(in), target :: dynamics

  call recv_from_atm_yac(istep, ice, partit, mesh)
  if (mod(istep-1, cpl_stride) == 0) &
     call exchange_oce_ice_yac(istep, ice, tracers, dynamics, partit, mesh)
  call send_to_atm_yac(ice, tracers, partit)

end subroutine exchange_atm_ice_yac

!-------------------------------------------------------------------------
!> Send the fields registered for the atmosphere.
subroutine send_to_atm_yac(ice, tracers, partit)

  implicit none
  type(t_ice)   , intent(inout), target :: ice
  type(t_tracer), intent(in),    target :: tracers
  type(t_partit), intent(inout), target :: partit
  !_____________________________________________________________________________
  integer                                            :: i
  real(kind=WP), dimension(:,:), allocatable , save  :: exchange
  logical                                            :: action
  real(kind=WP), dimension(:), pointer               :: a_ice, m_ice, m_snow

  a_ice            => ice%data(1)%values(:)
  m_ice            => ice%data(2)%values(:)
  m_snow           => ice%data(3)%values(:)

  if (.NOT. ALLOCATED(exchange)) then
     ALLOCATE(exchange(partit%myDim_nod2D, &
          MAX(MAXVAL(cpl_send_collection_size), MAXVAL(cpl_recv_collection_size))))
     exchange = 0
  end if
  do i=1,nsend
     if (cpl_send_partner(i) /= CPL_PARTNER_ATMOSPHERE) cycle
     exchange  =0.
     select case (trim(cpl_send(i)))
     case ('sst_feom')
        exchange(:,1)=tracers%data(1)%values(1, 1:partit%myDim_nod2d)+273.15 ! sea surface temperature [°K]
     case ('ocean_sea_ice_bundle')
        exchange(:,1) = m_ice(1:partit%myDim_nod2d)              ! ice thickness [m]
        exchange(:,2) = m_snow(1:partit%myDim_nod2d)             ! snow thickness
        exchange(:,3) = a_ice(1:partit%myDim_nod2d)              ! ice concentation [%]
     case default
        write(*,*) 'update_atm_forcing_yac: no handler for the sent field ', trim(cpl_send(i))
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
     end select
     call cpl_yac_send(i, exchange(:,1:cpl_send_collection_size(i)), action)
#ifdef VERBOSE
     if (partit%mype==0) write(*,*) 'SEND: field ', i, ' max val:', maxval(exchange), ' . ACTION? ', action
#endif
  enddo

end subroutine send_to_atm_yac

!-------------------------------------------------------------------------
!> Receive the fields registered from the atmosphere. With the sea ice in
!> FESIM, a copy of each atmospheric flux is kept as received, for
!> exchange_oce_ice_yac to forward.
subroutine recv_from_atm_yac(istep, ice, partit, mesh)

  implicit none
  integer,        intent(in)            :: istep
  type(t_ice)   , intent(inout), target :: ice
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  !_____________________________________________________________________________
  integer                  :: i, n
  real(kind=WP)            :: t1, t2, net
  real(kind=WP), dimension(:,:), allocatable , save  :: exchange
  logical                                          :: action
  logical                                          :: do_rotate_oce_wind=.false.
  logical                                          :: do_rotate_ice_wind=.false.
  !_____________________________________________________________________________
  ! pointer on necessary derived types
  real(kind=WP), dimension(:), pointer  :: stress_atmice_x, stress_atmice_y
  real(kind=WP), dimension(:), pointer  ::  oce_heat_flux, ice_heat_flux
  stress_atmice_x  => ice%stress_atmice_x(:)
  stress_atmice_y  => ice%stress_atmice_y(:)
  oce_heat_flux    => ice%atmcoupl%oce_flx_h(:)
  ice_heat_flux    => ice%atmcoupl%ice_flx_h(:)

  !_____________________________________________________________________________
  t1=MPI_Wtime()
  if (.NOT. ALLOCATED(exchange)) then
     ALLOCATE(exchange(partit%myDim_nod2D, &
          MAX(MAXVAL(cpl_send_collection_size), MAXVAL(cpl_recv_collection_size))))
     exchange = 0
  end if
  do i=1,nrecv
     if (cpl_recv_partner(i) /= CPL_PARTNER_ATMOSPHERE) cycle
     exchange =0.0
     CALL cpl_yac_recv (i, exchange(:,1:cpl_recv_collection_size(i)), action)
#ifdef VERBOSE
     if (partit%mype==0) then
        write(*,*) 'FESOM RECV: flux ', i, ', max val: ', maxval(exchange), ' . ACTION? ', action
     end if
#endif
     if (.not. action) cycle
     if (is_coupled_to_fesim) call keep_for_sea_ice(cpl_recv(i), exchange(:,1:cpl_recv_collection_size(i)))
     !Do not apply a correction at first time step!
     select case (trim(cpl_recv(i)))
     case ('taux')
        stress_atmoce_x(1:partit%myDim_nod2d) =  exchange(:,1)                    ! taux_oce
        call exchange_nod(stress_atmoce_x, partit)
        stress_atmice_x(1:partit%myDim_nod2d) =  exchange(:,2)                    ! taux_ice
        call exchange_nod(stress_atmice_x, partit)
        do_rotate_oce_wind=.true.
        do_rotate_ice_wind=.true.
     case ('tauy')
        stress_atmoce_y(1:partit%myDim_nod2d) =  exchange(:,1)                    ! tauy_oce
        call exchange_nod(stress_atmoce_y, partit)
        stress_atmice_y(1:partit%myDim_nod2d) =  exchange(:,2)                    ! tauy_ice
        call exchange_nod(stress_atmice_y, partit)
        do_rotate_oce_wind=.true.
        do_rotate_ice_wind=.true.
     case ('surface_fresh_water_flux')
        prec_rain(1:partit%myDim_nod2d)    =  exchange(:,1)/1000                   ! tot_prec ! kg m^(-2) s^(-1) -> m/s
        call exchange_nod(prec_rain, partit)
        prec_snow(1:partit%myDim_nod2d)    =  exchange(:,2)/1000                    ! snowfall ! kg m^(-2) s^(-1) -> m/s
        call exchange_nod(prec_snow, partit)
        evap_no_ifrac(1:partit%myDim_nod2d)     =  exchange(:,3)/1000               ! tot_evap ! kg m^(-2) s^(-1) -> m/s; change sign
        call exchange_nod(evap_no_ifrac, partit)
     case ('total_heat_flux')
        oce_heat_flux(1:partit%myDim_nod2d)     = exchange(:,2) + exchange(:,3) + exchange(:,4) ! heat_oce
        call exchange_nod(oce_heat_flux, partit)
        shortwave(1:partit%myDim_nod2d)         =  exchange(:,1)                ! heat_swr
        call exchange_nod(shortwave, partit)
     case ('atmosphere_sea_ice_bundle')
        ice_heat_flux(1:partit%myDim_nod2d)     = (exchange(:,2) + exchange(:,1)) ! heat_ice
        call exchange_nod(ice_heat_flux, partit)
!SL-- add river runoff ------------------------------
     case ('river_runoff')
        runoff(1:partit%myDim_nod2D) = exchange(:,1) * mesh%area_inv(1,1:partit%myDim_nod2D)
        call exchange_nod(runoff, partit)
        call integrate_nod(runoff, net, partit, mesh)
        if(partit%mype==0) write(*,*) 'RUNOFF CHECK:', net
!SL--------------------------------------------------
     case default
        write(*,*) 'update_atm_forcing_yac: no handler for the received field ', trim(cpl_recv(i))
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
     end select
  end do

  if ((do_rotate_oce_wind .AND. do_rotate_ice_wind) .AND. rotated_grid) then
     do n=1, partit%myDim_nod2D+partit%eDim_nod2D
        call vector_g2r(stress_atmoce_x(n), stress_atmoce_y(n), mesh%coord_nod2D(1, n), mesh%coord_nod2D(2, n), 0)
        call vector_g2r(stress_atmice_x(n), stress_atmice_y(n), mesh%coord_nod2D(1, n), mesh%coord_nod2D(2, n), 0)
     end do
     do_rotate_oce_wind=.false.
     do_rotate_ice_wind=.false.
  end if

  t2=MPI_Wtime()

#ifdef VERBOSE
  if (mod(istep,logfile_outfreq)==0 .and. partit%mype==0) then
     write(*,*) 'update forcing data took', t2-t1
  end if
#endif

end subroutine recv_from_atm_yac
#endif /* __cpl_yac */

! update_atm_forcing is compiled in every build. Under YAC it is the path
! taken when the ocean's only partner is the sea ice (FESOM + FESIM): the
! ocean reads its forcing files and runs its own bulk formulae.

subroutine update_atm_forcing(istep, ice, tracers, dynamics, partit, mesh)
  use fesom_profiler, only: fesom_profiler_start, fesom_profiler_end

  implicit none
  integer,        intent(in)            :: istep
  type(t_ice)   , intent(inout), target :: ice
  type(t_tracer), intent(in),    target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  type(t_dyn)   , intent(in), target :: dynamics
  !_____________________________________________________________________________
  integer		   :: i, itime,n2,n,nz,k,elem
  real(kind=WP)            :: i_coef, aux
  real(kind=WP)	           :: dux, dvy,tx,ty,tvol
  real(kind=WP)            :: t1, t2
  !---wiso-code
  integer                  :: nt1, nt2
  real(kind=WP), parameter :: zwisomin = 1.e-6_WP
  !---wiso-code-end
#if defined (__cpl_oasis)
  real(kind=WP)        				   :: flux_global(2), flux_local(2), eff_vol(2)
  real(kind=WP), dimension(:), allocatable , save  :: exchange
  real(kind=WP), dimension(:), allocatable , save  :: mask !, weight
  logical, save                                    :: firstcall=.true.
  logical                                          :: action
  logical                                          :: do_rotate_oce_wind=.false.
  logical                                          :: do_rotate_ice_wind=.false.
  INTEGER                                          :: my_global_rank, ierror
  INTEGER 					   :: status(MPI_STATUS_SIZE)
#endif
  !character(15)                         :: vari, filevari
  !character(4)                          :: fileyear
  !integer, parameter                    :: nci=192, ncj=94 ! T62 grid
  !real(kind=WP), dimension(nci,ncj)     :: array_nc, array_nc2,array_nc3,x
  !character(500)                        :: file
  !_____________________________________________________________________________
  ! pointer on necessary derived types
  real(kind=WP), dimension(:), pointer  :: u_ice, v_ice, u_w, v_w
  real(kind=WP), dimension(:), pointer  :: stress_atmice_x, stress_atmice_y
  real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
#if defined (__cpl_enabled)
  real(kind=WP), dimension(:), pointer  ::  oce_heat_flux, ice_heat_flux 
  real(kind=WP), dimension(:), pointer  ::  tmp_oce_heat_flux, tmp_ice_heat_flux 
#endif 
  real(kind=WP), dimension(:), pointer  :: ice_temp, ice_alb, enthalpyoffuse, runoff_liquid, runoff_solid
  real(kind=WP),               pointer  :: tmelt
  real(kind=WP), dimension(:,:,:), pointer :: UVnode
  logical                               :: l_ist
  real(kind=WP)              , pointer  :: rhoair
  u_ice            => ice%uice(:)
  v_ice            => ice%vice(:)
  u_w              => ice%srfoce_u(:)
  v_w              => ice%srfoce_v(:)
  stress_atmice_x  => ice%stress_atmice_x(:)
  stress_atmice_y  => ice%stress_atmice_y(:)
  a_ice            => ice%data(1)%values(:)
  m_ice            => ice%data(2)%values(:)
  m_snow           => ice%data(3)%values(:)
  tmelt            => ice%thermo%tmelt
  UVnode           => dynamics%uvnode(:,:,:)
  ! The ice surface temperature is only a tracer for the IFS-family partners.
  l_ist = ice%ist_itracer_idx > 0
  if (l_ist) ice_temp => ice%data(ice%ist_itracer_idx)%values(:)
#if defined (__cpl_enabled)
  ice_alb          => ice%atmcoupl%ice_alb(:)
  enthalpyoffuse   => ice%atmcoupl%enthalpyoffuse(:)
  runoff_liquid    => ice%atmcoupl%runoff_liquid(:)
  runoff_solid     => ice%atmcoupl%runoff_solid(:)
#endif  
#if defined (__cpl_enabled)
  oce_heat_flux    => ice%atmcoupl%oce_flx_h(:)
  ice_heat_flux    => ice%atmcoupl%ice_flx_h(:)
  tmp_oce_heat_flux=> ice%atmcoupl%tmpoce_flx_h(:)
  tmp_ice_heat_flux=> ice%atmcoupl%tmpice_flx_h(:)
#endif 
  rhoair           => ice%thermo%rhoair
  
  !_____________________________________________________________________________
  t1=MPI_Wtime()
#if defined (__cpl_oasis)
     if (firstcall) then
        allocate(exchange(partit%myDim_nod2D+partit%eDim_nod2D), mask(partit%myDim_nod2D+partit%eDim_nod2D))
        allocate(a2o_fcorr_stat(nrecv,6))
        a2o_fcorr_stat=0.
        exchange      =0.
        mask          =0.
        firstcall=.false.
     end if
     do i=1,nsend
         exchange  =0.
         ! The sent field set is the one define_cpl_fields registered for
         ! the partner; each field is filled by its name.
         select case (trim(cpl_send(i)))
         case ('sst_feom')
            ! sea surface temperature: OpenIFS takes [K], ECHAM6 [°C]
            if (is_coupled_to_oifs) then
               do n=1,partit%myDim_nod2D+partit%eDim_nod2D
                  exchange(n)=tracers%data(1)%values(1, n)+tmelt
               end do
            else
               do n=1,partit%myDim_nod2D+partit%eDim_nod2D
                  exchange(n)=tracers%data(1)%values(1, n)
               end do
            end if
         case ('sie_feom')
            exchange(:) = a_ice(:)                                  ! ice concentation [%]
         case ('snt_feom')
            exchange(:) = m_snow(:)                                 ! snow thickness
         case ('sit_feom')
            ! Effective (grid-mean = per-ice * concentration) sea-ice
            ! thickness. OpenIFS uses it for its ice-tile slab conduction
            ! (ECE_FESIM_GET_ICE_STATE divides by the received ice fraction
            ! under LNEMOLIMTHK). Weighted convention as ist/alb.
            exchange(:) = m_ice(:)                                  ! effective sea ice thickness
         case ('ist_feom')
            ! Concentration-weighted ice surface temperature (ist*a_ice).
            ! The o2a remap blend must be concentration-weighted so the
            ! atmosphere's prescribed ice-tile skin is representative of the
            ! ice actually present in the cell; PAIRED with
            ! ECE_CPL_NEMO_WEIGHTED_ICE=.true. (OIFS divides by the received
            ! ice fraction on ingest) and a weighted-convention rstos.nc.
            exchange(:) = ice_temp(:)*a_ice(:)                      ! ice surface temperature * concentration
         case ('sia_feom')
            exchange(:) = ice_alb(:)                                ! ice albedo
         case ('u_feom', 'v_feom')
            ! surface velocity, sent as eastward and northward component:
            ! rotate from the model grid to geographical coordinates
            do n=1,partit%myDim_nod2D+partit%eDim_nod2D
              tx = UVnode(1,1,n)
              ty = UVnode(2,1,n)
              if (rotated_grid) call vector_r2g(tx, ty, mesh%coord_nod2D(1, n), mesh%coord_nod2D(2, n), 0)
              if (trim(cpl_send(i)) == 'u_feom') then
                exchange(n) = tx
              else
                exchange(n) = ty
              end if
            end do
#if defined (__recom)
         case ('FCO2_feom')
            ! GloCO2flux_seaicemask is in [mmolCO2 m-2 s-1], need [kgCO2 m-2 s-1]
            ! Conversion: 1.0e-3_WP -> mol/s -> kg/s
            ! 1 mol CO2 = 44.0095 g/mol = 0.0440095 kg/mol (NIST 2018)
            ! *-1 for correct flux direction convention: oifs expects >0: downward; fesom: >0: upward
            exchange(:) = GloCO2flux_seaicemask(:) * 1.0e-3_WP * 0.0440095_WP * -1 ! [kgCO2 m-2 s-1]
#endif
!---wiso-code
         case ('o18w_oce')
            nt1 = index_wiso_tracers(1)
            nt2 = index_wiso_tracers(3)
            where (tracers%data(nt2)%values(1,:) > zwisomin)
              exchange(:) = (tracers%data(nt1)%values(1,:)/tracers%data(nt2)%values(1,:)/wiso_smow(1) - 1._WP)*1000._WP ! delta 18O of surface water
            end where
         case ('hdow_oce')
            nt1 = index_wiso_tracers(2)
            nt2 = index_wiso_tracers(3)
            where (tracers%data(nt2)%values(1,:) > zwisomin)
              exchange(:) = (tracers%data(nt1)%values(1,:)/tracers%data(nt2)%values(1,:)/wiso_smow(2) - 1._WP)*1000._WP ! delta D of surface water
            end where
         case ('o16w_oce')
            exchange(:) = 0.0_WP                          ! delta H216O of surface water is set to zero permill
         case ('o18i_oce')
            where (tr_arr_ice(:, 3) > zwisomin)
              exchange(:) = (tr_arr_ice(:, 1)/tr_arr_ice(:, 3)/wiso_smow(1) - 1._WP)*1000._WP ! delta 18O of sea ice
            end where
         case ('hdoi_oce')
            where (tr_arr_ice(:, 3) > zwisomin)
              exchange(:) = (tr_arr_ice(:, 2)/tr_arr_ice(:, 3)/wiso_smow(2) - 1._WP)*1000._WP ! delta 18O of sea ice
            end where
         case ('o16i_oce')
            exchange(:) = 0.0_WP                          ! delta H216O of sea ice is set to zero permill
!---wiso-code-end
         case default
            write(*,*) 'update_atm_forcing: no handler for the sent field ', trim(cpl_send(i))
            call par_abort(partit%MPI_COMM_FESOM, partit%mype)
         end select

#if defined(__recom) && defined(__usetp)
         if(partit%my_fesom_group == 0) then
#endif
         call cpl_oasis3mct_send(i, exchange, action, partit)
         ! Anchor for the implicit ice surface-temperature solve
         ! (ice_thermo_cpl.F90/ice_surftemp): remember the ist as ACTUALLY
         ! transmitted -- the temperature OIFS evaluates its ice-tile fluxes
         ! at for the coming coupling interval. `action` is only true on real
         ! OASIS transmissions, so this stays frozen between coupling events.
         if (trim(cpl_send(i)) == 'ist_feom' .and. action) &
            ice%atmcoupl%ist_ref(:) = exchange(:)
#if defined(__recom) && defined(__usetp)
         endif
#endif
      end do
#ifdef VERBOSE
      do i=1, nsend 
        if (partit%mype==0) write(*,*) 'SEND: field ', i, ' max val:', maxval(exchange), ' . ACTION? ', action 
      end do
#endif
      mask=1.
      do i=1,nrecv
         exchange =0.0
         call cpl_oasis3mct_recv (i, exchange, action, partit)
	 !if (.not. action) cycle
	 !Do not apply a correction at first time step!
         if (i==1 .and. action .and. istep/=1) call net_rec_from_atm(action, partit)
         ! The received field set is the one define_cpl_fields registered
         ! for the partner; each field is taken by its name.
         select case (trim(cpl_recv(i)))
         case ('taux_oce')
            if (.not. action) cycle
            stress_atmoce_x(:) =  exchange(:)                    ! taux_oce
            do_rotate_oce_wind=.true.
         case ('tauy_oce')
            if (.not. action) cycle
            stress_atmoce_y(:) =  exchange(:)                    ! tauy_oce
            do_rotate_oce_wind=.true.
         case ('taux_ico')
            if (.not. action) cycle
            stress_atmice_x(:) =  exchange(:)                    ! taux_ice
            do_rotate_ice_wind=.true.
         case ('tauy_ico')
            if (.not. action) cycle
            stress_atmice_y(:) =  exchange(:)                    ! tauy_ice
            do_rotate_ice_wind=.true.
         case ('prec_oce')
            if (action) then
                prec_rain(:)    =  exchange(:)	                  ! tot_prec
                mask=1.
                call force_flux_consv(prec_rain, mask, i, 0,action, partit, mesh)
            end if
         case ('snow_oce')
            if (action) then
                prec_snow(:)    =  exchange(:)                    ! snowfall
                mask=1.
                call force_flux_consv(prec_snow, mask,i,1,action, partit, mesh) ! Northern hemisphere
                call force_flux_consv(prec_snow, mask,i,2,action, partit, mesh) ! Southern Hemisphere
            end if
         case ('evap_oce')
            if (action) then
                evap_no_ifrac(:)     =  exchange(:)        	  ! tot_evap
                tmp_evap_no_ifrac(:) =  exchange(:) 		  ! to reset for flux
                                                              ! correction
            end if
            mask=1.-a_ice
            evap_no_ifrac(:)     =  tmp_evap_no_ifrac(:)
            call force_flux_consv(evap_no_ifrac,mask,i,0,action, partit, mesh)
         case ('subl_oce')
            if (action) then
                sublimation(:)       =  exchange(:)        	  ! tot_subl
                tmp_sublimation(:)   =  exchange(:) 		  ! to reset for flux
                                                              ! correction
            end if
            mask=a_ice
            sublimation(:)       =  tmp_sublimation(:)
            call force_flux_consv(sublimation,mask,i,1,action, partit, mesh) ! Northern hemisphere
            call force_flux_consv(sublimation,mask,i,2,action, partit, mesh) ! Southern Hemisphere
         case ('heat_oce')
            if (action) then
                oce_heat_flux(:)     =  exchange(:)        	  ! heat_oce
                tmp_oce_heat_flux(:) =  exchange(:) 		  ! to reset for flux
                                                              ! correction
            end if
            mask=1.-a_ice
            oce_heat_flux(:)     =  tmp_oce_heat_flux(:)
            call force_flux_consv(oce_heat_flux, mask, i, 0,action, partit, mesh)
         case ('heat_ico')
            if (action) then
                ice_heat_flux(:)     =  exchange(:)        	  ! heat_ice
                tmp_ice_heat_flux(:) =  exchange(:) 		  ! to reset for flux
                                                              ! correction
            end if
            mask=a_ice
            ice_heat_flux(:)     =  tmp_ice_heat_flux(:)
            call force_flux_consv(ice_heat_flux, mask, i, 1,action, partit, mesh) ! Northern hemisphere
            call force_flux_consv(ice_heat_flux, mask, i, 2,action, partit, mesh) ! Southern Hemisphere
         case ('heat_swo')
            if (action) then
                shortwave(:)         =  exchange(:)		        ! heat_swr
                tmp_shortwave(:)     =  exchange(:) 		    ! to reset for flux
                                                                ! correction
            end if
            mask=1.-a_ice
            shortwave(:)   =  tmp_shortwave(:)
            call force_flux_consv(shortwave, mask, i, 0,action, partit, mesh)
         case ('hydr_oce')
            if (action) then
                runoff(:)            =  exchange(:)        ! AWI-CM2: runoff, AWI-CM3: runoff + excess snow on glaciers
                mask=1.
                call force_flux_consv(runoff, mask, i, 0,action, partit, mesh)
            end if
         case ('calv_oce')
            if (action) then
                runoff_solid(:)      =  exchange(:)        ! solid water discharge (excess snow --> calving)
                enthalpyoffuse(:)            =  runoff_solid(:)*333.55*1000000.0        ! enthalpy of fusion of the solid water discharge from glaciers
                enthalpyoffuse = -min(enthalpyoffuse, 1000.0)
                runoff_liquid(:) = runoff(:)                     ! store previous liquid runoff
                runoff(:)            = runoff_liquid(:) + runoff_solid(:)                      ! Add solid runoff to the liquid runoff. Heatflux goes into enthalpyoffuse.

                mask=1.
                call force_flux_consv(runoff_liquid, mask, i, 0, action, partit, mesh)
                call force_flux_consv(runoff_solid, mask, i, 0, action, partit, mesh)
                call force_flux_consv(enthalpyoffuse, mask, i, 0, action, partit, mesh)
            end if
         case ('u10w_oce')
            if (action) then
                u_wind(:)                     = exchange(:)        ! zonal wind
            end if
            ! ECHAM6 sends the winds for the icebergs and corrects them like a flux
            if (is_coupled_to_echam) then
                mask=1.
                call force_flux_consv(u_wind, mask, i, 0, action, partit, mesh)
            end if
         case ('v10w_oce')
            if (action) then
                v_wind(:)                     = exchange(:)        ! meridional wind
            end if
            if (is_coupled_to_echam) then
                mask=1.
                call force_flux_consv(v_wind, mask, i, 0, action, partit, mesh)
            end if
         case ('tsk_ico')
            if (action) then
                ice%atmcoupl%ist_ref(:)       = exchange(:)        ! ice-tile skin temperature
            end if
         case ('mslp_oce')
            if (action) then
                ! Mean sea-level pressure [Pa]. Only the gas exchange reads it here:
                ! the pressure term of the momentum equation stays under l_mslp.
                press_air(:)                  = exchange(:)
            end if
#if defined (__recom)
         case ('XCO2_oce')
            if (action) then
                ! Convert mass mixing ratio (kg/kg) to ppm
                ! MW_CO2 = 44.0095 g/mol (NIST 2018)
                ! MW_dry_air = 28.9647 g/mol (standard atmosphere composition)
                x_co2atm(:) = exchange(:) * ((28.9647_WP/44.0095_WP)*1e6_WP)  ! [ppm]
            end if
#endif
!---wiso-code
         case ('w1_oce')
            if (action) then
                www3(:)         =  exchange(:)
            end if
            mask=1.
            call force_flux_consv(www3, mask, i, 0,action, partit, mesh)
         case ('w2_oce')
            if (action) then
                www1(:)         =  exchange(:)
            end if
            mask=1.
            call force_flux_consv(www1, mask, i, 0,action, partit, mesh)
         case ('w3_oce')
            if (action) then
            ! tot_prec_hdo over water: this variable includes (i) rain over open water and sea ice, (ii) snow and evap over open water,  (iii) river runoff
            www2(:)         =  exchange(:)
            end if
            mask=1.
            call force_flux_consv(www2, mask, i, 0,action, partit, mesh)
         case ('i1_oce')
            if (action) then
            ! snowfall_o16 over seaice: this variable includes snow and sublimation over sea ice areas
            iii3(:)         =  exchange(:)
            tmp_iii3(:)     =  exchange(:)                   ! to reset for flux correction
            end if
            mask=a_ice
            iii3(:)         =  tmp_iii3(:)
            call force_flux_consv(iii3,mask,i,1,action, partit, mesh) ! Northern hemisphere
            call force_flux_consv(iii3,mask,i,2,action, partit, mesh) ! Southern Hemisphere
         case ('i2_oce')
            if (action) then
            ! snowfall_o18 over seaice: this variable includes snow and sublimation over sea ice areas
            iii1(:)         =  exchange(:)
            tmp_iii1(:)     =  exchange(:)                   ! to reset for flux correction
            end if
            mask=a_ice
            iii1(:)         =  tmp_iii1(:)
            call force_flux_consv(iii1,mask,i,1,action, partit, mesh) ! Northern hemisphere
            call force_flux_consv(iii1,mask,i,2,action, partit, mesh) ! Southern Hemisphere
         case ('i3_oce')
            if (action) then
            ! snowfall_hdo over seaice: this variable includes snow and sublimation over sea ice areas
            iii2(:)         =  exchange(:)
            tmp_iii2(:)     =  exchange(:)                   ! to reset for flux correction
            end if
            mask=a_ice
            iii2(:)         =  tmp_iii2(:)
            call force_flux_consv(iii2,mask,i,1,action, partit, mesh) ! Northern hemisphere
            call force_flux_consv(iii2,mask,i,2,action, partit, mesh) ! Southern Hemisphere
!---wiso-code-end
         case default
            write(*,*) 'update_atm_forcing: no handler for the received field ', trim(cpl_recv(i))
            call par_abort(partit%MPI_COMM_FESOM, partit%mype)
         end select

#ifdef VERBOSE
	  if (partit%mype==0) then
		write(*,*) 'FESOM RECV: flux ', i, ', max val: ', maxval(exchange)
	  end if
#endif
      end do

    if ((do_rotate_oce_wind .AND. do_rotate_ice_wind) .AND. rotated_grid) then
        do n=1, partit%myDim_nod2D+partit%eDim_nod2D
            call vector_g2r(stress_atmoce_x(n), stress_atmoce_y(n), mesh%coord_nod2D(1, n), mesh%coord_nod2D(2, n), 0)
            call vector_g2r(stress_atmice_x(n), stress_atmice_y(n), mesh%coord_nod2D(1, n), mesh%coord_nod2D(2, n), 0)
        end do
        do_rotate_oce_wind=.false.
        do_rotate_ice_wind=.false.
    end if
#else
#ifndef __cpl_direct
  call fesom_profiler_start("forc_sbc_do")
  call sbc_do(partit, mesh)
  call fesom_profiler_end("forc_sbc_do")
!$OMP PARALLEL DO
  DO n=1, partit%myDim_nod2D+partit%eDim_nod2D
     u_wind(n)    = atmdata(i_xwind,n)
     v_wind(n)    = atmdata(i_ywind,n)
     if (l_xstre) then
        stress_atmoce_x(n)=atmdata(i_xstre, n)
     end if
     if (l_ystre) then
        stress_atmoce_y(n)=atmdata(i_ystre, n)
     end if
     shum(n)      = atmdata(i_humi ,n)
     longwave(n)  = atmdata(i_qlw  ,n)
     shortwave(n) = atmdata(i_qsr  ,n)
     Tair(n)      = atmdata(i_tair ,n)-273.15_WP
     prec_rain(n) = atmdata(i_prec ,n)/1000._WP
     prec_snow(n) = atmdata(i_snow ,n)/1000._WP
     if (l_mslp) then
        press_air(n) = atmdata(i_mslp ,n) ! unit should be Pa
     end if
  END DO
!$OMP END PARALLEL DO

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
            ! press_air is kept: a floating ice shelf passes the atmospheric pressure on to
            ! the water below, and zeroing it here made a pressure step at the ice-shelf front
            runoff(i)   = 0.0_WP
        end if
    end do
!$OMP END PARALLEL DO
  endif 
  ! second, compute exchange coefficients
  ! 1) drag coefficient 
  if(AOMIP_drag_coeff) then
     call fesom_profiler_start("forc_wind_drag")
     call cal_wind_drag_coeff(partit)
     call fesom_profiler_end("forc_wind_drag")
  end if
  ! 2) drag coeff. and heat exchange coeff. over ocean in case using ncar formulae
  if(ncar_bulk_formulae) then
!     cd_atm_oce_arr=0.0_WP
!     ch_atm_oce_arr=0.0_WP
!     ce_atm_oce_arr=0.0_WP
     call fesom_profiler_start("forc_ncar_fluxes")
     call ncar_ocean_fluxes_mode(ice, partit, mesh)
     call fesom_profiler_end("forc_ncar_fluxes")
  elseif(AOMIP_drag_coeff) then
     cd_atm_oce_arr=cd_atm_ice_arr
  end if
  ! third, compute wind stress
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dux, dvy, aux)

  do i=1,partit%myDim_nod2d+partit%eDim_nod2d   
     !__________________________________________________________________________
     if (mesh%ulevels_nod2d(i)>1) then
        stress_atmoce_x(i)=0.0_WP
        stress_atmoce_y(i)=0.0_WP
        stress_atmice_x(i)=0.0_WP
        stress_atmice_y(i)=0.0_WP
        cycle
     end if 
     
     !__________________________________________________________________________
     dux=u_wind(i)-(1.0_WP-Swind)*u_w(i) 
     dvy=v_wind(i)-(1.0_WP-Swind)*v_w(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     if ( .not. l_xstre) then
        stress_atmoce_x(i) = Cd_atm_oce_arr(i)*aux*dux
     end if
     if ( .not. l_ystre) then
        stress_atmoce_y(i) = Cd_atm_oce_arr(i)*aux*dvy
     end if
     !__________________________________________________________________________
     dux=u_wind(i)-u_ice(i) 
     dvy=v_wind(i)-v_ice(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     if (use_ice) then
        stress_atmice_x(i) = Cd_atm_ice_arr(i)*aux*dux
        stress_atmice_y(i) = Cd_atm_ice_arr(i)*aux*dvy
     end if
  end do
!$OMP END PARALLEL DO
  ! heat and fresh water fluxes are treated in i_therm and ice2ocean
#endif /* skip all in case of __cpl_direct */
#endif /* (__cpl_oasis) */

#if defined (__recom) /* consider in all cases */
  call sbc_do_recom(partit, mesh)
#endif

  t2=MPI_Wtime()

#ifdef VERBOSE
  if (mod(istep,logfile_outfreq)==0 .and. partit%mype==0) then
     write(*,*) 'update forcing data took', t2-t1
  end if
#endif

end subroutine update_atm_forcing

#if defined (__cpl_yac)
!-------------------------------------------------------------------------
!> Keep a copy of an atmospheric flux as received from ICON, for forwarding
!> to the sea ice (FESIM converts on receive).
subroutine keep_for_sea_ice(name, data)
  implicit none
  character(len=*), intent(in) :: name
  real(kind=WP),    intent(in) :: data(:,:)

  if (.not. allocated(fwd_taux)) call allocate_forwarded(size(data, 1))
  select case (trim(name))
  case ('taux')
     fwd_taux = data
  case ('tauy')
     fwd_tauy = data
  case ('surface_fresh_water_flux')
     fwd_fresh_water = data
  case ('total_heat_flux')
     fwd_heat_flux = data
  case ('atmosphere_sea_ice_bundle')
     fwd_atm_sea_ice_bundle = data
  end select   ! river_runoff is not forwarded as such: the ocean sends its runoff
end subroutine keep_for_sea_ice

subroutine allocate_forwarded(n)
  implicit none
  integer, intent(in) :: n
  allocate(fwd_taux(n, 2), fwd_tauy(n, 2), fwd_fresh_water(n, 3), &
           fwd_heat_flux(n, 4), fwd_atm_sea_ice_bundle(n, 2))
  fwd_taux               = 0.0_WP
  fwd_tauy               = 0.0_WP
  fwd_fresh_water        = 0.0_WP
  fwd_heat_flux          = 0.0_WP
  fwd_atm_sea_ice_bundle = 0.0_WP
end subroutine allocate_forwarded

!-------------------------------------------------------------------------
!> Ocean <-> sea-ice (FESIM) exchange over YAC, every cpl_stride steps. The
!> ocean receives the ice state and the ice->ocean fluxes, and sends its
!> surface state plus what it forwards from the atmosphere:
!>  - ICON: the fluxes as last received (exchange_atm_ice_yac);
!>  - IFS: the fluxes the IFS interface deposited in the ocean's arrays this
!>    step, already in FESOM units and rotated (stress_atmoce stays here: it
!>    drives the ocean). The ice surface temperature and albedo received back
!>    are relayed to IFS by the interface. The ocean also adds the open-water
!>    part of the atmospheric flux of this step itself (see below);
!>  - no atmosphere: the raw atmospheric state from the forcing files, which
!>    update_atm_forcing has just read.
!> The fields are the ones define_cpl_fields registered for the sea ice.
subroutine exchange_oce_ice_yac(istep, ice, tracers, dynamics, partit, mesh)
  implicit none
  integer,        intent(in)            :: istep
  type(t_ice)   , intent(inout), target :: ice
  type(t_tracer), intent(in),    target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
  type(t_dyn)   , intent(in),    target :: dynamics
  integer :: i, n, nod2
  integer :: k_drift   ! first component of the ice drift in sea_ice_bundle
  logical :: action
  real(kind=WP), dimension(:,:), allocatable, save :: exchange
  real(kind=WP), dimension(:),   allocatable, save :: flx_h_ice, flx_fw_ice   ! ice-related part of the received flux (IFS)
  real(kind=WP), dimension(:), pointer :: a_ice, m_ice, m_snow, ice_temp, ice_alb
  real(kind=WP), dimension(:), pointer :: oce_heat_flux, ice_heat_flux, enthalpyoffuse

  nod2   =  partit%myDim_nod2D
  a_ice  => ice%data(1)%values(:)
  m_ice  => ice%data(2)%values(:)
  m_snow => ice%data(3)%values(:)
  oce_heat_flux  => ice%atmcoupl%oce_flx_h(:)
  ice_heat_flux  => ice%atmcoupl%ice_flx_h(:)
#if defined (__cpl_enabled)
  enthalpyoffuse => ice%atmcoupl%enthalpyoffuse(:)
#endif
  if (is_coupled_to_ifs) then
     ! the ice surface tracer exists for the partners FESOM provides it to
     ice_temp => ice%data(ice%ist_itracer_idx)%values(:)
     ice_alb  => ice%atmcoupl%ice_alb(:)
  end if

  if (.NOT. ALLOCATED(exchange)) then
     ALLOCATE(exchange(nod2, MAX(MAXVAL(cpl_send_collection_size), MAXVAL(cpl_recv_collection_size))))
     exchange = 0
     ALLOCATE(flx_h_ice(nod2+partit%eDim_nod2D), flx_fw_ice(nod2+partit%eDim_nod2D))
     flx_h_ice  = 0.0_WP
     flx_fw_ice = 0.0_WP
  end if
  if (.not. allocated(fwd_taux)) call allocate_forwarded(nod2)

  ! ---- receive from the sea ice ------------------------------------------
  do i=1,nrecv
     if (cpl_recv_partner(i) /= CPL_PARTNER_SEA_ICE) cycle
     exchange = 0.0
     CALL cpl_yac_recv(i, exchange(:,1:cpl_recv_collection_size(i)), action)
     if (.not. action) cycle
     select case (trim(cpl_recv(i)))
     case ('sea_ice_bundle')
        m_ice(1:nod2)  = exchange(:,1)   ! ice thickness [m]
        m_snow(1:nod2) = exchange(:,2)   ! snow thickness [m]
        a_ice(1:nod2)  = exchange(:,3)   ! ice concentration
        call exchange_nod(m_ice, partit)
        call exchange_nod(m_snow, partit)
        call exchange_nod(a_ice, partit)
        k_drift = 4
        if (is_coupled_to_ifs) then
           ice_temp(1:nod2) = exchange(:,4)   ! ice surface temperature [K], relayed to IFS
           ice_alb(1:nod2)  = exchange(:,5)   ! ice albedo, relayed to IFS
           call exchange_nod(ice_temp, partit)
           call exchange_nod(ice_alb, partit)
           k_drift = 6
        end if
        ! ice drift: the ocean's under-ice mixing (TB04, oce_mo_conv) uses the ice speed
        ice%uice(1:nod2) = exchange(:,k_drift)     ! [m/s]
        ice%vice(1:nod2) = exchange(:,k_drift+1)   ! [m/s]
        call exchange_nod(ice%uice, partit)
        call exchange_nod(ice%vice, partit)
     case ('ice_to_ocean_stress')
        ! ice-ocean drag; oce_fluxes_mom uses it as received
        ice%stress_iceoce_x(1:nod2) = exchange(:,1)  ! [Pa]
        ice%stress_iceoce_y(1:nod2) = exchange(:,2)  ! [Pa]
        call exchange_nod(ice%stress_iceoce_x, partit)
        call exchange_nod(ice%stress_iceoce_y, partit)
     case ('ice_to_ocean_flux')
        ! net heat and freshwater flux; oce_fluxes applies heat_flux=-flx_h,
        ! water_flux=-flx_fw, as with the built-in sea ice
        ice%flx_h(1:nod2)  = exchange(:,1)  ! [W/m2]
        ice%flx_fw(1:nod2) = exchange(:,2)  ! [m/s]
        call exchange_nod(ice%flx_h, partit)
        call exchange_nod(ice%flx_fw, partit)
        flx_h_ice  = ice%flx_h
        flx_fw_ice = ice%flx_fw
     case ('ice_to_ocean_thermo')
        ! the terms of the ice thermodynamics that oce_fluxes integrates for the
        ! global freshwater balancing; the built-in sea ice sets them itself
        evaporation(1:nod2)            = exchange(:,1)  ! [m/s], negative up
        ice_sublimation(1:nod2)        = exchange(:,2)  ! [m/s]
        ice%thermo%thdgr(1:nod2)       = exchange(:,3)  ! ice growth rate  [m/s]
        ice%thermo%thdgrsn(1:nod2)     = exchange(:,4)  ! snow growth rate [m/s]
        ice%data(1)%values_old(1:nod2) = exchange(:,5)  ! a_ice before the thermodynamics
        real_salt_flux(1:nod2)         = exchange(:,6)  ! brine salt flux [m/s psu]
        call exchange_nod(evaporation, partit)
        call exchange_nod(ice_sublimation, partit)
        call exchange_nod(ice%thermo%thdgr, partit)
        call exchange_nod(ice%thermo%thdgrsn, partit)
        call exchange_nod(ice%data(1)%values_old, partit)
        call exchange_nod(real_salt_flux, partit)
     case default
        write(*,*) 'exchange_oce_ice_yac: no handler for the received field ', trim(cpl_recv(i))
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
     end select
  end do

  ! ---- send to the sea ice ------------------------------------------------
  do i=1,nsend
     if (cpl_send_partner(i) /= CPL_PARTNER_SEA_ICE) cycle
     exchange = 0.
     select case (trim(cpl_send(i)))
     case ('sst_feom_to_ice')
        ! the ocean surface state as ocean2ice prepares it for the built-in sea
        ! ice, in degC (a Kelvin round trip leaves 1e-14 K that flips the
        ! new-ice branch where the ocean sits at the freezing point)
        exchange(:,1) = ice%srfoce_temp(1:nod2)
     case ('ocean_to_ice_bundle')
        exchange(:,1) = ice%srfoce_salt(1:nod2)  ! SSS [psu]
        exchange(:,2) = ice%srfoce_ssh(1:nod2)   ! SSH [m]
     case ('ocean_to_ice_uv')
        exchange(:,1) = ice%srfoce_u(1:nod2)     ! surface velocity [m/s]
        exchange(:,2) = ice%srfoce_v(1:nod2)
     case ('runoff_to_ice')
        ! river runoff [m/s]: with the built-in sea ice it enters the ocean
        ! through the ice thermodynamics only, so the sea ice adds it there
        exchange(:,1) = runoff(1:nod2)
     ! ICON: the atmospheric fluxes as received
     case ('taux_to_ice')
        exchange(:,1:2) = fwd_taux
     case ('tauy_to_ice')
        exchange(:,1:2) = fwd_tauy
     case ('surface_fresh_water_flux_to_ice')
        exchange(:,1:3) = fwd_fresh_water
     case ('total_heat_flux_to_ice')
        exchange(:,1:4) = fwd_heat_flux
     case ('atmosphere_sea_ice_bundle_to_ice')
        exchange(:,1:2) = fwd_atm_sea_ice_bundle
     ! IFS: the fluxes the interface deposited this step, in FESOM units
     case ('atm_ice_flux_to_ice')
        exchange(:,1)  = ice%stress_atmice_x(1:nod2) ! [Pa]
        exchange(:,2)  = ice%stress_atmice_y(1:nod2) ! [Pa]
        exchange(:,3)  = oce_heat_flux(1:nod2)       ! [W/m2]
        exchange(:,4)  = ice_heat_flux(1:nod2)       ! [W/m2]
        exchange(:,5)  = shortwave(1:nod2)           ! [W/m2]
        exchange(:,6)  = prec_rain(1:nod2)           ! [m/s]
        exchange(:,7)  = prec_snow(1:nod2)           ! [m/s]
        exchange(:,8)  = evap_no_ifrac(1:nod2)       ! [m/s]
        exchange(:,9)  = sublimation(1:nod2)         ! [m/s]
        exchange(:,10) = enthalpyoffuse(1:nod2)      ! [W/m2]
     ! no atmosphere: the raw state from the forcing files (FESIM converts)
     case ('atm_state_to_ice')
        exchange(:,1) = atmdata(i_xwind, 1:nod2)   ! u_wind    [m/s]
        exchange(:,2) = atmdata(i_ywind, 1:nod2)   ! v_wind    [m/s]
        exchange(:,3) = atmdata(i_tair , 1:nod2)   ! t_air     [K]
        exchange(:,4) = atmdata(i_humi , 1:nod2)   ! shum      [kg/kg]
        exchange(:,5) = atmdata(i_qsr  , 1:nod2)   ! shortwave [W/m2] down
        exchange(:,6) = atmdata(i_qlw  , 1:nod2)   ! longwave  [W/m2] down
        exchange(:,7) = atmdata(i_prec , 1:nod2)   ! prec_rain [kg/m2/s]
        exchange(:,8) = atmdata(i_snow , 1:nod2)   ! prec_snow [kg/m2/s]
        exchange(:,9) = atmdata(i_mslp , 1:nod2)   ! mslp      [Pa]
     case default
        write(*,*) 'exchange_oce_ice_yac: no handler for the sent field ', trim(cpl_send(i))
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
     end select
     call cpl_yac_send(i, exchange(:,1:cpl_send_collection_size(i)), action)
  enddo

  ! IFS: open-water flux split. With the concurrent (lag-1) exchange the sea
  ! ice returns only the ice-related part of its heat and freshwater flux
  ! (thermodynamics_coupled subtracts (1-A)*a2ohf and the open-water P-E plus
  ! runoff, A = its concentration before the thermodynamics = values_old
  ! received above). The ocean adds the open-water part from the atmospheric
  ! fluxes of THIS step, so ice-free water is forced without the exchange
  ! lag; oce_fluxes then applies heat_flux = -flx_h, water_flux = -flx_fw.
  if (is_coupled_to_ifs) then
     do n = 1, nod2+partit%eDim_nod2D
        if (mesh%ulevels_nod2D(n) > 1) then
           ice%flx_h(n)  = flx_h_ice(n)
           ice%flx_fw(n) = flx_fw_ice(n)
        else
           ice%flx_h(n)  = flx_h_ice(n)  + (1.0_WP-ice%data(1)%values_old(n))*(oce_heat_flux(n)+shortwave(n)+enthalpyoffuse(n))
           ice%flx_fw(n) = flx_fw_ice(n) + (evap_no_ifrac(n)+prec_rain(n)+(1.0_WP-ice%data(1)%values_old(n))*prec_snow(n)+runoff(n)) &
                                           *ice%thermo%rhofwt/ice%thermo%rhowat
           ! the terms oce_fluxes integrates for the global freshwater balancing, from
           ! the SAME step's atmospheric fluxes, as the built-in thermodynamics does
           evaporation(n)     = evap_no_ifrac(n) + sublimation(n)
           ice_sublimation(n) = sublimation(n)
        end if
     end do
  end if

end subroutine exchange_oce_ice_yac
#endif /* __cpl_yac */

!
!------------------------------------------------------------------------------------
!
#if defined (__cpl_coupler)
!
!=================================================================
!
! Modifies the flux 'field2d' in order to conserve 
! the net fluxes.
! We distinguish between NH and SH fluxes (integer 'h') 
! where possible.
!
!================================================================= 
! History :
!  07-11  (D.Sidorenko,	AWI Germany) first routine
!  10-12  (T.Rackow, 	AWI Germany) code reordering and cleanup  
!-----------------------------------------------------------------
!
SUBROUTINE force_flux_consv(field2d, mask, n, h, do_stats, partit, mesh)

  IMPLICIT NONE
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit  
  real(kind=WP),  INTENT (INOUT) 	:: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP),  INTENT (IN)	        :: mask(partit%myDim_nod2D+partit%eDim_nod2D)
  INTEGER,        INTENT (IN)	        :: n
  INTEGER,        INTENT (IN)           :: h !hemisphere: 0=GL, 1=NH, 2=SH
  logical,        INTENT (IN)	        :: do_stats
  real(kind=WP)			        :: rmask(partit%myDim_nod2D+partit%eDim_nod2D)   
  real(kind=WP)		                :: weight(partit%myDim_nod2D+partit%eDim_nod2D)    
  real(kind=WP)                         :: flux_global(2), flux_local(2)
  real(kind=WP)                         :: eff_vol(2)


  ! Only the ECHAM coupling redistributes the flux residual. OIFS relies on
  ! OASIS3-MCT conservative remapping instead, and the YAC path never calls
  ! this routine at all.
  if (.not. is_coupled_to_echam) return

  if (mstep==1) then
	if (partit%mype == 0) write(*,*) 'Do not apply a correction at first time step for ', trim(cpl_recv(n))
	return
  end if

  
  !just keep NH or SH part in mask
  rmask=mask
  SELECT CASE(h)
	   CASE( 1 ) ; where (mesh%geo_coord_nod2D(2, :)< 0) !just NH 
	                     rmask=0. 
		       end where
	   CASE( 2 ) ; where (mesh%geo_coord_nod2D(2, :)>=0) !just SH 
	                     rmask=0. 
		       end where
	   CASE DEFAULT 
	   !keep global mask
  END SELECT
 
  !residual (net) fluxes; computes also oce_net_fluxes_*
  call compute_residual(field2d, rmask, n, partit, mesh)
  
#ifdef VERBOSE
  if (partit%mype == 0) then
  !atm net fluxes, oce net fluxes before modification
  write(*,'(3A,3e15.7)') 'atm NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			    atm_net_fluxes_north(n), atm_net_fluxes_south(n), 	&
			    atm_net_fluxes_north(n)+ atm_net_fluxes_south(n)
			    
  write(*,'(3A,3e15.7)') 'oce NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			    oce_net_fluxes_north(n), oce_net_fluxes_south(n), 	&
			    oce_net_fluxes_north(n)+ oce_net_fluxes_south(n)
  end if			    
#endif
  if (do_stats)	then
     if (h==0 .or. h==1) then
        a2o_fcorr_stat(n, 1) =a2o_fcorr_stat(n, 1)+atm_net_fluxes_north(n)
        a2o_fcorr_stat(n, 4) =a2o_fcorr_stat(n, 4)+oce_net_fluxes_north(n)
     end if
     if (h==0 .or. h==2) then
        a2o_fcorr_stat(n, 2) =a2o_fcorr_stat(n, 2)+atm_net_fluxes_south(n)
        a2o_fcorr_stat(n, 5) =a2o_fcorr_stat(n, 5)+oce_net_fluxes_south(n)
     end if     
  end if

  !integrate (masked) abs(field2d) to get positive weights
  call integrate_2D(flux_global, flux_local, eff_vol, abs(field2d), rmask, partit, mesh)
  
  !get weight pattern with integral 1
  if (abs(sum(flux_global))>1.e-10) then
    weight=abs(field2d/sum(flux_global))
  else
    !should rarely happen
    weight=1.0_WP / sum(eff_vol)
    write(*,*) 'Warning: Constant redistribution for flux ', trim(cpl_recv(n))
  end if
  
  !weight is still global 2D field, just keep NH or SH part
  where (rmask<1.e-10)
    weight=0.
  end where
  
  !redistribute the residual according to the mask
  SELECT CASE(h)
	   CASE( 0 ) ; field2d=field2d+weight*flux_correction_total(n) ! GL	       
	   CASE( 1 ) ; field2d=field2d+weight*flux_correction_north(n) ! NH
	   CASE( 2 ) ; field2d=field2d+weight*flux_correction_south(n) ! SH
  END SELECT
  
  !check conservation
  call integrate_2D(flux_global, flux_local, eff_vol, field2d, rmask, partit, mesh)
#ifdef VERBOSE
  if (partit%mype == 0) then
  write(*,'(3A,3e15.7)') 'oce NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			   flux_global(1), flux_global(2), sum(flux_global)
  end if
#endif
  
  !last flux			   
  if (n==nrecv .AND. partit%mype==0) write(*,*) 'Fluxes have been modified.'  
END SUBROUTINE force_flux_consv

!
! Compute the difference between the net fluxes seen by the atmosphere
! and ocean component (residual flux) for flux n.
!
SUBROUTINE compute_residual(field2d, mask, n, partit, mesh)

 
  IMPLICIT NONE
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit  
  real(kind=WP), INTENT(IN)   :: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP), INTENT(IN)   :: mask(partit%myDim_nod2D+partit%eDim_nod2D)  
  INTEGER,       INTENT(IN)   :: n
  
  real(kind=WP)               :: flux_global(2), flux_local(2)
  real(kind=WP)               :: eff_vol(2)

  !compute net flux (for flux n) on ocean side
  call integrate_2D(flux_global, flux_local, eff_vol, field2d, mask, partit, mesh)
  oce_net_fluxes_north(n)=flux_global(1)
  oce_net_fluxes_south(n)=flux_global(2)
  
  !compute the residual fluxes for NH, SH and Global
  flux_correction_north(n)= atm_net_fluxes_north(n) - oce_net_fluxes_north(n)
  flux_correction_south(n)= atm_net_fluxes_south(n) - oce_net_fluxes_south(n)
  flux_correction_total(n)= flux_correction_north(n) + flux_correction_south(n)
END SUBROUTINE compute_residual

! -field_2d (input) is any (partitioned) 2D field
! -flux_local  (returned) is the net local flux (for current pc)
! -flux_global (returned) is the communicated and summarized flux_local  
!
SUBROUTINE integrate_2D(flux_global, flux_local, eff_vol, field2d, mask, partit, mesh)
  IMPLICIT NONE
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(in),    target :: partit
  real(kind=WP), INTENT(OUT)  :: flux_global(2), flux_local(2)
  real(kind=WP), INTENT(OUT)  :: eff_vol(2)
  real(kind=WP), INTENT(IN)   :: field2d(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP), INTENT(IN)   :: mask(partit%myDim_nod2D   +partit%eDim_nod2D) 
   
  real(kind=WP)               :: eff_vol_local(2)
  integer                     :: MPIerr


  flux_local(1)=sum(mesh%lump2d_north*field2d(1:partit%myDim_nod2D)*mask(1:partit%myDim_nod2D))
  flux_local(2)=sum(mesh%lump2d_south*field2d(1:partit%myDim_nod2D)*mask(1:partit%myDim_nod2D))
  call MPI_AllREDUCE(flux_local, flux_global, 2, &
  		     MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, MPIerr)
		     
		     
  eff_vol_local(1)=sum(mesh%lump2d_north*mask(1:partit%myDim_nod2D))
  eff_vol_local(2)=sum(mesh%lump2d_south*mask(1:partit%myDim_nod2D))
  call MPI_AllREDUCE(eff_vol_local, eff_vol,  2, & 
  		     MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, MPIerr)
		     
END SUBROUTINE integrate_2D
!
!---------------------------------------------------------------------------------------------------
!  Receive atmospheric net fluxes (atm_net_fluxes_north and atm_net_fluxes_south)
!
!SUBROUTINE net_rec_from_atm(action)
!  
!
!  use g_forcing_arrays,	only : 	atm_net_fluxes_north, atm_net_fluxes_south, 	&
!  				oce_net_fluxes_north, oce_net_fluxes_south, 	&
!				flux_correction_north, flux_correction_south,	&
!				flux_correction_total
!  use cpl_driver,	only : 	nrecv
!  
!  IMPLICIT NONE
!
!  LOGICAL, INTENT(IN)   :: action
!
!  if (action) then
!  atm_net_fluxes_north=0.
!  atm_net_fluxes_south=0.
!
!  atm_net_fluxes_north=(/ 1.,2.,3.,4. , 0.5e7, -0.1e10, -0.3e8, 0.5e-13, -0.7e17, -0.4e-6, 0.5e5, 0.0/)
!  atm_net_fluxes_south=(/ 1.,2.,3.,4. , 0.4e7, 0.3e06, -0.3e8, -0.3,-0.7e17, -0.4e-6, 0.2e6, 0.0/)
!  end if
!END SUBROUTINE net_rec_from_atm
!
!
!---------------------------------------------------------------------------------------------------
!  Receieve atmospheric net fluxes (atm_net_fluxes_north and atm_net_fluxes_south)
!
SUBROUTINE net_rec_from_atm(action, partit)
!


  IMPLICIT NONE

  LOGICAL,        INTENT (IN)   		  :: action
  type(t_partit), intent(inout), target           :: partit
  INTEGER                                         :: my_global_rank, ierror
#if defined(__recom) && defined(__usetp)
  INTEGER                                         :: my_global_rank_test
#endif
  INTEGER                                         :: n  
  INTEGER 					  :: status(MPI_STATUS_SIZE,partit%npes) 
  INTEGER                                         :: request(2)
  real(kind=WP)                 		  :: aux(nrecv)
  ! Only the ECHAM coupling exchanges net fluxes out of band. OIFS relies on
  ! OASIS3-MCT conservative remapping and receives no net fluxes here.
  if (.not. is_coupled_to_echam) return
  ! NOTE (single precision): the MPI_DOUBLE_PRECISION calls below are a RAW
  ! FESOM<->atmosphere root exchange over MPI_COMM_WORLD (source_root/target_root),
  ! NOT routed through OASIS, so they are only used by the ECHAM/AWICM flux-correction
  ! path -- other partners return above and never reach them. They are left
  ! hardcoded double on purpose: atm_net_fluxes_* are real(kind=WP), but the buffer
  ! kind must match the ATMOSPHERE partner (double), not the local WP, so this is the
  ! one coupling spot that MUST NOT be switched to MPI_WP. Making FESOM single
  ! precision coupled to ECHAM would additionally require this exchange (and the
  ! atmosphere side) to agree on a precision -- out of scope for the OIFS SP work.

  if (action) then
     CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_global_rank, ierror)
     atm_net_fluxes_north=0.
     atm_net_fluxes_south=0.
#if defined(__recom) && defined(__usetp)
     my_global_rank_test = my_global_rank - (partit%my_fesom_group * partit%npes)
#endif

#if defined(__recom) && defined(__usetp)
! check for is root in group
     if (my_global_rank_test==target_root) then
        if(partit%my_fesom_group == 0) then
#else
     if (my_global_rank==target_root) then
#endif
        CALL MPI_IRecv(atm_net_fluxes_north(1), nrecv, MPI_WP, source_root, 111, MPI_COMM_WORLD, request(1), partit%MPIerr)
        CALL MPI_IRecv(atm_net_fluxes_south(1), nrecv, MPI_WP, source_root, 112, MPI_COMM_WORLD, request(2), partit%MPIerr)
        CALL MPI_Waitall(2, request, status, partit%MPIerr)
     end if

#if defined(__recom) && defined(__usetp)
        if(num_fesom_groups > 1) then
           call MPI_Bcast(atm_net_fluxes_north(1), nrecv, MPI_WP, 0, partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, partit%MPIerr)
           call MPI_Bcast(atm_net_fluxes_south(1), nrecv, MPI_WP, 0, partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, partit%MPIerr)
        end if
     end if ! (my_global_rank_test==target_root) then
#endif
  call MPI_AllREDUCE(atm_net_fluxes_north(1), aux, nrecv, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  atm_net_fluxes_north=aux
  call MPI_AllREDUCE(atm_net_fluxes_south(1), aux, nrecv, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  atm_net_fluxes_south=aux
  end if
END SUBROUTINE net_rec_from_atm
#endif

end module gen_forcing_couple_module

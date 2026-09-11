! Ice->ocean coupling, FESIM (sea-ice component) side.
!
! Only oce_fluxes_mom lives here: it forms the ice-ocean momentum drag into
! ice%stress_iceoce_x/y, which the ocean receives over YAC as
! `ice_to_ocean_stress`. It reads the ocean surface velocity from
! ice%srfoce_u/v -- i.e. the values FESIM received over YAC -- and never
! touches the ocean's own 3-D state.
!
! Removed 2026-09-11 (ocean-leftover audit): `ocean2ice` and `oce_fluxes`, both
! dead here. Their call sites in fesim_module were `!SL`-commented when the ice
! was split out, because in the decoupled design the ocean surface state arrives
! over YAC (ocean_coupling_interface -> update_atm_forcing_yac) and the ice->ocean
! heat/freshwater flux is sent as `ice_to_ocean_flux` from gen_forcing_couple,
! not applied to a local ocean here. Between them they carried ~690 lines of
! ocean-only code -- 3-D temp/salt/UV pointers, water-isotope and age-tracer
! handling, virtual salt flux, cavity and SPP paths -- none of it reachable.
! They remain in the ocean tree as `ice_oce_coupling.F90`, which is where they
! belong. Module `ocean2ice_interface` went with them.

module oce_fluxes_interface
    interface
        
        subroutine oce_fluxes_mom(ice, dynamics, partit, mesh)
        USE MOD_ICE
        USE MOD_DYN
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_dyn)   , intent(in)   , target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine oce_fluxes_mom
    end interface
end module oce_fluxes_interface

!
!
!_______________________________________________________________________________
! transmits the relevant fields from the ice to the ocean model
subroutine oce_fluxes_mom(ice, dynamics, partit, mesh)
    USE MOD_ICE
    USE MOD_DYN
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use o_ARRAYS
    USE g_CONFIG
    use g_comm_auto
    use cavity_interfaces    
#if defined (__icepack)
    use icedrv_main,   only: icepack_to_fesom
#endif
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_dyn)   , intent(in)   , target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: n, elem, elnodes(3),n1
    real(kind=WP)            :: aux
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice, a_ice, u_w, v_w
    real(kind=WP), dimension(:), pointer  :: stress_iceoce_x, stress_iceoce_y  
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice           => ice%uice(:)
    v_ice           => ice%vice(:)
    a_ice           => ice%data(1)%values(:)
    u_w             => ice%srfoce_u(:)
    v_w             => ice%srfoce_v(:)
    stress_iceoce_x => ice%stress_iceoce_x(:)
    stress_iceoce_y => ice%stress_iceoce_y(:)
  
    ! ==================
    ! momentum flux:
    ! ==================
    !___________________________________________________________________________

#if defined (__icepack)
     call icepack_to_fesom(nx_in=(myDim_nod2D+eDim_nod2D), &
                           aice_out=a_ice)
#endif
    !___________________________________________________________________________
    ! compute total surface stress (iceoce+atmoce) on nodes 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, elem, elnodes, n1, aux)
!$OMP DO
    do n=1,myDim_nod2D+eDim_nod2D   
        !_______________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(n)>1) cycle
        
        !_______________________________________________________________________
        if(a_ice(n)>0.001_WP) then
            aux=sqrt((u_ice(n)-u_w(n))**2+(v_ice(n)-v_w(n))**2)*density_0*ice%cd_oce_ice
            stress_iceoce_x(n) = aux * (u_ice(n)-u_w(n))
            stress_iceoce_y(n) = aux * (v_ice(n)-v_w(n))
        else
            stress_iceoce_x(n)=0.0_WP
            stress_iceoce_y(n)=0.0_WP
        end if
        
        stress_node_surf(1,n) = stress_iceoce_x(n)*a_ice(n) + stress_atmoce_x(n)*(1.0_WP-a_ice(n))
        stress_node_surf(2,n) = stress_iceoce_y(n)*a_ice(n) + stress_atmoce_y(n)*(1.0_WP-a_ice(n))
    end do
!$OMP END DO
    !___________________________________________________________________________
    ! compute total surface stress (iceoce+atmoce) on elements
!$OMP DO
    DO elem=1,myDim_elem2D
        !_______________________________________________________________________
        ! if cavity element skip it 
        if (ulevels(elem)>1) cycle
        
        !_______________________________________________________________________
        ! total surface stress (iceoce+atmoce) on elements 
        elnodes=elem2D_nodes(:,elem)
        
        !!PS stress_surf(1,elem)=sum(stress_iceoce_x(elnodes)*a_ice(elnodes) + &
        !!PS                         stress_atmoce_x(elnodes)*(1.0_WP-a_ice(elnodes)))/3.0_WP
        !!PS stress_surf(2,elem)=sum(stress_iceoce_y(elnodes)*a_ice(elnodes) + &
        !!PS                         stress_atmoce_y(elnodes)*(1.0_WP-a_ice(elnodes)))/3.0_WP
        stress_surf(1,elem)=sum(stress_node_surf(1,elnodes))/3.0_WP
        stress_surf(2,elem)=sum(stress_node_surf(2,elnodes))/3.0_WP

    END DO
!$OMP END DO
!$OMP END PARALLEL
    !___________________________________________________________________________
    if (use_cavity) call cavity_momentum_fluxes(dynamics, partit, mesh)
  
end subroutine oce_fluxes_mom

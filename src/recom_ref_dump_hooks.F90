!==============================================================================
! FESOM-side hooks of the recom_ref_dump reference dumps (see
! src/recom/recom_ref_dump.F90). Every routine returns immediately unless the
! current step is being dumped (FESOM_RECOM_DUMP_DIR / FESOM_RECOM_DUMP_STEPS),
! and only reads model arrays.
! Names: <stage>.<field>; per tracer <stage>.tr<ID>[.<field>]. Node fields are
! written over the owned nodes 1:myDim_nod2D (= all nodes on one rank),
! element fields over 1:myDim_elem2D.
!==============================================================================
module recom_ref_dump_hooks
#if defined(__recom)
    use recom_ref_dump
    use MOD_MESH
    use MOD_PARTIT
    use MOD_TRACER
    use MOD_DYN
    use MOD_ICE
    use o_ARRAYS, only: Kv, Ki, slope_tapered
    use g_forcing_arrays, only: press_air, u_wind, v_wind, shortwave
    use recom_glovar
    use recom_config, only: benthos_num
#if defined(__RECOM_WAVEBANDS)
    use REcoM_spectral
#endif
    implicit none
    private
    public :: rdh_step_begin, rdh_recom_in, rdh_recom_out, rdh_step_end, rdh_transport, &
              rdh_tracer, rdh_tracer_field, rdh_tracer_benthos, rdh_slp_force

    character(len=*), parameter :: fmt_id = '(a,a,i0)'

contains

    function trname(stage, id) result(s)
        character(len=*), intent(in) :: stage
        integer, intent(in) :: id
        character(len=64) :: s
        write(s, fmt_id) trim(stage), '.tr', id
    end function trname

    !--------------------------------------------------------------------------
    subroutine rdh_step_begin(mstep, yearnew, daynew, timenew, month, ndpyr, dt, partit, mesh)
        integer, intent(in) :: mstep, yearnew, daynew, month, ndpyr
        real(kind=WP), intent(in) :: timenew, dt
        type(t_partit), intent(in), target :: partit
        type(t_mesh), intent(in), target :: mesh
        integer :: nn, ne
        call rdump_init(partit%npes)
        if (.not. rdump_enabled) return
        call rdump_step(mstep, yearnew, daynew, real(timenew, 8), month, ndpyr, real(dt, 8), &
                        mesh%nod2D, mesh%elem2D, mesh%nl)
        if (.not. rdump_active()) return
        nn = partit%myDim_nod2D
        ne = partit%myDim_elem2D
        call rdump('mesh.gid_nod', partit%myList_nod2D(1:nn))
        call rdump('mesh.gid_elem', partit%myList_elem2D(1:ne))
        call rdump('mesh.ulevels_nod2D', mesh%ulevels_nod2D(1:nn))
        call rdump('mesh.nlevels_nod2D', mesh%nlevels_nod2D(1:nn))
        call rdump('mesh.ulevels', mesh%ulevels(1:ne))
        call rdump('mesh.nlevels', mesh%nlevels(1:ne))
        call rdump('mesh.areasvol', mesh%areasvol(:, 1:nn))
        call rdump('mesh.area', mesh%area(:, 1:nn))
        call rdump('mesh.elem_area', mesh%elem_area(1:ne))
        call rdump('mesh.geo_coord_nod2D', mesh%geo_coord_nod2D(:, 1:nn))
        call rdump('mesh.ocean_area', mesh%ocean_area)
    end subroutine rdh_step_begin

    !--------------------------------------------------------------------------
    subroutine dump_all_tracers(stage, tracers, partit)
        character(len=*), intent(in) :: stage
        type(t_tracer), intent(in), target :: tracers
        type(t_partit), intent(in), target :: partit
        integer :: tr
        do tr = 1, tracers%num_tracers
            call rdump(trname(stage, tracers%data(tr)%ID), tracers%data(tr)%values(:, 1:partit%myDim_nod2D))
        end do
    end subroutine dump_all_tracers

    !--------------------------------------------------------------------------
    ! Everything recom() reads, just before the call (fesom_module.F90).
    subroutine rdh_recom_in(ice, tracers, partit, mesh, dt, kappa, rad, daynew, month, mstep, &
                            ndpyr, yearold, timenew)
        type(t_ice), intent(in), target :: ice
        type(t_tracer), intent(in), target :: tracers
        type(t_partit), intent(in), target :: partit
        type(t_mesh), intent(in), target :: mesh
        real(kind=WP), intent(in) :: dt, kappa, rad, timenew
        integer, intent(in) :: daynew, month, mstep, ndpyr, yearold
        integer :: nn
        if (.not. rdump_active()) return
        nn = partit%myDim_nod2D
        call rdump('recom_in.dt', dt)
        call rdump('recom_in.kappa', kappa)
        call rdump('recom_in.rad', rad)
        call rdump('recom_in.timenew', timenew)
        call rdump('recom_in.calendar', (/daynew, month, mstep, ndpyr, yearold/))
        call rdump('recom_in.a_ice', ice%data(1)%values(1:nn))
        call rdump('recom_in.m_ice', ice%data(2)%values(1:nn))
        call rdump('recom_in.m_snow', ice%data(3)%values(1:nn))
        call rdump('recom_in.press_air', press_air(1:nn))
        call rdump('recom_in.u_wind', u_wind(1:nn))
        call rdump('recom_in.v_wind', v_wind(1:nn))
        call rdump('recom_in.shortwave', shortwave(1:nn))
        call rdump('recom_in.hnode', mesh%hnode(:, 1:nn))
        call rdump('recom_in.z_3d_n', mesh%z_3d_n(:, 1:nn))
        call rdump('recom_in.zbar_3d_n', mesh%zbar_3d_n(:, 1:nn))
        call dump_all_tracers('recom_in', tracers, partit)
        call rdump('recom_in.Benthos', Benthos(1:nn, :))
        call rdump('recom_in.GloHplus', GloHplus(1:nn))
        call rdump('recom_in.Alk_surf', Alk_surf(1:nn))
        call rdump('recom_in.GloFeDust', GloFeDust(1:nn))
        call rdump('recom_in.GloNDust', GloNDust(1:nn))
        call rdump('recom_in.AtmCO2', AtmCO2(:))
        call dump_carbonate('recom_in', nn)
    end subroutine rdh_recom_in

    subroutine dump_carbonate(stage, nn)
        character(len=*), intent(in) :: stage
        integer, intent(in) :: nn
        call rdump(stage//'.CO23D', CO23D(:, 1:nn))
        call rdump(stage//'.pH3D', pH3D(:, 1:nn))
        call rdump(stage//'.pCO23D', pCO23D(:, 1:nn))
        call rdump(stage//'.HCO33D', HCO33D(:, 1:nn))
        call rdump(stage//'.CO33D', CO33D(:, 1:nn))
        call rdump(stage//'.OmegaC3D', OmegaC3D(:, 1:nn))
        call rdump(stage//'.kspc3D', kspc3D(:, 1:nn))
        call rdump(stage//'.rhoSW3D', rhoSW3D(:, 1:nn))
    end subroutine dump_carbonate

    !--------------------------------------------------------------------------
    ! Everything recom() produces, just after the call.
    subroutine rdh_recom_out(tracers, partit)
        type(t_tracer), intent(in), target :: tracers
        type(t_partit), intent(in), target :: partit
        integer :: nn
        if (.not. rdump_active()) return
        nn = partit%myDim_nod2D
        call dump_all_tracers('recom_out', tracers, partit)
        call rdump('recom_out.Benthos', Benthos(1:nn, :))
        call rdump('recom_out.GlodecayBenthos', GlodecayBenthos(1:nn, :))
        call rdump('recom_out.GloHplus', GloHplus(1:nn))
        call dump_carbonate('recom_out', nn)
        call rdump('recom_out.relax_alk', relax_alk(1:nn))
        call rdump('recom_out.virtual_alk', virtual_alk(1:nn))
        call rdump('recom_out.AtmFeInput', AtmFeInput(1:nn))
        call rdump('recom_out.AtmNInput', AtmNInput(1:nn))
        call rdump('recom_out.GloPCO2surf', GloPCO2surf(1:nn))
        call rdump('recom_out.GlodPCO2surf', GlodPCO2surf(1:nn))
        call rdump('recom_out.GloCO2flux', GloCO2flux(1:nn))
        call rdump('recom_out.GloCO2flux_seaicemask', GloCO2flux_seaicemask(1:nn))
        call rdump('recom_out.GloO2flux', GloO2flux(1:nn))
        call rdump('recom_out.GloO2flux_seaicemask', GloO2flux_seaicemask(1:nn))
        call rdump('recom_out.x_co2atm', x_co2atm(1:nn))
        call rdump('recom_out.PAR3D', PAR3D(:, 1:nn))
        ! 2-D diagnostics (Diags = .true.); the fork (fesom2.6_recom_tra_diags_spectral) has per-PFT
        ! NPP/GPP/NNA/Chldeg for n, d, c, p and 3-D grazing by mesozooplankton per prey, no PistonVelocity,
        ! alphaCO2 or grazmacro/grazmicro arrays (REcoM 692f804 names kept where they exist)
        call rdump('recom_out.NPPn', NPPn(1:nn))
        call rdump('recom_out.NPPd', NPPd(1:nn))
        call rdump('recom_out.GPPn', GPPn(1:nn))
        call rdump('recom_out.GPPd', GPPd(1:nn))
        call rdump('recom_out.NNAn', NNAn(1:nn))
        call rdump('recom_out.NNAd', NNAd(1:nn))
        call rdump('recom_out.Chldegn', Chldegn(1:nn))
        call rdump('recom_out.Chldegd', Chldegd(1:nn))
        if (allocated(NPPc)) then
            call rdump('recom_out.NPPc', NPPc(1:nn))
            call rdump('recom_out.NPPp', NPPp(1:nn))
            call rdump('recom_out.GPPc', GPPc(1:nn))
            call rdump('recom_out.GPPp', GPPp(1:nn))
            call rdump('recom_out.NNAc', NNAc(1:nn))
            call rdump('recom_out.NNAp', NNAp(1:nn))
            call rdump('recom_out.Chldegc', Chldegc(1:nn))
            call rdump('recom_out.Chldegp', Chldegp(1:nn))
        end if
        ! 3-D diagnostics
        call rdump('recom_out.grazmeso_tot', grazmeso_tot(:, 1:nn))
        call rdump('recom_out.grazmeso_n', grazmeso_n(:, 1:nn))
        call rdump('recom_out.grazmeso_d', grazmeso_d(:, 1:nn))
        if (allocated(grazmeso_c)) call rdump('recom_out.grazmeso_c', grazmeso_c(:, 1:nn))
        call rdump('recom_out.aggn', aggn(:, 1:nn))
        call rdump('recom_out.aggd', aggd(:, 1:nn))
        call rdump('recom_out.docexn', docexn(:, 1:nn))
        call rdump('recom_out.docexd', docexd(:, 1:nn))
        call rdump('recom_out.respn', respn(:, 1:nn))
        call rdump('recom_out.respd', respd(:, 1:nn))
        call rdump('recom_out.NPPn3D', NPPn3D(:, 1:nn))
        call rdump('recom_out.NPPd3D', NPPd3D(:, 1:nn))
        call rdump('recom_out.calcdiss', calcdiss(:, 1:nn))
        call rdump('recom_out.calcif', calcif(:, 1:nn))
        call rdump('recom_out.respmeso', respmeso(:, 1:nn))
        call rdump('recom_out.respmacro', respmacro(:, 1:nn))
        call rdump('recom_out.respmicro', respmicro(:, 1:nn))
        call rdump('recom_out.TPhyCO2', TPhyCO2(:, 1:nn))
        call rdump('recom_out.TDiaCO2', TDiaCO2(:, 1:nn))
        call rdump('recom_out.TCphotLigLim_phyto', TCphotLigLim_phyto(:, 1:nn))
        call rdump('recom_out.TCphot_phyto', TCphot_phyto(:, 1:nn))
        call rdump('recom_out.TCphotLigLim_diatoms', TCphotLigLim_diatoms(:, 1:nn))
        call rdump('recom_out.TCphot_diatoms', TCphot_diatoms(:, 1:nn))
        call rdump('recom_out.TTemp_phyto', TTemp_phyto(:, 1:nn))
        call rdump('recom_out.TSi_assimDia', TSi_assimDia(:, 1:nn))
        if (allocated(NPPc3D)) then
            call rdump('recom_out.NPPc3D', NPPc3D(:, 1:nn))
            call rdump('recom_out.NPPp3D', NPPp3D(:, 1:nn))
            call rdump('recom_out.aggc', aggc(:, 1:nn))
            call rdump('recom_out.aggp', aggp(:, 1:nn))
            call rdump('recom_out.docexc', docexc(:, 1:nn))
            call rdump('recom_out.docexp', docexp(:, 1:nn))
            call rdump('recom_out.respc', respc(:, 1:nn))
            call rdump('recom_out.respp', respp(:, 1:nn))
            call rdump('recom_out.TCoccoCO2', TCoccoCO2(:, 1:nn))
            call rdump('recom_out.TPhaeoCO2', TPhaeoCO2(:, 1:nn))
            call rdump('recom_out.TCphotLigLim_cocco', TCphotLigLim_cocco(:, 1:nn))
            call rdump('recom_out.TCphot_cocco', TCphot_cocco(:, 1:nn))
            call rdump('recom_out.TCphotLigLim_phaeo', TCphotLigLim_phaeo(:, 1:nn))
            call rdump('recom_out.TCphot_phaeo', TCphot_phaeo(:, 1:nn))
        end if
#if defined(__RECOM_WAVEBANDS)
        ! the spectral light of this step (REcoM_spectral): the five planes per band, the surface fields
        if (allocated(Edz3D)) then
            call rdump('recom_out.Edz3D', Edz3D(:, 1:nn, :))
            call rdump('recom_out.Esz3D', Esz3D(:, 1:nn, :))
            call rdump('recom_out.Euz3D', Euz3D(:, 1:nn, :))
            call rdump('recom_out.Eutop3D', Eutop3D(:, 1:nn, :))
            call rdump('recom_out.Estop3D', Estop3D(:, 1:nn, :))
        end if
        if (allocated(Ed0_2D)) then
            call rdump('recom_out.Ed0_2D', Ed0_2D(1:nn, :))
            call rdump('recom_out.Es0_2D', Es0_2D(1:nn, :))
            call rdump('recom_out.Eu0_2D', Eu0_2D(1:nn, :))
        end if
        if (allocated(oasim_ed2D)) then
            call rdump('recom_out.oasim_ed2D', oasim_ed2D(1:nn, :))
            call rdump('recom_out.oasim_es2D', oasim_es2D(1:nn, :))
        end if
#endif
    end subroutine rdh_recom_out

    !--------------------------------------------------------------------------
    ! Transport fields of the step, at the start of solve_tracers_ale (before the
    ! bolus velocity is added to uv / w / w_e).
    subroutine rdh_transport(dynamics, partit, mesh)
        type(t_dyn), intent(in), target :: dynamics
        type(t_partit), intent(in), target :: partit
        type(t_mesh), intent(in), target :: mesh
        integer :: nn, ne
        if (.not. rdump_active()) return
        nn = partit%myDim_nod2D
        ne = partit%myDim_elem2D
        call rdump('transport.uv', dynamics%uv(:, :, 1:ne))
        call rdump('transport.w', dynamics%w(:, 1:nn))
        call rdump('transport.w_e', dynamics%w_e(:, 1:nn))
        call rdump('transport.w_i', dynamics%w_i(:, 1:nn))
        if (allocated(dynamics%fer_uv)) call rdump('transport.fer_uv', dynamics%fer_uv(:, :, 1:ne))
        if (allocated(dynamics%fer_w)) call rdump('transport.fer_w', dynamics%fer_w(:, 1:nn))
        call rdump('transport.hnode', mesh%hnode(:, 1:nn))
        call rdump('transport.hnode_new', mesh%hnode_new(:, 1:nn))
        call rdump('transport.helem', mesh%helem(:, 1:ne))
        call rdump('transport.zbar_3d_n', mesh%zbar_3d_n(:, 1:nn))
        call rdump('transport.z_3d_n', mesh%z_3d_n(:, 1:nn))
        call rdump('transport.Kv', Kv(:, 1:nn))
        if (allocated(Ki)) call rdump('transport.Ki', Ki(:, 1:nn))
        if (allocated(slope_tapered)) call rdump('transport.slope_tapered', slope_tapered(:, :, 1:nn))
    end subroutine rdh_transport

    !--------------------------------------------------------------------------
    ! One tracer's own fields inside the tracer loop: values (and optionally
    ! valuesold, valuesAB, del_ttf) at a named stage.
    subroutine rdh_tracer(stage, tr_num, tracers, partit, with_old, with_ab, with_del)
        character(len=*), intent(in) :: stage
        integer, intent(in) :: tr_num
        type(t_tracer), intent(in), target :: tracers
        type(t_partit), intent(in), target :: partit
        logical, intent(in) :: with_old, with_ab, with_del
        integer :: nn
        character(len=64) :: base
        if (.not. rdump_active()) return
        nn = partit%myDim_nod2D
        base = trname(stage, tracers%data(tr_num)%ID)
        call rdump(trim(base)//'.values', tracers%data(tr_num)%values(:, 1:nn))
        if (with_old) call rdump(trim(base)//'.valuesold', tracers%data(tr_num)%valuesold(1, :, 1:nn))
        if (with_ab) call rdump(trim(base)//'.valuesAB', tracers%data(tr_num)%valuesAB(:, 1:nn))
        if (with_del) call rdump(trim(base)//'.del_ttf', tracers%work%del_ttf(:, 1:nn))
    end subroutine rdh_tracer

    ! A REcoM work array of the tracer loop (dtr_bf, vert_sink, str_bf, ballasting).
    subroutine rdh_tracer_field(stage, field, id, a, nn)
        character(len=*), intent(in) :: stage, field
        integer, intent(in) :: id, nn
        real(kind=WP), intent(in) :: a(:, :)
        if (.not. rdump_active()) return
        call rdump(trim(trname(stage, id))//'.'//field, a(:, 1:nn))
    end subroutine rdh_tracer_field

    ! The benthos of the sinking tracer: Benthos_tr(:, :, tr_num), the per-tracer
    ! buffer summed into Benthos after the whole tracer loop (__usetp), and
    ! Benthos itself (unchanged inside the loop under __usetp).
    subroutine rdh_tracer_benthos(id, tr_num, nn)
        integer, intent(in) :: id, tr_num, nn
        if (.not. rdump_active()) return
#if defined(__usetp)
        call rdump(trim(trname('tra_sink', id))//'.Benthos_tr', Benthos_tr(1:nn, :, tr_num))
#endif
        call rdump(trim(trname('tra_sink', id))//'.Benthos', Benthos(1:nn, :))
    end subroutine rdh_tracer_benthos

    !--------------------------------------------------------------------------
    ! The SLP contribution of compute_vel_rhs (oce_ale_vel_rhs.F90, l_mslp):
    ! -grad(press_air/1000) per element [m/s2], the same arithmetic as the
    ! momentum code applies to p_air alone; the rhs adds it times elem_area on
    ! every layer.
    subroutine rdh_slp_force(l_mslp, partit, mesh)
        logical, intent(in) :: l_mslp
        type(t_partit), intent(in), target :: partit
        type(t_mesh), intent(in), target :: mesh
        integer :: elem, ne, elnodes(3)
        real(kind=WP) :: p_air(3), pre(3)
        real(kind=WP), allocatable :: force(:, :)
        if (.not. rdump_active()) return
        ne = partit%myDim_elem2D
        allocate(force(2, ne))
        do elem = 1, ne
            elnodes = mesh%elem2D_nodes(:, elem)
            if (l_mslp) then
                p_air = press_air(elnodes) / 1000
            else
                p_air = 0.0_WP
            end if
            pre = -p_air
            force(1, elem) = sum(mesh%gradient_sca(1:3, elem) * pre)
            force(2, elem) = sum(mesh%gradient_sca(4:6, elem) * pre)
        end do
        call rdump('vel_rhs.slp_force', force)
        deallocate(force)
    end subroutine rdh_slp_force

    !--------------------------------------------------------------------------
    ! State at the end of the step (after oce_timestep_ale).
    subroutine rdh_step_end(tracers, partit)
        type(t_tracer), intent(in), target :: tracers
        type(t_partit), intent(in), target :: partit
        integer :: nn
        if (.not. rdump_active()) return
        nn = partit%myDim_nod2D
        call dump_all_tracers('step_end', tracers, partit)
        call rdump('step_end.Benthos', Benthos(1:nn, :))
    end subroutine rdh_step_end
#endif
end module recom_ref_dump_hooks

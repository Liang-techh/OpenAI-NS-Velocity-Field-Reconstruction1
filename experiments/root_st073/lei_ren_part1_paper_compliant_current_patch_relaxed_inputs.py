"""Whole current implicit patch relaxed inputs, including partial-history gaps.

The fixed bump estimates are universal hypotheses, not admission of a legacy
field. Here they consume the checked current correlated source, implicit
coefficient family, inherited five primitives and unchanged analytic P0.
"""
import json
import math
from pathlib import Path
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_restore_relaxed_inputs as restore

packets=restore.packets;HERE,PREFIX,sha=restore.HERE,restore.PREFIX,restore.sha
NAME=PREFIX+'current_patch_relaxed_inputs.json'
RECEIPT=PREFIX+'current_patch_relaxed_inputs_check.json'
GATE='current_original_Rm_Rh_relaxed_generic_input_certified'
OPEN=restore.OPEN
EDGES=(49,51,59,61,69,71)


def patch_source_theorem():
    asts=packets.recovery.numeric.transport.SourceAST()
    wanted={
        'H':"f+x**c.mpf('.1')",'u':'am*H',
        'uy':"am*(fx*x+x**c.mpf('.1')/10)",
        'V':'z*4+g','Vy':'gx*x',
        'defects':'[a+b for a,b in zip(d,change)]',
        'mass':'z*4+defects[0]/x',
        'theta':"defects[2]+x**c.mpf('1.6')*c.mpf('.625')",
        'mixed':'(z*theta)*4+defects[1]',
        'energy':"defects[3]-x**c.mpf('1.2')*c.mpf(5)/12",
        'pressure_moment':"defects[4]+x**c.mpf('.2')*c.mpf('2.5')",
        'P':'p0+square(am)*pressure_moment',
        'terminal':'endpoints(x)[0]>=mp.mpf(71)/40'}
    nodes={name:asts.expression('actual_moment_patch','evaluate',name,wanted=value)
        for name,value in wanted.items()}
    x,Z,S=s.symbols('x Z Pstar',positive=True)
    f=s.Function('same_three_bump_f')(x,Z);g=s.Function('same_two_bump_g')(x,Z)
    am=S*s.exp(-s.Rational(3,5))/(1+Z*Z)
    env=dict(x=x,z=Z,am=am,f=f,fx=s.diff(f,x),g=g,gx=s.diff(g,x),
        c=type('ExactConstants',(),{'mpf':staticmethod(lambda v:s.Rational(str(v)))})())
    H=asts.evaluate(nodes['H'],env);env['H']=H
    u=asts.evaluate(nodes['u'],env);uy=asts.evaluate(nodes['uy'],env)
    V=asts.evaluate(nodes['V'],env);Vy=asts.evaluate(nodes['Vy'],env)
    a=1-2*uy/u;b=2*Vy/u
    checks={}
    def zero(name,value):
        if s.simplify(value)!=0:raise ArithmeticError('Original actual patch identity: '+name)
        checks[name]=True
    zero('original_same_full_swirl_and_radial_derivative',uy-x*s.diff(u,x))
    zero('original_same_full_axial_radial_derivative',Vy-x*s.diff(V,x))
    zero('original_full_signed_a',a-(s.Rational(4,5)+(f/5-2*x*s.diff(f,x))/H))
    zero('original_full_signed_b',b-2*x*s.diff(g,x)/u)
    zero('original_source_axial_zeta',s.diff(s.log(u),Z)+2*Z/(1+Z*Z)-s.diff(f,Z)/H)
    p1,p2,Q,N,R,L,J=s.symbols('p1 p2 Qangular N R L J',nonzero=True)
    bw=2*J*Vy/Q
    rules={p1:R*Q/L,p2:R*N/(L*u),J:N/u**2}
    zero('original_full_signed_bw_uses_axial_pressure_stress',(bw-b*p2/p1).subs(rules))
    zero('original_full_H0_uses_signed_p2',(p1-p2*b/a-p1*(a-bw)/a).subs(rules))
    # No gap can reset accumulated primitives: only local f,g vanish there.
    ds=s.symbols('d0:5');cs=s.symbols('partial_change0:5')
    defects=asts.evaluate(nodes['defects'],dict(d=ds,change=cs,zip=zip))
    zero('original_gap_preserves_nonzero_partial_defects',sum(defects)-sum(ds)-sum(cs))
    return dict(passed=True,original_patch_AST_bindings=asts.bindings,exact_identities=checks,
        full_signed_N_and_Q_theorem_inherited_from_same_five_ODEs=True,
        local_zero_bumps_do_not_zero_partial_history_or_pressure=True,
        terminal_closure_only_after_exact_last_support_edge=True,input_hashes=asts.hashes)


class CurrentPatchRelaxedInputs:
    def __init__(self,service=None):
        self.restore=restore.CurrentRestoreRelaxedInputs(service)
        self.service=self.restore.service;self.ctx=self.restore.ctx;self.family=self.restore.family
        self.conditions=[];self.rows={}
        record=json.loads((HERE/restore.RECEIPT).read_bytes())
        if not record['all_passed'] or not record[restore.GATE] or record['source_family']!=self.family:
            raise ValueError('Same checked current original Rz..Rm input required')
        self.service.bind_hashes(record['input_hashes']);self.service.bind_hashes({restore.RECEIPT:sha(restore.RECEIPT)})
        def load(stem,compliant=True):
            name=(PREFIX if compliant else 'lei_ren_part1_paper_')+stem+'.json'
            row=json.loads((HERE/name).read_bytes())
            self.service.bind_hashes(row['input_hashes']);self.service.bind_hashes({name:sha(name)})
            self.rows[stem]=row
            return row
        def yes(row,key,label):
            if row[key] is not True:raise ValueError('Actual patch attachment missing: '+label)
            self.conditions.append(label)
        def yesmap(row,key,label):
            if not row[key] or not all(row[key].values()):raise ValueError('Current original patch source differs: '+label)
            self.conditions.append(label)
        patch=load('current_actual_patch_background_tensor_check')
        if {k:patch[k] for k in self.family}!=self.family:raise ValueError('Current full patch source/pressure family differs')
        yes(patch,'all_passed','checked actual current full patch source')
        for key in ('current_actual_five_moment_patch_full_tensor_available',
            'current_actual_five_moment_patch_full_meridional_decomposition_available',
            'current_actual_patch_Rh_completed_tensor_join_certified',
            'current_actual_patch_six_internal_support_tensor_traces_certified'):
            yes(patch,key,'same original whole patch and source traces/'+key)
        theorem=patch['current_actual_patch_source_units_and_tensor_theorem']
        normalized=theorem['original_actual_patch_current_radius_normalization_theorem']
        yes(normalized,'passed','same current original full patch units and five primitive ODEs')
        for key in ('raw_five_history_ODE_m','raw_five_history_ODE_h','raw_five_history_ODE_k',
            'raw_five_history_ODE_e','raw_five_history_ODE_p','absolute_pressure_y_ODE'):
            yes(normalized['identities'],key,'original exact patch history/'+key)
        source=theorem['current_actual_patch_source_pressure_support_and_Rh_theorem']
        yes(source,'passed','same whole original current patch/Rh defining functions')
        yesmap(source,'live_original_callable_bindings','same actual current mixed/moment owners and unchanged algorithms')
        for key in ('current_terminal_closure_receipt_consumed',
            'current_Rh_open_neighborhood_uses_same_unique_implicit_full_weight_closure',
            'original_beta_flat_derivatives_0_through_4_and_continuous_partial_integrals_consumed',
            'same_analytic_P0_function_and_first_six_projection_not_hash_only',
            'interval_overlap_not_used_as_function_identity'):
            yes(source,key,'same current compact-support/P0/implicit function/'+key)
        if source['six_exact_rational_support_edges']!=[str(s.Rational(i,40)) for i in EDGES]:
            raise ValueError('Actual original three disjoint support groups changed')
        current=load('actual_feedback_moment_patch_check');report=load('actual_feedback_moment_patch')
        for row in (current,report):
            if (row['implicit_source_sha256']!=self.family['implicit_source_sha256']
                or row['datum_enclosure_sha256']!=self.family['datum_enclosure_sha256']):
                raise ValueError('Actual patch current implicit source and analytic P0 differ')
            for key in ('current_actual_moment_patch_installed','current_actual_patch_implicit_axial5_recomputed',
                'current_five_functional_terminal_identities_connected'):
                yes(row,key,'same unique current implicit coefficient source/'+key)
            yes(row['original_correlated_source_class_admission'],'current_C2_bound_controls_original_summed_C1_norm',
                'current true source C2 admits original summed C1 class; independent boxes are not norms')
        for key in ('current_Rm_source_parent_and_P0_retained','residual_zero_containment_not_used_as_closure_proof'):
            yes(current,key,'current defect/partial-history function/'+key)
        if current['current_source_tail_C1_caps_checked']!=5 or current['current_factored_source_decay_caps_checked']!=36:
            raise ValueError('All current exact tail and decay bounds required')
        self.conditions.append('all five current C1 tails and 36 factored decay bounds')
        if current['current_actual_defect_axial5_coefficients_checked']!=30 or current['current_implicit_control_axial5_coefficients_checked']!=30:
            raise ValueError('All original actual defects and unique implicit controls required')
        fixed=load('shared_bump_constants',compliant=False)
        yes(fixed,'fixed_CA_CQ_directed_bounds_certified','fixed original directed inverse/quadratic constants')
        yes(fixed,'fixed_CS_analytic_bound_certified_under_stated_inputs','universal original CS estimates, with current hypotheses below')
        if fixed['norm']!='sum_i(sup_Z |h_i|+sup_Z |h_i_Z|), Z in [-1,1]':
            raise ValueError('Fixed coefficient Banach norm differs')
        expected=dict(delta_upper='.001',epsilon0_upper='1e-6',Pstar_lower=1,
            coefficient_smallness='CS*||h||_1<=.01',entry_average='||Mz/R-4Z||_1<=2epsilon0',
            entry_pressure='||P||_1<=(Kp+100)Pstar^2',reference_velocity='u=Am*x^.1,V=4Z before correction')
        if fixed['CS_proof_hypotheses']!=expected or fixed['fixed_bump_radius']!='1/40' or fixed['centers']!=['5/4','3/2','7/4']:
            raise ValueError('Universal original patch proof hypotheses changed')
        admit=load('five_defect_admission')
        if (admit['implicit_source_sha256']!=self.family['implicit_source_sha256'] or admit['datum_enclosure_sha256']!=self.family['datum_enclosure_sha256']
            or admit['actual_five_defect_family_sha256']!=current['actual_five_defect_family_sha256']):
            raise ValueError('Current exact transported functional defects differ from analytic class')
        yes(admit,'complete_10_19_actual_functional_defect_test_certified','same original exact five defect functional norm')
        self.theorem=patch_source_theorem();self.service.bind_hashes(self.theorem['input_hashes'])
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
        self.proof=self.bounds()

    def bounds(self):
        c=self.ctx;read=lambda v:packets.interval(c,v)
        lower=lambda v:c.mpf(packets.recovery.endpoints(v)[0])
        upper=lambda v:c.mpf(packets.recovery.endpoints(v)[1])
        fixed=self.rows['shared_bump_constants'];admit=self.rows['five_defect_admission']
        report=self.rows['actual_feedback_moment_patch'];join=self.restore.reshape.rows['reference_join_bounds']
        norm=self.restore.reshape.inner.rows['physical_norm_family']
        eps=read(self.restore.reshape.inner.rows['K1_ledger']['epsilon0']);delta=c.mpf(self.service.data['delta'])
        e=read(admit['complete_actual_functional_e_upper'])
        t=read(report['inherited_coefficient_C1_t_upper'])
        CS=c.mpf(fixed['CS']);B=read(fixed['bump_sup_bound']);D=read(fixed['bump_first_derivative_sup_bound'])
        KN=c.mpf(join['KN']);Kp=c.mpf(join['pressure_Kp']);CSt=upper(CS*t)
        if fixed['KN']!=join['KN']:raise ValueError('Fixed/current KN normalization differs')
        inherited_N=read(join['N_after_Rz_over_Pstar_squared_upper'])
        if packets.recovery.endpoints(inherited_N)!=packets.recovery.endpoints(KN):
            raise ValueError('Same checked original KN transport bound at the Rm inlet required')
        if packets.recovery.endpoints(read(report['inherited_correlated_C1_e_upper']))!=packets.recovery.endpoints(e):
            raise ValueError('Current true functional source C1 norm differs')
        if packets.recovery.endpoints(t)[0]<packets.recovery.endpoints(2*fixed['CA']*e)[1]:
            raise ValueError('Stored current t does not bound the exact analytic 2CAe')
        # The original report used its own interval precision. Its CS product
        # must cover this recomputation, rather than equal its rounding bits.
        old=read(report['CS_times_C1_coefficient_bound']);product=CS*t
        if (packets.recovery.endpoints(old)[0]>packets.recovery.endpoints(product)[0]
            or packets.recovery.endpoints(old)[1]<packets.recovery.endpoints(product)[1]):
            raise ValueError('Stored current CS*t does not cover the same true coefficient bound')
        gates=dict(actual_delta_positive=delta,actual_delta_below_fixed_cap=c.mpf('.001')-delta,
            epsilon_below_fixed_cap=c.mpf('1e-6')-eps,
            actual_Pstar_above1_log_margin=self.service.saved('actual_patch').algebra.logs[1]/2,
            actual_Rm_above16_log_margin=read(norm['radius_log_margins']['Rm16_log_margin']),
            current_Banach_contraction_margin=c.mpf('.5')-4*fixed['CA']**2*fixed['CQ']*e,
            current_coefficient_t_star_margin=read(fixed['t_star'])-t,
            current_CS_smallness_margin=c.mpf('.01')-CSt,
            current_inherited_pressure_budget=Kp+100-read(join['normalized_join_pressure_C1_upper']),
            current_Rm_angular_Q_half_margin=self.restore.reshape.proof['angular_Q_after_one_log_unit_lower']-c.mpf('.5'),
            fixed_velocity_and_Vy_coefficient_budget=CS-2*D,
            fixed_angular_SQ_coefficient_budget=CS-read(fixed['elementary_coefficients']['SQ_source_error']))
        if any(packets.recovery.endpoints(v)[0]<=0 for v in gates.values()):
            raise ArithmeticError('Actual patch universal proof hypothesis unresolved')
        # Whole x[1,e] estimates follow from true C1 h, disjoint compact bumps,
        # positive partial unit mass and the same Rm histories. Quiet gaps keep
        # those partial histories. No independent Taylor-box norm is substituted.
        denominator=lower(1-B*t)
        umin=lower(c.exp(c.mpf('-.6'))*denominator/2)
        umax=upper(c.exp(c.mpf('-.6'))*(c.exp(c.mpf('.1'))+B*t))
        aerr=upper(read(fixed['elementary_coefficients']['a_error'])*t)
        amin,amax=c.mpf('.7'),c.mpf('.9')
        Vy=upper(2*D*t);V=upper(4+B*t)
        W=upper(3+4*delta+2*eps+CSt)
        # Original 9.32, now with the exact patch primitives and arbitrary Vy.
        SQ=lower(c.mpf('1.8')-delta/2-10*eps-CSt)
        # This unpatched join bound covers the full reference continuation to
        # Rh, not just its Rm inlet. Add the actual partial bump pressure change.
        pressure=upper(read(join['normalized_join_pressure_C1_upper'])+3*t)
        # Full signed N equation 9.35. Every energy/pressure and meridional
        # source term is bounded; N is never chosen from this norm cover.
        Nsource=upper(umax**2+(3+2*delta)*pressure+W*Vy+(1+delta)*(1+2*V)*V/2+(c.mpf('.5')+V)*V)
        Nnorm=2*KN;J=128*KN
        bS=upper(2*Vy/umin);bw=upper(512*KN*CSt)
        kappa=upper(amax+bS*bS/amin)
        cone=lower(8*(amin-bw)-2*amax)
        Hmargin=lower(cone/amax)
        margins=dict(positive_original_H_denominator=denominator-c.mpf('.99'),
            normalized_velocity_above_one_eighth=umin-c.mpf('.125'),normalized_velocity_below1=1-umax,
            actual_shear_error_below_point1=c.mpf('.1')-aerr,
            current_angular_SQ_above_seven_fifths=SQ-c.mpf('1.4'),
            current_angular_Q_half_crossing=SQ-(2-amin/2)/2,
            full_N_source_below_two_KN=2*KN-Nsource,
            original_signed_bw_below_point1=c.mpf('.1')-bw,
            original_kappa_below1=1-kappa,original_full_signed_cone_margin=cone,
            original_full_H0_minus2_lower=Hmargin)
        if any(packets.recovery.endpoints(v)[0]<=0 for v in margins.values()):
            raise ArithmeticError('Current whole original patch relaxed input failed')
        return dict(physical_domain='Rm<=R<=Rh=e*Rm=exp(-5)*Rref, x=R/Rm in[1,e], all Z[-1,1]',
            source_charts=['actual_patch'],actual_current_source_hypothesis_gates=gates,
            positive_directed_relaxed_margins=margins,current_true_source_C1_e_upper=e,
            current_unique_coefficient_C1_t_upper=t,current_CS_times_C1_coefficient_upper=CSt,
            original_normalized_swirl_lower=umin,original_normalized_swirl_upper=umax,
            original_angular_a_bounds=[amin,amax],original_signed_Vy_absolute_upper=Vy,
            original_signed_b_times_Pstar_absolute_upper=bS,current_absolute_pressure_C1_over_Pstar_squared_upper=pressure,
            inherited_Rm_N_over_Pstar_squared_upper=inherited_N,
            inherited_Rm_mean_C1_error_upper=2*eps,whole_original_W_absolute_upper=W,
            full_signed_N_source_over_Pstar_squared_upper=Nsource,whole_signed_N_over_Pstar_squared_upper=Nnorm,
            whole_signed_J_absolute_upper=J,original_signed_bw_absolute_upper=bw,
            original_kappa_upper=kappa,original_p1_uniform_lower=c.mpf(8),original_full_H0_minus2_uniform_lower=Hmargin,
            actual_entry_average_proof='same current positive mean transport of V110 and Mz110 then 4Z; summed C1 error<=2epsilon0; patch partial unit mass adds <=t',
            actual_Rm_N_source_attachment='checked restore receipt consumes original arbitrary-V_y five-ODE stress identities, same current pressure/N norm family and exact Rm source/P0 function join; its KN bound after Rz includes Rm',
            actual_pressure_proof='same original P0+Mp and exact inherited five histories; unpatched full reference-join C1 bound through Rh plus actual partial patch pressure change C1<=3t*Pstar^2',
            angular_Q_barrier='Q_y+(2-a/2)Q=SQ; actual Rm Q>=1/2, SQ>7/5 and a>=.7 preserve Q>=1/2 throughout supports and gaps',
            actual_p1_lower_proof='p1=R*Q/L>=Rm/2>=8; actual Rm>=16, 0<L<=1',
            full_N_barrier='N_y+N=Z*u^2+P_operator(P)-W*Vy-(1+delta)*(1-2ZV)*V/2-Hv*V_Z; same Rm |N|<=KN*S^2 and whole source<2KN*S^2 imply |N|<=2KN*S^2',
            signed_bw_proof='u/S>1/8 => |J=N/u^2|<=128KN; Q>=1/2 and |Vy|<=CS*t => |bw|=|2J*Vy/Q|<=512KN*CS*t',
            full_generic_relaxed_input='a in[.7,.9], kappa<1; p1(a-bw)-2a>=8(.7-|bw|)-1.8>0, H0=p1(a-bw)/a>2',
            conservative_full_signed_bw_bound_not_copied_from_legacy_256_estimate=True,
            full_signed_p2_pressure_energy_and_meridional_terms_retained=True,
            all_three_disjoint_support_groups_and_four_quiet_gaps_covered=True,
            quiet_gap_velocities_do_not_reset_current_partial_histories=True,
            terminal_functional_closure_only_at_or_after_last_support_edge='71/40',
            p1_p2_whole_path_norm_bounds_certified=False,strict_completed_tensor_cone_new_regions_admitted=0,
            source_radius_amplitude_and_width_not_materialized=True,**dict.fromkeys(OPEN,False))

    def query(self,x,Z=(-1,1)):
        c=self.ctx;v=c.mpf(x);z=c.mpf(Z)
        vl,vh=packets.recovery.endpoints(v);zl,zh=packets.recovery.endpoints(z)
        if (not all(mp.isfinite(a) for a in (vl,vh,zl,zh)) or vl<1
            or vh>packets.recovery.endpoints(c.exp(1))[1] or zl<-1 or zh>1):
            raise ValueError('Actual original patch x in[1,e], Z in[-1,1] required')
        return dict(chart='actual_patch',x=v,Z=z,source_family=self.family,analytic_source_certificate=self.proof,
            source_function_domain_not_saved_point_values=True,**{GATE:True},**dict.fromkeys(OPEN,False))


def run():
    owner=CurrentPatchRelaxedInputs();c=owner.ctx
    breaks=[c.mpf(1),*[c.mpf(i)/40 for i in EDGES],c.exp(1)]
    examples={str(i):owner.query((packets.recovery.endpoints(lo)[0],packets.recovery.endpoints(hi)[1]))
        for i,(lo,hi) in enumerate(zip(breaks,breaks[1:]))}
    result=dict(source_family=owner.family,current_source_function_attachment_conditions=owner.conditions,
        original_current_full_signed_patch_source_theorem=owner.theorem,
        current_whole_Rm_Rh_relaxed_input_and_bounds=owner.proof,analytic_support_and_quiet_gap_subbox_examples=examples,
        **{GATE:True},**dict.fromkeys(OPEN,False),strict_completed_tensor_cone_new_regions_admitted=0,
        source_graph_ancestor_constructors_called=False,
        scope='Current original Rm..Rh implicit five-moment patch relaxed input, all supports/partial-history gaps. No changed loop, whole norm scales, new repair/N, strict tensor or recursion admission.',
        input_hashes=owner.service.hashes)
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual current whole Rm..Rh implicit patch relaxed input PASS',flush=True)
    return result


if __name__=='__main__':run()

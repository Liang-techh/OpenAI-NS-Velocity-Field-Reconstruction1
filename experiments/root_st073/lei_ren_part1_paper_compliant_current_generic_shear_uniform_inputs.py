"""Actual whole generic-loop input margins and conservative logarithmic scales.

Original relaxed input includes the complete signed strong branch when
kappa>2. It is distinct from a strict closed-region tensor admission.
No microscopic lower bound, radius, amplitude, eta or inverse is resolved.
"""
import ast
import json
from pathlib import Path
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_loop_domain as domain

source=domain.source;bounds=source.bounds;packets=source.packets
HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
NAME=PREFIX+'current_generic_shear_uniform_inputs.json'
RECEIPT=PREFIX+'current_generic_shear_uniform_inputs_check.json'
GATE='current_original_whole_generic_input_margins_and_log_scales_certified'
OPEN=source.OPEN
REQUIRED={
 'current_inner_relaxed_inputs':'current_original_inner_Ra_R110_relaxed_generic_input_certified',
 'current_reshape_relaxed_inputs':'current_original_R110_Rz_relaxed_generic_input_certified',
 'current_restore_relaxed_inputs':'current_original_Rz_Rm_relaxed_generic_input_certified',
 'current_patch_relaxed_inputs':'current_original_Rm_Rh_relaxed_generic_input_certified',
 'current_O2_reference_slope_relaxed_cone':'current_original_Rh_O2_slope_relaxed_input_cone_certified',
 'current_O2_axial_relaxed_cone':'current_original_O2_axial_relaxed_input_cone_certified',
 'current_O2_relaxed_buffer_cone':'current_original_O2_buffer_relaxed_input_cone_certified',
 'current_O2_modified_taper_cone':'current_modified_O2_open_taper_signed_two_vector_cone_certified',
 'current_O3_transition_direction':'current_whole_O3_variable_transition_direction_bound_certified',
 'current_O3_power_cone':'current_whole_O3_power_signed_two_vector_cone_certified'}


def exact_normalized_source_theorem():
    """I/F and T/F differ by the original shear, including retained axial stress."""
    R,S,U,C=s.symbols('R S scalar_U axial_shape_C',positive=True)
    a,b,theta,axial=s.symbols('a b complete_inertial_theta_shape full_axial_shape',real=True)
    F=S*U*C/s.sqrt(2*R)
    It=S*s.sqrt(R/2)*U*theta;Iz=S**2*s.sqrt(R/2)*U**2*axial
    p1=It/F;p2=Iz/F;Ttheta=It-a*F;Tz=Iz+b*F
    tn,zn=Ttheta/F,Tz/F;kappa=a+b*b/a
    checks={}
    def zero(name,left,right):
        if s.cancel(left-right)!=0:raise ArithmeticError('Whole original input identity failed: '+name)
        checks[name]=True
    zero('same_inertial_p1',p1,R*theta/C)
    zero('same_full_retained_inertial_p2',p2,S*R*U*axial/C)
    zero('same_full_stress_theta',tn,p1-a)
    zero('same_full_stress_axial',zn,p2+b)
    H0=p1-b*p2/a;D=tn-b*zn/a;J=zn+b*tn/a
    zero('H0_minus2_is_D_plus_kappa_minus2',H0-2,D+kappa-2)
    zero('same_full_generic_transverse',J,p2+b*p1/a)
    theta_full=theta-a*C/R
    zero('b0_source_theta_after_radial_shear',tn,R*theta_full/C)
    zero('b0_full_stress_ratio_not_zero',Tz.subs(b,0)/Ttheta,S*U*axial/theta_full)
    zero('b0_strong_quadratic_relative_to_full_theta',
        (2*D**2-(kappa-2)*J**2).subs(b,0)/tn**2,
        2-(a-2)*(p2/tn)**2)
    asts=packets.recovery.numeric.transport.SourceAST()
    for method,target in (('actual_variable_transition_theorem','theta'),
        ('actual_variable_transition_theorem','full_axial'),
        ('actual_variable_transition_theorem','Pabs'),
        ('whole_current_transition_bounds','energy_log'),
        ('whole_current_transition_bounds','pressure_log'),
        ('whole_current_transition_bounds','direction')):
        asts.expression('current_O3_transition_direction_operator',method,target)
    asts.expression('current_O3_transition_direction_operator','whole_current_transition_bounds','direction',
        wanted='2*mu*(wE+wP)**2')
    for target in ('raw','theta','Pabs'):
        asts.expression('current_O2_modified_taper_cone','exact_theorem',target)
    return dict(passed=True,exact_identities=checks,original_complete_source_AST_bindings=asts.bindings,
        transition_coordinate='t=log(R/Rd) in[0,1], Rw=Rd*exp(1)',
        power_coordinate='t=log(R/Rw)=Tw*phase; first unit is phase[0,1/Tw]',
        scalar_U_not_the_full_Utheta_over_S='Utheta/S=scalar_U*C, C=1/(1+Z^2)',
        theta_floor_bounds_full_Ttheta_not_relabelled_as_exact_inertial_p1=True,
        ratio_budget_applies_to_Tz_over_Ttheta_with_radial_shear_in_denominator=True,
        no_axial_stress_dropped_when_local_Uz_zero=True,
        whole_closed_transition_strict_cone_not_claimed=True,input_hashes=asts.hashes)


def scale_formula_theorem():
    filename=PREFIX+'current_generic_shear_loop.py'
    tree=ast.parse((HERE/filename).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='GenericLoopScales')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    wanted={
        'self.d_star':'self.margin_min/4',
        'self.q_star':'c.sqrt(3/(2*self.a_min))',
        'self.B_star':'self.t0_abs_max+2*self.q_star+4*self.p2_abs_max*self.q_star**2/self.d_star',
        'self.J_star':'self.p1_abs_max*self.B_star+self.p2_abs_max',
        'self.eta':'min(c.mpf(1)/2,self.margin_min/8,self.margin_min**2/(8*(1+self.J_star**2)),self.boundary_kappa_excess_min/2)'}
    checks={}
    for target,expr in wanted.items():
        found=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        if len(found)!=1 or ast.dump(found[0])!=ast.dump(ast.parse(expr,mode='eval').body):
            raise ValueError('Original Section11 scale recipe changed: '+target)
        checks[target]=True
    return dict(passed=True,original_scale_recipe_AST_bindings=checks,
        conservative_upper_q_B_J_and_lower_eta_permitted=True,
        input_hashes={filename:sha(filename)})


class CurrentUniformInputs:
    def __init__(self,provider=None):
        self.domain=provider if provider is not None else domain.CurrentLoopDomain()
        self.service=self.domain.service;self.ctx=self.domain.ctx;self.family=self.domain.family;self.rows={}
        for stem,gate in REQUIRED.items():
            name=PREFIX+stem+'_check.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(gate):raise ValueError('Checked original margin prerequisite required: '+stem)
            family=receipt.get('source_family') or {k:receipt[k] for k in packets.FAMILY_KEYS}
            if family!=self.family:raise ValueError('Whole input source/pressure family differs: '+stem)
            self.service.bind_hashes(receipt['input_hashes']);self.service.bind_hashes({name:sha(name)})
            producer=PREFIX+stem+'.json';row=json.loads((HERE/producer).read_bytes())
            self.service.bind_hashes({producer:sha(producer)});self.rows[stem]=row
        receipt=json.loads((HERE/domain.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[domain.GATE] or receipt['source_family']!=self.family:
            raise ValueError('Checked actual two-sided loop domain required')
        self.service.bind_hashes(receipt['input_hashes']);self.service.bind_hashes({domain.RECEIPT:sha(domain.RECEIPT)})
        self.geometry=self.domain.domain();self.theorem=exact_normalized_source_theorem();self.scale_theorem=scale_formula_theorem()
        self.service.bind_hashes(self.theorem['input_hashes']);self.service.bind_hashes(self.scale_theorem['input_hashes'])
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
        self.previous=self.domain.provider.previous
        self.O3=json.loads((HERE/source.NAME).read_bytes())

    def margins(self):
        c=self.ctx;read=lambda v:packets.interval(c,v);endpoints=packets.recovery.endpoints
        lower=lambda v:c.mpf(endpoints(v)[0]);records={};conditions={}
        def logpositive(value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Positive current actual H0 margin required')
            return lower(c.ln(value))
        def admit(key,charts,logH,proof,branch):
            if not all(__import__('mpmath').isfinite(v) for v in endpoints(logH)):
                raise ArithmeticError('Finite actual H0 log lower required')
            records[key]=dict(source_charts=charts,log_H0_minus2_positive_lower=lower(logH),
                original_source_proof=proof,original_relaxed_branch=branch,
                full_pressure_energy_and_original_moments_retained=True)
        inner=self.rows['current_inner_relaxed_inputs']['current_whole_inner_relaxed_input_and_bounds']
        admit('inner',bounds.INNER,logpositive(read(inner['positive_directed_margins']['whole_inner_H0_minus2_lower'])),inner,
            'Original exact global exit source theorem supplies D,Q>0 where kappa>2; weak branch H0>2 throughout')
        shape=self.rows['current_reshape_relaxed_inputs']['current_whole_R110_Rz_relaxed_input_and_bounds']
        if not shape['whole_path_H0_minus2_strictly_above1']:raise ValueError('Same original whole angular barrier required')
        admit('reshape_reference',('reshape','inner_reference'),c.mpf(0),shape,'b=0, kappa=a<=.9<2, H0=p1>3')
        for stem,key,charts in (
            ('current_restore_relaxed_inputs','current_whole_Rz_Rm_relaxed_input_and_bounds',('axial_restore','restore_buffer')),
            ('current_patch_relaxed_inputs','current_whole_Rm_Rh_relaxed_input_and_bounds',('actual_patch',))):
            proof=self.rows[stem][key]
            admit(stem,charts,logpositive(read(proof['original_full_H0_minus2_uniform_lower'])),proof,'Original full signed H0 lower, kappa<1')
        ref=self.rows['current_O2_reference_slope_relaxed_cone']['whole_original_reference_slope_relaxed_cone']
        logp=read(ref['Itheta_over_F_log_lower'])
        if endpoints(logp-c.ln(3))[0]<=0:raise ArithmeticError('Reference/slope H0 lower must exceed3')
        conditions['reference_slope_log_p1_above_log3']=logp-c.ln(3)
        admit('O2_reference_slope',('Rh_reference','O2_slope'),c.mpf(0),ref,'b=0, kappa<=2; p1>3 gives H0-2>1')
        axial=self.rows['current_O2_axial_relaxed_cone'];base=axial['whole_original_axial_baseline'];proof=axial['whole_original_axial_relaxed_cone']
        axiallog=(self.service.data['logRref']+1+c.ln(read(base['canonical_theta_positive_reserve']))
            +c.ln(read(proof['full_D_over_theta_lower'])))
        if endpoints(axiallog)[0]<=0:raise ArithmeticError('Actual whole axial H0 margin must exceed1')
        conditions['axial_log_D_positive_lower']=lower(axiallog)
        admit('O2_axial',('O2_axial',),c.mpf(0),dict(baseline=base,signed_cone=proof),
            'H0-2=D+bs^2/2>=D; full signed D/theta lower; strict interior/nonzero-Z, relaxed flat edges/midplane')
        logRd=self.geometry['exact_log_Rw']-1
        for key,stem,basekey,start,selectors in (
            ('O2_buffer_first','current_O2_relaxed_buffer_cone','whole_original_buffer_baseline',-11,(0,9)),
            ('O2_buffer_last','current_O2_modified_taper_cone','original_whole_closed_O2_taper_baseline_bounds',-2,(9,11))):
            baseline=self.rows[stem][basekey]
            if not baseline.get('full_original_pressure_energy_moments_and_radial_sectors_retained',
                baseline.get('same_original_O2_pressure_energy_moments_and_all_stress_sectors',False)):
                raise ValueError('Use original buffer baseline, not modified taper output')
            logtheta=logRd+start+c.ln(read(baseline['canonical_theta_lower']))
            if endpoints(logtheta)[0]<=0:raise ArithmeticError('Actual pure buffer H0 margin must exceed1')
            conditions[key+'_full_Ttheta_over_F_log_lower']=lower(logtheta)
            admit(key,('O2_buffer',),c.mpf(0),dict(original_baseline=baseline,original_selector_domain=selectors,
                shared_offset_domain=(start,selectors[1]-11),modified_profile_and_frequency_ignored=True),
                'Original a=2,b=0,kappa=2; H0-2=Ttheta/F>1, exact same pre.axial source on both pieces')
        transition=self.rows['current_O3_transition_direction']['current_whole_variable_transition_bounds']
        transitionlog=read(transition['source_correlated_periodic_shear_loop']['original_stress_over_F_log_lower'])
        directional=read(transition['full_directional_expression_upper'])
        if endpoints(transitionlog)[0]<=0 or endpoints(2-directional)[0]<=0:
            raise ArithmeticError('Original full variable transition direction and stress reserve required')
        conditions['O3_transition_full_Ttheta_over_F_log_lower']=transitionlog
        conditions['O3_transition_full_signed_Q_over_Ttheta_squared_lower']=lower(2-directional)
        admit('O3_transition',('O3_slope_mu',),c.mpf(0),transition,
            'D=Ttheta/F>1, H0-2=D+2mu*sigma; Q/D^2>=2-direction>0. kappa=2 at0 relaxed, kappa>2 on(0,1] strict')
        power=self.rows['current_O3_power_cone']['current_whole_O3_power_correlated_bounds']
        powerlog=self.geometry['exact_log_Rw']+c.ln(read(power['theta_lower_constants']['theta_floor']))
        if endpoints(powerlog)[0]<=0:raise ArithmeticError('Actual first power unit H0 margin must exceed1')
        conditions['O3_power_full_Ttheta_over_F_log_lower']=lower(powerlog)
        conditions['O3_power_full_signed_Q_over_Ttheta_squared_lower']=lower(2-read(power['full_directional_expression_upper']))
        if endpoints(conditions['O3_power_full_signed_Q_over_Ttheta_squared_lower'])[0]<=0:
            raise ArithmeticError('Same whole original power signed cone required')
        admit('O3_power',('O3_power',),c.mpf(0),power,'b=0, kappa-2=2mu>0; full signed whole power D,Q bound')
        logm=c.mpf(min(endpoints(v['log_H0_minus2_positive_lower'])[0] for v in records.values()))
        return dict(whole_modification_domain='r_minus<=R<=r_plus, all Z[-1,1], original source only',
            chart_groups=records,source_margin_attachment_conditions=conditions,
            uniform_log_H0_minus2_positive_lower=logm,
            full11_unit_O2_buffer_composed_at_selector9_by_same_original_function=True,
            original_weak_and_strong_signed_branches_on_entire_modification_domain_certified=True,
            original_closed_transition_strict_cone_claimed=False,
            new_global_or_modified_cone_regions_admitted=False)

    def log_scales(self,margin):
        c=self.ctx;endpoints=packets.recovery.endpoints;read=lambda v:packets.interval(c,v)
        lower=lambda v:c.mpf(endpoints(v)[0]);upper=lambda v:c.mpf(endpoints(v)[1])
        prior=self.previous['current_original_source_log_bound_charts'];outer=self.O3['original_O3_quotient_log_norms']
        positive=self.previous['positive_noncore_source_denominator_theorem']['source_charts']
        loga=min(endpoints(read(v['log_actual_a_positive_lower']))[0] for v in positive.values())
        loga=min([loga,*[endpoints(read(v['actual_positive_denominator_theorem']['log_actual_a_positive_lower']))[0] for v in outer.values()]])
        loga=c.mpf(loga);maxima={};allcharts=list(positive)+list(outer)
        for key in ('t0','p1','p2'):
            values=[]
            for chart in allcharts:
                q=(prior[chart]['admitted_original_quotient_log_norms'] if chart in prior else outer[chart]['admitted_original_quotient_log_norms'])[key]['y0_Z0']
                if not q['exact_zero']:values.append(read(q['log_absolute_upper']))
            maxima[key]=None if not values else c.mpf(max(endpoints(v)[1] for v in values))
        logm=margin['uniform_log_H0_minus2_positive_lower'];logd=lower(logm-c.ln(4))
        logq=upper((c.ln(c.mpf(3)/2)-loga)/2)
        terms=[bounds.LogUpper(c,maxima['t0']),bounds.LogUpper(c,logq+c.ln(2)),
            bounds.LogUpper(c,None if maxima['p2'] is None else c.ln(4)+maxima['p2']+2*logq-logd)]
        logB=bounds.LogUpper.add(c,terms).log
        logJ=bounds.LogUpper.add(c,[bounds.LogUpper(c,None if maxima['p1'] is None else maxima['p1']+logB),
            bounds.LogUpper(c,maxima['p2'])]).log
        oneplusJ2=bounds.LogUpper.add(c,[bounds.LogUpper(c,c.mpf(0)),bounds.LogUpper(c,None if logJ is None else 2*logJ)]).log
        eta_constraints=dict(half=lower(-c.ln(2)),Hmargin_over8=lower(logm-c.ln(8)),
            quadratic=lower(2*logm-c.ln(8)-oneplusJ2),
            both_boundary_collars=read(self.geometry['allowed_eta_log_upper_for_automatic_constant_edges']))
        logeta=c.mpf(min(endpoints(v)[0] for v in eta_constraints.values()))
        return dict(whole_source_chart_inventory=allcharts,source_quotient_order='y0_Z0 for scale amplitudes; all exported y2/Z1 rows remain attached',
            logarithmic_selected_positive_lower_constants=dict(a_min=loga,margin_min=logm,d_star=logd,
                boundary_kappa_excess_min=self.geometry['both_boundary_kappa_excess_log_lower']),
            logarithmic_conservative_upper_constants=dict(t0_abs_max=maxima['t0'],p1_abs_max=maxima['p1'],p2_abs_max=maxima['p2'],
                q_star=logq,B_star=logB,J_star=logJ,one_plus_J_star_squared=oneplusJ2),
            selected_positive_eta_log=logeta,eta_required_log_upper_constraints=eta_constraints,
            constant_interpretation='Each selected positive constant is exp(its saved log); None upper is exact zero. No source exponentiation is performed.',
            eta_and_two_plus_eta_stay_formal_separate_sources=True,
            original_Section11_scale_recipe=self.scale_theorem,
            complete_modification_source_margin_and_norms_consumed=True,
            old_scalar_fixture_constants_not_consumed=True,
            phase_held_inverse_or_loop_jets_constructed=False,new_common_finite_N_selected=False)

    def run(self):
        margins=self.margins();scales=self.log_scales(margins)
        result=dict(source_family=self.family,exact_full_signed_source_normalization_theorem=self.theorem,
            whole_actual_original_generic_input_margin=margins,current_actual_logarithmic_loop_scales=scales,
            current_chosen_domain=self.geometry,
            **{GATE:True},**dict.fromkeys(OPEN,False),whole_modification_original_relaxed_input_certified=True,
            current_whole_source_loop_scales_in_log_form_certified=True,
            phase_held_loop_primitive_derivative_bounds_certified=False,
            original_outer_admission_twice_Rc_through_Rb_certified=False,
            source_graph_ancestor_constructors_called=False,
            scope='Actual original whole modification input H0 margin, correct weak/strong signed branches, full O2 buffer split, complete source quotient norm maxima and conservative Section11 logarithmic constants. No physical loop/inverse jets, changed histories, new repair/N, global cone or true recursion admission.',
            input_hashes=self.service.hashes)
        (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
        print('Actual whole original generic H0 margin and logarithmic loop scales PASS',flush=True)
        return result


def run():return CurrentUniformInputs().run()


if __name__=='__main__':run()

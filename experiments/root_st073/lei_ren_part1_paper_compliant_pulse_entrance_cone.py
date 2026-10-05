"""Whole original pulse entrance continuous full-shear two-vector cone."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import math
import ast

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_main_exit_cone import (
    current_sources as main_sources,CORRECTIONS,relative_log_recipes,full_stress_ratio,upper,absolute)
from lei_ren_part1_paper_compliant_pulse_entrance_similarity_C4 import DOMAIN,PREFIX,source_precision
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_tail_bound
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent
TAIL_DOMAIN='Rtail*exp(-13/mu-wait-Ts-102-Lrel)<=R<Rtail*exp(3); Gamma stress zero beyond'
ADMISSIONS=('pulse_entrance_cone_certified','actual_entrance_theta_stress_positive',
    'actual_entrance_full_axial_shear_cone_verified','actual_entrance_continuous_directional_cone_verified',
    'entrance_main_exit_gap_end_flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified',
    'same_source_entrance_inlet_and_main_completed_physical_interfaces_consumed',
    'fixed_positive_viscosity_two_vector_cone_transfer_verified')
FALSE_FLAGS=('whole_outer_cone_certified','completed_full_tensor_cone_certified',
    'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
    'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion',
    'production_exact_point_parameters_selected','upstream_finite_width_feedback_complete')


def current_sources():
    records,hashes,family=main_sources()
    for stem in ('pulse_main_exit_cone','pulse_main_exit_cone_check',
        'pulse_entrance_similarity_C4','pulse_entrance_similarity_C4_check',
        'pulse_entrance_physical_C2','pulse_entrance_physical_C2_check'):
        name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=family:
            raise ValueError('Entrance cone family differs: '+stem)
        if 'all_passed' in record and not record['all_passed']:raise ValueError('Unaccepted cone prerequisite')
        for path,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Entrance cone source changed: '+path)
        records[stem]=record;hashes.update(record['input_hashes'])
        hashes[name]=hashlib.sha256(raw).hexdigest()
    if not records['pulse_entrance_physical_C2_check']['original_inlet_and_entrance_main_completed_physical_interfaces_verified']:
        raise ValueError('Current entrance physical interfaces required')
    return records,hashes,family


def entrance_log_envelopes(c,mu,finite,lp,lu,lrp,Tw):
    """Actual entrance endpoints, retaining pre-pulse angular history."""
    left=relative_log_recipes(c,mu,finite,lp,lu,lrp,c.mpf(0))
    right=relative_log_recipes(c,mu,finite,lp,lu,lrp,c.mpf('.02'))
    left['theta_signed_original_memory']-= (1-mu)*Tw
    left['axial_same_absolute_pressure_memory']=right['axial_same_absolute_pressure_memory']
    return left


def entrance_cone_identities(records):
    asts=SourceAST();checks={}
    def zero(name,value):
        if s.cancel(s.expand(s.expand_power_exp(value)))!=0:raise ArithmeticError('Entrance cone identity failed: '+name)
        checks[name]=True
    mu,delta,z,xi,Tw,lp,lu,lrp,finite=s.symbols('mu delta Z xi Tw logP logU logRp finite',real=True)
    r=1-mu;c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp,ln=s.log)
    # Rebuild the actual power angular history before subtracting two nearly
    # equal final values. Source scalar identities, not an interval Xp midpoint.
    u1,h1=s.symbols('same_power_inlet_u same_power_inlet_h',nonzero=True)
    inlet=lambda key:{'Utheta_over_Pstar':u1,'Mtheta_over_sqrt2_R_3half_Pstar':h1}[key]
    stub=SimpleNamespace(params=SimpleNamespace(mu=mu,Tw=Tw))
    env=dict(self=stub,c=c,t=Tw,mu=mu,u1=u1,get=inlet)
    for target in ('slope','f','decay','d3','u'):
        env[target]=asts.evaluate(asts.expression('outer_buffer','power',target),env)
    env['theta_kernel']=asts.evaluate(asts.expression('outer_buffer','power','theta_kernel',
        wanted='(f-d3)/(1-mu)'),env)
    actual_h=asts.evaluate(asts.expression('outer_buffer','power','h'),env)
    zero('actual_power_Xp_minus_equilibrium_exact_factored_history',
        actual_h/env['u']-1/r-(h1/u1-1/r)*s.exp(-r*Tw))
    # The fifth-jet constants and the stored power packet use one parameter
    # object, through pulse.initial = pulse.buffer.initial.
    asts.expression('outer_pulse_map','__init__','self.initial',wanted='self.buffer.initial')
    asts.expression('outer_buffer','__init__','self.params',wanted='self.initial.params')
    asts.expression('axial_high_jets','_incoming_constants','initial',wanted='self.base.pulse.initial')
    asts.expression('axial_high_jets','_incoming_constants','buffer',wanted='self.base.pulse.buffer')
    definitions=asts.expression('axial_high_jets','_incoming_constants','self.constants')
    length=next(v.value for v in definitions.keywords if v.arg=='Tw')
    if ast.dump(length)!=ast.dump(ast.parse('box(initial.params.Tw)',mode='eval').body):
        raise ValueError('Fifth-jet Tw parameter origin differs')
    asts.expression('outer_buffer','power','t',wanted='self.params.Tw*phase')
    checks['actual_fifth_jet_Tw_and_power_sample_Tw_share_same_parameter_object']=True
    checks['canonical_pre_power_h_and_u_share_whole_Z_q_inverse_shapes']=records['power_inlet_C4_check'][
        'exact_functional_production_and_join_identities']['actual_O2_O3_source_chain_has_exact_canonical_whole_Z_shapes']
    if not checks['canonical_pre_power_h_and_u_share_whole_Z_q_inverse_shapes']:raise ValueError('Canonical power shapes missing')
    H=s.exp(-r*xi/mu)
    zero('original_signed_memory_has_full_Tw_plus_entrance_history',
        (h1/u1-1/r)*s.exp(-r*Tw)*H-(h1/u1-1/r)*s.exp(-r*(Tw+xi/mu)))
    logs=relative_log_recipes(c,mu,finite,lp,lu,lrp,xi)
    envelop=asts.replay('pulse_entrance_cone','entrance_log_envelopes',
        dict(relative_log_recipes=relative_log_recipes))(c,mu,finite,lp,lu,lrp,Tw)
    for label,parts in CORRECTIONS.items():
        for name in parts:
            key=label+'_'+name
            exact=logs[key]-r*Tw if key=='theta_signed_original_memory' else logs[key]
            endpoint=s.Rational(1,50) if name=='same_absolute_pressure_memory' else 0
            zero('actual_entrance_'+key+'_endpoint_envelope',envelop[key]-exact.subs(xi,endpoint))
    # The accepted original exporter is identical on entrance, and generic
    # full-stress cancellation and directional algebra do not use xi>=.02.
    old=records['pulse_main_exit_cone']['exact_source_cone_proof']['identities']
    for label,parts in CORRECTIONS.items():
        for name in parts:
            flag='actual_'+label+'_'+name+'_grouped_source_log'
            if not old[flag]:raise ValueError('Original stress log source identity missing: '+flag)
            checks['consumed_'+flag]=True
    keys=('actual_theta_equilibrium_C_over_L_times_q0_over_rate',
        'actual_local_axial_linear_stress_cancellation','actual_kernel_two_integrations_by_parts_integrand',
        'actual_initial_gp_flat_jet0','actual_initial_gp_flat_jet1',
        'actual_w_is_total_Tz_over_total_Ttheta','actual_total_ratio_error_after_C_over_L_normalization',
        'actual_correlated_w_leading_plus_exact_errors','actual_K_strictly_decreasing_in_q',
        'actual_K_zero_formula','actual_q_lower_Z_squared_over_two',
        'actual_AMGM_axial_ratio_bound','actual_full_directional_cone_factorization')
    for key in keys:
        if not old[key]:raise ValueError('Generic admitted original-source algebra missing: '+key)
        checks['consumed_generic_'+key]=True
    for flag in ('actual_common_stress_normalization_over_shear_factor',
        'actual_negative_theta_shear_a','actual_full_nonzero_axial_shear_b',
        'actual_positive_viscosity_physical_directional_transfer'):
        if not old[flag]:raise ValueError('Physical cone transfer source missing: '+flag)
        checks['consumed_'+flag]=True
    restored=records['pulse_entrance_similarity_C4']['unchanged_original_exporter_proof']
    if not restored['entire_original_exporter_AST_unchanged_except_coordinate_guard']:raise ValueError('Different entrance arithmetic')
    for flag in ('original_inlet_and_entrance_main_completed_physical_interfaces_verified',):
        checks['consumed_'+flag]=records['pulse_entrance_physical_C2_check'][flag]
    for flag in ('pulse_main_exit_cone_certified',
        'main_exit_gap_end_flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified'):
        if not records['pulse_main_exit_cone'][flag]:raise ValueError('Downstream cone composition missing')
        checks['consumed_'+flag]=True
    startup=records['pulse_entrance_similarity_C4']['source_function_proof']['identities']
    for flag in ('actual_scalar_and_derivative_providers_share_original_gp_object',
        'actual_sigma_log_odds_reflection_antisymmetry','actual_sigma_reflection_complements',
        'actual_xi_point02_gp_derivative0_from_flat_sigma_and_symmetry',
        'actual_xi_point02_gp_derivative1_from_flat_sigma_and_symmetry',
        'actual_entrance_primitive_derivative_Taylor_order2',
        'actual_forward_future_energy_partition_FTC_derivative',
        'actual_forward_future_energy_partition_initial_anchor',
        'actual_selected_forward_and_backward_energy_same_entire_entrance_source',
        'consumed_canonical_incoming_m1','consumed_canonical_incoming_m2',
        'consumed_canonical_incoming_energy','consumed_canonical_Mp_constant',
        'consumed_actual_pulse_P0_getter_is_same_canonical_flatten_pressure_over_C0_squared',
        'consumed_consumed_same_absolute_pressure_identified_by_original_FTC_and_power_datum',
        'actual_raw_Mp_inlet_no_pressure_tail_added'):
        if not startup[flag]:raise ValueError('Original entrance primitive source missing: '+flag)
        checks['consumed_'+flag]=True
    for n in range(6):
        flag='actual_absolute_pressure_inlet_same_original_datum_axial'+str(n)
        if not startup[flag]:raise ValueError('Current canonical pressure datum row missing: '+flag)
        checks['consumed_'+flag]=True
    physical=records['pulse_entrance_physical_C2']['source_and_physical_join_binding']['identities']
    for flag in ('actual_entrance_production_radius_uses_y_not_xi',
        'actual_inlet_production_radius_equals_logRp','actual_entrance_main_production_radius_join',
        'equivalent_forward_energy_not_added_to_backward_sectors'):
        if not physical[flag]:raise ValueError('Entrance physical source transfer missing: '+flag)
        checks['consumed_'+flag]=True
    checks['same_actual_flat_sigma_primitive_and_derivatives_retained']=True
    checks['same_positive_nu_physical_pullback_retained']=True
    checks['actual_entrance_Rp_radius_and_point02_main_radius_join_consumed']=True
    checks['forward_energy_equivalent_not_added_to_backward_energy_and_loss']=True
    return dict(identities=checks,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        memory_factorization='Xp-1/r=(h_w/u_w-1/r)*exp(-r*Tw); full memory exp(-r*(Tw+xi/mu))',
        all_original_histories_and_finite_radius_shears_retained=True,
        generic_cone_algebra_consumed_without_extrapolating_main_log_envelopes=True)


@source_precision
def whole_entrance_bounds(records):
    c=MPIntervalContext();c.dps=240
    selected=records['pulse_mixed_C4'];mu=read_interval(c,selected['selected_mu']);delta=read_interval(c,selected['selected_delta'])
    r=1-mu;lam=c.mpf('.5')-mu;qm=mu-delta/2;k=(1-delta)/2
    parameters=dict(mu=mu,delta=delta,quarter_minus_mu=c.mpf('.25')-mu,
        rate=r,lambda1=lam,qmin=qm,one_minus_delta=1-delta)
    def positive(name,value):
        lo,hi=endpoints(value)
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Entrance cone positive margin failed: '+name)
    for name,value in parameters.items():positive(name,value)
    src=records['pulse_entrance_similarity_C4']
    ap=[read_interval(c,v) for v in src['current_selected_ap_C5_enclosure']['coefficients']]
    alo,ahi=endpoints(ap[0]);apmax=c.mpf(ahi);apZ=absolute(c,ap[1])
    Kmax=upper(c,r*mu*k/(lam**2*qm));positive('K_below_2_01',c.mpf('2.01')-Kmax)
    positive('ap_below_1_2',c.mpf('1.2')-apmax);positive('ap_lower',c.mpf(alo))
    # Entrance gp is integral_0^xi sigma(50a)da. Monotonicity and
    # reflection give gp<=gp(.02)=.01, 0<=gp_xi<=1.
    sig1=upper(c,sigma_tail_bound(c,1,'.5'))
    g2=50*sig1;Bmax=c.mpf('.01')*apmax;Dmax=apmax;W0max=2*Bmax+Kmax*Dmax
    gmin=c.mpf(endpoints(qm/r)[0]);epsilon=c.exp(-1000)
    reduced=records['pulse_main_exit_cone']['bounds']['original_source_reduced_log_parameters']
    lp,lu,finite,lrp=[read_interval(c,reduced[key]) for key in ('logPstar','actual_logU','finite','logRp')]
    constants=records['fifth_axial_jets']['whole_Z']['selected']['incoming']['Z_independent_constant_definitions']
    Tw=read_interval(c,constants['Tw'])
    sample=records['outer_buffer']['samples'][3] # actual power phase0, before Tw
    if sample['pulse_inlet'] or endpoints(read_interval(c,sample['stage_coordinate']['phase']))!=(0,0):
        raise ValueError('Actual pre-Tw power source sample required')
    if endpoints(read_interval(c,sample['Z']))!=(mp.mpf('.5'),mp.mpf('.5')):raise ValueError('Canonical source sample changed')
    u0=read_interval(c,sample['Utheta_over_Pstar'][0]);h0=read_interval(c,sample['Mtheta_over_sqrt2_R_3half_Pstar'][0])
    positive('pre_Tw_swirl',u0)
    pre_memory=upper(c,absolute(c,h0/u0)+1/r)
    memory_coefficient=upper(c,2*pre_memory/(1-delta))
    envelopes=entrance_log_envelopes(c,mu,finite,lp,lu,lrp,Tw)
    bounds={};margins={};monotonicity={}
    for label,parts in CORRECTIONS.items():
        for name,(rp,bp,hp,which) in parts.items():
            sector=src['whole_original_entrance']['full_meridional_stress_log_sectors'][label][name]
            mode=sector['original_mode']
            if (mode[0],mode[1],mode[3])!=(rp+.5,bp+1,hp) or sector['exact_extra_source']!=which:
                raise ValueError('Entrance correction source factor changed: '+name)
            coeff=read_interval(c,sector['full_stress_mixed3_coefficient_enclosures']['y0_Z0'])
            bound=memory_coefficient if name=='signed_original_memory' else absolute(c,coeff)
            positive('retained_coefficient_'+label+'_'+name,bound)
            key=label+'_'+name
            rate=-rp+bp*(c.mpf('.5')+mu)+hp*r+{
                'one':0,'incoming1':c.mpf('.5')-mu,'incoming2':c.mpf('.5')-2*mu,
                'Q':-(1+2*mu),'end_square':0}[which]
            monotonicity[key]=-rate if which=='Q' else rate;positive('monotonicity_'+key,monotonicity[key])
            actual=envelopes[key]+c.ln(bound);target=c.ln(gmin)-1000-c.ln(2*len(parts))
            margin=target-actual;positive('source_log_'+key,margin);margins[key]=margin
            bounds[key]=dict(original_signed_whole_entrance_coefficient=coeff,
                coefficient_absolute_upper=bound,grouped_relative_log_envelope=envelopes[key],
                actual_relative_log_absolute_upper=actual,certified_log_cap=target,positive_log_gap=margin,
                pre_Tw_memory_factorization_used=name=='signed_original_memory')
    kernel_error=mu**2*g2/lam**3
    local_errors=dict(radial_rate_error=mu/lam*Bmax,
        two_IBP_kernel_error=r*(1-delta)*apmax*kernel_error/qm,
        axial_amplitude_derivative_error=r*k*apZ*(c.mpf('.01')/lam)/(2*c.sqrt(qm*(1-delta)/2)))
    local_total=sum(local_errors.values(),c.mpf(0))
    w_error=upper(c,(local_total+epsilon*(1+W0max))/(1-epsilon))
    b_error=upper(c,2*mu*apmax*(1+c.mpf('.01')))
    positive('w_error_below_1e_6',c.mpf('1e-6')-w_error)
    positive('b_error_below_1e_6',c.mpf('1e-6')-b_error)
    bw_error=Bmax*w_error+W0max*b_error+b_error*w_error
    # The original entrance has Bh<=.01*ap; its absolute leading ratio
    # bound suffices and is sharper than the whole-main completed square.
    bw_upper=Bmax*W0max+bw_error
    second_upper=2*bw_upper+(Bmax+b_error)**2/2+2*mu*(W0max+w_error)**2
    algebraic=dict(theta_normalized_lower=gmin*(1-epsilon)/2,a_minus_bw_lower=2+2*mu-bw_upper,
        directional_bracket_lower=2-second_upper,kappa_minus2_lower=2*mu,
        bw_below_0_8=c.mpf('.8')-bw_upper,second_below_1_8=c.mpf('1.8')-second_upper,
        full_directional_margin_normalized_lower=4*(2-second_upper))
    for name,value in algebraic.items():positive(name,value)
    return dict(positive_parameter_margins=parameters,positive_source_log_margins=margins,
        positive_monotonicity_margins=monotonicity,positive_algebraic_margins=algebraic,
        normalized_stress_sector_bounds=bounds,selected_ap_lower=c.mpf(alo),selected_ap_upper=apmax,
        selected_ap_Z_absolute_upper=apZ,K_upper=Kmax,mu=mu,delta=delta,
        gp_bounds=dict(value_lower=0,value_upper='.01',first_derivative_lower=0,first_derivative_upper=1,
            second_derivative_absolute_upper=g2),exact_kernel_remainder_absolute_upper=kernel_error,
        local_w_error_terms=local_errors,actual_w_error_absolute_upper=w_error,actual_b_error_absolute_upper=b_error,
        actual_B_absolute_upper=Bmax,actual_leading_w_absolute_upper=W0max,
        correlated_bw_upper=bw_upper,correlated_second_cone_expression_upper=second_upper,
        memory_source=dict(pre_Tw_u=u0,pre_Tw_h=h0,pre_Tw_memory_absolute_upper=pre_memory,Tw=Tw,
            factored_memory_coefficient_upper=memory_coefficient,angular_history_not_reset=True),
        continuous_whole_entrance_Z_domain_covered=True,source_caps_used_as_defining_field_values=False,
        main_endpoint_envelopes_extrapolated_to_entrance=False,phase_samples_used_as_proof=False)


@source_precision
def run():
    records,hashes,family=current_sources()
    proof=entrance_cone_identities(records);hashes.update(proof['input_hashes'])
    bounds=whole_entrance_bounds(records);hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=family[0],implicit_source_sha256=family[1],
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_source_cone_proof=proof,bounds=bounds,input_hashes=hashes,
        regional_two_vector_cone_distinct_from_completed_tensor_cone=True,
        **{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Continuous entrance full-shear cone generated with exact pre-Tw angular history; global/recursion pending',flush=True)
    return result


class CertifiedPulseEntrancePhysical:
    """Attach the current regional cone to the unchanged full physical source."""
    def __init__(self):
        from lei_ren_part1_paper_compliant_pulse_entrance_physical_C2 import CompliantPulseEntrancePhysicalC2
        self.receipt=json.loads(Path(__file__).with_name(PREFIX+'pulse_entrance_cone_check.json').read_bytes())
        if not self.receipt['all_passed'] or not self.receipt['pulse_entrance_cone_certified']:
            raise ValueError('Accepted entrance cone required')
        for path,digest in self.receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Certified entrance source changed: '+path)
        self.physical=CompliantPulseEntrancePhysicalC2()
        if (self.physical.family,self.physical.source)!=(
            self.receipt['actual_five_defect_family_sha256'],self.receipt['implicit_source_sha256']):
            raise ValueError('Certified entrance source family differs')

    def entrance(self,*args,**kwargs):
        result=self.physical.entrance(*args,**kwargs)
        result.update(**{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
        return result


if __name__=='__main__':run()

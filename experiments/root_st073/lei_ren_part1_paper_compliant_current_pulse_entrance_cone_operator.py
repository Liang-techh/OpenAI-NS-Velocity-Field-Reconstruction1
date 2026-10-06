"""Whole current entrance cone with its actual pre-pulse angular memory.

Only signed error sectors are bounded absolutely. The source tensor, its
incoming energy and absolute pressure are never replaced by these caps.
"""
import ast
import importlib
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_original_cone import (
    SourceAST,pack,encode,endpoints,source_precision,sha)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_operator import (
    CORRECTIONS,BASELINES,ORIGINAL_D_POWERS,upper,absolute)
from lei_ren_part1_paper_compliant_pulse_main_exit_cone import relative_log_recipes
from lei_ren_part1_paper_compliant_pulse_entrance_cone import entrance_log_envelopes
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_tail_bound
from lei_ren_part1_paper_compliant_current_pulse_entrance_incoming_background_tensor import (
    compiled_current_entrance_incoming_exporter,PREFIX)
from lei_ren_part1_paper_compliant_current_pulse_physical_assembly import (
    CompliantPressureDatum,CompliantOuterParameters,class_assignment)

DOMAIN=(0,'.02')


def current_entrance_parameter_bridge(owner):
    """Typed common Tw/mu definition for the current pre and native inlet.

    The selected future may carry a separate mu enclosure. Its pointer is
    not part of this incoming parameter proof.
    """
    pre=owner.physical.pre;native=owner.pulse.pulse.initial
    datums=(pre.datum,native.datum)
    if any(type(d) is not CompliantPressureDatum or type(d.parameters) is not CompliantOuterParameters
            or d.parameters.Md!='40' or d.parameters.precision!=160 for d in datums):
        raise ValueError('Pinned current parameter constructor and inputs required')
    for key in ('definition','source_sha','datum_sha','input_hashes'):
        if getattr(datums[0],key)!=getattr(datums[1],key):
            raise ValueError('Current pre/native parameter source differs: '+key)
    if pre.params is not pre.datum.parameters or native.params is not native.datum.parameters \
            or owner.pulse.pulse.buffer.params is not native.params:
        raise ValueError('Actual pre/native parameter constructor path differs')
    ctor=class_assignment('pressure_source','CompliantPressureDatum','__init__',
        'self.parameters','CompliantOuterParameters(Md,precision)')
    source=owner.proof['actual_both_production_O2_O3_function_equality']['parameter_source_bindings']
    if not source['passed']:raise ValueError('Actual Tw=-60*log_mu source formula required')
    return dict(typed_common_parameter_constructor=True,pinned_inputs=dict(Md='40',precision=160),
        constructor_source_assignment=ctor,actual_mu_Tw_defining_formula=source,
        same_exact_parameter_definition_and_source_hashes=True,
        pre_native_inlet_parameter_object_paths_identified=True,
        cross_instance_parameter_pointer_equality_not_required=True,
        selected_future_mu_enclosure_pointer_not_used_as_incoming_definition=True,
        input_hashes={**source['input_hashes'],PREFIX+'pressure_source.py':sha(PREFIX+'pressure_source.py')},passed=True)


def entrance_memory_theorem():
    """Replay the original power history and the entrance endpoint logs."""
    asts=SourceAST();checks={}
    def zero(name,value):
        if s.cancel(s.expand(s.expand_power_exp(value)))!=0:
            raise ArithmeticError('Current entrance identity failed: '+name)
        checks[name]=True
    mu,xi,Tw,lp,lu,lrp,finite=s.symbols('mu xi Tw logP logU logRp finite',real=True)
    r=1-mu;c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp,ln=s.log)
    u0,h0=s.symbols('actual_pre_power_u actual_pre_power_h',nonzero=True)
    get=lambda key:dict(Utheta_over_Pstar=u0,Mtheta_over_sqrt2_R_3half_Pstar=h0)[key]
    env=dict(self=SimpleNamespace(params=SimpleNamespace(mu=mu,Tw=Tw)),
        c=c,t=Tw,mu=mu,u1=u0,get=get)
    for target in ('slope','f','decay','d3','u'):
        env[target]=asts.evaluate(asts.expression('outer_buffer','power',target),env)
    env['theta_kernel']=asts.evaluate(asts.expression('outer_buffer','power','theta_kernel',
        wanted='(f-d3)/(1-mu)'),env)
    h=asts.evaluate(asts.expression('outer_buffer','power','h'),env)
    zero('actual_power_memory_factored_before_enclosure',h/env['u']-1/r-(h0/u0-1/r)*s.exp(-r*Tw))
    zero('actual_full_Tw_and_entrance_memory',
        (h0/u0-1/r)*s.exp(-r*Tw)*s.exp(-r*xi/mu)-(h0/u0-1/r)*s.exp(-r*(Tw+xi/mu)))
    constants=asts.expression('axial_high_jets','_incoming_constants','self.constants')
    actual=next(node.value for node in constants.keywords if node.arg=='Tw')
    if ast.dump(actual)!=ast.dump(ast.parse('box(initial.params.Tw)',mode='eval').body):
        raise ValueError('Actual high Tw source changed')
    asts.expression('pre_pulse_mixed_C4','power','t',wanted='self.params.Tw*phase')
    asts.expression('outer_buffer','power','t',wanted='self.params.Tw*phase')
    exact=relative_log_recipes(c,mu,finite,lp,lu,lrp,xi)
    envelopes=asts.replay('pulse_entrance_cone','entrance_log_envelopes',
        dict(relative_log_recipes=relative_log_recipes))(c,mu,finite,lp,lu,lrp,Tw)
    for label,parts in CORRECTIONS.items():
        for name in parts:
            key=label+'_'+name;value=exact[key]-r*Tw if name=='signed_original_memory' else exact[key]
            point=s.Rational(1,50) if name=='same_absolute_pressure_memory' else 0
            zero(key+'_actual_entrance_endpoint_envelope',envelopes[key]-value.subs(xi,point))
    return dict(identities=checks,input_hashes=asts.hashes,
        exact_memory='Xp-1/r=(h_pre/u_pre-1/r)*exp(-r*Tw)',
        historical_receipt_checks_loaded_or_promoted=False,passed=True)


def current_entrance_source_theorem(tailcone):
    tailcone.assert_graph();registry=tailcone.registry;owner=registry.owners['incoming'];owner.assert_graph()
    maincone=tailcone.flattencone.gapcone.maincone
    if not owner.acceptance_loaded or not owner.proof['passed'] or owner.main_tensor is not registry.owners['main']:
        raise ValueError('Same checked current entrance and main tensor required')
    if maincone.registry is not registry or not maincone.acceptance_loaded:
        raise ValueError('Same checked current main cone theorem required')
    production=owner.proof['actual_both_production_O2_O3_function_equality']
    live=production['actual_live_program_and_pressure_read_state']
    if not production['passed'] or not all(production['identities'].values()) or not live['passed']:
        raise ValueError('Current actual source function equality required')
    required=('actual_X_Z0_assignment_bound_to_exact_whole_Z_q_cancellation',
        'actual_live_high_and_pulse_ctor_source_constant_expressions_replayed',
        'actual_current_O3_power_phase1_equals_live_native_buffer_power_for_whole_Z',
        'actual_pre_and_native_pressure_datum_calls_replayed_on_same_defining_callable')
    if not all(production[k] for k in required):raise ValueError('Current U, X and absolute datum function witness missing')
    if not all(live['identities']['actual_pre_native_parameter_'+k] for k in ('Md','mu','Tw','yd','logPstar')):
        raise ValueError('Current shared defining parameter input missing')
    if not production['parameter_source_bindings']['passed'] or not all(production['current_source_graph'].values()):
        raise ValueError('Current parameter formula or native graph changed')
    parameter_bridge=current_entrance_parameter_bridge(owner)
    if not parameter_bridge['passed'] or not parameter_bridge['typed_common_parameter_constructor']:
        raise ValueError('Current typed parameter defining source bridge required')
    flat=owner.main_tensor.gap_tensor.end_tensor.flatten_power.flatten
    if flat.U is not owner.pulse.high.constants['U'] or flat.pulse is not owner.pulse:
        raise ValueError('Actual logU is not the current high constant source')
    exporter,extension=compiled_current_entrance_incoming_exporter()
    if owner.exporter.__code__!=exporter.__code__:
        raise ValueError('Current entrance source exporter differs')
    for name in exporter.__code__.co_names:
        if name in exporter.__globals__ and owner.exporter.__globals__.get(name) is not exporter.__globals__[name]:
            raise ValueError('Actual entrance exporter callable differs: '+name)
    module=importlib.import_module(PREFIX+'current_pulse_entrance_incoming_background_tensor')
    if owner.chart.__func__ is not module.CurrentPulseEntranceIncomingBackgroundTensor.chart:
        raise ValueError('Current entrance chart changed')
    for name in ('forward_entrance_energy_rows','gp_energy'):
        if owner.chart.__func__.__wrapped__.__globals__[name] is not getattr(module,name):
            raise ValueError('Original entrance integral changed')
    source=owner.proof['original_generic_entrance_source_theorem']
    physical=owner.proof['original_generic_entrance_physical_theorem']
    if not all(source['identities'].values()) or not all(physical['identities'].values()):
        raise ValueError('Current entrance energy and physical join proof missing')
    if not owner.proof['actual_incoming_entrance_full_tensor_source_function_identified'] \
            or not owner.proof['actual_entrance_main_same_source_function']:
        raise ValueError('Both entrance joins required')
    asts=SourceAST()
    asts.expression('current_pulse_entrance_incoming_background_tensor','chart','forward',
        wanted="forward_entrance_energy_rows(c,mu,xi,data['ap'],data['incoming_energy'],partial,Bh)")
    asts.method('pulse_entrance_similarity_C4','forward_entrance_energy_rows')
    pre_phase1=owner.physical.pre.power(0,1)
    if pre_phase1['chart']!='O3_power_to_Rp' or endpoints(pre_phase1['coverage_coordinate'])!=(1,1) \
            or endpoints(pre_phase1['Utheta_over_Pstar_axial5_coefficients'][0])[0]<=0:
        raise ValueError('Actual current pre-power phase1 U packet required')
    return dict(current_actual_incoming_function_proof=production,
        current_typed_parameter_defining_source_bridge=parameter_bridge,
        current_actual_pre_power_phase1_U_source_packet=pre_phase1,
        current_forward_backward_energy_and_absolute_pressure_theorem=source,
        current_full_incoming_entrance_and_entrance_main_physical_join_theorem=physical,
        current_generic_signed_kernel_shear_and_direction_algebra=maincone.theorem,
        current_actual_log_scale_source=maincone.bindings,
        actual_U_is_live_high_constant_and_pre_power_phase1_function=True,
        actual_Tw_uses_same_defining_recipe_and_inputs_not_parameter_pointer_equality=True,
        actual_Xpre_is_whole_Z_scalar_function_not_Z0_fit=True,
        all_full_incoming_energy_moments_pressure_radial_velocity_and_cross_terms_retained=True,
        forward_and_backward_energy_are_one_function_not_two_added_terms=True,
        original_exporter_unchanged_except_exact_reviewed_entrance_guard=True,
        both_current_function_joins_consumed_before_enclosure=True,
        input_hashes={**asts.hashes,**extension['input_hashes'],**production['input_hashes'],
            **parameter_bridge['input_hashes']},passed=True)


def validate_whole_entrance_view(tailcone,view):
    c=tailcone.ctx
    if tuple(view[k] for k in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')) \
            !=(tailcone.family,tailcone.source,tailcone.datum_sha):raise ValueError('Foreign current entrance source')
    raw=view['original_complete_view']
    if raw['chart']!='pulse_entrance' or endpoints(raw['Z'])!=endpoints(c.mpf([-1,1])) \
            or endpoints(raw['coordinate'])!=endpoints(c.mpf(list(DOMAIN))):
        raise ValueError('Whole current entrance Z[-1,1], xi[0,.02] required')
    for key in ('actual_full_stress_not_local_difference','all_nonzero_incoming_histories_and_radial_remainder_retained',
        'all_local_incoming_cross_and_square_radial_terms_retained','selected_forward_backward_energy_same_source_not_two_added_terms'):
        if not raw[key]:raise ValueError('Incomplete current entrance source: '+key)
    for key in ('current_incoming_moments','current_incoming_energy','current_signed_absolute_Rv_pressure',
        'current_forward_anchored_full_energy_rows','current_source_three_component_velocity_rows'):
        if key not in raw:raise ValueError('Current entrance history omitted: '+key)
    packet=raw['current_actual_source_stress_packet'];sectors=packet['full_meridional_stress_log_sectors']
    if endpoints(c.mpf(packet['exact_logD']))!=(0,0):raise ValueError('Entrance local D must equal one')
    for label in ('theta','axial'):
        if set(sectors[label])!=set(CORRECTIONS[label])|set(BASELINES[label]):
            raise ValueError('Current entrance signed sector omitted or duplicated')
        for name,part in sectors[label].items():
            rp,bp,hp,which=CORRECTIONS[label].get(name,(0,0,0,'one'))
            if tuple(part['mode'])!=(rp+.5,bp+1,ORIGINAL_D_POWERS[label][name],hp) or part['extra_source']!=which:
                raise ValueError('Current entrance signed source mode changed: '+name)
            row=raw['physical_cylindrical_stress_mixed3'][label][name]['r0_z0']
            if encode(pack(row['signed_coefficient']))!=encode(pack(part['full_stress_mixed3_coefficient_enclosures']['s0_Z0'])) \
                    or row['actual_source_log_parts']!=part['exact_source_log_parts']:
                raise ValueError('Actual physical signed coefficient/factor changed: '+name)
    return raw


@source_precision
def whole_current_entrance_bounds(tailcone,view,preTw):
    raw=validate_whole_entrance_view(tailcone,view);c=tailcone.ctx;owner=tailcone.registry.owners['incoming']
    mu=c.mpf(owner.pulse.mu);delta=c.mpf(owner.pulse.delta)
    r=1-mu;lam=c.mpf('.5')-mu;qm=mu-delta/2;k=(1-delta)/2
    def positive(name,value):
        lo,hi=endpoints(value)
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Current entrance positive bound failed: '+name)
    parameters=dict(mu=mu,delta=delta,quarter_minus_mu=c.mpf('.25')-mu,
        rate=r,one_minus_delta=1-delta,lambda1=lam,qmin=qm)
    for name,value in parameters.items():positive(name,value)
    ap=raw['current_selected_ap'];alo,ahi=endpoints(ap[0]);apmax=c.mpf(ahi);apZ=absolute(c,ap[1])
    Kmax=upper(c,r*mu*k/(lam**2*qm));positive('selected_ap_lower',c.mpf(alo))
    positive('K_below_2_01',c.mpf('2.01')-Kmax);positive('ap_below_1_2',c.mpf('1.2')-apmax)
    sig1=upper(c,sigma_tail_bound(c,1,'.5'));g2=50*sig1
    Bmax=c.mpf('.01')*apmax;Dmax=apmax;W0max=2*Bmax+Kmax*Dmax
    gmin=c.mpf(endpoints(qm/r)[0]);epsilon=c.exp(-1000)
    lp=c.mpf(owner.physical.logP);lu=c.ln(c.mpf(owner.main_tensor.gap_tensor.end_tensor.flatten_power.flatten.U))
    lrp=c.mpf(owner.physical.logRp);Tw=c.mpf(owner.physical.pre.params.Tw)
    saddle=owner.pulse.pulse.rows
    finite=-3*c.mpf(saddle['saddle_L'])+2*c.ln(c.mpf(saddle['saddle_u0']))-c.ln(6)/2-2*c.ln(mu)
    if preTw['chart']!='O3_power_to_Rp' or endpoints(preTw['coverage_coordinate'])!=(0,0) \
            or endpoints(preTw['Z'])!=(0,0) or preTw['pulse_inlet']:
        raise ValueError('Current actual pre-Tw power phase0 scalar source required')
    u0=c.mpf(endpoints(preTw['Utheta_over_Pstar_axial5_coefficients'][0]))
    h0=c.mpf(endpoints(preTw['actual_normalized_primitive_y_derivative_axial5']['h'][0][0]))
    positive('actual_pre_Tw_U',u0)
    pre_memory=upper(c,absolute(c,h0/u0)+1/r);memory_coefficient=upper(c,2*pre_memory/(1-delta))
    envelopes=entrance_log_envelopes(c,mu,finite,lp,lu,lrp,Tw)
    bounds={};margins={};monotonicity={}
    for label,parts in CORRECTIONS.items():
        for name,(rp,bp,hp,which) in parts.items():
            coeff=raw['current_actual_source_stress_packet']['full_meridional_stress_log_sectors'][label][name]['full_stress_mixed3_coefficient_enclosures']['s0_Z0']
            bound=memory_coefficient if name=='signed_original_memory' else absolute(c,coeff);key=label+'_'+name
            rate=-rp+bp*(c.mpf('.5')+mu)+hp*r+dict(one=0,incoming1=c.mpf('.5')-mu,
                incoming2=c.mpf('.5')-2*mu,Q=-(1+2*mu),end_square=0)[which]
            monotonicity[key]=-rate if which=='Q' else rate;positive('source_monotonicity_'+key,monotonicity[key])
            target=c.ln(gmin)-1000-c.ln(2*len(parts));value=envelopes[key]+c.ln(bound) if endpoints(bound)[1]>0 else None
            margin=target-value if value is not None else None
            if margin is not None:positive(key,margin);margins[key]=margin
            bounds[key]=dict(actual_signed_whole_entrance_coefficient=coeff,coefficient_absolute_upper=bound,
                grouped_relative_log_envelope=envelopes[key],actual_relative_log_absolute_upper=value,
                certified_log_cap=target,positive_log_gap=margin,structural_zero=endpoints(bound)==(0,0),
                pre_Tw_memory_factorization_used=name=='signed_original_memory',
                absolute_bound_only_for_signed_error_not_defining_stress=True)
    kernel_error=mu**2*g2/lam**3
    local_errors=dict(radial_rate_error=mu/lam*Bmax,two_IBP_kernel_error=r*(1-delta)*apmax*kernel_error/qm,
        axial_amplitude_derivative_error=r*k*apZ*(c.mpf('.01')/lam)/(2*c.sqrt(qm*(1-delta)/2)))
    w_error=upper(c,(sum(local_errors.values(),c.mpf(0))+epsilon*(1+W0max))/(1-epsilon))
    b_error=upper(c,2*mu*apmax*(1+c.mpf('.01')))
    shape=dict(selected_ap_lower=c.mpf(alo),selected_ap_below_1_2=c.mpf('1.2')-apmax,
        Kmax_below_2_01=c.mpf('2.01')-Kmax,w_error_below_one_millionth=c.mpf('1e-6')-w_error,
        b_error_below_one_millionth=c.mpf('1e-6')-b_error)
    for name,value in shape.items():positive(name,value)
    bw_error=Bmax*w_error+W0max*b_error+b_error*w_error;bw_upper=Bmax*W0max+bw_error
    second_upper=2*bw_upper+(Bmax+b_error)**2/2+2*mu*(W0max+w_error)**2
    algebraic=dict(theta_normalized_lower=gmin*(1-epsilon)/2,a_minus_bw_lower=2+2*mu-bw_upper,
        directional_bracket_lower=2-second_upper,kappa_minus2_lower=2*mu,
        bw_below_paper_0_8=c.mpf('.8')-bw_upper,second_below_paper_1_8=c.mpf('1.8')-second_upper,
        full_directional_margin_normalized_lower=4*(2-second_upper))
    for name,value in algebraic.items():positive(name,value)
    return dict(positive_parameter_margins=parameters,positive_log_margins=margins,
        positive_monotonicity_margins=monotonicity,positive_shape_and_error_margins=shape,positive_algebraic_margins=algebraic,
        normalized_stress_sector_bounds=bounds,mu=mu,delta=delta,current_selected_ap_C5_enclosure=ap,
        selected_ap_lower=c.mpf(alo),selected_ap_upper=apmax,selected_ap_Z_absolute_upper=apZ,K_upper=Kmax,
        gp_bounds=dict(value_lower=0,value_upper='.01',first_derivative_lower=0,first_derivative_upper=1,second_derivative_absolute_upper=g2),
        exact_kernel_remainder_absolute_upper=kernel_error,local_w_error_terms=local_errors,
        actual_w_error_absolute_upper=w_error,actual_b_error_absolute_upper=b_error,
        actual_B_absolute_upper=Bmax,actual_leading_w_absolute_upper=W0max,bw_absolute_error_upper=bw_error,
        correlated_bw_upper=bw_upper,correlated_second_cone_expression_upper=second_upper,
        memory_source=dict(pre_Tw_U=u0,pre_Tw_H=h0,pre_Tw_memory_absolute_upper=pre_memory,Tw=Tw,
            factored_memory_coefficient_upper=memory_coefficient,angular_history_not_reset=True,
            actual_pre_phase0_whole_Z_function_ratio_identified_before_Z0_enclosure=True),
        current_source_reduced_log_parameters=dict(logPstar=lp,actual_logU=lu,finite=finite,logRp=lrp),
        strict_current_whole_entrance_two_vector_cone=True,continuous_domain='xi[0,.02], all Z[-1,1]',
        xi0_included_with_full_nonzero_incoming_histories=True,all_fifteen_signed_sectors_retained=True,
        original_energy_cancellation_and_current_absolute_pressure_used=True,
        source_correlations_factored_before_interval_enclosure=True,
        phase_samples_used_as_proof=False,source_caps_used_as_defining_field_values=False)

"""Whole current gap cone with exact D0/D1/D2 and pressure correlations."""
import ast
import copy
import hashlib
import importlib
import mpmath as mp
import sympy as s
from types import SimpleNamespace
from lei_ren_part1_paper_compliant_current_original_cone import (
    HERE,PREFIX,SourceAST,pack,encode,endpoints,source_precision)
from lei_ren_part1_paper_compliant_pulse_gap_cone import FACTORS,relative_log_envelopes,absolute_bound

DOMAINS=dict(pulse_gap=(11,12),pulse_gap_end=(0,1))
LEGACY_STATEMENT_AST_SHA256=(
    'a3d2dd8315a4fb394e6e2960c237ef04f18b3335c64e339917fdbd31e432bfa5',
    '0816b1c83574450f887d3fac9cf8ac8a999567392c223c3a5ce8f50d2b59524e',
    'af0b16398dab5464b6c6746f062c75bcb44a3fc1da5cdb3e311e394211ba8e6c')


def generic_original_gap_cone_theorem():
    module=importlib.import_module(PREFIX+'pulse_gap_cone');asts=SourceAST()
    fn=copy.deepcopy(asts.method('pulse_gap_cone','gap_cone_identities'));fn.decorator_list=[]
    kept=[];removed=[];digests=[]
    for node in fn.body:
        if any(isinstance(v,ast.Name) and v.id in ('records','physical') for v in ast.walk(node)):
            removed.append(ast.unparse(node).splitlines()[0])
            digests.append(hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest())
        else:kept.append(node)
    if tuple(digests)!=LEGACY_STATEMENT_AST_SHA256:raise ValueError('Exact gap historical-check allowlist differs')
    fn.body=kept
    if any(isinstance(v,ast.Name) and v.id in ('records','physical') for v in ast.walk(fn)):
        raise ValueError('Historical gap admission retained in generic theorem')
    env=dict(vars(module))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original arbitrary-source gap cone; current composition separate>','exec'),env)
    result=env[fn.name](None)
    if not all(result['identities'].values()):raise ValueError('Original generic gap cone identity failed')
    result.update(input_hashes={**result['input_hashes'],**asts.hashes},
        omitted_historical_receipt_checks=removed,omitted_statement_AST_sha256=digests,
        historical_receipts_not_loaded_or_promoted=True)
    return result


def current_gap_cone_source_theorem(field,generic):
    field.assert_graph();owner=field.registry.owners['gap'];owner.assert_graph();proof=owner.proof
    if not owner.acceptance_loaded or not proof['passed'] or not proof['current_shared_outer_map_logscale_to_reduced_D0_AST_bridge']:
        raise ValueError('Checked current complete gap/logscale source required')
    if not all(proof[k] for k in ('same_native_complete_future_and_full_squared_beta_source',
        'current_absolute_pressure_source_is_checked_reduced_Rv_datum',
        'selected_source_functions_substituted_before_axial_differentiation',
        'current_full_stress_rows_identified_by_generic_splitter_with_bound_current_histories')):
        raise ValueError('Current selected full gap histories required')
    asts=SourceAST();mu,d,lp,lu,lrp,finite=s.symbols('mu distance logP logU logRp finite',real=True)
    c=SimpleNamespace(mpf=lambda x:x if isinstance(x,s.Expr) else s.Rational(str(x)),ln=s.log)
    native=SimpleNamespace(physical=SimpleNamespace(logRp=lrp,logP=lp),
        end_tensor=SimpleNamespace(flatten_power=SimpleNamespace(flatten=SimpleNamespace(U=s.exp(lu)))))
    env=dict(c=c,self=native,mu=mu,d=d,ss=-d/mu,finite=finite,pp=1+2*mu)
    logB=asts.evaluate(asts.expression('current_pulse_gap_background_tensor','chart','logB',
        wanted="dict(logPstar=self.physical.logP,actual_log_inlet_U=c.ln(c.mpf(self.end_tensor.flatten_power.flatten.U)),inverse_mu=(-13+d)/(2*mu),finite=-13+d)"),env)
    env['logB']=logB
    endB=asts.evaluate(asts.expression('current_pulse_gap_background_tensor','chart','logB',
        wanted="dict(logPstar=logB['logPstar'],actual_log_inlet_U=logB['actual_log_inlet_U'],inverse_mu=-13/(2*mu),finite=-13-(c.mpf('.5')+mu)*ss)"),env)
    if s.expand(sum(endB.values())-sum(logB.values()))!=0:raise ArithmeticError('Current reduced gap-end reference differs')
    logR=asts.evaluate(asts.expression('current_pulse_gap_background_tensor','chart','logR',
        wanted='self.physical.logRp+13/mu+ss'),env)
    D=asts.evaluate(asts.expression('current_pulse_gap_background_tensor','chart','D'),env)
    H=asts.evaluate(asts.expression('current_pulse_gap_background_tensor','chart','H'),env)
    Q=asts.evaluate(asts.expression('current_pulse_gap_background_tensor','chart','Q'),env)
    checks={'current_gap_and_reciprocal_gap_end_reference_logB_identical':True};rates={}
    for label,parts in FACTORS.items():
        for name,(rp,bp,dp,hp,which,pressure) in parts.items():
            key=label+'_'+name
            source=rp*logR+bp*sum(logB.values())+dp*D[int(which[-1])]+hp*H+(Q if pressure else 0)
            old=s.sympify(generic['exact_grouped_relative_log_recipes'][key],locals={
                'mu':mu,'distance':d,'logP':lp,'logU':lu,'logRp':lrp,'finite':finite})
            if s.expand(source-old)!=0:raise ArithmeticError('Current gap grouped recipe differs: '+key)
            checks[key+'_same_actual_current_grouped_log']=True;rate=s.diff(source,d);rates[key]=rate
            endpoint=4*mu if pressure else s.Integer(2)
            if s.simplify(source.subs(d,endpoint)-relative_log_envelopes(c,mu,finite,lp,lu,lrp)[key])!=0:
                raise ArithmeticError('Current exact endpoint envelope differs: '+key)
            checks[key+'_exact_correlated_endpoint']=True
    asts.expression('current_pulse_gap_background_tensor','chart','rows',
        wanted="gap_shapes(c,mu,delta,z,C,c.mpf(self.pulse.Xp),*data['M'],data['future'],data['J0'],data['P0'],d)")
    asts.expression('current_pulse_gap_background_tensor','chart','point',
        wanted="lift_physical_packet(c,packet,delta,rows['velocity'],lt,theta,nu)")
    asts.expression('pulse_gap_similarity_C4','gap_shapes','Bh',wanted='[zero for _ in range(5)]')
    asts.expression('pulse_gap_similarity_C4','gap_shapes','velocity',wanted='pulse_velocity_rows(c,delta,mu,z,C,Bh,m)')
    phase=s.symbols('phase',real=True)
    fn=asts.method('current_pulse_gap_background_tensor','chart')
    candidates=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
        and any(ast.unparse(t)=='d' for t in n.targets)
        and ast.dump(n.value)==ast.dump(ast.parse('4*mu+(1-4*mu)*(1-v)',mode='eval').body)]
    if len(candidates)!=1:raise ValueError('Current reciprocal distance recipe differs')
    reduced=asts.evaluate(candidates[0],dict(mu=mu,v=phase))
    if s.expand(reduced.subs(phase,0)-1)!=0 or s.expand(reduced.subs(phase,1)-4*mu)!=0:
        raise ArithmeticError('Exact current gap endpoint distance differs')
    from lei_ren_part1_paper_compliant_pulse_gap_similarity_C4 import gap_shapes
    from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import pulse_velocity_rows
    live=dict(current_gap_chart_uses_original_shapes=owner.chart.__func__.__wrapped__.__globals__['gap_shapes'] is gap_shapes,
        original_shapes_use_same_velocity_rows=gap_shapes.__globals__['pulse_velocity_rows'] is pulse_velocity_rows)
    if not all(live.values()):raise ValueError('Current gap shear/source callable differs')
    return dict(current_grouped_log_and_endpoint_identities=checks,current_source_log_distance_rates={k:str(v) for k,v in rates.items()},
        current_complete_gap_history_logscale_units_and_trace_theorem=proof,
        inherited_current_main_exit_complete_exit_gap_theorem=field.registry.owners['main'].proof,
        live_original_shape_and_velocity_bindings=live,
        exact_original_gap_domain='4*mu<=distance<=2; gap xi[11,12], gap-end phase[0,1]',
        actual_zero_axial_shear_from_original_Bh_zero=True,nonzero_radial_velocity_and_all_histories_retained=True,
        source_D1_D2_not_replaced_by_local_main_D_one=True,
        current_gap_coordinate_end_and_exit_gap_complete_tensor_attachments_consumed=True,
        input_hashes=asts.hashes,passed=True)


def validate_whole_gap_views(field,views):
    c=field.ctx
    if set(views)!=set(DOMAINS):raise ValueError('Both whole current gap charts required')
    actual={}
    for region,limits in DOMAINS.items():
        view=views[region]
        if (view['actual_five_defect_family_sha256'],view['implicit_source_sha256'],view['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
            raise ValueError('Foreign current gap source')
        raw=view['original_complete_view'];actual[region]=raw
        if raw['chart']!=region or endpoints(raw['Z'])!=endpoints(c.mpf([-1,1])) or endpoints(raw['coordinate'])!=endpoints(c.mpf(limits)):
            raise ValueError('Complete current gap domain required')
        velocity=raw['current_source_three_component_velocity_rows']
        if any(endpoints(row[n])!=(0,0) for row in velocity['axial'] for n in range(6)):
            raise ValueError('Actual gap axial velocity/shear is not structurally zero')
        sectors=raw['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
        for label,parts in FACTORS.items():
            if set(sectors[label])!=set(parts)|({'equilibrium'} if label=='theta' else set()):
                raise ValueError('Current complete signed gap sectors required')
            for name,part in sectors[label].items():
                rp,bp,dp,hp,which,pressure=parts.get(name,(0,0,0,0,'D0',False))
                if tuple(part['mode'])!=(rp+.5,bp+1,dp,hp) or part['selected_D_recipe']!=which \
                        or ('same_absolute_pressure_memory' in part['exact_source_log_parts'])!=pressure:
                    raise ValueError('Current original gap factor changed: '+name)
                row=raw['physical_cylindrical_stress_mixed3'][label][name]['r0_z0']
                if encode(pack(row['signed_coefficient']))!=encode(pack(part['full_stress_mixed3_coefficient_enclosures']['s0_Z0'])) \
                        or row['actual_source_log_parts']!=part['exact_source_log_parts']:
                    raise ValueError('Original signed current gap physical factor differs')
    return actual


@source_precision
def whole_current_gap_bounds(field,views):
    actual=validate_whole_gap_views(field,views);c=field.ctx;owner=field.registry.owners['gap']
    mu=c.mpf(owner.pulse.mu);delta=c.mpf(owner.pulse.delta);r=1-mu;b0=(1-delta)/2;g=(mu-delta/2)/r
    parameters=dict(mu=mu,delta=delta,quarter_minus_mu=c.mpf('.25')-mu,rate=r,one_minus_delta=1-delta,b0=b0,equilibrium_gap=g)
    def positive(name,value):
        lo,hi=endpoints(value)
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Current whole gap positive bound failed: '+name)
    for name,value in parameters.items():positive(name,value)
    g_lower=c.mpf(endpoints(g)[0]);theta_lower=g_lower/4
    eq_upper=c.mpf(endpoints((g+2*b0/r)/(1-delta))[1]);theta_upper=eq_upper+g_lower/8
    epsilon=c.exp(-1000);axial_cap=theta_lower*epsilon
    lp=c.mpf(owner.physical.logP);lu=c.ln(c.mpf(owner.end_tensor.flatten_power.flatten.U));lrp=c.mpf(owner.physical.logRp)
    saddle=owner.pulse.pulse.rows
    finite=-3*c.mpf(saddle['saddle_L'])+2*c.ln(c.mpf(saddle['saddle_u0']))-c.ln(6)/2-2*c.ln(mu)
    envelopes=relative_log_envelopes(c,mu,finite,lp,lu,lrp);bounds={};margins={};monotonicity={}
    for label,parts in FACTORS.items():
        for name,(rp,bp,dp,hp,which,pressure) in parts.items():
            key=label+'_'+name
            coeff=[actual[region]['current_actual_source_stress_packet']['full_meridional_stress_log_sectors'][label][name]
                ['full_stress_mixed3_coefficient_enclosures']['s0_Z0'] for region in DOMAINS]
            bound=c.mpf(max(endpoints(absolute_bound(c,v))[1] for v in coeff))
            target=c.ln(g_lower/64) if label=='theta' else c.ln(theta_lower)-1000-c.ln(4)
            value=envelopes[key]+c.ln(bound) if endpoints(bound)[1]>0 else None
            gap=target-value if value is not None else None
            if gap is not None:positive(key,gap);margins[key]=gap
            bounds[key]=dict(actual_signed_gap_and_gapend_coefficients=coeff,coefficient_absolute_upper=bound,
                relative_mode=(rp,bp,dp,hp),selected_D_recipe=which,pressure_memory_retained=pressure,
                exact_whole_gap_relative_log_envelope=envelopes[key],actual_relative_log_absolute_upper=value,
                certified_log_cap=target,positive_log_gap=gap,structural_zero=endpoints(bound)==(0,0),
                absolute_bound_only_for_error_around_signed_equilibrium=True)
    # Signed exact distance derivatives. Pressure memory decreases with d;
    # its maximum is at 4mu. The other sources increase and peak at d=2.
    rates=dict(theta_signed_original_memory=r/mu,theta_meridional_Mz_history=1/mu,
        theta_mixed_angular_axial_history=1/mu-1,theta_radial_shear=1/mu,
        axial_full_unperturbed_energy_and_pressure=1/(2*mu)+1,
        axial_same_absolute_pressure_memory=-(1/(2*mu)+1),
        axial_linear_axial_moment_history=1/(2*mu)-1,
        axial_selected_backward_energy_loss=1/(2*mu)+1)
    for name,value in rates.items():
        directed=-value if name=='axial_same_absolute_pressure_memory' else value
        positive(name,directed);monotonicity[name]=directed
    algebraic=dict(theta_lower=theta_lower,theta_upper=theta_upper,
        equilibrium_lower_slack=g_lower/2-4*g_lower/64-theta_lower,
        equilibrium_upper_slack=g_lower/8-4*g_lower/64,
        negative_dot_lower=theta_lower,cross_absolute_upper=axial_cap,kappa_minus2=2*mu,
        kappa_minus2_below2=2-2*mu,directional_margin=2*theta_lower**2-2*mu*axial_cap**2,
        normalized_directional_margin=1-mu*epsilon**2)
    for name,value in algebraic.items():positive(name,value)
    return dict(positive_parameter_margins=parameters,positive_log_margins=margins,
        positive_monotonicity_margins=monotonicity,positive_algebraic_margins=algebraic,
        signed_actual_log_distance_derivatives=rates,
        normalized_signed_stress_sector_bounds=bounds,current_reduced_log_parameters=dict(logPstar=lp,logU=lu,finite=finite,logRp=lrp),
        axial_uniform_absolute_upper=axial_cap,exact_axial_shear_ratio=c.mpf(0),exact_a_minus2=2*mu,
        normalized_cone='2*theta^2-2mu*axial^2>0; actual a=2+2mu,b=0',
        continuous_current_gap_and_gapend_domain_covered=True,
        all_nine_signed_sectors_per_region_retained=True,nonzero_radial_histories_and_remainder_retained=True,
        pressure_and_selected_full_beta_weight_sources_retained=True,
        exact_source_correlations_grouped_before_interval_enclosure=True,
        phase_samples_used_as_proof=False,source_caps_used_as_defining_field_values=False)

"""Correlated whole-main/exit cone on the actual admitted source graph.

The kernel is integrated by parts before taking its enclosure. Absolute
bounds apply only to the signed error sectors, not to the defining stress.
Historical admissions are not admissions for the current field.
"""
import ast
import copy
import hashlib
import importlib
import mpmath as mp
import sympy as s
from types import SimpleNamespace
from lei_ren_part1_paper_compliant_current_original_cone import (
    HERE,PREFIX,SourceAST,pack,encode,endpoints,source_precision)
from lei_ren_part1_paper_compliant_pulse_main_exit_cone import (
    CORRECTIONS,relative_log_recipes,relative_log_envelopes,upper,absolute)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_tail_bound

LEGACY_STATEMENT_AST_SHA256=(
    '79734fee33f38c57a6ad881902c5d9716693f8af256fa319b121784ea1386667',
    '5c82bd1b23b3e57fb2c947cc2f62c1d7216b06927bc8bb2f4f296098240dde8c',
    'd0ffceff2b618e9ed90e27cb9344b0947caf58db60a4240ac1bef51e3592460d',
    '9053cfb4dd54abdcfb97b094148b2b4d21b9022c1ea96a2d2a9e3e748976070c')
BASELINES=dict(theta=('equilibrium',),axial=('axial_transport','linear_axial_moment'))
ORIGINAL_D_POWERS=dict(theta=dict(equilibrium=0,signed_original_memory=0,meridional_transport=1,
    radial_shear=0,incoming_Mz_transport=1,incoming_Mtheta_z_transport=1),
    axial=dict(axial_transport=1,linear_axial_moment=1,full_energy_and_pressure=0,
    nonlinear_meridional_transport=2,axial_radial_shear=1,incoming_Mz_linear_axial_moment=1,
    incoming_Mz_nonlinear_meridional_transport=2,same_absolute_pressure_memory=0,selected_end_energy_loss=2))
DOMAINS=dict(pulse_main=('.02','10'),pulse_exit=('10','11'))


def generic_original_main_cone_theorem():
    """Retain the original mathematics, with exactly four old checks removed."""
    module=importlib.import_module(PREFIX+'pulse_main_exit_cone');asts=SourceAST()
    fn=copy.deepcopy(asts.method('pulse_main_exit_cone','main_cone_identities'))
    fn.decorator_list=[];kept=[];removed=[];digests=[]
    for node in fn.body:
        if any(isinstance(v,ast.Name) and v.id in ('records','flat') for v in ast.walk(node)):
            removed.append(ast.unparse(node).splitlines()[0])
            digests.append(hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest())
        else:kept.append(node)
    if tuple(digests)!=LEGACY_STATEMENT_AST_SHA256:
        raise ValueError('Historical check removal differs from exact statement allowlist')
    fn.body=kept
    if any(isinstance(v,ast.Name) and v.id in ('records','flat') for v in ast.walk(fn)):
        raise ValueError('Historical admission remains in generic theorem')
    env=dict(vars(module))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original main cone theorem; current admission supplied separately>','exec'),env)
    result=env[fn.name](None)
    if not all(result['identities'].values()):raise ValueError('Original main cone identity failed')
    result.update(input_hashes={**result['input_hashes'],**asts.hashes},
        omitted_historical_receipt_checks=removed,omitted_statement_AST_sha256=digests,
        historical_receipts_not_loaded_or_promoted=True,
        scope='Generic original signed-sector, kernel, shear and completed-square theorem only.')
    return result


def current_main_cone_source_theorem(field):
    field.assert_graph();owner=field.registry.owners['main'];owner.assert_graph()
    if not owner.acceptance_loaded or not owner.proof['passed']:
        raise ValueError('Checked current full main/exit owner required')
    needed=('current_C5_selected_ap_incoming_controls_future_algorithms_bound',
        'current_partial_forward_and_remaining_energy_defining_integrals_bound',
        'current_reduced_absolute_Rv_pressure_and_native_Pin_retained',
        'current_selected_two_row_inverse_and_energy_equation_instantiates_exit_gap_theorem',
        'all_local_incoming_and_square_meridional_cross_terms_retained')
    if not all(owner.proof[k] for k in needed):raise ValueError('Current complete-history composition required')
    gap=owner.gap_tensor
    if not gap.proof['current_shared_outer_map_logscale_to_reduced_D0_AST_bridge']:
        raise ValueError('Current reduced end factor source theorem required')
    selected=owner.history.selected
    if not (owner.pulse.pulse is selected.amplitude.pulse
            and selected.amplitude.log_end_scale is owner.pulse.pulse.logscale):
        raise ValueError('Current actual shared end logscale object differs')
    asts=SourceAST();mu,xi,lp,lu,lrp,finite=s.symbols('mu xi logP logU logRp finite',real=True)
    c=SimpleNamespace(mpf=lambda x:x if isinstance(x,s.Expr) else s.Rational(str(x)),ln=s.log)
    native=SimpleNamespace(physical=SimpleNamespace(logRp=lrp,logP=lp),
        gap_tensor=SimpleNamespace(end_tensor=SimpleNamespace(flatten_power=SimpleNamespace(
            flatten=SimpleNamespace(U=s.exp(lu))))),pulse=SimpleNamespace(logE=-1/mu+finite))
    env=dict(c=c,self=native,xi=xi,mu=mu)
    logR=asts.evaluate(asts.expression('current_pulse_main_exit_background_tensor','chart','logR',
        wanted='self.physical.logRp+xi/mu'),env)
    logB=asts.evaluate(asts.expression('current_pulse_main_exit_background_tensor','chart','logB'),env)
    logH=asts.evaluate(asts.expression('current_pulse_main_exit_background_tensor','chart','logH',
        wanted='-(1-mu)*xi/mu'),env)
    env['logD0']=asts.evaluate(asts.expression('current_pulse_main_exit_background_tensor','chart','logD0',
        wanted='c.mpf(self.pulse.logE)'),dict(env,c=SimpleNamespace(mpf=lambda x:x)))
    extras=asts.evaluate(asts.expression('current_pulse_main_exit_background_tensor','chart','extras'),env)
    recipes=relative_log_recipes(c,mu,finite,lp,lu,lrp,xi);checks={}
    for label,parts in CORRECTIONS.items():
        for name,(rp,bp,hp,which) in parts.items():
            key=label+'_'+name
            if s.expand(rp*logR+bp*sum(logB.values())+hp*logH+extras[which]-recipes[key])!=0:
                raise ArithmeticError('Current exact grouped log recipe differs: '+key)
            checks[key+'_current_grouped_log_recipe']=True
    asts.expression('current_pulse_main_exit_background_tensor','_data','ap',
        wanted="copy_jet(c,selected['selected_ap_Taylor'])")
    asts.expression('current_pulse_main_exit_background_tensor','chart','kernels',
        wanted="[partial_linear_kernel(c,mu,xi,c.mpf('.5')-i*mu,self.cells,self.window) for i in (1,2)]")
    asts.expression('current_pulse_main_exit_background_tensor','chart','rows',
        wanted="main_exit_shapes(c,mu,delta,z,C,c.mpf(self.pulse.Xp),data['ap'],data['incoming'],[v['enclosure'] for v in kernels],remaining,data['future'],data['J0'],data['P0'],xi)")
    return dict(current_grouped_log_identities=checks,current_actual_selected_history_theorem=owner.proof,
        current_actual_velocity_shear_and_lift_theorem=field.bindings,
        current_shared_end_logscale_AST_bridge=gap.proof,
        same_checked_current_selected_ap_incoming_energy_future_loss_and_absolute_pressure=True,
        original_partial_kernel_integral_not_replaced_by_point_or_truncation=True,
        all_signed_source_rows_retained=True,input_hashes=asts.hashes,passed=True)


def validate_whole_views(field,views):
    c=field.ctx;actual={}
    if set(views)!=set(DOMAINS):raise ValueError('Both whole original main/exit views required')
    for region,limits in DOMAINS.items():
        view=views[region]['original_complete_signed_tensor_view']
        if (view['actual_five_defect_family_sha256'],view['implicit_source_sha256'],view['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
            raise ValueError('Foreign current whole-box source')
        raw=view['original_complete_view']
        if raw['chart']!=region or endpoints(raw['Z'])!=endpoints(c.mpf([-1,1])) \
                or endpoints(raw['coordinate'])!=endpoints(c.mpf(list(limits))):
            raise ValueError('Whole original main/exit box required')
        packet=raw['current_actual_source_stress_packet'];sectors=packet['full_meridional_stress_log_sectors']
        if endpoints(c.mpf(packet['exact_logD']))!=(0,0):
            raise ValueError('Main/exit local D must be exactly one; end loss is a separate factor')
        for label in ('theta','axial'):
            if set(sectors[label])!=set(CORRECTIONS[label])|set(BASELINES[label]):
                raise ValueError('Signed source sector omitted or duplicated')
            for name,part in sectors[label].items():
                expected=CORRECTIONS[label].get(name,(0,0,0,'one'))
                rp,bp,hp,which=expected
                if tuple(part['mode'])!=(rp+.5,bp+1,ORIGINAL_D_POWERS[label][name],hp) or part['extra_source']!=which:
                    raise ValueError('Actual signed source mode changed: '+name)
                row=raw['physical_cylindrical_stress_mixed3'][label][name]['r0_z0']
                if encode(pack(row['signed_coefficient']))!=encode(pack(part['full_stress_mixed3_coefficient_enclosures']['s0_Z0'])) \
                        or row['actual_source_log_parts']!=part['exact_source_log_parts']:
                    raise ValueError('Source signed coefficient or physical factor differs')
        actual[region]=raw
    if encode(pack(actual['pulse_main']['current_selected_ap']))!=encode(pack(actual['pulse_exit']['current_selected_ap'])):
        raise ValueError('Main/exit must share the same actual whole-Z selected amplitude')
    return actual


@source_precision
def whole_current_main_exit_bounds(field,views):
    actual=validate_whole_views(field,views);c=field.ctx;owner=field.registry.owners['main']
    mu=c.mpf(owner.pulse.mu);delta=c.mpf(owner.pulse.delta)
    r=1-mu;lam=c.mpf('.5')-mu;qm=mu-delta/2;k=(1-delta)/2
    def positive(name,value):
        lo,hi=endpoints(value)
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Current whole main/exit positive bound failed: '+name)
    parameters=dict(mu=mu,delta=delta,quarter_minus_mu=c.mpf('.25')-mu,
        rate=r,one_minus_delta=1-delta,lambda1=lam,qmin=qm)
    for name,value in parameters.items():positive(name,value)
    ap=actual['pulse_main']['current_selected_ap'];alo,ahi=endpoints(ap[0])
    apmax=c.mpf(ahi);apZ=absolute(c,ap[1]);positive('selected_ap_lower',c.mpf(alo))
    Kmax=upper(c,r*mu*k/(lam**2*qm));positive('K_below_2_01',c.mpf('2.01')-Kmax)
    positive('ap_below_1_2',c.mpf('1.2')-apmax)
    sig1=upper(c,sigma_tail_bound(c,1,'.5'));sig2=upper(c,sigma_tail_bound(c,2,'.5'))
    g1=1+11*sig1;g2=52*sig1+11*sig2
    Bmax=11*apmax;Dabs=apmax*g1;Dupper=apmax;W0max=2*Bmax+Kmax*Dabs
    gmin=c.mpf(endpoints(qm/r)[0]);epsilon=c.exp(-1000)
    lp=c.mpf(owner.physical.logP)
    lu=c.ln(c.mpf(owner.gap_tensor.end_tensor.flatten_power.flatten.U));lrp=c.mpf(owner.physical.logRp)
    saddle=owner.pulse.pulse.rows
    finite=-3*c.mpf(saddle['saddle_L'])+2*c.ln(c.mpf(saddle['saddle_u0']))-c.ln(6)/2-2*c.ln(mu)
    envelopes=relative_log_envelopes(c,mu,finite,lp,lu,lrp)
    bounds={};margins={};monotonicity={}
    for label,parts in CORRECTIONS.items():
        for name,(rp,bp,hp,which) in parts.items():
            coeff=[actual[region]['current_actual_source_stress_packet']['full_meridional_stress_log_sectors'][label][name]
                ['full_stress_mixed3_coefficient_enclosures']['s0_Z0'] for region in DOMAINS]
            bound=c.mpf(max(endpoints(absolute(c,value))[1] for value in coeff));key=label+'_'+name
            magnitude=-rp+bp*(c.mpf('.5')+mu)+hp*(1-mu)+dict(one=0,incoming1=c.mpf('.5')-mu,
                incoming2=c.mpf('.5')-2*mu,Q=-(1+2*mu),end_square=0)[which]
            monotonicity[key]=-magnitude if which=='Q' else magnitude
            positive('source_monotonicity_'+key,monotonicity[key])
            target=c.ln(gmin)-1000-c.ln(2*len(parts))
            value=envelopes[key]+c.ln(bound) if endpoints(bound)[1]>0 else None
            margin=target-value if value is not None else None
            if margin is not None:positive(key,margin);margins[key]=margin
            bounds[key]=dict(actual_signed_main_and_exit_coefficients=coeff,coefficient_absolute_upper=bound,
                grouped_relative_log_envelope=envelopes[key],actual_relative_log_absolute_upper=value,
                certified_log_cap=target,positive_log_gap=margin,structural_zero=endpoints(bound)==(0,0),
                coefficient_absolute_bound_used_only_for_signed_error_not_defining_stress=True)
    kernel_error=mu**2*g2/lam**3
    local_errors=dict(radial_rate_error=mu/lam*Bmax,
        two_IBP_kernel_error=r*(1-delta)*apmax*kernel_error/qm,
        axial_amplitude_derivative_error=r*k*apZ*(11/lam)/(2*c.sqrt(qm*(1-delta)/2)))
    local_total=sum(local_errors.values(),c.mpf(0))
    w_error=upper(c,(local_total+epsilon*(1+W0max))/(1-epsilon))
    b_error=upper(c,2*mu*apmax*(g1+11))
    shape=dict(selected_ap_lower=c.mpf(alo),selected_ap_below_1_2=c.mpf('1.2')-apmax,
        Kmax_below_2_01=c.mpf('2.01')-Kmax,w_error_below_one_millionth=c.mpf('1e-6')-w_error,
        b_error_below_one_millionth=c.mpf('1e-6')-b_error)
    for name,value in shape.items():positive(name,value)
    A=Kmax*Dupper;bw_error=Bmax*w_error+W0max*b_error+b_error*w_error
    bw_upper=A**2/8+bw_error
    second_upper=2*A**2/7+2*bw_error+Bmax*b_error+b_error**2/2+2*mu*(W0max+w_error)**2
    algebraic=dict(theta_normalized_lower=gmin*(1-epsilon)/2,a_minus_bw_lower=2+2*mu-bw_upper,
        bw_below_paper_0_8=c.mpf('.8')-bw_upper,second_below_paper_1_8=c.mpf('1.8')-second_upper,
        directional_bracket_lower=2-second_upper,kappa_minus2_lower=2*mu,
        full_directional_margin_normalized_lower=4*(2-second_upper))
    for name,value in algebraic.items():positive(name,value)
    return dict(positive_parameter_margins=parameters,positive_log_margins=margins,
        positive_monotonicity_margins=monotonicity,positive_shape_and_error_margins=shape,
        positive_algebraic_margins=algebraic,normalized_stress_sector_bounds=bounds,
        mu=mu,delta=delta,current_selected_ap_C5_enclosure=ap,selected_ap_lower=c.mpf(alo),selected_ap_upper=apmax,
        selected_ap_Z_absolute_upper=apZ,K_upper=Kmax,
        gp_bounds=dict(value_lower=0,value_upper=11,first_derivative_upper=1,
            first_derivative_absolute_upper=g1,second_derivative_absolute_upper=g2),
        exact_kernel_remainder_absolute_upper=kernel_error,local_w_error_terms=local_errors,
        actual_w_error_absolute_upper=w_error,actual_b_error_absolute_upper=b_error,
        actual_B_absolute_upper=Bmax,actual_B_xi_upper=Dupper,actual_B_xi_absolute_upper=Dabs,
        actual_leading_w_absolute_upper=W0max,bw_absolute_error_upper=bw_error,
        correlated_bw_upper=bw_upper,correlated_second_cone_expression_upper=second_upper,
        current_source_reduced_log_parameters=dict(logPstar=lp,actual_logU=lu,finite=finite,logRp=lrp),
        strict_current_whole_main_exit_two_vector_cone=True,
        continuous_domain=dict(pulse_main='xi[.02,10], Z[-1,1]',pulse_exit='xi[10,11], Z[-1,1]'),
        full_order_one_axial_shear_and_falling_side_retained=True,
        all_fifteen_signed_theta_and_axial_sectors_retained=True,
        original_energy_cancellation_and_current_absolute_pressure_used=True,
        exact_source_correlations_grouped_before_interval_enclosure=True,
        phase_samples_used_as_proof=False,source_caps_used_as_defining_field_values=False)

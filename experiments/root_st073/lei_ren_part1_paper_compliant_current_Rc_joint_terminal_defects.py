"""Source-correlated finite-N Rc defects from the accepted complete route.

Reuses actual saved inverse jets. Forms k-a_eff*m before phase/branch
hulls throughout axial, buffer and transition. The slope inlet still has
only component enclosures. Two-tile covers are not a point/control oracle.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_buffer_period_transport as buffer
import lei_ren_part1_paper_compliant_current_generic_moment_repair_operator as repair

native,complete,rc,prior=buffer.native,buffer.complete,buffer.rc,buffer.prior
common,current,supported=buffer.common,buffer.current,buffer.supported
HERE,PREFIX,sha,encode,ep,iv=buffer.HERE,buffer.PREFIX,buffer.sha,buffer.encode,buffer.ep,buffer.iv
N=buffer.N;KEYS=buffer.KEYS;DROW=repair.ROWS[1]
NAME=PREFIX+'current_Rc_joint_terminal_defects.json'
RECEIPT=PREFIX+'current_Rc_joint_terminal_defects_check.json'
GATE='current_actual_fixed_N_source_correlated_Rc_five_terminal_defect_covers_executed'


def restored(rows,coordinates):
    return {key:complete.restore_half_source(value,coordinates) for key,value in rows.items()}


def positive_lower(v):
    if ep(v.coefficient)[0]<=0:raise ValueError('Original strictly positive amplitude/parameter required')
    return v.scale.evaluate()+v.ctx.ln(v.ctx.mpf(ep(v.coefficient)[0]))


def amplitude(coordinates,axial):
    """A=EY*exp(-7)*exp(-5mu/2), with tiny nonzero increment separate."""
    c=coordinates.ctx;t=coordinates.scalar(1)
    parent=axial['original_endpoint_parent_binding']
    mu_interval=iv(c,parent['original_saved_native_inputs']['original_mu']);mu=t.scalar(mu_interval)
    EY=[complete.restore_half_source(v,coordinates) for v in parent['original_Utheta_endpoint_values']['right']]
    epsilon=rc.density.density.factored_expm1(mu*c.mpf('-2.5'))
    base=EY[0]*c.exp(-7);baseZ=EY[1]*c.exp(-7)
    A=base+base*epsilon;AZ=baseZ+baseZ*epsilon
    z=c.mpf(axial['exact_Z_range']);ratio=t.scalar(-2*z/(1+z*z))
    return dict(A=A,A_Z=AZ,mu=mu,mu_interval=mu_interval,log_A_lower=positive_lower(A),log_mu_lower=positive_lower(mu),ratio=ratio,
        epsilon=epsilon,record=dict(original_EY=EY[0].record(),original_EY_Z=EY[1].record(),
        original_Rc_amplitude=A.record(),original_Rc_amplitude_Z=AZ.record(),original_mu=mu.record(),
        original_amplitude_ratio_Z=ratio.record(),original_nonzero_relative_mu_increment=epsilon.record(),
        exact_identity='A=EY*exp(-7-5*mu/2)=Pstar^(-1/2)*exp(-6/5-5*mu/2)/(1+Z^2)',
        buffer_width=11,transition_width=1,transition_J_at_one='1/2',quiet_power_width=2,
        exact_positive_mu_independent_of_Z=True,base_and_nonzero_increment_retained_separately=True,
        original_P0_not_part_of_amplitude_or_redefined=True))


def joint_density(E,EZ,V,VZ,primitive,epsilon,support=True,Ey=None):
    """epsilon=exp(-mu*remaining_sigma_integral)-1; a_eff=E*(1+epsilon).

    C=V*dE+dE*dV-E*epsilon*dV, evaluated before any branch hull.
    In O2 E_y=-E/2 and epsilon_y=0, giving the fixed-phase slow row.
    """
    c=E.ctx;factor=c.mpf(1)/N;x=primitive['A']*factor;z=primitive['B_over_Pstar']*factor
    xZ=primitive['A_Z']*factor;zZ=primitive['B_Z_over_Pstar']*factor
    if support:
        cover=object.__new__(supported.cutoff.NativeCutoffQCover)
        increment=cover.supported_expm1(x,N);exponential=cover.supported_exponential(x,N)
    else:
        increment=rc.density.density.factored_expm1(x)
        exponential=c.exp(rc.density.phase.bounded_value(x))
    dE=E*increment;dEZ=EZ*increment+E*exponential*xZ
    mismatch=-E*epsilon;mismatchZ=-EZ*epsilon
    C=V*dE+dE*z+mismatch*z
    CZ=VZ*dE+V*dEZ+dEZ*z+dE*zZ+mismatchZ*z+mismatch*zZ
    result=dict(C=C,C_Z=CZ)
    if Ey is not None:
        if not V.zero or not VZ.zero:raise ValueError('O2 buffer slow derivative requires original zero V')
        xy=primitive['A_y']*factor;zy=primitive['B_y_over_Pstar']*factor
        dEy=Ey*increment+E*exponential*xy
        result['C_y']=dEy*z+dE*zy-Ey*epsilon*z-E*epsilon*zy
    return result


def query_primitives(query,coordinates,key):
    return restored(query[key]['original_A_B_first_derivative_enclosures'],coordinates)


def joint_axial(coordinates,data,amp):
    c=coordinates.ctx;t=coordinates.scalar(1);zero=t.scalar(0);union=rc.density.local.same_source_union
    locals={};records=[];queries=0
    for row in data['actual_original_endpoint_source_queries']:
        roots=row['original_signed_root_source']['original_signed_roots']['E']
        E=complete.restore_half_source(roots['y0_Z0'],coordinates);EZ=complete.restore_half_source(roots['y0_Z1'],coordinates)
        phys=row['original_typed_native_source']['actual_original_ordinary_physical_rows']
        axial=phys['physical_velocity_pressure_y_derivative_Taylor']['Uz'][0]
        invPs=prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,(0,0,0,-2,0)),1,t.ledger)
        V=complete.restore_half_source(axial[0],coordinates)*invPs;VZ=complete.restore_half_source(axial[1],coordinates)*invPs
        local=row['actual_original_density_and_local_integrals'];qs=[]
        for index,q in enumerate(local['original_conditional_phase_inverse_and_density_queries']):
            primitive=query_primitives(q,coordinates,'actual_original_inverse_first_jet')
            den=joint_density(E,EZ,V,VZ,primitive,amp['epsilon'])
            qs.append(den);queries+=1
        den={k:union([q[k] for q in qs]) for k in ('C','C_Z')}
        mass=complete.restore_half_source(local['original_true_width_kernel_factors']['k']['mass'],coordinates)
        value=den['C']*mass;jet=den['C_Z']*mass
        locals[(row['endpoint'],row['label'])]=(value,jet)
        records.append(dict(endpoint=row['endpoint'],label=row['label'],original_saved_phase_query_count=len(qs),
            actual_joint_density_rows=[native.records(q) for q in qs],joint_phase_union=native.records(den),
            original_rate_three_halves_mass=mass.record(),local_C=value.record(),local_C_Z=jet.record()))
    # Rate 3/2 transports the global bracket k-A*exp(downstream/2)*m.
    inc=zero;incZ=zero;steps=[]
    total=iv(c,data['original_whole_axial_transport_with_genuine_incoming']['original_exact_total_log_radius_width'])
    for step in data['original_whole_axial_transport_with_genuine_incoming']['original_increasing_y_source_steps']:
        width=complete.restore_half_source(step['actual_true_geometry']['positive_true_ordinary_y_width']
            if 'positive_true_ordinary_y_width' in step['actual_true_geometry'] else step['actual_true_geometry']['positive_true_log_radius_width'],coordinates)
        geo=dict(width=width,regular=total if step['endpoint']=='middle' else c.mpf(0),scalar_cover=native.bounded(width),record=step['actual_true_geometry'])
        factor=rc.transfer.true_width_kernel(coordinates,geo,'3/2');value,jet=locals.get((step['endpoint'],step['label']),(zero,zero))
        inc=factor['decay']*inc+value;incZ=factor['decay']*incZ+jet
        steps.append(dict(endpoint=step['endpoint'],label=step['label'],original_decay=factor['decay'].record(),
            own_joint_source_exact_zero=(step['endpoint'],step['label']) not in locals,
            cumulative_joint_increment=inc.record(),cumulative_joint_Z_increment=incZ.record()))
    inlet=data['original_whole_axial_transport_with_genuine_incoming']['actual_genuine_slope_incoming']
    values=restored(inlet['actual_incoming_correction_C0'],coordinates);jets=restored(inlet['actual_incoming_correction_Z'],coordinates)
    left=data['original_endpoint_parent_binding']['original_Utheta_endpoint_values']['left']
    EL=complete.restore_half_source(left[0],coordinates);ELZ=complete.restore_half_source(left[1],coordinates)
    incoming=values['k']-EL*values['m']-EL*amp['epsilon']*values['m']
    incomingZ=jets['k']-ELZ*values['m']-EL*jets['m']-amp['epsilon']*(ELZ*values['m']+EL*jets['m'])
    decay=coordinates.decay(total,'3/2')
    return dict(C=decay*incoming+inc,C_Z=decay*incomingZ+incZ,
        slope=dict(C=decay*incoming,C_Z=decay*incomingZ),axial=dict(C=inc,C_Z=incZ),
        record=dict(source_cells=records,ordered_transport_steps=steps,
            slope_inlet_joint_from_component_covers=True,slope_inlet_function_correlation_not_recovered=True,
            original_genuine_2048_cell_slope_inlet_used=True,exact_whole_axial_attenuation=decay.record(),
            actual_original_source_joint_density_queries=queries,
            global_effective_amplitude_identity='a_eff(y)=E(y)*exp(-5*mu/2) on all O2 axial/buffer',
            exact_original_logPstar_identity='log(Pstar)=exp(Md)+11',
            original_nonzero_mu_mismatch_not_erased=True,branch_hull_after_joint_density=True))


def joint_buffer(coordinates,data,amp):
    c=coordinates.ctx;zero=coordinates.scalar(0);union=rc.density.local.same_source_union
    C=zero;CZ=zero;rows=[];count=0
    for row in data['actual_original_buffer_source_rows']:
        local=row['actual_period_source_and_integrals'];roots=restored(local['original_signed_stress_source']['original_signed_slow_roots']['E'],coordinates)
        bybin={i:[] for i in range(buffer.PHASE_BINS)};means={i:[] for i in range(buffer.PHASE_BINS//2)};query_records=[]
        for index,q in enumerate(local['actual_original_nonlinear_phase_density_queries']):
            primitive=query_primitives(q,coordinates,'original_conditional_first_jets')
            den=joint_density(roots['y0_Z0'],roots['y0_Z1'],zero,zero,primitive,amp['epsilon'],Ey=roots['y1_Z0'])
            bin_index=q['original_phase_bin'];bybin[bin_index].append(den)
            mean=buffer.reflected_density_mean(roots['y0_Z0'],primitive,N)['k'] if bin_index<buffer.PHASE_BINS//2 else None
            if mean is not None:means[bin_index].append(mean)
            query_records.append(dict(original_saved_query_index=index,original_phase_bin=bin_index,
                original_signed_u_branch=q['original_signed_u_branch'],joint_rows=native.records(den),
                frozen_joint_reflection_mean=None if mean is None else mean.record()));count+=1
        assert all(bybin.values()) and all(means.values())
        bins=[{k:union([d[k] for d in bybin[i]]) for k in ('C','C_Z','C_y')} for i in bybin]
        whole=union([d['C'] for d in bins]);whole_y=union([d['C_y'] for d in bins])
        average=sum((union(v)*(c.mpf(2)/buffer.PHASE_BINS) for v in means.values()),zero)
        slow_error=buffer.absolute_cover(whole_y)*(c.mpf(1)/N)
        weight_error=buffer.absolute_cover(whole)*(2*c.expm1(c.mpf('1.5')/N))
        factors=local['original_true_width_Duhamel_factors']['k'];mass=complete.restore_half_source(factors['mass'],coordinates);decay=complete.restore_half_source(factors['decay'],coordinates)
        val=mass*(average+buffer.symmetric_cover(slow_error+weight_error))
        weight=c.exp(c.mpf((-ep(c.mpf('1.5')/N)[1],ep(c.mpf('1.5')/N)[1])))*(c.mpf(1)/buffer.PHASE_BINS)
        jet=mass*sum((d['C_Z']*weight for d in bins),zero)
        C=decay*C+val;CZ=decay*CZ+jet
        rows.append(dict(original_exact_t_endpoints=row['original_exact_t_endpoints'],actual_joint_phase_queries=query_records,
            frozen_joint_reflection_mean=average.record(),full_joint_density_cover=whole.record(),
            fixed_phase_slow_joint_derivative_cover=whole_y.record(),full_slow_joint_error=slow_error.record(),
            period_joint_weight_error=weight_error.record(),local_joint_C=val.record(),local_joint_C_Z=jet.record(),
            original_rate_three_halves_mass=mass.record(),original_decay=decay.record(),
            original_Z_bin_weight=weight,cumulative_C=C.record(),cumulative_C_Z=CZ.record()))
    return dict(C=C,C_Z=CZ,record=dict(original_source_cells=rows,actual_joint_phase_queries=count,
        exact_complete_original_buffer_width=11,original_C0_reflection_mean_and_explicit_errors=True,
        reflection_cancels_linear_mu_mismatch_exactly=True,Z_direct_phase_bin_integration_no_cancellation_claim=True,
        mixed_yZ_derivative_not_assumed=True,nonlinear_joint_before_any_branch_or_phase_hull=True))


def joint_transition(coordinates,data,amp):
    c=coordinates.ctx;zero=coordinates.scalar(0);union=rc.density.local.same_source_union;C=zero;CZ=zero;rows=[];count=0
    for row in data['original_complete_prefix_source_rows']:
        if 'original_typed_native_source' not in row:
            width=complete.restore_half_source(row['actual_true_geometry']['positive_true_log_radius_width'],coordinates)
            geo=dict(width=width,regular=iv(c,row['actual_true_geometry']['regular_true_log_radius_width_cover']),scalar_cover=native.bounded(width),record=row['actual_true_geometry'])
            decay=rc.transfer.true_width_kernel(coordinates,geo,'3/2')['decay'];C=decay*C;CZ=decay*CZ
            rows.append(dict(label=row['label'],flat_source=True,original_decay=decay.record()));continue
        roots=restored(row['original_signed_root_source']['original_signed_roots']['E'],coordinates)
        J=complete.restore_half_source(row['original_typed_native_source']['original_safe_transition_kernels']['original_kernel_values']['J'],coordinates)
        epsilon=rc.density.density.factored_expm1((coordinates.scalar(c.mpf('2.5'))-J)*(-amp['mu']))
        queries=[]
        for q in row['actual_original_conditioned_first_jet_queries']:
            primitive=restored(q['original_A_B_first_derivative_enclosures'],coordinates)
            queries.append(joint_density(roots['y0_Z0'],roots['y0_Z1'],zero,zero,primitive,epsilon,support=False));count+=1
        local=row['actual_original_local_five_C0_Z_density_and_integral'];factors=local['original_true_width_kernel_factors']['k']
        mass=complete.restore_half_source(factors['mass'],coordinates);decay=complete.restore_half_source(factors['decay'],coordinates)
        val=mass*union([q['C'] for q in queries]);jet=mass*union([q['C_Z'] for q in queries])
        C=decay*C+val;CZ=decay*CZ+jet
        rows.append(dict(label=row['label'],original_remaining_sigma_integral=(coordinates.scalar(c.mpf('2.5'))-J).record(),
            original_nonzero_relative_mu_increment=epsilon.record(),actual_joint_phase_queries=[native.records(q) for q in queries],
            original_rate_three_halves_mass=mass.record(),original_decay=decay.record(),local_C=val.record(),local_C_Z=jet.record()))
    # Prefix ends at 2^-128. Remaining transition is flat, followed by 2 power units.
    quiet=coordinates.decay(c.mpf(3)-c.mpf(1)/2**128,'3/2')
    return dict(C=C*quiet,C_Z=CZ*quiet,record=dict(original_joint_source_cells=rows,actual_joint_phase_queries=count,
        exact_prefix_stop='2^-128',exact_remaining_flat_transition_plus_power_width='3-2^-128',
        original_flat_remaining_attenuation=quiet.record(),global_effective_amplitude_identity='a_eff(s)=E(s)*exp(-mu*(5/2-J_sigma(s)))',
        original_Z_independent_geometry_J_mu=True,original_power_own_source_exact_zero=True,
        original_nonzero_mu_mismatch_and_all_quadratic_terms_retained=True))


def normalize(values,jets,C,CZ,amp):
    A=amp['A'];ratio=amp['ratio'];den=A*A;logA=amp['log_A_lower'];out={};outZ={}
    for row,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
        divisor=A if degree==1 else den
        out[row]=values[key].positive_divide(divisor,degree*logA)
        outZ[row]=jets[key].positive_divide(divisor,degree*logA)-out[row]*ratio*degree
    divisor=den*amp['mu'];logden=2*logA+amp['log_mu_lower']
    out[DROW]=C.positive_divide(divisor,logden)
    outZ[DROW]=(CZ-C*ratio*2).positive_divide(divisor,logden)
    return dict(values={k:out[k] for k in repair.ROWS},Z_derivatives={k:outZ[k] for k in repair.ROWS})


def magnitude(v):
    return None if v.zero else v.record()['log_absolute_upper']


def summary(rows):
    return {key:dict(C0_log_absolute_upper=magnitude(v),Z_log_absolute_upper=magnitude(rows['Z_derivatives'][key]),
        signed_C0=v.record(),signed_Z=rows['Z_derivatives'][key].record()) for key,v in rows['values'].items()}


def component_source_parts(coordinates,bd,ad,td):
    """Keep pressure's rate-zero memory in the four source contributions."""
    raw=ad['original_whole_axial_transport_with_genuine_incoming'];op=raw['original_whole_axial_C0_Z_operator']
    incoming=raw['actual_genuine_slope_incoming'];v=restored(incoming['actual_incoming_correction_C0'],coordinates);z=restored(incoming['actual_incoming_correction_Z'],coordinates)
    ca=restored(op['incoming_C0_Z_decay_coefficients'],coordinates)
    buf=bd['actual_updated_Rd_Rc_correction_and_own_history']['original_whole_buffer_operator']
    cb=restored(buf['incoming_C0_Z_decay_coefficients'],coordinates)
    tr=td['original_Rd_to_Rc_operator'];ct=restored(tr['incoming_C0_Z_decay_coefficients'],coordinates)
    pair=lambda vv,zz,scale:dict(values={k:vv[k]*scale[k] for k in KEYS},Z_derivatives={k:zz[k]*scale[k] for k in KEYS})
    return dict(slope_inlet_component_cover=pair(v,z,{k:ca[k]*cb[k]*ct[k] for k in KEYS}),
        original_axial_joint_source=pair(restored(op['cumulative_signed_increment_C0_enclosures'],coordinates),restored(op['cumulative_signed_increment_Z_enclosures'],coordinates),{k:cb[k]*ct[k] for k in KEYS}),
        original_buffer_joint_source=pair(restored(buf['cumulative_signed_increment_C0_enclosures'],coordinates),restored(buf['cumulative_signed_increment_Z_enclosures'],coordinates),ct),
        original_transition_joint_source=pair(restored(tr['cumulative_signed_increment_C0_enclosures'],coordinates),restored(tr['cumulative_signed_increment_Z_enclosures'],coordinates),{k:coordinates.scalar(1) for k in KEYS}))


def linear_control_enclosures(coordinates,targets,amp):
    """First linear response only; do not admit the nonlinear repaired field."""
    c=coordinates.ctx;mu=amp['mu_interval']
    weights=repair.fresh_weights(c,mu,cells=256);matrix=repair.fresh_linear_inverse(c,mu,weights)
    zero=coordinates.scalar(0);h={};hZ={}
    for i,name in enumerate(repair.CONTROLS):
        h[name]=-sum((targets['values'][key]*matrix['inverse_enclosure'][i][j]*N for j,key in enumerate(repair.ROWS)),zero)
        hZ[name]=-sum((targets['Z_derivatives'][key]*matrix['inverse_enclosure'][i][j]*N for j,key in enumerate(repair.ROWS)),zero)
    return dict(original_exact_integral_matrix_weights=weights,original_exact_inverse_enclosure=matrix,
        first_linear_control_response_C0=native.records(h),first_linear_control_response_Z=native.records(hZ),
        original_control_equation='h_linear=-B_exact(mu)^(-1)*(N*r); nonlinear Q(mu,h)/N still omitted',
        original_five_distinct_control_directions_retained=True,
        exact_integral_weights_enclosed_not_selected=True,
        original_mu_not_replaced_by_zero=True,actual_nonlinear_control_functions_or_repair_admitted=False)


def execute_tile(coordinates,bd,ad,td):
    # Every branch must consume the same original selected parameter and
    # Z-independent geometry, not merely a numerically similar amplitude.
    assert bd['source_family']==ad['source_family']==td['source_family']==coordinates.family
    assert bd['candidate_N']==ad['candidate_N']==td['candidate_N']==N
    assert bd['exact_Z_range']==ad['exact_Z_range']==td['exact_Z_range']
    original_mu=iv(coordinates.ctx,ad['original_endpoint_parent_binding']['original_saved_native_inputs']['original_mu'])
    for row in td['original_complete_prefix_source_rows']:
        if 'original_typed_native_source' not in row:continue
        src=row['original_typed_native_source'];parameter=iv(coordinates.ctx,src['lossless_actual_native_parent_inputs']['original_mu'])
        assert parameter._mpi_==original_mu._mpi_, 'Canonical original mu mismatch'
        assert src['canonical_original_amplitude_source']['exact_source_logPstar_equals_exp_Md_plus_11']
        assert row['actual_true_geometry']['width_and_endpoints_independent_of_Z']
        assert src['lossless_actual_native_parent_inputs']['exact_Z_range']==ad['exact_Z_range']
    assert all(row['actual_period_source_and_integrals']['original_source']['exact_zero_Uz_and_all_its_slow_jets'] for row in bd['actual_original_buffer_source_rows'])
    amp=amplitude(coordinates,ad);ax=joint_axial(coordinates,ad,amp);buf=joint_buffer(coordinates,bd,amp);tr=joint_transition(coordinates,td,amp)
    decayBuffer=coordinates.decay(11,'3/2');decayTransition=coordinates.decay(3,'3/2')
    sources=dict(slope_inlet_component_cover={k:v*decayBuffer*decayTransition for k,v in ax['slope'].items()},
        original_axial_joint_source={k:v*decayBuffer*decayTransition for k,v in ax['axial'].items()},
        original_buffer_joint_source=dict(C=buf['C']*decayTransition,C_Z=buf['C_Z']*decayTransition),
        original_transition_joint_source=dict(C=tr['C'],C_Z=tr['C_Z']))
    zero=coordinates.scalar(0);C=sum((v['C'] for v in sources.values()),zero);CZ=sum((v['C_Z'] for v in sources.values()),zero)
    terminal=bd['actual_updated_Rd_Rc_correction_and_own_history']
    values=restored(terminal['actual_updated_Rc_correction_C0'],coordinates);jets=restored(terminal['actual_updated_Rc_correction_Z'],coordinates)
    targets=normalize(values,jets,C,CZ,amp)
    baselineC=values['k']-amp['A']*values['m'];baselineZ=jets['k']-amp['A_Z']*values['m']-amp['A']*jets['m']
    baseline=normalize(values,jets,baselineC,baselineZ,amp)
    parts=component_source_parts(coordinates,bd,ad,td)
    contributions={key:dict(joint_C=pair['C'].record(),joint_C_Z=pair['C_Z'].record(),
        five_target_diagnostics=summary(normalize(parts[key]['values'],parts[key]['Z_derivatives'],pair['C'],pair['C_Z'],amp))) for key,pair in sources.items()}
    return dict(source_family=coordinates.family,candidate_N=N,exact_Z_range=ad['exact_Z_range'],
        common_directed_coordinate_theorem=coordinates.record(),original_Rc_amplitude_and_parameter=amp['record'],
        actual_original_axial_joint_transport=ax['record'],actual_original_buffer_joint_transport=buf['record'],
        actual_original_transition_joint_transport=tr['record'],source_contributions_to_divided_joint_defect=contributions,
        actual_joint_terminal_C=C.record(),actual_joint_terminal_C_Z=CZ.record(),
        actual_signed_five_terminal_defect_C0=native.records(targets['values']),actual_signed_five_terminal_defect_Z=native.records(targets['Z_derivatives']),
        actual_N_scaled_targets_C0=native.records({key:value*N for key,value in targets['values'].items()}),
        actual_N_scaled_targets_Z=native.records({key:value*N for key,value in targets['Z_derivatives'].items()}),
        actual_five_target_diagnostics=summary(targets),independent_component_hull_baseline_diagnostics=summary(baseline),
        original_five_control_linear_response_enclosures=linear_control_enclosures(coordinates,targets,amp),
        original_five_rows=list(repair.ROWS),original_five_distinct_controls=list(repair.CONTROLS),
        exact_control_equation='B_exact(mu)*h + N*r(N,Z) + Q_exact(mu,h)/N = 0',
        continuous_Z_function_enclosures_not_selected_field_values=True,
        actual_original_phase_and_inverse_outputs_reused_without_rerunning_solvers=True,
        joint_density_formed_before_axial_buffer_transition_branch_hulls=True,
        source_correlated_slope_inlet_and_full_function_oracle_still_missing=True,
        same_canonical_original_mu_and_Z_independent_route_explicitly_checked=True,
        N_inverse_coefficient_split_not_assumed=True,original_P0_and_P0_Z_unchanged=True,
        actual_five_controls_or_functional_terminal_closure_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False))


def load_archive(archive):
    compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
    assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
    assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
    return json.loads(raw)


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();accepted=json.loads((HERE/buffer.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,buffer,accepted['source_family']);common.attach_receipt(hashes,repair,accepted['source_family'])
    c=MPIntervalContext();c.dps=240;archives=[]
    for archive in accepted['actual_original_buffer_period_transport_archives']:
        bd=load_archive(archive);ad=load_archive(bd['accepted_original_axial_archive']);td=load_archive(bd['accepted_original_transition_archive'])
        assert bd['exact_Z_range']==ad['exact_Z_range']==td['exact_Z_range']
        assert bd['candidate_N']==ad['candidate_N']==td['candidate_N']==N
        assert bd['source_family']==ad['source_family']==td['source_family']==accepted['source_family']
        coordinates=native.HalfPstarCoordinates(c,iv(c,bd['common_directed_coordinate_theorem']['common_log_bases'][1]),accepted['source_family'])
        result=execute_tile(coordinates,bd,ad,td);result['accepted_original_buffer_archive']=archive
        encoded=json.dumps(encode(result),indent=2).encode()+b'\n';compressed=gzip.compress(encoded,compresslevel=9,mtime=0)
        tag='positive' if Fraction(bd['exact_Z_range'][0])>0 else 'negative';name=PREFIX+'current_Rc_joint_terminal_defects_'+tag+'.json.gz'
        (HERE/name).write_bytes(compressed);hashes[name]=sha(name)
        archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(encoded),lossless_original_json_sha256=hashlib.sha256(encoded).hexdigest(),exact_Z_range=bd['exact_Z_range']))
        print('Actual source-correlated Rc five defects:',bd['exact_Z_range'],'executed',flush=True)
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=accepted['source_family'],candidate_N=N,
        actual_original_Rc_joint_terminal_defect_archives=archives,original_five_rows=list(repair.ROWS),original_controls=list(repair.CONTROLS),
        actual_axial_buffer_transition_joint_density_transport_and_Rc_normalization_executed=True,
        original_positive_mu_retained_no_zero_replacement=True,
        original_accepted_upstream_inverse_or_density_producers_not_rerun=True,
        source_correlated_slope_inlet_function_control_oracle_and_whole_Z_still_open=True,
        actual_five_controls_or_functional_terminal_closure_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Actual finite-N signed M,D,I,S,Cp C0/Z function covers on two strict-sign tiles. Original joint k-a_eff*m evaluated from saved inverse/source tuples before all axial, buffer and transition unions, with buffer nonlinear reflection/slow/weight errors. Genuine slope inlet still component-enclosed. P0 unchanged. No actual control function oracle, functional terminal closure, whole Z/axis/global N/heat/stress/recursion admission.')
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

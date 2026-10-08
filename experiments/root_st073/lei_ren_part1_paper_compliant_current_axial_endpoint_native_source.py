"""Original axial endpoint amplitude/history/stress, phase and five integrals.

Both end layers use the same accepted source family/P0 and original ODEs.
Microscopic exponential and weighted-source corrections remain formal.
"""
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_axial_endpoint_q_source as axial
import lei_ren_part1_paper_compliant_current_native_cutoff_density_oracle as supported

native=axial.native;complete=axial.complete;rc=native.rc;common=native.common;prior=native.prior
current=native.current;slow=native.slow;first=native.first
HERE,PREFIX,sha,encode,ep,iv=axial.HERE,axial.PREFIX,axial.sha,axial.encode,axial.ep,axial.iv
Taylor=native.ScaledTaylor;ZERO=(0,0);ORDERS=slow.ORDERS;N=complete.N;KEYS=tuple(rc.RATES)
NAME=PREFIX+'current_axial_endpoint_native_source.json'
RECEIPT=PREFIX+'current_axial_endpoint_native_source_check.json'
GATE='original_whole_axial_velocity_history_signed_density_C1_operator_and_genuine_incoming_executed'


def positive_prefix_integral(sigma,distance,rate):
    """Original integral of sigma or sigma^2 against exp(-rate*s)."""
    c=sigma.ctx
    # sigma is monotone with distance from the endpoint. Its whole prefix
    # range has lower zero and directed nonzero upper, never a selected value.
    lo,hi=ep(sigma.coefficient)
    cover=prior.ScaledEnclosure(sigma.scale,c.mpf((0,max(hi,0))),sigma.ledger)
    if rate==0:mass=distance
    elif rate==1:mass=-rc.density.density.factored_expm1(-distance)
    else:raise ValueError('Original endpoint prefix weights are rate 0 or 1')
    return cover*mass


class OriginalAxialEndpoint:
    def __init__(self,coordinates,coord,data,cells=64):
        self.coordinates=coordinates;self.coord=coord;self.c=c=coordinates.ctx;self.t=coordinates.scalar(1)
        self.family=coordinates.family;self.Z=data['exact_Z_range'];self.cells=cells;self.data=data
        inputs=data['original_source_rows'][0]['original_typed_native_source']['lossless_actual_native_parent_inputs']
        self.delta=iv(c,inputs['original_delta']);self.inputs=inputs
        self.P0=Taylor(self.t,[iv(c,v) for v in inputs['original_P0_axial_Taylor_coefficients']])
        self.z=Taylor(self.t,[c.mpf(self.Z),1,0,0,0,0]);self.qi=(1+self.z*self.z).reciprocal()
        self.invPs=prior.ScaledEnclosure(prior.FormalScale(self.t.scale.bases,(0,0,0,-2,0)),1,self.t.ledger)
        self.invP2=prior.ScaledEnclosure(prior.FormalScale(self.t.scale.bases,(0,0,0,-4,0)),1,self.t.ledger)
        self.Eleft=self.qi*c.exp(-c.mpf(1)/5)
        self.Eright=self.qi*prior.ScaledEnclosure(prior.FormalScale(self.t.scale.bases,(0,0,0,-1,0)),c.exp(c.mpf(29)/5),self.t.ledger)
        # The exact original slope endpoint quadrature consists of directed
        # rectangles, not point fits. No upstream owner or old producer runs.
        J,masses=native.original.slope_masses(c,c.mpf(1),cells)
        assert ep(J)==ep(c.mpf('.5'))
        h=self.qi*((c.mpf(5)/8+masses[0])*c.exp(-c.mpf('1.5')))
        self.left=dict(m=self.z*4,h=h,k=h*self.z*4,
            e=(self.z*self.z)*(16*self.invP2)-(self.qi*self.qi)*((c.mpf(5)/12+masses[2]/2)*c.exp(-1)),
            p=(self.qi*self.qi)*(c.mpf('2.5')+masses[1]/2))
        raw={k:Taylor(self.t,[iv(c,v) for v in rows[0]]) for k,rows in inputs['original_history_axial_Taylor_coefficients'].items()}
        self.right=dict(m=raw['m']*c.exp(11),h=raw['h']*c.exp(c.mpf('16.5'))-self.Eright*c.expm1(11),
            k=raw['k']*c.exp(c.mpf('16.5')),e=raw['e']*c.exp(11)+self.Eright*self.Eright*c.mpf('5.5'),
            p=raw['p']-self.Eright*self.Eright*((1-c.exp(-11))/2))
        self.parent_record=dict(source_family=self.family,exact_Z_range=self.Z,
            original_left_E_identity='exp(-1/5)/(1+Z**2)',original_right_E_identity='Pstar**(-1/2)*exp(29/5)/(1+Z**2)',
            actual_original_left_slope_J=J,actual_original_left_slope_masses=masses,original_left_slope_rectangle_cells=cells,
            original_left_history_values=native.records(self.left),original_right_history_values=native.records(self.right),
            original_right_Rd_raw_histories=native.records(raw),original_saved_native_inputs=inputs,
            original_separate_P0=native.records(self.P0),right_endpoint_inverted_original_11_unit_zero_V_ODEs=True,
            original_pre_background_not_pulse_correction=True,source_function_integrals_not_selected_values=True,
            original_Utheta_endpoint_values=native.records(dict(left=self.Eleft,right=self.Eright)))
        def factory(_,values):return Taylor(self.t,values)
        factory.variable=lambda _,value,order:Taylor(self.t,[value,1]+[0]*(order-1))
        self.phys,self.physproof=native.compile_original('pre_pulse_mixed_C4','physical_mixed',dict(
            IntervalTaylor=factory,square=lambda v:v*v,derivative=lambda v:v.derivative()))
    def source(self,geometry):
        c=self.c;t=self.t;z=self.z;rcover=geometry['cover'];side=geometry['side']
        sig=complete.closed_sigma_rows(t,rcover)
        whole_sig=complete.closed_sigma_rows(t,c.mpf((0,ep(rcover)[1])))[0]
        dy=complete.restore_half_source(encode(geometry['record']['original_nonzero_y_minus_endpoint']),self.coordinates)
        ledger=[];Ebase=self.Eleft if side=='left' else self.Eright
        E=Ebase*native.exponential_factor(dy*c.mpf('-.5'),ledger)
        y=geometry['y'];M=self.coord.M
        # Original turnoff D_y^k sigma(1-log(y)/M) through four. The
        # signed Stirling conversion acts once on phase-sigma derivatives.
        ordinary=[1-sig[0]]+[sig[j]*((-1)**(j+1)) for j in range(1,5)] if side=='left' else sig
        signed=((1,),(0,1),(0,-1,1),(0,2,-3,1),(0,-6,11,-6,1))
        B=[ordinary[0]]+[sum((ordinary[j]*((-1/M)**j)*signed[k][j] for j in range(1,k+1)),t.scalar(0))*(1/y**k) for k in range(1,5)]
        V=[z*(b*4) for b in B]
        integrals={}
        if side=='left':
            distance=dy;dm=native.exponential_factor(-distance,ledger)
            d3=native.exponential_factor(distance*c.mpf('-1.5'),ledger)
            positive_mass=rc.density.density.factored_expm1(distance)
            integrals['sigma_weighted']=positive_prefix_integral(whole_sig,distance,0)*native.exponential_factor(distance,ledger)
            integrals['two_sigma_minus_sigma_squared_weighted']=positive_prefix_integral(whole_sig*2,distance,0)*native.exponential_factor(distance,ledger)
            # Preserve cancellation of original left m=4Z before adding the
            # microscopic sigma correction. The other exact ODE formulas
            # retain every nonzero exponential correction separately.
            hist=dict(m=z*4-z*(integrals['sigma_weighted']*dm*4),
                h=self.left['h']*d3+Ebase*d3*positive_mass,
                k=self.left['k']*d3+Ebase*z*(positive_mass-integrals['sigma_weighted'])*d3*4,
                e=self.left['e']*dm+(z*z)*(16*self.invP2)*(positive_mass-integrals['two_sigma_minus_sigma_squared_weighted'])*dm
                    -Ebase*Ebase*(distance*dm*c.mpf('.5')),
                p=self.left['p']+Ebase*Ebase*(rc.density.density.factored_expm1(-distance)*c.mpf('-.5')))
        else:
            distance=-dy;d=native.exponential_factor(distance,ledger)
            d3=native.exponential_factor(distance*c.mpf('1.5'),ledger)
            integrals['JV']=z*(positive_prefix_integral(whole_sig,distance,1)*4)
            integrals['JUV']=self.Eright*z*(positive_prefix_integral(whole_sig,distance,1)*4)
            integrals['JV2']=(z*z)*(positive_prefix_integral(current.square(whole_sig),distance,1)*16)
            one_minus_decay=-rc.density.density.factored_expm1(-distance)
            hist=dict(m=(self.right['m']-integrals['JV'])*d,
                h=(self.right['h']-self.Eright*one_minus_decay)*d3,
                k=(self.right['k']-integrals['JUV'])*d3,
                e=(self.right['e']-integrals['JV2']*self.invP2+self.Eright*self.Eright*(distance*c.mpf('.5')))*d,
                p=self.right['p']-self.Eright*self.Eright*(rc.density.density.factored_expm1(distance)*c.mpf('.5')))
        class Context:
            def mpf(_,value):return native.ProtocolConstant(prior.FormalScale(t.scale.bases),c.mpf(value),t.ledger)
        logU=[z*0-c.mpf('.5')]+[z*0]*3
        physical=self.phys(Context(),self.Z,self.delta,E,logU,V,hist,self.P0,self.invP2)
        return dict(values=dict(u=E,logU=logU,V=V,history=hist,P0=self.P0,physical=physical,dy=dy),
            record=dict(source_family=self.family,exact_Z_range=self.Z,endpoint=side,
                original_true_geometry=geometry['record'],original_endpoint_parent_binding=self.parent_record,
                copied_original_physical_mixed_program=self.physproof,
                actual_original_Utheta_over_Pstar=native.records(E),actual_original_Uz_ordinary_y0_through_y4=native.records(V),
                original_five_history_source_values=native.records(hist),original_signed_weighted_turnoff_integrals=native.records(integrals),
                weighted_source_integrals_are_directed_prefix_enclosures_not_equal_to_majorants=True,
                left_two_sigma_minus_sigma_squared_uses_safe_nonnegative_two_sigma_majorant=True,
                original_profile_exponential_corrections=ledger,original_separate_P0=native.records(self.P0),
                actual_original_ordinary_physical_rows=native.records(physical),
                original_micro_y_increment_retained=not dy.zero,
                original_V_physical_units_before_common_Pstar_shift=True,
                ordinary_y_derivatives_and_all_five_history_ODEs_unmodified=True,
                actual_incoming_correction_not_replaced_by_pre_background=True))


def original_signed_roots(adapter,got,norm,qsource,geometry,logC):
    c=adapter.c;t=adapter.t;v=got['values'];phys=v['physical']
    u=phys['physical_velocity_pressure_y_derivative_Taylor']['Utheta_over_Pstar']
    V=phys['physical_velocity_pressure_y_derivative_Taylor']['Uz'];hist=phys['actual_normalized_primitive_y_derivative_axial5']
    P=[v['P0']+hist['p'][0]]+hist['p'][1:]
    fn,proof=native.compile_original('current_pre_pulse_stress_operator','raw_pre_stress_rows',dict(axial_derivative=lambda value:value.derivative()))
    parts=fn(c,adapter.delta,adapter.z,u,V,hist,P)
    Ps=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,2,0)),1,t.ledger)
    side=geometry['side'];power=20 if side=='left' else 22;offset=1 if side=='left' else -11
    radius_base=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,power,0),10*c.mpf(logC)+c.ln(110)+offset),1,t.ledger)
    exps=[];radius=radius_base*native.exponential_factor(v['dy'],exps);ns={}
    for label,omit in (('theta','variable_radial_shear'),('axial','axial_radial_shear')):
        raw=[u[j]*0 for j in range(3)]
        for name,part in parts[label].items():
            if name==omit:continue
            factor=1 if label=='theta' else adapter.invPs if part['mode'][1]==0 else Ps
            for j in range(3):raw[j]=raw[j]+part['shape'][j]*factor
        ns[label]=[sum((raw[i]*math.comb(j,i) for i in range(j+1)),raw[0]*0)*radius for j in range(3)]
    roots={name:{(j,k):rows[j][k]*math.factorial(k) for j,k in ORDERS}
        for name,rows in (('E',u),('p1_numerator',ns['theta']),('p2_numerator',ns['axial']))}
    low=roots['E'][ZERO].scale.evaluate()+c.ln(c.mpf(ep(roots['E'][ZERO].coefficient)[0]))
    roots['p1']=current.quotient_jet(roots.pop('p1_numerator'),roots['E'],low)
    roots['p2']=current.quotient_jet(roots.pop('p2_numerator'),roots['E'],low)
    roots.update(qsource['roots']);roots['b']=norm['b'];roots['t0']={o:-b*c.mpf('.5') for o,b in norm['b'].items()}
    return dict(roots=roots,q=qsource['rows'][ZERO],qrows=qsource['rows'],record=dict(
        source_family=adapter.family,copied_original_raw_stress_program=proof,
        original_signed_roots={name:{'y%d_Z%d'%o:v.record() for o,v in rows.items()} for name,rows in roots.items()},
        original_endpoint_radius=radius_base.record(),original_current_radius=radius.record(),original_radius_exponential_corrections=exps,
        original_Rref_Pstar_Cstar_relations_retained=True,original_b_and_t0_not_set_to_zero=True,
        actual_absolute_P0_plus_history_p_used_once=True,full_signed_p2_not_replaced_by_axial_velocity=True))


def actual_endpoint_phase(adapter,geometry,N=N):
    c=adapter.c;t=adapter.t;original=adapter.data['original_source_rows'][0]['actual_original_typed_phase']
    base=original['original_base_source_phase'];dy=complete.restore_half_source(encode(geometry['record']['original_nonzero_y_minus_endpoint']),adapter.coordinates)
    delta=-N*(adapter.coord.logP-1) if geometry['side']=='left' else c.mpf(0)
    boxes=[]
    for box in base['periodic_projection']['boxes']:
        boxes.extend(first.spatial.ordinary_mod_one(c,iv(c,box)+delta+native.bounded(dy*N))['boxes'])
    raw=base['actual_log_R_over_same_r_minus_enclosure'];old=raw['formal_positive_scale']
    assert old['source_exponents']==[0,0,0,0] and not old['radius_power']
    origin=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=iv(c,old['additional_log_interval'])),iv(c,raw['coefficient_interval']),t.ledger)
    regular_delta=-(adapter.coord.logP-1) if geometry['side']=='left' else c.mpf(-11)
    return dict(actual_phase_boxes=boxes,original_Rd_source_phase=base,
        original_endpoint_phase_relation='Rd_phase-N*(logPstar-1)' if geometry['side']=='left' else 'Rd_phase-11*N; 11*N integer',
        original_endpoint_log_radius_delta=regular_delta,original_micro_y_component=dy.record(),original_micro_N_y_component=(dy*N).record(),
        actual_original_log_R_over_r_minus=(origin+t.scalar(regular_delta)+dy).record(),
        actual_phase_is_N_times_original_log_radius=True,phase_independent_of_Z=True,
        exact_integer_cycles_collected_before_modulo=True,nonzero_micro_phase_not_selected_or_erased=not dy.zero)


def actual_density_integrals(adapter,got,roots,geometry,phase,dstar):
    c=adapter.c;t=adapter.t;src=dict(roots=roots['roots'],q=roots['q']);qrows=roots['qrows']
    u,branches,empty=supported.cover.signed_u_branches(src,dstar)
    support=object.__new__(supported.cutoff.NativeCutoffQCover)
    kernel,proof=supported.compile_supported_density(support)
    axial=got['values']['physical']['physical_velocity_pressure_y_derivative_Taylor']['Uz'][0]
    E,EZ=src['roots']['E'][ZERO],src['roots']['E'][(0,1)]
    V,VZ=axial[0]*adapter.invPs,axial[1]*adapter.invPs
    cells=[];den=[]
    for phi in phase['actual_phase_boxes']:
        local=[]
        for branch in branches:
            inverse=supported.cover.branch_first_jets(src,qrows,dstar,phi,branch)
            if inverse['values'] is None:raise ArithmeticError('Original axial inverse requires phase refinement')
            density=kernel(E,EZ,V,VZ,inverse['values'],N);local.append(density)
            cells.append(dict(actual_original_phase_box=phi,signed_u_branch=branch['name'],original_condition=branch['condition'],
                actual_original_inverse_first_jet=inverse['record'],
                actual_original_density_C0=native.records(density['kernels']),actual_original_density_Z=native.records(density['Z_derivatives'])))
        den.extend(local)
    if not den:raise ArithmeticError('Original phase union must be complete and nonempty')
    union=rc.density.local.same_source_union
    kernels={k:union([d['kernels'][k] for d in den]) for k in KEYS};jets={k:union([d['Z_derivatives'][k] for d in den]) for k in KEYS}
    width=geometry['physical_width'];geo=dict(width=width,regular=c.mpf(0),scalar_cover=native.bounded(width),record=geometry['record'])
    factors={k:rc.transfer.true_width_kernel(adapter.coordinates,geo,rate) for k,rate in rc.RATES.items()}
    values={k:kernels[k]*factors[k]['mass'] for k in KEYS};derivatives={k:jets[k]*factors[k]['mass'] for k in KEYS}
    return dict(geometry=geo,values=values,Z_derivatives=derivatives,record=dict(status='enclosed',
        original_unconditioned_signed_u=u.record(),original_signed_u_empty_branches=empty,
        original_conditional_phase_inverse_and_density_queries=cells,
        original_nonlinear_density_precedes_overlapping_branch_hull=True,
        complete_original_phase_union_C0=native.records(kernels),complete_original_phase_union_Z=native.records(jets),
        original_supported_exponential_AST_replacements=proof,
        original_A_over_N_and_A_Z_retained_exponential_bounds_only=True,
        original_true_ordinary_y_width=width.record(),physical_width_is_directed_bound_not_selected_exact_scalar=True,
        original_true_width_kernel_factors={k:dict(branch=v['branch'],mass=v['mass'].record(),decay=v['decay'].record()) for k,v in factors.items()},
        actual_original_local_five_C0_integral_contributions=native.records(values),actual_original_local_five_Z_integral_contributions=native.records(derivatives),
        original_common_Pstar_unit_once_and_physical_Jacobian_once=True,
        full_axial_incoming_operator_buffer_or_functional_targets_not_admitted=True))


def original_axial_endpoint_telescope():
    sy=axial.sy;H,M=sy.symbols('H M',positive=True)
    Q={'left':H+3*sy.log(H),'right':H+3*sy.log(H)+sy.exp(M)-1-2*M}
    rstar=sy.sqrt(2/(H-sy.log(H)-200));points=[]
    for side in ('left','right'):
        q=Q[side];distances=[sy.Integer(0),sy.Rational(3,4)*sy.sqrt(2/q)]
        distances += [sy.sqrt(2/(q+k)) for k in (8,0,-3,-sy.Rational(19,5),-sy.Rational(22,5))]
        distances.append(rstar)
        endpoint=[sy.exp(M*r) if side=='left' else sy.exp(M*(1-r)) for r in distances]
        points.append(endpoint)
    ordered=points[0]+list(reversed(points[1]))
    widths=[b-a for a,b in zip(ordered,ordered[1:])]
    assert sum(widths)==sy.exp(M)-1
    return dict(passed=True,original_ordered_ordinary_y_endpoints=[str(v) for v in ordered],
        exact_original_cell_width_expressions=[str(v) for v in widths],
        exact_original_total_log_radius_width='exp(Md)-1',
        adjacent_source_endpoints_collected_before_independent_directed_evaluation=True,
        original_right_endpoint_distance_order_reversed_for_increasing_y=True)


def genuine_slope_incoming(coordinates,Z,tail_report):
    target=tuple(Fraction(v) for v in Z)
    archive=next(a for a in tail_report['original_six_cell_tail_source_archives']
        if tuple(Fraction(v) for v in a['exact_Z_range'])==target)
    raw=gzip.decompress((HERE/archive['filename']).read_bytes())
    assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
    saved=json.loads(raw);parent=saved['genuine_primary_incoming_O2_record']
    assert parent['ordered_source_cells']==2048 and parent['candidate_N']==N
    original=rc.history.CommonSourceCoordinates(coordinates.ctx,coordinates.logP_squared,coordinates.family)
    restored=complete.tail.inlet_functions(parent,original)
    lift=lambda row:{k:coordinates.rebase(v,coordinates.family) for k,v in row.items()}
    return dict(values=lift(restored['values']),Z_derivatives=lift(restored['Z_derivatives']),
        record=dict(source_manifest=complete.tail.NAME,source_manifest_sha256=sha(complete.tail.NAME),
            lossless_actual_upstream_archive=archive,original_genuine_slope_exit_record=parent,
            original_route_endpoint='O2_slope(y=1)=O2_axial(p=0)',
            actual_incoming_correction_C0=native.records(lift(restored['values'])),
            actual_incoming_correction_Z=native.records(lift(restored['Z_derivatives'])),
            actual_incoming_is_genuine_2048_cell_slope_correction_not_background_or_Rd=True))


def whole_axial_operator(coord,data,adapter,locals_by_side,tail_report):
    coordinates=coord.coordinates;c=coord.c;family=adapter.family
    zero={k:coordinates.scalar(0) for k in KEYS}
    operator=rc.history.C1DuhamelOperator(coordinates);steps=[]
    connector_label,kind,left,right=axial.PLAN[-1]
    connectors={}
    for side in ('left','right'):
        geometry=coord.geometry(side,kind,left,right)
        norm=axial.normalized_excess(coord,data['exact_Z_range'],geometry)
        q=axial.original_axial_q(coord,norm['rows'])
        assert q['record']['branch']=='flat' and all(v.zero for v in q['rows'].values())
        connectors[side]=dict(geometry=dict(width=geometry['physical_width'],regular=c.mpf(0),
            scalar_cover=native.bounded(geometry['physical_width']),record=geometry['record']),values=zero,Z_derivatives=zero)
    # Keep both microscopic endpoint spans in the quiet width; the total
    # coefficient below is collected from exact source endpoints.
    x=coord.rstar*coord.M
    spans=dict(left=rc.density.density.factored_expm1(x),
        right=-rc.density.density.factored_expm1(-x)*c.exp(coord.M))
    total=c.expm1(coord.M);width=coordinates.scalar(total)-spans['left']-spans['right']
    if ep(width.coefficient)[0]<=0:raise ArithmeticError('Positive original whole axial quiet gap required')
    middle=dict(geometry=dict(width=width,regular=total,scalar_cover=native.bounded(width),
        record=dict(chart='O2_axial',original_endpoint_specification=['p=rstar','p=1-rstar'],
            exact_original_width='exp(Md*(1-rstar))-exp(Md*rstar)',
            original_nonzero_endpoint_y_spans=native.records(spans),positive_true_log_radius_width=width.record(),
            original_middle_flat_theorem=axial.middle_flat_proof(coord,data['exact_Z_range']),
            source_independent_of_Z=True)),values=zero,Z_derivatives=zero)
    ordered=[('left',row[0],local) for row,local in zip(axial.PLAN[:-1],locals_by_side['left'])]
    ordered += [('left',connector_label,connectors['left']),('middle','strict_flat_middle',middle),('right',connector_label,connectors['right'])]
    ordered += [('right',row[0],local) for row,local in reversed(list(zip(axial.PLAN[:-1],locals_by_side['right'])))]
    for side,label,local in ordered:
        rc.transfer.append_true_cell(operator,local['geometry'],local['values'],local['Z_derivatives'],family)
        steps.append(dict(endpoint=side,label=label,actual_true_geometry=local['geometry']['record'],
            local_signed_C0_increment=native.records(local['values']),local_signed_Z_increment=native.records(local['Z_derivatives']),
            cumulative_operator=operator.record(),zero_own_source_keeps_incoming_and_pressure_memory=side=='middle' or label==connector_label))
    operator.coefficients={k:coordinates.decay(total,r) for k,r in rc.RATES.items()}
    incoming=genuine_slope_incoming(coordinates,data['exact_Z_range'],tail_report)
    outgoing=operator.apply(incoming['values'],incoming['Z_derivatives'],family)
    original={k:adapter.right[k][0] for k in KEYS};originalZ={k:adapter.right[k][1] for k in KEYS}
    own={k:original[k]+outgoing['values'][k] for k in KEYS};ownZ={k:originalZ[k]+outgoing['Z_derivatives'][k] for k in KEYS}
    return dict(original_increasing_y_source_steps=steps,original_exact_endpoint_telescope=original_axial_endpoint_telescope(),
        original_whole_axial_C0_Z_operator=operator.record(),actual_genuine_slope_incoming=incoming['record'],
        original_exact_total_log_radius_width=total,actual_axial_exit_correction_C0=native.records(outgoing['values']),
        actual_axial_exit_correction_Z=native.records(outgoing['Z_derivatives']),
        actual_axial_exit_background_C0=native.records(original),actual_axial_exit_background_Z=native.records(originalZ),
        actual_axial_exit_own_history_C0=native.records(own),actual_axial_exit_own_history_Z=native.records(ownZ),
        actual_axial_exit_absolute_pressure=(adapter.P0[0]+own['p']).record(),
        actual_axial_exit_absolute_pressure_Z=(adapter.P0[1]+ownZ['p']).record(),
        P0_separate_and_added_once=True,entire_original_axial_phase_zero_to_one_covered=True,
        true_geometry_and_quiet_memory_not_replaced_by_endpoint_caps=True,
        original_full_buffer_Rc_targets_controls_global_N_axis_heat_stress_recursion_admitted=False)


def execute_tile(coord,data,tail_report):
    adapter=OriginalAxialEndpoint(coord.coordinates,coord,data);rows=[];locals_by_side={s:[] for s in ('left','right')}
    for side in ('left','right'):
        for label,kind,left,right in axial.PLAN[:-1]:
            geometry=coord.geometry(side,kind,left,right);norm=axial.normalized_excess(coord,data['exact_Z_range'],geometry)
            q=axial.original_axial_q(coord,norm['rows']);got=adapter.source(geometry)
            roots=original_signed_roots(adapter,got,norm,q,geometry,iv(coord.c,data['original_logC']))
            phase=actual_endpoint_phase(adapter,geometry)
            local=actual_density_integrals(adapter,got,roots,geometry,phase,iv(coord.c,data['original_dstar_log']))
            locals_by_side[side].append(local)
            rows.append(dict(endpoint=side,label=label,original_typed_native_source=got['record'],
                original_signed_root_source=roots['record'],actual_original_phase=phase,actual_original_density_and_local_integrals=local['record']))
            print('Original axial endpoint native source:',data['exact_Z_range'],side,label,'enclosed',flush=True)
    operator=whole_axial_operator(coord,data,adapter,locals_by_side,tail_report)
    return dict(source_family=adapter.family,candidate_N=N,exact_Z_range=data['exact_Z_range'],
        common_directed_coordinate_theorem=adapter.coordinates.record(),original_endpoint_parent_binding=adapter.parent_record,
        actual_original_endpoint_source_queries=rows,
        original_middle_flat_theorem=axial.middle_flat_proof(coord,data['exact_Z_range']),
        original_whole_axial_transport_with_genuine_incoming=operator,
        original_signed_history_pressure_and_inverse_sources_kept=True,
        complete_original_axial_operator_on_two_strict_sign_Z_tiles_executed=True,
        full_buffer_global_N_targets_controls_recursion_admitted=False)


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();accepted=json.loads((HERE/axial.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,axial,accepted['source_family']);common.attach_receipt(hashes,supported,accepted['source_family'])
    tail_report,_=common.attach_receipt(hashes,complete.tail,accepted['source_family'])
    c=MPIntervalContext();c.dps=240;archives=[]
    for archive in accepted['original_axial_q_source_archives']:
        raw=gzip.decompress((HERE/archive['filename']).read_bytes());assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
        qdata=json.loads(raw);parent=qdata['accepted_native_input_archive']
        original_raw=gzip.decompress((HERE/parent['filename']).read_bytes());assert hashlib.sha256(original_raw).hexdigest()==parent['lossless_original_json_sha256']
        data=json.loads(original_raw);coordinates=native.HalfPstarCoordinates(c,iv(c,qdata['saved_original_coordinate_theorem']['common_log_bases'][1]),accepted['source_family'])
        coord=axial.AxialCoordinates(coordinates,iv(c,qdata['original_parameter_sources']['original_eta_log']))
        got=execute_tile(coord,data,tail_report);got.update(accepted_original_q_source_archive=archive,accepted_original_native_parent_archive=parent,
            original_logC=data['original_logC'],original_dstar_log=data['original_dstar_log'],original_parameter_sources=coord.record())
        encoded=json.dumps(encode(got),indent=2).encode()+b'\n';compressed=gzip.compress(encoded,compresslevel=9,mtime=0)
        tag='positive' if Fraction(data['exact_Z_range'][0])>0 else 'negative';name=PREFIX+'current_axial_endpoint_native_source_'+tag+'.json.gz'
        (HERE/name).write_bytes(compressed);hashes[name]=sha(name)
        archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(encoded),lossless_original_json_sha256=hashlib.sha256(encoded).hexdigest(),exact_Z_range=data['exact_Z_range']))
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=accepted['source_family'],candidate_N=N,
        actual_original_axial_endpoint_native_source_archives=archives,actual_endpoint_queries=24,actual_local_C0_Z_integral_rows=240,
        same_original_Utheta_Uz_all_five_histories_P0_and_signed_stress_sources_executed=True,
        actual_radius_phase_signed_u_inverse_and_five_local_density_integrals_executed=True,
        accepted_native_owner_or_ancestor_producers_not_rerun=True,
        complete_original_axial_C0_Z_operator_with_genuine_slope_incoming_executed=True,
        full_buffer_global_N_targets_controls_recursion_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Complete original O2 axial phase [0,1] on two strict-sign Z tiles: original endpoint velocity/history/P0/stress, real N1024 phase and signed density integrals, both flat connectors and strict flat middle compose the five C0/Z affine operator with genuine accepted 2048-cell slope incoming. Actual axial-exit correction/own-history/pressure covers retained. Complete buffer/Rc update, full Z/axis, targets/controls/global N/heat/stress/recursion remain open.')
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

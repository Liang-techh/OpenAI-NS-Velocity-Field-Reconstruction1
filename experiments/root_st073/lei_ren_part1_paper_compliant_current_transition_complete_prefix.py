"""Closed original microscopic prefix and contiguous five C0/Z transport.

Uniform flat-endpoint inequalities replace singular interval divisions, not
the defining source. Actual native Rd parents and phase origins are reused
losslessly. The original transition correction has explicit genuine input.
"""
import ast
import copy
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import sympy as sy
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_transition_factored_native_source as native

previous=native.previous;rc=native.rc;common=native.common;prior=native.prior
first=native.first;current=native.current;slow=native.slow
support=previous.previous;tail=support.tail
HERE,PREFIX,sha,encode,ep,iv=native.HERE,native.PREFIX,native.sha,native.encode,native.ep,native.iv
N=native.N;ZERO=(0,0);KEYS=tuple(rc.RATES)
NAME=PREFIX+'current_transition_complete_prefix.json'
RECEIPT=PREFIX+'current_transition_complete_prefix_check.json'
GATE='original_closed_microscopic_prefix_and_contiguous_transition_C0_Z_operator_executed'


def poly_add(*polynomials):
    result={}
    for polynomial in polynomials:
        for power,value in polynomial.items():result[power]=result.get(power,0)+value
    return result


def poly_mul(a,b):
    result={}
    for i,x in a.items():
        for j,y in b.items():result[i+j]=result.get(i+j,0)+x*y
    return result


def weighted_flat_endpoint_rows(template,upper,anchor=None):
    """Original sigma/normalized-sigma derivatives on the whole [0,S]."""
    c=template.ctx;S=c.mpf(ep(c.mpf(upper))[1])
    if ep(S)[0]<=0 or ep(S)[1]>=ep(c.mpf('.1'))[0]:raise ValueError('Small positive closed original endpoint collar')
    # d/ds log(exp(L(s))*s^-m) = L1-m/s >= (2-m*s^2)/s^3.
    # Degree <=12, hence every weighted source has its supremum at S.
    if ep(2-12*S*S)[0]<=0:raise ArithmeticError('Original weighted-flat-source monotonicity')
    F=1/(1-S)**2;L=F-1/S**2
    if anchor is None:anchor=L
    odds=[{0:2/(1-S)**3,3:c.mpf(2)},
          {0:6/(1-S)**4,4:c.mpf(6)},
          {0:24/(1-S)**5,5:c.mpf(24)},
          {0:120/(1-S)**6,6:c.mpf(120)}]
    A,B,C,D=odds;A2=poly_mul(A,A);A3=poly_mul(A2,A);A4=poly_mul(A2,A2)
    polynomials=[{0:c.mpf(1)},A,poly_add(B,A2),
        poly_add(C,{k:3*v for k,v in poly_mul(A,B).items()},A3),
        poly_add(D,{k:4*v for k,v in poly_mul(A,C).items()},
            {k:3*v for k,v in poly_mul(B,B).items()},
            {k:6*v for k,v in poly_mul(A2,B).items()}, {k:2*v for k,v in A4.items()})]
    values=[]
    for order,polynomial in enumerate(polynomials):
        terms=[]
        for power,coefficient in polynomial.items():
            terms.append(prior.ScaledEnclosure(prior.FormalScale(template.scale.bases,
                offset=anchor-power*c.ln(S)),coefficient,template.ledger))
        magnitude=sum(terms,template.scalar(0))
        hi=max(abs(x) for x in ep(magnitude.coefficient))
        values.append(prior.ScaledEnclosure(magnitude.scale,c.mpf((0 if order<=1 else -hi,hi)),template.ledger))
    return dict(values=values,record=dict(original_closed_source_interval=[c.mpf(0),S],
        original_logistic_odds_upper=L,original_weighted_source_anchor=anchor,
        exact_endpoint_all_sigma_ordinary_derivatives_zero=True,
        weighted_source_monotonicity='L1-m/s >= (2-m*s^2)/s^3>0 for 0<s<=S, 0<=m<=12',
        positive_monotonicity_margin=2-12*S*S,
        original_logistic_derivative_absolute_polynomials=polynomials,
        ordinary_y0_through_y4_enclosures=[v.record() for v in values],
        logistic_identities_bound_only_not_source_redefinition=True,
        no_interval_division_by_zero_or_epsilon_endpoint=True,
        entire_source_interval_including_s_zero_enclosed=True))


def closed_sigma_rows(template,coordinate):
    c=template.ctx;lo,hi=ep(c.mpf(coordinate))
    if lo==hi==0:return [template.scalar(0)]*5
    if lo<0:raise ValueError('Original nonnegative radial coordinate')
    if lo>0:return native.original_sigma_rows(template,coordinate)
    return weighted_flat_endpoint_rows(template,c.mpf(hi))['values']


def source_programs():
    kernel,kproof=native.compile_original('current_transition_factored_native_source','safe_transition_kernels',
        dict(original_sigma_rows=closed_sigma_rows))
    tree=ast.parse(Path(native.__file__).read_text(encoding='utf8'))
    cls=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name=='FactoredOriginalTransition')
    fn=copy.deepcopy(next(node for node in cls.body if isinstance(node,ast.FunctionDef) and node.name=='source'))
    before=ast.dump(fn)
    assert isinstance(fn.body[2],ast.If) and ast.unparse(fn.body[2].test)=="geometry['kind'] == 'xi'"
    # Only the typed coordinate constructor changes. The original slope_mu,
    # physical operators and every inherited source assignment stay copied.
    fn.body[2]=ast.parse("s=geometry['original_s_source']").body[0]
    env=dict(vars(native));env.update(original_sigma_rows=closed_sigma_rows,safe_transition_kernels=kernel)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<same original source with closed typed coordinate>','exec'),env)
    proof=dict(source_module=Path(native.__file__).name,source_sha256=sha(Path(native.__file__).name),
        copied_source_wrapper_AST_sha256=hashlib.sha256(before.encode()).hexdigest(),
        only_source_coordinate_constructor_replaced=True,
        defining_native_slope_mu_and_physical_programs_unmodified=True,
        original_safe_kernel_program=kproof)
    return env['source'],proof


class CompletePrefixCoordinates(previous.OriginalTransitionCoordinates):
    def geometry(self,kind,left,right):
        if kind in ('xi','k'):
            got=super().geometry(kind,left,right)
            c=self.c;t=self.coordinates.scalar(1)
            if kind=='xi':
                interval=c.mpf((ep(self.cv(got['left']))[0],ep(self.cv(got['right']))[1]))
                got['original_s_source']=self.W*interval
            else:
                kval=c.mpf((ep(self.cv(got['right']))[0],ep(self.cv(got['left']))[1]))
                got['original_s_source']=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=-c.ln(self.D+kval)/2),1,t.ledger)
            return got
        if kind!='mixed':raise ValueError('Original xi/k or exact mixed endpoints required')
        a=Fraction(left);b=Fraction(right);c=self.c;t=self.coordinates.scalar(1)
        if not 0<a<1 or b<=0:raise ValueError('Original mixed bulk-to-seam interval')
        alpha=self.cv(a);beta=1/c.sqrt(1+self.cv(b)/self.D)
        if ep(beta-alpha)[0]<=0:raise ArithmeticError('Original positive mixed-coordinate width')
        width=self.W*(beta-alpha)
        ac=self.W_cover*alpha;bc=1/c.sqrt(self.D+self.cv(b))
        coordinate=c.mpf((ep(ac)[0],ep(bc)[1]))
        original_s=self.W*c.mpf((ep(alpha)[0],ep(beta)[1]))
        k_cover=c.mpf((ep(self.cv(b))[0],ep(self.D*(1/alpha**2-1))[1]))
        record=dict(chart='O3_slope_mu',source_family=self.coordinates.family,
            typed_coordinate_kind=kind,exact_normalized_endpoints=[str(a),str(b)],
            original_endpoint_specification=[dict(kind='xi',exact_multiplier=str(a),source_expression='s=W*xi'),
                dict(kind='k',exact_offset=str(b),source_expression='s=(D+k)**(-1/2)')],
            original_parameter_source=self.parameter_record,native_coordinate_box=coordinate,
            positive_true_log_radius_width=width.record(),original_mixed_s_component=original_s.record(),
            original_correlated_k_cover=k_cover,
            exact_width_identity='W*((1+k_right/D)^(-1/2)-xi_left)',
            regular_true_log_radius_width_cover=c.mpf(0),scalar_width_cover_used_only_for_directed_kernel_bounds=native.bounded(width),
            endpoint_difference_collected_symbolically_before_numeric_subtraction=True,
            integration_Jacobian_installed_once_in_true_width=True,ordinary_source_rows_are_y_derivatives_no_second_coordinate_Jacobian=True,
            width_and_endpoints_independent_of_Z=True,source_endpoints_not_scalar_parameter_selections=True)
        return dict(record=record,kind=kind,left=a,right=b,coordinate=coordinate,width=width,
            regular=c.mpf(0),scalar_cover=native.bounded(width),k_cover=k_cover,
            inverse_square=self.D+k_cover,original_s_source=original_s)

    def normalized_excess(self,template,geometry):
        c=self.c
        if geometry['kind']=='xi' and geometry['left']==0:
            S=c.mpf(ep(geometry['coordinate'])[1]);F=1/(1-S)**2
            anchor=c.ln(2)+self.D+F-1/S**2
            got=weighted_flat_endpoint_rows(template,S,anchor)
            rows={k:got['values'][k[0]] if not k[1] else template.scalar(0) for k in slow.ORDERS}
            return dict(rows=rows,record=dict(status='enclosed',branch='closed_original_endpoint',
                original_Delta_over_eta_source='2*mu*sigma(s)/eta',
                exact_s_zero_excess_and_all_slow_jets_zero=True,
                original_closed_endpoint_weighted_source=got['record'],
                original_normalized_excess_ordinary_rows={'y%d_Z%d'%k:v.record() for k,v in rows.items()},
                original_D_mu_eta_correlation_bound_with_entire_source_covers=True))
        if geometry['kind']!='mixed':return super().normalized_excess(template,geometry)
        s=geometry['coordinate'];F=1/(1-s)**2;k=geometry['k_cover']
        L=F-self.D-k;tail=template.bounded_exp(L);complement=1/(1+tail)
        G=c.ln(2)+F-k
        value=prior.ScaledEnclosure(prior.FormalScale(template.scale.bases,offset=G),complement,template.ledger)
        L1=2/(1-s)**3+2/s**3;L2=6/(1-s)**4-6/s**4
        ordinary=[value,value*complement*L1,value*complement*(L2+(2*complement-1)*L1**2)]
        rows={order:ordinary[order[0]] if not order[1] else template.scalar(0) for order in slow.ORDERS}
        return dict(rows=rows,record=dict(status='enclosed',branch='original_mixed_bulk_to_seam',
            original_source_identity='log(Delta/eta)=log(2)+F-k-log(1+exp(F-D-k))',
            original_correlated_k_cover=k,positive_original_logistic_denominator_inverse=complement,
            original_logistic_odds_cover=L,original_normalized_excess_log_numerator_cover=G,
            original_normalized_excess_ordinary_rows={'y%d_Z%d'%order:v.record() for order,v in rows.items()},
            shared_D_terms_cancelled_before_evaluation=True,ordinary_y_derivatives_no_mixed_coordinate_Jacobian=True))


def recorded_native_adapter(data,coordinates):
    c=coordinates.ctx;row=data['original_source_rows'][0]['original_typed_native_source']
    inputs=row['lossless_actual_native_parent_inputs'];Z=tuple(data['exact_Z_range'])
    parent=dict(actual_normalized_primitive_y_derivative_axial5={key:[IntervalTaylor(c,[iv(c,v) for v in vals]) for vals in rows]
        for key,rows in inputs['original_history_axial_Taylor_coefficients'].items()},
        original_P0_axial5_coefficients=[iv(c,v) for v in inputs['original_P0_axial_Taylor_coefficients']],
        log_Utheta_over_Pstar_base_source=iv(c,inputs['original_log_Utheta_over_Pstar_base_source']))
    adapter=object.__new__(native.FactoredOriginalTransition)
    adapter.c=c;adapter.coordinates=coordinates;adapter.template=coordinates.scalar(1)
    adapter.mu=iv(c,inputs['original_mu']);adapter.eta=iv(c,inputs['original_eta_log']);adapter.owner=SimpleNamespace(family=coordinates.family)
    adapter.proof=native.canonical_amplitude_proof();adapter.parents={Z:parent}
    adapter.pre=SimpleNamespace(params=SimpleNamespace(mu=adapter.mu),delta=iv(c,inputs['original_delta']))
    fn,proof=source_programs()
    adapter.source=lambda target,geometry,cells=16:fn(adapter,target,geometry,cells)
    adapter.closed_source_program=proof
    return adapter


def saved_original_phase(data,coordinates,geometry):
    c=coordinates.ctx;t=coordinates.scalar(1)
    actual=data['original_source_rows'][0]['actual_original_typed_phase'];base=actual['original_base_source_phase']
    if base['periodic_projection']['full_period']:boxes=[c.mpf((0,1))]
    else:
        boxes=[]
        for box in base['periodic_projection']['boxes']:
            boxes.extend(first.spatial.ordinary_mod_one(c,iv(c,box)+geometry['coordinate']*N)['boxes'])
    raw=base['actual_log_R_over_same_r_minus_enclosure'];scale=raw['formal_positive_scale']
    assert scale['source_exponents']==[0,0,0,0] and not scale['radius_power']
    origin=prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,offset=iv(c,scale['additional_log_interval'])),iv(c,raw['coefficient_interval']),t.ledger)
    micro=geometry['original_s_source']
    return dict(original_base_source_phase=base,typed_source_geometry=geometry['record'],
        original_microscopic_s_component=micro.record(),original_microscopic_N_s_component=(micro*N).record(),
        actual_original_log_R_over_r_minus=(origin+micro).record(),actual_phase_boxes=boxes,
        same_original_affine_radius_identity=actual['same_original_affine_radius_identity'],
        actual_phase_is_N_times_original_s_not_N_times_xi_or_k=True,nonzero_microscopic_phase_component_retained=not micro.zero,
        phase_independent_of_Z=True,phase_box_not_an_independently_selected_sample=True,
        genuine_native_origin_reused_with_no_midpoint_or_new_field=True)


def quiet_geometry(coordinates,coord,left_micro,right):
    c=coordinates.ctx;t=coordinates.scalar(1);r=Fraction(right);rr=c.mpf(r.numerator)/r.denominator
    width=t.scalar(rr)-left_micro
    if ep(width.coefficient)[0]<=0:raise ArithmeticError('Original positive support-to-quiet span')
    cover=native.bounded(width)
    record=dict(chart='O3_slope_mu',source_family=coordinates.family,
        original_endpoint_specification=[dict(kind='k',exact_offset='0',source_expression='s=W=D^(-1/2)'),str(r)],
        exact_width_expression='right-W',original_left_endpoint_nonzero_W=left_micro.record(),
        positive_true_log_radius_width=width.record(),regular_true_log_radius_width_cover=rr,
        scalar_width_cover_used_only_for_directed_kernel_bounds=cover,
        original_support_theorem=support.original_active_support_proof(coordinates,coord.mu,coord.eta_log),
        width_and_endpoints_independent_of_Z=True,endpoint_subtraction_not_used_to_erase_W=True)
    return dict(width=width,regular=rr,scalar_cover=cover,record=record)


PLAN=(('closed_endpoint','xi','0','1/4'),('bulk','xi','1/4','3/4'),
      ('bulk_to_seam','mixed','3/4','4'),('seam_outer','k','4','2'),
      ('seam_active','k','2','7/4'),('cutoff_seam','k','7/4','5/3'),('flat_seam','k','5/3','0'))


def original_endpoint_telescope():
    D=sy.symbols('D',positive=True);W=1/sy.sqrt(D);stop=sy.Rational(1,2**128)
    points=[sy.Integer(0),W/4,3*W/4,1/sy.sqrt(D+4),1/sy.sqrt(D+2),
        1/sy.sqrt(D+sy.Rational(7,4)),1/sy.sqrt(D+sy.Rational(5,3)),W,stop,sy.Integer(1)]
    widths=[b-a for a,b in zip(points,points[1:])]
    assert sy.simplify(sum(widths[:8])-stop)==0 and sy.simplify(sum(widths)-1)==0
    return dict(passed=True,original_parameter_symbol='D=log(mu)-log(eta)>0',
        original_ordered_source_endpoints=[str(p) for p in points],
        exact_source_cell_width_expressions=[str(w) for w in widths],
        exact_prefix_total_width=str(stop),exact_full_transition_total_width='1',
        exact_Rd_to_Rc_total_width='3',
        adjacent_source_endpoint_identities_not_scalar_endpoint_overlap_used=True,
        mixed_width_positive_if_D_exceeds='36/7',actual_width_order_checked_in_geometry=True)


def incoming_at_Rd(coordinates,Z,report):
    target=tuple(Fraction(v) for v in Z)
    row=next(r for r in report['actual_refined_O2_to_Rc_replays']
        if r['original_O2_ordered_source_cells']==2048 and tuple(Fraction(v) for v in r['exact_Z_range'])==target)
    assert row['candidate_N']==N
    before=row['six_original_cell_actual_histories'][2];transition=row['six_original_cell_actual_histories'][3]
    assert before['label']=='buffer_9_11' and transition['label']=='transition'
    assert before['actual_right_correction_C0']==transition['actual_inherited_C0']
    assert before['actual_right_correction_Z']==transition['actual_inherited_Z']
    values={k:coordinates.rebase(common.restore_common_source(v,rc.history.CommonSourceCoordinates(coordinates.ctx,coordinates.logP_squared,coordinates.family)),coordinates.family)
        for k,v in before['actual_right_correction_C0'].items()}
    jets={k:coordinates.rebase(common.restore_common_source(v,rc.history.CommonSourceCoordinates(coordinates.ctx,coordinates.logP_squared,coordinates.family)),coordinates.family)
        for k,v in before['actual_right_correction_Z'].items()}
    archive=row['actual_six_cell_source_binding']['archive']
    assert sha(archive)==row['actual_six_cell_source_binding']['sha256']
    source_data=json.loads(gzip.decompress((HERE/archive).read_bytes()))
    power=source_data['ordered_original_six_cell_source_records'][4:]
    assert len(power)==2 and all(r['original_whole_power_q_q_Z_and_density_C0_Z_exact_zero'] for r in power)
    return dict(values=values,Z_derivatives=jets,accepted_Rc_background=row['original_Rc_background_and_separate_P0_Z'],
        original_quiet_power_source_records=power,record=dict(source_manifest=tail.NAME,source_manifest_sha256=sha(tail.NAME),
        original_route_endpoint='Rd=O2_buffer(offset11)=O3_slope_mu(offset0)',exact_Z_range=list(Z),candidate_N=N,
        source_family=coordinates.family,ordered_source_cells=2048,source_cell_index=2,
        original_source_cell_record=before,following_original_transition_inherited_record=transition,
        actual_incoming_correction_C0=native.records(values),actual_incoming_correction_Z=native.records(jets),
        original_background_not_used_as_correction=True,actual_Rc_downstream_values_not_used_as_Rd_inlet=True,
        genuine_upstream_function_graph=row['genuine_slope_exit_binding'],
        actual_six_original_cell_source_binding=row['actual_six_cell_source_binding']))


def execute_tile(data,coordinates,tail_report):
    c=coordinates.ctx;adapter=recorded_native_adapter(data,coordinates)
    coord=CompletePrefixCoordinates(coordinates,adapter.mu,adapter.eta)
    operator=rc.history.C1DuhamelOperator(coordinates);records=[]
    seed=SimpleNamespace(core=SimpleNamespace(logC=iv(c,data['original_logC'])))
    for label,kind,left,right in PLAN:
        geometry=coord.geometry(kind,left,right);got=adapter.source(data['exact_Z_range'],geometry)
        norm=coord.normalized_excess(adapter.template,geometry)
        roots=native.signed_root_source(adapter,got,norm['rows'],adapter.eta,iv(c,data['original_positive_a_lower_log']),seed)
        phase=saved_original_phase(data,coordinates,geometry)
        local=native.original_local_density_integrals(adapter,got,roots,geometry,phase,iv(c,data['original_dstar_log']))
        if local['record']['status']!='enclosed':raise ArithmeticError('Original closed-prefix source unresolved: '+label)
        restore=lambda rec:restore_half_source(rec,coordinates)
        values={k:restore(v) for k,v in local['record']['original_local_five_signed_C0_Duhamel_contributions'].items()}
        jets={k:restore(v) for k,v in local['record']['original_local_five_signed_Z_Duhamel_contributions'].items()}
        rc.transfer.append_true_cell(operator,geometry,values,jets,coordinates.family)
        records.append(dict(label=label,actual_true_geometry=geometry['record'],
            original_typed_native_source=got['record'],original_closed_normalized_excess=norm['record'],
            original_signed_root_source=roots['record'],actual_original_typed_phase=phase,
            actual_original_conditioned_first_jet_queries=local['inverses'],
            actual_original_local_five_C0_Z_density_and_integral=local['record']))
        print('Original complete prefix:',data['exact_Z_range'],label,'enclosed',flush=True)
    quiet=quiet_geometry(coordinates,coord,coord.W,Fraction(1,2**128))
    zero={k:coordinates.scalar(0) for k in KEYS}
    rc.transfer.append_true_cell(operator,quiet,zero,zero,coordinates.family)
    records.append(dict(label='support_to_fixed_quiet_collar',actual_true_geometry=quiet['record'],
        original_q_and_all_required_slow_jets_exact_zero=True,
        all_five_original_own_density_C0_Z_and_integral_increments_exact_zero=True,
        incoming_history_not_zeroed_or_replaced=True))
    incoming=incoming_at_Rd(coordinates,data['exact_Z_range'],tail_report)
    # Endpoint expressions telescope to 2^-128 exactly. Collect this original
    # affine coefficient before independent directed endpoint uncertainty.
    stop=c.mpf(1)/2**128
    operator.coefficients={k:coordinates.decay(stop,r) for k,r in rc.RATES.items()}
    prefix=operator.apply(incoming['values'],incoming['Z_derivatives'],coordinates.family)
    prefix_record=dict(operator.record(),exact_original_total_width='2^-128',
        exact_source_endpoint_telescope_collected_before_coefficient_evaluation=True,
        actual_incoming= incoming['record'],actual_outgoing_correction_C0=native.records(prefix['values']),
        actual_outgoing_correction_Z=native.records(prefix['Z_derivatives']))
    # The checked original [2^-128,1] quiet theorem supplies zero source
    # increments and ordinary attenuation; all preceding memory is retained.
    collar=rc.history.C1DuhamelOperator(coordinates)
    collar.coefficients={k:coordinates.decay(c.mpf(1)-stop,r) for k,r in rc.RATES.items()}
    outgoing=collar.apply(prefix['values'],prefix['Z_derivatives'],coordinates.family)
    total=rc.history.C1DuhamelOperator(coordinates)
    total.coefficients={k:coordinates.decay(1,r) for k,r in rc.RATES.items()}
    total.increments={k:collar.coefficients[k]*operator.increments[k] for k in KEYS}
    total.Z_increments={k:collar.coefficients[k]*operator.Z_increments[k] for k in KEYS}
    recomposed=total.apply(incoming['values'],incoming['Z_derivatives'],coordinates.family)
    power=rc.history.C1DuhamelOperator(coordinates)
    power.coefficients={k:coordinates.decay(2,r) for k,r in rc.RATES.items()}
    at_Rc=power.apply(recomposed['values'],recomposed['Z_derivatives'],coordinates.family)
    final=rc.history.C1DuhamelOperator(coordinates)
    final.coefficients={k:coordinates.decay(3,r) for k,r in rc.RATES.items()}
    final.increments={k:power.coefficients[k]*total.increments[k] for k in KEYS}
    final.Z_increments={k:power.coefficients[k]*total.Z_increments[k] for k in KEYS}
    collected_Rc=final.apply(incoming['values'],incoming['Z_derivatives'],coordinates.family)
    background=incoming['accepted_Rc_background']
    original_coordinates=rc.history.CommonSourceCoordinates(c,coordinates.logP_squared,coordinates.family)
    lift=lambda rec:coordinates.rebase(common.restore_common_source(rec,original_coordinates),coordinates.family)
    originals={k:lift(v) for k,v in background['original_normalized_history_C0_enclosures'].items()}
    originalZ={k:lift(v) for k,v in background['original_normalized_history_Z_enclosures'].items()}
    own={k:originals[k]+collected_Rc['values'][k] for k in KEYS}
    ownZ={k:originalZ[k]+collected_Rc['Z_derivatives'][k] for k in KEYS}
    P0=lift(background['original_separate_P0_over_Pstar_squared']);P0Z=lift(background['original_separate_P0_Z_over_Pstar_squared'])
    return dict(source_family=coordinates.family,candidate_N=N,exact_Z_range=data['exact_Z_range'],
        reused_native_parent_source=Path(native.__file__).name,closed_typed_source_program=adapter.closed_source_program,
        common_directed_coordinate_theorem=coordinates.record(),original_complete_prefix_source_rows=records,
        original_endpoint_order_and_exact_width_telescope=original_endpoint_telescope(),
        original_prefix_operator_and_genuine_incoming=prefix_record,
        accepted_quiet_collar_source_manifest=support.NAME,accepted_quiet_collar_receipt=support.RECEIPT,
        original_complete_transition_operator=total.record(),
        actual_full_transition_outgoing_correction_C0=native.records(recomposed['values']),
        actual_full_transition_outgoing_correction_Z=native.records(recomposed['Z_derivatives']),
        separately_composed_outgoing_correction_C0=native.records(outgoing['values']),
        separately_composed_outgoing_correction_Z=native.records(outgoing['Z_derivatives']),
        original_quiet_power_source_records=incoming['original_quiet_power_source_records'],
        original_Rd_to_Rc_operator=final.record(),original_Rc_endpoint='Rc=Rw*exp(2)=Rd*exp(3)',
        actual_updated_Rc_correction_C0=native.records(collected_Rc['values']),
        actual_updated_Rc_correction_Z=native.records(collected_Rc['Z_derivatives']),
        separately_composed_Rc_correction_C0=native.records(at_Rc['values']),
        separately_composed_Rc_correction_Z=native.records(at_Rc['Z_derivatives']),
        original_Rc_background_and_separate_P0_Z=background,
        actual_updated_Rc_own_history_C0=native.records(own),actual_updated_Rc_own_history_Z=native.records(ownZ),
        actual_updated_Rc_absolute_pressure=native.records(P0+own['p']),
        actual_updated_Rc_absolute_pressure_Z=native.records(P0Z+ownZ['p']),
        actual_original_pressure_datum_kept_separate_and_added_once=True,
        full_original_transition_interval_zero_to_one_covered=True,actual_upstream_correction_not_reset=True,
        actual_Rc_targets_controls_or_global_N_or_recursion_admitted=False)


def restore_half_source(record,coordinates):
    c=coordinates.ctx;scale=record['formal_positive_scale']
    return prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,
        tuple(scale['source_exponents'])+(scale['radius_power'],),iv(c,scale['additional_log_interval'])),
        iv(c,record['coefficient_interval']),coordinates.ledger)


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();base=json.loads((HERE/native.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,native,base['source_family'])
    tail_report,_=common.attach_receipt(hashes,tail,base['source_family'])
    common.attach_receipt(hashes,support,base['source_family'])
    from mpmath.ctx_iv import MPIntervalContext
    c=MPIntervalContext();c.dps=240;archives=[]
    for archive in base['original_factored_native_source_archives']:
        raw=gzip.decompress((HERE/archive['filename']).read_bytes())
        assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
        data=json.loads(raw);coordinates=native.HalfPstarCoordinates(c,iv(c,data['common_directed_coordinate_theorem']['common_log_bases'][1]),base['source_family'])
        result=execute_tile(data,coordinates,tail_report)
        result.update(actual_native_source_input_archive=archive,
            original_logC=data['original_logC'],original_dstar_log=data['original_dstar_log'],original_positive_a_lower_log=data['original_positive_a_lower_log'])
        encoded=json.dumps(encode(result),indent=2).encode()+b'\n';compressed=gzip.compress(encoded,compresslevel=9,mtime=0)
        tag='positive' if Fraction(data['exact_Z_range'][0])>0 else 'negative'
        filename=PREFIX+'current_transition_complete_prefix_'+tag+'.json.gz'
        (HERE/filename).write_bytes(compressed);hashes[filename]=sha(filename)
        archives.append(dict(filename=filename,compressed_bytes=len(compressed),uncompressed_bytes=len(encoded),
            lossless_original_json_sha256=hashlib.sha256(encoded).hexdigest(),exact_Z_range=data['exact_Z_range']))
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=base['source_family'],candidate_N=N,
        actual_original_complete_transition_archives=archives,ordered_original_active_prefix_source_cells=len(PLAN),
        closed_original_endpoint_no_epsilon_gap=True,full_original_transition_C0_Z_operator_with_genuine_Rd_incoming=True,
        updated_actual_Rc_correction_and_own_history_C0_Z_covers_executed=True,
        native_upstream_owner_or_accepted_quadrature_producers_not_rerun=True,
        actual_full_axial_buffer_numerical_integrals_or_Rc_targets_controls_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,
        execution_seconds=time.monotonic()-began,
        scope='Original closed endpoint, mixed bulk-to-seam, all cutoff collars and quiet support compose complete transition [0,1] C0/Z operator, with genuine accepted Rd correction on two Z tiles at N1024. Exact quiet power updates actual Rc correction/background/pressure covers. No complete axial/buffer numerical source, terminal targets/controls/global N/heat/stress/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

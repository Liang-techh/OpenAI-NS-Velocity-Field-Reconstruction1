"""Same current Banach limit: actual implicit Z2 controls and C2 Picard tail.

The defining functions consume signed target handles. Quantitative target
caps certify their regularity and tails, and never become graph values.
The selected frequency remains the previously accepted repair-only integer.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C2_target_ranges as ranges

phase=ranges.phase;current=ranges.current;limit=phase.limit
source,controls=current.source,current.controls
HERE,PREFIX,sha,ep=ranges.HERE,ranges.PREFIX,ranges.sha,ranges.ep
NAME=PREFIX+'current_original_C2_limit_controls.json.gz'
RECEIPT=PREFIX+'current_original_C2_limit_controls_check.json'
GATES=('current_original_same_Banach_limit_C2_controls_installed',
    'current_original_same_repair_integer_quantitative_C2_Picard_tail_installed')
OPEN=ranges.OPEN


def encoded(value):
    if isinstance(value,source.C1Function):return current.encode_graph(value)
    if isinstance(value,dict):return {k:encoded(v) for k,v in value.items() if k!='graph'}
    if isinstance(value,(list,tuple)):return [encoded(v) for v in value]
    return phase.encoded(value)


def quadratic_ZZ(g,W,h):
    """DQ(h) h_ZZ + D2Q[h_Z,h_Z]; ordinary second derivatives."""
    first=controls.quadratic_action(g,W,[q.value for q in h],[q.ZZ for q in h])
    second=controls.quadratic_action(g,W,[q.Z for q in h],[q.Z for q in h])
    return [g.add(a,b) for a,b in zip(first,second)]


def picard_C2_map(g,W,N,d,h):
    overN=lambda q:g.quotient(q,N,'same selected exact positive integer N')
    Q=controls.quadratic(g,W,[q.value for q in h])
    QZ=controls.quadratic_action(g,W,[q.value for q in h],[q.Z for q in h])
    QZZ=quadratic_ZZ(g,W,h)
    rows=[]
    for key,nonlinear in (('value',Q),('Z',QZ),('ZZ',QZZ)):
        rhs=[g.add(getattr(q,key),overN(r)) for q,r in zip(d,nonlinear)]
        rows.append([g.neg(v) for v in controls.inverse_action(g,W,rhs)])
    return [phase.C2Function(*row) for row in zip(*rows)]


def quantitative_bounds(range_field,range_report):
    c=range_field.c;bd=ranges.Bounds(c)
    report=range_field.phase.report;frequency=report['actual_whole_Z_frequency_connection']
    conditions=frequency['repair_conditions']
    read=lambda key:ranges.LogUpper(c,current.packets.interval(c,conditions[key]['log_absolute_upper']))
    CA,CQ,rho=[read(key) for key in ('exact_integral_matrix_inverse_log_cap',
        'general_transformed_quadratic_C1_log_cap','formal_control_C1_ball_radius_log')]
    selection=frequency['selected_integer'];N=current.DoubleDyadicInteger(selection['J'])
    logN=N.log_interval(c);epsilon=ranges.LogUpper(c,-c.mpf(ep(logN)[0]))
    maximum=max(ep(current.packets.interval(c,row['actual_five_target_ZZ_maximum']['log_absolute_upper']))[1]
        for row in range_report['actual_four_Z_C2_target_transports'])
    D2=ranges.LogUpper(c,c.mpf(maximum))
    perturbation=bd.scale(CA*CQ*rho*epsilon,2)
    if ep(perturbation.log+c.ln(2))[1]>0:raise ValueError('Same original half-contraction required')
    inverse=bd.scale(CA,2)
    second= inverse*bd.sum(D2,bd.scale(CQ*bd.power(rho,2)*epsilon,2))
    # e2_(n+1) <= e2_n/2 + A*rho*2^-n, e2_0<=H2.
    forcing=CA*CQ*epsilon*bd.sum(bd.scale(second,2),bd.scale(rho,4))
    return dict(actual_same_repair_integer=selection,
        actual_log_N=logN,actual_inverse_epsilon_upper=epsilon.record(),
        original_B_inverse_upper=CA.record(),original_quadratic_upper=CQ.record(),
        original_C1_ball_radius_upper=rho.record(),actual_target_ZZ_global_upper=D2.record(),
        actual_J_preconditioned_perturbation_upper=perturbation.record(),
        actual_J_inverse_upper=inverse.record(),actual_limit_and_iterate_ZZ_upper=second.record(),
        actual_C2_error_forcing_upper=forcing.record(),
        actual_C2_error_forcing_times_rho_upper=(forcing*rho).record(),
        uniform_Z_domain=[-1,1],actual_N_not_reselected=True,
        inverse_theorem='J=B+DQ(h)/N; ||B^-1 DQ(h)/N||<=1/2 for every h in original C1 ball; ||J^-1||<=2CA',
        quadratic_theorem='Original positive integral coefficient absolute row sums give ||DQ(h)v||<=2CQ||h||||v|| and ||D2Q[v,w]||<=2CQ||v||||w||',
        C2_theorem='The signed target is C2 and J is invertible on the original unique C1 branch. IFT gives the same h* in C2, not a new branch.',
        second_bound='H2=2CA*(D2+2CQ*rho^2/N), valid for h*_ZZ and every exact C2 Picard iterate',
        error_recurrence='e2_(n+1)<=e2_n/2+A*rho*2^-n; A=CA*CQ*(2H2+4rho)/N; e2_0<=H2',
        actual_C2_tail='||h*_ZZ-h_n_ZZ||sup <= (H2+2*n*A*rho)*2^-n',
        C1_tail_unchanged='||h*-h_n||C1<=rho*2^-n',
        small_relative_tail_factors_kept_separate_from_large_logarithms=True,
        quadrature_phase_and_oracle_errors_not_included=True,
        physical_field_global_frequency_stress_and_recursion_not_admitted=True)


def factored_C2_tail(bounds,depth):
    if type(depth) is not int or not 0<=depth<=4096:raise ValueError('Exact tail depth in[0,4096] required')
    return dict(depth=depth,second_error_terms=[
        dict(log_cap=bounds['actual_limit_and_iterate_ZZ_upper'],exact_relative_numerator=1,
            exact_relative_denominator=2**depth),
        dict(log_cap=bounds['actual_C2_error_forcing_times_rho_upper'],exact_relative_numerator=2*depth,
            exact_relative_denominator=2**depth)],
        original_C1_error=dict(log_cap=bounds['original_C1_ball_radius_upper'],
            exact_relative_numerator=1,exact_relative_denominator=2**depth),
        second_rows_are_genuine_derivatives_of_exact_Picard_functions=True,
        numerical_field_evaluation_or_oracle_tail_not_installed=True)


def build(range_field,bounds):
    phase_field=range_field.phase;g=phase_field.built['graph']
    old=limit.exact_limit_adapter(phase_field.built,phase_field.report)
    W=old['repair_weights'];h1=old['exact_limit_controls'];d=list(range_field.target.target_functions().values())
    if tuple(range_field.target.target_functions())!=controls.ROWS:raise ValueError('Original five target row order required')
    overN=lambda q:g.quotient(q,old['N'],'same selected exact positive integer N')
    second=controls.quadratic_action(g,W,[q.Z for q in h1],[q.Z for q in h1])
    rhs=[g.neg(g.add(q.ZZ,overN(v))) for q,v in zip(d,second)]
    system=g.node('implicit_C2_ZZ_linear_system',matrix=current.encode_graph(old['implicit_Z_matrix']),
        rhs=current.encode_graph(rhs),limit_vector=old['exact_limit_vector'].node,
        original_C1_linear_system=old['implicit_Z_linear_system'].node,
        definition='(B+DQ(h*)/N)*h*_ZZ=-d_ZZ-D2Q[h*_Z,h*_Z]/N',
        inverse_and_C2_proof=phase.encoded(bounds),
        signed_target_C2_source=ranges.target.RECEIPT,signed_target_C2_source_sha256=sha(ranges.target.RECEIPT),
        quantitative_target_source=ranges.RECEIPT,quantitative_target_source_sha256=sha(ranges.RECEIPT),
        same_original_C1_limit_branch=True,numerical_value_not_installed=True)
    h=[phase.C2Function(q.value,q.Z,g.node('C2_Banach_limit_second_component',
        limit_vector=old['exact_limit_vector'].node,component=i,Z_order=2,
        implicit_ZZ_linear_system=system.node,genuine_ordinary_second_derivative=True,
        derivative_of_range_endpoint=False)) for i,q in enumerate(h1)]
    QZZ=quadratic_ZZ(g,W,h)
    residual=[g.add(q.ZZ,g.add(*(g.mul(b,v.ZZ) for b,v in zip(row,h))),overN(r))
        for q,row,r in zip(d,W['B'],QZZ)]
    zeros=[g.node('proved_C2_second_zero_identity',expression=v.node,
        theorem='twice differentiate same original signed residual and use actual implicit C2 system',
        implicit_ZZ_linear_system=system.node,zero=g.zero.node,row=row)
        for row,v in zip(controls.ROWS,residual)]
    alg=phase.C2Algebra(g);sequence=[[alg.fixed(0) for unused in range(5)]]
    for unused in range(3):sequence.append(picard_C2_map(g,W,old['N'],d,sequence[-1]))
    return dict(exact_C1_limit_adapter=old,exact_limit_controls_C2=h,
        actual_signed_targets_C2=dict(zip(controls.ROWS,d)),same_implicit_matrix=old['implicit_Z_matrix'],
        actual_implicit_ZZ_rhs=rhs,implicit_ZZ_linear_system=system,
        actual_second_control_residual=residual,actual_second_zero_theorem_nodes=zeros,
        genuine_C2_Picard_sequence=sequence,actual_quantitative_C2_bounds=bounds,
        existing_C1_value_and_Z_handles_retained=True,
        original_bump_weights_mu_and_N_Z_independent=True,
        C2_repair_band_fields_and_histories_not_yet_installed=True)


class CurrentC2LimitControls:
    def __init__(self,range_field=None,owner=None,require_checked=True):
        self.ranges=range_field if range_field is not None else ranges.CurrentC2TargetRanges(owner=owner)
        if not self.ranges.acceptance_loaded:raise ValueError('Accepted quantitative C2 target ranges required')
        self.identity=self.ranges.identity;self.hashes=dict(self.ranges.hashes)
        parent=json.loads((HERE/limit.RECEIPT).read_bytes())
        if not parent['all_passed'] or not parent[limit.GATE]:raise ValueError('Accepted same C1 limit band required')
        for name,digest in parent['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed original C1 limit source '+name)
        self.hashes.update(parent['input_hashes'])
        for name in (ranges.NAME,ranges.RECEIPT,limit.NAME,limit.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        report=json.loads(gzip.decompress((HERE/ranges.NAME).read_bytes()))
        if report['source_family']!=self.identity:raise ValueError('Same current C2 source family required')
        self.prefix=[dict(n) for n in self.ranges.phase.built['graph'].nodes]
        self.bounds=quantitative_bounds(self.ranges,report);self.functions=build(self.ranges,self.bounds)
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual C2 controls required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual C2 control source '+name)
            self.acceptance_loaded=True

    def control_functions(self):return self.functions['exact_limit_controls_C2']
    def tail(self,depth):return factored_C2_tail(self.bounds,depth)


def run(range_field=None,owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC2LimitControls(range_field=range_field,owner=owner,require_checked=False)
        report=dict(candidate_actual_C2_limit_controls_constructed=True,source_family=field.identity,
            original_C2_target_graph_prefix_length=len(field.prefix),
            exact_graph_nodes=field.ranges.phase.built['graph'].nodes,
            actual_C2_limit_controls_and_Picard_functions=encoded(field.functions),
            actual_quantitative_C2_bounds=phase.encoded(field.bounds),
            actual_C2_tail_examples=[field.tail(n) for n in (0,3,32)],
            source_bindings=dict(quadratic_ZZ=current.ast_binding(quadratic_ZZ),
                picard_C2_map=current.ast_binding(picard_C2_map),build=current.ast_binding(build),
                quantitative_bounds=current.ast_binding(quantitative_bounds),
                original_C1_limit_adapter=current.ast_binding(limit.exact_limit_adapter),
                original_exact_weights=current.ast_binding(controls.exact_weights),
                original_quadratic=current.ast_binding(controls.quadratic),
                original_quadratic_action=current.ast_binding(controls.quadratic_action),
                original_B_inverse_action=current.ast_binding(controls.inverse_action)),
            C2_repair_band_fields_and_histories_installed=False,
            current_numeric_point_field_oracle_installed=False,**dict.fromkeys(GATES+OPEN,False),
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(phase.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Same original repair-only Banach limit C2 controls and genuine C2 Picard tail constructed',flush=True)
    return field


if __name__=='__main__':run()

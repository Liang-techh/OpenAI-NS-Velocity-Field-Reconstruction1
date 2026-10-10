"""Same original C2 repaired limit: signed third controls and Picard tails."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C3_target_ranges as ranges
import lei_ren_part1_paper_compliant_current_original_C2_limit_controls as lower

phase,current,source,controls=ranges.phase,ranges.current,lower.source,lower.controls
target=ranges.target
HERE,PREFIX,sha,ep=ranges.HERE,ranges.PREFIX,ranges.sha,ranges.ep
NAME=PREFIX+'current_original_C3_limit_controls.json.gz'
RECEIPT=PREFIX+'current_original_C3_limit_controls_check.json'
GATES=('current_original_same_Banach_limit_C3_controls_installed',
    'current_original_same_repair_integer_quantitative_C3_Picard_tail_installed')
OPEN=ranges.OPEN


def encoded(value):
    if isinstance(value,target.C3Function):return target.encoded(value)
    if isinstance(value,dict):return {k:encoded(v) for k,v in value.items() if k!='graph'}
    if isinstance(value,(tuple,list)):return [encoded(v) for v in value]
    return lower.encoded(value)


def quadratic_ZZZ(g,W,h):
    first=controls.quadratic_action(g,W,[q.value for q in h],[q.ZZZ for q in h])
    cross=controls.quadratic_action(g,W,[q.Z for q in h],[q.ZZ for q in h])
    return [g.add(a,g.mul(g.constant(3),b)) for a,b in zip(first,cross)]


def picard_C3_map(g,W,N,d,h,old):
    nonlinear=quadratic_ZZZ(g,W,h)
    rhs=[g.add(q.ZZZ,g.quotient(v,N,'same selected exact positive integer N')) for q,v in zip(d,nonlinear)]
    third=[g.neg(v) for v in controls.inverse_action(g,W,rhs)]
    return [target.preserve(q,v) for q,v in zip(old,third)]


def quantitative_bounds(control_field,range_field,report):
    c=range_field.c;bd=ranges.Bounds(c);old=control_field.bounds
    read=lambda key:ranges.LogUpper(c,current.packets.interval(c,old[key]['log_absolute_upper']))
    CA,CQ,rho,H2,A=[read(k) for k in ('original_B_inverse_upper','original_quadratic_upper',
        'original_C1_ball_radius_upper','actual_limit_and_iterate_ZZ_upper','actual_C2_error_forcing_upper')]
    selection=old['actual_same_repair_integer'];N=current.DoubleDyadicInteger(selection['J']);logN=N.log_interval(c)
    epsilon=ranges.LogUpper(c,-c.mpf(ep(logN)[0]));inverse=read('actual_J_inverse_upper')
    maximum=max(ep(current.packets.interval(c,row['actual_five_target_ZZZ_maximum']['log_absolute_upper']))[1]
        for row in report['actual_four_Z_C3_target_transports'])
    D3=ranges.LogUpper(c,c.mpf(maximum));H3=inverse*bd.sum(D3,bd.scale(CQ*rho*H2*epsilon,6))
    forcingB=CA*CQ*rho*epsilon*bd.sum(bd.scale(H3,2),bd.scale(H2,12))
    forcingC=bd.scale(CA*CQ*A*bd.power(rho,2)*epsilon,12)
    return dict(actual_same_repair_integer=selection,actual_log_N=logN,actual_inverse_epsilon_upper=epsilon.record(),
        original_B_inverse_upper=CA.record(),original_quadratic_upper=CQ.record(),original_C1_ball_radius_upper=rho.record(),
        actual_J_inverse_upper=inverse.record(),actual_limit_and_iterate_ZZ_upper=H2.record(),
        actual_C2_error_forcing_upper=A.record(),actual_target_ZZZ_global_upper=D3.record(),
        actual_limit_and_iterate_ZZZ_upper=H3.record(),actual_C3_error_forcing_B_upper=forcingB.record(),
        actual_C3_error_forcing_C_upper=forcingC.record(),uniform_Z_domain=[-1,1],actual_N_not_reselected=True,
        same_C1_C2_branch_and_original_invertible_J_retained=True,
        third_bound='H3=2CA*(D3+6CQ*rho*H2/N), for same limit and every exact Picard third iterate',
        derivative_equation='J*h*_ZZZ=-d_ZZZ-3*D2Q[h*_Z,h*_ZZ]/N; J=B+DQ(h*)/N',
        C3_theorem='Signed original target is C3 and the same J is invertible: IFT differentiates the existing C2 branch.',
        error_recurrence='e3_(n+1)<=e3_n/2+(CA*CQ/N)*(2H3*e0_n+6H2*e1_n+6rho*e2_n)',
        lower_errors='e0_n,e1_n<=rho*2^-n; e2_n<=(H2+2n*A*rho)*2^-n',
        actual_C3_tail='e3_n<=(H3+2n*B+n(n-1)*C)*2^-n',
        tail_constants='B=CA*CQ*rho*(2H3+12H2)/N; C=12CA*CQ*A*rho^2/N',
        small_relative_tail_factors_kept_separate_from_large_logarithms=True,
        numerical_quadrature_phase_and_oracle_errors_not_included=True,
        physical_global_frequency_stress_and_recursion_not_admitted=True)


def factored_C3_tail(bounds,depth):
    if type(depth) is not int or not 0<=depth<=4096:raise ValueError('Exact tail depth in[0,4096] required')
    return dict(depth=depth,third_error_terms=[
        dict(log_cap=bounds[key],exact_relative_numerator=numerator,exact_relative_denominator=2**depth)
        for key,numerator in (('actual_limit_and_iterate_ZZZ_upper',1),
            ('actual_C3_error_forcing_B_upper',2*depth),('actual_C3_error_forcing_C_upper',depth*(depth-1)))],
        genuine_third_derivatives_of_same_exact_Picard_functions=True,
        numerical_field_or_oracle_tail_not_installed=True)


def build(field,bounds):
    parent=field.lower.functions;old=parent['exact_C1_limit_adapter'];g=field.ranges.phase.built['graph']
    W=old['repair_weights'];h2=field.lower.control_functions();d=list(field.ranges.target.target_functions().values())
    if tuple(field.ranges.target.target_functions())!=controls.ROWS:raise ValueError('Original row order required')
    for key,q in zip(controls.ROWS,d):
        original=parent['actual_signed_targets_C2'][key]
        if target.rows(q)[:3]!=tuple(vars(original).values()):raise ValueError('Original signed C2 target prefix changed')
    cross=controls.quadratic_action(g,W,[q.Z for q in h2],[q.ZZ for q in h2])
    rhs=[g.neg(g.add(q.ZZZ,g.quotient(g.mul(g.constant(3),v),old['N'],
        'same selected exact positive integer N'))) for q,v in zip(d,cross)]
    system=g.node('implicit_C3_ZZZ_linear_system',matrix=current.encode_graph(parent['same_implicit_matrix']),
        rhs=current.encode_graph(rhs),limit_vector=old['exact_limit_vector'].node,
        original_C1_linear_system=old['implicit_Z_linear_system'].node,
        original_C2_linear_system=parent['implicit_ZZ_linear_system'].node,
        definition='(B+DQ(h*)/N)*h*_ZZZ=-d_ZZZ-3*D2Q[h*_Z,h*_ZZ]/N',
        inverse_and_C3_proof=encoded(bounds),signed_C3_source=target.RECEIPT,signed_C3_source_sha256=sha(target.RECEIPT),
        quantitative_C3_source=ranges.RECEIPT,quantitative_C3_source_sha256=sha(ranges.RECEIPT),
        same_original_C2_limit_branch=True,numerical_value_not_installed=True)
    h=[target.preserve(q,g.node('C3_Banach_limit_third_component',limit_vector=old['exact_limit_vector'].node,
        component=i,Z_order=3,implicit_ZZZ_linear_system=system.node,genuine_ordinary_third_derivative=True,
        derivative_of_range_endpoint=False)) for i,q in enumerate(h2)]
    nonlinear=quadratic_ZZZ(g,W,h)
    residual=[g.add(q.ZZZ,g.add(*(g.mul(b,v.ZZZ) for b,v in zip(row,h))),
        g.quotient(r,old['N'],'same selected exact positive integer N')) for q,row,r in zip(d,W['B'],nonlinear)]
    zeros=[g.node('proved_C3_third_zero_identity',expression=v.node,implicit_ZZZ_linear_system=system.node,
        theorem='thrice differentiate the same original signed residual and use the actual C3 implicit system',
        zero=g.zero.node,row=row) for row,v in zip(controls.ROWS,residual)]
    lower_sequence=parent['genuine_C2_Picard_sequence'];sequence=[[target.preserve(q,g.zero) for q in lower_sequence[0]]]
    for n in range(3):sequence.append(picard_C3_map(g,W,old['N'],d,sequence[-1],lower_sequence[n+1]))
    return dict(exact_C2_limit_functions=parent,exact_limit_controls_C3=h,actual_signed_targets_C3=dict(zip(controls.ROWS,d)),
        same_implicit_matrix=parent['same_implicit_matrix'],actual_implicit_ZZZ_rhs=rhs,implicit_ZZZ_linear_system=system,
        actual_third_control_residual=residual,actual_third_zero_theorem_nodes=zeros,genuine_C3_Picard_sequence=sequence,
        actual_quantitative_C3_bounds=bounds,existing_C2_value_Z_ZZ_handles_retained=True,
        original_weights_mu_and_N_Z_independent=True,C3_repaired_band_fields_and_histories_not_yet_installed=True)


class CurrentC3LimitControls:
    def __init__(self,control_field=None,owner=None,require_checked=True):
        self.lower=control_field if control_field is not None else lower.CurrentC2LimitControls(owner=owner)
        if not self.lower.acceptance_loaded:raise ValueError('Accepted actual C2 limit controls required')
        # Append on the accepted C2 graph, preserving every original control handle.
        signed=target.CurrentC3SourceTargets(target_field=self.lower.ranges.target)
        self.ranges=ranges.CurrentC3TargetRanges(target_field=signed)
        self.identity=self.ranges.identity
        if self.identity!=self.lower.identity:raise ValueError('Same current C2/C3 source family required')
        self.hashes={**self.lower.hashes,**self.ranges.hashes}
        for name in (lower.NAME,lower.RECEIPT,ranges.NAME,ranges.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        report=json.loads(gzip.decompress((HERE/ranges.NAME).read_bytes()))
        if report['source_family']!=self.identity:raise ValueError('Same C3 target range source required')
        self.prefix=[dict(n) for n in self.ranges.phase.built['graph'].nodes]
        self.bounds=quantitative_bounds(self.lower,self.ranges,report);self.functions=build(self,self.bounds)
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual C3 controls required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual C3 control source '+name)
            self.acceptance_loaded=True

    def control_functions(self):return self.functions['exact_limit_controls_C3']
    def tail(self,depth):return factored_C3_tail(self.bounds,depth)


def source_bindings():
    return dict(quadratic_ZZZ=current.ast_binding(quadratic_ZZZ),picard_C3_map=current.ast_binding(picard_C3_map),
        build=current.ast_binding(build),quantitative_bounds=current.ast_binding(quantitative_bounds),
        factored_C3_tail=current.ast_binding(factored_C3_tail),
        original_quadratic=current.ast_binding(controls.quadratic),original_quadratic_action=current.ast_binding(controls.quadratic_action),
        original_inverse_action=current.ast_binding(controls.inverse_action),
        accepted_C2_limit_source=lower.RECEIPT,accepted_C2_limit_source_sha256=sha(lower.RECEIPT))


def run(control_field=None,owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC3LimitControls(control_field=control_field,owner=owner,require_checked=False)
        report=dict(candidate_actual_C3_limit_controls_constructed=True,source_family=field.identity,
            original_C3_target_and_C2_control_graph_prefix_length=len(field.prefix),
            exact_graph_nodes=field.ranges.phase.built['graph'].nodes,actual_C3_limit_controls_and_Picard_functions=encoded(field.functions),
            actual_quantitative_C3_bounds=encoded(field.bounds),actual_C3_tail_examples=[field.tail(n) for n in (0,3,32)],
            source_bindings=source_bindings(),C3_repaired_band_fields_and_histories_installed=False,
            current_numeric_point_field_oracle_installed=False,**dict.fromkeys(GATES+OPEN,False),
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Same original C2 limit: genuine third controls and C3 Picard tail constructed',flush=True)
    return field


if __name__=='__main__':run()

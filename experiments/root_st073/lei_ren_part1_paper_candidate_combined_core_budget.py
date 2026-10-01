"""Combine finite coefficient intervals with conditional analytic radial tails.

Reports the axial Taylor center and the whole scaled radial interval. The
fixed-point/gauge identity is checked as an explicit dependency, not inferred
from a small bound or agreement of sampled values.
"""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_candidate_core_evaluation import (
    coefficient_uncertainty_bound, scaled_mixed_derivative)
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
DEFAULT_STATE=('lei_ren_part1_paper_candidate_gauge_core_'
               'lei_ren_part1_paper_candidate_pressure_axis_jets_refined_state.json')


def run(state_name=DEFAULT_STATE, budget_name='lei_ren_part1_paper_candidate_core_tail_budget.json', output_name=None):
    path=HERE/state_name
    raw=path.read_bytes()
    state=json.loads(raw)
    pressure_name=state['target']['pressure_file']
    for name,digest in state['source_hashes'].items():
        filename=({'driver':state['target'].get('driver_file','lei_ren_part1_paper_candidate_gauge_core.py'),
                   'pressure_input':pressure_name}).get(name,name)
        if hashlib.sha256((HERE/filename).read_bytes()).hexdigest()!=digest:
            raise AssertionError('State dependency changed: '+filename)
    pressure=json.loads((HERE/pressure_name).read_text())
    for name,digest in pressure['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Pressure dependency changed: '+name)
    budget=json.loads((HERE/budget_name).read_text())
    for name,digest in budget['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Analytic tail dependency changed: '+name)
    if state['accepted_schedule_sha256']!=budget['accepted_schedule_sha256']:
        raise AssertionError('Analytic and finite core sources disagree')
    if mp.mpf(state['target']['Lambda'])!=mp.mpf(budget['Lambda']):
        raise AssertionError('Analytic and finite core Lambda disagree')
    identity_name='lei_ren_part1_paper_gauge_fixed_point_identity.json'
    identity=json.loads((HERE/identity_name).read_text())
    if not identity['normalized_equations_exactly_reproduce_core']:
        raise AssertionError('Gauge to fixed-point equation identity missing')
    for name,digest in identity['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Equation identity dependency changed: '+name)
    ctx=MPIntervalContext();ctx.dps=state['target']['precision']
    with mp.workdps(ctx.dps+40):
        def restore(row):
            return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                            mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
        A=[[restore(v) for v in row] for row in state['A_rows']]
        U=[[restore(v) for v in row] for row in state['Uz_rows']]
        degree=state['completed_radial_order']
        if len(A)!=degree+1 or len(U)!=degree+1:
            raise AssertionError('State radial row count inconsistent')
        h=restore(budget['Xh_parameter'])
        radius=ctx.mpf('4.1');target=ctx.mpf('1e-12')
        rows=[]
        all_pass=True
        for component,coeffs in (('Phi',A),('Psi',U)):
            norm=restore(budget['analytic_'+component+'_Xh_norm_upper'])
            for total in range(4):
                for i in range(total+1):
                    k=total-i
                    error=coefficient_uncertainty_bound(ctx,coeffs,
                        Lambda=state['target']['Lambda'],radius_upper=radius,
                        radial_order=i,axial_order=k,component=component)
                    finite=error['uniform_midpoint_error_upper']
                    tail=norm*tail_factor(ctx,degree=degree,radial_order=i,
                        axial_order=k,radius=radius,h=h)['tail_per_Xh_norm']
                    combined=finite+tail
                    passes=endpoints(combined)[1]<=endpoints(target)[0]
                    all_pass=all_pass and passes
                    exit_value=scaled_mixed_derivative(ctx,coeffs,
                        Lambda=state['target']['Lambda'],scaled_radius=radius,
                        radial_order=i,axial_order=k,component=component)
                    rows.append(dict(component=component,scaled_radial_order=i,
                        axial_order=k,finite_coefficient_uncertainty=finite,
                        conditional_infinite_radial_tail=tail,
                        conditional_combined_midpoint_error=combined,
                        finite_exit_derivative=exit_value,
                        normalized_error_target_pass=passes))
        result=dict(state_file=path.name,state_sha256=hashlib.sha256(raw).hexdigest(),
            input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                          for name in (budget_name,identity_name)},
            accepted_schedule_sha256=state['accepted_schedule_sha256'],
            Lambda=state['target']['Lambda'],Z=state['target']['Z'],
            completed_radial_order=degree,target_radial_degree=state['target']['radial_degree'],
            scaled_radial_interval=['0','4.1'],axial_scope='Taylor center only',
            normalized_mixed_C3_error_target=target,derivative_budgets=rows,
            all_conditional_normalized_error_budgets_pass=all_pass,
            exact_gauge_recursion_to_analytic_fixed_point_identity_audited=True,
            combined_bounds_conditional_on_that_identity=False,
            exact_analytic_core_error_enclosed_relative_to_accepted_datum=True,
            whole_axis_finite_core_certified=False,physical_coordinate_errors_enclosed=False,
            original_parameter_errors_enclosed=False,core_to_collar_matching_certified=False,
            temporal_recursion=False,full_NS_residual_target=False)
        output=HERE/output_name if output_name else Path(__file__).with_suffix('.json')
        output.write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
        print('Candidate combined radial order',degree,'all C3 budgets pass',all_pass,flush=True)
        for component in ('Phi','Psi'):
            worst=max(endpoints(r['conditional_combined_midpoint_error'])[1]
                      for r in rows if r['component']==component)
            print(component,'worst mixed C3 conditional error',mp.nstr(worst,18),flush=True)
        return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--state-file',default=DEFAULT_STATE)
    parser.add_argument('--tail-budget-file',default='lei_ren_part1_paper_candidate_core_tail_budget.json')
    parser.add_argument('--output-name')
    args=parser.parse_args()
    run(args.state_file,args.tail_budget_file,args.output_name)

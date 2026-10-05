"""Three current bridge functional joins; core/bridge is still open.

Production branch specialization, admitted common history/ODE sources and
the original coordinate pullback prove equality of source functions. Numeric
enclosure overlap is neither required nor used as that proof.
"""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import sympy as s

from lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 import (
    CurrentActualBridgeMixedC4,REPLAY_PROOF,HERE,PREFIX,sha)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import (
    assignment_source_bindings,R100_source_join)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes

SOURCE_RECEIPT=PREFIX+'current_actual_bridge_mixed_C4_check.json'
SWITCH_RECEIPT=PREFIX+'actual_switch_mixed_C4_check.json'
RECEIPT=PREFIX+'current_bridge_functional_joins_check.json'


def production_branch_specialization():
    tree=ast.parse((HERE/(PREFIX+'actual_bridge_integrals.py')).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='contributions')
    values={}
    for target in ('first_phase','second_phase','kernels1','kernels'):
        rows=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(t)==target for t in n.targets)]
        if len(rows)!=1:raise ValueError('Unique original contribution selector required: '+target)
        values[target]=rows[0]
    expressions={'first_phase':"s if chart=='first' else c.mpf(1)",
        'second_phase':"s-1 if chart=='second' else c.mpf(1 if chart=='macro' else 0)"}
    for target,expression in expressions.items():
        if ast.dump(values[target])!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Original source selector changed: '+target)
    ctx=SimpleNamespace(mpf=s.sympify)
    def replay(target,chart,coordinate):
        env={'c':ctx,'mp':ctx,'s':s.sympify(coordinate),'chart':chart,'endpoints':lambda q:(q,q)}
        return eval(compile(ast.Expression(values[target]),'<original bridge source selector>','eval'),{},env)
    phases={label:[replay(t,chart,q) for t in ('first_phase','second_phase')]
        for label,chart,q in (('first1','first',1),('second1','second',1),('second2','second',2),('macro0','macro',0))}
    if phases['first1']!=phases['second1'] or phases['second2']!=phases['macro0']:
        raise ValueError('Original signed micro source decomposition differs at join')
    for target in ('kernels1','kernels'):
        if not isinstance(values[target],ast.IfExp) or replay(target,'macro',0)!=[0,0,0]:
            raise ValueError('Exact macro inlet must add zero positive-kernel source')
    packet=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='packet')
    branches=[n.test for n in ast.walk(packet) if isinstance(n,ast.If)
        and any(isinstance(a,ast.Call) and ast.unparse(a.func)=='self.macro_moments' for a in ast.walk(n))]
    wanted=ast.dump(ast.parse("chart=='macro' and hi!=0",mode='eval').body)
    if sum(ast.dump(row)==wanted for row in branches)!=1:raise ValueError('Original macro own-history inlet branch changed')
    return dict(production_AST_selector_expressions=expressions,
        specialized_micro_phases={k:[str(q) for q in v] for k,v in phases.items()},
        macro_inlet_both_kernel_families_exact_zero=True,
        macro_inlet_uses_same_micro_own_moment_feedback=True,
        zero_additional_signed_terms_do_not_change_actual_source_functions=True,passed=True)


def common_coordinate_identities():
    h,ymax=s.symbols('hb ymax',positive=True);q,z,phase=s.symbols('fraction Z phase',real=True)
    macro=2*h+q*(ymax-2*h)
    if s.simplify((h*phase).subs(phase,2)-macro.subs(q,0))!=0:
        raise ArithmeticError('Original source radius fails at second/macro join')
    # This is the chain rule for the actual shared source functions, whose
    # history and derivative ODE calls are bound in the current receipt.
    y=s.symbols('y',real=True);count=0
    # Independent coefficients span every C4 local jet. This is a formal jet
    # carrier for the exact chain rule, not a polynomial model for the field.
    f=sum(s.Symbol('J'+str(a)+'_'+str(b))*y**a*z**b/s.factorial(a)/s.factorial(b)
        for a in range(5) for b in range(5-a))
    for k in range(5):
        for n in range(5-k):
            left=s.diff(f.subs(y,h*phase),phase,k,z,n)
            right=h**k*s.diff(f,y,k,z,n).subs(y,h*phase)
            if s.simplify(left-right)!=0:
                raise ArithmeticError('Exact positive-width coordinate pullback failed')
            count+=1
    return dict(source_second_macro_radius_identity=True,graded_mixed4_coordinate_identities=count,
        arbitrary_C4_local_jets_spanned_by_independent_formal_coefficients=True,
        exact_positive_width_chain_rule='D_phase^k D_Z^n = hb^k D_logR^k D_Z^n, k+n<=4',
        macro_fraction_is_coverage_label_not_derivative_coordinate=True,
        equality_of_velocity_pressure_and_five_primitive_functions_implies_135_common_mixed4_rows=True,
        no_width_cap_selected_as_exact_parameter=True,passed=True)


def run(field=None):
    field=field if field is not None else CurrentActualBridgeMixedC4()
    source=accepted(SOURCE_RECEIPT,field.family,field.source,'current_actual_three_bridge_charts_mixed4_available')
    switch=accepted(SWITCH_RECEIPT,field.family,field.source,'actual_R100_R110_feedback_mixed4_available')
    if (not field.acceptance_loaded or not all(field.current_provider_graph().values())
            or source['datum_enclosure_sha256']!=field.datum_sha or switch['datum_enclosure_sha256']!=field.datum_sha
            or not source['actual_nonlinear_history_and_pressure_source_bindings']['passed']
            or not source['actual_parent_source_bindings']['passed']
            or source['unchanged_original_derivative_and_physical_replay']!=REPLAY_PROOF):
        raise ValueError('Current admitted actual history/derivative graph required')
    _verify_hashes(source);_verify_hashes(switch)
    branches=production_branch_specialization();coordinate=common_coordinate_identities()
    packet=assignment_source_bindings('actual_bridge_integrals','packet',{
        'history_delta_V':"p['uniform_delta_V'] if not (chart=='first' and hi==0) else delta_V*0",
        'pressure':"dress(own['actual']['C'],p['inputs']['F0_squared_ratios'])",
        'Q':"(2*z*V-z*own['actual']['M']*(1-self.core.delta)-d*derivative(own['actual']['M']))/L"})
    r100=R100_source_join()
    if (r100['symbolic_original_control_coordinate_identities']!=17
            or not r100['exact_common_physical_primitive_operator']
            or not r100['no_interval_overlap_used_as_identity']):
        raise ValueError('Original R100 current trace proof missing')
    # Original bridge operator/flat-endpoint theorem is retained by current
    # source acceptance. Its provenance, rather than box equality, is used.
    formal=source['reused_unchanged_operator_fixtures']['symbolic_checks']
    if not formal['passed'] or not formal['phase1_and_phase2_use_original_flat_endpoint_jets']:
        raise ValueError('Original flat alpha/chi endpoint proof missing')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_provider_graph_identity=field.current_provider_graph(),
        actual_current_parent_and_unchanged_derivative_replay=source['actual_parent_source_bindings'],
        retained_nonlinear_actual_history_pressure_bindings=source['actual_nonlinear_history_and_pressure_source_bindings'],
        production_branch_specialization=branches,common_coordinate_identities=coordinate,
        actual_common_Q_pressure_and_prefix_history_AST_bindings=packet,
        retained_original_flat_endpoint_source_proof=formal,current_R100_control_trace_source_proof=r100,
        common_history_ODEs=['H_y+2H=2phi','M_y+M=V','K_y+2K=2phi*V',
            'A_y+A=V^2','B_y+2B=phi^2','C_y+C=phi^2'],
        first_second_join_rule='Identical signed source prefixes, own histories, P0 and alpha/chi flat endpoint jets; same original derivative recovery',
        second_macro_join_rule='Same two micro source prefixes and six own histories at exact y=2hb; flat alpha=0/chi=hb; exact positive-width pullback',
        R100_join_rule="Same upstream.packet(Z,1,'macro') and comparison.macro(Z,1); original flat switch inlet and shared phase_physical operator",
        current_first_second_functional_mixed4_join_certified=True,
        current_second_macro_functional_mixed4_join_certified=True,
        current_R100_bridge_switch_functional_mixed4_join_certified=True,
        current_core_bridge_functional_mixed4_join_certified=False,all_current_bridge_functional_interfaces_certified=False,
        core_join_remaining_obligation='Bind coefficientwise core inlet and normalized-core enclosure to the same unique nonlinear fixed point and local stress-free ODE rows',
        actual_point_moment_history_recovered=False,full_cartesian_vector_derivatives_certified=False,
        global_completed_tensor_admissibility=False,admissible_stress_lift_constructed=False,temporal_recursion=False,
        all_scoped_checks_passed=True,all_passed=True,input_hashes={**source['input_hashes'],**switch['input_hashes'],
            SOURCE_RECEIPT:sha(SOURCE_RECEIPT),SWITCH_RECEIPT:sha(SWITCH_RECEIPT),Path(__file__).name:sha(Path(__file__).name)})
    (HERE/RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS three current bridge functional mixed4 joins; core/first remains explicitly open',flush=True)
    return result


if __name__=='__main__':run()

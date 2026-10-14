"""Original source identities and independent C6/C7 transport calculus."""
import gzip
import json
import math
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions as current
import lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source_check as inlet_check
from lei_ren_part1_paper_interval_taylor import IntervalTaylor


class KernelAlgebra(inlet_check.KernelAlgebra):
    def at(self,index):
        if index not in self.cache and self.g.nodes[index]['operation']=='bound_variable':
            self.cache[index]=sy.Symbol(self.g.nodes[index]['name'])
        return super().at(index)


def free_variables(g,ref,cache):
    i=ref.node if isinstance(ref,current.source.FunctionRef) else ref
    if i in cache:return cache[i]
    row=g.nodes[i];op=row['operation'];visit=lambda j:free_variables(g,j,cache)
    if op=='bound_variable':out={row['name']}
    elif op in ('exact_rational','mathematical_positive_infinity','mathematical_pi'):out=set()
    elif op in ('original_source_parameter','current_original_source_parameter'):
        assert row['Z_independent'];out=set()
    elif op=='current_original_unique_analytic_core_component':out=visit(row['rho'])|visit(row['Z_variable'])
    elif op=='current_original_leading_function_recipe':
        assert row['quantity']=='P0' and row['coordinate_independent'];out=visit(row['Z_variable'])
    elif op=='current_original_axis_pressure_function_derivative':out=visit(row['Z_variable'])
    elif op in ('sum','product'):out=set().union(*(visit(j) for j in row['arguments']))
    elif op=='positive_quotient':out=visit(row['numerator'])|visit(row['denominator'])
    elif op in ('negative','analytic_unary','exact_original_source_flat_sigma'):out=visit(row['argument'])
    elif op=='definite_integral':out=(visit(row['integrand'])-{row['variable']})|visit(row['lower'])|visit(row['upper'])
    elif op=='function_substitution':
        original=visit(row['expression']);name=g.nodes[row['variable']]['name']
        out=(original-{name})|visit(row['value']) if name in original else original
    else:raise AssertionError('Unproved free-variable rule '+op)
    cache[i]=out;return out


def actual_source(owner):
    g=owner.g;r=current.rows;cache={};rho=g.nodes[owner.core.rho.node]['name']
    assert owner.inlet.acceptance_loaded
    for name,jet in owner.core7.items():
        assert len(r(jet))==8 and r(jet)[:7]==r(owner.inlet.core_functions[name])
    assert r(owner.core7['G'])[1:]==r(owner.core7['gradient'])[:7]
    assert all(rho not in free_variables(g,q,cache) for q in r(owner.core7['F0']))
    assert len(r(owner.P0))==8 and r(owner.P0)[:7]==r(owner.inlet.P0)
    assert g.nodes[owner.hb.node]['Z_independent']
    assert owner.bindings==current.bridge.bridge_control_bindings()
    counts={}
    for label,chart in (('micro',owner.micro),('macro',owner.macro)):
        assert all(len(r(q))==7 for q in chart['fields'].values())
        assert all(len(r(q))==6 for q in chart['generic'].values())
        assert all(len(r(q))==7 for q in chart['moments'].values())
        assert all(len(r(q))==7 for q in chart['radial'].values())
        for collection in ('actual','moments','fields','radial','generic'):
            for q in chart[collection].values():
                for ref in r(q):assert rho not in free_variables(g,ref,cache),(label,collection)
        counts[label]=sum(len(r(q)) for collection in ('fields','moments','radial','generic') for q in chart[collection].values())
    for chart in (owner.comparison_micro,owner.comparison_macro):
        assert all(len(r(q))==8 for q in chart['moments'].values())
        assert all(len(r(q))==7 for q in chart['directions'].values())
    assert owner.comparison_macro['phi'] is owner.comparison_micro_end['phi']
    assert owner.comparison_macro['U'] is owner.comparison_micro_end['U']
    # The original frozen-comparison B equilibrium is half the forcing.
    assert current.RATES==dict(H=2,M=1,K=2,A=1,B=2,C=1)
    z=owner.inlet.hierarchy.source.cell(('0.375','0.375'))
    assert owner.project_P0_seventh(owner.P0.rows[7],z)==owner.P0_rows(z)[7]*5040
    return dict(passed=True,actual_source_output_rows=counts,
        preserved_core_C6_before_actual_C7_comparison_source=True,
        G_higher_rows_follow_original_anchored_FTC=True,
        F0_is_Z_only_and_all_phase_bodies_have_no_free_core_rho=True,
        original_P0_C6_preserved_and_executable_Z7_projection=True,
        same_original_both_micro_controls_and_comparison_direction_source=True,
        known_comparison_C7_actual_background_C6_generic_C5=True,
        modified_finite_N_source_and_quantitative_C6_not_claimed=True)


def ordinary_jet(g,Z,body,order,phase=None):
    variables={'Z':g.symbol('Z')}
    if phase is not None:variables[str(phase)]=g.symbol(str(phase))
    return current.FunctionJet(tuple(current.source.expression(g,sy.diff(body,Z,k),variables) for k in range(order+1)))


def unary_fixture():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');body=2+Z**2+Z**7/31
    a=current.BridgeAlgebra(g,7);q=ordinary_jet(g,Z,body,7)
    result=dict(exp=a.exp(q),log=a.log(q,'fixture positive body'))
    alg=KernelAlgebra(g,{g.symbol('Z').node:sy.Rational(3,8)});count=0
    for name,jet in result.items():
        expected=getattr(sy,name)(body)
        for k,ref in enumerate(current.rows(jet)):
            assert sy.simplify(alg.at(ref.node)-sy.diff(expected,Z,k).subs(Z,sy.Rational(3,8)))==0,(name,k)
            count+=1
    return dict(passed=True,independent_ordinary_exp_log_rows=count,seventh_input_not_padded=True)


def direction_fixture():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');a=current.BridgeAlgebra(g,7)
    bodies=dict(phi=2+Z**2/10+Z**7/23,U=3+Z+Z**7/29,H=1+Z**3/37+Z**7/41,
        M=2+Z**2/17+Z**7/43,K=2+Z**4/31+Z**7/47,A=4+Z**2/13+Z**7/53,
        B=sy.Rational(1,2)+Z**2/19+Z**7/59,C=1+Z**3/23+Z**7/61,
        F0=2+Z**2/11+Z**7/67,P0=-3-Z**2/29-Z**7/71)
    jets={k:ordinary_jet(g,Z,v,7) for k,v in bodies.items()};delta=sy.Rational(1,5)
    new=current.comparison_directions(g,jets['phi'],jets['U'],{k:jets[k] for k in current.RATES},jets['F0'],jets['P0'],g.constant('1/5'))
    physical_swirl=current.BridgeAlgebra(g,6).mul(new['drive_swirl'],a.trim(a.mul(jets['F0'],jets['F0'])))
    class Exact:
        @staticmethod
        def mpf(value):return sy.Rational(value)
    c=Exact();point=sy.Rational(3,8)
    tj=lambda expr:IntervalTaylor(c,[sy.diff(expr,Z,k).subs(Z,point)/math.factorial(k) for k in range(8)])
    F=bodies['F0'];ratios=lambda f:[sy.diff(f,Z,k).subs(Z,point)/f.subs(Z,point) for k in range(8)]
    old=current.profiles.direction(c,point,delta,tj(bodies['phi']),tj(bodies['U']),
        current.bridge.named_moments({k:tj(bodies[k]) for k in current.RATES}),tj(bodies['P0']),ratios(F),ratios(F**2))
    alg=KernelAlgebra(g,{g.symbol('Z').node:point});count=0
    for name,jet in new.items():
        target=physical_swirl if name=='drive_swirl' else jet
        for k,ref in enumerate(current.rows(target)):
            factor=F.subs(Z,point)**2 if name=='drive_swirl' else 1
            assert sy.cancel(alg.at(ref.node)-old[name][k]*math.factorial(k)*factor)==0,(name,k)
            count+=1
    return dict(passed=True,actual_original_direction_method_ordinary_rows=count,
        genuine_seventh_mean_energy_P0_inputs=True,
        physical_swirl_F0_squared_dressing_and_external_scale_once=True,
        original_linear_quadratic_pressure_and_meridional_terms_retained=True)


def transport_fixture():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');S=sy.Symbol('fixture_transport_phase')
    a=current.BridgeAlgebra(g,6);phase=g.symbol(str(S));point=sy.Rational(3,8);jac=sy.Rational(2,7)
    phi=2+Z**2/17+Z**6/19;U0=3+Z+Z**6/23;F=2+Z**3/29+Z**6/31
    raw=dict(D_over_R=sy.Integer(0),drive_hydro=1+Z**6/7,drive_pressure=2+Z**3/11+Z**6/19,
        drive_swirl=3+Z**2/13+Z**6/23)
    jets={k:ordinary_jet(g,Z,v,6) for k,v in raw.items()};Pstar,R,chi=sy.Rational(3,2),sy.Rational(2,9),sy.Rational(1,3)
    got=current.actual_fields(a,jets,dict(phi=ordinary_jet(g,Z,phi,6),U=ordinary_jet(g,Z,U0,6)),
        phase,g.constant(str(jac)),g.constant(str(chi)),g.constant(str(R)),ordinary_jet(g,Z,F,6),
        g.constant(str(Pstar)),ordinary_jet(g,Z,phi,6),'fixture_transport')
    drive=R*raw['drive_hydro']+R*Pstar**2*raw['drive_pressure']+R**2*F**2*raw['drive_swirl']
    exact=dict(phi=phi,U=U0-S*jac*chi*drive,phi_y=sy.Integer(0),U_y=-chi*drive)
    alg=KernelAlgebra(g,{g.symbol('Z').node:point,phase.node:sy.Rational(1,2)});count=0
    for name,body in exact.items():
        for k,ref in enumerate(current.rows(got[name])):
            assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,k).subs({Z:point,S:sy.Rational(1,2)}))==0,(name,k)
            count+=1
    return dict(passed=True,independent_signed_transport_rows=count,
        nonunit_jacobian_chi_radius_Pstar_and_F0=True,
        full_axial_minus_sign_and_pressure_swirl_physical_scales_once=True)


def history_fixture():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');S=sy.Symbol('fixture_six_history_phase')
    phase=g.symbol(str(S));a=current.BridgeAlgebra(g,6);jac=sy.Rational(1,9);point=sy.Rational(3,8)
    phi=(2+Z**3/13+Z**6/19)*(1+S);U=(3+Z**2/17+Z**6/23)*(1-S/7)
    initial={k:1+(i+1)*Z+Z**6/(31+i) for i,k in enumerate(current.RATES)}
    got=current.histories(a,ordinary_jet(g,Z,phi,6,S),ordinary_jet(g,Z,U,6,S),
        {k:ordinary_jet(g,Z,v,6) for k,v in initial.items()},phase,g.constant(str(jac)),'fixture_six')
    densities=dict(H=2*phi,M=U,K=2*phi*U,A=U**2,B=phi**2,C=phi**2);t=sy.Symbol('fixture_physical_t')
    alg=KernelAlgebra(g,{g.symbol('Z').node:point});count=0
    for name,jet in got.items():
        rate=current.RATES[name]
        for k,ref in enumerate(current.rows(jet)):
            actual=alg.at(ref.node);density=sy.diff(densities[name],Z,k).subs(Z,point)
            incoming=sy.diff(initial[name],Z,k).subs(Z,point)
            assert sy.simplify(sy.diff(actual,S)-jac*(density-rate*actual))==0,(name,k,'ODE')
            assert sy.simplify(actual.subs(S,0)-incoming)==0,(name,k,'memory')
            count+=1
    return dict(passed=True,independent_six_history_ODE_and_inlet_rows=count,
        all_nonzero_incoming_histories_and_sixth_derivatives_retained=True,
        rates_H2_M1_K2_A1_B2_C1_and_B_forcing_not_halved=True,
        nonunit_physical_jacobian_applied_once=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['prescribed_bridge_C6_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks=dict(actual_source=actual_source(owner),unary_calculus=unary_fixture(),original_direction=direction_fixture())
    print('TRUE_COMPARISON_C7_DIRECTION_C6_CHECKED',flush=True)
    checks['signed_transport']=transport_fixture();checks['own_six_histories']=history_fixture()
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_PRESCRIBED_BRIDGE_C6_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses the retained checked core')

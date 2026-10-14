"""Scoped independent implicit calculus, signed density and memory checks."""
import json
import math
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions as current
import lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check as prior


class Algebra(prior.Algebra):
    def at(self,i):
        if i in self.cache:return self.cache[i]
        row=self.g.nodes[i];op=row['operation']
        if op=='exact_real_comparison':
            left,right=self.at(row['left']),self.at(row['right'])
            value={'<=':lambda:left<=right,'>=':lambda:left>=right}[row['operator']]()
        elif op=='logical_or':value=any(self.at(j) for j in row['arguments'])
        elif op=='original_lazy_flat_branch':
            value=self.at(row['flat_value']) if self.at(row['flat_predicate']) else self.at(row['active_body'])
        elif op=='substitute_original_inverse_angle':
            value=self.at(row['body']).subs(sy.Symbol(row['angle_variable']),self.at(row['inverse_angle']))
        elif op=='analytic_unary' and row['name'] in ('positive_sqrt','expm1'):
            x=self.at(row['argument']);value=sy.sqrt(x) if row['name']=='positive_sqrt' else sy.exp(x)-1
        elif op=='exact_original_source_flat_sigma':
            x=self.at(row['argument'])
            if x<=0:value=sy.Integer(0)
            elif x>=1:value=sy.Integer(1)
            else:
                e=sy.exp(1/(1-x)**2-1/x**2);value=e/(1+e)
        elif op=='mathematical_pi':value=sy.pi
        else:return super().at(i)
        self.cache[i]=value;return value


def free_variables(g,ref,cache):
    i=ref.node if isinstance(ref,current.source.FunctionRef) else ref
    if i in cache:return cache[i]
    row=g.nodes[i];op=row['operation'];visit=lambda j:free_variables(g,j,cache)
    if op=='exact_positive_integer_expression':out=set()
    elif op=='exact_real_comparison':out=visit(row['left'])|visit(row['right'])
    elif op=='logical_or':out=set().union(*(visit(j) for j in row['arguments']))
    elif op=='original_lazy_flat_branch':
        out=visit(row['active_body'])|visit(row['flat_value'])|visit(row['flat_predicate'])
    elif op=='original_monotone_phase_inverse':
        out=(visit(row['phase_function'])-{row['angle_variable']})|visit(row['target_phase'])
    elif op=='substitute_original_inverse_angle':
        original=visit(row['body']);out=(original-{row['angle_variable']})|visit(row['inverse_angle'])
    else:
        # Extend the existing lexical rules without losing scope at binders.
        if op=='bound_variable':out={row['name']}
        elif op in ('exact_rational','mathematical_pi','mathematical_positive_infinity'):out=set()
        elif op in ('original_source_parameter','current_original_source_parameter'):
            assert row['Z_independent'];out=set()
        elif op=='current_original_unique_analytic_core_component':out=visit(row['rho'])|visit(row['Z_variable'])
        elif op in ('current_original_leading_function_recipe','current_original_axis_pressure_function_derivative'):
            assert row.get('coordinate_independent',True);out=visit(row['Z_variable'])
        elif op in ('sum','product'):out=set().union(*(visit(j) for j in row['arguments']))
        elif op=='positive_quotient':out=visit(row['numerator'])|visit(row['denominator'])
        elif op in ('negative','analytic_unary','exact_original_source_flat_sigma'):out=visit(row['argument'])
        elif op=='definite_integral':out=(visit(row['integrand'])-{row['variable']})|visit(row['lower'])|visit(row['upper'])
        elif op=='function_substitution':
            original=visit(row['expression']);name=g.nodes[row['variable']]['name']
            out=(original-{name})|visit(row['value']) if name in original else original
        else:raise AssertionError('Unknown lexical loop source '+op)
    cache[i]=out;return out


def actual_source(owner):
    g=owner.g;r=current.rows;cache={};count=0;before=owner.zero_inlet
    assert owner.background.acceptance_loaded and tuple(owner.charts)==current.CHARTS
    assert g.nodes[owner.N.node]['one_common_N_for_all_charts_and_control_denominators']
    for ref,name in ((owner.N,None),(owner.eta,'eta'),(owner.dstar,'d_star'),(owner.sc,'sc')):
        assert ref.node<len(owner.background.prefix)
        assert g.nodes[ref.node]==owner.background.prefix[ref.node]
        assert owner.unique('exact_positive_integer_expression' if name is None else 'current_original_source_parameter',name)==ref
        if name is not None:
            row=g.nodes[ref.node]
            assert row['source_family']==owner.family and row['Z_independent']
            assert row['current_source_receipt_sha256']==current.sha(row['current_source_receipt'])
            if name in ('eta','d_star'):
                proof=row['canonical_parameter_proof']
                assert proof['source_family']==owner.family and proof['exact_original_parameter_tuples_equal']
    assert owner.bindings['original_loop']==current.parent.active.current.ast_binding(current.previous.connection.original_loop)
    assert owner.bindings['original_density_graph']==current.parent.active.current.ast_binding(current.original.build)
    assert all(q==g.zero for jet in before.values() for q in r(jet))
    for name,chart in owner.charts.items():
        native=g.nodes[chart['coordinate'].node]['name']
        assert chart['source_chart'] is owner.source_charts[name]
        assert chart['incoming'] is before
        assert chart['roots']['p2'] is chart['source_chart']['generic']['p2']
        assert all(len(r(jet))==7 for jet in chart['source_chart']['fields'].values())
        assert all(len(r(jet))==6 for jet in chart['source_chart']['generic'].values())
        assert free_variables(g,chart['actual_phase'],cache)<={native}
        assert 'Z' not in free_variables(g,chart['physical_measure'],cache)
        for group in (chart['roots'],chart['history'],chart['modified_fields'],
                      chart['finite_N']['five_signed_density_increments']):
            for jet in group.values():
                assert len(r(jet))==6
                for ref in r(jet):assert free_variables(g,ref,cache)<={'Z',native},(name,ref.node)
                count+=6
        for key in ('q','inverse','A','B_over_Pstar','A_phi','B_over_Pstar_phi'):
            assert len(r(chart['loop'][key]))==6
            for ref in r(chart['loop'][key]):assert free_variables(g,ref,cache)<={'Z',native}
        assert r(chart['modified_fields']['P0'])==r(owner.P0)[:6]
        assert chart['memory']['p']==g.one
        if name=='first_micro':
            assert chart['domain'][0]==g.mul(owner.sc,g.constant('1/2'))
            assert chart['loop']['original_collar_upper']==g.mul(g.constant('3/4'),owner.sc)
        else:assert chart['loop']['original_collar_upper'] is None
        for key,jet in chart['contributions'].items():
            for order,ref in enumerate(r(jet)):
                row=g.nodes[ref.node]
                assert row['operation']=='definite_integral' and row['ordinary_slow_Z_derivative_order']==order
                assert row['lower']==chart['domain'][0].node and row['upper']==chart['coordinate'].node
                assert row['physical_Jacobian_applied_once'] and row['predecessor_memory_never_reset']
                assert row['original_rate']==str(current.RATES[key])
        before=owner.endpoints[name]
    assert owner.Rm_correction is before
    return dict(passed=True,true_prescribed_charts=list(current.CHARTS),actual_source_history_density_rows=count,
        genuine_C6_inputs_and_C5_generic_roots=True,full_signed_p2_preserved_by_identity=True,
        same_original_N_eta_dstar_and_global_phase=True,all_ten_correction_predecessors_preserved=True,
        original_microscopic_collar_and_pressure_memory_retained=True,
        no_free_core_rho_angle_old_phase_or_integration_variable=True,
        all17_numeric_global_frequency_and_recursion_not_claimed=True)


def close(actual,expected,label):
    error=sy.N(actual-expected,85)
    assert abs(error)<sy.Float('1e-65',85),(label,error)


def jet(g,Z,body,order=5):
    return current.FunctionJet(tuple(current.source.expression(g,sy.diff(body,Z,k),{'Z':g.symbol('Z')})
        for k in range(order+1)))


def elementary_calculus():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');point=sy.Rational(1,7);a=current.LoopAlgebra(g,5)
    root=2+Z+Z**5/19;arg=sy.Rational(1,2)+Z/13+Z**2/17+Z**5/23
    sigma=sy.exp(1/(1-arg)**2-1/arg**2);sigma=sigma/(1+sigma)
    got=dict(sqrt=a.sqrt(jet(g,Z,root)),sigma=a.sigma(jet(g,Z,arg)))
    alg=Algebra(g,{g.symbol('Z').node:point});count=0
    for name,body in (('sqrt',sy.sqrt(root)),('sigma',sigma)):
        for k,ref in enumerate(current.rows(got[name])):
            close(alg.at(ref.node),sy.diff(body,Z,k).subs(Z,point),(name,k));count+=1
    for value in (sy.Rational(-1,3),sy.Rational(4,3)):
        flat=a.sigma(a.fixed(str(value)));alg=Algebra(g,{})
        assert alg.at(flat.value.node)==int(bool(value>1))
        assert all(alg.at(ref.node)==0 for ref in current.rows(flat)[1:])
    return dict(passed=True,independent_sqrt_and_interior_sigma_C5_rows=count,
        both_exact_flat_sigma_extensions_checked=True,numerical_calculus_fixture_only=True)


def implicit_calculus():
    """Closed explicit inverse exercises all total-degree cross coefficients."""
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');angle=g.symbol('fixture_angle');psi=sy.Symbol('fixture_angle')
    point=sy.Rational(3,8);phase=sy.Rational(7);body=(1+Z**2)*sy.exp(psi)+Z**3
    inverse=g.symbol('fixture_known_inverse');known=sy.log((phase-Z**3)/(1+Z**2));psi0=known.subs(Z,point)
    partials={(k,j):current.source.expression(g,sy.diff(body,Z,k,psi,j)/sy.factorial(k)/sy.factorial(j),
        {'Z':g.symbol('Z'),'fixture_angle':angle}) for k in range(6) for j in range(6-k)}
    got,delta,_=current.implicit_inverse(g,partials,inverse,angle)
    T1=Z**2*psi+Z*psi**3/3+psi**6/7
    T1partials={(k,j):current.source.expression(g,sy.diff(T1,Z,k,psi,j)/sy.factorial(k)/sy.factorial(j),
        {'Z':g.symbol('Z'),'fixture_angle':angle}) for k in range(6) for j in range(6-k)}
    composed=current.compose_inverse(g,T1partials,inverse,angle,delta)
    # Keep the free angle symbolic until at-inverse substitution.
    alg=Algebra(g,{g.symbol('Z').node:point,inverse.node:psi0})
    for k,ref in enumerate(current.rows(got)):close(alg.at(ref.node),sy.diff(known,Z,k).subs(Z,point),('inverse',k))
    for k,ref in enumerate(current.rows(composed)):
        close(alg.at(ref.node),sy.diff(T1.subs(psi,known),Z,k).subs(Z,point),('composition',k))
    # Check the defining equation independently on every solved coefficient.
    at={(k,j):alg.at(current.at_inverse(g,v,inverse,angle).node) for (k,j),v in partials.items()}
    powers=current.series_powers(g,delta,5)
    for n in range(1,6):
        residual=sum(at[k,j]*alg.at(powers[j][n-k].node)
            for k in range(n+1) for j in range(n-k+1))
        close(residual,sy.Integer(0),('implicit residual',n))
    return dict(passed=True,independent_explicit_inverse_C5_rows=6,
        independent_nonseparable_composition_C5_rows=6,implicit_equation_coefficients_checked=5,
        normalized_factorials_and_all_cross_terms_checked=True,numerical_calculus_fixture_only=True)


def finite_N_calculus():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');point=sy.Rational(1,7);N=sy.Integer(37)
    E=2+Z**2/13+Z**5/17;V=1+Z/19+Z**5/23;A=Z/29+Z**2/31+Z**5/37;B=-1+Z**3/41+Z**5/43
    roots={k:jet(g,Z,v) for k,v in dict(E=E,V=V).items()}
    loop={k:jet(g,Z,v) for k,v in dict(A=A,B_over_Pstar=B).items()}
    got=current.finite_N_functions(g,roots,loop,g.constant(37));r=E*(sy.exp(A/N)-1);s=B/N
    expected=dict(E_N=E*sy.exp(A/N),V_N=V+s,delta_E=r,delta_V=s)
    increments=dict(m=s,h=r,k=V*r+E*s+r*s,e=2*V*s+s*s-E*r-r*r/2,p=E*r+r*r/2)
    alg=Algebra(g,{g.symbol('Z').node:point});count=0
    for group,want in ((got,expected),(got['five_signed_density_increments'],increments)):
        for key,body in want.items():
            for k,ref in enumerate(current.rows(group[key])):
                close(alg.at(ref.node),sy.diff(body,Z,k).subs(Z,point),(key,k));count+=1
    return dict(passed=True,independent_nonzero_signed_velocity_and_density_C5_rows=count,
        original_nonunit_N_and_all_quadratic_cross_terms_checked=True,numerical_calculus_fixture_only=True)


def transport_calculus():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');s=g.symbol('fixture_phase');S=sy.Symbol('fixture_phase')
    variables={'Z':g.symbol('Z'),'fixture_phase':s};a=current.LoopAlgebra(g,5);lower=sy.Rational(1,3);jac=sy.Rational(3,2)
    density={name:(1+Z+Z**5/(11+i))*(1+S+S**2/(13+i)) for i,name in enumerate(current.RATES)}
    inlet={name:2+Z**2/(17+i)+Z**5/(19+i) for i,name in enumerate(current.RATES)}
    make=lambda body:current.FunctionJet(tuple(current.source.expression(g,sy.diff(body,Z,k),variables) for k in range(6)))
    got,_,_=current.transport(g,{k:make(v) for k,v in density.items()},{k:make(v) for k,v in inlet.items()},
        s,g.constant(str(lower)),g.constant(str(jac)),'fixture')
    alg=Algebra(g,{g.symbol('Z').node:Z,s.node:S});count=0
    for name,rate in current.RATES.items():
        for k,ref in enumerate(current.rows(got[name])):
            value=alg.at(ref.node)
            assert sy.simplify(sy.diff(value,S)-jac*(sy.diff(density[name],Z,k)-sy.Rational(rate.numerator,rate.denominator)*value))==0,(name,k)
            assert sy.simplify(value.subs(S,lower)-sy.diff(inlet[name],Z,k))==0,(name,k,'inlet')
            count+=1
    return dict(passed=True,independent_nonzero_Duhamel_ODE_and_inlet_C5_rows=count,
        all_five_original_rates_and_single_physical_measure_checked=True,
        nonzero_rate_zero_pressure_memory_checked=True,mathematical_fixture_only=True)


def lazy_flat_calculus():
    g=current.source.FunctionTransportGraph();a=current.LoopAlgebra(g,5);coordinate=g.symbol('fixture_coordinate')
    roots={key:a.fixed(value) for key,value in dict(a=2,b='1/5',t0='-1/10',E=3,V=1,p1=4,p2='-1/3',Delta=1).items()}
    count=0
    for label,delta,position in (('first_micro',0,'1/2'),('first_switch',1,0)):
        roots['Delta']=a.fixed(delta)
        loop=current.loop_functions(g,roots,g.constant('2/5'),g.constant('1/2'),g.constant('1/3'),label,coordinate,g.one)
        finite=current.finite_N_functions(g,roots,loop,g.constant(37));alg=Algebra(g,{coordinate.node:sy.Rational(position)})
        for jet in (loop['q'],loop['A'],loop['B_over_Pstar'],loop['A_phi'],loop['B_over_Pstar_phi'],
                    finite['delta_E'],finite['delta_V'],*finite['five_signed_density_increments'].values()):
            for ref in current.rows(jet):assert alg.at(ref.node)==0;count+=1
    return dict(passed=True,lazy_original_collar_and_shear_flat_rows=count,
        active_inverse_and_singular_sqrt_not_evaluated_on_flat_branch=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['prescribed_loop_C5_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks=dict(actual_source=actual_source(owner))
    print('TRUE_TEN_CHART_C5_SOURCE_AND_MEMORY_CHECKED',flush=True)
    for name,fn in (('elementary',elementary_calculus),('implicit_inverse',implicit_calculus),
                    ('finite_N',finite_N_calculus),('transport',transport_calculus),('lazy_flat',lazy_flat_calculus)):
        checks[name]=fn();print('TRUE_LOOP_C5_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_TEN_CHART_PRESCRIBED_LOOP_AND_TRANSPORT_C5_CHECKED',flush=True)
    return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses the actual prescribed source owner')

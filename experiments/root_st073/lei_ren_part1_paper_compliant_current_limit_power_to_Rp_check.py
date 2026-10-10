"""Independent actual graph calculus and accepted source prefix bindings."""
import gzip
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_limit_power_to_Rp as current


def graph_calculus():
    g=current.source.FunctionTransportGraph();pair=lambda k:current.source.C1Function(g.symbol(k),g.symbol(k+'_Z'))
    old={k:pair(k) for k in ('m','h','k','e','p')};A=pair('A');mu=g.symbol('mu');s=g.symbol('s')
    flow=current.power_flow(g,A,old,mu,s);frame=current.pulse_input(g,flow['fields']['E'],flow['own'],pair('P0'))
    symbols={}
    def read(ref):
        row=g.nodes[ref if type(ref) is int else ref.node];op=row['operation']
        if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
        if op=='bound_variable':return symbols.setdefault(row['name'],sy.Symbol(row['name'],real=True))
        if op=='sum':return sum(map(read,row['arguments']))
        if op=='product':return sy.prod(map(read,row['arguments']))
        if op=='negative':return -read(row['argument'])
        if op=='positive_quotient':return read(row['numerator'])/read(row['denominator'])
        if op=='analytic_unary':
            arg=read(row['argument'])
            if row['name']=='exp':return sy.exp(arg)
            if row['name']=='exprel':return sy.Integer(1) if arg==0 else (sy.exp(arg)-1)/arg
        raise AssertionError(row)
    sym=lambda k:sy.Symbol(k,real=True)
    a,mu,s=map(sym,('A','mu','s'));m,h,k,e,p=map(sym,('m','h','k','e','p'))
    ep=a*sy.exp(-(sy.Rational(1,2)+mu)*s)
    expected=dict(m=m*sy.exp(-s),h=h*sy.exp(-sy.Rational(3,2)*s)+a*(sy.exp(-(sy.Rational(1,2)+mu)*s)-sy.exp(-sy.Rational(3,2)*s))/(1-mu),
        k=k*sy.exp(-sy.Rational(3,2)*s),e=sy.exp(-s)*(e-a*a*(1-sy.exp(-2*mu*s))/(4*mu)),
        p=p+a*a*(1-sy.exp(-(1+2*mu)*s))/(2*(1+2*mu)))
    rows=[(flow['fields']['E'],ep)]+[(flow['own'][key],expected[key]) for key in expected]
    variables=(a,m,h,k,e,p);zrows=tuple(sym(str(q)+'_Z') for q in variables)
    for actual,expr in rows:
        assert sy.simplify(sy.expand_power_exp(read(actual.value)-expr))==0
        derivative=sum(sy.diff(expr,v)*vZ for v,vZ in zip(variables,zrows))
        assert sy.simplify(sy.expand_power_exp(read(actual.Z)-derivative))==0
    for key,pair in flow['own'].items():
        assert sy.simplify(sy.expand_power_exp(read(flow['own_y'][key].value)-sy.diff(expected[key],s)))==0,key
    frames=dict(u=ep,m1=expected['m']/ep,m2=expected['k']/ep**2,
        X=expected['h']/ep,energy=expected['e']/ep**2,Mp=expected['p'],P0=sym('P0'),pressure=sym('P0')+expected['p'])
    for key,expr in frames.items():
        assert sy.simplify(sy.expand_power_exp(read(frame[key].value)-expr))==0,key
        derivative=sum(sy.diff(expr,v)*vZ for v,vZ in zip(variables,zrows))+sy.diff(expr,sym('P0'))*sym('P0_Z')
        assert sy.simplify(sy.expand_power_exp(read(frame[key].Z)-derivative))==0,key+'_Z'
    return dict(actual_graph_all_five_semigroup_values_and_Z_checked=True,
        actual_radial_rows_match_independent_derivatives=True,
        actual_pulse_frame_values_and_amplitude_Z_quotient_rules_checked=True)


def run():
    began=time.monotonic();accepted,native,hashes=current.load_inputs()
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert report[current.GATE] and report['source_family']==accepted['source_family']
    for name,digest in report['input_hashes'].items():assert current.sha(name)==digest,name
    built=current.build(accepted);g=built['graph']
    assert g.nodes==report['exact_graph_nodes']
    assert current.current.current.current.encode_graph(built)==report['exact_power_to_Rp']
    assert g.nodes[:built['source_graph_prefix_length']]==accepted['exact_graph_nodes']
    assert built['source_graph_prefix_length']==6601
    assert built['P0'].value.node==accepted['exact_physical_recovery']['P0'][0]
    assert set(built['terminal_Rp_histories'])=={'m','h','k','e','p'}
    assert set(built['cartesian_velocity_pressure'])=={'u','v','w','p'}
    assert built['cartesian_velocity_pressure']['w']==g.zero
    assert built['recovery']['Q']!=g.zero
    propagation=g.nodes[built['relative_zero_propagation'].node]
    assert propagation['inlet_certificate_nodes']==accepted['exact_physical_recovery']['relative_terminal_zero_certificates']
    assert propagation['does_not_set_absolute_histories_to_zero'] and propagation['pressure_constant_not_reset']
    assert built['pulse_input_frame']['P0']==built['P0']
    assert built['pulse_input_frame']['common_m_k_already_divided_by_S']
    assert not built['native_bridge_requirements']['actual_Rp_native_frame_source_function_identified']
    assert not built['native_bridge_requirements']['exact_native_constructor_consumes_current_limit_frame']
    assert report['symbolic_checks']==current.symbolic_checks()
    calculus=graph_calculus()
    for key in ('actual_Rp_native_frame_source_function_identified','actual_native_O4_constructor_consumes_current_limit_frame',
        'actual_preheat_pressure_frame_identified','physical_original_exterior_five_targets_closed',
        'actual_numeric_point_source_oracle_installed','actual_numeric_controls_evaluated','actual_global_frequency_admitted',
        'full_recovered_velocity_pressure_and_heat_joins_admitted','higher_Z_velocity_jets_and_stress_cone_admitted',
        'actual_temporal_scale_recursion_installed',*current.current.current.outer.OPEN):assert report[key] is False,key
    hashes[current.NAME]=current.sha(current.NAME);hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=report['source_family'],
        exact_power_graph_nodes_checked=len(g.nodes),unchanged_original6601_node_source_prefix=True,
        source_recipe_and_native_Rp_acceptance_hashes_checked=True,independent_graph_calculus=calculus,
        actual_original_P0_and_all_absolute_incoming_histories_retained=True,
        same_parent_semigroup_and_2Rc_to_Rp_geometry_proved=True,
        native_constructor_source_equality_not_inferred_from_units_or_same_family=True,
        no_source_integration_or_ancestor_owner_replayed=True,
        actual_numerics_preheat_heat_high_jets_and_recursion_not_claimed=True,
        input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_LIMIT_POWER_TO_RP',len(g.nodes),'nodes; independent complete/Z flow and pulse units',flush=True)
    return receipt


if __name__=='__main__':run()

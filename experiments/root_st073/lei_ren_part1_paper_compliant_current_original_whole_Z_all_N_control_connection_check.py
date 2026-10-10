"""Check actual geometry, directed integer proof and the source-bound C1 map.

This does not evaluate the astronomical integer's point phase or controls.
Independent symbolic identities and current source recipe witnesses cover
the changed connection without replaying the existing 228 integrations.
"""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection as current


def forbidden(*args,**kwargs):raise AssertionError('Legacy 24-cell graph or fixed-N transport used')


def symbolic_checks():
    offsets,domains,x=current.exact_geometry();checks={}
    require=current.require
    for left,right in zip(current.CHARTS,current.CHARTS[1:]):
        difference=offsets[left].subs(x,domains[left][1])-offsets[right].subs(x,domains[right][0])
        require(sy.simplify(difference)==0,'Original radius seam failed: '+left+'/'+right)
        checks[left+'->'+right]=True
    require(sy.simplify(offsets['first_micro'].subs(x,domains['first_micro'][0]))==0,'Symbolic inlet phase is not zero')
    require(sy.diff(offsets['actual_patch'],x)==1/x,'Patch measure must be dx/x')
    require(sy.diff(offsets['O2_axial'],x)==40*sy.exp(40*x),'Axial physical Jacobian missing')
    for chart in ('restoration','postrestore','Rh_reference','O2_slope','O2_buffer','O3_transition','O3_power'):
        require(sy.diff(offsets[chart],x)==1,'Original unit log chart changed')
    # Check Poisson partial derivative and exact original primitive identities.
    r,psi=sy.symbols('r psi',real=True)
    D=1-2*r*sy.cos(psi)+r*r;nn=sy.cos(psi)-r
    require(sy.simplify(sy.diff(nn/D,r)-(-D+2*nn*nn)/D**2)==0,'Original w_r changed')
    u=sy.symbols('u',real=True)
    require(sy.simplify(sy.diff(u/sy.sqrt(1+u*u),u)-(1+u*u)**sy.Rational(-3,2))==0,'Original r_u changed')
    E,V,A,B,N=sy.symbols('E V A B N',nonzero=True)
    U=E*N*(sy.exp(A/N)-1)
    increments=dict(m=B,h=U,k=V*U+E*B+U*B/N,
        e=2*V*B+B*B/N-E*U-U*U/(2*N),p=E*U+U*U/(2*N))
    old=dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
    new={k:q.subs({E:E+U/N,V:V+B/N},simultaneous=True) for k,q in old.items()}
    for key in old:
        require(sy.expand(N*(new[key]-old[key])-increments[key])==0,'Normalized density identity failed: '+key)
    z=sy.symbols('Z');a=sy.Function('A')(z);h=sy.Function('H')(z);mu=sy.symbols('mu',positive=True)
    for degree in (1,2):
        require(sy.simplify(sy.diff(h/a**degree,z)-(sy.diff(h,z)-degree*sy.diff(a,z)*h/a)/a**degree)==0,
            'True relative target quotient Z rule failed')
    checks.update(symbolic_inlet_phase_zero=True,patch_Jacobian=True,axial_Jacobian=True,
        Poisson_w_r=True,Poisson_r_u=True,normalized_five_density_identities=True,true_relative_Z_quotients=True)
    return checks


def recipe_witnesses(parameter_proof):
    outer=current.current;owner=outer.WholeZAllNOuterRcFunctions();ends=outer.CELLS[0]
    assert current.parameter_identity(owner)==parameter_proof
    # Explicit source dispatch: each current provider keeps its own coordinate guard.
    switch=outer.previous.previous
    switch_view=SimpleNamespace(**vars(owner))
    switch_view.coordinate=lambda value:switch.WholeZAllNSwitchFunctions.coordinate(switch_view,value)
    long=outer.previous
    long_view=SimpleNamespace(**vars(owner))
    long_view.coordinate=lambda value,chart=None:long.WholeZAllNLongPatchFunctions.coordinate(long_view,value,chart)
    points=dict(first_micro=(1,2),second_micro=(3,2),frozen_macro=(1,2),
        first_switch=(1,2),second_switch=(1,2),post_power=(1,2),
        long_reshape=(1,2),reference=(1,2),restoration=(1,2),postrestore=(1,2),actual_patch=(3,2),
        Rh_reference=(-2,1),O2_slope=(1,2),O2_axial=(1,2),O2_buffer=(5,1),O3_transition=(1,2),O3_power=(2,1))
    witnesses=[]
    for chart in current.CHARTS:
        point=points[chart]
        if chart in current.CHARTS[:3]:packet=owner.owner.source_query(ends,chart,point)
        elif chart in current.CHARTS[3:6]:packet=switch.WholeZAllNSwitchFunctions.leading_packet(switch_view,ends,chart,point)
        elif chart in current.CHARTS[6:11]:packet=long.WholeZAllNLongPatchFunctions.leading_packet(long_view,ends,chart,point)
        else:packet=owner.leading_packet(ends,chart,point)
        op=owner.owner.owner(ends)
        assert packet['source_identity']==owner.identity and packet['candidate_N']==owner.N0
        assert packet['exact_common_P0_axial5'] is op.P0 and packet['actual_phase_Z_exact_zero']
        generic=packet['original_generic_source'];proof=packet['original_full_source_quotients'];roots=packet['original_roots']
        quantities={key:generic['common_velocity_'+key+'_axial5'][:2] for key in ('E','V')}
        quantities.update({key:[roots[key][(0,k)] for k in (0,1)] for key in ('a','t0','p2')})
        if chart=='long_reshape':
            assert proof['exact_b_and_t0_zero']
            quantities.update(b=[op.flow.scalar(0)]*2,Delta=proof['Delta_axial5'][:2])
        else:
            bkey='actual_nonzero_b_axial5' if chart in ('reference','restoration','postrestore','Rh_reference') else 'actual_b_axial5'
            quantities.update(b=proof[bkey][:2],Delta=proof['actual_Delta_axial5'][:2])
        assert outer.previous.equivalent_rows(quantities['E'],[roots['E'][(0,k)] for k in (0,1)])
        outer.core.relative.same_source(op.flow,[q for values in quantities.values() for q in values])
        if chart=='O3_power':assert all(q.zero for q in packet['original_q_C0_Z'].values())
        witnesses.append(dict(chart=chart,source_quantities=list(quantities),same_P0_and_actual_phase_Z_zero=True,
            b_and_Delta_from_actual_quotient_proof=True,field_values_not_read_from_caps=True))
    return witnesses


def run():
    began=time.monotonic();saved,hashes=current.load_inputs()
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert report[current.GATE] and report['source_family']==saved['source_family']
    for name,digest in report['input_hashes'].items():assert current.sha(name)==digest,name
    c=MPIntervalContext();c.dps=500
    with mp.workdps(540):
        chosen,frequency=current.frequency_connection(saved,c)
        assert current.current_range_containment(saved)==report['exact_current_source_uniform_range_containment']
        assert current.packets.encode(frequency)==report['actual_whole_Z_frequency_connection']
        assert chosen.J==1358356628656378313
        proof=frequency['selected_integer'];assert proof['log_N_lower_ge_required_upper'] and proof['N_ge_N0']
        for bad in (True,12.0,'12',-1):
            try:current.DoubleDyadicInteger(bad)
            except ValueError:pass
            else:raise AssertionError('Non-exact integer descriptor admitted')
        # Cover a small exact selection independently; do not materialize the live N.
        small=current.choose_integer(c,257,c.mpf(100));assert small.J==12
        assert current.ep(small.log_interval(c))[0]>=100
        manifest=current.recipe_manifest();assert manifest==report['exact_source_recipe_manifest']
        with patch.object(current.algebra,'exact_all_N_graph',forbidden),patch.object(current.source,'exact_radius_maps',forbidden):
            built=current.build_graph(saved['source_family'],saved['fixed_leading_input_N0'],chosen,manifest,
                report['canonical_original_source_parameter_proof'])
        assert built['graph'].nodes==report['exact_function_graph_nodes']
        assert current.encode_graph(built)==report['exact_current_integral_and_control_graph']
        windows=built['windows'];assert [w['chart'] for w in windows]==list(current.CHARTS)
        assert all(pair.value==built['graph'].zero and pair.Z==built['graph'].zero for pair in windows[0]['incoming'].values())
        for row in windows:
            assert row['memory']['p']==built['graph'].one
            phase=built['graph'].nodes[row['phase']]
            product=built['graph'].nodes[phase['argument']]
            assert phase['name']=='fractional_part' and built['N'].node in product['arguments']
            assert row['predecessor_memory_never_reset']
            partition=current.native_partitions()[row['chart']]
            assert [[p['exact_left'],p['exact_right']] for p in row['actual_current_source_cells']]==current.current.serialized(list(zip(partition,partition[1:])))
        for row in built['graph'].nodes:
            if row['operation']=='original_lazy_flat_branch' and row['original_collar_coordinate'] is not None:
                assert built['graph'].nodes[row['flat_predicate']]['operation']=='logical_or'
        for before,after in zip(windows,windows[1:]):assert before['outgoing']==after['incoming']
        assert windows[-1]['quiet_local_source_zero']
        assert all(p.value==built['graph'].zero and p.Z==built['graph'].zero for p in windows[-1]['contributions'].values())
        assert any(p.value!=built['graph'].zero for p in windows[-1]['outgoing'].values())
        assert len(built['finite_picard_sequence'])==4
        assert all(type(p.value) is current.source.FunctionRef for p in built['N_scaled_targets'].values())
        assert current.packets.encode([current.factored_tail(frequency['C1_limit'],n) for n in (0,3,32)])==report['exact_C1_tail_examples']
        symbolic=symbolic_checks();witnesses=recipe_witnesses(report['canonical_original_source_parameter_proof'])
    open_keys=(*current.current.OPEN,'actual_global_frequency_admitted','actual_numeric_controls_evaluated',
        'actual_terminal_Z_function_closure_installed','actual_power_repair_band_source_extended',
        'physical_original_exterior_five_targets_closed')
    for key in open_keys:assert report[key] is False,key
    input_hashes=dict(report['input_hashes']);input_hashes[current.NAME]=current.sha(current.NAME)
    input_hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=report['source_family'],
        exact_current_source_charts_checked=17,exact_radius_seams_symbolically_checked=16,
        exact_source_recipe_witnesses=witnesses,symbolic_identities=symbolic,
        exact_function_graph_nodes_checked=len(built['graph'].nodes),
        directed_integer_expression_and_actual_C1_contraction_checked=True,
        current_source_range_cells_already_bound=228,no_existing_source_integrations_replayed=True,
        numeric_phase_controls_and_global_frequency_not_claimed=True,
        input_hashes=input_hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_17_CHART_CONTROL_CONNECTION',len(built['graph'].nodes),'nodes; 16 seams; 17 source witnesses',flush=True)
    return receipt


if __name__=='__main__':run()

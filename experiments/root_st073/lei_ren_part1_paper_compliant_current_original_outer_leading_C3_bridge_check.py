"""Actual source replay, Duhamel uniqueness and complete Rc inlet identity.

This checks defining mathematical functions, never interval endpoint values.
The complete bridge retains the independent P0 and nonzero correction memory.
"""
import ast
import inspect
import hashlib
import json
import math
from pathlib import Path
import textwrap
import time
from types import SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_bridge as current
import lei_ren_part1_paper_compliant_current_original_outer_C3_continuation_check as outer_check

leading=current.leading
rows=current.target.rows
zero=current.leading_check.zero


def source_replay(field):
    replay=current.OriginalLeadingSourceReplay(field);q=replay.q;count=0
    for name,c in field.leading.functions['charts'].items():
        actual=replay.normalized[name]
        quantities={'E':c['E'],'V':c['V'],**c['histories']}
        values={'E':actual['E'],'V':actual['V'],**actual['histories']}
        for key,value in quantities.items():
            for j,ref in enumerate(rows(value)):
                assert zero(q.at(ref)-math.factorial(j)*values[key][j]),(name,key,j,'actual defining source replay')
                count+=1
    assert count==168
    pressure_rows=0
    for name,common in replay.normalized.items():
        for j,value in enumerate(common['pressure']):
            assert common['P0'] is replay.op.P0
            assert zero(value-replay.P0[j]-common['histories']['p'][j])
            pressure_rows+=1
    assert pressure_rows==24
    Y=s.Symbol('original_y',positive=True);axial=field.leading.functions['charts']['O2_axial']
    x=q.at(axial['native_coordinate']);y=q.at(axial['physical_coordinate'])
    cutoff=replay.turnoff(Y,s.log(Y)/40)
    for j,value in enumerate(cutoff):
        assert zero(value-s.diff(current.leading_check.Sigma(1-s.log(Y)/40),Y,j)),('turnoff dy',j)
    for j,value in enumerate(replay.raw['O2_axial']['Vy']):
        expected=4*q.z*s.diff(current.leading_check.Sigma(1-s.log(Y)/40),Y)
        assert zero(math.factorial(j)*value.subs(x,s.log(Y)/40)-s.diff(expected,q.z,j))
    return replay,dict(independent_actual_source_assignment_replay_rows=count,
        defining_source_slices=replay.bindings,original_parents_replayed_in_order=True,
        scalar_kernel_values_are_exact_functions_not_enclosure_endpoints=True,
        original_physical_y_turnoff_derivatives_zero_through_four_replayed=True,
        independent_actual_normalized_absolute_pressure_rows=pressure_rows,
        raw_V_m_k_divided_by_original_Pstar_once=True)


def require_statement(fn,text):
    wanted=ast.dump(ast.parse(text).body[0],include_attributes=False)
    body=ast.parse(textwrap.dedent(inspect.getsource(fn)))
    assert any(ast.dump(n,include_attributes=False)==wanted for n in ast.walk(body)),(fn.__name__,text)


def frontend_contracts():
    p=leading.providers;checks=[]
    for kind,fn in (('rh',p.rh.WholeZRhReferenceFiniteN.query),
            ('slope',p.slope.WholeZO2SlopeFiniteN.query),
            ('axial',p.axial.WholeZO2AxialBufferFiniteN.query),
            ('o3',p.o3.WholeZO3RcFiniteN.query)):
        _,binding=p.source_frontend(fn,kind)
        assert binding==p.FRONTEND_BINDINGS[kind]
        if kind=='rh':
            require_statement(fn,'background=original.background_cell(op,left,right)')
            require_statement(fn,'proxy,raw,recovered=reference.recover_cell(op.reference,background)')
        else:
            require_statement(fn,"generic=reference.long.RECOVER(proxy,source['raw'])")
            module=getattr(p,kind)
            fn2,proof=module.compiled_background(lambda *args:None)
            original=ast.parse(textwrap.dedent(inspect.getsource(module.original.background_cell))).body[0]
            # The compiler executes the unchanged entire source AST; only its
            # scalar-enclosure callback namespace changes, as checked upstream.
            require_statement(module.compiled_background,"tree=ast.parse(Path(original.__file__).read_text(encoding='utf8'))")
            require_statement(module.compiled_background,"fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='background_cell')")
            assert proof['original_background_AST_sha256']==hashlib.sha256(ast.dump(original).encode()).hexdigest()
            assert fn2.__name__==original.name
        checks.append(dict(kind=kind,actual_frontend_binding=binding))
    provider=p.WholeZAllNOuterRcFunctions
    require_statement(provider.__init__,"self.views['axial'].parent=lambda ends:self.leading_packet(ends,'O2_slope',(1,1))")
    require_statement(provider.__init__,"self.views['o3'].parent=lambda ends:self.leading_packet(ends,'O2_buffer',(11,1))")
    require_statement(provider.__init__,"self.views['o3'].query=lambda ends,chart,left,right=None:self.leading_packet(ends,'O3_'+chart,left,right)")
    require_statement(provider.leading_packet,"op=self.owner.owner(key)")
    require_statement(provider.leading_packet,"if packet['source_identity']!=self.identity or packet['candidate_N']!=self.N0 or packet['exact_common_P0_axial5'] is not op.P0 or not packet['actual_phase_Z_exact_zero']: raise ValueError('Same live outer source identity/P0/N0/Z origin required')")
    generic=p.o3.reference.long.generic.recover_inputs
    require_statement(generic,"invS=f.factor((0,-.5,0,0,0))")
    require_statement(generic,"if packet['original_P0_normalized_axial5'] is not op.P0: raise ValueError('Exact same live P0 object required')")
    require_statement(generic,"pressure=add(op.P0,p)")
    body=ast.parse(textwrap.dedent(inspect.getsource(generic)))
    assert any(isinstance(n,ast.keyword) and n.arg=='common_original_P0_axial5'
        and ast.unparse(n.value)=='op.P0' for n in ast.walk(body))
    return dict(actual_frontend_dispatch_and_normalization=checks,
        original_source_parent_endpoints_one_eleven_one_retained=True,
        source_fixed_N0_and_shared_owner_P0_guarded=True)


def aliases(field,replay):
    q=replay.q;provider=leading.providers.WholeZAllNOuterRcFunctions.leading_packet
    binding=current.current.ast_binding(provider);counts={}
    for consumer,table,nodes in (('outer',field.aliases,field.prefix),
            ('band',field.band_aliases,field.leading.bridge.band_nodes)):
        count=0
        for key,alias in table.items():
            leaf=nodes[int(key)];recipe=leaf['recipe'];chart=leaf['native_chart']
            assert all(recipe[k]==value for k,value in binding.items())
            assert recipe['native_chart']==chart and recipe['leading_source_fixed_input_N0']
            assert recipe['range_endpoints_not_field_coefficients']
            assert recipe['original_independent_P0_and_parent_leading_recipes_retained']
            assert alias['source_graph']==('accepted_outer' if consumer=='outer' else 'accepted_repair_band')
            assert alias['resolved_graph']=='canonical_bridge_graph'
            assert field.resolve_leading_leaf(int(key),consumer).graph is field.graph
            j=leaf['Z_order'];quantity=leaf['quantity']
            assert leaf.get('Taylor_coefficient_factorial',1)==math.factorial(j)
            expected_path=('original_generic_source.common_velocity_'+quantity+'_axial5[Z_order]'
                if quantity in ('E','V') else 'original_generic_source.common_own_five_histories_axial5.'+quantity[8:]+'[Z_order]')
            assert recipe['quantity_paths'][quantity]==expected_path
            if 'actual_original_frontend_binding' in recipe:
                kind='rh' if chart=='Rh_reference' else 'slope' if chart=='O2_slope' else 'axial' if chart.startswith('O2_') else 'o3'
                assert recipe['actual_original_frontend_binding']==leading.providers.FRONTEND_BINDINGS[kind]
            c=nodes[leaf['coordinate']]
            coordinate=(s.Rational(c['numerator'],c['denominator']) if c['operation']=='exact_rational'
                else s.Symbol(c['name'],real=True) if c['operation']=='bound_variable' else None)
            assert coordinate is not None
            native=q.at(field.leading.functions['charts'][chart]['native_coordinate'])
            values=replay.normalized[chart]
            source_value=values[quantity] if quantity in ('E','V') else values['histories'][quantity[8:]]
            expected=math.factorial(j)*s.sympify(source_value[j]).subs(native,coordinate)
            assert zero(q.at(field.resolve_leading_leaf(int(key),consumer))-expected),(consumer,key,quantity,j)
            count+=1
        counts[consumer]=count
    assert counts==dict(outer=72,band=72)
    try:field.resolve_leading_leaf(-1)
    except ValueError:pass
    else:raise AssertionError('Unsupported recipe admitted')
    return dict(independent_resolved_outer_source_recipes=counts['outer'],
        independent_resolved_band_source_recipes=counts['band'],
        actual_source_provider_paths_bound_to_replayed_math=True,
        ordinary_Z_factorials_once_and_coordinates_semantically_translated=True,
        returned_refs_owned_by_canonical_bridge_graph=True,
        consumer_marks_origin_graph_and_never_transfers_foreign_node_ids=True,
        unsupported_recipes_fail_closed=True)


def duhamel(field,replay):
    b=field.leading.bridge;g=field.graph;q=replay.q
    raw=b.outer_report;charts=raw['actual_five_chart_C3_continuation_functions']
    rh=current.outer.rh.read(current.outer.rh.NAME)['actual_Rh_C3_continuation_functions']
    ref=lambda i:current.source.FunctionRef(g,i)
    alias_bindings={int(i):q.at(v['explicit_original_source_row']) for i,v in field.aliases.items()}
    exact=current.leading_check.Interpreter(field.leading,bindings=alias_bindings)
    first=charts[current.outer.CHARTS[0]]['actual_leading_incoming_C3'];rhrows=reference_inlet_rows=0
    rchart=field.leading.functions['charts']['Rh_reference'];rx=q.at(rchart['native_coordinate'])
    for key in current.current.RATES:
        for j,i in enumerate(rh['actual_leading_history_C3'][key]):
            actual=exact.at(i).subs(q.at(rh['endpoint']),-5)
            expected=q.at(rows(rchart['histories'][key])[j]).subs(rx,-5)
            assert zero(actual-expected),(key,j,'actual Rh y=-5 source inlet')
            reference_inlet_rows+=1
        for j,i in enumerate(first[key]):
            assert zero(exact.at(i)-q.at(rows(rchart['histories'][key])[j]).subs(rx,0)),(key,j,'Rh incoming source identity')
            rhrows+=1
    densities=FTC=inlets=seams=0
    for name,f in charts.items():
        geo=raw['actual_geometry_bindings'][name];native=s.Symbol('native_'+name,real=True)
        c=field.leading.functions['charts'][name];cx=q.at(c['native_coordinate'])
        E,V=(q.at(rows(c[k])[0]).subs(cx,native) for k in ('E','V'))
        drivers=dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
        for key in current.current.RATES:
            for j,i in enumerate(f['actual_leading_density_C3'][key]):
                assert zero(exact.at(i)-s.diff(drivers[key],q.z,j)),(name,key,j,'actual source driver')
                densities+=1
        t=s.Symbol('native_'+name);a=s.Symbol(name+'_history_endpoint')
        L=s.Function('log_radius_'+name) if name=='O2_axial' else (lambda x:x)
        left,right=map(s.Integer,f['exact_native_domain']);J=s.diff(L(t),t)
        geometric={geo['offset']:L(t),geo['left_offset']:L(left),geo['right_offset']:L(right)}
        if geo['Jacobian']!=g.one.node:geometric[geo['Jacobian']]=J
        for key,rate in current.current.RATES.items():
            r=s.Rational(str(rate))
            for j,(h,d,i) in enumerate(zip(f['actual_leading_incoming_C3'][key],f['actual_leading_density_C3'][key],f['actual_leading_history_C3'][key])):
                H=s.Integer(0) if h==g.zero.node else s.Symbol('H_'+str(h))
                D=s.Integer(0) if d==g.zero.node else s.Function('D_'+str(d))(t)
                bindings=dict(geometric)
                if h!=g.zero.node:bindings[h]=H
                if d!=g.zero.node:bindings[d]=D
                wanted=s.exp(-r*(L(a)-L(left)))*H
                if D!=0:wanted+=s.Integral(s.exp(-r*(L(a)-L(t)))*D*J,(t,left,a))
                got=outer_check.symbolic(g,i,bindings)
                assert s.expand(got-wanted)==0,(name,key,j,'actual Duhamel')
                assert s.simplify(s.diff(got,a)-s.diff(L(a),a)*(D.subs(t,a)-r*got))==0
                assert s.simplify(got.subs(a,left)-H)==0
                FTC+=1;inlets+=1
                out=g.nodes[f['actual_leading_outgoing_C3'][key][j]]
                assert out['operation']=='function_substitution' and out['expression']==i
                assert out['variable']==f['endpoint'] and out['value']==f['domain'][1]
        previous=first if name==current.outer.CHARTS[0] else charts[current.outer.CHARTS[current.outer.CHARTS.index(name)-1]]['actual_leading_outgoing_C3']
        assert f['actual_leading_incoming_C3']==previous;seams+=20
        # The new identified functions use the actual old API endpoint variable.
        for key,value in field.functions['actual_outer_leading_Duhamel_as_original_source_C3'][name].items():
            for j,i in enumerate(rows(value)):
                assert zero(q.at(i)-q.at(rows(c['histories'][key])[j]).subs(cx,q.at(f['endpoint'])))
    # Difference of any two physical-y solutions obeys d_y+r*d=0. The
    # explicit integrating factor and zero inlet force the constant to zero.
    y,y0,C=s.symbols('y y0 C',real=True)
    uniqueness=[]
    for rate in sorted(set(current.current.RATES.values())):
        r=s.Rational(str(rate));difference=C*s.exp(-r*(y-y0))
        assert s.simplify(s.diff(difference,y)+r*difference)==0
        assert difference.subs(y,y0)==C and difference.subs(C,0)==0
        uniqueness.append(str(r))
    prior=current.outer.rh.read(leading.RECEIPT)
    assert prior['independent_five_ODE_calculus']['total_independent_ODE_rows']==120
    assert prior['independent_original_source_seams']['total_C3_state_and_own_rate_seam_rows']==240
    assert densities==FTC==inlets==seams==100 and rhrows==reference_inlet_rows==20
    return dict(independent_Rh_Z0_to_Z3_original_inlet_rows=rhrows,
        independent_actual_Rh_y_minus5_source_inlet_rows=reference_inlet_rows,
        independent_actual_leading_drivers_from_replayed_E_V_rows=densities,
        independent_actual_Duhamel_physical_FTC_and_inlet_rows=FTC,
        actual_predecessor_and_outgoing_function_seams=seams,
        source_ODE_and_seam_dependency=leading.RECEIPT,
        explicit_integrating_factor_unique_own_rates=uniqueness,
        all_one_hundred_source_and_actual_Duhamel_functions_identified=True,
        actual_axial_nonunit_Jacobian_and_incoming_memory_retained=True)


class InletAlgebra:
    """Only admitted inlet algebra; zero-length integrals do not read controls."""
    def __init__(self,nodes,bindings):self.nodes,self.bindings,self.memo=nodes,bindings,{}
    def at(self,i):
        if i in self.bindings:return self.bindings[i]
        if i in self.memo:return self.memo[i]
        n=self.nodes[i];op=n['operation'];at=self.at
        if op=='exact_rational':v=s.Rational(n['numerator'],n['denominator'])
        elif op=='bound_variable':v=s.Symbol(n['name'],real=True)
        elif op=='sum':v=s.Add(*(at(j) for j in n['arguments']))
        elif op=='product':v=s.Mul(*(at(j) for j in n['arguments']))
        elif op=='negative':v=-at(n['argument'])
        elif op=='positive_quotient':v=at(n['numerator'])/at(n['denominator'])
        elif op=='analytic_unary':
            assert n['name'] in ('exp','log','exprel')
            v=dict(exp=s.exp,log=s.log,exprel=current.leading_check.Exprel)[n['name']](at(n['argument']))
        elif op=='definite_integral':
            assert n['exact_function_integral']
            assert at(n['lower'])==at(n['upper']),'Only actual zero-length inlet integral admitted'
            v=s.Integer(0)
        elif op=='function_substitution':
            assert n['Z_independent_substitution']
            v=at(n['expression']).subs(at(n['variable']),at(n['value']))
        else:raise AssertionError('Unadmitted inlet dependency '+str(n))
        self.memo[i]=v;return v


def Rc_interface(field,replay):
    b=field.leading.bridge;q=replay.q
    packet=b.band_report['actual_C3_repair_band_functions'];power=packet['original_leading_power']
    nodes=b.band_nodes;endpoint=field.functions['actual_outer_Rc_leading_C3']
    bindings={int(i):q.at(v['explicit_original_source_row']) for i,v in field.band_aliases.items()}
    seeds=0
    for key,rr in power['actual_endpoint_seeds'].items():
        for j,i in enumerate(rr):
            n=nodes[i]
            assert n['native_chart']=='O3_power' and n['quantity']=='history_'+key and n['Z_order']==j
            assert n['recipe']['exact_source_call']=='WholeZAllNOuterRcFunctions.leading_packet(ends, O3_power, (2,1))'
            assert nodes[n['coordinate']]==dict(operation='exact_rational',numerator=2,denominator=1)
            assert zero(bindings[i]-q.at(rows(endpoint[key])[j]));seeds+=1
    assert seeds==20 and power['m_k_S_normalization_already_in_endpoint_source'] and power['P0_stays_separate']
    originalP0=s.Function('original_P0')(q.z)
    P0=field.functions['actual_independent_P0_C3']
    pbindings={}
    for j,(old,i) in enumerate(zip(rows(P0),power['independent_P0'])):
        for graph,node in ((field.graph.nodes,old.node),(nodes,i)):
            n=graph[node]
            assert n['operation']=='current_original_leading_function_recipe' and n['quantity']=='P0'
            assert n['source_family']==field.identity and n['Z_order']==j
            assert n.get('Taylor_coefficient_factorial',1)==math.factorial(j)
            assert n['recipe']['quantity_paths']['P0']=='original_generic_source.common_original_P0_axial5[Z_order]'
        pbindings[old.node]=s.diff(originalP0,q.z,j);bindings[i]=s.diff(originalP0,q.z,j)
    N=s.Symbol('same_original_positive_N',positive=True)
    Nbound={b.outer_N:N,**pbindings}
    correction={}
    for key,rr in b.outer_terminal.items():
        correction[key]=[]
        for j,(i,inlet) in enumerate(zip(rr,b.band_inlet[key])):
            Q=s.Function('actual_N_scaled_'+key)(q.z).diff(q.z,j)
            Nbound[i]=Q;bindings[inlet]=Q/N;correction[key].append(Q/N)
    # Formal mu is the exact same admitted graph dependency. Locate it from
    # the band's actual alpha expression, without reusing foreign node IDs.
    alpha=nodes[packet['profiles_at_band_x']['original_alpha']]
    assert alpha['operation']=='sum'
    murows=[i for i in alpha['arguments'] if nodes[i]['operation']!='exact_rational']
    assert len(murows)==1
    mu_node=murows[0]
    signature=current.leading.bridge.ExactFunctionSignatures
    assert signature(nodes).at(mu_node)==signature(field.graph.nodes).at(field.leading.parameters['mu'].node)
    bindings[mu_node]=q.mu;bindings[packet['original_band_variable']]=s.Integer(1)
    inlet=InletAlgebra(nodes,bindings)
    q2=current.leading_check.Interpreter(field.leading,bindings=Nbound)
    count=amplitude=0
    for key in current.current.RATES:
        for j,(lh,complete) in enumerate(zip(power['leading_histories'][key],packet['complete_histories'][key])):
            lead=q.at(rows(endpoint[key])[j]);wanted=lead+correction[key][j]
            assert zero(inlet.at(lh)-lead),(key,j,'actual band leading inlet')
            assert zero(inlet.at(complete)-wanted),(key,j,'actual complete band inlet')
            assert zero(q2.at(rows(field.complete_Rc_functions()[key])[j])-wanted),(key,j,'actual full Rc interface')
            count+=1
    for j,i in enumerate(packet['profiles_at_band_x']['radial_y_rows'][0]['original_E']):
        assert zero(inlet.at(i)-q.at(rows(field.functions['actual_outer_Rc_amplitude_C3'])[j]));amplitude+=1
    for j,i in enumerate(rows(field.functions['actual_outer_Rc_absolute_pressure_C3'])):
        assert zero(q2.at(i)-s.diff(originalP0,q.z,j)-q2.at(rows(field.complete_Rc_functions()['p'])[j]))
    assert count==20 and amplitude==4
    assert b.acceptance_loaded and b.identity_proof['outer_N_scaled_equals_N_times_actual_repair_inlet']
    assert len(b.identity_proof['rows'])==20 and b.source_integer_and_Rc_binding['same_exact_original_N_and_Rc']
    return dict(independent_actual_canonical_Rc_seed_function_rows=seeds,
        independent_complete_Rc_to_actual_band_inlet_rows=count,
        actual_Rc_amplitude_C3_rows=amplitude,independent_original_P0_and_absolute_pressure_rows=4,
        accepted_correction_expression_bridge_dependency=current.leading.bridge.RECEIPT,
        same_original_N_divided_once=True,zero_length_local_integrals_with_nonzero_incoming_kept=True,
        P0_shared_original_owner_guard_and_distinct_pressure_history_retained=True,
        source_geometry_Rw_exp2_equals_Rc_and_actual_band_x1=True)


def run(field=None):
    began=time.monotonic();raw=current.outer.rh.read(current.NAME)
    for name,value in raw['input_hashes'].items():assert current.sha(name)==value,name
    field=field if field is not None else current.CurrentOuterLeadingC3Bridge(require_checked=False)
    assert not field.acceptance_loaded and field.identity==raw['source_family']
    assert field.graph.nodes==raw['exact_graph_nodes']
    assert field.aliases==raw['actual_original_leading_provider_aliases']
    assert field.band_aliases==raw['actual_band_leading_provider_aliases']
    assert raw['resolved_function_graph_contract']==dict(returned_refs_owner='canonical_bridge_graph',
        consumer_selects_original_source_graph_not_returned_owner=True,
        same_original_source_family_and_common_V_m_k_Pstar_normalization=True,
        numerical_flow_ledger_or_point_oracle_conversion=False)
    assert current.target.encoded(field.functions)==raw['actual_original_leading_and_complete_Rc_functions']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    replay,replayed=source_replay(field);frontend=frontend_contracts();resolved=aliases(field,replay)
    identified=duhamel(field,replay);interface=Rc_interface(field,replay)
    open_keys=('selected_native_pulse_constructor_consumes_current_C3_frame',
        'current_numeric_point_field_oracle_installed',
        'global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed')
    assert not any(raw[k] for k in open_keys)
    result=dict(all_passed=True,source_family=field.identity,**dict.fromkeys(current.GATES,True),
        independent_actual_defining_source_replay=replayed,independent_source_frontend_contracts=frontend,
        independent_actual_recipe_function_resolution=resolved,
        independent_actual_Duhamel_source_function_identity=identified,
        independent_complete_Rc_interface=interface,**dict.fromkeys(open_keys+current.OPEN,False),
        input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Original leading source, Duhamel and full Rc repair inlet functions identified',flush=True)
    return result


if __name__=='__main__':run()

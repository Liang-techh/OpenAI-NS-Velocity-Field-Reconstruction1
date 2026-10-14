"""Actual source ODE transport from Rh through current Rp to the pulse."""
import copy
import json
from pathlib import Path
import time
from types import SimpleNamespace as NS
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source as current
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity_check as old_check
import lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check as compact_check


def zero(value):return old_check.zero(value)


class Expression:
    def __init__(self,g,bindings):self.g,self.bindings,self.cache=g,dict(bindings),{}
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings:return self.bindings[i]
        if i in self.cache:return self.cache[i]
        row=self.g.nodes[i];op=row['operation'];at=self.at
        if op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
        elif op=='bound_variable':value=sy.Symbol(row['name'],real=True)
        elif op=='sum':value=sy.Add(*(at(n) for n in row['arguments']))
        elif op=='product':value=sy.Mul(*(at(n) for n in row['arguments']))
        elif op=='negative':value=-at(row['argument'])
        elif op=='positive_quotient':value=at(row['numerator'])/at(row['denominator'])
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            value=at(row['expression']).subs(at(row['variable']),at(row['value']))
        elif op=='analytic_unary':
            x=at(row['argument']);name=row['name']
            if name=='exprel':value=sy.Integer(1) if x==0 else (sy.exp(x)-1)/x
            elif name in ('exp','log','sqrt'):value=getattr(sy,name)(x)
            else:raise AssertionError('Unbound actual Rp unary '+name)
        elif op=='exact_original_source_flat_sigma':value=old_check.Sigma(at(row['argument']))
        else:raise AssertionError('Unbound actual Rp source '+op+' '+str(i))
        self.cache[i]=value;return value


def bind(bindings,ids,body,Z):
    for n,i in enumerate(ids):
        value=sy.diff(body,Z,n)
        if i==0:assert value==0
        elif i==1:assert value==1
        else:bindings[i]=value


def old_context(owner):
    data=owner.inputs[current.old_frame.NAME];leading=owner.inputs[current.old_leading.NAME]
    assert leading['exact_graph_nodes']==data['exact_graph_nodes'][:len(leading['exact_graph_nodes'])]
    g=current.source.FunctionTransportGraph();g.nodes=copy.deepcopy(data['exact_graph_nodes'])
    g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
    parameters=owner.inputs[current.old_frame.quiet.NAME]['actual_C3_repaired_exit_power_and_pulse_functions']['parameters']
    field=NS(graph=g,contracts=current.old_leading.source_contracts(),
        parameters={k:current.source.FunctionRef(g,i) for k,i in parameters.items()})
    q=old_check.Interpreter(field);f=data['actual_Rp_native_frame_functions']
    Tw=sy.Symbol('Tw',positive=True);q.bindings[f['Tw']]=Tw
    datum=data['actual_common_P0_definition']
    assert datum['passed'] and datum['same_fourteen_continuous_atoms_and_analytic_flatten_function']
    assert datum['current_order6_first_six_and_native_order5_same_defining_rows']
    assert datum['common_P0_kept_independent_of_cumulative_pressure']
    P0=sy.sympify(datum['exact_analytic_function'],locals={'Z':q.z,'F_flat':sy.Function('F_flat')})
    for i,row in enumerate(g.nodes):
        if row['operation']=='current_original_leading_function_recipe' and row.get('quantity')=='P0':
            assert row['source_family']==owner.family
            q.bindings[i]=sy.diff(P0,q.z,row['ordinary_slow_Z_derivative_order'])
    return q,leading['actual_original_leading_C3_source_functions']['charts'],f,Tw,P0


def endpoint(q,body,x,value):
    out=body.subs(x,value)
    for i,fn in q.kernel_functions.items():out=out.subs(fn(q.at(q.g.nodes[i]['lower'])),0)
    return out


def source_inputs(owner):
    q=owner.inputs[current.quiet.NAME];outer=owner.inputs[current.outer.NAME]
    compact=owner.inputs[current.compact.NAME];g=owner.g
    assert q['incoming']==compact['terminal']
    assert q['relative_incoming']==compact['terminal']['correction']
    assert q['P0']==compact['P0'] and q['mu']==compact['mu']==outer['mu']
    assert q['N']==compact['N']==outer['N']
    assert q['Pstar']==compact['Pstar']==owner.seam.inlet.Pstar.node
    assert set(q['Rp_frame'])==set(current.FRAME_KEYS)
    assert current.previous.encode(owner.current_frame)==q['Rp_frame']
    assert all(len(current.previous.rows(jet))==6 for jet in owner.current_frame.values())
    for pressure,axis,moment in zip(q['Rp_frame']['pressure'],q['Rp_frame']['P0'],q['Rp_frame']['Mp']):
        assert owner.ref(pressure)==g.add(owner.ref(axis),owner.ref(moment))
    seed=outer['seed_identity'];Rh=current.parent.linear.read(current.outer.original_Rh.RECEIPT)
    branch=current.parent.linear.read(seed['true_C6_endpoint_receipt'])
    assert Rh['all_passed'] and Rh[current.outer.original_Rh.GATES[0]]
    assert Rh['current_full_weight_map_and_source_AST_checked'] and Rh['interval_overlap_not_used_as_function_identity']
    assert branch['all_passed'] and branch['same_exact_original_unique_branch_theorem_bound']
    assert seed['original_full_weight_function_identity']==current.outer.original_Rh.exact_join_proof()
    assert seed['original_normalized_history_ODE_binding']==current.parent.active.current.ast_binding(current.outer.original_leading.own_first_y)
    assert seed['propagated_original_histories_identified_by_linear_ODE_uniqueness']
    assert seed['same_exact_unique_branch_and_C6_implicit_derivatives']
    assert outer['initial_background']=={k:seed['actual_Rh_fields'][k] for k in current.outer.RATES}
    pc=owner.inputs[current.pulse3.NAME]
    original=current.source.FunctionTransportGraph();old=owner.inputs[current.old_frame.NAME]
    original.nodes=copy.deepcopy(old['exact_graph_nodes'])
    original.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(original.nodes)}
    transport=current.importer.GraphImport(g,original)
    frame=old['actual_Rp_native_frame_functions']['actual_identified_Rp_native_frame_C3']
    for key in current.FRAME_KEYS:
        assert [transport(current.source.FunctionRef(original,i)).node for i in frame[key]]==pc['native'][key]
    native5,evidence=current.pulse5.native_C5_functions(g,{key:current.pulse5.c3.C3Function(*(owner.ref(i) for i in ids)) for key,ids in pc['native'].items()})
    assert current.previous.encode(native5)==owner.inputs[current.pulse5.NAME]['native']
    assert q['P0'][:4]==pc['native']['P0']
    for n,i in enumerate(q['P0'][4:6],4):
        row=g.nodes[i]
        assert row['operation']=='current_original_axis_pressure_function_derivative'
        assert row['original_function']==q['P0'][0] and row['source_family']==owner.family
        assert row['ordinary_slow_Z_derivative_order']==n
        assert row['source_projection']==f'factorial({n})*P0_rows(Z).coefficients[{n}]'
    mu=Expression(g,{}).at(q['mu']);Tw=Expression(g,{}).at(q['Tw'])
    logmu=-sy.log(1000)-4*(sy.exp(40)+11)
    assert zero(mu-sy.exp(logmu)) and zero(Tw+60*logmu)
    return dict(passed=True,actual_Rp_C5_source_rows=48,
        true_Rh_full_weight_identity_and_C6_unique_branch_bound=True,
        actual_independent_fourteen_atom_P0_and_genuine_higher_pressure_provider=True,
        original_mu_and_Tw_recipes_not_interval_values=True,
        actual_complete_compact_exit_and_relative_only_correction_retained=True)


def duhamel_rows(owner,bg,bindings,Z,expected_incoming,label):
    """Check the actual graph bodies before invoking linear ODE uniqueness."""
    g=owner.g;X=Expression(g,bindings).at(bg['phase'])
    alg=Expression(g,bindings);lo=alg.at(bg['domain'][0]);hi=alg.at(bg['domain'][1])
    jac=alg.at(bg['physical_measure']) if 'physical_measure' in bg else 1/X
    offset=alg.at(bg['log_offset']) if 'log_offset' in bg else sy.log(X)
    assert zero(sy.diff(offset,X)-jac) and not offset.has(Z) and not jac.has(Z)
    for key,body in expected_incoming.items():
        for n,i in enumerate(bg['incoming'][key]):assert zero(alg.at(i)-sy.diff(body,Z,n)),(label,key,n,'incoming')
    fixtures=dict(bindings);E=sy.Function('input_E')(X)/(1+Z**2)
    V=sy.Integer(0) if bg['fields']['V'][0]==0 else sy.Function('input_V')(X)*Z
    bind(fixtures,bg['fields']['E'],E,Z);bind(fixtures,bg['fields']['V'],V,Z)
    calc=Expression(g,fixtures);checked=[]
    for key,rate in current.outer.RATES.items():
        rate=sy.Rational(str(rate));memory=calc.at(bg['memory'][key])
        assert zero(memory-sy.exp(-rate*(offset-offset.subs(X,lo))))
        for n,i in enumerate(bg['history_contributions'][key]):
            row=g.nodes[i];t=sy.Symbol(row['variable'],real=True)
            assert row['operation']=='definite_integral' and row['original_outer_Duhamel']
            assert row['ordinary_slow_Z_derivative_order']==n and row['physical_Jacobian_applied_once']
            assert sy.Rational(row['original_rate'])==rate
            assert zero(calc.at(row['lower'])-lo) and zero(calc.at(row['upper'])-X)
            et,vt=E.subs(X,t),V.subs(X,t)
            density=dict(m=vt,h=et,k=et*vt,e=vt**2-et**2/2,p=et**2/2)[key]
            wanted=jac.subs(X,t)*sy.exp(-rate*(offset-offset.subs(X,t)))*density
            assert zero(calc.at(row['integrand'])-sy.diff(wanted,Z,n)),(label,key,n,'density')
            assert owner.ref(bg['fields'][key][n])==g.add(g.mul(owner.ref(bg['incoming'][key][n]),
                owner.ref(bg['memory'][key])),owner.ref(i)),(label,key,n,'full history')
            checked.append(f'{label}_{key}_Z{n}')
    return checked


def outer_common_sources(owner):
    q,charts,f,Tw,P0=old_context(owner);data=owner.inputs[current.outer.NAME];g=owner.g;Z=q.z
    bindings={g.symbol('Z').node:Z,data['mu']:q.mu,owner.geometry['Pstar'].node:sy.exp(q.logP)}
    seed=data['seed_identity']['actual_Rh_fields'];E=sy.exp(-sy.Rational(1,2))/(1+Z**2);V=4*Z*sy.exp(-q.logP)
    Rh=dict(E=E,V=V,m=V,h=5*E/8,k=5*E*V/8,e=V**2-5*E**2/12,p=5*E**2/2)
    for key,body in Rh.items():bind(bindings,seed[key],body,Z)
    checked=[];profiles=[];inlet=Rh
    old_receipt=current.parent.linear.read(current.old_leading.RECEIPT)
    assert old_receipt['all_passed'] and all(old_receipt[k] for k in current.old_leading.GATES)
    assert old_receipt['independent_five_ODE_calculus']['all_five_own_rates_and_true_fundamental_theorem_of_calculus']
    assert old_receipt['independent_five_ODE_calculus']['total_independent_ODE_rows']==120
    assert old_receipt['independent_original_source_seams']['true_O3_x2_Rc_endpoint_rows']==20
    assert old_receipt['independent_original_native_geometry']['all_five_original_log_radius_seams']
    assert old_receipt['independent_original_native_geometry']['actual_axial_Jacobian_40_exp_40_native_once']
    for name,bg in data['background'].items():
        old=charts[name];X=q.at(old['native_coordinate']);bindings[bg['phase']]=X
        # Identify each actual sigma primitive by its full integrand/limits.
        for i,row in enumerate(g.nodes):
            if row['operation']=='definite_integral' and row.get('exact_original_sigma_integral') and row['upper']==bg['phase']:
                calc=Expression(g,bindings);t=sy.Symbol(row['variable'],real=True)
                assert calc.at(row['integrand'])==old_check.Sigma(t)
                assert calc.at(row['lower'])==0 and calc.at(row['upper'])==X
                bindings[i]=q.at(old['kernels']['J'])
        calc=Expression(g,bindings)
        for key in ('E','V'):
            body=q.at(old[key][0])
            for n,i in enumerate(bg['fields'][key]):
                assert zero(calc.at(i)-sy.diff(body,Z,n)),(name,key,n)
                profiles.append(f'{name}_{key}_Z{n}')
            bind(bindings,bg['fields'][key],body,Z)
        incoming={key:inlet[key] for key in current.outer.RATES}
        checked.extend(duhamel_rows(owner,bg,bindings,Z,incoming,name))
        # The accepted old bodies obey exactly these same five ODEs and
        # seams. Equal current profiles and inlets therefore identify the
        # defining histories, including all genuine differentiated integrals.
        lo=Expression(g,bindings).at(bg['domain'][0]);hi=Expression(g,bindings).at(bg['domain'][1])
        next_inlet={}
        for key,ids in old['histories'].items():
            body=q.at(ids[0]);assert zero(endpoint(q,body,X,lo)-incoming[key]),(name,key,'old inlet')
            bind(bindings,bg['fields'][key],body,Z)
            next_inlet[key]=endpoint(q,body,X,hi)
        next_inlet['E']=endpoint(q,q.at(old['E'][0]),X,hi)
        next_inlet['V']=endpoint(q,q.at(old['V'][0]),X,hi)
        inlet=next_inlet
        print('ACTUAL_RP_OUTER_COMMON_'+name.upper(),flush=True)
    for key,ids in owner.inputs[current.compact.NAME]['incoming_background'].items():
        calc=Expression(g,bindings)
        for n,i in enumerate(ids):assert zero(calc.at(i)-sy.diff(inlet[key],Z,n)),(key,n,'actual Rc inlet')
    return dict(passed=True,actual_profile_C6_rows=profiles,actual_current_Duhamel_C6_rows=checked,
        source_equalities_follow_same_actual_profiles_density_rates_inlets_and_geometry=True,
        O2_axial_exponential_Jacobian_retained_once=True,
        accepted_old_smooth_function_ODE_and_seam_theorem_consumed=True),dict(q=q,charts=charts,f=f,Tw=Tw,P0=P0,bindings=bindings,inlet=inlet)


def compact_quiet_common_sources(owner,state):
    g,ref=owner.g,owner.ref;q,Z=state['q'],state['q'].z;f,Tw,P0=state['f'],state['Tw'],state['P0']
    c=owner.inputs[current.compact.NAME];b=dict(state['bindings']);X=sy.Symbol('compact_x',positive=True)
    b[c['background']['phase']]=X;b[c['mu']]=q.mu;b[owner.geometry['Tw'].node]=Tw
    bind(b,c['P0'],P0,Z)
    calc=Expression(g,b)
    for n,i in enumerate(c['A']):assert zero(calc.at(i)-sy.diff(state['inlet']['E'],Z,n))
    power=f['original_power_at_Tw_C3'];Ec=q.at(power['profiles_y0_y1_y2_y3_y4_C3'][0]['E'][0])
    Hc={key:q.at(ids[0]) for key,ids in power['complete_histories_y0_y1_y2_y3_y4_C3'][0].items()}
    assert zero(Ec.subs(Tw,2)-state['inlet']['E'])
    for key in current.outer.RATES:assert zero(Hc[key].subs(Tw,2)-state['inlet'][key]),key
    expectedE=Ec.subs(Tw,2+sy.log(X));expectedH={key:body.subs(Tw,2+sy.log(X)) for key,body in Hc.items()}
    for key,body in (('E',expectedE),('V',sy.Integer(0))):
        for n,i in enumerate(c['background']['fields'][key]):assert zero(calc.at(i)-sy.diff(body,Z,n)),(key,n,'compact profile')
        bind(b,c['background']['fields'][key],body,Z)
    checked=duhamel_rows(owner,c['background'],b,Z,{key:state['inlet'][key] for key in current.outer.RATES},'compact_background')
    for key,body in expectedH.items():bind(b,c['background']['fields'][key],body,Z)
    import lei_ren_part1_paper_compliant_current_original_Rp_native_C5_terminal_projection as terminal_projection
    actual_terminal_proof=terminal_projection.prove(owner,Expression,zero)
    assert len(actual_terminal_proof['actual_integral_to_residual_C5_rows'])==30
    compact_receipt=current.parent.linear.read(current.compact.RECEIPT)
    assert compact_receipt['checks']['actual_profiles']['actual_original_raw_beta_end_flat_rows']==18
    assert compact_receipt['checks']['actual_inlet_and_endpoint']['terminal_integral_identity_uses_pinned_original_exact_bump_theorem']
    for key,i in c['terminal_certificates'].items():
        row=g.nodes[i]
        assert row['actual_endpoint_rows']==c['terminal']['correction'][key]
        assert row['equivalent_residual_rows']==c['terminal_rhs'][key]
        assert row['source_system']==c['system'] and row['relative_correction_only']
        assert row['absolute_leading_histories_and_pressure_not_zeroed']
    flat={c['background']['phase']:sy.Integer(2),c['N']:sy.Symbol('N',positive=True),c['mu']:q.mu,c['W']['normal']:sy.Integer(1),g.symbol('Z').node:Z}
    bind(flat,c['A'],state['inlet']['E'],Z)
    for j,ids in enumerate(c['h']):bind(flat,ids,sy.Function('actual_control_'+str(j))(Z),Z)
    flat_alg=compact_check.Algebra(g,flat)
    for triplet in c['profiles']['bumps']:
        for i in triplet:assert zero(flat_alg.at(i)),('actual beta at x2',i)
    for key in ('delta_E','delta_V'):
        for i in c['profiles']['radial_y_rows'][0][key]:assert zero(flat_alg.at(i)),key
    # Consume the certified relative zero, while retaining the full leading
    # absolute histories. The terminal modified field decomposition is checked.
    for key,body in {**expectedH,'E':expectedE,'V':sy.Integer(0),'P0':P0,'pressure':P0+expectedH['p']}.items():
        bind(b,c['background']['fields'][key],body,Z)
    for key,ids in c['terminal']['background_fields'].items():
        body=P0+expectedH['p'] if key=='pressure' else P0 if key=='P0' else expectedE if key=='E' else sy.Integer(0) if key=='V' else expectedH[key]
        for n,i in enumerate(ids):assert zero(Expression(g,b).at(i)-sy.diff(body.subs(X,2),Z,n)),(key,n,'background terminal')
    terminal={key:body.subs(Tw,2+sy.log(2)) for key,body in Hc.items()}
    terminal.update(E=Ec.subs(Tw,2+sy.log(2)),V=sy.Integer(0),P0=P0)
    terminal['pressure']=P0+terminal['p']
    # Check the actual source additions and actual terminal substitutions.
    # Linearity then combines the independently proved background and the
    # actual integral-to-branch terminal zero. No correction is bound to zero.
    terminal_rows=[]
    for key,ids in c['terminal']['modified_fields'].items():
        for n,i in enumerate(ids):
            row=g.nodes[i]
            assert row['operation']=='function_substitution' and row['variable']==c['background']['phase']
            assert row['value']==g.constant(2).node and row['Z_independent_substitution']
            assert row['expression']==c['modified_fields'][key][n]
            if key in current.outer.RATES:
                assert row['expression']==g.add(ref(c['background']['fields'][key][n]),ref(c['history'][key][n])).node
            elif key=='pressure':
                assert row['expression']==g.add(ref(c['P0'][n]),ref(c['modified_fields']['p'][n])).node
            elif key=='P0':assert row['expression']==c['P0'][n]
            elif key in ('E','V'):
                assert zero(flat_alg.at(i)-sy.diff(terminal[key],Z,n)),(key,n,'actual terminal flat profile')
            else:raise AssertionError('Unexpected actual compact terminal field')
            terminal_rows.append(f'{key}_Z{n}')
    quiet=owner.inputs[current.quiet.NAME];S=sy.Symbol(g.nodes[quiet['s']]['name'],real=True)
    native={key:q.at(f['actual_identified_Rp_native_frame_C3'][key][0]) for key in current.FRAME_KEYS}
    qb={g.symbol('Z').node:Z,quiet['s']:S,quiet['mu']:q.mu,quiet['Tw']:Tw}
    bind(qb,quiet['P0'],P0,Z)
    for key,body in terminal.items():bind(qb,quiet['incoming']['modified_fields'][key],body,Z)
    qb[quiet['s']]=Expression(g,qb).at(quiet['length'])
    calc=Expression(g,qb);rows=[]
    for key,jet in owner.current_frame.items():
        for n,ref in enumerate(current.previous.rows(jet)):
            assert zero(calc.at(ref)-sy.diff(native[key],Z,n)),(key,n,'actual current Rp native')
            rows.append(key+'_Z'+str(n))
    return dict(passed=True,actual_compact_background_Duhamel_C6_rows=checked,
        actual_relative_terminal_C5_rows=30,actual_terminal_velocity_flatness=True,
        actual_modified_terminal_decomposition_rows=terminal_rows,actual_eight_Rp_frame_C5_rows=rows,
        actual_terminal_integral_and_unique_branch_source_proof=actual_terminal_proof,
        complete_absolute_histories_and_independent_pressure_not_zeroed=True,
        same_finite_mu_power_semigroup_and_original_quiet_length=True),dict(native=native,Z=Z,q=q,Tw=Tw,P0=P0)


def callback_and_consumer(owner,proofs):
    assert proofs['source_inputs']['passed'] and proofs['outer_common_sources']['passed'] and proofs['compact_quiet']['passed']
    assert current.previous.encode(owner.incoming_source_rows())==owner.inputs[current.quiet.NAME]['Rp_frame']
    old=owner.inputs[current.pulse5.NAME];rewrite=current.PulseSourceSubstitution(owner.g,owner.input_substitution)
    count=0
    for group,new in ((old['selected_pulse'],owner.selected_pulse),(old['pressure'],owner.pressure)):
        for key,ids in group.items():
            for i,ref in zip(ids,current.previous.rows(new[key])):assert rewrite(owner.ref(i))==ref;count+=1
    for chart,fields in old['leading_fields'].items():
        for key,ids in fields.items():
            for i,ref in zip(ids,current.previous.rows(owner.leading_fields[chart][key])):assert rewrite(owner.ref(i))==ref;count+=1
    assert 0 not in owner.input_substitution and 1 not in owner.input_substitution
    assert owner.current_frame['Mp'].value!=owner.g.zero and owner.current_frame['P0'].value!=owner.g.zero
    Z=sy.Symbol('Z',real=True);mu,P=sy.symbols('mu Pstar',positive=True)
    u,axis,Mp=(sy.Function(name)(Z) for name in ('actual_native_u','actual_axis_P0','actual_Mp'))
    bindings={owner.g.symbol('Z').node:Z,owner.geometry['mu'].node:mu,owner.geometry['Pstar'].node:P}
    for key,body in (('u',u),('P0',axis),('Mp',Mp)):
        bind(bindings,[ref.node for ref in current.previous.rows(owner.current_frame[key])],body,Z)
    absolute=[]
    for chart,packet in owner.pulse_coordinates.items():
        x=sy.Symbol('query_'+chart,real=True);local=dict(bindings);local[packet['query'].node]=x
        y=x+13/mu if chart=='end' else x/mu
        expected=P**2*(axis+Mp+u**2*(1-sy.exp(-(1+2*mu)*y))/(2*(1+2*mu)))
        calc=Expression(owner.g,local)
        assert zero(calc.at(packet['y'])-y)
        for n,ref in enumerate(current.previous.rows(owner.absolute_pressure[chart])):
            assert zero(calc.at(ref)-sy.diff(expected,Z,n)),(chart,n,'absolute pressure')
            absolute.append(chart+'_Z'+str(n))
        assert zero(sy.diff(expected,x)-P**2*u**2*sy.exp(-(1+2*mu)*y)*sy.diff(y,x)/2)
    return dict(passed=True,actual_current_eight_function_source_callback_rows=48,
        actual_selected_coefficient_field_and_pressure_source_rows=count,
        actual_absolute_background_pressure_C3_rows=absolute,
        actual_current_Rp_coordinate_and_absolute_pressure_history_preserved=True,
        exact_expression_substitution_not_a_second_solution_or_numerical_center=True,
        unchanged_original_Gamma_phase_scales_and_once_physical_normalization=True,
        interval_numeric_callback_bounds_other_joins_and_global_n1_not_claimed=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['Rp_native_C5_pulse_source_bindings']==owner.snapshot
    checks={'source_inputs':source_inputs(owner)};print('RP_NATIVE_SOURCE_INPUTS_CHECKED',flush=True)
    checks['outer_common_sources'],state=outer_common_sources(owner)
    checks['compact_quiet'],native=compact_quiet_common_sources(owner,state)
    print('RP_NATIVE_ACTUAL_EIGHT_FRAME_IDENTITIES_CHECKED',flush=True)
    checks['callback_and_consumer']=callback_and_consumer(owner,checks)
    owner.assert_graph()
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ACTUAL_CURRENT_RP_NATIVE_C5_PULSE_SOURCE_CHECKED',flush=True);return rec


if __name__=='__main__':raise SystemExit('run(owner) checks the saved actual source graph')

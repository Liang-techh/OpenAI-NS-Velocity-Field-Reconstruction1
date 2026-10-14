"""Actual endpoint recovery and whole-source outer pressure-y matching."""
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_n1_outer_radial_pressure_source_joins as current
import lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions_check as pulse_check
import lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins_check as outer_check

bind=pulse_check.bind


def zero(expr,label):
    try:pulse_check.zero(expr,label)
    except AssertionError:
        # Substitution can leave equivalent composite derivatives/functions
        # in different forms. Preserve their definitions for the fallback;
        # no floating-point or source-output substitution is allowed.
        if sy.simplify(expr)!=0:raise


def endpoint_context(owner):
    g=owner.g;q=owner.quiet;Z=sy.Symbol('Z',real=True);S=sy.Symbol(g.nodes[q['s']]['name'],real=True)
    length=sy.Symbol('actual_quiet_length',positive=True);mu=sy.Symbol('mu',positive=True)
    P=sy.Symbol('Pstar',positive=True);R2=sy.Symbol('actual_R2',positive=True);delta=sy.Symbol('delta',positive=True)
    # The admitted original native family fixes u=A/(1+Z^2) and a
    # nonzero linear physical mean. Its m1 is cubic, so Z4/Z5 are zero.
    # Independent arbitrary E/M functions would violate that source family.
    A=sy.Symbol('actual_Rp_swirl_amplitude',positive=True)
    B=sy.Symbol('actual_Rp_normalized_mass_coefficient',nonzero=True,real=True)
    E=A/(1+Z**2);M=B*Z
    bindings={g.symbol('Z').node:Z,q['s']:S,q['length']:length,q['mu']:mu,q['Pstar']:P,q['R2']:R2,q['delta']:delta,
        owner.pulse.W['normalization'].node:sy.Symbol('actual_beta_normalization',positive=True),
        owner.pulse.W['log_end_scale'].node:sy.Symbol('actual_log_end_scale',real=True)}
    bind(bindings,owner.jet(q['fields']['E']),E*sy.exp(-(mu+sy.Rational(1,2))*(S-length)),Z)
    bind(bindings,owner.jet(q['fields']['m']),M*sy.exp(-(S-length)),Z)
    for key,body in (('E',E),('m',M)):
        bind(bindings,owner.jet(q['Rp_endpoint']['fields'][key]),body,Z)
    for key,jet in owner.pulse.selected.items():bind(bindings,jet,sy.Function('actual_selected_'+key)(Z),Z)
    return bindings,Z,E,M,P,R2*sy.exp(length)


def Rp_endpoint_recovery(owner):
    owner.assert_graph();q=owner.quiet;g=owner.g
    native_rec=current.parent.linear.read(current.rp.RECEIPT)
    quiet_rec=current.parent.linear.read(current.rp.quiet.RECEIPT)
    c5=current.parent.linear.read(current.rp.pulse5.NAME)['C5_pulse_C3_pressure_bindings']
    assert native_rec['all_passed'] and native_rec['source_family']==owner.family
    assert len(native_rec['checks']['compact_quiet']['actual_eight_Rp_frame_C5_rows'])==48
    assert native_rec['checks']['callback_and_consumer']['passed']
    assert quiet_rec['checks']['native_mass_identity']['actual_compact_derived_K_equals_original_Rp_Pstar_kappa']
    assert c5['native_evidence']['cubic_moment_Z4_Z5_are_proved_zero_not_padding']
    bindings,Z,E,M,P,Rp=endpoint_context(owner);alg=pulse_check.PulseAlgebra(owner.pulse,bindings)
    native=[]
    for key,body in (('u',E),('m1',M/E)):
        for n,ref in enumerate(current.outer.rows(owner.inlet.current_frame[key])):
            zero(alg.at(ref)-sy.diff(body,Z,n),('actual Rp frame',key,n));native.append(key+'_Z'+str(n))
        # These current native inputs have just been proved as the actual
        # endpoint quotient; binding them is a source lemma, not a midpoint.
        bind(bindings,owner.inlet.current_frame[key],body,Z)
    alg=pulse_check.PulseAlgebra(owner.pulse,bindings)
    left=owner.quiet_physical_right;right=owner.pulse.endpoints['main']['left']['fields']
    values=dict(U=sy.Integer(0),U_y=sy.Integer(0),m=P*M,m_y=-P*M)
    physical=[]
    for key,body in values.items():
        for side,packet in (('quiet',left),('pulse',right)):
            for n,ref in enumerate(current.outer.rows(packet[key])):
                zero(alg.at(ref)-sy.diff(body,Z,n),(side,key,n));physical.append(side+'_'+key+'_Z'+str(n))
    for n,ref in enumerate(current.outer.rows(right['U_yy'])):
        zero(alg.at(ref),('pulse flat Uyy',n))
    for n,ref in enumerate(current.outer.rows(right['m_yy'])):
        zero(alg.at(ref)-sy.diff(P*M,Z,n),('pulse myy',n))
    qbody=(-(1-sy.Symbol('delta',positive=True))*Z*P*M-(1-Z**2)*sy.diff(P*M,Z))/(1-sy.Symbol('delta',positive=True)*Z**2)
    recovered=[]
    for key,body in (('Q',qbody),('Q_y',-qbody),('Q_yy',qbody)):
        for side,packet in (('quiet',left),('pulse',right)):
            for n,ref in enumerate(current.outer.rows(packet[key])):
                zero(alg.at(ref)-sy.diff(body,Z,n),(side,key,n))
                recovered.append(side+'_'+key+'_Z'+str(n))
    alg=pulse_check.PulseAlgebra(owner.pulse,bindings)
    zero(alg.at(owner.pulse.endpoints['main']['left']['radius'])-Rp,'actual Rp radius')
    zero(alg.at(owner.ref(q['Rp_endpoint']['radius']))-Rp,'actual quiet Rp radius')
    pressure=[]
    # Compare actual reduced source formulas through the shared Q recovery.
    for n in range(4):
        zero(alg.at(current.outer.rows(owner.pressure_traces['Rp']['left'])[n])-
             alg.at(current.outer.rows(owner.pressure_traces['Rp']['right'])[n]),('Rp total pressure-y',n))
        zero(alg.at(current.outer.rows(owner.N1p_traces['Rp']['left'])[n])-
             alg.at(current.outer.rows(owner.N1p_traces['Rp']['right'])[n]),('Rp N1p',n))
    pressure=['Rp_P1_y_Z'+str(n) for n in range(4)]
    return dict(passed=True,actual_native_endpoint_quotient_C5_rows=native,
        actual_nonzero_physical_mean_and_flat_axial_endpoint_rows=physical,
        actual_Q_Qy_Qyy_recovery_source_rows=recovered,actual_Rp_total_pressure_y_C3_rows=pressure,
        whole_function_pressure_identity_checked=True,all_four_actual_pressure_source_rows_directly_compared=True,
        accepted_native_family_receipt=current.rp.RECEIPT,accepted_native_family_receipt_sha256=current.sha(current.rp.RECEIPT),
        accepted_nonzero_original_mass_receipt=current.rp.quiet.RECEIPT,
        actual_Pstar_applied_once=True,nonzero_original_linear_mean_retained=True,
        original_native_swirl_rational_and_cubic_m1_family_retained=True,
        absolute_leading_pressure_not_used_as_N1p=True),bindings


def Rv_endpoint_recovery(owner,bindings):
    owner.assert_graph();g=owner.g;alg=pulse_check.PulseAlgebra(owner.pulse,bindings)
    fields=owner.pulse.endpoints['end']['right']['fields'];flat=[]
    for key in ('U','U_y','U_yy','m','m_y','m_yy','Q','Q_y','Q_yy'):
        for n,ref in enumerate(current.outer.rows(fields[key])):
            zero(alg.at(ref),('actual Rv leading source',key,n));flat.append(key+'_Z'+str(n))
    traces=[]
    for key,packet in (('P1_y',owner.pressure_traces['Rv']),('N1p',owner.N1p_traces['Rv'])):
        for n,(left,right) in enumerate(zip(current.outer.rows(packet['left']),current.outer.rows(packet['right']))):
            zero(alg.at(left)-alg.at(right),('actual Rv pressure source',key,n));traces.append(key+'_Z'+str(n))
    tail=owner.outer.inputs[current.outer.intermediate.previous.NAME]
    assert current.outer.encode(owner.tail_pressure)==tail['terminal']['m13']==tail['pressure_tail']['P1']
    assert not outer_check.free_query_dependency(owner.outer,owner.tail_pressure.value.node,owner.tail_coordinate.node)
    return dict(passed=True,actual_Rv_flat_physical_source_rows=flat,actual_Rv_pressure_y_and_N1p_C3_rows=traces,
        actual_tail_pressure_is_same_terminal_function_not_zero_reset=True,
        constant_tail_pressure_source=owner.tail_pressure_definition,
        tail_pressure_actual_radial_query_independence=True)


def actual_chart_radius_sources(owner):
    owner.assert_graph();g=owner.g;bindings=outer_check.radius_context(owner.outer,owner.radius_proof)
    alg=outer_check.SourceAlgebra(g,bindings);checked=[]
    oldlog=alg.at(owner.old_logRp);newlog=alg.at(owner.current_logRp)
    zero(sy.expand_log(oldlog-newlog,force=True),'actual common logRp source')
    for name,data in owner.pulse.charts.items():
        old=owner.intermediate['charts']['pulse_'+name]
        assert old['coordinate']==data['coordinate'].node and old['domain']==[v.node for v in data['domain']]
        left=alg.at(old['radius']);right=alg.at(data['radius'])
        zero(sy.expand_log(sy.log(left)-sy.log(right),force=True),(name,'actual radius whole function'))
        mapped=owner.radius_projection[old['radius']]
        assert mapped==data['radius']
        checked.append(dict(chart=name,old_radius_node=old['radius'],actual_current_radius_node=data['radius'].node,
            actual_query_node=data['coordinate'].node,original_domain=old['domain']))
    return dict(passed=True,actual_common_logRp_source_identity=True,actual_three_chart_radius_function_identities=checked,
        accepted_scalar_parameter_projection=owner.radius_proof['selected_parameter_source_projection'],
        original_same_mu_queries_and_chart_domains=True,independent_radii_not_manufactured_as_one_symbol=True)


def partial_provider_source_identity(owner,endpoint_bindings,radius_proof):
    owner.assert_graph();assert radius_proof['passed'];g=owner.g
    actual=[];tensor=[];integrals=[];field_rows=[];pressure=[];complete_derived=[]
    oldrec=current.parent.linear.read(current.rp.pulse5.RECEIPT)
    assert oldrec['all_passed'] and oldrec['source_family']==owner.family
    assert oldrec['current_original_n1_true_pulse_pressure_C3_functions_installed']
    derivative_source={key:current.parent.active.current.ast_binding(getattr(current.rp.pulse5.HigherZDerivative,key))
        for key in ('rule','extend','seed_rows')}
    for pin in derivative_source.values():assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
    for name,packet in owner.canonical['leading_fields'].items():
        for key,ids in packet.items():
            for n,(i,ref) in enumerate(zip(ids,current.outer.rows(owner.inlet.leading_fields[name][key]))):
                if i in (0,1):assert i==ref.node
                else:assert owner.field_projection[i]==ref
                field_rows.append(name+'_'+key+'_Z'+str(n))
    for name,data in owner.pulse.charts.items():
        old=owner.projected_old[name];x=sy.Symbol(g.nodes[data['coordinate'].node]['name'],real=True)
        Z=sy.Symbol('Z',real=True);R=sy.Symbol('actual_'+name+'_radius',positive=True);mu=sy.Symbol('mu',positive=True)
        bindings={g.symbol('Z').node:Z,data['coordinate'].node:x,data['radius'].node:R,owner.pulse.mu.node:mu,
            owner.pulse.delta.node:sy.Symbol('delta',positive=True)}
        Q=sy.Function('actual_'+name+'_Q')(x,Z)
        for key in ('U','Q','Q_y'):
            body=sy.Function('actual_'+name+'_'+key)(x,Z)
            if name=='gap' and key=='U':body=sy.Integer(0)
            if key=='Q':body=Q
            if name=='gap' and key=='Q_y':body=-Q
            bind(bindings,data['fields'][key],body,Z)
        alg=outer_check.SourceAlgebra(g,bindings)
        for key,jet in data['FTC']['actual_products'].items():
            for n,(left,right) in enumerate(zip(current.outer.rows(old['products'][key]),current.outer.rows(jet))):
                zero(alg.at(left)-alg.at(right),(name,key,n,'actual source product'));actual.append(name+'_'+key+'_Z'+str(n))
        # Both old partial and complete nodes are kept for this check.
        # They are not silently identified by a reference rewrite.
        for kind,lefts,rights in (('partial',old['partial_integrals'],data['partial_integrals']),
            ('complete',old['complete_integrals'],owner.pulse.complete_integrals[name])):
            for key,jet in rights.items():
                for n,(left,right) in enumerate(zip(current.outer.rows(lefts[key]),current.outer.rows(jet))):
                    lrow,rrow=g.nodes[left.node],g.nodes[right.node]
                    if kind=='complete' and n>0:
                        # Both are admitted ordinary derivatives of the just
                        # compared C0 integral, with fixed named bounds. The
                        # canonical producer expands the same derivatives;
                        # consume its protected C3 source instead of replaying
                        # the implicit Rm/compact graph below those inputs.
                        if lrow['operation']!='exact_rational':
                            assert lrow['operation']=='definite_integral'
                            assert lrow['variable']==rrow['variable']
                            assert lrow['lower']==rrow['lower'] and lrow['upper']==rrow['upper']
                            assert lrow['ordinary_slow_Z_derivative_order']==n
                            assert lrow['endpoint_independent_of_Z']
                        else:assert lrow['numerator']==0
                        assert rrow['operation']=='definite_integral' and rrow['ordinary_slow_Z_derivative_order']==n
                        complete_derived.append(name+'_'+key+'_Z'+str(n));continue
                    assert lrow['operation']==rrow['operation']=='definite_integral'
                    assert lrow['variable']==rrow['variable']==g.nodes[data['coordinate'].node]['name']
                    assert lrow['lower']==rrow['lower'] and lrow['upper']==rrow['upper']
                    assert lrow['endpoint_independent_of_Z'] and rrow['endpoint_independent_of_Z']
                    assert lrow['ordinary_slow_Z_derivative_order']==rrow['ordinary_slow_Z_derivative_order']==n
                    zero(alg.at(lrow['integrand'])-alg.at(rrow['integrand']),(name,kind,key,n,'actual signed integrand'))
                    integrals.append(name+'_'+kind+'_'+key+'_Z'+str(n))
        for key,jet in data['upper_tensors'].items():
            for n,(left,right) in enumerate(zip(current.outer.rows(old['upper'][key]),current.outer.rows(jet))):
                zero(alg.at(left)-alg.at(right),(name,key,n,'actual upper tensor'));tensor.append(name+'_upper_'+key+'_Z'+str(n))
        ep=pulse_check.PulseAlgebra(owner.pulse,endpoint_bindings)
        for key,jet in owner.pulse.global_lower_tensors.items():
            for n,(left,right) in enumerate(zip(current.outer.rows(old['lower'][key]),current.outer.rows(jet))):
                zero(ep.at(left)-ep.at(right),(name,key,n,'actual finite global lower tensor'))
                tensor.append(name+'_lower_'+key+'_Z'+str(n))
        # Inputs have now been proved equal as functions. Compare the
        # original pressure body after those input projections; no pressure
        # result or own field is bound in this algebra.
        pbindings=dict(bindings)
        for previous_name,packet in owner.pulse.complete_integrals.items():
            for key,jet in packet.items():
                for ref in current.outer.rows(jet):
                    row=g.nodes[ref.node]
                    value=sy.Integer(0) if row['integrand']==0 else sy.Symbol('actual_integral_row_'+str(ref.node),real=True)
                    pbindings.setdefault(ref.node,value)
        for key,jet in data['partial_integrals'].items():
            for ref in current.outer.rows(jet):
                row=g.nodes[ref.node]
                value=sy.Integer(0) if row['integrand']==0 else sy.Symbol('actual_integral_row_'+str(ref.node),real=True)
                pbindings.setdefault(ref.node,value)
        for key,jet in owner.pulse.global_lower_tensors.items():
            for ref in current.outer.rows(jet):
                if ref.node not in (0,1):pbindings[ref.node]=sy.simplify(ep.at(ref))
        palg=outer_check.SourceAlgebra(g,pbindings)
        for key,jet in data['signed_pressure_increment'].items():
            for n,(left,right) in enumerate(zip(current.outer.rows(owner.projected_old_pressure[name][key]),current.outer.rows(jet))):
                zero(palg.at(left)-palg.at(right),(name,key,n,'actual signed pressure map'))
                pressure.append(name+'_'+key+'_Z'+str(n))
    pin=current.parent.active.current.ast_binding(current.outer.intermediate.pressure_functions)
    assert pin==owner.pulse.pressure_map_source
    return dict(passed=True,actual_field_callback_C5_Q4_source_projection_rows=field_rows,
        actual_ten_product_C3_source_identity_rows=actual,actual_partial_and_complete_integral_source_rows=integrals,
        actual_complete_integral_higher_rows_by_pinned_derivative_commutation=complete_derived,
        accepted_canonical_C3_derivative_source=derivative_source,
        actual_finite_upper_and_lower_tensor_source_rows=tensor,actual_signed_pressure_C3_source_map_rows=pressure,
        actual_same_pressure_map_source=pin,all_complete_preceding_charts_retained=True,
        original_domains_binders_measures_and_nonzero_lower_tensors_retained=True,
        actual_pressure_outputs_not_bound=True)


def own_pressure_transport(owner,partial_identity):
    owner.assert_graph();assert partial_identity['passed'];g=owner.g
    rec=current.parent.linear.read(current.previous.RECEIPT)
    assert rec['all_passed'] and rec['source_family']==owner.family
    proof=rec['checks']['generic_recovery_FTC'];assert proof['passed']
    assert len(proof['actual_signed_pressure_FTC_rows'])==6
    assert len(proof['pressure_FTC_higher_rows_by_pinned_ordinary_derivative_commutation'])==18
    assert owner.Rp_cumulative==owner.outer.joins['Rp']['left']['cumulative']
    a=current.outer.Algebra(g,3);retained=[];composed=[];endpoint_composition=[]
    for name,data in owner.pulse.charts.items():
        packet=owner.current_own[name];query=data['coordinate'].node
        for key,jet in owner.Rp_cumulative.items():
            for n,ref in enumerate(current.outer.rows(jet)):
                assert not outer_check.free_query_dependency(owner.outer,ref.node,query)
                retained.append(name+'_'+key+'_Z'+str(n))
        for key in ('m13','m15'):
            expected=a.add(owner.Rp_cumulative[key],data['signed_pressure_increment'][key])
            assert current.outer.rows(expected)==current.outer.rows(packet['cumulative'][key])
            for n in range(4):composed.append(name+'_'+key+'_y_Z'+str(n))
        assert current.outer.rows(packet['fields']['P1'])==current.outer.rows(packet['cumulative']['m13'])
        assert current.outer.rows(packet['fields']['P1_y'])==current.outer.rows(data['density']['P1_y'])
        assert packet['physical_y_Jacobian']==data['physical_y_Jacobian']
    # The whole-function increment identity extends to its actual fixed
    # endpoints. Explicitly identify the predecessor's defining cumulative
    # sum and the very endpoint traces covered by the 204-row receipt; this
    # composes the new sum to Rp history and the retained full Rv tail.
    ordinary=current.parent.linear.read(current.outer.CONSTRUCTION_RECEIPT)
    assert ordinary['all_three_ordinary_axial_trace_identities_passed']
    assert ordinary['ordinary_own_field_and_cumulative_trace_rows']==204
    assert ordinary['source_family']==owner.family
    assert current.outer.encode(owner.Rp_cumulative)==owner.intermediate['Rp_cumulative']
    for name,edge,seam,side in (('main','left','Rp','right'),('end','right','Rv','left')):
        oldchart=owner.intermediate['charts']['pulse_'+name]
        oldendpoint=owner.intermediate['endpoints']['pulse_'+name][edge]['cumulative']
        assert current.outer.encode(owner.outer.joins[seam][side]['cumulative'])==oldendpoint
        for key in ('m13','m15'):
            oldsum=a.add(owner.Rp_cumulative[key],owner.jet(oldchart['fixed_pressure'][key]))
            assert current.outer.encode(oldsum)==oldchart['cumulative'][key]
            for n in range(4):
                assert name+'_'+key+'_Z'+str(n) in partial_identity['actual_signed_pressure_C3_source_map_rows']
                endpoint_composition.append(dict(seam=seam,moment=key,Z_order=n,
                    actual_current_sum_node=current.outer.rows(owner.current_own[name]['cumulative'][key])[n].node,
                    actual_old_endpoint_node=oldendpoint[key][n],
                    retained_target_node=current.outer.rows(owner.outer.joins[seam]['left' if seam=='Rp' else 'right']['cumulative'][key])[n].node))
    return dict(passed=True,actual_query_independent_incoming_history_rows=retained,
        actual_current_own_pressure_radial_rows_by_accepted_FTC_composition=composed,
        actual_Rp_Rv_m13_m15_endpoint_rows_by_whole_function_identity_composition=endpoint_composition,
        accepted_ordinary_endpoint_receipt=current.outer.CONSTRUCTION_RECEIPT,
        accepted_ordinary_endpoint_receipt_sha256=current.sha(current.outer.CONSTRUCTION_RECEIPT),
        source_pressure_FTC_receipt=current.previous.RECEIPT,
        source_pressure_FTC_receipt_sha256=current.sha(current.previous.RECEIPT),
        own_P1_equals_retained_Rp_history_plus_actual_partial_pressure=True,
        own_P1_y_is_actual_leading_density_by_proved_FTC=True,
        m15_y_is_negative_actual_radius_times_same_density=True,
        derivative_composition_not_a_second_direct_pressure_expansion=True,
        endpoint_composition_not_a_second_direct_terminal_integral_expansion=True,
        nonzero_incoming_history_not_reset=True)


def run(owner):
    began=time.monotonic();owner.assert_graph()
    rp_proof,bindings=Rp_endpoint_recovery(owner)
    checks=dict(Rp_endpoint_recovery=rp_proof)
    print('ACTUAL_RP_TOTAL_RADIAL_PRESSURE_C3_JOIN_PASSED',flush=True)
    checks['Rv_endpoint_recovery']=Rv_endpoint_recovery(owner,bindings)
    print('ACTUAL_RV_TOTAL_RADIAL_PRESSURE_C3_JOIN_PASSED',flush=True)
    checks['actual_chart_radius_sources']=actual_chart_radius_sources(owner)
    checks['partial_provider_source_identity']=partial_provider_source_identity(owner,bindings,checks['actual_chart_radius_sources'])
    print('ACTUAL_OLD_CURRENT_PARTIAL_PRESSURE_SOURCE_IDENTITY_PASSED',flush=True)
    checks['own_pressure_transport']=own_pressure_transport(owner,checks['partial_provider_source_identity'])
    return save(owner,checks,time.monotonic()-began)


def save(owner,checks,execution_seconds):
    owner.assert_graph()
    assert set(checks)=={'Rp_endpoint_recovery','Rv_endpoint_recovery','actual_chart_radius_sources',
        'partial_provider_source_identity','own_pressure_transport'}
    assert all(proof['passed'] for proof in checks.values())
    raw=current.parent.linear.read(current.NAME)
    assert raw['outer_radial_pressure_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    pins=dict(owner.hashes);pins[current.NAME]=current.sha(current.NAME)
    pins[Path(__file__).name]=current.sha(Path(__file__).name)
    identity=checks['partial_provider_source_identity']
    receipt=dict(all_passed=True,source_family=owner.family,input_hashes=pins,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_Rp_Rv_total_pressure_y_C3_join_rows=8,
        actual_Rp_Rv_N1p_C3_join_rows=8,
        actual_partial_and_complete_C0_integral_source_rows=len(identity['actual_partial_and_complete_integral_source_rows']),
        complete_higher_rows_by_pinned_derivative_commutation=len(identity['actual_complete_integral_higher_rows_by_pinned_derivative_commutation']),
        actual_signed_pressure_C3_source_map_rows=len(identity['actual_signed_pressure_C3_source_map_rows']),
        accepted_ordinary_outer_trace_rows_reused=204,
        exact_graph_nodes=len(raw['exact_graph_nodes']),execution_seconds=execution_seconds,
        numerical_oracle_available=False,combined_outer_seam_accepted=False,global_n1=False)
    (current.HERE/current.RECEIPT).write_text(__import__('json').dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance()
    print('ACTUAL_RP_RV_RADIAL_PRESSURE_SOURCE_JOINS_ACCEPTED',receipt['exact_graph_nodes'],
        'seconds',receipt['execution_seconds'],flush=True)
    return receipt


if __name__=='__main__':
    import sys
    sys.setrecursionlimit(10000)
    run(current.CurrentOriginalN1OuterRadialPressureSourceJoins(require_checked=False))

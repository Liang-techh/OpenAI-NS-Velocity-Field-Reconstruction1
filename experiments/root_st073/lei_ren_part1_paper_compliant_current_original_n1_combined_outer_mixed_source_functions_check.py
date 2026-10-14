"""Compose actual finite own outer seams from pinned source theorems.

Every contracted row retains its actual left/right source node. Previously
proved primitive and pressure identities are consumed as theorems; the
new cumulative calculus and state recovery are checked here.
"""
from pathlib import Path
import json
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_n1_combined_outer_mixed_source_functions as current
import lei_ren_part1_paper_compliant_current_original_compact_higher_radial_pressure_C3_functions_check as compact_check
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions as support_source

outer_check=compact_check.outer_check
zero,bind=compact_check.zero,compact_check.bind
Algebra,rows,encode=current.outer.Algebra,current.outer.rows,current.outer.encode
PRESSURE=('P1','P1_y','P1_yy','P1_yyy')


def admitted(owner,module,receipt=None):
    receipt=receipt or module.RECEIPT
    rec=current.parent.linear.read(receipt);raw=current.parent.linear.read(module.NAME)
    assert rec['source_family']==raw['source_family']==owner.family
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    for name in (module.NAME,receipt,Path(module.__file__).name,Path(module.__file__).stem+'_check.py'):
        assert owner.hashes[name]==current.sha(name),name
    for name,digest in rec['input_hashes'].items():assert owner.hashes[name]==digest==current.sha(name),name
    if receipt==getattr(module,'CONSTRUCTION_RECEIPT',None):
        assert rec['all_three_ordinary_axial_trace_identities_passed'] and rec['ordinary_own_field_and_cumulative_trace_rows']==204
    else:
        assert rec['all_passed'] and all(rec[key] for key in module.GATES) and not any(rec[key] for key in module.OPEN)
    return rec


def source_contract(owner):
    owner.assert_graph();c=owner.compact
    compact_rec=admitted(owner,current.previous)
    high_rec=admitted(owner,current.previous.previous)
    radial_module=current.previous.previous.previous.previous
    radial_rec=admitted(owner,radial_module)
    ordinary=admitted(owner,current.outer,current.outer.CONSTRUCTION_RECEIPT)
    repaired=admitted(owner,current.outer.intermediate.previous)
    support_rec=admitted(owner,support_source)
    support_raw=current.parent.linear.read(support_source.NAME)['full_defect_linear_bindings']['reserved_background']
    assert support_rec['checks']['source']['actual_leading_local_Uz_zero_and_pure_swirl_support_guard']
    assert support_rec['checks']['source']['nonlocal_histories_and_fixed_N1p_retained']
    assert len(support_raw['actual_local_Uz_and_angular_increment_zero_rows'])==24
    assert support_raw['leading_swirl_equals_original_power_on_new_support']
    assert support_raw['nonlocal_mean_pressure_and_leading_N1p_not_zeroed']
    assert support_raw['exact_new_and_old_supports']
    assert owner.contract==current.CONTRACT and owner.aliases==current.ALIASES and current.ROWS_PER_SEAM==246
    assert set(owner.charts)=={'reserved','quiet','main','gap','end','tail'} and set(owner.joins)==set(current.SEAMS)
    for name,chart in owner.charts.items():
        for spec in owner.contract:
            packet=chart[spec['group']][spec['field']]
            assert len(packet)==spec['physical_y_order']+1
            assert all(len(rows(jet))==spec['ordinary_Z_order']+1 for jet in packet)
        expected=c.mu if name in ('main','gap') else owner.g.one
        assert chart['physical_y_Jacobian']==expected
        assert chart['fields']['P1']==chart['cumulative']['m13']
    pins=compact_check.contract(c)
    for function in (Algebra.at,Algebra.integral,current.weighted_product):
        pin=current.parent.active.current.ast_binding(function)
        assert owner.hashes[pin['module']]==pin['module_sha256'];pins[function.__name__]=pin
    transport=radial_rec['checks']['own_pressure_transport']
    assert transport['passed'] and transport['own_P1_equals_retained_Rp_history_plus_actual_partial_pressure']
    assert transport['own_P1_y_is_actual_leading_density_by_proved_FTC']
    assert transport['m15_y_is_negative_actual_radius_times_same_density']
    assert radial_rec['checks']['partial_provider_source_identity']['passed']
    assert compact_rec['checks']['actual_outside_own_radial_source']['passed']
    for name,packet in c.own_outside.items():
        assert all(ref==owner.g.zero for jet in packet['Utheta1'] for ref in rows(jet)),name
    assert compact_rec['checks']['actual_total_pressure_source']['passed']
    theta=compact_rec['checks']['actual_total_pressure_source']
    assert theta['theta0_actual_power_node']==compact_check.theta0_context(c)['actual_power_node']
    assert len(theta['theta0_actual_query_independent_amplitude_nodes'])==4
    radial_source=ordinary['checks']['actual_reserved_radial_source']
    assert ordinary['source_radial_calculus_passed'] and radial_source['passed']
    assert radial_source['actual_mass_FTC_M_y_equals_R_U_consumed']
    assert len(radial_source['actual_positive_order_Qy_Vy_C2_rows'])==6 and radial_source['own_one_plus_delta_recovery']
    assert repaired['checks']['actual_pressure_FTC']['passed']
    return dict(passed=True,ordinary_jet_and_integral_source_bindings=pins,
        contracted_rows_per_seam=246,contracted_total_rows=1230,
        accepted_compact_receipt=current.previous.RECEIPT,accepted_higher_receipt=current.previous.previous.RECEIPT,
        accepted_own_FTC_transport_receipt=radial_module.RECEIPT,accepted_ordinary_receipt=current.outer.CONSTRUCTION_RECEIPT,
        primitive_FTC_consumed_before_fixed_density_higher_derivative_lift=True,
        actual_own_patch_pressure_density_and_derivatives_zero_outside_reserved_support=True,
        actual_disjoint_support_identity_and_pure_swirl_power_receipt=support_source.RECEIPT,
        higher_support_derivatives_from_exact_local_function_identity_not_finite_row_count=True,
        arbitrary_order_numeric_bounds_global_n1_not_claimed=True)


def actual_cumulative_calculus(owner):
    c=owner.compact;g=owner.g;a=Algebra(g,3);b=c.repaired
    rec=admitted(owner,current.outer.intermediate.previous)
    assert rec['checks']['actual_densities']['passed']
    assert len(rec['checks']['actual_densities']['actual_selected_density_rows'])==20
    assert rec['checks']['actual_pressure_FTC']['fixed_leading_N1p_plus_patch_density_not_patch_alone']
    # The exact admitted patch integral has these density bodies and fixed
    # endpoints. Consume its FTC before differentiating the actual body.
    variable=b['patch_bodies']['variable'];history=[]
    for key in current.MOMENTS:
        for n,node in enumerate(b['patch_histories'][key]):
            row=g.nodes[node]
            assert row['operation']=='definite_integral' and row['variable']==g.nodes[variable]['name']
            assert row['lower']==0 and row['upper']==c.phase.node
            assert row['integrand']==b['patch_bodies']['bodies'][key][n]
            assert row['physical_measure']=='dR=R ds' and row['fixed_leading_N1p_not_subtracted']
            expected=a.at(c.own_outer.jet(b['patch_bodies']['per_ds'][key]),c.phase,c.own_outer.ref(variable))
            assert encode(expected)==b['patch_bodies']['bodies'][key]
            history.append(key+'_patch_FTC_Z'+str(n))
    Z=sy.Symbol('Z',real=True);Y=sy.Symbol('actual_reserved_y',real=True);Rc=sy.Symbol('Rc',positive=True)
    R=Rc*sy.exp(Y);mu=sy.Symbol('mu',positive=True)
    U,F=sy.Function('actual_Uz1')(Y,Z),sy.Function('actual_F1')(Y,Z)
    theta=sy.Function('actual_theta_amplitude')(Z)*sy.exp(-(sy.Rational(1,2)+mu)*Y)
    env={g.symbol('Z').node:Z,c.phase.node:Y,c.own_R.node:R,c.mu.node:mu}
    # These are input functions to the density product. Their ordinary and
    # physical derivative relations are the admitted compact source theorem.
    for packet,body in ((c.own_reserved['Uz1'],U),(c.own_reserved['F1'],F),(c.theta0_radial,theta)):
        for j,jet in enumerate(packet):bind(env,jet,sy.diff(body,Y,j),Z)
    alg=compact_check.CompactAlgebra(c,env);actual=[]
    for key,body in (('m12',2*R**2*F),('m14',R*sy.sqrt(2*R)*theta*U)):
        saved=c.own_outer.jet(b['patch_bodies']['per_ds'][key])
        for n,ref in enumerate(rows(saved)):
            zero(alg.at(ref)-sy.diff(body,Z,n),('actual reserved density body',key,n))
        for j,jet in enumerate(owner.charts['reserved']['cumulative'][key][1:],1):
            for n,ref in enumerate(rows(jet)):
                zero(alg.at(ref)-sy.diff(body,Y,j-1,Z,n),('actual cumulative Leibniz source',key,j,n))
                actual.append(key+'_y'+str(j)+'_Z'+str(n))
    # m11 uses the admitted exact mass FTC and m13 the admitted total
    # pressure FTC. m15 differentiates the same signed pressure density.
    reserved=owner.charts['reserved']
    assert reserved['cumulative']['m11']==c.own_reserved['M1'][:4]
    assert reserved['cumulative']['m13']==[c.own_pressure[key] for key in PRESSURE]
    pressure_rows=[];independent=[]
    for name,chart in owner.charts.items():
        R=chart['radius'];P=chart['fields']['P1'];mom=chart['cumulative'];query=chart['coordinate'].node
        assert mom['m13']==P
        for j,expected in enumerate((a.scale(P[1],g.neg(R)),a.scale(a.add(P[1],P[2]),g.neg(R)),
            a.scale(a.add(P[1],a.scale(P[2],2),P[3]),g.neg(R))),1):
            assert rows(mom['m15'][j])==rows(expected),(name,j)
            pressure_rows.extend(name+'_m15_y'+str(j)+'_Z'+str(n) for n in range(4))
        if name!='reserved':
            for key in ('m11','m12','m14'):
                for n,ref in enumerate(rows(mom[key][0])):
                    assert not outer_check.free_query_dependency(owner.own,ref.node,query),(name,key,n)
                    independent.append(name+'_'+key+'_Z'+str(n))
                assert all(ref==g.zero for jet in mom[key][1:] for ref in rows(jet))
        if name=='tail':
            for key in ('m13','m15'):
                for ref in rows(mom[key][0]):assert not outer_check.free_query_dependency(owner.own,ref.node,query)
    return dict(passed=True,actual_reserved_integral_FTC_rows=history,
        actual_reserved_m12_m14_y1_to_y3_C3_Leibniz_rows=actual,
        actual_all_chart_signed_m15_y1_to_y3_C3_rows=pressure_rows,
        actual_exterior_nonpressure_cumulative_query_independence_rows=independent,
        actual_fixed_plus_patch_total_pressure_FTC_not_just_leading_density=True,
        actual_reserved_mass_FTC_and_original_nonzero_inlets_retained=True,
        derivatives_of_R_and_sqrt_2R_and_theta0_included=True)


class QuietFTCAlgebra(outer_check.SourceAlgebra):
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i]
        if row['operation']=='definite_integral':
            assert row['endpoint_independent_of_Z'] and row['physical_Jacobian_applied_once']
            t=sy.Symbol(row['variable'],real=True)
            result=sy.Integral(self.at(row['integrand']),(t,self.at(row['lower']),self.at(row['upper'])))
            self.cache[i]=result;return result
        return super().at(i)


def actual_quiet_own_FTC_lift(owner):
    c=owner.compact;g=owner.g;q=owner.high.quiet['source'];a=Algebra(g,3)
    Z=sy.Symbol('Z',real=True);S=sy.Symbol(g.nodes[q['s']]['name'],real=True);R2=sy.Symbol('actual_R2',positive=True)
    R=R2*sy.exp(S);source=sy.Function('actual_reduced_S')(S,Z);Q=sy.Function('actual_Q')(S,Z)
    env={g.symbol('Z').node:Z,q['s']:S,q['R2']:R2}
    bind(env,owner.own.jet(q['S']),source,Z)
    for key,j in (('Q',0),('Q_y',1),('Q_yy',2)):
        bind(env,owner.own.jet(q['physical'][key]),sy.diff(Q,S,j),Z)
    for group in ('incoming_integrals','initial','incoming_pressure'):
        for key,ids in q[group].items():
            if not isinstance(ids,list):continue
            for node in ids:assert not outer_check.free_query_dependency(owner.own,node,q['s'])
            bind(env,owner.own.jet(ids),sy.Function('retained_quiet_'+group+'_'+key)(Z),Z)
    repaired=owner.own.inputs[current.outer.intermediate.previous.NAME]
    for key in current.MOMENTS:
        jet=owner.own.jet(repaired['band_endpoint'][key])
        for ref in rows(jet):assert not outer_check.free_query_dependency(owner.own,ref.node,q['s'])
        bind(env,jet,sy.Function('actual_band_endpoint_'+key)(Z),Z)
    alg=QuietFTCAlgebra(g,env);actual=[]
    Py=owner.charts['quiet']['fields']['P1'][1]
    for key,sign_weight in (('m13',sy.Integer(1)),('m15',-R)):
        expected=a.add(owner.own.jet(repaired['band_endpoint'][key]),
            a.sub(owner.own.jet(q['pressure'][key]),owner.own.jet(q['incoming_pressure'][key])))
        primitive=owner.charts['quiet']['cumulative'][key][0]
        assert rows(expected)==rows(primitive)
        for n,ref in enumerate(rows(primitive)):
            zero(sy.diff(alg.at(ref),S)-sign_weight*alg.at(rows(Py)[n]),('actual quiet own primitive FTC',key,n))
            actual.append(key+'_y_Z'+str(n))
    higher=admitted(owner,current.previous.previous)
    assert higher['checks']['actual_quiet_higher_source']['passed']
    for j,key in ((2,'P1_yy'),(3,'P1_yyy')):
        assert owner.charts['quiet']['fields']['P1'][j]==owner.high.quiet['higher_density'][key]
    return dict(passed=True,actual_quiet_own_pressure_and_weighted_pressure_FTC_rows=actual,
        actual_S0_S1_integrals_and_Q_Qy_finite_boundary_terms_differentiated=True,
        every_incoming_primitive_query_independent_before_FTC=True,
        higher_own_rows_lift_same_primitive_through_admitted_actual_density_derivatives=True)


def row_record(owner,seam,group,key,j,n,method,dependencies):
    packet=owner.joins[seam]
    return dict(seam=seam,group=group,field=key,physical_y_order=j,ordinary_Z_order=n,
        actual_left_node=rows(packet['left'][group][key][j])[n].node,
        actual_right_node=rows(packet['right'][group][key][j])[n].node,
        proof=method,dependencies=dependencies)


def cumulative_endpoint_lemma(owner,seam,key,ordinary,radial_rec,internal_rec):
    """Identify actual primitive traces before using them as recovery data."""
    g=owner.g;a=Algebra(g,3);radial=owner.radial
    definitions={'2Rc':('reserved','quiet'),'Rp':('quiet','main'),
        'main_gap':('main','gap'),'gap_end':('gap','end'),'Rv':('end','tail')}
    left_name,right_name=definitions[seam]
    internal=seam in ('main_gap','gap_end')
    for side,name,edge_index in (('left',left_name,1),('right',right_name,0)):
        chart=owner.charts[name];edge=chart['domain'][edge_index]
        actual=owner.joins[seam][side]['cumulative'][key][0]
        assert actual==a.at(chart['cumulative'][key][0],chart['coordinate'],edge)
        if name in ('reserved','quiet','tail'):
            expected=owner.own.joins[seam][side]['cumulative'][key]
            if name=='tail':
                assert all(not outer_check.free_query_dependency(owner.own,ref.node,chart['coordinate'].node) for ref in rows(expected))
                expected=a.at(expected,chart['coordinate'],g.zero)
            assert rows(actual)==rows(expected),(seam,side,key)
            continue
        packet=radial.current_own[name]
        assert chart['cumulative'][key][0]==packet['cumulative'][key]
        assert rows(actual)==rows(a.at(packet['cumulative'][key],packet['coordinate'],edge))
        if key not in ('m13','m15'):
            assert packet['cumulative'][key]==radial.Rp_cumulative[key]
            old=radial.intermediate['charts']['pulse_'+name]['cumulative'][key]
            assert encode(packet['cumulative'][key])==old
            assert all(not outer_check.free_query_dependency(owner.own,ref.node,chart['coordinate'].node) for ref in rows(packet['cumulative'][key]))
        else:
            increment=radial.pulse.charts[name]['signed_pressure_increment'][key]
            assert packet['cumulative'][key]==a.add(radial.Rp_cumulative[key],increment)
            if not internal:
                records=radial_rec['checks']['own_pressure_transport']['actual_Rp_Rv_m13_m15_endpoint_rows_by_whole_function_identity_composition']
                matches=[row for row in records if row['seam']==seam and row['moment']==key]
                assert len(matches)==4
                for n,row in enumerate(matches):
                    assert row['Z_order']==n and row['actual_current_sum_node']==rows(packet['cumulative'][key])[n].node
                    assert row['actual_old_endpoint_node']==rows(owner.own.joins[seam][side]['cumulative'][key])[n].node
                    other='left' if seam=='Rp' else 'right'
                    assert row['retained_target_node']==rows(owner.own.joins[seam][other]['cumulative'][key])[n].node
    if internal:
        if key in ('m13','m15'):
            proof=internal_rec['checks']['internal_primitive_composition'];assert proof['passed']
            labels=proof['actual_retained_m13_m15_C3_internal_pressure_primitive_rows']
            assert all(seam+'_'+key+'_Z'+str(n) in labels for n in range(4))
            assert proof['actual_common_Rp_cumulative_anchor_not_reset']
        else:
            assert owner.charts[left_name]['cumulative'][key][0]==owner.charts[right_name]['cumulative'][key][0]
    else:
        proof=ordinary['checks']['actual_own_'+seam+'_source_join'];assert proof['passed']
        assert all('cumulative_'+key+'_Z'+str(n) in proof['actual_own_'+seam+'_field_and_cumulative_trace_rows'] for n in range(4))
    return True


def actual_tail_calculus(owner):
    c=owner.compact;g=owner.g;tail=owner.own.inputs[current.outer.intermediate.previous.NAME]
    Z=sy.Symbol('Z',real=True);T=sy.Symbol('actual_tail_y',real=True)
    Lo,Lm=sy.symbols('actual_log_Rv actual_log_2Rc',real=True)
    env={g.symbol('Z').node:Z,owner.tail_coordinate.node:T,
        tail['tail_origins']['pressure_log_radius']:Lo,tail['tail_origins']['mass_log_radius']:Lm,
        c.delta.node:sy.Symbol('delta',positive=True)}
    expected_translation=g.add(owner.tail_coordinate,owner.own.ref(tail['tail_origins']['pressure_log_radius']),
        g.neg(owner.own.ref(tail['tail_origins']['mass_log_radius'])))
    assert owner.tail_mass_coordinate==expected_translation
    assert g.nodes[owner.tail_mass_coordinate.node]==dict(operation='sum',arguments=[owner.tail_coordinate.node,
        tail['tail_origins']['pressure_log_radius'],g.neg(owner.own.ref(tail['tail_origins']['mass_log_radius'])).node])
    for key,ids in tail['terminal'].items():
        jet=owner.own.jet(ids)
        for ref in rows(jet):assert not outer_check.free_query_dependency(owner.own,ref.node,owner.tail_coordinate.node)
        bind(env,jet,sy.Function('actual_terminal_'+key)(Z),Z)
    alg=compact_check.CompactAlgebra(c,env);chart=owner.charts['tail'];actual=[]
    zero(alg.at(chart['radius'])-sy.exp(Lo+T),'actual tail origin radius translation')
    for key,jets in chart['fields'].items():
        base=alg.at(jets[0].value)
        for j,jet in enumerate(jets):
            for n,ref in enumerate(rows(jet)):
                zero(alg.at(ref)-sy.diff(base,T,j,Z,n),('actual saved tail radial source',key,j,n))
                actual.append(key+'_y'+str(j)+'_Z'+str(n))
    return dict(passed=True,actual_saved_tail_radial_and_Z_source_rows=actual,
        actual_mass_coordinate_is_tail_y_plus_logRv_minus_log2Rc=True,
        nonzero_terminal_mass_and_pressure_functions_retained=True,
        actual_pressure_and_other_cumulative_query_independence_checked=True)


def actual_outer_state_joins(owner):
    c=owner.compact;g=owner.g;compact_rec=admitted(owner,current.previous)
    accepted=compact_rec['checks']['actual_2Rc_source_joins']['actual_own_2Rc_y0_to_y4_C3_QV_C2_rows']
    evidence=[]
    for key in current.STATE:
        for side in ('left','right'):
            assert owner.joins['2Rc'][side]['fields'][key]==c.own_2Rc_traces[side][key]
        for j,jet in enumerate(owner.joins['2Rc']['left']['fields'][key]):
            for n in range(len(rows(jet))):
                assert key+'_y'+str(j)+'_Z'+str(n) in accepted
                evidence.append(row_record(owner,'2Rc','fields',key,j,n,'accepted actual compact own higher trace',
                    [current.previous.RECEIPT]))
    # Abstract only the same retained cumulative mass input, after its
    # original source equality and both physical radii have been proved.
    construction=admitted(owner,current.outer,current.outer.CONSTRUCTION_RECEIPT)
    radius_proof=construction['checks']['actual_old_current_radius_sources']
    assert radius_proof['passed']
    Z=sy.Symbol('Z',real=True);mu=sy.Symbol('mu',positive=True)
    rad_rec=admitted(owner,current.previous.previous.previous.previous)
    internal_rec=admitted(owner,current.previous.previous.previous)
    for seam,left_name,right_name in (('Rp','quiet','main'),('main_gap','main','gap'),('gap_end','gap','end'),('Rv','end','tail')):
        cumulative_endpoint_lemma(owner,seam,'m11',construction,rad_rec,internal_rec)
        if seam in ('main_gap','gap_end'):
            internal_radius=internal_rec['checks']['leading_internal_sources']
            assert internal_radius['passed'] and seam in internal_radius['actual_internal_radius_identities']
        Rseam=sy.Symbol('actual_'+seam,positive=True);M=sy.Function('actual_'+seam+'_retained_mass')(Z)
        results={}
        for name,side in ((left_name,'left'),(right_name,'right')):
            chart=owner.charts[name];X=sy.Symbol(g.nodes[chart['coordinate'].node]['name'],real=True)
            edge=chart['domain'][1] if side=='left' else chart['domain'][0]
            # The radius recipe theorem covers old/new2Rc/Rp/Rv. The new
            # tail uses precisely the saved pressure-log origin plus y.
            env=outer_check.radius_context(owner.own,radius_proof);env[c.mu.node]=mu
            env[chart['coordinate'].node]=X
            edge_value=outer_check.SourceAlgebra(g,env).at(edge)
            jac=mu if name in ('main','gap') else sy.Integer(1)
            env[chart['radius'].node]=Rseam*sy.exp((X-edge_value)/jac)
            if name=='tail':
                tail=owner.own.inputs[current.outer.intermediate.previous.NAME]
                env[tail['tail_origins']['pressure_log_radius']]=sy.log(Rseam)
                env[tail['tail_origins']['mass_log_radius']]=sy.Symbol('actual_mass_log_origin',real=True)
            mass=chart['cumulative']['m11'][0]
            bind(env,mass,M,Z)
            alg=compact_check.CompactAlgebra(c,env)
            data=owner.at(name,edge);values={}
            for key in current.STATE:
                for j,jet in enumerate(data['fields'][key]):
                    for n,ref in enumerate(rows(jet)):values[(key,j,n)]=alg.at(ref)
            results[side]=values
        assert set(results['left'])==set(results['right'])
        for (key,j,n),left in results['left'].items():
            zero(left-results['right'][(key,j,n)],('actual exterior own mixed recovery join',seam,key,j,n))
            evidence.append(row_record(owner,seam,'fields',key,j,n,'actual recovery from same admitted mass and physical radius',
                [current.outer.CONSTRUCTION_RECEIPT,current.previous.RECEIPT]))
    assert len(evidence)==750
    return dict(passed=True,actual_state_mixed_join_rows=evidence,
        same_nonzero_cumulative_mass_input_not_mean_or_Q_outputs_bound=True,
        original_physical_radius_source_theorem_consumed=True,own_one_plus_delta_recovery_preserved=True)


def actual_pressure_and_cumulative_joins(owner):
    c=owner.compact;g=owner.g;radial=owner.radial;high=owner.high
    comp=admitted(owner,current.previous);old=admitted(owner,current.outer,current.outer.CONSTRUCTION_RECEIPT)
    rad=admitted(owner,current.previous.previous.previous.previous)
    higher=admitted(owner,current.previous.previous)
    internal=admitted(owner,current.previous.previous.previous)
    transport=rad['checks']['own_pressure_transport'];assert transport['passed']
    composition=transport['actual_Rp_Rv_m13_m15_endpoint_rows_by_whole_function_identity_composition']
    pressure=[];cumulative=[]
    for seam in current.SEAMS:
        # Value identities come from actual defining primitive endpoints.
        # For current pulse pressure the canonical/current whole-function
        # identity is explicitly consumed before its endpoint transport.
        for key in current.MOMENTS:
            cumulative_endpoint_lemma(owner,seam,key,old,rad,internal)
            for n in range(4):
                cumulative.append(row_record(owner,seam,'cumulative',key,0,n,'actual primitive endpoint theorem plus current-source transport',
                    [current.outer.CONSTRUCTION_RECEIPT,current.previous.previous.previous.previous.RECEIPT]))
        for j in range(4):
            if j==0:
                assert all(owner.joins[seam][side]['fields']['P1'][0]==owner.joins[seam][side]['cumulative']['m13'][0]
                    for side in ('left','right'))
                source='same retained cumulative m13 value theorem'
            else:
                if seam=='2Rc':
                    packet=owner.own.pressure_radial_traces[seam] if j==1 else c.pressure_2Rc_traces
                    expected={side:packet[side] if j==1 else packet[side][PRESSURE[j]] for side in ('left','right')}
                    theorem=current.outer.CONSTRUCTION_RECEIPT if j==1 else current.previous.RECEIPT
                    if j==1:assert len(old['checks']['actual_own_2Rc_source_join']['actual_2Rc_radial_pressure_trace_rows'])==4
                    else:assert len(comp['checks']['actual_2Rc_source_joins']['actual_total_own_2Rc_pressure_and_m15_C3_rows'])==16
                elif seam in ('Rp','Rv'):
                    packet=radial.pressure_traces[seam] if j==1 else high.outer_pressure_traces[seam]
                    expected={side:packet[side] if j==1 else packet[side][PRESSURE[j]] for side in ('left','right')}
                    theorem=current.previous.previous.previous.previous.RECEIPT if j==1 else current.previous.previous.RECEIPT
                    if j>1:
                        labels=higher['checks']['actual_outer_higher_pressure_joins']['actual_Rp_Rv_higher_pressure_C3_rows']
                        assert all(seam+'_'+PRESSURE[j]+'_Z'+str(n) in labels for n in range(4))
                else:
                    left_name,right_name=('main','gap') if seam=='main_gap' else ('gap','end')
                    if j==1:
                        expected=dict(left=high.pulse.endpoints[left_name]['right']['density']['P1_y'],
                            right=high.pulse.endpoints[right_name]['left']['density']['P1_y'])
                        labels=internal['checks']['leading_internal_sources']['actual_radial_pressure_and_weighted_density_C3_internal_rows']
                        theorem=current.previous.previous.previous.RECEIPT
                    else:
                        expected=dict(left=high.endpoints[left_name]['right']['higher_density'][PRESSURE[j]],
                            right=high.endpoints[right_name]['left']['higher_density'][PRESSURE[j]])
                        labels=higher['checks']['actual_internal_higher_pressure_joins']['actual_internal_higher_pressure_C3_rows']
                        theorem=current.previous.previous.RECEIPT
                    assert all(seam+'_'+PRESSURE[j]+'_Z'+str(n) in labels for n in range(4))
                for side in ('left','right'):
                    actual_rows=rows(owner.joins[seam][side]['fields']['P1'][j]);expected_rows=rows(expected[side])
                    assert len(actual_rows)==len(expected_rows)
                    for lr,rr in zip(actual_rows,expected_rows):
                        if lr!=rr:
                            # Endpoint substitution can retain a literal-zero
                            # wrapper. Check that exact recipe, not its ID.
                            alg0=outer_check.SourceAlgebra(g,{})
                            assert alg0.at(lr)==alg0.at(rr)==0,(seam,j,side,lr.node,rr.node)
                source='actual own primitive FTC then admitted fixed-density derivative join '+theorem
            for n in range(4):pressure.append(row_record(owner,seam,'fields','P1',j,n,source,
                [current.previous.RECEIPT,current.previous.previous.RECEIPT,current.previous.previous.previous.previous.RECEIPT]))
        # The state mass rows and pressure rows are already identified above.
        # m12/m14 are flat at2Rc by the actual guarded own bump source.
        # Elsewhere their exact primitives are query independent.
        if seam=='2Rc':
            env=compact_check.own_context(c)['bindings'];env.update(compact_check.theta0_context(c)['bindings'])
            env[c.phase.node]=sy.log(2);alg=compact_check.CompactAlgebra(c,env)
            for key in ('m12','m14'):
                for jet in owner.joins[seam]['left']['cumulative'][key][1:]:
                    for ref in rows(jet):zero(alg.at(ref),('actual flat cumulative2Rc density',key,ref.node))
        for key in current.MOMENTS:
            for j in range(1,4):
                if key=='m11':
                    assert all(owner.joins[seam][side]['cumulative'][key][j]==owner.joins[seam][side]['fields']['M1'][j]
                        for side in ('left','right'))
                    method='same admitted state mass derivative'
                elif key=='m13':method='same proved total pressure derivative'
                elif key=='m15':method='actual signed radius-weighted pressure derivative Leibniz composition'
                else:
                    if seam!='2Rc':
                        alg0=outer_check.SourceAlgebra(g,{})
                        for side in ('left','right'):
                            for ref in rows(owner.joins[seam][side]['cumulative'][key][j]):
                                zero(alg0.at(ref),('actual exterior zero density derivative',seam,key,j,side))
                    method='actual flat density at2Rc or query-independent exterior primitive'
                for n in range(4):cumulative.append(row_record(owner,seam,'cumulative',key,j,n,method,
                    ['actual_cumulative_calculus','actual_outer_state_joins','actual own pressure FTC join']))
    assert len(pressure)==80 and len(cumulative)==400
    return dict(passed=True,actual_total_own_pressure_join_rows=pressure,actual_cumulative_join_rows=cumulative,
        fixed_density_not_relabelled_without_actual_own_FTC=True,
        nonzero_primitive_values_and_all_finite_boundary_terms_retained=True,
        m15_derivative_joins_use_same_proved_radius_and_pressure_rows=True,
        actual_patch_density_flatness_consumed_at2Rc=True)


CHECKS={'source_contract','actual_cumulative_calculus','actual_quiet_own_FTC_lift','actual_tail_calculus','actual_outer_state_joins','actual_pressure_and_cumulative_joins'}


def save(owner,checks,execution_seconds):
    owner.assert_graph();assert set(checks)==CHECKS and all(proof['passed'] for proof in checks.values())
    evidence=checks['actual_outer_state_joins']['actual_state_mixed_join_rows']+checks['actual_pressure_and_cumulative_joins']['actual_total_own_pressure_join_rows']+checks['actual_pressure_and_cumulative_joins']['actual_cumulative_join_rows']
    expected={(seam,spec['group'],spec['field'],j,n) for seam in owner.joins for spec in owner.contract
        for j in range(spec['physical_y_order']+1) for n in range(spec['ordinary_Z_order']+1)}
    observed={(r['seam'],r['group'],r['field'],r['physical_y_order'],r['ordinary_Z_order']) for r in evidence}
    assert len(evidence)==len(observed)==1230 and observed==expected
    raw=current.parent.linear.read(current.NAME)
    assert raw['combined_outer_mixed_source_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        contracted_rows_per_seam=246,actual_contracted_mixed_join_rows=1230,
        exact_graph_nodes=len(raw['exact_graph_nodes']),finite_combined_outer_source_accepted=True,
        numerical_oracle_available=False,all_order_smoothness_quantitative_bounds=False,
        global_n1=False,
        execution_seconds=execution_seconds)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_FINITE_COMBINED_OUTER_1230_MIXED_ROWS_ACCEPTED',flush=True)
    return rec


def run(owner):
    began=time.monotonic();checks={}
    for fn in (source_contract,actual_cumulative_calculus,actual_quiet_own_FTC_lift,actual_tail_calculus,actual_outer_state_joins,actual_pressure_and_cumulative_joins):
        checks[fn.__name__]=fn(owner);print('COMBINED_OUTER_'+fn.__name__.upper()+'_PASSED',flush=True)
    return save(owner,checks,time.monotonic()-began)


if __name__=='__main__':
    import sys
    sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
    run(current.CurrentOriginalN1CombinedOuterMixedSourceFunctions(require_checked=False))

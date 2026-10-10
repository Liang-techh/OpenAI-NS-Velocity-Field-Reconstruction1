"""Independent nonunit-measure Duhamel/FTC and chained O2/O3 admission."""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_outer_C3_continuation as current
import lei_ren_part1_paper_compliant_current_original_Rh_C3_continuation_check as rh_check

symbolic=rh_check.symbolic


def encoded(value):return json.loads(json.dumps(current.target.encoded(current.ranges.record(value))))


def geometric_measure(field):
    offsets,domains,x=current.current.exact_geometry();before='Rh_reference';rows=[]
    for chart in current.CHARTS:
        left,right=domains[chart];p_right=domains[before][1]
        assert s.simplify(offsets[before].subs(x,p_right)-offsets[chart].subs(x,left))==0
        J=s.diff(offsets[chart],x)
        assert J==(40*s.exp(40*x) if chart=='O2_axial' else 1)
        width=s.simplify(offsets[chart].subs(x,right)-offsets[chart].subs(x,left))
        assert width==dict(O2_slope=1,O2_axial=s.exp(40)-1,O2_buffer=11,O3_transition=1,O3_power=2)[chart]
        # Change of variable dlogR=J dx is exact, including the exponential chart.
        assert s.simplify(s.integrate(J,(x,left,right))-width)==0
        w=field.windows[chart];geo=field.geometry[chart]
        assert w['incoming']==(field.rh.window if before=='Rh_reference' else field.windows[before])['outgoing']
        assert geo['phase']==field.graph.unary('fractional_part',field.graph.mul(field.N,
            current.source.FunctionRef(field.graph,geo['offset']))).node
        parts=w['original_physical_measure_and_partitions']
        assert parts[0]['exact_left']==[int(left),1] and parts[-1]['exact_right']==[int(right),1]
        assert all(a['exact_right']==b['exact_left'] for a,b in zip(parts,parts[1:]))
        rows.append(dict(chart=chart,exact_physical_width=str(width),exact_Jacobian=str(J),
            same_original_absolute_phase_and_predecessor=True))
        before=chart
    return dict(actual_five_chart_radius_measure_and_seams=rows,nonunit_axial_Jacobian_applied_once=True)


def history_FTC_and_outgoing(field):
    g=field.graph;partial_rows=partition_rows=leading_rows=alias_rows=endpoint_rows=0
    for chart in current.CHARTS:
        f=field.functions[chart];w=field.windows[chart];geo=field.geometry[chart]
        t=s.Symbol('native_'+chart);a=s.Symbol(chart+'_history_endpoint')
        # Unit-Jacobian charts contain the shared exact-one node. Do not replace
        # that global constant by an arbitrary derivative in other expressions.
        L=s.Function('log_radius_'+chart) if chart=='O2_axial' else (lambda q:q)
        left,right=map(s.Integer,field.domains[chart]);J=s.diff(L(t),t)
        base_bindings={geo['offset']:L(t),geo['left_offset']:L(left),geo['right_offset']:L(right),
            geo['width']:L(right)-L(left)}
        if geo['Jacobian']!=g.one.node:base_bindings[geo['Jacobian']]=J
        for kind,inlets,densities,partials,locals_ in (
            ('correction',f['source_incoming_N_scaled_C3'],f['source_density_N_scaled_C3'],f['partial_N_scaled_C3'],f['local_partial_N_scaled_C3']),
            ('leading',f['actual_leading_incoming_C3'],f['actual_leading_density_C3'],f['actual_leading_history_C3'],f['actual_leading_local_C3'])):
            for key,rate in current.current.RATES.items():
                r=s.Rational(str(rate));memory=s.exp(-r*(L(a)-L(left)))
                for j,(h,d,q,local) in enumerate(zip(*(current.target.rows(v[key]) for v in (inlets,densities,partials,locals_)))):
                    H=s.Integer(0) if h==g.zero else s.Symbol('H_'+str(h.node))
                    D=s.Integer(0) if d==g.zero else s.Function('D_'+str(d.node))(t)
                    bindings={**base_bindings}
                    if h!=g.zero:bindings[h.node]=H
                    if d!=g.zero:bindings[d.node]=D
                    expected=memory*H+s.Integral(s.exp(-r*(L(a)-L(t)))*D*J,(t,left,a)) if D!=0 else memory*H
                    got=symbolic(g,q.node,bindings)
                    assert s.expand(got-expected)==0,(chart,key,j,kind,'Duhamel')
                    assert s.simplify(s.diff(got,a)-s.diff(L(a),a)*(D.subs(t,a)-r*got))==0,(chart,key,j,kind,'physical FTC')
                    assert s.simplify(got.subs(a,left)-H)==0,(chart,key,j,kind,'inlet')
                    if d==g.zero:assert local==g.zero
                    else:
                        node=g.nodes[local.node]
                        assert node['variable']=='native_'+chart and node['lower']==geo['domain'][0]
                        assert node['upper']==f['endpoint'].node and node['ordinary_slow_Z_derivative_order']==j
                    if kind=='leading':leading_rows+=1
                    else:
                        partial_rows+=1
                        expected_out=s.exp(-r*(L(right)-L(left)))*H;cell_integrals=[]
                        if D!=0:
                            cell_integrals=[s.Integral(s.exp(-r*(L(right)-L(t)))*D*J,
                                (t,s.Rational(*p['exact_left']),s.Rational(*p['exact_right'])))
                                for p in w['original_physical_measure_and_partitions']]
                            expected_out+=sum(cell_integrals,s.Integer(0))
                        old_out=symbolic(g,w['outgoing'][key][j],bindings)
                        assert s.expand(old_out-expected_out)==0,(chart,key,j,'old partition')
                        # Explicitly compare the *new function* at the right
                        # endpoint with the old API alias after finite-integral
                        # additivity, rather than just checking the alias ID.
                        direct=s.expand(got.subs(a,right));full=s.Integer(0)
                        if D!=0:
                            full=s.expand(s.Integral(s.exp(-r*(L(right)-L(t)))*D*J,(t,left,right)))
                        normalized=s.expand(old_out)
                        expected_integrals=[s.expand(q) for q in cell_integrals]
                        assert normalized.atoms(s.Integral)==set(expected_integrals),(chart,key,j,'actual partition integrands')
                        collapsed=normalized.xreplace({q:s.Integer(0) for q in expected_integrals})+full
                        assert s.expand(direct-collapsed)==0,(chart,key,j,'direct substituted endpoint equals alias')
                        endpoint_rows+=1
                        partition_rows+=1
        assert encoded(field.correction_functions(chart,int(left)))==w['incoming']
        assert encoded(field.correction_functions(chart,int(right)))==w['outgoing'];alias_rows+=40
        for bad in (int(left)-1,int(right)+1):
            try:field.correction_functions(chart,bad)
            except ValueError:pass
            else:raise AssertionError('Invalid native extension admitted')
    assert partial_rows==leading_rows==partition_rows==endpoint_rows==100 and alias_rows==200
    return dict(independent_partial_correction_Duhamel_and_physical_FTC_rows=partial_rows,
        independent_leading_Duhamel_and_physical_FTC_rows=leading_rows,
        independently_checked_original_partitioned_outgoing_rows=partition_rows,
        independently_checked_new_partial_endpoint_equals_old_alias_rows=endpoint_rows,
        exact_incoming_and_outgoing_alias_rows=alias_rows,
        integral_additivity_on_original_partitions=True,physical_Jacobian_present_once=True,
        leading_and_correction_nonzero_inlets_retained=True)


def signed_complete_chain(field):
    g=field.graph;alg=current.target.C3Algebra(g);z=s.Symbol('Z');N=s.Symbol('N',positive=True)
    count=seam_rows=quiet_rows=leading_density_rows=0
    for chart in current.CHARTS:
        f=field.functions[chart];w=field.windows[chart];p=field.raw['actual_C3_primitive_functions'][chart]
        jet=lambda row:current.target.C3Function(*(current.source.FunctionRef(g,i) for i in row))
        inputs=[jet(p['source_roots'][key]) for key in ('E','V')]
        inputs.extend([alg.fixed(0) if w['F_N'] is None else jet(w['F_N']),jet(p['B'])])
        bindings={field.N.node:N};abstract=[]
        for name,q in zip(('E0','V0','F','B'),inputs):
            expr=s.Integer(0) if all(row==g.zero for row in current.target.rows(q)) else s.Function(chart+'_'+name)(z)
            for j,row in enumerate(current.target.rows(q)):
                if row!=g.zero:bindings[row.node]=s.diff(expr,z,j)
            abstract.append(expr)
        E0,V0,F,B=abstract;E,V=E0+F/N,V0+B/N
        leading=dict(m=V0,h=E0,k=E0*V0,e=V0**2-E0**2/2,p=E0**2/2)
        changed=dict(m=V,h=E,k=E*V,e=V**2-E**2/2,p=E**2/2)
        for key in current.current.RATES:
            for j,row in enumerate(current.target.rows(f['actual_leading_density_C3'][key])):
                assert s.expand(symbolic(g,row.node,bindings)-s.diff(leading[key],z,j))==0,(chart,key,j,'actual leading density from source roots')
                leading_density_rows+=1
            for j,row in enumerate(current.target.rows(f['source_density_N_scaled_C3'][key])):
                assert s.cancel(symbolic(g,row.node,bindings)-s.diff(N*(changed[key]-leading[key]),z,j))==0,(chart,key,j,'signed density')
                count+=1
        complete={}
        for key in current.current.RATES:
            L,H=(s.Function(chart+'_'+prefix+key)(z) for prefix in ('leading_','correction_'))
            for rows,expr in ((f['actual_leading_history_C3'][key],L),(f['partial_N_scaled_C3'][key],H)):
                for j,row in enumerate(current.target.rows(rows)):
                    if row!=g.zero:bindings[row.node]=s.diff(expr,z,j)
            complete[key]=L+H/N
            for j,row in enumerate(current.target.rows(f['actual_complete_history_C3'][key])):
                assert s.expand(symbolic(g,row.node,bindings)-s.diff(complete[key],z,j))==0;count+=1
        first=dict(m=V-complete['m'],h=E-3*complete['h']/2,k=E*V-3*complete['k']/2,
            e=V**2-E**2/2-complete['e'],p=E**2/2)
        for key,q in f['actual_complete_history_y_C3'].items():
            for j,row in enumerate(current.target.rows(q)):
                assert s.expand(symbolic(g,row.node,bindings)-s.diff(first[key],z,j))==0,(chart,key,j,'complete first-y')
                count+=1
        P0=s.Function('same_original_P0')(z)
        for j,row in enumerate(current.target.rows(f['actual_independent_P0_C3'])):bindings[row.node]=s.diff(P0,z,j)
        for j,row in enumerate(current.target.rows(f['actual_absolute_pressure_C3'])):
            assert s.expand(symbolic(g,row.node,bindings)-s.diff(P0+complete['p'],z,j))==0;count+=1
        assert encoded(f['actual_independent_P0_C3'])==encoded(field.rh.functions['actual_independent_P0_C3'])
        predecessor=field.initial_leading if chart==current.CHARTS[0] else field.functions[current.CHARTS[current.CHARTS.index(chart)-1]]['actual_leading_outgoing_C3']
        assert encoded(f['actual_leading_incoming_C3'])==encoded(predecessor);seam_rows+=20
        for key,q in f['actual_leading_outgoing_C3'].items():
            for row,partial in zip(current.target.rows(q),current.target.rows(f['actual_leading_history_C3'][key])):
                n=g.nodes[row.node]
                assert n['operation']=='function_substitution' and n['expression']==partial.node
                assert n['variable']==f['endpoint'].node and n['value']==f['domain'][1]
        if chart=='O3_power':
            assert F==B==0
            assert all(row==g.zero for q in f['source_density_N_scaled_C3'].values() for row in current.target.rows(q))
            assert any(q.value!=g.zero for q in f['source_incoming_N_scaled_C3'].values());quiet_rows=20
    assert count==320 and seam_rows==100 and quiet_rows==20 and leading_density_rows==100
    return dict(independent_signed_density_complete_first_y_and_pressure_rows=count,
        independent_leading_density_from_original_E_V_product_rule_rows=leading_density_rows,
        exact_leading_history_predecessor_alias_rows=seam_rows,quiet_local_zero_with_nonzero_incoming_rows=quiet_rows,
        original_same_P0_function_all_charts=True,leading_is_source_Duhamel_not_a_range_value=True)


def quantitative_ranges(field,raw):
    c=field.c;cells=0;windows=0
    for saved in raw['actual_four_Z_outer_C3_partial_ranges']:
        ends=tuple(saved['exact_Z_cell']);assert encoded(field.quantitative_chain(ends))==saved
        before=encoded(field.rh.quantitative_range(ends))
        for row in saved['actual_five_chart_C3_partial_ranges']:
            chart=row['chart'];windows+=1
            for kind in ('leading','N_scaled'):
                assert row['actual_'+kind+'_incoming_C3_bounds']==before['actual_'+kind+'_partial_history_C3_bounds']
                for key in current.current.RATES:
                    for bound,old in zip(row['actual_'+kind+'_partial_history_C3_bounds'][key],row['actual_'+kind+'_incoming_C3_bounds'][key]):
                        if old['exact_zero']:continue
                        assert current.ep(current.current.packets.interval(c,bound['log_absolute_upper']))[1]>=current.ep(current.current.packets.interval(c,old['log_absolute_upper']))[1]
            for cell in row['original_native_cells']:
                left,right=(s.Rational(*cell[key]) for key in ('exact_left','exact_right'))
                width=s.simplify(field.offsets[chart].subs(field.x,right)-field.offsets[chart].subs(field.x,left))
                assert cell['exact_physical_width']==str(width)
                actual=current.interval_constant(c,width)
                assert encoded(actual)==cell['physical_width'] and current.ep(actual)[0]>0
                for key,rate in current.current.RATES.items():
                    r=c.mpf(rate.numerator)/rate.denominator;mass=actual if rate==0 else -c.expm1(-r*actual)/r
                    assert encoded(mass)==cell['original_positive_own_rate_masses'][key] and current.ep(mass)[0]>0
                cells+=1
            assert row['range_covers_all_native_partial_endpoints'] and row['range_caps_not_function_values']
            assert row['actual_independent_P0_C3_bounds']==before['actual_independent_P0_C3_bounds']
            before=row
    assert cells==68 and windows==20
    assert all(len(rows)==4 for rows in field.P0_registry.values())
    return dict(actual_four_axial_cells=4,actual_five_chart_windows=windows,actual_native_cells=cells,
        true_radius_width_masses_and_nonunit_axial_measure_checked=True,
        previous_whole_partial_majorant_retained_without_memory_reset=True,
        all_partial_endpoints_covered_and_original_pressure_retained=True)


def run(field=None):
    began=time.monotonic();raw=current.rh.read(current.NAME)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    assert raw['candidate_actual_outer_C3_continuation_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    with mp.workdps(540):
        field=field if field is not None else current.CurrentOuterC3Continuation(require_checked=False)
        assert not field.acceptance_loaded and field.graph.nodes==raw['exact_graph_nodes']
        assert field.graph.nodes[:len(field.prefix)]==field.prefix
        assert encoded(field.functions)==raw['actual_five_chart_C3_continuation_functions']
        assert encoded(field.geometry_bindings)==raw['actual_geometry_bindings']
        assert encoded(field.P0_registry)==raw['original_same_P0_source_registry']
        assert raw['source_bindings']==current.source_bindings() and raw['source_family']==field.identity
        geometry=geometric_measure(field);histories=history_FTC_and_outgoing(field)
        signed=signed_complete_chain(field);bounds=quantitative_ranges(field,raw)
        assert not raw['final_selected_native_pulse_postpulse_collar_heat_registry_installed']
        assert not raw['current_numeric_point_field_oracle_installed']
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            independent_native_geometry_and_measure=geometry,independent_partial_histories_and_partition=histories,
            independent_signed_complete_history_chain=signed,actual_whole_native_directed_ranges=bounds,
            **dict.fromkeys(current.OPEN,False),current_numeric_point_field_oracle_installed=False,
            global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
            final_selected_native_pulse_postpulse_collar_heat_registry_installed=False,
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
            execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.target.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual outer C3 native-measure history chain, nonzero quiet memory and whole-domain bounds checks passed',flush=True)
    return result


if __name__=='__main__':run()

"""Role/context/fail-closed checks with new original factored source queries."""
import copy
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_Rc_factored_source_oracle as current

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha


def role_checks(owner,built,oracle):
    g=built['graph'];checked=0;expected={}
    for cell in built['cells']:
        if cell['source_flat_exact_zero']: continue
        view=owner.owner.views[cell['chart']]
        for key,pair in cell['contributions'].items():
            for order in ('value','Z'):
                integral=g.nodes[pair[order]];product=g.nodes[integral['integrand']]
                indices=[i for i in product['arguments'] if g.nodes[i]['operation']=='original_function_graph']
                if len(indices)!=1: raise AssertionError('Source function count differs')
                row=g.nodes[indices[0]]
                role='density_'+key+('_Z' if order=='Z' else '_C0')
                root=(view['five_signed_increment_rate_roots'][key] if order=='value' else
                    view['five_signed_increment_rate_first_derivatives']['Z'][key])
                if row['source_node']!=root or row['source_graph_namespace']!='function_graph_nodes':
                    raise AssertionError('Wrong full signed density namespace/root')
                if row['function_role']!=role or oracle.require_row(row)!=role:
                    raise AssertionError('Wrong original signed density dispatch')
                expected[indices[0]]=role;checked+=1
    for pairkey,role,order in (('value','Rc_E_C0','y0_Z0'),('Z','Rc_E_Z','y0_Z1')):
        index=getattr(built['amplitude'],pairkey).node;row=g.nodes[index]
        view=owner.owner.views['O3_power']['original_signed_input_graph']['jet_expression_dag']
        if row['source_graph_namespace']!='original_signed_input_graph.jet_expression_dag' or row['source_node']!=view['roots']['E'][order]:
            raise AssertionError('Wrong original amplitude DAG/derivative')
        if oracle.require_row(row)!=role: raise AssertionError('Wrong Rc amplitude role')
        expected[index]=role;checked+=1
    references=[i for i,row in enumerate(g.nodes) if row['operation']=='original_function_graph']
    if set(references)!=set(expected) or checked!=212:
        raise AssertionError('Untagged or missing original source function reference')
    return dict(passed=True,role_namespace_root_bindings=checked,
        integer_node_ID_without_graph_namespace_or_role_not_accepted=True,
        old_accepted_function_artifacts_unchanged=True)


def negative_checks(oracle,built,frame):
    g=built['graph'];row=next(q for q in g.nodes if q.get('chart')=='O2_slope' and q.get('function_role')=='density_m_C0')
    cases=[]
    for name,patch in (
        ('missing_role',{'function_role':None}),
        ('wrong_namespace',{'source_graph_namespace':'original_signed_input_graph.jet_expression_dag'}),
        ('wrong_source_node',{'source_node':row['source_node']+123456}),
        ('wrong_source_hash',{'graph_sha256':'different'}),
        ('unknown_role',{'function_role':'density_UNKNOWN_C0'})):
        try: oracle.require_row({**row,**patch})
        except (ValueError,KeyError): cases.append(name)
        else: raise AssertionError('Wrong original source accepted: '+name)
    unresolved=current.FactoredSourceFrame('O2_slope',frame.N,None,{'status':'requires_source_or_phase_refinement'})
    try: oracle.dispatch_function_range(row,unresolved)
    except ArithmeticError: cases.append('unresolved_no_zero_fallback')
    else: raise AssertionError('Unresolved source supplied a value')
    try: oracle.dispatch_function_range(row,current.FactoredSourceFrame('inner_reference',frame.N,frame.values,frame.record))
    except ValueError: cases.append('wrong_chart')
    else: raise AssertionError('Wrong chart supplied a value')
    # Rc amplitudes and density increments must not be confused even if two
    # IDs have the same numerical value in their independent graph namespaces.
    amplitude=g.nodes[built['amplitude'].value.node]
    try: oracle.require_row({**amplitude,'source_graph_namespace':'function_graph_nodes'})
    except ValueError: cases.append('amplitude_density_namespace_collision')
    else: raise AssertionError('Amplitude/density graph namespace collision admitted')
    try: current.transport.evaluate(built,built['targets']['M'].value,oracle=oracle,Z='.5',N=2048,ctx=oracle.ctx)
    except current.transport.SourceOracleRequired: cases.append('factored_oracle_not_plain_scalar_oracle')
    else: raise AssertionError('Factored ranges silently cast to scalar evaluator')
    return dict(passed=True,rejection_cases=cases)


@current.density.native.inlet.source_precision
def run(owner,report=None):
    began=time.monotonic();built=owner.build();oracle=current.NativeFactoredSourceOracle(owner,built)
    report=json.loads((HERE/current.NAME).read_bytes()) if report is None else report
    if not report[current.GATE] or report['source_family']!=oracle.source_family:
        raise AssertionError('Current original range producer required')
    if report['declared_density_queries']!=10 or report['enclosed_density_queries']!=10 or report['actual_role_dispatched_factored_ranges']!=82:
        raise AssertionError('Declared selected original queries not all resolved')
    if report['candidate_N_queries']!=[1024,2048]: raise AssertionError('Two distinct actual N queries required')
    for row in report['original_factored_source_query_records'].values():
        if row['status']!='enclosed' or row['saved_cover_midpoints_or_plain_scalar_cast_used']:
            raise AssertionError('Original range query replaced by scalar/cap selection')
        if not row['phase_derived_from_original_radius_not_supplied_as_free_angle'] or not row['all_phase_cells_unioned_with_source_functions_retained']:
            raise AssertionError('Different original phase/source convention')
        if set(row['signed_density_C0_ranges'])!=set(current.transport.RATES) or set(row['signed_density_Z_ranges'])!=set(current.transport.RATES):
            raise AssertionError('Missing signed density row')
        if 'separate_original_P0' not in row or 'separate_original_P0_Z' not in row:
            raise AssertionError('Original pressure datum merged or lost')
    c=oracle.ctx;ep=current.transport.current.ep
    box=c.mpf((ep(c.mpf('.13369999'))[0],ep(c.mpf('.13370001'))[1]))
    frame=oracle.density_frame(chart='O2_slope',Z=('.49','.51'),coordinate=box,N=2048)
    if frame.values is None: raise AssertionError('Declared original O2 radial/Z cell unresolved')
    rows={(q.get('chart'),q.get('function_role')):q for q in built['graph'].nodes if q['operation']=='original_function_graph'}
    count=0
    for key in current.transport.RATES:
        for order in ('C0','Z'):
            value=oracle.dispatch_function_range(rows['O2_slope','density_'+key+'_'+order],frame)
            if value.ctx is not c or value is not frame.values[order][key]:
                raise AssertionError('Dispatcher changed original range, context or units')
            count+=1
    # An actual new local integral query at N=2048, using the accepted
    # original signed source backend. This is a local interval result,
    # not full original Rc numerical transport or target/control values.
    integral=oracle.owner.contribution(Z=('.49','.51'),left='.13369999',right='.13370001',N=2048)
    if any(value.ctx is not c for group in ('contributions','Z_derivatives') for value in integral[group].values()):
        raise AssertionError('Original local integral context changed')
    if not any(not value.zero for value in integral['contributions'].values()):
        raise AssertionError('Nontrivial original source cell collapsed to zero')
    amplitude=oracle.amplitude_frame(Z=('.49','.51'),N=2048)
    for role in ('Rc_E_C0','Rc_E_Z'):
        value=oracle.dispatch_function_range(rows['O3_power',role],amplitude)
        if value.ctx is not c: raise AssertionError('Rc amplitude source context changed')
        count+=1
    broad=oracle.density_frame(chart='O2_slope',Z=(-1,1),coordinate=c.mpf((ep(c.mpf('.12'))[0],ep(c.mpf('.15'))[1])),N=2048)
    if broad.values is None:
        try: oracle.dispatch_function_range(rows['O2_slope','density_m_C0'],broad)
        except ArithmeticError: pass
        else: raise AssertionError('Actual unresolved broad source box silently evaluated')
    result=dict(all_passed=True,**{current.GATE:True},source_family=oracle.source_family,
        explicit_role_bindings=role_checks(owner,built,oracle),
        fail_closed_and_namespace_checks=negative_checks(oracle,built,frame),
        additional_original_cell_and_amplitude_dispatches=count,
        additional_original_O2_radial_Z_cell=frame.record,
        additional_actual_original_local_integral_C0_Z=integral['record'],
        additional_actual_original_local_integral_rows=10,
        actual_broad_O2_whole_Z_probe=broad.record,
        broad_probe_status=broad.record['status'],
        inherited_original_point_queries=10,inherited_actual_role_dispatches=82,
        numerical_original_source_ranges_queried=True,actual_original_local_integrals_queried=True,
        scalar_original_point_oracle_installed=False,all_17_chart_original_range_oracle_admitted=False,
        full_factored_function_graph_evaluator_installed=False,
        actual_original_full_route_numerical_integrals_evaluated=False,
        actual_five_controls_installed=False,actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(current.packets.OPEN,False),execution_seconds=time.monotonic()-began,
        input_hashes={**oracle.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='Role/namespace/context dispatch and fail-closed guards, selected genuine original ranges at two N values, one additional original radial/Z cell and ten new local integral C0/Z enclosures. Original backend source math is inherited unchanged; no whole-route target/control/commonN/closure/recursion admission.')
    (HERE/current.RECEIPT).write_text(json.dumps(current.packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original factored range oracle focused check PASS; broad O2 probe:',broad.record['status'],flush=True)
    return result

"""Bounded check of new actual power source, exact limit and terminal band."""
import ast
import copy
import gzip
import inspect
import json
from pathlib import Path
import textwrap
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_limit_repair_band as current


def source_witnesses():
    owner=current.outer.WholeZAllNOuterRcFunctions()
    service=current.ActualReservedPowerSource(owner);witnesses=[]
    for ends in current.outer.CELLS:
        op=owner.owner.owner(tuple(ends));f=op.flow
        inlet=service.query(ends,(1,1))
        packet=owner.leading_packet(ends,'O3_power',(2,1))
        old=packet['original_closed_O3_background']
        raw=lambda q:q['raw']['raw_current_radius_y_derivative_axial_coefficients']
        for group in ('histories','velocity'):
            assert raw(inlet)[group].keys()==raw(old)[group].keys()
            for key in raw(old)[group]:
                for newrows,oldrows in zip(raw(inlet)[group][key],raw(old)[group][key]):
                    assert current.outer.previous.equivalent_rows(newrows,oldrows),(ends,group,key)
        # Wide actual source cells retain all original histories, not merely
        # terminal amplitude. Their ranges are witnesses, never graph values.
        whole=service.query(ends,(1,1),(2,1));exit_=service.query(ends,(2,1))
        generic=packet['original_generic_source'];invS=f.factor((0,-.5,0,0,0))
        assert generic['common_original_P0_axial5'] is packet['exact_common_P0_axial5'] is op.P0
        for key in current.current.RATES:
            expected=[q*invS for q in raw(old)['histories'][key][0]] if key in ('m','k') else raw(old)['histories'][key][0]
            assert current.outer.previous.equivalent_rows(expected,generic['common_own_five_histories_axial5'][key]),key
        assert whole['exact_common_P0_axial5'] is exit_['exact_common_P0_axial5'] is op.P0
        current.outer.core.relative.same_source(f,[q for name in ('histories','velocity')
            for rows in raw(whole)[name].values() for row in rows for q in row])
        assert all(q.zero for rows in raw(whole)['velocity']['axial'] for q in rows)
        assert all(q.zero for q in whole['actual_collected_Delta_axial5'][1:])
        assert not whole['original_positive_formal_mu'].zero
        assert current.ep(service.reservation['actual_Tw'])[0]>current.ep(2+owner.c.ln(2))[1]
        assert current.ep(service.reservation['actual_quiet_power_log_margin'])[0]>=0
        assert not raw(exit_)['velocity']['theta'][0][0].zero
        witnesses.append(dict(exact_Z_cell=list(ends),exact_Rc_inlet_rows_identical=True,
            same_live_original_P0=True,whole_reserved_power_cell_queried=True,
            actual_band_exit_nonzero_amplitude_retained=True,
            original_background_and_radial_Z_axial5_rows_retained=True,
            actual_generic_history_leaf_paths_and_original_m_k_normalization_checked=True,
            exact_mu_quiet_power_and_Tw_reservation_checked=True,
            actual_source_rows_are_range_witnesses_not_graph_values=True))
    return witnesses


def run():
    began=time.monotonic();accepted,hashes=current.load_current()
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert report[current.GATE] and report['source_family']==accepted['source_family']
    for name,digest in report['input_hashes'].items():assert current.sha(name)==digest,name
    # Independent AST comparison: every original background statement from
    # key onward survives this explicit new coordinate adapter unchanged.
    old=ast.parse(textwrap.dedent(inspect.getsource(current.original.background_cell))).body[0]
    split=next(i for i,n in enumerate(old.body) if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='key' for t in n.targets))
    suffix=ast.Module(body=copy.deepcopy(old.body[split:]),type_ignores=[])
    assert current.hashlib.sha256(ast.dump(suffix).encode()).hexdigest()==report['original_power_background_source_binding']['unchanged_original_background_suffix_AST_sha256']
    for bad in (True,1,(0,1),(3,1),(1,0),(1.0,1),(2,1,1)):
        try:current.band_fraction(bad)
        except ValueError:pass
        else:raise AssertionError('Invalid reserved source coordinate admitted')
    try:current.original.fraction('power',(3,1))
    except ValueError:pass
    else:raise AssertionError('Historical original power guard was relaxed')
    built=current.build(accepted);g=built['graph']
    assert g.nodes==report['exact_graph_nodes']
    assert current.current.encode_graph(built)==report['exact_limit_repair_band']
    assert g.nodes[:built['original_graph_prefix_length']]==accepted['exact_function_graph_nodes']
    vector=g.nodes[built['exact_limit_vector'].node]
    assert vector['selected_integer']==accepted['actual_whole_Z_frequency_connection']['selected_integer']
    assert vector['shared_N']==built['N'].node
    assert vector['source_connection_report_sha256']==current.sha(current.current.NAME)
    assert len(built['exact_limit_controls'])==5
    linear_system=g.nodes[built['implicit_Z_linear_system'].node]
    assert linear_system['matrix']==current.current.encode_graph(built['implicit_Z_matrix'])
    assert linear_system['rhs']==current.current.encode_graph(built['implicit_Z_rhs'])
    for i,pair in enumerate(built['exact_limit_controls']):
        for k,q in enumerate((pair.value,pair.Z)):
            row=g.nodes[q.node]
            assert row['operation']=='C1_Banach_limit_component'
            assert row['component']==i and row['Z_order']==k and row['limit_vector']==built['exact_limit_vector'].node
            if k==1:assert row['implicit_Z_linear_system']==built['implicit_Z_linear_system'].node
    assert 'finite_picard_sequence' not in built and 'finite_terminal_residual_identities' not in built
    # Actual residual expressions remain in the graph; a theorem proves
    # their value/Z identity, rather than writing zero as an input field.
    for q,proof in zip(built['control_residual'],built['residual_zero_theorem_nodes']):
        assert q.value!=g.zero and g.nodes[proof.node]['expression_pair']==current.current.encode_graph(q)
    for key in current.current.RATES:
        pair=built['partial_band_histories'][key]
        assert pair.value!=g.zero
        proof=g.nodes[built['relative_terminal_zero_certificates'][key].node]
        assert proof['expression_pair']==current.current.encode_graph(built['relative_terminal_Duhamel_functions'][key])
        assert proof['equivalent_residual_pair']==current.current.encode_graph(built['exact_limit_terminal_residual_identities'][key])
        assert proof['relative_not_absolute_exterior']
    background=built['original_power_background']
    assert background['P0_not_reset_or_merged_into_p_history']
    assert background['V'].value==background['V'].Z==g.zero
    assert background['E']==background['corrected_fields']['original_E']
    fields=background['corrected_fields'];gmu=built['parameters']['mu']
    expected_y=g.c1add(g.c1scale(g.neg(g.add(g.constant('1/2'),gmu)),fields['original_E']),
        g.c1mul(built['amplitude'],fields['F_y']))
    assert fields['E_y']==expected_y
    assert background['leading_history_P0_source_binding']==current.leading_history_binding(accepted)
    assert len(background['leading_histories'])==len(background['complete_history_y_Z_pairs'])==5
    for key in ('E_y','V_y','F_y','G_y'):
        assert isinstance(background['corrected_fields'][key],current.source.C1Function)
    assert report['symbolic_theorems']==current.symbolic_theorems()
    with mp.workdps(540):witnesses=source_witnesses()
    for key in ('actual_numeric_controls_evaluated','actual_global_frequency_admitted',
        'physical_original_exterior_five_targets_closed','full_recovered_velocity_pressure_and_heat_joins_admitted',
        'actual_temporal_scale_recursion_installed',*current.outer.OPEN):assert report[key] is False,key
    hashes[current.NAME]=current.sha(current.NAME);hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=report['source_family'],
        exact_limit_repair_band_graph_nodes_checked=len(g.nodes),
        exact_source_limit_not_finite_iterate_checked=True,
        current17_exact_graph_prefix_reused_without_reintegration=True,
        original_background_AST_suffix_and_strict_reserved_coordinate_API_checked=True,
        four_Z_actual_reserved_power_source_witnesses=witnesses,
        original_five_ODE_semigroup_and_terminal_Z_identities_checked=True,
        original_bump_support_field_jet_joins_checked=True,
        relative_and_absolute_closure_gates_separate=True,
        input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_SOURCE_LIMIT_REPAIR_BAND',len(g.nodes),'nodes; four live source witnesses; relative C1 closure',flush=True)
    return receipt


if __name__=='__main__':run()

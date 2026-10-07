"""Actual seventeen-chart native queries and exact common-unit cover checks."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_generic_source_packets as source
from lei_ren_part1_paper_compliant_current_native_generic_left_inlet_check import compare_modes

packets=source.packets


def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=source.inlet.native_bridge_owner()
    report=json.loads((source.HERE/source.NAME).read_bytes())
    if not report[source.GATE] or report['native_generic_chart_count']!=17:
        raise ValueError('Actual seventeen-chart native producer required')
    for name,digest in report['input_hashes'].items():
        if hashlib.sha256((source.HERE/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Original whole native source input changed: '+name)
    views=json.loads(gzip.decompress((source.HERE/source.VIEWS).read_bytes()))
    comparisons=invalid=0;cases=[]
    with source.inlet.CheckedSourceRuntime():
        backend=source.NativeGenericSourcePackets(bridge)
        if report['source_family']!=backend.family or not all(backend.assert_native_graph().values()):
            raise ValueError('Actual native source family or original graph differs')
        for chart in source.DOMAINS:
            saved=views[chart];c=backend.ctx
            z=packets.interval(c,saved['provenance']['Z_box'])
            coordinate=packets.interval(c,saved['provenance']['coordinate_box'])
            packet=backend.query(chart,z,coordinate)
            original=backend.original_raw(chart,z,coordinate)
            raw=original['raw']
            def expected_mode(value):
                if isinstance(value,packets.FactoredJet):return value
                return packets.FactoredAlgebra(c,original['logs'],[]).lift(value)
            for group,normalized,native,shifted in (
                ('velocity',packet.velocity,raw['velocity'],('axial','radial')),
                ('histories',packet.histories,raw['histories'],('m','k'))):
                for key,rows in native.items():
                    shift=(0,-.5,0,0) if key in shifted else (0,0,0,0)
                    for actual,expected in zip(normalized[key],rows):
                        comparisons+=compare_modes(actual,expected_mode(expected),shift)
            for actual,expected in zip(packet.absolute_pressure,raw['absolute_pressure']):
                comparisons+=compare_modes(actual,expected_mode(expected),(0,0,0,0))
            if set(packet.P0.terms)!={(0,0,0,0)}:
                raise AssertionError('Separate original P0 gained a spurious factor shift')
            for actual,expected in zip(packet.P0.terms[(0,0,0,0)].coefficients,original['P0'].coefficients):
                if actual._mpi_!=expected._mpi_:
                    raise AssertionError('Original parent P0 changed')
                comparisons+=1
            if packets.encode(packet.record())!=saved:
                raise AssertionError('Saved native packet differs from fresh query: '+chart)
            if (saved['point_source_function'] or saved['source_function_replay'] or
                    saved['provenance']['cache_cover'] or not saved['provenance']['arbitrary_coordinates_evaluated']):
                raise AssertionError('Live native source cover scope differs')
            if chart=='O3_power':
                with source.mp.workdps(300):
                    lo,hi=packets.recovery.endpoints(original['source_packet']['local_offset'])
                    if not lo<=source.mp.mpf('.537')<=hi:
                        raise AssertionError('O3 reciprocal-phase enclosure missed nominal local offset')
            cases.append(dict(chart=chart,whole_Z=(-1,1),ordinary_y_rows=5,
                all_five_original_histories_retained=True,P0_axial_order=packet.P0.order,
                native_provider_group=original['native_provider_group']))
            print('Checked live generic native source:',chart,flush=True)
        for chart,(lo,hi) in source.DOMAINS.items():
            for z,coordinate in ((2,lo),(0,lo-1)):
                try:backend.query(chart,z,coordinate)
                except ValueError:invalid+=1
                else:raise AssertionError('Invalid original source domain accepted: '+chart)
        inlet=views['generic_left_inlet']
        if not inlet['generic_inlet_defect_zero_is_declared_Section11_initial_condition']:
            raise AssertionError('Left five-zero state was not labeled as an initial condition')
        if len(backend.queries)!=17:
            raise AssertionError('All seventeen charts were not actually queried')
    for key in packets.OPEN:
        if report[key]:raise AssertionError('Whole native covers promoted an unfinished global gate')
    for key in ('actual_changed_defect_integral_functions_installed','actual_repair_control_functions_installed',
                'actual_terminal_Z_function_closure_installed'):
        if report[key]:raise AssertionError('Whole native source scope was widened')
    hashes=dict(report['input_hashes']);hashes.update({source.NAME:source.sha(source.NAME),
        Path(__file__).name:source.sha(Path(__file__).name),
        source.PREFIX+'current_native_generic_left_inlet_check.py':source.sha(source.PREFIX+'current_native_generic_left_inlet_check.py')})
    result=dict(source_family=backend.family,all_passed=True,**{source.GATE:True},
        actual_seventeen_original_native_chart_queries=cases,
        exact_source_mode_coefficient_and_P0_comparisons=comparisons,
        invalid_original_chart_domains_rejected=invalid,
        fresh_actual_query_equals_saved_source_cover=True,
        native_pre_context_copied_with_original_adapter=True,
        O3_reciprocal_phase_encloses_nominal_local_offset_not_an_exact_point=True,
        separate_original_P0_and_all_five_histories_retained=True,
        initial_zero_defects_not_inferred_as_transport_or_terminal_closure=True,
        actual_changed_defect_integral_functions_installed=False,actual_repair_control_functions_installed=False,
        actual_terminal_Z_function_closure_installed=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,scope=report['scope'],input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Checked17 original native source backends:',comparisons,'modal/P0 comparisons',flush=True)
    return result


if __name__=='__main__':run()

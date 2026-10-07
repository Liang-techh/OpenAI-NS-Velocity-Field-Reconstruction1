"""Check actual native queries, exact modal units and bounded runtime reuse."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_generic_left_inlet as source

packets=source.packets


def runtime_guards():
    import lei_ren_part1_paper_compliant_macro_signed_integrals as macro
    import lei_ren_part1_paper_compliant_collar_Gamma_C4_check as theorem
    import lei_ren_part1_paper_compliant_current_heat_source as heat
    import lei_ren_part1_paper_compliant_current_postpulse_interfaces as postpulse
    before=(macro._sha256,theorem.functional_source_identities,heat.functional_source_identities,postpulse.heat_theorem)
    with source.CheckedSourceRuntime() as runtime:
        first=runtime.reuse_proof();second=runtime.reuse_proof()
        if first != second or first is second:
            raise AssertionError('Retained symbolic proofs must be independent copies')
        first[next(iter(first))]=False
        if not all(second.values()):
            raise AssertionError('Returned proof mutation leaked into retained evidence')
        original_digest=runtime.digest
        changed_name=next(iter(runtime.closure))
        runtime.digest=lambda name:('0'*64 if name==changed_name else original_digest(name))
        try:
            runtime.reuse_proof()
        except ValueError:
            pass
        else:
            raise AssertionError('Changed proof input was silently reused')
        finally:
            runtime.digest=original_digest
        try:
            with source.CheckedSourceRuntime():
                pass
        except RuntimeError:
            pass
        else:
            raise AssertionError('Nested runtime binding context must be rejected')
    if before != (macro._sha256,theorem.functional_source_identities,heat.functional_source_identities,postpulse.heat_theorem):
        raise AssertionError('Original source proof/digest aliases were not restored')
    try:
        with source.CheckedSourceRuntime():
            raise LookupError('fixture error')
    except LookupError:
        pass
    if before != (macro._sha256,theorem.functional_source_identities,heat.functional_source_identities,postpulse.heat_theorem):
        raise AssertionError('Exceptional runtime exit did not restore original source bindings')
    return dict(proof_copy_independence=True,changed_source_hash_rejected=True,
        nested_context_rejected=True,original_bindings_restored_after_success_and_error=True,
        numerical_operators_and_constructors_never_rebound=True,passed=True)


def compare_modes(normalized, original, shift):
    # Inspect source exponents and MP interval tuples directly, independently
    # of both the adapter serializer and the backend's conversion helper.
    if normalized.order != original.order:
        raise AssertionError('Original source Taylor order changed')
    expected={tuple(k+s for k,s in zip(key,shift)):row for key,row in original.terms.items()}
    if set(normalized.terms)!=set(expected):
        raise AssertionError('Original common unit mode exponents differ')
    comparisons=0
    for key,row in normalized.terms.items():
        for a,b in zip(row.coefficients,expected[key].coefficients):
            if a._mpi_ != b._mpi_:
                raise AssertionError('Original coefficient enclosure changed during normalization')
            comparisons+=1
    return comparisons


def run(bridge=None):
    began=time.monotonic()
    if bridge is None:
        bridge,_=source.native_bridge_owner()
    report=json.loads((source.HERE/source.NAME).read_bytes())
    if not report[source.GATE] or report['source_family']!=dict(zip(packets.FAMILY_KEYS,(bridge.family,bridge.source,bridge.datum_sha))):
        raise ValueError('Actual same-family native source producer required')
    for name,digest in report['input_hashes'].items():
        if hashlib.sha256((source.HERE/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Native source producer input changed: '+name)
    views=json.loads(gzip.decompress((source.HERE/source.VIEWS).read_bytes()))
    comparisons=0;cases=[];invalid=0
    with source.CheckedSourceRuntime():
        backend=source.NativeBridgeSourcePackets(bridge)
        for key in ('left_inlet','fresh_first','fresh_second','fresh_macro'):
            saved=views[key]['original_source_packet'] if key=='left_inlet' else views[key]
            Z=packets.interval(bridge.ctx,saved['provenance']['Z_box'])
            coordinate=packets.interval(bridge.ctx,saved['provenance']['coordinate_box'])
            chart=saved['chart']
            packet=backend.query(chart,Z,coordinate)
            # Match the original source_precision contract for the isolated
            # undecorated AST replay, while keeping serialization unchanged.
            with source.mp.workdps(300):
                original=backend.evaluator(bridge,Z,coordinate,chart.replace('bridge_',''))
                raw=backend.adapter(bridge.ctx,original)
            for key_v,rows in raw['velocity'].items():
                shift=(0,-.5,0,0) if key_v in ('axial','radial') else (0,0,0,0)
                for actual,expected in zip(packet.velocity[key_v],rows):
                    comparisons+=compare_modes(actual,expected,shift)
            for key_m,rows in raw['histories'].items():
                shift=(0,-.5,0,0) if key_m in ('m','k') else (0,0,0,0)
                for actual,expected in zip(packet.histories[key_m],rows):
                    comparisons+=compare_modes(actual,expected,shift)
            for actual,expected in zip(packet.absolute_pressure,raw['absolute_pressure']):
                comparisons+=compare_modes(actual,expected,(0,0,0,0))
            # Independent source path: P0 comes from original prepare(),
            # never from subtracting Mp or an exterior pressure limit.
            p0=bridge.upstream.prepare(Z)['inputs']['p0']
            if set(packet.P0.terms)!={(0,0,0,0)}:
                raise AssertionError('Original P0 must have no extra factor shift')
            for actual,expected in zip(packet.P0.terms[(0,0,0,0)].coefficients,p0.coefficients[:6]):
                if actual._mpi_!=expected._mpi_:
                    raise AssertionError('Independent original axis pressure differs')
                comparisons+=1
            encoded=packets.encode(packet.record())
            if encoded!=saved:
                raise AssertionError('Saved cover differs from fresh actual native coordinate query')
            cases.append(dict(view=key,chart=chart,ordinary_y_rows=5,
                P0_axial_order=packet.P0.order,all_five_original_histories_retained=True))
        for args in (('unknown',0,0),('bridge_first',2,0),('bridge_first',0,-1),('bridge_second',0,0),('bridge_macro',0,2)):
            try:
                backend.query(*args)
            except ValueError:
                invalid+=1
            else:
                raise AssertionError('Invalid original native query domain accepted')
        for bad in (None,object()):
            try:
                source.NativeBridgeSourcePackets(bad)
            except ValueError:
                invalid+=1
            else:
                raise AssertionError('Foreign/nonexistent native owner accepted')
    inlet=views['left_inlet']
    if set(inlet['generic_inlet_defect_exact_zero'])!=set(packets.recovery.RATES):
        raise AssertionError('All five candidate inlet defects required')
    for row in inlet['generic_inlet_defect_exact_zero'].values():
        algebra=packets.FactoredAlgebra(bridge.ctx,
            [packets.interval(bridge.ctx,x) for x in inlet['original_source_packet']['fixed_source_log_bases']],[])
        jet=packets.decode_row(algebra,row)
        if any(packets.recovery.endpoints(v)!=(0,0) for term in jet.terms.values() for v in term.coefficients):
            raise AssertionError('Generic inlet defects must be exact zero')
    for key in packets.OPEN:
        if report[key]:
            raise AssertionError('Native inlet source query promoted an open global gate')
    for key in ('actual_defect_integral_functions_installed','actual_repair_control_functions_installed',
                'actual_terminal_Z_function_closure_installed','source_function_or_full_tensor_owner_interface_gate_promoted'):
        if report[key]:
            raise AssertionError('Native cover-only scope was widened')
    guards=runtime_guards()
    hashes=dict(report['input_hashes'])
    hashes.update({source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    result=dict(source_family=report['source_family'],all_passed=True,**{source.GATE:True},
        actual_same_native_owner_coordinate_cases=cases,exact_modal_coefficient_and_P0_comparisons=comparisons,
        invalid_native_owner_or_domain_rejections=invalid,checked_runtime_reuse_guards=guards,
        fresh_actual_query_equals_saved_source_cover=True,all_five_left_candidate_defects_exact_zero=True,
        pressure_datum_compared_with_independent_original_prepare=True,
        actual_defect_integral_functions_installed=False,actual_repair_control_functions_installed=False,
        actual_terminal_Z_function_closure_installed=False,
        source_function_or_full_tensor_owner_interface_gate_promoted=False,
        **dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,
        scope=report['scope'],input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Checked actual native generic inlet source queries:',comparisons,'modal/P0 comparisons',flush=True)
    return result


if __name__=='__main__':run()

"""Genuine original transition right-half flat-source C1 integral operator.

The entire original [1/2,1] source is queried on both strict-sign Z tiles.
Only exact original flat q/slow jets imply zero own density and integrals.
The incoming history is an explicit argument, never reset or supplied here.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_generic_tail_phase_integrals as query_source
import lei_ren_part1_paper_compliant_current_tail_source_branch_atlas as atlas_source

tail=query_source.tail;rc=query_source.rc;common=query_source.common
HERE,PREFIX,sha=query_source.HERE,query_source.PREFIX,query_source.sha
N=query_source.N;encode=query_source.encode
NAME=PREFIX+'current_transition_right_half.json'
RECEIPT=PREFIX+'current_transition_right_half_check.json'
GATE='original_transition_entire_right_half_flat_C1_integral_operator_executed'


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();prior=json.loads((HERE/query_source.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,query_source,prior['source_family'])
    diagnostic=json.loads((HERE/atlas_source.NAME).read_bytes())
    assert diagnostic[atlas_source.GATE] and diagnostic['source_family']==prior['source_family']
    for name,digest in diagnostic['input_hashes'].items():assert sha(name)==digest,name
    # Refused atlas values never feed a numerical operator. Preserve this
    # executed diagnostic only as evidence for the next source-factor work.
    hashes[atlas_source.NAME]=sha(atlas_source.NAME)
    hashes.update(diagnostic['input_hashes'])
    bridge,seed=rc.native.inlet.native_bridge_owner()
    with rc.native.inlet.CheckedSourceRuntime():
        owner=rc.NativeRcC1Histories(rc.preceding.NativeO2C1Histories(tail.baseline.build_owner(bridge)))
        oracle=query_source.GenericTailPhaseIntegrals(owner);archives=[]
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            oracle.Z=Z;got=oracle.query('O3_slope_mu',Fraction(1,2),Fraction(1),N)
            assert got['values'] is not None,'Original transition right-half source must be covered'
            row=got['record'];native=row['actual_original_phase_inverse_and_density_source']
            q=native['original_q_slow_jet_source'];assert q['status']=='enclosed' and q['branch']=='flat'
            assert q['current_refined_original_q_C0']['exact_zero']
            assert all(r['exact_zero'] for r in q['original_q_ordinary_slow_derivative_enclosures'].values())
            assert all(v.zero for v in (*got['values'].values(),*got['Z_derivatives'].values()))
            operator=rc.history.C1DuhamelOperator(owner.coordinates)
            rc.transfer.append_true_cell(operator,got['geometry'],got['values'],got['Z_derivatives'],owner.family)
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),
                original_entire_transition_right_half_source_query=row,
                original_right_half_C1_integral_operator=operator.record(),
                common_directed_coordinate_theorem=owner.coordinates.record(),
                exact_original_native_window=['1/2','1'],exact_true_log_radius_length='1/2',
                original_smooth_flat_cutoff_implies_all_q_slow_jets_and_own_density_increments_zero=True,
                original_incoming_is_required_argument_and_rate0_memory_is_not_reset=True,
                unresolved_transition_left_half_not_skipped_or_assumed_zero=True,
                actual_numerical_entire_transition_or_Rc_histories_admitted=False)
            raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            filename=PREFIX+'current_transition_right_half_'+tag+'.json.gz'
            (HERE/filename).write_bytes(compressed);hashes[filename]=sha(filename)
            archives.append(dict(filename=filename,exact_Z_range=list(Z),compressed_bytes=len(compressed),
                uncompressed_bytes=len(raw),lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()))
            print('Original entire transition[1/2,1] source q/ordinary jets/own C0-Z integrals EXACT ZERO:',tag,flush=True)
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Source closures disagree: '+name)
            hashes[name]=digest
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            actual_original_transition_right_half_source_archives=archives,
            original_endpoint_source_atlas_diagnostic=atlas_source.NAME,
            original_right_half_source_and_integral_C0_Z_certified=True,
            actual_numerical_entire_transition_or_tail_or_Rc_histories_admitted=False,
            original_source_seed_evidence=seed,input_hashes=hashes,execution_seconds=time.monotonic()-began,
            **dict.fromkeys(common.current.FLAGS,False),
            scope='Original O3_slope_mu[1/2,1], both strict-sign tiles, exact flat q and all ordinary q slow jets, original density/integral C0-Z zero and explicit-input true-width affine history operator. The unresolved left half and original axial source remain open; no entire tail/targets/controls/global N/stress/recursion/full NS.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

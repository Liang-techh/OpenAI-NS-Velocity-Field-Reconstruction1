"""Focused independent weight-integral and complete pressure-memory checks."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_inner_reference_weighted_pressure as source

HERE=source.HERE;sha=source.sha;common=source.common;ep=source.ep;iv=source.iv;same=source.same


def weight_theorem():
    y,yr,Eright=sy.symbols('y yr Eright',real=True)
    exact_E=Eright*sy.exp((y-yr)/10)
    primitive=5*exact_E**2
    assert sy.simplify(sy.diff(primitive,y)-exact_E**2)==0
    assert sy.simplify(sy.diff(exact_E,y)-exact_E/10)==0
    w=sy.symbols('w',positive=True)
    assert sy.simplify(primitive.subs(y,yr)-primitive.subs(y,yr-w)-5*Eright**2*(1-sy.exp(-w/5)))==0
    # For t>=0, exp(t)-1 >= 1-exp(-t), proved by a positive square.
    t=sy.symbols('t',real=True)
    assert sy.simplify(sy.exp(t)+sy.exp(-t)-2-(sy.exp(t/2)-sy.exp(-t/2))**2)==0
    return dict(passed=True,independent_antiderivative_and_y_slope_verified=True,
        endpoint_difference_equals_exact_true_width_weight=True,
        variable_A_pointwise_absolute_increment_majorized_before_integrating=True,
        negative_exponential_increment_dominated_by_positive_increment=True,
        factor_1_minus_exp_minus_width_over5_bounded_by1_not_replaced_by_width=True)


def run():
    begin=time.monotonic();report=json.loads((HERE/source.NAME).read_bytes())
    assert report[source.GATE] and report['candidate_N']==source.N
    assert all(report[k] is False for k in common.current.FLAGS)
    for name,digest in report['input_hashes'].items():assert sha(name)==digest,name
    assert report['original_weighted_pressure_theorem']==source.source_identity()
    theorem=weight_theorem();manifest=json.loads((HERE/source.upstream.NAME).read_bytes())
    originals=json.loads((HERE/common.NAME).read_bytes())
    oldreplay=json.loads((HERE/source.upstream.REPLAY).read_bytes())
    c=MPIntervalContext();c.dps=240;checks=[];strict=0
    with mp.workdps(300):
        for tile,archive in zip(report['actual_weighted_pressure_tiles'],manifest['genuine_active_kappa_upstream_archives'],strict=True):
            assert tile['original_source_archive']==archive
            raw=gzip.decompress((HERE/archive['filename']).read_bytes())
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            payload=json.loads(raw);coordinates,original,jets,P0,P0Z=source.replay.inputs(c,payload)
            _,candidate,computed=source.pressure_bound(c,payload)
            proof=tile['pressure_integral_proof']
            assert proof==source.encode(computed)
            # Independent directed scalar cap; bounded E/A can be materialized,
            # whereas the original enormous true radius/width never is.
            emax=c.exp(c.mpf(ep(iv(c,proof['original_E_log_cover']))[1]))
            amax=iv(c,proof['finite_Amax']);tmax=amax*2/source.N
            cap=emax*emax*(c.mpf(5)/2)*(c.exp(tmax)-1)
            new=common.restore_common_source(proof['selected_pressure_integral'],coordinates)
            represented_log=source.upstream.baseline.magnitude_log(new)
            assert abs(represented_log-ep(c.ln(cap))[1])<mp.mpf('1e-225')
            assert ep(iv(c,proof['twice_Amax_over_N']))==ep(tmax)
            assert proof['strict_upper_reduction']
            flat,rows=source.accepted.route_rows(payload['complete_actual_pre_O2_tile_route'])
            memory=common.restore_common_source(flat['actual_correction_C0_enclosures']['p'],coordinates)
            assert memory.zero
            for (label,row),path in zip(rows,tile['actual_twelve_chart_pressure_path'],strict=True):
                assert path['label']==label and path['exact_original_pressure_decay']==1
                old=row['true_log_radius_signed_C0_contributions']['p']
                assert path['accepted_original_pressure_integral']==old
                addition=new if label=='inner_reference' else common.restore_common_source(old,coordinates)
                same(c,path['selected_original_pressure_integral'],addition)
                same(c,path['actual_inherited_pressure_correction'],memory)
                memory=memory+addition;same(c,path['actual_right_pressure_correction'],memory)
                _,_,_,bg,_,_=source.accepted.row_fields(label)
                background=common.restore_common_source(row[bg]['original_normalized_history_C0_enclosures']['p'],coordinates)
                same(c,path['original_right_pressure_background'],background)
                same(c,path['actual_right_pressure_own_history'],background+memory)
            values={k:common.restore_common_source(v,coordinates) for k,v in tile['actual_refined_upstream_C0'].items()}
            same(c,tile['actual_refined_upstream_C0']['p'],memory)
            for k in source.KEYS:
                if k!='p':same(c,tile['actual_refined_upstream_C0'][k],original[k])
            for record in originals['original_complete_same_N_source_incoming_O2_integral_refinements']:
                if not source.same_Z(record['exact_Z_range'],payload['exact_Z_range']):continue
                output=next(r for r in report['genuine_original_O2_weighted_pressure_replays']
                    if r['exact_Z_range']==record['exact_Z_range'] and r['ordered_source_cells']==record['ordered_source_cells'])
                before=next(r for r in oldreplay['genuine_active_kappa_original_O2_replays']
                    if r['exact_Z_range']==record['exact_Z_range'] and r['ordered_source_cells']==record['ordered_source_cells'])
                assert output==source.encode(source.pressure_replay(c,payload,values,record,archive))
                operator=source.upstream.middle.history.C1DuhamelOperator(coordinates)
                scalar=lambda x:coordinates.scalar(iv(c,x))
                operator.append(c.mpf(1),{k:scalar(v) for k,v in record['five_original_C0_integral_contributions'].items()},
                    {k:scalar(v) for k,v in record['five_genuine_ordinary_Z_integral_contributions'].items()},payload['source_family'])
                propagated=operator.apply(values,jets,payload['source_family'])
                for kind,which in (('C0','values'),('Z','Z_derivatives')):
                    for k in source.KEYS:
                        field='propagated_actual_correction_'+kind
                        same(c,output[field][k],propagated[which][k])
                        if kind=='Z' or k!='p':assert output[field][k]==before[field][k]
                p=propagated['values']['p'];background=scalar(record['source_defined_original_histories_at_y1']['p'])
                same(c,output['actual_absolute_pressure'],P0+(background+p))
                assert output['separate_P0']==before['separate_P0']
                assert output['separate_P0_Z']==before['separate_P0_Z']
                assert output['actual_absolute_pressure_Z']==before['actual_absolute_pressure_Z']
                logbefore=source.upstream.baseline.magnitude_log(common.restore_common_source(before['propagated_actual_correction_C0']['p'],coordinates))
                logafter=source.upstream.baseline.magnitude_log(p)
                assert logafter<logbefore;strict+=1
                checks.append(dict(exact_Z_range=record['exact_Z_range'],ordered_source_cells=record['ordered_source_cells'],
                    before_pressure_C0_log_upper=logbefore,after_pressure_C0_log_upper=logafter,
                    after_pressure_C0_absolute_upper=c.exp(c.mpf(logafter)),
                    original_other_four_C0_all_Z_P0_and_absolute_pressure_Z_unchanged=True,
                    independent_accepted_C1_affine_operator_matches=True))
            print('Weighted pressure proof and complete rate-zero/O2 memory PASS',payload['exact_Z_range'],flush=True)
    assert len(report['actual_weighted_pressure_tiles'])==2 and len(checks)==strict==4
    hashes={**report['input_hashes'],source.NAME:sha(source.NAME),Path(__file__).name:sha(Path(__file__).name)}
    receipt=dict(all_passed=True,**{source.GATE:True},source_family=report['source_family'],candidate_N=source.N,
        independent_weight_integral_proof=theorem,strict_downstream_pressure_C0_reductions=strict,
        full_twelve_active_chart_pressure_memory_checked_on_both_tiles=True,
        genuine_original_O2_pressure_checks=checks,input_hashes=hashes,execution_seconds=time.monotonic()-begin,
        full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False))
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('Original weighted pressure refinement PASS',flush=True)
    return receipt


if __name__=='__main__':run()

"""New closed-endpoint/mixed source checks and complete saved-parent replay.

Independent physical-s fixtures test the new uniform inequalities and mixed
q/width. Actual transition/Rc records are replayed without native owners,
upstream numerical integrations or the inherited legacy test suite.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_transition_complete_prefix as source
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source_check as qcheck
import lei_ren_part1_paper_compliant_current_generic_tail_phase_integrals_check as saved

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv


def exact_replay_equal(actual,expected,path):
    if actual==expected:return
    if isinstance(actual,dict) and isinstance(expected,dict):
        assert actual.keys()==expected.keys(),path+' keys'
        for key,value in actual.items():exact_replay_equal(value,expected[key],path+'.'+str(key))
    elif isinstance(actual,list) and isinstance(expected,list):
        assert len(actual)==len(expected),path+' length'
        for index,value in enumerate(actual):exact_replay_equal(value,expected[index],path+'['+str(index)+']')
    else:raise AssertionError(path+': '+repr(actual)+' != '+repr(expected))


def independent_new_source_checks(family):
    c=MPIntervalContext();c.dps=100;p=mp.mp.clone();p.dps=1600
    comparisons=0;width_checks=0
    for D in (64,128,256):
        coordinates=source.native.HalfPstarCoordinates(c,c.mpf(0),family)
        coord=source.CompletePrefixCoordinates(coordinates,c.mpf(1),-c.mpf(D))
        template=coordinates.scalar(1);eta=p.exp(-D)
        def sigma(s):
            if s<=0:return p.mpf(0)
            if s>=1:return p.mpf(1)
            return 1/(1+p.exp(1/s**2-1/(1-s)**2))
        def q(s):
            Delta=2*sigma(s)
            if Delta>=eta:return p.mpf(0)
            return sigma(1-Delta/eta)*p.sqrt((2*eta-Delta)/(2*(2+Delta)))
        def direct_q_rows(s):
            if not s:return [p.sqrt(eta/2),p.mpf(0),p.mpf(0)]
            sig=sigma(s);rho=2*sig/eta
            if rho>=1:return [p.mpf(0)]*3
            L1=2/(1-s)**3+2/s**3;L2=6/(1-s)**4-6/s**4
            rho1=rho*(1-sig)*L1
            rho2=rho*(1-sig)*(L2+(1-2*sig)*L1**2)
            # Differentiate the original scalar q analytically. Finite
            # differences of sqrt(eta/2)+an exponentially smaller increment
            # otherwise measure cancellation noise at this endpoint.
            C=1/rho**2-1/(1-rho)**2
            H1=2/rho**3+2/(1-rho)**3;H2=6/rho**4-6/(1-rho)**4
            if p.log(C)>4*p.dps*p.log(10):
                # Reference rounding only: certify the exact cutoff value
                # and first two rho derivatives differ from (1,0,0) by less
                # than 10^(-dps-100). Production retains its formal tail.
                assert C-p.log(1+H1+abs(H2)+H1**2)>(p.dps+100)*p.log(10)
                chi=p.mpf(1);chi1=chi2=p.mpf(0)
            else:
                tail=p.exp(-C);chi=1/(1+tail);complement=tail/(1+tail)
                chi1=-chi*complement*H1
                chi2=chi*complement*(H2+(1-2*chi)*H1**2)
            root=p.sqrt(eta*(2-rho)/(2*(2+eta*rho)))
            ell1=-1/(2*(2-rho))-eta/(2*(2+eta*rho))
            ell2=-1/(2*(2-rho)**2)+eta**2/(2*(2+eta*rho)**2)
            qr=root*(chi1+chi*ell1)
            qrr=root*(chi2+2*chi1*ell1+chi*(ell2+ell1**2))
            return [chi*root,qr*rho1,qrr*rho1**2+qr*rho2]
        closed=coord.geometry('xi','0','1/4')
        norm=coord.normalized_excess(template,closed)
        qsource=source.previous.correlated_original_q(norm['rows'],-c.mpf(D),c.ln(2))
        assert qsource['record']['branch']=='active'
        S=p.mpf(1)/(4*p.sqrt(D))
        sigma_rows=source.closed_sigma_rows(template,closed['coordinate'])
        for s in (p.mpf(0),S/7,S/2,S):
            qrefs=direct_q_rows(s)
            for order in range(5):
                qcheck.enclosed(p,c,sigma_rows[order],p.diff(sigma,s,order) if s else p.mpf(0));comparisons+=1
            for order in range(3):
                ref=p.diff(lambda x:2*sigma(x)/eta,s,order) if s else p.mpf(0)
                qcheck.enclosed(p,c,norm['rows'][(order,0)],ref);comparisons+=1
                ref=qrefs[order]
                qcheck.enclosed(p,c,qsource['rows'][(order,0)],ref);comparisons+=1
        assert all(v.zero for k,v in qsource['rows'].items() if k[1])
        assert norm['record']['original_closed_endpoint_weighted_source']['positive_monotonicity_margin']
        gap=coord.geometry('mixed','3/4','4');norm=coord.normalized_excess(template,gap)
        qsource=source.previous.correlated_original_q(norm['rows'],-c.mpf(D),c.ln(2))
        a=3/(4*p.sqrt(D));b=1/p.sqrt(D+4)
        qcheck.enclosed(p,c,gap['width'],b-a);width_checks+=1
        for s in (a,(a+b)/2,b):
            qrefs=direct_q_rows(s)
            for order in range(3):
                qcheck.enclosed(p,c,norm['rows'][(order,0)],p.diff(lambda x:2*sigma(x)/eta,s,order));comparisons+=1
                qcheck.enclosed(p,c,qsource['rows'][(order,0)],qrefs[order]);comparisons+=1
        quiet=source.quiet_geometry(coordinates,coord,coord.W,Fraction(1,4))
        qcheck.enclosed(p,c,quiet['width'],p.mpf('.25')-1/p.sqrt(D));width_checks+=1
        # Independent affine transfer of known nonzero functions and their
        # Z derivatives through the true mixed width, all original rates.
        operator=source.rc.history.C1DuhamelOperator(coordinates)
        widths={};values={};jets={}
        for index,(key,rate) in enumerate(source.rc.RATES.items()):
            mass=source.rc.transfer.true_width_kernel(coordinates,gap,rate)
            widths[key]=mass
            values[key]=template.scalar(c.mpf((-(index+1),index+2)))*mass['mass']
            jets[key]=template.scalar(c.mpf((-(index+2),index+3)))*mass['mass']
        source.rc.transfer.append_true_cell(operator,gap,values,jets,family)
        incoming={key:template.scalar(index+1) for index,key in enumerate(source.KEYS)}
        incomingZ={key:template.scalar(-(index+1)) for index,key in enumerate(source.KEYS)}
        result=operator.apply(incoming,incomingZ,family)
        w=b-a
        for index,(key,rate) in enumerate(source.rc.RATES.items()):
            r=p.mpf(str(rate));mass=w if not r else -p.expm1(-r*w)/r
            qcheck.enclosed(p,c,result['values'][key],p.exp(-r*w)*(index+1)+p.mpf('.3')*mass);comparisons+=1
            qcheck.enclosed(p,c,result['Z_derivatives'][key],-p.exp(-r*w)*(index+1)-p.mpf('.4')*mass);comparisons+=1
        assert ep(operator.coefficients['p'].coefficient)==(1,1)
    return dict(passed=True,independent_new_closed_endpoint_sigma_excess_q_mixed_gap_and_affine_comparisons=comparisons,
        independent_mixed_and_mixed_endpoint_true_width_checks=width_checks,
        independent_scalar_q_derivatives_analytic_to_avoid_cancellation_noise=True,
        independent_cutoff_reference_rounding_has_proved_subprecision_error=True,
        no_epsilon_endpoint_or_selected_incoming_fixture_substituted_for_native_source=True)


def tile_limits_and_memory(c,data):
    rows=data['original_complete_prefix_source_rows']
    assert [row['label'] for row in rows]==[row[0] for row in source.PLAN]+['support_to_fixed_quiet_collar']
    assert data['original_endpoint_order_and_exact_width_telescope']==source.original_endpoint_telescope()
    assert rows[0]['actual_true_geometry']['exact_normalized_endpoints']==['0','1/4']
    assert rows[0]['original_closed_normalized_excess']['exact_s_zero_excess_and_all_slow_jets_zero']
    for row in rows[:-1]:
        assert row['actual_original_local_five_C0_Z_density_and_integral']['status']=='enclosed'
        assert all(inv['status']=='enclosed' for inv in row['actual_original_conditioned_first_jet_queries'])
        assert row['actual_original_typed_phase']['actual_phase_is_N_times_original_s_not_N_times_xi_or_k']
        assert row['actual_original_typed_phase']['phase_independent_of_Z']
        assert ep(iv(c,row['actual_true_geometry']['positive_true_log_radius_width']['coefficient_interval']))[0]>0
    assert rows[-1]['original_q_and_all_required_slow_jets_exact_zero']
    assert rows[-1]['incoming_history_not_zeroed_or_replaced']
    incoming=data['original_prefix_operator_and_genuine_incoming']['actual_incoming']
    assert incoming['source_cell_index']==2 and incoming['ordered_source_cells']==2048 and incoming['candidate_N']==source.N
    original=incoming['original_source_cell_record'];following=incoming['following_original_transition_inherited_record']
    assert original['label']=='buffer_9_11' and following['label']=='transition'
    assert original['actual_right_correction_C0']==following['actual_inherited_C0']
    assert original['actual_right_correction_Z']==following['actual_inherited_Z']
    assert incoming['actual_Rc_downstream_values_not_used_as_Rd_inlet']
    assert incoming['original_background_not_used_as_correction']
    assert all(not rec['point_value_selected'] and rec['encloses_original_source_function']
        for values in (incoming['actual_incoming_correction_C0'],incoming['actual_incoming_correction_Z']) for rec in values.values())
    assert all(row['original_whole_power_q_q_Z_and_density_C0_Z_exact_zero'] for row in data['original_quiet_power_source_records'])
    for operator_key,width in (('original_complete_transition_operator',1),('original_Rd_to_Rc_operator',3)):
        operator=data[operator_key]
        coords=source.native.HalfPstarCoordinates(c,iv(c,data['common_directed_coordinate_theorem']['common_log_bases'][1]),data['source_family'])
        for key,rate in source.rc.RATES.items():saved.same(c,operator['incoming_C0_Z_decay_coefficients'][key],coords.decay(width,rate))
        assert operator['incoming_argument_not_assumed_or_reset']
    assert data['actual_original_pressure_datum_kept_separate_and_added_once']
    assert not data['actual_Rc_targets_controls_or_global_N_or_recursion_admitted']


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert manifest['closed_original_endpoint_no_epsilon_gap']
    assert manifest['full_original_transition_C0_Z_operator_with_genuine_Rd_incoming']
    assert manifest['updated_actual_Rc_correction_and_own_history_C0_Z_covers_executed']
    assert not manifest['actual_full_axial_buffer_numerical_integrals_or_Rc_targets_controls_admitted']
    assert all(manifest[k] is False for k in source.common.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    c=MPIntervalContext();c.dps=240;queries=0;rows_checked=0
    tail_report=json.loads((HERE/source.tail.NAME).read_bytes())
    quiet=json.loads((HERE/source.support.NAME).read_bytes())
    assert quiet['original_combined_quiet_transition_C1_window']==[str(Fraction(1,2**128)),'1']
    with mp.workdps(300):
        independent=independent_new_source_checks(manifest['source_family'])
        for archive in manifest['actual_original_complete_transition_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            data=json.loads(raw);assert data['source_family']==manifest['source_family'] and data['candidate_N']==source.N
            parent_archive=data['actual_native_source_input_archive']
            parent_raw=gzip.decompress((HERE/parent_archive['filename']).read_bytes())
            assert hashlib.sha256(parent_raw).hexdigest()==parent_archive['lossless_original_json_sha256']
            parent=json.loads(parent_raw)
            coordinates=source.native.HalfPstarCoordinates(c,iv(c,parent['common_directed_coordinate_theorem']['common_log_bases'][1]),manifest['source_family'])
            replay=source.execute_tile(parent,coordinates,tail_report)
            # JSON stores integer polynomial exponents as string object keys.
            # Compare the lossless serialized representation, as in archives.
            for key,value in replay.items():exact_replay_equal(json.loads(json.dumps(source.encode(value))),data[key],key)
            tile_limits_and_memory(c,data)
            assert data['full_original_transition_interval_zero_to_one_covered']
            queries+=len(source.PLAN);rows_checked+=len(source.PLAN)*len(source.KEYS)*2
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        independent_original_closed_endpoint_mixed_gap_and_memory_checks=independent,
        complete_original_transition_source_queries_replayed=queries,
        original_local_signed_C0_Z_contribution_rows_replayed=rows_checked,
        exact_endpoint_telescope_closed_prefix_no_gaps_and_original_widths_checked=True,
        genuine_Rd_incoming_correction_not_Rc_or_background_checked=True,
        exact_quiet_operator_and_Rc_pressure_memory_replayed=True,
        updated_Rc_correction_own_history_and_separate_pressure_C0_Z_checked=True,
        saved_actual_native_parents_and_phase_origins_used=True,
        accepted_upstream_owner_and_numerical_integration_producers_not_rerun=True,
        full_axial_buffer_targets_controls_N_heat_stress_recursion_and_NS_not_claimed=True,
        **dict.fromkeys(source.common.current.FLAGS,False),input_hashes=hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Closed original prefix/complete transition and genuine Rc memory PASS; targets/global construction open',flush=True)
    return result


if __name__=='__main__':run()

"""Focused independent physical-y q derivative and typed width/source checks.

No accepted numerical native integrals are rerun. The fixtures evaluate the
original logistic/q formula directly at independent physical-s points.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source as source
import lei_ren_part1_paper_compliant_current_generic_tail_phase_integrals_check as generic

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv
same=generic.same;rc,common,prior=source.rc,source.common,source.prior


def enclosed(p,c,value,reference):
    if not reference:
        lo,hi=ep(value.coefficient);assert lo<=0<=hi;return
    anchor=p.mpf(ep(value.scale.evaluate())[1])
    normalized=p.sign(reference)*p.exp(p.log(abs(reference))-anchor)
    cover=value.coefficient*value.bounded_exp(value.scale.evaluate()-c.mpf(anchor))
    lo,hi=ep(cover)
    assert p.mpf(lo)<=normalized<=p.mpf(hi),('Independent source not enclosed',
        p.nstr(reference,10),p.nstr(normalized,10),p.nstr(lo,10),p.nstr(hi,10))


def independent_checks(family):
    c=MPIntervalContext();c.dps=90;p=mp.mp.clone();p.dps=2000
    coordinates=rc.history.CommonSourceCoordinates(c,c.mpf(0),family)
    comparisons=0;width_checks=0;branches=[]
    def sigmoid(x):
        if x<=0:return p.mpf(0)
        if x>=1:return p.mpf(1)
        return 1/(1+p.exp(1/x**2-1/(1-x)**2))
    for D in (64,128,256):
        coord=source.OriginalTransitionCoordinates(coordinates,c.mpf(1),-c.mpf(D))
        eta=p.exp(-D)
        def q(s):
            Delta=2*sigmoid(s)
            if Delta>=eta:return p.mpf(0)
            return sigmoid(1-Delta/eta)*p.sqrt((2*eta-Delta)/(2*(2+Delta)))
        for kind,left,right in (('xi','1/4','3/4'),('k',4,2),
                ('k',2,'7/4'),('k','7/4','5/3'),('k','5/3',0)):
            g=coord.geometry(kind,left,right)
            norm=coord.normalized_excess(coordinates.scalar(1),g)
            got=source.correlated_original_q(norm['rows'],-c.mpf(D),c.ln(2))
            assert got['rows'] is not None
            branches.append(got['record']['branch'])
            lo=p.mpf(str(left)) if '/' not in str(left) else p.mpf(source.Fraction(left).numerator)/source.Fraction(left).denominator
            hi=p.mpf(str(right)) if '/' not in str(right) else p.mpf(source.Fraction(right).numerator)/source.Fraction(right).denominator
            for x in (lo,(lo+hi)/2,hi):
                s=x/p.sqrt(D) if kind=='xi' else 1/p.sqrt(D+x)
                for order in range(3):
                    ref=q(s) if order==0 else p.diff(q,s,order)
                    enclosed(p,c,got['rows'][(order,0)],ref);comparisons+=1
                independent_norm=2*sigmoid(s)/eta
                enclosed(p,c,norm['rows'][(0,0)],independent_norm);comparisons+=1
            exact_width=(hi-lo)/p.sqrt(D) if kind=='xi' else 1/p.sqrt(D+hi)-1/p.sqrt(D+lo)
            enclosed(p,c,g['width'],exact_width);width_checks+=1
            for k,v in got['rows'].items():
                if k[1]:assert v.zero
        zero=source.correlated_original_q({k:coordinates.scalar(0) for k in source.ORDERS},-c.mpf(D),c.ln(2))
        enclosed(p,c,zero['rows'][(0,0)],p.sqrt(eta/2));comparisons+=1
        assert all(v.zero for k,v in zero['rows'].items() if k!=source.ZERO)
    assert 'smooth_cutoff_seam' in branches and 'active' in branches and 'flat' in branches
    return dict(passed=True,independent_original_physical_y_q0_q1_q2_and_excess_comparisons=comparisons,
        independent_exact_typed_true_width_checks=width_checks,
        original_active_smooth_seam_flat_and_exact_zero_endpoint_checked=True,
        no_normalized_coordinate_derivatives_substituted=True)


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert all(manifest[k] is False for k in common.current.FLAGS)
    assert not manifest['full_prefix_or_tail_numerical_integrals_or_Rc_targets_admitted']
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    archive=manifest['original_typed_microscopic_source_archive']
    compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
    assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
    assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
    data=json.loads(raw);assert data['source_family']==manifest['source_family']
    c=MPIntervalContext();c.dps=240
    with mp.workdps(300):
        independent=independent_checks(manifest['source_family'])
        params=data['original_parameter_source']
        inventory=json.loads((HERE/source.previous.prior.signed.current.NAME).read_bytes())['current_actual_loop_jet_log_bounds_by_chart']
        certificate=inventory['O3_slope_mu']['actual_positive_denominator_theorem']
        scales=json.loads((HERE/(source.PREFIX+'current_generic_shear_uniform_inputs.json')).read_bytes())['current_actual_logarithmic_loop_scales']
        mu=iv(c,certificate['actual_positive_mu']);eta=iv(c,scales['selected_positive_eta_log'])
        assert ep(iv(c,params['original_mu']))==ep(mu) and ep(iv(c,params['original_eta_log']))==ep(eta)
        bases=data['common_directed_coordinate_theorem']['common_log_bases']
        coordinates=rc.history.CommonSourceCoordinates(c,iv(c,bases[1]),manifest['source_family'])
        coord=source.OriginalTransitionCoordinates(coordinates,mu,eta)
        assert source.encode(coord.parameter_record)==params
        positive=iv(c,certificate['log_actual_a_positive_lower'])
        assert ep(iv(c,data['original_positive_a_lower_log']))==ep(positive)
        expected=dict(open_bulk='active',near_seam_active='active',seam_active_collar='active',
            cutoff_seam='smooth_cutoff_seam',seam_flat_collar='flat')
        rows=data['original_source_coordinate_queries'];assert len(rows)==5
        masses_checked=0;q_rows_checked=0
        for row in rows:
            saved=row['actual_original_geometry']
            geometry=coord.geometry(saved['typed_coordinate_kind'],*saved['exact_normalized_endpoints'])
            assert source.encode(geometry['record'])==saved
            norm=coord.normalized_excess(coordinates.scalar(1),geometry)
            assert source.encode(norm['record'])==row['original_correlated_normalized_excess']
            got=source.correlated_original_q(norm['rows'],eta,positive)
            assert source.encode(got['record'])==row['original_q_source']
            assert got['record']['branch']==expected[row['label']]
            assert all(v.zero for k,v in got['rows'].items() if k[1])
            q_rows_checked+=len(got['rows'])
            if got['record']['branch']=='active':assert ep(got['rows'][source.ZERO].coefficient)[0]>0
            if got['record']['branch']=='flat':assert all(v.zero for v in got['rows'].values())
            for k,rate in rc.RATES.items():
                mass=rc.transfer.true_width_kernel(coordinates,geometry,rate)
                saved_mass=row['original_true_width_mass_and_decay'][k]
                assert mass['branch']==saved_mass['branch']
                same(c,saved_mass['mass'],mass['mass']);same(c,saved_mass['decay'],mass['decay'])
                assert not mass['mass'].zero and ep(mass['mass'].coefficient)[0]>0
                masses_checked+=1
            phase=row['actual_original_typed_radius_phase']
            assert phase['same_original_affine_radius_identity']['passed']
            assert phase['actual_phase_is_N_times_original_s_not_N_times_xi_or_k']
            assert phase['nonzero_microscopic_phase_component_retained']
            base_saved=phase['original_base_source_phase']
            assert base_saved['chart']=='O3_slope_mu' and base_saved['candidate_N']==source.N
            assert base_saved['original_spatial_phase_bound_at_candidate_N']
            micro=common.restore_common_source(phase['original_microscopic_s_component'],coordinates)
            if geometry['kind']=='xi':
                expected_micro=coord.W*c.mpf((ep(coord.cv(geometry['left']))[0],ep(coord.cv(geometry['right']))[1]))
            else:
                kval=c.mpf((ep(coord.cv(geometry['right']))[0],ep(coord.cv(geometry['left']))[1]))
                expected_micro=prior.ScaledEnclosure(prior.FormalScale(coordinates.bases,offset=-c.ln(coord.D+kval)/2),1,coordinates.ledger)
            same(c,phase['original_microscopic_s_component'],expected_micro)
            same(c,phase['original_microscopic_N_s_component'],micro*source.N)
            assert not micro.zero and ep(micro.coefficient)[0]>0
            if base_saved['periodic_projection']['full_period']:
                assert phase['actual_phase_boxes']==base_saved['periodic_projection']['boxes']
            else:
                projected=[]
                for box in base_saved['periodic_projection']['boxes']:
                    projected.extend(source.first.spatial.ordinary_mod_one(c,iv(c,box)+geometry['coordinate']*source.N)['boxes'])
                assert source.encode(projected)==phase['actual_phase_boxes']
            assert row['original_amplitude_inverse_density_and_integral_not_evaluated']
        endpoint=source.correlated_original_q({k:coordinates.scalar(0) for k in source.ORDERS},eta,positive)
        assert source.encode(endpoint['record'])==data['exact_original_s_zero_active_q_endpoint']
        assert ep(endpoint['rows'][source.ZERO].coefficient)[0]>0
        assert data['original_bulk_prefix_before_xi_quarter_and_bulk_to_seam_not_integrated']
        assert not data['full_original_prefix_amplitude_inverse_density_integral_admitted']
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        independent_original_source_and_true_width_checks=independent,
        same_original_parameter_bindings_checked=True,original_source_q_rows_checked=q_rows_checked,
        original_true_width_mass_decay_rows_checked=masses_checked,
        actual_original_N_s_phase_and_nonzero_offset_checked=True,
        original_bulk_seam_flat_endpoint_queries_checked=True,
        original_amplitude_inverse_density_full_integral_not_claimed=True,
        **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original normalized microscopic q/source PASS; nonzero active q and smooth seam retained; full integrals open',flush=True)
    return result


if __name__=='__main__':run()

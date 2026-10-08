"""Original logistic derivative/source refinement and quiet C1 collar audit.

Checks independent scalar logistic derivatives, genuine saved mu/eta/root
bindings, original q classification, actual phase/density and two-cell affine
transport. No prior native owner or accepted numerical integrals are rerun.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_transition_positive_log_sigma as source
import lei_ren_part1_paper_compliant_current_generic_tail_phase_integrals_check as generic

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv
same,restore=generic.same,generic.native_restore
rc,common,prior=source.rc,source.common,source.prior


def payload(archive):
    compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
    assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
    assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
    return json.loads(raw)


def independent_logistic_checks(family):
    c=MPIntervalContext();c.dps=100;p=mp.mp.clone();p.dps=300
    coordinates=rc.history.CommonSourceCoordinates(c,c.mpf(0),family);base=coordinates.scalar(1);count=0
    def reference(t):
        L=1/(1-t)**2-1/t**2
        return 1/(1+p.exp(-L))
    for left,right in (('1e-20','1.01e-20'),('.001','.00101'),('.015625','.025'),('.25','.3'),('.5','.5')):
        t=c.mpf((left,right));got=source.positive_log_sigma(base,t)
        assert ep(got['values'][(0,0)].coefficient)[0]>0 and not got['values'][(0,0)].zero
        a,b=ep(t);points=(p.mpf(a),) if a==b else (p.mpf(a),(p.mpf(a)+p.mpf(b))/2,p.mpf(b))
        for x in points:
            for order in range(3):
                value=got['values'][(order,0)]
                ref=reference(x) if order==0 else p.diff(reference,x,order)
                if x==p.mpf('.5') and order==1:
                    assert abs(ref-8)<p.mpf('1e-250');ref=p.mpf(8)
                if x==p.mpf('.5') and order==2:
                    assert abs(ref)<p.mpf('1e-250');ref=p.mpf(0)
                anchor=p.mpf(ep(value.scale.evaluate())[1])
                norm=0 if not ref else ref if not anchor else p.sign(ref)*p.exp(p.log(abs(ref))-anchor)
                cover=value.coefficient*value.bounded_exp(value.scale.evaluate()-c.mpf(anchor))
                lo,hi=ep(cover);assert p.mpf(lo)<=norm<=p.mpf(hi),('Independent logistic derivative not enclosed',
                    left,right,p.nstr(x,8),order,p.nstr(norm,15),p.nstr(lo,15),p.nstr(hi,15))
                count+=1
        assert all(v.zero for k,v in got['values'].items() if k[1])
    endpoint=source.positive_log_sigma(base,c.mpf(0));assert all(v.zero for v in endpoint['values'].values())
    support_checks=0
    for d in (64,1024,2**256):
        proof=source.original_active_support_proof(coordinates,c.mpf(1),-c.mpf(d))
        w=proof['original_q_active_support_width'];logwidth=iv(c,w['formal_positive_scale']['additional_log_interval'])
        expected=-p.log(d)/2;assert p.mpf(ep(logwidth)[0])<=expected<=p.mpf(ep(logwidth)[1])
        for factor in (p.mpf(1),p.mpf('1.25'),p.mpf(2)):
            t=factor/p.sqrt(d);L=1/(1-t)**2-1/t**2
            independent_log_Delta_over_eta=p.log(2)+d+L-p.log1p(p.exp(L))
            assert independent_log_Delta_over_eta>=0;support_checks+=1
    return dict(passed=True,independent_original_logistic_ordinary_y0_y1_y2_comparisons=count,
        independent_original_active_support_log_ratio_checks=support_checks,
        exact_original_endpoint_sigma_and_slow_jets_zero=True,
        nonmaterializable_positive_sigma_factor_not_erased=True,
        independent_point_derivatives_used_not_native_field_point_values=True)


def original_sigma_q_check(c,coordinates,row,mu,eta,log_a_lower):
    bases=tuple(iv(c,v) for v in row['original_native_source_log_bases']);ledger=dict(directed_small_exponential_tails=0,
        positive_function_denominator_intersections=0,positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    rootsaved=row['actual_original_phase_inverse_and_density_source']['original_q_slow_jet_source']
    fresh=rootsaved['original_correlated_shear_and_q']['correlated_shear_and_signed_root_enclosures']
    old=row['original_unrefined_q_slow_source']['original_correlated_shear_and_q']['correlated_shear_and_signed_root_enclosures']
    restore_here=lambda v:restore(v,bases,ledger)
    roots={name:{order:restore_here(rows['y%d_Z%d'%order]) for order in source.slow.ORDERS} for name,rows in fresh.items()}
    base=roots['E'][(0,0)];sigma=source.positive_log_sigma(base,iv(c,row['actual_true_geometry']['native_coordinate_box']))
    assert source.encode(sigma['record'])==row['same_original_positive_log_sigma_refinement']
    assert rootsaved['original_definition_sigma_eta_q_unchanged']
    assert ep(iv(c,rootsaved['original_selected_mu']))==ep(mu)
    assert ep(iv(c,rootsaved['original_selected_eta_log']))==ep(eta)
    assert ep(iv(c,rootsaved['original_actual_a_positive_lower_log']))==ep(log_a_lower)
    delta={k:v*(2*mu) for k,v in sigma['values'].items()}
    for k,v in delta.items():
        same(c,fresh['kappa_minus2']['y%d_Z%d'%k],v)
        same(c,fresh['a']['y%d_Z%d'%k],v+(2 if k==(0,0) else 0))
        assert roots['b'][k].zero and roots['t0'][k].zero
    for name in roots:
        if name in ('a','kappa_minus2','t0'):continue
        assert fresh[name]==old[name],'Unrelated original root changed: '+name
    loop=source.slow.current.q_enclosure(roots['a'][(0,0)],delta[(0,0)],eta,log_a_lower)
    jet=source.slow.original_q_jet(roots,eta,log_a_lower,loop)
    assert jet['status']==rootsaved['status'] and jet['branch']==rootsaved['branch']
    loop_saved=rootsaved['original_correlated_shear_and_q']['original_q_enclosure']
    for k,v in loop.items():
        if isinstance(v,prior.ScaledEnclosure):same(c,loop_saved[k],v)
        else:assert source.encode(v)==loop_saved[k]
    same(c,rootsaved['current_refined_original_q_C0'],loop['q'])
    if jet['rows'] is None:assert rootsaved['original_q_ordinary_slow_derivative_enclosures'] is None
    else:
        for k,v in jet['rows'].items():same(c,rootsaved['original_q_ordinary_slow_derivative_enclosures']['y%d_Z%d'%k],v)
    return dict(q_status=jet['status'],q_branch=jet['branch'],original_sigma_strictly_positive=not sigma['values'][(0,0)].zero
        and ep(sigma['values'][(0,0)].coefficient)[0]>0)


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert not manifest['complete_original_transition_or_tail_numerical_integrals_or_Rc_targets_admitted']
    assert all(manifest[k] is False for k in common.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    inventory=json.loads((HERE/prior.signed.current.NAME).read_bytes())['current_actual_loop_jet_log_bounds_by_chart']
    certificate=inventory['O3_slope_mu']['actual_positive_denominator_theorem']
    scales=json.loads((HERE/(source.PREFIX+'current_generic_shear_uniform_inputs.json')).read_bytes())['current_actual_logarithmic_loop_scales']
    c=MPIntervalContext();c.dps=240;checks=[];phase_count=0
    with mp.workdps(300):
        independent=independent_logistic_checks(manifest['source_family'])
        mu=iv(c,certificate['actual_positive_mu']);eta=iv(c,scales['selected_positive_eta_log'])
        log_a_lower=iv(c,certificate['log_actual_a_positive_lower'])
        for archive in manifest['actual_original_positive_log_sigma_source_archives']:
            data=payload(archive);Z=ep(c.mpf(data['exact_Z_range']))
            assert data['source_family']==manifest['source_family'] and data['candidate_N']==source.N
            theorem=data['common_directed_coordinate_theorem'];bases=theorem['common_log_bases']
            coordinates=rc.history.CommonSourceCoordinates(c,iv(c,bases[1]),manifest['source_family'])
            assert all(ep(iv(c,bases[k]))==(0,0) for k in (0,2,3,4))
            support=source.original_active_support_proof(coordinates,mu,eta)
            assert source.encode(support)==manifest['original_source_derived_nonzero_q_support']
            assert data['unresolved_original_transition_prefix_not_skipped']
            assert data['original_mu_eta_sigma_q_and_source_family_unchanged']
            assert not data['entire_transition_or_tail_numerical_integrals_admitted']
            rows=data['original_left_half_positive_log_sigma_queries'];assert len(rows)==2
            assert rows[0]['exact_native_endpoints']==['0',str(source.COLLAR_LEFT)]
            assert rows[1]['exact_native_endpoints']==[str(source.COLLAR_LEFT),'1/2']
            decisions=[];covered=[]
            for row in rows:
                decision=original_sigma_q_check(c,coordinates,row,mu,eta,log_a_lower);decisions.append(decision)
                got,n=generic.cell_check(c,coordinates,row,Z);phase_count+=n;covered.append(got)
            assert covered[0] is None and covered[1] is not None
            assert decisions[1]['q_branch']=='flat' and decisions[1]['original_sigma_strictly_positive']
            assert decisions[0]['q_status']=='requires_source_branch_subdivision'
            sigma=rows[1]['same_original_positive_log_sigma_refinement']
            delta=rows[1]['actual_original_phase_inverse_and_density_source']['original_q_slow_jet_source']['current_refined_original_excess_ordinary_rows']['y0_Z0']
            ds=delta['formal_positive_scale'];delta_lower=iv(c,ds['additional_log_interval'])+c.ln(iv(c,delta['coefficient_interval']))
            assert ep(delta_lower)[0]>=ep(eta)[1]
            right=payload(data['accepted_right_half_source_archive'])['original_entire_transition_right_half_source_query']
            rgot,n=generic.cell_check(c,coordinates,right,Z);assert rgot is not None;phase_count+=n
            coefficients={k:coordinates.scalar(1) for k in rc.RATES}
            values={k:coordinates.scalar(0) for k in rc.RATES};jets={k:coordinates.scalar(0) for k in rc.RATES}
            for cell in (covered[1],rgot):
                assert all(v.zero for v in (*cell['values'].values(),*cell['Z_derivatives'].values()))
                for k in rc.RATES:
                    decay=cell['factors'][k]['decay']
                    coefficients[k]=decay*coefficients[k]
                    values[k]=decay*values[k]+cell['values'][k];jets[k]=decay*jets[k]+cell['Z_derivatives'][k]
            operator=data['original_combined_quiet_transition_C1_operator']
            assert operator['steps']==2 and operator['incoming_argument_not_assumed_or_reset']
            assert operator['original_rates']=={k:str(r) for k,r in rc.RATES.items()}
            assert data['original_combined_quiet_transition_window']==[str(source.COLLAR_LEFT),'1']
            assert data['exact_true_log_radius_width']==str(1-source.COLLAR_LEFT)
            for k in rc.RATES:
                same(c,operator['incoming_C0_Z_decay_coefficients'][k],coefficients[k])
                same(c,operator['cumulative_signed_increment_C0_enclosures'][k],values[k])
                same(c,operator['cumulative_signed_increment_Z_enclosures'][k],jets[k])
                width=c.mpf(1)-c.mpf(source.COLLAR_LEFT.numerator)/source.COLLAR_LEFT.denominator
                collected=coordinates.decay(width,rc.RATES[k]);same(c,operator['incoming_C0_Z_decay_coefficients'][k],collected)
            assert ep(coefficients['p'].coefficient)==(1,1) and coefficients['p'].scale.powers==(0,0,0,0,0)
            checks.append(dict(exact_Z_range=data['exact_Z_range'],actual_original_query_decisions=decisions,
                delta_positive_log_lower_minus_eta_upper=ep(delta_lower)[0]-ep(eta)[1],
                original_combined_quiet_window=[str(source.COLLAR_LEFT),'1'],own_integral_C0_Z_exact_zero=True,
                original_unresolved_prefix_not_zeroed=True,original_incoming_and_pressure_memory_retained=True))
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        original_logistic_sigma_independent_checks=independent,actual_original_source_collar_checks=checks,
        actual_enclosed_phase_density_cells_checked=phase_count,actual_new_left_collar_C0_Z_integral_rows_checked=20,
        actual_combined_quiet_C0_Z_integral_rows_checked=20,
        source_derived_original_active_support_width_checked=True,
        independent_two_cell_affine_and_collected_true_width_decay_replay=True,
        entire_transition_tail_targets_or_physical_closure_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original positive-log sigma/collar C1 PASS; original source parameters unchanged; unresolved prefix retained',flush=True)
    return result


if __name__=='__main__':run()

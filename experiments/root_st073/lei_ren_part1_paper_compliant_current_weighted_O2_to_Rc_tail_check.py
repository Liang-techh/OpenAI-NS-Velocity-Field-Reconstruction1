"""Focused genuine source, true mass and complete injected tail-memory audit.

Saved native source covers are reconstructed; no native owner, original
upstream route or accepted O2 numerical integration is recomputed.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_weighted_O2_to_Rc_tail as source

common=source.common;HERE=source.HERE;sha=source.sha;ep=source.ep;iv=source.iv
same=source.inlet.same;native_restore=source.inlet.saved.native_restore


def source_check(c,coordinates,row,scales,Z):
    signed=row['original_full_box_signed_source_and_old_q_slow_status'];prov=signed['source_provenance']
    assert signed['chart']==row['chart']==prov['chart']
    assert prov['source_family']==coordinates.family and ep(iv(c,prov['Z_box']))==Z
    assert prov['native_provider_group'] in ('pre','O3') and 'amplitude_adapter' not in prov
    assert prov['one_original_seed_graph_asserted'] and not prov['cache_cover']
    assert prov['selector_derivatives_not_substituted'] and prov['numerical_caps_not_used_as_defining_function_values']
    assert prov['derivative_coordinate']=='ordinary y=log R'
    g=row['actual_true_cell_geometry']
    assert ep(iv(c,prov['coordinate_box']))==ep(iv(c,g['native_coordinate_box']))
    bases=(c.mpf(0),coordinates.logP_squared,c.mpf(0),c.mpf(0),iv(c,prov['logR_cover']))
    restore=lambda v:native_restore(v,bases,coordinates.ledger)
    roots={name:{(int(key[1]),int(key[-1])):restore(value) for key,value in rows.items()}
        for name,rows in signed['original_correlated_shear_and_q']['correlated_shear_and_signed_root_enclosures'].items()}
    positive=row['original_chart_uniform_a_positive_certificate']
    got=source.serial.whole_period_C1(roots,iv(c,scales['selected_positive_eta_log']),
        iv(c,positive['log_actual_a_positive_lower']),iv(c,scales['logarithmic_selected_positive_lower_constants']['d_star']))
    proof=row['original_periodic_parameter_C1_cover']
    assert source.encode(got['record'])==proof
    geometry=dict(width=common.restore_common_source(g['positive_true_log_radius_width'],coordinates),
        regular=iv(c,g['regular_true_log_radius_width_cover']),
        scalar_cover=iv(c,g['scalar_width_cover_used_only_for_directed_kernel_bounds']))
    assert g['width_and_endpoints_independent_of_Z'] and ep(geometry['width'].coefficient)[0]>0
    factors={};adds={};jets={}
    for k,rate in source.rc.RATES.items():
        factor=source.rc.transfer.true_width_kernel(coordinates,geometry,rate)
        saved=row['original_true_width_kernel_factors'][k]
        assert saved['branch']==factor['branch']
        same(c,saved['mass'],factor['mass']);same(c,saved['decay'],factor['decay'])
        rho=coordinates.rebase(restore(row['actual_signed_density_C0_covers'][k]),coordinates.family)
        rhoZ=coordinates.rebase(restore(row['actual_signed_density_Z_covers'][k]),coordinates.family)
        adds[k]=rho*factor['mass'];jets[k]=rhoZ*factor['mass'];factors[k]=factor
        same(c,row['actual_true_width_C0_contributions'][k],adds[k])
        same(c,row['actual_true_width_Z_contributions'][k],jets[k])
    if row['chart']=='O3_power':
        assert row['original_whole_power_q_q_Z_and_density_C0_Z_exact_zero']
        assert proof['original_full_box_q_and_q_Z_exact_zero_implies_primitive_C0_Z_exact_zero']
        assert all(v.zero for v in (*adds.values(),*jets.values()))
        assert roots['b'][(0,0)].zero and roots['b'][(0,1)].zero
        assert roots['a'][(0,1)].zero and roots['kappa_minus2'][(0,1)].zero
        assert positive['exact_a_minus2_source']=='2*mu'
        same(c,source.encode(roots['kappa_minus2'][(0,0)].record()),roots['kappa_minus2'][(0,0)].scalar(2*iv(c,positive['actual_positive_mu'])))
        assert row['quiet_local_density_does_not_reset_inherited_histories']
    assert factors['p']['decay'].scale.powers==(0,0,0,0,0)
    assert ep(factors['p']['decay'].coefficient)==(1,1)
    return geometry,adds,jets,factors


def run():
    begin=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert all(manifest[k] is False for k in common.current.FLAGS)
    assert not manifest['numerical_original_source_oracle_installed']
    assert not manifest['actual_original_numerical_tail_integrals_evaluated']
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    accepted=json.loads((HERE/source.inlet.NAME).read_bytes())
    scales=json.loads((HERE/(source.PREFIX+'current_generic_shear_uniform_inputs.json')).read_bytes())['current_actual_logarithmic_loop_scales']
    c=MPIntervalContext();c.dps=240;checks=[];quiet=0;serial_rows=0
    with mp.workdps(300):
        for archive in manifest['original_six_cell_tail_source_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            payload=json.loads(raw);rows=payload['ordered_original_six_cell_source_records']
            assert payload['source_requeried_on_actual_axial_tile_not_relabelled_whole_Z']
            assert tuple(r['original_tail_label'] for r in rows)==source.LABELS
            primary=payload['genuine_primary_incoming_O2_record']
            assert primary==next(r for r in accepted['genuine_original_O2_weighted_pressure_replays']
                if source.inlet.same_Z(r['exact_Z_range'],payload['exact_Z_range']) and r['ordered_source_cells']==2048)
            bg=rows[-1]['original_right_background_and_separate_P0_Z']
            bases=bg['common_directed_coordinate_theorem']['common_log_bases']
            assert all(ep(iv(c,bases[i]))==(0,0) for i in (0,2,3,4))
            coordinates=source.rc.history.CommonSourceCoordinates(c,iv(c,bases[1]),manifest['source_family'])
            Z=ep(c.mpf(payload['exact_Z_range']));original=source.inlet_functions(primary,coordinates)
            incoming=original;composite=source.rc.history.C1DuhamelOperator(coordinates)
            geometry_rows=[]
            for row in rows:
                assert row['source_family']==manifest['source_family'] and row['candidate_N']==source.N
                assert ep(iv(c,row['Z_box']))==Z
                geometry,adds,jets,factors=source_check(c,coordinates,row,scales,Z)
                if row['chart']=='O3_power':quiet+=10
                output={};outputZ={}
                for k in source.KEYS:
                    same(c,row['actual_inherited_correction_C0'][k],incoming['values'][k])
                    same(c,row['actual_inherited_correction_Z'][k],incoming['Z_derivatives'][k])
                    output[k]=factors[k]['decay']*incoming['values'][k]+adds[k]
                    outputZ[k]=factors[k]['decay']*incoming['Z_derivatives'][k]+jets[k]
                    same(c,row['actual_right_correction_C0'][k],output[k]);same(c,row['actual_right_correction_Z'][k],outputZ[k])
                    background=row['original_right_background_and_separate_P0_Z']
                    assert background['P0_not_merged_into_pressure_history']
                    C0=common.restore_common_source(background['original_normalized_history_C0_enclosures'][k],coordinates)
                    DZ=common.restore_common_source(background['original_normalized_history_Z_enclosures'][k],coordinates)
                    same(c,row['actual_right_own_history_C0'][k],C0+output[k])
                    same(c,row['actual_right_own_history_Z'][k],DZ+outputZ[k]);serial_rows+=6
                source.rc.transfer.append_true_cell(composite,geometry,adds,jets,manifest['source_family'])
                incoming=dict(values=output,Z_derivatives=outputZ);geometry_rows.append(geometry)
            assert composite.steps==6
            # Serial and collected interval evaluation may differ in rounding;
            # each independent source cover must enclose a common result.
            collected=composite.apply(original['values'],original['Z_derivatives'],manifest['source_family'])
            for kind in ('values','Z_derivatives'):
                for k in source.KEYS:
                    a=collected[kind][k];b=incoming[kind][k]
                    # Rebase to a common dominating upper coordinate and compare.
                    ref=max(ep(a.scale.evaluate())[1],ep(b.scale.evaluate())[1])
                    aa=a.coefficient*a.bounded_exp(a.scale.evaluate()-c.mpf(ref))
                    bb=b.coefficient*b.bounded_exp(b.scale.evaluate()-c.mpf(ref))
                    assert max(ep(aa)[0],ep(bb)[0])<=min(ep(aa)[1],ep(bb)[1])
            for input_record in accepted['genuine_original_O2_weighted_pressure_replays']:
                if not source.inlet.same_Z(input_record['exact_Z_range'],payload['exact_Z_range']):continue
                replay=next(r for r in manifest['actual_refined_O2_to_Rc_replays']
                    if r['exact_Z_range']==input_record['exact_Z_range'] and r['original_O2_ordered_source_cells']==input_record['ordered_source_cells'])
                regenerated=source.replay_saved(coordinates,input_record,rows)
                regenerated['actual_six_cell_source_binding']=dict(archive=archive['filename'],sha256=sha(archive['filename']),
                    source_is_genuine_original_same_N_tile_function_cover=True)
                assert source.encode(regenerated)==replay
                p=common.restore_common_source(replay['actual_Rc_correction_C0']['p'],coordinates)
                logp=source.baseline.magnitude_log(p)
                assert logp<0
                assert replay['original_Rc_endpoint']==dict(radius_expression='Rw*exp(2)',original_power_offset=2,original_power_phase_expression='2/Tw')
                checks.append(dict(exact_Z_range=input_record['exact_Z_range'],original_O2_ordered_source_cells=input_record['ordered_source_cells'],
                    actual_six_tail_cells_checked=True,actual_Rc_pressure_correction_log_upper=logp,
                    actual_Rc_pressure_correction_magnitude_upper=c.exp(c.mpf(logp)),
                    numerical_tail_targets_or_global_N_or_terminal_closure_not_claimed=True))
            for bad in (dict(primary,candidate_N=7),dict(primary,source_family={}),dict(primary,exact_y_window=['0','.5'])):
                try:source.inlet_functions(bad,coordinates)
                except ValueError:pass
                else:raise AssertionError('Changed source/N/domain input was admitted')
            print('Refined original six-cell source and actual Rc memory PASS',payload['exact_Z_range'],flush=True)
    assert len(checks)==4 and quiet==40 and serial_rows==360
    hashes={**manifest['input_hashes'],source.NAME:sha(source.NAME),Path(__file__).name:sha(Path(__file__).name)}
    receipt=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        original_twelve_tile_cells_source_primitive_geometry_mass_audited=True,
        full_inherited_correction_and_own_C0_Z_rows_checked=serial_rows,
        original_power_exact_zero_own_density_rows_checked=quiet,
        actual_refined_Rc_replays_checked=checks,input_hashes=hashes,execution_seconds=time.monotonic()-begin,
        source_density_program_and_prior_scalar_proofs_reused_not_reexecuted=True,
        numerical_original_source_oracle_installed=False,actual_original_numerical_tail_integrals_evaluated=False,
        full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False))
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('Original refined O2 to Rc six-cell continuation PASS',flush=True)
    return receipt


if __name__=='__main__':run()

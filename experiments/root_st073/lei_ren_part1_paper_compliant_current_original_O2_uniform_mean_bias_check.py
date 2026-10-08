"""Focused all-N mean-envelope calculus, source binding and own-rate check."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_uniform_mean_bias as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep;C0,Y,Z,YZ=current.mixed.ORDERS


def exact_interval(c,left,right):
    assert ep(saved.interval(c,left))==ep(saved.interval(c,right))


def exact_native(c,left,right):
    for row in (left,right):
        assert row['encloses_original_source_function'] and not row['point_value_selected']
    assert left['exact_zero']==right['exact_zero']
    ls,rs=left['formal_positive_scale'],right['formal_positive_scale']
    assert ls['source_exponents']==rs['source_exponents'] and ls['radius_power']==rs['radius_power']
    exact_interval(c,ls['additional_log_interval'],rs['additional_log_interval'])
    exact_interval(c,left['coefficient_interval'],right['coefficient_interval'])


def independent_uniform_pair_calculus(owner):
    """Finite compatible functions test formulas, not original-scale admission."""
    frame=owner.source.frame(64,0,branch='regular');a=frame.roots['q'].atlas;c=owner.ctx
    p=mp.mp.clone();p.dps=170;cases=0;comparisons=0
    for N in (160,257,1024,10**9):
        for v in ('-.21','.21'):
            for aa in ('-.5','0','1e-50','.5'):
                data=dict(E=('1.1','-.3'),V=(v,'.4'),A=(aa,'-.13'),B=('.53','.23'))
                jets={name:current.MixedJet(a,{order:a.scalar(c.mpf(rows[1] if order==Z else rows[0] if order==C0 else 0))
                    for order in current.mixed.ORDERS}) for name,rows in data.items()}
                zz=p.mpf('-.37')
                def f(name,z):
                    value,derivative=map(p.mpf,data[name]);return value+derivative*(z-zz)
                def density(E,V):return dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
                def pair(key,z):
                    E,V,A,B=(f(name,z) for name in ('E','V','A','B'))
                    plus=density(E*p.exp(A/N),V+B/N);minus=density(E*p.exp(-A/N),V-B/N)
                    return (plus[key]+minus[key])/2-density(E,V)[key]
                targets=[{key:N*N*pair(key,zz) for key in current.integrals.KEYS},
                    {key:N*N*p.diff(lambda z:pair(key,z),zz) for key in current.integrals.KEYS}]
                for physical in (False,True):
                    coefficients=current.native_coefficients(a,**jets,physical_amplitude=physical)
                    upper=lambda group,key:ep(coefficients[group][key].finite_interval())[1]
                    values=dict(m=c.mpf(0),h=c.mpf((0,upper('C0','h'))),
                        k=c.mpf((-upper('C0','k'),upper('C0','k'))),
                        e=c.mpf((-upper('C0','e_negative'),upper('C0','e_positive'))),
                        p=c.mpf((0,upper('C0','p'))))
                    derivatives={key:c.mpf((-upper('Z',key),upper('Z',key))) for key in current.integrals.KEYS}
                    for group,target in zip((values,derivatives),targets,strict=True):
                        for key in current.integrals.KEYS:
                            saved.contains(group[key],target[key],p.mpf('1e-120'));comparisons+=1
                assert targets[0]['m']==0 and targets[1]['m']==0
                assert targets[0]['h']>=0 and targets[0]['p']>=0;cases+=1
    theorem=current.uniform_amplitude_theorem()
    assert theorem['passed'] and theorem['exact_pair_ordinary_Z_product_rule_identities']==3
    return dict(passed=True,compatible_finite_function_cases=cases,
        independent_Ns=[160,257,1024,10**9],independent_N_squared_exponential_pair_C0_Z_comparisons=comparisons,
        physical_amplitude_boundary_zero_tiny_and_both_V_signs=True,
        tested_raw_native_derivative_caps_not_large_unit_normalized_fixtures=True,
        exact_ordinary_Z_product_rule_identities=3,finite_fixtures_not_original_source_replacements=True)


def actual_source_coefficients_and_masses(owner,manifest):
    c=owner.ctx;p=mp.mp.clone();p.dps=150
    report=manifest['actual_original_O2_all_N_mean_bias'];live=owner.integrate()
    encoded=current.base.encoded(live['report']);rows=report['actual_source_mean_coefficient_records']
    assert report['source_family']==owner.family and report['all_integer_N_lower']==160
    assert report['exact_y_window']==['0','1'] and report['exact_Z_range']==['-1','1']
    assert report['original_source_cells']==64 and report['original_source_predicate_frames']==192
    assert report['original_true_phase_input_cells']==768
    assert report['original_P0_datum_sha256']==owner.family['datum_enclosure_sha256']
    assert report['no_fixed_N_mean_density_or_integral_rescaling']
    assert report['inherited_histories_P0_and_pressure_zero_rate_memory_unchanged']
    exact_native(c,report['ordinary_Z_unit'],encoded['ordinary_Z_unit'])
    nativechecks=0;inputchecks=0;rangechecks=0;masschecks=0;parts=0
    for index,(row,replayed) in enumerate(zip(rows,encoded['actual_source_mean_coefficient_records'],strict=True)):
        assert row['exact_y_cell']==[str(s.Rational(index,64)),str(s.Rational(index+1,64))]
        assert row['branch_overlap_hulled_not_added']
        branches=row['same_source_native_all_N_mean_coefficient_branches']
        assert {b['branch'] for b in branches}==set(current.integrals.BRANCHES)
        for old,new in zip(branches,replayed['same_source_native_all_N_mean_coefficient_branches'],strict=True):
            frame=owner.source.frame(64,index,branch=old['branch'])
            assert old['source_family']==owner.family and old['exact_outer_Z_bounds']==frame.record['exact_outer_Z_bounds']
            assert old['actual_original_predicate']==frame.record['actual_original_source_predicate']
            assert ep(frame.roots['a'][C0].finite_interval())[1]<=2
            oldparts=old['original_phase_coefficient_records'];assert len(oldparts)==4
            for j,(part,again) in enumerate(zip(oldparts,new['original_phase_coefficient_records'],strict=True)):
                assert part['exact_true_phase_cell']==[str(s.Rational(j,8)),str(s.Rational(j+1,8))]
                assert part['exact_full_phase_weight']=='1/4' and part['source_family']==owner.family
                for flag in ('original_q_factors_not_floored_or_dropped_before_native_products',
                        'coefficient_envelopes_are_bounds_not_selected_mean_functions',
                        'candidate_N_density_outputs_not_rescaled','N_free_source_inputs_and_uniform_exp_floor_used'):
                    assert part[flag]
                for name,values in part['original_free_phase_function_inputs'].items():
                    for order,value in values.items():
                        exact_native(c,value,again['original_free_phase_function_inputs'][name][order]);inputchecks+=1
                for oldfamily,newfamily in zip(part['native_positive_coefficient_envelopes'],again['native_positive_coefficient_envelopes'],strict=True):
                    for group,values in oldfamily.items():
                        for key,value in values.items():exact_native(c,value,newfamily[group][key]);nativechecks+=1
                for family,other in zip(part['normalized_native_and_source_amplitude_coefficient_covers'],again['normalized_native_and_source_amplitude_coefficient_covers'],strict=True):
                    for group,recalculated in zip(family,other,strict=True):
                        for key,value in group.items():exact_interval(c,value,recalculated[key]);rangechecks+=1
                for group,other in zip(part['dual_same_function_C0_Z_coefficient_intersections'],again['dual_same_function_C0_Z_coefficient_intersections'],strict=True):
                    for key,value in group.items():exact_interval(c,value,other[key]);rangechecks+=1
                parts+=1
            for group,other in zip(old['normalized_full_mean_C0_Z_Nminus2_coefficients'],new['normalized_full_mean_C0_Z_Nminus2_coefficients'],strict=True):
                for key,value in group.items():exact_interval(c,value,other[key]);rangechecks+=1
        left,right=p.mpf(index)/64,p.mpf(index+1)/64
        for key,rate in current.integrals.five.RATES.items():
            rr=p.mpf(rate);target=right-left if not rr else (p.exp(-rr*(1-right))-p.exp(-rr*(1-left)))/rr
            mass=saved.interval(c,row['positive_original_own_rate_masses'][key]);assert ep(mass)[0]>0
            saved.contains(mass,target,p.mpf('1e-120'));masschecks+=1
        for name in ('normalized_C0_Z_Nminus2_coefficient_unions','normalized_Nminus2_mean_Duhamel_coefficients'):
            for group,other in zip(row[name],replayed[name],strict=True):
                for key,value in group.items():exact_interval(c,value,other[key]);rangechecks+=1
    for group,other in zip(report['normalized_all_N_mean_Duhamel_C0_Z_coefficients'],encoded['normalized_all_N_mean_Duhamel_C0_Z_coefficients'],strict=True):
        for key,value in group.items():exact_interval(c,value,other[key]);rangechecks+=1
    return live,dict(passed=True,actual_N_free_true_phase_source_parts=parts,
        source_C0_Z_native_input_rows=inputchecks,native_positive_coefficient_rows=nativechecks,
        actual_dual_range_phase_weight_branch_union_and_total_checks=rangechecks,
        independent_positive_own_rate_masses=masschecks,
        all_source_A_parameters_at_most_two=True,exact_source_domain_basis_predicates_revalidated=True,
        native_same_frame_products_before_unit_division=True)


def candidate_units_and_runtime_guards(owner,live,manifest):
    c=owner.ctx;queries=0;nativechecks=0;rejected=0
    coefficients=live['coefficients']
    for text,stored in manifest['candidate_mean_bias_queries'].items():
        N=int(text);actual=owner.at_candidate(live,N=N);factor=c.mpf(1)/(c.mpf(N)**2)
        assert stored['explicit_candidate_N']==N and stored['candidate_not_global_N_admission']
        assert stored['same_all_N_source_coefficient_bounds_used']
        for j,group in enumerate(stored['normalized_C0_Z_mean_contributions']):
            for key,value in group.items():
                exact_interval(c,value,coefficients[j][key]*factor);queries+=1
        encoded=current.base.encoded(actual['record'])
        for name in ('C0','ordinary_Z'):
            for key,value in stored[name].items():exact_native(c,value,encoded[name][key]);nativechecks+=1
    new=manifest['candidate_mean_bias_queries']['160']['normalized_C0_Z_mean_contributions']
    old=owner.saved['all_normalized_C0_Z_mean_bias_contributions']
    for j,group in enumerate(manifest['additional_N160_direct_and_all_N_mean_intersection']):
        for key,value in group.items():
            expected=current.intersect(c,saved.interval(c,new[j][key]),saved.interval(c,old[j][key]))
            exact_interval(c,value,expected);queries+=1
    for call in (lambda:owner.at_candidate(live,N=True),lambda:owner.at_candidate(live,N=159),
            lambda:owner.at_candidate(live,N=160.0),lambda:owner.at_candidate(dict(live),N=160),
            lambda:owner.at_candidate({},N=160),lambda:owner.at_candidate(live)):
        try:call()
        except (ValueError,TypeError):rejected+=1
    assert rejected==6
    frame=owner.source.frame(64,0,branch='regular')
    part=owner.saved['original_source_true_phase_mean_records'][0]['conditional_source_true_phase_mean_records'][0]['complete_original_half_phase_partition'][0]
    from dataclasses import replace
    for call in (lambda:owner.phase_inputs(replace(frame),part),
            lambda:owner.phase_inputs(frame,{**part,'source_index':1}),
            lambda:owner.phase_inputs(frame,{**part,'branch':'positive'})):
        try:call()
        except ValueError:rejected+=1
    assert rejected==9
    return dict(passed=True,explicit_candidate_Nminus2_and_N160_dual_queries=queries,
        native_C0_and_ordinary_Z_unit_queries=nativechecks,
        copied_missing_wrong_N_and_source_partition_guards=rejected,
        same_runtime_owner_integral_required=True,no_global_N_admission=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('sharp_O2_phase_averaging_installed','oscillatory_spatial_integral_remainder_enclosed',
        'actual_all_route_incoming_histories_installed','functional_terminal_identity_solved','current_whole_N_selected',
        'all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
        *current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2UniformMeanBias();checks={}
    with mp.workdps(owner.ctx.dps+40):
        checks['independent_all_N_pair_envelope_calculus']=independent_uniform_pair_calculus(owner)
        print('independent_all_N_pair_envelope_calculus PASS',flush=True)
        live,checks['actual_native_source_coefficients_and_Duhamel_masses']=actual_source_coefficients_and_masses(owner,manifest)
        print('actual_native_source_coefficients_and_Duhamel_masses PASS',flush=True)
        checks['all_N_queries_units_and_runtime_guards']=candidate_units_and_runtime_guards(owner,live,manifest)
        print('all_N_queries_units_and_runtime_guards PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Uniform original O2 all integer N>=160 five mean C0/Z N^-2 coefficient envelopes; physical amplitude, source q factors, same native inputs and own-rate means retained. Mean component only; spatial oscillatory remainder, actual all-route histories, terminal/global N and full recursion/corrected reconstruction remain open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Original O2 native all-N mean bias PASS',flush=True);return report


if __name__=='__main__':run()

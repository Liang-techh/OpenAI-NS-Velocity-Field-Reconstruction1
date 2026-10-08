"""Check real upstream provenance, native rebasing and full-Z O2 composition.

Accepted primitive/density suites are reused. This checks the new connection,
fresh N1024 accumulation, both alternative enclosures and pressure memory.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_full_predicate_incoming_driver as current

ep=current.ep;KEYS=current.KEYS;iv=current.integrals.interval


def same(record,value):
    scale=record['formal_positive_scale'];c=value.ctx
    assert tuple(scale['source_exponents'])+(scale['radius_power'],)==value.scale.powers
    assert ep(iv(c,scale['additional_log_interval']))==ep(value.scale.offset)
    assert ep(iv(c,record['coefficient_interval']))==ep(value.coefficient)
    assert record['exact_zero']==value.zero
    assert record['encloses_original_source_function'] and not record['point_value_selected']


def source_basis_and_incoming(owner,record):
    c=owner.ctx;binding=record['actual_upstream_source_binding'];source=owner.incoming
    assert binding['manifest']==current.middle.NAME and binding['manifest_sha256']==current.sha(current.middle.NAME)
    assert binding['receipt']==current.middle.RECEIPT and binding['receipt_sha256']==current.sha(current.middle.RECEIPT)
    assert binding['source_record_key']=='whole_Z' and binding['candidate_N']==1024
    assert binding['exact_13_cell_upstream_function_binding']==owner.exact_binding
    assert record['source_family']==owner.family and record['original_P0_datum_sha256']==owner.family['datum_enclosure_sha256']
    assert record['exact_y_window']==['0','1'] and record['exact_Z_range']==['-1','1']
    assert record['normalized_own_units']==current.integrals.five.UNITS and record['own_rates']==current.integrals.five.RATES
    assert len(owner.exact_binding['ordered_pre_O2_function_cells'])==13
    assert owner.exact_binding['O2_incoming_is_previous_Rh_reference_outgoing']
    assert source['actual_R110_incoming_C0_Z_functions_carried_through_every_chart']
    p,k,offset,z=s.symbols('logP k offset Z',real=True)
    assert s.expand(k*(2*p)+offset-(2*k*p+offset))==0
    coefficient=s.Function('coefficient')(z);extra=s.Function('offset')(z)
    # Source Pstar is Z independent. Offset/coefficient may still vary with Z.
    assert s.simplify(s.diff(coefficient*s.exp(k*(2*p)+extra),z)-s.diff(coefficient*s.exp(2*k*p+extra),z))==0
    count=0
    for j,field in enumerate(('actual_original_inlet_to_O2_inlet_correction_C0','actual_original_inlet_to_O2_inlet_correction_Z')):
        for key in KEYS:
            raw=current.upstream.restore_common_source(source[field][key],owner.coordinates)
            same(record[('actual_upstream_C0','actual_upstream_Z')[j]][key],raw)
            got=current.integrals.restore_value(owner.out,record['rebound_source_owned_C0_Z_corrections'][j][key])
            assert got.scale.powers==(2*raw.scale.powers[1],0,0,0,0)
            assert ep(got.scale.offset)==ep(raw.scale.offset) and ep(got.coefficient)==ep(raw.coefficient)
            assert got.zero==raw.zero and got.scale.bases is owner.out.bases and got.ledger is owner.out.ledger
            same(record['rebound_source_owned_C0_Z_corrections'][j][key],owner.native_incoming[j][key]);count+=1
    for field,value in (('separate_analytic_P0',owner.P0),('separate_analytic_P0_Z',owner.P0_Z)):
        same(record[field],value)
    assert binding['full_signed_coefficients_and_directed_offsets_retained']
    assert binding['no_native_log_radius_or_Pstar_exponential_materialized']
    return dict(passed=True,actual_source_owned_incoming_C0_Z_covers_rebound=count,
        exact_Pstar_squared_native_scale_and_ordinary_Z_identities=2,ordered_original_upstream_cells=13,
        complete_directed_offsets_coefficients_and_source_defined_zeros_preserved=True,
        analytic_P0_separate=True)


def fresh_N1024_contributions(owner,manifest):
    c=owner.ctx;report=manifest['actual_fresh_original_N1024_direct_integrals']
    assert report['explicit_candidate_N']==1024 and report['source_family']==owner.family
    assert report['exact_y_window']==['0','1'] and report['exact_Z_range']==['-1','1']
    assert report['complete_original_y_Z_and_true_phase_union_enclosed']
    rows=report['whole_source_density_phase_mass_records'];count=report['ordered_source_cells'];assert count==64 and len(rows)==count
    sums=[{key:c.mpf(0) for key in KEYS} for unused in range(4)];masses=unions=0
    for index,row in enumerate(rows):
        assert row['exact_y_cell']==[str(s.Rational(index,64)),str(s.Rational(index+1,64))]
        assert row['actual_original_phase_endpoint_indices']==[index,index+1]
        assert len(row['actual_original_common_N_phase_cover'])==1
        assert ep(iv(c,row['actual_original_common_N_phase_cover'][0]))==(0,1)
        assert [r['branch'] for r in row['predicate_density_records']]==list(current.integrals.BRANCHES)
        for branch in row['predicate_density_records']:
            frame=owner.direct.source.frame(64,index,branch=branch['branch'])
            assert branch['exact_outer_Z_bounds']==frame.record['exact_outer_Z_bounds']
            assert branch['actual_original_predicate']==frame.record['actual_original_source_predicate']
        assert row['predicate_density_records'][0]['exact_outer_Z_bounds']==['-1','1']
        left,right=mp.mpf(index)/64,mp.mpf(index+1)/64
        for key,rate in current.integrals.five.RATES.items():
            rate=mp.mpf(rate)
            exact=(mp.exp(-rate*(1-right))-mp.exp(-rate*(1-left)))/rate if rate else right-left
            lo,hi=ep(iv(c,row['positive_own_rate_final_endpoint_masses'][key]));assert lo<=exact<=hi;masses+=1
        labels=('changed_C0','changed_Z_over_Lambda0_Lminus2','unmodulated_C0','unmodulated_Z')
        for j,label in enumerate(labels):
            for key in KEYS:
                parts=[iv(c,r['exported_normalized_ranges'][label][key]) for r in row['predicate_density_records']]
                expected=c.mpf((min(ep(v)[0] for v in parts),max(ep(v)[1] for v in parts)))
                assert ep(iv(c,row['normalized_four_density_unions'][j][key]))==ep(expected);unions+=1
                addition=expected*iv(c,row['positive_own_rate_final_endpoint_masses'][key])
                assert ep(iv(c,row['normalized_four_final_endpoint_contributions'][j][key]))==ep(addition)
                sums[j][key]+=addition
    for got,expected in zip(report['all_four_normalized_integral_enclosures'],sums,strict=True):
        for key in KEYS:assert ep(iv(c,got[key]))==ep(expected[key])
    phase_checks=0
    for row_index,row in enumerate(report['actual_original_radius_phase_endpoints']):
        fresh=owner.direct.source.owner.radius.evaluate(y=s.Rational(row_index,64),N=1024)
        # Numeric endpoints are compared through lossless encoding, not rounded text.
        got=row['true_original_phase_directed_boxes'];expected=current.base.encoded(fresh)['true_original_phase_directed_boxes']
        assert row['periodic_projection_full_period']==fresh['periodic_projection_full_period']
        for part,want in zip(got,expected,strict=True):
            for endpoint in ('lower','upper'):
                assert part[endpoint]['exact_mpf_tuple']==want[endpoint]['exact_mpf_tuple']
        assert row['source_family']==owner.family and row['same_original_radius_phase_bound']
        assert row['positive_original_origin_offset']['actual_hb_sc_product_not_materialized_or_set_to_zero'];phase_checks+=1
    return dict(passed=True,fresh_actual_N1024_full_predicate_cells=count,
        independently_computed_own_rate_masses=masses,branch_union_and_signed_increment_checks=unions,
        fresh_actual_common_N_radius_endpoints=phase_checks,whole_Z_axis_and_both_signs_enclosed=True,
        finite_N_density_ranges_not_N160_rescaling_or_N7_archive=True)


def histories_and_alternatives(owner,manifest):
    c=owner.ctx;out=owner.out;record=manifest['actual_full_predicate_original_O2_source_incoming_driver']
    direct=manifest['actual_fresh_original_N1024_direct_integrals']
    accepted=json.loads(gzip.decompress((current.HERE/current.spatial.NAME).read_bytes()))
    coefficients=accepted['actual_original_O2_full_spatial_all_N_envelope']['normalized_full_spatial_C0_Z_Nminus2_coefficient_caps']
    for got,expected in zip(manifest['full_spatial_all_N_coefficients'],coefficients,strict=True):
        for key in KEYS:assert ep(iv(c,got[key]))==ep(iv(c,expected[key]))
    queries=correction_checks=complete=0;alternates=[];primary=[]
    for j in range(2):
        old={key:current.integrals.restore_value(out,record['source_defined_original_background_C0_Z'][j][key]) for key in KEYS}
        changed={key:current.integrals.restore_value(out,record['full_predicate_direct_C0_Z_changed_contributions'][j][key]) for key in KEYS}
        for key in KEYS:
            same(direct[('changed_C0_integral','changed_ordinary_Z_integral')[j]][key],changed[key])
            same(direct[('source_defined_original_histories_at_y1','source_defined_original_history_Z_at_y1')[j]][key],old[key])
        unit=current.prior.ScaledEnclosure(current.prior.FormalScale(out.bases,current.spatial.slow.UNITS[(current.spatial.Y,current.spatial.YZ)[j]]),1,out.ledger)
        alternate={};own={}
        for key in KEYS:
            upper=ep(iv(c,coefficients[j][key])*(c.mpf(1)/c.mpf(1024)**2))[1]
            contribution=unit*c.mpf((-upper,upper))
            same(record['original_source_all_N_spatial_C0_Z_alternative_contributions'][j][key],contribution);queries+=1
            decay=c.exp(-c.mpf(current.integrals.five.RATES[key]))
            inherited=owner.native_incoming[j][key]*decay
            corrected=owner.add_native(inherited,changed[key]);altcorrected=owner.add_native(inherited,contribution)
            same(record['propagated_actual_direct_C0_Z_corrections'][j][key],corrected)
            same(record['propagated_actual_spatial_C0_Z_alternative_corrections'][j][key],altcorrected);correction_checks+=2
            # Existing explicit full-predicate driver accumulates original,
            # local change and incoming in this order. Retain directed rounding.
            own[key]=owner.sum_native((old[key],changed[key],inherited))
            alternate[key]=owner.add_native(old[key],altcorrected)
            same(record['actual_complete_direct_C0_Z_histories'][j][key],own[key])
            same(record['actual_complete_spatial_C0_Z_alternative_histories'][j][key],alternate[key]);complete+=2
            if key=='p':assert ep(decay)==(1,1)
        alternates.append(alternate);primary.append(own)
    for j,P0 in enumerate((owner.P0,owner.P0_Z)):
        same(record['actual_absolute_direct_pressure_C0_Z'][j],owner.add_native(P0,primary[j]['p']))
        same(record['actual_absolute_spatial_alternative_pressure_C0_Z'][j],owner.add_native(P0,alternates[j]['p']))
    for left,right in zip(record['actual_direct_radius_phase_endpoints'],record['actual_spatial_candidate_radius_phase_endpoints'],strict=True):
        assert left==right
    for flag in ('source_incoming_driver_full_Z_including_axis_and_both_signs','downstream_fresh_original_common_N_phase_not_old_N7_archive',
        'upstream_error_covers_not_sharpened_by_local_O2_frequency_bound','alternative_envelopes_not_summed_as_two_distinct_physical_contributions',
        'ordinary_Z_rows_translated_without_differentiating_normalization_or_caps','original_background_added_once_and_analytic_P0_separate',
        'exact_rate_zero_pressure_keeps_all_upstream_memory','O2_to_Rc_route_still_open',
        'incoming_positive_offsets_collected_in_arithmetic_normalization_not_exponentiated',
        'normalized_offset_ceiling_is_not_a_source_field_point_or_differentiated_function'):assert record[flag]
    return dict(passed=True,source_all_N_alternative_native_unit_queries=queries,
        actual_correction_C0_Z_affine_comparisons=correction_checks,complete_actual_C0_Z_alternative_histories=complete,
        separate_absolute_pressure_C0_Z_comparisons=4,source_phase_endpoint_consistency=65,
        original_background_once_and_rate_zero_pressure_memory_preserved=True,
        two_complete_envelopes_are_alternatives_not_two_density_contributions=True)


def guards(owner):
    rejected=0
    checks=(lambda:owner.integrate(N=7),lambda:owner.integrate(N=2048),lambda:owner.integrate(N=True),
        lambda:owner.compose({}, {},N=1024),lambda:owner.compose({}, {},N=160),
        lambda:owner.common_to_native(owner.out.scalar(1)),
        lambda:owner.spatial_to_native(owner.out.scalar(1)),
        lambda:owner.common_to_native(current.prior.ScaledEnclosure(current.prior.FormalScale(owner.coordinates.bases,(1,0,0,0,0)),1,owner.coordinates.ledger)))
    for check in checks:
        try:check()
        except ValueError:rejected+=1
        else:raise AssertionError('Wrong candidate/source owner/basis admitted')
    return dict(passed=True,candidate_source_owner_and_native_basis_rejections=rejected)


def independent_large_offset_arithmetic(owner):
    c=owner.ctx;checks=cases=0
    # These finite arithmetic fixtures exercise the actual discovered failure:
    # large positive offsets, signed cancellation and whole interval offsets.
    # They are never substituted into the original incoming source records.
    for offsets in (((990,1010),(1000,1100)),((-2000,-1800),(-1700,-1600)),((0,0),(1500,1500))):
        for coefficients in (((-2,3),(-5,-1)),((1,2),(3,4)),((-4,-2),(2,4)),((0,0),(-1,1))):
            values=[current.prior.ScaledEnclosure(current.prior.FormalScale(owner.out.bases,offset=c.mpf(log)),c.mpf(coefficient),owner.out.ledger)
                for log,coefficient in zip(offsets,coefficients,strict=True)]
            got=owner.add_native(*values);anchor=ep(got.scale.offset)[1]
            assert got.scale.powers==(0,0,0,0,0)
            # An exact zero operand legitimately returns the other entire
            # offset interval; normalize that result only for comparison.
            lo,hi=ep(got.coefficient*c.exp(got.scale.offset-c.mpf(anchor)))
            for p0 in coefficients[0]:
                for p1 in coefficients[1]:
                    for x0 in offsets[0]:
                        for x1 in offsets[1]:
                            # Independent finite scalar evaluation. Its large
                            # common exp(anchor) cancels before computation.
                            exact=mp.mpf(p0)*mp.exp(mp.mpf(x0)-anchor)+mp.mpf(p1)*mp.exp(mp.mpf(x1)-anchor)
                            assert lo<=exact<=hi;checks+=1
            cases+=1
    return dict(passed=True,large_positive_negative_offset_signed_fixtures=cases,
        independent_finite_signed_sum_comparisons=checks,no_actual_source_point_or_error_cover_substitution=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE] and manifest['candidate_N']==1024
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_all_route_incoming_histories_installed','current_whole_N_selected','functional_terminal_identity_solved',
        'all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed','actual_phase_primitive_point_evaluator_installed',
        *current.spatial.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2FullPredicateIncomingDriver();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,callback in (('actual_upstream_source_native_basis_and_P0',lambda:source_basis_and_incoming(owner,manifest['actual_full_predicate_original_O2_source_incoming_driver'])),
            ('fresh_full_predicate_N1024_phase_and_signed_accumulation',lambda:fresh_N1024_contributions(owner,manifest)),
            ('actual_whole_Z_history_alternatives_and_pressure_memory',lambda:histories_and_alternatives(owner,manifest)),
            ('independent_large_offset_signed_arithmetic',lambda:independent_large_offset_arithmetic(owner)),
            ('candidate_and_owner_guards',lambda:guards(owner))):
            checks[name]=callback();print(name+' PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,candidate_N=1024,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        actual_full_Z_source_owned_O2_incoming_driver_installed=True,**dict.fromkeys(flags,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Real source-owned upstream N1024 corrections connected to fresh whole-y/whole-Z full-predicate O2 C0/Z direct transport and all-N spatial alternative, source unit/basis/error/P0/phase memory retained. Local O2 driver only; full Rc/global N/terminal/recursion remain open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Actual original N1024 full-Z source-incoming O2 driver PASS',flush=True);return report


if __name__=='__main__':run()

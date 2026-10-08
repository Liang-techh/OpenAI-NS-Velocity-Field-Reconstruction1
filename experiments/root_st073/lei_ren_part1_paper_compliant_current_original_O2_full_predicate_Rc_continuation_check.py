"""Focused actual-source six-cell continuation, width and pressure checks."""
import gzip
import hashlib
import json
from pathlib import Path
from fractions import Fraction
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_full_predicate_Rc_continuation as current

ep=current.ep;KEYS=current.KEYS;iv=current.driver.integrals.interval


def same(record,value):
    c=value.ctx;scale=record['formal_positive_scale']
    assert tuple(scale['source_exponents'])+(scale['radius_power'],)==value.scale.powers
    assert ep(iv(c,scale['additional_log_interval']))==ep(value.scale.offset)
    assert ep(iv(c,record['coefficient_interval']))==ep(value.coefficient)
    assert record['exact_zero']==value.zero and record['encloses_original_source_function'] and not record['point_value_selected']


def overlaps(left,right):
    assert left.scale.bases is right.scale.bases and left.ledger is right.ledger
    anchor=left.ctx.mpf(max(ep(left.scale.evaluate())[1],ep(right.scale.evaluate())[1]))
    a=left.coefficient*left.bounded_exp(left.scale.evaluate()-anchor)
    b=right.coefficient*right.bounded_exp(right.scale.evaluate()-anchor)
    return max(ep(a)[0],ep(b)[0])<=min(ep(a)[1],ep(b)[1])


def native_source_inheritance(owner,saved):
    binding=saved['accepted_full_predicate_O2_correction_binding']
    assert binding==current.encode(owner.inlet_binding)
    assert binding['same_live_original_pressure_parameter_owner']
    assert binding['source_manifest_sha256']==current.sha(current.driver.NAME)
    assert binding['source_receipt_sha256']==current.sha(current.driver.RECEIPT)
    assert binding['native_common_source_constructor_sha256']==current.sha(binding['native_common_source_constructor'])
    assert binding['inherited_corrections_not_original_plus_correction_histories']
    p,C,L,o,a=s.symbols('logP logC logL offset p',real=True)
    b,l=s.symbols('b l',real=True)
    assert s.expand(a*p+b*C+l*L+o-(a*s.Rational(1,2)*(2*p)+(o+b*C+l*L)))==0
    z=s.Symbol('Z',real=True);Q=s.Function('ordinary_Z_source_coefficient')(z);logL=s.Function('source_logL')(z)
    lhs=Q*s.exp(a*p+b*C+l*logL+o);rhs=Q*s.exp(a*s.Rational(1,2)*(2*p)+o+b*C+l*logL)
    assert s.simplify(lhs-rhs)==0
    rows=0;c=owner.ctx
    for alternative in range(2):
        for j,field in enumerate(('values','Z_derivatives')):
            for key in KEYS:
                native=owner.initial_native[alternative][j][key];common=owner.initial[alternative][field][key]
                same(saved['actual_rehydrated_native_C0_Z_correction_alternatives'][alternative][j][key],native)
                same(saved['actual_rebased_common_C0_Z_correction_alternatives'][alternative][j][key],common)
                assert common.scale.bases is owner.coordinates.bases and common.ledger is owner.coordinates.ledger
                p=native.scale.powers;expected=c.mpf(ep(native.scale.offset))
                for index in (1,2):
                    if p[index]:expected+=c.mpf(ep(owner.atlas.bases[index]))*p[index]
                assert common.scale.powers==(0,p[0]/2,0,0,0)
                assert ep(common.scale.offset)==ep(expected)
                original=c.mpf(ep(native.coefficient));assert ep(common.coefficient)==ep(original);rows+=1
    return dict(passed=True,actual_accepted_correction_C0_Z_native_and_common_covers=rows,
        exact_same_source_linear_factor_and_ordinary_Z_identities=2,
        original_pressure_parameter_owner_and_Pstar_factor_constructor_bound=True,
        accepted_covers_restored_not_source_points_or_numerical_source_owners=True,
        original_upstream_and_64_cell_O2_suites_not_replayed=True)


def independent_true_widths_and_masses(owner,saved):
    c=owner.ctx;rows=saved['actual_six_downstream_source_C0_Z_cell_records'];checks=0
    expected_widths=(c.expm1(40),c.mpf(9),c.mpf(2),c.mpf(1),c.mpf(1),c.mpf(1))
    for index,(row,width) in enumerate(zip(rows,expected_widths,strict=True)):
        geometry=row['actual_true_cell_geometry']
        assert geometry['width_and_endpoints_independent_of_Z'] and geometry['width_collected_before_huge_absolute_radius_endpoint_subtraction']
        assert ep(iv(c,geometry['regular_true_log_radius_width_cover']))==ep(width)
        assert not geometry['microscopic_true_log_radius_width_components']
        positive=current.driver.upstream.restore_common_source(geometry['positive_true_log_radius_width'],owner.coordinates)
        assert positive.scale.powers==(0,0,0,0,0) and ep(positive.scale.offset)==(0,0)
        assert ep(positive.coefficient)==ep(width)
        for key,rate in current.downstream.RATES.items():
            rate=Fraction(str(rate));rr=c.mpf(rate.numerator)/rate.denominator
            factor=row['original_true_width_kernel_factors'][key]
            decay=current.driver.upstream.restore_common_source(factor['decay'],owner.coordinates)
            mass=current.driver.upstream.restore_common_source(factor['mass'],owner.coordinates)
            assert ep(decay.coefficient)==(1,1) and decay.scale.powers==(0,0,0,0,0)
            assert ep(decay.scale.offset)==ep(-rr*width)
            assert mass.scale.powers==(0,0,0,0,0) and ep(mass.scale.offset)==(0,0)
            assert ep(mass.coefficient)[0]>0
            if not rate:
                assert factor['branch']=='exact_rate0' and ep(mass.coefficient)==ep(width)
            else:
                # Independent scalar integral, using the source width rather
                # than incorrectly assuming one native-coordinate unit.
                wl,wh=ep(width)
                for w in (wl,wh):
                    expected=-mp.expm1(-mp.mpf(rate.numerator)/rate.denominator*w)/(mp.mpf(rate.numerator)/rate.denominator)
                    assert ep(mass.coefficient)[0]<=expected<=ep(mass.coefficient)[1]
            checks+=1
        if index>=4:
            endpoints=geometry['original_endpoint_specification']
            assert endpoints==[dict(original_power_offset=index-4),dict(original_power_offset=index-3)]
    return dict(passed=True,independent_original_chart_widths=6,independent_own_rate_mass_and_decay_pairs=checks,
        axial_log_radius_width_is_exp40_minus1_not_one=True,power_offsets_have_unit_log_width_not_unit_phase=True)


def actual_six_cell_replay_and_memory(owner,saved):
    live=owner.integrate()['report'];assert current.encode(live)==saved
    rows=saved['actual_six_downstream_source_C0_Z_cell_records'];alts=saved['spatial_inlet_alternative_cell_records']
    assert saved['ordered_original_six_cell_route']==[dict(label=l,chart=ch,left=a,right=b) for l,ch,a,b in current.ROUTE]
    assert len(rows)==6 and len(alts)==6
    incoming=owner.initial[0];alternative=owner.initial[1];affine=ownchecks=quiet=pressure=0
    for step,((label,chart,left,right),row,altrow) in enumerate(zip(current.ROUTE,rows,alts,strict=True),1):
        assert row['source_family']==owner.family and row['candidate_N']==1024 and ep(iv(owner.ctx,row['Z_box']))==(-1,1)
        assert row['chart']==chart and row['ordered_continuation_step']==step
        assert row['new_full_predicate_original_O2_inherited_correction_used']
        proof=row['original_chart_uniform_a_positive_certificate']
        flag='whole_actual_source_positive_not_inferred_from_saved_denominator_box' if chart.startswith('O3_') else 'source_function_positivity_not_inferred_from_saved_box'
        assert proof[flag]
        assert row['original_right_background_and_separate_P0_Z']['P0_not_merged_into_pressure_history']
        values={k:current.driver.upstream.restore_common_source(v,owner.coordinates) for k,v in row['actual_true_width_C0_contributions'].items()}
        derivatives={k:current.driver.upstream.restore_common_source(v,owner.coordinates) for k,v in row['actual_true_width_Z_contributions'].items()}
        next0={};nextZ={};alt0={};altZ={}
        for key in KEYS:
            decay=current.driver.upstream.restore_common_source(row['original_true_width_kernel_factors'][key]['decay'],owner.coordinates)
            for j,(field,local) in enumerate((('values',values),('Z_derivatives',derivatives))):
                same(row[('actual_inherited_correction_C0','actual_inherited_correction_Z')[j]][key],incoming[field][key])
                got=decay*incoming[field][key]+local[key]
                altgot=decay*alternative[field][key]+local[key]
                same(row[('actual_right_correction_C0','actual_right_correction_Z')[j]][key],got)
                same(altrow['inherited_C0_Z_alternative'][j][key],alternative[field][key])
                same(altrow['corrected_C0_Z_alternative'][j][key],altgot);affine+=2
                background=current.driver.upstream.restore_common_source(row['original_right_background_and_separate_P0_Z'][('original_normalized_history_C0_enclosures','original_normalized_history_Z_enclosures')[j]][key],owner.coordinates)
                same(row[('actual_right_own_history_C0','actual_right_own_history_Z')[j]][key],background+got)
                same(altrow['own_C0_Z_alternative'][j][key],background+altgot);ownchecks+=2
                (next0,nextZ)[j][key]=got;(alt0,altZ)[j][key]=altgot
                if chart=='O3_power':
                    assert local[key].zero
                    assert row['actual_signed_density_C0_covers'][key]['exact_zero'] and row['actual_signed_density_Z_covers'][key]['exact_zero'];quiet+=1
                    if key=='p':
                        assert not incoming[field][key].zero and not alternative[field][key].zero
                        same(row[('actual_right_correction_C0','actual_right_correction_Z')[j]][key],incoming[field][key]);pressure+=1
        assert altrow['alternative_not_added_to_direct_history'] and altrow['identical_actual_source_cell_operator_used']
        if chart=='O3_power':
            primitive=row['original_periodic_parameter_C1_cover']
            assert primitive['original_full_box_q_and_q_Z_exact_zero_implies_primitive_C0_Z_exact_zero']
            assert row['original_whole_power_q_q_Z_density_C0_Z_exact_zero']
            assert row['original_actual_power_offsets']==[left,right]
        incoming=dict(values=next0,Z_derivatives=nextZ);alternative=dict(values=alt0,Z_derivatives=altZ)
        # The serialized step operators already bind original geometry; the
        # complete live composite below is compared to both serial alternatives.
    for j,(serial,initial) in enumerate(zip((incoming,alternative),owner.initial,strict=True)):
        for k in KEYS:
            for d,field in enumerate(('values','Z_derivatives')):
                same(saved[('actual_Rc_direct_correction_C0_Z','actual_Rc_spatial_alternative_correction_C0_Z')[j]][d][k],serial[field][k])
                value=current.driver.upstream.restore_common_source(saved['direct_and_spatial_alternative_composite_corrections'][j][d][k],owner.coordinates)
                assert overlaps(value,serial[field][k])
    endpoint=saved['original_Rc_endpoint'];Tw=owner.shell.transfer.geometry.binder.fixed['Tw']
    assert endpoint['original_power_offset']==2 and endpoint['original_power_phase_expression']=='2/Tw'
    assert ep(iv(owner.ctx,endpoint['actual_direct_endpoint_phase_cover']))==ep(2/Tw)
    assert endpoint['original_Rc_not_power_phase1']
    background=saved['original_Rc_background_and_separate_P0_Z']
    for j in range(2):
        P0=current.driver.upstream.restore_common_source(background[('original_separate_P0_over_Pstar_squared','original_separate_P0_Z_over_Pstar_squared')[j]],owner.coordinates)
        for field,pressurefield in (('actual_Rc_direct_own_history_C0_Z','absolute_Rc_direct_pressure_C0_Z'),
            ('actual_Rc_spatial_alternative_own_history_C0_Z','absolute_Rc_spatial_alternative_pressure_C0_Z')):
            value=current.driver.upstream.restore_common_source(saved[field][j]['p'],owner.coordinates)
            same(saved[pressurefield][j],P0+value)
    return dict(passed=True,fresh_original_downstream_whole_Z_source_cells=6,
        direct_and_alternative_C0_Z_affine_inheritance_comparisons=affine,
        original_background_once_C0_Z_comparisons=ownchecks,quiet_power_C0_Z_density_and_increment_rows=quiet,
        nonzero_inherited_quiet_pressure_memory_checks=pressure,complete_serial_composite_C0_Z_overlaps=20,
        separate_absolute_pressure_C0_Z_comparisons=4,original_Rc_power_offset_two_verified=True,
        accepted_upstream_and_original_O2_suites_not_rerun=True)


def guards(owner):
    rejected=0
    values=(lambda:owner.integrate(N=7),lambda:owner.integrate(N=2048),lambda:owner.integrate(N=True),
        lambda:owner.native_to_common(owner.coordinates.scalar(1)),
        lambda:owner.native_to_common(current.prior.ScaledEnclosure(current.prior.FormalScale(owner.atlas.bases,(0,0,1,1,0)),1,owner.atlas.ledger)),
        lambda:owner.native_to_common(current.prior.ScaledEnclosure(current.prior.FormalScale(owner.atlas.bases,(.5,0,0,0,0)),1,owner.atlas.ledger)))
    for callback in values:
        try:callback()
        except ValueError:rejected+=1
        else:raise AssertionError('Wrong candidate/native owner/uncollected source factors admitted')
    return dict(passed=True,candidate_owner_and_source_factor_guards=rejected)


@current.native.inlet.source_precision
def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE] and manifest['candidate_N']==1024
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_all_route_incoming_histories_installed','global_inlet_to_Rc_histories_admitted','functional_terminal_identity_solved',
        'current_whole_N_selected','all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',*current.downstream.packets.OPEN)
    assert all(manifest[flag] is False for flag in flags)
    bridge,_=current.native.inlet.native_bridge_owner();checks={}
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.OriginalO2FullPredicateRcContinuation(bridge);saved=manifest['actual_full_predicate_O2_correction_to_Rc_continuation']
        with mp.workdps(owner.ctx.dps+40):
            for label,callback in (('actual_accepted_source_cover_native_common_inheritance',lambda:native_source_inheritance(owner,saved)),
                ('independent_original_six_cell_true_width_mass_and_decay',lambda:independent_true_widths_and_masses(owner,saved)),
                ('actual_six_cell_source_replay_and_pressure_memory',lambda:actual_six_cell_replay_and_memory(owner,saved)),
                ('candidate_source_factor_guards',lambda:guards(owner))):
                checks[label]=callback();print(label+' PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,candidate_N=1024,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        actual_source_owned_full_predicate_to_Rc_C0_Z_covers_installed=True,**dict.fromkeys(flags,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Actual inherited original full-predicate O2 N1024 correction/whole-Z derivatives reach source-owned Rc through all six downstream cells. Source factors, actual geometry/masses, quiet power memory and separate P0 retained. Conservative function covers; no terminal/global N/source oracle/cone/recursion/full NS admission.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.encode(report),indent=2).encode()+b'\n')
    print('Actual original full-predicate O2 correction to Rc PASS',flush=True);return report


if __name__=='__main__':run()

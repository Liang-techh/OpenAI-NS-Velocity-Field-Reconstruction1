"""Focused saved-source density, true-mass and actual-period C1 integral audit.

Reconstructs the original density from its actual inverse/first-jet records
and independently replays directed integral transport. No native source owner
or preceding O2 quadrature is rerun.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_tail_phase_integrals as source
import lei_ren_part1_paper_compliant_current_weighted_O2_to_Rc_tail_check as accepted

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv
common=source.common;rc=source.rc;KEYS=source.KEYS
same,native_restore=accepted.same,accepted.native_restore


def cell_check(c,coordinates,row,Z):
    assert row['candidate_N']==source.N and row['source_family']==coordinates.family
    assert ep(c.mpf(row['exact_Z_range']))==Z
    assert row['actual_N_log_radius_phase_not_independent_phase_samples']
    assert row['whole_period_primitive_caps_not_used_as_A_B']
    assert row['original_ordinary_Z_density_derivatives_not_cap_derivatives']
    native=row['actual_original_phase_inverse_and_density_source'];prov=native['source_provenance']
    assert native['chart']==row['chart']==prov['chart'] and native['candidate_N']==source.N
    assert prov['source_family']==coordinates.family and ep(iv(c,prov['Z_box']))==Z
    assert prov['one_original_seed_graph_asserted'] and not prov['cache_cover']
    assert prov['numerical_caps_not_used_as_defining_function_values']
    phase=native['actual_original_radius_phase'];g=row['actual_true_geometry']
    assert phase['phase_independent_of_Z'] and phase['candidate_N']==source.N
    assert phase['source_radius_caps_not_consumed'] and phase['original_spatial_phase_bound_at_candidate_N']
    assert ep(iv(c,phase['coordinate_box']))==ep(iv(c,g['native_coordinate_box']))==ep(iv(c,prov['coordinate_box']))
    pieces=native['actual_spatial_signed_density_Z_cells']
    assert len(pieces)==len(phase['periodic_projection']['boxes'])
    if row['status']!='enclosed':
        assert not native['five_signed_density_Z_functions_installed']
        assert any(piece['status']!='enclosed' for piece in pieces)
        assert not any(k in row for k in ('original_local_C0_integral','original_local_Z_integral'))
        return None,0
    assert native['five_signed_density_Z_functions_installed']
    bases=tuple(iv(c,v) for v in row['original_native_source_log_bases'])
    assert ep(bases[1])==ep(coordinates.logP_squared)
    assert ep(bases[4])==ep(iv(c,prov['logR_cover']))
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    restore=lambda record:native_restore(record,bases,ledger)
    kernels=[];jets=[]
    for piece,phi in zip(pieces,phase['periodic_projection']['boxes'],strict=True):
        assert piece['status']=='enclosed' and piece['original_radius_phase_Z_derivative_exactly_zero']
        assert piece['same_original_factor_basis_and_ledger'] and piece['source_native_width_or_Pstar_conversion_not_reapplied']
        first=piece['original_phase_first_jet_source']
        assert not first['phase_is_independent_input_on_this_query'] and first['actual_spatial_fast_N_chain_installed']
        assert ep(iv(c,first['derived_spatial_fractional_phase_box']))==ep(iv(c,phi))
        primitives={k:restore(v) for k,v in first['original_A_B_first_derivative_enclosures'].items()}
        velocity={k:restore(v) for k,v in piece['original_and_candidate_normalized_velocity_Z_enclosures'].items()}
        # Accepted density differentiation is replayed on genuine stored
        # inverse outputs, never on hull selectors or period primitive caps.
        got=rc.density.density_Z_kernels(velocity['original_E'],velocity['original_E_Z'],
            velocity['original_V'],velocity['original_V_Z'],primitives,source.N)
        for k,v in got['velocities'].items():same(c,piece['original_and_candidate_normalized_velocity_Z_enclosures'][k],v)
        for k in KEYS:
            same(c,piece['five_original_signed_density_C0_enclosures'][k],got['kernels'][k])
            same(c,piece['five_original_signed_density_Z_derivative_enclosures'][k],got['Z_derivatives'][k])
        chain=first['original_total_spatial_first_derivative_enclosures']
        same(c,chain['A_Z_total'],primitives['A_Z']);same(c,chain['B_Z_total_over_Pstar'],primitives['B_Z_over_Pstar'])
        same(c,chain['A_y_total'],primitives['A_y']+source.N*primitives['A_phi'])
        same(c,chain['B_y_total_over_Pstar'],primitives['B_y_over_Pstar']+source.N*primitives['B_phi_over_Pstar'])
        if first['geometry']=='flat':
            assert all(v.zero for v in primitives.values()) and all(v.zero for v in (*got['kernels'].values(),*got['Z_derivatives'].values()))
        kernels.append(got['kernels']);jets.append(got['Z_derivatives'])
    geometry=dict(width=common.restore_common_source(g['positive_true_log_radius_width'],coordinates),
        regular=iv(c,g['regular_true_log_radius_width_cover']),scalar_cover=iv(c,g['scalar_width_cover_used_only_for_directed_kernel_bounds']))
    assert g['width_and_endpoints_independent_of_Z'] and row['true_log_radius_Jacobian_integrated_once']
    values={};derivatives={};factors={}
    for k,rate in rc.RATES.items():
        f=rc.transfer.true_width_kernel(coordinates,geometry,rate);saved=row['true_positive_kernel_factors'][k]
        assert saved['branch']==f['branch'];same(c,saved['mass'],f['mass']);same(c,saved['decay'],f['decay'])
        rho=rc.transfer.local.same_source_union([v[k] for v in kernels])
        rhoZ=rc.transfer.local.same_source_union([v[k] for v in jets])
        same(c,row['full_phase_union_signed_C0_density'][k],rho);same(c,row['full_phase_union_signed_Z_density'][k],rhoZ)
        values[k]=coordinates.rebase(rho,coordinates.family)*f['mass']
        derivatives[k]=coordinates.rebase(rhoZ,coordinates.family)*f['mass']
        same(c,row['original_local_C0_integral'][k],values[k]);same(c,row['original_local_Z_integral'][k],derivatives[k])
        factors[k]=f
    assert row['pressure_P0_not_added_to_density_or_integral']
    same(c,row['true_positive_kernel_factors']['p']['mass'],geometry['width'])
    return dict(geometry=geometry,values=values,Z_derivatives=derivatives,factors=factors),len(pieces)


def integral_check(c,coordinates,record,Z):
    assert record['status']=='enclosed' and record['chart']=='O2_buffer' and record['candidate_N']==source.N
    assert ep(c.mpf(record['exact_Z_range']))==Z
    lo,hi=map(Fraction,record['exact_native_endpoints']);count=record['ordered_source_cells']
    assert lo==5 and hi-lo==Fraction(1,source.N) and count in (8,32)
    assert len(record['queried_original_cells'])==count
    assert record['zero_operator_additive_identity_not_zero_actual_incoming']
    assert record['numerical_local_integral_enclosure_not_exact_field_point_or_whole_chart_closure']
    coefficients={k:coordinates.scalar(1) for k in KEYS}
    values={k:coordinates.scalar(0) for k in KEYS};jets={k:coordinates.scalar(0) for k in KEYS}
    phase_count=0
    for index,row in enumerate(record['queried_original_cells']):
        left=lo+(hi-lo)*index/count;right=lo+(hi-lo)*(index+1)/count
        assert list(map(Fraction,row['exact_native_endpoints']))==[left,right]
        got,n=cell_check(c,coordinates,row,Z);assert got is not None;phase_count+=n
        width=c.mpf(1)/(source.N*count)
        assert ep(got['geometry']['regular'])==ep(width)
        assert ep(got['geometry']['scalar_cover'])==ep(width)
        for k in KEYS:
            decay=got['factors'][k]['decay']
            coefficients[k]=decay*coefficients[k]
            values[k]=decay*values[k]+got['values'][k]
            jets[k]=decay*jets[k]+got['Z_derivatives'][k]
    operator=record['original_C1_integral_operator']
    assert operator['steps']==count and operator['original_rates']=={k:str(r) for k,r in rc.RATES.items()}
    assert operator['incoming_argument_not_assumed_or_reset'] and operator['quiet_cells_preserve_original_rate0_pressure_memory']
    for k in KEYS:
        same(c,operator['incoming_C0_Z_decay_coefficients'][k],coefficients[k])
        same(c,operator['cumulative_signed_increment_C0_enclosures'][k],values[k])
        same(c,operator['cumulative_signed_increment_Z_enclosures'][k],jets[k])
        same(c,record['five_original_C0_integrals'][k],values[k]);same(c,record['five_original_Z_integrals'][k],jets[k])
    assert ep(coefficients['p'].coefficient)==(1,1) and coefficients['p'].scale.powers==(0,0,0,0,0)
    return dict(values=values,Z_derivatives=jets),phase_count


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert not manifest['numerical_complete_tail_integrals_or_Rc_targets_or_closure_admitted']
    assert not manifest['full_Rc_closure_or_global_N_or_stress_or_recursion_claimed']
    assert all(manifest[k] is False for k in common.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    c=MPIntervalContext();c.dps=240;phase_count=0;queries=0;refinements=[]
    with mp.workdps(300):
        for archive in manifest['actual_original_tail_phase_density_source_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            payload=json.loads(raw);assert payload['source_family']==manifest['source_family'] and payload['candidate_N']==source.N
            assert not payload['full_axial_transition_or_complete_tail_numerical_integrals_admitted']
            Z=ep(c.mpf(payload['exact_Z_range']));rows=payload['full_original_tail_domain_phase_density_queries']
            assert tuple(r['original_tail_label'] for r in rows)==source.tail.LABELS
            first=next(r for r in rows if r['status']=='enclosed')
            coordinates=rc.history.CommonSourceCoordinates(c,iv(c,first['original_native_source_log_bases'][1]),manifest['source_family'])
            for row in rows:
                _,n=cell_check(c,coordinates,row,Z);phase_count+=n;queries+=1
            integrals=payload['actual_buffer_period_integral_refinements']
            assert [r['ordered_source_cells'] for r in integrals]==[8,32]
            results=[]
            for record in integrals:
                got,n=integral_check(c,coordinates,record,Z);phase_count+=n;queries+=record['ordered_source_cells'];results.append(got)
            comparisons={}
            for kind in ('values','Z_derivatives'):
                comparisons[kind]={}
                for k in KEYS:
                    a,b=(r[kind][k] for r in results)
                    assert a.scale.powers==b.scale.powers
                    anchor=max(ep(a.scale.offset)[1],ep(b.scale.offset)[1])
                    xa=a.coefficient*a.bounded_exp(a.scale.offset-anchor)
                    xb=b.coefficient*b.bounded_exp(b.scale.offset-anchor)
                    la,ua=ep(xa);lb,ub=ep(xb)
                    assert max(la,lb)<=min(ua,ub),'Disjoint true period enclosures: '+kind+' '+k
                    comparisons[kind][k]=dict(common_formal_offset=anchor,common_source_powers=list(a.scale.powers),
                        eight_cell_width=ua-la,thirty_two_cell_width=ub-lb,
                        width_ratio_32_over_8=None if ua==la else (ub-lb)/(ua-la),overlap=True)
            refinements.append(dict(exact_Z_range=payload['exact_Z_range'],actual_buffer_period_refinement_comparisons=comparisons))
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        actual_generic_source_cells_checked=queries,actual_enclosed_phase_inverse_density_C0_Z_cells_checked=phase_count,
        actual_period_integral_refinements_checked=4,actual_period_integral_C0_Z_rows_checked=40,
        true_period_length_and_single_Jacobian_checked=True,original_first_Z_jet_density_products_checked=True,
        independent_affine_integral_mass_and_downstream_decay_replay=True,refinement_comparisons=refinements,
        complete_tail_integrals_or_Rc_targets_or_physical_closure_admitted=False,**dict.fromkeys(common.current.FLAGS,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual generic tail phase/density and period C1 PASS:',queries,'source cells;',phase_count,'enclosed phase cells',flush=True)
    return result


if __name__=='__main__':run()

"""Independent scalar Z derivatives, source carriers and full-window checks."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_density_Z_integrals as current
import lei_ren_part1_paper_compliant_current_original_O2_positive_logq_cells_check as fixture
import lei_ren_part1_paper_compliant_current_generic_loop_point_Z as original

base=current.base;ep=current.ep;prev=current.current;KEYS=tuple(prev.five.RATES)
FLAGS=('numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed',
    'current_whole_N_selected',*base.point.source.inertial.profiles.loop.OPEN)


def iv(c,x):return prev.interval(c,x)


def finite_checks(c,graph):
    scales=original.loop.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max=0,p1_abs_max=8,p2_abs_max=5,dps=110)
    E,V=sy.symbols('E V');history=current.slow.base.point.source.inertial.profiles.loop
    expressions=prev.density.recovery.history_densities(E,V)
    jac={k:(sy.lambdify((E,V),sy.diff(x,E),'mpmath'),sy.lambdify((E,V),sy.diff(x,V),'mpmath')) for k,x in expressions.items()}
    derivatives=densities=coordinates=0;branches=set()
    with mp.workdps(c.dps+40):
        for p2 in ('0','.02','.2','-.2','2','-2'):
            reference=original.GenericLoopPointZ(scales,a='.8',b=0,p1=5,p2=p2,E='1.3',
                a_Z=0,b_Z=0,p2_Z='.11',E_Z='.12')
            part=fixture.finite_fixture(c,reference.base)
            for key,row in part['roots'].items():row[(0,1)]=part['kernel'].scalar(dict(p2='.11',E='.12',V='-.07').get(key,0))
            kernel=part['kernel']
            for phase in ('.137','.663'):
                expected=reference.evaluate(phase)
                C0=kernel.evaluate(phase,bits=80)['selected_inverse'];primitive=kernel.primitives(C0['coordinate_interval'],C0['chart'])
                for chart in (('psi','E') if kernel.geometry=='signed_Mobius' else ('psi',)):
                    inverse=kernel.inverse_bracket(phase,chart,80)
                    jets,proof=current.slow_values(part,inverse['coordinate_interval'],chart);branches.add(proof['branch'])
                    for key,target in (('psi_Z',expected.angle_Z),('A_Z_slow',expected.A_Z_slow),('B_Z_slow',expected.B_Z_slow)):
                        value=prev.bounded(jets[key])+c.mpf(('-1e-80','1e-80'))
                        assert fixture.contains(value,c.mpf(mp.nstr(target,115))),(p2,phase,chart,key)
                        derivatives+=1
                    values=prev.CellDensityGraph(graph,part,primitive,jets,7).outputs()['density_Z']
                    p=reference.ctx;e=p.mpf('1.3');v=p.mpf('.27');ez=p.mpf('.12');vz=p.mpf('-.07')
                    en=e*p.exp(expected.A/7);vn=v+expected.B/7
                    enz=p.exp(expected.A/7)*(ez+e*expected.A_Z_slow/7);vnz=vz+expected.B_Z_slow/7
                    for key,(fE,fV) in jac.items():
                        target=fE(en,vn)*enz+fV(en,vn)*vnz-fE(e,v)*ez-fV(e,v)*vz
                        assert fixture.contains(prev.bounded(values[key])+c.mpf(('-1e-80','1e-80')),c.mpf(mp.nstr(target,115))),(p2,phase,chart,key)
                        densities+=1
                    coordinates+=1
            for phase in ('0','.5','1'):
                jets,proof=current.slow_values(part,c.mpf(phase),'psi');assert all(v.zero for v in jets.values())
    assert {'conditioned_signed_psi','conditioned_signed_E','regular_small_r_Fourier','exact_midplane_nonzero_p2_Z'}<=branches
    return dict(passed=True,independent_original_quadrature_and_implicit_Z_components=derivatives,
        independent_original_history_Jacobian_density_Z_comparisons=densities,inverse_coordinate_cases=coordinates,
        derivative_branches=sorted(branches),reference_roundoff_allowed='1e-80',
        native_source_parameters_not_replaced_by_finite_fixtures=True,
        ordinary_source_derivatives_not_finite_differences=True)


def carrier_contract(c,record,source):
    original=source['original_source_function_coefficient_record']
    assert record['original_term_rows']==original['full_original_factored_coefficient_function_ranges']['p2'][1]
    assert record['exact_common_carrier_source_exponents']==[1,0,-1,0] and record['exact_common_carrier_radius_power']==1
    assert record['common_radial_factor_canceled_before_ratio'] and record['every_original_coefficient_and_positive_late_pressure_error_retained']
    got=record['correlated_same_original_source_root'];scale=got['formal_positive_scale']
    assert scale['source_exponents']==[1,0,-1,0] and scale['radius_power']==1
    assert ep(iv(c,scale['additional_log_interval']))==(0,0)
    assert got['encloses_original_source_function'] and not got['point_value_selected']
    assert record['original_wide_source_root']==original['native_source_root_enclosures']['p2']['y0_Z1']


def native_checks(c,manifest,owner,parent_manifest):
    totals={key:[] for key in KEYS};widths={key:[] for key in KEYS};cells=inverses=mass_checks=0
    with mp.workdps(c.dps+40):
        for result,old in zip(manifest['actual_original_five_density_Z_integral_refinements'],
                parent_manifest['actual_original_nonzero_Z_own_integral_refinements'],strict=True):
            count=result['ordered_source_cells'];rows=result['whole_original_density_and_density_Z_source_records']
            assert count==old['ordered_source_cells'] and len(rows)==count
            assert result['original_Z_exact']=='37/100' and result['exact_y_window']==['0','1'] and result['explicit_candidate_N']==7
            assert result['source_family']==manifest['source_family']==owner.family==old['source_family']
            assert result['own_rates']==prev.five.RATES and result['normalized_own_units']==prev.five.UNITS
            assert result['original_P0_datum_sha256']==owner.family['datum_enclosure_sha256'] and result['actual_incoming_histories_not_set_to_zero']
            assert result['complete_original_y_window_all_source_and_phase_pieces_integrated']
            assert not result['Z_interval_terminal_identities_or_global_controls_installed'] and not result['current_whole_N_selected']
            total0={key:c.mpf(0) for key in KEYS};total1={key:c.mpf(0) for key in KEYS}
            phase=owner.parent.parent.levels[count][1];p=mp.mp.clone();p.dps=c.dps+40
            for i,(row,before,phase_cell) in enumerate(zip(rows,old['whole_source_phase_density_records'],phase['whole_source_cells'],strict=True)):
                assert row['source']==before['source'] and row['source']['ordered_source_cell_index']==i
                assert row['five_signed_C0_density_hulls']==before['five_signed_normalized_density_hulls']
                assert row['five_C0_contributions']==before['five_signed_own_rate_contributions']
                assert row['positive_own_rate_final_endpoint_masses']==before['positive_own_rate_final_endpoint_masses']
                assert row['fixed_y_window_and_Z_independent_positive_masses']
                assert row['derivatives_integrated_from_original_graph_not_from_interval_endpoint_differences']
                pieces=row['actual_source_phase_held_Z_records'];phases=phase_cell['true_common_N_phase_boxes']
                assert len(pieces)==len(before['actual_native_inverse_phase_pieces'])
                for j,proof in enumerate(pieces):
                    assert proof['actual_true_phase_cover']==phases[j%len(phases)] and proof['source_piece_index']==j//len(phases)
                    assert proof['original_inverse']==before['actual_native_inverse_phase_pieces'][j]['selected_original_inverse']
                    derivative=proof['original_derivative_inverse'];target=iv(c,proof['actual_true_phase_cover'])
                    assert fixture.contains(iv(c,derivative['phase_image']),target)
                    assert derivative['chart'] in ('E','psi') and proof['density_Z_is_genuine_ordinary_source_derivative']
                    carrier_contract(c,proof['correlated_original_p2_Z_carrier'],row['source'])
                    contract=proof['original_slow_Z_contract']
                    assert contract['local_basis_order']==prev.LOCAL_ORDER and contract['q_Z_and_dstar_Z_exact_zero_source_identity']
                    assert contract['original_phase_held_Z_not_finite_difference'] and contract['source_Z_derivatives_not_interval_selector_derivatives']
                    assert all(v['encloses_original_source_function'] and not v['point_value_selected'] for v in proof['native_original_primitive_Z_enclosures'].values())
                    inverses+=1
                for key in KEYS:
                    h=iv(c,row['five_genuine_signed_ordinary_Z_density_hulls'][key]);m=iv(c,row['positive_own_rate_final_endpoint_masses'][key]);add=iv(c,row['five_ordinary_Z_contributions'][key])
                    assert ep(h*m)==ep(add) and all(mp.isfinite(v) for v in ep(h))
                    left,right=p.mpf(i)/count,p.mpf(i+1)/count;rate=p.mpf(prev.five.RATES[key])
                    exact=right-left if rate==0 else (p.exp(-rate*(1-right))-p.exp(-rate*(1-left)))/rate
                    assert 0<ep(m)[0]<=exact<=ep(m)[1]
                    total0[key]+=iv(c,row['five_C0_contributions'][key]);total1[key]+=add;mass_checks+=1
                cells+=1
            for key in KEYS:
                assert ep(total0[key])==ep(iv(c,result['five_original_C0_integral_contributions'][key]))==ep(iv(c,old['five_signed_original_nonzero_Z_own_rate_integral_contributions'][key]))
                got=iv(c,result['five_genuine_ordinary_Z_integral_contributions'][key]);assert ep(total1[key])==ep(got)
                totals[key].append(got);widths[key].append(ep(got)[1]-ep(got)[0])
        for key in KEYS:
            assert max(ep(v)[0] for v in totals[key])<=min(ep(v)[1] for v in totals[key])
        probes=manifest['native_whole_cell_derivative_probes'];formal_midplane=False
        for probe in probes:
            source=probe['source'];count,index=source['ordered_source_cell_level'],source['ordered_source_cell_index'];zl,zh=source['exact_Z_range']
            live=owner.query(count,index,Z_lower=zl,Z_upper=zh,phase=probe['fixed_phase'])
            assert base.encoded(live)==probe
            if not probe['source_pieces']:
                assert probe['mixed_signed_source_not_skipped'] and zl=='-1/100' and zh=='1/100';continue
            for piece in probe['source_pieces']:
                carrier_contract(c,piece['correlated_original_p2_Z_carrier'],source)
                assert piece['actual_original_density_Z_graph_executed']
                if zl==zh=='0':
                    assert any(x is None for x in piece['optional_finite_five_density_Z_ranges'].values())
                    assert piece['unmaterializable_original_Z_derivatives_retained_as_formal_sources'];formal_midplane=True
                    assert not piece['original_primitive_Z_enclosures']['A_Z_slow']['exact_zero']
        assert formal_midplane
    return dict(passed=True,whole_source_density_and_genuine_Z_integral_cells=cells,
        genuine_phase_held_Z_source_piece_records=inverses,independent_closed_form_positive_mass_checks=mass_checks,
        accepted_original_C0_integrals_exactly_preserved=True,
        all_five_derivative_refinements_have_common_intersection=True,
        derivative_refinement_interval_widths=widths,
        finest_genuine_five_ordinary_Z_contribution_intervals={k:v[-1] for k,v in totals.items()},
        original_extreme_midplane_derivatives_preserved_as_formal_sources=True,
        crossing_signed_source_remains_explicitly_unresolved=True,
        fixed_Z_derivative_enclosure_not_continuous_Z_terminal_matching=True)


def run():
    begin=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE] and manifest['fixed_nonzero_Z_C0_and_genuine_density_Z_full_window_enclosed']
    assert manifest['derivative_is_actual_phase_held_Z_not_finite_difference'] and all(manifest[k] is False for k in FLAGS)
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    assert current.identity_contract()['passed'] and manifest['exact_original_source_derivative_identities']['passed']
    parent_receipt=json.loads((current.HERE/prev.RECEIPT).read_bytes());saved=parent_receipt['compressed_producer_report']
    previous_raw=gzip.decompress((current.HERE/saved['filename']).read_bytes());assert hashlib.sha256(previous_raw).hexdigest()==saved['lossless_original_json_sha256']
    parent=json.loads(previous_raw);owner=current.OriginalO2DensityZIntegrals();c=owner.c
    scalar=finite_checks(c,owner.parent.owner.scales.graph)
    print('Original scalar quadrature and true density-Z Jacobians PASS',flush=True)
    native=native_checks(c,manifest,owner,parent)
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        independent_original_scalar_Z_and_density_derivative_checks=scalar,
        actual_native_genuine_density_Z_integral_contracts=native,
        compressed_producer_report=dict(filename=current.NAME,compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw),lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()),
        fixed_nonzero_Z_genuine_density_Z_integrals_enclosed=True,**dict.fromkeys(FLAGS,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-begin,
        scope='Complete original O2 C0 and genuine ordinary-Z five integral contribution enclosures at fixed Z37/100, candidate N7. Common p2_Z source carrier, original derivative E coordinate and fixed-phi identities checked; exact old C0 preserved. No continuous-Z functional matching, global controls/N, matched stress, recursion or full NS.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Genuine original O2 five density-Z integrals and C0 preservation PASS',flush=True);return result


if __name__=='__main__':run()

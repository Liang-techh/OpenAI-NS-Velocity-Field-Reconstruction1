"""Focused original-source, finite scalar and complete-window integral checks.

Directed native ranges certify the fixed-Z integral. Finite scalar fixtures
independently check the unchanged phase and equivalent B expression; they
are not replacements for the original extraordinarily small parameters.
No ancestor producer or defining source quadrature is reexecuted.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_positive_logq_cells as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original

base=current.base;ep=current.ep;KEYS=tuple(current.five.RATES)
FLAGS=('actual_changed_five_moment_integral_evaluated',
    'numerical_original_source_point_or_integral_oracle_installed',
    'actual_five_controls_installed','current_whole_N_selected',
    *base.point.source.inertial.profiles.loop.OPEN)


def interval(c,value):return current.interval(c,value)


def contains(a,b):
    lo,hi=ep(a);x,y=ep(b);return lo<=x and y<=hi


def same(a,b):return ep(a)==ep(b)


def fresh_ledger():
    return dict(directed_small_exponential_tails=0,
        positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)


def finite_fixture(c,reference,*,p2_cover=None,logu_cover=None):
    """Local source basis from the same original scalar defining parameters."""
    text=lambda x:mp.nstr(x,135)
    q=c.mpf(text(reference.q));p2=c.mpf(text(reference.p2))
    dstar=c.mpf(text(reference.scales.d_star))
    if p2_cover is not None:p2=c.mpf(p2_cover)
    pl,ph=ep(p2);sign=1 if pl>0 else -1 if ph<0 else 0
    if logu_cover is None:
        logu=c.mpf(0) if not sign else c.ln(abs(p2)*q/dstar)
    else:logu=c.mpf(logu_cover)
    bases=(c.mpf(0),c.mpf(0),c.mpf(0),logu,c.mpf(0));ledger=fresh_ledger()
    scalar=lambda v:base.prior.ScaledEnclosure(base.prior.FormalScale(bases),v,ledger)
    roots={key:{(0,0):scalar(value)} for key,value in dict(
        a=c.mpf(text(reference.a)),b=0,t0=0,E=c.mpf(text(reference.Utheta)),V='.27',p2=p2).items()}
    source_q=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,offset=c.ln(q)),1,ledger)
    source_u=scalar(0) if not sign else base.prior.ScaledEnclosure(
        base.prior.FormalScale(bases,(0,0,0,1,0)),sign,ledger)
    query=dict(roots=roots,q=source_q,original_u_source=source_u,ledger=ledger)
    kernel=current.PositiveLogQPhase(query,c.ln(dstar));query['kernel']=kernel
    return query


def original_scalar_at_coordinate(reference,x,chart):
    p=reference.ctx;fraction=p.mpf(mp.nstr(x,140))
    if chart=='psi' or fraction in (0,p.mpf('.5'),1):return 2*p.pi*fraction
    rho=reference.one_minus_abs_r
    plus,minus=(rho,2-rho) if reference.r>0 else (2-rho,rho)
    angle=2*p.atan2(plus*p.sin(p.pi*fraction),minus*p.cos(p.pi*fraction))
    return angle+2*p.pi if angle<0 else angle


def scalar_checks(c,graph):
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',
        boundary_kappa_excess_min='.019',t0_abs_max=0,p1_abs_max=8,p2_abs_max=5,dps=130)
    cases=0;phase_checks=0;primitive_checks=0;density_checks=0;branches=set()
    tolerance=c.mpf(('-1e-90','1e-90'))
    with mp.workdps(c.dps+40):
        for p2 in ('0','.02','-.02','.05','-.05','.2','-.2','2','-2'):
            reference=original.GenericShearLoop(scales,a='.8',b=0,p1=5,p2=p2,Utheta='1.3')
            query=finite_fixture(c,reference);kernel=query['kernel'];branches.add(kernel.geometry)
            for phase in ('0','.137','.337','.5','.663','.863','1'):
                got=kernel.evaluate(phase,bits=90);assert got['status']=='enclosed'
                selected=got['selected_inverse'];box=selected['coordinate_interval'];lo,hi=ep(box)
                left=reference.phase_at_angle(original_scalar_at_coordinate(reference,lo,selected['chart']))
                right=reference.phase_at_angle(original_scalar_at_coordinate(reference,hi,selected['chart']))
                target=reference.ctx.mpf(phase);allow=reference.ctx.mpf('1e-100')
                assert left<=target+allow and right>=target-allow,(p2,phase)
                scalar=reference.evaluate(phase);primitives=kernel.primitives(box,selected['chart'])
                for key,rkey in (('A','A'),('B_over_Pstar','B')):
                    assert contains(current.bounded(primitives[key])+tolerance,c.mpf(mp.nstr(scalar[rkey],135))),(p2,phase,key)
                    primitive_checks+=1
                if phase in ('0','.5','1'):assert all(v.zero for v in primitives.values())
                got_density=current.CellDensityGraph(graph,query,primitives,None,7).values()['densities']
                p=reference.ctx;E=reference.Utheta;V=p.mpf('.27')
                before=current.density.recovery.history_densities(E,V)
                after=current.density.recovery.history_densities(E*p.exp(scalar['A']/7),V+scalar['B']/7)
                for key in KEYS:
                    assert contains(current.bounded(got_density[key])+tolerance,c.mpf(mp.nstr(after[key]-before[key],135))),(p2,phase,key)
                    density_checks+=1
                phase_checks+=1
            cases+=1
        assert branches=={'small_r_series','signed_Mobius'}
        # Conditional log-u intervals retain the same original q and broad p2
        # ranges. Independent reference points verify both signs/coordinates.
        conditional_checks=0
        for sign in (-1,1):
            reference=original.GenericShearLoop(scales,a='.8',b=0,p1=5,p2=sign,Utheta='1.3')
            p2_cover=('.0001','5') if sign>0 else ('-5','-.0001')
            logfull=c.ln(c.mpf(('.0001','5'))*c.mpf(mp.nstr(reference.q,135))/c.mpf('.25'))
            lo,hi=ep(logfull);margin=c.ln(2)/100
            covers=[c.mpf((lo,ep(c.ln(c.mpf('.25'))-margin)[0])),
                    c.mpf((ep(c.ln(c.mpf('.125'))+margin)[1],hi))]
            for cover in covers:
                query=finite_fixture(c,reference,p2_cover=p2_cover,logu_cover=cover);kernel=query['kernel']
                assert kernel.geometry in ('small_r_series','signed_Mobius')
                for magnitude in ('.0001','.02','.05','.2','2','5'):
                    point=original.GenericShearLoop(scales,a='.8',b=0,p1=5,p2=sign*scales.ctx.mpf(magnitude),Utheta='1.3')
                    logvalue=point.ctx.ln(abs(point.u));cl,ch=ep(cover)
                    if not cl<=logvalue<=ch:continue
                    for phase in ('.137','.663'):
                        got=kernel.evaluate(phase,bits=70);selected=got['selected_inverse']
                        primitive=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                        target=point.evaluate(phase)
                        for key,rkey in (('A','A'),('B_over_Pstar','B')):
                            assert contains(current.bounded(primitive[key])+tolerance,c.mpf(mp.nstr(target[rkey],135)))
                        conditional_checks+=1
        assert conditional_checks>12
    return dict(passed=True,independent_original_scalar_cases=cases,
        directed_inverse_bracket_endpoint_comparisons=phase_checks,
        independent_original_A_B_comparisons=primitive_checks,
        independent_changed_minus_original_five_density_comparisons=density_checks,
        conditional_source_partition_scalar_comparisons=conditional_checks,
        opposite_p2_signs_regular_signed_and_exact_zero_checked=True,
        reference_roundoff_allowed='1e-90 primitive/density;1e-100 phase',
        finite_fixtures_not_native_source_or_integral_oracles=True)


def source_contract(c,row):
    contract=row['local_basis_contract'];q=row['native_positive_q']
    assert contract['order']==current.LOCAL_ORDER
    assert contract['original_root_fourth_powers_all_zero'] and contract['original_roots_unchanged_by_rebinding_zero_powers']
    assert contract['q_fourth_power_exactly_zero'] and contract['one_shared_basis_and_ledger_for_this_cell']
    assert contract['original_phase_math_unchanged_except_equivalent_u_source_binding']
    assert not contract['standard_packet_basis_compatible']
    assert contract['normalized_finite_density_ranges_exported_before_cross_cell_sum']
    assert not q['exact_zero'] and q['sign']=='positive' and not q['point_value_selected']
    assert ep(interval(c,q['coefficient_interval']))==(1,1)
    scale=q['formal_positive_scale'];assert scale['source_exponents']==[0,0,0,0] and scale['radius_power']==0
    assert same(interval(c,scale['additional_log_interval']),interval(c,row['whole_native_positive_logq_interval']))
    assert row['q_identity']['passed'] and row['original_positive_q_not_replaced_by_flat_or_field_point']
    assert not row['inverse_or_global_controls_installed']
    pieces=row['conditioned_source_piece_geometry'];split=row['overlapping_source_coordinate_split']
    assert all(p['branch'] in ('small_r_series','signed_Mobius') and not p['flat'] for p in pieces)
    if split:
        assert contract['strict_original_p2_sign'] in (-1,1)
        assert split['entire_native_logabsu_source_range_covered'] and split['conditional_original_u_equals_p2_q_over_dstar_on_each_piece']
        covers=[interval(c,v) for v in split['piece_logabsu_intervals']]
        full=interval(c,row['conditioned_geometry']['original_u']['log_absolute_upper'])
        assert ep(covers[0])[0]==ep(full)[0] and ep(covers[-1])[1]==ep(full)[1]
        assert len(covers)==len(pieces)==2 and ep(covers[0])[1]>=ep(covers[1])[0]
        assert ep(covers[0])[1]<ep(c.ln(c.mpf('.25')))[0]
        assert ep(covers[1])[0]>ep(c.ln(c.mpf('.125')))[1]
        assert [p['branch'] for p in pieces]==['small_r_series','signed_Mobius']
    return len(pieces),bool(split)


def native_checks(c,manifest,owner):
    mass_checks=0;source_checks=0;inverse_checks=0;split_cells=0
    refinements={key:[] for key in KEYS};widths={key:[] for key in KEYS}
    p=mp.mp.clone();p.dps=c.dps+40
    with mp.workdps(c.dps+40):
        for result in manifest['actual_original_nonzero_Z_own_integral_refinements']:
            count=result['ordered_source_cells'];rows=result['whole_source_phase_density_records']
            level,phase=owner.parent.levels[count];assert len(rows)==count
            assert result['source_family']==manifest['source_family']==owner.family
            assert result['original_Z_exact']=='37/100' and result['exact_y_window']==['0','1']
            assert result['explicit_candidate_N']==phase['explicit_candidate_N']==7
            assert result['own_rates']==current.five.RATES and result['normalized_own_units']==current.five.UNITS
            assert result['original_P0_datum_sha256']==owner.family['datum_enclosure_sha256']
            assert result['incoming_five_histories_and_P0_not_reset'] and result['contribution_only_no_default_actual_incoming_defects']
            assert result['complete_original_source_window_and_true_phase_unions_integrated']
            assert all(not result[k] for k in ('current_whole_N_selected','actual_five_controls_installed','Z_functional_terminal_matching_installed'))
            total={key:c.mpf(0) for key in KEYS}
            for index,(row,original_phase) in enumerate(zip(rows,phase['whole_source_cells'],strict=True)):
                source=row['source'];piece_count,split=source_contract(c,source)
                assert source['source_family']==owner.family and source['ordered_source_cell_level']==count
                assert source['ordered_source_cell_index']==index
                assert source['exact_y_cell']==original_phase['exact_y_cell']
                assert source['exact_Z_range']==['37/100','37/100']
                assert source['original_source_function_coefficient_record']['source_family']==owner.family
                node=owner.log_nodes(count);expected=current.ordered.hull(c,node[index+1],node[index])
                assert same(expected,interval(c,source['whole_native_positive_logq_interval']))
                assert same(interval(c,source['original_selected_eta_log']),owner.owner.scales.logs['eta'])
                inverse_rows=row['actual_native_inverse_phase_pieces']
                phases=original_phase['true_common_N_phase_boxes']
                assert len(inverse_rows)==piece_count*len(phases)
                for i,proof in enumerate(inverse_rows):
                    assert proof['source_coordinate_piece_index']==i//len(phases)
                    assert proof['true_source_phase']==phases[i%len(phases)]
                    assert proof['local_logabsu_basis_order']==current.LOCAL_ORDER
                    selected=proof['selected_original_inverse']
                    assert selected['bracket_proof'] in ('directed endpoint inequalities and exact strict source monotonicity','exact periodic/half-period symmetry')
                    assert contains(interval(c,selected['phase_image']),interval(c,proof['true_source_phase']))
                    assert all(v['encloses_original_source_function'] and not v['point_value_selected'] for v in proof['native_primitives'].values())
                    inverse_checks+=1
                assert row['native_B_original_V_and_all_cross_terms_retained'] and row['finite_enclosures_summed_after_leaving_local_u_factor_basis']
                left,right=p.mpf(index)/count,p.mpf(index+1)/count
                for key in KEYS:
                    weight=interval(c,row['positive_own_rate_final_endpoint_masses'][key]);wl,wh=ep(weight)
                    rate=p.mpf(current.five.RATES[key])
                    exact=right-left if rate==0 else (p.exp(-rate*(1-right))-p.exp(-rate*(1-left)))/rate
                    assert 0<wl<=exact<=wh
                    hull=interval(c,row['five_signed_normalized_density_hulls'][key])
                    contribution=interval(c,row['five_signed_own_rate_contributions'][key])
                    assert same(hull*weight,contribution)
                    total[key]+=contribution;mass_checks+=1
                source_checks+=1;split_cells+=split
            for key in KEYS:
                saved=interval(c,result['five_signed_original_nonzero_Z_own_rate_integral_contributions'][key])
                assert same(total[key],saved)
                refinements[key].append(saved);widths[key].append(ep(saved)[1]-ep(saved)[0])
        for key in KEYS:
            assert max(ep(v)[0] for v in refinements[key])<=min(ep(v)[1] for v in refinements[key])
            assert all(widths[key][i+1]<widths[key][i] for i in range(len(widths[key])-1))
        live_replays=0
        for record in manifest['native_whole_cell_positive_logq_queries']:
            count,index=record['ordered_source_cell_level'],record['ordered_source_cell_index']
            zl,zh=record['exact_Z_range']
            live=owner.query(count,index,Z_lower=zl,Z_upper=zh)
            assert base.encoded(live['record'])==record
            source_contract(c,record)
            if not live['pieces']:
                assert zl=='-1/100' and zh=='1/100'
                assert live['kernel'].geometry=='requires_signed_source_refinement'
                assert record['local_basis_contract']['strict_original_p2_sign']==0
                continue
            for part in live['pieces']:
                kernel=part['kernel']
                for phase in ('.137','.663'):
                    selected=kernel.evaluate(phase,bits=24)['selected_inverse']
                    lo,hi=ep(selected['coordinate_interval']);target=c.mpf(phase)
                    if lo>0:assert ep(kernel.phase(c.mpf(lo),selected['chart']))[1]<ep(target)[0]
                    if hi<1:assert ep(kernel.phase(c.mpf(hi),selected['chart']))[0]>ep(target)[1]
                    primitive=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                    assert all(v.scale.bases is part['q'].scale.bases and v.ledger is part['ledger'] for v in primitive.values())
                    assert all(all(mp.isfinite(x) for x in ep(current.bounded(v))) for v in primitive.values())
                    live_replays+=1
    return dict(passed=True,complete_original_window_source_cell_checks=source_checks,
        independent_positive_closed_form_Duhamel_mass_checks=mass_checks,
        actual_common_N_phase_union_inverse_records_checked=inverse_checks,
        native_conditional_small_signed_source_split_cells=split_cells,
        fresh_extreme_source_inverse_endpoint_and_primitive_checks=live_replays,
        common_refinement_intersection_and_decreasing_widths=True,
        refinement_interval_widths=widths,
        finest_five_normalized_contribution_intervals={k:v[-1] for k,v in refinements.items()},
        signed_Z_crossing_remains_explicitly_unresolved=True,
        no_density_Z_integral_or_functional_terminal_matching_claim=True)


def logarithmic_formula_checks(c):
    p=mp.mp.clone();p.dps=c.dps+40;comparisons=0
    with mp.workdps(c.dps+40):
        for eta in ('.01','1e-30'):
            for y in ('0','1/100','1/3','1/2','2/3','99/100','1'):
                rational=sy.Rational(y);point=p.mpf(int(rational.p))/int(rational.q)
                sigma=original.flat_step(p,point);a=p.mpf(4)/5+p.mpf(6)/5*sigma
                q2=(p.mpf(3)/5*(1-sigma)+p.mpf(eta))/a
                value=current.endpoint_logq(c,rational,c.ln(c.mpf(eta)))
                # Independent direct source arithmetic is used only at finite
                # scales; a declared reference tolerance covers its rounding.
                target=c.mpf(mp.nstr(p.ln(q2)/2,c.dps+20))
                assert contains(value+c.mpf(('-1e-150','1e-150')),target)
                comparisons+=1
        assert current.q_identity()['passed']
    return dict(passed=True,independent_finite_original_q_formula_comparisons=comparisons,
        original_endpoint_identity_and_exact_source_monotonicity_checked=True,
        original_constructor_AST_binding_guard_executed=True,
        unchanged_original_bounded_callback_method_code_objects=current.SPECIALIZED_METHODS,
        only_constructor_original_u_binding_and_equivalent_signed_B_specialized=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE] and manifest['one_fixed_nonzero_Z_complete_original_window_integral_enclosed']
    assert not manifest['standard_packet_basis_compatible'] and all(manifest[k] is False for k in FLAGS)
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.OriginalO2PositiveLogQCells();c=owner.c
    formula=logarithmic_formula_checks(c)
    print('Original positive q formula and AST binding PASS',flush=True)
    finite=scalar_checks(c,owner.owner.scales.graph)
    print('Independent original phase, A/B and signed five densities PASS',flush=True)
    native=native_checks(c,manifest,owner)
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        positive_native_logq_and_original_source_binding_checks=formula,
        independent_original_scalar_phase_AB_and_signed_density_checks=finite,
        complete_fixed_nonzero_Z_native_integral_contracts=native,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            uncompressed_bytes=len(raw),compressed_bytes=path.stat().st_size),
        one_fixed_nonzero_Z_complete_original_window_integral_enclosed=True,
        standard_packet_basis_compatible=False,**dict.fromkeys(FLAGS,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Whole-window native original O2 C0 five signed own-rate contribution enclosures at fixed Z37/100 and candidate N7. Positive logq, local conditional logabsu, original inverse/A and exact equivalent signed B checked. Not Z-functional matching, incoming chart assembly, global N/controls, stress, recursion or full NS.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Native positive logq, conditional original source and fixed nonzero-Z five integrals PASS',flush=True)
    return result


if __name__=='__main__':run()

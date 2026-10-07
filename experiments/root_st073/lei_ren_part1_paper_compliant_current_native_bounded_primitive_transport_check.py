"""Independent C0 primitive supports and original bounded-range transport."""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_bounded_primitive_transport as current
import lei_ren_part1_paper_compliant_current_native_cutoff_density_oracle_check as density_checks
import lei_ren_part1_paper_compliant_current_native_phase_first_jets_check as reference
import lei_ren_part1_paper_compliant_current_native_Rc_range_loss as loss

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets=current.packets;ep=current.ep;require=reference.require


def independent_support_checks(c):
    original=reference.original
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    p=scales.ctx;count=0
    for a in ('.8','1.8'):
        for b in ('-.2','.2'):
            aa,bb=p.mpf(a),p.mpf(b);Delta=aa+bb*bb/aa-2
            for p2_value in ('-.4','0','.4'):
                loop=original.GenericShearLoop(scales,a=aa,b=bb,p1=5,p2=p2_value,Utheta='1.3')
                for fraction in ('.01','.137','.49','.91','1'):
                    point=loop.at_angle(2*p.pi*p.mpf(fraction))
                    require(abs(point['A'])<=p.mpf(5)/4+p.mpf('1e-85'),'Independent original A support failed')
                    require(abs(point['B']/loop.Utheta)<=p.mpf(3)/2+p.mpf('1e-85'),'Independent original B/E support failed')
                    count+=2
    p2=mp.mp.clone();p2.dps=100;spectral=0
    for r in ('-.73','0','.73'):
        rr=p2.mpf(r);w=lambda x:(p2.cos(x)-rr)/(1-2*rr*p2.cos(x)+rr*rr)
        parts=[p2.pi*k/4 for k in range(9)]
        mean=p2.quad(w,parts)/(2*p2.pi);square=p2.quad(lambda x:w(x)**2,parts)/(2*p2.pi)
        require(abs(mean)<p2.mpf('1e-85') and abs(square-1/(2*(1-rr*rr)))<p2.mpf('1e-85'),
            'Independent signed/small Poisson mean-square normalization failed');spectral+=2
    bases=tuple(c.mpf(0) for _ in range(5));ledger=reference.qchecks.new_ledger()
    scalar=lambda value:current.cutoff.prior.ScaledEnclosure(current.cutoff.prior.FormalScale(bases),value,ledger)
    factor=lambda log,coefficient:current.cutoff.prior.ScaledEnclosure(current.cutoff.prior.FormalScale(bases,offset=c.mpf(log)),coefficient,ledger)
    tiny=factor('-1e6',(-1,1));bound=scalar((-1.25,1.25))
    got,row=current.signed_C0_support_range(tiny,bound)
    require(got is tiny and row['original_tighter_formal_range_retained'],'Tiny original A factor was replaced by scalar support cap')
    broad=factor('1e6',(-1,1));got,row=current.signed_C0_support_range(broad,bound)
    require(not row['original_tighter_formal_range_retained'] and ep(got.finite_interval())==ep(c.mpf((-1.25,1.25))),
        'Broad original A enclosure not support-bounded')
    E=factor('-1e6',1);B=E*factor('1e6',(-1,1));got,row=current.signed_C0_support_range(B,E*c.mpf((-1.5,1.5)))
    require(not got.zero and got.scale.powers==E.scale.powers and got.scale.offset._mpi_==E.scale.offset._mpi_,
        'B support product lost the original tiny E factor')
    zero=scalar(0);got,_=current.signed_C0_support_range(zero,bound)
    require(got is zero,'Exact original primitive zero changed')
    for sign in (-1,1):
        try:current.signed_C0_support_range(factor('1e6',sign),bound)
        except ArithmeticError:pass
        else:raise AssertionError('Contradictory sign-definite source range silently support-capped')
    return dict(passed=True,independent_original_A_B_support_comparisons=count,
        independent_signed_Poisson_mean_square_comparisons=spectral,
        tighter_tiny_A_original_object_and_B_E_factor_retained=True,contradictory_sign_definite_ranges_rejected=2,
        synthetic_fixtures_not_original_field_values=True,original_reference_precision_dps=100)


@current.native.inlet.source_precision
def run(role_owner,report,live):
    began=time.monotonic();owner,got=live
    require(report[current.GATE] and report['source_family']==owner.family,'Same original bounded primitive producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed bounded primitive prerequisite: '+name)
    require(report['enclosed_original_whole_cells']==24 and report['full24_original_C1_integral_range_transport_enclosed'],
        'Original full24 continuous source/integral ranges regressed')
    comparison=report['comparison_with_checked_original_fixed_N_target_ranges']
    require(all(comparison[key]['strict_absolute_upper_reduction'] for key in current.accepted.current.repair.ROWS),
        'All five original C0 target upper ranges must strictly improve')
    require(report['strict_target_absolute_upper_reductions']==8,'Actual target improvement count differs')
    support_rows=0;restricted=0;retained=0;live_rows=0
    for cell in got['cells']:
        for group in ('values','Z_derivatives'):
            for value in cell[group].values():
                require(value.ctx is owner.ctx and value.ledger is owner.coordinates.ledger,'Original transport context/ledger changed');live_rows+=1
        if cell['frame'] is None:continue
        source=cell['frame'].record['actual_original_spatial_source']
        for row in source['original_C0_A_B_source_support_ranges_before_nonlinear_density']:
            require(row['original_A_Z_and_B_Z_same_objects'] and row['B_support_retains_original_E_factor'] and
                row['original_E_E_Z_V_V_Z_and_all_cross_terms_retained'],'Original derivative or source cross row changed')
            for key in ('A','B_over_Pstar'):
                support_rows+=1
                if row[key]['original_tighter_formal_range_retained']:retained+=1
                else:restricted+=1
    for row in report['original_serial_cells']:
        if row['chart']=='O3_power':
            require(row['incoming_C0']['p']==row['outgoing_C0']['p'] and row['incoming_Z']['p']==row['outgoing_Z']['p'],
                'Original quiet pressure memory changed')
    require(not any(report.get(key) for key in packets.OPEN) and not report['useful_repair_contraction_or_actual_controls_established'],
        'Improved range enclosures cannot admit controls/global N/closure/recursion')
    inventory=json.loads((HERE/loss.NAME).read_bytes())
    require(inventory[loss.GATE] and inventory['source_family']==owner.family and inventory['original_continuous_cells']==24,
        'Same original dominant range-loss inventory required')
    require(inventory['original_A_enclosure_exceeds_known_support'],'Original support loss frontier missing')
    independent=independent_support_checks(owner.ctx)
    formulas=density_checks.independent_density_checks(owner.oracle)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=2048,
        actual_continuous_original_route_cells=24,actual_live_integral_C0_Z_rows=live_rows,
        branchwise_original_C0_support_ranges_checked=support_rows,source_theorem_restricted_C0_ranges=restricted,
        tighter_original_formal_C0_ranges_retained=retained,
        strict_original_C0_target_upper_reductions=5,strict_original_Z_target_upper_reductions=3,
        independent_original_primitive_support_checks=independent,independent_original_density_formula_checks=formulas,
        original_A_Z_B_Z_and_E_V_rows_unchanged=True,original_quiet_pressure_memory_preserved=True,
        completed_actual_controls_global_N_closure_or_recursion=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name),
            loss.NAME:sha(loss.NAME),Path(loss.__file__).name:sha(Path(loss.__file__).name)},
        scope='Original A/B C0 support proofs, signed/tiny factor guards and actual full24-cell source-function transport with five C0 and three Z strict target upper reductions. Original derivatives remain uncapped; controls, global N, closure and recursion remain open.')
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original C0 primitive support check PASS;5 C0 +3 Z target ranges improved',flush=True)
    return result

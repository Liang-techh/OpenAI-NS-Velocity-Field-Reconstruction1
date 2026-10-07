"""Independent rebase/affine transport identities and actual original replay."""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_C1_history_transfer as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks
import lei_ren_part1_paper_compliant_current_native_q_slow_jets_check as qchecks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets
require=checks.require;contains=checks.contains;prior=current.prior;density=current.density


def independent_transfer_checks(c):
    p=mp.mp.clone();p.dps=c.dps+30;family={key:'independent-fixture-'+key for key in packets.FAMILY_KEYS}
    coordinates=current.CommonSourceCoordinates(c,c.mpf('.22'),family);count=0
    for sign in (-1,1):
        old_bases=tuple(c.mpf(v) for v in ('-.3','.22','-.7','1.3','2.1'));ledger=qchecks.new_ledger()
        value=prior.ScaledEnclosure(prior.FormalScale(old_bases,(.5,-.5,1,-.5,.5),offset='.17'),sign,ledger)
        rebased=coordinates.rebase(value,family)
        reference=sign*p.exp(p.mpf('-.3')/2-p.mpf('.22')/2-p.mpf('.7')-p.mpf('1.3')/2+p.mpf('2.1')/2+p.mpf('.17'))
        require(contains(rebased.finite_interval(),c.mpf(p.nstr(reference,c.dps+20))),'Independent signed source rebase failed')
        require(rebased.scale.powers==(0,-.5,0,0,0),'Only same Pstar squared symbolic exponent may remain')
        require(rebased.scale.bases is coordinates.bases and rebased.ledger is coordinates.ledger,'Common rebase must attach one arithmetic basis/ledger')
        count+=1
    # Native bridge width bases carry a legacy context. Copy complete covers,
    # retaining the same accepted common-Pstar interval, rather than rejecting
    # or selecting one width endpoint.
    foreign=mp.ctx_iv.MPIntervalContext();foreign.dps=c.dps+20
    bases=(foreign.mpf(('-.31','-.29')),foreign.mpf(coordinates.logP_squared),foreign.mpf(0),foreign.mpf(0),foreign.mpf(0))
    value=prior.ScaledEnclosure(prior.FormalScale(bases,(1,0,0,0,0)),foreign.mpf(('1.1','1.2')),qchecks.new_ledger())
    copied=coordinates.rebase(value,family)
    require(copied.ctx is c and contains(copied.finite_interval(),c.mpf(value.finite_interval())),'Whole legacy-context source cover not retained')
    tiny=prior.ScaledEnclosure(prior.FormalScale(tuple(c.mpf(v) for v in ('-1e40','.22','0','0','0')),(1,0,0,0,0)),1,qchecks.new_ledger())
    require(not coordinates.rebase(tiny,family).zero,'Unmaterializable rebased source became zero')
    decay=coordinates.decay(c.mpf('1e40'),1)
    require(not decay.zero and ep(decay.coefficient)[0]>0,'Huge quiet-interval positive decay became zero')
    require(contains(coordinates.decay(c.mpf('1e40'),0).finite_interval(),c.mpf(1)),'Rate0 pressure memory was damped')
    for bad_family in ({**family,packets.FAMILY_KEYS[0]:'foreign'},{}):
        try:coordinates.rebase(tiny,bad_family)
        except ValueError:pass
        else:raise ArithmeticError('Foreign source family accepted by common rebase')
    bad_bases=tuple(c.mpf(v) for v in ('0','.23','0','0','0'))
    try:coordinates.rebase(prior.ScaledEnclosure(prior.FormalScale(bad_bases),1,qchecks.new_ledger()),family)
    except ValueError:pass
    else:raise ArithmeticError('Different common Pstar squared source cover accepted')
    # Actual affine source integration: two nonquiet cells and one quiet cell.
    # f_j(t,Z)=sign_j*(2+t)*(1+Z+Z²), with genuine nonzero incoming memory.
    operator=current.C1DuhamelOperator(coordinates);zbox=c.mpf(('.49','.51'));widths=('.13','.17','.11')
    left=p.mpf(0);analytic_parts={key:p.mpf(0) for key in current.RATES}
    coefficients={key:(-1 if index%2 else 1)*(index+1) for index,key in enumerate(current.RATES)}
    before_quiet=None
    for index,width in enumerate(widths):
        w=p.mpf(width);right=left+w;values={};jets={}
        for key,rate in current.RATES.items():
            r=p.mpf(rate.numerator)/rate.denominator
            mass=density.local.positive_kernel_mass(c,c.mpf(width),rate)
            if index==2:
                values[key]=coordinates.scalar(0);jets[key]=coordinates.scalar(0)
                analytic_parts[key]*=p.exp(-r*w)
            else:
                radial=c.mpf((p.nstr(2+left,c.dps+20),p.nstr(2+right,c.dps+20)))
                values[key]=coordinates.scalar(coefficients[key]*radial*(1+zbox+zbox**2)*mass)
                jets[key]=coordinates.scalar(coefficients[key]*radial*(1+2*zbox)*mass)
                exact_mass=w if not r else -p.expm1(-r*w)/r
                first_moment=w*w/2 if not r else (1-(1+r*w)*p.exp(-r*w))/(r*r)
                analytic_parts[key]=p.exp(-r*w)*analytic_parts[key]+coefficients[key]*((2+right)*exact_mass-first_moment)
        if index==2:before_quiet=(operator.increments['p'].record(),operator.Z_increments['p'].record())
        operator.append(c.mpf(width),values,jets,family);left=right
    require(operator.steps==3,'Three affine cells required')
    require(packets.encode(operator.increments['p'].record())==packets.encode(before_quiet[0]) and
        packets.encode(operator.Z_increments['p'].record())==packets.encode(before_quiet[1]),'Quiet pressure increment or Z memory reset')
    incoming={key:coordinates.scalar((index+2)*(1+zbox*c.mpf('.2'))) for index,key in enumerate(current.RATES)}
    incoming_Z={key:coordinates.scalar(c.mpf(index+2)*c.mpf('.2')) for index,key in enumerate(current.RATES)}
    got=operator.apply(incoming,incoming_Z,family);comparisons=0
    for point in ('.49','.5','.51'):
        z=p.mpf(point)
        for index,(key,rate) in enumerate(current.RATES.items()):
            decay=p.exp(-p.mpf(rate.numerator)/rate.denominator*left)
            references=(decay*(index+2)*(1+p.mpf('.2')*z)+analytic_parts[key]*(1+z+z*z),
                decay*(index+2)*p.mpf('.2')+analytic_parts[key]*(1+2*z))
            for value,reference in zip((got['values'][key],got['Z_derivatives'][key]),references):
                require(contains(value.finite_interval(),c.mpf(p.nstr(reference,c.dps+20))),'Independent signed affine C1 integration failed: '+key)
                comparisons+=1
    return dict(passed=True,independent_signed_formal_rebase_comparisons=count,
        independent_analytic_nonzero_incoming_piecewise_C0_Z_transfer_comparisons=comparisons,
        legacy_context_whole_cover_copy_checked=True,unmaterializable_rebase_and_quiet_decay_not_zeroed=True,
        rate0_pressure_C0_Z_memory_retained=True,foreign_family_and_common_Pstar_cover_rejected=True,
        scalar_fixtures_are_not_native_incoming_field_values=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['adjacent_radial_cell_count']==3 and saved['original_inlet_Z_query_count']==2,'Actual inlet and three adjacent cells required')
    require(not any(saved.get(k) for k in packets.OPEN),'Local affine transfer cannot complete global stages')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed native C1 transfer prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeC1HistoryTransfer(density.NativeDensityC1LocalIntegrals(density.first.NativePhaseFirstJets(density.slow.NativeQSlowJets(density.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge)))))))
        independent=independent_transfer_checks(owner.ctx);inlet_rows=0;serial_rows=0;regions={}
        for name,old in saved['actual_original_inlet_C1_history_records'].items():
            Z=packets.interval(owner.ctx,old['source_provenance']['Z_box']);got=owner.inlet(Z,saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual original inlet functions changed: '+name)
            require(got['geometry']['offset'].zero and len(got['geometry']['phase_boxes'])==1 and ep(got['geometry']['phase_boxes'][0])==(0,0),'Actual inlet must have exact zero original radius offset/phase')
            require(all(value.zero for value in got['initial_defect'].values()) and all(value.zero for value in got['initial_defect_Z'].values()),'Original checked inlet correction initial rows must be zero')
            functions=got['functions'];require(set(functions['originals'])==set(current.RATES) and set(functions['Z_derivatives'])==set(current.RATES),'Five original incoming C0/Z histories required')
            require(old['P0_not_merged_into_pressure_history'] and old['native_Pstar_width_or_radial_factor_not_reapplied'],'Incoming pressure/unit contract changed')
            for value in [*functions['originals'].values(),*functions['Z_derivatives'].values(),functions['P0'],functions['P0_Z'],functions['pressure'],functions['pressure_Z']]:
                require(value.scale.bases is owner.coordinates.bases and value.ledger is owner.coordinates.ledger,'Original incoming function common basis changed')
            require(not any(old.get(k) for k in packets.OPEN),'Inlet functions cannot complete global stages')
            inlet_rows+=14;regions[name]=dict(original_history_C0_Z_rows=10,separate_P0_C0_Z_rows=2,original_absolute_pressure_C0_Z_rows=2)
            print('Actual original inlet C1 memory checked:',name,flush=True)
        for name,old in saved['actual_adjacent_C1_serial_operator_records'].items():
            endpoints=[str(q['numerator'])+'/'+str(q['denominator']) for q in old['exact_adjacent_coordinate_partition']]
            got=owner.serial(Z=packets.interval(owner.ctx,old['Z_box']),endpoints=endpoints,N=saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual adjacent native C1 operator changed: '+name)
            require(Fraction(**old['exact_total_log_radius_width_fraction'])==Fraction(3,100),'Original expanded width.03 must be exact')
            require(got['operator'].steps==3 and len(old['actual_whole_cell_C1_signed_contributions'])==3,'Three true adjacent native cells required')
            for key,coefficient in got['operator'].coefficients.items():
                reference=owner.ctx.exp(-owner.ctx.mpf('0.03')*owner.ctx.mpf(current.RATES[key].numerator)/current.RATES[key].denominator)
                require(contains(coefficient.finite_interval(),reference),'Actual composed true-width decay changed: '+key)
            require(ep(got['operator'].coefficients['p'].finite_interval())==(1,1),'Native pressure decay must be exact1')
            require(old['actual_incoming_correction_at_local_left_not_computed_or_assumed_zero'] and old['missing_actual_inlet_to_local_gap_not_treated_as_quiet'],
                'Local operator cannot manufacture a global correction history')
            require(not old['global_inlet_to_Rc_histories_or_repair_admitted'] and not any(old.get(k) for k in packets.OPEN),'Local affine operator cannot admit global gates')
            require(old['original_right_endpoint_background_history_and_P0_Z']['P0_not_merged_into_pressure_history'],'Native original right background/P0 contract changed')
            serial_rows+=10;regions[name+'_serial']=dict(adjacent_cells=3,total_true_width='3/100',C0_Z_affine_increment_rows=10,
                full_phase_union_cells=sum(cell['original_whole_cell_signed_density_Z_source']['actual_original_radius_phase']['periodic_projection']['full_period'] for cell in old['actual_whole_cell_C1_signed_contributions']))
            print('Actual adjacent C1 signed operator checked:',name,flush=True)
        for points in (['.2','.1'],['.1','.1','.2'],['-.1','.2'],[owner.ctx.mpf('.1'),'.2']):
            try:owner.serial(Z=('.49','.51'),endpoints=points,N=1024)
            except ValueError:continue
            raise ArithmeticError('Invalid/non-exact adjacent native partition accepted')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=saved['candidate_N'],
        actual_inlet_C0_Z_and_pressure_rows_checked=inlet_rows,actual_inlet_Z_queries_checked=2,
        actual_adjacent_radial_cells_per_partition_checked=3,serial_Z_queries_checked=2,
        actual_cumulative_signed_C0_Z_affine_increment_rows_checked=serial_rows,
        actual_local_radial_coordinate_interval_checked=['.12','.15'],actual_true_width='3/100',
        expanded_radial_width_relative_to_previous_cell=1500000,
        independent_rebase_and_C1_Duhamel_checks=independent,regions=regions,
        original_incoming_memory_and_separate_P0_Z_retained=True,local_incoming_correction_not_invented=True,
        global_cumulative_history_or_Rc_or_repair_or_common_N_admission=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual original inlet/adjacent C1 PASS:',inlet_rows,'inlet rows;',serial_rows,'affine increment rows',flush=True)
    return result


if __name__=='__main__':run()

"""Check signed whole-cell integrals and actual native C0 function replay."""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_local_signed_integrals as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets
require=checks.require;contains=checks.contains


def independent_integration_checks(c):
    p=mp.mp.clone();p.dps=c.dps+20;count=0
    for rate in (Fraction(0),Fraction(1),Fraction(3,2)):
        r=p.mpf(rate.numerator)/rate.denominator
        for width in ('.00000002','.3'):
            w=p.mpf(width);mass=current.positive_kernel_mass(c,c.mpf(width),rate)
            true_mass=w if not rate else -p.expm1(-r*w)/r
            require(contains(mass,c.mpf(p.nstr(true_mass,c.dps+15))),'Original exact positive mass not enclosed')
            for sign in (-1,1):
                # Independent exact affine-source integration. f(s)=sign*(2+s).
                # u=w-s gives integral exp(-r*u)*(2+w-u)du.
                first=w*w/2 if not rate else (1-(1+r*w)*p.exp(-r*w))/(r*r)
                reference=sign*((2+w)*true_mass-first)
                box=c.mpf((2,ep(2+c.mpf(width))[1]))*sign
                enclosure=box*mass
                require(contains(enclosure,c.mpf(p.nstr(reference,c.dps+15))),'Signed analytic linear-source integral not enclosed')
                count+=1
    tiny=c.mpf(p.make_mpf((0,1,-10000,1)))
    mass=current.positive_kernel_mass(c,tiny,Fraction(3,2))
    require(ep(mass)[0]>0 and contains(mass/tiny,c.mpf(1)),'Microscopic positive cell mass dropped')
    require(ep(current.positive_kernel_mass(c,c.mpf('.3'),0))==ep(c.mpf('.3')),'Pressure kernel must have exact rate0 mass')
    require(ep(current.positive_kernel_mass(c,0,1))==(0,0),'Zero width mass must be exactly zero')
    for w,r in ((-1,1),(1,-1)):
        try:current.positive_kernel_mass(c,w,r)
        except ValueError:continue
        raise ArithmeticError('Invalid Duhamel mass accepted')
    bases=tuple(c.mpf(0) for _ in range(5));ledger={}
    values=[current.prior.ScaledEnclosure(current.prior.FormalScale(bases,offset=-2),(-2,-1),ledger),
            current.prior.ScaledEnclosure(current.prior.FormalScale(bases),(-3,-2),ledger)]
    hull=current.same_source_union(values).finite_interval()
    require(contains(hull,c.mpf(-3)) and contains(hull,-c.exp(c.mpf(-2))),'Signed phase-cell union hull lost a source range')
    require(ep(hull)[1]<0,'Uniformly negative phase-cell union incorrectly widened across zero')
    require(current.same_source_union([values[0]]) is values[0],'Single phase cell lost its original factored source')
    return dict(passed=True,independent_positive_mass_comparisons=6,independent_signed_affine_integral_comparisons=count,
        signed_phase_union_hull_checked=True,microscopic_positive_mass_retained=True,
        pressure_mass_rate_zero_exact=True,zero_width_mass_exact=True,invalid_mass_inputs_rejected=True,
        scalar_fixtures_not_native_field_values=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_radial_cell_count']==1 and saved['native_Z_query_count']==2,'Actual local native cell records required')
    require(not any(saved.get(k) for k in packets.OPEN),'Local integrals cannot complete global stages')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed local integral prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        binder=current.spatial.NativeSpatialPhase(current.native.NativeGenericSourcePackets(bridge))
        density=current.density.NativeCandidateDensities(current.phase.NativeConditionedPhase(current.spatial.current.NativeCorrelatedShearQ(current.prior.NativeSignedInputEnclosures(binder.native))))
        owner=current.NativeLocalSignedIntegrals(density,binder);independent=independent_integration_checks(owner.ctx)
        regions={};signs=dict(m='negative',h='negative',k='negative',e='positive',p='negative')
        for name,old in saved['actual_native_local_signed_integral_records'].items():
            Z=packets.interval(owner.ctx,old['Z_box']);got=owner.contribution(Z=Z,left='.13369999',right='.13370001',N=saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual local signed integral source changed: '+name)
            require(set(got['contributions'])==set(current.RATES),'Five original signed contributions required')
            require(got['record']['entire_radial_and_Z_cell_source_covers_used'] and got['record']['actual_phase_union_not_samples'],'Sample values cannot define actual integrals')
            require(got['record']['incoming_history_not_assumed_or_reset'] and got['record']['local_contributions_are_not_global_defect_histories'],'Local contribution overstates a global history')
            require(got['record']['C0_integral_function_enclosures_only'] and not got['record']['actual_C1_Z_derivatives_or_Rc_targets_installed'],'C0 source range cannot invent C1 Z/Rc targets')
            require(not any(got['record'].get(k) for k in packets.OPEN),'Local integral cannot admit global stages')
            require(Fraction(**got['record']['exact_log_radius_width_fraction'])==Fraction(1,50000000),'Original cell width not retained exactly')
            for key,value in got['contributions'].items():
                require(not value.zero and value.record()['sign']==signs[key],'Actual nonzero signed integral not retained: '+key)
                require(ep(got['record']['original_positive_Duhamel_kernel_mass_covers'][key])[0]>0,'Original positive kernel mass required')
            regions[name]=dict(Z_box=got['record']['Z_box'],whole_source_cell_covered=True,
                exact_coordinate_width='1/50000000',all_five_nonzero_integral_signs=signs)
            print('Actual local signed integral functions checked:',name,flush=True)
        for left,right,chart in (('.2','.1','O2_slope'),(owner.ctx.mpf('.1'),'.2','O2_slope'),('.1','.2','reshape')):
            try:owner.contribution(Z=('.5','.5'),left=left,right=right,N=1024,chart=chart)
            except ValueError:continue
            raise ArithmeticError('Unsupported or non-exact local coordinate request accepted')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},native_radial_cell_count=1,
        native_Z_query_count=2,signed_actual_local_integral_enclosures_checked=10,all_ten_signed_integrals_nonzero=True,
        actual_Z_interval_covered=['.49','.51'],candidate_N=saved['candidate_N'],regions=regions,
        independent_signed_integration_checks=independent,
        global_cumulative_history_or_C1_Z_or_Rc_or_common_N_admission=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual local signed Duhamel integrals PASS:10 nonzero enclosures',flush=True)
    return result


if __name__=='__main__':run()

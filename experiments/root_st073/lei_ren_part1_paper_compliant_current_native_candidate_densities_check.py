"""Native density replay and independent original velocity/moment formulas."""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_candidate_densities as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as phase_checks

phase=current.phase;HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets
require=phase_checks.require;contains=phase_checks.contains


def scalar_density_checks(c):
    original=phase_checks.original
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    count=0
    for a,b,p1,p2 in (('.8','.2','5','.2'),('.8','-.2','5','-.2'),('.8','.2','5','0'),('2.019','0','5','0')):
        ref=original.GenericShearLoop(scales,a=a,b=b,p1=p1,p2=p2,Utheta='1.3')
        source,scalar=phase_checks.fixture(c,a,b,p2,mp.nstr(ref.q,105))
        loop=phase.ConditionedPhase(source,c.ln(c.mpf(mp.nstr(scales.d_star,105))))
        for N in (1,1024):
            for phi in ('.137','.663'):
                got=current.candidate_at_phase(loop,scalar('-.3'),phi,N)
                candidate=ref.modulate(logR_offset=ref.ctx.mpf(phi)/N,N=N,Uz='-.3',slow_A_y=0,slow_B_y=0)
                moments=original.moment_increment_densities(ref.ctx,R='.5',Utheta='1.3',Uz='-.3',
                    delta_theta=candidate['delta_theta'],delta_z=candidate['delta_z'])
                # At R=.5, sqrt(2R)=1; the pressure dR density is twice
                # the dimensionless Duhamel kernel. S=1 in this fixture.
                reference=dict(m=moments['Mz'],h=moments['Mtheta'],k=moments['Mztheta'],e=moments['M2'],p=moments['Mp']/2)
                for key,rkey in (('deltaE','delta_theta'),('deltaV','delta_z'),('candidate_E','Utheta'),('candidate_V','Uz')):
                    covered=phase.bounded_value(got['values'][key])+c.mpf(('-1e-85','1e-85'))
                    require(contains(covered,c.mpf(mp.nstr(candidate[rkey],105))),'Original finite-N velocity not enclosed: '+key)
                for key,value in reference.items():
                    covered=phase.bounded_value(got['densities'][key])+c.mpf(('-1e-85','1e-85'))
                    require(contains(covered,c.mpf(mp.nstr(value,105))),'Original physical moment density not enclosed: '+key)
                count+=1
    require(count==16,'Expected signed/small/flat density reference cases')
    for N in (True,False,0,-1,'1024',1.5):
        try:current.candidate_integer(N)
        except ValueError:continue
        raise ArithmeticError('Invalid candidate integer accepted')
    # Even when x cannot be materialized, expm1(x) keeps x's own scale.
    source,scalar=phase_checks.fixture(c,2,0,(-2,3),1)
    tiny=current.prior.ScaledEnclosure(current.prior.FormalScale(source['q'].scale.bases,offset='-1e40'),1,source['q'].ledger)
    for sign in (-1,1):
        value=tiny*sign;expm1=current.factored_expm1(value)
        require(not expm1.zero and expm1.scale.record()==value.scale.record(),'Tiny expm1 factor was dropped')
        require(contains(expm1.coefficient/value.coefficient,c.mpf(1)),'Tiny exprel must enclose1')
    source['q']=tiny;loop=phase.ConditionedPhase(source,c.ln(c.mpf('.005')))
    got=current.candidate_at_phase(loop,scalar('-.3'),'.137',1024)
    require(not any(v.zero for v in got['densities'].values()),'Positive tiny q cannot erase all original signed density functions')
    return dict(passed=True,independent_original_finite_N_velocity_and_five_physical_density_comparisons=count,
        nonzero_original_axial_velocity_cross_terms_checked=True,opposite_signed_r_and_p2_zero_and_flat_checked=True,
        signed_tiny_expm1_scale_retained=True,tiny_five_density_kernels_not_zeroed=True,
        invalid_candidate_N_rejected=True,reference_roundoff_allowed='1e-85',fixtures_are_not_native_field_values=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_source_box_count']==5,'Five native candidate density boxes required')
    require(not any(saved.get(k) for k in packets.OPEN),'Candidate density backend cannot complete global stages')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed density prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeCandidateDensities(phase.NativeConditionedPhase(current.current.NativeCorrelatedShearQ(current.prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge)))))
        independent=scalar_density_checks(owner.ctx);regions={};count=0
        for chart,record in saved['actual_candidate_velocity_and_density_records'].items():
            p=record['original_source']['source_provenance'];Z=packets.interval(owner.ctx,p['Z_box']);coordinate=packets.interval(owner.ctx,p['coordinate_box'])
            source,loop,V=owner.source(chart,Z,coordinate)
            require(packets.encode(source['record'])==record['original_source'],'Native source changed: '+chart)
            require(packets.encode(V.record())==record['original_normalized_V'],'Native common-unit original V changed: '+chart)
            active=0
            for phi,old in record['candidate_queries'].items():
                got=current.candidate_at_phase(loop,V,phi,saved['candidate_N'])
                require(packets.encode(got['record'])==old,'Native candidate density changed: '+chart+' '+phi)
                require(got['record']['status']=='enclosed','Native candidate density unavailable')
                require(set(got['densities'])==set(current.RATES),'All five original signed kernels required')
                require(not any(got['record'].get(k) for k in packets.OPEN),'A local candidate cannot admit a global stage')
                require(got['record']['candidate_N_is_not_global_common_N_admission'] and not got['record']['spatial_phase_binding_installed'],'Candidate parameters cannot become original spatial phase')
                if phi in ('0','.5','1') or loop.flat:
                    require(all(v.zero for v in got['densities'].values()),'Flat/symmetric phase density increments must be exactly zero')
                else:active+=1
                if chart=='O2_buffer' and phi=='.137':
                    require(not got['values']['deltaE'].zero and not got['values']['deltaV'].zero,'Tiny buffer velocity increments incorrectly zeroed')
                    require(not any(v.zero for v in got['densities'].values()),'Tiny buffer signed density kernels incorrectly zeroed')
                count+=1
            regions[chart]=dict(candidate_queries=len(record['candidate_queries']),nontrivial_active_queries=active,
                original_V_enclosure_recovered_in_same_units=True)
            print('Live candidate velocity/five densities checked:',chart,flush=True)
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},
        native_source_box_count=len(regions),candidate_velocity_and_five_density_queries_checked=count,
        signed_density_kernel_enclosures_checked=count*5,candidate_N=saved['candidate_N'],
        independent_original_velocity_and_moment_density_checks=independent,regions=regions,
        spatial_phase_binding_or_moment_integral_or_common_N_admission=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Live native candidate velocities/five densities PASS:',count,'queries;',count*5,'kernel enclosures',flush=True)
    return result


if __name__=='__main__':run()

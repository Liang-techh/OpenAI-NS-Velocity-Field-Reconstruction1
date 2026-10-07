"""Original finite-N velocity increments and five signed density kernels.

The phase and N are explicit candidate parameters. Source functions are the
same native field used by the conditioned inverse. These kernels must still
be integrated with actual spatial phase, inlet histories and original lengths.
"""
import json
import operator
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_conditioned_phase as phase

current=phase.current;prior=phase.prior;packets=phase.packets;native=phase.native
HERE,PREFIX,sha=phase.HERE,phase.PREFIX,phase.sha;ep=phase.ep
NAME=PREFIX+'current_native_candidate_densities.json'
RECEIPT=PREFIX+'current_native_candidate_densities_check.json'
GATE='current_native_candidate_velocity_and_five_signed_density_C0_functions_executed'
RATES=dict(m=1,h='3/2',k='3/2',e=1,p=0)


def candidate_integer(N):
    if isinstance(N,bool):raise ValueError('Explicit positive integer candidate N required')
    try:N=operator.index(N)
    except TypeError as exc:raise ValueError('Explicit positive integer candidate N required') from exc
    if N<1:raise ValueError('Explicit positive integer candidate N required')
    return N


def factored_expm1(value):
    """expm1(x)=x*integral_0^1 exp(t*x)dt, with x's tiny factor retained."""
    if value.zero:return value
    c=value.ctx;x=phase.bounded_value(value);lo,hi=ep(x)
    if max(abs(lo),abs(hi))>1:raise ArithmeticError('Candidate N must keep this source exponent in [-1,1]')
    lower=ep(c.exp(c.mpf(min(lo,0))))[0];upper=ep(c.exp(c.mpf(max(hi,0))))[1]
    return value*c.mpf((lower,upper))


def signed_density_kernels(E,V,deltaE,deltaV):
    theta_cross=E*deltaE;theta_square=current.square(deltaE)
    axial_square=current.square(deltaV)
    return dict(m=deltaV,h=deltaE,
        k=V*deltaE+E*deltaV+deltaE*deltaV,
        e=V*deltaV*2+axial_square-theta_cross-theta_square*E.ctx.mpf('.5'),
        p=theta_cross+theta_square*E.ctx.mpf('.5'))


def candidate_at_phase(loop,V,phi,N):
    N=candidate_integer(N);result=loop.evaluate(phi)
    if result['status']!='enclosed':
        return dict(record=dict(status=result['status'],phase_result=result,candidate_N=N),values=None)
    selected=result['selected_inverse'];primitive=loop.primitives(selected['coordinate_interval'],selected['chart'])
    factor=loop.c.mpf(1)/N
    deltaE=loop.E*factored_expm1(primitive['A']*factor)
    deltaV=primitive['B_over_Pstar']*factor
    densities=signed_density_kernels(loop.E,V,deltaE,deltaV)
    values=dict(original_E=loop.E,original_V=V,deltaE=deltaE,deltaV=deltaV,
        candidate_E=loop.E+deltaE,candidate_V=V+deltaV)
    record=dict(status='enclosed',candidate_N=N,explicit_free_phase_parameter=result['phase'],
        original_phase_inverse_and_primitives=result,
        original_and_candidate_normalized_velocity_enclosures={k:v.record() for k,v in values.items()},
        five_signed_dimensionless_Duhamel_density_kernels={k:v.record() for k,v in densities.items()},
        original_Duhamel_rates=RATES,log_factored_expm1_preserves_tiny_source_increment=True,
        original_axial_velocity_and_all_cross_terms_retained=True,
        original_absolute_P0_unchanged=True,spatial_phase_binding_installed=False,
        candidate_N_is_not_global_common_N_admission=True,
        density_integrals_or_new_moment_histories_constructed=False,**dict.fromkeys(packets.OPEN,False))
    return dict(record=record,values=values,densities=densities,primitives=primitive)


class NativeCandidateDensities:
    def __init__(self,owner):
        if type(owner) is not phase.NativeConditionedPhase:raise ValueError('Same live conditioned phase owner required')
        admitted=json.loads((HERE/phase.RECEIPT).read_bytes())
        if not admitted['all_passed'] or admitted['source_family']!=owner.family:raise ValueError('Accepted same conditioned source family required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service
        self.service.bind_hashes(admitted['input_hashes']);self.service.bind_hashes({phase.RECEIPT:sha(phase.RECEIPT),Path(__file__).name:sha(Path(__file__).name)})
    def source(self,chart,Z,coordinate):
        source,loop=self.owner.query(chart,Z,coordinate);packet=source['packet']
        leaf=prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:packet.velocity['axial'][0]})
        V=self.owner.owner.owner.leaf(leaf,loop.E.scale.bases,loop.E.ledger)
        return source,loop,V


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeCandidateDensities(phase.NativeConditionedPhase(current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge)))))
        saved=json.loads((HERE/phase.NAME).read_bytes())['actual_C0_phase_inverse_and_primitive_records'];records={};N=1024
        for chart,old in saved.items():
            p=old['source']['source_provenance'];Z=packets.interval(owner.ctx,p['Z_box']);coordinate=packets.interval(owner.ctx,p['coordinate_box'])
            source,loop,V=owner.source(chart,Z,coordinate)
            values={phi:candidate_at_phase(loop,V,phi,N)['record'] for phi in old['candidate_phase_queries']}
            records[chart]=dict(original_source=source['record'],original_normalized_V=V.record(),candidate_queries=values)
            print('Live candidate velocity/five densities:',chart,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},native_source_box_count=len(records),candidate_N=N,
        actual_candidate_velocity_and_density_records=records,original_Duhamel_rates=RATES,
        explicit_candidate_phase_and_N_not_global_admission=True,spatial_phase_binding_installed=False,
        actual_changed_defect_integral_functions_installed=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='C0 actual source E_N,V_N and five signed dimensionless Duhamel density kernels on five original boxes at free candidate phases and N=1024. No actual spatial phase, moment histories/integrals, repair, common N/cone, recursion or full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

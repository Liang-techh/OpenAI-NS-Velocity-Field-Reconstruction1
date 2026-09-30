"""Bounded same-profile backward pressure for the source outer candidate.

Returns P/Pstar^2, avoiding huge amplitudes and radii. The omitted half-line
is bounded using s<=-1/2 beyond y=1, monotone Z-flattening, H<=1, and the
terminal heat normalization factor 1/(1-epsilon). No completed-core pressure
or Section 7 pressure-restoration claim follows from this temporary reference
extension. The source absolute constants are still uncertified.
"""
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
from numpy.polynomial.legendre import leggauss


def normalized_outer_pressure(schedule,Z,*,cutoff_y=64.,order=48):
    z=float(Z); end=Decimal(str(cutoff_y))
    if abs(z)>1 or not (0<schedule.mu<1 and 0<schedule.epsilon<1):
        raise ValueError('Require |Z|<=1, 0<mu<1 and 0<epsilon<1')
    if not (end>=1 and end<schedule.y_tail):
        raise ValueError('Cutoff must be between y=1 and the pre-heat tail')
    edges={Decimal(0),Decimal(1),end}
    for value in schedule.checkpoint_offsets.values():
        if 0<value<end: edges.add(value)
        if 0<value+1<end: edges.add(value+1)
    edges=sorted(edges)
    nodes,weights=leggauss(order)
    def normalized_squared(y):
        with localcontext() as context:
            context.prec=schedule.decimal_precision
            row=schedule.at_log_radius(schedule.logRref+Decimal(str(y)),z)
        return math.exp(2*float(row['log_angular_amplitude']-schedule.logPstar))
    integral=0.
    for a,b in zip(edges[:-1],edges[1:]):
        af,bf=float(a),float(b)
        integral+=(bf-af)/4*sum(float(w)*normalized_squared(
            af+(bf-af)*(float(x)+1)/2) for x,w in zip(nodes,weights))
    bound=.5*normalized_squared(float(end))/(1-float(schedule.epsilon))**2
    reference_inside=2.5/(1+z*z)**2
    return {'Z':z,'P_at_Rref_over_Pstar_squared':-integral,
            'reference_extension_P0_over_Pstar_squared':-reference_inside-integral,
            'omitted_outer_pressure_upper_bound_over_Pstar_squared':bound,
            'omitted_tail_only_pressure_interval_at_Rref':[-integral-bound,-integral],
            'quadrature_error_included_in_interval':False,
            'tail_bound_numeric_evaluation':'finite floating-point evaluation of the analytic bound formula',
            'cutoff_y':float(end),'order':order,
            'reference_extension_regular_at_axis':False,
            'completed_core_pressure_bound':False,'outer_moments_closed':False}


def run():
    from lei_ren_part1_paper_outer import PaperOuterSchedule
    schedule=PaperOuterSchedule(logPstar=14,logRref=10,delta='1e-32',Md='.5',
        c_mu='.001',c_delta='.001',c_epsilon='.01')
    rows=[]
    for z in (-.5,0.,.5):
        coarse=normalized_outer_pressure(schedule,z,order=32)
        fine=normalized_outer_pressure(schedule,z,order=64)
        rows.append({'coarse':coarse,'fine':fine,
                     'quadrature_absolute_difference':abs(coarse['P_at_Rref_over_Pstar_squared']-
                                                         fine['P_at_Rref_over_Pstar_squared'])})
    report={'source':'https://arxiv.org/html/2609.35406v1','equation':'6.6',
            'schedule':schedule.metadata(),'rows':rows,
            'scope':'Backward pressure of the uncorrected source outer plus temporary singular reference extension; normalized units, finite quadrature and analytic omitted-tail bound only.',
            'common_regular_core_constructed':False,'relaxed_cone_certified':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,default=str)+'\n',encoding='utf-8')
    print(json.dumps(rows),flush=True)
    return report


if __name__=='__main__': run()

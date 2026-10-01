"""Cancellation-safe angular stress and full-cell long-reshape relaxed cone.

Axial shear is exactly zero on this layer, so its relaxed cone margin depends
only on angular total stress. No axial inertial stress is invented or reset.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_long_reshape_field import (
    IntervalLongReshapeField,T_EXACT,RZ_LOG_EXACT,reshape_kernel,sigma_value_derivative)
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_jet,_pack
from lei_ren_part1_paper_interval_exit_switch_enclosure import zero_axial_source_cone
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent

def normalized_angular(c,z,delta,R,V,mass_source,rho,source_n,A,zeta,J,st):
    """Pure angular algebra; source_n uses physical Z derivatives / Ctheta0.

    mass_source is Mz0/R0 as a C1 jet; rho=R0/R. source_n keys are theta,
    theta_Z, mixed, mixed_Z. These derivatives include the source amplitude
    derivative, rather than the derivative of an independently divided ratio.
    """
    zv=z[0];d=1-zv*zv;L=1-delta*zv*zv
    mass=V+(mass_source-V)*rho
    H=-1+(1-delta)*zv*mass[0]+d*mass[1]
    source_combo=((1-delta/2)*source_n['theta']-(1-delta)*zv*source_n['theta_Z']/2
        -d*source_n['mixed_Z']+(2*delta-1)*zv*source_n['mixed'])*A
    C0=1-delta/2-d*V[1]+(2*delta-1)*zv*V[0]
    Cz=(1-delta)*zv/2+d*V[0]
    # Group the shared zeta*J before interval evaluation.
    increment_combo=J[0]*(C0-Cz*zeta)-Cz*J[1]
    inertial=R*(H+source_combo+increment_combo)/L
    return dict(I_theta_over_F=inertial,T_theta_over_F=inertial+st,
        S_theta_over_F=st,S_z_over_F=c.mpf(0),transport_over_R=H,
        source_angular_combo=source_combo,increment_angular_combo=increment_combo,
        coefficient_C0=C0,coefficient_Cz=Cz)

class IntervalLongReshapeAngularCone:
    def __init__(self,field=None):
        self.field=field if field is not None else IntervalLongReshapeField()
        c=self.field.calc.ctx
        with mp.workdps(self.field.calc.precision+60):
            state=self.field.source['normalized_exit_state']
            phi=restore_jet(c,state['phi'],order=1)
            theta=restore_jet(c,state['theta'],order=1)
            mixed=restore_jet(c,state['theta_z'],order=1)
            factor=self.field.calc.eps**2/(2*c.mpf(110)**2*phi[0])
            ell=self.field.calc.ell[0]
            self.n=dict(theta=theta[0]*factor,theta_Z=(theta[1]+ell*theta[0])*factor,
                mixed=mixed[0]*factor,mixed_Z=(mixed[1]+ell*mixed[0])*factor)
            self.mass_source=self.field.m0['z']/110
    def evaluate_cell(self,left,right):
        a,b=Fraction(left),Fraction(right)
        if not 0<=a<=b<=RZ_LOG_EXACT:raise ValueError('cell must lie between R110 and Rz')
        field=self.field;c=field.calc.ctx
        with mp.workdps(field.calc.precision+60):
            lo=c.mpf(a.numerator)/a.denominator;hi=c.mpf(b.numerator)/b.denominator
            y=c.mpf([endpoints(lo)[0],endpoints(hi)[1]])
            R=110*c.exp(y);rho=c.exp(-y)
            sig,ds=sigma_value_derivative(c,y/field.T)
            zeta=-((field.z*2)/(1+field.z*field.z))[0]+field.B[1]*(1-sig)
            J=reshape_kernel(c,c.mpf('1.6'),1,field.B,field.T,y)
            A=c.exp(-c.mpf('1.6')*y+field.B[0]*sig)
            st=-c.mpf('.8')-2*field.B[0]*ds/field.T
            data=normalized_angular(c,field.z,field.calc.delta,R,field.V,self.mass_source,
                rho,self.n,A,zeta,J,st)
            stress=dict(normalized_stress=data,cone=dict(status='not certified: source prerequisites'))
            certificate=zero_axial_source_cone(field.calc,stress,-st)
            return dict(left=str(a),right=str(b),normalized_angular_stress=data,cone=certificate,
                common_amplitude_cancelled=True,source_normalized_data_retained=True,
                pressure_and_axial_inertial_stress_unmodified=True,
                axial_inertial_stress_evaluated=False,whole_radial_cell=True)

def run():
    evaluator=IntervalLongReshapeAngularCone()
    points=[evaluator.evaluate_cell(y,y) for y in (0,1,T_EXACT/2,T_EXACT,RZ_LOG_EXACT)]
    cells=[evaluator.evaluate_cell(0,T_EXACT),evaluator.evaluate_cell(T_EXACT,RZ_LOG_EXACT)]
    complete=all(r['cone'].get('relaxed_cone_certified',False) for r in cells)
    result=dict(center_family=evaluator.field.calc.center_family,state_sha256=evaluator.field.calc.state_hash,
        points=points,cells=cells,full_R110_to_Rz_relaxed_cone_certified=complete,
        axial_shear_identically_zero=True,weak_margin_independent_of_axial_inertial_stress=True,
        physical_velocity_moments_and_pressure_changed=False,whole_axis=False,
        strong_admissibility_claimed=False,axial_restore_whole_path_cone_certified=False,
        original_parameter_remainders_enclosed=False,heat_exterior_matched=False,temporal_recursion=False)
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_long_reshape_field.py',
        evaluator.field.source_name,evaluator.field.defects_name,
        'lei_ren_part1_paper_interval_exit_switch_enclosure.py')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Normalized angular point cone:',[r['cone']['status'] for r in points],flush=True)
    print('Whole R110-to-Rz relaxed cone:',complete,'cell statuses:',[r['cone']['status'] for r in cells],flush=True)
    return result

if __name__=='__main__':run()

"""Actual reference/axial restore mixed logR/Z derivatives through4.

Postrestoration ends at Rm; the actual five-bump patch takes over there.
All positive amplitude sources and the inherited analytic P0 are retained.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_reference_restore_profiles import (
    CompliantReferenceRestoreProfiles,IntervalTaylor,relative_amplitude_jet)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets,positive_exp
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square,derivative
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def binomial_rate(rows,rate,k):
    return sum((rows[j]*(math.comb(k,j)*rate**(k-j)) for j in range(k+1)),rows[0]*0)


def source_mixed(c,Z,delta,E,alpha,centered,logu,p0,invP2):
    """Repeated original y=logR source equations and true product rules."""
    if len(alpha)!=5 or E.order!=5 or any(row.order!=5 for row in centered.values()):
        raise ValueError('Actual axial5 histories and cutoff y orders0..4 required')
    z=IntervalTaylor.variable(c,c.mpf(Z),5);zero=z*0;one=zero+1
    mismatch=[E*value for value in alpha];V=[z*4+mismatch[0]]+mismatch[1:]
    rows={name:[value] for name,value in centered.items()}
    for k in range(4):
        rows['mean_error'].append(mismatch[k]-rows['mean_error'][k])
        rows['angular_error'].append(rows['angular_error'][k]*c.mpf('-1.6'))
        rows['mixed_error'].append(mismatch[k]-rows['mixed_error'][k]*c.mpf('1.6'))
        quadratic=sum((mismatch[j]*mismatch[k-j]*math.comb(k,j) for j in range(k+1)),zero)
        rows['axial_square'].append(quadratic-rows['axial_square'][k])
        rows['swirl_error'].append(rows['swirl_error'][k]*c.mpf('-1.2'))
        rows['pressure_error'].append(rows['pressure_error'][k]*c.mpf('-.2'))
    H=[rows['angular_error'][k]+(c.mpf('.625') if k==0 else 0) for k in range(5)]
    K=[(z*H[k])*4+rows['mixed_error'][k] for k in range(5)]
    m=[rows['mean_error'][k]+(z*4 if k==0 else zero) for k in range(5)]
    A=[(z*rows['mean_error'][k])*8+rows['axial_square'][k]+(square(z)*16 if k==0 else zero) for k in range(5)]
    b=[rows['swirl_error'][k]+(c.mpf(5)/6 if k==0 else 0) for k in range(5)]
    p=[rows['pressure_error'][k]+(5 if k==0 else 0) for k in range(5)]
    L=1-square(z)*delta;d=1-square(z)
    Q=[(2*z*V[k]-(z*m[k])*(1-delta)-d*derivative(m[k]))/L for k in range(5)]
    amp=relative_amplitude_jet(logu);amp2=relative_amplitude_jet(logu,2)
    ratio0=positive_exp(c,2*logu[0]);ratio=amp2*ratio0
    physical=dict(Utheta_over_current_Utheta=[amp*c.mpf('.1')**k for k in range(5)],
        Uz=V,Ur_over_current_sqrt_R_over_2=[binomial_rate(Q,c.mpf('.5'),k) for k in range(5)],
        P_over_Pstar2=[p0+ratio*p[0]/2]+[ratio*binomial_rate(p,c.mpf('.2'),k)/2 for k in range(1,5)])
    primitives=dict(Mtheta_over_current_sqrt2_R_1p5_Utheta=[amp*binomial_rate(H,c.mpf('1.6'),k) for k in range(5)],
        Mtheta_z_over_current_sqrt2_R_1p5_Utheta=[amp*binomial_rate(K,c.mpf('1.6'),k) for k in range(5)],
        Mz_over_current_R=[binomial_rate(m,c.mpf(1),k) for k in range(5)],
        Mztheta_over_current_R_Pstar2=[binomial_rate(A,c.mpf(1),k)*invP2-ratio*binomial_rate(b,c.mpf('1.2'),k)/2 for k in range(5)],
        Mp_over_Pstar2=[ratio*binomial_rate(p,c.mpf('.2'),k)/2 for k in range(5)])
    grid=lambda values:{'y'+str(k)+'_Z'+str(n):row[n]*math.factorial(n)
        for k,row in enumerate(values) for n in range(5-k)}
    return dict(physical_velocity_pressure_y_Z_mixed4={name:grid(values) for name,values in physical.items()},
        physical_five_primitive_y_Z_mixed4={name:grid(values) for name,values in primitives.items()},
        physical_velocity_pressure_y_derivative_Taylor={name:[row.truncate(4-k) for k,row in enumerate(values)] for name,values in physical.items()},
        actual_centered_moment_y_derivative_axial5=rows,actual_Q_y_derivative_axial4=Q,
        actual_alpha_ordinary_y_derivatives=alpha,actual_Utheta_squared_over_Pstar_squared_enclosure=ratio0,
        exact_positive_swirl_source_log=logu[0],positive_exponential_cap_is_enclosure_only=True,
        radial_prefactors_differentiated_before_mixed_grid=True)


class CompliantReferenceRestoreMixedC4:
    def __init__(self):
        self.reference=CompliantReferenceRestoreProfiles();self.ctx=c=self.reference.ctx
        self.family=self.reference.family;self.source=self.reference.source;self.hashes=dict(self.reference.hashes)
        for part in ('reference_restore_profiles_check','actual_patch_mixed_C4_check'):
            name=PREFIX+part+'.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or receipt['actual_five_defect_family_sha256']!=self.family or receipt['implicit_source_sha256']!=self.source:
                raise ValueError('Actual reference/patch prerequisite receipt required')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Reference mixed source changed: '+path)
            self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.invP2=c.exp(-2*self.reference.core.logP)
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def packet(self,Z,parent,alpha,chart):
        c=self.ctx;jet=lambda name:IntervalTaylor(c,parent[name])
        centered={name:IntervalTaylor(c,row) for name,row in parent['actual_centered_moment_axial5_coefficients'].items()}
        E=jet('actual_E_V110_minus_4Z_axial5_coefficients');logu=jet('log_Utheta_over_Pstar_axial5_coefficients')
        p0=jet('pressure_axis_axial5_coefficients')
        result=source_mixed(c,Z,self.reference.core.delta,E,alpha,centered,logu,p0,self.invP2)
        return dict(Z=c.mpf(Z),chart=chart,**result,actual_inherited_axial5_packet=parent,
            exact_formal_prefactors=dict(Utheta='current positive Utheta; logUtheta=logPstar+source_log',
                Uz='1',Ur='current sqrt(R/2)',pressure='Pstar^2',y='logR; independent of Z'),
            actual_reference_restore_mixed4_available=True,actual_reference_moment_histories_and_P0_retained=True,
            Rz_restore_end_Rm_functional_joins_certified=True,Rsh_reshape_join_certified=False,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False)

    def reference_branch(self,Z,phase):
        c=self.ctx
        return self.packet(Z,self.reference.reference(Z,phase),[c.mpf(1)]+[c.mpf(0)]*4,'Rsh_to_Rz')

    def restoration(self,Z,t):
        c=self.ctx;cutoff=sigma_jets(c,c.mpf(t))
        alpha=[1-cutoff[0]]+[-cutoff[k]*math.factorial(k) for k in range(1,5)]
        return self.packet(Z,self.reference.restoration(Z,t),alpha,'original_axial_restoration')

    def postrestore(self,Z,offset):
        c=self.ctx;offset=c.mpf(offset)
        if endpoints(offset)[0]<-7 or endpoints(offset)[1]>-6:raise ValueError('Actual postrestore interval ends at Rm before patch')
        return self.packet(Z,self.reference.terminal(Z,offset),[c.mpf(0)]*5,'restore_exit_to_Rm')

    def report(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            whole_reference=self.reference_branch([-1,1],[0,1]),whole_restoration=self.restoration([-1,1],[0,1]),
            whole_postrestore=self.postrestore([-1,1],[-7,-6]),actual_Rsh_exit=self.reference_branch([-1,1],0),
            actual_Rz_reference_side=self.reference_branch([-1,1],1),actual_restore_inlet=self.restoration([-1,1],0),
            actual_restore_exit=self.restoration([-1,1],1),actual_postrestore_inlet=self.postrestore([-1,1],-7),
            actual_Rm_exit=self.postrestore([-1,1],-6),
            interior_packets=[self.restoration(z,t) for z,t in (('0','.25'),('.5','.5'),('0','.75'))],
            actual_reference_restore_mixed4_available=True,Rz_restore_end_Rm_functional_joins_certified=True,
            exact_flat_join_proof='original sigma derivatives vanish at0/1; same actual histories and physical primitive RHSs on both sides; postrestore terminates before the original actual patch',
            exact_positive_source_amplitudes_retained=True,postrestore_unpatched_interval_stops_at_Rm=True,
            Rsh_reshape_join_certified=False,full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantReferenceRestoreMixedC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual reference/restore through Rm mixed logR/Z derivatives through4 generated',flush=True)
    return result


if __name__=='__main__':run()

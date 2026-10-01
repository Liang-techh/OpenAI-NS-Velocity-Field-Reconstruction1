"""Controlled C1 axial restoration and reference continuation Rz..Rm.

Partial cutoff primitives use directed cell integration. Full restoration
weights reuse the same receipt that enters the five source defects.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_long_reshape_field import (
    IntervalLongReshapeField,RZ_LOG_EXACT,RM_LOG_EXACT,positive_decay_integral,sigma_value_derivative)
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value,_pack
from lei_ren_part1_paper_interval_repaired_reference_field import stress_values
from lei_ren_part1_paper_interval_exit_stress_enclosure import cone
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent

def partial_restore_integrals(c,t,cells=128):
    """Endpoint-normalized cutoff moments, for 0<=t<=1."""
    if cells<1 or endpoints(t)[0]<0 or endpoints(t)[1]>1:raise ValueError('invalid restoration quadrature domain')
    values={name:c.mpf(0) for name in ('K1','K16','K2')}
    step=t/cells
    for i in range(cells):
        s=c.mpf([endpoints(step*i)[0],endpoints(step*(i+1))[1]])
        cutoff=alpha_box(c,s+1,c.mpf(1))
        values['K1']+=c.exp(s)*cutoff*step
        values['K16']+=c.exp(c.mpf('1.6')*s)*cutoff*step
        values['K2']+=c.exp(s)*cutoff*cutoff*step
    return dict(K1=values['K1']*c.exp(-t),K16=values['K16']*c.exp(-c.mpf('1.6')*t),K2=values['K2']*c.exp(-t))

class IntervalAxialRestoreField:
    def __init__(self,reshape=None):
        self.reshape=reshape if reshape is not None else IntervalLongReshapeField()
        self.calc=self.reshape.calc;c=self.calc.ctx
        self.start=self.reshape.evaluate_log_offset(RZ_LOG_EXACT)
        self.z=self.reshape.z;self.g=self.reshape.g;self.P0=self.reshape.P0
        self.full={k:restore_value(c,v) for k,v in self.reshape.defects['restoration_constants'].items()}
    def evaluate_log_offset(self,y,cells=128):
        """Unified connecting-field API: y=log(R/110), 0<=y<=log(Rm/110)."""
        q=Fraction(y)
        if not 0<=q<=RM_LOG_EXACT:raise ValueError('connecting offset must lie between R110 and Rm')
        if q<=RZ_LOG_EXACT:return self.reshape.evaluate_log_offset(q)
        return self.evaluate_phase(q-RZ_LOG_EXACT,cells=cells)
    def evaluate_phase(self,t,cells=128):
        q=Fraction(t)
        if not 0<=q<=2:raise ValueError('require 0<=log(R/Rz)<=2')
        c=self.calc.ctx
        with mp.workdps(self.calc.precision+60):
            tt=c.mpf(q.numerator)/q.denominator;R=self.start['R']*c.exp(tt)
            if q>=1:sig,sigprime=c.mpf(1),c.mpf(0)
            else:sig,sigprime=sigma_value_derivative(c,tt)
            if q>=1:
                K=dict(K1=self.full['K1']*c.exp(-tt),K16=self.full['K16']*c.exp(-c.mpf('1.6')*tt),K2=self.full['K2']*c.exp(-tt))
            else:K=partial_restore_integrals(c,tt,cells)
            u=self.start['Utheta']*c.exp(tt/10);uy=u*c.mpf('.1')
            V=self.z*4+self.g*(1-sig);Vy=-self.g*sigprime
            theta=u*(c.sqrt(2)*R**c.mpf('1.5'))*positive_decay_integral(c,c.mpf('1.6'),tt)
            pressure=u*u/2*positive_decay_integral(c,c.mpf('.2'),tt)
            swirl=u*u*(R/2)*positive_decay_integral(c,c.mpf('1.2'),tt)
            deltaR=self.start['R']*(c.exp(tt)-1)
            mass=self.z*(4*deltaR)+self.g*(R*K['K1'])
            mixed=self.z*theta*4+self.g*u*(c.sqrt(2)*R**c.mpf('1.5')*K['K16'])
            axial=self.z*self.z*(16*deltaR)+self.z*self.g*(8*R*K['K1'])+self.g*self.g*(R*K['K2'])
            m0=self.start['physical_moments']
            moments=dict(z=m0['z']+mass,theta=m0['theta']+theta,theta_z=m0['theta_z']+mixed,
                z_theta=m0['z_theta']+axial-swirl,p=m0['p']+pressure)
            P=self.P0+moments['p'];root=c.sqrt(2*R);F=u/root
            stresses=stress_values(c,self.z,self.calc.delta,R,u,uy,V,Vy,moments,P)
            st=-c.mpf('.8');sz=2*Vy[0]/u[0]
            Itheta=stresses['stress']['I_theta'];Iz=stresses['stress']['I_z']
            normalized=dict(S_theta_over_F=st,S_z_over_F=sz,T_theta_over_F=Itheta/F[0]+st,T_z_over_F=Iz/F[0]+sz)
            stresses['normalized_stress']=normalized
            stresses['stress'].update(S_theta=F[0]*st,S_z=F[0]*sz,T_theta=Itheta+F[0]*st,T_z=Iz+F[0]*sz)
            stresses['cone']=cone(c,st,sz,normalized['T_theta_over_F'],normalized['T_z_over_F'])
            zv=self.z[0];L=1-self.calc.delta*zv*zv;mz=moments['z']
            Ur=(2*zv*R*V[0]-(1-self.calc.delta)*zv*mz[0]-(1-zv*zv)*mz[1])/(L*root)
            Ur_R=((1+self.calc.delta)*zv*V[0]+2*zv*Vy[0]-(1-zv*zv)*V[1])/(L*root)-Ur/(2*R)
            return dict(phase=str(q),region='axial_restoration' if q<1 else 'reference_before_five_bump_repair',
                R=R,F=F,Utheta=u,Utheta_y=uy,Uz=V,Uz_y=Vy,P=P,P0=self.P0,
                physical_moments=moments,Ur_value=Ur,Ur_R_value=Ur_R,restoration_integrals=K,
                whole_partial_integrals_enclosed=True,restoration_cells=cells,same_axis_pressure_preserved=True,
                moment_radial_rhs=dict(z=V,theta=u*root,theta_z=u*V*root,z_theta=V*V-u*u/2,p=u*u/(2*R)),
                Ur_Z_available=False,retained_axial_order=1,whole_axis=False,
                whole_path_cone_certified=False,temporal_recursion=False,**stresses)

def run():
    field=IntervalAxialRestoreField();samples=[field.evaluate_phase(t) for t in ('0','.5','1','2')]
    result=dict(center_family=field.calc.center_family,state_sha256=field.calc.state_hash,
        callable_axial_restore_and_reference_continuation=True,samples=samples,
        same_source_and_axis_pressure=True,source_defect_restoration_weights_reused=True,
        complete_R110_to_Rm_field_layers_installed=True,whole_path_cone_certified=False,
        higher_smoothness_or_Ur_Z_certified=False,original_parameter_remainders_enclosed=False,
        heat_exterior_matched=False,temporal_recursion=False)
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_long_reshape_field.py',
        'lei_ren_part1_paper_interval_functional_defects.json','lei_ren_part1_paper_interval_repaired_reference_field.py',
        'lei_ren_part1_paper_interval_comparison_enclosure.py')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Callable axial restoration installed; point cone statuses:',[s['cone']['status'] for s in samples],flush=True)
    return result

if __name__=='__main__':run()

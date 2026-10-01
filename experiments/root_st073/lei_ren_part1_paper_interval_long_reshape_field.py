"""Callable C1 long swirl reshape and reference continuation from R110.

Whole-kernel analytic bounds retain the original moments and axis pressure.
They can be broad during shaping: installing this field is not a cone proof.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_jet,_pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_functional_defects import logjet,abs_upper
from lei_ren_part1_paper_interval_repaired_reference_field import stress_values,validate_receipt
from lei_ren_part1_paper_interval_exit_stress_enclosure import cone
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
T_EXACT=Fraction('4e152')
RM_LOG_EXACT=10*(Fraction('5e151')+14)-6
RZ_LOG_EXACT=RM_LOG_EXACT-2

def positive_decay_integral(c,k,y):
    """Directed integral int_0^y exp(-k ell) dell, including cancellation."""
    if endpoints(k)[0]<=0 or endpoints(y)[0]<0:raise ValueError('positive rate and nonnegative y required')
    value=(1-c.exp(-k*y))/k
    return c.mpf([max(mp.mpf(0),endpoints(value)[0]),max(mp.mpf(0),endpoints(value)[1])])

def sigma_value_derivative(c,s):
    if endpoints(s)[1]<=0:return c.mpf(0),c.mpf(0)
    if endpoints(s)[0]>=1:return c.mpf(1),c.mpf(0)
    if endpoints(s)[0]<=0 or endpoints(s)[1]>=1:
        return 1-alpha_box(c,s+1,c.mpf(1)),c.mpf([0,32])
    sig=1-alpha_box(c,s+1,c.mpf(1))
    derivative=sig*(1-sig)*(2/s**3+2/(1-s)**3)
    # A universal analytic bound, independent of a numerical maximum search.
    lo,hi=endpoints(derivative)
    return sig,c.mpf([max(mp.mpf(0),lo),min(mp.mpf(32),hi)])

def reshape_kernel(c,k,m,B,T,y):
    """Enclose full endpoint-normalized kernel and first axial derivative.

    For B<=0 and 0<=sigma'<=32, exp(-(k+C)ell)<=integrand<=exp(-k ell),
    C=32*m*|B|/T. Its Z derivative is bounded by
    32*m*|B_Z|/T * min(y^2/2,1/k^2). No finite quadrature is used.
    """
    if endpoints(B[0])[1]>0:raise ValueError('negative shape amplitude required')
    C=abs_upper(c,B[0])*(32*m)/T
    lower=positive_decay_integral(c,k+C,y);upper=positive_decay_integral(c,k,y)
    val=c.mpf([endpoints(lower)[0],endpoints(upper)[1]])
    small=y*y/2;large=1/(k*k)
    moment_upper=c.mpf(min(endpoints(small)[1],endpoints(large)[1]))
    derivative_bound=abs_upper(c,B[1])*(32*m)*moment_upper/T
    bound=c.mpf([0,endpoints(derivative_bound)[1]])
    derivative=c.mpf([endpoints(-bound)[0],endpoints(bound)[1]])
    return IntervalTaylor(c,[val,derivative])

def physical_packet(c,z,delta,R,u,uy,V,moments,P0,angular_normalized_shear=None):
    root=c.sqrt(2*R);P=P0+moments['p'];zero=V*0
    zv=z[0];L=1-delta*zv*zv;mz=moments['z']
    Ur=(2*zv*R*V[0]-(1-delta)*zv*mz[0]-(1-zv*zv)*mz[1])/(L*root)
    Ur_R=((1+delta)*zv*V[0]-(1-zv*zv)*V[1])/(L*root)-Ur/(2*R)
    stresses=stress_values(c,z,delta,R,u,uy,V,zero,moments,P)
    if angular_normalized_shear is not None:
        # Cancel the common positive amplitude analytically. Its extremely
        # wide axial-family interval must not turn prescribed shear into zero.
        st=angular_normalized_shear;sz=c.mpf(0);F=u[0]/root
        Itheta=stresses['stress']['I_theta'];Iz=stresses['stress']['I_z']
        tt=Itheta/F+st;tz=Iz/F
        stresses['stress'].update(S_theta=F*st,S_z=sz,T_theta=Itheta+F*st,T_z=Iz)
        stresses['normalized_stress']=dict(S_theta_over_F=st,S_z_over_F=sz,T_theta_over_F=tt,T_z_over_F=tz)
        stresses['cone']=cone(c,st,sz,tt,tz)
    return dict(R=R,F=u/root,Utheta=u,Utheta_y=uy,Uz=V,Uz_y=zero,P=P,P0=P0,
        physical_moments=moments,Ur_value=Ur,Ur_R_value=Ur_R,
        moment_radial_rhs=dict(z=V,theta=u*root,theta_z=u*V*root,z_theta=V*V-u*u/2,p=u*u/(2*R)),
        Ur_Z_available=False,retained_axial_order=1,**stresses)

class IntervalLongReshapeField:
    def __init__(self,calc=None):
        self.calc=calc if calc is not None else IntervalComparisonJets(4)
        c=self.calc.ctx
        self.source_name='lei_ren_part1_paper_interval_exit_switch_enclosure.json'
        self.defects_name='lei_ren_part1_paper_interval_functional_defects.json'
        self.source=json.loads((HERE/self.source_name).read_text())
        self.defects=json.loads((HERE/self.defects_name).read_text())
        for raw,label in ((self.source,'switch'),(self.defects,'defects')):
            validate_receipt(raw,label)
            if raw['state_sha256']!=self.calc.state_hash:raise ValueError('source core mismatch')
        with mp.workdps(self.calc.precision+60):
            self.z=self.calc.z.truncate(1);self.P0=self.calc.p0.truncate(1)
            self.B=restore_jet(c,self.defects['B'],order=1)
            if endpoints(self.B[0])[1]>=0:raise ValueError('strictly negative B required')
            self.g=restore_jet(c,self.defects['g'],order=1);self.V=self.z*4+self.g
            self.m0={k:restore_jet(c,v,order=1) for k,v in self.source['physical_moments'].items()}
            self.u0=restore_jet(c,self.source['F'],order=1)*c.sqrt(220)
            self.I={k:restore_jet(c,v['jet'],order=1) for k,v in self.defects['flat_kernel_bounds'].items()}
            self.log_reference= -logjet(1+self.z*self.z)-c.mpf('5e151')
            self.T=c.mpf(T_EXACT.numerator)/T_EXACT.denominator
    def evaluate_phase(self,phase):
        q=Fraction(phase)
        if not 0<=q<=1:raise ValueError('reshape phase must lie in [0,1]')
        return self.evaluate_log_offset(q*T_EXACT)
    def evaluate_log_offset(self,y):
        q=Fraction(y)
        if not 0<=q<=RZ_LOG_EXACT:raise ValueError('require R110<=R<=Rz')
        c=self.calc.ctx
        with mp.workdps(self.calc.precision+60):
            yy=c.mpf(q.numerator)/q.denominator;R=110*c.exp(yy)
            if q>=T_EXACT:sig,sigprime=c.mpf(1),c.mpf(0)
            else:sig,sigprime=sigma_value_derivative(c,yy/self.T)
            if q==0:u=self.u0;uy=u*c.mpf('.1')
            else:
                logu=self.log_reference+yy/10+self.B*(1-sig)
                u=logu.exp();uy=u*(-self.B*sigprime/self.T+c.mpf('.1'))
            kernels={}
            for name,k,m in (('theta','1.6',1),('pressure','.2',2),('energy','1.2',2)):
                kk=c.mpf(k)
                if q>=T_EXACT:
                    kernels[name]=self.I[name]*c.exp(-kk*(yy-self.T))+positive_decay_integral(c,kk,yy)
                else:kernels[name]=reshape_kernel(c,kk,m,self.B,self.T,yy)
            theta=u*(c.sqrt(2)*R**c.mpf('1.5'))*kernels['theta']
            pressure=u*u*kernels['pressure']/2
            swirl=u*u*(R/2)*kernels['energy']
            mass=self.V*(110*(c.exp(yy)-1))
            moments=dict(z=self.m0['z']+mass,theta=self.m0['theta']+theta,
                theta_z=self.m0['theta_z']+self.V*theta,
                z_theta=self.m0['z_theta']+self.V*self.V*(110*(c.exp(yy)-1))-swirl,
                p=self.m0['p']+pressure)
            st=-c.mpf('.8')-2*self.B[0]*sigprime/self.T
            data=physical_packet(c,self.z,self.calc.delta,R,u,uy,self.V,moments,self.P0,st)
            data.update(log_radius_offset=str(q),region='long_swirl_reshape' if q<T_EXACT else 'reference_before_axial_restore',
                normalized_cumulative_kernels=kernels,shape_cutoff=sig,shape_cutoff_derivative=sigprime,
                full_kernel_and_first_axial_derivative_enclosed=True,same_axis_pressure_preserved=True,
                midpoint_projection=False,whole_axis=False,cone_certificate_claimed=False,
                whole_connecting_field_installed=False,axial_restore_installed=False,temporal_recursion=False)
            return data

def run():
    field=IntervalLongReshapeField()
    samples=[field.evaluate_log_offset(y) for y in (0,1,T_EXACT/2,T_EXACT,RZ_LOG_EXACT)]
    result=dict(center_family=field.calc.center_family,state_sha256=field.calc.state_hash,
        callable_long_reshape_and_reference_continuation=True,
        T=str(T_EXACT),Rm_log_offset=str(RM_LOG_EXACT),Rz_log_offset=str(RZ_LOG_EXACT),
        samples=samples,whole_kernel_enclosed=True,axis_pressure_preserved=True,
        whole_path_cone_certified=False,axial_restoration_complete=False,heat_exterior_matched=False,
        original_parameter_remainders_enclosed=False,temporal_recursion=False)
    names=(Path(__file__).name,field.source_name,field.defects_name,
        'lei_ren_part1_paper_interval_functional_defects.py','lei_ren_part1_paper_interval_repaired_reference_field.py',
        'lei_ren_part1_paper_interval_comparison_enclosure.py')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Callable long reshape and reference continuation installed; point cone statuses:',
        [s['cone']['status'] for s in samples],flush=True)
    return result

if __name__=='__main__':run()

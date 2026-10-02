"""Actual R110 long reshape: axial5 log field and full moment enclosures.

No exp(A), exp(logCstar), huge radius or tiny original amplitude is formed.
Actual moments are normalized by their CURRENT positive Utheta prefactors.
All inherited R110 histories and the analytic pressure datum remain present.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_inner_switch_profiles import (
    CompliantInnerSwitchProfiles, IntervalTaylor, MTH,MZ,MTHZ,MZT,MP)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    logarithm,derivative,square,norm,symmetric,coefficient_lists)
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets,positive_exp
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def bound_exp(c,logvalue,proofs):
    """Bound a factored positive magnitude; cap is never its source value."""
    logvalue=c.mpf(logvalue)
    if endpoints(logvalue)[1]<=-1000:
        proofs.append(dict(log_magnitude_upper=logvalue,log_cap=c.mpf(-1000),
                           exact_source_not_replaced=True))
    return positive_exp(c,logvalue)


def positive_decay(c,rate,y,proofs):
    return bound_exp(c,-c.mpf(rate)*c.mpf(y),proofs)


def decay_integral(c,rate,y,proofs):
    """Integral_0^y exp(-rate*t)dt with exact zero endpoint."""
    y=c.mpf(y);rate=c.mpf(rate)
    if endpoints(rate)[0]<=0 or endpoints(y)[0]<0:raise ValueError('Positive rate and y>=0 required')
    if endpoints(y)==(mp.mpf(0),mp.mpf(0)):return c.mpf(0)
    value=(1-positive_decay(c,rate,y,proofs))/rate
    return c.mpf([max(mp.mpf(0),endpoints(value)[0]),endpoints(value)[1]])


def normalized_amplitude_decay(c,logjet,rate_min,rate_max,y,proofs,initial=None):
    """Enclose exp(logjet), using admitted source rates for its zeroth term.

    Derivative Bell factors are combined with the source decay in LOGS before
    any cap. In particular a tiny exp(-y) is not first capped then multiplied
    by enormous axial derivative factors. At y=0 the exact inlet remains.
    """
    lo,hi=endpoints(c.mpf(y))
    if lo==hi==0:
        return initial if initial is not None else IntervalTaylor.constant(c,1,logjet.order)
    lower=positive_decay(c,rate_max,c.mpf(hi),proofs)
    upper=positive_decay(c,rate_min,c.mpf(lo),proofs)
    zeroth=c.mpf([endpoints(lower)[0],endpoints(upper)[1]])
    nonconstant=IntervalTaylor(c,[c.mpf(0)]+list(logjet.coefficients[1:])).exp()
    if initial is not None:nonconstant=initial*nonconstant
    out=[zeroth if initial is None else initial[0]*zeroth]
    for k in range(1,logjet.order+1):
        factor=norm(c,nonconstant[k])
        if endpoints(factor)[1]==0:out.append(c.mpf(0));continue
        bound=bound_exp(c,-c.mpf(rate_min)*c.mpf(lo)+c.ln(factor),proofs)
        out.append(symmetric(c,bound))
    return IntervalTaylor(c,out)


def backward_kernel(c,B,T,y,m,rate_min,rate_max,proofs):
    """Full, untruncated endpoint-normalized moment integral, axial5.

    K = integral_0^y exp(-k*t + m*B*(sigma((y-t)/T)-sigma(y/T)))dt.
    Source log-u slope in[.05,.15] gives the supplied positive rate bounds.
    |delta sigma|<=8*t/T. Differentiate the EXACT integrand in Z and bound
    its Bell-polynomial coefficient of t^p by full positive exponential
    moments p!/rate_min^(p+1), or the smaller finite-y polynomial mass.
    No quadrature cutoff, frozen sigma or infinite-y replacement defines K.
    """
    c0lo=decay_integral(c,rate_max,y,proofs)
    c0hi=decay_integral(c,rate_min,y,proofs)
    coefficients=[c.mpf([endpoints(c0lo)[0],endpoints(c0hi)[1]])]
    polynomials=[[c.mpf(1)]]
    upper_y=c.mpf(endpoints(c.mpf(y))[1]);rate_min=c.mpf(rate_min)
    for n in range(1,B.order+1):
        row=[c.mpf(0)]*(n+1)
        for k in range(1,n+1):
            factor=norm(c,B[k])*(8*m/T)*k/n
            for p,previous in enumerate(polynomials[n-k]):row[p+1]+=factor*previous
        polynomials.append(row)
        bound=c.mpf(0)
        for p in range(1,n+1):
            finite=upper_y**(p+1)/(p+1)
            complete=c.mpf(math.factorial(p))/rate_min**(p+1)
            mass=c.mpf(min(endpoints(finite)[1],endpoints(complete)[1]))
            bound+=row[p]*mass
        coefficients.append(symmetric(c,bound))
    return IntervalTaylor(c,coefficients)


def relative_amplitude_jet(logjet,m=1):
    return IntervalTaylor(logjet.ctx,[logjet.ctx.mpf(0)]+[q*m for q in logjet.coefficients[1:]]).exp()


class CompliantLongReshapeProfiles:
    def __init__(self):
        self.switch=CompliantInnerSwitchProfiles();self.core=self.switch.core
        self.ctx=c=self.switch.ctx;self.family=self.switch.family;self.source=self.switch.source
        self.hashes=dict(self.switch.hashes);self.records={};self.proofs=[];self.cache={}
        for part in ('inner_switch_profiles_check','reference_join_bounds','reference_join_check','flat_pulse_derivatives_check'):
            name=PREFIX+part+'.json';record=json.loads((HERE/name).read_bytes())
            for path,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Reshape prerequisite changed: '+path)
            self.records[part]=record;self.hashes.update(record['input_hashes'])
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        switch=self.records['inner_switch_profiles_check'];join=self.records['reference_join_bounds']
        if not (switch['all_passed'] and switch['actual_R110_inlet_axial5_available']
                and switch['actual_five_defect_family_sha256']==self.family and switch['implicit_source_sha256']==self.source
                and join['admitted_inner_parameter_family_sha256']==self.switch.bridge.records['K1_ledger']['admitted_inner_parameter_family_sha256']
                and join['uniform_Cstar_family_sha256']==self.core.records['physical_norm_family']['uniform_Cstar_family_sha256']
                and join['same_family_R110_Rh_relaxed_cone_analytically_certified']
                and join['same_preheat_axis_pressure_retained']
                and join['actual_moments_continuously_inherited']
                and self.records['reference_join_check']['reference_join_family_sha256']==join['reference_join_family_sha256']
                and self.records['flat_pulse_derivatives_check']['all_passed']):
            raise ValueError('Actual R110, same selected A/T/Cstar, pressure and cutoff required')
        self.A=read_interval(c,self.core.records['physical_norm_family']['A_upper']);self.T=400*self.A
        if endpoints(read_interval(c,join['B_C2_upper']))!=endpoints(2*self.A):
            raise ValueError('Actual B C2 theorem is not the selected 2A bound')
        for stored,expected in zip(join['shape_shear_interval'],('.7','.9')):
            lo,hi=endpoints(read_interval(c,stored));elo,ehi=endpoints(c.mpf(expected))
            if not lo<=elo<=ehi<=hi:raise ValueError('Original source shear certificate changed')
        step_name='lei_ren_part1_paper_shared_fixed_step_bound.json'
        step=json.loads((HERE/step_name).read_bytes())
        if not (step['global_derivative_bound_certified'] and step['global_derivative_upper']==8
                and hashlib.sha256((HERE/step_name).read_bytes()).hexdigest()==self.hashes[step_name]):
            raise ValueError('Original sigma derivative8 proof is missing')
        radius=read_interval(c,join['B_C2_upper'])*8/self.T
        self.q_source=c.mpf([endpoints(c.mpf('.1')-radius)[0],endpoints(c.mpf('.1')+radius)[1]])
        if endpoints(self.q_source)[0]<endpoints(c.mpf('.05'))[1] or endpoints(self.q_source)[1]>endpoints(c.mpf('.15'))[0]:
            raise ValueError('Kernel rate bounds do not follow from the actual source')
        self.logref=10*(self.core.logC+self.core.logP)
        if endpoints(self.logref-self.T)[0]<=8:
            raise ValueError('Original reshape-to-Rz ordering not admitted')
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def inputs(self,Z):
        c=self.ctx;Z=c.mpf(Z);key=Z._mpi_
        if key in self.cache:return self.cache[key]
        inlet=self.switch.inlet(Z);jet=lambda row:IntervalTaylor(c,row)
        original=jet(inlet['actual_R110_log_shape_axial5_coefficients'])
        # SAME-source C2 theorem bounds the actual B, independently of a
        # broad whole-Z rational interval. Higher B jets retain source data.
        B=IntervalTaylor(c,[intersection(c,q,c.mpf([-endpoints(2*self.A)[1]/math.factorial(k),
             endpoints(2*self.A)[1]/math.factorial(k)])) if k<=2 else q
             for k,q in enumerate(original.coefficients)])
        phi=jet(inlet['F_actual_over_F0_axial5_coefficients']);v=jet(inlet['Uz_actual_axial5_coefficients'])
        raw=inlet['actual_moment_shape_axial5_coefficients']
        moments={n:({k:jet(q) for k,q in row.items()} if n==MZT else jet(row)) for n,row in raw.items()}
        normalized=dict(theta=moments[MTH]/(phi*2),theta_z=moments[MTHZ]/(phi*2),
            pressure=moments[MP]/square(phi),swirl=moments[MZT]['swirl']/square(phi),
            mean=moments[MZ],axial=moments[MZT]['axial'])
        result=dict(inlet=inlet,B=B,phi=phi,v=v,moments=normalized,
            p0=jet(inlet['pressure_axis_axial5_coefficients']),Z=Z)
        self.cache[key]=result;return result

    def evaluate(self,Z,phase=None,y=None):
        """Entire original reshape in phase or exact R110-log offset."""
        c=self.ctx
        if (phase is None)==(y is None):raise ValueError('Supply exactly one phase or R110 log offset')
        if phase is None:
            y=c.mpf(y);phase=y/self.T
        else:
            phase=c.mpf(phase);y=self.T*phase
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Original reshape phase in[0,1] required')
        inp=self.inputs(Z);B=inp['B'];z=IntervalTaylor.variable(c,c.mpf(Z),5)
        cutoff=sigma_jets(c,phase)
        sig=intersection(c,cutoff[0],c.mpf([0,1]));ds=intersection(c,cutoff[1],c.mpf([0,8]))
        logq=logarithm(1+square(z))
        logu=B*(1-sig)-logq+(y/10-self.core.logC-self.core.logP)
        # Use the imported SAME-source shear bound for C0 only. The other
        # axial coefficients retain the actual B derivative source.
        a=B*(2*ds/self.T)+c.mpf('.8')
        a=IntervalTaylor(c,[intersection(c,a[0],c.mpf(['.7','.9']))]+list(a.coefficients[1:]))
        normalized_u=relative_amplitude_jet(logu)
        theta=positive_decay(c,1,y,self.proofs)
        mean=inp['moments']['mean']*theta+inp['v']*(1-theta)
        axial=inp['moments']['axial']*theta+square(inp['v'])*(1-theta)
        kernels={name:backward_kernel(c,B,self.T,y,m,low,high,self.proofs)
                 for name,m,low,high in (('theta',1,'1.55','1.65'),('pressure',2,'.1','.3'),('swirl',2,'1.1','1.3'))}
        decays={name:normalized_amplitude_decay(c,B*(sig*m)-y*k,low,high,y,self.proofs)
                for name,m,k,low,high in (('theta',1,'1.6','1.55','1.65'),('pressure',2,'.2','.1','.3'),('swirl',2,'1.2','1.1','1.3'))}
        # Convolve inherited histories with the derivative Bell factors
        # BEFORE bounding/capping any source exponential magnitude.
        inherited={name:normalized_amplitude_decay(c,B*(sig*m)-y*k,low,high,y,self.proofs,
                         initial=inp['moments'][name])
            for name,m,k,low,high in (('theta',1,'1.6','1.55','1.65'),('theta_z',1,'1.6','1.55','1.65'),
                                     ('pressure',2,'.2','.1','.3'),('swirl',2,'1.2','1.1','1.3'))}
        moment_shapes=dict(theta=inherited['theta']+kernels['theta'],
            theta_z=inherited['theta_z']+inp['v']*kernels['theta'],
            pressure=inherited['pressure']+kernels['pressure'],
            swirl=inherited['swirl']+kernels['swirl'],mean=mean,axial=axial)
        Q=(2*z*inp['v']-(z*mean)*(1-self.core.delta)-(1-square(z))*derivative(mean))/(1-square(z)*self.core.delta)
        amplitude=normalized_u;amplitude2=relative_amplitude_jet(logu,2)
        true_shapes=dict(theta=moment_shapes['theta']*amplitude,theta_z=moment_shapes['theta_z']*amplitude,
            pressure=moment_shapes['pressure']*amplitude2,swirl=moment_shapes['swirl']*amplitude2,
            mean=mean,axial=axial)
        result=dict(Z=c.mpf(Z),phase_enclosure=phase,log_R_over_110_enclosure=y,
            exact_source_radius='110*exp(y); y=400*Abar*phase',
            log_Utheta_over_Pstar_axial5_coefficients=list(logu.coefficients),
            Utheta_true_axial5_divided_by_current_Utheta=list(normalized_u.coefficients),
            Uz_actual_axial5_coefficients=list(inp['v'].coefficients),
            actual_Q_axial4_coefficients=list(Q.coefficients),
            pressure_axis_axial5_coefficients=list(inp['p0'].coefficients),
            angular_shear_axial5_coefficients=list(a.coefficients),
            actual_normalized_moment_shape_axial5_coefficients=coefficient_lists(moment_shapes),
            actual_true_moment_derivatives_divided_by_current_prefactor=coefficient_lists(true_shapes),
            full_backward_kernel_axial5_coefficients=coefficient_lists(kernels),
            inherited_R110_amplitude_decay_axial5_coefficients=coefficient_lists(decays),
            actual_inherited_R110_moment_contribution_axial5_coefficients=coefficient_lists(inherited),
            moment_prefactors=dict(theta='sqrt(2)*R^(3/2)*Utheta',theta_z='sqrt(2)*R^(3/2)*Utheta',
                mean='R',axial='R',swirl='R*Utheta^2/2',pressure='Utheta^2/2'),
            exact_physical_moment_definitions=dict(Mtheta='sqrt(2)*R^(3/2)*Utheta*theta_shape',
                Mtheta_z='sqrt(2)*R^(3/2)*Utheta*theta_z_shape',Mz='R*mean',
                Mztheta='R*axial-R*Utheta^2*swirl_shape/2',Mp='Utheta^2*pressure_shape/2',
                P='Pstar^2*original_P0_normalized+Mp',Ur='sqrt(R/2)*Q'),
            source_full_kernel='integral_0^y exp(-k*t+m*B*(sigma((y-t)/T)-sigma(y/T)))dt; (k,m)=(1.6,1),(.2,2),(1.2,2)',
            source_inherited_decay='exp(-k*y+m*B*sigma(y/T)); original nonzero R110 histories',
            source_axial_velocity='V(R,Z)=actual_V110(Z)',
            actual_R110_histories_retained=True,original_pressure_datum_retained=True,
            local_current_amplitude_normalization=True,moment_reset=False,
            numerical_kernel_truncation=False,newly_recomputed_point_coefficients=False,
            reshape_radial_mixed4_certified=False,reference_continuation_installed=False,
            all_radial_endpoint_joins_certified=False)
        return result

    def report(self):
        whole=self.evaluate([-1,1],phase=[0,1]);start=self.evaluate([-1,1],phase=0)
        end=self.evaluate([-1,1],phase=1)
        samples=[self.evaluate(z,phase=p) for z in ('0','.5') for p in ('.25','.5','.75')]
        samples+=[self.evaluate(z,y=1) for z in ('0','.5')]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            admitted_inner_parameter_family_sha256=self.switch.bridge.records['K1_ledger']['admitted_inner_parameter_family_sha256'],
            reference_join_family_sha256=self.records['reference_join_bounds']['reference_join_family_sha256'],
            Abar=self.A,T=self.T,exact_T_source='400 times selected A_upper; no new fit or smaller length',
            log_Rref_over_110=self.logref,exact_reference_offsets=dict(Rz=-8,restore_end=-7,Rm=-6,Rh=-5),
            whole_reshape=whole,actual_R110_inlet=start,actual_Rsh_exit=end,samples=samples,
            actual_B_axial5_coefficients=list(self.inputs([-1,1])['B'].coefficients),
            B_C2_bound_source=PREFIX+'reference_join_bounds.json; SAME actual110 B, theorem norm<=2Abar; intersection only through2',
            kernel_derivative_proof='|delta sigma|<=8*t/T; coefficient Bell polynomials; full positive moments int_0^infinity t^p exp(-r*t)dt=p!/r^(p+1), smaller finite-y mass allowed',
            kernel_rate_bounds=dict(theta=['1.55','1.65'],pressure=['.1','.3'],swirl=['1.1','1.3']),
            direct_same_source_log_u_slope_interval=self.q_source,
            sigma_derivative_bound_source='lei_ren_part1_paper_shared_fixed_step_bound.json; certified upper8',
            physical_radius_or_amplitude_materialized=False,cap_is_not_source=True,
            factored_exponential_cap_proofs=self.proofs,
            actual_long_reshape_log_velocity_axial5_available=True,
            actual_long_reshape_moment_axial5_enclosures_available=True,
            actual_long_reshape_radial_recovery_axial4_available=True,
            actual_Rsh_exit_axial5_available=True,
            reshape_radial_mixed4_certified=False,reference_continuation_installed=False,
            actual_five_moment_patch_connected=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantLongReshapeProfiles().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual R110 long reshape log field and full inherited moment axial5 enclosures generated',flush=True)
    return result


if __name__=='__main__':run()

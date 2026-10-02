"""Actual Rsh reference continuation, axial restore and prepatch five defects.

All moment differences are centered algebraically BEFORE interval bounds.
Formal radii use exact relative offsets. The original axis datum is kept.
These are axial5 enclosures, not full radial mixed4 or a connected repair.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_long_reshape_profiles import (
    CompliantLongReshapeProfiles,IntervalTaylor,relative_amplitude_jet)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    logarithm,square,derivative,norm,symmetric,coefficient_lists)
from lei_ren_part1_paper_compliant_frozen_comparison_field import axial_jet
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets,positive_exp
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def scale_history(c,initial,rate,gap,logcap,proofs):
    """Enclose exp(-rate*gap)*initial, all factors combined before capping."""
    gap=c.mpf(gap);rate=c.mpf(rate)
    if endpoints(gap)[0]<0:raise ValueError('Reference transport gap must be nonnegative')
    if endpoints(gap)==(mp.mpf(0),mp.mpf(0)):return initial
    out=[]
    for k,value in enumerate(initial.coefficients):
        absolute=norm(c,value)
        if endpoints(absolute)[1]==0:out.append(c.mpf(0));continue
        logvalue=-rate*gap+c.ln(absolute)
        if endpoints(logvalue)[1]<=endpoints(logcap)[0]:
            bound=c.mpf([0,endpoints(c.exp(logcap))[1]])
            proofs.append(dict(rate=rate,source_gap_enclosure=gap,input_absolute_upper=absolute,
                log_magnitude_upper=logvalue,log_cap=logcap,source_not_replaced=True))
        else:bound=positive_exp(c,logvalue)
        upper=endpoints(bound)[1]
        if endpoints(value)[0]>=0:out.append(c.mpf([0,upper]))
        elif endpoints(value)[1]<=0:out.append(c.mpf([-upper,0]))
        else:out.append(c.mpf([-upper,upper]))
    return IntervalTaylor(c,out)


def restoration_kernels(c,t,cells=128):
    """Full finite restoration integrals with exact positive cell weights.

    J(k,j;t)=int_0^t exp(-k*(t-s))*(1-sigma(s))^j ds.
    Interval t retains the complete uncertain endpoint segment separately.
    No midpoint sigma or omitted quadrature tail defines the kernels.
    """
    t=c.mpf(t);lo,hi=endpoints(t)
    if lo<0 or hi>1 or not isinstance(cells,int) or cells<1:raise ValueError('Restoration t in[0,1], positive cells required')
    values={name:c.mpf(0) for name in ('mean','mixed','square')}
    definitions=(('mean',1,1),('mixed',c.mpf('1.6'),1),('square',1,2))
    if lo:
        for i in range(cells):
            left=c.mpf(lo)*i/cells;right=c.mpf(lo)*(i+1)/cells
            box=c.mpf([endpoints(left)[0],endpoints(right)[1]])
            alpha=intersection(c,1-sigma_jets(c,box)[0],c.mpf([0,1]))
            for name,rate,power in definitions:
                weight=c.exp(-rate*(c.mpf(lo)-right))*(1-c.exp(-rate*(right-left)))/rate
                weight=c.mpf([max(mp.mpf(0),endpoints(weight)[0]),endpoints(weight)[1]])
                values[name]+=weight*alpha**power
    length=c.mpf(hi)-c.mpf(lo)
    alpha_lo=intersection(c,1-sigma_jets(c,c.mpf(lo))[0],c.mpf([0,1]))
    for name,rate,power in definitions:
        prefix_factor=c.mpf([endpoints(c.exp(-rate*length))[0],1])
        extension=(1-c.exp(-rate*length))*alpha_lo**power/rate
        values[name]=values[name]*prefix_factor+c.mpf([0,max(mp.mpf(0),endpoints(extension)[1])])
    return values


def restore_centered(c,initial,E,t,kernels):
    """Algebraically centered actual histories under V=4Z+E*(1-sigma)."""
    t=c.mpf(t)
    return dict(mean_error=initial['mean_error']*c.exp(-t)+E*kernels['mean'],
        angular_error=initial['angular_error']*c.exp(-c.mpf('1.6')*t),
        mixed_error=initial['mixed_error']*c.exp(-c.mpf('1.6')*t)+E*kernels['mixed'],
        axial_square=initial['axial_square']*c.exp(-t)+square(E)*kernels['square'],
        swirl_error=initial['swirl_error']*c.exp(-c.mpf('1.2')*t),
        pressure_error=initial['pressure_error']*c.exp(-c.mpf('.2')*t))


class CompliantReferenceRestoreProfiles:
    def __init__(self):
        self.reshape=CompliantLongReshapeProfiles();self.core=self.reshape.core
        self.ctx=c=self.reshape.ctx;self.family=self.reshape.family;self.source=self.reshape.source
        self.hashes=dict(self.reshape.hashes);self.cache={};self.kernel_cache={};self.proofs=[]
        name=PREFIX+'long_reshape_profiles_check.json';check=json.loads((HERE/name).read_bytes())
        for path,digest in check['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Reference restore prerequisite changed: '+path)
        if not (check['all_passed'] and check['actual_Rsh_exit_axial5_available']
                and check['actual_five_defect_family_sha256']==self.family and check['implicit_source_sha256']==self.source
                and check['actual_B_C2_shear_sigma_and_kernel_rate_theorems_directly_bound']):
            raise ValueError('Actual source-bound Rsh history is required')
        self.hashes.update(check['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        admitted_name=PREFIX+'five_defect_admission.json'
        self.defect_admission=json.loads((HERE/admitted_name).read_bytes())
        admitted=self.defect_admission
        if not (hashlib.sha256((HERE/admitted_name).read_bytes()).hexdigest()==self.hashes[admitted_name]
                and admitted['actual_five_defect_family_sha256']==self.family
                and admitted['implicit_source_sha256']==self.source
                and admitted['reference_join_family_sha256']==self.reshape.records['reference_join_bounds']['reference_join_family_sha256']
                and admitted['directed_fixed_cutoff_integrals_not_midpoint_fit']
                and admitted['actual_defect_functions_defined_by_source_primitives']):
            raise ValueError('Original source restoration integral admission changed')
        self.loggap=self.reshape.logref-self.reshape.T
        if endpoints(self.loggap)[0]<=8:raise ValueError('Rsh<Rz topology missing')
        self.logcap=-1000*c.ln(10)-2*self.core.logP
        self.tailcap=c.exp(self.logcap) # Compact exponent; no exp(logCstar).
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def inputs(self,Z):
        c=self.ctx;Z=c.mpf(Z);key=Z._mpi_
        if key in self.cache:return self.cache[key]
        parent=self.reshape.evaluate(Z,phase=1);jet=lambda row:IntervalTaylor(c,row)
        moments={n:jet(row) for n,row in parent['actual_normalized_moment_shape_axial5_coefficients'].items()}
        z=IntervalTaylor.variable(c,Z,5)
        centered=dict(mean_error=moments['mean']-z*4,
            angular_error=moments['theta']-c.mpf(5)/8,
            mixed_error=moments['theta_z']-(z*moments['theta'])*4,
            axial_square=moments['axial']-(z*moments['mean'])*8+square(z)*16,
            swirl_error=moments['swirl']-c.mpf(5)/6,pressure_error=moments['pressure']-5)
        # The exact centered axial primitive is int_0^R (V-4Z)^2 dR/R.
        row=centered['axial_square'];lower,upper=endpoints(row[0])
        centered['axial_square']=IntervalTaylor(c,[c.mpf([max(mp.mpf(0),lower),upper])]+list(row.coefficients[1:]))
        # Obtain E=V110-4Z without subtracting two broad baseline intervals.
        # Vcore=4Z+j+epsilon*Psi; actual bridge and first switch add their
        # admitted three drive-component caps each. The actual V function
        # is unchanged through the second switch and the whole reshape.
        source=self.core.normalized_jets(4,Z)
        E=axial_jet(c,source['Psi'])*self.core.epsilon+self.core.j
        self.reshape.switch.inputs(Z) # Recheck all actual width/product operands.
        extra=self.reshape.switch.bridge.cap*6
        E=E+IntervalTaylor(c,[symmetric(c,extra)]*6)
        ledger=self.reshape.switch.bridge.records['K1_ledger'];gate=self.reshape.switch.bridge.records['global_exit_certificate']
        rho=read_interval(c,ledger['rho_core_C2_bound'])+read_interval(c,gate['rho_bridge_C2_upper'])
        coefficients=[]
        for k,q in enumerate(E.coefficients):
            if k<=2:
                bound=symmetric(c,rho)/math.factorial(k)+(self.core.j if k==0 else 0)
                q=intersection(c,q,bound)
            coefficients.append(q)
        E=IntervalTaylor(c,coefficients)
        result=dict(parent=parent,z=z,E=E,rho=rho,source_centered_Rsh=centered,
                    original_axis_pressure=jet(parent['pressure_axis_axial5_coefficients']))
        self.cache[key]=result;return result

    def reference_centered(self,Z,gap):
        c=self.ctx;inp=self.inputs(Z);gap=c.mpf(gap)
        if endpoints(gap)==(mp.mpf(0),mp.mpf(0)):return inp['source_centered_Rsh']
        scale=lambda row,rate:scale_history(c,row,rate,gap,self.logcap,self.proofs)
        E=inp['E'];initial=inp['source_centered_Rsh']
        # Exact cancellation of the reference forcing occurs before bounds.
        d=scale(IntervalTaylor.constant(c,1,5),'1.6')[0]
        return dict(mean_error=E+scale(initial['mean_error']-E,1),
            angular_error=scale(initial['angular_error'],'1.6'),
            mixed_error=scale(initial['mixed_error'],'1.6')+E*((1-d)/c.mpf('1.6')),
            axial_square=square(E)+scale(initial['axial_square']-square(E),1),
            swirl_error=scale(initial['swirl_error'],'1.2'),
            pressure_error=scale(initial['pressure_error'],'.2'))

    def packet(self,Z,centered,V,logu,region,**metadata):
        c=self.ctx;inp=self.inputs(Z);z=inp['z']
        H=centered['angular_error']+c.mpf(5)/8
        shapes=dict(theta=H,theta_z=(z*H)*4+centered['mixed_error'],
            mean=z*4+centered['mean_error'],
            axial=square(z)*16+(z*centered['mean_error'])*8+centered['axial_square'],
            swirl=centered['swirl_error']+c.mpf(5)/6,pressure=centered['pressure_error']+5)
        m=shapes['mean'];d=1-square(z);L=1-square(z)*self.core.delta
        Q=(2*z*V-(z*m)*(1-self.core.delta)-d*derivative(m))/L
        amp=relative_amplitude_jet(logu);amp2=relative_amplitude_jet(logu,2)
        true=dict(theta=shapes['theta']*amp,theta_z=shapes['theta_z']*amp,mean=m,axial=shapes['axial'],
                  swirl=shapes['swirl']*amp2,pressure=shapes['pressure']*amp2)
        return dict(Z=c.mpf(Z),region=region,log_Utheta_over_Pstar_axial5_coefficients=list(logu.coefficients),
            Utheta_true_axial5_divided_by_current_Utheta=list(amp.coefficients),
            Uz_actual_axial5_coefficients=list(V.coefficients),actual_Q_axial4_coefficients=list(Q.coefficients),
            pressure_axis_axial5_coefficients=list(inp['original_axis_pressure'].coefficients),
            actual_centered_moment_axial5_coefficients=coefficient_lists(centered),
            actual_normalized_moment_shape_axial5_coefficients=coefficient_lists(shapes),
            actual_true_moment_derivatives_divided_by_current_prefactor=coefficient_lists(true),
            actual_E_V110_minus_4Z_axial5_coefficients=list(inp['E'].coefficients),
            exact_physical_definitions=dict(Mtheta='sqrt2*R^1.5*Utheta*theta_shape',
                Mtheta_z='sqrt2*R^1.5*Utheta*theta_z_shape',Mz='R*mean',
                Mztheta='R*axial-R*Utheta^2*swirl_shape/2',Mp='Utheta^2*pressure_shape/2',
                P='Pstar^2*original_axis_pressure+Mp',Ur='sqrt(R/2)*Q'),
            source_centered_definitions=dict(mean_error='m-4Z',angular_error='H-5/8',mixed_error='K-4Z*H',
                axial_square='A-8Z*m+16Z^2',swirl_error='b-5/6',pressure_error='p-5'),
            actual_histories_retained=True,original_axis_pressure_retained=True,moment_reset=False,
            physical_Rref_materialized=False,radial_mixed4_certified=False,
            actual_five_moment_patch_connected=False,**metadata)

    def reference(self,Z,phase):
        """Entire Rsh..Rz, phase*the exact formal log length."""
        c=self.ctx;phase=c.mpf(phase)
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Reference phase in[0,1] required')
        inp=self.inputs(Z);gap=(self.loggap-8)*phase
        centered=self.reference_centered(Z,gap);logq=logarithm(1+square(inp['z']))
        if endpoints(phase)==(mp.mpf(1),mp.mpf(1)):
            logu=-logq-c.mpf('.8')
        else:
            logu=-logq+(self.reshape.T/10-self.core.logC-self.core.logP+gap/10)
            # Monotone reference swirl cannot exceed its known Rz value.
            lo,hi=endpoints(logu[0]);upper=endpoints((-logq-c.mpf('.8'))[0])[1]
            logu=IntervalTaylor(c,[c.mpf([lo,min(hi,upper)])]+list(logu.coefficients[1:]))
        return self.packet(Z,centered,inp['z']*4+inp['E'],logu,'reference_before_axial_restore',
            phase=phase,source_log_gap_from_Rsh=gap,
            exact_source_radius='Rsh*exp(phase*(10*(logCstar+logPstar)-T-8))',
            exact_Rz_log_offset=-8 if endpoints(phase)==(mp.mpf(1),mp.mpf(1)) else None,
            axial_source='V=actualV110=4Z+E')

    def kernels(self,t):
        t=self.ctx.mpf(t);key=t._mpi_
        if key not in self.kernel_cache:
            values=restoration_kernels(self.ctx,t)
            if endpoints(t)==(mp.mpf(1),mp.mpf(1)):
                # Independent, previously admitted full scalar integrals of
                # the SAME original cutoff sharpen the identical endpoint
                # functions. Their positive dominant kernels do NOT erase
                # any negative inherited term in the actual histories.
                for name,part,rate in (('mean','linear_1',1),('mixed','linear_8_over_5','1.6'),('square','square_1',1)):
                    original=read_interval(self.ctx,self.defect_admission['directed_restoration_integrals'][part])
                    values[name]=intersection(self.ctx,values[name],original*self.ctx.exp(-self.ctx.mpf(rate)))
            self.kernel_cache[key]=values
        return self.kernel_cache[key]

    def restoration(self,Z,t):
        c=self.ctx;t=c.mpf(t)
        if endpoints(t)[0]<0 or endpoints(t)[1]>1:raise ValueError('Axial restoration t in[0,1] required')
        inp=self.inputs(Z);initial=self.reference_centered(Z,self.loggap-8);kernels=self.kernels(t)
        centered=restore_centered(c,initial,inp['E'],t,kernels)
        alpha=intersection(c,1-sigma_jets(c,t)[0],c.mpf([0,1]))
        V=inp['z']*4+inp['E']*alpha
        logu=-logarithm(1+square(inp['z']))+(-8+t)/10
        return self.packet(Z,centered,V,logu,'axial_restoration',restoration_t=t,
            exact_source_radius='Rz*exp(t), Rz=exp(-8)*Rref',exact_reference_log_offset=-8+t,
            original_restoration_kernels=kernels,axial_source='V=4Z+E*(1-sigma(t)), E=actualV110-4Z')

    def terminal(self,Z,offset):
        """Unpatched terminal reference-velocity branch, offsets[-7,-5]."""
        c=self.ctx;offset=c.mpf(offset)
        if endpoints(offset)[0]<-7 or endpoints(offset)[1]>-5:raise ValueError('Terminal offset in[-7,-5] required')
        inp=self.inputs(Z);initial=self.reference_centered(Z,self.loggap-8)
        restored=restore_centered(c,initial,inp['E'],c.mpf(1),self.kernels(1));length=offset+7
        rates=dict(mean_error=1,angular_error='1.6',mixed_error='1.6',axial_square=1,swirl_error='1.2',pressure_error='.2')
        centered={name:row*c.exp(-c.mpf(rates[name])*length) for name,row in restored.items()}
        logu=-logarithm(1+square(inp['z']))+offset/10
        return self.packet(Z,centered,inp['z']*4,logu,'unpatched_terminal_reference_velocity',
            exact_reference_log_offset=offset,exact_source_radius='exp(offset)*Rref',
            axial_source='V=4Z exactly; all actual moment differences remain',
            unpatched_profile=True)

    def defects(self,Z,offset=-5):
        c=self.ctx;packet=self.terminal(Z,offset);offset=c.mpf(offset);x=offset+6
        centered={n:IntervalTaylor(c,row) for n,row in packet['actual_centered_moment_axial5_coefficients'].items()}
        z=self.inputs(Z)['z'];invAm2=square(1+square(z))*c.exp(c.mpf('1.2')-2*self.core.logP)
        rows=[centered['mean_error']*c.exp(x),
            centered['mixed_error']*c.exp(c.mpf('1.6')*x),
            centered['angular_error']*c.exp(c.mpf('1.6')*x),
            (centered['axial_square']*invAm2)*c.exp(x)-centered['swirl_error']*(c.exp(c.mpf('1.2')*x)/2),
            centered['pressure_error']*(c.exp(c.mpf('.2')*x)/2)]
        return dict(Z=c.mpf(Z),exact_reference_offset=offset,
            actual_normalized_five_defect_axial5_coefficients=[list(row.coefficients) for row in rows],
            source_row_order=['Delta_z/Rm','(Delta_theta_z-4Z*Delta_theta)/(sqrt2*Rm^1.5*Am)',
                'Delta_theta/(sqrt2*Rm^1.5*Am)','(Delta_ztheta-8Z*Delta_z)/(Rm*Am^2)','Delta_p/Am^2'],
            exact_scale_source='Rm=exp(-6)*Rref, Am=exp(-.6)*Pstar/(1+Z^2)',
            baseline_reference_terms_cancelled_before_enclosure=True,
            source_moments=packet,actual_five_moment_patch_connected=False)

    def report(self):
        whole_reference=self.reference([-1,1],[0,1]);Rsh=self.reference([-1,1],0);Rz=self.reference([-1,1],1)
        whole_restore=self.restoration([-1,1],[0,1]);start=self.restoration([-1,1],0);end=self.restoration([-1,1],1)
        Rm=self.terminal([-1,1],-6);Rh=self.terminal([-1,1],-5)
        defects=self.defects([-1,1]);samples=[self.restoration(z,'.5') for z in ('0','.5')]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            reference_join_family_sha256=self.reshape.records['reference_join_bounds']['reference_join_family_sha256'],
            admitted_inner_parameter_family_sha256=self.reshape.switch.bridge.records['K1_ledger']['admitted_inner_parameter_family_sha256'],
            whole_reference=whole_reference,actual_Rsh_inlet=Rsh,actual_Rz_inlet=Rz,
            whole_axial_restore=whole_restore,restore_start=start,restore_end=end,unpatched_Rm=Rm,unpatched_Rh=Rh,
            actual_five_defects_at_Rh=defects,samples=samples,
            source_log_Rref_over_Rsh=self.loggap,exact_reference_offsets=dict(Rz=-8,restore_end=-7,Rm=-6,Rh=-5),
            centered_axial_source_proof=dict(core='Ecore=j+epsilon*Psi(4,Z), analytic fixed point',
                bridge_switch='Same actual source drive logs admit three caps per bridge/switch; total extra coefficient bound6cap',
                cap_coefficient=self.reshape.switch.bridge.cap*6,
                C2_error_about_j_upper=self.inputs([-1,1])['rho'],
                source_ledger=PREFIX+'K1_ledger.json',source_gate=PREFIX+'global_exit_certificate.json'),
            retained_source_transport_cap_proofs=self.proofs,positive_log_cap=self.logcap,positive_tail_cap=self.tailcap,
            restoration_integrals='int_0^t exp(-k*(t-s))*(1-sigma(s))^j ds; k/j=1/1,1.6/1,1/2',
            restoration_cells=128,original_sigma_and_complete_endpoint_segment_retained=True,
            identical_restoration_endpoint_refined_by_admitted_full_integrals=True,
            actual_reference_continuation_axial5_enclosures_available=True,
            actual_axial_restoration_axial5_enclosures_available=True,
            actual_unpatched_Rm_Rh_moment_axial5_enclosures_available=True,
            actual_five_defect_profile_axial5_enclosures_available=True,
            actual_radial_recovery_axial4_available=True,
            radial_mixed4_certified=False,all_radial_endpoint_joins_certified=False,
            actual_five_moment_patch_connected=False,newly_recomputed_point_coefficients=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantReferenceRestoreProfiles().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rsh reference/axial restore toRh and centered five-defect axial5 enclosures generated',flush=True)
    return result


if __name__=='__main__':run()

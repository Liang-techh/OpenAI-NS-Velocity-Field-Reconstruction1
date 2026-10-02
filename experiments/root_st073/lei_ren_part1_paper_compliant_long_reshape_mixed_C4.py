"""Actual full long reshape: physical logR/Z derivatives through four.

The original selected T and positive source amplitudes remain formal.
High radial rows come from physical primitive RHSs, not numerical kernels.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_long_reshape_profiles import (
    CompliantLongReshapeProfiles,IntervalTaylor,relative_amplitude_jet,bound_exp)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square,derivative
from lei_ren_part1_paper_compliant_inner_switch_profiles import MTH,MTHZ,MZ,MZT,MP
from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets,positive_exp
from lei_ren_part1_paper_compliant_reference_restore_mixed_C4 import binomial_rate
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def exponential_derivatives(log_derivatives):
    """Ordinary derivatives exp(L)/exp(L0), with L'(y)..L''''(y)."""
    one=log_derivatives[0]*0+1;out=[one]
    for n in range(4):
        out.append(sum((out[n-j]*log_derivatives[j]*math.comb(n,j) for j in range(n+1)),one*0))
    return out


def scaled_positive_source(c,logbase,factor,proofs):
    """Enclose exp(logbase)*factor AFTER combining large derivative factors.

    The cap is an enclosure, never a replacement for the exact positive
    source. Moderate finite fixtures retain direct interval multiplication.
    """
    if endpoints(logbase)[1]>-1000:return factor*positive_exp(c,logbase)
    values=[]
    for value in factor.coefficients:
        lo,hi=endpoints(value);size=max(abs(lo),abs(hi))
        if size==0:values.append(c.mpf(0));continue
        magnitude=bound_exp(c,logbase+c.ln(c.mpf(size)),proofs);upper=endpoints(magnitude)[1]
        values.append(c.mpf([0,upper] if lo>=0 else [-upper,0] if hi<=0 else [-upper,upper]))
    return IntervalTaylor(c,values)


def reshape_mixed(c,Z,delta,logu,log_y,V,shapes,p0,invP2,proofs=None):
    """Physical rows at a fixed current basepoint, ordinary y/Z units."""
    if proofs is None:proofs=[]
    if len(log_y)!=4 or any(row.order!=5 for row in log_y):
        raise ValueError('Original log-u y derivatives1..4 with axial5 required')
    z=IntervalTaylor.variable(c,c.mpf(Z),5);zero=z*0
    amp=relative_amplitude_jet(logu);amp2=relative_amplitude_jet(logu,2)
    shifted=lambda rate,m:exponential_derivatives([log_y[0]*m+rate]+[row*m for row in log_y[1:]])
    urows=shifted(0,1);theta_rhs=shifted(c.mpf('1.5'),1)
    swirl_rhs=shifted(1,2);pressure_rhs=shifted(0,2)
    mean=[shapes['mean']]+[(V-shapes['mean'])*((-1)**(k-1)) for k in range(1,5)]
    velocity=[V]+[zero]*4
    Q=[(2*z*velocity[k]-(z*mean[k])*(1-delta)-(1-square(z))*derivative(mean[k]))/(1-square(z)*delta) for k in range(5)]
    scaled=lambda row:scaled_positive_source(c,2*logu[0],amp2*row,proofs)
    pressure=[scaled(shapes['pressure'])/2]+[scaled(pressure_rhs[k-1])/2 for k in range(1,5)]
    physical=dict(Utheta_over_current_Utheta=[amp*row for row in urows],Uz=velocity,
        Ur_over_current_sqrt_R_over_2=[binomial_rate(Q,c.mpf('.5'),k) for k in range(5)],
        P_over_Pstar2=[p0+pressure[0]]+pressure[1:])
    primitives=dict(Mtheta_over_current_sqrt2_R_1p5_Utheta=[amp*shapes['theta']]+[amp*theta_rhs[k-1] for k in range(1,5)],
        Mtheta_z_over_current_sqrt2_R_1p5_Utheta=[amp*shapes['theta_z']]+[(amp*theta_rhs[k-1])*V for k in range(1,5)],
        Mz_over_current_R=[shapes['mean']]+[V]*4,
        Mztheta_over_current_R_Pstar2=[shapes['axial']*invP2-scaled(shapes['swirl'])/2]
            +[square(V)*invP2-scaled(swirl_rhs[k-1])/2 for k in range(1,5)],Mp_over_Pstar2=pressure)
    grid=lambda rows:{'y'+str(k)+'_Z'+str(n):row[n]*math.factorial(n) for k,row in enumerate(rows) for n in range(5-k)}
    return dict(physical_velocity_pressure_y_Z_mixed4={name:grid(rows) for name,rows in physical.items()},
        physical_five_primitive_y_Z_mixed4={name:grid(rows) for name,rows in primitives.items()},
        actual_log_Utheta_y_derivative_axial5=log_y,actual_Q_y_derivative_axial4=Q,
        exact_positive_swirl_source_log=logu[0],positive_source_cap_is_enclosure_only=True,
        large_derivative_factors_combined_before_positive_source_cap=True,
        physical_radial_prefactors_differentiated_before_mixed_grid=True)


class CompliantLongReshapeMixedC4:
    def __init__(self):
        self.reshape=CompliantLongReshapeProfiles();self.ctx=c=self.reshape.ctx
        self.family=self.reshape.family;self.source=self.reshape.source;self.hashes=dict(self.reshape.hashes)
        self.proofs=[]
        for part in ('long_reshape_profiles_check','reference_restore_mixed_C4_check'):
            name=PREFIX+part+'.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or receipt['actual_five_defect_family_sha256']!=self.family or receipt['implicit_source_sha256']!=self.source:
                raise ValueError('Same actual reshape/reference family and implicit source required')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Reshape mixed prerequisite changed: '+path)
            self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.invP2=c.exp(-2*self.reshape.core.logP)
        self.shared_axial_source=self.formal_axial_source()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def formal_axial_source(self):
        """Retain signed actual integrals; caps never define these functions.

        V and its centered E share the SAME two function references at all
        axial orders. Numerical directed enclosures remain separate from
        this exact source graph; no arbitrary cap representative is chosen.
        """
        switch=self.reshape.switch;bridge=switch.bridge
        actual_bridge=bridge.actual([-1,1],bridge.r/100)
        actual_switch=switch.phase([-1,1],1)
        sources=actual_bridge['source_integral_definitions']
        terms=dict(I_bridge=dict(variable='s=log(R/Ra)',bounds=['0','log(100/Ra)'],sign=-1,
            chi='1-(1-hb)*sigma(s/hb)',quotient='phi_actual(R,Z)/phi_bar(R,Z)',
            radius='R=Ra*exp(s)',drive_terms=['R*hydro','R*Pstar^2*pressure','R^2*F0^2*swirl'],
            exact_original_V_source=sources['V'],exact_original_chi_source=sources['chi']),
            I_first_switch=dict(variable='t',bounds=['0','1'],sign=-1,factor='hb^2',
            weight='1-sigma(t)',quotient='phi_actual(R,Z)/phi_bar(R,Z)',radius='R=100*exp(hb*t)',
            drive='same current-radius axial direction sum; not the angular Dbar',
            drive_terms=['R*hydro','R*Pstar^2*pressure','R^2*F0^2*swirl'],
            exact_original_V_source=actual_switch['exact_Uz_source']))
        bindings={name:self.hashes[name] for name in (PREFIX+'core_physical_field.py',
            PREFIX+'inner_bridge_profiles.py',PREFIX+'inner_switch_profiles.py',
            PREFIX+'reference_restore_profiles.py',PREFIX+'flat_pulse_derivatives.py')}
        namespace=hashlib.sha256(json.dumps(dict(family=self.family,implicit_source=self.source,
            integrals=terms,source_bindings=bindings),sort_keys=True).encode('utf8')).hexdigest()
        Zterm=dict(op='mul',args=[4,'Z']);jterm=dict(source='same selected j')
        psiterm=dict(op='mul',args=['same epsilon_core','same analytic Psi(rho=4,Z)'])
        bridge_ref=dict(formal_integral='I_bridge',shared_source=namespace)
        switch_ref=dict(formal_integral='I_first_switch',shared_source=namespace)
        v100=dict(op='sum',args=[Zterm,jterm,psiterm,bridge_ref])
        v110=dict(op='sum',args=v100['args']+[switch_ref])
        centered=dict(op='sum',args=v110['args'][1:])
        return dict(shared_source_namespace=namespace,formal_signed_integrals=terms,source_bindings=bindings,
            core_baseline='4Z+j+epsilon_core*Psi(4,Z)',V100=v100,V110=v110,E_V110_minus_4Z=centered,
            axial_Taylor_orders=list(range(6)),ordinary_derivative_factorials=[math.factorial(k) for k in range(6)],
            V110_enclosure_binding='actual switch.v_cover inherited unchanged through second switch, post power and Rsh',
            E_enclosure_binding='same V110 minus4Z; directed correlated core+six-drive bound only',
            caps_are_not_source_integrals=True,formal_integrals_numerically_reconstructed=False)

    def evaluate(self,Z,phase):
        c=self.ctx;parent=self.reshape.evaluate(Z,phase=phase);inp=self.reshape.inputs(Z)
        phase=c.mpf(phase);cutoff=sigma_jets(c,phase);B=inp['B'];T=self.reshape.T
        log_y=[B*(-cutoff[k]*math.factorial(k)/T**k)+(c.mpf('.1') if k==1 else 0) for k in range(1,5)]
        logu=IntervalTaylor(c,parent['log_Utheta_over_Pstar_axial5_coefficients'])
        shapes={name:IntervalTaylor(c,row) for name,row in parent['actual_normalized_moment_shape_axial5_coefficients'].items()}
        packet=reshape_mixed(c,Z,self.reshape.core.delta,logu,log_y,inp['v'],shapes,inp['p0'],self.invP2,self.proofs)
        return dict(Z=c.mpf(Z),phase=phase,**packet,actual_inherited_axial5_packet=parent,
            actual_axial_function_source=self.shared_axial_source['V110'],
            reference_centered_axial_function_source=self.shared_axial_source['E_V110_minus_4Z'],
            original_T=T,source_T='400*selected_Abar; all y derivatives keep original inverse-T factors',
            actual_R110_histories_and_P0_retained=True,full_backward_kernels_not_truncated=True,
            derivative_normalization='fixed current-basepoint physical prefactors; y=logR independent of Z',
            reshape_radial_mixed4_certified=True,Rsh_reference_mixed4_join_certified=True,
            R110_postswitch_power_reshape_local_mixed4_join_certified=True,
            microscopic_switch_mixed4_certified=False,full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False)

    def report(self):
        whole=self.evaluate([-1,1],[0,1]);start=self.evaluate([-1,1],0);end=self.evaluate([-1,1],1)
        samples=[self.evaluate(z,phase) for z,phase in (('0','.25'),('.5','.5'),('0','.75'))]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            whole_reshape=whole,actual_R110_inlet=start,actual_Rsh_exit=end,interior_packets=samples,
            whole_local_R109_R110_power=self.power_before_R110([-1,1],[109,110]),
            actual_R110_power_side=self.power_before_R110([-1,1],110),
            shared_exact_axial_source=self.shared_axial_source,
            factored_positive_source_cap_proofs=self.proofs,
            actual_long_reshape_all_mixed_derivatives_total_order_le4_available=True,
            Rsh_reference_mixed4_join_certified=True,R110_postswitch_power_reshape_local_mixed4_join_certified=True,
            microscopic_switch_mixed4_certified=False,full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False,
            input_hashes=self.hashes)

    def power_before_R110(self,Z,R):
        """Actual two-sided R110 neighborhood; excludes microscopic switches."""
        c=self.ctx;R=c.mpf(R)
        if endpoints(R)[0]<109 or endpoints(R)[1]>110:raise ValueError('Local original power neighborhood R109..R110 required')
        parent=self.reshape.switch.post(Z,R);inp=self.reshape.inputs(Z)
        phi=IntervalTaylor(c,parent['F_actual_over_F0_axial5_coefficients'])
        raw=parent['actual_moment_shape_axial5_coefficients']
        jet=lambda row:IntervalTaylor(c,row)
        shapes=dict(theta=jet(raw[MTH])/(phi*2),theta_z=jet(raw[MTHZ])/(phi*2),
            mean=jet(raw[MZ]),axial=jet(raw[MZT]['axial']),
            swirl=jet(raw[MZT]['swirl'])/square(phi),pressure=jet(raw[MP])/square(phi))
        z=IntervalTaylor.variable(c,c.mpf(Z),5)
        logu=inp['B']-logarithm(1+square(z))-self.reshape.core.logC-self.reshape.core.logP+c.ln(R/110)/10
        log_y=[IntervalTaylor.constant(c,'.1',5)]+[IntervalTaylor.constant(c,0,5)]*3
        result=reshape_mixed(c,Z,self.reshape.core.delta,logu,log_y,
            jet(parent['Uz_actual_axial5_coefficients']),shapes,inp['p0'],self.invP2,self.proofs)
        return dict(Z=c.mpf(Z),R=R,**result,actual_inherited_axial5_packet=parent,
            scope='original postswitch power R109..R110 only; microscopic switches excluded',
            original_power_and_R2_histories_retained=True,microscopic_switch_mixed4_certified=False,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False)


def run():
    with mp.workdps(280):result=CompliantLongReshapeMixedC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual original long reshape physical mixed4 generated',flush=True)
    return result


if __name__=='__main__':run()

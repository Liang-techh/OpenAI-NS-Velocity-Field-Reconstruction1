"""Original own mean and radial recovery with epsilon/Pstar sectors exposed.

The actual mean is recovered by its own FTC history from the genuine core
inlet. Collar increments receive analytic bounds with formal hb factors.
No subtraction of independently boxed raw V and M is used to expose epsilon.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_axial_epsilon_factor as axial
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=axial.op;domain=axial.domain;weighted=axial.weighted;core=axial.core;ep=axial.ep
HERE,PREFIX,sha=axial.HERE,axial.PREFIX,axial.sha
NAME=PREFIX+'current_original_n1_leading_mean_epsilon_factor.json.gz'
RECEIPT=PREFIX+'current_original_n1_leading_mean_epsilon_factor_check.json'
GATE='current_original_leading_own_mean_and_radial_Q_epsilon_factors_installed'
OPEN=axial.OPEN
ORDERS=(5,4,3)


def core_rows(algebra,physical,rho,Z):
    """Original Psi and its own radial mean, through total mixed order five."""
    c=algebra.ctx;axis=physical.axis_inputs(Z);z,u,d,L,H,p=(axis[k] for k in ('z','u','d','L','H','pressure'))
    hydro=((1-2*z*u)*u*((1+physical.delta)/2)+4*H)/(2*L)
    pressure=(d*core.plain_dz(p)-2*z*p*(1+physical.delta))/(2*L)
    slope=algebra.lift(hydro.truncate(5))+algebra.shift(algebra.lift(pressure.truncate(5)),(0,0,1,0))
    raw=physical.normalized_jets(rho,Z);key=domain.parent.parent.core_module.gridkey
    psi=[];mean=[];phi=[];axis_zero=ep(c.mpf(rho))[1]==0
    for i,n in enumerate(ORDERS):
        coefficients=[raw['nonlinear_correction_bounds'][key(i,k)]/math.factorial(k) for k in range(n+1)]
        err=algebra.lift(op.Jet(c,[c.mpf([-1,1])*v for v in coefficients]))
        av=err/(2 if i==0 else i+1)
        model=op.truncate(slope,n)*(rho if i==0 else 1 if i==1 else 0)
        psi.append(model+err);mean.append(model/2+av)
        phi.append(algebra.lift(op.Jet(c,[raw['Phi'][key(i,k)]/math.factorial(k) for k in range(n+1)])))
        if axis_zero and i==0:psi[i]=mean[i]=algebra.lift(0,n)
        if axis_zero and i==1:psi[i]=model;mean[i]=model/2
    return dict(Psi=psi,meanPsi=mean,Phi=phi,slope=slope,axis=axis)


def cauchy_signed_rows(bound,order,gap):
    """Ordinary derivative k! cancels the Taylor coefficient k! divisor."""
    c=bound.ctx
    return op.FactoredJet(bound.algebra,{p:op.Jet(c,[c.mpf([-1,1])*domain.upper(c,row[0]/gap**k)
                                                   for k in range(order+1)])
                                       for p,row in bound.terms.items()},order)


def complex_mean_bounds(algebra,physical,b,chi_tube):
    """Psi/meanPsi on eta/2, Qpsi on eta/4; original radial width powers."""
    c=algebra.ctx;eta=b['common_tubes']['outer_Z_radius'];z=1+eta;d=1+z*z
    L=b['denominator_bounds']['L_modulus_lower'];dt=abs(physical.delta)
    r=b['actual_inner_geometry']['real_rho_upper_bound'];cap=b['actual_inner_geometry']['numerical_cap_bounds_only']
    u=4*z+abs(physical.j);H=z*(1+dt)/2+d*u
    p=b['original_pressure']['normalized_modulus_upper'];pz=2*p/eta
    hydro=domain.upper(c,((1+2*z*u)*u*(1+dt)/2+4*H)/(2*L))
    pressure=domain.upper(c,(d*pz+2*z*p*(1+dt))/(2*L))
    lift=lambda value:algebra.lift(value,0)
    slope=lift(hydro)+algebra.shift(lift(pressure),(0,0,1,0))
    correction=domain.upper(c,4*physical.correction)
    pc1=b['actual_core_radial_bounds']['Phi'][1]['0']['coefficient']
    floor=b['denominator_bounds']['Phi_core_modulus_lower']
    phi=b['actual_collar_radial_bounds']['Phi'][0]['0']['coefficient']
    ratio=domain.upper(c,phi/floor);A=domain.upper(c,pc1/floor)
    sig=domain.sigma_jets(c,c.mpf([0,'.5']))
    sig1=domain.upper(c,abs(sig[1]));chir=algebra.width(lift(sig1/4),-1)
    core_psi=[slope*r+lift(correction),slope+lift(2*correction),lift(8*correction)]
    # Original chiSpatial is real and in[0,1] for the real radial collar;
    # it has no Z dependence. A 1+cap bound is additionally conservative.
    psi0=core_psi[0]+algebra.width(core_psi[1]*(ratio*(1+cap)*r/2))
    psi1=core_psi[1]*(ratio*(1+cap))
    psi2=(chir*core_psi[1]+core_psi[1]*(A*(1+cap)**2)+core_psi[2]*(1+cap))*ratio
    inlet_mean=slope*2+lift(correction)
    mean0=inlet_mean+algebra.width(psi0*(r/8))
    mean1=(psi0+mean0)/4;mean2=(psi1+2*mean1)/4
    psis=[psi0,psi1,psi2];means=[mean0,mean1,mean2]
    # Mean derivatives in Q are bounded by Cauchy eta/4 from eta/2.
    qpsi=[(P*(2*z)+M*((1+dt)*z)+M*(d/(eta/4)))/L for P,M in zip(psis,means)]
    uaxis=lift(u);qa=lift(((1+dt)*z*u+4*d)/L)
    eps=lambda v:algebra.shift(v,(0,0,0,-1))
    w=[uaxis+eps(psis[0]),eps(psis[1]),eps(psis[2])]
    mean=[uaxis+eps(means[0]),eps(means[1]),eps(means[2])]
    q=[qa+eps(qpsi[0]),eps(qpsi[1]),eps(qpsi[2])]
    # Point FTC increment bounds use the same original actual Phi ratio.
    psi_increment_per_hb=core_psi[1]*(ratio*(1+cap)*r)
    mean_increment_per_hb=psi0*(r/4)
    proof=dict(original_outer_eta=eta,Psi_and_meanPsi_bound_tube_radius=eta/2,
        Qpsi_and_Q_bound_tube_radius=eta/4,Q_mean_Z_Cauchy_gap=eta/4,
        point_increment_Taylor_Cauchy_radius=eta/4,actual_rho_lower=4,
        actual_mean_inlet_multiplier='4/rho<=1',actual_mean_integral_measure='drho=hb*rho*ds; phase<=1/2',
        actual_mean_integral_upper_factor='hb*rho_max/8',core_mean_nonlinear_value_factor=4,
        core_Psi_rho1_nonlinear_factor=8,core_Psi_rho2_nonlinear_factor=32,
        actual_Phi_over_core_modulus_upper=ratio,actual_core_logarithmic_derivative_modulus_upper=A,
        actual_chi_rho_inverse_hb_retained=True,meanPsi_rho2_has_no_inverse_hb=True,
        norm_bounds_not_selected_comparison_field=True,all_original_pressure_atoms_retained=True)
    return dict(Psi=psis,meanPsi=means,Qpsi=qpsi,Uz=w,mean=mean,Q=q,
                Psi_increment_per_hb=psi_increment_per_hb,mean_increment_per_hb=mean_increment_per_hb),proof


def actual_rows(algebra,physical,raw_rows,rho,phase,Z,uniform,eta):
    c=algebra.ctx;inlet=core_rows(algebra,physical,c.mpf(4),Z)
    extension=core_rows(algebra,physical,rho,Z)
    if ep(phase)[1]==0:
        return inlet
    # Each collar increment uses its own real integration history; the
    # entire analytic coefficient of each source sector has a Cauchy bound.
    psi_error=algebra.width(uniform['Psi_increment_per_hb']*phase)
    mean_error=algebra.width(uniform['mean_increment_per_hb']*phase)
    psi0=inlet['Psi'][0]+cauchy_signed_rows(psi_error,5,eta/4)
    mean0=inlet['meanPsi'][0]*(4/rho)+cauchy_signed_rows(mean_error,5,eta/4)
    phi0=raw_rows['Phi'][0];key=domain.parent.parent.core_module.gridkey
    raw_core=physical.normalized_jets(rho,Z)
    pc=[op.Jet(c,[raw_core['Phi'][key(i,k)]/math.factorial(k) for k in range(n+1)]) for i,n in enumerate(ORDERS)]
    ratio=phi0/pc[0];A=pc[1]/pc[0];Ar=pc[2]/pc[0]-A*A
    sigma=domain.sigma_jets(c,phase)
    chi=algebra.lift(1-sigma[0],4)+algebra.width(algebra.lift(sigma[0],4))
    chir=algebra.lift(sigma[1]/rho,3)-algebra.width(algebra.lift(sigma[1]/rho,3),-1)
    psi1=op.truncate(ratio*chi*extension['Psi'][1],4)
    psi2=op.truncate(ratio*(chir*extension['Psi'][1]+chi*(chi-1)*A*extension['Psi'][1]+chi*extension['Psi'][2]),3)
    mean1=op.truncate((psi0-mean0)/rho,4)
    mean2=op.truncate((psi1-2*mean1)/rho,3)
    phi=[phi0,op.truncate(phi0*chi*A,4),op.truncate(phi0*(chir*A+chi*chi*A*A+chi*Ar),3)]
    return dict(Psi=[psi0,psi1,psi2],meanPsi=[mean0,mean1,mean2],Phi=phi,
                slope=inlet['slope'],axis=inlet['axis'])


def factored_leading_rows(algebra,rows,physical):
    axis=rows['axis'];z=algebra.lift(axis['z']);u=algebra.lift(axis['u']);d=algebra.lift(axis['d']);L=algebra.lift(axis['L'])
    qa=(u*(z*(1+physical.delta))-d*op.dz(u))/L
    eps=lambda v:algebra.shift(v,(0,0,0,-1))
    psi,meanpsi=rows['Psi'],rows['meanPsi']
    qpsi=[(2*z*P-z*M*(1-physical.delta)-d*op.dz(M))/L for P,M in zip(psi,meanpsi)]
    w=[op.truncate(u,n)+eps(P) if i==0 else eps(P) for i,(P,n) in enumerate(zip(psi,ORDERS))]
    mean=[op.truncate(u,n)+eps(M) if i==0 else eps(M) for i,(M,n) in enumerate(zip(meanpsi,ORDERS))]
    q=[op.truncate(qa,n-1)+eps(op.truncate(Q,n-1)) if i==0 else eps(op.truncate(Q,n-1)) for i,(Q,n) in enumerate(zip(qpsi,ORDERS))]
    return dict(Uz=w,mean=mean,Q=q,Qaxis=qa,Qpsi=qpsi)


@dataclass(frozen=True,eq=False)
class OriginalLeadingMeanEpsilonPacket:
    leading_hierarchy_order:int=0


class CurrentOriginalN1LeadingMeanEpsilonFactor:
    @source_precision
    def __init__(self,accepted_axial,require_checked=True):
        if type(accepted_axial) is not axial.CurrentOriginalN1AxialEpsilonFactor or not accepted_axial.acceptance_loaded:
            raise ValueError('Accepted original axial epsilon field required')
        accepted_axial.assert_graph();self.axial=accepted_axial;self.inputs=accepted_axial.inputs
        self.local=accepted_axial.local;self.weighted=accepted_axial.weighted;self.physical=accepted_axial.physical;self.c=accepted_axial.c
        self.bound=accepted_axial.bound;self.family=copy.deepcopy(accepted_axial.family);self.hashes=dict(accepted_axial.hashes)
        for name in (axial.NAME,axial.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._axial=accepted_axial;self._identities=(self.inputs,self.local,self.weighted,self.physical,self.c)
        self._bound=op.snapshot(self.bound);self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Leading own-mean factor receipt differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Leading own-mean factor source changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.axial is not self._axial or (self.inputs,self.local,self.weighted,self.physical,self.c)!=self._identities
            or self.inputs is not self.axial.inputs or self.local is not self.axial.local or self.weighted is not self.axial.weighted
            or self.physical is not self.axial.physical or self.c is not self.axial.c or not self.axial.acceptance_loaded
            or self.family!=self.axial.family or self.hashes!=self._hashes or op.snapshot(self.bound)!=self._bound):
            raise ValueError('Leading own-mean factor graph changed')
        return self.axial.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate,chart='first'):
        self.assert_graph();began=time.monotonic();c=self.c;value=c.mpf(coordinate)
        if not isinstance(coordinate,(str,int,float)) and ep(value)[0]!=ep(value)[1]:raise ValueError('One exact selected endpoint required')
        lp=self.local.evaluate(Z,coordinate,chart);lv=self.local.report(lp);query=lv['query'];z=lv['Z']
        p=self.inputs.evaluate(z,query['actual_source_coordinate'],chart);sv=self.inputs.report(p);entry=self.inputs._packets[id(p)]
        algebra=entry[5];rho=sv['geometry'].get('rho',sv['geometry'].get('rho_enclosure_only'))
        uniform,proof=complex_mean_bounds(algebra,self.physical,self.bound,self.weighted.tube)
        rows=(core_rows(algebra,self.physical,rho,z) if chart=='core' else
              actual_rows(algebra,self.physical,entry[2],rho,query['actual_source_coordinate'],z,uniform,proof['original_outer_eta']))
        fields=factored_leading_rows(algebra,rows,self.physical)
        gradient=op.Jet(c,self.weighted.source.amplitude.gradient_jets(z))
        frows=[algebra.shift(phi*weighted.relative_amplitude(algebra,gradient,1,n),(0,1,0,0)) for phi,n in zip(rows['Phi'],ORDERS)]
        original_rows=dict(F=frows,Uz=fields['Uz'],Q=fields['Q'],mean=fields['mean'])
        zjet=algebra.lift(op.Jet.variable(c,z,5))
        forcing=op.known_forcing(algebra,original_rows,zjet,rho,self.physical.delta)
        result=dict(source_family=self.family,Z=z,query=query,geometry=sv['geometry'],
            leading_Psi_rho_rows=[op.export(v) for v in rows['Psi']],leading_own_meanPsi_rho_rows=[op.export(v) for v in rows['meanPsi']],
            leading_factored_Uz_rho_rows=[op.export(v) for v in fields['Uz']],leading_factored_own_mean_rho_rows=[op.export(v) for v in fields['mean']],
            leading_Qaxis=op.export(fields['Qaxis'],3),leading_Qpsi_rho_rows=[op.export(v) for v in fields['Qpsi']],
            leading_factored_Q_rho_rows=[op.export(v) for v in fields['Q']],
            genuine_n1_forcing_from_factored_leading_rows={k:op.export(v,2) for k,v in forcing.items()},
            uniform_complex_factored_bounds={k:[op.export(v,0) for v in uniform[k]] for k in ('Psi','meanPsi','Qpsi','Uz','mean','Q')},
            uniform_complex_bound_proof=proof,
            own_mean_identity='rho*meanPsi=4*meanPsi_core(4)+integral_4^rho Psi_actual drho',
            leading_radial_recovery_identity='Q=Qaxis+epsilon*Qpsi; Qpsi=[2Z*Psi-(1-delta)Z*meanPsi-d*meanPsi_Z]/L',
            exact_original_pressure_and_epsilon_sectors_retained=True,actual_own_mean_history_not_comparison_mean=True,
            independent_raw_M_minus_uaxis_subtraction_not_used=True,ordinary_mixed_total_order_five_available=True,
            collar_mean_increment_bounded_with_original_positive_hb=True,n1_axial_pressure_source_integrals_still_open=True,
            component_weighted_n1_tail_still_open=True,full_n1_function_or_repair_not_promoted=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalLeadingMeanEpsilonPacket();state=(result,rows,fields,forcing,uniform)
        self._packets[id(token)]=(token,result,lp,p,state,op.snapshot(state))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalLeadingMeanEpsilonPacket or entry is None or entry[0] is not packet:raise ValueError('Live original mean epsilon packet required')
        self.local.report(entry[2]);self.inputs.report(entry[3])
        if op.snapshot(entry[4])!=entry[5]:raise ValueError('Leading own-mean factor packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_axial):
    owner=CurrentOriginalN1LeadingMeanEpsilonFactor(accepted_axial,require_checked=False)
    packets={name:owner.evaluate('.371',coordinate,chart) for name,coordinate,chart in
             (('core_axis','0','core'),('core_exit','4','core'),('actual_first_inlet','0','first'),('selected_actual_inner_exit','1','first'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(p) for name,p in packets.items()},input_hashes=owner.hashes,
             **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('LEADING_OWN_MEAN_EPSILON_FACTORS_COMPUTED',len(packets),flush=True)
    return owner,packets

"""Natural n1 component scales, full feedback, and correlated axial seeds.

Pressure uses the exact viscous FTC primitive of the original Qpsi.
Axial seeds integrate the axis-only forcing exactly. All other source and
feedback terms receive analytic bounds; none are assigned zero.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_leading_mean_epsilon_factor as mean
import lei_ren_part1_paper_compliant_current_original_n1_weighted_feedback_bound as feedback
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=mean.op;domain=mean.domain;weighted=mean.weighted;core=mean.core;ep=mean.ep
HERE,PREFIX,sha=mean.HERE,mean.PREFIX,mean.sha
NAME=PREFIX+'current_original_n1_natural_component_enclosure.json.gz'
RECEIPT=PREFIX+'current_original_n1_natural_component_enclosure_check.json'
GATE='current_original_n1_natural_component_feedback_and_selected_axial_pressure_enclosures_installed'
OPEN=mean.OPEN
EXTRA=((0,0,0,0),(0,0,0,-1),(0,0,0,-1),(0,0,1,-1),(0,0,0,0),(0,0,0,-1))
WEIGHTS=tuple(tuple(a+b for a,b in zip(v,t)) for v,t in zip(weighted.WEIGHTS,EXTRA))
SCALES=dict(Uz1=(0,0,0,-1),K1=(0,0,0,-1),P1=(0,0,1,-1),
            d_x_Uz1=(0,0,0,-1),mean1=(0,0,0,-1),Q1=(0,0,0,-1),V1=(0,0,0,-2))
FLOOR='1e-320'


def conjugate_extra(algebra,matrices):
    """T is constant in x and Z, so there is no new gauge derivative."""
    out={}
    for name,rows in matrices.items():
        out[name]={}
        for key,v in rows.items():
            i,j=(int(s)-1 for s in key.split(','))
            out[name][key]=algebra.shift(v,tuple(b-a for a,b in zip(EXTRA[i],EXTRA[j])))
    return out


def axial_majorant(w,wr,wrr,wz,wzz,wrz,z,d,L,r,delta):
    """Triangle bound of the executable T_-2 T_(-1-delta) operator."""
    a=1+delta
    A=w*(a*z)+wz*d+wr*(2*z*r)
    Az=w*a+wz*((a+2)*z)+wzz*d+wr*(2*r)+wrz*(2*z*r)
    Ar=wr*((a+2)*z)+wrz*d+wrr*(2*z*r)
    return (A*(2*z/L)+(Az/L+A*(2*delta*z/(L*L)))*d+Ar*(2*z*r/L))/L


def natural_majorants(owner):
    """Matrix on eta/4, all initial sources on eta/8, same outer eta."""
    c=owner.c;b=owner.bound;p=owner.physical;eta=b['common_tubes']['outer_Z_radius']
    r=b['actual_inner_geometry']['real_rho_upper_bound'];a=c.sqrt(r);cap=b['actual_inner_geometry']['numerical_cap_bounds_only']
    z=1+eta;d=1+z*z;dt=abs(p.delta);L=b['denominator_bounds']['L_modulus_lower']
    Fmax=b['original_axis_amplitude']['complex_modulus_upper']
    algebra=op.N1LogAlgebra(c,(owner.weighted.source.logh,c.ln(Fmax),2*p.logP,p.logLambda),[])
    lift=lambda v:algebra.lift(v,0)
    eps=lambda v,n=1:algebra.shift(v,(0,0,0,-n))
    uniform,meanproof=mean.complex_mean_bounds(algebra,p,b,owner.weighted.tube)
    psiold,psiproof=mean.axial.complex_psi_bounds(algebra,p,b,owner.weighted.tube)
    u=4*z+abs(p.j);Haxis=z*(1+dt)/2+d*u
    pn=b['original_pressure']['normalized_modulus_upper'];pz=2*pn/eta
    slope=lift(((1+2*z*u)*u*(1+dt)/2+4*Haxis)/(2*L))+algebra.shift(lift((d*pz+2*z*pn*(1+dt))/(2*L)),(0,0,1,0))
    corr=domain.upper(c,4*p.correction)
    pc=[slope*r+lift(corr),slope+lift(2*corr),lift(8*corr)]
    mc=[slope*(r/2)+lift(corr),slope/2+lift(corr),lift(8*corr/3)]
    corefields=dict(Psi=pc,meanPsi=mc)
    collarfields=dict(Psi=uniform['Psi'],meanPsi=uniform['meanPsi'])
    qa=lift(((1+dt)*z*u+4*d)/L)
    def fields(rows):
        P,M=rows['Psi'],rows['meanPsi']
        Qp=[(v*(2*z)+m*((1+dt)*z)+m*(d/(eta/4)))/L for v,m in zip(P,M)]
        return dict(Psi=P,meanPsi=M,Qpsi=Qp,
                    Uz=[lift(u)+eps(P[0]),eps(P[1]),eps(P[2])],
                    Q=[qa+eps(Qp[0]),eps(Qp[1]),eps(Qp[2])])
    cf,bf=fields(corefields),fields(collarfields)
    def phi_rows(chart):
        return [domain.WidthBound(c,{int(k):v['coefficient'] for k,v in row.items()})
                for row in b['actual_'+chart+'_radial_bounds']['Phi']]
    phic,phib=phi_rows('core'),phi_rows('collar')
    phi=phic[0].numeric_upper(cap)+phib[0].numeric_upper(cap)
    pr=phic[1].numeric_upper(cap)+phib[1].numeric_upper(cap)
    P=[x+y for x,y in zip(cf['Psi'],bf['Psi'])]
    w=lift(u)+eps(P[0]);wr=eps(P[1]);wz=lift(4)+eps(P[0]/(eta/4))
    q=cf['Q'][0]+bf['Q'][0];E=w*z+lift(c.mpf('.5'));B=q+(lift(1)+w*(2*z))/L
    H=lift(z*(1+dt)/2)+w*d;Cf=phi+r*pr
    Ta=((2+dt)*z*phi+d*phi/(eta/4)+2*z*r*pr)/L
    Zw=(w*((1+dt)*z)+wz*d+wr*(2*z*r))/L
    m=d*psiproof['H_over_H_squared_plus_sigma_squared_modulus_upper']
    f2=lambda v:algebra.shift(lift(v),(0,2,0,0))
    old0={
        '4,1':f2(4*a*phi),'5,1':eps((q+E*((2+dt)/L))*2)+psiold['residual_gauge'],
        '5,2':eps(lift(2*Ta+2*(1+dt)*z*Cf/L),2)+eps(lift(2*m*phi)),
        '5,3':eps(lift(2*(1+dt)*z*Cf/L),2),'5,5':eps(B*a),
        '6,1':eps(f2(8*z*r*phi/L)),
        '6,2':eps((E*((1+dt)/L)+Zw+wr*((1+dt)*z*r/L))*2),
        '6,3':eps(wr*(2*(1+dt)*z*r/L)),'6,4':eps(lift(4*z/L)),'6,6':eps(B*a)}
    old1={'5,1':eps(H*(2/L)),'5,2':eps(lift(2*d*Cf/L),2),'5,3':eps(lift(2*d*Cf/L),2),
          '6,2':eps((H+wr*(d*r))*(2/L)),'6,3':eps(wr*(2*d*r/L)),'6,4':eps(lift(2*d/L))}
    mats=conjugate_extra(algebra,dict(C0entries=old0,C1entries=old1))
    C0=sum(mats['C0entries'].values(),lift(0));C1=sum(mats['C1entries'].values(),lift(0))
    axis_bound=axial_majorant(lift(u),lift(0),lift(0),lift(4),lift(0),lift(0),z,d,L,r,dt)
    def axial_source(rows):
        P=rows['Psi'];gap=eta/4
        return axial_majorant(P[0],P[1],P[2],P[0]/gap,P[0]*(2/gap**2),P[1]/gap,z,d,L,r,dt)
    def source(rows,ph):
        q0,q1,q2=rows['Q'];qz=q0/(eta/8);w0=rows['Uz'][0]
        inertial=((q0+qz*((1+dt)*z/2)+q1*r)/L+
                  q0*(q0/2+q1*r)+w0*(q0*(2*z)+qz*d+q1*(2*z*r))/L)/2
        viscous=rows['Qpsi'][1]*2+rows['Qpsi'][2]*r
        Np=inertial+viscous
        g4=algebra.shift(Np*(2*a),(0,0,-1,0))
        g4rest=algebra.shift(inertial*(2*a),(0,0,-1,0))
        g6rest=eps(axial_source(rows)*2)+eps(Np*(4*z*r/L))
        # This independently admitted g5 bound has the same normalization.
        # Keep every amplitude/gauge sector and only shrink the source tube.
        return dict(g4=g4,g4rest=g4rest,g6=axis_bound*2+g6rest,g6rest=g6rest)
    sc,sb=source(cf,phic),source(bf,phib)
    old=feedback.original_majorants(owner.feedback)
    def rebase(v):return op.FactoredJet(algebra,v.terms,v.order)
    ic={k:v*2 for k,v in sc.items()};ib={k:algebra.width(v)*(a/4) for k,v in sb.items()}
    ic['g5']=rebase(old['core_source_L1']['5']);ib['g5']=rebase(old['collar_source_L1']['5'])
    if any(p[0]<0 for v in ib.values() for p in v.terms):raise ValueError('Original collar inverse width not integrated')
    L1=sum((ic[k]+ib[k] for k in ('g4','g5','g6')),lift(0))
    cover=lambda v:feedback.directed_sector_cover(v,floor=FLOOR)
    source_cover=cover(L1);green=domain.upper(c,c.exp(2*a));M0=domain.upper(c,green*source_cover['modulus_upper'])
    c0,c1=cover(C0),cover(C1)
    tail=feedback.all_word_tail(c,M0,a,c0['modulus_upper'],c1['modulus_upper'],eta/16,0)
    tail.update(source_M0_tube_radius=eta/8,value_tube_radius=eta/16,
        ordinary_Z1_Z2_tube_radius=eta/32,ordinary_Z_final_Cauchy_gap=eta/32,
        ordinary_Z1_remainder_upper=domain.upper(c,tail['value_remainder_upper']/(eta/32)),
        ordinary_Z2_remainder_upper=domain.upper(c,2*tail['value_remainder_upper']/(eta/32)**2))
    pressure_rest=cover(ic['g4rest']+ib['g4rest'])
    axial_rest=cover(ic['g6rest']+ib['g6rest'])
    initial_errors=dict(P1=pressure_rest['modulus_upper'],axial=domain.upper(c,green*axial_rest['modulus_upper']))
    proof=dict(original_outer_eta=eta,matrix_bound_tube=eta/4,source_M0_tube=eta/8,
        original_hb_F0_Pstar_Lambda_parameters_unchanged=True,natural_weight_exponents=[list(p) for p in WEIGHTS],
        extra_T_is_x_and_Z_independent=True,no_new_weight_gauge=True,principal_Akin_Achi_unchanged=True,
        original_own_mean_and_Q_epsilon_sectors_used=True,Q_Z_source_Cauchy_gap=eta/8,
        Psi_Z_source_Cauchy_gap=eta/4,source_axis_Uz_Z1_exact_four=True,source_axis_Uz_Z2_exact_zero=True,
        all_original_pressure_sectors_retained=True,source_L1_sums_g4_g5_g6=True,
        original_core_dx_length=2,actual_collar_dx_measure='hb*x*ds/2; phase<=1/2',
        inverse_hb_cancelled_before_norm=True,Fmax_is_only_complex_norm_proxy=True,
        directed_floor_only_after_original_log_proof=FLOOR,
        mean_factor_proof=meanproof,remaining_gauge_proof=psiproof,
        remaining_pressure_inertial_source_not_set_zero=True,remaining_axial_Psi_and_pressure_sources_not_set_zero=True)
    return dict(C0=C0,C1=C1,**mats,C0cover=c0,C1cover=c1,
        core_L1=ic,collar_L1=ib,source_cover=source_cover,M0=M0,green=green,a=a,eta=eta,
        tail=tail,pressure_rest=pressure_rest,axial_rest=axial_rest,initial_errors=initial_errors,proof=proof)


def seeds(algebra,rows,fields,Z,rho,delta):
    c=algebra.ctx;z=algebra.lift(op.Jet.variable(c,Z,5));u=algebra.lift(rows['axis']['u']);zero=algebra.lift(0,5)
    C=op.truncate(-2*op.axial_second([u,zero,zero],z,rho,delta,-1-delta),2)
    U=C*(rho/4);K=C*(-rho/8);M=C*(rho/8)
    Q=op.truncate((U*(z*(1-delta))-K*(z*(1+delta))-op.dz(M)*(1-z*z))/(1-z*z*delta),1)
    P=algebra.shift(op.truncate(fields['Qpsi'][0]+fields['Qpsi'][1]*rho,2),(0,0,-1,0))
    return dict(Uz1=U,K1=K,P1=P,d_x_Uz1=C*(c.sqrt(rho)/2),mean1=M,Q1=Q,V1=Q*rho)


def natural_budget(value,error,k):
    result=core.integral_budget(value,error,k)
    result['relative_natural_component_width_certified']=result.pop('relative_initial_integral_width_certified')
    result['remaining_B0_B1_feedback_included']=True
    result['all_nonseed_sources_included']=True
    return result


@dataclass(frozen=True,eq=False)
class OriginalN1NaturalComponentPacket:
    kind:str='point'
    temporal_hierarchy_order:int=1


class CurrentOriginalN1NaturalComponentEnclosure:
    @source_precision
    def __init__(self,accepted_mean,accepted_feedback,require_checked=True):
        if (type(accepted_mean) is not mean.CurrentOriginalN1LeadingMeanEpsilonFactor or not accepted_mean.acceptance_loaded or
            type(accepted_feedback) is not feedback.CurrentOriginalN1WeightedFeedbackBound or not accepted_feedback.acceptance_loaded or
            accepted_mean.axial is not accepted_feedback.axial):raise ValueError('Same accepted original mean/feedback parents required')
        accepted_mean.assert_graph();accepted_feedback.assert_graph()
        self.mean=accepted_mean;self.feedback=accepted_feedback;self.c=accepted_mean.c;self.bound=accepted_mean.bound
        self.physical=accepted_mean.physical;self.weighted=accepted_mean.weighted;self.family=copy.deepcopy(accepted_mean.family)
        self.hashes=dict(accepted_mean.hashes)
        for name,digest in accepted_feedback.hashes.items():domain.parent.parent.original.bind(self.hashes,name,digest)
        for name in (mean.NAME,mean.RECEIPT,feedback.NAME,feedback.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._parents=(accepted_mean,accepted_feedback);self._identities=(self.c,self.bound,self.physical,self.weighted)
        self._bound=op.snapshot(self.bound);self._packets={};self._data=None;self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Natural component receipt differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Natural component input changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if ((self.mean,self.feedback)!=self._parents or (self.c,self.bound,self.physical,self.weighted)!=self._identities or
            self.mean.axial is not self.feedback.axial or not self.mean.acceptance_loaded or not self.feedback.acceptance_loaded or
            self.hashes!=self._hashes or self.family!=self.mean.family or self.family!=self.feedback.family or
            op.snapshot(self.bound)!=self._bound):raise ValueError('Natural component graph changed')
        self.mean.assert_graph();self.feedback.assert_graph()

    @source_precision
    def domain(self):
        self.assert_graph()
        if self._data is None:
            self._data=natural_majorants(self);self._data_snapshot=op.snapshot(self._data)
        if op.snapshot(self._data)!=self._data_snapshot:raise ValueError('Natural component bound changed')
        v=self._data
        result=dict(source_family=self.family,uniform_matrix_entry_sectors={k:{q:op.export(x,0) for q,x in v[k].items()} for k in ('C0entries','C1entries')},
            uniform_C0_cover=v['C0cover'],uniform_C1_cover=v['C1cover'],uniform_source_L1_cover=v['source_cover'],
            uniform_initial_M0=v['M0'],all_word_feedback_tail=v['tail'],uniform_proof=v['proof'],
            uniform_core_source_L1={k:op.export(x,0) for k,x in v['core_L1'].items()},
            uniform_collar_source_L1={k:op.export(x,0) for k,x in v['collar_L1'].items()},
            pressure_inertial_integral_cover=v['pressure_rest'],axial_nonaxis_integral_cover=v['axial_rest'],
            natural_initial_source_errors=v['initial_errors'],**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1NaturalComponentPacket('domain');self._packets[id(token)]=(token,result,None,None,op.snapshot(result))
        return token

    @source_precision
    def evaluate(self,Z,coordinate,chart='first'):
        self.assert_graph();began=time.monotonic();dp=self.domain();data=self._data;c=self.c
        mp=self.mean.evaluate(Z,coordinate,chart);mv=self.mean.report(mp);entry=self.mean._packets[id(mp)];rows,fields=entry[4][1:3]
        algebra=fields['Uz'][0].algebra;rho=mv['geometry'].get('rho',mv['geometry'].get('rho_enclosure_only'));z=mv['Z']
        values=seeds(algebra,rows,fields,z,rho,self.physical.delta)
        # Bounds are analytic on eta/16. Derivatives use the final eta/32
        # gap; factorials belong to ordinary derivatives, then divide for jets.
        eta=data['eta'];gap=eta/32;tail=data['tail']['value_remainder_upper']
        ev=domain.upper(c,data['initial_errors']['axial']+tail)
        pe=domain.upper(c,data['initial_errors']['P1']+tail)
        zmax=1+eta;dmax=1+zmax*zmax;L=self.bound['denominator_bounds']['L_modulus_lower']
        qe=domain.upper(c,(2*(1+abs(self.physical.delta))*zmax*ev+2*dmax*ev/gap)/L)
        norms=dict(Uz1=ev,K1=ev,P1=pe,d_x_Uz1=ev,mean1=2*ev,Q1=qe,V1=qe*rho)
        axis_zero=ep(rho)[1]==0
        errors={k:algebra.lift(0,v.order) if axis_zero else algebra.lift(op.Jet(c,[domain.upper(c,norms[k]/gap**j) for j in range(v.order+1)]))
                for k,v in values.items()}
        physical={k:algebra.shift(v,SCALES[k]) for k,v in values.items()}
        physical_errors={k:algebra.shift(v,SCALES[k]) for k,v in errors.items()}
        budgets={k:{'Z'+str(j):natural_budget(v,errors[k],j) for j in (0,1)} for k,v in values.items()}
        pb={k:{'Z'+str(j):natural_budget(v,physical_errors[k],j) for j in (0,1)} for k,v in physical.items()}
        result=dict(source_family=self.family,Z=z,query=mv['query'],geometry=mv['geometry'],
            normalized_signed_seed={k:op.export(v) for k,v in values.items()},normalized_total_error={k:op.export(v) for k,v in errors.items()},
            physical_signed_seed={k:op.export(v) for k,v in physical.items()},physical_total_error={k:op.export(v) for k,v in physical_errors.items()},
            normalized_component_budgets=budgets,physical_component_budgets=pb,
            passing_normalized_width_rows=sum(x['relative_natural_component_width_certified'] for v in budgets.values() for x in v.values()),
            passing_physical_width_rows=sum(x['relative_natural_component_width_certified'] for v in pb.values() for x in v.values()),
            requested_component_width_rows=14,exact_zero_axis_data=axis_zero,
            viscous_pressure_FTC_identity='Pvisc/(epsilon*Pstar^2)=(Qpsi(rho)-Qpsi(0)+rho*Qpsi_r)/Pstar^2; Qpsi(0)=0',
            axial_signed_seed_identity='C=-2*T_-2*T_(-1-delta)uaxis; U=C*rho/4; K=-C*rho/8; dxU=C*sqrt(rho)/2',
            radial_recovery_identity='Q1=[(1-delta)Z*Uz1-(1+delta)Z*K1-d*d_Z(Uz1+K1)]/L; V1=rho*epsilon*Q1',
            pressure_inertial_and_axial_nonaxis_sources_bounded_not_omitted=True,
            total_error_contains_all_six_row_feedback=True,original_own_mean_and_corrected_source_coordinate_used=True,
            axial_pressure_mean_Q_V_component_scales_Z_independent=True,physical_scale_derivative_gauge_not_added=True,
            original_zero_axis_regular_solution_used=axis_zero,radial_recovery_uses_own_n1_mean=True,
            point_widths_not_whole_Z_relative_claim=True,full_coefficient_function_and_equation_contracts_still_open=True,
            n1_own_five_moment_repair_and_n2_still_open=True,execution_seconds=time.monotonic()-began,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1NaturalComponentPacket();state=(values,errors,physical,physical_errors)
        self._packets[id(token)]=(token,result,mp,state,op.snapshot((result,state)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1NaturalComponentPacket or entry is None or entry[0] is not packet:raise ValueError('Live natural component packet required')
        if entry[2] is None:
            if op.snapshot(entry[1])!=entry[4] or op.snapshot(self._data)!=self._data_snapshot:raise ValueError('Natural bound packet changed')
        else:
            self.mean.report(entry[2])
            if op.snapshot((entry[1],entry[3]))!=entry[4]:raise ValueError('Natural component packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_mean,accepted_feedback):
    owner=CurrentOriginalN1NaturalComponentEnclosure(accepted_mean,accepted_feedback,require_checked=False)
    dp=owner.domain()
    packets={name:owner.evaluate('.371',value,chart) for name,value,chart in
             (('core_axis','0','core'),('core_exit','4','core'),('actual_first_inlet','0','first'),('selected_actual_inner_exit','1','first'))}
    raw=dict(source_family=owner.family,domain=owner.report(dp),samples={k:owner.report(v) for k,v in packets.items()},
             input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_NATURAL_COMPONENT_ENCLOSURES_COMPUTED',len(packets),flush=True)
    return owner,dp,packets

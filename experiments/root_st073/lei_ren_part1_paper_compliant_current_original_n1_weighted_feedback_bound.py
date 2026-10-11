"""Uniform residual matrix and all-word tail for the original n1 system.

The source amplitude is cancelled symbolically before boxing Phi derivatives.
Positive norm sectors use Fmax, explicitly distinguished from pointwise F0.
This admits feedback bounds; it does not publish a solved coefficient owner.
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
NAME=PREFIX+'current_original_n1_weighted_feedback_bound.json.gz'
RECEIPT=PREFIX+'current_original_n1_weighted_feedback_bound_check.json'
GATE='current_original_n1_full_weighted_residual_and_all_word_feedback_bounds_installed'
OPEN=axial.OPEN


def residual_matrix(algebra,phi,w,q,z,rho,delta,gradient,gauge,order=1):
    """Exact S-normalized matrices, with Akin and -2chi removed once."""
    c=algebra.ctx;x=c.sqrt(rho);d=1-z*z;L=1-z*z*delta
    f,fr=phi[:2];uz,ur=w[:2];H=z*((1-delta)/2)+d*uz
    E=z*uz-c.mpf('.5');B=q+(1-2*z*uz)/L;Cf=f+fr*rho
    Ta=(f*(z*(-2-delta))+op.dz(f)*d-fr*(z*(2*rho)))/L
    Zw=(uz*(z*(-1-delta))+op.dz(uz)*d-ur*(z*(2*rho)))/L
    m=algebra.lift(gradient)*d/L
    eps=lambda v,p=1:algebra.shift(v,(0,0,0,-p))
    f2=lambda v:algebra.shift(v*weighted.relative_amplitude(algebra,gradient,2,order),(0,2,0,0))
    B0={
        '4,1':f2(f*(4*x)),
        '5,1':eps(2*(q+E*(-2+delta)/L))+gauge,
        '5,2':eps(2*Ta-Cf*(z*(2*(-1+delta)))/L,2)-eps(2*m*f),
        '5,3':eps(-Cf*(z*(2*(1+delta)))/L,2),
        '5,5':eps(B*x),
        '6,1':eps(f2(-8*z*rho*f/L)),
        '6,2':eps(2*(E*(-1+delta)/L+Zw-ur*(z*((-1+delta)*rho))/L)),
        '6,3':eps(-ur*(z*(2*(1+delta)*rho))/L),
        '6,4':eps(-4*z/L),'6,6':eps(B*x)}
    B1={'5,1':eps(2*H/L),'5,2':eps(-2*d*Cf/L,2),'5,3':eps(-2*d*Cf/L,2),
        '6,2':eps(2*(H-d*rho*ur)/L),'6,3':eps(-2*d*rho*ur/L),'6,4':eps(2*d/L)}
    return {name:{k:op.truncate(v,order) for k,v in rows.items()} for name,rows in (('B0',B0),('B1',B1))}


def directed_sector_cover(value,k=0,floor='1e-240'):
    """Outward absolute cover; a tiny floor is used only after a log proof."""
    c=value.ctx;total=c.mpf(0);ledger=[];small=c.mpf(floor);cut=c.ln(small)
    if ep(small)[0]<=0 or ep(small)[1]>=1:raise ValueError('Positive small cover required')
    for powers,row in sorted(value.terms.items()):
        mag=domain.upper(c,abs(row[k]))
        if ep(mag)[1]==0:continue
        log=core.budget.log_ratio(c,value.algebra.logs,powers,(0,0,0,0),mag,c.mpf(1))
        if ep(log)[1]<=ep(cut)[0]:cover=domain.upper(c,small);rule='proved log<=log(floor)'
        elif ep(log)[1]<10000:cover=domain.upper(c,c.exp(log));rule='directed exp(log)'
        else:raise ValueError('Unresolved large norm sector')
        total+=cover
        ledger.append(dict(norm_exponents=list(powers),coefficient_modulus_upper=mag,
                           log_modulus_upper=log,cover=cover,rule=rule))
    return dict(modulus_upper=domain.upper(c,total),sectors=ledger,
                tiny_cover_floor=small,all_covers_have_directed_log_proof=True)


def width_norm(algebra,row,lambda_power=0):
    return op.FactoredJet(algebra,{(p,0,0,lambda_power):op.Jet.constant(algebra.ctx,v,0)
                                  for p,v in row.rows.items()},0)


def original_majorants(owner):
    """Matrix on eta/2; genuine normalized source L1 on eta/4."""
    c=owner.c;b=owner.bound;physical=owner.physical
    eta=b['common_tubes']['outer_Z_radius'];r=b['actual_inner_geometry']['real_rho_upper_bound'];a=c.sqrt(r)
    z=1+eta;d=1+z*z;dt=abs(physical.delta);L=b['denominator_bounds']['L_modulus_lower']
    Fmax=b['original_axis_amplitude']['complex_modulus_upper']
    # Fmax is a norm proxy, never a selected source amplitude.
    logs=(owner.weighted.source.logh,c.ln(Fmax),2*physical.logP,physical.logLambda)
    algebra=op.N1LogAlgebra(c,logs,[]);lift=lambda v:algebra.lift(v,0)
    eps=lambda v,p=1:algebra.shift(v,(0,0,0,-p))
    f2=lambda v:algebra.shift(lift(v),(0,2,0,0))
    def radial(label,chart):
        return [domain.WidthBound(c,{int(p):entry['coefficient'] for p,entry in row.items()})
                for row in b['actual_'+chart+'_radial_bounds'][label]]
    phic,phib=radial('Phi','core'),radial('Phi','collar')
    wc,wb=radial('Uz','core'),radial('Uz','collar')
    mc,mb=radial('mean','core'),radial('mean','collar')
    # Only radial orders 0 and 1 enter the matrix; neither has inverse hb.
    scalar=lambda rows,i:rows[i].numeric_upper(b['actual_inner_geometry']['numerical_cap_bounds_only'])
    phi=scalar(phic,0)+scalar(phib,0);pr=scalar(phic,1)+scalar(phib,1)
    w=scalar(wc,0)+scalar(wb,0);wr=scalar(wc,1)+scalar(wb,1);mean=scalar(mc,0)+scalar(mb,0)
    rz=eta/2;pz=phi/rz;wz=w/rz
    q=(2*z*w+(1+dt)*z*mean+d*mean/rz)/L
    Cf=phi+r*pr;Ta=((2+dt)*z*phi+d*pz+2*z*r*pr)/L
    Zw=((1+dt)*z*w+d*wz+2*z*r*wr)/L
    E=z*w+c.mpf('.5');B=q+(1+2*z*w)/L;H=z*(1+dt)/2+d*w
    uniform,psi_proof=axial.complex_psi_bounds(algebra,physical,b,owner.weighted.tube)
    gauge=op.truncate(uniform['residual_gauge'],0)
    K=psi_proof['H_over_H_squared_plus_sigma_squared_modulus_upper'];m=d*K
    C0entries={
        '4,1':f2(4*a*phi),
        '5,1':eps(lift(2*(q+(2+dt)*E/L)))+gauge,
        '5,2':eps(lift(2*Ta+2*(1+dt)*z*Cf/L),2)+eps(lift(2*m*phi)),
        '5,3':eps(lift(2*(1+dt)*z*Cf/L),2),'5,5':eps(lift(a*B)),
        '6,1':eps(f2(8*z*r*phi/L)),
        '6,2':eps(lift(2*((1+dt)*E/L+Zw+(1+dt)*z*r*wr/L))),
        '6,3':eps(lift(2*(1+dt)*z*r*wr/L)),'6,4':eps(lift(4*z/L)),'6,6':eps(lift(a*B))}
    C1entries={'5,1':eps(lift(2*H/L)),
               '5,2':eps(lift(2*d*Cf/L),2),'5,3':eps(lift(2*d*Cf/L),2),
               '6,2':eps(lift(2*(H+d*r*wr)/L)),
               '6,3':eps(lift(2*d*r*wr/L)),'6,4':eps(lift(2*d/L))}
    # Sum all entries: a conservative infinity norm with retained sectors.
    C0=sum(C0entries.values(),lift(0));C1=sum(C1entries.values(),lift(0))
    def source_bounds(phi_rows,w_rows,mean_rows):
        # This helper uses actual own mean, Cauchy Q on eta/2, and its
        # subsequent derivatives on eta/4. Fmax=1 supplies normalized Phi.
        h=domain.hierarchy_bounds(c,phi_rows,w_rows,mean_rows,c.mpf(1),physical.delta,eta,r,L)
        p0=phi_rows[0];p1=phi_rows[1]
        cross=((p0*((5+dt)*z)+p0*(2*d/(eta/4))+p1*(4*r*z))*(m/L)+p0*(d*m/(eta/4)/L))
        g5=width_norm(algebra,p0*(2*m*m))+width_norm(algebra,cross*2,-1)+width_norm(algebra,h['N1theta']*2,-2)
        g4=width_norm(algebra,h['N1p_inertial']*(2*a),-1)+width_norm(algebra,h['N1p_viscous_per_Lambda']*(2*a))
        g6=width_norm(algebra,h['N1z']*2,-1)+width_norm(algebra,h['N1p_inertial']*(4*z*r/L),-2)+width_norm(algebra,h['N1p_viscous_per_Lambda']*(4*z*r/L),-1)
        return {'4':g4,'5':g5,'6':g6}
    gc=source_bounds(phic,wc,mc);gb=source_bounds(phib,wb,mb)
    integrated_core={k:v*2 for k,v in gc.items()}
    integrated_collar={k:algebra.width(v)*(a/4) for k,v in gb.items()}
    if any(p[0]<0 for v in integrated_collar.values() for p in v.terms):raise ValueError('Actual width not cancelled')
    source_L1=sum(integrated_core.values(),lift(0))+sum(integrated_collar.values(),lift(0))
    source_cover=directed_sector_cover(source_L1);green=domain.upper(c,c.exp(2*a))
    M0=domain.upper(c,green*source_cover['modulus_upper'])
    proof=dict(original_eta=eta,matrix_bound_tube=eta/2,source_initial_bound_tube=eta/4,
        principal_Green_infinity_bound='exp(2*(x-s))',principal_row_norm_upper=c.mpf(2),
        regular_diagonal_exponents=[0,0,2,0,3,1],regular_zero_axis_data=True,
        normalized_m_modulus_upper_on_eta_over2=m,m_Z_modulus_upper_on_eta_over4=m/(eta/4),
        Phi_and_Uz_Z_Cauchy_gap_for_matrix=eta/2,source_derivative_Cauchy_gap=eta/4,
        norm_log_basis_names=['original_hb','uniform_Fmax_norm_proxy','original_Pstar_squared','original_Lambda'],
        norm_log_bases=list(logs),Fmax_is_uniform_complex_bound_not_point_amplitude=True,
        source_L1_uses_sum_of_rows_not_max=True,core_dx_length=2,actual_dx_measure='hb*x*ds/2; phase<=1/2',
        inverse_width_cancelled_before_source_norm=True,pressure_and_original_gauge_sectors_retained=True,
        selected_completed_inner_interval_contained_in_bounded_canonical_core_collar=True,
        inherited_psi_bound_proof=psi_proof)
    return dict(C0=C0,C1=C1,C0entries=C0entries,C1entries=C1entries,source_L1=source_L1,
                core_source_L1=integrated_core,collar_source_L1=integrated_collar,M0=M0,
                source_cover=source_cover,green=green,a=a,eta=eta,proof=proof)


def all_word_tail(c,M0,a,C0,C1,loss,N=0):
    """All derivative words: equal nested Cauchy gaps, then time simplex."""
    if isinstance(N,bool) or not isinstance(N,int) or N<0:raise ValueError('N>=0 required')
    M0,a,C0,C1,loss=(c.mpf(v) for v in (M0,a,C0,C1,loss))
    if any(ep(v)[0]<0 for v in (M0,C0,C1)) or ep(loss)[0]<=0:raise ValueError('Positive norm premises required')
    Q=domain.upper(c,c.e*a*(C0+C1/loss))
    if ep(Q)[1]>=1:raise ValueError('All-word feedback ratio must be below one')
    value=domain.upper(c,M0*c.exp(2*a)*Q**(N+1)/(1-Q))
    return dict(retained_indices=[0,N],ratio_upper=Q,value_remainder_upper=value,
                iteration_Cauchy_loss=loss,formula='M0*exp(2a)*Q^(N+1)/(1-Q); Q=e*a*(C0+C1/loss)',
                covers_all_derivative_words=True,old_nonadjacent_B1_word_count_not_reused=True,
                Green_and_coefficient_derivatives_handled_by_nested_Cauchy=True,
                ordered_radial_simplex_factor='a^k/k!',Green_factors_telescope=True)


@dataclass(frozen=True,eq=False)
class OriginalN1WeightedFeedbackBoundPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1WeightedFeedbackBound:
    @source_precision
    def __init__(self,accepted_axial,require_checked=True):
        if type(accepted_axial) is not axial.CurrentOriginalN1AxialEpsilonFactor or not accepted_axial.acceptance_loaded:
            raise ValueError('Accepted original axial epsilon owner required')
        accepted_axial.assert_graph();self.axial=accepted_axial;self.weighted=accepted_axial.weighted
        self.inputs=accepted_axial.inputs;self.local=accepted_axial.local;self.physical=accepted_axial.physical;self.c=accepted_axial.c
        self.bound=accepted_axial.bound;self.family=copy.deepcopy(accepted_axial.family);self.hashes=dict(accepted_axial.hashes)
        for name in (axial.NAME,axial.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._owner=accepted_axial;self._identities=(self.weighted,self.inputs,self.local,self.physical,self.c)
        self._bound=op.snapshot(self.bound);self.acceptance_loaded=False;self._packets={}
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Weighted feedback bound receipt differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Weighted feedback bound source changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.axial is not self._owner or (self.weighted,self.inputs,self.local,self.physical,self.c)!=self._identities
            or self.inputs is not self.axial.inputs or self.weighted is not self.axial.weighted
            or self.local is not self.axial.local or self.physical is not self.axial.physical or self.c is not self.axial.c
            or self.hashes!=self._hashes or self.family!=self.axial.family or op.snapshot(self.bound)!=self._bound
            or not self.axial.acceptance_loaded):raise ValueError('Weighted feedback bound graph changed')
        return self.axial.assert_graph()

    @source_precision
    def domain(self,N=0):
        self.assert_graph();began=time.monotonic();c=self.c;data=original_majorants(self)
        C0cover=directed_sector_cover(data['C0']);C1cover=directed_sector_cover(data['C1'])
        tail=all_word_tail(c,data['M0'],data['a'],C0cover['modulus_upper'],C1cover['modulus_upper'],data['eta']/8,N)
        tail.update(value_tube_radius=data['eta']/8,source_M0_tube_radius=data['eta']/4,
                    ordinary_Z1_tube_radius=data['eta']/16,ordinary_Z1_final_Cauchy_gap=data['eta']/16,
                    ordinary_Z1_remainder_upper=domain.upper(c,tail['value_remainder_upper']/(data['eta']/16)))
        result=dict(source_family=self.family,uniform_matrix_entry_norm_sectors={name:{k:op.export(v,0) for k,v in data[name].items()} for name in ('C0entries','C1entries')},
            uniform_residual_C0_norm_sectors=op.export(data['C0'],0),uniform_B1_norm_sectors=op.export(data['C1'],0),
            uniform_residual_C0_absolute_cover=C0cover,uniform_B1_absolute_cover=C1cover,
            uniform_core_source_L1_sectors={k:op.export(v,0) for k,v in data['core_source_L1'].items()},
            uniform_actual_collar_source_L1_sectors={k:op.export(v,0) for k,v in data['collar_source_L1'].items()},
            uniform_source_L1_cover=data['source_cover'],uniform_initial_resummed_source_M0=data['M0'],
            uniform_analytic_proof=data['proof'],all_word_feedback_tail=tail,
            remaining_feedback_bound_includes_all_six_rows=True,actual_feedback_iteration_prefix_not_computed=True,
            solved_coefficient_function_contracts_still_open=True,component_relative_accuracy_not_promoted=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1WeightedFeedbackBoundPacket();self._packets[id(token)]=(token,result,data,op.snapshot((result,data)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1WeightedFeedbackBoundPacket or entry is None or entry[0] is not packet:raise ValueError('Live weighted bound packet required')
        if op.snapshot((entry[1],entry[2]))!=entry[3]:raise ValueError('Weighted bound packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_axial):
    owner=CurrentOriginalN1WeightedFeedbackBound(accepted_axial,require_checked=False);packet=owner.domain()
    raw=dict(report=owner.report(packet),input_hashes=owner.hashes,source_family=owner.family,
             **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_WEIGHTED_FEEDBACK_BOUND_COMPUTED',flush=True)
    return owner,packet

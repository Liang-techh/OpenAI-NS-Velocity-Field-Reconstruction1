"""Genuine n1 operator functions on the core and ACTUAL first collar.

Original own field/moment functions are inputs. The comparison equals the
core on this collar, but is never substituted for the actual field. All
source width, amplitude and Lambda factors remain algebraic log sectors.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain as domain
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import FactoredAlgebra,FactoredJet
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as Jet
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=domain.HERE,domain.PREFIX,domain.sha
NAME=PREFIX+'current_original_n1_inner_operator.json.gz'
RECEIPT=PREFIX+'current_original_n1_inner_operator_check.json'
GATE='current_original_n1_actual_inner_operator_input_enclosures_installed'
OPEN=domain.OPEN
ep=domain.ep
SAMPLES=(('core_axis','core','0','.371'),('core_interior','core','2','.371'),
         ('first_inlet','first','0','.371'),('first_interior','first','.25','.371'),
         ('first_keep','first','.5','.371'))


class N1LogAlgebra(FactoredAlgebra):
    """Bases are (hb, F0_anchor, Pstar^2, Lambda), fixed during Z jets."""
    def lift(self,value,order=3):
        if isinstance(value,Jet) and value.ctx is not self.ctx:
            raise ValueError('Original interval coefficient context required')
        return super().lift(value,order)

    def resolve(self,*args,**kwargs):raise ValueError('N1 log sectors may not be materialized')
    def resolve_coefficient(self,*args,**kwargs):raise ValueError('N1 log sectors may not be materialized')


def truncate(value,order):
    if order>value.order:raise ValueError('Missing genuine source derivative')
    return FactoredJet(value.algebra,value.terms,order)


def dz(value):
    if value.order<1:raise ValueError('Missing genuine source Z derivative')
    return value.Zderivative()


def export(value,order=None):
    if order is not None:value=truncate(value,order)
    c=value.ctx
    return dict(ordinary_Z_derivative_order=value.order,terms=[dict(
        source_exponents=list(powers),
        ordinary_Z_derivatives=[math.factorial(k)*v for k,v in enumerate(jet.coefficients)],
        combined_source_log=sum((value.algebra.logs[i]*p for i,p in enumerate(powers) if p),c.mpf(0)))
        for powers,jet in sorted(value.terms.items())])


def fingerprint(value):
    if isinstance(value,FactoredJet):return export(value)
    if isinstance(value,Jet):return dict(coefficients=list(value.coefficients))
    if isinstance(value,dict):return {k:fingerprint(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [fingerprint(v) for v in value]
    return value


def snapshot(value):
    return copy.deepcopy(domain.parent.parent.first.encode(
        domain.parent.parent.original.serialized(fingerprint(value))))


def axial_second(rows,z,rho,delta,weight):
    """Exact Z_(a-1+delta)Z_a, ordinary rho rows, no Taylor-degree proxy."""
    d=1-z*z;L=1-z*z*delta
    f,fr,frr=rows;fz=dz(f);frz=dz(fr)
    A=f*(z*weight)+fz*d-fr*(z*(2*rho))
    Az=f*weight+fz*(z*(weight-2))+dz(fz)*d-fr*(2*rho)-frz*(z*(2*rho))
    Ar=fr*(z*(weight-2))+frz*d-frr*(z*(2*rho))
    return (A*(z*(weight-1+delta))/L
            +d*(Az/L+A*(z*(2*delta))/(L*L))-Ar*(z*(2*rho))/L)/L


def known_forcing(algebra,rows,z,rho,delta):
    f,w,q=rows['F'],rows['Uz'],rows['Q'];d=1-z*z;L=1-z*z*delta
    theta=-axial_second(f,z,rho,delta,-2-delta)
    axial=-axial_second(w,z,rho,delta,-1-delta)
    qr=q[1]*rho
    time=(q[0]+dz(q[0])*(z*((1-delta)/2))+qr)/L
    radial=q[0]*(q[0]/2+qr)
    transport=w[0]*(-q[0]*(z*2)+dz(q[0])*d-qr*(z*2))/L
    pressure=-(time+radial+transport)/2+algebra.shift(q[1]*2+q[2]*rho,(0,0,0,1))
    return dict(N1theta=truncate(theta,1),N1z=truncate(axial,1),N1p=truncate(pressure,1))


def matrix(algebra,rows,forcing,z,rho,delta):
    """Same six regular paper rows, now containing actual bridge sectors."""
    c=algebra.ctx;x=c.sqrt(rho);d=1-z*z;L=1-z*z*delta
    f,fr=rows['F'][:2];w,wr=rows['Uz'][:2];q=rows['Q'][0]
    H=z*((1-delta)/2)+d*w;B=q+(1-z*w*2)/L;E=z*w-c.mpf('.5')
    Cf=f+fr*rho
    Zf=(f*(z*(-2-delta))+dz(f)*d-fr*(z*(2*rho)))/L
    Zw=(w*(z*(-1-delta))+dz(w)*d-wr*(z*(2*rho)))/L
    eps=lambda value:algebra.shift(value,(0,0,0,-1))
    B0={};B1={}
    for i,j,sgn in ((1,5,1),(2,6,1),(3,6,-1)):B0[f'{i},{j}']=algebra.lift(sgn)
    B0['4,1']=eps(f*(4*x))
    for i,j,value in ((5,1,(q+E*(-2+delta)/L)*2),
        (5,2,Zf*2-Cf*(z*(2*(-1+delta)))/L),
        (5,3,-Cf*(z*(2*(1+delta)))/L),
        (6,1,eps(-f*(z*(8*rho))/L)),
        (6,2,(E*(-1+delta)/L+Zw-wr*(z*((-1+delta)*rho))/L)*2),
        (6,3,-wr*(z*(2*(1+delta)*rho))/L),(6,4,z*(-4)/L)):
        B0[f'{i},{j}']=eps(value)
    B0['5,5']=eps(B*x);B0['6,6']=eps(B*x)
    for i,j,value in ((5,1,H*2/L),(5,2,-Cf*d*2/L),(5,3,-Cf*d*2/L),
        (6,2,(H-wr*d*rho)*2/L),(6,3,-wr*d*(2*rho)/L),(6,4,d*2/L)):
        B1[f'{i},{j}']=eps(value)
    g={'4':eps(forcing['N1p']*(2*x)),
       '5':eps(forcing['N1theta']*2),
       '6':eps(forcing['N1z']*2)-eps(eps(forcing['N1p']*(z*(4*rho))/L))}
    return dict(B0={k:truncate(v,2) for k,v in B0.items()},
                B1={k:truncate(v,2) for k,v in B1.items()},g=g,D=(0,0,2,0,3,1))


def recover_Q(w,mean,z,delta):
    return (w*(z*2)-mean*(z*(1-delta))-dz(mean)*(1-z*z))/(1-z*z*delta)


def first_rows(algebra,core,raw,inputs,rho,phase,delta):
    """Actual FV/own moments + exact core comparison ODE, through total3."""
    c=algebra.ctx;lift=algebra.lift
    phi=lift(Jet(c,raw['actual_phi_axial5']))
    w0=lift(Jet(c,raw['actual_raw_V_axial5']))
    m0=lift(Jet(c,raw['actual_own_six_moments_axial5']['M']))
    pc=[Jet(c,[core['Phi'][domain.parent.parent.core_module.gridkey(i,k)]/math.factorial(k)
               for k in range(4-i)]) for i in range(3)]
    wc=[Jet(c,[core['Uz'][domain.parent.parent.core_module.gridkey(i,k)]/math.factorial(k)
               for k in range(4-i)]) for i in range(3)]
    A=pc[1]/pc[0];Ar=pc[2]/pc[0]-A*A
    sigma=domain.sigma_jets(c,phase)
    chi=lift(1-sigma[0])+algebra.width(sigma[0])
    chir=lift(sigma[1]/rho)-algebra.width(sigma[1]/rho,-1)
    pr=phi*chi*A
    prr=phi*(chir*A+chi*chi*A*A+chi*Ar)
    ratio=phi/pc[0]
    wr=ratio*chi*wc[1]
    wrr=ratio*(chir*wc[1]+chi*(chi-1)*A*wc[1]+chi*wc[2])
    mr=(w0-m0)/rho;mrr=(wr-mr*2)/rho
    z=lift(Jet.variable(c,raw['Z'],3))
    phirows=[phi,pr,prr];wrows=[w0,wr,wrr];mrows=[m0,mr,mrr]
    ratioF=Jet(c,[inputs['F0_ratios'][k]/math.factorial(k) for k in range(6)])
    frows=[algebra.shift(value*ratioF,(0,1,0,0)) for value in phirows]
    qrows=[recover_Q(w,m,z,delta) for w,m in zip(wrows,mrows)]
    C=lift(Jet(c,raw['actual_own_six_moments_axial5']['C']))
    ratioF2=Jet(c,[inputs['F0_squared_ratios'][k]/math.factorial(k) for k in range(6)])
    pressure=algebra.shift(inputs['p0'],(0,0,1,0))+algebra.shift(C*ratioF2*rho,(0,2,0,-1))
    return dict(F=frows,Uz=wrows,Q=qrows,mean=mrows,Phi=phirows,P0=pressure),z


def core_rows(algebra,profile,Z):
    c=algebra.ctx;g=profile['ordinary_mixed_profile_grids'];key=domain.parent.parent.core_module.gridkey
    def jets(label):
        return [algebra.lift(Jet(c,[g[label][key(i,k)]/math.factorial(k) for k in range(4-i)])) for i in range(3)]
    f=[algebra.shift(v,(0,1,0,0)) for v in jets(domain.parent.parent.core_module.F)]
    w=jets(domain.parent.parent.core_module.V);q=jets(domain.parent.parent.core_module.Q)
    mean=[algebra.lift(Jet(c,[profile['radial_average_Uz_derivatives'][key(i,k)]/math.factorial(k)
                            for k in range(4-i)])) for i in range(3)]
    pd=algebra.lift(Jet(c,[g[domain.parent.parent.core_module.PD][key(0,k)]/math.factorial(k) for k in range(4)]))
    pi=algebra.lift(Jet(c,[g[domain.parent.parent.core_module.PI][key(0,k)]/math.factorial(k) for k in range(4)]))
    return dict(F=f,Uz=w,Q=q,mean=mean,P0=algebra.shift(pd,(0,0,1,0))+algebra.shift(pi,(0,2,0,-1))),algebra.lift(Jet.variable(c,Z,3))


@dataclass(frozen=True,eq=False)
class OriginalN1InnerOperatorPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1InnerOperator:
    @source_precision
    def __init__(self,bounds,require_checked=True):
        if type(bounds) is not domain.CurrentOriginalN1InnerAnalyticDomain or not bounds.acceptance_loaded:
            raise ValueError('Accepted actual first-collar analytic bound owner required')
        bounds.assert_graph();self.bounds=bounds;self.source=bounds.source;self.n1=bounds.n1
        self.c=bounds.c;self.family=copy.deepcopy(bounds.family);self.hashes=dict(bounds.hashes)
        for name in (domain.NAME,domain.RECEIPT,Path(__file__).name,
                     PREFIX+'microswitch_mixed_C4.py'):self.hashes[name]=sha(name)
        self._owner=bounds;self._source=self.source;self._ctx=self.c;self._family=copy.deepcopy(self.family)
        self.frames={};self._frame_snapshots={};self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Actual inner operator receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Actual n1 operator input changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.bounds is not self._owner or self.source is not self._source or self.c is not self._ctx
                or self.source is not self.bounds.source or self.family!=self._family or self.family!=self.bounds.family
                or self.c is not self.bounds.c or not self.bounds.acceptance_loaded or self.hashes!=self._hashes):
            raise ValueError('Actual n1 operator source graph changed')
        return self.bounds.assert_graph()

    def frame(self,Z):
        z=self.c.mpf(Z)
        if ep(z)[0]<-1 or ep(z)[1]>1:raise ValueError('Original real Z domain required')
        key=z._mpi_
        if key not in self.frames:
            amplitude=self.source.amplitude.evaluate(z)
            logs=(self.source.logh,amplitude['logF0'],2*self.source.core.logP,self.source.core.logLambda)
            self.frames[key]=(z,N1LogAlgebra(self.c,logs,[]),self.source.inputs(z))
            self._frame_snapshots[key]=snapshot((z,logs,self.frames[key][2]))
        frame=self.frames[key]
        if frame[1].ctx is not self.c or snapshot((frame[0],frame[1].logs,frame[2]))!=self._frame_snapshots[key]:
            raise ValueError('Original n1 amplitude frame changed')
        return frame

    @source_precision
    def evaluate(self,Z,coordinate,chart='core'):
        self.assert_graph();c=self.c;z,algebra,inputs=self.frame(Z);value=c.mpf(coordinate)
        if chart=='core':
            if ep(value)[0]<0 or ep(value)[1]>4:raise ValueError('Actual core rho must lie in[0,4]')
            packet=self.n1.source.evaluate(value,z);profile=self.n1.source._fields[id(packet)][1]
            rows,zjet=core_rows(algebra,profile,z);rho=value;raw=None
            geometry=dict(chart='core',rho=rho,R='rho/Lambda',actual_source='canonical admitted analytic core')
        elif chart=='first':
            if ep(value)[0]<0 or ep(value)[1]>ep(c.mpf('.5'))[1]:raise ValueError('Actual first-collar phase in[0,1/2] required')
            raw=self.source.actual.packet(z,value,'first')
            radius=raw['radius_enclosure_only']*self.source.core.Lambda
            # Exact rho=4*exp(hb*s)>=4. Intersect an outward radius
            # enclosure; no rounded endpoint is selected as source radius.
            rho=(c.mpf(4) if ep(value)[1]==0 and raw['exact_core_inlet'] else
                 c.mpf([max(c.mpf(4).a, radius.a),min(c.mpf('4.1').b,radius.b)]))
            core=self.n1.source.core.normalized_jets(rho,z)
            rows,zjet=first_rows(algebra,core,raw,inputs,rho,value,self.source.core.delta)
            geometry=dict(chart='first',phase=value,rho_enclosure_only=rho,
                exact_rho='4*exp(hb*s)',R='Ra*exp(hb*s)',actual_source='original prescribed-shear FV and own cumulative moments',
                exact_width_not_materialized=True,comparison_only_used_for_known_core_direction=True)
        else:raise ValueError('Only actual core and first collar supported')
        forcing=known_forcing(algebra,rows,zjet,rho,self.source.core.delta)
        system=matrix(algebra,rows,forcing,zjet,rho,self.source.core.delta)
        result=dict(source_family=self.family,Z=z,geometry=geometry,
            source_log_basis_names=['hb','F0_anchor','Pstar_squared','Lambda'],source_log_bases=list(algebra.logs),
            ordinary_leading_rho_rows={name:[export(v) for v in rows[name]] for name in ('F','Uz','Q','mean')},
            independent_actual_pressure=export(rows['P0']),
            genuine_n1_known_forcing={name:export(v,1) for name,v in forcing.items()},
            regular_system=dict(D=list(system['D']),B0={key:export(v,2) for key,v in system['B0'].items()},
                B1={key:export(v,2) for key,v in system['B1'].items()},g={key:export(v,1) for key,v in system['g'].items()},
                equation='d_x W + D W/x = B0 W + B1 d_Z W + g',x='sqrt(rho)',
                unknowns=['F1','Uz1','K1','P1','d_x_F1','d_x_Uz1']),
            amplitude_anchor_fixed_during_jet_algebra=True,ordinary_Z_amplitude_derivatives_in_coefficients=True,
            source_width_and_Lambda_exponents_retained=True,actual_own_mean_not_comparison_mean=True,
            outputs_are_enclosures_of_original_functions_not_selected_point_values=True,
            actual_point_moment_history_recovered=False,
            radial_mean_rows_use_exact_FTC_of_original_own_mean=True,
            no_positive_order_solution_in_this_packet=True,**dict.fromkeys(OPEN,False))
        packet=OriginalN1InnerOperatorPacket()
        self._packets[id(packet)]=(packet,result,rows,forcing,system,algebra,
            snapshot(result),snapshot((rows,forcing,system)))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1InnerOperatorPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live actual n1 operator packet issued by this owner required')
        if snapshot(entry[1])!=entry[6] or snapshot((entry[2],entry[3],entry[4]))!=entry[7]:
            raise ValueError('Actual n1 operator packet changed')
        return dict(**copy.deepcopy(entry[1]),**{GATE:self.acceptance_loaded})


@source_precision
def run(bounds):
    began=time.monotonic();owner=CurrentOriginalN1InnerOperator(bounds,require_checked=False)
    packets={name:owner.evaluate(Z,coordinate,chart) for name,chart,coordinate,Z in SAMPLES}
    raw=dict(source_family=owner.family,samples={name:owner.report(v) for name,v in packets.items()},input_hashes=owner.hashes,
             **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    data=json.dumps(domain.parent.parent.first.encode(domain.parent.parent.original.serialized(raw)),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    print('ACTUAL_N1_INNER_OPERATOR_READY',len(packets),flush=True)
    return owner,packets

"""Genuine n=1 regular system and its coupled finite axis solution.

The hierarchy index is temporal n=1. Radial Taylor degree is a separate
index. Finite jets solve the coefficient equations at the axis, not the
full-interval PDE. No missing complex/tail estimate is replaced by a fit.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_core_hierarchy_source as parent
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as Jet
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=parent.HERE,parent.PREFIX,parent.sha
NAME=PREFIX+'current_original_n1_regular_system.json.gz'
RECEIPT=PREFIX+'current_original_n1_regular_system_check.json'
GATE='current_original_n1_regular_system_and_coupled_axis_jets_installed'
OPEN=parent.OPEN+('n1_full_inner_interval_solution_certified','n1_radial_Taylor_remainder_certified')
AXIS_SAMPLES={'whole_axis':['-1','1'],'offcenter':'.371','center':'0'}
SYSTEM_SAMPLES={'whole_core':(['0','4.1'],['-1','1']),
                'axis':('0','.371'),'inner':('2','.371'),'beyond_Ra':('4.1','.371')}
D=(0,0,2,0,3,1)


def derivative(jet):
    if jet.order<1:raise ValueError('An actual further source derivative is required')
    return Jet(jet.ctx,[(i+1)*jet[i+1] for i in range(jet.order)])


class Sectors:
    """Finite jets multiplied by formal A_anchor**p, never exp(logF0).

    A_anchor is constant during jet algebra. The original normalized F
    derivative grid already includes derivatives of the actual amplitude.
    Treating those rows as derivatives of Phi alone would be incorrect.
    """
    def __init__(self,c,rows,order=None):
        self.c=c
        self.order=min(v.order for v in rows.values()) if order is None else order
        self.rows={p:v.truncate(self.order) for p,v in rows.items()
                   if any(parent.ends(x)!=(0,0) for x in v.coefficients)}

    @classmethod
    def zero(cls,c,order):return cls(c,{},order)

    @classmethod
    def one(cls,c,value,order,p=0):
        jet=value if isinstance(value,Jet) else Jet.constant(c,value,order)
        return cls(c,{p:jet},order)

    def truncate(self,order):return Sectors(self.c,self.rows,min(order,self.order))

    def _pair(self,other):
        if not isinstance(other,Sectors):other=Sectors.one(self.c,other,self.order)
        if other.c is not self.c:raise ValueError('The original interval contexts must match')
        return self,other,min(self.order,other.order)

    def __add__(self,other):
        left,right,order=self._pair(other);zero=Jet.constant(self.c,0,order)
        return Sectors(self.c,{p:left.rows.get(p,zero)+right.rows.get(p,zero)
                               for p in left.rows.keys()|right.rows.keys()},order)

    __radd__=__add__
    def __neg__(self):return Sectors(self.c,{p:-v for p,v in self.rows.items()},self.order)
    def __sub__(self,other):return self+-other
    def __mul__(self,other):
        left,right,order=self._pair(other);rows={};zero=Jet.constant(self.c,0,order)
        for p,a in left.rows.items():
            for q,b in right.rows.items():rows[p+q]=rows.get(p+q,zero)+a*b
        return Sectors(self.c,rows,order)

    __rmul__=__mul__
    def __truediv__(self,other):
        if isinstance(other,Sectors):raise ValueError('Formal amplitude sectors may not be inverted')
        return self*(other.reciprocal() if isinstance(other,Jet) else self.c.mpf(1)/other)

    def dz(self):
        if self.order<1:raise ValueError('No unsupported Z derivative may be padded with zero')
        return Sectors(self.c,{p:derivative(v) for p,v in self.rows.items()},self.order-1)


def leading_radial_jets(c,grid):
    return [Jet(c,[grid[parent.core_module.gridkey(i,k)]/(math.factorial(i)*math.factorial(k))
                   for k in range(5-i)]) for i in range(5)]


def za_component(value,z,d,L,weight,radial_degree):
    return (z*value*(weight-2*radial_degree)+d*derivative(value))/L


def known_source_jets(c,grid,zvalue,delta,Lambda):
    """Rho Taylor coefficients of the genuine sources, total degree <=2."""
    z=Jet.variable(c,zvalue,4);d=1-z*z;L=1-z*z*delta
    f,w,v=(leading_radial_jets(c,grid[key]) for key in
           (parent.core_module.F,parent.core_module.V,parent.core_module.Q))
    theta=[];axial=[];pressure=[]
    for i in range(3):
        def second(rows,weight):
            first=za_component(rows[i],z,d,L,weight,i)
            return -za_component(first,z,d,L,weight-1+delta,i).truncate(2-i)
        theta.append(second(f,-2-delta));axial.append(second(w,-1-delta))
        omega=((1+i)*v[i]+z*derivative(v[i])*((1-delta)/2))/L
        for j in range(i+1):
            k=i-j
            omega+=v[j]*v[k]*(c.mpf(1)/2+k)
            omega+=w[j]*za_component(v[k],z,d,L,-2,k)
        pressure.append((-omega/2+v[i+1]*(Lambda*(i+1)*(i+2))).truncate(2-i))
    return f,w,v,theta,axial,pressure


def axis_solve(c,grid,zvalue,delta,Lambda,epsilon):
    """Solve n=1 axis equations to rho^3; first two rows include d_Z.

    At degree k the ordinary Z jets have order 3-k. The recovered v1=V1/R
    is available to rho^2. Its degree-three Z derivative is not present in
    the leading total-four source and is never assumed to be zero.
    """
    f,w,v,nt,nz,np=known_source_jets(c,grid,zvalue,delta,Lambda)
    z=Jet.variable(c,zvalue,4);d=1-z*z;L=1-z*z*delta
    b0,c0,b1,c1,p1,q1=-2-delta,-1-delta,-2+delta,-1+delta,-2,1+delta
    B=[v[i]+((1-2*z*w[i])/L if i==0 else -2*z*w[i]/L) for i in range(3)]
    H=[(z*((1-delta)/2)+d*w[i])/L if i==0 else d*w[i]/L for i in range(3)]
    E=[(z*w[i]-(c.mpf(1)/2 if i==0 else 0))/L for i in range(3)]
    Zf=[za_component(f[i],z,d,L,b0,i) for i in range(3)]
    Zw=[za_component(w[i],z,d,L,c0,i) for i in range(3)]
    F={0:Sectors.zero(c,3)};U={0:Sectors.zero(c,3)};P={0:Sectors.zero(c,3)}
    K={0:Sectors.zero(c,3)};Q={0:Sectors.zero(c,2)}
    for k in range(1,4):
        m=k-1;order=3-k
        rf=Sectors.one(c,nt[m],order,1);ru=Sectors.one(c,nz[m],order)
        rp=Sectors.one(c,np[m],order)
        for j in range(1,m+1):
            i=m-j
            rf+=F[j]*(j*B[i])+F[j].dz()*H[i]+F[j]*(v[i]+E[i]*b1)
            rf+=Sectors.one(c,(i+1)*f[i],order,1)*Q[j]
            rf+=Sectors.one(c,Zf[i],order,1)*U[j]
            ru+=U[j]*(j*B[i])+U[j].dz()*H[i]+U[j]*(E[i]*c1)+U[j]*Zw[i]
            if i>0:ru+=Q[j]*(i*w[i])
            rp+=2*Sectors.one(c,f[i],order,1)*F[j]
        if m>0:ru+=(P[m]*((p1-2*m)*z)+P[m].dz()*d)*(1/L)
        F[k]=(rf*(epsilon/(2*k*(k+1)))).truncate(order)
        U[k]=(ru*(epsilon/(2*k*k))).truncate(order)
        P[k]=(rp*(epsilon/k)).truncate(order)
        K[k]=U[k]*(-c.mpf(k)/(k+1))
        if k<3:
            Q[k]=(U[k]*(z*(2-q1/c.mpf(k+1)))-U[k].dz()*d/c.mpf(k+1))*(1/L)
    return dict(F1=F,Uz1=U,P1=P,K1=K,V1_over_R=Q),dict(N1theta=nt,N1z=nz,N1p=np)


def export_axis(c,rows,logF,logLambda):
    result={}
    for field,coefficients in rows.items():
        result[field]={}
        for k,value in coefficients.items():
            result[field][str(k)]={str(p):dict(
                rho_Taylor_coefficients=list(jet.coefficients),
                ordinary_Z_derivatives=[math.factorial(j)*v for j,v in enumerate(jet.coefficients)],
                log_amplitude_scale=p*logF,
                ordinary_R_derivative_at_axis=dict(coefficient=math.factorial(k)*jet[0],
                    log_scale=p*logF+k*logLambda),amplitude_power=p)
                for p,jet in value.rows.items()}
    return result


def regular_system(c,profile,forcing,delta,epsilon):
    """Paper (14.6)-(14.8), exactly transformed to x=sqrt(rho)."""
    g=profile['ordinary_mixed_profile_grids'];rho,z=profile['rho'],profile['Z'];x=c.sqrt(rho)
    at=lambda label,i,k:g[label][parent.core_module.gridkey(i,k)]
    f,w,v=at(parent.core_module.F,0,0),at(parent.core_module.V,0,0),at(parent.core_module.Q,0,0)
    fr,wr=at(parent.core_module.F,1,0),at(parent.core_module.V,1,0)
    fz,wz=at(parent.core_module.F,0,1),at(parent.core_module.V,0,1)
    L=1-delta*z*z;d=1-z*z;H=(1-delta)*z/2+d*w;B=v+(1-2*z*w)/L;E=z*w-c.mpf(1)/2
    b0,c0,b1,c1,p1,q1=-2-delta,-1-delta,-2+delta,-1+delta,-2,1+delta
    Cf=f+rho*fr;Zf=(b0*z*f+d*fz-2*z*rho*fr)/L;Zw=(c0*z*w+d*wz-2*z*rho*wr)/L
    A0={};A1={};rhs={}
    def put(matrix,i,j,p,value):matrix[str(i)+','+str(j)]={str(p):value}
    for i,j,sgn in ((1,5,1),(2,6,1),(3,6,-1)):put(A0,i,j,0,c.mpf(sgn))
    put(A0,4,1,1,4*epsilon*x*f)
    for i,j,p,value in ((5,1,0,2*(v+b1*E/L)),(5,2,1,2*Zf-2*c1*z*Cf/L),
          (5,3,1,-2*q1*z*Cf/L),(6,1,1,-8*epsilon*z*rho*f/L),
          (6,2,0,2*(c1*E/L+Zw-c1*z*rho*wr/L)),
          (6,3,0,-2*q1*z*rho*wr/L),(6,4,0,2*p1*z/L)):
        put(A0,i,j,p,epsilon*value)
    put(A0,5,5,0,epsilon*x*B);put(A0,6,6,0,epsilon*x*B)
    for i,j,p,value in ((5,1,0,2*H/L),(5,2,1,-2*d*Cf/L),(5,3,1,-2*d*Cf/L),
          (6,2,0,2*(H-d*rho*wr)/L),(6,3,0,-2*d*rho*wr/L),(6,4,0,2*d/L)):
        put(A1,i,j,p,epsilon*value)
    ntheta,nz,np=(forcing[key]['coefficient'] for key in ('N1theta','N1z','N1p'))
    rhs['4']={'0':2*epsilon*x*np};rhs['5']={'1':2*epsilon*ntheta}
    rhs['6']={'0':2*epsilon*nz-4*epsilon**2*z*rho*np/L}
    return dict(coordinate='x=sqrt(rho), rho=Lambda*R',unknowns=['F1','Uz1','K1','P1','d_x_F1','d_x_Uz1'],
        equation='d_x W + D W/x = B0 W + B1 d_Z W + g',D=list(D),B0=A0,B1=A1,g=rhs,
        recovery_q1='1+delta',recovery_K1='Mz1/R-Uz1',
        regular_inverse='G_x g_i=x*integral_0^1 t^D_i*g_i(tx,Z) dt',
        B1_maps_first_four_to_last_two_and_annihilates_last_two=True,
        separated_derivative_blocks_are_not_discarded=True,
        singular_axis_division_not_evaluated=True,pressure_2F0F1_retained=True)


@dataclass(frozen=True,eq=False)
class OriginalN1RegularPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1RegularSystem:
    @source_precision
    def __init__(self,source,require_checked=True):
        if type(source) is not parent.CurrentOriginalCoreHierarchySource or not source.acceptance_loaded:
            raise ValueError('Accepted live genuine original hierarchy source required')
        self.source=source;self.c=source.ctx;self.family=copy.deepcopy(source.family)
        self._source=source;self._context=source.ctx;self._family=copy.deepcopy(source.family)
        self.hashes=dict(source.hashes)
        for name in (parent.NAME,parent.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Genuine n1 finite-axis/system receipt or scope differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('N1 input changed '+name)
                parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hash_snapshot=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.source is not self._source or self.c is not self._context or self.c is not self.source.ctx
                or self.family!=self._family or self.family!=self.source.family or not self.source.acceptance_loaded
                or self.hashes!=self._hash_snapshot):
            raise ValueError('N1 owner source, interval context or admitted family changed')
        return self.source.assert_graph()

    def _issue(self,result):
        value=OriginalN1RegularPacket()
        self._packets[id(value)]=(value,result,copy.deepcopy(parent.first.encode(parent.original.serialized(result))))
        return value

    @source_precision
    def axis(self,Z):
        self.assert_graph()
        packet=self.source.evaluate('0',Z);view=self.source.report(packet)
        profile=self.source._fields[id(packet)][1];c=self.c;core=self.source.core
        rows,forcing=axis_solve(c,profile['ordinary_mixed_profile_grids'],profile['Z'],
                              core.delta,core.Lambda,core.epsilon)
        logF=view['genuine_order_one_known_forcing']['N1theta']['log_scale']
        return self._issue(dict(kind='finite_genuine_n1_axis_solution',Z=profile['Z'],
            temporal_hierarchy_order=1,radial_Taylor_coordinate='rho=Lambda*R',
            axis_coefficients=export_axis(c,rows,logF,core.logLambda),
            known_source_rho_Taylor_coefficients={key:[list(v.coefficients) for v in jets] for key,jets in forcing.items()},
            F1_Uz1_P1_K1_radial_degree=3,V1_over_R_radial_degree=2,
            Z_derivative_orders_by_radial_degree={'1':2,'2':1,'3':0},
            zero_axis_values=True,zero_x_derivative_axis_values=True,
            axis_values={key:'0' for key in ('F1','Uz1','K1','P1','V1','d_x_F1','d_x_Uz1')},
            original_pressure_coupling_in_P2_and_P3=True,
            original_pressure_coupling_backreaction_in_Uz1_degree3=True,
            amplitude_exp_not_materialized=True,positive_order_radial_jet_not_leading_radial_row=True,
            finite_axis_jets_only=True,whole_interval_polynomial_evaluation_forbidden=True,
            remaining_domain='Need common complex source bounds and full regular Volterra/PDE solve with controlled tail',
            **dict.fromkeys(OPEN,False)))

    @source_precision
    def system(self,rho,Z):
        self.assert_graph()
        packet=self.source.evaluate(rho,Z);view=self.source.report(packet)
        profile=self.source._fields[id(packet)][1];core=self.source.core
        result=regular_system(self.c,profile,view['genuine_order_one_known_forcing'],core.delta,core.epsilon)
        return self._issue(dict(kind='genuine_n1_full_regular_operator',rho=profile['rho'],Z=profile['Z'],
            amplitude_log=view['genuine_order_one_known_forcing']['N1theta']['log_scale'],
            formal_amplitude_scale='Every sector p means exp(p*logF0); never materialized',
            regular_system=result,**dict.fromkeys(OPEN,False)))

    @source_precision
    def report(self,value):
        entry=self._packets.get(id(value))
        if type(value) is not OriginalN1RegularPacket or entry is None or entry[0] is not value:
            raise ValueError('Live n1 packet issued by this owner required')
        self.assert_graph()
        if parent.first.encode(parent.original.serialized(entry[1]))!=entry[2]:raise ValueError('N1 packet changed')
        return dict(source_family=copy.deepcopy(self.family),**copy.deepcopy(entry[1]),**{GATE:self.acceptance_loaded})


@source_precision
def run(source):
    began=time.monotonic();owner=CurrentOriginalN1RegularSystem(source,require_checked=False)
    axes={name:owner.axis(z) for name,z in AXIS_SAMPLES.items()}
    systems={name:owner.system(r,z) for name,(r,z) in SYSTEM_SAMPLES.items()}
    raw=dict(source_family=owner.family,axis_samples={k:owner.report(v) for k,v in axes.items()},
        system_samples={k:owner.report(v) for k,v in systems.items()},input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    data=json.dumps(parent.first.encode(parent.original.serialized(raw)),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    print('GENUINE_N1_REGULAR_SYSTEM_AND_AXIS_READY',len(axes),len(systems),flush=True)
    return owner,axes,systems

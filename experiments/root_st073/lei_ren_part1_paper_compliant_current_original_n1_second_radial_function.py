"""Original leading rho3 sources and n1 second radial mixed Z1 jets.

Derivatives follow the actual prescribed-shear collar and six equations.
Axis limits are separate regular coefficients, never divisions by zero.
Finite limits do not certify a radial Taylor remainder or own repair.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_regular_series_function as parent
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=parent.op;natural=parent.natural;mean=natural.mean;weighted=parent.weighted
domain=parent.domain;ep=parent.ep;HERE,PREFIX,sha=parent.HERE,parent.PREFIX,parent.sha
NAME=PREFIX+'current_original_n1_second_radial_function.json.gz'
RECEIPT=PREFIX+'current_original_n1_second_radial_function_check.json'
GATE='current_original_n1_leading_rho3_sources_and_second_radial_mixed_Z1_installed'
OPEN=parent.OPEN


class RadialPair:
    """A value and its genuine ordinary radial derivative, both Z jets."""
    def __init__(self,value,radial):
        if value.algebra is not radial.algebra:raise ValueError('Same original amplitude algebra required')
        self.v=value;self.r=radial;self.algebra=value.algebra;self.ctx=value.ctx
        self.order=min(value.order,radial.order)

    def _pair(self,other):
        if isinstance(other,RadialPair):
            if other.algebra is not self.algebra:raise ValueError('Foreign radial pair')
            return other
        value=other if isinstance(other,op.FactoredJet) else self.algebra.lift(other,self.order)
        return RadialPair(value,self.algebra.lift(0,value.order))

    def __add__(self,other):
        other=self._pair(other);return RadialPair(self.v+other.v,self.r+other.r)
    __radd__=__add__
    def __neg__(self):return RadialPair(-self.v,-self.r)
    def __sub__(self,other):return self+-self._pair(other)
    def __rsub__(self,other):return self._pair(other)+-self
    def __mul__(self,other):
        other=self._pair(other);return RadialPair(self.v*other.v,self.r*other.v+self.v*other.r)
    __rmul__=__mul__
    def __truediv__(self,other):
        other=self._pair(other)
        return RadialPair(self.v/other.v,(self.r*other.v-self.v*other.r)/(other.v*other.v))
    def Zderivative(self):return RadialPair(op.dz(self.v),op.dz(self.r))


def shifted(pair,powers):
    return RadialPair(pair.algebra.shift(pair.v,powers),pair.algebra.shift(pair.r,powers))


def core_third_rows(algebra,physical,rho,Z):
    """Reuse the original total-order-five model and nonlinear bounds."""
    c=algebra.ctx;raw=physical.normalized_jets(rho,Z);key=domain.parent.parent.core_module.gridkey
    err=algebra.lift(op.Jet(c,[c.mpf([-1,1])*raw['nonlinear_correction_bounds'][key(3,k)]/math.factorial(k)
                              for k in range(3)]))
    phi=algebra.lift(op.Jet(c,[raw['Phi'][key(3,k)]/math.factorial(k) for k in range(3)]))
    return dict(Psi=err,meanPsi=err/4,Phi=phi)


def collar_third_rows(algebra,physical,rows,rho,phase,Z):
    """Differentiate the same original Phi/Psi prescribed-shear ODE."""
    c=algebra.ctx;extension=core_third_rows(algebra,physical,rho,Z)
    if ep(phase)[1]==0:return extension
    raw=physical.normalized_jets(rho,Z);key=domain.parent.parent.core_module.gridkey
    pc=[op.Jet(c,[raw['Phi'][key(i,k)]/math.factorial(k) for k in range(6-i)]) for i in range(4)]
    A=pc[1]/pc[0];Ar=pc[2]/pc[0]-A*A
    Arr=pc[3]/pc[0]-3*pc[1]*pc[2]/(pc[0]*pc[0])+2*A*A*A
    sigma=domain.sigma_jets(c,phase)
    chi=algebra.lift(1-sigma[0],4)+algebra.width(algebra.lift(sigma[0],4))
    chir=algebra.lift(sigma[1]/rho,3)-algebra.width(algebra.lift(sigma[1]/rho,3),-1)
    s1=algebra.lift(sigma[1]/rho**2,2);s2=algebra.lift(2*sigma[2]/rho**2,2)
    chirr=algebra.width(s1,-1)-s1-algebra.width(s2,-2)+algebra.width(s2,-1)
    C=chi*A;Cr=chir*A+chi*Ar;Crr=chirr*A+2*chir*Ar+chi*Arr
    phi3=op.truncate(rows['Phi'][0]*(C*C*C+3*C*Cr+Crr),2)
    core_rows=mean.core_rows(algebra,physical,rho,Z);E1,E2=core_rows['Psi'][1:3];E3=extension['Psi']
    ratio=rows['Phi'][0]/pc[0]
    psi3=op.truncate(ratio*(E1*(chirr+chir*(3*chi-2)*A+chi*(chi-1)*Ar+chi*(chi-1)*(chi-1)*A*A)
                               +E2*(2*chir+2*chi*(chi-1)*A)+chi*E3),2)
    mean3=op.truncate((rows['Psi'][2]-3*rows['meanPsi'][2])/rho,2)
    return dict(Psi=psi3,meanPsi=mean3,Phi=phi3,chi_rhorho=op.truncate(chirr,2))


def extend_leading_rows(owner,mean_packet):
    mv=owner.mean.report(mean_packet);entry=owner.mean._packets[id(mean_packet)]
    lower,fields=entry[4][1:3];algebra=fields['Uz'][0].algebra
    rho=mv['geometry'].get('rho',mv['geometry'].get('rho_enclosure_only'));Z=mv['Z']
    if mv['query']['chart']=='core':third=core_third_rows(algebra,owner.physical,rho,Z)
    else:
        phase=mv['query']['actual_source_coordinate']
        third=collar_third_rows(algebra,owner.physical,lower,rho,phase,Z)
    rows={k:[*lower[k],third[k]] for k in ('Psi','meanPsi','Phi')}
    z=algebra.lift(lower['axis']['z']);d=algebra.lift(lower['axis']['d']);L=algebra.lift(lower['axis']['L'])
    qpsi3=(rows['Psi'][3]*(2*z)-rows['meanPsi'][3]*(z*(1-owner.physical.delta))-op.dz(rows['meanPsi'][3])*d)/L
    fields={**fields,'Uz':[*fields['Uz'],algebra.shift(rows['Psi'][3],(0,0,0,-1))],
            'mean':[*fields['mean'],algebra.shift(rows['meanPsi'][3],(0,0,0,-1))],
            'Q':[*fields['Q'],algebra.shift(qpsi3,(0,0,0,-1))]}
    return rows,fields,third,mv


def source_rho_pairs(algebra,rows,fields,zvalue,rho,delta,gradient):
    c=algebra.ctx;lift=algebra.lift
    constant=lambda v,n=5:RadialPair(v if isinstance(v,op.FactoredJet) else lift(v,n),lift(0,n))
    z=constant(lift(op.Jet.variable(c,zvalue,5)));r=RadialPair(lift(rho,5),lift(1,5))
    d=1-z*z;L=1-z*z*delta
    pairs=lambda row:[RadialPair(row[i],row[i+1]) for i in range(3)]
    phi=pairs(rows['Phi']);w=pairs(fields['Uz']);q=pairs(fields['Q'])
    axial=-op.axial_second(w,z,r,delta,-1-delta)
    qr=q[1]*r
    time_term=(q[0]+op.dz(q[0])*(z*((1-delta)/2))+qr)/L
    radial=q[0]*(q[0]/2+qr)
    transport=w[0]*(-q[0]*(z*2)+op.dz(q[0])*d-qr*(z*2))/L
    Np=-(time_term+radial+transport)/2+shifted(q[1]*2+q[2]*r,(0,0,0,1))
    f,fr,_=phi;grad=constant(lift(gradient));m=grad*d/L;a=-2-delta;b=-3
    cross=m*(f*(z*(a+b))+op.dz(f)*(2*d)-fr*(z*(4*r)))/L+op.dz(m)*d*f/L
    g5=-2*m*m*f+shifted(2*cross,(0,0,0,-1))-shifted(2*op.axial_second(phi,z,r,delta,a),(0,0,0,-2))
    g6=2*axial-shifted(Np*z*(4*r)/L,(0,0,0,-1))
    return {'4_over_x':shifted(Np*2,(0,0,-1,0)),'5':g5,'6':g6}


def matrix_x_pairs(algebra,rows,fields,zvalue,rho,delta,gradient,chi):
    """Differentiate true B0/B1 in x, including explicit x and rho=x^2."""
    c=algebra.ctx;x=c.sqrt(rho);lift=algebra.lift
    constant=lambda v,n=5:RadialPair(v if isinstance(v,op.FactoredJet) else lift(v,n),lift(0,n))
    pair=lambda v,vr:RadialPair(v,vr*(2*x))
    z=constant(lift(op.Jet.variable(c,zvalue,5)));r=RadialPair(lift(rho,5),lift(2*x,5));xx=RadialPair(lift(x,5),lift(1,5))
    d=1-z*z;L=1-z*z*delta;grad=constant(lift(gradient))
    f,fr=pair(rows['Phi'][0],rows['Phi'][1]),pair(rows['Phi'][1],rows['Phi'][2])
    uz,ur=pair(fields['Uz'][0],fields['Uz'][1]),pair(fields['Uz'][1],fields['Uz'][2])
    q=pair(fields['Q'][0],fields['Q'][1]);psi=pair(rows['Psi'][0],rows['Psi'][1])
    H=z*((1-delta)/2)+d*uz;E=z*uz-c.mpf('.5');B=q+(1-2*z*uz)/L;Cf=f+fr*r
    Ta=(f*(z*(-2-delta))+op.dz(f)*d-fr*(z*(2*r)))/L
    Zw=(uz*(z*(-1-delta))+op.dz(uz)*d-ur*(z*(2*r)))/L;m=grad*d/L
    eps=lambda v,p=1:shifted(v,(0,0,0,-p))
    amplitude=constant(weighted.relative_amplitude(algebra,gradient,2,2))
    f2=lambda v:shifted(v*amplitude,(0,2,0,0))
    B0={'4,1':f2(f*(4*xx)),
        '5,1':eps(2*(q+E*(-2+delta)/L))+eps(-2*psi*d*grad/L),
        '5,2':eps(2*Ta-Cf*(z*(2*(-1+delta)))/L,2)-eps(2*m*f),
        '5,3':eps(-Cf*(z*(2*(1+delta)))/L,2),'5,5':eps(B*xx),
        '6,1':eps(f2(-8*z*r*f/L)),
        '6,2':eps(2*(E*(-1+delta)/L+Zw-ur*(z*(r*(-1+delta)))/L)),
        '6,3':eps(-ur*(z*(r*(2*(1+delta))))/L),'6,4':eps(-4*z/L),'6,6':eps(B*xx)}
    B1={'5,1':eps(2*H/L),'5,2':eps(-2*d*Cf/L,2),'5,3':eps(-2*d*Cf/L,2),
        '6,2':eps(2*(H-d*r*ur)/L),'6,3':eps(-2*d*r*ur/L),'6,4':eps(2*d/L)}
    matrices={}
    for label,entries in (('B0',B0),('B1',B1)):
        matrices[label]={}
        for key,value in entries.items():
            i,j=(int(v)-1 for v in key.split(','))
            matrices[label][key]=shifted(value,tuple(b-a for a,b in zip(natural.EXTRA[i],natural.EXTRA[j])))
    for key,sgn in (('1,5',1),('2,6',1),('3,6',-1)):matrices['B0'][key]=constant(sgn,2)
    matrices['B0']['5,1']=matrices['B0']['5,1']+constant(lift(-2*chi))
    return matrices


def axis_second_rho(algebra,system,srho,mx):
    """Even rho2 limits from the x^2 coefficient of the genuine system."""
    a=[system['g']['5']/8,system['g']['6']/4,-system['g']['6']/8,system['g4_over_x']/2]
    rhs={5:srho['5'],6:srho['6']}
    for row in (5,6):
        for column in range(1,5):
            key=f'{row},{column}'
            if key in system['B0']:rhs[row]+=system['B0'][key]*a[column-1]
            if key in system['B1']:rhs[row]+=system['B1'][key]*op.dz(a[column-1])
        rhs[row]+=mx['B0'][f'{row},{row}']*(2*a[0 if row==5 else 1])
    return [op.truncate(rhs[5]/12,1),op.truncate(rhs[6]/8,1),op.truncate(-rhs[6]/12,1),
            op.truncate((srho['4_over_x']+mx['B0']['4,1']*a[0])/2,1)]


def pair_records(sources,matrices):
    export=lambda v:dict(value=op.export(v.v),radial=op.export(v.r))
    return dict(sources={k:export(v) for k,v in sources.items()},
                matrices={label:{k:export(v) for k,v in entries.items()} for label,entries in matrices.items()})


def own_second_recovery(drr,first_recovery_dr,Z,rho,delta):
    algebra=drr[0].algebra;c=algebra.ctx;z=algebra.lift(op.Jet.variable(c,Z,1))
    d=1-z*z;L=1-z*z*delta;U,K=drr[1:3];M=U+K
    Q=(U*(z*(1-delta))-K*(z*(1+delta))-op.dz(M)*d)/L
    return dict(mean1=M,Q1=Q,V1=op.truncate(Q*rho+2*first_recovery_dr['Q1'],0))


@dataclass(frozen=True,eq=False)
class OriginalN1SecondRadialPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1SecondRadialFunction:
    @source_precision
    def __init__(self,accepted_regular,require_checked=True):
        if type(accepted_regular) is not parent.CurrentOriginalN1RegularSeriesFunction or not accepted_regular.acceptance_loaded:
            raise ValueError('Accepted original common-inner regular function required')
        accepted_regular.assert_graph();self.regular=accepted_regular;self.mean=accepted_regular.mean
        self.c=accepted_regular.c;self.physical=accepted_regular.physical;self.family=copy.deepcopy(accepted_regular.family)
        self.hashes=dict(accepted_regular.hashes)
        for name in (parent.NAME,parent.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._parents=(accepted_regular,self.mean,self.physical,self.c);self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Second radial receipt or scope differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Second radial source changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if ((self.regular,self.mean,self.physical,self.c)!=self._parents or self.mean is not self.regular.mean or
            self.physical is not self.regular.physical or self.hashes!=self._hashes or self.family!=self.regular.family):
            raise ValueError('Second radial source graph changed')
        return self.regular.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate,chart='core'):
        self.assert_graph();began=time.monotonic();c=self.c
        p=self.regular.evaluate(Z,coordinate,chart);pv=self.regular.report(p);pe=self.regular._packets[id(p)]
        W,dx,dr,_,first_recovery_dr,system=pe[5];algebra=W[0].algebra
        mean_packet=self.regular.natural._packets[id(pe[2])][2]
        rows,fields,third,mv=extend_leading_rows(self,mean_packet)
        x,rho=system['x'],system['rho'];axis=ep(x)[1]==0;z=mv['Z'];gradient=system['gradient']
        sources=source_rho_pairs(algebra,rows,fields,z,rho,self.physical.delta,gradient)
        source_rho={k:op.truncate(v.r,1) for k,v in sources.items()}
        matrices=matrix_x_pairs(algebra,rows,fields,z,rho,self.physical.delta,gradient,system['chi'])
        mx={label:{k:op.truncate(v.r,1) for k,v in entries.items()} for label,entries in matrices.items()}
        gx={'4':op.truncate(sources['4_over_x'].v+sources['4_over_x'].r*(2*rho),1),
            '5':source_rho['5']*(2*x),'6':source_rho['6']*(2*x)}
        if axis:
            dxx=[system['g']['5']/4,system['g']['6']/2,-system['g']['6']/4,system['g4_over_x'],
                 algebra.lift(0,1),algebra.lift(0,1)]
            dxx=[op.truncate(v,1) for v in dxx]
            drr=axis_second_rho(algebra,system,source_rho,mx)
        else:
            dxx=[]
            for i in range(6):
                value=op.truncate(W[i],1)*(system['D'][i]/rho)-op.truncate(dx[i],1)*(system['D'][i]/x)
                value+=gx.get(str(i+1),algebra.lift(0,1))
                for key,B in system['B0'].items():
                    row,column=(int(v)-1 for v in key.split(','))
                    if row==i:value+=mx['B0'][key]*op.truncate(W[column],1)+op.truncate(B,1)*op.truncate(dx[column],1)
                for key,B in system['B1'].items():
                    row,column=(int(v)-1 for v in key.split(','))
                    if row==i:value+=mx['B1'][key]*op.truncate(op.dz(W[column]),1)+op.truncate(B,1)*op.truncate(op.dz(dx[column]),1)
                dxx.append(op.truncate(value,1))
            drr=[op.truncate((vv-op.truncate(vx,1)/x)/(4*rho),1) for vv,vx in zip(dxx,dx)]
        physical_xx={name:weighted.transform(v,gradient,scale,1) for name,v,scale in zip(parent.FIELDS,dxx,natural.WEIGHTS)}
        physical_rr={name:weighted.transform(v,gradient,scale,1) for name,v,scale in zip(parent.FIELDS,drr,natural.WEIGHTS)}
        recovery2=own_second_recovery(drr,first_recovery_dr,z,rho,self.physical.delta)
        physical_recovery2={k:algebra.shift(v,natural.SCALES[k]) for k,v in recovery2.items()}
        result=dict(source_family=self.family,Z=z,query=pv['query'],geometry=pv['geometry'],
            leading_normalized_Phi_Psi_meanPsi_rho3={k:op.export(third[k],2) for k in ('Phi','Psi','meanPsi')},
            leading_physical_Q_rho3=op.export(fields['Q'][3],1),
            normalized_source_rho1={k:op.export(v,1) for k,v in source_rho.items()},
            normalized_B0_B1_x1={label:{k:op.export(v,1) for k,v in entries.items()} for label,entries in mx.items()},
            normalized_source_x1={k:op.export(v,1) for k,v in gx.items()},
            normalized_second_x_function_jets={name:op.export(v,1) for name,v in zip(parent.FIELDS,dxx)},
            normalized_second_rho_function_jets={name:op.export(v,1) for name,v in zip(parent.FIELDS,drr)},
            physical_second_x_function_jets={k:op.export(v,1) for k,v in physical_xx.items()},
            physical_second_rho_function_jets={k:op.export(v,1) for k,v in physical_rr.items()},
            normalized_own_mean_Q_V_rho2={k:op.export(v) for k,v in recovery2.items()},
            physical_own_mean_Q_V_rho2={k:op.export(v) for k,v in physical_recovery2.items()},
            source_rho3_Z2_uses_original_total_order_five_data=True,
            actual_collar_prescribed_shear_ODE_differentiated=True,
            original_hb_inverse_and_inverse_squared_sectors_retained=True,
            second_radial_jets_use_exact_six_equation_derivative=True,
            full_regular_function_and_feedback_enclosures_propagated=True,
            finite_axis_xx_and_even_rho2_limits_used=axis,
            odd_last_two_rho2_derivatives_at_axis_not_reported=axis,
            no_missing_source_grade_padded_with_zero=True,
            full_axis_Taylor_remainder_derivative_joins_own_repair_n2_still_open=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalN1SecondRadialPacket();state=(rows,fields,third,sources,matrices,source_rho,mx,gx,dxx,drr,recovery2)
        # Pair fingerprints consist solely of their genuine value/radial jets.
        frozen=op.snapshot((result,rows,fields,third,pair_records(sources,matrices),source_rho,mx,gx,dxx,drr,recovery2))
        self._packets[id(packet)]=(packet,result,p,state,frozen)
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1SecondRadialPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live original n1 second radial packet required')
        self.regular.report(entry[2]);rows,fields,third,sources,matrices,srho,mx,gx,dxx,drr,recovery2=entry[3]
        if op.snapshot((entry[1],rows,fields,third,pair_records(sources,matrices),srho,mx,gx,dxx,drr,recovery2))!=entry[4]:raise ValueError('Second radial packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_regular):
    owner=CurrentOriginalN1SecondRadialFunction(accepted_regular,require_checked=False)
    packets={name:owner.evaluate(z,r,chart) for name,z,r,chart in (
        ('core_axis','.371','0','core'),('core_interior','.371','1.5','core'),('core_exit','.371','4','core'),
        ('actual_first_inlet','.371','0','first'),('actual_first_interior','.371','.375','first'),
        ('selected_actual_inner_exit','.371','1','first'),('other_Z_core','-.4','2.5','core'))}
    raw=dict(source_family=owner.family,samples={k:owner.report(v) for k,v in packets.items()},
             input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_SECOND_RADIAL_FUNCTIONS_COMPUTED',len(packets),flush=True)
    return owner,packets

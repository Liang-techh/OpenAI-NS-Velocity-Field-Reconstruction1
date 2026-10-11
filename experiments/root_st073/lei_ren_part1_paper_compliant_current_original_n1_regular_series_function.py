"""One common-inner n1 function, regular principal inverse and equation jets.

The function is the convergent original Volterra series. Numeric queries
enclose its zeroth integral and all remaining words, retaining true core
and actual collar sources. Positive majorants never replace its matrices.
Higher radial jets and the full repaired-function admission remain open.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_natural_component_enclosure as natural
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=natural.op;domain=natural.domain;core=natural.core;weighted=natural.weighted;ep=natural.ep
feedback=natural.feedback;HERE,PREFIX,sha=natural.HERE,natural.PREFIX,natural.sha
NAME=PREFIX+'current_original_n1_regular_series_function.json.gz'
RECEIPT=PREFIX+'current_original_n1_regular_series_function_check.json'
GATE='current_original_n1_regular_common_inner_series_and_first_radial_function_enclosures_installed'
OPEN=natural.OPEN
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
POWERS=(2,2,2,2,1,1)
REGULAR_GREEN=64


def theta_source(algebra,phi,z,rho,delta,gradient,order=2):
    """Exact same normalized g5 identity, now through ordinary Z2."""
    f,fr,_=phi;d=1-z*z;L=1-z*z*delta;m=gradient*d/L
    a=-2-delta;b=-3
    cross=m*(f*(z*(a+b))+op.dz(f)*(2*d)-fr*(z*(4*rho)))/L+op.dz(m)*d*f/L
    second=op.axial_second(phi,z,rho,delta,a)
    return op.truncate(-2*m*m*f+algebra.shift(2*cross,(0,0,0,-1))+
                       algebra.shift(-2*second,(0,0,0,-2)),order)


def exact_operator(owner,mean_packet,order=2):
    """True source-bound functions; no WidthBound/norm proxy is used."""
    mv=owner.mean.report(mean_packet);entry=owner.mean._packets[id(mean_packet)]
    rows,fields,forcing=entry[4][1:4];algebra=fields['Uz'][0].algebra;c=owner.c
    rho=mv['geometry'].get('rho',mv['geometry'].get('rho_enclosure_only'));x=c.sqrt(rho)
    z=algebra.lift(op.Jet.variable(c,mv['Z'],5));d=1-z*z;L=1-z*z*owner.physical.delta
    gradient=op.Jet(c,owner.weighted.source.amplitude.gradient_jets(mv['Z']))
    grad=algebra.lift(gradient)
    gauge=algebra.shift(-2*rows['Psi'][0]*d*grad/L,(0,0,0,-1))
    residual=natural.conjugate_extra(algebra,feedback.residual_matrix(algebra,rows['Phi'],fields['Uz'],fields['Q'][0],
        z,rho,owner.physical.delta,gradient,gauge,order))
    B0=dict(residual['B0']);B1=residual['B1']
    for key,sgn in (('1,5',1),('2,6',1),('3,6',-1)):B0[key]=algebra.lift(sgn,order)
    axis=owner.physical.axis_inputs(mv['Z']);chi=axis['chi'].truncate(order)
    B0['5,1']=B0['5,1']+algebra.lift(-2*chi)
    g4_over_x=op.truncate(algebra.shift(forcing['N1p']*2,(0,0,-1,0)),order)
    g={'4':g4_over_x*x,
       '5':theta_source(algebra,rows['Phi'],z,rho,owner.physical.delta,grad,order),
       '6':op.truncate(forcing['N1z']*2-algebra.shift(forcing['N1p']*z*(4*rho)/L,(0,0,0,-1)),order)}
    return dict(B0=B0,B1=B1,C0=residual['B0'],g=g,g4_over_x=g4_over_x,
                D=(0,0,2,0,3,1),chi=chi,x=x,rho=rho,gradient=gradient,algebra=algebra)


def principal_kernel(c,chi,s,x,N=12,order=1):
    """6x3 regular inverse on source rows 4,5,6; tails make it infinite.

    Kchi=sum_j>=0(-2chi)^j T^j K0, T=(dxx+3dx/x)^-1.
    Pressure/axial kernels are exact; the swirl series has an explicit tail.
    """
    if order!=1:raise ValueError('This finite kernel adapter supplies value/Z1 tails only')
    if ep(s)[0]<0 or ep(x)[0]<=0 or ep(s)[1]>ep(x)[1]:raise ValueError('0<=s<=x and x>0 required')
    ratio=s/x;t=c.mpf([max(0,ep(ratio)[0]),min(1,ep(ratio)[1])])
    f,dx=weighted.potential_kernel(chi,t,x,N,order)
    kernels=weighted.initial.kernels(c,s,x)
    values=dict(F1=f,Uz1=kernels['Uz1'],K1=kernels['K1'],P1=kernels['P1'],
                d_x_F1=dx,d_x_Uz1=kernels['d_x_Uz1'])
    tails=weighted.potential_kernel_tail(c,domain.upper(c,x),domain.upper(c,abs(chi[1])),N)
    return values,tails


def regular_series_bound(owner):
    """All words in the axis-weighted space W_i/x^p_i."""
    dp=owner.natural.domain();nv=owner.natural.report(dp);data=owner.natural._data;c=owner.c
    a=data['a'];eta=data['eta'];q=a*a/2
    chi_tube=owner.weighted.tube
    if not chi_tube['chi_holomorphic'] or ep(chi_tube['chi_complex_modulus_upper'])[1]>1:
        raise ValueError('Original holomorphic potential bound not admitted')
    # Row 5, input power p=1: the combined j=0 weighted F/dxF
    # column is t^2*(1-t^2)/2+t^4 <= 1. For j>=1 its sum is
    # <= (j+1/2)*q^j/[j!*(j+1)!] <= q^j/(j-1)!.
    # Thus the whole column is <= 1+q*exp(q). Input power p=2
    # adds at most a. Row 4 gives t^2; row 6 gives at most
    # 1+1/(2e)+1/8, which is below the same uniform cover.
    column_bound=domain.upper(c,domain.upper(c,a)*(1+q*c.exp(q)))
    if ep(column_bound)[1]>=REGULAR_GREEN:raise ValueError('Regular weighted Green bound not admitted')
    # Recover pointwise core bounds from the actual source construction,
    # rather than treating an arbitrary L1 norm as a pointwise norm. In
    # natural_majorants(), core_L1 is exactly 2*the uniform source bound,
    # and row 4 is exactly 2*a*Np/Pstar^2. Hence g4/x has the bound below.
    # The pressure kernel integrates s*(g4/s); the other two source rows
    # have constant core bounds. Their zero-axis principal responses,
    # divided by x^p, are bounded by REGULAR_GREEN times this source sum
    # on x<=1. For x>=1 the already admitted unweighted M0 is sufficient.
    core_source=data['core_L1']
    near_axis_source={
        'g4_over_x':core_source['g4']/(2*a),
        'g5':core_source['g5']/2,
        'g6':core_source['g6']/2}
    near_axis_cover=feedback.directed_sector_cover(
        sum(near_axis_source.values(),core_source['g4'].algebra.lift(0,0)),floor=natural.FLOOR)
    if ep(near_axis_cover['modulus_upper'])[1]>ep(data['M0'])[0]:
        raise ValueError('Core pointwise regular source norm not dominated by M0')
    Mreg=domain.upper(c,REGULAR_GREEN*data['M0'])
    tail=feedback.all_word_tail(c,Mreg,a,REGULAR_GREEN*data['C0cover']['modulus_upper'],
        REGULAR_GREEN*data['C1cover']['modulus_upper'],eta/16,0)
    gap=eta/32
    tail.update(regular_axis_powers=list(POWERS),value_tube_radius=eta/16,
        original_holomorphic_chi_tube=copy.deepcopy(chi_tube),
        regular_Green_column_proof='row5 j0: t^2(1-t^2)/2+t^4<=1; j>=1: (j+1/2)/(j!*(j+1)!)<=1/(j-1)!; extra input s<=a; row4<=1; row6<=1+1/(2e)+1/8',
        core_pointwise_regular_source_cover=near_axis_cover,
        core_pointwise_regular_source_rows={k:op.export(v,0) for k,v in near_axis_source.items()},
        core_pointwise_bound_recovery='core_L1=2*core_uniform_source; g4_uniform=2*a*Np/Pstar^2; g4/x bound=core_L1[g4]/(2*a)',
        core_pointwise_regular_source_norm_dominated_by_M0=True,
        initial_weighted_norm_proof='x<=1: pressure integral of s*(g4/s), constant g5/g6 bounds and zero-axis kernels; x>=1: x^(-p)<=1 and admitted unweighted M0',
        derivative_tube_radius=eta/32,derivative_Cauchy_gap=gap,
        source_M0_tube_radius=eta/8,regular_initial_M0=Mreg,
        regular_solution_bound=domain.upper(c,Mreg+tail['value_remainder_upper']),
        regular_Green_column_majorant=column_bound,regular_Green_factor=REGULAR_GREEN,
        ordinary_Z1_remainder_upper=domain.upper(c,tail['value_remainder_upper']/gap),
        ordinary_Z2_remainder_upper=domain.upper(c,2*tail['value_remainder_upper']/gap**2),
        ordinary_Z3_remainder_upper=domain.upper(c,6*tail['value_remainder_upper']/gap**3))
    return dp,tail


def symmetric_error(value):
    c=value.ctx
    return op.FactoredJet(value.algebra,{p:op.Jet(c,[c.mpf([-1,1])*domain.upper(c,abs(v)) for v in row.coefficients])
                                       for p,row in value.terms.items()},value.order)


def extend_by_analytic_bound(value,order,bound,gap):
    """Uncomputed higher jets enclose the full function by Cauchy."""
    if value.order>=order:return op.truncate(value,order)
    c=value.ctx;terms={p:op.Jet(c,list(row.coefficients)+[c.mpf(0)]*(order-value.order)) for p,row in value.terms.items()}
    key=(0,0,0,0);coefficients=[c.mpf(0)]*(order+1)
    for k in range(value.order+1,order+1):coefficients[k]=c.mpf([-1,1])*domain.upper(c,bound/gap**k)
    addition=op.Jet(c,coefficients)
    terms[key]=terms[key]+addition if key in terms else addition
    return op.FactoredJet(value.algebra,terms,order)


def equation_radial_jets(system,W,order=2):
    """d_x W and mixed Z jets from the exact six-row equation."""
    algebra=system['algebra'];c=algebra.ctx;x=system['x']
    if ep(x)[1]==0:
        zero=algebra.lift(0,order)
        dx=[zero,zero,zero,zero,system['g']['5']/4,system['g']['6']/2]
        # Only the four even components have a finite rho derivative at 0.
        dr=[system['g']['5']/8,system['g']['6']/4,-system['g']['6']/8,system['g4_over_x']/2]
        return dx,dr
    dx=[]
    for i in range(6):
        value=system['g'].get(str(i+1),algebra.lift(0,order))-op.truncate(W[i],order)*(system['D'][i]/x)
        for key,B in system['B0'].items():
            row,column=(int(v)-1 for v in key.split(','))
            if row==i:value=value+B*op.truncate(W[column],order)
        for key,B in system['B1'].items():
            row,column=(int(v)-1 for v in key.split(','))
            if row==i:value=value+B*op.truncate(op.dz(W[column]),order)
        dx.append(op.truncate(value,order))
    return dx,[v/(2*x) for v in dx]


def recovered_fields(W,dr,Z,rho,delta,axis=False):
    algebra=W[0].algebra;c=algebra.ctx;z=algebra.lift(op.Jet.variable(c,Z,3));d=1-z*z;L=1-z*z*delta
    U,K=W[1],W[2];M=U+K
    Q=(U*(z*(1-delta))-K*(z*(1+delta))-op.dz(M)*d)/L
    mr=dr[1]+dr[2]
    qr=(dr[1]*(z*(1-delta))-dr[2]*(z*(1+delta))-op.dz(mr)*d)/L
    return dict(mean1=M,Q1=Q,V1=Q*rho),dict(mean1=mr,Q1=qr,V1=Q+qr*rho)


@dataclass(frozen=True,eq=False)
class OriginalN1RegularSeriesFunctionPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1RegularSeriesFunction:
    @source_precision
    def __init__(self,accepted_natural,require_checked=True):
        if type(accepted_natural) is not natural.CurrentOriginalN1NaturalComponentEnclosure or not accepted_natural.acceptance_loaded:
            raise ValueError('Accepted natural original source and full feedback owner required')
        accepted_natural.assert_graph();self.natural=accepted_natural;self.mean=accepted_natural.mean
        self.c=accepted_natural.c;self.physical=accepted_natural.physical;self.weighted=accepted_natural.weighted
        self.core=self.mean.axial.inner.core;self.collar=self.mean.axial.inner
        self.family=copy.deepcopy(accepted_natural.family);self.hashes=dict(accepted_natural.hashes)
        for name in (natural.NAME,natural.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._parents=(accepted_natural,self.mean,self.core,self.collar);self._identities=(self.c,self.physical,self.weighted)
        self.bound_packet,self.tail=regular_series_bound(self);self._tail=op.snapshot(self.tail)
        self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Regular series function receipt differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Regular series source changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if ((self.natural,self.mean,self.core,self.collar)!=self._parents or (self.c,self.physical,self.weighted)!=self._identities or
            self.mean is not self.natural.mean or self.core is not self.mean.axial.inner.core or self.collar is not self.mean.axial.inner or
            not self.natural.acceptance_loaded or self.hashes!=self._hashes or self.family!=self.natural.family or op.snapshot(self.tail)!=self._tail):
            raise ValueError('Regular series function graph changed')
        self.natural.report(self.bound_packet);return self.natural.assert_graph()

    @source_precision
    def contract(self,N=0):
        self.assert_graph();c=self.c
        if isinstance(N,bool) or not isinstance(N,int) or N<0:raise ValueError('Nonnegative integer word index required')
        data=self.natural._data;t=self.tail;gap=t['derivative_Cauchy_gap']
        term=domain.upper(c,t['regular_initial_M0']*data['green']*t['ratio_upper']**N)
        residual=domain.upper(c,data['a']**2*(data['C0cover']['modulus_upper']+data['C1cover']['modulus_upper']/gap)*term)
        return dict(source_family=self.family,exact_function_definition='W=sum_k>=0 V^k W0; W0=Gprincipal*g; V f=Gprincipal*(C0*f+B1*dZ(f))',
            exact_principal_source_columns=[4,5,6],exact_principal_unknowns=list(FIELDS),
            principal_kernel_identity='Gprincipal: [Kchi*h5, s*log(x/s)*h6, -s*(1-s^2/x^2)*h6/2, h4, dxKchi*h5, (s/x)*h6] integrated from 0 to x',
            exact_infinite_swirl_kernel='Kchi=sum_j>=0(-2chi)^j*T^j*K0; K0=s*(1-s^2/x^2)/2',
            principal_regular_space='W_i/x^p_i bounded; p=(2,2,2,2,1,1); even first four and odd last two at axis',
            exact_source_and_matrix_callbacks='exact_operator(original leading own-mean packet); norm proxies are separate',
            original_core_first_seam='same core inlet functions; one integral path 0..x; each V word inherits seam without reset',
            finite_equation_residual_identity='L(sum_0^N Wk)-g=-(C0*WN+B1*dZ(WN))',
            finite_equation_residual_upper_on_eta_over32=residual,finite_word_index=N,
            regular_series_bound=t,normalization='Wphysical=Snew*W; original F0 Z derivative applied once',
            numeric_centers_compute_only_zeroth_source=True,all_remaining_words_enclosed_not_set_zero=True,
            full_higher_radial_axis_and_join_contracts_still_open=True,n1_own_repair_n2_global_and_NS_still_open=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))

    @source_precision
    def evaluate(self,Z,coordinate,chart='core'):
        self.assert_graph();began=time.monotonic();c=self.c
        np=self.natural.evaluate(Z,coordinate,chart);nv=self.natural.report(np);ne=self.natural._packets[id(np)]
        mp=ne[2];system=exact_operator(self,mp);algebra=system['algebra'];rho=system['rho'];x=system['x'];z=nv['Z']
        if chart=='core':
            sp=self.core.evaluate(Z,coordinate);sv=self.core.report(sp);se=self.core._packets[id(sp)];source_values,source_errors=se[3:5];source_owner=self.core
        elif chart=='first':
            sp=self.collar.evaluate(Z,coordinate);sv=self.collar.report(sp);se=self.collar._packets[id(sp)];source_values,source_errors=se[6:8];source_owner=self.collar
        else:raise ValueError('Original core/selected actual first collar required')
        if any(v.algebra is not algebra for v in source_values.values()):raise ValueError('Original source amplitude frames differ')
        nvalues,nerrors=ne[3][:2];centers={**{k:nvalues[k] for k in ('Uz1','K1','P1','d_x_Uz1')},**source_values}
        errors={**{k:nerrors[k] for k in ('Uz1','K1','P1','d_x_Uz1')},**source_errors}
        W=[];total_errors={};gap=self.tail['derivative_Cauchy_gap']
        for name,power in zip(FIELDS,POWERS):
            shape=rho if power==2 else x
            extra=algebra.lift(op.Jet(c,[domain.upper(c,shape*self.tail['value_remainder_upper']/gap**k) for k in range(centers[name].order+1)]))
            total_errors[name]=errors[name]+extra
            enclosed=centers[name]+symmetric_error(total_errors[name])
            W.append(extend_by_analytic_bound(enclosed,3,shape*self.tail['regular_solution_bound'],gap))
        dx,dr=equation_radial_jets(system,W,2)
        recovery,recovery_dr=recovered_fields(W,dr,z,rho,self.physical.delta,ep(x)[1]==0)
        physical={name:weighted.transform(v,system['gradient'],scale,3) for name,v,scale in zip(FIELDS,W,natural.WEIGHTS)}
        physical_dx={name:weighted.transform(v,system['gradient'],scale,2) for name,v,scale in zip(FIELDS,dx,natural.WEIGHTS)}
        physical_dr={name:weighted.transform(v,system['gradient'],scale,2) for name,v,scale in zip(FIELDS,dr,natural.WEIGHTS)}
        physical_recovery={k:algebra.shift(v,natural.SCALES[k]) for k,v in recovery.items()}
        physical_recovery_dr={k:algebra.shift(v,natural.SCALES[k]) for k,v in recovery_dr.items()}
        # V1=epsilon^2*rho*Qnorm; its rho derivative already includes Qnorm.
        budgets={name:{'Z'+str(k):natural.natural_budget(centers[name],total_errors[name],k) for k in (0,1)} for name in FIELDS}
        result=dict(source_family=self.family,Z=z,query=nv['query'],geometry=nv['geometry'],
            normalized_function_jets={name:op.export(v) for name,v in zip(FIELDS,W)},
            normalized_signed_zeroth_centers={name:op.export(v) for name,v in centers.items()},
            normalized_center_total_errors={name:op.export(v) for name,v in total_errors.items()},
            normalized_d_x_function_jets={name:op.export(v) for name,v in zip(FIELDS,dx)},
            normalized_d_rho_function_jets={name:op.export(v) for name,v in zip(FIELDS,dr)},
            physical_function_jets={k:op.export(v) for k,v in physical.items()},
            physical_d_x_function_jets={k:op.export(v) for k,v in physical_dx.items()},
            physical_d_rho_function_jets={k:op.export(v) for k,v in physical_dr.items()},
            physical_own_mean_Q_V={k:op.export(v) for k,v in physical_recovery.items()},
            physical_own_mean_Q_V_rho={k:op.export(v) for k,v in physical_recovery_dr.items()},
            exact_normalized_equation={label:{key:op.export(v,2) for key,v in system[label].items()} for label in ('B0','B1','g')},
            component_budgets=budgets,passing_six_component_value_Z1_rows=sum(row['relative_natural_component_width_certified'] for rows in budgets.values() for row in rows.values()),
            same_original_common_inner_function_for_every_query=True,all_six_row_feedback_included=True,
            radial_derivatives_recovered_from_exact_equation_not_seed_differences=True,
            finite_axis_x_and_even_rho_limits_used=ep(x)[1]==0,
            odd_last_two_rho_derivatives_at_axis_not_reported=ep(x)[1]==0,
            ordinary_Z3_function_bound_available=True,first_radial_mixed_Z2_available=True,
            uncomputed_higher_Z_centers_use_symmetric_analytic_bound=True,
            original_source_F0_Z_gauge_applied_once_in_physical_output=True,
            full_higher_radial_axis_and_join_contracts_still_open=True,n1_own_repair_and_n2_still_open=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1RegularSeriesFunctionPacket()
        state=(W,dx,dr,recovery,recovery_dr,{k:v for k,v in system.items() if k!='algebra'})
        self._packets[id(token)]=(token,result,np,source_owner,sp,state,op.snapshot((result,state)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1RegularSeriesFunctionPacket or entry is None or entry[0] is not packet:raise ValueError('Live common-inner series function packet required')
        self.natural.report(entry[2]);entry[3].report(entry[4])
        if op.snapshot((entry[1],entry[5]))!=entry[6]:raise ValueError('Series function packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_natural):
    owner=CurrentOriginalN1RegularSeriesFunction(accepted_natural,require_checked=False)
    packets={name:owner.evaluate(z,r,chart) for name,z,r,chart in (
        ('core_axis','.371','0','core'),('core_interior','.371','1.5','core'),('core_exit','.371','4','core'),
        ('actual_first_inlet','.371','0','first'),('actual_first_interior','.371','.375','first'),
        ('selected_actual_inner_exit','.371','1','first'),('other_Z_core','-.4','2.5','core'))}
    raw=dict(source_family=owner.family,contract=owner.contract(),samples={k:owner.report(v) for k,v in packets.items()},
             input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_REGULAR_SERIES_FUNCTIONS_COMPUTED',len(packets),flush=True)
    return owner,packets

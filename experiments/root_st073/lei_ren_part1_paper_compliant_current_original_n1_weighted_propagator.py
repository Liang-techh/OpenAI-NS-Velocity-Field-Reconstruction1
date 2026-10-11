"""Original F0*Lambda-normalized n1 operator and resummed swirl potential.

The -2 chi gauge potential is solved in the regular radial propagator.
Actual source integrals preserve amplitude sectors and include a bounded
kernel-series tail. Remaining operator feedback is NOT solved here.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_component_budget as budget
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

leading=budget.leading;op=budget.op;domain=budget.domain;ep=budget.ep
initial=budget.coupled.initial
HERE,PREFIX,sha=budget.HERE,budget.PREFIX,budget.sha
NAME=PREFIX+'current_original_n1_weighted_propagator.json.gz'
RECEIPT=PREFIX+'current_original_n1_weighted_propagator_check.json'
GATE='current_original_n1_weighted_potential_propagator_and_source_integrals_computed'
OPEN=budget.OPEN
WEIGHTS=((0,1,0,1),(0,0,0,0),(0,0,0,0),(0,0,0,0),(0,1,0,1),(0,0,0,0))


def relative_amplitude(algebra,gradient,power,order):
    """Jet of F0(Z)^power/F0(Z0)^power, with Lambda kept formal."""
    if isinstance(power,bool) or not isinstance(power,int) or order>gradient.order:
        raise ValueError('Integer amplitude power and genuine gradient derivatives required')
    c=algebra.ctx;coefficients=[{0:c.mpf(1)}]
    for k in range(order):
        row={}
        for j in range(k+1):
            for p,v in coefficients[k-j].items():
                row[p+1]=row.get(p+1,c.mpf(0))-power*gradient[j]*v/(k+1)
        coefficients.append(row)
    terms={}
    for k,row in enumerate(coefficients):
        for p,v in row.items():
            if ep(v)==(0,0):continue
            key=(0,0,0,p)
            if key not in terms:terms[key]=[c.mpf(0)]*(order+1)
            terms[key][k]+=v
    return op.FactoredJet(algebra,{p:op.Jet(c,v) for p,v in terms.items()},order)


def transform(value,gradient,difference,order):
    value=op.truncate(value,order)
    return value.algebra.shift(value,tuple(difference))*relative_amplitude(value.algebra,gradient,difference[1],order)


def conjugate(system,gradient,order=2):
    """W=S What; B0*=S^-1 B0 S+S^-1 B1 S_Z, B1*=S^-1 B1 S."""
    algebra=next(iter(system['g'].values())).algebra
    ell=algebra.shift(algebra.lift(-op.Jet(algebra.ctx,gradient.coefficients[:order+1])),(0,0,0,1))
    result={label:{} for label in ('B0','B1','g')}
    for label in ('B0','B1'):
        for key,value in system[label].items():
            i,j=(int(v)-1 for v in key.split(','))
            diff=tuple(b-a for a,b in zip(WEIGHTS[i],WEIGHTS[j]))
            result[label][key]=transform(value,gradient,diff,order)
    for key,value in result['B1'].items():
        j=int(key.split(',')[1])-1
        if WEIGHTS[j][1]:
            result['B0'][key]=result['B0'].get(key,algebra.lift(0,order))+value*ell*WEIGHTS[j][1]
    for key,value in system['g'].items():
        result['g'][key]=transform(value,gradient,tuple(-p for p in WEIGHTS[int(key)-1]),order)
    return result


@lru_cache(None)
def kernel_polynomials(N):
    """K_j(x,s)=x^(2j+1) k_j(s/x), K_{j+1}=K0 composed with K_j.

    Exact rational integration of powers/logs, with no sampled fitting.
    """
    if isinstance(N,bool) or not isinstance(N,int) or not 0<=N<=64:
        raise ValueError('Bounded nonnegative potential-series cutoff required')
    rows=[{(1,0):Fraction(1,2),(3,0):Fraction(-1,2)}]
    def add(out,key,value):out[key]=out.get(key,Fraction(0))+value
    for j in range(N):
        out={}
        for (a,b),v in rows[-1].items():
            for m,sgn in ((2*j+3-a,1),(2*j+5-a,-1)):
                factor=v*sgn/2
                if b==0:
                    if m==0:add(out,(a,1),-factor)
                    else:
                        add(out,(a,0),factor/m);add(out,(a+m,0),-factor/m)
                elif b==1 and m!=0:
                    add(out,(a,1),factor/m);add(out,(a,0),factor/(m*m));add(out,(a+m,0),-factor/(m*m))
                else:raise ValueError('Unexpected resonant logarithm in radial kernel')
        rows.append({p:v for p,v in out.items() if v})
    derivatives=[]
    for j,row in enumerate(rows):
        out={}
        for (a,b),v in row.items():
            add(out,(a,b),v*(2*j+1-a))
            if b:add(out,(a,b-1),-b*v)
        derivatives.append({p:v for p,v in out.items() if v})
    return tuple(rows),tuple(derivatives)


def power_log_range(c,t,p,b):
    if ep(t)[0]<0 or ep(t)[1]>1 or p<1 or b not in (0,1):
        raise ValueError('Regular normalized kernel variable in[0,1] required')
    if b==0:return t**p
    if ep(t)[0]>0:return t**p*c.ln(t)
    # t^p log(t) has its minimum at exp(-1/p) and tends to zero at t=0.
    return -c.mpf([0,ep(domain.upper(c,1/(c.e*p)))[1]])


def polynomial_range(c,t,row):
    return sum((c.mpf(v.numerator)/v.denominator*power_log_range(c,t,p,b)
                for (p,b),v in row.items()),c.mpf(0))


def potential_kernel(chi,t,x,N=12,order=1):
    """Signed finite kernel jets; the separate tail covers all later j."""
    c=chi.ctx;rows,derivatives=kernel_polynomials(N);chi=chi.truncate(order)
    f=op.Jet.constant(c,0,order);dx=op.Jet.constant(c,0,order)
    for j,(row,derivative) in enumerate(zip(rows,derivatives)):
        scale=(-2*chi)**j
        f+=scale*x**(2*j+1)*polynomial_range(c,t,row)
        dx+=scale*x**(2*j)*polynomial_range(c,t,derivative)
    return f,dx


def potential_kernel_tail(c,a,chi_Z_norm,N=12):
    """Uniform pointwise kernel tails for real 0<=s<=x<=a and |chi|<=1.

    K_j=T^j K0 is positive and K0<=a/2. T^j1=x^(2j)/[4^j j!(j+1)!].
    Its x derivative is positive too. Differentiate chi^j directly for Z1.
    """
    kernel_polynomials(N)
    if ep(a)[0]<=0 or ep(chi_Z_norm)[0]<0:raise ValueError('Positive radius and derivative norm required')
    j=N+1;r=domain.upper(c,a*a/2);term=domain.upper(c,r**j/(math.factorial(j)*math.factorial(j+1)))
    q0=domain.upper(c,r/((j+1)*(j+2)));q1=domain.upper(c,r/(j*(j+2)))
    q2=domain.upper(c,r*(j+1)/(j*j*(j+2)))
    if any(ep(q)[1]>=1 for q in (q0,q1,q2)):raise ValueError('Kernel tail ratios not admitted')
    return dict(F1=op.Jet(c,[domain.upper(c,a*term/(2*(1-q0))),
                            domain.upper(c,a*chi_Z_norm*j*term/(2*(1-q1)))]),
                d_x_F1=op.Jet(c,[domain.upper(c,j*term/(1-q1)),
                                 domain.upper(c,chi_Z_norm*j*j*term/(1-q2))]),
                first_omitted_potential_index=j,ratio_bounds=[q0,q1,q2],
                bound_scope='potential-kernel truncation only; not remaining B0/B1 feedback',
                direct_chi_derivative_used_not_Cauchy_loss=True)


def chi_tube_certificate(c,sigma,j,delta,eta):
    # H=-4Z^3-jZ^2+(4+(1-delta)/2)Z+j on the original complex tube.
    b=domain.upper(c,eta*(12*(1+eta)**2+4*eta**2+2*abs(j)*(1+eta)+4+abs(1-delta)/2))
    gap=domain.lower(c,sigma*sigma-2*b*b)
    if ep(gap)[0]<=0:raise ValueError('Original complex chi tube not admitted')
    return dict(original_complex_Z_radius=eta,H_imaginary_modulus_upper=b,
                sigma_squared_minus_twice_H_imaginary_squared_lower=gap,
                chi_holomorphic=True,chi_complex_modulus_upper=c.mpf(1),
                proof='Re(H^2+sigma^2)>0; |H|^2<=Re(H^2+sigma^2) since 2(Im H)^2<sigma^2',
                original_tube_preserved_not_shrunk=True)


def absolute_jet(value,order=1):
    c=value.ctx
    return op.FactoredJet(value.algebra,{p:op.Jet(c,[domain.upper(c,abs(v)) for v in jet.coefficients[:order+1]])
                                      for p,jet in value.terms.items()},order)


def kernel_sector_budget(value,error,k):
    c=value.ctx;prefix=budget.component_budget(value,k,c.mpf(0));norm=c.mpf(0);terms=[]
    if prefix['prefix_exact_zero']:
        return dict(prefix=prefix,kernel_tail_exact_zero=not error.terms,
                    remaining_B0_B1_feedback_included=False)
    ref=prefix['reference_source_exponents'];M=prefix['directed_reference_coefficient_norm'];unresolved=False
    for p,jet in error.terms.items():
        v=domain.upper(c,abs(jet[k]))
        if ep(v)[1]==0:continue
        log=budget.log_ratio(c,value.algebra.logs,p,ref,v,M);upper,rule=budget.exp_cover(c,log)
        terms.append(dict(source_exponents=list(p),normalized_log_upper=log,norm_upper=upper,bound_rule=rule))
        if upper is None:unresolved=True
        else:norm+=upper
    return dict(prefix=prefix,kernel_tail_normalized_reference_norm_upper=None if unresolved else domain.upper(c,norm),
                kernel_tail_sector_contributions=terms,normalization_is_directed_source_bound_not_field_value=True,
                remaining_B0_B1_feedback_included=False,source_and_quadrature_uncertainty_not_separated=True)


@dataclass(frozen=True,eq=False)
class OriginalN1WeightedPropagatorPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1WeightedPropagator:
    @source_precision
    def __init__(self,accepted_budget,require_checked=True):
        if type(accepted_budget) is not budget.CurrentOriginalN1ComponentBudget or not accepted_budget.acceptance_loaded:
            raise ValueError('Accepted current component budget required')
        accepted_budget.assert_graph();self.budget=accepted_budget;self.local=accepted_budget.local
        self.inputs=self.local.inputs;self.c=self.local.c;self.source=self.local.source
        self.family=copy.deepcopy(self.local.family);self.hashes=dict(accepted_budget.hashes)
        for name in (budget.NAME,budget.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.physical=self.inputs.n1.source.core
        bound=self.inputs.bounds.report(self.inputs.bounds.domain())
        self.tube=chi_tube_certificate(self.c,self.physical.sigma,self.physical.j,self.physical.delta,
                                       bound['common_tubes']['outer_Z_radius'])
        self.a=domain.upper(self.c,self.c.sqrt(bound['actual_inner_geometry']['real_rho_upper_bound']))
        self._owner=accepted_budget;self._inputs=self.inputs;self._physical=self.physical;self._c=self.c
        self._state=op.snapshot((self.family,self.tube,self.a));self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Weighted propagator receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Weighted propagator dependency changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.budget is not self._owner or self.inputs is not self._inputs or self.physical is not self._physical
            or self.physical is not self.inputs.n1.source.core or self.c is not self._c or self.c is not self.local.c
            or self.local is not self.budget.local or not self.budget.acceptance_loaded
            or self.hashes!=self._hashes or self.family!=self.local.family
            or self._state!=op.snapshot((self.family,self.tube,self.a))):raise ValueError('Weighted source graph changed')
        return self.budget.assert_graph()

    def _operator(self,z,actual,chart):
        p=self.inputs.evaluate(z,actual,chart);view=self.inputs.report(p);entry=self.inputs._packets[id(p)]
        algebra=entry[5];gradient=op.Jet(self.c,self.source.amplitude.gradient_jets(z))
        axis=self.physical.axis_inputs(z);system=conjugate(entry[4],gradient,2)
        chi=op.Jet(self.c,axis['chi'].coefficients[:3]);potential=algebra.lift(-2*chi)
        residual=copy.copy(system['B0']);residual['5,1']=residual['5,1']-potential
        for key,sgn in (('1,5',1),('2,6',1),('3,6',-1)):
            residual[key]=residual[key]-algebra.lift(sgn,2)
        return p,view,system,gradient,chi,potential,residual

    @source_precision
    def evaluate(self,Z,coordinate,chart='first',core_cells=16,collar_cells=4,N=12):
        self.assert_graph();c=self.c
        endpoint=c.mpf(coordinate)
        if not isinstance(coordinate,(str,int,float)) and ep(endpoint)[0]!=ep(endpoint)[1]:
            raise ValueError('One exact outer endpoint required; only integration cells may be intervals')
        for count in (core_cells,collar_cells):
            if isinstance(count,bool) or not isinstance(count,int) or not 1<=count<=256:raise ValueError('Bounded positive partition required')
        kernel_polynomials(N)
        lp=self.local.evaluate(Z,coordinate,chart);lv=self.local.report(lp);q=lv['query'];z=lv['Z']
        p,view,system,gradient,chi,potential,residual=self._operator(z,q['actual_source_coordinate'],chart)
        algebra=self.inputs._packets[id(p)][5];rho=view['geometry'].get('rho',view['geometry'].get('rho_enclosure_only'))
        X=c.sqrt(rho);phase=q['actual_source_coordinate'] if chart=='first' else c.mpf(0)
        core_end=c.mpf(2) if chart=='first' else X
        sums={field:algebra.lift(0,1) for field in initial.FIELDS}
        errors={field:algebra.lift(0,1) for field in ('F1','d_x_F1')}
        kernel_tail=potential_kernel_tail(c,self.a,domain.upper(c,abs(chi[1])),N)
        ledger=[];began=time.monotonic()
        def cell_integral(packet,xcell,measure,collar):
            sv=self.inputs.report(packet);original=self.inputs._packets[id(packet)][4]
            sg=conjugate(original,gradient,1)['g']
            raw_t=xcell/X;t=c.mpf([max(0,ep(raw_t)[0]),min(1,ep(raw_t)[1])])
            f,dx=potential_kernel(chi,t,X,N,1)
            weights=initial.kernels(c,xcell,X)
            for field,row in (('F1','5'),('Uz1','6'),('K1','6'),('P1','4'),('d_x_F1','5'),('d_x_Uz1','6')):
                kernel=f if field=='F1' else dx if field=='d_x_F1' else weights[field]
                value=sg[row]*kernel*measure
                if collar:value=algebra.width(value)
                if any(p[0]<0 for p in value.terms):raise ValueError('Original inverse width not cancelled in source integral')
                sums[field]=sums[field]+op.truncate(value,1)
            for field in errors:
                value=absolute_jet(sg['5'])*kernel_tail[field]*measure
                if collar:value=algebra.width(value)
                if any(p[0]<0 for p in value.terms):raise ValueError('Inverse width not cancelled in kernel tail')
                errors[field]=errors[field]+op.truncate(value,1)
            ledger.append(dict(chart=sv['geometry']['chart'],rho_cell=sv['geometry'].get('rho',sv['geometry'].get('rho_enclosure_only')),
                               exact_ordered_cell_ratio_intersection=True))
        if ep(core_end)[1]>0:
            for j in range(core_cells):
                t=budget.coupled.partition(c,j,core_cells);xcell=core_end*t
                cell=self.inputs.evaluate(z,(rho*t*t if chart=='core' else 4*t*t),'core')
                cell_integral(cell,xcell,core_end/core_cells,False)
        if ep(phase)[1]>0:
            for j in range(collar_cells):
                s=phase*budget.coupled.partition(c,j,collar_cells);cell=self.inputs.evaluate(z,s,'first')
                sv=self.inputs.report(cell);xcell=c.sqrt(sv['geometry']['rho_enclosure_only'])
                cell_integral(cell,xcell,xcell*phase/(2*collar_cells),True)
        physical={field:transform(value,gradient,WEIGHTS[i],1) for i,(field,value) in enumerate(sums.items())}
        physical_errors={field:algebra.shift(absolute_jet(relative_amplitude(algebra,gradient,1,1))*value,WEIGHTS[initial.FIELDS.index(field)])
                         for field,value in errors.items()}
        diagnostics={field:{'Z'+str(k):kernel_sector_budget(sums[field],errors[field],k) for k in (0,1)}
                     for field in errors}
        result=dict(source_family=self.family,Z=z,query=q,
            normalization='W=diag(F0*Lambda,1,1,1,F0*Lambda,1)*What',source_weights=[list(p) for p in WEIGHTS],
            normalized_system={label:{key:op.export(v,2) for key,v in system[label].items()} for label in system},
            resummed_potential=op.export(potential,2),residual_B0={key:op.export(v,2) for key,v in residual.items()},
            resummed_kinematic_entries={'1,5':1,'2,6':1,'3,6':-1},
            remaining_feedback_equation='(B0star-Akin-Achi)*What+B1star*d_Z(What)',
            original_chi_axial_jet=list(chi.coefficients),original_chi_complex_tube=self.tube,
            potential_equation='f_xx+3*f_x/x+2*chi*f=h5; f(0)=f_x(0)=0',
            signed_kernel_formula='sum_j=0^N (-2chi)^j K0^(j+1); K0(x,s)=s*(1-s^2/x^2)/2',
            retained_potential_indices=[0,N],potential_kernel_tail={k:(list(v.coefficients) if isinstance(v,op.Jet) else v) for k,v in kernel_tail.items()},
            normalized_initial_integrals={field:op.export(v,1) for field,v in sums.items()},
            physical_initial_integrals={field:op.export(v,1) for field,v in physical.items()},
            normalized_kernel_tail_by_source_sector={field:op.export(v,1) for field,v in errors.items()},
            physical_kernel_tail_by_source_sector={field:op.export(v,1) for field,v in physical_errors.items()},
            normalized_kernel_sector_budgets=diagnostics,
            ordinary_Z_initial_derivative_order=1,integration_cells=ledger,
            finite_potential_kernel_tail_included_as_separate_source_sector_bound=True,
            known_F0_source_derivatives_not_regauged=True,unknown_weight_derivative_added_once=True,
            F0_or_source_width_not_materialized=True,initial_kinematic_and_gauge_potential_resummed=True,
            remaining_normalized_B0_B1_feedback_not_solved=True,full_sector_preserving_coupled_tail_still_open=True,
            component_relative_accuracy_not_promoted=True,old_coupled_prefix_not_replaced_by_initial_integral=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1WeightedPropagatorPacket();self._packets[id(token)]=(token,result,lp,p,sums,errors,op.snapshot((result,sums,errors)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1WeightedPropagatorPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live weighted source packet required')
        self.local.report(entry[2]);self.inputs.report(entry[3])
        if op.snapshot((entry[1],entry[4],entry[5]))!=entry[6]:raise ValueError('Weighted propagator packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_budget):
    owner=CurrentOriginalN1WeightedPropagator(accepted_budget,require_checked=False)
    packets={name:owner.evaluate('.371',coordinate,chart) for name,coordinate,chart in
             (('core_axis','0','core'),('selected_actual_inner_exit','1','first'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(p) for name,p in packets.items()},input_hashes=owner.hashes,
             **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_WEIGHTED_PROPAGATOR_SOURCE_INTEGRALS_COMPUTED',len(packets),flush=True)
    return owner,packets

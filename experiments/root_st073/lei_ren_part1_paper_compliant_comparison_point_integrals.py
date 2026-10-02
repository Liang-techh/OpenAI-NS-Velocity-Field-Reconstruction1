"""Signed original smoothed-comparison histories in exact formal hb factors.

Actual coupled core atoms feed both microscopic smoothing intervals and the
frozen macro. Second-order width coefficients retain a separately bounded
third-order remainder. Numerical caps constrain errors only, not the width.
This is the comparison field; actual prescribed-shear F/V remain separate.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_core_integral_atoms import CompliantCoreIntegralAtoms,bessel_radial_tail_coefficients
from lei_ren_part1_paper_compliant_bridge_mixed_C4 import CompliantBridgeMixedC4
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_compliant_core_physical_field import symmetric,intersection
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def upper(c,value):return c.mpf(endpoints(c.mpf(value))[1])


def jet_norm(jet,weight):
    return sum((upper(jet.ctx,abs(v))*weight**k for k,v in enumerate(jet.coefficients)),jet.ctx.mpf(0))


class WidthPolynomial:
    """Sum_{j=0}^2 hb^j a_j(Z) + hb^3 error, in a finite axial jet ring."""
    def __init__(self,rows,remainder,weight,cap):
        self.rows=list(rows);self.ctx=rows[0].ctx;self.weight=weight;self.cap=cap
        self.remainder=upper(self.ctx,remainder)
    def coerce(self,value):
        if isinstance(value,WidthPolynomial):
            if value.weight._mpi_!=self.weight._mpi_ or value.cap._mpi_!=self.cap._mpi_:
                raise ValueError('Different source norm/width cap')
            return value
        row=value if isinstance(value,IntervalTaylor) else self.rows[0]*0+value
        return WidthPolynomial([row,row*0,row*0],0,self.weight,self.cap)
    def __add__(self,value):
        other=self.coerce(value)
        return WidthPolynomial([a+b for a,b in zip(self.rows,other.rows)],self.remainder+other.remainder,self.weight,self.cap)
    __radd__=__add__
    def __neg__(self):return WidthPolynomial([-r for r in self.rows],self.remainder,self.weight,self.cap)
    def __sub__(self,value):return self+-self.coerce(value)
    def __rsub__(self,value):return self.coerce(value)+-self
    def __mul__(self,value):
        other=self.coerce(value);zero=self.rows[0]*0
        rows=[sum((self.rows[i]*other.rows[k-i] for i in range(k+1)),zero) for k in range(3)]
        left=[jet_norm(r,self.weight) for r in self.rows];right=[jet_norm(r,self.weight) for r in other.rows]
        discarded=sum((left[i]*right[j]*self.cap**(i+j-3) for i in range(3) for j in range(3) if i+j>=3),self.ctx.mpf(0))
        remainder=(discarded+self.remainder*sum((right[j]*self.cap**j for j in range(3)),self.ctx.mpf(0))
                   +other.remainder*sum((left[i]*self.cap**i for i in range(3)),self.ctx.mpf(0))
                   +self.cap**3*self.remainder*other.remainder)
        return WidthPolynomial(rows,remainder,self.weight,self.cap)
    __rmul__=__mul__
    def packet(self):
        return dict(signed_hb_power_axial_coefficients=[list(row.coefficients) for row in self.rows],
            omitted_hb_cubed_weighted_axial_jet_norm_upper=self.remainder,
            remainder_axial_coefficient_upper_per_hb_cubed=[self.remainder/self.weight**k for k in range(self.rows[0].order+1)],
            exact_source_representation='sum_{j=0}^2 hb^j a_j(Z)+hb^3 e(Z); sum_k |e_k|*weight^k<=remainder',
            width_not_materialized=True)


def smoothing_weights(c,phase,subintervals=256):
    """A0=int alpha, A1=int u alpha; directed local Simpson remainders."""
    s=c.mpf(phase);lo,hi=endpoints(s)
    if lo!=hi or lo<0 or hi>2:raise ValueError('A single normalized phase in[0,2] required')
    if hi<=1:return dict(A0=s,A1=s*s/2,J=s*s/2,quadrature_error_bounds=[c.mpf(0),c.mpf(0)],cells=0)
    d=s-1;step=d/subintervals;integrals=[c.mpf(0),c.mpf(0)];errors=[c.mpf(0),c.mpf(0)]
    for n in range(subintervals):
        a=n*step;b=(n+1)*step;m=(a+b)/2
        vals=[1-sigma_jets(c,x)[0] for x in (a,m,b)]
        whole=sigma_jets(c,c.mpf([endpoints(a)[0],endpoints(b)[1]]))
        ordinary=[v*math.factorial(k) for k,v in enumerate(whole.coefficients)]
        f4=[abs(ordinary[4]),abs((1+c.mpf([endpoints(a)[0],endpoints(b)[1]]))*ordinary[4]+4*ordinary[3])]
        for j in range(2):
            samples=[(1+x)**j*v for x,v in zip((a,m,b),vals)]
            integrals[j]+=step*(samples[0]+4*samples[1]+samples[2])/6
            errors[j]+=step**5*upper(c,f4[j])/2880
    A0=1+integrals[0]+symmetric(c,errors[0]);A1=c.mpf('.5')+integrals[1]+symmetric(c,errors[1])
    if hi==2:A0=c.mpf('1.5') # sigma(1-u)=1-sigma(u): exact endpoint symmetry.
    return dict(A0=A0,A1=A1,J=s*A0-A1,quadrature_error_bounds=errors,cells=subintervals)


class CompliantComparisonPointIntegrals:
    def __init__(self):
        self.atoms=CompliantCoreIntegralAtoms();self.core=self.atoms.core;self.ctx=c=self.atoms.ctx
        self.bridge=CompliantBridgeMixedC4();self.source=self.bridge.formal_comparison_source()
        if self.bridge.source!=self.core.source or self.bridge.family!=self.core.family:
            raise ValueError('Comparison and core must share the original family')
        self.hashes=dict(self.atoms.hashes);self.hashes.update(self.bridge.hashes)
        name=PREFIX+'core_integral_atoms_check.json';receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['actual_core_atoms_axial6_available']:
            raise ValueError('Actual coefficientwise core atom receipt required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Actual inlet source changed: '+path)
        self.hashes.update(receipt['input_hashes'])
        for path in (name,Path(__file__).name):self.hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
        self.cap=c.mpf('1e-180');self.weight=self.atoms.axial_weight
        if endpoints(self.bridge.logh)[1]>=endpoints(c.ln(self.cap))[0]:raise ValueError('Exact hb not bounded by remainder cap')
        if endpoints(4*c.exp(2*self.cap))[1]>=endpoints(c.mpf('4.1'))[0]:raise ValueError('Comparison leaves original core continuation')
        self.cache={};self.weight_cache={}

    def weights(self,phase):
        key=str(phase)
        if key not in self.weight_cache:self.weight_cache[key]=smoothing_weights(self.ctx,phase)
        return self.weight_cache[key]

    def source_profile_jet(self,packet,rho,radial_order,root):
        c=self.ctx;r=c.mpf(rho);i=radial_order;N=packet['radial_degree'];order=6
        phi=[];uz=[];rmax=c.mpf(endpoints(r)[1])
        # Explicit factorial-model coefficients through Z6, radial0..3.
        if root:model=[c.mpf(0)]*7
        else:
            z=IntervalTaylor.variable(c,packet['Z'],order);one=IntervalTaylor.constant(c,1,order)
            H=z*((1-self.core.delta)/2)+(one-z*z)*(4*z+self.core.j)
            square=list((H*H).coefficients);square[0]=H[0]**2
            chi=IntervalTaylor(c,square)/(IntervalTaylor(c,square)+self.core.sigma**2)
            model=bessel_radial_tail_coefficients(c,chi,N,6,rmax,i,self.weight)
        for k in range(7):
            correction=self.core.correction*tail_factor(c,degree=N,radial_order=i,axial_order=k,radius=rmax,h=self.core.h)['tail_per_Xh_norm']/math.factorial(k)
            values=[]
            for name,tail in (('A',model[k]+correction),('Uz',self.core.epsilon*correction)):
                finite=sum((packet['rows'][name][n][k]*(math.factorial(n)//math.factorial(n-i))*r**(n-i)
                            for n in range(i,N+1)),c.mpf(0))
                values.append(finite+symmetric(c,tail))
            phi.append(values[0]);uz.append(values[1])
        if i==0:
            phi[0]=intersection(c,phi[0],c.mpf([endpoints(self.core.phi_floor)[0],endpoints(self.core.phi_ceiling)[1]]))
        return IntervalTaylor(c,phi),IntervalTaylor(c,uz)

    def inlet(self,Z,root=False):
        key=('root' if root else self.ctx.mpf(Z)._mpi_)
        if key in self.cache:return self.cache[key]
        packet=self.atoms.field.build_root_rows(24,6) if root else self.atoms.rebuild.rebuild(Z,24,6)
        atoms=self.atoms.atoms_from_packet(packet,6,shared_root=root)
        c=self.ctx;local=[self.source_profile_jet(packet,4,i,root) for i in range(3)]
        phi0,V0=local[0];p1,v1=local[1];p2,v2=local[2]
        L1=4*p1/phi0;L2=(4*p1+16*p2)/phi0-L1*L1
        V1=4*v1;V2=4*v1+16*v2
        cover=[self.source_profile_jet(packet,[4,'4.1'],i,root) for i in range(4)]
        r=c.mpf([4,'4.1']);p=[v[0] for v in cover];v=[v[1] for v in cover]
        d1=p[1]*r/p[0];d2=(p[1]*r+p[2]*r*r)/p[0]
        L3=(p[1]*r+p[2]*(3*r*r)+p[3]*r**3)/p[0]-3*d2*d1+2*d1*d1*d1
        V3=v[1]*r+v[2]*(3*r*r)+v[3]*r**3
        result=dict(Z=packet['Z'],phi0=phi0,V0=V0,L1=L1,L2=L2,V1=V1,V2=V2,
            L3_norm=jet_norm(L3,self.weight),V3_norm=jet_norm(V3,self.weight),
            moments={name:IntervalTaylor(c,values) for name,values in atoms['actual_core_atom_axial_coefficients'].items()},
            original_core_atoms=atoms,exact_root=root)
        self.cache[key]=result;return result

    def fields(self,data,phase,A0,A1):
        c=self.ctx;s=c.mpf(phase);p=data['phi0'];a=data['L1']*A0;b=data['L2']*A1
        na=jet_norm(a,self.weight);nb=jet_norm(b,self.weight);Rg=data['L3_norm']*s**3/6
        e=self.cap*na+self.cap**2*nb+self.cap**3*Rg
        if endpoints(e)[1]>mp.mpf('.5'):
            raise ArithmeticError('Weighted exponential remainder too large: '+mp.nstr(endpoints(e)[1],18))
        Rphi=jet_norm(p,self.weight)*(Rg*c.exp(e)+na*nb+self.cap*nb*nb/2
            +(na+self.cap*nb)**3*c.exp(self.cap*na+self.cap**2*nb)/6)
        phi=WidthPolynomial([p,p*a,p*(b+a*a/2)],Rphi,self.weight,self.cap)
        V=WidthPolynomial([data['V0'],data['V1']*A0,data['V2']*A1],data['V3_norm']*s**3/6,self.weight,self.cap)
        return phi,V

    def micro(self,Z,phase,root=False):
        c=self.ctx;s=c.mpf(phase);weights=self.weights(phase);data=self.inlet(Z,root)
        phi,V=self.fields(data,s,weights['A0'],weights['A1'])
        uniform_phi,uniform_V=self.fields(data,s,c.mpf([0,endpoints(s)[1]]),c.mpf([0,endpoints(s*s/2)[1]]))
        rhs=dict(H=2*uniform_phi,M=uniform_V,K=2*uniform_phi*uniform_V,A=uniform_V*uniform_V,
                 B=uniform_phi*uniform_phi,C=uniform_phi*uniform_phi)
        p,v,L1,V1=data['phi0'],data['V0'],data['L1'],data['V1']
        r0=dict(H=2*p,M=v,K=2*p*v,A=v*v,B=p*p,C=p*p)
        r1=dict(H=2*p*L1,M=V1,K=2*(p*L1*v+p*V1),A=2*v*V1,B=2*p*p*L1,C=2*p*p*L1)
        moments={}
        for name,initial in data['moments'].items():
            rate=2 if name in ('H','K','B') else 1
            first=r0[name]-rate*initial
            second=r1[name]*weights['J']-rate*first*s*s/2
            exponential=c.exp(rate*self.cap*s)
            bound=(jet_norm(initial,self.weight)*(rate*s)**3/6
                +rate**2*s**3*jet_norm(r0[name],self.weight)/6
                +rate*s*s*jet_norm(rhs[name].rows[1],self.weight)/2
                +s*(jet_norm(rhs[name].rows[2],self.weight)+self.cap*rhs[name].remainder))*exponential
            moments[name]=WidthPolynomial([initial,first*s,second],bound,self.weight,self.cap)
        return dict(data=data,phi=phi,V=V,moments=moments,phase=s,weights=weights)

    def decay(self,rate,fraction,Y):
        c=self.ctx;q=c.mpf(fraction);constant=c.exp(-rate*q*Y);a=2*rate*q
        row=IntervalTaylor.constant(c,constant,6)
        return WidthPolynomial([row,row*a,row*a*a/2],upper(c,abs(constant)*a**3*c.exp(a*self.cap)/6),self.weight,self.cap)

    def macro(self,Z,fraction,root=False):
        c=self.ctx;q=c.mpf(fraction)
        if endpoints(q)[0]<0 or endpoints(q)[1]>1:raise ValueError('Macro fraction in[0,1] required')
        end=self.micro(Z,2,root);Y=c.ln(100/self.bridge.r);phi,V=end['phi'],end['V']
        if endpoints(Y-2*self.cap)[0]<=0:raise ValueError('Frozen macro length not positive')
        moments={}
        targets=dict(H=phi,M=V,K=phi*V,A=V*V,B=phi*phi*self.ctx.mpf('.5'),C=phi*phi)
        for name,initial in end['moments'].items():
            rate=2 if name in ('H','K','B') else 1;decay=self.decay(rate,q,Y)
            moments[name]=decay*initial+(1-decay)*targets[name]
        return dict(data=end['data'],phi=phi,V=V,moments=moments,fraction=q,
            exact_log_radius='y=2hb+fraction*(log(100/Ra)-2hb)',macro_log_length_without_width=Y,
            exact_endpoint_at_R100=endpoints(q)==(mp.mpf(1),mp.mpf(1)),
            initial_smoothing_history_retained=True)

    def packet(self,result,chart):
        return dict(Z=result['data']['Z'],chart=chart,coordinate=result.get('phase',result.get('fraction')),
            comparison_phi=result['phi'].packet(),comparison_raw_V=result['V'].packet(),
            comparison_own_six_moments={name:value.packet() for name,value in result['moments'].items()},
            comparison_source_namespace=self.source['source_namespace'],
            axial_jet_weight=self.weight,width_remainder_cap=self.cap,original_positive_width_log_enclosure=self.bridge.logh,
            normalized_smoothing_weights=result.get('weights'),
            exact_log_radius=result.get('exact_log_radius','y=hb*s'),
            exact_endpoint_at_R100=result.get('exact_endpoint_at_R100',False),
            actual_core_atom_inlet_used=True,source_width_materialized=False,cap_used_as_width_value=False,
            actual_prescribed_shear_field_substituted_by_comparison=False,
            original_signed_actual_bridge_integrals_resolved=False,temporal_recursion=False)

    def report(self):
        points={}
        for z,root in (('.5',False),('0',False),('exact_shared_root',True)):
            points[z]=dict(micro=[self.packet(self.micro(0 if root else z,s,root),'micro') for s in (0,1,'1.5',2)],
                           macro=[self.packet(self.macro(0 if root else z,q,root),'macro') for q in (0,'.5',1)])
            print('Signed original comparison histories: '+z+', micro0..2 and frozen macro',flush=True)
        return dict(actual_five_defect_family_sha256=self.core.family,implicit_source_sha256=self.core.source,
            datum_enclosure_sha256=self.core.datum.datum_sha,original_comparison_source=self.source,comparison_point_packets=points,
            actual_core_atom_comparison_histories_numerically_integrated=True,
            original_signed_actual_bridge_integrals_resolved=False,all_annular_source_values_resolved=False,
            full_point_physical_field_evaluation=False,measured_blowup_dynamics=False,
            physical_energy_integral_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False,
            input_hashes=self.hashes)


def run():
    with mp.workdps(400):result=CompliantComparisonPointIntegrals().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

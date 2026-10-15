"""Full original n1 core through native axial Taylor order two.

The principal response is integrated through Z2. Every feedback word is
enclosed from that principal origin; no unavailable Z2 nominal iterate is
invented. Leading-source, finite kernel and feedback tails stay separate.
"""
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _box, _bounds
from .original_first_collar_jets import bessel_jet_tail
from .original_n1_core_volterra_numeric import OriginalN1CoreVolterraNumeric, DifferenceJet


def forcing_difference_z2(A,U,mean,z,r,delta,epsilon,P2,g):
    """Exact source map with genuine leading mixed coefficients through Z4."""
    end=lambda v:v.cut(0,2)
    d,L=1-z*z,1-z*z*delta
    def T(weight,value):
        i,k=value.radial-1,value.axial-1
        cut=lambda q:q.cut(i,k)
        weight=value.constant(weight) if not isinstance(weight,DifferenceJet) else weight
        return (cut(weight)*cut(z)*cut(value)+cut(d)*value.dz().cut(i,k)
                -2*cut(z)*cut(r)*value.dr().cut(i,k))/cut(L)
    a=-2-delta;m=d*g/L;Ta=T(a,A)
    f=-2*end(m)*end(m)*end(A)
    f+=(end(T(-3,m*A))+end(m)*end(Ta))*(2*end(epsilon))
    f-=end(T(-3,Ta))*(2*end(epsilon)**2)
    u=-2*end(T(-2,T(-1-delta,U)))
    cut=lambda v:v.cut(2,3)
    Q=(2*cut(z)*cut(U)-cut(z)*cut(mean)*cut(1-delta)-cut(d)*mean.dz().cut(2,3))/cut(L)
    q,qr,qrr,qz=end(Q),end(Q.dr()),end(Q.dr().dr()),end(Q.dz())
    zz,rr,dd,LL,ee,pp=map(end,(z,r,d,L,epsilon,P2))
    flow=(q+zz*qz*end(1-delta)/2+rr*qr)/LL+q*(q/2+rr*qr)
    flow+=end(U)*(-2*zz*q+dd*qz-2*zz*rr*qr)/LL
    p=-flow/(2*pp)+(2*qr+rr*qrr)/(ee*pp)
    u-=4*zz*rr*p*(ee*pp)/LL
    return dict(f=f,u=u,p=p)


def kernel_tail_z2(iv,chi,radius,degree):
    """Positive factorial-moment tails of the actual finite Kchi kernel."""
    j=degree+1;r=iv.mpf(radius);q=r/2;x=iv.sqrt(r)
    first=q**j/(math.factorial(j)*math.factorial(j+1))
    moments=[]
    for m in range(4):
        ratio=q*iv.mpf(j+1)**(m-1)/(iv.mpf(j)**m*(j+2))
        if _bounds(mp.mp,ratio)[1]>=1:
            raise ValueError('Kernel factorial-moment tail does not contract')
        moments.append(j**m*first/(1-ratio))
    c1,c2=abs(chi[1]),abs(chi[2])
    return dict(f=[x*moments[0]/2,x*c1*moments[1]/2,
                   x*(c2*moments[1]+c1*c1*moments[2]/2)/2],
                fx=[moments[1],c1*moments[2],c2*moments[2]+c1*c1*moments[3]/2])


class OriginalN1CoreAxial2(OriginalN1CoreVolterraNumeric):
    def __init__(self,*,precision=85,quadrature_nodes=192,directory=None):
        super().__init__(precision=precision,quadrature_nodes=quadrature_nodes,directory=directory)
        for stem in ('compliant_current_original_n1_core_volterra_numeric',
                     'compliant_current_original_first_collar_jets'):
            path=self.core.directory/(PREFIX+stem+'.json')
            report=json.loads(path.read_text(encoding='utf-8'))
            if not report['checks']['passed'] or report['source_family']!=self.family:
                raise ValueError('Accepted same-family source/kernel prerequisites required')
            for name,digest in report['implementation_hashes'].items():
                if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed axial-source prerequisite: '+name)
                self.code_hashes[name]=digest
            self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        relative='src/openai_ns_reconstruction/original_n1_core_axial2.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        self._axial2_cache={}

    def _real_chi(self,data,Z,order=4):
        iv,J=self.iv,self.leading.taylor
        z=J.variable(iv,iv.mpf(Z),order);de,j=data['de'],self.leading.j_iv
        H=z*((1-de)/2)+(1-z*z)*(4*z+j)
        h2=J(iv,[H[0]**2]+list((H*H).coefficients[1:]))
        chi=h2/(h2+(j/500)**2)
        lo,hi=_bounds(self.ctx,chi[0])
        return J(iv,[iv.mpf([max(0,lo),min(1,hi)])]+list(chi.coefficients[1:]))

    def _source_tail2(self,data,Z,radius,degree):
        c,R=self.iv,self.R;rmax=c.mpf(radius)
        whole=c.mpf([0,_bounds(self.ctx,rmax)[1]])
        h=_box(c,self.leading.tube['Xh_parameter'])
        correction=_box(c,self.leading.major['scaled_map_size_upper'])/(1-_box(c,self.leading.major['scaled_map_Lipschitz_upper']))
        chi=self._real_chi(data,Z,4)
        model=[bessel_jet_tail(c,chi,degree,radius=rmax,radial_order=i) for i in range(3)]
        def error(i,k,angular):
            bound=correction*self.leading.tail_factor(c,degree=degree,radial_order=i,
                axial_order=k,radius=rmax,h=h)['tail_per_Xh_norm']/math.factorial(k)
            bound=bound+model[i][k] if angular else bound*data['e']
            bound/=math.factorial(i)
            high=_bounds(self.ctx,bound)[1]
            return c.mpf([-high,high])
        def field(key,mean=False):
            rows=data['rows'][key]
            if any(len(row)<5 for row in rows):raise ValueError('Actual leading Z4 coefficients required')
            base=[[sum((rows[n][k]*math.comb(n,i)*whole**(n-i)/(n+1 if mean else 1)
                for n in range(i,degree+1)),c.mpf(0)) for k in range(5)] for i in range(3)]
            err=[[error(i,k,key=='A') for k in range(5)] for i in range(3)]
            return DifferenceJet(R(base,2,4,c.mpf),R(err,2,4,c.mpf))
        def fixed(rows):
            value=R(rows,2,4,c.mpf);return DifferenceJet(value,value.constant(0))
        if len(data['fixed']['ell_Z_taylor'])<5:raise ValueError('Actual fixed gradient Z4 required')
        z,r=fixed([[c.mpf(Z),1]]),fixed([[whole],[1]])
        g=fixed([[-v for v in data['fixed']['ell_Z_taylor'][:5]]])
        out=forcing_difference_z2(field('A'),field('Uz'),field('Uz',True),z,r,
            fixed([[data['de']]]),fixed([[data['e']]]),fixed([[data['P2']]]),g)
        errors={name:[c.mpf(_bounds(self.ctx,abs(q))[1]) for q in value.error.rows[0]] for name,value in out.items()}
        sizes={name:[c.mpf(_bounds(self.ctx,abs(q))[1]) for q in value.base.rows[0]] for name,value in out.items()}
        return errors,sizes,chi

    def _prepare2(self,Z,degree,kernel):
        key=Z,degree,kernel
        if key not in self._axial2_cache:
            source,point_data=self._forcing_polynomial(Z,degree,False)
            boxes,data=self._forcing_polynomial(Z,degree,True)
            principal=self._principal(source,self._chi(point_data,Z,2),kernel)
            chi=self._real_chi(data,Z,2)
            box=self._principal(boxes,self.R([list(chi.coefficients)],0,2,self.iv.mpf),kernel)
            self._axial2_cache[key]=principal,box,data
        return self._axial2_cache[key]

    def evaluate(self,Z,rho,*,degree=18,kernel_degree=24):
        c,iv=self.ctx,self.iv
        if (type(Z) is not str or type(rho) is not str
                or not all(c.isfinite(c.mpf(v)) for v in (Z,rho))
                or not -1<c.mpf(Z)<1 or not 0<=c.mpf(rho)<=4
                or type(degree) is not int or not 8<=degree<=21
                or type(kernel_degree) is not int or not 8<=kernel_degree<=40):
            raise ValueError('Interior Z, rho0..4, source degree8..21 and kernel degree8..40 required')
        if c.mpf(rho)==0:
            zero=self._summary(c.mpf(0),iv.mpf(0))
            return dict(Z=Z,rho=rho,source_family=self.family,axial_Taylor_order=2,q_axial_Taylor_order=1,
                units='coefficient k = D_Z^k/k!; multiply value and error by 2 for ordinary Z2',
                normalized_original_n1={**{name:[dict(zero) for _ in range(3)] for name in self.fields},
                                       'mean_u':[dict(zero) for _ in range(3)],'q':[dict(zero) for _ in range(2)]},
                exact_regular_axis=True,full_n1_core_Z2_absolute_errors_available=True,
                core_only=True,q_Z2=False,collar_Z2_transport=False,temporal_n_ge_2_recursion=False)
        principal,box,data=self._prepare2(Z,degree,kernel_degree)
        point,boxes=self._values(principal,rho),self._values(box,rho)
        errors,sizes,chi=self._source_tail2(data,Z,rho,degree)
        r,x=iv.mpf(rho),iv.sqrt(iv.mpf(rho));q=r/2
        E=iv.exp(q);a=q*abs(chi[1]);b=q*abs(chi[2])
        expjet=[iv.mpf(1),a,b+a*a/2]
        conv=lambda left,right,k:sum((left[i]*right[k-i] for i in range(k+1)),iv.mpf(0))
        source=dict(f=[r*E*conv(expjet,errors['f'],k)/8 for k in range(3)],
            fx=[x*E*conv(expjet,errors['f'],k)/4 for k in range(3)],
            u=[r*v/4 for v in errors['u']],k=[r*v/8 for v in errors['u']],
            ux=[x*v/2 for v in errors['u']],p=[r*v for v in errors['p']])
        kt=kernel_tail_z2(iv,chi,rho,kernel_degree)
        kernel={name:[x*conv(kt[name],sizes['f'],k) for k in range(3)] for name in ('f','fx')}
        tails=[_box(iv,self.tail['value_remainder_upper']),_box(iv,self.tail['ordinary_Z1_remainder_upper']),
               _box(iv,self.tail['ordinary_Z2_remainder_upper'])/2]
        truth,ledger={},{}
        for name in self.fields:
            scale=x if name in ('fx','ux') else r
            truth[name],ledger[name]=[],[]
            for k,tail in enumerate(tails):
                feedback=scale*tail;kb=kernel[name][k] if name in kernel else iv.mpf(0)
                high=_bounds(c,source[name][k]+kb+feedback)[1]
                truth[name].append(boxes[name][k]+iv.mpf([-high,high]))
                ledger[name].append(dict(leading_source_tail=c.nstr(_bounds(c,source[name][k])[1],40),
                    finite_principal_kernel_tail=c.nstr(_bounds(c,kb)[1],40),
                    all_feedback_word_tail=c.nstr(_bounds(c,feedback)[1],40)))
        def recover(values,ctx,de):
            J=self.leading.taylor;z=J.variable(ctx,ctx.mpf(Z),2)
            u,k=J(ctx,values['u']),J(ctx,values['k']);mean=u+k
            mz=J(ctx,[(i+1)*mean[i+1] for i in range(2)])
            qr=(u*(2*z)-mean*(z*(1+de))-mz*(1-z*z))/(1-z*z*de)
            return list(mean.coefficients),list(qr.coefficients)
        pm,pq=recover(point,c,self.leading.delta);tm,tq=recover(truth,iv,self.leading.delta_iv)
        output={name:[self._summary(point[name][k],truth[name][k]) for k in range(3)] for name in self.fields}
        output['mean_u']=[self._summary(a,b) for a,b in zip(pm,tm)]
        output['q']=[self._summary(a,b) for a,b in zip(pq,tq)]
        return dict(Z=Z,rho=rho,source_family=self.family,axial_Taylor_order=2,q_axial_Taylor_order=1,
            units='coefficient k = D_Z^k/k!; multiply value and error by 2 for ordinary Z2',
            normalized_original_n1=output,error_ledger=ledger,leading_radial_degree=degree,
            principal_kernel_degree=kernel_degree,nominal_feedback_iterations=0,
            all_feedback_words_enclosed_from_principal_origin=True,
            finite_forcing_polynomial_integrated_analytically=True,
            full_n1_core_Z2_absolute_errors_available=True,
            limitations=dict(core_only=True,nominal_feedback_correction_computed=False,
                q_Z2=False,collar_Z2_transport=False,physical_mixed_derivatives=False,
                global_moment_repair=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False))

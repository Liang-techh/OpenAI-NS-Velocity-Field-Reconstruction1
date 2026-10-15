"""Original n1 core values from principal Volterra integrals and all-word bounds.

Finite forcing polynomials are integrated analytically, not interpreted as
a Taylor truncation of the unknown n1 solution. The error budget separately
covers the admitted leading-core tail, Kchi kernel tail and all feedback.
Nominal feedback iterates are computed, but the rigorous feedback enclosure
uses the accepted index-zero word tail, conservatively charging the nominal
correction against that same principal origin.
"""
import gzip
import hashlib
import json
import math

from .original_core_numeric import PREFIX, ROOT, _box, _bounds
from .original_n1_core_numeric import OriginalN1CoreNumeric, linear_rhs


class DifferenceJet:
    """A finite jet and an enclosure of its change, with exact product rule."""
    def __init__(self,base,error):
        self.base,self.error = base,error
        self.radial,self.axial = base.radial,base.axial

    def constant(self,value):
        return DifferenceJet(self.base.constant(value),self.base.constant(0))

    def coerce(self,value):
        return value if isinstance(value,DifferenceJet) else self.constant(value)

    def __add__(self,other):
        other = self.coerce(other)
        return DifferenceJet(self.base+other.base,self.error+other.error)

    __radd__ = __add__
    def __neg__(self):return DifferenceJet(-self.base,-self.error)
    def __sub__(self,other):return self+-self.coerce(other)
    def __rsub__(self,other):return self.coerce(other)+-self
    def __mul__(self,other):
        other = self.coerce(other)
        return DifferenceJet(self.base*other.base,
            self.base*other.error+self.error*other.base+self.error*other.error)
    __rmul__ = __mul__
    def __truediv__(self,other):
        other = self.coerce(other)
        if any(v != 0 for row in other.error.rows for v in row):
            raise ValueError("Only unchanged original denominators may be inverted")
        return DifferenceJet(self.base/other.base,self.error/other.base)
    def __pow__(self,power):
        if type(power) is not int or power < 0:raise ValueError("Nonnegative integer jet power required")
        result = self.constant(1)
        for _ in range(power):result = result*self
        return result
    def cut(self,radial,axial):
        def crop(v):
            if radial > v.radial or axial > v.axial:raise ValueError("Missing mixed derivative")
            return type(v)(v.rows,radial,axial,v.convert)
        return DifferenceJet(crop(self.base),crop(self.error))
    def dz(self):return DifferenceJet(self.base.dz(),self.error.dz())
    def dr(self):
        def diff(v):
            if v.radial < 1:raise ValueError("Missing radial derivative")
            return type(v)([[(i+1)*q for q in v.rows[i+1]] for i in range(v.radial)],
                v.radial-1,v.axial,v.convert)
        return DifferenceJet(diff(self.base),diff(self.error))


def forcing_difference(A,U,mean,z,r,delta,epsilon,P2,g):
    """Exact source map on mixed jets, propagating leading-tail differences."""
    end = lambda v:v.cut(0,1)
    d,L = 1-z*z,1-z*z*delta
    def T(weight,value):
        i,k = value.radial-1,value.axial-1
        cut = lambda v:v.cut(i,k)
        weight = value.constant(weight) if not isinstance(weight,DifferenceJet) else weight
        return (cut(weight)*cut(z)*cut(value)+cut(d)*value.dz().cut(i,k)
            -2*cut(z)*cut(r)*value.dr().cut(i,k))/cut(L)
    a = -2-delta
    m = d*g/L
    Ta = T(a,A)
    f = -2*end(m)*end(m)*end(A)
    f += (end(T(-3,m*A))+end(m)*end(Ta))*(2*end(epsilon))
    f -= end(T(-3,Ta))*(2*end(epsilon)**2)
    u = -2*end(T(-2,T(-1-delta,U)))
    cut = lambda v:v.cut(2,2)
    Q = (2*cut(z)*cut(U)-cut(z)*cut(mean)*cut(1-delta)-cut(d)*mean.dz())/cut(L)
    q,qr,qrr,qz = end(Q),end(Q.dr()),end(Q.dr().dr()),end(Q.dz())
    zz,rr,dd,LL,ee,pp = map(end,(z,r,d,L,epsilon,P2))
    S = (q+zz*qz*end(1-delta)/2+rr*qr)/LL+q*(q/2+rr*qr)
    S += end(U)*(-2*zz*q+dd*qz-2*zz*rr*qr)/LL
    p = -S/(2*pp)+(2*qr+rr*qrr)/(ee*pp)
    u -= 4*zz*rr*p*(ee*pp)/LL
    return dict(f=f,u=u,p=p)


class OriginalN1CoreVolterraNumeric:
    """Evaluated normalized original n1 core with a complete absolute error."""
    fields = ("f","u","k","p","fx","ux")
    def __init__(self,*,precision=60,quadrature_nodes=128,directory=None):
        self.core = OriginalN1CoreNumeric(precision=precision,quadrature_nodes=quadrature_nodes,directory=directory)
        self.leading = self.core.leading
        self.ctx,self.iv,self.R = self.core.ctx,self.core.iv,self.core.Rectangle
        self.precision,self.family = precision,self.core.family
        self.hashes,self.code_hashes = dict(self.core.hashes),dict(self.core.code_hashes)
        self.code_hashes['src/openai_ns_reconstruction/original_n1_core_numeric.py'] = hashlib.sha256(
            (ROOT/'src/openai_ns_reconstruction/original_n1_core_numeric.py').read_bytes()).hexdigest()
        stem = PREFIX+'compliant_current_original_n1_regular_series_function'
        receipt = json.loads((self.core.directory/(stem+'_check.json')).read_text(encoding='utf-8'))
        archive = self.core.directory/(stem+'.json.gz')
        if not receipt['all_passed'] or receipt['source_family'] != self.family:
            raise ValueError("Accepted same-family regular all-word bound required")
        for name,digest in receipt['input_hashes'].items():
            path = self.core.directory/name
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError("Original all-word source changed: "+name)
            self.hashes[name] = digest
        for name in (stem+'_check.json',stem+'.json.gz'):
            self.hashes[name] = hashlib.sha256((self.core.directory/name).read_bytes()).hexdigest()
        saved = json.loads(gzip.decompress(archive.read_bytes()))
        self.contract = saved['contract']
        if self.contract['finite_word_index'] != 0 or saved['source_family'] != self.family:
            raise ValueError("Index-zero principal-origin bound required")
        self.tail = self.contract['regular_series_bound']
        if self.tail['regular_axis_powers'] != [2,2,2,2,1,1]:
            raise ValueError("Natural six-state regular norm required")
        self.chi_certificate = self.tail['original_holomorphic_chi_tube']
        self._forcing_cache = {}
        self._solution_cache = {}

    def _leading(self,Z,degree,directed):
        c = self.iv if directed else self.ctx
        # One extra leading radial row leaves Z^4 in every retained row.
        fixed,rows,de,e = self.leading._seed(Z,degree+1,3,directed=directed)
        for n in range(degree+1):self.leading.advance(c,fixed,rows,n,c.mpf(Z),de,e)
        rows = {key:value[:degree+1] for key,value in rows.items()}
        lp = self.leading.lp_iv if directed else self.leading.logP
        return c,fixed,rows,de,e,c.exp(2*lp)

    def _forcing_polynomial(self,Z,degree,directed):
        cache_key = Z,degree,directed
        if cache_key in self._forcing_cache:return self._forcing_cache[cache_key]
        c,fixed,rows,de,e,P2 = self._leading(Z,degree,directed)
        radial,axial = 2*degree+1,2
        # This extension defines the finite leading *approximant*. Its
        # omitted true source is bounded separately by forcing_difference.
        leading = {key:[list(row[:axial+3]) for row in rows[key]]+
            [[c.mpf(0)]*(axial+3) for _ in range(radial+2-len(rows[key]))]
            for key in ('A','Uz')}
        zero = {key:[[c.mpf(0)]*(axial+2) for _ in range(radial+1)] for key in ('f','u','k','p')}
        rf,ru,rp = self.core._rhs(fixed,leading,zero,radial,axial,c.mpf(Z),de,e,P2,c)
        result = dict(f=rf,u=ru,p=rp),dict(c=c,fixed=fixed,rows=rows,de=de,e=e,P2=P2)
        self._forcing_cache[cache_key] = result
        return result

    def _chi(self,data,Z,order):
        c,de = data['c'],data['de']
        J = self.leading.taylor
        z = J.variable(c,c.mpf(Z),order)
        j = self.leading.j_iv if c is self.iv else self.leading.j
        H = z*((1-de)/2)+(1-z*z)*(4*z+j)
        chi = H*H/(H*H+(j/500)**2)
        return self.R([list(chi.coefficients)],0,order,c.mpf)

    def _principal(self,source,chi,kernel_degree):
        R,c = self.R,chi.convert
        n,K = source['f'].radial,source['f'].axial
        def integrate(value,field):
            rows = [[c(0)]*(K+1)]
            for i,row in enumerate(value.rows):
                factor = 4*(i+1)*(i+2) if field=='f' else 4*(i+1)**2 if field=='u' else i+1
                rows.append([v/factor for v in row])
            return R(rows,value.radial+1,K,c)
        angular = integrate(source['f'],'f')
        terms = [angular]
        for _ in range(kernel_degree):
            ch = R(chi.rows,angular.radial,K,c)
            angular = integrate(-2*ch*angular,'f')
            terms.append(angular)
        radial = terms[-1].radial
        f = sum((R(term.rows,radial,K,c) for term in terms),R([[0]],radial,K,c))
        u = integrate(source['u'],'u')
        p = integrate(source['p'],'p')
        u,p = R(u.rows,radial,K,c),R(p.rows,radial,K,c)
        k = R([[-c(i)*v/(i+1) for v in row] for i,row in enumerate(u.rows)],radial,K,c)
        return dict(f=f,u=u,k=k,p=p)

    def _residual_action(self,state,data,Z,chi):
        """Actual C0+B1*dZ polynomial action, with -2chi removed stably."""
        c,R = data['c'],self.R
        K = state['f'].axial-1
        N = len(data['rows']['A'])+state['f'].radial
        def ring(rows,order=K):return R(rows,N,order,c.mpf)
        z = ring([[c.mpf(Z),1]])
        de,e,P2 = (z.constant(data[key]) for key in ('de','e','P2'))
        d,L = 1-z*z,1-z*z*de
        A,U = (ring(data['rows'][key]) for key in ('A','Uz'))
        uz = ring([[(j+1)*row[j+1] for j in range(K+1)] for row in data['rows']['Uz']])
        az = ring([[(j+1)*row[j+1] for j in range(K+1)] for row in data['rows']['A']])
        mean = U.radial_average()
        meanz = uz.radial_average()
        Q = (2*z*U-z*mean*(1-de)-d*meanz)/L
        g = -ring([data['fixed']['ell_Z_taylor']])
        F2 = ring([data['fixed']['S_Z_taylor']])/(e*e)
        rho = ring([[0],[1]])
        v = {key:ring(value.rows) for key,value in state.items()}
        v.update({key+'Z':ring(value.dz().rows) for key,value in state.items()})
        f = v['f']
        # Cancel the principal term by the exact U0 identity, not by
        # subtracting near-equal floating point O(1) quantities.
        udiff = ring([[0]*(K+1)]+data['rows']['Uz'][1:])
        v.update(z=z,rho=rho,delta=de,epsilon=e,P2=P2,F2=F2,L=L,d=d,A=A,C=A+A.radial_euler(),g=g,
            E=z*U-c.mpf(1)/2,H=z*(1-de)/2+d*U,B=Q+(1-2*z*U)/L,Q=Q,
            Zf_without_gauge=((-2-de)*z*A+d*az-2*z*A.radial_euler())/L,
            Zw=((-1-de)*z*U+d*uz-2*z*U.radial_euler())/L,
            Euler_f=f.radial_euler(),Euler_u=v['u'].radial_euler(),Euler_U=U.radial_euler())
        v['f'] = z.constant(0)
        rf,ru,rp = linear_rhs(v)
        rf += (2*e*(Q+v['E']*(-2+de)/L)-2*d*udiff*g/L)*f
        ru -= 8*F2*z*rho*A*f/L
        rp += 2*F2/(e*P2)*A*f
        return dict(f=rf,u=ru,p=rp)

    def _source_tail(self,data,Z,radius,degree):
        """Difference of true forcing and the finite leading approximant."""
        c,R = self.iv,self.R
        rmax = c.mpf(radius)
        whole = c.mpf([0,_bounds(self.ctx,rmax)[1]])
        h = _box(c,self.leading.tube['Xh_parameter'])
        eta = _box(c,self.leading.tube['complex_tube_radius'])
        # Certify |chi|<=1 on this very leading-model complex tube before
        # using Cauchy on chi^n. This is not an n1 complex-radial theorem.
        j,de = self.leading.j_iv,data['de']
        b = eta*(12*(1+eta)**2+4*eta**2+2*j*(1+eta)+4+abs(1-de)/2)
        if _bounds(self.ctx,(j/500)**2-2*b*b)[0] <= 0:
            raise ValueError("Original leading chi tube not admitted")
        correction = _box(c,self.leading.major['scaled_map_size_upper'])/(1-_box(c,self.leading.major['scaled_map_Lipschitz_upper']))
        def error(i,k,angular):
            bound = correction*self.leading.tail_factor(c,degree=degree,radial_order=i,axial_order=k,
                radius=rmax,h=h)['tail_per_Xh_norm']
            if angular:
                n = degree+1
                first = rmax**(n-i)/(2**n*math.factorial(n-i)*math.factorial(n+1))
                ratio = rmax/(2*(n+1-i)*(n+2))
                bound += math.factorial(k)/(eta/2)**k*first/(1-ratio)
            else:bound *= data['e']
            high = _bounds(self.ctx,bound/(math.factorial(i)*math.factorial(k)))[1]
            return c.mpf([-high,high])
        def field(key,mean=False):
            rows = data['rows'][key]
            base = [[sum((rows[n][k]*math.comb(n,i)*whole**(n-i)/(n+1 if mean else 1)
                for n in range(i,degree+1)),c.mpf(0)) for k in range(4)] for i in range(3)]
            err = [[error(i,k,key=='A') for k in range(4)] for i in range(3)]
            return DifferenceJet(R(base,2,3,c.mpf),R(err,2,3,c.mpf))
        def fixed(rows):
            base = R(rows,2,3,c.mpf)
            return DifferenceJet(base,base.constant(0))
        z,r = fixed([[c.mpf(Z),1]]),fixed([[whole],[1]])
        g = fixed([[-v for v in data['fixed']['ell_Z_taylor'][:4]]])
        out = forcing_difference(field('A'),field('Uz'),field('Uz',True),z,r,
            fixed([[data['de']]]),fixed([[data['e']]]),fixed([[data['P2']]]),g)
        errors = {key:[c.mpf(_bounds(self.ctx,abs(v))[1]) for v in value.error.rows[0]] for key,value in out.items()}
        sizes = {key:[c.mpf(_bounds(self.ctx,abs(v))[1]) for v in value.base.rows[0]] for key,value in out.items()}
        return errors,sizes

    def _values(self,state,radius):
        c = state['f'].convert
        r = c(radius)
        result = {}
        for key,value in state.items():
            result[key] = [sum((row[k]*r**n for n,row in enumerate(value.rows)),c(0)) for k in range(value.axial+1)]
        for key,name in (('f','fx'),('u','ux')):
            result[name] = [sum((2*r**(c(1)/2)*n*row[k]*r**(n-1)
                for n,row in enumerate(state[key].rows) if n),c(0)) for k in range(state[key].axial+1)]
        return result

    def _prepare(self,Z,degree,kernel_degree,feedback_iterations):
        key = Z,degree,kernel_degree,feedback_iterations
        if key in self._solution_cache:return self._solution_cache[key]
        c = self.ctx
        source,nominal = self._forcing_polynomial(Z,degree,False)
        chi = self._chi(nominal,Z,2)
        principal = self._principal(source,chi,kernel_degree)
        increments = []
        if feedback_iterations:
            residual = self._residual_action(principal,nominal,Z,chi)
            increments.append(self._principal(residual,chi.truncate_axial(1),kernel_degree))
        N = max(v['f'].radial for v in [principal]+increments)
        total = {key:self.R(principal[key].rows,N,1,c.mpf) for key in ('f','u','k','p')}
        for increment in increments:
            for field in total:total[field] += self.R(increment[field].rows,N,1,c.mpf)
        source_box,directed = self._forcing_polynomial(Z,degree,True)
        chii = self._chi(directed,Z,2)
        principal_box = self._principal(source_box,chii,kernel_degree)
        result = total,increments,directed,chii,principal_box
        self._solution_cache[key] = result
        return result

    def evaluate(self,Z,rho,*,degree=16,kernel_degree=20,feedback_iterations=1):
        c,iv = self.ctx,self.iv
        if (type(Z) is not str or not -1<c.mpf(Z)<1 or type(rho) is not str or not 0<=c.mpf(rho)<=4
                or type(degree) is not int or not 4<=degree<=21 or type(kernel_degree) is not int
                or not 4<=kernel_degree<=32 or feedback_iterations not in (0,1)):
            raise ValueError("Interior decimal Z, rho in [0,4], degree4..21, kernel4..32, iterations0/1 required")
        if c.mpf(rho) == 0:
            zero = self._summary(c.mpf(0),iv.mpf(0))
            return dict(Z=Z,rho=rho,source_family=self.family,
                normalized_original_n1={**{key:[dict(zero),dict(zero)] for key in self.fields},
                    'mean_u':dict(zero),'q':dict(zero)},exact_regular_axis=True,
                original_full_core_n1_absolute_error_available=True,global_matched_n1=False,
                n0_plus_n1_physical_field=False,temporal_n_ge_2_recursion=False,full_corrected_NS_solution=False)
        total,increments,directed,chii,principal_box = self._prepare(Z,degree,kernel_degree,feedback_iterations)
        point = self._values(total,rho)
        nominal_increments = [{key:[c.nstr(v,self.precision) for v in values]
            for key,values in self._values(increment,rho).items()} for increment in increments]
        boxes = self._values(principal_box,rho)
        errors,sizes = self._source_tail(directed,Z,rho,degree)
        r,x = iv.mpf(rho),iv.sqrt(iv.mpf(rho))
        chiq = abs(chii.rows[0][1])
        expq = iv.exp(r/2)
        extra = dict(f=[r*expq*errors['f'][0]/8,
            r*expq*(errors['f'][1]+r*chiq*errors['f'][0]/2)/8],
            fx=[x*expq*errors['f'][0]/4,
                x*expq*(errors['f'][1]+r*chiq*errors['f'][0]/2)/4],
            u=[r*v/4 for v in errors['u']],k=[r*v/8 for v in errors['u']],
            ux=[x*v/2 for v in errors['u']],p=[r*v for v in errors['p']])
        j,q = kernel_degree+1,r/2
        term = q**j/(math.factorial(j)*math.factorial(j+1))
        q0,q1,q2 = q/((j+1)*(j+2)),q/(j*(j+2)),q*(j+1)/(j*j*(j+2))
        kt = dict(f=[x*term/(2*(1-q0)),x*chiq*j*term/(2*(1-q1))],
            fx=[j*term/(1-q1),chiq*j*j*term/(1-q2)])
        for key in ('f','fx'):
            extra[key][0] += kt[key][0]*x*sizes['f'][0]
            extra[key][1] += x*(kt[key][1]*sizes['f'][0]+kt[key][0]*sizes['f'][1])
        tail0 = _box(iv,self.tail['value_remainder_upper'])
        tail1 = _box(iv,self.tail['ordinary_Z1_remainder_upper'])
        truth = {}
        ledger = {}
        for key in self.fields:
            scale = x if key in ('fx','ux') else r
            truth[key],ledger[key] = [],[]
            for k,tail in enumerate((tail0,tail1)):
                feedback_error = tail*scale
                bound = extra[key][k]+feedback_error
                high = _bounds(c,bound)[1]
                truth[key].append(boxes[key][k]+iv.mpf([-high,high]))
                ledger[key].append(dict(source_and_kernel_tail_upper=c.nstr(_bounds(c,extra[key][k])[1],35),
                    full_feedback_word_tail_upper=c.nstr(_bounds(c,feedback_error)[1],35)))
        z,zi = c.mpf(Z),iv.mpf(Z)
        def recovery(values,zz,de):
            mean = values['u'][0]+values['k'][0]
            meanZ = values['u'][1]+values['k'][1]
            return mean,((1-de)*zz*values['u'][0]-(1+de)*zz*values['k'][0]-(1-zz*zz)*meanZ)/(1-de*zz*zz)
        pm,pq = recovery(point,z,self.leading.delta)
        tm,tq = recovery(truth,zi,self.leading.delta_iv)
        output = {}
        for key in self.fields:
            output[key] = [self._summary(point[key][k],truth[key][k]) for k in (0,1)]
        output['mean_u'] = self._summary(pm,tm)
        output['q'] = self._summary(pq,tq)
        return dict(Z=Z,rho=rho,leading_radial_degree=degree,potential_kernel_degree=kernel_degree,
            nominal_feedback_iterations=feedback_iterations,normalized_original_n1=output,
            nominal_feedback_increments=nominal_increments,
            error_ledger=ledger,accepted_all_word_ratio_upper=self.tail['ratio_upper']['upper'],
            source_family=self.family,original_full_core_n1_absolute_error_available=True,
            point_proxy_is_not_an_exact_source_oracle=True,
            all_feedback_bounded_from_principal_origin=True,
            nominal_feedback_correction_charged_in_distance_to_principal_enclosure=True,
            finite_forcing_polynomial_integrated_analytically=True,
            n1_radial_Taylor_tail_substituted=False,global_matched_n1=False,
            n0_plus_n1_physical_field=False,temporal_n_ge_2_recursion=False,full_corrected_NS_solution=False)

    def _summary(self,point,truth):
        c,iv = self.ctx,self.iv
        text = c.nstr(point,self.precision)
        error = max(abs(v) for v in _bounds(c,truth-iv.mpf(text)))
        return dict(approximation=text,absolute_error_upper=c.nstr(error*(1+c.mpf('1e-30')),40),
            relative_error_upper=None if not point else c.nstr(error/abs(point)*(1+c.mpf('1e-30')),35),
            error_applies_to_serialized_decimal=True,original_source_enclosure_exact_mpf_tuples=[list(v) for v in truth._mpi_])

"""Original O3 scalar integrals, with the tiny mu correction kept separate.

U_c=expm1(-c*mu*J)/mu stays order one. Its weighted integrals
are transported before multiplication by mu; no subtraction from one is used.
"""
from fractions import Fraction
import math

from .original_core_numeric import _bounds


class OriginalO3Scalar:
    names=('J','U1','U2','C1','C0','Cm')

    def __init__(self,ctx,iv,mu_point,mu_directed,*,cells=512):
        from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
        if type(cells) is not int or cells<16:raise ValueError('At least16 O3 scalar cells required')
        if not 0<_bounds(ctx,mu_directed)[0]<=_bounds(ctx,mu_directed)[1]<=ctx.mpf('.001'):
            raise ValueError('Original positive mu<=.001 required')
        self.ctx,self.iv,self.mu,self.mui=ctx,iv,mu_point,mu_directed
        self.cells,self.sigma_jets=cells,sigma_jets
        self.grid=[([ctx.mpf(0)]*6,[iv.mpf(0)]*6)];self.cache={};self.steps=0

    def _rational(self,q,ctx):return ctx.mpf(q.numerator)/q.denominator

    def _intersection(self,a,b):
        c,iv=self.ctx,self.iv;lo,hi=_bounds(c,a);ll,hh=_bounds(c,b)
        if max(lo,ll)>min(hi,hh):raise ArithmeticError('O3 scalar enclosure misses its source tube')
        return iv.mpf([max(lo,ll),min(hi,hh)])

    def _jets(self,ctx,state,sigma,exponentials,mu):
        exponential,negative_exponential=exponentials
        values=[[v] for v in state]
        for n in range(5):
            values[0].append(sigma[n]/(n+1))
            for index,c in ((1,1),(2,2)):
                rhs=sum((sigma[k]*((1 if n-k==0 else 0)+mu*values[index][n-k]) for k in range(n+1)),ctx.mpf(0))
                values[index].append(-c*rhs/(n+1))
            values[3].append(sum((exponential[k]*values[1][n-k] for k in range(n+1)),ctx.mpf(0))/(n+1))
            values[4].append(values[2][n]/(n+1))
            values[5].append(sum((negative_exponential[k]*values[2][n-k]
                for k in range(n+1)),ctx.mpf(0))/(n+1))
        return values

    def _step(self,left,width,state):
        c,iv=self.ctx,self.iv;p,t=state;h=self._rational(width,iv);hp=self._rational(width,c)
        x=self._rational(left,iv);box=iv.mpf([_bounds(c,x)[0],_bounds(c,x+h)[1]])
        tube=[]
        for i,cap in enumerate((1,1,2,_bounds(c,iv.exp(1))[1],2,2)):
            lo,hi=_bounds(c,t[i])
            if i==0:tube.append(iv.mpf([max(0,lo),min(1,_bounds(c,t[i]+h)[1])]))
            else:tube.append(iv.mpf([max(-cap,lo-_bounds(c,h*cap)[1]) if i<3 else lo-_bounds(c,h*cap)[1],min(0,hi)]))
        def coeff(q):
            return self.sigma_jets(iv,q),([iv.exp(q)/math.factorial(k) for k in range(5)],
                [iv.exp(-q)*(-1)**k/math.factorial(k) for k in range(5)])
        at,ex=coeff(x);whole,wx=coeff(box)
        mid=lambda v:sum(_bounds(c,v),c.mpf(0))/2
        pj=self._jets(c,p,list(map(mid,at)),tuple(list(map(mid,v)) for v in ex),self.mu)
        tj=self._jets(iv,t,at,ex,self.mui);rem=self._jets(iv,tube,whole,wx,self.mui)
        pp=[sum((v[k]*hp**k for k in range(5)),c.mpf(0)) for v in pj]
        tt=[sum((v[k]*h**k for k in range(5)),iv.mpf(0))+w[5]*h**5 for v,w in zip(tj,rem)]
        tt=[self._intersection(v,w) for v,w in zip(tt,tube)]
        # The defining expm1 identity supplies a second, correlated enclosure.
        for index,factor in ((1,1),(2,2)):
            tt[index]=self._intersection(tt[index],iv.expm1(-factor*self.mui*tt[0])/self.mui)
        self.steps+=1;return pp,tt

    def transition(self,value):
        if type(value) is not str:raise ValueError('Exact decimal O3 coordinate required')
        q=Fraction(value)
        if not 0<=q<=1:raise ValueError('Original O3 transition0<=x<=1 required')
        if q in self.cache:return self.cache[q]
        index=(q*self.cells).numerator//(q*self.cells).denominator
        while len(self.grid)<=index:
            left=Fraction(len(self.grid)-1,self.cells)
            self.grid.append(self._step(left,Fraction(1,self.cells),self.grid[-1]))
        p,t=self.grid[index];left=Fraction(index,self.cells)
        if q!=left:p,t=self._step(left,q-left,(p,t))
        p,t=list(p),list(t)
        if q>0:
            # The original flat sigmoid is monotone and strictly positive
            # inside (0,1). These integral comparisons retain signs even
            # when a Taylor remainder dominates the tiny flat-end value.
            iv,c=self.iv,self.ctx;x=self._rational(q,iv)
            sx=self.sigma_jets(iv,x)[0];sh=self.sigma_jets(iv,x/2)[0]
            t[0]=self._intersection(t[0],iv.mpf([_bounds(c,x*sh/2)[0],_bounds(c,x*sx)[1]]))
            U=[iv.expm1(-factor*self.mui*t[0])/self.mui for factor in (1,2)]
            at_threequarters=[iv.expm1(-factor*self.mui*(x*sh/4))/self.mui for factor in (1,2)]
            for index,lo,hi in ((3,x*iv.exp(x)*U[0],x*iv.exp(3*x/4)*at_threequarters[0]/4),
                    (4,x*U[1],x*at_threequarters[1]/4),(5,x*U[1],x*iv.exp(-x)*at_threequarters[1]/4)):
                t[index]=self._intersection(t[index],iv.mpf([_bounds(c,lo)[0],_bounds(c,hi)[1]]))
        if q==1:
            t[0]=self._intersection(t[0],self.iv.mpf('.5'));p[0]=self.ctx.mpf('.5')
        for index,factor in ((1,1),(2,2)):
            t[index]=self._intersection(t[index],self.iv.expm1(-factor*self.mui*t[0])/self.mui)
            p[index]=self.ctx.expm1(-factor*self.mu*p[0])/self.mu
        result=dict(point=dict(zip(self.names,p)),directed=dict(zip(self.names,t)),
            method='Taylor4 with ODE fifth-coefficient remainder, monotone sigmoid integral comparison; separate normalized expm1 corrections',
            uniform_cells=self.cells,endpoint_sigma_integral_half=q==1)
        self.cache[q]=result;return result

    def _moment(self,x,alpha,n,ctx):
        # Positive integral int_0^x exp(alpha*u)*(.5+u)^n du. Taylor
        # remainder of exp(alpha*u) is bounded over the whole [0,x].
        total=ctx.mpf(0);tail=ctx.mpf(0);length=100
        for k in range(n+1):
            coefficient=math.comb(n,k)*ctx.mpf('.5')**(n-k)
            series=sum(((alpha*x)**j/(math.factorial(j)*(k+j+1)) for j in range(length)),ctx.mpf(0))
            total+=coefficient*x**(k+1)*series
            tail+=coefficient*ctx.exp(abs(alpha)*x)*abs(alpha*x)**length*x**(k+1)/(math.factorial(length)*(k+length+1))
        if ctx is self.iv:
            upper=_bounds(self.ctx,tail)[1];total+=ctx.mpf([-upper,upper])
        return total

    def power(self,value):
        if type(value) is not str:raise ValueError('Exact decimal O3 power coordinate required')
        q=Fraction(value)
        if not 0<=q<=2:raise ValueError('Original O3 power0<=x<=2 required')
        inlet=self.transition('1');out={}
        for label,ctx,mu in (('point',self.ctx,self.mu),('directed',self.iv,self.mui)):
            x=ctx.mpf(value);old=inlet[label];s=dict(old);s['J']=ctx.mpf('.5')+x
            s['U1']=ctx.expm1(-mu*s['J'])/mu;s['U2']=ctx.expm1(-2*mu*s['J'])/mu
            for name,alpha,factor in (('C1',1,1),('C0',0,2),('Cm',-1,2)):
                moments=[self._moment(x,alpha,n,ctx) for n in range(1,5)]
                local=-factor*moments[0]+factor**2*mu*moments[1]/2-factor**3*mu**2*moments[2]/6
                if ctx is self.iv:
                    remainder=factor**4*mu**3*moments[3]/24
                    local+=ctx.mpf([0,_bounds(self.ctx,remainder)[1]])
                s[name]=old[name]+ctx.exp(alpha)*local
            out[label]=s
        out.update(method='original power_flow rewritten as normalized weighted expm1 integrals; explicit cubic-mu Taylor remainder',
            power_exp_series_terms=100,mu_remainder_order=3,transition_inlet_retained=True)
        return out

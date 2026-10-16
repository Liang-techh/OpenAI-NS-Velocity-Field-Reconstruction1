"""Directed Taylor integration of the original O2 slope scalar profiles.

Taylor coefficients include1/k!. A degree4 endpoint polynomial uses the
ODE-derived fifth coefficient over a whole-step value tube as remainder.
"""
from fractions import Fraction

from .original_core_numeric import _bounds


class OriginalSlopeScalar:
    names=('F','H','G','C')

    def __init__(self,ctx,iv,*,cells=512):
        if type(cells) is not int or cells<16:raise ValueError('At least16 scalar cells required')
        from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
        self.ctx,self.iv,self.cells,self.sigma_jets=ctx,iv,cells,sigma_jets
        p=[ctx.mpf(1),ctx.mpf(5)/8,ctx.mpf(5)/6,ctx.mpf(5)]
        t=[iv.mpf(1),iv.mpf(5)/8,iv.mpf(5)/6,iv.mpf(5)]
        self.grid=[(p,t)];self.cache={};self.steps=0

    def _rational(self,q,ctx):return ctx.mpf(q.numerator)/q.denominator

    def _jets(self,ctx,state,a):
        F,H,G,C=([v] for v in state)
        for n in range(5):
            F.append(sum((a[k]*F[n-k] for k in range(n+1)),ctx.mpf(0))/(n+1))
            H.append((F[n]-H[n]*ctx.mpf('1.5'))/(n+1))
            square=sum((F[k]*F[n-k] for k in range(n+1)),ctx.mpf(0))
            G.append((square-G[n])/(n+1));C.append(square/(n+1))
        return F,H,G,C

    def _intersection(self,a,b):
        c,iv=self.ctx,self.iv;lo,hi=_bounds(c,a);ll,hh=_bounds(c,b)
        if max(lo,ll)>min(hi,hh):raise ArithmeticError('Scalar Taylor enclosure misses its source value tube')
        return iv.mpf([max(lo,ll),min(hi,hh)])

    def _step(self,left,width,state):
        c,iv=self.ctx,self.iv;p,t=state;h=self._rational(width,iv);hp=self._rational(width,c)
        x=self._rational(left,iv);xl,xh=_bounds(c,x);wl,wh=_bounds(c,h)
        timebox=iv.mpf([xl,_bounds(c,x+h)[1]])
        # Differential comparison proves a tube for every intermediate time:
        # -.5<=a<=.1; all four states and the forcing F/F² are positive.
        lowF,highF=_bounds(c,t[0]);Fbox=iv.mpf([_bounds(c,t[0]*iv.exp(-h/2))[0],_bounds(c,t[0]*iv.exp(h/10))[1]])
        fhi=_bounds(c,Fbox)[1]
        tube=[Fbox]
        for initial,rate,source in ((t[1],iv.mpf('1.5'),iv.mpf(fhi)),(t[2],iv.mpf(1),iv.mpf(fhi)**2)):
            lo=_bounds(c,initial*iv.exp(-rate*h))[0]
            hi=_bounds(c,initial+source*h)[1]
            tube.append(iv.mpf([max(c.mpf(0),lo),hi]))
        tube.append(iv.mpf([_bounds(c,t[3])[0],_bounds(c,t[3]+iv.mpf(fhi)**2*h)[1]]))
        def coeff(box):
            sigma=self.sigma_jets(iv,box)
            return [iv.mpf('0.1')-sigma[0]*iv.mpf('0.6')]+[-sigma[k]*iv.mpf('0.6') for k in range(1,5)]
        at=coeff(x);whole=coeff(timebox)
        point_a=[sum(_bounds(c,v),c.mpf(0))/2 for v in at]
        pj=self._jets(c,p,point_a);tj=self._jets(iv,t,at);remainder=self._jets(iv,tube,whole)
        pp=[sum((q[k]*hp**k for k in range(5)),c.mpf(0)) for q in pj]
        tt=[sum((q[k]*h**k for k in range(5)),iv.mpf(0))+r[5]*h**5 for q,r in zip(tj,remainder)]
        tt=[self._intersection(q,b) for q,b in zip(tt,tube)]
        self.steps+=1
        return pp,tt

    def evaluate(self,value):
        if type(value) is not str:raise ValueError('Exact decimal slope coordinate required')
        q=Fraction(value)
        if not 0<=q<=1:raise ValueError('Original O2 slope0<=x<=1 required')
        if q in self.cache:return self.cache[q]
        index=(q*self.cells).numerator//(q*self.cells).denominator
        while len(self.grid)<=index:
            left=Fraction(len(self.grid)-1,self.cells)
            self.grid.append(self._step(left,Fraction(1,self.cells),self.grid[-1]))
        point,truth=self.grid[index];left=Fraction(index,self.cells)
        if q!=left:point,truth=self._step(left,q-left,(point,truth))
        point=list(point);truth=list(truth)
        if q==1:
            # sigma(1-x)=1-sigma(x) implies integral_0^1 sigma=1/2.
            # Hence F(1)=exp(1/10-3/10)=exp(-1/5), as a source identity.
            truth[0]=self._intersection(truth[0],self.iv.exp(-self.iv.mpf(1)/5))
            point[0]=self.ctx.exp(-self.ctx.mpf(1)/5)
        result=dict(point=dict(zip(self.names,point)),directed=dict(zip(self.names,truth)),
            method='degree4 Taylor + ODE fifth-coefficient interval remainder',uniform_cells=self.cells,
            endpoint_F_symmetry_identity_used=q==1,scalar_coefficients_are_Taylor_not_ordinary_derivatives=True)
        self.cache[q]=result;return result

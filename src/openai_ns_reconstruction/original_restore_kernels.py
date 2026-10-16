"""Positive scalar pressure kernels for the original axial restoration.

Point RK4/Simpson values are numerical proxies. Directed monotone cells,
including an uncertain endpoint segment, define the complete enclosures.
"""
from .original_core_numeric import _bounds


def powers(code):
    return code % 3, (code//3) % 3, code//9


MONOMIALS = tuple(i+3*j+9*k for i in range(3) for j in range(3-i)
                  for k in range(3-i-j))


class OriginalRestoreKernels:
    def __init__(self, reference, *, cells=128):
        if type(cells) is not int or cells < 4:
            raise ValueError('At least four positive restoration cells required')
        self.reference = reference
        self.ctx, self.iv = reference.ctx, reference.iv
        self.cells, self.cache = cells, {}

    def alpha(self, ctx, t):
        if ctx is self.iv:
            value = 1-self.reference.sigma_jets(ctx,t)[0]
            lo,hi = _bounds(self.ctx,value)
            return ctx.mpf([max(0,lo),min(1,hi)])
        if t <= 0: return ctx.mpf(1)
        if t >= 1: return ctx.mpf(0)
        odds = 1/(1-t)**2-1/t**2
        # All point tail truncations remain inside the independently
        # computed directed enclosure; they never redefine sigma.
        if odds > 10000: return ctx.mpf(0)
        if odds < -10000: return ctx.mpf(1)
        return 1/(1+ctx.exp(odds))

    def _point(self,t):
        c = self.ctx; h=t/self.cells; gamma=c.mpf(1)
        values = {code:c.mpf(0) for code in MONOMIALS}
        def step(a,g,h):
            k1=self.alpha(c,a)-g
            k2=self.alpha(c,a+h/2)-(g+h*k1/2)
            k3=self.alpha(c,a+h/2)-(g+h*k2/2)
            k4=self.alpha(c,a+h)-(g+h*k3)
            return g+h*(k1+2*k2+2*k3+k4)/6
        def densities(a,g):
            al=self.alpha(c,a);e=c.exp(-a);weight=c.exp(2*a)
            return {code:weight*al**i*g**j*e**k for code in MONOMIALS
                    for i,j,k in [powers(code)]}
        for n in range(self.cells):
            a=h*n; mid=step(a,gamma,h/2); end=step(a+h/2,mid,h/2)
            left,center,right=densities(a,gamma),densities(a+h/2,mid),densities(a+h,end)
            for code in values:
                values[code]+=h*(left[code]+4*center[code]+right[code])/6
            gamma=end
        return values,gamma

    def evaluate(self, point, truth):
        c,iv=self.ctx,self.iv;lo,hi=_bounds(c,truth)
        if lo < 0 or hi > 1: raise ValueError('Original restoration coordinate in[0,1] required')
        key=truth._mpi_
        if key in self.cache:return self.cache[key]
        values={code:iv.mpf(0) for code in MONOMIALS};gamma=iv.mpf(1)
        alpha0=self.alpha(iv,iv.mpf(0))
        cells=[]
        if lo:
            for n in range(self.cells):
                a=iv.mpf(lo)*n/self.cells;b=iv.mpf(lo)*(n+1)/self.cells;h=b-a
                ab=self.alpha(iv,b);e=iv.exp(-h)
                alow=_bounds(c,ab)[0];ahigh=_bounds(c,alpha0)[1]
                alpha=iv.mpf([alow,ahigh])
                end=e*gamma+(1-e)*alpha
                gl,gh=_bounds(c,end);oldhi=_bounds(c,gamma)[1]
                end=iv.mpf([max(alow,gl),min(1,oldhi,gh)])
                boxgamma=iv.mpf([_bounds(c,end)[0],oldhi])
                boxw=iv.mpf([_bounds(c,iv.exp(-b))[0],_bounds(c,iv.exp(-a))[1]])
                weight=iv.exp(2*a)*iv.expm1(2*h)/2
                for code in values:
                    i,j,k=powers(code);values[code]+=weight*alpha**i*boxgamma**j*boxw**k
                cells.append(dict(gamma_left=gamma,gamma_right=end,alpha=alpha,weight=weight))
                gamma,alpha0=end,ab
        # Exact mpf endpoint lo has a separately bounded [lo,hi] tail.
        # Integrands are nonnegative and <=exp(2s); gamma decreases.
        extension=iv.mpf(hi)-iv.mpf(lo)
        if hi > lo:
            weight=iv.exp(2*iv.mpf(lo))*iv.expm1(2*extension)/2
            at_hi=self.alpha(iv,iv.mpf(hi))
            e=iv.exp(-extension)
            end=e*gamma+(1-e)*iv.mpf([_bounds(c,at_hi)[0],_bounds(c,alpha0)[1]])
            for code in values:
                values[code]+=iv.mpf([0,_bounds(c,weight)[1]])
            gamma=iv.mpf([max(0,_bounds(c,end)[0]),_bounds(c,gamma)[1]])
        point_values,point_gamma=self._point(point)
        # Pure w modes have exact primitives. Preserve their exponential
        # correlation rather than accumulating avoidable cell roundoff.
        for code,rate in ((0,2),(9,1),(18,0)):
            values[code]=truth if rate==0 else iv.expm1(rate*truth)/rate
            point_values[code]=point if rate==0 else c.expm1(rate*point)/rate
        result=dict(point=point_values,truth=values,gamma_point=point_gamma,gamma_truth=gamma,
            cells=cells,uncertain_endpoint_segment=hi>lo,cell_count=self.cells if lo else 0,
            contract=dict(original_alpha_source=True,gamma_ODE_positive_recurrence=True,
                whole_cell_monotone_gamma=True,all_degree_two_scalar_monomials=True,
                positive_exp2_measure=True,complete_uncertain_endpoint=True,
                exact_pure_exponential_primitives=True,
                point_RK4_Simpson_is_only_proxy=True,interval_midpoint_not_source=True))
        self.cache[key]=result;return result

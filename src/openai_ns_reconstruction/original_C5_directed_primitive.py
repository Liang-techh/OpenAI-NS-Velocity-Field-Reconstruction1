"""Directed angle-coordinate evaluation of the prescribed C5 phase primitive.

Input coefficients are outward enclosures of the *same* root functions and
their ordinary mixed Z/y derivatives. This computational kernel does not
invent those enclosures, choose the huge-N phase, or evaluate the incoming
mean history. Actual-source installation is a separate, fail-closed gate.
"""
from __future__ import annotations

import ast
from fractions import Fraction
import math
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_finite_N_signed_source import load, sha, digest, require

PARENT='lei_ren_part1_paper_compliant_current_original_C5_seam_joined_history.json'
PARENT_SHA256='35fdc828c112e5fb09bf28a3f709b1d8a9817b54cd3a99fc11f76cfd92251d16'
SOURCE='lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
INDICES=tuple(sorted(((k,j) for k in range(4) for j in range(2)),key=lambda p:(sum(p),p)))


def interval(ctx,value):
    """Reject floating point/nearest-rounded inputs at the public boundary."""
    if hasattr(value,'_mpi_'):
        return ctx.mpf(value)
    if type(value) is int:
        return ctx.mpf(value)
    if isinstance(value,Fraction):
        return ctx.mpf(value.numerator)/value.denominator
    if isinstance(value,(list,tuple)) and len(value)==2:
        lo,hi=(interval(ctx,x) for x in value)
        require(lo.a<=hi.b,'Reversed input interval')
        return ctx.mpf([lo.a,hi.b])
    raise ValueError('Exact rational or outward interval input required')


class MixedJet:
    """Taylor coefficients in Z^0..3 y^0..1; public inputs are ordinary rows."""
    def __init__(self,ctx,values=None):
        self.ctx=ctx
        self.c={p:ctx.mpf(0) for p in INDICES}
        if values:
            for p,v in values.items():
                require(p in self.c,'Unsupported mixed derivative order')
                self.c[p]=interval(ctx,v)

    @classmethod
    def ordinary(cls,ctx,values):
        require(set(values)==set(INDICES),'All eight mixed ordinary root rows required')
        return cls(ctx,{p:interval(ctx,v)/math.factorial(p[0]) for p,v in values.items()})

    def other(self,value):
        if isinstance(value,MixedJet):
            require(value.ctx is self.ctx,'Mixed interval contexts')
            return value
        return MixedJet(self.ctx,{(0,0):value})

    def __add__(self,value):
        q=self.other(value);return MixedJet(self.ctx,{p:self.c[p]+q.c[p] for p in INDICES})

    __radd__=__add__

    def __neg__(self):
        return MixedJet(self.ctx,{p:-v for p,v in self.c.items()})

    def __sub__(self,value):return self+-self.other(value)
    def __rsub__(self,value):return self.other(value)+-self

    def __mul__(self,value):
        q=self.other(value)
        return MixedJet(self.ctx,{(k,j):sum((self.c[a,b]*q.c[k-a,j-b]
            for a in range(k+1) for b in range(j+1)),self.ctx.mpf(0)) for k,j in INDICES})

    __rmul__=__mul__

    def square(self):
        out=self*self
        out.c[0,0]=self.c[0,0]**2
        return out

    def inverse(self):
        require(self.c[0,0].a>0 or self.c[0,0].b<0,'Jet denominator contains zero')
        out=MixedJet(self.ctx,{(0,0):1/self.c[0,0]})
        for k,j in INDICES[1:]:
            known=sum((self.c[a,b]*out.c[k-a,j-b] for a in range(k+1)
                for b in range(j+1) if (a,b)!=(0,0)),self.ctx.mpf(0))
            out.c[k,j]=-known/self.c[0,0]
        return out

    def __truediv__(self,value):return self*self.other(value).inverse()
    def __rtruediv__(self,value):return self.other(value)*self.inverse()

    def sqrt(self):
        require(self.c[0,0].a>0,'Positive square root jet required')
        out=MixedJet(self.ctx,{(0,0):self.ctx.sqrt(self.c[0,0])})
        for k,j in INDICES[1:]:
            known=sum((out.c[a,b]*out.c[k-a,j-b] for a in range(k+1)
                for b in range(j+1) if (a,b) not in ((0,0),(k,j))),self.ctx.mpf(0))
            out.c[k,j]=(self.c[k,j]-known)/(2*out.c[0,0])
        return out

    def cos(self):
        # The ideal of zero-constant jets has nilpotence order five here.
        h=MixedJet(self.ctx,{p:v for p,v in self.c.items() if p!=(0,0)})
        cycle=(self.ctx.cos(self.c[0,0]),-self.ctx.sin(self.c[0,0]),
               -self.ctx.cos(self.c[0,0]),self.ctx.sin(self.c[0,0]))
        out=MixedJet(self.ctx);power=self.other(1)
        for n in range(5):
            out=out+power*(cycle[n%4]/math.factorial(n));power=power*h
        return out

    def ordinary_rows(self):
        return {p:v*math.factorial(p[0]) for p,v in self.c.items()}


def integrate_unit(body,end,panels):
    """Interval Riemann enclosure after h=end*s, including moving-limit jets.

    Each full cell is supplied to the body. No quadrature midpoint or
    estimated truncation error is used. Refinement tightens an enclosure;
    it does not turn its midpoint into the function value.
    """
    require(type(panels) is int and panels>=2,'At least two interval panels required')
    ctx=end.ctx;total=end.other(0)
    for n in range(panels):
        cell=ctx.mpf([interval(ctx,Fraction(n,panels)).a,
                      interval(ctx,Fraction(n+1,panels)).b])
        total=total+body(end*cell)/panels
    return end*total


class DirectedPrimitive:
    """C and fixed-phase C_y, through ordinary Z3, on an active loop.

    Roots E,a,t0,q,p2,dstar are supplied with all Z0..3/y0..1 rows.
    dstar may be constant but its zero rows must still be explicit. The
    source relation b=-a*t0 is retained before interval evaluation.
    """
    def __init__(self,roots,*,dps=35,panels=12,inverse_steps=16):
        require(set(roots)=={'E','a','t0','q','p2','dstar'},'Actual root tuple required')
        require(type(dps) is int and dps>=25,'Directed precision >=25 digits required')
        self.ctx=MPIntervalContext();self.ctx.dps=dps
        self.panels=panels;self.inverse_steps=inverse_steps
        self.roots={n:MixedJet.ordinary(self.ctx,v) for n,v in roots.items()}
        E,a,t0,q,p2,ds=(self.roots[n] for n in ('E','a','t0','q','p2','dstar'))
        require(E.c[0,0].a>0 and a.c[0,0].a>0 and ds.c[0,0].a>0,
            'Positive actual E,a,dstar required')
        u=p2*q/ds;hinv=(1+u.square()).sqrt().inverse()
        self.r=u*hinv;self.alpha=2*q*hinv
        require(self.r.c[0,0].a>-1 and self.r.c[0,0].b<1,
            'A strictly interior Poisson parameter is required')
        self.K=((1+t0.square()+2*q.square())*2*self.ctx.pi).inverse()
        self.twopi=2*self.ctx.pi

    def direction(self,h):
        c=(h*self.twopi).cos();r=self.r
        # Tighten only the constant coefficient by the exact correlated
        # inequality D>= (1-|r|)^2. Derivative coefficients stay unchanged.
        denominator=1-2*r*c+r*r
        rabs=max(abs(r.c[0,0].a),abs(r.c[0,0].b))
        floor=(1-rabs)**2
        denominator.c[0,0]=self.ctx.mpf([max(denominator.c[0,0].a,floor.a),
            denominator.c[0,0].b])
        return self.roots['t0']+self.alpha*(c-r)/denominator

    def T1_and_Phi(self,h):
        T=integrate_unit(self.direction,h,self.panels)*self.twopi
        T2=integrate_unit(lambda x:self.direction(x).square(),h,self.panels)*self.twopi
        return T,self.K*(h*self.twopi+T2)

    def lam_h(self,h):
        t=self.direction(h)
        return self.K*self.twopi*(1+t.square())

    def inverse_interval(self,phase):
        """Outward normalized-angle bracket for every supplied phase value."""
        ctx=self.ctx;p=interval(ctx,phase)
        require(p.a>=0 and p.b<=1,'Issued phase must lie in [0,1]')
        if p.a==p.b and any(p.a==ctx.mpf(n)/2 for n in (0,1,2)):
            # Reflection gives Phi(1/2)=1/2 for this exact source loop.
            return p
        if p.a==0 and p.b==1:return p
        lower=Fraction(0);right=Fraction(1)
        for _ in range(self.inverse_steps):
            mid=(lower+right)/2
            _,value=self.T1_and_Phi(MixedJet(ctx,{(0,0):mid}))
            if value.c[0,0].b<p.a:lower=mid
            else:right=mid
        left=Fraction(0);upper=Fraction(1)
        for _ in range(self.inverse_steps):
            mid=(left+upper)/2
            _,value=self.T1_and_Phi(MixedJet(ctx,{(0,0):mid}))
            if value.c[0,0].a>p.b:upper=mid
            else:left=mid
        return ctx.mpf([interval(ctx,lower).a,interval(ctx,upper).b])

    def inverse_jet(self,phase):
        ctx=self.ctx;p=interval(ctx,phase);base=self.inverse_interval(p)
        out=MixedJet(ctx,{(0,0):base})
        if p.a==p.b and any(p.a==ctx.mpf(n)/2 for n in (0,1,2)):
            return out
        lam=self.lam_h(MixedJet(ctx,{(0,0):base})).c[0,0]
        require(lam.a>0,'Directed inverse derivative lower bound required')
        for index in INDICES[1:]:
            _,phi=self.T1_and_Phi(out)
            out.c[index]=-phi.c[index]/lam
        return out

    def evaluate(self,phase):
        ctx=self.ctx;p=interval(ctx,phase)
        require(p.a>=0 and p.b<=1,'Issued phase must lie in [0,1]')
        if p.a==p.b and (p.a==0 or p.a==1):
            return MixedJet(ctx),dict(exact_periodic_gauge_zero=True,
                phase_selected=False,root_values_installed=False)
        h=self.inverse_jet(p)
        def body(angle):
            T,phi=self.T1_and_Phi(angle)
            B=self.roots['E']*self.roots['a']/2*(self.roots['t0']*phi-T/self.twopi)
            return B*self.lam_h(angle)
        result=integrate_unit(body,h,self.panels)
        return result,dict(exact_periodic_gauge_zero=False,phase_selected=False,
            input_phase_is_outward_interval=True,normalized_angle_bracket=str(h.c[0,0]),
            ordinary_Z0_Z3=True,fixed_phase_y_derivative=True,interval_Riemann=True,
            midpoint_or_error_cap_used_as_value=False,root_values_installed=False)


def _roots(**values):
    return {name:{p:values[name].get(p,0) for p in INDICES} for name in ('E','a','t0','q','p2','dstar')}


def interval_packet(value):
    return dict(display=str(value),exact_binary_endpoints=[list(v) for v in value._mpi_])


def numeric_controls():
    """New non-actual test loops: exact nonlinear half-phase C and mixed rows."""
    roots=_roots(E={(0,0):2,(1,0):2,(0,1):1,(1,1):1},a={(0,0):3},
        t0={},q={(0,0):Fraction(1,4),(1,0):Fraction(1,20),(0,1):Fraction(1,30)},
        p2={},dstar={(0,0):1})
    evaluator=DirectedPrimitive(roots,panels=12,inverse_steps=10)
    ctx=evaluator.ctx;computed,receipt=evaluator.evaluate(Fraction(1,2))
    q=evaluator.roots['q'];E=evaluator.roots['E'];a=evaluator.roots['a']
    # For t0=r=0, psi=pi at phi=1/2, and the sine integral is elementary.
    exact=-E*a*q*evaluator.K/ctx.pi*(1+4*q*q/3)
    contains={str(p):bool(computed.c[p].a<=exact.c[p].a and computed.c[p].b>=exact.c[p].b)
        for p in INDICES}
    zero0,_=evaluator.evaluate(0);zero1,_=evaluator.evaluate(1)
    gauge=all(v.a==0 and v.b==0 for jet in (zero0,zero1) for v in jet.c.values())
    # Independently test the formal inverse and a nonzero Poisson parameter.
    nonflat=_roots(E={(0,0):2},a={(0,0):3},t0={(0,0):Fraction(1,10)},
        q={(0,0):Fraction(1,4)},p2={(0,0):1},dstar={(0,0):1})
    nonflat['q'][1,0]=Fraction(1,20);nonflat['q'][0,1]=Fraction(1,30)
    other=DirectedPrimitive(nonflat,panels=8,inverse_steps=8)
    inv=other.inverse_jet(Fraction(1,3));_,ph=other.T1_and_Phi(inv)
    inverse_contains=bool(ph.c[0,0].a<=other.ctx.mpf(1)/3<=ph.c[0,0].b)
    inverse_rows_zero=all(v.a<=0<=v.b for p,v in ph.c.items() if p!=(0,0))
    checks=dict(nonlinear_mixed_rows_contain_exact=all(contains.values()),
        periodic_endpoint_gauge_exact=gauge,nonzero_r_inverse_contains_target=inverse_contains,
        nonconstant_inverse_mixed_residuals_contain_zero=inverse_rows_zero,
        all_eight_mixed_rows_present=len(computed.c)==8)
    require(all(checks.values()),'Directed primitive numeric kernel control failed')
    return dict(passed=True,checks=checks,ordinary_derivative_contains_exact=contains,
        half_phase_C0=interval_packet(computed.c[0,0]),half_phase_exact_C0=interval_packet(exact.c[0,0]),
        computed_ordinary_rows={str(p):interval_packet(v) for p,v in computed.ordinary_rows().items()},
        nonconstant_nonzero_r_inverse_rows={str(p):interval_packet(v) for p,v in inv.ordinary_rows().items()},
        controls_are_actual_source_values=False,controls_replay_accepted_quadrature=False,
        arithmetic='mpmath directed intervals; full-cell interval Riemann sums',receipt=receipt)


def build(root,parent):
    root=Path(root);directory=root/'experiments/root_st073'
    require(sha(directory/PARENT)==PARENT_SHA256 and digest(parent)==digest(load(directory/PARENT)),
        'Accepted joined-history source parent required')
    source=directory/SOURCE
    require(sha(source)==parent['input_hashes'][SOURCE],'Actual prescribed loop source changed')
    tree=ast.parse(source.read_text(encoding='utf-8'))
    functions={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
    loop=ast.get_source_segment(source.read_text(encoding='utf-8'),functions['loop_functions'])
    generic=ast.get_source_segment(source.read_text(encoding='utf-8'),functions['generic_roots'])
    compact=lambda x:''.join(x.split())
    patterns=('t=m.add(m.axial(roots[\'t0\']),m.mul(m.axial(alpha),m.div(numerator,den)))',
        'K=a.div(one,a.scale(ratio,twopi)',
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
        "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)")
    guards={v:compact(v) in compact(loop) for v in patterns}
    guards['b_equals_minus_a_t0']=all(compact(v) in compact(generic) for v in
        ('aa,bb=divide(C),divide(B)',"t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))"))
    require(all(guards.values()),'Actual loop angle or root correlation source changed')
    ports=parent['actual_C5_signed_history_IBP']['source_ports']
    controls=numeric_controls()
    return dict(passed=True,scope='Directed numerical kernel with actual formula binding; actual root/phase/history inputs still required',
        parent_report=PARENT,parent_sha256=PARENT_SHA256,source_pin=dict(name=SOURCE,sha256=sha(source)),
        source_AST_guards=guards,source_executed=False,
        actual_B_rows={chart:ports[chart]['source_B_rows'] for chart in ('O2_buffer','O3_transition')},
        actual_slow_B_y_rows={chart:ports[chart]['source_slow_B_y_rows'] for chart in ('O2_buffer','O3_transition')},
        formula='C(phi)=integral_0^h B(2*pi*s)*2*pi*K*(1+t(2*pi*s)^2)ds; Phi(2*pi*h)=phi',
        root_correlation='b=-a*t0 retained before boxing',
        mixed_orders=dict(ordinary_Z=[0,1,2,3],fixed_phase_y=[0,1]),
        slow_y_semantics='All root y rows are physical slow-y derivatives with issued phase held fixed; neither C_phi nor total fast-y derivative',
        root_box_semantics='One same-source mixed box throughout each nested integral; no independent b interval',
        inversion='outward monotone angle bracket, then rectangular formal implicit coefficients',
        quadrature='full-cell outward interval Riemann sums after moving-endpoint substitution',
        endpoint_gauge='C(0)=C(1)=0; not an individually zero actual issued phase',
        half_phase='Phi(pi)=1/2 for all actual root parameters by reflection',
        numeric_kernel_controls=controls,
        actual_root_mixed_enclosures_installed=False,actual_N_phase_enclosed=False,
        actual_signed_C_or_slow_integral_values_installed=False,actual_signed_H_installed=False,
        previous_H0_retained=True,old_errors_or_caps_modified=False,
        accepted_owner_fixture_recurrence_or_quadrature_replayed=False,
        missing_inputs=['Actual active/flat support certificate before computing q; flat q and all mixed rows exact zero, split branch-crossing cells',
            'Same-source mixed Z0..3/y0..1 root interval enclosures',
            'Actual global issued phase enclosure, or an explicitly uniform phase enclosure',
            'Directed evaluation of remaining slow physical-y integrals and retained incoming H0'])

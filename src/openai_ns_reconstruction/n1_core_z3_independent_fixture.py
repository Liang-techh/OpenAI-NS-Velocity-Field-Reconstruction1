"""Independent exact-rational checks for the new high-Z core computation."""
from __future__ import annotations

from fractions import Fraction
import math


def _fraction(bits):
    sign, mantissa, exponent, _ = bits
    value = Fraction(int(mantissa), 1)
    value = value * (1 << exponent) if exponent >= 0 else value / (1 << (-exponent))
    return -value if sign else value


def _contains(interval, rational):
    expected = Fraction(int(rational.p), int(rational.q))
    return _fraction(interval._mpi_[0]) <= expected <= _fraction(interval._mpi_[1])


def run_fixture(build, green, c):
    """Compare against direct symbolic differentiation, not source code jets."""
    import sympy as s

    r,z=s.symbols('r z', real=True)
    z0=s.Rational(1,4); r0=s.Integer(1); de=s.Rational(1,5); lam=s.Integer(3); j=s.Integer(2)
    d=1-z*z; L=1-de*z*z; H=z*(1-de)/2+d*(4*z+j); ell=-lam*L*H/(H*H+(j/500)**2)
    ell_d=[s.diff(ell,z,k).subs(z,z0) for k in range(5)]
    fder=[s.Integer(1)]
    for n in range(5):
        fder.append(sum(s.binomial(n,k)*ell_d[k]*fder[n-k] for k in range(n+1)))
    anchor=s.Rational(7,3)
    F0=anchor*sum(fder[n]*(z-z0)**n/s.factorial(n) for n in range(6))
    phi=1+r*z**5+r**2*z**3
    mean=1+r*z**5/2+r**2*z**3/3
    U=mean+r*s.diff(mean,r)
    F=F0*phi
    Q=(2*z*U-(1-de)*z*mean-d*s.diff(mean,z))/L
    S=(Q+(1-de)*z*s.diff(Q,z)/2+r*s.diff(Q,r))/L+Q*(Q/2+r*s.diff(Q,r))
    S+=U*(-2*z*Q+d*s.diff(Q,z)-2*z*r*s.diff(Q,r))/L
    Np=-S/2+lam*(2*s.diff(Q,r)+r*s.diff(Q,r,2))
    def axial(f,a):return (a*z*f+d*s.diff(f,z)-2*z*r*s.diff(f,r))/L
    Nt=-axial(axial(F,-2-de),-3)
    Nz=-axial(axial(U,-1-de),-2)
    gs={'4':2*s.sqrt(r)*Np/lam,'5':2*Nt/lam,
        '6':2*Nz/lam-4*z*r*Np/(lam**2*L)}
    expected={key:s.diff(expr,z,3).subs({r:r0,z:z0}) for key,expr in gs.items()}

    def interval(q):
        q=s.Rational(q);return c.mpf(int(q.p))/int(q.q)
    def grid(expr):
        return {f'rho{i}_Z{k}':interval(s.diff(expr,r,i,z,k).subs({r:r0,z:z0}))
                for i in range(3) for k in range(6-i)}
    F0_at_anchor=interval(anchor)
    result=build(c,interval(r0),interval(z0),interval(de),interval(lam),interval(j),
                 grid(phi),grid(U),grid(mean),c.log(F0_at_anchor))
    checks={}
    for row in ('4','5','6'):
        total=c.mpf(0)
        for term in result['g'+row+'_Z3']:
            total+=term['coefficient']*(c.mpf(3)**term['Lambda_power'])*(F0_at_anchor**term['F0_power'])
        checks['direct_symbolic_g'+row+'_Z3']=_contains(total,expected[row])
    checks['g5_has_current_F0_factor']=[term['F0_power'] for term in result['g5_Z3']]==[1]
    checks['g4_g6_have_no_F0_factor']=all(term['F0_power']==0 for row in ('g4_Z3','g6_Z3') for term in result[row])

    def axis_grid(expr):
        return {f'rho{i}_Z{k}':interval(s.diff(expr,r,i,z,k).subs({r:0,z:z0}))
                for i in range(3) for k in range(6-i)}
    axis=build(c,c.mpf(0),interval(z0),interval(de),interval(lam),interval(j),
               axis_grid(phi),axis_grid(U),axis_grid(mean),c.log(F0_at_anchor))
    checks['g4_zero_at_axis']=all(term['coefficient'].a==0 and term['coefficient'].b==0 for term in axis['g4_Z3'])

    cells=[]
    for index in range(8):
        left=c.mpf(index)/8;right=c.mpf(index+1)/8;t=c.make_mpf((left._mpi_[0],right._mpi_[1]))
        one={'F0_power':0,'Lambda_power':0,'coefficient':c.mpf(1)}
        cells.append({'t':t,'rho':4*t*t,'g':{'4':[one],'5':[one],'6':[one]}})
    integrated=green(c,cells)
    truth={'F1':s.Rational(1,2),'Uz1':s.Integer(1),'K1':s.Rational(-1,2),
           'P1':s.Integer(2),'d_x_F1':s.Rational(1,2),'d_x_Uz1':s.Integer(1)}
    for name, terms in integrated['fields'].items():
        total=c.mpf(0)
        for term in terms:total+=term['coefficient_interval']
        checks['constant_Green_'+name]=_contains(total,truth[name])
    for cell in cells:
        variable={'F0_power':0,'Lambda_power':0,'coefficient':cell['rho']}
        cell['g']={'4':[variable],'5':[variable],'6':[variable]}
    integrated=green(c,cells)
    truth={'F1':s.Rational(2,3),'Uz1':s.Integer(1),'K1':s.Rational(-2,3),
           'P1':s.Rational(8,3),'d_x_F1':s.Rational(4,3),'d_x_Uz1':s.Integer(2)}
    for name, terms in integrated['fields'].items():
        total=c.mpf(0)
        for term in terms:total+=term['coefficient_interval']
        checks['rho_Green_'+name]=_contains(total,truth[name])
    return {'passed':all(checks.values()),'checks':checks,
            'comparison':'exact rational independent symbolic differentiation and exact constant Green integrals',
            'fixture_is_not_actual_field':True}

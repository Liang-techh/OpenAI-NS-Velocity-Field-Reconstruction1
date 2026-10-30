"""Independent exact synthetic check of word1 feedback and Green integration."""
from __future__ import annotations
from fractions import Fraction


def _q(bits):
    sign,mantissa,exponent,_=bits
    a=Fraction(int(mantissa),1)
    a=a*(1<<exponent) if exponent>=0 else a/(1<<(-exponent))
    return -a if sign else a


def _contains(interval, rational):
    value=Fraction(int(rational.p),int(rational.q))
    return _q(interval._mpi_[0])<=value<=_q(interval._mpi_[1])


def run_fixture(build,c):
    import sympy as s
    C0=('4,1','5,1','5,2','5,3','5,5','6,1','6,2','6,3','6,4','6,6')
    B1=('5,1','5,2','5,3','6,2','6,3','6,4')
    term=lambda value:{'F0_power':0,'Lambda_power':0,'coefficient':c.mpf(value)}
    low={str(k):{row:[term(v)] for row in ('4','5','6')} for k,v in enumerate((1,2,3))}
    old=[];matrix=[]
    for i in range(8):
        left=c.mpf(i)/8;right=c.mpf(i+1)/8;t=c.make_mpf((left._mpi_[0],right._mpi_[1]))
        old.append({'g':{row:[term(4)] for row in ('4','5','6')}})
        a={str(k):{name:[term(v)] for name in C0} for k,v in enumerate((1,2,3))}
        b={str(k):{name:[term(v)] for name in B1} for k,v in enumerate((1,2,3))}
        matrix.append({'t':t,'rho':4*t*t,'C0':a,'B1':b})
    result=build(c,low,old,matrix)
    t=s.symbols('t',positive=True);x=2*t
    masses={1:x*x/8,2:x*x/4,3:-x*x/8,4:x,5:x/4,6:x/2}
    def h(row):
        ccols=[int(k.split(',')[1]) for k in C0 if int(k.split(',')[0])==row]
        bcols=[int(k.split(',')[1]) for k in B1 if int(k.split(',')[0])==row]
        return 14*sum(masses[v] for v in ccols)+22*sum(masses[v] for v in bcols)
    K=t*(1-t*t);kernels={'F1':K,'Uz1':2*t*s.log(1/t),'K1':-K,'P1':s.Integer(1),'d_x_F1':t**3,'d_x_Uz1':t}
    rows={'F1':5,'Uz1':6,'K1':6,'P1':4,'d_x_F1':5,'d_x_Uz1':6}
    checks={}
    for field,terms in result['W1_Z2_core_exit'].items():
        enclosure=c.mpf(0)
        for item in terms:enclosure+=item['coefficient']
        exact=s.integrate(2*kernels[field]*h(rows[field]),(t,0,1))
        checks['exact_word1_'+field]=_contains(enclosure,exact)
    return {'passed':all(checks.values()),'checks':checks,
        'fixture_is_not_actual_field':True,'method':'independent symbolic integration of polynomial feedback'}

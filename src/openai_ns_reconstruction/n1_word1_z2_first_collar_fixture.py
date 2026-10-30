"""Independent interval checks of the six exact homogeneous propagators."""
from __future__ import annotations
from fractions import Fraction


def _rational_endpoint(bits):
    sign,mantissa,exponent,_=bits
    value=Fraction(int(mantissa),1)
    value=value*(1<<exponent) if exponent>=0 else value/(1<<(-exponent))
    return -value if sign else value


def _contains_rational(interval,value):
    return _rational_endpoint(interval._mpi_[0])<=value<=_rational_endpoint(interval._mpi_[1])


def run_fixture(enclose,c):
    values={'F1':1,'Uz1':-2,'K1':3,'P1':-4,'d_x_F1':5,'d_x_Uz1':-6}
    core={k:[{'F0_power':1,'Lambda_power':-1,'coefficient':c.mpf(v)}] for k,v in values.items()}
    sc=c.mpf(1)/2048;hb=c.mpf(1)/4;one=c.mpf(1);s=sc/4;y=hb*s
    f0=c.mpf(7)/3;lam=c.mpf(3)
    result=enclose(c,core,sc=sc,q=one,feedback_Z2_upper=one)
    factor=c.exp(-y)
    direct={'F1':1+(1-factor)*5,'Uz1':-2+y*(-6),
        'K1':factor*3-(1-factor)*(-6),'P1':c.mpf(-4),
        'd_x_F1':c.exp(-3*y/2)*5,'d_x_Uz1':c.exp(-y/2)*(-6)}
    checks={}
    for field,terms in result['core_transport_signed_sectors'].items():
        total=c.mpf(0)
        for term in terms:
            total+=term['coefficient_interval']*(hb**term['hb_power'])*(f0**term['F0_power'])*(lam**term['Lambda_power'])
        expected=direct[field]*f0/lam
        if field=='P1':
            checks['propagator_'+field]=_contains_rational(total,Fraction(values[field]*7,9))
        elif field=='Uz1':
            checks['propagator_'+field]=_contains_rational(total,(Fraction(-2)-Fraction(6,32768))*Fraction(7,9))
        else:
            checks['propagator_'+field]=bool(total.a<=expected.a and total.b>=expected.b)
    zero=enclose(c,core,sc=sc,q=c.mpf(0),feedback_Z2_upper=one)
    for field,terms in zero['core_transport_signed_sectors'].items():
        total=c.mpf(0)
        for term in terms:
            total+=term['coefficient_interval']*(hb**term['hb_power'])*(f0**term['F0_power'])*(lam**term['Lambda_power'])
        checks['inlet_identity_'+field]=_contains_rational(total,Fraction(values[field]*7,9))
    dx=2*(c.exp(y/2)-1)
    for field,item in result['local_feedback_symmetric_radii'].items():
        radius=item['positive_coefficient_radius']*(hb**item['hb_power'])
        needed=dx**item['hb_power']
        checks['local_kernel_bound_'+field]=bool(radius.b>=needed.b)
    return {'passed':all(checks.values()),'checks':checks,'fixture_is_not_actual_field':True,
        'method':'direct exponential homogeneous solution and physical delta-x bounds'}

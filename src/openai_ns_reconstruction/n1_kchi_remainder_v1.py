"""Positive majorants for the Kchi kernel tail after j=1, ordinary Z0..2."""
from __future__ import annotations

def build(c, *, a, chi_z1, chi_z2):
    """Return pointwise kernel truncation bounds for 0<=s<=x<=a<=2.

    Caller certifies real |chi|<=1 and bounds on ordinary chi derivatives
    over the same Z domain. These are absolute errors, never signed fields.
    They do not include forcing, radial integration, or feedback-word tails.
    """
    if a.a <= 0 or a.b > 2: raise ValueError('0 < a <= 2 required')
    m1=abs(chi_z1).b; m2=abs(chi_z2).b
    r=a*a/2; term=r*r/12
    q0=r/12; q1=r/8; q2=3*r/16; q11=r/4; q21=3*r/8
    if any(q.b>=1 for q in (q0,q1,q2,q11,q21)): raise ValueError('tail ratio not below one')
    p0=a*term/(2*(1-q0)); p1=a*m1*2*term/(2*(1-q1)); p2=a*(m2*2*term/(1-q1)+m1*m1*2*term/(1-q11))/2
    d0=2*term/(1-q1); d1=m1*4*term/(1-q2); d2=m2*4*term/(1-q2)+m1*m1*4*term/(1-q21)
    def nonnegative(v): return c.make_mpf((c.mpf(0)._mpi_[0],v._mpi_[1]))
    return {'F1':tuple(nonnegative(x) for x in (p0,p1,p2)), 'd_x_F1':tuple(nonnegative(x) for x in (d0,d1,d2)),
      'scope':'pointwise potential-kernel remainder j>=2 only', 'ordinary_Z_orders':[0,1,2],
      'source_or_feedback_included':False,'input_abs_chi_bound':1,'first_omitted_index':2,'bounds_are_not_field_values':True}

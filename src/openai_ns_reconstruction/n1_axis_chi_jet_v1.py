"""Original axis chi and ordinary Z derivatives on a caller-owned interval."""
from __future__ import annotations

def build(c, *, z, delta, j):
    """Evaluate chi=H²/(H²+(j/500)²) and its ordinary first two jets.

    z is real; delta and j are Z-independent source intervals. Input source
    ownership and positivity of j belong to the caller. No prior producer is
    invoked and no point is selected from any interval.
    """
    if j.a <= 0:
        raise ValueError('j must have strictly positive lower endpoint')
    one=c.mpf(1); s2=(j/c.mpf(500))**2
    H=z*((one-delta)/2)+(one-z**2)*(4*z+j)
    Hp=(one-delta)/2+4-12*z**2-2*j*z
    Hpp=-24*z-2*j
    h2=H**2; D=h2+s2
    if D.a<=0:
        raise ValueError('positive chi denominator required')
    # chi is nonnegative and <= 1 independently of interval dependency.
    raw=h2/D
    lo=max(raw.a,c.mpf(0)); hi=min(raw.b,c.mpf(1))
    chi0=c.make_mpf((lo._mpi_[0],hi._mpi_[1]))
    chi1=2*s2*H*Hp/D**2
    chi2=2*s2*((Hp**2+H*Hpp)/D**2-4*h2*Hp**2/D**3)
    return {'chi':(chi0,chi1,chi2),'H':H,'H_Z':Hp,'H_ZZ':Hpp,'sigma_squared':s2,
      'ordinary_Z_derivatives':True,'Z_independent_delta_and_j_required':True,
      'source_values_selected':False,'point_value_computed':False}

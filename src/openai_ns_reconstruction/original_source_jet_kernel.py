"""New ordinary-Z source-jet arithmetic through order two.

These are spatial source derivatives.  They do not implement the separate
order-dependent temporal Navier-Stokes coefficient recursion.
"""
from mpmath.libmp import fzero, mpf_cmp


def _rows(values):
    if len(values) not in (2,3) or any(value is None for value in values):
        raise ValueError("Actual ordinary source rows through Z1 or Z2 are required; no zero padding")
    return values


def ordinary_product(left,right):
    _rows(left); _rows(right)
    if len(left)!=len(right):
        raise ValueError("Source jet orders must match explicitly")
    a,b = left[:2]; x,y = right[:2]
    result = [a*x,b*x+a*y]
    if len(left)==3:
        result.append(left[2]*x+2*b*y+a*right[2])
    return result


def ordinary_quotient(numerator,denominator):
    """Directed quotient with a strictly positive actual source denominator."""
    _rows(numerator); _rows(denominator)
    if len(numerator)!=len(denominator):
        raise ValueError("Source jet orders must match explicitly")
    a,b = numerator[:2]; x,y = denominator[:2]
    lower = x._mpi_[0] if hasattr(x,"_mpi_") else x._mpf_
    if lower[3] < 0 or mpf_cmp(lower,fzero) <= 0:
        raise ValueError("A positive source denominator enclosure is required")
    q0 = a/x
    q1 = (b-q0*y)/x
    result = [q0,q1]
    if len(numerator)==3:
        result.append((numerator[2]-2*q1*y-q0*denominator[2])/x)
    return result


def ordinary_exp(ctx,values):
    _rows(values); a,b = values[:2]
    value = ctx.exp(a)
    result = [value,value*b]
    if len(values)==3:
        result.append(value*(values[2]+b*b))
    return result


def transported_lower_dy(actual_phi,comparison_phi,comparison_psi_radial,rho):
    """Compute (actual Phi/comparison Phi)*comparison Psi_r*rho.

    Copying a result between numerical contexts requires exact endpoint
    copying only.  A macro Jacobian is not part of this source formula.
    """
    ratio = ordinary_quotient(actual_phi,comparison_phi)
    return ordinary_product(ordinary_product(ratio,comparison_psi_radial),rho)

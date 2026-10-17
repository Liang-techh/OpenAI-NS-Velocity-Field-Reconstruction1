"""Directed scalar algebra for two source-bound quiet mean-square increments.

All arguments are already in a shared scalar unit. Physical radii and
amplitudes remain outside this kernel as pinned logarithmic scale atoms.
"""
from __future__ import annotations

from mpmath.libmp import mpf_cmp


def square_interval(iv,value):
    a,b=value._mpi_;zero=iv.mpf(0)._mpi_[0]
    if mpf_cmp(a,zero)>=0:return value**2
    if mpf_cmp(b,zero)<=0:return (-value)**2
    top=max((-value).b,value.b)
    return iv.mpf([0,top.b])**2


def quiet_scalar_difference(iv,forward,canonical,length,*,source_bound,gap=None,gap_source_bound=False):
    """Enclose -length/4 * (forward²-canonical²), with exact ordinary Z1=0.

    A supplied gap must enclose the SAME signed source difference. It is
    never inferred from a midpoint or overlapping independent boxes.
    """
    if source_bound is not True:raise ValueError('Common-unit scalar source identity required')
    if length.a<=0:raise ValueError('Positive physical logarithmic length required')
    actual_gap=forward-canonical
    if gap is None:gap=actual_gap
    else:
        if gap_source_bound is not True:raise ValueError('Same-source signed-gap proof required')
        if gap.a<actual_gap.a or gap.b>actual_gap.b:
            raise ValueError('Signed gap enclosure must refine the source difference enclosure')
    factored=-(length/4)*gap*(forward+canonical)
    direct=-(length/4)*(square_interval(iv,forward)-square_interval(iv,canonical))
    lo=max(factored.a,direct.a);hi=min(factored.b,direct.b)
    if lo>hi:raise ArithmeticError('Same scalar expression enclosures are disjoint')
    return [iv.mpf([lo,hi]),iv.mpf(0)],{
        'identity':'-(Lq/4)*(Kf-Kc)*(Kf+Kc)',
        'ordinary_Z_orders':[0,1],'ordinary_Z1_exact_zero':True,
        'factored_and_direct_same_expression_intersection':True,
        'enclosure_overlap_used_as_equality_or_source_value':False,
        'physical_scale_materialized':False,'source_gap_selected_as_value':False}

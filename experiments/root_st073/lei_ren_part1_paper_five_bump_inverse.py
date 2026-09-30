"""Incremental Picard inverse of the paper's five-bump moment map.

Solve A h + Q(h,h) = -d, retaining every increment separately. First-Z
tangents propagate in the supplied coefficient ring. This is the analytic
five-moment inverse, not the temporal n-dependent background recursion.
"""
import operator
import mpmath as mp


def _sum_vectors(vectors):
    return tuple(sum(vector[i] for vector in vectors) for i in range(5))


def iterate_inverse(moment_map, defects, amplitude, *, steps=8, majorant=None):
    """Return finite Picard increments and cancellation-resistant residual.

steps counts nonlinear updates after the linear inverse. A supplied majorant
is conditional on its external uniform C1 defect bound and unenclosed weights.
No source bound is inferred from a point evaluation.
"""
    steps=operator.index(steps)
    if steps<0:raise ValueError('steps must be nonnegative')
    d=tuple(defects)
    if len(d)!=5:raise ValueError('five defect directions required')
    with mp.workdps(moment_map.precision):
        increments=[tuple(-v for v in moment_map.linear_inverse(d))]
        receipts=[]
        # At H0=-A^-1 d the residual is Q(H0,H0).
        residual=moment_map.bilinear(increments[0],increments[0],amplitude)
        for update in range(1,steps+1):
            step=tuple(-v for v in moment_map.linear_inverse(residual))
            previous=increments[:]
            increments.append(step)
            # R(Hk) = Q(Hk,Hk)-Q(Hk-1,Hk-1). Compute the
            # polarization directly, never subtract rounded large totals.
            cross=[tuple(2*v for v in moment_map.bilinear(old,step,amplitude))
                   for old in previous]
            cross.append(moment_map.bilinear(step,step,amplitude))
            residual=_sum_vectors(cross)
            receipts.append(dict(update=update,increment=step,residual=residual))
        conditional=None
        if majorant is not None:
            raw=majorant['raw']
            if not raw['contraction_condition_passed']:raise ValueError('Supplied contraction condition fails')
            # Starting at H0 in the contraction ball: error H0 <= q*r,
            # and each subsequent iteration multiplies by at most q.
            conditional=raw['contraction_radius']*raw['lipschitz_bound']**(steps+1)
        return dict(h=_sum_vectors(increments),increments=increments,
                    terminal_residual=residual,iterations=receipts,
                    conditional_C1_error_bound=conditional,
                    metadata=dict(nonlinear_updates=steps,
                        first_Z_propagated='if supplied as AxialDual',
                        increment_hierarchy_retained=True,
                        materialized_h_can_round_tiny_increments=True,
                        exact_infinite_inverse=False,functional_closure=False,
                        temporal_recursion=False,quadrature_enclosed=False,
                        uniform_source_bound_verified=False,
                        conditional_bound_scope='external majorant assumptions only'))

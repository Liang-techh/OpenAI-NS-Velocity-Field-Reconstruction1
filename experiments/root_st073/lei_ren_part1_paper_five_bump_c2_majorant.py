"""Conditional five-bump inverse bounds in a submultiplicative C2 norm.

On [-a,a], use ||f|| = sup|f| + sup|f'| + sup|f''|/2.
Leibniz proves ||fg|| <= ||f|| ||g||. No actual source norm is inferred.
The finite map quadrature and matrix constants are not enclosed here.
"""
import mpmath as mp
from lei_ren_part1_paper_five_bump_majorant import compute_majorant


def amplitude_inverse_square_c2_bound(Pstar, axial_radius):
    """Analytic upper bound for Am^-2=e^1.2(1+Z^2)^2/Pstar^2.

    This is the sum of three separately bounded derivative suprema,
    not the supremum of their pointwise sum. Scalars use caller precision.
    """
    p=mp.mpf(Pstar);a=mp.mpf(axial_radius)
    if not mp.isfinite(p) or p<1 or not mp.isfinite(a) or not 0<=a<1:
        raise ValueError('Require finite Pstar>=1 and 0<=axial_radius<1')
    scale=mp.exp(mp.mpf('1.2'))/(p*p)
    return scale*((1+a*a)**2+4*a*(1+a*a)+(4+12*a*a)/2)


def compute_c2_majorant(moment_map, uniform_defect_c2_bound, Pstar, *,
                        axial_radius, degree=3, inverse_steps=10):
    """Bound response/tail conditionally on a supplied interval defect norm.

    The defect vector norm is sum_i ||d_i||_C2. An isolated C2 jet at one Z
    cannot establish this input. Returned bounds inherit finite-map scope.
    """
    if not isinstance(inverse_steps,int) or isinstance(inverse_steps,bool) or inverse_steps<0:
        raise ValueError('inverse_steps must be a nonnegative integer')
    with mp.workdps(max(80,moment_map.precision+20)):
        a=mp.mpf(axial_radius)
        amplitude_bound=amplitude_inverse_square_c2_bound(Pstar,a)
        result=compute_majorant(moment_map,uniform_defect_c2_bound,Pstar,degree,
            amplitude_inverse_square_bound=amplitude_bound,
            norm_label='sum_i (sup|d_i| + sup|d_i_Z| + sup|d_i_ZZ|/2)')
        raw=result['raw'];q=raw['lipschitz_bound']
        # H0 lies in the contraction ball. The initial linear residual gives
        # ||H*-H0|| <= radius*q; k subsequent Picard updates give q^k.
        raw['incremental_inverse_tail_bound']=raw['contraction_radius']*q**(inverse_steps+1) if raw['contraction_condition_passed'] else mp.inf
        result['metadata'].update(axial_radius=mp.nstr(a,60),inverse_steps=inverse_steps,
            amplitude_C2_bound=mp.nstr(amplitude_bound,60),
            incremental_inverse_tail_bound=mp.nstr(raw['incremental_inverse_tail_bound'],60),
            interval_source_norm_verified=False,source_error_enclosed=False,
            amplitude_formula='Am=Pstar*exp(-0.6)/(1+Z^2)',
            norm_product_proof='Leibniz with second derivative factorial weight 1/2')
        return result

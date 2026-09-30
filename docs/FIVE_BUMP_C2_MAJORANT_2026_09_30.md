# Conditional C2 control for the five-bump inverse

The finite five-bump inverse now has a conditional interval C2 majorant.
Use the scalar norm N(f)=sup|f|+sup|f_Z|+sup|f_ZZ|/2 on [-a,a], a<1,
and the sum of these norms over five vector components. Leibniz's rule gives
N(fg)<=N(f)N(g): the cross term in the second derivative has coefficient
one after division by 2 and is included in the product of norms.

The common amplitude is Am=Pstar exp(-0.6)/(1+Z^2). Its inverse square is
an explicit polynomial, so a valid real interval bound is

    N(Am^-2) <= exp(1.2)/Pstar^2 *
       [(1+a^2)^2 + 4a(1+a^2) + (4+12a^2)/2].

This replaces the C1 pressure-amplitude bound in the quadratic tensor
estimate. The other quadratic terms have constant coefficients. The fixed
matrix l1 inverse bound CA and bilinear bound CQ therefore give the same
Banach-algebra contraction formulas, with e=sum_i N(d_i):

    e <= 1/(8 CA^2 CQ), radius=2 CA e, q=4 CA^2 CQ e <= 1/2.

Starting at H0=-A^-1 d, the distance to the fixed point is bounded by
radius*q. After k subsequent Picard updates a conservative tail is
radius*q^(k+1). This bounds all three derivatives in the declared norm,
conditionally on a valid uniform defect bound and valid operator constants.
The existing Catalan coefficient-series tail is also returned.

## Implementation and checks

`lei_ren_part1_paper_five_bump_c2_majorant.py` implements the amplitude formula
and wraps the shared majorant with an explicit norm contract. The shared
function retains the default C1 behavior and its historical receipt keys.

`lei_ren_part1_paper_five_bump_c2_majorant_fixture.py` uses the analytically
bounded family d_i=c_i(1+Z^2) on [-.8,.8]. The defect norm is calculated
from polynomial derivative suprema, independently of the sampled inverse.
At four axial points the 10-update value/first/second errors relative to a
40-update reference lie below the conditional tail 1.15111e-33. The reference
is a finite comparison, not an independently certified exact inverse.
The contraction boundary gives q=.5 and twice the threshold is rejected.
The original C1 zero/boundary/failed smoke checks also pass.

## Scope

No actual common-source defect interval norm has yet been enclosed. Finite
bump quadrature and matrix constants are not enclosed. Floating-point
evaluation of these formulas is not a directed-rounding certificate.
Thus this API provides the missing C2 conditional algebra and tail contract,
not a proof of actual uniform closure. It must not be fed a pointwise jet
norm and described as an interval certificate. Relative-flat row closure,
endpoint strips and source numerical remainders remain separate requirements.

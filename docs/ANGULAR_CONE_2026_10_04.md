# Whole original angular two-vector stress cone — 2026-10-04

## Accepted increment

The original angular repair s in [-4,0], Z in [-1,1] has the original two-vector stress cone against the SAME complete A/E/P and analytic absolute pressure used by the accepted entry and following tail. Original selected axial coefficients, two compact beta pulses, source velocities, radius, Ts/Lrel and physical viscosity units remain.

The continuous proof covers the whole angular interval, both supports and their crossings. It does not use phase samples or the raw whole-theta interval as positivity evidence.

Accepted regional source/physical joins compose the nonzero two-vector cone tail

~~~
Rtail*exp(-wait-Ts-6) <= R < Rtail*exp(3)
~~~

Gamma stress is exactly zero separately beyond the collar. The preceding-power region, its left physical/stress join, original 100-unit flatten, completed full-tensor/global cone, independent global flatness/energy, actual recursion and final corrected NS residual are still open.

Files:

- experiments/root_st073/lei_ren_part1_paper_compliant_angular_cone.py
- experiments/root_st073/lei_ren_part1_paper_compliant_angular_cone.json
- experiments/root_st073/lei_ren_part1_paper_compliant_angular_cone_check.py
- experiments/root_st073/lei_ren_part1_paper_compliant_angular_cone_check.json

Controller stage: angularcone.

## Stable native cumulative history

Write r=1-mu, k=1-a, b=(1-delta)/2, delta=2*a, q=-wait-Ts-2+s, G=KR*exp((a-mu)*s), K=G*F and F=1+h. The native cumulative angular moment is A=G*N, where

~~~
N = 1/r + (Xf-1/r)*exp(-r*(Lrel+s))
    +sum_j d_j(Z)*exp(-r*(s-center_j))*A_j(s)
~~~

Here Xf is the actual flatten exit, center_j is -3 or -1, and A_j(s) is the original normalized positive beta cumulative integral with its actual clamped endpoint. The selected d_j(Z) are the source physical_coefficient_Taylor functions.

The inertial theta numerator is G*W with

~~~
W = k*N-b*Z*N_Z-F
  = (mu-a)/r
    +exp(-r*(Lrel+s))*(k*(Xf-1/r)-b*Z*Xf_Z)
    +sum_j exp(-r*(s-center_j))*A_j(s)*(k*d_j-b*Z*d_j_Z)
    -h
~~~

The equilibrium surplus (mu-a)/r is retained before interval subtraction. The actual native source ASTs replay this identity. Entry's different sigmoid ODE comparison is not reused on the angular interval.

Positive normalized beta density and actual clamping give 0<=A_j(s)<=the SAME full A weight. On the original whole domain, exp(-r*(Lrel+s))<=exp(-r*(Lrel-4)), and exp(-r*(s-center_j))<=exp(r*(4+center_j)). The latter is deliberately conservative and remains valid before, inside and after each support.

Whole-Z source jets bound Xf, d_j and d_j_Z. Native ordinary-beta rows bound |h| and |F_s| directly; the implementation does not compute tiny h by subtracting F-1. Their sum gives a uniform positive W lower bound. W_lower/mu is approximately 0.999999999953212 for the unchanged selected source family.

## Actual shear and positive theta

The original source gives

~~~
kappa-2 = 2*mu-2*F_s/F
shear/G = 2*S*exp(-q)*(F_s-(1+mu)*F)
~~~

F has a positive whole-domain enclosure. The source bounds establish mu>|F_s|/Fmin and (1+mu)*Fmin>|F_s|. Therefore the actual kappa-2 is positive and source shear is strictly negative. S=1/Rtail is strictly positive by the original source definition; a cap whose lower endpoint is zero is only used for an upper bound, never to define the field or strict sign.

The largest inverse-radius factor occurs at s=-4, giving exp(wait+Ts+6). Let Sigma_max bound the absolute shear/G there. Since a<mu and s<=0, G>=KR; and 0<L=1-delta*Z^2<=1. Thus

~~~
Ctheta >= KR*(W_lower-Sigma_max) > 0
~~~

This is the positivity proof. The raw whole-theta hull crosses zero because it loses shared unit/history correlations, and is retained only as a conservative signed enclosure.

## Directional inequality and physical composition

The actual source B decreases with s. Bmax occurs at q=-wait-Ts-6. Its native logarithmic source parts regroup the common Ts/wait terms into -k*Ts+6*bh before enclosure. Each other part and the complete source log sum are AST-verified; no numerical B cap is selected as a field.

Using the actual positive variable-rate upper m_upper=2*mu+2*|F_s|/Fmin and the SAME whole absolute Cz upper, the checker verifies

~~~
log(m_upper)+2*log(Bmax)+2*log(Cz_upper/Ctheta_lower) < log(2)
~~~

Together with negative source shear and positive theta, this proves the original two-vector cone continuously. Common positive nu/lambda factors preserve it under the accepted physical map. The completed diagonal is not included in this regional two-vector assertion; full-tensor/global gates remain false.

Current angular K4/full moments/stress3/pressure4 and completed physical stress/diagonal/divergence/remainder right joins identify the SAME entry field. The already accepted entry/power/exit/waiting/collar cone then composes without another tail replay.

## Focused verification and next work

Producer/checker PASS: 24 positive source margins and 374 current input hashes. Eleven exact identities and 15 source-bridge facts are replayed; whole-domain bounds, actual variable shear and source Bmax are recomputed. No phase grid or shortened domain is used.

CertifiedAngularPhysical attaches regional cone metadata to the unchanged angular physical method and preserves its generally nonzero axial-viscosity remainder. An actual interior call at Z=.5, s=-3, log_tau=-2, theta=.3, nu=.7 passes with the regional scope retained. Global flatness and physical energy are not inferred from this regional result.

NEXT: SAME complete moments/stress/pressure through preceding original power and its left angular joins; original 100-unit flatten and its regional physical/cone companions; finite-width bridge feedback; independent global flat/volume/energy bounds; actual n-dependent recursion, oscillatory cancellation and final corrected Cartesian residual.

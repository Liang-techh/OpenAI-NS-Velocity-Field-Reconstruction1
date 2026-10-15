# Actual first-collar radial derivatives and core seam — 2026-10-15

Subsequent progress: [physical first gradient and vorticity](CURRENT_ORIGINAL_FIRST_COLLAR_GRADIENT_2026_10_15.md)
are now evaluated on the admitted first-half collar with original
coordinate/F0/delta factors and a core-boundary first-derivative seam.
Physical second derivatives and complete Cartesian C2 joins remain open.

Subsequent progress: [full original first-half collar Z2, q_Z and q_rho_Z](CURRENT_ORIGINAL_FIRST_COLLAR_AXIAL2_2026_10_15.md)
are now evaluated through a separate native API on 0<=s<=1/2. Physical
mixed-derivative assembly and complete Cartesian C2 joins remain open.

Subsequent progress: [full n1 core Z2 and q_Z numerical data](CURRENT_ORIGINAL_N1_CORE_AXIAL2_2026_10_15.md)
are now available through a separate core-only API. Transporting them
through the actual collar and complete physical mixed derivatives remain open.

## Evaluated derivatives

`OriginalFirstCollarDerivatives.evaluate(Z,s)` evaluates native radial
derivatives of the actual leading and full n1 field on `|Z|<1, 0<=s<=1/2`.
It uses the original coupled leading source and enclosed full n1 solution.
Second radial values are recovered from the local differential equation;
no nominal unknown radial polynomial or finite difference is differentiated.

```python
from openai_ns_reconstruction.original_first_collar_derivatives import OriginalFirstCollarDerivatives
field = OriginalFirstCollarDerivatives()
packet = field.evaluate('.3', '.5')
packet['leading_radial_derivatives']['Phi_rhorho']
packet['n1_radial_derivatives']['f_rhorho']
packet['n1_radial_derivatives']['p_rho']
```

All radial derivatives below are ordinary derivatives in rho, where
`rho=4exp(hb*s)`, `R=epsilon*rho`, and the computational radial variable is
`x=sqrt(rho)`. This x is not a physical Cartesian coordinate.

| Field | Radial orders exported | Available native Z Taylor orders |
| --- | --- | --- |
| Leading Phi, Psi, own Mpsi | 1, 2 | 0, 1 |
| Full leading Q | 1 | 0, 1 |
| Natural n1 f, u, k, own mean | 1 | 0, 1 |
| Natural n1 f, u, k | 2 | 0 |
| Natural n1 p | 1 | 0 |
| Recovered natural n1 q | 1 | 0 |

Here `Psi=(Uz_leading-(4Z+j))/epsilon`. The n1 frame is
`f=epsilon*F1/F0`, `u=Uz1/epsilon`, `k=K1/epsilon`,
`p=P1/(epsilon*Pstar^2)`, `fx=D_x f`, `ux=D_x u`.
Native axial coefficient k means `D_Z^k/k!`. Leading Q_rhorho, pressure
second-radial values, n1 q_Z and mixed derivatives requiring a numerical
W_ZZ are not supplied. This is not a complete physical Cartesian C2 API.

## Explicit inverse-width sector

The actual smooth control is `chi=1-sigma(s)+hb*sigma(s)`. Its exact
radial derivative is split as

`chi_rho = -sigma_prime(s)/(hb*rho) + sigma_prime(s)/rho`.

Both signs and the finite positive term are retained. With
`A=Phi_core_rho/Phi_core`, the leading derivatives obey

```
Phi_rho = Phi*chi*A
Phi_rhorho = Phi*(chi^2*A^2 + chi*A_rho + sigma_prime*A/rho)
          + hb^-1*(-Phi*sigma_prime*A/rho)

Psi_rho = (Phi/Phi_core)*chi*Psi_core_rho
Psi_rhorho = (Phi/Phi_core)*[
    chi*(chi-1)*A*Psi_core_rho + chi*Psi_core_rhorho
    + sigma_prime*Psi_core_rho/rho]
    + hb^-1*[-(Phi/Phi_core)*sigma_prime*Psi_core_rho/rho]
```

Each returned derivative is `sum(hb^power * coefficient)` with exact
directed coefficient endpoints and point-coefficient errors. The powers
are -1 and/or 0. Coefficients themselves retain finite-width interval
dependence through actual chi, actual/comparison values and rho. Thus
this is an enclosure of a factored identity, not a coefficientwise formal
Laurent polynomial or a source field with hb set to zero.

The point branch uses the inlet/rho=4 and zero-width coefficient proxy;
the error branch encloses the finite original-width source. The log(hb)
source is included separately. No hb inverse is exponentiated and no
width enclosure endpoint/midpoint is selected as the original field.

## Full n1 derivative recovery

Let `(rf,ru,rp)=linear_rhs(W,W_Z)` in the accepted natural frame, and let
`gf,gu,Np` be the actual known source. The operator applies the F0 gauge
once and retains original nonzero S, pressure and coupled feedback.

```
f_rho = fx/(2x)
u_rho = ux/(2x)
k_rho = -u_rho-k/rho
mean_u_rho = -k/rho

f_rhorho = (rf+gf-4fx/x)/(4rho)
u_rhorho = (ru+gu-2ux/x)/(4rho)
k_rhorho = -u_rhorho+u_rho/rho+2k/rho^2
p_rho = rp+Np

q_rho = [Z*ux/x + (1+delta)*Z*k/rho + (1-Z^2)*k_Z/rho]/L
```

The transport forcing packet named p is `2x*Np`, whereas the independent
core forcing polynomial's p row is Np. The derivative adapter divides
the transport packet by 2x before pressure recovery. This conversion is
part of the independent seam check.

The first-radial f/u/k/mean Z1 fields need the available fx_Z/ux_Z/k_Z.
Second radial values need W and W_Z, which are already enclosed. Their
next Z derivatives would need W_ZZ and are intentionally not claimed.
The mean and q identities use the same own k field; no fitted radial
velocity is introduced.

## Core seam

At s=0, the original flat pulse has sigma=sigma_prime=0. The inverse-width
coefficients vanish exactly, actual chi=1, actual/comparison values share
the original core endpoint, and the derivative equations reduce to the
core equations. The finite `+sigma_prime/rho` sector also vanishes.

The driver compares the actual phase-zero fx/ux known source and
`p/(2x)` source to the independent existing core forcing polynomial,
including its separate true-leading-source tail. It also compares leading
radial derivatives to the differentiated original inlet. That polynomial
is a core-seam reference only; it is never used as a finite-collar source.

These source/chain identities and the common nonzero state support the
native radial seam. They do not establish all physical mixed derivative
joins, the remaining collars or a global numerical C2 bound.

## Executed evidence

The accepted run exits 0 in **168.395 seconds**, producing ten original
radial derivative packets: Z=-0.3,0,0.3 at s=0,0.25,0.5, plus the offgrid
point Z=0.17,s=0.37. It passes:

- seven independent chain/kinematic/recovery/width-split identities;
- 33 independent core-known-source and leading radial seam comparisons;
- 24 exactly zero inverse-width coefficient slots at the core seam;
- 105 complete coefficient-error comparisons across precision 60 and 85;
- four unsupported-input rejections and 20 distinct original-width log proofs.

At Z=0.3,s=0.5, the normalized n1 pressure radial derivative is
`p_rho=hb^-1*a_minus1+a_zero`, with sampled coefficients

`a_minus1 = -10.72531587157666670210`, error below `8.246e-20`,

`a_zero = 10.79840156241311573778`, error below `8.302e-20`.

The sampled inverse-width coefficient relative errors are below 3.185e-43
for leading Phi_rhorho and n1 f_rhorho, and 7.688e-21 for n1 u_rhorho and
p_rho. These are coefficient errors, not a bound on an ordinary materialized
full derivative or physical relative error.

The 708070-byte report pins 1478 source inputs and 11 implementation files.
Fresh saved-report inspection exits 0. The existing read-only Luna/max
worker confirms the product/chain identities, pressure factor, n1 operator
and q recovery. The table above limits each exported derivative separately;
no leading Q_rhorho or n1 second-pressure derivative is implied.

Driver: `scripts/run_original_first_collar_derivatives.py`.
Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_collar_derivatives.json`.

## Remaining work

Tasks 14/15 remain IN PROGRESS; the full goal remains ACTIVE / INCOMPLETE.

1. Extend full n1 numerical axial data to the orders consumed by q_Z and
   physical mixed/second spatial derivatives, with their own source/kernel/
   feedback errors. Do not silently pad missing W_ZZ with zero.
2. Finish quantitative physical derivative joins and propagate evaluated
   data through the remaining admitted/extended collar charts.
3. Apply source-owned cutoff and own five-moment repair; complete global
   numerical matching and physical coordinate/chart assembly.
4. Complete exact heat exterior, finite energy, admissible stress and
   flat remainder, genuine temporal n>=2 recursion and its own repair,
   oscillatory corrections and independent corrected Cartesian NS checks.

Large inverse-width derivatives describe the fixed original thin spatial
layer. They are not a measured NS residual, a temporal blow-up test or a
proof of higher-order scale recursion.

# Original steep power: complete-future stress and pressure — 2026-10-04

Update: the regional physical decomposition and whole original Ts cone have subsequently passed. Read STEEP_POWER_PHYSICAL_2026_10_04.md and STEEP_POWER_CONE_2026_10_04.md for current scope and the next steep-entry dependency. The next-work list below records the dependencies identified at similarity-stage acceptance.

The original Ts-long steep-power region now recovers the same complete angular, energy and absolute-pressure moments as the accepted steep-exit/waiting/collar/Gamma continuation. The companion supplies signed similarity stress through mixed order 3 and absolute pressure through mixed order 4. Actual formula identities join its right endpoint to steep-exit phase 0. Original coefficients, velocities and pressure datum are retained.

This completes the steep-power similarity companion. Its physical stress/remainder transfer and whole-domain directional cone remain unfinished, as do steep entry and the preceding repair/power/flatten regions. True n-dependent coefficient recursion is not implemented by this increment.

## Original source and complete moments

Domain: original steep_power phase in [0,1], source Z in [-1,1]. The original map has

```
y = Ts*phase
ell = Ts*(1-phase)
q = log(R/Rtail) = -wait-1-ell
Ts = 4*(log(2)-log(delta))
a = delta/2, k = 1-a, p = 1+delta, bh = 1/2+a
KQ = (1-epsilon)*exp(k/2)
K = KQ*exp(k*ell)
B = Ev0*theta_base*exp(-bh*q)
C = (Ev0/Pstar)^2*theta_base^2
A = K*X_original
E = 2*K^2*energy_original
P = -pressure_absolute*exp(p*q)/C
```

The actual steep-exit phase-0 endpoint supplies AQ/EQ/PQ and their axial derivatives. Its future already includes the original waiting/collar and full Gamma tail; no new terminal data or finite-cutoff future replaces it. Exact source factors are distinguished from interval enclosures, including Ev2 and the tiny positive S=1/Rtail.

The source-normalized equations and their closed backward solutions are

```
A_y = K-k*A
E_y = delta*E-K^2
P_y = p*P-K^2/2
K_y = -k*K

A = exp(k*ell)*(AQ-KQ*ell)
E = K^2/2 + exp(-delta*ell)*(EQ-KQ^2/2)
P = K^2/6 + exp(-p*ell)*(PQ-KQ^2/6)
```

The exact rates delta+2*k=2 and p+2*k=3 give the denominators 2 and 6. Defect rows use expm1 and exact unit-baseline cancellation before enclosure. Ordinary y derivatives through order 4 feed the accepted stress recovery.

## Original angular correlation and endpoint binding

At the enormous original Ts, independently enclosing AQ/KQ and then subtracting Ts discards the small inlet X. The implementation instead retains the exact source identity

```
X = XS+Ts*phase
X_y = 1
X_yy = X_yyy = X_yyyy = 0
```

XS is the actual original source function, not a chosen interval endpoint. To prove that this is the same full angular future, the bridge replays actual XS/XQ/XT/waiting/exit expressions and the child backward-moment formulas. It consumes the accepted common zero-angular-history closure.

Both primitives use the same actual sigma callable. Source AST guards bind the two nextJ updates, J(1)=1/2, the backward nextf update, its zero-length terminal integral, and the original logistic return exp(odds)/(1+exp(odds)). Actual odds reflection gives sigma(v)+sigma(1-v)=1. FTC and the common endpoint then prove

```
f(v) = J(v)-v+1/2
exp(k*v)*(K(v)-1) = KQ*exp(k*J(v))-exp(k*v)
AQ/KQ = XQ = XS+Ts
```

These are source-functional equalities. Interval overlap is not used as a proof. The whole-Z inlet XS enclosure is approximately [1.491282615002014,1.506043547510141]. Preserving this correlation restores a positive whole-domain normalized theta-stress lower bound, approximately 0.8099880972790025. This bound alone is not a directional-cone certificate.

## Stress, pressure and joins

Positive source factors are Qtheta=sqrt(R/2)*B and Qz=sqrt(R/2)*B^2. Homogeneous angular/energy/pressure modes and opposing unit modes are combined analytically before interval evaluation. The original source shear gives kappa-2=2 on this power region; source S>0 proves the corresponding strict shear sign independently of any enclosing lower endpoint.

Absolute pressure retains its complete datum. With PQ=1/(2*p)+PdQ,

```
pressure_absolute = -C*exp(p*(wait+1))*(PQ+KQ^2*expm1(3*ell)/6)
```

The implementation retains original forward-history aliases and their pressure derivatives. The equivalent full-future representation supplies the companion stress and pressure rows.

Actual source-row ASTs at phase 1 and steep-exit phase 0, evaluated on arbitrary terminal axial functions, prove 5 K, 20 stress and 15 pressure identities. The stress join is mixed3 and the pressure join mixed4. Selected meridional terminal histories and the original zero-density FTC route retain the exact meridional primitive/velocity zeros.

## Reproduction and evidence

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage steeppowerstress
```

The focused producer and checker are the steep_power_stress_C3.py and steep_power_stress_C3_check.py companions. Current acceptance requires all_passed, matching family/source, every current input hash, the new endpoint/source proofs, true local gates and false physical/cone/global/recursion gates.

The checker covers 80 finite signed actual stress coefficients, 60 pressure coefficients and 96 exact meridional zero rows, including whole original phase/Z enclosures. One bounded independent moderate exponential-power fixture checks 6 finite original integrals, 30 complete-moment derivatives, 20 normalized quotient derivatives, 40 native unnormalized stress mixed3 coefficients and 30 absolute-pressure mixed4 coefficients at tolerance 1e-60. It uses a separately convergent full future and does not certify the actual Gamma source or global NS accuracy. Actual source history is admitted through the source bridge separately.

## Concrete next work

1. Own new steep_power_physical_C2.py/.json and checker companions. Reuse the unchanged general-K collar physical mapper with an original steep_power dispatch adapter and q=-wait-1-Ts*(1-phase). Preserve viscosity, lambda, amplitude, pressure datum and all source coefficients.
2. Complete Ttheta_theta=r*partial_z(Tz), reduced divergence and generally nonzero Etheta=-nu*partial_zz(utheta). Derive the pure-power homogeneous cancellations before enclosure. Do not import the waiting region's zero-remainder assertion.
3. Bind the actual source K_y=-k*K and B source factor. Useful reduced theta modes have rates -1.5 and -2.5; the axial pressure homogeneous divergence cancels after the B^2 factor. Prove the identities and endpoint map rather than relying on this planning formula alone.
4. Use one bounded native Cartesian fixture for any new physical operator terms; compare the complete tensor, pressure and retained nonzero remainder. Admit the physical right join only from common source functions, factors and their derivative identities.
5. Own steep_power_cone.py/.json and checker. Use the entire original Ts and Z[-1,1], kappa-2=2, actual full moments and exact positive S. Bound the original B maximum at q=-wait-1-Ts using common source log parts and cancel correlated Ts/wait terms before enclosure. Prove both directional inequalities; the positive theta box alone is insufficient.
6. Compose the power cone with the accepted steep-exit/waiting/collar tail only after physical/source joins and margins pass. Keep the full-tensor/global cone gate false.
7. Extend these same complete moments back through original steep entry, then angular repair, following power and the original 100-unit flatten. Recover signed stress/absolute-pressure derivatives, all functional joins and each regional physical/cone proof separately.
8. Finish actual Ra-to-R100 finite-width feedback and controlled omitted terms, independent global flat/physical-volume estimates and required-domain energy. Then establish common leading inputs before solving n=1 and n>=2 recovery equations, independent moment repairs and oscillatory stress cancellation.

Global admissible stress, global flat remainder, physical energy, true coefficient recursion, oscillatory correction and full corrected Cartesian NS validation remain unfinished.

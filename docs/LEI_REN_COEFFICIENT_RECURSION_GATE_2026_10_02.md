# Lei–Ren coefficient recursion: equations and admission gate — 2026-10-02

This note records the actual Part I coefficient hierarchy and the evidence
needed before a numerical `n=1` solve can be called paper-faithful. The source
is the local Lei–Ren v2 text `work_paper_cache/lei_ren_part1.txt`, equations
(13.15)–(13.27), Assumption 14.1, Lemma 14.2, and the regular system (14.19).
The corresponding equation locations are plaintext lines 21713–21811,
22273–22295, 22300–22383, 22619–22655, 23213–23266, and 24075–24111.

## Exact order dependence

For each order `n`, set

```text
e_n = 2*n*delta
b_n = -2-delta+e_n
c_n = -1-delta+e_n
p_n = -2-2*delta+e_n
```

The profile operators are

```text
T_a G = L^(-1) * (-a*G/2 + (1-delta)*Z*G_Z/2 + R*G_R)
Z_a G = L^(-1) * (a*Z*G + d*G_Z - 2*Z*R*G_R)
Z^[2]_a G = Z_(a-1+delta)(Z_a G)
L = 1-delta*Z^2, d = 1-Z^2.
```

For `n >= 1`, the coupled equations are

```text
2*(R*F_n,RR + 2*F_n,R)
  = T_(b_n) F_n
    + sum_(i+j=n) [V_i*(F_j,R + F_j/R) + Uz_i*Z_(b_j) F_j]
    - Z^[2]_(b_(n-1)) F_(n-1)

2*(R*Uz_n,RR + Uz_n,R)
  = T_(c_n) Uz_n
    + sum_(i+j=n) [V_i*Uz_j,R + Uz_i*Z_(c_j) Uz_j]
    + Z_(p_n) P_n
    - Z^[2]_(c_(n-1)) Uz_(n-1)

P_n,R = sum_(i+j=n) F_i*F_j - Omega_(n-1)/(2*R)

Omega_k = T_(e_k) V_k
    + sum_(i+j=k) [V_i*(V_j,R - V_j/(2*R)) + Uz_i*Z_(e_j) V_j]
    - 2*R*V_k,RR - Z^[2]_(e_(k-1)) V_(k-1).
```

Negative-index profiles are zero. At `n=1`, the known axial-viscosity
forcings are `-Z^[2]_(b_0)F_0` and `-Z^[2]_(c_0)Uz_0`; the pressure equation
still contains the unknown linear term `2*F_0*F_1`. Solve the three profiles
together with zero axis data. `Omega_0/(R)` must be evaluated through a
regular axis factorization, not by floating-point division near `R=0`.

## What qualifies as recursion

Lei–Ren Assumption 14.1 requires one fixed `R_in > R_a` shared by every order,
with `F_0, Uz_0, P_0, V_0/R` and the needed radial derivatives smooth on
`[0,R_in]` and holomorphic in a common neighborhood of the real `Z` interval.
Lemma 14.2 then gives the regular zero-axis order-`n` solution on that same
interval. The regular first-order realization uses
`xi=sqrt(R)`, `K_n=M^z_n/R-Uz_n`, and
`Y_n=(F_n,Uz_n,K_n,P_n,partial_xi F_n,partial_xi Uz_n)`.

After solving `n=1`, its actual five order-one moments must be repaired as
functions of `Z`, including the renormalized angular-viscosity moment and the
two stress-integral cancellations. Only then may those corrected profiles
source `n=2`. The same common radial interval and genuinely order-dependent
recovery coefficients must be retained. A radial Taylor series, coordinate
rescaling, or pressure-row-only implementation is not this recursion.

## Repository evidence and current gate

- `background_recurrence.py` implements an older OpenAI pressure/
  `Omega_k/X` row; its own docstring says it is not the coupled solver.
- `paper_core_series.py` and `lei_ren_part1_paper_core_recursion.py` are local
  radial Taylor constructions, not the temporal coefficient hierarchy.
- The v2 compliant core and shared five-moment repair provide substantial
  candidate-local leading-profile evidence. The current physical-field receipt
  still flags core/annulus interfaces, full-field smoothness, whole-background
  admissible stress lift, and `temporal_recursion` as incomplete.
- `lei_ren_part1_first_order_sources.py` is an older v1 extended-field source
  diagnostic. It computes sampled `Omega_0` and viscosity inputs but explicitly
  reports `first_coefficient_solved=false`; it cannot certify a v2 `n=1` solve.

Therefore no positive-order coefficient profile is admitted yet. The next
faithful implementation must first bind the v2 leading profile and its
`V_0/R` regular jets on a single common interval, complete the required leading
stress/cone and flat-remainder gates, and then solve the coupled v2 `n=1`
system. Preserve the distinction between a local candidate computation and
the theorem-level assumptions used for the all-order induction.

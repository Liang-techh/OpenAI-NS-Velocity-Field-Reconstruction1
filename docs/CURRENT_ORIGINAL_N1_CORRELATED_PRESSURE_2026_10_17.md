# Correlated actual leading-pressure increment

Tasks14/15 and the full long-term goal remain **ACTIVE / INCOMPLETE**. This stage computes a correlated enclosure of the actual first-collar leading N1p increment at the saved Z interval near0.371 and phase s=sc/4, including ordinary Z0 andZ1. Both P1_W0 increments are strictly positive; both differently weighted m15 leading-N1p contributions are strictly negative. Full n1 pressure, full C_inner, all-Z coverage and independent terminal matching remain open.

## One actual inlet function and exact mean removal

The new evaluator uses the accepted actual_first_inlet core Psi_rho jets, the checked same-source mean archive and the native transport identity. It does not subtract the separately boxed inlet/exit Qpsi values. The selected phase is the exact sc/4 used by selected_actual_inner_exit, not the unrelated historical s=1/2 sample.

Use the same physical conventions as the accepted parent:

```
Q = Qaxis(Z) + epsilon*b(rho,Z), b=Qpsi
U = Uaxis(Z) + epsilon*Psi(rho,Z), epsilon*Lambda=1
L = 1-delta*Z^2, d=1-Z^2
N1p = -S/2 + Lambda*(2*Q_rho+rho*Q_rhorho)
P0 = integral epsilon*N1p d-rho
J15_Np = -integral epsilon^2*rho*N1p d-rho.
```

P0 is the local W0 pressure increment, not the full feedback-dependent n1 pressure coefficient. J15_Np is only the existing leading-N1p source component of C_inner (canonical handles53161..53164), not a new ledger term. P1/m13 and m15 remain distinct data.

The actual own mean m obeys m+rho*m_rho=Psi, and

```
b = [2Z*Psi-(1-delta)Z*m-d*m_Z]/L
b+rho*b_rho = [2Z*rho*Psi_rho+(1+delta)Z*Psi-d*Psi_Z]/L
Delta(rho^2*b_rho)
 = [2Z*Delta(rho^2*Psi_rho)-(1-delta)Z*J-d*J_Z]/L
J = integral_4^rhoB rho*Psi_rho d-rho.
```

These are exact same-function FTC identities. Shared mean and Psi-value inlet constants cancel before boxing; the uncertainty of the single actual core inlet derivative pA(Z)=Psi_core,rho(4,Z) remains. Seven fresh symbolic checks cover the mean removal, m15 primitive through Z1, regular factors and flat-term signs. The separate sigma monotonicity proof checks two additional symbolic identities and positivity.

## Correlated radial transport with the original flat atom

The source identity on this first window is

```
w = Psi_rho = chi*R*Psi_core,rho(rho,Z)
chi=1-(1-hb)*sigma(s), R(0)=1
R_s=hb*rho*(chi-1)*(Phi_core,rho/Phi_core)*R.
```

Here sigma is the original leading-shear logistic switch. It is distinct from the own-velocity cutoff eta, whose selected endpoint value is zero. The exact nonzero atom Sigma=sigma(s) is retained as a fifth source factor:

```
odds=(1-s)^(-2)-s^(-2)
Sigma=exp(odds)*kappa, kappa=1/(1+exp(odds)).
```

The enormous negative odds is never exponentiated. Directed bounds enclose kappa and log(Sigma), prove Sigma>0, and use s*exp(-1000) only as a positive upper bound. The code checks the pinned source's logistic derivative and left zero branch. Its derivative is strictly positive on0<x<1: after x=y/(1+y), E=exp(odds)>0, it equals 2*(1+y)^3*(1+y^(-3))*E/(1+E)^2. Consequently0<=sigma(t)<=Sigma for0<=t<=s; this is required before bounding the transport integrals.

The accepted radial factors are retained:

```
rhoB-4=hb*Cd, Cd=4*s*exprel(hb*s)
integral_4^rhoB rho d-rho=hb*Cw, Cw=16*s*exprel(2*hb*s)
K=[(3+delta)Z*pA-d*pA_Z]/(2L), V=2Z*pA/L
P_regular=2*epsilon*hb*Cd*K
J_regular=-2*epsilon^2*hb*Cw*K
P_flat=-epsilon*Sigma*(1-hb)*rhoB*V
J_flat=epsilon^2*Sigma*(1-hb)*rhoB^2*V.
```

The signed axis inertia and its remaining S1/S2 volume radii are inherited from the accepted pressure-primitive stage. They are not replaced by bound midpoints. Ordinary Z derivatives use the actual pA jets through order2. Z-dependent F0 source factors are forbidden from the constant-factor arithmetic.

## Certified error and actual numerical result

B1 is the accepted complex core first-derivative cover. C is the saved nonlinear contraction size divided by1-q, and B2=32*C bounds the core radial second derivative using the accepted radius1/2 Cauchy domain. The ratio estimate gives |R|<2. With rhoMax=4.1 and the accepted core logarithmic-derivative bound Ac,

```
T=2*Cd*B2+2*s*rhoMax*Ac*B1*Sigma
|wB-chiB*pA| <= hb*T
Epsi=DeltaPsi-pA*DeltaRho
|Epsi| <= hb*Cd*Sigma*B1+hb^2*Cd*T
EJ=J-pA*integral rho d-rho
|EJ| <= hb*Cw*Sigma*B1+hb^2*Cw*T.
```

The factor2 on Cd*B2 includes the |R|<2 cover and was strengthened following independent review. Spending eta/8 for the required first Z derivative and eta/16 for the whole error's Z1 Cauchy bound includes derivatives of the complete error function and denominators. No complex-rho analyticity of the flat switch is asserted; its transport integration is along real rho.

New1800-bit directed arithmetic keeps hb, Pstar^2, Lambda and Sigma as source factors. After division by the stated positive physical scale, identical factors cancel before their logarithms are boxed. Each of34 remaining signed/error monomials is bounded in absolute value by2^-1024. Their sum34*2^-1024 (about1.892e-307) is a radius, never a signed function value. The uncertainty of the dominant actual inlet coefficient is retained.

The following rounded decimal displays summarize the directed enclosures; exact interval endpoints are serialized as mpf tuples in the report. The normalizing scales are hb*Pstar^2*epsilon*s for P0 and hb*Pstar^2*epsilon^2*s for J15_Np. They are Z-independent positive source scales.

| ordinary Z order | P0 / scale | J15_Np / scale | relative interval width upper |
|---|---:|---:|---:|
|0|approximately6.688947443079413668758621|approximately-26.75578977231765467503449|7.788e-29|
|1|approximately134.4941558700555848611644|approximately-537.9766234802223394446576|2.093e-29|

These certify the actual local leading-N1p integrals and their signs at the saved Z interval. They do not supply the missing cut-velocity/coupling source, complete C_inner, all-Z function or terminal conditions.

## Run, independent review and preservation

Production exits0 in 3.741 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_correlated_pressure_increment.json` is 10,436,545 bytes, SHA256 `895e772d41b0f0756648aa66c08c788e2cd93c096b50dfb7cd62d56a9499e3ee`, with 1986 input pins and 328 implementation pins. Focused Python compilation, saved-only inspection and the scoped raw-byte index audit pass. No accepted ancestor owner, source fixture, coefficient recurrence, old norm producer or old quadrature was replayed.

A fixed luna_max_worker performed the bounded independent mathematical review, identified the factor2 and explicit monotonicity gate, and accepted their resolution. The requested child configuration was GPT-5.6 Luna/max with its original fork_turns:none brief. The tool does not expose independent runtime model metadata, so none is claimed. Root integrated the corrections and final source gates.

Every inherited source/error packet remains equal at JSON-value level except stage/hash/contract/check/time metadata. Rc178/quiet179/raw213/reduced183, six signed core/lower contributions, signed slots19/20/41/42 and all finite-N volume caps remain unchanged. Target/quiet remain consumed once. This new correlated representation replaces no old ledger entry and removes no old error cap. It supplies a usable alternative evaluation of one existing component.

## Next executable work

Use this accepted correlated leading-N1p component when forming the same-source signed C_inner. Evaluate the remaining raw cut-velocity and2*F_lead*F1 coupling source with the same physical scaling and core-to-collar history; refine remaining inertia only where its error matters. Preserve each original source exactly once. Then extend the axial coverage required for functional claims. Do not repeat the endpoint-box subtraction or old accepted builders.

The independent reserved own_swirl/negative complete pulse m15 comparison and outer inlet p2/phase/H0/signed C/I1/I2/H work remain open. Accepted O2 E is not rebuilt. Functional all-five/common-N closure, exact heat/energy, admissible stress/flat remainder, true temporal n>=2 recursion, mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The1e-3 full-field target is not applied to this leading component.

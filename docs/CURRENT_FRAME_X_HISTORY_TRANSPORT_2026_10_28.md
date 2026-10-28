# Conditional ordinary-Z transport for the X histories — 2026-10-28

Goal ACTIVE / INCOMPLETE. A new single-entry leaf bound the original Gamma-history source and its saved receipt to the previously recorded `X`, `Xf` and `Xt` ports. It provides callable exact jet algebra and defining integrands through ordinary Z2. Eleven symbolic checks passed. It did not evaluate a source function or integral, traverse an expression graph, invoke an ancestor, or remove an error term.

Let `X` denote the inlet angular-history function, and set `r=1-mu`, `k=1-delta/2`, `L=-30 log(mu)`, `e=1/r`, `q=(1+Z²)/2`, `C=exp(-13r/mu)`, and `P=exp(-100r)/q`. With the original smooth switch `sigma`, define the **unevaluated finite integral**

```
J_f(Z) = integral_0^100 exp(r v) q(Z)^sigma(v/100) dv
D(Z)   = C X(Z) + (1-C)e + J_f(Z)
Xf(Z)  = P(Z) D(Z)
Xt(Z)  = A Xf(Z) + B
A      = exp(-k/2-r/2-rL)
B      = exp(-k/2) [exp(-r/2)((1-exp(-rL))e + I_in) + T_s + I_out].
```

These follow from the source definitions of `Xv`, the flattening primitive, `Xrel`, `Xin`, `Xs` and `Xt`, using `sigma(1)=1`. `I_in`, `I_out` and `T_s` are retained source-defined terms, not zeros or interval centers. `J_f` is the flattening integral; it is distinct from the unit steep primitive also named `J` in the source, and from the pulse moment integrals used elsewhere. The formula assumes `mu>0`, `r!=0`, source-parameter independence from Z and the original switch semantics.

For supplied same-source ordinary jets `X_j`, `J_j`, `P_j`, the implemented algebra uses

```
D0=C X0+(1-C)e+J0;  D1=C X1+J1;  D2=C X2+J2
Xf0=P0 D0;  Xf1=P1 D0+P0 D1;  Xf2=P2 D0+2 P1 D1+P0 D2
Xt0=A Xf0+B;  Xt1=A Xf1;  Xt2=A Xf2.
```

The finite-integral derivative integrands are `exp(r v) sigma Z q^(sigma-1)` and `exp(r v) sigma [q^(sigma-1)+(sigma-1)Z²q^(sigma-2)]`, with `sigma=sigma(v/100)`. Their formal first and second derivatives were checked symbolically. On a compact Z interval, the stated original smooth switch is bounded and `q>=1/2`; the finite integration interval then supplies a uniform integrable bound, so differentiation under this **finite** integral follows under those source-function assumptions. This is distinct from the later improper Gamma integral, whose transport is not proved here.

At the axis, `P1(0)=J1(0)=0`, so `Xf1(0)=2 exp(-100r) C X1(0)` and `Xt1(0)=A Xf1(0)` under the same assumptions. These are derivatives evaluated at the axis, not the ordinary derivative of the fixed scalar `Xf(0)` or `Xt(0)`. Their second derivatives retain both `P2 D0` and `J2`; they must not be dropped.

The callable algebra is conditional: the saved current references `[265772,265776,265782]` are ports mapped from the older X rows, but their actual common C2 function and ordinary derivative identities have not been proved. Nor have capture-free inlet rewriting, the actual axis composition, improper-Gamma differentiation, X4/X5 compatibility, or the physical frame seam been established. No numeric source values, certified error reduction, global five-moment closure or later physical layer is claimed. Existing errors, nonzero histories, independent reference and once-only reserved contribution remain; `P1` is distinct from `m15`. Runtime model metadata was unavailable (`null`).

Next: prove the current Rp source's `X=H/E` port and its ordinary Z1/Z2 rows are one C2 jet, then bind those actual rows to this transport and the axis compositions. Separately resolve whether current X4/X5 are semantically zero before relying on the old constant-X C5 assumption; do not infer function values from node IDs.

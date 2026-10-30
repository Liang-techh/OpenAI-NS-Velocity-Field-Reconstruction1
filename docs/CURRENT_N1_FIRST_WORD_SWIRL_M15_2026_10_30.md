# Same-source first-word swirl m15 contribution, 2026-10-30

The saved finite first word `W1` and the actual leading swirl source now give a directed-interval enclosure of its **first-collar swirl contribution** to physical m15 and ordinary `Z1`. This is one finite component under the parent source assumptions; it is not the full `W0+W1` prefix, pressure history, reserved own-swirl response, complete m15, terminal discrepancy or global five-moment closure. It is not appended to the existing 212-term ledger.

## Physical source and F0 connection

For raw physical `F1_W1`, `epsilon=Lambda^-1`, `s=(sc/4)q`, `rho=4 exp(hb*s)` and the original own-velocity cutoff `eta(q)=1-sigma(4q-2)`, the physical contribution is

```text
I_F,W1 = -2*epsilon^2*hb*(sc/4)
         * integral_0^1 rho(q)^2*F_lead(q,Z)*F1_W1(q,Z)*eta(q) dq.
```

Ordinary `Z1` uses `F_lead_Z*F1_W1 + F_lead*F1_W1_Z`; the physical cutoff and radial measure are Z-independent under the authenticated geometry. The negative sign and `epsilon^2` come from the original signed m15 density, not from a pressure endpoint or the reserved own-swirl response. The exact cutoff mass is `5/8`; it is used as a positive weight on whole-phase interval products, not as a point evaluation of the integrand. The radial multiplier is conservatively enclosed by `exp(2 hb s) in [1,1/(1-sc/4)]`; the resulting coefficient factor is `-5*sc` with external formal powers `hb^1 Lambda^-2`.

The radial receipt supplies ordinary derivatives of normalized `Phi`; physical leading swirl is `F_lead=F0*Phi`. The new source-factored connection applies the amplitude derivative exactly once:

```text
F_lead_Z = F0*(Phi_Z + Lambda*ell_coeff*Phi)
ell_coeff = -L*H/(H^2+(j/500)^2)
L=1-delta*Z^2, H=Z*(1-delta)/2+(1-Z^2)*(4Z+j).
```

`F0` and `Lambda` remain formal atoms, and `ell_coeff` is a directed interval, not a selected function value. Core and local `W1` branches remain separate. Grouping only performs outward interval addition within the same powers and uncertainty classes; it proves no same-source cancellation or sign of the physical field.

## Run and limits

The one-entry consumer authenticated the joined first word, leading radial source, prior axial contribution, core/transfer and source definitions, including the original negative swirl density and amplitude connection. It retained all 807 input IDs. For Z0 the core/local branches contain 64/158 grouped terms from 352/274 products; for Z1 they contain 128/316 groups from 1,408/822 products. The core-cell, source-uncertainty, Cauchy-phase, cutoff, mean-phase and F0-connection interval classes remain represented. These are interval sectors, not independent samples.

An exact-rational synthetic connection fixture passed 7 checks; an independent nonconstant synthetic swirl fixture passed 16 checks including direct high-precision cutoff integrals, physical powers and provenance. The latter's quadrature is a numerical cross-check, not a certified quadrature bound. The actual product enclosure uses whole-phase directed boxes and can be wide.

The successful receipt is `experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz`, SHA-256 `a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4`; the one-entry result SHA-256 is `5955cd84a89d5708392b993ae743b89ee62e262ec7280395da26d615349d35d0`. Its `F0_connection` index lists several operand labels in a single `path` field, which is not a traversable nested path. The metadata-only correction `experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json` (SHA-256 `f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4`) records the exact separate core and transfer paths; it changes no numerical result and reruns no calculation.

The earlier W1 axial component, Z2, Cauchy alternative, all-word tail and all historical m15 terms and errors remain separate. The saved receipts do not contain complete first-collar `W0_Z0/Z1` sectors: the former W0 local Green terms were internal intermediates and not exported. Next, form or fuse the missing same-source local W0 contribution with the physical m15 integrals without rerunning the old Z2 program; retain the separately required pressure history and endpoint. The corrected-reserved/current-pulse owner and measure bridge is still open; `P1` is not m15.

# Same-source first-word axial m15 contribution, 2026-10-30

The saved finite `W1` and actual leading axial velocity now give a directed-interval enclosure of **one component** of the physical first-collar m15 increment and its ordinary `Z1`: the axial coupling due to `W1`. This is conditional on the parent actual-source enclosures. It is not the full finite prefix `W0+W1`, the pressure or swirl contribution, a complete m15 function, a terminal discrepancy, or five-moment closure. It has not been appended to the existing cumulative ledger.

## Physical formula and cutoff

The joined `Uz1` is the raw physical n1 axial port, so with `epsilon=Lambda^-1`, `s=(sc/4)q`, `rho=4 exp(hb*s)` and `R=epsilon*rho`, its first-collar axial contribution is

```text
I_W1 = epsilon*hb*(sc/4) integral_0^1 2*rho*U_lead*Uz1_W1*eta(q) dq
eta(q) = 1 - sigma(4q-2).
```

The ordinary Z1 uses `U_lead_Z*Uz1_W1 + U_lead*Uz1_W1_Z`. The saved geometry makes `epsilon`, `hb`, `sc`, the phase and cutoff Z-independent. This velocity cutoff `eta` is distinct from the leading prescribed-shear cutoff and does not turn the pressure history off. The authenticated exact flat sigmoid obeys `sigma(t)+sigma(1-t)=1`, so `integral_0^1 eta dq=5/8` exactly. This identity is used only as a positive integration weight on **whole-phase interval enclosures**; it does not replace the integrand with an inlet value or select the midpoint of an interval.

The factor `rho/4=exp(hb*s)` is enclosed by `[1,1/(1-sc/8)]` using the admitted `0<hb<1/2`. Each resulting sector carries the formal external powers `hb^1 Lambda^-1` and the exact coefficient factor `5*sc/4`. The actual leading `U` and its ordinary derivative come from the same-source radial receipt; the core-transport and local branches of `W1` remain separate. Each product and interval addition is outward directed. Equal-power/error-class sectors may be grouped, but no sign, endpoint or same-source cancellation follows from that grouping.

## Run and limits

The actual one-entry consumer authenticated the joined finite-word and leading radial receipts, exact source atoms and phase, and the original physical measure, velocity cutoff and exact flat-sigma source definitions. It retained all 1,118 distinct input IDs with paths back to the parent receipts. For each of Z0 and Z1 the output has 128 core-branch and 231 local-branch groups; before grouping, Z0 contains 480 and 394 products, and Z1 contains 1,440 and 788 products. The core-cell, unclassified core source uncertainty, Cauchy-phase, mean-phase and cutoff error classes remain represented. These are interval terms, not independent samples or actual function values.

An independent nonconstant synthetic profile passed 16 checks, including direct high-precision cutoff integrals for both orders and branches, the cutoff mass, formal physical powers and error provenance. The fixture's quadrature is a numerical cross-check, not a certified quadrature bound. The actual range calculation uses the positive-weight interval enclosure; it may be wide, and no physical sign is inferred.

The successful receipt is `experiments/root_st073/n1_first_word_axial_m15_v1_receipt.json.gz`, SHA-256 `02acdcb22d9b0f03b4201f502388814d54d49937124ebb2288e66b9706770ab5`; the one-entry result SHA-256 is `a6fb10ca1b30d0dc81e0bba3b05eeb92ac95ef035a0110efe910a0c83a9a507d`. No ancestor program or previous Z2 calculation was rerun. The prior Z2 result, Cauchy alternative, all-word tail, 212 cumulative terms, six signed terms, reserved contribution and inherited errors remain separate and unchanged.

Next, construct the corresponding `W0` axial contribution and the physical swirl and pressure/endpoint contributions at ordinary Z0/Z1, keeping the exact physical scales and cutoff/pressure histories distinct. Only then can a nonduplicating finite-prefix m15 interface be assembled. The independent corrected-reserved/current-pulse source-owner, measure and moving-endpoint bridge remains open; `P1` is not m15.

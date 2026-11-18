# ST073 signed lower-collar m15 radial FTC — 2026-11-18

This checkpoint adds a function-level representation of the signed lower-collar pressure term. It uses the checked canonical source archive and keeps the old forward enclosures separate. The [previous reserved m15 C1 handoff](CURRENT_RESERVED_ACTUAL_M15_C1_2026_11_03.md) and the [Lei–Ren Part II auxiliary operators](CURRENT_LEI_REN_AUXILIARY_TORUS_2026_11_18.md) address different stages of the construction.

## Source identity

Use `rho(s) = 4 exp(hb s)` and `q = ratio`, with `q(0,Z)=1`. The saved log-ratio ODE gives `q_s = hb rho (chi-1) (Phi_core,rho/Phi_core) q`. Write `T = rho² q Psi_core,rho` and `T0 = 16 Psi_core,rho(4,Z)`. The signed singular coefficient of the lower-collar pressure satisfies

```
A  = (2Z/L)(1-hb) T
A0 = (2Z/L) T0
K  = rho² q [2 Psi_core,rho + rho Psi_core,rhorho
              + rho (chi-1) (Phi_core,rho/Phi_core) Psi_core,rho]
(A-A0)/hb = -A0 + (1-hb)(2Z/L) ∫₀ˢ K(t,Z) dt
```

The FTC for `T_s = hb K` uses the same `T0`, retaining `-A0` as the inlet term. [`original_collar_signed_ftc.py`](../src/openai_ns_reconstruction/original_collar_signed_ftc.py) defines the ordinary `Z0/Z1` pair with the quotient and product rules, keeping the physical factor `epsilon² hb sigma′` outside that normalized coefficient. The resulting signed primitive is its integral from `0` to `s`. The principal `epsilon² A0 sigma` and the existing core pressure history are separately referenced. Here `s=(sc/4)q_chart` with upper endpoint `sc/2`. The `ratio` variable above and `q_chart` describe different quantities.

## Actual checked run

The frozen consumer authenticated the raw lower-collar archive before decoding, checked the matching acceptance receipt, and bound the restored radial core source rows for `Phi`, `Phi_rho`, `Psi_rho`, and `Psi_rhorho` with their ordinary Z derivatives. It used the saved ratio pair `[35940,35941]`, the `rho` node `35851`, and the checked core `m15` primitive `[53125,53126]` as the incoming history reference. No historical producer or fixture was run.

| New local graph port | Ordinary `Z0/Z1` handles |
| --- | --- |
| `(A-A0)/hb` | `[145,146]` |
| Physical signed integrand | `[151,154]` |
| Signed primitive | `[157,158]` |
| Principal plus signed primitive | `[165,166]` |

The [passing v2 receipt](../experiments/root_st073/original_collar_signed_ftc_v2_receipt.json.gz) contains the 167-node expression graph and inputs. Raw SHA-256 `88c36dde29a4cc3eb38cb0800686dd6731159f20f127ff8fc2402c397331c26f`. [The first receipt](../experiments/root_st073/original_collar_signed_ftc_receipt.json.gz) records a runner lookup error in the source of the core incoming-history handle. This stopped before the leaf was called. The corrected reader ran as a separate frozen entry.

## Remaining task

Evaluate or certify this signed primitive using the same actual C1 source evaluator for the core radial derivative, the ratio and the lower pressure partition. Then prove the corresponding relation to the legacy finite-`N` inner/lower approximations and their error leaves before reducing any ledger. The new result defines an exact signed source expression. Its physical value, cross-owner identity to the forward 212 rows, independent terminal match, and global five-moment closure remain open. The independent negative full-pulse target and the once-only reserved patch remain the terminal comparison to perform.

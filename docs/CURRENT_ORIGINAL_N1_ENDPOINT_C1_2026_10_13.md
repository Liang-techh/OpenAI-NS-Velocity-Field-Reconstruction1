# Original n1 finite pressure endpoint C1 source — 2026-10-13

The actual lower `first_micro/sc_half` and upper `O3_power` pressure endpoints now have checked first-Z derivative data on all four whole axial cells. This removes the missing finite-endpoint source adapter from the n1 pressure C1 work. Complete pre-Rp and Rin-to-infinity derivative accuracy, actual n1 repair coefficients, repaired own M/V/P, five terminal identities and temporal recursion remain open.

The endpoint source uses the retained checked all-chart physical owner. It takes the actual `V_y` source rows through order three, the genuine `V[2]` Taylor row multiplied by 2 for `U_ZZ`, and `V_y[1]` for `U_yZ`. The common physical `Pstar` is applied once. Incoming finite-N increments vanish at the lower source-owned collar, while original nonlocal histories remain present. No magnitude cap is differentiated or selected as a signed field value.

For `L=1-delta*Z^2`, the original endpoint identities are

```text
qplus = q+q_y = ((1+delta)*Z*U+2*Z*U_y-(1-Z^2)*U_Z)/L
qplus_Z = ((1+delta)*U+(3+delta)*Z*U_Z+2*U_y
           +2*Z*U_yZ-(1-Z^2)*U_ZZ)/L
           +2*delta*Z*((1+delta)*Z*U+2*Z*U_y-(1-Z^2)*U_Z)/L^2
q_yZ = qplus_Z-q_Z
```

The upper endpoint has exact `U=U_y=U_Z=0`, hence `qplus=qplus_Z=0`, but its finite-N `q` and `q_Z` and their nonlocal histories are retained. In particular `q_y=-q` and `q_yZ=-q_Z`. The signed viscous boundary operators are `B(q+q_y)` for m13 and `-B(R*q_y)` for m15; the actual endpoint radius is Z independent.

The quartet in `experiments/root_st073` is `lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1` with `.py`, `.json.gz`, `_check.py` and `_check.json`. Its local gate is `current_original_n1_preRc_endpoint_C1_installed`. Production and independent checking covered all four Z cells, including a separate nonunit finite-radius fixture and the public query at caller precision 15. Root checked all 1506 receipt input pins against the saved bytes before publication. The global, repair and full-pressure gates remain false.

Reproduce this new adapter with `python scripts/run_original_n1_preRc_endpoint_C1.py`. This restores checked owners and runs the new producer/checker; it does not execute ancestor producers or their suites. Keep a retained owner in a live driver for subsequent source extraction.

- [x] Recover actual lower U_yZ/U_ZZ source rows.
- [x] Differentiate the finite quotient including its L derivative.
- [x] Preserve upper q_yZ=-q_Z and signed weighted boundary operators.
- [x] Check all four axial cells, fixed finite radii, one Pstar scale and public precision restoration.
- [ ] Replace the excessive pre-Rp m_ZZZ bound using the same original primitive source and selected N.
- [ ] Propagate that genuine derivative bound through q_ZZ; retain all previous-chart histories at seams.
- [ ] Consume the new q_ZZ and this endpoint source in differentiated pre-Rp reduced-S integrals.
- [ ] Keep compact-band source derivative terms, prove quiet-band Z independence and count each radial piece once.
- [ ] Compose inner-cut, lower-collar, pre-Rp and checked outer-pulse C1 data in fixed original pressure units.
- [ ] Bind genuine signed full five-defect C1 functions to the installed n1 linear repair operator.
- [ ] Compute actual alpha/beta functions and recover the repaired family's own M/V/P.
- [ ] Establish all five terminal moment identities before implementing genuine n2/temporal recursion.

A scratch Fourier L2 estimate for the same C3 primitive is promising, but has no accepted production/pressure-integration gate yet. Do not promote it to full derivative accuracy or a recursion result.

Full original goal **ACTIVE / INCOMPLETE**. See [the outer-pulse C1 checkpoint](CURRENT_ORIGINAL_N1_PULSE_C1_2026_10_13.md) and [the physical n1 linear operator](CURRENT_ORIGINAL_N1_LINEAR_OPERATOR_2026_10_13.md) for the already accepted adjoining inputs.

# Compliant inner, outer and exact heat family (F37/F38)

The distinct `c_epsilon=.001` family now continues from the admitted analytic
core through the reference join, functional five-moment repair, outer
candidate and exact Gamma heat component. The ordered inner and outer
rebuilds finish successfully. This completes F37-A/B/C and F38-A/B, but not
the coupled outer angular/pressure repairs, actual axial pulse selection,
or a globally matched background.

## Source and interval-envelope transfer

All new modules use the prefix `lei_ren_part1_paper_compliant_`. Legacy `.01`
sources and receipts remain unchanged. The new source is
`5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae`;
its pressure-datum enclosure is
`dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010`.
Source epsilon is `.001*delta`; core epsilon is `1/Lambda`. They are distinct.
The actual delta is retained; the analytic-tube upper bound is not used as
a replacement parameter.

`compliant_seed_inclusion` proves containment of fourteen pressure mass
boxes, the complex flatten envelope and all 148 scaled pressure input jets
in the completed old directed recurrence inputs. The exact stored units are
reproduced as `exp(2*logPstar-logLambda)*normalized_P0_jet`. The same delta,
mu, logPstar, j, unique axis anchor and F0 definition are retained.

`compliant_Cstar_envelope_transfer` also proves containment of all 148
`S=epsilon_core^2*F0^2` jets and the selected h box. It explicitly checks
that the new Kp is dominated by the old envelope (both 17 in P0/Pstar^2
units), retains KN=18000, and links the F36 finite-core sensitivity receipt.
Independent new analytic existence/uniqueness identifies the new solution
inside these input envelopes. The old local receipts remain old-source
interval envelopes on Z in [.49,.51]; they are not relabelled point solutions.
The whole-axis exit proof is separately recomputed analytically.

## Inner and functional matching results

- Normalized core Phi lower bound: 0.26538107638889.
- Combined exit C2 bound: 3.4956700879212e-22.
- Mixed C3 radial-tail bound at degree 144: 6.8421059309919e-14.
  The physical epsilon prefactor is essential; a raw Psi target is not proved.
- Frozen angular Q upper: 0.31845729166667.
- K1: 192922259423111113597; all thirteen width gates pass.
- Whole-axis constant-power exit SQ lower: 1.7604999996667.
- Reference-join SQ lower: 1.6114999994723; restoration |bw| upper:
  0.0012799288928393.
- Whole-Z five-defect C1 bound: 4.9017390830378e-23, or
  0.00035055790161485 times e_star.
- Five-bump coefficient C1 bound: 3.37239648913e-20; corrected |bw|
  upper: 8.9742822813401e-6.

The new five-moment family is
`3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894`.
The repair is a smooth implicit functional closure on Z in [-1,1], with
callable coefficient and partial-moment enclosures. It is not yet a complete
pointwise physical velocity evaluator. Reference-join checks include two
independent stress identities and seven interfaces; repair checks include
thirteen independent integrand, primitive and stress identities.

## Outer and exact heat components

The new outer inlet, slope/buffer continuation, trial axial pulse map,
angular flatten/steep/waiting candidate and heat component all consume this
same five-moment family. The .001 waiting root is recomputed from its own
source. Its refined enclosure width remains approximately .0063099; the
legacy .01 waiting root is not reused.

Checks cover twenty inlet identities/twelve interfaces, ten buffer
identities/eight interfaces, twelve pulse-map identities/five fixtures,
seven angular identities/twenty interfaces, and nine Gamma/PDE/heat-moment
identities. The pulse map is a trial c1(a,Z), c2(a,Z) recovery. The actual
energy-selected ap remains undetermined.

The exact heat factor is defined by the positive Gamma integral, with a
finite directed remainder. The formal inverse radius S=1/Rtail is retained;
exp(-1000) is only a proved upper bound. Neither H nor its positive deficit
is replaced by its limiting value. Angular, pressure and energy defects
carry their distinct physical prefactors, including S or a*S.

The Section 7 numerical gate now follows from the exact correlated source
rational epsilon/delta=1/1000, rather than independently dividing rounded
intervals at a threshold. This proves the numerical epsilon requirement;
it does not prove all Section 7 hypotheses. Incoming moments are not yet
matched to the heat component.

## Reproduction

From the repository root, using Python with the existing mpmath/SymPy
dependencies:

```text
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage inner
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage outer
```

Both ordered runs passed on 2026-10-01. `--stage source` reproduces F36;
`--stage all` runs the full source/inner/outer sequence, and `--list` prints
the ordered modules. Each receipt retains source hashes and narrow scope.

## Next dependencies, in order

1. Derive the actual angular terminal defect as a function of Z relative to
   Z=0, preserving flatten-history correlation and the mu^30 suppression.
   Combine it with the exact Gamma heat defect in a common Rrel normalization.
2. Recover two angular bump coefficients from the angular and pressure
   moments, including the pressure quadratic term. Use a uniform inverse
   and a smooth functional branch; restore the prescribed preheat P0 rather
   than redefine it from the unmatched heat field.
3. Compute the entire corrected future swirl energy, including steep,
   waiting, collar and infinite heat-tail contributions. Select ap and
   install its actual c1/c2 pulse, preserving positive future-energy targets.
4. Transport all five primitives through the corrected outer field, recover
   Ur and pressure from those primitives, and prove the terminal identities
   as functions of Z, interface smoothness, finite energy and whole-outer cone.
5. Build admissible stress and flat remainder; implement the genuine
   n-dependent coefficient recursion and smooth divergence-preserving sum.
   Oscillatory correction and independent full Cartesian residual checks
   follow only after those dependencies hold.

No global admissible lift, flat remainder, temporal recursion, oscillatory
stress cancellation or full NS residual completion is claimed here.

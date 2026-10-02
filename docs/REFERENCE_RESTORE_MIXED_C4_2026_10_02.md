# Actual reference and axial restoration mixed derivatives through four

The accepted actual Rsh histories now supply logR/Z mixed velocity, pressure
and all five physical primitive derivatives of total order at most four.
Coverage is the complete reference continuation Rsh..Rz, original axial
restoration Rz..eRz, and postrestoration eRz..Rm. The unpatched provider stops
at Rm; the accepted actual five-bump patch takes over there.

This extends the SAME selected core, Cstar, bridge/switch/Rsh source, actual
five-defect family and analytic P0. It does not reconstruct numerical point
coefficients or install a complete Cartesian field.

## Actual equations and physical derivative units

Let y=logR and E=V110-4Z. The actual axial velocity is raw V=4Z+E*alpha,
without normalization by Pstar. Alpha is one before restoration, the original
1-sigma(t) during restoration, and zero afterward. Its ordinary y derivatives
through four come from the original flat cutoff, not an artificial width.

The six actual centered histories obey

```text
(D+1) em = E*alpha
(D+8/5) eh = 0
(D+8/5) ek = E*alpha
(D+1) aa = (E*alpha)^2
(D+6/5) eb = 0
(D+1/5) ep = 0
```

Repeated differentiation retains every quadratic convolution and actual
incoming history. Five physical primitives and pressure use their original
R, R^(3/2)*Utheta and Utheta^2 factors. Ur uses (D+1/2)^k Q, because the
physical sqrt(R/2) factor must be differentiated before forming the grid.
Axial Taylor coefficients are converted to ordinary derivatives using n!.

Utheta and angular primitive rows are divided by the FIXED current-basepoint
Utheta, not by the moving function Utheta(Z). Relative amplitude Taylor jets
therefore remain in the Z derivatives. Positive amplitude sources are
retained as exact formal logarithms; exponential caps are enclosures only.
The axial squared-velocity term in Mztheta uses Pstar^-2, while the original
patch centered-energy inlet uses Am^-2. Original P0 remains unchanged.

## Functional interfaces

Rz reference/restoration and restoration-end/postrestoration packets reuse
the identical actual histories. Original sigma endpoint derivatives vanish.
The independent checker requires exact serialized endpoint grid equality.

For the Rm join, exponential terminal transport cancels the original Rh
defect normalizations algebraically. The actual inlet defects are

```text
d = (em, ek, eh, aa/Am^2-eb/2, ep/2) at Rm
```

These give the exact five patch inlet histories, including the -5/12 axial
energy and 5/2 pressure baselines. The original bump sources vanish on an
open neighborhood of Rm, so identical primitives, P0, swirl and primitive
RHSs imply identical mixed derivatives through four. Interval overlap for
Uz/Ur/P is only a consistency diagnostic, not the functional join proof.

The Rsh/reshape mixed interface remains uncertified until the actual reshape
mixed derivatives are supplied.

## Evidence and reproduction

- 720 actual velocity/pressure and 900 actual five-primitive mixed bounds
  cover twelve whole-interval, endpoint and interior packets.
- 135 independent derivatives compare against closed physical primitive
  integrals with nonconstant cutoff and nonzero histories. This finite
  fixture checks formulas; it does not admit the actual source.
- Five original primitive RHS identities, five physical radial prefactor
  identities, five Rh-to-Rm defect identities and five actual Rm inlet
  identities are checked symbolically.
- Ten flat cutoff endpoint values/derivatives, exact Rz/restoration-end
  packets and SAME actual source/family hashes are checked.
- The existing read-only GPT-5.6 Luna/max worker reviewed the formulas and
  identified source-hash and Rm proof bindings; both were included.

Run the new stage only with accepted prerequisites:

```text
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage restoremixed
```

The complete ordered pipeline has 96 modules. Producer and independent
checker are `lei_ren_part1_paper_compliant_reference_restore_mixed_C4.py`
and its `_check.py`; each writes its own source-bound JSON receipt.

## Remaining work

Actual long reshape, both microscopic switches and prescribed-shear bridge
still need mixed radial derivatives and their remaining functional joins.
Then build the inner/pre-O3 dispatcher and complete Cartesian assembly.
Required-domain energy, admissible divergence-form stress, independent flat
remainder, actual n-dependent recursion and oscillatory correction remain
unfinished. Full-field/stress/global-energy/temporal gates remain false.
Original unlocalized whole-space kinetic energy remains infinite.

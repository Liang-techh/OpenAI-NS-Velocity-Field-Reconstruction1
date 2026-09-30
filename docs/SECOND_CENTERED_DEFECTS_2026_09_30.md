# Common-source second derivatives through five-moment repair

## Implemented chain

The actual second-axial R110 source now feeds all 69 labeled centered defect
parts, the flat shape kernels, the incremental nonlinear five-bump inverse,
and the reference-plus-correction similarity field. Each stage retains value,
first Z derivative and second Z derivative in the P9/W2 pressure/width ring.
The same analytic preheat datum P0 and its two derivatives are preserved.

The flat shape composer accepts an explicit B_ZZ. Its second derivative uses
both I^(n+2)(B0) B_Z^2 and I^(n+1)(B0) B_ZZ, with the same delta^n/n!
weights as the retained expansion. The actual order N=11 therefore requires
kernel derivatives through order 13. Omitting B_ZZ preserves the original
first-order interface; missing second derivatives are not guessed.

The inverse retains the initial linear solve and ten incremental nonlinear
updates, with all three axial slots. The corrected field transforms each
centered source part back into physical moment coordinates and adds the
reference and full nonlinear partial bump map: 71 labeled physical parts.
It recovers Ur and Ur_Z from the fixed-R radial moment identity. Ur_ZZ is
unavailable because that identity requires third moment derivatives.

## Evidence and reproducibility

- `lei_ren_part1_paper_flat_shape_second_fixture.py`: independent finite-Z
  second derivative check, maximum relative discrepancy approximately 7.99e-20.
- `lei_ren_part1_paper_second_centered_component_defects_fixture.py`: scalar
  direct density integration at shifted Z, five-point second derivative
  stencil with step 1e-4. Maximum relative discrepancy among nonzero second
  derivatives is 1.75147e-14. Row 1 is exactly affine and uses a zero-derivative
  check instead. All labeled second parts reassemble exactly in working precision.
- `lei_ren_part1_paper_second_five_bump_field_fixture.py`: independent shifted
  scalar field checks for fields and moments, plus Ur_Z; maximum scaled
  discrepancy 9.64652e-32.
- `lei_ren_part1_paper_second_centered_component_defects_check.py`: actual
  Z=.3, 473 digits, P9/W2, 69 source parts. Maximum nominal second-Z terminal
  inverse residual after ten updates is 5.85616e-209.
- `lei_ren_part1_paper_second_five_bump_field_check.py`: actual Z=.3, x=1.25,
  all 71 physical parts and their derivatives recorded. The recovered Ur
  agrees with the shared stress-module recovery with zero retained-atom
  discrepancy; P0 and P0_ZZ are preserved.

Each script writes its adjacent JSON receipt. The actual checks require local
internally generated second-source caches. Generate the common source using
the commands in SECOND_AXIAL_PROPAGATION_2026_09_30.md before running the
centered check; that check writes the cache consumed by the field check.
The caches are restricted to their declared Z=.3 source and are not portable
evidence of another axial point.

## Remaining requirements

These are derivatives of a finite numerical construction, not certified
derivatives of the exact profile. Core truncation, angular/preheat quadrature,
RK propagation, restoration and flat quadrature still lack enclosed errors.
Neither the single actual point nor the fixtures establish a uniform Z bound.
The absolute inverse residual does not prove relative closure of extremely
flat rows or bound the infinite inverse tail. Retained pressure/width atoms
and labels must remain separate when reporting these quantities.

Next establish source and defect C1/C2 bounds on a declared compact axial
interval, propagate all numerical remainder bounds, and use these to justify
the inverse contraction and relative-flat hierarchy. Extending to |Z|<1
requires endpoint control. Finite-energy radial matching, the exact heat
exterior, admissible stress in all regions, temporal n-dependent recursion
and oscillatory correction remain open.

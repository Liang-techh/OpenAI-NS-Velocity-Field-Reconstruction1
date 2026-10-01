# Functional reference targets and the directed physical repair backend

The completed Lambda120 core remains an axial-center calculation. This update
implements functional reference targets and the directed physical algebra needed
to feed the existing paper five-bump map. It does not fit an alternative profile
or claim that the actual transition has been matched.

## Reference endpoint targets

`lei_ren_part1_paper_reference_endpoint_targets.py` implements (9.3) from the
authoritative accepted schedule. The canonical parameter hash is checked, so
obsolete reference receipts cannot silently supply Rh or delta. The API
`reference_targets(ctx, Z, order=3)` accepts a real center or interval of centers
in `[-1,1]` and returns ordinary axial Taylor coefficients.

The physical moment factors remain separate:

| Moment | Radial factor | Factored target |
| --- | --- | --- |
| theta | Rh sqrt(2Rh) | 5 uh / 8 |
| z | Rh | 4Z |
| theta_z | Rh sqrt(2Rh) | 4Z times the factored theta target |
| z_theta | Rh | 16Z² - 5 uh² / 12 |
| p | 1 | 5 uh² / 2 |

Here `uh=exp(13.5)/(1+Z²)` and `logRh=accepted_logRref-5`.
`physical_moment_jets` materializes the factors in arbitrary precision when
needed. The factored representation avoids converting the enormous radius to
ordinary floating point. The core scaled radius `s=4` is never used as Rh.

Verification covers80 independently differentiated coefficients at four centers,
their enclosure by the whole-axis target jets, and physical scale conversion.
These are target-function checks, not a whole-axis core or matching certificate.

## Exact change of the five physical moments

`lei_ren_part1_paper_five_bump_moment_map.py` accepts the actual F/Uz and an
integration backend for three angular functions f_i and two axial functions g_j.
The backend supplies `int R^weight product(factors) dR` on the common compact
band. It must provide the actual support, source and integration error data.

For `dF=sum h_i f_i` and `dU=sum h_j g_j`, the map keeps exactly:

- `dMtheta=2 int R dF` and `dMz=int dU`.
- `dMthetaz=2 int R (Uz dF + F dU + dF dU)`.
- `dMztheta=int (2 Uz dU + dU² - 2R F dF - R dF²)`.
- `dMp=int (2F dF + dF²)`.

`assemble_map` returns A and a symmetric Q; `evaluate_map` evaluates the full
double sum `A*h+Q(h,h)`. `control_jacobian` returns its exact derivative. Axial
Taylor/interval coefficients pass through the same algebra, including every
quadratic angular/axial cross term.

This is an independent physical identity backend for the existing
`lei_ren_part1_paper_five_bump_map.py`, which already implements normalized
(10.8) with the paper's compact bumps. The existing map has finite quadrature
data and does not enclose quadrature error. It remains the source of the actual
paper supports; the new polynomial fixture does not replace those supports.

## Connecting actual defects to the existing paper map

The existing paper control order is `(c1,c2,xi1,xi2,xi3)`; the physical assembler
uses `(xi1,xi2,xi3,c1,c2)`. Keep this distinction when connecting the data.
`physical_defects_to_paper_rows` transforms the actual five defects into the
existing row order `(z,z_weighted,theta,z_theta,p)`:

`d1=Delta_z/Rm`;
`d3=Delta_theta/(sqrt(2)*Rm^(3/2)*Am)`;
`d2=Delta_thetaz/(sqrt(2)*Rm^(3/2)*Am)-4Z*d3`;
`d4=(Delta_ztheta-8Z*Delta_z)/(Rm*Am²)`;
`d5=Delta_p/Am²`.

Am and Z retain their axial jets; freezing Am would produce incorrect derivative
data. `paper_rows_to_physical_changes` implements the inverse transformation.
The equation to solve remains the paper's `A*h+Q_Z(h,h)=-d(Z)`.

The independent fixture checks20 changed-field moment coefficients,100 control
Jacobian coefficients and20 physical/normalized conversion coefficients. Actual
transition defects, certified compact-bump integration and the repair solution
are still required. No functional five-moment closure is asserted.

Run both new receipts:

```powershell
python experiments/root_st073/lei_ren_part1_paper_reference_endpoint_targets.py
python experiments/root_st073/lei_ren_part1_paper_five_bump_moment_map.py
```

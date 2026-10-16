# Original leading implicit Rm patch - 2026-10-15

`OriginalRmPatchNumeric.coefficients(Z)` solves the original leading five-bump map from the actual signed unmodified Rm defect. `evaluate(Z,x)` recovers its leading fields, pressure, own moments and radial velocity coefficient for1<=x=R/Rm<=e; use x="Rh" for the exact source endpoint. All source-owned finite-N incoming histories and the positive-order n1 primitive remain separately attached. They have not yet been transported through this numerical patch.

## Original control coordinates and existence

Let s=10^10*t^2*Pstar^-2 and S=diag(t,t,s,s,s), with the original fixed t from the accepted repair. These are normalization coordinates, not new physical parameters. Write h=S*u. The directed calculation uses the exact conjugated map S^-1*(L*S*u+Q(S*u,alpha)+d). Its L and gg weights remain unchanged; fg, ff and ff_over_x acquire the positive angular unit s. The normalized axial nonlinear coefficient is exp(6/5)*(1+Z^2)^2/10^10. No angular floor is used as a divisor.

A strict self-map and contraction on the original normalized unit box enclose the actual leading root and its first Z derivative. The second Taylor coefficient is solved with the SAME actual base Jacobian:

```
J*h2 = -d2 - Q2(h0,h1)
```

All known lower-order convolutions and amplitude derivatives are retained. No missing Z3..5 defect rows are filled with zeros. The ordinary point root is a numerical proxy computed using full-support numerical bump integrals; its complete error is measured against the directed source/control enclosure. Midpoint data is used only to choose the fixed preconditioner for the interval proof.

## Partial moments, pressure and shared radial recovery

With d_i(x) the actual incoming defects plus the original partial bump-map change,

```
U = 4Z + axial_bumps
M = 4Z + d0(x)/x
theta = d2(x) + (5/8)*x^(8/5)
mixed = 4Z*theta + d1(x)
energy = d3(x) - (5/12)*x^(6/5)
pressure_moment = d4(x) + (5/2)*x^(1/5)
```

The generic S=Pstar histories are m=M/Pstar, h=am*x^(-3/2)*theta, k=am*x^(-3/2)*mixed/Pstar, e=(8Z*M-16Z^2)/Pstar^2+am^2*energy/x, and p=am^2*pressure_moment. Both the mass-defect1/x and centered-energy1/x are essential. The independent source P0/Pstar^2 is added to p exactly once.

The angular profile is Utheta/Pstar=am*(x^(1/10)+angular_bumps). Its native y=log(x) derivative and the axial y derivative are exported. Meridional recovery uses the same complete mass primitive:

```
Q = (2Z*U-(1-delta)*Z*M-(1-Z^2)*M_Z)/(1-delta*Z^2)
Ur = sqrt(R/2)*Q
```

The numerical packet supplies Q through Z1; similarity-to-physical time factors and Cartesian assembly are a later step. Native x is independent of Z, so M_Z includes the derivative4 as well as the true defect derivative.

For x>=71/40 all three original bump supports are exhausted. The code retains the unrefined full-weight map diagnostic, checks its enclosure, and uses the same unique implicit-function equation to refine the five LEADING terminal defects to exact zero. This is a functional identity of the coupled map; it is not a fitted endpoint reset and does not declare the finite-N moments closed.

## Focused execution receipt

```
python scripts/run_original_Rm_patch_numeric.py --write-report
python scripts/run_original_Rm_patch_numeric.py --inspect-report
```

The full focused run exits0 in176.607s; saved-report inspection exits0. Fourteen actual packets cover Z=-.3,0,.3,.17,.23 at x=1,Rh, plus x=1.25,1.5,1.75,1.8 at Z=.3. The last axial point is evaluated from the live source rather than an old saved sample. The maximum actual strict self-map contraction is7.3388740271418854543e-6, below1.

Checks include17 manufactured-known-coupled-control/Z2/contraction checks,34 independent raw five-density/radial-recovery integral coefficients,90 leading terminal five-identity coefficients,75 Rm inlet coefficients,38 leading/finite-N history distinctions and42 unchanged-P0/once-only pressure coefficients. The independent radial test integrates the original own-history ODE densities, including the distinct h/k3/2 rates, p rate0 and physical dx/x. It does not check closure by fitting sampled endpoints.

Read-only Luna/max review found no material algebra/proof/unit blocker. The new numerical receipt certifies the queried source branches and derivatives; it does not claim a new whole-axis numerical cover. The accepted whole-axis original-source branch theorem remains pinned separately.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_patch_numeric.json`,2716678 bytes, SHA256 `2144b8a14406dadc435e856452ee2761ba0f244624ff835baf26edb979199900`. It pins1760 source inputs and47 implementations. Accepted ancestor source/report bytes remain unchanged.

- [x] Actual original leading control Z2 and patch moment/pressure/radial-recovery reader, with focused source-current checks.
- [ ] Original finite-N density/history transport through this leading patch.
- [ ] Fixed-time physical cylindrical/Cartesian patch assembly and outer continuation.

## Remaining work

Continue the original finite-N signed density/history source through this fixed leading patch using dx/x, from the retained nonzero inlet. Re-solving the leading controls on modified diagnostic defects would change the accepted construction. The subsequent reserved repair remains a separate layer.

Then assemble the patch's fixed-time physical cylindrical/Cartesian values and extend outer/Rh matching. Higher-Z numerical defect/control rows, total radial/fast-phase derivatives, pressure joins and physical C2 remain open. No new whole-axis numerical certificate, global common-N admission, exact heat/energy/tails, admissible stress/flat remainder, genuine n>=2 recursion/smooth sum, oscillatory cancellation or corrected NS is claimed. Tasks14/15 and the full long-term goal remain active and incomplete.

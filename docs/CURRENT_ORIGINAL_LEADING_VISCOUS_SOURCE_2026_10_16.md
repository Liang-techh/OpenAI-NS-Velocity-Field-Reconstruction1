# Original signed leading-viscous source interfaces

Tasks 14/15 remain **IN PROGRESS**; the full goal remains **ACTIVE / INCOMPLETE**. This stage implements signed source operators, independent leading FTC recipes, 19 macro local coefficient evaluations, 38 first-switch endpoint evaluations and three phase-one source-integral C0 evaluations. It does not complete source-integral C1 or the five-moment terminal problem. The raw 213 and correlated 183 ledgers are inherited unchanged by their exact parent hash.

## Independent leading functions and physical scales

The independent high-radial background registry supplies the leading normalized functions. Modified physical Qy/RQy ports are not leading ports and are not aliased here. There are no archived leading Qy output nodes; four local endpoint recipes remain explicitly unregistered as graph nodes.

| Chart | V Z0..2 | m Z0..2 | Vy Z0..2 | Independent mean-y Z0..2 |
| --- | --- | --- | --- | --- |
| second_micro | 69298,69299,69300 | 69305,69306,69307 | 69540,69541,69542 | 121876,121877,121878 |
| frozen_macro | 72142,72143,72144 | 72149,72150,72151 | 72382,72383,72384 | 133790,133791,133792 |
| first_switch | 87584,87585,87586 | 87591,87592,87593 | 87824,87825,87826 | 141159,141160,141161 |

All nine independent mean-y rows have the actual graph form sum(V,-m). All 27 normalized input rows retain their original inverse-Pstar factor node 69245, whose denominator is Pstar node 33572. Physical leading pressure therefore restores Pstar once. Physical y=log(R/110), d_y=(1/J)d_native and R_y=R retain the actual chart Jacobians. The macro J is Y-2hbB; the switch J is hbS. Native query node 48 is value 2, not node 2 (value 40).

Let L=1-delta*Z^2 and use ordinary axial derivatives. The independent pressure recipe is

```
q_y^L = [2Z Vy-(1-delta)Z(V-m)-(1-Z^2)(V_Z-m_Z)]/L,
H = R*Pstar*q_y^L.
```

The local H recipes use macro0, macro1, first0 and first1 scopes. The accepted independent V3/m3/Vy2 seam theorem and common physical-radius theorem give Hmacro1=Hfirst0 as source functions. Their proof is separate from the already-resolved finite-N D=R(Qy_modified-Qy_leading) category.

## Signed FTC and mean-ODE partition

The pinned physical density/FTC source contains Py=...+Qy+Qyy, N1p=Py/R and m15_y=-R*Py. Its viscous part is -d_y(R Qy). Hence

```
Vmacro+Vfirst = (Hmacro0-Hmacro1)+(Hfirst0-Hfirst1)
              = Hmacro0-Hfirst1.
```

For physical B=U_y^L, define E=2Z*R*B/L and C(B)=[(1-delta)Z*B+(1-Z^2)B_Z]/L. The explicit independent mean ODE gives

```
Hlower-Hupper = Elower-Eupper + integral R*C(B) dy.
```

This derives the upper negative sign, lower positive sign and positive source primitive. Absolute mean histories remain in H. The source control makes Uy zero at switch phase one; it does not make H, Qy or the mean zero. The new independent fixture includes this nonzero-history counterexample.

Raw 23/24 are positive bounds for actual-versus-frozen-proxy H discrepancies. Raw 43/44/46 contain signed endpoint/source-primitive approximations; raw 45/47 are radius/weighted-primitive remainder bounds. The FTC does not turn these caps into signed residual functions. No old term is removed, zeroed or replaced in this stage.

## New evaluated interfaces

The shared ordinary-derivative kernel computes 2Z*B/L and its Z1 row from signed B/B_Z. It separately implements C(B) C0 and requires B_ZZ for C(B) Z1; missing rows raise rather than being zero-padded. Inputs retain their saved numerical approximations and exact directed source enclosures. Recomputed approximations are never labeled exact source values, and cap/log atoms are rejected as signed coefficient inputs.

All 19 archived macro queries now have new local source coefficient C0/C1 rows. They are local q source evaluations; they are not the missing transported second-micro lower endpoint. The old interval-context transfer copies exact endpoints and applies no extra macro J factor.

All 38 first-switch upper/lower endpoint coefficients were recomputed. The 32 nonempty cases agree with independently saved signed point decimals at fixed relative tolerance 1e-70 (maximum observed endpoint relative difference about 4.57e-85). Six empty increment comparisons are explicitly skipped because the old producer stores an empty interval increment as zero. This does not zero its local source. Source-error radii and interval overlap do not supply the point-comparison tolerance; interval overlap is informational only.

The numerical endpoint rows carry R100*epsilon times hbS (upper) or hbB (lower). Their R100 reference is distinct from actual R(s); raw 45 retains the upper-radius difference. The original modal radial dependence inside B remains present. Physical R occurs once as the external moment weight, epsilon=Lambda^-1 remains factored, and hbB/hbS remain distinct source nodes.

The pinned pulse_masses and pulse_control_masses ASTs both define mass(0)=0 and mass(1)=1/2 exactly. No ancestor mass function was called. At phase one, the new kernel computes C(B) C0 times the exact half, with three independent raw46 C0 comparisons passing at the same point-only tolerance. These are signed integral-component approximations with source enclosures and the retained weighted remainder; exact mass does not make the whole source integral exact. Empty source-integral increments retain their own exact zero scope.

## Actual execution and open inputs

New production exits 0 in 5.452 s; serialized inspection exits 0 without constructing an owner, fixture or ancestor. The new report is 1,697,875 bytes, SHA256 `d1aa3644ea8158fe78dccb210a0a75a4c81ad257b8504758cb5797845b8003d4`, with 1957 input pins and 210 implementation pins. It references the unchanged 56 MB parent by SHA rather than duplicating its ledger. The 34 independent controls pass; maximum supporting difference is 5.7149369564113749111e-101. Tests differentiate nonconstant polynomial/rational source operators independently and verify the signed FTC with an analytic nonzero mean history. They do not validate an original Poisson inverse or Navier-Stokes residual.

The first-switch packet deliberately keeps complete=False and full-scope passed=False. Its advertised endpoint and phase-one C0 interfaces pass independently. There are one macro lower-dy blocked record and 35 switch blocked records. The source graph contains Vy_ZZ as a function, but the factored numerical B_ZZ row is not serialized. Interior mass values and weighted_primitive_remainder_per_hb are also not yet recovered. These absences are explicit; the saved raw46 output is never inverted to manufacture its missing inputs.

The inherited raw213/reduced183, independent negative complete-pulse m15 target, six signed inner/lower contributions and reserved-once entry have exact digest/count witnesses at all three saved Z queries. All prior error budgets remain active. The largest remaining compact/Rm finite-N volume bounds are inherited unchanged and are not signed residual values. Total signed m15, functional five-moment terminal closure, exact heat/finite energy, admissible stress/flat remainder, genuine temporal n>=2 recursion/summation, mean/oscillatory corrections and independent corrected Cartesian NS remain incomplete.

## Next executable task

Implement a bounded source-jet evaluator/export for the factored numeric B_ZZ row from the same original leading generator; graph Vy[2] is already defined, but that is not a numerical coefficient. Recover second-micro lower_dy from its actual/comparison Phi/Psi radial rows and rho. Then recover interior mass and weighted primitive source/remainder rows, compute source-integral C1 and compare signed coupled contributions with the independent m15 target. Keep the source function, numerical approximation and certified remainder separate. Do not replay accepted ancestor constructors/producers/calculations, use cap centers as values, delete source errors from an FTC formula alone, or claim five-moment closure from sampled comparisons.

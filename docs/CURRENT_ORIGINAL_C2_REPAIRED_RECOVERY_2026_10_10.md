# Actual repaired radial and inertial C2 recovery (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [e92a8ab5](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e92a8ab54c692fc9d94baa34d6406968cc67531f) consumes the accepted C3 repaired-band E/V, full five histories and independent P0 to install actual radial Q/Q_y, all four signed inertial sectors, generic shear/inertial numerators and quotients, and original cylindrical velocity/pressure functions through ordinary Z2. Strict same-source L/E/C positivity is now proved on all four axial cells and the whole reserved x[1,2] band. Directed recovery magnitude ranges are installed. Higher radial stress derivatives, Rh/absolute heat matching, the complete physical Cartesian field and temporal recursion remain open.

## Actual recovered functions

The source inputs are the corrected complete histories, including their original nonzero incoming memories. The ordinary derivative shift consumes the accepted C3 row directly: derivative(q) has rows (q_Z,q_ZZ,q_ZZZ). It introduces no additional Taylor factorial and uses no leading-only replacement for a repaired high row. Existing normalized E/V/history/P0 handles are retained exactly.

With z=Z, d=1-z^2, L=1-delta*z^2, and P=P0+p, the unchanged source recovery uses

```
transport=(1-delta)*z*m + d*m_Z
Q=(2*z*V-transport)/L
Q_y=(2*z*V_y-(1-delta)*z*m_y-d*(m_y)_Z)/L

theta_linear=(-E+(1-delta/2)*h-(1-delta)*z*h_Z/2)/L
theta_quadratic=((2*delta-1)*z*k-d*k_Z+E*transport)/L
axial_linear=(-V+(1-delta)*(m-z*m_Z)/2)/L
axial_quadratic=(V*transport+2*delta*z*e-d*e_Z+2*(1+delta)*z*P-d*P_Z)/L
```

The same physical source radius and scale restore It=R*(theta_linear+Pstar*theta_quadratic) and Iz=R*(axial_linear+Pstar*axial_quadratic). The signed pressure and every meridional, linear and quadratic term are retained. C=E-2E_y and B=2V_y supply the full source numerators den=C*E, kap=C^2+B^2, excess=kap-2den, H=C*It-B*Iz, D=H-kap and J=C*Iz+B*It.

Actual C2 quotients a=C/E, b=B/E, t0=-B/C, p1=It/E, p2=Iz/E, kap/den, D/den, J/den and Delta=kap/den-2 are installed. The complete generic discriminant (2*D^2*den-excess*J^2)/den^3 is a signed function with a directed magnitude bound. Its cone sign is not admitted here.

The original cylindrical source factors are

```
Utheta=Pstar*E
Uz=Pstar*V
Ur=Pstar*sqrt(R/2)*Q
Ur_y=Pstar*sqrt(R/2)*(Q_y+Q/2)
Pi=Pstar^2*(P0+p)
```

These are source-coordinate functions before the physical time and Cartesian map. They are neither a complete global dispatcher nor a numerical point oracle. The exact normalized radial divergence identity follows from m_y=V-m:

```
L*(Q+Q_y)=(1+delta)*Z*V-(1-Z^2)*V_Z+2*Z*V_y
```

## Strict denominator margins and parameter identity

The actual positive terminal amplitude A factors both E=A*(power+F) and C=A*((1+2alpha)*power+F-2F_y), where alpha=1/2+mu and power=x^-alpha. The unchanged selected repair N gives strict whole-band |F|,|F_y|<1/8; 1/2<=alpha<2/3 and x[1,2] imply power>1/2. Therefore

```
E/A>3/8>1/4
C/A>5/8>1/2
L>=1-delta>0
```

Directed lower logarithms are log(E)>=log(A_min)-log4 and log(C)>=log(A_min)-log2. The same actual A is factored before bounding, preserving its correlation. No interval endpoint is selected as a function value. These margins admit the local positive quotients, not the signed generic stress cone.

The pressure-source family fixes Md=40, logPstar=exp(40)+11 and delta=min(1e-200,exp(-4logPstar-30)). Every live source cell binds the reference logP exactly to flow.logs[1]/2, proves the exponential branch is below 1e-200, and checks exact tuple equality of the live delta with that branch in the same context. The original reference constructor and whole-Z source-to-reference owner are AST-bound. Thus the function graph and denominator proof use one exact source delta; overlap of unrelated enclosures is not used to establish identity.

## Radius, units and directed ranges

The correct source radius is R=Rc*x with Rc=op.Rm_factor*exp(exp40+20). The formal original_Rc_offset already includes the Rm log offset, so it is not exponentiated and multiplied by Rm a second time. The relative Rm radius offset is logPstar+9+log(x); its identity with the original absolute graph offset is independently checked. Large radii and field magnitudes stay formal or logarithmic.

Common V/m/k were already normalized by the accepted band source; this recovery does not divide them by Pstar again. P0 remains separate from p. Directed bounds use the accepted complete third history/first-y rows, original P0 thirds, actual R/Pstar source witnesses, full product/quotient recurrences, and the strict positive lower bounds above. Function definitions and magnitude caps remain distinct.

The dedicated evidence passed:

- 108 independent signed ordinary value/Z/ZZ rows covering Q/Q_y, pressure, all inertial sectors/numerators/quotients/discriminant and cylindrical source factors.
- 108 directed finite polynomial-jet comparisons against independently differentiated formulas, including pressure, quotient cross terms and full signed sectors.
- 78 exact input derivative aliases from the accepted C3 repaired E/V, histories, first-y rows and independent P0.
- All 4 actual axial cells, 16 radius/scale/delta source bindings, exact live delta/logP tuple guards and the original absolute/relative radius identity.
- Same selected-N whole-band positive margins and the exact own-moment radial divergence identity.

The accepted default constructor succeeded on the reused live owner. The exact Git-index audit passed all 1201 receipt dependencies. The sole reused static reviewer remained GPT-5.6 Luna / max; it found and helped close the delta provenance guard, with no remaining material issue. Producer gates remain false; the separate checker admits only the repaired C2 recovery, local positive denominators and directed magnitude bounds.

```python
from lei_ren_part1_paper_compliant_current_original_C2_repaired_recovery import CurrentC2RepairedRecovery

recovery = CurrentC2RepairedRecovery()  # reuse the live accepted band/owner when available
Q, Q_y = recovery.radial_functions()   # ordinary value/Z/ZZ
cylindrical = recovery.cylindrical_functions()  # Utheta, Uz, Ur, Ur_y, Pi
```

Installed quartet: `lei_ren_part1_paper_compliant_current_original_C2_repaired_recovery.py`, `.json.gz`, `_check.py`, `_check.json`.

## Completed bounded tasks

- [x] **CURRENT-C3-REPAIR-BAND-PROFILES** Actual y0/y1/y2, Z0/Z1/Z2/Z3 normalized corrected velocity functions.
- [x] **CURRENT-C3-REPAIRED-HISTORIES** Signed C3 partial/complete histories, independent P0, first-y rows and relative terminal identities.
- [x] **CURRENT-PATCHED-RECOVERY-SOURCE** Actual repaired radial Q/Q_y, full signed inertial sectors/numerators/quotients and cylindrical source factors through ordinary Z2.
- [x] **CURRENT-PATCHED-RECOVERY-POSITIVE-DENOMINATORS** Whole same-N repaired band strict L/E/C margins, same-source delta/logP identity and actual Rc*x radius.
- [x] **CURRENT-PATCHED-RECOVERY-C2-RANGES** Whole-band directed magnitude bounds for the full repaired recovery, pressure and physical source factors on all four axial cells.

## Next tasks, in dependency order

- [ ] **CURRENT-PATCHED-RADIAL-Y4-CONTRACT** Own a companion producer/checker for corrected velocity y0..4 and mixed Z0..3. Bind original beta third/fourth y derivatives and compact endpoint flatness to the existing original beta source; include the e^-y and ell/J0 factors before differentiation. Generate the required five ordinary radial rows for the unchanged `current_generic_shear_moment_recovery.py::field_rows`, with exact y0..2 aliases preserved. Repeated history y derivatives must come from signed own-rate FTC, and P0 must have zero positive-y rows.
- [ ] **CURRENT-PATCHED-MIXED-RECOVERY** Use those actual five profile rows and complete moment FTC rows to recover radial y0..4 and inertial y0..3 through available Z2. Differentiate sqrt(R/2) and R before applying physical derivative rows. Keep explicit derivative-order contracts; do not claim recovered Z3 from C3 histories. Obtain all actual mixed source bounds before local stress admission.
- [ ] **CURRENT-PATCHED-RECOVERY-Z3-PREREQUISITE** If stress assembly needs recovered Q/inertial Z3, build genuine fourth-Z source/phase/target/limit/history/P0 rows and their quantitative bounds first. Preserve the existing branch and all lower handles; no leading-only fourth append.
- [ ] **CURRENT-C2-RH-CORRECTION-JOIN** Continue actual repaired second histories through Rh with the same phase/radius and nonzero memory. Reuse the accepted leading join separately from the correction continuation. Prove interface values/derivatives and source geometry; do not use overlap as functional equality.
- [ ] **CURRENT-C2-ABSOLUTE-FUTURE-INTEGRALS** Recover source-owned absolute exterior five identities through Z2 with true heat amplitude, independent P0, future-integral FTC/Gamma tails and original cancellations. Retain the exact final reference power and distinguish it from the current mu power. Relative band zeros alone do not close absolute exterior moments.
- [ ] **CURRENT-C2-PRESSURE-HEAT-ASSEMBLY** Assemble actual signed pressure restoration and exact/controlled heat exterior, source-derived axis/interface regularity and finite-energy radial tail bounds.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Supply all chart/interface mixed derivatives needed by global physical velocity and stress; record exact orders and one shared source family. The current recovery is local to the reserved band.
- [ ] **CURRENT-PHYSICAL-DISPATCHER-AND-ORACLE** Assemble one core-to-heat similarity/physical time/Cartesian dispatcher returning u/v/w and axis limits. Combine phase/Picard/integral errors; select one finite frequency satisfying physical/matching/stress/numeric constraints. Current N remains repair-only.
- [ ] **CURRENT-STRESS-CONE** Use complete absolute histories and required mixed derivatives to recover divergence-form stress and flat remainder separately. Prove signed cone margins over exit/matching/pulse/end/flatten/collar and quantify volume norms and scale dependence. The new discriminant is a function, not a positive-cone receipt.
- [ ] **CURRENT-TEMPORAL-RECURSION** Implement original n-dependent coefficient recovery and independent moment repair on the common core domain, curl-preserving truncation and finite-order/smooth-sum remainder. Spatial derivative recovery is a prerequisite, not temporal scale recursion.
- [ ] **CURRENT-OSCILLATORY-CANCELLATION** Install both original pulse families after accepted background stress/recursion, with averaged quadratic momentum-flux cancellation, common frequency hierarchy and flat forcing remainder.
- [ ] **CURRENT-DYNAMICS-AND-CARTESIAN-RESIDUAL** Measure contraction, axial aspect ratio, swirl/vorticity amplification, true material winding, interscale recurrence and finite energy on the completed time-dependent field. Independently evaluate Cartesian divergence and the full forced NS residual after oscillatory correction.

Read this handoff first. Continue the radial y0..4 contract and C2 Rh/absolute heat matching; completed C3-band and C2 radial recovery tasks remain complete. Full reconstruction stays ACTIVE / INCOMPLETE.

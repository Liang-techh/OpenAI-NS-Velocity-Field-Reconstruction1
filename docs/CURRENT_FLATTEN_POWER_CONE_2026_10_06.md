# Whole current flatten and outer-power signed cone

Implementation and scoped whole current flatten/power cone receipt: commit [caa9f74f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/caa9f74f851391267adb0ffaacfdf98f30e5a8d0).

## Constructed result

**F57C-cone1b-flatten is complete.** The same current source graph now satisfies the strict original two-vector cone on flatten offset[0,100] and outer-power phase[0,1], for every Z[-1,1]. Together with main, exit, gap, gap-end and pulse-end, seven consecutive regions have whole-domain source-function cone bounds. The other 26 registry regions and the global completed-tensor admission remain open.

The implementation preserves full forward X, variable axial K and K_Z, current complete remaining energy, absolute pressure and all tensor/divergence/Cartesian remainder records. It consumes the current pulse-end/flatten, flatten/power and power/angular functional tensor joins. No new matching data or independently chosen velocity/pressure is introduced.

This advances admissible-stress stage 4. It does not implement actual oscillatory correction families or higher-order temporal coefficient recursion. The leading origin remainder still has its previously certified nonflat term.

## Correlated source bounds

Let a=delta/2, r=1-mu, k=1-a, b=(1-2a)/2, g=(mu-a)/r>0, f=(1+Z^2)/2, rho=log(f), and j=2Z^2/(1+Z^2). Original sigma is evaluated at t/100; sigma_t is its ordinary t derivative.

For flatten, F=f^sigma, N=F*X and A=K*X are the actual current functions. The full angular inertial numerator divided by K is

`W=((k+b*j)*N-b*Z*N_Z)/F-1`.

Set V=W-g. The original radial ODE gives

`V_t+(r+rho*sigma_t)*V=b*j*(1-sigma)-rho*sigma_t*(1+g)>=0`.

The signed current pulse memory is retained through its exact scalar recipe

`Xv=1/r+(Xp-1/r)*exp(-13*r/mu)`.

The directed bound on its contribution is `2*k*(abs(Xp)+1/r)*exp(-13*r/mu)`. It is proved smaller than `g_lower*exp(-1000)` before forming the positive uniform W lower. The broad Xv interval near 1 is never subtracted from a rounded equilibrium to infer a sign.

For outer power, with u=(Lrel-4)*phase,

`W=g+Hf(Z)*exp(-r*u)`.

The same normalized flatten integral gives the strict uniform lower

`Hf>=exp(-100*r)*(b*expm1(50*r)-2*k)/(2*r)>0`.

This uses original sigma<=1/2 on the first half, f in[1/2,1], the exponential/log tangent inequalities, and the correlated axial integrand. Independent Xf and Xf_Z boxes are not used to prove Hf positivity.

The complete theta source is `K*(W/L+2/R*(K_t/K-1-a))`, L=1-2aZ^2. The lower bound uses 0<L<=1 and W>0 on the entire Z domain. The exact inverse radius is positive and equals S_exact*exp(-q)=exp(-logR). A cap relative to g bounds its negative shear contribution; the Scap enclosure is not a defining field value. Minimum K follows from F/f=f^(sigma-1)>=1 and a-mu<0.

Current axial velocity/shear is structurally zero by the selected terminal histories and full FTC propagation. Axial stress remains nonzero and uses full E/P. Its whole current signed interval supplies an absolute upper bound, retaining the positive Qz factor during physical conversion.

Current source logs already contain the exact KR cancellation. With logU from the same current inlet, the uniform B maxima are

`flatten: logPstar+logU-13/(2mu)+(a-mu)*(100+Lrel)-13-log(2)`;

`power: logPstar+logU-13/(2mu)-100*(1/2+mu)+(a-mu)*Lrel-13-log(2)`.

There is no second KR multiplication. The final comparison is `2*Ctheta^2-(vs-2)*B^2*Cz^2>0`. Positive constant-viscosity and lambda factors preserve this cone. Full pressure remains `-positive_B_squared*P`; it is not confused with the cone pressure coordinate.

## Reproducible entry points and evidence

Use prefix `experiments/root_st073/lei_ren_part1_paper_compliant_`:

- `current_flatten_power_cone_operator.py`: 26 generic source algebra identities, 12 current source/history/log/radius identities, 32 strict directed margins, and an eight-interval continuous sigma derivative cover.
- `current_flatten_power_cone.py`: `CurrentFlattenPowerCone(gapcone=checked_current_gap_cone)`. The constructor obtains both whole views itself. `native(...)` admits flatten/power separately and preserves nested earlier pulse admissions.
- `current_flatten_power_cone.json` / `_check.json`: actual current source composition, bounds, scope and hashes. Focused rejection cases cover foreign source, partial domain, missing pressure and incomplete history.
- `current_flatten_power_cone_views.json.gz`: deterministic gzip with both complete unpruned source views. Decode with `json.loads(gzip.decompress(path.read_bytes()))`.
- Focused controller stage: `currentflattenpowercone`. Reuse the warm graph; avoid cold reconstruction and inherited producer reruns when the source is unchanged.

Producer, checker, focused controller, checked regional API and inherited scope passed. Read-only GPT-5.6 Luna/max review identified the explicit A=K*X binding; it is now bound to the actual full-row AST and its flatten W identity. Input hashes matched all793 dependency files in both the working tree and Git index before publication. Graph inventory remains 33 regions / 32 adjacent / 14 internal tensor traces, with primitive atlas 14 / 8.

## Next executable tasks

- [x] **F57C-cone1b-flatten:** complete actual pressure/energy/forward X/K_Z, signed current pulse memory, correlated Hf lower, exact shear and both whole domains with all three source joins.
- [ ] **F57C-cone1b-angular-1:** start from checked `registry.owners['angular']` on original s[-4,0], Z[-1,1]. Retain both selected beta controls and full quadratic terms. Derive the source-correlated angular inertial comparison and bind actual ordinary s rates. Do not fit three velocity components or substitute a cap endpoint as a control.
- [ ] **F57C-cone1b-angular-2:** bound all signed angular/axial stress contributions using current full E/P. Preserve the four internal support traces and power/angular and angular/entry joins. If a raw interval crosses zero, separate a signed main term and correlated errors before subdividing. Report a real sign failure separately from an inconclusive bound.
- [ ] **F57C-cone1b-tail-1:** current steep entry/power/exit and waiting, using their checked source maps, complete moments and original rate changes. Port generic cone identities onto current objects; recompute actual margins and consume current tensor joins.
- [ ] **F57C-cone1b-tail-2:** current heat collar and exact Gamma exterior. Keep absolute pressure and higher-pressure identities. Prove strict cone where stress is nonzero; treat the zero exterior and limiting uniform edge direction separately.
- [ ] **F57C-cone1b-entrance:** whole pulse entrance, O3 power/slope and O2/reference regions. Retain real incoming moments, pressure and variable log slopes; publish any sign failure with its exact source sector and domain.
- [ ] **F57C-cone1b-inner:** actual patch/restore/reshape/switch/bridge/core. Where required, construct the original periodic shear loop, finite uniform N and independent five-moment repair on reserved intervals. A source-family change requires dependent tensors, joins and receipts to be rebuilt.
- [ ] **F57C-cone1c:** construct two actual homogeneous oscillatory pulse families, their covariance integrals and finite errors, positive squared amplitudes, flat edge weights and smooth square-root extension, then signed linear lift/mean correction and averaged quadratic cancellation. Reference covariance algebra alone is insufficient.
- [ ] **F57C-recursion:** actual n=1 and distinct n>=2 equations and moment repair, cancellation/absorption of nonflat leading origin E, finite-order estimates and smooth sum. Then resolve physical u/v/w, independently validate corrected Cartesian NS residual and prescribed-domain energy, and measure scale-recursion/material winding.

Global cone/lift/waves, completed flatness, corrected NS/energy and genuine temporal recursion gates remain false. The long-term goal remains active.

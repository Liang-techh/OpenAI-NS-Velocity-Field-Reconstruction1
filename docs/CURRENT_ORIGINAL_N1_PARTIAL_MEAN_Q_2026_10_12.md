# Interior own mean, radial recovery and reduced pressure-tail route

The actual finite-N correction now has a **free-endpoint own-mean function through ordinary Z3**, its log-radius derivative through Z3, and the recovered meridional Q increment through Z2. Directed whole-cell ranges execute on the two named micro cells. The original finite N, common velocity unit S=Pstar, phase and predecessor history are retained. Full reconstruction remains **ACTIVE / INCOMPLETE**.

Read this before [the finite-N velocity/endpoint checkpoint](CURRENT_ORIGINAL_N1_MODIFIED_LEADING_C3_2026_10_12.md). This replaces the endpoint-only mean limitation for the new micro-cell queries. It also changes the critical pressure-tail route: integration by parts removes the need for internal radial/phase velocity derivatives in the two integrated forcing moments. Finite endpoint q_R remains necessary.

## Moving endpoint, original history and physical units

At fixed Z, let y be the original log-radius offset, t its native chart coordinate, and x a distinct free native query. The signed correction is

```text
delta_m(x) = [exp(-(y(x)-y_left))*incoming_N_scaled_m
              + integral_left^x exp(-(y(x)-y(t)))*B(t)*J(t) dt]/N
delta_m_y  = delta_V-delta_m
delta_Q    = [2Z*delta_V-(1-delta)*Z*delta_m
              -(1-Z^2)*delta_m_Z]/(1-delta*Z^2).
```

J is the original dy/dt and appears once. N is the original selected DoubleDyadicInteger, divided once and never materialized. The incoming mean is retained even where local B vanishes. The integration coordinate and free endpoint are separate graph symbols; a chart-right history is not substituted for every interior query.

All four mean rows are ordinary Z derivatives. Q through Z2 consumes delta_m through Z3. The original kernel, endpoints, Jacobian and phase are Z independent. The exact delta recipe is `exp(-4*(exp(40)+11)-30)`.

The numerical envelope on a target cell starts with its original incoming N-scaled Z3 magnitude and adds B's whole-source-cell caps times positive original dy widths up to the target right endpoint. The rate-one kernel and incoming decay are at most one for preceding source points. This bounds every free query inside the target cell. The envelope then uses the actual selected log(N), not a materialized integer. It is an absolute magnitude enclosure; it is neither a signed point value nor a relative-width certificate.

The signed functions cover all 17 original charts. Numerical acceptance here covers only Z in [0,1/2], first_micro native [1/4,1/2] and second_micro native [5/4,3/2], through the admitted micro range API. Physical increments are S*delta_m and S*delta_Q. Baseline mean and correction must still be combined before interpreting them as the complete physical field.

## Reduced pressure-tail integrals

Use physical R, y=log(R), z=Z, a=1-delta, d=1-z^2 and L=1-delta*z^2. For the leading physical axial velocity U and own mean m=Mz0/R,

```text
m_y=U-m
L*q=2z*U-a*z*m-d*m_Z
S=(q+a*z*q_Z/2+q_y)/L + q*(q/2+q_y)
  + U*(-2z*q+d*q_Z-2z*q_y)/L
N1p=-S/2+2*q_R+R*q_RR.
```

Here I0(f)=integral f dR, I1(f)=integral R*f dR and B(f)=f(upper)-f(lower). Exact integration by parts gives

```text
I0(S)=B(R*q)/L + a*z*I0(q_Z)/(2L) + B(R*q^2)/2
      +d*I0(U*q_Z)/L -2z^2*B(R*U^2)/L^2
      -2delta*z^2*I0(U^2)/L^2 +2z*d*I0(U*U_Z)/L^2

I1(S)=B(R^2*q)/L -I1(q)/L +a*z*I1(q_Z)/(2L)
      +B(R^2*q^2)/2 -I1(q^2)/2 +d*I1(U*q_Z)/L
      -2z^2*B(R^2*U^2)/L^2 +2a*z^2*I1(U^2)/L^2
      +2z*d*I1(U*U_Z)/L^2

integral N1p dR  = -I0(S)/2 + B(q+R*q_R)
integral R*N1p dR= -I1(S)/2 + B(R^2*q_R).
```

Thus the ten scalar moment integrands need only U, U_Z, q and q_Z. Their ordinary Z1 rows need U through Z2, own mean through Z3 and q through Z2. Internal radial/phase derivatives are no longer prerequisites for these integrals. They remain relevant for pointwise stress, derivative joins and finite-endpoint q_R.

Every moment and boundary in these formulas is in physical R units. For rho=Lambda*R, dR=d(rho)/Lambda, weighted measure R*dR=d(rho)*rho/Lambda^2 and q_R=Lambda*q_rho. Neither radius factors nor common S may be restored twice. Upper boundary terms may vanish only after a proof for the actual completed leading support; the n1 velocity cutoff does not prove leading q or q_R vanish there.

## Evidence and callable entry

The independent exact fixture checks 11 rows: mean Z0..3, mean_y Z0..3 and Q Z0..2, with a nonzero incoming mean, oscillatory density and a moving endpoint. A separate symbolic arbitrary-smooth-mean identity verifies both reduced S formulas and the two finite viscous boundary identities. The read-only GPT-5.6 Luna/max reviewer independently confirmed the Jacobian/N multiplicity, ordinary Z grades, incoming history and the whole-cell envelope conditions. Source, delta, scope, packet mutation and foreign-packet checks preserve the narrow contract.

```python
owner = CurrentOriginalN1PartialMeanQ(accepted_modified_leading_C3)
packet = owner.range_cell(('0', '.5'), 'first_micro', (1, 4), (1, 2))
view = owner.report(packet)
```

Module: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q.py`.

Gate: `current_original_finite_N_partial_own_mean_Z3_mean_y_Z3_and_Q_increment_Z2_installed`.

Identity helper: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_radial_free_pressure_integrals.py`. Its physical formulas are verified; it does not itself supply completed numerical pressure tails.

## Ordered tasks for continuing agents

- [x] Keep the original selected finite N, ordinary Z convention, delta and common S unchanged.
- [x] Define the signed partial Duhamel mean with a separate free endpoint in every original chart.
- [x] Retain predecessor history and the genuine original native Jacobian once.
- [x] Export mean/mean_y through Z3 and the meridional Q increment through Z2.
- [x] Execute genuine whole-cell micro envelopes, distinct from outgoing endpoint histories.
- [x] Verify the reduced physical S/RS and finite viscous endpoint identities independently.
- [ ] Bind baseline `common_own_five_histories_axial5['m']` and `common_velocity_V_axial5` from the same live leading source, through ordinary Z3. These inputs already exist in the generic source packets; do not replace baseline history by zero.
- [ ] Combine baseline and finite-N correction in their common S=Pstar unit, then restore S once. Recover complete physical q through Z2 with the same delta/P0/source family.
- [ ] Define and integrate the four unweighted and six R-weighted nonradial moment integrands above. Use physical dR=R*J*dt and R*dR=R^2*J*dt; keep Z1 rows and true whole-cell source ranges.
- [ ] Extend numerical partial-mean bounds to every actual leading chart needed by the tail. Preserve original native partitions, positive widths, fixed-Z geometry and all quiet-chart incoming histories.
- [ ] Accumulate the ten moments from the actual original Rin to the genuine leading source-support endpoint. Use the unchanged background before r_minus and actual finite-N leading field after it.
- [ ] Supply lower finite-endpoint q/q_R from the unchanged original Rin source. Retain the upper terms until actual completed leading support and q=q_R=0 are bound; local loop flatness is not leading radial support.
- [ ] Apply the reduced physical pressure identities to obtain both complete fixed forcing integrals, then append them with m13 plus the unweighted tail and m15 minus the weighted tail.
- [ ] Bind actual positive I2/Astar and the normalized 2+3 disjoint repair bumps. Compute alpha=-Bz^-1*(m11,m14/Astar), beta=-Btheta^-1*(m12,m13/(2*Astar),-m15/Astar), with complete signed defects.
- [ ] Reconstruct the repaired n1 own M/V/P from repaired F1/Uz1 and verify all five terminal identities, exact own recovery and pressure compatibility.
- [ ] Finish original full-inner axis remainder/derivative joins and whole-domain n1 contracts, then genuine n2 forcing, finite-order remainder and smooth temporal summation.
- [ ] Complete admissible stress, oscillatory averaged stress cancellation, flat remainder and independent corrected NS residual diagnostics at their actual scope.

All full-inner/solved/repair/completed-leading-support/n2/temporal-recursion/global-stress/flat-remainder/full-corrected-NS gates remain false. A local pressure-integrand or moment-envelope interface must not be reported as a completed pressure-tail closure.

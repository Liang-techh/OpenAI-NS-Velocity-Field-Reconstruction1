# Original waiting: full-moment stress and physical zero remainder

The original O7 waiting region now has a callable source-bound stress and
absolute-pressure companion, followed by a physical divergence-form tensor.
Both focused checkers pass with current source hashes. The native waiting
velocity and pressure datum are retained; no new profile or pressure tail is
fitted. This advances leading-background stress recovery, not coefficient
recursion or a full corrected Navier–Stokes solution.

## Region and actual endpoint data

Use the original selected waiting length and phase:

\[
q=\log(R/R_{\rm tail})=\mathrm{wait}\,(\mathrm{phase}-1),
\qquad 0\le\mathrm{phase}\le1,
\qquad -\mathrm{wait}\le q\le0.
\]

Set \(a=\delta/2\), \(k=1-a\), \(p=1+\delta\),
\(b=(1-\delta)/2\), \(d=1-Z^2\), \(L=1-\delta Z^2\),
\(\epsilon=.001\delta\), \(K_0=1-\epsilon\),
\(Q_0=K_0^2-1=-2\epsilon+\epsilon^2\).
The actual source is \(B=E_{v0}\theta_{\rm base}e^{-(1/2+a)q}\),
\(R=R_{\rm tail}e^q\), \(U_\theta=BK_0\).

The endpoint defects \(A_{d0},E_{d0},P_{d0}\) come from the accepted collar at
offset zero, including its complete sigma/phi/Gamma future integrals to
infinity. The new implementation evaluates the collar only at zero. It does
not call any collar provider outside its admitted domain.

The accepted angular repair, positive half-normalized future energy and
absolute-pressure closure identify these endpoint moments with the original
waiting history. The selected pulse terminal meridional moments, the actual
`packet -> _packet -> flatten_mixed` zero-velocity route and the zero-density
FTC preserve \(M_z=M_{\theta z}=U_r=U_z=0\) throughout waiting.

## Backward transport and baseline cancellation

With \(A=1/k+A_d\), \(E=1/\delta+E_d\),
\(P=1/(2p)+P_d\), the original FTC equations are

\[
A_q=K_0-kA,\qquad E_q=\delta E-K_0^2,
\qquad P_q=pP-K_0^2/2.
\]

Their exact backward solutions are

\[
A_d=-\epsilon/k+(A_{d0}+\epsilon/k)e^{-kq},
\]
\[
E_d=Q_0/\delta+(E_{d0}-Q_0/\delta)e^{\delta q},\qquad
P_d=Q_0/(2p)+(P_{d0}-Q_0/(2p))e^{pq}.
\]

Unit baselines cancel analytically before interval enclosure. Define

\[
N_\theta=kA_{d0}+\epsilon-bZ\partial_ZA_{d0},
\]
\[
N_E=\delta ZE_{d0}-ZQ_0-d\partial_ZE_{d0}/2,
\quad
N_P=-2pZP_{d0}+ZQ_0+d\partial_ZP_{d0}.
\]

The normalized original stresses are

\[
C_\theta=\frac{N_\theta e^{-kq}}L
-2S(1+a)K_0e^{-q},\qquad
C_z=\frac{N_Ee^{\delta q}+N_Pe^{pq}}L.
\]

Here \(S=1/R_{\rm tail}>0\) is the exact source; `[0,S_cap]` only encloses it.
The numerical box does not itself prove a strict shear sign, so the new
waiting shear-sign/cone flag remains false. The exact strength margin
\(\kappa-2=\delta>0\) is available.

Mixed-three stress rows include the positive factor derivatives:
\(Q_\theta=\sqrt{R/2}B\) and \(Q_z=\sqrt{R/2}B^2\).
The angular inertial/shear radial rates are \(-1,-(1+a)\); the axial
energy/pressure rates are \(-1/2,+1/2\). Original forward moment rows remain
accessible beside the equivalent full-future representation.

The original absolute pressure in Rtail reference units is

\[
\frac{P_{\rm profile}}{P_\ast^2}
=-\mathsf C\left[\frac{K_0^2}{2p}e^{-pq}
+P_{d0}-\frac{Q_0}{2p}\right],\qquad
\mathsf C=\frac{E_{v0}^2\theta_{\rm base}^2}{P_\ast^2}.
\]

This retains the full future datum; it is equivalent to the original forward
pressure by the admitted endpoint closure and the same density FTC. Runtime
`Ev2` is an enclosure of the exact amplitude source, not a replacement value.

## Physical decomposition

The accepted general-K physical transfer is reused through a separate phase
adapter. It calls the original waiting chart at its original phase and uses
constant waiting K rows. No interval quotient is used to infer a new phase
from an independently enclosed negative offset. Existing admitted source
files and their receipts remain unchanged. The production waiting-radius AST
is also evaluated symbolically to prove
\(\log R(\mathrm{phase})-\log R(1)=\mathrm{wait}(\mathrm{phase}-1)=q\).

The physical map applies on \(\tau>0,r>0,|Z|<1\); source endpoints
\(Z=\pm1\) represent physical infinity limits. The physical stress prefactor
is \(\nu\lambda^{-2-\delta}\). Complete the
symmetric tensor with \(T_{r\theta}=T_\theta\), \(T_{rz}=T_z\),
\(T_{\theta\theta}=r\partial_zT_z\); other entries are zero. The completion
cancels radial tensor divergence.

The physical swirl is exactly

\[
u_\theta=\sqrt\nu\lambda^{-1-\delta}B K_0
=\nu^{1+a}E_{v0}\theta_{\rm base}K_0
(2R_{\rm tail})^{1/2+a}r^{-1-\delta}.
\]

It is independent of physical z and time. Hence
\(E_\theta=-\nu\partial_{zz}u_\theta=0\) exactly throughout waiting.
All regional remainder components vanish, while the stress and background
momentum residual generally do not:

\[
\mathcal R_B=-\nabla\cdot T_B+E_B,
\qquad E_B=0\quad\text{in waiting}.
\]

Angular inertial stress and axial energy stress have homogeneous divergences
that cancel before enclosure. The reduced divergence grids retain the
angular shear and axial pressure contributions, including their source
radius factor derivatives. Pressure axial/time derivatives retain their
full datum; only the proved pure-power velocity rows are set to exact zero.

## Evidence and reproduction

Run from the repository root:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage waitingphysical
```

This runs the two new producers and their checkers in dependency order.
`waitingstress` runs just the similarity-stress/pressure pair.

- Similarity checker: 80 finite signed stress rows, 60 pressure rows,
  96 exact zero meridional-history rows, 72 actual flat endpoint K derivative
  checks and 234 current input hashes.
- Actual waiting/collar formula ASTs are replayed on arbitrary full terminal
  axial functions: 20 mixed stress and 15 mixed pressure join identities.
  Sample interval overlap is not used as functional proof.
- Independent moderate transport fixture: full convergent future integrals,
  60 original stress derivatives, 45 pressure derivatives and 27 moment FTC
  checks. Its analytic future is a test model, not the actual Gamma source.
  Actual Gamma histories are consumed separately through accepted receipts.
- Physical checker: 124 finite signed rows, 72 exact zeros and 314 current
  input hashes. Native implicit-coordinate Cartesian differentiation at
  viscosities .01 and .7 checks all six decomposition components, including
  the full completed tensor. Errors are below 1.3e-61. Another 24 mixed-two
  reduced-divergence comparisons pass; both nonzero stresses are exercised.
  These moderate errors are not the project's full NS residual.

## Remaining work

The waiting/collar stress mixed-three and pressure mixed-four joins are
admitted. The steep-exit/waiting stress join is still false until the
preceding stress representation and its common primitive route are replayed.
The waiting cone, all-region stress cone, global flat remainder and volume
bounds, physical energy, finite-width upstream bridge feedback, coupled
coefficient recursion, oscillatory correction and full corrected NS residual
remain unfinished. A zero regional remainder does not close any global gate.

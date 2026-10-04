# Whole original waiting two-vector cone

The original waiting region now has a whole-domain admissible two-vector
stress-cone proof, using its actual full-future endpoint moments and original
source amplitude. This composes with the accepted whole collar proof at their
functional common join. No velocity, pressure, repair coefficient or previous
source receipt is changed.

## Source domain and mode argument

Let \(q=\log(R/R_{\rm tail})=\mathrm{wait}(\mathrm{phase}-1)\),
\(0\le\mathrm{phase}\le1\), \(-\mathrm{wait}\le q\le0\).
Use \(a=\delta/2\), \(k=1-a>0\), \(p=1+\delta\),
\(K_0=1-\epsilon>0\), and exact \(S=1/R_{\rm tail}>0\).
The entire source Z domain is [-1,1]. Finite physical points have |Z|<1;
source endpoints are physical infinity limits.

The accepted full moments give the coefficient functions at the collar
endpoint:

\[
I_0(Z)=\frac{kA_{d0}+\epsilon-bZ\partial_ZA_{d0}}{1-\delta Z^2},
\]
\[
J_E(Z)=\frac{\delta ZE_{d0}-ZQ_0-(1-Z^2)\partial_ZE_{d0}/2}
{1-\delta Z^2},
\quad
J_P(Z)=\frac{-2pZP_{d0}+ZQ_0+(1-Z^2)\partial_ZP_{d0}}
{1-\delta Z^2}.
\]

Here \(Q_0=-2\epsilon+\epsilon^2\). These are obtained from the SAME
complete sigma/phi/Gamma future moments, not a new terminal value selection.
The source-proven waiting stresses are

\[
c_\theta=I_0e^{-kq}-2S(1+a)K_0e^{-q},
\qquad c_z=J_Ee^{\delta q}+J_Pe^{pq}.
\]

For \(v=-q\in[0,\mathrm{wait}]\),

\[
c_\theta=e^{kv}\left[I_0-2S(1+a)K_0e^{av}\right].
\]

Thus the whole-region sufficient margin is

\[
\gamma_\theta=I_{0,\min}
-2S_{\rm cap}(1+a)K_{0,\max}e^{a\,\mathrm{wait}}>0,
\qquad c_\theta\ge\gamma_\theta.
\]

This explicitly controls the faster-growing negative shear. Endpoint
positivity alone would not suffice. The actual whole-Z source enclosure
gives \(\gamma_\theta/\epsilon\) a lower bound approximately 2.6788.

Both axial modes contract backwards, so

\[
|c_z|\le C_z=|J_E|_{\max}+|J_P|_{\max},
\qquad C_z/\epsilon\lesssim5.16953.
\]

The absolute-mode sum uses directed interval arithmetic. Numerical bounds
are used only for inequalities; no cap endpoint defines a field value.

## Strict source signs and directional inequality

The exact source S is positive because it is the exponential of minus the
original finite log radius. Therefore the waiting source shear is strictly
negative, although its `[0,S_cap]` enclosure contains zero. The new proof
does not infer this strict sign from an enclosure endpoint.

Constant K0 gives \(\kappa-2=\delta>0\). Since
\(Q_z=Q_\theta B\), \(S_z=0\), \(c_\theta>0\), and the angular shear is
negative, the two-vector dot-product sign is correct. The remaining
directional condition reduces to

\[
2c_\theta^2-\delta B^2c_z^2>0.
\]

The actual B increases backwards to its original waiting inlet:

\[
B_{\max}=E_{v0}\theta_{\rm base}e^{(1/2+a)\mathrm{wait}}
=E_{v0}\theta_T/K_0.
\]

The production log-factor function and the owned inlet cancellation are
executed symbolically before interval enclosure. They verify Bmax is the
original source B at q=-wait, with no materialized amplitude cap.
The whole-domain directed test is

\[
\log\delta+2\log B_{\max}
+2\log(C_z/\gamma_\theta)<\log2.
\]

It passes with the actual source logs. Fixed positive viscosity and physical
lambda factors preserve the two-vector cone. The completed tensor diagonal
is outside this two-vector condition.

## Accepted companions and joined scope

`CertifiedWaitingPhysical().waiting(Z,phase,...)` consumes the accepted
receipt, verifies all its input hashes and family/source identities, and
adds regional cone metadata to the unchanged physical source packet. It
does not project stress components or alter pressure/velocity.

The waiting/collar mixed-three stress and mixed-four pressure joins are
already proved by actual formula AST replay on arbitrary full terminal
functions. The new cone proof consumes those joins and the actual production
phase/radius identity. It composes the two regional cone proofs on

\[
R_t\le R<R_{\rm tail}e^3.
\]

Stress vanishes exactly at the Gamma join and in the accepted heat exterior.
Those zeros are separate from the strict inequalities on the nonzero tail.
The physical waiting remainder remains exactly zero as previously proved.

Run from the repository root:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage waitingcone
```

This runs only the new cone producer/checker against current admitted
waiting/physical/collar prerequisites. The checker recomputes the entire
source Z enclosure and positive margins, replays the exact mode and source
log identities, and covers every waiting phase through the functional
theorem. No additional phase sampling or repeated Cartesian fixture is
needed for this cone increment. A bounded read-only review accepted the
mode bounds, strict source sign, log factors and cone normalization.

## Remaining work

The steep-exit/waiting stress join and preceding steep/entry/angular/power/
flatten stress-cone regions remain open. This is a joined outer-tail proof,
not the whole outer domain or the global completed tensor. Independent
global flat remainder and volume bounds, required-domain energy, upstream
finite-width bridge feedback, coupled coefficient recursion, oscillatory
correction and the full corrected NS residual remain unfinished.

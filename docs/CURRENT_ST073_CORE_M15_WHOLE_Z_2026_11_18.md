# ST073 whole-axis core m15 and the signed six-source inlet

The first two actual core contributions and the sum of all six actual core and lower-collar contributions now have directed C0 and ordinary Z1 enclosures in 128 closed cells covering \(Z\in[-1,1]\). The six-term sum is the actual \(q=2\) inlet, at \(s=s_c/2\). Its unit is \(\epsilon^2P_\star^2\). The first core axial term is also available in its own archived unit \(\epsilon^2\).

## Core axial source and numerical errors

The previously bound axis source \(C(Z)=-2T_{-2}T_{-1-\delta}(4Z+j)\) has regular raw axial response \(U_z^{(1)}/\epsilon=C(Z)\rho/4\) for constant axis forcing. Its contribution to the core axial moment is \(4(4Z+j)C(Z)\) in \(\epsilon^2\) units, since \(2\int_0^4\rho/4\,d\rho=4\). The new calculation extracts the signed product from the saved collar axis model using its positive geometry weight and applies this core weight. It retains separate bounds for the nonaxis leading velocity, the nonaxis zeroth-word source, the first feedback word and the remaining words. The old recurrence and quadrature were left in their saved state.

## Shared physical pressure endpoint

The true part of the pressure density involving \(N_{1p}\) satisfies
\[
-\epsilon^2\rho N_{1p}
=\partial_\rho(-\epsilon^2\rho^2 b_\rho)+\epsilon^2\rho S_{\rm flow}/2,\qquad
Q=Q_{\rm axis}+\epsilon b .
\]
Analytic axis regularity removes the lower endpoint at \(\rho=0\). The radius-four endpoint is \(E_0=16\epsilon^2 b_\rho(4,Z)\). The lower-collar graph retains the identical inlet \(\Psi,\Psi_\rho,\operatorname{mean}\Psi\), since its ratio and \(\chi\) start at one. Both sides use the same own-mean identity to obtain the radial mean derivative. This proves the endpoint match for the function and its ordinary Z derivative.

The new graph stores a single shared \(E_0\) at local nodes [51,52], and the actual upper \(E_S=\epsilon^2\rho(S)^2 b_\rho(S,Z)\) at [61,64], with \(S=s_c/2\). Hence the core term \(-E_0\) and collar term \(E_0-E_S\) combine to \(-E_S\), before any interval subtraction of the same endpoint. The inertial and raw swirl integrals keep separate error bounds. The graph also exposes core axial [90,92] and core pressure [95,96] through the original radius-four moment source.

For the core endpoint enclosure itself, write \(\overline\Psi\) for the true radial average. Then \(\overline\Psi_\rho(4)-\Psi_\rho(4)/2=\int_0^1 t[\Psi_\rho(4t)-\Psi_\rho(4)]\,dt\). Its absolute bound and each ordinary Z derivative through order two use \( (2/3)\sup_{\rho\in[0,4]}|\partial_\rho^2\Psi| \) with the matching Z derivative. The accepted \(X_h\) embedding supplies this curvature of the nonlinear core. The computation combines the common \(\Psi_\rho(4)/2\) pieces algebraically and keeps the correction explicit. The split also keeps signed axis inertia, nonaxis inertia and the unsigned bound on the actual \(F_{\rm lead}F_1\) integral.

## Correlated pressure comparison

The first receipt gave containment for ten of twelve core comparisons with the old archives. At both off-axis inputs the actual pressure value interval intersected the archived interval and extended beyond its upper edge. The second receipt factors the common actual normalized \(\beta=2\) mass \(M_2\) before any pressure derivative products. This independently evaluates a single rational mass coefficient for each ordinary order. Late pressure masses, hydrodynamic and nonlinear inlet terms, the mean correction, inertia and swirl keep their own bounded contributions.

With this correlation the six first-row comparisons and the six second-row comparisons all lie inside the archived directed intervals. At \(Z=\pm0.3\), the tightened second-row value minus the archived finite decimal is in \([2.12145,2.13447]\times10^{-26}\) in \(\epsilon^2P_\star^2\) units. The narrower result comes from algebraic reuse of the actual mass. The archive comparison itself allows only numerical inclusion. The old binary rounding relation remains open.

The whole \(q=2\) six-source sum resolves value sign on 124 cells and derivative sign on 114. At \(Z=0.371\) its value and derivative are approximately \(-13.377894886158827\) and \(-268.98831174011117\) in the common pressure unit. The three preserved forward ledgers each contain 212 terms and retain all prior error budgets. The new inlet is an evaluation of the existing six actual sources.

## Verification and continuation

The [whole-axis core receipt](../experiments/root_st073/original_core_m15_whole_z_receipt.json.gz) exited zero, SHA-256 244e06ec4ddbd710f836a4d473f6a00e8b7584d6b3a3bde3e334be0842a1de30. The [shared mass receipt](../experiments/root_st073/original_core_pressure_mass_pair_receipt.json.gz) exited zero, SHA-256 4b432d1cd0c6d9372e13c3b800340896a2798f29f9893686cbbb590a34b45eff. Both used only frozen archived inputs and fresh source algebra.

The next source join is the genuine first-micro volume plus its upper \(RQ_y\) endpoint. Use the shared outgoing \(E_S\) at its lower boundary, and leave the separate active and compact errors in place until their own sources can be paired. The independent negative complete-pulse comparison and the once-only reserved own-swirl contribution remain separate from the \(q=2\) inlet.

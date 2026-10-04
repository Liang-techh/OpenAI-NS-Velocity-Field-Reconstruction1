# Original entry physical stress and retained remainder — 2026-10-04

CompliantSteepEntryPhysicalC2 maps the original sigmoid entry and its SAME
complete-future stress/absolute pressure through the unchanged general-K
physical transfer. The original field, coefficients, sigmoid and datum are
preserved. The entry regional decomposition and right physical connection
are now accepted; entry cone and global flatness remain open.

## Reproduce

    python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage steepentryphysical

The ordered stage runs the entry similarity producer/checker, then the new
steep_entry_physical_C2.py producer and steep_entry_physical_C2_check.py.
The callable is CompliantSteepEntryPhysicalC2.steep_in(Z,phase,log_tau,theta,viscosity).
Phase is the original entry t, so dq/dt=1, not a Ts-scaled power phase.

## Original physical map and tensor

The source radius is q=log(R/Rtail)=-wait-Ts-2+t. Dispatch keeps the original
steep_entry chart and consumes the accepted companion pressure through both
native mixed/Taylor packet routes. An ephemeral heat reference supplies the
actual entry K rows to the general mapper without evaluating a collar at a
negative collar offset.

The physical map retains

    lambda^2-lambda^(2*delta)*z^2/nu=tau,
    R=r^2/(2*nu*lambda^2),
    u_phys=sqrt(nu)*u_source, p_phys=nu*p_source.

The two original stress profiles carry nu*lambda^(-2-delta), with spatial
derivative and completion factors treated separately. The symmetric tensor
has Trtheta=Ttheta, Trz=Tz and Ttheta_theta=r*partial_z(Tz); other entries
are zero. Completing this diagonal cancels radial stress divergence.

The regional physical identity is

    partial_t u + (u dot grad)u + grad p - nu*Delta u = -div(T) + E,
    E=(0,-nu*partial_zz(utheta),0) in cylindrical components.

The remainder is generally nonzero. It is retained in the Cartesian
decomposition instead of being removed or relabelled globally flat.

## Exact divergence reduction

Let a=delta/2, k=1-a, p=1+delta, L=1-delta*Z^2 and d=1-Z^2. Original entry
K is independent of Z and has K_q/K=a-mu-(1-mu)*sigma. The full moments
satisfy A_q=K-kA, E_q=delta*E-K^2 and P_q=pP-K^2/2.

The homogeneous A and E terms cancel analytically, leaving

    Dtheta=-K_q/L+2*S*exp(-q)*(K_qq-p*K_q+a*(1+a)*K),
    Dz=(d*P_Z-2*Z*(p*P-K^2/2))/L.

With Pd=P-1/(2*p) and Qd=K^2-1, Dz is evaluated as

    Dz=(d*Pd_Z-2*p*Z*Pd+Z*Qd)/L.

Thus unit baselines cancel before enclosure while retaining the full
absolute-pressure history. Full-factor derivative rows use rates
-(a+1/2) for Dtheta and -p for Dz. The latter rows through order two are

    Dz, Z*(K^2)_q/L, Z*((K^2)_qq-p*(K^2)_q)/L.

Actual divergence-function AST replay proves these rows and mixed axial
derivatives. Physical mixed2 divergence then uses the original implicit map.
The same general-K viscosity operator gives

    Etheta coefficient = 2*F*K_q/L^3-4*Z^2*K_qq/L^2,
    F=1-2*(1-delta)*Z^2-delta^2*Z^4.

Its mixed2 derivatives need K through order four, exactly the admitted
original entry shape derivatives. Positive source factors remain logarithmic
and correlated; caps do not define replacement fields.

## Source and physical connections

Actual radius AST replay proves outer_angular offset0 equals entry t0,
and entry t1 equals steep_power phase0. Current native C4 left field and
pressure identities are consumed. This finishes the native left radius/
field/pressure connection; upstream angular stress is still unfinished.

At entry t1, the original flat sigmoid has sigma=1 and higher derivatives
zero. Actual theta and logarithmic rates equal the native power inlet.
Actual flatten_mixed source replay proves45 velocity mixed4 identities.
Both companion pressure packet routes consume the already admitted mixed4
pressure connection. The original entry/power K4 and stress mixed3 joins,
common radius/map/source factors and shared completed mapper give the right
physical stress/diagonal/divergence/remainder connection. Remainder operators
are the same source AST and require no derivatives beyond the admitted K4.

## Accepted evidence and remaining scope

Focused producer/checker PASS with338 current hashes:156 finite signed
physical rows and40 exact zeros. Stress is available through mixed3; diagonal,
divergence and axial-viscosity remainder through mixed2. The whole original
entry/Z source box is included, with finite logarithmic source bounds.

One independent fixture differentiates the implicit physical coordinates and
every entry of the native completed Cartesian tensor. It uses a moderate
nonconstant polynomial K and complete analytic future in the original entry
radius chart. At nu=.01 and .7 it verifies6 Cartesian decomposition equations,
two divergence zeros, two nonzero remainders,24 reduced-divergence mixed2
comparisons and12 remainder mixed2 comparisons, all within1e-36. Actual
sigmoid/Gamma history is consumed separately through the current source
receipts; this fixture is not a global actual-field NS accuracy result.

Entry cone, upstream angular stress, global flat/volume/energy estimates,
true n-dependent recursion, oscillatory correction and full corrected NS
validation remain unfinished. No cone coverage is extended by this stage.
Next: prove the entire original entry cone with variable kappa-2 in[2*mu,2],
actual source B logs and common moments, then compose the power tail. Continue
the upstream angular/power/flatten physical stress and cone afterward.

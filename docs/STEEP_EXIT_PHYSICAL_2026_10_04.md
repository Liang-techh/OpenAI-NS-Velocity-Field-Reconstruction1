# Original steep exit: physical stress and nonzero remainder — 2026-10-04

The steep-exit companion now transfers the original full-moment stresses into physical coordinates, completes the divergence-form tensor, and retains the regional axial-viscosity remainder. Its physical momentum decomposition and right waiting join are accepted. The steep-exit two-vector cone and independent global flatness/energy estimates remain unfinished.

## Actual field and physical map

CompliantSteepExitPhysicalC2.steep_out(Z,phase,log_tau,theta,viscosity) evaluates the original steep_exit chart at its original phase in [0,1]. The unchanged source relation is q=-wait-1+phase, relative to the same Rtail as the accepted waiting/collar/Gamma chain.

The physical map is

```
lambda^2-lambda^(2delta)*z^2/nu=tau
R=r^2/(2nu*lambda^2)
u_phys=sqrt(nu)*u_source
p_phys=nu*p_source
T_phys=nu*lambda^(-2-delta)*T_source
```

Ephemeral shape/stress/assembly/dispatch adapters reuse the unchanged accepted CompliantCollarPhysicalC2.collar general-K algebra. The original heat/collar source is never evaluated at negative q. The new pressure companion supplies the same full absolute datum to the original physical field assembly. Exact Ev0, source radius, original velocities and forward histories remain retained.

The symmetric completed tensor has Trtheta=Ttheta, Trz=Tz, Ttheta_theta=r*partial_z(Tz), with all other independent cylindrical entries zero. The completion cancels radial tensor divergence. Physical velocity divergence, radial momentum, radial remainder and axial remainder are exactly zero in this pure-swirl region.

## Divergence without large homogeneous cancellations

The source has K=K(q), independent of Z. Write L=1-delta*Z^2, d=1-Z^2, p=1+delta, a=delta/2 and Kj=partial_q^j K. Full angular and energy homogeneous moments cancel analytically before interval evaluation:

```
Dtheta=-K1/L
       +2S*exp(-q)*(K2-p*K1+a*(1+a)*K0_current)
Dz=(d*P_Z-2Z*(p*P-K0_current^2/2))/L
```

Here K0_current denotes the current zeroth K row, not the constant waiting value Ktail=1-epsilon. The external factors are Qtheta*sqrt(2/R)=B and Qz*sqrt(2/R)=B^2. Their q rates are -1/2-a and -p. The physical divergence scale is sqrt(nu)*lambda^(-3-delta)*Qalpha*sqrt(2/R).

The stationary pressure defect is separated too:

```
Q0=Ktail^2-1
Qplus(q)=K(q)^2-Ktail^2
Iplus(t)=integral_t^1 exp(-p*(v-t))*Qplus(v)dv
Isigned(t)=Q0*(1-exp(-p*(1-t)))/p+Iplus(t)
Np_tail=-2p*Z*Pd_tail+Z*Q0+d*partial_Z(Pd_tail)

full-factor axial row0=(Np_tail*exp(p*q)+Z*(Qplus-p*Iplus))/L
full-factor axial row1=Z*partial_q(K^2)/L
full-factor axial row2=Z*(partial_qq(K^2)-p*partial_q(K^2))/L
```

Both producer integrands are replayed from their actual ASTs. Their sigmoid callable, cumulative primitive, correlated distances, cell counts and endpoint overrides are the same. Exact integrand splitting and the stationary exponential integral prove the relation between full signed pressure and positive fluctuation integrals. Their numerical enclosure boxes need not be identical. The second derivative follows by applying the full-factor operator twice with Iplus_q=p*Iplus-Qplus.

## Regional remainder and waiting join

The complete angular momentum retains

```
Etheta=-nu*partial_zz(utheta)
Er=Ez=0
```

For source K independent of Z, its normalized coefficient is

```
2*(1-2*(1-delta)*Z^2-delta^2*Z^4)*K1/L^3
  -4*Z^2*K2/L^2
```

The exact physical factor is sqrt(nu)*lambda^(delta-3)*B. Mixed remainder derivatives through order 2 use K0 through K4 and the accepted physical operators. The generally nonzero interior remainder is preserved.

At phase 1, the original sigmoid is flat and K1 through K4 vanish. All six mixed2 remainder coefficients vanish exactly, and axial/time velocity derivatives become the waiting pure-power zero jets. The stress itself and its divergence generally remain nonzero there. Same source stress mixed3/pressure mixed4 joins, identical physical maps/factors and the flat-K operator identities establish the physical right join. No steep-power left join is claimed.

## Reproduction and evidence

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage steepexitphysical
```

The stage lists stress producer/checker followed by physical producer/checker. The focused physical producer/checker passed with 322 current checker input hashes:

- 147 finite signed physical rows and 49 exact structural/derivative zeros in the stored samples and whole-domain enclosures.
- Six exact mixed2 remainder zeros at the actual waiting endpoint.
- Direct original-sigma quadrature checks of the new positive pressure kernel at two interior phases and its exact endpoint zeros.
- Six independent native Cartesian momentum decomposition components, including the full tensor diagonal and moving cylindrical basis, at nu=.01 and .7.
- Both fixtures exercise nonzero axial viscosity.
- 24 native mixed2 divergence comparisons and 12 native mixed2 remainder comparisons against the production coefficient helpers.

The bounded Cartesian fixture uses a moderate nonconstant polynomial K with flat derivatives at its waiting endpoint, a constant waiting segment and complete convergent analytic future moments. It differentiates the implicit coordinates and every completed tensor component natively. Actual sigma/Gamma history is consumed separately through source-bound receipts. Its tolerance is 1e-36 and Cartesian decomposition errors are below 1.4e-60; these fixture errors are not a global NS residual certificate.

## Next dependency

Prove the original steep-exit whole-domain two-vector cone against the same full moments and actual B source logs. Use kappa-2=delta+2*k*(1-sigma), lying in [delta,2], without substituting waiting's constant delta. Establish theta sign, shear dot-product sign and the directional inequality on the entire phase/Z domain, then compose the accepted right join with waiting/collar.

Continue full stress/pressure through steep power and entry, then angular/power/flatten. Independently finish upstream finite-width bridge feedback, global flat/volume estimates and required-domain energy before admitting common leading data, genuine coefficient recursion, oscillatory correction or full corrected NS completion.

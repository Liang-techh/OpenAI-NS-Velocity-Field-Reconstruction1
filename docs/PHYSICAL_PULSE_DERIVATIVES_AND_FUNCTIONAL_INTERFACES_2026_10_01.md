# Actual physical pulse derivative map and functional interfaces

The same compliant leading pulse now has physical cylindrical component
derivatives in r and z through total order four, using the paper's actual
implicit lambda(t,z) and anisotropic velocity/pressure scales. A complete
six-chart bound ledger keeps huge positive radii/norms in logarithmic form.
The internal main/gap/end primitive formulas are now identified by exact
selected-source functional equations, rather than finite interval overlap.

This completes F38-D3c1g2a1b (physical bound ledger) and g2a2 (internal
functional chart identities). Quantitative flat velocity/moment comparison
bounds, two-sided O3/O5 high-order joining data and full post-pulse C4 are
still pending. No energy integral, stress cone, Cartesian vector residual,
temporal recursion or blow-up acceptance is claimed.

## Original physical coordinates

The supplied v2 paper Sections2.2/2.3 and Lemma2.1 give

    R=r²/(2lambda²), Z=z/lambda^(1-delta),
    lambda²-lambda^(2delta)z²=tau=1-physical_time>0,
    d=1-Z², L=1-delta*Z².

The physical components are lambda^beta times the leading profiles, with
beta=-1-delta for angular/axial velocity, beta=-1 for radial velocity, and
beta=-2-2delta for pressure. Existing normalized pulse jets differentiate
the actual profile components in y=log(R),Z. They are not already physical
r/z derivatives: lambda depends on physical z.

For every a+b<=4, the new exact map is

    dr^a dz^b(lambda^beta G)
      =lambda^(beta-a+b(delta-1))*(2/R)^(a/2)*H_ab,
    H_a0=prod_(j=0..a-1)(Dy-j/2)G,
    H_a,b+1=((beta+b(delta-1))*Z*H_ab+d*DZ(H_ab)-2Z*Dy(H_ab))/L.

The apparent -a from the radial lambda factor cancels when axial
differentiation acts on R^(-a/2). Retaining an extra -a in H would be wrong.
Every operator is a finite exact rational combination of the existing
15 true mixed derivatives. Ordinary Taylor coefficients are not treated
as derivative values.

## Uniform bound ledger and large scales

The six whole-Z source boxes cover entrance, main, exit, both gap charts
and end. Each actual mixed derivative operator is enclosed over these
entire boxes. Physical scales use the unchanged

    Rref=110*(Cstar*Pstar)^10,
    Rp=Rref*exp(yd+1+Tw),
    t_local=log(R/Rp).

Absolute radial origins, inverse-mu decay and finite offsets remain
separate source terms. No exp(logCstar), absolute Rv, or super-small end
exponential is materialized. A logarithmic upper bound represents a finite
formal positive norm; it is not a numerical value of the norm.

For angular/axial profiles the combined radial prefactor after a r
derivatives is

    Pstar*(2/Rp)^(a/2)*exp(-(.5+mu+a/2)*t_local).

For radial profiles it is

    Pstar*2^((a-1)/2)*Rp^((1-a)/2)*exp(-(mu+a/2)*t_local).

For pressure it is Pstar²*(2/Rp)^(a/2)*exp(-a*t_local/2).
The factored derivative jets already include the profile product rules
for sqrt(R/2) and inlet swirl, so those derivatives are not added again.
All rates are nonpositive. Source chart lower times are respectively
0, .02/mu,10/mu,11/mu,12/mu and13/mu-4. Inverse-mu and finite pieces are
kept separate when evaluating the log bounds.

Since lambda>=sqrt(tau), every negative physical exponent supplies an
upper factor tau^(gamma/2), gamma=beta-a+b(delta-1). Receipts give all
derivative bounds at log(tau)=-1,-10,-100; the formula supports any finite
log(tau). For finite physical z, Z lies in(-1,1); closed-Z interval boxes
also bound the limiting endpoints. These are pointwise component derivative
suprema, not kinetic energy or volume L2 norms. Cartesian vector derivatives
also need differentiation of the cylindrical basis.

## Exact functional chart composition

Let lambda_i=.5-i*mu. The full main integral is

    I_i=int_0^(11/mu) exp(lambda_i*v)*gp(mu*v)dv.

The actual selected end equations imply

    sum_j Aij*cj=-exp(-13lambda_i/mu)*(mi0+ap*I_i),
    Aij=exp(lambda_i*center_j)*Wji.

The certificate translates the PRODUCTION two-row inverse, common
log-end scale and incoming row normalization, then proves both equations
through the exact divided difference A2=A1+mu*D. Numerical caps are not
substituted into these identities. This yields the main/gap function
equality at xi11 for arbitrary smooth axial histories, including axial
derivatives through five.

The actual selected energy equation uses the necessary Rp/Rv factor:

    ap²*Kpulse+mu*exp(-26)*end_energy
      =(1-exp(-26))/4-mu*e0+mu*exp(-26)*ev.

It turns the forward main energy into the identical backward gap energy.
The substitutions xi=13+mu*s, D=-mu*s identify the two gap moment logs,
pressure time/decay, angular history and energy expressions. At s=-4,
both end bumps have full future supports, so their weights equal the gap
weights. Disjoint supports make quadratic cross terms exactly zero.
At s=0 the future linear supports are empty while ev remains positive.

At the pulse entrance, the original inlet values, angular history and
P0+Mp pressure are identified with the zero-length main integral. Repeated
y derivative agreement inside O4 follows from the same primitive ODEs,
common boundary functions and flat forcing derivatives. Independent
two-sided high-order O3/O5 field providers and quantitative velocity/moment
flat-comparison ledgers are still required for full pulse acceptance.

## Evidence and reproduction

Under experiments/root_st073, with prefix lei_ren_part1_paper_:

- compliant_pulse_physical_bounds.py/.json and _check.py/.json;
- compliant_pulse_interface_certificate.py/.json;
- compliant_reconstruction.py --stage physicaljets;
- compliant_reconstruction.py --stage interfacejets.

There are62 ordered modules. Sixteen symbolic chain-rule identities,
180 independently differentiated physical implicit-coordinate values,
180 independent radius/time scale identities and1080 actual full-chart
derivative/time-sector log bounds pass. The physical fixture solves the
implicit lambda equation directly in physical z with nonzero delta;
it does not reuse the production derivative operator.

Thirty-two source-bound functional interface identities pass, including
production inverse and log normalization bindings, true divided-difference
weights, selected energy/units, coordinate changes and full compact future
supports. A read-only GPT-5.6 Luna/max worker independently derived the
interface equations and confirmed the physical derivative recurrence.

Earlier source modules and receipts remain unchanged. New receipts bind
the entire existing source graph. Full pulse/outer C4, physical energy,
stress/flat remainder and actual n-dependent coefficient recursion remain
false. The next action is quantitative flat velocity/moment bounds and
two-sided high-order joins, followed by post-pulse derivative propagation.

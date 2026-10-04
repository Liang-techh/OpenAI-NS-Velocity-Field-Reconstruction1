# Actual collar heat-part cone and flat join direction - 2026-10-04

`collar_heat_cone` proves the admissible two-vector cone on the entire actual
source domain `1 <= offset < 3`, `-1 <= Z <= 1`. This is the heat part of the
existing collar, with its original cutoff, moments and positive source
scales. The sigmoid transition `0 <= offset <= 1` remains open, so the
whole-collar and global cone flags remain false.

## Exact source and normalization

The production sigmoid is identically one on this domain. With `y=3-offset`,
the SAME admitted source is

```
phi=exp(-4/y^2)
xi=2*(1-Z^2)*S*exp(-offset), exact S=1/Rtail>0
K=(1-eps*phi)*H_delta(xi)
eps=.001*delta>0, a=delta/2, k=1-a, b=(1-delta)/2
L=1-delta*Z^2
```

The full Gamma source has exactly zero stress. Subtract it from the same
complete angular, energy and pressure integrals, whose differences are
supported only before offset3. With `F=integral_offset^3 exp(k*(q-offset))
phi(q)/phi(offset)*H(xi(q)) dq`, the actual normalized angular stress is

```
ctheta/(eps*phi)=(H+k*F-b*Z*F_Z)/L
                +2*S*exp(-offset)*[(8/y^3+1+a)*H-H_t]
```

No field parameter is selected from an interval cap. The AST bridge binds
the production sigma/phi/K composition, canonical full Gamma source,
epsilon relation and admitted common-moment/Gamma closure. The actual B
amplitude logs come from the same selected inlet/pressure source.

## Whole-domain source bounds

The full positive Gamma expectation gives
`hmin=1-2*a*(1+a)*Scap <= H <= 1`,
`abs(Z*H_Z)<=4*a*(1+a)*Scap`, and
`0<=H_t<=2*a*(1+a)*Scap`, where H_t is the offset derivative,
`H_t=-xi*H_xi`, rather than the derivative in y=3-offset.
The sharp cubic bound `abs(Z)*(1-Z^2)<=1/2` also gives
`(1-Z^2)*abs(H_Z)<=2*a*(1+a)*Scap`, and hence
`(1-Z^2)*abs(Q_Z)<=8*eps*phi*a*(1+a)*Scap` for the squared-density
difference Q. Both canonical chain rules are verified symbolically.
The actual source verifies
`k*hmin-4*b*a*(1+a)*Scap>0` and
`(1+a)*hmin-2*a*(1+a)*Scap>0`. Therefore
`ctheta>=eps*phi*hmin>0` throughout the open domain.

For `0<=u<y`, the exact nonnegative identity

```
4*((y-u)^(-2)-y^(-2))-8*u/y^3
  =4*u^2*(3*y-2*u)/(y^3*(y-u)^2)
```

implies `integral_0^y phi(offset+u)/phi(offset) du <= y^3/8`.
The SAME squared-density difference is
`(-2*eps*phi+eps^2*phi^2)*H^2`. Regrouping its energy/pressure contribution
gives `abs(cz)<=eps*phi*y^3*Cz`, where
`Cz=(1+2*delta+4*a*(1+a)*Scap)/(4*(1-delta))`.

Since `Sz=0` and the already accepted shear margin is positive, the two
remaining cone inequalities are `ctheta>0` and
`2*ctheta^2-(kappa-2)*B^2*cz^2>0`. Here `K_y>0`, so
`0<kappa-2<=delta`. The actual source-log receipt proves

```
log(delta)+2*log(B)+2*log(8*Cz/hmin) < log(2)
```

uniformly on the full source domain. This bounds the normalized directional
margin; the unnormalized stress tends to zero at the join.

## Endpoint and remaining work

The uniform bound `abs(Tz/Ttheta)<=B*Cz*y^3/hmin` proves that the direction
tends to positive e_theta. The sharper lower angular term also gives
`abs(Tz/Ttheta)<=Bmax*Cz*y^6/(16*S*exp(-3)*hmin)` for fixed source parameters.
This matches the paper's exponents in (5.23)-(5.24): angular stress is
`phi*y^(-3)*btheta`, axial stress is `phi*y^3*bz`, and their ratio is O(y^6).
At offset3 the admitted common moments and flat source jets give exact zero
stress. Strict cone inequalities are imposed only before that endpoint.

This cone concerns `(Ttheta,Tz)` against `(Stheta,Sz)`. The completed diagonal
`T_theta_theta=r*partial_z(Tz)` does not enter it; no cone property for the
entire completed tensor is claimed. Physical momentum acceptance remains
separate, as do global flatness, energy, recursion and oscillatory correction.

Run the producer and focused source/proof checker with stage `collarheatcone`.
Next close the actual sigmoid interval `[0,1]` using the original full K and
common future moments; preserve this full-domain heat proof at their shared
offset1 boundary. Do not replace it with a finite-grid sign check or shrink
the original collar to a convenient width.

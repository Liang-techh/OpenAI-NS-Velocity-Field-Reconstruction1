# Actual R110 long angular reshape and inherited moments

The accepted ACTUAL R110 inlet now feeds the original long angular reshape
(4.38)/(9.30). The new provider encloses log(Utheta/Pstar), normalized true
Utheta derivatives, actual Uz, pressure datum and all five cumulative moment
histories through axial order5, with mean-recovered radial velocity through
axial4. It covers the entire original phase[0,1], including the actual Rsh
exit. It does not yet implement the following reference/axial continuation.

The source length remains the previously selected T=400*Abar, where Abar is
the admitted A_upper for the selected physical norm family. No shorter
length, alternative inlet, fitted B, finite core or comparison history is
substituted. These are interval enclosures of original analytic profiles,
not newly recomputed point coefficients or full mixed spatial derivatives.

## Exact log field and inherited normalization

Let y=log(R/110), s=y/T, u=Utheta, u1=actual Utheta110, v1=actual Uz110,
Qz=1+Z^2 and B=log(Cstar*u1*Qz). The source is

```
log(u/Pstar) = y/10-logCstar-logPstar-log(Qz)+(1-sigma(s))*B,
V(R,Z) = v1(Z),
a = 4/5 + 2*B*sigma'(s)/T.
```

B is taken directly from the actual110 log jet, after the exact logCstar
cancellation in the previous stage. Its coefficients through2 are
intersected with the accepted SAME-source theorem ||B||C2<=2*Abar. This
theorem does not replace the actual source or assert higher derivative
bounds; coefficients3..5 retain the analytic inlet bounds. The previously
admitted sigma derivative bound8 gives a in[.7,.9] and u_y/u in[.05,.15].

The current-amplitude normalization avoids forming impossible positive or
tiny source numbers:

```
Mtheta   = sqrt(2)*R^(3/2)*u*H
Mtheta_z = sqrt(2)*R^(3/2)*u*K
Mz       = R*m
Mztheta  = R*A - R*u^2*b/2
Mp       = u^2*p/2
P        = Pstar^2*original_P0_normalized + Mp.
```

H,K,m,A,b,p are normalized shapes. They are distinct from their true
derivatives divided by the current basepoint prefactor. Both are returned;
the latter include the exact Bell factors of the log-u jet.

From the previous stage's F0-normalized R110 shapes H110,K110,m110,A110,
B110,C110 and phi110, the normalized starting histories are

```
H0=H110/(2*phi110), K0=K110/(2*phi110),
b0=B110/phi110^2, p0=C110/phi110^2,
m0=m110, A0=A110.
```

This p0 is the incoming normalized pressure MOMENT. The original axis
pressure datum is stored separately and never reset. At y=0 the new
provider retains every normalized incoming interval exactly.

## Full, untruncated radial integrations

For the angular, pressure and swirl-energy moments respectively, let
(k,multiplicity)=(1.6,1),(.2,2),(1.2,2). Their endpoint-normalized kernels are

```
I(k,multiplicity;y,Z) = integral_0^y
 exp[-k*t + multiplicity*B(Z)*(sigma(y/T)-sigma((y-t)/T))]dt.
```

The corresponding inherited-decay factor is

```
D(k,multiplicity;y,Z) = exp[-k*y + multiplicity*B(Z)*sigma(y/T)].
```

The actual histories satisfy

```
H = H0*Dtheta + Itheta
K = K0*Dtheta + v1*Itheta
b = b0*Dswirl + Iswirl
p = p0*Dpressure + Ipressure
m = exp(-y)*m0 + (1-exp(-y))*v1
A = exp(-y)*A0 + (1-exp(-y))*v1^2.
```

Every incoming term remains in its exact source expression. Normalization
does not erase it even when its numerical enclosure containszero.

The original log-u slope yields positive kernel rate ranges[1.55,1.65],
[.1,.3] and[1.1,1.3]. Zeroth-order bounds use the corresponding exact
positive decay integrals. For axial derivatives,
|sigma((y-t)/T)-sigma(y/T)|<=8*t/T. The Taylor Bell coefficients are positive
polynomial majorants in t; their integrals are bounded by p!/rate^(p+1), or
the smaller finite-y mass y^(p+1)/(p+1). The full infinite integral is an
UPPER BOUND only: it does not define or truncate the original finite kernel.
Both signs of B are supported.

Inherited moment Taylor coefficients are convolved with the decay's Bell
factors before combining their norms with the exact source decay in logs.
Only then may a sufficiently small product be enclosed by exp(-1000).
This cap never defines the original amplitude. No enormous radius, exp(A),
exp(logCstar), selected F0 or physical u is materialized.

## Verification and reproduction

Run `lei_ren_part1_paper_compliant_reconstruction.py --stage reshapeprofiles`
from experiments/root_st073, with accepted prerequisites. The complete
ordered pipeline now has88 modules. The source and checker are
`lei_ren_part1_paper_compliant_long_reshape_profiles.py` and
`lei_ren_part1_paper_compliant_long_reshape_profiles_check.py`, each with a
separate source-bound JSON receipt.

Checks pass426 actual factored exponential-cap inequalities,385 profile/log
bounds and1254 normalized moment/kernel bounds. Exact checks preserve the
selected A/T, original pressure datum, all R110 normalized histories and
constant axial velocity. Independent full quadrature with both signs of B
passes144 kernel axial derivatives,144 nonzero-history decay derivatives
and96 normalized amplitude derivatives. Seven symbolic identities verify
the five PHYSICAL primitive RHSs and reference endpoint/radius compatibility.
Finite fixtures do not admit the actual source.

A read-only GPT-5.6 Luna/max audit verified the source interpolation,
physical moment factors, inherited decays and kernel derivative bounds.
Its request for direct proof bindings is addressed: B_C2_upper=2*Abar,
actual continuously inherited histories, certified sigma derivative8,
selected shear interval and recomputed log-u/kernel rates are checked
explicitly by both constructor and checker.

At s=1 the exact log field becomes
log(u/Pstar)=T/10-logCstar-logPstar-log(1+Z^2), equal to the reference at the
same formal Rsh=110exp(T). This log identity is not a claim that all radial
derivative interfaces or five reference moments are already certified.

## Next dependencies and open scope

Continue actual Rsh moments through the reference branch to Rz=e^-8*Rref,
where Rref=110*(Cstar*Pstar)^10. Use formal log-radius differences factored
BEFORE numerical enclosure and exact offsets -8,-7,-6,-5. Do not recover
these offsets by subtracting huge rounded absolute logs. Restore the
actual axial mismatch using the original sigma on[Rz,eRz], then continue
toRh. Derive the actual five defects for the existing repair on[Rm,2Rm].
Matching prescribed Utheta/Uz alone does not match their accumulated
moments, pressure or recovered Ur.

Bridge/switch/reshape phase-radial mixed4, full high-order core/annular joins,
complete physical field, local/terminal-energy domains, admissible stress
lift, independent flat remainder, actual n-dependent recursion and
oscillatory correction remain incomplete. The original unlocalized source
has infinite whole-space physical kinetic energy. This stage does not
change that source or insert a cutoff. Full-field/cone/energy/temporal gates
remainfalse.

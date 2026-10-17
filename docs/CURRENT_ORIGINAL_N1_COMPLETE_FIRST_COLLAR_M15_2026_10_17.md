# Complete actual first-collar m15 increment

Tasks14/15 and the full long-term goal remain **ACTIVE / INCOMPLETE**. This stage adds the actual cut-velocity and swirl-coupling contributions, with computed source errors, to the accepted leading-N1p contribution. The resulting complete selected first-collar m15 increment and its ordinary Z1 derivative are strictly negative at the saved Z interval near0.371. It is an increment over rho4 to the selected inner endpoint, not the core integral constant or full C_inner, and not an all-Z terminal identity.

## Source and physical radial measure

The canonical own-cutoff definition gives

```
S = sc/4, rho(s)=4*exp(hb*s), d-rho=hb*rho*ds
q=s/S, eta(s)=1-sigma(4q-2)
F1_cut=eta*F1_raw, Uz1_cut=eta*Uz1_raw
D15_rho=2*epsilon*Ulead*eta*Uz1_raw
         -epsilon^2*rho*(2*Flead*eta*F1_raw+N1p_physical).
```

N1p is uncut and in physical units. The preceding accepted pressure packet already encloses its integral with the weight -epsilon^2*rho. This stage adds exactly the other two integrals over0<=s<=S. It does not add the leading-N1p contribution twice, reset the core inlet, or multiply pressure/own mean by eta. Original leading-shear sigma(s) and the independent own-velocity eta(s) remain distinct. P1/m13 is not m15.

The exact C_inner source handles remain53161..53164. Raw first-chart F1 handles are43324..43327 and Uz1 handles43328..43331. The full source relation is

```
C_inner,15 = M15_core(4,Z) + DeltaM15_first_collar.
```

This stage evaluates only DeltaM15_first_collar. M15_core is not appended here, and the old six signed core/lower ledger items remain unchanged.

## Actual saved jets and F0 connection

The accepted selected-word packet supplies the signed core-entry finite prefix W0+W1 and its ordinary Z0/Z1 derivatives. Its inherited source intervals and core quadrature errors stay in every coefficient interval. Finite words are used as approximations to the actual all-word raw field, with local and omitted-word radii kept separately.

The accepted inner_operator archive supplies actual leading F inlet jets in samples.first_inlet.ordinary_leading_rho_rows.F[0]. Geometry is phase0, rho4. The checked receipt certifies the nonconstant-amplitude ordinary Z derivative convention. Its four source log bases are identical to the selected/coupled raw packet's bases. The saved F derivative coefficients already include the F0 derivative connection; it is not added again. The actual leading U inlet jets come from the checked leading_mean_epsilon_factor archive.

These are enclosures of the same original nonlinear fields, not signed values chosen from norm caps. The F inlet coefficient begins0.2829799868805425 times the exact positive F0_anchor factor. The source formula is F_lead=F0*Phi_actual; the raw F1 physical scaling is already present in the selected prefix. Multiplication adds the leading F0 sector once and retains all epsilon factors.

No Bessel/core owner is rerun. A bounded worker initially overlooked the stored leading-F rows; inspection of the actual archive found them and the worker verified the correction and receipt.

## Integral approximation and separate errors

The original cutoff has exact unweighted mass

```
integral_0^S eta(s) ds = (5/8)*S.
```

Let F_A,U_A be actual leading inlet enclosures and f_A,u_A the raw inlet finite-prefix enclosures. The signed inlet-product integral approximations are

```
I_U_model = 5*epsilon*hb*S*u_A*U_A
I_F_model = -20*epsilon^2*hb*S*f_A*F_A.
```

These models use actual saved signed coefficients, not midpoint/radius values. The following errors enclose their difference from the true integrals.

For each raw field and ordinary Z order separately,

```
E_raw = (sc/8)*abs(prefix_d_x_field_A)
        + local_prefix_radius(q=1)
        + uniform_raw_allword_tail.
```

The homogeneous displacement bound is sc/8. It is distinct from the looser source-kernel bound dx=sc/4 used in the saved local-radii formulas. Those positive formulas are monotone in q, so their q1 values are uniform on the entire selected window. One all-word tail at the evaluation point is added per field; the finite prefix is never declared the exact inlet function. Ordinary Z1 uses its own saved radius, not a formal derivative of an upper bound.

The leading-field variations are controlled on a common complex Z tube while rho/s remain real. The accepted inner analytic-domain archive provides the uniform complex amplitude cap Fmax and the actual Phi bound:

```
Fmax = original_axis_amplitude.complex_modulus_upper
     = about8.362428613635886e-817812180069139571
Phi_bound = about189253.6086780747
B_F = Fmax*Phi_bound.
```

Fmax is a bound on the physical F0 throughout the tube, not the center F0_anchor. The accepted source proves logCstar>Lambda*Gbar+2*logLambda+1000, hence log|F0|<=-2*logLambda-1000. The native ODE, positive core-Phi floor and real0<=sigma<=1 give

```
|DeltaF_lead| <= hb*Cd*Ac*B_F
|DeltaU_lead| <= epsilon*hb*Cd*actual_Psi_rho_bound
rhoB-4=hb*Cd, rhoB^2-16=2*hb*Cw.
```

Here Ac is the admitted complex core Phi_rho/Phi modulus bound. The actual native Psi_rho bound includes its ratio/chi cover. Cauchy with gap eta/16 bounds the ordinary Z1 of each whole difference; only Z is complexified. No complex-rho continuation through a flat switch endpoint is used.

For inlet approximations A,B with symmetric variations EA,EB, the code evaluates

```
E_AB = abs(A)*EB+abs(B)*EA+EA*EB.
```

Its ordinary Z1 jet includes all six product-rule error terms. For rho^k, k1 or2, the additional geometry radius is bounded using the exact hb*Cd or2*hb*Cw difference, retaining its extra hb factor. Seven fresh symbolic checks cover these product identities, both inlet-mass coefficients, both radial geometry identities and positivity by real FTC.

## Actual numerical result and scope

New1800-bit directed arithmetic computes both signed models and all their error radii. Before merging with pressure, every one of44 axial signed/error monomials and19 swirl signed/error monomials is separately bounded in normalized absolute value by2^-1024. The corresponding normalized absolute upper bounds are44*2^-1024 (about2.448e-307) and19*2^-1024 (about1.057e-307). These are bounds only, never signed component values.

After adding the accepted same-source leading-N1p packet, the complete collar increment is enclosed as follows. Values below are rounded decimal displays of the directed intervals; exact endpoints are serialized mpf tuples.

| ordinary Z order | DeltaM15 / (hb*Pstar^2*epsilon^2*S) | relative interval width upper |
|---|---:|---:|
|0|approximately-26.75578977231765467503449|7.788e-29|
|1|approximately-537.9766234802223394446576|2.093e-29|

The dominant actual coefficient uncertainty remains in the interval. Each of74 remaining merged signed/error monomials is bounded by2^-1024 after identical source factors cancel before log boxing. This is a certified actual local integral enclosure, not a claim that its tiny omitted components vanish.

Production exits0 in 2.976 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_complete_first_collar_m15.json` is 11,109,344 bytes, SHA256 `333d9afeb4d16542b922aecb41101875eac12c3256f980bce1a13f1574a167bd`, with 1987 input pins and 330 implementation pins. Focused Python compilation, seven fresh symbolic checks, saved-only report inspection and scoped raw-byte index checks pass. A fixed Luna/max worker independently reviewed the new physical factors, error algebra, tube bounds and scope; root integrated its required gates. Tools do not expose independent runtime model metadata, so none is claimed.

Every inherited source/error packet remains equal at JSON-value level except stage/hash/contract/check/time metadata. Rc178/quiet179/raw213/reduced183, the six signed core/lower terms, signed slots19/20/41/42, all finite-N volume caps and target/quiet consumed-once state are unchanged. No old error cap is deleted, no ancestor owner/source fixture/recurrence/norm builder/old quadrature runs, and no new item is appended to the old six-term ledger.

## Next executable work

Recover the same-source actual raw M15_core(4,Z) value with a complete signed approximation/error budget, then add this accepted complete collar increment exactly once to form signed C_inner. Distinguish forward normalized core fixtures from genuine all-word raw core data, and retain physical epsilon/F0/Pstar factors and ordinary Z1. Extend the axial coverage required for functional claims before any terminal five-moment conclusion.

The accepted core endpoint ports53125/53126 substitute cumulative core functions52601/52602 at rho4; they have no installed numerical values. The full endpoint ports53161/53162 already contain the collar and must not be added to a separate collar increment. The old core_m15_numeric archive (SHA256 ba2462a9a31259ff20ef82701335f6bfc44698d0bb5af3aebfe2d54a94005f8f) contains Z=-.3,0,.3 only, so its nearest .3 value cannot substitute for .371. Its useful exact scale identity is Mcore=epsilon^2*[2*integral Ulead*u+Pstar^2*(integral p-4*p(4))], with u=Uz1/epsilon and p=P1/(epsilon*Pstar^2). The missing next operation is a new source adapter at .371 with the complete C0/C1 error, not a replay of the old core owner.

Independent reserved own_swirl/negative complete pulse comparison, outer p2/phase/H0/signed C/I1/I2/H, all-five/common-N functional closure, exact heat/energy, admissible stress/flat remainder, true temporal n>=2 recursion, mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The full-field1e-3 threshold is not applied to this leading source component. The long-term goal remains active.

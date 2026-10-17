# Actual signed A-minus-A0 and a correlated q2 inlet

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. At the same saved Z interval near .371, the actual signed remainder epsilon^2*integral_0^(sc/2)(A-A0)*d-sigma and its ordinary Z derivative are now enclosed with strict negative sign. A new typed q2 inlet provider exposes one shared source packet to new first_micro and preRp consumers. The old numerical owners and six-row archives are unchanged; downstream first_micro increments, full forward Rp and all-Z matching remain open.

## Actual source, scale and positive measure

For S=sc/2, rho(s)=4*exp(hb*s), L=1-delta*Z^2, pcore=Psi_core,rho(rho(s),Z), pA=Psi_core,rho(4,Z), and R=Phi_actual/Phi_core,

```
A(s) = 2*rho(s)^2*Z*(1-hb)*R(s)*pcore(s)/L
A0   = 32*Z*pA/L
J6   = epsilon^2*integral_0^S (A(s)-A0)*sigma_prime(s) ds.
```

The accepted local primitive nodes are [81,94], principal A0*sigma nodes [83,96], with physical epsilon squared outside exactly once. Radial-order-one source checks exclude radial-order-zero Psi. The current native first-window identity covers q in [0,2], uses the actual pointwise quotient and has no extra exponential after the evolved field. The exact S is the selected sc/2, not a midpoint of a cover.

Both sigma(S)>0 and d-sigma>0 are certified from the original monotone flat source, sigma(0)=0, and sigma is ordinary-Z independent. Sigma(S) stays an exact positive source factor with its original log interval. It is neither the own velocity cutoff sigma(4q-2) nor a selected vanishing cap.

## Signed model and complete variation radius

The exact algebra is

```
A-A0 = -hb*A0 + (1-hb)*(2Z/L)*B
B = (rho^2-16)*R*pcore + 16*(R-1)*pcore + 16*(pcore-pA).
```

Thus the signed integral model is -hb*epsilon^2*A0*sigma(S). The actual variation is retained as an error radius, with the same admitted complex Z tube:

```
|B|/hb <= 4*Cw*B1 + 32*S*rho_max*Ac*sigma(S)*B1 + 16*Cd*B2
Cd = (rho(S)-4)/hb; Cw = (rho(S)^2-16)/(2hb)
B1 >= |Psi_core,rho|; B2 >= |Psi_core,rhorho|; Ac >= |Phi_core,rho/Phi_core|.
```

The ratio modulus is at most 2; its deviation has the same-source bound |R-1|/hb<=2*S*rho_max*Ac*sigma(S). The code checks the actual nonzero core denominator, complex L modulus lower bound, common outer/Psi/Qpsi tubes, accepted ratio exponent margin and actual core inlet rows. Multiplication by 2*|Z|max/Lmin gives Rbound. The integrated error is at most hb*epsilon^2*Rbound*sigma(S), because the positive d-sigma mass is sigma(S). The ordinary Z1 error bounds the entire error function by Cauchy with gap eta/16; a cap is not differentiated as though it were the function.

The directed normalized enclosures have the following rounded displays:

| ordinary Z order | J6 / (hb*epsilon^2*Pstar^2*sigma(sc/2)) | relative interval width upper |
|---|---:|---:|
| 0 | approximately -28.381823699387874303 | 1.644e-29 |
| 1 | approximately -90.110255333306329434 | 2.717e-29 |

Every error remains a radius around a signed model; no old norm bound, midpoint or overlap supplies a signed value. Exact directed endpoints and formal physical scale factors are stored in the report. No hyper-exponential scale is materialized.

## Six correlated source entries

The new inlet has these six function entries:

1. core axial source;
2. core pressure and finite-boundary source, including its unsigned swirl error;
3. cutoff axial source;
4. direct leading pressure + cutoff swirl - principal - signed remainder;
5. principal epsilon^2*A0*sigma(S);
6. actual signed J6.

The fourth identity is bound to the accepted same-leading pressure partition. Principal and remainder are shared function leaves, including their signed error functions. Adding the six source weights removes both principal and remainder before boxing, leaving exactly core + cutoff velocity + direct leading pressure[0,S]. The aggregate is the accepted actual q2 enclosure from the signed-lower report. A negative error radius is never substituted for a signed term; independently summing all component error radii would lose this correlation and is explicitly disallowed.

Each leaf records whether it uses the original five-factor or extended six-factor scale basis. No factor list is silently truncated. The old sixth name represented a width-error enclosure; the new sixth name is `first_lower_signed_A_minus_A0`, so it cannot masquerade as the old unsigned error slot. Old six numerical rows and all their finite-N errors remain archived and active.

## Executable shared inlet

`ActualQ2Inlet(report).first_micro_inlet('.371')` and `.preRp_inner_lower('0.371')` return caller copies of the same checked actual q2 source packet through two explicit new consumer interfaces. They check the hard-pinned accepted parent, the current provider source hash, actual source family/sc tuple, aggregate equality and the shared source-weight reduction. The exact endpoint tuple is S=sc/2; the q baseline sc/4 is stored separately, and both are checked against the accepted parent. A digest binds the complete per-order physical source/error leaves, which accompany every query. Foreign axial queries are rejected; returned caller copies cannot mutate the provider. The query keeps its input decimal spelling while the admitted rational is exactly 371/1000.

This is a new correlated format, not a drop-in invocation of the old six-row numerical owner. The old first_micro and preRp classes are not modified or executed. A subsequent actual first_micro increment must consume this inlet and preserve its common source/error leaves; this stage supplies the cumulative initial moment at sc/2 and does not invent the downstream increment.

## Run and next work

Production exits 0 in 2.698 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_correlated_q2_inlet.json` is 14,463,158 bytes, SHA256 `420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13`, with 1991 input and 338 implementation pins. Focused compilation, six fresh exact identities, directed variation arithmetic, source-weight cancellation, parent immutability, shared-provider queries, foreign-Z/caller-mutation rejection, saved-only inspection and scoped raw index checks pass. Fixed Luna/max workers checked source ports and independently reviewed the signed primitive and correlated inlet. Tool responses expose no independent runtime model metadata; none is claimed.

No accepted owner, fixture, recurrence, norm builder or old quadrature is replayed. Inherited packets remain equal at JSON-value level except stage/hash/contract/check/time metadata. Original six terms, Rc178/quiet179/raw213/reduced183, signed slots 19/20/41/42, all finite-N errors and target/quiet consumed-once contracts remain unchanged. The finite old prefix is not declared equal to the infinite raw function.

Next consume this actual q2 source in a new first_micro continuation. Bind the lower R*Q_y endpoint to the canonical/native sc/2 endpoint using the accepted zero-correction inlet theorem and actual V_y source; keep the nonzero accumulated n1 pressure/mean. Contract only proved same-function internal endpoints, then evaluate signed increment sources with ordinary Z1 and full finite-N error provenance. Do not rerun the old numeric inlet or append a seventh inner term.

Full first_micro propagation, full forward Rp, all-five/common-N functional terminal closure, exact heat/finite energy, admissible stress/flat remainder, true temporal n>=2 recursion/summation, realizable mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. Goal ACTIVE.

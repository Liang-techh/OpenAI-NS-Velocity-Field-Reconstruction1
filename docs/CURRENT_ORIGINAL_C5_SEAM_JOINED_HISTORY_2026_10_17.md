# Same-source C5 phase seam and joined signed mean history

Tasks 14/15 remain **IN PROGRESS**. The long-term goal remains **ACTIVE / INCOMPLETE**. This stage certifies the O2-buffer upper/O3-transition lower phase-primitive seam as a function of independent issued phase, removes its paired internal endpoint from the signed transport expression, and applies the known pressure functional to the joined inlet. Signed numerical inlet values, full n1 function equality and global terminal closure remain open.

## Actual seam, with distinct endpoint substitutions

The actual O2 buffer coordinate232895 is evaluated at source node4=11. The actual O3 transition coordinate241536 is evaluated at source node0=0. The latter is the lower endpoint; node1=1 is the distinct outgoing transition endpoint. Physical log-radius offsets6+exp(40)+11 and17+exp(40)+0 agree exactly, both physical Jacobians equal1, and the same Rm factor gives the same physical radius.

The transition source starts from the actual buffer outgoing E. Its extra exponent -x/2-mu*J(x) vanishes at x=0 because J has equal integral endpoints. The flat sigma branch is zero there; E_y=-E/2 and V=V_y=0 agree at the seam. The original Duhamel construction starts with the buffer outgoing histories and zero-width lower integrals. The shared P0, leading fields, physical radius and parameters therefore feed the same generic recovery and the same prescribed root inputs.

The new source-only certificate binds those facts to saved graph nodes and pinned definitions. Endpoint normalization explicitly checks20 history-field rows,8 radial rows and the radius, then20 generic theta/axial/p2 rows and32 root rows. Buffer endpoint aliases are validated against their actual wrappers; moving-field aliases require native11, while global inputs retain their actual source-node identity. It proves equality of the loop as a function of arbitrary independent phase, with consistent renaming of the angle and its bound integral variable, rather than equating unrelated actual-phase node IDs or matching a sample. The phase inverse uses its original monotone branch. The universal anchored primitive C(y,Z,phi)=integral_0^phi B(y,Z,theta)dtheta consequently agrees at the seam. Ordinary Z0..Z3 follow from the same function identity and the accepted C5 source regularity.

Both sides also use the identical unwrapped global phase N*(log(R/Ra)-phase_origin), with actual N node18, Ra node65604 and phase_origin node24954. Equal endpoint radius gives equal unwrapped phase, and hence the same periodic issued phase after any wrap. The result identifies two endpoint functions; it does not set either endpoint individually to zero.

## Joined signed transport and pressure application

Let M1 and M2 be the actual rate-one memory factors of the buffer and transition, I1/I2 their signed slow integrals from the accepted IBP expressions, and H0=N*delta-m at the original buffer inlet. C_in is the initial buffer primitive and C_out the outgoing transition primitive. The two nested updates reduce exactly to

```
H_joined = M2*M1*H0
         + (C_out-M2*M1*C_in-M2*I1-I2)/N.
```

The removed pair is M2*(C_buffer_at11-C_transition_at0)/N and is zero by the same-source seam certificate. The initial signed history, both outer endpoints and both signed slow integrals remain. All four ordinary Z rows are carried. Separate moving-query and actual transition-endpoint1 expressions are retained. The new graph references the immutable accepted197-node IBP graph through an explicit nested graph namespace; it does not overwrite that graph or confuse its node IDs with the original C5 archive.

New universal algebra checks compare the nested and contracted formulas through ordinary Z0..Z3, retain the incoming response and both slow integrals, and recover the exact nonzero signed difference when the two seam values are unequal. No positive bound is subtracted and no interval center is chosen.

The new endpoint H0..H3 are then supplied to the previously accepted pure O3 pressure kernel as a new source application. With physical Q0=-4*mstar*exp(-12), independent of Z, this preserves

```
qhat = -[Z*(1-delta)*H+(1-Z^2)*H_Z]/(1-delta*Z^2),
Delta m15 = R0^2*((Pstar/N)*A+(Pstar/N)^2*B),
A = Z*(1-delta)*qhat_Z*expm1(x)/(4*(1-delta*Z^2))-Q0*qhat*x/2,
B = -qhat^2*x/4.
```

Ordinary Z1 consumes H_Z3 and remains explicit. H already includes its N scaling; no extra N is added. The background cross term and both nonlinear factors are retained. The pressure kernel's old controls, coefficient interval evaluation and source builder are not rerun. This graph has 210 local nodes. Its functions still require signed phase/slow-integral evaluation; it does not contain evaluated numerical H or a numerical pressure increment.

## Exact weighted complex Green kernel, with actual bounds still open

A separate new kernel handles v=(F1/x^2,Uz1/x^2,K1/x^2,P1/x^2,F1x/x,Uz1x/x). It reads the pinned raw matrix AST without executing that builder. All16 off-kinematic entries satisfy x^(p_j-q_i)B_ij=rho*A_ij with even nonnegative residual x degree, p=(2,2,2,2,1,1), q=(1,0,0) on source rows4/5/6. This removes axis divisions by analytic continuation, conditional on the actual leading fields being holomorphic in rho.

The six exact source Green integral constants are1/8,1/4,1/8,1/2,1/4,1/2. For the operator acting on rho*A*v along s=rho*u, its six kernels are u(1-u)/4, u*log(1/u)/4, -u(1-u)/4, u/2, u^2/2 and u/2. Their absolute values are at most1/2 on0<=u<=1. These kernels already include the three kinematic rows; no additional I+J factor is applied.

For genuine source/operator sup bounds, set Msrc=max(||g4/x||/2,||g5||/4,||g6||/2) and C_a=max row-sum sup|A^a_ij|. A total ordinary-Z Cauchy gap sigma is divided into sigma/k at word order k. Ordered complex radial paths and k!>= (k/e)^k give

```
word_k <= Msrc*(R/2)^k*(C0+k*C1/sigma)^k/k!,
beta = e*R/2*(C0+C1/sigma),
||v|| <= Msrc/(1-beta) <= 2*Msrc   if beta<=1/2.
```

New symbolic controls verify the six integrals, transformed kernels, logarithmic-kernel maximum and factorial integral identity. At R=9/2, r=4, the weighted degree-N Taylor factor is9*(8/9)^(N+1). A raw base-field degree16 prefix gives only a weighted degree15 prefix because W1..W4=rho*v. Raw derivative rows remain x times a rho-holomorphic row; their tail therefore also carries sqrt(r) when converting a weighted remainder back to W5/W6.

This is an exact conditional kernel and source-operator parity result. No actual complex source sup norms, outward solution norm or radial remainder is installed. An earlier approximate candidate mixed raw and weighted norms and used rounding without certification; it was rejected and replaced, and none of its numerical estimates appear as accepted bounds. Actual coefficient holomorphy and nonvanishing denominators on the common tube, Xh-to-complex-sup embedding and directed source/A0/A1/Msrc/beta arithmetic remain required. The kernel permits exact rational or outward interval arithmetic for word evaluation; it does not certify a supplied numerical input merely by accepting it.

## Actual execution and preserved state

Production exits0 in 3.584 s. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C5_seam_joined_history.json` is 4,805,892 bytes, SHA256 `35fdc828c112e5fb09bf28a3f709b1d8a9817b54cd3a99fc11f76cfd92251d16`. It pins 1975 input files and 291 implementation files. Saved-only inspection, focused compilation and raw staged-byte audit support this publication. The production does not execute accepted owners, constructors, old fixtures, numerical core recurrences or quadratures.

All Rc178 and quiet179 views and inherited raw213/reduced183 packets remain equal to the accepted parent, including every active error, six signed core/lower contributions, signed leading leaves19/20/41/42, remaining O3 lower endpoint, all volume caps, external viscous contraction and pressure throughRv. The target and quiet pair were already consumed once. P1/m13 remains distinct from m15. Original uncommitted files are preserved.

Bounded workers use the fixed luna_max_worker route, initially fork_turns:none, requested GPT-5.6 Luna/max. Tool responses provide task/nickname rather than independent runtime-model metadata. Root owns integration, the contraction formula, ordinary-Z controls, pressure application and publication. No Astra child is requested.

## Remaining work

Implement directed signed evaluation of the remaining phase primitives and slow integrals, starting from the retained prior history. The new seam identity may simplify that evaluator; it does not justify deleting any existing unsigned volume or finite-N error cap. Use function definitions, numerical approximations and certified error bounds as separate objects.

For the core, a complex-radial regularized n1 source/operator/Green and solution bound at some R>4 is still required before applying a radial Taylor remainder. The actual smooth collar is not that analytic continuation. A Volterra word tail on real rho<=4.1 cannot replace this bound. Raw W5/W6 are odd in x=sqrt(rho); use divided-by-x or axis-weighted coordinates for rho holomorphy.

Full core response equality and the q-dependent collar bridge still precede any removal of the six signed core/lower terms. Global functional all-five/common-N closure, exact heat/energy, admissible stress/flat remainder, true temporal n>=2 recursion and smooth summation, realizable mean/oscillatory correction and independent corrected Cartesian NS remain incomplete. The final corrected-NS1e-3 threshold is not applied here.

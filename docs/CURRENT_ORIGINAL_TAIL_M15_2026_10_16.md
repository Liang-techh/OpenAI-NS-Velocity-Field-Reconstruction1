# Original compact/quiet seams and zero relative viscous tail

Tasks14/15 remain **IN PROGRESS**; the full goal remains **ACTIVE / INCOMPLETE**. The two actual O3/compact/quiet source seams and a separately checked compact mass-row identity now prove D(2Rc)=D(Rp)=0, where D=R*(Qy_modified-Qy_original_leading). Absolute mean and pressure histories are retained. The previous long-to-Rc signed viscous chain plus compact sector210 therefore becomes D(long native0). Raw213 and the independent negative complete-pulse m15 target remain verbatim. A separate correlated view has189 terms, replacing190, with188 other entries unchanged. This is not global five-moment or numeric Rp closure.

## Source coordinates and physical recovery

| Seam | Left source coordinate/query node/value | Right source coordinate/query node/value |
| --- | --- | --- |
|O3 power upper Rc -> compact lower Rc|250215 /48 /2|260181 /1 /1|
|compact upper2Rc -> quiet lower2Rc|260181 /48 /2|263344 /0 /0|

Graph handle48 is the exact number2. With physical y=log(R/110), the O3 source has J=1, compact uses x=R/Rc and J=1/x, and quiet uses the physical log-offset s with J=1. The original O3 global phase is preserved where present. Compact repair and quiet power flow retain the same original N and parameter context; they are not assigned a new fast phase.

The saved numeric quiet frontend uses a fraction q, not source s: s=L*q, L=Tw-2-log(2). Its q=1 packet is Rp, corresponding to source endpoint node263346, not s=1. The actual quiet frontend geometry AST/hash is bound to this conversion and its Jacobian L. Source quiet J=1 and numeric fractional Jacobian L are kept distinct; physical dy is applied once. Rp radius is node263349.

Thirty recovered physical trace rows bind U3,Uy2,m3,Q3,Qy2,RQy2 at both seams. An additional modified m_Z3 input is retained to recover Q_Z2; Qy_Z1 uses U/m through Z2 and Uy through Z1. Original leading V3,m3,Vy2 are independently bound at each seam, including quiet's separate leading_flow and inlet.leading_fields. No modified mean is reused under an original-leading name.

Compact inlet and terminal packets do not contain a physical subgroup. Their U/Uy/m/Q/Qy endpoint evaluations are explicitly local recipes on actual full-chart function nodes and actual native substitution scopes. These recipes are never registered as archived graph IDs. The actual compact RQy endpoint nodes are263155..263158 and262914..262917. O3 and quiet actual physical endpoint wrappers remain explicitly bound. The shared positive pressure denominator L=1-delta*Z^2, original Pstar factor and physical R factor are retained once.

## Independent compact M-row identity

The new bounded algebra checker treats the actual source jets of A, incoming relative mass, N and h as opaque functions. It admits only exact rational graph operations, together with the exact constant identity exp(-log(2))=1/2. For ordinary Z orders0 through3 it independently verifies all12 identities:

```
N * incoming_relative_m = A * target_M,
residual_M = target_M + h0 + h1,
terminal_relative_m = (A * residual_M)/(2*N).
```

All products use the ordinary derivative binomial coefficients, including1,3,3,1 at Z3. The actual mass row of B is[1,1,0,0,0] and its quadratic term is zero. The five value leaves belong to the same actual unique Banach system259187. Their ordinary derivatives match the actual shared-Jacobian recurrence; the M-row known right-hand side is exactly the corresponding target_M derivative. The accepted same-branch invertibility/IFT premise makes residual_M zero as a function. No finite iterate, midpoint, interval endpoint, control-zero choice or numeric cap supplies its value.

The actual terminal certificate node is tied to the actual relative mass endpoint rows262678..262683, terminal RHS263189..263194, the same source system and the pinned original full-bump endpoint theorem/AST. This gives zero relative mass through the required ordinary derivatives. Original compact beta support is flat at both ends; modified and leading U/Uy coincide there. Applying the original pressure formula gives

```
Qy=[2Z Uy-(1-delta)Z(U-m)-(1-Z^2)(Uz-mz)]/(1-delta Z^2),
DeltaU=DeltaUy=Delta_m=Delta_mZ=Delta_mZZ=0,
DeltaQy=DeltaQyZ=0; D(2Rc)=D_Z(2Rc)=0.
```

These are source identities on the admitted axial branch, not conclusions from the three saved Z queries. The mass-branch parameter derivatives are not the still-missing temporal n>=2 coefficient recursion.

## Quiet transport, signed chain and preserved errors

Actual quiet source profiles give U=Uy=0 while the original finite-mu swirl profile E remains in the source flow. Both modified and leading physical means satisfy m_y=-m and R_y=R. Hence D_y=0 and the proved zero at2Rc propagates to Rp, including ordinary Z1. The absolute quiet mean remains the original K*Z/R, with its source K retained; Q=-K/R and Qy=K/R are not set to zero.

The source chain is

```
[Dlong0-DRc] + [DcompactRc-Dcompact2Rc]
 = Dlong0-Dcompact2Rc = Dlong0,
DRc=DcompactRc; Dcompact2Rc=Dquiet0=DRp=0.
```

Only the old external replacement and compact finite-N viscous sector210 are replaced. The active error budget is now a directed positive union of the two saved own-long0 C0/Z1 radius-weighted pressure bounds. The full old202 and compact210 aggregates remain in raw/provenance data and are not active budget components. Their removal is justified by exact signed source identities; no positive cap is subtracted. The numeric adapter retains a safe full210 fallback when the zero certificate is absent.

Exact full dictionaries verify compact lower incoming207, complete upper211, quiet lower incoming211 and complete upper212 packets, their zero lower increments, and the same selected N. Reserved212 is retained once. All188 other sectors remain verbatim, including compact mean-square207, endpoint208, volume209, quiet mean-square211 and all prior numerical/width/volume errors. The six signed inner/lower contributions and independent target/reference*A*Rv^2 remain unchanged. P1/m13 and m15 remain distinct data.

## Actual execution and limits

New production exits0 in12.02s; saved inspection exits0 without constructing any owner or invoking the fixture. Report `lei_ren_part1_paper_compliant_current_original_tail_m15_numeric.json` has11,238,040 bytes, SHA256 `a394d3b538e420936e382b36b63dc7405b7a098b57a557c668deba7a2cc60c9a`, 1953 input pins and183 implementation pins. One prepublication attempt rejected a full-function versus endpoint leading namespace mismatch; the actual scopes were then aligned before any report was accepted. No accepted ancestor constructor, producer, checker or calculation was replayed.

The independent N=17 synthetic local-germ fixture passes all12 controls. Maximum supporting difference is0.0; its nonzero incoming/leading means, varying compact Jacobian, physical R factor, ordinary Z1, genuine leading subtraction and separate relative-terminal-zero case are retained. The toy cutoff checks the derivative orders used in the fixture; original beta flatness comes from the separate source theorem. This is not an original Poisson inverse evaluation or an NS residual.

The largest remaining cap is the surviving long0 replacement: normalized log about5.13045578980e+408906090034569676 at Z=0 and3.46460448428e+408906090034569677 at Z=+/-.3. These remain enormous upper bounds, not signed observations or useful terminal accuracy. The terminal relative viscous function is zero, but the full cumulative m15 difference is not thereby zero. Original signed inlet values, numeric Rp/global all-five/common-N closure, exact heat/finite energy, admissible stress/flat remainder, true temporal n>=2/summation, realizable mean/oscillatory correction and corrected Cartesian NS remain incomplete.

## Next executable task

Connect the surviving own-long0 physical trace to the preceding post_power upper endpoint using actual source germs and a same-source numerical packet binding. Prioritize a batch of still-uncontracted active finite-N viscous sectors27/50/68/87 where their interfaces can be proved. Macro leading viscous23/24 and first-switch leading/source-integral43..47 are different terms and cannot be deleted as finite-N endpoint errors. Preserve raw213,189-view/188 other sectors and the own-long0 budget until the upstream signed endpoint identity and new surviving inlet bound are both established. Do not rerun accepted ancestors or reinterpret cap overlap as equality.

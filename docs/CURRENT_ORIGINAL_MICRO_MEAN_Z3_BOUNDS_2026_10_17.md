# Actual micro C/Cy and complete physical mean Z3 bounds

Tasks 14/15 and the long-term goal remain **ACTIVE / INCOMPLETE**. The
source-owned C and Cy primitives and both moving convolutions now have
directed ordinary Z0..Z3 norm caps on the original micro charts for
Z in [37/100,93/250]. The complete physical mean correction and its total
physical-y derivative are each bounded by 2^-1024 through ordinary Z3.
This is a local correction-error statement, not an NS residual threshold.

## Functions, derivatives and physical scale

The accepted graph functions remain authoritative. Write B0 for the source
field named B_over_Pstar, with physical B=Pstar*B0. The new implementation
consumes the saved ordinary derivative caps of B0 and fixed-phase slow B0_y,
including all E product-rule cross terms. With target phase phi in [0,1],

```
C(y,Z,phi)  = integral_0^phi B0(y,Z,theta) dtheta
Cy(y,Z,phi) = integral_0^phi B0_slow_y(y,Z,theta) dtheta.
```

Thus the phase integration length is at most 1. Each physical Pstar*C_Zk and
Pstar*Cy_Zk cap uses the corresponding saved B0 or slow B0_y cap normalized by frozen
E0, multiplied by the already accepted uniform upper bound on Pstar*E.
E0 denotes the actual positive local E, frozen only for differentiation;
it is not a chosen value. These phase primitives commute with ordinary Z
derivatives because the actual phase and integration endpoints are
Z-independent. There are no hidden factorials.

C(0)=C(1)=0 is an accepted identity of parameter functions. Differentiating
it at fixed phase gives the corresponding Cy endpoint identity. Periodic
extension across integer fast-phase cuts therefore retains the source
transport identity; no huge-N phase is selected or materialized.

## Joined moving history

The first original chart is s in [sc/2,1], the second is s in [1,2], and
y=hb*s=log(R/Ra). The source-owned left extension has zero correction.
For either primitive f=C or Cy, the first convolution is

```
Jf_first(s) = integral_(sc/2)^s hb*exp(-hb*(s-t))*f_first(t) dt.
Jf_second_joint(s) = exp(-hb*(s-1))*Jf_first(1)
                    + integral_1^s hb*exp(-hb*(s-t))*f_second(t) dt.
```

The second formula retains the first-chart incoming history. Multiplying
the memory exponential by the first kernel gives exactly the full joined
kernel; it is not an interval-overlap argument. The positive mass is
1-exp(-hb*length), at most hb for either native width<=1 and at most 2hb for
the joined history. The hb Jacobian is applied once. Source-proved first
Cin=0 is used only at the original left collar. Neither second history is
reset, and C at the outer endpoint is retained.

The physical correction is bounded as the complete expression

```
delta_m = Pstar/N^2 * (C - JC_joint - JCy_joint).
delta_m_total_y + delta_m = Pstar*B0/N.
```

The second line is the accepted source transport identity, with
C_total_y=Cy+N*B0, rather than a slow-y substitution. It provides the total
y derivative bound without requiring a second slow-y derivative. Both
identities hold through ordinary Z3 using the same common N and Pstar.

## Certified bounds and retained limitations

Write Csel for the selected positive logC. Both charts satisfy
log||delta_m_Zk|| < -4500*Csel and
log||delta_m_total_y_Zk|| < -600*Csel for k=0..3. Each is also strictly
below log(2^-1024). Directed records, including physical scales, are in the
report; the following are rounded upper-endpoint displays divided by Csel.

| ordinary Z order | first mean log-cap/Csel | second mean log-cap/Csel | first total-y mean log-cap/Csel | second total-y mean log-cap/Csel |
|---|---:|---:|---:|---:|
| 0 | -10702.2466839868503580702 | -10702.2466839868503580702 | -5877.62334199342517903509 | -5977.62334199342517903509 |
| 1 | -8682.24668398685035807018 | -8682.24668398685035807018 | -4669.62334199342517903509 | -4968.62334199342517903509 |
| 2 | -6662.24668398685035807018 | -6662.24668398685035807018 | -2649.62334199342517903509 | -3147.62334199342517903509 |
| 3 | -4642.24668398685035807018 | -4642.24668398685035807018 | -629.623341993425179035091 | -1326.62334199342517903509 |

The identical first/second mean bound at Z3 comes from the retained first
Cy history. It does not assert equality of those signed function values.
The report also separates the known C part and the retained Cy part.
The norms give symmetric error enclosures around zero; the actual mean
function is not set equal to zero. No cap, midpoint or interval overlap is
used as a signed value.

All old 212 forward terms, six signed inner/lower sources, endpoint
identities and error packets remain unchanged. The ancestor source packet's
historical source-only flags remain unchanged; this new packet supplies the
norm bounds that were previously missing. Actual signed mean values,
finite-N velocity derivative errors, QZZ, outer RQy, independent terminal
matching and global five-moment closure remain open. Heat/energy,
stress/flat remainder, true n-dependent recursion/summation, realizable
stress corrections and independent corrected Cartesian NS validation also
remain incomplete.

## Execution and next work

New production exits 0 in 3.114 seconds. Report `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_micro_mean_Z3_bounds.json`
is 17,772,325 bytes, SHA256 `557195be5893c17e9398371017b3b4b7a52b6a0eae473b6befc19cafe84b8aef`, with 2003 input and
371 implementation pins. Saved-only inspection,
same-source AST guards, joined-kernel and physical-transport controls,
directed norm checks and independent bounded Luna/max math review pass.
Separate runtime model metadata is unavailable. No accepted producer,
owner, norm computation, recurrence, quadrature or seam checker was replayed.

Implementation: [mean Z3 provider](../src/openai_ns_reconstruction/original_micro_mean_Z3_bounds.py),
[runner](../scripts/run_original_micro_mean_Z3_bounds.py).

```
python scripts/run_original_micro_mean_Z3_bounds.py --inspect-report
```

NEXT bound delta_U=Pstar*B0/N through ordinary Z2 and apply the exact linear
Q recovery to delta_U and delta_m_Z0..Z3, producing a delta_Q_Z0..Z2 error
provider. For the outer RQy boundary, retain delta_m_y=delta_U-delta_m and
delta_U_y=Pstar*(slow_B_y/N+B0_phi); the N phase chain is already canceled
in the B0_phi term. Then bound the six C1 I1 defects and outer boundary,
retaining any non-small bounds, before signed micro integration and terminal
reserved-pulse functional matching. Goal ACTIVE.

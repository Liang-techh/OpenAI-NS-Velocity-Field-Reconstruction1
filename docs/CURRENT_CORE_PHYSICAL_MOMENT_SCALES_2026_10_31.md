# Saved core inlet moments with physical scale factors — 2026-10-31

The checked core-interior archive already contains the six inlet moments, their ordinary derivatives and radial-tail bounds. This component reads those saved rows directly: it does not rerun their integral, coefficient recurrence, fixed-point solver, provider, or earlier checker. It adds a physical-scale and ordinary-Z1 interface while keeping the saved finite terms and tails separate.

At fixed normalized radius \(\rho=4\), with \(F_0>0\), \(\epsilon>0\), and \(\epsilon\) independent of \(Z\), the six physical radial integrals are

\[
(I_V,I_{V^2},I_{RF},I_{RFV},I_{RF^2},I_{F^2})=
(4\epsilon M,4\epsilon A,8\epsilon^2F_0H,8\epsilon^2F_0K,
16\epsilon^2F_0^2B,4\epsilon F_0^2C).
\]

For any entry written as \(c\epsilon^pF_0^kJ\), the Z1 coefficient is \(c(J_Z+kLJ)\), where the saved relative derivative \(L=F_{0,Z}/F_0\) is a directed interval. The separate inherited absolute-tail bound is \(c(e_{J_Z}+k\sup|L|e_J)\); none is cancelled. The output gives the finite coefficient enclosure, the tail cap and their widened enclosure, as well as the exact symbolic powers \((p,k)\). Neither \(F_0\) nor \(\epsilon\) is assigned a numerical representative. If either radius or \(\epsilon\) depends on \(Z\), this fixed-radius/constant-parameter formula does not apply without the missing chain terms. P1 and m15 remain different data; no Pstar factor is introduced into these normalized moments.

The actual input was the authenticated prior core-interior archive's `inlet` view, covering \(Z\in[-1,1]\) at \(\rho=4\), and its independent `fresh_inlet` view at \(Z=0.537\). For both, all six scale rows and their ordinary Z1 rows returned with the inherited tails. The archived Q4 Phi component declarations, the core-interior archive and the fixed-point check all carry the same three family/source/datum hashes. This is a match of saved source records, not a proof of capture-safe Q4 substitution or full archived-constructor lineage. The sample's `I_RF2` Z1 coefficient enclosure is strictly negative, conditional on those saved enclosures and the positive symbolic factors; it is not a global stress or terminal-matching claim.

Several saved endpoints have binary exponents of enormous magnitude. The replacement consumer keeps endpoints as signed mantissa/exponent pairs and performs each addition and multiplication with outward rounding to at most 256 significant binary bits. It does not expand a power with an enormous exponent or flush a nonzero number to zero. In particular, an unresolved cancellation produces a wider interval; it is not used to infer a zero or a sign. Nineteen fresh synthetic checks passed, including negative log derivative, a nonzero \(F_0^2\) tail, huge exponent gaps, near cancellation, tiny terms, source mismatch and missing-error rejection. Source-family and certificate checks are data-derived; arithmetic soundness is conditional on the saved source bounds and the implementation assumptions recorded in the receipt.

The first, separately frozen v1 entry passed its 15 synthetic scale checks and five rounding checks, but stopped before producing any view because a saved endpoint exceeded its fixed Fraction expansion limit. Its complete failure, entry, source manifest and input hashes are preserved in the v1 failure receipt. The v2 consumer is a new frozen entry, with new synthetic inputs and a different exponent-preserving arithmetic implementation; the v1 entry was not replayed.

## Remaining work

The intervals can be extremely wide and are not actual function values. The \(\epsilon,F_0\) amplitude factors remain symbolic, and the six entries include H and K even though those labels were not among the 13 previously inspected Q4 frontiers. Q4's nested binder, the opaque core substitutions and the actual source-to-terminal function identity remain unproved. Neither interval overlap nor the sign of a single sample proves a five-moment terminal identity. Preserve all 212 accumulated terms, the six signed inner/lower terms, the reserved own-swirl contribution exactly once, the independent negative pulse, ordinary Z1, physical scales and all existing errors.

The next executable task is to connect the saved physical-scale function records and their exact source family to the Q4 nested substitutions without re-entering prior consumers, and then establish the same-source signed endpoint identities required before any inherited error can be removed. Admissible stress, flat remainder, true time recursion, oscillatory correction and independent Cartesian NS validation remain open.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| v1 failed receipt, raw gzip | 6fb9214c37f0d82ac996fb2a37cbc41047b46a43aecb3926a0ecbd42e35ee9fe |
| v1 complete failed result | fd6fed257cc09ab28b7d0b4d57d971080c5928c1596e53072723f4e5be58c078 |
| v2 receipt, raw gzip | e068e2d66c61aa7b7c41e210e97e91c2943253f949c4eb462721a7f4f2a505f8 |
| v2 frozen source manifest | d95b1f8eef34ef2e03c69e8e766097d1e0891b16835582de0a62e14e440d727d |
| v2 exclusive entry | 746c8f4f70c2b02e05fd82c5b2bcc2eaf5cfbb4f9cada75b75b08fc13aeeec71 |
| v2 complete result | 993344b80959b4ba7b8304c28e1edef9e985d55067bfee24767012458b07e7e2 |

> **Superseded: physical V_y scaling error.** The v1 evaluator and synthetic fixture below multiply the already physical raw background V_y by Pstar. That is incorrect; see the source normalization recorded in CURRENT_ORIGINAL_N1_FIRST_MICRO_SEAM_2026_10_17.md. The v1 numerical formulas and fixture results are not accepted as physical flux evidence. Preserve the frozen receipt as incident evidence. See the later correction handoff for the replacement.

# Signed lower source-packet flux evaluator — 2026-10-31

A new conditional evaluator accepts an already-held live first-micro lower source view and computes signed physical R Q_y and its ordinary Z1 enclosure, without replacing source quantities by magnitude caps or selected interval points. It is an implementation with synthetic validation only: no real source provider or original ancestor was executed in this entry, and no new numerical value or endpoint identity is claimed.

The input path follows the existing source_rows view: V, original m history, physical R, Pstar, and the source packet's actual_background_source.V_y. It checks live context, formal basis and ledger identity; direct object links back to the packet; the exact sc_half singleton query and common P0; and op.Pstar. It retains m and its derivatives. The caller must additionally supply source-qualified, packet-bound conditions for the requested Z cell, fixed R/Pstar/delta, the original history, and positivity of L. These conditions are not independently established merely because the caller supplies them. The evaluator also checks positivity of the supplied interval denominator directly.

The implementation converts the second Taylor rows to ordinary derivatives once and forms:

    L = 1 - delta Z^2
    A = U - m
    Ny = 2 Z U_y - (1-delta) Z A - (1-Z^2) A_Z
    Ny_Z = 2 U_y + 2 Z U_yZ - (1-delta)(A + Z A_Z)
           + 2 Z A_Z - (1-Z^2) A_ZZ
    R Q_y = R Ny / L
    d_Z(R Q_y) = R (Ny_Z L + 2 delta Z Ny) / L^2

Here U=Pstar V and U_y=Pstar V_y, with Pstar applied once to the unscaled V_y source. Radius and scale remain formal source enclosures; the result is not a point value. The ordinary Z1 expression assumes R, Pstar and delta are Z-independent on the requested cell. It does not prove that native initial ports [129567,129568] and corrected lower ports [53450,53462] are the same source function, nor that their independently computed interval enclosures may be subtracted as a correlated exact zero.

## Validation and scope

After freezing the source and input SHA manifest and recording the exclusive entry, one new synthetic fixture invocation returned all 17 checks true. The fixture checks exact rational enclosure arithmetic, preservation of nonzero m, ordinary Z1, wrong basis/query and missing-condition rejection. An independent finite difference of synthetic polynomials checks the derivative and the second-row factorial. Its candidate N=1024 is synthetic and has no bearing on the original candidate N. No real source evaluation, physical terminal matching, or error reduction occurred.

The source packet certifies the local finite-N flat collar, not vanishing background histories. The 212 ledger terms, six signed core/lower terms, separate negative complete-pulse target, and every previous error remain unchanged. The source-level reserved construction already defines one own_swirl slot and an independent pulse target; do not add another slot or pulse.

Next, build a new, frozen same-run bridge that obtains the native Q4 lower source and corrected lower packet from qualified owners, explicitly connects U, m, U_y, R, Pstar, delta, N, predecessor history and ordinary Z jets, and verifies both are the same defining source before asserting C1 cancellation. The present evaluator provides the signed packet-side calculation, but shared basis/ledger identity alone does not establish that cross-owner theorem. Then address the Rm-to-Rp physical germ and whole-Z reserved/pulse comparison.

## Reproducibility

| Artifact | SHA-256 |
| --- | --- |
| Receipt, raw gzip | d2975cf73f798ceef9989f9f79357e5094c887c91dc15d5a23bd63b2a2aa8ee5 |
| Frozen source manifest | 74b0da303ed8efdad35140a8c5356236438656b31872aba947d468a9c27503bc |
| Exclusive entry | b9236f3f38cca2dfa56c484aa7535b88c58c0362caf14df1d663951275699af1 |
| Complete result | 9aacc756d9e981d30e758bf56532be6134dc086d7b365803526aa264ec55d828 |
| Evaluator | 84926c3a437e56966e4f34e8e3d7a97f478f752df24995c7bbe6aadccc79c03d |
| Fixture | 5b5b31c95ec07c1bef5ff57847e25c9b2ca56a2ec21f65a5b499404ab6461b20 |
| Runner | 8c120c5a5b69da1c17746cb2b699560723f2d61ce59a03ba21dca0411207e685 |


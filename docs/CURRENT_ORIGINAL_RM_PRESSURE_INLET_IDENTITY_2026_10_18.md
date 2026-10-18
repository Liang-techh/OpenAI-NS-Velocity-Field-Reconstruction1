# Actual Rm pressure inlet and relative patch pressure source

Goal and Tasks 14/15 remain **ACTIVE / INCOMPLETE**. A new exact-source
consumer reduces the actual current-Rm and patch-left pressure expressions
through ordinary Z0/Z1. It cancels their identical S1 memory and initial
RQy function, leaving the signed physical-germ defect explicitly. It also
constructs the same-source pressure increment from Rm to Rh. No new pressure
or velocity numerical enclosure is computed, and no RQy equality is assumed.

## Actual source identity

The component is `original_Rm_pressure_inlet_identity.py`, with root runner
`scripts/run_original_Rm_pressure_inlet_identity.py`. The root authenticates
the saved Q4 graph (183,700 nodes), true leading patch (196,042 nodes), and
finite-N patch (205,147 nodes). The two smaller graphs are exactly equal to
the corresponding prefixes of the finite-N graph, including node metadata.
No source builder, owner, ancestor checker or numerical provider runs.

The current adopted pressure source retains Q4 Rm RQy roots183696/183697.
The patch-left boundary roots are204932/204933. Their common source family
and graph prefix are established provenance; they do not establish equality
between these distinct boundary functions.

Write Phi15 for the PRESSURE PRIMITIVE component

    Phi15(R) = S1(R)/2 - RQy(R) + RQy(initial).

This is the saved `pressure.m15` component. It is not the entire cumulative
m15 with all n1 channels, and it is not P1 or loop-root q. The new reducer
uses actual graph operands, exact rational arithmetic and scoped endpoint
substitutions. An integral is removed only when its lower and upper bounds
reduce to the identical source expression. A source alias must have the
same exact substitution expression, variable, value and Z-independent guard.

At patch x=1 the local S1 integral is exactly empty. The saved nonzero
incoming S1 roots183680/183681 and initial RQy129567/129568 are the same
functions used by the current Q4 pressure component. Thus

    Phi15_current(Rm) = S1_current/2 - B_current + B_initial
    Phi15_patch(Rm)   = S1_current/2 - B_patch   + B_initial
    Phi15_patch(Rm)-Phi15_current(Rm) = B_current-B_patch.

The last line is proved for actual source functions through ordinary Z1,
not for assigned sample values. The two B atoms remain independent with
coefficients +1 and -1. Neither a tiny normalized enclosure nor a sign cap
is selected as B or used to prove that their difference vanishes.

The actual patch Rh endpoint uses pressure roots205135/205136, S1 roots
205127/205128, and RQy roots205143/205144. The new exact increment is

    Phi15_patch(Rh)-Phi15_patch(Rm)
      = local_S1(Rh)/2 - B_patch(Rh) + B_patch(Rm).

The inherited S1 and initial boundary cancel as identical source functions
before any boxing. The local S1(Rh) remains an unevaluated actual definite
integral, identified by its exact Rh substitution descriptor. It is never
replaced by a bound or a computed center. This increment is a usable signed
function-source expression for later common-Z evaluation; it does not yet
transport the complete current m15 or configure the legacy forward chain.

The physical measures remain those of the saved patch: R=Rm*x and dy=dx/x;
the S1 integral uses Rm^2*x*dx exactly once. All expressions above use
physical source units. No extra Pstar, Rm, normalization, Z factorial or
F0 derivative is introduced. A later division by a frozen Rm^2 or common
unit must remain separate from this physical source identity.

## Validation and retained state

The new expression reducer ran once and was checkpointed before saved inspection.
Recorded execution seconds: {"new_actual_Rm_pressure_inlet_identity": 0.849186399980681}.
It traversed45 actual nodes, retained20 alias witnesses, and reduced exactly
the two empty local integrals204689/204690. No physical source value was
evaluated by that reduction.
The new source and identity checkpoints are authenticated in the report.
The saved-only inspector reads source bytes/checkpoints and preserves the
parent absolute-Rm and corrected own-swirl packets exactly; it does not
replay the new reducer. New Python passed LF/single-EOF/AST preflight.
The new synthetic case checks independent rational source arithmetic,
empty-integral cancellation, surviving boundary defect, metadata preservation
and invalid-edge/unsupported-operation guards. Its saved fixture is pinned
to the reviewed module; completed fixtures are not repeated.

Source preparation first corrected the parent's nested source-family path.
An import attempt then stopped before the operator because the dynamic
module was not registered for Python dataclasses. A subsequent binding-layout
validation stopped before reducer construction: the direct binding maps
also have an unrelated provenance key named `bindings`. The root now
supplies explicit archive wrappers `modified_Q4_pressure_C3_bindings` and
`Rm_patch_finite_N_C5_bindings`. These pre-reduction failures are preserved
in the report. The existing source checkpoint was reused after the loader
and wrapper fixes. An entry marker prevents replay after reduction starts;
the accepted mathematical expression reduction was not repeated.

The new module reduces pressure expressions, not the old numerical212-term
ledger. All212 forward terms, six structural channels (five signed plus
core-swirl error-only), every finite-N/width error, all old sample rows and
the single reserved registration remain unchanged. The corrected own-swirl
C1 source graph still has53,580 nodes; its numeric response remains cached
only at exact Z=371/1000,x=2. No old Z=-0.3,0,0.3 row is added to the current
Rm tube[0.37,0.372]. No numerical full-m15 terminal value is produced here.

The old finite-N inlet checker substituted common toy bodies into the
actual formula inputs. Its useful generic identities do not provide the
missing normalization certificate for the actual modified endpoint seeds.
The new source reducer does not inherit an old fingerprint's seam-closed
claim. It identifies precisely the remaining boundary function instead.

Three bounded tasks used fresh fixed `luna_max_worker` roles for the Rm
port audit, downstream source audit, and new module/fixture work; another
existing Luna worker performed read-only review. Tool results expose role
and task names, not an independent runtime model attestation.

## Remaining physical join

The current Rm->patch modified-field germ is still unproved. Actual seed
aliases are needed for E,V,m,h,k,e,p and E_y,V_y, with the required ordinary
Z and native derivatives, followed by the same global loop phase/inverse
and physical U,U_y,m,m_y,Q_y source map. Leading mean-only matching and
same family/N are insufficient. In particular the explicit current/patch
RQy difference above must remain until that source theorem is supplied.

There are further downstream gates. The finite-N Rh boundary205143/205144
and outer Rh-reference inlet have separate source traces. The saved native
join receipt does not prove their endpoint equality. Compact and quiet
carry their incoming function sources, but the lower active quiet seam is
also explicitly open. A common source-prefix chain through265,856 nodes
does not erase those physical-germ obligations.

NEXT prove the actual Rm modified endpoint seed normalization and first
native jets, then derive the pressure-germ equality from the shared loop
and physical map. Consume this turn's pressure defect and relative increment
without rerunning them. Evaluate new, admitted source consumers on a common
Z domain before combining numerical packets. Preserve every unresolved
signed boundary and all old errors. Functional five-moment closure, exact
heat exterior, admissible stress/flat remainder, true n-dependent recursion,
realizable corrections and independent Cartesian NS validation remain open.

Report experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_pressure_inlet_identity.json: 2345628 bytes,
SHA256 2a95b27e9851c225592c11e4481e3f8f4c92af317972d6521e45b29f1145e18b.

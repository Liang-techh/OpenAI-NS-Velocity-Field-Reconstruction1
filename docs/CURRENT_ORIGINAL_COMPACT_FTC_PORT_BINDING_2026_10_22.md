# Current original compact FTC port binding (FAILED, 2026-10-22)

This tag review records a failed source binding with preserved partial proof evidence. The four source recipe gates passed individually (native C1, native local S1, old C1 IBP, and endpoint/m15), while the aggregate source FTC transport remained false because the old Q and Q_y recovery residuals did not vanish. The saved result and six lossless artifacts remain byte exact; publication is metadata-only.

Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_compact_ftc_port_binding.json` (SHA256 `cc15d8dadfc3cd269f13a64f4ed48e5451be7ae4405af372ddcab95bb1ae2c70`). It retains 1891 source pins, 176 checkpoints, and a 50-artifact failed union.

The saved attempt performed 59 bounded source checks: 54 zero residuals and 5 nonzero residuals, specifically three old Q and two old Q_y residuals. Fourteen of sixteen proof groups passed; all nine native groups and the old I1/S1, endpoint, and m15 groups passed. PortCodec controls passed 12 of 12.

The semantic failure is a producer-label interpretation issue: the authenticated compact_functions producer allocates U 30709–30712, M 30713–30716, U_y 30717–30720, and M_y 30721–30724, while the consumer swapped the middle U_y/M labels. All 630 literal product guards were authenticated; no old-source physical bug is claimed.

Next task: Create a fresh unentered cached-witness semantic-port reconciliation. Authenticate compact_functions producer field order and source q/q_y paths, correct the U_y/M labels, and reuse the 54 passed plus 5 failed residual witnesses in a fresh five-difference proof only. Do not replay entered codec, proof, or normalizer code; do not evaluate source functions numerically or delete retained caps.

# Coupled stress probe: mismatch retained

Run lei_ren_part1_paper_continuous_bundle_stress_check.py for one shared prepared Z=.3 candidate. The two points are after axial support. The receipt records actual moments, pressure, analytic velocity jets, inertial/shear/total stress, and stress/shear ratios. Finiteness is the only runtime acceptance check, alongside presence of all five moments. Large ratios fail the intended admissible-stress matching; they are not a measured NS momentum residual.

Run lei_ren_part1_paper_continuous_pressure_tail_defect.py for a cheap audit of the saved output. It identifies the materialized pressure transport coefficient, without reconstructing suppressed subdominant terms. This audit is based on rounded receipts and does not establish pressure closure.

The probe now retains separate theta transport/mean/mixed and axial transport/linear/quadratic/pressure contributions, plus the exact coupled B_theta target. Summing the separated pieces replays the original inertial stresses at both points. At saved precision, axial pressure divided by I_z is 1 at both points. I_theta/abs(theta_transport) is about +0.0989369 at Rv+50 and -8.373795e-16 at Rtail+4; Rv+50 must not be assumed to have completed the flattening stage. Subdominant mixed terms are retained but tiny cancellation is not certified. No moment is replaced by the target.

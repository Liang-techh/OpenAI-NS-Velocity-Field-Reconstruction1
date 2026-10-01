"""Reproducible partial-core diagnostic, not a certified comparison."""
import hashlib
import json
from pathlib import Path
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    c = IntervalComparisonJets(4, diagnostic_partial=True)
    start = c.evaluate('0')
    end = c.evaluate('.01')
    driver = c.exit_driver_jets('.01', None, '1')
    frozen = c.evaluate('.03')
    if end['outside_analytic_core_radial_domain'] or not frozen['outside_analytic_core_radial_domain']:
        raise AssertionError('Frozen continuation core-certificate scope mislabeled')
    partial_rejected = False
    if c.degree < 124:
        try:
            IntervalComparisonJets(4)
        except ValueError as exc:
            if 'Production comparison requires' not in str(exc):
                raise
            partial_rejected = True
        if not partial_rejected:
            raise AssertionError('Partial core accepted for production')
    scalar_rejected = False
    try:
        c.evaluate(0, '.5')
    except ValueError:
        scalar_rejected = True
    if not scalar_rejected:
        raise AssertionError('Interval family mislabeled as scalar slice')
    here = Path(__file__).parent
    out = dict(completed_degree=c.degree, state_sha256=c.state_hash,
        center_family=c.center_family, steps=4, diagnostic_partial=True,
        production_partial_rejected=partial_rejected, scalar_slice_rejected=scalar_rejected,
        outside_domain_frozen_continuation_labeled=True,
        comparison_driver_only=True, stress_derivatives_supplied=False,
        start_D=list(start['D'].coefficients), endpoint_D=list(end['D'].coefficients),
        endpoint_I_z=list(end['I_z'].coefficients),
        driver_A=list(driver['A'].coefficients), driver_B=list(driver['B'].coefficients),
        analytic_tail_propagated=False, ODE_discretization_error_enclosed=False,
        temporal_recursion=False, whole_axis_transition_generated=False,
        source_hashes={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in
            ('lei_ren_part1_paper_interval_comparison_jets.py',
             'lei_ren_part1_paper_candidate_comparison_jets.py',
             'lei_ren_part1_paper_candidate_shared_inlet.py',
             'lei_ren_part1_paper_interval_comparison_jets_check.py')})
    (here/'lei_ren_part1_paper_interval_comparison_jets_diagnostic.json').write_text(
        json.dumps(encode(out),indent=2)+'\n', encoding='utf-8')
    print('Fresh interval comparison diagnostic at degree', c.degree,
          '; partial production and scalar slice rejected:', partial_rejected, scalar_rejected, flush=True)
    return out


if __name__ == '__main__':
    run()

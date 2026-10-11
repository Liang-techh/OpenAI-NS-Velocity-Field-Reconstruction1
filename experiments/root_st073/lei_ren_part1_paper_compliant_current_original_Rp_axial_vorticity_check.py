"""Independent cylindrical curl law and actual original signed mu result."""
import copy
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_axial_vorticity as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:
        fn()
    except (ValueError, TypeError, ArithmeticError, KeyError):
        return True
    raise AssertionError('Invalid original axial vorticity source admitted')


def independent_curl_identity():
    r, theta, R, lam, mu, delta, A = s.symbols('r theta R lambda mu delta A', positive=True)
    vr, vt = s.Function('vr')(r), s.Function('vt')(r)
    ux, uy = vr*s.cos(theta)-vt*s.sin(theta), vr*s.sin(theta)+vt*s.cos(theta)
    dx = lambda value: s.cos(theta)*s.diff(value, r)-s.sin(theta)/r*s.diff(value, theta)
    dy = lambda value: s.sin(theta)*s.diff(value, r)+s.cos(theta)/r*s.diff(value, theta)
    assert s.trigsimp(s.expand(dx(uy)-dy(ux)-s.diff(vt, r)-vt/r)) == 0
    theta_source = A*R**(-s.Rational(1, 2)-mu)
    physical_swirl = lam**(-1-delta)*theta_source
    physical_radius = lam*s.sqrt(2*R)
    curl = s.diff(physical_radius*physical_swirl, R)/s.diff(physical_radius, R)/physical_radius
    expected = -mu*s.sqrt(2/R)*lam**(-2-delta)*theta_source
    assert s.simplify(curl-expected) == 0
    assert s.simplify(s.diff(theta_source, R)*R+theta_source/2+mu*theta_source) == 0
    return dict(Cartesian_to_cylindrical_axial_curl=True, original_theta_source_rate=True,
        original_physical_radial_and_lambda_factors=True, exact_mu_retained_before_subtraction=True)


@source_precision
def run(owner, values, fields):
    began = time.monotonic()
    raw = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.OPEN)
    assert raw['source_family'] == owner.family_record and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name] == digest for name, digest in owner.before.hashes.items())
    symbolic = independent_curl_identity()
    c = MPIntervalContext(); c.dps = 800
    observations = {}
    for name, value in values.items():
        field = fields[name]
        report = owner.report(value)
        record = owner.before._require(field)
        issued = owner.arithmetic._require(value)
        reader = owner.before._reader(record)
        bases = [term for component, index, q in current.differential.standard_forms()['vorticity_z']
            for term in record['entries'][component, index]['raw'].terms
            if term.source_label == current.correlated.physical.mixed.UT and term.source_row.derivative == (0, 0)]
        assert bases and all(term.source_row.coefficients[0]._mpi_ == bases[0].source_row.coefficients[0]._mpi_ for term in bases)
        base = bases[0]
        source = c.mpf(base.source_row.coefficients[0])
        independent = -c.sqrt(2)*source
        actual = owner.arithmetic._coefficient(issued)
        assert current.differential.product.contains(actual, independent)
        log = owner.graph.add(*(ref for _, ref in base.log_scale_parts), owner.product.logmu)
        assert not reader.polynomial(owner.graph.sub(log, issued['scale']))
        old = owner.before._linear(record, current.differential.standard_forms()['vorticity_z'])
        old_scale = current.box.pulse.radius.FunctionRef(owner.graph, old['exact_reference_log_scale_function'])
        base_scale = owner.graph.add(*(ref for _, ref in base.log_scale_parts))
        assert not reader.polynomial(owner.graph.sub(base_scale, old_scale))
        refined_mu = c.mpf(owner.product.inclusion['mu']['fresh_box'])
        assert current.differential.product.contains(old['common_scale_coefficient_enclosure'], independent*refined_mu)
        assert current.ends(source)[0] > 0 and current.ends(actual)[1] < 0
        assert issued['unit'] == 'velocity/length**1' and report['value']['sign'] == -1
        assert report['value']['proved_nonzero'] and report['original_mu_retained_before_enclosure']
        assert report['near_equal_derivatives_not_subtracted'] and report['exact_radial_primitive_cancellations'] == 2
        assert not report['value']['ordinary_numeric_materialized'] and report['actual_time_growth_not_measured']
        assert not report['unrestricted_physical_point_API'] and not report['full_certified_physical_accuracy']
        assert report['physical_accuracy']['ordinary_numeric_relative_width_satisfied']
        assert report['value']['original_absolute_width_target_satisfied']
        assert report['value']['exact_source_inverse_mu']['coefficient']['numerator'] == 0
        assert all(not report[key] for key in current.OPEN)
        assert current.ends(owner.arithmetic.ratio(value, value)['ordinary_numeric_ratio_enclosure']) == (1, 1)
        observations[name] = dict(sign=-1, proved_nonzero=True, physical_unit=issued['unit'],
            exact_source_inverse_mu_coefficient=0, total_relative_width_target_satisfied=True,
            exact_retained_logmu_function=owner.product.logmu.node,
            physical_accuracy=report['physical_accuracy'], ordinary_numeric_materialized=False,
            time_growth_recursion_and_global_scope_open=True)
    value = next(iter(values.values())); field = next(iter(fields.values()))
    invalid = dict(copied_value=rejected(lambda: owner.report(copy.copy(value))),
        copied_field=rejected(lambda: owner.evaluate(copy.copy(field))))
    saved = owner.refined.pulse.main
    try:
        owner.refined.pulse.main = owner.refined.pulse.end
        invalid['substituted_actual_refined_source_method'] = rejected(owner.assert_graph)
    finally:
        owner.refined.pulse.main = saved
    saved = owner.refined.pulse.mu
    try:
        owner.refined.pulse.mu = owner.refined.ctx.mpf(saved)*2
        invalid['changed_actual_mu_source_enclosure'] = rejected(owner.assert_graph)
    finally:
        owner.refined.pulse.mu = saved
    proof = owner._issued[id(value)][3]['exact_source_law_proof']
    saved = copy.deepcopy(owner.graph.nodes[proof])
    try:
        owner.graph.nodes[proof]['identity'] = 'omega_z=0'
        invalid['changed_defining_curl_law'] = rejected(lambda: owner.report(value))
    finally:
        owner.graph.nodes[proof] = saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes = dict(owner.hashes)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(all_passed=True, source_family=owner.family_record, actual_theta_source_law=owner.law,
        independent_symbolic_curl_identities=symbolic, independent_800_digit_theta_coefficient_inclusion=True,
        recovered_mu_curl_included_in_original_unresolved_enclosure=True,
        actual_original_axial_vorticity_observations=observations, invalid_inputs_rejected=invalid,
        original_mu_pressure_scales_and_physical_units_retained=True,
        no_time_growth_stress_flat_recursion_or_global_claim=True,
        input_hashes=hashes, execution_seconds=time.monotonic()-began,
        **dict.fromkeys(current.GATES, True), **dict.fromkeys(current.OPEN, False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt), indent=2)+'\n', encoding='utf8', newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_AXIAL_VORTICITY exact original mu and signed curl', flush=True)
    return receipt

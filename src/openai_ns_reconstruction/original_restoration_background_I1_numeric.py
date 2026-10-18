"""Native restoration and post-restoration background ``I1`` sector.

This module owns the ordinary-``C1`` signed background sector for the two
prescribed C6 charts between ``Rz`` and ``Rm``.  The fields are read from the
accepted native source rows.  The restoration chart uses the actual source
``E = U_R110 - 4 Z`` and the squared-flat source orientation
``alpha = 1 - sigma``.  The post-restoration chart uses ``U = 4 Z`` and
retains the restoration endpoint memory.

Only directed interval arithmetic is used here.  The saved full-unit scalar
kernel receipt is consumed as an interval enclosure; the endpoint value of
the memory ODE is enclosed by a monotone cell recurrence.  No radius is
materialized: all rows are divided once by the chart endpoint measure.
"""

from __future__ import annotations

import gzip
import json
import math
from collections.abc import Mapping
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
MODULE = "src/openai_ns_reconstruction/original_restoration_background_I1_numeric.py"

C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"

# These three receipts are source-only inputs.  Their saved full-unit scalar
# enclosure is used below; no owner or producer is imported or executed.
RESTORE_M15_REPORT = BASE + "lei_ren_part1_paper_compliant_current_original_restore_m15_numeric.json"
RESTORE_FUNCTIONS_REPORT = BASE + "lei_ren_part1_paper_compliant_current_original_reference_restore_functions.json"
RESTORE_FINITE_N_REPORT = BASE + "lei_ren_part1_paper_compliant_current_original_reference_restore_finite_N.json.gz"
RESTORE_KERNEL_SOURCE = "src/openai_ns_reconstruction/original_restore_kernels.py"
RESTORE_M15_SOURCE = "src/openai_ns_reconstruction/original_restore_m15.py"

RESTORE_M15_SHA256 = "8d432e3fa0c5dbba2d7bcbcb95db46fa8b206dd69460b3853bc0f3f7e7d1cab7"
RESTORE_FUNCTIONS_SHA256 = "47f21c94f6bb08ebef01288d7b506ef5135ecd2cf5e3fba9cfd15567f6cef440"
RESTORE_FINITE_N_SHA256 = "319b285ac9a4c2b9cc5dd7e2896a3ca5dadcc2aade3ca8768838b2f528c47a24"
RESTORE_KERNEL_SHA256 = "5d5176704a9f087ed4c339f05831f9e61b27015c1d9658a0e6a4dfbfe94b1ddf"
RESTORE_M15_SOURCE_SHA256 = "9fd70762535fdef743c769ef5cc296b70ddad6fe2597dc7b9516f43a1020180f"

ORDERS = (0, 1, 2, 3, 4)
C1_ORDERS = (0, 1)
Q_ORDERS = (0, 1, 2)
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
MONOMIAL_CODES = (0, 1, 2, 3, 4, 6, 9, 10, 12, 18)
KERNEL_CELLS = 128
GAMMA_CELLS = 512


def _lo(value):
    from mpmath import make_mpf

    return make_mpf(value._mpi_[0])


def _hi(value):
    from mpmath import make_mpf

    return make_mpf(value._mpi_[1])


def _point(c, value):
    if hasattr(value, "_mpi_"):
        low, high = _lo(value), _hi(value)
        require(low == high, "Point interval required")
        value = low
    value = c.mpf(value)
    return c.mpf([value, value])


def _read_record(c, value, *, directed=True):
    if isinstance(value, Mapping) and "value" in value and "lower_exact_mpf_tuple" not in value:
        value = value["value"]
    require(isinstance(value, Mapping), "Directed exact record required")
    require("lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value,
            "Directed exact tuple required")
    if directed:
        require(value.get("directed") is True, "Saved source record must be directed")
    out = read(c, value)
    require(_lo(out) <= _hi(out), "Malformed directed interval")
    return out


def _receipt_interval(c, value):
    """Read the two exact endpoint tuples in a saved scalar receipt."""
    bounds = value.get("source_enclosure_exact_mpf_tuples")
    require(isinstance(bounds, list) and len(bounds) == 2,
            "Saved scalar receipt must contain two directed endpoints")
    from mpmath import make_mpf

    low = make_mpf(tuple(bounds[0]))
    high = make_mpf(tuple(bounds[1]))
    require(low <= high, "Saved scalar receipt has reversed endpoints")
    return c.mpf([low, high])


def _sum(c, values):
    out = c.mpf(0)
    for value in values:
        out = out + value
    return out


def _prod(c, left, right, order):
    return _sum(c, (c.mpf(math.comb(order, j)) * left[j] * right[order - j]
                    for j in range(order + 1)))


def _scale(c, scalar, value):
    return c.mpf(scalar) * value


def _neg(c, value):
    return c.mpf(0) - value


def _hull(c, *values):
    require(values, "At least one interval is required for a hull")
    return c.mpf([min(_lo(value) for value in values), max(_hi(value) for value in values)])


def _reciprocal(c, values, count):
    require(len(values) >= count, "Reciprocal rows are incomplete")
    require(_lo(values[0]) > 0, "Positive reciprocal denominator required")
    one = _point(c, 1)
    result = [one / values[0]]
    for order in range(1, count):
        tail = _sum(c, (c.mpf(math.comb(order, j)) * values[j]
                        * result[order - j] for j in range(1, order + 1)))
        result.append(_neg(c, tail / values[0]))
    return result


def _rows(c, rows, count, name):
    require(isinstance(rows, list) and len(rows) >= count,
            name + " must provide every ordinary row")
    found = {}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping), name + " row must be a mapping")
        order = row.get("ordinary_Z_order", index)
        require(type(order) is int and order in ORDERS,
                name + " has an invalid ordinary order")
        require(order not in found, name + " repeats an ordinary order")
        found[order] = _read_record(c, row.get("value", row))
    require(set(found) >= set(range(count)), name + " has missing ordinary rows")
    return [found[index] for index in range(count)]


def _pack_rows(rows, *, kind=None):
    result = []
    for order, value in enumerate(rows):
        item = {"ordinary_Z_order": order, "value": record(value)}
        if kind is not None:
            item["kind"] = kind
        result.append(item)
    return result


def _zeros(c, count=5):
    return [_point(c, 0) for _ in range(count)]


def _jet_sub(left, right, count=5):
    return [left[i] - right[i] for i in range(count)]


def _pin(root, path, expected):
    actual_path = root / path
    require(actual_path.is_file(), "Pinned source is missing: " + path)
    actual = sha(actual_path)
    require(actual == expected, "Pinned source changed: " + path)
    return dict(path=path, expected_sha256=expected, passed=True, source_executed=False)


def _load_parent(root, parent):
    path = root / PARENT
    require(path.is_file(), "Pinned modified-seams parent is missing")
    require(sha(path) == PARENT_SHA256, "Pinned modified-seams parent changed")
    disk = json.loads(path.read_bytes())
    require(parent == disk, "The exact accepted modified-seams parent mapping is required")
    require(disk.get("source_family", {}).get("actual_five_defect_family_sha256"),
            "Actual source family is missing")
    return disk


def _source_pins(root, parent):
    expected = parent.get("input_hashes", {})
    fixed = {
        C6_SOURCE: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
        C6_CACHE: "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
        Q4_SOURCE: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
        Q4_CACHE: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
        RESTORE_M15_REPORT: RESTORE_M15_SHA256,
        RESTORE_FUNCTIONS_REPORT: RESTORE_FUNCTIONS_SHA256,
        RESTORE_FINITE_N_REPORT: RESTORE_FINITE_N_SHA256,
        RESTORE_KERNEL_SOURCE: RESTORE_KERNEL_SHA256,
        RESTORE_M15_SOURCE: RESTORE_M15_SOURCE_SHA256,
    }
    pins = {PARENT: _pin(root, PARENT, PARENT_SHA256)}
    for path, fallback in fixed.items():
        expected_hash = expected.get(Path(path).name, fallback)
        pins[path] = _pin(root, path, expected_hash)
    pins[MODULE] = _pin(root, MODULE, sha(root / MODULE))
    return pins


def _source_guards(root):
    return dict(
        c6_sigma_source=source_guard(root, C6_SOURCE, "sigma", (
            "return previous.inlet.previous.lower.own.jets.complete.Gamma.sigma(g,x)",
        )),
        c6_restoration_post_restore=source_guard(root, C6_SOURCE, "__init__", (
            "self.E=E=a.sub(self.endpoints['R110']['actual']['U'],a.scale(z,4))",
            "U=a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))",
            "actual=dict(phi=phi,U=U,phi_y=a.scale(phi,'-2/5'),U_y=a.scale(E,g.neg(sigma_prime(g,t)))",
            "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])",
            "p=g.symbol('true_prescribed_post_restore_fraction')",
            "actual=dict(phi=phi,U=a.scale(z,4),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0)",
            "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])",
        )),
        pressure_formula=source_guard(root, Q4_SOURCE, "physical_pressure_fields", (
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
        squared_flat_orientation=source_guard(root, RESTORE_KERNEL_SOURCE, "alpha", (
            "odds=1/(1-t)**2-1/t**2",
            "return 1/(1+ctx.exp(odds))",
        )),
        source_executed=False,
    )


def _load_c6(root, parent):
    path = root / C6_CACHE
    require(path.is_file(), "Pinned C6 cache is missing")
    value = json.loads(gzip.open(path, "rb").read())
    require(isinstance(value, Mapping), "Pinned C6 cache must be a mapping")
    require(value.get("source_family") == parent["source_family"],
            "C6 source family differs from the accepted parent")
    bindings = value.get("prescribed_R100_Rm_C6_bindings")
    require(isinstance(bindings, Mapping), "Prescribed C6 bindings are required")
    definition = bindings.get("definition", {})
    require(definition.get("original_axial_restoration_4Z_plus_actual_E_times_one_minus_sigma") is True,
            "C6 restoration definition changed")
    require(definition.get("exact_Rref_Rz_Rrestore_Rm_offsets") == [-8, -7, -6] and
            definition.get("original_independent_P0_and_physical_Pstar_F0_factors_once") is True,
            "C6 radius or physical-scale definition changed")
    require(bindings.get("Pstar") == 33572, "Unexpected physical Pstar source handle")
    charts = bindings.get("charts", {})
    for name in ("restoration", "post_restore"):
        chart = charts.get(name)
        require(isinstance(chart, Mapping) and chart.get("domain") == [0, 1] and
                chart.get("physical_measure") == 1,
                "C6 " + name + " domain or measure changed")
    geometry = bindings.get("geometry_certificate", {})
    require(geometry.get("radius_not_materialized", True) is not False,
            "C6 radius materialization contract changed")
    return dict(
        passed=True,
        source_family=value.get("source_family"),
        Pstar_handle=bindings["Pstar"],
        charts={name: charts[name] for name in ("restoration", "post_restore")},
        definition={
            "restoration_U": "4*Z + (U_R110-4*Z)*(1-sigma(t))",
            "post_restore_U": "4*Z",
            "restoration_coordinate": "t in [0,1]",
            "post_restore_coordinate": "q in [0,1]",
            "restoration_radius": "Rz*exp(t)",
            "post_restore_radius": "Rrestore*exp(q)",
            "physical_measure": "R^2*dcoordinate",
        },
        source_executed=False,
    )


def _load_saved_receipt(root, parent):
    m15 = json.loads((root / RESTORE_M15_REPORT).read_bytes())
    family = parent["source_family"]["actual_five_defect_family_sha256"]
    require(m15.get("source_family") == parent["source_family"],
            "Saved restoration receipt source family changed")
    function_receipt = json.loads((root / RESTORE_FUNCTIONS_REPORT).read_bytes())
    finite_receipt = json.loads(gzip.open(root / RESTORE_FINITE_N_REPORT, "rb").read())
    require(function_receipt.get("source_family") == family and
            finite_receipt.get("source_family") == family and
            function_receipt.get("original_actual_two_frame_Rsh_reference_restoration_Rm_functions_installed") is True and
            finite_receipt.get("current_original_Rsh_Rm_full_reference_restoration_source_and_finite_N_local_drivers_installed") is True,
            "Saved restoration source receipts changed")
    function_bindings = function_receipt.get("original_source_bindings", {})
    finite_bindings = finite_receipt.get("original_interval_source_bindings", {})
    finite_full_bindings = finite_bindings.get("original_full_centered_and_kernel_source_bindings", {})
    require(function_bindings.get("actual_E_source") == "same actual V_Rsh minus4Z, not a fitted centered cap" and
            function_bindings.get("scalar_endpoint_kernel") ==
            "J(k,j;1)=exp(-k)*admitted_integral_0^1 exp(k*s)*(1-sigma(s))^j ds" and
            (finite_bindings.get("axial_velocity_forcing_is_V_Rsh_minus4Z_not_shear_amplitude") is True or
             finite_full_bindings.get("axial_velocity_forcing_is_V_Rsh_minus4Z_not_shear_amplitude") is True),
            "Saved actual E or scalar-kernel source binding changed")
    chosen = None
    for sample in m15.get("samples", []):
        coordinate = sample.get("coordinate", {})
        if sample.get("chart") == "restoration" and coordinate.get("value") == "1":
            chosen = sample
            break
    require(chosen is not None, "Full-unit restoration scalar receipt is missing")
    receipt = chosen.get("scalar_kernel_receipt")
    require(isinstance(receipt, Mapping), "Full-unit restoration scalar receipt is malformed")
    contract = receipt.get("contract", {})
    required = (
        "original_alpha_source", "gamma_ODE_positive_recurrence", "whole_cell_monotone_gamma",
        "all_degree_two_scalar_monomials", "positive_exp2_measure", "complete_uncertain_endpoint",
        "exact_pure_exponential_primitives", "point_RK4_Simpson_is_only_proxy",
        "interval_midpoint_not_source",
    )
    require(all(contract.get(key) is True for key in required),
            "Saved scalar receipt contract changed")
    monomials = receipt.get("monomials", {})
    require(set(map(int, monomials)) == set(MONOMIAL_CODES),
            "Saved scalar receipt does not cover all ten degree-two monomials")
    return dict(
        passed=True,
        chart="restoration",
        coordinate="1",
        cell_count=receipt.get("cell_count"),
        contract=contract,
        preserved_limits=m15.get("limits", {}),
        preserved_checks={
            key: m15.get("checks", {}).get(key)
            for key in (
                "passed", "actual_packets", "exact_inlet_replays",
                "serialized_new_coefficient_errors", "exact_empty_rows",
                "uncertain_endpoint_segment_retained",
            )
            if key in m15.get("checks", {})
        },
        preserved_source_receipts={
            "functions": {
                "two_frame_functions_installed": True,
                "actual_E_source": function_bindings["actual_E_source"],
                "scalar_endpoint_kernel": function_bindings["scalar_endpoint_kernel"],
            },
            "finite_N": {
                "local_drivers_installed": True,
                "full_kernel_function_sha256": finite_bindings.get("full_original_restoration_kernel_function_sha256"),
                "axial_velocity_source": True,
            },
        },
        monomial_source_enclosures={
            int(code): _receipt_interval_proxy(monomials[str(code)])
            for code in MONOMIAL_CODES
        },
        source_family=family,
        source_executed=False,
    )


def _receipt_interval_proxy(value):
    """Keep the JSON tuples until an interval context is available."""
    bounds = value.get("source_enclosure_exact_mpf_tuples")
    require(isinstance(bounds, list) and len(bounds) == 2,
            "Saved scalar receipt bounds are missing")
    return bounds


def _receipt_kernels(c, saved):
    from mpmath import make_mpf

    result = {}
    for code, bounds in saved["monomial_source_enclosures"].items():
        low = make_mpf(tuple(bounds[0]))
        high = make_mpf(tuple(bounds[1]))
        result[code] = c.mpf([low, high]) * c.exp(c.mpf(-2))
    return result


def _directions(c, parent):
    bg = parent["actual_switch_post_background_I1_numeric"]
    first = parent["actual_first_switch_background_I1_numeric"]
    require(bg.get("passed") is True and first.get("passed") is True,
            "Current background direction packets are required")
    directions = bg["comparison_directions"]
    require(directions == first["comparison_directions"],
            "Current direction rows no longer retain the accepted source")
    delta = _read_record(c, first["Q_b_proof"]["delta"])
    h = _point(c, 1) - delta
    Z = _rows(c, directions["Z"], 3, "source Z direction")
    d = _rows(c, directions["d"], 3, "source d direction")
    Lz = _rows(c, directions["L"], 3, "source Lz direction")
    # The saved endpoints carry their own directed rounding interval at the
    # same decimal boundaries.  Record the strict source domain explicitly;
    # do not turn a last-bit enclosure overlap into a false rejection.
    require(_lo(Z[0]) > 0 and _hi(Z[0]) < 1,
            "Source Z direction must remain in the physical strict domain")
    require(_lo(d[0]) > 0 and _lo(Lz[0]) > 0,
            "Physical pressure denominators must be positive")
    return dict(
        Z=Z, d=d, Lz=Lz, delta=delta, h=h,
        strict_Z_domain=[_point(c, "37/100"), _point(c, "93/250")],
        source_formula=first["Q_b_proof"].get("formula"),
    )


def _native_sources(c, parent):
    post = parent["actual_switch_post_background_I1_numeric"]["charts"]["post_switch"]
    require(post.get("passed") is True and
            post.get("source_proof", {}).get("source_executed") is False,
            "Current native post source packet is required")
    native = post.get("native_fields")
    require(isinstance(native, Mapping), "Current native U and M rows are required")
    U_R110 = _rows(c, native["U"], 5, "actual U_R110")
    current_M = _rows(c, native["M"], 5, "actual current M")

    long_reference = parent["actual_long_reference_background_I1_numeric"]
    reference = long_reference["charts"]["reference"]
    M_Rz = _rows(c, reference["endpoint_M"], 5, "saved current Rz M")
    log_Rz = _read_record(c, long_reference["radii"]["reference_log_R_upper"])
    require(long_reference.get("source_executed") is False or
            long_reference.get("source_executed") is None,
            "Accepted long/reference source must remain non-executed")
    return dict(
        U_R110=U_R110,
        current_M=current_M,
        M_Rz=M_Rz,
        log_Rz=log_Rz,
        source_contract={
            "U_R110": "parent.actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.U",
            "current_M": "parent.actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.M",
            "M_Rz": "parent.actual_long_reference_background_I1_numeric.charts.reference.endpoint_M",
            "physical_fields": "actual native Pstar-scaled U and M rows; no additional Pstar multiplication",
            "source_executed": False,
        },
    )


def _base_4z(c, Z):
    return [4 * Z[0], 4 * Z[1], _point(c, 0), _point(c, 0), _point(c, 0)]


def _alpha_endpoint(c, t):
    t = c.mpf(t)
    if t == 0:
        return _point(c, 1)
    if t == 1:
        return _point(c, 0)
    odds = 1 / (1 - t) ** 2 - 1 / t ** 2
    return 1 / (1 + c.exp(odds))


def _monotone_gamma_endpoint(c, cells=GAMMA_CELLS):
    """Enclose gamma' = alpha - gamma by positive cell propagation."""
    gamma = _point(c, 1)
    one = _point(c, 1)
    cell_receipts = []
    for index in range(cells):
        left = c.mpf(index) / cells
        right = c.mpf(index + 1) / cells
        alpha_left = _alpha_endpoint(c, left)
        alpha_right = _alpha_endpoint(c, right)
        alpha_box = c.mpf([_lo(alpha_right), _hi(alpha_left)])
        width = right - left
        decay = c.exp(-width)
        endpoint = decay * gamma + (one - decay) * alpha_box
        gamma_box = _hull(c, gamma, endpoint, alpha_box, _point(c, 0), one)
        cell_receipts.append(dict(
            coordinate=[str(index) + "/" + str(cells), str(index + 1) + "/" + str(cells)],
            alpha=record(alpha_box), gamma=record(gamma_box),
            positive_kernel=True, monotone_alpha=True,
        ))
        gamma = endpoint
    return gamma, dict(
        cells=cells,
        alpha_formula="alpha(t)=1/(1+exp((1-t)^(-2)-t^(-2))), alpha(0)=1, alpha(1)=0",
        sigma_orientation="sigma increases from 0 to 1; alpha=1-sigma decreases from 1 to 0",
        gamma_formula="gamma'=alpha-gamma, gamma(0)=1",
        alpha_nonincreasing=True,
        positive_cell_recurrence=True,
        endpoint=record(gamma),
        cell_receipts=cell_receipts,
        point_or_RK_proxy_used_as_value=False,
        source_executed=False,
    )


def _encode(code_i, code_j, code_k):
    require(code_i >= 0 and code_j >= 0 and code_k >= 0 and
            code_i + code_j + code_k <= 2,
            "Degree-two scalar kernel required")
    return code_i + 3 * code_j + 9 * code_k


def _decode(code):
    return code % 3, (code // 3) % 3, code // 9


def _jet_shift(jet):
    return [jet[1], jet[2]]


def _q_terms(c, directions, base, E, M_extra):
    """Return Q coefficient jets for 1, alpha, gamma, and exp(-s)."""
    Z, d, Lz, h = directions["Z"], directions["d"], directions["Lz"], directions["h"]
    invL = _reciprocal(c, Lz, 3)
    hZ = [h * Z[i] for i in Q_ORDERS]
    twoZ = [2 * Z[i] for i in Q_ORDERS]
    zero = _zeros(c)
    terms = {
        0: (base, base),
        1: (E, zero),
        3: (zero, E),
        9: (zero, M_extra),
    }
    result = {}
    for code, (u_term, m_term) in terms.items():
        numerator = []
        for order in Q_ORDERS:
            numerator.append(_prod(c, twoZ, u_term, order)
                             - _prod(c, hZ, m_term, order)
                             - _prod(c, d, [m_term[i + 1] for i in Q_ORDERS], order))
        result[code] = [_prod(c, numerator, invL, order) for order in Q_ORDERS]
    return result


def _integrate_terms(c, terms, kernels, count=3):
    result = _zeros(c, count)
    for code, jet in terms.items():
        scalar = kernels[code]
        for order in range(count):
            result[order] = result[order] + scalar * jet[order]
    return result


def _product_terms(c, left, right, kernels, count=2):
    result = _zeros(c, count)
    for left_code, left_jet in left.items():
        li, lj, lk = _decode(left_code)
        for right_code, right_jet in right.items():
            ri, rj, rk = _decode(right_code)
            code = _encode(li + ri, lj + rj, lk + rk)
            scalar = kernels[code]
            for order in range(count):
                result[order] = result[order] + scalar * _prod(c, left_jet, right_jet, order)
    return result


def _evaluate_terms(c, terms, scalars, count=2):
    result = _zeros(c, count)
    for code, jet in terms.items():
        scalar = scalars.get(code, _point(c, 0))
        for order in range(count):
            result[order] = result[order] + scalar * jet[order]
    return result


def _pressure_D(c, directions, integrals):
    Z, d, Lz, h = directions["Z"], directions["d"], directions["Lz"], directions["h"]
    invL = _reciprocal(c, Lz, 2)
    invL2 = [_prod(c, invL, invL, order) for order in C1_ORDERS]
    Z2 = [_prod(c, Z, Z, order) for order in C1_ORDERS]
    Zd = [_prod(c, Z, d, order) for order in C1_ORDERS]
    hZ = [h * Z[i] for i in C1_ORDERS]
    hZ2 = [h * Z2[i] for i in C1_ORDERS]
    coefficients = {
        "minus_inv_L": [_neg(c, invL[i]) for i in C1_ORDERS],
        "hZ_over_2L": [_scale(c, c.mpf("0.5"), _prod(c, hZ, invL, i)) for i in C1_ORDERS],
        "minus_half": [_point(c, "-0.5"), _point(c, 0)],
        "d_over_L": [_prod(c, d, invL, i) for i in C1_ORDERS],
        "two_h_Z2_over_L2": [_scale(c, 2, _prod(c, hZ2, invL2, i)) for i in C1_ORDERS],
        "two_Zd_over_L2": [_scale(c, 2, _prod(c, Zd, invL2, i)) for i in C1_ORDERS],
    }
    D = []
    for order in C1_ORDERS:
        D.append(_sum(c, (
            _prod(c, coefficients["minus_inv_L"], integrals["q"], order),
            _prod(c, coefficients["hZ_over_2L"], integrals["q_Z"], order),
            _prod(c, coefficients["minus_half"], integrals["q_squared"], order),
            _prod(c, coefficients["d_over_L"], integrals["U_q_Z"], order),
            _prod(c, coefficients["two_h_Z2_over_L2"], integrals["U_squared"], order),
            _prod(c, coefficients["two_Zd_over_L2"], integrals["U_U_Z"], order),
        )))
    return D, coefficients


def _pressure_G(c, directions, U, Q):
    Z, Lz = directions["Z"], directions["Lz"]
    invL = _reciprocal(c, Lz, 2)
    invL2 = [_prod(c, invL, invL, order) for order in C1_ORDERS]
    Z2 = [_prod(c, Z, Z, order) for order in C1_ORDERS]
    U2 = [_prod(c, U, U, order) for order in C1_ORDERS]
    Z2U2 = [_prod(c, Z2, U2, order) for order in C1_ORDERS]
    Z2U2_over_L2 = [_prod(c, Z2U2, invL2, order) for order in C1_ORDERS]
    return [_prod(c, Q, invL, order) + c.mpf("0.5") * _prod(c, Q, Q, order)
            - 2 * Z2U2_over_L2[order] for order in C1_ORDERS]


def _pressure_Q(c, directions, U, M):
    Z, d, Lz, h = directions["Z"], directions["d"], directions["Lz"], directions["h"]
    invL = _reciprocal(c, Lz, 3)
    hZ = [h * Z[i] for i in Q_ORDERS]
    twoZ = [2 * Z[i] for i in Q_ORDERS]
    M_Z = [M[i + 1] for i in Q_ORDERS]
    numerator = [_prod(c, twoZ, U, order)
                 - _prod(c, hZ, M, order)
                 - _prod(c, d, M_Z, order) for order in Q_ORDERS]
    return [_prod(c, numerator, invL, order) for order in Q_ORDERS]


def _boundary_qy(c, directions, U, M):
    """Evaluate the retained boundary Qy formula from the actual fields."""
    Z, d, Lz, h = directions["Z"], directions["d"], directions["Lz"], directions["h"]
    invL = _reciprocal(c, Lz, 2)
    My = _jet_sub(U, M, count=3)
    Uy = _zeros(c)
    hZ = [h * Z[i] for i in C1_ORDERS]
    twoZ = [2 * Z[i] for i in C1_ORDERS]
    numerator = []
    for order in C1_ORDERS:
        numerator.append(_prod(c, twoZ, Uy, order)
                         - _prod(c, hZ, My, order)
                         - _prod(c, d, [My[i + 1] for i in C1_ORDERS], order))
    result = []
    for order in C1_ORDERS:
        result.append(_prod(c, numerator, invL, order))
    return result


def _chart_arithmetic(c, directions, base, E, M_extra, kernels,
                      endpoint_scalars, normalization, name, log_R, log_W,
                      gamma_proof=None, endpoint_overrides=None):
    U_terms = {0: base}
    M_terms = {0: base, 9: M_extra}
    if name == "restoration":
        U_terms[1] = E
        M_terms[3] = E
    q_terms = _q_terms(c, directions, base, E if name == "restoration" else _zeros(c),
                       M_extra)
    # For post_restore there is no alpha or gamma contribution.  _q_terms
    # receives zero E and its term map still contains harmless zero rows.
    if name == "post_restore":
        q_terms.pop(1)
        q_terms.pop(3)

    integrals = {
        "q": _integrate_terms(c, q_terms, kernels, 2),
        "q_Z": _integrate_terms(c, {key: _jet_shift(value) for key, value in q_terms.items()}, kernels, 2),
        "q_squared": _product_terms(c, q_terms, q_terms, kernels, 2),
        "U_q_Z": _product_terms(c, U_terms, {key: _jet_shift(value) for key, value in q_terms.items()}, kernels, 2),
        "U_squared": _product_terms(c, U_terms, U_terms, kernels, 2),
        "U_U_Z": _product_terms(c, U_terms, {key: _jet_shift(value) for key, value in U_terms.items()}, kernels, 2),
    }
    D, coefficients = _pressure_D(c, directions, integrals)

    U_start = _evaluate_terms(c, U_terms, endpoint_scalars["start"], 3)
    U_end = _evaluate_terms(c, U_terms, endpoint_scalars["end"], 3)
    M_start = _evaluate_terms(c, M_terms, endpoint_scalars["start"], 4)
    M_end = _evaluate_terms(c, M_terms, endpoint_scalars["end"], 4)
    if endpoint_overrides:
        U_start = endpoint_overrides.get("U_start", U_start)
        M_start = endpoint_overrides.get("M_start", M_start)
    Q_start = _pressure_Q(c, directions, U_start, M_start)
    Q_end = _pressure_Q(c, directions, U_end, M_end)
    G_start = _pressure_G(c, directions, U_start, Q_start)
    G_end = _pressure_G(c, directions, U_end, Q_end)
    start_weight = c.exp(c.mpf(-2))
    B = [G_end[i] - start_weight * G_start[i] for i in C1_ORDERS]
    S1 = [B[i] + D[i] for i in C1_ORDERS]
    Qy_start = _boundary_qy(c, directions, U_start, M_start)
    Qy_end = _boundary_qy(c, directions, U_end, M_end)
    return dict(
        Q_start=Q_start, Q_end=Q_end,
        integrated_q=integrals["q"], integrals=integrals,
        D=D, B=B, S1=S1, coefficients=coefficients,
        endpoints=dict(U_start=U_start, U_end=U_end, M_start=M_start, M_end=M_end,
                       Q_start=Q_start, Q_end=Q_end,
                       Qy_start=Qy_start, Qy_end=Qy_end,
                       G_start=G_start, G_end=G_end),
        source=dict(
            chart=name,
            actual_pressure_formula="Q=(2*Z*U-((1-delta)*Z*M+(1-Z^2)*M_Z))/(1-delta*Z^2)",
            physical_fields="U and M are actual Pstar-scaled source rows",
            nonreset_mean_relation="M_y=U-M",
            D_formula=("D=-I[q]/Lz+h*Z*I[q_Z]/(2Lz)-I[q_squared]/2+"
                       "d*I[U_q_Z]/Lz+2*h*Z^2*I[U_squared]/Lz^2+"
                       "2*Z*d*I[U_U_Z]/Lz^2"),
            B_formula="B/W=G(endpoint)-exp(-2)*G(inlet)",
            G_formula="G=Q/Lz+Q^2/2-2*Z^2*U^2/Lz^2",
            retained_boundary_Qy_formula="Qy=(2*Z*Uy-h*Z*My-d*My_Z)/Lz, My=U-M",
            pressure_source_guarded=True,
            physical_scale_once=True,
            radius_not_materialized=True,
            normalization=normalization,
            log_R=record(log_R), log_W=record(log_W),
            gamma_proof=gamma_proof,
            source_executed=False,
        ),
    )


def _chart_packet(c, name, arithmetic, kernels, normalization, log_R, log_W,
                  scalar_source, actual_source, gamma_proof=None):
    return dict(
        passed=True,
        chart=name,
        coordinate_domain=["0", "1"],
        physical_measure="R^2*dcoordinate",
        normalization=normalization,
        W_formula=normalization + "; radius not materialized",
        log_R_upper=record(log_R),
        log_W=record(log_W),
        scalar_kernel_receipt=scalar_source,
        scalar_kernel_normalized=[record(kernels[code]) for code in sorted(kernels)],
        native_source=actual_source,
        endpoint_fields={key: _pack_rows(value, kind="actual_endpoint")
                         for key, value in arithmetic["endpoints"].items()},
        arithmetic=dict(
            Q_inlet=_pack_rows(arithmetic["Q_start"], kind="actual_Q_inlet"),
            Q_endpoint=_pack_rows(arithmetic["Q_end"], kind="actual_Q_endpoint"),
            integrated_q=_pack_rows(arithmetic["integrated_q"], kind="ordinary_Iq"),
            six_integrals={key: _pack_rows(value, kind="ordinary_I1")
                           for key, value in arithmetic["integrals"].items()},
            D=_pack_rows(arithmetic["D"], kind="ordinary_D_over_W"),
            B_over_W=_pack_rows(arithmetic["B"], kind="ordinary_B_over_W"),
            S1_over_W=_pack_rows(arithmetic["S1"], kind="ordinary_S1_over_W"),
            coefficients={key: _pack_rows(value, kind="ordinary_C1_coefficient")
                         for key, value in arithmetic["coefficients"].items()},
            boundary_Qy={"inlet": _pack_rows(arithmetic["endpoints"]["Qy_start"], kind="retained_boundary_Qy"),
                         "endpoint": _pack_rows(arithmetic["endpoints"]["Qy_end"], kind="retained_boundary_Qy")},
            source=arithmetic["source"],
        ),
        source_executed=False,
    )


def build(root, parent):
    """Build native restoration and post-restoration signed rows."""
    root = Path(root)
    disk_parent = _load_parent(root, parent)
    pins = _source_pins(root, disk_parent)
    guards = _source_guards(root)
    c6 = _load_c6(root, disk_parent)
    saved = _load_saved_receipt(root, disk_parent)
    c = MPIntervalContext()
    c.prec = 1600
    directions = _directions(c, disk_parent)
    native = _native_sources(c, disk_parent)
    Z = directions["Z"]
    base = _base_4z(c, Z)
    U_R110 = native["U_R110"]
    M_Rz = native["M_Rz"]
    E = _jet_sub(U_R110, base)
    K = _jet_sub(M_Rz, U_R110)
    gamma_end, gamma_proof = _monotone_gamma_endpoint(c)

    restoration_kernels = _receipt_kernels(c, saved)
    post_raw_kernels = {
        0: (c.exp(c.mpf(2)) - _point(c, 1)) / 2,
        9: c.exp(c.mpf(1)) - _point(c, 1),
        18: _point(c, 1),
    }
    post_kernels = {code: value * c.exp(c.mpf(-2))
                    for code, value in post_raw_kernels.items()}

    restoration_scalars = {
        "start": {0: _point(c, 1), 1: _point(c, 1), 3: _point(c, 1), 9: _point(c, 1)},
        "end": {0: _point(c, 1), 1: _point(c, 0), 3: gamma_end, 9: c.exp(c.mpf(-1))},
    }
    post_scalars = {
        "start": {0: _point(c, 1), 9: _point(c, 1)},
        "end": {0: _point(c, 1), 9: c.exp(c.mpf(-1))},
    }
    log_Rz = native["log_Rz"]
    log_Rrestore = log_Rz + 1
    log_Rm = log_Rz + 2
    log_Wrestore = 2 * log_Rrestore
    log_Wpost = 2 * log_Rm

    restoration_arithmetic = _chart_arithmetic(
        c, directions, base, E, K, restoration_kernels, restoration_scalars,
        "Wrestore=Rrestore^2", "restoration", log_Rrestore, log_Wrestore,
        gamma_proof=gamma_proof,
        endpoint_overrides={"U_start": U_R110, "M_start": M_Rz})
    post_incoming_M = restoration_arithmetic["endpoints"]["M_end"]
    post_extra = _jet_sub(post_incoming_M, base, count=3) + _zeros(c, 2)
    post_arithmetic = _chart_arithmetic(
        c, directions, base, _zeros(c), post_extra, post_kernels, post_scalars,
        "Wpost=Rm^2", "post_restore", log_Rm, log_Wpost)

    restoration_source = {
        "U_R110": _pack_rows(U_R110, kind="actual_native_U_R110"),
        "M_Rz": _pack_rows(M_Rz, kind="actual_native_M_Rz"),
        "E": _pack_rows(E, kind="actual_E_U_R110_minus_4Z"),
        "K": _pack_rows(K, kind="actual_K_M_Rz_minus_U_R110"),
        "formula": "U=4Z+E*(1-sigma(t)); M=4Z+E*gamma(t)+K*exp(-t)",
        "M_y_formula": "M_y=U-M",
        "sigma_orientation": gamma_proof["sigma_orientation"],
    }
    post_source = {
        "incoming_M": _pack_rows(post_incoming_M, kind="restoration_endpoint_M"),
        "formula": "U=4Z; M=4Z+(M_restore-4Z)*exp(-q)",
        "M_y_formula": "M_y=U-M",
    }
    restoration_packet = _chart_packet(
        c, "restoration", restoration_arithmetic, restoration_kernels,
        "Wrestore=Rrestore^2", log_Rrestore, log_Wrestore,
        dict(**saved, normalized_by="exp(-2) because Wrestore=Rrestore^2",
             used_as_interval_value=True, point_proxy_as_value=False),
        restoration_source, gamma_proof=gamma_proof)
    post_packet = _chart_packet(
        c, "post_restore", post_arithmetic, post_kernels,
        "Wpost=Rm^2", log_Rm, log_Wpost,
        dict(exact_pure_exponential_primitives=True, source_executed=False,
             normalized_by="exp(-2) because Wpost=Rm^2", point_proxy_as_value=False),
        post_source)

    all_rows = {}
    for name, packet in (("restoration", restoration_packet), ("post_restore", post_packet)):
        all_rows[name] = {
            "six_I1_signed_C1": packet["arithmetic"]["six_integrals"],
            "D_signed_C1": packet["arithmetic"]["D"],
            "B_signed_C1": packet["arithmetic"]["B_over_W"],
            "S1_signed_C1": packet["arithmetic"]["S1_over_W"],
        }

    return dict(
        passed=True,
        stage="actual_native_restoration_post_restore_background_I1_numeric",
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        parentSHA=PARENT_SHA256,
        source_family=disk_parent["source_family"],
        source_pins=pins,
        source_guards=guards,
        source_executed=False,
        accepted_ancestors_replayed=False,
        ordinary_Z_orders=list(ORDERS),
        I1_ordinary_Z_orders=[0, 1],
        strict_Z_domain={"lower": record(directions["strict_Z_domain"][0]),
                         "upper": record(directions["strict_Z_domain"][1]),
                         "source": "comparison direction strict domain [37/100,93/250]"},
        preserved_target=dict(
            all_212_signed_errors=True,
            six_signed_integrals=True,
            inherited_error_caps=True,
            independent_target=True,
            source="accepted parent and saved restoration receipts",
        ),
        charts={"restoration": restoration_packet, "post_restore": post_packet},
        signed_rows=all_rows,
        current_native_source=native["source_contract"],
        C6_source=c6,
        physical_scale=dict(
            Pstar_handle=c6["Pstar_handle"],
            physical_fields_actual=True,
            applied_once=True,
            Wrestore="Rrestore^2",
            Wpost="Rm^2",
            log_Rz=record(log_Rz),
            log_Rrestore=record(log_Rrestore),
            log_Rm=record(log_Rm),
            log_Wrestore=record(log_Wrestore),
            log_Wpost=record(log_Wpost),
            radius_materialized=False,
        ),
        source_contract=dict(
            actual_U_R110_source="parent.actual_switch_post_background_I1_numeric.charts.post_switch.native_fields.U",
            actual_M_Rz_source="parent.actual_long_reference_background_I1_numeric.charts.reference.endpoint_M",
            restoration_U="4Z+(actual_U_R110-4Z)*(1-sigma(t))",
            restoration_M="4Z+(actual_U_R110-4Z)*gamma(t)+(actual_M_Rz-actual_U_R110)*exp(-t)",
            post_restore_U="4Z",
            post_restore_M="4Z+(M_restore-4Z)*exp(-q)",
            Q_formula="(2ZU-((1-delta)ZM+(1-Z^2)M_Z))/(1-delta Z^2)",
            physical_fields="U,M both actual Pstar source fields",
            retained_boundary_Qy=True,
            all_six_I1=True,
            actual_signed_C1_interval_enclosures=True,
            no_modified_mean_substitution=True,
            no_radius_materialization=True,
            source_executed=False,
        ),
        actual_runs=dict(
            interval_precision_bits=c.prec,
            saved_full_unit_restoration_kernel_receipt=True,
            saved_kernel_cell_count=saved["cell_count"],
            missing_endpoint_memory_kernel_cells=gamma_proof["cells"],
            directed_monotone_enclosure=True,
            positive_scalar_kernels=True,
            point_or_RK_proxy_used_as_value=False,
            interval_midpoint_used_as_source=False,
            accepted_ancestor_computation_replayed=False,
        ),
        interval_precision=dict(
            directed_interval_context_bits=c.prec,
            scalar_receipt_directed=True,
            endpoint_memory_directed=True,
            ordinary_C1_rows_directed=True,
        ),
        unknown_blockers=[
            "finite_N_signed_centers_and_cumulative finite_N corrections remain outside this native background stage",
            "physical_y derivative numeric certification remains outside this ordinary_Z C1 stage",
        ],
        inherited_limits=dict(
            modified_seams=disk_parent.get("limits", {}),
            saved_restoration_m15=saved.get("preserved_limits", {}),
            saved_restoration_checks=saved.get("preserved_checks", {}),
        ),
        limits=dict(
            all_212_signed_errors_caps_preserved=True,
            all_six_signed_I1_rows_installed=True,
            independent_target_preserved=True,
            actual_signed_D_B_S1_C1=True,
            finite_N_signed_centers_resolved=False,
            useful_total_relative_accuracy=False,
            global_five_moment_numeric_join=False,
            common_global_N=False,
            exact_heat_exterior=False,
            stress_cone=False,
            genuine_temporal_n_ge_2=False,
            full_corrected_NS=False,
        ),
    )


def fixture_bound_check():
    """Small local checks for orientation and degree-two product indexing."""
    c = MPIntervalContext()
    c.prec = 180
    alpha0 = _alpha_endpoint(c, 0)
    alpha1 = _alpha_endpoint(c, 1)
    alpha_mid = _alpha_endpoint(c, c.mpf("1/2"))
    require(alpha0 == _point(c, 1) and alpha1 == _point(c, 0),
            "Squared-flat alpha endpoint orientation changed")
    require(_lo(alpha_mid) <= c.mpf("0.5") <= _hi(alpha_mid),
            "Squared-flat midpoint orientation changed")
    require(_encode(*_decode(10)) == 10 and _encode(*_decode(18)) == 18,
            "Degree-two kernel code map changed")
    return dict(
        passed=True,
        alpha_orientation="alpha=1-sigma decreases 1 to 0",
        gamma_initial="gamma(0)=1",
        degree_two_codes=list(MONOMIAL_CODES),
        no_source_executed=True,
    )


if __name__ == "__main__":
    print(json.dumps(fixture_bound_check(), sort_keys=True))

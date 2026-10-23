"""Standalone symbolic graph construction for the native volume error.

This module only builds a fresh integer DAG.  It does not import a source
producer, evaluate a source function, run a proof, normalize a bound, or
perform an integral.  The outer consumer owns source admission and cap flags.
"""

from fractions import Fraction


_MISSING = object()

# Archive identities are explicit source contracts, separate from source-row
# extraction.  They are never counted as default source-contract rows.
_NATIVE_ARCHIVE_RAW_CONTRACT = (
    "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"
)
_OLD_ARCHIVE_RAW_CONTRACT = (
    "3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce"
)
_NATIVE_ARCHIVE_PATH_CONTRACT = (
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
    "true_compact_C5_functions.json.gz"
)
_OLD_ARCHIVE_PATH_CONTRACT = (
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
    "n1_preRp_moment_functions.json.gz"
)


def _freeze(value):
    if isinstance(value, dict):
        return tuple(sorted((str(key), _freeze(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _fraction(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, str):
        return Fraction(value)
    if isinstance(value, dict):
        return Fraction(int(value["numerator"]), int(value["denominator"]))
    raise TypeError("graph constants must be exact integer or Fraction values")


def _fraction_json(value):
    value = _fraction(value)
    return {"numerator": int(value.numerator), "denominator": int(value.denominator)}


def _path_get(record, path):
    value = record
    for part in path.split(".") if isinstance(path, str) else path:
        if not isinstance(value, dict) or part not in value:
            return _MISSING
        value = value[part]
    return value


def _pick(record, paths, default):
    for path in paths:
        value = _path_get(record, path)
        if value is not _MISSING and value is not None:
            return value, path
    return default, None


def _row_int(value):
    if isinstance(value, dict):
        for key in ("node", "row", "source_row", "value"):
            if key in value and isinstance(value[key], (int, str)):
                value = value[key]
                break
    if isinstance(value, bool):
        raise TypeError("boolean is not a source row")
    return int(value)


def _rows(value, count, label):
    if isinstance(value, dict) and "rows" in value:
        value = value["rows"]
    if not isinstance(value, (list, tuple)) or len(value) != count:
        raise ValueError("%s must contain exactly %d ordinary-Z rows" % (label, count))
    return [_row_int(item) for item in value]


def _node_int(value, label):
    if isinstance(value, dict):
        for key in ("node", "row", "source_row", "value"):
            if key in value:
                value = value[key]
                break
    return _row_int(value)


class _DAG:
    """Small deterministic integer DAG with JSON-safe node records."""

    def __init__(self):
        self.nodes = []
        self._interned = {}
        self.counters = {
            "add": 0,
            "sub": 0,
            "mul": 0,
            "div": 0,
            "neg": 0,
            "source_ref": 0,
            "cached_ref": 0,
            "bound_variable": 0,
            "formal_integral": 0,
            "sparse_linear": 0,
            "c1_product_rules": 0,
            "c1_quotient_rules": 0,
            "c2_product_rules": 0,
        }

    def _intern(self, key, record):
        old = self._interned.get(key)
        if old is not None:
            return old
        node_id = len(self.nodes)
        record = dict(record)
        record["id"] = node_id
        self.nodes.append(record)
        self._interned[key] = node_id
        operation = record.get("op")
        if operation in self.counters:
            self.counters[operation] += 1
        return node_id

    def const(self, value):
        value = _fraction(value)
        key = ("const", value.numerator, value.denominator)
        return self._intern(
            key,
            {"op": "const", "value": _fraction_json(value)},
        )

    def source(self, row, chart, label, order):
        row = _row_int(row)
        key = ("source_ref", row, chart, label, int(order))
        return self._intern(
            key,
            {
                "op": "source_ref",
                "source_row": row,
                "chart": chart,
                "label": label,
                "ordinary_Z_order": int(order),
            },
        )

    def cached(self, row, label, order):
        row = _row_int(row)
        key = ("cached_ref", row, label, int(order))
        return self._intern(
            key,
            {
                "op": "cached_ref",
                "cached_node": row,
                "label": label,
                "ordinary_Z_order": int(order),
                "replay": False,
            },
        )

    def local(self, name):
        key = ("bound_variable", name)
        return self._intern(
            key,
            {"op": "bound_variable", "name": name, "ordinary_Z_derivative": 0},
        )

    def op(self, operation, args, **metadata):
        args = tuple(int(item) for item in args)
        frozen_metadata = _freeze(metadata)
        key = (operation, args, frozen_metadata)
        record = {"op": operation, "args": list(args)}
        record.update(metadata)
        return self._intern(key, record)


def _c1(value, dz):
    return {"value": int(value), "dz": int(dz)}


def _c2(value, dz, dzz):
    return {"value": int(value), "dz": int(dz), "dzz": int(dzz)}


def _c1_const(graph, value):
    zero = graph.const(0)
    return _c1(graph.const(value), zero)


def _c2_const(graph, value):
    zero = graph.const(0)
    node = graph.const(value)
    return _c2(node, zero, zero)


def _c1_add(graph, left, right):
    return _c1(
        graph.op("add", [left["value"], right["value"]]),
        graph.op("add", [left["dz"], right["dz"]]),
    )


def _c1_sub(graph, left, right):
    return _c1(
        graph.op("sub", [left["value"], right["value"]]),
        graph.op("sub", [left["dz"], right["dz"]]),
    )


def _c1_neg(graph, value):
    return _c1(
        graph.op("neg", [value["value"]]),
        graph.op("neg", [value["dz"]]),
    )


def _c1_mul(graph, left, right):
    graph.counters["c1_product_rules"] += 1
    value = graph.op("mul", [left["value"], right["value"]])
    dz = graph.op(
        "add",
        [
            graph.op("mul", [left["dz"], right["value"]]),
            graph.op("mul", [left["value"], right["dz"]]),
        ],
    )
    return _c1(value, dz)


def _c1_div(graph, numerator, denominator):
    graph.counters["c1_quotient_rules"] += 1
    value = graph.op("div", [numerator["value"], denominator["value"]])
    denominator_sq = graph.op("mul", [denominator["value"], denominator["value"]])
    top = graph.op(
        "sub",
        [
            graph.op("mul", [numerator["dz"], denominator["value"]]),
            graph.op("mul", [numerator["value"], denominator["dz"]]),
        ],
    )
    return _c1(value, graph.op("div", [top, denominator_sq]))


def _c1_sum(graph, values):
    result = _c1_const(graph, 0)
    for value in values:
        result = _c1_add(graph, result, value)
    return result


def _c2_sub(graph, left, right):
    return _c2(
        graph.op("sub", [left["value"], right["value"]]),
        graph.op("sub", [left["dz"], right["dz"]]),
        graph.op("sub", [left["dzz"], right["dzz"]]),
    )


def _c2_mul(graph, left, right):
    graph.counters["c2_product_rules"] += 1
    value = graph.op("mul", [left["value"], right["value"]])
    dz = graph.op(
        "add",
        [
            graph.op("mul", [left["dz"], right["value"]]),
            graph.op("mul", [left["value"], right["dz"]]),
        ],
    )
    twice = graph.const(2)
    cross = graph.op("mul", [twice, graph.op("mul", [left["dz"], right["dz"]])])
    dzz = graph.op(
        "add",
        [
            graph.op("mul", [left["dzz"], right["value"]]),
            graph.op(
                "add",
                [cross, graph.op("mul", [left["value"], right["dzz"]])],
            ),
        ],
    )
    return _c2(value, dz, dzz)


def _c2_source(graph, rows, chart, label):
    return _c2(
        graph.source(rows[0], chart, label, 0),
        graph.source(rows[1], chart, label, 1),
        graph.source(rows[2], chart, label, 2),
    )


def _c1_geometry_source(graph, row, chart, label, derivative):
    return _c1(
        graph.source(row, chart, label, 0),
        graph.const(derivative),
    )


def _extract_port_rows(record, paths, default, count, label):
    value, path = _pick(record, paths, default)
    return _rows(value, count, label), path, path is None


def _extract_node(record, paths, default, label):
    value, path = _pick(record, paths, default)
    return _node_int(value, label), path, path is None


def _source_c2(graph, rows, chart, label):
    return _c2_source(graph, rows, chart, label)


def _scaled_c2(graph, scale, normalized):
    return _c2_mul(graph, scale, normalized)


def _chart_geometry(graph, chart, rc_row, delta_row, z_row, coordinate_row, x_node):
    # Keep the local-coordinate binding as an integer DAG id at the caller,
    # while accepting the explicit C1 lift required by the pullback rules.
    # This guard also keeps the helper total for direct static callers.
    if isinstance(x_node, int):
        x_node = _c1(x_node, graph.const(0))
    one = _c1_const(graph, 1)
    two = _c1_const(graph, 2)
    z = _c1_geometry_source(graph, z_row, chart, "Z", 1)
    delta = _c1_geometry_source(graph, delta_row, chart, "delta", 0)
    rc = _c1_geometry_source(graph, rc_row, chart, "Rc", 0)
    z2 = _c1_mul(graph, z, z)
    d = _c1_sub(graph, one, z2)
    h = _c1_sub(graph, one, delta)
    L = _c1_sub(graph, one, _c1_mul(graph, delta, z2))
    L2 = _c1_mul(graph, L, L)
    radius = _c1_mul(graph, rc, x_node)
    mu = _c1_div(graph, one, x_node)
    radius2 = _c1_mul(graph, radius, radius)
    weight = _c1_mul(graph, radius2, mu)
    coefficients = {
        "q": _c1_div(graph, _c1_const(graph, Fraction(-1, 2)), L),
        "qZ": _c1_div(
            graph,
            _c1_mul(graph, h, z),
            _c1_mul(graph, _c1_const(graph, 4), L),
        ),
        "q2": _c1_const(graph, Fraction(-1, 4)),
        "UqZ": _c1_div(graph, d, _c1_mul(graph, _c1_const(graph, 2), L)),
        "U2": _c1_div(graph, _c1_mul(graph, h, z2), L2),
        "UUZ": _c1_div(graph, _c1_mul(graph, z, d), L2),
    }
    return {
        "chart": chart,
        "Z": z,
        "delta": delta,
        "Rc": rc,
        "d": d,
        "h": h,
        "L": L,
        "L2": L2,
        "radius": radius,
        "mu": mu,
        "weight": weight,
        "coefficients": coefficients,
        "coordinate_source": graph.source(
            coordinate_row, chart, "archive_coordinate_handle", 0
        ),
        "x": x_node,
        "two": two,
    }


def _chart_v(graph, geometry, fields):
    U = fields["U"]
    Q = fields["Q"]
    # qz is a C1 jet whose value is Q_Z and whose derivative is Q_ZZ.
    qz = _c1(Q["dz"], Q["dzz"])
    u = _c1(U["value"], U["dz"])
    uz = _c1(U["dz"], U["dzz"])
    q = _c1(Q["value"], Q["dz"])
    products = {
        "Q": q,
        "QZ": qz,
        "Q2": _c1_mul(graph, q, q),
        "UQZ": _c1_mul(graph, u, qz),
        "U2": _c1_mul(graph, u, u),
        "UUZ": _c1_mul(graph, u, uz),
    }
    coeffs = geometry["coefficients"]
    terms = {
        "q_linear": _c1_mul(graph, coeffs["q"], products["Q"]),
        "z_qz_linear": _c1_mul(graph, coeffs["qZ"], products["QZ"]),
        "q_square": _c1_mul(graph, coeffs["q2"], products["Q2"]),
        "u_qz": _c1_mul(graph, coeffs["UqZ"], products["UQZ"]),
        "u_square": _c1_mul(graph, coeffs["U2"], products["U2"]),
        "u_uz": _c1_mul(graph, coeffs["UUZ"], products["UUZ"]),
    }
    return {
        "products": products,
        "terms": terms,
        "v": _c1_sum(graph, list(terms.values())),
    }


def _formal_integral(graph, value, label):
    metadata = {
        "label": label,
        "lower": 1,
        "upper": 2,
        "measure": "dx",
        "derivative_under_integral": True,
    }
    value_metadata = dict(metadata)
    value_metadata["ordinary_Z_order"] = 0
    derivative_metadata = dict(metadata)
    derivative_metadata["ordinary_Z_order"] = 1
    return _c1(
        graph.op("formal_integral", [value["value"]], **value_metadata),
        graph.op("formal_integral", [value["dz"]], **derivative_metadata),
    )


def _cached_c1(graph, endpoint, signed_record):
    value, value_path = _pick(
        endpoint,
        ["value", "value_node", "C1.value"],
        _MISSING,
    )
    dz, dz_path = _pick(
        endpoint,
        ["dz", "dz_node", "C1.dz"],
        _MISSING,
    )
    if value is _MISSING or dz is _MISSING:
        raise ValueError(
            "cached endpoint requires actual value and ordinary-Z derivative nodes"
        )
    return _c1(
        graph.cached(value, "endpoint_difference_x2_minus_x1", 0),
        graph.cached(dz, "endpoint_difference_x2_minus_x1", 1),
    ), {"value_path": value_path, "dz_path": dz_path, "replay": False}


def _sparse_c1(graph, atoms, coefficients):
    names = list(atoms)
    coeffs = [_fraction(value) for value in coefficients]
    terms = []
    value_args = []
    dz_args = []
    for name, coefficient in zip(names, coeffs):
        terms.append({"atom": name, "coefficient": _fraction_json(coefficient)})
        value_args.append(atoms[name]["value"])
        dz_args.append(atoms[name]["dz"])
    metadata = {
        "atom_order": names,
        "atom_terms": terms,
        "all_atom_arguments_retained": True,
        "zero_coefficient_atoms_retained": True,
        "exact_merge": True,
        "integration_linearity_only": True,
    }
    value = graph.op("sparse_linear", value_args, **metadata)
    dz = graph.op("sparse_linear", dz_args, **metadata)
    return _c1(value, dz), terms


def _as_c1_dict(value):
    return {"value": int(value["value"]), "dz": int(value["dz"])}


def build_volume_source_graph(native_record, signed_record):
    """Build the native N/BG and registered-old signed volume source graph.

    The function consumes source row metadata only.  It returns a JSON-safe
    mapping of integer DAG nodes and C1/C2 node pairs.  Cap rows, midpoint
    values, source evaluation, quadrature, proofs, and old checkers are not
    inputs to this builder.
    """
    graph = _DAG()

    native_u_rows, native_u_path, native_u_default = _extract_port_rows(
        native_record,
        [
            "field_ports.actual_physical.U.rows",
            "field_ports.actual_physical.U_C2",
            "native.field_ports.actual_physical.U.rows",
            "true_compact_C5_bindings.physical.U",
            "native.true_compact_C5_bindings.physical.U",
            "physical.U",
        ],
        [262269, 262270, 262271],
        3,
        "native physical U",
    )
    native_q_rows, native_q_path, native_q_default = _extract_port_rows(
        native_record,
        [
            "field_ports.actual_physical.Q.rows",
            "field_ports.actual_physical.Q_C2",
            "native.field_ports.actual_physical.Q.rows",
            "true_compact_C5_bindings.physical.Q",
            "native.true_compact_C5_bindings.physical.Q",
            "physical.Q",
        ],
        [262364, 262368, 262374],
        3,
        "native physical Q",
    )
    bg_v_rows, bg_v_path, bg_v_default = _extract_port_rows(
        native_record,
        [
            "field_ports.pure_background_leading.normalized_V_jet.rows",
            "native.field_ports.pure_background_leading.normalized_V_jet.rows",
            "true_compact_C5_bindings.background.fields.V",
            "native.true_compact_C5_bindings.background.fields.V",
            "background.fields.V",
        ],
        [0, 0, 0],
        3,
        "normalized background V",
    )
    bg_q_rows, bg_q_path, bg_q_default = _extract_port_rows(
        native_record,
        [
            "field_ports.pure_background_leading.Q_b.rows",
            "field_ports.pure_background_leading.Q.rows",
            "native.field_ports.pure_background_leading.Q_b.rows",
            "native.field_ports.pure_background_leading.Q.rows",
            "true_compact_C5_bindings.background.generic.Q",
            "native.true_compact_C5_bindings.background.generic.Q",
            "background.generic.Q",
        ],
        [260923, 260927, 260933],
        3,
        "normalized background generic Q",
    )

    native_rc, native_rc_path, native_rc_default = _extract_node(
        native_record,
        [
            "native_view.geometry.Rc",
            "native_view.geometry.Rc.node",
            "geometry.Rc.node",
            "native.geometry.Rc.node",
            "true_compact_C5_bindings.Rc",
            "native.true_compact_C5_bindings.Rc",
            "Rc",
        ],
        251066,
        "native Rc",
    )
    native_delta, native_delta_path, native_delta_default = _extract_node(
        native_record,
        [
            "native_view.geometry.delta",
            "native_view.geometry.delta.node",
            "geometry.delta.node",
            "native.geometry.delta.node",
            "true_compact_C5_bindings.delta",
            "native.true_compact_C5_bindings.delta",
            "delta",
        ],
        25220,
        "native delta",
    )
    native_z, native_z_path, native_z_default = _extract_node(
        native_record,
        [
            "native_view.geometry.Z",
            "native_view.geometry.Z.node",
            "geometry.Z.node",
            "native.geometry.Z.node",
            "true_compact_C5_bindings.Z",
            "native.true_compact_C5_bindings.Z",
            "Z",
        ],
        20,
        "native Z",
    )
    native_coordinate, native_coordinate_path, native_coordinate_default = _extract_node(
        native_record,
        [
            "native_view.geometry.local_coordinate",
            "native_view.geometry.local_coordinate.node",
            "geometry.local_coordinate.node",
            "native.geometry.local_coordinate.node",
            "true_compact_C5_bindings.coordinate",
            "native.true_compact_C5_bindings.coordinate",
            "coordinate",
        ],
        260181,
        "native coordinate",
    )
    pstar, pstar_path, pstar_default = _extract_node(
        native_record,
        [
            "native_view.geometry.Pstar",
            "native_view.geometry.Pstar.node",
            "geometry.Pstar.node",
            "native.geometry.Pstar.node",
            "true_compact_C5_bindings.Pstar",
            "native.true_compact_C5_bindings.Pstar",
            "Pstar",
        ],
        33572,
        "native Pstar",
    )

    old_source, old_source_path = _pick(
        signed_record,
        ["component_result.source_interface.old", "source_interface.old", "old"],
        {},
    )
    old_u_rows, old_u_path, old_u_default = _extract_port_rows(
        old_source,
        ["U", "U_C2", "physical_U"],
        [30709, 30710, 30711],
        3,
        "old physical U",
    )
    old_q_rows, old_q_path, old_q_default = _extract_port_rows(
        old_source,
        ["Q", "Q_C2", "physical_Q"],
        [30758, 30762, 30768],
        3,
        "old physical Q",
    )
    old_rc, old_rc_path, old_rc_default = _extract_node(
        old_source,
        ["Rc_node", "Rc"],
        30828,
        "old Rc",
    )
    old_delta, old_delta_path, old_delta_default = _extract_node(
        old_source,
        ["delta_node", "delta"],
        25220,
        "old delta",
    )
    old_z, old_z_path, old_z_default = _extract_node(
        old_source,
        ["Z_node", "Z"],
        20,
        "old Z",
    )
    old_coordinate, old_coordinate_path, old_coordinate_default = _extract_node(
        old_source,
        ["coordinate_node", "coordinate"],
        5556,
        "old coordinate",
    )

    endpoint, endpoint_path = _pick(
        signed_record,
        [
            "component_result.signed_compact_graph.endpoint_difference_x2_minus_x1",
            "signed_compact_graph.endpoint_difference_x2_minus_x1",
        ],
        _MISSING,
    )
    if endpoint is _MISSING or not isinstance(endpoint, dict):
        raise ValueError(
            "signed endpoint_difference_x2_minus_x1 must be an actual C1 mapping"
        )

    x = graph.local("spb_local_common_x")
    x_jet = _c1(x, graph.const(0))
    native_geometry = _chart_geometry(
        graph,
        "native",
        native_rc,
        native_delta,
        native_z,
        native_coordinate,
        x_jet,
    )
    old_geometry = _chart_geometry(
        graph,
        "old",
        old_rc,
        old_delta,
        old_z,
        old_coordinate,
        x_jet,
    )

    native_fields = {
        "U": _source_c2(graph, native_u_rows, "native", "physical_U"),
        "Q": _source_c2(graph, native_q_rows, "native", "physical_Q"),
    }
    pstar_c2 = _c2(
        graph.source(pstar, "native", "Pstar", 0),
        graph.const(0),
        graph.const(0),
    )
    bg_v_norm = _source_c2(graph, bg_v_rows, "native_BG", "normalized_V")
    bg_q_norm = _source_c2(graph, bg_q_rows, "native_BG", "normalized_Q")
    bg_fields = {
        "U": _scaled_c2(graph, pstar_c2, bg_v_norm),
        "Q": _scaled_c2(graph, pstar_c2, bg_q_norm),
    }
    old_fields = {
        "U": _source_c2(graph, old_u_rows, "old", "physical_U"),
        "Q": _source_c2(graph, old_q_rows, "old", "physical_Q"),
    }

    native_v = _chart_v(graph, native_geometry, native_fields)
    bg_v = _chart_v(graph, native_geometry, bg_fields)
    old_v = _chart_v(graph, old_geometry, old_fields)

    six_names = [
        "q_linear",
        "z_qz_linear",
        "q_square",
        "u_qz",
        "u_square",
        "u_uz",
    ]
    six_difference = {}
    for name in six_names:
        six_difference[name] = _as_c1_dict(
            _c1_sub(graph, native_v["terms"][name], bg_v["terms"][name])
        )

    eU = _c2_sub(graph, native_fields["U"], bg_fields["U"])
    eQ = _c2_sub(graph, native_fields["Q"], bg_fields["Q"])

    vN = native_v["v"]
    vBG = bg_v["v"]
    vO = old_v["v"]
    WN = native_geometry["weight"]
    WO = old_geometry["weight"]
    atom_WNvN = _c1_mul(graph, WN, vN)
    atom_WNvBG = _c1_mul(graph, WN, vBG)
    atom_WOvO = _c1_mul(graph, WO, vO)
    weighted_vol_NBG = _c1_sub(graph, atom_WNvN, atom_WNvBG)
    weighted_delta_v = _c1_sub(graph, atom_WNvN, atom_WOvO)
    reference_offset = _c1_sub(graph, atom_WNvBG, atom_WOvO)
    geometry_offset = _c1_mul(
        graph,
        _c1_sub(graph, WN, WO),
        vBG,
    )
    field_offset = _c1_mul(graph, WO, _c1_sub(graph, vBG, vO))
    sparse_atoms = {
        "WNvN": atom_WNvN,
        "WNvBG": atom_WNvBG,
        "WOvO": atom_WOvO,
    }
    evol_input_coefficients = [
        Fraction(1),
        Fraction(-1),
        Fraction(0),
    ]
    offset_input_coefficients = [
        Fraction(0),
        Fraction(1),
        Fraction(-1),
    ]
    total_coefficients = [
        evol_coefficient + offset_coefficient
        for evol_coefficient, offset_coefficient in zip(
            evol_input_coefficients,
            offset_input_coefficients,
        )
    ]
    response_density_sparse, response_atom_terms = _sparse_c1(
        graph,
        sparse_atoms,
        total_coefficients,
    )
    Evol = _formal_integral(graph, weighted_vol_NBG, "Evol_NBG")
    offset = _formal_integral(graph, reference_offset, "offset_NBG_to_O")
    cached_endpoint, endpoint_metadata = _cached_c1(graph, endpoint, signed_record)
    response = _c1_add(graph, _c1_add(graph, Evol, offset), cached_endpoint)

    source_default_flags = {
        "native_U": native_u_default,
        "native_Q": native_q_default,
        "background_V": bg_v_default,
        "background_Q": bg_q_default,
        "native_Rc": native_rc_default,
        "native_delta": native_delta_default,
        "native_Z": native_z_default,
        "native_coordinate": native_coordinate_default,
        "Pstar": pstar_default,
        "old_U": old_u_default,
        "old_Q": old_q_default,
        "old_Rc": old_rc_default,
        "old_delta": old_delta_default,
        "old_Z": old_z_default,
        "old_coordinate": old_coordinate_default,
    }
    default_rows = [
        name for name, used_default in source_default_flags.items() if used_default
    ]
    if default_rows:
        raise ValueError(
            "default source contract rows are forbidden: " + ", ".join(default_rows)
        )

    native_archive_raw, native_archive_path = _pick(
        native_record,
        [
            "native_view.native_archive.raw_sha256",
            "native_archive.raw_sha256",
            "archive.raw_sha256",
            "raw_sha256",
        ],
        _NATIVE_ARCHIVE_RAW_CONTRACT,
    )
    old_archive_raw, old_archive_path = _pick(
        old_source,
        ["raw_sha256", "archive_raw_sha256"],
        _OLD_ARCHIVE_RAW_CONTRACT,
    )
    native_archive_name, native_archive_name_path = _pick(
        native_record,
        [
            "native_view.native_archive.path",
            "native_archive.path",
            "archive.path",
            "archive_path",
        ],
        _NATIVE_ARCHIVE_PATH_CONTRACT,
    )
    old_archive_name, old_archive_name_path = _pick(
        old_source,
        ["archive", "archive_path"],
        _OLD_ARCHIVE_PATH_CONTRACT,
    )

    def jet(value):
        return {"value": int(value["value"]), "dz": int(value["dz"])}

    def c2jet(value):
        return {
            "value": int(value["value"]),
            "dz": int(value["dz"]),
            "dzz": int(value["dzz"]),
        }

    source_refs = {
        "native": {
            "U_C2": native_u_rows,
            "Q_C2": native_q_rows,
            "Rc": native_rc,
            "delta": native_delta,
            "Z": native_z,
            "coordinate": native_coordinate,
            "Pstar": pstar,
            "source_path": native_u_path or native_q_path or "contract_defaults",
        },
        "native_BG": {
            "normalized_V_C2": bg_v_rows,
            "normalized_Q_C2": bg_q_rows,
            "physical_lift_scale": pstar,
            "additional_N_factor": False,
            "physical_lift": [
                {"ordinary_Z_order": order, "scale_node": pstar, "normalized_V_row": row}
                for order, row in enumerate(bg_v_rows)
            ],
            "physical_Q_lift": [
                {"ordinary_Z_order": order, "scale_node": pstar, "normalized_Q_row": row}
                for order, row in enumerate(bg_q_rows)
            ],
            "source_path": bg_v_path or bg_q_path or "contract_defaults",
        },
        "old": {
            "U_C2": old_u_rows,
            "Q_C2": old_q_rows,
            "Rc": old_rc,
            "delta": old_delta,
            "Z": old_z,
            "coordinate": old_coordinate,
            "source_path": old_source_path or old_u_path or old_q_path or "contract_defaults",
        },
    }

    counts = {
        "integer_dag_nodes": len(graph.nodes),
        "source_ref_nodes": graph.counters["source_ref"],
        "cached_ref_nodes": graph.counters["cached_ref"],
        "ordinary_Z_orders": [0, 1],
        "source_field_jet_orders": [0, 1, 2],
        "c1_product_rules": graph.counters["c1_product_rules"],
        "c1_quotient_rules": graph.counters["c1_quotient_rules"],
        "c2_product_rules": graph.counters["c2_product_rules"],
        "mul_nodes_total": graph.counters["mul"],
        "div_nodes_total": graph.counters["div"],
        "formal_integral_nodes": graph.counters["formal_integral"],
        "sparse_linear_nodes": graph.counters["sparse_linear"],
        "new_normalizations": 0,
        "original_source_evaluations": 0,
        "numeric_evaluations": 0,
        "quadrature_calls": 0,
        "proof_calls": 0,
        "old_checker_calls": 0,
    }

    return {
        "schema": "compact_volume_error_graph_v1",
        "source_only": True,
        "graph": {
            "nodes": graph.nodes,
            "node_count": len(graph.nodes),
            "initial_constant_values_are_exact": True,
            "integer_node_ids": True,
            "no_source_graph_replay": True,
        },
        "source_refs": source_refs,
        "source_identity": {
            "native_archive_raw_sha256": str(native_archive_raw),
            "old_archive_raw_sha256": str(old_archive_raw),
            "native_archive_path": str(native_archive_name),
            "old_archive_path": str(old_archive_name),
            "native_source_namespace": "true_compact_C5_bindings",
            "native_background_namespace": "true_compact_C5_bindings.background",
            "old_source_namespace": "signed_record.component_result.source_interface.old",
            "native_and_old_coordinate_handles_not_used_as_common_x": True,
            "common_x": "spb_local_common_x",
            "common_x_binding_node": x,
            "common_x_derivative_node": x_jet["dz"],
        },
        "pointwise": {
            "vN": jet(vN),
            "vBG": jet(vBG),
            "vO": jet(vO),
            "native_geometry_reused_for_N_BG": True,
            "old_geometry_namespace_separate": True,
            "native_fields": {
                "U": c2jet(native_fields["U"]),
                "Q": c2jet(native_fields["Q"]),
            },
            "background_physical_fields": {
                "U": c2jet(bg_fields["U"]),
                "Q": c2jet(bg_fields["Q"]),
            },
            "old_fields": {
                "U": c2jet(old_fields["U"]),
                "Q": c2jet(old_fields["Q"]),
            },
        },
        "pointwise_vN": jet(vN),
        "pointwise_vBG": jet(vBG),
        "pointwise_vO": jet(vO),
        "geometry": {
            "native": {
                "Z": jet(native_geometry["Z"]),
                "delta": jet(native_geometry["delta"]),
                "Rc": jet(native_geometry["Rc"]),
                "d": jet(native_geometry["d"]),
                "h": jet(native_geometry["h"]),
                "L": jet(native_geometry["L"]),
                "L2": jet(native_geometry["L2"]),
                "coefficients": {
                    name: jet(value)
                    for name, value in native_geometry["coefficients"].items()
                },
                "radius": jet(native_geometry["radius"]),
                "mu": jet(native_geometry["mu"]),
                "weight": jet(native_geometry["weight"]),
            },
            "native_BG": {
                "shares_native_geometry": True,
                "coefficient_namespace": "native",
            },
            "old": {
                "Z": jet(old_geometry["Z"]),
                "delta": jet(old_geometry["delta"]),
                "Rc": jet(old_geometry["Rc"]),
                "d": jet(old_geometry["d"]),
                "h": jet(old_geometry["h"]),
                "L": jet(old_geometry["L"]),
                "L2": jet(old_geometry["L2"]),
                "coefficients": {
                    name: jet(value)
                    for name, value in old_geometry["coefficients"].items()
                },
                "radius": jet(old_geometry["radius"]),
                "mu": jet(old_geometry["mu"]),
                "weight": jet(old_geometry["weight"]),
            },
        },
        "field_errors": {
            "eU": c2jet(eU),
            "eQ": c2jet(eQ),
            "ordinary_Z_orders": [0, 1, 2],
            "background_physical_lift": "Pstar*normalized jet once at each ordinary-Z order; no extra N",
        },
        "six_correlated_source_difference_channels": {
            "orientation": "actual native N product minus native BG product",
            "ordinary_Z_orders": [0, 1],
            "cached_generic_product_algebra_reused": True,
            "new_proof_or_normalization": False,
            "channels": six_difference,
        },
        "weights": {
            "WN": jet(WN),
            "WO": jet(WO),
            "native_R": jet(native_geometry["radius"]),
            "old_R": jet(old_geometry["radius"]),
            "native_mu": jet(native_geometry["mu"]),
            "old_mu": jet(old_geometry["mu"]),
            "native_and_old_Rc_not_aliased": True,
        },
        "weighted_vol_NBG": jet(weighted_vol_NBG),
        "weighted_Delta_v": jet(weighted_delta_v),
        "reference_offset": jet(reference_offset),
        "geometry_offset": jet(geometry_offset),
        "field_offset": jet(field_offset),
        "pointwise_atoms": {
            "WNvN": jet(atom_WNvN),
            "WNvBG": jet(atom_WNvBG),
            "WOvO": jet(atom_WOvO),
            "assembly_inputs": {
                "Evol_density": [
                    {
                        "atom": name,
                        "coefficient": _fraction_json(coefficient),
                    }
                    for name, coefficient in zip(
                        sparse_atoms,
                        evol_input_coefficients,
                    )
                ],
                "offset_density": [
                    {
                        "atom": name,
                        "coefficient": _fraction_json(coefficient),
                    }
                    for name, coefficient in zip(
                        sparse_atoms,
                        offset_input_coefficients,
                    )
                ],
            },
            "computed_total_coefficients": [
                _fraction_json(coefficient) for coefficient in total_coefficients
            ],
            "computed_total_coefficients_numeric": [
                int(coefficient) for coefficient in total_coefficients
            ],
            "signed_sparse_coefficients": response_atom_terms,
            "signed_sparse_atom_order": ["WNvN", "WNvBG", "WOvO"],
            "signed_sparse_C0_C1_coefficients": [
                _fraction_json(coefficient) for coefficient in total_coefficients
            ],
            "signed_sparse_C0_C1_coefficients_numeric": [
                int(coefficient) for coefficient in total_coefficients
            ],
            "response_density_sparse": jet(response_density_sparse),
            "integration_linearity_only": True,
        },
        "integrals": {
            "Evol": jet(Evol),
            "offset": jet(offset),
            "bounds": [1, 2],
            "measure": "dx",
            "numerical_value": False,
        },
        "Evol": jet(Evol),
        "offset": jet(offset),
        "cached_endpoint": {
            "endpoint_difference_x2_minus_x1": jet(cached_endpoint),
            "source": endpoint_metadata,
            "cached_dag_replayed": False,
            "signed_record_path": endpoint_path,
        },
        "cachedendpoint": jet(cached_endpoint),
        "response": {
            "value": int(response["value"]),
            "dz": int(response["dz"]),
            "formula": "Evol + offset + cached endpoint_difference_x2_minus_x1",
            "ordinary_Z_scope": "value and first derivative",
        },
        "contracts": {
            "native_source_archive_raw_sha256": str(native_archive_raw),
            "old_source_archive_raw_sha256": str(old_archive_raw),
            "native_Pstar_Z_independent": True,
            "native_Rc_Z_independent": True,
            "native_delta_Z_independent": True,
            "old_Rc_Z_independent": True,
            "old_delta_Z_independent": True,
            "L_native_nonzero": "caller/source premise; not proved here",
            "L_old_nonzero": "caller/source premise; not proved here",
            "ordinary_Z_derivative_is_not_temporal": True,
            "cap_consumed": False,
            "midpoint_consumed": False,
            "default_source_contract_rows_used": source_default_flags,
            "source_paths": {
                "native_U": native_u_path,
                "native_Q": native_q_path,
                "background_V": bg_v_path,
                "background_Q": bg_q_path,
                "old_source_interface": old_source_path,
                "cached_endpoint": endpoint_path,
            },
        },
        "counts": counts,
        "nonclaims": [
            "No cap or midpoint enters the graph.",
            "No native equals old or native N equals BG zero claim is made.",
            "No numerical source value, quadrature, proof, old checker, or normalization is run.",
            "The cached endpoint is referenced as an external C1 node and its DAG is not replayed.",
            "The generic six-channel product algebra is consumed as an accepted contract and not re-proved.",
        ],
    }
# Root-owned source contracts and archived cap provenance. No ancestor code is
# imported, evaluated, normalized, or used to choose a signed function value.
_BASE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_'
_NATIVE = _BASE + 'true_compact_C5_functions.json.gz'
_SIGNED = _BASE + 'signed_pullback_difference.json'
_MACRO = _BASE + 'macro_endpoint_m15_numeric.json'
_QUIET = _BASE + 'compact_quiet_m15_numeric.json'
_NATIVE_RAW = '3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269'
_OLD_RAW = '3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce'
_NORMALIZATION = _BASE + 'modified_Q4_pressure_C3_functions.py'
_NORMALIZATION_RAW = '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7'


def _need(condition, message):
    if not condition:
        raise ValueError(message)


def _source_contract(source):
    records = source['records']
    native, signed, macro, quiet = (records[x] for x in [_NATIVE, _SIGNED, _MACRO, _QUIET])
    parent = source['parent_report']
    _need(parent['accepted_identity_flag'] is True and parent['source_ftc_transport_defined'] is True,
          'Cached conditional FTC source theorem is not accepted')
    b, nodes = native['true_compact_C5_bindings'], native['exact_graph_nodes']
    interface = signed['component_result']['source_interface']
    family = native['source_family']
    _need(family == signed['component_result']['source_family'] == macro['source_family'] == quiet['source_family'],
          'Actual construction family fingerprints differ')
    _need(source['record_metadata'][_NATIVE]['raw_sha256'] == _NATIVE_RAW, 'Native source record changed')
    _need(source['record_metadata'][_NORMALIZATION]['raw_sha256'] == _NORMALIZATION_RAW,
          'Physical/generic normalization producer changed')
    _need(b['physical']['U'][:3] == [262269, 262270, 262271] and b['physical']['Q'][:3] == [262364, 262368, 262374],
          'Actual physical compact ports changed')
    _need(b['background']['generic']['Q'][:3] == [260923, 260927, 260933], 'Normalized background Q ports changed')
    _need(b['background']['fields']['V'][:3] == [0, 0, 0] and nodes[0] == {'operation': 'exact_rational', 'numerator': 0, 'denominator': 1},
          'Background normalized V zero source changed')
    _need(nodes[260862]['operation'] == 'product' and nodes[260862]['arguments'] == [251066, 260181],
          'Native physical radius pullback changed')
    _need(b['Rc'] == 251066 and b['delta'] == 25220 and b['Pstar'] == 33572 and b['background']['phase'] == 260181
          and interface['native']['Z_node'] == 20, 'Native geometry source contract changed')
    _need(interface['native']['U'] == [262269, 262270, 262271] and interface['native']['Q'] == [262364, 262368, 262374],
          'Native graph / cached FTC field namespace differs')
    for k, row_id in enumerate([262269, 262270, 262271]):
        row = nodes[row_id]
        _need(row['operation'] == 'product' and row['arguments'] == [[260375, 260378, 260382][k], 33572],
              'Actual U physical scale no longer appears exactly once')
    for row_id in [260923, 260927, 260933]:
        _need(nodes[row_id]['operation'] == 'positive_quotient' and nodes[row_id]['denominator'] == 25272,
              'Background normalized Q quotient ancestry changed')
    n = nodes[b['N']]
    _need(n['type'] == 'DoubleDyadicInteger' and n['J'] == 1358356628656378313 and n['Z_independent'] is True,
          'Native exact selected integer changed')
    for name, digest in [('true_compact_C5_functions.json.gz', _NATIVE_RAW), ('n1_preRp_moment_functions.json.gz', _OLD_RAW)]:
        _need(quiet['input_hashes'][_BASE.split('/')[-1] + name] == digest,
              'Retained compact cap does not directly pin the selected source archive')
    return {'passed': True, 'same_family_fingerprint': True, 'family': family,
            'selected_integer': n, 'selected_integer_materialized': False,
            'native_actual_U_C2': b['physical']['U'][:3], 'native_actual_Q_C2': b['physical']['Q'][:3],
            'background_normalized_V_C2': b['background']['fields']['V'][:3],
            'background_normalized_Q_C2': b['background']['generic']['Q'][:3],
            'background_normalized_m_C3': b['background']['fields']['m'][:4],
            'background_physical_lift': 'U_b=Pstar*V_b; Q_b=Pstar*q_b, once, Pstar_Z=0',
            'background_Q_is_not_already_physical': True, 'additional_N_factor': False,
            'normalization_producer_raw_sha256': _NORMALIZATION_RAW,
            'measure': 'Rc_N^2*x dx on x in [1,2]; old Rc_O retained independently',
            'native_numeric_source_archive_direct_pin': True,
            'physical_lift_basis': 'generic_functions recovers q=N(V,m)/L; physical_pressure_fields scales V,m by Pstar before the same linear Q recovery',
            'old_cap_to_signed_function_identity_proved': False}


def _cap_provenance(records, n):
    macro, quiet = records[_MACRO], records[_QUIET]
    snapshots = []
    for sample in macro['samples']:
        z = sample['Z']
        matches = [x for x in quiet['samples'] if x['Z'] == z and x['chart'] == 'compact' and str(x['native_coordinate']) == '2']
        _need(len(matches) == 1, 'Unique compact upper cap sample missing')
        compact = matches[0]
        raw, reduced = sample['raw_m15_C0_C1_terms'], sample['correlated_reduced_m15_C0_C1_terms']
        _need(len(raw) == 213 and len(reduced) == 183 and len(compact['m15_C0_C1_terms']) == 211,
              'Raw/reduced/compact packet cardinalities changed')
        term = raw[209]
        _need(term['name'] == 'compact_finite_N_volume_error' and term == reduced[180] == compact['m15_C0_C1_terms'][209],
              'Archived volume cap row does not have an exact three-way provenance match')
        _need(sample['selected_N_evidence'] == n == compact['actual_selected_integer'], 'Selected integer differs across packets')
        binding = term['actual_function_binding']
        _need(binding['chart'] == 'compact' and binding['incoming_terms_retained_once'] == 207 and binding['error_enclosure_only'] is True,
              'Volume cap was relabeled as an evaluated field')
        _need(all(x.get('source_value_selected') is False and x.get('positive_error_bound_atom_only') is True for x in term['log_scale_terms']),
              'Volume cap log atom no longer carries the error-only contract')
        _need(raw[212]['name'] == 'reserved_n1_own_swirl_once', 'Reserved contribution is not registered exactly once at the retained slot')
        snapshots.append({'Z_sample': z, 'ordinary_Z_orders': [0, 1], 'raw_selector': '$.samples[Z=' + z + '].raw_m15_C0_C1_terms[209]',
                          'reduced_selector': '$.samples[Z=' + z + '].correlated_reduced_m15_C0_C1_terms[180]',
                          'raw_equals_reduced_equals_compact': True, 'cap_row_exact_copy': term,
                          'directed_volume_logs_exact_copy': compact['finite_N_error_logs_in_physical_units']['volume'],
                          'signed_function_value_selected': False})
    _need(len(snapshots) == 3, 'Expected three archived cap samples')
    return {'passed': True, 'samples': snapshots, 'selected_N_metadata_matched': True,
            'same_family_is_not_a_function_identity': True, 'numeric_centers_used_as_source_values': False,
            'cap_to_source_bridge_proved': False, 'all_original_caps_retained': True,
            'unresolved_bridge': 'Bind the numeric baseline mean and complete-prefix error definitions to the native background m_C3 and the newly lifted physical Q_b; row/N/archive provenance alone is insufficient.'}


def derive(source, journal):
    """Define the actual native-minus-background signed source, preserving caps."""
    result = {'schema': 'original_compact_volume_error_source_packet', 'component': 'compact_volume_error_source_packet',
              'passed': False, 'metadata_only': False, 'source_only': True, 'source_packet_definition_only': True,
              'input_digest': source['input_digest'], 'original_source_evaluations': 0, 'runtime_model_metadata': None,
              'source_record_hashes': {name: meta['raw_sha256'] for name, meta in source['record_metadata'].items()},
              'source_ftc_transport_defined': source['parent_report']['source_ftc_transport_defined'],
              'source_packet_defined': False, 'cap_to_source_identity_proved': False,
              'fresh_source_normalizations': 0, 'cached_formal_proof_calls': 0,
              'retained_construction_data': source['parent_report']['retained_construction_data'],
              'retained_integral_term_inventory': source['parent_report']['retained_integral_term_inventory'],
              'source_inspection_audit': {'old_source_consumer_calls': 0, 'old_proof_functions_called': 0,
                                          'old_operator_calls': 0, 'old_source_normalizer_calls': 0,
                                          'original_source_evaluations': 0, 'quadratures': 0,
                                          'temporal_recurrence_calls': 0, 'numeric_field_evaluations': 0,
                                          'new_source_expression_builder_calls': 0}}
    result.update({x: False for x in source['parent_report']['false_flags']})
    try:
        contract = _source_contract(source)
        result['source_contract'] = contract
        journal({'event': 'actual_background_port_contract_bound', 'normalization': 'BG generic Q lifted by Pstar once'})
        cap = _cap_provenance(source['records'], contract['selected_integer'])
        result['retained_cap_provenance'] = cap
        journal({'event': 'three_archived_cap_rows_joined', 'raw_row': 209, 'reduced_row': 180, 'selected_signed_value': False})
        native_view = dict(source['records'][_NATIVE])
        native_view['geometry'] = {'Rc': {'node': 251066}, 'delta': {'node': 25220}, 'Pstar': {'node': 33572},
                                   'Z': {'node': 20}, 'local_coordinate': {'node': 260181}}
        native_view['native_archive'] = {'path': _NATIVE, 'raw_sha256': _NATIVE_RAW}
        result['source_inspection_audit']['new_source_expression_builder_calls'] = 1
        graph = build_volume_source_graph(native_view, source['records'][_SIGNED])
        result['signed_volume_source_graph'] = graph
        result['source_packet_defined'] = True
        result['accepted_scope'] = 'Source-bound C1 signed function definition and exact archived cap provenance; cap-to-source inequality bridge remains open'
        result['formal_preconditions'] = {'ordinary_Z': 'Z_Z=1; delta_Z=Rc_Z=Pstar_Z=x_Z=0',
                                           'L_native_and_old_nonzero_required': True,
                                           'common_source_domain_required': True, 'physical_domain_admission': False,
                                           'cached_FTC_radial_jet_and_mixed_derivative_contracts_retained': True}
        normalized_upper_mass = (Fraction(1) - Fraction(1, 4)) / 2
        normalized_lower_mass = (Fraction(4) - Fraction(1)) / 2
        _need(Fraction(4) * normalized_upper_mass == normalized_lower_mass, 'Physical compact mass normalization failed')
        result['measure_normalization'] = {'upper_radius_R2_mass': str(normalized_upper_mass), 'lower_Rc2_mass': str(normalized_lower_mass),
                                            'exact_rational_check': '4*(3/8)=3/2', 'exact_rational_check_passed': True,
                                            'new_exact_rational_mass_checks': 1,
                                            'bounds': ['1', '2'], 'no_huge_radius_materialized': True,
                                            'cap_inequality_not_recertified': True}
        result['remaining_bridges'] = [cap['unresolved_bridge'],
                                        'Compute the signed native-background volume source and the distinct background-to-old reference offset with certified numerical errors.',
                                        'Complete the ledger-wide source correspondence before replacing any inherited quantitative cap.']
        result['passed'] = True
        journal({'event': 'signed_source_packet_defined', 'source_evaluations': 0, 'cap_to_source_identity': False})
    except Exception as error:
        import traceback
        result['passed'] = False
        result['failure'] = {'type': type(error).__name__, 'message': str(error), 'traceback': traceback.format_exc()}
        journal({'event': 'new_source_packet_failure', 'type': type(error).__name__, 'message': str(error)})
    return result

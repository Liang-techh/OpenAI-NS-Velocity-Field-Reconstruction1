"""Fresh cached-C2 volume graph construction; intentionally unentered."""
from fractions import Fraction

_NATIVE_ARCHIVE = "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"
_OLD_ARCHIVE = "3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce"
_MISSING = object()

def _f(x):
    if isinstance(x, Fraction): return x
    if isinstance(x, int): return Fraction(x, 1)
    if isinstance(x, str): return Fraction(x)
    if isinstance(x, dict): return Fraction(int(x["numerator"]), int(x["denominator"]))
    raise TypeError("exact rational required")

def _fj(x):
    x = _f(x); return {"numerator": x.numerator, "denominator": x.denominator}

def _get(obj, path):
    for part in path.split("."):
        if not isinstance(obj, dict) or part not in obj: return _MISSING
        obj = obj[part]
    return obj

def _need(obj, *paths):
    for path in paths:
        value = _get(obj, path)
        if value is not _MISSING and value is not None: return value
    raise KeyError("missing explicit input: " + ", ".join(paths))

def _node(x, label):
    if isinstance(x, dict):
        for key in ("node", "row", "source_row", "value"):
            if key in x: x = x[key]; break
    if isinstance(x, bool): raise TypeError(label + " cannot be boolean")
    return int(x)

def _rows(x, label):
    if isinstance(x, dict) and "rows" in x: x = x["rows"]
    if not isinstance(x, (list, tuple)) or len(x) != 3:
        raise ValueError(label + " requires exactly three C2 rows")
    return [_node(v, label) for v in x]

class _DAG:
    def __init__(self):
        self.nodes, self.index = [], {}
        self.counts = {"add": 0, "sub": 0, "mul": 0, "div": 0,
                       "neg": 0, "source_ref": 0, "cached_ref": 0,
                       "bound_variable": 0, "formal_integral": 0,
                       "sparse_linear": 0, "c1_product_rules": 0,
                       "c1_quotient_rules": 0, "c2_product_rules": 0}
    def _put(self, key, rec):
        if key in self.index: return self.index[key]
        i = len(self.nodes); rec = dict(rec); rec["id"] = i
        self.nodes.append(rec); self.index[key] = i
        if rec.get("op") in self.counts: self.counts[rec["op"]] += 1
        return i
    def const(self, x):
        x = _f(x); return self._put(("const", x.numerator, x.denominator),
            {"op": "const", "value": _fj(x)})
    def source(self, archive, chart, node, label, order):
        node, order = int(node), int(order)
        return self._put(("source_ref", archive, chart, node, label, order), {
            "op": "source_ref", "archive_raw_sha256": str(archive),
            "chart": chart, "source_node": node, "label": label,
            "ordinary_Z_order": order})
    def cached(self, node, label, order):
        node, order = int(node), int(order)
        return self._put(("cached_ref", node, label, order), {
            "op": "cached_ref", "cached_node": node, "label": label,
            "ordinary_Z_order": order, "replay": False})
    def local(self, name):
        return self._put(("bound_variable", name), {
            "op": "bound_variable", "name": name, "ordinary_Z_derivative": 0})
    def op(self, op, args, **meta):
        args = tuple(int(a) for a in args)
        key = (op, args, tuple(sorted((k, repr(v)) for k, v in meta.items())))
        rec = {"op": op, "args": list(args)}; rec.update(meta)
        return self._put(key, rec)

def _c1(v, dz): return {"value": int(v), "dz": int(dz)}
def _c2(v, dz, dzz): return {"value": int(v), "dz": int(dz), "dzz": int(dzz)}
def _c1k(g, x): return _c1(g.const(x), g.const(0))
def _c2k(g, x): return _c2(g.const(x), g.const(0), g.const(0))
def _a(g, x, y): return _c1(g.op("add", [x["value"], y["value"]]), g.op("add", [x["dz"], y["dz"]]))
def _s(g, x, y): return _c1(g.op("sub", [x["value"], y["value"]]), g.op("sub", [x["dz"], y["dz"]]))
def _n(g, x): return _c1(g.op("neg", [x["value"]]), g.op("neg", [x["dz"]]))
def _m(g, x, y):
    g.counts["c1_product_rules"] += 1
    return _c1(g.op("mul", [x["value"], y["value"]]), g.op("add", [g.op("mul", [x["dz"], y["value"]]), g.op("mul", [x["value"], y["dz"]])]))
def _d(g, x, y):
    g.counts["c1_quotient_rules"] += 1
    den = g.op("mul", [y["value"], y["value"]])
    top = g.op("sub", [g.op("mul", [x["dz"], y["value"]]), g.op("mul", [x["value"], y["dz"]])])
    return _c1(g.op("div", [x["value"], y["value"]]), g.op("div", [top, den]))
def _sum(g, xs):
    out = _c1k(g, 0)
    for x in xs: out = _a(g, out, x)
    return out
def _c2m(g, x, y):
    g.counts["c2_product_rules"] += 1
    v = g.op("mul", [x["value"], y["value"]])
    dz = g.op("add", [g.op("mul", [x["dz"], y["value"]]), g.op("mul", [x["value"], y["dz"]])])
    cross = g.op("mul", [g.const(2), g.op("mul", [x["dz"], y["dz"]])])
    zz = g.op("add", [g.op("mul", [x["dzz"], y["value"]]), g.op("add", [cross, g.op("mul", [x["value"], y["dzz"]])])])
    return _c2(v, dz, zz)
def _src_c2(g, rows, archive, chart, label):
    return _c2(*(g.source(archive, chart, row, label, i) for i, row in enumerate(rows)))
def _geom(g, archive, chart, geom, x):
    one = _c1k(g, 1); z = _c1(g.source(archive, chart, geom["Z"], "Z", 0), g.const(1))
    delta = _c1(g.source(archive, chart, geom["delta"], "delta", 0), g.const(0))
    rc = _c1(g.source(archive, chart, geom["Rc"], "Rc", 0), g.const(0))
    z2 = _m(g, z, z); d = _s(g, one, z2); h = _s(g, one, delta)
    L = _s(g, one, _m(g, delta, z2)); L2 = _m(g, L, L)
    radius = _m(g, rc, x); mu = _d(g, one, x); W = _m(g, _m(g, radius, radius), mu)
    c = {"q": _d(g, _c1k(g, Fraction(-1, 2)), L),
         "qZ": _d(g, _m(g, h, z), _m(g, _c1k(g, 4), L)),
         "q2": _c1k(g, Fraction(-1, 4)),
         "UQZ": _d(g, d, _m(g, _c1k(g, 2), L)),
         "U2": _d(g, _m(g, h, z2), L2), "UUZ": _d(g, _m(g, z, d), L2)}
    return {"Z": z, "delta": delta, "Rc": rc, "L": L, "d": d, "h": h,
            "weight": W, "coefficients": c,
            "coordinate_source": g.source(archive, chart, geom["local_coordinate"], "archive_coordinate_handle", 0)}
def _v(g, geom, fields):
    U, Q = fields["U"], fields["Q"]; u = _c1(U["value"], U["dz"]); uz = _c1(U["dz"], U["dzz"])
    q = _c1(Q["value"], Q["dz"]); qz = _c1(Q["dz"], Q["dzz"]); c = geom["coefficients"]
    products = {"Q": q, "QZ": qz, "Q2": _m(g, q, q), "UQZ": _m(g, u, qz), "U2": _m(g, u, u), "UUZ": _m(g, u, uz)}
    coefficient_key = {"Q": "q", "QZ": "qZ", "Q2": "q2", "UQZ": "UQZ", "U2": "U2", "UUZ": "UUZ"}
    terms = {name: _m(g, c[coefficient_key[name]], products[name]) for name in products}
    return {"products": products, "terms": terms, "v": _sum(g, list(terms.values()))}
def _sub2(g, x, y): return _c2(g.op("sub", [x["value"], y["value"]]), g.op("sub", [x["dz"], y["dz"]]), g.op("sub", [x["dzz"], y["dzz"]]))
def _integral(g, x, label):
    meta = {"label": label, "lower": 1, "upper": 2, "measure": "dx"}
    value_meta = dict(meta, ordinary_Z_order=0, derivative_under_integral=False)
    derivative_meta = dict(meta, ordinary_Z_order=1, derivative_under_integral=True)
    return _c1(g.op("formal_integral", [x["value"]], **value_meta), g.op("formal_integral", [x["dz"]], **derivative_meta))
def _jet(x): return {"value": int(x["value"]), "dz": int(x["dz"])}
def _jet2(x): return {"value": int(x["value"]), "dz": int(x["dz"]), "dzz": int(x["dzz"])}
def _sparse(g, atoms, coeffs):
    names = list(atoms); coeffs = [_f(x) for x in coeffs]
    meta = {"atom_order": names, "atom_terms": [{"atom": n, "coefficient": _fj(c)} for n, c in zip(names, coeffs)], "all_atom_arguments_retained": True, "exact_merge": True, "integration_linearity_only": True}
    return _c1(g.op("sparse_linear", [atoms[n]["value"] for n in names], **meta), g.op("sparse_linear", [atoms[n]["dz"] for n in names], **meta))

def build_cached_C2_volume_source(saved_contract, signed_component):
    """Build a fresh graph from explicit saved three-row C2 ports only.

    `saved_contract.native_geometry` must explicitly contain Rc, delta, Z,
    local_coordinate, and Pstar node ids. No source-row or higher-order fallback
    is accepted. `signed_component` must contain source_interface.old and the
    saved signed_compact_graph endpoint C1 pair.
    """
    if not isinstance(saved_contract, dict) or not isinstance(signed_component, dict): raise TypeError("explicit mappings required")
    required = ("native_actual_U_C2", "native_actual_Q_C2", "background_normalized_V_C2", "background_normalized_Q_C2")
    rows = {k: _rows(saved_contract[k], k) for k in required}
    ng = _need(saved_contract, "native_geometry")
    geom = {k: _node(_need(ng, k), "native_geometry." + k) for k in ("Rc", "delta", "Z", "local_coordinate", "Pstar")}
    old = _need(signed_component, "source_interface.old")
    old_rows = {"U": _rows(_need(old, "U"), "old U"), "Q": _rows(_need(old, "Q"), "old Q")}
    oldg = {"Rc": _node(_need(old, "Rc_node"), "old Rc"), "delta": _node(_need(old, "delta_node"), "old delta"), "Z": _node(_need(old, "Z_node"), "old Z"), "local_coordinate": _node(_need(old, "coordinate_node"), "old coordinate")}
    endpoint = _need(signed_component, "signed_compact_graph.endpoint_difference_x2_minus_x1")
    ev, ed = _node(_need(endpoint, "value"), "endpoint value"), _node(_need(endpoint, "dz"), "endpoint dz")
    g = _DAG(); x = g.local("cached_C2_common_x"); xj = _c1(x, g.const(0))
    ngm = _geom(g, _NATIVE_ARCHIVE, "native", geom, xj); ogm = _geom(g, _OLD_ARCHIVE, "old", oldg, xj)
    nfields = {"U": _src_c2(g, rows["native_actual_U_C2"], _NATIVE_ARCHIVE, "native", "physical_U"), "Q": _src_c2(g, rows["native_actual_Q_C2"], _NATIVE_ARCHIVE, "native", "physical_Q")}
    pstar = _c2(g.source(_NATIVE_ARCHIVE, "native", geom["Pstar"], "Pstar", 0), g.const(0), g.const(0))
    bfields = {"U": _c2m(g, pstar, _src_c2(g, rows["background_normalized_V_C2"], _NATIVE_ARCHIVE, "native_BG", "normalized_V")), "Q": _c2m(g, pstar, _src_c2(g, rows["background_normalized_Q_C2"], _NATIVE_ARCHIVE, "native_BG", "normalized_Q"))}
    ofields = {"U": _src_c2(g, old_rows["U"], _OLD_ARCHIVE, "old", "physical_U"), "Q": _src_c2(g, old_rows["Q"], _OLD_ARCHIVE, "old", "physical_Q")}
    nv, bv, ov = _v(g, ngm, nfields), _v(g, ngm, bfields), _v(g, ogm, ofields)
    WN, WO = ngm["weight"], ogm["weight"]; atoms = {"WNvN": _m(g, WN, nv["v"]), "WNvBG": _m(g, WN, bv["v"]), "WOvO": _m(g, WO, ov["v"])}
    NBG = _s(g, atoms["WNvN"], atoms["WNvBG"]); offset = _s(g, atoms["WNvBG"], atoms["WOvO"])
    weighted_delta = _s(g, atoms["WNvN"], atoms["WOvO"])
    Evol, Off = _integral(g, NBG, "Evol_NBG"), _integral(g, offset, "reference_offset")
    endpoint_c1 = _c1(g.cached(ev, "endpoint_difference_x2_minus_x1", 0), g.cached(ed, "endpoint_difference_x2_minus_x1", 1))
    response = _a(g, _a(g, Evol, Off), endpoint_c1)
    six = {n: _s(g, nv["terms"][n], bv["terms"][n]) for n in nv["terms"]}
    evol_coefficients = [_f(1), _f(-1), _f(0)]
    offset_coefficients = [_f(0), _f(1), _f(-1)]
    total_coefficients = [a + b for a, b in zip(evol_coefficients, offset_coefficients)]
    sparse = _sparse(g, atoms, total_coefficients)
    eU = _sub2(g, nfields["U"], bfields["U"])
    eQ = _sub2(g, nfields["Q"], bfields["Q"])
    weighted_six = {n: _m(g, WN, six[n]) for n in six}
    geometry_offset = _m(g, _s(g, WN, WO), bv["v"])
    field_offset = _m(g, WO, _s(g, bv["v"], ov["v"]))
    graph_count = len(g.nodes)
    count_snapshot = dict(g.counts, new_builder_calls=1, original_source_evaluations=0,
                          numeric_evaluations=0, quadrature_calls=0,
                          normalization_calls=0, old_proof_calls=0)
    return {
        "schema": "compact_cached_C2_volume_graph_v1",
        "source_only": True,
        "graph_constructed": True,
        "graph": {"nodes": g.nodes, "node_count": graph_count},
        "source_refs": {
            "native_archive_raw_sha256": _NATIVE_ARCHIVE,
            "old_archive_raw_sha256": _OLD_ARCHIVE,
            "native_rows": rows,
            "old_rows": old_rows,
            "native_geometry": geom,
            "old_geometry": oldg,
        },
        "pointwise": {
            "vN": _jet(nv["v"]),
            "vBG": _jet(bv["v"]),
            "vO": _jet(ov["v"]),
            "eU": _jet2(eU),
            "eQ": _jet2(eQ),
            "six_weighted_N_minus_BG_term_differences": {
                n: _jet(weighted_six[n]) for n in weighted_six
            },
        },
        "weighted_vol_NBG": _jet(NBG),
        "weighted_Delta_v": _jet(weighted_delta),
        "reference_offset": _jet(offset),
        "geometry_offset": _jet(geometry_offset),
        "field_offset": _jet(field_offset),
        "pointwise_atoms": {
            "atoms": {n: _jet(v) for n, v in atoms.items()},
            "Evol_coefficients": [_fj(x) for x in evol_coefficients],
            "offset_coefficients": [_fj(x) for x in offset_coefficients],
            "computed_total_coefficients": [_fj(x) for x in total_coefficients],
            "response_density_sparse": _jet(sparse),
        },
        "integrals": {
            "Evol": _jet(Evol),
            "reference_offset": _jet(Off),
            "bounds": [1, 2],
            "measure": "dx",
        },
        "cached_endpoint": {"value": ev, "dz": ed, "replay": False},
        "response": _jet(response),
        "counts": count_snapshot,
    }
def derive(source, journal):
    """Consume cached C2 ports once; archived enclosure data never enters DAG."""
    saved = source['parent_result']
    result = {
        'schema': 'original_compact_cached_C2_volume_source', 'component': 'compact_cached_C2_volume_source',
        'passed': False, 'source_only': True, 'metadata_only': False, 'source_packet_definition_only': True,
        'source_packet_defined': False, 'input_digest': source['input_digest'], 'runtime_model_metadata': None,
        'source_ftc_transport_defined': True, 'cap_to_source_identity_proved': False,
        'original_source_evaluations': 0, 'cached_formal_proof_calls': 0, 'fresh_source_normalizations': 0,
        'retained_construction_data': saved['retained_construction_data'],
        'retained_integral_term_inventory': saved['retained_integral_term_inventory'],
        'source_inspection_audit': {'new_C2_expression_builder_calls': 0, 'old_source_consumer_calls': 0,
                                    'old_proof_functions_called': 0, 'old_operator_calls': 0,
                                    'old_source_normalizer_calls': 0, 'quadratures': 0,
                                    'numeric_field_evaluations': 0, 'temporal_recurrence_calls': 0},
    }
    result.update({key: False for key in source['parent_report']['false_flags']})
    try:
        contract = saved['source_contract']
        cap = saved['retained_cap_provenance']
        if contract['passed'] is not True or cap['passed'] is not True or saved['source_ftc_transport_defined'] is not True:
            raise ValueError('Authenticated cached partial source/cap contracts are not available')
        # Only metadata authentication and C2 interface binding are repeated;
        # no original source theorem, 100 MB ledger join, or old builder runs.
        native = source['signed_component']['source_interface']['native']
        expected = {
            'native_actual_U_C2': [262269, 262270, 262271],
            'native_actual_Q_C2': [262364, 262368, 262374],
            'background_normalized_V_C2': [0, 0, 0],
            'background_normalized_Q_C2': [260923, 260927, 260933],
        }
        for name, rows in expected.items():
            if contract[name] != rows or len(contract[name]) != 3:
                raise ValueError('Explicit cached C2 port binding changed: ' + name)
        if native['U'] != expected['native_actual_U_C2'] or native['Q'] != expected['native_actual_Q_C2']:
            raise ValueError('Cached C2 fields differ from the accepted signed source interface')
        view = dict(contract)
        view['native_geometry'] = {'Rc': native['Rc_node'], 'delta': native['delta_node'], 'Z': native['Z_node'],
                                   'local_coordinate': native['coordinate_node'], 'Pstar': 33572}
        if view['native_geometry'] != {'Rc': 251066, 'delta': 25220, 'Z': 20, 'local_coordinate': 260181, 'Pstar': 33572}:
            raise ValueError('Explicit cached native geometry changed')
        result['cached_source_contract'] = contract
        result['retained_cap_provenance'] = cap
        result['cached_provenance_reuse'] = {'source_contract_passed': True, 'cap_provenance_passed': True,
                                             'cap_rows_recomputed': 0, 'selected_N_recomputed': 0,
                                             'cached_proof_functions_called': 0,
                                             'parent_failure_result_sha256': '809ae44b56ee4428da44708630a72a38f0538102671b9ab05d5f404ac4cffb5d'}
        journal({'event': 'cached_C2_interface_bound', 'higher_order_archive_arrays_consumed': False})
        result['source_inspection_audit']['new_C2_expression_builder_calls'] = 1
        graph = build_cached_C2_volume_source(view, source['signed_component'])
        result['signed_volume_source_graph'] = graph
        mass_upper = (Fraction(1) - Fraction(1, 4)) / 2
        mass_lower = (Fraction(4) - Fraction(1)) / 2
        if 4 * mass_upper != mass_lower:
            raise ValueError('Compact physical measure normalization failed')
        result['measure_normalization'] = {'upper_R2_mass': str(mass_upper), 'lower_Rc2_mass': str(mass_lower),
                                           'exact_mass_identity_passed': True, 'new_exact_rational_checks': 1,
                                           'physical_radius_materialized': False, 'cap_bound_recertified': False}
        result['formal_preconditions'] = {'ordinary_Z': 'Z_Z=1; delta_Z=Rc_Z=Pstar_Z=x_Z=0',
                                           'L_native_and_old_nonzero_required': True,
                                           'source_common_domain_required': True,
                                           'cached_FTC_radial_jet_and_mixed_derivative_contracts_required': True,
                                           'physical_domain_admission': False}
        result['source_packet_defined'] = True
        result['passed'] = True
        result['accepted_scope'] = 'C1 signed source definition on explicit cached C2 ports; numerical cap-to-source bridge remains open'
        result['remaining_bridge'] = cap['unresolved_bridge']
        journal({'event': 'new_cached_C2_graph_constructed', 'node_count': graph['graph']['node_count'],
                 'source_value_evaluations': 0, 'cap_identity_proved': False})
    except Exception as error:
        import traceback
        result['passed'] = False
        result['failure'] = {'type': type(error).__name__, 'message': str(error), 'traceback': traceback.format_exc()}
        journal({'event': 'cached_C2_source_failure', 'type': type(error).__name__, 'message': str(error)})
    return result

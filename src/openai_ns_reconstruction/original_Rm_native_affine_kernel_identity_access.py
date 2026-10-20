"""Source-only affine identities for the native C2/C5 transport kernels.

This module consumes an already authenticated graph and source descriptors.  It
does not import a producer, execute a graph operator, evaluate a numerical
parameter, or call an old checker.  The proof is deliberately limited to the
real-coordinate affine forms: a C2 offset is differentiated in its graph, a
C5 radius exponent is differentiated in its graph, and the two exact slopes
are compared after an explicit typed coordinate pullback.  The reference
chart uses the separately authenticated source projection from the accepted
cross-projection receipt; no interval or name equality is used as a source
identity.
"""

from __future__ import annotations

import ast
import base64
import binascii
from collections.abc import Mapping, Sequence
from fractions import Fraction
import gzip
import hashlib
import json
import textwrap
from typing import Any


SCHEMA = "original_Rm_native_affine_kernel_identity_access_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
PARENT_REPORT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
POSITIVE_REPORT_SHA256 = "7c88052348a78ab7e26728f492296fb57858cead1794b9d8cb6e57ad873a7b9b"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
CROSS_PROJECTION_RECEIPT = PREFIX + "n1_outer_own_source_joins_construction_check.json"
CROSS_PROJECTION_RECEIPT_SHA256 = "b67073446ac25da8736cae9c09ef963cafe3474f730a3c4d709f4cd47b6b83a6"
PHYSICAL_NORM_SOURCE = "lei_ren_part1_paper_compliant_physical_norm_family.py"
PHYSICAL_NORM_DATUM = "lei_ren_part1_paper_compliant_physical_norm_family.json"
N1_PARAMETER_IDENTITY_RECEIPT = PREFIX + "n1_core_fixed_point_functions_check.json"
SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
EXPECTED_NODE6_FIXED_LEADING_INPUT_N0_SHA256 = "42df26c8f2262256dbfc33da6ff555ebaf18f895d8014212f1304dabb680d6d9"
EXPECTED_NODE35560_PROVIDER = {
    "module": PREFIX + "n1_core_fixed_point_functions.py",
    "module_sha256": "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    "function": "CurrentOriginalN1CoreFixedPointFunctions.attach_oracle",
    "AST_sha256": "90c40b6f5111f7104975096c30cbc23657ee6317819ad9279c396b7cdc4bcf9f",
}

# These are the source files named by the accepted node-6/node-35560
# projection.  They remain ordinary authenticated records; no source is
# imported or executed by this module.
LOGC_PROJECTION_HASHES = {
    CROSS_PROJECTION_RECEIPT: CROSS_PROJECTION_RECEIPT_SHA256,
    PHYSICAL_NORM_SOURCE: "89c01b43baee8e89a4060f195f54cdf6aa2c65acb606457197ed28ef56dcf0e8",
    PHYSICAL_NORM_DATUM: "e8b97192b8d026a8fa2b22848bb07bae38bc876261d909fb70419b4898efb7e5",
    PREFIX + "whole_Z_long_Rm_finite_N.py": "dd5830925de6befe08254e4954cb25e4a3ec1d5dec448cf6684e3c6fade87c2d",
    PREFIX + "reference_restore_functions.py": "a86da4ef8791c65ce678a067a999f1f24945ac7c1ca0673f6bafc8c8969702b7",
    PREFIX + "whole_Z_long_Rm_signed_functions.py": "071f7978fb9c5e29b91d11763c0bbf4beaedc19a7f2d6b133ed1d97a9b2ae033",
    PREFIX + "reshape_terminal_kernels.py": "a846be160da9a9e4e14296d5a47da235ce899101191ca23f1da7290b4888b5bc",
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point.py": "dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f",
    "lei_ren_part1_paper_compliant_core_physical_field.py": "19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54",
    PREFIX + "whole_Z_all_N_outer_Rc_functions_check.json": "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158",
}
EXPECTED_LOGC_ASSIGNMENTS = (
    {
        "module": PREFIX + "whole_Z_long_Rm_finite_N.py",
        "module_sha256": "dd5830925de6befe08254e4954cb25e4a3ec1d5dec448cf6684e3c6fade87c2d",
        "class_name": "WholeZLongRmFiniteN",
        "function": "__init__",
        "AST_sha256": "ccbcff7ea37243b26ef62fca3331407b22577129d369fd7e83fdee46366930d1",
        "actual_assignments": {
            "params": "json.loads((HERE / kernels.PARAMS).read_bytes())",
            "self.logC": "read(c, params['selected_logCstar'])",
        },
    },
    {
        "module": PREFIX + "reshape_terminal_kernels.py",
        "module_sha256": "a846be160da9a9e4e14296d5a47da235ce899101191ca23f1da7290b4888b5bc",
        "class_name": None,
        "function": None,
        "AST_sha256": "7718b28f605bbe5abe9c57e14a6bd619eceb2f6e84b3cc34bf3c29f4c94b6188",
        "actual_assignments": {"PARAMS": "PREFIX + 'physical_norm_family.json'"},
    },
    {
        "module": PREFIX + "whole_Z_long_Rm_finite_N.py",
        "module_sha256": "dd5830925de6befe08254e4954cb25e4a3ec1d5dec448cf6684e3c6fade87c2d",
        "class_name": "WholeZLongRmFiniteN",
        "function": "owner",
        "AST_sha256": "8139c3859c2b0ba51090c7c9e6770dd5c247ecb4e39070b71b5c424dd80eb8f2",
        "actual_assignments": {
            "ref": "background.ActualReferenceRestoreFunctions(f, z, source.core.delta, terminal['actual_terminal_normalized_six_history_shapes'], fields['V'], P0, self.T, self.logC, self.provider)",
        },
    },
    {
        "module": PREFIX + "reference_restore_functions.py",
        "module_sha256": "a86da4ef8791c65ce678a067a999f1f24945ac7c1ca0673f6bafc8c8969702b7",
        "class_name": "ActualReferenceRestoreFunctions",
        "function": "__init__",
        "AST_sha256": "0219b5224273ec6ee3e22d8d524637cdabadf3e8ada9061c7bae7d85153b6380",
        "actual_assignments": {"self.logC": "c.mpf(logC)"},
    },
    {
        "module": PREFIX + "whole_Z_long_Rm_signed_functions.py",
        "module_sha256": "071f7978fb9c5e29b91d11763c0bbf4beaedc19a7f2d6b133ed1d97a9b2ae033",
        "class_name": "CurrentLongRadiusPhase",
        "function": "__init__",
        "AST_sha256": "2fe00a905b473d920d60a8c4bf0de29fcf492d2ca9187bef8909dca2797da160",
        "actual_assignments": {"self.logC": "long_owner.reference.logC"},
    },
    {
        "module": "lei_ren_part1_paper_compliant_current_core_common_fixed_point.py",
        "module_sha256": "dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f",
        "class_name": "CurrentCoreCommonFixedPoint",
        "function": "__init__",
        "AST_sha256": "9a30ab6417a3b7dab13ae5afcc193ff290fb22de8082beaff219b48bdc4c033f",
        "actual_assignments": {"self.core": "self.source.core"},
    },
    {
        "module": "lei_ren_part1_paper_compliant_core_physical_field.py",
        "module_sha256": "19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54",
        "class_name": "CompliantCorePhysicalField",
        "function": "__init__",
        "AST_sha256": "97db8ff5f248a806fdf30ba5fb76c92542ae9b787986f4174cf193c128afd4df",
        "actual_assignments": {
            "norm": "self.records['physical_norm_family']",
            "self.logC": "read(norm, 'selected_logCstar')",
        },
    },
)

# Every item is an accepted parent input pin.  Gzip records are authenticated
# from their compressed bytes before they are decoded.
RECORD_HASHES = {
    PREFIX + "whole_Z_all_N_control_connection.py": "aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d",
    PREFIX + "whole_Z_all_N_control_connection_check.py": "e1980b87bafe334264fede244219ec62b2c27b0bd63d7a05bbf4e0c2c993fb1a",
    PREFIX + "whole_Z_all_N_control_connection_check.json": "f451f291643653baa912a672ebdd914eb4f668cf3316bb37e3f7eea5b4310abb",
    PREFIX + "whole_Z_all_N_control_connection.json.gz": "588f411e9122b9a76e6a42efc2b5754137f6af66a9cae42f8d6a7c2c59920de8",
    PREFIX + "whole_Z_all_N_outer_Rc_functions.py": "9d73659a1c76163dddee6653af53f56ebffd7045986d871dc559ba1e053620b7",
    PREFIX + "whole_Z_all_N_outer_Rc_functions_check.py": "9fc3a4921ab1a7c633f3dd60acd175b65803871ab0cc7963e1b964d90942fda9",
    PREFIX + "whole_Z_all_N_outer_Rc_functions_check.json": "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158",
    PREFIX + "whole_Z_all_N_outer_Rc_functions.json.gz": "93364a44e86b79dfbe0a942f9d6e6157a0eb0b2e149dae58833cfc50a41d72d0",
    PREFIX + "C2_density_transport.py": "8bfa32c87c73c5bedf3056c7171707ae934dafd950bd0df874dbd3b3450379d1",
    PREFIX + "C2_density_transport_check.json": "c4d2c62dd67e9dd7e171fabbdb9ba18308e420c5e0c98d01147c708840b8f583",
    PREFIX + "C2_density_transport.json.gz": "0f347824bea798ea692eca7967d7acadba49fe4336989bea1a5c86c922260877",
    PREFIX + "prescribed_loop_C5_functions.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    PREFIX + "prescribed_loop_C5_functions_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    PREFIX + "prescribed_loop_C5_functions.json.gz": "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    PREFIX + "prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    PREFIX + "prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    PREFIX + "prescribed_bridge_C6_functions.json.gz": "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92",
    PREFIX + "prescribed_R100_Rm_C6_functions.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    PREFIX + "prescribed_R100_Rm_C6_functions_check.py": "8ad5e2edc1daba05763b44dd031ac1671d74f30b1733fc80570873e2ca12089a",
    PREFIX + "prescribed_R100_Rm_C6_functions_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    PREFIX + "prescribed_R100_Rm_C6_functions.json.gz": "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
    PREFIX + "n1_lower_collar_functions.py": "d00dfe285735844eabc985f99d94190caa1cfd187f668eb14a2686055cf175fa",
    PREFIX + "n1_lower_collar_functions_check.json": "1a7b6b370b5ecf7d6b7e97001e0311e9525698ecf3389fc3da608af308a13ba9",
    PREFIX + "n1_lower_collar_functions.json.gz": "b2f983ab934875ae9aeb24b292277911dd5647e3a23943ddcf2b3b6e82595cb5",
    PREFIX + "n1_core_fixed_point_functions.py": "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    PREFIX + "n1_core_fixed_point_functions_check.json": "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead",
    PREFIX + "n1_core_fixed_point_functions.json.gz": "999aefa1e4e285dff1ebc29e1bfc4c2fb58e5507acd01a0386fbd52dfdd41b43",
}
RECORD_HASHES.update(LOGC_PROJECTION_HASHES)

CHARTS = (
    "first_micro",
    "second_micro",
    "frozen_macro",
    "first_switch",
    "second_switch",
    "post_power",
    "long_reshape",
    "reference",
    "restoration",
)

# These handles are the frozen exact C2 prefix handles.  The C5 graph is an
# extension of that prefix, so the IDs remain valid in the full graph.
C2 = {
    "first_micro": dict(coordinate=23, lower=29, upper=1, offset=27, jacobian=8),
    "second_micro": dict(coordinate=328, lower=1, upper=48, offset=330, jacobian=8),
    "frozen_macro": dict(coordinate=638, lower=0, upper=1, offset=648, jacobian=654),
    "first_switch": dict(coordinate=935, lower=0, upper=1, offset=937, jacobian=9),
    "second_switch": dict(coordinate=1219, lower=0, upper=1, offset=1222, jacobian=9),
    "post_power": dict(coordinate=1505, lower=0, upper=1, offset=1512, jacobian=1517),
    "long_reshape": dict(coordinate=1796, lower=0, upper=1, offset=1798, jacobian=7),
    "reference": dict(coordinate=2080, lower=0, upper=1, offset=2088, jacobian=2086),
    "restoration": dict(coordinate=2377, lower=0, upper=1, offset=2378, jacobian=1),
}

# The saved C2 source descriptor records the final right-offset handle of each
# authenticated physical partition.  These are source-descriptor handles, not
# the compact graph ``offset`` handles above; the endpoint binding is checked
# against the compressed descriptor before any C2/C5 identity is admitted.
EXPECTED_C2_RIGHT_OFFSET_ENDPOINTS = {
    "first_micro": 31,
    "second_micro": 332,
    "frozen_macro": 650,
    "first_switch": 938,
    "second_switch": 1224,
    "post_power": 1513,
    "long_reshape": 1799,
    "reference": 2092,
    "restoration": 2380,
}

# C5 chart source coordinates/radii/measures/global offsets.  These are graph
# handles, not names used as function identity; every handle is checked below.
C5 = {
    "first_micro": dict(coordinate=66988, variable="original_prescribed_C6_micro_s", lower=97295, upper=1,
                         radius=68644, measure=8, global_offset=97299, base=65604),
    "second_micro": dict(coordinate=66988, variable="original_prescribed_C6_micro_s", lower=1, upper=48,
                          radius=68644, measure=8, global_offset=97299, base=65604),
    "frozen_macro": dict(coordinate=70265, variable="original_prescribed_C6_macro_fraction", lower=0, upper=1,
                          radius=71448, measure=66987, global_offset=101907, base=65604),
    "first_switch": dict(coordinate=73190, variable="true_prescribed_first_switch_s", lower=0, upper=1,
                          radius=85905, measure=9, global_offset=104535, base=31709),
    "second_switch": dict(coordinate=75949, variable="true_prescribed_second_switch_t", lower=0, upper=1,
                           radius=88600, measure=9, global_offset=107163, base=31709),
    "post_power": dict(coordinate=78472, variable="true_prescribed_post_switch_fraction", lower=0, upper=1,
                        radius=91107, measure=91103, global_offset=109572, base=31709),
    "long_reshape": dict(coordinate=79966, variable="true_prescribed_long_reshape_fraction", lower=0, upper=1,
                          radius=79970, measure=7, global_offset=111980, base=79968),
    "reference": dict(coordinate=81520, variable="true_prescribed_reference_fraction", lower=0, upper=1,
                       radius=81524, measure=81519, global_offset=114388, base=79968),
    "restoration": dict(coordinate=83010, variable="true_prescribed_restoration_t", lower=0, upper=1,
                         radius=83014, measure=1, global_offset=116796, base=79968),
}

EXPECTED_C2_VARIABLES = {
    "first_micro": "native_first_micro",
    "second_micro": "native_second_micro",
    "frozen_macro": "native_frozen_macro",
    "first_switch": "native_first_switch",
    "second_switch": "native_second_switch",
    "post_power": "native_post_power",
    "long_reshape": "native_long_reshape",
    "reference": "native_reference",
    "restoration": "native_restoration",
}

# Handle aliases are used only after the source-parameter row has been checked.
PARAMETER_HANDLES = {5: "logP", 6: "logC", 7: "T", 8: "hbB", 9: "hbS", 10: "sc", 35560: "true_core_logC"}


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8")
    ).hexdigest()


def _gate(gates: list[dict[str, Any]], name: str, **details: Any) -> None:
    item = {"gate": name}
    item.update(details)
    gates.append(item)


def _raw_bytes(record: Mapping[str, Any]) -> bytes | None:
    encoded = record.get("raw_bytes_base64")
    if isinstance(encoded, str):
        try:
            return base64.b64decode(encoded, validate=True)
        except (ValueError, binascii.Error):
            return None
    for key in ("raw_bytes", "source_bytes", "gzip_bytes"):
        value = record.get(key)
        if isinstance(value, bytes):
            return value
        if isinstance(value, bytearray):
            return bytes(value)
    return None


def _read_record(source_records: Mapping[str, Any], name: str, gates: list[dict[str, Any]]) -> dict[str, Any] | None:
    expected = RECORD_HASHES[name]
    record = source_records.get(name) if isinstance(source_records, Mapping) else None
    if not isinstance(record, Mapping):
        _gate(gates, "source_record_missing", record=name)
        return None
    expected_path = "experiments/root_st073/" + name
    if record.get("path") != expected_path or record.get("sha256") != expected:
        _gate(gates, "source_record_pin_changed", record=name, expected_path=expected_path, expected_sha256=expected)
        return None
    if name.endswith(".gz"):
        raw = _raw_bytes(record)
        if raw is None or hashlib.sha256(raw).hexdigest() != expected:
            _gate(gates, "gzip_raw_hash_changed", record=name, expected_sha256=expected)
            return None
        if not isinstance(record.get("source_text"), str) or not isinstance(record.get("decoded_sha256"), str):
            _gate(gates, "gzip_decoded_descriptor_missing", record=name)
            return None
        try:
            decoded = gzip.decompress(raw)
            text = decoded.decode("utf-8")
        except (OSError, EOFError, UnicodeDecodeError) as error:
            _gate(gates, "gzip_decode_failed", record=name, error=str(error))
            return None
        if record.get("source_text") != text or record.get("decoded_sha256") != hashlib.sha256(decoded).hexdigest():
            _gate(gates, "gzip_decoded_text_changed", record=name)
            return None
        return {"name": name, "text": text, "decoded": decoded, "raw_sha256": expected}
    text = record.get("source_text")
    if not isinstance(text, str) or hashlib.sha256(text.encode("utf-8")).hexdigest() != expected:
        _gate(gates, "source_text_hash_changed", record=name, expected_sha256=expected)
        return None
    return {"name": name, "text": text, "decoded": text.encode("utf-8"), "raw_sha256": None}


def _ast_sha256(node: ast.AST) -> str:
    return hashlib.sha256(ast.dump(node, include_attributes=False).encode("utf-8")).hexdigest()


def _ast_source_segment_sha256(text: str, node: ast.AST | None) -> str | None:
    if node is None:
        return None
    segment = ast.get_source_segment(text, node)
    if segment is None:
        return None
    try:
        normalized = textwrap.dedent(segment)
        return hashlib.sha256(ast.dump(ast.parse(normalized)).encode("utf-8")).hexdigest()
    except SyntaxError:
        return None


def _ast_path(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _ast_path(node.value)
        return f"{parent}.{node.attr}" if parent is not None else None
    return None


def _ast_call(node: ast.AST, function: str) -> ast.Call | None:
    if not isinstance(node, ast.Call) or _ast_path(node.func) != function:
        return None
    return node


def _ast_assignment_value(tree: ast.AST, target: str) -> ast.AST | None:
    def locate(assigned_target: ast.AST, assigned_value: ast.AST) -> ast.AST | None:
        if _ast_path(assigned_target) == target:
            return assigned_value
        if isinstance(assigned_target, (ast.Tuple, ast.List)) and isinstance(assigned_value, (ast.Tuple, ast.List)):
            if len(assigned_target.elts) != len(assigned_value.elts):
                return None
            for target_item, value_item in zip(assigned_target.elts, assigned_value.elts, strict=True):
                found = locate(target_item, value_item)
                if found is not None:
                    return found
        return None

    for item in ast.walk(tree):
        if isinstance(item, ast.Assign):
            for assigned_target in item.targets:
                found = locate(assigned_target, item.value)
                if found is not None:
                    return found
        elif isinstance(item, ast.AnnAssign) and _ast_path(item.target) == target:
            return item.value
    return None


def _ast_definition(tree: ast.AST, class_name: str | None, function_name: str | None) -> ast.AST | None:
    if class_name is None and function_name is None:
        return tree
    if class_name is None or function_name is None:
        return None
    owner = next((item for item in ast.walk(tree) if isinstance(item, ast.ClassDef) and item.name == class_name), None)
    if owner is None:
        return None
    return next((item for item in ast.walk(owner) if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == function_name), None)


def _check_c5_transport_operands(text: str, gates: list[dict[str, Any]]) -> dict[str, Any]:
    """Authenticate the actual C5 kernel operands from source AST only."""

    name = PREFIX + "prescribed_loop_C5_functions.py"
    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError as error:
        _gate(gates, "C5_transport_ast_invalid", record=name, error=str(error))
        return {"passed": False, "record": name, "error": str(error)}
    functions = [item for item in ast.walk(tree) if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == "transport"]
    if len(functions) != 1:
        _gate(gates, "C5_transport_definition_changed", record=name, count=len(functions))
        return {"passed": False, "record": name, "function_count": len(functions)}
    function = functions[0]
    kernel_value = _ast_assignment_value(function, "kernel")
    kernel_conditional_guard = False
    if isinstance(kernel_value, ast.IfExp):
        kernel_conditional_guard = (
            isinstance(kernel_value.test, ast.Name) and kernel_value.test.id == "rate"
            and _ast_path(kernel_value.orelse) == "g.one"
        )
        kernel_value = kernel_value.body
    width_value = _ast_assignment_value(function, "width")
    kernel_unary = _ast_call(kernel_value, "g.unary") if kernel_value is not None else None
    kernel_mul = None
    if kernel_unary is not None and len(kernel_unary.args) == 2:
        kernel_mul = _ast_call(kernel_unary.args[1], "g.mul")
    body_value = _ast_assignment_value(function, "body")
    body_scale = _ast_call(body_value, "a.scale") if body_value is not None else None
    body_kernel = None
    if body_scale is not None and len(body_scale.args) == 2:
        body_kernel = _ast_call(body_scale.args[1], "g.mul")
    decay_value = _ast_assignment_value(function, "decay")
    decay_conditional_guard = False
    if isinstance(decay_value, ast.IfExp):
        decay_conditional_guard = (
            isinstance(decay_value.test, ast.Name) and decay_value.test.id == "rate"
            and _ast_path(decay_value.orelse) == "g.one"
        )
        decay_value = decay_value.body
    decay_unary = _ast_call(decay_value, "g.unary") if decay_value is not None else None
    decay_mul = None
    if decay_unary is not None and len(decay_unary.args) == 2:
        decay_mul = _ast_call(decay_unary.args[1], "g.mul")
    integral_calls = [
        item for item in ast.walk(function)
        if isinstance(item, ast.Call) and _ast_path(item.func) == "a.integral"
    ]
    operand_checks = {
        "kernel_conditional_guard": kernel_conditional_guard,
        "kernel_exp": kernel_unary is not None and len(kernel_unary.args) == 2
            and isinstance(kernel_unary.args[0], ast.Constant) and kernel_unary.args[0].value == "exp",
        "kernel_mul_arity": kernel_mul is not None and len(kernel_mul.args) == 3,
        "kernel_rate_operand": False,
        "kernel_jacobian_operand": False,
        "kernel_coordinate_minus_transport_time": False,
        "width_jacobian_operand": False,
        "width_coordinate_minus_lower": False,
        "body_gmul_jacobian_kernel": (
            body_kernel is not None and len(body_kernel.args) == 2
            and isinstance(body_kernel.args[0], ast.Name) and body_kernel.args[0].id == "jacobian"
            and isinstance(body_kernel.args[1], ast.Name) and body_kernel.args[1].id == "kernel"
        ),
        "integral_body_binding": bool(integral_calls) and all(
            len(call.args) >= 1 and isinstance(call.args[0], ast.Name) and call.args[0].id == "body"
            for call in integral_calls
        ),
        "decay_conditional_guard": decay_conditional_guard,
        "decay_exp": decay_unary is not None and len(decay_unary.args) == 2
            and isinstance(decay_unary.args[0], ast.Constant) and decay_unary.args[0].value == "exp",
        "decay_width_operand": False,
        "integral_has_physical_jacobian_guard": False,
        "integral_has_predecessor_memory_guard": False,
    }
    if kernel_mul is not None and len(kernel_mul.args) == 3:
        first, second, third = kernel_mul.args
        rate_call = _ast_call(first, "g.constant")
        operand_checks["kernel_rate_operand"] = (
            rate_call is not None and len(rate_call.args) == 1 and isinstance(rate_call.args[0], ast.UnaryOp)
            and isinstance(rate_call.args[0].op, ast.USub) and isinstance(rate_call.args[0].operand, ast.Name)
            and rate_call.args[0].operand.id == "rate"
        )
        operand_checks["kernel_jacobian_operand"] = isinstance(second, ast.Name) and second.id == "jacobian"
        difference = _ast_call(third, "g.sub")
        operand_checks["kernel_coordinate_minus_transport_time"] = (
            difference is not None and len(difference.args) == 2
            and all(isinstance(arg, ast.Name) for arg in difference.args)
            and difference.args[0].id == "coordinate" and difference.args[1].id == "t"
        )
    width_call = _ast_call(width_value, "g.mul") if width_value is not None else None
    if width_call is not None and len(width_call.args) == 2:
        operand_checks["width_jacobian_operand"] = isinstance(width_call.args[0], ast.Name) and width_call.args[0].id == "jacobian"
        difference = _ast_call(width_call.args[1], "g.sub")
        operand_checks["width_coordinate_minus_lower"] = (
            difference is not None and len(difference.args) == 2
            and all(isinstance(arg, ast.Name) for arg in difference.args)
            and difference.args[0].id == "coordinate" and difference.args[1].id == "lower"
        )
    if decay_mul is not None and len(decay_mul.args) == 2:
        decay_rate, decay_width = decay_mul.args
        decay_rate_call = _ast_call(decay_rate, "g.constant")
        operand_checks["decay_width_operand"] = isinstance(decay_width, ast.Name) and decay_width.id == "width"
        operand_checks["decay_rate_operand"] = (
            decay_rate_call is not None and len(decay_rate_call.args) == 1
            and isinstance(decay_rate_call.args[0], ast.UnaryOp)
            and isinstance(decay_rate_call.args[0].op, ast.USub)
            and isinstance(decay_rate_call.args[0].operand, ast.Name)
            and decay_rate_call.args[0].operand.id == "rate"
        )
    else:
        operand_checks["decay_rate_operand"] = False
    if integral_calls:
        operand_checks["integral_has_physical_jacobian_guard"] = all(
            isinstance(
                {kw.arg: kw.value for kw in call.keywords}.get("physical_Jacobian_applied_once"),
                ast.Constant,
            )
            and {kw.arg: kw.value for kw in call.keywords}["physical_Jacobian_applied_once"].value is True
            for call in integral_calls
        )
        operand_checks["integral_has_predecessor_memory_guard"] = all(
            isinstance(
                {kw.arg: kw.value for kw in call.keywords}.get("predecessor_memory_never_reset"),
                ast.Constant,
            )
            and {kw.arg: kw.value for kw in call.keywords}["predecessor_memory_never_reset"].value is True
            for call in integral_calls
        )
    passed = all(operand_checks.values())
    if not passed:
        _gate(gates, "C5_transport_operands_changed", record=name, checks=operand_checks)
    return {
        "passed": passed,
        "record": name,
        "function": "transport",
        "AST_sha256": _ast_sha256(function),
        "operand_checks": operand_checks,
        "formula": "exp(g.mul(g.constant(-rate), jacobian, g.sub(coordinate, t)))",
        "source_executed": False,
    }


def _ast_ref_field(node: ast.AST, owner: str, field: str) -> bool:
    reference = _ast_call(node, "ref")
    if reference is None or len(reference.args) != 1:
        return False
    item = reference.args[0]
    return (
        isinstance(item, ast.Subscript)
        and isinstance(item.value, ast.Name)
        and item.value.id == owner
        and isinstance(item.slice, ast.Constant)
        and item.slice.value == field
    )


def _check_c2_right_offset_source_descriptor(
    source_text: str, descriptor_text: str, gates: list[dict[str, Any]],
) -> dict[str, Any]:
    """Bind C2 kernel endpoints to the authenticated saved source descriptor."""

    source_name = PREFIX + "C2_density_transport.py"
    descriptor_name = PREFIX + "C2_density_transport.json.gz"
    result: dict[str, Any] = {
        "passed": False,
        "source_record": source_name,
        "descriptor_record": descriptor_name,
        "source_executed": False,
        "source_ast_checks": {},
        "descriptor_checks": {},
        "endpoint_handles": {},
    }
    source_checks: dict[str, bool] = {}
    try:
        tree = ast.parse(source_text, filename=source_name)
    except SyntaxError as error:
        _gate(gates, "C2_right_offset_descriptor_changed", record=source_name, error=str(error))
        result["error"] = str(error)
        return result
    builds = [
        item for item in ast.walk(tree)
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == "build"
    ]
    source_checks["single_build_definition"] = len(builds) == 1
    build = builds[0] if len(builds) == 1 else None
    kernel_value = _ast_assignment_value(build, "kernel") if build is not None else None
    kernel_exp = _ast_call(kernel_value, "g.unary") if kernel_value is not None else None
    kernel_neg = None
    kernel_mul = None
    if kernel_exp is not None and len(kernel_exp.args) == 2:
        kernel_neg = _ast_call(kernel_exp.args[1], "g.neg")
    if kernel_neg is not None and len(kernel_neg.args) == 1:
        kernel_mul = _ast_call(kernel_neg.args[0], "g.mul")
    source_checks["kernel_exp"] = (
        kernel_exp is not None and len(kernel_exp.args) == 2
        and isinstance(kernel_exp.args[0], ast.Constant) and kernel_exp.args[0].value == "exp"
    )
    source_checks["kernel_neg"] = kernel_neg is not None and len(kernel_neg.args) == 1
    source_checks["kernel_mul_rate_difference"] = kernel_mul is not None and len(kernel_mul.args) == 2
    source_checks["kernel_rate_operand"] = False
    source_checks["kernel_right_offset_operand"] = False
    source_checks["kernel_offset_operand"] = False
    if kernel_mul is not None and len(kernel_mul.args) == 2:
        rate = _ast_call(kernel_mul.args[0], "g.constant")
        difference = _ast_call(kernel_mul.args[1], "g.sub")
        source_checks["kernel_rate_operand"] = (
            rate is not None and len(rate.args) == 1
            and isinstance(rate.args[0], ast.Name) and rate.args[0].id == "rate"
        )
        source_checks["kernel_right_offset_operand"] = (
            difference is not None and len(difference.args) == 2
            and _ast_ref_field(difference.args[0], "old", "right_offset")
        )
        source_checks["kernel_offset_operand"] = (
            difference is not None and len(difference.args) == 2
            and _ast_ref_field(difference.args[1], "old", "offset")
        )
    integral_calls = [
        item for item in ast.walk(build) if isinstance(item, ast.Call)
        and _ast_path(item.func) == "current.controls.integral"
    ] if build is not None else []
    source_checks["integral_calls_present"] = bool(integral_calls)
    body_bindings: list[bool] = []
    endpoint_bindings: list[bool] = []
    for call in integral_calls:
        body = _ast_call(call.args[1], "g.mul") if len(call.args) >= 2 else None
        body_bindings.append(
            body is not None and len(body.args) == 3
            and isinstance(body.args[0], ast.Name) and body.args[0].id == "kernel"
            and isinstance(body.args[1], ast.Name) and body.args[1].id == "q"
            and _ast_ref_field(body.args[2], "old", "Jacobian")
        )
        endpoint_bindings.append(
            len(call.args) >= 5
            and _ast_ref_field(call.args[3], "p", "lower")
            and _ast_ref_field(call.args[4], "p", "upper")
        )
    source_checks["integral_body_kernel_q_Jacobian"] = bool(body_bindings) and all(body_bindings)
    source_checks["integral_endpoint_bindings"] = bool(endpoint_bindings) and all(endpoint_bindings)
    result["source_ast_checks"] = source_checks

    descriptor_checks: dict[str, bool] = {}
    try:
        payload = json.loads(descriptor_text)
    except (TypeError, json.JSONDecodeError) as error:
        _gate(gates, "C2_right_offset_descriptor_changed", record=descriptor_name, error=str(error))
        result["error"] = str(error)
        return result
    transport = payload.get("actual_C2_density_transport_and_target_functions") if isinstance(payload, Mapping) else None
    windows = transport.get("windows") if isinstance(transport, Mapping) else None
    descriptor_checks["windows_list"] = isinstance(windows, list)
    if isinstance(windows, list):
        for chart, expected in EXPECTED_C2_RIGHT_OFFSET_ENDPOINTS.items():
            matches = [item for item in windows if isinstance(item, Mapping) and item.get("chart") == chart]
            descriptor_checks[f"{chart}.single_window"] = len(matches) == 1
            pieces = matches[0].get("original_physical_measure_and_partitions") if len(matches) == 1 else None
            descriptor_checks[f"{chart}.source_cells"] = (
                isinstance(pieces, list) and bool(pieces) and all(isinstance(piece, Mapping) for piece in pieces)
            )
            endpoint = pieces[-1].get("right_offset") if descriptor_checks[f"{chart}.source_cells"] else None
            descriptor_checks[f"{chart}.right_offset_endpoint"] = endpoint == expected
            result["endpoint_handles"][chart] = endpoint
    result["descriptor_checks"] = descriptor_checks
    result["passed"] = all(source_checks.values()) and all(descriptor_checks.values())
    if not result["passed"]:
        _gate(gates, "C2_right_offset_descriptor_changed", record=descriptor_name,
              source_ast_checks=source_checks, descriptor_checks=descriptor_checks)
    return result


def _source_evidence(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for name in RECORD_HASHES:
        row = _read_record(source_records, name, gates)
        if row is not None:
            records[name] = row

    source_tokens = {
        PREFIX + "whole_Z_all_N_control_connection.py": (
            "def exact_geometry", "P,C,T,B,S,sc,x=sy.symbols", "Y=sy.log(25)+1000+4*P", "frozen_macro", "post_power", "restoration",
        ),
        PREFIX + "C2_density_transport.py": (
            "def build", "g.unary('exp'", "right_offset", "old['right_offset']", "ref(old['offset'])",
            "ref(old['Jacobian'])", "actual_current_source_cells", "current.controls.integral",
        ),
        PREFIX + "prescribed_loop_C5_functions.py": (
            "def transport", "physical_measure", "predecessor_memory_never_reset", "offset=g.sub",
        ),
        PREFIX + "prescribed_bridge_C6_functions.py": (
            "self.hb=hb=", "self.Ra=Ra=", "self.Y=g.sub", "R=g.mul(Ra", "Rmacro=g.mul(Ra", "physical_measure=hb", "physical_measure=length",
        ),
        PREFIX + "prescribed_R100_Rm_C6_functions.py": (
            "native_parameter('hbS'", "native_parameter('T'", "self.post_length=g.sub", "self.logref=logref=g.mul", "self.reference_length=length=g.sub",
        ),
        PREFIX + "n1_lower_collar_functions.py": (
            "n['name']=='hbB'", "exact_binary_scalar", "physical_Ra_is_4_over_Lambda_not_4", "parameters={**self.core.parameters",
        ),
        PREFIX + "n1_core_fixed_point_functions.py": (
            "name='true_core_'+key", "logLambda=g.add", "Lambda,epsilon", "common.core.logC",
        ),
    }
    declaration_names = {
        PREFIX + "whole_Z_all_N_control_connection.py": ("exact_geometry", "build_graph"),
        PREFIX + "C2_density_transport.py": ("build",),
        PREFIX + "prescribed_loop_C5_functions.py": ("transport", "CurrentOriginalPrescribedLoopC5Functions"),
        PREFIX + "prescribed_bridge_C6_functions.py": ("CurrentOriginalPrescribedBridgeC6Functions",),
        PREFIX + "prescribed_R100_Rm_C6_functions.py": ("CurrentOriginalPrescribedR100RmC6Functions",),
        PREFIX + "n1_lower_collar_functions.py": ("CurrentOriginalN1LowerCollarFunctions",),
        PREFIX + "n1_core_fixed_point_functions.py": ("CurrentOriginalN1CoreFixedPointFunctions",),
    }
    declaration_lines: dict[str, dict[str, int]] = {}
    for name, wanted in declaration_names.items():
        row = records.get(name)
        if row is None:
            continue
        try:
            tree = ast.parse(row["text"], filename=name)
        except SyntaxError as error:
            _gate(gates, "source_ast_invalid", record=name, error=str(error))
            continue
        found: dict[str, int] = {}
        for item in ast.walk(tree):
            if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and item.name in wanted:
                found[item.name] = item.lineno
        for item in wanted:
            if item not in found:
                _gate(gates, "source_declaration_missing", record=name, declaration=item)
        declaration_lines[name] = found
    for name, wanted in source_tokens.items():
        row = records.get(name)
        if row is None:
            continue
        missing = [token for token in wanted if token not in row["text"]]
        if missing:
            _gate(gates, "source_formula_token_missing", record=name, tokens=missing)
    c5_transport = None
    c5_name = PREFIX + "prescribed_loop_C5_functions.py"
    if c5_name in records:
        c5_transport = _check_c5_transport_operands(records[c5_name]["text"], gates)
    c2_right_offset = None
    c2_name = PREFIX + "C2_density_transport.py"
    c2_descriptor_name = PREFIX + "C2_density_transport.json.gz"
    if c2_name in records and c2_descriptor_name in records:
        c2_right_offset = _check_c2_right_offset_source_descriptor(
            records[c2_name]["text"], records[c2_descriptor_name]["text"], gates,
        )
    receipts: dict[str, bool] = {}
    for name, row in records.items():
        if not name.endswith("_check.json") or name == CROSS_PROJECTION_RECEIPT:
            continue
        try:
            decoded = json.loads(row["text"])
        except json.JSONDecodeError as error:
            _gate(gates, "receipt_json_invalid", record=name, error=str(error))
            receipts[name] = False
            continue
        receipts[name] = decoded.get("all_passed") is True
        if not receipts[name]:
            _gate(gates, "receipt_not_accepted", record=name)
    return records, {"declarations": declaration_lines, "receipts": receipts, "token_checks": source_tokens,
                     "C5_transport_operands": c5_transport,
                     "C2_right_offset_source_descriptor": c2_right_offset}


def _check_logc_source_projection(records: Mapping[str, Mapping[str, Any]], gates: list[dict[str, Any]]) -> dict[str, Any]:
    """Check the accepted old-node/new-node source projection before aliasing."""

    result: dict[str, Any] = {
        "passed": False,
        "receipt": {"path": "experiments/root_st073/" + CROSS_PROJECTION_RECEIPT,
                    "sha256": CROSS_PROJECTION_RECEIPT_SHA256},
        "source_family": dict(SOURCE_FAMILY),
        "canonical_field": "selected_logCstar",
        "canonical_datum": PHYSICAL_NORM_DATUM,
        "canonical_recipe": "selected_logCstar=b.hi(2*max(explicit logarithmic lower endpoints))",
        "source_assignments": [],
        "input_pins": {},
        "parameter_identity_receipt": {},
        "ast_operand_dataflow": {},
        "scope": {
            "source_only": True,
            "source_executed": False,
            "saved_selected_logC_zero_width_supplemental": True,
            "interval_or_box_equality_used_as_identity": False,
            "numeric_value_selected": False,
            "operator_calls": [],
        },
    }
    receipt_row = records.get(CROSS_PROJECTION_RECEIPT)
    if receipt_row is None:
        _gate(gates, "logC_projection_receipt_missing", record=CROSS_PROJECTION_RECEIPT)
        return result
    try:
        receipt = json.loads(receipt_row["text"])
    except (KeyError, TypeError, json.JSONDecodeError) as error:
        _gate(gates, "logC_projection_receipt_invalid", record=CROSS_PROJECTION_RECEIPT, error=str(error))
        return result
    if receipt.get("source_family") != SOURCE_FAMILY:
        _gate(gates, "logC_projection_source_family_changed")
    old = receipt.get("checks", {}).get("actual_old_current_radius_sources", {})
    projection = old.get("selected_parameter_source_projection") if isinstance(old, Mapping) else None
    if not isinstance(projection, Mapping) or old.get("passed") is not True:
        _gate(gates, "logC_projection_receipt_not_accepted")
        return result
    expected_scalars = {
        "passed": True,
        "old_parameter_node": 6,
        "true_parameter_node": 35560,
        "canonical_datum": PHYSICAL_NORM_DATUM,
        "canonical_datum_sha256": LOGC_PROJECTION_HASHES[PHYSICAL_NORM_DATUM],
        "canonical_datum_key": "selected_logCstar",
        "accepted_true_core_parameter_identity_consumed": True,
        "same_selected_scalar_recipe_not_interval_overlap_or_new_common_radius": True,
    }
    scalar_checks = {key: projection.get(key) == value for key, value in expected_scalars.items()}
    for key, passed in scalar_checks.items():
        if not passed:
            _gate(gates, "logC_projection_receipt_scalar_changed", field=key,
                  expected=expected_scalars[key], actual=projection.get(key))
    assignments = projection.get("actual_source_assignments")
    assignments_match = assignments == list(EXPECTED_LOGC_ASSIGNMENTS)
    if not assignments_match:
        _gate(gates, "logC_projection_assignments_changed")
    result["source_assignments"] = list(assignments) if isinstance(assignments, list) else []
    result["assignment_scalars"] = scalar_checks

    receipt_inputs = receipt.get("input_hashes")
    pins_ok = isinstance(receipt_inputs, Mapping)
    for name, expected in LOGC_PROJECTION_HASHES.items():
        if name == CROSS_PROJECTION_RECEIPT:
            continue
        actual = receipt_inputs.get(name) if isinstance(receipt_inputs, Mapping) else None
        result["input_pins"][name] = {"sha256": expected, "receipt_sha256": actual, "passed": actual == expected}
        if actual != expected:
            pins_ok = False
            _gate(gates, "logC_projection_input_pin_changed", record=name, expected_sha256=expected,
                  actual_sha256=actual)

    datum: Mapping[str, Any] | None = None
    datum_row = records.get(PHYSICAL_NORM_DATUM)
    try:
        datum_value = json.loads(datum_row["text"]) if datum_row is not None else None
        datum = datum_value if isinstance(datum_value, Mapping) else None
    except (KeyError, TypeError, json.JSONDecodeError) as error:
        _gate(gates, "physical_norm_datum_invalid", error=str(error))
    definition = datum.get("definition") if isinstance(datum, Mapping) else None
    selected = datum.get("selected_logCstar") if isinstance(datum, Mapping) else None
    selected_definition = definition.get("selected_logC") if isinstance(definition, Mapping) else None
    datum_checks = {
        "canonical_field_present": isinstance(selected, Mapping),
        "definition_reuses_same_field": selected_definition == selected,
        "selection_recipe": isinstance(definition, Mapping)
        and definition.get("selection") == "logCstar=2*max(explicit logarithmic lower bounds)",
        "zero_width_supplemental": isinstance(selected, Mapping)
        and selected.get("width") == "0.0"
        and selected.get("lower") == selected.get("upper")
        and selected.get("lower_exact_mpf_tuple") == selected.get("upper_exact_mpf_tuple"),
    }
    for field, passed in datum_checks.items():
        if not passed:
            _gate(gates, "physical_norm_canonical_field_changed", field=field)
    result["canonical_datum_checks"] = datum_checks

    parameter_receipt_row = records.get(N1_PARAMETER_IDENTITY_RECEIPT)
    parameter_identity_checks = {
        "receipt_all_passed": False,
        "source_family": False,
        "actual_oracle_passed": False,
        "logC_matches_original_selected_source": False,
        "logC_leaf_uses_this_exact_source": False,
        "selected_j_logC_and_anchor_are_from_this_actual_core": False,
        "interval_overlap_not_used_as_a_function_identity": False,
    }
    try:
        parameter_receipt = json.loads(parameter_receipt_row["text"]) if parameter_receipt_row is not None else None
        exact_identity = (
            parameter_receipt.get("checks", {}).get("actual_oracle", {})
            .get("exact_parameter_and_native_P0_identity", {})
            if isinstance(parameter_receipt, Mapping) else {}
        )
        scalar_identity = exact_identity.get("actual_source_scalar_and_unique_anchor_identity", {})
        parameter_identity_checks = {
            "receipt_all_passed": parameter_receipt.get("all_passed") is True,
            "source_family": parameter_receipt.get("source_family") == SOURCE_FAMILY,
            "actual_oracle_passed": exact_identity.get("passed") is True,
            "logC_matches_original_selected_source": scalar_identity.get("logC_matches_original_selected_source") is True,
            "logC_leaf_uses_this_exact_source": scalar_identity.get("logC_leaf_uses_this_exact_source") is True,
            "selected_j_logC_and_anchor_are_from_this_actual_core": exact_identity.get("selected_j_logC_and_anchor_are_from_this_actual_core") is True,
            "interval_overlap_not_used_as_a_function_identity": exact_identity.get("interval_overlap_not_used_as_a_function_identity") is True,
        }
    except (KeyError, TypeError, json.JSONDecodeError) as error:
        _gate(gates, "logC_parameter_identity_receipt_invalid", record=N1_PARAMETER_IDENTITY_RECEIPT, error=str(error))
    for field, passed in parameter_identity_checks.items():
        if not passed:
            _gate(gates, "logC_parameter_identity_receipt_changed", field=field)
    result["parameter_identity_receipt"] = {
        "path": "experiments/root_st073/" + N1_PARAMETER_IDENTITY_RECEIPT,
        "sha256": RECORD_HASHES.get(N1_PARAMETER_IDENTITY_RECEIPT),
        "checks": parameter_identity_checks,
    }

    source_ast_checks: list[dict[str, Any]] = []
    source_ast_ok = True
    for assignment in EXPECTED_LOGC_ASSIGNMENTS:
        name = assignment["module"]
        row = records.get(name)
        check: dict[str, Any] = {"module": name, "class_name": assignment["class_name"],
                                  "function": assignment["function"], "passed": False}
        if row is None:
            _gate(gates, "logC_projection_source_missing", record=name)
            source_ast_checks.append(check)
            source_ast_ok = False
            continue
        try:
            tree = ast.parse(row["text"], filename=name)
            definition_node = _ast_definition(tree, assignment["class_name"], assignment["function"])
            actual_ast = _ast_sha256(definition_node) if definition_node is not None else None
            assignment_checks: dict[str, bool] = {}
            for target, expression in assignment["actual_assignments"].items():
                actual_value = _ast_assignment_value(definition_node, target) if definition_node is not None else None
                expected_value = ast.parse(expression, mode="eval").body
                assignment_checks[target] = actual_value is not None and _ast_sha256(actual_value) == _ast_sha256(expected_value)
            check.update(expected_AST_sha256=assignment["AST_sha256"], AST_sha256=actual_ast,
                         AST_hash_matches=actual_ast == assignment["AST_sha256"],
                         operand_checks=assignment_checks)
            check["passed"] = check["AST_hash_matches"] and all(assignment_checks.values())
        except (KeyError, SyntaxError, TypeError, ValueError) as error:
            check["error"] = str(error)
        if not check["passed"]:
            source_ast_ok = False
            _gate(gates, "logC_projection_AST_operand_changed", record=name,
                  class_name=assignment["class_name"], function=assignment["function"])
        source_ast_checks.append(check)

    provider_name = EXPECTED_NODE35560_PROVIDER["module"]
    provider_row = records.get(provider_name)
    provider_check: dict[str, Any] = {"module": provider_name, "function": EXPECTED_NODE35560_PROVIDER["function"],
                                      "passed": False}
    if provider_row is None:
        _gate(gates, "logC_projection_provider_missing", record=provider_name)
    else:
        try:
            provider_tree = ast.parse(provider_row["text"], filename=provider_name)
            owner, function_name = EXPECTED_NODE35560_PROVIDER["function"].split(".", 1)
            provider_class = next((item for item in ast.walk(provider_tree)
                                   if isinstance(item, ast.ClassDef) and item.name == owner), None)
            provider_function = next((item for item in ast.walk(provider_class) if isinstance(item, ast.FunctionDef)
                                      and item.name == function_name), None) if provider_class is not None else None
            provider_init = next((item for item in ast.walk(provider_class)
                                  if isinstance(item, ast.FunctionDef) and item.name == "__init__"), None) if provider_class is not None else None
            provider_operand = False
            if provider_init is not None:
                for call in ast.walk(provider_init):
                    if isinstance(call, ast.Call) and _ast_path(call.func) == "g.node":
                        source_provider = next((keyword.value for keyword in call.keywords
                                                if keyword.arg == "source_provider"), None)
                        if (_ast_call(source_provider, "ast") is not None
                                and len(_ast_call(source_provider, "ast").args) == 1
                                and _ast_path(_ast_call(source_provider, "ast").args[0]) == "self.attach_oracle"):
                            provider_operand = True
                            break
            provider_ast = _ast_source_segment_sha256(provider_row["text"], provider_function)
            provider_check.update(
                expected_AST_sha256=EXPECTED_NODE35560_PROVIDER["AST_sha256"],
                source_python_AST_sha256=provider_ast,
                AST_hash_matches=provider_ast == EXPECTED_NODE35560_PROVIDER["AST_sha256"],
                source_function_present=provider_function is not None,
                source_provider_operand=provider_operand,
                AST_descriptor_pinned_in_graph=True,
            )
            provider_check["passed"] = provider_check["AST_hash_matches"] and provider_operand
        except (SyntaxError, TypeError, ValueError) as error:
            provider_check["error"] = str(error)
        if not provider_check["passed"]:
            _gate(gates, "logC_projection_provider_AST_changed", record=provider_name)
    source_ast_checks.append(provider_check)

    signed_name = PREFIX + "whole_Z_long_Rm_signed_functions.py"
    signed_row = records.get(signed_name)
    radius_checks = {
        "RmRadiusPhase_base_constructor": False,
        "exact_source_logC_singleton": False,
        "exact_T_logC_singleton": False,
    }
    if signed_row is not None:
        try:
            signed_tree = ast.parse(signed_row["text"], filename=signed_name)
            phase = _ast_definition(signed_tree, "CurrentLongRadiusPhase", "__init__")
            base_value = _ast_assignment_value(phase, "self.base") if phase is not None else None
            base_call = _ast_call(base_value, "backend.phase.RmRadiusPhase") if base_value is not None else None
            radius_checks["RmRadiusPhase_base_constructor"] = (
                base_call is not None and len(base_call.args) == 3
                and isinstance(base_call.args[1], ast.Subscript)
                and _ast_path(base_call.args[1].value) == "identity"
                and isinstance(base_call.args[2], ast.Name) and base_call.args[2].id == "identity"
            )
            radius_checks["exact_source_logC_singleton"] = "self.logC._mpi_ != self.base.logC._mpi_" in signed_row["text"]
            radius_checks["exact_T_logC_singleton"] = "(self.T, self.logC)" in signed_row["text"] and "Actual source-selected singleton T/logC required" in signed_row["text"]
        except (SyntaxError, TypeError, ValueError):
            pass
    for field, passed in radius_checks.items():
        if not passed:
            _gate(gates, "RmRadiusPhase_source_guard_changed", field=field)

    norm_source_checks: dict[str, bool] = {
        "LogBounds_hi_upper_endpoint": False,
        "selected_logC_recipe": False,
    }
    norm_row = records.get(PHYSICAL_NORM_SOURCE)
    if norm_row is not None:
        try:
            norm_tree = ast.parse(norm_row["text"], filename=PHYSICAL_NORM_SOURCE)
            hi = _ast_definition(norm_tree, "LogBounds", "hi")
            hi_value = next((item.value for item in ast.walk(hi) if isinstance(item, ast.Return)), None) if hi is not None else None
            expected_hi = ast.parse("self.c.mpf(endpoints(self.c.mpf(value))[1])", mode="eval").body
            norm_source_checks["LogBounds_hi_upper_endpoint"] = hi_value is not None and _ast_sha256(hi_value) == _ast_sha256(expected_hi)
            recipe_value = _ast_assignment_value(norm_tree, "selected_logC")
            expected_recipe = ast.parse("b.hi(2*max(endpoints(v)[1] for v in restrictions.values()))", mode="eval").body
            norm_source_checks["selected_logC_recipe"] = recipe_value is not None and _ast_sha256(recipe_value) == _ast_sha256(expected_recipe)
        except (SyntaxError, TypeError, ValueError):
            pass
    for field, passed in norm_source_checks.items():
        if not passed:
            _gate(gates, "physical_norm_source_recipe_changed", field=field)

    result["source_ast_checks"] = source_ast_checks
    result["ast_operand_dataflow"] = {
        "C2_source_choice_chain": "LongRmFiniteN.selected_logCstar -> ActualReferenceRestoreFunctions -> CurrentLongRadiusPhase",
        "C5_source_choice_chain": "common.core.logC -> physical_norm_family.selected_logCstar",
        "C5_parameter_identity": "n1_core_fixed_point_functions_check.actual_source_scalar_and_unique_anchor_identity",
        "RmRadiusPhase": radius_checks,
        "physical_norm_family": norm_source_checks,
    }
    result["passed"] = (
        all(scalar_checks.values()) and assignments_match and pins_ok and all(datum_checks.values())
        and all(parameter_identity_checks.values()) and source_ast_ok and provider_check["passed"]
        and all(radius_checks.values()) and all(norm_source_checks.values())
    )
    if not result["passed"]:
        _gate(gates, "logC_source_projection_unproved")
    return result


def _node(nodes: Sequence[Any], index: int, gates: list[dict[str, Any]], label: str) -> Mapping[str, Any] | None:
    if type(index) is not int or index < 0 or index >= len(nodes):
        _gate(gates, "graph_handle_out_of_range", label=label, handle=index)
        return None
    row = nodes[index]
    if not isinstance(row, Mapping):
        _gate(gates, "graph_node_not_mapping", label=label, handle=index)
        return None
    return row


def _graph_parameter_source_guard(nodes: Sequence[Any], gates: list[dict[str, Any]]) -> dict[str, Any]:
    """Require exact source metadata for the two graph leaves before aliasing."""

    result: dict[str, Any] = {
        "passed": False,
        "source_family": dict(SOURCE_FAMILY),
        "node6": {},
        "node35560": {},
        "canonical_alias": "node35560 true_core_logC -> node6 logC",
    }
    node6 = _node(nodes, 6, gates, "C2 logC source parameter")
    node35560 = _node(nodes, 35560, gates, "C5 true_core_logC source parameter")
    if node6 is None or node35560 is None:
        _gate(gates, "logC_source_graph_metadata_missing")
        return result
    node6_expected = {
        "operation": "current_original_source_parameter",
        "name": "logC",
        "definition": "same selected_logCstar and CurrentLongRadiusPhase exact singleton tuple",
        "quantity_is_exact_original_parameter_not_range_endpoint": True,
        "Z_independent": True,
        "current_source_receipt": PREFIX + "whole_Z_all_N_outer_Rc_functions_check.json",
        "current_source_receipt_sha256": "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158",
    }
    node35560_expected = {
        "operation": "current_original_source_parameter",
        "name": "true_core_logC",
        "Z_independent": True,
        "exact_fixed_source_parameter": True,
        "source_path": "common.core.logC",
        "parameter_enclosure_not_selected_value": True,
        "anchor_definition": None,
    }
    checks: dict[str, bool] = {}
    for label, row, expected in (("node6", node6, node6_expected), ("node35560", node35560, node35560_expected)):
        for key, value in expected.items():
            checks[f"{label}.{key}"] = row.get(key) == value
            if not checks[f"{label}.{key}"]:
                _gate(gates, "logC_source_graph_attribute_changed", node=6 if label == "node6" else 35560,
                      attribute=key, expected=value, actual=row.get(key))
        checks[f"{label}.source_family"] = row.get("source_family") == SOURCE_FAMILY
        if not checks[f"{label}.source_family"]:
            _gate(gates, "logC_source_graph_family_changed", node=6 if label == "node6" else 35560)
    fixed_input = node6.get("fixed_leading_input_N0")
    checks["node6.fixed_leading_input_N0_exact"] = (
        type(fixed_input) is int
        and hashlib.sha256(str(fixed_input).encode("utf-8")).hexdigest() == EXPECTED_NODE6_FIXED_LEADING_INPUT_N0_SHA256
    )
    if not checks["node6.fixed_leading_input_N0_exact"]:
        _gate(gates, "logC_source_graph_fixed_input_changed", node=6)
    provider = node35560.get("source_provider")
    checks["node35560.source_provider"] = provider == EXPECTED_NODE35560_PROVIDER
    if not checks["node35560.source_provider"]:
        _gate(gates, "logC_source_graph_provider_changed", node=35560, expected=EXPECTED_NODE35560_PROVIDER,
              actual=provider)
    result["node6"] = {key: node6.get(key) for key in node6_expected}
    result["node6"]["fixed_leading_input_N0_sha256"] = hashlib.sha256(str(fixed_input).encode("utf-8")).hexdigest() if type(fixed_input) is int else None
    result["node35560"] = {key: node35560.get(key) for key in node35560_expected}
    result["node35560"]["source_provider"] = provider
    result["checks"] = checks
    result["passed"] = all(checks.values())
    if not result["passed"]:
        _gate(gates, "logC_source_graph_metadata_unproved")
    return result


def _arguments(row: Mapping[str, Any]) -> list[int] | None:
    values = row.get("arguments", row.get("args"))
    if not isinstance(values, list) or any(type(item) is not int for item in values):
        return None
    return values


def _same_args(values: Sequence[int] | None, expected: Sequence[int]) -> bool:
    return values is not None and list(values) == list(expected)


# Immutable symbolic expression constructors.  Fractions are exact and no
# decimal or interval value is ever produced.
Expr = tuple[Any, ...]
ZERO: Expr = ("const", 0, 1)
ONE: Expr = ("const", 1, 1)


def _const(value: Fraction | int) -> Expr:
    q = value if isinstance(value, Fraction) else Fraction(value)
    return ("const", q.numerator, q.denominator)


def _is_const(value: Expr) -> bool:
    return value and value[0] == "const"


def _fraction(value: Expr) -> Fraction:
    return Fraction(value[1], value[2])


def _sort_key(value: Expr) -> str:
    return repr(value)


Monomial = tuple[Expr, ...]
Polynomial = dict[Monomial, Fraction]


def _poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, Fraction(0)) + coefficient
        if result[monomial] == 0:
            del result[monomial]
    return result


def _poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(sorted(left_monomial + right_monomial, key=_sort_key))
            coefficient = left_coefficient * right_coefficient
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
            if result[monomial] == 0:
                del result[monomial]
    return result


def _canonical_factor(value: Expr) -> Expr:
    """Canonicalize one opaque polynomial factor without evaluating it."""

    tag = value[0]
    if tag == "const":
        return _const(_fraction(value))
    if tag == "atom":
        return value
    if tag == "log":
        return _log(value[1])
    if tag == "exp":
        inner = _from_polynomial(_polynomial(value[1]))
        return (tag, inner)
    if tag == "quot":
        numerator = _from_polynomial(_polynomial(value[1]))
        denominator = _from_polynomial(_polynomial(value[2]))
        if _is_const(numerator) and _is_const(denominator) and _fraction(denominator) != 0:
            return _const(_fraction(numerator) / _fraction(denominator))
        return numerator if denominator == ONE else ("quot", numerator, denominator)
    if tag in ("add", "mul", "neg"):
        return _from_polynomial(_polynomial(value))
    raise ValueError(f"unadmitted symbolic expression tag {tag!r}")


def _polynomial(value: Expr) -> Polynomial:
    """Convert an expression to an exact polynomial over opaque factors."""

    tag = value[0]
    if tag == "const":
        coefficient = _fraction(value)
        return {} if coefficient == 0 else {(): coefficient}
    if tag == "add":
        result: Polynomial = {}
        for item in value[1]:
            result = _poly_add(result, _polynomial(item))
        return result
    if tag == "mul":
        result: Polynomial = {(): Fraction(1)}
        for item in value[1]:
            result = _poly_mul(result, _polynomial(item))
            if not result:
                break
        return result
    if tag == "neg":
        return {monomial: -coefficient for monomial, coefficient in _polynomial(value[1]).items()}
    if tag == "log":
        canonical = _log(value[1])
        if canonical[0] != "log":
            return _polynomial(canonical)
        return {(canonical,): Fraction(1)}
    factor = _canonical_factor(value)
    return {(factor,): Fraction(1)}


def _monomial_expr(monomial: Monomial, coefficient: Fraction) -> Expr:
    factors = tuple(sorted(monomial, key=_sort_key))
    parts: list[Expr] = []
    if coefficient != 1:
        parts.append(_const(coefficient))
    parts.extend(factors)
    if not parts:
        return ONE
    if len(parts) == 1:
        return parts[0]
    return ("mul", tuple(parts))


def _from_polynomial(polynomial: Polynomial) -> Expr:
    terms = [
        _monomial_expr(monomial, coefficient)
        for monomial, coefficient in polynomial.items()
        if coefficient != 0
    ]
    if not terms:
        return ZERO
    terms.sort(key=_sort_key)
    if len(terms) == 1:
        return terms[0]
    return ("add", tuple(terms))


def _add(*values: Expr) -> Expr:
    polynomial: Polynomial = {}
    for value in values:
        polynomial = _poly_add(polynomial, _polynomial(value))
    return _from_polynomial(polynomial)


def _mul(*values: Expr) -> Expr:
    polynomial: Polynomial = {(): Fraction(1)}
    for value in values:
        polynomial = _poly_mul(polynomial, _polynomial(value))
        if not polynomial:
            return ZERO
    return _from_polynomial(polynomial)


def _neg(value: Expr) -> Expr:
    return _from_polynomial({monomial: -coefficient for monomial, coefficient in _polynomial(value).items()})


def _sub(left: Expr, right: Expr) -> Expr:
    return _add(left, _neg(right))


def _integer_prime_factors(value: int) -> dict[int, int]:
    factors: dict[int, int] = {}
    remainder = abs(value)
    divisor = 2
    while divisor * divisor <= remainder:
        while remainder % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            remainder //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if remainder > 1:
        factors[remainder] = factors.get(remainder, 0) + 1
    return factors


def _rational_log(value: Fraction) -> Expr:
    """Represent a positive rational logarithm by exact prime-log factors."""

    if value <= 0:
        return ("log", _const(value))
    if value == 1:
        return ZERO
    exponents = _integer_prime_factors(value.numerator)
    for prime, exponent in _integer_prime_factors(value.denominator).items():
        exponents[prime] = exponents.get(prime, 0) - exponent
        if exponents[prime] == 0:
            del exponents[prime]
    terms: list[Expr] = []
    for prime, exponent in sorted(exponents.items()):
        atom = _atom(f"log({prime})")
        terms.append(atom if exponent == 1 else _mul(_const(exponent), atom))
    return _add(*terms)


def _log(value: Expr) -> Expr:
    inner = _from_polynomial(_polynomial(value))
    if _is_const(inner):
        rational = _fraction(inner)
        if rational > 0:
            return _rational_log(rational)
    # All admitted source parameters are real. exp(a)>0 and each admitted
    # rational factor below is positive, so these are real logarithm laws.
    if inner[0] == "exp":
        return inner[1]
    if inner[0] == "mul" and all(
        factor[0] == "exp" or (_is_const(factor) and _fraction(factor) > 0)
        for factor in inner[1]
    ):
        return _add(*(_log(factor) for factor in inner[1]))
    return ("log", inner)


def _exp(value: Expr) -> Expr:
    return _canonical_factor(("exp", value))


def _quotient(left: Expr, right: Expr) -> Expr:
    return _canonical_factor(("quot", left, right))


def _atom(name: str) -> Expr:
    return ("atom", name)


def _render(value: Expr) -> str:
    tag = value[0]
    if tag == "const":
        q = Fraction(value[1], value[2])
        return str(q.numerator) if q.denominator == 1 else f"({q.numerator}/{q.denominator})"
    if tag == "atom":
        return str(value[1])
    if tag == "add":
        return " + ".join(_render(item) for item in value[1])
    if tag == "mul":
        return " * ".join(_render(item) for item in value[1])
    if tag == "log":
        return f"log({_render(value[1])})"
    if tag == "exp":
        return f"exp({_render(value[1])})"
    if tag == "quot":
        return f"({_render(value[1])})/({_render(value[2])})"
    return repr(value)


def _parameter_atom(nodes: Sequence[Any], index: int, gates: list[dict[str, Any]]) -> Expr | None:
    row = _node(nodes, index, gates, "parameter")
    if row is None:
        return None
    name = PARAMETER_HANDLES.get(index)
    if index == 5:
        if row.get("operation") != "sum" or not _same_args(_arguments(row), [3, 4]):
            _gate(gates, "logP_handle_shape_changed", handle=index)
            return None
        return _atom("logP")
    if row.get("operation") != "current_original_source_parameter" or row.get("name") != name:
        _gate(gates, "source_parameter_handle_changed", handle=index, expected=name,
              operation=row.get("operation"), actual_name=row.get("name"))
        return None
    return _atom(name)


def _expr(nodes: Sequence[Any], index: int, gates: list[dict[str, Any]], memo: dict[int, Expr] | None = None) -> Expr:
    if memo is None:
        memo = {}
    if index in memo:
        return memo[index]
    if index in PARAMETER_HANDLES:
        value = _parameter_atom(nodes, index, gates)
        if value is None:
            raise ValueError(f"parameter handle {index} is not pinned")
        memo[index] = value
        return value
    row = _node(nodes, index, gates, "expression")
    if row is None:
        raise ValueError(f"graph expression {index} unavailable")
    operation = row.get("operation")
    if operation == "exact_rational":
        numerator, denominator = row.get("numerator"), row.get("denominator")
        if type(numerator) is not int or type(denominator) is not int or denominator == 0:
            raise ValueError(f"non-scalar exact rational {index}")
        value = _const(Fraction(numerator, denominator))
    elif operation in ("sum", "add"):
        arguments = _arguments(row)
        if arguments is None:
            raise ValueError(f"sum arguments changed {index}")
        value = _add(*(_expr(nodes, item, gates, memo) for item in arguments))
    elif operation in ("product", "mul"):
        arguments = _arguments(row)
        if arguments is None:
            raise ValueError(f"product arguments changed {index}")
        value = _mul(*(_expr(nodes, item, gates, memo) for item in arguments))
    elif operation == "negative":
        argument = row.get("argument")
        if type(argument) is not int:
            raise ValueError(f"negative argument changed {index}")
        value = _neg(_expr(nodes, argument, gates, memo))
    elif operation == "analytic_unary":
        argument = row.get("argument")
        if type(argument) is not int:
            raise ValueError(f"analytic argument changed {index}")
        argument_value = _expr(nodes, argument, gates, memo)
        name = row.get("name")
        if name == "log":
            value = _log(argument_value)
        elif name == "exp":
            value = _exp(argument_value)
        else:
            raise ValueError(f"unadmitted analytic operator {name!r}")
    elif operation in ("positive_quotient", "quotient"):
        numerator, denominator = row.get("numerator"), row.get("denominator")
        if type(numerator) is not int or type(denominator) is not int:
            raise ValueError(f"quotient operands changed {index}")
        value = _quotient(_expr(nodes, numerator, gates, memo), _expr(nodes, denominator, gates, memo))
    elif operation == "current_original_source_parameter":
        name = row.get("name")
        if not isinstance(name, str):
            raise ValueError(f"source parameter name changed {index}")
        value = _atom(name)
    else:
        raise ValueError(f"unadmitted graph operation {operation!r} at {index}")
    memo[index] = value
    return value


def _affine(
    nodes: Sequence[Any], index: int, coordinate: int, gates: list[dict[str, Any]],
    memo: dict[tuple[int, int], tuple[Expr, Expr]] | None = None,
) -> tuple[Expr, Expr]:
    """Return exact ``slope, intercept`` for an admitted graph expression."""

    if memo is None:
        memo = {}
    key = (index, coordinate)
    if key in memo:
        return memo[key]
    if index in PARAMETER_HANDLES:
        result = (ZERO, _expr(nodes, index, gates))
        memo[key] = result
        return result
    if index == coordinate:
        result = (ONE, ZERO)
        memo[key] = result
        return result
    row = _node(nodes, index, gates, "affine")
    if row is None:
        raise ValueError(f"affine node {index} unavailable")
    operation = row.get("operation")
    if operation in ("bound_variable", "symbol"):
        raise ValueError(f"unrelated bound variable in affine expression {index}")
    if operation in ("exact_rational", "current_original_source_parameter"):
        result = (ZERO, _expr(nodes, index, gates))
        memo[key] = result
        return result
    if operation == "analytic_unary":
        argument = row.get("argument")
        if type(argument) is not int:
            raise ValueError(f"analytic argument changed {index}")
        argument_slope, _ = _affine(nodes, argument, coordinate, gates, memo)
        if argument_slope != ZERO:
            raise ValueError(f"coordinate-dependent analytic unary at {index}")
        result = (ZERO, _expr(nodes, index, gates))
        memo[key] = result
        return result
    if operation in ("sum", "add"):
        arguments = _arguments(row)
        if arguments is None:
            raise ValueError(f"sum arguments changed {index}")
        slopes, intercepts = zip(*(_affine(nodes, item, coordinate, gates, memo) for item in arguments), strict=False)
        result = (_add(*slopes), _add(*intercepts))
        memo[key] = result
        return result
    if operation == "negative":
        argument = row.get("argument")
        if type(argument) is not int:
            raise ValueError(f"negative argument changed {index}")
        slope, intercept = _affine(nodes, argument, coordinate, gates, memo)
        result = (_neg(slope), _neg(intercept))
        memo[key] = result
        return result
    if operation in ("product", "mul"):
        arguments = _arguments(row)
        if arguments is None:
            raise ValueError(f"product arguments changed {index}")
        slope, intercept = ZERO, ONE
        for item in arguments:
            item_slope, item_intercept = _affine(nodes, item, coordinate, gates, memo)
            if slope != ZERO and item_slope != ZERO:
                raise ValueError(f"non-affine product at {index}")
            slope, intercept = (
                _add(_mul(slope, item_intercept), _mul(intercept, item_slope)),
                _mul(intercept, item_intercept),
            )
        result = (slope, intercept)
        memo[key] = result
        return result
    if operation in ("positive_quotient", "quotient"):
        numerator, denominator = row.get("numerator"), row.get("denominator")
        if type(numerator) is not int or type(denominator) is not int:
            raise ValueError(f"quotient operands changed {index}")
        n_slope, n_intercept = _affine(nodes, numerator, coordinate, gates, memo)
        d_slope, d_intercept = _affine(nodes, denominator, coordinate, gates, memo)
        if d_slope != ZERO:
            raise ValueError(f"variable denominator at {index}")
        result = (_quotient(n_slope, d_intercept), _quotient(n_intercept, d_intercept))
        memo[key] = result
        return result
    raise ValueError(f"unadmitted affine graph operation {operation!r} at {index}")


def _bound_name(nodes: Sequence[Any], index: int, gates: list[dict[str, Any]], label: str) -> str | None:
    row = _node(nodes, index, gates, label)
    if row is None:
        return None
    if row.get("operation") != "bound_variable" or not isinstance(row.get("name"), str):
        _gate(gates, "coordinate_not_typed_bound_variable", chart=label, handle=index)
        return None
    return row["name"]


def _check_ra_recipe(nodes: Sequence[Any], gates: list[dict[str, Any]]) -> bool:
    ok = True
    expected = {
        2: ("exact_rational", 40, 1),
        13: ("exact_rational", 4, 1),
        28: ("exact_rational", 1, 2),
        48: ("exact_rational", 2, 1),
        643: ("exact_rational", 1000, 1),
        31709: ("exact_rational", 100, 1),
        79968: ("exact_rational", 110, 1),
    }
    for index, (operation, numerator, denominator) in expected.items():
        row = _node(nodes, index, gates, "Ra recipe")
        if row is None or row.get("operation") != operation or row.get("numerator") != numerator or row.get("denominator") != denominator:
            _gate(gates, "Ra_scalar_recipe_changed", handle=index)
            ok = False
    row3 = _node(nodes, 3, gates, "logP exp recipe")
    row4 = _node(nodes, 4, gates, "logP scalar recipe")
    row5 = _node(nodes, 5, gates, "logP")
    if row3 is None or row3.get("operation") != "analytic_unary" or row3.get("name") != "exp" or row3.get("argument") != 2:
        _gate(gates, "logP_exp_recipe_changed", handle=3)
        ok = False
    if row4 is None or row4.get("operation") != "exact_rational" or row4.get("numerator") != 11 or row4.get("denominator") != 1:
        _gate(gates, "logP_scalar_recipe_changed", handle=4)
        ok = False
    if row5 is None or row5.get("operation") not in ("sum", "add") or not _same_args(_arguments(row5), [3, 4]):
        _gate(gates, "logP_source_expression_changed", handle=5)
        ok = False
    row14 = _node(nodes, 14, gates, "4logP")
    row35562 = _node(nodes, 35562, gates, "logLambda")
    row35564 = _node(nodes, 35564, gates, "negative_logLambda")
    row35565 = _node(nodes, 35565, gates, "epsilon")
    row65604 = _node(nodes, 65604, gates, "Ra")
    if row14 is None or row5 is None or row14.get("operation") != "product" or not _same_args(_arguments(row14), [13, 5]):
        _gate(gates, "4logP_recipe_changed")
        ok = False
    if row35562 is None or row35562.get("operation") != "sum" or not _same_args(_arguments(row35562), [14, 643]):
        _gate(gates, "logLambda_recipe_changed")
        ok = False
    if row35564 is None or row35564.get("operation") != "negative" or row35564.get("argument") != 35562:
        _gate(gates, "negative_logLambda_recipe_changed")
        ok = False
    if row35565 is None or row35565.get("operation") != "analytic_unary" or row35565.get("name") != "exp" or row35565.get("argument") != 35564:
        _gate(gates, "epsilon_recipe_changed")
        ok = False
    if row65604 is None or row65604.get("operation") != "product" or not _same_args(_arguments(row65604), [13, 35565]):
        _gate(gates, "Ra_recipe_changed")
        ok = False
    phase = _node(nodes, 24954, gates, "common offset origin")
    if phase is None or phase.get("operation") != "product" or not _same_args(_arguments(phase), [8, 10, 28]):
        _gate(gates, "common_offset_origin_recipe_changed")
        ok = False
    return ok


def _global_shape(nodes: Sequence[Any], chart: str, data: Mapping[str, Any], gates: list[dict[str, Any]]) -> bool:
    row = _node(nodes, int(data["global_offset"]), gates, f"{chart} global offset")
    if row is None or row.get("operation") not in ("sum", "add"):
        _gate(gates, "C5_global_offset_shape_changed", chart=chart)
        return False
    arguments = _arguments(row)
    if arguments is None or len(arguments) != 2:
        _gate(gates, "C5_global_offset_arity_changed", chart=chart)
        return False
    radius = int(data["radius"])
    positive_log = False
    phase_negative = False
    for argument in arguments:
        item = _node(nodes, argument, gates, f"{chart} global offset term")
        if item is None:
            continue
        if item.get("operation") == "analytic_unary" and item.get("name") == "log":
            quotient = item.get("argument")
            qrow = _node(nodes, quotient, gates, f"{chart} radius quotient") if type(quotient) is int else None
            if qrow is not None and qrow.get("operation") == "positive_quotient" and qrow.get("numerator") == radius and qrow.get("denominator") == 65604:
                positive_log = True
        if item.get("operation") == "negative" and item.get("argument") == 24954:
            phase_negative = True
    if not (positive_log and phase_negative):
        _gate(gates, "C5_global_offset_operands_changed", chart=chart)
    return positive_log and phase_negative


def _radius_slope(
    nodes: Sequence[Any], chart: str, data: Mapping[str, Any], gates: list[dict[str, Any]],
) -> tuple[Expr, Expr, bool]:
    radius = int(data["radius"])
    coordinate = int(data["coordinate"])
    row = _node(nodes, radius, gates, f"{chart} radius")
    if row is None or row.get("operation") not in ("product", "mul"):
        _gate(gates, "C5_radius_shape_changed", chart=chart)
        return ZERO, ZERO, False
    arguments = _arguments(row)
    if arguments is None:
        _gate(gates, "C5_radius_arguments_changed", chart=chart)
        return ZERO, ZERO, False
    exp_arguments = []
    base_arguments = []
    for argument in arguments:
        child = _node(nodes, argument, gates, f"{chart} radius child")
        if child is not None and child.get("operation") == "analytic_unary" and child.get("name") == "exp":
            exp_arguments.append(argument)
        else:
            base_arguments.append(argument)
    if len(exp_arguments) != 1:
        _gate(gates, "C5_radius_exp_factor_changed", chart=chart)
        return ZERO, ZERO, False
    try:
        exponent_node = _node(nodes, exp_arguments[0], gates, f"{chart} radius exponential")
        exponent_argument = exponent_node.get("argument") if exponent_node is not None else None
        if type(exponent_argument) is not int:
            raise ValueError("radius exponential argument unavailable")
        slope, intercept = _affine(nodes, exponent_argument, coordinate, gates)
        for base in base_arguments:
            base_slope, _ = _affine(nodes, base, coordinate, gates)
            if base_slope != ZERO:
                raise ValueError("radius base depends on coordinate")
    except ValueError as error:
        _gate(gates, "C5_radius_not_affine_exponential", chart=chart, error=str(error))
        return ZERO, ZERO, False
    expected_base = int(data["base"])
    if base_arguments != [expected_base]:
        _gate(gates, "C5_radius_base_changed", chart=chart, expected=expected_base, actual=base_arguments)
        return slope, intercept, False
    return slope, intercept, True


def _expected_forms(c_symbol: str) -> dict[str, tuple[Expr, Expr, Expr]]:
    p = _atom("logP")
    c = _atom(c_symbol)
    b = _atom("hbB")
    s = _atom("hbS")
    t = _atom("T")
    sc = _atom("sc")
    y = _add(_log(_const(25)), _const(1000), _mul(_const(4), p))
    q = _log(_const(Fraction(11, 10)))
    phase = _mul(_const(Fraction(1, 2)), b, sc)
    g = _add(_mul(_const(10), c), _mul(_const(10), p), _neg(t), _const(-8))
    return {
        "first_micro": (b, _neg(phase), b),
        "second_micro": (b, _neg(phase), b),
        "frozen_macro": (_sub(y, _mul(_const(2), b)), _sub(_mul(_const(2), b), phase), _sub(y, _mul(_const(2), b))),
        "first_switch": (s, _sub(y, phase), s),
        "second_switch": (s, _add(y, s, _neg(phase)), s),
        "post_power": (_sub(q, _mul(_const(2), s)), _add(y, _mul(_const(2), s), _neg(phase)), _sub(q, _mul(_const(2), s))),
        "long_reshape": (t, _add(y, q, _neg(phase)), t),
        "reference": (g, _add(y, q, t, _neg(phase)), g),
        "restoration": (ONE, _add(y, q, _mul(_const(10), _add(c, p)), _const(-8), _neg(phase)), ONE),
    }


def _alias_logc(value: Expr, enabled: bool) -> Expr:
    """Apply the authenticated node 35560 -> node 6 source alias."""

    if not enabled:
        return value
    tag = value[0]
    if tag == "atom":
        return _atom("logC") if value[1] == "true_core_logC" else value
    if tag == "add":
        return _add(*(_alias_logc(item, enabled) for item in value[1]))
    if tag == "mul":
        return _mul(*(_alias_logc(item, enabled) for item in value[1]))
    if tag == "neg":
        return _neg(_alias_logc(value[1], enabled))
    if tag == "log":
        return _log(_alias_logc(value[1], enabled))
    if tag == "exp":
        return _exp(_alias_logc(value[1], enabled))
    if tag == "quot":
        return _quotient(_alias_logc(value[1], enabled), _alias_logc(value[2], enabled))
    return value


def _check_frozen_macro_measure_tree(nodes: Sequence[Any], gates: list[dict[str, Any]]) -> bool:
    """Require the saved C5 frozen-macro measure tree by exact node handles."""

    checks: dict[str, bool] = {}
    root = _node(nodes, 66987, gates, "frozen macro measure root")
    root_args = _arguments(root) if root is not None else None
    checks["66987=sum[66984,66986]"] = (
        root is not None and root.get("operation") in ("sum", "add") and _same_args(root_args, [66984, 66986])
    )
    first = _node(nodes, 66984, gates, "frozen macro measure first sum")
    first_args = _arguments(first) if first is not None else None
    checks["66984=sum[66981,66983]"] = (
        first is not None and first.get("operation") in ("sum", "add") and _same_args(first_args, [66981, 66983])
    )
    log_31709 = _node(nodes, 66981, gates, "frozen macro log 31709")
    checks["66981=log(31709)"] = (
        log_31709 is not None and log_31709.get("operation") == "analytic_unary"
        and log_31709.get("name") == "log" and log_31709.get("argument") == 31709
    )
    negative_log_ra = _node(nodes, 66983, gates, "frozen macro negative log Ra")
    checks["66983=-66982"] = (
        negative_log_ra is not None and negative_log_ra.get("operation") == "negative"
        and negative_log_ra.get("argument") == 66982
    )
    log_ra = _node(nodes, 66982, gates, "frozen macro log Ra")
    checks["66982=log(65604)"] = (
        log_ra is not None and log_ra.get("operation") == "analytic_unary"
        and log_ra.get("name") == "log" and log_ra.get("argument") == 65604
    )
    negative_width = _node(nodes, 66986, gates, "frozen macro negative width")
    checks["66986=-66985"] = (
        negative_width is not None and negative_width.get("operation") == "negative"
        and negative_width.get("argument") == 66985
    )
    width = _node(nodes, 66985, gates, "frozen macro width")
    width_args = _arguments(width) if width is not None else None
    checks["66985=48*8"] = (
        width is not None and width.get("operation") in ("product", "mul") and _same_args(width_args, [48, 8])
    )
    passed = all(checks.values())
    if not passed:
        _gate(gates, "frozen_macro_measure_tree_changed", checks=checks)
    return passed


def _c5_measure_normal_form(
    nodes: Sequence[Any], chart: str, data: Mapping[str, Any], expected: tuple[Expr, Expr, Expr],
    gates: list[dict[str, Any]], ra_ok: bool,
) -> tuple[Expr, bool]:
    measure = int(data["measure"])
    # These physical measures are exact shared handles in the source graph.
    if chart in ("first_micro", "second_micro", "first_switch", "second_switch", "long_reshape", "restoration"):
        try:
            value = _expr(nodes, measure, gates)
        except ValueError as error:
            _gate(gates, "C5_measure_expression_unavailable", chart=chart, error=str(error))
            return ZERO, False
        return value, value == expected[0]
    if chart == "frozen_macro":
        tree_ok = _check_frozen_macro_measure_tree(nodes, gates)
        if ra_ok and tree_ok:
            return expected[0], True
        _gate(gates, "frozen_macro_measure_normalization_unproved", chart=chart)
        return ZERO, False
    if chart == "post_power":
        try:
            value = _expr(nodes, measure, gates)
        except ValueError as error:
            _gate(gates, "C5_post_power_measure_unavailable", chart=chart, error=str(error))
            return ZERO, False
        return value, value == expected[0]
    if chart == "reference":
        try:
            value = _expr(nodes, measure, gates)
        except ValueError as error:
            _gate(gates, "C5_reference_measure_unavailable", chart=chart, error=str(error))
            return ZERO, False
        return value, value == expected[0]
    raise ValueError(chart)


def _domain_expr(nodes: Sequence[Any], index: int, gates: list[dict[str, Any]]) -> Expr:
    return _expr(nodes, index, gates)


def _row_result(
    nodes: Sequence[Any], chart: str, source_ok: bool, logc_alias_ok: bool, ra_ok: bool,
    gates: list[dict[str, Any]],
) -> dict[str, Any]:
    c2, c5 = C2[chart], C5[chart]
    expected_c2 = _expected_forms("logC")[chart]
    expected_c5 = _expected_forms("true_core_logC")[chart]
    row_gates: list[dict[str, Any]] = []
    before = len(gates)
    # Use one shared gate list for deterministic root diagnostics, then keep a
    # chart-local projection in the returned row.
    try:
        c2_name = _bound_name(nodes, c2["coordinate"], gates, chart + " C2 coordinate")
        c5_name = _bound_name(nodes, c5["coordinate"], gates, chart + " C5 coordinate")
        c2_slope, c2_intercept = _affine(nodes, c2["offset"], c2["coordinate"], gates)
        c2_j = _expr(nodes, c2["jacobian"], gates)
        c5_slope, c5_exponent_intercept, radius_ok = _radius_slope(nodes, chart, c5, gates)
        c5_j, c5_measure_ok = _c5_measure_normal_form(nodes, chart, c5, expected_c5, gates, ra_ok)
        c2_offset_ok = (c2_slope, c2_intercept) == (expected_c2[0], expected_c2[1])
        c2_j_ok = c2_j == expected_c2[2]
        c2_derivative_ok = c2_offset_ok and c2_slope == c2_j and c2_j_ok
        c5_derivative_ok = radius_ok and c5_measure_ok and c5_slope == c5_j
        c2_domain = (_domain_expr(nodes, c2["lower"], gates), _domain_expr(nodes, c2["upper"], gates))
        c5_domain = (_domain_expr(nodes, c5["lower"], gates), _domain_expr(nodes, c5["upper"], gates))
        domain_ok = c2_domain == c5_domain
        variable_ok = c2_name == EXPECTED_C2_VARIABLES[chart] and c5_name == c5["variable"]
        global_shape_ok = _global_shape(nodes, chart, c5, gates)
        # The C5 offset is log(base/Ra) + exponent - hbB*sc/2.  The base
        # normalization is exact: Ra gives zero, 100 gives Y, and 110 gives Y+Q.
        base_term = ZERO if c5["base"] == 65604 else (
            _add(_log(_const(100)), _neg(_log(_const(4))), _const(1000), _mul(_const(4), _atom("logP")))
            if c5["base"] == 31709 else
            _add(_log(_const(100)), _neg(_log(_const(4))), _const(1000), _mul(_const(4), _atom("logP")), _log(_const(Fraction(11, 10))))
        )
        phase = _mul(_const(Fraction(1, 2)), _atom("hbB"), _atom("sc"))
        c5_intercept = _sub(_add(base_term, c5_exponent_intercept), phase)
        c5_offset_ok = global_shape_ok and (c5_slope, c5_intercept) == (expected_c5[0], expected_c5[1])
        expected_c5_aliased = tuple(_alias_logc(item, logc_alias_ok) for item in expected_c5)
        c5_j_aliased = _alias_logc(c5_j, logc_alias_ok)
        actual_j_equal = c2_j == c5_j_aliased
        global_ok = c2_offset_ok and c5_offset_ok and (expected_c2[1] == expected_c5_aliased[1])
        cross_ok = source_ok and logc_alias_ok and c2_derivative_ok and c5_derivative_ok and domain_ok and variable_ok and actual_j_equal
        exact_forms = {
            "c2_offset": {"slope": _render(c2_slope), "intercept": _render(c2_intercept)},
            "c2_jacobian": _render(c2_j),
            "c5_radius_exponent": {"slope": _render(c5_slope), "intercept": _render(c5_exponent_intercept)},
            "c5_physical_measure": _render(c5_j),
            "expected_c2": {"slope": _render(expected_c2[0]), "intercept": _render(expected_c2[1]), "jacobian": _render(expected_c2[2])},
            "expected_c5": {"slope": _render(expected_c5[0]), "intercept": _render(expected_c5[1]), "jacobian": _render(expected_c5[2])},
            "expected_c5_after_logC_alias": {"slope": _render(expected_c5_aliased[0]), "intercept": _render(expected_c5_aliased[1]), "jacobian": _render(expected_c5_aliased[2])},
            "c5_global_offset": {"slope": _render(c5_slope), "intercept": _render(c5_intercept)},
        }
        row = {
            "chart": chart,
            "field": "m",
            "c2_handles": {"coordinate": c2["coordinate"], "lower": c2["lower"], "upper": c2["upper"], "offset": c2["offset"], "Jacobian": c2["jacobian"]},
            "c5_handles": {"coordinate": c5["coordinate"], "lower": c5["lower"], "upper": c5["upper"], "radius": c5["radius"], "physical_measure": c5["measure"], "global_offset": c5["global_offset"]},
            "typed_alpha_rename": {"c2": {"handle": c2["coordinate"], "name": c2_name}, "c5": {"handle": c5["coordinate"], "name": c5_name}, "real_coordinate": "x", "domain_equal": domain_ok},
            "pullback_domains": {"c2": [_render(c2_domain[0]), _render(c2_domain[1])], "c5": [_render(c5_domain[0]), _render(c5_domain[1])]},
            "exact_forms": exact_forms,
            "offset_derivative_proved_on_both": c2_derivative_ok and c5_derivative_ok,
            "actual_c2_J_equals_actual_c5_J": actual_j_equal,
            "memory_width_identity_proved": cross_ok,
            "cross_archive_kernel_identity_proved": cross_ok,
            "global_offset_identity_proved": global_ok,
            "kernel_forms": {"C2": "exp[-rate*(right_offset-offset(t))]", "C5": "exp[-rate*J*(coordinate-t)]"},
            "source_formula_pins": {"C2_exact_geometry": True, "C2_kernel": True, "C5_transport": True, "Ra_recipe": ra_ok,
                                    "logC_source_projection": source_ok, "logC_alias_applied": logc_alias_ok},
            "admitted": cross_ok,
        }
        # A row-local gate projection is useful to reviewers, while the root
        # gate list retains deterministic all-source diagnostics.
        row["diagnostic_gate_count_before"] = before
        row["diagnostic_gate_count_after"] = len(gates)
        return row
    except (KeyError, TypeError, ValueError, IndexError) as error:
        _gate(gates, "chart_affine_normalization_failed", chart=chart, error=str(error))
        return {
            "chart": chart,
            "field": "m",
            "c2_handles": dict(c2),
            "c5_handles": dict(c5),
            "offset_derivative_proved_on_both": False,
            "actual_c2_J_equals_actual_c5_J": False,
            "memory_width_identity_proved": False,
            "cross_archive_kernel_identity_proved": False,
            "global_offset_identity_proved": False,
            "kernel_forms": {"C2": "exp[-rate*(right_offset-offset(t))]", "C5": "exp[-rate*J*(coordinate-t)]"},
            "admitted": False,
            "error": str(error),
        }


def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Derive exact affine kernel identities from saved source and graph data."""

    gates: list[dict[str, Any]] = []
    node_list_schema_ok = type(nodes) is list and all(isinstance(row, Mapping) for row in nodes)
    if not node_list_schema_ok:
        _gate(gates, "graph_node_list_schema_changed", expected="list[Mapping]",
              actual_type=type(nodes).__name__)
    validated_nodes: Sequence[Any] = nodes if node_list_schema_ok else []
    try:
        graph_digest = _digest(nodes) if node_list_schema_ok else ""
    except (TypeError, ValueError) as error:
        graph_digest = ""
        _gate(gates, "graph_node_list_digest_unavailable", error=str(error))
    graph_ok = node_list_schema_ok and len(nodes) == EXPECTED_NODE_COUNT and graph_digest == EXPECTED_NODES_DIGEST
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "expected_graph_digest_changed", expected=EXPECTED_NODES_DIGEST, actual=expected_nodes_digest)
    if not graph_ok:
        _gate(gates, "canonical_graph_changed", expected=EXPECTED_NODES_DIGEST, actual=graph_digest, node_count=len(nodes) if isinstance(nodes, Sequence) else None)

    records, source_meta = _source_evidence(source_records, gates)
    logc_projection = _check_logc_source_projection(records, gates)
    graph_parameter = _graph_parameter_source_guard(validated_nodes, gates)
    c2_descriptor = source_meta.get("C2_right_offset_source_descriptor")
    c2_descriptor_ok = isinstance(c2_descriptor, Mapping) and c2_descriptor.get("passed") is True
    expected_graph_digest_ok = expected_nodes_digest == EXPECTED_NODES_DIGEST
    # This flag is computed before chart arithmetic so it represents only the
    # authenticated source/graph preconditions.  A single changed pin,
    # schema, node attribute, or projection assignment keeps every row out.
    source_ok = (
        all(name in records for name in RECORD_HASHES)
        and c2_descriptor_ok
        and logc_projection["passed"]
        and graph_parameter["passed"]
        and expected_graph_digest_ok
        and graph_ok
        and not gates
    )
    logc_alias_ok = logc_projection["passed"] and graph_parameter["passed"]
    ra_ok = _check_ra_recipe(validated_nodes, gates) if graph_ok else False
    rows = [_row_result(validated_nodes, chart, source_ok, logc_alias_ok, ra_ok, gates) for chart in CHARTS]
    admitted_rows = [row["chart"] for row in rows if row.get("cross_archive_kernel_identity_proved") is True]
    unresolved_rows = [row["chart"] for row in rows if row.get("cross_archive_kernel_identity_proved") is not True]

    # This is the only expected open source identity in the frozen graph.  Do
    # not replace it with an interval/name/owner equality.
    reference_row = next(row for row in rows if row["chart"] == "reference")
    if reference_row.get("actual_c2_J_equals_actual_c5_J") is not True:
        _gate(gates, "reference_logC_source_leaf_alias_unbound",
              required="an accepted source identity tying C2 logC handle 6 to C5 true_core_logC handle 35560")

    # Keep structural source diagnostics separate from the mathematical scope.
    source_gate_names = {"source_record_missing", "source_record_pin_changed", "source_text_hash_changed", "gzip_raw_hash_changed",
                         "gzip_decoded_descriptor_missing", "gzip_decode_failed", "gzip_decoded_text_changed", "source_ast_invalid",
                         "source_declaration_missing", "source_formula_token_missing", "receipt_json_invalid", "receipt_not_accepted",
                         "C2_right_offset_descriptor_changed"}
    structural_source_ok = not any(item["gate"] in source_gate_names for item in gates)
    global_offset_identity_proved = all(row.get("global_offset_identity_proved") is True for row in rows)
    all_nine_rows_structurally_supported = (
        len(admitted_rows) == len(CHARTS)
        and not unresolved_rows
        and global_offset_identity_proved
        and logc_projection["passed"]
        and graph_parameter["passed"]
        and not gates
    )
    # Keep this admission predicate deliberately conjunctive.  In particular,
    # an incorrect caller-supplied digest or any diagnostic gate can never be
    # masked by a coincidentally complete row list.
    passed = (
        expected_nodes_digest == EXPECTED_NODES_DIGEST
        and expected_graph_digest_ok
        and graph_ok
        and structural_source_ok
        and c2_descriptor_ok
        and logc_projection["passed"]
        and graph_parameter["passed"]
        and len(admitted_rows) == len(CHARTS)
        and not unresolved_rows
        and global_offset_identity_proved
        and not gates
    )
    return {
        "schema": SCHEMA,
        "source_only": True,
        "nodes_digest": graph_digest,
        "expected_nodes_digest": expected_nodes_digest,
        "node_count": len(nodes) if isinstance(nodes, Sequence) else None,
        "passed": passed,
        "unresolved_gates": gates,
        "rows": rows,
        "admitted_rows": admitted_rows,
        "unresolved_rows": unresolved_rows,
        "offset_derivative_equals_jacobian": all(row.get("offset_derivative_proved_on_both") is True for row in rows),
        "cross_archive_kernel_equality": len(admitted_rows) == len(CHARTS),
        "cross_archive_kernel_identity_proved": len(admitted_rows) == len(CHARTS),
        "global_offset_identity": global_offset_identity_proved,
        "global_offset_identity_proved": global_offset_identity_proved,
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "physical_Q_y_Z_seam": False,
        "source_parameter_aliases": {"logP": 5, "logC": 6, "T": 7, "hbB": 8, "hbS": 9, "sc": 10, "true_core_logC": 35560},
        "ra_recipe": {"Ra": 65604, "logLambda": 35562, "four_logP": 14, "formula": "4*exp(-(4*logP+1000))", "proved": ra_ok},
        "source_evidence": source_meta,
        "logC_source_projection_guard": logc_projection,
        "logC_graph_parameter_source_guard": graph_parameter,
        "all_nine_rows_structurally_supported": all_nine_rows_structurally_supported,
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": {name: hashlib.sha256(row["decoded"]).hexdigest() for name, row in records.items()},
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "operator_calls": [],
        "producer_provider_builder_checker_replayed": False,
        "full_goal_complete": False,
    }

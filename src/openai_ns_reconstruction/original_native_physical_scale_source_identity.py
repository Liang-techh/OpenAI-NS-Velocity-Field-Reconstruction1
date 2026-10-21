"""Bind exact physical scalar sources and contract the saved terminal defect.

Only frozen source syntax, literal graph rows and saved source certificates are
consumed. No archived constructor, parameter evaluator or normalizer is run.
All aliases concern underlying exact source functions, not interval values.
"""
from __future__ import annotations

import ast
from fractions import Fraction
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_'
ORIG = PREFIX + 'current_original_'
PARENT = ORIG + 'native_physical_terminal_log_repair.json'
PHYSICAL = ORIG + 'native_physical_pulse_source.json'
MASTER = ORIG + 'Rp_native_C5_pulse_source.json.gz'
OWN = ORIG + 'n1_pulse_own_fields.json.gz'
QUIET_CHECK = ORIG + 'true_quiet_Rp_C5_pressure_C3_functions_check.json'
QUIET_SOURCE = ORIG + 'true_quiet_Rp_C5_pressure_C3_functions.py'
QUIET_CHECK_SOURCE = ORIG + 'true_quiet_Rp_C5_pressure_C3_functions_check.py'
CORE_CHECK = ORIG + 'n1_core_fixed_point_functions_check.json'
LOCAL_CHECK = ORIG + 'n1_outer_own_source_joins_construction_check.json'
NORM = PREFIX + 'physical_norm_family.json'
PRESSURE = PREFIX + 'pressure_source.json'
TRANSFER = PREFIX + 'core_transfer.json'
PROVIDER = ORIG + 'n1_all_chart_physical_pressure_fields.py'
PATCH = ORIG + 'Rm_defect_patch_inverse.py'
MACRO = ORIG + 'bridge_macro_functions.py'
WHOLE = ORIG + 'whole_Z_bridge_source.py'
SCALED = PREFIX + 'current_native_signed_input_enclosures.py'
FRAME = ORIG + 'C3_Rp_native_frame.json.gz'
FRAME_CHECK = ORIG + 'C3_Rp_native_frame_check.json'
FRAME_SOURCE = ORIG + 'C3_Rp_native_frame.py'
FRAME_CHECK_SOURCE = ORIG + 'C3_Rp_native_frame_check.py'
PRE = ORIG + 'n1_preRp_moment_functions.json.gz'
PRE_SOURCE = ORIG + 'n1_preRp_moment_functions.py'
PRIMITIVE = ORIG + 'n1_quiet_pressure_primitives.json'
PRIMITIVE_CHECK = ORIG + 'n1_quiet_pressure_primitives_check.json'
PRIMITIVE_SOURCE = ORIG + 'n1_quiet_pressure_primitives.py'
SCHEMA = 'original_native_physical_scale_source_identity_v1'
STAGE_FLAGS = (
    'native_logC_original_control_identity_proved',
    'native_original_absolute_Rc_Rp_identity_proved',
    'old_preRp_Pstar_source_identity_proved',
    'native_original_quiet_K_identity_proved',
    'physical_quiet_Z0_source_identity_proved',
    'current_physical_complete_pulse_C3_identity_proved',
    'terminal_difference_reduced_to_active_compact',
)
FALSE_FLAGS = (
    'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1',
    'global_five_moment_closure', 'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed',
    'full_goal_complete', 'old212_adapter_source_identity_proved',
    'independent_terminal_function_match_complete',
    'native_vs_registered_preRp_m15_identity_proved',
)
BASE_GRAPH = '5cdd0c8fd6a758b9e176bf67bd29ea2f576c4ddae2c7cb63d438ea1ca8388c1e'
PHYSICAL_GRAPH = '29d2257340a9bb9399448dfe0ab6bd2cd31e1815d3646741a242dac506a0a84c'
CORRECTED_GRAPH = 'bb8c47eb39b4bfce588a7ffe4311d07acf23f7098d8f5d49ee07859ea2525834'
QUIET_GRAPH = '15d51cd83dd82f41eba990979ace934bd5c8d5c339cd516b9fc96bc20a64af56'
RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_native_physical_terminal_log_repair.json': '9566d63baaf92856b8210e1ce07824a5b12c97b6af7a641bfa4c5cab36f5f040', 'lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json': '8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b', 'lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz': '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe', 'lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.json.gz': '2c0f83a9ac36f4ea7970b7938bac31d36af7c19df03e9c1584a6a11aee14713b', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions_check.json': '3cb8b2e7c48ee62ce2b3dba9c6b0d68b901010d6a154b18d466d3221e2c8e7c0', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.py': '0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions_check.py': 'f8f9b12bb3c157810d4fc124274c150fcae29ade80db1c676be4233af9b7d57d', 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json': '45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead', 'lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins_construction_check.json': 'b67073446ac25da8736cae9c09ef963cafe3474f730a3c4d709f4cd47b6b83a6', 'lei_ren_part1_paper_compliant_physical_norm_family.json': 'e8b97192b8d026a8fa2b22848bb07bae38bc876261d909fb70419b4898efb7e5', 'lei_ren_part1_paper_compliant_pressure_source.json': '936a30c1d5dd8df266cdf4464e5fe8cd542dbec240f593caae58bb993eea4f9e', 'lei_ren_part1_paper_compliant_core_transfer.json': '4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497', 'lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields.py': 'df9d282e3e1b5c7c74d5101c6e37046b6e28138d40b08ca05024dbe790c59f28', 'lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse.py': 'c88ac849f55f64aa1b338b0eb91c55d561006e495e2256f848a27ade07e379df', 'lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py': '74511ee00a655184021b0516478cffb7c8fdca830b5dc01ab44dfde2e9ad648b', 'lei_ren_part1_paper_compliant_current_native_signed_input_enclosures.py': 'e5e322e5194e5b52cd6e39be79d38a8561265f0d024532bdf020b9003ef428b0', 'lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source.py': 'df8b204d023c1d8762da16d07190601ebd1e11447077b01a734bba8baa595a6a', 'lei_ren_part1_paper_compliant_current_original_C2_phase_primitives.py': '84f5da3eb64096f82df88dd828b22b09bce18192a03925d91cfee0b193902d2c', 'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_source_functions.py': 'f0f6058bf77e52d9f5f04db3259ca9f737dc0fdb9281e2650b7a6ae5ca3f25fe', 'lei_ren_part1_paper_compliant_current_original_whole_Z_signed_loop_functions.py': '2d534e7108cbd518f0938aa50e9ba802c10c2b0a68dc65c70b2e05eec2508779', 'lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py': '83e5a8bdae75d20a261d24dd22d30cc7cb5e603e8fa2f0133958ca880b5e76db', 'lei_ren_part1_paper_compliant_current_original_whole_Z_O2_axial_buffer_finite_N.py': '2a832a91a5d67247ff17f18a3e6a46cfcfe61298554390f8ab3ecff118a10f5c', 'lei_ren_part1_paper_compliant_current_original_whole_Z_O2_slope_finite_N.py': 'dee791efd33a8976409a5f478f50ce9519a6a62390316ef3099d0ebb963863f8', 'lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_reference_finite_N.py': '137c196b18ae71047ca32ba4d8c2fe4f48c53d61c97207d4798eb58f1b1487ed', 'lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_finite_N.py': '587c279ccc6337c7a6c7805d67e3075232630dc1563fad1f91c463052d521f44', 'lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_mixed4.py': 'b8be390a9fef7af9df83cee84a84680fdddfe731d0365d79d4cdf39c9c0688d3', 'lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_centered_inverse.py': 'b0044165763821001978135ff6601a93930a2660165d507721d23a9d9bc15a96', 'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py': 'd71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59', 'lei_ren_part1_paper_compliant_core_physical_field.py': '19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54', 'lei_ren_part1_paper_compliant_current_original_whole_Z_long_Rm_finite_N.py': 'dd5830925de6befe08254e4954cb25e4a3ec1d5dec448cf6684e3c6fade87c2d', 'lei_ren_part1_paper_compliant_current_original_reshape_terminal_kernels.py': 'a846be160da9a9e4e14296d5a47da235ce899101191ca23f1da7290b4888b5bc', 'lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.json.gz': '859b59d422104d184092bcf60f8afa95267335ce87e8b41f2c271c18130e9a28', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce', 'lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame_check.json': 'd34d0e4f90b131e68bbde16caad33029ff89b75174741ad02d03d0a6e135d9e8', 'lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_primitives.json': 'e77b0e182b18b261e145db5c0ecbbb7652f1e1219170b2971c4028d9b02a7852', 'lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_primitives_check.json': '8316c8bf06f3f524ef6c57aa7313f4b5546399f5bd59f2bc209a806cf71bdf0d', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py': '661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e', 'lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_primitives.py': '8e024c3c1047c60e08bc249c99211d407c2865cb1d7b6de377ab5408e42db732', 'lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.py': 'cb7bc1bc1801e04c198f40458dd88246434e74318f500260e8a9043c49398674', 'lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame_check.py': '5baaac7050500641aa1a7f1481f2fc3a2994af9d4566bcb96956f5bfd8ebd06c'}
DECODED_RECORD_HASHES = {MASTER: '56137195254ae71d7eedfccc82ac7b5c048c8f3acd8de87fedfee5c6725fd7b1', OWN: '8cfee4eb2ad845c98df7cb8e7dafac6a32df473c56f02eef56751902293715a3', FRAME: '1d563553e3c7e6f9115a024a52e61f3b308d822104aaf557f3d235c54f26738d', PRE: '348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64'}

# These are source-level ports.  The integer values are graph handles, never
# evaluated scalar values.  In particular, the Pstar port is the provider's
# function projection and the quiet K port is a source identity certificate.
SOURCE_ALIASES = {
    'logC': {
        'source_record': NORM,
        'json_pointer': '/selected_logCstar',
        'native_node': 35560,
        'old_node': 6,
        'identity_scope': 'selected source scalar from the canonical datum',
        'interval_overlap_is_not_identity': True,
    },
    'Rc': {
        'source_graph_digest': BASE_GRAPH,
        'native_node': 251066,
        'old_node': 30828,
        'identity_scope': 'absolute-radius log recipe only',
        'adapter_alias': False,
    },
    'Rp': {
        'source_graph_digest': BASE_GRAPH,
        'native_node': 263349,
        'old_node': 31162,
        'identity_scope': 'absolute-radius log recipe only',
        'adapter_alias': False,
    },
    'Pstar': {
        'source_graph_digest': BASE_GRAPH,
        'provider_node': 25203,
        'underlying_formula_node': 33572,
        'raw_quantity_path': 'outer.owner.owner(ends).Pstar',
        'identity_scope': 'universal source function over admissible ends',
        'materialized_value_alias': False,
    },
    'quietK': {
        'source_graph_digest': QUIET_GRAPH,
        'native_node': 265180,
        'native_Rp_node': 263349,
        'Pstar_node': 33572,
        'old_kappa_node': 31204,
        'identity_scope': 'native quiet mass source after explicit frame bridge',
        'normalized_pulse_m1_is_not_kappa': True,
    },
}

SOURCE_IDENTITIES = {
    'pulseC3congruence': {
        'kind': 'ordinary_Z_function_identity',
        'ordinary_Z_orders': [0, 1, 2, 3],
        'formula': 'current physical complete C3 pulse equals old frame-only pulse and saved old pulse by parent bridge',
        'geometry_delta_C3': '0',
        'saved_frame_only_to_old_pulse_bridge': True,
        'source_only': True,
    },
    'quietZ0cancellation': {
        'kind': 'ordinary_Z_function_identity',
        'ordinary_Z_order': 0,
        'formula': 'native quiet primitive minus old quiet primitive equals 0',
        'primitive': '-K^2*length/4',
        'source_only': True,
    },
    'signedremainingnative2Rc-correctedactive-oldcompact': {
        'kind': 'signed_source_function_decomposition',
        'coefficients': {'native_2Rc': '1', 'corrected_active': '-1', 'old_compact': '-1'},
        'formula': 'native_2Rc - corrected_active - old_compact',
        'zero_claimed': False,
        'source_only': True,
    },
}

# External static receipts are pinned as metadata.  They are not opened by
# derive; the caller has already authenticated their raw bytes.  The linear
# receipt is the explicit cross-graph bridge needed for the quiet K identity.
STATIC_RECEIPT_HASHES = {
    'logC_exact_source_binding_static_2026_10_21.json': '68a721177f9d076a59228ff75aa0f0ddde0ffee915541f0916e953cc837b4bb2',
    'pstar_universal_source_binding_static_2026_10_21.json': '2e622656ed538a5fdf04eca293a92f019e27f651229ea9ee580375cd8f5e7d45',
    'pstar_definition_static_scan_2026_10_20.json': 'c09ab9ef4b8b02f52bc65fb877dd2c504574ea5c967ec12fbf7b17f64e821252',
    'quiet_kappa_same_source_static_2026_10_21.json': 'f9d93b8063a1c0bbec7b3eb64651a9e39abbe2526b70c7a29f38fbcd37b2e710',
    'quiet_kappa_linear_source_static_2026_10_21.json': 'eb93dcf81b702b2bacd3ea5ecd1ff38b88058732e655a970742968eb387ef09d',
}


RADIUS_LOG_DIFFERENCE = [
    {'monomial': ['native_core_logC'], 'coefficient': '10'},
    {'monomial': ['old_selected_logCstar'], 'coefficient': '-10'},
]

RETAINED_COUNTS = {
    'cumulative_terms': 212,
    'raw_signed_terms': 213,
    'reduced_terms': 183,
}

QUIET_DEFECT_GRAPH = '9a09124e300a64bc3d955ba1541ac09ae7f4b072f62ae6025a86a34794c01edb'


def _digest(value):
    """Digest JSON metadata using the runner's canonical representation."""
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
        allow_nan=False).encode('utf-8')).hexdigest()


def _dict(value):
    return value if isinstance(value, dict) else {}


def _list(value):
    return value if isinstance(value, list) else []


def _json_record(source_records, name):
    record = _dict(source_records.get(name))
    text = record.get('source_text')
    if not isinstance(text, str):
        return {}
    try:
        value = json.loads(text)
    except (TypeError, ValueError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _ast_record(source_records, name):
    record = _dict(source_records.get(name))
    text = record.get('source_text')
    if not isinstance(text, str):
        return None
    try:
        return ast.parse(text)
    except (SyntaxError, ValueError, TypeError):
        return None


def _function(tree, name):
    if tree is None:
        return None
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    return None


def _method(tree, class_name, method_name):
    if tree is None:
        return None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == method_name:
                    return child
    return None


def _ast_text(node):
    if node is None:
        return ''
    try:
        return ast.unparse(node)
    except (AttributeError, ValueError):
        return ''


def _ast_has(node, *needles):
    text = _ast_text(node)
    return all(needle in text for needle in needles)


def _node(nodes, handle):
    if not isinstance(handle, int) or handle < 0 or handle >= len(nodes):
        return None
    value = nodes[handle]
    return value if isinstance(value, dict) else None


def _field(obj, *keys, default=None):
    value = obj
    for key in keys:
        if not isinstance(value, dict):
            return default
        value = value.get(key, default)
    return value


def _all_true(mapping, keys):
    return all(mapping.get(key) is True for key in keys)


class _SourceCertificateGraph:
    """A small exact integer certificate graph over source-function ports."""

    def __init__(self):
        self.nodes = []
        self.forms = {}
        self.ports = {}

    def append(self, operation, arguments=(), form=None, **metadata):
        ref = SCHEMA + ':new' + str(len(self.nodes))
        row = {'id': ref, 'operation': operation, 'arguments': list(arguments)}
        row.update(metadata)
        self.nodes.append(row)
        if form is None:
            self.forms[ref] = {ref: Fraction(1)}
        else:
            self.forms[ref] = {key: value for key, value in form.items() if value}
        return ref

    def port(self, namespace, source_node, order, source_graph_digest):
        key = (namespace, str(source_node), int(order))
        if key not in self.ports:
            self.ports[key] = self.append(
                'source_function_port', [],
                namespace=namespace, source_node=source_node,
                source_graph_digest=source_graph_digest,
                ordinary_Z_order=int(order), bound_or_midpoint=False,
                physical_m15=True)
        return self.ports[key]

    def signed(self, terms, **metadata):
        form = {}
        arguments = []
        coefficients = []
        for coefficient, ref in terms:
            arguments.append(ref)
            coefficients.append(str(Fraction(coefficient)))
            for atom, value in self.forms.get(ref, {}).items():
                form[atom] = form.get(atom, Fraction(0)) + Fraction(coefficient) * value
        return self.append('signed_source_function_sum', arguments, form,
                           coefficients=coefficients, **metadata)

    def identity(self, source, target, admitted, identity, **metadata):
        form = self.forms.get(target, {}) if admitted else self.forms.get(source, {})
        return self.append(
            'qualified_source_function_identity', [source, target], form,
            admitted=bool(admitted), identity=identity, **metadata)

    def compose_identity(self, source, middle, target, admitted, identity, **metadata):
        """Record a source identity followed by an already saved alias."""
        form = self.forms.get(target, {}) if admitted else self.forms.get(source, {})
        return self.append(
            'qualified_source_function_identity_composition',
            [source, middle, target], form,
            admitted=bool(admitted), identity=identity, **metadata)

    def parent_atom(self, parent_atom_id, source_ref, coefficient, **metadata):
        """Import one canonical parent atom while retaining its source role."""
        return self.append(
            'imported_saved_parent_canonical_atom', [source_ref],
            self.forms.get(source_ref, {}),
            parent_atom=parent_atom_id, parent_coefficient=str(coefficient),
            **metadata)

    def contract_parent(self, parent_atoms, substitutions, admitted, **metadata):
        """Contract saved parent atoms through qualified source substitutions."""
        form = {}
        arguments = list(parent_atoms)
        for coefficient, ref in substitutions:
            arguments.append(ref)
            for atom, value in self.forms.get(ref, {}).items():
                form[atom] = form.get(atom, Fraction(0)) + Fraction(coefficient) * value
        return self.append(
            'contract_saved_parent_canonical_terminal', arguments, form,
            admitted=bool(admitted), substitutions=[
                {'coefficient': str(Fraction(coefficient)), 'ref': ref}
                for coefficient, ref in substitutions], **metadata)

    def cancel(self, source, admitted, identity, **metadata):
        form = {} if admitted else self.forms.get(source, {})
        return self.append(
            'qualified_source_function_cancellation', [source], form,
            admitted=bool(admitted), identity=identity, **metadata)

    def normal(self, ref):
        return {
            atom: str(value)
            for atom, value in sorted(self.forms.get(ref, {}).items())
            if value
        }


def derive(source_records, *, expected_input_digest):
    """Derive the source-only physical scale certificate once.

    Every operation below consumes authenticated JSON text, literal graph rows,
    or AST syntax.  No provider, constructor, checker, normalizer, evaluator,
    quadrature routine, or archived graph reducer is imported or called.
    """
    records = source_records if isinstance(source_records, dict) else {}
    checks = {}

    def check(key, condition):
        checks[key] = bool(condition)
        return checks[key]

    check('input_digest', _digest(records) == expected_input_digest)

    authenticated = True
    for name, expected in RECORD_HASHES.items():
        record = _dict(records.get(name))
        decoded_expected = DECODED_RECORD_HASHES.get(name, expected)
        valid = (
            record.get('raw_authenticated_before_decode') is True
            and record.get('raw_sha256') == expected
            and record.get('decoded_sha256') == decoded_expected
            and isinstance(record.get('source_text'), str)
            and hashlib.sha256(record.get('source_text', '').encode('utf-8')).hexdigest() == decoded_expected
        )
        authenticated = authenticated and valid
        check('authenticated:' + name, valid)
    check('all_31_pinned_inputs_authenticated', authenticated)
    check('all_pinned_inputs_authenticated', authenticated)

    data = {}
    trees = {}
    for name in RECORD_HASHES:
        if name.endswith('.json') or name.endswith('.json.gz'):
            data[name] = _json_record(records, name)
        elif name.endswith('.py'):
            trees[name] = _ast_record(records, name)

    parent = _dict(data.get(PARENT))
    parent_component = _dict(parent.get('component_result'))
    physical = _dict(data.get(PHYSICAL))
    physical_component = _dict(physical.get('component_result'))
    master = _dict(data.get(MASTER))
    master_nodes = _list(master.get('exact_graph_nodes'))
    own = _dict(data.get(OWN))
    own_nodes = _list(own.get('exact_graph_nodes'))
    quiet_check = _dict(data.get(QUIET_CHECK))
    quiet_checks = _dict(quiet_check.get('checks'))
    core_check = _dict(data.get(CORE_CHECK))
    core_checks = _dict(core_check.get('checks'))
    local_check = _dict(data.get(LOCAL_CHECK))
    local_checks = _dict(local_check.get('checks'))
    norm = _dict(data.get(NORM))
    pressure = _dict(data.get(PRESSURE))
    transfer = _dict(data.get(TRANSFER))
    frame = _dict(data.get(FRAME))
    frame_check = _dict(data.get(FRAME_CHECK))
    pre = _dict(data.get(PRE))
    primitive = _dict(data.get(PRIMITIVE))
    primitive_check = _dict(data.get(PRIMITIVE_CHECK))
    pre_tree = trees.get(PRE_SOURCE)
    primitive_tree = trees.get(PRIMITIVE_SOURCE)
    primitive_init_fn = _method(
        primitive_tree, 'CurrentOriginalN1QuietPressurePrimitives', '__init__')

    family = parent.get('source_family')
    if not isinstance(family, dict):
        family = physical.get('source_family')
    if not isinstance(family, dict):
        family = master.get('source_family')
    if not isinstance(family, dict):
        family = {}
    family_records = (parent, physical, master, own, quiet_check, core_check, local_check)
    family_present = [row.get('source_family') for row in family_records
                      if isinstance(row.get('source_family'), dict)]
    check('same_source_family', bool(family) and all(item == family for item in family_present))
    check('source_family_has_actual_defect_implicit_and_datum', set(family) >= {
        'actual_five_defect_family_sha256', 'implicit_source_sha256', 'datum_enclosure_sha256'})
    check('source_family_pressure_and_transfer_match',
          pressure.get('compliant_source', {}).get('implicit_source_sha256') == family.get('implicit_source_sha256')
          and pressure.get('compliant_source', {}).get('datum_enclosure_sha256') == family.get('datum_enclosure_sha256')
          and transfer.get('implicit_source_sha256') == family.get('implicit_source_sha256')
          and transfer.get('datum_enclosure_sha256') == family.get('datum_enclosure_sha256'))

    retained = parent.get('retained_construction_data')
    if not isinstance(retained, dict):
        retained = parent_component.get('retained_construction_data')
    retained = retained if isinstance(retained, dict) else {}
    check('retained_construction_counts_212_213_183',
          all(retained.get(key) == value for key, value in RETAINED_COUNTS.items()))
    check('retained_old_error_data_and_channels',
          retained.get('six_signed_channels_retained') is True
          and retained.get('all_old_error_data_retained') is True
          and retained.get('inherited_errors_deleted') is False)
    check('reserved_contribution_registered_once',
          retained.get('reserved_contribution_registered_once') is True
          and retained.get('accepted_physical_baseline_commit') == '8f36f19e3c524b424a9286c46cc8e0b855f75307')
    check('legacy_numeric_chain_unchanged',
          parent_component.get('legacy_numeric_chain_reconfigured') is False
          and parent_component.get('old_error_terms_removed') is False
          and parent_component.get('extra_reserved_contribution_added') is False)

    parent_audit = _dict(parent_component.get('source_inspection_audit'))
    physical_audit = _dict(physical_component.get('source_inspection_audit'))
    check('parent_saved_source_only', parent_component.get('source_only') is True)
    check('parent_saved_no_old_or_numeric_consumption',
          parent_audit.get('old_definition_graphs_rebuilt') is False
          and parent_audit.get('old_operators_replayed') is False
          and parent_audit.get('numeric_values_computed') is False
          and parent_audit.get('module_imports_of_ancestor_constructors') == 0)
    check('physical_saved_source_only', physical.get('source_only') is True)
    check('physical_saved_no_old_or_numeric_consumption',
          physical_audit.get('old_operators_replayed') is False
          and physical_audit.get('numeric_values_computed') is False
          and physical_audit.get('accepted_ancestor_operator_calls') == 0)

    # LOGC: consume the saved local projection and the saved core oracle.  The
    # local report's full seam acceptance remains deliberately out of scope.
    selected_projection = _dict(_field(local_checks, 'actual_old_current_radius_sources',
                                       'selected_parameter_source_projection'))
    core_oracle = _dict(_field(core_checks, 'actual_oracle'))
    core_parameter = _dict(core_oracle.get('exact_parameter_and_native_P0_identity'))
    core_source_identity = _dict(core_parameter.get('actual_source_scalar_and_unique_anchor_identity'))
    logc_projection = (
        selected_projection.get('passed') is True
        and selected_projection.get('old_parameter_node') == 6
        and selected_projection.get('true_parameter_node') == 35560
        and selected_projection.get('canonical_datum') == NORM
        and selected_projection.get('canonical_datum_sha256') == RECORD_HASHES[NORM]
        and selected_projection.get('canonical_datum_key') == 'selected_logCstar'
        and selected_projection.get('accepted_true_core_parameter_identity_consumed') is True
        and selected_projection.get('same_selected_scalar_recipe_not_interval_overlap_or_new_common_radius') is True
    )
    check('local_saved_selected_logC_source_projection', logc_projection)
    check('local_full_parent_seam_acceptance_not_consumed',
          local_check.get('full_seam_acceptance_passed') is False)
    check('core_saved_parameter_and_native_P0_oracle',
          core_oracle.get('passed') is True
          and core_parameter.get('passed') is True
          and core_parameter.get('selected_j_logC_and_anchor_are_from_this_actual_core') is True
          and core_parameter.get('interval_overlap_not_used_as_a_function_identity') is True)
    check('core_saved_logC_and_anchor_exact_source',
          core_source_identity.get('logC_matches_original_selected_source') is True
          and core_source_identity.get('logC_leaf_uses_this_exact_source') is True
          and core_source_identity.get('j_leaf_uses_this_exact_source') is True
          and core_source_identity.get('anchor_leaf_uses_this_exact_source') is True)
    logc_stage = all(checks[key] for key in (
        'local_saved_selected_logC_source_projection',
        'local_full_parent_seam_acceptance_not_consumed',
        'core_saved_parameter_and_native_P0_oracle',
        'core_saved_logC_and_anchor_exact_source'))
    check('native_logC_source_identity_complete', logc_stage)

    # Absolute radius aliases use the exact saved log recipes.  Independent
    # interval endpoints are intentionally never compared here.
    radius = _dict(parent_component.get('absolute_radius_source_audit'))
    check('absolute_radius_saved_master_graph',
          radius.get('source_archive') == MASTER
          and radius.get('source_graph_digest') == BASE_GRAPH
          and radius.get('native_Rc') == 251066
          and radius.get('old_Rc') == 30828
          and radius.get('native_Rp') == 263349
          and radius.get('old_Rp') == 31162)
    check('absolute_radius_exact_logRc_difference',
          radius.get('logRc_difference') == RADIUS_LOG_DIFFERENCE)
    check('absolute_radius_exact_logRp_difference',
          radius.get('logRp_difference') == RADIUS_LOG_DIFFERENCE)
    check('absolute_radius_quiet_length_difference_empty', radius.get('quiet_length_difference') == [])
    check('absolute_radius_saved_endpoint_and_positive_scope',
          radius.get('substituted_endpoint_handle') == 48
          and radius.get('substituted_endpoint_exact_value') == '2'
          and radius.get('exact_positive_radius_log_laws_only') is True
          and radius.get('source_values_or_bounds_evaluated') is False
          and radius.get('old_math_consumers_called') == 0)
    radius_stage = logc_stage and all(checks[key] for key in (
        'absolute_radius_saved_master_graph',
        'absolute_radius_exact_logRc_difference',
        'absolute_radius_exact_logRp_difference',
        'absolute_radius_quiet_length_difference_empty',
        'absolute_radius_saved_endpoint_and_positive_scope'))
    check('native_absolute_Rc_Rp_source_identity_complete', radius_stage)

    # Pstar: inspect the provider chain and the exact graph recipe.  This
    # proves a function projection for every admissible ends tuple, rather than
    # identifying a runtime enclosure object or a huge materialized value.
    provider_tree = trees.get(PROVIDER)
    patch_tree = trees.get(PATCH)
    macro_tree = trees.get(MACRO)
    whole_tree = trees.get(WHOLE)
    scaled_tree = trees.get(SCALED)
    provider_fn = _function(provider_tree, 'physical_source_scale')
    patch_init = _method(patch_tree, '_ActualRmLeadingPatch', '__init__')
    macro_init = _method(macro_tree, 'MacroFlow', '__init__')
    macro_factor = _method(macro_tree, 'MacroFlow', 'factor')
    whole_init = _method(whole_tree, 'OriginalWholeZBridgeSource', '__init__')
    formal_init = _method(scaled_tree, 'FormalScale', '__init__')
    scaled_init = _method(scaled_tree, 'ScaledEnclosure', '__init__')
    formal_evaluate = _method(scaled_tree, 'FormalScale', 'evaluate')
    finite_interval = _method(scaled_tree, 'ScaledEnclosure', 'finite_interval')
    check('Pstar_provider_ast_projection',
          _ast_has(provider_fn, 'owner.modified.ranges.phase.outer.owner.owner(tuple(ends))',
                   "_same_source(op.Pstar, owner.c, op.flow, 'Pstar')",
                   'return op.Pstar')
          and _ast_has(_function(provider_tree, '_same_source'),
                       'value.ctx is not c',
                       'value.scale.bases is not flow.logs',
                       'value.ledger is not flow.ledger', 'raise ValueError'))
    check('Pstar_patch_ast_half_power',
          _ast_has(patch_init, 'self.flow = f = reference.flow',
                   'self.Pstar = f.factor((0, 0.5, 0, 0, 0))'))
    check('Pstar_macro_ast_basis_and_factor',
          _ast_has(macro_init, 'self.logs = tuple((c.mpf(v) for v in (logh, logP2, logF02, logRa, 0)))')
          and _ast_has(macro_factor,
                       'prior.FormalScale(self.logs, powers, offset)',
                       'prior.ScaledEnclosure', 'offset=0'))
    check('Pstar_formal_and_scaled_constructors_preserve_source',
          _ast_has(formal_init, 'self.bases = bases', 'self.powers = tuple(powers)',
                   'self.offset = self.ctx.mpf(offset)')
          and _ast_has(scaled_init, 'self.coefficient = self.ctx.mpf(coefficient)',
                       'self.scale = scale', 'self.ledger = ledger', 'if magnitude',
                       'self.scale = FormalScale')
          and 'prior.ScaledEnclosure(prior.FormalScale(self.logs, powers, offset), 1, self.ledger)' in
              _ast_text(macro_factor))
    check('Pstar_whole_bridge_ast_core_transfer',
          _ast_has(whole_init, "logP=scalar(major, 'logPstar')")
          and _ast_has(_function(whole_tree, 'owner'),
                       'macro.MacroFlow(c, self.logh, 2 * self.core.logP,')
          and core_parameter.get('exact_graph_logP_delta_logLambda_epsilon_recipes_verified') is True
          and _field(core_parameter, 'actual_core_encloses_same_exact_recipes', 'logP') is True)
    check('Pstar_scaled_enclosure_ast_is_formal',
          _ast_has(formal_evaluate, 'return sum((b * p for b, p in zip(self.bases, self.powers) if p), self.offset)')
          and _ast_has(finite_interval, 'self.coefficient * self.ctx.exp(self.scale.evaluate())'))

    implicit_definition = _dict(_field(pressure, 'compliant_source', 'implicit_source_definition'))
    check('Pstar_pressure_source_recipe',
          implicit_definition.get('Md') == '40'
          and implicit_definition.get('logPstar') == 'exp(Md)+11')
    check('Pstar_pressure_and_transfer_source_family',
          _field(pressure, 'compliant_source', 'implicit_source_sha256') == family.get('implicit_source_sha256')
          and _field(pressure, 'compliant_source', 'datum_enclosure_sha256') == family.get('datum_enclosure_sha256')
          and transfer.get('logPstar') is not None)

    master_pstar = _node(master_nodes, 25203)
    master_md = _node(master_nodes, 2)
    master_exp_md = _node(master_nodes, 3)
    master_plus = _node(master_nodes, 4)
    master_logp = _node(master_nodes, 5)
    master_pstar_formula = _node(master_nodes, 33572)
    check('Pstar_master_graph_exact_recipe',
          master_md == {'operation': 'exact_rational', 'numerator': 40, 'denominator': 1}
          and master_exp_md == {'operation': 'analytic_unary', 'name': 'exp', 'argument': 2}
          and master_plus == {'operation': 'exact_rational', 'numerator': 11, 'denominator': 1}
          and master_logp == {'operation': 'sum', 'arguments': [3, 4]}
          and master_pstar_formula == {'operation': 'analytic_unary', 'name': 'exp', 'argument': 5})
    check('Pstar_master_provider_function_port',
          isinstance(master_pstar, dict)
          and master_pstar.get('operation') == 'current_original_physical_source_scale'
          and master_pstar.get('name') == 'Pstar'
          and master_pstar.get('raw_quantity_path') == 'outer.owner.owner(ends).Pstar'
          and master_pstar.get('Z_independent') is True
          and master_pstar.get('physical_common_unit') is True
          and master_pstar.get('defining_quantity_not_a_range_value') is True
          and master_pstar.get('same_live_basis_context_and_ledger') is True
          and _field(master_pstar, 'source_rows_provider', 'module') == PROVIDER
          and _field(master_pstar, 'source_rows_provider', 'module_sha256') == RECORD_HASHES[PROVIDER]
          and _field(master_pstar, 'source_rows_provider', 'function') == 'physical_source_scale'
          and _field(master_pstar, 'source_rows_provider', 'AST_sha256') ==
              '55164d3a3db57726ae947afbc5fb4f97ed9f2c656a1ea57e38289e115c50a9a2'
          and _field(master_pstar, 'source_projection_binding', 'module') == PROVIDER
          and _field(master_pstar, 'source_projection_binding', 'module_sha256') == RECORD_HASHES[PROVIDER]
          and _field(master_pstar, 'source_projection_binding', 'AST_sha256') ==
              '55164d3a3db57726ae947afbc5fb4f97ed9f2c656a1ea57e38289e115c50a9a2')
    physical_scale = _dict(_field(physical_component, 'physical_scale_receipt'))
    check('Pstar_saved_physical_receipt',
          physical_scale.get('Pstar') == 33572
          and physical_scale.get('old_pulse_Pstar') == 33572
          and physical_scale.get('common_exact_expression_port') is True
          and physical_scale.get('extra_Pstar_factor_applied') is False
          and physical_scale.get('old_preRp_physical_scale_projection_25203_not_aliased') is True
          and physical_scale.get('all_three_charts_and_finite_endpoints_retained') is True
          and physical_scale.get('m15_paper_sign_retained') is True
          and physical_scale.get('target_negative_applied_once') is True
          and physical_scale.get('positive_actual_Rp_and_original_common_source_domain_required') is True)
    pstar_stage = all(checks[key] for key in (
        'Pstar_provider_ast_projection',
        'Pstar_patch_ast_half_power',
        'Pstar_macro_ast_basis_and_factor',
        'Pstar_whole_bridge_ast_core_transfer',
        'Pstar_scaled_enclosure_ast_is_formal',
        'Pstar_pressure_source_recipe',
        'Pstar_pressure_and_transfer_source_family',
        'Pstar_master_graph_exact_recipe',
        'Pstar_master_provider_function_port',
        'Pstar_saved_physical_receipt',
        'Pstar_formal_and_scaled_constructors_preserve_source'))
    check('old_preRp_Pstar_source_identity_complete', pstar_stage)

    # Quiet K: the saved native mass theorem is combined with the explicit
    # frame bridge.  The bridge uses raw history m[0], raw ordinary derivative
    # m[1], and the saved GraphImport 26263 -> 31204 mapping.  The normalized
    # pulse field called m1 is deliberately excluded from this identity.
    native_mass = _dict(quiet_checks.get('native_mass_identity'))
    master_k = _node(master_nodes, 265180)
    master_kappa2 = _node(master_nodes, 265179)
    master_linear_mass = _node(master_nodes, 265181)
    frame_functions = _dict(frame.get('actual_Rp_native_frame_functions'))
    frame_histories = _dict(frame_functions.get('actual_identified_Rp_histories_C3'))
    frame_m_rows = _list(frame_histories.get('m'))
    pre_quiet = _dict(_field(pre, 'defining_bindings', 'quiet'))
    pre_mapping = _dict(pre_quiet.get('original_native_frame_import_mapping'))
    primitive_identity = _dict(primitive.get('actual_source_identity'))
    primitive_flags = _dict(primitive_check.get('checks'))
    pre_input_hashes = _dict(pre.get('input_hashes'))
    primitive_input_hashes = _dict(primitive.get('input_hashes'))
    pre_quiet_fn = _function(pre_tree, 'quiet_functions')
    primitive_identity_fn = _method(
        primitive_tree, 'CurrentOriginalN1QuietPressurePrimitives', 'actual_source_identity')
    primitive_assert_fn = _method(
        primitive_tree, 'CurrentOriginalN1QuietPressurePrimitives', 'assert_source')
    check('quiet_saved_check_all_passed', quiet_check.get('all_passed') is True)
    check('quiet_saved_native_mass_identity',
          native_mass.get('passed') is True
          and native_mass.get('actual_compact_derived_K_equals_original_Rp_Pstar_kappa') is True
          and native_mass.get('actual_Rh_source_identity_and_same_original_native_mass_FunctionRef_used') is True
          and native_mass.get('common_Rp_Rc_ratio_and_finite_mu_Tw_recipe_bound') is True
          and native_mass.get('no_numeric_kernel_value_or_range_endpoint_used') is True)
    check('quiet_native_K_graph_recipe',
          master_k == {'operation': 'product', 'arguments': [262918, 265179]}
          and master_kappa2 == {'operation': 'product', 'arguments': [33572, 265178]}
          and isinstance(master_linear_mass, dict)
          and master_linear_mass.get('operation') == 'proved_original_quiet_linear_mass_source_identity'
          and master_linear_mass.get('physical_scale') == 33572
          and master_linear_mass.get('scalar_kappa2') == 265179
          and master_linear_mass.get('constant_K') == 265180
          and master_linear_mass.get('source_constant_not_a_saved_range_endpoint') is True
          and _node(master_nodes, 31205) == {
              'operation': 'product', 'arguments': [31162, 25203, 31204]})
    check('quiet_kappa_linear_frame_bridge',
          frame.get('source_family') == family
          and frame_m_rows == [26262, 26263, 26264, 26264]
          and pre.get('source_family') == family
          and pre_quiet.get('kappa') == 31204
          and pre_quiet.get('K') == 31205
          and pre_mapping.get('26263') == 31204
          and len(pre_mapping) == 12173
          and all(pre_input_hashes.get(name) == RECORD_HASHES[name] for name in (
              FRAME_SOURCE, FRAME, FRAME_CHECK_SOURCE, FRAME_CHECK,
              PRIMITIVE_SOURCE, PRIMITIVE, PRIMITIVE_CHECK))
          and all(primitive_input_hashes.get(name) == RECORD_HASHES[name] for name in (
              FRAME_SOURCE, FRAME, FRAME_CHECK_SOURCE, FRAME_CHECK,
              PRIMITIVE_SOURCE))
          and frame_check.get('all_passed') is True
          and frame_check.get('current_original_complete_C3_chain_to_Rp_native_frame_identity_installed') is True
          and _field(frame_check, 'independent_actual_complete_chain_and_Rp_frame',
                     'genuine_compact_band_and_complete_exit_consumed') is True
          and _ast_has(pre_quiet_fn, 'GraphImport',
                       "raw['actual_Rp_native_frame_functions']",
                       "f['actual_identified_Rp_histories_C3']['m'][1]",
                       'kappa = ref')
          and _ast_has(primitive_init_fn,
                       "self.mean_common = s.cancel(parsed['u'] * parsed['m1'])",
                       'self.kappa = s.cancel(self.mean_common / self.Z)'))
    check('quiet_kappa_frame_and_primitive_receipts',
          primitive.get('source_family') == family
          and primitive_identity.get('actual_frame_strings') ==
              frame.get('actual_analytic_native_inlet_functions')
          and 'Z' in primitive_identity.get('mean_common', '')
          and 'Z' not in primitive_identity.get('kappa', '')
          and primitive_identity.get('kappa_Z_derivative') == '0'
          and primitive_identity.get('q_common_plus_kappa') == '0'
          and primitive.get('same_original_endpoint_function_nodes') is not None
          and primitive_identity.get('exact_Z_independent_kappa') is True
          and primitive_identity.get('exact_q_common_at_Rp_is_minus_kappa') is True
          and primitive_check.get('all_passed') is True
          and _all_true(primitive_flags, (
              'actual_native_frame_mean_and_q_source_identities_verified',
              'same_original_endpoint_function_nodes_and_raw_frame_strings_bound',
              'formal_K_Rp_Pstar_kappa_keeps_physical_R_units'))
          and _ast_has(primitive_identity_fn, 's.simplify(s.diff(self.kappa, self.Z))',
                       'self.q_common_U0 + self.kappa')
          and _ast_has(primitive_assert_fn, 'self.frame_raw', 'self.quiet_raw',
                       'self.endpoint_nodes'))
    check('quiet_kappa_analytic_extension_from_linear_source',
          frame_m_rows[:2] == [26262, 26263]
          and 'Z' in primitive_identity.get('mean_common', '')
          and 'Z' not in primitive_identity.get('kappa', '')
          and primitive_identity.get('kappa_Z_derivative') == '0'
          and primitive_identity.get('exact_Z_independent_kappa') is True
          and _ast_has(primitive_identity_fn, 's.diff(self.kappa, self.Z)')
          and _ast_has(_function(trees.get(FRAME_CHECK_SOURCE), 'native_inlet'),
                       "M = s.cancel(H['m'] / z)",
                       'U = s.cancel(E * Q)',
                       'C1=M / U',
                       'z not in value.free_symbols',
                       "u=U / Q, m1=constants['C1'] * z * Q",
                       "zero(value - frame[key])")
          and _field(frame_check, 'independent_actual_six_row_native_inlet',
                     'current_complete_Rp_qshapes_identified_from_actual_leading_functions') is True)
    kappa_stage = radius_stage and all(checks[key] for key in (
        'quiet_saved_check_all_passed',
        'quiet_saved_native_mass_identity',
        'quiet_native_K_graph_recipe',
        'quiet_kappa_linear_frame_bridge',
        'quiet_kappa_frame_and_primitive_receipts',
        'quiet_kappa_analytic_extension_from_linear_source'))
    check('native_original_quiet_K_source_identity_complete', pstar_stage and kappa_stage)

    # Current physical pulse congruence.  The saved C5 qualification supplies
    # the frame identity; the own-field graph supplies only the radius
    # subexpressions needed to verify the three chart bindings.
    physical_prereq = _dict(physical_component.get('source_prerequisite_audit'))
    physical_prereq_keys = (
        'accepted_parent_and_open_geometry_gate',
        'authenticated_master_and_function_prefixes',
        'actual_current_frame_and_geometry',
        'whole_function_native_frame_qualification',
        'complete_C3_source_qualification',
        'common_pulse_Pstar_is_an_exact_port',
        'old_radius_ports_are_actual_expressions',
        'source_family',
        'all_physical_radius_ports_are_current',
        'old_physical_geometry_absent_outside_current_ports',
        'actual_current_radius_reaches_m15',
        'current_input_ports_unchanged',
        'geometry_binding_Z_independent',
    )
    check('physical_saved_source_prerequisites',
          _all_true(physical_prereq, physical_prereq_keys))
    physical_bindings = _dict(physical_component.get('physical_pulse_bindings'))
    master_bindings = _dict(master.get('Rp_native_C5_pulse_source_bindings'))
    current_frame = _dict(master_bindings.get('current_frame'))
    saved_current_frame = _dict(physical_bindings.get('current_frame'))
    check('physical_current_frame_matches_master',
          bool(current_frame) and current_frame == saved_current_frame)
    current_geometry = _dict(master_bindings.get('geometry'))
    saved_geometry = _dict(physical_bindings.get('geometry'))
    check('physical_geometry_matches_master',
          current_geometry == saved_geometry
          and current_geometry.get('N') == 18
          and current_geometry.get('mu') == 17
          and current_geometry.get('Pstar') == 33572
          and current_geometry.get('delta') == 25220
          and current_geometry.get('Rp') == 263349
          and current_geometry.get('R2') == 262918
          and current_geometry.get('Tw') == 205147
          and current_geometry.get('length') == 263346)
    master_coordinates = _dict(master_bindings.get('pulse_coordinates'))
    saved_coordinates = _dict(physical_bindings.get('coordinates'))
    check('physical_coordinates_match_master', master_coordinates == saved_coordinates)
    check('physical_query_y_coordinates_are_exact', all(
        _field(master_coordinates, chart, 'query') == expected_query
        and _field(master_coordinates, chart, 'y') == expected_y
        and _field(master_coordinates, chart, 'original_mu_and_actual_current_Rp_used') is True
        for chart, expected_query, expected_y in (
            ('main', 33578, 33579), ('gap', 33829, 33830), ('end', 34005, 34006))))

    own_fields = _dict(own.get('own_main_gap_end_fields'))
    own_rows_ok = True
    own_expected = {
        'main': {'query': 33578, 'y': 33579, 'log_radius': 33804, 'radius': 33805,
                 'domain_coordinate': 'xi=mu*log(R/Rp)'},
        'gap': {'query': 33829, 'y': 33830, 'log_radius': 33981, 'radius': 33982,
                'domain_coordinate': 'xi=mu*log(R/Rp)'},
        'end': {'query': 34005, 'y': 34006, 'log_radius': 34271, 'radius': 34272,
                'domain_coordinate': 's=log(R/Rv)'},
    }
    for chart, expected in own_expected.items():
        field = _dict(own_fields.get(chart))
        own_rows_ok = own_rows_ok and all(
            (field.get(key) == value if key != 'domain_coordinate'
             else _field(field, 'domain', 'coordinate') == value)
            for key, value in expected.items())
        own_rows_ok = own_rows_ok and field.get('physical_Pstar_applied_exactly_once') is True
        own_rows_ok = own_rows_ok and field.get('pressure_is_not_n1_pressure_forcing') is True
        own_rows_ok = own_rows_ok and field.get('Q_is_auxiliary_Ur_over_sqrt_R_half') is True
    check('physical_own_chart_fields_exact', own_rows_ok)
    check('physical_own_radius_literal_rows',
          _node(own_nodes, 33804) == {'operation': 'sum', 'arguments': [31161, 33579]}
          and _node(own_nodes, 33805) == {'operation': 'analytic_unary', 'name': 'exp', 'argument': 33804}
          and _node(own_nodes, 33981) == {'operation': 'sum', 'arguments': [31161, 33830]}
          and _node(own_nodes, 33982) == {'operation': 'analytic_unary', 'name': 'exp', 'argument': 33981}
          and _node(own_nodes, 34271) == {'operation': 'sum', 'arguments': [31161, 34006]}
          and _node(own_nodes, 34272) == {'operation': 'analytic_unary', 'name': 'exp', 'argument': 34271})
    check('physical_own_log_radius_literal_rows',
          _node(own_nodes, 31161) == {'operation': 'sum', 'arguments': [31157, 31160]})
    check('physical_current_radius_literal_products',
          _node(master_nodes, 280140) == {'operation': 'analytic_unary', 'name': 'exp', 'argument': 33579}
          and _node(master_nodes, 280141) == {'operation': 'product', 'arguments': [263349, 280140]}
          and _node(master_nodes, 280154) == {'operation': 'analytic_unary', 'name': 'exp', 'argument': 33830}
          and _node(master_nodes, 280155) == {'operation': 'product', 'arguments': [263349, 280154]}
          and _node(master_nodes, 280168) == {'operation': 'analytic_unary', 'name': 'exp', 'argument': 34006}
          and _node(master_nodes, 280169) == {'operation': 'product', 'arguments': [263349, 280168]})
    geometry_substitution = _dict(physical_bindings.get('geometry_substitution'))
    check('physical_geometry_substitution_exact',
          geometry_substitution == {
              '31161': 280182, '31162': 263349, '33805': 280141,
              '33982': 280155, '34272': 280169})
    substitution = _dict(physical_component.get('substitution_receipt'))
    check('physical_saved_substitution_is_Z_independent',
          substitution.get('source_coefficients_or_field_values_evaluated') is False
          and substitution.get('ordinary_Z_preserved_by_exact_C5_inlet_and_Z_independent_geometry') is True
          and substitution.get('old_geometry_outside_current_ports_present') is False)
    check('physical_source_audit_accepts_current_frame',
          physical.get('current_native_physical_pulse_target_defined') is True
          and physical_component.get('actual_native_pulse_geometry_source_bound') is True
          and physical_component.get('actual_native_physical_pulse_C3_functions_defined') is True
          and physical_component.get('ordinary_Z_derivatives_preserved_under_geometry_binding') is True
          and physical_component.get('independent_complete_physical_pulse_m15_target_defined') is True)
    check('physical_universal_C3_template_controls_measures_endpoints',
          set(_dict(physical_bindings.get('integrals'))) >= {'main', 'gap', 'end'}
          and set(_dict(physical_bindings.get('endpoints'))) >= {'main', 'gap', 'end'}
          and set(_dict(physical_bindings.get('pieces'))) >= {'main', 'gap', 'end'}
          and set(_dict(physical_bindings.get('leading_fields'))) >= {'main', 'gap', 'end'}
          and len(_list(physical_bindings.get('global_endpoints'))) == 2
          and _field(physical_bindings, 'pressure', 'm15') == [287417, 287491, 287663, 287956]
          and len(_list(physical_bindings.get('negative_complete_m15_target'))) == 4
          and len(_list(physical_bindings.get('pulse_geometry_correction_m15'))) == 4
          and physical_prereq.get('whole_function_native_frame_qualification') is True
          and physical_prereq.get('complete_C3_source_qualification') is True
          and physical_prereq.get('current_input_ports_unchanged') is True)

    parent_interface = _dict(parent_component.get('source_interface'))
    parent_physical = _dict(parent_interface.get('physical_pulse_source'))
    parent_native = _dict(parent_interface.get('native_base_source'))
    parent_quiet = _dict(parent_interface.get('quiet_source'))
    parent_corrected = _dict(parent_interface.get('corrected_source'))
    check('physical_old_and_new_pulse_ports_saved',
          parent_physical.get('graph_digest') == PHYSICAL_GRAPH
          and parent_physical.get('pressure_m15') == [287417, 287491, 287663, 287956]
          and parent_native.get('graph_digest') == BASE_GRAPH
          and parent_native.get('old_frame_only') == [279426, 279542, 279770, 280139]
          and parent_native.get('old_pulse') == [35556, 35557])
    check('parent_terminal_source_roles_saved',
          parent_quiet.get('graph_digest') == QUIET_GRAPH
          and parent_quiet.get('native_Rc') == [258845, 258846]
          and parent_quiet.get('native_2Rc') == [262906, 262907]
          and parent_quiet.get('old_compact') == [31155, 31156]
          and parent_corrected.get('graph_digest') == CORRECTED_GRAPH
          and parent_corrected.get('corrected_active') == [53466, 53470])
    check('physical_all_three_chart_endpoints_retained',
          _field(physical_component, 'physical_scale_receipt', 'all_three_charts_and_finite_endpoints_retained') is True
          and len(_field(physical_component, 'physical_pulse_bindings', 'coordinates', default={})) == 3)
    check('physical_shared_query_domain_and_derivative_scope',
          _field(master_bindings, 'definition', 'actual_current_Rp_eight_C5_jets_are_direct_callback_inputs') is True
          and _field(master_bindings, 'definition', 'original_mu_phase_Gamma_kernels_and_nonzero_histories_retained') is True
          and _field(master_bindings, 'definition', 'numerical_interval_selected_owner_not_created_or_modified') is True)
    pulse_stage = all(checks[key] for key in (
        'physical_saved_source_prerequisites',
        'physical_current_frame_matches_master',
        'physical_geometry_matches_master',
        'physical_coordinates_match_master',
        'physical_query_y_coordinates_are_exact',
        'physical_own_chart_fields_exact',
        'physical_own_radius_literal_rows',
        'physical_own_log_radius_literal_rows',
        'physical_current_radius_literal_products',
        'physical_geometry_substitution_exact',
        'physical_saved_substitution_is_Z_independent',
        'physical_source_audit_accepts_current_frame',
        'physical_universal_C3_template_controls_measures_endpoints',
        'physical_old_and_new_pulse_ports_saved',
        'physical_all_three_chart_endpoints_retained',
        'physical_shared_query_domain_and_derivative_scope'))
    check('current_physical_complete_pulse_C3_source_identity_complete', pulse_stage)

    # The quiet source has the closed primitive -K^2*length/4.  The saved
    # parent already proves its ordinary-Z1 derivative is zero; Z0 requires
    # the newly bound K and absolute-radius source identities.
    quiet_geometry = _dict(quiet_checks.get('actual_source_geometry'))
    quiet_pressure = _dict(quiet_checks.get('actual_physical_pressure'))
    check('quiet_source_geometry_saved',
          quiet_geometry.get('passed') is True
          and quiet_geometry.get('exact_2Rc_s_Tw_minus2_minus_log2_Rp_geometry') is True
          and quiet_geometry.get('original_constant_K_source_Z_independent') is True
          and quiet_geometry.get('true_compact_exit_and_all_original_nonzero_histories_preserved') is True)
    check('quiet_source_pressure_saved',
          quiet_pressure.get('passed') is True
          and quiet_pressure.get('actual_general_pressure_equals_closed_original_quiet_primitive_C3_rows') == 8
          and quiet_pressure.get('finite_nonzero_initial_pressure_and_original_rminus_boundary_retained') is True
          and quiet_pressure.get('physical_R_and_Rweighted_signs_checked') is True)
    check('quiet_parent_length_identity_saved',
          parent_component.get('native_old_quiet_length_identity_proved') is True
          and radius.get('quiet_length_difference') == [])
    comparison_identities = _list(parent_component.get('comparison_identities'))
    saved_z1_quiet = any(
        row.get('ordinary_Z_order') == 1
        and row.get('quiet_defect_Z1_zero_proved') is True
        and row.get('quiet_derivative_removed_only_by_saved_proof_Z1', True) is not False
        for row in comparison_identities if isinstance(row, dict))
    check('quiet_parent_Z1_zero_saved', saved_z1_quiet)
    quiet_stage = all(checks[key] for key in (
        'native_original_quiet_K_source_identity_complete',
        'quiet_source_geometry_saved',
        'quiet_source_pressure_saved',
        'quiet_parent_length_identity_saved',
        'quiet_parent_Z1_zero_saved'))
    check('physical_quiet_Z0_source_identity_complete', quiet_stage)

    # Read, but do not rebuild, the parent terminal normal forms.  The new
    # graph below contracts only the two qualified identities and leaves the
    # three signed source terms visible.
    parent_graph = _dict(parent_component.get('definition_graph'))
    parent_outputs = _list(parent_graph.get('outputs'))
    parent_normal_forms = []
    parent_nodes_by_id = {
        row.get('id'): row for row in _list(parent_graph.get('nodes'))
        if isinstance(row, dict) and isinstance(row.get('id'), str)
    }
    expected_parent_atoms = {
        0: {
            'original_native_physical_terminal_log_repair_v1:s12': '1',
            'original_native_physical_terminal_log_repair_v1:s13': '-1',
            'original_native_physical_terminal_log_repair_v1:s14': '1',
            'original_native_physical_terminal_log_repair_v1:s4': '-1',
            'original_native_physical_terminal_log_repair_v1:s5': '1',
            'original_native_physical_terminal_log_repair_v1:s9': '-1',
        },
        1: {
            'original_native_physical_terminal_log_repair_v1:s49': '-1',
            'original_native_physical_terminal_log_repair_v1:s50': '1',
            'original_native_physical_terminal_log_repair_v1:s54': '-1',
            'original_native_physical_terminal_log_repair_v1:s57': '1',
            'original_native_physical_terminal_log_repair_v1:s58': '-1',
        },
    }
    parent_normal_forms_ok = (
        parent_graph.get('namespace') == 'original_native_physical_terminal_log_repair_v1'
        and parent_graph.get('node_count') == 90
        and parent_graph.get('nodes_digest') == 'bbd07075dd4f07953a30ffff7f2c0674d1ae94ff541ef29755267113d6837c5c'
        and parent_graph.get('proof_graph_not_numeric_evaluator') is True
        and len(parent_outputs) == 2)
    for order in (0, 1):
        row = next((item for item in parent_outputs
                    if isinstance(item, dict) and item.get('ordinary_Z_order') == order), None)
        row = row if isinstance(row, dict) else {}
        terminal = row.get('canonical_terminal_difference')
        geometry = row.get('canonical_geometry_decomposition')
        parent_normal_forms.append({
            'ordinary_Z_order': order,
            'canonical_terminal_difference': terminal if isinstance(terminal, dict) else {},
            'canonical_geometry_decomposition': geometry if isinstance(geometry, dict) else {},
            'exact_saved_equality': terminal == geometry,
        })
        parent_normal_forms_ok = (
            parent_normal_forms_ok and terminal == geometry
            and terminal == expected_parent_atoms[order])
    check('parent_saved_terminal_normal_forms_read_only', parent_normal_forms_ok)
    check('parent_saved_canonical_atoms_are_source_ports', all(
        isinstance(parent_nodes_by_id.get(atom), dict)
        and parent_nodes_by_id[atom].get('operation') == 'authenticated_source_function_port'
        for atom in {
            atom for mapping in expected_parent_atoms.values() for atom in mapping
        }))
    check('parent_saved_frame_only_to_old_pulse_identity', all(
        isinstance(row.get('frame_only_to_old_pulse_difference'), str)
        and _dict(parent_nodes_by_id.get(
            row.get('frame_only_to_old_pulse_difference'))).get('identity') ==
            'frame_only_pulse-old_pulse'
        for row in parent_outputs if isinstance(row, dict)))
    check('parent_saved_terminal_decomposition_accepted',
          parent_component.get('passed') is True
          and parent_component.get('signed_physical_terminal_decomposition_proved') is True
          and parent_component.get('absolute_radius_log_recipe_difference_proved') is True)

    # Build a fresh source-port certificate graph.  Its atoms are handles from
    # the authenticated saved graphs; no saved graph is rebuilt or reduced.
    certificate = _SourceCertificateGraph()
    terminal_outputs = []
    parent_contractions = []
    parent_atom_role_maps = []
    physical_ports = [287417, 287491]
    old_frame_only_ports = [279426, 279542]
    old_pulse_ports = [35556, 35557]
    native_2rc_ports = [262906, 262907]
    native_rc_ports = [258845, 258846]
    corrected_ports = [53466, 53470]
    old_compact_ports = [31155, 31156]
    quiet_defect_ports = ['s19', 's41']
    for order in (0, 1):
        current_pulse = certificate.port(
            'current_physical_pulse_C3', physical_ports[order], order, PHYSICAL_GRAPH)
        old_frame_only = certificate.port(
            'saved_old_frame_only_pulse', old_frame_only_ports[order], order, BASE_GRAPH)
        old_pulse = certificate.port(
            'saved_old_pulse', old_pulse_ports[order], order, BASE_GRAPH)
        pulse_to_frame = certificate.identity(
            current_pulse, old_frame_only, pulse_stage,
            'current_complete_physical_pulse_equals_old_frame_only_pulse',
            ordinary_Z_orders=[0, 1, 2, 3],
            physical_geometry_delta_C3='0',
            pressure_ports=[287417, 287491, 287663, 287956])
        frame_to_old = certificate.identity(
            old_frame_only, old_pulse,
            checks.get('parent_saved_frame_only_to_old_pulse_identity') is True,
            'saved_frame_only_pulse_equals_saved_old_pulse',
            parent_graph_namespace=parent_graph.get('namespace'),
            parent_output='frame_only_to_old_pulse_difference',
            source_nodes={
                'old_frame_only': old_frame_only_ports[order],
                'old_pulse': old_pulse_ports[order],
            })
        pulse_identity = certificate.compose_identity(
            current_pulse, old_frame_only, old_pulse,
            pulse_stage and checks.get('parent_saved_frame_only_to_old_pulse_identity') is True,
            'current_complete_physical_pulse_equals_saved_old_pulse_by_frame_only_bridge',
            first_identity='current_complete_physical_pulse_equals_old_frame_only_pulse',
            second_identity='saved_frame_only_pulse_equals_saved_old_pulse',
            source_map={
                'current_physical': physical_ports[order],
                'old_frame_only': old_frame_only_ports[order],
                'old_pulse': old_pulse_ports[order],
            })
        pulse_difference = certificate.signed(
            [(1, pulse_identity), (-1, old_pulse)],
            ordinary_Z_orders=[0, 1, 2, 3],
            source_identity='pulseC3congruence',
            saved_frame_only_bridge=frame_to_old)
        quiet_defect = certificate.port(
            'saved_quiet_defect', quiet_defect_ports[order], order, QUIET_DEFECT_GRAPH)
        quiet_cancellation = certificate.cancel(
            quiet_defect,
            quiet_stage if order == 0 else checks.get('quiet_parent_Z1_zero_saved') is True,
            ('native_quiet_Z0_primitive_minus_old_quiet_primitive'
             if order == 0 else 'saved_quiet_Z1_zero'),
            primitive='-K^2*length/4', quiet_length_difference=[],
            ordinary_Z_order=order)
        native_2rc = certificate.port(
            'native_quiet_2Rc', native_2rc_ports[order], order, QUIET_GRAPH)
        native_rc = certificate.port(
            'native_quiet_Rc', native_rc_ports[order], order, QUIET_GRAPH)
        corrected_active = certificate.port(
            'corrected_active_physical', corrected_ports[order], order, CORRECTED_GRAPH)
        old_compact = certificate.port(
            'old_compact_physical', old_compact_ports[order], order, QUIET_GRAPH)
        source_refs_by_node = {
            physical_ports[order]: current_pulse,
            old_frame_only_ports[order]: old_frame_only,
            old_pulse_ports[order]: old_pulse,
            native_2rc_ports[order]: native_2rc,
            native_rc_ports[order]: native_rc,
            corrected_ports[order]: corrected_active,
            old_compact_ports[order]: old_compact,
        }
        parent_atom_refs = []
        parent_substitutions = []
        parent_source_map = {}
        parent_atom_role_map = {}
        parent_atom_map_valid = True
        for parent_atom, coefficient in expected_parent_atoms[order].items():
            parent_node = _dict(parent_nodes_by_id.get(parent_atom))
            source_node = parent_node.get('source_node')
            if parent_node.get('role') == 'accepted_quiet_output':
                mapped_ref = quiet_defect
                source_node_valid = str(source_node).endswith(
                    quiet_defect_ports[order])
            else:
                mapped_ref = source_refs_by_node.get(source_node)
                source_node_valid = mapped_ref is not None
            expected_digest = (
                QUIET_DEFECT_GRAPH if parent_node.get('role') == 'accepted_quiet_output'
                else PHYSICAL_GRAPH if source_node == physical_ports[order]
                else CORRECTED_GRAPH if source_node == corrected_ports[order]
                else BASE_GRAPH if source_node in (
                    old_frame_only_ports[order], old_pulse_ports[order])
                else QUIET_GRAPH if source_node in (
                    native_2rc_ports[order], native_rc_ports[order], old_compact_ports[order])
                else None)
            parent_atom_map_valid = (
                parent_atom_map_valid
                and parent_node.get('operation') == 'authenticated_source_function_port'
                and source_node_valid
                and parent_node.get('ordinary_Z_order') == order
                and expected_digest is not None
                and parent_node.get('source_graph_digest') == expected_digest)
            if mapped_ref is None:
                mapped_ref = certificate.append(
                    'unresolved_saved_parent_source_port', [],
                    parent_atom=parent_atom, ordinary_Z_order=order)
            parent_ref = certificate.parent_atom(
                parent_atom, mapped_ref, coefficient,
                ordinary_Z_order=order,
                source_node=source_node,
                source_role=parent_node.get('role'),
                source_graph_digest=parent_node.get('source_graph_digest'))
            parent_atom_refs.append(parent_ref)
            parent_source_map[parent_atom] = {
                'source_node': source_node,
                'source_role': parent_node.get('role'),
                'mapped_certificate_ref': mapped_ref,
                'coefficient': coefficient,
            }
            parent_atom_role_map[parent_atom] = parent_node.get('accepted_field')
            if mapped_ref == current_pulse:
                parent_substitutions.append((coefficient, pulse_identity))
            elif mapped_ref == quiet_defect:
                parent_substitutions.append((coefficient, quiet_cancellation))
            else:
                parent_substitutions.append((coefficient, mapped_ref))
        parent_contracted = certificate.contract_parent(
            parent_atom_refs, parent_substitutions,
            parent_atom_map_valid and pulse_stage and quiet_stage,
            ordinary_Z_order=order,
            parent_graph_namespace=parent_graph.get('namespace'),
            parent_canonical_atoms=expected_parent_atoms[order],
            parent_source_map=parent_source_map,
            parent_source_fields=parent_atom_role_map,
            qualified_aliases={
                'physical_to_old_pulse': pulse_identity,
                'quiet_defect_to_zero': quiet_cancellation,
            },
            source_identity='saved_parent_canonical_terminal_contracted_once')
        parent_contractions.append(parent_contracted)
        parent_atom_role_maps.append(parent_source_map)
        remaining = parent_contracted
        active = certificate.signed(
            [(1, native_rc), (-1, corrected_active)],
            ordinary_Z_order=order,
            source_identity='active_nativeRc_minus_correctedactive',
            zero_claimed=False)
        native_compact = certificate.signed(
            [(1, native_2rc), (-1, native_rc)],
            ordinary_Z_order=order,
            source_identity='native2Rc_minus_nativeRc',
            zero_claimed=False)
        compact = certificate.signed(
            [(1, native_compact), (-1, old_compact)],
            ordinary_Z_order=order,
            source_identity='native_compact_minus_oldcompact',
            zero_claimed=False)
        active_compact = certificate.signed(
            [(1, active), (1, compact)],
            ordinary_Z_order=order,
            source_identity='active_plus_compact_decomposition',
            zero_claimed=False)
        qualified_terminal = certificate.signed(
            [(1, remaining), (1, quiet_cancellation), (1, pulse_difference)],
            ordinary_Z_order=order,
            source_identity='qualified_saved_terminal_difference',
            zero_claimed=False)
        terminal_outputs.append({
            'ordinary_Z_order': order,
            'canonical_terminal_difference': {
                'source_terms': [
                    {'source': 'native_2Rc', 'coefficient': '1', 'source_node': native_2rc_ports[order]},
                    {'source': 'corrected_active', 'coefficient': '-1', 'source_node': corrected_ports[order]},
                    {'source': 'old_compact', 'coefficient': '-1', 'source_node': old_compact_ports[order]},
                ],
                'normal_form': certificate.normal(remaining),
                'zero_claimed': False,
                'graph_ref': remaining,
            },
            'canonical_geometry_decomposition': {
                'active': certificate.normal(active),
                'compact': certificate.normal(compact),
                'normal_form': certificate.normal(active_compact),
                'zero_claimed': False,
                'graph_ref': active_compact,
            },
            'qualified_terminal_difference': {
                'normal_form': certificate.normal(qualified_terminal),
                'graph_ref': qualified_terminal,
                'all_qualified_identities_admitted': bool(pulse_stage and quiet_stage),
            },
            'saved_parent_canonical_contract': {
                'graph_ref': parent_contracted,
                'ordinary_Z_order': order,
                'parent_atoms': expected_parent_atoms[order],
                'source_map': parent_atom_role_maps[order],
                'normal_form': certificate.normal(parent_contracted),
            },
            'physical_pulse_to_frame_only_graph_ref': pulse_to_frame,
            'saved_frame_only_to_old_pulse_graph_ref': frame_to_old,
            'physical_pulse_identity_graph_ref': pulse_identity,
            'quiet_Z0_cancellation_graph_ref': quiet_cancellation,
            'active_difference_graph_ref': active,
            'compact_difference_graph_ref': compact,
            'exact_saved_parent_normal_form': parent_normal_forms[order],
        })

    check('fresh_parent_canonical_source_map', all(
        len(mapping) == len(expected_parent_atoms[order])
        and all(isinstance(item, dict) and item.get('mapped_certificate_ref')
                for item in mapping.values())
        for order, mapping in enumerate(parent_atom_role_maps)))
    check('fresh_parent_canonical_terminal_contracted_once', all(
        _dict(next((node for node in certificate.nodes
                    if node.get('id') == ref), {})).get('operation') ==
        'contract_saved_parent_canonical_terminal'
        and _dict(next((node for node in certificate.nodes
                        if node.get('id') == ref), {})).get('admitted') is True
        for ref in parent_contractions))
    check('contracted_parent_equals_active_compact_signed_form', all(
        row['canonical_terminal_difference']['normal_form'] ==
        row['canonical_geometry_decomposition']['normal_form']
        and len(row['canonical_terminal_difference']['normal_form']) == 3
        for row in terminal_outputs))

    terminal_stage = (
        parent_normal_forms_ok
        and checks.get('parent_saved_terminal_decomposition_accepted') is True
        and pulse_stage and quiet_stage
        and all(row['canonical_terminal_difference']['zero_claimed'] is False
                for row in terminal_outputs)
        and all(row['canonical_geometry_decomposition']['zero_claimed'] is False
                for row in terminal_outputs))
    check('fresh_source_certificate_graph_complete',
          certificate.nodes != [] and len(terminal_outputs) == 2)
    check('fresh_source_certificate_graph_preserves_three_signed_terms',
          all(row['canonical_terminal_difference']['source_terms'] == [
              {'source': 'native_2Rc', 'coefficient': '1', 'source_node': native_2rc_ports[row['ordinary_Z_order']]},
              {'source': 'corrected_active', 'coefficient': '-1', 'source_node': corrected_ports[row['ordinary_Z_order']]},
              {'source': 'old_compact', 'coefficient': '-1', 'source_node': old_compact_ports[row['ordinary_Z_order']]},
          ] for row in terminal_outputs))
    terminal_stage = terminal_stage and checks.get(
        'fresh_source_certificate_graph_complete') is True
    terminal_stage = terminal_stage and checks.get(
        'fresh_source_certificate_graph_preserves_three_signed_terms') is True
    terminal_stage = terminal_stage and all(checks.values())
    check('terminal_source_reduction_premises_complete', terminal_stage)

    stage_values = {
        'native_logC_original_control_identity_proved': logc_stage,
        'native_original_absolute_Rc_Rp_identity_proved': radius_stage,
        'old_preRp_Pstar_source_identity_proved': pstar_stage,
        'native_original_quiet_K_identity_proved': pstar_stage and kappa_stage,
        'physical_quiet_Z0_source_identity_proved': quiet_stage,
        'current_physical_complete_pulse_C3_identity_proved': pulse_stage,
        'terminal_difference_reduced_to_active_compact': terminal_stage,
    }
    unresolved = sorted(key for key, value in checks.items() if value is not True)
    source_interface = {
        'required_alias_keys': ['logC', 'Rc', 'Rp', 'Pstar', 'quietK'],
        'required_identity_keys': [
            'pulseC3congruence', 'quietZ0cancellation',
            'signedremainingnative2Rc-correctedactive-oldcompact'],
        'source_aliases': SOURCE_ALIASES,
        'source_identities': SOURCE_IDENTITIES,
        'remaining_physical_C1_ports': {
            'main': {'query': 33578, 'y': 33579, 'log_radius': 33804, 'radius': 33805},
            'gap': {'query': 33829, 'y': 33830, 'log_radius': 33981, 'radius': 33982},
            'end': {'query': 34005, 'y': 34006, 'log_radius': 34271, 'radius': 34272},
            'pressure_m15': [287417, 287491, 287663, 287956],
            'geometry_delta_C3': [290917, 290919, 290921, 290923],
        },
        'physical_pulse_source_ports': {
            'current_pressure_m15': [287417, 287491, 287663, 287956],
            'old_frame_only': [279426, 279542, 279770, 280139],
            'old_pulse': [35556, 35557],
            'negative_complete_m15_target': _list(
                physical_bindings.get('negative_complete_m15_target')),
            'pulse_geometry_correction_m15': _list(
                physical_bindings.get('pulse_geometry_correction_m15')),
            'ordinary_Z_orders': [0, 1, 2, 3],
            'geometry_delta_C3': '0',
        },
        'saved_parent_canonical_source_map': {
            str(order): parent_atom_role_maps[order] for order in (0, 1)
        },
    }
    definition_graph = {
        'namespace': SCHEMA + '_certificate_graph',
        'nodes': certificate.nodes,
        'node_count': len(certificate.nodes),
        'nodes_digest': _digest(certificate.nodes),
        'outputs': terminal_outputs,
        'proof_graph_not_numeric_evaluator': True,
        'source_graph_digests': {
            'native_base': BASE_GRAPH,
            'physical_pulse': PHYSICAL_GRAPH,
            'quiet': QUIET_GRAPH,
            'corrected': CORRECTED_GRAPH,
            'quiet_defect': QUIET_DEFECT_GRAPH,
        },
    }
    result = {
        'schema': SCHEMA,
        'passed': all(stage_values.values()) and not unresolved,
        'source_only': True,
        'ordinary_Z_orders': [0, 1],
        'ordinary_Z_not_temporal_n': True,
        'old_error_terms_removed': False,
        'new_unentered_expression_simulations': 0,
        'legacy_numeric_chain_reconfigured': False,
        'pulse_C3_ordinary_Z_orders': [0, 1, 2, 3],
        'source_family': family,
        'retained_construction_data': retained,
        'retained_construction_counts': RETAINED_COUNTS,
        'source_prerequisite_audit': checks,
        'definition_graph': definition_graph,
        'source_interface': source_interface,
        'source_inspection_audit': {
            'source_only': True,
            'accepted_ancestor_operator_calls': 0,
            'old_operators_replayed': False,
            'numeric_values_computed': False,
            'new_source_consumer_calls': 1,
            'ordinary_Z_not_temporal_n': True,
            'old_error_terms_removed': False,
            'new_unentered_expression_simulations': 0,
            'legacy_numeric_chain_reconfigured': False,
        },
        'source_witnesses': {
            'pstar_provider': _dict(master_pstar.get('source_rows_provider')),
            'kappa_frame': {
                'frame_record': FRAME,
                'frame_source_record': FRAME_SOURCE,
                'frame_check_record': FRAME_CHECK,
                'pre_record': PRE,
                'pre_source_record': PRE_SOURCE,
                'primitive_record': PRIMITIVE,
                'primitive_check_record': PRIMITIVE_CHECK,
                'primitive_source_record': PRIMITIVE_SOURCE,
                'frame_history_m_nodes': frame_m_rows,
                'pre_kappa_node': pre_quiet.get('kappa'),
                'pre_K_node': pre_quiet.get('K'),
                'pre_mapping_size': len(pre_mapping),
                'linear_mass_source_expression': primitive_identity.get('mean_common'),
                'quotient_kappa_source_expression': primitive_identity.get('kappa'),
                'quotient_rule_source_ast': "self.kappa = s.cancel(self.mean_common / self.Z)",
                'ordinary_derivative_source_expression': primitive_identity.get(
                    'kappa_Z_derivative'),
                'analytic_Z0_extension_admitted': checks.get(
                    'quiet_kappa_analytic_extension_from_linear_source') is True,
                'authenticated_raw_hashes': {
                    name: RECORD_HASHES[name] for name in (
                        FRAME, FRAME_SOURCE, FRAME_CHECK, PRE, PRE_SOURCE,
                        PRIMITIVE, PRIMITIVE_CHECK, PRIMITIVE_SOURCE)
                },
                'decoded_hashes': {
                    name: DECODED_RECORD_HASHES[name] for name in (FRAME, PRE)
                },
            },
        },
        'SOURCE_ALIASES': SOURCE_ALIASES,
        'SOURCE_IDENTITIES': SOURCE_IDENTITIES,
        'source_aliases': SOURCE_ALIASES,
        'source_identities': SOURCE_IDENTITIES,
        'parent_saved_terminal_normal_forms': parent_normal_forms,
        'native_original_quiet_K_identity_proved': pstar_stage and kappa_stage,
        'physical_quiet_Z0_source_identity_proved': quiet_stage,
        'old_physical_geometry_alias_proved': radius_stage,
        'old_physical_geometry_alias_scope': 'absolute Rc/Rp source recipes only',
        'new_physical_pulse_equals_old_pulse_proved': pulse_stage,
        'terminal_difference_reduced_to_active_compact': terminal_stage,
        'unresolved_gates': unresolved,
        'runtime_model_metadata': None,
        'source_audit': {
            'new_source_consumer_calls': 1,
            'old_math_consumers_called': 0,
            'old_definition_graphs_rebuilt': False,
            'accepted_mathematical_operators_replayed': False,
            'source_values_or_bounds_evaluated': False,
            'numerical_values_computed': False,
            'legacy_numeric_chain_unchanged': True,
        },
        'input_record_hashes': dict(RECORD_HASHES),
        'decoded_record_hashes': dict(DECODED_RECORD_HASHES),
        'static_receipt_metadata': dict(STATIC_RECEIPT_HASHES),
    }
    for flag, value in stage_values.items():
        result[flag] = bool(value)
    for flag in FALSE_FLAGS:
        result[flag] = False
    return result

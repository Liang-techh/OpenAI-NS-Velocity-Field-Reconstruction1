"""Bind one new lower-collar radial FTC expression to the archived source."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
COMPONENT = 'original_collar_signed_ftc_v2'
PREFIX = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_'
ARCHIVE = PREFIX + 'n1_lower_collar_functions.json.gz'
CHECK = PREFIX + 'n1_lower_collar_functions_check.json'
CORE = PREFIX + 'core_collar_pressure_identity.json'
PINS = {
    ARCHIVE: 'b2f983ab934875ae9aeb24b292277911dd5647e3a23943ddcf2b3b6e82595cb5',
    CHECK: '1a7b6b370b5ecf7d6b7e97001e0311e9525698ecf3389fc3da608af308a13ba9',
    CORE: 'ae91b0d72d0415e8146f7ee9926d1a35ddbedffff6b89eac1f26558dbb5b44c7',
}
LEAF = 'src/openai_ns_reconstruction/original_collar_signed_ftc.py'
GRAPH = 'src/openai_ns_reconstruction/original_signed_density_expression.py'
RUNNER = 'scripts/run_original_collar_signed_ftc_v2_once.py'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_new(path, result):
    with path.open('x', encoding='utf-8', newline='\n') as out:
        json.dump(result, out, sort_keys=True, separators=(',', ':'), allow_nan=False)
        out.write('\n')


def consume(inputs, result):
    from openai_ns_reconstruction.original_collar_signed_ftc import build
    from openai_ns_reconstruction.original_signed_density_expression import SignedDensityExpressionGraph

    lower = json.loads(gzip.decompress(inputs[ARCHIVE]))
    checked = json.loads(inputs[CHECK])
    core_report = json.loads(inputs[CORE])
    core = core_report['actual_core_pressure_m15_primitive']
    core_port = core_report['six_source_identity_progress']['canonical_core_m15_port']
    require(checked['all_passed'] and checked['source_family'] == lower['source_family']
            and checked['input_hashes'][Path(ARCHIVE).name] == PINS[ARCHIVE]
            and checked['checks']['actual_source']['checked_true_core_nodes_only']
            and checked['checks']['ratio']['rho4_ratio_inlet_exact_one'], 'actual core-to-collar receipt mismatch')
    require(core['actual_symbolic_primitive_identity'] and core_port == [53125, 53126],
            'core pressure occurrence differs')

    saved = lower['lower_collar_bindings']
    fields, jets, parameters = (saved[x] for x in ('fields', 'core_radial_rows', 'parameters'))
    nodes = lower['exact_graph_nodes']
    s_id, r_id = fields['query'], fields['rho']
    require((s_id, r_id, parameters['hb'], parameters['epsilon']) == (35848, 35851, 8, 35565),
            'radial source parameters changed')
    require(fields['ratio'][:2] == [35940, 35941] and nodes[35940]['argument'] == fields['logratio'][0]
            and nodes[fields['logratio'][0]]['original_ratio_ODE'] is True,
            'ratio is not the genuine log-Volterra source')
    require(nodes[35949]['argument'] == nodes[35953]['argument'] == s_id and nodes[35953]['name'] == 'original_flat_sigma_prime',
            'sigma and radial derivative do not share the historical source')
    require(jets['Psi'][1][:2] == [35715, 35833], 'kernel radial derivative history differs')
    require(all(nodes[x]['rho'] == 35558 and nodes[x]['Z_variable'] == 20
                and nodes[x]['quantity'] == name and nodes[x]['radial_order'] == order
                and nodes[x]['ordinary_slow_Z_derivative_order'] == k
                for name, order in (('Phi', 0), ('Phi', 1), ('Psi', 1), ('Psi', 2))
                for k, x in enumerate(jets[name][order][:2])), 'core r/Z derivative order mismatch')

    g = SignedDensityExpressionGraph()
    def ref(index):
        return g.reference(Path(ARCHIVE).name, PINS[ARCHIVE], index, source_nodes=nodes)

    def as_pair(ids):
        return tuple(ref(x) for x in ids[:2])

    def fixed(index):
        return ref(index), g.zero

    def at_rho(ids, point):
        return tuple(g.insert({'operation': 'function_substitution', 'expression': ref(index).node,
                               'variable': ref(35558).node, 'value': point.node}) for index in ids[:2])

    f, f_z = at_rho(jets['Phi'][0], ref(r_id))
    p, p_z = at_rho(jets['Phi'][1], ref(r_id))
    ports = {
        'Z': (ref(20), g.one), 'L': as_pair([30382, 30381]),
        'hb': fixed(parameters['hb']), 'epsilon': fixed(parameters['epsilon']), 'rho': fixed(r_id),
        'ratio': as_pair(fields['ratio']), 'chi': fixed(fields['chi']), 'sigma': fixed(35949), 'sigma_prime': fixed(35953),
        'phi_r_over_phi': (p/f, p_z/f - p*f_z/f**2),
        'psi_r': at_rho(jets['Psi'][1], ref(r_id)), 'psi_rr': at_rho(jets['Psi'][2], ref(r_id)),
        'psi_r_at_4': at_rho(jets['Psi'][1], ref(13)),
    }
    result['actual_input_handles'] = {'query': s_id, 'rho': r_id, 'phi_rho_orders_0_1': [x[:2] for x in jets['Phi'][:2]],
                                      'psi_rho_orders_1_2': [x[:2] for x in jets['Psi'][1:3]], 'ratio': fields['ratio'][:2]}
    history = {'report': CORE, 'sha256': PINS[CORE], 'canonical_core_C1_port': core_port}
    result['actual_FTC'] = build(source_ports=ports, limits={'lower': 0, 'upper': parameters['sc_half']},
                                 incoming_history=history, graph=g, source_coordinate=ref(s_id))
    require(result['actual_FTC']['passed'], str(result['actual_FTC'].get('gaps')))
    result['forward212_same_source_proved'] = False
    result['legacy_errors_removed'] = False
    result['m15_value_evaluated'] = False


def main():
    path = Path(sys.argv[1]).resolve()
    manifest_raw = path.read_bytes()
    manifest = json.loads(manifest_raw)
    files = manifest['files']
    require(manifest['component'] == COMPONENT and manifest['one_entry_only'] is True
            and set(files) == {*PINS, LEAF, GRAPH, RUNNER}, 'wrong frozen inputs')
    require(all(files[f] == sha for f, sha in PINS.items()), 'historical pins changed')
    raw = {}
    for name, sha in files.items():
        raw[name] = (ROOT / name).read_bytes()
        require(digest(raw[name]) == sha, f'changed raw input {name}')
    brief_raw = Path(manifest['brief_path']).read_bytes()
    brief = json.loads(brief_raw)
    require(brief['phase'] == COMPONENT + '_pending' and brief['source_checkpoint']['sha256'] == digest(manifest_raw),
            'missing pending source checkpoint')
    entry = {'component': COMPONENT, 'manifest_sha256': digest(manifest_raw), 'brief_sha256': digest(brief_raw)}
    write_new(path.parent / 'entry.json', entry)
    sys.path.insert(0, str(ROOT / 'src'))
    result = {**entry, 'raw_input_sha256': files, 'historical_producer_executed': False}
    try:
        consume(raw, result)
        result['passed'] = True
    except Exception as error:
        result.update(passed=False, error=str(error), traceback=traceback.format_exc())
    write_new(path.parent / 'result.json', result)
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())

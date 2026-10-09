"""Live whole-Z microscopic and macro source/history functions through R100.

Fresh source cells feed the original MicroFunctions and ActualMacroMoments.
The accepted analytic Dbar theorem supplies a logarithmic lower bound when
the raw interval extension crosses zero. It never selects a source value.
The micro endpoint fields and all six histories are the macro incoming.
"""
import ast
import gzip
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source as source
import lei_ren_part1_paper_compliant_current_original_micro_functions as micro
import lei_ren_part1_paper_compliant_current_original_first_switch_functions as first
import lei_ren_part1_paper_compliant_current_original_bridge_macro_moments as moments
import lei_ren_part1_paper_compliant_current_original_R100_endpoint as endpoint

macro = micro.downstream
switch = macro.downstream
fields, parameters, ep = micro.fields, micro.parameters, source.ep
HERE, PREFIX, sha, read, bind = source.HERE, source.PREFIX, source.sha, source.read, source.bind
NAME = PREFIX+'current_original_whole_Z_micro_macro.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_micro_macro_check.json'
GATE = 'current_original_whole_Z_live_micro_macro_six_history_functions_through_R100'
OPEN = source.OPEN
MICRO_CELLS = (('first_micro', '0', '.5'), ('first_micro', '.5', '1'),
               ('second_micro', '1', '1.5'), ('second_micro', '1.5', '2'))
MACRO_CELLS = (((0,1),(1,2)), ((1,2),(1,1)))


def replay_positive_callback(module, name, callback, expected_role):
    """Original arithmetic unchanged; attach a same-source analytic proof."""
    tree = ast.parse(Path(module.__file__).read_text(encoding='utf8'))
    fn = next(node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)
              and node.name == name and any(isinstance(n, ast.Call)
              and ast.unparse(n.func) == 'parameters.positive_source' for n in ast.walk(node)))
    calls = [node for node in ast.walk(fn) if isinstance(node, ast.Call)
             and ast.unparse(node.func) == 'parameters.positive_source']
    if len(calls) != 1 or len(calls[0].args) != 3 or not isinstance(calls[0].args[2], ast.Constant) \
            or calls[0].args[2].value != expected_role:
        raise ValueError('Original Dbar proof callback changed: '+name)
    env = dict(vars(module))
    env['parameters'] = SimpleNamespace(positive_source=callback)
    exec(compile(ast.Module(body=[fn], type_ignores=[]), module.__file__, 'exec'), env)
    return env[name], dict(original_method_AST_sha256=hashlib.sha256(
        ast.dump(fn).encode('utf8')).hexdigest(), unchanged_original_arithmetic=True,
        only_positive_source_proof_callback_supplied=True, exact_role=expected_role)


class WholeZMicroMacro:
    def __init__(self, dps=500):
        self.source = source.OriginalWholeZBridgeSource(dps)
        self.c = c = self.source.c
        self.hashes = dict(self.source.hashes)
        self.identity = dict(actual_five_defect_family_sha256=self.source.family,
            implicit_source_sha256=self.source.source, datum_enclosure_sha256=self.source.datum)
        for stem, gate in (
            ('current_original_whole_Z_bridge_source', source.GATE),
            ('current_inner_relaxed_inputs', 'current_original_inner_Ra_R110_relaxed_generic_input_certified')):
            name = PREFIX+stem+'_check.json'
            receipt = json.loads((HERE/name).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(gate):
                raise ValueError('Accepted whole-Z same-source prerequisite required: '+stem)
            identity = receipt.get('source_family') or {key: receipt[key] for key in self.identity}
            if identity != self.identity:
                raise ValueError('Whole-Z source/pressure/leading family differs')
            for path, digest in receipt['input_hashes'].items():
                bind(self.hashes, path, digest)
            bind(self.hashes, name, sha(name))
        self.records = {}
        for stem in ('current_inner_relaxed_inputs', 'current_inner_exit_strict_collar',
                     'global_exit_certificate', 'K1_ledger'):
            name = PREFIX+stem+'.json'
            row = json.loads((HERE/name).read_bytes())
            for path, digest in row.get('input_hashes', {}).items():
                bind(self.hashes, path, digest)
            bind(self.hashes, name, sha(name))
            self.records[stem] = row
        norm = self.source.records['physical_norm_family']
        ledger = self.records['K1_ledger']
        exit_ = self.records['global_exit_certificate']
        inner = self.records['current_inner_relaxed_inputs']
        collar = self.records['current_inner_exit_strict_collar']
        if inner['source_family'] != self.identity:
            raise ValueError('Same original whole inner input required')
        for key in ('base_analytic_core_family_sha256', 'uniform_Cstar_family_sha256'):
            if norm[key] != ledger[key] or norm[key] != exit_[key]:
                raise ValueError('Original positive-Dbar analytic family differs')
        if (not exit_['actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified']
                or not ledger['h_b_equals_epsilon_b_by_definition']
                or not all(ledger['smallness_checks'].values())
                or exit_['real_axial_domain'] != ['-1','1']
                or exit_['physical_radial_domain'] != ['Ra','110']):
            raise ValueError('Original complete analytic Dbar domain/smallness required')
        if exit_['proof']['comparison'] != '9.14 inputs, K1 ledger and shared width giveDbar>=1/(2K),|q|<=2K,Hbar>=2+2gamma,Hbar<=K^10; Dbar>=3.5 on100..110':
            raise ValueError('Exact original positive-Dbar theorem changed')
        attachment = collar['cross_graph_amplitude_packet_pressure_width_attachment']
        for key in ('base_analytic_core_family_sha256', 'uniform_Cstar_family_sha256',
                    'admitted_inner_parameter_family_sha256'):
            if attachment[key] != ledger[key] or attachment[key] != exit_[key]:
                raise ValueError('Actual source theorem attachment differs')
        if collar['exact_current_exit_source_attachment']['admitted_global_exit_prescription'] != exit_['field_prescription']:
            raise ValueError('Same original field equations required')
        self.logK = read(c, norm['selected_log_K_upper'])
        self.logh = read(c, ledger['shared_positive_width_log_enclosure'])
        if self.logh._mpi_ != self.source.logh._mpi_ or ep(self.logK-c.ln(1000000))[0] <= 0:
            raise ValueError('Same source hb and full Kmax bound required')
        self.logDlower = c.mpf(ep(-c.ln(2)-self.logK)[0])
        self.positive_theorem = dict(original_domain='Z[-1,1], Ra<R<=110',
            original_bound='Dbar>=1/(2Kmax)', log_Kmax_upper=self.logK,
            log_Dbar_positive_lower=self.logDlower, source_width_log=self.logh,
            source_identity=self.identity, inner_receipt_sha256=sha(PREFIX+'current_inner_relaxed_inputs_check.json'),
            global_exit_sha256=sha(PREFIX+'global_exit_certificate.json'),
            K1_ledger_sha256=sha(PREFIX+'K1_ledger.json'),
            physical_norm_family_sha256=sha(PREFIX+'physical_norm_family.json'),
            lower_exponential_not_materialized=True, inherited_from_analytic_leading_source_not_N1024=True)
        self.micro_evaluate, micro_binding = replay_positive_callback(micro, 'evaluate',
            self.positive_Dbar, 'actual_micro_comparison_Dbar')
        self.macro_evaluate, macro_binding = replay_positive_callback(macro, 'macro_cell',
            self.positive_Dbar, 'actual_frozen_macro_Dbar')
        self.bindings = dict(micro=micro_binding, macro=macro_binding,
            exact_source_controls=micro.source_bindings(), endpoint=endpoint.source_bindings(),
            original_micro_exit_is_macro_incoming_by_source_definition=True,
            frame_packet_or_N1024_inlet_not_used=True)
        self.owners = {}
        for module in (micro, first, moments, endpoint, macro, switch, parameters):
            name = Path(module.__file__).name
            bind(self.hashes, name, sha(name))
        bind(self.hashes, Path(__file__).name, sha(Path(__file__).name))

    def positive_Dbar(self, flow, row, role):
        """Bind the theorem to this exact original Dbar row, not a cap."""
        if role not in ('actual_micro_comparison_Dbar', 'actual_frozen_macro_Dbar'):
            raise ValueError('Only the original comparison Dbar theorem is admitted')
        parameters.same_source(flow, [row])
        if flow.c is not self.c or flow.logs[0]._mpi_ != self.logh._mpi_:
            raise ValueError('Original same-current source width/context required')
        upper = ep(row.scale.evaluate()+self.c.ln(self.c.mpf(ep(row.coefficient)[1])))[1] \
            if ep(row.coefficient)[1] > 0 else -mp.inf
        if not mp.isfinite(upper) or ep(self.logDlower)[0] > upper:
            raise ValueError('Original Dbar theorem contradicts the raw source cover')
        positive = row.positive_intersection(self.logDlower)
        proof = parameters.positive_source(flow, positive, role)
        return {**proof, 'original_raw_source_row':row,
            'theorem_backed_positive_intersection':True,
            'original_analytic_theorem':self.positive_theorem,
            'source_log_lower':self.logDlower}

    def owner(self, ends):
        Z = self.source.cell(ends)
        if Z._mpi_ in self.owners:
            return self.owners[Z._mpi_]
        flow, source_proof = self.source.owner(ends)
        c = self.c
        prepared = self.source.actual.prepare(Z)
        # Scalar series needs only flow and context, not an R100 constructor.
        series = first.FirstSwitchFunctions.__new__(first.FirstSwitchFunctions)
        series.flow, series.c = flow, c
        decoder = micro.HydratedComparison.__new__(micro.HydratedComparison)
        decoder.flow, decoder.ctx, decoder.weight = flow, c, flow.weight
        decoder.cap, decoder.Z, decoder.delta = self.source.cap, Z, self.source.core.delta
        decoder.data, decoder.inputs = prepared['data'], self.source.inputs(Z)
        decoder.comparison, decoder.cache = decoder, {}
        decoder.proof = dict(live_original_interval_inlet=True, source_frame_values_not_used=True,
            source_L3_norm_upper=decoder.data['L3_norm'], source_V3_norm_upper=decoder.data['V3_norm'])
        original_micro = micro.MicroFunctions(flow, series, decoder,
            source.bridge._encode(source_proof['fresh_actual_core']))
        # Use the actual current micro source at the seam; retain its signed
        # width terms directly instead of subtracting broad V and V0 boxes.
        at = original_micro.prefix('second_micro', c.mpf(2))
        inlet, weights = original_micro.histories('second_micro', c.mpf(2))
        flow.ell_in = at['ell']
        flow.deltaV_in = flow.add(*at['deltaV'].values())
        actual_macro = moments.ActualMacroMoments(flow, inlet)
        inputs = decoder.inputs
        p0 = inputs['p0'].truncate(5)
        ratios = [source.IntervalTaylor(c, [inputs[key][n]/math.factorial(n)
                  for n in range(6)]) for key in ('F0_ratios','F0_squared_ratios')]
        r100 = actual_macro.evaluate((1,1))
        axis = endpoint.recover_endpoint(flow, Z, decoder.delta, p0, *ratios,
            r100['same_original_macro_field_functions'], r100['actual_six_moment_functions'])
        reference = SimpleNamespace(flow=flow, c=c, z=source.IntervalTaylor.variable(c,Z,5),
                                    delta=decoder.delta, P0=axis['pressure_axis_over_Pstar_squared_axial5'])
        reference.zrows = flow.jet(reference.z)
        result = SimpleNamespace(flow=flow, c=c, Z=Z, series=series,
            micro=original_micro, macro=actual_macro, axis=axis, reference=reference,
            source_proof=source_proof, micro_inlet_weights=weights, actual_R100=r100)
        self.owners[Z._mpi_] = result
        return result

    def query(self, ends, chart, left, right):
        op = self.owner(ends)
        if chart in micro.CHARTS:
            value = self.micro_evaluate(op.micro, chart, left, right)
        elif chart == 'frozen_macro':
            value = self.macro_evaluate(op.macro, op.series, left, right)
            a = value['correlated_a_axial5']
            a[0] = a[0].positive_intersection(self.logh+self.logDlower)
        else:
            raise ValueError('Original micro or frozen macro chart required')
        return value

    def generic(self, ends, chart, left, right):
        op = self.owner(ends)
        background = self.query(ends, chart, left, right)
        proxy, raw, recovered, a = switch.recover_source(op.reference, op.axis, background)
        recovered.pop('source_frame_conditional_on_same_current_R100_background')
        recovered['same_original_whole_Z_micro_macro_source']=True
        return dict(background=background, full_signed_generic_source=recovered,
                    original_raw_source=raw, proxy=proxy, a=a)

    def seam(self, ends):
        op = self.owner(ends)
        back = self.query(ends, 'second_micro', 2, 2)
        next_ = self.query(ends, 'frozen_macro', (0,1), (0,1))
        rows = 0
        for name in ('phi', 'V'):
            for left, right in zip(back['fields'][name], next_['fields'][name]):
                delta = left-right
                lo, hi = ep(delta.coefficient)
                if not delta.zero and not lo <= 0 <= hi:
                    raise ValueError('Same original field seam has inconsistent source covers')
                rows += 1
        for name in moments.RATES:
            assert next_['histories'][name] is op.macro.inlet[name]
            for left, right in zip(back['histories'][name], next_['histories'][name]):
                delta = left-right
                lo, hi = ep(delta.coefficient)
                if not delta.zero and not lo <= 0 <= hi:
                    raise ValueError('Same original history seam has inconsistent source covers')
                rows += 1
        return dict(overlap_rows=rows,
            actual_micro_histories_are_the_live_macro_incoming=True,
            join_identity_from_same_source_definitions_and_history_ODE_uniqueness=True,
            overlaps_are_consistency_evidence_not_the_functional_join_proof=True)


def exported(value):
    """Consumer source rows only; exclude executable proxy objects."""
    return dict(source_background=value['background'],
        full_signed_generic_source=value['full_signed_generic_source'],
        original_raw_source=value['original_raw_source'])


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZMicroMacro()
        cells = []
        for ends in source.CELLS:
            op = owner.owner(ends)
            queries = [exported(owner.generic(ends, chart, left, right))
                       for chart, left, right in MICRO_CELLS]
            queries += [exported(owner.generic(ends, 'frozen_macro', left, right))
                        for left, right in MACRO_CELLS]
            cells.append(source.serialized(dict(exact_Z_cell=list(ends),
                whole_radial_source_cells=queries, actual_R100=op.actual_R100,
                actual_physical_R100=op.axis, exact_common_P0=op.reference.P0,
                live_micro_macro_seam=owner.seam(ends))))
            print('Live whole-Z micro/macro and six histories: '+str(ends), flush=True)
        result = dict(**{GATE:True}, source_family=owner.identity,
            exact_Z_domain=['-1','1'], exact_Z_partition=[list(cell) for cell in source.CELLS],
            source_cells=cells, original_positive_Dbar_theorem=owner.positive_theorem,
            original_source_bindings=owner.bindings,
            current_actual_micro_endpoint_fields_and_histories_are_macro_incoming=True,
            complete_original_signed_micro_macro_generic_inputs_available=True,
            no_ancestor_constructors_or_frame_packet_replay=True,
            **dict.fromkeys(OPEN,False), input_hashes=owner.hashes,
            execution_seconds=time.monotonic()-began,
            scope='Whole-Z actual microscopic and macro background source functions, '
                  'six own histories and full signed generic input graph through R100. '
                  'Original analytic Dbar lower is source-bound and logarithmic. '
                  'Current finite-N density integration, downstream whole-Z repair and recursion remain open.')
        with (HERE/NAME).open('wb') as target:
            with gzip.GzipFile(filename='',fileobj=target,mode='wb',mtime=0) as output:
                output.write(json.dumps(source.bridge._encode(result),indent=2).encode()+b'\n')
    print('Live whole-Z micro/macro sources and histories generated',flush=True)
    return result


if __name__ == '__main__':
    run()

"""Live whole-Z leading bridge inputs, without frame or ancestor reconstruction.

The admitted fixed point, pressure atoms and root disks are hydrated once.
Every requested real Z cell then rebuilds the coupled radial coefficients,
integrates the six core atoms and evaluates the original comparison and
prescribed-shear integrals. These are leading source functions, not finite-N
correction histories and not the temporal coefficient recursion.
"""
import gzip
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_bridge_macro_functions as macro
import lei_ren_part1_paper_compliant_actual_bridge_integrals as bridge
import lei_ren_part1_paper_compliant_core_coefficient_rebuild as rebuild_module
import lei_ren_part1_paper_compliant_anchored_axis_amplitude as amplitude_module
import lei_ren_part1_paper_compliant_core_physical_field as core_module
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

HERE, PREFIX, sha = macro.HERE, macro.PREFIX, macro.sha
ep, read, bind = macro.ep, macro.previous.read_interval, macro.previous.bind
NAME = PREFIX + 'current_original_whole_Z_bridge_source.json.gz'
RECEIPT = PREFIX + 'current_original_whole_Z_bridge_source_check.json'
GATE = 'current_original_whole_Z_live_leading_bridge_source_functions_installed'
CELLS = (('-1', '-.5'), ('-.5', '0'), ('0', '.5'), ('.5', '1'))
OPEN = ('current_finite_N_whole_Z_correction_prefix_through_Rc',
        'whole_Z_five_moment_terminal_repair_closed',
        'global_completed_tensor_admissibility', 'exact_heat_exterior_installed',
        'temporal_recursion', 'full_corrected_Navier_Stokes_solution')


def restore(c, row):
    """Decode full certified disks, never their provisional centers."""
    if isinstance(row, dict):
        if 'lower' in row and 'upper' in row:
            return read(c, row)
        if set(row) == {'real', 'imag'}:
            return c.mpc(restore(c, row['real']), restore(c, row['imag']))
        return {key: restore(c, value) for key, value in row.items()}
    if isinstance(row, list):
        return [restore(c, value) for value in row]
    return row


def serialized(row):
    if isinstance(row, IntervalTaylor):
        return serialized(list(row.coefficients))
    if isinstance(row, bridge.WidthPolynomial):
        return serialized(row.packet())
    if isinstance(row, dict):
        return {key: serialized(value) for key, value in row.items()}
    if isinstance(row, (tuple, list)):
        return [serialized(value) for value in row]
    return macro.serialized(row)


class OriginalWholeZBridgeSource:
    """An interval function API over [-1,1], not a saved-label adapter."""
    def __init__(self, dps=500):
        # This admission owner only reads checked data; its frame owner is
        # deliberately never called. None of the old source constructors run.
        self.admission = admission = macro.OriginalBridgeMacroFunctions(dps)
        self.c = c = admission.c
        self.hashes = dict(admission.hashes)
        self.source, self.family, self.datum = admission.source, admission.family, admission.datum
        self.records = dict(admission.records)
        for stem in ('core_uniform_bounds', 'physical_norm_family',
                     'core_coefficient_rebuild_check', 'core_integral_atoms_check',
                     'core_physical_field_check', 'core_transfer_check'):
            name = PREFIX + stem + '.json'
            row = json.loads((HERE / name).read_bytes())
            for path, digest in row.get('input_hashes', {}).items():
                bind(self.hashes, path, digest)
            bind(self.hashes, name, sha(name))
            self.records[stem] = row
            if stem.endswith('_check') and not row.get('all_passed'):
                raise ValueError('Accepted leading source prerequisite required: ' + stem)
            if row.get('implicit_source_sha256', self.source) != self.source:
                raise ValueError('Leading source identity differs: ' + stem)
        for stem in ('shared_analytic_tube',):
            name = 'lei_ren_part1_paper_' + stem + '.json'
            row = json.loads((HERE / name).read_bytes())
            for path, digest in row.get('input_hashes', {}).items():
                bind(self.hashes, path, digest)
            bind(self.hashes, name, sha(name))
            self.records[stem] = row
        major = self.records['core_transfer']
        norm = self.records['physical_norm_family']
        uniform = self.records['core_uniform_bounds']
        tube = self.records['shared_analytic_tube']
        analytic = major['analytic_core_family_sha256']
        if (norm['base_analytic_core_family_sha256'] != analytic
                or uniform['analytic_core_family_sha256'] != analytic
                or tube['analytic_core_family_sha256'] != analytic
                or norm['datum_enclosure_sha256'] != self.datum
                or not major['contraction_proved']
                or not norm['uniform_analytic_fixed_point_admission_extended']
                or not uniform['full_real_axis_normalized_swirl_positive']):
            raise ValueError('One admitted positive selected-Cstar analytic family required')
        scalar = lambda record, key: read(c, record[key])
        ps = admission.records['pressure_source']['compliant_source']
        delta = c.exp(scalar(ps['parameter_bounds'], 'log_delta'))
        q = scalar(major, 'scaled_map_Lipschitz_upper')
        if ep(q)[1] >= 1:
            raise ValueError('Original analytic contraction required')
        core = SimpleNamespace(ctx=c, source=self.source, family=self.family,
            delta=delta, logLambda=scalar(major, 'logLambda'),
            Lambda=scalar(major, 'Lambda'), epsilon=scalar(major, 'epsilon'),
            logP=scalar(major, 'logPstar'), logC=scalar(norm, 'selected_logCstar'),
            j=scalar(major, 'required_j'), h=scalar(tube, 'Xh_parameter'),
            correction=scalar(major, 'scaled_map_size_upper') / (1-q),
            phi_floor=scalar(uniform, 'actual_Phi_lower'),
            phi_ceiling=scalar(uniform, 'actual_Phi_upper'),
            Gbar=scalar(major, 'Gupper_in_logC_definition'),
            datum=admission.pressure, records=self.records, hashes=self.hashes)
        core.sigma = core.j / 500
        self.core = core
        rebuild = rebuild_module.CompliantCoreCoefficientRebuild.__new__(
            rebuild_module.CompliantCoreCoefficientRebuild)
        rebuild.core, rebuild.ctx, rebuild.hashes = core, c, self.hashes
        rebuild.eta = scalar(tube, 'complex_tube_radius')
        rebuild.phi_norm = scalar(major, 'Phi_ball_norm_upper')
        rebuild.psi_norm = scalar(major, 'Psi_ball_norm_upper')
        if ep(core.logC - (core.Lambda*core.Gbar + 2*core.logLambda + 1000))[0] <= 0:
            raise ValueError('Same admitted complex amplitude guard required')
        rebuild.S_bound = c.exp(-6*core.logLambda - 2000)
        if ep(rebuild.S_bound)[0] <= 0:
            raise ValueError('Nonzero positive nonlinear swirl source bound required')
        rebuild.model_tail_cache = {}
        self.rebuild = rebuild
        amplitude = amplitude_module.CompliantAnchoredAxisAmplitude.__new__(
            amplitude_module.CompliantAnchoredAxisAmplitude)
        amplitude.core, amplitude.ctx, amplitude.rebuild = core, c, rebuild
        amplitude.hashes = self.hashes
        amplitude.params = dict(j=core.j, delta=core.delta, sigma=core.sigma,
                                a=(9-core.delta)/2)
        amplitude.radius = c.mpf('1e-180')
        amplitude.roots = restore(c, admission.records['anchored_axis_amplitude']['root_certificates'])
        amplitude.rows, amplitude.cache = {}, {}
        if not amplitude.roots['all_six_poles_certified']:
            raise ValueError('Full certified source pole disks required')
        self.amplitude = amplitude
        original = admission.records['actual_bridge_integrals']
        self.weight = weight = scalar(original, 'axial_jet_weight')
        atoms = bridge.BridgeCoreAtoms.__new__(bridge.BridgeCoreAtoms)
        atoms.core, atoms.ctx, atoms.rebuild = core, c, rebuild
        atoms.axial_weight, atoms.cache, atoms.hashes = weight, {}, self.hashes
        self.atoms = atoms
        logRa = c.ln(4) - core.logLambda
        self.logRa = logRa
        self.logh = scalar(original, 'source_log_hb_enclosure')
        self.cap = c.mpf('1e-50000')
        if ep(self.logh)[1] >= ep(c.ln(self.cap))[0]:
            raise ValueError('Same true source hb must satisfy whole-axis remainder cap')
        geometry = SimpleNamespace(r=c.exp(logRa), logh=self.logh,
                                   inputs=self.inputs, source=self.source, family=self.family)
        comparison = bridge.BridgeComparison.__new__(bridge.BridgeComparison)
        comparison.atoms, comparison.core, comparison.ctx = atoms, core, c
        comparison.bridge = SimpleNamespace(r=geometry.r, logh=self.logh, bridge=geometry)
        comparison.source = admission.records['comparison_point_integrals']['original_comparison_source']
        comparison.weight, comparison.cap = weight, self.cap
        comparison.cache, comparison.weight_cache, comparison.hashes = {}, {}, self.hashes
        self.comparison = comparison
        actual = bridge.CompliantActualBridgeIntegrals.__new__(bridge.CompliantActualBridgeIntegrals)
        actual.comparison, actual.bridge, actual.core, actual.ctx = comparison, geometry, core, c
        actual.direction_context = SimpleNamespace(ctx=c, delta=delta)
        actual.weight, actual.cap, actual.logh = weight, self.cap, self.logh
        actual.Y = c.ln(100)-logRa
        actual.pressure_log = 2*core.logP
        # This is a bounding log enclosure for the original positive function,
        # not a constant substitute for G(Z); each flow below uses actual G.
        actual.swirl_log = c.mpf([ep(-2*core.logC-2*core.Lambda*core.Gbar)[0],
                                 ep(-2*core.logC)[1]])
        actual.hashes, actual.bindings = self.hashes, admission.bindings
        actual.proofs, actual.cache, actual.pulse_cache = [], {}, {}
        self.actual = actual
        self.cache, self.input_cache = {}, {}
        for module in (macro, bridge, rebuild_module, amplitude_module, core_module):
            name = Path(module.__file__).name
            bind(self.hashes, name, sha(name))
        bind(self.hashes, Path(__file__).name, sha(Path(__file__).name))

    def cell(self, ends):
        c = self.c
        if not isinstance(ends, (tuple, list)) or len(ends) != 2:
            raise ValueError('Explicit ordered Z cell endpoints required')
        left, right = c.mpf(ends[0]), c.mpf(ends[1])
        if ep(left)[0] != ep(left)[1] or ep(right)[0] != ep(right)[1]:
            raise ValueError('Exact cell endpoints required')
        lo, hi = ep(left)[0], ep(right)[1]
        if not -1 <= lo <= hi <= 1:
            raise ValueError('Source cell must lie in [-1,1]')
        return c.mpf([lo, hi])

    def inputs(self, Z):
        c = self.c
        z = c.mpf(Z)
        if ep(z)[0] < -1 or ep(z)[1] > 1:
            raise ValueError('Source Z domain required')
        if z._mpi_ in self.input_cache:
            return self.input_cache[z._mpi_]
        gradient = self.amplitude.gradient_jets(z, 5)
        def ratios(multiplier):
            rows = [c.mpf(1)]
            for n in range(1, 7):
                rows.append(sum((-multiplier*self.core.Lambda*gradient[j]*rows[n-1-j]
                                 for j in range(n)), c.mpf(0))/n)
            return [value*math.factorial(n) for n, value in enumerate(rows)]
        pressure = self.admission.pressure.normalized_jets(ep(z), 6)
        result = dict(p0=IntervalTaylor(c, [c.mpf(ep(value)) for value in
                      pressure['normalized_pressure_coefficients']]),
                      F0_ratios=ratios(1), F0_squared_ratios=ratios(2))
        self.input_cache[z._mpi_] = result
        return result

    def owner(self, ends):
        """Build one live MacroFlow from fresh implicit core and micro inputs."""
        Z = self.cell(ends)
        if Z._mpi_ in self.cache:
            return self.cache[Z._mpi_]
        c = self.c
        amplitude = self.amplitude.evaluate(Z)
        prepared = self.actual.prepare(Z)
        core = self.actual.packet(Z, 0, 'first')
        angular, axial = self.actual.contributions(prepared, 'second', c.mpf(2))
        flow = macro.MacroFlow(c, self.logh, 2*self.core.logP,
            2*amplitude['logF0'], self.logRa, c.ln(100)-self.logRa, self.weight)
        def terms(rows, part=None):
            pieces = []
            for row in rows:
                jet = IntervalTaylor(c, row['normalized_axial_coefficients'])
                pieces.append(flow.term(jet, row['width_power'], part))
                powers = (row['width_power']+1, int(part == 'pressure'),
                          int(part == 'swirl'), 0, 0)
                pieces.append(flow.error_rows(flow.factor(powers) *
                    row['nonlinear_prefix_error_weighted_norm_per_next_hb']))
            return flow.add(*pieces)
        ell = terms(angular)
        dV = flow.add(*[terms(axial[part], part) for part in macro.PARTS])
        modes = prepared['modes']
        flow.set_sources([row.truncate(5) for row in modes['D_over_R']],
            {part: [row.truncate(5) for row in modes['drive_'+part]]
             for part in macro.PARTS}, prepared['quotient'], ell,
            prepared['data']['phi0'].truncate(5), prepared['data']['V0'].truncate(5), dV)
        proof = dict(Z=Z, implicit_source_sha256=self.source,
            actual_five_defect_family_sha256=self.family, datum_enclosure_sha256=self.datum,
            comparison_source_namespace=self.comparison.source['source_namespace'],
            fresh_coupled_radial_degree=24, fresh_axial_depth=6,
            fresh_actual_core=core, actual_anchored_amplitude=amplitude,
            comparison_micro_endpoint=prepared['end'],
            original_micro_inlet_angular_terms=angular, original_micro_inlet_axial_terms=axial,
            actual_G_not_Gbar_used_in_flow=True,
            positive_nonlinear_swirl_bound_retained=self.rebuild.S_bound,
            original_pressure_Z0_through_Z6=list(self.inputs(Z)['p0'].coefficients),
            saved_label_packets_or_old_finite_radial_rows_used=False,
            archived_N1024_correction_rows_reused=False,
            original_ancestor_constructors_or_producers_executed=False,
            source_interval_is_not_interpolation=True)
        self.cache[Z._mpi_] = flow, proof
        return flow, proof

    def evaluate(self, ends):
        flow, proof = self.owner(ends)
        Z = proof['Z']
        value = flow.evaluate((1, 1))
        actual_R100 = self.actual.packet(Z, 1, 'macro')
        return dict(source_proof=proof, live_macro_R100=value,
            actual_leading_R100_own_histories=actual_R100,
            one_live_source_basis_and_ledger=True,
            positive_formal_F0_squared_retained=True,
            leading_source_not_finite_N_incoming=True,
            **dict.fromkeys(OPEN, False))


def run():
    began = time.monotonic()
    with mp.workdps(540):
        provider = OriginalWholeZBridgeSource()
        packets = []
        for ends in CELLS:
            packets.append(serialized(provider.evaluate(ends)))
            print('Fresh whole-Z original leading bridge cell: '+str(ends), flush=True)
        result = dict(**{GATE: True}, implicit_source_sha256=provider.source,
            actual_five_defect_family_sha256=provider.family, datum_enclosure_sha256=provider.datum,
            exact_Z_domain=['-1', '1'], exact_Z_partition=[list(cell) for cell in CELLS],
            source_cells=packets, interval_function_API_accepts_arbitrary_source_subcells=True,
            core_rows_rebuilt_not_loaded_from_frames=True,
            original_signed_micro_and_macro_leading_source_integrals_retained=True,
            leading_actual_six_history_enclosures_available_on_whole_Z=True,
            no_old_ancestor_constructors_or_producers_executed=True,
            **dict.fromkeys(OPEN, False), input_hashes=provider.hashes,
            execution_seconds=time.monotonic()-began,
            scope='Live original leading source functions over the complete real axial domain. '
                  'Fresh coupled radial rows, differentiated nonlinear tails, analytic P0, '
                  'certified anchored amplitude and original prescribed-shear integrals. '
                  'Current finite-N correction transport and terminal five-moment repair remain open.')
        payload = json.dumps(bridge._encode(result), indent=2).encode('utf8')+b'\n'
        with (HERE / NAME).open('wb') as target:
            with gzip.GzipFile(filename='', fileobj=target, mode='wb', mtime=0) as output:
                output.write(payload)
    print('Live whole-Z bridge source provider generated', flush=True)
    return result


if __name__ == '__main__':
    run()

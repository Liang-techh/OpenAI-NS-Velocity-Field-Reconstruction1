"""Genuine signed inverse functions and C1 transport on the actual Rm patch.

The original x=R/Rm function and dx/x measure are retained. The existing
correction-only Rm incoming is transported, separately from the background.
Full inner-to-Rc defect functions and terminal repair controls remain open.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_signed_loop_functions as previous

current = previous.current
patch = current.upstream.upstream.upstream.upstream
HERE, PREFIX, sha, bind, ep = previous.HERE, previous.PREFIX, previous.sha, previous.bind, previous.ep
CELLS, RATES, C0, Z, OPEN = previous.CELLS, previous.RATES, previous.C0, previous.Z, previous.OPEN
PARTITION = patch.PARTITION
NAME = PREFIX + 'current_original_whole_Z_Rm_Rh_signed_functions.json.gz'
RECEIPT = PREFIX + 'current_original_whole_Z_Rm_Rh_signed_functions_check.json'
GATE = 'current_original_whole_Z_actual_Rm_patch_signed_functions_and_C1_transport_through_Rh_installed'
serialized, encode = previous.serialized, previous.encode


def coordinate_record(value):
    return 'Rh' if value == 'Rh' else list(value)


def function_digest(value):
    return hashlib.sha256(json.dumps(encode(serialized(value)), sort_keys=True,
        separators=(',', ':')).encode()).hexdigest()


class WholeZRmRhSignedFunctions(previous.WholeZSignedLoopFunctions):
    def __init__(self, dps=500):
        super().__init__(dps)
        checked = json.loads((HERE / previous.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(previous.GATE) \
                or checked['source_family'] != self.identity or checked['candidate_N'] != self.N:
            raise ValueError('Checked current whole-Z signed inverse function backend required')
        for name, digest in checked['input_hashes'].items():
            bind(self.hashes, name, digest)
        bind(self.hashes, previous.RECEIPT, sha(previous.RECEIPT))
        bind(self.hashes, Path(previous.__file__).name, sha(Path(previous.__file__).name))
        bind(self.hashes, Path(__file__).name, sha(Path(__file__).name))
        self.patch_source = self.source.source.source.source.source
        if type(self.patch_source) is not patch.WholeZRmPatchFiniteN:
            raise ValueError('Actual current whole-Z Rm patch owner chain required')
        self.patch_radius = {}

    def source_query(self, ends, chart, left, right=None):
        if chart == 'actual_patch':
            return self.patch_source.query(ends, left, None if right == left else right)
        return super().source_query(ends, chart, left, right)

    def radius_query(self, ends, chart, left, right=None):
        if chart != 'actual_patch':
            return super().radius_query(ends, chart, left, right)
        right = left if right is None else right
        key = tuple(ends)
        if key not in self.patch_radius:
            self.patch_radius[key] = previous.phase.RmRadiusPhase(
                self.owner(ends), self.identity['actual_five_defect_family_sha256'], self.identity)
        mapper = self.patch_radius[key]
        if left == right:
            if left != 'Rh':
                previous.exact(left)
            result = mapper.point(left, self.N, decimal_digits=80)
            full = False
        else:
            if left == 'Rh':
                raise ValueError('Actual patch cell must start below Rh')
            previous.exact(left)
            if right != 'Rh':
                previous.exact(right)
            result = mapper.cell(left, right, self.N)
            full = result['whole_period_cover']
        return dict(source_identity=self.identity, candidate_N=self.N, actual_chart=chart,
            exact_left=coordinate_record(left), exact_right=coordinate_record(right),
            actual_original_Rm_patch_phase_binding=result, actual_phase_boxes=result['phase_boxes'],
            full_period=full, actual_phase_Z_exact_zero=True,
            actual_same_source_Rm_radius=self.owner(ends).Rm_factor,
            global_phase_not_restarted=True, phase_not_supplied_as_free_angle=True,
            adaptive_analytic_interval_digits=max(260, 80 + (self.N.bit_length()*30103+99999)//100000 + 60),
            exact_original_global_phase='frac(N*(logRm+log(x)-logRa-hb*s_c/2))',
            actual_radial_measure='dy=dx/x', symbolic_Rh_endpoint='x=exp(1)')

    def query(self, ends, chart, left, right=None):
        result = super().query(ends, chart, left, right)
        if chart == 'actual_patch':
            result['exact_left'] = coordinate_record(left)
            result['exact_right'] = coordinate_record(left if right is None else right)
        return result

    def local_integral(self, ends, chart, left, right):
        if chart != 'actual_patch':
            return super().local_integral(ends, chart, left, right)
        source = self.query(ends, chart, left, right)
        op = self.owner(ends)
        f, c = op.flow, op.c
        width, _, _, _, _ = patch.patch_weights(f, left, right, next(iter(RATES.values())))
        if ep(width)[0] <= 0:
            raise ValueError('Strictly positive actual logarithmic patch interval required')
        density = source['actual_signed_nonlinear_density_C0_Z']
        integrals, masses = {}, {}
        for name, rate in RATES.items():
            mass = current.original.incoming.weighted.terminal.local.positive_kernel_mass(c, width, rate)
            integrals[name] = [density[part][name]*mass for part in ('kernels', 'Z_derivatives')]
            masses[name] = mass
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            actual_chart=chart, exact_left=coordinate_record(left), exact_right=coordinate_record(right),
            exact_common_P0_axial5=op.P0, actual_source_function=source,
            actual_physical_log_width=width, actual_positive_own_rate_masses=masses,
            actual_signed_local_integral_C0_Z=integrals,
            integral_recipe='int_left^right exp(-rate*log(x_right/x))*signed_density(x,Z,N) dx/x',
            live_inverse_function_extension_under_Z_independent_integral=True,
            directed_rectangle_function_enclosure_not_selected_integral_value=True,
            physical_Jacobian_applied_once=True, incoming_correction_not_supplied_or_reset=True,
            actual_full_prefix_C1_defect_functions_installed=False,
            actual_terminal_controls_installed=False, **dict.fromkeys(OPEN, False))

    def transport(self, ends):
        op = self.owner(ends)
        f, c = op.flow, op.c
        inlet = self.patch_source.inlet(ends)
        incoming = inlet['actual_Rm_correction_C0_Z']
        local = {name: [f.scalar(0), f.scalar(0)] for name in RATES}
        cell_records = []
        for left, right in zip(PARTITION, PARTITION[1:]):
            integral = self.local_integral(ends, 'actual_patch', left, right)
            factors, contribution = {}, {}
            for name, rate in RATES.items():
                width, suffix, mass, decay, tail = patch.patch_weights(f, left, right, rate)
                contribution[name] = [row*tail for row in integral['actual_signed_local_integral_C0_Z'][name]]
                for i in range(2):
                    local[name][i] += contribution[name][i]
                factors[name] = dict(physical_log_width=width, suffix_log_width=suffix,
                    original_positive_full_mass=mass, original_cell_decay=decay,
                    original_suffix_decay=tail, own_rate=str(rate))
            source = integral['actual_source_function']
            cell_records.append(dict(exact_left=coordinate_record(left), exact_right=coordinate_record(right),
                current_signed_function_and_integral_sha256=function_digest(integral),
                actual_signed_primitive_C0_Z_phi=source['actual_signed_primitive_C0_Z_phi'],
                actual_signed_nonlinear_density_C0_Z=source['actual_signed_nonlinear_density_C0_Z'],
                actual_signed_local_integral_C0_Z=integral['actual_signed_local_integral_C0_Z'],
                actual_signed_contribution_at_Rh_C0_Z=contribution, original_own_rate_factors=factors,
                signed_inverse_alternative_count=len(source['actual_signed_inverse_function_alternatives']),
                actual_global_phase_full_period_cover=source['actual_source_geometry']['full_period'],
                full_function_packet_available_from_live_oracle_and_replay_digest=True))
            print('Current signed Rm/Rh integral: '+str(ends)+' '+str(left)+' '+str(right), flush=True)
        memory = {name: patch.long.kernel_weight(f, c.mpf(1), c.mpf(0), rate)[1]
            for name, rate in RATES.items()}
        inherited = {name: [row*memory[name] for row in incoming[name]] for name in RATES}
        outgoing = {name: [inherited[name][i]+local[name][i] for i in range(2)] for name in RATES}
        terminal = self.query(ends, 'actual_patch', 'Rh')
        background = {name: list(rows[:2]) for name, rows in terminal['actual_original_source_packet']
            ['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete = {name: [background[name][i]+outgoing[name][i] for i in range(2)] for name in RATES}
        return dict(source_identity=self.identity, exact_Z_cell=list(ends), candidate_N=self.N,
            exact_common_P0_axial5=op.P0, exact_Rm_radius=op.Rm_factor,
            exact_Rh_radius=op.Rm_factor*c.exp(1), actual_Rm_incoming_binding=inlet,
            exact_total_log_width=1, exact_length_identity='sum log(x_right/x_left)=log(exp(1))=1',
            actual_signed_patch_cells=cell_records, actual_Rm_correction_incoming_C0_Z=incoming,
            original_incoming_own_rate_memory=memory, actual_retained_incoming_C0_Z=inherited,
            actual_signed_local_driver_at_Rh_C0_Z=local, actual_Rh_signed_correction_C0_Z=outgoing,
            actual_original_Rh_background_C0_Z=background, actual_Rh_signed_complete_own_history_C0_Z=complete,
            actual_original_Rh_signed_source_function=terminal,
            current_source_signed_inverse_used_in_every_cell=True,
            interval_caps_not_substituted_for_defining_functions=True,
            correction_only_incoming_distinct_from_background=True,
            actual_nonzero_incoming_memory_retained=True, actual_rate_zero_pressure_memory_one=True,
            actual_full_prefix_C1_defect_functions_installed=False,
            actual_terminal_controls_installed=False, **dict.fromkeys(OPEN, False))


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZRmRhSignedFunctions()
        points, transports = [], []
        for ends in CELLS:
            for point in ((1, 1), (3, 2), 'Rh'):
                points.append(serialized(owner.query(ends, 'actual_patch', point)))
            transports.append(serialized(owner.transport(ends)))
        report = dict(**{GATE: True}, source_family=owner.identity, candidate_N=owner.N,
            exact_Z_partition=CELLS, exact_patch_partition=serialized(PARTITION),
            actual_signed_patch_point_functions=points, actual_signed_patch_transports=transports,
            actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False, **dict.fromkeys(OPEN, False), input_hashes=owner.hashes,
            execution_seconds=time.monotonic()-began,
            scope='Current whole-Z actual Rm/Rh patch inverse-defined signed C0/Z/phase functions '
                'and local C1 integral transport with dy=dx/x, true suffix kernels and retained same-N '
                'correction-only incoming. Full inner-to-Rc defect functions, terminal controls and '
                'global N/heat/cone/temporal recursion remain unfinished.')
        (HERE / NAME).write_bytes(gzip.compress((json.dumps(encode(report),
            separators=(',', ':'))+'\n').encode(), compresslevel=6, mtime=0))
    print('Current whole-Z signed Rm patch functions and genuine C1 transport reach Rh', flush=True)
    return report


if __name__ == '__main__':
    run()

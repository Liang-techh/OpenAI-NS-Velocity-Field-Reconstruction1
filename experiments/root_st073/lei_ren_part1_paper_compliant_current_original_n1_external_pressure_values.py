"""Compose checked original pressure values on the whole [Rin,infinity).

The lower collar keeps its signed source models and positive scale ratio.
No ancestor field is rebuilt and no affine approximation selects a field.
Inner cut moments and full derivative accuracy remain separate open layers.
"""
import copy
import gzip
import hashlib
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source as upper
import lei_ren_part1_paper_compliant_current_original_n1_unchanged_axial_pressure_models as lower
import lei_ren_part1_paper_compliant_current_original_Rp_amplitude_binding as amplitude

HERE, PREFIX, sha, ep = upper.HERE, upper.PREFIX, upper.sha, upper.ep
record, interval, source_precision = upper.record, upper.interval, upper.source_precision
NAME = PREFIX+'current_original_n1_external_pressure_values.json.gz'
RECEIPT = PREFIX+'current_original_n1_external_pressure_values_check.json'
GATE = 'current_original_n1_Rin_to_infinity_pressure_values_installed'
OPEN = upper.OPEN
NORM = PREFIX+'physical_norm_family.json'
CORE = PREFIX+'core_transfer.json'


def snapshot(value):
    return hashlib.sha256(json.dumps(record(value), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def absolute_upper(c, value):
    return c.mpf(max(abs(x) for x in ep(value)))


def checked_bundle():
    """Check each inherited dependency byte once across the three archives."""
    hashes = {}; archives = []
    for module, gates in ((upper, (upper.GATE,)), (lower, (lower.GATE,)),
                          (amplitude, amplitude.GATES)):
        contents = (HERE/module.NAME).read_bytes()
        raw = json.loads(gzip.decompress(contents) if module.NAME.endswith('.gz') else contents)
        receipt = json.loads((HERE/module.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or any(not receipt.get(gate) for gate in gates):
            raise ValueError('Accepted original pressure/amplitude archive required '+module.NAME)
        if raw['source_family'] != receipt['source_family']:
            raise ValueError('Archive source family differs '+module.NAME)
        if any(receipt.get(key, False) for key in module.OPEN):
            raise ValueError('Inherited open physical gate differs '+module.NAME)
        for name, digest in receipt['input_hashes'].items():
            if name in hashes and hashes[name] != digest:
                raise ValueError('Conflicting original source '+name)
            hashes[name] = digest
        for name, digest in raw.get('input_hashes', {}).items():
            if hashes.get(name) != digest:
                raise ValueError('Raw original source binding differs '+name)
        for name in (module.NAME, module.RECEIPT, Path(module.__file__).name):
            hashes[name] = sha(name)
        archives.append(raw)
    if any(raw['source_family'] != archives[0]['source_family'] for raw in archives):
        raise ValueError('Same original cut/collar/pulse family required')
    if not all(name in hashes for name in (NORM, CORE)):
        raise ValueError('Original core and selected Cstar must already be bound')
    for name, digest in hashes.items():
        if sha(name) != digest:
            raise ValueError('Changed inherited original source '+name)
    hashes[Path(__file__).name] = sha(Path(__file__).name)
    return archives, hashes


def ratio_log(c, power, loghb, logP, logLambda, logRc, logUp, M):
    """Log(hb*Pstar^2*Lambda^-p*M / (Phi*A*Rv^p)).

    Original logE and Rv share mu. Combine their inverse-mu terms exactly:
    13/2 + 1 - 13*p = 15/2 - 13*p; never subtract unrelated huge boxes.
    """
    if power not in (1, 2) or ep(M)[0] <= 0:
        raise ValueError('Positive original pressure scale required')
    logmu = -4*logP-c.ln(1000)
    invmu = c.exp(-logmu)
    phi = c.mpf('10.99')*c.sqrt(2*c.pi)*(c.exp(5)-c.exp(3))
    two_thirds = c.mpf(2)/3
    return (loghb+logP-power*(logLambda+logRc)-logUp+c.ln(M)-c.ln(phi)+13+2*power
            +(60*power+c.mpf(4)/3)*logmu+(c.mpf('7.5')-13*power)*invmu
            +3*c.exp(-two_thirds*(c.ln(4)+logmu))-two_thirds*c.ln(4)+c.ln(6)/2)


def scales(c, upper_data, lower_data, amplitude_data):
    core = json.loads((HERE/CORE).read_bytes())
    norm = json.loads((HERE/NORM).read_bytes())
    loghb, logP2, logLambda = [interval(c, row) for row in lower_data['original_source_log_bases']]
    logP = interval(c, core['logPstar'])
    if not (ep(logP2)[0] <= ep(2*logP)[1] and ep(logP2)[1] >= ep(2*logP)[0]):
        raise ValueError('Lower pressure and original Pstar differ')
    if interval(c, core['logLambda'])._mpi_ != logLambda._mpi_:
        raise ValueError('Lower pressure and original Lambda differ')
    logRc = c.ln(110)+10*interval(c, norm['selected_logCstar'])+11*c.exp(40)+124
    logUp = interval(c, amplitude_data['directed_logU0'])
    units = upper_data['original_physical_units']
    if (units['source_A'] != 'Pstar*U0*exp(logE-13/(2*mu)-13)' or
        units['source_Rv_formula'] != 'Rp*exp(13/mu)' or
        units['source_logRp_formula'] != 'logRc-2-60*log(mu)'):
        raise ValueError('Accepted original pulse unit identity required')
    geometry = upper_data['source_data']['geometry']
    if geometry['original_log_r_minus_formula'] != 'log(4)-4*logP-1000+hbB*sc/2':
        raise ValueError('Accepted common original rminus geometry required')
    rows = {}
    for key, power in (('m13', 1), ('m15', 2)):
        normalization = lower_data['normalizations'][key]
        if normalization['reference_source_exponents'] != [1, 0, 1, -power]:
            raise ValueError('Actual collar source exponents changed')
        M = interval(c, normalization['directed_reference_coefficient_norm'])
        logarithm = ratio_log(c, power, loghb, logP, logLambda, logRc, logUp, M)
        if ep(logarithm)[1] >= -1000:
            raise ValueError('Original positive lower scale ratio not covered')
        rows[key] = dict(original_source_exponents=[1, 0, 1, -power],
                         coefficient_norm=M, exact_positive_ratio_log=logarithm,
                         positive_ratio_enclosure=c.mpf([0, ep(c.exp(-1000))[1]]),
                         positive_ratio_is_bound_not_selected_zero=True)
    return dict(fields=rows, loghb=loghb, logP=logP, logLambda=logLambda, logRc=logRc,
                logUp=logUp, original_basis_names=['hb', 'F0_anchor', 'Pstar_squared', 'Lambda'],
                common_rminus=dict(lower_phase='sc/2', upper_phase='sc/2',
                    Rin_phase='sc/4', same_original_family=True, shared_endpoint_counted_once=True),
                exact_log_scale_formula='log(hb*Pstar^2*Lambda^(-p)*M/(Phi*A*Rv^p))',
                inverse_mu_coefficient='15/2-13*p', no_physical_scale_materialization=True)


def compose_models(c, upper_data, lower_data, scale_data):
    """Refine the accepted outer grid to the lower grid; retain both errors."""
    result = []
    for i, cell in enumerate(lower_data['cells']):
        left, right = [c.mpf(x) for x in cell['exact_Z_cell']]
        center = c.mpf(cell['exact_center']); half = (right-left)/2
        outer_index = i//2
        parent = upper_data['models'][outer_index]
        if ep(left)[0] < ep(interval(c, parent['left']))[0] or ep(right)[1] > ep(interval(c, parent['right']))[1]:
            raise ValueError('Exact common model partition required')
        shift = center-interval(c, parent['center']); fields = {}
        for key in ('m13', 'm15'):
            old = parent['fields'][key]; local = cell['models'][key]
            ratio = scale_data['fields'][key]['positive_ratio_enclosure']
            a, b = [interval(c, value) for value in local['affine_coefficients_Z0_Z1']]
            A = interval(c, old['center_intercept'])+interval(c, old['center_slope'])*shift+ratio*a
            B = interval(c, old['center_slope'])+ratio*b
            lower_remainder = ratio*interval(c, local['Taylor_remainder_abs_upper'])
            remainder = interval(c, old['linear_Taylor_remainder'])+lower_remainder
            source_error = interval(c, old['full_nonzero_source_error'])
            ambiguity = (A.b-A.a)/2+(B.b-B.a)*half/2
            error = absolute_upper(c, remainder+source_error+ambiguity)
            fields[key] = dict(center_intercept=A, center_slope=B,
                inherited_outer_Taylor_remainder=interval(c, old['linear_Taylor_remainder']),
                lower_Taylor_remainder=lower_remainder, linear_Taylor_remainder=remainder,
                coefficient_ambiguity=ambiguity, full_nonzero_source_error=source_error,
                absolute_error=error, absolute_accuracy_certified=ep(error)[1]<c.mpf('.001'),
                lower_signed_coefficients_retained=True,
                full_pressure_derivative_accuracy_not_claimed=True)
        result.append(dict(left=left, right=right, center=center, half_width=half,
            lower_model_index=i, upper_model_index=outer_index, fields=fields,
            approximation='interval affine enclosure; does not select the physical field'))
    if len(result) != 512 or len(upper_data['models']) != 256:
        raise ValueError('Accepted common full-Z model counts changed')
    return result


class CurrentOriginalN1ExternalPressureValues:
    @source_precision
    def __init__(self, require_checked=True):
        self.c = MPIntervalContext(); self.c.dps = 540
        archives, self.hashes = checked_bundle()
        self.upper, self.lower, self.amplitude = archives
        self.family = copy.deepcopy(self.upper['source_family'])
        self.scale_data = scales(self.c, *archives)
        self.models = compose_models(self.c, self.upper, self.lower, self.scale_data)
        self.acceptance_loaded = False
        self.bindings = {name:upper.pulse.end.coefficients.function_binding(fn) for name, fn in (
            ('same_source_unit_ratio', ratio_log), ('original_scale_binding', scales),
            ('common_Z_model_composition', compose_models), ('checked_archive_bundle', checked_bundle))}
        if require_checked:
            receipt = json.loads((HERE/RECEIPT).read_bytes())
            raw = json.loads(gzip.decompress((HERE/NAME).read_bytes()))
            if (not receipt.get('all_passed') or not receipt.get(GATE) or
                receipt['source_family'] != self.family or any(receipt[k] for k in OPEN)):
                raise ValueError('Checked complete external pressure values required')
            for name, digest in receipt['input_hashes'].items():
                if sha(name) != digest: raise ValueError('Changed complete pressure source '+name)
            if raw['models'] != record(self.models) or raw['source_scales'] != record(self.scale_data):
                raise ValueError('Accepted complete external pressure models differ')
            self.hashes.update(receipt['input_hashes']); self.acceptance_loaded = True
        self._packets = {}
        self._state = self.state()

    def state(self):
        return snapshot((self.family, self.hashes, self.upper, self.lower, self.amplitude,
                         self.scale_data, self.models, self.bindings))

    def assert_graph(self):
        if self.state() != self._state or sha(Path(__file__).name) != self.hashes[Path(__file__).name]:
            raise ValueError('Complete external pressure source changed')

    def _make_view(self, Z):
        c = self.c; z = c.mpf(Z)
        if ep(z)[0] < -1 or ep(z)[1] > 1: raise ValueError('Original real Z domain [-1,1] required')
        values = {key:[] for key in ('m13','m15')}; active = []
        for i, model in enumerate(self.models):
            lo = max(ep(z)[0], ep(model['left'])[0]); hi = min(ep(z)[1], ep(model['right'])[1])
            if lo > hi: continue
            active.append(i); offset = c.mpf([lo,hi])-model['center']
            for key, field in model['fields'].items():
                error = absolute_upper(c, field['linear_Taylor_remainder']+field['full_nonzero_source_error'])
                values[key].append(field['center_intercept']+field['center_slope']*offset+c.mpf([-ep(error)[1],ep(error)[1]]))
        if not active: raise ValueError('Missing original pressure model coverage')
        hulls = {key:c.mpf([min(ep(v)[0] for v in rows),max(ep(v)[1] for v in rows)]) for key,rows in values.items()}
        return dict(Z=z, source_family=self.family, radial_scope='Rin to infinity',
            full_pressure_values=hulls, active_models=active, original_physical_units=self.upper['original_physical_units'],
            inner_cut_pressure_moments_not_included=True, no_affine_midpoint_selected_as_field=True,
            positive_lower_scale_ratio_retained=True, full_pressure_derivative_accuracy_not_claimed=True,
            **{GATE:self.acceptance_loaded}, **dict.fromkeys(OPEN,False))

    @source_precision
    def evaluate(self, Z):
        self.assert_graph(); view = self._make_view(Z); packet = object()
        self._packets[id(packet)] = (packet, view, snapshot(view)); return packet

    @source_precision
    def report(self, packet):
        self.assert_graph(); entry = self._packets.get(id(packet))
        if entry is None or entry[0] is not packet or entry[2] != snapshot(entry[1]):
            raise ValueError('Live complete external pressure packet required')
        return upper.pulse.end.copy_view(self.c, entry[1])


@source_precision
def run():
    began = time.monotonic(); owner = CurrentOriginalN1ExternalPressureValues(False)
    owner.assert_graph()
    samples = {z:owner._make_view(z) for z in ('-1','-.639','0','.211','.371','.639','1')}
    owner.assert_graph()
    raw = dict(source_family=owner.family, models=record(owner.models), samples=record(samples),
        source_scales=record(owner.scale_data), original_physical_units=owner.upper['original_physical_units'],
        source_bindings=owner.bindings, input_hashes=owner.hashes,
        radial_scope='Rin to infinity', inner_cut_pressure_moments_not_included=True,
        ancestor_math_suites_not_rerun=True, inherited_models_reused=True,
        **{GATE:False}, **dict.fromkeys(OPEN,False), execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_COMPLETE_EXTERNAL_PRESSURE_VALUES_COMPUTED',len(owner.models),flush=True)
    return owner, raw

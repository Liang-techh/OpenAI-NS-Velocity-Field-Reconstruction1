"""Actual original 2Rc-to-Rp quiet velocity/history/pressure enclosures."""
import copy
from fractions import Fraction
import gzip
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_core_velocity import OriginalCoreVelocity
from .original_compact_n1_pressure import OriginalCompactN1PressureVelocity
from .original_quiet_scalar import OriginalQuietScalar


class OriginalQuietNumeric:
    def __init__(self, *, precision=85, amplitude_precision=110, directory=None):
        self.compact = OriginalCompactN1PressureVelocity(precision=precision,
                    amplitude_precision=amplitude_precision, directory=directory)
        f = self.compact
        self.ctx, self.iv, self.J = f.ctx, f.iv, f.J
        self.family, self.directory = f.family, f.numeric.directory
        self.amplitude, self.owner = f.amplitude, f.numeric.owner
        self.outer = f.numeric.outer
        self.hashes, self.code_hashes = dict(f.hashes), dict(f.code_hashes)
        path = self.directory / (PREFIX + 'compliant_current_original_compact_n1_pressure_numeric.json')
        accepted = json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family'] != self.family:
            raise ValueError('Accepted original compact velocity and P1 inlet required')
        for name, digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Changed actual compact pressure implementation: '+name)
            self.code_hashes[name] = digest
        self.hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.inlets = {p['Z']: p for p in accepted['samples'] if p['native_compact_x'] == '2'}
        specifications = (
            ('current_original_C3_power_to_Rp', ('current_original_actual_repaired_exit_C3_quiet_power_and_pulse_frame_installed',)),
            ('current_original_true_quiet_Rp_C5_pressure_C3_functions', ('current_original_true_2Rc_to_Rp_quiet_C5_pressure_C3_functions_installed',
                'current_original_true_quiet_K_matches_original_native_source_installed')),
            ('current_original_n1_quiet_pressure_source', ('current_original_n1_quiet_source_scaled_signed_pressure_moments_installed',)),
            ('current_original_n1_repaired_own_C3_terminal_functions', ('current_original_n1_mass_pressure_C3_QV_C2_tail_functions_installed',)))
        raw = {}
        for stem, gates in specifications:
            receipt_path = self.directory / (PREFIX+'compliant_'+stem+'_check.json')
            receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or not all(receipt[k] for k in gates) or receipt['source_family'] != self.family:
                raise ValueError('Checked original quiet source required: '+stem)
            for name, digest in receipt['input_hashes'].items():
                if name in self.hashes:
                    if self.hashes[name] != digest: raise ValueError('Conflicting original quiet source: '+name)
                elif hashlib.sha256((self.directory/name).read_bytes()).hexdigest() != digest:
                    raise ValueError('Changed original quiet source: '+name)
                self.hashes[name] = digest
            self.hashes[receipt_path.name] = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
            p = self.directory / (PREFIX+'compliant_'+stem+'.json.gz')
            self.hashes[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
            raw[stem] = json.loads(gzip.decompress(p.read_bytes()))
            if raw[stem]['source_family'] != self.family: raise ValueError('Original quiet source family changed')
        self.binding = raw['current_original_true_quiet_Rp_C5_pressure_C3_functions']['true_quiet_Rp_C5_C3_bindings']
        self.tail = raw['current_original_n1_repaired_own_C3_terminal_functions']['repaired_own_C3_terminal_bindings']
        if not self.tail['definition']['mass_mean_K_C3_QV_C2_zero_after_2Rc']:
            raise ValueError('Actual n1 mass and velocity zero-tail identity required')
        quiet = raw['current_original_n1_quiet_pressure_source']['actual_quiet_source']
        if not quiet['source_scaled_Z1_exact_zero']: raise ValueError('Original quiet scalar K must be Z independent')
        self.logk = f.reserved._decode(quiet['K']['log_coefficient'])
        self.length = -60*self.outer.logmu-2-self.ctx.log(2)
        self.lengthi = -60*self.outer.logmui-2-self.iv.log(2)
        if _bounds(self.ctx, self.lengthi)[0] <= 0 or _bounds(self.ctx, self.outer.mui*self.lengthi)[1] >= self.ctx.mpf('.1'):
            raise ValueError('Original positive quiet length and small positive mu*length required')
        self.scalar = OriginalQuietScalar(self.ctx, self.iv, self.outer.mu, self.outer.mui)
        self._packets = {}
        for name in ('src/openai_ns_reconstruction/original_quiet_scalar.py', 'src/openai_ns_reconstruction/original_quiet_numeric.py'):
            self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _rows(self, point, truth): return self.compact.numeric._rows(point, truth)
    def _jet(self, rows, directed): return self.compact.numeric._jet(rows, directed)
    def _row(self, point, truth): return self.compact._summary((point, truth))

    def _flow(self, inlet, s, ctx):
        directed = ctx is self.iv
        J = self.J
        initial = {n: self._jet(v, directed).truncate(2) for n, v in inlet['generic_leading_fields'].items()}
        A = initial['E']; AA = A*A
        scalar = self.scalar.evaluate(s, ctx)
        zero = J.constant(ctx, 0, 2)
        base = dict(E=A*scalar['decay_half'], V=zero, m=initial['m']*scalar['decay1'],
                    k=initial['k']*scalar['decay_three_half'],
                    h=initial['h']*scalar['decay_three_half']+A*(scalar['decay_half']*scalar['mass1']),
                    e=initial['e']*scalar['decay1']-AA*(scalar['decay1']*s/2),
                    p=initial['p']+AA*(scalar['mass1']/2),
                    pressure=initial['pressure']+AA*(scalar['mass1']/2))
        correction = dict(E=base['E']*scalar['E_mu_divided'],
                          h=A*(scalar['decay_half']*scalar['h_mu_divided']),
                          e=-AA*(scalar['decay1']*scalar['e_mass_mu_divided']/2),
                          p=AA*(scalar['p_mass_mu_divided']/2))
        correction['pressure'] = correction['p']
        mu = self.outer.mui if directed else self.outer.mu
        full = {n: v+correction[n]*mu if n in correction else v for n, v in base.items()}
        return base, correction, full, scalar

    def evaluate(self, Z, q):
        c, iv, J = self.ctx, self.iv, self.J
        if type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1:
            raise ValueError('Interior original decimal Z required')
        if type(q) is not str or not c.isfinite(c.mpf(q)) or not 0<=Fraction(q)<=1:
            raise ValueError('Original quiet fraction q in[0,1] required')
        key = Z, q
        if key in self._packets: return self._packets[key]
        if Z not in self.inlets: self.inlets[Z] = self.compact.evaluate(Z, '2')
        inlet = self.inlets[Z]
        sp, st = self.length*c.mpf(q), self.lengthi*iv.mpf(q)
        basep, corp, fullp, scalarp = self._flow(inlet, sp, c)
        baset, cort, fullt, scalart = self._flow(inlet, st, iv)
        # The exact original relative compact histories vanish at2Rc. Full
        # leading histories and the independent datum are not set to zero.
        zero = J.constant(c, 0, 2); zeroi = J.constant(iv, 0, 2)
        Qp = self._jet(inlet['physical_similarity']['Ur_over_sqrt_R_over2'], False)*scalarp['decay1']
        Qt = self._jet(inlet['physical_similarity']['Ur_over_sqrt_R_over2'], True)*scalart['decay1']
        Pp, Pt = c.exp(self.amplitude.lp), iv.exp(self.amplitude.lp_iv)
        simp = dict(Uz=zero, Utheta_over_Pstar=basep['E'], Mz_over_R=fullp['m']*Pp, Ur_over_sqrt_R_over2=Qp)
        simt = dict(Uz=zeroi, Utheta_over_Pstar=baset['E'], Mz_over_R=fullt['m']*Pt, Ur_over_sqrt_R_over2=Qt)
        pressure = copy.deepcopy(inlet['reserved_n1_pressure'])
        pressure.update(native_compact_x=None, native_quiet_fraction=q,
                        incoming_pressure_is_actual_2Rc_value=True, fixed_quiet_source_C1_exact_zero=True)
        pressure['incoming_compact_cumulative_true_function'] = pressure.pop('cumulative_true_function')
        pressure['current_pressure_function_recipe'] = dict(
            incoming_compact_function=pressure['incoming_compact_cumulative_true_function'],
            incoming_evaluation='original compact s=log2',
            plus_quiet_increment='K^2/4*(1/(2Rc)-1/(2Rc*exp(s)))',
            quiet_coordinate='s=log(R/(2Rc))', K_source=self.binding['K'])
        pressure['limits'].update(beyond_2Rc=True, beyond_Rp=False)
        quiet_atoms = copy.deepcopy(inlet['coordinate']['global_y_log_terms'][:4]) + [self._row(c.log(110), iv.log(110)),
                        self._row(2*sum(_bounds(c, self.logk))/2-c.log(8), 2*self.logk-iv.log(8))]
        # K^2/4 * (1/(2Rc)-1/(2Rc*exp(s))) = Rc*k^2/8*(1-exp(-s)).
        pressure['P1_C0_C1_terms'].append(dict(name='quiet_fixed_leading', coefficient_C0_C1=[
            self._row(scalarp['mass1'], scalart['mass1']), self._row(c.mpf(0), iv.mpf(0))],
            log_scale_terms=quiet_atoms, actual_function_binding=dict(source=PREFIX+'compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz',
                function=self.binding['closed_pressure'], fixed_N1p='K^2/(4R^2)', physical_dR_applied_once=True,
                component='local quiet increment from2Rc in original closed-pressure formula',
                lower_radius='2Rc', primitive='Rc*k^2/8*(1-exp(-s))')))
        # Produce native pulse data from the same current fields. This is an
        # enclosure; admission to the independently selected pulse is separate.
        frames = []
        for fields, ctx, directed in ((fullp, c, False), (fullt, iv, True)):
            E, E2 = fields['E'], fields['E']*fields['E']
            p0 = self._jet(inlet['independent_P0_over_Pstar_squared'], directed).truncate(2)
            frames.append(dict(u=E, m1=fields['m']/E, m2=fields['k']/E2, X=fields['h']/E,
                               energy=fields['e']/E2, Mp=fields['p'], P0=p0, pressure=fields['pressure']))
        out = dict(Z=Z, native_quiet_fraction=q, source_family=self.family,
            selected_integer=inlet['selected_integer'], chart='actual_original_quiet_2Rc_Rp',
            quiet_log_radius=self._row(sp, st), quiet_length=self._row(self.length, self.lengthi),
            coordinate=dict(global_y_log_terms=copy.deepcopy(inlet['coordinate']['global_y_log_terms'])+[self._row(sp, st)]),
            generic_leading_fields={n:self._rows(fullp[n], fullt[n]) for n in fullp},
            generic_base_fields={n:self._rows(basep[n], baset[n]) for n in basep},
            quiet_mu_divided_fields={n:self._rows(corp[n], cort[n]) for n in corp},
            physical_similarity={n:self._rows(simp[n], simt[n]) for n in simp},
            finite_N_history_errors={n:self._rows(zero, zeroi) for n in ('m','h','k','e','p')},
            finite_N_local_errors={n:self._rows(zero, zeroi) for n in ('E','V')},
            independent_P0_over_Pstar_squared=inlet['independent_P0_over_Pstar_squared'],
            reserved_n1_pressure=pressure, native_pulse_input_frame={n:self._rows(frames[0][n], frames[1][n]) for n in frames[0]},
            normalized_mu_profiles={n:self._row(scalarp[n], scalart[n]) for n in scalarp},
            original_mu=copy.deepcopy(inlet['original_mu']),
            source_bindings=dict(flow=self.binding['flow'], K=self.binding['K'], linear_mass_theorem=self.binding['linear_mass_theorem'],
                                 quiet_pressure=self.binding['closed_pressure'], n1_zero_tail=self.tail['definition']),
            contract=dict(actual_compact_2Rc_inlet_consumed=True, full_leading_histories_not_reset=True,
                fixed_original_N_retained=True, independent_P0_added_once=True, positive_mu_terms_retained_separately=True,
                no_long_exponential_power_series=True, only_relative_finite_N_histories_zero=True,
                n1_velocity_zero_from_same_original_source_theorem=True, incoming_P1_not_zeroed=True,
                quiet_pressure_integrated_from_2Rc=True),
            limits=dict(velocity_and_pressure_values=True, leading_Z2=True, P1_C1=True,
                global_physical_gradient=False, resolved_signed_P1_centers=False, relative_pressure_cancellation_accuracy=False,
                Rp_selected_pulse_numerical_join=False, beyond_Rp=False, common_global_N=False,
                exact_heat_exterior=False, stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))
        self._packets[key] = out
        return out


class OriginalQuietVelocity:
    def __init__(self, **kwargs):
        self.numeric = OriginalQuietNumeric(**kwargs)
        f = self.numeric.compact
        self.ctx, self.iv, self.J, self.amplitude = f.ctx, f.iv, f.J, f.amplitude
        self.family = f.family
        self.hashes, self.code_hashes = dict(self.numeric.hashes), dict(self.numeric.code_hashes)
        self._summary, self._pair, self._restore, self._rotate = f._summary, f._pair, f._restore, f._rotate

    def evaluate(self, Z, q, *, log_tau='-1', scaled_time=None, azimuth='0'):
        p = self.numeric.evaluate(Z, q); c, iv = self.ctx, self.iv
        if type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)): raise ValueError('Finite decimal azimuth required')
        ell, de, time = OriginalCoreVelocity._time(self, Z, log_tau, scaled_time)
        y = [self._pair(r) for r in p['coordinate']['global_y_log_terms']]
        def scale(v, k): return v[0]*k, v[1]*k
        def logs(a, b, P=0, radius=0, constant=None, mu=False):
            out = [self._summary(scale(v, a)) for v in ell if a]+[self._summary(scale(v, b)) for v in de if b]
            if P: out.append(self._summary(scale((self.amplitude.lp, self.amplitude.lp_iv), P)))
            if radius:
                out += [self._summary(scale(v, c.mpf(radius))) for v in y]
                out.append(self._summary((c.mpf(radius)*c.log(110), iv.mpf(radius)*iv.log(110))))
            if constant is not None: out.append(self._summary(constant))
            if mu: out.append(copy.deepcopy(p['original_mu']['log_mu']))
            return out
        def term(name, row, atoms, order=0):
            return dict(name=name, temporal_order=order, coefficient=copy.deepcopy(row), log_scale_terms=atoms)
        zero = self._summary((c.mpf(0), iv.mpf(0)))
        sim = p['physical_similarity']; mu = p['quiet_mu_divided_fields']
        radial = [term('quiet_background_radial', sim['Ur_over_sqrt_R_over2'][0],
                       logs(-1, 0, radius='.5', constant=(-c.log(2)/2, -iv.log(2)/2))),
                  term('n1_quiet_radial_zero', zero, [], 1)]
        swirl = [term('quiet_background_swirl_base', sim['Utheta_over_Pstar'][0], logs(-1, -1, P=1)),
                 term('quiet_background_swirl_mu', mu['E'][0], logs(-1, -1, P=1, mu=True)),
                 term('n1_quiet_swirl_zero', zero, [], 1)]
        axial = [term('quiet_background_axial_zero', zero, []), term('n1_quiet_axial_zero', zero, [], 1)]
        stream = [term('quiet_background_streamfunction', sim['Mz_over_R'][0], logs(1, -1, radius=1)),
                  term('n1_quiet_mass_zero', zero, [], 1)]
        pressure = [term('quiet_background_pressure_base', p['generic_base_fields']['pressure'][0], logs(-2, -2, P=2)),
                    term('quiet_background_pressure_mu', mu['pressure'][0], logs(-2, -2, P=2, mu=True))]
        for source in p['reserved_n1_pressure']['P1_C0_C1_terms']:
            pressure.append(term('n1_pressure_'+source['name'], source['coefficient_C0_C1'][0],
                                 logs(-2, 0)+copy.deepcopy(source['log_scale_terms']), 1))
        angle, anglei = c.mpf(azimuth), iv.mpf(azimuth)
        cs, sn = (c.cos(angle), iv.cos(anglei)), (c.sin(angle), iv.sin(anglei))
        radius = term('physical_radius', self._summary((c.mpf(1), iv.mpf(1))),
                      logs(1, 0, radius='.5', constant=(c.log(2)/2, iv.log(2)/2)), None)
        physical = dict(source_family=self.family, time=time, native=dict(chart=p['chart'], Z=Z, q=q,
            s_definition='q*(-60logmu-2-log2)', R_definition='2Rc*exp(s)', physical_viscosity=1),
            coordinates=dict(r=radius, x=self._rotate(radius, cs), y=self._rotate(radius, sn),
                             z=term('physical_z', self._summary((c.mpf(Z), iv.mpf(Z))), logs(1, -1), None),
                             t=dict(kind='one_minus_exp_of_log_tau', definition=time)),
            cylindrical=dict(ur=radial, utheta=swirl, uz=axial),
            Cartesian=dict(u=[self._rotate(v, cs) for v in radial]+[self._rotate(v, (-sn[0], -sn[1])) for v in swirl],
                           v=[self._rotate(v, sn) for v in radial]+[self._rotate(v, cs) for v in swirl], w=copy.deepcopy(axial)),
            streamfunction=stream, pressure=pressure, coefficient_orders_included=[0, 1], pressure_orders_included=[0, 1],
            shared_streamfunction_meridional_recovery=True, full_n1_pressure_included=True,
            n1_velocity_zero_from_source_identity=True, leading_pressure_includes_same_P0_once=True,
            positive_mu_physical_terms_separate=True, n1_pressure_time_weight='lambda^(-2)',
            physical_exponentials_materialized=False, useful_combined_pressure_relative_error_certified=False,
            global_field_available=False, genuine_n_ge_2_temporal_recursion=False, full_corrected_NS_solution=False)
        result = copy.deepcopy(p); result['physical'] = physical
        return result

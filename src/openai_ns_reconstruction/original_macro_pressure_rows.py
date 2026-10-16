"""Original macro pressure jets with actual radius and incoming mean.

Background numerical proxies, directed source errors, and selected-N log
errors remain separate. Small physical sectors also retain their epsilon
factor, so adding them to the affine background cannot erase their data.
"""
import copy
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_macro_numeric import OriginalMacroNumeric
from .original_second_micro_m15 import OriginalSecondMicroPressureRows
from .original_first_micro_pressure_rows import pressure_coefficients


class OriginalMacroPressureRows(OriginalSecondMicroPressureRows):
    def __init__(self, *, micro=None, **kwargs):
        if micro is not None and kwargs:
            raise ValueError('Supply retained macro source or constructor options')
        micro = micro if micro is not None else OriginalMacroNumeric(**kwargs)
        if not isinstance(micro, OriginalMacroNumeric):
            raise ValueError('Actual original macro source required')
        super().__init__(micro=micro)
        path = self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop = json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart']['bridge_macro']
        if (self.loop['source_family'] != self.family or self.loop['chart'] != 'bridge_macro'
                or self.loop['primitive_units']['B'] != 'original B/Pstar; E=Utheta/Pstar'):
            raise ValueError('Original bridge_macro derivative source and Pstar units required')
        A = self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if A is not None and _bounds(self.audit, micro.logN)[0] < _bounds(self.audit, A)[1]:
            raise ValueError('Selected original N must admit macro A/N')
        self._macro_inputs = {}
        name = 'src/openai_ns_reconstruction/original_macro_pressure_rows.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def geometry(self, q):
        """Combine physical epsilon/rho logs before exponentiation."""
        c, iv, o = self.ctx, self.iv, self.o
        qp = sum(_bounds(c, q))/2
        tiny = o._width_error(2*(1-q))
        logR = iv.log(100)-self.micro.Y*(1-q)+iv.mpf([0, _bounds(self.audit, tiny)[1]])
        logRp = c.log(100)-self.micro.Ypoint*(1-qp)
        length = self.micro.Y-iv.mpf([0, _bounds(self.audit, 2*o.hcap)[1]])
        theta = iv.exp(-length*q)
        if q._mpi_ == iv.mpf(0)._mpi_:
            theta = iv.mpf(1)
        lo, hi = _bounds(self.audit, 1-theta)
        mass = iv.mpf([max(self.audit.mpf(0), lo), min(self.audit.mpf(1), hi)])
        return dict(logR=logR, logR_point=logRp, length=length, theta=theta, mass=mass,
                    theta_point=c.exp(-self.micro.Ypoint*qp))

    def source_Uy(self, prepared, geometry, order=5):
        """The actual prescribed macro drive, before its one hb factor.

        Comparison modes determine the source. Actual phi/barphi has the
        exp(ell) factor; the accepted weighted-jet norm of ell is <=6 cap.
        Thus ||exp(ell)-1|| <=6 cap exp(6 cap). No radial derivative of an
        integral cap is used as a source value.
        """
        o, iv, J, c = self.o, self.iv, self.J, self.ctx
        b = o.bridge
        module = __import__(PREFIX+'compliant_actual_bridge_integrals', fromlist=['weighted_modes'])
        R = b.ctx.mpf([mp.make_mpf(v) for v in geometry['logR']._mpi_])
        R = b.ctx.exp(R)
        total = J.constant(iv, 0, order)
        norm = iv.mpf(0)
        parts = {}
        with mp.workdps(280):
            for name, power, scale in (('hydro', 1, b.ctx.mpf(0)),
                                       ('pressure', 1, b.pressure_log),
                                       ('swirl', 2, b.swirl_log)):
                modes = [prepared['R0']**j*R**(power-j) for j in range(3)]
                row = -prepared['quotient']*module.weighted_modes(prepared['modes']['drive_'+name], modes)
                term = b.source_term(row, 1, scale)
                normalized = J(iv, [o._scaled(v, term, divide_eps=True)
                                    for v in row.coefficients[:order+1]])
                parts[name] = normalized
                total += normalized
                norm += o._norm(normalized)
        nonlinear = norm*o._transfer(6*b.cap*b.ctx.exp(6*b.cap))
        truth = total+o._norm_error_jet(nonlinear, order)
        point = J(c, [sum(_bounds(c, v))/2 for v in total.coefficients])
        # The signed factored coefficient is retained above. This separate
        # finite numerical enclosure bounds hb*coefficient without choosing hb.
        physical_normalized = o._norm_error_jet(o._width_error(o._norm(truth)), order)
        return dict(point=point, truth=truth, normalized_value_truth=physical_normalized,
                    normalized_value_point=J.constant(c, 0, order), parts=parts,
                    nonlinear_coefficient_error_norm=nonlinear)

    def inputs(self, Z, q):
        key = Z, q._mpi_
        if key in self._macro_inputs:
            return self._macro_inputs[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        order = 5
        # Obtain genuine Z0..5 second-micro source jets, before the older
        # public velocity packet truncates its mean to order2.
        incoming = super().inputs(Z, iv.mpf(2))
        inlet = incoming['inlet']
        a = inlet['truth']['values']
        prepared = o._prepared[Z]
        _, Bpsi = self.micro._macro_increments(prepared, _bounds(self.audit, q)[1], order)
        geo = self.geometry(q)
        point = dict(Psi=inlet['point']['values']['Psi'])
        point['Mpsi'] = (inlet['point']['values']['Mpsi']*geo['theta_point']
                        +point['Psi']*(1-geo['theta_point']))
        truth = dict(Psi=a['Psi']+o._norm_error_jet(o._width_error(Bpsi), order))
        truth['Mpsi'] = (incoming['actual']['Mpsi']*geo['theta']+a['Psi']*geo['mass']
                        +o._norm_error_jet(self.micro._uniform_Psi_increment(Z, prepared)*geo['mass'], order))
        uy = self.source_Uy(prepared, geo, order)
        point['Psi_y'], truth['Psi_y'] = uy['normalized_value_point'], uy['normalized_value_truth']
        result = dict(point=point, truth=truth, geometry=geo, source_Uy=uy, incoming=incoming)
        self._macro_inputs[key] = result
        return result

    def background(self, Z, q):
        key = Z, q._mpi_
        if key in self._backgrounds:
            return self._backgrounds[key]
        v = self.inputs(Z, q)
        def recover(ctx, values, delta, eps, axis_j):
            J = self.J
            z = J.variable(ctx, ctx.mpf(Z), 5)
            d, L = 1-z*z, 1-z*z*delta
            u0 = 4*z+axis_j
            psi, mean, psiy = (values[k] for k in ('Psi', 'Mpsi', 'Psi_y'))
            qpsi = (psi*(2*z)-mean*(z*(1-delta))-self.o._derivative(mean)*d)/L
            my = psi-mean
            qpsiy = (psiy*(2*z)-my*(z*(1-delta))-self.o._derivative(my)*d)/L
            qaxis = (u0*(z*(1+delta))-4*d)/L
            return dict(U=u0+psi*eps, m=u0+mean*eps, Q=qaxis+qpsi*eps,
                        U_y=psiy*eps, Q_y=qpsiy*eps,
                        Psi=psi, Mpsi=mean, Qpsi=qpsi, Psi_y=psiy, Qpsi_y=qpsiy)
        o = self.o
        point = recover(self.ctx, v['point'], o.core.delta, o.core.epsilon, o.core.j)
        truth = recover(self.iv, v['truth'], o.core.delta_iv, o.core.epsilon_iv, o.core.j_iv)
        result = point, truth, v['geometry']
        self._backgrounds[key] = result
        return result

    def evaluate(self, Z, coordinate, *, selected_sc_multiple=False):
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        if (type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1 < c.mpf(Z) < 1
                or selected_sc_multiple is not False):
            raise ValueError('Interior decimal Z and ordinary macro coordinate required')
        q, _, recipe = self.micro._coordinate(coordinate, False)
        key = Z, coordinate
        if key in self._packets:
            return self._packets[key]
        local, prefix, frames = self.micro._covers(Z, q)
        # All distinct axial frames survive, including shared chart endpoints.
        cover, seen = [], set()
        for cell in [p[0] for p in prefix]+local:
            if id(cell) not in seen:
                seen.add(id(cell)); cover.append(cell)
        _, B, hist = self.micro._correction_logs(cover, prefix, False)
        lp, ll = o.core.lp_iv, o.core.ll_iv
        errors = dict(U=[None if v is None else v-ll for v in B],
                      m=[None if v is None else v+lp for v in hist['m']])
        errors['U_Z'], errors['m_Z'] = errors['U'][1:], errors['m'][1:]
        slow, fast = self.loop['slow_phase_held_log_bounds']['B'], self.loop['phase_and_first_slow_phase_log_bounds']
        slow_y = [self.read_log(slow[k]) for k in ('y1_Z0', 'y1_Z1')]
        fast_y = [self.read_log(fast[k]) for k in ('B_phi', 'B_Z_phi')]
        slow_y = [None if v is None else v+lp-self.micro.logN for v in slow_y]
        fast_y = [None if v is None else v+lp for v in fast_y]
        errors['U_y'] = [self.sum_logs(v) for v in zip(slow_y, fast_y)]
        coeff = pressure_coefficients(J, iv, Z, o.core.delta_iv, 2)
        errors['Q'] = self.combine(coeff['Q'], errors, 2)
        errors['Q_y'] = self.combine(coeff['Q_y'], errors, 1)
        points, truths, geo = self.background(Z, q)
        def ordinary(name, order):
            return [o._summary(points[name][k]*math.factorial(k), truths[name][k]*math.factorial(k))
                    for k in range(order+1)]
        fields = {name: dict(signed_background_ordinary_Z=ordinary(name, order),
                            finite_N_absolute_error_log_bounds=self.log_rows(errors[name][:order+1]),
                            finite_N_log_bound_is_not_selected_signed_value=True)
                  for name, order in (('U', 2), ('m', 3), ('Q', 2), ('U_y', 1), ('Q_y', 1))}
        normalized = {name: ordinary(name, order) for name, order in
                      (('Psi', 2), ('Mpsi', 3), ('Qpsi', 2), ('Psi_y', 1), ('Qpsi_y', 1))}
        v = self.inputs(Z, q)
        uy_source = v['source_Uy']
        hb = o._summary(sum(_bounds(c, o.logh))/2, o.logh)
        hb['source_log_atom_only'] = True
        out = dict(Z=Z, coordinate=recipe, source_family=self.family, physical_pressure_rows=fields,
            normalized_background_ordinary_Z=normalized,
            normalized_sector_units=dict(U='U0+epsilon*Psi', m='U0+epsilon*Mpsi',
                Q='Qaxis+epsilon*Qpsi', U_y='epsilon*Psi_y', Q_y='epsilon*Qpsi_y',
                derivative_coordinate='ordinary fixed-native-q Z derivatives'),
            leading_U_y_factored=dict(coefficient_ordinary_Z=[o._summary(
                uy_source['point'][k]*math.factorial(k), uy_source['truth'][k]*math.factorial(k)) for k in (0, 1)],
                log_scale_terms=[o._summary(-o.core.logLambda, -ll), hb],
                definition='U_y=epsilon*hb*(original normalized macro drive plus bounded exp(ell) source error)',
                raw_numerical_zero_proxy_does_not_zero_factored_source=True,
                actual_source='-hb*(Phi/barPhi)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl)',
                comparison_modes_not_replaced_by_core_radial_derivatives=True),
            R_Q_y_endpoint=dict(signed_background_ordinary_Z=normalized['Qpsi_y'],
                log_physical_factor=o._summary(geo['logR_point']-o.core.logLambda, geo['logR']-ll),
                finite_N_coefficient_error_log_bounds=self.log_rows(errors['Q_y']),
                finite_N_log_physical_factor=o._summary(geo['logR_point'], geo['logR']),
                definition='background R*epsilon*Qpsi_y; finite-N R*deltaQ_y; these factors differ by epsilon'),
            geometry=dict(log_physical_radius=o._summary(geo['logR_point'], geo['logR']),
                actual_log_length=o._summary(self.micro.Ypoint, geo['length']), original_log_hb=hb,
                R_definition='4epsilon*exp(2hb+(Y-2hb)*q)', L_definition='Y-2hb',
                native_measure='dy=(Y-2hb)dq; R*dR=R^2*(Y-2hb)dq',
                R100_endpoint_exact=q._mpi_ == iv.mpf(1)._mpi_,
                defining_hb_not_selected_from_cap=True, epsilon_rho_logs_combined_first=True),
            selected_integer=copy.deepcopy(self.micro.integer), original_sc_tuple=self.t.sc_tuple,
            finite_N_axial_radial_error_components=dict(slow_B_y_over_N=self.log_rows(slow_y),
                fast_B_phi=self.log_rows(fast_y), fast_Z_row_source='B_Z_phi',
                formula='Pstar*(B_y/N+B_phi); mixed Pstar*(B_yZ/N+B_Z_phi)', extra_hb_factor=False),
            source_cover=dict(exact_Z_cells=frames, local_source_cells=len(local), prefix_source_cells=len(prefix),
                uniform_prefix_source_cells=len(cover), generic_derivative_chart='bridge_macro',
                cumulative_history_charts=['first_micro', 'second_micro', 'frozen_macro'],
                macro_prefix_measure_bound_only='Y*dq; defining physical measure is (Y-2hb)dq',
                phase='frac(N*(2hb+(Y-2hb)*q-hb*original_sc/2))',
                incoming_second_mean_not_zeroed=True, all_axial_frames_retained=True),
            contract=dict(actual_selected_N_not_candidate_257=True, original_Pstar_once=True,
                genuine_mean_Z3_before_error_serialization=True, mean_Z3_error_kept_logarithmic=True,
                Q_ZZ_uses_mean_ZZZ_error=True, no_missing_rows_zero_padded=True,
                full_second_micro_mean_inlet_retained=True, actual_large_Y_mean_relaxation=True,
                mean_y_equals_Psi_minus_Mpsi=True, genuine_macro_modal_Uy_source=True,
                nonlinear_exp_ell_source_error_retained=True, slow_B_y_over_N=True,
                fast_B_phi_not_divided_by_N=True, fast_mixed_B_Z_phi_retained=True,
                no_extra_hb_in_fast_correction=True, actual_macro_geometry_and_phase=True,
                finite_R_Qy_and_ordinary_Z1_retained=True, normalized_epsilon_sectors_retained=True),
            limits=dict(actual_macro_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors=True,
                resolved_finite_N_signed_centers=False, finite_total_relative_accuracy=False,
                macro_signed_m15_integrated=False, independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False, common_global_N=False, exact_heat_exterior=False,
                stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))
        self._packets[key] = out
        return out

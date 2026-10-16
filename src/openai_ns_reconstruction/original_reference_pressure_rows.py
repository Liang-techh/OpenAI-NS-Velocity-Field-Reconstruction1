"""Genuine pressure sectors through reference, axial restore and Rm buffer."""
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_long_pressure_rows import OriginalLongPressureRows
from .original_macro_pressure_rows import OriginalMacroPressureRows
from .original_r110_numeric import OriginalSecondSwitchNumeric
from .original_reference_numeric import OriginalReferenceNumeric


class OriginalReferencePressureRows(OriginalLongPressureRows):
    OWNERS = dict(reference='inner_reference', restoration='axial_restore', postrestore='restore_buffer')

    def __init__(self, *, reference=None, chart='reference', long_rows=None, previous_rows=None, **kwargs):
        if reference is not None and kwargs:
            raise ValueError('Supply original reference source or constructor options')
        if chart not in self.OWNERS: raise ValueError('Original reference-family chart required')
        self.reference = reference if reference is not None else OriginalReferenceNumeric(**kwargs)
        self.chart, self.derivative_chart = chart, self.OWNERS[chart]
        self.long = self.reference.previous
        view = object.__new__(OriginalSecondSwitchNumeric)
        view.__dict__ = dict(self.reference.owner.__dict__)
        view._micro_cache = {}; view._upgraded = set(); view.chart = chart
        view._covers = lambda Z, q: self.reference._covers(Z, self.chart, q)
        OriginalMacroPressureRows.__init__(self, micro=view)
        self._locations = {}; view._coordinate = self._coordinate
        self.long_rows = long_rows if long_rows is not None else OriginalLongPressureRows(long=self.long)
        self.previous_rows = previous_rows if previous_rows is not None else self.long_rows
        if self.long_rows.ctx is not self.ctx or self.previous_rows.iv is not self.iv:
            raise ValueError('Shared exact source/Taylor contexts required')
        if chart != 'reference' and previous_rows is None:
            raise ValueError('Explicit preceding genuine pressure source required')
        self._reference_inputs = {}; self._next = None; self._mean_kernels = {}
        self.hashes.update(self.reference.hashes); self.code_hashes.update(self.reference.code_hashes)
        self.hashes.update(self.long_rows.hashes); self.code_hashes.update(self.long_rows.code_hashes)
        path = self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop = json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart'][self.derivative_chart]
        if self.loop['source_family'] != self.family or self.loop['primitive_units']['B'] != 'original B/Pstar; E=Utheta/Pstar':
            raise ValueError('Original chart-specific derivative owner and Pstar units required')
        bound = self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if bound is not None and _bounds(self.audit, self.micro.logN)[0] < _bounds(self.audit, bound)[1]:
            raise ValueError('Original selected N must admit local A/N')
        for suffix in ('long_m15_numeric', 'reference_velocity'):
            path = self.directory/(PREFIX+'compliant_current_original_'+suffix+'.json')
            raw = path.read_bytes(); data = json.loads(raw)
            if not data['checks']['passed'] or data['source_family'] != self.family:
                raise ValueError('Accepted original reference pressure prerequisites required')
            for base, group, target in ((self.directory, data['input_hashes'], self.hashes),
                                         (ROOT, data['implementation_hashes'], self.code_hashes)):
                for name, digest in group.items():
                    if target.get(name, digest) != digest or hashlib.sha256((base/name).read_bytes()).hexdigest() != digest:
                        raise ValueError('Changed original reference pressure source: '+name)
                    target[name] = digest
            self.hashes[path.name] = hashlib.sha256(raw).hexdigest()
        name = 'src/openai_ns_reconstruction/original_reference_pressure_rows.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):
        if self.chart == 'postrestore': raise ValueError('Original pressure source already at Rm')
        if self._next is None:
            chart = 'restoration' if self.chart == 'reference' else 'postrestore'
            self._next = OriginalReferencePressureRows(reference=self.reference, chart=chart,
                long_rows=self.long_rows, previous_rows=self)
        return self._next

    def at_reference_gap(self, Z, gap):
        return self.evaluate(Z, 'g:'+gap)

    def _coordinate(self, coordinate, selected=False):
        c, iv = self.ctx, self.iv
        if selected or type(coordinate) is not str:
            raise ValueError('Original decimal chart coordinate required')
        native = coordinate.startswith('g:')
        raw = coordinate[2:] if native else coordinate
        if not c.isfinite(c.mpf(raw)) or (native and self.chart != 'reference'):
            raise ValueError('Finite native coordinate; g:gap belongs to reference only')
        qp = c.mpf(raw); q = iv.mpf(raw)
        if native:
            gp, gap = qp, q; q = gap/self.reference.length
        else:
            gp, gap = self.reference.lengthp*qp, self.reference.length*q
        lo, hi = _bounds(self.audit, q)
        if lo < 0 or hi > 1: raise ValueError('Original chart coordinate in[0,1] required')
        recipe = dict(kind='exact_reference_gap' if native else 'exact_decimal_'+self.chart, value=raw)
        self._locations[q._mpi_] = (gp, gap, recipe) if self.chart == 'reference' else (qp, q, recipe)
        return q, False, recipe

    def geometry(self, q):
        c, iv, ref = self.ctx, self.iv, self.reference
        post = self.chart == 'postrestore'
        if self.chart == 'reference':
            default = (ref.lengthp*sum(_bounds(c, q))/2, ref.length*q, None)
            gp, gap, recipe = self._locations.get(q._mpi_, default)
            qp = gp/ref.lengthp
            if recipe and recipe['kind'] == 'exact_reference_gap':
                terms = [(ref.Tp, ref.T), (gp, gap)]
            else:
                terms = [(ref.Tp*(1-qp), ref.T*(1-q)), (10*ref.logCp*qp, 10*ref.logC*q),
                         (10*ref.logPp*qp, 10*ref.logP*q), (-8*qp, -8*q)]
            length = ref.length
        else:
            gp, gap, recipe = self._locations.get(q._mpi_, (sum(_bounds(c, q))/2, q, None))
            offset = -7 if post else -8
            terms = [(10*ref.logCp, 10*ref.logC), (10*ref.logPp, 10*ref.logP), (offset+gp, offset+gap)]
            length = iv.mpf(1)
        proof = []
        decay = ref.previous.kernels.positive_decay(iv, 1, gap, proof)
        lo, hi = _bounds(self.audit, 1-decay)
        return dict(logR_point=c.log(110)+sum((a for a,b in terms), c.mpf(0)),
            logR=iv.log(110)+sum((b for a,b in terms), iv.mpf(0)), length=length,
            log_offset_terms=terms, local_point=gp, local=gap, mass=iv.mpf([max(0,lo), min(1,hi)]),
            gap_decay_proofs=proof, recipe=recipe)

    def inputs(self, Z, q):
        key = Z, q._mpi_
        if key in self._reference_inputs: return self._reference_inputs[key]
        c, iv, J, o, ref = self.ctx, self.iv, self.J, self.o, self.reference
        incoming = self.previous_rows.inputs(Z, iv.mpf(1)); geo = self.geometry(q)
        gp, gap = geo['local_point'], geo['local']; proof = []
        point, truth = {}, {}
        if self.chart == 'reference':
            pair = ref._decay(incoming['point']['Mpsi']-incoming['point']['Psi'],
                incoming['truth']['Mpsi']-incoming['truth']['Psi'], 1, gp, gap, proof, 'pressure_reference_mean_memory')
            for ctx, values, initial, delta in ((c, point, incoming['point'], pair[0]),
                                                 (iv, truth, incoming['truth'], pair[1])):
                z = J.variable(ctx, ctx.mpf(Z), 5)
                values.update(Psi=initial['Psi'], Mpsi=initial['Mpsi'] if not gp else initial['Psi']+delta,
                    Psi_y=J.constant(ctx, 0, 5), U0=4*z+(o.core.j if ctx is c else o.core.j_iv),
                    M0=4*z+(o.core.j if ctx is c else o.core.j_iv), U0_y=J.constant(ctx, 0, 5))
        elif self.chart == 'restoration':
            ap, at = ref._alpha(gp, gap)
            if gap._mpi_ not in self._mean_kernels:
                # The point proxy may round gp; the original interval source
                # must use the exact query enclosure, including long decimals.
                proxy_coordinate = str(gp)
                kp, kt = ref.kernels(proxy_coordinate)['mean']
                if iv.mpf(proxy_coordinate)._mpi_ != gap._mpi_:
                    kt = ref.source.restoration_kernels(iv, gap, cells=128)['mean']
                self._mean_kernels[gap._mpi_] = kp, kt
            kp, kt = self._mean_kernels[gap._mpi_]
            if gp == 0 or gp == 1: sp = c.mpf(0)
            else: sp = ap*(1-ap)*(2/(1-gp)**3+2/gp**3)
            st = ref.sigma_jets(iv, gap)[1]
            for ctx, values, initial, alpha, sigma1, kernel, local, j in (
                (c, point, incoming['point'], ap, sp, kp, gp, o.core.j),
                (iv, truth, incoming['truth'], at, st, kt, gap, o.core.j_iv)):
                z = J.variable(ctx, ctx.mpf(Z), 5); decay = ctx.exp(-local)
                values.update(Psi=initial['Psi']*alpha, Mpsi=initial['Mpsi']*decay+initial['Psi']*kernel,
                    Psi_y=initial['Psi']*(-sigma1), U0=4*z+j*alpha, M0=4*z+j*(decay+kernel),
                    U0_y=J.constant(ctx, -j*sigma1, 5))
        else:
            for ctx, values, initial, local in ((c, point, incoming['point'], gp), (iv, truth, incoming['truth'], gap)):
                z = J.variable(ctx, ctx.mpf(Z), 5); decay = ctx.exp(-local)
                values.update(Psi=J.constant(ctx, 0, 5),
                    Mpsi=initial['Mpsi'] if not gp else initial['Mpsi']*decay,
                    Psi_y=J.constant(ctx, 0, 5), U0=4*z,
                    M0=initial['M0'] if not gp else 4*z+(initial['M0']-4*z)*decay,
                    U0_y=J.constant(ctx, 0, 5))
        source = dict(point=point['Psi_y'], truth=truth['Psi_y'], parts={},
            normalized_value_point=point['Psi_y'], normalized_value_truth=truth['Psi_y'], nonlinear_coefficient_error_norm=iv.mpf(0))
        result = dict(point=point, truth=truth, geometry=geo, source_Uy=source, incoming=incoming, mean_memory_proofs=proof)
        self._reference_inputs[key] = result
        return result

    def recover_values(self, Z, ctx, values):
        o, J = self.o, self.J
        de, eps = (o.core.delta, o.core.epsilon) if ctx is self.ctx else (o.core.delta_iv, o.core.epsilon_iv)
        z = J.variable(ctx, ctx.mpf(Z), 5); A, D = 1-z*z*de, 1-z*z
        def recover(u, m, uy):
            q = (u*(2*z)-m*(z*(1-de))-o._derivative(m)*D)/A
            qy = (uy*(2*z)-(u-m)*(z*(1-de))-o._derivative(u-m)*D)/A
            return q, qy
        q0, qy0 = recover(values['U0'], values['M0'], values['U0_y'])
        q1, qy1 = recover(values['Psi'], values['Mpsi'], values['Psi_y'])
        return dict(U=values['U0']+values['Psi']*eps, m=values['M0']+values['Mpsi']*eps,
            Q=q0+q1*eps, U_y=values['U0_y']+values['Psi_y']*eps, Q_y=qy0+qy1*eps,
            Psi=values['Psi'], Mpsi=values['Mpsi'], Psi_y=values['Psi_y'], Qpsi=q1, Qpsi_y=qy1,
            U0=values['U0'], M0=values['M0'], Q0=q0, U0_y=values['U0_y'], Q0_y=qy0)

    def background(self, Z, q):
        key = Z, q._mpi_
        if key not in self._backgrounds:
            v = self.inputs(Z, q)
            self._backgrounds[key] = self.recover_values(Z,self.ctx,v['point']), self.recover_values(Z,self.iv,v['truth']), v['geometry']
        return self._backgrounds[key]

    def uniform_background(self, Z, q):
        iv, J = self.iv, self.J; v = self.inputs(Z, q); initial = v['incoming']['truth']
        z = J.variable(iv, iv.mpf(Z), 5); j = self.o.core.j_iv
        if self.chart == 'reference':
            values = dict(v['truth']); values['Mpsi'] = initial['Mpsi']+(initial['Psi']-initial['Mpsi'])*iv.mpf([0,1])
        elif self.chart == 'restoration':
            e = iv.mpf([_bounds(self.audit, iv.exp(-q))[0],1]); a = iv.mpf([0,1]); k = iv.mpf([0,_bounds(self.audit,1-iv.exp(-q))[1]])
            values = dict(Psi=initial['Psi']*a, Mpsi=initial['Mpsi']*e+initial['Psi']*k,
                Psi_y=v['truth']['Psi_y'], U0=4*z+j*a, M0=4*z+j*(e+k), U0_y=v['truth']['U0_y'])
        else:
            values = dict(v['truth']); e = iv.mpf([_bounds(self.audit,iv.exp(-q))[0],1])
            values['M0'] = 4*z+(initial['M0']-4*z)*e; values['Mpsi'] = initial['Mpsi']*e
        return self.recover_values(Z, iv, values)

    def evaluate(self, Z, coordinate, *, selected_sc_multiple=False):
        out = OriginalMacroPressureRows.evaluate(self,Z,coordinate,selected_sc_multiple=selected_sc_multiple)
        q,_,_ = self._coordinate(coordinate,False); v = self.inputs(Z,q); geo=v['geometry']; o=self.o
        out['leading_U_y_factored'] = dict(definition='U_y=U0_y+epsilon*Psi_y; j baseline is separate',
            baseline_ordinary_Z=o._jet_summary(v['point']['U0_y'].truncate(1),v['truth']['U0_y'].truncate(1)),
            epsilon_coefficient_ordinary_Z=o._jet_summary(v['point']['Psi_y'].truncate(1),v['truth']['Psi_y'].truncate(1)),
            epsilon_log_atom=o._summary(-o.core.logLambda,-o.core.ll_iv), zero_only_on_constant_U_source=self.chart!='restoration')
        out['baseline_pressure_rows']={name:o._jet_summary(self.background(Z,q)[0][name].truncate(2),self.background(Z,q)[1][name].truncate(2))
            for name in ('U0','M0','Q0','U0_y','Q0_y')}
        out['baseline_pressure_row_convention']='Taylor coefficients D_Z^k/k!; physical_pressure_rows are ordinary derivatives'
        out['normalized_sector_units'].update(m='M0+epsilon*Mpsi', Q='Q0+epsilon*Qpsi',
            U_y='U0_y+epsilon*Psi_y', Q_y='Q0_y+epsilon*Qpsi_y')
        endpoint = out['R_Q_y_endpoint']
        epsilon_rows = endpoint.pop('signed_background_ordinary_Z', out['normalized_background_ordinary_Z']['Qpsi_y'])
        endpoint.update(baseline_ordinary_Z=o._jet_summary(self.background(Z,q)[0]['Q0_y'].truncate(1),self.background(Z,q)[1]['Q0_y'].truncate(1)),
            baseline_log_physical_factor=o._summary(geo['logR_point'],geo['logR']),
            epsilon_coefficient_ordinary_Z=epsilon_rows,
            definition='R*Q0_y + R*epsilon*Qpsi_y; finite-N R*deltaQy has no extra epsilon')
        out['geometry']=dict(log_physical_radius=out['geometry']['log_physical_radius'],
            atomic_log_offset_terms=[o._summary(a,b) for a,b in geo['log_offset_terms']],
            actual_native_length=o._summary(self.reference.lengthp if self.chart=='reference' else self.ctx.mpf(1),geo['length']),
            reference_T_singleton_exact_mpf_tuple=self.long.T_tuple,
            local_source_coordinate=o._summary(geo['local_point'],geo['local']),
            native_measure='Lref dq' if self.chart=='reference' else 'dt or dp once',
            unit_offsets_not_merged_into_huge_radius=True, physical_radius_materialized=False,
            exact_Rm_endpoint=self.chart=='postrestore' and q._mpi_==self.iv.mpf(1)._mpi_,
            mean_memory_proofs=self.long._serializable_proof(v['mean_memory_proofs']))
        charts=['first_micro','second_micro','frozen_macro','first_switch','second_switch','post_power','long_reshape']
        charts+=list(self.reference.charts[:self.reference.charts.index(self.chart)+1])
        out['source_cover'].update(generic_derivative_chart=self.derivative_chart,cumulative_history_charts=charts,
            generic_source_domain=['-7','-6'] if self.chart=='postrestore' else ['0','1'],
            wrapper_coordinate_domain=['0','1'],phase='frac(N*(log(rho/4)-hbB*original_sc/2))')
        out['source_cover'].pop('macro_prefix_measure_bound_only',None)
        for name in ('genuine_macro_modal_Uy_source','actual_macro_geometry_and_phase','actual_large_Y_mean_relaxation'):
            out['contract'].pop(name,None)
        out['contract'].update(genuine_long_fivejet_pressure_inlet=True,actual_j_baseline_restoration_retained=True,
            full_mean_ODE_through_Z5=True,correct_local_generic_derivative_owner=True,
            both_baseline_and_epsilon_Uy_source_retained=True,atomic_reference_offsets_retained=True,
            all_upstream_and_current_prefix_measures_retained=True,mean_Z3_uncertainty_not_zero_padded=True)
        out['limits'].pop('actual_macro_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors',None)
        out['limits'].pop('macro_signed_m15_integrated',None)
        out['limits'].update(actual_reference_restore_post_pressure_rows=True,reference_family_m15_integrated=False)
        return out

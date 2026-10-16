"""Complete first-micro pressure transport with uniform prefix errors.

Original quarter-slice implementations and their receipts remain unchanged.
This extension covers every original first-micro source cell and handles
the exact sigma=1 endpoint without forming an infinite log-odds atom.
"""
import copy
import hashlib
import math

import mpmath as mp

from .original_core_numeric import ROOT,_bounds
from .original_first_micro_pressure_rows import OriginalFirstMicroPressureRows,pressure_coefficients
from .original_first_micro_m15 import OriginalFirstMicroM15
from .original_n1_inner_m15_numeric import normalized_pressure_source


class OriginalFirstMicroFullPressureRows(OriginalFirstMicroPressureRows):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        expected=[(['sc_half'],[1,4]),([1,4],[1,2]),([1,2],[3,4]),([3,4],[1,1])]
        for frame in self.micro.ranges:
            chart=next(v for v in frame['actual_17_chart_C3_transports'] if v['chart']=='first_micro')
            if [(v['exact_left'],v['exact_right']) for v in chart['cells']]!=expected:
                raise ValueError('Complete original first-micro partition required')
        name='src/openai_ns_reconstruction/original_first_micro_full_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    @staticmethod
    def prefix_cover(local,prefix):
        """Keep all axial-frame cells; deduplicate only identical objects."""
        cells=[];seen=set()
        for cell in [v for v,width in prefix]+list(local):
            if id(cell) not in seen:seen.add(id(cell));cells.append(cell)
        return cells

    def evaluate(self,Z,coordinate,*,selected_sc_multiple=False):
        c,iv=self.ctx,self.iv
        if (type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1
                or type(selected_sc_multiple) is not bool):
            raise ValueError('Interior decimal Z and boolean original-coordinate selector required')
        s,flat,recipe=self.micro._coordinate(coordinate,selected_sc_multiple)
        key=Z,coordinate,selected_sc_multiple
        if key in self._packets:return self._packets[key]
        local,prefix,frames=self.micro._covers(Z,s)
        cover=self.prefix_cover(local,prefix)
        _,B,hist=self.micro._correction_logs(cover,prefix,flat)
        lp,ll=self.o.core.lp_iv,self.o.core.ll_iv
        errors=dict(U=[None if q is None else q-ll for q in B],
                    m=[None if q is None else q+lp for q in hist['m']])
        errors['U_Z'],errors['m_Z']=errors['U'][1:],errors['m'][1:]
        slow=self.loop['slow_phase_held_log_bounds']['B'];fast=self.loop['phase_and_first_slow_phase_log_bounds']
        slow_y=[self.read_log(slow[k]) for k in ('y1_Z0','y1_Z1')]
        phase_y=[self.read_log(fast[k]) for k in ('B_phi','B_Z_phi')]
        slow_y=[None if flat or v is None else v+lp-self.micro.logN for v in slow_y]
        phase_y=[None if flat or v is None else v+lp for v in phase_y]
        errors['U_y']=[self.sum_logs(v) for v in zip(slow_y,phase_y)]
        coeff=pressure_coefficients(self.J,iv,Z,self.o.core.delta_iv,2)
        errors['Q']=self.combine(coeff['Q'],errors,2);errors['Q_y']=self.combine(coeff['Q_y'],errors,1)
        points,truths,rho=self.background(Z,s)
        fields={}
        for name,order in (('U',2),('m',3),('Q',2),('U_y',1),('Q_y',1)):
            rows=[self.o._summary(points[name][k]*math.factorial(k),truths[name][k]*math.factorial(k))
                  for k in range(order+1)]
            fields[name]=dict(signed_background_ordinary_Z=rows,
                finite_N_absolute_error_log_bounds=self.log_rows(errors[name][:order+1]),
                total_error='background absolute_error_upper plus exp(finite-N log bound)',
                finite_N_log_bound_is_not_selected_signed_value=True)
        endpoint=dict(signed_background_ordinary_Z=[self.o._summary(points['Q_y'][k]*4,truths['Q_y'][k]*rho[0]) for k in (0,1)],
            log_physical_factor=self.o._summary(-self.o.core.logLambda,-ll),
            finite_N_coefficient_error_log_bounds=self.log_rows([None if v is None else v+iv.log(rho[0]) for v in errors['Q_y']]),
            definition='R*Q_y=epsilon*rho*Q_y; epsilon multiplies value and errors once')
        source_cover=dict(exact_Z_cells=frames,local_source_cells=len(local),prefix_source_cells=len(prefix),
            uniform_prefix_source_cells=len(cover),exact_source_boundaries=[(v['exact_left'],v['exact_right']) for v in cover],
            exact_flat_collar=flat,generic_derivative_chart='bridge_first',
            actual_radius='R=4epsilon*exp(hb*s)',native_measure='dy=hb*ds',phase='frac(N*hb*(s-sc/2))',
            original_sc_tuple=self.t.sc_tuple,correction_errors_cover_full_original_prefix=True,
            all_axial_source_frames_retained=True,local_caps_not_misused_as_whole_prefix=True)
        out=dict(Z=Z,coordinate=recipe,source_family=self.family,physical_pressure_rows=fields,
            R_Q_y_endpoint=endpoint,selected_integer=copy.deepcopy(self.micro.integer),source_cover=source_cover,
            finite_N_axial_radial_error_components=dict(slow_B_y_over_N=self.log_rows(slow_y),fast_B_phi=self.log_rows(phase_y),
                fast_Z_row_source='B_Z_phi',formula='Pstar*(B_y/N+B_phi); mixed Pstar*(B_yZ/N+B_Z_phi)',extra_hb_factor=False),
            contract=dict(actual_selected_N_not_candidate_257=True,original_Pstar_applied_once=True,
                slow_B_y_divided_by_N=True,fast_B_phi_not_divided_by_N=True,fast_mixed_B_Z_phi_retained=True,
                ordinary_Z_factorial_once=True,mean_correction_history_from_rminus_retained=True,
                mean_y_equals_U_minus_m=True,Q_ZZ_includes_mean_ZZZ_error=True,
                original_first_micro_phase_and_absolute_radius_preserved=True,
                all_original_source_cells_before_query_in_uniform_error_cover=True),
            limits=dict(actual_full_first_micro_pressure_rows_with_factored_errors=True,
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                second_micro_pressure_integrated=False,independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out
        return out


class OriginalFirstMicroFullM15(OriginalFirstMicroM15):
    def __init__(self,*,cells=8,rows=None,**kwargs):
        if rows is not None and kwargs:raise ValueError('Supply retained full rows or constructor options')
        rows=rows if rows is not None else OriginalFirstMicroFullPressureRows(**kwargs)
        if not isinstance(rows,OriginalFirstMicroFullPressureRows):
            raise ValueError('Actual full-prefix first-micro rows required')
        super().__init__(cells=cells,rows=rows)
        name='src/openai_ns_reconstruction/original_first_micro_full_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def sigma_atom(self,s):
        # On the upper half retain sigma directly as the coefficient. At
        # s=1 it is exactly1; no division by (1-s)^2 is performed.
        if _bounds(self.audit,s)[0]>=self.audit.mpf('.5'):
            with mp.workdps(280):sigma=self.t.sigma_jets(self.iv,s)[0]
            atom=self.o._summary(self.ctx.mpf(0),self.iv.mpf(0))
            atom.update(original_nonzero_sigma_log_atom=True,upper_half_sigma_direct_coefficient=True)
            point=sum(_bounds(self.ctx,sigma))/2
            return atom,point,sigma
        return super().sigma_atom(s)

    def regular_integral(self,Z,s):
        key=Z,s._mpi_
        if key in self._integrals:return self._integrals[key]
        c,iv,J=self.ctx,self.iv,self.J;p=self.prepare(Z,s)
        value=J.constant(c,0,1);truth=J.constant(iv,0,1)
        ends=[self.s0]+[iv.mpf(v) for v in ('.25','.5','.75') if self.audit.mpf(v)<_bounds(self.audit,s)[0]]+[s]
        for lo,hi in zip(ends,ends[1:]):
            width=(hi-lo)/self.cells
            for k in range(self.cells):
                left,right=lo+width*k,lo+width*(k+1)
                box=iv.mpf([_bounds(self.audit,left)[0],_bounds(self.audit,right)[1]])
                fp,ft=self.frame(p,Z,box,point=True),self.frame(p,Z,box,point=False)
                dp=-(normalized_pressure_source(fp).rows[0]*fp['rho']*fp['rho']).truncate(1)
                dt=-(normalized_pressure_source(ft).rows[0]*ft['rho']*ft['rho']).truncate(1)
                value+=dp*(sum(_bounds(c,width))/2);truth+=dt*width
        self._integrals[key]=value,truth
        return value,truth

    def evaluate(self,Z,coordinate,*,selected_sc_multiple=False):
        out=super().evaluate(Z,coordinate,selected_sc_multiple=selected_sc_multiple)
        source=self.rows.evaluate(Z,coordinate,selected_sc_multiple=selected_sc_multiple)
        out['finite_N_full_prefix_cover']=copy.deepcopy(source['source_cover'])
        out['original_first_micro_partition']=['sc/2','1/4','1/2','3/4','1']
        out['contract'].update(all_traversed_finite_N_source_cells_retained=True,
            regular_cells_split_at_original_quarters=True,upper_sigma_endpoint_avoids_log_odds_division_by_zero=True)
        out['limits']['complete_first_micro_m15_with_factored_finite_N_errors']=True
        out['limits']['second_micro_m15_integrated']=False
        return out

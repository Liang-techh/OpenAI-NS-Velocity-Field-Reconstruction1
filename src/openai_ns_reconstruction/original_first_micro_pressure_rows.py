"""First-micro physical pressure rows with selected-N factored errors.

Signed leading approximations and directed errors are separate from the
unknown signed finite-N correction. Fast phase derivatives are retained;
large error exponentials, the selected integer and its angle stay formal.
"""
import copy
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_micro_numeric import OriginalFirstMicroNumeric
from .original_n1_cutoff_numeric import stable_sigma_point
from .original_n1_inner_m15_numeric import leading_regular_frame


def pressure_coefficients(J,ctx,Z,delta,order):
    """Exact physical leading Q and Q_y recovery, as axial Taylor jets."""
    z=J.variable(ctx,ctx.mpf(Z),order)
    L=1-z*z*delta;d=1-z*z
    return dict(Q={'U':2*z/L,'m':-z*(1-delta)/L,'m_Z':-d/L},
        Q_y={'U_y':2*z/L,'U':-z*(1-delta)/L,'m':z*(1-delta)/L,
             'U_Z':-d/L,'m_Z':d/L})


class OriginalFirstMicroPressureRows:
    def __init__(self,*,micro=None,**kwargs):
        if micro is not None and kwargs:raise ValueError('Supply retained micro or constructor options')
        self.micro=micro if micro is not None else OriginalFirstMicroNumeric(**kwargs)
        self.t,self.o=self.micro.transport,self.micro.leading
        self.ctx,self.iv,self.J=self.t.ctx,self.t.iv,self.t.Jet
        self.audit=mp.mp.clone();self.audit.dps=self.iv.dps+30
        self.directory,self.family=self.o.core.directory,self.micro.family
        self.hashes,self.code_hashes=dict(self.micro.hashes),dict(self.micro.code_hashes)
        name=PREFIX+'compliant_current_generic_shear_loop_jet_bounds'
        path=self.directory/(name+'.json');receipt=self.directory/(name+'_check.json')
        data=json.loads(path.read_bytes());checked=json.loads(receipt.read_bytes())
        gate='current_actual_phase_held_inverse_and_loop_y2_yZ_log_bounds_certified'
        if (data['source_family']!=self.family or checked['source_family']!=self.family
                or not checked['all_passed'] or not checked[gate] or not data[gate]
                or not checked['B_bounds_in_original_common_Pstar_unit']):
            raise ValueError('Same original physical slow/fast derivative bounds required')
        for key,digest in checked['input_hashes'].items():
            if (self.hashes.get(key,digest)!=digest
                    or hashlib.sha256((self.directory/key).read_bytes()).hexdigest()!=digest):
                raise ValueError('Changed original loop derivative input: '+key)
            self.hashes[key]=digest
        for p in (path,receipt):self.hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        self.loop=data['current_actual_loop_jet_log_bounds_by_chart']['bridge_first']
        if (self.loop['source_family']!=self.family or self.loop['chart']!='bridge_first'
                or self.loop['primitive_units']['B']!='original B/Pstar; E=Utheta/Pstar'
                or not self.loop['primitive_units']['Pstar_is_global_constant_under_slow_derivatives']):
            raise ValueError('Actual first-micro bridge_first and B/Pstar units required')
        # Bind the actual selected integer; no candidate_N=257 projection.
        A=self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if A is not None and _bounds(self.ctx,self.micro.logN)[0]<_bounds(self.ctx,A)[1]:
            raise ValueError('Actual selected N must satisfy |A/N|<=1')
        for frame in self.micro.ranges:
            chart=next(v for v in frame['actual_17_chart_C3_transports'] if v['chart']=='first_micro')
            if not all(v['exact_zero'] for v in chart['actual_N_scaled_incoming_C3']['m']):
                raise ValueError('Source-defined zero correction inlet at rminus required')
        self._backgrounds,self._packets={},{}
        for name in ('src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py',
                     'src/openai_ns_reconstruction/original_first_micro_pressure_rows.py'):
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def read_log(self,row):
        return None if row['exact_zero'] else self.micro._decode(row['log_absolute_upper'])

    def sum_logs(self,values):
        """Directed positive sum bound; never exponentiate a huge log."""
        values=[v for v in values if v is not None]
        if not values:return None
        high=max(_bounds(self.ctx,v)[1] for v in values)
        return self.iv.mpf(high)+self.iv.log(len(values))

    def scaled_log(self,value,coefficient):
        if value is None:return None
        # Source interval endpoints may carry more bits than the point
        # context. Do not round their absolute upper endpoint inward.
        cap=max(abs(x) for x in _bounds(self.audit,coefficient))
        return None if not cap else value+self.iv.log(self.iv.mpf(cap))

    def log_rows(self,values):return [self.micro._source_summary(v) for v in values]

    def combine(self,coefficients,errors,order):
        out=[]
        for k in range(order+1):
            pieces=[]
            for name,coef in coefficients.items():
                for j in range(k+1):
                    ordinary=coef[j]*(math.factorial(j)*math.comb(k,j))
                    pieces.append(self.scaled_log(errors[name][k-j],ordinary))
            out.append(self.sum_logs(pieces))
        return out

    def background(self,Z,s):
        key=Z,s._mpi_
        if key in self._backgrounds:return self._backgrounds[key]
        c,iv,J=self.ctx,self.iv,self.J;upper=_bounds(c,s)[1]
        v=self.t._inputs(Z,upper,order=5)
        with mp.workdps(280):sigma=self.t.sigma_jets(iv,s)[0]
        inlet=v['inlet']['point'];z=J.variable(c,c.mpf(Z),5)
        point=leading_regular_frame(inlet['values'],inlet['radial'],z,J.constant(c,4,5),
            J.constant(c,2,5),self.o.core.delta,self.o.core.epsilon,c.exp(2*self.o.core.logP),
            J.constant(c,stable_sigma_point(c,sum(_bounds(c,s))/2),5),c.mpf(0),self.o.core.j)
        truth=leading_regular_frame(v['actual'],v['comparison'],v['z'],v['rho'],v['x'],
            v['de'],v['e'],v['P2'],J.constant(iv,sigma,5),
            iv.mpf([0,_bounds(c,self.o._width_error(iv.mpf(1)))[1]]),self.o.core.j_iv)
        def extract(frame):
            return dict(U=frame['U'][0].rows[0],m=frame['M'][0].rows[0],
                Q=frame['Q'][0].rows[0],U_y=frame['U'][1].rows[0]*frame['rho'],
                Q_y=frame['Q'][1].rows[0]*frame['rho'])
        result=extract(point),extract(truth),v['rho']
        self._backgrounds[key]=result
        return result

    def evaluate(self,Z,coordinate,*,selected_sc_multiple=False):
        c,iv,J=self.ctx,self.iv,self.J
        if (type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1
                or type(selected_sc_multiple) is not bool):
            raise ValueError('Interior decimal Z and boolean original-coordinate selector required')
        s,flat,recipe=self.micro._coordinate(coordinate,selected_sc_multiple)
        if _bounds(c,s)[1]>c.mpf('.25'):
            raise ValueError('This missing first-micro pressure slice ends at1/4')
        key=Z,coordinate,selected_sc_multiple
        if key in self._packets:return self._packets[key]
        local,prefix,frames=self.micro._covers(Z,s)
        _,B,hist=self.micro._correction_logs(local,prefix,flat)
        lp,ll=self.o.core.lp_iv,self.o.core.ll_iv
        errors=dict(U=[None if q is None else q-ll for q in B],
                    m=[None if q is None else q+lp for q in hist['m']])
        errors['U_Z'],errors['m_Z']=errors['U'][1:],errors['m'][1:]
        slow=self.loop['slow_phase_held_log_bounds']['B']
        fast=self.loop['phase_and_first_slow_phase_log_bounds']
        slow_y=[self.read_log(slow[k]) for k in ('y1_Z0','y1_Z1')]
        phase_y=[self.read_log(fast[k]) for k in ('B_phi','B_Z_phi')]
        slow_y=[None if flat or v is None else v+lp-self.micro.logN for v in slow_y]
        phase_y=[None if flat or v is None else v+lp for v in phase_y]
        errors['U_y']=[self.sum_logs(v) for v in zip(slow_y,phase_y)]
        coeff=pressure_coefficients(J,iv,Z,self.o.core.delta_iv,2)
        errors['Q']=self.combine(coeff['Q'],errors,2)
        errors['Q_y']=self.combine(coeff['Q_y'],errors,1)
        points,truths,rho=self.background(Z,s)
        fields={}
        for name,order in (('U',2),('m',3),('Q',2),('U_y',1),('Q_y',1)):
            # Production Taylor coefficients become ordinary derivatives once.
            rows=[self.o._summary(points[name][k]*math.factorial(k),truths[name][k]*math.factorial(k))
                  for k in range(order+1)]
            fields[name]=dict(signed_background_ordinary_Z=rows,
                finite_N_absolute_error_log_bounds=self.log_rows(errors[name][:order+1]),
                approximation='signed_background_ordinary_Z[*].approximation',
                total_error='background absolute_error_upper plus exp(finite_N log bound); exact_zero means no finite-N term',
                finite_N_log_bound_is_not_selected_signed_value=True)
        rpoint=c.mpf(4)
        endpoint=dict(signed_background_ordinary_Z=[self.o._summary(points['Q_y'][k]*rpoint,
            truths['Q_y'][k]*rho[0]) for k in (0,1)],
            log_physical_factor=self.o._summary(-self.o.core.logLambda,-ll),
            finite_N_coefficient_error_log_bounds=self.log_rows([
                None if q is None else q+iv.log(rho[0]) for q in errors['Q_y']]),
            definition='R*Q_y=epsilon*rho*Q_y; common epsilon factor multiplies background and finite-N error once')
        out=dict(Z=Z,coordinate=recipe,source_family=self.family,physical_pressure_rows=fields,
            R_Q_y_endpoint=endpoint,selected_integer=copy.deepcopy(self.micro.integer),
            selected_logN=self.micro._source_summary(self.micro.logN),
            finite_N_axial_radial_error_components=dict(slow_B_y_over_N=self.log_rows(slow_y),
                fast_B_phi=self.log_rows(phase_y),fast_Z_row_source='B_Z_phi',
                formula='delta U_y=Pstar*(B_y/N+B_phi); delta U_yZ=Pstar*(B_yZ/N+B_Z_phi)',
                phase_y='N',phase_Z='0',extra_hb_factor=False),
            source_cover=dict(exact_Z_cells=frames,local_source_cells=len(local),prefix_source_cells=len(prefix),
                exact_flat_collar=flat,generic_derivative_chart='bridge_first',
                actual_radius='R=4epsilon*exp(hb*s)',native_measure='dy=hb*ds',
                phase='frac(N*hb*(s-sc/2))',original_sc_tuple=self.t.sc_tuple),
            contract=dict(actual_selected_N_not_candidate_257=True,original_Pstar_applied_once=True,
                slow_B_y_divided_by_N=True,fast_B_phi_not_divided_by_N=True,
                fast_mixed_B_Z_phi_retained=True,ordinary_Z_factorial_once=True,
                mean_correction_history_from_rminus_retained=True,mean_y_equals_U_minus_m=True,
                Q_ZZ_includes_mean_ZZZ_error=True,large_error_exponentials_not_materialized=True,
                original_first_micro_phase_and_absolute_radius_preserved=True),
            limits=dict(actual_first_micro_pressure_rows_with_factored_errors=True,
                resolved_finite_N_signed_centers=False,finite_total_relative_error=False,
                signed_first_micro_m15_integrated=False,independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out
        return out

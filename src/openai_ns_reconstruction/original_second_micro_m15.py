"""Second-micro cumulative pressure with constant-hb leading controls.

The complete first-window moment and finite-N histories are retained.
Only the pressure increment on native [1,s] is added. Both modified-field
pressure endpoints remain in its error; the lower one is not a flat inlet.
"""
import copy
from fractions import Fraction
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_second_micro_numeric import OriginalSecondMicroNumeric
from .original_first_micro_full_m15 import OriginalFirstMicroFullPressureRows,OriginalFirstMicroFullM15
from .original_first_micro_m15 import OriginalFirstMicroM15
from .original_n1_inner_m15_numeric import leading_regular_frame,normalized_pressure_source


class OriginalSecondMicroPressureRows(OriginalFirstMicroFullPressureRows):
    def __init__(self,*,micro=None,**kwargs):
        if micro is not None and kwargs:raise ValueError('Supply second source or constructor options')
        micro=micro if micro is not None else OriginalSecondMicroNumeric(**kwargs)
        if not isinstance(micro,OriginalSecondMicroNumeric):raise ValueError('Actual original second-micro source required')
        super().__init__(micro=micro)
        path=self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop=json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart']['bridge_second']
        if (self.loop['source_family']!=self.family or self.loop['chart']!='bridge_second'
                or self.loop['primitive_units']['B']!='original B/Pstar; E=Utheta/Pstar'):
            raise ValueError('Original bridge_second derivative source/units required')
        A=self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if A is not None and _bounds(self.audit,self.micro.logN)[0]<_bounds(self.audit,A)[1]:
            raise ValueError('Actual selected N must admit second-chart A/N')
        expected=[([1,1],[5,4]),([5,4],[3,2]),([3,2],[7,4]),([7,4],[2,1])]
        for frame in self.micro.ranges:
            chart=next(v for v in frame['actual_17_chart_C3_transports'] if v['chart']=='second_micro')
            if [(v['exact_left'],v['exact_right']) for v in chart['cells']]!=expected:
                raise ValueError('Complete original second-micro partition required')
        self._second_inputs={}
        name='src/openai_ns_reconstruction/original_second_micro_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def inputs(self,Z,s):
        key=Z,s._mpi_
        if key in self._second_inputs:return self._second_inputs[key]
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;order=5
        inlet=o.inlet_jets(Z,degree=24,order=order)
        if Z not in o._prepared:
            with mp.workdps(280):o._prepared[Z]=o.bridge.prepare(Z)
        upper=_bounds(self.audit,s)[1]
        # These are cumulative FIELD changes: first window plus second
        # suffix. The moment integral below still runs only from1 to s.
        _,Bell,_,Bpsi=self.micro._second_increments(o._prepared[Z],upper,order)
        a=inlet['truth']['values']
        delta=dict(Phi=o._norm(a['Phi'])*Bell*iv.exp(o._width_error(Bell)),Psi=Bpsi,
            Mpsi=iv.mpf(upper)*o._norm(a['Psi']-a['Mpsi'])+o._width_error(iv.mpf(upper)*Bpsi))
        actual={name:a[name]+o._norm_error_jet(o._width_error(v),order) for name,v in delta.items()}
        comparison={name:[inlet['truth']['radial'][name][i]+o._norm_error_jet(
            o._width_error(self.t._core_comparison_changes(name,i,upper,order)),order) for i in range(3)]
            for name in ('Phi','Psi')}
        rho=J.constant(iv,iv.mpf([4,_bounds(self.audit,iv.mpf(4)+
            o._width_error(4*iv.mpf(upper)*iv.exp(o.hcap*iv.mpf(upper))))[1]]),order)
        out=dict(inlet=inlet,actual=actual,comparison=comparison,delta=delta,
            z=J.variable(iv,iv.mpf(Z),order),rho=rho,x=J.constant(iv,iv.sqrt(rho[0]),order),
            de=o.core.delta_iv,e=o.core.epsilon_iv,P2=iv.exp(2*o.core.lp_iv),
            prefix='full first control window plus second suffix s-1')
        self._second_inputs[key]=out
        return out

    def frame(self,Z,s,*,point):
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;v=self.inputs(Z,s)
        if point:
            inp=v['inlet']['point']
            return leading_regular_frame(inp['values'],inp['radial'],J.variable(c,c.mpf(Z),5),
                J.constant(c,4,5),J.constant(c,2,5),o.core.delta,o.core.epsilon,
                c.exp(2*o.core.logP),J.constant(c,1,5),c.mpf(0),o.core.j)
        hb=iv.mpf([0,_bounds(self.audit,o._width_error(iv.mpf(1)))[1]])
        # effective sigma=1 is a CONTROL switch: chi=hb constant. No
        # sigma derivative exists in this chart, so all radial rows are full.
        return leading_regular_frame(v['actual'],v['comparison'],v['z'],v['rho'],v['x'],
            v['de'],v['e'],v['P2'],J.constant(iv,1,5),hb,o.core.j_iv)

    def background(self,Z,s):
        key=Z,s._mpi_
        if key not in self._backgrounds:
            p,t=self.frame(Z,s,point=True),self.frame(Z,s,point=False)
            def extract(f):
                return dict(U=f['U'][0].rows[0],m=f['M'][0].rows[0],Q=f['Q'][0].rows[0],
                    U_y=f['U'][1].rows[0]*f['rho'],Q_y=f['Q'][1].rows[0]*f['rho'])
            self._backgrounds[key]=extract(p),extract(t),t['rho']
        return self._backgrounds[key]

    def evaluate(self,Z,coordinate,*,selected_sc_multiple=False):
        out=super().evaluate(Z,coordinate,selected_sc_multiple=selected_sc_multiple)
        out['source_cover'].update(generic_derivative_chart='bridge_second',native_domain=['1','2'],
            cumulative_history_charts=['first_micro','second_micro'],first_window_memory_retained=True)
        out['contract'].pop('original_first_micro_phase_and_absolute_radius_preserved',None)
        out['contract'].update(original_shared_micro_phase_and_absolute_radius_preserved=True,
            leading_chi_equals_hb_constant=True,first_window_field_error_retained_once=True)
        out['limits'].pop('actual_full_first_micro_pressure_rows_with_factored_errors',None)
        out['limits']['actual_second_micro_pressure_rows_with_factored_errors']=True
        return out


class OriginalSecondMicroM15(OriginalFirstMicroM15):
    def __init__(self,*,rows=None,**kwargs):
        if rows is not None and kwargs:raise ValueError('Supply retained second rows or constructor options')
        self.rows=rows if rows is not None else OriginalSecondMicroPressureRows(**kwargs)
        if not isinstance(self.rows,OriginalSecondMicroPressureRows):raise ValueError('Actual second pressure rows required')
        self.t,self.o=self.rows.t,self.rows.o
        self.ctx,self.iv,self.J,self.audit=self.rows.ctx,self.rows.iv,self.rows.J,self.rows.audit
        self.directory,self.family=self.rows.directory,self.rows.family
        self.hashes,self.code_hashes=dict(self.rows.hashes),dict(self.rows.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_first_micro_full_m15_numeric.json')
        raw=path.read_bytes();d=json.loads(raw)
        if (d['source_family']!=self.family or not d['checks']['passed']
                or not d['limits']['full_first_micro_leading_signed_m15_with_complete_factored_finite_N_errors']):
            raise ValueError('Accepted complete original first-micro inlet required')
        for root,group,target in ((self.directory,d['input_hashes'],self.hashes),(ROOT,d['implementation_hashes'],self.code_hashes)):
            for name,h in group.items():
                if target.get(name,h)!=h or hashlib.sha256((root/name).read_bytes()).hexdigest()!=h:
                    raise ValueError('Changed or conflicting second-micro inlet: '+name)
                target[name]=h
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest();self.inlet_report=path.name
        self.inlets={Fraction(p['Z']):p for p in d['samples']
            if p['coordinate']==dict(kind='exact_decimal_phase',value='1') and p['cells_per_piece']==8}
        self._first_engine=None;self._integrals,self._packets={},{}
        name='src/openai_ns_reconstruction/original_second_micro_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            if self._first_engine is None:self._first_engine=OriginalFirstMicroFullM15(cells=8,precision=self.o.precision,directory=self.directory)
            self.inlets[key]=self._first_engine.evaluate(Z,'1')
            self.hashes.update(self._first_engine.hashes);self.code_hashes.update(self._first_engine.code_hashes)
        p=self.inlets[key]
        if (p['source_family']!=self.family or Fraction(p['Z'])!=key
                or p['original_sc_tuple']!=self.t.sc_tuple or len(p['m15_C0_C1_terms'])!=13):
            raise ValueError('Same-Z complete thirteen-sector first-micro inlet required')
        return p

    def regular_integral(self,Z,s):
        key=Z,s._mpi_
        if key in self._integrals:return self._integrals[key]
        c,iv,J=self.ctx,self.iv,self.J
        # The current fields already have the full first-window history.
        # Their constant-control density bounds the entire second prefix.
        p,t=self.rows.frame(Z,s,point=True),self.rows.frame(Z,s,point=False)
        dp=-(normalized_pressure_source(p).rows[0]*p['rho']*p['rho']).truncate(1)
        dt=-(normalized_pressure_source(t).rows[0]*t['rho']*t['rho']).truncate(1)
        # Integrate the constant signed proxy exactly. The actual varying
        # source is enclosed uniformly over [1,s], so its entire variation
        # error is charged in dt. Artificially subdividing this same box
        # would not improve it or constitute a convergence experiment.
        width=s-1
        value=dp*(sum(_bounds(c,width))/2);truth=dt*width
        self._integrals[key]=value,truth
        return value,truth

    def evaluate(self,Z,coordinate):
        packet=self.rows.evaluate(Z,coordinate);s,_,recipe=self.rows.micro._coordinate(coordinate,False)
        key=Z,coordinate
        if key in self._packets:return self._packets[key]
        c,iv,J=self.ctx,self.iv,self.J;incoming=self.inlet(Z);terms=copy.deepcopy(incoming['m15_C0_C1_terms'])
        empty=s._mpi_==iv.mpf(1)._mpi_
        epsilon=self.o._summary(-self.o.core.logLambda,-self.o.core.ll_iv)
        eps2=self.o._summary(-2*self.o.core.logLambda,-2*self.o.core.ll_iv)
        P2=self.o._summary(2*self.o.core.logP,2*self.o.core.lp_iv)
        h=self.o._summary(sum(_bounds(c,self.o.logh))/2,self.o.logh);h['source_log_atom_only']=True
        zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1)
        def append(name,p,t,atoms,binding):
            terms.append(dict(name=name,coefficient_C0_C1=self.o._jet_summary(p.truncate(1),t.truncate(1)),
                log_scale_terms=copy.deepcopy(atoms),actual_function_binding=binding))
        binding=dict(interval='actual second_micro [1,S]',leading_control='chi=hb constant',
            physical_radius='R=4epsilon*exp(hb*s)',measure='dy=hb ds; R*dR=epsilon^2*rho^2*hb ds',
            first_window_moment_not_added_again=True,first_window_field_and_mean_history_retained=True)
        lower_packet=self.rows.evaluate(Z,'1')
        lower=self.finite_N_errors(Z,iv.mpf(1),lower_packet)
        if empty:
            pp,tt=zero,zeroi;error={name:[None,None] for name in ('endpoint','volume','viscous')}
        else:
            pp,tt=self.regular_integral(Z,s)
            upper=self.finite_N_errors(Z,s,packet)
            error=dict(volume=upper['volume'])
            for name in ('endpoint','viscous'):
                error[name]=[self.rows.sum_logs((upper[name][k],lower[name][k])) for k in (0,1)]
        append('second_micro_regular_pressure',pp,tt,[eps2,P2,h],dict(binding,density='-rho^2*Np ds',
            no_sigma_prime_or_inverse_hb_row=True,complete_directed_whole_prefix_density=True))
        for name,atoms in (('endpoint',[eps2]),('volume',[eps2,h]),('viscous',[epsilon])):
            bound=self.rows.sum_logs(error[name]);truth=zeroi
            if bound is not None:
                top=iv.mpf(_bounds(self.audit,bound)[1]);atom=self.o._summary(_bounds(c,top)[1],top)
                atom.update(positive_error_bound_atom_only=True,source_value_selected=False)
                atoms=atoms+[atom]
                if name=='volume':
                    length=iv.log(s-1);atoms=atoms+[self.o._summary(sum(_bounds(c,length))/2,length)]
                truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
            append('second_micro_finite_N_'+name+'_error',zero,truth,atoms,
                dict(binding,error_enclosure_only=True,both_upper_and_lower_pressure_endpoints_retained=True,
                     nonzero_first_window_correction_not_declared_flat=True,quadratic_cross_errors_retained=True,
                     ordinary_C0_C1_bounds_not_derivatives_of_caps=True))
        out=dict(Z=Z,coordinate=recipe,source_family=self.family,m15_C0_C1_terms=terms,
            incoming_first_micro_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),inlet_source_report=self.inlet_report,
            integration_method='constant signed density proxy integrated exactly; complete uniform source error',
            original_sc_tuple=self.t.sc_tuple,exact_empty_increment=empty,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']),
            second_micro_full_prefix_cover=copy.deepcopy(packet['source_cover']),
            finite_N_error_logs_before_physical_measures={k:self.rows.log_rows(v) for k,v in error.items()},
            incoming_finite_N_endpoint_error_logs={k:self.rows.log_rows(lower[k]) for k in ('endpoint','viscous')},
            contract=dict(complete_first_window_m15_retained_once=True,only_second_suffix_moment_integrated=True,
                cumulative_first_plus_second_field_and_mean_errors=True,actual_constant_hb_control=True,
                no_new_sigma_prime_singular_increment=True,previous_singular_pressure_history_not_reset=True,
                both_finite_N_pressure_endpoints_and_Z1_retained=True,nonzero_lower_correction_history_retained=True,
                original_selected_N_phase_and_absolute_radius_preserved=True,physical_hb_ds_and_RdR_once=True,
                all_source_cells_and_axial_frames_in_uniform_prefix_error=True,raw_n1_transport_not_extended=True),
            limits=dict(second_micro_leading_signed_m15_with_complete_factored_finite_N_errors=True,
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                remaining_active_compact_signed_pressure=False,independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out
        return out

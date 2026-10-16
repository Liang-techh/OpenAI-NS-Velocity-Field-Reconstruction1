"""Original pressure through O2 buffer and both O3 charts, with full C3 errors."""
import copy
from fractions import Fraction
import hashlib
import json
import math
from types import SimpleNamespace

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_axial_pressure_rows import OriginalAxialPressureRows
from .original_outer_pressure_rows import OriginalOuterPressureRows
from .original_outer_buffer_numeric import OriginalOuterBufferNumeric
from .original_outer_o3_numeric import OriginalOuterO3Numeric
from .original_Rm_patch_numeric import encode


class OriginalPostAxialPressureRows(OriginalAxialPressureRows):
    DOMAINS={'O2_buffer':(Fraction(0),Fraction(11)),
        'O3_transition':(Fraction(0),Fraction(1)),'O3_power':(Fraction(0),Fraction(2))}
    ORIGINS={'O2_buffer':0,'O3_transition':11,'O3_power':12}
    PREVIOUS={'O2_buffer':('O2_axial','1'),'O3_transition':('O2_buffer','11'),'O3_power':('O3_transition','1')}
    LOOP_CHARTS={'O2_buffer':'O2_buffer','O3_transition':'O3_slope_mu','O3_power':'O3_power'}
    history_logs=OriginalOuterPressureRows.history_logs

    def __init__(self,*,chart='O2_buffer',previous=None,**kwargs):
        if chart not in self.DOMAINS:raise ValueError('Original post-axial chart required')
        if previous is not None and kwargs:raise ValueError('Supply predecessor or constructor options')
        if previous is None:
            if chart!='O2_buffer':raise ValueError('Complete original predecessor required')
            previous=OriginalAxialPressureRows(**kwargs)
        expected,self.incoming_coordinate=self.PREVIOUS[chart]
        if previous.chart!=expected:raise ValueError('Contiguous original post-axial predecessor required')
        self.chart,self.previous=chart,previous
        for name in ('ctx','iv','J','audit','directory','family','o','t','micro','owner','reference','cover_owner','patch_rows'):
            setattr(self,name,getattr(previous,name))
        self.hashes,self.code_hashes=dict(previous.hashes),dict(previous.code_hashes)
        if chart=='O2_buffer':
            self.axial=previous
            path=self.directory/(PREFIX+'compliant_current_original_axial_m15_numeric.json');raw=path.read_bytes();data=json.loads(raw)
            if not data['checks']['passed'] or data['source_family']!=self.family:raise ValueError('Accepted actual axial pressure required')
            for base,group,target in ((self.directory,data['input_hashes'],self.hashes),(ROOT,data['implementation_hashes'],self.code_hashes)):
                for name,digest in group.items():
                    if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest or target.get(name,digest)!=digest:
                        raise ValueError('Changed original post-axial prerequisite: '+name)
                    target[name]=digest
            self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
            endpoints=[p for p in data['pressure_row_samples'] if p['native_coordinate']=='1']
            if not endpoints:raise ValueError('Positive original axial terminal mean required')
            self.mean_source_summary=copy.deepcopy(endpoints[0]['leading_scalars']['mean'])
            if any(p['leading_scalars']['mean']!=self.mean_source_summary for p in endpoints):
                raise ValueError('Original scalar mean must be independent of Z')
            self.mean_source_proof=copy.deepcopy(endpoints[0]['positive_mean_memory_receipt'])
            self.mean_point=self.ctx.mpf(self.mean_source_summary['approximation'])
            self.mean_truth=self.iv.mpf([self.ctx.make_mpf(tuple(v)) for v in self.mean_source_summary['source_enclosure_exact_mpf_tuples']])
            if _bounds(self.audit,self.mean_truth)[0]<=0:raise ValueError('Original nonzero terminal memory must be preserved')
            self.axial_inlets={Fraction(p['Z']):p for p in data['samples'] if p['native_coordinate']=='1'}
        else:
            for name in ('axial','mean_source_summary','mean_source_proof','mean_point','mean_truth','axial_inlets'):
                setattr(self,name,getattr(previous,name))
        path=self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        loop_chart=self.LOOP_CHARTS[chart]
        self.loop=json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart'][loop_chart]
        if self.loop['chart']!=loop_chart or self.loop['source_family']!=self.family or self.loop['primitive_units']['B']!='original B/Pstar; E=Utheta/Pstar':
            raise ValueError('Actual post-axial slow/fast pressure units required')
        A=self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if A is not None and _bounds(self.audit,self.owner.logN)[0]<_bounds(self.audit,A)[1]:raise ValueError('Original N must admit A/N')
        cls=OriginalOuterBufferNumeric if chart=='O2_buffer' else OriginalOuterO3Numeric
        utility=object.__new__(cls)
        for name in ('ctx','iv','family','owner'):setattr(utility,name,getattr(self,name))
        self.utility=utility if chart=='O2_buffer' else SimpleNamespace(_cover=lambda Z,q:utility._cover(Z,q,self.chart))
        self._histories={};self._incoming={};self._backgrounds={};self._packets={};self._errors={};self._scalars={};self._uniforms={};self._next=None
        for name in ('original_post_axial_pressure_rows.py','original_outer_buffer_numeric.py','original_outer_o3_numeric.py'):
            name='src/openai_ns_reconstruction/'+name;self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):
        if self.chart=='O3_power':raise ValueError('Next actual source is original compact continuation')
        if self._next is None:
            chart='O3_transition' if self.chart=='O2_buffer' else 'O3_power'
            self._next=OriginalPostAxialPressureRows(chart=chart,previous=self)
        return self._next

    def geometry(self,value):
        q,p,t=self.coordinate(value);c,iv=self.ctx,self.iv;offset=self.ORIGINS[self.chart]
        terms=list(self.axial.geometry('1')['atomic_log_radius_terms'])+[(c.mpf(offset),iv.mpf(offset)),(p,t)]
        return dict(native=q,coordinate_point=p,coordinate=t,local_point=p,local=t,
            physical_jacobian_point=c.mpf(1),physical_jacobian=iv.mpf(1),
            atomic_log_radius_terms=terms,logR_point=sum((a for a,b in terms),c.mpf(0)),
            logR=sum((b for a,b in terms),iv.mpf(0)),global_post_axial_offset=offset)

    def scalars(self,value):
        if value not in self._scalars:
            q,p,t=self.coordinate(value);offset=self.ORIGINS[self.chart]
            self._scalars[value]=dict(point=dict(beta=self.ctx.mpf(0),beta_y=self.ctx.mpf(0),mean=self.mean_point*self.ctx.exp(-offset-p)),
                truth=dict(beta=self.iv.mpf(0),beta_y=self.iv.mpf(0),mean=self.mean_truth*self.iv.exp(-offset-t)))
        return self._scalars[value]

    def uniform_background(self,Z,value):
        key=Z,value
        if key not in self._uniforms:
            end=self.scalars(value)['truth']['mean'];start=self.scalars('0')['truth']['mean']
            mean=self.iv.mpf([max(0,_bounds(self.audit,end)[0]),_bounds(self.audit,start)[1]])
            self._uniforms[key]=self.fields(Z,self.iv,self.iv.mpf(0),mean,self.iv.mpf(0))
        return self._uniforms[key]

    def evaluate(self,Z,value):
        key=Z,value
        if key in self._packets:return self._packets[key]
        if type(Z) is not str or not self.ctx.isfinite(self.ctx.mpf(Z)) or not -1<self.ctx.mpf(Z)<1:raise ValueError('Interior decimal Z required')
        p,t,geo=self.background(Z,value);h=self.history_logs(Z,value);errors,components=self.error_data(Z,value)
        fields={name:dict(signed_background_ordinary_Z=[self.o._summary(p[name][k]*math.factorial(k),t[name][k]*math.factorial(k)) for k in range(order+1)],
            finite_N_absolute_error_log_bounds=self.log_rows(errors[name][:order+1]),finite_N_log_bound_is_not_selected_signed_value=True)
            for name,order in (('U',2),('m',3),('Q',2),('U_y',1),('Q_y',1))}
        v=self.scalars(value)
        out=dict(Z=Z,chart=self.chart,native_coordinate=value,source_family=self.family,selected_integer=copy.deepcopy(self.owner.integer),
            physical_pressure_rows=fields,leading_mean=self.o._summary(v['point']['mean'],v['truth']['mean']),
            positive_original_axial_mean_source=copy.deepcopy(self.mean_source_summary),positive_original_axial_mean_proof=copy.deepcopy(self.mean_source_proof),
            finite_N_history_error_logs={name:self.log_rows(values) for name,values in h['history'].items()},
            finite_N_incoming_history_error_logs={name:self.log_rows(values) for name,values in self.incoming_history_logs(Z).items()},
            finite_N_prefix_history_error_logs={name:self.log_rows(values) for name,values in h['prefix_history'].items()},
            finite_N_axial_radial_error_components=components,physical_transport_steps=encode(h['steps']),
            geometry=dict(atomic_log_radius_terms=[self.o._summary(a,b) for a,b in geo['atomic_log_radius_terms']],
                physical_local_log_distance=self.o._summary(geo['local_point'],geo['local']),
                global_post_axial_origin=self.ORIGINS[self.chart],physical_Jacobian=1,native_coordinate_is_physical_log_offset=True,
                physical_radius_materialized=False,exact_start=geo['native']==0),
            source_cover=dict(exact_Z_cells=h['frames'],unique_native_cells=len(h['groups']),local_source_cells=len(h['local']),
                generic_derivative_chart=self.LOOP_CHARTS[self.chart],transport_chart=self.chart,
                native_domain=[str(v) for v in self.DOMAINS[self.chart]],all_upstream_histories_retained=True),
            contract=dict(original_U_zero_M_4Z_positive_decay=True,exact_rational_pressure_recovery=True,
                original_Q_and_Qy_nonzero_opposites=True,analytic_mean_Z3_and_QZZ=True,
                all_four_history_rows_retained=True,actual_chart_slow_fast_units=True,
                chart_and_predecessor_endpoint_budgets=True,physical_native_measure_once=True,
                O3_mu_angular_sources_retained_in_original_finite_N_bounds=True,no_false_zero_mean=True),
            limits=dict(actual_post_axial_pressure_with_complete_absolute_errors=True,
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

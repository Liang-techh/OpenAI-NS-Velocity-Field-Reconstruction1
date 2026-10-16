"""Original cumulative pressure across buffer/O3 with exact leading cancellation."""
import copy
from fractions import Fraction
import hashlib

from .original_core_numeric import ROOT,_bounds
from .original_axial_m15 import OriginalAxialM15
from .original_post_axial_pressure_rows import OriginalPostAxialPressureRows


class OriginalPostAxialM15(OriginalAxialM15):
    COUNTS={'O2_buffer':195,'O3_transition':199,'O3_power':203}

    def __init__(self,*,rows=None,previous=None,**kwargs):
        if rows is not None and kwargs:raise ValueError('Supply original rows or constructor options')
        self.rows=rows if rows is not None else OriginalPostAxialPressureRows(**kwargs)
        self.chart=self.rows.chart;self.previous=previous
        if self.chart!='O2_buffer' and (previous is None or previous.chart!=self.rows.PREVIOUS[self.chart][0]):
            raise ValueError('Complete actual post-axial moment predecessor required')
        for name in ('t','o','ctx','iv','J','audit','directory','family'):setattr(self,name,getattr(self.rows,name))
        self.hashes,self.code_hashes=dict(self.rows.hashes),dict(self.rows.code_hashes)
        if previous is not None:self.hashes.update(previous.hashes);self.code_hashes.update(previous.code_hashes)
        self.inlets=dict(self.rows.axial_inlets) if self.chart=='O2_buffer' else {}
        self._axial_engine=None;self._packets={};self._next=None
        name='src/openai_ns_reconstruction/original_post_axial_m15.py';self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):
        if self.chart=='O3_power':raise ValueError('Next pressure source is the original compact interval')
        if self._next is None:self._next=OriginalPostAxialM15(rows=self.rows.next_field(),previous=self)
        return self._next

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            if self.previous is None:
                if self._axial_engine is None:self._axial_engine=OriginalAxialM15(rows=self.rows.axial)
                source,coordinate=self._axial_engine,'1'
            else:source,coordinate=self.previous,self.rows.incoming_coordinate
            self.inlets[key]=source.evaluate(Z,coordinate)
            self.hashes.update(source.hashes);self.code_hashes.update(source.code_hashes)
        return self.inlets[key]

    def evaluate(self,Z,value):
        key=Z,value
        if key in self._packets:return self._packets[key]
        packet=self.rows.evaluate(Z,value);incoming=self.inlet(Z);count=self.COUNTS[self.chart]
        if len(incoming['m15_C0_C1_terms'])!=count or incoming['actual_selected_integer']!=packet['selected_integer']:
            raise ValueError('Complete same-N original pressure inlet required')
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;geo=self.rows.geometry(value);empty=geo['native']==0
        zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1);terms=copy.deepcopy(incoming['m15_C0_C1_terms'])
        # The global axial inlet is common to all charts. At an O3 lower
        # endpoint Rlower^2 is larger but its mean squared is smaller by the
        # reciprocal factor; do not mix Rlower with the global mstar.
        anchor=self.rows.axial.geometry('1')
        atoms=[o._summary(2*a,2*b) for a,b in anchor['atomic_log_radius_terms']]
        mp,mt=self.rows.mean_point,self.rows.mean_truth
        coefficientp=J.constant(c,-4*mp*mp*geo['local_point'],1)
        coefficientt=J.constant(iv,-4*mt*mt*geo['local'],1)
        def append(name,p,t,scale,**binding):
            terms.append(dict(name=self.chart+'_'+name,coefficient_C0_C1=o._jet_summary(zero if empty else p,zeroi if empty else t),
                log_scale_terms=[] if empty else copy.deepcopy(scale),
                actual_function_binding=dict(chart=self.chart,incoming_terms_retained_once=count,**binding)))
        append('exact_mean_square_integral',coefficientp,coefficientt,atoms,
            exact_primitive='-4*R_axial_end^2*mstar^2*local_length',
            common_source_function='original_positive_axial_terminal_mean',global_mean_and_radius_anchor=True,
            leading_linear_pressure_terms_cancel_exactly=True,half_S1_already_included=True,
            both_leading_RQy_endpoints_nonzero_and_equal=True,no_new_free_parameter=True)
        if empty:finite={name:[None,None] for name in ('endpoint','volume','viscous')};budget=None
        else:finite,budget=self.finite_errors(Z,value,packet)
        for name in ('endpoint','volume','viscous'):
            bound=self.rows.sum_logs(finite[name]);truth=zeroi;scale=[]
            if bound is not None:
                top=iv.mpf(_bounds(self.audit,bound)[1]);atom=o._summary(_bounds(c,top)[1],top)
                atom.update(positive_error_bound_atom_only=True,source_value_selected=False);scale.append(atom)
                truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
            append('finite_N_'+name+'_error',zero,truth,scale,error_enclosure_only=True,
                both_actual_endpoints=True,full_prefix_C3_and_quadratic_errors=True)
        out=dict(Z=Z,chart=self.chart,native_coordinate=value,source_family=self.family,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']),m15_C0_C1_terms=terms,
            incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']),original_global_mean_source=copy.deepcopy(self.rows.mean_source_summary),
            leading_RQy_endpoint_difference_exact_zero=True,leading_RQy_endpoints_individually_nonzero=True,
            finite_N_error_logs_in_physical_units={name:self.rows.log_rows(v) for name,v in finite.items()},
            actual_lower_finite_N_endpoint_budget=budget,
            contract=dict(complete_incoming_cumulative_state_once=True,original_positive_mean_not_reset=True,
                exact_half_S1_primitive=True,global_radius_and_global_mean_paired=True,
                no_extra_O3_exp22_or_exp24_factor=True,nonzero_RQy_cancels_as_same_source_function=True,
                leading_viscous_cancellation_does_not_reset_errors=True,both_finite_N_endpoints=True,
                full_prefix_not_final_decayed_history=True,actual_mu_angular_error_sources_retained=True,
                physical_native_measure_once=True,no_huge_radius_materialized=True),
            limits=dict(original_post_axial_cumulative_m15_with_complete_absolute_errors=True,
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                independent_reserved_Rp_seam=False,global_five_moment_numeric_join=False,common_global_N=False,
                exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

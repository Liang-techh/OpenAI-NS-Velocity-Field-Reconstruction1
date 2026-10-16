"""Complete Rc pressure inlet through actual compact repair and exact quiet flow."""
import copy
from fractions import Fraction
import hashlib

from .original_core_numeric import ROOT,_bounds
from .original_post_axial_m15 import OriginalPostAxialM15
from .original_compact_pressure_rows import OriginalCompactPressureRows


class OriginalCompactQuietM15(OriginalPostAxialM15):
    COUNTS={'compact':207,'quiet':211}

    def __init__(self,*,rows=None,previous=None,**kwargs):
        if rows is not None and kwargs:raise ValueError('Supply compact rows or constructor options')
        self.rows=rows if rows is not None else OriginalCompactPressureRows(**kwargs)
        self.chart=self.rows.chart;self.previous=previous
        if self.chart=='quiet' and (previous is None or previous.chart!='compact'):raise ValueError('Complete original 2Rc moment predecessor required')
        for name in ('t','o','ctx','iv','J','audit','directory','family'):setattr(self,name,getattr(self.rows,name))
        self.hashes,self.code_hashes=dict(self.rows.hashes),dict(self.rows.code_hashes)
        if previous is not None:self.hashes.update(previous.hashes);self.code_hashes.update(previous.code_hashes)
        self.inlets=dict(self.rows.rc_inlets) if self.chart=='compact' else {}
        self._outer_engine=None;self._packets={};self._next=None
        name='src/openai_ns_reconstruction/original_compact_quiet_m15.py';self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):
        if self.chart!='compact':raise ValueError('Next is independent Rp pressure matching')
        if self._next is None:self._next=OriginalCompactQuietM15(rows=self.rows.next_field(),previous=self)
        return self._next

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            if self.previous is not None:source,coordinate=self.previous,'2'
            else:
                if self._outer_engine is None:
                    power=self.rows.previous;transition=power.previous;buffer=transition.previous
                    b=OriginalPostAxialM15(rows=buffer)
                    t=OriginalPostAxialM15(rows=transition,previous=b)
                    self._outer_engine=OriginalPostAxialM15(rows=power,previous=t)
                source,coordinate=self._outer_engine,'2'
            self.inlets[key]=source.evaluate(Z,coordinate);self.hashes.update(source.hashes);self.code_hashes.update(source.code_hashes)
        return self.inlets[key]

    def finite_errors(self,Z,value,packet):
        rows=self.rows;_,base,geo=rows.background(Z,value);_,lower_base,lower_geo=rows.background(Z,rows.START)
        incoming,_=rows.error_data(Z,rows.START);upper=self.correction_densities(Z,packet,base)
        lower=self.correction_densities(Z,self.packet_for_errors(incoming),lower_base)
        prefix,_=rows.error_data(Z,value,prefix=True)
        volume=self.correction_densities(Z,self.packet_for_errors(prefix),rows.uniform_background(Z,value))['volume']
        result={}
        for name,power in (('endpoint',2),('viscous',1)):
            result[name]=[rows.sum_logs((None if upper[name][k] is None else upper[name][k]+power*geo['logR'],
                None if lower[name][k] is None else lower[name][k]+power*lower_geo['logR'])) for k in (0,1)]
        mass=-self.iv.expm1(-2*geo['local'])/2
        result['volume']=[None if v is None else v+2*geo['logR']+self.iv.log(mass) for v in volume]
        return result,dict(endpoint=rows.log_rows(lower['endpoint']),viscous=rows.log_rows(lower['viscous']),
            prefix_error_rows={name:rows.log_rows(v) for name,v in prefix.items() if name in ('U','m','Q','U_y','Q_y')},
            normalized_positive_volume_mass=self.o._summary(-self.ctx.expm1(-2*geo['local_point'])/2,mass),
            both_actual_endpoints=True,normalization='Rupper²; physical dy once')

    def evaluate(self,Z,value):
        key=Z,value
        if key in self._packets:return self._packets[key]
        packet=self.rows.evaluate(Z,value);incoming=self.inlet(Z);count=self.COUNTS[self.chart]
        if len(incoming['m15_C0_C1_terms'])!=count or incoming['actual_selected_integer']!=packet['selected_integer']:
            raise ValueError('Complete same-N compact/quiet cumulative inlet required')
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;geo=self.rows.geometry(value);empty=geo['exact_start']
        zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1);terms=copy.deepcopy(incoming['m15_C0_C1_terms'])
        anchor=self.rows.axial.geometry('1');atoms=[o._summary(2*a,2*b) for a,b in anchor['atomic_log_radius_terms']]
        coefficientp=J.constant(c,-4*self.rows.mean_point**2*geo['local_point'],1)
        coefficientt=J.constant(iv,-4*self.rows.mean_truth**2*geo['local'],1)
        def append(name,p,t,scale,**binding):
            terms.append(dict(name=self.chart+'_'+name,coefficient_C0_C1=o._jet_summary(zero if empty else p,zeroi if empty else t),
                log_scale_terms=[] if empty else copy.deepcopy(scale),actual_function_binding=dict(chart=self.chart,
                    incoming_terms_retained_once=count,**binding)))
        append('exact_original_mean_square_integral',coefficientp,coefficientt,atoms,
            exact_primitive='-4*R_axial_end²*mstar²*physical_local_log_length',common_source_function='original_positive_axial_terminal_mean',
            half_S1_already_included=True,global_radius_and_global_mean_paired=True,
            leading_RQy_endpoints_nonzero_and_equal=True,original_quiet_K_identity_used=self.chart=='quiet')
        budget=None
        if self.chart=='compact' and not empty:finite,budget=self.finite_errors(Z,value,packet)
        else:finite={name:[None,None] for name in ('endpoint','volume','viscous')}
        if self.chart=='compact':
            for name in ('endpoint','volume','viscous'):
                bound=self.rows.sum_logs(finite[name]);truth=zeroi;scale=[]
                if bound is not None:
                    top=iv.mpf(_bounds(self.audit,bound)[1]);atom=o._summary(_bounds(c,top)[1],top)
                    atom.update(positive_error_bound_atom_only=True,source_value_selected=False);scale.append(atom)
                    truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
                append('finite_N_'+name+'_error',zero,truth,scale,error_enclosure_only=True,
                    actual_supported_compact_control_sources=True,full_prefix_C3_and_quadratic_errors=True)
        else:
            errors,_=self.rows.error_data(Z,value)
            if any(v is not None for name in ('U','m','Q','U_y','Q_y') for v in errors[name]):raise ValueError('Exact full quiet source correction must vanish')
        out=dict(Z=Z,chart=self.chart,native_coordinate=value,source_family=self.family,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']),m15_C0_C1_terms=terms,
            incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']),finite_N_error_logs_in_physical_units={name:self.rows.log_rows(v) for name,v in finite.items()},
            actual_lower_finite_N_endpoint_budget=budget,
            contract=dict(complete_incoming_cumulative_state_once=True,exact_half_S1_leading_primitive=True,
                positive_original_mean_and_absolute_pressure_not_reset=True,global_mean_radius_anchor=True,
                compact_error_caps_not_selected_signed_centers=True,physical_log_measure_once=True,
                compact_terminal_zero_is_relative_source_identity=True,quiet_full_Qy_equals_minus_Q_from_source=True,
                earlier_cumulative_errors_retained_in_quiet=True,no_huge_radius_materialized=True),
            limits=dict(original_compact_quiet_forward_m15_with_absolute_errors=True,
                exact_quiet_increment_from_full_source=self.chart=='quiet',resolved_signed_compact_finite_N_centers=False,
                useful_total_relative_accuracy=False,independent_reserved_Rp_seam=False,global_five_moment_numeric_join=False,
                common_global_N=False,exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

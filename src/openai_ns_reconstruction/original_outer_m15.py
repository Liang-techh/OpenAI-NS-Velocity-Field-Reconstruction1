"""Continue original cumulative m15 through Rh-reference and O2 angular slope."""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_outer_pressure_rows import OriginalOuterPressureRows
from .original_macro_m15 import OriginalMacroM15,half_S1_coefficients,densities,combine


class OriginalOuterM15(OriginalMacroM15):
    def __init__(self,*,rows=None,previous=None,**kwargs):
        if rows is not None and kwargs:raise ValueError('Supply retained original outer rows or constructor options')
        self.rows=rows if rows is not None else OriginalOuterPressureRows(**kwargs)
        self.chart=self.rows.chart;self.previous=previous
        if self.chart=='O2_slope' and (previous is None or previous.chart!='Rh_reference'):
            raise ValueError('Complete original outer reference moment predecessor required')
        for name in ('t','o','ctx','iv','J','audit','directory','family'):setattr(self,name,getattr(self.rows,name))
        self.hashes,self.code_hashes=dict(self.rows.hashes),dict(self.rows.code_hashes)
        if previous is not None:self.hashes.update(previous.hashes);self.code_hashes.update(previous.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_Rm_m15_numeric.json');raw=path.read_bytes();data=json.loads(raw)
        if not data['checks']['passed'] or data['source_family']!=self.family:
            raise ValueError('Accepted complete168-term Rh pressure moment required')
        for base,group,target in ((self.directory,data['input_hashes'],self.hashes),(ROOT,data['implementation_hashes'],self.code_hashes)):
            for name,digest in group.items():
                if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest or target.get(name,digest)!=digest:
                    raise ValueError('Changed original outer moment prerequisite: '+name)
                target[name]=digest
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        self.inlets={Fraction(p['Z']):p for p in data['samples'] if p['native_x']=='Rh'} if self.chart=='Rh_reference' else {}
        self._patch_engine=None;self._packets={};self._densities={};self._next=None
        name='src/openai_ns_reconstruction/original_outer_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):
        if self.chart!='Rh_reference':raise ValueError('O2 axial turnoff needs its genuine changing-source mean')
        if self._next is None:self._next=OriginalOuterM15(rows=self.rows.next_field(),previous=self)
        return self._next

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            if self.previous is not None:source=self.previous;coordinate='0'
            else:
                if self._patch_engine is None:
                    from .original_Rm_m15 import OriginalRmM15
                    self._patch_engine=OriginalRmM15(rows=self.rows.patch_rows)
                source=self._patch_engine;coordinate='Rh'
            self.inlets[key]=source.evaluate(Z,coordinate)
            self.hashes.update(source.hashes);self.code_hashes.update(source.code_hashes)
        return self.inlets[key]

    def base_density(self,Z):
        if Z not in self._densities:
            result={};left=str(self.rows.DOMAINS[self.chart][0])
            point,truth,_=self.rows.background(Z,left)
            for label,ctx,v,de in (('point',self.ctx,point,self.o.core.delta),('truth',self.iv,truth,self.o.core.delta_iv)):
                bc,vc=half_S1_coefficients(self.J,ctx,Z,de)
                values=densities({(0,0):v['U'].truncate(2)},{(0,0):v['Q']},self.J)
                result[label]=2*combine(bc,values)[(0,0)]+combine(vc,values)[(0,0)]
            self._densities[Z]=result
        return self._densities[Z]

    def packet_for_errors(self,errors):
        return dict(physical_pressure_rows={name:dict(finite_N_absolute_error_log_bounds=self.rows.log_rows(errors[name]))
                    for name in ('U','Q','Q_y')})

    def finite_errors(self,Z,value,packet):
        rows=self.rows;_,base,geo=rows.background(Z,value);left=str(rows.DOMAINS[self.chart][0])
        _,lower_base,lower_geo=rows.background(Z,left);incoming,_=rows.error_data(Z,left)
        upper=self.correction_densities(Z,packet,base)
        lower=self.correction_densities(Z,self.packet_for_errors(incoming),lower_base)
        prefix,_=rows.error_data(Z,value,prefix=True)
        volume=self.correction_densities(Z,self.packet_for_errors(prefix),rows.uniform_background(Z,value))['volume']
        result={}
        for name,power in (('endpoint',2),('viscous',1)):
            result[name]=[rows.sum_logs((None if upper[name][k] is None else upper[name][k]+power*geo['logR'],
                None if lower[name][k] is None else lower[name][k]+power*lower_geo['logR'])) for k in (0,1)]
        mass=self.iv.expm1(2*geo['local'])/2
        result['volume']=[None if v is None else v+2*lower_geo['logR']+self.iv.log(mass) for v in volume]
        return result,dict(endpoint=rows.log_rows(lower['endpoint']),viscous=rows.log_rows(lower['viscous']),
            log_radius=self.o._summary(lower_geo['logR_point'],lower_geo['logR']),
            prefix_error_rows={name:rows.log_rows(v) for name,v in prefix.items() if name in ('U','m','Q','U_y','Q_y')},
            exact_positive_native_mass=self.o._summary(self.ctx.expm1(2*geo['local_point'])/2,mass))

    def evaluate(self,Z,value):
        packet=self.rows.evaluate(Z,value);key=Z,value
        if key in self._packets:return self._packets[key]
        incoming=self.inlet(Z);count=168 if self.chart=='Rh_reference' else 172
        if len(incoming['m15_C0_C1_terms'])!=count or incoming['actual_selected_integer']!=packet['selected_integer']:
            raise ValueError('Complete same-N incoming original outer pressure moment required')
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;geo=self.rows.geometry(value);left=self.rows.DOMAINS[self.chart][0]
        empty=geo['native']==left;zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1)
        terms=copy.deepcopy(incoming['m15_C0_C1_terms']);base=self.base_density(Z)
        kp,kt=c.expm1(2*geo['local_point'])/2,iv.expm1(2*geo['local'])/2
        ref=self.rows.reference;originp=c.mpf(left.numerator)/left.denominator;origint=iv.mpf(left.numerator)/left.denominator
        atoms=[o._summary(2*c.log(110),2*iv.log(110)),o._summary(20*ref.logCp,20*ref.logC),
               o._summary(20*ref.logPp,20*ref.logP),o._summary(2*originp,2*origint)]
        def append(name,p,t,scale,**binding):
            terms.append(dict(name=self.chart+'_'+name,coefficient_C0_C1=o._jet_summary(p.truncate(1),t.truncate(1)),
                log_scale_terms=scale,actual_function_binding=dict(chart=self.chart,incoming_terms_retained_once=count,
                    original_constant_axial_pressure_source=True,**binding)))
        append('constant_background_integral',zero if empty else base['point']*kp,
            zeroi if empty else base['truth']*kt,atoms,exact_primitive='Rlower²*(2F+G)*expm1(2length)/2',
            correlated_source_function='original_U_M_4Z_half_S1',native_physical_Jacobian=1,
            same_angular_dependent_finite_N_sources_separate=True,no_large_endpoint_subtraction=True)
        if not empty:finite,lower=self.finite_errors(Z,value,packet)
        else:finite={name:[None,None] for name in ('endpoint','volume','viscous')};lower=None
        for name in ('endpoint','volume','viscous'):
            bound=self.rows.sum_logs(finite[name]);truth=zeroi;scale=[]
            if bound is not None:
                top=iv.mpf(_bounds(self.audit,bound)[1]);atom=o._summary(_bounds(c,top)[1],top)
                atom.update(positive_error_bound_atom_only=True,source_value_selected=False);scale.append(atom)
                truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
            append('finite_N_'+name+'_error',zero,truth,scale,error_enclosure_only=True,both_nonflat_endpoints=True,
                full_prefix_C3_history_and_local_slow_fast_bounds=True,no_epsilon_or_hb_scale=True)
        out=dict(Z=Z,chart=self.chart,native_coordinate=value,source_family=self.family,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']),m15_C0_C1_terms=terms,
            incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']),full_prefix_cover=copy.deepcopy(packet['source_cover']),
            positive_native_kernel=o._summary(kp,kt),leading_signed_RQy_endpoints_exact_zero=True,
            finite_N_error_logs_in_physical_units={name:self.rows.log_rows(v) for name,v in finite.items()},
            actual_lower_finite_N_endpoint_budget=lower,
            contract=dict(complete_incoming_cumulative_state_once=True,exact_original_U_M_4Z_identity=True,
                analytic_mean_Z3_and_QZZ=True,exact_background_positive_native_integral=True,
                original_half_S1_and_RQy_sign=True,leading_Qy_zero_does_not_reset_errors=True,
                angular_effects_in_all_finite_N_sources_retained=True,both_finite_N_endpoints=True,
                full_prefix_not_final_decayed_history=True,new_chart_local_slow_fast_bounds_at_seams=True,
                native_measure_once=True,unit_radius_offsets_separate=True),
            limits=dict(outer_reference_cumulative_signed_m15_with_complete_absolute_errors=True,
                O2_slope_cumulative_signed_m15_with_complete_absolute_errors=self.chart=='O2_slope',
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                independent_preRp_signed_seam=False,global_five_moment_numeric_join=False,common_global_N=False,
                exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

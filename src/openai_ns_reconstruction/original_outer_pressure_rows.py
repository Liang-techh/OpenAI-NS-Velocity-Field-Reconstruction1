"""Genuine pressure/mean rows on the original Rh-reference and O2-slope charts."""
import copy
from fractions import Fraction
import hashlib
import json
import math

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_pressure_rows import OriginalRmPressureRows
from .original_Rm_finite_N_numeric import RATES
from .original_outer_reference_numeric import OriginalOuterReferenceNumeric
from .original_outer_slope_numeric import OriginalOuterSlopeNumeric
from .original_first_micro_pressure_rows import OriginalFirstMicroPressureRows,pressure_coefficients


class OriginalOuterPressureRows(OriginalFirstMicroPressureRows):
    DOMAINS={'Rh_reference':(Fraction(-5),Fraction(0)),'O2_slope':(Fraction(0),Fraction(1))}

    def __init__(self,*,patch_rows=None,chart='Rh_reference',previous=None,**kwargs):
        if patch_rows is not None and kwargs:raise ValueError('Supply actual patch source or constructor options')
        if chart not in self.DOMAINS:raise ValueError('Original constant-axial outer chart required')
        if chart=='O2_slope' and (previous is None or previous.chart!='Rh_reference'):
            raise ValueError('Complete original Rh_reference predecessor required')
        self.patch_rows=patch_rows if patch_rows is not None else OriginalRmPressureRows(**kwargs)
        p=self.patch_rows
        for name in ('ctx','iv','J','audit','directory','family','o','t','micro','owner','reference','cover_owner'):
            setattr(self,name,getattr(p,name))
        self.chart=chart;self.previous=p if previous is None else previous
        self.incoming_coordinate='Rh' if chart=='Rh_reference' else '0'
        self.hashes,self.code_hashes=dict(p.hashes),dict(p.code_hashes)
        if previous is not None:self.hashes.update(previous.hashes);self.code_hashes.update(previous.code_hashes)
        for suffix,gate in (
            ('compliant_current_original_true_outer_C6_C5_functions_check.json','current_original_true_six_outer_charts_finite_N_C5_functions_installed'),
            ('compliant_current_original_whole_Z_Rh_functional_join_check.json','current_whole_Z_repaired_leading_Rh_reference_function_identity_checked')):
            path=self.directory/(PREFIX+suffix);raw=path.read_bytes();data=json.loads(raw)
            if not data['all_passed'] or not data[gate] or data['source_family']!=self.family:
                raise ValueError('Same original smooth Rh and outer source identities required')
            for name,digest in data['input_hashes'].items():
                if hashlib.sha256((self.directory/name).read_bytes()).hexdigest()!=digest or self.hashes.get(name,digest)!=digest:
                    raise ValueError('Changed original outer pressure source: '+name)
                self.hashes[name]=digest
            self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        path=self.directory/(PREFIX+'compliant_current_original_Rm_pressure_rows_numeric.json');raw=path.read_bytes();data=json.loads(raw)
        if not data['checks']['passed'] or data['source_family']!=self.family:raise ValueError('Accepted genuine C3 Rh pressure source required')
        for base,group,target in ((self.directory,data['input_hashes'],self.hashes),(ROOT,data['implementation_hashes'],self.code_hashes)):
            for name,digest in group.items():
                if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest or target.get(name,digest)!=digest:
                    raise ValueError('Changed genuine C3 Rh prerequisite: '+name)
                target[name]=digest
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        path=self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop=json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart'][chart]
        if self.loop['source_family']!=self.family or self.loop['chart']!=chart or self.loop['primitive_units']['B']!='original B/Pstar; E=Utheta/Pstar':
            raise ValueError('Exact original outer slow/fast derivative chart and units required')
        A=self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if A is not None and _bounds(self.audit,self.owner.logN)[0]<_bounds(self.audit,A)[1]:
            raise ValueError('Original selected N must admit local A/N')
        cls=OriginalOuterReferenceNumeric if chart=='Rh_reference' else OriginalOuterSlopeNumeric
        self.utility=object.__new__(cls)
        for name in ('ctx','iv','family','owner'):setattr(self.utility,name,getattr(self,name))
        self._histories={};self._incoming={};self._backgrounds={};self._packets={};self._errors={};self._next=None
        for name in ('original_outer_pressure_rows.py','original_outer_reference_numeric.py','original_outer_slope_numeric.py'):
            name='src/openai_ns_reconstruction/'+name;self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):
        if self.chart!='Rh_reference':raise ValueError('The next O2 axial source changes axial velocity')
        if self._next is None:self._next=OriginalOuterPressureRows(patch_rows=self.patch_rows,chart='O2_slope',previous=self)
        return self._next

    def coordinate(self,value):
        c,iv=self.ctx,self.iv
        if type(value) is not str or not c.isfinite(c.mpf(value)):raise ValueError('Finite exact decimal outer coordinate required')
        q=Fraction(value);left,right=self.DOMAINS[self.chart]
        if not left<=q<=right:raise ValueError('Original outer native domain exceeded')
        return q,c.mpf(value),iv.mpf(value)

    def geometry(self,value):
        q,p,t=self.coordinate(value);left,_=self.DOMAINS[self.chart];c,iv,ref=self.ctx,self.iv,self.reference
        terms=[(c.log(110),iv.log(110)),(10*ref.logCp,10*ref.logC),(10*ref.logPp,10*ref.logP),(p,t)]
        return dict(native=q,coordinate_point=p,coordinate=t,local_point=p-self.ctx.mpf(left.numerator)/left.denominator,
            local=t-self.iv.mpf(left.numerator)/left.denominator,atomic_log_radius_terms=terms,
            logR_point=sum((a for a,b in terms),c.mpf(0)),logR=sum((b for a,b in terms),iv.mpf(0)))

    def background(self,Z,value):
        key=Z,value
        if key not in self._backgrounds:
            geo=self.geometry(value);sides=[]
            for directed,ctx in ((False,self.ctx),(True,self.iv)):
                J=self.J;z=J.variable(ctx,ctx.mpf(Z),3);U=4*z;de=self.o.core.delta_iv if directed else self.o.core.delta
                Q=((U*(z*(1+de))-4*(1-z*z))/(1-z*z*de)).truncate(2)
                sides.append(dict(U=U,m=U,Q=Q,U_y=J.constant(ctx,0,1),Q_y=J.constant(ctx,0,1)))
            self._backgrounds[key]=*sides,geo
        return self._backgrounds[key]

    def uniform_background(self,Z,value):return self.background(Z,value)[1]

    def incoming_history_logs(self,Z):
        if Z not in self._incoming:
            raw=self.previous.history_logs(Z,self.incoming_coordinate)['history']
            if any(len(v)!=4 for v in raw.values()):raise ValueError('Genuine four-row outer history inlet required')
            self._incoming[Z]={name:list(v) for name,v in raw.items()}
        return self._incoming[Z]

    def history_logs(self,Z,value):
        key=Z,value
        if key in self._histories:return self._histories[key]
        q,_,_=self.coordinate(value);iv=self.iv;frames,groups,local=self.utility._cover(Z,q)
        hist={name:list(v) for name,v in self.incoming_history_logs(Z).items()};prefix=copy.copy({name:list(v) for name,v in hist.items()})
        steps=[];visited=[]
        for left,right,cells in groups:
            if q<=left:continue
            end=min(q,right);length=end-left;width=iv.mpf(length.numerator)/length.denominator;weights={}
            visited.extend(cells)
            for name,rate in RATES.items():
                r=iv.mpf(rate);mass=width if r==0 else -iv.expm1(-r*width)/r;out=[];sup=[]
                for k in range(4):
                    source=self.cover_owner._log_cover([p['actual_density_ranges']['N_scaled'][name][k] for p in cells])
                    driver=None if source is None else source-self.owner.logN+iv.log(mass)
                    memory=None if hist[name][k] is None else hist[name][k]-r*width
                    out.append(self.sum_logs((memory,driver)));sup.append(self.sum_logs((prefix[name][k],driver)))
                hist[name],prefix[name]=out,sup;weights[name]=dict(rate=rate,decay_log=-r*width,positive_mass=mass)
            steps.append(dict(exact_left=[left.numerator,left.denominator],exact_end=[end.numerator,end.denominator],
                physical_log_width=width,weights=weights,overlapping_Z_covers=len(cells)))
        out=dict(history=hist,prefix_history=prefix,frames=frames,groups=groups,local=local,
            visited_cells=visited if visited else local,steps=steps)
        self._histories[key]=out;return out

    def recover_errors(self,Z,hist,cells,loop):
        lp=self.o.core.lp_iv;B=[self.cover_owner._log_cover([cell['actual_source_and_primitive_range']['B'][k] for cell in cells]) for k in range(4)]
        shift=lambda v,amount:None if v is None else v+amount
        errors=dict(U=[shift(v,lp-self.owner.logN) for v in B],m=[shift(v,lp) for v in hist['m']])
        slow=loop['slow_phase_held_log_bounds']['B'];fast=loop['phase_and_first_slow_phase_log_bounds']
        sy=[shift(self.read_log(slow[key]),lp-self.owner.logN) for key in ('y1_Z0','y1_Z1')]
        fy=[shift(self.read_log(fast[key]),lp) for key in ('B_phi','B_Z_phi')]
        errors['U_y']=[self.sum_logs(pair) for pair in zip(sy,fy)]
        return errors,dict(slow_B_y_over_N=self.log_rows(sy),fast_B_phi=self.log_rows(fy),extra_hb_factor=False,extra_x_factor=False)

    def finish_errors(self,Z,errors):
        errors['U_Z'],errors['m_Z']=errors['U'][1:],errors['m'][1:]
        coeff=pressure_coefficients(self.J,self.iv,Z,self.o.core.delta_iv,2)
        errors['Q']=self.combine(coeff['Q'],errors,2);errors['Q_y']=self.combine(coeff['Q_y'],errors,1)
        return errors

    def incoming_error_data(self,Z):
        if isinstance(self.previous,OriginalOuterPressureRows):return self.previous.error_data(Z,self.incoming_coordinate)
        p=self.patch_rows;raw=p.history_logs(Z,'Rh')
        errors,components=self.recover_errors(Z,raw['history'],raw['local'],p.loop)
        return self.finish_errors(Z,errors),components

    def error_data(self,Z,value,prefix=False):
        key=Z,value,prefix
        if key in self._errors:return self._errors[key]
        q,_,_=self.coordinate(value);left,_=self.DOMAINS[self.chart]
        h=self.history_logs(Z,value)
        errors,components=self.recover_errors(Z,h['prefix_history' if prefix else 'history'],
            h['visited_cells' if prefix else 'local'],self.loop)
        if prefix or q==left:
            # Keep the new chart's own slow/fast enclosure even at its inlet.
            # Preserve the predecessor budget by union; caps are not values.
            old,_=self.incoming_error_data(Z)
            for name in ('U','U_y'):
                errors[name]=[self.sum_logs((a,b)) for a,b in zip(errors[name],old[name])]
            components['additional_incoming_Uy_log_bounds']=self.log_rows(old['U_y'])
            components['new_chart_local_slow_fast_bounds_not_identified_with_predecessor']=True
        out=self.finish_errors(Z,errors),components
        self._errors[key]=out;return out

    def evaluate(self,Z,value):
        key=Z,value
        if key in self._packets:return self._packets[key]
        if type(Z) is not str or not self.ctx.isfinite(self.ctx.mpf(Z)) or not -1<self.ctx.mpf(Z)<1:
            raise ValueError('Interior exact decimal Z required')
        p,t,geo=self.background(Z,value);h=self.history_logs(Z,value);errors,components=self.error_data(Z,value)
        fields={name:dict(signed_background_ordinary_Z=[self.o._summary(p[name][k]*math.factorial(k),t[name][k]*math.factorial(k)) for k in range(order+1)],
            finite_N_absolute_error_log_bounds=self.log_rows(errors[name][:order+1]),finite_N_log_bound_is_not_selected_signed_value=True)
            for name,order in (('U',2),('m',3),('Q',2),('U_y',1),('Q_y',1))}
        out=dict(Z=Z,chart=self.chart,native_coordinate=value,source_family=self.family,selected_integer=copy.deepcopy(self.owner.integer),
            physical_pressure_rows=fields,finite_N_history_error_logs={name:self.log_rows(v) for name,v in h['history'].items()},
            finite_N_incoming_history_error_logs={name:self.log_rows(v) for name,v in self.incoming_history_logs(Z).items()},
            finite_N_prefix_history_error_logs={name:self.log_rows(v) for name,v in h['prefix_history'].items()},
            finite_N_axial_radial_error_components=components,
            physical_transport_steps=self.reference.previous._serializable_proof(h['steps']),
            geometry=dict(atomic_log_radius_terms=[self.o._summary(a,b) for a,b in geo['atomic_log_radius_terms']],
                actual_native_length=self.o._summary(geo['local_point'],geo['local']),physical_radius_materialized=False,
                native_coordinate_is_physical_log_offset=True,exact_start=geo['native']==self.DOMAINS[self.chart][0]),
            source_cover=dict(exact_Z_cells=h['frames'],unique_native_cells=len(h['groups']),local_source_cells=len(h['local']),
                generic_derivative_chart=self.chart,native_domain=[str(v) for v in self.DOMAINS[self.chart]],
                all_upstream_histories_retained=True,overlapping_Z_cells_do_not_duplicate_width=True),
            contract=dict(original_U_equal_M_equal_4Z_identity=True,genuine_analytic_mean_Z3_zero=True,
                rational_Q_Z2_from_actual_source=True,leading_Uy_Qy_zero_does_not_reset_errors=True,
                full_C3_incoming_and_local_history_rows=True,native_physical_width_once=True,
                incoming_N_factor_not_reapplied=True,new_driver_N_factor_once=True,
                original_Pstar_once=True,slow_fast_phase_derivative_units=True,
                shared_boundary_local_errors_reused=True,point_derivatives_not_error_cap_derivatives=True),
            limits=dict(original_outer_pressure_UQ_Z2_mean_Z3_UyQy_Z1_with_complete_absolute_errors=True,
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

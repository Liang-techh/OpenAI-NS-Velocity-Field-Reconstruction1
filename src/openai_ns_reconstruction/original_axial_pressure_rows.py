"""Actual O2 axial pressure from the flat cutoff and nonzero positive mean memory."""
import copy
from fractions import Fraction
import hashlib
import json
import math

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_outer_pressure_rows import OriginalOuterPressureRows
from .original_outer_axial_numeric import OriginalOuterAxialNumeric
from .original_axial_rate_memory import OriginalAxialRateMemory
from .original_Rm_finite_N_numeric import RATES
from .original_Rm_patch_numeric import encode


class OriginalAxialPressureRows(OriginalOuterPressureRows):
    DOMAINS={'O2_axial':(Fraction(0),Fraction(1))}

    def __init__(self,*,previous=None,**kwargs):
        if previous is not None and kwargs:raise ValueError('Supply original slope pressure or constructor options')
        self.previous=previous if previous is not None else OriginalOuterPressureRows(**kwargs).next_field()
        if self.previous.chart!='O2_slope':raise ValueError('Complete original O2 slope pressure predecessor required')
        p=self.previous
        for name in ('ctx','iv','J','audit','directory','family','o','t','micro','owner','reference','cover_owner','patch_rows'):
            setattr(self,name,getattr(p,name))
        self.chart='O2_axial';self.incoming_coordinate='1'
        self.hashes,self.code_hashes=dict(p.hashes),dict(p.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_outer_m15_numeric.json');raw=path.read_bytes();data=json.loads(raw)
        if not data['checks']['passed'] or data['source_family']!=self.family:raise ValueError('Accepted original slope C3 pressure state required')
        for base,group,target in ((self.directory,data['input_hashes'],self.hashes),(ROOT,data['implementation_hashes'],self.code_hashes)):
            for name,digest in group.items():
                if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest or target.get(name,digest)!=digest:
                    raise ValueError('Changed original axial pressure source: '+name)
                target[name]=digest
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        path=self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop=json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart'][self.chart]
        if self.loop['chart']!=self.chart or self.loop['source_family']!=self.family or self.loop['primitive_units']['B']!='original B/Pstar; E=Utheta/Pstar':
            raise ValueError('Actual physical-y O2 axial loop derivatives required')
        A=self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if A is not None and _bounds(self.audit,self.owner.logN)[0]<_bounds(self.audit,A)[1]:raise ValueError('Original N must admit A/N')
        self.utility=object.__new__(OriginalOuterAxialNumeric)
        for name in ('ctx','iv','family','owner'):setattr(self.utility,name,getattr(self,name))
        self.memory=OriginalAxialRateMemory(self.ctx,self.iv,self.J)
        self._histories={};self._incoming={};self._backgrounds={};self._packets={};self._errors={};self._scalars={};self._uniforms={}
        for name in ('original_axial_pressure_rows.py','original_axial_memory.py','original_axial_rate_memory.py','original_outer_axial_numeric.py'):
            name='src/openai_ns_reconstruction/'+name;self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):raise ValueError('Continue the actual O2 buffer with the nonzero axial mean endpoint')

    def geometry(self,value):
        q,sp,st=self.coordinate(value);c,iv,ref=self.ctx,self.iv,self.reference
        dp,dt=c.expm1(40*sp),iv.expm1(40*st)
        terms=[(c.log(110),iv.log(110)),(10*ref.logCp,10*ref.logC),(10*ref.logPp,10*ref.logP),
               (c.mpf(1),iv.mpf(1)),(dp,dt)]
        return dict(native=q,coordinate_point=sp,coordinate=st,local_point=dp,local=dt,
            physical_jacobian_point=40*c.exp(40*sp),physical_jacobian=40*iv.exp(40*st),
            atomic_log_radius_terms=terms,logR_point=sum((a for a,b in terms),c.mpf(0)),
            logR=sum((b for a,b in terms),iv.mpf(0)))

    def scalars(self,value):
        if value not in self._scalars:
            q,sp,st=self.coordinate(value);memory=self.memory.evaluate(value,1);sides=[]
            for directed,ctx,s in ((False,self.ctx,sp),(True,self.iv,st)):
                if q==0:beta,by=ctx.mpf(1),ctx.mpf(0)
                elif q==1:beta,by=ctx.mpf(0),ctx.mpf(0)
                else:
                    logjet=self.memory._log_beta(ctx,self.J.variable(ctx,1-s,4));jet=logjet.exp()
                    beta=jet[0];by=-jet[1]/(40*ctx.exp(40*s))
                sides.append(dict(beta=beta,beta_y=by,mean=memory['directed' if directed else 'point']))
            self._scalars[value]=dict(point=sides[0],truth=sides[1],memory=memory)
        return self._scalars[value]

    def fields(self,Z,ctx,beta,mean,beta_y):
        z=self.J.variable(ctx,ctx.mpf(Z),3);de=self.o.core.delta_iv if ctx is self.iv else self.o.core.delta
        A=1-z*z*de
        # Exact analytic cancellation in the original mean recovery:
        # Q=8 Z² beta/A -4 mean, with beta/mean independent of Z.
        U=4*z*beta;M=4*z*mean;Uy=4*z*beta_y
        Q=(z*z*(8*beta)/A-4*mean).truncate(2)
        Qy=(z*z*(8*beta_y)/A-4*(beta-mean)).truncate(2)
        return dict(U=U,m=M,Q=Q,U_y=Uy,Q_y=Qy)

    def background(self,Z,value):
        key=Z,value
        if key not in self._backgrounds:
            scalars=self.scalars(value);sides=[]
            for label,ctx in (('point',self.ctx),('truth',self.iv)):
                v=scalars[label];sides.append(self.fields(Z,ctx,v['beta'],v['mean'],v['beta_y']))
            self._backgrounds[key]=*sides,self.geometry(value)
        return self._backgrounds[key]

    def uniform_background(self,Z,value):
        key=Z,value
        if key not in self._uniforms:
            v=self.scalars(value)['truth'];lo=lambda row:max(0,_bounds(self.audit,row)[0])
            # beta and its positive relaxation mean decrease from1. Only
            # U/Q are used by the volume bound; no false uniform Uy=0 claim.
            beta=self.iv.mpf([lo(v['beta']),1]);mean=self.iv.mpf([lo(v['mean']),1])
            z=self.J.variable(self.iv,self.iv.mpf(Z),3);A=1-z*z*self.o.core.delta_iv
            self._uniforms[key]=dict(U=4*z*beta,m=4*z*mean,Q=(z*z*(8*beta)/A-4*mean).truncate(2))
        return self._uniforms[key]

    def history_logs(self,Z,value):
        key=Z,value
        if key in self._histories:return self._histories[key]
        q,_,_=self.coordinate(value);iv=self.iv;frames,groups,local=self.utility._cover(Z,q)
        hist={name:list(v) for name,v in self.incoming_history_logs(Z).items()};prefix={name:list(v) for name,v in hist.items()}
        steps=[];visited=[];rational=lambda f:iv.mpf(f.numerator)/f.denominator
        for left,right,cells in groups:
            if q<=left:continue
            end=min(q,right);width=iv.exp(40*rational(left))*iv.expm1(40*rational(end-left));visited.extend(cells);weights={}
            for name,rate in RATES.items():
                r=iv.mpf(rate);mass=width if r==0 else -iv.expm1(-r*width)/r;out=[];sup=[]
                for k in range(4):
                    source=self.cover_owner._log_cover([cell['actual_density_ranges']['N_scaled'][name][k] for cell in cells])
                    driver=None if source is None else source-self.owner.logN+iv.log(mass)
                    memory=None if hist[name][k] is None else hist[name][k]-r*width
                    out.append(self.sum_logs((memory,driver)));sup.append(self.sum_logs((prefix[name][k],driver)))
                hist[name],prefix[name]=out,sup;weights[name]=dict(rate=rate,decay_log=-r*width,positive_mass=mass)
            steps.append(dict(exact_left=[left.numerator,left.denominator],exact_end=[end.numerator,end.denominator],
                physical_log_width=width,weights=weights,overlapping_Z_covers=len(cells),
                exact_measure='dd=40exp(40s)ds once'))
        result=dict(history=hist,prefix_history=prefix,frames=frames,groups=groups,local=local,
            visited_cells=visited if visited else local,steps=steps)
        self._histories[key]=result;return result

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
            physical_pressure_rows=fields,leading_scalars={name:self.o._summary(v['point'][name],v['truth'][name]) for name in ('beta','beta_y','mean')},
            positive_mean_memory_receipt=encode(v['memory']['proof']),
            finite_N_history_error_logs={name:self.log_rows(values) for name,values in h['history'].items()},
            finite_N_incoming_history_error_logs={name:self.log_rows(values) for name,values in self.incoming_history_logs(Z).items()},
            finite_N_prefix_history_error_logs={name:self.log_rows(values) for name,values in h['prefix_history'].items()},
            finite_N_axial_radial_error_components=components,physical_transport_steps=encode(h['steps']),
            geometry=dict(atomic_log_radius_terms=[self.o._summary(a,b) for a,b in geo['atomic_log_radius_terms']],
                physical_log_distance=self.o._summary(geo['local_point'],geo['local']),
                physical_Jacobian=self.o._summary(geo['physical_jacobian_point'],geo['physical_jacobian']),
                native_coordinate_is_physical_log_offset=False,physical_radius_materialized=False,exact_start=geo['native']==0),
            source_cover=dict(exact_Z_cells=h['frames'],unique_native_cells=len(h['groups']),local_source_cells=len(h['local']),
                generic_derivative_chart='O2_axial',native_domain=['0','1'],all_upstream_histories_retained=True),
            contract=dict(actual_U_4Z_beta_and_M_4Z_positive_memory=True,genuine_analytic_mean_Z3_zero=True,
                exact_rational_QZZ_and_actual_QyZ=True,cutoff_derivative_divided_by_physical_Jacobian=True,
                local_terminal_zero_does_not_reset_mean=True,terminal_Qy_equals_four_mean=True,
                full_C3_incoming_and_local_histories=True,exponential_physical_width_once=True,
                slow_fast_bounds_in_physical_y_units=True,new_chart_and_predecessor_endpoint_budgets=True,
                no_extra_N_epsilon_hb_or_native_Jacobian_on_stored_By=True),
            limits=dict(actual_O2_axial_pressure_UQ_Z2_mean_Z3_UyQy_Z1_with_complete_absolute_errors=True,
                changing_axial_m15_integrated=False,resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

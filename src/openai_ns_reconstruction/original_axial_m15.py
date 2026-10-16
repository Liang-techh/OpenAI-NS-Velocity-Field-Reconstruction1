"""Original changing O2 axial pressure moment with positive own-mean kernels."""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_outer_m15 import OriginalOuterM15
from .original_axial_pressure_rows import OriginalAxialPressureRows
from .original_axial_rate_memory import OriginalAxialRateMemory
from .original_Rm_patch_numeric import encode


class OriginalAxialM15(OriginalOuterM15):
    def __init__(self,*,rows=None,**kwargs):
        if rows is not None and kwargs:raise ValueError('Supply axial pressure rows or constructor options')
        self.rows=rows if rows is not None else OriginalAxialPressureRows(**kwargs)
        self.chart=self.rows.chart
        for name in ('t','o','ctx','iv','J','audit','directory','family'):setattr(self,name,getattr(self.rows,name))
        self.hashes,self.code_hashes=dict(self.rows.hashes),dict(self.rows.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_outer_m15_numeric.json')
        data=json.loads(path.read_bytes())
        self.inlets={Fraction(p['Z']):p for p in data['samples'] if p['chart']=='O2_slope' and p['native_coordinate']=='1'}
        if not self.inlets or not data['checks']['passed']:raise ValueError('Accepted complete slope pressure inlet required')
        self._slope_engine=None;self._packets={};self._kernels={}
        self.rate2=OriginalAxialRateMemory(self.ctx,self.iv,self.J,kernel_rate=2)
        name='src/openai_ns_reconstruction/original_axial_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):raise ValueError('Continue O2 buffer with the complete nonzero own-mean endpoint')

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            if self._slope_engine is None:
                reference=OriginalOuterM15(rows=self.rows.previous.previous)
                self._slope_engine=OriginalOuterM15(rows=self.rows.previous,previous=reference)
            self.inlets[key]=self._slope_engine.evaluate(Z,'1')
            self.hashes.update(self._slope_engine.hashes);self.code_hashes.update(self._slope_engine.code_hashes)
        return self.inlets[key]

    def coefficients(self,Z,ctx):
        z=self.J.variable(ctx,ctx.mpf(Z),1)
        de=self.o.core.delta_iv if ctx is self.iv else self.o.core.delta
        A=1-z*z*de;D=1-z*z
        F=dict(b=4*z*z/A**2,bm=-16*z*z/A,m2=self.J.constant(ctx,4,1),m=-2/A)
        P=z**4*(de*de+de)-z*z*(2*de+3)+3
        G=dict(b2=16*z*z*P/A**3,bm=16*z*z/A,m2=self.J.constant(ctx,-4,1),
            # Preserve the genuine delta-sized source; subtracting its
            # unsimplified O(1) pieces loses it at original parameters.
            b=-4*z*z*D*de/A**3,m=2/A)
        return F,G,dict(by=-8*z*z/A,b=self.J.constant(ctx,4,1),m=self.J.constant(ctx,-4,1))

    def kernels(self,value):
        if value not in self._kernels:
            self._kernels[value]=dict(beta=self.rate2.evaluate(value,1),beta2=self.rate2.evaluate(value,2),
                mean2=self.rows.memory.squared_mean_integral(value))
        return self._kernels[value]

    def finite_errors(self,Z,value,packet):
        rows=self.rows;_,base,geo=rows.background(Z,value);_,lower_base,lower_geo=rows.background(Z,'0')
        incoming,_=rows.error_data(Z,'0');upper=self.correction_densities(Z,packet,base)
        lower=self.correction_densities(Z,self.packet_for_errors(incoming),lower_base)
        prefix,_=rows.error_data(Z,value,prefix=True)
        volume=self.correction_densities(Z,self.packet_for_errors(prefix),rows.uniform_background(Z,value))['volume']
        result={}
        for name,power in (('endpoint',2),('viscous',1)):
            result[name]=[rows.sum_logs((None if upper[name][k] is None else upper[name][k]+power*geo['logR'],
                None if lower[name][k] is None else lower[name][k]+power*lower_geo['logR'])) for k in (0,1)]
        # integral R^2 dy=Rupper^2*(1-exp(-2D))/2; no huge radius.
        mass=-self.iv.expm1(-2*geo['local'])/2
        result['volume']=[None if v is None else v+2*geo['logR']+self.iv.log(mass) for v in volume]
        return result,dict(endpoint=rows.log_rows(lower['endpoint']),viscous=rows.log_rows(lower['viscous']),
            prefix_error_rows={name:rows.log_rows(v) for name,v in prefix.items() if name in ('U','m','Q','U_y','Q_y')},
            normalized_positive_volume_mass=self.o._summary(-self.ctx.expm1(-2*geo['local_point'])/2,mass),
            normalization='Rupper^2; physical dd once')

    def evaluate(self,Z,value):
        key=Z,value
        if key in self._packets:return self._packets[key]
        packet=self.rows.evaluate(Z,value);incoming=self.inlet(Z);count=176
        if len(incoming['m15_C0_C1_terms'])!=count or incoming['actual_selected_integer']!=packet['selected_integer']:
            raise ValueError('Complete same-N 176-term slope pressure inlet required')
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;geo=self.rows.geometry(value);empty=geo['native']==0
        zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1)
        cp,ct=self.coefficients(Z,c),self.coefficients(Z,iv);sc=self.rows.scalars(value)
        p,t=sc['point'],sc['truth'];terms=copy.deepcopy(incoming['m15_C0_C1_terms'])
        lower_geo=self.rows.geometry('0')
        def atoms(geometry,power):return [o._summary(power*a,power*b) for a,b in geometry['atomic_log_radius_terms']]
        lower,upper,viscous=atoms(lower_geo,2),atoms(geo,2),atoms(geo,1)
        def append(name,pv,tv,scale,**binding):
            terms.append(dict(name='O2_axial_'+name,
                coefficient_C0_C1=o._jet_summary(zero if empty else pv,zeroi if empty else tv),
                log_scale_terms=[] if empty else copy.deepcopy(scale),
                actual_function_binding=dict(chart=self.chart,incoming_terms_retained_once=count,
                    original_changing_axial_pressure_source=True,**binding)))
        for mode,scalar in (('b','beta'),('bm',None),('m2',None),('m','mean')):
            sp=p['beta']*p['mean'] if mode=='bm' else p['mean']**2 if mode=='m2' else p[scalar]
            st=t['beta']*t['mean'] if mode=='bm' else t['mean']**2 if mode=='m2' else t[scalar]
            append('F_upper_'+mode,cp[0][mode]*sp,ct[0][mode]*st,upper,sector='Rupper^2*F')
        append('F_lower',-sum(cp[0].values(),zero),-sum(ct[0].values(),zeroi),lower,sector='-Rlower^2*F(0)')
        kernels=None if empty else self.kernels(value)
        def source_integral(name,mode,coefficientp,coefficientt):
            if empty:append(name,zero,zeroi,[],sector='exact empty beta integral');return
            d=kernels[mode]['distribution'];sp=d['normalized_source_point'];st=d['normalized_source_integral']
            extra=[o._summary(d['anchor'],iv.mpf(d['anchor'])),o._summary(c.log(d['scale']),iv.log(iv.mpf(d['scale'])))]
            append(name,coefficientp*sp,coefficientt*st,upper+extra,
                sector='positive original integral exp(2d)*beta^p dd',source_kernel=mode,
                homogeneous_exponential_not_subtracted=True,source_function='original_flat_axial_beta')
        source_integral('G_b2','beta2',cp[1]['b2'],ct[1]['b2'])
        append('G_bm_upper',cp[1]['bm']*(p['mean']**2/2),ct[1]['bm']*(t['mean']**2/2),upper,
            exact_reduction='I_bm=(exp(2D)*m(D)^2-1)/2')
        append('G_bm_lower',-cp[1]['bm']/2,-ct[1]['bm']/2,lower,sector='same exact I_bm lower endpoint')
        if empty:append('G_m2',zero,zeroi,[],sector='exact empty own-mean-square integral')
        else:
            k=kernels['mean2']
            # Keep 2*anchor and 2*log(scale) separate from 2D to avoid
            # cancellation when constructing the physical exponent atoms.
            d=self.rows.memory.evaluate(value,1)['distribution']
            extra=[o._summary(2*d['anchor'],2*iv.mpf(d['anchor'])),o._summary(2*c.log(d['scale']),2*iv.log(iv.mpf(d['scale'])))]
            append('G_m2',cp[1]['m2']*k['point'],ct[1]['m2']*k['truth'],upper+extra,
                sector='positive squared own-mean integral',source_function='original_positive_own_mean',
                exact_positive_future_CDF=True)
        source_integral('G_b','beta',cp[1]['b'],ct[1]['b'])
        append('G_m_upper',cp[1]['m']*p['mean'],ct[1]['m']*t['mean'],upper,
            exact_reduction='I_m=exp(2D)*m(D)-1-I_b')
        append('G_m_lower',-cp[1]['m'],-ct[1]['m'],lower,sector='same exact I_m lower endpoint')
        source_integral('G_m_minus_Ib','beta',-cp[1]['m'],-ct[1]['m'])
        for mode,scalar in (('by','beta_y'),('b','beta'),('m','mean')):
            append('minus_RQy_upper_'+mode,cp[2][mode]*p[scalar],ct[2][mode]*t[scalar],viscous,
                sector='-Rupper*Qy',terminal_mean_not_zero=True,lower_leading_Qy_exact_zero=True)
        if empty:finite={name:[None,None] for name in ('endpoint','volume','viscous')};budget=None
        else:finite,budget=self.finite_errors(Z,value,packet)
        for name in ('endpoint','volume','viscous'):
            bound=self.rows.sum_logs(finite[name]);truth=zeroi;scale=[]
            if bound is not None:
                top=iv.mpf(_bounds(self.audit,bound)[1]);atom=o._summary(_bounds(c,top)[1],top)
                atom.update(positive_error_bound_atom_only=True,source_value_selected=False);scale.append(atom)
                truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
            append('finite_N_'+name+'_error',zero,truth,scale,error_enclosure_only=True,
                both_actual_endpoints=True,full_prefix_C3_history_and_local_slow_fast_bounds=True)
        proofs=None if empty else dict(beta_rate2=encode(kernels['beta']['proof']),
            beta2_rate2=encode(kernels['beta2']['proof']),mean2=encode(kernels['mean2']['proof']))
        out=dict(Z=Z,chart=self.chart,native_coordinate=value,source_family=self.family,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']),m15_C0_C1_terms=terms,
            incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']),positive_kernel_proofs=proofs,
            finite_N_error_logs_in_physical_units={name:self.rows.log_rows(v) for name,v in finite.items()},
            actual_lower_finite_N_endpoint_budget=budget,
            leading_upper_Qy=copy.deepcopy(packet['physical_pressure_rows']['Q_y']['signed_background_ordinary_Z']),
            contract=dict(complete_176_term_inlet_once=True,original_changing_U_M_pressure=True,
                original_half_S1_once=True,exact_I_bm_and_I_m_reductions=True,
                positive_source_aware_beta_and_mean_square_integrals=True,delta_sized_G_b_retained=True,
                nonzero_terminal_RQy_retained=True,both_finite_N_endpoints=True,
                full_prefix_C3_errors_not_final_decayed_history=True,physical_measure_once=True,
                no_huge_radius_materialized=True,all_quadrature_errors_in_coefficients=True),
            limits=dict(original_O2_axial_signed_cumulative_m15_with_absolute_errors=True,
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                independent_reserved_Rp_seam=False,global_five_moment_numeric_join=False,
                common_global_N=False,exact_heat_exterior=False,stress_cone=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

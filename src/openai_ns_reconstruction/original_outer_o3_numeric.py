"""Actual original O3 transition/power through Rc, retaining tiny mu sectors."""
from fractions import Fraction
import hashlib
import json
import math

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_finite_N_numeric import RATES,transport_step
from .original_Rm_patch_numeric import encode
from .original_Rm_patch_velocity import OriginalRmPatchVelocity
from .original_outer_buffer_numeric import OriginalOuterBufferNumeric,OriginalOuterBufferVelocity
from .original_o3_scalar import OriginalO3Scalar


class OriginalOuterO3Numeric:
    lengths={'O3_transition':1,'O3_power':2}

    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,buffer=None):
        self.buffer=buffer or OriginalOuterBufferNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        a=self.buffer;self.patch=a.patch;self.ctx,self.iv,self.J=a.ctx,a.iv,a.J
        self.family=a.family;self.owner=a.owner;self.directory=a.directory
        self.hashes=dict(a.hashes);self.code_hashes=dict(a.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_outer_buffer_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:raise ValueError('Accepted same-original buffer endpoint required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed buffer implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._inlets={p['Z']:p for p in accepted['samples'] if p['native_buffer']=='11'};self._packets={}
        c,iv=self.ctx,self.iv
        self.logmu=c.log(c.mpf('.001'))-4*self.patch.leading.defect.logPp
        self.logmui=iv.log(iv.mpf('.001'))-4*self.patch.leading.defect.logP
        self.mu,self.mui=c.exp(self.logmu),iv.exp(self.logmui)
        if not _bounds(c,self.logmui)[1]<-1000 or not _bounds(c,-60*self.logmui)[0]>2+_bounds(c,iv.log(2))[1]:
            raise ValueError('Original tiny positive mu and Rc reservation required')
        self.scalar=OriginalO3Scalar(c,iv,self.mu,self.mui)
        for name in ('src/openai_ns_reconstruction/original_o3_scalar.py','src/openai_ns_reconstruction/original_outer_o3_numeric.py'):
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _rows(self,p,t):return self.buffer._rows(p,t)
    def _jet(self,rows,directed):return self.buffer._jet(rows,directed)

    def _cover(self,Z,q,chart):
        c,iv=self.ctx,self.iv;zl,zh=_bounds(c,iv.mpf(Z));frames=[];groups={};local=[]
        for frame in self.owner.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            if frame['source_identity']!=self.family or not frame['uniform_for_all_integer_N_ge_N0'] or frame['epsilon_domain']!='[0,1/N0]':
                raise ValueError('Same-original all-N O3 ranges required')
            frames.append(frame['exact_Z_cell'])
            source=next(p for p in frame['actual_17_chart_C3_transports'] if p['chart']==chart)
            for cell in source['cells']:
                left,right=Fraction(*cell['exact_left']),Fraction(*cell['exact_right'])
                groups.setdefault((left,right),[]).append(cell)
                if left<=q<=right:local.append(cell)
        ordered=sorted(groups)
        if (not local or not ordered or ordered[0][0]!=0 or ordered[-1][1]!=self.lengths[chart]
            or any(a[1]!=b[0] for a,b in zip(ordered,ordered[1:]))):raise ValueError('Complete original O3 source cover required')
        return frames,[(a,b,groups[a,b]) for a,b in ordered],local

    def _leading(self,Z,x,chart,buffer_inlet,scalar,ctx):
        c,iv=self.ctx,self.iv;directed=ctx is iv
        # The mu=0 reference is only an algebraic frame. The full original
        # field includes every separately enclosed mu-scaled correction below.
        # Exact decimal strings are retained by Fraction in the traversal;
        # directed arithmetic constructs the reference coordinate directly.
        t=ctx.mpf(x)+(1 if chart=='O3_power' else 0)
        axial=dict(original_slope_scalar_inlet=buffer_inlet['original_slope_scalar_inlet'],
            axial_profiles=buffer_inlet['original_axial_profile_inlet'])
        fields,sectors,Q=self.buffer._leading(Z,ctx.mpf(11)+t,axial,ctx)
        E0=self._jet(buffer_inlet['generic_leading_fields']['E'],directed)
        s=scalar['directed' if directed else 'point'];mu=self.mui if directed else self.mu
        correction=dict(E=fields['E']*(mu*s['U1']),h=E0*(ctx.exp(-ctx.mpf('1.5')*t)*mu*s['C1']),
            e=-(E0*E0)*(ctx.exp(-t)*mu*s['C0']/2),p=E0*E0*(mu*s['Cm']/2))
        for n,v in correction.items():fields[n]+=v;sectors['mu_correction_'+n]=v
        return fields,sectors,Q

    def evaluate(self,Z,chart,x):
        c,iv,J=self.ctx,self.iv,self.J
        if chart not in self.lengths or type(x) is not str or not c.isfinite(c.mpf(x)):raise ValueError('Original O3 chart and finite decimal coordinate required')
        q=Fraction(x)
        if not 0<=q<=self.lengths[chart]:raise ValueError('Original O3 native domain exceeded')
        key=(Z,chart,x)
        if key in self._packets:return self._packets[key]
        if Z not in self._inlets:self._inlets[Z]=self.buffer.evaluate(Z,'11')
        buffer_inlet=self._inlets[Z]
        inlet=buffer_inlet if chart=='O3_transition' else self.evaluate(Z,'O3_transition','1')
        scalar=self.scalar.transition(x) if chart=='O3_transition' else self.scalar.power(x)
        frames,groups,local=self._cover(Z,q,chart);proof=[];steps=[]
        history={n:self._jet(v,True) for n,v in inlet['finite_N_history_errors'].items()}
        for left,right,cells in groups:
            if q<=left:continue
            end=min(q,right);width=iv.mpf((end-left).numerator)/(end-left).denominator;weights={}
            for name,rate in RATES.items():
                coefficients=[]
                for k in range(3):
                    log=self.patch._log_cover([v['actual_density_ranges']['N_scaled'][name][k] for v in cells])
                    cap=self.patch._small(None if log is None else log-self.owner.logN,proof,chart+'_density_'+name+'_Z'+str(k))/math.factorial(k)
                    coefficients.append(self.patch.leading.reference.leading._symmetric(cap))
                history[name],decay,mass=transport_step(iv,history[name],J(iv,coefficients),width,rate)
                weights[name]=dict(rate=rate,decay=decay,mass=mass)
            steps.append(dict(exact_left=[left.numerator,left.denominator],exact_end=[end.numerator,end.denominator],
                physical_log_width=width,own_weights=weights,overlapping_Z_covers=len(cells),physical_measure='dx once',
                incoming_already_contains_one_over_N=True,new_driver_divided_by_N_once=True))
        capsA=[];capsB=[]
        for k in range(3):
            for label,target in (('A',capsA),('B',capsB)):
                log=self.patch._log_cover([v['actual_source_and_primitive_range'][label][k] for v in local])
                target.append(self.patch._small(None if log is None else log-self.owner.logN,proof,chart+'_local_'+label+'_Z'+str(k))/math.factorial(k))
        if _bounds(c,capsA[0])[1]>=c.mpf('.25'):raise ArithmeticError('Small original O3 exponent required')
        point,sectorsp,Qp=self._leading(Z,x,chart,buffer_inlet,scalar,c)
        truth,sectorst,Qt=self._leading(Z,x,chart,buffer_inlet,scalar,iv)
        p0p=self._jet(inlet['independent_P0_over_Pstar_squared'],False).truncate(2)
        p0t=self._jet(inlet['independent_P0_over_Pstar_squared'],True).truncate(2)
        point['pressure']=point['p']+p0p;truth['pressure']=truth['p']+p0t
        sectorsp['P0_plus_p_constant']=p0p+sectorsp['p_constant'];sectorst['P0_plus_p_constant']=p0t+sectorst['p_constant']
        symmetric=self.patch.leading.reference.leading._symmetric
        dE=truth['E']*J(iv,[symmetric(v) for v in (2*capsA[0],2*capsA[1],2*capsA[2]+capsA[1]*capsA[1])])
        dV=J(iv,[symmetric(v) for v in capsB])
        # The transition endpoint shares its exact flat local source with power.
        if q==0 and chart=='O3_transition':dE,dV=(self._jet(inlet['finite_N_local_errors'][n],True) for n in ('E','V'))
        modified={n:truth[n]+history[n] for n in RATES}
        modified.update(E=truth['E']+dE,V=truth['V']+dV,pressure=truth['pressure']+history['p'])
        if q==0 and chart=='O3_transition':
            for n in ('E','V'):modified[n]=self._jet(inlet['generic_modified_fields'][n],True)
        Pp,Pt=c.exp(self.patch.leading.defect.logPp),iv.exp(self.patch.leading.defect.logP)
        z=J.variable(iv,iv.mpf(Z),2);dm=history['m'];dmz=J(iv,[dm[1],2*dm[2]])
        de=iv.mpf(self.patch.leading.reference.leading.core.delta_iv)
        dQ=((z*dV*2-z*dm*(1-de)-(1-z*z)*dmz)/(1-z*z*de)*Pt).truncate(1)
        sp=dict(Uz=point['V']*Pp,Utheta_over_Pstar=point['E'],Mz_over_R=point['m']*Pp,Ur_over_sqrt_R_over2=Qp)
        st=dict(Uz=modified['V']*Pt,Utheta_over_Pstar=modified['E'],Mz_over_R=modified['m']*Pt,Ur_over_sqrt_R_over2=Qt+dQ)
        zero=J.constant(c,0,2);ref=self.patch.leading.reference;offset=11 if chart=='O3_transition' else 12
        packet=dict(Z=Z,chart='actual_original_'+chart,native_O3=x,native_O3_chart=chart,native_x='exp('+str(6+offset)+'+exp40+x)',source_family=self.family,
            selected_integer=inlet['selected_integer'],generic_leading_fields={n:self._rows(point[n],truth[n]) for n in point},
            generic_modified_fields={n:self._rows(point[n],modified[n]) for n in modified},physical_similarity={n:self._rows(sp[n],st[n]) for n in sp},
            leading_small_sectors={n:self._rows(sectorsp[n],sectorst[n]) for n in sectorsp},leading_Q_exact_Z_independent=self._rows(Qp,Qt),
            normalized_mu_profiles={n:ref._summary(scalar['point'][n],scalar['directed'][n]) for n in scalar['point']},
            scalar_method={k:v for k,v in scalar.items() if k not in ('point','directed')},
            original_mu=dict(log_mu=ref._summary(self.logmu,self.logmui),mu=ref._summary(self.mu,self.mui),recipe='mu=exp(-4logPstar)/1000',
                strictly_positive=True,kept_separate_from_order_one_sums=True),
            finite_N_history_errors={n:self._rows(zero,v) for n,v in history.items()},finite_N_local_errors=dict(E=self._rows(zero,dE),V=self._rows(zero,dV)),
            finite_N_nonzero_incoming_errors=inlet['finite_N_history_errors'],physical_transport_steps=encode(steps),
            independent_P0_over_Pstar_squared=inlet['independent_P0_over_Pstar_squared'],positive_order_n1_own_history=inlet['positive_order_n1_own_history'],
            coordinate=dict(global_y_log_terms=[ref._summary(10*ref.logCp,10*ref.logC),
                ref._summary(10*self.patch.leading.defect.logPp,10*self.patch.leading.defect.logP),ref._summary(c.exp(40),iv.exp(40)),
                ref._summary(c.mpf(offset),iv.mpf(offset)),ref._summary(c.mpf(x),iv.mpf(x))]),
            source_cover=dict(exact_Z_cells=frames,unique_native_cells=len(groups),local_source_cells=len(local),original_all_N_ge_N0=True,phase_not_selected=True),
            small_error_proofs=proof,contract=dict(original_N_retained=True,pressure_memory_exactly_one=True,P0_added_once=True,
                original_mu_correction_not_dropped=True,mu_zero_reference_is_only_an_algebraic_frame=True,
                actual_finite_N_history_inlets_not_reset=True,quiet_power_local_zero_does_not_zero_history=True,
                local_E_V_are_source_values_not_additively_transported=True,leading_Q_from_actual_mass_identity=True),
            limits=dict(physical_Cartesian_assembly=False,resolved_signed_finite_N_centers=False,total_radial_or_fast_phase_derivatives=False,
                beyond_Rc=False,reduced_m13_m15_integrals_continued=False,relative_pressure_cancellation_accuracy=False,
                admissible_stress_cone_margin=False,global_n1_high_order_joins=False,global_N_admission=False,global_field=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=packet;return packet


class OriginalOuterO3Velocity(OriginalOuterBufferVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.numeric=OriginalOuterO3Numeric(buffer=self.numeric)
        self.hashes.update(self.numeric.hashes);self.code_hashes.update(self.numeric.code_hashes)
        self.n1_support=dict(self.n1_support);self.n1_support.pop('lower_to_buffer_order')
        self.n1_support.update(lower_to_O3_order='rminus<R110<Rm<Rh<Rref<Rbuffer_end<=R<=Rc; reserved repair starts at Rc',
            O3_order='Rw is transition endpoint; Rc=Rw*exp(2), power x<=2',
            exact_all17_preRc_homogeneous_n1_source_used=True,continuation_beyond_Rc_admitted=False)

    def evaluate(self,Z,chart,x,*,log_tau='-1',scaled_time=None,azimuth='0'):
        result=OriginalRmPatchVelocity._assemble(self,self.numeric.evaluate(Z,chart,x),Z,log_tau,scaled_time,azimuth)
        physical=result['physical'];offset=11 if chart=='O3_transition' else 12
        definition='log(R/110)=10logCstar+10logPstar+exp40+'+str(offset)+'+x'
        physical['native'].update(chart='actual_original_'+chart,O3_x=x,R_definition='110*exp('+definition.split('=')[1]+')',
            x_definition='R/Rm=exp('+str(6+offset)+'+exp40+x)')
        def rewrite(v):
            if isinstance(v,dict):
                if v.get('kind')=='original_Rm_geometry_y_log':v.update(kind='original_O3_geometry_y_log',y_definition=definition)
                for q in v.values():rewrite(q)
            elif isinstance(v,list):
                for q in v:rewrite(q)
        rewrite(physical);return result

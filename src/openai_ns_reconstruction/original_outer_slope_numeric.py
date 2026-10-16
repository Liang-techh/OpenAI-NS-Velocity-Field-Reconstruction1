"""Original O2 slope values and physical sums with directed profile integration."""
from fractions import Fraction
import hashlib
import json
import math

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_finite_N_numeric import RATES,transport_step
from .original_Rm_patch_numeric import encode
from .original_Rm_patch_velocity import OriginalRmPatchVelocity
from .original_outer_reference_numeric import OriginalOuterReferenceNumeric
from .original_outer_reference_velocity import OriginalOuterReferenceVelocity
from .original_slope_scalar import OriginalSlopeScalar


class OriginalOuterSlopeNumeric:
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,reference=None,scalar_cells=512):
        self.reference=reference or OriginalOuterReferenceNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        r=self.reference;self.patch=r.patch;self.ctx,self.iv,self.J=r.ctx,r.iv,r.J
        self.family=r.family;self.owner=r.owner;self.directory=r.directory
        self.hashes=dict(r.hashes);self.code_hashes=dict(r.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_outer_reference_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:raise ValueError('Accepted same-original Rref inlet required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed Rref implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._inlets={p['Z']:p for p in accepted['samples'] if p['native_offset']=='0'};self._packets={}
        self.scalar=OriginalSlopeScalar(self.ctx,self.iv,cells=scalar_cells)
        for name in ('src/openai_ns_reconstruction/original_slope_scalar.py',
                     'src/openai_ns_reconstruction/original_outer_slope_numeric.py'):
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        for suffix in ('compliant_flat_pulse_derivatives.py','compliant_flat_pulse_derivatives_check.json',
                       'compliant_current_original_true_outer_C6_C5_functions.py'):
            path=self.directory/(PREFIX+suffix);digest=hashlib.sha256(path.read_bytes()).hexdigest()
            if self.hashes.get(path.name,digest)!=digest:raise ValueError('Changed original slope/sigma source')
            self.hashes[path.name]=digest

    def _rows(self,p,t):return self.reference._rows(p,t)
    def _jet(self,rows,directed):return self.reference._jet(rows,directed)

    def _cover(self,Z,q):
        c,iv=self.ctx,self.iv;zl,zh=_bounds(c,iv.mpf(Z));frames=[];groups={};local=[]
        for frame in self.owner.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            if frame['source_identity']!=self.family or not frame['uniform_for_all_integer_N_ge_N0'] or frame['epsilon_domain']!='[0,1/N0]':
                raise ValueError('Same-original all-N slope source ranges required')
            frames.append(frame['exact_Z_cell'])
            chart=next(p for p in frame['actual_17_chart_C3_transports'] if p['chart']=='O2_slope')
            for cell in chart['cells']:
                left,right=Fraction(*cell['exact_left']),Fraction(*cell['exact_right'])
                groups.setdefault((left,right),[]).append(cell)
                if left<=q<=right:local.append(cell)
        ordered=sorted(groups)
        if (not local or not ordered or ordered[0][0]!=0 or ordered[-1][1]!=1
            or any(a[1]!=b[0] for a,b in zip(ordered,ordered[1:]))):raise ValueError('Complete contiguous O2 slope source cover required')
        return frames,[(a,b,groups[a,b]) for a,b in ordered],local

    def _leading(self,Z,state,ctx):
        J=self.J;z=J.variable(ctx,ctx.mpf(Z),2);A=(1+z*z).reciprocal()
        P=ctx.exp(self.patch.leading.defect.logPp if ctx is self.ctx else self.patch.leading.defect.logP)
        B=z*(4/P);F,H,G,C=(state[n] for n in ('F','H','G','C'))
        E=A*F;h=A*H
        return dict(E=E,V=B,m=B,h=h,k=h*B,e=B*B-A*A*(G/2),p=A*A*(C/2))

    def evaluate(self,Z,x):
        c,iv,J=self.ctx,self.iv,self.J
        if type(x) is not str or not c.isfinite(c.mpf(x)):raise ValueError('Finite decimal O2 slope coordinate required')
        q=Fraction(x)
        if not 0<=q<=1:raise ValueError('Original O2 slope0<=x<=1 required')
        key=(Z,x)
        if key in self._packets:return self._packets[key]
        if Z not in self._inlets:self._inlets[Z]=self.reference.evaluate(Z,'0')
        inlet=self._inlets[Z];state=self.scalar.evaluate(x);frames,groups,local=self._cover(Z,q);proof=[];steps=[]
        history={n:self._jet(v,True) for n,v in inlet['finite_N_history_errors'].items()}
        for left,right,cells in groups:
            if q<=left:continue
            end=min(q,right);width=iv.mpf((end-left).numerator)/(end-left).denominator;weights={}
            for name,rate in RATES.items():
                coefficients=[]
                for k in range(3):
                    log=self.patch._log_cover([v['actual_density_ranges']['N_scaled'][name][k] for v in cells])
                    cap=self.patch._small(None if log is None else log-self.owner.logN,proof,'slope_density_'+name+'_Z'+str(k))/math.factorial(k)
                    coefficients.append(self.patch.leading.reference.leading._symmetric(cap))
                history[name],decay,mass=transport_step(iv,history[name],J(iv,coefficients),width,rate)
                weights[name]=dict(rate=rate,decay=decay,mass=mass)
            steps.append(dict(exact_left=[left.numerator,left.denominator],exact_end=[end.numerator,end.denominator],
                physical_log_width=width,own_weights=weights,overlapping_Z_covers=len(cells),
                physical_measure='d(O2_slope_x) once',incoming_already_contains_one_over_N=True,new_driver_divided_by_N_once=True))
        A=[];B=[]
        for k in range(3):
            for label,target in (('A',A),('B',B)):
                log=self.patch._log_cover([v['actual_source_and_primitive_range'][label][k] for v in local])
                target.append(self.patch._small(None if log is None else log-self.owner.logN,proof,'slope_local_'+label+'_Z'+str(k))/math.factorial(k))
        if _bounds(c,A[0])[1]>=c.mpf('.25'):raise ArithmeticError('Small original slope exponent required')
        symmetric=self.patch.leading.reference.leading._symmetric
        point=self._leading(Z,state['point'],c);truth=self._leading(Z,state['directed'],iv)
        p0p=self._jet(inlet['independent_P0_over_Pstar_squared'],False).truncate(2)
        p0t=self._jet(inlet['independent_P0_over_Pstar_squared'],True).truncate(2)
        point['pressure']=point['p']+p0p;truth['pressure']=truth['p']+p0t
        dE=truth['E']*J(iv,[symmetric(v) for v in (2*A[0],2*A[1],2*A[2]+A[1]*A[1])])
        dV=J(iv,[symmetric(v) for v in B])
        if q==0:dE,dV=(self._jet(inlet['finite_N_local_errors'][n],True) for n in ('E','V'))
        modified={n:truth[n]+history[n] for n in RATES}
        modified.update(E=truth['E']+dE,V=truth['V']+dV,pressure=truth['pressure']+history['p'])
        if q==0:
            for n in ('E','V'):modified[n]=self._jet(inlet['generic_modified_fields'][n],True)
        def similarity(fields,ctx):
            z=J.variable(ctx,ctx.mpf(Z),2);m=fields['m'];mz=J(ctx,[m[1],2*m[2]])
            P=ctx.exp(self.patch.leading.defect.logPp if ctx is c else self.patch.leading.defect.logP)
            de=ctx.mpf(self.patch.leading.reference.leading.core.delta if ctx is c else self.patch.leading.reference.leading.core.delta_iv)
            Q=((z*fields['V']*2-z*m*(1-de)-(1-z*z)*mz)/(1-z*z*de)*P).truncate(1)
            return dict(Uz=fields['V']*P,Utheta_over_Pstar=fields['E'],Mz_over_R=m*P,Ur_over_sqrt_R_over2=Q)
        sp,st=similarity(point,c),similarity(modified,iv);zero=J.constant(c,0,2);ref=self.patch.leading.reference
        packet=dict(Z=Z,chart='actual_original_O2_slope',native_slope=x,native_x='exp(6+slope_x)',source_family=self.family,
            selected_integer=inlet['selected_integer'],generic_leading_fields={n:self._rows(point[n],truth[n]) for n in point},
            generic_modified_fields={n:self._rows(point[n],modified[n]) for n in modified},
            physical_similarity={n:self._rows(sp[n],st[n]) for n in sp},
            finite_N_history_errors={n:self._rows(zero,v) for n,v in history.items()},finite_N_local_errors=dict(E=self._rows(zero,dE),V=self._rows(zero,dV)),
            finite_N_nonzero_Rref_incoming_errors=inlet['finite_N_history_errors'],physical_transport_steps=encode(steps),
            independent_P0_over_Pstar_squared=inlet['independent_P0_over_Pstar_squared'],positive_order_n1_own_history=inlet['positive_order_n1_own_history'],
            scalar_profiles={n:ref._summary(state['point'][n],state['directed'][n]) for n in state['point']},
            scalar_integration={k:v for k,v in state.items() if k not in ('point','directed')},
            coordinate=dict(slope_x=x,global_y_log_terms=[ref._summary(10*ref.logCp,10*ref.logC),
                ref._summary(10*self.patch.leading.defect.logPp,10*self.patch.leading.defect.logP),ref._summary(c.mpf(x),iv.mpf(x))]),
            source_cover=dict(exact_Z_cells=frames,unique_native_cells=len(groups),local_source_cells=len(local),original_all_N_ge_N0=True,phase_not_selected=True),
            small_error_proofs=proof,contract=dict(original_N_retained=True,pressure_memory_exactly_one=True,P0_added_once=True,
                leading_history_initials_from_exact_Rref_identity=True,actual_finite_N_history_inlets_not_reset=True,
                local_E_V_are_source_values_not_additively_transported=True,exact_Rref_local_E_V_reused=q==0),
            limits=dict(physical_Cartesian_assembly=False,resolved_signed_finite_N_centers=False,total_radial_or_fast_phase_derivatives=False,
                beyond_O2_slope=False,reduced_m13_m15_integrals_continued=False,global_n1_high_order_joins=False,
                global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=packet;return packet


class OriginalOuterSlopeVelocity(OriginalOuterReferenceVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,scalar_cells=512):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.numeric=OriginalOuterSlopeNumeric(reference=self.numeric,scalar_cells=scalar_cells)
        self.hashes.update(self.numeric.hashes);self.code_hashes.update(self.numeric.code_hashes)
        if _bounds(self.ctx,self.amplitude.lp_iv)[0]<=0:raise ValueError('Positive logPstar required before reserved repair')
        self.n1_support=dict(self.n1_support)
        self.n1_support.pop('lower_to_reference_order')
        self.n1_support.update(lower_to_slope_order='rminus<R110<Rm<Rh<Rref<=R<=e*Rref<Rw<Rc',
            slope_outer_order='log(Rw/(e*Rref))=logPstar>0',
            positive_log_Rw_over_slope_end=self._summary((self.amplitude.lp,self.amplitude.lp_iv)))

    def evaluate(self,Z,x,*,log_tau='-1',scaled_time=None,azimuth='0'):
        result=OriginalRmPatchVelocity._assemble(self,self.numeric.evaluate(Z,x),Z,log_tau,scaled_time,azimuth)
        physical=result['physical'];physical['native'].update(chart='actual_original_O2_slope',slope_x=x,
            R_definition='110*exp(10logCstar+10logPstar+slope_x)',x_definition='R/Rm=exp(6+slope_x)')
        def rewrite(v):
            if isinstance(v,dict):
                if v.get('kind')=='original_Rm_geometry_y_log':v.update(kind='original_O2_slope_geometry_y_log',y_definition='log(R/110)=10logCstar+10logPstar+slope_x')
                for q in v.values():rewrite(q)
            elif isinstance(v,list):
                for q in v:rewrite(q)
        rewrite(physical);return result

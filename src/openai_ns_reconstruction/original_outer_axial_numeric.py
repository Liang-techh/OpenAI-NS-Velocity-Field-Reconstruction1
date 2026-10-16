"""Original O2 axial turn-off with positive own-memory and finite-N errors."""
from fractions import Fraction
import hashlib
import json
import math

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_axial_memory import OriginalAxialMemory
from .original_Rm_finite_N_numeric import RATES,transport_step
from .original_Rm_patch_numeric import encode
from .original_Rm_patch_velocity import OriginalRmPatchVelocity
from .original_outer_slope_numeric import OriginalOuterSlopeNumeric,OriginalOuterSlopeVelocity


class OriginalOuterAxialNumeric:
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,slope=None):
        self.slope=slope or OriginalOuterSlopeNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        p=self.slope;self.patch=p.patch;self.ctx,self.iv,self.J=p.ctx,p.iv,p.J
        self.family=p.family;self.owner=p.owner;self.directory=p.directory
        self.hashes=dict(p.hashes);self.code_hashes=dict(p.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_outer_slope_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:raise ValueError('Accepted same-original O2 slope endpoint required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed O2 slope implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._inlets={p['Z']:p for p in accepted['samples'] if p['native_slope']=='1'};self._packets={}
        self.memory=OriginalAxialMemory(self.ctx,self.iv,self.J)
        for name in ('src/openai_ns_reconstruction/original_axial_memory.py','src/openai_ns_reconstruction/original_outer_axial_numeric.py'):
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _rows(self,p,t):return self.slope._rows(p,t)
    def _jet(self,rows,directed):return self.slope._jet(rows,directed)

    def _cover(self,Z,q):
        c,iv=self.ctx,self.iv;zl,zh=_bounds(c,iv.mpf(Z));frames=[];groups={};local=[]
        for frame in self.owner.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            if frame['source_identity']!=self.family or not frame['uniform_for_all_integer_N_ge_N0'] or frame['epsilon_domain']!='[0,1/N0]':
                raise ValueError('Same-original all-N axial source ranges required')
            frames.append(frame['exact_Z_cell'])
            chart=next(p for p in frame['actual_17_chart_C3_transports'] if p['chart']=='O2_axial')
            for cell in chart['cells']:
                left,right=Fraction(*cell['exact_left']),Fraction(*cell['exact_right'])
                groups.setdefault((left,right),[]).append(cell)
                if left<=q<=right:local.append(cell)
        ordered=sorted(groups)
        if (not local or not ordered or ordered[0][0]!=0 or ordered[-1][1]!=1
            or any(a[1]!=b[0] for a,b in zip(ordered,ordered[1:]))):raise ValueError('Complete contiguous O2 axial cover required')
        return frames,[(a,b,groups[a,b]) for a,b in ordered],local

    def _leading(self,Z,s,inlet,memories,ctx):
        c,iv,J=self.ctx,self.iv,self.J;directed=ctx is iv
        restore=lambda r:self.patch.leading.reference.previous._restore(r) if directed else c.mpf(r['approximation'])
        F1,H1,G1,C1=(restore(inlet['scalar_profiles'][n]) for n in ('F','H','G','C'))
        z=J.variable(ctx,ctx.mpf(Z),2);A=(1+z*z).reciprocal()
        P=ctx.exp(self.patch.leading.defect.logP if directed else self.patch.leading.defect.logPp);B=z*(4/P)
        d=ctx.expm1(40*ctx.mpf(s));key='directed' if directed else 'point'
        M,W=memories[1][key],memories[2][key]
        if ctx.mpf(s)==0:beta=ctx.mpf(1)
        elif ctx.mpf(s)==1:beta=ctx.mpf(0)
        else:beta=ctx.exp(self.memory._log_beta(ctx,J.constant(ctx,1-ctx.mpf(s),4))[0])
        F=F1*ctx.exp(-d/2)
        H=F+(H1-F1)*ctx.exp(-ctx.mpf('1.5')*d)
        K=F*(M+(H1/F1-1)*ctx.exp(-d))
        G=(G1+F1*F1*d)*ctx.exp(-d)
        C=C1+F1*F1*(-ctx.expm1(-d))
        fields=dict(E=A*F,V=B*beta,m=B*M,h=A*H,k=A*B*K,e=B*B*W-A*A*(G/2),p=A*A*(C/2))
        return fields,dict(beta=beta,M=M,W=W,F=F,H=H,K=K,G=G,C=C,physical_log_distance=d)

    def evaluate(self,Z,s):
        c,iv,J=self.ctx,self.iv,self.J
        if type(s) is not str or not c.isfinite(c.mpf(s)):raise ValueError('Finite decimal O2 axial coordinate required')
        q=Fraction(s)
        if not 0<=q<=1:raise ValueError('Original O2 axial0<=s<=1 required')
        key=(Z,s)
        if key in self._packets:return self._packets[key]
        if Z not in self._inlets:self._inlets[Z]=self.slope.evaluate(Z,'1')
        inlet=self._inlets[Z];memories={k:self.memory.evaluate(s,k) for k in (1,2)}
        frames,groups,local=self._cover(Z,q);proof=[];steps=[]
        history={n:self._jet(v,True) for n,v in inlet['finite_N_history_errors'].items()}
        rational=lambda v:iv.mpf(v.numerator)/v.denominator
        for left,right,cells in groups:
            if q<=left:continue
            end=min(q,right)
            # dd=40exp(40s)ds. Integrate the actual physical logarithmic
            # width exactly once; a native width is not a physical width.
            width=iv.exp(40*rational(left))*iv.expm1(40*rational(end-left));weights={}
            for name,rate in RATES.items():
                coefficients=[]
                for k in range(3):
                    log=self.patch._log_cover([v['actual_density_ranges']['N_scaled'][name][k] for v in cells])
                    cap=self.patch._small(None if log is None else log-self.owner.logN,proof,'axial_density_'+name+'_Z'+str(k))/math.factorial(k)
                    coefficients.append(self.patch.leading.reference.leading._symmetric(cap))
                history[name],decay,mass=transport_step(iv,history[name],J(iv,coefficients),width,rate)
                weights[name]=dict(rate=rate,decay=decay,mass=mass)
            steps.append(dict(exact_native_left=[left.numerator,left.denominator],exact_native_end=[end.numerator,end.denominator],
                physical_log_width=width,own_weights=weights,overlapping_Z_covers=len(cells),
                physical_measure='dd=40exp(40s)ds once',incoming_already_contains_one_over_N=True,new_driver_divided_by_N_once=True))
        capsA=[];capsB=[]
        for k in range(3):
            for label,target in (('A',capsA),('B',capsB)):
                log=self.patch._log_cover([v['actual_source_and_primitive_range'][label][k] for v in local])
                target.append(self.patch._small(None if log is None else log-self.owner.logN,proof,'axial_local_'+label+'_Z'+str(k))/math.factorial(k))
        if _bounds(c,capsA[0])[1]>=c.mpf('.25'):raise ArithmeticError('Small original axial exponent required')
        point,profilep=self._leading(Z,s,inlet,memories,c);truth,profilet=self._leading(Z,s,inlet,memories,iv)
        p0p=self._jet(inlet['independent_P0_over_Pstar_squared'],False).truncate(2)
        p0t=self._jet(inlet['independent_P0_over_Pstar_squared'],True).truncate(2)
        point['pressure']=point['p']+p0p;truth['pressure']=truth['p']+p0t
        symmetric=self.patch.leading.reference.leading._symmetric
        dE=truth['E']*J(iv,[symmetric(v) for v in (2*capsA[0],2*capsA[1],2*capsA[2]+capsA[1]*capsA[1])])
        dV=J(iv,[symmetric(v) for v in capsB])
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
        packet=dict(Z=Z,chart='actual_original_O2_axial',native_axial=s,native_x='exp(6+exp(40s))',source_family=self.family,
            selected_integer=inlet['selected_integer'],generic_leading_fields={n:self._rows(point[n],truth[n]) for n in point},
            generic_modified_fields={n:self._rows(point[n],modified[n]) for n in modified},physical_similarity={n:self._rows(sp[n],st[n]) for n in sp},
            finite_N_history_errors={n:self._rows(zero,v) for n,v in history.items()},finite_N_local_errors=dict(E=self._rows(zero,dE),V=self._rows(zero,dV)),
            finite_N_nonzero_slope_incoming_errors=inlet['finite_N_history_errors'],physical_transport_steps=encode(steps),
            independent_P0_over_Pstar_squared=inlet['independent_P0_over_Pstar_squared'],positive_order_n1_own_history=inlet['positive_order_n1_own_history'],
            axial_profiles={n:ref._summary(profilep[n],profilet[n]) for n in profilep},
            memory_proofs={str(k):encode(memories[k]['proof']) for k in memories},
            original_slope_scalar_inlet=inlet['scalar_profiles'],coordinate=dict(axial_s=s,
                global_y_log_terms=[ref._summary(10*ref.logCp,10*ref.logC),ref._summary(10*self.patch.leading.defect.logPp,10*self.patch.leading.defect.logP),
                    ref._summary(c.exp(40*c.mpf(s)),iv.exp(40*iv.mpf(s)))]),
            source_cover=dict(exact_Z_cells=frames,unique_native_cells=len(groups),local_source_cells=len(local),original_all_N_ge_N0=True,phase_not_selected=True),
            small_error_proofs=proof,contract=dict(original_N_retained=True,pressure_memory_exactly_one=True,P0_added_once=True,
                actual_leading_and_finite_N_inlet_memories_not_reset=True,resonant_swirl_square_history_retained=True,
                local_E_V_are_source_values_not_additively_transported=True,exact_slope_local_E_V_reused=q==0,
                local_terminal_axial_zero_does_not_zero_mass_or_square_histories=True),
            limits=dict(physical_Cartesian_assembly=False,resolved_signed_finite_N_centers=False,total_radial_or_fast_phase_derivatives=False,
                beyond_O2_axial=False,reduced_m13_m15_integrals_continued=False,global_n1_high_order_joins=False,
                global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=packet;return packet


class OriginalOuterAxialVelocity(OriginalOuterSlopeVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.numeric=OriginalOuterAxialNumeric(slope=self.numeric)
        self.hashes.update(self.numeric.hashes);self.code_hashes.update(self.numeric.code_hashes)
        gap=self.amplitude.lp_iv+1-self.iv.exp(40)
        if _bounds(self.ctx,gap)[0]<=0:raise ValueError('Axial endpoint must precede Rw and Rc')
        self.n1_support=dict(self.n1_support);self.n1_support.pop('lower_to_slope_order')
        self.n1_support.update(lower_to_axial_order='rminus<R110<Rm<Rh<Rref<e*Rref<=R<=Rref*exp(exp40)<Rw<Rc',
            axial_outer_order='log(Rw/Raxial_end)=logPstar+1-exp40=12 from original logPstar=exp40+11',
            positive_log_Rw_over_axial_end=self._summary((self.amplitude.lp+1-self.ctx.exp(40),gap)))

    def evaluate(self,Z,s,*,log_tau='-1',scaled_time=None,azimuth='0'):
        result=OriginalRmPatchVelocity._assemble(self,self.numeric.evaluate(Z,s),Z,log_tau,scaled_time,azimuth)
        physical=result['physical'];physical['native'].update(chart='actual_original_O2_axial',axial_s=s,
            R_definition='110*exp(10logCstar+10logPstar+exp(40s))',x_definition='R/Rm=exp(6+exp(40s))')
        def rewrite(v):
            if isinstance(v,dict):
                if v.get('kind')=='original_Rm_geometry_y_log':v.update(kind='original_O2_axial_geometry_y_log',y_definition='log(R/110)=10logCstar+10logPstar+exp(40s)')
                for q in v.values():rewrite(q)
            elif isinstance(v,list):
                for q in v:rewrite(q)
        rewrite(physical);return result

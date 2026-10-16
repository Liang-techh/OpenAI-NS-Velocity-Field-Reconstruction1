"""Original O2 buffer: exact memories, retained small sectors, selected N."""
from fractions import Fraction
import hashlib
import json
import math

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_finite_N_numeric import RATES,transport_step
from .original_Rm_patch_numeric import encode
from .original_Rm_patch_velocity import OriginalRmPatchVelocity
from .original_outer_axial_numeric import OriginalOuterAxialNumeric,OriginalOuterAxialVelocity


class OriginalOuterBufferNumeric:
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,axial=None):
        self.axial=axial or OriginalOuterAxialNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        a=self.axial;self.patch=a.patch;self.ctx,self.iv,self.J=a.ctx,a.iv,a.J
        self.family=a.family;self.owner=a.owner;self.directory=a.directory
        self.hashes=dict(a.hashes);self.code_hashes=dict(a.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_outer_axial_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:raise ValueError('Accepted same-original axial endpoint required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed axial implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._inlets={p['Z']:p for p in accepted['samples'] if p['native_axial']=='1'};self._packets={}
        name='src/openai_ns_reconstruction/original_outer_buffer_numeric.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _rows(self,p,t):return self.axial._rows(p,t)
    def _jet(self,rows,directed):return self.axial._jet(rows,directed)

    def _cover(self,Z,q):
        c,iv=self.ctx,self.iv;zl,zh=_bounds(c,iv.mpf(Z));frames=[];groups={};local=[]
        for frame in self.owner.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            if frame['source_identity']!=self.family or not frame['uniform_for_all_integer_N_ge_N0'] or frame['epsilon_domain']!='[0,1/N0]':
                raise ValueError('Same-original all-N buffer ranges required')
            frames.append(frame['exact_Z_cell'])
            chart=next(p for p in frame['actual_17_chart_C3_transports'] if p['chart']=='O2_buffer')
            for cell in chart['cells']:
                left,right=Fraction(*cell['exact_left']),Fraction(*cell['exact_right'])
                groups.setdefault((left,right),[]).append(cell)
                if left<=q<=right:local.append(cell)
        ordered=sorted(groups)
        if (not local or not ordered or ordered[0][0]!=0 or ordered[-1][1]!=11
            or any(a[1]!=b[0] for a,b in zip(ordered,ordered[1:]))):raise ValueError('Complete original buffer source cover required')
        return frames,[(a,b,groups[a,b]) for a,b in ordered],local

    def _leading(self,Z,b,inlet,ctx):
        c,iv,J=self.ctx,self.iv,self.J;directed=ctx is iv
        restore=lambda r:self.patch.leading.reference.previous._restore(r) if directed else c.mpf(r['approximation'])
        F1,H1,G1,C1=(restore(inlet['original_slope_scalar_inlet'][n]) for n in ('F','H','G','C'))
        M,W,D=(restore(inlet['axial_profiles'][n]) for n in ('M','W','physical_log_distance'))
        z=J.variable(ctx,ctx.mpf(Z),2);A=(1+z*z).reciprocal()
        P=ctx.exp(self.patch.leading.defect.logP if directed else self.patch.leading.defect.logPp);B=z*(4/P)
        x=ctx.mpf(b);decay=ctx.exp(-x);half=ctx.exp(-x/2);threehalf=ctx.exp(-ctx.mpf('1.5')*x)
        F0=F1*ctx.exp(-D/2);offset=(H1-F1)*ctx.exp(-ctx.mpf('1.5')*D)
        E=A*(F0*half);zero=J.constant(ctx,0,2)
        sectors=dict(h_minus_E=A*(offset*threehalf),k_mass_history=A*B*(F0*M*threehalf),
            k_initial_offset=A*B*(offset*threehalf),e_axial_square_history=B*B*(W*decay),
            e_incoming_swirl_square=-(A*A)*((G1+F1*F1*D)*ctx.exp(-D)*decay/2),
            e_buffer_swirl_square_forcing=-(A*A)*(F1*F1*ctx.exp(-D)*x*decay/2),
            p_constant=A*A*((C1+F1*F1)/2),p_tail=-(A*A)*(F1*F1*ctx.exp(-D)*decay/2),
            p_buffer_increment=A*A*(F1*F1*ctx.exp(-D)*(-ctx.expm1(-x))/2))
        fields=dict(E=E,V=zero,m=B*(M*decay),h=E+sectors['h_minus_E'],
            k=sectors['k_mass_history']+sectors['k_initial_offset'],
            e=sectors['e_axial_square_history']+sectors['e_incoming_swirl_square']+sectors['e_buffer_swirl_square_forcing'],
            p=sectors['p_constant']+sectors['p_tail'])
        # V=0 and m=4Z*M*exp(-b)/P give Q=-4M*exp(-b) identically.
        Q=J.constant(ctx,-4*M*decay,1)
        return fields,sectors,Q

    def evaluate(self,Z,b):
        c,iv,J=self.ctx,self.iv,self.J
        if type(b) is not str or not c.isfinite(c.mpf(b)):raise ValueError('Finite decimal buffer coordinate required')
        q=Fraction(b)
        if not 0<=q<=11:raise ValueError('Original O2 buffer0<=b<=11 required')
        key=(Z,b)
        if key in self._packets:return self._packets[key]
        if Z not in self._inlets:self._inlets[Z]=self.axial.evaluate(Z,'1')
        inlet=self._inlets[Z];frames,groups,local=self._cover(Z,q);proof=[];steps=[]
        history={n:self._jet(v,True) for n,v in inlet['finite_N_history_errors'].items()}
        for left,right,cells in groups:
            if q<=left:continue
            end=min(q,right);width=iv.mpf((end-left).numerator)/(end-left).denominator;weights={}
            for name,rate in RATES.items():
                coefficients=[]
                for k in range(3):
                    log=self.patch._log_cover([v['actual_density_ranges']['N_scaled'][name][k] for v in cells])
                    cap=self.patch._small(None if log is None else log-self.owner.logN,proof,'buffer_density_'+name+'_Z'+str(k))/math.factorial(k)
                    coefficients.append(self.patch.leading.reference.leading._symmetric(cap))
                history[name],decay,mass=transport_step(iv,history[name],J(iv,coefficients),width,rate)
                weights[name]=dict(rate=rate,decay=decay,mass=mass)
            steps.append(dict(exact_left=[left.numerator,left.denominator],exact_end=[end.numerator,end.denominator],
                physical_log_width=width,own_weights=weights,overlapping_Z_covers=len(cells),physical_measure='db once',
                incoming_already_contains_one_over_N=True,new_driver_divided_by_N_once=True))
        capsA=[];capsB=[]
        for k in range(3):
            for label,target in (('A',capsA),('B',capsB)):
                log=self.patch._log_cover([v['actual_source_and_primitive_range'][label][k] for v in local])
                target.append(self.patch._small(None if log is None else log-self.owner.logN,proof,'buffer_local_'+label+'_Z'+str(k))/math.factorial(k))
        if _bounds(c,capsA[0])[1]>=c.mpf('.25'):raise ArithmeticError('Small original buffer exponent required')
        point,sectorsp,Qp=self._leading(Z,b,inlet,c);truth,sectorst,Qt=self._leading(Z,b,inlet,iv)
        p0p=self._jet(inlet['independent_P0_over_Pstar_squared'],False).truncate(2)
        p0t=self._jet(inlet['independent_P0_over_Pstar_squared'],True).truncate(2)
        point['pressure']=point['p']+p0p;truth['pressure']=truth['p']+p0t
        sectorsp['P0_plus_p_constant']=p0p+sectorsp['p_constant'];sectorst['P0_plus_p_constant']=p0t+sectorst['p_constant']
        symmetric=self.patch.leading.reference.leading._symmetric
        dE=truth['E']*J(iv,[symmetric(v) for v in (2*capsA[0],2*capsA[1],2*capsA[2]+capsA[1]*capsA[1])])
        dV=J(iv,[symmetric(v) for v in capsB])
        if q==0:dE,dV=(self._jet(inlet['finite_N_local_errors'][n],True) for n in ('E','V'))
        modified={n:truth[n]+history[n] for n in RATES}
        modified.update(E=truth['E']+dE,V=truth['V']+dV,pressure=truth['pressure']+history['p'])
        if q==0:
            for n in ('E','V'):modified[n]=self._jet(inlet['generic_modified_fields'][n],True)
        Pp,Pt=c.exp(self.patch.leading.defect.logPp),iv.exp(self.patch.leading.defect.logP)
        z=J.variable(iv,iv.mpf(Z),2);dm=history['m'];dmz=J(iv,[dm[1],2*dm[2]])
        de=iv.mpf(self.patch.leading.reference.leading.core.delta_iv)
        dQ=((z*dV*2-z*dm*(1-de)-(1-z*z)*dmz)/(1-z*z*de)*Pt).truncate(1)
        sp=dict(Uz=point['V']*Pp,Utheta_over_Pstar=point['E'],Mz_over_R=point['m']*Pp,Ur_over_sqrt_R_over2=Qp)
        st=dict(Uz=modified['V']*Pt,Utheta_over_Pstar=modified['E'],Mz_over_R=modified['m']*Pt,Ur_over_sqrt_R_over2=Qt+dQ)
        zero=J.constant(c,0,2);ref=self.patch.leading.reference
        packet=dict(Z=Z,chart='actual_original_O2_buffer',native_buffer=b,native_x='exp(6+exp40+b)',source_family=self.family,
            selected_integer=inlet['selected_integer'],generic_leading_fields={n:self._rows(point[n],truth[n]) for n in point},
            generic_modified_fields={n:self._rows(point[n],modified[n]) for n in modified},physical_similarity={n:self._rows(sp[n],st[n]) for n in sp},
            leading_small_sectors={n:self._rows(sectorsp[n],sectorst[n]) for n in sectorsp},
            leading_Q_exact_Z_independent=self._rows(Qp,Qt),finite_N_history_errors={n:self._rows(zero,v) for n,v in history.items()},
            finite_N_local_errors=dict(E=self._rows(zero,dE),V=self._rows(zero,dV)),finite_N_nonzero_axial_incoming_errors=inlet['finite_N_history_errors'],
            physical_transport_steps=encode(steps),independent_P0_over_Pstar_squared=inlet['independent_P0_over_Pstar_squared'],
            positive_order_n1_own_history=inlet['positive_order_n1_own_history'],original_slope_scalar_inlet=inlet['original_slope_scalar_inlet'],
            original_axial_profile_inlet=inlet['axial_profiles'],coordinate=dict(buffer_b=b,
                global_y_log_terms=[ref._summary(10*ref.logCp,10*ref.logC),ref._summary(10*self.patch.leading.defect.logPp,10*self.patch.leading.defect.logP),
                    ref._summary(c.exp(40),iv.exp(40)),ref._summary(c.mpf(b),iv.mpf(b))]),
            source_cover=dict(exact_Z_cells=frames,unique_native_cells=len(groups),local_source_cells=len(local),original_all_N_ge_N0=True,phase_not_selected=True),
            small_error_proofs=proof,contract=dict(original_N_retained=True,pressure_memory_exactly_one=True,P0_added_once=True,
                leading_homogeneous_forced_sectors_and_resonance_retained=True,tiny_pressure_increment_exported_before_large_sum=True,
                P0_plus_pressure_constant_not_asserted_zero=True,leading_Q_minus_four_mass_scalar_from_streamfunction_identity=True,
                actual_finite_N_history_inlets_not_reset=True,local_E_V_are_source_values_not_additively_transported=True),
            limits=dict(physical_Cartesian_assembly=False,resolved_signed_finite_N_centers=False,total_radial_or_fast_phase_derivatives=False,
                beyond_O2_buffer=False,reduced_m13_m15_integrals_continued=False,relative_pressure_cancellation_accuracy=False,
                global_n1_high_order_joins=False,global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=packet;return packet


class OriginalOuterBufferVelocity(OriginalOuterAxialVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.numeric=OriginalOuterBufferNumeric(axial=self.numeric)
        self.hashes.update(self.numeric.hashes);self.code_hashes.update(self.numeric.code_hashes)
        gap=self.amplitude.lp_iv-self.iv.exp(40)-10
        if _bounds(self.ctx,gap)[0]<=0:raise ValueError('Buffer endpoint must precede Rw and Rc')
        self.n1_support=dict(self.n1_support);self.n1_support.pop('lower_to_axial_order')
        self.n1_support.update(lower_to_buffer_order='rminus<R110<Rm<Rh<Rref<Raxial_end<=R<=Rbuffer_end<Rw<Rc',
            buffer_outer_order='log(Rw/Rbuffer_end)=logPstar-exp40-10=1 from original logPstar=exp40+11',
            positive_log_Rw_over_buffer_end=self._summary((self.amplitude.lp-self.ctx.exp(40)-10,gap)))

    def evaluate(self,Z,b,*,log_tau='-1',scaled_time=None,azimuth='0'):
        result=OriginalRmPatchVelocity._assemble(self,self.numeric.evaluate(Z,b),Z,log_tau,scaled_time,azimuth)
        physical=result['physical'];physical['native'].update(chart='actual_original_O2_buffer',buffer_b=b,
            R_definition='110*exp(10logCstar+10logPstar+exp40+b)',x_definition='R/Rm=exp(6+exp40+b)')
        def rewrite(v):
            if isinstance(v,dict):
                if v.get('kind')=='original_Rm_geometry_y_log':v.update(kind='original_O2_buffer_geometry_y_log',y_definition='log(R/110)=10logCstar+10logPstar+exp40+b')
                for q in v.values():rewrite(q)
            elif isinstance(v,list):
                for q in v:rewrite(q)
        rewrite(physical);return result

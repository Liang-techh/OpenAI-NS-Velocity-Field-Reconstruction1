"""Same original selected-N Rh-to-Rref values, with own history transport.

Leading fields are the source's smooth reference identity. Signed finite-N
centers remain unresolved; full original all-N caps enclose their deviations.
"""
import hashlib
import json
import math
from fractions import Fraction

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_finite_N_numeric import OriginalRmFiniteNNumeric,RATES,transport_step
from .original_Rm_patch_numeric import encode


class OriginalOuterReferenceNumeric:
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,patch=None):
        self.patch=patch or OriginalRmFiniteNNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        p=self.patch;self.ctx,self.iv,self.J=p.ctx,p.iv,p.J
        self.family=p.family;self.owner=p.owner;self.directory=p.leading.directory
        self.hashes=dict(p.hashes);self.code_hashes=dict(p.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_Rm_finite_N_numeric.json')
        saved=json.loads(path.read_text(encoding='utf-8'))
        if not saved['checks']['passed'] or saved['source_family']!=self.family:raise ValueError('Accepted actual selected-N Rh inlet required')
        for name,digest in saved['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed Rh implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._inlets={v['Z']:v for v in saved['samples'] if v['native_x']=='Rh'};self._packets={}
        for suffix,gate in (
            ('compliant_current_original_true_outer_C6_C5_functions_check.json','current_original_true_six_outer_charts_finite_N_C5_functions_installed'),
            ('compliant_current_original_whole_Z_Rh_functional_join_check.json','current_whole_Z_repaired_leading_Rh_reference_function_identity_checked')):
            path=self.directory/(PREFIX+suffix);receipt=json.loads(path.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or not receipt[gate] or receipt['source_family']!=self.family:
                raise ValueError('Same-source outer reference identity required')
            for name,digest in receipt['input_hashes'].items():
                if hashlib.sha256((self.directory/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed outer reference source: '+name)
                if self.hashes.get(name,digest)!=digest:raise ValueError('Inconsistent outer source pin')
                self.hashes[name]=digest
            self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        name='src/openai_ns_reconstruction/original_outer_reference_numeric.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _rows(self,p,t):return self.patch._rows(p,t)
    def _jet(self,rows,directed):return self.patch._jet(rows,directed)

    def _cover(self,Z,offset):
        c,iv=self.ctx,self.iv;lo,hi=_bounds(c,iv.mpf(Z));groups={};frames=[];local=[]
        for frame in self.owner.ranges:
            left,right=map(c.mpf,frame['exact_Z_cell'])
            if right<lo or left>hi:continue
            if frame['source_identity']!=self.family or not frame['uniform_for_all_integer_N_ge_N0'] or frame['epsilon_domain']!='[0,1/N0]':
                raise ValueError('Original uniform all-N reference range required')
            frames.append(frame['exact_Z_cell'])
            chart=next(q for q in frame['actual_17_chart_C3_transports'] if q['chart']=='Rh_reference')
            for cell in chart['cells']:
                a,b=Fraction(*cell['exact_left']),Fraction(*cell['exact_right'])
                groups.setdefault((a,b),[]).append(cell)
                if a<=offset<=b:local.append(cell)
        ordered=sorted(groups)
        if (not frames or not local or ordered[0][0]!=-5 or ordered[-1][1]!=0
            or any(a[1]!=b[0] for a,b in zip(ordered,ordered[1:]))):
            raise ValueError('Complete contiguous original reference range required')
        return frames,[(a,b,groups[a,b]) for a,b in ordered],local

    def _leading(self,Z,offset,ctx):
        J=self.J;z=J.variable(ctx,ctx.mpf(Z),2)
        P=ctx.exp(self.patch.leading.defect.logPp if ctx is self.ctx else self.patch.leading.defect.logP)
        E=(1+z*z).reciprocal()*ctx.exp(ctx.mpf(offset)/10);V=z*(4/P)
        h=E*(ctx.mpf(5)/8);p=E*E*(ctx.mpf(5)/2)
        return dict(E=E,V=V,m=V,h=h,k=h*V,e=V*V-E*E*(ctx.mpf(5)/12),p=p)

    def evaluate(self,Z,offset):
        if type(offset) is not str:raise ValueError('Exact decimal reference offset required')
        c,iv,J=self.ctx,self.iv,self.J
        if not c.isfinite(c.mpf(offset)):raise ValueError('Finite reference offset required')
        query=Fraction(offset)
        if not -5<=query<=0:raise ValueError('Original Rh/Rref offsets[-5,0] required')
        key=(Z,offset)
        if key in self._packets:return self._packets[key]
        if Z not in self._inlets:self._inlets[Z]=self.patch.evaluate(Z,'Rh')
        inlet=self._inlets[Z];frames,groups,local=self._cover(Z,query);proof=[];steps=[]
        incoming={n:self._jet(rows,True) for n,rows in inlet['finite_N_history_errors'].items()}
        history=dict(incoming)
        for left,right,cells in groups:
            if query<=left:continue
            end=min(query,right);width=iv.mpf((end-left).numerator)/(end-left).denominator;weights={}
            for name,rate in RATES.items():
                bounds=[]
                for k in range(3):
                    log=self.patch._log_cover([cell['actual_density_ranges']['N_scaled'][name][k] for cell in cells])
                    bound=self.patch._small(None if log is None else log-self.owner.logN,proof,'reference_density_'+name+'_Z'+str(k))/math.factorial(k)
                    bounds.append(self.patch.leading.reference.leading._symmetric(bound))
                history[name],decay,mass=transport_step(iv,history[name],J(iv,bounds),width,rate)
                weights[name]=dict(rate=rate,decay=decay,mass=mass)
            steps.append(dict(exact_left=[left.numerator,left.denominator],exact_end=[end.numerator,end.denominator],
                physical_log_width=width,own_weights=weights,overlapping_Z_covers=len(cells),
                physical_measure='d(reference_offset) once',incoming_already_contains_one_over_N=True,new_driver_divided_by_N_once=True))
        A=[];B=[]
        for k in range(3):
            for label,target in (('A',A),('B',B)):
                log=self.patch._log_cover([cell['actual_source_and_primitive_range'][label][k] for cell in local])
                target.append(self.patch._small(None if log is None else log-self.owner.logN,proof,'reference_local_'+label+'_Z'+str(k))/math.factorial(k))
        if _bounds(c,A[0])[1]>=c.mpf('.25'):raise ArithmeticError('Small original reference exponent required')
        symmetric=self.patch.leading.reference.leading._symmetric
        exp_error=J(iv,[symmetric(q) for q in (2*A[0],2*A[1],2*A[2]+A[1]*A[1])])
        dV=J(iv,[symmetric(q) for q in B])
        point=self._leading(Z,offset,c);truth=self._leading(Z,offset,iv)
        p0p=self._jet(inlet['independent_P0_over_Pstar_squared'],False).truncate(2)
        p0t=self._jet(inlet['independent_P0_over_Pstar_squared'],True).truncate(2)
        point['pressure']=point['p']+p0p;truth['pressure']=truth['p']+p0t
        dE=truth['E']*exp_error
        if query==-5:
            # E_N/V_N are local algebraic source functions, not histories.
            # At the shared endpoint consume their actual inlet enclosures;
            # away from it evaluate the same source with local range caps.
            dE=self._jet(inlet['finite_N_local_errors']['E'],True)
            dV=self._jet(inlet['finite_N_local_errors']['V'],True)
        modified={n:truth[n]+history[n] for n in RATES}
        modified.update(E=truth['E']+dE,V=truth['V']+dV,pressure=truth['pressure']+history['p'])
        if query==-5:
            for n in ('E','V'):modified[n]=self._jet(inlet['generic_modified_fields'][n],True)
        def similarity(fields,ctx):
            z=J.variable(ctx,ctx.mpf(Z),2)
            P=ctx.exp(self.patch.leading.defect.logPp if ctx is c else self.patch.leading.defect.logP)
            delta=ctx.mpf(self.patch.leading.reference.leading.core.delta if ctx is c else self.patch.leading.reference.leading.core.delta_iv)
            m=fields['m'];mz=J(ctx,[m[1],2*m[2]])
            Q=((z*fields['V']*2-z*m*(1-delta)-(1-z*z)*mz)/(1-z*z*delta)*P).truncate(1)
            return dict(Uz=fields['V']*P,Utheta_over_Pstar=fields['E'],Ur_over_sqrt_R_over2=Q,Mz_over_R=m*P)
        sp,st=similarity(point,c),similarity(modified,iv)
        zero=J.constant(c,0,2);ref=self.patch.leading.reference
        packet=dict(Z=Z,chart='actual_original_Rh_reference',native_offset=offset,native_x='exp(6+offset)',source_family=self.family,
            selected_integer=inlet['selected_integer'],generic_leading_fields={n:self._rows(point[n],truth[n]) for n in point},
            generic_modified_fields={n:self._rows(point[n],modified[n]) for n in modified},
            physical_similarity={n:self._rows(sp[n],st[n]) for n in sp},
            finite_N_history_errors={n:self._rows(zero,v) for n,v in history.items()},
            finite_N_local_errors=dict(E=self._rows(zero,dE),V=self._rows(zero,dV)),
            finite_N_nonzero_Rh_incoming_errors=inlet['finite_N_history_errors'],physical_transport_steps=encode(steps),
            independent_P0_over_Pstar_squared=inlet['independent_P0_over_Pstar_squared'],positive_order_n1_own_history=inlet['positive_order_n1_own_history'],
            coordinate=dict(reference_offset=offset,global_y_log_terms=[ref._summary(10*ref.logCp,10*ref.logC),
                ref._summary(10*self.patch.leading.defect.logPp,10*self.patch.leading.defect.logP),ref._summary(c.mpf(offset),iv.mpf(offset))]),
            source_cover=dict(exact_Z_cells=frames,unique_native_cells=len(groups),local_source_cells=len(local),
                original_all_N_ge_N0=True,original_phase_preserved=True,phase_not_selected=True),small_error_proofs=proof,
            contract=dict(leading_smooth_function_identity_from_same_original_Rh_repair=True,own_histories_not_reset=True,
                local_E_V_are_source_values_not_additively_transported=True,exact_Rh_local_E_V_enclosures_reused=query==-5,
                original_N_retained=True,pressure_memory_exactly_one=True,P0_added_once=True,zero_point_error_centers_not_source_zeros=True),
            limits=dict(physical_Cartesian_assembly=False,resolved_signed_finite_N_centers=False,total_radial_or_fast_phase_derivatives=False,
                beyond_Rref=False,reduced_m13_m15_integrals_continued=False,
                global_n1_high_order_joins=False,global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=packet;return packet

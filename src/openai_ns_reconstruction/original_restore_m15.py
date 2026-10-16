"""Changing-source restoration and postrestore cumulative m15 to Rm."""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_macro_m15 import half_S1_coefficients,densities,combine
from .original_reference_m15 import OriginalReferenceM15
from .original_restore_kernels import OriginalRestoreKernels,powers


def pressure_polynomials(J,ctx,Z,delta,U,M):
    """Exact half-S1 F/G coefficients; key=(epsilon degree, scalar code)."""
    z=J.variable(ctx,ctx.mpf(Z),3);zero=J.constant(ctx,0,3);A,D=1-z*z*delta,1-z*z
    derivative=lambda row:J(ctx,[(k+1)*row[k+1] for k in range(row.order)])
    Q={}
    for key in set(U)|set(M):
        u,m=U.get(key,zero),M.get(key,zero)
        Q[key]=((u*(2*z)-m*(z*(1-delta))-derivative(m)*D)/A).truncate(2)
    U={key:row.truncate(2) for key,row in U.items()}
    bc,vc=half_S1_coefficients(J,ctx,Z,delta);rho=densities(U,Q,J)
    return dict(F=combine(bc,rho),G=combine(vc,rho))


class OriginalRestoreM15(OriginalReferenceM15):
    def __init__(self, *, reference_engine=None, chart='restoration', previous=None, cells=128, **kwargs):
        if reference_engine is not None and kwargs:
            raise ValueError('Supply reference moment source or constructor options')
        if chart not in ('restoration','postrestore'):
            raise ValueError('Original restoration or postrestore chart required')
        self.reference_engine=reference_engine if reference_engine is not None else OriginalReferenceM15(**kwargs)
        base=self.reference_engine
        for name in ('t','o','ctx','iv','J','audit','directory','family'):
            setattr(self,name,getattr(base,name))
        self.previous=base if chart=='restoration' else previous
        if self.previous is None or (chart=='postrestore' and self.previous.chart!='restoration'):
            raise ValueError('Full original restoration moment inlet required')
        self.rows=self.previous.rows.next_field();self.chart=chart
        self.hashes,self.code_hashes=dict(base.hashes),dict(base.code_hashes)
        self.hashes.update(self.rows.hashes);self.code_hashes.update(self.rows.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_reference_m15_numeric.json')
        raw=path.read_bytes();data=json.loads(raw)
        if (not data['checks']['passed'] or data['source_family']!=self.family
                or not data['limits']['reference_cumulative_signed_m15_with_complete_factored_errors']):
            raise ValueError('Accepted complete124-term reference moment required')
        for root,group,target in ((self.directory,data['input_hashes'],self.hashes),(ROOT,data['implementation_hashes'],self.code_hashes)):
            for name,digest in group.items():
                if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed original restoration prerequisite: '+name)
                target[name]=digest
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        self.inlets={Fraction(p['Z']):p for p in data['samples']
            if p['coordinate']==dict(kind='exact_decimal_reference',value='1')} if chart=='restoration' else {}
        self.kernels=OriginalRestoreKernels(self.rows.reference,cells=cells)
        self._packets={};self._polynomials={};self._next=None
        for name in ('original_restore_kernels.py','original_restore_m15.py'):
            name='src/openai_ns_reconstruction/'+name
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def next_field(self):
        if self.chart!='restoration':raise ValueError('Already at the original Rm endpoint')
        if self._next is None:
            self._next=OriginalRestoreM15(reference_engine=self.reference_engine,chart='postrestore',previous=self,cells=self.kernels.cells)
        return self._next

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            self.inlets[key]=self.previous.evaluate(Z,'1')
            self.hashes.update(self.previous.hashes);self.code_hashes.update(self.previous.code_hashes)
        return self.inlets[key]

    def polynomials(self,Z):
        if Z in self._polynomials:return self._polynomials[Z]
        incoming=self.rows.inputs(Z,self.iv.mpf(0))['incoming']
        proof=[];memory=None
        if self.chart=='restoration':
            ref=self.rows.reference;long_incoming=self.rows.previous_rows.inputs(Z,self.iv.mpf(1))['incoming']
            memory=ref._decay(long_incoming['point']['Mpsi']-long_incoming['point']['Psi'],
                long_incoming['truth']['Mpsi']-long_incoming['truth']['Psi'],1,ref.lengthp,ref.length,
                proof,'restoration_correlated_reference_mean_difference')
        result={}
        for index,(ctx,values,de,j) in enumerate(((self.ctx,incoming['point'],self.o.core.delta,self.o.core.j),
                                                 (self.iv,incoming['truth'],self.o.core.delta_iv,self.o.core.j_iv))):
            J=self.J;z=J.variable(ctx,ctx.mpf(Z),3)
            if self.chart=='restoration':
                psi=values['Psi'].truncate(3)
                U={(0,0):4*z,(0,1):J.constant(ctx,j,3),(1,1):psi}
                M={(0,0):4*z,(0,3):J.constant(ctx,j,3),(1,3):psi,(1,9):memory[index].truncate(3)}
            else:
                U={(0,0):4*z}
                M={(0,0):4*z,(0,9):values['M0'].truncate(3)-4*z,(1,9):values['Mpsi'].truncate(3)}
            result['point' if index==0 else 'truth']=pressure_polynomials(J,ctx,Z,de,U,M)
        result['reference_memory_proofs']=proof;self._polynomials[Z]=result;return result

    def integral(self,Z,q):
        c,iv,J=self.ctx,self.iv,self.J;poly=self.polynomials(Z);geo=self.rows.geometry(q)
        t=geo['local_point'];source=geo['local']
        if self.chart=='restoration':
            kernels=self.kernels.evaluate(t,source);kp,kt=kernels['point'],kernels['truth']
        else:
            kp,kt={},{}
            for key in poly['point']['G']:
                _,code=key;i,j,k=powers(code)
                if i or j:raise ValueError('Postrestore has only actual decaying mean modes')
                rate=2-k
                kp[code]=t if rate==0 else c.expm1(rate*t)/rate
                kt[code]=source if rate==0 else iv.expm1(rate*source)/rate
            kernels=dict(point=kp,truth=kt,cell_count=0,uncertain_endpoint_segment=False,
                contract=dict(exact_postrestore_exponential_resonant_kernels=True))
        result={}
        for name,ctx,ks in (('point',c,kp),('truth',iv,kt)):
            out={d:J.constant(ctx,0,1) for d in range(3)}
            for (d,code),row in poly[name]['G'].items():out[d]+=row*ks[code]
            result[name]=out
        result['kernels']=kernels;return result

    def finite_errors(self,Z,q,packet):
        iv=self.iv;_,base,geo=self.rows.background(Z,q);owner=self.previous.rows
        lower_packet=owner.evaluate(Z,'1');_,lower_base,lower_geo=owner.background(Z,iv.mpf(1))
        upper=self.correction_densities(Z,packet,base);lower=self.correction_densities(Z,lower_packet,lower_base)
        volume=self.correction_densities(Z,packet,self.rows.uniform_background(Z,q))['volume']
        result={}
        for name,power in (('endpoint',2),('viscous',1)):
            result[name]=[self.rows.sum_logs((None if upper[name][k] is None else upper[name][k]+power*geo['logR'],
                None if lower[name][k] is None else lower[name][k]+power*lower_geo['logR'])) for k in (0,1)]
        mass=iv.mpf(min(_bounds(self.audit,q)[1],self.audit.mpf('.5')))
        result['volume']=[None if v is None else v+2*geo['logR']+iv.log(mass) for v in volume]
        return result,dict(endpoint=self.rows.log_rows(lower['endpoint']),viscous=self.rows.log_rows(lower['viscous']),
            log_radius=self.o._summary(lower_geo['logR_point'],lower_geo['logR']))

    def evaluate(self,Z,coordinate):
        packet=self.rows.evaluate(Z,coordinate);q,_,recipe=self.rows._coordinate(coordinate,False)
        key=Z,coordinate
        if key in self._packets:return self._packets[key]
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;incoming=self.inlet(Z);count=124 if self.chart=='restoration' else 140
        if len(incoming['m15_C0_C1_terms'])!=count or incoming['actual_selected_integer']!=packet['selected_integer']:
            raise ValueError('Complete same-N incoming moment state required')
        empty=q._mpi_==iv.mpf(0)._mpi_;terms=copy.deepcopy(incoming['m15_C0_C1_terms'])
        zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1);geo=self.rows.geometry(q)
        upper=self.rows.background(Z,q);lower=self.rows.background(Z,iv.mpf(0))
        startoffset=-8 if self.chart=='restoration' else -7
        def atoms(power,degree,local=0):
            ref=self.rows.reference
            return [o._summary(power*c.log(110)-degree*o.core.logLambda,power*iv.log(110)-degree*o.core.ll_iv),
                    o._summary(power*10*ref.logCp,power*10*ref.logC),
                    o._summary(power*10*ref.logPp,power*10*ref.logP),
                    o._summary(c.mpf(power*startoffset),iv.mpf(power*startoffset)),
                    o._summary(power*(geo['local_point'] if local else 0),power*(geo['local'] if local else iv.mpf(0)))]
        def append(name,p,t,scale,**binding):
            terms.append(dict(name=self.chart+'_'+name,coefficient_C0_C1=o._jet_summary(p.truncate(1),t.truncate(1)),
                log_scale_terms=scale,actual_function_binding=dict(chart=self.chart,incoming_terms_retained_once=count,
                    epsilon_degree_kept_before_scaling=True,original_native_unit_interval=True,**binding)))
        def endpoint_F(ctx,values):
            bc,_=half_S1_coefficients(J,ctx,Z,o.core.delta if ctx is c else o.core.delta_iv)
            U={(0,0):values['U0'].truncate(2),(1,0):values['Psi'].truncate(2)}
            Q={(0,0):values['Q0'].truncate(2),(1,0):values['Qpsi'].truncate(2)}
            return combine(bc,densities(U,Q,J))
        if not empty:
            integral=self.integral(Z,q);finite,lower_budget=self.finite_errors(Z,q,packet)
        else:
            integral=None;finite={name:[None,None] for name in ('endpoint','volume','viscous')};lower_budget=None
        for label,values,sign,local in (('upper',upper,-1,1),('lower',lower,1,0)):
            # F endpoints have the opposite orientation to -[R Q_y].
            fp,ft=endpoint_F(c,values[0]),endpoint_F(iv,values[1])
            for d in range(3):
                append('R2F_'+label+'_epsilon_'+str(d),zero if empty else fp.get((d,0),zero)*(-sign),
                    zeroi if empty else ft.get((d,0),zeroi)*(-sign),atoms(2,d,local),
                    exact_endpoint_identity='integral R2*(2F+Fy)=[R2F]',signed_endpoint_orientation=-sign)
        for d in range(3):
            append('G_integral_epsilon_'+str(d),zero if empty else integral['point'][d],
                zeroi if empty else integral['truth'][d],atoms(2,d),
                positive_scalar_kernels=True,quadrature_error_in_coefficient_enclosure=True,
                no_extra_half_on_existing_half_S1_coefficients=True)
        for label,values,sign,local in (('upper',upper,-1,1),('lower',lower,1,0)):
            for d,name in ((0,'Q0_y'),(1,'Qpsi_y')):
                append('RQy_'+label+'_epsilon_'+str(d),zero if empty else values[0][name]*sign,
                    zeroi if empty else values[1][name]*sign,atoms(1,d,local),
                    signed_endpoint_orientation=sign,actual_changing_source_endpoint=True)
        for name in ('endpoint','volume','viscous'):
            bound=self.rows.sum_logs(finite[name]);scale=[];truth=zeroi
            if bound is not None:
                top=iv.mpf(_bounds(self.audit,bound)[1]);atom=o._summary(_bounds(c,top)[1],top)
                atom.update(positive_error_bound_atom_only=True,source_value_selected=False);scale.append(atom)
                truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
            append('finite_N_'+name+'_error',zero,truth,scale,error_enclosure_only=True,both_nonflat_endpoints=True)
        kreceipt=None
        if integral is not None:
            ks=integral['kernels'];kreceipt=dict(cell_count=ks['cell_count'],uncertain_endpoint_segment=ks['uncertain_endpoint_segment'],
                contract=ks['contract'],monomials={str(code):o._summary(ks['point'][code],value) for code,value in ks['truth'].items()})
        out=dict(Z=Z,chart=self.chart,coordinate=recipe,source_family=self.family,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']),m15_C0_C1_terms=terms,
            incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']),full_prefix_cover=copy.deepcopy(packet['source_cover']),
            scalar_kernel_receipt=kreceipt,finite_N_error_logs_in_physical_units={name:self.rows.log_rows(v) for name,v in finite.items()},
            actual_lower_finite_N_endpoint_budget=lower_budget,
            reference_mean_memory_proofs=self.rows.long._serializable_proof(self.polynomials(Z)['reference_memory_proofs']),
            contract=dict(complete_incoming_moment_once=True,changing_j_baseline_retained=True,epsilon_degrees_separate=True,
                genuine_meanZ3_QZZ=True,exact_F_endpoint_identity=True,actual_G_integral=True,
                both_signed_RQy_endpoints=True,both_finite_N_endpoints=True,positive_native_measure_once=True,
                full_source_and_quadrature_enclosures=True,unit_offsets_separate_from_huge_radius=True,
                no_epsilon_on_finite_N_pressure_error=True),
            limits=dict(restoration_cumulative_signed_m15_with_complete_absolute_errors=True,
                postrestore_cumulative_signed_m15_to_Rm=self.chart=='postrestore',
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                independent_preRp_signed_seam=False,global_five_moment_numeric_join=False,common_global_N=False,
                exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

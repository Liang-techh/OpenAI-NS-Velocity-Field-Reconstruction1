"""Actual Rm..Rh signed pressure moment, preserving the 156-term inlet."""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_pressure_rows import OriginalRmPressureRows
from .original_Rm_moment_kernels import OriginalRmMomentKernels,powers
from .original_macro_m15 import OriginalMacroM15,half_S1_coefficients,densities,combine
from .original_restore_m15 import pressure_polynomials
from .original_first_micro_pressure_rows import pressure_coefficients


class OriginalRmM15(OriginalMacroM15):
    def __init__(self,*,rows=None,cells=32,partial_cells=128,**kwargs):
        if rows is not None and kwargs:raise ValueError('Supply actual patch rows or constructor options')
        self.rows=rows if rows is not None else OriginalRmPressureRows(**kwargs)
        for name in ('t','o','ctx','iv','J','audit','directory','family'):
            setattr(self,name,getattr(self.rows,name))
        self.hashes,self.code_hashes=dict(self.rows.hashes),dict(self.rows.code_hashes)
        self.inlets={}
        for suffix in ('restore_m15_numeric','Rm_pressure_rows_numeric'):
            path=self.directory/(PREFIX+'compliant_current_original_'+suffix+'.json');raw=path.read_bytes();data=json.loads(raw)
            if not data['checks']['passed'] or data['source_family']!=self.family:raise ValueError('Accepted same-family pressure inlet/rows required')
            for base,group,target in ((self.directory,data['input_hashes'],self.hashes),(ROOT,data['implementation_hashes'],self.code_hashes)):
                for name,digest in group.items():
                    if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest or target.get(name,digest)!=digest:
                        raise ValueError('Changed active pressure prerequisite: '+name)
                    target[name]=digest
            self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
            if suffix=='restore_m15_numeric':
                self.inlets={Fraction(p['Z']):p for p in data['samples'] if p['chart']=='postrestore'
                    and p['coordinate']==dict(kind='exact_decimal_postrestore',value='1')}
        self.kernels=OriginalRmMomentKernels(self.rows,cells=cells,partial_cells=partial_cells)
        self.code_hashes.update(self.kernels.code_hashes)
        self._inlet_engine=None;self._packets={};self._polynomials={};self._prefix_errors={}
        name='src/openai_ns_reconstruction/original_Rm_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            from .original_restore_m15 import OriginalRestoreM15
            from .original_reference_m15 import OriginalReferenceM15
            from .original_reference_pressure_rows import OriginalReferencePressureRows
            if self._inlet_engine is None:
                reference=OriginalReferenceM15(rows=OriginalReferencePressureRows(reference=self.rows.reference))
                self._inlet_engine=OriginalRestoreM15(reference_engine=reference).next_field()
            self.inlets[key]=self._inlet_engine.evaluate(Z,'1')
            self.hashes.update(self._inlet_engine.hashes);self.code_hashes.update(self._inlet_engine.code_hashes)
        return self.inlets[key]

    def polynomials(self,Z):
        if Z in self._polynomials:return self._polynomials[Z]
        data=self.rows.controls.coefficients(Z);result={}
        for label,ctx,de in (('point',self.ctx,self.o.core.delta),('truth',self.iv,self.o.core.delta_iv)):
            z=self.J.variable(ctx,ctx.mpf(Z),3);h=data['normalized_'+label]
            U={(0,0):4*z,(1,1):h[0],(1,3):h[1]}
            # Exact mass-row closure: d0+h0+h1=0 on the same unique branch.
            M={(0,0):4*z,(1,9):-h[0],(1,27):-h[1]}
            result[label]=pressure_polynomials(self.J,ctx,Z,de,U,M)
        result['mass_row_diagnostic']=data['defects_truth'][0]+data['normalized_truth'][0]+data['normalized_truth'][1]
        self._polynomials[Z]=result;return result

    def normalized_fields(self,Z,x,directed):
        label='truth' if directed else 'point';ctx=self.iv if directed else self.ctx
        geo=self.rows.geometry(x);xx=geo['x' if directed else 'x_point'];q=None if x=='Rh' else Fraction(x)
        h=self.rows.controls.coefficients(Z)['normalized_'+label];J=self.J;z=J.variable(ctx,ctx.mpf(Z),3)
        b=[];tails=[]
        for i in range(2):
            right=self.kernels.centers[i]+Fraction(1,40)
            left=self.kernels.centers[i]-Fraction(1,40)
            if q is None or q>=right:beta,bx,tail=ctx.mpf(0),ctx.mpf(0),ctx.mpf(0)
            elif q<=left:beta,bx,tail=ctx.mpf(0),ctx.mpf(0),ctx.mpf(1)
            else:
                beta,bx=self.kernels.beta(ctx,i,xx,directed);tail=self.kernels.tail(i,q)[int(directed)]
            b.append((beta,bx));tails.append(tail)
        U=h[0]*b[0][0]+h[1]*b[1][0];M=-(h[0]*tails[0]+h[1]*tails[1])/xx
        Uy=(h[0]*b[0][1]+h[1]*b[1][1])*xx
        de=self.o.core.delta_iv if directed else self.o.core.delta;A,D=1-z*z*de,1-z*z
        dz=lambda row:J(ctx,[(k+1)*row[k+1] for k in range(row.order)])
        Q=((U*(2*z)-M*(z*(1-de))-dz(M)*D)/A).truncate(2)
        Qy=((Uy*(2*z)-(U-M)*(z*(1-de))-dz(U-M)*D)/A).truncate(1)
        base=4*z;Q0=((base*(z*(1+de))-dz(base)*D)/A).truncate(2)
        return dict(U0=base,Q0=Q0,U=U,Q=Q,Q_y=Qy,m=M,U_y=Uy)

    def endpoint_F(self,Z,x,directed):
        ctx=self.iv if directed else self.ctx;v=self.normalized_fields(Z,x,directed)
        bc,_=half_S1_coefficients(self.J,ctx,Z,self.o.core.delta_iv if directed else self.o.core.delta)
        return combine(bc,densities({(0,0):v['U0'].truncate(2),(1,0):v['U'].truncate(2)},
            {(0,0):v['Q0'],(1,0):v['Q']},self.J)),v

    def integral(self,Z,x):
        poly=self.polynomials(Z);ks=self.kernels.evaluate(x);result={}
        for label,ctx in (('point',self.ctx),('truth',self.iv)):
            out={d:self.J.constant(ctx,0,1) for d in range(3)}
            for (d,code),row in poly[label]['G'].items():out[d]+=row*ks[label][code]
            result[label]=out
        result['kernels']=ks;return result

    def prefix_error_packet(self,Z,x):
        key=Z,x
        if key in self._prefix_errors:return self._prefix_errors[key]
        rows=self.rows;h=rows.history_logs(Z,x);iv=self.iv;hist={name:list(v) for name,v in rows.incoming_history_logs(Z).items()}
        lookup={json.dumps([g['exact_left'],g['exact_right']]):g for g in h['groups']};cells=[]
        for step in h['steps']:
            group=lookup[json.dumps([step['exact_left'],step['exact_right']])];cells.extend(group['cells'])
            for name in hist:
                mass=step['weights'][name]['positive_mass']
                for k in range(4):
                    source=rows.cover_owner._log_cover([p['actual_density_ranges']['N_scaled'][name][k] for p in group['cells']])
                    added=None if source is None else source-rows.owner.logN+iv.log(mass)
                    hist[name][k]=rows.sum_logs((hist[name][k],added))
        if not cells:cells=h['local']
        B=[rows.cover_owner._log_cover([cell['actual_source_and_primitive_range']['B'][k] for cell in cells]) for k in range(4)]
        lp=self.o.core.lp_iv
        errors=dict(U=[None if v is None else v+lp-rows.owner.logN for v in B],
                    m=[None if v is None else v+lp for v in hist['m']])
        errors['U_Z'],errors['m_Z']=errors['U'][1:],errors['m'][1:]
        slow=rows.loop['slow_phase_held_log_bounds']['B'];fast=rows.loop['phase_and_first_slow_phase_log_bounds']
        shift=lambda value,offset:None if value is None else value+offset
        errors['U_y']=[rows.sum_logs((shift(rows.read_log(slow[s]),lp-rows.owner.logN),shift(rows.read_log(fast[f]),lp)))
            for s,f in (('y1_Z0','B_phi'),('y1_Z1','B_Z_phi'))]
        coeff=pressure_coefficients(self.J,iv,Z,self.o.core.delta_iv,2)
        errors['Q']=rows.combine(coeff['Q'],errors,2);errors['Q_y']=rows.combine(coeff['Q_y'],errors,1)
        result=dict(physical_pressure_rows={name:dict(finite_N_absolute_error_log_bounds=rows.log_rows(errors[name]))
                    for name in ('U','Q','Q_y')},
            finite_N_prefix_history_logs={name:rows.log_rows(v) for name,v in hist.items()},
            contract=dict(incoming_decay_omitted_only_for_prefix_supremum=True,all_positive_driver_masses_retained=True,
                local_B_bounds_over_full_prefix=True,slow_fast_phase_bounds_over_whole_chart=True,
                rate_zero_mass_is_native_width=True,bounds_not_selected_signed_values=True))
        self._prefix_errors[key]=result;return result

    def finite_errors(self,Z,x,packet):
        rows=self.rows;_,up,geo=rows.background(Z,x);lowpacket=rows.evaluate(Z,'1');_,low,lowgeo=rows.background(Z,'1')
        upper=self.correction_densities(Z,packet,up);lower=self.correction_densities(Z,lowpacket,low)
        prefix=self.prefix_error_packet(Z,x)
        volume=self.correction_densities(Z,prefix,rows.uniform_background(Z,x))['volume'];result={}
        for name,power in (('endpoint',2),('viscous',1)):
            result[name]=[rows.sum_logs((None if upper[name][k] is None else upper[name][k]+power*geo['logR'],
                None if lower[name][k] is None else lower[name][k]+power*lowgeo['logR'])) for k in (0,1)]
        mass=(geo['x']*geo['x']-1)/2
        result['volume']=[None if v is None else v+2*lowgeo['logR']+self.iv.log(mass) for v in volume]
        return result,dict(endpoint=rows.log_rows(lower['endpoint']),viscous=rows.log_rows(lower['viscous']),
            log_radius=self.o._summary(lowgeo['logR_point'],lowgeo['logR']),prefix_error_packet=prefix)

    def evaluate(self,Z,x):
        packet=self.rows.evaluate(Z,x);key=Z,x
        if key in self._packets:return self._packets[key]
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;incoming=self.inlet(Z)
        if len(incoming['m15_C0_C1_terms'])!=156 or incoming['actual_selected_integer']!=packet['selected_integer']:
            raise ValueError('Complete same-N156-term original Rm inlet required')
        empty=x=='1';terms=copy.deepcopy(incoming['m15_C0_C1_terms']);geo=self.rows.geometry(x)
        zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1);ref=self.rows.reference;patch=self.rows.patch
        def atoms(power,degree,upper=False):
            return [o._summary(power*c.log(110),power*iv.log(110)),
                o._summary(power*10*ref.logCp,power*10*ref.logC),
                o._summary(power*10*ref.logPp,power*10*ref.logP),
                o._summary(-6*power*c.mpf(1),-6*power*iv.mpf(1)),
                o._summary(degree*c.log(patch.tp),degree*iv.log(patch.t)),
                o._summary(power*(geo['y_point'] if upper else 0),power*(geo['y'] if upper else iv.mpf(0)))]
        def append(name,p,t,scale,**binding):
            terms.append(dict(name='Rm_patch_'+name,coefficient_C0_C1=o._jet_summary(p.truncate(1),t.truncate(1)),
                log_scale_terms=scale,actual_function_binding=dict(chart='actual_patch',incoming156_terms_retained_once=True,
                    scale_degree_is_original_fixed_axial_t_not_epsilon=True,**binding)))
        upperp,up= self.endpoint_F(Z,x,False);uppert,ut=self.endpoint_F(Z,x,True)
        lowerp,lp= self.endpoint_F(Z,'1',False);lowert,lt=self.endpoint_F(Z,'1',True)
        integral=self.integral(Z,x);poly=self.polynomials(Z);ks=integral['kernels']
        if not empty:finite,lower_budget=self.finite_errors(Z,x,packet)
        else:finite={name:[None,None] for name in ('endpoint','volume','viscous')};lower_budget=None
        bp=(2*poly['point']['F'][(0,0)]+poly['point']['G'][(0,0)])*ks['point'][0]
        bt=(2*poly['truth']['F'][(0,0)]+poly['truth']['G'][(0,0)])*ks['truth'][0]
        append('constant_background_integral',zero if empty else bp,zeroi if empty else bt,atoms(2,0),
            exact_native_primitive='(2F0+G0)*(x²-1)/2',no_large_endpoint_subtraction=True)
        for degree in (1,2):
            for label,sign,p,t,upper in (('upper',1,upperp,uppert,True),('lower',-1,lowerp,lowert,False)):
                append('R2F_'+label+'_t'+str(degree),zero if empty else p[(degree,0)]*sign,
                    zeroi if empty else t[(degree,0)]*sign,atoms(2,degree,upper),signed_endpoint_orientation=sign)
            append('G_integral_t'+str(degree),zero if empty else integral['point'][degree],
                zeroi if empty else integral['truth'][degree],atoms(2,degree),
                actual_positive_scalar_kernels=True,quadrature_error_in_coefficient_enclosure=True,
                no_extra_half_on_existing_half_S1_coefficients=True)
        for label,p,t,sign,upper in (('upper',up,ut,-1,True),('lower',lp,lt,1,False)):
            append('RQy_'+label+'_t1',zero if empty else p['Q_y']*sign,zeroi if empty else t['Q_y']*sign,
                atoms(1,1,upper),signed_endpoint_orientation=sign,actual_mean_ODE_pressure_endpoint=True)
        for name in ('endpoint','volume','viscous'):
            bound=self.rows.sum_logs(finite[name]);truth=zeroi;scale=[]
            if bound is not None:
                top=iv.mpf(_bounds(self.audit,bound)[1]);atom=o._summary(_bounds(c,top)[1],top)
                atom.update(positive_error_bound_atom_only=True,source_value_selected=False);scale.append(atom)
                truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
            append('finite_N_'+name+'_error',zero,truth,scale,error_enclosure_only=True,
                both_nonflat_endpoints=True,volume_uses_complete_prefix_not_final_decayed_error=True)
        out=dict(Z=Z,chart='actual_patch',native_x=x,source_family=self.family,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']),m15_C0_C1_terms=terms,
            incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']),full_prefix_cover=copy.deepcopy(packet['source_cover']),
            scalar_kernel_receipt=dict(cells_per_support=self.kernels.cells,partial_CDF_cells=self.kernels.partial_cells,
                contract=ks['contract'],kernels={str(code):o._summary(ks['point'][code],value) for code,value in ks['truth'].items()}),
            finite_N_error_logs_in_physical_units={name:self.rows.log_rows(v) for name,v in finite.items()},
            actual_lower_finite_N_endpoint_budget=lower_budget,
            exact_mass_row_branch_diagnostic=o._jet_summary(J.constant(c,0,3),poly['mass_row_diagnostic']),
            contract=dict(complete156term_inlet_once=True,actual_partial_patch_pressure_integrated=True,
                original_mass_constraint_from_same_unique_branch=True,control_scale_degrees_separate=True,
                exact_constant_background_integral=True,exact_F_endpoint_identity=True,actual_G_integral=True,
                positive_dx_over_x_and_RdR_once=True,both_signed_RQy_endpoints=True,both_finite_N_endpoints=True,
                complete_prefix_errors_not_final_endpoint_caps=True,full_source_and_quadrature_enclosures=True,
                no_extra_epsilon_N_hb_or_half=True,huge_physical_radius_not_materialized=True),
            limits=dict(Rm_Rh_cumulative_signed_m15_with_complete_absolute_errors=True,
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                independent_preRp_signed_seam=False,global_five_moment_numeric_join=False,common_global_N=False,
                exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

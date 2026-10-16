"""Physical original selected-N Rm/Rh leading+n1 Cartesian velocity sums.

The geometry y means log(R/110), not logR. The inherited constant n1
primitive has no radius factor in its physical streamfunction.
"""
import copy
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_core_velocity import OriginalCoreVelocity
from .original_long_velocity import OriginalLongVelocity
from .original_Rm_finite_N_numeric import OriginalRmFiniteNNumeric


class OriginalRmPatchVelocity(OriginalLongVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.numeric=OriginalRmFiniteNNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.amplitude=self.numeric.owner.amplitude;self.ctx,self.iv=self.amplitude.ctx,self.amplitude.iv
        self.J=self.numeric.J;self.family=self.numeric.family
        self.hashes=dict(self.numeric.hashes);self.code_hashes=dict(self.numeric.code_hashes)
        directory=self.numeric.leading.directory
        path=directory/(PREFIX+'compliant_current_original_Rm_finite_N_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:
            raise ValueError('Accepted actual finite-N patch numerical source required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed finite-N patch implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._packets={(v['Z'],v['native_x']):v for v in accepted['samples']}
        # The actual lower cutoff is at rho=4 exp(hb*sc/2). Source hb<1,
        # 0<sc<=1, Lambda>1 imply Rminus<4exp(1/2)<110. The accepted
        # positive T and L, then restoration lengths1+1, place Rm beyond110.
        transport=self.numeric.owner.transport
        sc=self.iv.mpf(mp.make_mpf(tuple(transport.sc_tuple)))
        logh=self.numeric.owner.leading.logh;logLambda=self.amplitude.ll_iv
        if not (0<_bounds(self.ctx,sc)[0]<=_bounds(self.ctx,sc)[1]<=1
                and _bounds(self.ctx,logh)[1]<0 and _bounds(self.ctx,logLambda)[0]>0):
            raise ValueError('Original lower-cutoff domain ordering required')
        ref=self.numeric.leading.reference
        if _bounds(self.ctx,ref.previous.T)[0]<=0 or _bounds(self.ctx,ref.length)[0]<=0:
            raise ValueError('Original positive long/reference lengths required')
        # Consume checked source functions for this homogeneous local branch.
        # This does not promote their still-open global high-order field joins.
        receipts=(('compliant_current_generic_shear_loop_domain_check.json',
                   'current_original_generic_loop_two_sided_collars_and_repair_geometry_certified'),
                  ('compliant_current_original_n1_intermediate_own_C3_chart_functions_check.json',
                   'current_original_n1_lower_and_17_preRc_own_C3_QV_C2_chart_functions_installed'),
                  ('compliant_current_original_n1_reserved_linear_operator_check.json',
                   'current_original_n1_reserved_normalized_linear_repair_operator_installed'))
        for suffix,gate in receipts:
            receipt_path=directory/(PREFIX+suffix)
            receipt=json.loads(receipt_path.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or not receipt[gate] or receipt['source_family']!=self.family:
                raise ValueError('Checked original homogeneous n1 source required: '+suffix)
            for name,digest in receipt['input_hashes'].items():
                if hashlib.sha256((directory/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed homogeneous n1 source: '+name)
                if self.hashes.get(name,digest)!=digest:raise ValueError('Inconsistent n1 source pin')
                self.hashes[name]=digest
            self.hashes[receipt_path.name]=hashlib.sha256(receipt_path.read_bytes()).hexdigest()
        gap=self.amplitude.lp_iv+6
        if _bounds(self.ctx,gap)[0]<=0:raise ValueError('Rh must precede Rw and Rc')
        self.n1_support=dict(cutoff_phase_at_rminus='sc/2; twice original inner phase_in=sc/4',
            cutoff_zero_fraction='3/4',source_hb_lt1=True,source_sc_in_0_1=True,Lambda_gt1=True,
            Rminus_upper='4*exp(1/2)/Lambda <110',original_T_and_reference_L_positive=True,
            lower_to_patch_order='rminus < R110 < Rm <= R <= Rh; before later Rc reserved repair',
            outer_order='log(Rw/Rh)=logPstar+6>0; Rc=Rw*exp(2); reserved supports are beyond Rc',
            positive_log_Rw_over_Rh=self._summary((self.amplitude.lp+6,gap)),
            primitive_ODE='dJ/dlog(rho)=rho*u1; local u1=0 after the exact cutoff and before later repair',
            own_mean_identity='mean1_y=-mean1, hence (rho*mean1)_y=0',
            local_homogeneous_branch_only=True,global_high_order_n1_field_joins_admitted=False,
            constant_primitive_is_source_continuation_not_data_copy=True)
        for name in ('src/openai_ns_reconstruction/original_n1_cutoff_numeric.py',
                     'src/openai_ns_reconstruction/original_long_velocity.py',
                     'src/openai_ns_reconstruction/original_core_velocity.py'):
            digest=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            if self.code_hashes.get(name,digest)!=digest:raise ValueError('Changed physical/cutoff source')
            self.code_hashes[name]=digest
        name='src/openai_ns_reconstruction/original_Rm_patch_velocity.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _assemble(self,packet,Z,log_tau,scaled_time,azimuth):
        c,iv,J=self.ctx,self.iv,self.J
        if type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)):raise ValueError('Finite decimal azimuth required')
        if packet['source_family']!=self.family:raise ValueError('Same original physical source required')
        history=packet['positive_order_n1_own_history'];terminal=history['actual_terminal_primitive']
        if not history['f_and_u_exact_zero'] or not terminal['exact_same_terminal_history_used']:
            raise ValueError('Actual post-cutoff homogeneous n1 history required')
        ell,de,time=OriginalCoreVelocity._time(self,Z,log_tau,scaled_time)
        ll=self.amplitude.ll,self.amplitude.ll_iv;lp=self.amplitude.lp,self.amplitude.lp_iv
        yterms=[self._pair(row) for row in packet['coordinate']['global_y_log_terms']]
        primitive=terminal['actual_primitive']
        pp=J(c,[self._pair(r)[0] for r in primitive]);pt=J(iv,[self._pair(r)[1] for r in primitive])
        def recover(j,ctx,delta):
            z=J.variable(ctx,ctx.mpf(Z),j.order)
            jz=J(ctx,[(i+1)*j[i+1] for i in range(j.order)])
            return (-z*j*(1+delta)-(1-z*z)*jz)/(1-z*z*delta)
        qp,qt=recover(pp,c,self.amplitude.delta),recover(pt,iv,self.amplitude.delta_iv)
        jp,jt=(pp[0],pt[0]),(qp[0],qt[0])
        def scaled(pair,k):return pair[0]*k,pair[1]*k
        def logs(a,b,L=0,P=0,constant=None,yfactor=None):
            out=[self._summary(scaled(q,a)) for q in ell if a]
            out.extend(self._summary(scaled(q,b)) for q in de if b)
            if L:out.append(self._summary(scaled(ll,L)))
            if P:out.append(self._summary(scaled(lp,P)))
            if constant is not None:out.append(self._summary(constant))
            if yfactor is not None:
                for i,q in enumerate(yterms):
                    row=self._summary(scaled(q,c.mpf(yfactor)))
                    row.update(kind='original_Rm_geometry_y_log',component=i,multiplier=yfactor,
                        y_definition='log(R/110)=10logCstar+10logPstar-6+log(x)',unit_offsets_kept_separate=True)
                    out.append(row)
            return out
        def term(name,order,pair,weights):
            return dict(name=name,temporal_order=order,coefficient=self._summary(pair),log_scale_terms=weights)
        zero=lambda name:term(name,1,(c.mpf(0),iv.mpf(0)),[])
        coefficients={n:self._pair(packet['physical_similarity'][n][0]) for n in
            ('Uz','Utheta_over_Pstar','Ur_over_sqrt_R_over2','Mz_over_R')}
        radial=[term('leading_radial',0,coefficients['Ur_over_sqrt_R_over2'],
                    logs(-1,0,constant=(c.ln(55)/2,iv.ln(55)/2),yfactor='.5')),
            term('n1_radial',1,jt,logs(-1,2,L=-2,constant=(-c.ln(220)/2,-iv.ln(220)/2),yfactor='-.5'))]
        swirl=[term('leading_swirl',0,coefficients['Utheta_over_Pstar'],logs(-1,-1,P=1)),zero('n1_swirl')]
        axial=[term('leading_axial',0,coefficients['Uz'],logs(-1,-1)),zero('n1_axial')]
        stream=[term('leading_mass',0,coefficients['Mz_over_R'],
                    logs(1,-1,constant=(c.ln(110),iv.ln(110)),yfactor='1')),
            term('n1_mass',1,jp,logs(1,1,L=-2))]
        radius=term('physical_radius',None,(c.mpf(1),iv.mpf(1)),
            logs(1,0,constant=(c.ln(220)/2,iv.ln(220)/2),yfactor='.5'))
        z,zi=c.mpf(Z),iv.mpf(Z);angle,anglei=c.mpf(azimuth),iv.mpf(azimuth)
        cosine,sine=(c.cos(angle),iv.cos(anglei)),(c.sin(angle),iv.sin(anglei))
        negsine=-sine[0],-sine[1]
        leading_pressure=term('leading_pressure',0,self._pair(packet['generic_modified_fields']['pressure'][0]),logs(-2,-2,P=2))
        physical=dict(source_family=self.family,time=time,
            native=dict(chart='actual_finite_N_Rm_patch',Z=Z,x=packet['native_x'],azimuth=azimuth,
                y_definition='log(R/110)',R_definition='110*exp(10logCstar+10logPstar-6)*x',
                geometry_atoms=packet['coordinate']['global_y_log_terms'],physical_viscosity=1,
                original_N_and_construction_parameters_retained=True),
            coordinates=dict(r=radius,x=self._rotate(radius,cosine),y=self._rotate(radius,sine),
                z=term('physical_z',None,(z,zi),logs(1,-1)),t=dict(kind='one_minus_exp_of_log_tau',definition=time)),
            cylindrical=dict(ur=radial,utheta=swirl,uz=axial),
            Cartesian=dict(u=[self._rotate(r,cosine) for r in radial]+[self._rotate(r,negsine) for r in swirl],
                v=[self._rotate(r,sine) for r in radial]+[self._rotate(r,cosine) for r in swirl],w=copy.deepcopy(axial)),
            streamfunction=stream,leading_pressure=leading_pressure,
            n1_own_history=dict(actual_terminal_primitive=terminal,primitive_coefficient=self._summary(jp),
                Q_primitive_coefficient=self._summary(jt),constant_primitive_support_admission=self.n1_support,
                local_f_u_exact_zero=True,primitive_log_radius_derivative_exact_zero=True,
                mean_u_definition='J/(Lambda*R)',radial_identity='Q_J*lambda^(-1+2delta)*Lambda^-2/sqrt(2R)',
                streamfunction_identity='J*lambda^(1+delta)*Lambda^-2',
                geometry_cancelled_before_enclosure=True,streamfunction_has_no_radius_atom=True),
            original_selected_integer=packet['selected_integer'],coefficient_orders_included=[0,1],
            shared_streamfunction_meridional_recovery=True,actual_selected_N_errors_included=True,
            absolute_physical_exponentials_materialized=False,native_coordinates_required=True,
            direct_physical_inverse_available=False,physical_gradients_or_vorticity_available=False,
            leading_pressure_includes_same_P0_once=True,full_n1_pressure_included=False,
            whole_global_temporal_sum_admitted=False,global_field_available=False,genuine_n_ge_2_temporal_recursion=False,
            full_corrected_NS_solution=False,value_semantics='sum coefficient*exp(sum log_scale_terms), with full directed errors')
        result=copy.deepcopy(packet);result['physical']=physical
        result['limits'].update(physical_Cartesian_assembly=True,physical_pressure_leading_only=True,
            physical_gradient=False,direct_physical_inverse=False)
        return result

    def evaluate(self,Z,x,*,log_tau='-1',scaled_time=None,azimuth='0'):
        key=(Z,x)
        if key not in self._packets:self._packets[key]=self.numeric.evaluate(Z,x)
        return self._assemble(self._packets[key],Z,log_tau,scaled_time,azimuth)

    def at_log_radius(self,*args,**kwargs):raise ValueError('Use native x=R/Rm, or x=Rh')

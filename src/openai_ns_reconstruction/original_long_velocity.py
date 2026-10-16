"""Physical original finite-long leading+n1 velocity as bounded log sums.

The exact long radius never becomes an ordinary float. Retained n1 primitive
and radius factors cancel before interval arithmetic, not after separately
boxing enormous positive and negative logarithms.
"""
import copy
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_core_velocity import OriginalCoreVelocity
from .original_long_numeric import OriginalLongNumeric


class OriginalLongVelocity:
    """Native (Z,q,angle,time) -> physical Cartesian u,v,w factored sums."""
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.numeric=OriginalLongNumeric(precision=precision,
            amplitude_precision=amplitude_precision,directory=directory)
        self.amplitude=self.numeric.owner.amplitude
        self.ctx,self.iv=self.amplitude.ctx,self.amplitude.iv
        self.J=self.numeric.J;self.family=self.numeric.family
        self.hashes=dict(self.numeric.hashes);self.code_hashes=dict(self.numeric.code_hashes)
        path=self.numeric.leading.core.directory/(PREFIX+'compliant_current_original_long_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:
            raise ValueError('Accepted original finite-long numerical reader required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed accepted long numerical implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        name='src/openai_ns_reconstruction/original_long_velocity.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _restore(self,row):
        return self.iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']])

    def _pair(self,row):return self.ctx.mpf(row['approximation']),self._restore(row)

    def _summary(self,pair):return self.amplitude.summarize(*pair)

    def _rotate(self,row,pair):
        out=copy.deepcopy(row);p,t=self._pair(row['coefficient'])
        out['coefficient']=self._summary((p*pair[0],t*pair[1]))
        return out

    def _assemble(self,packet,Z,log_tau,scaled_time,azimuth):
        c,iv,J=self.ctx,self.iv,self.J
        if type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)):
            raise ValueError('Finite decimal azimuth required')
        if packet['source_family']!=self.family:raise ValueError('One physical/numerical source family required')
        ell,de,time=OriginalCoreVelocity._time(self,Z,log_tau,scaled_time)
        ll=self.amplitude.ll,self.amplitude.ll_iv
        yp,y=self._pair(packet['coordinate']['y'])
        z,zi=c.mpf(Z),iv.mpf(Z)
        u0=4*z+self.amplitude.j,4*zi+self.amplitude.j_iv
        qaxis=(((1+self.amplitude.delta)*z*u0[0]-4*(1-z*z))/(1-self.amplitude.delta*z*z),
            ((1+self.amplitude.delta_iv)*zi*u0[1]-4*(1-zi*zi))/(1-self.amplitude.delta_iv*zi*zi))
        rows=packet['normalized_modified_leading']
        psi,mean,qpsi=(self._pair(rows[n][0]) for n in ('Psi','Mpsi','Qpsi'))
        terminal=packet['positive_order_n1_own_history']['actual_terminal_primitive']
        primitive=terminal['actual_primitive']
        pp=J(c,[self._pair(row)[0] for row in primitive])
        pt=J(iv,[self._pair(row)[1] for row in primitive])
        def recover(prim,ctx,delta):
            zz=J.variable(ctx,ctx.mpf(Z),prim.order)
            derivative=J(ctx,[(k+1)*prim[k+1] for k in range(prim.order)])
            return (-zz*prim*(1+delta)-(1-zz*zz)*derivative)/(1-zz*zz*delta)
        qp,qt=recover(pp,c,self.amplitude.delta),recover(pt,iv,self.amplitude.delta_iv)
        masspair,radialpair=(pp[0],pt[0]),(qp[0],qt[0])
        def scaled(pair,k):return pair[0]*k,pair[1]*k
        def logs(a,b,L=0,constant=None,yfactor=None):
            out=[self._summary(scaled(q,a)) for q in ell if a]
            out.extend(self._summary(scaled(q,b)) for q in de if b)
            if L:out.append(self._summary(scaled(ll,L)))
            if constant is not None:out.append(self._summary(constant))
            if yfactor is not None:
                row=self._summary((yp*c.mpf(yfactor),y*iv.mpf(yfactor)))
                row.update(kind='original_finite_long_y_log',multiplier=yfactor,
                    definition='multiplier*T*q; T is original finite singleton',
                    original_T_singleton_exact_mpf_tuple=packet['original_T_singleton_exact_mpf_tuple'],
                    same_y_source_as_profile=True)
                out.append(row)
            return out
        def term(name,order,pair,exponents):
            return dict(name=name,temporal_order=order,coefficient=self._summary(pair),log_scale_terms=exponents)
        def zero(name):return term(name,1,(c.mpf(0),iv.mpf(0)),[])
        # Physical radial geometry is sqrt(R/2). Lambda cancels against
        # rho=Lambda*R analytically; the epsilon correction stays separate.
        half55=c.ln(55)/2,iv.ln(55)/2
        radial=[term('leading_radial_axis',0,qaxis,logs(-1,0,constant=half55,yfactor='.5')),
            term('leading_radial_correction',0,qpsi,logs(-1,0,-1,half55,'.5')),
            term('n1_radial',1,radialpair,logs(-1,2,-2,
                (-c.ln(220)/2,-iv.ln(220)/2),'-.5'))]
        # Full current F0*sqrt(2R)*Phi_N is already in this source log.
        swirl_log=copy.deepcopy(packet['physical_similarity_log_Utheta'][0])
        swirl_log.update(kind='original_long_similarity_swirl_log',
            F0_current_Phi_and_sqrt_2R_included_once=True)
        swirl=[term('leading_swirl',0,(c.mpf(1),iv.mpf(1)),
            [swirl_log]+logs(-1,-1)),zero('n1_swirl')]
        axial=[term('leading_axis',0,u0,logs(-1,-1)),
            term('leading_axial_correction',0,psi,logs(-1,-1,-1)),zero('n1_axial')]
        log110=c.ln(110),iv.ln(110)
        stream=[term('leading_mass_axis',0,u0,logs(1,-1,constant=log110,yfactor='1')),
            term('leading_mass_correction',0,mean,logs(1,-1,-1,log110,'1')),
            term('n1_mass',1,masspair,logs(1,1,-2))]
        # n1 radial coefficient is Q_J; its -y/2 is the combined radius and
        # memory factor. n1 mass has exactly no y dependence.
        radius=term('physical_radius',None,(c.mpf(1),iv.mpf(1)),
            logs(1,0,constant=(c.ln(220)/2,iv.ln(220)/2),yfactor='.5'))
        angle,anglei=c.mpf(azimuth),iv.mpf(azimuth)
        cosine,sine=(c.cos(angle),iv.cos(anglei)),(c.sin(angle),iv.sin(anglei))
        negsine=-sine[0],-sine[1]
        physical=dict(source_family=self.family,time=time,
            native=dict(chart='actual_modified_long_reshape',Z=Z,coordinate=packet['coordinate'],azimuth=azimuth,
                R_definition='110*exp(y)',rho_definition='110*Lambda*exp(y)',
                physical_radius_definition='lambda*sqrt(220)*exp(y/2)',
                original_finite_T_and_selected_N_retained=True),
            coordinates=dict(r=radius,x=self._rotate(radius,cosine),y=self._rotate(radius,sine),
                z=term('physical_z',None,(z,zi),logs(1,-1)),
                t=dict(kind='one_minus_exp_of_log_tau',definition=time)),
            cylindrical=dict(ur=radial,utheta=swirl,uz=axial),
            Cartesian=dict(u=[self._rotate(row,cosine) for row in radial]+[self._rotate(row,negsine) for row in swirl],
                v=[self._rotate(row,sine) for row in radial]+[self._rotate(row,cosine) for row in swirl],w=copy.deepcopy(axial)),
            streamfunction=stream,
            n1_own_history=dict(actual_terminal_primitive=terminal,
                primitive_coefficient=self._summary(masspair),Q_primitive_coefficient=self._summary(radialpair),
                mean_u_definition='J/(110*Lambda)*exp(-y)',
                radial_identity='Q_J*lambda^(-1+2delta)*Lambda^(-2)/sqrt(220)*exp(-y/2)',
                streamfunction_identity='J*lambda^(1+delta)*Lambda^(-2)',
                current_radius_and_nonzero_primitive_retained=True,
                y_cancellation_before_enclosure=True,f_and_u_exact_zero=True),
            original_selected_integer=packet['selected_integer'],
            coefficient_orders_included=[0,1],shared_streamfunction_meridional_recovery=True,
            actual_modified_leading_errors_included=True,
            absolute_F0_or_physical_velocity_exponentials_materialized=False,
            original_collar_geometry_kept_in_separate_log_terms=True,
            useful_combined_physical_relative_error_certified=False,
            native_coordinates_required=True,direct_Cartesian_inverse_available=False,
            physical_gradient_or_vorticity_available=False,
            whole_global_temporal_sum_admitted=False,global_field_available=False,
            genuine_n_ge_2_temporal_recursion=False,full_corrected_NS_solution=False,
            value_semantics='sum coefficient*exp(sum log_scale_terms); full directed coefficient/log errors, canonical phase unresolved')
        result=copy.deepcopy(packet);result['physical']=physical
        result['limits']['full_physical_vector_assembled']=True
        result['limits'].update(direct_physical_inverse=False,physical_gradient=False)
        return result

    def evaluate(self,Z,q,*,log_tau='-1',scaled_time=None,azimuth='0'):
        return self._assemble(self.numeric.evaluate(Z,q),Z,log_tau,scaled_time,azimuth)

    def at_log_radius(self,Z,y,*,log_tau='-1',scaled_time=None,azimuth='0'):
        return self._assemble(self.numeric.at_log_radius(Z,y),Z,log_tau,scaled_time,azimuth)

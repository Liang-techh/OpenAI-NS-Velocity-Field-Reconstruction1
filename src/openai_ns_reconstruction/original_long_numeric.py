"""Finite original long reshape in logarithms and current-amplitude units.

Ordinary numerical centers are proxies with directed source errors. No huge
radius, tiny Phi or finite-N fast phase is materialized. All six leading
inlet histories remain separate from finite-N and positive-order n1 data.
"""
import gzip
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_r110_numeric import OriginalSecondSwitchNumeric


def gradient_jet(J, ctx, z, j, delta, width, order):
    """The accepted anchored primitive's exact rational derivative."""
    Z=J.variable(ctx,z,order)
    H=-4*Z*Z*Z-Z*Z*j+Z*(9-delta)/2+j
    denominator=H*H+width*width
    # Preserve the positive real square before reciprocal interval algebra.
    denominator=J(ctx,[H[0]**2+width*width]+list(denominator.coefficients[1:]))
    return (1-Z*Z*delta)*H/denominator


class OriginalLongNumeric:
    """evaluate(Z,q) or at_log_radius(Z,y), with y=T*q, 0<=q<=1."""
    def __init__(self, *, precision=85, amplitude_precision=110, directory=None):
        self.owner=owner=OriginalSecondSwitchNumeric(precision=precision,
            amplitude_precision=amplitude_precision,directory=directory)
        self.inlet=owner.post_power_field()
        self.leading=owner.leading;self.transport=owner.transport
        self.ctx,self.iv=self.transport.ctx,self.transport.iv
        self.J=self.leading.Jet;self.family=owner.family
        self.hashes=dict(owner.hashes);self.code_hashes=dict(owner.code_hashes)
        self._inputs={};self._terminal={}
        path=self.leading.core.directory/(PREFIX+'compliant_current_original_r110_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:
            raise ValueError('Accepted actual original R110 numerical inlet required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed accepted R110 implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        stem=PREFIX+'compliant_current_original_prescribed_R100_Rm_C6_functions'
        receipt=self.leading.core.directory/(stem+'_check.json')
        check=json.loads(receipt.read_text(encoding='utf-8'))
        if not check['all_passed'] or check['source_family']!=self.family:
            raise ValueError('Accepted same-family finite-long source required')
        for name,digest in check['input_hashes'].items():
            if self.hashes.get(name)!=digest:
                if hashlib.sha256((self.leading.core.directory/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed finite-long prerequisite: '+name)
            self.hashes[name]=digest
        self.hashes[receipt.name]=hashlib.sha256(receipt.read_bytes()).hexdigest()
        archive=self.leading.core.directory/(stem+'.json.gz')
        raw=archive.read_bytes();digest=hashlib.sha256(raw).hexdigest()
        if self.hashes.get(archive.name,digest)!=digest:raise ValueError('Changed exact T archive')
        self.hashes[archive.name]=digest
        geometry=json.loads(gzip.decompress(raw))['prescribed_R100_Rm_C6_bindings']['geometry_certificate']
        selected=geometry['actual_T_singleton_tuple']
        if selected['lower_exact_mpf_tuple']!=selected['upper_exact_mpf_tuple']:
            raise ValueError('Original finite T must be one exact singleton')
        self.T_tuple=selected['lower_exact_mpf_tuple']
        self.Tp=self.ctx.make_mpf(tuple(self.T_tuple))
        self.T=self.iv.mpf([mp.make_mpf(tuple(self.T_tuple))]*2)
        if _bounds(self.ctx,self.T)[0]<=0:raise ValueError('Positive original finite T required')
        import lei_ren_part1_paper_compliant_long_reshape_profiles as kernels
        from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm
        self.kernels=kernels;self.logarithm=logarithm
        # Reuse certified finite kernels, never the old comparison inlet.
        for name in (PREFIX+'compliant_long_reshape_profiles.py',
                     PREFIX+'compliant_inner_bridge_profiles.py',
                     PREFIX+'compliant_anchored_axis_amplitude.py',
                     PREFIX+'compliant_anchored_axis_amplitude_check.py',
                     PREFIX+'shared_fixed_step_bound.json'):
            p=self.leading.core.directory/name;actual=hashlib.sha256(p.read_bytes()).hexdigest()
            if self.hashes.get(name,actual)!=actual:raise ValueError('Changed analytic long helper: '+name)
            self.hashes[name]=actual
        step=json.loads((self.leading.core.directory/(PREFIX+'shared_fixed_step_bound.json')).read_text(encoding='utf-8'))
        if not step['global_derivative_bound_certified'] or step['global_derivative_upper']!=8:
            raise ValueError('Original sigma derivative bound8 required')
        name='src/openai_ns_reconstruction/original_long_numeric.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _restore(self,row):
        return self.iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']])

    def _summary(self,p,t):return self.leading._summary(p,t)

    @staticmethod
    def _serializable_proof(value):
        if hasattr(value,'_mpi_'):
            return dict(lower_exact_mpf_tuple=list(value._mpi_[0]),
                        upper_exact_mpf_tuple=list(value._mpi_[1]))
        if isinstance(value,dict):return {k:OriginalLongNumeric._serializable_proof(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):return [OriginalLongNumeric._serializable_proof(v) for v in value]
        return value

    def inputs(self,Z):
        if Z in self._inputs:return self._inputs[Z]
        c,iv,J=self.ctx,self.iv,self.J
        if type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1:
            raise ValueError('Interior decimal axial Z required')
        point,truth,_=self.inlet._background(Z,iv.mpf(1))
        sixp,sixt=self.inlet._own_six(Z,iv.mpf(1))
        amplitude=self.owner.amplitude.global_amplitude(Z)
        Gp0=c.mpf(amplitude['G']['approximation']);Gt0=self._restore(amplitude['G'])
        source=self.owner.amplitude
        params=((c,c.mpf(Z),c.mpf(source.j),c.mpf(source.delta),c.mpf(source.sigma)),
                (iv,iv.mpf(Z),iv.mpf(source.j_iv),iv.mpf(source.delta_iv),iv.mpf(source.sigma_iv)))
        out=[]
        for directed,(ctx,z,j,de,width) in enumerate(params):
            grad=gradient_jet(J,ctx,z,j,de,width,2)
            G=J(ctx,[Gt0 if directed else Gp0]+[grad[k-1]/k for k in range(1,4)])
            values=truth if directed else point
            zjet=J.variable(ctx,z,3);logshape=self.logarithm(1+zjet*zjet)
            lp=self.logarithm(values['Phi'])
            Lambda=iv.mpf(source.Lambda_iv) if directed else c.mpf(source.Lambda)
            B=-G*Lambda+lp+ctx.ln(220)/2+logshape
            out.append(dict(values=values,six=sixt if directed else sixp,G=G,
                            logPhi110=lp,logshape=logshape,B=B))
        self._inputs[Z]=out
        return out

    def _cutoff_times(self,point,truth,qp,q,proof):
        """Multiply by sigma and its complement before any flat-tail cap."""
        c,iv,J=self.ctx,self.iv,self.J
        ql,qh=_bounds(c,q)
        zero_p=J.constant(c,0,point.order);zero_t=J.constant(iv,0,truth.order)
        if ql==qh==0:return zero_p,zero_t,point,truth
        if ql==qh==1:return point,truth,zero_p,zero_t
        odds=1/(1-q)**2-1/q**2
        low,high=_bounds(c,odds)
        if high<-10000 or low>10000:
            small_log=odds if high<0 else -odds
            errors=[]
            for k in range(truth.order+1):
                magnitude=self.leading._abs(truth[k])
                log=None if _bounds(c,magnitude)[1]==0 else iv.ln(magnitude)+small_log
                errors.append(self.leading._symmetric(self.owner._cap(log,proof,'flat_sigma_product_Z'+str(k))))
            tail=J(iv,errors)
            if high<0:return zero_p,tail,point,truth-tail
            return point,truth-tail,zero_p,tail
        # Stable logistic; only a bounded exponential is materialized.
        op=1/(1-qp)**2-1/qp**2
        if low>=0:
            e=iv.exp(-odds);sig=1/(1+e);comp=e/(1+e)
            ep=c.exp(-op);sp=1/(1+ep);cp=ep/(1+ep)
        else:
            e=iv.exp(odds);sig=e/(1+e);comp=1/(1+e)
            ep=c.exp(op);sp=ep/(1+ep);cp=1/(1+ep)
        return point*sp,truth*sig,point*cp,truth*comp

    def _covers(self,Z,q):
        o=self.owner;c,iv=self.ctx,self.iv
        _,prefix,frames=self.inlet._covers(Z,iv.mpf(1))
        zl,zh=_bounds(c,iv.mpf(Z));ql,qh=_bounds(c,q);local=[]
        for frame in o.ranges:
            a,b=map(c.mpf,frame['exact_Z_cell'])
            if b<zl or a>zh:continue
            current=next(w for w in frame['actual_17_chart_C3_transports'] if w['chart']=='long_reshape')
            for cell in current['cells']:
                ll,lh=_bounds(c,o._endpoint(cell['exact_left']))
                rl,rh=_bounds(c,o._endpoint(cell['exact_right']))
                if not rh<ql and not ll>qh:local.append(cell)
                if min(qh,rh)>ll:
                    width=iv.mpf(_bounds(c,iv.mpf(min(qh,rh))-iv.mpf(ll))[1])
                    prefix.append((cell,width,iv.ln(self.T),'long_reshape'))
        if not local:raise ValueError('Complete actual long source cover required')
        return local,prefix,frames

    def _six(self,inp,sigB,qp,q,yp,y,proof):
        c,iv,J=self.ctx,self.iv,self.J
        p,t=inp;pp,tt={},{};kernels={};rate_proofs={}
        initialp=dict(H_over_Phi=p['six']['H']/p['values']['Phi'],
            Kpsi_over_Phi=p['six']['Kpsi']/p['values']['Phi'],
            B_over_Phi2=p['six']['B']/(p['values']['Phi']*p['values']['Phi']),
            C_over_Phi2=p['six']['C']/(p['values']['Phi']*p['values']['Phi']))
        initialt=dict(H_over_Phi=t['six']['H']/t['values']['Phi'],
            Kpsi_over_Phi=t['six']['Kpsi']/t['values']['Phi'],
            B_over_Phi2=t['six']['B']/(t['values']['Phi']*t['values']['Phi']),
            C_over_Phi2=t['six']['C']/(t['values']['Phi']*t['values']['Phi']))
        for k,m,names in (('1.6',1,('H_over_Phi','Kpsi_over_Phi')),
                           ('1.2',2,('B_over_Phi2',)),('.2',2,('C_over_Phi2',))):
            rate=iv.mpf(k);radius=8*m*self.leading._abs(t['B'][0])/self.T
            low,high=rate-radius,rate+radius
            if _bounds(c,low)[0]<=0:raise ArithmeticError('Actual local kernel rate not positive')
            # sigma(0)=0 and |sigma'|<=8 on the entire finite window imply
            # |m*B0*sigma(q)| <= radius*y. These PERTURBED rates enclose
            # the full zeroth exponential, including the B0 term, BEFORE
            # normalized_amplitude_decay applies any exponential cap.
            rate_proofs[k]=dict(m=m,rate_min=low,rate_max=high,radius=radius,
                global_sigma_inequality='abs(sigma(q))<=8*q=8*y/T',
                full_zeroth_log_inequality='-rate_max*y<=-k*y+m*B0*sigma(q)<=-rate_min*y',
                scalar_B0_factor_included_before_cap=True)
            kernel=self.kernels.backward_kernel(iv,t['B'].truncate(2),self.T,y,m,low,high,proof)
            # This finite constant-rate center is only a numerical proxy;
            # the full variable-sigma finite kernel encloses its error.
            kp=-c.expm1(-c.mpf(k)*yp)/c.mpf(k) if yp<10000 else 1/c.mpf(k)
            kpoint=J.constant(c,kp,2)
            kernels[k]=self.leading._jet_summary(kpoint,kernel)
            for name in names:
                logs=sigB[1].truncate(2)*m-y*rate
                it=self.kernels.normalized_amplitude_decay(iv,logs,low,high,y,proof,
                    initial=initialt[name].truncate(2))
                if yp<10000:
                    ip=initialp[name].truncate(2)*(sigB[0].truncate(2)*m-c.mpf(k)*yp).exp()
                else:ip=J.constant(c,0,2)
                sourcep=(2*p['values']['Psi'].truncate(2) if name=='Kpsi_over_Phi' else
                         J.constant(c,2 if name=='H_over_Phi' else 1,2))
                sourcet=(2*t['values']['Psi'].truncate(2) if name=='Kpsi_over_Phi' else
                         J.constant(iv,2 if name=='H_over_Phi' else 1,2))
                pp[name]=ip+sourcep*kpoint;tt[name]=it+sourcet*kernel
        # Combine incoming coefficient with exp(-y) in logs before a cap.
        for name,targetp,targett in (('Mpsi',p['values']['Psi'],t['values']['Psi']),
                ('Apsi',p['values']['Psi']*p['values']['Psi'],t['values']['Psi']*t['values']['Psi'])):
            dp=(p['six'][name]-targetp).truncate(2);dt=(t['six'][name]-targett).truncate(2)
            if yp<10000:
                pp[name]=targetp.truncate(2)+dp*c.exp(-yp)
                tt[name]=targett.truncate(2)+dt*iv.exp(-y)
            else:
                errors=[]
                for k in range(3):
                    mag=self.leading._abs(dt[k]);log=None if _bounds(c,mag)[1]==0 else iv.ln(mag)-y
                    errors.append(self.leading._symmetric(self.owner._cap(log,proof,name+'_inlet_memory_Z'+str(k))))
                pp[name]=targetp.truncate(2);tt[name]=targett.truncate(2)+J(iv,errors)
        return pp,tt,kernels,rate_proofs

    def evaluate(self,Z,q):return self._evaluate(Z,q,None)

    def at_log_radius(self,Z,y):return self._evaluate(Z,None,y)

    def _evaluate(self,Z,qvalue,yvalue):
        c,iv,J=self.ctx,self.iv,self.J;o=self.owner;lead=self.leading
        raw=qvalue if yvalue is None else yvalue
        if type(raw) is not str or not c.isfinite(c.mpf(raw)):raise ValueError('Finite decimal coordinate required')
        if yvalue is None:
            qp,q=c.mpf(raw),iv.mpf(raw);yp,y=self.Tp*qp,self.T*q
        else:
            yp,y=c.mpf(raw),iv.mpf(raw);qp,q=yp/self.Tp,y/self.T
        lo,hi=_bounds(c,q)
        if lo<0 or hi>1:raise ValueError('Original finite long coordinate in[0,1] required')
        inp=self.inputs(Z);p,t=inp;proof=[]
        sbp,sbt,cbp,cbt=self._cutoff_times(p['B'],t['B'],qp,q,proof)
        lp=p['logPhi110']-sbp-c.mpf(2)*yp/5
        lt=t['logPhi110']-sbt-iv.mpf(2)*y/5
        if lo==hi==1:
            # Cancel the inlet logPhi analytically at the exact endpoint.
            lp=p['G']*c.mpf(o.amplitude.Lambda)-c.ln(220)/2-p['logshape']-c.mpf(2)*yp/5
            lt=t['G']*iv.mpf(o.amplitude.Lambda_iv)-iv.ln(220)/2-t['logshape']-iv.mpf(2)*y/5
        logup=cbp-p['logshape']+yp/10-c.mpf(o.amplitude.logC)
        logut=cbt-t['logshape']+y/10-iv.mpf(o.amplitude.logC_iv)
        local,prefix,frames=self._covers(Z,q)
        A,B,hist=o._correction_logs(local,prefix,False)
        aj=J(iv,[lead._symmetric(o._cap(A[k],proof,'long_A_over_N_Z'+str(k))/math.factorial(k)) for k in range(4)])
        bj=J(iv,[lead._symmetric(o._cap(B[k],proof,'long_delta_Psi_Z'+str(k))/math.factorial(k)) for k in range(4)])
        means=[None if v is None else v+lead.core.lp_iv+lead.core.ll_iv for v in hist['m']]
        mj=J(iv,[lead._symmetric(o._cap(means[k],proof,'long_delta_Mpsi_Z'+str(k))/math.factorial(k)) for k in range(3)])
        sixp,sixt,kernels,rate_proofs=self._six(inp,(sbp,sbt),qp,q,yp,y,proof)
        point=dict(logPhi=lp,Psi=p['values']['Psi'],Mpsi=sixp['Mpsi'])
        truth=dict(logPhi=lt+aj,Psi=t['values']['Psi']+bj,Mpsi=sixt['Mpsi']+mj)
        def recover(values,ctx,de):
            z=J.variable(ctx,ctx.mpf(Z),2);mean=values['Mpsi'];psi=values['Psi'].truncate(2)
            return ((2*z*psi-z*mean*(1-de)-(1-z*z)*lead._derivative(mean))/(1-z*z*de)).truncate(1)
        point['Qpsi']=recover(point,c,lead.core.delta);truth['Qpsi']=recover(truth,iv,lead.core.delta_iv)
        if Z not in self._terminal:self._terminal[Z]=self.transport.lower_collar(Z,'2')['normalized_own_primitive']
        return dict(Z=Z,source_family=self.family,chart='long_reshape',
            coordinate=dict(kind='exact_decimal_q' if yvalue is None else 'exact_decimal_R110_log_offset',value=raw,
                q=self._summary(qp,q),y=self._summary(yp,y),physical_measure='T*dq; y=T*q'),
            original_T_singleton_exact_mpf_tuple=self.T_tuple,
            anchored_G=lead._jet_summary(p['G'],t['G']),
            unmodified_log_shape_B=lead._jet_summary(p['B'],t['B']),
            normalized_modified_leading={n:lead._jet_summary(point[n],truth[n]) for n in point},
            physical_similarity_log_Utheta=lead._jet_summary(logup,logut+aj),
            log_geometry=dict(R_definition='110*exp(y)',rho_definition='110*Lambda*exp(y)',
                log_R=self._summary(c.ln(110)+yp,iv.ln(110)+y),
                log_rho=self._summary(c.ln(110)+c.mpf(o.amplitude.ll)+yp,iv.ln(110)+iv.mpf(o.amplitude.ll_iv)+y),
                radius_or_amplitude_exponentials_materialized=False),
            unmodified_six_current_amplitude_units={n:lead._jet_summary(sixp[n],sixt[n]) for n in sixp},
            finite_backward_kernels=kernels,
            finite_decay_rate_proofs=self._serializable_proof(rate_proofs),
            six_history_contract=dict(H='Phi*(H_over_Phi)',Kpsi='Phi*(Kpsi_over_Phi)',
                B='Phi^2*(B_over_Phi2)',C='Phi^2*(C_over_Phi2); undressed',Mpsi='Mpsi',Apsi='Apsi',
                finite_kernel='int_0^y exp(-k*t+m*B*(sigma(q)-sigma(q-t/T)))dt',
                all_nonzero_actual_R110_inlets_retained=True,infinite_window_replacement=False,
                kernel_centers_are_finite_constant_rate_proxies_with_full_errors=True,
                finite_N_corrections_separate=True,positive_order_n1_history_separate=True),
            selected_integer=dict(recipe=o.integer,original_N_not_reselected=True,
                log_N_exact_mpf_tuples=[list(v) for v in o.logN._mpi_],global_high_derivative_stress_heat_admission=False),
            source_cover=dict(exact_Z_cells=frames,local_source_cells=len(local),prefix_source_cells=len(prefix),
                original_fast_phase='frac(N*(Y+log(11/10)+y-hbB*sc/2))',fast_phase_not_chosen=True,
                physical_current_prefix_measure='T*dq once',all_incoming_five_corrections_retained=True),
            original_five_correction_history_bounds={n:[o._source_summary(v) for v in rows] for n,rows in hist.items()},
            positive_order_n1_own_history=dict(actual_terminal_primitive=self._terminal[Z],
                mean_u_definition='actual_terminal_primitive/(110*Lambda) * exp(-y)',
                q_definition='[-(1+delta)*Z*mean_u-(1-Z^2)*dZ(mean_u)]/(1-delta*Z^2)',
                actual_primitive_retained=True,f_and_u_exact_zero=True,exponential_not_materialized=True),
            small_error_proofs=self._serializable_proof(proof),
            semantics='Actual same-family functions enclosed about numerical proxies; no source midpoint or new fast phase. Taylor Zk=D_Z^k/k!.',
            limits=dict(logPhi_Psi_Z_order=3,own_mean_Z_order=2,Qpsi_Z_order=1,unmodified_six_Z_order=2,
                full_physical_vector_assembled=False,modified_all_six_merged=False,pressure_exported=False,
                resolved_signed_correction_integrals=False,resolved_fast_phase=False,total_radial_derivatives=False,
                mean_Z3_and_QZZ=False,later_charts=False,global_N_admission=False,global_field=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))

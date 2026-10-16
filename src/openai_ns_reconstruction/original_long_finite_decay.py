"""Extend original long values across representable ordinary log distances.

The accepted old reader stays unchanged and hash-verifiable. Direct ordinary
exponentials retain the whole incoming coefficient; truly enormous distances
still require the original combined-log cap proof.
"""
import hashlib
import json

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_long_numeric import OriginalLongNumeric
from .original_long_velocity import OriginalLongVelocity


DIRECT_LOG_DISTANCE_LIMIT='1e30'


class OriginalLongFiniteDecayNumeric(OriginalLongNumeric):
    """Same exact source and selected N, with ordinary-distance histories."""
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        name='src/openai_ns_reconstruction/original_long_finite_decay.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _six(self,inp,sigB,qp,q,yp,y,proof):
        c,iv,J=self.ctx,self.iv,self.J
        direct=_bounds(c,y)[1]<=c.mpf(DIRECT_LOG_DISTANCE_LIMIT)
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
            rate_proofs[k]=dict(m=m,rate_min=low,rate_max=high,radius=radius,
                global_sigma_inequality='abs(sigma(q))<=8*q=8*y/T',
                full_zeroth_log_inequality='-rate_max*y<=-k*y+m*B0*sigma(q)<=-rate_min*y',
                scalar_B0_factor_included_before_cap=True,
                inherited_decay_evaluation='direct finite exponential' if direct else 'combined-log enclosure')
            kernel=self.kernels.backward_kernel(iv,t['B'].truncate(2),self.T,y,m,low,high,proof)
            kp=-c.expm1(-c.mpf(k)*yp)/c.mpf(k) if direct else 1/c.mpf(k)
            kpoint=J.constant(c,kp,2);kernels[k]=self.leading._jet_summary(kpoint,kernel)
            for name in names:
                logs=sigB[1].truncate(2)*m-y*rate
                if direct:
                    # The exact finite log Taylor jet is representable here.
                    # Do not cap exp(-y) before multiplying an enormous
                    # normalized incoming coefficient or its derivatives.
                    it=initialt[name].truncate(2)*logs.exp()
                    ip=initialp[name].truncate(2)*(sigB[0].truncate(2)*m-c.mpf(k)*yp).exp()
                else:
                    it=self.kernels.normalized_amplitude_decay(iv,logs,low,high,y,proof,
                        initial=initialt[name].truncate(2))
                    ip=J.constant(c,0,2)
                sourcep=(2*p['values']['Psi'].truncate(2) if name=='Kpsi_over_Phi' else
                         J.constant(c,2 if name=='H_over_Phi' else 1,2))
                sourcet=(2*t['values']['Psi'].truncate(2) if name=='Kpsi_over_Phi' else
                         J.constant(iv,2 if name=='H_over_Phi' else 1,2))
                pp[name]=ip+sourcep*kpoint;tt[name]=it+sourcet*kernel
        for name,targetp,targett in (('Mpsi',p['values']['Psi'],t['values']['Psi']),
                ('Apsi',p['values']['Psi']*p['values']['Psi'],t['values']['Psi']*t['values']['Psi'])):
            dp=(p['six'][name]-targetp).truncate(2);dt=(t['six'][name]-targett).truncate(2)
            if direct:
                pp[name]=targetp.truncate(2)+dp*c.exp(-yp)
                tt[name]=targett.truncate(2)+dt*iv.exp(-y)
            else:
                errors=[]
                for k in range(3):
                    mag=self.leading._abs(dt[k]);log=None if _bounds(c,mag)[1]==0 else iv.ln(mag)-y
                    errors.append(self.leading._symmetric(self.owner._cap(log,proof,name+'_inlet_memory_Z'+str(k))))
                pp[name]=targetp.truncate(2);tt[name]=targett.truncate(2)+J(iv,errors)
        return pp,tt,kernels,rate_proofs

    def _evaluate(self,Z,qvalue,yvalue):
        result=super()._evaluate(Z,qvalue,yvalue)
        _,truth=self._pair_y(result['coordinate']['y'])
        result['finite_decay_numeric_extension']=dict(
            direct_log_distance_limit=DIRECT_LOG_DISTANCE_LIMIT,
            direct_finite_exponential_used=_bounds(self.ctx,truth)[1]<=self.ctx.mpf(DIRECT_LOG_DISTANCE_LIMIT),
            actual_original_T_N_cutoff_G_and_inlets_unchanged=True,
            old_absolute_cap_not_weakened=True,
            truly_enormous_distances_require_combined_log_proof=True)
        return result

    def _pair_y(self,row):return self.ctx.mpf(row['approximation']),self._restore(row)


class OriginalLongFiniteDecayVelocity(OriginalLongVelocity):
    """Physical value API using the extended original finite-long reader."""
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.numeric=OriginalLongFiniteDecayNumeric(precision=precision,
            amplitude_precision=amplitude_precision,directory=directory)
        self.amplitude=self.numeric.owner.amplitude;self.ctx,self.iv=self.amplitude.ctx,self.amplitude.iv
        self.J=self.numeric.J;self.family=self.numeric.family
        self.hashes=dict(self.numeric.hashes);self.code_hashes=dict(self.numeric.code_hashes)
        path=self.numeric.leading.core.directory/(PREFIX+'compliant_current_original_long_velocity.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:
            raise ValueError('Accepted physical original-long assembly required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed physical source: '+name)
            self.code_hashes[name]=digest
        for name,digest in accepted['input_hashes'].items():
            if self.hashes.get(name)!=digest:
                if hashlib.sha256((self.numeric.leading.core.directory/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed physical input: '+name)
                self.hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()

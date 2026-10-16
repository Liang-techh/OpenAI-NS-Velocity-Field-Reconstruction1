"""Physical reference/restoration values through the actual unpatched Rm."""
import copy
import hashlib
import json

from .original_core_numeric import ROOT,PREFIX
from .original_long_velocity import OriginalLongVelocity
from .original_reference_numeric import OriginalReferenceNumeric


class OriginalReferenceVelocity(OriginalLongVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.numeric=OriginalReferenceNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.amplitude=self.numeric.owner.amplitude;self.ctx,self.iv=self.amplitude.ctx,self.amplitude.iv
        self.J=self.numeric.J;self.family=self.numeric.family
        self.hashes=dict(self.numeric.hashes);self.code_hashes=dict(self.numeric.code_hashes)
        path=self.numeric.leading.core.directory/(PREFIX+'compliant_current_original_long_velocity.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:
            raise ValueError('Accepted physical long-chart assembly required')
        for name,digest in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed physical prerequisite: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        name='src/openai_ns_reconstruction/original_reference_velocity.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _assemble(self,packet,Z,log_tau,scaled_time,azimuth):
        result=super()._assemble(packet,Z,log_tau,scaled_time,azimuth)
        out=result['physical'];c,iv=self.ctx,self.iv
        base=packet['baseline_axial_histories'];alpha=self._pair(base['alpha']);gamma=self._pair(base['gamma1'])
        z,zi=c.mpf(Z),iv.mpf(Z);j,ji=self.amplitude.j,self.amplitude.j_iv
        velocity=(4*z+j*alpha[0],4*zi+ji*alpha[1])
        mean=(4*z+j*gamma[0],4*zi+ji*gamma[1])
        def radial(zz,delta,u,m):
            return (2*zz*u-(1-delta)*zz*m-4*(1-zz*zz))/(1-delta*zz*zz)
        qbase=(radial(z,self.amplitude.delta,velocity[0],mean[0]),
               radial(zi,self.amplitude.delta_iv,velocity[1],mean[1]))
        replacements={'leading_axis':velocity,'leading_radial_axis':qbase,'leading_mass_axis':mean}
        groups=[*out['cylindrical'].values(),out['streamfunction']]
        groups.append([out['coordinates']['r']])
        for group in groups:
            for row in group:
                if row['name'] in replacements:row['coefficient']=self._summary(replacements[row['name']])
                logs=[]
                for log in row['log_scale_terms']:
                    if log.get('kind')=='original_finite_long_y_log':
                        factor=log['multiplier']
                        for n,component in enumerate(packet['coordinate']['y_log_terms']):
                            pp,tt=self._pair(component)
                            item=self._summary((pp*c.mpf(factor),tt*iv.mpf(factor)))
                            item.update(kind='original_reference_geometry_log',component=n,multiplier=factor,
                                source_chart=packet['chart'],unit_offsets_kept_separate=True)
                            logs.append(item)
                    else:logs.append(log)
                row['log_scale_terms']=logs
        angle,anglei=c.mpf(azimuth),iv.mpf(azimuth)
        cosine,sine=(c.cos(angle),iv.cos(anglei)),(c.sin(angle),iv.sin(anglei))
        negsine=-sine[0],-sine[1];ur,ut,uz=(out['cylindrical'][n] for n in ('ur','utheta','uz'))
        out['Cartesian']=dict(u=[self._rotate(r,cosine) for r in ur]+[self._rotate(r,negsine) for r in ut],
            v=[self._rotate(r,sine) for r in ur]+[self._rotate(r,cosine) for r in ut],w=copy.deepcopy(uz))
        out['coordinates']['x']=self._rotate(out['coordinates']['r'],cosine)
        out['coordinates']['y']=self._rotate(out['coordinates']['r'],sine)
        out['native'].update(chart='actual_modified_'+packet['chart'],
            y_definition='T*(1-q)+10*logC*q+10*logP*q-8q for reference; 10logC+10logP+offset for restoration/postrestore',
            unit_radius_offsets_atomic=True)
        out['baseline_axial_history_contract']=packet['baseline_contract']
        out['baseline_coefficients']=dict(U=self._summary(velocity),M=self._summary(mean),Q=self._summary(qbase))
        out['n1_own_history']['radial_identity']='Q_J*lambda^(-1+2delta)*Lambda^-2/sqrt220*exp(-global_y/2); atomic source geometry'
        result['limits'].update(reference_restoration_postrestore_physical_values=True,actual_patch_connected=False)
        return result

    def reference(self,Z,q,**kwargs):return self._assemble(self.numeric.reference(Z,q),Z,**self._options(kwargs))
    def at_reference_gap(self,Z,gap,**kwargs):return self._assemble(self.numeric.at_reference_gap(Z,gap),Z,**self._options(kwargs))
    def restoration(self,Z,t,**kwargs):return self._assemble(self.numeric.restoration(Z,t),Z,**self._options(kwargs))
    def postrestore(self,Z,p,**kwargs):return self._assemble(self.numeric.postrestore(Z,p),Z,**self._options(kwargs))

    @staticmethod
    def _options(options):
        unexpected=set(options)-{'log_tau','scaled_time','azimuth'}
        if unexpected:raise TypeError('Unsupported physical options: '+str(unexpected))
        return dict(log_tau=options.get('log_tau','-1'),scaled_time=options.get('scaled_time'),azimuth=options.get('azimuth','0'))

    def evaluate(self,Z,q,**kwargs):return self.reference(Z,q,**kwargs)
    def at_log_radius(self,*args,**kwargs):raise ValueError('Use at_reference_gap for distance from actual Rsh')

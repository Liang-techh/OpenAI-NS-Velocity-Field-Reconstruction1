"""Original Rh/Rref physical velocity with separate finite-N and n1 histories."""
import hashlib

from .original_core_numeric import ROOT
from .original_Rm_patch_velocity import OriginalRmPatchVelocity
from .original_outer_reference_numeric import OriginalOuterReferenceNumeric


class OriginalOuterReferenceVelocity(OriginalRmPatchVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.numeric=OriginalOuterReferenceNumeric(patch=self.numeric)
        self.hashes.update(self.numeric.hashes);self.code_hashes.update(self.numeric.code_hashes)
        self.n1_support=dict(self.n1_support)
        self.n1_support.pop('lower_to_patch_order')
        self.n1_support.update(lower_to_reference_order='rminus<R110<Rm<Rh<=R<=Rref<Rw<Rc',
            reference_outer_order='log(Rw/Rref)=logPstar+1>0; reserved supports are beyond Rc',
            positive_log_Rw_over_Rref=self._summary((self.amplitude.lp+1,self.amplitude.lp_iv+1)))
        name='src/openai_ns_reconstruction/original_outer_reference_velocity.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def evaluate(self,Z,offset,*,log_tau='-1',scaled_time=None,azimuth='0'):
        packet=self.numeric.evaluate(Z,offset)
        result=super()._assemble(packet,Z,log_tau,scaled_time,azimuth)
        physical=result['physical'];native=physical['native']
        native.update(chart='actual_original_Rh_reference',reference_offset=offset,
            R_definition='110*exp(10logCstar+10logPstar+reference_offset)',
            x_definition='R/Rm=exp(6+reference_offset)')
        # Only the presentation metadata changes; the shared physical algebra
        # consumes this packet's actual log(R/110) atoms without an x lookup.
        def rewrite(obj):
            if isinstance(obj,dict):
                if obj.get('kind')=='original_Rm_geometry_y_log':
                    obj.update(kind='original_outer_reference_geometry_y_log',y_definition='log(R/110)=10logCstar+10logPstar+reference_offset')
                for value in obj.values():rewrite(value)
            elif isinstance(obj,list):
                for value in obj:rewrite(value)
        rewrite(physical)
        return result

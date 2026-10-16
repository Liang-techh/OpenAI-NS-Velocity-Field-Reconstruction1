"""Original compact temporal orders0/1 with actual reserved n1 mass recovery."""
import copy
import hashlib
import json

from .original_core_numeric import ROOT,PREFIX
from .original_core_velocity import OriginalCoreVelocity
from .original_compact_background_velocity import OriginalCompactBackgroundVelocity
from .original_reserved_n1_numeric import OriginalReservedN1Numeric


class OriginalCompactN1Velocity(OriginalCompactBackgroundVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,cells=256):
        super().__init__(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        path=self.numeric.directory/(PREFIX+'compliant_current_original_compact_numeric.json');accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:raise ValueError('Accepted actual compact background required')
        for n,h in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h:raise ValueError('Changed compact background implementation: '+n)
            self.code_hashes[n]=h
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self.reserved=OriginalReservedN1Numeric(self.numeric,cells=cells)
        self.hashes.update(self.reserved.hashes);self.code_hashes.update(self.reserved.code_hashes)
        for name in ('src/openai_ns_reconstruction/original_compact_n1_velocity.py','src/openai_ns_reconstruction/original_n1_physical_source.py'):
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def evaluate(self,Z,x,*,log_tau='-1',scaled_time=None,azimuth='0'):
        result=super().evaluate(Z,x,log_tau=log_tau,scaled_time=scaled_time,azimuth=azimuth)
        n1=self.reserved.evaluate(Z,x,self.numeric.evaluate(Z,x));c,iv=self.ctx,self.iv
        ell,de,_=OriginalCoreVelocity._time(self,Z,log_tau,scaled_time)
        y=[self._pair(v) for v in result['coordinate']['global_y_log_terms']]
        def scale(v,k):return v[0]*k,v[1]*k
        def logs(a,b):
            out=[self._summary(scale(v,a)) for v in ell if a]
            out.extend(self._summary(scale(v,b)) for v in de if b)
            return out
        def geometry(multiplier,local=True):
            rows=[self._summary(scale(v,c.mpf(multiplier))) for v in (y if local else y[:4])]
            rows.append(self._summary((c.log(110)*c.mpf(multiplier),iv.log(110)*iv.mpf(multiplier))))
            return rows
        def term(name,source,weights):
            return dict(name=name,temporal_order=1,coefficient=self._summary(self._pair(source)),log_scale_terms=weights)
        axial=term('n1_reserved_axial',n1['Uz1_times_Rc_C0_C1'][0],logs(-1,1)+geometry('-1',False))
        swirl=term('n1_reserved_swirl',n1['Utheta1_normalized_C0_C1'][0],logs(-1,1)+copy.deepcopy(n1['Utheta1_log_unit_atoms']))
        radial=term('n1_reserved_radial',n1['R_times_Q1'][0],logs(-1,2)+geometry('-.5')+
            [self._summary((-c.log(2)/2,-iv.log(2)/2))])
        stream=term('n1_reserved_streamfunction',n1['M1_C0_C1'][0],logs(1,1))
        p=result['physical'];p['cylindrical']['uz'].append(axial);p['cylindrical']['utheta'].append(swirl);p['cylindrical']['ur'].append(radial)
        p['streamfunction'].append(stream)
        angle,anglei=c.mpf(azimuth),iv.mpf(azimuth);cs=c.cos(angle),iv.cos(anglei);sn=c.sin(angle),iv.sin(anglei)
        p['Cartesian']=dict(u=[self._rotate(v,cs) for v in p['cylindrical']['ur']]+[self._rotate(v,(-sn[0],-sn[1])) for v in p['cylindrical']['utheta']],
            v=[self._rotate(v,sn) for v in p['cylindrical']['ur']]+[self._rotate(v,cs) for v in p['cylindrical']['utheta']],w=copy.deepcopy(p['cylindrical']['uz']))
        p.update(coefficient_orders_included=[0,1],compact_repaired_n1_field_included=True,
            n1_physical_source_units_already_applied=True,n1_extra_Lambda_squared_factor_applied=False,
            n1_swirl_is_Utheta1_not_F1=True,n1_mass_and_radial_from_same_cumulative_streamfunction=True,
            full_n1_pressure_included=False,n1_pressure_available=False,n1_numeric_Q1_Z_derivative_available=False)
        result['reserved_n1']=n1;result['limits'].update(physical_temporal_orders=[0,1],compact_n1_numeric_fields=True,
            compact_n1_mass_C1=True,compact_n1_Q1_value=True,compact_n1_Q1_Z=False,n1_pressure=False)
        return result

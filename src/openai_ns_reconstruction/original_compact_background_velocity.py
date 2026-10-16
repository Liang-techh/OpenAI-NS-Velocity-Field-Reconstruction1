"""Physical temporal-order-zero field in the original compact repair band.

The independently repaired n1 chart is deliberately absent here. It is not
continued by copying its pre-Rc homogeneous primitive.
"""
import copy
import hashlib

from .original_core_numeric import ROOT
from .original_core_velocity import OriginalCoreVelocity
from .original_long_velocity import OriginalLongVelocity
from .original_compact_numeric import OriginalCompactNumeric


class OriginalCompactBackgroundVelocity(OriginalLongVelocity):
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.numeric=OriginalCompactNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.amplitude=self.numeric.owner.amplitude;self.ctx,self.iv=self.amplitude.ctx,self.amplitude.iv
        self.J=self.numeric.J;self.family=self.numeric.family
        self.hashes=dict(self.numeric.hashes);self.code_hashes=dict(self.numeric.code_hashes)
        name='src/openai_ns_reconstruction/original_compact_background_velocity.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def evaluate(self,Z,x,*,log_tau='-1',scaled_time=None,azimuth='0'):
        p=self.numeric.evaluate(Z,x);c,iv=self.ctx,self.iv
        if type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)):raise ValueError('Finite decimal azimuth required')
        ell,de,time=OriginalCoreVelocity._time(self,Z,log_tau,scaled_time)
        y=[self._pair(r) for r in p['coordinate']['global_y_log_terms']]
        lp=self.amplitude.lp,self.amplitude.lp_iv
        def scale(pair,k):return pair[0]*k,pair[1]*k
        def logs(a,b,P=0,constant=None,yfactor=None):
            out=[self._summary(scale(v,a)) for v in ell if a]
            out.extend(self._summary(scale(v,b)) for v in de if b)
            if P:out.append(self._summary(scale(lp,P)))
            if constant is not None:out.append(self._summary(constant))
            if yfactor is not None:
                for i,v in enumerate(y):
                    row=self._summary(scale(v,c.mpf(yfactor)));row.update(kind='original_compact_geometry_y_log',component=i,multiplier=yfactor,
                        y_definition='log(R/110)=10logCstar+10logPstar+exp40+14+log(x)',unit_offsets_kept_separate=True)
                    out.append(row)
            return out
        def term(name,pair,weights):return dict(name=name,temporal_order=0,coefficient=self._summary(pair),log_scale_terms=weights)
        rows={n:self._pair(p['physical_similarity'][n][0]) for n in ('Uz','Utheta_over_Pstar','Ur_over_sqrt_R_over2','Mz_over_R')}
        radial=term('background_radial',rows['Ur_over_sqrt_R_over2'],logs(-1,0,constant=(c.log(55)/2,iv.log(55)/2),yfactor='.5'))
        swirl=term('background_swirl',rows['Utheta_over_Pstar'],logs(-1,-1,P=1))
        axial=term('background_axial',rows['Uz'],logs(-1,-1))
        stream=term('background_streamfunction',rows['Mz_over_R'],logs(1,-1,constant=(c.log(110),iv.log(110)),yfactor='1'))
        radius=term('physical_radius',(c.mpf(1),iv.mpf(1)),logs(1,0,constant=(c.log(220)/2,iv.log(220)/2),yfactor='.5'))
        angle,anglei=c.mpf(azimuth),iv.mpf(azimuth);cs=c.cos(angle),iv.cos(anglei);sn=c.sin(angle),iv.sin(anglei)
        pressure=term('background_pressure',self._pair(p['generic_modified_fields']['pressure'][0]),logs(-2,-2,P=2))
        physical=dict(source_family=self.family,time=time,native=dict(chart='actual_original_compact_Rc_2Rc',Z=Z,x=x,azimuth=azimuth,
                R_definition='110*exp(10logCstar+10logPstar+exp40+14)*x',physical_viscosity=1,original_N_and_construction_parameters_retained=True),
            coordinates=dict(r=radius,x=self._rotate(radius,cs),y=self._rotate(radius,sn),
                z=term('physical_z',(c.mpf(Z),iv.mpf(Z)),logs(1,-1)),t=dict(kind='one_minus_exp_of_log_tau',definition=time)),
            cylindrical=dict(ur=[radial],utheta=[swirl],uz=[axial]),
            Cartesian=dict(u=[self._rotate(radial,cs),self._rotate(swirl,(-sn[0],-sn[1]))],
                v=[self._rotate(radial,sn),self._rotate(swirl,cs)],w=[copy.deepcopy(axial)]),
            streamfunction=[stream],leading_pressure=pressure,original_selected_integer=p['selected_integer'],coefficient_orders_included=[0],
            shared_streamfunction_meridional_recovery=True,actual_selected_N_errors_included=True,
            absolute_physical_exponentials_materialized=False,native_coordinates_required=True,
            compact_repaired_n1_field_included=False,preRc_n1_homogeneous_history_not_extrapolated=True,
            direct_physical_inverse_available=False,physical_gradients_or_vorticity_available=False,
            leading_pressure_includes_same_P0_once=True,full_n1_pressure_included=False,
            whole_global_temporal_sum_admitted=False,global_field_available=False,genuine_n_ge_2_temporal_recursion=False,
            full_corrected_NS_solution=False,value_semantics='sum coefficient*exp(sum log_scale_terms), with full directed errors')
        result=copy.deepcopy(p);result['physical']=physical
        result['limits'].update(physical_Cartesian_assembly=True,physical_temporal_orders=[0],physical_pressure_leading_only=True)
        return result

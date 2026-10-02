"""Independent O.3 pure-power continuation immediately left of Rp.

Canonical incoming functions are rebuilt from their source constants,
not copied pulse derivative grids. Original pressure is evaluated from the
analytic datum. Local y=log(R/Rp) lies in [-1,0], wholly inside the actual
long power buffer. All derivatives have leading-profile scope only.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import transport_mixed
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


class CompliantPowerInletC4:
    def __init__(self):
        self.ctx=c=MPIntervalContext(); c.dps=240
        name=PREFIX+'compliant_pulse_flat_comparison_check.json'
        receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['same_nonzero_interface_histories_preserved']:
            raise ValueError('Admitted same-history flat differences required')
        self.hashes=dict(receipt['input_hashes'])
        for source,digest in self.hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Power inlet source changed: '+source)
        self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.family=receipt['actual_five_defect_family_sha256']; self.source=receipt['implicit_source_sha256']
        fifth=json.loads((HERE/(PREFIX+'compliant_fifth_axial_jets.json')).read_bytes())
        pulse=json.loads((HERE/(PREFIX+'compliant_pulse_mixed_C4.json')).read_bytes())
        buffer_name=PREFIX+'compliant_outer_buffer.json'
        buffer=json.loads((HERE/buffer_name).read_bytes())
        if any(r['actual_five_defect_family_sha256']!=self.family or r['implicit_source_sha256']!=self.source for r in (fifth,pulse,buffer)):
            raise ValueError('Power inlet source/family mismatch')
        self.constants={k:read_interval(c,v) for k,v in
            fifth['whole_Z']['selected']['incoming']['Z_independent_constant_definitions'].items()
            if isinstance(v,dict) and 'lower' in v}
        self.mu=read_interval(c,pulse['selected_mu']); self.delta=read_interval(c,pulse['selected_delta'])
        self.rate=1-self.mu; self.prate=1+2*self.mu
        if endpoints(self.constants['Tw'])[0]<=1:raise ValueError('Unit left neighborhood is outside the actual power buffer')
        sample=buffer['samples'][-1]
        if not sample['pulse_inlet'] or endpoints(read_interval(c,sample['Z']))!=(mp.mpf('.5'),mp.mpf('.5')):
            raise ValueError('Source-identified final O3 power sample required')
        # Exact canonical q^-1 and q^-2 shapes transport unchanged through
        # O3. Recover their Z-independent constants, not a fitted profile.
        q=c.mpf('1.25')
        self.inlet_H=read_interval(c,sample['Mtheta_over_sqrt2_R_3half_Pstar'][0])*q
        self.inlet_P=read_interval(c,sample['Mp_over_Pstar_squared'][0])*q*q
        self.Xp=self.inlet_H/self.constants['U']
        self.datum=CompliantPressureDatum('40',precision=160)
        if self.datum.source_sha!=self.source or self.datum.datum_sha!=buffer['datum_enclosure_sha256']:
            raise ValueError('Original analytic pressure datum mismatch')
        self.hashes[buffer_name]=hashlib.sha256((HERE/buffer_name).read_bytes()).hexdigest()
        self.hashes.update(self.datum.input_hashes)
        self.radial=CompliantPulseRadialC4.__new__(CompliantPulseRadialC4)
        self.radial.ctx=c; self.radial.mu=self.mu; self.radial.delta=self.delta; self.radial.prate=self.prate
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def incoming(self,Z):
        c=self.ctx; Z=c.mpf(Z)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Z in[-1,1] required')
        z=IntervalTaylor(c,[Z,1,0,0,0,0])
        q=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]); invq=q.reciprocal(); k=self.constants
        u=invq*k['U']; m=[z*q*k[key] for key in ('C1','C2')]
        e=z*z*q*q*k['C_E']+k['C0']
        raw=self.datum.normalized_jets(Z,5)['normalized_pressure_coefficients']
        p0=IntervalTaylor(c,[c.mpf(endpoints(v)) for v in raw])
        return dict(u=u,m1=m[0],m2=m[1],energy=e,X=IntervalTaylor.constant(c,self.Xp,5),
            P0=p0,Mp=invq*invq*self.inlet_P)

    def power(self,Z,distance):
        c=self.ctx; d=c.mpf(distance)
        if endpoints(d)[0]<0 or endpoints(d)[1]>1:raise ValueError('Left local distance in[0,1] required')
        data=self.incoming(Z); u=data['u']; zero=u*0
        # y=-d. The primitive differences are cancellation-free, including
        # the actual sub-precision positive mu. All histories are retained.
        m1=data['m1']*c.exp((c.mpf('.5')-self.mu)*d)
        m2=data['m2']*c.exp((c.mpf('.5')-2*self.mu)*d)
        e=data['energy']*c.exp(-2*self.mu*d)+decay_integral(c,2*self.mu,d)/2
        X=IntervalTaylor.constant(c,1/self.rate,5)+(data['X']-1/self.rate)*c.exp(self.rate*d)
        kernel=c.exp(self.prate*d)*decay_integral(c,self.prate,d)
        Mp=data['Mp']-u*u*(kernel/2)
        pressure=dict(Mp_over_Pstar_squared=Mp,P0_over_Pstar_squared=data['P0'],
            P_over_Pstar_squared=Mp+data['P0'],P_y_over_Pstar_squared=u*u*(c.exp(self.prate*d)/2),
            original_analytic_datum_preserved=True)
        point=dict(Mz_over_R_Utheta=m1,Mtheta_z_over_sqrt2_R_3half_Utheta_squared=m2,
            Mtheta_over_sqrt2_R_3half_Utheta=X,Mztheta_over_R_Utheta_squared=e,pressure=pressure)
        mixed=transport_mixed(self.radial,c.mpf(Z),point,[zero for _ in range(5)],u)
        point.update(mixed,Z=c.mpf(Z),left_distance=d,y=-d,
            canonical_incoming_Taylor=data,
            velocity_common_radial_log_factors=dict(theta_axial=(c.mpf('.5')+self.mu)*d,radial=self.mu*d),
            coordinate_scope='y=log(R/Rp), paper Z; mixed derivatives of physical profile components with common radial factors separated',
            radial_velocity_not_reset_when_axial_input_is_zero=True,
            independently_evaluated_O3_power_high_mixed_derivatives=True,
            source_shapes_exact_not_fitted=True,full_pulse_C4_installed=False,
            full_outer_C4_certified=False,physical_energy_integral_certified=False,
            whole_outer_cone_certified=False,temporal_recursion=False)
        return point

    def report(self):
        with mp.workdps(270):
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                datum_enclosure_sha256=self.datum.datum_sha,
                samples=[self.power(z,d) for z in ('-1','0','.5','1') for d in ('1','.01','0')],
                whole_Z_left_box=self.power([-1,1],[0,1]),whole_Z_inlet=self.power([-1,1],0),
                selected_mu=self.mu,selected_delta=self.delta,
                exact_constant_definitions=dict(Hp='q(Z)*O3.buffer.power(Z,1).Mtheta',
                    Pin='q(Z)^2*O3.buffer.power(Z,1).Mp',Xp='Hp/U',
                    incoming='u=U/q; mi=Ci*Z*q; e=C0+C_E*Z^2*q^2; Mp=Pin/q^2'),
                independently_evaluated_O3_power_high_mixed_derivatives=True,
                two_sided_O3_pulse_join_certified=False,
                full_pulse_C4_installed=False,full_outer_C4_certified=False,
                physical_energy_integral_certified=False,whole_outer_cone_certified=False,
                temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantPowerInletC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Original O3 left power neighborhood: independent canonical C5 primitives and mixed velocity/pressure order<=4 generated',flush=True)
    return result


if __name__=='__main__':run()

"""Implicit derivatives of the continuous axial correction equations.

The matrix/pulse/Gram atoms and mu are fixed in Z for this source schedule.
Callers must supply real incoming-row and energy-target derivatives. The
standalone diagnostic uses a declared input direction, not measured Z jets.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log


class ContinuousAxialAlgebra:
    def __init__(self,source,receipt):
        self.precision=receipt['algebra_precision']
        with mp.workdps(self.precision):
            self.mu=mp.mpf(receipt['input_mu'])
            self.M=[[mp.mpf(x) for x in row] for row in receipt['linear_matrix']]
            lam1=mp.mpf('.5')-self.mu;lam2=mp.mpf('.5')-2*self.mu
            B1=self.M[0][0]*mp.exp(3*lam1)
            B2=self.M[1][0]*mp.exp(3*lam2)
            self.det=-2*mp.exp(-2+6*self.mu)*mp.sinh(self.mu)*B1*B2
            self.p=[mp.exp(mp.mpf(x['log_normalized_pulse_integral'])) for x in receipt['pulse_rows']]
            inputs=source['axial'].get('linear_rhs_inputs')
            base=(inputs['base'] if inputs is not None else
                  [source['incoming']['row_normalization'][key]
                   for key in ('scaled_base_m1','scaled_base_m2')])
            self.base=[from_signed_log(x) for x in base]
            target=source.get('energy_target',source['axial'].get('energy_target'))
            if target is None:raise ValueError('Actual energy target must be retained')
            self.target=mp.mpf(target)
            self.K=[mp.mpf(x) for x in receipt['K_bump']]
            self.Kp=mp.mpf(receipt['K_p'])
            self.v=self.affine(self.p)

    def affine(self,rhs):
        return [(-rhs[0]*self.M[1][1]+self.M[0][1]*rhs[1])/self.det,
                (-self.M[0][0]*rhs[1]+rhs[0]*self.M[1][0])/self.det]

    def solve(self,base=None,target=None):
        with mp.workdps(self.precision):
            base=self.base if base is None else base
            target=self.target if target is None else target
            u=self.affine(base);v=self.v
            quadratic=self.Kp+self.mu*sum(k*x*x for k,x in zip(self.K,v))
            linear=2*self.mu*sum(k*x*y for k,x,y in zip(self.K,u,v))
            constant=self.mu*sum(k*x*x for k,x in zip(self.K,u))-target
            disc=linear*linear-4*quadratic*constant
            if constant>=0 or disc<=0:raise ArithmeticError('No positive energy branch')
            a=-2*constant/(linear+mp.sqrt(disc))
            return a,[x+a*y for x,y in zip(u,v)]

    def tangent(self,base_Z,target_Z):
        """Differentiate M c+b+a p=0, Kp a^2+mu sum(K c^2)=target."""
        with mp.workdps(self.precision):
            a,c=self.solve();u_Z=self.affine(base_Z)
            denominator=2*self.Kp*a+2*self.mu*sum(k*x*y for k,x,y in zip(self.K,c,self.v))
            if denominator<=0:raise ArithmeticError('Energy branch derivative is singular')
            a_Z=(target_Z-2*self.mu*sum(k*x*y for k,x,y in zip(self.K,c,u_Z)))/denominator
            c_Z=[x+y*a_Z for x,y in zip(u_Z,self.v)]
            linear_replay=[]
            for i in range(2):
                terms=[self.M[i][j]*c_Z[j] for j in range(2)]+[base_Z[i],a_Z*self.p[i]]
                scale=sum(abs(x) for x in terms)
                linear_replay.append(mp.nstr(abs(sum(terms))/scale,40) if scale else '0')
            energy_terms=[2*self.Kp*a*a_Z,
                2*self.mu*sum(k*x*y for k,x,y in zip(self.K,c,c_Z)),-target_Z]
            scale=sum(abs(x) for x in energy_terms)
            return dict(a_Z=a_Z,c_Z=c_Z,linear_tangent_relative_replay=linear_replay,
                energy_tangent_relative_replay=mp.nstr(abs(sum(energy_terms))/scale,40) if scale else '0')


def run():
    folder=Path(__file__).parent
    source=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    receipt=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    algebra=ContinuousAxialAlgebra(source,receipt)
    with mp.workdps(algebra.precision):
        direction=algebra.base;target_direction=algebra.target/10
        jet=algebra.tangent(direction,target_direction)
        refinements=[]
        for h in (mp.mpf('1e-3'),mp.mpf('5e-4')):
            values={i:algebra.solve([b+i*h*d for b,d in zip(algebra.base,direction)],
                                   algebra.target+i*h*target_direction) for i in (-2,-1,1,2)}
            numerical=[(values[-2][0]-8*values[-1][0]+8*values[1][0]-values[2][0])/(12*h)]
            numerical += [(values[-2][1][j]-8*values[-1][1][j]+8*values[1][1][j]-values[2][1][j])/(12*h) for j in range(2)]
            analytic=[jet['a_Z']]+jet['c_Z']
            refinements.append(dict(step=mp.nstr(h,20),relative_errors=[mp.nstr(abs(x/y-1),40) for x,y in zip(numerical,analytic)]))
        report=dict(input_direction='base derivative = actual base; target derivative = actual target / 10; a declared algebra direction, not physical Z',
            a_directional_derivative=signed_log(jet['a_Z'],algebra.precision),
            c_directional_derivatives=[signed_log(x,algebra.precision) for x in jet['c_Z']],
            linear_tangent_relative_replay=jet['linear_tangent_relative_replay'],
            energy_tangent_relative_replay=jet['energy_tangent_relative_replay'],
            independent_directional_refinement=refinements,
            actual_incoming_Z_jets_installed=False,finite_energy_certified=False,
            fixed_continuous_atoms=True,quadrature_enclosure_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('linear_tangent_relative_replay',
        'energy_tangent_relative_replay','independent_directional_refinement')}),flush=True)
    return report


if __name__=='__main__':run()

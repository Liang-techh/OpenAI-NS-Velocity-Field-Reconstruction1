"""Finite normalized core plus same-source analytic tails at scaled exit.

Local axial family only. Provides callable jets for subsequent comparison.
No physical amplitude selection, repaired inlet or full residual claim.
"""
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_candidate_gauge_core import _unpack_rows,_unpack_vector
from lei_ren_part1_paper_candidate_shared_inlet import finite_moments
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
STEM='lei_ren_part1_paper_shared_interval_core_Z049_Z051'

class LogarithmicCoreExit:
    def __init__(self):
        self.receipt=json.loads((HERE/(STEM+'.json')).read_bytes())
        if not self.receipt['completed_target']:raise ValueError('Fresh selected-degree finite core unfinished')
        state_path=HERE/self.receipt['state_file']
        if hashlib.sha256(state_path.read_bytes()).hexdigest()!=self.receipt['state_sha256']:
            raise ValueError('core state changed')
        self.state=json.loads(state_path.read_bytes())
        seed_name='lei_ren_part1_paper_shared_core_seed_check.json'
        seed=json.loads((HERE/seed_name).read_bytes())
        if not seed['seed_acceptance_passed']:raise ValueError('seed acceptance missing')
        fixed_sha=hashlib.sha256(json.dumps(self.state['fixed'],sort_keys=True).encode()).hexdigest()
        if fixed_sha!=seed['fixed_seed_sha256']:raise ValueError('accepted seed changed')
        for n,d in seed['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:raise ValueError('seed dependency changed: '+n)
        for n,d in self.receipt['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:raise ValueError('dependency changed: '+n)
        self.tail=json.loads((HERE/'lei_ren_part1_paper_shared_core_tail_admission.json').read_bytes())
        identity=self.state['identity']
        if self.tail['implicit_source_sha256']!=identity['implicit_source_sha256'] or self.tail['datum_enclosure_sha256']!=identity['datum_enclosure_sha256']:
            raise ValueError('tail source changed')
        if self.tail['radial_degree']!=identity['target_radial_degree']:raise ValueError('tail degree changed')
        if not identity['source_j_eta_tol_relation_verified'] or self.tail['analytic_core_family_sha256']!=identity['analytic_core_family_sha256'] or seed['analytic_core_family_sha256']!=identity['analytic_core_family_sha256']:
            raise ValueError('New-j core/tail/seed family mismatch')
        self.ctx=c=MPIntervalContext();c.dps=identity['precision']
        with mp.workdps(c.dps+40):
            self.rows={n:_unpack_rows(c,rs) for n,rs in self.state['rows'].items()}
            self.fixed={n:_unpack_vector(c,rs) for n,rs in self.state['fixed'].items()}
        self.input_hashes={**self.receipt['input_hashes'],self.receipt['state_file']:self.receipt['state_sha256'],
            STEM+'.json':hashlib.sha256((HERE/(STEM+'.json')).read_bytes()).hexdigest(),
            seed_name:hashlib.sha256((HERE/seed_name).read_bytes()).hexdigest(),
            'lei_ren_part1_paper_candidate_shared_inlet.py':hashlib.sha256((HERE/'lei_ren_part1_paper_candidate_shared_inlet.py').read_bytes()).hexdigest(),
            Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

    def mixed(self,r,i,k):
        c=self.ctx;r=c.mpf(r)
        if i<0 or k<0 or i+k>3:raise ValueError('mixed derivative total order <=3 required')
        if endpoints(r)[0]<0 or endpoints(r)[1]>endpoints(c.mpf('4.1'))[1]:raise ValueError('scaled radius outside [0,4.1]')
        with mp.workdps(c.dps+40):
            finite={name:sum((row[k]*math.factorial(k)*
                (math.factorial(n)//math.factorial(n-i))*r**(n-i)
                for n,row in enumerate(rows) if n>=i),c.mpf(0))
                for name,rows in self.rows.items()}
            tails=next(row for row in self.tail['tail_rows'] if row['scaled_radial_order']==i and row['axial_order']==k)
            full={}
            for name,key in (('A','Phi_tail'),('Uz','Uz_tail')):
                b=read_interval(c,tails[key]);upper=endpoints(b)[1]
                full[name]=finite[name]+c.mpf([-upper,upper])
            return dict(finite=finite,full_Phi=full['A'],full_Uz=full['Uz'],
                Phi_tail=read_interval(c,tails['Phi_tail']),Uz_tail=read_interval(c,tails['Uz_tail']))

    def core_moments(self,r):
        """Full scaled moment and epsilon-P jets through axial order two.

        Uniform value/axial tail bounds are integrated over [0,r]. They
        bound each ordinary Taylor coefficient, not one midpoint sample.
        """
        c=self.ctx
        with mp.workdps(c.dps+40):
            r=c.mpf(r)
            if endpoints(r)[0]<0 or endpoints(r)[1]>endpoints(c.mpf('4.1'))[1]:
                raise ValueError('radius outside analytic tail domain')
            jet=lambda values:IntervalTaylor(c,values[:3])
            phi=[jet(row) for row in self.rows['A']]
            u=[jet(row) for row in self.rows['Uz']]
            finite=finite_moments(phi,u,r)
            def absjet(rows):
                return IntervalTaylor(c,[sum((c.mpf(max(abs(v) for v in endpoints(row[k])))*r**n
                    for n,row in enumerate(rows)),c.mpf(0)) for k in range(3)])
            def tailjet(name):
                out=[]
                for k in range(3):
                    record=next(t for t in self.tail['tail_rows'] if t['scaled_radial_order']==0 and t['axial_order']==k)
                    bound=read_interval(c,record[name])/math.factorial(k)
                    out.append(c.mpf([0,endpoints(bound)[1]]))
                return IntervalTaylor(c,out)
            bp,bu=absjet(self.rows['A']),absjet(self.rows['Uz'])
            ep,eu=tailjet('Phi_tail'),tailjet('Uz_tail')
            pp=bp*ep*2+ep*ep;pu=bp*eu+bu*ep+ep*eu;uu=bu*eu*2+eu*eu
            errors=dict(theta=ep*r*r,z=eu*r,theta_z=pu*r*r,p=pp*r,
                u_squared=uu*r,weighted_phi_squared=pp*r*r/2)
            def symmetric(e):
                return IntervalTaylor(c,[c.mpf([-endpoints(x)[1],endpoints(x)[1]]) for x in e.coefficients])
            full={name:finite[name]+symmetric(errors[name]) for name in finite}
            pressure=jet(self.fixed['P0_Z_taylor'])+jet(self.fixed['S_Z_taylor'])*full['p']
            return dict(full_moment_coefficients={n:list(v.coefficients) for n,v in full.items()},
                moment_error_coefficient_bounds={n:list(v.coefficients) for n,v in errors.items()},
                restored_epsilon_pressure_coefficients=list(pressure.coefficients),
                pressure_identity='epsilon P = epsilon P0 + epsilon^2 F0^2 integral_0^r Phi^2 ds',
                moment_units='normalized r=Lambda R integrals; amplitude and Lambda factors must be restored',
                terminal_five_moment_constraints_solved=False)

    def report(self):
        c=self.ctx
        with mp.workdps(c.dps+40):
            cases=[]
            for r in ('4','4.1'):
                for total in range(4):
                    for i in range(total+1):
                        k=total-i
                        cases.append(dict(scaled_radius=r,scaled_radial_order=i,axial_order=k,**self.mixed(r,i,k)))
            positive_exits={r:endpoints(next(x for x in cases if x['scaled_radius']==r and x['scaled_radial_order']==0 and x['axial_order']==0)['full_Phi'])[0]>0 for r in ('4','4.1')}
            return dict(input_hashes=self.input_hashes,identity=self.state['identity'],
                normalized_swirl_shape_positive_at_exits=positive_exits,
                local_axial_family=['.49','.51'],finite_radial_degree=self.state['identity']['target_radial_degree'],mixed_C3_enclosures=cases,
                scaled_core_moments={r:self.core_moments(r) for r in ('4','4.1')},
                finite_and_same_source_tail_combined=True,
                saved_pressure_rows_only_finite=True,
                analytic_pressure_primitive_with_controlled_moment_tail_available=True,
                physical_swirl_amplitude_remains_implicit=True,
                actual_velocity_inlet_or_five_moment_repair_completed=False,
                physical_R_derivatives_certified=False,temporal_recursion=False,
                full_NS_residual_certified=False)

def run():
    e=LogarithmicCoreExit();result=e.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    print('Local Md40 finite+analytic-tail exit enclosures:',len(result['mixed_C3_enclosures']))
    return result

if __name__=='__main__':run()

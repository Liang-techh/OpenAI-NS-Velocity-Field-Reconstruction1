"""Perturbed interval view of the completed local 144-order core.

Each view is old coefficient enclosure plus a proved positive sensitivity
cap. This is an enclosure of the NEW analytic solution, not its recomputed
point coefficients and not a changed source label on the old state file.
"""
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_candidate_gauge_core import _unpack_rows, _pack_vector
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
STATE='lei_ren_part1_paper_shared_interval_core_Z049_Z051_state.json'
CORE='lei_ren_part1_paper_compliant_core_transfer.json'
CHECK='lei_ren_part1_paper_compliant_core_transfer_check.json'


class CompliantFiniteCore:
    def __init__(self):
        self.ctx=c=MPIntervalContext();c.dps=160
        self.major=json.loads((HERE/CORE).read_bytes())
        self.check=json.loads((HERE/CHECK).read_bytes())
        self.state=json.loads((HERE/STATE).read_bytes())
        self.hashes={}
        for record in (self.major,self.check,self.state):
            for name,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Finite transfer dependency changed: '+name)
                self.hashes[name]=digest
        for name in (CORE,CHECK,STATE,Path(__file__).name):
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        tr=self.major['core_transfer'];identity=self.state['identity']
        if (not self.check['all_passed'] or not tr['analytic_core_perturbation_transfer_proved']
            or identity['implicit_source_sha256']!=tr['old_implicit_source_sha256']
            or identity['analytic_core_family_sha256']!=self.major['analytic_core_family_sha256']):
            raise ValueError('Old/new analytic source binding missing')
        self.degree=identity['target_radial_degree']
        if (self.state['seed_metadata']['logC_definition']!=self.major['logC_definition']
            or self.state['seed_metadata']['implicit_positive_F0_definition']!=
                'exp(-logC-Lambda G), G primitive anchored at H root in [-j,0]'
            or identity['required_j']!=self.major['required_j']):
            raise ValueError('Shared F0 definition, primitive branch or j differs')
        if self.degree!=144 or self.state['completed_radial_order']!=self.degree:
            raise ValueError('Expected completed 144-order local core')
        if identity['center_Z']!=['.49','.51'] or identity['row_units']!={
            'A':'physical_A[n]/Lambda**n', 'Uz':'physical_Uz[n]/Lambda**n',
            'P':'epsilon*physical_P[n]/Lambda**n'}:
            raise ValueError('Coefficient domain/units changed')
        self.rows={name:_unpack_rows(c,rs) for name,rs in self.state['rows'].items()}
        if set(self.rows)!={'A','Uz','P'} or any(
            len(rs)!=self.degree+1 or any(len(row)!=self.degree+4-n for n,row in enumerate(rs))
            for rs in self.rows.values()):
            raise ValueError('Old coefficient layout incomplete')
        tube=json.loads((HERE/'lei_ren_part1_paper_shared_analytic_tube.json').read_bytes())
        with mp.workdps(c.dps+40):
            self.h=read_interval(c,tube['Xh_parameter'])
            self.eps=read_interval(c,self.major['epsilon'])
            oldeps=read_interval(c,identity['epsilon'])
            if endpoints(self.eps)!=endpoints(oldeps):
                raise ValueError('Old/new Lambda normalization differs')
            self.norms=dict(A=read_interval(c,tr['solution_difference_Xh_upper']),
                Uz=read_interval(c,tr['physical_Uz_difference_Xh_upper']),
                P=self.eps*read_interval(c,tr['physical_pressure_difference_Xh_upper']))
            self.axis_pressure_cap=self.eps*read_interval(c,self.major['pressure_perturbation']['physical_pressure_difference_abs_upper'])

    def cap(self,field,n,m):
        if field not in self.rows or not isinstance(n,int) or not isinstance(m,int):
            raise ValueError('Field and integer coefficient indices required')
        if n<0 or n>self.degree or m<0 or m>=len(self.rows[field][n]):
            raise ValueError('Coefficient outside stored finite jet layout')
        c=self.ctx
        if n==0:
            # Same exact axis velocity/swirl; only pressure value changes.
            return self.axis_pressure_cap if field=='P' and m==0 else c.mpf(0)
        with mp.workdps(c.dps+40):
            weight=c.mpf(math.comb(n+m,m))/(20**n*self.h**m*(n+1)**2*(m+1)**2)
            return self.norms[field]*weight

    def coefficient(self,field,n,m):
        c=self.ctx
        with mp.workdps(c.dps+40):
            b=endpoints(self.cap(field,n,m))[1]
            return self.rows[field][n][m]+c.mpf([-b,b])

    def row(self,field,n):
        return [self.coefficient(field,n,m) for m in range(len(self.rows[field][n]))]

    def report(self):
        c=self.ctx;tr=self.major['core_transfer']
        with mp.workdps(c.dps+40):
            samples=[]
            for field in self.rows:
                for n in (0,1,72,144):
                    for m in (0,3):
                        value=self.coefficient(field,n,m)
                        old=self.rows[field][n][m]
                        if not endpoints(value)[0]<=endpoints(old)[0]<=endpoints(old)[1]<=endpoints(value)[1]:
                            raise AssertionError('Perturbed coefficient view fails to contain prior enclosure')
                        samples.append(dict(field=field,radial_order=n,axial_order=m,
                            old_enclosure=_pack_vector([old])[0],
                            positive_difference_cap=_pack_vector([self.cap(field,n,m)])[0],
                            new_enclosure=_pack_vector([value])[0]))
            counts={name:sum(map(len,rs)) for name,rs in self.rows.items()}
            return dict(old_source_sha256=tr['old_implicit_source_sha256'],
                new_source_sha256=tr['new_implicit_source_sha256'],
                new_datum_enclosure_sha256=self.major['datum_enclosure_sha256'],
                analytic_core_family_sha256=self.major['analytic_core_family_sha256'],
                local_axis_centers=['.49','.51'],radial_degree=144,
                coefficient_count_per_field=counts,total_supported_coefficients=sum(counts.values()),
                normalized_difference_norms=self.norms,axis_scaled_pressure_difference_cap=self.axis_pressure_cap,
                coefficient_cap_formula='norm[field]*binomial(n+m,m)/(20^n*h^m*(n+1)^2*(m+1)^2)',
                source_state_unchanged=True,representation='callable old interval + symmetric sensitivity cap',
                caps_retained_separately_even_when_old_interval_width_dominates=True,
                same_axis_A_Uz_and_pressure_derivatives_preserved=True,
                new_finite_coefficient_enclosure_view_available=True,
                new_point_coefficients_recomputed=False,old_coefficients_relabelled=False,
                sample_enclosures=samples,
                full_axis_point_evaluator_available=False,
                downstream_exit_moment_cone_transfer_completed=False,
                actual_heat_pressure_change_included=False,temporal_recursion=False,
                input_hashes=self.hashes)


def run():
    core=CompliantFiniteCore();result=core.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    print('New local finite core enclosure view:',result['total_supported_coefficients'],
          'coefficients over 144 radial orders; old state preserved',flush=True)
    return result


if __name__=='__main__':run()

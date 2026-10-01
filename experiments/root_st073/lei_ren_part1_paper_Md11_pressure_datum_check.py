"""Independent exact stage, schedule, axial-series and datum dependency checks."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_Md11_pressure_datum import pressure_jets,BETA2
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_candidate_pressure_axis_jets import _interval_from_exact
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    here=Path(__file__).parent;name='lei_ren_part1_paper_Md11_pressure_datum.json'
    datum=json.loads((here/name).read_bytes());c=MPIntervalContext();c.dps=300
    for dep,digest in datum['input_hashes'].items():
        if hashlib.sha256((here/dep).read_bytes()).hexdigest()!=digest:raise AssertionError('changed dependency: '+dep)
    inputs=datum['new_schedule_inputs']
    sha=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    if sha!=datum['new_schedule_sha256']:raise AssertionError('new schedule digest mismatch')
    if sha=='736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4':
        raise AssertionError('old schedule relabeled as new')
    contains=0
    with mp.workdps(340):
        md=mp.mpf(inputs['Md']);lp=mp.mpf(inputs['logPstar']);yd=mp.exp(md)+11
        if not md>1 or not lp>yd-1:raise AssertionError('explicit paper parameter inequality failed')
        # Independent exact angular pressure primitive on [1,yd]. Md affects
        # the end but not the slope -1/2 or the primitive itself.
        # The enclosures concern the stored Decimal schedule. Its exp(Md)
        # rounding is explicitly outside scope; do not demand that its
        # 300-digit integration interval enclose an unrounded new parameter.
        stored_yd=mp.mpf(datum['stored_y_d'])
        exact=mp.mpf('.5')*(mp.exp(mp.mpf('-.4'))-mp.exp(mp.mpf('.6')-stored_yd))
        lo,hi=endpoints(read_interval(c,datum['stages']['axial_turnoff']['mass']))
        if not lo<=exact<=hi:raise AssertionError('exact axial-turnoff pressure mass excluded')
        contains+=1
        if len(datum['stages'])!=14 or sum(r['beta']==2 for r in datum['stages'].values())!=6:
            raise AssertionError('wrong full-stage partition')
        sum2=c.mpf(0);sum0=c.mpf(0)
        for stage,row in datum['stages'].items():
            mass=read_interval(c,row['mass'])
            if endpoints(mass)[0]<0:raise AssertionError('negative pressure mass')
            if row['beta']==2:
                if stage not in BETA2:raise AssertionError('wrong beta2 stage')
                sum2+=mass
            elif row['beta']==0:sum0+=mass
        for field,total in (('fixed_beta2_mass_upper_normalized',sum2),('fixed_beta0_mass_upper_normalized',sum0)):
            saved=_interval_from_exact(c,datum[field])
            if saved._mpi_!=total._mpi_:raise AssertionError('mass sum inconsistent')
        # Independent exact q^-2 axis Taylor series: even coefficient
        # (-1)^k(k+1), odd coefficients zero. Use contained midpoint fixed
        # masses, with zero flatten term (included in its positive envelope).
        mass2=_interval_from_exact(c,datum['fixed_beta2_mass_upper_normalized'])
        mass0=_interval_from_exact(c,datum['fixed_beta0_mass_upper_normalized'])
        m2=sum(endpoints(mass2))/2;m0=sum(endpoints(mass0))/2;scale=mp.exp(2*lp)
        jets=pressure_jets(c,0,12,datum)
        for n,box in enumerate(jets['physical_pressure_coefficients']):
            a=mp.mpf(0) if n%2 else mp.mpf((-1)**(n//2)*(n//2+1))
            value=-scale*(m2*a+(m0 if n==0 else 0))
            lo,hi=endpoints(box)
            if not lo<=value<=hi:raise AssertionError(('exact q-series coefficient excluded',n))
            contains+=1
        old=json.loads((here/'lei_ren_part1_paper_candidate_pressure_axis_jets_refined.json').read_bytes())
        oldm=_interval_from_exact(c,old['fixed_beta2_mass_upper_normalized'])
        difference=mass2-oldm
        checks=dict(passed=True,new_schedule_sha256=sha,all_14_nonnegative_true_masses_present=True,
            exact_axial_turnoff_primitive_contained=True,exact_axis_series_coefficients_contained=13,
            total_independent_exact_containments=contains,stage_mass_sums_consistent=True,
            explicit_Md_greater_than_one=True,explicit_logPstar_greater_than_Td=True,
            stored_y_d_minus_unrounded_y_d=stored_yd-yd,
            normalized_beta2_mass_difference_from_old=difference,
            new_core_generated=False,waiting_root_certified=False,
            unknown_paper_constants_verified=False,
            input_hashes={n:hashlib.sha256((here/n).read_bytes()).hexdigest() for n in
                (name,Path(__file__).name,'lei_ren_part1_paper_Md11_pressure_datum.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(checks),indent=2)+'\n',encoding='utf-8')
    print('New Md1.1 datum checks passed;14 stages;14 exact containments; schedule fingerprint and dependencies',flush=True)
    return checks

if __name__=='__main__':run()

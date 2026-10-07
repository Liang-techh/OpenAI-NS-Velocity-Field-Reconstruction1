"""Actual C0 signed Duhamel contributions from whole native O2-slope cells.

The existing source/phase interval kernels are integrated with the positive
original recovery weights. No sampled density defines an integral. These
are local contributions, before global inlet histories, C1 Z jets or repair.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_spatial_phase as spatial

density=spatial.density;phase=spatial.phase;prior=spatial.prior;native=spatial.native;packets=spatial.packets
HERE,PREFIX,sha=spatial.HERE,spatial.PREFIX,spatial.sha;ep=spatial.ep
NAME=PREFIX+'current_native_local_signed_integrals.json'
RECEIPT=PREFIX+'current_native_local_signed_integrals_check.json'
GATE='current_actual_native_C0_local_signed_Duhamel_contributions_executed'
RATES={key:Fraction(value) for key,value in density.RATES.items()}


def positive_kernel_mass(c,width,rate):
    """Enclose integral_0^width exp(-rate*t)dt without 1-exp cancellation."""
    w=c.mpf(width);r=Fraction(rate)
    if ep(w)[0]<0 or any(not mp.isfinite(v) for v in ep(w)) or r<0:
        raise ValueError('Finite nonnegative source width and recovery rate required')
    if not r:return w
    return w*packets.recovery.exp_average(c,-w*c.mpf(r.numerator)/r.denominator)


def same_source_union(values):
    """Hull a phase-cell union in one shared original formal log basis."""
    if not values:raise ValueError('Nonempty actual phase-cell union required')
    first=values[0];c=first.ctx
    for value in values:first.coerce(value)
    if len(values)==1:return first
    ref=max(values,key=lambda value:ep(value.scale.evaluate())[1]).scale
    threshold=2*c.dps*mp.log(10)
    if any(ep((value.scale-ref).evaluate())[1]>threshold for value in values if not value.zero):
        # Only an arithmetic coordinate: not a source-selected field scale.
        ref=prior.FormalScale(first.scale.bases,offset=c.mpf(max(ep(v.scale.evaluate())[1] for v in values)))
        first.ledger['directed_independent_log_rescalings']+=1
    coefficients=[c.mpf(0) if value.zero else value.coefficient*value.bounded_exp((value.scale-ref).evaluate()) for value in values]
    return prior.ScaledEnclosure(ref,c.mpf((min(ep(v)[0] for v in coefficients),max(ep(v)[1] for v in coefficients))),first.ledger)


class NativeLocalSignedIntegrals:
    def __init__(self,owner,binder):
        if type(owner) is not density.NativeCandidateDensities or type(binder) is not spatial.NativeSpatialPhase:
            raise ValueError('Actual existing native density and radius owners required')
        if owner.family!=binder.family or owner.owner.owner.owner.native.seed is not binder.seed:
            raise ValueError('One same-family original native seed required')
        receipt=json.loads((HERE/spatial.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[spatial.GATE] or receipt['source_family']!=binder.family:
            raise ValueError('Accepted actual source radius and spatial phase required')
        if {key:Fraction(value) for key,value in packets.recovery.RATES.items()}!=RATES:
            raise ValueError('Same original five recovery rates required')
        self.owner=owner;self.binder=binder;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (spatial.RECEIPT,
            PREFIX+'current_generic_shear_moment_recovery.py',PREFIX+'current_O3_independent_repair_operator.py',Path(__file__).name)})

    def contribution(self,*,Z,left,right,N,chart='O2_slope'):
        # Only this chart has the installed exact coordinate-length identity
        # dy=dcoordinate. Other charts need their own true Jacobian integration.
        if chart!='O2_slope':raise ValueError('Current local integral scope is the original O2_slope chart')
        lo=spatial.exact_coordinate(left);hi=spatial.exact_coordinate(right)
        if lo is None or hi is None or not 0<=lo<hi<=1:
            raise ValueError('Explicit exact ordered O2_slope endpoints in[0,1] required')
        c=self.ctx;cv=lambda q:c.mpf(q.numerator)/q.denominator
        coordinate=c.mpf((ep(cv(lo))[0],ep(cv(hi))[1]));width=cv(hi-lo)
        geometry=self.binder.query(chart,Z,coordinate,N)
        source,loop,V=self.owner.source(chart,Z,geometry['raw']['coordinate'])
        candidates=[];cells=[]
        for box in geometry['phase_boxes']:
            got=density.candidate_at_phase(loop,V,box,N)
            if got['values'] is None:raise ArithmeticError('Refine actual source/phase cell before integration')
            candidates.append(got);cells.append(spatial.spatial_candidate_record(got,geometry['record']))
        kernels={key:same_source_union([got['densities'][key] for got in candidates]) for key in RATES}
        masses={key:positive_kernel_mass(c,width,rate) for key,rate in RATES.items()}
        contributions={key:kernels[key]*masses[key] for key in RATES}
        record=dict(source_family=self.family,chart=chart,Z_box=geometry['raw']['Z'],
            exact_native_coordinate_endpoints=[spatial.fractional_record(lo),spatial.fractional_record(hi)],
            exact_log_radius_width_fraction=spatial.fractional_record(hi-lo),log_radius_width_enclosure=width,
            original_coordinate_identity='dy=dcoordinate on original O2_slope; dy/dcoordinate=1',
            candidate_N=N,actual_whole_cell_spatial_phase=geometry['record'],original_native_source=source['record'],
            actual_phase_union_candidate_cells=cells,whole_cell_signed_density_kernel_covers={k:v.record() for k,v in kernels.items()},
            original_positive_Duhamel_kernel_mass_covers=masses,
            five_actual_local_signed_Duhamel_contribution_enclosures={k:v.record() for k,v in contributions.items()},
            actual_integrals='I_j(Z)=integral_left^right exp(-lambda_j*(right-s))*f_j(s,Z) ds',
            recovery_rates={key:str(value) for key,value in RATES.items()},
            entire_radial_and_Z_cell_source_covers_used=True,actual_phase_union_not_samples=True,
            positive_kernel_signed_range_integration=True,no_extra_N_or_radius_or_velocity_unit_factor=True,
            original_P0_same_common_unit_unchanged=geometry['raw']['P0'],
            pressure_rate_zero_and_memory_not_reset=True,local_contributions_are_not_global_defect_histories=True,
            incoming_history_not_assumed_or_reset=True,C0_integral_function_enclosures_only=True,
            actual_C1_Z_derivatives_or_Rc_targets_installed=False,global_common_N_admitted=False,
            **dict.fromkeys(packets.OPEN,False))
        return dict(record=record,contributions=contributions,kernels=kernels,width=width,source=source,geometry=geometry)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        binder=spatial.NativeSpatialPhase(native.NativeGenericSourcePackets(bridge))
        density_owner=density.NativeCandidateDensities(phase.NativeConditionedPhase(spatial.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(binder.native))))
        owner=NativeLocalSignedIntegrals(density_owner,binder);N=1024;records={}
        for name,Z in (('Z_half',('.5','.5')),('Z_interval',('.49','.51'))):
            got=owner.contribution(Z=Z,left='.13369999',right='.13370001',N=N);records[name]=got['record']
            print('Actual whole-cell signed integrals:',name,{k:v.record()['sign'] for k,v in got['contributions'].items()},flush=True)
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=N,native_radial_cell_count=1,
        native_Z_query_count=len(records),signed_local_integral_enclosures_count=len(records)*5,
        actual_native_local_signed_integral_records=records,
        actual_C0_signed_local_integral_functions_installed=True,global_cumulative_histories_or_C1_Z_jets_or_repair_installed=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Actual signed C0 local Duhamel contributions over original O2_slope coordinate [.13369999,.13370001] at Z=.5 and whole Z interval [.49,.51], using true spatial phase at candidate N=1024. No global incoming histories/C1 Z/Rc targets/repair/common N/cone/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

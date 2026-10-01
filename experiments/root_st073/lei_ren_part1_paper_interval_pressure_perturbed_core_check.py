"""Propagate the accepted pressure-integral error to the finite collar endpoint.

This is conditional on the stored axis and finite model. Original parameter,
adapter-generation, infinite core and omitted-pressure/width errors are open.
"""
import json
import time
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_pressure_perturbed_core import regenerate_core,encode_receipt_metadata
from lei_ren_part1_paper_interval_collar_core_inlet import build_inlet_from_interval_core,encode_snapshot
from lei_ren_part1_paper_interval_collar_core_inlet_check import encode_result
from lei_ren_part1_paper_interval_axial_second_jet import IntervalAxialSecondJet
from lei_ren_part1_paper_interval_pressure_width_jet import IntervalPressureWidthJet
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_collar_second_width_moments import endpoint_integrated_switch,second_width_coefficients
from lei_ren_part1_paper_collar_first_width_interval_endpoint import first_width_coefficients,encode_coefficients
from lei_ren_part1_paper_collar_width_physical_endpoint import physical_endpoint


def error_bounds(jet):
    """Sum atom error bounds at pressure=1 separately for each width order.

Width powers remain distinct; this does not bound a full physical field.
"""
    slots={name:getattr(jet,name) for name in ('value','tangent','second')} if isinstance(jet,IntervalAxialSecondJet) else {'value':jet}
    result={}
    for slot,ring in slots.items():
        rows={}
        for w in range(3):
            size=ring.ctx.mpf(0)
            for (p,width),pair in ring.atoms.items():
                if width!=w:continue
                lo,hi=endpoints(pair.difference);size+=ring.ctx.mpf(max(abs(lo),abs(hi)))
            rows[str(w)]=endpoints(size)[1]
        result[slot]=rows
    return result


def run():
    started=time.monotonic()
    with mp.workdps(520):
        data=regenerate_core()
        print('pressure-perturbed degree18 core regenerated',flush=True)
        ctx=data['ctx'];core=data['core'];source=data['raw_source']
        inlet=build_inlet_from_interval_core(core,ctx,data['Lambda'],data['delta'])
        print('paired pressure errors propagated to inlet',flush=True)
        profile,_=accepted_profile();switch=endpoint_integrated_switch(ScheduleEndpointEnclosures(profile.schedule))
        K=ctx.mpf(list(endpoints(switch['integrated_switch_primitive'])))
        first=first_width_coefficients(inlet);second=second_width_coefficients(inlet,K)
        physical=physical_endpoint(inlet,first,second)
        audit={}
        for name in ('F','Uz','P'):
            count=0;outside=0;relative=mp.mpf(0)
            for n,row in enumerate(core[name]):
                for k,jet in enumerate(row[:4]):
                    for p in range(10):
                        stored=source['coefficients'][name][n][k].atoms.get(p,mp.mpf(0))
                        lo,hi=endpoints(jet.component(p,0).nominal)
                        distance=max(lo-stored,stored-hi,mp.mpf(0))
                        count+=1;outside+=distance!=0
                        if stored:relative=max(relative,distance/abs(stored))
            audit[name]=dict(comparisons=count,legacy_values_outside_nominal_intervals=outside,
                             maximum_relative_distance=relative,
                             legacy_generation_roundoff_not_enclosed=True)
        if physical['P0'] is not inlet['P0']:raise AssertionError('Pressure datum replaced')
        fields={name:error_bounds(inlet[name]) for name in ('F','Uz','P','I_theta','I_z','D','D_R','I_z_R')}
        endpoint={name:error_bounds(physical[name]) for name in ('F','Uz','P','Ur','Ur_Z')}
        moment_errors={name:error_bounds(jet) for name,jet in physical['moments'].items()}
        report=dict(error_provenance=encode_receipt_metadata(data['receipt_metadata']),
            conditional_pressure_integral_error_propagated=True,
            pressure_prefix_and_post_bounds_applied_separately=True,
            pressure_subset_bound_conservatism='each subset receives full combined bound; at pressure=1 seed bound can double',
            inlet=encode_snapshot(inlet),first_width_coefficients=encode_coefficients(first),
            second_width_coefficients=encode_coefficients(second),physical_endpoint=encode_result(physical),
            inlet_pressure_error_bounds=encode(fields),endpoint_pressure_error_bounds_by_width_order=encode(endpoint),
            five_endpoint_moment_pressure_error_bounds_by_width_order=encode(moment_errors),
            legacy_coefficient_nominal_audit=encode(audit),original_pressure_datum_preserved=True,
            pressure_parameter_order=9,width_order=2,radial_degree=18,Z='.3',precision=473,
            axis_or_adapter_generation_roundoff_enclosed=False,original_parameter_errors_enclosed=False,
            infinite_core_remainder_enclosed=False,omitted_pressure_orders_enclosed=False,
            omitted_width_orders_enclosed=False,full_collar_ODE_error_enclosed=False,
            global_axial_domain_certified=False,temporal_recursion=False,elapsed_seconds=time.monotonic()-started)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode_result(report),indent=2)+'\n')
        for name in ('Uz','P','D'):
            print(name,'inlet conditional pressure error',mp.nstr(fields[name]['value']['0'],16),flush=True)
        print('conditional endpoint receipt saved, seconds',round(time.monotonic()-started,2),flush=True)


if __name__=='__main__':run()

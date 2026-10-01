"""Actual accepted finite-core collar endpoint interval receipt at Z=.3.

No expensive source core or spatial collar chain is rebuilt. Stored source
coefficients are exact inputs here; their approximation errors remain open.
"""
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from lei_ren_part1_paper_interval_collar_core_inlet import build_inlet,encode_snapshot,DEFAULT_CACHE,_read_cache
from lei_ren_part1_paper_interval_axial_second_jet import IntervalAxialSecondJet
from lei_ren_part1_paper_interval_pressure_width_jet import IntervalPressureWidthJet
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_collar_second_width_moments import endpoint_integrated_switch,second_width_coefficients
from lei_ren_part1_paper_collar_first_width_interval_endpoint import first_width_coefficients,encode_coefficients
from lei_ren_part1_paper_collar_width_physical_endpoint import physical_endpoint
from lei_ren_part1_paper_pressure_width_second_axial_comparison import SecondAxialPressureWidthComparison
from lei_ren_part1_paper_component_pressure_core import evaluate_component_core_coefficients
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


def encode_ring(ring):
    return {str(key):dict(nominal=encode(pair.nominal),perturbation=encode(pair.difference))
            for key,pair in ring.atoms.items()}


def encode_result(value):
    if isinstance(value,IntervalAxialSecondJet):
        return {slot:encode_ring(getattr(value,slot)) for slot in ('value','tangent','second')}
    if isinstance(value,IntervalPressureWidthJet):return encode_ring(value)
    if isinstance(value,dict):return {key:encode_result(item) for key,item in value.items() if key!='ctx'}
    return encode(value)


def run():
    begin=time.monotonic()
    with mp.workdps(520):
        inlet=build_inlet(precision=473)
        print('actual directed finite-core inlet available',flush=True)
        profile,_=accepted_profile();e=ScheduleEndpointEnclosures(profile.schedule)
        switch=endpoint_integrated_switch(e)
        K=inlet['ctx'].mpf(list(endpoints(switch['integrated_switch_primitive'])))
        first=first_width_coefficients(inlet);second=second_width_coefficients(inlet,K)
        physical=physical_endpoint(inlet,first,second)
        if physical['P0'] is not inlet['P0']:raise AssertionError('Pressure datum replaced')
        if physical['Ur_ZZ'] is not None:raise AssertionError('Unknown third derivative fabricated')
        # Compare to the original low-precision inlet, preserving discrepancies
        # instead of expanding the directed intervals to force inclusion.
        raw,_,_=_read_cache(DEFAULT_CACHE)
        # Independent radial differentiation of scalar, pressure-power-zero
        # core data checks the derivative formulas used by dynamic g2/u2.
        scalar_core=dict(raw['coefficients'])
        for name in ('F','Uz','P'):
            scalar_core[name]=[[entry.atoms.get(0,mp.mpf(0)) for entry in row]
                               for row in raw['coefficients'][name]]
        Ra=mp.mpf(4)/raw['Lambda'];z=mp.mpf('.3');delta=raw['delta']
        def scalar_drivers(t):
            r=Ra*t;fields=evaluate_component_core_coefficients(scalar_core,r,z,delta)
            root=mp.sqrt(2*r)
            stress=evaluate_mp_stress(0,z,delta,Utheta=root*fields['F'],Uz=fields['Uz'],
                Utheta_y=root*(r*fields['F_R']+fields['F']/2),Utheta_Z=root*fields['F_Z'],
                Uz_y=r*fields['Uz_R'],Uz_Z=fields['Uz_Z'],moments=fields['moments'],
                moments_Z=fields['moments_Z'],P=fields['P'],P_Z=fields['P_Z'],
                precision=mp.mp.dps,radius_override=r)
            return dict(I_theta=stress['I_theta'],I_z=stress['I_z'],D=stress['I_theta']/fields['F'])
        radial_audit=[]
        for name in ('I_theta','I_z','D'):
            oracle=mp.diff(lambda t:scalar_drivers(t)[name],1)
            interval=(inlet['R']*inlet[name+'_R']).value.component(0,0).nominal
            lo,hi=endpoints(interval)
            if not lo<=oracle<=hi:raise AssertionError((name,'radial derivative oracle outside interval'))
            radial_audit.append(dict(field=name,derivative_convention='R*d/dR',
                                     independent_oracle=oracle,directed_interval=interval))
        with mp.workdps(260):
            axis=SimpleNamespace(precision=260,Lambda=raw['Lambda'],delta=raw['delta'])
            def coefficients(z):
                if mp.mpf(str(z))!=mp.mpf('.3'):raise ValueError('Only Z=.3 cached')
                return raw['coefficients']
            component=SimpleNamespace(axis=axis,precision=260,coefficients=coefficients)
            comparison=SecondAxialPressureWidthComparison(
                dict(precision=260,axis=axis,component_pressure_core=component),
                h_b=mp.exp(-100-100*mp.mpf('1e152')),pressure_order=9,width_order=2)
            legacy=comparison.evaluate(0,'.3')
        audit=[];outside=0
        for name in ('F','Uz','P','I_theta','I_z','D'):
            for slot in ('value','tangent','second'):
                for p in range(10):
                    stored=getattr(legacy[name],slot).atoms.get((p,0),mp.mpf(0))
                    interval=getattr(inlet[name],slot).component(p,0).nominal
                    lo,hi=endpoints(interval);contained=lo<=stored<=hi
                    outside+=not contained
                    audit.append(dict(field=name,Z_slot=slot,pressure_power=p,
                        legacy_value=stored,directed_interval=interval,
                        legacy_value_contained=contained,
                        distance_to_interval=max(lo-stored,stored-hi,mp.mpf(0))))
        # Scalar display values only; exact full atom intervals are below.
        display={name:mp.nstr(sum(endpoints(jet.value.component(0,0).nominal))/2,18)
                 for name,jet in second.items()}
        report=dict(precision=473,Z='.3',accepted_cache_sha256=inlet['cache_sha256'],
            inlet=encode_snapshot(inlet),weighted_switch=encode(switch),
            first_width_coefficients=encode_coefficients(first),
            second_width_coefficients=encode_coefficients(second),
            physical_endpoint=encode_result(physical),
            second_width_pressure_power_zero_display=display,
            legacy_comparisons=encode(audit),legacy_comparison_count=len(audit),
            legacy_values_outside_high_precision_intervals=outside,
            independent_radial_driver_derivative_checks=encode(radial_audit),
            actual_finite_core_endpoint_evaluated=True,
            all_eight_states_at_width_orders_one_and_two_evaluated=True,
            all_five_endpoint_moments_evaluated=True,original_pressure_datum_preserved=True,
            source_coefficient_errors_enclosed=False,core_series_remainder_enclosed=False,
            omitted_width_orders_enclosed=False,full_ODE_error_enclosed=False,
            uniform_Z_certified=False,temporal_recursion=False,
            expensive_source_or_collar_chain_rebuilt=False,elapsed_seconds=time.monotonic()-begin)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print('eight second-width coefficients',display,flush=True)
        print('legacy comparisons',len(audit),'outside',outside,'seconds',round(time.monotonic()-begin,2),flush=True)


if __name__=='__main__':run()

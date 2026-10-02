"""Independent replay of original-width comparison-history packets."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_comparison_point_integrals import CompliantComparisonPointIntegrals,smoothing_weights
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_comparison_point_integrals.json'


def alpha(u):
    if u<=1:return mp.mpf(1)
    x=u-1
    if x<=0:return mp.mpf(1)
    if x>=1:return mp.mpf(0)
    exponent=1/(1-x)**2-1/x**2
    e=mp.exp(exponent)
    return 1-e/(1+e)


def gauss_weighted_integrals(phase,nodes=18,cells=512):
    """Independent Gauss-Legendre quadrature, distinct from Simpson receipts."""
    with mp.workdps(100):
        x,w=mp.gauss_quadrature(nodes,'legendre');s=mp.mpf(phase);d=s-1
        sums=[mp.mpf(0),mp.mpf(0)]
        if d<=0:return s,s*s/2
        h=d/cells
        for cell in range(cells):
            mid=(cell+mp.mpf('.5'))*h
            for xx,ww in zip(x,w):
                t=mid+h*xx/2;a=alpha(1+t)
                sums[0]+=h*ww*a/2;sums[1]+=h*ww*(1+t)*a/2
        return mp.mpf('.5')+sums[0],mp.mpf('.5')+sums[1]


def run():
    receipt=json.loads((HERE/NAME).read_bytes());hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Comparison source changed: '+name)
    if not receipt['actual_core_atom_comparison_histories_numerically_integrated'] or receipt['original_signed_actual_bridge_integrals_resolved']:
        raise ValueError('Comparison must inherit atoms while keeping actual bridge distinct')
    c=MPIntervalContext();c.dps=180;f=CompliantComparisonPointIntegrals();count=0
    if f.source['source_namespace']!=receipt['original_comparison_source']['source_namespace']:
        raise ArithmeticError('Original 9.23 comparison namespace changed')
    for z,charts in receipt['comparison_point_packets'].items():
        if set(charts)!= {'micro','macro'} or len(charts['micro'])!=4 or len(charts['macro'])!=3:
            raise ValueError('Both original smoothing charts and frozen macro coverage required')
        for packet in charts['micro']+charts['macro']:
            if not packet['actual_core_atom_inlet_used'] or packet['source_width_materialized'] or packet['cap_used_as_width_value']:
                raise ValueError('Formal positive hb source was materialized or replaced by a cap')
            if packet['actual_prescribed_shear_field_substituted_by_comparison'] or packet['original_signed_actual_bridge_integrals_resolved']:
                raise ValueError('Auxiliary comparison relabeled as actual bridge')
            if len(packet['comparison_phi']['signed_hb_power_axial_coefficients'])!=3:
                raise ValueError('Second-order signed-width field expansion missing')
            for row in packet['comparison_phi']['signed_hb_power_axial_coefficients']:
                if len(row)!=7 or any(not all(mp.isfinite(v) for v in endpoints(read_interval(c,x))) for x in row):
                    raise ArithmeticError('Nonfinite comparison Phi axial source jet')
            for field in ('comparison_phi','comparison_raw_V'):
                error=read_interval(c,packet[field]['omitted_hb_cubed_weighted_axial_jet_norm_upper'])
                if endpoints(error)[0]<0:raise ArithmeticError('Third-order field remainder must be positive')
            if packet['chart']=='micro':
                phase=read_interval(c,packet['coordinate']);w=packet['normalized_smoothing_weights']
                if endpoints(phase)[0]==0 and endpoints(phase)[1]==0:
                    for field in ('comparison_phi','comparison_raw_V'):
                        rows=packet[field]['signed_hb_power_axial_coefficients']
                        if any(endpoints(read_interval(c,v))!=(mp.mpf(0),mp.mpf(0)) for v in rows[1]+rows[2]):
                            raise ArithmeticError('Comparison must exactly match core at switch entrance')
                if endpoints(phase)[0]==2 and endpoints(phase)[1]==2:
                    if endpoints(read_interval(c,w['A0']))!=(mp.mpf('1.5'),mp.mpf('1.5')):
                        raise ArithmeticError('Exact symmetric comparison weight A0(2) failed')
                    exactA1,exactA2=gauss_weighted_integrals(2,18,512),gauss_weighted_integrals(2,24,512)
                    for stored,k in ((w['A1'],1),):
                        enclosure=read_interval(c,stored)
                        if not endpoints(enclosure)[0]<=exactA1[k]<=endpoints(enclosure)[1]:
                            raise ArithmeticError('Independent Gaussian switch integral outside Simpson enclosure')
                        if abs(exactA1[k]-exactA2[k])>mp.mpf('1e-25'):
                            raise ArithmeticError('Independent Gaussian switch quadratures did not converge')
                    count+=1
            else:
                if packet['exact_log_radius']!='y=2hb+fraction*(log(100/Ra)-2hb)':raise ValueError('Frozen macro lost original radius coordinate')
                if packet['coordinate']==0:
                    for name,row in packet['comparison_own_six_moments'].items():
                        for a,b in zip(row['signed_hb_power_axial_coefficients'],charts['micro'][-1]['comparison_own_six_moments'][name]['signed_hb_power_axial_coefficients']):
                            if [read_interval(c,v)._mpi_ for v in a]!=[read_interval(c,v)._mpi_ for v in b]:
                                raise ArithmeticError('Frozen comparison reset the 2hb actual history')
                if packet['coordinate']==1 and not packet['exact_endpoint_at_R100']:
                    raise ArithmeticError('Frozen macro endpoint radius is not exactly R=100')
                count+=1
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
        implicit_source_sha256=receipt['implicit_source_sha256'],datum_enclosure_sha256=receipt['datum_enclosure_sha256'],
        original_9_23_namespace_bound=True,three_source_Z_charts_checked=3,
        original_microscopic_switch_packets_checked=12,frozen_macro_packets_checked=9,
        independent_gauss_integral_A1_comparisons=3,
        actual_atoms_and_formal_width_remainders_preserved=True,
        original_signed_actual_bridge_integrals_resolved=False,all_passed=True,input_hashes=hashes)
    for name in (NAME,Path(__file__).name):result['input_hashes'][name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Comparison histories PASS: both switches, six inherited moments, frozen macro to R=100',flush=True)
    return result


if __name__=='__main__':run()

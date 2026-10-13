"""Check genuine row-one composition and full pressure-moment C1 coverage."""
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_n1_complete_pressure_C1 as current


def signed_source_fixture():
    c = MPIntervalContext()
    c.dps = 90
    results = []
    for ends in (['1','3'],['-3','-1'],['-2','3'],['0','0']):
        coefficient = c.mpf(ends)
        cover,error = current.signed_small_cover(c,coefficient,c.mpf(-1004))
        actual = coefficient*c.exp(-1005)
        if current.ep(cover)[0] > current.ep(actual)[0] or current.ep(cover)[1] < current.ep(actual)[1]:
            raise AssertionError('Signed nonzero physical derivative source was not covered')
        if ends == ['0','0']:
            if current.ep(error) != (mp.mpf(0),mp.mpf(0)):
                raise AssertionError('Actual exact-zero derivative changed')
        elif current.ep(error)[0] <= 0:
            raise AssertionError('A bounded derivative was silently set to zero')
        results.append(ends)
    return dict(positive_negative_mixed_and_exact_zero=results,
                independent_nonzero_physical_scale=True)


@current.source_precision
def run(owner=None):
    began = time.monotonic()
    owner = owner or current.CurrentOriginalN1CompletePressureC1(False)
    raw = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    if raw['source_family'] != owner.family or any(raw.get(k) for k in current.OPEN):
        raise AssertionError('Complete pressure source family or scope changed')
    for name,digest in raw['input_hashes'].items():
        if current.sha(name) != digest:
            raise AssertionError('Changed complete pressure source '+name)
    data = current.record(current.complete_source_data(owner.c,**owner.parents))
    if data != raw['source_data']:
        raise AssertionError('Actual cutoff/collar first-derivative archive differs')
    if len(data['inner_cut_cells']) != 4:
        raise AssertionError('Complete inner-cut Z cover missing')
    for cell in data['inner_cut_cells']:
        for key,field in cell['fields'].items():
            for term in field['terms']:
                actual = current.interval(owner.c,term['actual_ordinary_Z1_coefficient'])
                original = current.interval(owner.c,term['original_source_term']['ordinary_Z_derivatives'][1])
                if actual._mpi_ != original._mpi_:
                    raise AssertionError('The actual first derivative source was replaced')
                if not term['no_cap_or_affine_error_differentiation']:
                    raise AssertionError('A value error was differentiated')
    grid = owner.C1_grid()
    maxima = {}
    if len(grid) != 4096:
        raise AssertionError('Complete original C1 grid missing')
    for key in ('m13','m15'):
        maxima[key] = max(current.ep(v['fields'][key]['total_error'])[1] for v in grid)
        if maxima[key] >= current.ep(owner.c.mpf('.001'))[0]:
            raise AssertionError('Complete pressure moment derivative accuracy exceeds 1e-3')
        for row in grid:
            if not row['fields'][key]['absolute_accuracy_certified']:
                raise AssertionError('A full derivative cell is unaccepted')
            if current.ep(row['fields'][key]['inner_and_collar_nonzero_error'])[0] <= 0:
                raise AssertionError('The actual nonzero lower source disappeared')
    for Z in ('-1','-.5','0','.5','1',['-.511','-.489']):
        view = owner.derivative_enclosure(Z)
        if view['radial_scope'] != 'original inner cutoff plus Rin to infinity':
            raise AssertionError('Full moment domain missing at a source/grid seam')
    before = mp.mp.dps
    with mp.workdps(15):
        owner.derivative_enclosure('.371')
        if mp.mp.dps != 15:
            raise AssertionError('Public full C1 query leaked source precision')
    if before != mp.mp.dps:
        raise AssertionError('Public full C1 query changed caller precision')
    fixture = signed_source_fixture()
    pins = dict(owner.hashes)
    for name in (current.NAME,Path(current.__file__).name,Path(__file__).name):
        pins[name] = current.sha(name)
    receipt = dict(all_passed=True,source_family=owner.family,
        checks=dict(genuine_inner_Z1_terms=True,unchanged_collar_genuine_full_Z1_range=True,
            fixed_positive_original_units=True,nonzero_source_errors=True,
            full_four_Z_cut_cover=True,full_derivative_cells=4096,moment_rows=8192,
            grid_and_source_seam_queries=True,default_precision_restored=True,
            independent_signed_source_fixture=fixture,
            no_actual_function_or_repair_coefficient_selected=True),
        maximum_normalized_derivative_error=maxima,input_hashes=pins,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.record(receipt),separators=(',',':'))+'\n').encode())
    print('PASS_COMPLETE_PRESSURE_MOMENT_C1',len(grid),flush=True)
    return receipt


if __name__=='__main__':
    run()

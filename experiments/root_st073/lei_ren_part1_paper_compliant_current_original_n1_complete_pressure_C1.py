"""Complete original m13/m15 first-Z derivative enclosures before n1 repair.

Compose the genuine signed inner cutoff rows, unchanged Rin-to-r_minus
collar derivative ranges, and checked r_minus-to-infinity C1 source. Fixed
positive original moment units and nonzero source allowances are retained.
Numerical enclosures do not select the actual signed defect functions.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_C1 as upper
import lei_ren_part1_paper_compliant_current_original_n1_external_pressure_values as external
import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values as defects

HERE, PREFIX, sha, ep = upper.HERE, upper.PREFIX, upper.sha, upper.ep
record, interval, source_precision = upper.record, upper.interval, upper.source_precision
NAME = PREFIX+'current_original_n1_complete_pressure_C1.json.gz'
RECEIPT = PREFIX+'current_original_n1_complete_pressure_C1_check.json'
GATE = 'current_original_n1_complete_pressure_moment_C1_accuracy_installed'
OPEN = tuple(dict.fromkeys((*upper.OPEN,*defects.OPEN)))


def signed_small_cover(c, coefficient, log_upper):
    """Use the true coefficient sign and a directed nonzero source allowance."""
    lo, hi = ep(coefficient)
    if lo == hi == 0:
        return c.mpf(0), c.mpf(0)
    if ep(log_upper)[1] >= -1000:
        raise ValueError('Actual inner/collar derivative exceeds the original source budget')
    e = c.exp(-1000)
    cover = c.mpf([0,ep(e)[1]]) if lo >= 0 else (
        c.mpf([-ep(e)[1],0]) if hi <= 0 else c.mpf([-ep(e)[1],ep(e)[1]]))
    return cover,e


def complete_source_data(c, lower, external_data, defect_data):
    scales = external_data['source_scales']
    if not lower['Z2_is_half_the_ordinary_second_derivative'] or not lower['original_physical_measures_and_m15_sign_retained']:
        raise ValueError('Actual unchanged collar derivative recipe required')
    collar = {}
    for key in ('m13','m15'):
        normalization = lower['normalizations'][key]
        if normalization['source_basis_Z_independent_indices'] != [0,2,3]:
            raise ValueError('Fixed original collar units must be Z independent')
        derivative = interval(c,normalization['full_domain_first_derivative_range'])
        norm = external.absolute_upper(c,derivative)
        ratio = interval(c,scales['fields'][key]['exact_positive_ratio_log'])
        logarithm = ratio+c.ln(norm) if ep(norm)[1] else None
        cover,error = (signed_small_cover(c,derivative,logarithm) if logarithm is not None else (c.mpf(0),c.mpf(0)))
        collar[key] = dict(actual_full_Z_derivative=derivative,correlated_positive_ratio_log=ratio,
            normalized_log_absolute_upper=logarithm,normalized_signed_derivative=cover,absolute_upper=error,
            actual_first_derivative_not_C0_model_error=True,original_m15_sign_retained=True)
    cells = []
    for source_cell in defect_data['source_cells']:
        ends = tuple(source_cell['ends'])
        if ends not in upper.CELLS:
            raise ValueError('Original four whole-Z cutoff cells required')
        fields = {}
        for key,power in (('m13',1),('m15',2)):
            old = source_cell['inner_pressure'][key]
            actual = old['original_signed_Z1_sectors']
            if actual['ordinary_Z_derivative_order'] < 1:
                raise ValueError('Genuine signed inner cutoff first derivative required')
            unit = defects.pressure_unit_log(c,scales,power)
            if unit._mpi_ != interval(c,old['unit_log'])._mpi_:
                raise ValueError('Inner cut and exterior pressure moment units differ')
            terms = []
            for source_term in actual['terms']:
                rows = source_term['ordinary_Z_derivatives']
                if len(rows) < 2:
                    raise ValueError('Actual cutoff Z1 source row missing')
                coefficient = interval(c,rows[1])
                norm = external.absolute_upper(c,coefficient)
                logarithm = interval(c,source_term['combined_source_log'])-unit+c.ln(norm) if ep(norm)[1] else None
                cover,error = (signed_small_cover(c,coefficient,logarithm) if logarithm is not None else (c.mpf(0),c.mpf(0)))
                terms.append(dict(original_source_term=source_term,actual_ordinary_Z1_coefficient=coefficient,
                    normalized_log_absolute_upper=logarithm,normalized_signed_derivative=cover,absolute_upper=error,
                    exact_source_derivative_zero=ep(norm)[1]==0,no_cap_or_affine_error_differentiation=True))
            value = sum((v['normalized_signed_derivative'] for v in terms),c.mpf(0))
            bound = sum((v['absolute_upper'] for v in terms),c.mpf(0))
            fields[key] = dict(unit_log=unit,terms=terms,normalized_signed_derivative=value,absolute_upper=bound,
                actual_cutoff_source_Z1_retained=True)
        cells.append(dict(ends=ends,fields=fields))
    if {tuple(v['ends']) for v in cells} != set(upper.CELLS):
        raise ValueError('Complete inner-cut C1 axial cover required')
    errors = {key:collar[key]['absolute_upper']+max((v['fields'][key]['absolute_upper'] for v in cells),key=lambda v:ep(v)[1])
              for key in collar}
    return dict(collar=collar,inner_cut_cells=cells,additional_nonzero_errors=errors,
        radial_scope='original inner cutoff plus Rin to infinity',
        source_family=defect_data['source_family'],same_original_pressure_units=True,
        signed_inner_and_collar_rows_not_selected_as_field_values=True)


class CurrentOriginalN1CompletePressureC1:
    @source_precision
    def __init__(self,require_checked=True):
        self.upper = upper.CurrentOriginalN1PreRpPressureC1(True)
        self.c,self.family = self.upper.c,copy.deepcopy(self.upper.family)
        self.hashes = dict(self.upper.hashes)
        self.hashes[upper.RECEIPT] = sha(upper.RECEIPT)
        self.parents = dict(lower=upper.load_checked(external.lower,self.hashes,self.family),
            external_data=upper.load_checked(external,self.hashes,self.family),
            defect_data=upper.load_checked(defects,self.hashes,self.family))
        self.data = record(complete_source_data(self.c,**self.parents))
        self.extra = {key:interval(self.c,value) for key,value in self.data['additional_nonzero_errors'].items()}
        self.acceptance_loaded = False
        if require_checked:
            receipt = json.loads((HERE/RECEIPT).read_bytes())
            if not (receipt.get('all_passed') and receipt.get(GATE) and receipt['source_family']==self.family
                    and not any(receipt.get(k) for k in OPEN)):
                raise ValueError('Checked complete pressure-moment C1 receipt required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:
                    raise ValueError('Changed complete pressure C1 source '+name)
            raw = json.loads(gzip.decompress((HERE/NAME).read_bytes()))
            if raw['source_data'] != self.data:
                raise ValueError('Complete pressure C1 source rows differ')
            self.hashes.update(receipt['input_hashes'])
            self.acceptance_loaded = True
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        self._snapshot = json.dumps(record((self.data,self.family,self.hashes)),sort_keys=True)

    @source_precision
    def assert_graph(self):
        self.upper.assert_graph()
        if self._snapshot != json.dumps(record((self.data,self.family,self.hashes)),sort_keys=True):
            raise ValueError('Complete original pressure C1 source changed')

    @source_precision
    def C1_grid(self):
        self.assert_graph()
        grid = self.upper.upper_C1_grid()
        for model in grid:
            for key,field in model['fields'].items():
                field['inner_and_collar_nonzero_error'] = self.extra[key]
                field['total_error'] = interval(self.c,field['total_error'])+self.extra[key]
                field['absolute_accuracy_certified'] = ep(field['total_error'])[1] < ep(self.c.mpf('.001'))[0]
            model['radial_scope'] = 'original inner cutoff plus Rin to infinity'
        return grid

    @source_precision
    def derivative_enclosure(self,Z):
        self.assert_graph()
        previous = self.upper.derivative_enclosure(Z)
        lo,hi = ep(previous['Z'])
        values = {key:[] for key in self.extra}
        for cell in self.data['inner_cut_cells']:
            a,b = cell['ends']
            if ep(self.c.mpf(b))[1] < lo or ep(self.c.mpf(a))[0] > hi:
                continue
            for key in values:
                lower = interval(self.c,self.data['collar'][key]['normalized_signed_derivative'])
                inner = interval(self.c,cell['fields'][key]['normalized_signed_derivative'])
                values[key].append(previous['normalized_pressure_Z'][key]+lower+inner)
        if any(not rows for rows in values.values()):
            raise ValueError('Original complete cutoff derivative domain not covered')
        return dict(Z=previous['Z'],source_family=self.family,
            normalized_pressure_moment_Z={key:self.c.mpf([min(ep(v)[0] for v in rows),max(ep(v)[1] for v in rows)])
                                          for key,rows in values.items()},
            radial_scope='original inner cutoff plus Rin to infinity',
            original_signed_inner_and_collar_sources_retained=True,
            numerical_enclosures_not_signed_defect_function_selection=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))


@source_precision
def run():
    began = time.monotonic()
    owner = CurrentOriginalN1CompletePressureC1(False)
    grid = owner.C1_grid()
    maximum = {key:max(ep(v['fields'][key]['total_error'])[1] for v in grid) for key in owner.extra}
    raw = dict(source_family=owner.family,source_data=owner.data,
        C1_composition=dict(parent=upper.NAME,parent_sha256=sha(upper.NAME),grid_count=len(grid),
            added_inner_collar_error=owner.extra,maximum_normalized_derivative_error=maximum,
            radial_scope='original inner cutoff plus Rin to infinity',original_positive_units_fixed=True),
        samples={z:record(owner.derivative_enclosure(z)) for z in ('-1','-.639','0','.211','.371','.639','1')},
        source_bindings=dict(inner='actual cutoff ordinary_Z_derivatives[1]',
            lower='unchanged collar full_domain_first_derivative_range',
            units='five_defect_values.pressure_unit_log and external.source_scales',
            upper=upper.NAME),input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(record(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('COMPLETE_PRESSURE_MOMENT_C1_COMPOSED',len(grid),flush=True)
    return owner,raw


if __name__=='__main__':
    run()

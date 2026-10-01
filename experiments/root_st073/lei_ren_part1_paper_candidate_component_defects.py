"""Candidate-specific center moment defects from persisted actual R110 data.

Formal component algebra and finite quadrature only. These are not global
defect functions, integral certificates, or solved five-bump controls.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_centered_component_defects import CenteredComponentDefects

HERE = Path(__file__).resolve().parent
INPUT = 'lei_ren_part1_paper_candidate_transition_analytic.json'


class CandidateEndpointSource:
    def __init__(self):
        raw = (HERE/INPUT).read_bytes()
        self.input_hash = hashlib.sha256(raw).hexdigest()
        self.receipt = json.loads(raw)
        for name,digest in self.receipt['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Actual transition input changed: '+name)
        self.precision = 260
        with mp.workdps(500):
            # Preserve the canonical decimal across downstream MP contexts;
            # comparing differently rounded MP values can falsely reject it.
            self.delta = '1e-200'
            self.center = mp.mpf(self.receipt['center'])
            def unpack(record):
                if 'exact_mpf_tuple' in record:
                    return mp.make_mpf(tuple(record['exact_mpf_tuple']))
                return {k:unpack(v) for k,v in record.items()}
            self.endpoint = unpack(self.receipt['R110_refined_endpoint_exact_mpf'])
            self.endpoint.update(F_Z=self.endpoint['FZ'],Uz_Z=self.endpoint['UZ'],
                                 P_Z=self.endpoint['PZ'],moments_Z=self.endpoint['momentsZ'])
        self.shared_parameters = dict(logCstar='5e151',logPstar='14',delta='1e-200',Lambda='1e120')

    def evaluate_R(self,R,Z):
        with mp.workdps(500):
            if abs(mp.mpf(R)-110) > mp.mpf('1e-255'):
                raise ValueError('Persisted source only supplies the R110 endpoint')
            if abs(mp.mpf(Z)-self.center) > mp.mpf('1e-255'):
                raise ValueError('Persisted source only supplies its actual axial center')
            return self.endpoint


def packed_atoms(jet):
    return [dict(pressure_order=p,width_order=w,exact_mpf_tuple=list(v._mpf_))
            for (p,w),v in sorted(jet.atoms.items())]


def run():
    source = CandidateEndpointSource()
    centered = CenteredComponentDefects(source,order=16,restore_order=32)
    packet = centered.evaluate(source.center)
    with mp.workdps(centered.work_precision):
        rows = {str(row):packed_atoms(value) for row,value in packet['d'].items()}
        tangents = {str(row):packed_atoms(value) for row,value in packet['d_Z'].items()}
        nominal = {str(row):mp.nstr(value.evaluate(1,1),70) for row,value in packet['d'].items()}
        parts = {name:{str(row):packed_atoms(value) for row,value in data.items()}
                 for name,data in packet['parts'].items()}
        dependencies = (INPUT,'lei_ren_part1_paper_centered_component_defects.py',
                        'lei_ren_part1_paper_axial_dual.py','lei_ren_part1_paper_pressure_width_jet.py')
        report = dict(center='.3',Lambda='1e120',actual_R110_source_used=True,
                      source_transition_sha256=source.input_hash,
                      candidate_core_state_sha256=source.receipt['state_sha256'],
                      d_atoms=rows,d_Z_atoms=tangents,component_parts_atoms=parts,
                      nominal_row_values=nominal,metadata=packet['metadata'],
                      source_pressure_and_amplitude_fitted=False,
                      source_endpoint_is_nominal_fixed_parameter_projection=True,
                      source_pressure_parameter_remainder_enclosed=False,
                      functional_d_Z_domain_certified=False,
                      quadrature_error_enclosed=False,five_bump_controls_solved=False,
                      terminal_functional_five_moment_closure=False,
                      input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                    for name in dependencies},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Candidate actual-endpoint five defect rows and analytic first axial tangents generated; no closure claim',flush=True)
        return report


if __name__ == '__main__':
    run()

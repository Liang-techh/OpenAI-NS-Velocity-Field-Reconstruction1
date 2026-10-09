"""Actual whole-Z leading Rm..Rh patch functions on closed radial cells."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_centered_inverse as upstream
import lei_ren_part1_paper_compliant_current_original_Rm_patch_mixed4_cells as mixed

HERE,PREFIX,sha,read,bind,ep=upstream.HERE,upstream.PREFIX,upstream.sha,upstream.read,upstream.bind,upstream.ep
NAME=PREFIX+'current_original_whole_Z_Rm_patch_mixed4.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_Rm_patch_mixed4_check.json'
GATE='current_original_whole_Z_leading_Rm_Rh_closed_cell_mixed4_functions_installed'
OPEN=upstream.OPEN
PARTITION=((1,1),(49,40),(5,4),(51,40),(59,40),(3,2),(61,40),(69,40),(7,4),(71,40),(2,1),'Rh')


def serialized(value):return upstream.upstream.reference.serialized(upstream.leading.pack(value))


class WholeZRmPatchMixed4:
    def __init__(self,dps=500):
        self.source=upstream.WholeZRmCenteredInverse(dps);self.c=self.source.c
        self.identity=self.source.identity;self.hashes=dict(self.source.hashes);self.N=self.source.N
        checked=json.loads((HERE/upstream.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(upstream.GATE) or checked['source_family']!=self.identity \
                or not checked.get('same_common_unit_box_proves_unique_whole_Z_smooth_leading_family'):
            raise ValueError('Checked current whole-Z smooth leading controls required')
        for name,digest in checked['input_hashes'].items():bind(self.hashes,name,digest)
        flat=json.loads((HERE/mixed.FLAT).read_bytes())
        if not all(flat.get(gate) for gate in ('all_passed','original_radial_shape_derivatives_C4_available',
            'quantitative_flat_support_majorants_available','support_crossing_derivatives_available')):
            raise ValueError('Checked original beta flat derivative and crossing envelopes required')
        if flat['actual_five_defect_family_sha256']!=self.identity['actual_five_defect_family_sha256'] \
                or flat['implicit_source_sha256']!=self.identity['implicit_source_sha256']:
            raise ValueError('Same original beta source and whole-Z control family required')
        for name,digest in flat['input_hashes'].items():bind(self.hashes,name,digest)
        self.bindings=mixed.source_bindings();self.owners={};self.cache={}
        for name in (upstream.RECEIPT,mixed.FLAT,Path(mixed.__file__).name,Path(mixed.original.__file__).name,
            'lei_ren_part1_paper_compliant_flat_pulse_derivatives.py',
            'lei_ren_part1_paper_compliant_current_patch_stress_operator.py',Path(__file__).name):bind(self.hashes,name,sha(name))

    def owner(self,ends):
        key=tuple(ends)
        if key not in self.owners:
            source=self.source.owner(ends)
            self.owners[key]=mixed._MixedPatchOwner(source.patch)
        return self.owners[key]

    def query(self,ends,left,right=None):
        op=self.owner(ends);key=(tuple(ends),left,right)
        if key in self.cache:return self.cache[key]
        value=dict(op.evaluate(left) if right is None else op.cell(left,right))
        # The old pure backend's combined local-frame flag is replaced by
        # explicit wrapper admissions; finite-N density remains a new layer.
        value.pop('source_frames_conditional_on_same_accepted_Rm_inlet')
        value.pop('whole_axis_function_provider_or_finite_N_Rc_patch_installed')
        value.update(source_identity=self.identity,exact_Z_cell=list(ends),
            whole_Z_original_leading_patch_function_provider_installed=True,
            whole_Z_leading_control_common_box_admitted=True,
            actual_patch_finite_N_density_oracle_installed=False,
            genuine_finite_N_correction_incoming_not_used_as_leading_background=True)
        self.cache[key]=value;return value


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZRmPatchMixed4();cells=[]
        for ends in upstream.source.CELLS:
            op=owner.owner(ends);radial=[]
            for left,right in zip(PARTITION,PARTITION[1:]):
                radial.append(serialized(owner.query(ends,left,right)))
                print('Whole-Z actual mixed4 patch cell: '+str(ends)+' '+str(left)+' '+str(right),flush=True)
            cells.append(dict(exact_Z_cell=list(ends),source_identity=owner.identity,
                source_radial_cells=radial,actual_Rm_function=serialized(owner.query(ends,(1,1))),
                actual_Rh_function=serialized(owner.query(ends,'Rh')),
                exact_common_P0_axial5=serialized(op.op.P0),exact_Rm_radius=serialized(op.op.Rm_factor),
                exact_Rh_radius=serialized(op.op.Rm_factor*owner.c.exp(1))))
        report=dict(**{GATE:True},source_family=owner.identity,exact_Z_domain=['-1','1'],
            exact_Z_partition=[list(ends) for ends in upstream.source.CELLS],
            exact_x_domain=['1','exp(1)'],exact_x_partition=serialized(PARTITION),source_cells=cells,
            original_source_bindings=owner.bindings,whole_Z_leading_Rm_implicit_controls_solved=True,
            full_radial_mixed4_leading_patch_installed=True,
            current_same_N_finite_correction_after_Rm_installed=False,
            candidate_N_not_used_for_leading_patch=owner.N,**dict.fromkeys(OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual whole-Z original leading patch velocities/pressure/five primitives '
                  'and x/y/R mixed derivatives through total order4 on all closed radial cells Rm..Rh. '
                  'Genuine same-N post-Rm correction transport, final Rc repair, cone/heat and n-recursion stay open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(upstream.source.bridge._encode(report),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z actual leading Rm..Rh closed-cell mixed4 atlas generated',flush=True);return report


if __name__=='__main__':run()

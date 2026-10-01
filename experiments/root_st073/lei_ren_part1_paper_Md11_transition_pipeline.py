"""Fresh pressure/core -> comparison -> R110 exit -> five-moment inverse.

Generic equations are reused with explicit new-source arguments. No frozen
old field/defect/control receipt is substituted. Strong cone coverage and
whole-axis matching remain separate gates.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_Md11_comparison_jets import Md11ComparisonJets,STEM,TAIL
from lei_ren_part1_paper_interval_comparison_cells import build_cells
from lei_ren_part1_paper_interval_exit_bridge_enclosure import integrate_exit
from lei_ren_part1_paper_interval_exit_continuation_enclosure import continue_exit,_pack
from lei_ren_part1_paper_interval_exit_switch_enclosure import integrate_switch
from lei_ren_part1_paper_interval_functional_defects import assemble
from lei_ren_part1_paper_interval_five_bump_inverse import certify,weights
from lei_ren_part1_paper_interval_repaired_reference_field import validate_receipt
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_Md11_'


def run():
    calc=Md11ComparisonJets(4);c=calc.ctx
    names=(Path(__file__).name,'lei_ren_part1_paper_Md11_comparison_jets.py',STEM+'_state.json',TAIL,
        'lei_ren_part1_paper_interval_comparison_cells.py','lei_ren_part1_paper_interval_exit_bridge_enclosure.py',
        'lei_ren_part1_paper_interval_exit_continuation_enclosure.py','lei_ren_part1_paper_interval_exit_switch_enclosure.py',
        'lei_ren_part1_paper_interval_functional_defects.py','lei_ren_part1_paper_interval_five_bump_inverse.py')
    hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    produced={}
    def save(label,data):
        # Legacy generic producers emit their old default source label even
        # when their equations receive a supplied calculator. Replace that
        # label with the actual source identity before writing any receipt.
        data.pop('accepted_schedule_sha256',None)
        data.update(new_schedule_sha256=calc.new_schedule_sha256,source_pressure_Md='1.1',
            state_sha256=calc.state_hash,center_family=calc.center_family,input_hashes=dict(hashes),
            all_inputs_from_new_pressure_core=True,old_field_receipts_reused=False)
        filename=PREFIX+label+'.json'
        (HERE/filename).write_text(json.dumps(encode(_pack(data)),indent=2)+'\n',encoding='utf-8')
        produced[label]=filename;hashes[filename]=hashlib.sha256((HERE/filename).read_bytes()).hexdigest()
        print('Fresh Md1.1 stage saved:',label,flush=True)
    with mp.workdps(calc.precision+60):
        comparison=build_cells(calc,cells_per_half=16);save('comparison_cells',comparison)
        bridge=integrate_exit(calc,comparison,cells_per_half=16);save('exit_bridge',bridge)
        continuation=continue_exit(calc,bridge,comparison['endpoint'],comparison['initial_phi'],target_R='100')
        save('exit_continuation',continuation)
        switch=integrate_switch(calc,continuation,cells_per_half=16);save('exit_switch',switch)
        defects=assemble(calc,switch,bridge,continuation,restore_cells=512);save('functional_defects',defects)
        wname='lei_ren_part1_paper_bump_integral_enclosures_check.json'
        wraw=json.loads((HERE/wname).read_bytes());validate_receipt(wraw,'dimensionless bump weights')
        hashes[wname]=hashlib.sha256((HERE/wname).read_bytes()).hexdigest()
        z=calc.z.truncate(1);Am=(1+z*z).reciprocal()*c.exp(c.mpf('13.4'))
        d=[defects['rows'][str(i)] for i in range(1,6)]
        inverse=certify(c,weights(c,wraw),d,(Am*Am).reciprocal())
        if not inverse['certified']:raise ValueError('new uniform implicit five-moment inverse failed')
        save('five_bump_inverse',inverse)
        summary=dict(new_schedule_sha256=calc.new_schedule_sha256,state_sha256=calc.state_hash,
            center_family=calc.center_family,completed_core_degree=calc.degree,artifacts=produced,
            new_core_to_R110_exit_and_functional_moment_inverse_generated=True,
            local_uniform_implicit_C1_family_exists=inverse['uniform_implicit_C1_family_exists'],
            new_whole_transition_cone_certified=False,new_heat_exterior_matched=False,
            original_parameter_errors_enclosed=False,whole_axis=False,temporal_recursion=False,
            input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(summary)),indent=2)+'\n',encoding='utf-8')
        print('Fresh Md1.1 transition/inverse pipeline complete',flush=True)
        return summary

if __name__=='__main__':run()

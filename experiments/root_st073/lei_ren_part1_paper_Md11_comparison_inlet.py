"""Controlled fresh-source core inlet and first comparison transition."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_Md11_comparison_jets import Md11ComparisonJets,STEM,TAIL
from lei_ren_part1_paper_Md11_pressure_datum import pressure_jets
from lei_ren_part1_paper_candidate_gauge_core import _unpack_fixed
from lei_ren_part1_paper_interval_comparison_enclosure import integrate,core_box
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    calc=Md11ComparisonJets(4);c=calc.ctx
    state_name=STEM+'_state.json';state=json.loads((HERE/state_name).read_bytes())
    datum=json.loads((HERE/'lei_ren_part1_paper_Md11_pressure_datum.json').read_bytes())
    with mp.workdps(calc.precision+40):
        pressure=pressure_jets(c,c.mpf(calc.center_family),127,datum)
        fixed=_unpack_fixed(c,state['fixed'])
        if any(a._mpi_!=b._mpi_ for a,b in zip(fixed['P0_Z_taylor'],pressure['physical_pressure_coefficients'])):
            raise ValueError('fresh128 pressure coefficients not installed into core')
        old=json.loads((HERE/'lei_ren_part1_paper_candidate_interval_core_Z049_Z051_state.json').read_bytes())
        oldfixed=_unpack_fixed(c,old['fixed'])
        changed=sum(a._mpi_!=b._mpi_ for a,b in zip(fixed['P0_Z_taylor'],oldfixed['P0_Z_taylor']))
        if changed==0:raise ValueError('new core pressure coincides with old saved datum; investigate source')
        core=core_box(calc,c.mpf(4));comparison=integrate(calc,cells=16)
        result=dict(new_schedule_sha256=calc.new_schedule_sha256,state_sha256=calc.state_hash,
            completed_radial_order=calc.degree,center_family=calc.center_family,
            all128_pressure_coefficients_match_new_datum=True,pressure_coefficients_changed_from_old=changed,
            analytic_tail_gate=calc.tail,controlled_core_inlet=core,controlled_comparison=comparison,
            fresh_source_core_and_comparison_generated=True,old_core_coefficients_reused=False,
            new_exit_bridge_generated=False,new_five_moment_inverse_generated=False,
            whole_axis_finite_core=False,whole_outer_profile_matched=False,
            full_NS_residual_validated=False,temporal_recursion=False)
        names=(Path(__file__).name,'lei_ren_part1_paper_Md11_comparison_jets.py',state_name,TAIL,
            'lei_ren_part1_paper_Md11_pressure_datum.json','lei_ren_part1_paper_Md11_pressure_datum.py',
            'lei_ren_part1_paper_interval_comparison_enclosure.py','lei_ren_part1_paper_interval_comparison_jets.py',
            'lei_ren_part1_paper_candidate_interval_core_Z049_Z051_state.json')
        result['input_hashes']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names}
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
        print('Fresh Md1.1 core/inlet generated; degree',calc.degree,'changed pressure rows',changed,
            'controlled first comparison16cells',flush=True)
        return result

if __name__=='__main__':run()

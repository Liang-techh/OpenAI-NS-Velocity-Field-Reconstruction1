"""Same-source exact meridional forcing support at the true Rv endpoint.

The actual selected pulse has an empty backward end-bump tail near Rv.
The original flatten_mixed factory gives Uz=Ur=0 in every future chart.
Angular velocity/pressure and their heat tails are not declared zero.
"""
import ast
import copy
import gzip
import inspect
import json
from pathlib import Path
import textwrap

import lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields as physical
import lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius as radius
import lei_ren_part1_paper_compliant_flatten_mixed_C4 as mixed
import lei_ren_part1_paper_compliant_axial_pulse_field as axial
import lei_ren_part1_paper_compliant_power_angular_C4 as power
import lei_ren_part1_paper_compliant_steep_waiting_C4 as steep
import lei_ren_part1_paper_compliant_collar_Gamma_C4 as heat
import lei_ren_part1_paper_compliant_current_original_n1_inner_operator as operator

HERE,PREFIX,sha=physical.HERE,physical.PREFIX,physical.sha
NAME=PREFIX+'current_original_n1_pressure_forcing_support.json'
RECEIPT=PREFIX+'current_original_n1_pressure_forcing_support_check.json'
GATE='current_original_n1_known_pressure_forcing_zero_from_Rv_to_infinity_certified'
OPEN=physical.OPEN
FUTURE=('flatten','outer_power','outer_angular','steep_entry','steep_power',
        'steep_exit','waiting','heat_collar','heat_exterior')


def read(name):
    data=(HERE/name).read_bytes()
    return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)


def load_inputs():
    reports={};hashes={};family=None
    for module in (physical,radius.post.selected,radius.post,radius):
        receipt=read(module.RECEIPT)
        gates=module.GATES if hasattr(module,'GATES') else (module.GATE,)
        if not receipt['all_passed'] or not all(receipt[k] for k in gates):
            raise ValueError('Accepted current pressure-support source required '+module.RECEIPT)
        if family is None:family=receipt['source_family']
        if receipt['source_family']!=family:raise ValueError('Same original leading source family required')
        for name,digest in receipt['input_hashes'].items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Conflicting pressure-support source '+name)
            hashes[name]=digest
        for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):hashes[name]=sha(name)
        reports[module.NAME]=read(module.NAME)
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    for name,digest in hashes.items():
        if sha(name)!=digest:raise ValueError('Pressure-support source changed '+name)
    return family,reports,hashes


def source_proof():
    """Bind literal zero field definitions and actual future caller chain."""
    def assignment(fn,target):
        tree=ast.parse(textwrap.dedent(inspect.getsource(fn)))
        values=[node.value for node in ast.walk(tree) if isinstance(node,ast.Assign)
                and any(ast.unparse(item)==target for item in node.targets)]
        if len(values)!=1:raise ValueError('Unique actual source assignment required '+target)
        return ast.unparse(values[0])
    if assignment(mixed.flatten_mixed,'zero')!='theta * 0':raise ValueError('Actual zero field definition changed')
    rows=assignment(mixed.flatten_mixed,'rows')
    expected='{UZ: [zero for _ in range(5)], UT: theta_rows, UR: [zero for _ in range(5)], P: prows}'
    if ast.dump(ast.parse(rows,mode='eval'))!=ast.dump(ast.parse(expected,mode='eval')):
        raise ValueError('Actual flatten/future velocity derivative rows changed')
    caller=assignment(power.CompliantPowerAngularC4._packet,'mixed')
    if not caller.startswith('flatten_mixed('):raise ValueError('Actual outer packet must use original mixed factory')
    if assignment(steep.CompliantSteepWaitingC4.packet,'result').split('(')[0]!='self.outer._packet':
        raise ValueError('Actual steep packet must retain outer mixed factory')
    if assignment(heat.CompliantCollarGammaC4.packet,'out').split('(')[0]!='self.outer._packet':
        raise ValueError('Actual heat packet must retain outer mixed factory')
    if assignment(axial.CompliantAxialPulseField.end,'ell')!="c.mpf('.15')":
        raise ValueError('Actual end bump width changed')
    binding=physical.current.ast_binding
    return dict(actual_Uz_and_Ur_y0_to_y4_definitions_exact_zero=True,
        actual_mixed_zero_definition=rows,actual_future_outer_call=caller,
        actual_end_bump_centers=[-3,-1],actual_end_bump_half_width='3/20',
        open_empty_end_bump_collar='-17/20 < log(R/Rv) <= 0',
        end_backward_linear_mean_zero_by_empty_source_support_not_numeric_cancellation=True,
        source_bindings={key:binding(fn) for key,fn in {
            'mixed_factory':mixed.flatten_mixed,'native_end':axial.CompliantAxialPulseField.end,
            'backward_weights':axial.backward_bump_weights,'outer_packet':power.CompliantPowerAngularC4._packet,
            'steep_packet':steep.CompliantSteepWaitingC4.packet,
            'heat_packet':heat.CompliantCollarGammaC4.packet,
            'heat_collar':heat.CompliantCollarGammaC4.collar,
            'heat_exterior':heat.CompliantCollarGammaC4.exterior,
            'original_n1_known_forcing':operator.known_forcing,
            'radius_maps':radius.CurrentOriginalRpSegmentedRadius._map}.items()})


class CurrentOriginalN1PressureForcingSupport:
    """Receipt-bound function theorem; no enormous endpoint is materialized."""
    def __init__(self,require_checked=True):
        self.family,self.reports,self.hashes=load_inputs();self.proof=source_proof()
        data=self.reports[radius.NAME]
        self.nodes=copy.deepcopy(data['exact_function_graph_nodes'])
        self.maps=copy.deepcopy(data['current_absolute_radius_maps'])
        if not set(FUTURE)<=set(self.maps):raise ValueError('Complete actual future radius routes required')
        self.endpoint=dict(definition='Rv=exp(original_logRp+13/original_mu)',
            original_logRp_node=data['original_logRp_node'],
            original_pulse_term_node=self.maps['flatten']['pulse_term'],
            exact_radius_map=self.maps['flatten'],native_flatten_coordinate=0,
            source_family=self.family,angular_and_pressure_not_zeroed=True)
        self.acceptance_loaded=False
        if require_checked:
            receipt=read(RECEIPT)
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Pressure forcing support receipt/scope differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Pressure-support input changed '+name)
                self.hashes[name]=digest
            self.acceptance_loaded=True
        self._frozen=physical.fingerprint((self.family,self.hashes,self.proof,self.nodes,self.maps,self.endpoint))

    def assert_graph(self):
        if physical.fingerprint((self.family,self.hashes,self.proof,self.nodes,self.maps,self.endpoint))!=self._frozen:
            raise ValueError('Same-source pressure-support theorem changed')

    def future_zero(self,chart):
        self.assert_graph()
        if chart not in FUTURE:raise ValueError('Only actual original future charts beyond Rv have this theorem')
        return dict(source_family=copy.deepcopy(self.family),chart=chart,
            whole_original_chart=True,exact_zero_ordinary_Z_rows={key:[0,0] for key in ('U','q','q_R','S','N1p')},
            exact_zero_unweighted_and_R_weighted_forcing_integrals=[0,0],
            actual_endpoint=copy.deepcopy(self.endpoint),
            full_velocity_support_or_angular_heat_pressure_closure_implied=False,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))


def run():
    owner=CurrentOriginalN1PressureForcingSupport(require_checked=False)
    raw=dict(source_family=owner.family,proof=owner.proof,endpoint=owner.endpoint,
        future_charts={chart:owner.future_zero(chart) for chart in FUTURE},
        input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes((json.dumps(raw,separators=(',',':'))+'\n').encode())
    return owner

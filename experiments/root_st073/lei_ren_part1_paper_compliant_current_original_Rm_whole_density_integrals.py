"""Actual Rm active-to-terminal five-density C0/Z integral operators.

The source atlas follows the exact angular-control support boundaries.
Resolved contributions compose with the accepted terminal operator. Any
unresolved signed-source cell remains an explicit unknown integral, never
a zero contribution; the finite-N Rm inlet also stays unsupplied.
"""
import ast
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_terminal_density_integrals as previous

fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_whole_density_integrals.json.gz'
RECEIPT=PREFIX+'current_original_Rm_whole_density_integrals_check.json'
GATE='original_actual_Rm_active_source_atlas_and_partial_C0_Z_Duhamel_contributions_installed'
RATES=previous.RATES
ACTIVE_PARTITION=((1,1),(49,40),(51,40),(59,40),(61,40),(69,40),(71,40))
JOIN=(71,40)


def compile_whole_patch_integrator():
    """Keep the accepted integral graph; extend only its chart lower bound."""
    tree=ast.parse(Path(previous.__file__).read_text(encoding='utf8'))
    nodes=[copy.deepcopy(next(n for n in tree.body if getattr(n,'name',None)==name))
           for name in ('coordinate','integrate_source_cell')]
    edits=[]
    class Chart(ast.NodeTransformer):
        def visit_Compare(self,node):
            self.generic_visit(node)
            if ast.dump(node)==ast.dump(ast.parse('q<Fraction(71,40)',mode='eval').body):
                node.comparators[0]=ast.parse('Fraction(1,1)',mode='eval').body
                edits.append('source chart lower endpoint 71/40 -> 1')
            return node
        def visit_Constant(self,node):
            if node.value=='Actual terminal support x>=71/40 required':
                return ast.copy_location(ast.Constant('Actual Rm support x>=1 required'),node)
            return node
    nodes[0]=Chart().visit(nodes[0])
    if edits!=['source chart lower endpoint 71/40 -> 1']:
        raise ValueError('Accepted exact source-chart bound changed')
    original=next(n for n in tree.body if getattr(n,'name',None)=='integrate_source_cell')
    if ast.dump(nodes[1])!=ast.dump(original):
        raise ValueError('Original signed integration graph must be unchanged')
    env=dict(vars(previous))
    exec(compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),
                 '<accepted integral; exact Rm chart [1,e]>','exec'),env)
    return env['coordinate'],env['integrate_source_cell'],dict(
        original_signed_integral_defining_AST_unchanged=True,
        exact_chart_lower_endpoint_edit=edits,
        same_source_family_P0_radius_phase_and_width_guards_retained=True,
        original_own_rates_masses_and_suffixes_unchanged=True)


coordinate,integrate_source_cell,INTEGRATION_BINDING=compile_whole_patch_integrator()


def validate_active_partition(c,partition):
    if not isinstance(partition,(tuple,list)) or len(partition)<2:
        raise ValueError('Complete ordered actual active source partition required')
    if partition[0]!=(1,1) or partition[-1]!=JOIN:
        raise ValueError('Actual active partition must cover exactly [1,71/40]')
    values=[coordinate(c,v) for v in partition]
    if any(ep(values[j+1])[0]<=ep(values[j])[1] for j in range(len(values)-1)):
        raise ValueError('Strictly ordered exact active source cells required')
    edges={Fraction(*v) for v in partition}
    if not {Fraction(*v) for v in ACTIVE_PARTITION}.issubset(edges):
        raise ValueError('All six actual angular-control support edges required')
    return values


def midpoint(left,right):
    value=(Fraction(*left)+Fraction(*right))/2
    return value.numerator,value.denominator


class OriginalRmWholeDensityIntegrals:
    mode='actual_Rm_active_source_atlas_resolved_C0_Z_contributions_and_explicit_unresolved_source_cells'
    def __init__(self,dps=500,require_checked=True):
        self.terminal=previous.OriginalRmTerminalDensityIntegrals(dps)
        self.phase=self.terminal.upstream;self.c=self.terminal.c
        self.family=self.terminal.family;self.hashes=dict(self.terminal.hashes);self.cache={}
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual Rm whole-source integral receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def contribution(self,label,partition=ACTIVE_PARTITION,N=257,max_depth=3):
        N=previous.previous.candidate_N(N);validate_active_partition(self.c,partition)
        if type(max_depth) is not int or not 0<=max_depth<=16:
            raise ValueError('Bounded source refinement depth in [0,16] required')
        key=(label,tuple(partition),N,max_depth)
        if key in self.cache:return self.cache[key]
        op=self.phase.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        with mp.workdps(c.dps+40):
            cells=[];unresolved=[];atlas=[];refinements=[]
            def enclosed(left,right,depth):
                source=self.phase.spatial_cell(label,left,right,N)
                if source['actual_original_spatial_Z_density_interface_installed']:
                    cell=integrate_source_cell(f,source,left,right,JOIN,
                        source_family=self.family,P0=op.P0,Rm_factor=op.Rm_factor)
                    cells.append(cell);atlas.append(dict(source_geometry=cell['source_geometry'],
                        density_integral_enclosed=True,resolved_cell_index=len(cells)-1))
                    return
                if depth==max_depth:
                    packet=dict(source_geometry=previous.exact_geometry(left,right),source_family=self.family,
                        actual_source_bound_unresolved_phase_density=source,
                        exact_common_P0_axial5=op.P0,candidate_N=N,
                        missing_integral_formula='integral_left^right exp(-lambda*log(71/40/x))*f_j(x,Z)*dx/x',
                        unresolved_integral_not_zero_or_replaced_by_a_sample=True)
                    unresolved.append(packet);atlas.append(dict(source_geometry=packet['source_geometry'],
                        density_integral_enclosed=False,unresolved_cell_index=len(unresolved)-1))
                    return
                middle=midpoint(left,right)
                refinements.append(dict(exact_left=left,exact_right=right,exact_midpoint=middle,
                    reason='original phase primitive needs narrower actual radial source enclosure'))
                enclosed(left,middle,depth+1);enclosed(middle,right,depth+1)
            for left,right in zip(partition,partition[1:]):enclosed(left,right,0)
            active={name:[sum((cell['actual_cell_to_target_integral_C0_Z'][name][n]
                              for cell in cells),f.scalar(0)) for n in range(2)] for name in RATES}
            tail=self.terminal.contribution(label,N=N)
            if tail['exact_common_P0_axial5'] is not op.P0:
                raise ValueError('Same actual active/terminal P0 source required')
            tail_width=tail['full_original_logarithmic_interval_width']
            whole=previous.affine_transport(f,active,tail['actual_terminal_local_defect_integral_C0_Z'],tail_width)
            active_width=c.ln(coordinate(c,JOIN));full_width=c.mpf(1)
            inlet=self.phase.upstream.upstream.evaluate(label,(1,1))
            if inlet['common_original_P0_axial5'] is not op.P0:
                raise ValueError('Same actual Rm leading inlet/P0 source required')
            actual_partition=[tuple(atlas[0]['source_geometry']['exact_left'])]
            actual_partition.extend(tuple(cell['source_geometry']['exact_right']) for cell in atlas)
            for j,cell in enumerate(atlas):
                if cell['source_geometry']!=previous.exact_geometry(actual_partition[j],actual_partition[j+1]):
                    raise ValueError('Complete adjacent actual active source atlas required')
            if actual_partition[0]!=(1,1) or actual_partition[-1]!=JOIN:
                raise ValueError('Actual active atlas lost source endpoint coverage')
            packet=dict(source_family=self.family,source_frame=label,candidate_N=N,
                requested_active_partition=list(partition),exact_actual_active_partition=actual_partition,
                actual_source_refinements=refinements,complete_active_source_atlas=atlas,
                actual_active_cells_to_join=cells,actual_unresolved_active_cells=unresolved,
                actual_active_resolved_cell_integral_C0_Z=active,
                exact_active_terminal_join=JOIN,
                actual_terminal_local_defect_integral_C0_Z=tail['actual_terminal_local_defect_integral_C0_Z'],
                actual_terminal_source_partition=tail['exact_radial_partition'],
                original_terminal_incoming_to_Rh_decays=tail['original_incoming_to_Rh_decays'],
                actual_whole_patch_known_cell_contribution_C0_Z=whole,
                complete_whole_patch_density_integrals_installed=not unresolved,
                complete_active_density_integrals_installed=not unresolved,
                unknown_source_integral_count=len(unresolved),
                whole_patch_history_recipe='deltaH_j(Rh)=exp(-lambda_j)*deltaH_j(Rm)+known_j+exp(-lambda_j*log(e/(71/40)))*sum(unresolved_I_j)',
                exact_common_P0_axial5=op.P0,full_original_logarithmic_interval_width=full_width,
                actual_active_logarithmic_interval_width=active_width,
                actual_terminal_logarithmic_interval_width=tail_width,
                actual_leading_Rm_memory=inlet['common_own_five_histories_axial5'],
                actual_leading_join_memory=tail['actual_leading_inlet_memory'],
                actual_leading_Rh_memory=tail['actual_leading_Rh_memory'],
                finite_N_Rm_incoming_defect_is_unsupplied_affine_argument=True,
                known_cell_integral_enclosures_are_bounds_not_selected_values=True,
                unresolved_integrals_are_explicit_unknown_source_terms=True,
                active_and_terminal_operators_use_same_actual_source_algebra=True,
                pressure_incoming_and_join_decay_exactly_one=True,
                exact_original_support_edges_retained=True,
                two_axial_source_frames_not_whole_Z_provider=True,
                finite_N_complete_prefix_corrected_Rh_or_moment_closure_not_claimed=True,
                **dict.fromkeys(fields.previous.OPEN,False))
            self.cache[key]=packet;return packet

    def transport_supplied_incoming(self,label,incoming,partition=ACTIVE_PARTITION,N=257):
        if incoming is None:raise ValueError('Real finite-N correction at Rm cannot be silently reset')
        packet=self.contribution(label,partition,N)
        if not packet['complete_whole_patch_density_integrals_installed']:
            raise ValueError('Unresolved actual active source integrals cannot be silently omitted')
        op=self.phase.upstream.upstream.upstream.owner(label).op
        with mp.workdps(self.c.dps+40):
            rows=previous.affine_transport(op.flow,incoming,packet['actual_whole_patch_known_cell_contribution_C0_Z'],
                                          packet['full_original_logarithmic_interval_width'])
        return dict(source_family=self.family,source_frame=label,exact_common_P0_axial5=packet['exact_common_P0_axial5'],
            conditional_whole_patch_transported_defect_C0_Z=rows,
            supplied_Rm_incoming_enclosure_not_complete_upstream_prefix_proof=True,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalRmWholeDensityIntegrals(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=owner.contribution(label)
        print('Actual Rm active source atlas',label,len(frames[label]['actual_active_cells_to_join']),
              'resolved cells',len(frames[label]['actual_unresolved_active_cells']),'unresolved cells',flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_integral_chart_binding=INTEGRATION_BINDING,frames=fields.serialized(frames),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),
                                      compresslevel=6,mtime=0));return report


if __name__=='__main__':run()

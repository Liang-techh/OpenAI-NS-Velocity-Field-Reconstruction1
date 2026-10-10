"""Actual current C2 density, own-rate transport and five target functions.

Differentiates the signed defining functions at fixed original radial
coordinate and global phase. Directed C2 target ranges, quantitative C2
repair tails and a point-field oracle are separate, uninstalled tasks.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C2_phase_primitives as phase

current=phase.current;source=current.source;HERE,PREFIX,sha=phase.HERE,phase.PREFIX,phase.sha
NAME=PREFIX+'current_original_C2_density_transport.json.gz'
RECEIPT=PREFIX+'current_original_C2_density_transport_check.json'
GATES=('current_original_17_chart_C2_signed_density_functions_installed',
    'current_original_17_chart_C2_own_rate_transport_and_five_target_functions_installed')
OPEN=phase.OPEN


def coefficient_C2(g,E,V,A,B,N):
    """Full exp/exprel, including both second-order mixed products."""
    alg=phase.C2Algebra(g);invN=g.quotient(g.one,N,'same exact positive integer N')
    old,F1=current.algebra.coefficient_pairs(g,source.C1Function(E.value,E.Z),
        source.C1Function(V.value,V.Z),source.C1Function(A.value,A.Z),source.C1Function(B.value,B.Z),N)
    arg=g.mul(A.value,invN);exprel=g.unary('exprel',arg);exponential=g.unary('exp',arg)
    FZZ=g.add(g.mul(E.ZZ,A.value,exprel),g.mul(g.constant(2),E.Z,A.Z,exponential),
        g.mul(E.value,exponential,g.add(A.ZZ,g.mul(A.Z,A.Z,invN))))
    F=phase.C2Function(F1.value,F1.Z,FZZ)
    EF,VF,EB,VB=[alg.mul(a,b) for a,b in ((E,F),(V,F),(E,B),(V,B))]
    F2,B2=alg.mul(F,F),alg.mul(B,B);zero=alg.fixed(0)
    first=dict(m=B,h=F,k=alg.add(VF,EB),e=alg.add(alg.scale(VB,2),alg.neg(EF)),p=EF)
    second=dict(m=zero,h=zero,k=alg.mul(F,B),e=alg.add(B2,alg.neg(alg.scale(F2,'1/2'))),p=alg.scale(F2,'1/2'))
    coefficients={order:{key:phase.C2Function(old[order][key].value,old[order][key].Z,q.ZZ)
        for key,q in rows.items()} for order,rows in ((-1,first),(-2,second))}
    return coefficients,F


def build(field):
    """Retain exact C1 window handles; append only genuine second rows."""
    built=field.built;g=built['graph'];alg=phase.C2Algebra(g);N=built['N'];ref=lambda i:source.FunctionRef(g,i)
    history={key:alg.fixed(0) for key in current.RATES};windows=[]
    for chart in phase.CHARTS:
        old=field.windows[chart];primitive=field.functions[chart];roots=primitive['source_roots']
        if old['quiet_local_source_zero']:
            density={key:alg.fixed(0) for key in current.RATES};coefficients=None;F=None
        else:
            coefficients,F=coefficient_C2(g,roots['E'],roots['V'],primitive['A'],primitive['B'],N)
            invN=alg.fixed(g.quotient(g.one,N,'same exact positive integer N'))
            density={key:alg.add(coefficients[-1][key],alg.mul(invN,coefficients[-2][key])) for key in current.RATES}
            density={key:phase.C2Function(*map(ref,old['density'][key]),q.ZZ) for key,q in density.items()}
        before=history;after={};contributions={}
        for key,rate in current.RATES.items():
            kernel=g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(ref(old['right_offset']),ref(old['offset'])))))
            pieces=old['actual_current_source_cells'];q=density[key].ZZ
            zz=g.zero if q==g.zero else g.add(*(current.controls.integral(g,g.mul(kernel,q,ref(old['Jacobian'])),
                'native_'+chart,ref(p['lower']),ref(p['upper']),
                measure='original dlogR: Jacobian exactly once',exact_native_left=p['exact_left'],
                exact_native_right=p['exact_right'],actual_current_chart=chart,
                N_scaled_density=True,ordinary_slow_Z_derivative_order=2,
                integration_endpoints_kernel_radius_and_global_phase_Z_independent=True) for p in pieces))
            contributions[key]=phase.C2Function(*map(ref,old['contributions'][key]),zz)
            after[key]=phase.C2Function(*map(ref,old['outgoing'][key]),
                g.add(g.mul(ref(old['memory'][key]),before[key].ZZ),zz))
        history=after
        windows.append(dict(chart=chart,signed_coefficients=coefficients,F_N=F,density=density,
            incoming=before,contributions=contributions,outgoing=history,
            original_own_rate_memory=old['memory'],original_physical_measure_and_partitions=pieces,
            predecessor_ZZ_memory_retained=True,quiet_local_zero_does_not_reset_incoming=True))
    oldA=built['amplitude'];definition=dict(g.nodes[oldA.value.node]);definition['Z_order']=2
    definition.update(ordinary_slow_Z_derivative_order=2,Taylor_coefficient_factorial=2,
        source_projection_binding=current.ast_binding(phase.CurrentC2PhasePrimitives.source_packet),
        source_projection='same actual terminal O3_power E Taylor coefficient[2] times 2',
        derivative_of_range_endpoint=False)
    amplitude=phase.C2Function(oldA.value,oldA.Z,g.node(definition.pop('operation'),**definition))
    A2=alg.mul(amplitude,amplitude);mu=alg.fixed(built['parameters']['mu']);targets={}
    for row,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
        targets[row]=alg.div(history[key],amplitude if degree==1 else A2)
    joint=alg.add(history['k'],alg.neg(alg.mul(amplitude,history['m'])))
    targets[current.controls.ROWS[1]]=alg.div(joint,alg.mul(mu,A2))
    targets={key:phase.C2Function(built['N_scaled_targets'][key].value,built['N_scaled_targets'][key].Z,
        targets[key].ZZ) for key in current.controls.ROWS}
    return dict(windows=windows,terminal_N_scaled_histories=history,actual_terminal_amplitude=amplitude,
        actual_joint_angular_numerator=joint,actual_N_scaled_five_targets=targets,
        original_C1_target_prefix_unchanged=True,global_phase_not_restarted=True,
        all_original_native_partitions_and_own_rates_retained=True,
        differentiation_under_integral='Each original chart is compact with Z-independent endpoints/Jacobian/kernel/phase. '
            'The same source-backed C2 primitives and their strictly positive denominators are continuous on compact source domains; '
            'differentiate the defining finite integral twice before taking ranges.',
        quantitative_current_C2_target_ranges_installed=False,actual_C2_repaired_limit_controls_installed=False)


class CurrentC2DensityTransport:
    def __init__(self,owner=None,phase_field=None,require_checked=True):
        self.phase=phase_field if phase_field is not None else phase.CurrentC2PhasePrimitives(owner=owner)
        if not self.phase.acceptance_loaded:raise ValueError('Accepted actual 17-chart C2 phase source required')
        self.prefix=[dict(n) for n in self.phase.built['graph'].nodes];self.identity=self.phase.identity
        self.functions=build(self.phase);self.hashes=dict(self.phase.hashes)
        for name in (phase.RECEIPT,phase.NAME,Path(__file__).name):self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked C2 target source required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed C2 transport source: '+name)
            self.acceptance_loaded=True

    def target_functions(self):return self.functions['actual_N_scaled_five_targets']


def run(owner=None,phase_field=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC2DensityTransport(owner=owner,phase_field=phase_field,require_checked=False)
        terminal=field.phase.source_packet(('0','.5'),'O3_power',(2,1))
        report=dict(candidate_C2_transport_constructed=True,source_family=field.identity,
            original_phase_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.phase.built['graph'].nodes,
            actual_C2_density_transport_and_target_functions=phase.encoded(field.functions),
            actual_terminal_source_packet=phase.encoded(terminal),
            actual_source_bindings=dict(coefficient_C2=current.ast_binding(coefficient_C2),
                original_C1_coefficient_functions=current.ast_binding(current.algebra.coefficient_pairs),
                original_C1_target_normalization=current.ast_binding(current.algebra.normalize_history),
                actual_transport_and_normalization=current.ast_binding(build)),
            quantitative_current_C2_target_ranges_installed=False,actual_C2_repaired_limit_controls_installed=False,
            **dict.fromkeys(GATES+OPEN,False),input_hashes=field.hashes,
            execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(phase.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    return field


if __name__=='__main__':run()

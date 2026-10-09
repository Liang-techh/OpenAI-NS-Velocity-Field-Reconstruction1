"""Live whole-Z original switches with the genuine same-N R100 incoming.

Six-jet background source and C0/Z finite-N cumulative corrections have
different roles. Only the background enters the original switch owners;
the correction incoming is carried by its own true Volterra memory.
"""
import ast
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import time

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_R100_finite_N as upstream

current,switch=upstream.current,upstream.switch
first,second=current.first,switch.second
original,phase,primitives,bounds,parameters=upstream.original,upstream.phase,upstream.primitives,upstream.bounds,upstream.parameters
HERE,PREFIX,sha,read,bind,ep=upstream.HERE,upstream.PREFIX,upstream.sha,upstream.read,upstream.bind,upstream.ep
RATES,C0,Z=upstream.RATES,upstream.C0,upstream.Z
NAME=PREFIX+'current_original_whole_Z_switch_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_switch_finite_N_check.json'
GATE='current_original_whole_Z_same_N_R100_R110_switch_source_and_actual_history_transport_enclosed'
OPEN=upstream.OPEN
PARTITION=((0,1),(1,2),(1,1))
CHARTS=('first_switch','second_switch','post_power')


def decode(f,row):
    if isinstance(row,dict):
        if 'formal_positive_scale' in row:return first.endpoint.restore_row(f,row)
        if 'lower_exact_mpf_tuple' in row:return read(f.c,row)
        return {key:decode(f,value) for key,value in row.items()}
    if isinstance(row,list):return [decode(f,value) for value in row]
    return row


def replay_prefix(callback):
    tree=ast.parse(Path(switch.__file__).read_text(encoding='utf8'))
    fn=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=='prefix_cell')
    calls=[node for node in ast.walk(fn) if isinstance(node,ast.Call)
           and isinstance(node.func,ast.Name) and node.func.id=='positive_jet']
    if len(calls)!=2 or {node.args[2].value for node in calls}!={'current_first_Dbar','current_second_Dbar'}:
        raise ValueError('Original switch Dbar proof sites changed')
    env={**vars(switch),'positive_jet':callback}
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<unchanged original switch prefix; source-bound Dbar proof>','exec'),env)
    return env['prefix_cell'],dict(original_prefix_AST_sha256=hashlib.sha256(ast.dump(fn).encode()).hexdigest(),
        unchanged_original_arithmetic=True,only_original_positive_Dbar_proof_callback_supplied=True)


class WholeZSwitchFiniteN:
    def __init__(self,dps=500):
        self.background=current.WholeZMicroMacro(dps)
        self.c=c=self.background.c;self.identity=self.background.identity
        self.hashes=dict(self.background.hashes)
        receipt=json.loads((HERE/upstream.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(upstream.GATE) or receipt['source_family']!=self.identity:
            raise ValueError('Checked same-source whole-Z finite-N R100 incoming required')
        for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,upstream.RECEIPT,sha(upstream.RECEIPT))
        self.saved=json.loads(gzip.decompress((HERE/upstream.NAME).read_bytes()))
        self.N=phase.candidate_N(self.saved['candidate_N'])
        if self.saved['source_family']!=self.identity or self.N!=receipt['candidate_N']:
            raise ValueError('Same actual source and candidate N required')
        self.eta_log=read(c,self.saved['source_owned_left_inlet_proof']['selected_eta_log'])
        self.dstar_log=read(c,self.saved['source_owned_left_inlet_proof']['selected_dstar_log'])
        self.saved_cells={tuple(row['exact_Z_cell']):row for row in self.saved['source_cells']}
        if set(self.saved_cells)!=set(current.source.CELLS):raise ValueError('Complete original whole-Z incoming cover required')
        self.bindings=switch.source_bindings()
        for name,digest in self.bindings['original_first_second_and_post_switch_control_bindings']['input_hashes'].items():
            bind(self.hashes,name,digest)
        comparison=self.background.records['global_exit_certificate']['proof']['comparison']
        if 'Dbar>=3.5on100..110' not in comparison.replace(' ',''):
            raise ValueError('Same original stronger switch-domain Dbar theorem required')
        self.logDlower=c.ln(c.mpf('3.5'))
        self.positive_theorem=dict(original_domain='Z[-1,1],100<=R<=110',original_bound='Dbar>=3.5',
            log_Dbar_positive_lower=self.logDlower,source_identity=self.identity,
            original_mode_attachment=self.bindings,source_width_log=self.background.logh,
            global_exit_sha256=sha(PREFIX+'global_exit_certificate.json'),
            source_function_not_selected_cap=True)
        self.prefix,self.prefix_binding=replay_prefix(self.positive_Dbar)
        self.owners={};self.cache={}
        for module in (upstream,switch,first,second):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def positive_Dbar(self,f,rows,role):
        if role not in ('current_first_Dbar','current_second_Dbar'):
            raise ValueError('Only genuine switch comparison Dbar roles admitted')
        if not any(op.flow is f for op in self.owners.values()) or f.c is not self.c:
            raise ValueError('Live original whole-Z mode owner required')
        parameters.same_source(f,rows)
        raw=rows[0];hi=ep(raw.coefficient)[1]
        upper=ep(raw.scale.evaluate()+self.c.ln(self.c.mpf(hi)))[1] if hi>0 else -mp.inf
        if upper<ep(self.logDlower)[0]:raise ValueError('Original raw Dbar cover contradicts switch theorem')
        refined=raw.positive_intersection(ep(self.logDlower)[0])
        return {**parameters.positive_source(f,refined,role),
            'source_log_lower':self.c.mpf(ep(self.logDlower)[0]),'original_raw_source_row':raw,
            'original_analytic_theorem':self.positive_theorem,'theorem_backed_positive_intersection':True}

    def owner(self,ends):
        key=tuple(ends)
        if key in self.owners:return self.owners[key]
        if key not in self.saved_cells:raise ValueError('Actual checked whole-Z incoming cell required')
        op=self.background.owner(ends);f,c=op.flow,op.c
        if ep(c.ln(c.mpf(11)/10)-2*f.hupper)[0]<=0:
            raise ValueError('Both actual switch windows must lie below110')
        first_op=first.FirstSwitchFunctions(f,op.axis['normalized_actual_fields'],op.axis['normalized_actual_own_six_moments'])
        at=first_op.evaluate((1,1))
        second_op=second.SecondSwitchFunctions(first_op,at['actual_fields'],at['actual_six_histories'])
        incoming=self.saved_cells[key]
        saved_P0=decode(f,incoming['exact_common_P0_axial5'])
        if len(saved_P0)!=len(op.reference.P0) or any(
            live.coefficient._mpi_!=saved.coefficient._mpi_
            or live.scale.powers!=saved.scale.powers or live.scale.offset._mpi_!=saved.scale.offset._mpi_
            or live.zero!=saved.zero for live,saved in zip(op.reference.P0,saved_P0)):
            raise ValueError('Same exact live original analytic P0 tuples required')
        result=SimpleNamespace(flow=f,c=c,background=op,first=first_op,second=second_op,
            incoming=decode(f,incoming['actual_current_R100_correction_C0_Z']),
            incoming_background=decode(f,incoming['actual_original_R100_background_C0_Z']))
        for name,rows in result.incoming.items():
            if len(rows)!=2:raise ValueError('Genuine finite-N C0/Z correction pair required')
            parameters.same_source(f,rows)
        self.owners[key]=result;return result

    def query(self,ends,chart,left,right=None):
        if chart not in CHARTS:raise ValueError('Actual original switch chart required')
        right=left if right is None else right;key=(tuple(ends),chart,left,right)
        if key in self.cache:return self.cache[key]
        op=self.owner(ends);f,c=op.flow,op.c
        l,r=switch.long.fraction(left),switch.long.fraction(right)
        if r<l:raise ValueError('Ordered original source phase cell required')
        t=switch.downstream.scalar_hull(c,c.mpf(l.numerator)/l.denominator,c.mpf(r.numerator)/r.denominator)
        source=(switch.post_cell(op.first,op.second,t) if chart=='post_power' else
                self.prefix(op.first,op.second,chart,left,right))
        if chart!='post_power':
            # The original a identity uses the same Dbar modes. Refine its
            # value by that identity while preserving every derivative.
            lower=self.logDlower+f.logs[0]
            if chart=='second_switch':lower=c.mpf(min(ep(lower)[0],ep(c.ln(c.mpf('.8')))[0]))
            raw_a=source['correlated_a_axial5'][0]
            source['correlated_a_axial5'][0]=raw_a.positive_intersection(ep(lower)[0])
            source['phi_y']=f.scale(f.multiply(source['correlated_a_axial5'],source['fields']['phi']),-c.mpf('.5'))
            source['original_a_value_before_source_identity_intersection']=raw_a
            source['original_a_positive_identity_lower']=lower
        proxy,raw,recovered,a=switch.recover_source(op.background.reference,op.background.axis,source)
        roots,qr,proof=switch.general_quotients(proxy,recovered,a,self.eta_log)
        got=primitives.all_u_primitive_bounds(f,roots,qr,self.dstar_log,c.mpf([0,1]))
        if qr[C0].zero and qr[Z].zero:
            got['values']={name:f.scalar(0) for name in got['values']}
            got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
        original_A=got['values']['A']
        log_A_upper=None if original_A.zero else ep(original_A.record()['log_absolute_upper'])[1]
        if log_A_upper is not None and ep(c.ln(self.N))[0]<=log_A_upper:
            raise ValueError('Same incoming N fails actual switch exponent budget at '+str(key)
                             +'; logAupper='+mp.nstr(log_A_upper,18))
        got=upstream.bounded_exponent_cover(f,got,self.N)
        E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
        density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],self.N)
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,chart=chart,
            exact_left=list(left),exact_right=list(right),exact_common_P0_axial5=op.background.reference.P0,
            actual_background_source=source,original_raw_source=raw,original_generic_source=recovered,
            original_full_source_quotients=proof,original_roots=roots['roots'],original_q_C0_Z=qr,
            original_primitive_values=got['values'],original_primitive_proof=got['record'],
            actual_original_A_log_absolute_upper=None if log_A_upper is None else c.mpf(log_A_upper),
            original_signed_five_density_C0_Z=density,
            actual_phase_definition='frac(N*(log100+physical_log_radius_over100-logRa-hb*s_c/2))',
            actual_phase_Z_exact_zero=True,actual_phase_full_period_cover=c.mpf([0,1]),
            phase_average_cancellation_not_claimed=True,real_same_N_R100_correction_supplied=True)
        self.cache[key]=result;return result

    def contribution(self,ends):
        op=self.owner(ends);f,c=op.flow,op.c
        incoming={name:list(rows) for name,rows in op.incoming.items()}
        windows=[]
        for chart in CHARTS:
            total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
            for left,right in zip(PARTITION,PARTITION[1:]):
                query=self.query(ends,chart,left,right);rows={};weight_rows={}
                density=query['original_signed_five_density_C0_Z']
                for name,rate in RATES.items():
                    mass,decay,tail=switch.own_weights(op.first,chart,left,right,rate)
                    parameters.positive_source(f,mass,'actual_whole_Z_switch_Duhamel_mass')
                    pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*tail)
                          for part in ('kernels','Z_derivatives')]
                    for n,row in enumerate(pair):total[name][n]+=row
                    rows[name]=pair;weight_rows[name]=dict(true_full_mass=mass,true_cell_decay=decay,
                        true_suffix_decay=tail,own_rate=str(rate),true_radial_measure_applied_once=True)
                cells.append(dict(source=query,signed_cell_driver_C0_Z=rows,own_rate_weights=weight_rows))
            memory={name:switch.own_weights(op.first,chart,(0,1),(1,1),rate)[1] for name,rate in RATES.items()}
            inherited={name:[row*memory[name] for row in incoming[name]] for name in RATES}
            outgoing={name:[inherited[name][n]+total[name][n] for n in range(2)] for name in RATES}
            windows.append(dict(chart=chart,source_cells=cells,actual_incoming_correction_C0_Z=incoming,
                incoming_own_rate_memory=memory,actual_local_driver_C0_Z=total,
                actual_retained_incoming_C0_Z=inherited,actual_exit_correction_C0_Z=outgoing))
            incoming=outgoing
            print('Whole-Z same-N R100..R110 window: '+str(ends)+' '+chart,flush=True)
        terminal=self.query(ends,'post_power',(1,1))
        background={name:list(rows[:2]) for name,rows in
                    terminal['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete={name:[background[name][n]+incoming[name][n] for n in range(2)] for name in RATES}
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.background.reference.P0,
            actual_R100_correction_incoming_C0_Z=op.incoming,
            actual_source_windows=windows,actual_R110_correction_C0_Z=incoming,
            actual_original_R110_background_C0_Z=background,actual_R110_complete_own_history_C0_Z=complete,
            actual_R110_background_six_history_source=terminal['actual_background_source']['histories'],
            actual_R110_background_field_source=terminal['actual_background_source']['fields'],
            actual_R110_candidate_velocity_C0_Z=terminal['original_signed_five_density_C0_Z']['velocities'],
            R110_candidate_velocity_is_same_N_full_phase_enclosure=True,
            exact_R110_radius=f.scalar(110),total_log_radius_length=c.ln(c.mpf(11)/10),
            exact_total_length_identity='hb+hb+(log(1.1)-2hb)=log(1.1)',
            exact_total_memory_identity='exp(-rate*log(1.1)); pressure rate0 memory1',
            actual_R100_incoming_not_zeroed_or_replaced_by_background=True,
            higher_background_jets_not_fabricated_from_C0_Z_correction=True,
            **dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZSwitchFiniteN()
        cells=[original.serialized(owner.contribution(ends)) for ends in current.source.CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_domain=['-1','1'],exact_Z_partition=[list(cell) for cell in current.source.CELLS],
            source_cells=cells,original_source_bindings=owner.bindings,
            original_prefix_replay_binding=owner.prefix_binding,original_switch_positive_Dbar_theorem=owner.positive_theorem,
            same_live_source_basis_P0_and_current_N_with_real_R100_incoming=True,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual whole-Z original switch/post-power source functions and same-N finite correction '
                  'R100..R110 with real incoming, five signed C0/Z drivers and true own-rate memory. '
                  'Whole-Z downstream Rm/Rref/Rc, terminal repair/global N/cone/heat and recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(current.source.bridge._encode(report),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z genuine same-N finite correction prefix reaches R110',flush=True)
    return report


if __name__=='__main__':run()

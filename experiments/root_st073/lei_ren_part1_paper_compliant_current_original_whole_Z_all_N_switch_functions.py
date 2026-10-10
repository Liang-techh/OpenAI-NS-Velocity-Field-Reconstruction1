"""Original all-N switch functions and normalized R100 -> R110 transport.

Only the hash-bound all-N bridge incoming is consumed. The switch leading
source, inverse functions, physical measures and P0 remain unchanged.
"""
from fractions import Fraction
import ast
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_source_functions as previous
import lei_ren_part1_paper_compliant_current_original_whole_Z_switch_signed_functions as signed_switch

HERE,PREFIX,sha,bind,ep=previous.HERE,previous.PREFIX,previous.sha,previous.bind,previous.ep
CELLS,RATES,C0,Z,OPEN=previous.CELLS,previous.RATES,previous.C0,previous.Z,previous.OPEN
Pair,add,scale,record=previous.Pair,previous.add,previous.scale,previous.record
original=signed_switch.original
CHARTS,PARTITION=signed_switch.CHARTS,signed_switch.PARTITION
NAME=PREFIX+'current_original_whole_Z_all_N_switch_functions.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_all_N_switch_functions_check.json'
GATE='current_original_whole_Z_full_phase_all_N_switch_C1_functions_and_normalized_R110_transport_installed'


def source_binding():
    functions=(original.prefix_cell,original.post_cell,original.recover_source,
        original.general_quotients,original.own_weights,signed_switch.inverse_identity_branch)
    return dict(original_function_AST_sha256={fn.__module__+'.'+fn.__name__:
        hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(fn)))).encode()).hexdigest()
        for fn in functions},unchanged_source_function_entry_points=True,
        original_switch_inverse_identity='A=a/2*(phi-psi_fraction)',
        genuine_original_inverse_Z_phi_derivative_rows_retained=True,
        actual_global_phase='frac(N*(log100+physical_log_radius_over100-logRa-hb*s_c/2))',
        original_radius_maps=dict(first_switch='log100+h*s',second_switch='log100+h*(1+s)',
            post_power='log100+2*h+(log(11/10)-2*h)*t'),
        full_phase_is_outer_parameter_cover_not_claimed_physical_image=True)


class WholeZAllNSwitchFunctions(previous.WholeZAllNBridgeSourceFunctions):
    def __init__(self,dps=500):
        super().__init__(dps)
        receipt=json.loads((HERE/previous.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(previous.GATE) \
                or receipt['source_family']!=self.identity or receipt['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Checked same-input all-N bridge incoming required')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed all-N bridge prerequisite: '+name)
            bind(self.hashes,name,digest)
        bind(self.hashes,previous.RECEIPT,sha(previous.RECEIPT))
        saved=json.loads(gzip.decompress((HERE/previous.NAME).read_bytes()))
        if not saved.get(previous.GATE) or saved['source_family']!=self.identity \
                or saved['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Same-source normalized bridge report required')
        self.incoming_rows={tuple(row['exact_Z_cell']):row for row in saved['actual_four_Z_all_N_bridge_transports']}
        if set(self.incoming_rows)!=set(CELLS):raise ValueError('Complete four-Z all-N incoming required')
        self.switch_source=self.owner.switch_source
        self.switch_binding=source_binding()
        self.prefix_binding=self.switch_source.prefix_binding
        for module in (signed_switch,signed_switch.switch,original):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.packet_cache={}

    def coordinate(self,value):
        q=previous.backend.exact(value)
        if not 0<=q<=1:raise ValueError('Actual switch coordinate in [0,1] required')
        return (q.numerator,q.denominator)

    def leading_packet(self,ends,chart,left,right=None):
        key=tuple(ends)
        if key not in CELLS or chart not in CHARTS:raise ValueError('Admitted actual switch source domain required')
        l=self.coordinate(left);r=self.coordinate(left if right is None else right)
        if Fraction(*r)<Fraction(*l):raise ValueError('Ordered actual switch source coordinates required')
        cachekey=(key,chart,l,r)
        if cachekey in self.packet_cache:return self.packet_cache[cachekey]
        # The source-bound switch owner registers the same live flow for its
        # Dbar callback. Its saved raw correction is never consumed here.
        op=self.owner.switch_owner(key);live=self.owner.owner(key);f,c=op.flow,op.c
        low,high=Fraction(*l),Fraction(*r)
        t=original.downstream.scalar_hull(c,c.mpf(low.numerator)/low.denominator,
            c.mpf(high.numerator)/high.denominator)
        source=original.post_cell(op.first,op.second,t) if chart=='post_power' else \
            self.switch_source.prefix(op.first,op.second,chart,l,r)
        if chart!='post_power':
            lower=self.switch_source.logDlower+f.logs[0]
            if chart=='second_switch':lower=c.mpf(min(ep(lower)[0],ep(c.ln(c.mpf('.8')))[0]))
            raw_a=source['correlated_a_axial5'][0]
            source['correlated_a_axial5'][0]=raw_a.positive_intersection(ep(lower)[0])
            source['phi_y']=f.scale(f.multiply(source['correlated_a_axial5'],source['fields']['phi']),-c.mpf('.5'))
            source['original_a_value_before_source_identity_intersection']=raw_a
            source['original_a_positive_identity_lower']=lower
        proxy,raw,recovered,a=original.recover_source(op.background.reference,op.background.axis,source)
        roots,qr,proof=original.general_quotients(proxy,recovered,a,self.switch_source.eta_log)
        packet=dict(source_identity=self.identity,candidate_N=self.N0,exact_Z_cell=list(key),chart=chart,
            exact_left=list(l),exact_right=list(r),exact_common_P0_axial5=live.P0,
            original_switch_independent_P0_axial5=op.background.reference.P0,
            actual_background_source=source,original_raw_source=raw,original_generic_source=recovered,
            original_full_source_quotients=proof,original_roots=roots['roots'],original_q_C0_Z=qr,
            original_eta_log=self.switch_source.eta_log,actual_phase_Z_exact_zero=True,
            identical_source_P0_tuples_rebound_to_live_descendant_object=True,
            fixed_N_density_and_fixed_N_phase_provider_not_called=True,
            first_micro_collar_not_inherited=True)
        self.packet_cache[cachekey]=packet
        return packet

    def primitive_cell(self,ends,chart,left,right=None):
        key=tuple(ends);l=self.coordinate(left);r=self.coordinate(left if right is None else right)
        cachekey=(key,chart,l,r)
        if cachekey in self.cache:return self.cache[cachekey]
        packet=self.leading_packet(key,chart,l,r);op=self.owner.owner(key)
        if packet['source_identity']!=self.identity or packet['candidate_N']!=self.N0 \
                or packet['exact_common_P0_axial5'] is not op.P0:
            raise ValueError('Same actual switch family, prefix N0 and independent P0 required')
        qr=packet['original_q_C0_Z']
        source=dict(q=qr[C0],roots=packet['original_roots'])
        previous.relative.same_source(op.flow,[v for row in source['roots'].values() for v in row.values()]+list(qr.values()))
        _,branches,empty=previous.backend.signed.signed_u_branches(source,self.switch_source.dstar_log)
        alternatives=[]
        def phase_piece(a,b,depth=0):
            phi=self.c.mpf((ep(self.c.mpf(a.numerator)/a.denominator)[0],
                ep(self.c.mpf(b.numerator)/b.denominator)[1]))
            got=[signed_switch.inverse_identity_branch(source,qr,self.switch_source.dstar_log,phi,branch) for branch in branches]
            if any(row['values'] is None for row in got):
                if depth>=5:raise ArithmeticError('Actual full-phase switch inverse needs source/phase refinement')
                middle=(a+b)/2;phase_piece(a,middle,depth+1);phase_piece(middle,b,depth+1);return
            alternatives.extend(dict(exact_phase_left=[a.numerator,a.denominator],
                exact_phase_right=[b.numerator,b.denominator],full_phase_parameter_box=phi,
                branch=branch['name'],primitive_C0_Z_phi=row['values'],original_inverse_proof=row['record'])
                for branch,row in zip(branches,got))
        for i in range(4):phase_piece(Fraction(i,4),Fraction(i+1,4))
        E,V=packet['original_generic_source']['common_velocity_E_axial5'],packet['original_generic_source']['common_velocity_V_axial5']
        result=dict(source_identity=self.identity,source_input_fixed_prefix_N0=self.N0,exact_Z_cell=list(key),
            actual_chart=chart,exact_left=list(l),exact_right=list(r),exact_common_P0_axial5=op.P0,
            original_source_packet=packet,original_E_V_C0_Z=dict(E=Pair(*E[:2]),V=Pair(*V[:2])),
            full_phase_inverse_function_alternatives=alternatives,empty_signed_branches=empty,
            explicit_full_phase_parameter_cover=[0,1],full_cover_contains_original_frac_N_radius_for_every_integer_N=True,
            phase_is_a_range_parameter_not_a_selected_spatial_angle=True,phase_Z_exact_zero=True,
            original_roots_q_velocity_P0_widths_unchanged=True,first_micro_collar_not_inherited=True)
        self.cache[cachekey]=result
        return result

    def normalized_incoming(self,ends):
        op=self.owner.owner(ends);saved=self.incoming_rows[tuple(ends)];f=op.flow
        if saved['source_identity']!=self.identity or saved['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Same fixed leading input and normalized incoming family required')
        decode=lambda row:previous.relative.restore_row(f,row)
        if len(op.P0)!=len(saved['exact_common_P0_axial5']) or any(
                not previous.relative.exact_row(a,decode(b)) for a,b in zip(op.P0,saved['exact_common_P0_axial5'])):
            raise ValueError('Exact same bridge/switch independent P0 required')
        incoming={name:Pair(*[decode(row) for row in rows]) for name,rows in
            saved['actual_uniform_N_scaled_R100_correction_C0_Z'].items()}
        if set(incoming)!=set(RATES) or all(q.value.zero and q.Z.zero for q in incoming.values()):
            raise ValueError('Actual nonzero normalized five-moment R100 incoming required')
        return incoming

    def transport(self,ends):
        op=self.owner.owner(ends);native=self.owner.switch_owner(ends);f=op.flow
        incoming=self.normalized_incoming(ends);initial=dict(incoming);windows=[]
        for chart in CHARTS:
            local={name:Pair(f.scalar(0),f.scalar(0)) for name in RATES};cells=[]
            for left,right in zip(PARTITION,PARTITION[1:]):
                query=self.query(ends,chart,left,right);contributions={};weights={}
                for name,rate in RATES.items():
                    mass,decay,suffix=original.own_weights(native.first,chart,left,right,rate)
                    contribution=scale(query['N_scaled_density_C0_Z'][name],mass*suffix)
                    local[name]=add(local[name],contribution);contributions[name]=contribution
                    weights[name]=dict(original_mass=mass,cell_decay=decay,suffix_decay=suffix)
                packet=query['source']['original_source_packet']
                cells.append(dict(exact_left=list(left),exact_right=list(right),
                    full_phase_inverse_function_alternatives=query['source']['full_phase_inverse_function_alternatives'],
                    original_E_V_C0_Z={k:record(v) for k,v in query['source']['original_E_V_C0_Z'].items()},
                    original_q_C0_Z=packet['original_q_C0_Z'],original_q_Z_scope=packet['original_full_source_quotients']['original_q_Z_scope'],
                    original_positive_chart_log_length=packet['actual_background_source']['window_length'],
                    epsilon=query['epsilon'],N_scaled_density_C0_Z={k:record(v) for k,v in query['N_scaled_density_C0_Z'].items()},
                    actual_N_scaled_contribution_C0_Z={k:record(v) for k,v in contributions.items()},
                    original_own_rate_weights=weights,physical_Jacobian_applied_once=True,
                    all_N_phase_cover_not_fixed_N_phase_reuse=True,first_micro_collar_not_inherited=True))
            memory={name:original.own_weights(native.first,chart,(0,1),(1,1),rate)[1] for name,rate in RATES.items()}
            outgoing={name:add(scale(incoming[name],memory[name]),local[name]) for name in RATES}
            windows.append(dict(actual_chart=chart,actual_all_N_source_cells=cells,
                actual_N_scaled_incoming_C0_Z={k:record(v) for k,v in incoming.items()},
                original_incoming_own_rate_memory=memory,actual_N_scaled_local_C0_Z={k:record(v) for k,v in local.items()},
                actual_N_scaled_outgoing_C0_Z={k:record(v) for k,v in outgoing.items()}))
            incoming=outgoing
            print('Actual all-N normalized switch transport: '+str(ends)+' '+chart,flush=True)
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),fixed_leading_input_N0=self.N0,
            exact_common_P0_axial5=op.P0,epsilon=previous.epsilon_cover(self.c,self.N0),
            actual_uniform_N_scaled_R100_incoming_C0_Z={k:record(v) for k,v in initial.items()},
            actual_all_N_switch_windows=windows,actual_uniform_N_scaled_R110_correction_C0_Z={k:record(v) for k,v in incoming.items()},
            exact_R110_radius=f.scalar(110),exact_total_log_width_identity='h+h+(log(11/10)-2*h)=log(11/10)',
            actual_nonzero_normalized_R100_incoming_not_reset=True,normalized_incoming_not_passed_to_raw_correction_inlet=True,
            original_positive_widths_and_rate_zero_pressure_memory_retained=True,
            fixed_N_local_drivers_not_rescaled_or_hydrated=True,actual_uniform_N_scaled_Rc_targets_installed=False,
            actual_global_frequency_admitted=False,**dict.fromkeys(OPEN,False))


def serialized(value):
    if isinstance(value,Pair):return record(value)
    if isinstance(value,dict):return {('y%d_Z%d'%k if isinstance(k,tuple) else k):serialized(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serialized(v) for v in value]
    return previous.current.serialized(value)


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZAllNSwitchFunctions()
        rows=[owner.transport(ends) for ends in CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,fixed_leading_input_N0=owner.N0,
            exact_Z_partition=CELLS,current_switch_charts=CHARTS,exact_native_partition=PARTITION,
            original_switch_source_binding=owner.switch_binding,original_prefix_replay_binding=owner.prefix_binding,
            actual_four_Z_all_N_switch_transports=rows,actual_full_phase_all_N_switch_driver_functions_installed=True,
            actual_uniform_N_scaled_Rc_targets_installed=False,actual_all_regions_all_N_adapter_installed=False,
            actual_N_squared_averaging_cancellation_proved=False,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual unchanged whole-Z switch leading source and original signed inverse C1 functions '
                'over full phase and epsilon=1/N. Fresh physical-measure local drivers propagate the accepted '
                'nonzero normalized R100 incoming to R110. No uniform Rc target, global frequency, N^-2 '
                'averaging, solved controls, heat/stress, temporal recursion or full NS admission.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(previous.current.encode(serialized(report)),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()

"""Same-input all-N long/reference and analytic-patch correction functions.

Fresh defining source functions transport the accepted normalized R110
incoming through Rm and Rh. Fixed-N local densities/inlets are not used.
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
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_switch_functions as previous
import lei_ren_part1_paper_compliant_current_original_whole_Z_long_Rm_signed_functions as long_signed
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_finite_N as patch_source

HERE,PREFIX,sha,bind,ep=previous.HERE,previous.PREFIX,previous.sha,previous.bind,previous.ep
CELLS,RATES,C0,Z,OPEN=previous.CELLS,previous.RATES,previous.C0,previous.Z,previous.OPEN
Pair,add,scale,record,serialized=previous.Pair,previous.add,previous.scale,previous.record,previous.serialized
core=previous.previous
LONG_CHARTS,LONG_PARTITION=long_signed.CHARTS,long_signed.PARTITION
PATCH_PARTITION=patch_source.PARTITION
CHARTS=(*LONG_CHARTS,'actual_patch')
NAME=PREFIX+'current_original_whole_Z_all_N_long_patch_functions.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_all_N_long_patch_functions_check.json'
GATE='current_original_whole_Z_full_phase_all_N_long_patch_C1_functions_and_normalized_Rh_transport_installed'


def source_binding():
    functions=(patch_source.recover_patch,patch_source.switch.general_quotients,
        patch_source.patch_weights,long_signed.source_module.long.kernel_weight,
        previous.signed_switch.inverse_identity_branch)
    return dict(original_long_source_frontend=long_signed.FRONTEND_BINDING,
        original_patch_function_AST_sha256={fn.__module__+'.'+fn.__name__:
            hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(fn)))).encode()).hexdigest()
            for fn in functions},pure_patch_mixed4_leading_provider_used=True,
        actual_long_partition_source=Path(long_signed.source_module.__file__).name,
        original_generic_inverse_identity='A=a/2*(phi-psi_fraction)',
        same_identity_already_used_by_accepted_long_and_full_signed_patch_query=True,
        original_global_phase='frac(N*(logR-logRa-hb*s_c/2))',
        original_long_radius_maps=dict(long_reshape='log110+T*t',reference='log110+T+gap*t',
            restoration='log110+10*(logC+logP)-8+t',postrestore='log110+10*(logC+logP)-7+t'),
        original_patch_radius='Rm*x',original_patch_measure='dy=dx/x',symbolic_Rh='x=exp(1)',
        full_phase_is_outer_parameter_cover_not_claimed_physical_image=True)


def equivalent_rows(actual,expected):
    if len(actual)!=len(expected):return False
    return all(core.relative.exact_row(a,b) and a.ctx is b.ctx and a.scale.bases is b.scale.bases
        and a.ledger is b.ledger for a,b in zip(actual,expected))


class WholeZAllNLongPatchFunctions(previous.WholeZAllNSwitchFunctions):
    def __init__(self,dps=500):
        super().__init__(dps)
        receipt=json.loads((HERE/previous.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(previous.GATE) \
                or receipt['source_family']!=self.identity or receipt['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Checked same-input all-N R110 incoming required')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed all-N switch prerequisite: '+name)
            bind(self.hashes,name,digest)
        bind(self.hashes,previous.RECEIPT,sha(previous.RECEIPT))
        saved=json.loads(gzip.decompress((HERE/previous.NAME).read_bytes()))
        if not saved.get(previous.GATE) or saved['source_family']!=self.identity \
                or saved['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Same source all-N switch report required')
        self.R110_rows={tuple(row['exact_Z_cell']):row for row in saved['actual_four_Z_all_N_switch_transports']}
        if set(self.R110_rows)!=set(CELLS):raise ValueError('Complete four-Z normalized R110 cover required')
        self.long_patch_binding=source_binding()
        for module in (long_signed,patch_source,patch_source.upstream):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def coordinate(self,value,chart=None):
        if chart=='actual_patch' and value=='Rh':return 'Rh'
        q=core.backend.exact(value)
        if chart=='actual_patch':
            if q<1:raise ValueError('Original patch coordinate x>=1 required')
        elif not 0<=q<=1:raise ValueError('Original long native coordinate in [0,1] required')
        return (q.numerator,q.denominator)

    def leading_packet(self,ends,chart,left,right=None):
        key=tuple(ends)
        if key not in CELLS or chart not in CHARTS:raise ValueError('Admitted original long/patch domain required')
        l=self.coordinate(left,chart);r=self.coordinate(left if right is None else right,chart)
        if l=='Rh' and r!='Rh' or r!='Rh' and l!='Rh' and Fraction(*r)<Fraction(*l):
            raise ValueError('Ordered original long/patch coordinates required')
        cachekey=(key,chart,l,r)
        if cachekey in self.packet_cache:return self.packet_cache[cachekey]
        live=self.owner.owner(key)
        if chart in LONG_CHARTS:
            self.owner.long_owner(key)
            packet=dict(self.owner.source_query(key,chart,l,r))
        else:
            source=self.owner.patch_source
            op=source.owner(key)
            if op.flow is not live.flow or op.c is not live.c or not equivalent_rows(op.P0,live.P0) \
                    or not equivalent_rows([op.Rm_factor],[live.Rm_factor]):
                raise ValueError('Same live patch flow/context/P0/Rm radius required')
            leading=source.source.query(key,l,None if l==r else r)
            if leading['source_identity']!=self.identity or leading['original_P0_normalized_axial5'] is not op.P0 \
                    or leading['actual_patch_finite_N_density_oracle_installed'] is not False:
                raise ValueError('Same original pure leading patch source required')
            proxy,raw,recovered,a,shear=patch_source.recover_patch(op,leading)
            roots,qr,proof=patch_source.switch.general_quotients(proxy,recovered,a,source.eta_log)
            packet=dict(source_identity=self.identity,candidate_N=self.N0,exact_Z_cell=list(key),chart=chart,
                exact_left=l,exact_right=r,exact_common_P0_axial5=live.P0,
                original_patch_independent_P0_axial5=op.P0,actual_background_source=leading,
                original_raw_source=raw,original_generic_source=recovered,original_roots=roots['roots'],
                original_q_C0_Z=qr,original_full_source_quotients=proof,actual_correlated_shear_source=shear,
                original_window_length=self.c.mpf(1),actual_phase_Z_exact_zero=True,
                identical_source_P0_tuples_rebound_to_live_descendant_object=True)
        if packet['source_identity']!=self.identity or packet['candidate_N']!=self.N0 \
                or packet['exact_common_P0_axial5'] is not live.P0 or not packet['actual_phase_Z_exact_zero']:
            raise ValueError('Same fixed leading-input source and independent live P0 required')
        packet.update(fixed_N_density_and_phase_and_inlet_not_called=True,bridge_collar_not_inherited=True)
        self.packet_cache[cachekey]=packet
        return packet

    def primitive_cell(self,ends,chart,left,right=None):
        l=self.coordinate(left,chart);r=self.coordinate(left if right is None else right,chart)
        key=(tuple(ends),chart,l,r)
        if key in self.cache:return self.cache[key]
        packet=self.leading_packet(ends,chart,l,r);op=self.owner.owner(ends);qr=packet['original_q_C0_Z']
        source=dict(q=qr[C0],roots=packet['original_roots'])
        core.relative.same_source(op.flow,[v for row in source['roots'].values() for v in row.values()]+list(qr.values()))
        dstar=self.owner.long_source.dstar_log if chart in LONG_CHARTS else self.owner.patch_source.dstar_log
        _,branches,empty=core.backend.signed.signed_u_branches(source,dstar)
        alternatives=[]
        def phase_piece(a,b,depth=0):
            phi=self.c.mpf((ep(self.c.mpf(a.numerator)/a.denominator)[0],
                ep(self.c.mpf(b.numerator)/b.denominator)[1]))
            got=[previous.signed_switch.inverse_identity_branch(source,qr,dstar,phi,branch) for branch in branches]
            if any(row['values'] is None for row in got):
                if depth>=5:raise ArithmeticError('Actual long/patch full-phase inverse needs source/phase refinement')
                middle=(a+b)/2;phase_piece(a,middle,depth+1);phase_piece(middle,b,depth+1);return
            alternatives.extend(dict(exact_phase_left=[a.numerator,a.denominator],
                exact_phase_right=[b.numerator,b.denominator],full_phase_parameter_box=phi,
                branch=branch['name'],primitive_C0_Z_phi=row['values'],original_inverse_proof=row['record'])
                for branch,row in zip(branches,got))
        for i in range(4):phase_piece(Fraction(i,4),Fraction(i+1,4))
        E,V=packet['original_generic_source']['common_velocity_E_axial5'],packet['original_generic_source']['common_velocity_V_axial5']
        result=dict(source_identity=self.identity,source_input_fixed_prefix_N0=self.N0,exact_Z_cell=list(ends),
            actual_chart=chart,exact_left=l,exact_right=r,exact_common_P0_axial5=op.P0,original_source_packet=packet,
            original_E_V_C0_Z=dict(E=Pair(*E[:2]),V=Pair(*V[:2])),full_phase_inverse_function_alternatives=alternatives,
            empty_signed_branches=empty,explicit_full_phase_parameter_cover=[0,1],
            full_cover_contains_original_frac_N_radius_for_every_integer_N=True,
            phase_is_a_range_parameter_not_a_selected_spatial_angle=True,phase_Z_exact_zero=True,
            original_roots_q_velocity_P0_widths_unchanged=True,bridge_collar_not_inherited=True)
        self.cache[key]=result
        return result

    def window_length(self,ends,chart):
        if chart=='actual_patch':return self.c.mpf(1)
        op=self.owner.long_owner(ends);source=self.owner.long_source
        if source.T._mpi_!=(400*source.A)._mpi_ or source.T._mpi_!=op.long.T._mpi_ \
                or source.T._mpi_!=op.reference.T._mpi_:
            raise ValueError('Original defining T=400*A and live long/reference identity required')
        if op.reference.gap._mpi_!=(10*(op.reference.logC+op.flow.logs[1]/2)-source.T-8)._mpi_:
            raise ValueError('Exact original positive reference gap identity required')
        length=source.T if chart=='long_reshape' else op.reference.gap if chart=='reference' else self.c.mpf(1)
        if chart not in LONG_CHARTS or ep(length)[0]<=0:raise ValueError('Positive original long chart length required')
        return length

    def weights(self,ends,chart,left,right,rate):
        f=self.owner.owner(ends).flow;c=self.c
        if chart=='actual_patch':
            width,suffix,mass,decay,tail=patch_source.patch_weights(f,left,right,rate)
        else:
            l,r=Fraction(*left),Fraction(*right);length=self.window_length(ends,chart)
            width=length*(c.mpf((r-l).numerator)/(r-l).denominator)
            suffix=length*(c.mpf((1-r).numerator)/(1-r).denominator)
            mass,decay,tail=patch_source.long.kernel_weight(f,width,suffix,rate)
        if ep(width)[0]<=0:raise ValueError('Strictly positive original log-radius cell required')
        return dict(original_mass=mass,cell_decay=decay,suffix_decay=tail,
            physical_log_width=width,suffix_log_width=suffix)

    def normalized_R110_incoming(self,ends):
        row=self.R110_rows[tuple(ends)];op=self.owner.owner(ends);f=op.flow
        if row['source_identity']!=self.identity or row['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Same-input normalized R110 correction required')
        decode=lambda value:core.relative.restore_row(f,value)
        if not equivalent_rows(op.P0,[decode(v) for v in row['exact_common_P0_axial5']]):
            raise ValueError('Exact same normalized R110/long P0 required')
        incoming={name:Pair(*[decode(v) for v in values]) for name,values in
            row['actual_uniform_N_scaled_R110_correction_C0_Z'].items()}
        if set(incoming)!=set(RATES) or all(q.value.zero and q.Z.zero for q in incoming.values()):
            raise ValueError('Original nonzero normalized five-moment incoming required')
        return incoming

    def transport(self,ends):
        op=self.owner.owner(ends);f=op.flow;incoming=self.normalized_R110_incoming(ends);initial=dict(incoming)
        windows=[];Rm=None
        for chart in CHARTS:
            partition=PATCH_PARTITION if chart=='actual_patch' else LONG_PARTITION
            length=self.window_length(ends,chart)
            local={name:Pair(f.scalar(0),f.scalar(0)) for name in RATES};cells=[]
            for left,right in zip(partition,partition[1:]):
                query=self.query(ends,chart,left,right);contributions={};weights={}
                for name,rate in RATES.items():
                    weight=self.weights(ends,chart,left,right,rate)
                    contribution=scale(query['N_scaled_density_C0_Z'][name],weight['original_mass']*weight['suffix_decay'])
                    local[name]=add(local[name],contribution);contributions[name]=contribution;weights[name]=weight
                cells.append(dict(exact_left=left,exact_right=right,epsilon=query['epsilon'],
                    full_phase_inverse_function_alternatives=query['source']['full_phase_inverse_function_alternatives'],
                    original_E_V_C0_Z={k:record(v) for k,v in query['source']['original_E_V_C0_Z'].items()},
                    N_scaled_density_C0_Z={k:record(v) for k,v in query['N_scaled_density_C0_Z'].items()},
                    actual_N_scaled_contribution_C0_Z={k:record(v) for k,v in contributions.items()},
                    original_own_rate_weights=weights,physical_Jacobian_applied_once=True,
                    all_N_phase_cover_not_fixed_N_phase_reuse=True,bridge_collar_not_inherited=True))
                print('Actual all-N '+chart+' cell: '+str(ends)+' '+str(left)+' '+str(right),flush=True)
            memory={name:patch_source.long.kernel_weight(f,length,self.c.mpf(0),rate)[1] for name,rate in RATES.items()}
            outgoing={name:add(scale(incoming[name],memory[name]),local[name]) for name in RATES}
            windows.append(dict(actual_chart=chart,actual_original_window_length=length,actual_all_N_source_cells=cells,
                actual_N_scaled_incoming_C0_Z={k:record(v) for k,v in incoming.items()},original_incoming_own_rate_memory=memory,
                actual_N_scaled_local_C0_Z={k:record(v) for k,v in local.items()},
                actual_N_scaled_outgoing_C0_Z={k:record(v) for k,v in outgoing.items()}))
            incoming=outgoing
            if chart=='postrestore':Rm=dict(incoming)
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),fixed_leading_input_N0=self.N0,
            exact_common_P0_axial5=op.P0,epsilon=core.epsilon_cover(self.c,self.N0),exact_Rm_radius=op.Rm_factor,
            exact_Rh_radius=op.Rm_factor*self.c.exp(1),actual_uniform_N_scaled_R110_incoming_C0_Z={k:record(v) for k,v in initial.items()},
            actual_all_N_long_patch_windows=windows,actual_uniform_N_scaled_Rm_correction_C0_Z={k:record(v) for k,v in Rm.items()},
            actual_uniform_N_scaled_Rh_correction_C0_Z={k:record(v) for k,v in incoming.items()},
            exact_long_total_log_length_identity='T+gap+1+1=10*(logC+logP)-6',
            exact_patch_total_log_length_identity='sum log(x_right/x_left)=log(exp(1))=1',
            normalized_correction_only_incoming_background_P0_distinct=True,
            fixed_N_local_densities_inlets_and_phase_not_used=True,actual_uniform_N_scaled_Rc_targets_installed=False,
            actual_global_frequency_admitted=False,**dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZAllNLongPatchFunctions()
        rows=[owner.transport(ends) for ends in CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,fixed_leading_input_N0=owner.N0,
            exact_Z_partition=CELLS,current_long_patch_charts=CHARTS,exact_long_native_partition=LONG_PARTITION,
            exact_patch_partition=PATCH_PARTITION,original_long_patch_source_binding=owner.long_patch_binding,
            actual_four_Z_all_N_long_patch_transports=rows,actual_full_phase_all_N_long_patch_driver_functions_installed=True,
            actual_uniform_N_scaled_Rc_targets_installed=False,actual_all_regions_all_N_adapter_installed=False,
            actual_N_squared_averaging_cancellation_proved=False,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual unchanged whole-Z long/reference and analytic patch source functions over full phase '
                'and epsilon=1/N. Original T/gap/1/1 and dx/x measures carry the accepted correction-only '
                'normalized R110 incoming through Rm and Rh. No uniform Rc target, global frequency, N^-2 '
                'averaging, solved controls, heat/stress, temporal recursion or full NS admission.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(core.current.encode(serialized(report)),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()

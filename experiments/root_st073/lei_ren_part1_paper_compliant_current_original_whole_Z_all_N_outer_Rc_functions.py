"""Pure all-N outer source functions, normalized Rh -> Rc and Rc targets.

Original leading source frontends stop before support/density evaluation.
Only checked N-scaled correction histories are transported; the fixed
leading background and its independent pressure datum remain distinct.
"""
from fractions import Fraction
import ast
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_long_patch_functions as previous
import lei_ren_part1_paper_compliant_current_original_whole_Z_full_signed_transport as outer

HERE,PREFIX,sha,bind,ep=previous.HERE,previous.PREFIX,previous.sha,previous.bind,previous.ep
CELLS,RATES,C0,Z,OPEN=previous.CELLS,previous.RATES,previous.C0,previous.Z,previous.OPEN
Pair,add,scale,record,serialized=previous.Pair,previous.add,previous.scale,previous.record,previous.serialized
core=previous.core
rh,slope,axial,o3=outer.rh,outer.slope,outer.axial,outer.o3
CHARTS,PARTITIONS=outer.OUTER_CHARTS,outer.OUTER_PARTITIONS
NAME=PREFIX+'current_original_whole_Z_all_N_outer_Rc_functions.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_all_N_outer_Rc_functions_check.json'
GATE='current_original_whole_Z_full_phase_all_N_outer_C1_functions_normalized_Rc_and_relative_targets_installed'


def source_frontend(fn,kind):
    """Retain every original leading-source statement before all-u support."""
    node=ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
    stops=[i for i,statement in enumerate(node.body) if isinstance(statement,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='got' for t in statement.targets)]
    if not stops or ast.unparse(node.body[stops[0]].value.func)!='reference.primitives.all_u_primitive_bounds':
        raise ValueError('Original outer source/support boundary changed: '+kind)
    prefix=node.body[:stops[0]]
    proof=dict(original_query_module=Path(inspect.getfile(fn)).name,
        original_query_sha256=sha(Path(inspect.getfile(fn)).name),
        exact_unchanged_leading_frontend_AST_sha256=hashlib.sha256(ast.dump(ast.Module(
            body=prefix,type_ignores=[])).encode()).hexdigest(),
        stop_before_all_u_support_and_fixed_N_density=True,
        original_parent_leading_source_recipe_not_changed=True)
    source_key={'rh':'actual_background_source','slope':'original_closed_O2_slope_background',
        'axial':'original_closed_O2_axial_buffer_background','o3':'original_closed_O3_background'}[kind]
    generic='recovered' if kind=='rh' else 'generic'
    background='background' if kind=='rh' else 'source'
    tail=ast.parse(f'''result=dict(source_identity=self.identity,candidate_N=self.N,
        exact_Z_cell=list(ends),exact_left=list(left),exact_right=list(right),
        exact_common_P0_axial5=op.P0,original_generic_source={generic},
        original_full_source_quotients=proof,original_roots=roots['roots'],
        original_q_C0_Z=qr,actual_phase_Z_exact_zero=True,
        actual_original_source_frontend_binding=FRONTEND_BINDING,
        {source_key}={background})
self.cache[key]=result
return result
''').body
    node.name='pure_original_'+kind+'_frontend';node.body=prefix+tail
    scope={**fn.__globals__,'FRONTEND_BINDING':proof}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),
        '<unchanged pure original '+kind+' leading frontend>','exec'),scope)
    return scope[node.name],proof


FRONTENDS={};FRONTEND_BINDINGS={}
for kind,fn in (('rh',rh.WholeZRhReferenceFiniteN.query),('slope',slope.WholeZO2SlopeFiniteN.query),
        ('axial',axial.WholeZO2AxialBufferFiniteN.query),('o3',o3.WholeZO3RcFiniteN.query)):
    FRONTENDS[kind],FRONTEND_BINDINGS[kind]=source_frontend(fn,kind)


def normalized_targets(f,histories,amplitude,mu,logmu):
    """Direct targets of H=N*D, with no second multiplication by N."""
    if set(histories)!=set(RATES) or len(amplitude)!=2 or any(not isinstance(v,Pair) for v in histories.values()):
        raise ValueError('Five normalized correction C0/Z functions required')
    A,AZ=amplitude
    core.relative.same_source(f,[v for pair in histories.values() for v in (pair.value,pair.Z)]+[A,AZ,mu])
    if ep(A.coefficient)[0]<=0 or ep(mu.coefficient)[0]<=0:
        raise ValueError('Original strictly positive terminal A_rc and mu required')
    c=f.c;logA=c.mpf(ep(A.scale.evaluate()+c.ln(A.coefficient))[0])
    rho=AZ.positive_divide(A,logA);den=A*A;targets={}
    for out,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
        H=histories[key];divisor=A if degree==1 else den
        targets[out]=Pair(H.value.positive_divide(divisor,degree*logA),
            (H.Z-H.value*rho*degree).positive_divide(divisor,degree*logA))
    numerator=Pair(histories['k'].value-A*histories['m'].value,
        histories['k'].Z-AZ*histories['m'].value-A*histories['m'].Z)
    divided=core.relative.repair.ROWS[1]
    targets[divided]=Pair(numerator.value.positive_divide(den*mu,2*logA+logmu),
        (numerator.Z-numerator.value*rho*2).positive_divide(den*mu,2*logA+logmu))
    if tuple(targets)!=('M','I','S','Cp',divided):raise ValueError('Original relative normalization rows changed')
    repair=core.relative.repair
    caps={name:repair.LogUpper.add(c,[repair.LogUpper(c,None if v.zero else v.record()['log_absolute_upper'])
        for v in (targets[name].value,targets[name].Z)]) for name in repair.ROWS}
    nonzero=[cap for cap in caps.values() if cap.log is not None]
    maximum=repair.LogUpper(c,c.mpf(max(ep(cap.log)[1] for cap in nonzero))) if nonzero else repair.LogUpper(c,None)
    return dict(actual_uniform_N_scaled_target_C0_Z={name:record(targets[name]) for name in repair.ROWS},
        actual_joint_N_scaled_axial_numerator_C0_Z=record(numerator),actual_positive_A_rc_C0_Z=list(amplitude),
        actual_positive_mu=mu,actual_positive_A_rc_log_lower=logA,actual_positive_mu_log=logmu,
        actual_A_rc_Z_over_A_rc=rho,
        target_C1_caps=dict(transformed_N_scaled_target_C1_caps={name:cap.record() for name,cap in caps.items()},
            whole_target_C1_cap=maximum.record(),norm='max_i(sup|N*r_i|+sup|N*r_i_Z|)',
            uniform_for_all_integer_N_ge_N0=True,epsilon_domain='[0,1/N0]'),
        joint_numerator_formed_before_actual_mu_division=True,
        correction_only_complete_minus_identical_background_identity_used=True,
        saved_fixed_N_totals_or_targets_not_rescaled=True,
        interval_hulls_do_not_prove_joint_cancellation=True)


class WholeZAllNOuterRcFunctions(previous.WholeZAllNLongPatchFunctions):
    def __init__(self,dps=500):
        super().__init__(dps)
        receipt=json.loads((HERE/previous.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(previous.GATE) \
                or receipt['source_family']!=self.identity or receipt['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Checked same-input all-N Rh incoming required')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed all-N Rh prerequisite: '+name)
            bind(self.hashes,name,digest)
        bind(self.hashes,previous.RECEIPT,sha(previous.RECEIPT))
        saved=json.loads(gzip.decompress((HERE/previous.NAME).read_bytes()))
        if not saved.get(previous.GATE) or saved['source_family']!=self.identity \
                or saved['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Same-input normalized Rh evaluation required')
        self.Rh_rows={tuple(row['exact_Z_cell']):row for row in saved['actual_four_Z_all_N_long_patch_transports']}
        if set(self.Rh_rows)!=set(CELLS):raise ValueError('Complete four-Z normalized Rh cover required')
        self.views={}
        sources={'o3':self.owner.source,'axial':self.owner.source.source,
            'slope':self.owner.source.source.source,'rh':self.owner.source.source.source.source}
        for kind,source in sources.items():
            self.views[kind]=SimpleNamespace(**{**vars(source),'cache':{},'owner':self.owner.owner})
        self.views['axial'].parent=lambda ends:self.leading_packet(ends,'O2_slope',(1,1))
        self.views['o3'].parent=lambda ends:self.leading_packet(ends,'O2_buffer',(11,1))
        self.views['o3'].query=lambda ends,chart,left,right=None:self.leading_packet(ends,'O3_'+chart,left,right)
        self.outer_binding=dict(original_pure_source_frontends=FRONTEND_BINDINGS,
            pure_parent_chain='axial/buffer <- slope(1); O3 transition <- buffer(11); power <- transition(1)',
            original_global_phase='frac(N*(logRm+original_outer_offset-logRa-hb*s_c/2))',
            outer_offsets=dict(Rh_reference='6+t',O2_slope='6+t',O2_axial='6+exp(40*t)',
                O2_buffer='6+exp40+t',O3_transition='17+exp40+t',O3_power='18+exp40+t'),
            original_physical_measures_and_pressure_rate_zero_retained=True,
            full_phase_is_outer_parameter_cover_not_claimed_physical_image=True,
            relative_target_reference=core.relative.REFERENCE)
        self.logmu=self.owner.source.logmu
        admitted,Tw,margin=o3.source_mu_admission(self.c,self.owner.source.eta_log)
        if admitted._mpi_!=self.logmu._mpi_ or Tw._mpi_!=self.owner.source.Tw._mpi_ \
                or margin._mpi_!=self.owner.source.power_margin._mpi_:
            raise ValueError('Exact same original O3 positive mu/quiet power admission required')
        for module in (outer,rh,slope,axial,o3,rh.original,slope.original,axial.original,o3.original):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def coordinate(self,value,chart=None):
        if chart not in CHARTS:raise ValueError('Original outer chart required')
        q=core.backend.exact(value)
        lo,hi={'Rh_reference':(-5,0),'O2_slope':(0,1),'O2_axial':(0,1),
            'O2_buffer':(0,11),'O3_transition':(0,1),'O3_power':(0,2)}[chart]
        if not lo<=q<=hi:raise ValueError('Original outer native domain required')
        return (q.numerator,q.denominator)

    def leading_packet(self,ends,chart,left,right=None):
        key=tuple(ends);l=self.coordinate(left,chart);r=self.coordinate(left if right is None else right,chart)
        if key not in CELLS or Fraction(*r)<Fraction(*l):raise ValueError('Ordered admitted outer source cell required')
        cachekey=(key,chart,l,r)
        if cachekey in self.packet_cache:return self.packet_cache[cachekey]
        kind='rh' if chart=='Rh_reference' else 'slope' if chart=='O2_slope' else 'axial' if chart.startswith('O2_') else 'o3'
        view=self.views[kind]
        if kind in ('rh','slope'):packet=FRONTENDS[kind](view,key,l,r)
        else:packet=FRONTENDS[kind](view,key,chart.split('_',1)[1],l,r)
        op=self.owner.owner(key)
        if packet['source_identity']!=self.identity or packet['candidate_N']!=self.N0 \
                or packet['exact_common_P0_axial5'] is not op.P0 or not packet['actual_phase_Z_exact_zero']:
            raise ValueError('Same live outer source identity/P0/N0/Z origin required')
        if chart=='O3_power' and (not packet['original_q_C0_Z'][C0].zero or not packet['original_q_C0_Z'][Z].zero):
            raise ValueError('Original admitted O3 quiet power must be exactly flat')
        packet.update(fixed_N_density_support_phase_and_inlet_not_called=True,bridge_collar_not_inherited=True,
            actual_chart=chart,original_power_flat_admission=self.owner.source.power_margin if chart=='O3_power' else None)
        self.packet_cache[cachekey]=packet
        return packet

    def primitive_cell(self,ends,chart,left,right=None):
        # The parent implementation dynamically calls this class's pure
        # leading_packet/coordinate, and uses the generic inverse identity.
        # All outer charts use exactly the same original dstar definition.
        if self.owner.patch_source.dstar_log._mpi_!=self.owner.source.dstar_log._mpi_:
            raise ValueError('Same original inverse dstar required across outer source')
        return super().primitive_cell(ends,chart,left,right)

    def weights(self,ends,chart,left,right,rate):
        f=self.owner.owner(ends).flow
        width,suffix,mass,decay,tail=outer.WholeZFullSignedTransport.outer_weights(self.owner,ends,chart,left,right,rate)
        if ep(width)[0]<=0 or ep(mass.coefficient)[0]<=0:raise ValueError('Positive actual physical outer mass/width required')
        return dict(original_mass=mass,cell_decay=decay,suffix_decay=tail,
            physical_log_width=width,suffix_log_width=suffix)

    def normalized_Rh_incoming(self,ends):
        row=self.Rh_rows[tuple(ends)];op=self.owner.owner(ends);f=op.flow
        if row['source_identity']!=self.identity or row['fixed_leading_input_N0']!=self.N0:
            raise ValueError('Same-input normalized Rh correction required')
        decode=lambda value:core.relative.restore_row(f,value)
        if not previous.equivalent_rows(op.P0,[decode(v) for v in row['exact_common_P0_axial5']]) \
                or not previous.equivalent_rows([op.Rm_factor,op.Rm_factor*self.c.exp(1)],
                    [decode(row['exact_Rm_radius']),decode(row['exact_Rh_radius'])]):
            raise ValueError('Exact same live outer P0/Rm/Rh required')
        incoming={name:Pair(*[decode(v) for v in values]) for name,values in
            row['actual_uniform_N_scaled_Rh_correction_C0_Z'].items()}
        if set(incoming)!=set(RATES) or all(v.value.zero and v.Z.zero for v in incoming.values()):
            raise ValueError('Genuine nonzero correction-only normalized Rh incoming required')
        return incoming

    def transport(self,ends):
        op=self.owner.owner(ends);f=op.flow;incoming=self.normalized_Rh_incoming(ends);initial=dict(incoming);windows=[]
        for chart in CHARTS:
            partition=PARTITIONS[chart];local={name:Pair(f.scalar(0),f.scalar(0)) for name in RATES};cells=[]
            for left,right in zip(partition,partition[1:]):
                query=self.query(ends,chart,left,right);weights={};contributions={}
                for name,rate in RATES.items():
                    weight=self.weights(ends,chart,left,right,rate)
                    contribution=scale(query['N_scaled_density_C0_Z'][name],weight['original_mass']*weight['suffix_decay'])
                    local[name]=add(local[name],contribution);weights[name]=weight;contributions[name]=contribution
                cells.append(dict(exact_left=left,exact_right=right,epsilon=query['epsilon'],
                    full_phase_inverse_function_alternatives=query['source']['full_phase_inverse_function_alternatives'],
                    original_E_V_C0_Z={k:record(v) for k,v in query['source']['original_E_V_C0_Z'].items()},
                    N_scaled_density_C0_Z={k:record(v) for k,v in query['N_scaled_density_C0_Z'].items()},
                    actual_N_scaled_contribution_C0_Z={k:record(v) for k,v in contributions.items()},
                    original_own_rate_weights=weights,physical_Jacobian_applied_once=True,
                    all_N_phase_cover_not_fixed_N_phase_reuse=True,bridge_collar_not_inherited=True))
                print('Actual all-N outer '+chart+' cell: '+str(ends)+' '+str(left)+' '+str(right),flush=True)
            memory={name:self.weights(ends,chart,partition[0],partition[-1],rate)['cell_decay'] for name,rate in RATES.items()}
            outgoing={name:add(scale(incoming[name],memory[name]),local[name]) for name in RATES}
            windows.append(dict(actual_chart=chart,actual_all_N_source_cells=cells,
                actual_N_scaled_incoming_C0_Z={k:record(v) for k,v in incoming.items()},
                original_incoming_own_rate_memory=memory,actual_N_scaled_local_C0_Z={k:record(v) for k,v in local.items()},
                actual_N_scaled_outgoing_C0_Z={k:record(v) for k,v in outgoing.items()}))
            incoming=outgoing
        terminal=self.leading_packet(ends,'O3_power',PARTITIONS['O3_power'][-1])
        amplitude=[terminal['original_roots']['E'][part] for part in (C0,Z)]
        mu=terminal['original_closed_O3_background']['original_positive_formal_mu']
        if not previous.equivalent_rows([mu],[f.factor((0,0,0,0,0),self.logmu)]):
            raise ValueError('Terminal source mu differs from exact admitted O3 parameter')
        # Keep the actual O3 offset representation. The equivalent
        # Pstar*exp(9) radius has different formal powers; raw tuple
        # equality is not a valid test between those representations.
        if self.owner.source.parameters['exact_logPstar']!='exp(40)+11':
            raise ValueError('Same original Rc/Pstar defining parameter identity required')
        radius=op.Rm_factor*f.factor((0,0,0,0,0),self.c.exp(40)+18+self.c.mpf(2))
        actual_radius=terminal['original_closed_O3_background']['actual_source_radius']
        if not core.relative.exact_row(radius,actual_radius):raise ValueError('Original Rc radius source identity changed')
        target=normalized_targets(f,incoming,amplitude,mu,self.logmu)
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),fixed_leading_input_N0=self.N0,
            epsilon=core.epsilon_cover(self.c,self.N0),exact_common_P0_axial5=op.P0,
            exact_Rm_radius=op.Rm_factor,exact_Rh_radius=op.Rm_factor*self.c.exp(1),exact_Rc_radius=radius,
            exact_Rc_radius_definition='logRc=logRm+exp(40)+20=logRm+logPstar+9',
            actual_uniform_N_scaled_Rh_incoming_C0_Z={k:record(v) for k,v in initial.items()},
            actual_all_N_outer_windows=windows,actual_uniform_N_scaled_Rc_correction_C0_Z={k:record(v) for k,v in incoming.items()},
            actual_original_Rc_background_C0_Z={name:rows[:2] for name,rows in
                terminal['original_generic_source']['common_own_five_histories_axial5'].items()},
            actual_positive_Rc_normalized_amplitude_C0_Z=amplitude,exact_source_log_mu=self.logmu,
            original_quiet_power_flat_log_margin=self.owner.source.power_margin,**target,
            normalized_correction_only_incoming_background_P0_distinct=True,
            actual_uniform_N_scaled_Rc_targets_installed=True,actual_global_frequency_admitted=False,
            **dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZAllNOuterRcFunctions();rows=[owner.transport(ends) for ends in CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,fixed_leading_input_N0=owner.N0,
            exact_Z_partition=CELLS,current_outer_charts=CHARTS,exact_outer_native_partitions=PARTITIONS,
            original_outer_source_binding=owner.outer_binding,actual_four_Z_all_N_outer_Rc_transports=rows,
            exact_quotient_Z_theorem=core.relative.exact_quotient_theorem(),
            actual_full_phase_all_N_outer_driver_functions_installed=True,
            actual_uniform_N_scaled_Rc_targets_installed=True,actual_all_regions_all_N_adapter_installed=True,
            actual_N_squared_averaging_cancellation_proved=False,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Same leading input and full-phase epsilon family of outer C1 functions carries checked normalized '
                'Rh correction through six original physical windows to Rc and directly constructs the five '
                'N-scaled relative targets with joint axial normalization. No fixed-N totals/targets are rescaled. '
                'Common N selection, contraction/control solution, N^-2 averaging, input exterior/heat, stress, '
                'higher jets, genuine temporal recursion, pulses and full corrected NS remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(core.current.encode(serialized(report)),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()

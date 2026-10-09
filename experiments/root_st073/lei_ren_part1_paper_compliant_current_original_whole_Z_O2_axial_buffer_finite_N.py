"""Whole-Z original O2 axial/buffer functions and genuine same-N Rd prefix."""
import ast
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_O2_slope_finite_N as upstream
import lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N as original

reference,switch=original.reference,original.switch
HERE,PREFIX,sha,bind,ep=upstream.HERE,upstream.PREFIX,upstream.sha,upstream.bind,upstream.ep
RATES,C0,Z,OPEN,CELLS=upstream.RATES,upstream.C0,upstream.Z,upstream.OPEN,upstream.CELLS
PARTITIONS,MD,KERNEL_CELLS,KERNEL_WINDOW=original.PARTITIONS,original.MD,original.KERNEL_CELLS,original.KERNEL_WINDOW
NAME=PREFIX+'current_original_whole_Z_O2_axial_buffer_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_O2_axial_buffer_finite_N_check.json'
GATE='current_original_whole_Z_same_N_actual_O2_axial_buffer_source_and_correction_through_Rd_enclosed'
serialized,encode=upstream.serialized,upstream.encode


def compiled_background(callback):
    """Original background AST unchanged; cache only its pure scalar kernels."""
    tree=ast.parse(Path(original.__file__).read_text(encoding='utf8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='background_cell')
    calls=[n.func.attr for n in ast.walk(fn) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
        and isinstance(n.func.value,ast.Name) and n.func.value.id=='kernels']
    if calls!=['turnoff_kernels','turnoff_kernels']:
        raise ValueError('Original turnoff kernel call sites changed')
    binding=dict(original_background_AST_sha256=hashlib.sha256(ast.dump(fn).encode()).hexdigest(),
        unchanged_original_source_arithmetic=True,only_original_pure_scalar_repeated_evaluation_cached=True,
        cache_key_binds_context_precision_exact_y_interval_Md_partition_window=True)
    namespace={**vars(original),'kernels':SimpleNamespace(turnoff_kernels=callback)}
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<unchanged original axial source with kernel memoization>','exec'),namespace)
    return namespace['background_cell'],binding


def exact_row_equal(live,old):
    return live.coefficient._mpi_==old.coefficient._mpi_ and live.scale.powers==old.scale.powers \
        and live.scale.offset._mpi_==old.scale.offset._mpi_


class WholeZO2AxialBufferFiniteN:
    def __init__(self,dps=500):
        self.source=upstream.WholeZO2SlopeFiniteN(dps)
        self.c=self.source.c;self.identity=self.source.identity;self.N=self.source.N
        self.hashes=dict(self.source.hashes)
        receipt=json.loads((HERE/upstream.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(upstream.GATE) or receipt['source_family']!=self.identity \
                or receipt['candidate_N']!=self.N or not receipt.get('genuine_Rref_correction_only_incoming_true_memory_and_slope_exit_checked'):
            raise ValueError('Checked same-source whole-Z slope-exit correction required')
        for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,upstream.RECEIPT,sha(upstream.RECEIPT))
        self.saved=json.loads(gzip.decompress((HERE/upstream.NAME).read_bytes()))
        if not self.saved[upstream.GATE] or self.saved['candidate_N']!=self.N or self.saved['source_family']!=self.identity:
            raise ValueError('Same current whole-Z slope-exit packet and unchanged N required')
        self.saved_cells={tuple(row['exact_Z_cell']):row for row in self.saved['source_cells']}
        if set(self.saved_cells)!=set(CELLS):raise ValueError('Complete whole-Z slope-exit correction cover required')
        self.parameters=original.parameter_binding(self.identity['actual_five_defect_family_sha256'],self.hashes)
        if self.parameters['exact_parameter_source_family']!=self.identity:
            raise ValueError('Same whole-Z implicit definition and pressure datum family required')
        self.eta_log=self.source.eta_log;self.dstar_log=self.source.dstar_log;self.sc=self.source.sc
        self.kernel_cache={};self.kernel_evaluations=0;self.kernel_hits=0
        self.background,self.kernel_binding=compiled_background(self.kernels)
        self.bindings=dict(original_source=original.source_bindings(),kernel_backend=self.kernel_binding,
            actual_global_phase='frac(N*(logRm+6+y-logRa-hb*s_c/2))',
            actual_current_source_identity_not_N257_label_plumbing=True)
        self.cache={};self.inlets={};self.parents={};self.joins={}
        for module in (upstream,original,original.kernels,original.previous.original,
                original.parameter_source,reference,reference.long,switch,reference.parameters,
                reference.primitives,reference.phase,reference.bounds):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def kernels(self,c,y,md,cells,window):
        if c is not self.c:raise ValueError('Same live scalar interval context required')
        key=(c.dps,y._mpi_,md,cells,window)
        if key not in self.kernel_cache:
            self.kernel_cache[key]=original.kernels.turnoff_kernels(c,y,md,cells,window);self.kernel_evaluations+=1
        else:self.kernel_hits+=1
        return self.kernel_cache[key]

    def owner(self,ends):return self.source.owner(ends)

    def inlet(self,ends):
        key=tuple(ends)
        if key in self.inlets:return self.inlets[key]
        op=self.owner(ends);f,c=op.flow,op.c;stored=self.saved_cells[key]
        if stored['source_identity']!=self.identity or stored['candidate_N']!=self.N \
                or not stored['full_same_N_correction_prefix_through_O2_slope_exit_installed']:
            raise ValueError('Genuine current same-N slope-exit source required')
        decode=lambda row:upstream.upstream.upstream.prefix.upstream.decode(f,row)
        P0=decode(stored['exact_common_P0_axial5'])
        if len(P0)!=6 or any(not exact_row_equal(live,old) for live,old in zip(op.P0,P0)):
            raise ValueError('Same independent original P0 source tuples required')
        if not exact_row_equal(op.Rm_factor,decode(stored['exact_Rm_radius'])) \
                or not exact_row_equal(op.Rm_factor*c.exp(7),decode(stored['exact_O2_slope_exit_radius'])):
            raise ValueError('Same exact current Rm/slope-exit radii required')
        incoming=decode(stored['actual_O2_slope_exit_correction_C0_Z'])
        if set(incoming)!=set(RATES) or any(len(rows)!=2 for rows in incoming.values()):
            raise ValueError('Actual slope-exit correction-only five C0/Z rows required')
        reference.parameters.same_source(f,[row for rows in incoming.values() for row in rows])
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.P0,exact_Rm_radius=op.Rm_factor,exact_O2_slope_exit_radius=op.Rm_factor*c.exp(7),
            actual_O2_slope_exit_correction_C0_Z=incoming,incoming_source_report=upstream.NAME,
            correction_only_not_complete_own_history=True,no_old_N257_label_incoming_or_phase_receipt_transplanted=True)
        self.inlets[key]=result;return result

    def parent(self,ends):
        key=tuple(ends)
        if key not in self.parents:self.parents[key]=self.source.query(ends,(1,1))
        return self.parents[key]

    def query(self,ends,chart,left,right=None):
        right=left if right is None else right;key=(tuple(ends),chart,left,right)
        original.coordinate(self.c,chart,left,right)
        if key in self.cache:return self.cache[key]
        op=self.owner(ends);f,c=op.flow,op.c
        source=self.background(op,self.parent(ends),chart,left,right);a=source['actual_a_axial5']
        E=source['raw']['raw_current_radius_y_derivative_axial_coefficients']['velocity']['theta'][0]
        proxy=SimpleNamespace(flow=f,c=c,reference=op.reference,zrows=op.zrows,P0=op.P0,
            Pstar=op.Pstar,source_radius=source['actual_source_radius'],correlated_C=f.multiply(a,E))
        generic=reference.long.RECOVER(proxy,source['raw'])
        generic.pop('source_frame_conditional_on_same_actual_Rm_inlet')
        generic.update(source_from_same_current_whole_Z_original_O2_axial_buffer_background=True,
            leading_background_not_finite_N_correction=True)
        roots,qr,proof=original.CRITICAL_QUOTIENTS(proxy,generic,a,self.eta_log)
        got=reference.primitives.all_u_primitive_bounds(f,roots,qr,self.dstar_log,c.mpf([0,1]))
        if qr[C0].zero and qr[Z].zero:
            got['values']={name:f.scalar(0) for name in got['values']}
            got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
        A=got['values']['A'];logA=None if A.zero else ep(A.record()['log_absolute_upper'])[1]
        if logA is not None and logA>=ep(c.ln(self.N))[0]:
            raise ValueError('Same N fails actual O2 axial/buffer A/N budget at '+str(key))
        got=upstream.upstream.upstream.prefix.upstream.upstream.bounded_exponent_cover(f,got,self.N)
        V=generic['common_velocity_V_axial5']
        density=reference.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],self.N)
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_left=list(left),exact_right=list(right),exact_common_P0_axial5=op.P0,actual_chart=chart,
            original_closed_O2_axial_buffer_background=source,original_generic_source=generic,
            original_full_source_quotients=proof,original_roots=roots['roots'],original_q_C0_Z=qr,
            original_primitive_values=got['values'],original_primitive_proof=got['record'],
            actual_original_A_log_absolute_upper=None if logA is None else c.mpf(logA),actual_A_below_same_candidate_N=True,
            original_signed_five_density_C0_Z=density,
            actual_global_radius_phase='frac(N*(logRm+6+y-logRa-hb*s_c/2))',
            actual_phase_origin_s_c=self.sc,actual_phase_Z_exact_zero=True,
            actual_radius_phase_full_period_cover=c.mpf([0,1]),phase_not_restarted_at_O2_slope_exit=True,
            phase_average_cancellation_not_claimed=True,
            actual_physical_log_measure='dy=40*exp(40*phase)*dphase' if chart=='axial' else 'dy=dtau',
            signed_axis_no_absZ_or_positive_rho_division=True,
            exact_critical_a2_and_collected_Delta_b_squared_over2=True,
            actual_high_order_finite_N_correction_jets_installed=False)
        self.cache[key]=result;return result

    def seam(self,ends,left,right):
        key=(tuple(ends),left,right)
        if key in self.joins:return self.joins[key]
        op=self.owner(ends)
        before=self.parent(ends) if left=='slope_exit' else self.query(ends,'axial',(1,1))
        after=self.query(ends,'axial',(0,1)) if right=='axial_inlet' else self.query(ends,'buffer',(0,1))
        a,b=before['original_generic_source'],after['original_generic_source'];count=0
        pairs=[(a['common_velocity_E_axial5'],b['common_velocity_E_axial5']),
            (a['common_velocity_V_axial5'],b['common_velocity_V_axial5'])]
        pairs.extend((a['common_own_five_histories_axial5'][name],b['common_own_five_histories_axial5'][name]) for name in RATES)
        for old,new in pairs:
            for x,y in zip(old,new):
                difference=x-y
                if not difference.zero and not ep(difference.coefficient)[0]<=0<=ep(difference.coefficient)[1]:
                    raise ValueError('Defining whole-Z axial/buffer source covers do not meet')
                count+=1
        result=dict(exact_original_function_join_by_closed_kernel_and_decay_identities=True,
            directed_source_overlap_consistency_rows=count,overlap_only_consistency_not_functional_identity_proof=True,
            same_original_P0_and_slope_history_seed=a['common_original_P0_axial5'] is b['common_original_P0_axial5'] is op.P0,
            exact_radius_and_global_phase_join=True)
        self.joins[key]=result;return result

    def contribution(self,ends):
        op=self.owner(ends);f,c=op.flow,op.c;inlet=self.inlet(ends)
        initial=inlet['actual_O2_slope_exit_correction_C0_Z'];current=initial;charts={}
        joins=dict(slope_axial=self.seam(ends,'slope_exit','axial_inlet'))
        for chart,partition in PARTITIONS.items():
            if chart=='buffer':joins['axial_buffer']=self.seam(ends,'axial_exit','buffer_inlet')
            local={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
            for left,right in zip(partition,partition[1:]):
                source=self.query(ends,chart,left,right);density=source['original_signed_five_density_C0_Z'];drivers={};weights={}
                for name,rate in RATES.items():
                    width,suffix,mass,decay=original.physical_weights(f,chart,left,right,rate)
                    reference.parameters.positive_source(f,f.scalar(mass),'actual_O2_'+chart+'_Duhamel_mass')
                    pair=[reference.bounds.symmetric(f,reference.bounds.magnitude(f,density[part][name])*mass*decay)
                        for part in ('kernels','Z_derivatives')]
                    for n,row in enumerate(pair):local[name][n]+=row
                    drivers[name]=pair;weights[name]=dict(true_full_mass=mass,true_suffix_decay=decay,
                        own_rate=str(rate),true_radial_measure_applied_once=True)
                cells.append(dict(source=source,actual_log_width=width,actual_suffix_to_chart_exit=suffix,
                    signed_cell_driver_C0_Z=drivers,own_rate_weights=weights))
                print('Whole-Z actual same-N O2 '+chart+' driver: '+str(ends)+' '+str(left)+' '+str(right),flush=True)
            width=c.expm1(MD) if chart=='axial' else c.mpf(11)
            memory={name:f.scalar(1) if not rate else f.factor((0,0,0,0,0),-width*c.mpf(rate.numerator)/rate.denominator)
                for name,rate in RATES.items()}
            inherited={name:[row*memory[name] for row in current[name]] for name in RATES}
            output=original.formal_affine_transport(f,current,local,width)
            terminal=self.query(ends,chart,partition[-1])
            background={name:list(rows[:2]) for name,rows in terminal['original_generic_source']['common_own_five_histories_axial5'].items()}
            complete={name:[background[name][n]+output[name][n] for n in range(2)] for name in RATES}
            charts[chart]=dict(actual_chart_incoming_C0_Z=current,incoming_own_rate_memory=memory,
                actual_retained_incoming_C0_Z=inherited,actual_chart_local_driver_C0_Z=local,
                actual_chart_exit_correction_C0_Z=output,actual_original_chart_exit_background_C0_Z=background,
                actual_chart_exit_complete_own_history_C0_Z=complete,actual_source_cells=cells,
                exact_source_partition=partition,exact_total_log_width=width,pressure_memory_exactly_one=True,
                actual_original_chart_exit_source=terminal)
            current=output
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.P0,exact_Rm_radius=op.Rm_factor,exact_O2_slope_exit_radius=op.Rm_factor*c.exp(7),
            exact_O2_axial_exit_radius=op.Rm_factor*f.factor((0,0,0,0,0),6+c.exp(MD)),
            exact_Rd_radius=op.Rm_factor*f.factor((0,0,0,0,0),17+c.exp(MD)),
            actual_slope_exit_incoming_binding=inlet,actual_O2_slope_exit_correction_incoming_C0_Z=initial,
            actual_O2_axial_buffer_charts=charts,typed_actual_O2_axial_buffer_source_joins=joins,
            actual_Rd_correction_C0_Z=current,
            actual_original_Rd_background_C0_Z=charts['buffer']['actual_original_chart_exit_background_C0_Z'],
            actual_Rd_complete_own_history_C0_Z=charts['buffer']['actual_chart_exit_complete_own_history_C0_Z'],
            actual_original_Rd_source=charts['buffer']['actual_original_chart_exit_source'],
            full_same_N_correction_prefix_through_Rd_installed=True,
            original_turnoff_kernels_nonzero_far_tail_and_buffer_history_retained=True,
            original_P0_background_and_finite_N_correction_kept_distinct=True,**dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZO2AxialBufferFiniteN();cells=[serialized(owner.contribution(ends)) for ends in CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_domain=['-1','1'],exact_Z_partition=[list(ends) for ends in CELLS],
            exact_axial_selector_domain=['0','1'],exact_buffer_tau_domain=['0','11'],
            exact_chart_partitions=PARTITIONS,source_cells=cells,original_source_bindings=owner.bindings,
            exact_original_parameter_binding=owner.parameters,original_eta_log=owner.eta_log,original_dstar_log=owner.dstar_log,
            original_kernel_partition=KERNEL_CELLS,original_kernel_window=KERNEL_WINDOW,
            original_kernel_cache_evaluations=owner.kernel_evaluations,original_kernel_cache_hits=owner.kernel_hits,
            full_same_N_correction_prefix_through_Rd_installed=True,
            actual_high_order_finite_N_correction_jets_installed=False,**dict.fromkeys(OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual whole-Z original O2 axial turnoff and11-unit buffer functions, critical a2 '
                  'Delta=b^2/2 and genuine same-N slope-exit correction through Rd. O3/Rc transport, '
                  'high correction jets/repair, global cone/heat/N and n-recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z genuine same-N correction prefix reaches Rd',flush=True);return report


if __name__=='__main__':run()

"""Whole-Z original O2 slope functions and genuine same-N correction prefix."""
import ast
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_reference_finite_N as upstream
import lei_ren_part1_paper_compliant_current_original_O2_slope_finite_N as original

reference,switch=original.reference,original.switch
HERE,PREFIX,sha,read,bind,ep=upstream.HERE,upstream.PREFIX,upstream.sha,upstream.read,upstream.bind,upstream.ep
RATES,C0,Z,OPEN=upstream.RATES,upstream.C0,upstream.Z,upstream.OPEN
CELLS=upstream.upstream.upstream.upstream.source.CELLS
PARTITION=original.PARTITION
SCALAR_CELLS=original.SCALAR_CELLS
NAME=PREFIX+'current_original_whole_Z_O2_slope_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_O2_slope_finite_N_check.json'
GATE='current_original_whole_Z_same_N_actual_O2_slope_source_and_correction_through_slope_exit_enclosed'


def serialized(value):return upstream.serialized(value)


def encode(value):return upstream.upstream.upstream.upstream.source.bridge._encode(value)


def compiled_background(callback):
    """Unchanged original source AST, with memoized original scalar helper."""
    tree=ast.parse(Path(original.__file__).read_text(encoding='utf8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='background_cell')
    calls=[n.func.attr for n in ast.walk(fn) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
        and isinstance(n.func.value,ast.Name) and n.func.value.id=='original']
    if sorted(calls)!=['sigma_jets','slope_masses']:
        raise ValueError('Original slope scalar helper call sites changed')
    binding=dict(original_background_AST_sha256=hashlib.sha256(ast.dump(fn).encode()).hexdigest(),
        unchanged_original_source_arithmetic=True,only_original_pure_scalar_repeated_evaluation_cached=True,
        cache_key_binds_context_precision_exact_y_interval_and_scalar_partition=True)
    namespace={**vars(original),'original':SimpleNamespace(slope_masses=callback,
        sigma_jets=original.original.sigma_jets)}
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<unchanged original slope with scalar memoization>','exec'),namespace)
    return namespace['background_cell'],binding


def slope_weights(f,left,right,rate):
    c=f.c;l,r=original.fraction(left),original.fraction(right)
    if l>=r:raise ValueError('Strictly ordered original slope transport cell required')
    width=c.mpf((r-l).numerator)/(r-l).denominator;suffix=1-c.mpf(r.numerator)/r.denominator
    mass,decay,tail=reference.long.kernel_weight(f,width,suffix,rate)
    return width,suffix,mass,decay,tail


class WholeZO2SlopeFiniteN:
    def __init__(self,dps=500):
        self.source=upstream.WholeZRhReferenceFiniteN(dps)
        self.c=self.source.c;self.identity=self.source.identity;self.N=self.source.N
        self.hashes=dict(self.source.hashes)
        receipt=json.loads((HERE/upstream.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(upstream.GATE) or receipt['source_family']!=self.identity \
                or receipt['candidate_N']!=self.N or not receipt.get('genuine_Rh_correction_only_incoming_true_memory_and_Rref_exit_checked'):
            raise ValueError('Checked same-source whole-Z Rref correction incoming required')
        for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,upstream.RECEIPT,sha(upstream.RECEIPT))
        self.saved=json.loads(gzip.decompress((HERE/upstream.NAME).read_bytes()))
        if not self.saved[upstream.GATE] or self.saved['candidate_N']!=self.N or self.saved['source_family']!=self.identity:
            raise ValueError('Same current whole-Z Rref packet and unchanged N required')
        self.saved_cells={tuple(row['exact_Z_cell']):row for row in self.saved['source_cells']}
        if set(self.saved_cells)!=set(CELLS):raise ValueError('Complete original whole-Z Rref correction cover required')
        self.eta_log=self.source.eta_log;self.dstar_log=self.source.dstar_log;self.sc=self.source.sc
        self.scalar_cache={};self.scalar_evaluations=0;self.scalar_hits=0
        self.background,self.scalar_binding=compiled_background(self.scalars)
        self.bindings=dict(original_source=original.source_bindings(),scalar_backend=self.scalar_binding,
            actual_global_phase='frac(N*(logRm+6+y-logRa-hb*s_c/2))',
            actual_current_source_identity_not_N257_label_plumbing=True)
        self.cache={};self.inlets={};self.joins={}
        for module in (upstream,original,original.original,reference,reference.long,switch,
                reference.parameters,reference.primitives,reference.phase,reference.bounds):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,'lei_ren_part1_paper_interval_outer_slope_field.py',sha('lei_ren_part1_paper_interval_outer_slope_field.py'))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def scalars(self,c,y,cells):
        if c is not self.c:raise ValueError('Same live scalar interval context required')
        key=(c.dps,y._mpi_,cells)
        if key not in self.scalar_cache:
            self.scalar_cache[key]=original.original.slope_masses(c,y,cells);self.scalar_evaluations+=1
        else:self.scalar_hits+=1
        return self.scalar_cache[key]

    def owner(self,ends):return self.source.owner(ends)

    def inlet(self,ends):
        key=tuple(ends)
        if key in self.inlets:return self.inlets[key]
        op=self.owner(ends);f,c=op.flow,op.c;stored=self.saved_cells[key]
        if stored['source_identity']!=self.identity or stored['candidate_N']!=self.N \
                or not stored['full_same_N_correction_prefix_through_Rref_installed']:
            raise ValueError('Genuine current same-N Rref source required')
        decode=lambda row:upstream.upstream.prefix.upstream.decode(f,row)
        P0=decode(stored['exact_common_P0_axial5'])
        def equal(live,old):
            return live.coefficient._mpi_==old.coefficient._mpi_ and live.scale.powers==old.scale.powers \
                and live.scale.offset._mpi_==old.scale.offset._mpi_
        if len(P0)!=6 or any(not equal(live,old) for live,old in zip(op.P0,P0)):
            raise ValueError('Same independent original P0 source tuples required')
        if not equal(op.Rm_factor,decode(stored['exact_Rm_radius'])) \
                or not equal(op.Rm_factor*c.exp(6),decode(stored['exact_Rref_radius'])):
            raise ValueError('Same exact current Rm/Rref radii required')
        incoming=decode(stored['actual_Rref_correction_C0_Z'])
        if set(incoming)!=set(RATES) or any(len(rows)!=2 for rows in incoming.values()):
            raise ValueError('Actual Rref correction-only five C0/Z rows required')
        reference.parameters.same_source(f,[row for rows in incoming.values() for row in rows])
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.P0,exact_Rm_radius=op.Rm_factor,exact_Rref_radius=op.Rm_factor*c.exp(6),
            actual_Rref_correction_C0_Z=incoming,incoming_source_report=upstream.NAME,
            correction_only_not_complete_own_history=True,no_old_N257_label_incoming_or_phase_receipt_transplanted=True)
        self.inlets[key]=result;return result

    def query(self,ends,left,right=None):
        right=left if right is None else right;key=(tuple(ends),left,right)
        if key in self.cache:return self.cache[key]
        op=self.owner(ends);f,c=op.flow,op.c
        source=self.background(op,left,right,SCALAR_CELLS);a=source['actual_a_axial5']
        E=source['raw']['raw_current_radius_y_derivative_axial_coefficients']['velocity']['theta'][0]
        proxy=SimpleNamespace(flow=f,c=c,reference=op.reference,zrows=op.zrows,P0=op.P0,
            Pstar=op.Pstar,source_radius=source['actual_source_radius'],correlated_C=f.multiply(a,E))
        generic=reference.long.RECOVER(proxy,source['raw'])
        generic.pop('source_frame_conditional_on_same_actual_Rm_inlet')
        generic.update(source_from_same_current_whole_Z_original_O2_slope_background=True,
            leading_background_not_finite_N_correction=True)
        roots,qr,proof=switch.general_quotients(proxy,generic,a,self.eta_log)
        got=reference.primitives.all_u_primitive_bounds(f,roots,qr,self.dstar_log,c.mpf([0,1]))
        if qr[C0].zero and qr[Z].zero:
            got['values']={name:f.scalar(0) for name in got['values']}
            got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
        A=got['values']['A'];logA=None if A.zero else ep(A.record()['log_absolute_upper'])[1]
        if logA is not None and logA>=ep(c.ln(self.N))[0]:
            raise ValueError('Same N fails actual O2 slope A/N budget at '+str(key))
        got=upstream.upstream.prefix.upstream.upstream.bounded_exponent_cover(f,got,self.N)
        V=generic['common_velocity_V_axial5']
        density=reference.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],self.N)
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_left=list(left),exact_right=list(right),exact_common_P0_axial5=op.P0,
            original_closed_O2_slope_background=source,original_generic_source=generic,
            original_full_source_quotients=proof,original_roots=roots['roots'],original_q_C0_Z=qr,
            original_primitive_values=got['values'],original_primitive_proof=got['record'],
            actual_original_A_log_absolute_upper=None if logA is None else c.mpf(logA),actual_A_below_same_candidate_N=True,
            original_signed_five_density_C0_Z=density,
            actual_global_radius_phase='frac(N*(logRm+6+y-logRa-hb*s_c/2))',
            actual_phase_origin_s_c=self.sc,actual_phase_Z_exact_zero=True,
            actual_radius_phase_full_period_cover=c.mpf([0,1]),phase_not_restarted_at_Rref=True,
            phase_average_cancellation_not_claimed=True,actual_physical_log_measure='dy=d(y)',
            signed_axis_no_absZ_or_positive_rho_division=True,
            actual_high_order_finite_N_correction_jets_installed=False)
        self.cache[key]=result;return result

    def seam(self,ends):
        key=tuple(ends)
        if key in self.joins:return self.joins[key]
        op=self.owner(ends)
        background=upstream.original.background_cell(op,(0,1),(0,1))
        _,_,ref=reference.recover_cell(op.reference,background)
        slope=self.query(ends,(0,1))['original_generic_source'];count=0
        pairs=[(ref['common_velocity_E_axial5'],slope['common_velocity_E_axial5']),
            (ref['common_velocity_V_axial5'],slope['common_velocity_V_axial5'])]
        pairs.extend((ref['common_own_five_histories_axial5'][name],slope['common_own_five_histories_axial5'][name]) for name in RATES)
        for before,after in pairs:
            for a,b in zip(before,after):
                difference=a-b
                if not difference.zero and not ep(difference.coefficient)[0]<=0<=ep(difference.coefficient)[1]:
                    raise ValueError('Defining whole-Z Rref/slope source covers do not meet')
                count+=1
        result=dict(exact_original_function_join_by_J_and_masses_zero_at_y0=True,
            same_original_Rref_P0_and_five_history_seed=ref['common_original_P0_axial5'] is op.P0,
            exact_radius_and_phase_join='Rref*exp(0)=Rref; slope y=0 equals reference offset0',
            directed_source_overlap_consistency_rows=count,overlap_only_consistency_not_functional_identity_proof=True)
        self.joins[key]=result;return result

    def contribution(self,ends):
        op=self.owner(ends);f,c=op.flow,op.c;inlet=self.inlet(ends);join=self.seam(ends)
        incoming=inlet['actual_Rref_correction_C0_Z'];total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
        for left,right in zip(PARTITION,PARTITION[1:]):
            source=self.query(ends,left,right);density=source['original_signed_five_density_C0_Z'];rows={};weights={}
            for name,rate in RATES.items():
                width,suffix,mass,decay,tail=slope_weights(f,left,right,rate)
                reference.parameters.positive_source(f,mass,'actual_O2_slope_Duhamel_mass')
                pair=[reference.bounds.symmetric(f,reference.bounds.magnitude(f,density[part][name])*mass*tail)
                    for part in ('kernels','Z_derivatives')]
                for n,row in enumerate(pair):total[name][n]+=row
                rows[name]=pair;weights[name]=dict(true_full_mass=mass,true_cell_decay=decay,
                    true_suffix_decay=tail,own_rate=str(rate),true_radial_measure_applied_once=True)
            cells.append(dict(source=source,actual_log_width=width,actual_suffix_to_slope_exit=suffix,
                signed_cell_driver_C0_Z=rows,own_rate_weights=weights))
            print('Whole-Z actual same-N O2 slope driver: '+str(ends)+' '+str(left)+' '+str(right),flush=True)
        memory={name:reference.long.kernel_weight(f,c.mpf(1),c.mpf(0),rate)[1] for name,rate in RATES.items()}
        inherited={name:[row*memory[name] for row in incoming[name]] for name in RATES}
        outgoing={name:[inherited[name][n]+total[name][n] for n in range(2)] for name in RATES}
        terminal=self.query(ends,(1,1))
        background={name:list(rows[:2]) for name,rows in terminal['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete={name:[background[name][n]+outgoing[name][n] for n in range(2)] for name in RATES}
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.P0,exact_Rm_radius=op.Rm_factor,exact_Rref_radius=op.Rm_factor*c.exp(6),
            exact_O2_slope_exit_radius=op.Rm_factor*c.exp(7),actual_Rref_incoming_binding=inlet,
            exact_total_log_length=1,actual_source_cells=cells,typed_actual_Rref_slope_source_join=join,
            actual_Rref_correction_incoming_C0_Z=incoming,incoming_own_rate_memory=memory,
            actual_local_driver_C0_Z=total,actual_retained_incoming_C0_Z=inherited,
            actual_O2_slope_exit_correction_C0_Z=outgoing,actual_original_O2_slope_exit_background_C0_Z=background,
            actual_O2_slope_exit_complete_own_history_C0_Z=complete,actual_original_O2_slope_exit_source=terminal,
            full_same_N_correction_prefix_through_O2_slope_exit_installed=True,
            original_P0_background_and_finite_N_correction_kept_distinct=True,**dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZO2SlopeFiniteN();cells=[serialized(owner.contribution(ends)) for ends in CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_domain=['-1','1'],exact_Z_partition=[list(ends) for ends in CELLS],
            exact_slope_y_domain=['0','1'],exact_slope_partition=PARTITION,source_cells=cells,
            original_source_bindings=owner.bindings,original_eta_log=owner.eta_log,original_dstar_log=owner.dstar_log,
            original_scalar_partition=SCALAR_CELLS,original_scalar_cache_evaluations=owner.scalar_evaluations,
            original_scalar_cache_hits=owner.scalar_hits,full_same_N_correction_prefix_through_O2_slope_exit_installed=True,
            actual_high_order_finite_N_correction_jets_installed=False,**dict.fromkeys(OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual whole-Z original O2 slope scalar history functions, variable shear with exact '
                  'critical endpoint and genuine same-N Rref correction through slope exit. Axial/buffer '
                  'and O3/Rc transport, high correction jets/repair, global cone/heat/N and n-recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z genuine same-N correction prefix reaches O2 slope exit',flush=True);return report


if __name__=='__main__':run()

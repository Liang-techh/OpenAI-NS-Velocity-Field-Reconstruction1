"""Original whole-Z Rh reference functions and genuine same-N prefix to Rref."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_finite_N as upstream
import lei_ren_part1_paper_compliant_current_original_Rh_reference_finite_N as original

reference=original.reference
HERE,PREFIX,sha,read,bind,ep=upstream.HERE,upstream.PREFIX,upstream.sha,upstream.read,upstream.bind,upstream.ep
RATES,C0,Z,OPEN=upstream.RATES,upstream.C0,upstream.Z,upstream.OPEN
PARTITION=original.PARTITION
NAME=PREFIX+'current_original_whole_Z_Rh_reference_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_Rh_reference_finite_N_check.json'
GATE='current_original_whole_Z_same_N_actual_Rh_reference_source_and_correction_through_Rref_enclosed'


def serialized(value):return upstream.serialized(value)


def reference_weights(f,left,right,rate):
    c=f.c;l,r=original.fraction(left),original.fraction(right)
    if l>=r:raise ValueError('Strictly ordered original reference transport cell required')
    width=c.mpf((r-l).numerator)/(r-l).denominator;suffix=-c.mpf(r.numerator)/r.denominator
    mass,decay,tail=upstream.long.kernel_weight(f,width,suffix,rate)
    return width,suffix,mass,decay,tail


class WholeZRhReferenceFiniteN:
    def __init__(self,dps=500):
        self.source=upstream.WholeZRmPatchFiniteN(dps)
        self.c=self.source.c;self.identity=self.source.identity;self.N=self.source.N
        self.hashes=dict(self.source.hashes)
        receipt=json.loads((HERE/upstream.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(upstream.GATE) or receipt['source_family']!=self.identity \
                or receipt['candidate_N']!=self.N or not receipt.get('genuine_Rm_correction_only_incoming_true_own_rate_memory_and_Rh_exit_checked'):
            raise ValueError('Checked genuine same-source whole-Z Rh correction incoming required')
        for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,upstream.RECEIPT,sha(upstream.RECEIPT))
        self.saved=json.loads(gzip.decompress((HERE/upstream.NAME).read_bytes()))
        if not self.saved[upstream.GATE] or self.saved['candidate_N']!=self.N or self.saved['source_family']!=self.identity:
            raise ValueError('Same current whole-Z Rh packet and unchanged N required')
        self.saved_cells={tuple(row['exact_Z_cell']):row for row in self.saved['source_cells']}
        if set(self.saved_cells)!=set(upstream.upstream.upstream.source.CELLS):
            raise ValueError('Complete original whole-Z Rh correction cover required')
        self.eta_log=self.source.eta_log;self.dstar_log=self.source.dstar_log;self.sc=self.source.sc
        self.bindings=original.source_bindings()
        # This old metadata describes a different native formal rebase. Its
        # phase recipe is not transported to the current whole-Z algebra.
        self.bindings.pop('exact_reference_phase')
        self.bindings.update(actual_current_global_radius_phase='frac(N*(logRm+6+offset-logRa-hb*s_c/2))',
            original_global_positive_phase_origin_retained=True,
            source_admission='current whole-Z unique leading patch full-weight terminal identities',
            same_original_reference_source_with_actual_whole_Z_Rh_inlet=True)
        self.cache={};self.inlets={};self.joins={}
        for module in (upstream,original,reference,reference.background,reference.long,reference.parameters,
                       reference.primitives,reference.phase,reference.bounds):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def owner(self,ends):return self.source.owner(ends)

    def inlet(self,ends):
        key=tuple(ends)
        if key in self.inlets:return self.inlets[key]
        op=self.owner(ends);f,c=op.flow,op.c;stored=self.saved_cells[key]
        if stored['source_identity']!=self.identity or stored['candidate_N']!=self.N \
                or not stored['full_same_N_correction_prefix_through_Rh_installed']:
            raise ValueError('Genuine current same-N Rh source required')
        decode=lambda row:upstream.prefix.upstream.decode(f,row)
        P0=decode(stored['exact_common_P0_axial5'])
        def equal(live,old):
            return live.coefficient._mpi_==old.coefficient._mpi_ and live.scale.powers==old.scale.powers \
                and live.scale.offset._mpi_==old.scale.offset._mpi_
        if len(P0)!=6 or any(not equal(live,old) for live,old in zip(op.P0,P0)):
            raise ValueError('Same independent original P0 source tuples required')
        if not equal(op.Rm_factor,decode(stored['exact_Rm_radius'])) \
                or not equal(op.Rm_factor*c.exp(1),decode(stored['exact_Rh_radius'])):
            raise ValueError('Same exact current Rm and Rh radii required')
        incoming=decode(stored['actual_Rh_correction_C0_Z'])
        if set(incoming)!=set(RATES) or any(len(rows)!=2 for rows in incoming.values()):
            raise ValueError('Actual Rh correction-only five C0/Z rows required')
        upstream.parameters.same_source(f,[row for rows in incoming.values() for row in rows])
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.P0,exact_Rm_radius=op.Rm_factor,exact_Rh_radius=op.Rm_factor*c.exp(1),
            actual_Rh_correction_C0_Z=incoming,incoming_source_report=upstream.NAME,
            correction_only_not_complete_own_history=True,original_positive_phase_origin_s_c=self.sc,
            no_old_N257_label_incoming_or_phase_receipt_transplanted=True)
        self.inlets[key]=result;return result

    def query(self,ends,left,right=None):
        right=left if right is None else right;key=(tuple(ends),left,right)
        if key in self.cache:return self.cache[key]
        op=self.owner(ends);f,c=op.flow,op.c
        background=original.background_cell(op,left,right)
        if background['original_P0_normalized_axial5'] is not op.P0:
            raise ValueError('Same live independent reference P0 required')
        proxy,raw,recovered=reference.recover_cell(op.reference,background)
        recovered.pop('source_frame_conditional_on_same_current_Rsh_background')
        recovered.update(source_from_same_current_whole_Z_repaired_Rh_background=True,
            leading_background_not_finite_N_correction=True)
        roots,qr,proof=reference.quotient_cell(proxy,recovered,self.eta_log)
        got=reference.primitives.all_u_primitive_bounds(f,roots,qr,self.dstar_log,c.mpf([0,1]))
        if qr[C0].zero and qr[Z].zero:
            got['values']={name:f.scalar(0) for name in got['values']}
            got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
        A=got['values']['A'];logA=None if A.zero else ep(A.record()['log_absolute_upper'])[1]
        if logA is not None and logA>=ep(c.ln(self.N))[0]:
            raise ValueError('Same N fails actual reference A/N budget at '+str(key))
        got=upstream.prefix.upstream.upstream.bounded_exponent_cover(f,got,self.N)
        E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
        density=reference.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],self.N)
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_left=list(left),exact_right=list(right),exact_common_P0_axial5=op.P0,
            actual_background_source=background,original_raw_source=raw,original_generic_source=recovered,
            original_full_source_quotients=proof,original_roots=roots['roots'],original_q_C0_Z=qr,
            original_primitive_values=got['values'],original_primitive_proof=got['record'],
            actual_original_A_log_absolute_upper=None if logA is None else c.mpf(logA),
            actual_A_below_same_candidate_N=True,original_signed_five_density_C0_Z=density,
            actual_global_radius_phase='frac(N*(logRm+6+offset-logRa-hb*s_c/2))',
            actual_phase_origin_s_c=self.sc,actual_phase_Z_exact_zero=True,
            actual_radius_phase_full_period_cover=c.mpf([0,1]),phase_not_restarted_at_Rh=True,
            phase_average_cancellation_not_claimed=True,
            reference_physical_log_measure='dy=d(offset)',
            original_reference_a_four_fifths_only_after_patch_terminal=True,
            signed_axis_no_absZ_or_positive_rho_division=True,
            actual_high_order_finite_N_correction_jets_installed=False)
        self.cache[key]=result;return result

    def seam(self,ends):
        key=tuple(ends)
        if key in self.joins:return self.joins[key]
        op=self.owner(ends);f=op.flow
        leading_packet=self.source.source.query(ends,'Rh')
        _,_,leading,_,_=upstream.recover_patch(op,leading_packet)
        tail=self.query(ends,(-5,1))['original_generic_source'];count=0
        pairs=[(leading['common_velocity_E_axial5'],tail['common_velocity_E_axial5']),
            (leading['common_velocity_V_axial5'],tail['common_velocity_V_axial5'])]
        pairs.extend((leading['common_own_five_histories_axial5'][name],tail['common_own_five_histories_axial5'][name]) for name in RATES)
        for before,after in pairs:
            for a,b in zip(before,after):
                difference=a-b
                if not difference.zero and not ep(difference.coefficient)[0]<=0<=ep(difference.coefficient)[1]:
                    raise ValueError('Defining whole-Z leading Rh/reference source covers do not meet')
                count+=1
        result=dict(source_function_identity_from_original_full_weight_unique_leading_map=True,
            exact_radius_identity='Rm*e=Rm*exp(6-5)=Rref*exp(-5)',
            exact_angular_amplitude_identity='exp(-.6)*e^.1=exp(-.5)',
            same_five_closed_histories_and_original_pressure=True,
            shared_original_P0=leading['common_original_P0_axial5'] is op.P0,
            directed_source_overlap_consistency_rows=count,
            overlap_only_consistency_not_functional_identity_proof=True)
        self.joins[key]=result;return result

    def contribution(self,ends):
        op=self.owner(ends);f,c=op.flow,op.c;inlet=self.inlet(ends);join=self.seam(ends)
        incoming=inlet['actual_Rh_correction_C0_Z']
        total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
        for left,right in zip(PARTITION,PARTITION[1:]):
            source=self.query(ends,left,right);density=source['original_signed_five_density_C0_Z']
            rows={};weights={}
            for name,rate in RATES.items():
                width,suffix,mass,decay,tail=reference_weights(f,left,right,rate)
                reference.parameters.positive_source(f,mass,'actual_Rh_reference_Duhamel_mass')
                pair=[reference.bounds.symmetric(f,reference.bounds.magnitude(f,density[part][name])*mass*tail)
                    for part in ('kernels','Z_derivatives')]
                for n,row in enumerate(pair):total[name][n]+=row
                rows[name]=pair;weights[name]=dict(true_full_mass=mass,true_cell_decay=decay,
                    true_suffix_decay=tail,own_rate=str(rate),true_radial_measure_applied_once=True)
            cells.append(dict(source=source,actual_log_width=width,actual_suffix_to_Rref=suffix,
                signed_cell_driver_C0_Z=rows,own_rate_weights=weights))
            print('Whole-Z actual same-N Rh reference driver: '+str(ends)+' '+str(left)+' '+str(right),flush=True)
        memory={name:upstream.long.kernel_weight(f,c.mpf(5),c.mpf(0),rate)[1] for name,rate in RATES.items()}
        inherited={name:[row*memory[name] for row in incoming[name]] for name in RATES}
        outgoing={name:[inherited[name][n]+total[name][n] for n in range(2)] for name in RATES}
        terminal=self.query(ends,(0,1))
        background={name:list(rows[:2]) for name,rows in terminal['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete={name:[background[name][n]+outgoing[name][n] for n in range(2)] for name in RATES}
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.P0,exact_Rm_radius=op.Rm_factor,exact_Rh_radius=op.Rm_factor*c.exp(1),
            exact_Rref_radius=op.Rm_factor*c.exp(6),actual_Rh_incoming_binding=inlet,
            exact_total_log_length=5,actual_source_cells=cells,typed_actual_Rh_reference_source_join=join,
            actual_Rh_correction_incoming_C0_Z=incoming,incoming_own_rate_memory=memory,
            actual_local_driver_C0_Z=total,actual_retained_incoming_C0_Z=inherited,
            actual_Rref_correction_C0_Z=outgoing,actual_original_Rref_background_C0_Z=background,
            actual_Rref_complete_own_history_C0_Z=complete,actual_original_Rref_source=terminal,
            full_same_N_correction_prefix_through_Rref_installed=True,
            original_P0_background_and_finite_N_correction_kept_distinct=True,
            **dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZRhReferenceFiniteN()
        cells=[serialized(owner.contribution(ends)) for ends in upstream.upstream.upstream.source.CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_domain=['-1','1'],exact_Z_partition=[list(ends) for ends in upstream.upstream.upstream.source.CELLS],
            exact_reference_offset_domain=['-5','0'],exact_reference_partition=PARTITION,source_cells=cells,
            original_source_bindings=owner.bindings,original_eta_log=owner.eta_log,original_dstar_log=owner.dstar_log,
            full_same_N_correction_prefix_through_Rref_installed=True,
            actual_high_order_finite_N_correction_jets_installed=False,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual whole-Z original repaired Rh reference functions, full signed C0/Z '
                  'drivers and genuine unchanged-N incoming to Rref. Later O2/O3/Rc corrections, '
                  'high correction jets/final repair, global N/cone/heat and n-recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(upstream.upstream.upstream.source.bridge._encode(report),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z genuine same-N correction prefix reaches Rref',flush=True);return report


if __name__=='__main__':run()

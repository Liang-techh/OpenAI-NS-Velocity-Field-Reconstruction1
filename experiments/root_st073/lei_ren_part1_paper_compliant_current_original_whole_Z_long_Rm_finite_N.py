"""Whole-Z actual long/reference source and same-N R110 correction to Rm."""
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_switch_finite_N as upstream

current=upstream.current;long=upstream.switch.long;reference=upstream.switch.downstream
background=reference.background;endpoint=background.previous;kernels=endpoint.kernels
HERE,PREFIX,sha,read,bind,ep=upstream.HERE,upstream.PREFIX,upstream.sha,upstream.read,upstream.bind,upstream.ep
phase,primitives,bounds,parameters=upstream.phase,upstream.primitives,upstream.bounds,upstream.parameters
RATES,C0,Z=upstream.RATES,upstream.C0,upstream.Z
NAME=PREFIX+'current_original_whole_Z_long_Rm_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_long_Rm_finite_N_check.json'
GATE='current_original_whole_Z_same_N_actual_long_reference_source_and_correction_through_Rm_enclosed'
OPEN=upstream.OPEN
PARTITION=upstream.PARTITION
CHARTS=('long_reshape','reference','restoration','postrestore')


class WholeZLongRmFiniteN:
    def __init__(self,dps=500):
        self.source=upstream.WholeZSwitchFiniteN(dps)
        self.c=c=self.source.c;self.identity=self.source.identity;self.N=self.source.N
        self.hashes=dict(self.source.hashes)
        receipt=json.loads((HERE/upstream.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(upstream.GATE) or receipt['source_family']!=self.identity:
            raise ValueError('Checked same-source whole-Z R110 input required')
        for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,upstream.RECEIPT,sha(upstream.RECEIPT))
        self.saved=json.loads(gzip.decompress((HERE/upstream.NAME).read_bytes()))
        if self.saved['candidate_N']!=self.N or self.saved['source_family']!=self.identity:
            raise ValueError('Same genuine N and current R110 source required')
        self.saved_cells={tuple(row['exact_Z_cell']):row for row in self.saved['source_cells']}
        if set(self.saved_cells)!=set(current.source.CELLS):raise ValueError('Complete whole-Z R110 cover required')
        params=json.loads((HERE/kernels.PARAMS).read_bytes())
        self.A=read(c,params['A_upper']);self.T=400*self.A;self.logC=read(c,params['selected_logCstar'])
        join_name=PREFIX+'reference_join_bounds.json';join=json.loads((HERE/join_name).read_bytes())
        if read(c,join['B_C2_upper'])._mpi_!=(2*self.A)._mpi_:
            raise ValueError('Same original B norm theorem and frozen T required')
        kernel_receipt=json.loads((HERE/kernels.RECEIPT).read_bytes())
        if not kernel_receipt.get('all_passed') or not kernel_receipt.get(kernels.GATE) \
                or kernel_receipt['source_family']!=self.identity['actual_five_defect_family_sha256']:
            raise ValueError('Checked actual finite terminal kernel backend required')
        for name,digest in kernel_receipt['input_hashes'].items():bind(self.hashes,name,digest)
        admitted=json.loads((HERE/background.ADMISSION).read_bytes())
        checked=json.loads((HERE/background.ADMISSION_CHECK).read_bytes())
        if not checked.get('all_passed') or not checked.get('original_restoration_integrals_directly_source_bound'):
            raise ValueError('Original complete restoration scalar integrals required')
        for key,want in self.identity.items():
            if admitted[key]!=want:raise ValueError('Same original restoration source/datum/family required')
        for name,digest in checked['input_hashes'].items():bind(self.hashes,name,digest)
        self.provider=background.RestorationKernelProvider(c,{name:read(c,row)
            for name,row in admitted['directed_restoration_integrals'].items()})
        self.bindings=dict(long_recovery=long.RECOVERY_BINDING,
            reference=reference.source_bindings(),endpoint=endpoint.source_bindings())
        self.eta_log=self.source.eta_log;self.dstar_log=self.source.dstar_log
        self.owners={};self.cache={}
        for name in (kernels.PARAMS,join_name,kernels.RECEIPT,background.ADMISSION,background.ADMISSION_CHECK):
            bind(self.hashes,name,sha(name))
        for module in (upstream,long,reference,background,endpoint,kernels):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def owner(self,ends):
        key=tuple(ends)
        if key in self.owners:return self.owners[key]
        if key not in self.saved_cells:raise ValueError('Actual checked whole-Z R110 cell required')
        source=self.source.background.source
        f,proof=source.owner(ends);c=self.c;z=proof['Z'];incoming=self.saved_cells[key]
        fields=upstream.decode(f,incoming['actual_R110_background_field_source'])
        histories=upstream.decode(f,incoming['actual_R110_background_six_history_source'])
        p0=current.source.IntervalTaylor(c,proof['original_pressure_Z0_through_Z6'][:6])
        P0=f.jet(p0);savedP0=upstream.decode(f,incoming['exact_common_P0_axial5'])
        if any(live.coefficient._mpi_!=old.coefficient._mpi_ or live.scale.powers!=old.scale.powers
               or live.scale.offset._mpi_!=old.scale.offset._mpi_ for live,old in zip(P0,savedP0)):
            raise ValueError('Same exact original pressure datum tuples required')
        # B is the original anchored log-shape, distinct from history B.
        amplitude=proof['actual_anchored_amplitude']
        if f.logs[2]._mpi_!=(2*amplitude['logF0'])._mpi_:
            raise ValueError('Same actual anchored logF0 source required')
        gradient=source.amplitude.gradient_jets(z,5)
        # gradient_jets is the Taylor series of G', not of G itself.
        G=current.source.IntervalTaylor(c,[amplitude['G']]
            +[gradient[n-1]/n for n in range(1,6)])
        phi=current.source.IntervalTaylor(c,[f.ordinary_cover(row) for row in fields['phi']])
        if ep(phi[0])[0]<=0:raise ValueError('Positive actual R110 angular source required')
        zj=current.source.IntervalTaylor.variable(c,z,5)
        B=G*(-source.core.Lambda)+endpoint.logarithm(phi)+endpoint.logarithm(1+zj*zj)+c.ln(220)/2
        for n in range(3):
            if max(abs(v) for v in ep(B[n]))>ep(2*self.A/math.factorial(n))[1]:
                raise ValueError('Actual whole-Z B C2 source bound needs refinement: '+str(ends)+' order'+str(n))
        long_op=endpoint.ActualLongReshapeEndpoint(f,z,source.core.delta,fields['phi'],fields['V'],
            histories,B,self.T,self.logC,p0)
        terminal=long_op.evaluate()
        ref=background.ActualReferenceRestoreFunctions(f,z,source.core.delta,
            terminal['actual_terminal_normalized_six_history_shapes'],fields['V'],P0,
            self.T,self.logC,self.provider)
        correction=upstream.decode(f,incoming['actual_R110_correction_C0_Z'])
        result=SimpleNamespace(flow=f,c=c,Z=z,long=long_op,reference=ref,terminal=terminal,
            R110_fields=fields,R110_histories=histories,incoming=correction,source_proof=proof,
            original_B_source=dict(actual_anchored_G=G,actual_log_phi=endpoint.logarithm(phi),
                actual_log_shape_B=B,original_Lambda=source.core.Lambda,
                B_distinct_from_history_B=True,actual_B_C2_bound_verified_without_field_cap=True))
        self.owners[key]=result
        print('Whole-Z actual long/reference owner: '+str(ends),flush=True)
        return result

    def query(self,ends,chart,left,right=None):
        if chart not in CHARTS:raise ValueError('Original long/reference source chart required')
        right=left if right is None else right;key=(tuple(ends),chart,left,right)
        if key in self.cache:return self.cache[key]
        owner=self.owner(ends);f,c=owner.flow,owner.c
        l,r=long.fraction(left),long.fraction(right)
        if r<l:raise ValueError('Ordered original source phase cell required')
        t=reference.scalar_hull(c,c.mpf(l.numerator)/l.denominator,c.mpf(r.numerator)/r.denominator)
        if chart=='long_reshape':
            y=self.T*t;radius=f.factor((0,0,0,0,0),c.ln(110)+y)
            proxy,raw,recovered,a,source=long.source_packet(owner.long,owner.reference,radius,t,y)
            roots,qr,proof=long.source_quotients(proxy,recovered,a,self.eta_log)
            length=self.T
        else:
            source=reference.background_cell(owner.reference,chart,t)
            proxy,raw,recovered=reference.recover_cell(owner.reference,source)
            roots,qr,proof=reference.quotient_cell(proxy,recovered,self.eta_log)
            length=source['original_window_length']
        got=primitives.all_u_primitive_bounds(f,roots,qr,self.dstar_log,c.mpf([0,1]))
        if qr[C0].zero and qr[Z].zero:
            got['values']={name:f.scalar(0) for name in got['values']}
            got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
        A=got['values']['A'];logA=None if A.zero else ep(A.record()['log_absolute_upper'])[1]
        if logA is not None and ep(c.ln(self.N))[0]<=logA:
            raise ValueError('Same R110 N fails actual downstream exponent budget at '+str(key)+'; logAupper='+mp.nstr(logA,18))
        got=upstream.upstream.bounded_exponent_cover(f,got,self.N)
        E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
        density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],self.N)
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            chart=chart,exact_left=list(left),exact_right=list(right),original_window_length=length,
            exact_common_P0_axial5=owner.reference.P0,actual_background_source=source,
            original_raw_source=raw,original_generic_source=recovered,original_full_source_quotients=proof,
            original_q_C0_Z=qr,original_primitive_values=got['values'],original_primitive_proof=got['record'],
            actual_original_A_log_absolute_upper=None if logA is None else c.mpf(logA),
            original_signed_five_density_C0_Z=density,
            actual_radius_phase_definition='frac(N*(logR-logRa-hb*s_c/2))',
            actual_radius_phase_Z_exact_zero=True,actual_radius_phase_full_period_cover=c.mpf([0,1]),
            phase_average_cancellation_not_claimed=True,actual_same_N_R110_incoming_available=True)
        self.cache[key]=result;return result

    def contribution(self,ends):
        op=self.owner(ends);f,c=op.flow,op.c
        incoming={name:list(rows) for name,rows in op.incoming.items()};windows=[]
        for chart in CHARTS:
            total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
            length=self.T if chart=='long_reshape' else op.reference.gap if chart=='reference' else c.mpf(1)
            for left,right in zip(PARTITION,PARTITION[1:]):
                query=self.query(ends,chart,left,right);density=query['original_signed_five_density_C0_Z']
                l,r=long.fraction(left),long.fraction(right)
                width=length*(c.mpf((r-l).numerator)/(r-l).denominator)
                suffix=length*(c.mpf((1-r).numerator)/(1-r).denominator)
                rows={};weights={}
                for name,rate in RATES.items():
                    mass,decay,tail=long.kernel_weight(f,width,suffix,rate)
                    parameters.positive_source(f,mass,'actual_whole_Z_downstream_Duhamel_mass')
                    pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*tail)
                          for part in ('kernels','Z_derivatives')]
                    for n,row in enumerate(pair):total[name][n]+=row
                    rows[name]=pair;weights[name]=dict(true_full_mass=mass,true_cell_decay=decay,
                        true_suffix_decay=tail,own_rate=str(rate),true_radial_measure_applied_once=True)
                cells.append(dict(source=query,signed_cell_driver_C0_Z=rows,own_rate_weights=weights))
            memory={name:long.kernel_weight(f,length,c.mpf(0),rate)[1] for name,rate in RATES.items()}
            inherited={name:[row*memory[name] for row in incoming[name]] for name in RATES}
            outgoing={name:[inherited[name][n]+total[name][n] for n in range(2)] for name in RATES}
            windows.append(dict(chart=chart,actual_window_length=length,source_cells=cells,
                actual_incoming_correction_C0_Z=incoming,incoming_own_rate_memory=memory,
                actual_local_driver_C0_Z=total,actual_retained_incoming_C0_Z=inherited,
                actual_exit_correction_C0_Z=outgoing))
            incoming=outgoing
            print('Whole-Z same-N downstream window: '+str(ends)+' '+chart,flush=True)
        terminal=self.query(ends,'postrestore',(1,1))
        original={name:list(rows[:2]) for name,rows in terminal['original_generic_source']['common_own_five_histories_axial5'].items()}
        complete={name:[original[name][n]+incoming[name][n] for n in range(2)] for name in RATES}
        Rm=f.factor((0,5,0,0,0),10*self.logC+c.ln(110)-6)
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            exact_common_P0_axial5=op.reference.P0,actual_R110_correction_incoming_C0_Z=op.incoming,
            actual_original_anchored_B_source=op.original_B_source,
            actual_original_finite_long_endpoint=op.terminal,
            actual_source_windows=windows,actual_Rsh_correction_C0_Z=windows[0]['actual_exit_correction_C0_Z'],
            actual_Rm_correction_C0_Z=incoming,actual_original_Rm_background_C0_Z=original,
            actual_Rm_complete_own_history_C0_Z=complete,exact_Rm_radius=Rm,
            actual_original_Rm_background_function=op.reference.postrestore((-6,1)),
            exact_total_length_identity='T+gap+1+1=10*(logC+logP)-6',
            original_full_window_and_nonzero_incoming_memory_retained=True,
            actual_R110_background_not_replaced_by_finite_correction=True,
            whole_Z_leading_Rm_implicit_controls_solved=False,**dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZLongRmFiniteN()
        cells=[reference.serialized(owner.contribution(ends)) for ends in current.source.CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_domain=['-1','1'],exact_Z_partition=[list(ends) for ends in current.source.CELLS],
            source_cells=cells,original_source_bindings=owner.bindings,
            original_fixed_T=owner.T,original_fixed_logCstar=owner.logC,
            whole_Z_leading_Rm_implicit_controls_solved=False,**dict.fromkeys(OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual whole-Z original long/reference/restoration source and genuine same-N R110 '
                  'correction incoming through Rm. Background and correction remain distinct. '
                  'Whole-Z Rm implicit patch controls, Rref/Rc repair, sharpness/global N/cone/heat and recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(current.source.bridge._encode(report),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z genuine same-N correction prefix reaches Rm',flush=True)
    return report


if __name__=='__main__':run()

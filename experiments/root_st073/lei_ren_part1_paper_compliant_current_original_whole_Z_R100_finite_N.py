"""Source-owned whole-Z current finite-N correction histories to R100.

The defining unchanged-left collar supplies zero correction at r_minus,
while the original background histories are retained. All subsequent five
density increments use the live micro/macro source and their true own-rate
Volterra weights. No archived finite-N incoming or point interpolation.
"""
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_micro_macro as current
import lei_ren_part1_paper_compliant_current_original_micro_finite_N as original

HERE,PREFIX,sha,read,bind,ep=current.HERE,current.PREFIX,current.sha,current.read,current.bind,current.ep
macro,switch=current.macro,current.switch
phase,primitives,bounds,parameters=macro.phase,macro.primitives,macro.bounds,macro.parameters
RATES,C0,Z=macro.RATES,macro.C0,macro.Z
NAME=PREFIX+'current_original_whole_Z_R100_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_R100_finite_N_check.json'
GATE='current_original_source_owned_whole_Z_finite_N_prefix_through_R100_enclosed'
OPEN=current.OPEN


def bounded_exponent_cover(f,got,N):
    """Outward A cover bounded by this candidate N; preserve the A_Z source.

    The legacy 1000 cutoff is a materialization guard, not a source
    inequality. Here exp(log|A|) is materialized only after proving it is
    below the explicit at-most-4096-bit integer N.
    """
    row=got['values']['A']
    if row.zero:return got
    lower,upper=ep(row.scale.evaluate())
    if lower < -1000 <= upper:
        log_upper=ep(row.record()['log_absolute_upper'])[1]
        if not mp.isfinite(log_upper) or log_upper>=ep(f.c.ln(N))[0]:
            raise ValueError('Strict actual A/N exponent budget required')
        cap=f.c.exp(f.c.mpf(log_upper))
        got['record']['original_A_C0_before_outward_exponent_cover']=row.record()
        got['values']['A']=f.scalar(f.c.mpf([-ep(cap)[1],ep(cap)[1]]))
        got['record']['bounded_A_C0_cover_only_not_source_or_A_Z_replacement']=True
        got['record']['density_A_C0_outward_enclosure']=got['values']['A'].record()
        got['record']['materialized_A_cap_below_same_explicit_candidate_N']=True
    return got


class WholeZR100FiniteN:
    def __init__(self,dps=500):
        self.background=current.WholeZMicroMacro(dps)
        self.c=c=self.background.c
        self.identity=self.background.identity
        self.hashes=dict(self.background.hashes)
        self.records={}
        for stem,gate in (
            ('current_original_whole_Z_micro_macro',current.GATE),
            ('current_generic_shear_uniform_inputs','current_original_whole_generic_input_margins_and_log_scales_certified'),
            ('current_generic_shear_loop_domain','current_original_generic_loop_two_sided_collars_and_repair_geometry_certified')):
            checked_name=PREFIX+stem+'_check.json'
            receipt=json.loads((HERE/checked_name).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(gate) or receipt['source_family']!=self.identity:
                raise ValueError('Same-source current finite-N prerequisite required: '+stem)
            for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
            bind(self.hashes,checked_name,sha(checked_name))
            if stem!='current_original_whole_Z_micro_macro':
                name=PREFIX+stem+'.json'
                self.records[stem]=json.loads((HERE/name).read_bytes())
                bind(self.hashes,name,sha(name))
        scales=self.records['current_generic_shear_uniform_inputs']['current_actual_logarithmic_loop_scales']
        domain=self.records['current_generic_shear_loop_domain']['current_original_generic_loop_domain']
        self.eta_log=read(c,scales['selected_positive_eta_log'])
        self.dstar_log=read(c,scales['logarithmic_selected_positive_lower_constants']['d_star'])
        collar=self.background.records['current_inner_exit_strict_collar']['explicit_current_inner_exit_strict_collar']
        self.sc=read(c,collar['selected_first_phase_endpoint'])
        if ep(self.sc)[0]!=ep(self.sc)[1] or not 0<ep(self.sc)[0]<mp.mpf('.25'):
            raise ValueError('Exact original positive compact phase origin below1/4 required')
        if read(c,collar['exact_source_width_log_enclosure'])._mpi_!=self.background.logh._mpi_:
            raise ValueError('Same exact original source width required')
        excess=read(c,domain['left_collar_kappa_minus2_lower'])
        if ep(excess)[0]<=0 or ep(self.eta_log)[1]>ep(c.ln(excess)-c.ln(2))[0] \
                or ep(self.eta_log)[1]>ep(read(c,domain['allowed_eta_log_upper_for_automatic_constant_edges']))[0]:
            raise ValueError('Same selected eta must prove the actual flat left join')
        self.inlet_proof=dict(source_owned_definition='original background unchanged before r_minus',
            exact_radius='r_minus=Ra*exp(hb*s_c/2)',selected_positive_s_c=self.sc,
            original_left_kappa_excess=excess,selected_eta_log=self.eta_log,
            selected_dstar_log=self.dstar_log,
            exact_zero_correction_follows_from_defining_Duhamel_lower_bound=True,
            actual_background_incoming_histories_remain_nonzero=True,
            q_A_B_and_Z_rows_flat_on_same_source_two_sided_left_collar=True,
            no_saved_frame_or_N1024_initial_correction_used=True)
        self.cache={};self.prepared={}
        for module in (original,phase,primitives,bounds,parameters):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def support(self,chart,left,right):
        c=self.c
        if chart=='frozen_macro':return 'original_modified_right'
        if chart not in current.micro.CHARTS:raise ValueError('Original source chart required')
        lo,hi=ep(c.mpf(left))[0],ep(c.mpf(right))[1]
        if chart=='first_micro' and hi<=ep(self.sc/2)[0]:return 'unchanged_left'
        if chart=='first_micro' and lo>=ep(self.sc/4)[0] and hi<=ep(self.sc*3/4)[0]:return 'original_flat_left_collar'
        if chart=='first_micro' and lo<ep(self.sc/2)[1]<hi:return 'unchanged_modified_union'
        return 'original_modified_right'

    def prepare(self,ends,chart,left,right):
        op=self.background.owner(ends);f,c=op.flow,op.c
        key=(op.Z._mpi_,chart,str(left),str(right))
        if key in self.prepared:return self.prepared[key]
        value=self.background.generic(ends,chart,left,right)
        r=value['full_signed_generic_source']
        support=self.support(chart,left,right)
        if support in ('unchanged_left','original_flat_left_collar'):
            q={C0:f.scalar(0),Z:f.scalar(0)}
            got=dict(values={key:f.scalar(0) for key in phase.OUTPUTS},
                record=dict(exact_flat_original_A_B_and_Z_zero_by_source_owned_left_collar=True,
                            active_inverse_not_evaluated_on_original_flat_collar=True))
            proof=dict(source_owned_flat_left_collar=True,original_q_C0_Z_exact_zero=True)
        else:
            roots,q,proof=switch.general_quotients(value['proxy'],r,value['a'],self.eta_log)
            got=primitives.all_u_primitive_bounds(f,roots,q,self.dstar_log,c.mpf([0,1]))
            if q[C0].zero and q[Z].zero:
                got['values']={key:f.scalar(0) for key in got['values']}
                got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
        self.prepared[key]=dict(value=value,q=q,got=got,proof=proof,support=support)
        return self.prepared[key]

    def query(self,ends,chart,left,right,N=257):
        N=phase.candidate_N(N)
        op=self.background.owner(ends);f,c=op.flow,op.c
        key=(op.Z._mpi_,chart,str(left),str(right),N)
        if key in self.cache:return self.cache[key]
        saved=self.prepare(ends,chart,left,right)
        value,q,proof,support=(saved[key] for key in ('value','q','proof','support'))
        got=dict(values=dict(saved['got']['values']),record=dict(saved['got']['record']))
        got=bounded_exponent_cover(f,got,N)
        r=value['full_signed_generic_source']
        E,V=r['common_velocity_E_axial5'],r['common_velocity_V_axial5']
        density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],N)
        source=value['background']
        unwrapped=(source['exact_log_radius_over_R0']+f.h*(2-self.sc/2)
                   if chart=='frozen_macro' else source['exact_log_radius_over_Ra']-f.h*(self.sc/2))
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=N,
            original_background_source=source,common_original_P0=op.reference.P0,
            full_signed_generic_source=r,original_full_quotient_proof=proof,
            original_q_C0_Z=q,original_primitive_values=got['values'],original_primitive_proof=got['record'],
            original_signed_five_density_C0_Z=density,source_owned_modification_support=support,
            actual_unwrapped_log_radius_phase=unwrapped,actual_phase_definition='frac(N*log(R/r_minus))',
            actual_phase_Z_exact_zero=True,full_period_phase_cover=c.mpf([0,1]),
            phase_average_cancellation_not_claimed=True)
        self.cache[key]=result
        return result

    def partitions(self):
        c=self.c
        return (('first_micro',(self.sc/2,c.mpf('.25'),c.mpf('.5'),c.mpf('.75'),c.mpf(1)),c.mpf(1)),
                ('second_micro',(c.mpf(1),c.mpf('1.25'),c.mpf('1.5'),c.mpf('1.75'),c.mpf(2)),c.mpf(2)),
                ('frozen_macro',((0,1),(1,2),(1,1)),None))

    def select_N(self):
        """One fresh dyadic N from actual whole-prefix A bounds; no transplant."""
        c=self.c;maximum=-mp.inf;requirements=[]
        for ends in current.source.CELLS:
            for chart,partition,_ in self.partitions():
                for left,right in zip(partition,partition[1:]):
                    got=self.prepare(ends,chart,left,right)['got'];A=got['values']['A']
                    upper=None if A.zero else ep(A.record()['log_absolute_upper'])[1]
                    if upper is not None:
                        if not mp.isfinite(upper):
                            raise ValueError('Whole-prefix A bound needs correlated refinement: '
                                +str((ends,chart,str(left),str(right)))+'; logAupper='+mp.nstr(upper,18))
                        maximum=max(maximum,upper)
                    requirements.append(dict(exact_Z_cell=list(ends),chart=chart,left=left,right=right,
                        exact_A_zero=A.zero,log_A_absolute_upper=None if upper is None else c.mpf(upper)))
            print('Whole-Z original prefix primitive requirements: '+str(ends),flush=True)
        power=9 if maximum==-mp.inf else max(9,int(mp.ceil(maximum/mp.log(2)))+1)
        if power>=4096:
            raise ValueError('Whole-prefix needs '+str(power+1)+' N bits; refine the source cells '
                             'before exceeding the current 4096-bit candidate interface')
        N=phase.candidate_N(2**power)
        margin=None if maximum==-mp.inf else c.ln(N)-c.mpf(maximum)
        if margin is not None and ep(margin)[0]<=0:raise ValueError('Strict common exponent budget required')
        return N,dict(candidate_N=N,exact_dyadic_power=power,whole_prefix_A_requirements=requirements,
            maximum_log_A_absolute_upper=None if maximum==-mp.inf else c.mpf(maximum),
            log_N_minus_whole_prefix_log_A_upper=margin,
            actual_whole_prefix_A_over_N_below_one=True,
            N_selected_from_live_source_not_archived_N1024_or_point_frames=True,
            source_N257_not_reused_if_its_exponent_budget_fails=True,
            global_all_region_N_cone_or_terminal_repair_admission=False)

    def contribution(self,ends,N=257):
        N=phase.candidate_N(N)
        op=self.background.owner(ends);f,c=op.flow,op.c
        start=self.sc/2
        inlet=self.query(ends,'first_micro',start,start,N)
        density=inlet['original_signed_five_density_C0_Z']
        assert all(row.zero for part in ('kernels','Z_derivatives') for row in density[part].values())
        zero={name:[f.scalar(0),f.scalar(0)] for name in RATES}
        incoming={name:list(rows) for name,rows in zero.items()}
        windows=[]
        for chart,partition,end in self.partitions():
            total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
            for left,right in zip(partition,partition[1:]):
                source=self.query(ends,chart,left,right,N);rows={};weight_rows={}
                for name,rate in RATES.items():
                    weights=(macro.weights(op.series,left,right,rate) if chart=='frozen_macro' else
                             original.weights(op.series,left,right,end,rate))
                    mass,decay,tail=weights
                    parameters.positive_source(f,mass,'current_whole_Z_true_cell_Duhamel_mass')
                    density=source['original_signed_five_density_C0_Z']
                    pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*tail)
                          for part in ('kernels','Z_derivatives')]
                    for n,row in enumerate(pair):total[name][n]+=row
                    rows[name]=pair;weight_rows[name]=dict(true_full_mass=mass,true_cell_decay=decay,
                        true_suffix_decay=tail,own_rate=str(rate),true_radial_measure_applied_once=True)
                cells.append(dict(source=source,signed_cell_driver_C0_Z=rows,own_rate_weights=weight_rows))
            memory={name:(macro.weights(op.series,partition[0],partition[-1],rate)[1]
                          if chart=='frozen_macro' else
                          original.weights(op.series,partition[0],partition[-1],end,rate)[1]) for name,rate in RATES.items()}
            inherited={name:[row*memory[name] for row in incoming[name]] for name in RATES}
            outgoing={name:[inherited[name][n]+total[name][n] for n in range(2)] for name in RATES}
            windows.append(dict(chart=chart,source_cells=cells,actual_incoming_correction_C0_Z=incoming,
                incoming_own_rate_memory=memory,actual_local_driver_C0_Z=total,
                actual_retained_incoming_C0_Z=inherited,actual_exit_correction_C0_Z=outgoing))
            incoming=outgoing
            print('Whole-Z current N'+str(N)+' correction window: '+str(ends)+' '+chart,flush=True)
        terminal=self.background.generic(ends,'frozen_macro',(1,1),(1,1))
        original_histories=terminal['full_signed_generic_source']['common_own_five_histories_axial5']
        background={name:list(rows[:2]) for name,rows in original_histories.items()}
        own={name:[background[name][n]+incoming[name][n] for n in range(2)] for name in RATES}
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=N,
            exact_common_P0_axial5=op.reference.P0,actual_source_owned_initial_condition=self.inlet_proof,
            actual_current_left_background=inlet['full_signed_generic_source'],
            actual_current_inlet_correction_C0_Z=zero,actual_source_windows=windows,
            actual_current_R0_correction_C0_Z=windows[1]['actual_exit_correction_C0_Z'],
            actual_current_R100_correction_C0_Z=incoming,actual_original_R100_background_C0_Z=background,
            actual_current_R100_complete_own_history_C0_Z=own,exact_R100_radius=f.scalar(100),
            actual_micro_macro_source_join=self.background.seam(ends),
            genuine_current_finite_N_incoming_propagated_not_reset=True,
            candidate_N_is_not_global_N_or_cone_admission=True,**dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZR100FiniteN()
        N,selection=owner.select_N()
        cells=[original.serialized(owner.contribution(ends,N)) for ends in current.source.CELLS]
        result=dict(**{GATE:True},source_family=owner.identity,candidate_N=N,
            current_whole_prefix_N_selection=original.serialized(selection),
            exact_Z_domain=['-1','1'],exact_Z_partition=[list(cell) for cell in current.source.CELLS],
            source_cells=cells,source_owned_left_inlet_proof=owner.inlet_proof,
            same_current_N_and_live_implicit_source_used_throughout=True,
            no_saved_frame_N1024_boundary_or_two_point_interpolation_used=True,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Genuine whole-Z fresh common finite-N source-owned correction histories r_minus..R0..R100, '
                  'all five C0/Z signed density integrals and true own-rate memory. '
                  'R100 background, correction and complete-own incoming exported separately. '
                  'Downstream whole-Z Rm/Rref/Rc repair, sharpness, global N/cone/heat and recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(current.source.bridge._encode(result),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Whole-Z genuine fresh common finite-N correction prefix reaches R100',flush=True)
    return result


if __name__=='__main__':run()

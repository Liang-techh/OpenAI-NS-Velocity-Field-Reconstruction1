"""Original 11-unit buffer: actual nonlinear period covers and Rd/Rc memory.

Frozen C0 means use exact reflection, with explicit full slow-envelope and
Duhamel-weight errors. Z integrates original nonlinear phase-bin covers.
No accepted numerical owner/ancestor producer is constructed or rerun.
"""
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_axial_endpoint_native_source as axial

native=axial.native;complete=axial.complete;rc=axial.rc;prior=axial.prior
common=axial.common;current=axial.current;supported=axial.supported;first=axial.first
Taylor=axial.Taylor;HERE,PREFIX,sha,encode,ep,iv=axial.HERE,axial.PREFIX,axial.sha,axial.encode,axial.ep,axial.iv
ZERO=axial.ZERO;ORDERS=axial.ORDERS;KEYS=axial.KEYS;N=axial.N
NAME=PREFIX+'current_buffer_period_transport.json';RECEIPT=PREFIX+'current_buffer_period_transport_check.json'
GATE='original_complete_buffer_actual_period_density_C1_transport_and_new_Rc_executed'
CELL_COUNT=22;PHASE_BINS=16


def restore_taylor(rows,coordinates):
    return Taylor(coordinates.scalar(1),[complete.restore_half_source(v,coordinates) for v in rows])


class OriginalBuffer:
    def __init__(self,coordinates,data):
        self.coordinates=coordinates;self.c=c=coordinates.ctx;self.t=coordinates.scalar(1)
        self.family=coordinates.family;self.Z=data['exact_Z_range'];self.data=data
        parent=data['original_endpoint_parent_binding'];inputs=parent['original_saved_native_inputs']
        self.delta=iv(c,inputs['original_delta']);self.P0=restore_taylor(parent['original_separate_P0'],coordinates)
        self.z=Taylor(self.t,[c.mpf(self.Z),1,0,0,0,0])
        self.EY=restore_taylor(parent['original_Utheta_endpoint_values']['right'],coordinates)
        self.initial={k:restore_taylor(v,coordinates) for k,v in parent['original_right_history_values'].items()}
        self.invPs=prior.ScaledEnclosure(prior.FormalScale(self.t.scale.bases,(0,0,0,-2,0)),1,self.t.ledger)
        self.invP2=prior.ScaledEnclosure(prior.FormalScale(self.t.scale.bases,(0,0,0,-4,0)),1,self.t.ledger)
        self.eta_log=iv(c,data['original_parameter_sources']['original_eta_log'])
        self.q=prior.ScaledEnclosure(prior.FormalScale(self.t.scale.bases,offset=self.eta_log/2-c.ln(2)/2),1,self.t.ledger)
        self.qrows={o:self.q if o==ZERO else self.t.scalar(0) for o in ORDERS}
        self.logC=iv(c,data['original_logC']);self.dstar=iv(c,data['original_dstar_log'])
        def factory(_,values):return Taylor(self.t,values)
        factory.variable=lambda _,value,order:Taylor(self.t,[value,1]+[0]*(order-1))
        self.phys,self.physproof=native.compile_original('pre_pulse_mixed_C4','physical_mixed',dict(
            IntervalTaylor=factory,square=lambda v:v*v,derivative=lambda v:v.derivative()))
        self.stress,self.stressproof=native.compile_original('current_pre_pulse_stress_operator','raw_pre_stress_rows',dict(axial_derivative=lambda v:v.derivative()))
        support=object.__new__(supported.cutoff.NativeCutoffQCover)
        self.kernel,self.kernelproof=supported.compile_supported_density(support)

    def source(self,box):
        c=self.c;t=c.mpf(box);old=self.initial;EY=self.EY
        d=c.exp(-t);d3=c.exp(-c.mpf('1.5')*t);root=c.exp(-t/2);zero=self.z*0
        E=EY*root
        hist=dict(m=old['m']*d,h=old['h']*d3+EY*(d3*c.expm1(t)),k=old['k']*d3,
            e=(old['e']-EY*EY*(t/2))*d,p=old['p']+EY*EY*(-c.expm1(-t)/2))
        class Context:
            def mpf(_,value):return native.ProtocolConstant(prior.FormalScale(self.t.scale.bases),c.mpf(value),self.t.ledger)
        physical=self.phys(Context(),self.Z,self.delta,E,[zero-c.mpf('.5')]+[zero]*3,[zero]*5,hist,self.P0,self.invP2)
        return dict(E=E,history=hist,physical=physical,record=dict(original_forward_buffer_t=box,
            original_EY=native.records(EY),original_velocity_amplitude=native.records(E),
            original_history_C0_Z_Taylor=native.records(hist),original_separate_P0=native.records(self.P0),
            original_physical_mixed4=native.records(physical),copied_original_physical_program=self.physproof,
            exact_zero_Uz_and_all_its_slow_jets=True,all_five_histories_follow_original_ODEs=True,
            original_background_not_pulse_incoming=True))

    def signed_roots(self,got,box):
        c=self.c;t=self.t;phys=got['physical']
        U=phys['physical_velocity_pressure_y_derivative_Taylor']['Utheta_over_Pstar']
        V=phys['physical_velocity_pressure_y_derivative_Taylor']['Uz'];hist=phys['actual_normalized_primitive_y_derivative_axial5']
        P=[self.P0+hist['p'][0]]+hist['p'][1:]
        parts=self.stress(c,self.delta,self.z,U,V,hist,P)
        Ps=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,2,0)),1,t.ledger)
        radius=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,22,0),10*self.logC+c.ln(110)-11+c.mpf(box)),1,t.ledger)
        numerator={}
        for label,omit in (('theta','variable_radial_shear'),('axial','axial_radial_shear')):
            rows=[U[j]*0 for j in range(3)]
            for name,part in parts[label].items():
                if name==omit:continue
                factor=1 if label=='theta' else self.invPs if part['mode'][1]==0 else Ps
                for j in range(3):rows[j]=rows[j]+part['shape'][j]*factor
            numerator[label]=[sum((rows[i]*math.comb(j,i) for i in range(j+1)),rows[0]*0)*radius for j in range(3)]
        roots={name:{(j,k):rows[j][k]*math.factorial(k) for j,k in ORDERS}
            for name,rows in (('E',U),('p1_numerator',numerator['theta']),('p2_numerator',numerator['axial']))}
        low=roots['E'][ZERO].scale.evaluate()+c.ln(c.mpf(ep(roots['E'][ZERO].coefficient)[0]))
        roots['p1']=current.quotient_jet(roots.pop('p1_numerator'),roots['E'],low)
        roots['p2']=current.quotient_jet(roots.pop('p2_numerator'),roots['E'],low)
        for name in ('a','b','t0','kappa_minus2'):
            roots[name]={o:t.scalar(2 if name=='a' and o==ZERO else 0) for o in ORDERS}
        return dict(roots=roots,q=self.q,record=dict(original_signed_slow_roots={name:{'y%d_Z%d'%o:v.record() for o,v in rows.items()} for name,rows in roots.items()},
            original_radius=radius.record(),original_radius_identity='R=Rd*exp(t-11)=110*Cstar^10*Pstar^11*exp(t-11)',
            original_a2_b0_t0zero_constant_q=True,original_nonzero_q=self.q.record(),
            copied_original_signed_stress_program=self.stressproof,original_P0_plus_p_added_once=True,
            common_Pstar_units_not_reapplied=True))


def absolute_cover(v):
    c=v.ctx;lo,hi=ep(v.coefficient)
    return prior.ScaledEnclosure(v.scale,c.mpf(max(abs(lo),abs(hi))),v.ledger)


def symmetric_cover(v):
    c=v.ctx;lo,hi=ep(v.coefficient);bound=max(abs(lo),abs(hi))
    return prior.ScaledEnclosure(v.scale,c.mpf((-bound,bound)),v.ledger)


def reflected_density_mean(E,primitives,N):
    c=E.ctx;A=primitives['A']*(c.mpf(1)/N);B=primitives['B_over_Pstar']*(c.mpf(1)/N)
    cap=c.mpf(5)/(4*N)
    # cosh(x)-1=x^2*integral_0^1(1-s)*cosh(s*x)ds. The exact
    # nonzero x^2 remains formal; only the analytic factor is bounded.
    cm1=lambda x,a:current.square(x)*c.mpf(('.5',ep(c.exp(a)/2)[1]))
    sinh=A*c.mpf((1,ep(c.exp(cap))[1]))
    p=current.square(E)*cm1(A*2,cap*2)*c.mpf('.5')
    return dict(m=E.scalar(0),h=E*cm1(A,cap),k=E*B*sinh,
        e=current.square(B)-p,p=p)


def phase_origin(adapter):
    row=adapter.data['actual_original_endpoint_source_queries'][0]['actual_original_phase']
    base=row['original_Rd_source_phase'];assert row['phase_independent_of_Z']
    return dict(original_saved_Rd_phase=base,original_buffer_phase_identity='phi(t)=Rd_phase+N*(t-11) mod1',
        exact_whole_buffer_cycles=11*N,phase_independent_of_Z=True,
        actual_phase_origin_not_selected=True,all_possible_origin_shifts_retained_in_weight_bounds=True)


def integrate_cell(adapter,left,right,bins=PHASE_BINS):
    c=adapter.c;coordinates=adapter.coordinates;box=c.mpf((left,right));width=c.mpf(right)-c.mpf(left)
    cycles=width*N
    if ep(cycles)[0]!=ep(cycles)[1] or ep(cycles)[0]%1 or ep(width)[0]<=0:
        raise ValueError('Exact positive integer-period original source cell required')
    got=adapter.source(box);source=adapter.signed_roots(got,box)
    u,branches,empty=supported.cover.signed_u_branches(source,adapter.dstar)
    phase_cells=[];densities=[];slow_rows=[];means=[];z_bins=[]
    for index in range(bins):
        phi=c.mpf((c.mpf(index)/bins,c.mpf(index+1)/bins));local=[];slow_local=[];mean_local=[]
        for branch in branches:
            inverse=supported.cover.branch_first_jets(source,adapter.qrows,adapter.dstar,phi,branch)
            if inverse['values'] is None:raise ArithmeticError('Buffer phase inverse requires refinement: '+str((left,right,index)))
            values=inverse['values'];E=source['roots']['E'][ZERO]
            EZ=source['roots']['E'][(0,1)];Ey=source['roots']['E'][(1,0)];zero=E.scalar(0)
            density=adapter.kernel(E,EZ,zero,zero,values,N)
            yvalues=dict(values,A_Z=values['A_y'],B_Z_over_Pstar=values['B_y_over_Pstar'])
            slow=adapter.kernel(E,Ey,zero,zero,yvalues,N)['Z_derivatives']
            mean=reflected_density_mean(E,values,N) if index<bins//2 else None
            local.append(density);slow_local.append(slow)
            if mean is not None:mean_local.append(mean)
            phase_cells.append(dict(original_phase_bin=index,actual_complete_period_phase_box=phi,
                original_signed_u_branch=branch['name'],original_conditional_first_jets=inverse['record'],
                actual_nonlinear_five_C0_density=native.records(density['kernels']),actual_nonlinear_five_Z_density=native.records(density['Z_derivatives']),
                actual_complete_fixed_phase_slow_y_density=native.records(slow),
                actual_reflected_C0_pair_density=None if mean is None else native.records(mean)))
        union=rc.density.local.same_source_union
        densities.append({k:union([d['kernels'][k] for d in local]) for k in KEYS})
        slow_rows.append({k:union([d[k] for d in slow_local]) for k in KEYS})
        z_bins.append({k:union([d['Z_derivatives'][k] for d in local]) for k in KEYS})
        if mean_local:means.append({k:union([d[k] for d in mean_local]) for k in KEYS})
    union=rc.density.local.same_source_union;zero=coordinates.scalar(0)
    whole={k:union([d[k] for d in densities]) for k in KEYS}
    whole_y={k:union([d[k] for d in slow_rows]) for k in KEYS}
    averaged={k:sum((d[k]*(c.mpf(2)/bins) for d in means),zero) for k in KEYS}
    geometry=dict(width=coordinates.scalar(width),regular=width,scalar_cover=width,
        record=dict(chart='O2_buffer',original_forward_t_endpoints=[str(left),str(right)],
            positive_true_log_radius_width=coordinates.scalar(width).record(),ordinary_log_radius_Jacobian=1,
            geometry_and_phase_origin_independent_of_Z=True,exact_original_complete_periods=cycles))
    factors={k:rc.transfer.true_width_kernel(coordinates,geometry,r) for k,r in rc.RATES.items()}
    values={};jets={};errors={};z_weights={}
    for key,rate in rc.RATES.items():
        r=Fraction(rate);rr=c.mpf(r.numerator)/r.denominator
        slow_error=absolute_cover(whole_y[key])*(c.mpf(1)/N)
        weight_error=absolute_cover(whole[key])*(2*c.expm1(rr/N))
        error=slow_error+weight_error;errors[key]=dict(full_slow_envelope_error=slow_error.record(),period_weight_error=weight_error.record())
        values[key]=factors[key]['mass']*(averaged[key]+symmetric_cover(error))
        # For any actual phase origin, each complete period spends exactly
        # 1/bins of its width in each bin. Duhamel weight varies by at most
        # exp(rate/N) within it, including wrapped bins. No origin selected.
        weight=c.exp(c.mpf((-ep(rr/N)[1],ep(rr/N)[1])))*(c.mpf(1)/bins)
        jets[key]=factors[key]['mass']*sum((d[key]*weight for d in z_bins),zero)
        z_weights[key]=dict(per_phase_bin_fraction_cover=weight,per_period_origin_uniform_ratio='exp(-rate/N) <= bin_mass/(period_mass/bins) <= exp(rate/N)')
    return dict(geometry=geometry,values=values,Z_derivatives=jets,record=dict(status='enclosed',original_source=got['record'],
        original_signed_stress_source=source['record'],actual_original_phase_binding=phase_origin(adapter),
        actual_source_unconditioned_signed_u=u.record(),empty_signed_u_branches=empty,
        actual_original_nonlinear_phase_density_queries=phase_cells,original_supported_density_AST_replacements=adapter.kernelproof,
        fixed_slow_reflection_mean_C0_cover=native.records(averaged),complete_fixed_phase_density_C0_cover=native.records(whole),
        complete_fixed_phase_slow_y_density_cover=native.records(whole_y),original_per_rate_period_and_slow_errors=errors,
        original_Z_phase_bin_weight_theorem=z_weights,
        actual_local_five_C0_integrals=native.records(values),actual_local_five_Z_integrals=native.records(jets),
        original_true_width_Duhamel_factors={k:dict(mass=v['mass'].record(),decay=v['decay'].record(),branch=v['branch']) for k,v in factors.items()},
        original_true_geometry=geometry['record'],nonlinear_evaluation_before_signed_branch_hull=True,
        C0_reflection_mean_is_frozen_slow_with_explicit_errors=True,
        Z_direct_density_phase_integrals_not_exact_reflection_cancellation=True,
        primitive_support_caps_not_used_as_field_values=True))


def propagate_to_Rc(coordinates,data,operator,transition):
    family=coordinates.family
    transport=data['original_whole_axial_transport_with_genuine_incoming']
    restore=lambda rec:complete.restore_half_source(rec,coordinates)
    incoming={k:restore(v) for k,v in transport['actual_axial_exit_correction_C0'].items()}
    incomingZ={k:restore(v) for k,v in transport['actual_axial_exit_correction_Z'].items()}
    atRd=operator.apply(incoming,incomingZ,family)
    final=rc.history.C1DuhamelOperator(coordinates);saved=transition['original_Rd_to_Rc_operator']
    final.coefficients={k:restore(v) for k,v in saved['incoming_C0_Z_decay_coefficients'].items()}
    final.increments={k:restore(v) for k,v in saved['cumulative_signed_increment_C0_enclosures'].items()}
    final.Z_increments={k:restore(v) for k,v in saved['cumulative_signed_increment_Z_enclosures'].items()}
    atRc=final.apply(atRd['values'],atRd['Z_derivatives'],family)
    old=rc.history.CommonSourceCoordinates(coordinates.ctx,coordinates.logP_squared,family)
    lift=lambda rec:coordinates.rebase(common.restore_common_source(rec,old),family)
    background=transition['original_Rc_background_and_separate_P0_Z']
    original={k:lift(v) for k,v in background['original_normalized_history_C0_enclosures'].items()}
    originalZ={k:lift(v) for k,v in background['original_normalized_history_Z_enclosures'].items()}
    P0=lift(background['original_separate_P0_over_Pstar_squared']);P0Z=lift(background['original_separate_P0_Z_over_Pstar_squared'])
    own={k:original[k]+atRc['values'][k] for k in KEYS};ownZ={k:originalZ[k]+atRc['Z_derivatives'][k] for k in KEYS}
    return dict(actual_axial_exit_incoming_C0=native.records(incoming),actual_axial_exit_incoming_Z=native.records(incomingZ),
        original_whole_buffer_operator=operator.record(),actual_updated_Rd_correction_C0=native.records(atRd['values']),actual_updated_Rd_correction_Z=native.records(atRd['Z_derivatives']),
        accepted_original_Rd_to_Rc_operator=saved,actual_updated_Rc_correction_C0=native.records(atRc['values']),actual_updated_Rc_correction_Z=native.records(atRc['Z_derivatives']),
        actual_updated_Rc_own_history_C0=native.records(own),actual_updated_Rc_own_history_Z=native.records(ownZ),
        actual_updated_Rc_absolute_pressure=(P0+own['p']).record(),actual_updated_Rc_absolute_pressure_Z=(P0Z+ownZ['p']).record(),
        original_Rc_background_and_separate_P0_Z=background,
        new_actual_axial_exit_not_historical_coarse_tail_incoming_used=True,
        accepted_transition_density_not_rerun_only_operator_applied=True,original_P0_added_once_and_rate0_pressure_memory_retained=True,
        actual_terminal_defects_controls_full_Z_axis_global_N_stress_recursion_admitted=False)


def execute_tile(coordinates,data,transition,cells=CELL_COUNT,bins=PHASE_BINS):
    adapter=OriginalBuffer(coordinates,data);operator=rc.history.C1DuhamelOperator(coordinates);rows=[];c=coordinates.ctx
    for index in range(cells):
        left=Fraction(11*index,cells);right=Fraction(11*(index+1),cells)
        cv=lambda v:c.mpf(v.numerator)/v.denominator
        local=integrate_cell(adapter,cv(left),cv(right),bins)
        rc.transfer.append_true_cell(operator,local['geometry'],local['values'],local['Z_derivatives'],coordinates.family)
        rows.append(dict(original_exact_t_endpoints=[str(left),str(right)],actual_period_source_and_integrals=local['record']))
        print('Original buffer period source:',data['exact_Z_range'],str(left),str(right),'enclosed',flush=True)
    operator.coefficients={k:coordinates.decay(11,r) for k,r in rc.RATES.items()}
    result=propagate_to_Rc(coordinates,data,operator,transition)
    return dict(source_family=coordinates.family,candidate_N=N,exact_Z_range=data['exact_Z_range'],
        original_buffer_parameter_and_phase_binding=phase_origin(adapter),ordinary_width_is_11_not_exp_M=True,
        source_cells=cells,full_period_phase_bins=bins,actual_original_buffer_source_rows=rows,
        common_directed_coordinate_theorem=coordinates.record(),complete_original_buffer_C0_Z_source_and_operator_executed=True,
        original_nonlinear_C0_period_reflection_with_slow_and_weight_errors=True,
        actual_updated_Rd_Rc_correction_and_own_history=result,
        full_original_target_controls_N_axis_heat_stress_recursion_admitted=False)


@rc.native.inlet.source_precision
def run():
    began=time.monotonic();accepted=json.loads((HERE/axial.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,axial,accepted['source_family'])
    transition,_=common.attach_receipt(hashes,complete,accepted['source_family'])
    c=MPIntervalContext();c.dps=240;archives=[]
    for archive in accepted['actual_original_axial_endpoint_native_source_archives']:
        raw=gzip.decompress((HERE/archive['filename']).read_bytes());assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
        data=json.loads(raw);ta=next(a for a in transition['actual_original_complete_transition_archives'] if a['exact_Z_range']==data['exact_Z_range'])
        traw=gzip.decompress((HERE/ta['filename']).read_bytes());assert hashlib.sha256(traw).hexdigest()==ta['lossless_original_json_sha256']
        coordinates=native.HalfPstarCoordinates(c,iv(c,data['common_directed_coordinate_theorem']['common_log_bases'][1]),accepted['source_family'])
        result=execute_tile(coordinates,data,json.loads(traw));result.update(accepted_original_axial_archive=archive,accepted_original_transition_archive=ta)
        encoded=json.dumps(encode(result),indent=2).encode()+b'\n';compressed=gzip.compress(encoded,compresslevel=9,mtime=0)
        tag='positive' if Fraction(data['exact_Z_range'][0])>0 else 'negative';name=PREFIX+'current_buffer_period_transport_'+tag+'.json.gz'
        (HERE/name).write_bytes(compressed);hashes[name]=sha(name)
        archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(encoded),lossless_original_json_sha256=hashlib.sha256(encoded).hexdigest(),exact_Z_range=data['exact_Z_range']))
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=accepted['source_family'],candidate_N=N,actual_original_buffer_period_transport_archives=archives,
        actual_original_source_cells=2*CELL_COUNT,actual_phase_bins=PHASE_BINS,exact_periods_per_buffer=11*N,
        original_whole_buffer_C0_Z_integrals_and_new_Rd_Rc_histories_executed=True,
        original_C0_nonlinear_period_reflection_with_explicit_slow_weight_errors=True,
        Z_period_bin_weight_covers_not_exact_cancellation=True,
        original_accepted_upstream_or_transition_numerical_producers_not_rerun=True,
        actual_functional_targets_controls_full_Z_axis_N_heat_stress_recursion_admitted=False,
        **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Entire original 11-unit buffer on two strict-sign tiles at N1024, actual nonlinear half-period reflection C0 averages with full slow/weight errors and direct nonlinear Z phase-bin integrals, genuine new axial exit to Rd then accepted full transition/power to updated Rc correction/own-history/pressure. Functional terminal defects/controls/full Z/axis/global N/heat/stress/recursion remain open.')
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

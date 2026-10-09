"""Actual O3 slope-mu source and real current finite-N prefix through Rc.

Positive mu/Delta remain formal source parameters. Directed scalar covers
only enclose their amplitude/kernel arithmetic; no parameter is selected
from a cover. Quiet power follows from the same-source eta comparison.
"""
import ast
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N as previous
import lei_ren_part1_paper_compliant_outer_buffer as kernels

fields,base,ep=previous.fields,previous.base,previous.ep
incoming,reference,switch=previous.incoming,previous.reference,previous.switch
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
RATES,C0,Z=previous.RATES,previous.C0,previous.Z
PARTITIONS={'transition':((0,1),(1,4),(1,2),(3,4),(1,1)),
            'power':((0,1),(1,1),(2,1))}
SCALAR_CELLS=128
NAME=PREFIX+'current_original_O3_Rc_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_O3_Rc_finite_N_check.json'
GATE='current_original_source_bound_finite_N_O3_prefix_through_Rc_enclosed'


def fraction(chart,value):
    if chart not in PARTITIONS or not isinstance(value,tuple) or len(value)!=2 or any(type(v) is not int for v in value) or value[1]<=0:
        raise ValueError('Exact original O3 local log-radius coordinate required')
    q=Fraction(*value)
    if not 0<=q<=(1 if chart=='transition' else 2):raise ValueError('O3 transition[0,1] or Rw..Rc power[0,2] required')
    return q


def nonnegative_Delta_q(a,Delta,eta_log,log_a_lower):
    """Original q, sharpened by a>=2,Delta>=0 on its active subdomain."""
    if ep(Delta.coefficient)[0]<0:raise ValueError('Actual nonnegative O3 Delta source required')
    original=switch.parameters.original
    result=original.q_enclosure(a,Delta,eta_log,log_a_lower)
    if result['branch']=='flat' or Delta.zero:return result
    # On Delta<eta, 0<=q^2<=eta/a<=eta/2; the flat subdomain is0.
    # This source theorem avoids using a huge independent gamma cover.
    eta=a.scalar(1);eta=type(eta)(type(eta.scale)(eta.scale.bases,offset=eta_log),1,a.ledger)
    cap=original.nonnegative_sqrt(eta*mp.mpf('.5'))
    result['q']=type(cap)(cap.scale,a.ctx.mpf([0,ep(cap.coefficient)[1]]),a.ledger)
    result['original_nonnegative_Delta_q_cap_intersection']=dict(
        exact_source_a_ge2=True,exact_source_Delta_ge0=True,
        active_source_q_squared_upper='eta/2',flat_source_exact_q_zero=True,
        positive_mu_not_replaced_by_arithmetic_cover=True)
    return result


def compile_correlated_O3_quotients():
    tree=ast.parse(Path(switch.__file__).read_text(encoding='utf8'))
    fn=copy.deepcopy(next(n for n in tree.body if getattr(n,'name',None)=='general_quotients'))
    wanted=ast.dump(ast.parse('Delta=parameters.add(f,kappa,[f.scalar(-2)]+[f.scalar(0)]*5)').body[0])
    sites=[n for n in ast.walk(fn) if isinstance(n,ast.Assign) and ast.dump(n)==wanted]
    if len(sites)!=1:raise ValueError('Original generic shear excess source changed')
    sites[0].value=ast.Name(id='source_Delta',ctx=ast.Load())
    calls=[n for n in ast.walk(fn) if isinstance(n,ast.Call) and ast.unparse(n.func)=='parameters.original.q_enclosure']
    if len(calls)!=1:raise ValueError('Original q callback changed')
    calls[0].func=ast.Name(id='nonnegative_Delta_q',ctx=ast.Load())
    fn.args.args.append(ast.arg(arg='source_Delta'));fn.name='O3_quotients'
    env=dict(vars(switch),nonnegative_Delta_q=nonnegative_Delta_q)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original O3 collected positive excess>','exec'),env)
    return env[fn.name]


O3_QUOTIENTS=compile_correlated_O3_quotients()


def source_bindings():
    original=previous.previous.original
    tree=ast.parse(Path(original.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if getattr(n,'name',None)=='CompliantPrePulseMixedC4')
    recipes={
        'slope_mu':("factor=c.exp(-t/2-mu*K['J'])",
            "hist=dict(m=old['m']*d,h=old['h']*d3+u1*(d3*K['theta']),k=old['k']*d3,e=old['e']*d-square(u1)*(d*K['energy']/2),p=old['p']+square(u1)*(K['pressure']/2))"),
        'power':("f=c.exp((-c.mpf('.5')-mu)*t)",
            "hist=dict(m=old['m']*d,h=old['h']*d3+u1*theta,k=old['k']*d3,e=old['e']*d-square(u1)*(d*decay_integral(c,2*mu,t)/2),p=old['p']+square(u1)*(decay_integral(c,1+2*mu,t)/2))")}
    for name,statements in recipes.items():
        fn=next(n for n in cls.body if getattr(n,'name',None)==name)
        for statement in statements:
            wanted=ast.dump(ast.parse(statement).body[0])
            if sum(ast.dump(n)==wanted for n in ast.walk(fn) if isinstance(n,ast.Assign))!=1:
                raise ValueError('Original O3 source changed: '+name+'.'+statement)
    return dict(original_O3_transition_power_assignments=recipes,
        exact_log_mu='log(.001)-4*(exp(40)+11)',exact_Tw='-60*log_mu',
        actual_Rd_background_and_finite_N_incoming_supplied_separately=True,
        exact_collected_Delta_transition='2*mu*sigma(t)',exact_collected_Delta_power='2*mu',
        whole_transition_kernels='positive endpoint hulls of original same-mu integral functions',
        scalar_mu_cover_is_not_source_mu=True,exact_positive_formal_mu_before_arithmetic_with2=True,
        original_nonnegative_Delta_active_q_bound='q^2<=eta/2 from a>=2,Delta>=0',
        exact_radius='Rd*exp(t); power Rw*exp(t), Rc=Rw*exp2',
        actual_power_offsets_not_phase_or_Tw_widths=True,
        original_full_signed_recovery=reference.long.generic.source_bindings())


def transition_kernel_hull(c,t,mu_cover,cells=SCALAR_CELLS):
    """Positive integrals are monotone in endpoint t for the same exact mu.

    Evaluate the pure directed helper at fixed endpoint boxes, rather than
    losing the common t in a=t*i/n,b=t*(i+1)/n on a wide interval.
    """
    lo,hi=ep(t)
    if lo<0 or hi>1 or type(cells) is not int or not 16<=cells<=4096:
        raise ValueError('Original O3 transition endpoint interval and bounded cell count required')
    left=kernels.transition_kernels(c,c.mpf(lo),mu_cover,cells)
    right=left if lo==hi else kernels.transition_kernels(c,c.mpf(hi),mu_cover,cells)
    result={key:c.mpf([max(mp.mpf(0),ep(left[key])[0]),ep(right[key])[1]]) for key in ('J','theta','energy','pressure')}
    result.update(cells=cells,positive_same_source_endpoint_monotonicity_used=True,
        actual_mu_enclosure_not_parameter_selection=True)
    return result


def positive_decay_mass(c,k,t):
    """int_0^t exp(-k*s)ds, allowing a zero lower bound for tiny k."""
    if ep(k)[0]<0 or ep(t)[0]<0:raise ValueError('Nonnegative decay source and offset required')
    if ep(t)==(0,0):return c.mpf(0)
    if ep(k)[0]>0 and ep(k*t)[0]>mp.mpf('.001'):
        mass=-c.expm1(-k*t)/k
        return c.mpf([max(mp.mpf(0),ep(mass)[0]),ep(mass)[1]])
    return c.mpf([max(mp.mpf(0),ep(t*c.exp(-k*t))[0]),ep(t)[1]])


def background_cell(op,parent,chart,left,right,logmu):
    lo,hi=fraction(chart,left),fraction(chart,right)
    if hi<lo:raise ValueError('Ordered actual O3 source cell required')
    f,c=op.flow,op.c;cv=lambda q:c.mpf(q.numerator)/q.denominator
    t=c.mpf([ep(cv(lo))[0],ep(cv(hi))[1]])
    key='original_closed_O2_axial_buffer_background' if chart=='transition' else 'original_closed_O3_background'
    raw=parent[key]['raw']['raw_current_radius_y_derivative_axial_coefficients']
    if parent['exact_common_P0_axial5'] is not op.P0:raise ValueError('Same live original O3 parent P0 required')
    u1=raw['velocity']['theta'][0];old={name:rows[0] for name,rows in raw['histories'].items()}
    mu=f.factor((0,0,0,0,0),logmu);mu_cover=f.ordinary_cover(mu)
    if ep(mu_cover)[0]<0 or ep(mu_cover)[1]>=mp.mpf('.001'):raise ValueError('Valid original small positive mu enclosure required')
    zero=[f.scalar(0)]*6;d=c.exp(-t);d3=c.exp(-c.mpf('1.5')*t)
    if chart=='transition':
        K=transition_kernel_hull(c,t,mu_cover);sig=previous.previous.original.sigma_jets(c,t)[0]
        factor=c.exp(-t/2-mu_cover*K['J']);theta=d3*K['theta'];energy=d*K['energy'];pressure=K['pressure']
    else:
        sig=c.mpf(1);factor=c.exp(-t/2-mu_cover*t)
        # Original (factor-d3)/(1-mu), written as a positive exprel.
        average=incoming.weighted.terminal.local.packets.recovery.exp_average(c,(1-mu_cover)*t)
        theta=d3*t*average;energy=d*positive_decay_mass(c,2*mu_cover,t)
        pressure=positive_decay_mass(c,1+2*mu_cover,t)
        K=dict(theta=theta,energy=energy,pressure=pressure,
            original_positive_exprel_theta_identity=True,actual_original_power_offset=t)
    E=f.scale(u1,factor);E2=f.multiply(u1,u1);V,Vy=zero,zero
    excess=mu*(2*sig);Delta=[excess]+zero[1:];a=[f.scalar(2)+excess]+zero[1:]
    Ey=f.add(f.scale(E,-c.mpf('.5')),f.scale(E,-mu*sig))
    hist=dict(m=f.scale(old['m'],d),h=f.add(f.scale(old['h'],d3),f.scale(u1,theta)),
        k=f.scale(old['k'],d3),e=f.add(f.scale(old['e'],d),f.scale(E2,-energy/2)),
        p=f.add(old['p'],f.scale(E2,pressure/2)))
    dy=dict(m=f.scale(hist['m'],-1),h=f.add(E,f.scale(hist['h'],-c.mpf('1.5'))),
        k=f.scale(hist['k'],-c.mpf('1.5')),
        e=f.add(f.scale(hist['e'],-1),f.scale(f.multiply(E,E),-c.mpf('.5'))),
        p=f.scale(f.multiply(E,E),c.mpf('.5')))
    radius=op.Rm_factor*f.factor((0,0,0,0,0),c.exp(40)+(17 if chart=='transition' else 18)+t)
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=dict(histories={name:[row,dy[name]] for name,row in hist.items()},
        velocity=dict(theta=[E,Ey],axial=[V,Vy])),original_P0_normalized_axial5=op.P0,
        geometry=dict(chart='O3_'+chart,exact_source_coordinate_cell=(left,right),actual_physical_radius=radius,
            actual_local_log_radius_coordinate=t,phase_and_radius_Z_independent=True))
    return dict(raw=packet,actual_a_axial5=a,actual_collected_Delta_axial5=Delta,
        original_positive_formal_mu=mu,actual_mu_arithmetic_enclosure_only=mu_cover,
        original_transition_or_power_kernels=K,original_sigma=sig,
        exact_source_Rm_factor=op.Rm_factor,actual_source_radius=radius,actual_log_mu=logmu,
        exact_source_parameter_not_selected_from_cover=True,
        original_all_five_incoming_background_histories_retained=True,
        original_field_and_histories_are_closed_source_functions=True)


def guard_incoming(family,label,N,op,packet):
    switch.guard_saved_source(family,label,N,op.P0,packet)
    if packet.get('genuine_current_finite_N_prefix_through_Rd_boundary_enclosures_supplied') is not True:
        raise ValueError('Genuine current finite-N Rd input required')
    radius=switch.canonical_expression([op.Rm_factor*op.flow.factor((0,0,0,0,0),17+op.c.exp(40))])[0]
    if incoming.canonical_saved(packet['exact_source_Rd_factor'])!=radius:
        raise ValueError('Same actual canonical Rd radius required')
    incoming.row_packet(packet,'actual_current_Rd_correction_C0_Z')


class OriginalO3RcFiniteN:
    mode='current_original_O3_transition_quiet_power_and_actual_finite_N_Rd_to_Rc_prefix'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalO2AxialBufferFiniteN(dps)
        self.c=self.upstream.c;self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes)
        self.bindings=source_bindings();self.cache={};self.cells={};self.parents={}
        self.loop_parameters=self.upstream.upstream.upstream.upstream.upstream.parameters
        self.logmu=self.c.ln(self.c.mpf('.001'))-4*(self.c.exp(40)+11)
        self.Tw=-60*self.logmu
        if ep(self.logmu)[1]>=-1000:raise ValueError('Original positive mu must admit its directed tiny arithmetic cover')
        self.power_margin=self.logmu+self.c.ln(2)-self.loop_parameters.eta_log
        if ep(self.power_margin)[0]<0:raise ValueError('Original O3 2mu>=eta quiet power not proved')
        if ep(self.Tw)[0]<=2+ep(self.c.ln(2))[1]:raise ValueError('Actual Rc/repair reservation lies outside original power chart')
        self.saved=json.loads(gzip.decompress((HERE/previous.NAME).read_bytes()))
        if not self.saved[previous.GATE] or self.saved['source_family']!=self.family or self.saved['candidate_N']!=257:
            raise ValueError('Accepted same-current actual Rd prefix required')
        for name,digest in self.saved['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
        for module in (kernels,previous.previous.original,switch):fields.previous.bind(self.hashes,Path(module.__file__).name,sha(Path(module.__file__).name))
        for name in (previous.NAME,Path(__file__).name):fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or not checked.get(GATE) or checked['source_family']!=self.family:
                raise ValueError('Accepted actual O3 Rc finite-N receipt required')
            for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def owner(self,label):return self.upstream.owner(label)

    def query(self,label,chart,left,right,N=257):
        N=reference.phase.candidate_N(N)
        if N!=257 or label not in ('0','.5'):raise ValueError('Same current frames0,.5,N257 required')
        lo,hi=fraction(chart,left),fraction(chart,right);key=(label,chart,lo,hi,N)
        if key in self.cells:return self.cells[key]
        if hi<lo:raise ValueError('Ordered actual O3 source cell required')
        op=self.owner(label);f,c=op.flow,op.c
        with mp.workdps(c.dps+40):
            if label not in self.parents:self.parents[label]=self.upstream.query(label,'buffer',(11,1),(11,1))
            parent=self.parents[label] if chart=='transition' else self.query(label,'transition',(1,1),(1,1),N)
            source=background_cell(op,parent,chart,left,right,self.logmu);a=source['actual_a_axial5']
            E=source['raw']['raw_current_radius_y_derivative_axial_coefficients']['velocity']['theta'][0]
            proxy=SimpleNamespace(flow=f,c=c,reference=op.reference,zrows=op.zrows,P0=op.P0,
                Pstar=op.Pstar,source_radius=source['actual_source_radius'],correlated_C=f.multiply(a,E))
            generic=reference.long.RECOVER(proxy,source['raw']);generic.pop('source_frame_conditional_on_same_actual_Rm_inlet')
            generic['source_from_same_current_original_O3_background']=True
            roots,qr,quotients=O3_QUOTIENTS(proxy,generic,a,self.loop_parameters.eta_log,source['actual_collected_Delta_axial5'])
            got=reference.primitives.all_u_primitive_bounds(f,roots,qr,self.loop_parameters.dstar_log,c.mpf([0,1]))
            if qr[C0].zero and qr[Z].zero:
                got['values']={name:f.scalar(0) for name in got['values']}
                got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
            if chart=='power' and (not qr[C0].zero or not qr[Z].zero):raise ValueError('Original admitted power source must be exactly flat')
            V=generic['common_velocity_V_axial5']
            density=reference.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],N)
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=op.P0,
                original_closed_O3_background=source,original_generic_source=generic,
                original_full_signed_quotients=quotients,original_roots=roots['roots'],original_q_rows=qr,
                original_primitive_values=got['values'],original_primitive_proof=got['record'],
                original_signed_five_density_C0_Z=density,exact_source_radius=source['actual_source_radius'],
                actual_compact_source_coordinate_cell=(left,right),actual_chart=chart,
                actual_phase_full_period_cover=c.mpf([0,1]),actual_phase_Z_exact_zero=True,
                actual_global_radius_phase='frac(N*(logRd+t-logRa-hb*s_c/2))' if chart=='transition' else 'frac(N*(logRw+t-logRa-hb*s_c/2))',
                actual_physical_log_measure='dy=dt; local log-radius offset, not power phase',
                actual_same_source_power_flat_log_margin=self.power_margin if chart=='power' else None,
                positive_formal_mu_and_correlated_Delta_before_critical_addition=True,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cells[key]=result;return result

    def seam(self,label,join):
        self.query(label,'transition',(0,1),(0,1))
        left=self.parents[label] if join=='Rd' else self.query(label,'transition',(1,1),(1,1))
        right=self.query(label,'transition',(0,1),(0,1)) if join=='Rd' else self.query(label,'power',(0,1),(0,1))
        left,right=left['original_generic_source'],right['original_generic_source'];rows=0
        pairs=[(left['common_velocity_E_axial5'],right['common_velocity_E_axial5']),
            (left['common_velocity_V_axial5'],right['common_velocity_V_axial5'])]
        pairs.extend((left['common_own_five_histories_axial5'][k],right['common_own_five_histories_axial5'][k]) for k in RATES)
        for a,b in pairs:
            for x,y in zip(a,b):
                d=x-y
                if not d.zero and not ep(d.coefficient)[0]<=0<=ep(d.coefficient)[1]:raise ValueError('Original common O3 seam source covers must meet')
                rows+=1
        return dict(exact_function_join_by_zero_local_integral_and_decay_one=True,
            source_overlap_consistency_rows=rows,overlap_not_function_identity_proof=True,
            same_original_P0_and_five_history_seed=True)

    def contribution(self,label,N=257):
        N=reference.phase.candidate_N(N)
        if N!=257 or label not in ('0','.5'):raise ValueError('Same supplied current frame/frequency required')
        if (label,N) in self.cache:return self.cache[(label,N)]
        op=self.owner(label);f,c=op.flow,op.c;accepted=self.saved['frames'][label]
        with mp.workdps(c.dps+40):
            guard_incoming(self.family,label,N,op,accepted)
            restore=lambda row:switch.first.endpoint.restore_row(f,row)
            initial={name:[restore(row) for row in values] for name,values in incoming.row_packet(accepted,'actual_current_Rd_correction_C0_Z').items()}
            joins={name:self.seam(label,name) for name in ('Rd','Rw')};current=initial;charts={}
            for chart,partition in PARTITIONS.items():
                local={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[];end=1 if chart=='transition' else 2
                for left,right in zip(partition,partition[1:]):
                    source=self.query(label,chart,left,right,N);density=source['original_signed_five_density_C0_Z'];drivers={};weights={}
                    l,r=fraction(chart,left),fraction(chart,right);cv=lambda q:c.mpf(q.numerator)/q.denominator
                    width=cv(r-l);suffix=cv(end-r)
                    for name,rate in RATES.items():
                        mass=incoming.weighted.terminal.local.positive_kernel_mass(c,width,rate)
                        decay=f.scalar(1) if not rate else f.factor((0,0,0,0,0),-suffix*c.mpf(rate.numerator)/rate.denominator)
                        pair=[reference.bounds.symmetric(f,reference.bounds.magnitude(f,density[part][name])*mass*decay) for part in ('kernels','Z_derivatives')]
                        for n,row in enumerate(pair):local[name][n]+=row
                        drivers[name]=pair;weights[name]=dict(actual_own_rate_mass=mass,actual_suffix_to_chart_exit=decay)
                    cells.append(dict(source=source,actual_cell_to_chart_exit_driver_C0_Z=drivers,
                        actual_original_log_width=width,actual_original_suffix=suffix,original_own_rate_weights=weights))
                    print('Actual current O3 finite-N source cell',label,chart,left,right,flush=True)
                output=previous.formal_affine_transport(f,current,local,c.mpf(end))
                charts[chart]=dict(actual_current_chart_incoming_C0_Z=current,actual_chart_local_driver_C0_Z=local,
                    actual_current_chart_exit_C0_Z=output,actual_source_cells=cells,exact_source_partition=partition,
                    full_actual_logarithmic_width=c.mpf(end),pressure_memory_exactly_one=True)
                current=output
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=op.P0,
                exact_source_Rd_factor=op.Rm_factor*f.factor((0,0,0,0,0),17+c.exp(40)),
                exact_source_Rw_factor=op.Rm_factor*f.factor((0,0,0,0,0),18+c.exp(40)),
                exact_source_Rc_factor=op.Rm_factor*f.factor((0,0,0,0,0),20+c.exp(40)),
                actual_current_Rd_incoming_correction_C0_Z=initial,actual_O3_charts=charts,
                actual_current_Rc_correction_C0_Z=current,actual_Rc_background_source=self.query(label,'power',(2,1),(2,1)),
                typed_actual_O3_source_joins=joins,exact_source_log_mu=self.logmu,exact_source_Tw=self.Tw,
                original_same_source_power_flat_log_margin=self.power_margin,
                genuine_current_finite_N_prefix_through_Rc_boundary_enclosures_supplied=True,
                original_O3_transition_and_power_background_histories_retained=True,
                quiet_power_does_not_reset_real_current_incoming=True,
                original_P0_background_and_actual_correction_kept_separate=True,
                no_old_O3_owner_N1024_tail_or_source_replay=True,
                source_enclosures_not_chosen_field_values_or_functional_closure=True,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cache[(label,N)]=result;return result


def run():
    began=time.monotonic();owner=OriginalO3RcFiniteN(require_checked=False)
    frames={label:owner.contribution(label) for label in ('0','.5')}
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_source_bindings=owner.bindings,frames=reference.serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Actual current finite-N prefix reaches original Rc through slope-mu and proved quiet power',flush=True);return report


if __name__=='__main__':run()

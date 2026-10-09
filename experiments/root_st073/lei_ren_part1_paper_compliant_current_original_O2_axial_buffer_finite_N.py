"""Original O2 axial turnoff/buffer and real same-current N257 incoming.

Large exponential factors remain formal. Whole source cells are deliberately
conservative; their bounds are not a selected field or terminal closure.
Only pure cutoff/kernel helpers are reused, never old source owners.
"""
import ast
import copy
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_slope_finite_N as previous
import lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame as parameter_source
import lei_ren_part1_paper_compliant_outer_initial as kernels

fields,base,ep=previous.fields,previous.base,previous.ep
incoming,reference,switch=previous.incoming,previous.reference,previous.switch
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
RATES,C0,Z=previous.RATES,previous.C0,previous.Z
PARTITIONS={'axial':((0,1),(1,4),(1,2),(3,4),(1,1)),
            'buffer':((0,1),(1,1),(5,1),(11,1))}
KERNEL_CELLS=128
KERNEL_WINDOW=800
MD=40
NAME=PREFIX+'current_original_O2_axial_buffer_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_O2_axial_buffer_finite_N_check.json'
GATE='current_original_source_bound_finite_N_O2_axial_buffer_prefix_through_Rd_enclosed'


def fraction(chart,value):
    if chart not in PARTITIONS or not isinstance(value,tuple) or len(value)!=2 or any(type(v) is not int for v in value) or value[1]<=0:
        raise ValueError('Exact rational original axial phase or buffer offset required')
    q=Fraction(*value)
    if not 0<=q<=(1 if chart=='axial' else 11):
        raise ValueError('Original axial phase[0,1] or buffer offset[0,11] required')
    return q


def compile_correlated_critical_quotients():
    """Unchanged signed recovery, collecting kappa-2=b^2/2 before addition.

    Here a is the exact source constant2. Adding a tiny b^2/2 to2 and
    subtracting2 numerically would erase the defining positive shear excess.
    """
    tree=ast.parse(Path(switch.__file__).read_text(encoding='utf8'))
    fn=copy.deepcopy(next(n for n in tree.body if getattr(n,'name',None)=='general_quotients'))
    wanted=ast.dump(ast.parse('Delta=parameters.add(f,kappa,[f.scalar(-2)]+[f.scalar(0)]*5)').body[0])
    sites=[n for n in ast.walk(fn) if isinstance(n,ast.Assign) and ast.dump(n)==wanted]
    if len(sites)!=1:raise ValueError('Original generic critical-shear subtraction changed')
    sites[0].value=ast.parse('parameters.quotient(f,square,a,pa)',mode='eval').body
    fn.name='critical_quotients';env=dict(vars(switch))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original critical a2 collected excess>','exec'),env)
    return env[fn.name]


CRITICAL_QUOTIENTS=compile_correlated_critical_quotients()


def source_bindings():
    tree=ast.parse(Path(previous.original.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if getattr(n,'name',None)=='CompliantPrePulseMixedC4')
    fn=next(n for n in cls.body if getattr(n,'name',None)=='axial')
    assignments=('y=c.exp(md*phase)','B=turnoff_derivatives(c,y,md,phase)',
        'y=c.exp(md)+selector',"B=[c.mpf(0)]*5",'t=y-1',
        "hist=dict(m=old['m']*d+z*(4*K['B_mass']),h=old['h']*d3+u1*(root-d3),k=old['k']*d3+u1*z*(4*root*K['B_mass']),e=old['e']*d+square(z)*(16*self.invP2*K['B_squared_mass'])-square(u1)*(t*d/2),p=old['p']+square(u1)*((1-d)/2))")
    for statement in assignments:
        wanted=ast.dump(ast.parse(statement).body[0])
        if sum(ast.dump(n)==wanted for n in ast.walk(fn) if isinstance(n,ast.Assign))!=1:
            raise ValueError('Original axial/buffer source changed: '+statement)
    return dict(original_axial_buffer_assignments=assignments,actual_Md=MD,
        exact_parameter_logPstar='exp(40)+11',actual_slope_exit_background_and_correction_supplied_separately=True,
        exact_a_equals2_and_E_y_equals_minus_E_over2=True,
        exact_collected_Delta='b^2/2; no rounded (2+b^2/2)-2',
        unchanged_original_generic_quotients_except_algebraically_collected_Delta=True,
        buffer_kernel_identity='K_j(exp(Md)+tau)=exp(-tau)*K_j(exp(Md))',
        source_decays_and_radius_remain_positive_formal_exponentials=True,
        exact_global_phase='frac(N*(logRref+y-logRa-hb*s_c/2))',
        actual_log_measure='dy=Md*exp(Md*phase)*dphase; buffer dy=dtau',
        physical_width='exp(Md*right)-exp(Md*left), never right-left phase width',
        original_all_five_signed_generic_recovery=reference.long.generic.source_bindings())


def parameter_binding(family,hashes):
    """Hydrate checked exact definitions without constructing a point owner."""
    record=json.loads((HERE/parameter_source.RECEIPT).read_bytes())
    if not record.get('all_passed') or not record.get(parameter_source.GATE) or record['source_family']['actual_five_defect_family_sha256']!=family or not record['source_definition_binding']['passed']:
        raise ValueError('Checked exact same-current O2 parameter definitions required')
    for name,digest in record['input_hashes'].items():fields.previous.bind(hashes,name,digest)
    for name in (parameter_source.RECEIPT,Path(parameter_source.__file__).name):
        fields.previous.bind(hashes,name,sha(name))
    source=json.loads((HERE/parameter_source.PRESSURE_SOURCE).read_bytes())['compliant_source']
    definition=source['implicit_source_definition']
    canonical=lambda obj:hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()
    if canonical(definition)!=record['source_family']['implicit_source_sha256'] or source['datum_enclosure_sha256']!=record['source_family']['datum_enclosure_sha256']:
        raise ValueError('Exact parameter definition and pressure datum family differ')
    if definition['Md']!='40' or definition['logPstar']!='exp(Md)+11' or definition['c_mu']!='.001':
        raise ValueError('Original current Md40/logPstar/mu definitions required')
    return dict(source_family=family,exact_parameter_source_family=record['source_family'],actual_Md=40,exact_logPstar='exp(40)+11',
        selected_physical_family_sha256=record['selected_physical_family_sha256'],
        same_original_implicit_parameter_definition=True,no_old_parameter_or_point_owner_constructed=True)


def coordinate(c,chart,left,right):
    lo,hi=fraction(chart,left),fraction(chart,right)
    if hi<lo:raise ValueError('Ordered actual source cell required')
    l=c.mpf(lo.numerator)/lo.denominator;r=c.mpf(hi.numerator)/hi.denominator
    selector=c.mpf([ep(l)[0],ep(r)[1]])
    y=c.exp(MD*selector) if chart=='axial' else c.exp(MD)+selector
    return selector,y


def background_cell(op,parent,chart,left,right):
    f,c=op.flow,op.c;selector,y=coordinate(c,chart,left,right)
    raw=parent['original_closed_O2_slope_background']['raw']['raw_current_radius_y_derivative_axial_coefficients']
    u1=raw['velocity']['theta'][0];old={key:rows[0] for key,rows in raw['histories'].items()}
    if parent['exact_common_P0_axial5'] is not op.P0:raise ValueError('Same live slope-exit P0 required')
    B=previous.original.turnoff_derivatives(c,y,MD,selector) if chart=='axial' else [c.mpf(0)]*5
    if chart=='axial':K=kernels.turnoff_kernels(c,y,MD,KERNEL_CELLS,KERNEL_WINDOW)
    else:
        # The original B vanishes past exp(Md). Keep its incoming kernel
        # tail; the eleven-unit buffer multiplies it by exp(-tau).
        seed=kernels.turnoff_kernels(c,c.exp(MD),MD,KERNEL_CELLS,KERNEL_WINDOW)
        K={key:seed[key]*c.exp(-selector) for key in ('B_mass','B_squared_mass','retained_far_tail')}
        K.update(cells=KERNEL_CELLS,window=KERNEL_WINDOW,buffer_exact_retained_kernel_memory=True)
    t=y-1;formal=lambda log:f.factor((0,0,0,0,0),log)
    d,root,d3=formal(-t),formal(-t/2),formal(-c.mpf('1.5')*t)
    E=f.scale(u1,root);E2=f.multiply(u1,u1);z=op.zrows;z2=f.multiply(z,z)
    V=f.scale(z,4*B[0]);Vy=f.scale(z,4*B[1]);zero=[f.scalar(0)]*6
    hist=dict(m=f.add(f.scale(old['m'],d),f.scale(z,4*K['B_mass'])),
        h=f.add(f.scale(old['h'],d3),f.scale(u1,root-d3)),
        k=f.add(f.scale(old['k'],d3),f.scale(f.multiply(u1,z),root*(4*K['B_mass']))),
        e=f.add(f.scale(old['e'],d),f.scale(z2,f.factor((0,-1,0,0,0))*(16*K['B_squared_mass'])),f.scale(E2,d*(-t/2))),
        p=f.add(old['p'],f.scale(E2,(f.scalar(1)-d)*c.mpf('.5'))))
    Ey=f.scale(E,-c.mpf('.5'));a=[f.scalar(2)]+zero[1:]
    dy=dict(m=f.add(V,f.scale(hist['m'],-1)),h=f.add(E,f.scale(hist['h'],-c.mpf('1.5'))),
        k=f.add(f.multiply(E,V),f.scale(hist['k'],-c.mpf('1.5'))),
        e=f.add(f.scale(f.multiply(V,V),f.factor((0,-1,0,0,0))),f.scale(hist['e'],-1),f.scale(f.multiply(E,E),-c.mpf('.5'))),
        p=f.scale(f.multiply(E,E),c.mpf('.5')))
    radius=op.Rm_factor*formal(6+y)
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=dict(histories={name:[row,dy[name]] for name,row in hist.items()},
        velocity=dict(theta=[E,Ey],axial=[V,Vy])),original_P0_normalized_axial5=op.P0,
        geometry=dict(chart='O2_'+chart,exact_source_coordinate_cell=(left,right),actual_physical_radius=radius,
            actual_compact_source_coordinate=y,actual_source_selector=selector,phase_and_radius_Z_independent=True))
    return dict(raw=packet,actual_a_axial5=a,original_turnoff_kernels=K,original_cutoff_ordinary_y_derivatives=B,
        exact_source_Rm_factor=op.Rm_factor,actual_source_radius=radius,actual_compact_y=y,
        original_slope_exit_seed_is_live_current_source=True,positive_formal_history_decays=dict(mean=d,theta=d3,amplitude=root),
        original_field_and_histories_are_closed_source_functions=True,
        buffer_kernel_is_true_turnoff_exit_memory_not_zero=True)


def guard_incoming(family,label,N,op,packet):
    switch.guard_saved_source(family,label,N,op.P0,packet)
    if packet.get('genuine_current_finite_N_prefix_through_O2_slope_exit_supplied') is not True:
        raise ValueError('Genuine current finite-N O2 slope-exit input required')
    radius=switch.canonical_expression([op.Rm_factor*op.c.exp(7)])[0]
    if incoming.canonical_saved(packet['exact_source_O2_slope_exit_factor'])!=radius:
        raise ValueError('Same actual canonical O2 slope-exit radius required')
    incoming.row_packet(packet,'actual_current_O2_slope_exit_correction_C0_Z')


def physical_weights(f,chart,left,right,rate):
    """Exact endpoints/Jacobian, stable positive mass and formal suffix."""
    c=f.c;l,r=fraction(chart,left),fraction(chart,right)
    if r<=l:raise ValueError('Positive actual source interval required')
    cv=lambda q:c.mpf(q.numerator)/q.denominator
    if chart=='axial':
        yl,yr=c.exp(MD*cv(l)),c.exp(MD*cv(r))
        # Collect the width before enclosures, avoiding endpoint subtraction.
        width=yl*c.expm1(MD*cv(r-l));suffix=yr*c.expm1(MD*cv(1-r))
    else:width=cv(r-l);suffix=cv(11-r)
    mass=incoming.weighted.terminal.local.positive_kernel_mass(c,width,rate)
    decay=f.scalar(1) if not rate else f.factor((0,0,0,0,0),-suffix*c.mpf(rate.numerator)/rate.denominator)
    return width,suffix,mass,decay


def formal_affine_transport(f,initial,local,width):
    if set(initial)!=set(RATES) or set(local)!=set(RATES):raise ValueError('All five real incoming/local rows required')
    if ep(width)[0]<0:raise ValueError('Nonnegative actual physical logarithmic width required')
    c=f.c;result={}
    for name,rate in RATES.items():
        if len(initial[name])!=2 or len(local[name])!=2:raise ValueError('Real C0 and ordinary Z pairs required')
        decay=f.scalar(1) if not rate else f.factor((0,0,0,0,0),-width*c.mpf(rate.numerator)/rate.denominator)
        result[name]=[row*decay+local[name][n] for n,row in enumerate(initial[name])]
    return result


class OriginalO2AxialBufferFiniteN:
    mode='current_original_O2_axial_buffer_whole_source_and_actual_finite_N_slope_exit_incoming'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalO2SlopeFiniteN(dps)
        self.c=self.upstream.c;self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes)
        self.bindings=source_bindings();self.parameters=parameter_binding(self.family,self.hashes);self.cache={};self.cells={};self.parents={}
        self.saved=json.loads(gzip.decompress((HERE/previous.NAME).read_bytes()))
        if not self.saved[previous.GATE] or self.saved['source_family']!=self.family or self.saved['candidate_N']!=257:
            raise ValueError('Accepted same-current actual slope-exit prefix required')
        for name,digest in self.saved['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
        for module in (kernels,previous.original,switch):fields.previous.bind(self.hashes,Path(module.__file__).name,sha(Path(module.__file__).name))
        for name in (previous.NAME,Path(__file__).name):fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or not checked.get(GATE) or checked['source_family']!=self.family:
                raise ValueError('Accepted actual O2 axial/buffer finite-N receipt required')
            for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def owner(self,label):return self.upstream.owner(label)

    def query(self,label,chart,left,right,N=257):
        N=reference.phase.candidate_N(N)
        if N!=257 or label not in ('0','.5'):raise ValueError('Same current frames0,.5,N257 required')
        lo,hi=fraction(chart,left),fraction(chart,right);key=(label,chart,lo,hi,N)
        if key in self.cells:return self.cells[key]
        op=self.owner(label);f,c=op.flow,op.c
        with mp.workdps(c.dps+40):
            if label not in self.parents:self.parents[label]=self.upstream.query(label,(1,1),(1,1))
            source=background_cell(op,self.parents[label],chart,left,right);a=source['actual_a_axial5']
            E=source['raw']['raw_current_radius_y_derivative_axial_coefficients']['velocity']['theta'][0]
            proxy=SimpleNamespace(flow=f,c=c,reference=op.reference,zrows=op.zrows,P0=op.P0,
                Pstar=op.Pstar,source_radius=source['actual_source_radius'],correlated_C=f.multiply(a,E))
            generic=reference.long.RECOVER(proxy,source['raw'])
            generic.pop('source_frame_conditional_on_same_actual_Rm_inlet')
            generic['source_from_same_current_original_O2_axial_buffer_background']=True
            parameters=self.upstream.upstream.upstream.upstream.parameters
            roots,qr,quotients=CRITICAL_QUOTIENTS(proxy,generic,a,parameters.eta_log)
            got=reference.primitives.all_u_primitive_bounds(f,roots,qr,parameters.dstar_log,c.mpf([0,1]))
            if qr[C0].zero and qr[Z].zero:
                got['values']={name:f.scalar(0) for name in got['values']}
                got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
            V=generic['common_velocity_V_axial5']
            density=reference.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],N)
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=op.P0,
                original_closed_O2_axial_buffer_background=source,original_generic_source=generic,
                original_full_signed_quotients=quotients,original_roots=roots['roots'],original_q_rows=qr,
                original_primitive_values=got['values'],original_primitive_proof=got['record'],
                original_signed_five_density_C0_Z=density,exact_source_radius=source['actual_source_radius'],
                actual_compact_source_coordinate_cell=(left,right),actual_chart=chart,
                actual_global_radius_phase='frac(N*(logRref+y-logRa-hb*s_c/2))',
                actual_phase_full_period_cover=c.mpf([0,1]),actual_phase_Z_exact_zero=True,
                actual_physical_log_measure='dy=40*exp(40*phase)*dphase' if chart=='axial' else 'dy=dtau',
                exact_a2_collected_positive_Delta_and_axis_Z_inertia_retained=True,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cells[key]=result;return result

    def seam(self,label,left,right):
        op=self.owner(label);rows=0
        l=self.parents[label] if left=='slope_exit' else self.query(label,'axial',(1,1),(1,1))
        r=self.query(label,'axial',(0,1),(0,1)) if right=='axial_inlet' else self.query(label,'buffer',(0,1),(0,1))
        l,r=l['original_generic_source'],r['original_generic_source']
        pairs=[(l['common_velocity_E_axial5'],r['common_velocity_E_axial5']),
            (l['common_velocity_V_axial5'],r['common_velocity_V_axial5'])]
        pairs.extend((l['common_own_five_histories_axial5'][k],r['common_own_five_histories_axial5'][k]) for k in RATES)
        for a,b in pairs:
            for x,y in zip(a,b):
                d=x-y
                if not d.zero and not ep(d.coefficient)[0]<=0<=ep(d.coefficient)[1]:
                    raise ValueError('Original shared axial/buffer source covers must meet')
                rows+=1
        return dict(exact_original_function_join_by_closed_kernel_and_decay_identities=True,
            source_overlap_consistency_rows=rows,overlap_not_function_identity_proof=True,
            same_original_background_P0_and_history_seed=True)

    def contribution(self,label,N=257):
        N=reference.phase.candidate_N(N)
        if N!=257 or label not in ('0','.5'):raise ValueError('Same supplied current frame/frequency required')
        if (label,N) in self.cache:return self.cache[(label,N)]
        op=self.owner(label);f,c=op.flow,op.c;accepted=self.saved['frames'][label]
        with mp.workdps(c.dps+40):
            guard_incoming(self.family,label,N,op,accepted)
            restore=lambda row:switch.first.endpoint.restore_row(f,row)
            initial={name:[restore(row) for row in values] for name,values in incoming.row_packet(accepted,'actual_current_O2_slope_exit_correction_C0_Z').items()}
            self.query(label,'axial',(0,1),(0,1),N)
            joins=dict(slope_axial=self.seam(label,'slope_exit','axial_inlet'))
            current=initial;charts={}
            for chart,partition in PARTITIONS.items():
                local={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
                if chart=='buffer':joins['axial_buffer']=self.seam(label,'axial_exit','buffer_inlet')
                for left,right in zip(partition,partition[1:]):
                    source=self.query(label,chart,left,right,N);density=source['original_signed_five_density_C0_Z'];drivers={};weights={}
                    for name,rate in RATES.items():
                        width,suffix,mass,decay=physical_weights(f,chart,left,right,rate)
                        pair=[reference.bounds.symmetric(f,reference.bounds.magnitude(f,density[part][name])*mass*decay) for part in ('kernels','Z_derivatives')]
                        for n,row in enumerate(pair):local[name][n]+=row
                        drivers[name]=pair;weights[name]=dict(actual_own_rate_mass=mass,actual_suffix_to_chart_exit=decay)
                    cells.append(dict(source=source,actual_cell_to_chart_exit_driver_C0_Z=drivers,
                        actual_original_log_width=width,actual_original_suffix=suffix,original_own_rate_weights=weights))
                    print('Actual current O2 axial/buffer finite-N cell',label,chart,left,right,flush=True)
                width=c.expm1(MD) if chart=='axial' else c.mpf(11)
                output=formal_affine_transport(f,current,local,width)
                charts[chart]=dict(actual_current_chart_incoming_C0_Z=current,actual_chart_local_driver_C0_Z=local,
                    actual_current_chart_exit_C0_Z=output,actual_source_cells=cells,exact_source_partition=partition,
                    full_actual_logarithmic_width=width,pressure_memory_exactly_one=True)
                current=output
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=op.P0,
                exact_source_O2_slope_exit_factor=op.Rm_factor*c.exp(7),
                exact_source_O2_axial_exit_factor=op.Rm_factor*f.factor((0,0,0,0,0),6+c.exp(MD)),
                exact_source_Rd_factor=op.Rm_factor*f.factor((0,0,0,0,0),17+c.exp(MD)),
                actual_current_O2_slope_exit_incoming_correction_C0_Z=initial,actual_O2_axial_buffer_charts=charts,
                actual_current_Rd_correction_C0_Z=current,actual_Rd_background_source=self.query(label,'buffer',(11,1),(11,1)),
                typed_actual_O2_axial_buffer_source_joins=joins,
                genuine_current_finite_N_prefix_through_Rd_boundary_enclosures_supplied=True,
                original_turnoff_kernels_nonzero_far_tail_and_buffer_history_retained=True,
                exact_critical_a2_and_positive_eta_axis_q_retained=True,
                original_P0_background_and_actual_correction_kept_separate=True,
                no_old_O2_owner_or_source_replay=True,
                source_enclosures_not_chosen_field_values_or_functional_closure=True,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cache[(label,N)]=result;return result


def run():
    began=time.monotonic();owner=OriginalO2AxialBufferFiniteN(require_checked=False)
    frames={label:owner.contribution(label) for label in ('0','.5')}
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_source_bindings=owner.bindings,exact_original_parameter_binding=owner.parameters,
        frames=reference.serialized(frames),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Actual current finite-N prefix reaches original Rd after axial turnoff and11-unit buffer',flush=True);return report


if __name__=='__main__':run()

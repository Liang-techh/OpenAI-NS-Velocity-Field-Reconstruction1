"""Actual ordinary Z3 source, same inverse, signed density and five targets.

Preserves every accepted C2 handle. Defining integral functions are separate
from quantitative C3 ranges, third repaired controls and numerical fields.
"""
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_C2_density_transport as parent

phase,current,source=parent.phase,parent.current,parent.source
HERE,PREFIX,sha,ep=parent.HERE,parent.PREFIX,parent.sha,current.ep
NAME=PREFIX+'current_original_C3_source_targets.json.gz'
RECEIPT=PREFIX+'current_original_C3_source_targets_check.json'
GATES=('current_original_actual_17_chart_C3_source_and_same_inverse_functions_installed',
    'current_original_actual_17_chart_C3_density_transport_and_five_target_functions_installed')
OPEN=parent.OPEN


@dataclass(frozen=True)
class C3Function:
    value: object
    Z: object
    ZZ: object
    ZZZ: object


def preserve(old,third):return C3Function(old.value,old.Z,old.ZZ,third)
def rows(q):return (q.value,q.Z,q.ZZ,q.ZZZ)


def encoded(value):
    if isinstance(value,C3Function):return [q.node for q in rows(value)]
    if isinstance(value,dict):return {k:encoded(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encoded(v) for v in value]
    return phase.encoded(value)


def flat_sigma_third(c,x):return 6*phase.flat_source.sigma_jets(c,x)[3]


class C3Algebra:
    def __init__(self,g):self.g=g
    def fixed(self,value):
        q=self.g.constant(value) if isinstance(value,(str,int,Fraction)) else value
        return C3Function(q,self.g.zero,self.g.zero,self.g.zero)
    def add(self,*values):return C3Function(*(self.g.add(*(rows(q)[j] for q in values)) for j in range(4)))
    def neg(self,q):return C3Function(*(self.g.neg(v) for v in rows(q)))
    def mul(self,a,b):
        g=self.g;aa,bb=rows(a),rows(b)
        return C3Function(*(g.add(*(g.mul(g.constant(math.comb(j,k)),aa[k],bb[j-k]) for k in range(j+1))) for j in range(4)))
    def scale(self,a,k):return self.mul(a,self.fixed(k))
    def div(self,n,d):
        g=self.g;nr,dr=rows(n),rows(d);out=[]
        for j in range(4):
            remainder=g.add(*(g.mul(g.constant(math.comb(j,k)),dr[k],out[j-k]) for k in range(1,j+1)))
            out.append(g.quotient(g.sub(nr[j],remainder),d.value,'same actual positive source denominator'))
        return C3Function(*out)
    def sqrt(self,q):
        g=self.g;v=g.unary('positive_sqrt',q.value);out=[v];den=g.mul(g.constant(2),v)
        for j in range(1,4):
            remainder=g.add(*(g.mul(g.constant(math.comb(j,k)),out[k],out[j-k]) for k in range(1,j)))
            out.append(g.quotient(g.sub(rows(q)[j],remainder),den,'actual positive square root'))
        return C3Function(*out)
    def sigma(self,q):
        g=self.g;first=g.unary('original_flat_sigma_prime',q.value)
        second=g.node('analytic_unary',name='original_flat_sigma_second',argument=q.value.node,
            defining_source=current.ast_binding(phase.flat_sigma_second),
            derivative_source=current.ast_binding(phase.flat_source.sigma_jets),
            ordinary_derivative='2*original_sigma_Taylor_coefficient[2]',
            range_enclosure_is_not_selected_function_value=True)
        third=g.node('analytic_unary',name='original_flat_sigma_third',argument=q.value.node,
            defining_source=current.ast_binding(flat_sigma_third),
            derivative_source=current.ast_binding(phase.flat_source.sigma_jets),
            ordinary_derivative='6*original_sigma_Taylor_coefficient[3]',
            range_enclosure_is_not_selected_function_value=True)
        return C3Function(g.unary('original_flat_sigma',q.value),g.mul(first,q.Z),
            g.add(g.mul(second,q.Z,q.Z),g.mul(first,q.ZZ)),
            g.add(g.mul(third,q.Z,q.Z,q.Z),g.mul(g.constant(3),second,q.Z,q.ZZ),g.mul(first,q.ZZZ)))
    def integral(self,q,angle):
        g=self.g
        return C3Function(*(current.controls.integral(g,v,angle,g.zero,g.symbol(angle),
            measure='original angle dpsi at fixed slow source parameters',
            ordinary_slow_Z_derivative_order=j,
            angle_endpoint_is_fixed_for_partial_Z_derivative=True,
            same_original_Z_independent_global_phase=True) for j,v in enumerate(rows(q))))


def primitive_C3(field,chart):
    g=field.phase.built['graph'];alg=C3Algebra(g);old=field.phase.functions[chart];window=field.phase.windows[chart]
    roots={}
    for name,q in old['source_roots'].items():
        definition=dict(g.nodes[q.value.node]);definition['Z_order']=3
        definition.update(ordinary_slow_Z_derivative_order=3,Taylor_coefficient_factorial=6,
            source_projection_binding=current.ast_binding(CurrentC3SourceTargets.source_packet),
            source_projection='6*same actual current source Taylor coefficient[3]',
            derivative_of_range_endpoint=False,C3_quantity_paths=g.nodes[q.ZZ.node]['C2_quantity_paths'])
        roots[name]=preserve(q,g.node(definition.pop('operation'),**definition))
    if window['loop'] is None:
        return dict(chart=chart,source_roots=roots,A=alg.fixed(0),B=alg.fixed(0),
            exact_flat_from_original_power_admission=True)
    parameters={n['name']:source.FunctionRef(g,i) for i,n in enumerate(g.nodes)
        if n['operation']=='current_original_source_parameter'}
    eta,dstar=[alg.fixed(parameters[k]) for k in ('eta','d_star')]
    a,b,E,p2,t0,D=[roots[k] for k in ('a','b','E','p2','t0','Delta')]
    qold=g.nodes[window['loop']['q'][0]]
    def flat(v):
        attrs={key:value for key,value in qold.items() if key not in ('operation','active_body')}
        attrs.update(ordinary_slow_Z_derivative_order=3,not_additional_Taylor_conversion=True)
        return g.node('original_lazy_flat_branch',active_body=v.node,**attrs)
    gamma=alg.add(alg.scale(eta,2),alg.neg(D));cut=alg.add(alg.fixed(1),alg.neg(alg.div(D,eta)))
    rawq=alg.mul(alg.sigma(cut),alg.sqrt(alg.div(gamma,alg.scale(a,2))))
    q=preserve(old['q'],flat(rawq.ZZZ))
    u=alg.div(alg.mul(p2,q),dstar);hinv=alg.div(alg.fixed(1),alg.sqrt(alg.add(alg.fixed(1),alg.mul(u,u))))
    r=alg.mul(u,hinv);alpha=alg.scale(alg.mul(q,hinv),2)
    angle='loop_angle_'+chart;psi=g.symbol(angle);cp=g.unary('cos',psi);sp=g.unary('sin',psi)
    nn=alg.add(alg.fixed(cp),alg.neg(r))
    den=alg.add(alg.fixed(1),alg.neg(alg.scale(r,g.mul(g.constant(2),cp))),alg.mul(r,r))
    w=alg.div(nn,den);t=alg.add(t0,alg.mul(alpha,w))
    pi=g.node('mathematical_pi');twopi=g.mul(g.constant(2),pi)
    K=alg.div(alg.fixed(1),alg.scale(alg.add(alg.fixed(1),alg.mul(t0,t0),alg.scale(alg.mul(q,q),2)),twopi))
    T1=alg.integral(t,angle);P=alg.add(alg.fixed(psi),alg.integral(alg.mul(t,t),angle));Phi=alg.mul(K,P)
    denpsi=alg.scale(r,g.mul(g.constant(2),sp));denpsipsi=alg.scale(r,g.mul(g.constant(2),cp))
    wpsi=alg.div(alg.add(alg.fixed(g.neg(sp)),alg.neg(alg.mul(w,denpsi))),den)
    wpsipsi=alg.div(alg.add(alg.fixed(g.neg(cp)),alg.neg(alg.scale(alg.mul(wpsi,denpsi),2)),alg.neg(alg.mul(w,denpsipsi))),den)
    tpsi=alg.mul(alpha,wpsi);tpsipsi=alg.mul(alpha,wpsipsi)
    phi_psi=alg.mul(K,alg.add(alg.fixed(1),alg.mul(t,t)))
    phi_psipsi=alg.scale(alg.mul(K,alg.mul(t,tpsi)),2)
    phi_psipsipsi=g.mul(g.constant(2),K.value,g.add(g.mul(tpsi.value,tpsi.value),g.mul(t.value,tpsipsi.value)))
    inverse=old['inverse'];p1,p2inverse=old['psi_Z'],old['psi_ZZ']
    def at(v):return g.node('substitute_original_inverse_angle',body=v.node,inverse_angle=inverse.node,angle_variable=angle)
    numerator=g.add(at(Phi.ZZZ),g.mul(g.constant(3),at(phi_psi.ZZ),p1),
        g.mul(g.constant(3),at(phi_psipsi.Z),p1,p1),g.mul(at(phi_psipsipsi),p1,p1,p1),
        g.mul(g.constant(3),g.add(at(phi_psi.Z),g.mul(at(phi_psipsi.value),p1)),p2inverse))
    p3=g.neg(g.quotient(numerator,at(phi_psi.value),'same strictly positive inverse Jacobian'))
    composed=C3Function(at(T1.value),g.add(at(T1.Z),g.mul(at(t.value),p1)),
        g.add(at(T1.ZZ),g.mul(g.constant(2),at(t.Z),p1),g.mul(at(tpsi.value),p1,p1),g.mul(at(t.value),p2inverse)),
        g.add(at(T1.ZZZ),g.mul(g.constant(3),at(t.ZZ),p1),g.mul(g.constant(3),at(tpsi.Z),p1,p1),
            g.mul(at(tpsipsi.value),p1,p1,p1),g.mul(g.constant(3),at(t.Z),p2inverse),
            g.mul(g.constant(3),at(tpsi.value),p1,p2inverse),g.mul(at(t.value),p3)))
    phase_ref=source.FunctionRef(g,window['phase'])
    chi=C3Function(g.sub(phase_ref,g.quotient(inverse,twopi,'2*pi positive')),
        *[g.neg(g.quotient(v,twopi,'2*pi positive')) for v in (p1,p2inverse,p3)])
    Araw=alg.scale(alg.mul(a,chi),'1/2')
    M=alg.add(alg.neg(alg.div(alg.mul(a,composed),alg.fixed(twopi))),alg.neg(alg.scale(b,phase_ref)))
    Braw=alg.scale(alg.mul(E,M),'1/2')
    return dict(chart=chart,source_roots=roots,q=q,inverse=inverse,psi_Z=p1,psi_ZZ=p2inverse,psi_ZZZ=p3,
        Phi_C3=Phi,Phi_psi_C3=phi_psi,Phi_psipsi_C3=phi_psipsi,Phi_psipsipsi=phi_psipsipsi,
        original_t_C3=t,original_tpsi_C3=tpsi,original_tpsipsi_C3=tpsipsi,fixed_angle_T1_C3=T1,composed_T1_C3=composed,
        A=preserve(old['A'],flat(Araw.ZZZ)),B=preserve(old['B'],flat(Braw.ZZZ)),
        current_C2_prefix_and_original_inverse_retained=True,original_flat_collar_predicates_retained=True,
        phase_and_radius_Z_ZZ_ZZZ_exact_zero=True,phase_ZZZ_exact_zero=True)


def coefficient_C3(g,E,V,A,B,N,old):
    alg=C3Algebra(g);invN=g.quotient(g.one,N,'same exact positive integer N')
    arg=g.mul(A.value,invN);exponential=g.unary('exp',arg);exprel=g.unary('exprel',arg)
    F3=g.add(g.mul(E.ZZZ,A.value,exprel),g.mul(g.constant(3),E.ZZ,A.Z,exponential),
        g.mul(g.constant(3),E.Z,exponential,g.add(A.ZZ,g.mul(A.Z,A.Z,invN))),
        g.mul(E.value,exponential,g.add(A.ZZZ,g.mul(g.constant(3),A.Z,A.ZZ,invN),g.mul(A.Z,A.Z,A.Z,invN,invN))))
    F=preserve(old['F_N'],F3)
    EF,VF,EB,VB=[alg.mul(a,b) for a,b in ((E,F),(V,F),(E,B),(V,B))]
    F2,B2=alg.mul(F,F),alg.mul(B,B);zero=alg.fixed(0)
    first=dict(m=B,h=F,k=alg.add(VF,EB),e=alg.add(alg.scale(VB,2),alg.neg(EF)),p=EF)
    second=dict(m=zero,h=zero,k=alg.mul(F,B),e=alg.add(B2,alg.neg(alg.scale(F2,'1/2'))),p=alg.scale(F2,'1/2'))
    coefficients={order:{key:preserve(old['signed_coefficients'][order][key],q.ZZZ) for key,q in values.items()}
        for order,values in ((-1,first),(-2,second))}
    return coefficients,F


def build_transport(field):
    g=field.phase.built['graph'];alg=C3Algebra(g);N=field.phase.built['N'];ref=lambda i:source.FunctionRef(g,i)
    history={key:alg.fixed(0) for key in current.RATES};windows=[]
    for chart,old in zip(phase.CHARTS,field.parent.functions['windows']):
        base=field.phase.windows[chart];primitive=field.primitives[chart]
        if base['quiet_local_source_zero']:
            density={key:alg.fixed(0) for key in current.RATES};coefficients=None;F=None
        else:
            roots=primitive['source_roots'];coefficients,F=coefficient_C3(g,roots['E'],roots['V'],primitive['A'],primitive['B'],N,old)
            epsilon=alg.fixed(g.quotient(g.one,N,'same exact positive integer N'))
            density={key:preserve(old['density'][key],alg.add(coefficients[-1][key],alg.mul(epsilon,coefficients[-2][key])).ZZZ)
                for key in current.RATES}
        before=history;after={};contributions={}
        for key,rate in current.RATES.items():
            kernel=g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(ref(base['right_offset']),ref(base['offset'])))))
            q=density[key].ZZZ
            third=g.zero if q==g.zero else g.add(*(current.controls.integral(g,g.mul(kernel,q,ref(base['Jacobian'])),
                'native_'+chart,ref(p['lower']),ref(p['upper']),measure='original dlogR: Jacobian exactly once',
                exact_native_left=p['exact_left'],exact_native_right=p['exact_right'],actual_current_chart=chart,
                N_scaled_density=True,ordinary_slow_Z_derivative_order=3,
                integration_endpoints_kernel_radius_and_global_phase_Z_independent=True)
                for p in base['actual_current_source_cells']))
            contributions[key]=preserve(old['contributions'][key],third)
            after[key]=preserve(old['outgoing'][key],g.add(g.mul(ref(base['memory'][key]),before[key].ZZZ),third))
        history=after;windows.append(dict(chart=chart,signed_coefficients=coefficients,F_N=F,density=density,
            incoming=before,contributions=contributions,outgoing=history,original_own_rate_memory=base['memory'],
            original_physical_measure_and_partitions=base['actual_current_source_cells'],predecessor_ZZZ_memory_retained=True))
    oldA=field.parent.functions['actual_terminal_amplitude'];definition=dict(g.nodes[oldA.value.node]);definition['Z_order']=3
    definition.update(ordinary_slow_Z_derivative_order=3,Taylor_coefficient_factorial=6,
        source_projection_binding=current.ast_binding(CurrentC3SourceTargets.source_packet),
        source_projection='6*same actual terminal O3_power E Taylor coefficient[3]',derivative_of_range_endpoint=False)
    amplitude=preserve(oldA,g.node(definition.pop('operation'),**definition));AA=alg.mul(amplitude,amplitude)
    mu=alg.fixed(field.phase.built['parameters']['mu']);targets={}
    for row,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
        targets[row]=alg.div(history[key],amplitude if degree==1 else AA)
    joint=alg.add(history['k'],alg.neg(alg.mul(amplitude,history['m'])))
    targets[current.controls.ROWS[1]]=alg.div(joint,alg.mul(mu,AA))
    targets={key:preserve(field.parent.target_functions()[key],targets[key].ZZZ) for key in current.controls.ROWS}
    return dict(windows=windows,terminal_N_scaled_histories=history,actual_terminal_amplitude=amplitude,
        actual_joint_angular_numerator=joint,actual_N_scaled_five_targets=targets,
        original_C2_signed_function_prefix_retained=True,quantitative_C3_ranges_installed=False,
        actual_C3_repaired_limit_controls_installed=False)


class CurrentC3SourceTargets:
    def __init__(self,target_field=None,owner=None,require_checked=True):
        self.parent=target_field if target_field is not None else parent.CurrentC2DensityTransport(owner=owner)
        if not self.parent.acceptance_loaded:raise ValueError('Accepted original C2 source/target functions required')
        self.phase=self.parent.phase;self.identity=self.parent.identity;self.cache={}
        self.prefix=[dict(n) for n in self.phase.built['graph'].nodes];self.hashes=dict(self.parent.hashes)
        for name in (parent.NAME,parent.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.primitives={chart:primitive_C3(self,chart) for chart in phase.CHARTS}
        self.functions=build_transport(self);self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual C3 source/targets required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual C3 source '+name)
            self.acceptance_loaded=True

    def source_packet(self,ends,chart,left,right=None):
        key=(tuple(ends),chart,left,right)
        if key not in self.cache:
            lower=self.phase.source_packet(ends,chart,left,right)
            packet=self.phase.leading_source(ends,chart,left,right);op=self.phase.outer.owner.owner(ends);f=op.flow
            raw=packet['original_generic_source'];proof=packet['original_full_source_quotients']
            def get(*names):
                for name in names:
                    if name in proof:return proof[name]
                raise ValueError('Actual third source Taylor arrays absent '+str(names))
            zero=[f.scalar(0)]*6
            values=dict(E=raw['common_velocity_E_axial5'],V=raw['common_velocity_V_axial5'],
                a=get('actual_a_axial5','actual_correlated_a_axial5','full_original_shear_a_axial5'),
                b=zero if proof.get('exact_b_and_t0_zero') else get('actual_b_axial5','actual_nonzero_b_axial5'),
                t0=zero if proof.get('exact_b_and_t0_zero') else get('actual_t0_axial5'),
                p2=get('full_signed_p2_axial4','exact_full_p2_axial4'),Delta=get('actual_Delta_axial5','Delta_axial5'))
            projected={}
            for name,row in values.items():
                if len(row)<4 or any(q.ctx is not self.phase.outer.c or q.scale.bases is not f.logs or q.ledger is not f.ledger for q in row):
                    raise ValueError('Actual Taylor3 in same live context/basis/ledger required')
                projected[name]=[row[0],row[1],2*row[2],6*row[3]]
                if not current.current.previous.equivalent_rows(projected[name][:3],lower['actual_raw_root_ordinary_Z2'][name]):
                    raise ValueError('Actual C2 source prefix changed '+name)
            if (packet['source_identity']!=self.identity or packet['exact_common_P0_axial5'] is not op.P0
                or packet.get('actual_phase_Z_exact_zero') is not True):
                raise ValueError('Same actual independent P0 required')
            # Generic recovery retains the same P0 rows but may copy the list.
            if (len(op.P0)<4 or
                not current.current.previous.equivalent_rows(raw['common_original_P0_axial5'],op.P0) or
                any(q.ctx is not self.phase.outer.c or q.scale.bases is not f.logs or q.ledger is not f.ledger for q in op.P0)):
                raise ValueError('Same independent P0 Taylor3 rows/context/basis/ledger required')
            self.cache[key]=dict(source_family=self.identity,exact_Z_cell=list(ends),chart=chart,
                actual_raw_root_ordinary_Z3=projected,actual_independent_P0_ordinary_Z3=[op.P0[0],op.P0[1],2*op.P0[2],6*op.P0[3]],
                actual_C3_source_function_handles=encoded(self.primitives[chart]['source_roots']),
                third_Taylor_factorial_applied_once=True,same_live_basis_context_and_ledger=True,
                live_original_phase_Z_exact_zero=packet['actual_phase_Z_exact_zero'],
                original_C2_source_rows_retained=True,source_radius_and_global_phase_Z_independent=True,
                source_ranges_not_function_values=True)
        return self.cache[key]

    def target_functions(self):return self.functions['actual_N_scaled_five_targets']


def run(target_field=None,owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC3SourceTargets(target_field=target_field,owner=owner,require_checked=False)
        coordinates={chart:(0,1) for chart in phase.CHARTS}
        coordinates.update(first_micro='inlet',second_micro=(1,1),actual_patch=(5,4),Rh_reference=(-5,1),O3_power=(2,1))
        packets=[field.source_packet(ends,chart,coordinates[chart]) for ends in current.current.CELLS for chart in phase.CHARTS]
        partitions=current.native_partitions()
        cells=[field.source_packet(ends,chart,left,right) for ends in current.current.CELLS for chart in phase.CHARTS
            for left,right in zip(partitions[chart],partitions[chart][1:])]
        report=dict(candidate_actual_C3_source_targets_constructed=True,source_family=field.identity,
            original_C2_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.phase.built['graph'].nodes,
            actual_C3_primitive_functions=encoded(field.primitives),actual_C3_density_transport_targets=encoded(field.functions),
            actual_four_Z_17_chart_third_source_packets=encoded(packets),
            actual_four_Z_full_native_partition_third_source_packets=encoded(cells),
            source_bindings=dict(algebra=current.ast_binding(C3Algebra),primitive_C3=current.ast_binding(primitive_C3),
                coefficient_C3=current.ast_binding(coefficient_C3),build_transport=current.ast_binding(build_transport),
                source_projection=current.ast_binding(CurrentC3SourceTargets.source_packet),
                flat_sigma_third=current.ast_binding(flat_sigma_third),original_flat_sigma_jets=current.ast_binding(phase.flat_source.sigma_jets)),
            regularity='Same original C2 source/inverse branch. Source Taylor3 is retained; original sigma is C-infinity flat. '
                'a,active gamma,Poisson D and Phi_psi stay positive on the same compact admitted domains. '
                'Ordinary third differentiation therefore preserves the same branch and finite integral functions; no quantitative C3 caps here.',
            quantitative_C3_ranges_installed=False,actual_C3_repaired_limit_controls_installed=False,
            current_numeric_point_field_oracle_installed=False,**dict.fromkeys(GATES+OPEN,False),
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Actual C3 source, same inverse, signed density, own-rate transport and five target functions constructed',flush=True)
    return field


if __name__=='__main__':run()

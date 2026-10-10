"""Independent ordinary cubic calculus and actual source/transport replay.

Checks function definitions only. No magnitude cap is a signed value, and
this receipt does not install C3 limit controls or a physical point oracle.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C3_source_targets as current

phase,parent=current.phase,current.parent


def symbolic(g,node,bindings):
    if node in bindings:return bindings[node]
    n=g.nodes[node];op=n['operation'];at=lambda i:symbolic(g,i,bindings)
    if op=='exact_rational':return s.Rational(n['numerator'],n['denominator'])
    if op=='mathematical_pi':return s.pi
    if op=='sum':return s.Add(*(at(i) for i in n['arguments']))
    if op=='product':return s.Mul(*(at(i) for i in n['arguments']))
    if op=='negative':return -at(n['argument'])
    if op=='positive_quotient':return at(n['numerator'])/at(n['denominator'])
    if op=='original_lazy_flat_branch':return at(n['active_body'])
    if op=='analytic_unary':
        v=at(n['argument']);name=n['name']
        if name=='exp':return s.exp(v)
        if name=='exprel':return (s.exp(v)-1)/v
        if name=='positive_sqrt':return s.sqrt(v)
        sigma_names=('original_flat_sigma','original_flat_sigma_prime',
            'original_flat_sigma_second','original_flat_sigma_third')
        if name in sigma_names:
            x=s.Symbol('sigma_argument');f=s.Function('sigma')(x)
            return s.diff(f,x,sigma_names.index(name)).subs(x,v)
    raise ValueError('Unbound actual C3 symbolic function '+str(n))


def independent_algebra():
    g=current.source.FunctionTransportGraph();alg=current.C3Algebra(g);z=s.Symbol('Z');bindings={}
    def source(name):
        value=s.Function(name)(z);handles=[]
        for j in range(4):
            q=g.symbol(name+'_'+str(j));bindings[q.node]=s.diff(value,z,j);handles.append(q)
        return current.C3Function(*handles),value
    a,av=source('a');b,bv=source('b');count=0
    for q,want in ((alg.mul(a,b),av*bv),(alg.div(a,b),av/bv),
                   (alg.sqrt(a),s.sqrt(av)),(alg.sigma(a),s.Function('sigma')(av))):
        for j,h in enumerate(current.rows(q)):
            assert s.simplify(symbolic(g,h.node,bindings)-s.diff(want,z,j))==0
            count+=1
    return dict(independent_ordinary_product_quotient_sqrt_sigma_rows=count,
        cubic_Leibniz_coefficients_and_sigma_chain_checked=True)


def independent_inverse_and_primitives(field):
    """Compare actual DAG with unrestricted bivariate polynomial local jets."""
    g=field.phase.built['graph'];f=field.primitives['reference'];z,x=s.symbols('Z psi')
    p1,p2,p3=s.symbols('psi1 psi2 psi3');angle=p1*z+p2*z*z/2+p3*z**3/6
    def jet(label):
        coeff={(i,j):s.Symbol(label+str(i)+str(j)) for i in range(4) for j in range(4-i)}
        value=s.Add(*(v*x**i*z**j/(s.factorial(i)*s.factorial(j)) for (i,j),v in coeff.items()))
        return coeff,value
    pc,Phi=jet('phi');hc,H=jet('H');bindings={f['psi_Z'].node:p1,f['psi_ZZ'].node:p2}
    body={}
    for key,i in (('Phi_C3',0),('Phi_psi_C3',1),('Phi_psipsi_C3',2),
                  ('fixed_angle_T1_C3',0),('original_t_C3',1),('original_tpsi_C3',2),('original_tpsipsi_C3',3)):
        table=pc if key.startswith('Phi') else hc
        for j,q in enumerate(current.rows(f[key])):
            body[q.node]=table.get((i,j),s.Symbol(key+'_'+str(j)))
    body[f['Phi_psipsipsi'].node]=pc[3,0]
    for node,n in enumerate(g.nodes):
        if n['operation']=='substitute_original_inverse_angle' and n['inverse_angle']==f['inverse'].node and n['body'] in body:
            bindings[node]=body[n['body']]
    third=s.diff(Phi.subs(x,angle),z,3).subs(z,0)
    want=-s.expand(third-pc[1,0]*p3)/pc[1,0]
    actual=symbolic(g,f['psi_ZZZ'].node,bindings)
    assert s.simplify(actual-want)==0
    bindings[f['psi_ZZZ'].node]=p3
    composed=H.subs(x,angle);T1=f['composed_T1_C3']
    for j,q in enumerate(current.rows(T1)):
        assert s.simplify(symbolic(g,q.node,bindings)-s.diff(composed,z,j).subs(z,0))==0
    # Bind the retained C2 source rows and test the actual cubic A/B bodies.
    values={}
    for name,q in f['source_roots'].items():
        polynomial=sum(s.Symbol(name+str(j))*z**j/s.factorial(j) for j in range(4))
        values[name]=polynomial
        for j,handle in enumerate(current.rows(q)):bindings[handle.node]=s.diff(polynomial,z,j).subs(z,0)
    phase_value,psi0=s.symbols('phase psi0')
    bindings[field.phase.windows['reference']['phase']]=phase_value
    bindings[f['inverse'].node]=psi0
    psi=psi0+angle;T=H.subs(x,angle)
    A=values['a']*(phase_value-psi/(2*s.pi))/2
    M=-values['a']*T/(2*s.pi)-values['b']*phase_value
    B=values['E']*M/2
    for key,value in (('A',A),('B',B)):
        assert s.simplify(symbolic(g,f[key].ZZZ.node,bindings)-s.diff(value,z,3).subs(z,0))==0
    assert f['phase_ZZZ_exact_zero']
    return dict(unrestricted_local_cubic_inverse_identity_checked=True,
        fixed_endpoint_and_composed_T1_rows_checked=4,actual_A_B_third_bodies_checked=2,
        same_C2_inverse_and_Z_independent_global_phase=True,
        scope='exact local jet identities, not a numerical current point-field oracle')


def independent_density():
    g=current.source.FunctionTransportGraph();z=s.Symbol('Z');N=s.Symbol('N',positive=True);bindings={}
    def source(name):
        value=s.Function(name)(z);handles=[]
        for j in range(4):
            q=g.symbol(name+'_'+str(j));bindings[q.node]=s.diff(value,z,j);handles.append(q)
        return current.C3Function(*handles),value
    E,e=source('E');V,v=source('V');A,a=source('A');B,b=source('B')
    n=g.symbol('N');bindings[n.node]=N
    lower=lambda q:phase.C2Function(*current.rows(q)[:3])
    old_pairs,oldF=parent.coefficient_C2(g,lower(E),lower(V),lower(A),lower(B),n)
    pairs,F=current.coefficient_C3(g,E,V,A,B,n,dict(F_N=oldF,signed_coefficients=old_pairs))
    full=N*e*(s.exp(a/N)-1);zero=s.Integer(0)
    expected={-1:dict(m=b,h=full,k=v*full+e*b,e=2*v*b-e*full,p=e*full),
        -2:dict(m=zero,h=zero,k=full*b,e=b*b-full**2/2,p=full**2/2)}
    assert s.simplify(symbolic(g,F.ZZZ.node,bindings)-s.diff(full,z,3))==0
    for order,values in pairs.items():
        for key,q in values.items():
            assert current.rows(q)[:3]==tuple(vars(old_pairs[order][key]).values())
            assert s.simplify(symbolic(g,q.ZZZ.node,bindings)-s.diff(expected[order][key],z,3))==0
    return dict(independent_full_exponential_third_derivative_checked=True,
        independent_signed_third_density_coefficients_checked=10,
        both_inverse_N_orders_and_original_C2_handles_retained=True)


def independent_normalization(field):
    z=s.Symbol('Z');mu=s.Symbol('mu',positive=True);g=field.phase.built['graph'];bindings={}
    H={name:s.Function(name)(z) for name in current.current.RATES};A=s.Function('Am')(z)
    for name,q in field.functions['terminal_N_scaled_histories'].items():
        for j,handle in enumerate(current.rows(q)):bindings[handle.node]=s.diff(H[name],z,j)
    for j,handle in enumerate(current.rows(field.functions['actual_terminal_amplitude'])):
        bindings[handle.node]=s.diff(A,z,j)
    bindings[field.phase.built['parameters']['mu'].node]=mu
    expected={'M':H['m']/A,'I':H['h']/A,'S':H['e']/A**2,'Cp':H['p']/A**2,
        current.current.controls.ROWS[1]:(H['k']-A*H['m'])/(mu*A**2)}
    for name,q in field.target_functions().items():
        assert s.simplify(symbolic(g,q.ZZZ.node,bindings)-s.diff(expected[name],z,3))==0
        assert current.rows(q)[:3]==tuple(vars(field.parent.target_functions()[name]).values())
    return dict(independent_third_five_target_normalizations_checked=5,
        joint_angular_numerator_and_actual_amplitude_derivatives_retained=True)


def source_bindings():
    return dict(algebra=current.current.ast_binding(current.C3Algebra),
        primitive_C3=current.current.ast_binding(current.primitive_C3),
        coefficient_C3=current.current.ast_binding(current.coefficient_C3),
        build_transport=current.current.ast_binding(current.build_transport),
        source_projection=current.current.ast_binding(current.CurrentC3SourceTargets.source_packet),
        flat_sigma_third=current.current.ast_binding(current.flat_sigma_third),
        original_flat_sigma_jets=current.current.ast_binding(phase.flat_source.sigma_jets))


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_C3_source_targets_constructed']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC3SourceTargets(require_checked=False)
        assert not field.acceptance_loaded and field.parent.acceptance_loaded
        assert raw['source_family']==field.identity
        g=field.phase.built['graph'];ref=lambda node:current.source.FunctionRef(g,node)
        assert g.nodes[:len(field.prefix)]==field.prefix and g.nodes==raw['exact_graph_nodes']
        assert current.encoded(field.primitives)==raw['actual_C3_primitive_functions']
        assert json.loads(json.dumps(current.encoded(field.functions)))==raw['actual_C3_density_transport_targets']
        assert source_bindings()==raw['source_bindings']
        counts=dict(actual_Z_cells=4,live_chart_source_packets=0,actual_ordinary_third_root_rows=0,
            actual_independent_P0_ordinary_rows=0,same_original_inverse_handles=0,
            actual_own_rate_steps=0,third_radial_integrals=0,quiet_nonzero_third_memories=0,
            representative_chart_packets=0,full_native_partition_packets=0,third_fixed_angle_integrals=0)
        coordinates={chart:(0,1) for chart in phase.CHARTS}
        coordinates.update(first_micro='inlet',second_micro=(1,1),actual_patch=(5,4),Rh_reference=(-5,1),O3_power=(2,1))
        representative=[(ends,chart,coordinates[chart],None) for ends in current.current.current.CELLS for chart in phase.CHARTS]
        partitions=current.current.native_partitions()
        cells=[(ends,chart,left,right) for ends in current.current.current.CELLS for chart in phase.CHARTS
            for left,right in zip(partitions[chart],partitions[chart][1:])]
        assert len(representative)==len(raw['actual_four_Z_17_chart_third_source_packets'])==68
        assert len(cells)==len(raw['actual_four_Z_full_native_partition_third_source_packets'])==228
        queries=list(zip(representative,raw['actual_four_Z_17_chart_third_source_packets']))
        queries+=list(zip(cells,raw['actual_four_Z_full_native_partition_third_source_packets']))
        for (ends,chart,left,right),stored in queries:
            packet=field.source_packet(ends,chart,left,right);assert current.encoded(packet)==stored
            owner=field.phase.outer.owner.owner(ends);f=owner.flow
            original=field.phase.leading_source(ends,chart,left,right);proof=original['original_full_source_quotients']
            assert original['actual_phase_Z_exact_zero'] and packet['live_original_phase_Z_exact_zero']
            generic=original['original_generic_source'];get=lambda *names:next(proof[k] for k in names if k in proof)
            zero=[f.scalar(0)]*6
            actual=dict(E=generic['common_velocity_E_axial5'],V=generic['common_velocity_V_axial5'],
                a=get('actual_a_axial5','actual_correlated_a_axial5','full_original_shear_a_axial5'),
                b=zero if proof.get('exact_b_and_t0_zero') else get('actual_b_axial5','actual_nonzero_b_axial5'),
                t0=zero if proof.get('exact_b_and_t0_zero') else get('actual_t0_axial5'),
                p2=get('full_signed_p2_axial4','exact_full_p2_axial4'),Delta=get('actual_Delta_axial5','Delta_axial5'))
            for name,rows in actual.items():
                assert len(rows)>=4
                assert all(q.ctx is field.phase.outer.c and q.scale.bases is f.logs and q.ledger is f.ledger for q in rows)
                assert current.current.current.previous.equivalent_rows(packet['actual_raw_root_ordinary_Z3'][name],
                    [rows[0],rows[1],2*rows[2],6*rows[3]])
                q=field.primitives[chart]['source_roots'][name];old=field.phase.functions[chart]['source_roots'][name]
                assert current.rows(q)[:3]==tuple(vars(old).values())
                n=g.nodes[q.ZZZ.node]
                assert n['native_chart']==chart and n['quantity']==name and n['Z_order']==3
                assert n['Taylor_coefficient_factorial']==6 and not n['derivative_of_range_endpoint']
                assert n['source_family']==field.identity
                assert n['source_projection_binding']==current.current.ast_binding(current.CurrentC3SourceTargets.source_packet)
                counts['actual_ordinary_third_root_rows']+=1
            assert original['exact_common_P0_axial5'] is owner.P0
            assert current.current.current.previous.equivalent_rows(packet['actual_independent_P0_ordinary_Z3'],
                [owner.P0[0],owner.P0[1],2*owner.P0[2],6*owner.P0[3]])
            counts['actual_independent_P0_ordinary_rows']+=4;counts['live_chart_source_packets']+=1
            counts['representative_chart_packets' if right is None else 'full_native_partition_packets']+=1
        assert counts['representative_chart_packets']==68 and counts['full_native_partition_packets']==228
        assert counts['live_chart_source_packets']==296 and counts['actual_ordinary_third_root_rows']==2072
        for chart in phase.CHARTS:
            f=field.primitives[chart];old=field.phase.functions[chart];window=field.phase.windows[chart]
            for key in ('A','B'):assert current.rows(f[key])[:3]==tuple(vars(old[key]).values())
            if window['loop'] is None:
                assert all(q==g.zero for key in ('A','B') for q in current.rows(f[key]))
            else:
                assert f['inverse']==old['inverse'] and f['psi_Z']==old['psi_Z'] and f['psi_ZZ']==old['psi_ZZ']
                assert f['phase_ZZZ_exact_zero'] and f['phase_and_radius_Z_ZZ_ZZZ_exact_zero']
                previous=g.nodes[window['loop']['q'][0]]
                for key in ('q','A','B'):
                    assert current.rows(f[key])[:3]==tuple(vars(old[key]).values())
                    n=g.nodes[f[key].ZZZ.node]
                    assert n['ordinary_slow_Z_derivative_order']==3 and n['not_additional_Taylor_conversion']
                    for attr in ('flat_predicate','flat_value','Delta','eta'):assert n[attr]==previous[attr]
                counts['same_original_inverse_handles']+=1
        for n in g.nodes[len(field.prefix):]:
            if n['operation']=='definite_integral' and n['variable'].startswith('loop_angle_'):
                assert n['angle_endpoint_is_fixed_for_partial_Z_derivative'] and n['same_original_Z_independent_global_phase']
                assert n['ordinary_slow_Z_derivative_order'] in range(4)
                if n['ordinary_slow_Z_derivative_order']==3:counts['third_fixed_angle_integrals']+=1
        assert counts['third_fixed_angle_integrals']==32
        for chart,window,old in zip(phase.CHARTS,field.functions['windows'],field.parent.functions['windows']):
            base=field.phase.windows[chart];assert window['chart']==chart
            assert window['original_physical_measure_and_partitions']==base['actual_current_source_cells']
            assert window['original_own_rate_memory']==base['memory']
            for key,rate in current.current.RATES.items():
                q=window['density'][key];before=window['incoming'][key];after=window['outgoing'][key]
                contribution=window['contributions'][key]
                for name in ('density','contributions','outgoing'):
                    assert current.rows(window[name][key])[:3]==tuple(vars(old[name][key]).values())
                assert after.ZZZ==g.add(g.mul(ref(base['memory'][key]),before.ZZZ),contribution.ZZZ)
                if base['quiet_local_source_zero']:
                    assert q.ZZZ==g.zero and contribution.ZZZ==g.zero and before.ZZZ!=g.zero and after.ZZZ!=g.zero
                    counts['quiet_nonzero_third_memories']+=1
                else:
                    n=g.nodes[contribution.ZZZ.node];ids=n['arguments'] if n['operation']=='sum' else [contribution.ZZZ.node]
                    assert len(ids)==len(base['actual_current_source_cells'])
                    kernel=g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(ref(base['right_offset']),ref(base['offset'])))))
                    for node,p in zip(ids,base['actual_current_source_cells']):
                        integral=g.nodes[node]
                        assert integral['operation']=='definite_integral' and integral['variable']=='native_'+chart
                        assert integral['integrand']==g.mul(kernel,q.ZZZ,ref(base['Jacobian'])).node
                        assert integral['lower']==p['lower'] and integral['upper']==p['upper']
                        assert integral['ordinary_slow_Z_derivative_order']==3
                        assert integral['integration_endpoints_kernel_radius_and_global_phase_Z_independent']
                        assert integral['measure']=='original dlogR: Jacobian exactly once'
                        counts['third_radial_integrals']+=1
                counts['actual_own_rate_steps']+=1
        assert counts['actual_own_rate_steps']==85 and counts['third_radial_integrals']==275
        assert counts['quiet_nonzero_third_memories']==5 and counts['same_original_inverse_handles']==16
        amplitude=field.functions['actual_terminal_amplitude'];n=g.nodes[amplitude.ZZZ.node]
        assert n['native_chart']=='O3_power' and n['quantity']=='E' and n['Z_order']==3
        assert g.nodes[n['coordinate']]==dict(operation='exact_rational',numerator=2,denominator=1)
        assert n['Taylor_coefficient_factorial']==6 and not n['derivative_of_range_endpoint']
        assert n['source_projection_binding']==current.current.ast_binding(current.CurrentC3SourceTargets.source_packet)
        algebra=independent_algebra();inverse=independent_inverse_and_primitives(field)
        density=independent_density();normalization=independent_normalization(field)
        z,y,lam=s.symbols('Z y lambda');u=s.Symbol('s');H0=s.Function('H0')(z);D=s.Function('D')
        H=s.exp(-lam*y)*H0+s.Integral(s.exp(-lam*(y-u))*D(u,z),(u,0,y))
        assert s.simplify(s.diff(H,z,3).subs(y,0)-s.diff(H0,z,3))==0
        assert s.simplify(s.diff(H,z,3,y)+lam*s.diff(H,z,3)-s.diff(D(y,z),z,3))==0
    c=MPIntervalContext();c.dps=50
    for x in (0,1,-1,2):assert current.ep(current.flat_sigma_third(c,c.mpf(x)))==(0,0)
    for x in ('.2','.5','.8'):
        with mp.workdps(70):
            sigma=lambda z:1/(1+mp.exp(-(1/(1-z)**2-1/z**2)))
            want=mp.diff(sigma,mp.mpf(x),3)
            lo,hi=current.ep(current.flat_sigma_third(c,c.mpf(x)));assert lo<=want<=hi
    result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
        actual_replay_counts=counts,independent_cubic_algebra=algebra,
        independent_same_inverse_and_primitive_calculus=inverse,
        independent_signed_density_calculus=density,independent_target_normalization=normalization,
        third_history_FTC_and_inlet_identities_checked=True,original_C2_function_prefix_retained=True,
        original_flat_sigma_third_source_provider_checked=True,
        quantitative_C3_ranges_installed=False,actual_C3_repaired_limit_controls_installed=False,
        current_numeric_point_field_oracle_installed=False,**dict.fromkeys(current.OPEN,False),
        input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual C3 source/inverse/signed-density/five-target function checks passed',flush=True)
    return result


if __name__=='__main__':run()

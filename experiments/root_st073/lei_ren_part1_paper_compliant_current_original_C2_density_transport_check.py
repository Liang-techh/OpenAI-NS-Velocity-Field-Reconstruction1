"""Independent C2 signed-density/normalization and actual transport checks."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C2_density_transport as current

phase=current.phase


def symbolic(g,node,overrides):
    if node in overrides:return overrides[node]
    n=g.nodes[node];op=n['operation'];at=lambda i:symbolic(g,i,overrides)
    if op=='exact_rational':return s.Rational(n['numerator'],n['denominator'])
    if op=='sum':return s.Add(*(at(i) for i in n['arguments']))
    if op=='product':return s.Mul(*(at(i) for i in n['arguments']))
    if op=='negative':return -at(n['argument'])
    if op=='positive_quotient':return at(n['numerator'])/at(n['denominator'])
    if op=='analytic_unary':
        v=at(n['argument'])
        if n['name']=='exp':return s.exp(v)
        if n['name']=='exprel':return (s.exp(v)-1)/v
    raise ValueError('Unbound C2 symbolic function '+str(n))


def independent_density_proof():
    g=current.source.FunctionTransportGraph();z=s.Symbol('Z');N=s.Symbol('N',positive=True);bindings={}
    def source(name):
        f=s.Function(name)(z);handles=[]
        for order in range(3):
            handle=g.symbol(name+'_'+str(order));bindings[handle.node]=s.diff(f,z,order);handles.append(handle)
        return phase.C2Function(*handles),f
    E,e=source('E');V,v=source('V');A,a=source('A');B,b=source('B')
    n=g.symbol('N');bindings[n.node]=N
    coefficients,F=current.coefficient_C2(g,E,V,A,B,n)
    independent=N*e*(s.exp(a/N)-1)
    zero=s.Integer(0)
    expected={-1:dict(m=b,h=independent,k=v*independent+e*b,e=2*v*b-e*independent,p=e*independent),
        -2:dict(m=zero,h=zero,k=independent*b,e=b*b-independent**2/2,p=independent**2/2)}
    count=0
    for order,rows in coefficients.items():
        for name,f in rows.items():
            for derivative,handle in enumerate((f.value,f.Z,f.ZZ)):
                assert s.simplify(symbolic(g,handle.node,bindings)-s.diff(expected[order][name],z,derivative))==0
                count+=1
    assert s.simplify(symbolic(g,F.ZZ.node,bindings)-s.diff(independent,z,2))==0
    return dict(independent_signed_density_rows=count,full_exponential_second_derivative_checked=True,
        both_mixed_product_terms_checked=True,exprel_removable_zero_identity_preserved=True)


def independent_normalization(field):
    z=s.Symbol('Z');mu=s.Symbol('mu',positive=True);g=field.phase.built['graph'];bindings={}
    H={name:s.Function(name)(z) for name in current.current.RATES};A=s.Function('Am')(z)
    for name,q in field.functions['terminal_N_scaled_histories'].items():
        for derivative,handle in enumerate((q.value,q.Z,q.ZZ)):bindings[handle.node]=s.diff(H[name],z,derivative)
    q=field.functions['actual_terminal_amplitude']
    for derivative,handle in enumerate((q.value,q.Z,q.ZZ)):bindings[handle.node]=s.diff(A,z,derivative)
    bindings[field.phase.built['parameters']['mu'].node]=mu
    expected={'M':H['m']/A,'I':H['h']/A,'S':H['e']/A**2,'Cp':H['p']/A**2,
        current.current.controls.ROWS[1]:(H['k']-A*H['m'])/(mu*A**2)}
    for name,target in field.target_functions().items():
        assert s.simplify(symbolic(g,target.ZZ.node,bindings)-s.diff(expected[name],z,2))==0
        old=field.phase.built['N_scaled_targets'][name]
        assert target.value==old.value and target.Z==old.Z
    return dict(independent_second_target_normalizations=5,
        joint_angular_numerator_and_positive_amplitude_derivatives_checked=True)


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_C2_transport_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC2DensityTransport(require_checked=False)
        assert not field.acceptance_loaded and field.phase.acceptance_loaded
        assert raw['source_family']==field.identity
        g=field.phase.built['graph'];ref=lambda node:current.source.FunctionRef(g,node)
        assert g.nodes[:len(field.prefix)]==field.prefix
        assert g.nodes==raw['exact_graph_nodes']
        assert json.loads(json.dumps(phase.encoded(field.functions)))==raw['actual_C2_density_transport_and_target_functions']
        assert raw['actual_source_bindings']==dict(coefficient_C2=current.current.ast_binding(current.coefficient_C2),
            original_C1_coefficient_functions=current.current.ast_binding(current.current.algebra.coefficient_pairs),
            original_C1_target_normalization=current.current.ast_binding(current.current.algebra.normalize_history),
            actual_transport_and_normalization=current.current.ast_binding(current.build))
        density=independent_density_proof();normalization=independent_normalization(field)
        counts=dict(actual_chart_histories=0,actual_original_own_rate_steps=0,
            second_row_radial_integrals=0,quiet_nonzero_predecessor_memories=0)
        for chart,window in zip(phase.CHARTS,field.functions['windows']):
            old=field.phase.windows[chart];assert window['chart']==chart
            assert window['original_physical_measure_and_partitions']==old['actual_current_source_cells']
            assert window['original_own_rate_memory']==old['memory']
            for name,rate in current.current.RATES.items():
                q=window['density'][name];before=window['incoming'][name];after=window['outgoing'][name]
                contribution=window['contributions'][name]
                assert after.value.node==old['outgoing'][name][0] and after.Z.node==old['outgoing'][name][1]
                assert after.ZZ==g.add(g.mul(ref(old['memory'][name]),before.ZZ),contribution.ZZ)
                if old['quiet_local_source_zero']:
                    assert q.ZZ==g.zero and contribution.ZZ==g.zero
                    assert before.ZZ!=g.zero and after.ZZ!=g.zero
                    counts['quiet_nonzero_predecessor_memories']+=1
                else:
                    assert q.value.node==old['density'][name][0] and q.Z.node==old['density'][name][1]
                    n=g.nodes[contribution.ZZ.node]
                    ids=n['arguments'] if n['operation']=='sum' else [contribution.ZZ.node]
                    assert len(ids)==len(old['actual_current_source_cells'])
                    kernel=g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(ref(old['right_offset']),ref(old['offset'])))))
                    for node,p in zip(ids,old['actual_current_source_cells']):
                        integral=g.nodes[node]
                        assert integral['operation']=='definite_integral' and integral['variable']=='native_'+chart
                        assert integral['integrand']==g.mul(kernel,q.ZZ,ref(old['Jacobian'])).node
                        assert integral['lower']==p['lower'] and integral['upper']==p['upper']
                        assert integral['ordinary_slow_Z_derivative_order']==2
                        assert integral['integration_endpoints_kernel_radius_and_global_phase_Z_independent']
                        assert integral['measure']=='original dlogR: Jacobian exactly once'
                        counts['second_row_radial_integrals']+=1
                counts['actual_original_own_rate_steps']+=1
            counts['actual_chart_histories']+=1
        z,y,lam=s.symbols('Z y lambda');q=s.Symbol('s');H0=s.Function('H0')(z);D=s.Function('D')
        H=s.exp(-lam*y)*H0+s.Integral(s.exp(-lam*(y-q))*D(q,z),(q,0,y))
        assert s.simplify(s.diff(H,z,2).subs(y,0)-s.diff(H0,z,2))==0
        assert s.simplify(s.diff(H,z,2,y)+lam*s.diff(H,z,2)-s.diff(D(y,z),z,2))==0
        terminal=field.phase.source_packet(('0','.5'),'O3_power',(2,1))
        assert phase.encoded(terminal)==raw['actual_terminal_source_packet']
        amplitude=field.functions['actual_terminal_amplitude'];n=g.nodes[amplitude.ZZ.node]
        assert n['native_chart']=='O3_power' and n['quantity']=='E' and n['Z_order']==2
        assert g.nodes[n['coordinate']] == dict(operation='exact_rational',numerator=2,denominator=1)
        assert n['source_family']==field.identity and not n['derivative_of_range_endpoint']
        assert n['source_projection_binding']==current.current.ast_binding(phase.CurrentC2PhasePrimitives.source_packet)
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            actual_five_target_ZZ_function_representations_installed=True,
            independent_density_proof=density,independent_target_normalization=normalization,
            actual_transport_replay_counts=counts,original_C1_function_prefix_retained=True,
            same_actual_terminal_amplitude_P0_source_and_global_phase_retained=True,
            current_C2_target_functions_are_enclosure_free_definitions=True,
            quantitative_current_C2_target_ranges_installed=False,actual_C2_repaired_limit_controls_installed=False,
            current_numeric_point_field_oracle_installed=False,**dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),
                Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(phase.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Current original C2 signed density, transport and five target function checks passed',flush=True)
    return result


if __name__=='__main__':run()

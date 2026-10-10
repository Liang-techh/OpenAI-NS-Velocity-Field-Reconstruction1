"""Source-bound C2 projection checks and independent local calculus checks.

The finite scalar fixture tests calculus, not the astronomical current
field's numerical oracle. Actual source replay separately binds all seventeen
charts to the same checked owner, pressure, Taylor basis and ledger.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_C2_phase_primitives as current


def local_calculus_check(field):
    """Evaluate the actual appended DAG against independently solved loops."""
    g=field.built['graph'];chart='reference';functions=field.functions[chart]
    old=field.windows[chart]['loop'];phase_id=field.windows[chart]['phase']
    def sigma(x):
        if x<=0:return mp.mpf(0)
        if x>=1:return mp.mpf(1)
        odds=1/(1-x)**2-1/x**2
        return 1/(1+mp.exp(-odds))
    def roots(z):
        a=mp.mpf('2.1')+mp.mpf('.05')*z+mp.mpf('.03')*z*z
        b=mp.mpf('.1')+mp.mpf('.02')*z+mp.mpf('.01')*z*z
        return dict(a=a,b=b,E=mp.mpf('1.2')+mp.mpf('.1')*z+mp.mpf('.07')*z*z,
            V=3*z,t0=-b/a,p2=mp.mpf('.12')+mp.mpf('.03')*z,Delta=a+b*b/a-2)
    eta,dstar,phase=mp.mpf('.2'),mp.mpf('.7'),mp.mpf('.27')
    def independent(z):
        v=roots(z);q=sigma(1-v['Delta']/eta)*mp.sqrt((2*eta-v['Delta'])/(2*v['a']))
        u=v['p2']*q/dstar;hinv=1/mp.sqrt(1+u*u);r=u*hinv;alpha=2*q*hinv
        def t(x):return v['t0']+alpha*(mp.cos(x)-r)/(1-2*r*mp.cos(x)+r*r)
        K=1/(2*mp.pi*(1+v['t0']**2+2*q*q))
        def Phi(x):return K*(x+mp.quad(lambda angle:t(angle)**2,[0,x]))
        inverse=mp.findroot(lambda x:Phi(x)-phase,(mp.mpf('.5'),mp.mpf('3')))
        T1=mp.quad(t,[0,inverse]);M=-v['a']*T1/(2*mp.pi)-v['b']*phase
        return dict(q=q,psi=inverse,A=v['a']*(phase-inverse/(2*mp.pi))/2,B=v['E']*M/2)
    z=mp.mpf('.15');truth=independent(z);cache={}
    source_data={(name,order):mp.diff(lambda zz:roots(zz)[name],z,order)
        for name in roots(z) for order in range(3)}
    def ev(node,angles=None):
        angles={} if angles is None else angles
        key=(node,tuple(angles.items()))
        if key in cache:return cache[key]
        n=g.nodes[node];op=n['operation'];at=lambda i:ev(i,angles)
        if node==phase_id:v=phase
        elif op=='exact_rational':v=mp.mpf(n['numerator'])/n['denominator']
        elif op=='mathematical_pi':v=mp.pi
        elif op=='current_original_source_parameter':
            v={'eta':eta,'d_star':dstar}[n['name']]
        elif op=='current_original_leading_function_recipe':
            v=source_data[(n['quantity'],n['Z_order'])]
        elif op=='sum':v=sum(at(i) for i in n['arguments'])
        elif op=='product':v=mp.fprod(at(i) for i in n['arguments'])
        elif op=='negative':v=-at(n['argument'])
        elif op=='positive_quotient':v=at(n['numerator'])/at(n['denominator'])
        elif op=='bound_variable':v=angles[n['name']]
        elif op=='analytic_unary':
            x=at(n['argument']);name=n['name']
            if name=='original_flat_sigma':v=sigma(x)
            elif name=='original_flat_sigma_prime':v=mp.diff(sigma,x)
            elif name=='original_flat_sigma_second':v=mp.diff(sigma,x,2)
            else:v={'cos':mp.cos,'sin':mp.sin,'positive_sqrt':mp.sqrt}[name](x)
        elif op=='original_lazy_flat_branch':
            v=at(n['flat_value']) if at(n['Delta'])>=at(n['eta']) else at(n['active_body'])
        elif op=='original_monotone_phase_inverse':v=truth['psi']
        elif op=='substitute_original_inverse_angle':
            v=ev(n['body'],{**angles,n['angle_variable']:at(n['inverse_angle'])})
        elif op=='definite_integral':
            v=mp.quad(lambda x:ev(n['integrand'],{**angles,n['variable']:x}),[at(n['lower']),at(n['upper'])])
        else:raise ValueError('Unsupported local fixture operation '+op)
        cache[key]=v;return v
    errors={}
    for name in ('q','A','B'):
        f=functions[name]
        for order,handle in enumerate((f.value,f.Z,f.ZZ)):
            want=mp.diff(lambda zz:independent(zz)[name],z,order)
            errors[name+'_'+str(order)]=abs(ev(handle.node)-want)
    for order,name in ((1,'psi_Z'),(2,'psi_ZZ')):
        want=mp.diff(lambda zz:independent(zz)['psi'],z,order)
        errors[name]=abs(ev(functions[name].node)-want)
    assert max(errors.values())<mp.mpf('1e-24'),errors
    assert functions['A'].value.node==old['A'][0] and functions['A'].Z.node==old['A'][1]
    assert functions['B'].value.node==old['B'][0] and functions['B'].Z.node==old['B'][1]
    return dict(scope='finite local calculus fixture only; not current point-field oracle',
        genuine_current_numeric_oracle_installed=False,checked_derivative_rows=len(errors),
        maximum_absolute_error=str(max(errors.values())),source_fixture_pressure_not_current_P0=True)


def run(field=None):
    began=time.monotonic()
    raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_C2_source_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC2PhasePrimitives(require_checked=False)
        assert not field.acceptance_loaded and field.identity==raw['source_family']
        assert field.prefix==field.built['graph'].nodes[:len(field.prefix)]
        assert current.encoded(current.source_bindings())==raw['actual_source_bindings']
        assert field.built['graph'].nodes==raw['exact_graph_nodes']
        assert current.encoded(field.functions)==raw['current_C2_phase_functions']
        assert current.regularity_proof()==raw['actual_current_regularity_proof']
        counts=dict(live_chart_source_packets=0,actual_ordinary_Z2_root_rows=0,
            independent_pressure_rows=0,same_original_inverse_handles=0,domain_rejections=0)
        coordinates={chart:(0,1) for chart in current.CHARTS}
        coordinates.update(first_micro='inlet',second_micro=(1,1),actual_patch=(5,4),Rh_reference=(-5,1))
        for chart,stored in zip(current.CHARTS,raw['actual_current_Z2_source_packets']):
            packet=field.source_packet(('0','.5'),chart,coordinates[chart])
            assert current.encoded(packet)==stored
            op=field.outer.owner.owner(('0','.5'))
            assert packet['exact_common_P0_axial5'] is op.P0
            # A fresh direct source query ensures the projection is bound to
            # the retained live Taylor arrays, not serialized endpoints.
            source=field.leading_source(('0','.5'),chart,coordinates[chart])
            f=op.flow;proof=source['original_full_source_quotients'];generic=source['original_generic_source']
            def get(*names):return next(proof[k] for k in names if k in proof)
            zero=[f.scalar(0)]*6
            actual=dict(E=generic['common_velocity_E_axial5'],V=generic['common_velocity_V_axial5'],
                a=get('actual_a_axial5','actual_correlated_a_axial5','full_original_shear_a_axial5'),
                b=zero if proof.get('exact_b_and_t0_zero') else get('actual_b_axial5','actual_nonzero_b_axial5'),
                t0=zero if proof.get('exact_b_and_t0_zero') else get('actual_t0_axial5'),
                p2=get('full_signed_p2_axial4','exact_full_p2_axial4'),Delta=get('actual_Delta_axial5','Delta_axial5'))
            for name,rows in actual.items():
                projected=packet['actual_raw_root_ordinary_Z2'][name]
                assert all(row.ctx is field.outer.c and row.scale.bases is f.logs and row.ledger is f.ledger for row in rows)
                assert current.current.current.previous.equivalent_rows(projected,[rows[0],rows[1],rows[2]*2])
                handle=field.functions[chart]['source_roots'][name].ZZ
                n=field.built['graph'].nodes[handle.node]
                assert n['native_chart']==chart and n['quantity']==name and n['Z_order']==2
                assert n['Taylor_coefficient_factorial']==2 and not n['derivative_of_range_endpoint']
                assert n['source_family']==field.identity
                assert n['source_projection_binding']==current.current.ast_binding(current.CurrentC2PhasePrimitives.source_packet)
                counts['actual_ordinary_Z2_root_rows']+=1
            assert current.current.current.previous.equivalent_rows(packet['actual_independent_P0_ordinary_Z2'],
                [op.P0[0],op.P0[1],op.P0[2]*2])
            counts['independent_pressure_rows']+=3
            functions=field.functions[chart];window=field.windows[chart]
            if window['loop'] is None:
                assert window['quiet_source_proof'] and window['quiet_local_source_zero']
                assert functions['exact_flat_from_original_power_admission']
                assert all(v==field.built['graph'].zero for name in ('A','B') for v in vars(functions[name]).values())
            else:
                assert functions['inverse'].node==window['loop']['inverse']
                assert functions['A'].value.node==window['loop']['A'][0]
                assert functions['A'].Z.node==window['loop']['A'][1]
                assert functions['B'].value.node==window['loop']['B'][0]
                assert functions['B'].Z.node==window['loop']['B'][1]
                old=field.built['graph'].nodes[window['loop']['q'][0]]
                for name in ('A','B'):
                    new=field.built['graph'].nodes[functions[name].ZZ.node]
                    assert new['flat_predicate']==old['flat_predicate'] and new['flat_value']==old['flat_value']
                counts['same_original_inverse_handles']+=1
            counts['live_chart_source_packets']+=1
        for ends,chart,coordinate in ((('0','.5'),'invalid',(0,1)),(('0','0'),'reference',(0,1)),
                                      (('0','.5'),'reference',(-1,1)),(('0','.5'),'first_switch',(2,1))):
            try:field.source_packet(ends,chart,coordinate)
            except ValueError:counts['domain_rejections']+=1
            else:raise AssertionError('Unadmitted C2 source domain must reject')
    with mp.workdps(40):calculus=local_calculus_check(field)
    c=MPIntervalContext();c.dps=50
    for x in (0,1,-1,2):assert current.current.ep(current.flat_sigma_second(c,c.mpf(x)))==(0,0)
    for x in ('.2','.5','.8'):
        with mp.workdps(70):
            xx=mp.mpf(x);sigma=lambda z:1/(1+mp.exp(-(1/(1-z)**2-1/z**2)))
            exact=mp.mpf(0) if x=='.5' else mp.diff(sigma,xx,2)
            lo,hi=current.current.ep(current.flat_sigma_second(c,c.mpf(x)))
            assert lo<=exact<=hi
    result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
        source_charts=list(current.CHARTS),replay_counts=counts,local_calculus_check=calculus,
        original_C1_graph_prefix_and_inverse_handles_preserved=True,
        ordinary_source_second_derivatives_use_actual_Taylor_rows=True,
        same_independent_P0_and_source_basis_ledger_checked=True,
        original_flat_sigma_second_derivative_source_bound_and_checked=True,
        upstream_six_chart_C2_provider_installed=True,actual_five_target_ZZ_installed=False,
        actual_C2_repaired_limit_controls_installed=False,**dict.fromkeys(current.OPEN,False),
        input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Current original seventeen-chart C2 source/phase primitive checks passed',flush=True)
    return result


if __name__=='__main__':run()

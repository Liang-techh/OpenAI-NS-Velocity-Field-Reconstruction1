"""Own histories, all signed moments, physical units and original flat jets."""
import json
from pathlib import Path
import time
from types import SimpleNamespace
import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check as algebra_check


def own_history_fixture():
    g=current.source.FunctionTransportGraph();z=algebra_check.symbol('Z');r=algebra_check.symbol('cut_test_rho');s=algebra_check.symbol('test_s')
    variables={str(q):g.symbol(str(q)) for q in (z,r,s)}
    bind=lambda q:current.source.expression(g,sy.sympify(q),variables)
    jet=lambda q:current.c3.C3Function(*(bind(sy.diff(q,z,k)) for k in range(4)))
    hb=sy.Rational(1,5);eps=sy.Rational(1,9);de=sy.Rational(1,5);rs=4*sy.exp(hb*s)
    f=1+z+r/7;u=2+z+r/3;Np=2+z*z+(1-z)*r
    F=(1+z*z)*r+z*r*r;U=(2-z+z**3)*r+(1+z*z)*r*r
    m=(2-z+z**3)*r/2+(1+z*z)*r*r/3
    t=sy.Symbol('independent_cutoff_t',nonnegative=True)
    primitive=lambda q,lo,hi:sy.integrate(q.subs(r,t),(t,lo,hi))
    P=eps*primitive(2*f*F+Np,0,r)
    chi=1-(r-4)**3/37;Fc=chi*F;Uc=chi*U
    M=primitive(U,0,4)+primitive(Uc,4,r);mc=M/r
    Pc=P.subs(r,4)+eps*primitive(2*f*Fc+Np,4,r)
    radial_rows={};leading={};systems={};functions={};cuts={}
    for chart,radius in (('core',r),('first',rs)):
        sub=lambda q:q.subs(r,radius)
        radial_rows[chart]={name:[jet(sub(sy.diff(q,r,i))) for i in range(3)]
            for name,q in (('F1',F),('Uz1',U),('mean1',m),('K1',m-U))}
        leading[chart]=dict(F_rho_rows=[jet(sub(f))],U_rho_rows=[jet(sub(u))])
        systems[chart]=dict(rho=bind(radius),N1p=jet(sub(Np)))
        functions[chart]=dict(P1=jet(sub(P)))
        cuts[chart]=([jet(1),jet(0),jet(0)] if chart=='core' else [jet(sub(sy.diff(chi,r,i))) for i in range(3)])
    source=SimpleNamespace(parameters=dict(epsilon=bind(eps),hb=bind(hb),delta=bind(de)),
        coordinates=dict(core=variables[str(r)],first=variables[str(s)]),systems=systems,functions=functions,leading=leading)
    fields,integrands,cumulative,bodies=current.own_cutoff_fields(g,source,radial_rows,cuts)
    densities=lambda F1,U1:dict(m11=eps*U1,m12=2*eps**2*r*F1,m13=eps*(2*f*F1+Np),
        m14=2*eps**2*r*(u*F1+U1*f),m15=2*eps*u*U1-eps**2*r*(2*f*F1+Np))
    core_densities=densities(F,U);first_densities=densities(Fc,Uc)
    cache={};count=0
    def check(jet_value,want,label,order=3):
        nonlocal count
        refs=(jet_value.value,jet_value.Z,jet_value.ZZ,jet_value.ZZZ) if order==3 else (jet_value.value,jet_value.Z,jet_value.ZZ)
        for k,ref in enumerate(refs):
            got=algebra_check.symbolic(g,ref,cache)
            assert sy.simplify(got-sy.diff(want,z,k))==0,(label,k)
            count+=1
    for chart,expected in (('core',dict(F1=F,Uz1=U,mean1=m,K1=m-U)),
        ('first',dict(F1=Fc,Uz1=Uc,mean1=mc,K1=mc-Uc))):
        radius=r if chart=='core' else rs
        for name,q in expected.items():
            for i in range(3):check(fields[chart]['rows'][name][i],sy.diff(q,r,i).subs(r,radius),(chart,name,i))
        check(fields[chart]['radial_primitive_M1'],(eps*r*expected['mean1']).subs(r,radius),(chart,'M1'))
        check(fields[chart]['pressure_P1'],(P if chart=='core' else Pc).subs(r,radius),(chart,'P1'))
        Q=(2*z*expected['Uz1']-(1+de)*z*expected['mean1']-(1-z*z)*sy.diff(expected['mean1'],z))/(1-de*z*z)
        for name,q in (('Q1',Q),('V1',eps*r*Q)):
            for i in range(3):check(fields[chart][name][i],sy.diff(q,r,i).subs(r,radius),(chart,name,i),2)
        for name in current.MOMENTS:
            want=(primitive(core_densities[name],0,r) if chart=='core' else
                primitive(core_densities[name],0,4)+primitive(first_densities[name],4,r))
            check(cumulative[chart][name],want.subs(r,radius),(chart,name))
    # The unchanged leading pressure source must survive a zero velocity
    # cutoff; three velocity moments then have zero incremental density.
    zero=current.c3.C3Algebra(g).fixed(0)
    zeros=dict(F1=[zero]*3,Uz1=[zero]*3)
    tail=current.signed_rho_integrands(g,leading['core'],zeros,systems['core'],bind(eps))
    check(tail['m13'],eps*Np,'zero-cut m13')
    check(tail['m15'],-eps**2*r*Np,'zero-cut m15')
    for name in ('m11','m12','m14'):check(tail[name],sy.Integer(0),('zero-cut',name))
    return dict(passed=True,mathematical_fixture_only=True,exact_signed_C3_and_radial_identities=count,
        independent_nonconstant_cutoff_own_histories=True,nonunit_epsilon_hb_delta=True,
        true_core_pressure_and_mean_inlets_retained=True,all_five_signed_physical_moments=True,
        surviving_leading_pressure_source_after_velocity_zero=True)


def original_flat_fixture():
    g=current.source.FunctionTransportGraph();s=g.symbol('flat_test_s');r=g.symbol('flat_test_rho');c=g.constant
    phase=sy.Rational(1,2);hb=sy.Rational(1,8)
    rows,geometry=current.original_cutoff_rows(g,s,r,dict(sc_quarter=c('1/2'),hb=c('1/8')))
    ctx=current.inner.collar.core.common.MPIntervalContext();ctx.dps=90
    ends=current.inner.collar.core.common.endpoints
    def evaluate(ref,frame,cache):
        if ref.node in cache:return cache[ref.node]
        n=g.nodes[ref.node];op=n['operation'];get=lambda i:evaluate(current.source.FunctionRef(g,i),frame,cache)
        if op=='exact_rational':v=ctx.mpf(n['numerator'])/n['denominator']
        elif op=='bound_variable':v=frame[n['name']]
        elif op=='sum':v=sum((get(i) for i in n['arguments']),ctx.mpf(0))
        elif op=='product':
            v=ctx.mpf(1)
            for i in n['arguments']:v*=get(i)
        elif op=='negative':v=-get(n['argument'])
        elif op=='positive_quotient':v=get(n['numerator'])/get(n['denominator'])
        elif op=='exact_original_source_flat_sigma':v=current.own.c2.flat_source.sigma_jets(ctx,get(n['argument']))[0]
        elif op=='analytic_unary':
            q=get(n['argument']);sigma=current.own.c2.flat_source.sigma_jets(ctx,q)
            assert n['name'] in ('original_flat_sigma_prime','original_flat_sigma_second')
            v=sigma[1] if n['name'].endswith('prime') else 2*sigma[2]
        else:raise AssertionError('Unexpected flat fixture operation '+op)
        cache[ref.node]=v;return v
    worst=mp.mpf(0)
    with mp.workdps(90):
        ph=mp.mpf(1)/2;b=mp.mpf(1)/8
        def cutoff(rho):
            q=4*mp.log(rho/4)/(ph*b)-2
            if q<=0:return mp.mpf(1)
            if q>=1:return mp.mpf(0)
            return 1-1/(1+mp.exp(1/q**2-1/(1-q)**2))
        for fraction in ('0','.5','.5625','.625','.75','1'):
            q=mp.mpf(fraction);sv=ph*q;rv=4*mp.exp(b*sv)
            frame={'flat_test_s':ctx.mpf(sv),'flat_test_rho':ctx.mpf(rv)}
            for i,row in enumerate(rows):
                box=evaluate(row.value,frame,{})
                if q<=mp.mpf('.5'):want=mp.mpf(1) if i==0 else mp.mpf(0)
                elif q>=mp.mpf('.75'):want=mp.mpf(0)
                else:want=mp.diff(cutoff,rv,i)
                lo,hi=ends(box);error=max(abs(lo-want),abs(hi-want));worst=max(worst,error)
                assert error<mp.mpf('1e-65'),(fraction,i,error)
        # q=.75 is the geometric mean in physical R, with the same epsilon.
        Rkeep=4*mp.exp(b*ph/2);Rzero=4*mp.exp(b*ph*3/4);Rin=4*mp.exp(b*ph)
        assert abs(Rzero**2-Rkeep*Rin)<mp.mpf('1e-80')
    return dict(passed=True,mathematical_fixture_only=True,original_flat_source_and_ordinary_rho_jets=True,
        ordinary_sigma_second_factor_two_checked=True,keep_zero_support_branches_checked=True,
        physical_zero_radius_is_geometric_mean=True,max_fixture_absolute_difference=str(worst))


def source_fixture(owner):
    assert set(owner.endpoint)==set(current.MOMENTS)
    assert owner.inner.definition['common_inner_endpoint_is_Rin_not_rminus']
    assert owner.definition['mean_K_Q_V_pressure_recovered_from_own_history']
    assert owner.definition['leading_pressure_forcing_unchanged']
    assert not owner.definition['full_five_defects_and_linear_repair_installed']
    g=owner.g;s=owner.inner.coordinates['first'];end=owner.parameters['sc_quarter']
    for jet in owner.endpoint.values():
        for ref in current.c3.rows(jet):
            n=g.nodes[ref.node]
            assert n['operation']=='function_substitution' and n['variable']==s.node and n['value']==end.node
    assert owner.fields['core']['rows']['mean1']==owner.radial.rows['core']['mean1']
    assert owner.fields['core']['rows']['K1']==owner.radial.rows['core']['K1']
    return dict(passed=True,true_regular_solution_and_axis_jets_used=True,
        actual_flat_cutoff_and_nonzero_original_sc_hb_used=True,
        five_signed_inner_endpoint_functions_not_caps=True,
        original_common_inner_Rin_scope_preserved=True,full_five_repair_not_claimed=True)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1InnerCutoffFunctions(require_checked=False)
    raw=current.parent.linear.read(current.NAME);owner.assert_graph()
    assert raw['exact_graph_nodes']==owner.g.nodes and raw['inner_cutoff_bindings']==owner.snapshot
    checks=dict(source=source_fixture(owner));print('TRUE_N1_CUTOFF_SOURCE_CHECKED',flush=True)
    checks['own_history']=own_history_fixture();print('SIGNED_INNER_CUTOFF_HISTORIES_CHECKED',flush=True)
    checks['original_flat']=original_flat_fixture();print('ORIGINAL_FLAT_RHO_JETS_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,**{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        original_numeric_n1_solution_oracle_executed=False,full_five_defects_and_linear_repair_installed=False,
        numerical_original_n1_point_accuracy_claimed=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_TRUE_INNER_CUTOFF_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()

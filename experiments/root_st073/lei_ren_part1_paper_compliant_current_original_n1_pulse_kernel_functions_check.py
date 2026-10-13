"""Finite-mu direct-integral and positive-branch C1 selector checks."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions as current


def evaluate(g,ref,bindings=None):
    bindings = bindings or {}
    cache = {}
    def at(i,env):
        if i in env:return env[i]
        key=(i,tuple(sorted((k,str(v)) for k,v in env.items())))
        if key in cache:return cache[key]
        row=g.nodes[i];op=row['operation'];child=lambda j:at(j,env)
        if op=='exact_rational':value=mp.mpf(row['numerator'])/row['denominator']
        elif op=='sum':value=sum(map(child,row['arguments']),mp.mpf(0))
        elif op=='product':value=mp.fprod(map(child,row['arguments']))
        elif op=='negative':value=-child(row['argument'])
        elif op=='positive_quotient':value=child(row['numerator'])/child(row['denominator'])
        elif op=='analytic_unary':
            x=child(row['argument'])
            value=(mp.expm1(x)/x if x else mp.mpf(1)) if row['name']=='exprel' else getattr(mp,row['name'])(x)
        elif op=='compact_raw_beta':
            x=child(row['argument']);value=mp.exp(-1/(1-x*x)) if abs(x)<1 else mp.mpf(0)
        elif op=='definite_integral':
            variable=next(j for j,r in enumerate(g.nodes) if r==dict(operation='bound_variable',name=row['variable']))
            lo,hi=child(row['lower']),child(row['upper'])
            points=[lo,0,hi] if lo<0<hi else [lo,hi]
            value=mp.quad(lambda x:at(row['integrand'],{**env,variable:x}),points)
        else:raise ValueError('No fixture oracle for '+op)
        cache[key]=value;return value
    return at(ref.node,bindings)


def direct_basis_fixture():
    with mp.workdps(50):
        g=current.source.FunctionTransportGraph();mu=mp.mpf('.07')
        W=current.kernel_functions(g,g.constant('7/100'))
        ell=mp.mpf('.15');beta=lambda r:mp.exp(-1/(1-r*r)) if abs(r)<1 else mp.mpf(0)
        norm=mp.quad(beta,[-1,0,1]);lam1,lam2=mp.mpf('.5')-mu,mp.mpf('.5')-2*mu
        A=[];B=[];D=[]
        for center in (-3,-1):
            direct=lambda l:mp.quad(lambda r:beta(r)*mp.exp(l*(center+ell*r)),[-1,0,1])/norm
            first,second=direct(lam1),direct(lam2)
            A.append(first);B.append(second);D.append((second-first)/mu)
        det=A[0]*D[1]-A[1]*D[0]
        assert det<0
        for name,want in (('first_row',A),('divided_row',D)):
            for ref,v in zip(W[name],want):assert abs(evaluate(g,ref)-v)<mp.mpf('1e-40')
        assert abs(evaluate(g,W['negative_divided_determinant'])-det)<mp.mpf('1e-40')
        gram=mp.quad(lambda r:beta(r)**2*mp.exp(-2*mu*ell*r),[-1,0,1])/(ell*norm**2)
        assert abs(evaluate(g,W['energy_gram'])-gram)<mp.mpf('1e-40')
        for ref,c in zip(W['original_end_energy_weights'],(-3,-1)):
            assert abs(evaluate(g,ref)-mp.exp(-26-2*mu*c)*gram)<mp.mpf('1e-40')
        J1,J2,m=sy.symbols('J1 J2 mu',positive=True)
        A1,A2=sy.exp(-3*(sy.Rational(1,2)-m))*J1,sy.exp(-(sy.Rational(1,2)-m))*J1
        B1,B2=sy.exp(-3*(sy.Rational(1,2)-2*m))*J2,sy.exp(-(sy.Rational(1,2)-2*m))*J2
        identity=A1*(B2-A2)/m-A2*(B1-A1)/m
        exact=-sy.exp(-2+5*m)*J1*J2*(sy.exp(2*m)-1)/m
        assert sy.simplify(identity-exact)==0
    return dict(passed=True,mathematical_fixture_only=True,original_mu_not_substituted=True,
        nonzero_fixture_mu='.07',ell='.15',direct_two_rows_and_divided_rows=True,
        negative_determinant_identity_for_every_positive_mu=True,nonzero_end_energy_weights_checked=True)


def selector_fixture():
    with mp.workdps(50):
        g=current.source.FunctionTransportGraph();z=g.symbol('Z')
        fixed=g.constant;Z=sy.Symbol('Z');b={str(Z):z}
        expression=lambda f:current.source.expression(g,sy.sympify(f),b)
        pair=lambda f:current.source.C1Function(expression(f),expression(sy.diff(f,Z)))
        W=dict(first_row=list(map(fixed,(1,2))),divided_row=list(map(fixed,(3,5))),
            positive_determinant_absolute=fixed(1),mu=fixed('7/100'),
            incoming_factors=[fixed('2/5'),fixed('3/10')],
            scaled_main_rows=[fixed('3/10'),fixed('1/2')],
            fixed_pulse_energy=fixed('1/10'),exact_nonzero_end_energy=[fixed('1/50'),fixed('1/40')])
        incoming=[pair(1+Z**2/10),pair(2-Z/9)];energy=pair(sy.Rational(1,5)+Z**2/5)
        # Independently solve the negative-determinant matrix and choose a
        # nonconstant positive amplitude. The future here is a mathematical
        # test function, explicitly not the project's full future energy.
        mu=sy.Rational(7,100);M1=sy.Rational(2,5)*(1+Z**2/10);M2=sy.Rational(3,10)*(2-Z/9)
        inverse=lambda r1,r2:[2*r2-5*r1,3*r1-r2]
        u=inverse(-mu*M1,-(M2-M1));v=inverse(-mu*sy.Rational(3,10),-sy.Rational(1,5))
        ap=1+Z/40;C=[p+ap*q for p,q in zip(u,v)]
        base=(1-sy.exp(-26))/4
        future=sy.Rational(1,10)*ap**2+C[0]**2/50+C[1]**2/40+mu*(sy.Rational(1,5)+Z**2/5)-base
        selected,_=current.select_functions(g,W,incoming,energy,pair(future))
        for value in ('-.6','0','.7'):
            env={z.node:mp.mpf(value)}
            for name,want in zip(('ap','C1','C2'),[ap]+C):
                for order,ref in enumerate((selected[name].value,selected[name].Z)):
                    expected=mp.mpf(str(sy.N(sy.diff(want,Z,order).subs(Z,sy.Rational(value)),48)))
                    assert abs(evaluate(g,ref,env)-expected)<mp.mpf('1e-40')
        try:current.select_functions(g,W,[1,2],energy,pair(future))
        except TypeError:pass
        else:raise AssertionError('Numeric incoming caps accepted')
    return dict(passed=True,mathematical_fixture_only=True,full_original_future_not_installed=True,
        nonconstant_positive_amplitude_and_two_signed_coefficients=True,
        true_first_Z_rows_checked=True,nonzero_end_energy_preserved=True,numeric_caps_rejected=True)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1PulseKernelFunctions(require_checked=False)
    owner.assert_graph()
    raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['exact_graph_nodes']==owner.g.nodes
    assert raw['exact_pulse_kernel_functions']==current.parent.encode(owner.functions)
    assert raw['actual_native_incoming_C1_functions']==current.parent.encode(owner.incoming)
    assert set(owner.incoming)=={'u','m1','m2','energy'}
    for p in owner.incoming.values():owner.g.ids([p.value,p.Z])
    assert raw['full_future_C1_function_still_required']
    assert not any(raw[k] for k in current.OPEN)
    fixtures=dict(direct_basis=direct_basis_fixture(),positive_selector=selector_fixture())
    fixture_frame=current.parent.linear.read(current.parent.quiet.FRAME_NAME)
    native=fixture_frame['actual_Rp_native_frame_functions']['actual_identified_Rp_native_frame_C3']
    original=current.parent.linear.realization.restore_graph(fixture_frame['exact_graph_nodes'])
    importer=current.parent.GraphImport(owner.g,original)
    expected={key:current.source.C1Function(*(importer(current.source.FunctionRef(original,i))
        for i in native[key][:2])) for key in owner.incoming}
    assert current.parent.encode(expected)==current.parent.encode(owner.incoming)
    fixtures['actual_native_incoming']=dict(passed=True,
        same_original_u_m1_m2_energy_C1_nodes=True,native_common_Pstar_normalization_preserved=True,
        genuine_functions_not_interval_constants=True)
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=fixtures,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_EXACT_PULSE_KERNEL_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()

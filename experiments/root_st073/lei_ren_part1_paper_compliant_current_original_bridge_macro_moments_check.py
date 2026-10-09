"""Independent full nonlinear six-history quadratures and source contracts."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_bridge_macro_moments as current
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields=current.fields;ep=current.ep


def finite(row):return row.coefficient*row.ctx.exp(row.scale.evaluate()) if not row.zero else row.ctx.mpf(0)


def contains(row,value):
    lo,hi=ep(finite(row));assert lo<=value<=hi,(mp.nstr(lo,20),mp.nstr(value,20),mp.nstr(hi,20))


def kernels():
    c=MPIntervalContext();c.dps=90
    logRa=c.ln(c.mpf('.7'));Y=c.ln(100)-logRa
    f=fields.MacroFlow(c,c.ln(c.mpf('.0002')),0,0,logRa,Y,c.mpf('.1'))
    poly=current.ExponentialPolynomial.constant(f,[f.scalar(1)]+[f.scalar(0)]*5)
    S=c.mpf('3.3');R0=f.scalar(c.mpf('.7'));R1=f.scalar(c.mpf('.7')*c.exp(S))
    count=primitives=0
    with mp.workdps(120):
        for rate in (1,2):
            for a in (-4,-2,-1,0,3):
                for k in (0,1,3,6):
                    expected=mp.quad(lambda t:mp.exp(-rate*(mp.mpf('3.3')-t)+a*t)*t**k,[0,mp.mpf('1.65'),mp.mpf('3.3')])
                    contains(poly.mass(rate,a,k,S,R0,R1),expected);count+=1
        for a in (-2,0,3):
            for k in (0,1,3,6):
                term=current.ExponentialPolynomial(f,{(a,k):[f.scalar(1)]+[f.scalar(0)]*5})
                result=term.primitive().evaluate(S,R0,R1)
                expected=mp.quad(lambda t:mp.exp(a*t)*t**k,[0,mp.mpf('1.65'),mp.mpf('3.3')])
                contains(result[0],expected);primitives+=1
    # The finite closed formula must satisfy the defining ODE/FTC at every
    # exponent/degree used by the production source polynomial algebra.
    x=sy.symbols('x',nonnegative=True);identities=0
    for rate in (1,2):
        for a in range(-6,9):
            for k in range(7):
                b=a+rate
                if b==0:F=sy.exp(-rate*x)*x**(k+1)/sy.Integer(k+1)
                else:
                    F=sy.exp(a*x)*sum((-1)**j*sy.factorial(k)*x**(k-j)/(sy.factorial(k-j)*sy.Integer(b)**(j+1)) for j in range(k+1))
                    F-=sy.exp(-rate*x)*(-1)**k*sy.factorial(k)/sy.Integer(b)**(k+1)
                assert sy.simplify(sy.diff(F,x)+rate*F-x**k*sy.exp(a*x))==0
                assert sy.simplify(F.subs(x,0))==0;identities+=1
    return dict(passed=True,independent_resonant_mass_quadratures=count,
        independent_zero_inlet_primitive_quadratures=primitives,exact_Volterra_ODE_and_zero_inlet_identities=identities)


def complete_nonlinear_histories():
    c=MPIntervalContext();c.dps=90
    z=sy.symbols('z');a=sy.symbols('a1:6')
    expression=sy.series(sy.exp(sum(a[n-1]*z**n for n in range(1,6))),z,0,6).removeO()
    bell=[sy.lambdify(a,expression.coeff(z,n),'mpmath') for n in range(6)]
    histories=radial=fieldchecks=0;records=[]
    with mp.workdps(120):
        def product(a,b):return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(6)]
        def plus(a,b):return [a[n]+b[n] for n in range(6)]
        def scale(a,v):return [x*v for x in a]
        def jet(values):return IntervalTaylor(c,[c.mpf(str(v)) for v in values])
        for sign,htext in ((1,'.0002'),(-1,'.0004')):
            h=mp.mpf(htext);Ra=mp.mpf('.7');logRa=c.ln(c.mpf('.7'));Y=c.ln(100)-logRa
            f=fields.MacroFlow(c,c.ln(c.mpf(htext)),c.mpf('.3'),c.mpf('-.2'),logRa,Y,c.mpf('.1'))
            d=[sign*mp.mpf('.006'),mp.mpf('.001'),mp.mpf('-.0007'),mp.mpf('.0004'),mp.mpf('-.0002'),mp.mpf('.0001')]
            q=list(map(mp.mpf,('1.1','.03','-.02','.01','-.004','.002')))
            ph=list(map(mp.mpf,('.8','.02','-.01','.005','-.002','.001')))
            v0=list(map(mp.mpf,('.3','.04','-.02','.01','-.004','.002')))
            vin=scale(list(map(mp.mpf,('.02','.001','-.0006','.0003','-.0001','.00005'))),h)
            lin=scale(list(map(mp.mpf,('.01','.002','-.001','.0005','-.0002','.0001'))),h)
            kh=list(map(mp.mpf,('.04','.004','-.002','.001','-.0004','.0002')))
            kp=list(map(mp.mpf,('.03','.003','-.0015','.0008','-.0003','.00015')))
            ks=list(map(mp.mpf,('.02','.002','-.001','.0005','-.0002','.0001')))
            zero=jet([0]*6)
            drives={'hydro':[jet(product(kh,d)),zero,zero],
                'pressure':[jet(product(kp,d)),zero,zero],
                'swirl':[jet(scale(product(ks,product(d,d)),h/4)),zero,zero]}
            f.set_sources([jet(d),zero,zero],drives,jet(q),f.jet(jet(lin)),jet(ph),jet(v0),f.jet(jet(vin)))
            initial={name:list(map(mp.mpf,('0.71','.01','-.006','.003','-.001','.0005'))) for name in current.RATES}
            owner=current.ActualMacroMoments(f,{name:f.jet(jet(row)) for name,row in initial.items()})
            eta=scale(d,h/2);etaR0=scale(eta,Ra*mp.exp(2*h));one=[mp.mpf(1)]+[mp.mpf(0)]*5
            qi=product(q,plus(kh,scale(kp,mp.exp(mp.mpf('.3')))));qs=scale(product(q,ks),mp.exp(mp.mpf('-.2')))
            ein=[mp.exp(lin[0])*fn(*lin[1:]) for fn in bell]
            Vin=plus(v0,vin)
            for fraction in ((1,2),(1,1)):
                result=owner.evaluate(fraction);R0=Ra*mp.exp(2*h)
                S=mp.mpf(fraction[0])/fraction[1]*(mp.log(100/Ra)-2*h)
                cache={}
                def sources(t):
                    key=t._mpf_
                    if key in cache:return cache[key]
                    R=R0*mp.exp(t);ell=plus(lin,scale(eta,-(R-R0)))
                    E=[mp.exp(ell[0])*fn(*ell[1:]) for fn in bell]
                    Phi=product(ph,E)
                    velocity=plus(Vin,scale(product(qi,plus(E,scale(ein,-1))),2))
                    swirl=plus(product(plus(one,etaR0),ein),scale(product(plus(one,scale(eta,R)),E),-1))
                    velocity=plus(velocity,scale(product(qs,swirl),-1))
                    src=dict(H=scale(Phi,2),M=velocity,K=scale(product(Phi,velocity),2),
                        A=product(velocity,velocity),B=product(Phi,Phi),C=product(Phi,Phi))
                    cache[key]=(src,Phi,velocity);return cache[key]
                cuts=[0,S/3,2*S/3,S];endpoint=sources(S)
                for name,rate in current.RATES.items():
                    for n in range(6):
                        integral=mp.quad(lambda t:mp.exp(-rate*(S-t))*sources(t)[0][name][n],cuts)
                        value=mp.exp(-rate*S)*initial[name][n]+integral
                        contains(result['actual_six_moment_functions'][name][n],value);histories+=1
                        derivative=endpoint[0][name][n]-rate*value
                        contains(result['actual_radial_ODE_derivative_functions'][name][n],derivative);radial+=1
                for key,index in (('phi',1),('V',2)):
                    for n in range(6):contains(result['same_original_macro_field_functions'][key][n],endpoint[index][n]);fieldchecks+=1
                records.append(dict(sign=sign,h=htext,fraction=list(fraction),complete_nonlinear_FV_sources_used=True,
                    nonzero_swirl_pressure_and_inlet=True,diagnostic_functions_only=True))
    return dict(passed=True,independent_complete_actual_history_Taylor_integrals=histories,
        independent_radial_ODE_derivative_comparisons=radial,independent_higher_expansion_field_comparisons=fieldchecks,
        records=records)


def native_checks(owner):
    c=owner.c;histories=derivatives=nonzero_errors=nonzero_decays=joins=0;degrees=set()
    for label in ('0','.5'):
        model,proof=owner.owner(label);f=model.flow
        assert proof['actual_G_not_Gbar_used'] and proof['original_micro_inlet_errors_retained']
        for q in ((0,1),(1,2),(1,1)):
            result=model.evaluate(q)
            assert result['original_volterra_ODEs_preserved'] and result['comparison_moments_not_substituted']
            assert result['analytic_P0_remains_separate'] and result['source_basis_and_micro_inlet_errors_retained']
            for name,rate in current.RATES.items():
                evidence=result['complete_integral_evidence'][name]
                assert evidence['rate']==rate and evidence['source_products_formed_before_enclosure']
                assert evidence['original_finite_interval_not_shortened']
                assert not evidence['original_incoming_decay'].zero;nonzero_decays+=1
                degrees.add(evidence['max_polynomial_degree'])
                if q[0]:
                    assert not evidence['complete_integral_error_norm'].zero;nonzero_errors+=1
                for row in result['actual_six_moment_functions'][name]:
                    assert row.scale.bases is f.logs and not row.record()['point_value_selected'];histories+=1
                for row in result['actual_radial_ODE_derivative_functions'][name]:
                    assert row.scale.bases is f.logs and not row.record()['point_value_selected'];derivatives+=1
                if not q[0]:
                    assert [v.record() for v in result['actual_six_moment_functions'][name]]==[v.record() for v in model.inlet[name]];joins+=1
            if q==(1,1):assert ep(finite(result['geometry']['R1']))==(100,100)
            assert result['kernel_coordinate']=='t=log(R/R0); dR/R=dt; original positive kernel applied once'
    return dict(passed=True,genuine_actual_history_Taylor_records=histories,genuine_radial_ODE_Taylor_records=derivatives,
        nonzero_complete_nonlinear_integral_error_norms=nonzero_errors,nonzero_actual_incoming_decay_factors=nonzero_decays,
        unchanged_actual_micro_inlet_joins=joins,observed_source_polynomial_degrees=sorted(degrees),
        B_source_phi_squared_with_rate2_preserved=True,
        micro_functions_whole_Z_switches_controls_and_recursion_remain_open=True)


def guards(owner):
    model,_=owner.owner('0');f=model.flow;count=0
    def reject(fn):
        nonlocal count
        try:fn()
        except (ValueError,ArithmeticError):count+=1;return
        raise AssertionError('Invalid moment source request accepted')
    reject(lambda:current.ActualMacroMoments(f,{}))
    reject(lambda:current.ExponentialPolynomial(f,{(0,-1):f.phi0}))
    reject(lambda:current.ExponentialPolynomial(f,{(True,0):f.phi0}))
    reject(lambda:current.ExponentialPolynomial(f,{(0,0):f.phi0[:4]}))
    foreign=owner.owner('.5')[0].flow
    reject(lambda:current.ExponentialPolynomial.constant(f,f.phi0)+current.ExponentialPolynomial.constant(foreign,foreign.phi0))
    for q in ((-1,1),(2,1),(1,0)):reject(lambda q=q:model.evaluate(q))
    reject(lambda:current.ExponentialPolynomial.constant(f,f.phi0).mass(3,0,0,1,f.scalar(1),f.scalar(2)))
    reject(lambda:owner.owner('anchor'))
    assert count==10
    return dict(passed=True,invalid_history_source_polynomial_rate_and_geometry_requests_rejected=count)


def run():
    began=time.monotonic();owner=current.OriginalBridgeMacroMoments()
    producer=json.loads((current.HERE/current.NAME).read_bytes())
    assert producer[current.GATE] and producer['source_family']==owner.family
    for name in fields.previous.OPEN:assert producer[name] is False
    for path,digest in producer['input_hashes'].items():assert current.sha(path)==digest
    mass=kernels();print('Exact polynomial primitive/Volterra kernels PASS',flush=True)
    refs=complete_nonlinear_histories();print('Independent complete nonlinear six-history integrals PASS',flush=True)
    native=native_checks(owner);bad=guards(owner)
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,exact_signed_polynomial_kernels=mass,
        independent_complete_nonlinear_actual_histories=refs,genuine_native_source_histories=native,guards=bad,
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Original actual macro six moment functions and radial ODEs on two native source frames, full nonlinear source errors and actual micro inlet retained. Full bridge, whole Z, switches and higher reconstruction remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original macro six histories and nonzero errors PASS',flush=True);return result


if __name__=='__main__':run()

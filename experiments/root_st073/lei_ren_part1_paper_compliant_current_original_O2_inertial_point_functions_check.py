"""Independent physical I/F and pressure-density jet checks for O2 functions.

Finite reference units and a manufactured pressure test differentiation only.
They never select the original source parameters or replace original P0.
"""
import json
from pathlib import Path
import time
import sympy as s
from sympy.core.function import AppliedUndef
import lei_ren_part1_paper_compliant_current_original_O2_inertial_point_functions as point
import lei_ren_part1_paper_compliant_pulse_end_physical_C2 as physical


def physical_reference(owner):
    t=owner.template;y,z,R,ps,delta=(t[k] for k in ('y','z','R','Pstar','delta'))
    E=t['E'];V=4*z;history=t['histories'];root=s.sqrt(2*R)
    moments=dict(theta=s.sqrt(2)*R**s.Rational(3,2)*ps*history['h'],
        z=R*history['m'],theta_z=s.sqrt(2)*R**s.Rational(3,2)*ps*history['k'],
        z_theta=R*ps**2*history['e'],p=ps**2*history['p'])
    P=ps**2*(t['P0']+history['p'])
    actual,hashes=physical.original_full_stress(dict(a=ps*E,b=V,
        ay=ps*s.diff(E,y),by=0,az=ps*s.diff(E,z),bz=4,dt=delta,z=z,R=R,
        root=root,L=1-delta*z*z,d=1-z*z,m=moments,
        mz={key:s.diff(value,z) for key,value in moments.items()},p=P,pz=s.diff(P,z)))
    checks={}
    for key,reference in (('p1',actual['Itheta']/(ps*E/root)),('p2',actual['Iz']/(ps*E/root))):
        for name,expected in ((key,reference),(key+'_Z',s.diff(reference,z))):
            assert s.cancel(s.expand(t[name]-expected))==0,name
            checks[name+'_equals_original_full_physical_I_over_F']=True
    assert t['p1'].has(R) and not t['p1'].has(ps,t['P0'])
    assert t['p2'].has(ps,t['P0'],s.diff(t['P0'],z))
    assert t['p2_Z'].has(s.diff(t['P0'],z,2))
    checks['absolute_pressure_and_required_P0_ZZ_not_discarded']=True
    hashes[Path(physical.__file__).name]=point.sha(Path(physical.__file__).name)
    return dict(passed=True,identities=checks,input_hashes=hashes)


def pressure_density_references(owner):
    operator=owner.pressure;p=operator.partition;z,t,q=p['z'],p['t'],p['q']
    comparisons=0
    assert len(operator.original_densities)==14
    assert not operator.raw_waiting_root.has(z)
    for name,density in operator.original_densities.items():
        assert all(not s.sympify(bound).has(z) for bound in p['domains'][name])
        if name in point.pressure.BETA2:
            factors=(-4*z/q,(-4+20*z*z)/q**2)
        elif name in point.pressure.BETA0:
            factors=(0,0)
        else:
            assert name=='z_flatten'
            sigma=p['sigma'](t/100)
            factors=(4*z*(sigma-1)/q,
                4*(sigma-1)/q+8*z*z*(sigma-1)*(2*sigma-3)/q**2)
        for order,factor in enumerate(factors,1):
            quotient=s.diff(density,z,order)/density-factor
            assert s.simplify(s.expand_log(s.expand_power_exp(quotient),force=True))==0,(name,order)
            comparisons+=1
    raw0=-sum((s.Integral(density,(t,*p['domains'][name]))
        for name,density in operator.original_densities.items()),s.Integer(0))
    assert owner.pressure_Z_jet(0)==raw0.subs(p['W'],operator.raw_waiting_root).xreplace({p['delta']:owner.template['delta']})
    def integral_rows(expression):
        # Differentiation may pull a Z-only coefficient outside an integral.
        # Compare the linear integral operators after lifting it back inside.
        rows={}
        for term in s.Add.make_args(s.expand(expression,deep=False)):
            factors=s.Mul.make_args(term)
            integrals=[factor for factor in factors if isinstance(factor,s.Integral)]
            if not integrals:
                assert term==0
                continue
            assert len(integrals)==1
            integral=integrals[0];coefficient=s.Mul(*(factor for factor in factors if factor is not integral))
            assert not coefficient.has(p['t'])
            rows[integral.limits]=rows.get(integral.limits,0)+coefficient*integral.function
        return rows
    for order in (1,2):
        differentiated=integral_rows(s.diff(owner.pressure_Z_jet(0),z,order))
        generated=integral_rows(owner.pressure_Z_jet(order))
        for limits in set(differentiated)|set(generated):
            difference=differentiated.get(limits,0)-generated.get(limits,0)
            assert s.simplify(s.expand_log(s.expand_power_exp(difference),force=True))==0,(order,limits)
    assert all(not owner.pressure_Z_jet(k).has(p['W']) for k in range(3))
    assert all(not owner.pressure_Z_jet(k).has(p['delta']) for k in range(3))
    assert owner.pressure_Z_jet(0).has(owner.template['delta'])
    return dict(passed=True,original_stage_count=14,independent_density_Z_jet_comparisons=comparisons,
        fixed_Z_independent_bounds_and_original_raw_waiting_root=True,
        same_delta_symbol_in_pressure_and_inertial_functions=True,
        normalized_pressure_sign_and_ordinary_integral_derivatives_preserved=True,
        pressure_enclosures_not_consumed_as_singleton_values=True)


def point_function_references(owner):
    t=owner.template;z=t['z'];got=owner.functions('.53',resolve_pressure=True)
    assert got['pressure_jet_derivatives_from_original_fourteen_stage_integral']
    assert got['p2'][0].has(s.Integral) and got['p2'][1].has(s.Integral)
    for pair in (got['p1'],got['p2']):
        for row in pair:
            assert not any(atom.func.__name__=='original_normalized_P0' for atom in row.atoms(AppliedUndef))
    c=owner.profiles.ctx
    sigma=point.profiles.loop.flat_step(c,owner.profiles.radial('.53')['y'])
    sign,mantissa,exponent,_=sigma._mpf_
    exact_binary=s.Integer((-1 if sign else 1)*mantissa)*s.Integer(2)**exponent
    assert got['a']==s.Rational(4,5)+s.Rational(6,5)*exact_binary
    assert got['b']==0 and got['approximate_radial_coefficients'] and not got['numerical_error_certified']
    compact=owner.functions('.53',resolve_pressure=False)
    # Manufactured pressure is an algebra/derivative probe, never original P0.
    probe=-3/(1+z*z)**2-s.Rational(1,1000)*(1+z**4)
    substitution={t['P0']:probe,s.diff(t['P0'],z):s.diff(probe,z),s.diff(t['P0'],z,2):s.diff(probe,z,2)}
    comparisons=0;maximum=s.Float(0)
    for R,ps,delta in ((2,7,s.Rational(1,10)),(5,13,s.Rational(1,100))):
        units={t['R']:R,t['Pstar']:ps,t['delta']:delta}
        for key in ('p1','p2'):
            value,derivative=(row.xreplace(substitution).subs(units) for row in compact[key])
            value_fn=s.lambdify(z,value,'mpmath');derivative_fn=s.lambdify(z,derivative,'mpmath')
            import mpmath as mp
            with mp.workdps(50):
                step=mp.mpf('1e-8')
                for coordinate in ('-.7','0','.37','.8'):
                    x=mp.mpf(coordinate)
                    finite_difference=(value_fn(x+step)-value_fn(x-step))/(2*step)
                    error=abs(finite_difference-derivative_fn(x));assert error<mp.mpf('2e-12'),(key,x,error)
                    maximum=max(maximum,s.Float(str(error),30));comparisons+=1
    rejected=0
    for function in (lambda:owner.pressure_at(2),lambda:owner.pressure_at(s.oo),
        lambda:owner.pressure_at(s.I),lambda:owner.pressure_Z_jet(-1),
        lambda:owner.pressure_Z_jet(3),lambda:owner.pressure_Z_jet(True),lambda:owner.functions(-1)):
        try:function()
        except ValueError:rejected+=1
    assert rejected==7
    return dict(passed=True,exact_original_integrals_in_returned_p2_and_Z=True,
        unresolved_P0_function_placeholders_eliminated=True,
        original_a_b_returned_separately=True,finite_difference_Z_comparisons=comparisons,
        maximum_reference_Z_error=str(maximum),reference_units_R_Pstar_delta=[[2,7,'.1'],[5,13,'.01']],
        manufactured_pressure_and_finite_units_are_reference_probes_only=True,
        domain_and_derivative_order_rejections=rejected)


def run():
    began=time.monotonic();owner=point.OriginalO2InertialPointFunctions()
    manifest=json.loads((point.HERE/point.NAME).read_bytes())
    assert manifest[point.GATE] and manifest['source_family']==owner.family
    for name,digest in manifest['input_hashes'].items():assert point.sha(name)==digest,name
    physical_checks=physical_reference(owner)
    flags=('original_p1_p2_scalar_point_values_installed','native_source_scale_and_phase_point_oracle_installed',
        'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed',
        'certified_actual_fixed_point_tail_installed','actual_terminal_Z_function_closure_installed',
        'current_whole_N_selected',*point.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    result=dict(all_passed=True,**{point.GATE:True},source_family=owner.family,
        independent_original_full_physical_reference=physical_checks,
        original_pressure_density_jet_references=pressure_density_references(owner),
        point_function_references=point_function_references(owner),
        **{key:False for key in flags},
        input_hashes={**manifest['input_hashes'],**physical_checks['input_hashes'],
            point.NAME:point.sha(point.NAME),Path(__file__).name:point.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original O2 full inertial expression/Z functions and exact original P0 jets only. Approximate radial coefficients, formal original scales and unevaluated exact pressure integrals remain; native numerical oracle, phase, certified errors, installed controls and global recursion are open.')
    (point.HERE/point.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Original O2 full inertial functions: physical I/F, pressure jets and point Z checks PASS',flush=True)
    return result


if __name__=='__main__':run()

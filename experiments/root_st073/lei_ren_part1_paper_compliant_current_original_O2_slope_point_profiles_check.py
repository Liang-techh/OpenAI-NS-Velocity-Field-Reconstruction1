"""Independent ODE integration checks for the original O2 point coefficients."""
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_O2_slope_point_profiles as point


def original_ODE_reference(y,z,Pstar,steps=1200):
    """Independent double RK4 of J and the five value/Z moment equations.

    Pstar is an explicit reference unit used only here to check each formal
    exponent and incoming history. It is not the native chosen parameter.
    """
    C=1/(1+z*z);CZ=-2*z*C*C;V=4*z/Pstar;VZ=4/Pstar
    h0=C*5/8;hz0=CZ*5/8
    state=[0,V,h0,V*h0,V*V-C*C*5/12,C*C*2.5,
        VZ,hz0,VZ*h0+V*hz0,2*V*VZ-C*CZ*5/6,5*C*CZ]
    def rhs(x,state):
        if x<=0:sigma=0
        elif x>=1:sigma=1
        else:
            odds=1/(1-x)**2-1/x**2
            e=math.exp(-abs(odds));sigma=e/(1+e) if odds<=0 else 1/(1+e)
        J,m,h,k,e,p,mz,hz,kz,ez,pz=state
        f=math.exp(x/10-.6*J);E=C*f;EZ=CZ*f
        return [sigma,V-m,E-1.5*h,E*V-1.5*k,V*V-E*E/2-e,E*E/2,
            VZ-mz,EZ-1.5*hz,EZ*V+E*VZ-1.5*kz,2*V*VZ-E*EZ-ez,E*EZ]
    dt=y/steps
    for i in range(steps):
        x=i*dt;k1=rhs(x,state);k2=rhs(x+dt/2,[v+dt*q/2 for v,q in zip(state,k1)])
        k3=rhs(x+dt/2,[v+dt*q/2 for v,q in zip(state,k2)]);k4=rhs(x+dt,[v+dt*q for v,q in zip(state,k3)])
        state=[v+dt*(a+2*b+2*c+d)/6 for v,a,b,c,d in zip(state,k1,k2,k3,k4)]
    return state


def references(owner,manifest):
    c=owner.ctx;errors=[];points=[];dispatch_count=0
    graph=json.loads((point.HERE/point.ALLN).read_bytes())['exact_function_graph_nodes']
    rows={r['function_role']:r for r in graph if r.get('chart')=='O2_slope' and r.get('function_role') in owner.roles}
    assert len(rows)==4
    assert owner.leaf_identity['passed']
    before=len(owner.radial_cache)
    owner.background(Z='.37',y='.271')
    assert len(owner.radial_cache)==before
    evaluate=lambda v,p:sum(coefficient*c.mpf(p)**power for power,coefficient in v.terms)
    for y in ('.137','.53','.91'):
        for z in ('-.6','0','.7'):
            got=owner.evaluate(Z=z,y=y)
            for pstar in (7,13):
                ref=original_ODE_reference(float(y),float(z),pstar)
                assert abs(float(got['original_defining_integrals']['J'])-ref[0])<2e-11
                for i,key in enumerate(('m','h','k','e','p')):
                    for order,expected in (('value',ref[i+1]),('Z',ref[i+6])):
                        actual=evaluate(getattr(got['histories'][key],order),pstar)
                        error=abs(actual-c.mpf(expected));assert error<c.mpf('2e-11'),(y,z,pstar,key,order,error)
                        errors.append(error)
            for role,row in rows.items():
                actual=owner.dispatch_original_background(row,coordinate=y,Z=z)
                expected=getattr(got['E' if '_E_' in role else 'V'],'Z' if role.endswith('_Z') else 'value')
                assert actual==expected;dispatch_count+=1
            C=1/(1+c.mpf(z)**2)
            assert abs(got['E'].Z.unscaled_scalar()+2*c.mpf(z)*C*got['E'].value.unscaled_scalar())<c.mpf('1e-45')
            assert got['a'].Z.unscaled_scalar()==got['b'].value.unscaled_scalar()==got['b'].Z.unscaled_scalar()==0
            assert got['absolute_pressure']['original_analytic_P0_reference']['datum_enclosure_sha256']==owner.family['datum_enclosure_sha256']
            assert not got['absolute_pressure']['pressure_datum_reset']
            points.append(dict(y=got['y'],Z=got['Z'],E=got['E'].record(),V=got['V'].record(),
                original_point_coefficients=True,quadrature_errors_not_certified=True))
    # Both original source endpoints; same reference inlet and slope=-1/2 outlet.
    left=owner.evaluate(Z='.37',y=0);right=owner.evaluate(Z='.37',y=1)
    assert left['a'].value.unscaled_scalar()==c.mpf('.8')
    assert right['a'].value.unscaled_scalar()==2 and right['log_E_y']==-c.mpf('.5')
    assert owner.J(0)==0 and owner.J(1)==c.mpf('.5')
    assert left['histories']['h'].value.unscaled_scalar()==left['E'].value.unscaled_scalar()*c.mpf('.625')
    rejected=0
    for function in (lambda:owner.evaluate(Z=2,y='.5'),lambda:owner.evaluate(Z=0,y=-1),
        lambda:owner.evaluate(Z='nan',y='.5'),lambda:left['V'].Z.unscaled_scalar(),
        lambda:owner.dispatch_original_background({**rows['all_N_original_E_C0'],'source_node':-1},coordinate='.5',Z=0),
        lambda:owner.dispatch_original_background({**rows['all_N_original_E_C0'],'function_role':'all_N_periodic_A_C0'},coordinate='.5',Z=0)):
        try:function()
        except ValueError:rejected+=1
    assert rejected==6
    assert not manifest['quadrature_and_roundoff_certified'] and not manifest['full_original_O2_source_point_provider_installed']
    return dict(passed=True,native_original_profile_point_queries=len(points),
        independent_moment_value_Z_ODE_reference_comparisons=len(errors),
        independent_ODE_reference_Pstar_units=[7,13],reference_units_not_native_parameter_selection=True,
        maximum_absolute_independent_ODE_reference_error=max(errors),original_E_V_role_dispatches=dispatch_count,
        original_profile_endpoints_checked=2,wrong_source_scope_domain_or_scalar_cast_rejected=rejected,
        same_original_P0_reference_and_nonzero_axial_factor_retained=True,point_queries=points)


def run():
    began=time.monotonic();owner=point.OriginalO2SlopePointProfiles()
    manifest=json.loads((point.HERE/point.NAME).read_bytes())
    assert manifest[point.GATE] and manifest['source_family']==owner.family
    for name,digest in manifest['input_hashes'].items():assert point.sha(name)==digest,name
    result=dict(all_passed=True,**{point.GATE:True},source_family=owner.family,
        original_recipe_binding=owner.recipe,independent_original_profile_references=references(owner,manifest),
        original_E_V_graph_leaf_identity=owner.leaf_identity,
        **{key:manifest[key] for key in ('original_p1_p2_and_phase_provider_installed',
            'full_original_O2_source_point_provider_installed','quadrature_and_roundoff_certified',
            'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed',
            'certified_actual_fixed_point_tail_installed','actual_terminal_Z_function_closure_installed','current_whole_N_selected',*point.loop.OPEN)},
        input_hashes={**manifest['input_hashes'],point.NAME:point.sha(point.NAME),Path(__file__).name:point.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original O2 profile coefficients and E/V background source roles only. Independent finite-unit ODE references test factor powers and histories; original global Pstar, phase, p1/p2, P0 point functions, certified errors and field are not installed.')
    (point.HERE/point.RECEIPT).write_text(json.dumps(point.loop.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original O2 point profiles: independent moment ODE/units/Z and source-role dispatch PASS',flush=True)
    return result


if __name__=='__main__':run()

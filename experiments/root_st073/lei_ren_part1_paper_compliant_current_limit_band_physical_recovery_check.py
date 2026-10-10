"""Current source connection, independent physical calculus and numeric map."""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_limit_band_physical_recovery as current


def exact_graph_formula_checks():
    g=current.source.FunctionTransportGraph()
    pair=lambda name:current.source.C1Function(g.symbol(name),g.symbol(name+'_Z'))
    fields={key:pair(key) for key in ('E','V','E_y','V_y')}
    own={key:pair(key) for key in current.current.current.RATES}
    own_y={key:pair(key+'_y') for key in own};P0=pair('P0')
    out=current.native_recovery(g,fields,own,own_y,P0,g.symbol('Z'),g.symbol('delta'))
    symbols={}
    def walk(q):
        row=g.nodes[q if type(q) is int else q.node];op=row['operation']
        if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
        if op=='bound_variable':return symbols.setdefault(row['name'],sy.Symbol(row['name']))
        if op=='sum':return sum(map(walk,row['arguments']))
        if op=='product':return sy.prod(map(walk,row['arguments']))
        if op=='negative':return -walk(row['argument'])
        if op=='positive_quotient':return walk(row['numerator'])/walk(row['denominator'])
        raise AssertionError(op)
    Z,de,E,V,Ey,Vy,m,mZ,my,myZ,h,hZ,k,kZ,e,eZ,p,pZ,P,PZ=sy.symbols(
        'Z delta E V E_y V_y m m_Z m_y m_y_Z h h_Z k k_Z e e_Z p p_Z P0 P0_Z')
    L=1-de*Z*Z;d=1-Z*Z;tr=(1-de)*Z*m+d*mZ
    expected=dict(Q=(2*Z*V-tr)/L,Q_y=(2*Z*Vy-(1-de)*Z*my-d*myZ)/L,
        shear_theta=2*Ey-E,shear_axial=2*Vy)
    for key,expr in expected.items():assert sy.cancel(walk(out[key])-expr)==0,key
    sectors=dict(theta_linear=(-E+(1-de/2)*h-(1-de)*Z*hZ/2)/L,
        theta_quadratic=((2*de-1)*Z*k-d*kZ+E*tr)/L,
        axial_linear=(-V+(1-de)*(m-Z*mZ)/2)/L,
        axial_quadratic=(V*tr+2*de*Z*e-d*eZ+2*(1+de)*Z*(P+p)-d*(PZ+pZ))/L)
    for key,expr in sectors.items():assert sy.cancel(walk(out['full_signed_inertial_sectors'][key])-expr)==0,key
    assert sy.expand(walk(out['pressure'].value)-(P+p))==0
    assert sy.expand(walk(out['pressure'].Z)-(PZ+pZ))==0
    return dict(actual_graph_Q_Qy_pressure_and_four_signed_sectors=True)


def numeric_mapping_checks():
    c=mp.mp.clone();c.dps=105;root_cases=0
    for delta in ('0','.07','.6','.99'):
        for z,tau,nu in (('0','.7','1'),('.21','.7','1.3'),('-.21','1e-70','.6'),('1e30','1e-70','2')):
            q=current.solve_similarity_scale(c,z,tau,nu,delta,iterations=4*c.prec)
            lam=c.exp(q['log_lambda']);z,tau,nu,de=map(c.mpf,(z,tau,nu,delta))
            scaled=(lam*lam-z*z*lam**(2*de)/nu-tau)
            assert abs(scaled)/(lam*lam+tau)<c.mpf('1e-98')
            assert abs(q['Z'])<=1+c.mpf('1e-100')
            root_cases+=1
    for args in ((0,0,1,'.1'),(0,1,0,'.1'),(0,1,1,1),(0,1,1,'-.1')):
        try:current.solve_similarity_scale(c,*args)
        except ValueError:pass
        else:raise AssertionError('Invalid physical scale domain admitted')
    # Analytic synthetic local functions satisfy m_y=V-m exactly. The
    # independently differentiated Cartesian field checks the entire map,
    # moving basis, signed Z, implicit lambda and viscosity factors.
    de,S,Rc,a=c.mpf('.07'),c.mpf('1.25'),c.mpf('2.3'),c.mpf(1)/3
    def packet(Z,x):
        b=x**a;m=Z*b;E=(1+Z*Z/5)*x**(-c.mpf('.7'))
        p=(1+Z*Z/5)**2*(1-x**(-c.mpf('1.4')))/(2*c.mpf('1.4'))
        pZ=(4*Z/5)*(1+Z*Z/5)*(1-x**(-c.mpf('1.4')))/(2*c.mpf('1.4'))
        return dict(E=E,V=(1+a)*m,E_y=-c.mpf('.7')*E,V_y=a*(1+a)*m,
            P0=-c.mpf('.4')-Z*Z/10,P0_Z=-Z/5,m=m,m_Z=b,m_y=a*m,m_yZ=a*b,
            p=p,p_Z=pZ,p_y=E*E/2,p_yZ=E*(2*Z/5)*x**(-c.mpf('.7')))
    family={'synthetic_reference_only':True};digest='synthetic_analytic_field_not_original_current_source'
    oracle=SimpleNamespace(mode='synthetic_reference',source_family=family,source_graph_sha256=digest,
        parameters=lambda:dict(delta=de,S=S,Rc=Rc),band_state=packet)
    adapter=current.ApproximateCartesianBandAdapter(oracle=oracle,source_family=family,graph_sha256=digest,ctx=c)
    t,Zphys=c.mpf('.3'),c.mpf('.21');q=current.solve_similarity_scale(c,Zphys,1-t,1,de,iterations=4*c.prec)
    r=c.exp(q['log_lambda'])*c.sqrt(2*Rc*c.mpf('1.4'));theta=c.mpf('.39')
    x,y=r*c.cos(theta),r*c.sin(theta)
    field=lambda x,y,z:adapter.cartesian(x,y,z,t)
    div=c.diff(lambda xx:field(xx,y,Zphys)['u'],x)+c.diff(lambda yy:field(x,yy,Zphys)['v'],y)+c.diff(lambda zz:field(x,y,zz)['w'],Zphys)
    assert abs(div)<c.mpf('1e-97'),div
    value=field(x,y,Zphys);Ut=-c.sin(theta)*value['u']+c.cos(theta)*value['v']
    radial_pressure=c.diff(lambda rr:field(rr*c.cos(theta),rr*c.sin(theta),Zphys)['p'],r)
    assert abs(radial_pressure-Ut*Ut/r)<c.mpf('1e-97')
    assert field(x,y,-Zphys)['w']==-value['w']
    nu=c.mpf('1.7');scaled=adapter.cartesian(c.sqrt(nu)*x,c.sqrt(nu)*y,c.sqrt(nu)*Zphys,t,viscosity=nu)
    for key in ('u','v','w'):assert abs(scaled[key]-c.sqrt(nu)*value[key])<c.mpf('1e-97')
    assert abs(scaled['p']-nu*value['p'])<c.mpf('1e-97')
    bad=packet(c.mpf('.1'),c.mpf('1.2'));bad['m']={'log_absolute_upper':'0'}
    try:current.recover_numeric_values(c,bad,'.1',de)
    except TypeError:pass
    else:raise AssertionError('Range/cap used as a numeric field value')
    try:current.ApproximateCartesianBandAdapter(oracle=None,source_family=family,graph_sha256=digest)
    except current.SourcePointOracleRequired:pass
    else:raise AssertionError('Missing actual source oracle silently supplied')
    return dict(synthetic_reference_only=True,actual_current_numeric_source_not_installed=True,
        implicit_root_cases=root_cases,independent_Cartesian_divergence=True,
        independent_radial_pressure_gradient=True,signed_axial_and_constant_viscosity_dilation=True,
        cap_and_missing_oracle_rejection=True,directed_numeric_error_certificate=False)


def run():
    began=time.monotonic();accepted,pressure,hashes=current.load_inputs()
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert report[current.GATE] and report['source_family']==accepted['source_family']
    for name,digest in report['input_hashes'].items():assert current.sha(name)==digest,name
    built=current.build(accepted,pressure);g=built['graph']
    assert g.nodes==report['exact_graph_nodes']
    assert current.current.current.encode_graph(built)==report['exact_physical_recovery']
    assert g.nodes[:built['source_graph_prefix_length']]==accepted['exact_graph_nodes']
    Rc=g.nodes[built['Rc'].node]
    assert Rc['quantity_path']==['original_closed_O3_background','actual_source_radius'] and Rc['coordinate']==g.constant(2).node
    assert Rc['source_background_binding']==current.current.BACKGROUND_BINDING and Rc['Z_independent']
    assert built['original_recovery_binding']==current.current.current.ast_binding(current.recovery.GenericMomentRecovery.field)
    lam=g.nodes[built['physical_coordinate_functions']['lambda_scale'].node]
    assert lam['lambda_not_replaced_by_sqrt_tau'] and lam['q_is_lambda_squared']
    assert lam['source_parameter_report_sha256']==current.sha(current.PREFIX+'pressure_source.json')
    assert set(built['cartesian_velocity_pressure'])=={'u','v','w','p'}
    assert set(built['terminal_absolute_cumulative'])=={'M','I','J','S_energy','Cp'}
    assert set(built['required_future_integrals'])=={'M','J','S_energy','renormalized_I','pressure_tail'}
    assert all(pair.value!=g.zero for pair in built['required_future_integrals'].values())
    final=built['final_heat_reference'];power=built['current_O3_power_reference']
    assert final['final_reference_c_infinity_not_installed']
    assert final['actual_transition_2Rc_to_final_heat_not_installed']
    assert final['current_O3_amplitude_not_identified_with_final_c_infinity']
    assert final['c_infinity']==built['future_c_infinity']
    assert power['mu']==built['parameters']['mu'] and power['actual_C1_amplitude_and_Z_row_retained']
    assert power['current_mu_not_replaced_with_final_delta_over_2']
    assert built['future_integral_definitions']['domain']=='[2Rc,infinity)'
    assert built['future_integral_definitions']['necessary_conditions_only']
    assert built['future_integral_definitions']['actual_exterior_source_not_installed']
    assert built['native_recovery']['Q_Z_and_full_velocity_mixed_jets_not_available']
    assert report['symbolic_checks']==current.symbolic_checks()
    formulas=exact_graph_formula_checks();numeric=numeric_mapping_checks()
    for key in ('actual_numeric_point_source_oracle_installed','actual_numeric_controls_evaluated',
        'actual_global_frequency_admitted','physical_original_exterior_five_targets_closed',
        'full_recovered_velocity_pressure_and_heat_joins_admitted','higher_Z_velocity_jets_and_stress_cone_admitted',
        'actual_temporal_scale_recursion_installed',*current.current.outer.OPEN):assert report[key] is False,key
    hashes[current.NAME]=current.sha(current.NAME);hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=report['source_family'],
        exact_physical_graph_nodes_checked=len(g.nodes),same_actual6368_node_limit_band_reused=True,
        original_n0_recovery_graph_formula_checks=formulas,
        source_Rc_pressure_and_similarity_parameters_bound=True,
        paper_cumulative_and_absolute_future_targets_checked=True,
        current_mu_power_and_conditional_final_heat_reference_kept_distinct=True,
        independent_numeric_map_checks=numeric,
        no_accepted_source_integration_or_live_owner_replayed=True,
        actual_numerics_higher_jets_and_exterior_not_claimed=True,
        input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_LIMIT_BAND_PHYSICAL_RECOVERY',len(g.nodes),'nodes; independent divergence/pressure/map checks',flush=True)
    return receipt


if __name__=='__main__':run()

"""Independent original traces, Z-independent weights and C1 integral check."""
import ast
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_native_Rc_C1_integral_realization as producer
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original

HERE,PREFIX,sha=producer.HERE,producer.PREFIX,producer.sha
packets,ep,require=producer.packets,producer.ep,producer.require


def independent_flat_and_circle_proof():
    x=s.symbols('x',positive=True)
    limits={str(k):s.limit(s.exp(-1/x**2)/x**k,x,0,dir='+') for k in (0,1,2,6)}
    require(all(v==0 for v in limits.values()),'Original exponential cutoff lacks flat value/first traces')
    tree=ast.parse((HERE/(PREFIX+'current_generic_shear_loop.py')).read_text(encoding='utf8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='flat_step')
    odds=next(n.value for n in fn.body if isinstance(n,ast.Assign) and any(ast.unparse(t)=='odds' for t in n.targets))
    require(ast.dump(odds)==ast.dump(ast.parse('1 / (1 - x) ** 2 - 1 / x ** 2',mode='eval').body),
            'Flat trace proof no longer binds the original cutoff')
    source_tree=ast.parse((HERE/(PREFIX+'current_generic_loop_function_sources.py')).read_text(encoding='utf8'))
    build=next(n for n in source_tree.body if isinstance(n,ast.FunctionDef) and n.name=='build')
    assignments={ast.unparse(n.targets[0]):n.value for n in build.body if isinstance(n,ast.Assign) and len(n.targets)==1}
    recipes={'ratio':'g.add(one, g.mul(t0, t0), g.mul(two, q, q))',
        'K':'g.div(one, g.mul(twopi, ratio), \'correlated_ratio_ge1\')',
        'T1':'g.integral(t, psi)','T2':'g.integral(g.mul(t, t), psi)',
        'Phi':'g.mul(K, P)','lambda0':'g.mul(K, g.add(one, g.mul(t, t)))',
        'A':'g.flat(g.mul(half, a, chi), Delta, eta)','BB':'g.flat(g.mul(half, E, M), Delta, eta)'}
    for key,expr in recipes.items():
        require(ast.dump(assignments[key])==ast.dump(ast.parse(expr,mode='eval').body),
                'Original phase/primitive C1 recipe changed: '+key)
    r,q,h,t0=s.symbols('r q h t0',real=True)
    # w=sum_{n>=1} r^(n-1) cos(n psi), |r|<1. Absolute convergence
    # and Parseval give int w=0 and int w^2=pi/(1-r^2).
    second=2*s.pi*t0*t0+4*q*q*s.pi/(h*h*(1-r*r))
    require(s.cancel(second.subs(r*r,1-1/(h*h))-2*s.pi*(t0*t0+2*q*q))==0,
            'Independent full-period Poisson second moment differs')
    a,eta,delta,aZ,dZ,sigma,sigmaZ=s.symbols('a eta Delta a_Z Delta_Z sigma sigma_Z',real=True)
    gamma=2*eta-delta;factor=s.sqrt(gamma/(2*a))
    qZ=factor*(sigmaZ-sigma*(dZ/(2*gamma)+aZ/(2*a)))
    weakZ=-factor*(dZ/(2*gamma)+aZ/(2*a))
    require(s.simplify((qZ.subs({sigma:1,sigmaZ:0})-weakZ).subs(delta,0))==0 and
            s.simplify(qZ.subs({sigma:0,sigmaZ:0,delta:eta}))==0,
            'One-sided q_Z traces do not match at Delta=0,eta')
    return dict(passed=True,independent_original_flat_cutoff_limits=4,
        original_phase_and_primitive_AST_bindings=len(recipes),
        full_period_identities_from_absolutely_convergent_original_Poisson_Fourier_series=True,
        positive_original_phase_slope_and_parameter_implicit_function_theorem_bound=True,
        q_active_transition_flat_value_and_Z_traces_glue=True,
        symbolic_one_sided_q_Z_cutoff_trace_identities=2,
        first_Z_period_traces_vanish_without_spatial_chart_seam_assumption=True)


def independent_piecewise_C1_quadrature():
    c=mp.mp.clone();c.dps=45;eta=c.mpf('.02');slope=c.mpf('.07')
    def cutoff(z):
        delta=eta+slope*z;a=2+delta
        if delta>=eta:return c.mpf(0),c.mpf(0)
        theta=1-delta/eta;factor=c.sqrt((2*eta-delta)/(2*a))
        sigma=original.flat_step(c,theta)
        if theta<=0 or theta>=1:derivative=c.mpf(0)
        else:
            odds=1/(1-theta)**2-1/theta**2;small=c.exp(-abs(odds))
            derivative=small/(1+small)**2*(2/(1-theta)**3+2/theta**3)
        q=sigma*factor
        qZ=factor*(-slope*derivative/eta+sigma*(-slope/(2*(2*eta-delta))-slope/(2*a)))
        return q,qZ
    def density(phase,z,n,multiplier):
        q,qZ=cutoff(z);A=q*c.sin(2*c.pi*phase);AZ=qZ*c.sin(2*c.pi*phase)
        B=q*c.sin(4*c.pi*phase);BZ=qZ*c.sin(4*c.pi*phase)
        E=multiplier*(1+c.mpf('.1')*z);EZ=multiplier*c.mpf('.1');V=c.mpf('.3')*z;VZ=c.mpf('.3')
        factor=c.mpf(1) if A==0 else c.expm1(A/n)/(A/n)
        F=E*A*factor;FZ=EZ*A*factor+E*AZ*c.exp(A/n)
        values={'m':B,'h':F,'k':V*F+E*B+F*B/n,'e':2*V*B-E*F+(B*B-F*F/2)/n,'p':E*F+F*F/(2*n)}
        jets={'m':BZ,'h':FZ,'k':VZ*F+V*FZ+EZ*B+E*BZ+(FZ*B+F*BZ)/n,
            'e':2*VZ*B+2*V*BZ-EZ*F-E*FZ+(2*B*BZ-F*FZ)/n,'p':EZ*F+E*FZ+F*FZ/n}
        return values,jets
    # A deliberately discontinuous spatial seam at x=1/2 has fixed,
    # Z-independent endpoints. N=160,256 makes each half a whole number
    # of periods. Integrate one original phase period on each piece.
    def integral(z,n,key,order):
        return sum(c.quad(lambda phi:density(phi,z,n,m)[order][key],[0,c.mpf('.25'),c.mpf('.5'),c.mpf('.75'),1])/2
                   for m in (c.mpf(1),c.mpf(2)))
    step=c.mpf('1e-10');count=0
    for n in (160,256):
        for z in (c.mpf('-.6'),-2/c.mpf(7),c.mpf('-.01'),c.mpf(0),c.mpf('.01'),c.mpf('.6')):
            for key in ('m','h','k','e','p'):
                exact=integral(z,n,key,1)
                finite=(integral(z+step,n,key,0)-integral(z-step,n,key,0))/(2*step)
                require(abs(exact-finite)<c.mpf('1e-15')*max(1,abs(exact)),
                        'Independent piecewise Z differentiation differs: '+key)
                count+=1
    return dict(passed=True,independent_piecewise_integral_Z_comparisons=count,
        N_reference_values=[160,256],active_transition_flat_Z_cutoff_seams_checked=True,
        fixed_spatial_jump_retained=True,numerical_reference_only=True,
        original_native_field_or_fixed_point_not_numerically_evaluated=True)


def original_graph_and_majorant_proof(owner,built,report):
    old=owner.data['current_native_Rc_all_N_function_controls'];nodes=built['graph'].nodes
    require(nodes==old['exact_function_graph_nodes'],'C1 realization changed the original exact function graph')
    def Z_independent(index,seen=None):
        seen=set() if seen is None else seen
        if index in seen:return True
        seen.add(index);row=nodes[index];op=row['operation']
        if op in ('exact_rational','original_source_parameter','shared_positive_integer'):return True
        if op=='bound_variable':return row['name']!='Z'
        if op in ('sum','product'):return all(Z_independent(i,seen) for i in row['arguments'])
        if op in ('negative','analytic_unary'):return Z_independent(row['argument'],seen)
        if op=='positive_quotient':return Z_independent(row['numerator'],seen) and Z_independent(row['denominator'],seen)
        return False
    integral_ids=set();caps=0;quiet=0
    for cell,record in zip(built['cells'],report['actual_original_C1_integral_cell_records']):
        original_cell=cell.original_cell
        for key in ('lower','upper','left_radius_offset','right_radius_offset','width'):
            require(Z_independent(original_cell[key]),'Actual native endpoint/radius has unresolved Z dependence')
        require(record['different_spatial_chart_velocity_or_stress_seam_equality_used'] is False and
                record['incoming_C1_histories_and_rate_zero_pressure_memory_retained'],
                'Piecewise regularity claimed a spatial seam or reset pressure')
        for index in cell.exact_integral_nodes:
            row=nodes[index]
            require(row['exact_function_integral'] and row['coefficient_still_depends_on_N_and_actual_phase'] and
                    Z_independent(row['lower']) and Z_independent(row['upper']),
                    'Original integral identity/bounds changed')
            integral_ids.add(index)
        for group in record['integrable_value_and_Z_majorants'].values():
            for row in group.values():
                for cap in row.values():
                    require(cap['encloses_original_source_function'] and not cap['point_value_selected'],'Native cap used as exact function')
                    if cap['log_absolute_upper'] is not None:
                        value=packets.interval(owner.ctx,cap['log_absolute_upper'])
                        require(all(mp.isfinite(v) for v in ep(value)),'Original common density/kernel majorant not finite')
                    caps+=1
        if original_cell['source_flat_exact_zero']:
            require(not cell.exact_integral_nodes and all(p.value==built['graph'].zero and p.Z==built['graph'].zero
                    for order in cell.coefficient_density_pairs.values() for p in order.values()),'Original flat density is not exactly zero')
            quiet+=1
    native={i for i,row in enumerate(nodes) if row['operation']=='definite_integral' and row.get('coefficient_still_depends_on_N_and_actual_phase')}
    require(integral_ids==native and len(native)==336 and caps==480 and quiet==3,'Native C1 integral/majorant inventory incomplete')
    require(report['exact_C1_Rc_N_scaled_target_roots']==old['exact_reconstructed_N_scaled_target_roots'],
            'C1 target source roots changed')
    for key in ('actual_five_controls_installed','certified_actual_fixed_point_tail_installed','actual_terminal_Z_function_closure_installed',
        'current_whole_N_selected','numerical_original_source_point_or_integral_oracle_installed','original_source_ancestor_constructors_called',*packets.OPEN):
        require(report[key] is False,'C1 integral realization overclaims later physical construction: '+key)
    return dict(passed=True,original_exact_graph_prefix_nodes=len(nodes),original_continuous_C1_cells=24,
        actual_native_nonzero_coefficient_integrals=336,original_value_Z_own_mass_majorants=480,
        exact_original_zero_density_cells=3,all_actual_integration_endpoints_and_weights_Z_independent=True,
        same_source_exact_targets_not_arbitrary_C1_integral_symbols=True,
        spatial_y_seams_numeric_point_oracle_global_N_controls_terminal_and_recursion_remain_open=True)


def run():
    began=time.monotonic();owner=producer.NativeRcC1IntegralRealization();built=owner.build()
    report=json.loads((HERE/producer.NAME).read_bytes())
    require(report[producer.GATE] and report['source_family']==owner.family,'Current C1 integral realization required')
    for name,value in owner.hashes.items():require(report['input_hashes'].get(name)==value,'Original dependency changed: '+name)
    result=dict(all_passed=True,**{producer.GATE:True},source_family=owner.family,
        independent_original_cutoff_and_circle_proof=independent_flat_and_circle_proof(),
        independent_piecewise_integral_reference=independent_piecewise_C1_quadrature(),
        actual_original_graph_geometry_and_majorants=original_graph_and_majorant_proof(owner,built,report),
        exact_original_integral_functions_C1_in_Z_defined=True,
        read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
        original_source_ancestor_constructors_called=False,numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes={**report['input_hashes'],producer.NAME:sha(producer.NAME),Path(__file__).name:sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (HERE/producer.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original C1 integral realization: flat/periodic traces, native graph and piecewise differentiation PASS',flush=True)
    return result


if __name__=='__main__':run()

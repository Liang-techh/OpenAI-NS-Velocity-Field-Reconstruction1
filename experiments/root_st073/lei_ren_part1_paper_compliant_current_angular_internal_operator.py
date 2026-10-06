"""Two-sided full angular source/operator traces, before any interval bounds."""
import ast
import functools
import math
from types import SimpleNamespace,FunctionType
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import shifted_rows
from lei_ren_part1_paper_compliant_current_pulse_gap_background_tensor import canonical_tensor_groups
from lei_ren_part1_paper_compliant_current_pulse_support_interfaces import weighted_support_FTC_theorem


@functools.lru_cache(maxsize=1)
def two_sided_full_angular_trace_theorem():
    """Replay both complete nonzero sides after source-bound FTC limits.

    The prior weighted source theorem supplies beta jets zero and shared
    full/empty cumulative values. Inside beta and remaining kernels are
    kept symbolic until that exact limit, independently of the outside input.
    """
    asts=SourceAST();v,z,a,mu=s.symbols('source_s source_Z a mu',real=True)
    support_limit=weighted_support_FTC_theorem()
    if not support_limit['passed'] or not all(support_limit['identities'].values()):
        raise ValueError('Original weighted FTC and flat beta limits required')
    asts.method('current_pulse_support_interfaces','weighted_support_FTC_theorem')
    asts.method('current_angular_support_interfaces','source_edge_integrals')
    delta=2*a;beta=-2-delta
    c=SimpleNamespace(mpf=s.sympify,exp=s.exp,ln=s.log,cos=s.cos,sin=s.sin)
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,delta=delta,k=1-a,prate=1+delta,
        S=s.Symbol('same_exact_inverse_R',positive=True))
    outer=SimpleNamespace(mu=mu,prate=1+2*mu)
    Xpast=s.Function('same_nonzero_flatten_plus_accumulated_angular_trace')(z)
    Epost=s.Function('same_nonzero_complete_post_energy')(z)
    Ppost=s.Function('same_nonzero_complete_post_pressure_with_P0')(z)
    JE=s.Function('same_exact_remaining_E_F_trace')(z)
    JP=s.Function('same_exact_remaining_B_D_trace')(z)
    inside_beta=s.Function('inside_local_normalized_beta')(v,z)
    dj=s.Function('same_selected_C5_coefficient')(z)
    inside_JE=s.Function('inside_original_remaining_E_F')(v,z)
    inside_JP=s.Function('inside_original_remaining_B_D')(v,z)
    stub=SimpleNamespace(constant=lambda ctx,value,order:s.sympify(value),variable=lambda ctx,value,order:value)
    def decay(ctx,p,length):return (1-s.exp(-p*length))/p
    env=dict(IntervalTaylor=stub,copy_jet=lambda ctx,value:value,math=math,
        product_rows=product_rows,decay_integral=decay)
    moment=asts.replay('current_angular_background_stress','normalized_full_moment_rows',env)
    F=1+dj*inside_beta
    native_inside=dict(swirl_factor_one_plus_h_Taylor=F,angular_Taylor=Xpast/F,
        actual_angular_bump_y_derivatives=[s.diff(F,v,j) for j in range(5)])
    native_outside=dict(swirl_factor_one_plus_h_Taylor=s.Integer(1),angular_Taylor=Xpast,
        actual_angular_bump_y_derivatives=[s.Integer(1)]+[s.Integer(0)]*4)
    inside=moment(heat,outer,native_inside,dict(post=Epost),dict(energy=inside_JE,pressure=inside_JP),
        dict(normalized_full_post_pressure=Ppost),v,c)
    outside=moment(heat,outer,native_outside,dict(post=Epost),dict(energy=JE,pressure=JP),
        dict(normalized_full_post_pressure=Ppost),v,c)
    def exact_limit(expression):
        replacements={inside_beta:s.Integer(0),inside_JE:JE,inside_JP:JP}
        for derivative in expression.atoms(s.Derivative):
            orders=dict(derivative.variable_count)
            if derivative.expr==inside_beta:replacements[derivative]=s.Integer(0)
            elif derivative.expr in (inside_JE,inside_JP):
                trace=JE if derivative.expr==inside_JE else JP
                replacements[derivative]=s.Integer(0) if orders.get(v,0) else s.diff(trace,z,orders.get(z,0))
        return s.cancel(expression.xreplace(replacements))
    left={name:[exact_limit(value) for value in rows] for name,rows in inside.items()}
    right={name:[s.cancel(value) for value in rows] for name,rows in outside.items()}
    moments={}
    for name in ('A','E','P','K'):
        for j in range(5):
            residual=s.cancel(left[name][j]-right[name][j])
            if residual!=0:raise ValueError('Two-sided full angular normalized moment source differs')
            for k in range(5-j):
                if s.diff(residual,z,k)!=0:raise ValueError('Two-sided full mixed4 source trace differs')
                moments[name+'/y'+str(j)+'_Z'+str(k)]=True
    # Feed both complete nonzero sources into the actual original stress AST.
    stress_env=dict(IntervalTaylor=stub,axial_derivative=lambda value:s.diff(value,z),
        shifted_rows=shifted_rows)
    stress_program=asts.replay('collar_stress_C3','collar_stress_rows',stress_env)
    class SourceJet:
        def __init__(self,value):self.value=value
        def __getitem__(self,n):return s.diff(self.value,z,n)/math.factorial(n)
    grid_env=dict(math=math)
    ordinary_grid=asts.replay('pulse_end_physical_C2','ordinary_grid',grid_env)
    from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
    from lei_ren_part1_paper_compliant_collar_physical_C2 import collar_velocity_bracket
    # Replay the physical loops and Cartesian completion verbatim. Capture
    # signed defining fields, not logarithmic upper enclosures or source caps.
    original=asts.method('current_angular_background_stress','angular')
    begin=next(i for i,node in enumerate(original.body) if isinstance(node,ast.Assign) and ast.unparse(node.targets[0])=='ps')
    end=next(i for i,node in enumerate(original.body) if isinstance(node,ast.Assign) and ast.unparse(node.targets[0])=='pressure_rows')
    statements=original.body[begin:end]
    if not statements or not any(isinstance(node,ast.For) for node in statements):raise ValueError('Actual full physical angular loops changed')
    lam,nu=s.symbols('same_lambda same_nu',positive=True);lt,lr,angle=s.symbols('same_log_tau same_log_R same_angle',real=True)
    parts={name:dict(original_same_log=s.Symbol('actual_'+name+'_log',real=True)) for name in ('B','Qtheta','Qz','completed_diagonal')}
    actual_exponent=asts.expression('collar_physical_C2','physical_source_row','exponent')
    actual_combined=asts.expression('collar_physical_C2','physical_source_row','combined')
    asts.expression('collar_physical_C2','scale_row',"result['signed_coefficient']",wanted="row['signed_coefficient']*multiplier")
    def signed_row(ctx,coefficient,logs,gamma,logtau,viscosity,N,radial_log=0,nu_base=1):
        local=dict(c=ctx,parts=logs,spatial_order=N,radial_log=radial_log,nu_base=nu_base)
        exponent=asts.evaluate(actual_exponent,local);combined=asts.evaluate(actual_combined,local)
        return dict(formal_signed_field=coefficient*s.exp(combined)*lam**gamma*viscosity**exponent)
    def signed_scale(ctx,row,factor):return dict(formal_signed_field=factor*row['formal_signed_field'])
    def physical(rows):
        stress=stress_program(heat,dict(K_rows=rows['K']),dict(angular_defect_rows=rows['A'],
            energy_defect_rows=rows['E'],pressure_defect_rows=rows['P'],K_defect_rows=rows['K']),z,s.Integer(0))
        grids={name:ordinary_grid([SourceJet(value) for value in stress[name]],3) for name in ('theta','axial','theta_inertial','theta_shear')}
        K=[SourceJet(value) for value in rows['K']]
        local=dict(self=SimpleNamespace(heat=heat),c=c,Z=z,theta=angle,nu=nu,lt=lt,logR=lr,
            beta=beta,parts=parts,grids=grids,K=K,math=math,physical_bracket=physical_bracket,
            collar_velocity_bracket=collar_velocity_bracket,physical_source_row=signed_row,scale_row=signed_scale)
        exec(compile(ast.fix_missing_locations(ast.Module(body=statements,type_ignores=[])),
            '<actual-full-angular-two-sided-physical-source>','exec'),local)
        view=dict(physical_cylindrical_stress_mixed3={name:dict(angular_full=grid) for name,grid in local['ps'].items()},
            physical_cylindrical_stress_divergence_mixed2={name:dict(angular_full=grid) for name,grid in local['div'].items()},
            completed_theta_theta_stress_mixed2=dict(angular_full=local['diagonal']),
            physical_completed_stress_tensor_cartesian=local['tensor'],
            physical_completed_stress_divergence_cartesian=local['divcart'],physical_remainder_cartesian=local['ecart'],
            physical_momentum_residual_decomposition_cartesian=local['residual'])
        zero={key:signed_scale(c,row,s.Integer(0)) for key,row in local['error'].items()}
        view['physical_three_component_remainder_mixed2']=dict(radial=dict(exact_zero=zero),
            theta=dict(angular_full=local['error']),axial=dict(exact_zero=zero))
        return canonical_tensor_groups(view)
    lhs=physical(left);rhs=physical(right);checks={};nonzero=0
    if set(lhs)!=set(rhs) or len(lhs)!=71:raise ValueError('Complete both-sided canonical tensor source required')
    for key in lhs:
        x=sum((row['formal_signed_field'] for row in lhs[key]),s.Integer(0))
        y=sum((row['formal_signed_field'] for row in rhs[key]),s.Integer(0))
        residual=s.cancel(x-y)
        if residual!=0:raise ValueError('Two-sided full physical angular tensor source differs: '+key)
        nonzero+=int(x!=0)
        checks[key]=True
    if not nonzero:raise ValueError('Full source replaced by zero support difference')
    for stem,method in (('pulse_physical_bounds','physical_bracket'),('collar_physical_C2','collar_velocity_bracket'),
        ('collar_physical_C2','physical_source_row'),('collar_physical_C2','scale_row'),
        ('collar_Gamma_C4','product_rows')):asts.method(stem,method)
    return dict(original_two_sided_full_moment_mixed4_identities=moments,
        full_moment_identity_count=len(moments),original_two_sided_full_physical_canonical_group_identities=checks,
        full_physical_group_identity_count=len(checks),nonzero_common_formal_group_count=nonzero,
        same_original_moment_stress_and_physical_loop_AST_replayed_on_both_complete_sides=True,
        retained_nonzero_X_post_energy_post_pressure_P0_and_remaining_quadratic_histories=True,
        exact_inside_beta_and_remaining_weight_FTC_limits_before_operators=True,
        original_lambda_nu_radius_amplitude_and_Cartesian_factors_identical=True,
        formal_signed_source_not_interval_upper_or_inverse_radius_cap=True,
        checked_original_weighted_FTC_and_flat_beta_limits=support_limit,
        input_hashes=asts.hashes,passed=True)

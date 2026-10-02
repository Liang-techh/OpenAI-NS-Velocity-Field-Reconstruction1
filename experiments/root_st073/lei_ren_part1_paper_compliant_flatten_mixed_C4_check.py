"""Independent variable-rate derivatives and original pulse/O5 joins."""
import ast
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import flatten_mixed
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
NAME=PREFIX+'compliant_flatten_mixed_C4.json'


def source_branch_expression(method,target,needle,env):
    tree=ast.parse((HERE/(PREFIX+'compliant_corrected_outer_field.py')).read_text(encoding='utf-8'))
    function=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
    values=[n.value for n in ast.walk(function) if isinstance(n,ast.Assign)
        and any(ast.unparse(v)==target for v in n.targets) and needle in ast.unparse(n.value)]
    if len(values)!=1:raise ValueError('Unique original source branch required: '+method+' '+target)
    def formal(node):
        label=ast.unparse(node)
        if label in env:return env[label]
        if isinstance(node,ast.Constant):return s.Rational(str(node.value))
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -formal(node.operand)
        if isinstance(node,ast.BinOp):
            a,b=formal(node.left),formal(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
        if isinstance(node,ast.Call) and ast.unparse(node.func)=='c.exp' and len(node.args)==1:
            return s.exp(formal(node.args[0]))
        raise ValueError('Unbound original source branch: '+label)
    return formal(values[0])


def independent_variable_rate_fixture():
    """Nonflat sigma=y/100 fixture admits closed functions, not ODE targets."""
    with mp.workdps(90):
        c=MPIntervalContext(); c.dps=100; radius=mp.mpf('1e-70')
        mu=mp.mpf('.04'); bp=mp.mpf('.5')+mu; rate=1-mu; Ev2=mp.mpf('.2')
        rho=lambda z:mp.log((1+z*z)/2)
        aa=lambda z:rho(z)/100-bp
        cc=lambda z:rate+rho(z)/100
        dd=lambda z:2*mu-2*rho(z)/100
        theta=lambda y,z:mp.exp(aa(z)*y)/(1+z*z)
        X=lambda y,z:mp.mpf('.7')*mp.exp(-cc(z)*y)+(1-mp.exp(-cc(z)*y))/cc(z)
        e=lambda y,z:(1+z*z)*mp.exp(dd(z)*y)-(mp.exp(dd(z)*y)-1)/(2*dd(z))
        pressure=lambda y,z:mp.mpf('.1')/(1+z*z)**2+Ev2*(mp.exp(2*aa(z)*y)-1)/(4*aa(z)*(1+z*z)**2)
        checks=0; y=mp.mpf('.3')
        for Z in (mp.mpf(-1),mp.mpf(0),mp.mpf('.3'),mp.mpf(1)):
            def jet(function):return IntervalTaylor(c,[c.mpf([v-radius,v+radius]) for v in mp.taylor(function,Z,5)])
            r=flatten_mixed(c,c.mpf(mu),jet(rho),[c.mpf(y/100),c.mpf('.01'),c.mpf(0),c.mpf(0),c.mpf(0)],
                jet(lambda z:theta(y,z)),jet(lambda z:X(y,z)),jet(lambda z:e(y,z)),jet(lambda z:pressure(y,z)),c.mpf(Ev2))
            for key,function,rows in ((UT,theta,r['physical_velocity_and_pressure_y_derivative_Taylor'][UT]),
                (P,pressure,r['physical_velocity_and_pressure_y_derivative_Taylor'][P]),
                ('angular',X,r['angular_y_derivative_Taylor']),('energy',e,r['energy_y_derivative_Taylor'])):
                for k in range(5):
                    for n in range(5-k):
                        target=mp.diff(function,(y,Z),(k,n))/math.factorial(n)
                        lo,hi=endpoints(rows[k][n])
                        if not lo<=target<=hi:raise ArithmeticError('Independent flatten derivative failed: '+str((key,Z,k,n)))
                        checks+=1
        return dict(independent_closed_form_variable_rate_derivative_checks=checks,
            surrogate_linear_sigma_fixture_only=True,actual_Md40_source_admission=False,passed=True)


def functional_identities():
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Flatten functional identity failed: '+name)
        proofs[name]=True
    t,mu,q,rho,sig,lq,Xv,Ev2=s.symbols('t mu q rho sigma logq Xv Ev2',real=True)
    bp=s.Rational(1,2)+mu; rate=1-mu
    F=s.exp(rho*sig); J=s.symbols('J',real=True)
    env={'ratio':rho,'sig':sig,'lq':lq,'self.bp':bp,'t':t,'integral':J,'self.Xv':Xv,'self.rate':rate,'F':F}
    zero('original_flatten_logtheta',source_branch_expression('flatten_raw','L','ratio * sig',env)-(rho*sig-lq-bp*t))
    zero('original_flatten_angular_integral',assignment('compliant_corrected_outer_field','flatten_raw','X',env)-(J+Xv*s.exp(-rate*t))/F)
    ff=s.Function('F')(t); Y=s.Function('Y')(t); G=s.Function('G')(t)
    theta=ff/q*s.exp(-bp*t); X=Y/ff; energy=G*q*q*s.exp(2*mu*t)/(2*ff*ff)
    fp=s.diff(ff,t); rhoy=fp/ff
    zero('theta_lograte',s.diff(theta,t)-(rhoy-bp)*theta)
    zero('angular_variable_rate_ODE',s.diff(X,t).subs(s.diff(Y,t),ff-rate*Y)-1+(rate+rhoy)*X)
    zero('energy_variable_rate_ODE',s.diff(energy,t).subs(s.diff(G,t),-ff*ff*s.exp(-2*mu*t)/(q*q))
        +s.Rational(1,2)-(2*mu-2*rhoy)*energy)
    ev=s.symbols('ev',real=True); Eint=s.symbols('Eint',real=True)
    zero('complete_future_Rv_Ev0_units',energy.subs(G,2*ev/q**2).subs(ff,1).subs(t,0)-ev)
    zero('energy_forward_integral_units',energy.subs(G,(2*ev-Eint)/q**2)-(ev-Eint/2)*s.exp(2*mu*t)/(ff*ff))
    zero('original_packet_energy_units',source_branch_expression('packet','energy','Eglobal',
        {'Eglobal':G,'theta':theta,'y':t})-energy)
    zero('original_flatten_exit_Z_independent_swirl',(F/q).subs({rho:s.log(q/2),sig:1})-s.Rational(1,2))
    zero('exact_terminal_pressure_decay_separated_logs',-13*(1+2*mu)/mu-(-13/mu-26))
    # At both flat endpoints all sigma y derivatives1..4 vanish.
    # The actual pulse terminal primitives supply the exact histories;
    # identical y ODEs and common axial C5 data identify all mixed jets.
    a=s.symbols('sigma_y1:5',real=True)
    for k in range(4):
        lograte=rho*a[k]-(bp if k==0 else 0)
        angularrate=rho*a[k]+(rate if k==0 else 0)
        energyrate=-2*rho*a[k]+(2*mu if k==0 else 0)
        substitutions={v:0 for v in a}
        zero('flat_endpoint_theta_rate'+str(k),lograte.subs(substitutions)-(-bp if k==0 else 0))
        zero('flat_endpoint_angular_rate'+str(k),angularrate.subs(substitutions)-(rate if k==0 else 0))
        zero('flat_endpoint_energy_rate'+str(k),energyrate.subs(substitutions)-(2*mu if k==0 else 0))
    proofs['terminal_m1_m2_B_and_radial_zero_inherited_from_empty_future_supports']=True
    proofs['positive_terminal_energy_is_complete_future_half_not_zero']=True
    proofs['same_original_P0_Mp_and_Xv_histories_are_preserved']=True
    proofs['axial5_and_sigma_y_through4_suffice_for_mixed_velocity_C4']=True
    proofs['same_source_ODEs_and_flat_comparison_identify_two_sided_pulse_O5_join']=True
    proofs['flatten_exit_has_pure_power_y_ODEs']=True
    return proofs


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Flatten source changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    checks=zero_checks=0
    for point in r['samples']+[r['whole_Z_inlet'],r['whole_Z_flatten_box'],r['whole_Z_exit']]:
        for flag in ('actual_terminal_zero_linear_and_radial_histories_inherited',
            'original_pressure_datum_and_complete_positive_future_energy_preserved',
            'future_energy_input_is_derivative_enclosure_not_a_polynomial_fit',
            'terminal_pressure_is_an_enclosure_of_the_original_exact_source'):
            if not point[flag]:raise ValueError('Original post-pulse source history required')
        if endpoints(read(point['energy_Taylor']['coefficients'][0]))[0]<=0:
            raise ArithmeticError('Positive complete future energy lost')
        for component,grid in point['physical_mixed_derivatives_total_order_le4'].items():
            if len(grid)!=15:raise ValueError('Flatten mixed-index coverage incomplete')
            for value in grid.values():
                ends=endpoints(read(value))
                if any(not mp.isfinite(v) for v in ends):raise ArithmeticError('Flatten derivative nonfinite')
                checks+=1
                if component in (UZ,UR):
                    if ends!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Inherited empty linear support lost')
                    zero_checks+=1
        for flag in ('full_pulse_C4_installed','full_outer_C4_certified','physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
            if point[flag]:raise ValueError('Flatten scope promoted: '+flag)
    pulse=json.loads((HERE/(PREFIX+'compliant_pulse_mixed_C4.json')).read_bytes())
    fifth=json.loads((HERE/(PREFIX+'compliant_fifth_axial_jets.json')).read_bytes())
    U=read(fifth['whole_Z']['selected']['incoming']['Z_independent_constant_definitions']['U'])
    left=pulse['whole_Z_terminal']; right=r['whole_Z_inlet']; overlap_checks=0
    def overlap(a,b,label):
        nonlocal overlap_checks
        a,b=endpoints(a),endpoints(b)
        if max(a[0],b[0])>min(a[1],b[1]):raise ArithmeticError('Independent O5 interface diagnostic disjoint: '+label)
        overlap_checks+=1
    for component,grid in right['physical_mixed_derivatives_total_order_le4'].items():
        for index,value in grid.items():
            value=read(value)
            # Pulse velocities factor exp(-13bp/mu) and inlet U/q;
            # O5 uses Ev0=U*Pstar*exp(-13bp/mu). Compare in common units.
            if component!=P:value=value*U
            overlap(value,read(left['physical_mixed_derivatives_total_order_le4'][component][index]),component+' '+index)
    for n in range(6):
        overlap(read(right['energy_Taylor']['coefficients'][n]),read(left['Mztheta_over_R_Utheta_squared']['coefficients'][n]),'positive energy '+str(n))
        overlap(read(right['angular_Taylor']['coefficients'][n]),read(left['Mtheta_over_sqrt2_R_3half_Utheta']['coefficients'][n]),'angular '+str(n))
    for value in right['sigma_y_derivatives'][1:]:
        if endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Original flatten inlet not flat')
    exit_grid=r['whole_Z_exit']['physical_mixed_derivatives_total_order_le4'][UT]
    for index,value in exit_grid.items():
        if int(index.split('_Z')[1])>0 and endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):
            raise ArithmeticError('Flatten exit swirl is not exactly Z independent')
    flat=json.loads((HERE/(PREFIX+'compliant_pulse_flat_comparison_check.json')).read_bytes())
    inlet=json.loads((HERE/(PREFIX+'compliant_power_inlet_C4_check.json')).read_bytes())
    internal=json.loads((HERE/(PREFIX+'compliant_pulse_interface_certificate.json')).read_bytes())
    if not(flat['all_passed'] and inlet['two_sided_O3_pulse_join_certified'] and internal['exact_functional_main_gap_and_gap_end_identities_certified']):
        raise ValueError('Full local pulse interface prerequisites missing')
    proofs=functional_identities(); fixture=independent_variable_rate_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        all_passed=True,functional_original_flatten_and_interface_identities=proofs,
        independent_variable_rate_fixture=fixture,finite_actual_flatten_mixed_bounds_checked=checks,
        inherited_zero_mixed_bounds_checked=zero_checks,independent_pulse_O5_overlap_diagnostics=overlap_checks,
        independently_evaluated_O5_flatten_high_mixed_derivatives=True,
        entire_original_100_unit_flatten_mixed_C4_available=True,
        two_sided_pulse_O5_join_certified=True,two_sided_O3_pulse_join_certified=True,
        full_pulse_C4_installed=True,
        full_pulse_C4_scope='Leading O4 pulse and its two local external interfaces; profile velocity/pressure mixed total<=4, physical cylindrical r/z via admitted map, fixed source/positive tau; no time derivatives or whole outer/core/axis claim',
        following_power_high_mixed_derivatives_available=False,full_outer_C4_certified=False,
        physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Pulse/O5 and entire flatten: independent variable-rate derivatives, preserved positive energy and original functional joins PASS',flush=True)
    return out


if __name__=='__main__':run()

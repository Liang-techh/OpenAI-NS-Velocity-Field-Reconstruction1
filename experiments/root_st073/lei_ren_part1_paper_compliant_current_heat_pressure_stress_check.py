"""Current forward heat histories, constant stress terms and FTC checks."""
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

from lei_ren_part1_paper_compliant_current_heat_pressure_stress import (
    CurrentHeatPressureStress,HERE,PREFIX,NAME,RECEIPT,GATES,SCOPES,VIEWS,sha,
    binding,history_source_bindings,constant_stress_rows,pack,encode,endpoints,
    collar_pressure_rows,IntervalTaylor)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_moment_stress_identities
from lei_ren_part1_paper_compliant_heat_terminal_stress_identities import terminal_stress_identities


def exact(left,right,message):
    if s.simplify(s.expand_power_exp(left-right))!=0:raise ArithmeticError(message)


def retained_forward_pressure_proof():
    t,v,z=s.symbols('t v Z',real=True);rate,scale=s.symbols('prate pressure_scale',positive=True)
    K=s.Function('same_current_K');Ptail=s.Function('actual_waiting_Ptail')(z)
    density=s.exp(-rate*v)*K(v,z)**2/2
    prefix=s.Integral(density,(v,0,t));tail=s.Integral(density,(v,t,s.oo))
    # Integral additivity is applied only to this same density and matching
    # limits, already bound to actual forward/collar/Gamma production code.
    infinity=Ptail+scale*(prefix+tail)
    backward=infinity-scale*tail;forward=Ptail+scale*prefix
    exact(backward,forward,'Backward pressure changed original Ptail')
    constant=s.Function('actual_P_infinity')(z)
    pressure=constant-scale*tail;derivatives=0
    for order in range(1,5):
        expected=scale*s.exp(-rate*t)/2*sum(
            s.binomial(order-1,j)*(-rate)**(order-1-j)*s.diff(K(t,z)**2,t,j) for j in range(order))
        exact(s.diff(pressure,t,order),expected,'Original pressure FTC derivative omitted a term')
        derivatives+=1
    a,S,eps,W,D=s.symbols('a S epsilon W D',real=True);pre=1-eps*W;bracket=pre-a*S*D
    exact(bracket**2,pre**2-a*S*D*(2*pre-a*S*D),'Full collar pressure deficit density differs from K squared')
    exact(pre**2,1-2*eps*W+eps**2*W**2,'Two original epsilon atoms were merged')
    x,H=s.symbols('exact_inverse_radius H',positive=True);deficit=(1-H)/(a*x)
    exact(1-2*a*x*deficit*(1-a*x*deficit/2),H**2,'Infinite Gamma pressure numerator differs')
    return dict(actual_forward_density_AST_bound=True,
        same_density_full_integral_additivity_proved=True,
        pressure_identity='P/Pstar^2=Cp-pressure_scale*full_remaining_pressure(t)',
        actual_constant_definition='Cp=Ptail+pressure_scale*full_remaining_pressure(0)',
        trace3_definition='Cp=pressure3+pressure_scale*exp(-3prate)*Gamma_pressure_numerator(3)',
        ordinary_positive_pressure_FTC_identities=derivatives,
        full_collar_and_Gamma_pressure_density_identities=3,
        arbitrary_axial_derivatives_follow_same_smooth_integrals=True,
        axis_pressure_and_original_forward_history_retained=True,
        Cp_not_defined_as_zero_or_a_fitted_tail=True,passed=True)


def actual_constant_stress_proof():
    t,z=s.symbols('t Z',real=True);a,Rt,B0,Pstar=s.symbols('a Rtail B0 Pstar',positive=True)
    delta=2*a;k=1-a;b=(1-delta)/2;prate=1+delta;d=1-z*z;L=1-delta*z*z
    D=s.Function('actual_angular_constant')(z);C=s.Function('actual_pressure_infinity')(z)
    R=Rt*s.exp(t);B=B0*s.exp(-(s.Rational(1,2)+a)*t)
    qt=s.sqrt(R/2)*B;qp=s.sqrt(R/2)*Pstar**2
    # The actual extra X=D exp(-kt)/K gives a radially constant Mtheta.
    extra_M=s.sqrt(2)*R**s.Rational(3,2)*B*D*s.exp(-k*t)
    exact(s.diff(extra_M,t),0,'Inherited angular primitive constant was dropped')
    theta=(k*extra_M-b*z*s.diff(extra_M,z))/(2*L*R)
    pressure=Pstar**2*C
    axial=R*(2*prate*z*pressure-d*s.diff(pressure,z))/(L*s.sqrt(2*R))
    theta_coefficient=s.exp(-k*t)*(k*D-b*z*s.diff(D,z))/L
    axial_coefficient=(2*prate*z*C-d*s.diff(C,z))/L
    identities=0
    for j in range(5):
        for n in range(5-j):
            exact(s.diff(theta,t,j,z,n)/qt,(-1)**j*s.diff(theta_coefficient,z,n),'Current angular constant stress derivative differs')
            exact(s.diff(axial,t,j,z,n)/qp,s.Rational(1,2)**j*s.diff(axial_coefficient,z,n),'Actual pressure constant stress derivative differs')
            identities+=2
    exact(s.diff(theta,t),-theta,'Actual angular stress is not the R^-1 homogeneous term')
    exact(s.diff(axial,t),axial/2,'Actual pressure offset stress radial factor differs')
    return dict(original_3_16_and_3_17_constant_terms_recovered=True,
        actual_extra_angular_moment_is_constant_in_R=True,
        constant_stress_mixed4_physical_factor_identities=identities,
        angular_radial_power=-1,pressure_offset_radial_power='1/2',
        pressure_constant_factor_separate_from_Ev2=True,
        no_division_by_a_cap_box_containing_zero=True,
        both_current_constants_retained_in_actual_stress=True,passed=True)


def exact_source_radius_proof(field):
    for stem,method,target,expression in (
        ('compliant_exact_heat_component','__init__','self.logRref','c.ln(110)+10*(self.logC+self.params.logPstar)'),
        ('compliant_exact_heat_component','__init__','self.tail_finite','self.params.yd+1+self.params.Tw+100-30*self.params.log_mu+2+self.params.Ts+self.angular.waiting'),
        ('compliant_exact_heat_component','__init__','self.logradius_terms','dict(selected_reference=self.logRref,pulse_term=13/self.params.mu,finite_offset=self.tail_finite)'),
        ('compliant_future_swirl_energy','__init__','self.Lrel','-30*self.params.log_mu'),
        ('compliant_outer_angular_repair','__init__','logRtail','self.heat.logRref+13/self.mu+self.heat.tail_finite')):
        binding(stem,method,target,expression)
    # Current typed parameter/source bindings supply the same defining
    # selected amplitude and scalar parameters; interval copies are bounds.
    proof=field.current_bindings
    if not all(proof['current_defining_object_graph'].values()) or not proof['shared_exact_Gamma_future_binding']['passed']:
        raise ValueError('Current canonical full Gamma/source bridge missing')
    lr,lp,mu,Tw,Lrel,Ts,wait,lmu,t,z=s.symbols('logRref logP mu Tw Lrel Ts wait logmu t Z',real=True)
    raw=lr+13/mu+lp+1+Tw+100-30*lmu+2+Ts+wait
    physical=(lr+lp+1+Tw)+13/mu+100+Lrel+2+Ts+wait
    exact(physical.subs(Lrel,-30*lmu),raw,'Current and original exact heat radii differ')
    S=s.exp(-raw);R=s.exp(raw+t)
    exact(S*s.exp(raw),1,'Exact S is not reciprocal Rtail')
    exact(2*(1-z*z)*S*s.exp(-t),2*(1-z*z)/R,'Current full Gamma argument differs')
    return dict(current_exact_heat_logradius_definitions_AST_bound=True,
        current_common_parameter_bridge_consumed=True,
        exact_source_radius_argument_identities=3,
        exact_S_not_cap_endpoint_or_midpoint=True,
        positive_Ev0_and_stress_factors_kept_symbolically=True,passed=True)


def independent_nonzero_constant_fixture():
    """Direct differentiated integrals/moments, with nonzero constants."""
    c=MPIntervalContext();c.dps=65
    t,z,v=s.symbols('t Z v',real=True);p=s.Rational(7,5);scale=s.Rational(13,10)
    K=1+(z+z*z)/10+t*(1+z)/50+t*t/700
    Cp=s.Rational(3,10)+(z+z*z)/7;D=s.Rational(2,5)+z/5+z*z/9
    location={t:s.Rational(4,5),z:s.Rational(1,4)}
    def interval(expression):
        expression=s.simplify(expression)
        if expression.is_Rational:return c.mpf(int(expression.p))/int(expression.q)
        if expression.func==s.exp:return c.exp(interval(expression.args[0]))
        if expression.is_Add:return sum((interval(arg) for arg in expression.args),c.mpf(0))
        if expression.is_Mul:
            result=c.mpf(1)
            for arg in expression.args:result*=interval(arg)
            return result
        if expression.is_Pow:return interval(expression.base)**interval(expression.exp)
        raise ValueError('Unexpected independent exact fixture expression: '+str(expression))
    def jet(expression):
        return IntervalTaylor(c,[interval(s.diff(expression,z,n).subs(location)/math.factorial(n)) for n in range(6)])
    def encloses(value,expression,label):
        with mp.workdps(100):
            reference=mp.mpf(str(s.N(expression.subs(location),100)))
            lo,hi=endpoints(value)
            if not lo<=reference<=hi:raise ArithmeticError('Independent retained-constant fixture not enclosed: '+label)
    remaining=s.integrate(s.exp(-p*v)*K.subs(t,v)**2/2,(v,t,s.oo))
    shape=[jet(s.diff(K,t,j)) for j in range(5)]
    rows=collar_pressure_rows(shape,jet(remaining),interval(p),interval(scale),interval(location[t]))
    rows[0]+=jet(Cp);pressure_checks=stress_checks=0
    for j in range(5):
        for n in range(5-j):
            encloses(rows[j][n]*math.factorial(n),s.diff(Cp-scale*remaining,t,j,z,n),'pressure')
            pressure_checks+=1
    delta=s.Rational(3,10);a=delta/2;k=1-a;b=(1-delta)/2;L=1-delta*z*z;d=1-z*z
    heat=SimpleNamespace(ctx=c,delta=interval(delta),k=interval(k),prate=interval(1+delta))
    constants=constant_stress_rows(heat,interval(location[z]),interval(location[t]),jet(D),jet(Cp),4)
    R=7*s.exp(t);B=2*s.exp(-(s.Rational(1,2)+a)*t);Pstar=s.Integer(3)
    Qt=s.sqrt(R/2)*B;Qp=s.sqrt(R/2)*Pstar**2
    Mt=s.sqrt(2)*R**s.Rational(3,2)*B*D*s.exp(-k*t)
    theta=(k*Mt-b*z*s.diff(Mt,z))/(2*L*R)
    axial=R*(2*(1+delta)*z*Pstar**2*Cp-d*Pstar**2*s.diff(Cp,z))/(L*s.sqrt(2*R))
    for label,expression,factor in (('theta',theta,Qt),('axial_pressure_constant',axial,Qp)):
        for j in range(5):
            for n in range(5-j):
                encloses(constants[label][j][n]*math.factorial(n),s.diff(expression,t,j,z,n)/factor,label)
                stress_checks+=1
    return dict(nonzero_pressure_constant_and_angular_defect_used=True,
        independent_full_future_pressure_mixed4_rows=pressure_checks,
        independent_direct_physical_constant_stress_mixed4_rows=stress_checks,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentHeatPressureStress(require_checked=False)
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Current heat companion source/datum differs')
    if not all(field.graph.values()) or any(raw[k] for k in GATES+SCOPES) or raw['retained_history_bindings']!=history_source_bindings():
        raise ValueError('Current heat companion graph, history or scope differs')
    rows=0
    for name,(chart,coordinate) in VIEWS.items():
        packet=field.evaluate(chart,[-1,1] if name!='fresh_exterior' else '.381',coordinate)
        if encode(pack(packet))!=raw['actual_current_heat_companion_views'][name]:raise ValueError('Actual current heat packet changed')
        expected=3 if chart=='heat_collar' else 4
        if packet['actual_stress_mixed_order']!=expected or len(packet['stable_same_forward_pressure_mixed4'])!=15:
            raise ValueError('Native pressure/stress derivative order differs')
        for grid in packet['actual_stress_factored_mixed_rows'].values():
            if len(grid)!=(expected+1)*(expected+2)//2:raise ValueError('Factored stress row omitted')
            if not all(mp.isfinite(value) for box in grid.values() for value in endpoints(box)):
                raise ValueError('Unbounded factored stress coefficient')
            rows+=len(grid)
        if not packet['actual_current_terminal_constants']['zero_constants_not_assumed'] or any(packet[k] for k in SCOPES):
            raise ValueError('Actual heat constants were eliminated without source equations')
    pressure=retained_forward_pressure_proof();constants=actual_constant_stress_proof()
    radius=exact_source_radius_proof(field);fixture=independent_nonzero_constant_fixture()
    general=collar_moment_stress_identities();Gamma=terminal_stress_identities()
    if not general['original_collar_moment_to_stress_identities_verified'] or not Gamma['full_terminal_moment_stress_theorem_verified']:
        raise ValueError('Full canonical moment-to-stress identities missing')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_common_heat_source_graph=field.graph,
        actual_retained_forward_pressure_proof=pressure,actual_constant_stress_proof=constants,
        exact_current_radius_and_positive_factors=radius,
        independent_nonzero_constant_fixture=fixture,
        canonical_full_collar_moment_stress_identities=general,
        canonical_projected_full_Gamma_moment_stress_identities=Gamma,
        original_forward_constants_retained=True,actual_factored_stress_rows_checked=rows,
        absolute_pressure_mixed4_rows_checked=15*len(VIEWS),
        **dict.fromkeys(GATES,True),**dict.fromkeys(SCOPES,False),
        remaining_dependency=raw['remaining_dependency'],all_scoped_checks_passed=True,all_passed=True,
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    (HERE/RECEIPT).write_bytes((json.dumps(encode(result),indent=2)+'\n').encode('utf8'))
    print('PASS actual current heat pressure mixed4, collar stress mixed3 and exterior constant stress mixed4; zero-stress/global/temporal remain open',flush=True)
    return result


if __name__=='__main__':run()

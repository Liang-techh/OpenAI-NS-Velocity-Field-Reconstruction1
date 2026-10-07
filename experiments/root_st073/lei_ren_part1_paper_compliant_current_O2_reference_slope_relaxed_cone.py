"""Full original Rh/reference and slope relaxed cone before shear repair.

At bs=0 and kappa<=2 the paper requires I_theta/F>2, not simply
T_theta>0. The actual transported angular deficit supplies this bound.
"""
import ast
import functools
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_patch_physical_numeric_bounds as numeric
import lei_ren_part1_paper_compliant_current_O3_actual_parameter_majorants as parameters
from lei_ren_part1_paper_compliant_current_O2_modified_taper_cone import OPEN
from lei_ren_part1_paper_compliant_current_original_cone_operator import original_cone_theorem
from lei_ren_part1_paper_compliant_current_O3_power_cone_operator import original_correlated_O2_O3_theorem
from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 import slope_masses

HERE,PREFIX,sha=numeric.HERE,numeric.PREFIX,numeric.sha
NAME=PREFIX+'current_O2_reference_slope_relaxed_cone.json'
RECEIPT=PREFIX+'current_O2_reference_slope_relaxed_cone_check.json'
TENSOR=PREFIX+'current_O2_background_tensor.json.gz'
TENSOR_CHECK=PREFIX+'current_O2_background_tensor_check.json'
DOMAINS={'Rh_reference':(-5,0),'O2_slope':(0,1)}
endpoints=numeric.transport.endpoints
read=numeric.transport.read_interval


def relaxed_cone_margins(c,a,bs,theta_over_F,axial_over_F,vs_minus2=None):
    """Paper (3.23), signed normalized units; crossing boxes sufficient only."""
    a,bs,theta,axial=(c.mpf(v) for v in (a,bs,theta_over_F,axial_over_F))
    if not all(mp.isfinite(v) for value in (a,bs,theta,axial) for v in endpoints(value)):
        raise ValueError('Finite signed normalized source cone inputs required')
    if endpoints(a)[0]<=0:return dict(admitted=False,status='positive_source_shear_unresolved')
    km=c.mpf(vs_minus2) if vs_minus2 is not None else a-2+bs**2/a
    if not all(mp.isfinite(v) for v in endpoints(km)):raise ValueError('Finite exact source kappa-2 required')
    D=theta-bs*axial/a;J=axial+bs*theta/a;lo,hi=endpoints(km)
    if hi<=0:
        branch='kappa_le2';margins=dict(signed_direction=D,stronger_signed_direction=D+km)
    elif lo>0:
        branch='kappa_gt2';margins=dict(signed_direction=D,signed_quadratic=2*D**2-km*J**2)
    else:
        branch='crossing_box_sufficient';negative=c.mpf([min(lo,mp.mpf(0)),0]);positive=c.mpf([0,max(hi,mp.mpf(0))])
        margins=dict(signed_direction=D,stronger_signed_direction=D+negative,
            signed_quadratic=2*D**2-positive*J**2)
    passed=all(endpoints(value)[0]>0 for value in margins.values())
    return dict(admitted=passed,status='relaxed_input_cone' if passed else 'directed_box_unresolved_or_failed',
        branch=branch,source_vs_minus2=km,margins=margins,
        strict_source_shear_for_entire_box=lo>0,global_target_admission=False)


@functools.lru_cache(maxsize=1)
def exact_theorem():
    asts=numeric.transport.SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(s.expand_power_exp(a-b)))!=0:raise ArithmeticError('Original slope relaxed identity: '+name)
        checks[name]=True
    z,y,delta=s.symbols('Z actual_logR_offset delta',real=True)
    C=1/(1+z*z);U,X,M=s.symbols('actual_scalar_U actual_h_over_U actual_scalar_M',real=True)
    L=1-delta*z*z;A=(1-delta/2)*C+(1-delta)*z*z*C*C
    B=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
    J,H,E,P=s.symbols('actual_J actual_angular_mass actual_energy_mass actual_pressure_mass',real=True)
    ctx=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp)
    env=dict(c=ctx,y=y,qi=C,z=z,self=SimpleNamespace(invP2=s.Symbol('actual_invPstar2',positive=True)),
        J=J,mass=(H,P,E),square=lambda v:v*v)
    rh_u=asts.evaluate(asts.expression('pre_pulse_mixed_C4','reference','u',wanted='qi*c.exp(y/10)'),env)
    rh_h=asts.evaluate(asts.expression('pre_pulse_mixed_C4','reference','h',wanted="u*c.mpf('.625')"),{**env,'u':rh_u})
    zero('actual_reference_h_over_U_is5_over8',rh_h/rh_u,s.Rational(5,8))
    factor=asts.evaluate(asts.expression('pre_pulse_mixed_C4','slope','factor',wanted="c.exp(y/10-c.mpf('.6')*J)"),env)
    h=asts.evaluate(asts.expression('pre_pulse_mixed_C4','slope','h',wanted="qi*(c.mpf('.625')+mass[0])*c.exp(-c.mpf('1.5')*y)"),env)
    AA=s.Rational(8,5)*y-s.Rational(3,5)*J
    zero('actual_slope_X_defining_integral',h/(C*factor),s.exp(-AA)*(s.Rational(5,8)+H))
    asts.expression('pre_pulse_mixed_C4','slope','(J, mass)',wanted='slope_masses(c,y,self.cells)')
    asts.method('pre_pulse_mixed_C4','slope_masses')
    asts.expression('pre_pulse_mixed_C4','slope','V',wanted='z*4')
    asts.expression('pre_pulse_mixed_C4','reference','V',wanted='z*4')
    asts.expression('pre_pulse_mixed_C4','slope','hist',wanted="dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(qi)*(c.mpf(5)/12+mass[2]/2)*c.exp(-y),p=square(qi)*(c.mpf('2.5')+mass[1]/2))")
    integral_name='lei_ren_part1_paper_interval_outer_slope_field.py'
    tree=ast.parse((HERE/integral_name).read_text(encoding='utf8'))
    integral=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='transition_integrals')
    assignments={ast.unparse(t):n.value for n in ast.walk(integral) if isinstance(n,ast.Assign) for t in n.targets}
    rates=asts.evaluate(assignments['rates'],dict(c=ctx))
    if rates!=(s.Rational(8,5),s.Rational(1,5),s.Rational(6,5)):raise ValueError('Actual slope integral rates differ')
    if asts.evaluate(assignments['powers'],{})!=(1,2,2):raise ValueError('Actual slope mass powers differ')
    expected=ast.parse("dy*c.exp(rate*scell-c.mpf('.6')*power*jcell)",mode='eval').body
    if not any(isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='masses[k]' and ast.dump(n.value)==ast.dump(expected) for n in ast.walk(integral)):
        raise ValueError('Original defining slope mass enclosure density changed')
    asts.hashes[integral_name]=sha(integral_name)
    checks['same_angular_mass_integrand_exp_A_bound_to_original_source']=True
    jf,hf=s.Function('actual_J')(y),s.Function('actual_H')(y);sig=s.Function('actual_sigma')(y)
    af=s.Rational(8,5)*y-s.Rational(3,5)*jf;xf=s.exp(-af)*(s.Rational(5,8)+hf)
    rules={s.diff(jf,y):sig,s.diff(hf,y):s.exp(af)}
    zero('actual_slope_X_ODE',s.diff(xf,y).subs(rules),1-(s.Rational(8,5)-s.Rational(3,5)*sig)*xf)
    zero('actual_slope_positive_deficit_integrating_factor_density',
        s.diff(s.exp(af)*(1-xf),y).subs(rules),s.Rational(3,5)*(1-sig)*s.exp(af))
    zero('actual_slope_positive_X_minus5_over8_density',
        s.diff(s.exp(af)*(xf-s.Rational(5,8)),y).subs(rules),s.Rational(3,8)*sig*s.exp(af))
    checks['same_source_deficit_at0_is3_over8_and_lower_ge3_over8_exp_minus8_over5']=True
    # Native leading theta is replayed from the complete arbitrary source
    # program; the axial stress, energy and absolute pressure are retained.
    from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
    from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
    rawop=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',dict(
        axial_derivative=lambda v:s.diff(v,z),shifted_rows=shifted_rows,product_rows=product_rows))
    lam=s.Symbol('actual_log_U_slope',real=True);pressure=s.Function('actual_absolute_pressure')(z)
    energy=s.Function('actual_full_energy')(z)
    u=[U*C,lam*U*C]+[s.Integer(0)]*3;V=[4*z]+[s.Integer(0)]*4
    history=dict(m=[4*z]+[s.Integer(0)]*4,h=[U*X*C]+[s.Integer(0)]*4,
        k=[4*z*U*X*C]+[s.Integer(0)]*4,e=[energy]+[s.Integer(0)]*4)
    raw=rawop(ctx,delta,z,u,V,history,[pressure]+[s.Integer(0)]*4)
    theta=(A*X-C)/L+4*C+4*B*X
    zero('actual_complete_inertial_theta_native0',sum(p['shape'][0] for name,p in raw['theta'].items() if name!='variable_radial_shear')/U,theta)
    D=s.Symbol('same_actual_angular_deficit',real=True)
    zero('actual_positive_deficit_full_theta_split',theta.subs(X,1-D),(A-C)/L+4*(C+B)+(-A/L-4*B)*D)
    x=s.Symbol('Z_squared',real=True);cx=(1-delta)/4
    initial=(A-C)/L
    zero('actual_reference_shape_lower_remainder',initial-cx*z*z+delta/(2*(1-delta)),
        (1-delta)*z*z*(C*C/L-s.Rational(1,4))+delta*(1/(1-delta)-C/L)/2)
    zero('actual_nonnegative_geometry_C_over_L_gap',1/(1-delta)-C/L,
        (delta+z*z*(1-2*delta)+delta*z*z*(1-z*z))/((1-delta)*(1+z*z)*L))
    F,R,Ps=s.symbols('same_F R Pstar',positive=True);a=1-2*lam
    inertial=Ps*s.sqrt(R/2)*U*theta;shear=-a*Ps*U*C/s.sqrt(2*R)
    actualF=Ps*U*C/s.sqrt(2*R)
    zero('actual_full_theta_over_shear_scale',(inertial+shear)/actualF,R*theta/C-a)
    zero('actual_relaxed_stronger_direction_gap',(inertial+shear)/actualF-(2-a),R*theta/C-2)
    zero('actual_reference_shear_a',a.subs(lam,s.Rational(1,10)),s.Rational(4,5))
    zero('actual_slope_shear_a',a.subs(lam,s.Rational(1,10)-s.Rational(3,5)*sig),s.Rational(4,5)+s.Rational(6,5)*sig)
    aa=s.Symbol('positive_source_a',positive=True);bb,tn,zn=s.symbols('signed_bs theta_over_F axial_over_F',real=True)
    kk=aa+bb*bb/aa;dot=F**2*(-aa*tn+bb*zn);cross=-F**2*(bb*tn+aa*zn)
    zero('generic_paper_relaxed_strength_in_normalized_signed_units',(-dot+F*(-F*aa)*(2-kk))/(F**2*aa),tn-bb*zn/aa+kk-2)
    zero('generic_paper_quadratic_in_normalized_signed_units',(2*dot**2-(kk-2)*cross**2)/(F**4*aa**2),
        2*(tn-bb*zn/aa)**2-(kk-2)*(zn+bb*tn/aa)**2)
    geometry=original_correlated_O2_O3_theorem();cone=original_cone_theorem()
    for key in ('actual_nonnegative_M_coefficient','actual_deficit_shape_numerator_lower_positive_polynomial',
        'actual_deficit_shape_denominator_below_four'):
        if not geometry['identities'].get(key):raise ValueError('Shared nonnegative whole-Z geometry missing')
        checks['consumed_shared_geometry/'+key]=True
    return dict(passed=True,identities=checks,source_domains=DOMAINS,
        positive_deficit_formula='D=exp(-A)*(3/8+3/5*integral(1-sigma)*exp(A)); A=8y/5-3J/5',
        source_uniform_deficit_lower='3/8*exp(-8/5)',
        native_leading_order_only=True,higher_source_rows_not_claimed_from_native0_symbolic_fixture=True,
        full_axial_stress_pressure_energy_are_retained_in_checked_source_packets=True,
        relaxed_condition_is_D_over_F_gt2_minus_kappa_not_only_theta_positive=True,
        input_hashes={**asts.hashes,**geometry['input_hashes'],**cone['input_hashes'],Path(__file__).name:sha(Path(__file__).name)})


class CurrentOriginalReferenceSlopeRelaxedCone:
    def __init__(self,require_checked=True):
        self.data=numeric.inputs();self.ctx=self.data['ctx'];self.delta=self.data['delta'];self.theorem=exact_theorem()
        self.family=self.data['accepted']['source_family'];self.hashes={**self.data['hashes'],**self.theorem['input_hashes']}
        receipt=json.loads((HERE/TENSOR_CHECK).read_bytes())
        if not receipt['all_passed'] or not receipt['current_actual_four_O2_completed_tensor_joins_certified']:
            raise ValueError('Checked original O2 full tensor/source joins required')
        for key,value in self.family.items():
            if receipt.get(key)!=value:raise ValueError('Foreign O2 tensor source family')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed original O2 tensor input: '+name)
        if receipt['input_hashes'].get(TENSOR)!=sha(TENSOR):raise ValueError('Original O2 packet not receipt-bound')
        self.tensor=json.loads(gzip.decompress((HERE/TENSOR).read_bytes()))
        for chart,domain in DOMAINS.items():
            view=self.tensor['current_actual_O2_tensor_views'][chart+'_whole'];pre=view['actual_upstream_original_pre_O2_source']
            if read(self.ctx,view['coverage_coordinate'])._mpi_!=self.ctx.mpf(domain)._mpi_ or read(self.ctx,pre['Z'])._mpi_!=self.ctx.mpf([-1,1])._mpi_:
                raise ValueError('Whole original reference/slope source packet required')
            if not view['actual_full_stress_not_local_difference'] or not pre['actual_five_histories_and_analytic_pressure_retained']:
                raise ValueError('Full original stress/history/pressure source required')
        self.hashes.update({**receipt['input_hashes'],TENSOR_CHECK:sha(TENSOR_CHECK)})
        self.proof=self.prove()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed']:raise ValueError('Checked whole original reference/slope relaxed cone required')
            for name,digest in record['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed original relaxed direction source: '+name)

    def prove(self):
        c=self.ctx;delta=self.delta;positive={}
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Reference/slope relaxed inequality unresolved: '+name)
            positive[name]=value
        require('delta_above_or_equal_zero_with_positive_gap_to1_over8',c.mpf(1)/8-delta)
        if endpoints(delta)[0]<0:raise ValueError('Nonnegative original delta required')
        Dmin=c.mpf(3)/8*c.exp(-c.mpf(8)/5);cD=(c.mpf(15)/8-14*delta)/4
        cX=(1-delta)/4;dcap=delta/(2*(1-delta));floor=cD*Dmin-dcap
        require('actual_angular_deficit_lower',Dmin);require('positive_whole_Z_deficit_coefficient',cD)
        require('full_canonical_inertial_theta_floor',floor)
        logratio=self.data['logRref']-5+c.ln(floor)
        require('actual_radius_full_theta_over_F_log_gap',logratio-1000)
        gap=c.exp(1000)-2;require('stronger_relaxed_signed_direction_reserve',gap)
        return dict(whole_source_domains=DOMAINS,whole_Z_domain=(-1,1),positive_margins=positive,
            source_angular_deficit_lower=Dmin,source_cD=cD,source_cX=cX,canonical_inertial_theta_lower=floor,
            exact_Itheta_over_F='R*Theta/C',Itheta_over_F_log_lower=logratio,
            stronger_relaxed_D_over_F_minus2_plus_kappa_lower=gap,
            original_Ttheta_over_F_lower=gap,actual_bs=0,actual_kappa_range=(c.mpf(4)/5,c.mpf(2)),
            full_original_pressure_energy_and_all_five_histories_retained=True,
            actual_source_radius_and_selected_family_used=True,
            current_original_Rh_O2_slope_relaxed_input_cone_certified=True,
            current_original_Rh_O2_slope_strict_cone_certified=False,**{key:False for key in OPEN})

    def query(self,chart,coordinate):
        if chart not in DOMAINS:raise ValueError('Original Rh_reference or O2_slope chart required')
        c=self.ctx;y=c.mpf(coordinate);lo,hi=endpoints(y);domain=DOMAINS[chart]
        if lo<domain[0] or hi>domain[1] or not all(mp.isfinite(v) for v in (lo,hi)):
            raise ValueError('Finite original reference/slope coordinate subset of declared domain required')
        if chart=='Rh_reference':U=c.exp(y/10);X=c.mpf(5)/8
        else:
            J,mass=slope_masses(c,y,256);U=c.exp(y/10-c.mpf(3)/5*J)
            X=(c.mpf(5)/8+mass[0])*c.exp(-c.mpf(3)/2*y)/U
            xl,xh=endpoints(X);yl,yh=endpoints(1-self.proof['source_angular_deficit_lower'])
            X=c.mpf([max(xl,endpoints(c.mpf(5)/8)[0]),min(xh,yh)])
        return dict(chart=chart,coordinate=y,native_logR=self.data['logRref']+y,source_family=self.family,
            actual_original_scalar_U_enclosure=U,actual_original_h_over_U_enclosure=X,
            same_source_correlated_X_bounds_applied=chart=='O2_slope',
            whole_signed_relaxed_input_certificate=self.proof,
            full_original_source_tensor_packet=TENSOR,current_original_Rh_O2_slope_relaxed_input_cone_certified=True,
            complete_modified_signed_cone_certified_for_entire_query_box=False,**{key:False for key in OPEN})


def run():
    field=CurrentOriginalReferenceSlopeRelaxedCone(require_checked=False)
    result=dict(source_family=field.family,exact_reference_slope_source_relaxed_cone_theorem=field.theorem,
        whole_original_reference_slope_relaxed_cone=field.proof,
        checked_full_source_tensor_and_four_seam_theorem=field.tensor['current_actual_O2_tensor_and_four_joins_theorem'],
        examples={name:field.query(chart,y) for name,chart,y in (
            ('whole_reference','Rh_reference',(-5,0)),('Rh','Rh_reference',-5),('reference_slope_left','Rh_reference',0),
            ('whole_slope','O2_slope',(0,1)),('reference_slope_right','O2_slope',0),('slope_interior','O2_slope','.537'),('slope_axial_seam','O2_slope',1))},
        current_original_Rh_O2_slope_relaxed_input_cone_certified=True,
        current_original_Rh_O2_slope_strict_cone_certified=False,**{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Original Rh/reference and O2 slope stronger relaxed input cone generated',flush=True)
    return result


if __name__=='__main__':run()

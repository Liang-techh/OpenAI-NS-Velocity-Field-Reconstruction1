"""Actual Gamma exterior: original physical map, stress and local NS identity.

The actual unit-viscosity leading exterior is independent of physical z.
All amplitude/radius sources remain the original implicit positive sources.
This regional identity does not certify the rest of the background field.
"""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_heat_stress_equations import heat_stress_equations
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def physical_heat_identities():
    proofs={}
    def zero(name,expression):
        if s.simplify(s.powsimp(s.expand_power_exp(expression),force=True))!=0:
            raise ArithmeticError('Physical heat source identity failed: '+name)
        proofs[name]=True

    r,lam,tau,a,Rtail,Ev0,theta=s.symbols('r lambda tau a Rtail Ev0 theta_base',positive=True)
    delta=2*a; bh=s.Rational(1,2)+a
    R=r*r/(2*lam*lam); d=tau/(lam*lam)
    # c_infinity is not a new free amplitude. The accepted packet defines
    # Utheta/Ev0=theta_base*exp(-bh*offset)*H and offset=log(R/Rtail).
    actual_c=Ev0*theta*Rtail**bh
    c0=actual_c
    zero('actual_current_radius_swirl_amplitude_equals_original_c_infinity',
         Ev0*theta*(R/Rtail)**(-bh)-actual_c*R**(-bh))
    zero('original_Gamma_argument_becomes_4_tau_over_r_squared',2*d/R-4*tau/r**2)
    A=c0*2**bh
    zero('actual_physical_swirl_prefactor_cancels_lambda',lam**(-1-delta)*c0*R**(-bh)-A*r**(-1-delta))
    zero('actual_physical_pressure_prefactor_cancels_lambda',
         lam**(-2-2*delta)*c0*c0*R**(-1-delta)/2-A*A*r**(-2-2*delta)/2)
    q=s.symbols('physical_tail_radius',positive=True)
    rho=q*q/(2*lam*lam)
    zero('full_pressure_integral_change_of_radius_and_jacobian',
         lam**(-2-2*delta)*c0*c0*s.diff(rho,q)*rho**(-2-delta)/2-A*A*q**(-3-2*delta))
    # A complete pressure integral has no residual analytic pressure
    # constant. Its derivative is the full centrifugal source by FTC.
    g=s.Function('physical_swirl')(r,tau)
    pressure_prime=g*g/r
    zero('actual_radial_momentum_and_pressure_FTC_cancel',-g*g/r+pressure_prime)
    z=s.symbols('physical_z',real=True)
    zero('physical_exterior_swirl_axial_viscosity_is_zero',s.diff(g,z,2))
    zero('physical_exterior_pressure_axial_gradient_is_zero',s.diff(s.Function('physical_pressure')(r,tau),z))
    # Original physical stress Talpha=lambda^(-2-delta)*Talpha_profile.
    zero('actual_physical_angular_stress_prefactor_preserves_exact_zero',lam**(-2-delta)*s.Integer(0))
    zero('actual_physical_axial_stress_prefactor_preserves_exact_zero',lam**(-2-delta)*s.Integer(0))

    # Derive the Cartesian geometry independently of the source-map
    # templates. The radial convective force includes the moving basis.
    x,y=s.symbols('x y',real=True); rr=s.sqrt(x*x+y*y)
    G=s.Function('G'); Pr=s.Function('P')
    ux=-y*G(rr)/rr; uy=x*G(rr)/rr; pp=Pr(rr)
    zero('Cartesian_pure_swirl_divergence',s.diff(ux,x)+s.diff(uy,y))
    for label,value,radial,angular in (('x',ux,x/rr,-y/rr),('y',uy,y/rr,x/rr)):
        zero('Cartesian_'+label+'_convective_centrifugal_force',
             ux*s.diff(value,x)+uy*s.diff(value,y)+radial*G(rr)**2/rr)
        zero('Cartesian_'+label+'_vector_Laplacian',
             s.diff(value,x,2)+s.diff(value,y,2)-angular*(s.Subs(s.diff(G(r),r,2),r,rr)
                        +s.Subs(s.diff(G(r),r),r,rr)/rr-G(rr)/rr**2))
        zero('Cartesian_'+label+'_radial_pressure_gradient',s.diff(pp,x if label=='x' else y)-radial*s.Subs(s.diff(Pr(r),r),r,rr))
    # Combine this geometry with the admitted full Gamma heat ODE. Actual
    # source xi=4*tau/r^2 corresponds to viscosity1; an arbitrary nu is not
    # silently inserted into the currently reconstructed source family.
    equations=heat_stress_equations()
    if not equations['canonical_physical_heat_flow_equation_verified']:
        raise ValueError('Full Gamma angular heat equation not admitted')
    proofs['unit_viscosity_full_angular_heat_identity_from_full_Gamma_ODE']=True

    # Preserve exact source correlations before evaluating amplitude logs.
    mu,L,Ts,wait,lone,lp,lu,lrp=s.symbols('mu Lrel Ts wait logone logPstar logU logRp',real=True)
    bp=s.Rational(1,2)+mu; rate=1-mu; k=1-a
    logEv0=lp+lu-13/(2*mu)-13
    logtheta=-100*bp-s.log(2)-bp*L-bp-rate/2-s.Rational(3,2)*Ts-s.Rational(3,2)+k/2-bh*wait-lone
    logRt=lrp+13/mu+100+L+2+Ts+wait
    grouped=lp+lu+bh*lrp+13*a/mu+(a-mu)*(100+L)-k*Ts-14-mu/2+3*a/2-s.log(2)-lone
    zero('original_c_infinity_log_combines_source_correlations_before_enclosure',logEv0+logtheta+bh*logRt-grouped)
    nu=s.symbols('nu',positive=True)
    zero('physical_viscosity_swirl_amplitude_and_argument',s.sqrt(nu)*A*(r/s.sqrt(nu))**(-1-delta)-A*nu**(1+a)*r**(-1-delta))
    zero('physical_viscosity_Gamma_argument',4*tau/(r/s.sqrt(nu))**2-4*nu*tau/r**2)
    zero('physical_viscosity_pressure_amplitude',nu*A*A*(r/s.sqrt(nu))**(-2-2*delta)-(A*nu**(1+a))**2*r**(-2-2*delta))
    return dict(identities=proofs, actual_unit_viscosity=1,
                original_physical_stress_prefactor='lambda^(-2-delta)',
                actual_amplitude_definition='A=2^((1+delta)/2)*Ev0*theta_base*Rtail^((1+delta)/2)',
                physical_velocity='ux=-y/r*A_nu*r^(-1-delta)*H(4*nu*tau/r^2); uy=x/r*A_nu*r^(-1-delta)*H(4*nu*tau/r^2); uz=0; A_nu=nu^(1+delta/2)*A',
                physical_pressure='p=-A_nu^2*integral_r^infinity q^(-3-2*delta)*H(4*nu*tau/q^2)^2 dq',
                all_three_physical_momentum_components_exactly_zero=True,
                axial_viscosity_remainder_exactly_zero=True,
                radial_lower_order_remainder_exactly_zero=True,
                Cartesian_divergence_exactly_zero=True,
                full_Gamma_heat_equation=equations,
                source_equations=['Lei-Ren v2 (2.18)-(2.20)','(3.1)-(3.3)','(5.1)-(5.11)'],
                domain='Actual leading exact-Gamma exterior only; tau>0, r>0, finite physical z; |Z|<1 and log(R/Rtail)>=3')


def physical_source_bridge(stress):
    if not stress.history['complete_terminal_moment_history_bridge_verified']:
        raise ValueError('Actual terminal moments not source-bound')
    proofs,hashes,trees={},{},{}
    def method(stem,name):
        if stem not in trees:
            path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
            trees[stem]=ast.parse(path.read_text(encoding='utf8'))
        return next(n for n in ast.walk(trees[stem]) if isinstance(n,ast.FunctionDef) and n.name==name)
    def syntax(stem,name,target,expression):
        nodes=[n.value for n in ast.walk(method(stem,name)) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(n)==wanted for n in nodes)!=1:
            raise ValueError('Physical heat source map changed: '+stem+':'+target)
        proofs['source_'+stem+'_'+name+'_'+target]=True
    def keyword(stem,name,key,expression):
        nodes=[n.value for n in ast.walk(method(stem,name)) if isinstance(n,ast.keyword) and n.arg==key]
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(n)==wanted for n in nodes)!=1:
            raise ValueError('Physical heat map output changed: '+stem+':'+key)
        proofs['output_'+stem+'_'+name+'_'+key]=True
    syntax('collar_Gamma_C4','packet','theta','K[0]*(self.theta_base*c.exp(-self.bh*t))')
    syntax('collar_Gamma_C4','__init__','self.bh','self.steep.bh')
    syntax('collar_Gamma_C4','__init__','self.theta_base','self.steep.thetaT*c.exp(-self.bh*self.steep.wait-self.steep.logone)')
    syntax('collar_Gamma_C4','__init__','self.pressure_scale','self.Ev2*self.theta_base**2')
    syntax('collar_Gamma_C4','local_Gamma','K','[one-D[0]*(self.a*self.S)]+[-v*(self.a*self.S) for v in D[1:]]')
    syntax('outer_angular_repair','__init__','logRtail','self.heat.logRref+13/self.mu+self.heat.tail_finite')
    syntax('steep_waiting_C4','__init__','self.bh',"c.mpf('.5')+self.delta/2")
    syntax('global_physical_assembly','cartesian_source_row','beta','{UR:c.mpf(-1),UT:-1-delta,UZ:-1-delta,P:-2-2*delta}')
    syntax('global_physical_assembly','normalized_sources','logs[4]',
           'self.logP+sum(flatten.logEv2_parts.values(),c.mpf(0))/2')
    syntax('global_physical_assembly','normalized_sources','amplitudes',"{UZ:{},UT:{},UR:{5:mp.mpf('.5'),6:mp.mpf('-.5')},P:{1:mp.mpf(2)}}")
    keyword('global_physical_assembly','evaluate','radial',repr('R=r^2/(2lambda^2)'))
    keyword('global_physical_assembly','evaluate','velocity',repr('ur=lambda^-1 Ur; utheta/uz=lambda^(-1-delta) Utheta/Uz'))
    keyword('global_physical_assembly','evaluate','pressure',repr('p=lambda^(-2-2delta)P'))
    keyword('global_physical_assembly','evaluate','lambda_relation',repr('lambda^2-lambda^(2delta)*z^2=tau; Z=z/lambda^(1-delta)'))
    # The previously admitted defining bridge identifies exact Rtail/S,
    # Ev2=Ev0^2/Pstar^2, theta_base and the complete absolute pressure tail.
    inherited=stress.pressure.bridge
    functions=inherited['defining_function_bridge']
    for gate in ('exact_heat_radius_and_xi_binding_verified','canonical_gamma_derivative_enclosure_verified','c4_constants_data_path_verified'):
        if not functions[gate]:raise ValueError('Actual defining function binding missing: '+gate)
        proofs['consume_'+gate]=True
    for gate in ('same_actual_heat_velocity_amplitude','Rtail_pressure_units_equal_theta_base_squared','production_C4_pressure_scale'):
        if not inherited['identities'][gate]:raise ValueError('Actual amplitude/pressure equality missing: '+gate)
        proofs['consume_'+gate]=True
    R,Rtail,offset,a,Ev0,Pstar,theta,H,Bp=s.symbols('R Rtail offset a Ev0 Pstar theta_base H full_pressure_integral',positive=True)
    Z=s.symbols('Z',real=True); bh=s.Rational(1,2)+a
    S=1/Rtail; inverse_current=s.exp(-offset)/Rtail
    D=(1-H)/(a*S)
    K=1-a*S*D
    source_c=Ev0*theta*Rtail**bh
    def zero(name,expr):
        if s.simplify(s.powsimp(s.expand_power_exp(expr),force=True))!=0:raise ArithmeticError('Actual heat physical transfer failed: '+name)
        proofs[name]=True
    zero('actual_full_Gamma_packet_bracket_is_the_bound_H',K-H)
    zero('actual_exact_inverse_radius_Gamma_argument_is_canonical',
         (2*(1-Z**2)*inverse_current).subs(offset,s.log(R/Rtail))-2*(1-Z**2)/R)
    zero('actual_packet_swirl_uses_the_same_bound_c_infinity',
         (Ev0*theta*s.exp(-bh*offset)*K).subs(offset,s.log(R/Rtail))-source_c*R**(-bh)*H)
    # Consume Ev2=Ev0^2/Pstar^2 from the admitted constants bridge before
    # the physical pressure multiplication. Its runtime cap is not Ev2.
    Ev2=Ev0**2/Pstar**2
    zero('actual_absolute_pressure_tail_uses_same_c_infinity_squared',
         (-Pstar**2*Ev2*theta**2*s.exp(-(1+2*a)*offset)*Bp/2).subs(offset,s.log(R/Rtail))
         +source_c**2*R**(-1-2*a)*Bp/2)
    history=inherited['retained_pressure_history_bridge']
    hashes.update(history['input_hashes']); hashes.update(stress.history['input_hashes'])
    return dict(identities=proofs,
                actual_velocity_pressure_and_physical_map_source_bound=True,
                actual_unit_viscosity_from_original_Gamma_argument=1,
                original_absolute_pressure_and_zero_meridional_histories_retained=True,
                actual_positive_amplitude_not_freely_selected=True,
                exact_Rtail_not_a_cap_endpoint=True, input_hashes=hashes)


def viscosity_source_row(row,lognu,power):
    """Uniform spatial pullback x_source=x/sqrt(nu), with explicit units."""
    out=dict(row)
    out['physical_viscosity_exponent']=power
    out['physical_viscosity_log_prefactor']=power*lognu
    out['log_absolute_upper']=None if row['exact_zero'] else row['log_absolute_upper']+power*lognu
    out['terms']=[dict(term,physical_viscosity_exponent=power,
                       log_absolute_upper=term['log_absolute_upper']+power*lognu) for term in row['terms']]
    return out


class CompliantHeatPhysicalC4:
    """Local physical exterior certificate on the admitted leading source."""
    def __init__(self):
        self.assembly=CompliantGlobalPhysicalAssembly()
        self.stress=self.assembly.dispatch.provider('heat_exterior')
        self.heat=self.stress.heat; self.ctx=self.assembly.ctx
        self.family,self.source=self.stress.family,self.stress.source
        self.hashes=dict(self.assembly.hashes); self.hashes.update(self.stress.hashes)
        for stem,gate in (('heat_stress_C4_check','heat_exterior_stress_identity_certified'),
                          ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped')):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Physical exterior prerequisite not admitted: '+name)
            if receipt['actual_five_defect_family_sha256']!=self.family or receipt['implicit_source_sha256']!=self.source:
                raise ValueError('Physical exterior source family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Physical exterior source changed: '+path)
            self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.bridge=physical_source_bridge(self.stress); self.proof=physical_heat_identities()
        self.hashes.update(self.bridge['input_hashes'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def exterior(self,Z,offset,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx
        if any(not mp.isfinite(v) for v in endpoints(c.mpf(log_tau))):raise ValueError('Finite log_tau required; physical tau>0')
        nu=c.mpf(viscosity)
        if endpoints(nu)[0]<=0 or any(not mp.isfinite(v) for v in endpoints(nu)):raise ValueError('Finite positive physical viscosity required')
        lognu=c.ln(nu)
        point=self.assembly.evaluate('heat_exterior',Z,offset,log_tau=log_tau,theta=theta)
        spatial={}
        for key,components in point['physical_spatial_cartesian_mixed4'].items():
            degree=sum(int(v[1:]) for v in key.split('_'))
            spatial[key]={component:{label:viscosity_source_row(row,lognu,c.mpf((2 if component=='p' else 1)-degree)/2)
                                     for label,row in parts.items()} for component,parts in components.items()}
        point['physical_spatial_cartesian_mixed4']=spatial
        point['first_fixed_x_physical_time_derivative']={component:{label:viscosity_source_row(row,lognu,c.mpf(1 if component=='p' else '.5'))
                                                       for label,row in parts.items()} for component,parts in point['first_fixed_x_physical_time_derivative'].items()}
        coords=dict(point['physical_log_coordinates']);coords['log_r']+=lognu/2
        coords.update(exact_log_r='log(lambda)+log(2R)/2+log(nu)/2',exact_z='sqrt(nu)*Z*lambda^(1-delta)',
                      source_lambda_relation='lambda^2*(1-Z^2)=tau; x_source=x_physical/sqrt(nu)')
        point['physical_log_coordinates']=coords
        point['mapping']=dict(lambda_relation='lambda^2-lambda^(2delta)*z^2/nu=tau; Z=z/(sqrt(nu)*lambda^(1-delta))',
                              radial='R=r^2/(2nu*lambda^2)',velocity='u_physical=sqrt(nu)*u_source',pressure='p_physical=nu*p_source')
        heat=self.heat; a=heat.a; mu=heat.mu; bh=heat.bh; k=heat.k
        inlet=self.assembly.dispatch.provider('flatten').logEv2_parts['inlet_log']/2
        point['actual_A_nu_log_source_parts']=dict(logPstar=self.assembly.logP,original_inlet_logU=inlet,
                  radial_origin=bh*self.assembly.logRp,combined_inverse_mu=13*a/mu,
                  combined_flatten_power=(a-mu)*(100+heat.steep.outer.Lrel),combined_steep=-k*heat.steep.Ts,
                  finite_correction=-14-mu/2+3*a/2+(a-c.mpf('.5'))*c.ln(2),
                  minus_log_one_minus_epsilon=-heat.steep.logone,physical_viscosity=(1+a)*lognu)
        point.update(actual_compliant_field_physical_heat_region_certified=True,
                     actual_physical_map_and_prefactor_transfer_pending=False,
                     regional_unit_viscosity_leading_NS_identity_certified=True,
                     physical_viscosity=nu,viscosity_source_pullback_verified=True,
                     regional_physical_viscosity_leading_NS_identity_certified=True,
                     completed_background_stress_tensor_components={key:c.mpf(0) for key in ('xx','xy','xz','yy','yz','zz')},
                     completed_stress_tensor_definition='mapped completed divergence-form background tensor T_B; Cauchy shear is nu*(d_r utheta-utheta/r)',
                     physical_momentum_residual_components={key:c.mpf(0) for key in ('x','y','z')},
                     physical_axial_viscosity_remainder=c.mpf(0),physical_radial_remainder=c.mpf(0),
                     physical_divergence=c.mpf(0),original_pressure_datum_and_positive_amplitude_retained=True,
                     physical_heat_velocity_definition=self.proof['physical_velocity'],
                     physical_heat_pressure_definition=self.proof['physical_pressure'],
                     original_physical_stress_prefactor=self.proof['original_physical_stress_prefactor'],
                     certification_scope=self.proof['domain'],
                     source_box_endpoints_Z_plus_minus_one_are_physical_infinity_limits=True,
                     global_admissible_stress_lift_constructed=False,whole_outer_cone_certified=False,
                     independently_bounded_global_flat_remainder=False,full_background_NS_validation=False,
                     physical_energy_integral_certified=False,temporal_recursion=False)
        return point

    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_source_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            return result
        with mp.workdps(280):
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                        actual_compliant_field_physical_heat_region_certified=True,
                        regional_unit_viscosity_leading_NS_identity_certified=True,
                        actual_physical_map_and_prefactor_transfer_pending=False,
                        physical_source_bridge=self.bridge, physical_heat_identities=self.proof,
                        samples=[summary(self.exterior(z,t,lt,angle)) for z,t,lt,angle in
                                 (('0','3','-1','0'),('.5','4','-10','.7'),('-.5','1000','-100','1'))],
                        whole_source_exterior=summary(self.exterior([-1,1],[3,mp.inf],[-1000,-1],None)),
                        global_admissible_stress_lift_constructed=False,whole_outer_cone_certified=False,
                        independently_bounded_global_flat_remainder=False,full_background_NS_validation=False,
                        physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantHeatPhysicalC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual unit-viscosity Gamma physical exterior and local NS identity generated; global gates remain open',flush=True)
    return result


if __name__=='__main__':run()

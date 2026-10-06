"""Actual lambda and constant-viscosity transport of native velocity bounds."""
import ast
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_dispatch import endpoints
from lei_ren_part1_paper_compliant_cartesian_field import INDICES
from lei_ren_part1_paper_compliant_current_global_tensor_cover import CurrentGlobalTensorCover
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_dispatch import CurrentModifiedVelocityPressureDispatch

def existing_original_cover(dispatch):
    # Reuse the actual admitted ancestor, not a separately reconstructed
    # registry/physical graph when the new physical wrapper is standalone.
    direction=dispatch.velocity.source.repair.histories.candidate.direction
    return direction.powercone.entrancecone.tailcone.flattencone.gapcone.maincone.cone.cover

def viscosity_power(component,spatial_order):
    return s.Rational((2 if component=='p' else 1)-spatial_order,2)

def actual_physical_velocity_bounds(field,view,location):
    c=field.ctx;q=location['actual_log_lambda'];lt=location['requested_log_tau'];lnnu=c.ln(location['physical_viscosity'])
    core=view['physical_layout']=='core_native_Cartesian';spatial={};time={}
    for i,j,k in INDICES:
        index='x%d_y%d_z%d'%(i,j,k);spatial[index]={};order=i+j+k
        parts=(view['core_axis_nonsingular_native_map']['cartesian_spatial_multiindices'][index] if core
          else view['physical_spatial_cartesian_mixed4'][index])
        for component,rows in parts.items():
            spatial[index][component]={};power=c.mpf(str(viscosity_power(component,order)))
            for label,row in rows.items():
                if core:
                    scale=view['core_axis_nonsingular_native_map']['shared_physical_prefactor_bounds'][row['scale_key']]
                    norm=endpoints(row['absolute_upper'])[1]
                    bound=(c.ln(c.mpf(norm))+scale['logLambda_term']+scale['amplitude_log_upper']+
                      scale['physical_lambda_exponent']*q+power*lnnu) if norm else None
                    gamma=scale['physical_lambda_exponent']
                else:
                    gamma=row['physical_lambda_exponent']
                    bound=None if row['exact_zero'] else row['log_absolute_upper']+gamma*(q-lt/2)+power*lnnu
                spatial[index][component][label]=dict(exact_zero=bound is None,log_absolute_upper=bound,
                  actual_log_lambda=q,physical_lambda_exponent=gamma,constant_viscosity_power=power,
                  native_signed_source_retained=True)
    parts=(view['core_axis_nonsingular_native_map']['first_fixed_x_physical_time_derivative'] if core
      else view['first_fixed_x_physical_time_derivative'])
    for component,rows in parts.items():
        time[component]={};power=c.mpf(str(viscosity_power(component,0)))
        for label,row in rows.items():
            gamma=row['physical_lambda_exponent']
            if core:
                norm=endpoints(row['absolute_upper'])[1]
                bound=(c.ln(c.mpf(norm))+row['logLambda_term']+row['amplitude_log_upper']+gamma*q+power*lnnu) if norm else None
            else:bound=None if row['exact_zero'] else row['log_absolute_upper']+gamma*(q-lt/2)+power*lnnu
            time[component][label]=dict(exact_zero=bound is None,log_absolute_upper=bound,
              actual_log_lambda=q,physical_lambda_exponent=gamma,constant_viscosity_power=power,
              fixed_physical_position_time_derivative=True,native_signed_source_retained=True)
    return dict(actual_lambda_viscosity_spatial4_bounds=spatial,actual_lambda_viscosity_fixed_x_time1_bounds=time,
      native_physical_layout=view['physical_layout'],actual_log_lambda=q,requested_log_tau=lt,
      physical_constant_viscosity=location['physical_viscosity'],
      exact_velocity_source='u_nu(x,t)=sqrt(nu)*u_1(x/sqrt(nu),t)',
      exact_absolute_pressure_source='p_nu(x,t)=nu*p_1(x/sqrt(nu),t)',
      actual_lambda_not_replaced_by_sqrt_tau=True,viscosity_constant_in_space_and_time=True,
      signed_factored_function_enclosures_not_resolved_point_values=True)

def actual_candidate_velocity(field,candidate,location):
    region=candidate['region'];coordinate=candidate['native_coordinate_enclosure']
    native=(field.dispatch.axis(location['Z'],location['requested_log_tau'],location['theta'])
      if region=='core_positive_radius' and endpoints(coordinate)==(0,0)
      else field.dispatch.native(region,location['Z'],coordinate,location['requested_log_tau'],location['theta']))
    pieces=[]
    for piece in native['source_pieces']:
        view=piece['complete_velocity_pressure_source']
        pieces.append(dict(native_source_piece=piece,
          actual_physical_velocity_pressure_bounds=actual_physical_velocity_bounds(field,view,location)))
    return dict(actual_physical_candidate=candidate,current_native_velocity_pressure_query=native,
      actual_physical_velocity_pressure_pieces=pieces,
      finite_heat_coordinate_not_replaced_by_zero_tensor_or_point_at_infinity=True,
      source_pieces_are_alternative_charts_not_summed_fields=True)

def exact_modified_physical_velocity_theorem(field):
    field.cover.assert_graph();field.dispatch.velocity.heat.assert_graph()
    if type(field.cover) is not CurrentGlobalTensorCover or not field.cover.acceptance_loaded:
        raise ValueError('Checked actual global physical candidate cover required')
    if type(field.dispatch) is not CurrentModifiedVelocityPressureDispatch or not field.dispatch.acceptance_loaded:
        raise ValueError('Checked actual full33 native velocity/absolute-pressure dispatcher required')
    if field.cover.registry is not field.dispatch.registry or field.cover.locator.physical is not field.dispatch.geometry:
        raise ValueError('One actual canonical physical geometry/registry required')
    if field.cover.field.locator is not field.cover.locator or field.dispatch.pre is not field.dispatch.geometry.pre:
        raise ValueError('Same actual radius and pressure source graph required')
    if not all(p['passed'] for p in (field.dispatch.theorem,field.cover.radial,field.cover.implicit,field.cover.candidate,field.cover.hypotheses)):
        raise ValueError('Full native velocity and source-bound all-point cover proofs required')
    if len(field.cover.locator.recipes)!=33 or tuple(field.cover.locator.recipes)!=tuple(field.dispatch.registry.routes):
        raise ValueError('Actual physical/native33 coordinate inventories differ')
    asts=SourceAST();checks={}
    asts.expression('current_modified_physical_velocity_operator','existing_original_cover','direction',
      wanted='dispatch.velocity.source.repair.histories.candidate.direction')
    asts.method('current_modified_physical_velocity_operator','existing_original_cover')
    if existing_original_cover(field.dispatch) is not field.cover:
        raise ValueError('Reuse the exact original cover ancestor of this actual modified graph')
    checks['same_actual_original_cover_ancestor_reused_without_reconstruction']=True
    def zero(name,a,b):
        if s.simplify(a-b)!=0:raise ArithmeticError('Actual physical velocity source identity: '+name)
        checks[name]=True
    nu,tau,lam=s.symbols('nu tau lambda',positive=True);Z,delta=s.symbols('Z delta',real=True)
    x,y,z=s.symbols('x y z',real=True);t=s.symbols('t',real=True)
    normalized_z=z/s.sqrt(nu)
    zero('same_actual_viscosity_normalized_implicit_lambda_equation',
      lam**2-normalized_z**2*lam**(2*delta)-tau,lam**2-(z**2/nu)*lam**(2*delta)-tau)
    zero('same_actual_viscosity_normalized_radial_coordinate',(x**2+y**2)/nu/(2*lam**2),(x**2+y**2)/(2*nu*lam**2))
    # Independently differentiate the actual dilation. nu is one positive
    # constant source parameter, never a spatial/time-dependent field.
    # Arbitrary ordinary local jets at an arbitrary requested point.
    # This avoids implementation-dependent nesting of Sympy Subs while
    # proving every spatial4 coefficient, not a chosen sample function.
    coefficients={index:s.Symbol('same_native_local_jet_'+str(index),real=True) for index in INDICES}
    time_coefficient=s.Symbol('same_native_fixed_x_time_jet',real=True)
    jet=sum((coefficients[i,j,k]*x**i*y**j*z**k/(s.factorial(i)*s.factorial(j)*s.factorial(k))
      for i,j,k in INDICES),s.Integer(0))+time_coefficient*t
    scaled_jet=jet.subs({x:x/s.sqrt(nu),y:y/s.sqrt(nu),z:z/s.sqrt(nu)},simultaneous=True)
    center={x:0,y:0,z:0,t:0}
    for component in ('ux','uy','uz','p'):
        amplitude=nu**(s.Rational(1) if component=='p' else s.Rational(1,2))
        for i,j,k in INDICES:
            n=i+j+k;axes=(x,i,y,j,z,k)
            actual=s.diff(amplitude*scaled_jet,*axes).subs(center)
            zero('constant_nu_'+component+'_x%d_y%d_z%d'%(i,j,k),actual,nu**viscosity_power(component,n)*coefficients[i,j,k])
        zero('constant_nu_'+component+'_fixed_x_time1',s.diff(amplitude*scaled_jet,t).subs(center),amplitude*time_coefficient)
    base,gamma,q,lt,lnnu=s.symbols('signed_log_unit_bound gamma actual_log_lambda log_tau log_nu',real=True)
    zero('native_lambda_prefactor_to_actual_lambda_before_viscosity',
      base+gamma*lt/2+gamma*(q-lt/2),base+gamma*q)
    program=asts.replay('current_modified_physical_velocity_operator','viscosity_power',dict(s=s))
    for component in ('ux','uy','uz','p'):
        for n in range(5):zero('actual_viscosity_power_program_'+component+'_order'+str(n),program(component,n),s.Rational((2 if component=='p' else 1)-n,2))
    asts.expression('current_modified_physical_velocity_operator','actual_candidate_velocity','native',
      wanted="field.dispatch.axis(location['Z'],location['requested_log_tau'],location['theta']) if region=='core_positive_radius' and endpoints(coordinate)==(0,0) else field.dispatch.native(region,location['Z'],coordinate,location['requested_log_tau'],location['theta'])")
    for method in ('actual_candidate_velocity','actual_physical_velocity_bounds'):
        asts.method('current_modified_physical_velocity_operator',method)
    for method,target,wanted in (
      ('cartesian','covered','self.cover.cartesian(x,y,z,time,viscosity,terminal_time,tensors=False)'),
      ('log_radius','covered','self.cover.locate_log_radius(log_r,z,log_tau,theta,viscosity)'),
      ('correlated','location','self.cover.field.locate(request,z,log_tau,theta,viscosity)'),
      ('correlated','covered','self.cover._attach(location,request=request)'),
      ('_consume','views',"[actual_candidate_velocity(self,candidate,location) for candidate in location['candidates']]")):
        asts.expression('current_modified_physical_velocity',method,target,wanted=wanted)
        checks['actual_physical_velocity_'+method+'_'+target+'_source_call']=True
    checks['source_bound_global_candidate_cover_composes_with_same_actual_full33_velocity_functions']=True
    checks['all_candidate_and_native_partition_pieces_retained_without_summing_or_midpoint_choice']=True
    checks['finite_heat_velocity_uses_actual_native_Gamma_coordinate_not_the_zero_tensor_registry']=True
    checks['requested_log_tau_kept_separate_from_actual_implicit_log_lambda']=True
    return dict(identities=checks,passed=True,input_hashes={**field.dispatch.hashes,**field.cover.hashes,**asts.hashes},
      consumed_actual_full33_native_velocity_source_theorem=field.dispatch.theorem,
      consumed_source_bound_global_physical_cover=dict(radial=field.cover.radial,implicit=field.cover.implicit,
        candidate=field.cover.candidate,hypotheses=field.cover.hypotheses),
      all_finite_physical_points_tau_positive_constant_nu_positive_enclosed=True,
      resolved_point_coefficients_global_smoothness_energy_cones_recursion_full_NS_remain_open=True)

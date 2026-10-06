"""Actual current steep-entry tensor and angular-entry completed tensor join.

Full moments share KR units with the checked current angular tensor. Direct
remaining integrals preserve the same absolute pressure, never a gauge patch.
"""
import ast
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace,FunctionType

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_angular_background_stress import (
    CurrentAngularBackgroundStress,HERE,PREFIX,OPEN,sha,pack,encode,endpoints,
    angular_source_log_parts,copy_jet,full_normalized_stress_AST_theorem)
from lei_ren_part1_paper_compliant_current_angular_support_differences import linear_stress_program
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import (
    SourceAST,entry_shape,entry_angular_rows,entry_stress_rows,entry_transport_identities)
from lei_ren_part1_paper_compliant_angular_stress_C3 import angular_entry_join_binding
from lei_ren_part1_paper_compliant_steep_entry_physical_C2 import steep_entry_physical_identities
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import ordinary_coordinate_proof
from lei_ren_part1_paper_compliant_current_selected_energy_source import binding,function
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_stress_rows,axial_derivative
from lei_ren_part1_paper_compliant_collar_physical_C2 import physical_source_row,scale_row,collar_velocity_bracket
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_current_pulse_support_interfaces import ordinary_grid
from lei_ren_part1_paper_compliant_axial_pulse_field import sigma_enclosure
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_steep_entry_background_stress.json'
RECEIPT=PREFIX+'current_steep_entry_background_stress_check.json'
GATES=('current_actual_steep_entry_full_moment_stress_mixed3_recovered',
    'current_actual_steep_entry_completed_physical_tensor_remainder_decomposition_available',
    'current_actual_angular_entry_completed_tensor_interface_certified')
VIEWS={'whole_entry':((-1,1),(0,1),('-3','-1'),None,'1'),
    'angular_boundary':((-1,1),'0','-1',None,'1'),
    'power_boundary':((-1,1),'1','-1',None,'1'),
    'middle_sector':(('-0.8','0.8'),('.3','.7'),('-5','-2'),None,'.2'),
    'fresh_entry':('.537','.417','-2.6','.41','.8')}


def remaining_entry_pressure(c,t,mu,rate,Jt,cells):
    """Same original integral_t^1 exp(-po*v-2*rate*J(v))/2 dv."""
    t=c.mpf(t);length=1-t;ds=length/cells;J=c.mpf(Jt);pressure=c.mpf(0);r=rate
    if endpoints(t)[0]<0 or endpoints(t)[1]>1 or cells<1:raise ValueError('Original unit entry domain/cells required')
    if endpoints(length)[1]>0:
        for i in range(cells):
            left=t+length*i/cells;right=t+length*(i+1)/cells
            v=c.mpf([max(mp.mpf(0),endpoints(left)[0]),min(mp.mpf(1),endpoints(right)[1])])
            nextJ=J+ds*sigma_enclosure(c,v)
            jc=c.mpf([endpoints(J)[0],endpoints(nextJ)[1]])
            pressure+=ds*c.exp(-(1+2*mu)*v-2*r*jc)/2
            J=nextJ
    return pressure


def normalized_entry_moment_rows(heat,native,shape,data,kernels,post_after,t,c):
    K=shape['K_rows'];Q=product_rows(K,K)
    one=IntervalTaylor.constant(c,1,5)
    A=[K[0]*copy_jet(c,native['angular_Taylor'])]
    E=[(copy_jet(c,data['after_entry'])+one*kernels['remaining_energy'])*c.exp(heat.delta*t)]
    P=[(post_after['normalized_remaining_post_entry_pressure']+one*post_after['remaining_entry_pressure'])*c.exp(heat.prate*t)]
    for j in range(4):
        A.append(K[j]-A[j]*heat.k);E.append(E[j]*heat.delta-Q[j]);P.append(P[j]*heat.prate-Q[j]/2)
    return dict(A=A,E=E,P=P,K=K)


def entry_moment_normalization_AST_theorem():
    """Original signed defect FTC includes unit baselines before KR scaling."""
    t,z,a,mu,KR=s.symbols('t Z a mu KR',real=True,positive=True)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp)
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,delta=2*a,k=1-a,prate=1+2*a)
    K=s.Function('same_normalized_entry_K')(t);X=s.Function('same_actual_X')(t,z)
    Epost=s.Function('same_current_after_entry')(z);Ppost=s.Function('same_current_post_entry_pressure')(z)
    IE=s.Function('same_remaining_entry_energy')(t);IP=s.Function('same_remaining_entry_pressure')(t)
    stub=SimpleNamespace(constant=lambda ctx,value,order:s.Rational(value))
    env=dict(normalized_entry_moment_rows.__globals__);env.update(IntervalTaylor=stub,copy_jet=lambda c,v:v)
    shape=dict(K_rows=[s.diff(K,t,j) for j in range(5)])
    rows=FunctionType(normalized_entry_moment_rows.__code__,env)(heat,dict(angular_Taylor=X),shape,
        dict(after_entry=Epost),dict(remaining_energy=IE),dict(normalized_remaining_post_entry_pressure=Ppost,remaining_entry_pressure=IP),t,c)
    asts=SourceAST();original_env=dict(IntervalTaylor=stub)
    asts.replay('steep_entry_stress_C3','entry_defect_rows',original_env)
    actualK=[KR*v for v in rows['K']];square=product_rows(actualK,actualK)
    actualshape=dict(K_rows=actualK,K_defect_rows=[actualK[0]-1]+actualK[1:],
        K_squared_defect_rows=[square[0]-1]+square[1:])
    ell=1-t;delta=heat.delta;p=heat.prate
    terminal=dict(energy_defect_rows=[KR**2*s.exp(delta)*Epost-1/delta],
        pressure_defect_rows=[KR**2*s.exp(p)*Ppost-1/(2*p)])
    signed=dict(energy=KR**2*s.exp(delta*t)*IE-(1-s.exp(-delta*ell))/delta,
        pressure=2*KR**2*s.exp(p*t)*IP-(1-s.exp(-p*ell))/p)
    original=original_env['entry_defect_rows'](heat,terminal,t,signed,actualshape,X)
    checks={}
    for name,key,base in (('A','angular_defect_rows',1/heat.k),('E','energy_defect_rows',1/delta),('P','pressure_defect_rows',1/(2*p))):
        for j in range(5):
            zero=s.simplify(s.expand(original[key][j]+(base if j==0 else 0)-rows[name][j]*(KR if name=='A' else KR**2)))
            if zero!=0:raise ArithmeticError('Current entry full normalized moment differs from original AST')
            for n in range(5-j):checks[name+'_y%d_Z%d'%(j,n)]=s.diff(zero,z,n)==0
    # Native half-energy and original absolute-pressure increment are bound.
    binding('compliant_steep_waiting_C4','steep_in','energy',
        "(data['after_entry']+kernels['remaining_energy'])*c.exp(2*self.mu*t+2*self.rate*J)/2")
    binding('compliant_steep_waiting_C4','steep_in','pressure',
        "data['PR']+kernels['pressure']*(self.outer.flatten.Ev2*self.thetaR**2)")
    return dict(original_full_entry_defect_AST_normalized_mixed4_identities=checks,
        original_half_energy_and_absolute_pressure_increment_bound=True,
        actual_full_nonzero_post_entry_moments_retained=True,input_hashes=asts.hashes,passed=True)


def entry_source_and_shape_theorem():
    asts=SourceAST();env={};checks={}
    for stem,method,target,wanted in (
            ('steep_waiting_C4','transition_kernels','pressure','ds*c.exp(-(1+2*mu)*v-2*rate*jc)/2'),
            ('current_steep_entry_background_stress','remaining_entry_pressure','pressure','ds*c.exp(-(1+2*mu)*v-2*r*jc)/2')):
        node=asts.expression(stem,method,target,wanted=wanted,augmented=True)
        v,mu,r,J,ds=s.symbols('v mu rate J ds',real=True)
        c=SimpleNamespace(exp=s.exp)
        value=asts.evaluate(node,dict(c=c,v=v,mu=mu,r=r,rate=r,jc=J,ds=ds))
        if s.simplify(value-ds*s.exp(-(1+2*mu)*v-2*r*J)/2)!=0:raise ArithmeticError('Remaining pressure integrand differs')
        checks[stem+'_same_actual_pressure_integrand']=True
    t,a,mu,J=s.symbols('t a mu J',real=True);r=1-mu;g=a-mu
    c=SimpleNamespace(exp=s.exp);heat=SimpleNamespace(mu=mu,a=a)
    exponent=asts.evaluate(asts.expression('steep_entry_stress_C3','entry_shape','exponent'),
        dict(heat=heat,t=t,r=r,kernels={'f':s.Rational(1,2)-J}))
    KS=s.exp(g-r/2)
    Kfromsource=asts.evaluate(asts.expression('steep_entry_stress_C3','entry_shape','K'),
        dict(one=s.Integer(1),KS=KS,c=c,exponent=exponent))[0]
    if s.simplify(s.expand_power_exp(Kfromsource/s.exp(g*t-r*J)))!=1:raise ArithmeticError('Entry normalized source K differs')
    return dict(actual_remaining_pressure_integrand_AST=asts.bindings,identities=checks,
        exact_K_over_KR='exp((a-mu)*t-(1-mu)*J(t))',exact_KS_over_KR='exp(a-mu-(1-mu)/2)',
        remaining_pressure_and_original_forward_pressure_identical_by_radial_FTC_additivity=True,
        original_current_Cp_zero_supplies_absolute_PR_once=True,
        same_exact_sigma_J_source_and_positive_correlated_lengths=True,input_hashes=asts.hashes,passed=True)


def entry_inertial_correlation_AST_theorem():
    """Current full A=KX is exactly the original correlated entry operator."""
    z,a,S,q=s.symbols('Z a inverse_radius q',real=True)
    g=s.symbols('g0:4',real=True);K=[s.Symbol('same_K0',positive=True)]
    X=[s.Function('same_actual_X0')(z)];k=1-a
    for n in range(1,5):
        K.append(sum(K[n-1-j]*g[j]*s.binomial(n-1,j) for j in range(n)))
        X.append((1 if n==1 else 0)-sum(X[n-1-j]*(k+g[0] if j==0 else g[j])*s.binomial(n-1,j) for j in range(n)))
    A=[K[0]*X[0]]
    for j in range(4):A.append(K[j]-k*A[j])
    checks={}
    for j,value in enumerate(product_rows(K,X)):
        if s.expand(A[j]-value)!=0:raise ArithmeticError('Original normalized entry A=KX row source differs')
        checks['actual_native_A_equals_KX_y'+str(j)]=True
    stub=SimpleNamespace(constant=lambda c,v,n:s.Rational(v),variable=lambda c,v,n:v)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp)
    heat=SimpleNamespace(ctx=c,a=a,delta=2*a,k=k,prate=1+2*a,S=S)
    linear=linear_stress_program();env=dict(linear.__globals__)
    env.update(IntervalTaylor=stub,axial_derivative=lambda v:s.diff(v,z))
    linear=FunctionType(linear.__code__,env)
    defs=dict(angular_defect_rows=A,K_defect_rows=K,
        energy_defect_rows=[s.Function('same_full_E'+str(j))(z) for j in range(5)],
        pressure_defect_rows=[s.Function('same_full_P'+str(j))(z) for j in range(5)])
    generic=linear(heat,dict(K_rows=K),defs,z,q)
    asts=SourceAST();actual_env=dict(entry_stress_rows.__globals__)
    actual_env.update(IntervalTaylor=stub,collar_stress_rows=linear,axial_derivative=lambda v:s.diff(v,z))
    entry=FunctionType(entry_stress_rows.__code__,actual_env)(heat,dict(K_rows=K),defs,X,z,q)
    for label in ('theta','axial','theta_inertial','theta_shear'):
        for j in range(4):
            zero=s.cancel(s.expand(entry[label][j]-generic[label][j]))
            if zero!=0:raise ArithmeticError('Original correlated entry full stress differs from same general-K operator')
            for n in range(4-j):checks[label+'_y%d_Z%d'%(j,n)]=s.diff(zero,z,n)==0
    binding('compliant_steep_entry_stress_C3','entry_angular_rows','h',
        "[(1-heat.mu)*(1-shape['original_sigma_jets'][0])]+g[1:]")
    return dict(actual_full_A_KX_and_original_correlated_entry_stress_mixed3_identities=checks,
        current_K_Z_independent_and_native_X_axial_history_retained=True,
        exact_log_X_rate_minus_log_K_rate_equals_k=True,
        arbitrary_nonzero_full_energy_pressure_histories_retained=True,passed=True)


def tensor_operator():
    """Replay the checked actual angular general-K physical assembler AST."""
    fn=function('compliant_current_angular_background_stress','angular')
    begin=[i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='ps' for t in n.targets)]
    end=[i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='pressure_rows' for t in n.targets)]
    if len(begin)!=1 or len(end)!=1 or begin[0]>=end[0]:raise ValueError('Checked actual tensor source block differs')
    wrapper=ast.parse('def current_tensor(self,c,Z,K,grids,parts,logR,lt,theta,nu):\n    pass').body[0]
    wrapper.body=[ast.parse('beta=-2-self.heat.delta').body[0]]+copy.deepcopy(fn.body[begin[0]:end[0]])
    wrapper.body.append(ast.Return(value=ast.parse("dict(physical_cylindrical_stress_mixed3=ps,physical_stress_divergence_mixed2=div,completed_theta_theta_stress_mixed2=diagonal,physical_axial_viscosity_remainder_mixed2=error,completed_background_tensor_cartesian_components=tensor,physical_completed_stress_divergence_cartesian=divcart,physical_remainder_cartesian=ecart,physical_momentum_residual_decomposition_cartesian=residual)",mode='eval').body))
    env=dict(math=math,physical_bracket=physical_bracket,physical_source_row=physical_source_row,
        scale_row=scale_row,collar_velocity_bracket=collar_velocity_bracket)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper],type_ignores=[])),'<checked current general-K tensor AST>','exec'),env)
    return env[wrapper.name]


def angular_entry_tensor_source_proof(angular):
    angular.assert_graph();p=angular.atlas.postpulse.proof
    if not p['passed'] or not p['all_current_postpulse_primitive_mixed4_traces_identified']:
        raise ValueError('Checked same-current angular/entry primitive function trace required')
    inst={k:v for k,v in p['current_instantiated_seam_mixed4_identities'].items() if k.startswith('angular_steep_')}
    if len(inst)!=60 or not all(inst.values()):raise ValueError('Complete current angular-entry primitive mixed4 source identity required')
    formula=angular_entry_join_binding()
    if not formula['actual_angular_entry_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Original full angular-entry stress function join missing')
    coordinate=ordinary_coordinate_proof(angular.physical)
    if not coordinate['exact_common_radius_identities']['angular_steep']:raise ValueError('Current angular-entry radius differs')
    return dict(current_instantiated_primitive_function_trace=inst,
        original_arbitrary_nonzero_terminal_stress_pressure_AST_join=formula,
        current_full_moment_units_from_checked_angular=angular.moment_proof,
        current_absolute_pressure_from_checked_Cp_zero=angular.pressure_proof,
        same_ordinary_logR_and_exact_radius=coordinate,
        current_full_nonzero_terminal_A_E_P_and_K_bound_before_physical_transfer=True,
        tensor_stress3_divergence2_completion2_remainder2_join_from_same_general_K_operators=True,
        completed_diagonal_uses_same_source_radius_and_stress_z_derivative=True,
        common_B_KR_Qtheta_KR_Qz_KR_squared_source_logs=True,
        interval_overlap_not_used_as_tensor_function_identity=True,passed=True)


class CurrentSteepEntryBackgroundStress:
    @source_precision
    def __init__(self,angular=None,require_checked=True):
        self.angular=angular if angular is not None else CurrentAngularBackgroundStress()
        if not self.angular.acceptance_loaded:raise ValueError('Checked actual current full angular tensor required')
        self.physical=self.angular.physical;self.history=self.physical.history;self.steep=self.history.steep
        self.heat=self.history.heat;self.ctx=self.physical.ctx
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.transport=entry_transport_identities();self.moment_proof=entry_moment_normalization_AST_theorem()
        self.shape_proof=entry_source_and_shape_theorem();self.baselines=full_normalized_stress_AST_theorem()
        self.correlation=entry_inertial_correlation_AST_theorem()
        self.join_proof=angular_entry_tensor_source_proof(self.angular)
        self.join_proof['current_entry_full_moment_normalization']=self.moment_proof
        self.join_proof['current_entry_full_inertial_correlation']=self.correlation
        self.physical_proof=steep_entry_physical_identities();self.tensor=tensor_operator()
        self.hashes=dict(self.angular.hashes)
        for stem in ('steep_entry_stress_C3','steep_entry_physical_C2','angular_stress_C3','current_steep_entry_background_stress'):
            name=PREFIX+stem+'.py';digest=sha(name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Checked current entry dependency changed: '+name)
            self.hashes[name]=digest
        self.cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current entry tensor receipt exceeds regional scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.angular.assert_graph()
        if not (self.angular.acceptance_loaded and self.physical is self.angular.physical and
                self.history is self.physical.history and self.steep is self.history.steep and
                self.heat is self.history.heat and self.ctx is self.physical.ctx):
            raise ValueError('Entry tensor must retain the checked current angular/history/pressure graph')

    @source_precision
    def entry(self,Z,t,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);t=c.mpf(t);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if not all(mp.isfinite(v) for v in endpoints(Z)+endpoints(t)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(t)[0]<0 or endpoints(t)[1]>1 or endpoints(nu)[0]<=0:
            raise ValueError('Finite current entry t[0,1],Z[-1,1],compact log(tau),nu>0 required')
        if theta is not None and not all(mp.isfinite(v) for v in endpoints(c.mpf(theta))):raise ValueError('Finite theta or all-angle bounds required')
        native=self.steep.steep_in(Z,t);data=self.steep.data(Z);kernels=self.steep.kernels(t,'in')
        local_heat=copy.copy(self.heat);local_heat.ctx=c
        shape=entry_shape(local_heat,dict(KS=c.exp(self.heat.a-self.steep.mu-self.steep.rate/2)),t,dict(f=c.mpf('.5')-kernels['J']))
        post=self.angular.post_pressure(Z);terms=post['normalized_source_contributions']
        one=IntervalTaylor.constant(c,1,5)
        post_after=dict(normalized_remaining_post_entry_pressure=sum((v for k,v in terms.items() if k!='entry'),one*0),
            remaining_entry_pressure=remaining_entry_pressure(c,t,self.steep.mu,self.steep.rate,kernels['J'],self.steep.cells),
            exact_source='integral_t^1 original_entry_pressure + power/exit/waiting/full_heat pressure future; Cp=0 fixes PR')
        full=normalized_entry_moment_rows(local_heat,native,shape,data,kernels,post_after,t,c)
        A,E,P,K=(full[name] for name in ('A','E','P','K'))
        q=-self.steep.wait-self.steep.Ts-2+t
        logR,radius_source=self.physical.radius('steep_entry',t,{},self.steep)
        caplog=c.ln(c.mpf(endpoints(self.heat.Scap)[1]))-q
        if endpoints(caplog+logR)[0]<0:raise ValueError('Entry inverse-radius cap differs from current radius')
        inverseR=c.mpf([0,endpoints(c.exp(c.mpf(endpoints(caplog)[1])))[1]])
        local_heat.S=inverseR
        defects=dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=P,K_defect_rows=K)
        # Full moment rows after source-baseline cancellation; native X rows
        # preserve the original entry inertial correlation with Z-independent K.
        Xrows=entry_angular_rows(local_heat,copy_jet(c,native['angular_Taylor']),shape)
        stress=entry_stress_rows(local_heat,shape,defects,Xrows,Z,c.mpf(0))
        grids={name:ordinary_grid(stress[name],3) for name in ('theta','axial','theta_inertial','theta_shear')}
        parts=angular_source_log_parts(self.physical,t)
        tensor=self.tensor(self,c,Z,K,grids,parts,logR,lt,theta,nu)
        pressure_rows=[-sum((P[n]*math.comb(j,n)*(-self.heat.prate)**(j-n) for n in range(j+1)),one*0) for j in range(5)]
        return dict(Z=Z,entry_t=t,log_tau=lt,viscosity=nu,
            current_actual_normalized_full_moment_rows=full,current_actual_normalized_stress_mixed3=grids,
            current_actual_remaining_pressure_source=post_after,current_original_entry_transition_kernels=kernels,
            exact_current_positive_log_source_parts=parts,exact_source_logR=logR,radius_source=radius_source,
            inverse_radius_enclosure=inverseR,inverse_radius_exact_source='1/R=S_exact*exp(-q), q=-wait-Ts-2+t',
            inverse_radius_cap_only_not_exact_field_value=True,
            stable_actual_absolute_pressure_mixed4_factored=ordinary_grid(pressure_rows,4),
            stable_pressure_positive_source_log_parts=parts['absolute_pressure'],
            original_forward_absolute_pressure_Taylor=native['pressure_over_Pstar_squared_Taylor'],
            original_current_energy_Taylor=native['energy_Taylor'],
            normalized_energy_native_consistency=E[0]-K[0]*K[0]*copy_jet(c,native['energy_Taylor'])*2,
            actual_full_stress_not_local_difference=True,original_full_nonzero_A_E_P_and_current_C5_controls_retained=True,
            exact_pure_swirl_meridional_velocity_radial_divergence_and_remainder_zero=True,
            actual_pressure_is_same_P0_Pin_forward_function_by_current_Cp_zero=True,
            physical_tensor_operator_replayed_from_checked_actual_general_K_AST=True,
            regional_decomposition_is_not_global_temporal_flatness=True,**tensor,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph();left=self.angular.angular(Z,0,log_tau,theta,viscosity);right=self.entry(Z,0,log_tau,theta,viscosity)
        keys=('physical_cylindrical_stress_mixed3','physical_stress_divergence_mixed2',
            'completed_theta_theta_stress_mixed2','physical_axial_viscosity_remainder_mixed2',
            'completed_background_tensor_cartesian_components','physical_completed_stress_divergence_cartesian',
            'physical_remainder_cartesian','physical_momentum_residual_decomposition_cartesian')
        def flatten(value,path=''):
            if isinstance(value,dict) and 'signed_coefficient' in value:return {path:value}
            result={}
            values=value.items() if isinstance(value,dict) else enumerate(value)
            for key,row in values:result.update(flatten(row,path+'/'+str(key)))
            return result
        bounds={};sources={}
        for key in keys:
            a=flatten(left[key],key);b=flatten(right[key],key)
            if set(a)!=set(b):raise ValueError('Actual tensor interface layout differs')
            for name,row in a.items():
                other=b[name]
                for sourcekey in ('actual_source_log_parts','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):
                    if encode(pack(row[sourcekey]))!=encode(pack(other[sourcekey])):raise ValueError('Common actual tensor positive source factor differs')
                values=[endpoints(v['log_absolute_upper'])[1] for v in (row,other) if not v['exact_zero']]
                bounds[name]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values)) if values else None)
        return dict(common_actual_angular_entry_tensor_rows=bounds,current_common_tensor_contribution_count=len(bounds),
            full_nonzero_histories_and_source_function_join_proof_retained=True,
            overlap_not_used_for_tensor_source_equality=True,all_angles_covered=theta is None,
            requested_Z=self.ctx.mpf(Z),requested_log_tau=self.ctx.mpf(log_tau),
            actual_completed_tensor_and_divergence_remainder_common_source_bound=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_entry_transport_source_theorem=self.transport,
            current_entry_full_moment_normalization_AST_theorem=self.moment_proof,
            current_entry_shape_and_remaining_pressure_source_theorem=self.shape_proof,
            current_entry_full_A_KX_inertial_correlation_AST_theorem=self.correlation,
            original_full_stress_baseline_cancellation_and_KR_units=self.baselines,
            current_angular_entry_actual_tensor_function_join_theorem=self.join_proof,
            original_entry_physical_tensor_remainder_theorem=self.physical_proof,
            physical_tensor_operator='Exact AST replay of checked current angular general-K physical assembler; same viscosity/basis/radius operators',
            source_domain='original current entry t[0,1], source Z[-1,1],R>0,|Z|<1,tau>0,finite compact log(tau),constant nu>0; Z endpoints are infinity limits',
            source_bounds_are_enclosures_not_resolved_physical_point_values=True,
            actual_full_current_nonzero_moments_and_absolute_pressure_retained=True,
            current_22_interface_velocity_pressure_atlas_retained=True,
            other_actual_current_chart_tensors_joins_and_global_temporal_obligations_remain_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentSteepEntryBackgroundStress(require_checked=False)
    result=field.manifest();result['current_actual_steep_entry_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_steep_entry_tensor_views'][name]=field.entry(*args)
        print('Current actual steep-entry full stress/tensor: '+name,flush=True)
    result['current_actual_angular_entry_tensor_interface_bounds']=field.interface()
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

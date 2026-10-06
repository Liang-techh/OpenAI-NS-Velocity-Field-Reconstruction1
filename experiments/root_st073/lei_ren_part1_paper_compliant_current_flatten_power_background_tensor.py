"""Actual current flatten/power tensors, preserving the complete future.

The original variable axial flatten K and forward X are retained. Quadratic
moments use direct remaining integrals in the checked common KR units.
"""
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace, FunctionType

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_heat_background_tensor import (
    CurrentHeatBackgroundTensor,HERE,PREFIX,OPEN,sha,pack,encode,endpoints,
    copy_jet,angular_source_log_parts,current_full_moment_rows,TENSOR_KEYS,
    source_precision,accepted,_verify_hashes,SourceAST,ordinary_grid)
from lei_ren_part1_paper_compliant_flatten_stress_C3 import (
    flatten_shape,flatten_remaining_kernels,flatten_transport_identities,flatten_power_join_binding)
from lei_ren_part1_paper_compliant_outer_power_stress_C3 import (
    outer_power_shape,outer_power_transport_identities,outer_power_angular_join_binding)
from lei_ren_part1_paper_compliant_current_angular_support_differences import linear_stress_program
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import ordinary_coordinate_proof
from lei_ren_part1_paper_compliant_current_selected_energy_source import binding
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

NAME=PREFIX+'current_flatten_power_background_tensor.json'
RECEIPT=PREFIX+'current_flatten_power_background_tensor_check.json'
CHARTS=('flatten','outer_power')
SEAMS=('flatten_power','power_angular')
GATES=('current_actual_flatten_outer_power_full_moment_tensor_available',
    'current_actual_flatten_outer_power_physical_decomposition_available',
    'current_actual_flatten_power_angular_two_tensor_joins_certified')
VIEWS={'flatten_whole':('flatten',(-1,1),(0,100),('-3','-1'),None,'1'),
    'flatten_inlet':('flatten',(-1,1),'0','-1',None,'1'),
    'flatten_exit':('flatten',(-1,1),'100','-1',None,'1'),
    'flatten_fresh':('flatten','.537','41.7','-2.6','.41','.8'),
    'power_whole':('outer_power',(-1,1),(0,1),('-3','-1'),None,'1'),
    'power_inlet':('outer_power',(-1,1),'0','-1',None,'1'),
    'power_exit':('outer_power',(-1,1),'1','-1',None,'1'),
    'power_fresh':('outer_power','.731','.417','-2.6','.41','.8')}


def full_power_rows(heat,K,X,Eright,Pright,offset,c):
    """Remaining quadratic history from angular s=-4, ordinary logR rows."""
    one=IntervalTaylor.constant(c,1,5);square=K[0]*K[0]
    E0=Eright*c.exp(heat.delta*offset)+square*decay_integral(c,2*heat.mu,-offset)
    P0=Pright*c.exp(heat.prate*offset)+square*(decay_integral(c,1+2*heat.mu,-offset)/2)
    return current_full_moment_rows(heat,K,K[0]*X,E0,P0)


def full_flatten_rows(heat,K,X,Eright,Pright,Kright,offset,kernels,c):
    """Use F/f correlated suffix kernels, without subtracting past energy."""
    factor=Kright*Kright*c.exp(2*(heat.a-heat.mu)*offset)
    E0=Eright*c.exp(heat.delta*offset)+kernels['energy']*factor
    P0=Pright*c.exp(heat.prate*offset)+kernels['pressure']*factor
    return current_full_moment_rows(heat,K,K[0]*X,E0,P0)


def full_normalization_AST_theorem():
    """Replay both original defect implementations on arbitrary full data."""
    a,mu,z,q,KR,Kr=s.symbols('a mu Z offset KR normalized_Kright',positive=True)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp,expm1=lambda v:s.exp(v)-1)
    h=SimpleNamespace(ctx=c,a=a,mu=mu,delta=2*a,k=1-a,prate=1+2*a)
    stub=SimpleNamespace(constant=lambda c,v,n:s.Rational(v))
    K=[s.Function('same_normalized_K'+str(j))(z) for j in range(5)]
    X=s.Function('same_native_forward_X')(z)
    Er=s.Function('same_full_right_E')(z);Pr=s.Function('same_full_right_P')(z)
    kernels={k:s.Function('same_suffix_'+k)(z) for k in ('energy','pressure')}
    checks={};asts=SourceAST()
    for stem,method,formula in (('flatten_stress_C3','flatten_defect_rows',full_flatten_rows),
            ('outer_power_stress_C3','outer_power_defect_rows',full_power_rows)):
        env=dict(formula.__globals__)
        env.update(IntervalTaylor=stub,decay_integral=lambda c,p,t:(1-s.exp(-p*t))/p)
        formula=FunctionType(formula.__code__,env)
        full=(formula(h,K,X,Er,Pr,Kr,q,kernels,c) if stem=='flatten_stress_C3'
            else formula(h,K,X,Er,Pr,q,c))
        actualK=[KR*v for v in K];Q=[]
        for j in range(5):Q.append(sum(actualK[n]*actualK[j-n]*s.binomial(j,n) for n in range(j+1)))
        shape=dict(K_rows=actualK,K_defect_rows=[actualK[0]-1]+actualK[1:],
            K_squared_defect_rows=[Q[0]-1]+Q[1:])
        terminal=dict(Kright=KR*Kr,energy_defect_rows=[KR**2*Er-1/h.delta],
            pressure_defect_rows=[KR**2*Pr-1/(2*h.prate)])
        original=asts.replay(stem,method,dict(IntervalTaylor=stub))
        result=(original(h,terminal,q,shape,X,kernels) if stem=='flatten_stress_C3'
            else original(h,terminal,q,shape,X))
        for name,key,base,unit in (('A','angular_defect_rows',1/h.k,KR),
                ('E','energy_defect_rows',1/h.delta,KR**2),('P','pressure_defect_rows',1/(2*h.prate),KR**2)):
            for j in range(5):
                residual=result[key][j]+(base if j==0 else 0)-unit*full[name][j]
                if stem=='outer_power_stress_C3':
                    # The power K has the original constant logarithmic rate.
                    residual=residual.subs(K[0],Kr*s.exp((a-mu)*q))
                residual=s.cancel(s.expand(residual))
                if residual!=0:raise ArithmeticError('Current full flatten/power AST normalization differs')
                for n in range(5-j):checks[stem+'_'+name+'_y%d_Z%d'%(j,n)]=s.diff(residual,z,n)==0
    return dict(original_full_defect_AST_KR_normalization_mixed4_identities=checks,
        arbitrary_full_nonzero_axial_functions_used=True,actual_K_Z_dependence_not_discarded=True,
        positive_remaining_integrals_before_enclosure=True,input_hashes=asts.hashes,passed=True)


def current_source_units_proof(heat_tensor):
    heat_tensor.assert_graph();angular=heat_tensor.o7.entry_source.angular
    asts=SourceAST();a,mu,L,y,t,F,z=s.symbols('a mu Lrel y t F Z',real=True)
    bp=s.Rational(1,2)+mu;bh=s.Rational(1,2)+a;d=a-mu;f=(1+z*z)/2
    c=SimpleNamespace(exp=s.exp);native=SimpleNamespace(bp=bp)
    h=SimpleNamespace(ctx=c,a=a,mu=mu)
    stub=SimpleNamespace(constant=lambda c,v,n:s.Rational(v))
    sj=[s.Function('original_sigma_t'+str(j))(t) for j in range(5)]
    env=dict(flatten_shape.__globals__)
    env.update(IntervalTaylor=stub,log_taylor=s.log,flatten_f=lambda c,z:(1+z*z)/2)
    actual_flatten=asts.replay('flatten_stress_C3','flatten_shape',env)(h,
        dict(Kright=s.exp(-d*L)),t-100,
        dict(Z=z,F_Taylor=F,sigma_y_derivatives=sj,actual_original_flatten_right_endpoint=False))
    actual_power=asts.replay('outer_power_stress_C3','outer_power_shape',
        dict(outer_power_shape.__globals__,IntervalTaylor=stub))(h,dict(Kright=s.exp(-4*d)),y-L+4)
    thetaR=s.exp(-bp*(100+L))/2
    flatten_theta=asts.evaluate(asts.expression('flatten_mixed_C4','flatten','theta',wanted='F/q*c.exp(-self.bp*t)'),
        dict(F=F,q=2*f,c=c,self=native,t=t))
    power_theta=asts.evaluate(asts.expression('power_angular_C4','power','theta'),
        dict(one=1,c=c,self=native,y=y))
    checks={}
    for name,theta,v,K in (('flatten',flatten_theta,-L+t-100,s.exp(d*(-L+t-100))*F/f),
            ('outer_power',power_theta,-L+y,s.exp(d*(-L+y)))):
        residual=s.simplify(s.expand_power_exp(theta/(thetaR*s.exp(-bh*v))/K))
        if residual!=1:raise ArithmeticError('Actual native theta/current K unit differs')
        checks[name+'_actual_theta_K_over_KR_identified']=True
        actual=(actual_flatten if name=='flatten' else actual_power)['K_rows'][0]
        if s.simplify(s.expand_power_exp(actual/K))!=1:raise ArithmeticError('Actual original shape/current theta K differs')
        checks[name+'_original_shape_K_source_identified']=True
    rate=actual_flatten['log_K_rate_rows'][0]
    if s.simplify(rate-(d+s.log(f)*sj[1]))!=0:raise ArithmeticError('Original flatten variable axial/radial log-rate differs')
    checks['flatten_variable_K_rate_and_full_axial_dependence_identified']=True
    asts.expression('current_flatten_power_background_tensor','chart','Kright',wanted='c.exp(-(h.a-h.mu)*c.mpf(self.outer.Lrel))')
    asts.expression('current_flatten_power_background_tensor','power_full','shape',wanted='outer_power_shape(h,dict(Kright=c.exp(-4*d)),offset)')
    binding('compliant_flatten_mixed_C4','flatten','Fc','(rho*sigma_jets(c,v/100)[0]).exp()')
    asts.expression('flatten_stress_C3','flatten_remaining_kernels','density',
        wanted='(rho*(2*(sig-1))).exp()')
    binding('compliant_power_angular_C4','power','X',
        'one/self.rate+(data[\'flatten_exit_X\']-1/self.rate)*c.exp(-self.rate*y)')
    checks['flatten_same_original_F_squared_density_divided_by_endpoint_f_squared']=True
    return dict(native_theta_and_suffix_density_AST_bindings=asts.bindings,identities=checks,
        checked_actual_complete_energy_transfer=angular.atlas.postpulse.proof['current_complete_energy_substitution'],
        checked_current_pressure_and_full_future=angular.pressure_proof,
        checked_current_closed_exterior=heat_tensor.units,
        same_common_KR_units=True,pressure_sign='absolute p=-positive_B_squared_KR_squared*P',
        native_forward_X_preserved_without_backward_exponential_amplification=True,
        full_quadratic_future_transferred_by_same_density_FTC_additivity=True,
        flatten_actual_variable_axial_K_and_K_Z_retained=True,
        no_tiny_amplitude_cap_division_or_radius_materialization=True,input_hashes=asts.hashes,passed=True)


def joined_source_proofs(heat_tensor,units,normalization):
    angular=heat_tensor.o7.entry_source.angular;post=angular.atlas.postpulse.proof
    coordinate=ordinary_coordinate_proof(heat_tensor.physical)
    formulas=(flatten_power_join_binding(),outer_power_angular_join_binding());result={}
    for seam,formula in zip(SEAMS,formulas):
        primitive={k:v for k,v in post['current_instantiated_seam_mixed4_identities'].items() if k.startswith(seam+'_')}
        if len(primitive)!=60 or not all(primitive.values()) or not coordinate['exact_common_radius_identities'][seam]:
            raise ValueError('Complete current flatten/power source seam missing')
        flags=[v for k,v in formula.items() if k.endswith('_verified')]
        if not flags or not all(flags):raise ValueError('Original full stress/pressure function join missing')
        result[seam]=dict(current_primitive_mixed4_function_identities=primitive,
            original_arbitrary_full_terminal_stress_pressure_AST_join=formula,
            same_common_KR_full_moment_normalization=normalization,same_complete_energy_pressure_source=units,
            same_ordinary_logR_and_source_radius=coordinate,
            same_general_K_completed_tensor_stress3_divergence2_completion2_remainder2=True,
            source_function_equality_precedes_common_bounds=True,interval_overlap_not_used_for_source_equality=True,passed=True)
    return result


class CurrentFlattenPowerBackgroundTensor:
    @source_precision
    def __init__(self,heat_tensor=None,require_checked=True):
        self.heat_tensor=heat_tensor if heat_tensor is not None else CurrentHeatBackgroundTensor()
        if not self.heat_tensor.acceptance_loaded:raise ValueError('Checked actual current heat tensor chain required')
        self.physical=self.heat_tensor.physical;self.history=self.physical.history;self.heat=self.history.heat
        self.flatten=self.history.flatten;self.outer=self.history.outer;self.angular=self.heat_tensor.o7.entry_source.angular
        self.ctx=self.physical.ctx;self.family=self.heat_tensor.family;self.source=self.heat_tensor.source;self.datum_sha=self.heat_tensor.datum_sha
        self.normalization=full_normalization_AST_theorem();self.units=current_source_units_proof(self.heat_tensor)
        self.transport=dict(flatten=flatten_transport_identities(),outer_power=outer_power_transport_identities())
        self.joins=joined_source_proofs(self.heat_tensor,self.units,self.normalization)
        self.tensor=self.heat_tensor.tensor;self.stress=linear_stress_program();self.cache={}
        self.hashes=dict(self.heat_tensor.hashes)
        for stem in ('flatten_stress_C3','outer_power_stress_C3','current_flatten_power_background_tensor'):
            name=PREFIX+stem+'.py';digest=sha(name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Current checked flatten/power dependency changed')
            self.hashes[name]=digest
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Current flatten/power tensor receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.heat_tensor.assert_graph()
        if not (self.heat_tensor.acceptance_loaded and self.physical is self.heat_tensor.physical and self.history is self.physical.history and self.heat is self.history.heat and self.flatten is self.history.flatten and self.outer is self.history.outer and self.angular is self.heat_tensor.o7.entry_source.angular and self.angular.acceptance_loaded and self.ctx is self.physical.ctx):
            raise ValueError('Current flatten/power tensor must retain same checked complete physical graph')

    def right_data(self,Z):
        key=endpoints(Z)
        if key not in self.cache:
            view=self.angular.angular(Z,-4)
            full=view['current_actual_normalized_full_moment_rows']
            self.cache[key]=dict(E=full['E'][0],P=full['P'][0])
        return self.cache[key]

    def power_full(self,Z,x,h):
        c=self.ctx;one=IntervalTaylor.constant(c,1,5);d=h.a-h.mu;length=c.mpf(self.outer.Lrel)-4
        y=length*x;offset=-length*(1-x);v=-4+offset
        native=self.outer.power(Z,x);right=self.right_data(Z)
        shape=outer_power_shape(h,dict(Kright=c.exp(-4*d)),offset)
        full=full_power_rows(h,shape['K_rows'],copy_jet(c,native['angular_Taylor']),right['E'],right['P'],offset,c)
        return full,shape,native,v,dict(remaining_power_length=-offset,complete_angular_right_E=right['E'],complete_angular_right_P=right['P'])

    @source_precision
    def chart(self,chart,Z,x,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);x=c.mpf(x);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        upper=100 if chart=='flatten' else 1
        if chart not in CHARTS or not all(mp.isfinite(v) for v in endpoints(Z)+endpoints(x)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(x)[0]<0 or endpoints(x)[1]>upper or endpoints(nu)[0]<=0:
            raise ValueError('Original flatten t[0,100]/power phase[0,1],Z[-1,1],finite logtau and nu>0 required')
        if theta is not None and not all(mp.isfinite(v) for v in endpoints(c.mpf(theta))):raise ValueError('Finite theta or all-angle bounds required')
        h=copy.copy(self.heat);h.ctx=c;h.mu=c.mpf(self.outer.mu)
        if chart=='outer_power':full,shape,native,v,remaining=self.power_full(Z,x,h)
        else:
            native=self.flatten.flatten(Z,x);right,_,_,_,_=self.power_full(Z,c.mpf(0),h)
            offset=x-100;v=-c.mpf(self.outer.Lrel)+offset;Kright=c.exp(-(h.a-h.mu)*c.mpf(self.outer.Lrel))
            packet=dict(Z=Z,F_Taylor=copy_jet(c,native['F_Taylor']),sigma_y_derivatives=[c.mpf(v) for v in native['sigma_y_derivatives']],
                actual_original_flatten_right_endpoint=endpoints(x)==(mp.mpf(100),mp.mpf(100)))
            shape=flatten_shape(h,dict(Kright=Kright),offset,packet)
            kernels=flatten_remaining_kernels(c,h.mu,Z,x,self.flatten.cells)
            full=full_flatten_rows(h,shape['K_rows'],copy_jet(c,native['angular_Taylor']),right['E'][0],right['P'][0],Kright,offset,kernels,c)
            remaining=dict(original_direct_flatten_suffix=kernels,complete_power_inlet_E=right['E'][0],complete_power_inlet_P=right['P'][0])
        K=full['K'];q=-c.mpf(self.history.steep.wait)-c.mpf(self.history.steep.Ts)-2+v
        provider=self.flatten if chart=='flatten' else self.outer
        logR,radius_source=self.physical.radius(chart,x,{},provider)
        caplog=c.ln(c.mpf(endpoints(self.heat.Scap)[1]))-q
        if endpoints(caplog+logR)[0]<0:raise ValueError('Current flatten/power inverse-radius cap disagrees with source radius')
        h.S=c.mpf([0,endpoints(c.exp(c.mpf(endpoints(caplog)[1])))[1]])
        defects=dict(angular_defect_rows=full['A'],energy_defect_rows=full['E'],pressure_defect_rows=full['P'],K_defect_rows=K)
        # Checked exact constant-baseline cancellation and KR homogeneity.
        # The general operator retains K_Z through d_Z A; no independent
        # axial-K approximation or zero-history shortcut is applied.
        stress=self.stress(h,shape,defects,Z,c.mpf(0))
        grids={name:ordinary_grid(stress[name],3) for name in ('theta','axial','theta_inertial','theta_shear')}
        parts=angular_source_log_parts(self.physical,v);tensor=self.tensor(self,c,Z,K,grids,parts,logR,lt,theta,nu)
        one=IntervalTaylor.constant(c,1,5);P=full['P']
        pressure_rows=[-sum((P[n]*math.comb(j,n)*(-h.prate)**(j-n) for n in range(j+1)),one*0) for j in range(5)]
        return dict(chart=chart,Z=Z,original_coordinate=x,ordinary_logR_offset_from_Rrel=v,log_tau=lt,viscosity=nu,
            current_actual_normalized_full_moment_rows=full,current_actual_normalized_stress_mixed3=grids,
            stable_complete_remaining_quadratic_sources=remaining,
            exact_current_positive_log_source_parts=parts,exact_source_logR=logR,radius_source=radius_source,
            stable_actual_absolute_pressure_mixed4_factored=ordinary_grid(pressure_rows,4),stable_pressure_positive_source_log_parts=parts['absolute_pressure'],
            original_forward_absolute_pressure_Taylor=(native['pressure']['P_over_Pstar_squared'] if chart=='flatten' else native['pressure_over_Pstar_squared_Taylor']),
            original_current_energy_Taylor=native['energy_Taylor'],
            normalized_energy_native_consistency=full['E'][0]-K[0]*K[0]*copy_jet(c,native['energy_Taylor'])*2,
            actual_full_stress_not_local_difference=True,actual_variable_axial_K_retained=chart=='flatten',
            current_complete_pressure_energy_and_forward_X_retained=True,ordinary_logR_derivatives_not_phase_derivatives=True,
            source_bounds_not_resolved_physical_point_values=True,regional_result_not_global_temporal_flatness=True,
            **tensor,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,seam,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph()
        if seam not in SEAMS or not self.joins[seam]['passed']:raise ValueError('Checked current flatten/power function seam required')
        left=self.chart('flatten' if seam=='flatten_power' else 'outer_power',Z,100 if seam=='flatten_power' else 1,log_tau,theta,viscosity)
        right=self.chart('outer_power',Z,0,log_tau,theta,viscosity) if seam=='flatten_power' else self.angular.angular(Z,-4,log_tau,theta,viscosity)
        def flatten(value,path):
            if isinstance(value,dict) and 'signed_coefficient' in value:return {path:value}
            result={}
            for key,row in (value.items() if isinstance(value,dict) else enumerate(value)):result.update(flatten(row,path+'/'+str(key)))
            return result
        bounds={}
        for key in TENSOR_KEYS:
            a=flatten(left[key],key);b=flatten(right[key],key)
            if set(a)!=set(b):raise ValueError('Common actual flatten/power tensor layout differs')
            for name,row in a.items():
                other=b[name]
                for sk in ('physical_lambda_exponent','physical_viscosity_exponent'):
                    if encode(pack(row[sk]))!=encode(pack(other[sk])):raise ValueError('Common tensor physical exponents differ')
                values=[endpoints(v['log_absolute_upper'])[1] for v in (row,other) if not v['exact_zero']]
                bounds[name]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values)) if values else None)
        return dict(seam=seam,common_actual_tensor_rows=bounds,current_common_tensor_contribution_count=len(bounds),
            overlap_not_used_for_tensor_source_equality=True,actual_completed_tensor_and_divergence_remainder_common_source_bound=True,
            same_actual_complete_current_KR_units_histories_and_radius=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_full_flatten_power_normalization_AST_theorem=self.normalization,current_flatten_power_source_units_theorem=self.units,
            current_flatten_power_full_remaining_transport_theorems=self.transport,current_two_flatten_power_completed_tensor_function_join_theorems=self.joins,
            checked_original_full_stress_baseline_KR_homogeneity=self.heat_tensor.o7.entry_source.baselines,
            checked_original_general_K_physical_tensor_decomposition=self.heat_tensor.physical_proof,
            actual_current_tensor_regions_available=['flatten','outer_power']+self.heat_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=8,current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
            source_domain='flatten t[0,100],outer power phase[0,1],y=(Lrel-4)*phase;Z[-1,1],R>0,|Z|<1,tau>0,finite logtau,constant nu>0',
            variable_flatten_axial_K_and_complete_nonzero_A_E_P_preserved=True,
            remaining_pulse_core_global_cone_lift_NS_temporal_remainder_energy_points_n_recursion_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentFlattenPowerBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_flatten_power_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_flatten_power_tensor_views'][name]=field.chart(*args)
        print('Current actual flatten/power tensor: '+name,flush=True)
    result['current_actual_two_flatten_power_tensor_interface_bounds']={seam:field.interface(seam) for seam in SEAMS}
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

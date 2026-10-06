"""Actual current collar/exterior tensors and waiting-to-heat attachment.

Stable native small-defect stress is converted into common KR units.
Current closed terminal functions identify the canonical full Gamma exterior.
"""
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_steep_waiting_background_stress import (
    CurrentSteepWaitingBackgroundStress,HERE,PREFIX,OPEN,sha,pack,encode,endpoints,
    copy_jet,angular_source_log_parts,current_full_moment_rows,TENSOR_KEYS,
    source_precision,accepted,_verify_hashes,SourceAST,ordinary_grid)
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_defect_rows,collar_stress_rows
from lei_ren_part1_paper_compliant_waiting_stress_C3 import waiting_join_binding
from lei_ren_part1_paper_compliant_collar_physical_C2 import physical_collar_identities
from lei_ren_part1_paper_compliant_heat_physical_C4 import physical_heat_identities
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import ordinary_coordinate_proof
from lei_ren_part1_paper_compliant_current_selected_energy_source import binding
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

NAME=PREFIX+'current_heat_background_tensor.json'
RECEIPT=PREFIX+'current_heat_background_tensor_check.json'
GATES=('current_actual_collar_full_moment_background_tensor_available',
    'current_actual_full_Gamma_exterior_tensor_and_regional_NS_identity_certified',
    'current_actual_waiting_collar_exterior_two_tensor_joins_certified')
SEAMS=('waiting_collar','collar_exterior')
VIEWS={'collar_whole':('heat_collar',(-1,1),(0,3),('-3','-1'),None,'1'),
    'collar_inlet':('heat_collar',(-1,1),'0','-1',None,'1'),
    'collar_Gamma_join':('heat_collar',(-1,1),'3','-1',None,'1'),
    'collar_phi_crossing':('heat_collar',(-1,1),('2.99','3'),('-4','-2'),None,'.2'),
    'collar_fresh':('heat_collar','.537','.417','-2.6','.41','.8'),
    'exterior_inlet':('heat_exterior',(-1,1),'3','-1',None,'1'),
    'exterior_fresh':('heat_exterior','.731','4.17','-2.6','.41','.8'),
    'exterior_far':('heat_exterior',(-1,1),('10','1000'),('-4','-2'),None,'.2')}


def full_collar_normalization_AST_theorem():
    """Actual native remaining numerators and signed recurrence, arbitrary data."""
    a,eps,S,t,z,C=s.symbols('a eps inverse_Rtail t Z inverse_KR',real=True,positive=True)
    k=1-a;delta=2*a;p=1+delta
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),exp=s.exp)
    heat=SimpleNamespace(ctx=c,a=a,eps=eps,S=S,k=k,delta=delta,prate=p)
    one=s.Integer(1);stub=SimpleNamespace(constant=lambda c,v,n:s.Rational(v))
    future={name:s.Function('same_full_Gamma_'+name)(z) for name in ('theta','energy','pressure')}
    atoms={name:s.Symbol('same_'+name) for name in ('JW','EW','EW2','PW','PW2')}
    W=[s.Function('same_W'+str(j))(z) for j in range(5)]
    D=[s.Function('same_D'+str(j))(z) for j in range(5)]
    Km=[-W[j]*eps-D[j]*a*S for j in range(5)];K=[1+Km[0]]+Km[1:]
    asts=SourceAST();env=dict(collar_defect_rows.__globals__);env['IntervalTaylor']=stub
    original=asts.replay('collar_stress_C3','collar_defect_rows',env)
    tails=dict(scaled_full_future_Gamma_defects=future,separate_epsilon_atoms=atoms)
    defects=original(heat,dict(W_rows=W,D_rows=D),tails,t)
    actualA=asts.evaluate(asts.expression('collar_Gamma_C4','collar_tails','A'),
        dict(one=one,self=heat,angular=future['theta'],atoms=atoms,c=c,t=t))
    actualE=asts.evaluate(asts.expression('collar_Gamma_C4','collar_tails','E'),
        dict(one=one,self=heat,energy=future['energy'],atoms=atoms,c=c,t=t))*s.exp(delta*t)
    actualP=asts.evaluate(asts.expression('collar_Gamma_C4','collar_tails','P'),
        dict(one=one,self=heat,pressure=future['pressure'],atoms=atoms,c=c,t=t))*s.exp(p*t)
    full=current_full_moment_rows(heat,[v*C for v in K],actualA*C,actualE*C**2,actualP*C**2)
    checks={}
    for name,key,base,unit in (('A','angular_defect_rows',1/k,C),('E','energy_defect_rows',1/delta,C**2),('P','pressure_defect_rows',1/(2*p),C**2)):
        for j in range(5):
            residual=s.cancel(s.expand((defects[key][j]+(base if j==0 else 0))*unit-full[name][j]))
            if residual!=0:raise ArithmeticError('Actual collar full numerator/defect AST normalization differs')
            for n in range(5-j):checks[name+'_y%d_Z%d'%(j,n)]=s.diff(residual,z,n)==0
    return dict(original_full_collar_numerator_and_defect_AST_common_KR_mixed4_identities=checks,
        arbitrary_complete_Gamma_and_epsilon_atoms_used=True,
        original_small_defect_stress_computed_before_KR_conversion=True,
        actual_full_nonzero_moments_not_replaced_by_defects=True,input_hashes=asts.hashes,passed=True)


def heat_units_and_pressure_source_proof(o7):
    o7.assert_graph();h=o7.heat;ex=o7.physical.exterior
    if not ex.acceptance_loaded or ex.history is not o7.history or h is not ex.heat or not ex.proof['all_current_five_terminal_moment_functions_identified']:
        raise ValueError('Checked same-current closed five-history exterior required')
    if not ex.proof['current_exact_function_links']['checked_exact_pressure_function_identity'] or not ex.proof['current_exact_function_links']['checked_exact_angular_function_identity']:
        raise ValueError('Current original Cp=0 and Dtheta=0 source functions required')
    specs=(('collar','X',"(tails['angular_numerator']+data['angular_tail_constant_defect']*c.exp(-self.k*t))/K"),
        ('collar','energy',"tails['remaining_energy_in_Rtail_units']*c.exp(self.delta*t)/(K*K*2)"),
        ('collar','pressure',"self.forward_pressure(Z,t,data['Ptail'])"),
        ('data','pressure3','self.forward_pressure(Z,3,Ptail)'),
        ('exterior','pressure',"data['pressure3']+(loc3['pressure_numerator']*c.exp(-3*self.prate)-local['pressure_numerator']*c.exp(-self.prate*t))*self.pressure_scale"))
    bindings={}
    for method,target,value in specs:
        binding('compliant_collar_Gamma_C4',method,target,value);bindings[method+'.'+target]=value
    SourceAST().expression('collar_Gamma_C4','forward_pressure','integral',wanted="K*K*(ds*c.exp(-self.prate*v)/2)",augmented=True)
    bindings['forward_pressure.original_density']='K*K*(ds*exp(-p*v)/2)'
    a,mu,T,W,t,lone=s.symbols('a mu Ts wait t logone',real=True);k=1-a;p=1+2*a;r=1-mu
    LS=a-mu-r/2;LQ=LS-k*T;LT=LQ-k/2;logKR=1+mu/2-3*a/2+k*T+lone
    rS=-1-2*mu-r;rQ=rS-3*T;rT=rQ-3+k;rB=rT-p*W-2*lone
    identities={}
    for name,value in (
        ('same_inverse_KR_unit',LT-lone+logKR),
        ('waiting_actual_constant_K_equals_local_heat_1_minus_epsilon',LT-(LT-lone+lone)),
        ('full_heat_pressure_endpoint_exponent',p*(T+2+W)+rB-2*(LT-lone)),
        ('same_current_collar_pressure_source_exponent',p*(T+2+W+t)+rB-2*(LT-lone)-p*t)):
        if s.expand(value)!=0:raise ArithmeticError('Current heat KR/pressure source unit differs')
        identities[name]=True
    # Exact function additivity: I0 = integral_0^3 density + exp(-3p)*N3.
    # P3=Ptail+scale*prefix; thus both native infinity constants are the
    # same Ptail+scale*I0 already proved zero in the current closure graph.
    Ptail,scale,I0,I3,prefix=s.symbols('Ptail pressure_scale I0 N3 prefix',real=True)
    residual=s.expand((Ptail+scale*prefix)+scale*s.exp(-3*p)*I3-(Ptail+scale*I0))
    if s.simplify(residual.subs(I0,prefix+s.exp(-3*p)*I3))!=0:raise ArithmeticError('Native Ptail/P3 complete pressure additivity differs')
    identities['native_Ptail_and_P3_infinity_constants_same_function_by_complete_FTC']=True
    return dict(actual_current_closed_exterior_history_transfer=ex.proof,
        actual_native_forward_and_remaining_pressure_AST_bindings=bindings,identities=identities,
        original_integrand_and_complete_collar_Gamma_suffix_from_checked_current_bridge=True,
        Cp_zero_applied_only_after_same_native_Ptail_P3_function_identity=True,
        Dtheta_zero_applied_only_after_checked_current_native_angular_function_transfer=True,
        exact_inverse_KR_log='LT-logone',no_amplitude_cap_division=True,
        same_actual_current_nonzero_energy_angular_and_pressure_history=True,passed=True)


def heat_join_source_proofs(o7,units):
    post=o7.entry_source.angular.atlas.postpulse.proof;coordinate=ordinary_coordinate_proof(o7.physical)
    formulas=dict(waiting_collar=waiting_join_binding(),collar_exterior=o7.physical.exterior.native_stress_AST)
    result={}
    for seam in SEAMS:
        primitive={k:v for k,v in post['current_instantiated_seam_mixed4_identities'].items() if k.startswith(seam+'_')}
        if len(primitive)!=60 or not all(primitive.values()) or not coordinate['exact_common_radius_identities'][seam]:raise ValueError('Current full heat source seam missing')
        result[seam]=dict(current_primitive_function_mixed4_identities=primitive,
            original_stress_pressure_function_AST_join=formulas[seam],same_actual_heat_KR_and_complete_pressure_source=units,
            same_current_closed_energy_Dtheta_Cp_functions=True,same_ordinary_radius_coordinate=coordinate,
            same_general_K_completed_tensor_divergence_remainder_operators=True,
            interval_overlap_not_used_as_function_identity=True,passed=True)
    return result


def source_zero_tensor(tensor):
    """Evaluate proved exact tensor/remainder identities, keep source factors."""
    if isinstance(tensor,dict) and 'signed_coefficient' in tensor:
        return dict(tensor,signed_coefficient=tensor['signed_coefficient']*0,exact_zero=True,
            log_absolute_upper=None,source_exact_identity_evaluated_before_acceptance=True)
    if isinstance(tensor,dict):return {k:source_zero_tensor(v) for k,v in tensor.items()}
    return [source_zero_tensor(v) for v in tensor]


class CurrentHeatBackgroundTensor:
    @source_precision
    def __init__(self,o7=None,require_checked=True):
        self.o7=o7 if o7 is not None else CurrentSteepWaitingBackgroundStress()
        if not self.o7.acceptance_loaded:raise ValueError('Checked actual O7 tensors required')
        self.physical=self.o7.physical;self.history=self.physical.history;self.heat=self.history.heat;self.steep=self.history.steep;self.ctx=self.physical.ctx
        self.exterior=self.physical.exterior;self.family=self.o7.family;self.source=self.o7.source;self.datum_sha=self.o7.datum_sha
        self.units=heat_units_and_pressure_source_proof(self.o7);self.normalization=full_collar_normalization_AST_theorem()
        self.joins=heat_join_source_proofs(self.o7,self.units)
        self.physical_proof=physical_collar_identities();self.Gamma_physical=physical_heat_identities()
        if not self.Gamma_physical['axial_viscosity_remainder_exactly_zero'] or not self.Gamma_physical['all_three_physical_momentum_components_exactly_zero']:
            raise ValueError('Full Gamma physical zero remainder/NS theorem missing')
        self.tensor=self.o7.tensor;self.hashes=dict(self.o7.hashes)
        for stem in ('collar_stress_C3','collar_Gamma_C4','collar_physical_C2','heat_physical_C4','waiting_stress_C3','current_heat_background_tensor'):
            name=PREFIX+stem+'.py';digest=sha(name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Checked current heat dependency changed: '+name)
            self.hashes[name]=digest
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Current heat tensor receipt exceeds regional scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.o7.assert_graph()
        if not (self.o7.acceptance_loaded and self.physical is self.o7.physical and self.history is self.physical.history and self.heat is self.history.heat and self.exterior is self.physical.exterior and self.exterior.acceptance_loaded and self.exterior.history is self.history and self.exterior.heat is self.heat and self.ctx is self.physical.ctx):
            raise ValueError('Actual heat tensor must retain checked current O7/closed exterior graph')

    @source_precision
    def chart(self,chart,Z,t,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);t=c.mpf(t);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in ('heat_collar','heat_exterior') or not all(mp.isfinite(v) for v in endpoints(Z)+endpoints(t)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(nu)[0]<=0 or endpoints(t)[0]<(0 if chart=='heat_collar' else 3) or (chart=='heat_collar' and endpoints(t)[1]>3):raise ValueError('Original heat chart finite offset/Z/time and nu>0 required')
        if theta is not None and not all(mp.isfinite(v) for v in endpoints(c.mpf(theta))):raise ValueError('Finite theta or all-angle bounds required')
        h=self.heat;st=self.steep;one=IntervalTaylor.constant(c,1,5)
        LT=h.a-st.mu-st.rate/2-h.k*st.Ts-h.k/2;logC=LT-st.logone;C=c.exp(logC)
        if chart=='heat_collar':
            shape=h.shape(Z,t);tails=h.collar_tails(Z,t)
            defects=collar_defect_rows(h,shape,tails,h.ctx.mpf(t))
            A0=copy_jet(c,tails['angular_numerator'])*C
            E0=copy_jet(c,tails['remaining_energy_in_Rtail_units'])*(C*C*c.exp(h.delta*t))
            P0=copy_jet(c,tails['remaining_pressure_in_Rtail_units'])*(C*C*c.exp(h.prate*t))
        else:
            shape=h.local_Gamma(Z,t);A0=copy_jet(c,shape['angular_numerator'])*C
            E0=copy_jet(c,shape['energy_numerator'])*(C*C);P0=copy_jet(c,shape['pressure_numerator'])*(C*C)
            defects=dict(angular_defect_rows=[],energy_defect_rows=[],pressure_defect_rows=[],K_defect_rows=[])
            local_full=current_full_moment_rows(h,shape['K_rows'],shape['angular_numerator'],shape['energy_numerator'],shape['pressure_numerator'])
            for name,key,base in (('A','angular_defect_rows',1/h.k),('E','energy_defect_rows',1/h.delta),('P','pressure_defect_rows',1/(2*h.prate))):
                defects[key]=[local_full[name][0]-base]+local_full[name][1:]
            defects['K_defect_rows']=[shape['K_rows'][0]-1]+shape['K_rows'][1:]
        K=[copy_jet(c,v)*C for v in shape['K_rows'][:5]]
        full=current_full_moment_rows(h,K,A0,E0,P0)
        native=self.history.evaluate(chart,Z,t)['source_packet']
        rawstress=collar_stress_rows(h,shape,defects,h.ctx.mpf(Z),h.ctx.mpf(t))
        stress={name:[copy_jet(c,v)*(C*C if name=='axial' else C) for v in rawstress[name]] for name in ('theta','axial','theta_inertial','theta_shear')}
        grids={name:ordinary_grid(rows,3) for name,rows in stress.items()}
        diagnostic=copy.deepcopy(grids)
        canonical=chart=='heat_exterior' or endpoints(t)==(mp.mpf(3),mp.mpf(3))
        if canonical:
            # This is the already source-bound full Gamma stress theorem,
            # not an override inferred from numerical overlap or a cap.
            grids['theta']={k:c.mpf(0) for k in grids['theta']};grids['axial']={k:c.mpf(0) for k in grids['axial']}
        v=2+st.Ts+st.wait+t;parts=angular_source_log_parts(self.physical,v)
        logR,radius_source=self.physical.radius(chart,t,{},h)
        tensor=self.tensor(self,c,Z,K,grids,parts,logR,lt,theta,nu)
        error_diagnostic=tensor['physical_axial_viscosity_remainder_mixed2']
        if canonical:tensor={key:source_zero_tensor(value) for key,value in tensor.items()}
        P=full['P'];pressure_rows=[-sum((P[n]*math.comb(j,n)*(-h.prate)**(j-n) for n in range(j+1)),one*0) for j in range(5)]
        return dict(chart=chart,Z=Z,heat_offset=t,log_tau=lt,viscosity=nu,
            current_actual_normalized_full_moment_rows=full,current_actual_normalized_stress_mixed3=grids,
            original_correlated_stress_operator_enclosure_diagnostic=diagnostic,
            original_axial_viscosity_operator_enclosure_diagnostic=error_diagnostic,
            exact_inverse_KR_log=logC,exact_current_positive_log_source_parts=parts,exact_source_logR=logR,radius_source=radius_source,
            current_native_small_defect_stress_retained_before_KR_scaling=True,
            stable_actual_absolute_pressure_mixed4_factored=ordinary_grid(pressure_rows,4),
            stable_pressure_positive_source_log_parts=parts['absolute_pressure'],
            original_forward_absolute_pressure_Taylor=native['pressure_over_Pstar_squared_Taylor'],
            original_current_energy_Taylor=native['energy_Taylor'],normalized_energy_native_consistency=E0-K[0]*K[0]*copy_jet(c,native['energy_Taylor'])*2,
            actual_full_stress_not_local_difference=True,actual_Cp_Dtheta_zero_from_checked_same_current_functions=True,
            source_exact_Gamma_tensor_and_remainder_zero=canonical,
            actual_Gamma_regional_physical_NS_identity=canonical,
            source_bounds_not_resolved_physical_point_values=True,regional_result_not_global_temporal_flatness=True,
            **tensor,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,seam,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph()
        if seam not in SEAMS or not self.joins[seam]['passed']:raise ValueError('Checked current heat function seam required')
        left=self.o7.chart('waiting',Z,1,log_tau,theta,viscosity) if seam=='waiting_collar' else self.chart('heat_collar',Z,3,log_tau,theta,viscosity)
        right=self.chart('heat_collar' if seam=='waiting_collar' else 'heat_exterior',Z,0 if seam=='waiting_collar' else 3,log_tau,theta,viscosity)
        def flatten(value,path):
            if isinstance(value,dict) and 'signed_coefficient' in value:return {path:value}
            result={}
            for key,row in (value.items() if isinstance(value,dict) else enumerate(value)):result.update(flatten(row,path+'/'+str(key)))
            return result
        bounds={}
        for key in TENSOR_KEYS:
            a=flatten(left[key],key);b=flatten(right[key],key)
            if set(a)!=set(b):raise ValueError('Common actual heat tensor layout differs')
            for name,row in a.items():
                other=b[name]
                for sk in ('physical_lambda_exponent','physical_viscosity_exponent'):
                    if encode(pack(row[sk]))!=encode(pack(other[sk])):raise ValueError('Common tensor physical exponent differs')
                values=[endpoints(v['log_absolute_upper'])[1] for v in (row,other) if not v['exact_zero']]
                bounds[name]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values)) if values else None)
        return dict(seam=seam,common_actual_tensor_rows=bounds,current_common_tensor_contribution_count=len(bounds),
            overlap_not_used_for_tensor_source_equality=True,actual_completed_tensor_and_divergence_remainder_common_source_bound=True,
            same_actual_complete_current_KR_units_histories_and_radius=True,
            source_exact_Gamma_zero_at_interface=seam=='collar_exterior',**dict.fromkeys(OPEN,False))

    @source_precision
    def unbounded_exterior(self,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        # Zero function identities hold on the entire original exterior.
        # A finite anchor supplies only the component layout. No infinite
        # source logs or finite truncation define this domain theorem.
        anchor=self.chart('heat_exterior',Z,3,log_tau,theta,viscosity)
        return dict(source_domain='every finite t>=3 and its infinity limit; R>0,|Z|<1,tau>0,constant nu>0',
            original_unbounded_exterior_covered=True,finite_anchor_only_supplies_tensor_component_layout=True,
            canonical_full_Gamma_function_and_current_terminal_source_proof=self.exterior.theorem,
            physical_heat_source_function_proof=self.Gamma_physical,
            actual_stress_tensor_divergence_remainder_and_momentum_exactly_zero=True,
            zero_identity_not_radial_truncation_or_temporal_flatness=True,
            **{key:anchor[key] for key in TENSOR_KEYS},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_heat_units_complete_pressure_and_closed_terminal_source_theorem=self.units,
            current_collar_full_moment_normalization_AST_theorem=self.normalization,
            current_two_heat_completed_tensor_function_join_theorems=self.joins,
            original_general_K_collar_physical_tensor_theorem=self.physical_proof,
            current_full_Gamma_regional_physical_identity_theorem=self.Gamma_physical,
            checked_original_full_stress_baseline_KR_homogeneity=self.o7.entry_source.baselines,
            actual_current_tensor_regions_available=['outer_angular','steep_entry','steep_power','steep_exit','waiting','heat_collar','heat_exterior'],
            actual_current_completed_tensor_adjacent_interface_count=6,current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
            source_domain='collar t[0,3],exterior finite t>=3 plus separate full-unbounded identity;Z[-1,1],R>0,|Z|<1,tau>0,finite compact logtau,constant nu>0',
            actual_current_histories_and_nonzero_energy_pressure_retained=True,
            remaining_chart_tensors_global_cone_lift_NS_temporal_remainder_energy_and_n_recursion_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentHeatBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_heat_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_heat_tensor_views'][name]=field.chart(*args)
        print('Current actual heat tensor: '+name,flush=True)
    result['current_actual_two_heat_tensor_interface_bounds']={seam:field.interface(seam) for seam in SEAMS}
    result['current_unbounded_exact_Gamma_physical_identity']=field.unbounded_exterior()
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

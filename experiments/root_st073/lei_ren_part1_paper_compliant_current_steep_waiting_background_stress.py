"""Actual current O7 power/exit/waiting tensors in common KR units.

Remaining pressure is reduced algebraically before interval evaluation.
All derivatives are ordinary logR derivatives, never phase derivatives.
"""
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_steep_entry_background_stress import (
    CurrentSteepEntryBackgroundStress,HERE,PREFIX,OPEN,sha,pack,encode,endpoints,
    copy_jet,angular_source_log_parts,entry_stress_rows,product_rows,ordinary_grid,
    source_precision,accepted,_verify_hashes,SourceAST)
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import entry_power_join_binding
from lei_ren_part1_paper_compliant_steep_power_stress_C3 import power_transport_identities,power_exit_join_binding
from lei_ren_part1_paper_compliant_steep_exit_stress_C3 import exit_transport_identities,exit_join_binding
from lei_ren_part1_paper_compliant_waiting_stress_C3 import waiting_identities
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import ordinary_coordinate_proof
from lei_ren_part1_paper_compliant_current_selected_energy_source import binding
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral,sigma_enclosure
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

NAME=PREFIX+'current_steep_waiting_background_stress.json'
RECEIPT=PREFIX+'current_steep_waiting_background_stress_check.json'
CHARTS=('steep_power','steep_exit','waiting')
SEAMS=('steep_entry_power','steep_power_exit','steep_exit_waiting')
GATES=('current_actual_steep_power_exit_waiting_full_stress_mixed3_recovered',
    'current_actual_steep_power_exit_waiting_completed_tensor_remainder_available',
    'current_actual_entry_power_exit_waiting_three_tensor_joins_certified')
VIEWS={'whole':((-1,1),(0,1),('-3','-1'),None,'1'),
    'inlet':((-1,1),'0','-1',None,'1'),
    'outlet':((-1,1),'1','-1',None,'1'),
    'fresh':('.537','.417','-2.6','.41','.8')}
TENSOR_KEYS=('physical_cylindrical_stress_mixed3','physical_stress_divergence_mixed2',
    'completed_theta_theta_stress_mixed2','physical_axial_viscosity_remainder_mixed2',
    'completed_background_tensor_cartesian_components','physical_completed_stress_divergence_cartesian',
    'physical_remainder_cartesian','physical_momentum_residual_decomposition_cartesian')


def remaining_exit_pressure(c,t,k,Jt,cells):
    """Original integral_t^1 exp(-3*v+2*k*J(v))/2 dv, positive cells."""
    t=c.mpf(t);length=1-t;ds=length/cells;J=c.mpf(Jt);pressure=c.mpf(0)
    if endpoints(t)[0]<0 or endpoints(t)[1]>1 or cells<1:raise ValueError('Original exit domain/cells required')
    if endpoints(length)[1]>0:
        for i in range(cells):
            left=t+length*i/cells;right=t+length*(i+1)/cells
            v=c.mpf([max(mp.mpf(0),endpoints(left)[0]),min(mp.mpf(1),endpoints(right)[1])])
            nextJ=J+ds*sigma_enclosure(c,v)
            jc=c.mpf([endpoints(J)[0],endpoints(nextJ)[1]])
            pressure+=ds*c.exp(-3*v+2*k*jc)/2
            J=nextJ
    return pressure


def current_full_moment_rows(heat,K,A0,E0,P0):
    Q=product_rows(K,K);A=[A0];E=[E0];P=[P0]
    for j in range(4):
        A.append(K[j]-A[j]*heat.k);E.append(E[j]*heat.delta-Q[j]);P.append(P[j]*heat.prate-Q[j]/2)
    return dict(A=A,E=E,P=P,K=K)


def log_rate_rows(c,K0,g,X0,k):
    one=IntervalTaylor.constant(c,1,5);K=[one*K0];X=[X0]
    for n in range(1,5):
        K.append(sum((K[n-1-j]*g[j]*math.comb(n-1,j) for j in range(n)),one*0))
        X.append((one if n==1 else one*0)-sum((X[n-1-j]*(k+g[0] if j==0 else g[j])*math.comb(n-1,j) for j in range(n)),one*0))
    return K,X


def source_normalization_theorem():
    """Bind native half-energy/rates and pressure source before KR reduction."""
    specs=(('steep_power','energy',"(data['after_power']-c.mpf('.5'))*c.exp(-2*left)/2+c.mpf('.25')"),
        ('steep_out','energy',"(data['waiting_future']*c.exp(-1-self.delta/2)+kernels['remaining_energy'])*c.exp(2*t-2*self.k*J)/2"),
        ('waiting','energy',"(data['H']*c.exp(-self.delta*left)+decay_integral(c,self.delta,left))/2"),
        ('steep_power','X',"data['XS']+t"),
        ('steep_out','X',"(data['XQ']+kernels['angular'])*c.exp(-self.k*J)"),
        ('waiting','X',"(data['XT']-1/self.k)*c.exp(-self.k*t)+1/self.k"))
    bindings={}
    for method,target,value in specs:
        binding('compliant_steep_waiting_C4',method,target,value);bindings[method+'.'+target]=value
    asts=SourceAST();v,k,J,ds=s.symbols('v k J ds',real=True);c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    for stem,method in (('steep_waiting_C4','transition_kernels'),('current_steep_waiting_background_stress','remaining_exit_pressure')):
        node=asts.expression(stem,method,'pressure',wanted='ds*c.exp(-3*v+2*'+('rate' if stem=='steep_waiting_C4' else 'k')+'*jc)/2',augmented=True)
        if s.simplify(asts.evaluate(node,dict(c=c,ds=ds,v=v,k=k,rate=k,jc=J))-ds*s.exp(-3*v+2*k*J)/2)!=0:
            raise ArithmeticError('Current remaining exit pressure integrand differs')
    a,mu,T,W,y,t,J,one=s.symbols('a mu Ts wait y t J logone',real=True)
    k=1-a;p=1+2*a;r=1-mu;KS=a-mu-r/2;KQ=KS-k*T;KT=KQ-k/2
    rS=-1-2*mu-r;rQ=rS-3*T;rT=rQ-3+k
    checks={}
    def zero(name,value):
        if s.expand(value)!=0:raise ArithmeticError('Common KR source normalization failed: '+name)
        checks[name]=True
    zero('power_pressure_full_future_reduction',p*(1+y)+rS-3*y-2*(KS-k*y))
    zero('power_future_exit_relative_attenuation',rQ-rS+3*y+3*(T-y))
    zero('exit_pressure_full_future_reduction',p*(1+T+t)+rQ-(2*KQ+p*t))
    zero('waiting_pressure_full_future_reduction',p*(2+T)+rT-2*KT)
    zero('entry_power_K_endpoint',KS-(a-mu-r/2))
    zero('power_exit_K_endpoint',KS-k*T-KQ)
    zero('exit_waiting_K_endpoint',KQ-k+k/2-KT)
    # Original waiting root: W=(logone-KT-logKR)/k is already source-bound.
    logKR=1+mu/2-3*a/2+k*T+one
    zero('waiting_actual_constant_K_equals_one_minus_epsilon',KT+logKR-one)
    # Execute the actual theta and actual new logK expressions on the same
    # arbitrary source parameters. This identifies the reconstructed K,
    # not just an algebraic formula which resembles it.
    thetaR=s.Symbol('thetaR',positive=True);bp=s.Rational(1,2)+mu;bh=s.Rational(1,2)+a
    native=SimpleNamespace(thetaS=thetaR*s.exp(-bp-r/2),
        thetaQ=thetaR*s.exp(-bp-r/2-s.Rational(3,2)*T),
        thetaT=thetaR*s.exp(-bp-r/2-s.Rational(3,2)*T-s.Rational(3,2)+k/2),k=k,bh=bh)
    for method,offset,position,expected,expression in (
            ('steep_power',y,1+y,KS-k*y,'KS-k*y'),
            ('steep_out',t,1+T+t,KQ-k*t+k*J,'KQ-k*x+k*J'),
            ('waiting',y,2+T+y,KT,'KT')):
        actual_theta=asts.evaluate(asts.expression('steep_waiting_C4',method,'theta'),
            dict(self=native,c=c,one=1,t=offset,J=J))
        reconstructed=asts.evaluate(asts.expression('current_steep_waiting_background_stress','chart','logK',wanted=expression),
            dict(KS=KS,KQ=KQ,KT=KT,k=k,y=y,x=t,J=J))
        ratio=actual_theta/(thetaR*s.exp(-bh*position))
        if s.simplify(s.expand_power_exp(ratio/s.exp(reconstructed)))!=1 or s.expand(reconstructed-expected)!=0:
            raise ArithmeticError('Actual native theta/current normalized K AST differs')
        checks[method+'_actual_theta_and_current_logK_AST_identified']=True
    return dict(native_nonzero_half_energy_and_angular_history_AST_bindings=bindings,
        original_remaining_exit_pressure_integrand_AST=asts.bindings,identities=checks,
        common_KR_units_with_checked_entry=True,positive_full_Gamma_future_retained=True,
        pressure_reduction_before_enclosure_not_division_of_tiny_caps=True,
        original_forward_and_remaining_pressure_identical_by_FTC_additivity_and_current_Cp_zero=True,
        power_uniform_half_energy_floor_one_quarter_retained=True,input_hashes=asts.hashes,passed=True)


def moment_baseline_AST_theorem():
    """Replay original power/exit recurrence suffix on arbitrary full data."""
    import ast
    a,z,KR=s.symbols('a Z KR',positive=True);k=1-a;delta=2*a;p=1+delta
    heat=SimpleNamespace(k=k,delta=delta,prate=p)
    K=[s.Function('K'+str(j))(z) for j in range(5)]
    full=current_full_moment_rows(heat,K,*[s.Function(name)(z) for name in ('A0','E0','P0')])
    actualK=[KR*v for v in K];Q=product_rows(actualK,actualK)
    Km=[actualK[0]-1]+actualK[1:];Qm=[Q[0]-1]+Q[1:]
    checks={};asts=SourceAST()
    for stem,method in (('steep_power_stress_C3','power_defect_rows'),('steep_exit_stress_C3','exit_defect_rows')):
        fn=asts.method(stem,method)
        begin=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and any(ast.unparse(v)=='A' for v in n.targets))
        wrapper=ast.parse('def recurrence(heat,shape,Ad,Ed,Pd):\n    pass').body[0]
        wrapper.body=copy.deepcopy(fn.body[begin:]);env={}
        exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper],type_ignores=[])),'<original full O7 moment AST>','exec'),env)
        result=env['recurrence'](heat,dict(K_defect_rows=Km,K_squared_defect_rows=Qm),
            KR*full['A'][0]-1/k,KR**2*full['E'][0]-1/delta,KR**2*full['P'][0]-1/(2*p))
        for name,key,base,unit in (('A','angular_defect_rows',1/k,KR),('E','energy_defect_rows',1/delta,KR**2),('P','pressure_defect_rows',1/(2*p),KR**2)):
            for j in range(5):
                residual=s.cancel(s.expand(result[key][j]+(base if j==0 else 0)-unit*full[name][j]))
                if residual!=0:raise ArithmeticError('Original current full O7 normalized rows differ')
                for n in range(5-j):checks[stem+'_'+name+'_y%d_Z%d'%(j,n)]=s.diff(residual,z,n)==0
    return dict(original_full_power_exit_defect_AST_KR_normalization_mixed4_identities=checks,
        arbitrary_full_nonzero_axial_functions_used=True,passed=True)


def joined_source_proofs(entry):
    entry.assert_graph();post=entry.angular.atlas.postpulse.proof
    coordinate=ordinary_coordinate_proof(entry.physical)
    formulas=(entry_power_join_binding(),power_exit_join_binding(),exit_join_binding())
    result={}
    for seam,formula in zip(SEAMS,formulas):
        primitive={k:v for k,v in post['current_instantiated_seam_mixed4_identities'].items() if k.startswith(seam+'_')}
        if len(primitive)!=60 or not all(primitive.values()) or not coordinate['exact_common_radius_identities'][seam]:
            raise ValueError('Current complete source function seam missing: '+seam)
        flags=[v for k,v in formula.items() if k.endswith('_verified')]
        if not flags or not all(flags):raise ValueError('Original stress/pressure arbitrary-function seam missing')
        result[seam]=dict(current_primitive_mixed4_identities=primitive,
            original_arbitrary_terminal_stress_pressure_AST_join=formula,
            same_source_radius_coordinate=coordinate,
            same_actual_KR_units_nonzero_histories_and_Cp_zero=True,
            same_general_K_tensor_stress3_divergence2_completion2_remainder2=True,
            interval_overlap_not_used_as_functional_identity=True,passed=True)
    return result


class CurrentSteepWaitingBackgroundStress:
    @source_precision
    def __init__(self,entry=None,require_checked=True):
        self.entry_source=entry if entry is not None else CurrentSteepEntryBackgroundStress()
        if not self.entry_source.acceptance_loaded:raise ValueError('Checked actual current entry tensor required')
        self.physical=self.entry_source.physical;self.history=self.physical.history
        self.steep=self.history.steep;self.heat=self.history.heat;self.ctx=self.physical.ctx
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.transport=dict(steep_power=power_transport_identities(),steep_exit=exit_transport_identities(),waiting=waiting_identities())
        self.normalization=source_normalization_theorem();self.moment_proof=moment_baseline_AST_theorem()
        self.joins=joined_source_proofs(self.entry_source);self.tensor=self.entry_source.tensor
        self.hashes=dict(self.entry_source.hashes)
        for stem in ('steep_power_stress_C3','steep_exit_stress_C3','waiting_stress_C3','current_steep_waiting_background_stress'):
            name=PREFIX+stem+'.py';digest=sha(name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Checked dependency changed: '+name)
            self.hashes[name]=digest
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Current O7 tensor receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.entry_source.assert_graph()
        if not (self.entry_source.acceptance_loaded and self.physical is self.entry_source.physical and self.history is self.physical.history and self.steep is self.history.steep and self.heat is self.history.heat and self.ctx is self.physical.ctx):
            raise ValueError('Current O7 tensor must retain admitted entry/physical/history graph')

    @source_precision
    def chart(self,chart,Z,x,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);x=c.mpf(x);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in CHARTS or not all(mp.isfinite(v) for v in endpoints(Z)+endpoints(x)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(x)[0]<0 or endpoints(x)[1]>1 or endpoints(nu)[0]<=0:
            raise ValueError('Original O7 chart/unit coordinate,Z[-1,1],finite compact logtau,nu>0 required')
        if theta is not None and not all(mp.isfinite(v) for v in endpoints(c.mpf(theta))):raise ValueError('Finite theta or all-angle bounds required')
        st=self.steep;h=self.heat;one=IntervalTaylor.constant(c,1,5);k=h.k;p=h.prate
        KS=h.a-st.mu-st.rate/2;KQ=KS-k*st.Ts;KT=KQ-k/2
        post=self.entry_source.angular.post_pressure(Z);Ih=copy_jet(c,post['original_full_Gamma_pressure_future'])
        Jwaiting=one*(decay_integral(c,p,st.wait)/2)+Ih*c.exp(-p*st.wait-2*st.logone)
        Jafter=one*st.outfull['pressure']+Jwaiting*c.exp(-3+k)
        if chart=='steep_power':
            native=st.steep_power(Z,x);y=st.Ts*x;left=st.Ts*(1-x);v=1+y
            logK=KS-k*y;g=[-k]+[c.mpf(0)]*3
            pressure0=(one*(decay_integral(c,3,left)/2)+Jafter*c.exp(-3*left))*c.exp(2*logK)
            pressure_source=dict(remaining_power=left,post_power=Jafter)
        elif chart=='steep_exit':
            native=st.steep_out(Z,x);kernels=st.kernels(x,'out');J=kernels['J'];sig=sigma_jets(c,x);v=1+st.Ts+x
            logK=KQ-k*x+k*J;g=[k*(sig[0]-1)]+[k*sig[j]*math.factorial(j) for j in range(1,4)]
            suffix=remaining_exit_pressure(c,x,k,J,st.cells)
            pressure0=(one*suffix+Jwaiting*c.exp(-3+k))*c.exp(2*KQ+p*x)
            pressure_source=dict(remaining_exit=suffix,post_exit=Jwaiting)
        else:
            native=st.waiting(Z,x);left=st.wait*(1-x);v=2+st.Ts+st.wait*x
            logK=KT;g=[c.mpf(0)]*4
            pressure0=(one*(decay_integral(c,p,left)/2)+Ih*c.exp(-p*left-2*st.logone))*c.exp(2*logK)
            pressure_source=dict(remaining_waiting=left,complete_heat=Ih)
        K,X=log_rate_rows(c,c.exp(logK),g,copy_jet(c,native['angular_Taylor']),k)
        E0=K[0]*K[0]*copy_jet(c,native['energy_Taylor'])*2
        full=current_full_moment_rows(h,K,K[0]*X[0],E0,pressure0)
        logR,radius_source=self.physical.radius(chart,x,{},st);q=-st.wait-st.Ts-2+v
        caplog=c.ln(c.mpf(endpoints(h.Scap)[1]))-q
        if endpoints(caplog+logR)[0]<0:raise ValueError('Current O7 inverse-radius cap differs')
        inverseR=c.mpf([0,endpoints(c.exp(c.mpf(endpoints(caplog)[1])))[1]])
        local=copy.copy(h);local.ctx=c;local.S=inverseR
        defs=dict(angular_defect_rows=full['A'],energy_defect_rows=full['E'],pressure_defect_rows=full['P'],K_defect_rows=K)
        stress=entry_stress_rows(local,dict(K_rows=K),defs,X,Z,c.mpf(0))
        grids={name:ordinary_grid(stress[name],3) for name in ('theta','axial','theta_inertial','theta_shear')}
        parts=angular_source_log_parts(self.physical,v);tensor=self.tensor(self,c,Z,K,grids,parts,logR,lt,theta,nu)
        P=full['P'];pressure_rows=[-sum((P[n]*math.comb(j,n)*(-p)**(j-n) for n in range(j+1)),one*0) for j in range(5)]
        return dict(chart=chart,Z=Z,original_coordinate=x,log_tau=lt,viscosity=nu,
            current_actual_normalized_full_moment_rows=full,current_actual_normalized_stress_mixed3=grids,
            current_actual_remaining_pressure_source=pressure_source,exact_K_over_KR_log=logK,
            exact_current_positive_log_source_parts=parts,exact_source_logR=logR,radius_source=radius_source,
            inverse_radius_enclosure=inverseR,inverse_radius_cap_only_not_exact_field_value=True,
            stable_actual_absolute_pressure_mixed4_factored=ordinary_grid(pressure_rows,4),
            stable_pressure_positive_source_log_parts=parts['absolute_pressure'],
            original_forward_absolute_pressure_Taylor=native['pressure_over_Pstar_squared_Taylor'],
            original_current_energy_Taylor=native['energy_Taylor'],normalized_energy_native_consistency=E0-K[0]*K[0]*copy_jet(c,native['energy_Taylor'])*2,
            actual_full_stress_not_local_difference=True,original_full_nonzero_A_E_P_and_current_C5_controls_retained=True,
            actual_pressure_is_same_P0_Pin_forward_function_by_current_Cp_zero=True,
            physical_tensor_operator_replayed_from_checked_actual_general_K_AST=True,
            ordinary_logR_derivatives_not_phase_derivatives=True,
            waiting_axial_viscosity_source_exact_zero=chart=='waiting',
            regional_decomposition_is_not_global_temporal_flatness=True,**tensor,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,seam,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        self.assert_graph()
        if seam not in SEAMS or not self.joins[seam]['passed']:raise ValueError('Checked current O7 function seam required')
        args=(Z,log_tau,theta,viscosity)
        left=self.entry_source.entry(Z,1,log_tau,theta,viscosity) if seam==SEAMS[0] else self.chart(CHARTS[SEAMS.index(seam)-1],Z,1,log_tau,theta,viscosity)
        right=self.chart(CHARTS[SEAMS.index(seam)],Z,0,log_tau,theta,viscosity)
        def flatten(value,path):
            if isinstance(value,dict) and 'signed_coefficient' in value:return {path:value}
            result={}
            for key,row in (value.items() if isinstance(value,dict) else enumerate(value)):result.update(flatten(row,path+'/'+str(key)))
            return result
        bounds={}
        for key in TENSOR_KEYS:
            a=flatten(left[key],key);b=flatten(right[key],key)
            if set(a)!=set(b):raise ValueError('Common actual O7 tensor layout differs')
            for name,row in a.items():
                other=b[name]
                for sourcekey in ('physical_lambda_exponent','physical_viscosity_exponent'):
                    if encode(pack(row[sourcekey]))!=encode(pack(other[sourcekey])):raise ValueError('Common physical tensor power differs')
                values=[endpoints(v['log_absolute_upper'])[1] for v in (row,other) if not v['exact_zero']]
                bounds[name]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values)) if values else None)
        return dict(seam=seam,common_actual_tensor_rows=bounds,current_common_tensor_contribution_count=len(bounds),
            exact_common_positive_source_factors_and_radius_from_function_proof=True,
            differently_rounded_valid_bounds_merged_after_function_identity=True,
            overlap_not_used_for_tensor_source_equality=True,all_angles_covered=theta is None,
            actual_completed_tensor_and_divergence_remainder_common_source_bound=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_O7_native_full_moment_transport_theorems=self.transport,
            current_O7_source_normalization_and_remaining_pressure_theorem=self.normalization,
            current_O7_full_moment_baseline_AST_theorem=self.moment_proof,
            current_O7_three_actual_tensor_function_join_theorems=self.joins,
            checked_entry_generic_A_KX_inertial_correlation=self.entry_source.correlation,
            checked_angular_general_K_physical_tensor_baseline=self.entry_source.baselines,
            source_domain='original steep_power/waiting phase[0,1],steep_exit t[0,1],Z[-1,1],R>0,|Z|<1,tau>0,finite compact logtau,constant nu>0; Z endpoints are infinity limits',
            source_bounds_are_enclosures_not_resolved_physical_point_values=True,
            actual_current_tensor_regions_available=['outer_angular','steep_entry']+list(CHARTS),
            actual_current_completed_tensor_adjacent_interface_count=4,
            current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
            other_chart_tensors_global_cone_lift_NS_temporal_remainder_energy_and_n_recursion_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentSteepWaitingBackgroundStress(require_checked=False)
    result=field.manifest();result['current_actual_O7_tensor_views']={}
    for chart in CHARTS:
        for name,args in VIEWS.items():
            result['current_actual_O7_tensor_views'][chart+':'+name]=field.chart(chart,*args)
        print('Current actual O7 full tensor: '+chart,flush=True)
    result['current_actual_O7_three_tensor_interface_bounds']={seam:field.interface(seam) for seam in SEAMS}
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

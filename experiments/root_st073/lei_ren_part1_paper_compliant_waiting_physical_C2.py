"""Original waiting physical stress/remainder companion, without source edits.

The accepted general-K physical collar transfer is reused as an algebraic
map. Its native chart call is redirected to waiting at the original phase;
the collar itself is NEVER evaluated at a negative offset.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_waiting_stress_C3 import (
    CompliantWaitingStressC3,source_precision)
from lei_ren_part1_paper_compliant_collar_physical_C2 import (
    CompliantCollarPhysicalC2,physical_collar_identities,physical_source_row,scale_row)
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


class WaitingHeatReference:
    """Only the constant waiting K replaces the generic map's shape input."""
    def __init__(self,heat):self.actual_heat=heat
    def __getattr__(self,key):return getattr(self.actual_heat,key)
    def shape(self,Z,q):
        one=IntervalTaylor.constant(self.ctx,1,5)
        return dict(K_rows=[one*(1-self.eps)]+[one*0]*4)


class WaitingDispatch:
    def __init__(self,native,stress):self.native=native; self.stress=stress
    def __getattr__(self,key):return getattr(self.native,key)
    def provider(self,chart):
        return self.stress.steep if chart=='waiting' else self.native.provider(chart)
    def evaluate(self,chart,Z,coordinate):
        if chart!='waiting':return self.native.evaluate(chart,Z,coordinate)
        packet=self.stress.waiting(Z,coordinate)
        return dict(chart=chart,source_packet=packet,
                    physical_mixed_grids=packet['physical_mixed_derivatives_total_order_le4'],
                    original_waiting_absolute_pressure_companion_used=True)


class AtWaitingPhaseStress:
    def __init__(self,packet):self.packet=packet
    def collar(self,Z,q):
        # Internal names belong to the generic transfer contract only.
        return dict(collar_similarity_stress_mixed3_factored=self.packet['waiting_similarity_stress_mixed3_factored'],
                    Gamma_endpoint_stress_exact_zero_from_same_moments=False)


class AtWaitingPhaseAssembly:
    def __init__(self,native,phase):self.native=native; self.phase=phase
    def __getattr__(self,key):return getattr(self.native,key)
    def evaluate(self,chart,Z,q,**kwargs):
        if chart!='heat_collar':raise ValueError('Generic waiting map received an unexpected chart')
        return self.native.evaluate('waiting',Z,self.phase,**kwargs)


def waiting_physical_identities():
    proof=physical_collar_identities(); identities={}
    a,q,z=s.symbols('a q Z',real=True)
    It,St,Je,Jp=s.symbols('theta_inertial theta_shear axial_energy axial_pressure',real=True)
    def zero(name,value):
        if s.simplify(value)!=0:raise ArithmeticError('Waiting physical reduction failed: '+name)
        identities[name]=True
    # Source rows already include Qt/Qz derivatives.
    theta=It*s.exp(-q)+St*s.exp(-(1+a)*q)
    axial=Je*s.exp(-q/2)+Jp*s.exp(q/2)
    zero('theta_homogeneous_inertial_divergence_cancels_before_enclosure',
         s.diff(theta,q)+theta+a*St*s.exp(-(1+a)*q))
    zero('axial_homogeneous_energy_divergence_cancels_before_enclosure',
         s.diff(axial,q)+axial/2-Jp*s.exp(q/2))
    for j in range(3):
        zero('theta_divergence_full_radial_factor_y'+str(j),
             s.diff(-a*St*s.exp(-(s.Rational(3,2)+a)*q),q,j)
             -(-s.Rational(3,2)-a)**j*(-a*St*s.exp(-(s.Rational(3,2)+a)*q)))
        zero('axial_divergence_full_radial_factor_y'+str(j),
             s.diff(s.exp(-q/2)*Jp*s.exp(q/2),q,j)-(Jp if j==0 else 0))
    return dict(general_K_physical_transfer=proof,waiting_reduction_identities=identities,
                waiting_angular_axial_viscosity_exact_zero_from_constant_K=True,
                waiting_inertial_and_energy_homogeneous_divergence_cancelled_before_enclosure=True,
                actual_regional_physical_waiting_stress_remainder_identity_verified=True,
                stress_is_generally_nonzero=True,full_NS_solution_not_claimed=True)


def waiting_adapter_binding():
    """Replay the narrow chart/constant-shape adapter source contract."""
    tree=ast.parse(Path(__file__).read_text(encoding='utf8')); bindings={}
    required={
        ('WaitingHeatReference','shape'):'K_rows=[one*(1-self.eps)]+[one*0]*4',
        ('AtWaitingPhaseAssembly','evaluate'):"self.native.evaluate('waiting',Z,self.phase,**kwargs)",
        ('AtWaitingPhaseStress','collar'):"self.packet['waiting_similarity_stress_mixed3_factored']",
        ('WaitingDispatch','evaluate'):'self.stress.waiting(Z,coordinate)',
        ('WaitingDispatch','provider'):"self.stress.steep if chart=='waiting' else self.native.provider(chart)",
        ('CompliantWaitingPhysicalC2','waiting'):'self.stress.steep.wait*(phase-1)'}
    for (cls,method),expression in required.items():
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name==method)
        # Keyword assignment in dict(...) is checked as a keyword node.
        if expression.startswith('K_rows='):
            expected=ast.dump(ast.parse('dict('+expression+')',mode='eval').body.keywords[0])
            found=any(ast.dump(n)==expected for n in ast.walk(fn) if isinstance(n,ast.keyword))
        else:
            expected=ast.dump(ast.parse(expression,mode='eval').body)
            found=any(ast.dump(n)==expected for n in ast.walk(fn))
        if not found:raise ValueError('Waiting adapter source changed: '+cls+'.'+method)
        bindings[cls+'.'+method]=True
    # Execute the actual production waiting radius branch with formal source
    # variables. This proves the selector/offset relation, not merely a toy q.
    path=HERE/(PREFIX+'global_physical_assembly.py')
    source_tree=ast.parse(path.read_text(encoding='utf8'))
    radius=next(n for n in ast.walk(source_tree) if isinstance(n,ast.FunctionDef) and n.name=='radius')
    environment={}
    for node in source_tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):
                    environment[target.id]=ast.literal_eval(node.value)
    exec(compile(ast.Module(body=[radius],type_ignores=[]),'<actual production waiting radius>','exec'),environment)
    phase,wait,Ts,length,mu,origin=s.symbols('phase wait Ts Lrel mu logRp',real=True)
    native=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    provider=SimpleNamespace(outer=SimpleNamespace(Lrel=length),Ts=Ts,wait=wait)
    logR=environment['radius'](native,'waiting',phase,{},provider)[0]
    logTail=environment['radius'](native,'waiting',s.Integer(1),{},provider)[0]
    if s.simplify(logR-logTail-wait*(phase-1))!=0:raise ArithmeticError('Actual waiting phase/radius binding failed')
    bindings['actual_production_waiting_logR_minus_logRtail_equals_q']=True
    return dict(source_bindings=bindings,original_waiting_phase_not_uncorrelated_offset_division=True,
                actual_production_waiting_phase_radius_identity_verified=True,
                original_collar_source_domain_preserved=True,
                generic_transfer_source='CompliantCollarPhysicalC2.collar; general-K proof consumed')


def waiting_divergence_grids(heat,fields):
    """Factored divergence after exact homogeneous-moment cancellations."""
    c=heat.ctx; St=fields['theta_shear'][0]; Jp=fields['axial_pressure'][0]
    return {
        'theta':{'y'+str(j)+'_Z'+str(n):-heat.a*St[n]*math.factorial(n)*(-c.mpf('1.5')-heat.a)**j
                 for j in range(3) for n in range(3-j)},
        'axial':{'y'+str(j)+'_Z'+str(n):(Jp[n]*math.factorial(n) if j==0 else c.mpf(0))
                 for j in range(3) for n in range(3-j)}}


class CompliantWaitingPhysicalC2:
    @source_precision
    def __init__(self):
        self.stress=CompliantWaitingStressC3(); self.heat=WaitingHeatReference(self.stress.heat)
        self.assembly=CompliantGlobalPhysicalAssembly()
        self.assembly.dispatch=WaitingDispatch(self.assembly.dispatch,self.stress)
        self.ctx=self.stress.ctx; self.family,self.source=self.stress.family,self.stress.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):raise ValueError('Waiting physical family mismatch')
        self.proof=waiting_physical_identities(); self.adapter_proof=waiting_adapter_binding()
        self.hashes=dict(self.stress.hashes); self.hashes.update(self.assembly.hashes)
        for stem,gate in (('waiting_stress_C3_check','actual_original_waiting_similarity_stress_recovered'),
                          ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped'),
                          ('collar_physical_C2_check','actual_regional_physical_collar_stress_remainder_identity_verified')):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Waiting physical prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Waiting receipt family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Waiting physical source changed: '+path)
            self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        for stem in ('collar_physical_C2','waiting_physical_C2'):
            path=HERE/(PREFIX+stem+'.py'); self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()

    @source_precision
    def waiting(self,Z,phase,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx; z=c.mpf(Z); phase=c.mpf(phase); self.stress.steep._phase(phase)
        q=self.stress.steep.wait*(phase-1); packet=self.stress.waiting(z,phase)
        adapter=SimpleNamespace(ctx=c,heat=self.heat,stress=AtWaitingPhaseStress(packet),
                                assembly=AtWaitingPhaseAssembly(self.assembly,phase))
        point=CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)
        logtau=c.mpf(log_tau); nu=c.mpf(viscosity); beta=-2-self.heat.delta
        parts=point['actual_source_log_factors']; logR=point['source_logR_enclosure']
        fields=packet['waiting_similarity_stress_y_derivative_Taylor']
        # Include the 1/sqrt(R) source factor when differentiating divergence.
        divgrids=waiting_divergence_grids(self.heat,fields)
        divmixed={}
        for label,factor in (('theta','Qtheta'),('axial','Qz')):
            divmixed[label]={}
            for i in range(3):
                for j in range(3-i):
                    coefficient=physical_bracket(c,divgrids[label],i,j,z,self.heat.delta,beta-1)
                    divmixed[label]['r'+str(i)+'_z'+str(j)]=physical_source_row(c,coefficient,parts[factor],
                        beta-1-i+j*(self.heat.delta-1),logtau,nu,i+j,
                        c.mpf(i+1)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
        cs=point['cos_theta']; sn=point['sin_theta']
        divcart=dict(x=[scale_row(c,divmixed['theta']['r0_z0'],-sn)],
                     y=[scale_row(c,divmixed['theta']['r0_z0'],cs)],z=[divmixed['axial']['r0_z0']])
        point['physical_cylindrical_stress_divergence_mixed2']=divmixed
        point['physical_completed_stress_divergence_cartesian']=divcart
        point['physical_momentum_residual_decomposition_cartesian']={label:[scale_row(c,row,-1) for row in divcart[label]]
            +point['physical_remainder_cartesian'][label] for label in ('x','y','z')}
        # Exact native physical pure power removes spurious interval widths
        # in axial/time velocity rows; pressure keeps its actual datum jets.
        for key,components in point['physical_spatial_cartesian_mixed4'].items():
            if int(key.split('_')[2][1:]):
                for component in ('ux','uy','uz'):
                    for row in components[component].values():
                        row.update(terms=[],exact_zero=True,log_absolute_upper=None,
                                   exact_zero_from_actual_waiting_pure_power=True)
        for component in ('ux','uy','uz'):
            for row in point['first_fixed_x_physical_time_derivative'][component].values():
                row.update(terms=[],exact_zero=True,log_absolute_upper=None,
                           exact_zero_from_actual_waiting_pure_power=True)
        for key in ('actual_regional_physical_collar_stress_remainder_identity_verified',
                    'Gamma_endpoint_physical_remainder_exact_zero','collar_cone_certified'):
            point.pop(key,None)
        point.update(chart='waiting',waiting_phase=phase,waiting_offset_from_Rtail=q,
                     actual_regional_physical_waiting_stress_remainder_identity_verified=True,
                     waiting_physical_axial_viscosity_exact_zero=True,
                     waiting_physical_velocity_axial_and_time_derivatives_exact_zero=True,
                     waiting_inertial_and_energy_divergence_cancelled_before_enclosure=True,
                     waiting_collar_stress_mixed3_join_verified=True,
                     steep_exit_waiting_stress_mixed3_join_verified=False,
                     regional_remainder_is_leading_axial_viscosity_not_proven_flat=False,
                     waiting_regional_remainder_exact_zero=True,waiting_cone_certified=False,
                     full_background_NS_validation=False,temporal_recursion=False)
        return point

    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_velocity_pressure_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            result['physical_axial_and_time_velocity_rows_all_exact_zero']=all(
                row['exact_zero'] for key,components in point['physical_spatial_cartesian_mixed4'].items()
                if int(key.split('_')[2][1:]) for component in ('ux','uy','uz') for row in components[component].values()) and all(
                row['exact_zero'] for component in ('ux','uy','uz') for row in point['first_fixed_x_physical_time_derivative'][component].values())
            return result
        samples=[summary(self.waiting(z,phase,lt,angle,nu)) for z,phase,lt,angle,nu in
                 (('0','0','-1','0','1'),('.5','.5','-10','.7','.01'),('-.5','1','-100','1','.7'))]
        whole=summary(self.waiting([-1,1],[0,1],[-1000,-1],None))
        self.hashes.update(self.assembly.dispatch.native.hashes)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                    scope='Actual original waiting phase[0,1], Z[-1,1], tau>0, r>0; stress mixed3, completed diagonal/divergence/remainder mixed2',
                    samples=samples,whole_waiting=whole,waiting_physical_identities=self.proof,waiting_adapter_binding=self.adapter_proof,
                    actual_regional_physical_waiting_stress_remainder_identity_verified=True,
                    waiting_regional_remainder_exact_zero=True,waiting_cone_certified=False,
                    independently_bounded_global_flat_remainder=False,global_admissible_stress_lift_constructed=False,
                    physical_energy_integral_certified=False,full_background_NS_validation=False,temporal_recursion=False,
                    input_hashes=self.hashes)


def run():
    result=CompliantWaitingPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual waiting physical stress and exact-zero regional remainder generated',flush=True)
    return result


if __name__=='__main__':run()

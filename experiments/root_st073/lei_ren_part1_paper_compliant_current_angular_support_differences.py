"""Current local angular stress3/error2 differences with factored full moments.

Same nonzero boundary histories; only the local beta input is removed.
The enormous common KR is an exact log source, never a materialized cap.
This is local interface flatness, not the global temporal flat remainder.
"""
import ast
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace,FunctionType

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_angular_support_interfaces import (
    CurrentAngularSupportInterfaces,EDGES,OPEN,HERE,PREFIX,sha,pack,encode,endpoints,
    statement,function,binding,BASE)
from lei_ren_part1_paper_compliant_current_pulse_support_interfaces import (
    beta_tail_bound,symmetric_jet,ordinary_grid)
from lei_ren_part1_paper_compliant_angular_stress_C3 import angular_transport_identities
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_stress_rows,collar_moment_stress_identities
from lei_ren_part1_paper_compliant_collar_physical_C2 import (
    collar_velocity_bracket,collar_velocity_operator_coefficients,stress_source_log_parts,physical_source_row)
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_angular_support_differences.json'
RECEIPT=PREFIX+'current_angular_support_differences_check.json'
GATES=('current_angular_support_local_stress3_error2_flat_difference_bounds_available',)
DISTANCES=('.01','.000001','0')


def exact_current_KR_source(physical):
    physical.assert_graph();c=physical.ctx;heat=physical.history.heat;steep=physical.history.steep
    bindings={}
    for target,value in (('self.thetaR','c.exp(-100*self.bp-self.bp*self.outer.Lrel)/2'),
        ('self.thetaS','self.thetaR*c.exp(-self.bp-self.rate/2)'),
        ('self.thetaQ',"self.thetaS*c.exp(-c.mpf('1.5')*self.Ts)"),
        ('self.thetaT',"self.thetaQ*c.exp(-c.mpf('1.5')+self.k/2)")):
        class_assignment('current_steep_waiting_source','CurrentSteepWaitingC4','__init__',target,value)
        bindings[target]=value
    class_assignment('current_heat_source','CurrentCollarGammaC4','__init__','self.theta_base',
        'self.steep.thetaT*c.exp(-self.bh*self.steep.wait-self.steep.logone)')
    if heat.steep is not steep or steep.outer is not physical.history.outer:
        raise ValueError('Current angular/steep/heat normalization graph differs')
    if endpoints(steep.Ts)!=endpoints(c.mpf(endpoints(steep.future.params.Ts))):
        raise ValueError('Current KR uses a different steep length')
    if endpoints(steep.logone)!=endpoints(c.mpf(endpoints(steep.future.angular.waiting_logone))):
        raise ValueError('Current KR uses a different waiting normalization')
    text="B=dict(common,inverse_mu=-13/(2*mu),flatten_power=-bp*(100+heat.steep.outer.Lrel),steep=-c.mpf('1.5')*heat.steep.Ts,waiting_and_current=-bh*(heat.steep.wait+offset),finite=-15-mu/2-a/2-c.ln(2))"
    bindings['original_physical_B_log_parts']=statement('compliant_collar_physical_C2','stress_source_log_parts',text)
    mu,a,T,W,l=s.symbols('mu a Ts wait log_one_minus_epsilon',real=True)
    bp=s.Rational(1,2)+mu;bh=s.Rational(1,2)+a;r=1-mu;k=1-a
    rel=-bp-r/2-s.Rational(3,2)*T-s.Rational(3,2)+k/2-bh*W-l
    source=-rel-bh*(W+T+2);reduced=1+mu/2-3*a/2+k*T+l
    if s.expand(source-reduced)!=0:raise ArithmeticError('Current exact KR logarithm reduction failed')
    q=y=s.symbols('source_s',real=True)-W-T-2
    logEv0,logthetaR=s.symbols('exact_logEv0 exact_logthetaR',real=True)
    BminusEv0=logthetaR-bp-r/2-s.Rational(3,2)*T-s.Rational(3,2)+k/2-bh*(W+q)-l
    if s.expand(BminusEv0+reduced+(a-mu)*(q+W+T+2)-(logthetaR-bp*(q+W+T+2)))!=0:
        raise ArithmeticError('Current native angular amplitude and tail stress normalization differ')
    logKR=1+steep.mu/2-3*heat.a/2+heat.k*steep.Ts+steep.logone
    if not all(mp.isfinite(v) for v in endpoints(logKR)):raise ValueError('Finite current exact KR log required')
    return dict(actual_current_normalization_AST=bindings,
        exact_KR_definition='thetaR/(theta_base*exp(-bh*qR)), qR=-wait-Ts-2',
        exact_KR_log='1+mu/2-3*a/2+(1-a)*Ts+log(1-epsilon)',logKR=logKR,
        exact_theta_cancellation='B*KR*exp((a-mu)*s)=Ev0*thetaR*exp(-(.5+mu)*s)',
        no_division_by_theta_base_cap_or_materialization_of_KR=True,passed=True)


def linear_stress_program():
    """Original stress operator, excluding its unrelated K1/K0 diagnostic."""
    fn=copy.deepcopy(function('compliant_collar_stress_C3','collar_stress_rows'))
    ret=fn.body[-1]
    if not isinstance(ret,ast.Return) or not isinstance(ret.value,ast.Call):raise ValueError('Original stress return changed')
    removed=[v for v in ret.value.keywords if v.arg=='shear_strength_margin']
    if len(removed)!=1 or ast.dump(removed[0].value)!=ast.dump(ast.parse('K[1]/K[0]*(-2)+heat.delta',mode='eval').body):
        raise ValueError('Original non-difference shear diagnostic changed')
    ret.value.keywords=[v for v in ret.value.keywords if v.arg!='shear_strength_margin']
    env=dict(collar_stress_rows.__globals__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original-linear-current-angular-stress>','exec'),env)
    return env[fn.name]


def local_difference_theorem():
    mu,a,y=s.symbols('mu a s',real=True);h=s.Function('local_beta_input')(y)
    Kref=s.exp((a-mu)*y);actual=Kref*(1+h)
    checks={}
    expressions={'linear_K':actual-Kref-Kref*h,
        'quadratic_K':actual**2-Kref**2-Kref**2*(2*h+h*h)}
    for label,value in expressions.items():
        if s.expand(value)!=0:raise ArithmeticError('Current angular signed input difference failed')
        checks[label]=True
    dK=[s.Symbol('deltaK'+str(j)) for j in range(5)]
    dQ=[s.Symbol('deltaQ'+str(j)) for j in range(5)]
    A=[s.Symbol('same_boundary_deltaA')];E=[s.Symbol('same_boundary_deltaE')];P=[s.Symbol('same_boundary_deltaP')]
    for j in range(4):A.append(dK[j]-(1-a)*A[j]);E.append(2*a*E[j]-dQ[j]);P.append((1+2*a)*P[j]-dQ[j]/2)
    subs={A[0]:0,E[0]:0,P[0]:0,**dict.fromkeys(dK+dQ,0)}
    for label,rows in (('A',A),('E',E),('P',P)):
        for j,value in enumerate(rows):
            if s.expand(value.subs(subs))!=0:raise ArithmeticError('Nonzero angular difference endpoint')
            checks[label+'_endpoint_order'+str(j)]=True
    return dict(identities=checks,
        exact_difference_equations=dict(A='deltaA_y=deltaK-k*deltaA',E='deltaE_y=delta*deltaE-deltaQ',
            P='deltaP_y=prate*deltaP-deltaQ/2'),
        normalized_differences='deltaK/KR, deltaA/KR; deltaQ/KR^2, deltaE/KR^2, deltaP/KR^2',
        both_orientations_integral_bound='abs(deltaMoment) <= h*exp(abs(rate)*h)*sup(abs(input difference))',
        exact_source_boundary_differences_zero_actual_histories_not_zeroed=True,
        linear_stress_and_axial_viscosity_operators_preserve_difference_flatness=True,
        full_nonzero_histories_pressure_datum_and_complete_future_retained=True,
        local_flat_limit='For every fixed current positive source factor, beta derivatives through4 and their vanishing-interval integrals tend to zero as h->0; this is not temporal tau-flatness',passed=True)


def original_operator_difference_theorem():
    """Actual AST superposition with arbitrary nonzero reference histories."""
    y,z,a,S=s.symbols('source_y Z a exact_inverse_Rtail',real=True)
    ctx=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),exp=s.exp)
    heat=SimpleNamespace(ctx=ctx,a=a,delta=2*a,k=1-a,prate=1+2*a,S=S)
    program=linear_stress_program();env=dict(program.__globals__)
    env['IntervalTaylor']=SimpleNamespace(variable=lambda c,Z,n:Z)
    env['axial_derivative']=lambda value:s.diff(value,z)
    program=FunctionType(program.__code__,env)
    refs={name:s.Function('nonzero_reference_'+name)(y,z) for name in ('A','E','P','K')}
    diffs={name:s.Function('current_difference_'+name)(y,z) for name in refs}
    def evaluate(values):
        rows={name:[s.diff(value,y,j) for j in range(5)] for name,value in values.items()}
        return program(heat,dict(K_rows=rows['K']),dict(angular_defect_rows=rows['A'],
            energy_defect_rows=rows['E'],pressure_defect_rows=rows['P'],K_defect_rows=rows['K']),z,y)
    actual=evaluate({name:refs[name]+diffs[name] for name in refs});reference=evaluate(refs);difference=evaluate(diffs)
    proofs={}
    for label in ('theta','axial','theta_inertial','theta_shear'):
        for j in range(4):
            zero=s.expand(actual[label][j]-reference[label][j]-difference[label][j])
            if zero!=0:raise ArithmeticError('Original stress difference dropped nonzero history: '+label)
            for n in range(4-j):proofs[label+'_y%d_Z%d'%(j,n)]=s.diff(zero,z,n)==0
    remainder={}
    for i in range(3):
        for j in range(3-i):
            operator=collar_velocity_operator_coefficients(i,j+2)
            if any(sum(key)>4 for key in operator):raise ValueError('Remainder difference exceeds mixed4')
            aa=sum(value*s.Symbol('actual_row%d_%d'%key) for key,value in operator.items())
            rr=sum(value*s.Symbol('reference_row%d_%d'%key) for key,value in operator.items())
            dd=sum(value*(s.Symbol('actual_row%d_%d'%key)-s.Symbol('reference_row%d_%d'%key)) for key,value in operator.items())
            if s.expand(aa-rr-dd)!=0:raise ArithmeticError('Original axial viscosity difference is not linear')
            remainder['r%d_z%d'%(i,j)]=True
    return dict(actual_AST_stress_mixed3_superposition_identities=proofs,
        original_axial_viscosity_mixed2_superposition_identities=remainder,
        arbitrary_nonzero_reference_A_E_P_K_retained=True,passed=True)


class CurrentAngularSupportDifferences:
    @source_precision
    def __init__(self,support=None,require_checked=True):
        self.support=support if support is not None else CurrentAngularSupportInterfaces()
        if not self.support.acceptance_loaded:raise ValueError('Checked current eight internal source traces required')
        self.physical=self.support.physical;self.outer=self.support.outer;self.heat=self.physical.history.heat;self.ctx=self.outer.ctx
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.KR=exact_current_KR_source(self.physical);self.stress_operator=linear_stress_program()
        self.transport=angular_transport_identities();self.units=collar_moment_stress_identities();self.proof=local_difference_theorem()
        self.operator_proof=original_operator_difference_theorem()
        self.hashes=dict(self.support.hashes)
        for stem in ('angular_stress_C3','collar_stress_C3','collar_physical_C2','collar_Gamma_C4','pulse_physical_bounds'):
            name=PREFIX+stem+'.py'
            if name in self.hashes and self.hashes[name]!=sha(name):raise ValueError('Current angular difference operator source changed')
            self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or any(receipt[k] for k in OPEN):raise ValueError('Angular difference receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def interface(self,edge,h,Z=(-1,1),log_tau='-1',viscosity='1'):
        if edge not in EDGES:raise ValueError('Unknown current angular support edge')
        self.physical.assert_graph();c=self.ctx;h=c.mpf(h);Z=c.mpf(Z);nu=c.mpf(viscosity);logtau=c.mpf(log_tau)
        if endpoints(h)[0]<0 or endpoints(h)[1]>endpoints(c.mpf('.01'))[1]:raise ValueError('Local h in[0,.01] required')
        if endpoints(nu)[0]<=0 or not all(mp.isfinite(v) for v in endpoints(nu)+endpoints(logtau)):raise ValueError('Finite positive nu and finite log(tau) required')
        controls=self.outer.data(Z)['coeff'];zero=controls[0]*0;ell=c.mpf('.15');pos=s.Rational(edge['exact_edge'])
        center=c.mpf(str(pos.p))/pos.q
        coordinate=c.mpf([endpoints(center-h)[0],endpoints(center+h)[1]])
        g=self.heat.a-self.outer.mu;factor=c.exp(g*coordinate)
        H=[zero]*5 if endpoints(h)[1]==0 else [symmetric_jet(controls[edge['row']]*(beta_tail_bound(c,j,2*h/ell)/(ell**(j+1)*self.outer.normalization))) for j in range(5)]
        Q=[H[j]*2+product_rows(H,H)[j] for j in range(5)]
        dK=[sum((H[j]*math.comb(n,j)*g**(n-j) for j in range(n+1)),zero)*factor for n in range(5)]
        dQ=[sum((Q[j]*math.comb(n,j)*(2*g)**(n-j) for j in range(n+1)),zero)*factor**2 for n in range(5)]
        A=[symmetric_jet(dK[0]*(h*c.exp(self.heat.k*h)))]
        E=[symmetric_jet(dQ[0]*(h*c.exp(self.heat.delta*h)))]
        P=[symmetric_jet(dQ[0]*(h*c.exp(self.heat.prate*h)/2))]
        for j in range(4):A.append(dK[j]-A[j]*self.heat.k);E.append(E[j]*self.heat.delta-dQ[j]);P.append(P[j]*self.heat.prate-dQ[j]/2)
        q=-self.heat.steep.wait-self.heat.steep.Ts-2+coordinate
        defects=dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=P,K_defect_rows=dK)
        stress=self.stress_operator(self.heat,dict(K_rows=dK),defects,Z,q)
        grids={name:ordinary_grid(stress[name],3) for name in ('theta','axial','theta_inertial','theta_shear')}
        error={'r%d_z%d'%(i,j):-collar_velocity_bracket(c,dK,i,j+2,Z,self.heat.delta) for i in range(3) for j in range(3-i)}
        # Original source log prefactors and physical operators; KR/KR^2
        # remain common source logs and are never substituted as field caps.
        logR,radius_source=BASE.radius(self.physical,'outer_angular',coordinate,{},self.outer)
        parts=stress_source_log_parts(self.physical,self.heat,q)
        parts={key:{**value,'current_exact_KR_log':self.KR['logKR']*(2 if key=='Qz' else 1)} for key,value in parts.items() if key in ('B','Qtheta','Qz')}
        beta=-2-self.heat.delta;physical_stress={}
        for label,f in (('theta','Qtheta'),('axial','Qz')):
            physical_stress[label]={}
            for i in range(4):
                for j in range(4-i):
                    bracket=physical_bracket(c,grids[label],i,j,Z,self.heat.delta,beta)
                    physical_stress[label]['r%d_z%d'%(i,j)]=physical_source_row(c,bracket,parts[f],beta-i+j*(self.heat.delta-1),logtau,nu,i+j,c.mpf(i)/2*(c.ln(2)-logR))
        physical_error={}
        for i in range(3):
            for j in range(3-i):
                physical_error['r%d_z%d'%(i,j)]=physical_source_row(c,error['r%d_z%d'%(i,j)],parts['B'],-1-self.heat.delta-i+(j+2)*(self.heat.delta-1),logtau,nu,i+j,c.mpf(i)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
        return dict(edge=edge,h=h,Z=Z,exact_source_edge=edge['exact_edge'],diagnostic_source_coordinate_enclosure=coordinate,
            local_beta_difference_rows=H,normalized_K_difference_rows=dK,normalized_quadratic_K_difference_rows=dQ,
            same_boundary_normalized_full_moment_difference_rows=dict(A=A,E=E,P=P),
            normalized_similarity_stress_difference_mixed3=grids,
            normalized_axial_viscosity_remainder_difference_mixed2=error,
            original_physical_stress_difference_mixed3_log_bounds=physical_stress,
            original_physical_axial_viscosity_difference_mixed2_log_bounds=physical_error,
            exact_common_source_factor=self.KR,radius_source=radius_source,
            current_nonzero_actual_histories_and_absolute_pressure_preserved=True,
            h_zero_means_difference_zero_not_actual_field_zero=True,
            temporal_flat_remainder_not_inferred=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_exact_KR_source=self.KR,
            original_angular_full_moment_transport_theorem=self.transport,original_paper_stress_units=self.units,
            current_angular_local_difference_source_theorem=self.proof,
            current_original_stress_and_viscosity_difference_operator_theorem=self.operator_proof,
            original_linear_stress_AST_sha256=sha(PREFIX+'collar_stress_C3.py'),
            original_physical_axial_viscosity_operator_sha256=sha(PREFIX+'collar_physical_C2.py'),
            only_shear_strength_ratio_diagnostic_omitted_from_difference_replay=True,
            current_source_identified_adjacent_interface_count=14,current_source_identified_internal_support_count=8,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentAngularSupportDifferences(require_checked=False)
    result=field.manifest();result['current_local_angular_support_difference_views']=[]
    for edge in EDGES:
        for h in DISTANCES:result['current_local_angular_support_difference_views'].append(field.interface(edge,h))
        print('Current angular local stress/error generated: '+edge['exact_edge'],flush=True)
    result['input_hashes']=field.hashes
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

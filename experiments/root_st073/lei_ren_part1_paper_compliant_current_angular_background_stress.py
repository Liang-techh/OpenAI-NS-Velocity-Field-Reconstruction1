"""Actual current angular background stress, stable full pressure and tensor.

Normalized full moments retain nonzero histories. Unit baselines cancel in
the original stress AST before KR/KR^2 factors are separated into source logs.
Only the original angular region is admitted, not a global cone or flatness.
"""
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace,FunctionType

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_interface_atlas import (
    CurrentInterfaceAtlas,HERE,PREFIX,OPEN,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_angular_support_differences import exact_current_KR_source,linear_stress_program
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import statement
from lei_ren_part1_paper_compliant_current_selected_energy_source import binding
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_current_full_exterior_stress import current_full_history_transfer
from lei_ren_part1_paper_compliant_angular_stress_C3 import angular_future_changes,angular_transport_identities,SourceAST
from lei_ren_part1_paper_compliant_future_energy_high_jets import copy_jet
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_stress_rows,collar_moment_stress_identities
from lei_ren_part1_paper_compliant_angular_physical_C2 import angular_physical_identities
from lei_ren_part1_paper_compliant_collar_physical_C2 import physical_source_row,scale_row,collar_velocity_bracket
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_current_pulse_support_interfaces import ordinary_grid
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_angular_background_stress.json'
RECEIPT=PREFIX+'current_angular_background_stress_check.json'
GATES=('current_actual_angular_full_moment_stress_mixed3_recovered',
    'current_actual_angular_absolute_pressure_remaining_integral_identified',
    'current_actual_angular_completed_physical_tensor_remainder_decomposition_available')
VIEWS={'whole_angular':((-1,1),(-4,0),('-3','-1'),None,'1'),
    'power_boundary':((-1,1),'-4','-1',None,'1'),
    'entry_boundary':((-1,1),'0','-1',None,'1'),
    'first_support':((-1,1),('-3.16','-2.84'),'-1',None,'1'),
    'second_support':((-1,1),('-1.16','-.84'),'-1',None,'1'),
    'fresh_first':('.537','-2.97','-2.6','.41','.8'),
    'fresh_gap':('-.419','-2.1','-4.3','.7','.03')}


def normalized_full_moment_rows(heat,outer,native,data,future,post,v,c):
    """Directed copies of same native functions in one Taylor context."""
    one=IntervalTaylor.constant(c,1,5)
    F=[copy_jet(c,native['swirl_factor_one_plus_h_Taylor'])]+[
        copy_jet(c,row) for row in native['actual_angular_bump_y_derivatives'][1:]]
    g=heat.a-outer.mu;fac=c.exp(g*v)
    K=[sum((F[n]*math.comb(j,n)*g**(j-n) for n in range(j+1)),one*0)*fac for j in range(5)]
    Q=product_rows(K,K);A=[K[0]*copy_jet(c,native['angular_Taylor'])]
    E=[(copy_jet(c,data['post'])+copy_jet(c,future['energy'])+
        one*(decay_integral(c,2*outer.mu,-v)*c.exp(-2*outer.mu*v)))*c.exp(heat.delta*v)]
    P=[(copy_jet(c,post['normalized_full_post_pressure'])+copy_jet(c,future['pressure'])+
        one*(decay_integral(c,outer.prate,-v)*c.exp(-outer.prate*v)/2))*c.exp(heat.prate*v)]
    for j in range(4):
        A.append(K[j]-A[j]*heat.k);E.append(E[j]*heat.delta-Q[j]);P.append(P[j]*heat.prate-Q[j]/2)
    return dict(A=A,E=E,P=P,K=K)


def full_normalized_moment_AST_theorem():
    """Replay original full defect rows, retaining every nonzero datum."""
    v,z,a,mu,KR=s.symbols('s Z a mu KR',real=True,positive=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),exp=s.exp,
        expm1=lambda value:s.exp(value)-1)
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,delta=2*a,k=1-a,prate=1+2*a)
    outer=SimpleNamespace(mu=mu,prate=1+2*mu)
    F=s.Function('same_current_F')(v,z);X=s.Function('same_actual_X')(v,z)
    ER=s.Function('same_full_post_energy')(z);PR=s.Function('same_full_post_pressure')(z)
    JE=s.Function('same_remaining_beta_energy')(v,z);JP=s.Function('same_remaining_beta_pressure')(v,z)
    native=dict(swirl_factor_one_plus_h_Taylor=F,angular_Taylor=X,
        actual_angular_bump_y_derivatives=[F]+[s.diff(F,v,j) for j in range(1,5)])
    stub=SimpleNamespace(constant=lambda ctx,value,order:s.Rational(value))
    env=dict(normalized_full_moment_rows.__globals__);env.update(IntervalTaylor=stub,
        copy_jet=lambda ctx,value:value,decay_integral=lambda ctx,p,length:(1-s.exp(-p*length))/p)
    rows=FunctionType(normalized_full_moment_rows.__code__,env)(heat,outer,native,
        dict(post=ER),dict(energy=JE,pressure=JP),dict(normalized_full_post_pressure=PR),v,c)
    asts=SourceAST();original_env=dict(IntervalTaylor=stub,math=math,
        decay_integral=env['decay_integral'])
    asts.replay('angular_stress_C3','angular_defect_rows',original_env)
    shape=dict(K_rows=[KR*row for row in rows['K']],
        K_defect_rows=[KR*rows['K'][0]-1]+[KR*row for row in rows['K'][1:]])
    square=product_rows(shape['K_rows'],shape['K_rows'])
    shape['K_squared_defect_rows']=[square[0]-1]+square[1:]
    original=original_env['angular_defect_rows'](heat,
        dict(KR=KR,energy_defect_rows=[KR**2*ER-1/heat.delta],
            pressure_defect_rows=[KR**2*PR-1/(2*heat.prate)]),v,shape,X,dict(energy=JE,pressure=JP))
    checks={}
    for name,key,base in (('A','angular_defect_rows',1/heat.k),
            ('E','energy_defect_rows',1/heat.delta),('P','pressure_defect_rows',1/(2*heat.prate))):
        for j in range(5):
            zero=s.cancel(s.expand(original[key][j]+(base if j==0 else 0)-
                rows[name][j]*(KR if name=='A' else KR**2)))
            if zero!=0:raise ArithmeticError('Current full normalized moment differs from original AST: '+name+str(j))
            for n in range(5-j):checks[name+'_y%d_Z%d'%(j,n)]=s.diff(zero,z,n)==0
    binding('compliant_power_angular_C4','angular','energy',
        "((data['post']+futureE)*c.exp(2*self.mu*s)+decay_integral(c,2*self.mu,-s))/(F[0]*F[0]*2)")
    return dict(original_full_defect_AST_normalized_mixed4_identities=checks,
        original_signed_quadratic_beta_future_and_nonzero_post_moments_retained=True,
        original_native_energy_normalization_bound=True,input_hashes=asts.hashes,passed=True)


def full_normalized_stress_AST_theorem():
    """Retain then cancel exact unit baselines before normalized full rows."""
    y,z,a,S,KR=s.symbols('y Z a inverse_radius KR',real=True,positive=True)
    ctx=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp)
    heat=SimpleNamespace(ctx=ctx,a=a,delta=2*a,k=1-a,prate=1+2*a,S=S)
    program=linear_stress_program();env=dict(program.__globals__)
    env['IntervalTaylor']=SimpleNamespace(variable=lambda c,Z,n:Z)
    env['axial_derivative']=lambda v:s.diff(v,z)
    program=FunctionType(program.__code__,env)
    full={name:s.Function('actual_normalized_'+name)(y,z) for name in ('A','E','P','K')}
    def evaluate(values,defects):
        rows={name:[s.diff(value,y,j) for j in range(5)] for name,value in values.items()}
        return program(heat,dict(K_rows=rows['K']),dict(angular_defect_rows=defects['A'],
            energy_defect_rows=defects['E'],pressure_defect_rows=defects['P'],K_defect_rows=defects['K']),z,y)
    actual={name:value*(KR if name in ('A','K') else KR**2) for name,value in full.items()}
    actual_rows={name:[s.diff(value,y,j) for j in range(5)] for name,value in actual.items()}
    defect_rows={name:list(rows) for name,rows in actual_rows.items()}
    for name,baseline in (('A',1/heat.k),('E',1/heat.delta),('P',1/(2*heat.prate)),('K',1)):
        defect_rows[name][0]-=baseline
    original=evaluate(actual,defect_rows)
    normalized=evaluate(full,{name:[s.diff(value,y,j) for j in range(5)] for name,value in full.items()})
    checks={}
    for name in ('theta','axial','theta_inertial','theta_shear'):
        for j in range(4):
            zero=s.cancel(s.expand(original[name][j]-normalized[name][j]*(KR**2 if name=='axial' else KR)))
            if zero!=0:raise ArithmeticError('Original unit baseline or KR homogeneity differs')
            for n in range(4-j):checks[name+'_y%d_Z%d'%(j,n)]=s.diff(zero,z,n)==0
    return dict(actual_original_stress_mixed3_baseline_cancellation_and_KR_homogeneity=checks,
        original_defects='A-1/k, E-1/delta, P-1/(2prate), K-1 retained before exact AST cancellation',
        resulting_source_units='Qtheta*KR for theta; Qz*KR^2 for axial',
        actual_nonzero_full_A_E_P_K_functions_retained=True,passed=True)


def stable_post_pressure_source_proof(physical):
    history=physical.history;physical.assert_graph();transfer=current_full_history_transfer(history)
    if not transfer['passed'] or not all(transfer['current_exact_function_links'].values()):
        raise ValueError('Current pressure function link missing')
    binding('compliant_current_steep_waiting_source','data','PS',
        "PR+self.infull['pressure']*(self.outer.flatten.Ev2*self.thetaR**2)")
    binding('compliant_current_steep_waiting_source','data','PQ',
        'PS+decay_integral(c,3,self.Ts)*(self.outer.flatten.Ev2*self.thetaS**2/2)')
    binding('compliant_current_steep_waiting_source','data','PT',
        "PQ+self.outfull['pressure']*(self.outer.flatten.Ev2*self.thetaQ**2)")
    binding('compliant_steep_waiting_C4','waiting','pressure',
        "data['PT']+decay_integral(c,1+self.delta,t)*(self.outer.flatten.Ev2*self.thetaT**2/2)")
    binding('compliant_current_angular_terminal_closure','terminal_constants','pressure_infinity',
        "terminal['Ptail']+tail0['remaining_pressure_in_Rtail_units']*terminal['pressure_scale']")
    rS,rQ,rT,rB,Ii,I3,Io,Iw,Ih,C0,PR=s.symbols('rS rQ rT rB Iin I3 Iout Iwait Iheat C0 actual_PR',real=True)
    post=Ii+rS*I3/2+rQ*Io+rT*Iw/2+rB*Ih
    env={'PR':PR,"self.infull['pressure']":Ii,'self.outer.flatten.Ev2':C0,
        'self.thetaR':s.Integer(1),'self.thetaS':s.sqrt(rS),'self.thetaQ':s.sqrt(rQ),
        'self.thetaT':s.sqrt(rT),'decay_integral(c, 3, self.Ts)':I3,
        'decay_integral(c, 1 + self.delta, t)':Iw}
    env['PS']=assignment('compliant_current_steep_waiting_source','data','PS',env)
    env['PQ']=assignment('compliant_current_steep_waiting_source','data','PQ',env)
    env["self.outfull['pressure']"]=Io
    PT=assignment('compliant_current_steep_waiting_source','data','PT',env)
    env["data['PT']"]=PT
    tail=assignment('compliant_steep_waiting_C4','waiting','pressure',env)
    env.update({"terminal['Ptail']":tail,"terminal['pressure_scale']":C0*rB,
        "tail0['remaining_pressure_in_Rtail_units']":Ih})
    Cp=assignment('compliant_current_angular_terminal_closure','terminal_constants','pressure_infinity',env)
    zero=s.expand(Cp-PR-C0*post)
    if zero!=0:raise ArithmeticError('Current whole pressure primitive does not telescope')
    checks={}
    z=s.Symbol('Z')
    for n in range(6):
        if s.diff(zero,z,n)!=0:raise ArithmeticError('Absolute pressure terminal axial identity failed')
        checks['actual_PR_from_current_Cp_zero_Z'+str(n)]=True
    return dict(current_full_history_and_original_pressure_function_transfer=transfer,
        current_native_pressure_increment_AST_and_Cp_definition_replayed=True,
        exact_AST_telescoping_identity='Cp=PR+Ev2*thetaR^2*Jpost; checked current Cp=0 implies PR=-Ev2*thetaR^2*Jpost',
        exact_normalized_post_pressure='Iin+rS*I3(Ts)/2+rQ*Iout+rT*Iprate(wait)/2+rB*remaining_heat_pressure(0)',
        exact_ratio_logs=dict(rS='-2bp-rate',rQ='-2bp-rate-3Ts',rT='-2bp-rate-3Ts-3+k',
            rB='-2bp-rate-3Ts-3+k-2bh*wait-2logone'),
        identities=checks,current_Cp_zero_is_original_analytic_P0_identity=True,
        heat_prate='1+delta',angular_beta_pressure_rate='1+2mu',
        no_pressure_gauge_shift_or_division_by_Ev2_theta_base_caps=True,passed=True)


def angular_source_log_parts(physical,offset):
    """Original collar parts * KR factors, cancelled analytically first."""
    c=physical.ctx;heat=physical.history.heat;mu=heat.mu;a=heat.a;bp=c.mpf('.5')+mu
    common=dict(logPstar=physical.logP,actual_inlet_logU=physical.flatten.logEv2_parts['inlet_log']/2)
    B=dict(common,inverse_mu=-13/(2*mu),flatten_power=-bp*(100+heat.outer.Lrel),
        angular_offset=-heat.bh*offset,finite=-13-c.ln(2))
    Qt=dict(common,radial_origin=physical.logRp/2,flatten_power=-mu*(100+heat.outer.Lrel),
        angular_offset=-a*offset,finite=-13-c.mpf('1.5')*c.ln(2))
    Qz={**{'Qtheta_'+k:v for k,v in Qt.items()},**{'B_'+k:v for k,v in B.items()}}
    diag={k:2*v for k,v in Qt.items()};diag['log2']=c.ln(2)
    pressure={k:2*v for k,v in B.items()}
    return dict(B=B,Qtheta=Qt,Qz=Qz,completed_diagonal=diag,absolute_pressure=pressure)


def source_log_cancellation_theorem():
    lp,lu,rp,mu,a,L,T,W,l,v=s.symbols('logP logU logRp mu a Lrel Ts wait logone s',real=True)
    bp=s.Rational(1,2)+mu;bh=s.Rational(1,2)+a;q=-W-T-2+v
    KR=1+mu/2-3*a/2+(1-a)*T+l
    oldB=lp+lu-l-13/(2*mu)-bp*(100+L)-s.Rational(3,2)*T-bh*(W+q)-15-mu/2-a/2-s.log(2)
    oldQt=lp+lu-l+rp/2-mu*(100+L)-T-a*(W+q)-14-mu/2-a/2-s.Rational(3,2)*s.log(2)
    newB=lp+lu-13/(2*mu)-bp*(100+L)-bh*v-13-s.log(2)
    newQt=lp+lu+rp/2-mu*(100+L)-a*v-13-s.Rational(3,2)*s.log(2)
    pairs={'B_KR':(oldB+KR,newB),'Qtheta_KR':(oldQt+KR,newQt),
        'Qz_KR_squared':(oldQt+oldB+2*KR,newQt+newB),
        'diagonal_KR_squared':(2*oldQt+s.log(2)+2*KR,2*newQt+s.log(2))}
    checks={name:s.expand(left-right)==0 for name,(left,right) in pairs.items()}
    if not all(checks.values()):raise ArithmeticError('Current exact angular positive source log cancellation failed')
    return dict(identities=checks,KR_not_materialized=True,current_Ts_wait_and_logone_cancel_before_enclosure=True,passed=True)


class CurrentAngularBackgroundStress:
    @source_precision
    def __init__(self,atlas=None,require_checked=True):
        self.atlas=atlas if atlas is not None else CurrentInterfaceAtlas()
        if not self.atlas.acceptance_loaded:raise ValueError('Checked current interface atlas required')
        self.physical=self.atlas.physical;self.history=self.physical.history
        self.outer=self.history.outer;self.heat=self.history.heat;self.ctx=self.physical.ctx
        self.pressure=self.history.selected.pressure
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.KR=exact_current_KR_source(self.physical);self.units=collar_moment_stress_identities()
        self.transport=angular_transport_identities();self.moment_proof=full_normalized_moment_AST_theorem()
        self.baselines=full_normalized_stress_AST_theorem()
        self.pressure_proof=stable_post_pressure_source_proof(self.physical)
        self.physical_proof=angular_physical_identities();self.log_proof=source_log_cancellation_theorem()
        self.hashes=dict(self.atlas.hashes)
        for stem in ('angular_stress_C3','angular_physical_C2','collar_stress_C3','collar_physical_C2','current_full_exterior_stress'):
            name=PREFIX+stem+'.py'
            if name in self.hashes and self.hashes[name]!=sha(name):raise ValueError('Actual stress operator changed')
            self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current actual angular tensor admission exceeds regional scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.atlas.assert_graph();self.physical.assert_graph()
        if not (self.physical is self.atlas.physical and self.history is self.physical.history and
                self.outer is self.history.outer and self.heat is self.history.heat and
                self.pressure is self.history.selected.pressure and self.ctx is self.physical.ctx):
            raise ValueError('Actual angular stress must retain the same current physical/history/pressure owners')
        if not self.pressure.acceptance_loaded:raise ValueError('Original checked current absolute pressure required')
        steep=self.history.steep
        if steep.infull is not steep.kernels(1,'in') or steep.outfull is not steep.kernels(1,'out'):
            raise ValueError('Current pressure transition integrals do not retain their native source owners')

    @source_precision
    def post_pressure(self,Z):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);key=tuple(endpoints(Z))
        if not all(mp.isfinite(x) for x in key) or key[0]<-1 or key[1]>1:
            raise ValueError('Finite source Z[-1,1] required')
        if key not in self.cache:
            steep=self.history.steep;heat=self.heat
            if steep.infull is not steep.kernels(1,'in') or steep.outfull is not steep.kernels(1,'out'):
                raise ValueError('Current pressure transition integrals do not retain their native source owners')
            rS=-2*steep.bp-steep.rate;rQ=rS-3*steep.Ts;rT=rQ-3+steep.k
            rB=rT-2*steep.bh*steep.wait-2*steep.logone
            logs=dict(rS=rS,rQ=rQ,rT=rT,rB=rB);factors={k:c.exp(v) for k,v in logs.items()}
            one=IntervalTaylor.constant(c,1,5);tails=heat.collar_tails(Z,0)
            terms=dict(entry=one*steep.infull['pressure'],
                power=one*(factors['rS']*decay_integral(c,3,steep.Ts)/2),
                exit=one*(factors['rQ']*steep.outfull['pressure']),
                waiting=one*(factors['rT']*decay_integral(c,heat.prate,steep.wait)/2),
                complete_heat=copy_jet(c,tails['remaining_pressure_in_Rtail_units'])*factors['rB'])
            self.cache[key]=dict(normalized_full_post_pressure=sum(terms.values(),one*0),
                exact_positive_ratio_logs=logs,normalized_source_contributions=terms,
                original_full_Gamma_pressure_future=tails['remaining_pressure_in_Rtail_units'])
        return self.cache[key]

    @source_precision
    def angular(self,Z,offset,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);v=c.mpf(offset);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(v)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(v)[0]<-4 or endpoints(v)[1]>0:
            raise ValueError('Actual angular region Z[-1,1], s[-4,0] required')
        if endpoints(nu)[0]<=0 or not all(mp.isfinite(x) for x in endpoints(nu)+endpoints(lt)):
            raise ValueError('Finite compact log(tau) and finite constant nu>0 required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):
            raise ValueError('Finite theta or all-angle enclosure required')
        native=self.outer.angular(Z,v);data=self.outer.data(Z);future=angular_future_changes(self.outer,data,v)
        post=self.post_pressure(Z)
        full_rows=normalized_full_moment_rows(self.heat,self.outer,native,data,future,post,v,c)
        A,E,P,K=(full_rows[name] for name in ('A','E','P','K'));one=IntervalTaylor.constant(c,1,5)
        q=-self.history.steep.wait-self.history.steep.Ts-2+v
        logR,radius_source=self.physical.radius('outer_angular',v,{},self.outer)
        # Enclose exact 1/R from the same S*exp(-q) source, combining its
        # cap log first. This is a bound, never an exact inverse-radius value.
        caplog=c.ln(c.mpf(endpoints(self.heat.Scap)[1]))-q
        if endpoints(caplog+logR)[0]<0:raise ValueError('Angular inverse-radius cap does not enclose the same exact radius')
        cap=c.exp(c.mpf(endpoints(caplog)[1]));inverseR=c.mpf([0,endpoints(cap)[1]])
        local_heat=copy.copy(self.heat);local_heat.S=inverseR;local_heat.ctx=c
        full=dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=P,K_defect_rows=K)
        stress=collar_stress_rows(local_heat,dict(K_rows=K),full,Z,c.mpf(0))
        grids={name:ordinary_grid(stress[name],3) for name in ('theta','axial','theta_inertial','theta_shear')}
        parts=angular_source_log_parts(self.physical,v);beta=-2-self.heat.delta
        ps={};div={}
        for name,factor in (('theta','Qtheta'),('axial','Qz')):
            ps[name]={};div[name]={}
            for i in range(4):
                for j in range(4-i):
                    bracket=physical_bracket(c,grids[name],i,j,Z,self.heat.delta,beta)
                    ps[name]['r%d_z%d'%(i,j)]=physical_source_row(c,bracket,parts[factor],beta-i+j*(self.heat.delta-1),lt,nu,i+j,c.mpf(i)/2*(c.ln(2)-logR))
            for i in range(3):
                for j in range(3-i):
                    bracket=physical_bracket(c,grids[name],i+1,j,Z,self.heat.delta,beta)
                    for m in range(i+1):bracket+=physical_bracket(c,grids[name],i-m,j,Z,self.heat.delta,beta)*((2 if name=='theta' else 1)*math.comb(i,m)*(-1)**m*math.factorial(m)/c.mpf(2)**(m+1))
                    div[name]['r%d_z%d'%(i,j)]=physical_source_row(c,bracket,parts[factor],beta-1-i+j*(self.heat.delta-1),lt,nu,i+j,c.mpf(i+1)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
        diagonal={};error={}
        for i in range(3):
            for j in range(3-i):
                bracket=physical_bracket(c,grids['axial'],i,j+1,Z,self.heat.delta,beta)
                if i:bracket+=physical_bracket(c,grids['axial'],i-1,j+1,Z,self.heat.delta,beta)*(c.mpf(i)/2)
                diagonal['r%d_z%d'%(i,j)]=physical_source_row(c,bracket,parts['completed_diagonal'],beta+self.heat.delta-i+j*(self.heat.delta-1),lt,nu,i+j,c.mpf(i)/2*(c.ln(2)-logR))
                error['r%d_z%d'%(i,j)]=physical_source_row(c,-collar_velocity_bracket(c,K,i,j+2,Z,self.heat.delta),parts['B'],-1-self.heat.delta-i+(j+2)*(self.heat.delta-1),lt,nu,i+j,c.mpf(i)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
        cs=sn=c.mpf([-1,1]) if theta is None else None
        if theta is not None:cs=c.cos(c.mpf(theta));sn=c.sin(c.mpf(theta))
        tt=ps['theta']['r0_z0'];tz=ps['axial']['r0_z0'];dd=diagonal['r0_z0'];er=error['r0_z0']
        zero=physical_source_row(c,c.mpf(0),parts['B'],-1-self.heat.delta,lt,nu,0)
        tensor=dict(xx=[scale_row(c,tt,-2*cs*sn),scale_row(c,dd,sn*sn)],
            xy=[scale_row(c,tt,cs*cs-sn*sn),scale_row(c,dd,-cs*sn)],xz=[scale_row(c,tz,cs)],
            yy=[scale_row(c,tt,2*cs*sn),scale_row(c,dd,cs*cs)],yz=[scale_row(c,tz,sn)],zz=[zero])
        divcart=dict(x=[scale_row(c,div['theta']['r0_z0'],-sn)],y=[scale_row(c,div['theta']['r0_z0'],cs)],z=[div['axial']['r0_z0']])
        ecart=dict(x=[scale_row(c,er,-sn)],y=[scale_row(c,er,cs)],z=[zero])
        residual={name:[scale_row(c,row,-1) for row in divcart[name]]+ecart[name] for name in ('x','y','z')}
        pressure_rows=[-sum((P[n]*math.comb(j,n)*(-self.heat.prate)**(j-n) for n in range(j+1)),one*0) for j in range(5)]
        return dict(Z=Z,angular_offset=v,log_tau=lt,viscosity=nu,
            current_actual_normalized_full_moment_rows=dict(A=A,E=E,P=P,K=K),
            current_normalized_post_pressure_source=post,current_signed_remaining_angular_quadratic_kernels=future,
            current_actual_normalized_stress_mixed3=grids,
            physical_cylindrical_stress_mixed3=ps,physical_stress_divergence_mixed2=div,
            completed_theta_theta_stress_mixed2=diagonal,physical_axial_viscosity_remainder_mixed2=error,
            completed_background_tensor_cartesian_components=tensor,
            physical_completed_stress_divergence_cartesian=divcart,physical_remainder_cartesian=ecart,
            physical_momentum_residual_decomposition_cartesian=residual,
            stable_actual_absolute_pressure_mixed4_factored=ordinary_grid(pressure_rows,4),
            original_forward_absolute_pressure_Taylor=native['pressure_over_Pstar_squared_Taylor'],
            original_current_energy_Taylor=native['energy_Taylor'],normalized_energy_native_consistency=E[0]-K[0]*K[0]*copy_jet(c,native['energy_Taylor'])*2,
            exact_current_positive_log_source_parts=parts,exact_source_logR=logR,radius_source=radius_source,
            inverse_radius_enclosure=inverseR,inverse_radius_exact_source='1/R=S_exact*exp(-q), q=-wait-Ts-2+s',
            inverse_radius_cap_only_not_exact_field_value=True,
            exact_pure_swirl_meridional_velocity_radial_divergence_and_remainder_zero=True,
            actual_pressure_is_same_P0_Pin_forward_function_by_current_Cp_zero=True,
            actual_full_stress_not_local_difference=True,original_full_nonzero_A_E_P_and_signed_C5_controls_retained=True,
            regional_decomposition_is_not_global_temporal_flatness=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_exact_KR_source=self.KR,
            original_full_moment_stress_units=self.units,current_angular_transport_source_theorem=self.transport,
            current_full_moment_normalization_AST_theorem=self.moment_proof,
            current_full_moment_baseline_and_normalization_AST_theorem=self.baselines,
            current_original_absolute_pressure_future_source_theorem=self.pressure_proof,
            original_general_K_physical_tensor_remainder_theorem=self.physical_proof,
            exact_current_angular_source_log_cancellation=self.log_proof,
            domain='current original angular s[-4,0], source Z[-1,1], R>0, finite compact log(tau), constant nu>0; Z endpoints are infinity limits',
            scoped_actual_tensor='Trtheta,T rz and Ttheta_theta=r*partial_z(Trz), symmetric; others zero',
            scoped_actual_remainder='Er=Ez=0, Etheta=-nu*partial_zz(utheta)',
            tensor_cone_and_global_flatness_not_inferred=True,current_complete_22_interface_atlas_retained=True,
            other_current_actual_chart_tensors_and_tensor_joins_remain_open=True,input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentAngularBackgroundStress(require_checked=False)
    result=field.manifest();result['current_actual_angular_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_angular_tensor_views'][name]=field.angular(*args)
        print('Current actual angular full stress/tensor: '+name,flush=True)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

"""Independent functional-interface fixture and current source receipt check."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_end_flatten_join import report,DOMAIN,FALSE_FLAGS
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import pulse_coefficients
from lei_ren_part1_paper_compliant_flatten_stress_C3 import flatten_shape
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_stress_rows,shifted_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_heat_pressure_C4 import pressure_y_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_physical_interface_fixture():
    """Direct original full moment formulas, with unequal reference rates."""
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=95
        mu=mp.mpf('.2');a=mp.mpf('.06');delta=2*a;r=1-mu
        bp=mp.mpf('.5')+mu;bh=mp.mpf('.5')+a;p=1+2*a;pp=1+2*mu
        R0=mp.mpf('2.8');B0=mp.mpf('.45');C0=mp.mpf('.7')
        Xp=mp.mpf('.3');H0=mp.mpf('.08');D=mp.mpf('.07');q0=mp.mpf('-3')
        Xv=1/r+(Xp-1/r)*H0
        C=lambda z:1/(1+z*z)
        Cz=lambda z:-2*z/(1+z*z)**2
        future=lambda z:mp.mpf('.8')+mp.mpf('.1')*z*z+mp.mpf('.02')*z**3
        future_z=lambda z:mp.mpf('.2')*z+mp.mpf('.06')*z*z
        pressure0=lambda z:-mp.mpf('.35')*C(z)**2-mp.mpf('.04')*z*z
        pressure0_z=lambda z:-mp.mpf('.7')*C(z)*Cz(z)-mp.mpf('.08')*z
        datum=lambda z:mp.mpf('.12')+mp.mpf('.02')*z*z
        def direct(v,z):
            R=R0*mp.exp(v);B=B0*mp.exp(-bp*v);Ut=B*C(z)
            X=1/r+(Xp-1/r)*H0*mp.exp(-r*v)
            e=future(z)*mp.exp(2*mu*v)-mp.expm1(2*mu*v)/(4*mu)
            ez=future_z(z)*mp.exp(2*mu*v)
            kernel=-mp.expm1(-pp*v)/(2*pp)
            pressure=B0**2*(pressure0(z)+C(z)**2*kernel)
            pressure_z=B0**2*(pressure0_z(z)+2*C(z)*Cz(z)*kernel)
            Mt=mp.sqrt(2)*R**mp.mpf('1.5')*Ut*X
            Mt_z=mp.sqrt(2)*R**mp.mpf('1.5')*B*Cz(z)*X
            Me=R*Ut*Ut*e;Me_z=R*B*B*(2*C(z)*Cz(z)*e+C(z)**2*ez)
            L=1-delta*z*z;d=1-z*z;b=(1-delta)/2
            It=((1-a)*Mt-b*z*Mt_z-R*mp.sqrt(2*R)*Ut)/(2*L*R)
            Iz=(2*delta*z*Me-d*Me_z+R*(2*p*z*pressure-d*pressure_z))/(L*mp.sqrt(2*R))
            shear=(-2*bp*Ut-Ut)/mp.sqrt(2*R)
            return dict(Utheta=Ut,Mtheta=Mt,Mztheta=Me,pressure=pressure,
                Mp=pressure-datum(z),theta=It+shear,axial=Iz)
        checks=0;misses=[];sides={}
        def compare(label,rows,zv,total,side):
            nonlocal checks
            for j in range(total+1):
                for n in range(total-j+1):
                    actual=rows[j][n]*mp.factorial(n)
                    expected=mp.diff(lambda z:mp.diff(lambda v:direct(v,z)[label],mp.mpf(0),j),zv,n)
                    lo,hi=endpoints(actual);miss=max(lo-expected,expected-hi,mp.mpf(0))
                    if miss>mp.mpf('1e-55'):
                        raise ArithmeticError('Independent interface fixture failed: '+str((side,label,zv,j,n,miss)))
                    checks+=1;misses.append(str(miss));sides[side]=sides.get(side,0)+1
        for zv in map(mp.mpf,('0','.5','-.7')):
            iv=lambda value:c.mpf(value)
            jet=lambda fn:IntervalTaylor(c,mp.taylor(fn,zv,5))
            one=IntervalTaylor.constant(c,1,5);zero=one*0
            z=IntervalTaylor.variable(c,iv(zv),5);CJ=jet(C);fu=jet(future);PJ=jet(pressure0)
            m=iv(mu);av=iv(a);dv=2*av;kv=1-av;pv=1+dv;rv=1-m
            heat=SimpleNamespace(ctx=c,a=av,mu=m,k=kv,delta=dv,prate=pv,S=iv(mp.exp(q0)/R0))
            shape=flatten_shape(heat,dict(Kright=iv(C0*mp.exp(100*(a-mu))/2)),iv(-100),
                dict(Z=iv(zv),F_Taylor=one,sigma_y_derivatives=[iv(0)]*5,actual_original_flatten_right_endpoint=False))
            K=shape['K_rows'];Q=shape['K_squared_defect_rows']
            fullA=K[0]*iv(Xv);fullE=CJ*CJ*fu*iv(2*C0**2);fullP=-PJ*iv(C0**2)
            A=[fullA-one/kv];E=[fullE-one/dv];Pr=[fullP-one/(2*pv)]
            for j in range(4):
                A.append(shape['K_defect_rows'][j]-A[j]*kv)
                E.append(E[j]*dv-Q[j]);Pr.append(Pr[j]*pv-Q[j]/2)
            defects=dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
                K_defect_rows=shape['K_defect_rows'],K_squared_defect_rows=Q)
            flat=collar_stress_rows(heat,shape,defects,iv(zv),iv(q0))
            e=[fu];P=[PJ];X=[one*iv(Xv)]
            for j in range(4):
                e.append(e[j]*(2*m)-(one/2 if j==0 else zero))
                P.append(P[j]*(1+2*m)+(CJ*CJ/2 if j==0 else zero))
                X.append((one if j==0 else zero)-X[j]*rv)
            pulse=pulse_coefficients(dv,m,z,CJ,iv(Xp),[zero]*5,[zero]*5,[zero]*5,e,[zero]*5,P)
            for label in ('theta','axial'):
                rows=[zero for _ in range(4)]
                for part in pulse[label].values():
                    rp,bpow,dp,hp=part['mode']
                    factor=iv(R0)**rp*iv(B0)**bpow*iv(D)**dp*iv(H0)**hp/c.sqrt(2)
                    for j in range(4):rows[j]+=part['full_derivative_rows'][j]*factor
                compare(label,rows,zv,3,'pulse')
                factor=mp.sqrt(R0/2)*(B0/C0 if label=='theta' else B0**2/C0**2)
                compare(label,[v*iv(factor) for v in flat[label]],zv,3,'flatten')
            rawP=pressure_y_rows(K,fullP,pv,iv(B0**2/C0**2*mp.exp(p*q0)),iv(q0))
            compare('pressure',rawP,zv,4,'flatten')
            compare('pressure',[v*iv(B0**2) for v in shifted_rows(P,-iv(pp),4)],zv,4,'pulse')
            datumJ=jet(datum)
            for side,rows in (
                ('pulse',[v*iv(B0**2) for v in shifted_rows(P,-iv(pp),4)]),
                ('flatten',rawP)):
                compare('Mp',[rows[0]-datumJ]+rows[1:],zv,4,side)
            compare('Utheta',[v*iv(B0/C0) for v in shifted_rows(K,-iv(bh),4)],zv,4,'flatten')
            compare('Mtheta',[v*iv(mp.sqrt(2)*R0**mp.mpf('1.5')*B0/C0)
                for v in shifted_rows([fullA]+A[1:],kv,4)],zv,4,'flatten')
            compare('Mztheta',[v*iv(R0*B0**2/(2*C0**2))
                for v in shifted_rows([fullE]+E[1:],-dv,4)],zv,4,'flatten')
        return dict(passed=True,independent_original_physical_mixed_derivative_checks=checks,
            checks_by_side=sides,tolerance='1e-55',misses=misses,
            unequal_reference_rates=True,nonzero_signed_angular_memory=True,
            analytic_datum_retained_in_Mp=True,actual_project_cone_or_NS_validation=False)


def run():
    name=PREFIX+'pulse_end_flatten_join.json';record=json.loads((HERE/name).read_bytes())
    for source,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Interface source changed: '+source)
    fresh=encode(pack(report()))
    if fresh!=record:raise ValueError('Current functional join report differs')
    if record['domain']!=DOMAIN or not record['pulse_end_flatten_full_moment_stress_pressure_functional_join_verified']:
        raise ValueError('Actual complete similarity interface not established')
    for flag in FALSE_FLAGS+('source_caps_used_as_defining_field_values','interval_overlap_used_as_join_proof'):
        if record[flag]:raise ValueError('Interface scope overclaimed: '+flag)
    fixture=independent_physical_interface_fixture()
    hashes=dict(record['input_hashes']);hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
        implicit_source_sha256=record['implicit_source_sha256'],domain=DOMAIN,
        pulse_end_flatten_full_moment_stress_pressure_functional_join_verified=True,
        current_actual_endpoint_source_recomputed=True,
        exact_functional_join_identities=len(record['functional_join']['identities']),
        actual_local_inlet_and_datum_source_identities=len(record['local_inlet_and_datum_source_binding']['identities']),
        actual_endpoint_source_identities=len(record['source_endpoint_binding']['identities']),
        independent_physical_interface_fixture=fixture,input_hashes=hashes,
        source_caps_used_as_defining_field_values=False,interval_overlap_used_as_join_proof=False,
        **{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS actual pulse-end / flatten full moment mixed4, stress3 and absolute pressure4 functional join; physical/cone pending',flush=True)
    return result


if __name__=='__main__':run()

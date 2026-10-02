"""Independent implicit physical coordinate derivatives and log ledger gates."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    PulsePhysicalBounds,physical_operators,physical_bracket,UZ,UT,UR,P,ZSYM,DSYM,BSYM)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_pulse_physical_bounds.json'


def operator_identities():
    z,delta,beta=ZSYM,DSYM,BSYM; y=s.symbols('y',real=True); f=s.Function('G')(y,z)
    op=physical_operators(); checks={}
    def apply(row):return sum(v*s.diff(f,y,k,z,n) for (k,n),v in row.items())
    for (a,b),row in op.items():
        h=apply(row)
        if b>0:
            previous=apply(op[(a,b-1)])
            expected=((beta+(b-1)*(delta-1))*z*previous+(1-z*z)*s.diff(previous,z)-2*z*s.diff(previous,y))/(1-delta*z*z)
            if s.simplify(h-expected)!=0:raise ArithmeticError('Physical axial operator recurrence failed')
        elif a>0:
            previous=apply(op[(a-1,0)])
            if s.simplify(h-(s.diff(previous,y)-s.Rational(a-1,2)*previous))!=0:
                raise ArithmeticError('Physical radial lowering factor failed')
        checks['physical_r'+str(a)+'_z'+str(b)+'_operator']=True
    # The crucial cancellation of the radial weight in the first axial
    # chain rule is checked separately, rather than accepting a missing -a.
    a=s.symbols('a',integer=True,nonnegative=True); H=s.Function('H')(y,z)
    scaled=s.exp(-a*y/2)*H
    full=(beta-a)*z*scaled+(1-z*z)*s.diff(scaled,z)-2*z*s.diff(scaled,y)
    reduced=s.exp(-a*y/2)*(beta*z*H+(1-z*z)*s.diff(H,z)-2*z*s.diff(H,y))
    if s.simplify(full-reduced)!=0:raise ArithmeticError('Radial weight/axial lambda cancellation failed')
    checks['radial_weight_cancels_axial_minus_a']=True
    return checks


def independent_physical_coordinates():
    with mp.workdps(90):
        c=MPIntervalContext(); c.dps=110; tol=mp.mpf('1e-65')
        delta=mp.mpf('.03'); tau=mp.mpf('.7'); y0=mp.mpf('.2'); checks={}
        G=lambda y,z:mp.exp(mp.mpf('.17')*y)*(1+mp.mpf('.2')*z+mp.mpf('.3')*z*z)+mp.exp(-mp.mpf('.07')*y)*mp.sin(z)
        for Z in (mp.mpf('-.4'),mp.mpf(0),mp.mpf('.6')):
            lam=mp.sqrt(tau/(1-Z*Z)); R=mp.exp(y0); rr=lam*mp.sqrt(2*R); zz=lam**(1-delta)*Z
            grid={}
            for k in range(5):
                for n in range(5-k):
                    v=mp.diff(G,(y0,Z),(k,n)); grid['y'+str(k)+'_Z'+str(n)]=c.mpf([v-tol,v+tol])
            for label,beta in ((UZ,-1-delta),(UT,-1-delta),(UR,mp.mpf(-1)),(P,-2-2*delta)):
                def physical(r,z):
                    # An independent positive implicit root in PHYSICAL
                    # z, then the actual paper coordinates. No operator
                    # or production coordinate derivative is used here.
                    root=mp.findroot(lambda l:l*l-l**(2*delta)*z*z-tau,
                        (mp.sqrt(tau+z*z),mp.sqrt(tau+z*z)+mp.mpf('.1')),
                        tol=mp.eps*16,verify=True)
                    coordinate=z/root**(1-delta)
                    radius=r*r/(2*root*root)
                    return root**beta*G(mp.log(radius),coordinate)
                for (a,b) in physical_operators():
                    target=mp.diff(physical,(rr,zz),(a,b))
                    bracket=physical_bracket(c,grid,a,b,c.mpf(Z),c.mpf(delta),c.mpf(beta))
                    exponent=c.mpf(beta)-a+b*(c.mpf(delta)-1)
                    # Half-integer radial power must be evaluated by sqrt;
                    # the production map itself stores the exact log scale.
                    prefactor=c.exp(exponent*c.ln(c.mpf(lam)))*c.sqrt(c.mpf(2)/c.mpf(R))**a
                    value=bracket*prefactor; lo,hi=endpoints(value)
                    if not lo<=target<=hi:raise ArithmeticError('Independent physical implicit derivative failed: '+str((label,Z,a,b)))
                    checks[label+'_Z'+str(Z)+'_r'+str(a)+'_z'+str(b)]=True
        return dict(checks=checks,implicit_physical_coordinate_derivative_count=len(checks),
            finite_parameter_fixture_only=True,actual_Md40_source_admission=False,passed=True)


def independent_scale_fixture():
    with mp.workdps(85):
        c=MPIntervalContext(); c.dps=100
        field=PulsePhysicalBounds.__new__(PulsePhysicalBounds); field.ctx=c
        field.mu=c.mpf('.04'); field.delta=c.mpf('.03'); field.logP=c.mpf('2.3')
        field.logRp_parts=dict(logCstar=c.mpf('3.7'),logPstar=c.mpf('23'),finite_outer_offset=c.mpf('4.1'))
        logRp=sum(field.logRp_parts.values(),c.mpf(0)); count=0
        for label in (UZ,UT,UR,P):
            for a,b in physical_operators():
                for chi,offset in (('.02','0'),('12','0'),('13','-4')):
                    scale=field.component_scale(label,a,b,chi,offset,'-10')
                    actual=sum((v for key,v in scale.items() if key.endswith('_term')),c.mpf(0))
                    beta=scale['original_profile_lambda_exponent']; gamma=scale['physical_lambda_exponent']
                    t=c.mpf(chi)/field.mu+c.mpf(offset); bp=c.mpf('.5')+field.mu
                    if label==UR:expected=field.logP+(1-a)*logRp/2+(a-1)*c.ln(2)/2-(field.mu+c.mpf(a)/2)*t-gamma*5
                    elif label==P:expected=2*field.logP-a*logRp/2+a*c.ln(2)/2-c.mpf(a)*t/2-gamma*5
                    else:expected=field.logP-a*logRp/2+a*c.ln(2)/2-(bp+c.mpf(a)/2)*t-gamma*5
                    if max(endpoints(actual)[0],endpoints(expected)[0])>min(endpoints(actual)[1],endpoints(expected)[1]):
                        raise ArithmeticError('Physical radial/log scale mismatch')
                    count+=1
        return dict(independent_scale_equivalence_count=count,passed=True,actual_source_admission=False)


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Physical ledger source changed: '+name)
    c=MPIntervalContext(); c.dps=260; read=lambda v:read_interval(c,v)
    count=0
    for chart,components in r['all_chart_physical_derivative_log_bounds'].items():
        for label,rows in components.items():
            if len(rows)!=15:raise ValueError('Physical mixed derivative coverage missing')
            for row in rows.values():
                if any(not mp.isfinite(v) for v in endpoints(read(row['factored_bracket_absolute_upper']))):raise ArithmeticError('Bracket norm not finite')
                for logtau,sector in row['sectors'].items():
                    if endpoints(read(sector['physical_lambda_exponent']))[1]>=0:raise ArithmeticError('Invalid lambda lower-bound direction')
                    if sector['derivative_exactly_zero']:
                        if sector['physical_derivative_log_upper'] is not None:raise ArithmeticError('Exact zero bound replaced')
                    else:
                        if any(not mp.isfinite(v) for v in endpoints(read(sector['physical_derivative_log_upper']))):raise ArithmeticError('Actual physical log bound not finite')
                        parts=[read(v) for key,v in sector.items() if key.endswith('_term')]+[read(sector['factored_bracket_log_upper'])]
                        total=sum(parts,c.mpf(0)); admitted=read(sector['physical_derivative_log_upper'])
                        if max(endpoints(total)[0],endpoints(admitted)[0])>min(endpoints(total)[1],endpoints(admitted)[1]):raise ArithmeticError('Segmented log bound lost')
                    count+=1
    for key in ('full_pulse_C4_installed','full_outer_C4_certified','full_cartesian_vector_derivatives_certified','physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
        if r[key]:raise ValueError('Physical derivative scope overclaimed: '+key)
    proof=operator_identities(); fixture=independent_physical_coordinates(); scales=independent_scale_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        symbolic_physical_chain_rule_checks=proof,independent_physical_coordinate_fixture=fixture,independent_physical_scale_fixture=scales,
        actual_chart_derivative_sector_log_bounds_checked=count,
        all_cylindrical_r_z_derivatives_total_order_le4_mapped=True,
        complete_pulse_physical_spatial_supremum_log_ledger_available=True,all_passed=True,
        full_pulse_C4_installed=False,full_outer_C4_certified=False,full_cartesian_vector_derivatives_certified=False,
        physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Physical pulse map: exact chain rules, independent implicit coordinate derivatives, radial/time scales and all-chart log bounds PASS',flush=True)
    return out


if __name__=='__main__':run()

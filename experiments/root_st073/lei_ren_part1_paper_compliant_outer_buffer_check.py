"""Independent cumulative-ODE and interface checks for the O.3 buffer."""
# Recomputed for the distinct compliant pressure/moment family.
# Formula origin: lei_ren_part1_paper_shared_outer_buffer_check.py; old .01 receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_outer_buffer import SharedOuterBuffer,decay_integral
from lei_ren_part1_paper_compliant_outer_initial_check import overlaps
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def symbolic():
    t,mu=s.symbols('t mu',positive=True);u0,m0,h0,k0,e0,p0=s.symbols('u0 m0 h0 k0 e0 p0')
    J,K1,K2,K3=[s.Function(n)(t) for n in ('J','K1','K2','K3')];sig=s.Function('sig')(t)
    f=s.exp(-t/2-mu*J);d=s.exp(-t);d3=s.exp(-3*t/2);u=u0*f
    moments=[m0*d,h0*d3+u0*d3*K1,k0*d3,e0*d-u0**2*d*K2/2,p0+u0**2*K3/2]
    rules={s.diff(J,t):sig,s.diff(K1,t):s.exp(t-mu*J),
           s.diff(K2,t):s.exp(-2*mu*J),s.diff(K3,t):s.exp(-t-2*mu*J)}
    def test(u,values,rules):
        m,h,k,e,p=values;rhs=[-m,u-3*h/2,-3*k/2,-u**2/2-e,u**2/2]
        residuals=[s.simplify((s.diff(v,t)-r).subs(rules)) for v,r in zip(values,rhs)]
        if any(v!=0 for v in residuals):raise ArithmeticError('Outer buffer moment ODE failed: '+str(residuals))
        return [str(v) for v in residuals]
    transition=test(u,moments,rules)
    f=s.exp((-s.Rational(1,2)-mu)*t);u=u0*f
    Ktheta=(f-d3)/(1-mu);D2=(1-s.exp(-2*mu*t))/(2*mu);Dp=(1-s.exp(-(1+2*mu)*t))/(1+2*mu)
    power=test(u,[m0*d,h0*d3+u0*Ktheta,k0*d3,e0*d-u0**2*d*D2/2,p0+u0**2*Dp/2],{})
    return dict(slope_mu_five_ODE_residuals=transition,pure_power_five_ODE_residuals=power)


def run():
    result=symbolic();field=SharedOuterBuffer();c=field.ctx
    name=PREFIX+'compliant_outer_buffer.json';raw=json.loads((HERE/name).read_bytes())
    for source,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Buffer dependency changed: '+source)
    keys=('Utheta_over_Pstar','Utheta_y_over_Pstar','Uz','Uz_y','Mz_over_R',
          'Mtheta_over_sqrt2_R_3half_Pstar','Mtheta_z_over_sqrt2_R_3half_Pstar',
          'Mztheta_over_R_Pstar_squared','Mp_over_Pstar_squared','P_over_Pstar_squared')
    interfaces=[]
    with mp.workdps(210):
        for Z in ('-1','0','.5','1'):
            for label,left,right in [('d',field.initial.axial(Z,buffer_offset=11),field.slope_mu(Z,0)),
                                     ('w',field.slope_mu(Z,1),field.power(Z,0))]:
                checks={k:all(overlaps(a,b) for a,b in zip(left[k],right[k])) for k in keys}
                checks['Ur']=overlaps(left['Ur_over_sqrt_R_over_2'],right['Ur_over_sqrt_R_over_2'])
                if not all(checks.values()):raise ArithmeticError('O3 interface mismatch '+label+' '+Z+str(checks))
                interfaces.append(dict(Z=Z,interface=label,checks=checks))
        for phase in ('0','.5','1'):
            packet=field.power([-1,1],phase,cells=64)
            if not all(len(packet[k])==2 for k in keys):raise ArithmeticError('Whole-axis C1 buffer API failed')
        inlet=field.power('.5',1)
        if not inlet['pulse_inlet'] or endpoints(inlet['Uz'][0])!=(mp.mpf(0),mp.mpf(0)):
            raise ArithmeticError('Pulse inlet velocity/interface mismatch')
        if endpoints(inlet['Mz_over_R'][0])[1]<=0:raise ArithmeticError('Axial history reset in buffer')
        for k,t in (('.2','.0001'),('.5','2'),('1e-100','1e20')):
            enclosure=decay_integral(c,c.mpf(k),c.mpf(t))
            kk=mp.mpf(k);tt=mp.mpf(t);actual=-mp.expm1(-kk*tt)/kk
            if not endpoints(enclosure)[0]<=actual<=endpoints(enclosure)[1]:
                raise ArithmeticError('Cancellation-free decay integral enclosure failed')
    result.update(total_independent_symbolic_identities=10,interface_checks=interfaces,
        value_and_first_axial_derivative_interfaces_checked=True,
        whole_axis_C1_buffer_API_checked=True,retained_axial_history_at_pulse_inlet_checked=True,
        cancellation_free_decay_integral_examples_checked=True,
        actual_five_defect_family_sha256=field.initial.family,
        axial_pulse_and_outer_moment_repairs_installed=False,heat_exterior_matched=False,
        whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=dict(raw['input_hashes']))
    for source in (Path(__file__).name,name):
        result['input_hashes'][source]=hashlib.sha256((HERE/source).read_bytes()).hexdigest()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('O.3 buffer: 10 independent moment identities, 8 interfaces, whole-axis C1 and pulse inlet history PASS',flush=True)
    return result


if __name__=='__main__':run()

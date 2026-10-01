"""Independent Gamma heat/PDE, moments and retained-deficit checks."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_shared_exact_heat_component import SharedExactHeatComponent
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def symbolic():
    a,xi,v=s.symbols('a xi v',positive=True)
    w=s.exp(-v)*v**a;f=w*(1+xi*v)**(-a)
    ode=xi**2*s.diff(f,xi,2)+(1+2*(1+a)*xi)*s.diff(f,xi)+a*(1+a)*f
    boundary=a*s.diff(s.exp(-v)*v**(a+1)*(1+xi*v)**(-a-1),v)
    gamma=s.simplify(ode-boundary)
    r,H,Hx,Hxx=s.symbols('r H Hx Hxx',positive=True);b=a+s.Rational(1,2)
    ut=-4*Hx/r**2
    ur=(-2*b*H-2*xi*Hx)/r
    urr=(2*b*(2*b+1)*H+(8*b+6)*xi*Hx+4*xi**2*Hxx)/r**2
    heat=s.simplify(ut-urr-ur/r+H/r**2+4*(xi**2*Hxx+(1+(2+2*a)*xi)*Hx+a*(1+a)*H)/r**2)
    Rt,c,t,K=s.symbols('Rt c t K',positive=True);R=Rt*s.exp(t);U=c*R**(-b)*K
    moments=[s.simplify(s.sqrt(2*R)*U*R/(s.sqrt(2)*c*Rt**(1-a))-s.exp((1-a)*t)*K),
             s.simplify(U**2/2/(c*c*Rt**(-1-2*a))-s.exp(-(1+2*a)*t)*K**2/2),
             s.simplify(R*U**2/(c*c*Rt**(-2*a))-s.exp(-2*a*t)*K**2)]
    K0,SS,CC,Dhat,sig=s.symbols('K0 S C Dhat sig')
    actual=K0-a*SS*sig*CC*Dhat
    square=s.simplify((K0*K0-actual*actual)/(a*SS)-sig*CC*Dhat*(2*K0-a*SS*sig*CC*Dhat))
    # Differentiate the SAME source Dhat, not independent velocity fits.
    Dz=s.symbols('Dz')
    square_z=s.simplify(s.diff((K0*K0-actual*actual)/(a*SS),Dhat)*Dz-sig*CC*Dz*(2*K0-2*a*SS*sig*CC*Dhat))
    kk=s.symbols('k',positive=True);uraw=s.symbols('uraw');radius=s.symbols('radius',positive=True)
    angular_ftc=s.simplify(s.sqrt(2)*c*radius**(-a)+s.sqrt(2*radius)*(uraw-c*radius**(-b))-s.sqrt(2*radius)*uraw)
    eps,W=s.symbols('epsilon W')
    preheat=s.expand((1-eps*W)**2-(1-2*eps*W+eps**2*W**2))
    residuals=[gamma,heat]+moments+[square,square_z,angular_ftc,preheat]
    if any(v!=0 for v in residuals):raise ArithmeticError('Independent exact heat identity failed '+str(residuals))
    return dict(Gamma_ODE_boundary_identity=str(gamma),physical_angular_heat_PDE_factor_identity=str(heat),
        physical_three_moment_normalizations=[str(v) for v in moments],
        same_source_square_deficit_identity=str(square),same_source_square_derivative_identity=str(square_z),
        renormalized_angular_backward_FTC_identity=str(angular_ftc),separate_epsilon_correction_identity=str(preheat),
        total_symbolic_identities=len(residuals))


def run():
    result=symbolic();field=SharedExactHeatComponent();c=field.ctx
    name=PREFIX+'shared_exact_heat_component.json';raw=json.loads((HERE/name).read_bytes())
    for source,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Heat component dependency changed: '+source)
    if raw['paper_Section7_c_epsilon_gate'] or raw['heat_exterior_matched_to_incoming_five_moments']:
        raise ValueError('Current source threshold/global matching must remain unresolved')
    with mp.workdps(210):
        for Z in ('-1','0','.5','1'):
            h=field.deficit(Z,3);D=h['deficit_scaled_Taylor']
            if abs(mp.mpf(Z))==1:
                if endpoints(D[0])!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Heat deficit axis-endpoint value not zero')
            elif endpoints(D[0])[0]<=0:raise ArithmeticError('True positive heat deficit lost')
            defects=field.future_defects(Z,cells=256)
            for key in ('angular_heat_difference_scaled','pressure_heat_difference_scaled','swirl_energy_heat_difference_scaled'):
                if abs(mp.mpf(Z))==1:
                    if endpoints(defects[key][0])!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Endpoint heat defect not zero')
                elif endpoints(defects[key][0])[0]<=0:raise ArithmeticError('Positive future heat defect lost')
        for Z in ([-1,1],'.5'):
            if field.deficit(Z,4)['deficit_scaled_Taylor'].order!=1:raise ArithmeticError('Heat C1 jet lost')
        begin=field.profile('.5',0);end=field.profile('.5',3)
        if endpoints(begin['actual_bracket_deficit_scaled'][0])!=(mp.mpf(0),mp.mpf(0)):
            raise ArithmeticError('Heat collar did not match raw waiting inlet')
        if not endpoints(end['actual_bracket_deficit_scaled'][0])[0]>0:
            raise ArithmeticError('Exact heat deficit discarded at exterior start')
        # Independent finite Gamma evaluations verify the inequalities away
        # from the formal selected scale; the source itself uses proof bounds.
        finite=[]
        for aa,xx in (('.01','.001'),('.1','.05'),('.25','.1')):
            av=mp.mpf(aa);xv=mp.mpf(xx)
            H=mp.quad(lambda v:mp.exp(-v)*v**av*(1+xv*v)**(-av),[0,1,mp.inf])/mp.gamma(1+av)
            D=1-H;B=(1+av)**2*(2+av)
            if not av*(1+av)*xv-av*B*xv*xv/2<=D<=av*(1+av)*xv:
                raise ArithmeticError('Independent finite Gamma deficit inequality failed')
            finite.append(dict(a=aa,xi=xx,finite_deficit_bounds_contain_Gamma_integral=True))
    result.update(Gamma_boundary_terms_vanish_at_zero_and_infinity=True,
        positive_deficits_and_endpoint_zeros_checked=True,whole_axis_C1_heat_deficit_supported=True,
        heat_collar_endpoint_sources_checked=True,finite_Gamma_reference_examples=finite,
        exact_selected_heat_component_specified=True,actual_five_defect_family_sha256=field.angular.initial.family,
        current_source_Section7_c_epsilon_gate=False,source_not_changed_or_relabelled=True,
        incoming_five_moment_heat_matching_completed=False,actual_ap_selected=False,
        whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=dict(raw['input_hashes']))
    for source in (Path(__file__).name,name):
        result['input_hashes'][source]=hashlib.sha256((HERE/source).read_bytes()).hexdigest()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Exact heat: 9 Gamma/PDE/moment identities, positive deficits and collar endpoint sources PASS; threshold/global match unresolved',flush=True)
    return result


if __name__=='__main__':run()

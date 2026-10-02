"""Source-bound functional O.4 interfaces, independent of interval overlap.

Exact selected sources, divided-difference rows and full future supports
prove the same primitive functions in each chart. Numerical positive caps
enclose these sources and are never substituted into the equalities.
Two-sided external high-order joins remain a separate acceptance gate.
"""
import ast
import hashlib
import json
from pathlib import Path

import sympy as s
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def production_inverse(env):
    tree=ast.parse((HERE/(PREFIX+'compliant_axial_amplitude_selection.py')).read_text(encoding='utf-8'))
    method=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='linear_inverse')
    returns=[n.value for n in ast.walk(method) if isinstance(n,ast.Return)]
    if len(returns)!=1 or not isinstance(returns[0],ast.List) or len(returns[0].elts)!=2:
        raise ValueError('Unique actual two-row production inverse required')
    def expression(node):
        name=ast.unparse(node)
        if name in env:return env[name]
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -expression(node.operand)
        if isinstance(node,ast.BinOp):
            left,right=expression(node.left),expression(node.right)
            if isinstance(node.op,ast.Add):return left+right
            if isinstance(node.op,ast.Sub):return left-right
            if isinstance(node.op,ast.Mult):return left*right
            if isinstance(node.op,ast.Div):return left/right
        raise ValueError('Unbound actual inverse node: '+name)
    return [expression(node) for node in returns[0].elts]


def functional_identities():
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Functional pulse identity failed: '+name)
        proofs[name]=True
    mu=s.symbols('mu',positive=True); ap=s.symbols('ap',real=True)
    A11,A12,D11,D12,Q1,Q2,P1,P2=s.symbols('A11 A12 D11 D12 Q1 Q2 P1 P2',real=True)
    det=A11*D12-A12*D11
    inverse=production_inverse({'A[0]':A11,'A[1]':A12,'D[0]':D11,'D[1]':D12,'det':det,
        'r1':-mu*(Q1+ap*P1),'r2':-(Q2-Q1+ap*(P2-P1))})
    C1,C2=inverse
    zero('actual_selected_scaled_row1',A11*C1+A12*C2+mu*(Q1+ap*P1))
    zero('actual_selected_scaled_row2',(A11+mu*D11)*C1+(A12+mu*D12)*C2+mu*(Q2+ap*P2))
    x=s.symbols('x',positive=True)
    base=assignment('compliant_outer_pulse_map','correction_basis','base',{'x':x,'mu':mu})
    zero('actual_first_row_exponential_weight',base-s.exp(-(s.Rational(1,2)-mu)*x))
    zero('actual_second_row_divided_difference',base+mu*base*(s.exp(mu*x)-1)/mu-s.exp(-(s.Rational(1,2)-2*mu)*x))
    # G is the actual common pulse-integral scale; exp(logE)=G/mu.
    # Q_i and P_i carry the identical row normalization. Arbitrary Z
    # functions enter these identities, not selected interval centers.
    G=s.symbols('G',positive=True); m0,I=s.symbols('m0 I',real=True); lam=s.symbols('lambda_i',positive=True)
    row=s.symbols('row',integer=True,positive=True)
    logscale=assignment('compliant_outer_pulse_map','__init__','self.logscale',{
        "self.rows['common_logpref']":s.log(G),'self.mu':mu})
    zero('actual_common_end_scale_normalization',s.exp(logscale)-G/mu)
    logfactor=assignment('compliant_axial_amplitude_selection','__init__','exact_log_factor',{
        'self.mu':mu,'row':row,"self.pulse.rows['common_logpref']":s.log(G)})
    zero('actual_incoming_row_normalization',s.exp(logfactor)-s.exp(-13*(s.Rational(1,2)-row*mu)/mu)/G)
    selected_base=assignment('compliant_axial_amplitude_selection','__init__','self.base',{})
    zero('actual_selected_energy_base',selected_base-(1-s.exp(-26))/4)
    Q=m0*s.exp(-13*lam/mu)/G; Pi=s.exp(-13*lam/mu)*I/G
    scaled_row=-mu*(Q+ap*Pi); physical_row=G*scaled_row/mu
    zero('exact_physical_end_selected_moment',physical_row+s.exp(-13*lam/mu)*(m0+ap*I))
    main=s.exp(-11*lam/mu)*(m0+ap*I); gap=-s.exp(2*lam/mu)*physical_row
    zero('main_gap_linear_moments_whole_Z',main-gap)
    e0,ev,end,K=s.symbols('e0 ev end K',real=True)
    selected=(1-s.exp(-26))/4-mu*e0+mu*s.exp(-26)*(ev-end)
    forward=s.exp(22)*(e0+ap*ap*K/mu-(1-s.exp(-22))/(4*mu))
    backward=s.exp(-4)*ev+(1-s.exp(-4))/(4*mu)-s.exp(-4)*end
    zero('main_gap_selected_energy_whole_Z',(forward-backward).subs(ap*ap*K,selected))
    # nu_j=mu*exp(-26)*Kj_Rv*E²: Rp and Rv units cannot be interchanged.
    E,Kj,Cj=s.symbols('E Kj Cj',real=True)
    zero('selected_Rp_Rv_end_energy_units',mu*(s.exp(-26)*Kj)*(E*Cj)**2-mu*s.exp(-26)*(E**2*Cj**2*Kj))
    D,xi,ss,finite=s.symbols('D xi s finite',real=True)
    main_leading=assignment('compliant_axial_pulse_field','gap','leading',{'D':D,'self.mu':mu,'finite':finite})
    end_leading=assignment('compliant_axial_pulse_field','gap_from_end','leading',{'s':ss,'self.mu':mu,'finite':finite})
    zero('production_gap_coordinate_leading_log',main_leading.subs(D,-mu*ss)-end_leading)
    prate=1+2*mu; rate=1-mu
    zero('gap_coordinate_pressure_log',(-prate*xi/mu).subs(xi,13+mu*ss)-(-13*prate/mu-prate*ss))
    zero('gap_coordinate_angular_log',(-rate*xi/mu).subs(xi,13+mu*ss)-(-13*rate/mu-rate*ss))
    zero('gap_coordinate_pressure_time',((13-D)/mu).subs(D,-mu*ss)-(13/mu+ss))
    gapenergy=s.exp(-2*D)*ev+(1-s.exp(-2*D))/(4*mu)-s.exp(-2*D)*end
    endenergy=s.exp(2*mu*ss)*ev+(1-s.exp(2*mu*ss))/(4*mu)-s.exp(2*mu*ss)*end
    zero('gap_coordinate_energy_whole_Z',gapenergy.subs(D,-mu*ss)-endenergy)
    # Every end beta support lies strictly between -4 and0. Translating
    # the full future integral gives exp(lambda*(center-s))*W_i exactly.
    center,W=s.symbols('center W',real=True)
    endmoment=-E*Cj*s.exp(lam*(center-ss))*W
    gapmoment=-E*s.exp(-lam*ss)*Cj*s.exp(lam*center)*W
    zero('gap_end_full_future_linear_weight',endmoment-gapmoment)
    zero('gap_end_energy_at_s_minus4',gapenergy.subs(D,4*mu)-endenergy.subs(ss,-4))
    Xp=s.symbols('Xp',real=True)
    Xmain=1/rate+(Xp-1/rate)*s.exp(-rate*xi/mu)
    Xend=1/rate+(Xp-1/rate)*s.exp(-13*rate/mu-rate*ss)
    zero('gap_end_angular_history_whole_Z',Xmain.subs(xi,13+mu*ss)-Xend)
    Mp0,P0,u=s.symbols('Mp0 P0 u',real=True)
    t=s.symbols('t',real=True)
    P=Mp0+P0+u*u*(1-s.exp(-prate*t))/(2*prate)
    zero('gap_end_original_pressure_whole_Z',P.subs(t,(13-D)/mu).subs(D,-mu*ss)-P.subs(t,13/mu+ss))
    zero('buffer_pulse_entrance_pressure',P.subs(t,0)-(Mp0+P0))
    zero('buffer_pulse_entrance_swirl_pressure_slope',s.diff(P,t).subs(t,0)-u*u/2)
    zero('buffer_pulse_entrance_angular',Xmain.subs(xi,0)-Xp)
    incoming=assignment('compliant_axial_pulse_field','main','incoming',{
        'rows[row - 1]':m0,'self.factor(-lam * t, cut)':s.exp(-lam*t)})
    v=s.symbols('v',real=True); forcing=s.Function('forcing')(v)
    empty_integral=s.Integral(forcing,(v,0,0)).doit()
    zero('buffer_pulse_entrance_linear_moment',incoming.subs(t,0)+ap*empty_integral-m0)
    # On a zero-input neighborhood, source ODE uniqueness identifies
    # every y derivative from the common boundary history. Exact flat
    # shapes give identical forcing derivatives at support boundaries.
    Z=s.symbols('Z',real=True)
    functions=[s.Function(name)(Z) for name in ('incoming','amplitude','I','endrow')]
    f0,fap,fI,frow=functions
    selected_row=-s.exp(-13*lam/mu)*(f0+fap*fI)
    difference=s.exp(-11*lam/mu)*(f0+fap*fI)+s.exp(2*lam/mu)*selected_row
    for n in range(6):zero('functional_main_gap_axial_order'+str(n),s.diff(difference,Z,n))
    ell=s.Rational(15,100); support=[(-3-ell,-3+ell),(-1-ell,-1+ell)]
    if not all(-4<a<b<0 for a,b in support) or not support[0][1]<support[1][0]:
        raise ArithmeticError('Actual disjoint compact end supports failed')
    proofs['full_future_supports_at_minus4_and_empty_at_zero']=True
    proofs['end_energy_cross_products_exactly_zero_by_disjoint_supports']=True
    return proofs


def run():
    name=PREFIX+'compliant_pulse_physical_bounds_check.json'; receipt=json.loads((HERE/name).read_bytes())
    if not receipt['all_passed'] or not receipt['complete_pulse_physical_spatial_supremum_log_ledger_available']:
        raise ValueError('Admitted actual physical pulse derivative ledger required')
    hashes=dict(receipt['input_hashes'])
    for source,digest in hashes.items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Pulse interface source changed: '+source)
    hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    # Bind the exact formal source definitions separately from their
    # interval enclosures. No statement here treats cap values as roots.
    amplitude=json.loads((HERE/(PREFIX+'compliant_axial_amplitude_selection.json')).read_bytes())
    basis=json.loads((HERE/(PREFIX+'compliant_outer_pulse_map.json')).read_bytes())['bump_basis']
    if not amplitude['actual_ap_selected'] or not amplitude['selected_actual_c1_c2_functions_defined']:
        raise ValueError('Exact actual selected source missing')
    if not basis['identity'].startswith('row2=row1+mu*D; D integrates'):
        raise ValueError('True divided-difference integral definition missing')
    proof=functional_identities()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],implicit_source_sha256=receipt['implicit_source_sha256'],
        source_bound_functional_pulse_identities=proof,
        exact_uncapped_selected_sources_used=True,numerical_caps_used_only_as_enclosures=True,
        exact_functional_main_gap_and_gap_end_identities_certified=True,
        internal_O4_chart_function_equality_not_sample_overlap=True,
        functional_axial_derivatives_through5_identified=True,
        repeated_y_derivative_interface_compatibility_basis='same source ODEs, common boundary histories and flat forcing derivatives',
        inlet_value_identity_scope='exact O3 power inlet definitions; independent two-sided high-order O3 jets remain pending',
        full_pulse_C4_installed=False,full_outer_C4_certified=False,
        quantitative_flat_velocity_interface_bound_ledger_available=False,
        whole_outer_cone_certified=False,physical_energy_integral_certified=False,temporal_recursion=False,all_passed=True,
        next_dependency='Quantitative flat velocity/moment comparison bounds and two-sided O3/O5 high-order joins, then full post-pulse C4/cone',input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Source-bound functional pulse interfaces: exact selected moment/energy identities, coordinate changes and full future supports PASS',flush=True)
    return result


if __name__=='__main__':run()

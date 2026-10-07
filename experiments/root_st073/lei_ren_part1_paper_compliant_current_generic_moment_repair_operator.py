"""New Rc..2Rc five-bump functional inverse for arbitrary generic defects.

The actual incoming source functions remain a separate prerequisite. This
operator proves a uniform C1(Z) implicit inverse under their whole-source
log bounds. It does not replace a defect by a midpoint or assume the old
special O3 axial correlation. Exact integral weights remain functions;
their directed enclosures supply constants only.
"""
import ast
import json
from pathlib import Path
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_five_defect_bounds as defects
from lei_ren_part1_paper_compliant_outer_pulse_map import raw_beta

HERE,PREFIX,sha=defects.HERE,defects.PREFIX,defects.sha
packets=defects.packets;LogUpper=defects.LogUpper
NAME=PREFIX+'current_generic_moment_repair_operator.json'
RECEIPT=PREFIX+'current_generic_moment_repair_operator_check.json'
GATE='current_generic_Rc_2Rc_five_bump_C1_functional_inverse_log_conditions_certified'
OPEN=defects.OPEN
ROWS=('M','D=(J-M)/mu','I','S','Cp')
CONTROLS=('axial0','axial2','swirl0','swirl1','swirl2')


def upper(c,x):return c.mpf(packets.recovery.endpoints(x)[1])


def exp_average(c,x):
    """Enclose exact integral_0^1 exp(t*x)dt, including x=0."""
    lo,hi=packets.recovery.endpoints(x)
    return c.exp(c.mpf([min(0,lo),max(0,hi)]))


def fresh_weights(c,mu,cells=256):
    """Fresh ln2 geometry; no old repair coefficient or normalizer value."""
    ep=packets.recovery.endpoints;L=c.ln(2);ell=L/40
    centers=[L/5,L/2,4*L/5]
    normal=c.mpf(0);H={key:c.mpf(0) for key in ('I','S','Cp','axial')}
    gram={key:c.mpf(0) for key in ('cross','energy','pressure')}
    divided=[c.mpf(0),c.mpf(0)]
    cells_y=[c.mpf([ep(c.mpf(-1)+c.mpf(2*i)/cells)[0],ep(c.mpf(-1)+c.mpf(2*(i+1))/cells)[1]]) for i in range(cells)]
    raw=[raw_beta(c,y) for y in cells_y]
    for beta in raw:normal+=beta*2/cells
    if ep(normal)[0]<=0:raise ArithmeticError('Exact raw bump normalization needs positive enclosure')
    for y,beta in zip(cells_y,raw):
        w=ell*y;mass=beta*(c.mpf(2)/cells)/normal
        square=beta*beta*(c.mpf(2)/cells)/(ell*normal*normal)
        for key,power in (('I',c.mpf('.5')),('S',-c.mpf('.5')-mu),('Cp',-c.mpf('1.5')-mu),('axial',-mu)):
            H[key]+=mass*c.exp(power*w)
        for key,power in (('cross',c.mpf('-.5')),('energy',c.mpf(-1)),('pressure',c.mpf(-2))):
            gram[key]+=square*c.exp(power*w)
        for j,center in enumerate((centers[0],centers[2])):
            t=center+w;divided[j]+=mass*(-t)*exp_average(c,-mu*t)
    if any(ep(value)[0]<=0 for value in (*H.values(),*gram.values())):
        raise ArithmeticError('Positive new integral weights required')
    return dict(log_band_length=L,centers=centers,radius=ell,raw_normalization=normal,
        H=H,gram=gram,divided_axial_rows=divided,cells=cells,
        exact_first_axial_row=[1,1],same_exact_raw_integral_normalization=True,
        interval_weights_are_enclosures_not_chosen_defining_values=True)


def fresh_linear_inverse(c,mu,W):
    ci=W['centers'];H=W['H'];gap=ci[2]-ci[0];step=ci[1]-ci[0]
    D=W['divided_axial_rows']
    detA=-gap*H['axial']*c.exp(-mu*ci[0])*exp_average(c,-mu*gap)
    inverseA=[[D[1]/detA,-1/detA],[-D[0]/detA,1/detA]]
    powers=[c.mpf('.5'),-c.mpf('.5')-mu,-c.mpf('1.5')-mu];signs=[1,-1,1]
    keys=('I','S','Cp')
    E=[[signs[i]*H[keys[i]]*c.exp(powers[i]*center) for center in ci] for i in range(3)]
    factors=[signs[i]*H[keys[i]]*c.exp(powers[i]*ci[0]) for i in range(3)]
    detE=factors[0]*factors[1]*factors[2]
    # Difference is kept correlated before enclosure; no rounded nodes
    # are subtracted and mu is not removed from the defining map.
    for i,j in ((0,1),(0,2),(1,2)):
        detE*=c.exp(powers[i]*step)*(powers[j]-powers[i])*step*exp_average(c,(powers[j]-powers[i])*step)
    ep=packets.recovery.endpoints
    if ep(detA)[1]>=0 or ep(detE)[0]<=0:raise ArithmeticError('New divided axial/swirl determinant signs unresolved')
    inverseE=[]
    for i in range(3):
        row=[]
        for j in range(3):
            ii=[k for k in range(3) if k!=j];jj=[k for k in range(3) if k!=i]
            numerator=(E[ii[0]][jj[0]]*E[ii[1]][jj[1]]-E[ii[0]][jj[1]]*E[ii[1]][jj[0]])*((-1)**(i+j))
            row.append(numerator/detE)
        inverseE.append(row)
    B=[[c.mpf(1),c.mpf(1)]+[c.mpf(0)]*3,D+[c.mpf(0)]*3]
    B += [[c.mpf(0)]*2+row for row in E]
    inv=[[c.mpf(0)]*5 for _ in range(5)]
    for i in range(2):inv[i][:2]=inverseA[i]
    for i in range(3):inv[i+2][2:]=inverseE[i]
    norm=upper(c,c.mpf(max(ep(sum((upper(c,abs(v)) for v in row),c.mpf(0)))[1] for row in inv)))
    gram=W['gram']
    return dict(linear_enclosure=B,inverse_enclosure=inv,inverse_infinity_norm_upper=norm,
        divided_axial_determinant=detA,swirl_determinant=detE,
        cross_weights=[gram['cross']*c.exp(-ci[i]/2) for i in (0,2)],
        energy_weights=[gram['energy']*c.exp(-center) for center in ci],
        pressure_weights=[gram['pressure']*c.exp(-2*center) for center in ci],
        actual_inverse_is_of_exact_integral_matrix=True,
        source_mu_cover_includes_zero_only_for_analytic_extension=True,
        exact_first_two_rows_are_M_and_divided_J_minus_M=True)


def transformed_quadratic(c,mu,N,matrix,h):
    """General repair controls F=e*g/N, G=a*g/N; no mu prefactor."""
    a0,a2,e0,e1,e2=h;energy=matrix['energy_weights'];pressure=matrix['pressure_weights'];zero=a0*0
    return [zero,(matrix['cross_weights'][0]*a0*e0+matrix['cross_weights'][1]*a2*e2)/(mu*N),
        zero,(energy[0]*a0*a0+energy[2]*a2*a2-(energy[0]*e0*e0+energy[1]*e1*e1+energy[2]*e2*e2)/2)/N,
        (pressure[0]*e0*e0+pressure[1]*e1*e1+pressure[2]*e2*e2)/(2*N)]


def target_C1_bounds(c,inlet,raw,log_A_lower,log_mu_lower):
    """Arbitrary incoming dJ-dM; its inverse mu is never canceled by a claim."""
    C=lambda v:LogUpper.constant(c,v);read=lambda row:LogUpper(c,None if row['exact_zero'] else packets.interval(c,row['log_absolute_upper']))
    logAZ=read(raw['E']['y0_Z1'])
    ratio=logAZ.divide_positive(log_A_lower)
    polys={}
    for key,power in (('m',1),('h',1),('k',2),('e',2),('p',2)):
        value=defects.read_poly(c,inlet[key]['value']);Z=defects.read_poly(c,inlet[key]['Z'])
        scale=LogUpper(c,-power*log_A_lower)
        polys[key]=(value+Z+value.scale(ratio*C(power))).scale(scale)
    # N*d has only powers0,-1; for N>=1 take the sum of its actual
    # log coefficient bounds without evaluating N or its exponential.
    def bounded(poly):
        if any(power>=0 for power in poly.terms):raise ArithmeticError('N-scaled incoming defect bound needs negative original powers')
        return LogUpper.add(c,list(poly.terms.values()))
    normal={key:bounded(poly) for key,poly in polys.items()}
    rowcaps=[normal['m'],LogUpper.add(c,[normal['k'],normal['m']]).divide_positive(log_mu_lower),
             normal['h'],normal['e'],normal['p']]
    ep=packets.recovery.endpoints
    maximum=LogUpper(c,c.mpf(max(ep(cap.log)[1] for cap in rowcaps if cap.log is not None)))
    return dict(normalized_N_scaled_inlet_C1_cap={key:cap.record() for key,cap in normal.items()},
        transformed_N_scaled_target_C1_caps={key:cap.record() for key,cap in zip(ROWS,rowcaps)},
        whole_target_C1_cap=maximum.record(),original_A_Z_over_A_log_cap=ratio.record(),
        norm='max_i(sup_Z|v_i|+sup_Z|d_Z v_i|), all Z in[-1,1]',
        quotient_C1_formula='||(D/A^p)||C1 <= A_lower^-p*(|D|+|D_Z|+p*|A_Z/A|*|D|)',
        generic_axial_difference_not_assumed_O_mu=True,
        exact_transformed_target='N*(D_m/A,(D_k/A^2-D_m/A)/mu,D_h/A,D_e/A^2,D_p/A^2)',
        actual_target_signed_functions_not_replaced_by_caps=True)


def contraction_log_conditions(c,target,matrix,weights,log_mu_lower,source_logN_lower):
    C=lambda value:LogUpper.constant(c,value)
    D=LogUpper(c,packets.interval(c,target['whole_target_C1_cap']['log_absolute_upper']))
    CA=C(matrix['inverse_infinity_norm_upper'])
    rows=[C(sum(matrix['cross_weights'])).divide_positive(log_mu_lower),
          C(matrix['energy_weights'][0]+matrix['energy_weights'][2]+sum(matrix['energy_weights'])/2),
          C(sum(matrix['pressure_weights'])/2)]
    ep=packets.recovery.endpoints
    CQ=LogUpper(c,c.mpf(max(ep(cap.log)[1] for cap in rows)))
    radius=D*CA*C(2);threshold=CA*CQ*radius*C(4)
    # Disjoint compact bumps give |F|<=g_max*rho/N. Require the swirl
    # correction to be below half the original power floor2^(-2/3).
    gmax=C(c.exp(-1)/(weights['radius']*weights['raw_normalization']))
    positivity=radius*gmax*C(2)*LogUpper(c,c.mpf(2)/3*c.ln(2))
    logN=c.mpf(max(0,ep(source_logN_lower)[1],ep(threshold.log)[1],ep(positivity.log)[1]))
    return dict(actual_target_C1_log_cap=D.record(),exact_integral_matrix_inverse_log_cap=CA.record(),
        general_transformed_quadratic_C1_log_cap=CQ.record(),
        formal_control_C1_ball_radius_log=radius.record(),
        repair_sufficient_common_log_N_lower=logN,
        compact_bump_profile_log_upper=gmax.record(),
        swirl_positivity_sufficient_log_N_lower=positivity.log,
        sufficient_recipe='rho=2*CA*D; N>=max(1, source_expm1_N,4*CA*CQ*rho,2*g_max*rho*2^(2/3))',
        contraction_at_most='1/2',image_radius_at_most='3*rho/4',
        exact_implicit_function='h(Z)=-B(mu)^-1*(d_scaled(Z)+Q(mu,h(Z))/N)',
        unique_C1_control_functions_in_certified_ball_for_any_compatible_actual_C1_target=True,
        repaired_swirl_positive_lower_relative_to_A='2^(-2/3)/2',
        B_mu_Z_independent_and_C1_product_norm_submultiplicative=True,
        no_old_special_correlated_defect_or_old_finite_N_used=True,
        actual_finite_integer_N_not_selected=True,
        actual_signed_defect_functions_and_control_functions_not_evaluated=True)


def exact_theorem():
    mu,L,x,A,S,Rc,N=s.symbols('mu L x A S Rc N',positive=True)
    m,h,k,e,p,F,G=s.symbols('D_m D_h D_k D_e D_p F G',real=True);checks={}
    alpha=s.Rational(1,2)+mu
    density={'m':G,'h':s.sqrt(x)*F,'k':s.sqrt(x)*(x**(-alpha)+F)*G,
             'e':G*G-x**(-alpha)*F-F*F/2,'p':(x**(-alpha)*F+F*F/2)/x}
    # Convert each normalized cumulative ODE to its physical moment
    # integral; the original E=Ac/Pstar and V=0 on this power band.
    in_units={'m':A,'h':A,'k':A*A,'e':A*A,'p':A*A}
    rates={'m':1,'h':s.Rational(3,2),'k':s.Rational(3,2),'e':1,'p':0}
    delta=packets.recovery.increment_densities(A*x**(-alpha),s.Integer(0),A*F,A*G)
    for key in density:
        actual=x**(rates[key]-1)*delta[key]/in_units[key]
        if s.simplify(s.expand(actual-density[key]))!=0:raise ArithmeticError('New general repair physical density differs: '+key)
        checks['normalized_physical_density_'+key]=True
    beta,t=s.symbols('same_normalized_beta t',real=True)
    if s.simplify(s.sqrt(x)*x**(-alpha)-x**(-mu))!=0:raise ArithmeticError('Axial row linear weight differs')
    checks['new_axial_J_linear_weight']=True
    # Exact divided row keeps every generic cross term /mu, unlike the
    # old special sqrt(mu),mu scaled source/controls.
    Jlinear=x**(-mu)*G
    divided=(density['k']-density['m'])/mu
    expected=(Jlinear-G)/mu+s.sqrt(x)*F*G/mu
    if s.simplify(s.expand(divided-expected))!=0:raise ArithmeticError('Generic divided row omitted inverse mu')
    checks['general_divided_J_minus_M_cross_over_mu']=True
    rho,CA,CQ,D=s.symbols('rho CA CQ D',positive=True)
    r0=2*CA*D;N0=4*CA*CQ*r0
    if s.simplify(2*CA*CQ*r0/N0-s.Rational(1,2))!=0 or s.simplify(CA*D+CA*CQ*r0*r0/N0-3*r0/4)!=0:
        raise ArithmeticError('C1 Banach image/contraction recipe differs')
    checks['C1_contraction_half_and_image_three_quarters']=True
    u,v,w=s.symbols('u v w');V=s.Matrix([[1,r,r*r] for r in (u,v,w)])
    if s.expand(V.det()-(v-u)*(w-u)*(w-v))!=0:raise ArithmeticError('Fresh swirl Vandermonde determinant differs')
    checks['fresh_equal_log_center_Vandermonde_determinant']=True
    Hax,c0,c2=s.symbols('same_Hax c0 c2',real=True)
    d0=(Hax*s.exp(-mu*c0)-1)/mu;d2=(Hax*s.exp(-mu*c2)-1)/mu
    det=Hax*s.exp(-mu*c0)*(s.exp(-mu*(c2-c0))-1)/mu
    if s.simplify(d2-d0-det)!=0:raise ArithmeticError('Fresh divided axial determinant differs')
    checks['fresh_divided_axial_determinant_analytic_extension']=True
    tree=ast.parse((HERE/(PREFIX+'outer_pulse_map.py')).read_text(encoding='utf8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='raw_beta')
    if not any(isinstance(n,ast.Assign) and any(ast.unparse(q)=='lower' for q in n.targets)
       and ast.dump(n.value)==ast.dump(ast.parse('c.exp(-1/(1-c.mpf(far)**2)) if far<1 else c.mpf(0)',mode='eval').body) for n in ast.walk(fn)):
        raise ValueError('Original exact raw bump changed')
    checks['same_exact_raw_beta_definition_AST_bound']=True
    return dict(passed=True,identities=checks,
        exact_bump_normalization='J0=integral_-1^1 raw_beta(t)dt; beta_ell(w)=raw_beta(w/ell)/(ell*J0)',
        exact_profile='E=A(Z)*(x^(-1/2-mu)+sum e_i(Z)g_i(x)/N), V=A(Z)*sum a_i(Z)g_i(x)/N',
        exact_terminal_condition='d+(M,J,I,S,Cp)=0; normalized moment propagation multiplies both inlet and integral by2^-rate',
        actual_original_incoming_histories_and_P0_unchanged=True,
        new_defect_functions_must_be_actual_same_family_C1_functions=True)


class CurrentGenericMomentRepairOperator:
    def __init__(self):
        self.ctx=MPIntervalContext();self.ctx.dps=240
        receipt=json.loads((HERE/defects.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[defects.GATE] or any(receipt.get(k) for k in OPEN):
            raise ValueError('Checked current whole five-defect envelope required')
        self.hashes=dict(receipt['input_hashes'])
        for filename,digest in self.hashes.items():
            if sha(filename)!=digest:raise ValueError('Changed generic repair prerequisite: '+filename)
        self.hashes[defects.RECEIPT]=sha(defects.RECEIPT)
        self.incoming=json.loads((HERE/defects.NAME).read_bytes());self.family=self.incoming['source_family']
        name=PREFIX+'current_generic_shear_O3_sources.json';self.O3=json.loads((HERE/name).read_bytes())
        if self.O3['source_family']!=self.family:raise ValueError('Same Ac/P0/original power source required')
        check=PREFIX+'current_generic_shear_O3_sources_check.json';O3check=json.loads((HERE/check).read_bytes())
        if not O3check['all_passed'] or O3check['source_family']!=self.family or any(O3check.get(k) for k in OPEN) or not all(O3check.get(k) for k in (
            'current_original_O3_common_unit_packets_and_quotient_log_bounds_certified',
            'current_original_O3_right_edge_and_reserved_repair_geometry_certified')):
            raise ValueError('Checked original O3 raw packets and Rc reservation required')
        self.hashes[name]=sha(name);self.hashes[check]=sha(check)
        self.reserve=self.O3['right_edge_and_new_repair_reservation']
        if self.reserve['original_reserved_profile']!='Utheta=Ac_theta(Z)*(R/Rc)^(-1/2-mu), Uz=0':
            raise ValueError('Same original Rc power amplitude function required')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.hashes[PREFIX+'outer_pulse_map.py']=sha(PREFIX+'outer_pulse_map.py')
        self.hashes[PREFIX+'current_generic_left_inlet.py']=sha(PREFIX+'current_generic_left_inlet.py')
        self.theorem=exact_theorem()

    def compute(self):
        c=self.ctx;ep=packets.recovery.endpoints;read=lambda row:packets.interval(c,row)
        positive=self.O3['original_O3_quotient_log_norms']['O3_power']['actual_positive_denominator_theorem']
        mu=read(positive['actual_positive_mu'])
        if not 0<ep(mu)[0]<=ep(mu)[1]<ep(c.mpf(1)/6)[0]:raise ArithmeticError('Same original positive mu<1/6 required')
        log_mu_lower=c.ln(c.mpf(ep(mu)[0]));cover=c.mpf([0,ep(c.mpf(1)/6)[1]])
        weights=fresh_weights(c,cover);matrix=fresh_linear_inverse(c,cover,weights)
        logA=read(self.reserve['positive_Ac_over_S_log_lower'])
        target=target_C1_bounds(c,self.incoming['actual_new_repair_inlet_Rc_pre_repair_defect_majorants'],
            self.O3['original_O3_quotient_log_norms']['O3_power']['ordinary_mixed_source_log_norms'],logA,log_mu_lower)
        certificate=contraction_log_conditions(c,target,matrix,weights,log_mu_lower,read(self.incoming['required_positive_log_N_lower']))
        return dict(source_family=self.family,**{GATE:True},**dict.fromkeys(OPEN,False),
            new_repair_geometry=dict(R_left='Rc',R_right='2Rc',x='R/Rc in[1,2]',log_length=weights['log_band_length'],
                log_centers=weights['centers'],log_radius=weights['radius'],
                all_bump_supports_strictly_inside_original_power_band=True,
                original_Ac_amplitude_source=self.reserve['original_Ac_definition'],
                positive_Ac_over_Pstar_log_lower=logA,actual_positive_mu_log_lower=log_mu_lower,
                original_actual_mu_not_replaced_by_uniform_cover=True,
                actual_mu_cover_for_constants_only=cover),
            fresh_exact_weight_definitions_and_enclosures=weights,
            fresh_divided_linear_inverse_and_enclosures=matrix,
            actual_generic_defect_target_C1_bounds=target,
            actual_generic_repair_C1_log_contraction_conditions=certificate,
            exact_new_general_five_moment_repair_theorem=self.theorem,
            row_order=ROWS,control_order=CONTROLS,
            explicit_original_left_inlet_query_adapter=dict(module=PREFIX+'current_generic_left_inlet.py',
                callable='original_left_inlet(existing_checked_bridge_owner,Z=...)',
                original_inlet_coordinate='bridge_first at s_c/2, r_minus=Ra*exp(hb*s_c/2)',
                no_owner_constructor_or_saved_cover_fallback=True,
                source_rows_remain_genuine_query_covers_not_selected_values=True,
                actual_existing_owner_success_path_exercised=False),
            implicit_function_contract=dict(unknowns='five C1(Z) functions h=(a0,a2,e0,e1,e2)',
                equation='B_exact(mu)*h+d_scaled(Z)+Q_exact(mu,h)/N=0',
                exact_weights_are_raw_bump_integrals_not_enclosure_endpoints=True,
                incoming_signed_generic_defect_functions_required=True,
                old_special_O3_defects_or_axial_correlation_not_used=True,
                actual_controls_installed=False,actual_terminal_Z_function_closure_installed=False),
            sufficient_mixed4_velocity_or_mixed3_stress_admitted=False,
            current_whole_N_selected=False,actual_changed_history_integration_or_repair_evaluated=False,
            source_graph_ancestor_constructors_called=False,
            scope='New(Rc,2Rc) five-bump map and uniform divided inverse for arbitrary same-family incoming generic defects, actual C1 target log bounds and sufficient log-N Banach conditions. Conditional functional operator only: actual source replay/seams/integrals/controls, high derivatives, terminal field closure/global cone/recursion remain open.',
            input_hashes=self.hashes)

    def run(self):
        row=self.compute();(HERE/NAME).write_text(json.dumps(packets.encode(row),indent=2)+'\n',encoding='utf8')
        print('New Rc..2Rc generic five-bump C1 inverse/log conditions defined; arbitrary axial defect difference retained',flush=True)
        return row


def run():return CurrentGenericMomentRepairOperator().run()


if __name__=='__main__':run()

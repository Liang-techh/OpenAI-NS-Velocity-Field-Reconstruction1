"""Fresh repair geometry, exact nonlinear density fixtures and C1 scope."""
import json
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_generic_moment_repair_operator as source
from lei_ren_part1_paper_compliant_current_generic_left_inlet import original_left_inlet


def independent_fixture(manifest):
    c=mp.mp.clone();c.dps=65;L=c.ln(2);ell=L/40;centers=[L/5,L/2,4*L/5]
    beta=lambda t:c.exp(-1/(1-t*t)) if abs(t)<1 else c.mpf(0)
    J0=c.quad(beta,[-1,0,1]);comparisons=0;cases=[]
    ep=source.packets.recovery.endpoints;iv=source.MPIntervalContext();iv.dps=100
    cover=lambda row:source.packets.interval(iv,row)
    enclosed=manifest['fresh_divided_linear_inverse_and_enclosures']
    for mu_text,N_text in (('.04','1e12'),('1e-8','1e32')):
        mu=c.mpf(mu_text);N=c.mpf(N_text);alpha=c.mpf('.5')+mu
        quad=lambda f:c.quad(f,[-1,0,1])
        H={key:quad(lambda u,p=p:beta(u)/J0*c.exp(p*ell*u)) for key,p in
            (('I',c.mpf('.5')),('S',-alpha),('Cp',-alpha-1),('axial',-mu))}
        divided=[quad(lambda u,center=center:beta(u)/J0*c.expm1(-mu*(center+ell*u))/mu) for center in (centers[0],centers[2])]
        gram={key:quad(lambda u,p=p:beta(u)**2/(ell*J0**2)*c.exp(p*ell*u)) for key,p in
            (('cross',c.mpf('-.5')),('energy',c.mpf(-1)),('pressure',c.mpf(-2)))}
        B=c.matrix(5);B[0,0]=1;B[0,1]=1;B[1,0]=divided[0];B[1,1]=divided[1]
        for row,key,power,sign in ((2,'I',c.mpf('.5'),1),(3,'S',-alpha,-1),(4,'Cp',-alpha-1,1)):
            for j,center in enumerate(centers):B[row,j+2]=sign*H[key]*c.exp(power*center)
        inverse=B**-1
        matrix=dict(cross_weights=[gram['cross']*c.exp(-centers[i]/2) for i in (0,2)],
            energy_weights=[gram['energy']*c.exp(-center) for center in centers],
            pressure_weights=[gram['pressure']*c.exp(-2*center) for center in centers])
        for i in range(5):
            for j in range(5):
                lo,hi=ep(cover(enclosed['linear_enclosure'][i][j]))
                if not lo<=B[i,j]<=hi:raise ArithmeticError('Fresh exact scalar matrix not enclosed')
                lo,hi=ep(cover(enclosed['inverse_enclosure'][i][j]))
                if not lo<=inverse[i,j]<=hi:raise ArithmeticError('Fresh exact inverse not enclosed')
                comparisons+=2
        def original_target(Z):
            # N-scaled ordinary dJ and dM are independent affine data.
            # Their difference is not divided away using any old source law.
            return [c.mpf('.01')+c.mpf('.003')*Z,-c.mpf('.02')+c.mpf('.001')*Z,
                c.mpf('.01')-c.mpf('.002')*Z,-c.mpf('.04')+c.mpf('.003')*Z,c.mpf('.002')+c.mpf('.001')*Z]
        def solve(Z):
            d=original_target(Z);target=c.matrix([d[0],(d[1]-d[0])/mu,*d[2:]])
            h=-inverse*target
            for _ in range(20):
                h=-inverse*(target+c.matrix(source.transformed_quadratic(c,mu,N,matrix,list(h))))
            residual=B*h+target+c.matrix(source.transformed_quadratic(c,mu,N,matrix,list(h)))
            if max(abs(v) for v in residual)>c.mpf('1e-45')*(1+max(abs(v) for v in target)):
                raise ArithmeticError('Modest arbitrary-defect nonlinear inverse did not converge')
            return h,target
        for Z in (c.mpf('-.6'),c.mpf('.37')):
            h,target=solve(Z);actual=B*h+c.matrix(source.transformed_quadratic(c,mu,N,matrix,list(h)))
            flux=[c.mpf(0)]*5
            # Independently integrate the original full signed density
            # differences; no response matrix or quadratic formula used.
            for i,center in enumerate(centers):
                def integrand(u,key):
                    t=center+ell*u;x=c.exp(t);g=beta(u)/(ell*J0*x)
                    G=(h[0] if i==0 else h[1] if i==2 else c.mpf(0))*g/N
                    F=h[i+2]*g/N;E0=x**(-alpha);dx=ell*x
                    f=[G,c.sqrt(x)*(E0+F)*G,c.sqrt(x)*F,
                        G*G-E0*F-F*F/2,(E0*F+F*F/2)/x][key]
                    return N*f*dx
                for key in range(5):flux[key]+=quad(lambda u,key=key:integrand(u,key))
            transformed=c.matrix([flux[0],(flux[1]-flux[0])/mu,*flux[2:]])
            for value,expected in zip(transformed,actual):
                if abs(value-expected)>c.mpf('1e-45')*(1+abs(expected)):
                    raise ArithmeticError('New operator differs from full physical density integrals')
                comparisons+=1
            for value,expected in zip(transformed,-target):
                if abs(value-expected)>c.mpf('1e-45')*(1+abs(expected)):
                    raise ArithmeticError('Independent five terminal density conditions failed')
                comparisons+=1
            # Numerical Z derivative checks the implicit function rather
            # than only a collection of independent point fits.
            jac=c.matrix([[c.diff(lambda v:source.transformed_quadratic(c,mu,N,matrix,
                [v if k==j else h[k] for k in range(5)])[i],h[j]) for j in range(5)] for i in range(5)])
            td=c.matrix([c.mpf('.003'),-c.mpf('.002')/mu,-c.mpf('.002'),c.mpf('.003'),c.mpf('.001')])
            hz=-(B+jac)**-1*td;step=c.mpf('1e-4')
            hm2=solve(Z-2*step)[0];hm1=solve(Z-step)[0];hp1=solve(Z+step)[0];hp2=solve(Z+2*step)[0]
            fd=(-hp2+8*hp1-8*hm1+hm2)/(12*step)
            for value,expected in zip(hz,fd):
                if abs(value-expected)>c.mpf('1e-15')*(1+abs(expected)):raise ArithmeticError('Implicit C1 controls differ from Z derivative')
                comparisons+=1
            cases.append(dict(mu=mu_text,N=N_text,Z=c.nstr(Z,5),arbitrary_dJ_minus_dM_without_old_correlation=True))
    return dict(comparisons=comparisons,cases=cases,
        full_signed_density_integrals_vs_map_allowance='1e-45 relative/absolute, modest fixtures only',
        implicit_Z_derivative_fourth_order_difference_allowance='1e-15, step1e-4',
        original_source_scale_or_defect_functions_not_materialized=True,
        fixture_integral_weights_are_independent_raw_beta_quadratures=True,
        old_special_O3_defect_or_control_scaling_not_used=True)


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes());provider=source.CurrentGenericMomentRepairOperator()
    if manifest!=source.packets.encode(provider.compute()):raise ValueError('New functional inverse differs from actual source attachment')
    if not manifest[source.GATE] or any(manifest[k] for k in source.OPEN):raise ArithmeticError('Conditional repair promoted global completion')
    contract=manifest['implicit_function_contract'];certificate=manifest['actual_generic_repair_C1_log_contraction_conditions']
    if contract['actual_controls_installed'] or contract['actual_terminal_Z_function_closure_installed']:
        raise ArithmeticError('Actual incoming functions/controls still missing')
    if not manifest['actual_generic_defect_target_C1_bounds']['generic_axial_difference_not_assumed_O_mu']:
        raise ArithmeticError('Generic target got the old special axial correlation')
    if certificate['actual_finite_integer_N_not_selected'] is not True or certificate['no_old_special_correlated_defect_or_old_finite_N_used'] is not True:
        raise ArithmeticError('Repair-only threshold was mistaken for whole-source finite N')
    c=provider.ctx;ep=source.packets.recovery.endpoints;read=lambda row:source.packets.interval(c,row)
    geometry=manifest['new_repair_geometry'];L=read(geometry['log_length']);ell=read(geometry['log_radius']);ci=[read(row) for row in geometry['log_centers']]
    if ep(ci[0]-ell)[0]<=0 or ep(L-ci[-1]-ell)[0]<=0 or any(ep(ci[i+1]-ci[i]-2*ell)[0]<=0 for i in range(2)):
        raise ArithmeticError('New bumps overlap or reach either endpoint')
    norm=read(certificate['exact_integral_matrix_inverse_log_cap']['log_absolute_upper'])
    Q=read(certificate['general_transformed_quadratic_C1_log_cap']['log_absolute_upper'])
    rho=read(certificate['formal_control_C1_ball_radius_log']['log_absolute_upper'])
    minimum=read(certificate['repair_sufficient_common_log_N_lower'])
    if ep(minimum-(norm+Q+rho+c.ln(4)))[0]<0:
        # Comparing two directed enclosures of the same expression may
        # cross zero by rounding. Compare certified upper endpoints.
        if ep(minimum)[0]<ep(norm+Q+rho+c.ln(4))[1]:raise ArithmeticError('Repair contraction log threshold insufficient')
    if ep(minimum)[0]<ep(read(certificate['swirl_positivity_sufficient_log_N_lower']))[1]:raise ArithmeticError('Repair swirl positivity threshold insufficient')
    guards=0
    for owner in (None,object()):
        try:original_left_inlet(owner)
        except ValueError:guards+=1
        else:raise ArithmeticError('Left-inlet adapter accepted missing/foreign owner')
    fixture=independent_fixture(manifest)
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=source.sha(source.NAME);hashes[Path(__file__).name]=source.sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=provider.family,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
        exact_new_general_repair_identities=len(provider.theorem['identities']),
        fresh_ln2_band_disjoint_flat_bump_geometry=True,uniform_positive_mu_divided_axial_and_swirl_inverse=True,
        arbitrary_generic_axial_difference_inverse_mu_and_quadratic_inverse_mu_retained=True,
        actual_A_Z_over_A_and_pressure_P0_preserved=True,
        C1_Banach_image_contraction_and_swirl_positivity_log_conditions=True,
        independent_general_repair_fixture=fixture,
        explicit_left_inlet_missing_or_foreign_owner_guards=guards,
        actual_left_inlet_existing_owner_success_path_exercised=False,
        actual_signed_incoming_defect_functions_or_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        sufficient_mixed4_velocity_or_mixed3_stress_admitted=False,
        source_graph_ancestor_constructors_called=False,input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('New Rc..2Rc arbitrary-generic five-bump C1 inverse PASS: fresh matrix, full signed density and implicit Z kernels',flush=True)
    return result


if __name__=='__main__':run()

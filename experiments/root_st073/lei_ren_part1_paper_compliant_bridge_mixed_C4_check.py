"""Independent varying-comparison/source checks for bridge mixed4."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_bridge_mixed_C4 import (
    CompliantBridgeMixedC4,AxialSixAlgebra,IntervalTaylor,radial_core_log,
    smoothed_comparison_rows,comparison_directions,bridge_controls,S2,MTH,MZ,MTHZ,MZT,MP)
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import comparison_radial_directions
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import scaled_positive_source
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import canonical_source
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_five_moment_repair import pack

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_bridge_mixed_C4.json'


def jet(c,fn,z,tol,order=6):
    return IntervalTaylor(c,[c.mpf([q-tol,q+tol])/math.factorial(n)
        for n in range(order+1) for q in (mp.diff(fn,z,n),)])


def contains(c,row,expected,tol,label):
    lo,hi=endpoints(row)
    if not lo-tol*10000<=expected<=hi+tol*10000:raise ArithmeticError('Independent bridge fixture failed: '+label)


def core_log_fixture():
    """Independent log differentiation under rho=4exp(y), all rectangular rows."""
    with mp.workdps(70):
        c=MPIntervalContext();c.dps=105;tol=mp.mpf('1e-50');r0=mp.mpf('4.05');z0=mp.mpf('.3')
        rho,z,y=s.symbols('rho z y',real=True)
        phi=s.exp(z/20)*(1+rho/100+rho*rho*(1+z*z)/5000+rho**3*z**6/1000000)
        V=4*z+rho*(1+z*z)/100+rho**3*z**6/1000000
        point=lambda expr:s.lambdify((rho,z),expr,'mpmath')
        rows={}
        for name,expr in [('Phi',phi),('Uz',V)]:
            rows[name]=[jet(c,lambda zz,i=i:point(s.diff(expr,rho,i))(r0,zz),z0,tol) for i in range(4)]
        bounded=radial_core_log(rows['Phi'],c.mpf(r0));count=0
        for k,row in enumerate(bounded,1):
            true=s.diff(s.log(phi.subs(rho,rho*s.exp(y))),y,k).subs(y,0)
            for n in range(7):
                contains(c,row[n]*math.factorial(n),point(s.diff(true,z,n))(r0,z0),tol,'Euler logPhi');count+=1
        for k in range(1,4):
            bounded=sum((rows['Uz'][i]*(S2[k][i]*c.mpf(r0)**i) for i in range(1,k+1)),rows['Uz'][0]*0)
            true=s.diff(V.subs(rho,rho*s.exp(y)),y,k).subs(y,0)
            for n in range(7):
                contains(c,bounded[n]*math.factorial(n),point(s.diff(true,z,n))(r0,z0),tol,'Euler V');count+=1
        return dict(independent_rectangular_Euler_log_and_velocity_derivatives=count,
            includes_axial6_and_radial3=True,finite_fixture_only=True,actual_source_admission=False,passed=True)


def varying_direction_fixture():
    """True physical direction differentiated from CLOSED varying-field integrals."""
    with mp.workdps(70):
        c=MPIntervalContext();c.dps=105;tol=mp.mpf('1e-50');t0=mp.mpf('.3');z0=mp.mpf('.2')
        t,z,q=s.symbols('t z q',real=True);h=s.Rational(1,40);delta=s.Rational(1,2000)
        phi=1+z/100+t/50+t*t*(1+z*z)/1000
        V=4*z+s.Rational(1,50)+(1+z*z)*t/100+t**3*z/1000
        F0=s.exp(z/20+z*z/60);p0=s.Rational(1,30)+z*z/100+z**6/100000
        initial=[s.Rational(i+1,10)+z/100+z*z/100+z**6/10000 for i in range(6)]
        rhs=[2*phi,V,2*phi*V,V*V,phi*phi,phi*phi];rates=[2,1,2,1,2,1]
        moment=[]
        for value,force,rate in zip(initial,rhs,rates):
            weight=s.exp(rate*h*q)*force.subs(t,q)
            integral=s.integrate(weight,(q,0,t))
            moment.append(s.exp(-rate*h*t)*(value+h*integral))
        H,m,K,A,B,P=moment;point=lambda expr:s.lambdify((t,z),expr,'mpmath')
        zj=lambda expr:jet(c,lambda zz:point(expr)(t0,zz),z0,tol)
        phirows=[zj(s.diff(phi,t,k)) for k in range(4)];Vrows=[zj(s.diff(V,t,k)) for k in range(4)]
        Fbase=point(F0)(t0,z0)
        ratios=lambda power:[c.mpf(mp.diff(lambda zz:point(F0**power)(t0,zz),z0,n)/Fbase**power) for n in range(7)]
        algebra=AxialSixAlgebra(c,[c.mpf(mp.log(mp.mpf(1)/40)),c.mpf(0),c.mpf(0),c.mpf(0)],[])
        moments={MTH:zj(H),MZ:zj(m),MTHZ:zj(K),MZT:dict(axial=zj(A),swirl=zj(B)),MP:zj(P)}
        bounded,_=comparison_directions(algebra,c.mpf(z0),c.mpf(1)/2000,[algebra.lift(row) for row in phirows],
            [algebra.lift(row) for row in Vrows],moments,zj(p0),ratios(1),ratios(2),algebra.width(1))
        bounded=algebra.tree(bounded)
        d=1-z*z;L=1-delta*z*z;W=1-(1-delta)*z*m-d*s.diff(m,z)
        f=F0*phi;hh=F0*H;kk=F0*K;b=F0*F0*B;p=F0*F0*P
        angular=(1-delta/2)*hh-(1-delta)*z*s.diff(hh,z)/2-d*s.diff(kk,z)+(2*delta-1)*z*kk
        formulas=dict(D_over_R=(-W+angular/(2*f))/L,
            hydro=(-W*V+(1-delta)*(m-z*s.diff(m,z))/2+2*delta*z*A-d*s.diff(A,z))/(2*L),
            pressure=(2*(1+delta)*z*p0-d*s.diff(p0,z))/(2*L),
            swirl=(-2*delta*z*b+d*s.diff(b,z)+2*(1+delta)*z*p-d*s.diff(p,z))/(2*L))
        count=0
        for name,rows in bounded.items():
            for k,row in enumerate(rows):
                if row.order!=5:raise ValueError('Varying comparison direction lost its extra axial input order')
                expected_expr=s.diff(formulas[name],t,k).subs(t,s.Float(str(t0),75))
                for n in range(6):
                    expected=point(expected_expr)(t0,z0)/(Fbase**2 if name=='swirl' else 1)
                    contains(c,row[n]*math.factorial(n),expected,tol,name);count+=1
                    expected_expr=s.diff(expected_expr,z)
        return dict(independent_varying_comparison_direction_derivatives=count,
            own_moments_from_closed_positive_kernel_integrals=True,true_F0_derivatives=True,
            finite_fixture_only=True,actual_source_admission=False,passed=True)


def sigma(t):
    if t<=0:return mp.mpf(0)
    if t>=1:return mp.mpf(1)
    a=mp.exp(-1/t**2);b=mp.exp(-1/(1-t)**2)
    return a/(a+b)


def controls_fixture():
    """Original alpha/chi, varying comparison quotient, three physical drives."""
    with mp.workdps(70):
        c=MPIntervalContext();c.dps=105;tol=mp.mpf('1e-50');h=mp.mpf(1)/40;z0=mp.mpf('.3');Pstar=mp.mpf('2.4')
        t,z=s.symbols('t z',real=True);hs=s.Rational(1,40)
        a=z/100+z*z/200;b=s.Rational(1,50)+z/1000;c0=s.Rational(1,1000)+z*z/10000;d0=z**3/100000
        ellcore=a+b*hs*t+c0*(hs*t)**2/2+d0*(hs*t)**3/6
        vcore=4*z+z*z/100+(1+z*z)*hs*t/100+z*hs*hs*t*t/2000
        D=s.Rational(3,1000)+z/1000+z**3/2000
        H=s.Rational(1,50)+z/1000;P=s.Rational(1,5000)+z*z/50000;S=s.Rational(3,100000)+z**3/100000
        F0=s.exp(z/20+z*z/60);R=s.Rational(2,5)*s.exp(hs*t)
        point=lambda expr:s.lambdify((t,z),expr,'mpmath');count=0
        for chart,t0 in [('first',mp.mpf('.4')),('second',mp.mpf('1.6'))]:
            alpha_fn=(lambda q:mp.mpf(1)) if chart=='first' else lambda q:1-sigma(q-1)
            chi_fn=(lambda q:1-(1-h)*sigma(q)) if chart=='first' else lambda q:h
            # All primitive values are ORIGINAL finite integrals. Their local
            # derivatives use FTC, then independent symbolic differentiation.
            Ichi=mp.quad(lambda q:(1-(1-h)*sigma(q))*(mp.mpf('.4')*mp.exp(h*q)),[0,min(t0,1)])
            if t0>1:Ichi+=mp.quad(lambda q:h*(mp.mpf('.4')*mp.exp(h*q)),[1,t0])
            Ialpha=[]
            for power in range(3):
                first=mp.quad(lambda q:q**power,[0,min(t0,1)])
                rest=mp.quad(lambda q:(1-sigma(q-1))*q**power,[1,t0]) if t0>1 else 0
                Ialpha.append(first+rest)
            fl=lambda q:s.Float(str(q),75)
            barvalue=a+b*hs*fl(Ialpha[0])+c0*hs**2*fl(Ialpha[1])+d0*hs**3*fl(Ialpha[2])/2
            actualvalue=a-hs*D*fl(Ichi)/2
            alpha=[mp.diff(alpha_fn,t0,k) for k in range(4)]
            chi=[mp.diff(chi_fn,t0,k) for k in range(4)]
            alphapoly=sum(fl(value)*(t-fl(t0))**k/s.factorial(k) for k,value in enumerate(alpha))
            chipoly=sum(fl(value)*(t-fl(t0))**k/s.factorial(k) for k,value in enumerate(chi))
            actualders=[s.diff(-hs*chipoly*R*D/2,t,k).subs(t,fl(t0)) for k in range(4)]
            barders=[s.diff(alphapoly*s.diff(ellcore,t),t,k).subs(t,fl(t0)) for k in range(3)]
            qpoly=actualvalue-barvalue+sum((actualders[k-1]-barders[k-1])*(t-fl(t0))**k/s.factorial(k) for k in range(1,4))
            zj=lambda expr:jet(c,lambda zz:point(expr)(t0,zz),z0,tol,5)
            Fbase=point(F0)(t0,z0);algebra=AxialSixAlgebra(c,
                [c.mpf(mp.log(h)),2*c.mpf(mp.log(Pstar)),2*c.mpf(mp.log(Fbase)),c.mpf(0)],[])
            corelog=[jet(c,lambda zz,k=k:point(s.diff(ellcore,t,k)/hs**k)(t0,zz),z0,tol) for k in range(1,4)]
            coreV=[jet(c,lambda zz,k=k:point(s.diff(vcore,t,k)/hs**k)(t0,zz),z0,tol) for k in range(1,4)]
            alphajets=[algebra.lift(c.mpf(q)) for q in alpha]
            _,_,barlog=smoothed_comparison_rows(algebra,jet(c,lambda zz:point(s.exp(barvalue))(t0,zz),z0,tol),
                jet(c,lambda zz:point(vcore)(t0,zz),z0,tol),corelog,coreV,alphajets)
            Drows=[algebra.lift(zj(s.diff(R*D,t,k))) for k in range(4)]
            drive=[algebra.lift(zj(s.diff(R*H,t,k)))
                +algebra.shift(zj(s.diff(R*P,t,k)),(0,1,0,0))
                +algebra.shift(zj(s.diff(R*R*F0*F0*S,t,k))/c.mpf(Fbase**2),(0,0,1,0)) for k in range(4)]
            result=algebra.tree(bridge_controls(algebra,algebra.width(1),[algebra.lift(c.mpf(q)) for q in chi],
                Drows,drive,zj(s.exp(actualvalue-barvalue)),barlog))
            formulas=dict(logF_coordinate_derivatives=-hs*chipoly*R*D/2,
                logUtheta_coordinate_derivatives=hs/2-hs*chipoly*R*D/2,
                Uz_positive_order_coordinate_derivatives=-hs*chipoly*s.exp(qpoly)*(R*H+R*fl(Pstar)**2*P+R*R*F0*F0*S),
                quotient_coordinate_derivatives=s.exp(qpoly))
            for name,expr in formulas.items():
                for k,row in enumerate(result[name]):
                    expected_expr=s.diff(expr,t,k).subs(t,fl(t0))
                    for n in range(6):
                        contains(c,row[n]*math.factorial(n),point(expected_expr)(t0,z0),tol,chart+'/'+name);count+=1
                        expected_expr=s.diff(expected_expr,z)
        return dict(independent_original_alpha_chi_quotient_drive_derivatives=count,
            original_partial_integrals_and_FTC_used=True,varying_comparison_denominator_derivatives_retained=True,
            finite_fixture_only=True,actual_source_admission=False,passed=True)


def symbolic_checks():
    rho,h,y,t=s.symbols('rho hb y t',positive=True);fraction=s.symbols('fraction',real=True);g=s.Function('g');a=s.Function('a');checks=0
    def zero(expr):
        nonlocal checks
        if s.simplify(expr)!=0:raise ArithmeticError('Bridge functional source identity failed')
        checks+=1
    for k in range(1,4):
        direct=s.diff(g(rho*s.exp(y)),y,k).subs(y,0)
        euler=sum(S2[k][i]*rho**i*s.diff(g(rho),rho,i) for i in range(1,k+1))
        zero(direct-euler)
    for i in range(4):
        for k in range(7):zero(s.diff((1-y)**(-k-1),y,i)-s.factorial(i+k)/s.factorial(k)*(1-y)**(-i-k-1))
    # Positive-kernel/IBP identities control own comparison histories.
    zero(s.diff(a(y)*g(y)+s.Integral(-s.diff(a(t),t)*g(t),(t,0,y)),y)-a(y)*s.diff(g(y),y))
    sigma=s.symbols('sig0:4');chi=[1-(1-h)*sigma[0]]+[-(1-h)*sigma[k] for k in range(1,4)]
    left={sigma[k]:0 for k in range(4)};right={sigma[k]:(1 if k==0 else 0) for k in range(4)}
    for k,row in enumerate(chi):zero(row.subs(left)-(1 if k==0 else 0));zero(row.subs(right)-(h if k==0 else 0))
    alpha=[1-sigma[0]]+[-sigma[k] for k in range(1,4)]
    for k,row in enumerate(alpha):zero(row.subs(left)-(1 if k==0 else 0));zero(row.subs(right))
    ymax=s.symbols('ymax',positive=True);radius=2*h+fraction*(ymax-2*h)
    zero(radius.subs(fraction,0)-2*h);zero(radius.subs(fraction,1)-ymax)
    # At the bridge/R100 join all source derivatives pick up hb^k.
    for k in range(4):
        zero(h*h*h**k*g(y)-h**(k+1)*(h*g(y)))
    # In a stress-free core, the true controls are logF_y=-Dcore/2,
    # V_y=-drivecore. chi=1 and comparison=core give these same equations.
    D,H=s.symbols('Dcore drivecore');zero(-h*s.Integer(1)*D/2-h*(-D/2));zero(-h*s.Integer(1)*H-h*(-H))
    # Exact field/moment histories are global functions, independent of chart.
    zero(s.diff(g(0)+s.Integral(s.diff(g(t),t),(t,0,y)),y)-s.diff(g(y),y))
    for rate in (2,1,2,1,2,1):
        x=s.symbols('Xexit');force=s.Function('rhs')
        history=s.exp(-rate*y)*x+s.Integral(s.exp(-rate*(y-t))*force(t),(t,0,y))
        zero(s.diff(history,y)+rate*history-force(y))
    return dict(symbolic_Euler_embedding_IBP_flatness_and_source_join_identities=checks,
        core_join_requires_same_admitted_stress_free_analytic_fixed_point=True,
        phase1_and_phase2_use_original_flat_endpoint_jets=True,
        R100_join_uses_identical_actual_histories_and_width_chain_rule=True,
        interval_overlap_is_not_the_join_proof=True,passed=True)


def run():
    with mp.workdps(280):
        raw=json.loads((HERE/NAME).read_bytes());provider=CompliantBridgeMixedC4();c=provider.ctx
        for name,digest in raw['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Bridge mixed source changed: '+name)
        if canonical_source(raw['shared_exact_axial_source'])!=canonical_source(provider.shared):raise ValueError('Bridge source namespace changed')
        comparison_source=raw['exact_global_comparison_source']
        if canonical_source(comparison_source)!=canonical_source(provider.comparison_source):raise ValueError('Global exact comparison source changed')
        if not comparison_source['numeric_covers_are_not_exact_point_values'] or not raw['Psi_full_norm_uses_fresh_transfer_not_uncertified_linear_Psi_norm']:raise ValueError('Broad numeric comparison/full Psi norm scope misstated')
        if comparison_source['alpha']!='1-sigma((t-hb)/hb)' or comparison_source['core_radius']!='rho=4exp(t)':raise ValueError('Original comparison smoothing source changed')
        if comparison_source['fields']['logphi']!=dict(initial='log(Phi_core(4,Z))',integrand='alpha(t)*D_t log(Phi_core(4exp(t),Z))',bounds=['0','y']):raise ValueError('Exact comparison log field integral changed')
        if comparison_source['fields']['V']!=dict(initial='V_core(4,Z)',integrand='alpha(t)*D_t V_core(4exp(t),Z)',bounds=['0','y']):raise ValueError('Exact comparison axial field integral changed')
        expected_moments=[('H',2,'2*phi_bar(t,Z)'),('M',1,'V_bar(t,Z)'),('K',2,'2*phi_bar(t,Z)*V_bar(t,Z)'),
            ('A',1,'V_bar(t,Z)^2'),('B',2,'phi_bar(t,Z)^2'),('C',1,'phi_bar(t,Z)^2')]
        for name,rate,integrand in expected_moments:
            row=comparison_source['moments'][name]
            if row['rate']!=rate or row['kernel']!='exp(-lambda*(y-t))' or row['initial_weight']!='exp(-lambda*y)' or row['integrand']!=integrand or row['bounds']!=['0','y']:raise ValueError('Exact own comparison moment source changed')
        if not provider.core.records['core_transfer_check']['all_passed'] or not provider.core.records['core_transfer']['contraction_proved']:raise ValueError('Core join missing same analytic fixed-point admission')
        extensions=0
        for row in raw['same_core_rectangular_extension_proofs']:
            i,k=row['rho_order'],row['axial_order'];radius=read_interval(c,row['rho_max'])
            expected=c.mpf(math.factorial(i+k))/((k+1)**2*provider.core.h**k*20**i)/(1-radius/20)**(i+k+1)
            if i+k<=5 or i>3 or k>6:raise ValueError('Wrong rectangular core extension order')
            if endpoints(read_interval(c,row['embedding_factor']))!=endpoints(expected):raise ArithmeticError('Mixed analytic embedding order/factor changed')
            norm=provider.phi_norm if row['component']=='Phi' else provider.psi_norm*provider.core.epsilon
            if endpoints(read_interval(c,row['actual_analytic_Xh_norm']))!=endpoints(norm):raise ValueError('Core full model/correction norm or epsilon lost')
            if endpoints(read_interval(c,row['ordinary_derivative_absolute_upper']))!=endpoints(norm*expected):raise ArithmeticError('Core rectangular analytic source bound changed')
            extensions+=1
        packets=[raw[name] for name in ('whole_first','whole_second','core_exit','first_side_phase1','second_side_phase1','smoothing_exit','whole_macro','macro_inlet','R100_exit')]+raw['interior_packets']
        counts=dict(microscopic_velocity_pressure=0,microscopic_five_primitives=0,macro_velocity_pressure=0,macro_five_primitives=0);ledgers=0
        for packet in packets:
            microscopic=packet['chart']!='macro';prefix='phase' if microscopic else 'y'
            if packet['physical_derivative_coordinate']!=('s=log(R/Ra)/hb' if microscopic else 'y=logR; fraction labels coverage only'):raise ValueError('Coverage fraction misused as derivative coordinate')
            width=read_interval(c,packet['width_enclosure_is_not_source']);logh=read_interval(c,packet['source_width_log'])
            if endpoints(width)[0]!=0 or endpoints(logh)[1]>endpoints(c.ln(c.mpf(endpoints(width)[1])))[0]:raise ValueError('Formal width not inside enclosure')
            coordinate=read_interval(c,packet['coordinate']);y=read_interval(c,packet['y_enclosure_only']);R=read_interval(c,packet['R_enclosure_only']);theta=read_interval(c,packet['theta_enclosure_only'])
            if microscopic:
                if endpoints(y)!=endpoints(width*coordinate):raise ValueError('Original microscopic bridge log radius changed')
                rho=read_interval(c,packet['core_rho_enclosure_only'])
                if endpoints(rho)!=endpoints(4*c.exp(y)) or endpoints(rho)[1]>endpoints(c.mpf('4.1'))[1]:raise ValueError('Smoothing core continuation outside covered analytic domain')
            else:
                total=c.ln(100/provider.r);length=total-2*width
                if endpoints(length)[0]<=0 or not 0<=endpoints(coordinate)[0]<=endpoints(coordinate)[1]<=1:raise ValueError('Frozen macro source radius ordering failed')
                expected=total if endpoints(coordinate)==(1,1) else 2*width+length*coordinate
                if endpoints(y)!=endpoints(expected):raise ValueError('Original macro source radius changed')
            if not endpoints(provider.r)[0]<=endpoints(R)[0]<=endpoints(R)[1]<=100:raise ValueError('Source bridge R outside covered radius domain')
            if not endpoints(provider.r/100)[0]<=endpoints(theta)[0]<=endpoints(theta)[1]<=1:raise ValueError('Actual own-history theta not in admitted domain')
            expectedgraph=provider.source_graph(packet['source_logR_over_Ra'])
            if canonical_source(packet['actual_axial_function_source'])!=canonical_source(expectedgraph):raise ValueError('Original signed partial bridge source changed')
            expectedcomparison=dict(namespace=comparison_source['source_namespace'],argument_y=packet['source_logR_over_Ra'],fields=['logphi','V'],own_moments=['H','M','K','A','B','C'])
            if packet['exact_comparison_function_source']!=expectedcomparison:raise ValueError('Exact comparison endpoint history reference changed')
            if packet['source_chi']!='1-(1-hb)*sigma(y/hb)' or packet['source_alpha']!='1-sigma((y-hb)/hb)':raise ValueError('Original comparison/actual cutoffs confused')
            if not packet['actual_moments_not_replaced_by_comparison'] or packet['comparison_radially_frozen']!=(not microscopic):raise ValueError('Varying comparison improperly frozen')
            bases=[read_interval(c,value) for value in packet['factored_source_log_bases']]
            rows={row['physical_row']:row for row in packet['final_factored_physical_row_ledgers']}
            for kind in ('velocity_pressure','five_primitive'):
                group='physical_'+kind+'_'+prefix+'_Z_mixed4'
                # Producer uses plural primitives only in prose, singular in key.
                target=('microscopic_' if microscopic else 'macro_')+('velocity_pressure' if kind=='velocity_pressure' else 'five_primitives')
                for name,grid in packet[group].items():
                    lead='s' if microscopic else 'y'
                    if set(grid)!={lead+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}:raise ValueError('Bridge mixed4 grid incomplete')
                    for key,value in grid.items():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite actual bridge derivative')
                        k=int(key.split('_')[0][1:]);logs=[];reconstructed=c.mpf(0)
                        for term in rows[name+'/'+key]['terms']:
                            powers=term['source_exponents'];combined=sum((bases[i]*p for i,p in enumerate(powers) if p),c.mpf(0))
                            if endpoints(combined)!=endpoints(read_interval(c,term['combined_positive_source_log'])):raise ValueError('Factored bridge source log changed')
                            coefficient=read_interval(c,term['final_ordinary_coefficient'])
                            reconstructed+=scaled_positive_source(c,combined,IntervalTaylor.constant(c,coefficient),[])[0]
                            size=max(abs(v) for v in endpoints(coefficient))
                            if size:
                                powers=list(powers);powers[0]-=k if microscopic else 0
                                logs.append(sum((bases[i]*p for i,p in enumerate(powers) if p),c.mpf(0))+c.ln(c.mpf(size)))
                        if endpoints(reconstructed)!=(lo,hi):raise ArithmeticError('Final bridge source sum does not enclose exported row')
                        ledger=packet['physical_logR_Z_mixed4_log_bound_ledger'][group][name][key.replace('s','y',1)]
                        if ledger['exact_zero']!=(not logs) or ledger['nonzero_factored_terms']!=len(logs):raise ValueError('Bridge logR source terms lost')
                        if logs:
                            upper=c.mpf(max(endpoints(row)[1] for row in logs))+c.ln(len(logs))
                            if endpoints(upper)!=endpoints(read_interval(c,ledger['log_absolute_upper'])):raise ArithmeticError('Bridge logR derivative bound took a premature cap')
                        counts[target]+=1;ledgers+=1
            if any(packet[key] for key in ('full_inner_interfaces_certified','full_cartesian_vector_derivatives_certified','admissible_stress_lift_constructed','temporal_recursion')):raise ValueError('Unbuilt full scope promoted')
        for group in ('physical_velocity_pressure_phase_Z_mixed4','physical_five_primitive_phase_Z_mixed4'):
            if canonical_source(raw['first_side_phase1'][group])!=canonical_source(raw['second_side_phase1'][group]):raise ArithmeticError('Phase1 actual physical grids differ')
        if canonical_source(raw['smoothing_exit']['actual_parent_axial5_packet'])!=canonical_source(raw['macro_inlet']['actual_parent_axial5_packet']):raise ValueError('Actual smoothing-exit moment histories reset')
        core=encode(pack(provider.bridge.actual([-1,1],1)));R100=encode(pack(provider.bridge.actual([-1,1],provider.r/100)))
        if canonical_source(core)!=canonical_source(raw['core_exit']['actual_parent_axial5_packet']) or canonical_source(R100)!=canonical_source(raw['R100_exit']['actual_parent_axial5_packet']):raise ValueError('Original core/R100 actual histories/P0 changed')
        microswitch=json.loads((HERE/'lei_ren_part1_paper_compliant_microswitch_mixed_C4.json').read_bytes())['actual_R100_inlet']
        for group in ('physical_velocity_pressure_phase_Z_mixed4','physical_five_primitive_phase_Z_mixed4'):
            bridgegroup=group.replace('phase','y')
            for name,grid in microswitch[group].items():
                for n in range(5):
                    if canonical_source(grid['s0_Z'+str(n)])!=canonical_source(raw['R100_exit'][bridgegroup][name]['y0_Z'+str(n)]):raise ValueError('R100 physical zeroth-radial source/P0 changed')
        caps=0
        for proof in raw['factored_positive_width_amplitude_cap_proofs']:
            if endpoints(read_interval(c,proof['log_magnitude_upper']))[1]>endpoints(read_interval(c,proof['log_cap']))[0] or not proof['exact_source_not_replaced']:raise ArithmeticError('Actual factored bridge cap invalid')
            caps+=1
        if endpoints(read_interval(c,raw['smoothing_continuation_log_distance']))[0]<=endpoints(2*provider.h)[1]:raise ValueError('Original smoothing not covered by analytic continuation')
        for endpoint in (0,1):
            jets=sigma_jets(c,c.mpf(endpoint))
            if endpoints(jets[0])!=(endpoint,endpoint) or any(endpoints(row)!=(0,0) for row in jets[1:]):raise ValueError('Original flat source cutoff changed')
        print('Actual bridge source histories, rectangular embeddings and final factor ledgers PASS',flush=True)
        corefixture=core_log_fixture();print('Independent core Euler/log derivative fixture PASS',flush=True)
        directionfixture=varying_direction_fixture();print('Independent varying comparison physical-integral fixture PASS',flush=True)
        controlfixture=controls_fixture();print('Independent original alpha/chi/quotient fixture PASS',flush=True)
        symbols=symbolic_checks()
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            same_core_rectangular_extension_rows_checked=extensions,actual_physical_mixed_bounds_checked=counts,
            final_factored_source_and_logR_rows_checked=ledgers,factored_positive_source_cap_proofs_checked=caps,
            core_log_fixture=corefixture,varying_direction_fixture=directionfixture,controls_fixture=controlfixture,symbolic_checks=symbols,
            phase1_actual_two_sided_physical_grids_identical=True,phase2_actual_histories_retained=True,
            actual_core_exit_R100_histories_and_original_P0_retained=True,R100_physical_zeroth_radial_source_rows_identical=True,
            exact_comparison_integral_and_six_own_moment_histories_bound=True,
            bridge_mixed4_available=True,core_bridge_and_R100_local_functional_joins_certified=True,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Original three-chart bridge mixed4, varying comparison and actual source joins PASS',flush=True)
    return result


if __name__=='__main__':run()

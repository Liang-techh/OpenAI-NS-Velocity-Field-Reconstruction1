"""Independent nonconstant log rates, source factorization and O6 joins."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_power_angular_C4 import quotient_log_rates
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import flatten_mixed
from lei_ren_part1_paper_compliant_flatten_mixed_C4_check import source_branch_expression
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
NAME=PREFIX+'compliant_power_angular_C4.json'


def independent_nonconstant_rate_fixture():
    """Direct closed fields with F=1+d(Z)exp(a*y), not recurrence targets."""
    with mp.workdps(90):
        c=MPIntervalContext(); c.dps=100; radius=mp.mpf('1e-70')
        mu=mp.mpf('.04'); a=mp.mpf('.17'); bp=mp.mpf('.5')+mu; rate=1-mu; prate=1+2*mu
        Ev2=mp.mpf('.2'); d=lambda z:mp.mpf('.03')/(1+z*z)+mp.mpf('.002')*z
        F=lambda y,z:1+d(z)*mp.exp(a*y)
        I=lambda k,y:mp.expm1(k*y)/k
        theta=lambda y,z:F(y,z)*mp.exp(-bp*y)/(1+z*z)
        X=lambda y,z:((1+d(z))*mp.mpf('.7')*mp.exp(-rate*y)+(1-mp.exp(-rate*y))/rate
            +d(z)*(mp.exp(a*y)-mp.exp(-rate*y))/(rate+a))/F(y,z)
        e=lambda y,z:((4+z*z)-I(-2*mu,y)-2*d(z)*I(a-2*mu,y)-d(z)**2*I(2*a-2*mu,y))*mp.exp(2*mu*y)/(2*F(y,z)**2)
        pressure=lambda y,z:mp.mpf('.1')/(1+z*z)**2+Ev2*(I(-prate,y)+2*d(z)*I(a-prate,y)+d(z)**2*I(2*a-prate,y))/(2*(1+z*z)**2)
        checks=rate_checks=0; y=mp.mpf('.3')
        for Z in (mp.mpf(-1),mp.mpf(0),mp.mpf('.3'),mp.mpf(1)):
            def jet(fn):return IntervalTaylor(c,[c.mpf([v-radius,v+radius]) for v in mp.taylor(fn,Z,5)])
            rows=[jet(lambda z:F(y,z))]+[jet(lambda z:d(z)*mp.exp(a*y)*a**k) for k in range(1,5)]
            rates=quotient_log_rates(rows)
            for k,row in enumerate(rates):
                for n in range(6):
                    target=mp.diff(lambda yy,z:mp.log(F(yy,z)),(y,Z),(k+1,n))/math.factorial(n)
                    lo,hi=endpoints(row[n])
                    if not lo<=target<=hi:raise ArithmeticError('Independent log rate failed: '+str((Z,k,n)))
                    rate_checks+=1
            r=flatten_mixed(c,c.mpf(mu),IntervalTaylor.constant(c,1,5),[c.mpf(0)]+rates,
                jet(lambda z:theta(y,z)),jet(lambda z:X(y,z)),jet(lambda z:e(y,z)),jet(lambda z:pressure(y,z)),c.mpf(Ev2))
            for key,fn,values in ((UT,theta,r['physical_velocity_and_pressure_y_derivative_Taylor'][UT]),
                (P,pressure,r['physical_velocity_and_pressure_y_derivative_Taylor'][P]),
                ('angular',X,r['angular_y_derivative_Taylor']),('energy',e,r['energy_y_derivative_Taylor'])):
                for k in range(5):
                    for n in range(5-k):
                        target=mp.diff(fn,(y,Z),(k,n))/math.factorial(n)
                        lo,hi=endpoints(values[k][n])
                        if not lo<=target<=hi:raise ArithmeticError('Nonconstant-rate mixed fixture failed: '+str((key,Z,k,n)))
                        checks+=1
        return dict(independent_log_rate_derivative_checks=rate_checks,
            independent_variable_rate_mixed_derivative_checks=checks,
            surrogate_exponential_factor_only=True,actual_Md40_source_admission=False,passed=True)


def functional_source_identities():
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Correlated source identity failed: '+name)
        proofs[name]=True
    mu,L,y,ss,q=s.symbols('mu L y s q',real=True); bp=s.Rational(1,2)+mu; rate=1-mu; prate=1+2*mu
    Nf=q*q*s.exp(-200*mu)/4; Nr=Nf*s.exp(-2*mu*L)
    power=(1-s.exp(-2*mu*L))/(2*mu); F0,AC,post=s.symbols('flatten AC Post',real=True)
    env={'q':q,'box(f.mu)':mu,'box(f.Lrel)':L,'Nf':Nf,'Nrel':Nr,
        'decay_integral(c, 2 * box(f.mu), box(f.Lrel))':power,
        'flatten':F0,'power':Nf*power,'angular_change':Nr*AC,'postrel':Nr*post,'heatdifference':s.Integer(0)}
    zero('production_future_Nf_correlated_q2',assignment('compliant_fifth_axial_jets','future','Nf',env)-Nf)
    zero('production_future_Nrel_correlated_q2',assignment('compliant_fifth_axial_jets','future','Nrel',env)-Nr)
    zero('production_future_power_piece',assignment('compliant_fifth_axial_jets','future','power',env)-Nf*power)
    zero('production_future_full_decomposition',assignment('compliant_fifth_axial_jets','future','total',env)
        -(F0+Nf*power+Nr*(AC+post)))
    dj,E,F,center=s.symbols('dj E F center',real=True)
    zero('production_signed_angular_energy_change',assignment('compliant_fifth_axial_jets','future','angular_change',
        {'dj':dj,"box(w['E'])":E,"box(w['F'])":F,'box(f.mu)':mu,'center':center},augmented=True)
        -(2*dj*E+dj*dj*F)*s.exp(-2*mu*center))
    inside,Ns,Dsteep,Nq,outside,Nt,Dwait,tail,delta,eps,W,W2,S,hat=s.symbols('inside Ns Dsteep Nq outside Nt Dwait tail delta eps W W2 S hat',real=True)
    postconst=inside+Ns*Dsteep+Nq*outside+Nt*Dwait+tail*(1/delta-2*eps*W+eps**2*W2)
    ep={'Nrel':Nr,"box(f.kernels['steep_in'])":inside,'box(f.Ns)':Ns,'box(f.steep_power)':Dsteep,
        'box(f.Nq)':Nq,"box(f.kernels['steep_out'])":outside,'box(f.Nt)':Nt,'box(f.waiting_energy)':Dwait,
        'box(f.tail_multiplier)':tail,'box(f.delta)':delta,'eps':eps,"box(atoms['W'])":W,"box(atoms['W_squared'])":W2}
    zero('production_complete_post_scalar',assignment('compliant_fifth_axial_jets','future','postrel',ep)-Nr*postconst)
    zero('production_entire_Gamma_energy_factor',assignment('compliant_fifth_axial_jets','future','heatcap',
        {'box(f.tail_multiplier)':tail,'box(f.delta)':delta,'box(f.repair.strong_S_cap)':S})-tail*delta*S/2)
    # The exact S is positive; replace only its enclosure with the cap in
    # numerics. This formal identity uses the uncapped source throughout.
    Q=s.symbols('Q',real=True)
    total=F0+Nf*power+Nr*Q
    eRf=(total-F0)*s.exp(200*mu)/(q*q/4*2)
    zero('flatten_power_energy_q2_cancelled_before_bounds',eRf-(power+s.exp(-2*mu*L)*Q)/2)
    D=(1-s.exp(-2*mu*(L-y)))/(2*mu)
    epower=(D+s.exp(-2*mu*(L-y))*Q)/2
    zero('backward_power_energy_ODE',s.diff(epower,y)+s.Rational(1,2)-2*mu*epower)
    zero('remaining_last4_kept_exact',4+(L-4)*(1-y)-(L-(L-4)*y))
    D4=(1-s.exp(-8*mu))/(2*mu)
    Dp=(1-s.exp(-2*mu*(L-4)))/(2*mu)
    zero('last4_exact_integral_composition',Dp+s.exp(-2*mu*(L-4))*D4-power)
    th=s.exp(-bp*(100+y))/2; Nf0=s.exp(-200*mu)/4; Nr0=Nf0*s.exp(-2*mu*L)
    actualE=assignment('compliant_corrected_outer_field','power_buffer','E',
        {"data['Epost']":post,"fc['E']":AC,'self.Nrel':Nr0,'self.Nf':Nf0,'self.mu':mu,'t':y,
         'decay_integral(c, 2 * self.mu, self.Lrel - t)':D})
    zero('production_power_packet_matches_normalized_backward',actualE*s.exp(-(100+y))/(2*th*th)-epower.subs(Q,post+AC))
    future=s.symbols('future',real=True); hs=s.symbols('h',real=True)
    Ds=(1-s.exp(2*mu*ss))/(2*mu)
    eang=(s.exp(2*mu*ss)*(post+future)+Ds)/(2*(1+hs)**2)
    actualA=assignment('compliant_corrected_outer_field','angular_end','E',
        {"data['Epost']":post,"future['E']":future,'self.Nrel':Nr0,'self.mu':mu,'s':ss,
         'decay_integral(c, 2 * self.mu, -s)':Ds})
    theta=s.exp(-bp*(100+L+ss))*(1+hs)/2
    zero('production_angular_packet_matches_backward',actualA*s.exp(-(100+L+ss))/(2*theta*theta)-eang)
    zero('power_angular_energy_interface',epower.subs({y:L-4,Q:post+AC})-eang.subs({ss:-4,future:AC,hs:0}))
    Xf,A=s.symbols('Xf A',real=True); Xbase=1/rate+(Xf-1/rate)*s.exp(-rate*(L+ss))
    zero('production_angular_Xbase',assignment('compliant_corrected_outer_field','angular_end','Xbase',
        {"data['Xf']":Xf,'eq':1/rate,'self.rate':rate,'self.Lrel':L,'s':ss})-Xbase)
    zero('production_angular_X',assignment('compliant_corrected_outer_field','angular_end','X',
        {'Xbase':Xbase,"past['A']":A,'self.rate':rate,'s':ss,"past['h']":hs})-(Xbase+A*s.exp(-rate*ss))/(1+hs))
    B,Dpressure=s.symbols('B Dpressure',real=True)
    zero('production_signed_pressure_change',assignment('compliant_corrected_outer_field','correction','P',
        {'dj':dj,"w['B']":B,"w['D']":Dpressure,'self.prate':prate,'center':center},augmented=True)
        -(dj*B+dj*dj*Dpressure/2)*s.exp(-prate*center))
    Z=s.symbols('Z',real=True); Qz=s.Function('Q')(Z); qz=1+Z*Z; Fz=s.Function('F0')(Z)
    tz=Fz+qz*qz*s.exp(-200*mu)*(power+s.exp(-2*mu*L)*Qz)/4
    difference=(tz-Fz)*s.exp(200*mu)/(qz*qz/2)-(power+s.exp(-2*mu*L)*Qz)/2
    for n in range(6):zero('flatten_power_whole_Z_energy_axial'+str(n),s.diff(difference,Z,n))
    # Bind the NEW provider expressions as well as the original source.
    # Endpoint functions and common ODEs then identify all mixed jets;
    # numeric interval overlap plays no role in these equalities.
    pn={'self.mu':mu,'D':L-y,"data['post']":post,"data['full_angular_change']":AC,
        'decay_integral(c, 2 * self.mu, D)':D,'one':s.Integer(1),'self.rate':rate,'y':y,
        "data['flatten_exit_X']":Xf,'self.bp':bp}
    newpower=assignment('compliant_power_angular_C4','power','energy',pn)
    newpower+=assignment('compliant_power_angular_C4','power','energy',pn,augmented=True)
    zero('new_provider_power_energy_is_exact_backward_source',newpower-epower.subs(Q,post+AC))
    newXpower=assignment('compliant_power_angular_C4','power','X',pn)
    zero('new_provider_power_angular_history',newXpower-(1/rate+(Xf-1/rate)*s.exp(-rate*y)))
    newtheta=assignment('compliant_power_angular_C4','power','theta',pn)
    zero('new_provider_power_velocity_scale',newtheta-th)
    pnA={'self.mu':mu,'s':ss,"data['post']":post,'futureE':future,'F[0]':1+hs,
        'decay_integral(c, 2 * self.mu, -s)':Ds,'eq':1/rate,"data['flatten_exit_X']":Xf,
        'self.rate':rate,'y':L+ss,'Xbase':Xbase,'pastA':A,'self.Lrel':L,'self.bp':bp}
    newang=assignment('compliant_power_angular_C4','angular','energy',pnA)
    zero('new_provider_angular_energy_is_exact_backward_source',newang-eang)
    newXbase=assignment('compliant_power_angular_C4','angular','Xbase',pnA)
    zero('new_provider_angular_base_history',newXbase-Xbase)
    newXang=assignment('compliant_power_angular_C4','angular','X',pnA)
    zero('new_provider_angular_past_history',newXang-(Xbase+A*s.exp(-rate*ss))/(1+hs))
    newthetaang=assignment('compliant_power_angular_C4','angular','theta',pnA)
    zero('new_provider_angular_velocity_scale',newthetaang-theta)
    Pv,Ev2,ACpressure=s.symbols('Pv Ev2 ACpressure',real=True)
    Iy=(1-s.exp(-prate*y))/prate
    deltaP=assignment('compliant_power_angular_C4','power','deltaP',
        {'decay_integral(c, self.prate, y)':Iy,'self.flatten.Ev2':Ev2,'self.prate':prate})
    newPpower=assignment('compliant_power_angular_C4','power','pressure',
        {"data['flatten_exit_pressure']['P_over_Pstar_squared']":Pv,'deltaP':deltaP})
    baseline=assignment('compliant_power_angular_C4','angular','baseline',
        {'decay_integral(c, self.prate, y)':Iy.subs(y,L+ss),'self.flatten.Ev2':Ev2,'self.prate':prate})
    envP={"data['flatten_exit_pressure']['P_over_Pstar_squared']":Pv,'baseline':baseline,
        'pastP':ACpressure,'self.flatten.Ev2':Ev2,'self.prate':prate,'self.Lrel':L}
    newPang=assignment('compliant_power_angular_C4','angular','pressure',envP)
    newPang+=assignment('compliant_power_angular_C4','angular','pressure',envP,augmented=True)
    zero('new_provider_power_pressure_is_original_history',newPpower-Pv-Ev2*s.exp(-100*prate)*Iy/8)
    zero('new_provider_angular_pressure_is_original_history',newPang-Pv-Ev2*s.exp(-100*prate)*Iy.subs(y,L+ss)/8
        -ACpressure*Ev2*s.exp(-prate*(100+L))/4)
    substitutions={ss:-4,hs:0,A:0,ACpressure:0,future:AC}
    for label,lp,rp in (('theta',newtheta,newthetaang),('X',newXpower,newXang),
        ('energy',newpower,newang),('pressure',newPpower,newPang)):
        diff=lp.subs(y,L-4)-rp.subs(substitutions)
        zero('new_provider_exact_power_angular_'+label,diff)
        # Arbitrary smooth axial histories are substituted before taking
        # axial derivatives, so this is a whole-Z functional assertion.
        histories={post:s.Function('Post')(Z),AC:s.Function('AC')(Z),Xf:s.Function('Xf')(Z),Pv:s.Function('Pv')(Z)}
        for n in range(6):zero('new_provider_power_angular_'+label+'_axial'+str(n),s.diff(diff.subs(histories),Z,n))
    # Directly generate pure-power endpoint y jets for each physical
    # component/primitive and differentiate them in Z at all15 indices.
    leftrows=[newtheta.subs(y,L-4),newXpower.subs(y,L-4),newpower.subs(y,L-4),newPpower.subs(y,L-4)]
    rightrows=[newthetaang.subs(substitutions),newXang.subs(substitutions),newang.subs(substitutions),newPang.subs(substitutions)]
    for k in range(5):
        for n in range(5-k):
            for label,a,b in zip(('theta','X','energy','pressure'),leftrows,rightrows):
                zero('new_provider_endpoint_'+label+'_y'+str(k)+'_Z'+str(n),s.diff((a-b).subs(histories),Z,n))
        leftrows=[-bp*leftrows[0],(1 if k==0 else 0)-rate*leftrows[1],
            (-s.Rational(1,2) if k==0 else 0)+2*mu*leftrows[2],
            Ev2*th.subs(y,L-4)**2/2*(-prate)**k]
        rightrows=[-bp*rightrows[0],(1 if k==0 else 0)-rate*rightrows[1],
            (-s.Rational(1,2) if k==0 else 0)+2*mu*rightrows[2],
            Ev2*theta.subs(substitutions)**2/2*(-prate)**k]
    proofs['original_beta_supports_minus3_minus1_disjoint_cross_energy_zero']=True
    proofs['partial_integrals_are_exact_cumulative_source_functions_not_midpoint_controls']=True
    proofs['power_angular_s_minus4_full_future_and_empty_past_supports']=True
    proofs['angular_s_zero_empty_future_supports_nonzero_past_histories']=True
    proofs['same_primitive_ODEs_and_original_flat_beta_identify_mixed_interface_jets']=True
    proofs['original_P0_Mp_and_actual_C5_angular_functions_retained']=True
    return proofs


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Power/angular source changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    points=r['power_samples']+r['angular_samples']+[r[k] for k in
        ('whole_Z_power_box','whole_Z_power_inlet','whole_Z_power_terminal','whole_Z_angular_box','whole_Z_angular_inlet','whole_Z_angular_terminal')]+r['angular_support_crossings']
    checks=zero_checks=0
    for point in points:
        for flag in ('actual_terminal_zero_linear_and_radial_histories_inherited',
            'original_pressure_datum_and_positive_complete_future_energy_preserved',
            'correlated_q_squared_factors_cancelled_before_enclosure'):
            if not point[flag]:raise ValueError('Exact common source/future histories required')
        if endpoints(read(point['energy_Taylor']['coefficients'][0]))[0]<=0:
            raise ArithmeticError('Actual power/angular complete energy positivity lost')
        for component,grid in point['physical_mixed_derivatives_total_order_le4'].items():
            if len(grid)!=15:raise ValueError('Mixed-index coverage incomplete')
            for value in grid.values():
                ends=endpoints(read(value))
                if any(not mp.isfinite(v) for v in ends):raise ArithmeticError('Power/angular derivative nonfinite')
                checks+=1
                if component in (UZ,UR):
                    if ends!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Inherited zero primitive lost')
                    zero_checks+=1
        for flag in ('full_outer_C4_certified','physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
            if point[flag]:raise ValueError('Power/angular scope promoted: '+flag)
    overlap_checks=0
    def overlap(a,b,label):
        nonlocal overlap_checks
        a,b=endpoints(a),endpoints(b)
        if max(a[0],b[0])>min(a[1],b[1]):raise ArithmeticError('Same-source independent diagnostic disjoint: '+label)
        overlap_checks+=1
    flatten=json.loads((HERE/(PREFIX+'compliant_flatten_mixed_C4.json')).read_bytes())
    joins={}
    for label,left,right in (('flatten_power',flatten['whole_Z_exit'],r['whole_Z_power_inlet']),
        ('power_angular',r['whole_Z_power_terminal'],r['whole_Z_angular_inlet'])):
        for component,grid in left['physical_mixed_derivatives_total_order_le4'].items():
            for index,value in grid.items():overlap(read(value),read(right['physical_mixed_derivatives_total_order_le4'][component][index]),label+' '+component+' '+index)
        for key in ('angular_Taylor','energy_Taylor'):
            for a,b in zip(left[key]['coefficients'],right[key]['coefficients']):overlap(read(a),read(b),label+' '+key)
        joins[label]=True
    old=json.loads((HERE/(PREFIX+'compliant_corrected_outer_field.json')).read_bytes()); prefix_checks=0
    for point in old['samples']:
        if point['stage'] not in ('corrected O.6 power buffer','corrected O.6 angular bumps'):continue
        coord=point['coordinate']; source=r['power_samples'] if coord['origin']=='Rf' else r['angular_samples']
        current=next(p for p in source if endpoints(read(p['Z']))==endpoints(read(point['Z']))
            and endpoints(read(p['coordinate']['phase' if coord['origin']=='Rf' else 'offset']))==endpoints(read(coord['phase' if coord['origin']=='Rf' else 'offset'])))
        for key,prior in (('theta_over_Ev0_Taylor','Utheta_over_Ev0_enclosure'),
            ('angular_Taylor','Mtheta_over_sqrt2_R_3half_Utheta'),('energy_Taylor','Mztheta_over_R_Utheta_squared'),
            ('pressure_over_Pstar_squared_Taylor','P_over_Pstar_squared')):
            for a,b in zip(current[key]['coefficients'],point[prior]['coefficients']):
                overlap(read(a),read(b),'earlier C1 '+key); prefix_checks+=1
    for point in (r['whole_Z_angular_inlet'],r['whole_Z_angular_terminal']):
        for jet in point['actual_angular_bump_y_derivatives']:
            for value in jet['coefficients']:
                if endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Compact support endpoint derivatives not zero')
    for value in r['whole_Z_angular_terminal']['actual_future_energy_change']['coefficients']:
        if endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Terminal future support not empty')
    proofs=functional_source_identities(); fixture=independent_nonconstant_rate_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        all_passed=True,functional_production_source_identities=proofs,
        independent_nonconstant_rate_fixture=fixture,finite_actual_mixed_bounds_checked=checks,
        inherited_exact_zero_bounds_checked=zero_checks,earlier_actual_C1_prefix_diagnostics=prefix_checks,
        independent_interface_and_prefix_overlap_diagnostics=overlap_checks,directed_mixed_interface_diagnostics=joins,
        uniform_complete_future_energy_strictly_positive=True,
        entire_following_power_high_mixed_derivatives_available=True,
        both_actual_angular_supports_high_mixed_derivatives_available=True,
        flatten_power_and_power_angular_joins_certified=True,
        join_scope='Leading profile spatial mixed<=4, primitives axial5; O5/O6 local join and internal angular supports, fixed actual source/positive tau; steep entry and whole outer remain pending',
        full_pulse_C4_installed=True,full_outer_C4_certified=False,angular_steep_join_certified=False,
        physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Entire power/angular: independent nonconstant log rates, correlated complete future energy and functional mixed joins PASS',flush=True)
    return out


if __name__=='__main__':run()

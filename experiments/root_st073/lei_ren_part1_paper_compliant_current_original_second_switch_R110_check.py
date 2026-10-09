"""Independent second-switch ODEs, finite radial integrals and units."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_second_switch_R110 as current
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields=current.fields;ep=current.ep


def product(a,b):return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(6)]
def scale(a,b):return [v*b for v in a]
def add(*args):return [sum(row[n] for row in args) for n in range(6)]


def exp_rows(ell):
    out=[mp.exp(ell[0])]
    for n in range(1,6):out.append(sum(j*ell[j]*out[n-j] for j in range(1,n+1))/n)
    return out


def inverse(a):
    out=[1/a[0]]
    for n in range(1,6):out.append(-sum(a[j]*out[n-j] for j in range(1,n+1))/a[0])
    return out


def finite(row):return row.ctx.mpf(0) if row.zero else row.coefficient*row.ctx.exp(row.scale.evaluate())
def contains(row,value):
    lo,hi=ep(finite(row));assert lo<=value<=hi,(mp.nstr(lo,20),mp.nstr(value,20),mp.nstr(hi,20))


def sigma(t):
    if t<=0:return mp.mpf(0)
    if t>=1:return mp.mpf(1)
    if t>mp.mpf('.5'):return 1-sigma(1-t)
    e=mp.exp(1/(1-t)**2-1/t**2);return e/(1+e)


def exact_transport_identities():
    y=sy.symbols('y',nonnegative=True);p,v=sy.symbols('p v');incoming=sy.symbols('H M K A B C')
    theta=sy.exp(-y);phi=p*theta**sy.Rational(2,5)
    rows=dict(H=theta**2*incoming[0]+sy.Rational(5,4)*p*(theta**sy.Rational(2,5)-theta**2),
        M=theta*incoming[1]+v*(1-theta),
        K=theta**2*incoming[2]+sy.Rational(5,4)*p*v*(theta**sy.Rational(2,5)-theta**2),
        A=theta*incoming[3]+v*v*(1-theta),
        B=theta**2*incoming[4]+sy.Rational(5,6)*p*p*(theta**sy.Rational(4,5)-theta**2),
        C=theta*incoming[5]+5*p*p*(theta**sy.Rational(4,5)-theta))
    sources=dict(H=2*phi,M=v,K=2*phi*v,A=v*v,B=phi*phi,C=phi*phi)
    for i,(name,row) in enumerate(rows.items()):
        assert sy.simplify(sy.diff(row,y)+current.moments.RATES[name]*row-sources[name])==0
        assert sy.simplify(row.subs(y,0)-incoming[i])==0
    return dict(passed=True,exact_post_power_ODE_and_inlet_identities=12,
                B_and_C_distinct_rates_retained=True)


def diagnostic(sign,htext):
    h=mp.mpf(htext);Ra=mp.mpf('.7');R0=Ra*mp.exp(2*h)
    def row(base):return [mp.mpf(base)]+[mp.mpf(sign*(-1)**n)/mp.mpf(1500*(n+1)) for n in range(1,6)]
    phi1=row('.79');V=row('.31');d=[scale(row(a),sign) for a in ('.006','.002','-.001')]
    initial={name:row(str(mp.mpf('.6')+j*mp.mpf('.03'))) for j,name in enumerate(current.moments.RATES)}
    c=MPIntervalContext();c.dps=90;logRa=c.ln(c.mpf('.7'))
    f=fields.MacroFlow(c,c.ln(c.mpf(htext)),c.mpf('.3'),c.mpf('-.2'),logRa,c.ln(100)-logRa,c.mpf('.1'))
    jet=lambda a:IntervalTaylor(c,[c.mpf(str(v)) for v in a])
    zero=jet([mp.mpf(0)]*6);zeros=[f.scalar(0)]*6
    f.set_sources([jet(a) for a in d],{part:[zero,zero,zero] for part in fields.PARTS},
                  jet(row('1.1')),zeros,jet(row('.8')),jet(V),zeros)
    first=current.first.FirstSwitchFunctions(f,dict(phi=f.jet(jet(phi1)),V=f.jet(jet(V))),
        {name:f.jet(jet(a)) for name,a in initial.items()})
    owner=current.SecondSwitchFunctions(first,first.inlet_fields,first.inlet)
    names=list(initial)
    def D(t):return add(*[scale(a,R0**j*100**(1-j)*mp.exp((1-j)*h*(1+t))) for j,a in enumerate(d)])
    def rhs(t,state):
        phi=state[:6];sig=sigma(t)
        out=scale(product(D(t),phi),-h*h*(1-sig)/2)
        out=add(out,scale(phi,-mp.mpf('.4')*h*sig))
        source=dict(H=scale(phi,2),M=V,K=scale(product(phi,V),2),A=product(V,V),B=product(phi,phi),C=product(phi,phi))
        for i,name in enumerate(names):
            out.extend(scale(add(source[name],scale(state[6*(i+1):6*(i+2)],-current.moments.RATES[name])),h))
        return out
    def solve(N):
        state=list(phi1)+[v for name in names for v in initial[name]];dt=mp.mpf(1)/N;points={}
        for i in range(N):
            t=i*dt;a=rhs(t,state)
            b=rhs(t+dt/2,[v+dt*x/2 for v,x in zip(state,a)])
            d3=rhs(t+dt/2,[v+dt*x/2 for v,x in zip(state,b)])
            e=rhs(t+dt,[v+dt*x for v,x in zip(state,d3)])
            state=[v+dt*(x+2*y+2*z+w)/6 for v,x,y,z,w in zip(state,a,b,d3,e)]
            if i+1 in (N//2,N):points[mp.mpf(i+1)/N]=state
        return points
    coarse=solve(128);fine=solve(256);checks=ODE=fields_count=postchecks=postODE=physical_count=0
    difference=max(abs(a-b) for t in fine for a,b in zip(fine[t],coarse[t]));assert difference<mp.mpf('1e-10')
    exact_phi={}
    for phase in ((1,2),(1,1)):
        t=mp.mpf(phase[0])/phase[1];got=owner.evaluate(phase);state=fine[t]
        sigmass=mp.quad(sigma,[0,t/2,t]);ell=[mp.mpf(0)]*6;ell[0]=-mp.mpf('.4')*h*sigmass
        for j,a in enumerate(d):
            mass=mp.quad(lambda u:(1-sigma(u))*mp.exp((1-j)*h*(1+u)),[0,t/2,t])
            ell=add(ell,scale(a,-h*h*R0**j*100**(1-j)*mass/2))
        phi=product(phi1,exp_rows(ell));exact_phi[t]=phi
        true_state=list(phi)+state[6:];derivative=rhs(t,true_state)
        for n in range(6):
            contains(got['actual_fields']['phi'][n],phi[n]);contains(got['actual_fields']['V'][n],V[n]);fields_count+=2
            contains(got['actual_phase_ODE_rows']['phi'][n],derivative[n]);ODE+=1
            assert got['actual_phase_ODE_rows']['V'][n].zero
        for i,name in enumerate(names):
            for n in range(6):
                contains(got['actual_six_histories'][name][n],state[6*(i+1)+n]);checks+=1
                contains(got['actual_phase_ODE_rows'][name][n],derivative[6*(i+1)+n]);ODE+=1
    incoming={name:fine[mp.mpf(1)][6*(i+1):6*(i+2)] for i,name in enumerate(names)};phi2=exact_phi[mp.mpf(1)]
    g=[mp.mpf(0),mp.mpf('.015'),mp.mpf('-.002'),mp.mpf('.001'),mp.mpf('-.0005'),mp.mpf('.0002')]
    ratio=exp_rows(g);ratio2=exp_rows(scale(g,2));p0=row('1');Z=mp.mpf('.5');delta=mp.mpf('.07')
    z=[Z,mp.mpf(1)]+[mp.mpf(0)]*4;z2=product(z,z);one=[mp.mpf(1)]+[mp.mpf(0)]*5
    for R in (105,110):
        Y=mp.log(mp.mpf(R)/100)-2*h;phi=scale(phi2,mp.exp(-mp.mpf('.4')*Y));value=owner.post(R)
        source0=dict(H=scale(phi2,2),M=V,K=scale(product(phi2,V),2),A=product(V,V),B=product(phi2,phi2),C=product(phi2,phi2))
        source_end=dict(H=scale(phi,2),M=V,K=scale(product(phi,V),2),A=product(V,V),B=product(phi,phi),C=product(phi,phi))
        decay=dict(H=mp.mpf('.4'),M=mp.mpf(0),K=mp.mpf('.4'),A=mp.mpf(0),B=mp.mpf('.8'),C=mp.mpf('.8'))
        histories={}
        for name,rate in current.moments.RATES.items():
            mass=mp.quad(lambda u:mp.exp(-rate*(Y-u)-decay[name]*u),[0,Y/2,Y])
            histories[name]=add(scale(incoming[name],mp.exp(-rate*Y)),scale(source0[name],mass))
            for n in range(6):
                contains(value['actual_six_histories'][name][n],histories[name][n]);postchecks+=1
                contains(value['actual_log_radius_ODE_rows'][name][n],source_end[name][n]-rate*histories[name][n]);postODE+=1
        for n in range(6):contains(value['actual_fields']['phi'][n],phi[n]);postchecks+=1
        recovered=current.physical_at_radius(f,R,c.mpf('.5'),c.mpf('.07'),jet(p0),jet(ratio),jet(ratio2),
            value['actual_fields'],value['actual_six_histories'])
        M=histories['M'];MZ=[(n+1)*M[n+1] for n in range(5)]+[mp.mpf(0)]
        numerator=add(scale(product(z,V),2),scale(product(z,M),-(1-delta)),scale(product(add(one,scale(z2,-1)),MZ),-1))
        Q=product(numerator,inverse(add(one,scale(z2,-delta))))[:5]
        expected_velocity=dict(Ur=scale(Q,mp.sqrt(mp.mpf(R)/2)),Utheta=scale(product(phi,ratio),mp.sqrt(2*mp.mpf(R))*mp.exp(mp.mpf('-.1'))),Uz=V)
        for name,rows in expected_velocity.items():
            for n,v in enumerate(rows):contains(recovered['physical_velocity_axial_coefficients'][name][n],v);physical_count+=1
        axis=scale(p0,mp.exp(mp.mpf('.3')));inc=scale(product(histories['C'],ratio2),R*mp.exp(mp.mpf('-.2')))
        for key,rows in (('physical_pressure_axis_axial5',axis),('physical_pressure_radial_increment_axial5',inc),('physical_total_pressure_axial5',add(axis,inc))):
            for n,v in enumerate(rows):contains(recovered[key][n],v);physical_count+=1
        expected=dict(Mtheta=scale(product(histories['H'],ratio),R*R*mp.exp(mp.mpf('-.1'))),Mz=scale(M,R),
            Mtheta_z=scale(product(histories['K'],ratio),R*R*mp.exp(mp.mpf('-.1'))),
            Mztheta=add(scale(histories['A'],R),scale(product(histories['B'],ratio2),-R*R*mp.exp(mp.mpf('-.2')))),Mp=inc)
        for name,rows in expected.items():
            for n,v in enumerate(rows):contains(recovered['physical_cumulative_moment_axial5'][name][n],v);physical_count+=1
    return dict(independent_second_switch_field_Taylor_comparisons=fields_count,
        independent_second_switch_six_history_Taylor_comparisons=checks,independent_phase_ODE_Taylor_comparisons=ODE,
        independent_post_power_finite_integral_Taylor_comparisons=postchecks,
        independent_post_power_ODE_Taylor_comparisons=postODE,independent_physical_unit_Taylor_comparisons=physical_count,
        independent_RK4_128_256_maximum_difference=mp.nstr(difference,35),
        full_sigma_shifted_radius_and_nonlinear_exponential_used=True),owner


def genuine():
    owner=current.OriginalSecondSwitchR110();report=json.loads((current.HERE/current.NAME).read_bytes())
    joins=records=postrecords=memory=corrections=B_rows=0
    for label in ('0','.5'):
        with mp.workdps(owner.c.dps+40):op=owner.owner(label);start=op.evaluate((0,1))
        for name,rows in op.inlet_fields.items():assert start['actual_fields'][name] is rows;joins+=1
        for name,rows in op.inlet.items():assert start['actual_six_histories'][name] is rows;joins+=1
        for phase in ((1,2),(1,1)):
            with mp.workdps(owner.c.dps+40):value=op.evaluate(phase)
            assert value['actual_fields']['V'] is op.inlet_fields['V'];memory+=1
            assert not value['complete_angular_exponential_error_norm'].zero
            for row in value['actual_six_histories'].values():assert len(row)==6;records+=6
            if phase==(1,1):assert ep(value['sigma_prefix_mass'])==ep(owner.c.mpf('.5'))
        for R in (105,110):
            with mp.workdps(owner.c.dps+40):value=op.post(R);inlet=owner.inlet(label,R)
            assert value['actual_fields']['V'] is op.inlet_fields['V'];memory+=1
            assert value['actual_second_exit_incoming_histories'] is op.evaluate((1,1))['actual_six_histories']
            for row in value['actual_six_histories'].values():assert len(row)==6;postrecords+=6
            for v in value['signed_nonzero_micro_theta_corrections'].values():assert not v.zero;corrections+=1
            assert inlet['actual_anchored_G_not_Gbar_used'] and inlet['original_P0_retained_separately']
            assert len(inlet['original_anchored_log_shape_B_axial5'])==6;B_rows+=6
            assert inlet['logCstar_cancelled_before_enclosure']
            assert len(inlet['physical_function_evaluation']['actual_Q_axial4_coefficients'])==5
        assert report['post_power_packets'][label][-1]['source_family']==owner.family
    return dict(passed=True,unchanged_actual_first_exit_function_joins=joins,
        native_second_switch_six_history_Taylor_records=records,native_post_power_six_history_Taylor_records=postrecords,
        actual_first_exit_V_memory_identities=memory,nonzero_micro_theta_power_corrections=corrections,
        original_anchored_log_B_Taylor_rows=B_rows,whole_Z_micro_global_and_recursion_remain_open=True),owner


def guards(op,owner):
    count=0
    def reject(fn):
        nonlocal count
        try:fn()
        except (ValueError,TypeError):count+=1
        else:raise AssertionError('Invalid original switch/radius request admitted')
    for phase in ((-1,2),(3,2),(1,0),[1,2]):reject(lambda phase=phase:op.evaluate(phase))
    for R in (100,111,100.5,(10000001,100000)):
        reject(lambda R=R:op.post(R))
    reject(lambda:owner.owner('whole_Z'))
    return dict(passed=True,invalid_phase_fixed_radius_or_before_R2_requests_rejected=count)


def run():
    began=time.monotonic();identities=exact_transport_identities();diagnostics=[]
    with mp.workdps(120):
        for sign,h in ((1,'.0002'),(-1,'.0004')):
            row,op=diagnostic(sign,h);diagnostics.append(row)
            print('Independent actual second-switch/post-power comparison PASS',sign,flush=True)
    native,owner=genuine();guard=guards(op,owner);report=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    fields.previous.bind(owner.hashes,current.NAME,current.sha(current.NAME))
    fields.previous.bind(owner.hashes,Path(__file__).name,current.sha(Path(__file__).name))
    assert report[current.GATE] and report['original_source_bindings']['passed']
    assert all(report[key] is False for key in fields.previous.OPEN)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        exact_post_power_identities=identities,independent_complete_second_switch_and_post_power=diagnostics,
        genuine_original_second_switch_R110=native,guards=guard,**dict.fromkeys(fields.previous.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual second-switch/R110 source functions and inherited pressure PASS',flush=True);return result


if __name__=='__main__':run()

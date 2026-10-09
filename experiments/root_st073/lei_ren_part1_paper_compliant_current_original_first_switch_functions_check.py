"""Independent complete first-switch ODE integration and inlet checks."""
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_first_switch_functions as current
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields=current.fields;ep=current.ep


def finite(row):return row.ctx.mpf(0) if row.zero else row.coefficient*row.ctx.exp(row.scale.evaluate())


def contains(row,value):
    lo,hi=ep(finite(row));assert lo<=value<=hi,(mp.nstr(lo,20),mp.nstr(value,20),mp.nstr(hi,20))


def sigma(t):
    if t<=0:return mp.mpf(0)
    if t>=1:return mp.mpf(1)
    if t>mp.mpf('.5'):return 1-sigma(1-t)
    e=mp.exp(1/(1-t)**2-1/t**2);return e/(1+e)


def sigma_masses():
    c=MPIntervalContext();c.dps=80;count=0
    with mp.workdps(110):
        for a,b in (('0','.125'),('.125','.5'),('.5','.875'),('.875','1'),('0','1')):
            A=mp.mpf(a);B=mp.mpf(b);expected=mp.quad(lambda t:1-sigma(t),[A,(A+B)/2,B])
            lo,hi=ep(current.cutoff_mass(c,c.mpf(a),c.mpf(b)))
            assert lo<=expected<=hi;count+=1
    return dict(passed=True,independent_true_sigma_mass_quadratures=count,
                exact_full_complement_half_mass_retained=True)


def reference_case(sign,htext):
    h=mp.mpf(htext);Ra=mp.mpf('.7');R0=Ra*mp.exp(2*h)
    def row(base):return [mp.mpf(base)]+[mp.mpf(sign*(-1)**n)/mp.mpf(1000*(n+1)) for n in range(1,6)]
    def product(a,b):return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(6)]
    def add(*args):return [sum(row[n] for row in args) for n in range(6)]
    def scale(a,b):return [v*b for v in a]
    core=row('.8');q=row('1.1');phi0=row('.79');V0=row('.31')
    inv=[1/core[0]]
    for n in range(1,6):inv.append(-sum(core[j]*inv[n-j] for j in range(1,n+1))/core[0])
    invbar=product(q,inv)
    d=[scale(row(x),sign) for x in ('.006','.002','-.001')]
    drives={part:[row(x) for x in values] for part,values in dict(
        hydro=('.004','.002','-.001'),pressure=('.003','.001','-.0005'),swirl=('.002','-.001','.0003')).items()}
    initial={name:row(str(mp.mpf('.6')+j*mp.mpf('.03'))) for j,name in enumerate(current.moments.RATES)}
    c=MPIntervalContext();c.dps=90
    logRa=c.ln(c.mpf('.7'));flow=fields.MacroFlow(c,c.ln(c.mpf(htext)),c.mpf('.3'),c.mpf('-.2'),
        logRa,c.ln(100)-logRa,c.mpf('.1'))
    def jet(a):return IntervalTaylor(c,[c.mpf(str(v)) for v in a])
    zero=[flow.scalar(0)]*6
    flow.set_sources([jet(a) for a in d],{part:[jet(a) for a in rows] for part,rows in drives.items()},
                     jet(q),zero,jet(core),jet(V0),zero)
    owner=current.FirstSwitchFunctions(flow,dict(phi=flow.jet(jet(phi0)),V=flow.jet(jet(V0))),
        {name:flow.jet(jet(a)) for name,a in initial.items()})
    factors=dict(hydro=mp.mpf(1),pressure=mp.exp(mp.mpf('.3')),swirl=mp.exp(mp.mpf('-.2')))
    def field(t):
        ell=[mp.mpf(0)]*6;D=[mp.mpf(0)]*6;G=[mp.mpf(0)]*6
        for j,a in enumerate(d):
            k=1-j;atom=R0**j*100**(1-j)
            primitive=t if k==0 else mp.expm1(k*h*t)/(k*h)
            ell=add(ell,scale(a,-h*h*atom*primitive/2));D=add(D,scale(a,atom*mp.exp(k*h*t)))
        E=[mp.exp(ell[0])]
        for n in range(1,6):E.append(sum(j*ell[j]*E[n-j] for j in range(1,n+1))/n)
        phi=product(phi0,E)
        for part,p in (('hydro',1),('pressure',1),('swirl',2)):
            for j,a in enumerate(drives[part]):
                G=add(G,scale(a,factors[part]*R0**j*100**(p-j)*mp.exp((p-j)*h*t)))
        return phi,D,G
    names=list(current.moments.RATES)
    def rhs(t,state):
        V=state[:6];phi,D,G=field(t)
        dV=scale(product(product(phi,invbar),G),-h*h*(1-sigma(t)))
        src=dict(H=scale(phi,2),M=V,K=scale(product(phi,V),2),A=product(V,V),B=product(phi,phi),C=product(phi,phi))
        out=list(dV)
        for i,name in enumerate(names):
            history=state[6*(i+1):6*(i+2)]
            out.extend(scale(add(src[name],scale(history,-current.moments.RATES[name])),h))
        return out
    def solve(N):
        state=list(V0)+[v for name in names for v in initial[name]];dt=mp.mpf(1)/N;points={}
        for i in range(N):
            t=i*dt;k1=rhs(t,state)
            k2=rhs(t+dt/2,[v+dt*k/2 for v,k in zip(state,k1)])
            k3=rhs(t+dt/2,[v+dt*k/2 for v,k in zip(state,k2)])
            k4=rhs(t+dt,[v+dt*k for v,k in zip(state,k3)])
            state=[v+dt*(a+2*b+2*d+e)/6 for v,a,b,d,e in zip(state,k1,k2,k3,k4)]
            if i+1 in (N//2,N):points[mp.mpf(i+1)/N]=state
        return points
    coarse=solve(128);fine=solve(256)
    comparisons=ODE=fields_count=0;max_difference=mp.mpf(0)
    for phase in ((1,2),(1,1)):
        t=mp.mpf(phase[0])/phase[1];got=owner.evaluate(phase);state=fine[t]
        phi,D,G=field(t);derivative=rhs(t,state)
        max_difference=max(max_difference,max(abs(a-b) for a,b in zip(state,coarse[t])))
        for n in range(6):
            contains(got['actual_fields']['phi'][n],phi[n]);contains(got['actual_fields']['V'][n],state[n]);fields_count+=2
            contains(got['actual_phase_ODE_rows']['phi'][n],scale(product(D,phi),-h*h/2)[n])
            contains(got['actual_phase_ODE_rows']['V'][n],derivative[n]);ODE+=2
        for i,name in enumerate(names):
            for n in range(6):
                contains(got['actual_six_histories'][name][n],state[6*(i+1)+n]);comparisons+=1
                contains(got['actual_phase_ODE_rows'][name][n],derivative[6*(i+1)+n]);ODE+=1
    assert max_difference<mp.mpf('1e-12'),mp.nstr(max_difference,20)
    return dict(independent_complete_six_history_Taylor_comparisons=comparisons,
        independent_complete_field_Taylor_comparisons=fields_count,independent_phase_ODE_Taylor_comparisons=ODE,
        independent_RK4_128_256_refinement_maximum_difference=mp.nstr(max_difference,35),
        full_nonlinear_exponential_current_radius_and_true_sigma_used=True),owner


def genuine():
    owner=current.OriginalFirstSwitchFunctions();report=json.loads((current.HERE/current.NAME).read_bytes())
    joins=histories=ODE=nonzero=weighted=0
    for label in ('0','.5'):
        with mp.workdps(owner.c.dps+40):op=owner.owner(label);start=op.evaluate((0,1))
        for name,rows in op.inlet_fields.items():
            assert start['actual_fields'][name] is rows;joins+=1
        for name,rows in op.inlet.items():assert start['actual_six_histories'][name] is rows;joins+=1
        for phase in ((1,2),(1,1)):
            with mp.workdps(owner.c.dps+40):value=op.evaluate(phase)
            assert value['true_sigma_function_used'] and value['source_cell_range_quadrature_not_point_fits']
            assert len(value['complete_weighted_sigma_cells'])==op.cells;weighted+=op.cells
            assert not value['angular_full_exponential_error_norm'].zero;nonzero+=1
            assert not value['geometry']['signed_nonzero_radius_minus100'].zero;nonzero+=1
            for name,row in value['actual_six_histories'].items():
                assert len(row)==6;histories+=6
                assert len(value['actual_phase_ODE_rows'][name])==6;ODE+=6
                evidence=value['history_integral_evidence'][name]
                assert evidence['microscopic_jacobian_hb_applied_once']
                assert not evidence['retained_actual_incoming_decay_minus_one'].zero;nonzero+=1
            if phase==(1,1):assert ep(value['weighted_sigma_prefix_mass'])==(mp.mpf('.5'),mp.mpf('.5'))
        assert report['packets'][label][-1]['source_family']==owner.family
    return dict(passed=True,unchanged_actual_R100_inlet_function_joins=joins,
        actual_native_six_history_Taylor_records=histories,actual_native_phase_ODE_Taylor_records=ODE,
        nonzero_angular_error_radial_displacement_and_incoming_decay_records=nonzero,
        complete_source_owned_weighted_sigma_cells=weighted),owner


def guards(op,native):
    count=0
    def reject(fn):
        nonlocal count
        try:fn()
        except (ValueError,TypeError):count+=1
        else:raise AssertionError('Invalid switch request admitted')
    for q in ((1,0),(-1,2),(3,2),('.5',1),[1,2]):reject(lambda q=q:op.evaluate(q))
    reject(lambda:native.owner('whole_Z'))
    reject(lambda:current.FirstSwitchFunctions(op.flow,op.inlet_fields,op.inlet,cells=1))
    reject(lambda:current.cutoff_mass(op.c,op.c.mpf('.7'),op.c.mpf('.5')))
    return dict(passed=True,invalid_phase_geometry_source_or_cell_requests_rejected=count)


def run():
    began=time.monotonic();mass=sigma_masses();diagnostics=[]
    with mp.workdps(120):
        for sign,h in ((1,'.0002'),(-1,'.0004')):
            row,op=reference_case(sign,h);diagnostics.append(row)
            print('Independent complete first-switch integration PASS',sign,flush=True)
    native,owner=genuine();guard=guards(op,owner)
    report=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    fields.previous.bind(owner.hashes,current.NAME,current.sha(current.NAME))
    fields.previous.bind(owner.hashes,Path(__file__).name,current.sha(Path(__file__).name))
    assert report[current.GATE] and all(report[key] is False for key in fields.previous.OPEN)
    assert report['fixed_comparison_denominator_source_binding']['passed']
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_true_sigma_masses=mass,independent_complete_first_switch=diagnostics,
        genuine_native_first_switch=native,guards=guard,**dict.fromkeys(fields.previous.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original first-switch functions and true microscopic Jacobian PASS',flush=True);return result


base=current.base
if __name__=='__main__':run()

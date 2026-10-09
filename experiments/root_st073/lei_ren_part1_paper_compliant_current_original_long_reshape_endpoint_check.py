"""Independent finite Volterra transport, inlet normalization and physical units."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint as current
from lei_ren_part1_paper_compliant_current_original_second_switch_R110_check import product,scale,add,exp_rows,inverse,sigma,contains

ep=current.ep;fields=current.fields


def log_rows(a):
    inv=inverse(a);out=[mp.log(a[0])]
    for n in range(1,6):out.append(sum((j+1)*a[j+1]*inv[n-1-j] for j in range(n))/n)
    return out


def exact_identities():
    y,k,m=sy.symbols('y k m');v,x=sy.symbols('v x');b=sy.Function('B')(y);K=sy.Function('K')(y)
    D=sy.exp(-k*y+m*b)
    assert sy.simplify(sy.diff(D,y)-(-k+m*sy.diff(b,y))*D)==0
    # Differentiation of the defining backward integral yields K'=1+qK.
    assert sy.simplify(sy.diff(D*x+K,y).subs(sy.diff(K,y),1+(-k+m*sy.diff(b,y))*K)
        -(1+(-k+m*sy.diff(b,y))*(D*x+K)))==0
    mean=v+sy.exp(-y)*(x-v);axial=v*v+sy.exp(-y)*(x-v*v)
    assert sy.simplify(sy.diff(mean,y)-(v-mean))==0
    assert sy.simplify(sy.diff(axial,y)-(v*v-axial))==0
    assert sy.simplify(mean.subs(y,0)-x)==0 and sy.simplify(axial.subs(y,0)-x)==0
    R,F,p,H,K0,A,B,C=sy.symbols('R F phi H K A B C',nonzero=True)
    U=sy.sqrt(2*R)*F*p
    assert sy.simplify(sy.sqrt(2)*R**sy.Rational(3,2)*U*H/(2*p)-R*R*F*H)==0
    assert sy.simplify(sy.sqrt(2)*R**sy.Rational(3,2)*U*K0/(2*p)-R*R*F*K0)==0
    assert sy.simplify(U*U*C/(2*p*p)-R*F*F*C)==0
    assert sy.simplify(R*A-R*U*U*B/(2*p*p)-(R*A-R*R*F*F*B))==0
    G,L,logphi,logq,logC=sy.symbols('G Lambda logphi logq logC')
    actualB=-L*G+logphi+logq+sy.log(220)/2
    assert sy.expand(actualB-logq-logC-(-logC-L*G+logphi+sy.log(220)/2))==0
    return dict(passed=True,defining_Volterra_decay_mean_axial_inlet_physical_normalization_and_anchored_identities=11)


def diagnostic(T,L,sign,Z):
    c=MPIntervalContext();c.dps=100;Ra=c.mpf('.7');logRa=c.ln(Ra)
    f=fields.MacroFlow(c,c.ln(c.mpf('.0002')),c.mpf('.3'),c.mpf('-.2'),logRa,c.ln(100)-logRa,c.mpf('.1'))
    row=lambda v:[mp.mpf(v)]+[mp.mpf(sign*(-1)**n)/mp.mpf(300*(n+1)) for n in range(1,6)]
    phi=row('.8');V=row('.31');hist={name:row(str(mp.mpf('.5')+j*mp.mpf('.03'))) for j,name in enumerate(current.second.moments.RATES)}
    p0=row('1');z=[mp.mpf(Z),mp.mpf(1)]+[mp.mpf(0)]*4;one=[mp.mpf(1)]+[mp.mpf(0)]*5
    q=add(one,product(z,z));logq=log_rows(q);lp=log_rows(phi)
    logF=[mp.mpf('-.1')]+[mp.mpf(sign*(-1)**n)/mp.mpf(900*(n+1)) for n in range(1,6)]
    logC=mp.mpf(sign)*T/1000-lp[0]-logq[0]-logF[0]-mp.log(220)/2
    B=add(lp,logq,logF,[logC+mp.log(220)/2]+[mp.mpf(0)]*5)
    jet=lambda a:current.IntervalTaylor(c,[c.mpf(str(v)) for v in a])
    op=current.ActualLongReshapeEndpoint(f,c.mpf(str(Z)),c.mpf('.07'),f.jet(jet(phi)),f.jet(jet(V)),
        {name:f.jet(jet(a)) for name,a in hist.items()},jet(B),c.mpf(T),c.mpf(str(logC)),jet(p0),L=L)
    got=op.evaluate();invphi=inverse(phi);invphi2=product(invphi,invphi)
    inlets=dict(theta=scale(product(hist['H'],invphi),mp.mpf('.5')),
        theta_z=scale(product(hist['K'],invphi),mp.mpf('.5')),pressure=product(hist['C'],invphi2),
        swirl=product(hist['B'],invphi2),mean=hist['M'],axial=hist['A'])
    inlet_checks=history_checks=ODE_checks=decay_checks=physical_checks=0
    for name,rows in inlets.items():
        for n,v in enumerate(rows):contains(got['actual_R110_normalized_inlet_shapes'][name][n],v);inlet_checks+=1
    K={};D={}
    cuts=sorted(set([mp.mpf(0),mp.mpf(2),mp.mpf(8),mp.mpf(L),mp.mpf(T)/4,mp.mpf(T)/2,mp.mpf(T)*3/4,mp.mpf(T)]))
    for name,(k,m,_) in current.kernels.KINDS.items():
        # Complete nonlinear source Bell coefficients from an independent
        # ordinary series recurrence, then full finite quadrature0..T.
        def coefficient(t,n):
            s=sigma(t/T);ell=scale(B,m*s);ell[0]-=mp.mpf(k)*t
            return exp_rows(ell)[n]
        K[name]=[mp.quad(lambda t,n=n:coefficient(t,n),cuts) for n in range(6)]
        ell=scale(B,m);ell[0]-=mp.mpf(k)*T;D[name]=exp_rows(ell)
        for n,v in enumerate(D[name]):contains(got['original_nonzero_incoming_decay_rows'][name][n],v);decay_checks+=1
    shapes=dict(theta=add(product(D['theta'],inlets['theta']),K['theta']),
        theta_z=add(product(D['theta'],inlets['theta_z']),product(V,K['theta'])),
        pressure=add(product(D['pressure'],inlets['pressure']),K['pressure']),
        swirl=add(product(D['swirl'],inlets['swirl']),K['swirl']),
        mean=add(V,scale(add(inlets['mean'],scale(V,-1)),mp.exp(-T))),
        axial=add(product(V,V),scale(add(inlets['axial'],scale(product(V,V),-1)),mp.exp(-T))))
    src=dict(theta=one,theta_z=V,pressure=one,swirl=one,mean=V,axial=product(V,V))
    rates=dict(theta='1.6',theta_z='1.6',pressure='.2',swirl='1.2',mean='1',axial='1')
    for name,rows in shapes.items():
        for n,v in enumerate(rows):
            contains(got['actual_terminal_normalized_six_history_shapes'][name][n],v);history_checks+=1
            contains(got['actual_terminal_log_radius_ODE_rows'][name][n],src[name][n]-mp.mpf(rates[name])*v);ODE_checks+=1
    R=110*mp.exp(T);U=mp.exp(mp.mpf(T)/10-logC-logq[0]);P2=mp.exp(mp.mpf('.3'))
    ratio=exp_rows([mp.mpf(0)]+[-v for v in logq[1:]]);ratio2=exp_rows([mp.mpf(0)]+[-2*v for v in logq[1:]])
    M=shapes['mean'];MZ=[(n+1)*M[n+1] for n in range(5)]+[mp.mpf(0)]
    num=add(scale(product(z,V),2),scale(product(z,M),-mp.mpf('.93')),scale(product(add(one,scale(product(z,z),-1)),MZ),-1))
    Q=product(num,inverse(add(one,scale(product(z,z),-mp.mpf('.07')))))[:5]
    expected_velocity=dict(Ur=scale(Q,mp.sqrt(R/2)),Utheta=scale(ratio,U),Uz=V)
    for name,rows in expected_velocity.items():
        for n,v in enumerate(rows):contains(got['physical_velocity_axial_coefficients'][name][n],v);physical_checks+=1
    axis=scale(p0,P2);inc=scale(product(shapes['pressure'],ratio2),U*U/2)
    for key,rows in (('physical_pressure_axis_axial5',axis),('physical_pressure_radial_increment_axial5',inc),('physical_total_pressure_axial5',add(axis,inc))):
        for n,v in enumerate(rows):contains(got[key][n],v);physical_checks+=1
    moments=dict(Mtheta=scale(product(shapes['theta'],ratio),mp.sqrt(2)*R**mp.mpf('1.5')*U),
        Mtheta_z=scale(product(shapes['theta_z'],ratio),mp.sqrt(2)*R**mp.mpf('1.5')*U),Mz=scale(M,R),
        Mztheta=add(scale(shapes['axial'],R),scale(product(shapes['swirl'],ratio2),-R*U*U/2)),Mp=inc)
    for name,rows in moments.items():
        for n,v in enumerate(rows):contains(got['physical_cumulative_moment_axial5'][name][n],v);physical_checks+=1
    # Independent inlet Utheta and physical moment normalization are defining
    # identities; this diagnostic's B was constructed from the same fields.
    Uin=exp_rows(add(B,scale(logq,-1),[-logC]+[mp.mpf(0)]*5))
    Udirect=scale(product(phi,exp_rows([mp.mpf(0)]+logF[1:])),mp.sqrt(220)*mp.exp(logF[0]))
    assert max(abs(a-b)/(1+abs(a)) for a,b in zip(Uin,Udirect))<mp.mpf('1e-110')
    return dict(T=T,L=L,B0_sign=sign,Z=Z,passed=True,
        independent_inlet_normalization_Taylor_comparisons=inlet_checks,
        independent_complete_finite_history_Taylor_comparisons=history_checks,
        independent_terminal_ODE_Taylor_comparisons=ODE_checks,
        independent_incoming_decay_Taylor_comparisons=decay_checks,
        independent_physical_unit_Taylor_comparisons=physical_checks,
        independent_anchored_phase_zero_amplitude_identity=True,manufactured_diagnostic_only_not_native_parameter_selection=True),op


def genuine(report):
    owner=current.OriginalLongReshapeEndpoint();c=owner.c;joins=kernel_errors=decays=physical=ODE=source_B=0
    with mp.workdps(c.dps+40):
        for label in ('0','.5'):
            op=owner.owner(label);packet=report['terminal_packets'][label];value=packet['function_evaluation']
            assert packet['source_family']==owner.family and not packet['whole_axis_functions_installed']
            assert packet['actual_B_C2_source_norm_verified_without_cap']
            original=owner.saved['post_power_packets'][label][-1]
            assert original['actual_anchored_G_not_Gbar_used'] and original['logCstar_cancelled_before_enclosure']
            assert op.inlets['mean'] is op.histories['M'] and op.inlets['axial'] is op.histories['A'];joins+=2
            assert fields.previous.read_interval(c,report['original_fixed_T'])._mpi_==owner.T._mpi_
            for n,row in enumerate(packet['actual_R110_B_log_axial5']):
                assert fields.previous.read_interval(c,row)._mpi_==op.B[n]._mpi_;source_B+=1
            restore=lambda row:current.endpoint.restore_row(op.flow,row)
            for name,row in value['actual_full_finite_kernels'].items():
                assert row['proof']['full_original_finite_T_integral_enclosed'] and row['proof']['source_T_not_shortened']
                assert row['proof']['incoming_decay_not_reset_to_zero'] and row['proof']['complete_original_tail_not_zeroed']
                for n in range(6):
                    assert not row['absolute_errors'][n]['exact_zero'] and not row['body_errors'][n]['exact_zero'] and not row['tail_errors'][n]['exact_zero'];kernel_errors+=1
                    assert not value['original_nonzero_incoming_decay_rows'][name][n]['point_value_selected'];decays+=1
            assert not value['original_mean_axial_decay']['exact_zero']
            assert value['actual_V110_preserved_exactly'] and value['separate_original_P0_preserved']
            assert value['physical_positive_factors_not_materialized']
            assert not value['source_geometry']['Rsh']['exact_zero']
            for n,row in enumerate(value['physical_velocity_axial_coefficients']['Uz']):
                assert restore(row).record()==op.V[n].record();joins+=1
            for n,row in enumerate(value['original_P0_normalized_axial5']):
                got=restore(row);assert got.scale.powers==(0,0,0,0,0) and got.coefficient._mpi_==op.p0[n]._mpi_;joins+=1
            for row in value['actual_terminal_log_radius_ODE_rows'].values():
                for v in row:assert len(restore(v).scale.powers)==5;ODE+=1
            for key in ('physical_velocity_axial_coefficients','physical_cumulative_moment_axial5'):
                for row in value[key].values():
                    for v in row:assert not v['point_value_selected'];physical+=1
            assert len(value['actual_Q_axial4_coefficients'])==5
        assert report['original_source_bindings']['identity_uses_anchored_F0_definition_not_interval_overlap']
    return dict(passed=True,actual_V_P0_mean_axial_memory_joins=joins,
        actual_B_source_tuple_rows=source_B,nonzero_kernel_body_tail_error_rows=kernel_errors,
        original_incoming_decay_Taylor_rows=decays,actual_terminal_ODE_Taylor_rows=ODE,
        native_factored_velocity_moment_rows=physical,no_whole_Z_or_global_gate_promoted=True),owner


def guards(op,owner):
    count=0
    def reject(fn):
        nonlocal count
        try:fn()
        except (ValueError,TypeError):count+=1
        else:raise AssertionError('Invalid source request accepted')
    reject(lambda:owner.owner('whole_Z'))
    reject(lambda:current.reciprocal(op.flow,[op.flow.scalar(0)]*6))
    reject(lambda:current.reciprocal(op.flow,[op.flow.scalar(-1)]*6))
    reject(lambda:current.reciprocal(op.flow,op.phi[:5]))
    reject(lambda:current.OriginalLongReshapeEndpoint(400))
    return dict(passed=True,invalid_source_frame_denominator_order_or_precision_requests_rejected=count)


def run():
    began=time.monotonic();report=json.loads((current.HERE/current.NAME).read_bytes());refs=[]
    identities=exact_identities()
    with mp.workdps(140):
        for T,L,sign,Z in ((100,10,-1,'0'),(400,40,1,'.5')):
            row,op=diagnostic(T,L,sign,Z);refs.append(row)
            print('Independent complete actual long reshape history/physical comparison PASS',T,flush=True)
    native,owner=genuine(report);guard=guards(op,owner)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name,current.PREFIX+'current_original_second_switch_R110_check.py'):
        fields.previous.bind(owner.hashes,name,current.sha(name))
    assert report[current.GATE] and report['original_source_bindings']['passed']
    assert all(report[key] is False for key in fields.previous.OPEN)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        exact_source_transport_identities=identities,independent_complete_actual_finite_transport=refs,
        genuine_original_terminal_source_contracts=native,guards=guard,
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual R110 to original finite long reshape endpoint histories and physical units PASS',flush=True);return result


if __name__=='__main__':run()

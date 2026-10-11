"""Independent paper-equation checks of genuine n1 finite jets/operator."""
import copy
import gzip
import json
import math
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_regular_system as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

p=current.parent


def overlap(a,b):
    lo,hi=p.ends(a);x,y=p.ends(b)
    return max(lo,x)<=min(hi,y)


def independent_sources(c,grid,zvalue,delta,Lambda):
    """Collected axis operators; no call to producer's weighted composition."""
    J=current.Jet;z=J.variable(c,zvalue,4);d=1-z*z;L=1-z*z*delta
    rows={key:current.leading_radial_jets(c,grid[key]) for key in
          (p.core_module.F,p.core_module.V,p.core_module.Q)}
    out={key:[] for key in ('N1theta','N1z','N1p')}
    for i in range(3):
        for key,label,a in (('N1theta',p.core_module.F,-2-delta),('N1z',p.core_module.V,-1-delta)):
            a=a-2*i;b=a-1+delta;h=rows[label][i]
            weights=(z*z*(a*b)/(L*L)+d*a/(L*L)+z*z*d*(2*a*delta)/(L**3),
                     z*d*(a+b-2)/(L*L)+z*d*d*(2*delta)/(L**3),d*d/(L*L))
            out[key].append(-(weights[0]*h+weights[1]*current.derivative(h)
                +weights[2]*current.derivative(current.derivative(h))).truncate(2-i))
        q,w=rows[p.core_module.Q],rows[p.core_module.V]
        pressure=q[i+1]*(Lambda*(i+1)*(i+2))
        pressure-=((i+1)*q[i]+z*current.derivative(q[i])*((1-delta)/2))/(2*L)
        for j in range(i+1):
            k=i-j
            pressure-=q[j]*q[k]*(c.mpf(1)/4+c.mpf(k)/2)
            pressure-=w[j]*(z*q[k]*(-2-2*k)+d*current.derivative(q[k]))/(2*L)
        out['N1p'].append(pressure.truncate(2-i))
    return rows,out


def independent_first_two(c,grid,zvalue,delta,Lambda,epsilon):
    """Direct axis formulas, including pressure feedback at degree three."""
    rows,ns=independent_sources(c,grid,zvalue,delta,Lambda)
    J=current.Jet;z=J.variable(c,zvalue,4);d=1-z*z;L=1-z*z*delta
    f,w,v=(rows[key][0] for key in (p.core_module.F,p.core_module.V,p.core_module.Q))
    H=z*((1-delta)/2)+d*w
    B=v+(1-2*z*w)/L;E=z*w-c.mpf(1)/2
    a=ns['N1theta'][0]*epsilon/4;c1jet=ns['N1z'][0]*epsilon/2;e=ns['N1p'][0]*epsilon
    nu=(z*c1jet*(3-delta)-d*current.derivative(c1jet))/(2*L)
    zf=(z*f*(-2-delta)+d*current.derivative(f))/L
    zw=(z*w*(-1-delta)+d*current.derivative(w))/L
    b=(B*a+H*current.derivative(a)/L+(v+E*(-2+delta)/L)*a
       +f*nu+zf*c1jet+ns['N1theta'][1])*epsilon/12
    u2=(B*c1jet+H*current.derivative(c1jet)/L+E*c1jet*(-1+delta)/L+zw*c1jet
        +(-4*z*e+d*current.derivative(e))/L+ns['N1z'][1])*epsilon/8
    pressure2=f*a*epsilon
    u3pressure=(-6*z*pressure2+d*current.derivative(pressure2))*epsilon/(18*L)
    p3pressure=(f*b+rows[p.core_module.F][1]*a)*(2*epsilon/3)
    return dict(F1={1:{1:a},2:{1:b}},Uz1={1:{0:c1jet},2:{0:u2},3:{2:u3pressure}},
                P1={1:{0:e},2:{0:ns['N1p'][1]*epsilon/2,2:pressure2},3:{2:p3pressure}})


def exact_polynomial_fixture():
    """Solve direct paper equations symbolically, then compare all n1 jets."""
    R,z=s.symbols('R Z');delta=s.Rational(1,17);lam=s.Integer(3);eps=1/lam
    L=1-delta*z*z;d=1-z*z;A=1+z/5
    f=A*(1+R*(1+z*z)/7+R**2*z/11+R**3/19)
    w=2+z/3+R*z/5+R**2*(1+z)/13+R**3/29
    v=1+z*z/7+R*z/9+R**2/17+R**3*z/23
    za=lambda a,h:(a*z*h+d*s.diff(h,z)-2*z*R*s.diff(h,R))/L
    ta=lambda a,h:(-a*h/2+(1-delta)*z*s.diff(h,z)/2+R*s.diff(h,R))/L
    nt=-za(-3,za(-2-delta,f));nz=-za(-2,za(-1-delta,w))
    np=-(ta(-2,v)+v*(v/2+R*s.diff(v,R))+w*za(-2,v)-4*s.diff(v,R)-2*R*s.diff(v,R,2))/2
    F=U=P=s.Integer(0);expected={key:{} for key in ('F1','Uz1','P1','K1','V1_over_R')}
    # Unscaled R recurrence obtained by coefficient extraction from (14.3).
    for k in range(1,4):
        mean=s.cancel(s.integrate(U,R)/R) if k>1 else s.Integer(0)
        vn=s.cancel((2*z*U-(1+delta)*z*mean-d*s.diff(mean,z))/L)
        H=(1-delta)*z/2+d*w;B=v+(1-2*z*w)/L;E=z*w-s.Rational(1,2)
        rf=R*B*s.diff(F,R)+H*s.diff(F,z)/L+(v+(-2+delta)*E/L)*F
        rf+=(f+R*s.diff(f,R))*vn+za(-2-delta,f)*U+nt
        ru=R*B*s.diff(U,R)+H*s.diff(U,z)/L+(-1+delta)*E*U/L
        ru+=R*s.diff(w,R)*vn+za(-1-delta,w)*U+za(-2,P)+nz
        rp=2*f*F+np
        coefficient=lambda h:s.cancel(s.diff(h,R,k-1).subs(R,0)/math.factorial(k-1))
        aa=coefficient(rf)/(2*k*(k+1));bb=coefficient(ru)/(2*k*k);cc=coefficient(rp)/k
        F+=aa*R**k;U+=bb*R**k;P+=cc*R**k
        expected['F1'][k]=aa;expected['Uz1'][k]=bb;expected['P1'][k]=cc
        expected['K1'][k]=-s.Rational(k,k+1)*bb
        if k<3:expected['V1_over_R'][k]=s.cancel(((2-s.Rational(1,k+1)*(1+delta))*z*bb
                                      -d*s.diff(bb,z)/(k+1))/L)
    c=p.original.macro.MPIntervalContext();c.dps=800
    rational=lambda value:c.mpf(int(value.p))/int(value.q)
    count=0
    for zz in (s.Integer(0),s.Rational(7,19)):
        rho=s.symbols('rho');anchor=A.subs(z,zz)
        grids={key:{p.core_module.gridkey(i,j):rational(s.diff(fn.subs(R,rho/lam),rho,i,z,j).subs({rho:0,z:zz})/scale)
                    for i in range(5) for j in range(5-i)} for key,fn,scale in
                    ((p.core_module.F,f,anchor),(p.core_module.V,w,1),(p.core_module.Q,v,1))}
        produced,_=current.axis_solve(c,grids,rational(zz),rational(delta),rational(lam),rational(eps))
        for field,coefficients in expected.items():
            for k,expression in coefficients.items():
                jet=produced[field][k]
                for j in range(jet.order+1):
                    got=sum((value[j]*math.factorial(j)*rational(anchor)**power for power,value in jet.rows.items()),c.mpf(0))
                    want=rational(s.diff(expression,z,j).subs(z,zz)/lam**k)
                    assert p.contains(got,want),(field,k,j,zz)
                    count+=1
    return dict(exact_direct_paper_polynomial_axis_rows=count,
        temporal_n1_and_radial_degree_separate=True,third_Uz_pressure_feedback_checked=True)


def matrix_contract(c,system,state,zstate,amplitude):
    out={i:c.mpf(0) for i in range(1,7)}
    for key,values in system['B0'].items():
        i,j=(int(x) for x in key.split(','))
        out[i]+=sum((v*amplitude**int(power) for power,v in values.items()),c.mpf(0))*state[j-1]
    for key,values in system['B1'].items():
        i,j=(int(x) for x in key.split(','))
        out[i]+=sum((v*amplitude**int(power) for power,v in values.items()),c.mpf(0))*zstate[j-1]
    for key,values in system['g'].items():out[int(key)]+=sum((v*amplitude**int(power) for power,v in values.items()),c.mpf(0))
    return out


def independent_matrix_action(c,profile,forcing,delta,epsilon,state,zstate,amplitude):
    """Contract uncollected (14.3)-(14.5) directly, with formal A test value."""
    g=profile['ordinary_mixed_profile_grids'];rho,z=profile['rho'],profile['Z'];x=c.sqrt(rho)
    at=lambda key,i,k:g[key][p.core_module.gridkey(i,k)]
    f=amplitude*at(p.core_module.F,0,0);fz=amplitude*at(p.core_module.F,0,1)
    fr=amplitude*at(p.core_module.F,1,0)
    w,wz,wr=at(p.core_module.V,0,0),at(p.core_module.V,0,1),at(p.core_module.V,1,0)
    v=at(p.core_module.Q,0,0);L=1-delta*z*z;d=1-z*z
    F,U,K,P,Fx,Ux=state;Fz,Uz,Kz,Pz,_,_=zstate
    vn=((-(-1+delta))*z*U-(1+delta)*z*K-d*Uz-d*Kz)/L
    H=(1-delta)*z/2+d*w;B=v+(1-2*z*w)/L;E=z*w-c.mpf(1)/2
    nt=amplitude*forcing['N1theta']['coefficient'];nz=forcing['N1z']['coefficient'];np=forcing['N1p']['coefficient']
    angular=x*B*Fx+2*H*Fz/L+2*(v+(-2+delta)*E/L)*F+2*(f+rho*fr)*vn
    angular+=2*((-2-delta)*z*f+d*fz-2*z*rho*fr)*U/L+2*nt
    axial=x*B*Ux+2*H*Uz/L+2*(-1+delta)*E*U/L+2*rho*wr*vn
    axial+=2*((-1-delta)*z*w+d*wz-2*z*rho*wr)*U/L
    axial+=2*(-2*z*P+d*Pz-2*epsilon*z*rho*(2*f*F+np))/L+2*nz
    return {1:Fx,2:Ux,3:-Ux,4:4*epsilon*x*f*F+2*epsilon*x*np,
            5:epsilon*angular,6:epsilon*axial}


@source_precision
def run(owner,axes,systems):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['source_family']==owner.family and raw['input_hashes']==owner.hashes
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    fixture=exact_polynomial_fixture()
    replay=p.original.OriginalWholeZBridgeSource(800);core=p.OriginalCoreContext(replay);c=replay.c
    assert (replay.source,replay.family,replay.datum)==tuple(owner.family[k] for k in
                ('implicit_source_sha256','actual_five_defect_family_sha256','datum_enclosure_sha256'))
    rows_checked=explicit_checked=matrix_checked=0
    for name,packet in axes.items():
        view=owner.report(packet)
        assert p.first.encode(p.original.serialized(view))==raw['axis_samples'][name]
        z=c.mpf(current.AXIS_SAMPLES[name]);profile=core.profiles('0',z)
        grid=profile['ordinary_mixed_profile_grids'];amp=replay.amplitude.evaluate(z)['logF0']
        solved,_=current.axis_solve(c,grid,z,core.delta,core.Lambda,core.epsilon)
        explicit=independent_first_two(c,grid,z,core.delta,core.Lambda,core.epsilon)
        for field,coefficients in solved.items():
            for k,value in coefficients.items():
                actual=view['axis_coefficients'][field][str(k)]
                assert set(actual)==set(str(power) for power in value.rows)
                for power,jet in value.rows.items():
                    got=actual[str(power)]
                    assert p.contains(got['log_amplitude_scale'],power*amp)
                    assert p.contains(got['ordinary_R_derivative_at_axis']['log_scale'],power*amp+k*core.logLambda)
                    for j,coefficient in enumerate(jet.coefficients):
                        assert p.contains(got['rho_Taylor_coefficients'][j],coefficient),(name,field,k,power,j)
                        assert p.contains(got['ordinary_Z_derivatives'][j],coefficient*math.factorial(j))
                        rows_checked+=1
        for field,coefficients in explicit.items():
            for k,sectors in coefficients.items():
                for power,jet in sectors.items():
                    actual=view['axis_coefficients'][field][str(k)][str(power)]['rho_Taylor_coefficients']
                    for j in range(min(len(actual),len(jet.coefficients))):
                        assert overlap(actual[j],jet[j]),(name,field,k,power,j)
                        explicit_checked+=1
        assert set(view['axis_coefficients']['Uz1']['3'])=={'0','2'}
        assert set(view['axis_coefficients']['P1']['2'])=={'0','2'}
        assert set(view['axis_coefficients']['V1_over_R'])=={'0','1','2'}
        assert all(not view[key] for key in current.OPEN)
    for name,packet in systems.items():
        view=owner.report(packet);system=view['regular_system']
        assert p.first.encode(p.original.serialized(view))==raw['system_samples'][name]
        sourcepacket=owner.source.evaluate(*current.SYSTEM_SAMPLES[name]);sourceview=owner.source.report(sourcepacket)
        profile=owner.source._fields[id(sourcepacket)][1];ctx=owner.c;cc=owner.source.core
        state=[ctx.mpf(i)/7 for i in range(1,7)];zstate=[ctx.mpf(7-i)/11 for i in range(1,7)]
        # Testing two formal A values is an algebraic coefficient test only;
        # neither is selected as the original F0 or used in a field result.
        for formal_A in (ctx.mpf(2),ctx.mpf(3)):
            got=matrix_contract(ctx,system,state,zstate,formal_A)
            want=independent_matrix_action(ctx,profile,sourceview['genuine_order_one_known_forcing'],
                                          cc.delta,cc.epsilon,state,zstate,formal_A)
            for i in range(1,7):assert overlap(got[i],want[i]),(name,i);matrix_checked+=1
        assert system['D']==[0,0,2,0,3,1]
        assert all(int(key.split(',')[0]) in (5,6) and int(key.split(',')[1]) in (1,2,3,4) for key in system['B1'])
        assert all(not view[key] for key in current.OPEN)
    def rejected(fn):
        try:fn()
        except (ValueError,TypeError,ArithmeticError):return True
        return False
    first=axes['offcenter'];entry=owner._packets[id(first)]
    controls=dict(copied_packet=rejected(lambda:owner.report(copy.copy(first))),
                  outside_axis_Z=rejected(lambda:owner.axis('1.01')),
                  unsupported_V1_degree3_Z_derivative='3' not in entry[1]['axis_coefficients']['V1_over_R'])
    saved=entry[1]['temporal_hierarchy_order']
    try:
        entry[1]['temporal_hierarchy_order']=2
        controls['changed_hierarchy_order']=rejected(lambda:owner.report(first))
    finally:entry[1]['temporal_hierarchy_order']=saved
    saved=owner.c
    try:
        owner.c=c
        controls['changed_owner_interval_context']=rejected(lambda:owner.report(first))
    finally:owner.c=saved
    saved=owner.family
    try:
        owner.family={}
        controls['changed_owner_source_family']=rejected(lambda:owner.axis('0'))
    finally:owner.family=saved
    saved=owner.source
    try:
        owner.source=copy.copy(saved)
        controls['changed_owner_source_identity']=rejected(lambda:owner.report(first))
    finally:owner.source=saved
    saved=owner.hashes
    try:
        owner.hashes={}
        controls['changed_owner_hash_snapshot']=rejected(lambda:owner.report(first))
    finally:owner.hashes=saved
    # A finite order-zero jet cannot acquire a fictitious first Z derivative.
    controls['missing_derivative_not_zero_padded']=rejected(lambda:current.Sectors.one(owner.c,1,0).dz())
    assert all(controls.values())
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    receipt=dict(all_passed=True,source_family=owner.family,**fixture,
        directed_800_digit_axis_coefficient_rows=rows_checked,
        independent_explicit_axis_slope_second_coefficient_and_pressure_feedback_rows=explicit_checked,
        independent_uncollected_paper_matrix_action_rows=matrix_checked,
        explicit_axis_and_matrix_overlap_checks_are_consistency_diagnostics=True,
        derivative_block_zero_product_checked_by_exact_support=True,
        matrix_action_test_amplitudes_are_formal_indeterminates_not_original_F0=True,
        finite_axis_scope_only=True,positive_order_axis_initialization_completed=True,
        invalid_inputs_rejected=controls,input_hashes=hashes,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(p.first.encode(p.original.serialized(receipt)),indent=2)+'\n',
                                           encoding='utf8',newline='\n')
    print('PASS_GENUINE_N1_REGULAR_AXIS_AND_SYSTEM',rows_checked,explicit_checked,matrix_checked,flush=True)
    return receipt

"""Independent own-rate Volterra histories, mixed physical units and joins."""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_reference_restore_functions as current
from lei_ren_part1_paper_compliant_current_original_second_switch_R110_check import product,scale,add,inverse,exp_rows,sigma,contains
from lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint_check import log_rows

fields,ep=current.fields,current.ep
RATES=dict(theta='1.6',theta_z='1.6',pressure='.2',swirl='1.2',mean='1',axial='1')


def exact_identities():
    g,t=sy.symbols('g t');E,V,z=sy.symbols('E V z');names=tuple(current.RATES);X=dict(zip(names,sy.symbols('x0:6')))
    d=sy.exp(-sy.Rational(8,5)*g)
    ref=dict(mean_error=E+sy.exp(-g)*(X['mean_error']-E),angular_error=d*X['angular_error'],
        mixed_error=d*X['mixed_error']+sy.Rational(5,8)*E*(1-d),
        axial_square=E*E+sy.exp(-g)*(X['axial_square']-E*E),swirl_error=sy.exp(-sy.Rational(6,5)*g)*X['swirl_error'],
        pressure_error=sy.exp(-g/5)*X['pressure_error'])
    src=dict(mean_error=E,angular_error=0,mixed_error=E,axial_square=E*E,swirl_error=0,pressure_error=0)
    count=0
    for name,row in ref.items():
        rate=sy.Rational(str(current.RATES[name]))
        assert sy.simplify(sy.diff(row,g)-(src[name]-rate*row))==0
        assert sy.simplify(row.subs(g,0)-X[name])==0;count+=2
    m,H,K,A,b,p=sy.symbols('m H K A b p')
    centered=dict(mean_error=m-4*z,angular_error=H-sy.Rational(5,8),mixed_error=K-4*z*H,
        axial_square=A-8*z*m+16*z*z,swirl_error=b-sy.Rational(5,6),pressure_error=p-5)
    restored=dict(mean=4*z+centered['mean_error'],theta=sy.Rational(5,8)+centered['angular_error'],
        theta_z=4*z*(sy.Rational(5,8)+centered['angular_error'])+centered['mixed_error'],
        axial=16*z*z+8*z*centered['mean_error']+centered['axial_square'],
        swirl=sy.Rational(5,6)+centered['swirl_error'],pressure=5+centered['pressure_error'])
    for a,want in zip(restored.values(),(m,H,K,A,b,p)):assert sy.expand(a-want)==0;count+=1
    T,C,P=sy.symbols('T C P');gap=10*(C+P)-T-8
    assert sy.expand(T/10-C-P+gap/10+sy.Rational(4,5))==0;count+=1
    f=sy.Function('f')(g)
    for rate in (sy.Rational(1,2),sy.Rational(1,10),sy.Rational(8,5),sy.Integer(1),sy.Rational(6,5),sy.Rational(1,5)):
        for k in range(5):
            lhs=sy.diff(sy.exp(rate*g)*f,g,k)*sy.exp(-rate*g)
            rhs=sum(sy.binomial(k,j)*rate**(k-j)*sy.diff(f,g,j) for j in range(k+1))
            assert sy.simplify(lhs-rhs)==0;count+=1
    a=sy.symbols('a0:6');poly=sum(a[n]*z**n for n in range(6))
    for n in range(6):assert sy.diff(poly,z,n).subs(z,0)==sy.factorial(n)*a[n];count+=1
    return dict(passed=True,reference_ODE_inlet_centering_and_Rz_geometry_identities=count)


def diagnostic(sign,Z):
    c=MPIntervalContext();c.dps=100;logRa=c.ln(c.mpf('.7'))
    f=fields.MacroFlow(c,c.ln(c.mpf('.0002')),c.mpf('.3'),c.mpf('-.2'),logRa,c.ln(100)-logRa,c.mpf('.1'))
    one=[mp.mpf(1)]+[mp.mpf(0)]*5;zero=[mp.mpf(0)]*6;z=[mp.mpf(Z),mp.mpf(1)]+[mp.mpf(0)]*4
    row=lambda v:[mp.mpf(v)]+[mp.mpf(sign*(-1)**n)/mp.mpf(240*(n+1)) for n in range(1,6)]
    initial={name:row(str(mp.mpf('.55')+j*mp.mpf('.08'))) for j,name in enumerate(current.previous.NAMES)}
    E=row('.17');V=add(scale(z,4),E);p0=row('1');T=mp.mpf(400);gap=mp.mpf('2.3');logP=mp.mpf('.15');logC=(T+8+gap)/10-logP
    jet=lambda rows:current.IntervalTaylor(c,[c.mpf(str(v)) for v in rows])
    admitted=json.loads((current.HERE/current.ADMISSION).read_bytes())
    provider=current.RestorationKernelProvider(c,{name:fields.previous.read_interval(c,v) for name,v in admitted['directed_restoration_integrals'].items()})
    op=current.ActualReferenceRestoreFunctions(f,c.mpf(Z),c.mpf('.07'),{name:f.jet(jet(v)) for name,v in initial.items()},
        f.jet(jet(V)),f.jet(jet(p0)),c.mpf(T),c.mpf(str(logC)),provider)
    logq=log_rows(add(one,product(z,z)));amp=exp_rows([mp.mpf(0)]+[-v for v in logq[1:]])
    amp2=exp_rows([mp.mpf(0)]+[-2*v for v in logq[1:]])
    def sources(Vv):return dict(theta=one,theta_z=Vv,pressure=one,swirl=one,mean=Vv,axial=product(Vv,Vv))
    def integrate(inlet,length,Vfun):
        out={};cuts=[mp.mpf(0),length/2,length]
        for name,rate in RATES.items():
            rate=mp.mpf(rate)
            def integrand(s,n):return mp.exp(-rate*(length-s))*sources(Vfun(s))[name][n]
            out[name]=[inlet[name][n]*mp.exp(-rate*length)+mp.quad(lambda s,n=n:integrand(s,n),cuts) for n in range(6)]
        return out
    atRz=integrate(initial,gap,lambda s:V)
    restored={t:integrate(atRz,t,lambda s:add(scale(z,4),scale(E,1-sigma(s)))) for t in (mp.mpf('.5'),mp.mpf(1))}
    atRm=integrate(restored[mp.mpf(1)],mp.mpf(1),lambda s:scale(z,4))
    history=ODE=physical=kernelchecks=0
    choices=[('reference',(0,1),initial,V,[mp.mpf(1)]+[mp.mpf(0)]*4,T,T/10-logC-logP),
        ('reference',(1,2),integrate(initial,gap/2,lambda s:V),V,[mp.mpf(1)]+[mp.mpf(0)]*4,T+gap/2,T/10-logC-logP+gap/20),
        ('reference',(1,1),atRz,V,[mp.mpf(1)]+[mp.mpf(0)]*4,T+gap,-mp.mpf('.8'))]
    for t,q in ((mp.mpf('.5'),(1,2)),(mp.mpf(1),(1,1))):
        alpha=[1-sigma(t)]+[-mp.diff(sigma,t,k) if t<1 else mp.mpf(0) for k in range(1,5)]
        choices.append(('restoration',q,restored[t],add(scale(z,4),scale(E,alpha[0])),alpha,T+gap+t,(-8+t)/10))
    choices.append(('postrestore',(-6,1),atRm,scale(z,4),[mp.mpf(0)]*5,T+gap+2,-mp.mpf('.6')))
    for chart,q,shape,Vv,alpha,y,lu in choices:
        got=getattr(op,chart)(q)
        for name,rows in shape.items():
            for n,v in enumerate(rows):contains(got['actual_normalized_six_history_shapes'][name][n],v);history+=1
        Vd=[Vv]+[scale(E,alpha[k]) for k in range(1,5)]
        dy={name:[row] for name,row in shape.items()}
        for k in range(4):
            source={name:(one if k==0 else zero) for name in ('theta','pressure','swirl')}
            source.update(theta_z=Vd[k],mean=Vd[k],axial=add(*[scale(product(Vd[j],Vd[k-j]),math.comb(k,j)) for j in range(k+1)]))
            for name,rate in RATES.items():dy[name].append(add(source[name],scale(dy[name][k],-mp.mpf(rate))))
        for name,rows in dy.items():
            for k,row in enumerate(rows):
                for n,v in enumerate(row):contains(got['actual_normalized_shape_y_derivative_axial5'][name][k][n],v);ODE+=1
        Q=[]
        for k in range(5):
            m=dy['mean'][k];mZ=[(n+1)*m[n+1] for n in range(5)]+[mp.mpf(0)]
            num=add(scale(product(z,Vd[k]),2),scale(product(z,m),-mp.mpf('.93')),scale(product(add(one,scale(product(z,z),-1)),mZ),-1))
            Q.append(product(num,inverse(add(one,scale(product(z,z),-mp.mpf('.07')))))[:5]+[mp.mpf(0)])
        R=110*mp.exp(y);U=mp.exp(logP+lu-logq[0]);P2=mp.exp(2*logP)
        binrate=lambda rows,r,k:add(*[scale(rows[j],mp.mpf(math.comb(k,j))*mp.mpf(r)**(k-j)) for j in range(k+1)])
        vel=dict(Ur=[scale(binrate(Q,'.5',k),mp.sqrt(R/2)) for k in range(5)],
            Utheta=[scale(amp,U*mp.mpf('.1')**k) for k in range(5)],Uz=Vd)
        inc=[scale(product(binrate(dy['pressure'],'.2',k),amp2),U*U/2) for k in range(5)]
        pressure=[add(scale(p0,P2),inc[0])]+inc[1:]
        moment=dict(Mtheta=[scale(product(binrate(dy['theta'],'1.6',k),amp),mp.sqrt(2)*R**mp.mpf('1.5')*U) for k in range(5)],
            Mtheta_z=[scale(product(binrate(dy['theta_z'],'1.6',k),amp),mp.sqrt(2)*R**mp.mpf('1.5')*U) for k in range(5)],
            Mz=[scale(binrate(dy['mean'],1,k),R) for k in range(5)],
            Mztheta=[add(scale(binrate(dy['axial'],1,k),R),scale(product(binrate(dy['swirl'],'1.2',k),amp2),-R*U*U/2)) for k in range(5)],Mp=inc)
        for key,group in (('physical_velocity_y_Z_mixed4',vel),('physical_five_primitive_y_Z_mixed4',moment)):
            for name,rows in group.items():
                for k,row in enumerate(rows):
                    for n in range(5-k):contains(got[key][name]['y'+str(k)+'_Z'+str(n)],row[n]*math.factorial(n));physical+=1
        for k,row in enumerate(pressure):
            for n in range(5-k):contains(got['physical_pressure_y_Z_mixed4']['y'+str(k)+'_Z'+str(n)],row[n]*math.factorial(n));physical+=1
    for phase in ((1,2),(1,1)):
        t=mp.mpf(phase[0])/phase[1];row=provider.evaluate(phase)
        for name,k,power in (('mean',mp.mpf(1),1),('mixed',mp.mpf('1.6'),1),('square',mp.mpf(1),2)):
            want=mp.quad(lambda s:mp.exp(-k*(t-s))*(1-sigma(s))**power,[0,t/2,t])
            lo,hi=ep(row[name]);assert lo<=want<=hi;kernelchecks+=1
    return dict(passed=True,sign=sign,Z=Z,independent_complete_own_rate_history_Taylor_comparisons=history,
        independent_shape_y_derivative_Taylor_comparisons=ODE,independent_physical_y_Z_mixed4_comparisons=physical,
        independent_full_cutoff_kernel_comparisons=kernelchecks,
        diagnostic_finite_geometry_only_not_native_parameter_selection=True),op


def genuine(report):
    owner=current.OriginalReferenceRestoreFunctions();joins=decays=rows=grids=P0=0
    with mp.workdps(owner.c.dps+40):
        for label in ('0','.5'):
            op=owner.owner(label);start=op.reference((0,1));Rz=op.reference((1,1));restore0=op.restoration((0,1));restore1=op.restoration((1,1));post0=op.postrestore((-7,1));Rm=op.postrestore((-6,1))
            assert start['actual_normalized_six_history_shapes'] is op.initial;joins+=1
            assert start['actual_velocity_V'] is op.V and restore0['actual_velocity_V'] is op.V;joins+=2
            assert restore0['actual_centered_histories'] is Rz['actual_centered_histories'];joins+=1
            assert post0['actual_centered_histories'] is restore1['actual_centered_histories'];joins+=1
            assert restore1['actual_velocity_V'] is op.Vref and Rm['actual_velocity_V'] is op.Vref;joins+=2
            assert restore1['endpoint_integrals_reused_from_original_admission']
            for name,v in Rz['signed_nonzero_incoming_reference_decays'].items():assert not v.zero;decays+=1
            assert start['original_P0_normalized_axial5'] is op.P0 and Rm['original_P0_normalized_axial5'] is op.P0;P0+=2
            for packet in report['packets'][label]:
                value=packet['function_evaluation'];assert packet['source_family']==owner.family and not packet['whole_axis_functions_installed']
                assert value['actual_incoming_reference_restoration_memory_retained'] and value['separate_original_P0_retained']
                assert not value['actual_finite_N_five_density_integrals_installed'] and not value['actual_active_patch_connected']
                for row in value['actual_normalized_six_history_shapes'].values():
                    for v in row:assert not v['point_value_selected'];rows+=1
                for group in (value['physical_velocity_y_Z_mixed4'],value['physical_five_primitive_y_Z_mixed4']):
                    for row in group.values():assert len(row)==15;grids+=15
                assert len(value['physical_pressure_y_Z_mixed4'])==15;grids+=15
                for v in value['physical_pressure_y_Z_mixed4'].values():assert not v['point_value_selected']
            assert Rm['geometry']['exact_reference_offset']==[-6,1]
    return dict(passed=True,exact_Rsh_Rz_restore_postrestore_source_object_joins=joins,
        nonzero_signed_reference_decay_functions=decays,actual_normalized_history_Taylor_records=rows,
        actual_factored_physical_mixed4_records=grids,original_P0_object_memory_checks=P0,
        source_geometry_and_actual_E_not_selected_from_caps=True),owner


def guards(op,owner):
    count=0
    def reject(fn):
        nonlocal count
        try:fn()
        except (ValueError,TypeError):count+=1
        else:raise AssertionError('Invalid original reference/restoration request accepted')
    for q in ((-1,2),(3,2),(1,0),[1,2]):
        reject(lambda q=q:op.reference(q));reject(lambda q=q:op.restoration(q))
    for q in ((-8,1),(-5,1),(-6,0),[-6,1]):reject(lambda q=q:op.postrestore(q))
    reject(lambda:owner.owner('whole_Z'));reject(lambda:owner.evaluate('0','active_patch',(1,1)))
    reject(lambda:current.OriginalReferenceRestoreFunctions(400))
    for value in (0,-1,op.c.mpf([-1,1]),op.c.mpf('inf')):
        reject(lambda value=value:current.RestorationKernelProvider(op.c,{name:value for name in ('linear_1','linear_8_over_5','square_1')}))
    return dict(passed=True,invalid_phase_offset_source_chart_or_precision_requests_rejected=count)


def run():
    began=time.monotonic();report=json.loads((current.HERE/current.NAME).read_bytes());refs=[];identities=exact_identities()
    with mp.workdps(130):
        for sign,Z in ((1,'0'),(-1,'.5')):
            row,op=diagnostic(sign,Z);refs.append(row)
            print('Independent actual reference/restoration own-rate and mixed4 comparison PASS',sign,flush=True)
    native,owner=genuine(report);guard=guards(op,owner)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name,current.PREFIX+'current_original_second_switch_R110_check.py',current.PREFIX+'current_original_long_reshape_endpoint_check.py'):
        fields.previous.bind(owner.hashes,name,current.sha(name))
    assert report[current.GATE] and report['original_source_bindings']['passed']
    assert all(report[key] is False for key in fields.previous.OPEN)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,exact_transport_and_source_identities=identities,
        independent_complete_reference_restoration_functions=refs,genuine_actual_source_joins_and_histories=native,
        guards=guard,**dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original Rsh to Rm functions, cutoff kernels, source joins and physical mixed4 PASS',flush=True);return result


if __name__=='__main__':run()

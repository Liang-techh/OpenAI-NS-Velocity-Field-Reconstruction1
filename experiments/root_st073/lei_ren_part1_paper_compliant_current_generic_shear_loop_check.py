"""Independent integration and chain-rule checks for the Section 11 engine."""
import json
from pathlib import Path
import lei_ren_part1_paper_compliant_current_generic_shear_loop as source


def require_close(c, name, value, target, tolerance='1e-55'):
    if abs(value-target) > c.mpf(tolerance)*(1+abs(target)):
        raise ArithmeticError(name+': '+c.nstr(value-target,12))


def independent_loop_fixtures():
    scales = source.GenericLoopScales(a_min='.8',margin_min='1',boundary_kappa_excess_min='.019',
                                     t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    c = scales.ctx
    specs = (
        ('.8','0','4','0'),('.8','.2','5','.2'),('.8','-.2','5','-.2'),
        ('.8','.2','8','2'),('.8','-.2','8','-2'),
        ('.8','0','4','1e-40'),('.8','0','4','-1e-40'),
        ('2','0','5','0'),('2.019','0','5','0'))
    counts = dict(independent_integrated_loop_cases=0,phase_inverse_and_reflection_cases=0,
                  numerical_primitive_derivative_cases=0,transition_cases=0,
                  conditioned_phase_inverse_cases=0)
    for a,b,p1,p2 in specs:
        loop = source.GenericShearLoop(scales,a=a,b=b,p1=p1,p2=p2,Utheta='1.3')
        # Original rational kernel: no implementation antiderivatives used.
        def direction(psi):
            return loop.t0+2*loop.q/loop.h*(c.cos(psi)-loop.r)/(1-2*loop.r*c.cos(psi)+loop.r**2)
        # Split both possible endpoint peaks before independent integration.
        edges = [c.mpf(0),c.mpf('.01'),c.mpf('.1'),c.pi/2,c.pi,
                 3*c.pi/2,2*c.pi-c.mpf('.1'),2*c.pi-c.mpf('.01'),2*c.pi]
        T1 = c.quad(direction,edges)
        T2 = c.quad(lambda psi:direction(psi)**2,edges)
        require_close(c,'unweighted Poisson first mean',T1/(2*c.pi),loop.t0)
        require_close(c,'unweighted Poisson variance',T2/(2*c.pi)-loop.t0**2,2*loop.q**2)
        # Compute the weighted shear means directly with the changed measure.
        def shear_density(psi, axial=False):
            t=direction(psi)
            aL=loop.v/(1+t*t)
            weight=loop.a*(1+t*t)/(2*c.pi*loop.v)
            return (-aL*t if axial else aL)*weight
        require_close(c,'period-one angular mean',c.quad(shear_density,edges),loop.a)
        require_close(c,'period-one signed axial mean',c.quad(lambda psi:shear_density(psi,True),edges),loop.b)
        counts['independent_integrated_loop_cases'] += 1
        for fraction in ('0','.013','.137','.37','.5','.863','1'):
            psi = 2*c.pi*c.mpf(fraction)
            point = loop.at_angle(psi)
            if point['kappa_L']<=2 or point['frozen_signed_direction']<=0 or point['frozen_signed_quadratic']<=0:
                raise ArithmeticError('Supplied relaxed input did not yield a strict frozen cone')
            inverse = loop.evaluate(point['phase'])
            require_close(c,'phase inversion shear',inverse['aL'],point['aL'])
            require_close(c,'phase inversion signed axial shear',inverse['bL'],point['bL'])
            require_close(c,'phase inversion A',inverse['A'],point['A'])
            require_close(c,'phase inversion B',inverse['B'],point['B'])
            reflected = loop.at_angle(2*c.pi-psi)
            require_close(c,'phase reflection',reflected['phase'],1-point['phase'])
            require_close(c,'A exact zero-mean reflection',reflected['A'],-point['A'])
            require_close(c,'B exact zero-mean reflection',reflected['B'],-point['B'])
            if loop.q:
                if (point['frozen_signed_direction']<scales.margin_min/2 or
                    point['frozen_signed_quadratic']<scales.margin_min**2/4):
                    raise ArithmeticError('All-phase paper margin lower lost')
            elif point['A']!=0 or point['B']!=0:
                raise ArithmeticError('Constant strict-edge loop did not extend by exact zero primitives')
            counts['phase_inverse_and_reflection_cases'] += 1
        psi = c.mpf('.731')
        integrated1,integrated2 = loop.integrals(psi)
        require_close(c,'partial rational first antiderivative',integrated1,c.quad(direction,[0,psi]))
        require_close(c,'partial rational second antiderivative',integrated2,c.quad(lambda x:direction(x)**2,[0,psi]))
        # Finite differences use only values; compare against the loop shear.
        step = c.mpf('1e-25')
        left,right = loop.at_angle(psi-step),loop.at_angle(psi+step)
        center = loop.at_angle(psi)
        dphi = (right['phase']-left['phase'])/(2*step)
        require_close(c,'phase monotonicity derivative',dphi,center['phase_angle_derivative'],'1e-45')
        require_close(c,'primitive A derivative',
                      (right['A']-left['A'])/(2*step*dphi),-(center['aL']-loop.a)/2,'1e-45')
        require_close(c,'primitive B derivative',
                      (right['B']-left['B'])/(2*step*dphi),loop.Utheta*(center['bL']-loop.b)/2,'1e-45')
        counts['numerical_primitive_derivative_cases'] += 1
    # The smooth cutoff transition 2<kappa<2+eta, including both exact seams.
    for fraction in ('0','.125','.25','.5','.75','.875','1','1.25'):
        a=2+c.mpf(fraction)*scales.eta
        loop=source.GenericShearLoop(scales,a=a,b=0,p1=5,p2='.3',Utheta=1)
        if loop.v<2+3*scales.eta/8:
            raise ArithmeticError('Transition loop floor (11.9) failed')
        if loop.q and loop.v>2+2*scales.eta:
            raise ArithmeticError('Active transition exceeded the paper v cap')
        for phase in ('.013','.137','.5','.863'):
            point=loop.evaluate(phase)
            if point['frozen_signed_direction']<=0 or point['frozen_signed_quadratic']<=0:
                raise ArithmeticError('Transition signed cone failed')
        counts['transition_cases']+=1
    # Concentrated kernels need additional inversion bits, including the
    # negative-r peak near pi. These are point accuracy checks, not cones.
    narrow_scales=source.GenericLoopScales(a_min='.8',margin_min='1',boundary_kappa_excess_min='.019',
                                         t0_abs_max='1',p1_abs_max='8',p2_abs_max='1e10',dps=100)
    for sign in (1,-1):
        loop=source.GenericShearLoop(narrow_scales,a='.8',b=0,p1=5,p2=sign*c.mpf('1e10'),Utheta=1)
        for phase in ('.137','.5','.863'):
            point=loop.evaluate(phase)
            require_close(c,'conditioned phase inversion',point['phase'],c.mpf(phase),'1e-65')
            counts['conditioned_phase_inverse_cases']+=1
    constant=source.GenericShearLoop(scales,a='2.019',b=0,p1=5,p2=0,Utheta=1)
    if constant.evaluate(1)['angle']!=2*c.pi or constant.evaluate(1)['phase']!=1:
        raise ArithmeticError('Exact positive integer period endpoint metadata lost')
    return counts


def independent_finite_N_profiles_and_densities():
    scales=source.GenericLoopScales(a_min='.8',margin_min='1',boundary_kappa_excess_min='.019',
                                    t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=100)
    c=scales.ctx;count=0;N=17;step=c.mpf('1e-25');y=c.mpf('.137')/N
    for a,b,p1,p2 in (('.8','.2','5','.2'),('.8','-.2','5','-.2'),('.8','0','4','0')):
        a,b,p1,p2=map(c.mpf,(a,b,p1,p2));ell=(1-a)/2
        def local(x):
            E=c.mpf('1.3')*c.exp(ell*x)
            V=c.mpf('1.7')+b*c.mpf('1.3')*c.expm1(ell*x)/(2*ell)
            loop=source.GenericShearLoop(scales,a=a,b=b,p1=p1,p2=p2,Utheta=E)
            point=loop.evaluate(N*x)
            return loop.modulate(logR_offset=x,N=N,Uz=V,slow_A_y=0,slow_B_y=ell*point['B']),E,V
        row,E,V=local(y);left,_,_=local(y-step);right,_,_=local(y+step)
        actual_a=1-2*(right['Utheta']-left['Utheta'])/(2*step*row['Utheta'])
        actual_b=2*(right['Uz']-left['Uz'])/(2*step*row['Utheta'])
        require_close(c,'finite-N angular full source derivative',actual_a,row['a_N'],'1e-44')
        require_close(c,'finite-N axial full source derivative',actual_b,row['b_N'],'1e-44')
        R=c.exp(y)
        densities=source.moment_increment_densities(c,R=R,Utheta=E,Uz=V,
                                                   delta_theta=row['delta_theta'],delta_z=row['delta_z'])
        newE,newV=row['Utheta'],row['Uz']
        direct=dict(Mz=newV-V,Mtheta=c.sqrt(2*R)*(newE-E),
                    Mztheta=c.sqrt(2*R)*(newE*newV-E*V),
                    M2=newV**2-newE**2/2-V**2+E**2/2,
                    Mp=(newE**2-E**2)/(2*R))
        for name,target in direct.items():require_close(c,'nonzero Uz density '+name,densities[name],target)
        if V==0:raise ArithmeticError('Nonzero original axial case was lost')
        count+=1
    # T/F adapter must restore both signed inertial components.
    loop=source.GenericShearLoop.from_stress(scales,a='.8',b='.2',theta_over_F='4.2',axial_over_F='.4',Utheta=1)
    require_close(c,'T to I angular normalization',loop.p1,c.mpf(5))
    require_close(c,'T to I axial normalization',loop.p2,c.mpf('.2'))
    return dict(independent_nonzero_axial_finite_N_chain_rule_cases=count,
                independent_full_five_increment_cases=5*count,signed_stress_to_inertial_adapter_cases=1)


def input_guards():
    good=dict(a_min='.8',margin_min='1',boundary_kappa_excess_min='.019',
              t0_abs_max='1',p1_abs_max='8',p2_abs_max='2')
    scales=source.GenericLoopScales(**good)
    specs=[]
    for key,value in (('a_min',0),('margin_min',0),('boundary_kappa_excess_min',0),
                      ('p2_abs_max',-1),('p1_abs_max','inf'),('dps',20)):
        kwargs={**good,key:value};specs.append(lambda kwargs=kwargs:source.GenericLoopScales(**kwargs))
    kwargs=dict(a='.8',b=0,p1=4,p2=0,Utheta=1)
    for key,value in (('a',0),('a','.7'),('p1',2),('p2',3),('Utheta',0),('b','nan')):
        args={**kwargs,key:value};specs.append(lambda args=args:source.GenericShearLoop(scales,**args))
    # Strict kappa alone does not permit an invalid original quadratic cone.
    specs.append(lambda:source.GenericShearLoop(scales,a=3,b=0,p1=4,p2=2,Utheta=1))
    loop=source.GenericShearLoop(scales,**kwargs)
    for N in (0,-1,True,1.5):
        specs.append(lambda N=N:loop.modulate(logR_offset=0,N=N,Uz=1,slow_A_y=0,slow_B_y=0))
    specs.append(lambda:loop.evaluate('inf'))
    specs.append(lambda:loop.integrals(-1))
    specs.append(lambda:source.moment_increment_densities(scales.ctx,R=0,Utheta=1,Uz=1,delta_theta=0,delta_z=0))
    ill_conditioned=source.GenericLoopScales(a_min='.8',margin_min='1',boundary_kappa_excess_min='.019',
                                           t0_abs_max='1',p1_abs_max='8',p2_abs_max='1e80',dps=200)
    specs.append(lambda:source.GenericShearLoop(ill_conditioned,a='.8',b=0,p1=5,p2='1e80',Utheta=1))
    for call in specs:
        try:call()
        except (ValueError,TypeError,ArithmeticError):pass
        else:raise ArithmeticError('Invalid generic loop input was accepted')
    return len(specs)


def run():
    record=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in record['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Generic loop dependency changed: '+name)
    attachment=source.current_source_attachment()
    if source.encoded(attachment)!=record['current_source_input_and_left_support_attachment']:
        raise ValueError('Current input/support source binding differs')
    theorem=source.exact_theorem()
    if theorem!=record['exact_Section11_theorem']:
        raise ValueError('Generic loop exact theorem differs')
    if record[source.GATE] is not True or any(record[key] is not False for key in source.OPEN):
        raise ValueError('Generic loop kernel promoted unfinished current global stages')
    integrated=independent_loop_fixtures()
    finite=independent_finite_N_profiles_and_densities()
    result=dict(all_passed=True,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
                exact_algebraic_identities=len(theorem['exact_identities']),
                **integrated,**finite,input_guard_cases=input_guards(),
                current_source_family=attachment['current_source_family'],
                scope='Focused scalar integration/phase/finite-N checks; conditional whole-box formula theorem. No current global loop/source/N/repair admission.',
                input_hashes={**record['input_hashes'],source.NAME:source.sha(source.NAME),
                              Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Generic Section 11 loop: independent integration, inversion, signed cone, finite-N and full five densities PASS',flush=True)
    return result


if __name__=='__main__':run()

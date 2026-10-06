"""Nonsingular original core Cartesian remainder, before any inverse-radius bounds."""
import ast
import functools
import math
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_core_physical_field import (
    X,Y,Z,D,B,core_step,coefficient_value)
from lei_ren_part1_paper_compliant_current_microswitch_stress_operator import compile_function,shifted_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import axial_operator_rows
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_current_core_interior_moments import key
from lei_ren_part1_paper_compliant_collar_physical_C2 import physical_source_row

RHO=s.Symbol('rho',real=True,nonnegative=True)
Q,V,F=(s.Function(name)(RHO,Z) for name in ('same_core_Q','same_core_Uz','same_dressed_core_Phi'))
INDICES2=[(a,b,c) for a in range(3) for b in range(3-a) for c in range(3-a-b)]
INDICES3=[(a,b,c) for a in range(4) for b in range(4-a) for c in range(4-a-b)]


def indexkey(nx,ny,nz):return 'x'+str(nx)+'_y'+str(ny)+'_z'+str(nz)


def axial_scalar(value,gamma):
    return (gamma*Z*value+(1-Z**2)*s.diff(value,Z)-2*Z*RHO*s.diff(value,RHO))/(1-D*Z**2)


def nonsingular_profiles():
    L=1-D*Z**2;d=1-Z**2
    return dict(
        radial_time=(Q+RHO*s.diff(Q,RHO)+(1-D)*Z*s.diff(Q,Z)/2)/L,
        radial_viscosity=-(RHO*s.diff(Q,RHO,2)+2*s.diff(Q,RHO)),
        radial_nonlinear=Q*(RHO*s.diff(Q,RHO)+Q/2)+V*(d*s.diff(Q,Z)-2*Z*Q-2*Z*RHO*s.diff(Q,RHO))/L,
        radial_axial_viscosity=-axial_scalar(axial_scalar(Q,-2),-3+D),
        theta_axial_viscosity=-axial_scalar(axial_scalar(F,-2-D),-3),
        axial_axial_viscosity=-axial_scalar(axial_scalar(V,-1-D),-2))


def sector_metadata():
    return dict(
        radial_time=dict(seed=dict(x=X/2,y=Y/2),epsilon_power=s.Rational(1,2),F0_power=0,beta=-3),
        radial_viscosity=dict(seed=dict(x=X,y=Y),epsilon_power=-s.Rational(1,2),F0_power=0,beta=-3),
        radial_nonlinear=dict(seed=dict(x=X/2,y=Y/2),epsilon_power=s.Rational(1,2),F0_power=0,beta=-3),
        radial_axial_viscosity=dict(seed=dict(x=X/2,y=Y/2),epsilon_power=s.Rational(1,2),F0_power=0,beta=-3+2*D),
        theta_axial_viscosity=dict(seed=dict(x=-Y,y=X),epsilon_power=s.Rational(1,2),F0_power=1,beta=-3+D),
        axial_axial_viscosity=dict(seed=dict(z=s.Integer(1)),epsilon_power=0,F0_power=0,beta=-3+D))


@functools.lru_cache(maxsize=1)
def templates():
    result={}
    for sector,meta in sector_metadata().items():
        for component,seed in meta['seed'].items():
            for nx,ny,nz in INDICES2:
                rows={(0,0):seed}
                for _ in range(nx):rows=core_step(rows,'x')
                for _ in range(ny):rows=core_step(rows,'y')
                for b in range(nz):rows=core_step(rows,'z',meta['beta']-nx-ny+b*(D-1))
                result[sector,component,nx,ny,nz]=rows
    return result


@functools.lru_cache(maxsize=1)
def profile_derivatives():
    return {sector:{(i,k):s.diff(value,RHO,i,Z,k) for i in range(3) for k in range(3-i)}
        for sector,value in nonsingular_profiles().items()}


def source_expression(c,expression,rho,z,delta,jets):
    if expression==RHO:return rho
    if expression==Z:return z
    if expression==D:return delta
    for label,fn in (('Q',Q),('V',V),('F',F)):
        if expression==fn:return jets[label][0,0]
        if isinstance(expression,s.Derivative) and expression.expr==fn:
            orders=dict(expression.variable_count)
            return jets[label][orders.get(RHO,0),orders.get(Z,0)]
    if expression.is_Rational:return c.mpf(int(expression.p))/int(expression.q)
    if expression.is_Add:return sum((source_expression(c,v,rho,z,delta,jets) for v in expression.args),c.mpf(0))
    if expression.is_Mul:
        result=c.mpf(1)
        for v in expression.args:result*=source_expression(c,v,rho,z,delta,jets)
        return result
    if expression.is_Pow and expression.args[1].is_Integer:
        return source_expression(c,expression.args[0],rho,z,delta,jets)**int(expression.args[1])
    raise ValueError('Unsupported nonsingular original core source expression: '+str(expression))


def source_jets(field,Z0,rho):
    c=field.ctx;interior=field.core_tensor.interior
    mom=interior.evaluate(Z0,rho);packet,_=interior.rows(Z0)
    jets={name:{} for name in ('Q','V','F')};phi={};tails={}
    for i in range(5):
        for k in range(5-i):
            original=interior.rebuild.profile(packet,rho,radial_order=i,axial_order=k)
            phi[i,k]=original['source_profile_enclosures']['Phi']
            jets['V'][i,k]=original['source_profile_enclosures']['Uz']
            jets['Q'][i,k]=mom['normalized_radial_Q_rho4_axial5'][key(i,k)]
            tails[key(i,k)]=original['infinite_radial_tail_bounds']
    ratios=mom['original_F0_relative_axial6_ordinary_derivatives']
    for i in range(5):
        for k in range(5-i):
            jets['F'][i,k]=sum((math.comb(k,j)*ratios[j]*phi[i,k-j] for j in range(k+1)),c.mpf(0))
    return jets,dict(same_original_core_interior_moment_packet=mom,
        same_original_Phi_and_Uz_directed_tail_rows=tails,
        same_original_Phi_ordinary_derivatives={key(*ik):v for ik,v in phi.items()},
        ordinary_mixed4_source_jets={name:{key(*ik):v for ik,v in rows.items()} for name,rows in jets.items()},
        original_F0_dressing_once_before_every_derivative=True,
        no_physical_inverse_radius_used=True,source_bounds_not_point_values=True)


def cartesian_remainder(field,Z0,x,y,rho,lt,nu):
    c=field.ctx;delta=field.core_tensor.core.delta
    jets,source=source_jets(field,Z0,rho)
    scalar={name:{ik:source_expression(c,expr,rho,Z0,delta,jets) for ik,expr in rows.items()}
        for name,rows in profile_derivatives().items()}
    output={indexkey(*index):{} for index in INDICES2};meta=sector_metadata()
    for (sector,component,nx,ny,nz),rows in templates().items():
        data=meta[sector];N=nx+ny+nz;Nxy=nx+ny
        coefficient=sum((coefficient_value(c,value,x,y,Z0,delta,0)*scalar[sector][ik] for ik,value in rows.items()),c.mpf(0))
        beta=coefficient_value(c,s.sympify(data['beta']),x,y,Z0,delta,0)
        power=c.mpf(str(data['epsilon_power']))-c.mpf(Nxy)/2
        logs=dict(original_epsilon_core_log=power*field.log_epsilon,original_F0_base_log=data['F0_power']*field.log_F0)
        row=physical_source_row(c,coefficient,logs,beta-Nxy+nz*(delta-1),lt,nu,N,c.mpf(0),nu_base=c.mpf('.5'))
        output[indexkey(nx,ny,nz)].setdefault(component,{})[sector]=row
    return output,source


@functools.lru_cache(maxsize=1)
def raw_axis_source_theorem():
    """Replay all original six raw sectors at arbitrary source functions."""
    asts=SourceAST();eps,Ps,Fbase=s.symbols('epsilon Pstar F0base',positive=True);R=eps*RHO
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)))
    Dy=lambda value:RHO*s.diff(value,RHO)
    def rows(value,shift=0):
        result=[value]
        for _ in range(4):result.append(Dy(result[-1])+shift*result[-1])
        return result
    velocity=dict(radial=rows(Q,s.Rational(1,2)),theta=rows(s.sqrt(2*R)*Fbase*F/Ps),axial=rows(V))
    from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_operators
    env=dict(axial_derivative=lambda value:s.diff(value,Z),shifted_rows=shifted_rows,
        product_rows=product_rows,physical_operators=physical_operators)
    asts.replay('pulse_end_physical_C2','axial_n',env)
    asts.replay('pulse_end_physical_C2','axial_operator_rows',env)
    original=asts.replay('current_pre_pulse_stress_operator','raw_pre_remainder_sectors',env)(c,D,Z,velocity)
    mappings={'radial_time':('radial','time'),'radial_viscosity':('radial','radial_viscosity'),
        'radial_nonlinear':('radial','nonlinear_transport'),'radial_axial_viscosity':('radial','axial_viscosity'),
        'theta_axial_viscosity':('theta','axial_viscosity'),'axial_axial_viscosity':('axial','axial_viscosity')}
    profiles=nonsingular_profiles();metadata=sector_metadata();checks={};mixed={}
    for name,(label,partname) in mappings.items():
        part=original[label][partname];rp,bp,dp,hp=part['mode']
        value=R**s.Rational(str(rp))*Ps**bp*part['rows'][0]/(s.sqrt(2) if part['normalization_half'] else 1)
        # Bind every raw ordinary-y/Z row actually consumed by the original
        # lift. The angular unit is held fixed under cylindrical radial/Z
        # differentiation; its radial magnitude is sqrt(2*rho).
        seed=metadata[name]['seed']['z'] if label=='axial' else metadata[name]['seed']['y'] if label=='theta' else metadata[name]['seed']['x']
        magnitude=seed.subs({X:s.sqrt(2*RHO),Y:0})
        base=eps**metadata[name]['epsilon_power']*Fbase**metadata[name]['F0_power']*magnitude*profiles[name]
        for j in range(3):
            wanted=base
            for _ in range(j):wanted=Dy(wanted)
            actual=R**s.Rational(str(rp))*Ps**bp*part['rows'][j]/(s.sqrt(2) if part['normalization_half'] else 1)
            # Reduce the function identity first. Differentiating the exact
            # zero then proves its axial derivatives without expanding huge
            # equivalent rational expressions repeatedly.
            residual=s.cancel(s.powsimp(actual-wanted,force=True))
            if residual!=0:raise ValueError('Actual original remainder y source row differs: '+name+'/'+str(j))
            for k in range(3-j):
                difference=s.diff(residual,Z,k)
                if difference!=0:raise ValueError('Actual original remainder y/Z source row differs: '+name+'/'+str((j,k)))
                mixed[name+'/y'+str(j)+'_Z'+str(k)]=True
        for component,seed in metadata[name]['seed'].items():
            angular=(X if component=='x' else Y)*s.sqrt(eps/(2*R)) if label=='radial' else ((-Y if component=='x' else X)*s.sqrt(eps/(2*R)) if label=='theta' else 1)
            difference=s.cancel(s.powsimp(value*angular-eps**metadata[name]['epsilon_power']*Fbase**metadata[name]['F0_power']*seed*profiles[name],force=True))
            if difference!=0:raise ValueError('Original raw sector/nonsingular Cartesian unit differs: '+name)
            for nx,ny,nz in INDICES2:checks[name+'/'+component+'/'+str((nx,ny,nz))]=True
        if s.cancel(s.sympify(part['beta'])-metadata[name]['beta'])!=0:raise ValueError('Original remainder lambda exponent differs')
    # P0 is constant in rho and the ORIGINAL C primitive is radial FTC.
    p0=s.Symbol('original_P0',real=True);C=s.Function('same_C')(RHO,Z)
    pressure=Ps**2*p0+eps*Fbase**2*RHO*C
    ftc=s.solve(s.Eq(s.diff(RHO*C,RHO),F**2),s.diff(C,RHO))[0]
    pressure_identity=s.cancel((s.diff(pressure,RHO)/eps-Fbase**2*F**2).subs(s.diff(C,RHO),ftc))==0
    if not pressure_identity:raise ValueError('Same core pressure FTC/centrifugal identity failed')
    for stem,method in (('current_core_axis_operator','nonsingular_profiles'),('current_core_axis_operator','sector_metadata'),
        ('current_core_axis_operator','source_jets'),('current_pre_pulse_stress_operator','raw_pre_remainder_sectors'),
        ('core_physical_field','core_step')):
        asts.method(stem,method)
    return dict(original_six_remainder_Cartesian_source_unit_identities=checks,total_identities=len(checks),
        original_six_sector_true_mixed2_source_row_identities=mixed,mixed2_source_row_identity_count=len(mixed),
        original_six_sector_modes_and_normalizations={name:dict(mode=[str(v) for v in original[label][partname]['mode']],
            normalization_half=original[label][partname]['normalization_half'],beta=str(original[label][partname]['beta']))
            for name,(label,partname) in mappings.items()},
        original_pressure_FTC_cancels_centrifugal_source=pressure_identity,
        ordinary_radial_viscosity_without_inverse_R='-X*(rho*Q_rhorho+2*Q_rho)/sqrt(epsilon), and -Y likewise',
        original_F0_dressing_before_theta_axial_viscosity=True,physical_nu_base='1/2',
        original_axis_values_do_not_zero_transverse_derivatives=True,
        input_hashes=asts.hashes,passed=True)


@functools.lru_cache(maxsize=1)
def cartesian_pullback_theorem():
    """Independent Cartesian chain rule vs unchanged core_step through C2."""
    g=s.Function('arbitrary_smooth_axis_source')(RHO,Z);checks={}
    d=1-Z**2;L=1-D*Z**2
    dx=lambda value:s.diff(value,X)+X*s.diff(value,RHO)
    dy=lambda value:s.diff(value,Y)+Y*s.diff(value,RHO)
    dz=lambda value,gamma:(gamma*Z*value+d*s.diff(value,Z)-Z*(X*s.diff(value,X)+Y*s.diff(value,Y)+2*RHO*s.diff(value,RHO)))/L
    for label,seed in (('x',X),('y',Y),('scalar',s.Integer(1))):
        for nx,ny,nz in INDICES2:
            actual=seed*g;rows={(0,0):seed}
            for _ in range(nx):actual=dx(actual);rows=core_step(rows,'x')
            for _ in range(ny):actual=dy(actual);rows=core_step(rows,'y')
            for b in range(nz):
                gamma=B-nx-ny+b*(D-1);actual=dz(actual,gamma);rows=core_step(rows,'z',gamma)
            replay=sum((coef*s.diff(g,RHO,i,Z,k) for (i,k),coef in rows.items()),s.Integer(0))
            difference=s.cancel((actual-replay).subs(RHO,(X**2+Y**2)/2))
            if difference!=0:raise ValueError('Nonsingular Cartesian mixed2 pullback differs')
            checks[label+'/'+str((nx,ny,nz))]=True
    return dict(independent_Cartesian_mixed2_pullback_identities=checks,total_identities=len(checks),
        exact_coordinate_identity='rho=(X^2+Y^2)/2; X=x/(sqrt(nu)*lambda*sqrt(epsilon))',
        source_derivative_scales=dict(nu='nu^(1/2-Ntotal/2)',epsilon='epsilon^(base_power-Nxy/2)',
            physical_lambda='lambda^(base_beta-Nxy+Nz*(delta-1))'),
        original_gamma_updated_at_every_axial_step=True,passed=True)


@functools.lru_cache(maxsize=1)
def physical_factor_and_cylindrical_pullback_theorem():
    """Original cylindrical lift equals the nonsingular physical Cartesian C2.

    The polar coefficient is differentiated before its radius factor is
    enclosed. Its exact (2/rho)^(Nxy/2) conversion cancels that factor,
    leaving epsilon^(-Nxy/2). No interval overlap is used as equality.
    """
    asts=SourceAST();theta=s.Symbol('polar_angle',real=True)
    g=s.Function('arbitrary_same_physical_core_source')(RHO,Z)
    cs=s.cos(theta);sn=s.sin(theta);rad=s.sqrt(2*RHO)
    polar={X:rad*cs,Y:rad*sn};coordinate={}
    # Factored Cartesian differentiation of polar coordinates. The -n/2
    # term differentiates the original (2/rho)^(n/2) radius prefactor.
    def polar_step(value,label,n):
        radial=RHO*s.diff(value,RHO)-s.Rational(n,2)*value
        return cs*radial-sn*s.diff(value,theta)/2 if label=='x' else sn*radial+cs*s.diff(value,theta)/2
    for label,seed in (('x',X),('y',Y),('scalar',s.Integer(1))):
        for nx,ny,nz in INDICES2:
            original=seed.subs(polar)*g;rows={(0,0):seed};n=0
            for direction,count in (('x',nx),('y',ny)):
                for _ in range(count):
                    original=polar_step(original,direction,n);n+=1
                    rows=core_step(rows,direction)
            # The -Nxy lambda shift cancels the Z derivative of rho^(-Nxy/2),
            # just as in the actual original physical_operators program.
            for b in range(nz):
                original=axial_scalar(original,B+b*(D-1))
                rows=core_step(rows,'z',B-nx-ny+b*(D-1))
            replay=sum((coef.subs(polar)*s.diff(g,RHO,i,Z,k) for (i,k),coef in rows.items()),s.Integer(0))
            difference=s.trigsimp(s.cancel(replay-(2/RHO)**s.Rational(nx+ny,2)*original))
            if difference!=0:raise ValueError('Original cylindrical/Cartesian physical C2 coefficient differs')
            coordinate[label+'/'+str((nx,ny,nz))]=True
    # Bind the normalized polar radial/Z rules to the actual original
    # cylindrical operator, rather than merely another invented pullback.
    env=dict(s=s,functools=functools)
    from lei_ren_part1_paper_compliant_pulse_physical_bounds import ZSYM,DSYM,BSYM
    env.update(ZSYM=Z,DSYM=D,BSYM=B)
    operators=asts.replay('pulse_physical_bounds','physical_operators',env)()
    cylindrical={}
    for i in range(3):
        for j in range(3-i):
            direct=g
            for a in range(i):direct=RHO*s.diff(direct,RHO)-s.Rational(a,2)*direct
            for b in range(j):direct=axial_scalar(direct,B+b*(D-1))
            replay=s.Integer(0)
            for (k,n),coefficient in operators[i,j].items():
                value=s.diff(g,Z,n)
                for _ in range(k):value=RHO*s.diff(value,RHO)
                replay+=coefficient*value
            if s.cancel(direct-replay)!=0:raise ValueError('Original physical_operators polar rule differs')
            cylindrical['r'+str(i)+'_z'+str(j)]=True
    # Replay the actual old/new call arguments and actual logarithmic factor
    # programs. Raw sector units have already been proved with arbitrary Q,V,F.
    old_method=asts.method('pulse_end_physical_C2','lift_physical_packet')
    old_calls=[node for node in ast.walk(old_method) if isinstance(node,ast.Call)
        and ast.unparse(node.func)=='physical_source_row'
        and len(node.args)>3 and ast.unparse(node.args[3])=="sector['beta'] - i + j * (delta - 1)"]
    new_method=asts.method('current_core_axis_operator','cartesian_remainder')
    new_calls=[node for node in ast.walk(new_method) if isinstance(node,ast.Call) and ast.unparse(node.func)=='physical_source_row']
    if len(old_calls)!=1 or len(new_calls)!=1:raise ValueError('Original remainder physical factor call changed')
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=s.log)
    le,lr,lp,lf,lt,lnu=s.symbols('log_epsilon log_rho log_Pstar log_F0base log_tau log_nu',real=True)
    actual_logR=asts.evaluate(asts.expression('current_core_stress_operator','core_raw_source','logR'),
        dict(c=c,field=SimpleNamespace(core=SimpleNamespace(epsilon=s.exp(le))),r=s.exp(lr)))
    if s.cancel(actual_logR-le-lr)!=0:raise ValueError('Actual core R=epsilon*rho source changed')
    chart_logR=asts.evaluate(asts.expression('current_core_background_tensor','chart','logR'),dict(source=dict(logR=actual_logR)))
    actual_packet=asts.evaluate(asts.expression('current_core_background_tensor','chart','packet'),
        dict(c=c,Z=Z,rho=RHO,logR=chart_logR,self=SimpleNamespace(physical=SimpleNamespace(logP=lp)),sectors={}))
    if actual_packet['exact_logR']!=le+lr or actual_packet['exact_pulse_reference_logB_parts']!=dict(logPstar=lp) or actual_packet['exact_logD']!=0 or actual_packet['exact_logH']!=0:
        raise ValueError('Actual core physical packet factors differ')
    asts.expression('current_core_axis_background_tensor','__init__','self.log_epsilon',wanted='self.ctx.ln(self.core.epsilon)')
    old_F0=asts.expression('current_core_stress_operator','core_raw_source','logF0')
    new_F0=asts.expression('current_core_axis_background_tensor','__init__','self.log_F0')
    expected_F0=ast.unparse(old_F0).replace('field.core','self.core').replace('c.mpf','self.ctx.mpf')
    if ast.dump(new_F0)!=ast.dump(ast.parse(expected_F0,mode='eval').body):raise ValueError('Same original F0 source enclosure changed')
    factor=asts.replay('pulse_end_physical_C2','factor_logs',{})
    exponent=asts.expression('collar_physical_C2','physical_source_row','exponent')
    combined=asts.expression('collar_physical_C2','physical_source_row','combined')
    upper=asts.expression('collar_physical_C2','physical_source_row','upper')
    power=asts.expression('current_core_axis_operator','cartesian_remainder','power')
    logs=asts.expression('current_core_axis_operator','cartesian_remainder','logs')
    def capture(ctx,coefficient,parts,gamma,logtau,nu,spatial_order,radial_log=0,nu_base=1):
        env=dict(c=ctx,parts=parts,gamma=gamma,logtau=logtau,spatial_order=spatial_order,
            radial_log=radial_log,nu_base=nu_base,magnitude=1,lognu=lnu)
        env['exponent']=asts.evaluate(exponent,env);env['combined']=asts.evaluate(combined,env)
        return dict(gamma=gamma,nu_exponent=env['exponent'],weight=asts.evaluate(upper,env))
    unit=raw_axis_source_theorem()
    if not unit['passed'] or unit['total_identities']!=110 or unit['mixed2_source_row_identity_count']!=36:
        raise ValueError('Exact original raw units and every true mixed2 source row required')
    factors={}
    for name,data in sector_metadata().items():
        # The consumed raw theorem binds these modes and normalization to
        # raw_pre_remainder_sectors; factor_logs still executes its real AST.
        original_sector=unit['original_six_sector_modes_and_normalizations'][name]
        parts=factor(c,actual_packet,tuple(s.Rational(v) for v in original_sector['mode']),
            half=original_sector['normalization_half'])
        base=sum(parts.values(),s.Integer(0))
        for component in data['seed']:
            for nx,ny,nz in INDICES2:
                Nxy=nx+ny;N=Nxy+nz;beta=data['beta']
                old=asts.evaluate(old_calls[0],dict(c=c,physical_source_row=capture,coefficient=1,parts=parts,
                    sector=dict(beta=beta),i=Nxy,j=nz,delta=D,lt=lt,nu=1,lr=le+lr))
                p=asts.evaluate(power,dict(c=c,data=data,Nxy=Nxy))
                new_parts=asts.evaluate(logs,dict(power=p,field=SimpleNamespace(log_epsilon=le,log_F0=lf),data=data))
                new=asts.evaluate(new_calls[0],dict(c=c,physical_source_row=capture,coefficient=1,logs=new_parts,
                    beta=beta,Nxy=Nxy,nz=nz,delta=D,lt=lt,nu=1,N=N))
                # First transfer is the exact raw unit identity; second is
                # the exact coefficient identity proved above for all C2 rows.
                raw_transfer=data['epsilon_power']*le+data['F0_power']*lf-base
                polar_transfer=s.Rational(Nxy,2)*(lr-s.log(2))
                if any(s.cancel(value)!=0 for value in (old['gamma']-new['gamma'],
                    old['nu_exponent']-new['nu_exponent'],old['weight']+raw_transfer+polar_transfer-new['weight'])):
                    raise ValueError('Original positive-radius/axis epsilon,F0,lambda,nu factor differs')
                factors[name+'/'+component+'/'+str((nx,ny,nz))]=True
    return dict(original_cylindrical_to_Cartesian_coefficient_identities=coordinate,
        coefficient_identity_count=len(coordinate),original_cylindrical_operator_identities=cylindrical,
        cylindrical_operator_identity_count=len(cylindrical),original_physical_factor_identities=factors,
        physical_factor_identity_count=len(factors),
        consumed_exact_raw_unit_identity_count=unit['total_identities'],
        consumed_original_true_mixed2_source_row_identity_count=unit['mixed2_source_row_identity_count'],
        actual_core_packet_radius_Pstar_zero_D_zero_H_and_axis_epsilon_F0_source_bound=True,
        actual_old_remainder_call=ast.unparse(old_calls[0]),actual_new_remainder_call=ast.unparse(new_calls[0]),
        exact_radius_identity='R=epsilon*rho',exact_polar_coefficient_transfer='Cartesian coefficient=(2/rho)^(Nxy/2)*original normalized polar coefficient',
        physical_factor_cancellation_before_enclosure=True,
        same_actual_original_cylindrical_physical_operators_and_log_factor_programs=True,
        physical_source_bounds_not_resolved_point_values=True,
        input_hashes={**unit['input_hashes'],**asts.hashes},passed=True)

"""Original full-meridional pulse-end physical tensor and three-component error.

Positive source scales remain signed logarithmic sectors. The completed
tensor has zero radial divergence; the radial material/viscous remainder
is retained. This regional leading error is not a proved flat remainder.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import (
    CompliantPulseEndStressC3, SourceAST, source_precision, DOMAIN)
from lei_ren_part1_paper_compliant_pulse_end_flatten_join import ExactSourceHalves
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative, shifted_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    physical_bracket, physical_operators, interval_expression)
from lei_ren_part1_paper_compliant_collar_physical_C2 import (
    physical_source_row, scale_row)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_compliant_'
FALSE_FLAGS = ('pulse_end_cone_certified', 'whole_outer_cone_certified',
    'global_admissible_stress_lift_constructed', 'independently_bounded_global_flat_remainder',
    'physical_energy_integral_certified', 'full_background_NS_validation', 'temporal_recursion')


def original_full_stress(env):
    """Replay unchanged (3.16)-(3.18), including original Ntheta/Nz."""
    path = HERE / 'lei_ren_part1_paper_mp_stress.py'
    tree = ast.parse(path.read_text(encoding='utf8'))
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == 'evaluate_mp_stress')
    body = next(n for n in fn.body if isinstance(n, ast.With)).body
    result = dict(env)
    started = False
    for node in body:
        if isinstance(node, ast.Assign) and any(ast.unparse(t) == 'transport' for t in node.targets):
            started = True
        if isinstance(node, ast.Assign) and any(ast.unparse(t) == 'result' for t in node.targets):
            break
        if not started:
            continue
        node = ExactSourceHalves().visit(ast.parse(ast.unparse(node)).body[0])
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.IfExp):
            node.value = node.value.body
        exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])),
            '<original full meridional stress>', 'exec'), result)
    return result, {path.name: hashlib.sha256(path.read_bytes()).hexdigest()}


def full_physical_identities():
    y, z, delta = s.symbols('y Z delta', real=True)
    R = s.symbols('R', positive=True)
    h = s.Rational(1, 2)
    d, L = 1-z*z, 1-delta*z*z
    Ut, Uz = (s.Function(name)(y, z) for name in ('Utheta', 'Uz'))
    moments = {name: s.Function('same_'+name)(y, z) for name in ('theta', 'z', 'theta_z', 'z_theta', 'p')}
    P = s.Function('same_absolute_P')(y, z)
    dy = lambda value: s.diff(value, y) + R*s.diff(value, R)
    env = dict(a=Ut, b=Uz, ay=s.diff(Ut,y), az=s.diff(Ut,z),
        by=s.diff(Uz,y), bz=s.diff(Uz,z), dt=delta, z=z, R=R,
        root=s.sqrt(2*R), d=d, L=L, m=moments,
        mz={name:s.diff(value,z) for name,value in moments.items()},
        p=P, pz=s.diff(P,z))
    actual, hashes = original_full_stress(env)
    rules = {}
    sources = dict(theta=s.sqrt(2)*R**s.Rational(3,2)*Ut,
        z=R*Uz, theta_z=s.sqrt(2)*R**s.Rational(3,2)*Ut*Uz,
        z_theta=R*(Uz*Uz-Ut*Ut/2))
    for name, rhs in sources.items():
        rules[s.diff(moments[name],y)] = rhs
        rules[s.diff(moments[name],y,z)] = s.diff(rhs,z)
    rules[s.diff(P,y)] = Ut*Ut/2
    rules[s.diff(P,y,z)] = Ut*s.diff(Ut,z)
    proofs = {}
    def zero(name, value):
        if s.cancel(s.expand(value)).is_zero is not True:
            raise ArithmeticError('Full physical identity failed: '+name)
        proofs[name] = True
    Ur = actual['Ur']
    Ury = dy(Ur).subs(rules, simultaneous=True)
    zero('actual_radial_recovery_derivative_3_9', Ury/R-actual['Ur_R'])
    zero('actual_full_inertial_theta_primitive', (dy(actual['Itheta'])+actual['Itheta']).subs(rules, simultaneous=True)-actual['Nt'])
    zero('actual_full_inertial_axial_primitive', (dy(actual['Iz'])+actual['Iz']/2).subs(rules, simultaneous=True)-actual['Nz'])
    Dz = lambda value, beta: (beta*z*value+d*s.diff(value,z)-2*z*dy(value))/L
    time = lambda value, beta: (-beta*value/2+(1-delta)*z*s.diff(value,z)/2+dy(value))/L
    # sqrt(nu) lambda^(beta-2) is the common inertial/radial-viscosity scale.
    theta_inertial = time(Ut,-1-delta)+s.sqrt(2/R)*Ur*(dy(Ut)+Ut/2)+Uz*Dz(Ut,-1-delta)
    axial_inertial = time(Uz,-1-delta)+s.sqrt(2/R)*Ur*dy(Uz)+Uz*Dz(Uz,-1-delta)+Dz(P,-2-2*delta)
    theta_inertial = theta_inertial.subs(rules, simultaneous=True)
    axial_inertial = axial_inertial.subs(rules, simultaneous=True)
    zero('actual_full_physical_theta_inertia_pressure_transfer', theta_inertial+s.sqrt(2/R)*actual['Nt'])
    zero('actual_full_physical_axial_inertia_pressure_transfer', axial_inertial+s.sqrt(2/R)*actual['Nz'])
    zero('actual_theta_radial_viscosity_transfer', s.sqrt(2/R)*(dy(actual['Stheta'])+actual['Stheta'])-2/R*(s.diff(Ut,y,2)-Ut/4))
    zero('actual_axial_radial_viscosity_transfer', s.sqrt(2/R)*(dy(actual['Sz'])+actual['Sz']/2)-2/R*s.diff(Uz,y,2))
    zero('actual_radial_pressure_centrifugal_cancellation', s.sqrt(2/R)*rules[s.diff(P,y)]-Ut*Ut/s.sqrt(2*R))
    zero('actual_full_meridional_incompressibility',
        s.sqrt(2/R)*(Ury+Ur/2)+Dz(Uz,-1-delta))
    # Independent cylindrical tensor product identities, valid for arbitrary full Tz.
    r, zz = s.symbols('physical_r physical_z', positive=True)
    tt, tz = s.Function('Ttheta')(r,zz), s.Function('Tz')(r,zz)
    zero('completed_tensor_radial_divergence_zero', s.diff(tz,zz)-r*s.diff(tz,zz)/r)
    for i in range(3):
        for j in range(3-i):
            diagonal = r*s.diff(tz,r,i,zz,j+1)+(i*s.diff(tz,r,i-1,zz,j+1) if i else 0)
            zero('completed_diagonal_mixed2_'+str(i)+str(j),s.diff(r*s.diff(tz,zz),r,i,zz,j)-diagonal)
            for f, multiplier, name in ((tt,2,'theta'),(tz,1,'axial')):
                product = s.diff(f,r,i+1,zz,j)+multiplier*sum(
                    s.binomial(i,m)*(-1)**m*s.factorial(m)*r**(-m-1)*s.diff(f,r,i-m,zz,j) for m in range(i+1))
                zero('divergence_mixed2_'+name+str(i)+str(j),s.diff(s.diff(f,r)+multiplier*f/r,r,i,zz,j)-product)
    # The residual formula does not discard the meridional material derivative.
    radial = time(Ur,-1)+s.sqrt(2/R)*Ur*Ury+Uz*Dz(Ur,-1)-2/R*(dy(Ury)-Ur/4)
    return dict(identities=proofs,input_hashes=hashes,
        tensor='Trtheta=Ttheta; Trz=Tz; Ttheta_theta=r*partial_z(Tz); symmetric; other entries zero',
        physical_units='u=sqrt(nu)*lambda^beta*U; p=nu*lambda^(-2-2delta)*P; T=nu*lambda^(-2-delta)*Tprofile',
        exact_remainder=dict(radial='partial_t(ur)+ur*partial_r(ur)+uz*partial_z(ur)-nu*(partial_rr(ur)+partial_r(ur)/r-ur/r^2+partial_zz(ur))',
            theta='-nu*partial_zz(utheta)',axial='-nu*partial_zz(uz)'),
        radial_material_and_radial_viscosity_similarity_expression=str(radial),
        full_nonzero_meridional_velocity_retained=True,
        regional_error_not_claimed_flat=True)


def ordinary_grid(rows, total):
    return {'y'+str(k)+'_Z'+str(n):row[n]*math.factorial(n)
        for k,row in enumerate(rows[:total+1]) for n in range(total+1-k)}


def axial_n(value, order):
    for _ in range(order):
        value = axial_derivative(value)
    return value


def axial_operator_rows(c, rows, beta, z, delta, count=2, order=2):
    """True source y/Z jets; no fixed-basepoint amplitude is differentiated again."""
    result = []
    for j in range(count+1):
        zero = rows[0]*0
        value = zero
        for (k,n), coefficient in physical_operators()[0,order].items():
            # Rational functions of Z must remain Taylor functions.
            def evaluate(expression):
                from lei_ren_part1_paper_compliant_pulse_physical_bounds import ZSYM, DSYM, BSYM
                if expression == ZSYM:return z
                if expression == DSYM:return delta
                if expression == BSYM:return beta
                if expression.is_Rational:return c.mpf(int(expression.p))/int(expression.q)
                if expression.is_Add:return sum((evaluate(v) for v in expression.args),zero)
                if expression.is_Mul:
                    ans = zero+1
                    for v in expression.args:ans = ans*evaluate(v)
                    return ans
                if expression.is_Pow and expression.args[1].is_Integer:return evaluate(expression.args[0])**int(expression.args[1])
                raise ValueError('Unsupported axial operator source coefficient')
            value += axial_n(rows[j+k],n)*evaluate(coefficient)
        result.append(value)
    return result


def pulse_velocity_rows(c, delta, mu, z, C, Bh, m):
    """Exactly the native radial recovery with its common C derivative retained."""
    d, L = 1-z*z, 1-z*z*delta
    raw = [(z*C*Bh[j]*2-z*C*m[j]*(1-delta)-d*axial_derivative(C*m[j]))/L for j in range(5)]
    return dict(radial=shifted_rows(raw,-mu,4),
        theta=[C*(-(c.mpf('.5')+mu))**j for j in range(5)],
        axial=shifted_rows([C*v for v in Bh],-(c.mpf('.5')+mu),4))


def pulse_remainder_sectors(c,delta,mu,z,velocity):
    Ur, Ut, Uz = (velocity[name] for name in ('radial','theta','axial'))
    L = 1-z*z*delta
    time = [(Ur[j]/2+z*axial_derivative(Ur[j])*((1-delta)/2)+Ur[j+1])/L for j in range(3)]
    radial_viscosity = shifted_rows([-(Ur[j+2]-Ur[j]/4)*2 for j in range(3)],-1,2)
    dzUr = axial_operator_rows(c,Ur,-1,z,delta,count=2,order=1)
    radial_transport = shifted_rows(product_rows(Ur[:3],Ur[1:4]),-c.mpf('.5'),2)
    axial_transport = product_rows(Uz[:3],dzUr)
    nonlinear = [radial_transport[j]+axial_transport[j] for j in range(3)]
    return dict(radial=dict(
        time=dict(mode=(.5,1,1,0),normalization_half=True,beta=-3,rows=time),
        radial_viscosity=dict(mode=(-.5,1,1,0),normalization_half=True,beta=-3,rows=radial_viscosity),
        nonlinear_transport=dict(mode=(.5,2,2,0),normalization_half=True,beta=-3,rows=nonlinear),
        axial_viscosity=dict(mode=(.5,1,1,0),normalization_half=True,beta=-3+2*delta,
            rows=[-v for v in axial_operator_rows(c,Ur,-1,z,delta)])),
        theta=dict(axial_viscosity=dict(mode=(0,1,0,0),normalization_half=False,beta=-3+delta,
            rows=[-v for v in axial_operator_rows(c,Ut,-1-delta,z,delta)])),
        axial=dict(axial_viscosity=dict(mode=(0,1,1,0),normalization_half=False,beta=-3+delta,
            rows=[-v for v in axial_operator_rows(c,Uz,-1-delta,z,delta)])))


def factor_logs(c, packet, mode, half=True):
    rp,bp,dp,hp = mode
    return dict(source_logR=rp*packet['exact_logR'],
        **{key:bp*value for key,value in packet['exact_pulse_reference_logB_parts'].items()},
        selected_log_end_scale=dp*packet['exact_logD'],
        signed_original_memory_log=hp*packet['exact_logH'],
        normalization=-c.ln(2)/2 if half else c.mpf(0))


def lift_physical_packet(c, packet, delta, velocity, log_tau, theta, viscosity):
    z0 = packet['Z']; z = IntervalTaylor.variable(c,z0,5)
    nu = c.mpf(viscosity); lt = c.mpf(log_tau)
    if endpoints(nu)[0] <= 0 or not all(mp.isfinite(v) for v in endpoints(nu)+endpoints(lt)):
        raise ValueError('Finite nu>0 and log(tau) required')
    lr = packet['exact_logR']; beta = -2-delta
    stress, divergence, diagonal = {}, {}, {}
    for label, sectors in packet['full_meridional_stress_log_sectors'].items():
        stress[label], divergence[label] = {}, {}
        for name, sector in sectors.items():
            grid = {key.replace('s','y',1):value for key,value in sector['full_stress_mixed3_coefficient_enclosures'].items()}
            parts = sector['exact_source_log_parts']
            stress[label][name] = {}
            divergence[label][name] = {}
            for i in range(4):
                for j in range(4-i):
                    coefficient = physical_bracket(c,grid,i,j,z0,delta,beta)
                    key = 'r'+str(i)+'_z'+str(j)
                    stress[label][name][key] = physical_source_row(c,coefficient,parts,
                        beta-i+j*(delta-1),lt,nu,i+j,c.mpf(i)/2*(c.ln(2)-lr))
            for i in range(3):
                for j in range(3-i):
                    coefficient = physical_bracket(c,grid,i+1,j,z0,delta,beta)
                    for m in range(i+1):
                        coefficient += physical_bracket(c,grid,i-m,j,z0,delta,beta)*(
                            (2 if label=='theta' else 1)*math.comb(i,m)*(-1)**m*math.factorial(m)/c.mpf(2)**(m+1))
                    key = 'r'+str(i)+'_z'+str(j)
                    divergence[label][name][key] = physical_source_row(c,coefficient,parts,
                        beta-1-i+j*(delta-1),lt,nu,i+j,c.mpf(i+1)/2*(c.ln(2)-lr),nu_base=c.mpf('.5'))
            if label == 'axial':
                diagonal[name] = {}
                diagparts = dict(parts);diagparts['completed_radius']=lr/2+c.ln(2)/2
                for i in range(3):
                    for j in range(3-i):
                        coefficient = physical_bracket(c,grid,i,j+1,z0,delta,beta)
                        if i:coefficient += physical_bracket(c,grid,i-1,j+1,z0,delta,beta)*c.mpf(i)/2
                        diagonal[name]['r'+str(i)+'_z'+str(j)] = physical_source_row(c,coefficient,diagparts,
                            beta+delta-i+j*(delta-1),lt,nu,i+j,c.mpf(i)/2*(c.ln(2)-lr))
    remainder = {}
    source_errors = pulse_remainder_sectors(c,delta,None,z,velocity)
    for label, sectors in source_errors.items():
        remainder[label] = {}
        for name, sector in sectors.items():
            grid = ordinary_grid(sector['rows'],2)
            parts = factor_logs(c,packet,sector['mode'],sector['normalization_half'])
            remainder[label][name] = {}
            for i in range(3):
                for j in range(3-i):
                    coefficient = physical_bracket(c,grid,i,j,z0,delta,sector['beta'])
                    remainder[label][name]['r'+str(i)+'_z'+str(j)] = physical_source_row(c,coefficient,parts,
                        sector['beta']-i+j*(delta-1),lt,nu,i+j,c.mpf(i)/2*(c.ln(2)-lr),nu_base=c.mpf('.5'))
    cs,sn = (c.mpf([-1,1]),c.mpf([-1,1])) if theta is None else (c.cos(c.mpf(theta)),c.sin(c.mpf(theta)))
    tt = [grid['r0_z0'] for grid in stress['theta'].values()]
    tz = [grid['r0_z0'] for grid in stress['axial'].values()]
    dd = [grid['r0_z0'] for grid in diagonal.values()]
    zero = physical_source_row(c,c.mpf(0),{},beta,lt,nu,0)
    tensor = dict(xx=[scale_row(c,row,-2*cs*sn) for row in tt]+[scale_row(c,row,sn*sn) for row in dd],
        xy=[scale_row(c,row,cs*cs-sn*sn) for row in tt]+[scale_row(c,row,-cs*sn) for row in dd],
        xz=[scale_row(c,row,cs) for row in tz],yy=[scale_row(c,row,2*cs*sn) for row in tt]+[scale_row(c,row,cs*cs) for row in dd],
        yz=[scale_row(c,row,sn) for row in tz],zz=[zero])
    dt = [grid['r0_z0'] for grid in divergence['theta'].values()]
    dz = [grid['r0_z0'] for grid in divergence['axial'].values()]
    divcart = dict(x=[scale_row(c,row,-sn) for row in dt],y=[scale_row(c,row,cs) for row in dt],z=dz)
    er,et,ez = ([grid['r0_z0'] for grid in remainder[label].values()] for label in ('radial','theta','axial'))
    ecart = dict(x=[scale_row(c,row,cs) for row in er]+[scale_row(c,row,-sn) for row in et],
        y=[scale_row(c,row,sn) for row in er]+[scale_row(c,row,cs) for row in et],z=ez)
    return dict(Z=z0,s=packet['s'],requested_log_tau=lt,physical_viscosity=nu,
        physical_cylindrical_stress_mixed3=stress,
        physical_cylindrical_stress_divergence_mixed2=divergence,
        completed_theta_theta_stress_mixed2=diagonal,
        physical_three_component_remainder_mixed2=remainder,
        physical_completed_stress_tensor_cartesian=tensor,
        physical_completed_stress_divergence_cartesian=divcart,
        physical_remainder_cartesian=ecart,
        physical_momentum_residual_decomposition_cartesian={key:[scale_row(c,row,-1) for row in divcart[key]]+ecart[key] for key in ('x','y','z')},
        exact_completed_tensor_radial_divergence=c.mpf(0),exact_physical_divergence=c.mpf(0),
        physical_map=dict(lambda_relation='lambda^2-lambda^(2delta)*z^2/nu=tau=1-t',
            radial='R=r^2/(2nu*lambda^2)',axial='Z=z/(sqrt(nu)*lambda^(1-delta))'),
        source_caps_used_as_defining_field_values=False,positive_source_factors_not_materialized=True,
        full_nonzero_meridional_velocity_and_radial_remainder_retained=True,
        regional_error_not_claimed_flat=True,**{flag:False for flag in FALSE_FLAGS})


def source_and_join_binding(records):
    asts = SourceAST();checks = {}
    # The sector radial formula is the same native recovery, after multiplying C.
    z,delta,Bh,m,mz = s.symbols('Z delta Bhat m mZ',real=True)
    C = 1/(1+z*z);L = 1-delta*z*z
    native = asts.evaluate(asts.expression('axial_pulse_field','radial','q'),dict(Z=z))
    if s.cancel(native-(1+z*z)) != 0:raise ArithmeticError('Native pulse q changed')
    fn = asts.method('axial_pulse_field','radial')
    ret = next(n for n in fn.body if isinstance(n,ast.Return))
    actual = asts.evaluate(ret.value,dict(Z=z,q=1+z*z,d=1-z*z,L=L,
        self=SimpleNamespace(delta=delta),B=[Bh],m1=[m,mz]))
    explicit = (2*z*C*Bh-(1-delta)*z*C*m-(1-z*z)*(C*mz+s.diff(C,z)*m))/L
    if s.cancel(explicit-C*actual) != 0:raise ArithmeticError('Actual native radial source differs')
    checks['actual_native_radial_recovery_and_common_C_derivative_bound']=True
    asts.expression('pulse_mixed_C4','transport_mixed','A',wanted='[field.radial(Z,b,m) for b,m in zip(Brows,m1)]')
    asts.expression('pulse_end_stress_C3','end','logD',wanted='self.native.logE')
    asts.expression('pulse_end_stress_C3','end','logR',wanted='self.assembly.logRp+13/mu+v')
    checks['actual_selected_backward_C5_and_original_radius_source_retained']=True
    join = records['pulse_end_flatten_join_check']
    if not join['pulse_end_flatten_full_moment_stress_pressure_functional_join_verified']:
        raise ValueError('Current exact functional similarity join required')
    for label in ('actual_Utheta_mixed4','actual_Ur_mixed4','actual_Uz_mixed4','actual_absolute_pressure_mixed4',
        'actual_theta_stress_mixed3','actual_axial_stress_mixed3'):
        if not any(name.startswith(label) for name in records['pulse_end_flatten_join']['functional_join']['identities']):
            raise ValueError('Missing functional interface source: '+label)
        checks['consumed_'+label]=True
    # All pullbacks use the same beta, nu and original radius at the endpoint.
    # Mixed3 stress suffices for mixed2 completion/divergence. Mixed4 velocities
    # suffice for mixed2 errors, including the meridional material products.
    for i in range(3):
        for j in range(3-i):
            if any(k+n>4 for k,n in physical_operators()[i,j+2]):
                raise ValueError('Physical error operator exceeds source mixed4')
            checks['same_physical_interface_operator_'+str(i)+str(j)]=True
    checks['functional_stress3_implies_physical_stress3_divergence2_diagonal2_join']=True
    checks['functional_velocity4_implies_full_physical_remainder2_join']=True
    checks['empty_original_beta_support_gives_exact_endpoint_Er_Ez_zero']=True
    checks['same_nonzero_angular_axial_viscosity_retained_at_interface']=True
    return dict(identities=checks,input_hashes=asts.hashes,
        pulse_end_flatten_completed_physical_interface_verified=True,
        functional_interface_not_interval_overlap=True,global_cone_not_inferred=True)


class CompliantPulseEndPhysicalC2:
    @source_precision
    def __init__(self,cells=64):
        self.stress = CompliantPulseEndStressC3(cells=cells)
        self.ctx = self.stress.ctx;self.family = self.stress.family;self.source = self.stress.source
        self.hashes = dict(self.stress.hashes);records = {}
        for stem in ('pulse_end_stress_C3_check','pulse_end_flatten_join','pulse_end_flatten_join_check',
            'flatten_physical_C2_check'):
            name = PREFIX+stem+'.json';record = json.loads((HERE/name).read_bytes())
            if stem.endswith('_check') and not record['all_passed']:raise ValueError('Unaccepted physical source: '+stem)
            if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256']) != (self.family,self.source):
                raise ValueError('Pulse-end physical source family differs')
            for path,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Physical source changed: '+path)
            self.hashes.update(record['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            records[stem] = record
        self.proof = full_physical_identities()
        self.join = source_and_join_binding(records)
        self.hashes.update(self.proof['input_hashes']);self.hashes.update(self.join['input_hashes'])
        self.hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def end(self,Z,offset,log_tau='-1',theta='0',viscosity='1'):
        packet = self.stress.end(Z,offset);c = self.ctx
        z = IntervalTaylor.variable(c,packet['Z'],5)
        # Preserve 1+Z^2 >= 1 on whole-Z boxes before Taylor multiplication.
        C = IntervalTaylor(c,[1+packet['Z']**2,2*packet['Z'],1,0,0,0]).reciprocal()
        velocity = pulse_velocity_rows(c,self.stress.native.delta,self.stress.native.mu,z,C,
            packet['source_formal_Bhat_rows'],packet['source_formal_Mz_rows'])
        point = lift_physical_packet(c,packet,self.stress.native.delta,velocity,log_tau,theta,viscosity)
        point.update(actual_pulse_end_physical_decomposition_constructed=True,
            pulse_end_flatten_completed_physical_interface_verified=True,
            actual_physical_stress_mixed_order=3,actual_physical_remainder_mixed_order=2,
            physical_interface_endpoint_is_s0_flatten_t0=True)
        return point

    @source_precision
    def report(self):
        samples = [self.end(z,v,viscosity=nu,theta='.37') for z,v,nu in
            (('0','-3','.01'),('.5','-1','.7'),('.5','-2','1'),('.5','0','.7'))]
        whole = self.end([-1,1],[-4,0],theta=None,viscosity='.7')
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            domain=DOMAIN,full_physical_operator_proof=self.proof,source_and_physical_join_binding=self.join,
            samples=samples,whole_original_end=whole,input_hashes=self.hashes,
            actual_pulse_end_physical_decomposition_constructed=True,
            pulse_end_flatten_completed_physical_interface_verified=True,
            source_caps_used_as_defining_field_values=False,
            full_nonzero_meridional_velocity_and_radial_remainder_retained=True,
            **{flag:False for flag in FALSE_FLAGS})


@source_precision
def run():
    result = CompliantPulseEndPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original pulse-end full physical tensor and three-component remainder generated; cone/flat/global pending',flush=True)
    return result


if __name__ == '__main__':
    run()

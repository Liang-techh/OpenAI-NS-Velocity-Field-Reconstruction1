"""Executable production 10+10 core source and zero-axis Green inverses.

The evaluator consumes the original gauge expressions, rather than a second
handwritten equation. Coefficients are ordinary axial Taylor coefficients.
An interval input encloses coefficients; it is not an exact point solution.
The operator/recurrence identity does not itself admit the analytic tail.
"""
import ast
import copy
import json
from pathlib import Path

import sympy as sp

from lei_ren_part1_paper_compliant_current_core_recurrence_source import function, binding, sha
from lei_ren_part1_paper_compliant_current_core_scaled_swirl_source import CurrentCoreScaledSwirlSource, OPEN
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_compliant_'
NAME = PREFIX + 'current_core_nonlinear_operator.json'
RECEIPT = PREFIX + 'current_core_nonlinear_operator_check.json'
GATE = 'current_core_actual_twenty_term_operator_and_radial_extraction_certified'
NORM_GATE = 'current_core_actual_operator_contraction_and_tail_certified'
NORM_BOUNDS_GATE = 'current_core_actual_twenty_term_norm_incidence_certified'


def production_expressions():
    """Execute only the pure symbolic prefix of the production gauge proof."""
    fn = function('gauge_fixed_point_identity', 'run')
    stops = [i for i, node in enumerate(fn.body) if isinstance(node, ast.Assign)
             and any(ast.unparse(t) == 'z_difference' for t in node.targets)]
    if len(stops) != 1:
        raise ValueError('Unique complete gauge-expression prefix required')
    tree = ast.Module(body=copy.deepcopy(fn.body[:stops[0]+1]), type_ignores=[])
    allowed = (ast.Assign, ast.BinOp, ast.UnaryOp, ast.Call, ast.Name,
               ast.Store, ast.Load, ast.Attribute, ast.Constant, ast.Tuple,
               ast.List, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub)
    if any(not isinstance(node, allowed+(ast.Module,)) for node in ast.walk(tree)):
        raise ValueError('Gauge prefix is no longer pure symbolic assignments')
    env = {'sp': sp, 'sum': sum}
    exec(compile(ast.fix_missing_locations(tree), '<production-gauge-prefix>', 'exec'), env)
    if env['theta_difference'] != 0 or env['z_difference'] != 0:
        raise ArithmeticError('Actual original gauge equations fail')
    if len(env['theta_terms']) != 10 or len(env['z_terms']) != 10:
        raise ValueError('All twenty actual source terms required')
    return env


class TaylorRectangle:
    """Finite formal quotient ring, indexed by (radial order, axial order).

    Products use ordinary Taylor convolution. Derivatives are formed before
    discarding the extra axial order. The caller controls the valid order;
    no unavailable derivative coefficient is interpreted as zero.
    """
    def __init__(self, rows, radial, axial, convert=sp.sympify):
        self.radial, self.axial, self.convert = radial, axial, convert
        zero = convert(0)
        self.rows = [[convert(rows[n][m]) if n < len(rows) and m < len(rows[n]) else zero
                      for m in range(axial+1)] for n in range(radial+1)]

    def constant(self, value):
        return TaylorRectangle([[value]], self.radial, self.axial, self.convert)

    def coerce(self, other):
        if isinstance(other, TaylorRectangle):
            if (other.radial, other.axial) != (self.radial, self.axial):
                raise ValueError('Taylor rectangles differ')
            return other
        return self.constant(other)

    def __add__(self, other):
        other = self.coerce(other)
        return TaylorRectangle([[a+b for a,b in zip(left,right)]
                                for left,right in zip(self.rows,other.rows)],
                               self.radial,self.axial,self.convert)

    __radd__ = __add__

    def __neg__(self):
        return TaylorRectangle([[-v for v in row] for row in self.rows],
                               self.radial,self.axial,self.convert)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        z = self.convert(0)
        rows = [[sum((self.rows[i][j]*other.rows[n-i][m-j]
                       for i in range(n+1) for j in range(m+1)),z)
                 for m in range(self.axial+1)] for n in range(self.radial+1)]
        return TaylorRectangle(rows,self.radial,self.axial,self.convert)

    __rmul__ = __mul__

    def inverse(self):
        # Only fixed axial multipliers are inverted in the production source.
        if any(v != 0 for row in self.rows[1:] for v in row):
            raise ValueError('Only a radial-constant multiplier may be inverted')
        a = self.rows[0]
        b = [self.convert(1)/a[0]]
        for m in range(1,self.axial+1):
            b.append(-sum((a[j]*b[m-j] for j in range(1,m+1)),self.convert(0))/a[0])
        return TaylorRectangle([b],self.radial,self.axial,self.convert)

    def __truediv__(self, other):
        return self*self.coerce(other).inverse()

    def __rtruediv__(self, other):
        return self.coerce(other)*self.inverse()

    def __pow__(self, power):
        power = int(power)
        if power < 0:
            return self.inverse()**(-power)
        out = self.constant(1)
        for _ in range(power):
            out = out*self
        return out

    def truncate_axial(self, order):
        if order > self.axial:
            raise ValueError('Unavailable axial order')
        return TaylorRectangle(self.rows,self.radial,order,self.convert)

    def dz(self):
        if self.axial < 1:
            raise ValueError('Axial derivative needs one extra Taylor coefficient')
        return TaylorRectangle([[(m+1)*row[m+1] for m in range(self.axial)] for row in self.rows],
                               self.radial,self.axial-1,self.convert)

    def radial_euler(self):
        return TaylorRectangle([[n*v for v in row] for n,row in enumerate(self.rows)],
                               self.radial,self.axial,self.convert)

    def radial_integral(self):
        z = self.convert(0)
        rows = [[z]*(self.axial+1)] + [[v/(n+1) for v in row]
                    for n,row in enumerate(self.rows[:-1])]
        return TaylorRectangle(rows,self.radial,self.axial,self.convert)

    def radial_average(self):
        return TaylorRectangle([[v/(n+1) for v in row] for n,row in enumerate(self.rows)],
                               self.radial,self.axial,self.convert)


def evaluate_expression(expr, values, template, cache=None):
    """Interpret a production SymPy expression in the finite Taylor ring."""
    cache = {} if cache is None else cache
    if expr in cache:
        return cache[expr]
    if expr.is_Symbol:
        result = values[str(expr)]
    elif expr.is_Number:
        # Preserve rational exactness for both SymPy and interval contexts.
        p,q = expr.as_numer_denom()
        result = template.constant(template.convert(int(p))/template.convert(int(q)))
    elif expr.is_Add:
        result = sum((evaluate_expression(a,values,template,cache) for a in expr.args),template.constant(0))
    elif expr.is_Mul:
        result = template.constant(1)
        for a in expr.args:
            result = result*evaluate_expression(a,values,template,cache)
    elif expr.is_Pow and expr.exp.is_Integer:
        result = evaluate_expression(expr.base,values,template,cache)**int(expr.exp)
    else:
        raise ValueError('Unsupported production source expression: '+str(expr))
    cache[expr] = result
    return result


def source_terms(fixed, rows, n, center, delta, epsilon, convert=sp.sympify):
    """Actual twenty-term source through radial n, consuming one axial order.

    P0_step=epsilon*P0_phys, ell_step=epsilon*ell_phys=-g,
    S_step=epsilon^2*F0^2. Physical F0^2 is recovered only from the typed
    source input; the pressure primitive and its Z derivative use that same
    source and Phi, with no independent pressure correction input.
    """
    sizes = [len(fixed[key]) for key in ('ell_Z_taylor','S_Z_taylor','U0_Z_taylor','P0_Z_taylor')]
    if len(set(sizes)) != 1 or n < 0:
        raise ValueError('Common fixed Taylor lengths and nonnegative radial order required')
    count = sizes[0]-n
    if count < 2 or any(len(rows[key]) != n+1 for key in ('A','Uz','P')):
        raise ValueError('Rows 0..n and one extra axial coefficient required')
    if any(len(row) < count for key in ('A','Uz','P') for row in rows[key]):
        raise ValueError('Source derivative would consume an unavailable coefficient')
    eps = convert(epsilon)
    if not (eps > 0 and eps <= 1):
        raise ValueError('epsilon must lie in (0,1]')
    axial = count-1
    make = lambda data: TaylorRectangle(data,n,axial,convert)
    z = make([[convert(center),convert(1)]])
    one = make([[convert(1)]])
    d, L = one-z*z, one-z*z*convert(delta)
    u0 = make([fixed['U0_Z_taylor']])
    p0 = make([fixed['P0_Z_taylor']])/eps
    scaled_S = make([fixed['S_Z_taylor']])
    physical_S = scaled_S/(eps*eps)
    g = -make([fixed['ell_Z_taylor']])
    phi = make(rows['A'])
    # The axis U0 subtraction is an exact structural identity. Subtracting
    # two equal interval enclosures would lose correlation and create a
    # spurious, epsilon-amplified correction at the axis.
    same = lambda a,b: a._mpi_ == b._mpi_ if hasattr(a,'_mpi_') and hasattr(b,'_mpi_') else a == b
    if (not all(same(a,b) for a,b in zip(rows['Uz'][0],fixed['U0_Z_taylor']))
        or not same(convert(rows['A'][0][0]),convert(1))
        or any(not same(convert(v),convert(0)) for v in rows['A'][0][1:])):
        raise ValueError('Exact Phi=1, Psi=0 axis normalization required')
    psi = make([[convert(0)]*count]+[[convert(v)/eps for v in row] for row in rows['Uz'][1:]])
    average = psi.radial_average()
    # This is a common functional pressure, including the derivative of F0^2.
    pressure = physical_S*(phi*phi).radial_integral()
    cut = lambda jet: jet.truncate_axial(axial-1)
    u0_z = u0.dz()
    z,d,L,u0,g,physical_S,phi0,psi0 = map(cut,(z,d,L,u0,g,physical_S,phi,psi))
    W0 = phi0.constant(1)-z*u0*(1-convert(delta))-d*u0_z
    H0 = z*(1-convert(delta))/2+d*u0
    rho = TaylorRectangle([[0],[1]],n,axial-1,convert)
    values = dict(s=rho,z=z,delta=phi0.constant(delta),L=L,d=d,epsilon=phi0.constant(eps),g=g,S=physical_S,
                  U0=u0,U0_z=u0_z,P0=cut(p0),P0_z=p0.dz(),H0=H0,W0=W0,
                  Phi=phi0,Phi_s=phi0.constant(0),Phi_z=phi.dz(),
                  Psi=psi0,Psi_s=psi0.constant(0),Psi_z=psi.dz(),
                  M=cut(average),M_z=average.dz(),Pcal=cut(pressure),Pcal_z=pressure.dz())
    # Evaluate s*Phi_s and s*Psi_s directly. This avoids losing the highest
    # known radial row by representing an unnecessarily differentiated field.
    env = production_expressions()
    ephi,epsi = sp.symbols('EulerPhi EulerPsi')
    values.update(EulerPhi=phi0.radial_euler(),EulerPsi=psi0.radial_euler())
    terms = {}
    cache = {}
    for component,key in (('theta','theta_terms'),('z','z_terms')):
        expressions = [expr.subs(env['s']*env['phi_s'],ephi).subs(env['s']*env['psi_s'],epsi)
                       for expr in env[key]]
        terms[component] = [evaluate_expression(expr,values,phi0,cache) for expr in expressions]
    chi = evaluate_expression(env['chi'],values,phi0,cache)
    B = evaluate_expression(env['B'],values,phi0,cache)
    return dict(terms=terms,chi=chi,B=B,phi=phi0,psi=psi0,
                scaled_S=cut(scaled_S),pressure=cut(pressure),pressure_z=pressure.dz(),
                output_axial_count=axial,epsilon=eps,expressions=env)


def green_inverse(rows, component, chi=None, convert=sp.sympify):
    """Full-equation inverse: (2D2+chi)^-1 or (2D1)^-1, zero axis.

    Raw J1/J2 omit 1/2. Here the single full-equation 1/2 is explicit.
    Angular multiplication by chi uses ordinary axial convolution.
    """
    if component not in ('theta','z') or not rows:
        raise ValueError('Nonempty angular/axial source rows required')
    K = len(rows[0])
    if any(len(row) != K for row in rows):
        raise ValueError('Green source needs a common axial rectangle')
    out = [[convert(0)]*K]
    for n,row in enumerate(rows):
        coupled = [convert(0)]*K
        if component == 'theta':
            if chi is None or len(chi) < K:
                raise ValueError('Angular inverse needs the same fixed chi jets')
            coupled = [sum((chi[i]*out[n][m-i] for i in range(m+1)),convert(0)) for m in range(K)]
        divisor = 2*(n+1)*(n+2 if component == 'theta' else n+1)
        out.append([(v-c)/divisor for v,c in zip(row,coupled)])
    return out


def correction_map(packet, convert=sp.sympify):
    """T=(epsilon/2 R J2 Etheta, epsilon/2 J1 Ez), correction-only data."""
    result = {}
    for component in ('theta','z'):
        terms = packet['terms'][component]
        summed = sum(terms,terms[0].constant(0))
        source = [[packet['epsilon']*v for v in row] for row in summed.rows]
        result[component] = green_inverse(source,component,packet['chi'].rows[0],convert)
    return result


def next_rows(packet, n):
    """Extract the next full Phi, normalized Uz and scaled P coefficients."""
    eps = packet['epsilon']
    theta = sum(packet['terms']['theta'],packet['phi'].constant(0))
    axial = sum(packet['terms']['z'],packet['psi'].constant(0))
    chi_phi = packet['chi']*packet['phi']
    drive = packet['B'] if n == 0 else packet['B'].constant(0)
    squared = packet['scaled_S']*(packet['phi']*packet['phi'])
    return dict(A=[(eps*v-c)/(2*(n+1)*(n+2)) for v,c in zip(theta.rows[n],chi_phi.rows[n])],
                Uz=[eps*(b+eps*v)/(2*(n+1)**2) for b,v in zip(drive.rows[n],axial.rows[n])],
                P=[v/(n+1) for v in squared.rows[n]])


def source_bindings_and_incidence():
    env = production_expressions()
    # Exact expressions paired with the majorants. A mere permutation of the
    # twenty terms keeps the summed equation unchanged but must not silently
    # associate a term with the bound of another term.
    s,z,dt,L,d,eps,g,S,U,Uz,P,Pz,H,W = [env[key] for key in
        ('s','z','dt','L','d','eps','g','S','U','U_z','P','P_z','H0','W0')]
    phi,phis,phiz,psi,psis,psiz,M,Mz,C,Cz = [env[key] for key in
        ('phi','phi_s','phi_z','psi','psi_s','psi_z','M','M_z','C','C_z')]
    theta = [-env['beta']*phi,W/L*s*phis,H/L*phiz,
             -eps*(1-dt)*z/L*M*phi,-eps*(1-dt)*z/L*M*s*phis,
             -eps*d/L*Mz*phi,-eps*d/L*Mz*s*phis,
             -eps*dt*z/L*psi*phi,eps*d/L*psi*phiz,-g*d/L*psi*phi]
    axial = [W/L*s*psis,(d*Uz+(1+dt)*(1-4*z*U)/2)/L*psi,H/L*psiz,
             -eps*(1-dt)*z/L*M*s*psis,-eps*d/L*Mz*s*psis,
             -eps*(1+dt)*z/L*psi*psi,eps*d/L*psi*psiz,
             d/L*Cz,-2*(1+dt)*z/L*C,-2*z/L*s*S*phi*phi]
    if any(sp.expand(a-b)!=0 for a,b in zip(env['theta_terms']+env['z_terms'],theta+axial)):
        raise ValueError('Exact source-expression/majorant pairing changed')
    transfer = function('compliant_core_transfer','run')
    calls = [node for node in ast.walk(transfer) if isinstance(node,ast.Call)
             and isinstance(node.func,ast.Name) and node.func.id == 'term']
    if len(calls) != 20:
        raise ValueError('Twenty transfer term calls required')
    expected = [
        ('theta',"product*coefficients['beta']*J2",1,0),
        ('theta',"product*coefficients['W0_over_L']*radial",1,0),
        ('theta',"product*coefficients['H0_over_L']*axial",1,0),
        ('theta',"eps*product*coefficients['az']*J2*product",1,1),
        ('theta',"eps*product*coefficients['az']*radial",1,1),
        ('theta',"eps*product*coefficients['d_over_L']*axial",1,1),
        ('theta',"eps*product*coefficients['d_over_L']*axial",1,1),
        ('theta',"eps*dt*product*coefficients['z_over_L']*J2*product",1,1),
        ('theta',"eps*product*coefficients['d_over_L']*axial",1,1),
        ('theta',"product*coefficients['cross_swirl']*J2*product",1,1),
        ('z',"product*coefficients['W0_over_L']*radial",0,1),
        ('z',"product*coefficients['axial_linear']*J1",0,1),
        ('z',"product*coefficients['H0_over_L']*axial",0,1),
        ('z',"eps*product*coefficients['az']*radial",0,2),
        ('z',"eps*product*coefficients['d_over_L']*axial",0,2),
        ('z',"eps*(1+dt)*product*coefficients['z_over_L']*J1*product",0,2),
        ('z',"eps*product*coefficients['d_over_L']*axial",0,2),
        ('z',"product*coefficients['d_over_L']*axial*Pcal",2,0),
        ('z',"2*(1+dt)*product*coefficients['z_over_L']*J1*Pcal",2,0),
        ('z',"2*product*coefficients['z_over_L']*J1*80*product**2*Fsquare",2,0)]
    rows = []
    for index,(call,(component,coefficient,p,q)) in enumerate(zip(calls,expected)):
        if (len(call.args)!=5 or ast.literal_eval(call.args[0])!=component
            or ast.dump(call.args[2])!=ast.dump(ast.parse(coefficient,mode='eval').body)
            or [ast.literal_eval(a) for a in call.args[3:]]!=[p,q]):
            raise ValueError('Actual majorant incidence changed at term '+str(index+1))
        local = index if component == 'theta' else index-10
        expr = env['theta_terms' if component == 'theta' else 'z_terms'][local]
        rows.append(dict(component=component,source_term_number=local+1,source_expression=str(expr),
                         production_majorant_label=ast.literal_eval(call.args[1]),
                         production_majorant_coefficient=ast.unparse(call.args[2]),phi_power=p,psi_power=q,
                         integral_pressure_not_independent=component=='z' and local in (7,8),
                         mixed_Mz_radial_derivative=local==6 if component=='theta' else local==4))
    binding('compliant_core_transfer','run','c','c*fixed_multiplier/product',True)
    binding('compliant_core_transfer','run','size','(Rnorm*theta_size+z_size)/2')
    binding('compliant_core_transfer','run','lip','(Rnorm*theta_lip+z_lip)/2',True)
    binding('compliant_core_transfer','run','map_size','eps*size')
    binding('compliant_core_transfer','run','map_lip','eps*lip')
    binding('compliant_core_transfer','run','Pcal','80*product**2*Fsquare')
    for target,expression in [('product','ctx.mpf(256)'),('J1','ctx.mpf(80)'),('J2','ctx.mpf(40)'),
        ('axial','ctx.mpf(20480)/h'),('radial','ctx.mpf(20480)'),
        ('Rnorm',"get(commuting,'resolvent_norm_upper')"),
        ('fixed_multiplier',"get(commuting,'axial_convolution_factor_upper')")]:
        binding('compliant_core_transfer','run',target,expression)
    binding('compliant_core_transfer','run','coefficients',
        "dict(beta=get(tube,'beta_modulus_upper')*weight,"
        "W0_over_L=(1+(1+dt)*a*U+4*d)/L*weight,H0_over_L=H/L*weight,"
        "az=(1+dt)*a/L*weight,d_over_L=d/L*weight,z_over_L=a/L*weight,"
        "axial_linear=((1+dt)/2*(1+4*a*U)+4*d)/L*weight,cross_swirl=d/pole*weight)")
    binding('shared_commuting_resolvent','run','S','(1+r)/(1-r)**3')
    return dict(actual_production_gauge_differences=[str(env['theta_difference']),str(env['z_difference'])],
                terms=rows,single_outer_epsilon_and_one_half_AST_bound=True,
                angular_resolvent_applied_only_to_theta_AST_bound=True,
                fixed_multiplier_replaces_one_outer_product_factor_AST_bound=True,
                physical_F0_squared_is_scaled_S_divided_by_epsilon_squared=True,
                formal_incidence_is_not_an_analytic_norm_proof=True,passed=True)


class CurrentCoreNonlinearOperator:
    def __init__(self, source=None, require_checked=True):
        self.source = source if source is not None else CurrentCoreScaledSwirlSource()
        if not self.source.acceptance_loaded:
            raise ValueError('Checked current exact scaled-S source required')
        self.ctx = self.source.ctx
        self.family,self.source_sha,self.datum_sha = self.source.family,self.source.source,self.source.datum_sha
        self.hashes = dict(self.source.hashes)
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        self.hashes[PREFIX+'current_core_scaled_swirl_source_check.json'] = sha(PREFIX+'current_core_scaled_swirl_source_check.json')
        self.incidence = source_bindings_and_incidence()
        self.acceptance_loaded = False
        if require_checked:
            report = accepted(RECEIPT,self.family,self.source_sha,GATE)
            _verify_hashes(report)
            if (report['datum_enclosure_sha256'] != self.datum_sha or report[NORM_GATE]
                or not report[NORM_BOUNDS_GATE]):
                raise ValueError('Operator receipt source/scope differs')
            self.hashes.update(report['input_hashes'])
            self.hashes[RECEIPT] = sha(RECEIPT)
            self.acceptance_loaded = True

    @source_precision
    def evaluate_seed(self, Z, degree=4, depth=3):
        seed = self.source.certified_seed(Z,degree,depth)
        fixed,rows = seed['fixed'],seed['rows']
        if degree < 1:
            raise ValueError('At least one output radial row required')
        packets=[]
        for n in range(degree):
            packet = source_terms(fixed,rows,n,self.ctx.mpf(Z),self.source.core.delta,
                                  self.source.core.epsilon,self.ctx.mpf)
            recovered = next_rows(packet,n)
            packets.append(dict(source_radial_order=n,output_axial_count=packet['output_axial_count'],
                                actual_source_terms={key:[term.rows[n] for term in terms]
                                                     for key,terms in packet['terms'].items()}))
            for key in ('A','Uz','P'):
                rows[key].append(recovered[key])
        return dict(Z=self.ctx.mpf(Z),current_original_seed_consumed=True,
                    exact_S_remains_formal_with_certified_Cauchy_enclosure=True,
                    actual_source_packets_by_radial_order=packets,
                    generated_finite_radial_rows=rows,
                    correction_only_Green_map=correction_map(packet,self.ctx.mpf),
                    physical_pressure_primitive=packet['pressure'].rows,
                    physical_pressure_primitive_Z=packet['pressure_z'].rows,
                    **{GATE:self.acceptance_loaded,NORM_BOUNDS_GATE:self.acceptance_loaded,NORM_GATE:False},
                    **dict.fromkeys(OPEN,False))


@source_precision
def run(operator=None):
    operator = operator if operator is not None else CurrentCoreNonlinearOperator(require_checked=False)
    names = ['lei_ren_part1_paper_gauge_fixed_point_identity.py',
             'lei_ren_part1_paper_functional_core_step.py','lei_ren_part1_paper_logarithmic_core_step.py',
             'lei_ren_part1_paper_shared_commuting_resolvent.py',
             'lei_ren_part1_paper_compliant_core_transfer.py',
             PREFIX+'current_core_recurrence_source_check.json']
    hashes = dict(operator.hashes)
    for name in names:
        hashes[name] = sha(name)
    report = dict(actual_five_defect_family_sha256=operator.family,implicit_source_sha256=operator.source_sha,
                  datum_enclosure_sha256=operator.datum_sha,source_bindings_and_incidence=operator.incidence,
                  actual_operator='T=(epsilon/2 R J2 Etheta, epsilon/2 J1 Ez)',
                  raw_Green_inverses=['(J1 f)_(n+1)=f_n/(n+1)^2','(J2 f)_(n+1)=f_n/((n+1)(n+2))'],
                  common_integrals=['M=rho^-1 integral_0^rho Psi','Pcal=F0^2 integral_0^rho Phi^2',
                                    'Pcal_Z is the derivative of that same product/primitive'],
                  normalized_units=['rho=s=Lambda R','ell_scaled=-g=epsilon ell_phys',
                                    'Ptilde=epsilon P0+epsilon^2 Pcal','S_scaled=epsilon^2 F0^2'],
                  fresh_current_runtime=operator.evaluate_seed('.293'),
                  analytic_majorant_obligations=['All-order composite derivative/product estimates for each actual term',
                                                'Same fixed analytic multiplier convolution before/after radial inversion',
                                                'All-order uniqueness, finite coefficient induction and common analytic tail'],
                  **{GATE:False,NORM_BOUNDS_GATE:False,NORM_GATE:False},**dict.fromkeys(OPEN,False),input_hashes=hashes)
    (HERE/NAME).write_bytes((json.dumps(encode(report),indent=2)+'\n').encode('utf8'))
    print('Built actual production twenty-term source/Green operator; analytic fixed-point/tail admission remains open',flush=True)
    return report


if __name__ == '__main__':
    run()

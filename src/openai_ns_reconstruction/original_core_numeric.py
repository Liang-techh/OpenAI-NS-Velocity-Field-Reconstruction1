"""Numerical approximations of the original core, paired with directed errors.

The point calculation uses the production coupled radial recurrence. The
independent interval calculation retains every true pressure atom and the
nonzero implicit swirl source. Its admitted analytic tail bounds the distance
from the point approximation to the original nonlinear core. No midpoint of
a source box is selected as an original value.
"""
from fractions import Fraction
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "experiments/root_st073"
PREFIX = "lei_ren_part1_paper_"


def _box(ctx, raw):
    return ctx.mpf([mp.make_mpf(tuple(raw[side+"_exact_mpf_tuple"])) for side in ("lower", "upper")])


def _bounds(ctx, box):
    return tuple(ctx.make_mpf(row) for row in box._mpi_)


def _exact_fraction(raw):
    sign,mantissa,exponent,_ = raw
    value = Fraction((-1 if sign else 1)*mantissa)
    return value*2**exponent if exponent >= 0 else value/Fraction(2**(-exponent))


def _chebyshev_value(ctx, coefficients, x):
    b1 = b2 = ctx.mpf(0)
    for coefficient in reversed(coefficients[1:]):
        b0 = coefficient+2*x*b1-b2
        b2,b1 = b1,b0
    return coefficients[0]+x*b1-b2


def original_early_measure(ctx, nodes=128):
    """Compute the original integral, using a numerical cutoff primitive.

    The polynomial is quadrature machinery, not a new defining pressure.
    Exact-original error is obtained from the directed source enclosure.
    """
    if type(nodes) is not int or not 32 <= nodes <= 512:
        raise ValueError("Use 32..512 Chebyshev panels")
    sigma = []
    for j in range(nodes+1):
        y = (1+ctx.cos(ctx.pi*j/nodes))/2
        if j in (0,nodes):
            value = ctx.mpf(1 if j == 0 else 0)
        else:
            # The same exp(-1/y^2) flat cutoff as the original schedule.
            value = 1/(1+ctx.exp(1/y**2-1/(1-y)**2))
        sigma.append(value)
    coefficients = []
    for k in range(nodes+1):
        total = (sigma[0]+(-1)**k*sigma[-1])/2
        total += ctx.fsum(sigma[j]*ctx.cos(ctx.pi*k*j/nodes) for j in range(1,nodes))
        coefficients.append(total*2/nodes/(2 if k in (0,nodes) else 1))
    primitive = [ctx.mpf(0)]*(nodes+2)
    primitive[1] += coefficients[0]/2
    primitive[2] += coefficients[1]/8
    for k in range(2,nodes+1):
        primitive[k+1] += coefficients[k]/(4*(k+1))
        primitive[k-1] -= coefficients[k]/(4*(k-1))
    primitive[0] -= ctx.fsum(coefficient*(-1)**k for k,coefficient in enumerate(primitive))
    integrand = lambda y: ctx.exp(y/5-ctx.mpf(6)/5*_chebyshev_value(ctx,primitive,2*y-1))/2
    value = ctx.quad(integrand, [0,ctx.mpf(1)/4,ctx.mpf(1)/2,ctx.mpf(3)/4,1])
    return value


class OriginalCoreNumeric:
    """Compute nonzero normalized core approximations with same-source bounds.

    Physical F0 remains logarithmic/implicit. Output is the leading core
    profile, not n1, not an assembled spacetime velocity field or recursion.
    """
    def __init__(self, *, precision=60, quadrature_nodes=128, directory=None):
        if type(precision) is not int or not 35 <= precision <= 120:
            raise ValueError("Working precision must be 35..120 decimal digits")
        self.directory = Path(directory) if directory else DIRECTORY
        self.ctx = c = mp.mp.clone()
        c.dps = precision
        self.iv = iv = MPIntervalContext()
        iv.dps = precision+30
        self.precision = precision
        self.hashes = {}
        def read(stem):
            name = PREFIX+stem+".json"
            path = self.directory/name
            self.hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
            return json.loads(path.read_text(encoding="utf-8"))
        accepted = read("compliant_current_core_common_fixed_point_check")
        if not accepted.get("all_passed"):
            raise ValueError("Accepted original nonlinear fixed point required")
        for name,digest in accepted["input_hashes"].items():
            if hashlib.sha256((self.directory/name).read_bytes()).hexdigest() != digest:
                raise ValueError("Changed original numerical source: "+name)
            self.hashes[name] = digest
        self.family = {key:accepted[key] for key in (
            "actual_five_defect_family_sha256", "implicit_source_sha256", "datum_enclosure_sha256")}
        current_n1 = read("compliant_current_original_Rm_outer_mixed_source_joins_check")
        if current_n1["source_family"] != self.family:
            raise ValueError("Numerical leading core must share the current n1 source family")
        self.major = read("compliant_core_transfer")
        self.tube = read("shared_analytic_tube")
        self.linear = read("shared_linear_resolvent")
        fixed = read("shared_bump_constants")
        pressure = read("compliant_pressure_source")["compliant_source"]
        if (pressure["implicit_source_sha256"] != self.family["implicit_source_sha256"]
                or pressure["datum_enclosure_sha256"] != self.family["datum_enclosure_sha256"]
                or self.major["datum_enclosure_sha256"] != self.family["datum_enclosure_sha256"]
                or not self.major["contraction_proved"] or pressure["stage_count_total"] != 14):
            raise ValueError("Same complete pressure and nonlinear core family required")
        # The active minimum branch is exactly rational; recompute it instead
        # of selecting either endpoint or a midpoint of required_j.
        CA,CQ,CS,KN = (fixed[key] for key in ("CA","CQ","CS","KN"))
        if any(type(value) is not int or value <= 0 for value in (CA,CQ,CS,KN)):
            raise ValueError("Original fixed integer constants required")
        j = Fraction(1,16000000*KN*CS*CA)
        exact_j = iv.mpf(j.numerator)/j.denominator
        if (_bounds(c,exact_j)[1] >= _bounds(c,_box(iv,fixed["epsilon0"])/32)[0]
                or j >= Fraction(1,6400*CA*CA*CQ)):
            raise ValueError("Original minimum branch changed")
        lo,hi = (_exact_fraction(self.major["required_j"][side+"_exact_mpf_tuple"])
                 for side in ("lower","upper"))
        if not lo <= j <= hi:
            raise ValueError("Recomputed exact j differs from original source")
        self.j_fraction = j
        self.j = c.mpf(j.numerator)/j.denominator
        self.j_iv = exact_j
        self.logP = c.exp(40)+11
        self.logLambda = 4*self.logP+1000
        self.delta = c.exp(-4*self.logP-30)
        self.epsilon = c.exp(-self.logLambda)
        self.lp_iv = iv.exp(40)+11
        self.ll_iv = 4*self.lp_iv+1000
        self.delta_iv = iv.exp(-4*self.lp_iv-30)
        self.epsilon_iv = iv.exp(-self.ll_iv)
        if not 0 < _bounds(c,self.delta_iv)[1] < c.mpf("1e-200"):
            raise ValueError("Pinned Md40 source requires the exponential delta minimum branch")
        self.mass_point = original_early_measure(c,quadrature_nodes)
        # Extending the turnoff integral to infinity is an approximation.
        # The interval branch below retains its actual truncated value and
        # all nonzero late atoms; none are declared equal to zero.
        self.M2_point = c.mpf(5)/2+self.mass_point+c.exp(-c.mpf(2)/5)/2
        self.pressure = pressure
        if str(self.directory) not in sys.path:
            sys.path.insert(0,str(self.directory))
        self.taylor = importlib.import_module(PREFIX+"interval_taylor").IntervalTaylor
        self.q_jets = importlib.import_module(PREFIX+"candidate_pressure_function").q_jets
        self.advance = importlib.import_module(PREFIX+"logarithmic_core_step").advance_scaled_one
        self.initial_rows = importlib.import_module(PREFIX+"functional_core_step").initial_rows
        self.tail_factor = importlib.import_module(PREFIX+"analytic_radial_tail").tail_factor
        self.quadrature_nodes = quadrature_nodes

    def _seed(self, Z, degree, depth, *, directed):
        c = self.iv if directed else self.ctx
        z0 = c.mpf(Z)
        n = degree+depth
        variable = self.taylor.variable(c,z0,n)
        one = self.taylor.constant(c,1,n)
        j = self.j_iv if directed else self.j
        delta = self.delta_iv if directed else self.delta
        eps = self.epsilon_iv if directed else self.epsilon
        ll = self.ll_iv if directed else self.logLambda
        lp = self.lp_iv if directed else self.logP
        U = 4*variable+j
        L = one-variable*variable*delta
        H = variable*(1-delta)/2+(one-variable*variable)*U
        gradient = L*H/(H*H+(j/500)**2)
        qs = self.q_jets(c,z0,n)
        if directed:
            M2 = _box(c,self.pressure["fixed_beta2_mass_normalized"])
            M0 = _box(c,self.pressure["fixed_beta0_mass_normalized"])
            flat = _box(c,self.pressure["stages"]["z_flatten"]["mass"])
            bound = _box(c,self.pressure["flatten_complex_mass_upper"])
            radius = c.mpf(".25")
            pressure = []
            for k,q in enumerate(qs):
                remainder = flat if k == 0 else c.mpf(0) if Z == "0" and k%2 else (
                    c.mpf([-_bounds(self.ctx,bound/radius**k)[1], _bounds(self.ctx,bound/radius**k)[1]]))
                pressure.append(-(M2*q+(M0 if k == 0 else 0)+remainder))
            S_bound = c.exp(-6*ll-2000)
            eta = _box(c,self.tube["complex_tube_radius"])
            S = []
            for k in range(n+1):
                high = _bounds(self.ctx,S_bound/(eta/2)**k)[1]
                S.append(c.mpf([0 if k == 0 else -high, high]))
        else:
            pressure = [-self.M2_point*q for q in qs]
            S = [c.mpf(0)]*(n+1)  # numerical proxy; actual S is enclosed above
        fixed = dict(ell_Z_taylor=[-v for v in gradient.coefficients], S_Z_taylor=S,
            U0_Z_taylor=list(U.coefficients), P0_Z_taylor=[c.exp(2*lp-ll)*v for v in pressure])
        return fixed,self.initial_rows(c,fixed,z0,degree,required_depth=depth),delta,eps

    def coefficients(self, Z, *, degree=18, depth=1):
        if (type(Z) is not str or not -1 < self.ctx.mpf(Z) < 1
                or type(degree) is not int or not 4 <= degree <= 24
                or type(depth) is not int or not 1 <= depth <= 3):
            raise ValueError("Use an interior decimal Z string, degree4..24 and depth1..3")
        result = dict(Z=Z,degree=degree,depth=depth)
        for directed,label in ((False,"point"),(True,"directed")):
            c = self.iv if directed else self.ctx
            fixed,rows,delta,eps = self._seed(Z,degree,depth,directed=directed)
            for n in range(degree):
                self.advance(c,fixed,rows,n,c.mpf(Z),delta,eps)
            result[label] = rows
        return result

    def evaluate(self, packet, rho):
        c,iv = self.ctx,self.iv
        if type(rho) is not str or not 0 <= c.mpf(rho) <= c.mpf("4.1"):
            raise ValueError("Use a decimal rho string in [0,4.1]")
        r = c.mpf(rho)
        ri = iv.mpf(rho)
        N = packet["degree"]
        eps,epsi = self.epsilon,self.epsilon_iv
        h = _box(iv,self.tube["Xh_parameter"])
        values = {}
        finite_points,truths = {},{}
        correction = (_box(iv,self.major["scaled_map_size_upper"])
            /(1-_box(iv,self.major["scaled_map_Lipschitz_upper"])))
        def summarize(point,truth,tail):
            point_text = c.nstr(point,self.precision)
            error = max(abs(v) for v in _bounds(c,truth-iv.mpf(point_text)))
            text = lambda v:c.nstr(v,min(45,self.precision))
            upper_text = lambda v:text(v*(1+c.mpf("1e-30")))
            return dict(approximation=point_text,absolute_error_upper=upper_text(error),
                raw_directed_error_exact_mpf_tuple=list(error._mpf_),
                displayed_error_has_relative_safety_padding="1e-30",
                error_applies_to_serialized_decimal_approximation=True,
                relative_error_upper=None if point == 0 else upper_text(error/abs(point)),
                source_enclosure=[text(v) for v in _bounds(c,truth)],
                source_enclosure_display_is_rounded=True,
                source_enclosure_exact_mpf_tuples=[list(v) for v in truth._mpi_],
                original_analytic_radial_tail_upper=upper_text(_bounds(c,tail)[1]))
        for field in ("Phi","Psi","mean_Psi"):
            def finite(rows, radius, epsilon, zero):
                key = "A" if field == "Phi" else "Uz"
                first = 0 if field == "Phi" else 1
                return sum((rows[key][n][0]*radius**n/(epsilon if field != "Phi" else 1)
                    /(n+1 if field == "mean_Psi" else 1) for n in range(first,N+1)),zero)
            point = finite(packet["point"],r,eps,c.mpf(0))
            finite_box = finite(packet["directed"],ri,epsi,iv.mpf(0))
            # Subtract the admitted analytic model before bounding the tail.
            # Psi's model is linear in rho, so every row above N belongs to
            # the small nonlinear correction. A full Psi-ball tail would
            # lose its enormous correlation with the pressure-driven model.
            if r == 0:
                tail = iv.mpf(0)
            else:
                tail = correction*self.tail_factor(iv,degree=N,radial_order=0,axial_order=0,radius=ri,h=h)["tail_per_Xh_norm"]
                if field == "Phi":
                    q = ri/2
                    tail += (q**(N+1)/(math.factorial(N+1)*math.factorial(N+2))
                        /(1-q/((N+2)*(N+3))))
            bound = _bounds(c,tail)[1]
            truth = finite_box+iv.mpf([-bound,bound])
            # Directed subtraction covers the serialized numerical proxy,
            # whether or not it lies in the actual finite-prefix enclosure.
            finite_points[field],truths[field] = point,truth
            values[field] = summarize(point,truth,tail)
        # Recover actual leading axial and radial-recovery profile values.
        # The mean's Z derivative is consumed at fixed original rho; epsilon
        # is a source constant, so the source FTC and recovery are unchanged.
        zm,ziv = c.mpf(packet["Z"]),iv.mpf(packet["Z"])
        U0,U0iv = 4*zm+self.j,4*ziv+self.j_iv
        U = U0+eps*finite_points["Psi"]
        M = U0+eps*finite_points["mean_Psi"]
        Uiv = U0iv+epsi*truths["Psi"]
        Miv = U0iv+epsi*truths["mean_Psi"]
        meanZ = sum((packet["point"]["Uz"][n][1]*r**n/(eps*(n+1)) for n in range(1,N+1)),c.mpf(0))
        meanZiv = sum((packet["directed"]["Uz"][n][1]*ri**n/(epsi*(n+1)) for n in range(1,N+1)),iv.mpf(0))
        ztail = (iv.mpf(0) if r == 0 else correction*self.tail_factor(iv,degree=N,
            radial_order=0,axial_order=1,radius=ri,h=h)["tail_per_Xh_norm"])
        high = _bounds(c,ztail)[1]
        meanZiv += iv.mpf([-high,high])
        Q = (2*zm*U-(1-self.delta)*zm*M-(1-zm*zm)*(4+eps*meanZ))/(1-self.delta*zm*zm)
        Qiv = (2*ziv*Uiv-(1-self.delta_iv)*ziv*Miv-(1-ziv*ziv)*(4+epsi*meanZiv))/(1-self.delta_iv*ziv*ziv)
        values["Uz"] = summarize(U,Uiv,iv.mpf(0))
        values["Q"] = summarize(Q,Qiv,iv.mpf(0))
        pressure_point = -self.M2_point/(1+zm*zm)**2
        pressure_truth = -(_box(iv,self.pressure["fixed_beta2_mass_normalized"])/(1+ziv*ziv)**2
            +_box(iv,self.pressure["fixed_beta0_mass_normalized"])
            +_box(iv,self.pressure["stages"]["z_flatten"]["mass"]))
        values["P0_over_Pstar_squared"] = summarize(pressure_point,pressure_truth,iv.mpf(0))
        return dict(Z=packet["Z"],rho=rho,radial_degree=N,profiles=values,
            U0=summarize(U0,U0iv,iv.mpf(0)),
            physical_axial_recovery="Uz=U0+epsilon*Psi; retain the sum factored to avoid cancellation",
            physical_swirl_recovery="Utheta=sqrt(2*epsilon*rho)*F0*Phi; original F0 still implicit",
            numerical_approximation_of_original_normalized_core=True,
            point_uses_midpoint_of_source_boxes=False, all_twenty_production_terms_called=True,
            implicit_swirl_approximated_by_zero_only_in_point_proxy=True,
            actual_nonzero_swirl_and_all14_pressure_atoms_retained_in_error_branch=True,
            full_original_physical_velocity_evaluated=False, temporal_recursion=False)

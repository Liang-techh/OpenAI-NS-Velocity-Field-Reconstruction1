"""Original anchored amplitude in numerical logarithms and resolved peak coordinates.

The selected log Cstar is an exact binary singleton. Neither exp(-logCstar)
nor exp(-Lambda*G) is materialized. Actual certified pole disks provide the
global primitive enclosure; exact-root coordinates resolve its narrow peak.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import DIRECTORY, ROOT, PREFIX, _box, _bounds, _exact_fraction


class OriginalCoreLogAmplitude:
    def __init__(self, *, precision=90, directory=None):
        if type(precision) is not int or not 60 <= precision <= 120:
            raise ValueError("Use 60..120 decimal digits")
        self.directory = Path(directory) if directory else DIRECTORY
        self.ctx = c = mp.mp.clone()
        c.dps = precision
        self.iv = iv = MPIntervalContext()
        iv.dps = precision+40
        self.precision = precision
        self.hashes = {}
        def read(stem):
            name = PREFIX+stem+".json"
            path = self.directory/name
            self.hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
            return json.loads(path.read_text(encoding="utf-8"))
        self.numeric = read("compliant_current_original_core_numeric_approximation")
        if not self.numeric["checks"]["passed"]:
            raise ValueError("Checked same-source numerical core required")
        self.family = self.numeric["source_family"]
        for name,digest in self.numeric["implementation_hashes"].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
                raise ValueError("Changed numerical core implementation: "+name)
        certificate = read("compliant_anchored_axis_amplitude_check")
        self.archive = read("compliant_anchored_axis_amplitude")
        for record in (certificate,self.archive):
            if any(record[key] != value for key,value in self.family.items()):
                raise ValueError("Same original amplitude/core source family required")
        if not (certificate["all_passed"] and certificate["actual_anchored_G_resolved_with_directed_error"]):
            raise ValueError("Accepted original anchored primitive required")
        for record in (self.numeric,certificate):
            for name,digest in record["input_hashes"].items():
                if (hashlib.sha256((self.directory/name).read_bytes()).hexdigest() != digest
                        or self.hashes.get(name,digest) != digest):
                    raise ValueError("Changed original amplitude dependency: "+name)
                self.hashes[name] = digest
        fixed = read("shared_bump_constants")
        CA,CQ,CS,KN = (fixed[key] for key in ("CA","CQ","CS","KN"))
        j = Fraction(1,16000000*KN*CS*CA)
        if str(j) != self.numeric["levels"][-1]["exact_original_j_fraction"]:
            raise ValueError("Exact original j recipe differs")
        bounds = (_exact_fraction(fixed["required_j"][side+"_exact_mpf_tuple"]) for side in ("lower","upper"))
        low,high = bounds
        if not low <= j <= high or j >= Fraction(1,6400*CA*CA*CQ):
            raise ValueError("Original active j branch not retained")
        self.j,self.j_iv = c.mpf(j.numerator)/j.denominator,iv.mpf(j.numerator)/j.denominator
        self.lp,self.lp_iv = c.exp(40)+11,iv.exp(40)+11
        self.ll,self.ll_iv = 4*self.lp+1000,4*self.lp_iv+1000
        self.delta,self.delta_iv = c.exp(-4*self.lp-30),iv.exp(-4*self.lp_iv-30)
        if not _bounds(c,self.delta_iv)[1] < c.mpf("1e-200"):
            raise ValueError("The pinned source requires the exponential delta minimum")
        self.eps,self.eps_iv = c.exp(-self.ll),iv.exp(-self.ll_iv)
        self.Lambda,self.Lambda_iv = c.exp(self.ll),iv.exp(self.ll_iv)
        self.sigma,self.sigma_iv = self.j/500,self.j_iv/500
        norm = read("compliant_physical_norm_family")
        selected = norm["selected_logCstar"]
        self.selected_logC_tuple = selected["lower_exact_mpf_tuple"]
        if self.selected_logC_tuple != selected["upper_exact_mpf_tuple"]:
            raise ValueError("Canonical selected logCstar must be an exact singleton")
        self.logC = c.make_mpf(tuple(self.selected_logC_tuple))
        self.logC_iv = _box(iv,selected)
        self.roots = self.archive["root_certificates"]
        if not all(self.roots[key] for key in ("unique_real_anchor_certified","all_six_poles_certified",
                "uniform_interval_parameter_root_certificates","G_nonnegative_on_real_domain")):
            raise ValueError("Original uniform anchor and six-pole certificates required")
        self.anchor_iv = _box(iv,self.roots["anchor_interval"])
        self.anchor = c.re(self._root(-self.j/((9-self.delta)/2),0))
        a = (9-self.delta)/2
        plus = [self._root(seed,c.j*self.sigma) for seed in (-c.sqrt(a/4),-self.j/a,c.sqrt(a/4))]
        self.poles = [q for root in plus for q in (root,c.conj(root))]
        self._core_backend = None
        self._core_rows = {}

    def _H(self, z, directed=False):
        j,de = (self.j_iv,self.delta_iv) if directed else (self.j,self.delta)
        return -4*z**3-j*z*z+(9-de)*z/2+j

    def _hp(self,z,directed=False):
        j,de = (self.j_iv,self.delta_iv) if directed else (self.j,self.delta)
        return -12*z*z-2*j*z+(9-de)/2

    def _L(self,z,directed=False):
        return 1-(self.delta_iv if directed else self.delta)*z*z

    def _root(self,seed,target):
        c = self.ctx
        z = c.mpc(seed)
        for _ in range(50):
            step = (self._H(z)-target)/self._hp(z)
            z -= step
            if abs(step) < c.eps*max(1,abs(z)):
                return z
        raise ArithmeticError("Original cubic nominal root did not converge")

    def _complex_box(self,raw):
        return self.iv.mpc(_box(self.iv,raw["real"]),_box(self.iv,raw["imag"]))

    def summarize(self,point,truth):
        c,iv = self.ctx,self.iv
        decimal = c.nstr(point,self.precision)
        error = max(abs(v) for v in _bounds(c,truth-iv.mpf(decimal)))
        return dict(approximation=decimal,
            absolute_error_upper=c.nstr(error*(1+c.mpf("1e-30")),45),
            raw_directed_error_exact_mpf_tuple=list(error._mpf_),
            source_enclosure_exact_mpf_tuples=[list(v) for v in truth._mpi_],
            error_applies_to_serialized_decimal=True)

    def global_amplitude(self,Z):
        c,iv = self.ctx,self.iv
        anchor = Z == "anchor"
        if type(Z) is not str or not anchor and not -1 <= c.mpf(Z) <= 1:
            raise ValueError("Use a decimal Z string in [-1,1] or the exact anchor")
        if anchor:
            point,truth = c.mpf(0),iv.mpf(0)
        else:
            z,zi = c.mpf(Z),iv.mpf(Z)
            point = c.mpc(0)
            for root in self.poles:
                point += self._L(root)/(2*self._hp(root))*(c.log(z-root)-c.log(self.anchor-root))
            point = c.re(point)
            total = iv.mpc(0)
            for row in self.roots["poles"]:
                root = self._complex_box(row["disk"])
                numerator,denominator = iv.mpc(zi,0)-root,iv.mpc(self.anchor_iv,0)-root
                for argument in (numerator,denominator):
                    low,high = _bounds(c,argument.imag)
                    if not (low > 0 or high < 0):
                        raise ArithmeticError("Certified log path crossed its imaginary half-plane")
                total += self._L(root,True)/(2*self._hp(root,True))*(iv.log(numerator)-iv.log(denominator))
            low,high = _bounds(c,total.real)
            ilow,ihigh = _bounds(c,total.imag)
            if not ilow <= 0 <= ihigh or high < 0:
                raise ArithmeticError("Original real nonnegative primitive enclosure failed")
            truth = iv.mpf([max(c.mpf(0),low),high])
        offset_point,offset_truth = -self.Lambda*point,-self.Lambda_iv*truth
        return dict(Z=Z,anchor=self.summarize(self.anchor,self.anchor_iv),
            G=self.summarize(point,truth),log_F0_relative_to_anchor=self.summarize(offset_point,offset_truth),
            positive_F0=dict(kind="positive_exponential_of_separate_terms",
                exact_anchor_log_term=dict(negative_of_exact_mpf_tuple=self.selected_logC_tuple),
                relative_log_term=self.summarize(offset_point,offset_truth)),
            log_S_definition="-2 logLambda - 2 selected_logCstar - 2 Lambda G",
            selected_logC_is_original_singleton=True, six_original_certified_poles_used=True,
            source_midpoints_selected=False, absolute_F0_materialized=False,
            physical_velocity_materialized=False, temporal_recursion=False)

    def resolved_peak(self,xi):
        """Resolve Z=a+b*xi, b=sigma/sqrt(Lambda), retaining exact H(a)=0."""
        c,iv = self.ctx,self.iv
        if type(xi) is not str or not -4 <= c.mpf(xi) <= 4:
            raise ValueError("Use a decimal xi string in [-4,4]")
        values = []
        for directed in (False,True):
            ctx = iv if directed else c
            a = self.anchor_iv if directed else self.anchor
            j,de = (self.j_iv,self.delta_iv) if directed else (self.j,self.delta)
            eps = self.eps_iv if directed else self.eps
            sigma = self.sigma_iv if directed else self.sigma
            s = ctx.mpf(xi)
            b = sigma*ctx.sqrt(eps)
            hp = (9-de)/2-12*a*a-2*j*a
            q1,La = -12*a-j,1-de*a*a
            coefficients = [La*hp,b*(La*q1-2*de*a*hp),
                b*b*(-4*La-2*de*a*q1-de*hp),
                b**3*(8*de*a-de*q1),4*de*b**4]
            polynomial = sum((v*s**(k+2)/(k+2) for k,v in enumerate(coefficients)),ctx.mpf(0))
            if directed:
                offset = b*4
                qchange = abs(q1)*offset+4*offset**2
                Lchange = abs(de)*(2*abs(a)*offset+offset**2)
                if _bounds(c,hp-qchange)[0] <= 0 or _bounds(c,La-Lchange)[0] <= 0:
                    raise ArithmeticError("Original peak numerator lost positivity")
                remainder = eps*(abs(La)+Lchange)*(abs(hp)+qchange)**3*abs(s)**4/4
                high = _bounds(c,remainder)[1]
                raw = polynomial-iv.mpf([0,high])
                lo,hi = _bounds(c,raw)
                K = iv.mpf([max(c.mpf(0),lo),hi]) if c.mpf(xi) else iv.mpf(0)
                values.append((K,iv.exp(-K),b,remainder))
            else:
                values.append((polynomial,c.exp(-polynomial),b))
        (point,relative,b), (truth,relative_truth,btruth,remainder) = values
        return dict(xi=xi,coordinate="Z=exact_shared_anchor + positive_width_b*xi",
            anchor=self.summarize(self.anchor,self.anchor_iv),
            width_b=self.summarize(b,btruth),root_offset_not_added_to_anchor=True,
            scaled_G_Lambda_times_G=self.summarize(point,truth),
            F0_relative_to_exact_anchor=self.summarize(relative,relative_truth),
            denominator_remainder_upper=c.nstr(_bounds(c,remainder)[1]*(1+c.mpf("1e-30")),45),
            exact_H_at_anchor_zero_used=True,all_original_numerator_terms_retained=True,
            exact_source_denominator_retained_in_remainder=True,
            absolute_log_F0_terms=dict(exact_anchor_negative_logC_tuple=self.selected_logC_tuple,
                finite_negative_Lambda_G=self.summarize(-point,-truth)),
            absolute_F0_materialized=False, measured_blowup_dynamics=False, temporal_recursion=False)

    def leading_velocity_terms(self,Z,rho,*,log_tau="-1",azimuth="0"):
        """Leading physical core velocity as signed exponential sums.

        Both coefficient and logarithmic uncertainties are explicit. This
        is neither a floating-point vector nor the n0+n1 assembled field.
        """
        c,iv = self.ctx,self.iv
        if (type(Z) is not str or not -1 < c.mpf(Z) < 1 or type(rho) is not str
                or not 0 <= c.mpf(rho) <= 4):
            raise ValueError("Original interior decimal core coordinates required")
        if type(log_tau) is not str or not c.isfinite(c.mpf(log_tau)) or c.mpf(log_tau) > 0:
            raise ValueError("Use finite log(tau)<=0")
        if type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)):
            raise ValueError("Use a finite decimal azimuth")
        matches = [p for p in self.numeric["levels"][-1]["samples"] if p["Z"] == Z and p["rho"] == rho]
        if len(matches) == 1:
            sample = matches[0]
            numeric_origin = "saved_same_source_numerical_approximation"
        else:
            # Continue the original numerical function away from the saved
            # demonstration points; do not turn a sample table into a field.
            if self._core_backend is None:
                from .original_core_numeric import OriginalCoreNumeric
                self._core_backend = OriginalCoreNumeric(precision=max(70,self.precision),
                    quadrature_nodes=192,directory=self.directory)
                if self._core_backend.family != self.family:
                    raise ValueError("On-demand core evaluator changed the source family")
            if Z not in self._core_rows:
                self._core_rows[Z] = self._core_backend.coefficients(Z,degree=18)
            sample = self._core_backend.evaluate(self._core_rows[Z],rho)
            numeric_origin = "on_demand_same_source_production_recurrence"
        z,r,phi = c.mpf(Z),c.mpf(rho),c.mpf(azimuth)
        zi,ri,phii = iv.mpf(Z),iv.mpf(rho),iv.mpf(azimuth)
        ell = (c.mpf(log_tau)-c.log(1-z*z))/2
        elliv = (iv.mpf(log_tau)-iv.ln(1-zi*zi))/2
        amplitude = self.global_amplitude(Z)
        def profile(name):
            row = sample["profiles"][name]
            box = iv.mpf([mp.make_mpf(tuple(q)) for q in row["source_enclosure_exact_mpf_tuples"]])
            return c.mpf(row["approximation"]),box
        Q,Qiv = profile("Q")
        Phi,Phiiv = profile("Phi")
        Uz,Uziv = profile("Uz")
        def term(point,box,logs):
            return dict(coefficient=self.summarize(point,box),log_scale_terms=logs)
        axial = term(Uz,Uziv,[self.summarize(-(1+self.delta)*ell,-(1+self.delta_iv)*elliv)])
        if not r:
            radial = swirl = term(c.mpf(0),iv.mpf(0),[])
            radius = term(c.mpf(0),iv.mpf(0),[])
        else:
            radial_log = (c.log(r)-self.ll-c.log(2))/2-ell
            radial_iv = (iv.ln(ri)-self.ll_iv-iv.ln(2))/2-elliv
            swirl_log = (c.log(r)-self.ll+c.log(2))/2-(1+self.delta)*ell
            swirl_iv = (iv.ln(ri)-self.ll_iv+iv.ln(2))/2-(1+self.delta_iv)*elliv
            radial = term(Q,Qiv,[self.summarize(radial_log,radial_iv)])
            swirl = term(Phi,Phiiv,[dict(negative_of_exact_mpf_tuple=self.selected_logC_tuple),
                amplitude["log_F0_relative_to_anchor"],self.summarize(swirl_log,swirl_iv)])
            log_radius = ell+(c.log(2*r)-self.ll)/2
            log_radius_iv = elliv+(iv.ln(2*ri)-self.ll_iv)/2
            radius = term(c.mpf(1),iv.mpf(1),[self.summarize(log_radius,log_radius_iv)])
        def rotate(packet,point,box):
            row = packet["coefficient"]
            bounds = iv.mpf([mp.make_mpf(tuple(q)) for q in row["source_enclosure_exact_mpf_tuples"]])
            return term(c.mpf(row["approximation"])*point,bounds*box,packet["log_scale_terms"])
        return dict(Z=Z,rho=rho,log_tau=log_tau,azimuth=azimuth,
            coordinates=dict(r=radius,x=rotate(radius,c.cos(phi),iv.cos(phii)),
                y=rotate(radius,c.sin(phi),iv.sin(phii)),
                z=term(z,zi,[self.summarize((1-self.delta)*ell,(1-self.delta_iv)*elliv)]),
                t=dict(kind="one_minus_positive_exponential",exact_log_tau=log_tau)),
            cylindrical=dict(ur=radial,utheta=swirl,uz=axial),
            Cartesian=dict(u=[rotate(radial,c.cos(phi),iv.cos(phii)),rotate(swirl,-c.sin(phi),-iv.sin(phii))],
                v=[rotate(radial,c.sin(phi),iv.sin(phii)),rotate(swirl,c.cos(phi),iv.cos(phii))],w=[axial]),
            value_semantics="sum coefficient*exp(sum log_scale_terms), with original source enclosures",
            source_family=self.family,leading_core_only=True,coefficient_orders_included=[0],
            numerical_core_origin=numeric_origin,native_core_domain="0<=rho<=4, |Z|<1, log(tau)<=0",
            structurally_shared_mean_recovery=True,absolute_exponentials_materialized=False,
            numerical_float_velocity_vector_available=False,global_field_available=False,
            temporal_recursion=False,full_corrected_NS_solution=False)

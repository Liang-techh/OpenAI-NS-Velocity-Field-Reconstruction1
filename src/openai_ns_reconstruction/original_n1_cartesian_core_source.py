"""Axis-regular Cartesian n1 core sources; no original point oracle is implied."""
from .original_n1_physical_source import OriginalN1PhysicalSource
from .original_n1_source_dispatcher import _rational


class OriginalN1CartesianCoreSource(OriginalN1PhysicalSource):
    """Use physical x,y and paper Z,tau, including the original regular axis.

    The disk restriction rho=(x*x+y*y)/(2*lambda**2*epsilon)<=4 is part
    of the function domain. This class does not select an outer chart or
    invert physical z to Z. All original profile inputs remain unevaluated.
    """

    def __init__(self, source=None):
        super().__init__(source)
        inner = self.source._inner_inputs["rin"]["inputs"]
        self.core = inner["inner_cutoff_bindings"]["fields"]["core"]
        self.volterra = inner["inner_Volterra_bindings"]
        self.rho_variable = self.volterra["coordinates"]["core"]
        self.epsilon = self.volterra["parameters"]["epsilon"]
        self.axis_zero_sources = set()
        self.axis_source_evidence = []
        for quantity in ("F1", "Uz1", "K1", "P1"):
            for order, node in enumerate(self.volterra["functions"]["core"][quantity]):
                row = self.graph.nodes[node]
                expected = dict(operation="current_original_unique_regular_n1_component",
                    quantity=quantity, chart="core", coordinate=self.rho_variable,
                    Z_variable=self.Z_variable, coordinate_domain_lower=0,
                    ordinary_slow_Z_derivative_order=order, zero_axis_regular_power=2,
                    unique_regular_analytic_solution=True,
                    all_words_retained_in_analytic_limit=True)
                if any(row.get(key) != value for key, value in expected.items()):
                    raise ValueError("Same original zero-axis regular n1 solution required")
                self.axis_zero_sources.add(node)
                self.axis_source_evidence.append(dict(quantity=quantity, Z_order=order,
                    node=node, zero_axis_power_in_sqrt_rho=2))
        for quantity in ("F1", "Uz1", "K1", "P1"):
            actual = (self.core["pressure_P1"] if quantity == "P1"
                      else self.core["rows"][quantity][0])
            if actual != self.volterra["functions"]["core"][quantity]:
                raise ValueError("Cutoff core must retain the same guarded regular source: " + quantity)
        self.axis_derived_zeros = {}
        for quantity, rows in (("mean1", self.core["rows"]["mean1"][0]),
                               ("Q1", self.core["Q1"][0])):
            for order, node in enumerate(rows):
                if not self._zero_axis_algebra(node):
                    raise ValueError("Unproved original axis limit: " + quantity)
                self.axis_derived_zeros[f"{quantity}_Z{order}"] = node
        self.axis_derivative_evidence = {f"{quantity}_rho{order}": self._axis_derivative(quantity, order)
            for quantity in ("F1", "Uz1", "mean1", "K1") for order in (1, 2)}
        self.axial_rho_axis = self.axis_derivative_evidence["Uz1_rho1"]

    def _zero_axis_algebra(self, node):
        if node in self.axis_zero_sources or node == 0:
            return True
        row = self.graph.nodes[node]
        op = row["operation"]
        if op == "sum":
            return all(self._zero_axis_algebra(i) for i in row["arguments"])
        if op == "product":
            return any(self._zero_axis_algebra(i) for i in row["arguments"])
        if op == "negative":
            return self._zero_axis_algebra(row["argument"])
        if op == "positive_quotient":
            return self._zero_axis_algebra(row["numerator"])
        return False

    def _axis_derivative(self, quantity, order):
        node = self.core["rows"][quantity][order][0]
        row = self.graph.nodes[node]
        if (row.get("operation") != "regular_even_radial_axis_continuation"
                or row.get("rho") != self.rho_variable
                or row.get("quantity") != quantity
                or row.get("ordinary_rho_derivative_order") != order
                or not row.get("axis_value_is_exact_equation_limit")):
            raise ValueError("Original equation-derived even axis derivative required")
        return row["axis_expression"]

    def cartesian_core_function(self, *, tau=None, Z=None):
        """Export x,y,Z,tau on the original closed core disk, including x=y=0."""
        return self._function(tau=tau, Z=Z, axis=False)

    def axis_function(self, *, tau=None, Z=None):
        """Exact n1 axis values and the original transverse axial curvature."""
        return self._function(tau=tau, Z=Z, axis=True)

    def _function(self, *, tau, Z, axis):
        qt = None if tau is None else _rational(tau)
        qz = None if Z is None else _rational(Z)
        if qt is not None and not 0 < qt <= 1:
            raise ValueError("Normalized time requires 0 < tau <= 1")
        if qz is not None and not -1 < qz < 1:
            raise ValueError("Original interior requires -1 < Z < 1")
        g = self.graph
        tau_node = g.node("bound_variable", name="original_physical_tau") if qt is None else g.constant(qt)
        zeta = self.Z_variable if qz is None else g.constant(qz)
        bind_z = lambda node: node if qz is None else g.at(node, self.Z_variable, zeta)
        d = g.add(1, g.neg(g.mul(zeta, zeta)))
        lam = g.unary("sqrt", g.quotient(tau_node, d, "tau>0 and |Z|<1"))
        power = lambda exponent: g.unary("exp", g.mul(exponent, g.unary("log", lam)))
        x = 0 if axis else g.node("bound_variable", name="original_physical_x")
        y = 0 if axis else g.node("bound_variable", name="original_physical_y")
        R = g.quotient(g.add(g.mul(x, x), g.mul(y, y)), g.mul(g.constant(2), lam, lam), "lambda>0")
        rho = g.quotient(R, self.epsilon, "same original epsilon>0")
        raw = dict(F1=self.core["rows"]["F1"][0][0], Uz1=self.core["rows"]["Uz1"][0][0],
            mean1=self.core["rows"]["mean1"][0][0], Q1=self.core["Q1"][0][0],
            P1=self.core["pressure_P1"][0])
        # rho depends on Z at fixed physical x,y,tau. Do not reuse the native
        # query helper's Z-independent-substitution annotation for this map.
        profiles = {key: 0 if axis else bind_z(g.node("function_substitution",
            expression=node, variable=self.rho_variable, value=rho,
            Z_independent_substitution=False, physical_coordinate_pullback=True)) for key, node in raw.items()}
        de = self.delta
        radial = g.mul(g.constant("1/2"), power(g.add(g.constant(-2), g.mul(g.constant(2), de))), profiles["Q1"])
        swirl = g.mul(power(g.add(g.constant(-2), de)), profiles["F1"])
        u = g.add(g.mul(radial, x), g.neg(g.mul(swirl, y)))
        v = g.add(g.mul(radial, y), g.mul(swirl, x))
        w = g.mul(power(g.add(g.constant(-1), de)), profiles["Uz1"])
        p = g.mul(power(g.constant(-2)), profiles["P1"])
        psi = g.mul(power(g.add(1, de)), R, profiles["mean1"])
        curvature = g.mul(power(g.add(g.constant(-3), de)),
                          g.quotient(bind_z(self.axial_rho_axis), self.epsilon, "original epsilon>0"))
        return dict(chart="inner_core_axis" if axis else "inner_core_cartesian",
            coordinates=dict(x=x, y=y, z=g.mul(power(g.add(1, g.neg(de))), zeta),
                t=g.add(1, g.neg(tau_node)), R=R, rho=rho, Z=zeta, tau=tau_node, lambda_source=lam),
            source_inputs=dict(raw_core=raw, bound_profiles=profiles, delta=de,
                epsilon=self.epsilon, original_axial_rho_axis=self.axial_rho_axis),
            transverse_velocity_factors=dict(ur_over_r=radial, utheta_over_r=swirl),
            velocity=dict(u=u, v=v, w=w), pressure=p, streamfunction=psi,
            axis_limits=dict(u=0, v=0, w=0, pressure=0, streamfunction=0,
                all_first_Cartesian_velocity_derivatives=0,
                axial_transverse_hessian=dict(xx=curvature, yy=curvature, xy=0)),
            parameter_domain=dict(tau="0 < tau <= 1", Z="-1 < Z < 1",
                spatial="x=y=0" if axis else "0 <= rho=(x*x+y*y)/(2*lambda^2*epsilon) <= 4"),
            definition=dict(coefficient_order=1, normalized_viscosity=1, normalized_blowup_time=1,
                source_only=True, axis_included=True, physical_xy_coordinates=True,
                transverse_formula_has_no_radius_division=True, no_azimuth_required=True,
                same_original_streamfunction_recovery=True,
                axis_values_from_original_regular_solution=True,
                free_time_variable=qt is None, source_domain_must_be_honored=True,
                physical_z_inverse_available=False, automatic_outer_chart_selection=False,
                includes_leading_background=False, includes_higher_temporal_coefficients=False,
                original_numeric_field_values_evaluated=False, numerical_oracle_available=False,
                globally_joined_physical_field=False, full_corrected_Navier_Stokes_solution=False))

    def snapshot(self, packet):
        payload = super().snapshot(packet)
        payload.update(axis_source_evidence=self.axis_source_evidence,
                       axis_derived_zero_sources=self.axis_derived_zeros,
                       axis_derivative_evidence=self.axis_derivative_evidence)
        return payload

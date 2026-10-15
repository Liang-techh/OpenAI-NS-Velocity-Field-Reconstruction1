"""Regional physical n=1 velocity/pressure sources in paper coordinates.

This exports an unevaluated parametric source field, with viscosity normalized
to one. It is a single temporal coefficient, not the complete corrected NS
solution. Automatic inversion from arbitrary Cartesian points is separate.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
from types import SimpleNamespace

from .original_n1_source_dispatcher import OriginalN1SourceDispatcher, _rational, _sha


class OriginalN1PhysicalSource:
    """Recover radial/axial velocity from one original streamfunction."""
    def __init__(self, source=None):
        self.source = source if source is not None else OriginalN1SourceDispatcher().extend_inner_high_sources()
        self.graph = self.source.graph
        self.delta = self.source._inner_inputs["rin"]["inputs"]["inner_Volterra_bindings"]["parameters"]["delta"]
        component = self.source._inner_inputs["rin"]["inputs"]["inner_Volterra_bindings"]["functions"]["core"]["F1"][0]
        self.Z_variable = self.graph.nodes[component]["Z_variable"]
        self.prefix_nodes = (self.source.source_construction_end
                             if getattr(self.source,"inner_high_source_restored",False) else self.graph.archived_nodes)
        self.mass_identity_evidence = []
        for name, chart in self.source.charts.items():
            for field,moment in zip(chart["fields"]["M1"][0],chart["cumulative"]["m11"][0]):
                method = "identical_source_node"
                if field != moment:
                    row = self.graph.nodes[field]
                    if (row["operation"] != "function_substitution" or row["expression"] != moment
                            or self._may_depend(moment,row["variable"])):
                        raise ValueError("The same original cumulative mass source is required")
                    method = "capture_free_substitution_of_same_terminal_mass"
                self.mass_identity_evidence.append(dict(chart=name,field_node=field,moment_node=moment,method=method))

    def _may_depend(self, node, variable):
        """Conservative occurrence proof on the append-only original DAG.

        A previously created source cannot reference a later bound variable.
        This proves the tail coordinate translation leaves terminal mass fixed.
        Unknown later operations fail closed rather than erasing dependence.
        """
        if node == variable:
            return True
        if node < variable:
            return False
        row = self.graph.nodes[node]
        op = row["operation"]
        if op in ("sum","product"):
            args = row["arguments"]
        elif op in ("negative","analytic_unary"):
            args = [row["argument"]]
        elif op == "positive_quotient":
            args = [row["numerator"],row["denominator"]]
        elif op == "function_substitution":
            args = ([row["expression"]] if row["variable"] != variable else [])+[row["value"]]
        elif op == "definite_integral":
            args = [row["integrand"],row["lower"],row["upper"]]
        elif op in ("exact_rational","bound_variable"):
            return False
        else:
            raise ValueError("Unproved terminal-mass variable independence")
        return any(self._may_depend(i,variable) for i in args)

    def regional_function(self,chart,*,tau=None,Z=None,azimuth=None):
        """Export a whole native interval with free spatial/time parameters."""
        return self.query(chart,tau=tau,Z=Z,azimuth=azimuth,_regional_function=True)

    def query(self, chart, *, fraction=None, tail_log_offset=None, tau="1/2", Z=None, azimuth=None,
              _regional_function=False):
        """Use native chart parameters, tau=1-t, |Z|<1 and angle in radians.

        Omit Z or azimuth to keep that variable free in the returned source
        expressions. Exact rational parameters avoid hidden point rounding.
        """
        qtau = None if tau is None and _regional_function else _rational(tau)
        if qtau is not None and not 0 < qtau <= 1:
            raise ValueError("Normalized time requires 0 < tau=1-t <= 1")
        qZ = None if Z is None else _rational(Z)
        if qZ is not None and not -1 < qZ < 1:
            raise ValueError("Physical interior requires -1 < Z < 1")
        qphi = None if azimuth is None else _rational(azimuth)
        # Domain/axis checks are performed before adding physical source nodes.
        if _regional_function:
            if fraction is not None or tail_log_offset is not None or chart not in self.source.charts:
                raise ValueError("Choose a known complete original native chart")
            region = self.source.charts[chart]
            def query(field,physical_y_order=0,ordinary_Z_order=0):
                return SimpleNamespace(bound_node=region["fields"][field][physical_y_order][ordinary_Z_order],
                    radius_node=region["radius"],exact_coordinate_parameter=None)
        else:
            def query(field,**orders):
                return self.source.query(chart,field,fraction=fraction,tail_log_offset=tail_log_offset,**orders)
        mass = query("M1")
        my = query("M1",physical_y_order=1)
        mz = query("M1",ordinary_Z_order=1)
        theta = query("Utheta1")
        axial = query("Uz1")
        auxiliary = query("Q1")
        pressure = query("P1")
        g = self.graph
        zeta = self.Z_variable if qZ is None else g.constant(qZ)
        angle = (g.node("bound_variable",name="original_physical_azimuth") if qphi is None else g.constant(qphi))
        bind = lambda i: i if qZ is None else g.at(i,self.Z_variable,zeta)
        R = bind(mass.radius_node)
        M, My, MZ, Ut, Uz, Q, P = (bind(q.bound_node) for q in (mass,my,mz,theta,axial,auxiliary,pressure))
        tau_node = (g.node("bound_variable",name="original_physical_tau") if qtau is None else g.constant(qtau))
        de = self.delta
        d = g.add(1,g.neg(g.mul(zeta,zeta)))
        L = g.add(1,g.neg(g.mul(de,zeta,zeta)))
        lam = g.unary("sqrt",g.quotient(tau_node,d,"tau>0 and |Z|<1"))
        power = lambda exponent: g.unary("exp",g.mul(exponent,g.unary("log",lam)))
        one_plus = g.add(1,de)
        radial_square_root = g.unary("sqrt",g.mul(g.constant(2),R))
        r = g.mul(lam,radial_square_root)
        z = g.mul(power(g.add(1,g.neg(de))),zeta)
        cos, sin = g.unary("cos",angle),g.unary("sin",angle)
        axial_weight = power(g.add(g.constant(-1),de))
        radial_weight = power(g.add(g.constant(-1),g.mul(g.constant(2),de)))
        psi = g.mul(power(one_plus),M)
        # psi=lambda^(1+delta) M1. Using the actual M_y and M_Z rows
        # makes both meridional components derivatives of this same source.
        radial_numerator = g.add(g.mul(g.constant(2),zeta,My),
            g.neg(g.mul(one_plus,zeta,M)),g.neg(g.mul(d,MZ)))
        ur = g.quotient(g.mul(power(g.mul(g.constant(2),de)),radial_numerator),
                        g.mul(r,L),"r>0, 0<delta<1 and |Z|<1")
        uz = g.mul(axial_weight,g.quotient(My,R,"original positive regional R"))
        ut = g.mul(axial_weight,Ut)
        u = g.add(g.mul(ur,cos),g.neg(g.mul(ut,sin)))
        v = g.add(g.mul(ur,sin),g.mul(ut,cos))
        p = g.mul(power(g.constant(-2)),P)
        return dict(chart=chart,exact_native_parameter=mass.exact_coordinate_parameter,
            exact_tau=None if qtau is None else str(qtau),exact_Z=None if qZ is None else str(qZ),
            exact_azimuth=None if qphi is None else str(qphi),
            coordinates=dict(x=g.mul(r,cos),y=g.mul(r,sin),z=z,t=g.add(1,g.neg(tau_node)),
                             r=r,R=R,Z=zeta,lambda_source=lam,tau=tau_node,azimuth=angle),
            source_inputs=dict(M1=M,M1_y=My,M1_Z=MZ,Utheta1=Ut,Uz1=Uz,Q1=Q,P1=P,delta=de),
            streamfunction=psi,velocity=dict(u=u,v=v,w=uz,ur=ur,utheta=ut,uz=uz),pressure=p,
            parameter_domain=dict(native_coordinate=self.source.charts[chart]["coordinate"],
                native_domain=self.source.charts[chart]["domain"],
                native_lower_open=self.source.charts[chart]["lower_open"],
                tau="0 < tau <= 1",Z="-1 < Z < 1",azimuth="all real angles (radians)"),
            direct_profile_recovery=dict(uz=g.mul(axial_weight,Uz),
                ur=g.mul(radial_weight,g.unary("sqrt",g.mul(g.constant(Fraction(1,2)),R)),Q)),
            definition=dict(coefficient_order=1,normalized_viscosity=1,normalized_blowup_time=1,
                source_only=True,parametric_physical_coordinates=True,
                complete_regional_spatial_source=_regional_function,free_time_variable=qtau is None,
                M1_is_original_cumulative_m11=True,
                radial_axial_velocity_from_one_streamfunction=True,
                divergence_structure="u_r=-(1/r) psi_z; u_z=(1/r) psi_r",
                includes_leading_background=False,includes_higher_temporal_coefficients=False,
                includes_oscillatory_correction=False,axis_included=False,
                inverse_Cartesian_chart_selection_available=False,numerical_oracle_available=False,
                globally_joined_physical_field=False,full_corrected_Navier_Stokes_solution=False))

    def snapshot(self, packet):
        """Small query suffix relative to the checked original source chain."""
        base = self.source.source_archives["Rm_outer_mixed_source_joins"]
        chain = [dict(archive=base["archive"],sha256=_sha(self.source.directory/base["archive"]))]
        if getattr(self.source,"inner_high_source_restored",False):
            name = "lei_ren_part1_paper_compliant_current_original_n1_inner_high_source_functions.json.gz"
            chain.append(dict(archive=name,sha256=_sha(self.source.directory/name)))
        return dict(source_family=self.source.source_family,source_graph_chain=chain,
            mass_identity_evidence=self.mass_identity_evidence,
            source_archives=self.source.source_archives,
            inner_high_source_restored=getattr(self.source,"inner_high_source_restored",False),
            graph_prefix_nodes=self.prefix_nodes,query_graph_suffix=self.graph.nodes[self.prefix_nodes:],
            physical_source=packet)

    @classmethod
    def restore_export(cls,path,source=None):
        """Restore checked physical source roots and their small query suffix."""
        path = Path(path)
        receipt_path = path.with_name(path.name.removesuffix(".json.gz")+"_check.json")
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if not receipt["checks"]["passed"] or _sha(path) != receipt["archive_sha256"]:
            raise ValueError("Unchanged checked physical source export required")
        root = Path(__file__).resolve().parents[2]
        for name,digest in receipt["implementation_hashes"].items():
            if _sha(root/name) != digest:
                raise ValueError("Changed physical source implementation: "+name)
        with gzip.open(path,"rt",encoding="utf-8") as stream:
            raw = json.load(stream)
        adapter = cls(source)
        if (raw["source_family"] != adapter.source.source_family
                or raw["implementation_hashes"] != receipt["implementation_hashes"]
                or raw["graph_prefix_nodes"] != adapter.prefix_nodes
                or len(adapter.graph.nodes) != adapter.prefix_nodes
                or raw["mass_identity_evidence"] != adapter.mass_identity_evidence):
            raise ValueError("Same checked original source prefix required before physical queries")
        for item in raw["source_graph_chain"]:
            if _sha(adapter.source.directory/item["archive"]) != item["sha256"]:
                raise ValueError("Changed original physical source graph chain")
        adapter.graph.nodes.extend(raw["query_graph_suffix"])
        for i in range(adapter.prefix_nodes,len(adapter.graph.nodes)):
            adapter.graph._keys[adapter.graph._key(adapter.graph.nodes[i])] = i
        return adapter,raw["physical_source"]

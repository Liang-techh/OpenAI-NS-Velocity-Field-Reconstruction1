"""Query accepted original n=1 regional source DAGs without replaying owners.

These are source expressions, not evaluated velocity samples. In particular,
the inner finite derivative contract is weaker than the outer one. The axis,
numerical oracles and a uniform global high derivative contract remain open.
"""
from dataclasses import asdict, dataclass
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path


PREFIX = "lei_ren_part1_paper_compliant_current_original_"
INPUTS = (
    ("Rm_outer_mixed_source_joins", "Rm_outer_mixed_source_join_bindings"),
    ("prescribed_high_radial_source_functions", "prescribed_high_radial_source_bindings"),
    ("n1_combined_outer_mixed_source_functions", "combined_outer_mixed_source_bindings"),
    ("n1_Rin_source_join", "Rin_source_join_bindings"),
    ("n1_intermediate_own_C3_chart_functions", "intermediate_own_C3_chart_bindings"),
    ("lower_higher_radial_Rin_join", "lower_higher_radial_Rin_bindings"),
)
STATE = ("Uz1", "Utheta1", "F1", "M1", "mean1", "K1", "Q1", "V1")
MOMENTS = ("m11", "m12", "m13", "m14", "m15")


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise TypeError("Use an exact rational int, string or Fraction")
    return Fraction(value)


def _sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


class _QueryGraph:
    """An immutable archived prefix with a separate append-only query suffix.

    Only the small suffix is interned. Identical expressions in the archived
    prefix may acquire a new query handle; original node numbers never change.
    """
    def __init__(self, nodes):
        self.nodes = nodes
        self.archived_nodes = len(nodes)
        self._keys = {}
        for i, value in enumerate((0, 1)):
            row = dict(operation="exact_rational", numerator=value, denominator=1)
            if nodes[i] != row:
                raise ValueError("Original exact zero/one prefix required")
            self._keys[self._key(row)] = i

    @staticmethod
    def _key(row):
        return json.dumps(row, sort_keys=True, separators=(",", ":"))

    def node(self, operation, **attributes):
        row = dict(operation=operation, **attributes)
        key = self._key(row)
        if key not in self._keys:
            self._keys[key] = len(self.nodes)
            self.nodes.append(row)
        return self._keys[key]

    def constant(self, value):
        q = Fraction(value)
        return self.node("exact_rational", numerator=q.numerator, denominator=q.denominator)

    def add(self, *arguments):
        ids = [i for i in arguments if i != 0]
        return 0 if not ids else ids[0] if len(ids) == 1 else self.node("sum", arguments=ids)

    def mul(self, *arguments):
        if 0 in arguments:
            return 0
        ids = [i for i in arguments if i != 1]
        return 1 if not ids else ids[0] if len(ids) == 1 else self.node("product", arguments=ids)

    def neg(self, value):
        return 0 if value == 0 else self.node("negative", argument=value)

    def unary(self, name, argument):
        return self.node("analytic_unary", name=name, argument=argument)

    def at(self, expression, variable, value):
        return self.node("function_substitution", expression=expression, variable=variable,
                         value=value, Z_independent_substitution=True)


@dataclass(frozen=True)
class SourceQuery:
    chart: str
    group: str
    field: str
    physical_y_order: int
    ordinary_Z_order: int
    exact_coordinate_parameter: str
    source_node: int
    bound_node: int
    native_query_node: int
    radius_node: int
    numerical_oracle_available: bool = False

    def as_dict(self):
        return asdict(self)


class OriginalN1SourceDispatcher:
    """One entry point for regional own fields and all five cumulative moments.

    A finite chart is addressed by a fraction of its *original native domain*.
    The tail is addressed by log(R/Rv). Physical y derivatives are not native
    derivatives. Unsupported derivative requests fail before adding graph nodes.
    """
    def __init__(self, source_directory=None):
        self.directory = Path(source_directory) if source_directory else (
            Path(__file__).resolve().parents[2] / "experiments" / "root_st073")
        self.charts = {}
        self.input_hashes = {}
        self.source_archives = {}
        self.source_family = None
        self.graph = None
        packets = {}
        for stem, binding in INPUTS:
            name = PREFIX + stem + ".json.gz"
            receipt_name = PREFIX + stem + "_check.json"
            receipt = json.loads((self.directory / receipt_name).read_text(encoding="utf-8"))
            if not receipt.get("all_passed") or not receipt.get("checks") or not all(
                proof.get("passed") for proof in receipt["checks"].values()
            ):
                raise ValueError("Accepted source receipt required: " + receipt_name)
            pins = dict(receipt["input_hashes"])
            if pins.get(name) != _sha(self.directory / name):
                raise ValueError("Changed accepted archive: " + name)
            for pinned_name, digest in pins.items():
                if pinned_name in self.input_hashes and self.input_hashes[pinned_name] != digest:
                    raise ValueError("Inconsistent ancestor pin: " + pinned_name)
                self.input_hashes[pinned_name] = digest
            self.input_hashes[receipt_name] = _sha(self.directory / receipt_name)
            with gzip.open(self.directory / name, "rt", encoding="utf-8") as stream:
                raw = json.load(stream)
            nodes = raw["exact_graph_nodes"]
            recorded_count = receipt.get("exact_graph_nodes")
            if (recorded_count is not None and len(nodes) != recorded_count
                    or raw["source_family"] != receipt["source_family"]):
                raise ValueError("Archive/receipt identity mismatch: " + name)
            if self.graph is None:
                self.graph = _QueryGraph(nodes)
                self.source_family = raw["source_family"]
            elif raw["source_family"] != self.source_family or len(nodes) > self.graph.archived_nodes or any(
                row != self.graph.nodes[i] for i, row in enumerate(nodes)
            ):
                raise ValueError("Source is not a prefix of the current accepted graph: " + name)
            expected_pins = dict(pins)
            expected_pins.pop(name)
            expected_pins.pop(PREFIX + stem + "_check.py")
            if expected_pins != raw["input_hashes"]:
                raise ValueError("Archive/receipt input pins disagree: " + name)
            packets[stem] = raw[binding]
            self.source_archives[stem] = dict(archive=name, receipt=receipt_name, binding=binding)
            del raw, nodes
        # Integrity checks only; no ancestor mathematical suite is executed.
        for name, digest in self.input_hashes.items():
            if _sha(self.directory / name) != digest:
                raise ValueError("Changed accepted source input: " + name)
        for stem in ("prescribed_high_radial_source_functions", "Rm_outer_mixed_source_joins"):
            for name, chart in packets[stem]["charts"].items():
                self._install(name, chart, chart["own"], stem)
        stem = "n1_combined_outer_mixed_source_functions"
        for name, chart in packets[stem]["charts"].items():
            self._install(name, chart, chart, stem)
        rin = packets["n1_Rin_source_join"]
        lower = packets["n1_intermediate_own_C3_chart_functions"]["charts"]["lower"]
        self._install("inner_first", dict(coordinate=rin["coordinate"], radius=rin["radius"],
                      domain=[0, rin["endpoint"]]), rin["left"], "n1_Rin_source_join")
        self._install("lower", lower, rin["right"], "n1_Rin_source_join")
        self._install_inner_core(rin)
        high = packets["lower_higher_radial_Rin_join"]
        self._install("inner_flat", dict(coordinate=rin["coordinate"], radius=rin["radius"],
                      domain=high["flat_domain"]), high["left"], "lower_higher_radial_Rin_join")
        self.coverage = dict(regional_source_dispatcher=True, original_preRc_charts=17,
            axis_included=False, axis_regular_Cartesian_recovery=False,
            uniform_global_y4_pressure_y3_contract=False, global_n1_numerical_field=False,
            numerical_oracle_available=False, physical_coordinate_velocity_export=False,
            temporal_recursion=False, full_corrected_Navier_Stokes_solution=False,
            unresolved=["axis and Cartesian regularity", "inner high radial derivatives outside the flat collar",
                        "complete lower high derivative interval contract", "original numerical evaluation/bounds",
                        "physical recovery and genuine n>=2 temporal recursion"])
        self._validate_rows()
        self.source_construction_end = len(self.graph.nodes)

    def _install(self, name, geometry, own, source):
        self.charts[name] = dict(coordinate=geometry["coordinate"], radius=geometry["radius"],
            domain=list(geometry["domain"]), fields=own["fields"], cumulative=own["cumulative"],
            source=source, lower_open=False)

    def _install_inner_core(self, rin):
        """Lift the saved core rho rows with D_y=rho D_rho, without division at rho=0."""
        g = self.graph
        inner = rin["inputs"]["inner_Volterra_bindings"]
        cut = rin["inputs"]["inner_cutoff_bindings"]
        x = inner["coordinates"]["core"]
        rho = inner["systems"]["core"]["rho"]
        radius = g.mul(inner["parameters"]["epsilon"], rho)
        raw = cut["fields"]["core"]
        native = raw["rows"]
        scale = lambda jet, coefficient: [g.mul(i, coefficient) for i in jet]
        add = lambda *jets: [g.add(*items) for items in zip(*jets)]
        rho2 = g.mul(rho, rho)
        fields = {}
        for key in ("F1", "Uz1", "mean1", "K1"):
            v, r, rr = native[key]
            fields[key] = [v, scale(r, rho), add(scale(r, rho), scale(rr, rho2))]
        mean = native["mean1"]
        fields["M1"] = [raw["radial_primitive_M1"],
            scale(add(mean[0], scale(mean[1], rho)), radius),
            scale(add(mean[0], scale(mean[1], g.mul(g.constant(3), rho)), scale(mean[2], rho2)), radius)]
        F = native["F1"]
        root = g.unary("sqrt", g.mul(g.constant(2), radius))
        fields["Utheta1"] = [scale(F[0], root),
            scale(add(scale(F[0], g.constant("1/2")), scale(F[1], rho)), root),
            scale(add(scale(F[0], g.constant("1/4")), scale(F[1], g.mul(g.constant(2), rho)),
                      scale(F[2], rho2)), root)]
        for key in ("Q1", "V1"):
            v, r, rr = raw[key]
            fields[key] = [v, scale(r, rho), add(scale(r, rho), scale(rr, rho2))]
        fields["P1"] = [raw["pressure_P1"], scale(raw["pressure_rho"], rho)]
        cumulative = {key: [value, scale(cut["integrands"]["core"][key], rho)]
                      for key, value in cut["cumulative"]["core"].items()}
        self._install("inner_core", dict(coordinate=x, radius=radius, domain=[0, g.constant(4)]),
                      dict(fields=fields, cumulative=cumulative), "n1_Rin_source_join")
        self.charts["inner_core"]["lower_open"] = True

    def _validate_rows(self):
        for name, chart in self.charts.items():
            if set(chart["fields"]) != set(STATE) | {"P1"} or set(chart["cumulative"]) != set(MOMENTS):
                raise ValueError("Incomplete own source fields: " + name)
            for group in ("fields", "cumulative"):
                for field, jets in chart[group].items():
                    zrows = 3 if field in ("Q1", "V1") else 4
                    if not jets or any(len(jet) != zrows for jet in jets):
                        raise ValueError("Wrong ordinary-Z contract: " + name + ":" + field)
                    if any(type(i) is not int or not 0 <= i < len(self.graph.nodes) for jet in jets for i in jet):
                        raise ValueError("Unbound source node")

    def contract(self, chart):
        if chart not in self.charts:
            raise ValueError("Unknown original chart: " + str(chart))
        packet = self.charts[chart]
        return {group: {field: dict(physical_y_order=len(jets)-1, ordinary_Z_order=len(jets[0])-1)
                        for field, jets in packet[group].items()} for group in ("fields", "cumulative")}

    def query(self, chart, field, *, fraction=None, tail_log_offset=None, group="fields",
              physical_y_order=0, ordinary_Z_order=0):
        if chart not in self.charts:
            raise ValueError("Unknown original chart: " + str(chart))
        packet = self.charts[chart]
        if group not in ("fields", "cumulative") or field not in packet[group]:
            raise ValueError("Choose an actual own field or cumulative moment")
        for order in (physical_y_order, ordinary_Z_order):
            if type(order) is not int or order < 0:
                raise TypeError("Derivative orders must be nonnegative integers")
        jets = packet[group][field]
        if physical_y_order >= len(jets) or ordinary_Z_order >= len(jets[physical_y_order]):
            raise ValueError("Unsupported derivative for this region; inspect contract(chart)")
        if chart == "tail":
            if fraction is not None or tail_log_offset is None:
                raise ValueError("Tail requires only tail_log_offset=log(R/Rv)")
            q = _rational(tail_log_offset)
            if q < 0:
                raise ValueError("Tail log(R/Rv) must be nonnegative")
        else:
            if fraction is None or tail_log_offset is not None:
                raise ValueError("Finite charts require only fraction of the native domain")
            q = _rational(fraction)
            if not 0 <= q <= 1 or q == 0 and packet["lower_open"]:
                raise ValueError("Uncovered coordinate (the axis is not admitted)")
        g = self.graph
        if chart == "tail":
            value = g.constant(q)
        else:
            lo, hi = packet["domain"]
            value = lo if q == 0 else hi if q == 1 else g.add(lo, g.mul(g.constant(q), g.add(hi, g.neg(lo))))
        node = jets[physical_y_order][ordinary_Z_order]
        return SourceQuery(chart, group, field, physical_y_order, ordinary_Z_order, str(q), node,
            g.at(node, packet["coordinate"], value), value,
            g.at(packet["radius"], packet["coordinate"], value))

    def manifest(self):
        return dict(source_family=self.source_family, source_archives=self.source_archives,
            input_hashes=self.input_hashes, archived_graph_nodes=self.graph.archived_nodes,
            charts={name: dict(domain=chart["domain"], lower_open=chart["lower_open"],
                    coordinate=chart["coordinate"], radius=chart["radius"], source=chart["source"],
                    contract=self.contract(name)) for name, chart in self.charts.items()},
            source_construction_suffix=self.graph.nodes[self.graph.archived_nodes:self.source_construction_end],
            coverage=self.coverage)

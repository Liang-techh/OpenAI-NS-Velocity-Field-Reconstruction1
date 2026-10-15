"""Export/check a regional physical n1 source coefficient (not numeric samples)."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from openai_ns_reconstruction.original_n1_physical_source import OriginalN1PhysicalSource

STEM = "lei_ren_part1_paper_compliant_current_original_n1_parametric_physical_source"
DIRECTORY = ROOT/"experiments"/"root_st073"


class ScalarAlgebra:
    def __init__(self,graph,held):
        import sympy as s
        self.s,self.g,self.cache = s,graph,dict(held)

    def at(self,i):
        if i in self.cache:
            return self.cache[i]
        s = self.s
        row = self.g.nodes[i]
        op = row["operation"]
        if op == "exact_rational":
            out = s.Rational(row["numerator"],row["denominator"])
        elif op == "sum":
            out = s.Add(*(self.at(j) for j in row["arguments"]))
        elif op == "product":
            out = s.Mul(*(self.at(j) for j in row["arguments"]))
        elif op == "negative":
            out = -self.at(row["argument"])
        elif op == "positive_quotient":
            out = self.at(row["numerator"])/self.at(row["denominator"])
        elif op == "analytic_unary":
            out = {"exp":s.exp,"log":s.log,"sqrt":s.sqrt,"cos":s.cos,"sin":s.sin}[row["name"]](self.at(row["argument"]))
        else:
            raise ValueError("An original source input was not bound in the independent fixture: "+op)
        self.cache[i] = out
        return out


def check(adapter):
    import sympy as s
    source,g = adapter.source,adapter.graph
    packets = []
    mass_rows = 0
    regional_functions = 0
    for chart in source.charts:
        args = dict(tail_log_offset="1") if chart == "tail" else dict(fraction="1/2")
        packet = adapter.query(chart,tau="1/2",Z="1/4",azimuth="1/3",**args)
        assert packet["velocity"]["w"] == packet["velocity"]["uz"]
        assert packet["definition"]["coefficient_order"] == 1
        assert not packet["definition"]["numerical_oracle_available"]
        assert packet["definition"]["M1_is_original_cumulative_m11"]
        assert packet["source_inputs"]["P1"] != 0
        mass_rows += len(source.charts[chart]["fields"]["M1"][0])
        packets.append(packet)
        whole = adapter.regional_function(chart)
        assert whole["definition"]["complete_regional_spatial_source"]
        assert whole["definition"]["free_time_variable"]
        assert whole["source_inputs"]["M1"] == source.charts[chart]["fields"]["M1"][0][0]
        assert whole["parameter_domain"]["native_coordinate"] == source.charts[chart]["coordinate"]
        regional_functions += 1
    # Interpret actual export nodes on arbitrary independent source functions.
    # These are algebra fixtures, never substituted original point values.
    packet = adapter.query("inner_first",fraction="5/8",tau="1/2")
    R,lam = s.symbols("R lambda",positive=True)
    Z,phi,de = s.symbols("Z phi delta",real=True)
    M = s.Function("original_mass")(R,Z)
    Ut = s.Function("original_Utheta1")(R,Z)
    P = s.Function("original_P1")(R,Z)
    My,MZ = R*s.diff(M,R),s.diff(M,Z)
    inputs = packet["source_inputs"]
    held = {packet["coordinates"]["R"]:R,packet["coordinates"]["lambda_source"]:lam,
        packet["coordinates"]["Z"]:Z,packet["coordinates"]["azimuth"]:phi,inputs["delta"]:de,
        inputs["M1"]:M,inputs["M1_y"]:My,inputs["M1_Z"]:MZ,inputs["Utheta1"]:Ut,inputs["P1"]:P}
    a = ScalarAlgebra(g,held)
    d,L = 1-Z**2,1-de*Z**2
    r = lam*s.sqrt(2*R)
    ur = lam**(2*de)/(r*L)*(2*Z*My-(1+de)*Z*M-d*MZ)
    uz = lam**(-1+de)*s.diff(M,R)
    ut = lam**(-1+de)*Ut
    wanted = dict(u=ur*s.cos(phi)-ut*s.sin(phi),v=ur*s.sin(phi)+ut*s.cos(phi),w=uz,
                  ur=ur,utheta=ut,uz=uz)
    identities = 0
    for field,value in wanted.items():
        assert s.simplify(a.at(packet["velocity"][field])-value) == 0,field
        identities += 1
    for node,value in ((packet["streamfunction"],lam**(1+de)*M),(packet["pressure"],lam**-2*P),
                       (packet["coordinates"]["r"],r),(packet["coordinates"]["z"],lam**(1-de)*Z),
                       (packet["coordinates"]["x"],r*s.cos(phi)),(packet["coordinates"]["y"],r*s.sin(phi))):
        assert s.simplify(a.at(node)-value) == 0
        identities += 1
    # Independent physical-coordinate chain rule at fixed tau and fixed r.
    # lambda_Z/lambda=Z/d; R_Z=-2ZR/d; z_Z=lambda^(1-delta)L/d.
    psi = lam**(1+de)*M
    dz = lambda f: d/(lam**(1-de)*L)*(s.diff(f,Z)+lam*Z/d*s.diff(f,lam)-2*Z*R/d*s.diff(f,R))
    dr = lambda f: r/lam**2*s.diff(f,R)
    assert s.simplify(-dz(psi)/r-ur) == 0
    assert s.simplify(dr(psi)/r-uz) == 0
    assert s.simplify(dr(r*ur)/r+dz(uz)) == 0
    # Recover the printed Q formula independently from the streamfunction.
    U,K = s.symbols("Uz1 K1",real=True)
    mean = U+K
    meanZ = s.Symbol("mean1_Z")
    ownQ = ((1-de)*Z*U-(1+de)*Z*K-d*meanZ)/L
    expectedQ = (2*Z*U-(1+de)*Z*mean-d*meanZ)/L
    assert s.cancel(ownQ-expectedQ) == 0
    # Physical normalization and cutoff-independent domain guards.
    rejected = 0
    invalid = [dict(tau="0"),dict(tau="-1"),dict(tau="3/2"),dict(tau=0.5),
               dict(Z="1"),dict(Z="-1"),dict(Z="3/2"),dict(Z=0.25),dict(azimuth=True)]
    for kwargs in invalid:
        before = len(g.nodes)
        try:
            adapter.query("inner_core",fraction="1/2",**kwargs)
        except (TypeError,ValueError):
            pass
        else:
            raise AssertionError("Invalid physical source parameters accepted")
        assert len(g.nodes) == before
        rejected += 1
    try:
        adapter.query("inner_core",fraction="0")
    except ValueError:
        rejected += 1
    else:
        raise AssertionError("Unproved physical axis admitted")
    return dict(passed=True,regional_physical_packets=len(packets),actual_mass_cumulative_source_rows=mass_rows,
        complete_regional_spacetime_source_exports=regional_functions,
        direct_mass_source_identities=sum(e["method"] == "identical_source_node" for e in adapter.mass_identity_evidence),
        capture_free_tail_mass_source_identities=sum(e["method"] != "identical_source_node" for e in adapter.mass_identity_evidence),
        actual_export_algebra_identities=identities,independent_streamfunction_chain_identities=2,
        independent_axisymmetric_divergence_identity=True,own_Q_radial_scaling_identity=True,
        rejected_invalid_or_axis_requests=rejected,mathematical_source_fixtures_only=True,
        original_numeric_field_values_evaluated=False,global_interface_germs_proved=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check",action="store_true")
    p.add_argument("--write-report",action="store_true")
    p.add_argument("--restore-source",action="store_true")
    p.add_argument("--regional-function",action="store_true")
    p.add_argument("--chart",default="inner_first")
    p.add_argument("--fraction",default="5/8")
    p.add_argument("--tail-log-offset")
    p.add_argument("--tau")
    p.add_argument("--Z")
    p.add_argument("--azimuth")
    args = p.parse_args()
    if args.write_report and not args.check:
        p.error("--write-report requires --check")
    began = time.monotonic()
    if args.restore_source:
        if args.check or args.write_report:
            p.error("Physical export recovery is separate from construction/checking")
        adapter,packet = OriginalN1PhysicalSource.restore_export(DIRECTORY/(STEM+".json.gz"))
        roots = list(packet["velocity"].values())+[packet["pressure"],packet["streamfunction"]]
        assert all(type(i) is int and 0 <= i < len(adapter.graph.nodes) for i in roots)
        print(json.dumps(dict(source_only_restored=True,full_owner_restored=False,
            actual_physical_source_roots=len(roots),query_suffix_nodes=len(adapter.graph.nodes)-adapter.prefix_nodes,
            graph_prefix_nodes=adapter.prefix_nodes,definition=packet["definition"],
            execution_seconds=round(time.monotonic()-began,3)),indent=2))
        return
    adapter = OriginalN1PhysicalSource()
    if args.regional_function:
        packet = adapter.regional_function(args.chart,tau=args.tau,Z=args.Z,azimuth=args.azimuth)
    else:
        packet = adapter.query(args.chart,fraction=None if args.chart == "tail" else args.fraction,
            tail_log_offset=args.tail_log_offset,tau=args.tau or "1/2",Z=args.Z,azimuth=args.azimuth)
    payload = adapter.snapshot(packet)
    report = dict(source_only=True,physical_source=packet,
                  graph_prefix_nodes=adapter.prefix_nodes,query_suffix_nodes=len(payload["query_graph_suffix"]))
    if args.check:
        report["checks"] = check(adapter)
    if args.write_report:
        files = [Path(__file__),ROOT/"src/openai_ns_reconstruction/original_n1_physical_source.py"]
        pins = {str(path.relative_to(ROOT)).replace("\\","/"):hashlib.sha256(path.read_bytes()).hexdigest() for path in files}
        payload["implementation_hashes"] = pins
        path = DIRECTORY/(STEM+".json.gz")
        path.write_bytes(gzip.compress((json.dumps(payload,separators=(",",":"))+"\n").encode(),mtime=0))
        with gzip.open(path,"rt",encoding="utf-8") as f:
            assert json.load(f) == payload
        report.update(implementation_hashes=pins,archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      source_graph_chain=payload["source_graph_chain"],lossless_source_export_roundtrip=True)
    report["execution_seconds"] = round(time.monotonic()-began,3)
    if args.write_report:
        (DIRECTORY/(STEM+"_check.json")).write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()

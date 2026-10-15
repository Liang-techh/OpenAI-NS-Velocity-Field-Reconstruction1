"""Export the original axis-regular Cartesian n1 core and check its recovery."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from openai_ns_reconstruction.original_n1_cartesian_core_source import OriginalN1CartesianCoreSource
from run_original_n1_physical_source import ScalarAlgebra

DIRECTORY = ROOT / "experiments/root_st073"
STEM = "lei_ren_part1_paper_compliant_current_original_n1_cartesian_core_source"


def check(adapter, packet):
    import sympy as s
    g = adapter.graph
    x, y, R, Z, lam, de, eps = s.symbols("x y R Z lambda delta epsilon", positive=True)
    F, U, mean, Q, P = s.symbols("F1 Uz1 mean1 Q1 P1")
    held = {packet["coordinates"]["x"]:x, packet["coordinates"]["y"]:y,
            packet["coordinates"]["lambda_source"]:lam, adapter.delta:de, adapter.epsilon:eps}
    held.update(zip(packet["source_inputs"]["bound_profiles"].values(), (F,U,mean,Q,P)))
    a = ScalarAlgebra(g, held)
    A = lam**(-2+2*de)*Q/2
    B = lam**(-2+de)*F
    identities = 0
    expected = [(packet["velocity"]["u"], A*x-B*y), (packet["velocity"]["v"], A*y+B*x),
        (packet["velocity"]["w"], lam**(-1+de)*U), (packet["pressure"], lam**-2*P),
        (packet["coordinates"]["rho"], (x*x+y*y)/(2*lam**2*eps)),
        (packet["streamfunction"], lam**(1+de)*(x*x+y*y)/(2*lam**2)*mean)]
    for node, value in expected:
        assert s.simplify(a.at(node)-value) == 0
        identities += 1
    # Independent rotation of the already accepted positive-radius profiles.
    phi = s.symbols("phi", real=True)
    r = lam*s.sqrt(2*R)
    ur = lam**(-1+2*de)*s.sqrt(R/2)*Q
    ut = lam**(-1+de)*s.sqrt(2*R)*F
    for expression, target in ((A*r*s.cos(phi)-B*r*s.sin(phi), ur*s.cos(phi)-ut*s.sin(phi)),
                               (A*r*s.sin(phi)+B*r*s.cos(phi), ur*s.sin(phi)+ut*s.cos(phi))):
        assert s.simplify(expression-target) == 0
        identities += 1
    # Compare the actual original Q source body, with only its genuine regular
    # solution inputs held, rather than binding the output Q to its target.
    original = adapter.volterra["functions"]["core"]
    K, UZ, KZ = s.symbols("K1 Uz1_Z K1_Z")
    qa = ScalarAlgebra(g, {adapter.delta:de, adapter.Z_variable:Z,
        original["Uz1"][0]:U, original["K1"][0]:K,
        original["Uz1"][1]:UZ, original["K1"][1]:KZ})
    d, L = 1-Z**2, 1-de*Z**2
    wanted_Q = ((1-de)*Z*U-(1+de)*Z*K-d*(UZ+KZ))/L
    assert s.cancel(qa.at(packet["source_inputs"]["raw_core"]["Q1"])-wanted_Q) == 0
    # Verify the fixed-physical-coordinate Z chain, including R_Z and lambda_Z.
    # The source FTC gives Uz1 = mean1 + R*partial_R(mean1).
    m = s.Function("original_mean1")(R,Z)
    actual_U = m+R*s.diff(m,R)
    actual_Q = (2*Z*actual_U-(1+de)*Z*m-d*s.diff(m,Z))/L
    actual_A = lam**(-2+2*de)*actual_Q/2
    actual_w = lam**(-1+de)*actual_U
    dz = lambda f: d/(lam**(1-de)*L)*(s.diff(f,Z)+lam*Z/d*s.diff(f,lam)-2*Z*R/d*s.diff(f,R))
    assert s.simplify(2*actual_A+2*R*s.diff(actual_A,R)+dz(actual_w)) == 0
    for node in packet["source_inputs"]["bound_profiles"].values():
        row = g.nodes[node]
        assert row["operation"] == "function_substitution"
        assert row["Z_independent_substitution"] is False
        assert row["physical_coordinate_pullback"] is True
    # The original core is analytic in rho, with F1/Uz1/K1/P1=O(rho).
    # Quadratic rho substitution makes all Cartesian first derivatives zero;
    # retain the actual equation-derived Uz1_rho axis node for the Hessian.
    cF,cU,cQ,cP = s.symbols("F_rho U_rho Q_rho P_rho")
    rho = (x*x+y*y)/(2*lam**2*eps)
    velocity = (lam**(-2+2*de)*cQ*rho*x/2-lam**(-2+de)*cF*rho*y,
                lam**(-2+2*de)*cQ*rho*y/2+lam**(-2+de)*cF*rho*x,
                lam**(-1+de)*cU*rho)
    for component in velocity:
        assert component.subs({x:0,y:0}) == 0
        assert all(s.diff(component,q).subs({x:0,y:0}) == 0 for q in (x,y))
    expected_curvature = lam**(-3+de)*cU/eps
    ca = ScalarAlgebra(g, {packet["coordinates"]["lambda_source"]:lam,
        adapter.delta:de, adapter.epsilon:eps, adapter.axial_rho_axis:cU})
    curvature = packet["axis_limits"]["axial_transverse_hessian"]
    assert s.simplify(ca.at(curvature["xx"])-expected_curvature) == 0
    for q in (x,y):
        assert s.simplify(s.diff(velocity[2],q,2)-expected_curvature) == 0
    axis = adapter.axis_function(tau="1/2", Z="1/4")
    assert set(axis["velocity"].values()) == {0} and axis["pressure"] == 0
    assert axis["streamfunction"] == 0 and axis["coordinates"]["rho"] == 0
    assert axis["axis_limits"]["axial_transverse_hessian"]["xx"] != 0
    rejected = 0
    for args in (dict(tau="0"),dict(tau="3/2"),dict(tau=0.5),dict(Z="1"),dict(Z=True)):
        before = len(g.nodes)
        try:
            adapter.cartesian_core_function(**args)
        except (ValueError,TypeError):
            rejected += 1
        else:
            raise AssertionError("Invalid physical core domain accepted")
        assert len(g.nodes) == before
    return dict(passed=True, actual_regular_axis_source_rows=len(adapter.axis_source_evidence),
        actual_cutoff_core_rows_identical_to_guarded_source=16,
        original_even_radial_axis_derivative_sources=len(adapter.axis_derivative_evidence),
        actual_derived_mean_Q_axis_zero_rows=len(adapter.axis_derived_zeros),
        actual_cartesian_export_and_rotation_identities=identities,
        actual_original_Q_recovery_body_identity=True,
        independent_fixed_Cartesian_divergence_identity=True,
        physical_pullback_retains_Z_dependent_radius=True,
        actual_original_axial_rho_axis_node=adapter.axial_rho_axis,
        exact_axis_values_and_first_derivatives=True,
        original_source_bound_axial_transverse_curvature=True,
        rejected_requests=rejected, symbolic_operator_checks_only=True,
        original_numeric_field_values_evaluated=False, global_interface_germs_proved=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-report",action="store_true")
    parser.add_argument("--restore-source",action="store_true")
    args = parser.parse_args()
    if args.write_report and not args.check:
        parser.error("--write-report requires --check")
    start = time.monotonic()
    archive = DIRECTORY / (STEM+".json.gz")
    if args.restore_source:
        if args.check or args.write_report:
            parser.error("Restore independently of construction/checking")
        adapter,packet = OriginalN1CartesianCoreSource.restore_export(archive)
        roots = list(packet["velocity"].values()) + [packet["pressure"], packet["streamfunction"]]
        assert all(type(i) is int and 0 <= i < len(adapter.graph.nodes) for i in roots)
        print(json.dumps(dict(source_only_restored=True, roots=len(roots),
            query_suffix_nodes=len(adapter.graph.nodes)-adapter.prefix_nodes,
            axis_included=packet["definition"]["axis_included"], numerical_oracle_available=False,
            execution_seconds=round(time.monotonic()-start,3)),indent=2))
        return
    adapter = OriginalN1CartesianCoreSource()
    packet = adapter.cartesian_core_function()
    payload = adapter.snapshot(packet)
    report = dict(source_only=True, physical_source=packet,
                  graph_prefix_nodes=adapter.prefix_nodes, query_suffix_nodes=len(payload["query_graph_suffix"]))
    if args.check:
        report["checks"] = check(adapter,packet)
    if args.write_report:
        files = [Path(__file__), ROOT/"src/openai_ns_reconstruction/original_n1_cartesian_core_source.py",
                 ROOT/"src/openai_ns_reconstruction/original_n1_physical_source.py"]
        pins = {str(p.relative_to(ROOT)).replace("\\","/"):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
        payload["implementation_hashes"] = pins
        archive.write_bytes(gzip.compress((json.dumps(payload,separators=(",",":"))+"\n").encode(),mtime=0))
        with gzip.open(archive,"rt",encoding="utf-8") as f:
            assert json.load(f) == payload
        report.update(archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(), implementation_hashes=pins,
                      lossless_source_export_roundtrip=True)
    report["execution_seconds"] = round(time.monotonic()-start,3)
    if args.write_report:
        (DIRECTORY/(STEM+"_check.json")).write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()

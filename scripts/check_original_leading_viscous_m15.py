"""Independent differentiation checks for the new signed viscous operators."""
from pathlib import Path
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_leading_viscous_kernel import (
    decode_pair, encode_pair, endpoint_coefficient, mean_source_coefficient,
    mean_source_value,
)


def run_fixture():
    c = mp.mp.clone(); c.dps = 100
    iv = MPIntervalContext(); iv.dps = 120
    delta = c.mpf(".02")
    controls = []
    differences = []
    def passed(name, condition):
        if not condition:
            raise AssertionError(name)
        controls.append({"name": name, "passed": True})
    def b(z):
        return c.mpf("1.2") + c.mpf(".9")*z - c.mpf(".7")*z*z + c.mpf(".4")*z*z*z
    def endpoint(z):
        return 2*z*b(z)/(1-delta*z*z)
    def mean(z):
        return ((1-delta)*z*b(z)+(1-z*z)*c.diff(b,z))/(1-delta*z*z)
    rows = []
    for text in ("-.3", "0", ".23"):
        z = c.mpf(text)
        inputs = [c.diff(b,z,n) for n in range(3)]
        direct = endpoint_coefficient(z, delta, inputs)
        recovered = mean_source_coefficient(z, delta, inputs)
        for label, values, oracle in (("endpoint", direct, endpoint), ("mean", recovered, mean)):
            for order in (0,1):
                expected = c.diff(oracle,z,order)
                difference = abs(values[order]-expected)
                differences.append(difference)
                passed(label+"_"+text+"_ordinary_Z"+str(order), difference < c.mpf("1e-90"))
        passed("mean_C0_two_inputs_"+text,
               abs(mean_source_value(z,delta,inputs[:2])-recovered[0]) < c.mpf("1e-90"))
        interval_inputs = [iv.mpf([c.nstr(value-c.mpf("1e-95"),c.dps),
                                   c.nstr(value+c.mpf("1e-95"),c.dps)]) for value in inputs]
        enclosed = endpoint_coefficient(iv.mpf(text),iv.mpf(".02"),interval_inputs)
        passed("endpoint_directed_"+text, all(
            iv.mpf(c.nstr(value,c.dps)) in box for value,box in zip(direct,enclosed)))
        rows.append({"Z":text,"endpoint_C0_C1":[c.nstr(v,30) for v in direct],
                     "mean_C0_C1":[c.nstr(v,30) for v in recovered]})
    z = c.mpf(".23"); inputs = [c.diff(b,z,n) for n in range(3)]
    passed("upper_lower_sign", endpoint_coefficient(z,delta,inputs,sign=-1)
           == [-v for v in endpoint_coefficient(z,delta,inputs)])
    passed("missing_Bz_rejected", _raises(lambda: endpoint_coefficient(z,delta,inputs[:1])))
    passed("missing_Bzz_rejected", _raises(lambda: mean_source_coefficient(z,delta,inputs[:2])))
    passed("Taylor_second_coefficient_not_ordinary", abs(
        mean_source_coefficient(z,delta,[inputs[0],inputs[1],inputs[2]/2])[1]
        - c.diff(mean,z)) > c.mpf(".01"))
    passed("derivative_of_rational_factor_required", abs(
        2*z*inputs[1]/(1-delta*z*z)-c.diff(endpoint,z)) > c.mpf(".01"))
    box = iv.mpf(["1.1","1.3"])
    summary = encode_pair(c,iv,c.mpf("1.2"),box)
    point, restored = decode_pair(c,iv,summary)
    passed("exact_enclosure_roundtrip", point==c.mpf("1.2") and restored._mpi_==box._mpi_)
    passed("point_not_promoted", summary["approximation_is_exact_source_value"] is False)
    passed("cap_input_rejected", _raises(lambda: decode_pair(c,iv,dict(summary,positive_error_bound_atom_only=True))))
    passed("log_atom_rejected", _raises(lambda: decode_pair(c,iv,dict(summary,source_log_atom_only=True))))
    # An independent analytic mean history verifies the signed FTC partition.
    # U=a+b*y^2 and m=a+b*(y^2-2y+2)+history*exp(-y) solve m_y=U-m.
    lower,upper = c.mpf(".2"),c.mpf(".8")
    def amplitude(q):
        return c.mpf(".7")-q/5+3*q*q*q/10
    def history(q):
        return 2-q+q*q
    def h(y,q):
        source = 2*y*amplitude(q)
        n = amplitude(q)*(2*y-2)-history(q)*c.exp(-y)
        nz = c.diff(amplitude,q)*(2*y-2)-c.diff(history,q)*c.exp(-y)
        return c.exp(y)*(2*q*source-(1-delta)*q*n-(1-q*q)*nz)/(1-delta*q*q)
    def e(y,q):
        return c.exp(y)*2*q*(2*y*amplitude(q))/(1-delta*q*q)
    def source_primitive(q):
        factor = 2*(c.exp(upper)*(upper-1)-c.exp(lower)*(lower-1))
        return factor*mean_source_value(q,delta,[amplitude(q),c.diff(amplitude,q)])
    def direct_interval(q):
        return h(lower,q)-h(upper,q)
    def partition(q):
        return e(lower,q)-e(upper,q)+source_primitive(q)
    for text in ("-.3","0",".23"):
        z = c.mpf(text)
        for order in (0,1):
            difference = abs(c.diff(direct_interval,z,order)-c.diff(partition,z,order))
            differences.append(difference)
            passed("signed_mean_FTC_"+text+"_ordinary_Z"+str(order),difference<c.mpf("1e-90"))
    history_H = mean_source_value(c.mpf(".23"),delta,[history(c.mpf(".23")),c.diff(history,c.mpf(".23"))])
    passed("zero_Uy_does_not_zero_RQy_or_mean_history",abs(history_H)>c.mpf(".01"))
    return {"passed":True,"controls":controls,"control_count":len(controls),"samples":rows,
            "maximum_supporting_difference":c.nstr(max(differences),20),
            "independent_oracle":"mpmath differentiation of nonconstant polynomial source and rational pressure",
            "ordinary_Z1":True,"source_values_of_original_field_evaluated":False,
            "ancestor_constructor_producer_checker_fixture_run":False,"NS_residual":False}


def _raises(function):
    try:
        function()
    except ValueError:
        return True
    return False


if __name__ == "__main__":
    import json
    print(json.dumps(run_fixture(),indent=2))

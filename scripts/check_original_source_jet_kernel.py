"""Independent analytic checks for the new ordinary source-jet arithmetic."""
from pathlib import Path
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

from openai_ns_reconstruction.original_source_jet_kernel import (
    ordinary_exp, ordinary_product, ordinary_quotient, transported_lower_dy,
)


def _raises(function):
    try:
        function()
    except ValueError:
        return True
    return False


def run_fixture():
    c = mp.mp.clone(); c.dps = 100
    iv = MPIntervalContext(); iv.dps = 120
    controls = []; differences = []
    def passed(name,condition):
        if not condition:
            raise AssertionError(name)
        controls.append({"name":name,"passed":True})
    def a(z):
        return c.exp(z/3)+z*z
    def b(z):
        return 2+z/7+z**3/5
    def p(z):
        return c.sin(z)+c.mpf(".7")
    def rho(z):
        return 4+z*z/10
    samples = []
    for text in ("-.3","0",".23"):
        z = c.mpf(text)
        jets = [[c.diff(f,z,k) for k in range(3)] for f in (a,b,p,rho)]
        operations = (
            ("product",ordinary_product(jets[0],jets[1]),lambda q:a(q)*b(q)),
            ("quotient",ordinary_quotient(jets[0],jets[1]),lambda q:a(q)/b(q)),
            ("exponential",ordinary_exp(c,jets[2]),lambda q:c.exp(p(q))),
            ("transported_lower_dy",transported_lower_dy(*jets),lambda q:a(q)/b(q)*p(q)*rho(q)),
        )
        for label,result,oracle in operations:
            for k,value in enumerate(result):
                difference = abs(value-c.diff(oracle,z,k))
                differences.append(difference)
                passed(label+"_"+text+"_ordinary_Z"+str(k),difference<c.mpf("1e-90"))
            passed(label+"_Z1_order_supported_"+text,
                   (ordinary_product(jets[0][:2],jets[1][:2]) if label=="product" else
                    ordinary_quotient(jets[0][:2],jets[1][:2]) if label=="quotient" else
                    ordinary_exp(c,jets[2][:2]) if label=="exponential" else
                    transported_lower_dy(*(row[:2] for row in jets)))==result[:2])
        boxes = [[iv.mpf([c.nstr(v-c.mpf("1e-95"),100),
                             c.nstr(v+c.mpf("1e-95"),100)]) for v in row] for row in jets]
        truth = transported_lower_dy(*boxes)
        point = transported_lower_dy(*jets)
        passed("directed_transport_"+text,all(iv.mpf(c.nstr(v,100)) in box for v,box in zip(point,truth)))
        samples.append({"Z":text,"lower_dy_ordinary_Z0_Z1_Z2":[c.nstr(v,30) for v in point]})
    row = [c.mpf(1),c.mpf(2),c.mpf(3)]
    passed("missing_derivative_rejected",_raises(lambda:ordinary_product(row,[c.mpf(1)])))
    passed("None_row_rejected",_raises(lambda:ordinary_product(row,[c.mpf(1),None,c.mpf(2)])))
    passed("unequal_orders_rejected",_raises(lambda:ordinary_quotient(row,row[:2])))
    passed("zero_denominator_rejected",_raises(lambda:ordinary_quotient(row,[c.mpf(0),c.mpf(1),c.mpf(2)])))
    passed("negative_denominator_rejected",_raises(lambda:ordinary_quotient(row,[-c.mpf(1),c.mpf(1),c.mpf(2)])))
    passed("denominator_straddling_zero_rejected",_raises(lambda:ordinary_quotient(
        [iv.mpf(v) for v in row],[iv.mpf([-1,1]),iv.mpf(1),iv.mpf(2)])))
    z = c.mpf(".23"); jets = [[c.diff(f,z,k) for k in range(3)] for f in (a,b,p,rho)]
    actual = transported_lower_dy(*jets)
    passed("no_extra_macro_J",abs(actual[0]*17-c.diff(lambda q:a(q)/b(q)*p(q)*rho(q),z,0))>1)
    halved = [row[:] for row in jets]; halved[0][2] /= 2
    passed("Taylor_Z2_not_ordinary_Z2",abs(transported_lower_dy(*halved)[2]-actual[2])>c.mpf(".01"))
    from openai_ns_reconstruction.original_lower_dy_source import _symmetric, _upper_abs
    fine = MPIntervalContext(); fine.dps = 180
    narrow = fine.mpf(1)+fine.mpf(2)**-500
    raw = _upper_abs(c,narrow)
    passed("cap_upper_not_rounded_in_point_context",mpf_cmp(raw._mpf_,narrow._mpi_[1])>=0)
    symmetric = _symmetric(fine,c,narrow)
    passed("symmetric_cap_keeps_directed_high_precision",mpf_cmp(symmetric._mpi_[1],narrow._mpi_[1])>=0)
    return {"passed":True,"control_count":len(controls),"controls":controls,"samples":samples,
            "maximum_supporting_difference":c.nstr(max(differences),25),
            "independent_oracle":"analytic mpmath differentiation of nonconstant exp, polynomial and sin functions",
            "ancestor_constructor_producer_checker_fixture_run":False,
            "original_source_rows_evaluated":False,"temporal_n_recursion":False,"NS_residual":False}


if __name__=="__main__":
    import json
    print(json.dumps(run_fixture(),indent=2))

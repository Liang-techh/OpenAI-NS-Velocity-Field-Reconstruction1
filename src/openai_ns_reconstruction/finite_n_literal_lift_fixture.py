"""New scoped literal-lift verification; no density or ancestor math."""
from fractions import Fraction
from openai_ns_reconstruction.original_finite_N_density_kernel import ScalarOperations, coerce_ordinary_jet

def run_fixture():
    calls=[]
    def deny(*args):
        raise AssertionError("arithmetic must not be entered")
    def constant(value):
        calls.append(value)
        return ("constant",value)
    ops=ScalarOperations(deny,deny,deny,deny,deny,constant)
    value=Fraction(3,7)
    numeric=coerce_ordinary_jet(value,operations=ops,name="literal")
    checks={"explicit_numeric_literal_and_zeros_use_constant": calls==[value,0,0] and numeric.value==("constant",value) and numeric.dz==("constant",0) and numeric.dzz==("constant",0) and numeric.operations is ops}
    calls.clear()
    token=object()
    opaque=coerce_ordinary_jet(token,operations=ops,name="opaque_source")
    checks["explicit_opaque_source_preserved_zero_rows_converted"]=opaque.value is token and calls==[0,0] and opaque.operations is ops
    return {"schema":"finite_n_literal_lift_v2","checks":checks,"passed":all(checks.values()),"source_only":True,"old_math_calls":0,"runtime_model_metadata":None}

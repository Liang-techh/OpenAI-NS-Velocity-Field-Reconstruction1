"""Current Rsh functional mixed4 join; accepted source artifacts stay unchanged.

The constant-power continuation below computes the boundary jet only.
It never replaces the finite long-reshape neighborhood or selects a point
from a production source enclosure.
"""
import ast
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import (
    CompliantActualReferenceRestoreMixedC4, current_E_source_bindings,
    accepted, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import canonical_source
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

GROUPS=("physical_velocity_pressure_y_Z_mixed4","physical_five_primitive_y_Z_mixed4")
SCOPES=("current_actual_moment_patch_installed","full_implicit_leading_inputs_recomputed",
        "actual_point_moment_history_recovered","global_completed_tensor_admissibility",
        "full_inner_interfaces_certified","full_cartesian_vector_derivatives_certified",
        "admissible_stress_lift_constructed","temporal_recursion")


def return_call_binding(module, method, expected):
    tree=ast.parse((HERE/(PREFIX+module+".py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
    want=ast.dump(ast.parse(expected,mode="eval").body)
    if sum(isinstance(n,ast.Return) and n.value is not None and ast.dump(n.value)==want
           for n in ast.walk(fn))!=1:
        raise ValueError("Original return source changed: "+module+"."+method)
    return True



def expression_call_binding(module, method, expected):
    tree=ast.parse((HERE/(PREFIX+module+".py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
    want=ast.dump(ast.parse(expected,mode="eval").body)
    if sum(isinstance(n,ast.Call) and ast.dump(n)==want for n in ast.walk(fn))!=1:
        raise ValueError("Original source equation changed: "+module+"."+method)
    return True


def keyword_binding(module, method, keyword, expression):
    tree=ast.parse((HERE/(PREFIX+module+".py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
    wanted=ast.dump(ast.parse(expression,mode="eval").body)
    values=[kw.value for n in ast.walk(fn) if isinstance(n,ast.Call)
            for kw in n.keywords if kw.arg==keyword]
    if len(values)!=1 or ast.dump(values[0])!=wanted:
        raise ValueError("Original source keyword changed: "+module+"."+method+"."+keyword)
    return True


def class_assignment(module, class_name, method, target, expression):
    tree=ast.parse((HERE/(PREFIX+module+".py")).read_text(encoding="utf8"))
    cls=next(n for n in ast.walk(tree) if isinstance(n,ast.ClassDef) and n.name==class_name)
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==method)
    wanted=ast.dump(ast.parse(expression,mode="eval").body)
    values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(v)==target for v in n.targets)]
    if len(values)!=1 or ast.dump(values[0])!=wanted:
        raise ValueError("Original class source assignment changed: "+module+"."+target)
    return True


def source_bindings():
    left=assignment_source_bindings("long_reshape_profiles","evaluate",{
        "y":"self.T*phase",
        "logu":"B*(1-sig)-logq+(y/10-self.core.logC-self.core.logP)",
        "moment_shapes":"dict(theta=inherited['theta']+kernels['theta'],theta_z=inherited['theta_z']+inp['v']*kernels['theta'],pressure=inherited['pressure']+kernels['pressure'],swirl=inherited['swirl']+kernels['swirl'],mean=mean,axial=axial)",
    })
    log_y=assignment_source_bindings("long_reshape_mixed_C4","evaluate",{
        "log_y":"[B*(-cutoff[k]*math.factorial(k)/T**k)+(c.mpf('.1') if k==1 else 0) for k in range(1,5)]",
        "packet":"reshape_mixed(c,Z,self.reshape.core.delta,logu,log_y,inp['v'],shapes,inp['p0'],self.invP2,self.proofs)",
    })
    centered=assignment_source_bindings("actual_reference_restore_mixed_C4","inputs",{
        "parent":"self.reshape.evaluate(Z,phase=1)",
        "moments":"{n:jet(row) for n,row in parent['actual_normalized_moment_shape_axial5_coefficients'].items()}",
        "centered":"dict(mean_error=moments['mean']-z*4,angular_error=moments['theta']-c.mpf(5)/8,mixed_error=moments['theta_z']-(z*moments['theta'])*4,axial_square=moments['axial']-(z*moments['mean'])*8+square(z)*16,swirl_error=moments['swirl']-c.mpf(5)/6,pressure_error=moments['pressure']-5)",
    })
    reference=assignment_source_bindings("reference_restore_profiles","reference",{
        "gap":"(self.loggap-8)*phase",
        "centered":"self.reference_centered(Z,gap)",
        "logu":"-logq+(self.reshape.T/10-self.core.logC-self.core.logP+gap/10)",
    })
    shapes=assignment_source_bindings("reference_restore_profiles","packet",{
        "shapes":"dict(theta=H,theta_z=(z*H)*4+centered['mixed_error'],mean=z*4+centered['mean_error'],axial=square(z)*16+(z*centered['mean_error'])*8+centered['axial_square'],swirl=centered['swirl_error']+c.mpf(5)/6,pressure=centered['pressure_error']+5)",
    })
    reference_branch=return_call_binding("reference_restore_mixed_C4","reference_branch",
        "self.packet(Z,self.reference.reference(Z,phase),[c.mpf(1)]+[c.mpf(0)]*4,'Rsh_to_Rz')")
    mixed=assignment_source_bindings("reference_restore_mixed_C4","packet",{
        "p0":"jet('pressure_axis_axial5_coefficients')",
        "result":"source_mixed(c,Z,self.reference.core.delta,E,alpha,centered,logu,p0,self.invP2)",
    })
    left_physical=assignment_source_bindings("long_reshape_mixed_C4","reshape_mixed",{
        "mean":"[shapes['mean']]+[(V-shapes['mean'])*((-1)**(k-1)) for k in range(1,5)]",
        "velocity":"[V]+[zero]*4",
        "Q":"[(2*z*velocity[k]-(z*mean[k])*(1-delta)-(1-square(z))*derivative(mean[k]))/(1-square(z)*delta) for k in range(5)]",
        "physical":"dict(Utheta_over_current_Utheta=[amp*row for row in urows],Uz=velocity,Ur_over_current_sqrt_R_over_2=[binomial_rate(Q,c.mpf('.5'),k) for k in range(5)],P_over_Pstar2=[p0+pressure[0]]+pressure[1:])",
        "primitives":"dict(Mtheta_over_current_sqrt2_R_1p5_Utheta=[amp*shapes['theta']]+[amp*theta_rhs[k-1] for k in range(1,5)],Mtheta_z_over_current_sqrt2_R_1p5_Utheta=[amp*shapes['theta_z']]+[(amp*theta_rhs[k-1])*V for k in range(1,5)],Mz_over_current_R=[shapes['mean']]+[V]*4,Mztheta_over_current_R_Pstar2=[shapes['axial']*invP2-scaled(shapes['swirl'])/2]+[square(V)*invP2-scaled(swirl_rhs[k-1])/2 for k in range(1,5)],Mp_over_Pstar2=pressure)",
        "grid":"lambda rows:{'y'+str(k)+'_Z'+str(n):row[n]*math.factorial(n) for k,row in enumerate(rows) for n in range(5-k)}",
    })
    right_physical=assignment_source_bindings("reference_restore_mixed_C4","source_mixed",{
        "mismatch":"[E*value for value in alpha]",
        "V":"[z*4+mismatch[0]]+mismatch[1:]",
        "physical":"dict(Utheta_over_current_Utheta=[amp*c.mpf('.1')**k for k in range(5)],Uz=V,Ur_over_current_sqrt_R_over_2=[binomial_rate(Q,c.mpf('.5'),k) for k in range(5)],P_over_Pstar2=[p0+ratio*p[0]/2]+[ratio*binomial_rate(p,c.mpf('.2'),k)/2 for k in range(1,5)])",
        "primitives":"dict(Mtheta_over_current_sqrt2_R_1p5_Utheta=[amp*binomial_rate(H,c.mpf('1.6'),k) for k in range(5)],Mtheta_z_over_current_sqrt2_R_1p5_Utheta=[amp*binomial_rate(K,c.mpf('1.6'),k) for k in range(5)],Mz_over_current_R=[binomial_rate(m,c.mpf(1),k) for k in range(5)],Mztheta_over_current_R_Pstar2=[binomial_rate(A,c.mpf(1),k)*invP2-ratio*binomial_rate(b,c.mpf('1.2'),k)/2 for k in range(5)],Mp_over_Pstar2=[ratio*binomial_rate(p,c.mpf('.2'),k)/2 for k in range(5)])",
        "grid":"lambda values:{'y'+str(k)+'_Z'+str(n):row[n]*math.factorial(n) for k,row in enumerate(values) for n in range(5-k)}",
    })
    rates={name:expression_call_binding("reference_restore_mixed_C4","source_mixed",expression)
        for name,expression in {
            "mean":"rows['mean_error'].append(mismatch[k]-rows['mean_error'][k])",
            "theta":"rows['angular_error'].append(rows['angular_error'][k]*c.mpf('-1.6'))",
            "theta_z":"rows['mixed_error'].append(mismatch[k]-rows['mixed_error'][k]*c.mpf('1.6'))",
            "axial":"rows['axial_square'].append(quadratic-rows['axial_square'][k])",
            "swirl":"rows['swirl_error'].append(rows['swirl_error'][k]*c.mpf('-1.2'))",
            "pressure":"rows['pressure_error'].append(rows['pressure_error'][k]*c.mpf('-.2'))",
        }.items()}
    clamp=assignment_source_bindings("actual_reference_restore_mixed_C4","inputs",{
        "centered['axial_square']":"IntervalTaylor(c,[c.mpf([max(mp.mpf(0),lower),upper])]+list(row.coefficients[1:]))",
    })
    return dict(original_long_log_source=left,original_long_log_y_units=log_y,
        current_six_centered_Rsh_coordinates=centered,original_reference_log_source=reference,
        original_reference_shape_recovery=shapes,original_constant_alpha_branch=reference_branch,
        original_reference_pressure_and_operator=mixed,
        original_left_physical_rows_and_normalization=left_physical,
        original_right_physical_rows_and_normalization=right_physical,
        original_six_centered_source_equations=rates,
        original_same_parent_axis_pressure_keyword=keyword_binding("actual_reference_restore_mixed_C4","inputs","original_axis_pressure","jet(parent['pressure_axis_axial5_coefficients'])"),
        original_left_axis_pressure_keyword=keyword_binding("long_reshape_profiles","inputs","p0","jet(inlet['pressure_axis_axial5_coefficients'])"),
        left_invP2_source=class_assignment("actual_long_reshape_mixed_C4","CompliantActualLongReshapeMixedC4","__init__","self.invP2","c.exp(-2*self.reshape.core.logP)"),
        right_invP2_source=class_assignment("actual_reference_restore_mixed_C4","CompliantActualReferenceRestoreMixedC4","__init__","self.invP2","c.exp(-2*self.reference.core.logP)"),
        actual_reference_Rsh_phase_zero_invocation=keyword_binding("reference_restore_mixed_C4","report","actual_Rsh_exit","self.reference_branch([-1,1],0)"),
        axial_square_enclosure_refinement_AST=clamp,
        centered_function_identity_established_before_enclosure_clamp=True,
        clamp_is_not_a_source_function_assignment=True)


def boundary_identities():
    """Arbitrary smooth axial data: boundary extension, not left finite field."""
    y,z=s.symbols("y z",real=True)
    T,R0,Pstar=s.symbols("T Rsh Pstar",positive=True)
    delta,C,P=s.symbols("delta logC logP",real=True)
    V,U,P0,B=[s.Function(n)(z) for n in ("V","U","P0","B")]
    H,K,m,A,b,p=[s.Function(n)(z) for n in ("H0","K0","m0","A0","b0","p0shape")]
    E=V-4*z
    em=m-4*z;eh=H-s.Rational(5,8);ek=K-4*z*H
    aa=A-8*z*m+16*z*z;eb=b-s.Rational(5,6);ep=p-5
    exp=s.exp
    left=[
        s.Rational(5,8)+(H-s.Rational(5,8))*exp(-8*y/5),
        s.Rational(5,8)*V+(K-s.Rational(5,8)*V)*exp(-8*y/5),
        V+(m-V)*exp(-y),V**2+(A-V**2)*exp(-y),
        s.Rational(5,6)+(b-s.Rational(5,6))*exp(-6*y/5),
        5+(p-5)*exp(-y/5)]
    eh_y=eh*exp(-8*y/5)
    em_y=E+(em-E)*exp(-y)
    ek_y=ek*exp(-8*y/5)+E*(1-exp(-8*y/5))/s.Rational(8,5)
    aa_y=E**2+(aa-E**2)*exp(-y)
    right=[s.Rational(5,8)+eh_y,4*z*(s.Rational(5,8)+eh_y)+ek_y,
        4*z+em_y,16*z*z+8*z*em_y+aa_y,
        s.Rational(5,6)+eb*exp(-6*y/5),5+ep*exp(-y/5)]
    differences=[s.simplify(a-bb) for a,bb in zip(left,right)]
    if any(q!=0 for q in differences):
        raise ArithmeticError("Same Rsh moment function failed centered coordinate identity")
    # The flat cutoff's ordinary phase derivatives d1..d4 vanish at1.
    # T is a source constant independent of Z; 1/T factors stay in the jet.
    ds=s.symbols("sigma1:5",real=True)
    cutoff=1+sum(ds[k-1]*y**k/(T**k*s.factorial(k)) for k in range(1,5))
    left_log=B*(1-cutoff)-s.log(1+z*z)+(T+y)/10-C-P
    right_log=-s.log(1+z*z)+T/10-C-P+y/10
    logrows=[]
    for k in range(5):
        for n in range(5-k):
            q=s.simplify(s.diff(left_log-right_log,y,k,z,n).subs(y,0).subs(dict.fromkeys(ds,0)))
            if q!=0:
                raise ArithmeticError("Flat original log-amplitude boundary jet mismatch")
            logrows.append("y"+str(k)+"_Z"+str(n))
    def physical(shapes):
        hh,kk,mm,aaa,bb,pp=shapes
        radius=R0*exp(y);u=U*exp(y/10)
        Q=(2*z*V-(1-delta)*z*mm-(1-z*z)*s.diff(mm,z))/(1-delta*z*z)
        return dict(Utheta=u,Uz=V,Ur=s.sqrt(radius/2)*Q,
            P_over_Pstar2=P0+u*u*pp/(2*Pstar**2),
            Mtheta=s.sqrt(2)*radius**s.Rational(3,2)*u*hh,
            Mtheta_z=s.sqrt(2)*radius**s.Rational(3,2)*u*kk,
            Mz=radius*mm,Mztheta_over_Pstar2=radius*(aaa-u*u*bb/2)/Pstar**2,
            Mp_over_Pstar2=u*u*pp/(2*Pstar**2))
    fields_left=physical(left);fields_right=physical(right);rows={}
    for name,expression in fields_left.items():
        difference=s.simplify(expression-fields_right[name])
        if difference!=0:
            raise ArithmeticError("Same physical source function differs: "+name)
        rows[name]=[]
        for k in range(5):
            for n in range(5-k):
                if s.diff(difference,y,k,z,n).subs(y,0)!=0:
                    raise ArithmeticError("Physical Rsh mixed row differs")
                rows[name].append("y"+str(k)+"_Z"+str(n))
    return dict(arbitrary_smooth_axial_inlet_functions=True,
        exact_six_centered_history_function_identities=6,
        exact_original_flat_log_mixed4_boundary_rows=len(logrows),
        physical_boundary_source_identities=list(fields_left),
        physical_mixed4_rows_implied_by_exact_function_identities=rows,
        total_physical_mixed4_rows_implied=sum(map(len,rows.values())),
        original_source_rates=["8/5","8/5","1","1","6/5","1/5"],
        original_source_forcings=["1","V","V","V^2","1","1"],
        pressure_is_Pstar_squared_times_same_axis_datum_plus_Mp=True,
        Q_keeps_axial_derivative_of_same_mean_and_original_delta=True,
        constant_power_extension_used_only_for_boundary_jet=True,
        finite_left_neighborhood_replaced=False,
        same_fixed_Rsh_Utheta_Pstar_normalization_on_both_sides=True,
        normalization_is_frozen_after_physical_differentiation=True,
        interval_overlap_is_not_functional_proof=True,passed=True)


@source_precision
def build():
    provider=CompliantActualReferenceRestoreMixedC4()
    left=provider.long_mixed;right=provider.reference
    family,source=provider.family,provider.source
    hashes=dict(provider.hashes);records={}
    specifications=(
        ("actual_long_reshape_mixed_C4_check","current_actual_long_reshape_mixed4_available"),
        ("actual_reference_restore_mixed_C4_check","current_actual_reference_restore_mixed4_available"),
    )
    for part,flag in specifications:
        name=PREFIX+part+".json";record=accepted(name,family,source,flag)
        if record["datum_enclosure_sha256"]!=right.core.datum.datum_sha:
            raise ValueError("Current Rsh pressure datum differs")
        hashes.update(record["input_hashes"]);hashes[name]=sha(name);records[part]=record
    long_name=PREFIX+"actual_long_reshape_mixed_C4.json"
    ref_name=PREFIX+"actual_reference_restore_mixed_C4.json"
    long=json.loads((HERE/long_name).read_bytes());ref=json.loads((HERE/ref_name).read_bytes())
    _verify_hashes(long);_verify_hashes(ref)
    hashes[long_name]=sha(long_name);hashes[ref_name]=sha(ref_name)
    if (right.reshape is not left.reshape or right.history is not left.history
            or right.core is not left.reshape.core or provider.invP2._mpi_!=left.invP2._mpi_
            or provider.shared_axial_source!=left.shared_axial_source):
        raise ValueError("Rsh two-sided provider/physical normalization identity lost")
    if (long["shared_exact_axial_source"]!=ref["shared_exact_axial_source"]
            or ref["shared_exact_axial_source"]!=provider.shared_axial_source
            or ref["current_centered_E_source_bindings"]!=current_E_source_bindings(provider)):
        raise ValueError("Current Rsh V and centered E source graph differs")
    graph=provider.shared_axial_source
    if graph["E_V110_minus_4Z"]["args"]!=graph["V110"]["args"][1:]:
        raise ValueError("Current exact E is not V110 minus exact4Z")
    parent=long["actual_Rsh_exit"]["actual_inherited_axial5_packet"]
    if canonical_source(parent)!=canonical_source(ref["current_actual_Rsh_parent"]):
        raise ValueError("Current Rsh reference parent is not current long-reshape parent")
    rp=ref["actual_Rsh_exit"]["actual_inherited_axial5_packet"]
    if canonical_source(parent["pressure_axis_axial5_coefficients"])!=canonical_source(rp["pressure_axis_axial5_coefficients"]):
        raise ValueError("Current physical pressure datum transfer differs")
    long_generic=left.original_check;ref_generic=provider.original_check
    if (not long_generic["symbolic_checks"]["passed"] or not long_generic["independent_physical_fixture"]["passed"]
            or not ref_generic["structural_checks"]["passed"] or not ref_generic["independent_physical_fixture"]["passed"]):
        raise ValueError("Hash-current original physical operator/fixture admission missing")
    jets=sigma_jets(provider.ctx,provider.ctx.mpf(1))
    if endpoints(jets[0])!=(mp.mpf(1),mp.mpf(1)) or any(endpoints(jets[k])!=(mp.mpf(0),mp.mpf(0)) for k in range(1,5)):
        raise ArithmeticError("Original sigma phase-one flat C4 jet changed")
    if endpoints(left.reshape.T)!=endpoints(400*left.reshape.A):
        raise ValueError("Original source T=400A changed")
    native_flags=dict(long=long["Rsh_reference_mixed4_join_certified"],
        reference=ref["current_Rsh_source_functional_join_certified"])
    if any(native_flags.values()):
        raise ValueError("Dedicated new source-join receipt expected; native artifacts stay scoped")
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    return dict(actual_five_defect_family_sha256=family,implicit_source_sha256=source,
        datum_enclosure_sha256=right.core.datum.datum_sha,
        shared_exact_axial_source_namespace=graph["shared_source_namespace"],
        current_Rsh_source_bindings=source_bindings(),boundary_source_proof=boundary_identities(),
        current_two_sided_provider_and_core_identity_verified=True,
        current_Rsh_parent_and_axis_pressure_identical=True,
        current_V_E_common_signed_source_graph_verified=True,
        unchanged_physical_operator_fixtures_reused_by_current_hashes=True,
        generic_fixture_receipts=[PREFIX+"long_reshape_mixed_C4_check.json",PREFIX+"reference_restore_mixed_C4_check.json"],
        original_flat_cutoff_endpoint_derivatives_checked=5,
        native_historical_Rsh_flags_preserved=native_flags,
        current_Rsh_source_functional_join_certified=True,
        **dict.fromkeys(SCOPES,False),input_hashes=hashes)


def run():
    result=build()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("Current Rsh exact source mixed4 join generated: 6 histories, 9 fields, 135 implied rows",flush=True)
    return result


if __name__=="__main__":
    run()

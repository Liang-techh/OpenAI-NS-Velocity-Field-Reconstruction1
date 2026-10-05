"""Current correlated E and reference/restoration mixed4 through Rm.

E is bounded from its common core and actual source increments, never by
subtracting independent V110/4Z boxes or selecting cap representatives.
Original transport, restoration kernels and physical source equations stay.
"""
import ast
import json
import math
from pathlib import Path

import mpmath as mp

import lei_ren_part1_paper_compliant_reference_restore_profiles as profile_module
import lei_ren_part1_paper_compliant_reference_restore_mixed_C4 as mixed_module
from lei_ren_part1_paper_compliant_actual_long_reshape_mixed_C4 import (
    CompliantActualLongReshapeMixedC4, accepted, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square, symmetric
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_five_moment_repair import pack

PROFILE=profile_module.CompliantReferenceRestoreProfiles
MIXED=mixed_module.CompliantReferenceRestoreMixedC4



def current_centered_core(history,Z):
    """Same fresh core recurrence with 4Z removed before interval summation."""
    comparison=history.upstream.comparison
    c=comparison.ctx;core=comparison.core
    packet=comparison.atoms.rebuild.rebuild(Z,24,6)
    baseline=IntervalTaylor.variable(c,packet["Z"],5)*4+core.j
    if len(packet["rows"]["Uz"][0])<6:
        raise ValueError("Fresh core baseline axial5 required")
    for k in range(6):
        if c.mpf(packet["rows"]["Uz"][0][k])._mpi_!=baseline[k]._mpi_:
            raise ValueError("Fresh core row zero is not the exact affine baseline")
    centered_packet=dict(packet)
    centered_rows=dict(packet["rows"])
    centered_rows["Uz"]=[list(row) for row in packet["rows"]["Uz"]]
    centered_rows["Uz"][0]=[c.mpf(core.j)]+[c.mpf(0)]*(len(centered_rows["Uz"][0])-1)
    centered_packet["rows"]=centered_rows
    E=comparison.source_profile_jet(centered_packet,4,0,False)[1].truncate(5)
    return dict(E=E,baseline_rows_checked=6,
        source="current fresh radial recurrence; exact 4Z removed from radial row0 before finite sum",
        tail_source="same current source_profile_jet nonlinear epsilon*Psi tail, unchanged",
        defining_field_not_changed=True)


def fresh_core_source_bindings():
    seed=assignment_source_bindings("core_coefficient_rebuild","seed",{"u":"4*z+self.core.j"})
    inlet=assignment_source_bindings("comparison_point_integrals","inlet",{
        "packet":"self.atoms.field.build_root_rows(24,6) if root else self.atoms.rebuild.rebuild(Z,24,6)",
        "local":"[self.source_profile_jet(packet,4,i,root) for i in range(3)]",
        "(phi0, V0)":"local[0]",
    })
    tree=ast.parse((HERE/"lei_ren_part1_paper_functional_core_step.py").read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=="initial_rows")
    dictionaries=[n for n in ast.walk(fn) if isinstance(n,ast.Dict)
                  and any(isinstance(k,ast.Constant) and k.value=="Uz" for k in n.keys)]
    if len(dictionaries)!=1:
        raise ValueError("Unique original fresh Uz row-zero assignment required")
    node=dictionaries[0]
    value=next(v for k,v in zip(node.keys,node.values) if isinstance(k,ast.Constant) and k.value=="Uz")
    expected=ast.parse('[_coerce_vector(ctx,fixed["U0_Z_taylor"])]',mode="eval").body
    if ast.dump(value)!=ast.dump(expected):
        raise ValueError("Fresh baseline row-zero definition changed")
    seed_tree=ast.parse((HERE/(PREFIX+"core_coefficient_rebuild.py")).read_text(encoding="utf8"))
    seed_fn=next(n for n in ast.walk(seed_tree) if isinstance(n,ast.FunctionDef) and n.name=="seed")
    values=[kw.value for n in ast.walk(seed_fn) if isinstance(n,ast.Call) for kw in n.keywords if kw.arg=="U0_Z_taylor"]
    if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse("list(u.coefficients)",mode="eval").body):
        raise ValueError("Fresh fixed U0 no longer comes from the same affine source")
    centered=assignment_source_bindings("actual_reference_restore_mixed_C4","current_centered_core",{
        "comparison":"history.upstream.comparison",
        "packet":"comparison.atoms.rebuild.rebuild(Z,24,6)",
        "centered_packet":"dict(packet)",
        "centered_rows":"dict(packet['rows'])",
        "centered_rows['Uz']":"[list(row) for row in packet['rows']['Uz']]",
        "centered_rows['Uz'][0]":"[c.mpf(core.j)]+[c.mpf(0)]*(len(centered_rows['Uz'][0])-1)",
        "centered_packet['rows']":"centered_rows",
        "E":"comparison.source_profile_jet(centered_packet,4,0,False)[1].truncate(5)",
    })
    return dict(seed_AST_bindings=seed,inlet_AST_bindings=inlet,
        centered_current_packet_AST_bindings=centered,
        fixed_U0_and_initial_rows_AST_bound=True,
        centered_core_uses_same_current_finite_rows_and_tail=True)


def current_E_source_bindings(provider):
    bridge=assignment_source_bindings("actual_bridge_integrals","packet",{
        "V0":"p['data']['V0'].truncate(5)",
        "V":"V0+delta_V",
    })
    switch=assignment_source_bindings("inner_switch_profiles","inputs",{
        "v":"jet(incoming['Uz_actual_axial5_coefficients'])",
        "vc":"v+IntervalTaylor(c,inc)",
    })
    recipe=assignment_source_bindings("actual_reference_restore_mixed_C4","inputs",{
        "core_packet":"current_centered_core(self.history,Z)",
        "core_E":"core_packet['E']",
        "incoming":"self.history.upstream.packet(Z,1,'macro')",
        "bridge_delta":"IntervalTaylor(c,incoming['actual_delta_V_axial5'])",
        "first_delta":"IntervalTaylor(c,switch['velocity_increment'])",
        "E":"core_E+bridge_delta+first_delta",
    })
    restore=assignment_source_bindings("reference_restore_profiles","restoration",{
        "V":"inp['z']*4+inp['E']*alpha",
        "centered":"restore_centered(c,initial,inp['E'],t,kernels)",
    })
    checks={name:getattr(CompliantActualReferenceRestoreProfiles,name) is getattr(PROFILE,name)
            for name in ("reference_centered","packet","reference","kernels","restoration","terminal","defects")}
    checks.update({name:getattr(CompliantActualReferenceRestoreMixedC4,name) is getattr(MIXED,name)
                   for name in ("packet","reference_branch","restoration","postrestore")})
    if not all(checks.values()):
        raise ValueError("Original reference/restoration algorithms changed")
    return dict(current_bridge_increment_AST_bindings=bridge,
        original_first_switch_increment_AST_bindings=switch,
        correlated_current_E_recipe_AST_bindings=recipe,
        original_restore_source_AST_bindings=restore,
        unchanged_original_callables=checks,
        current_fresh_core_source_bindings=fresh_core_source_bindings(),
        same_core_bridge_first_switch_defining_E_source=True)


class CompliantActualReferenceRestoreProfiles(PROFILE):
    @source_precision
    def __init__(self,long_mixed):
        super().__init__()
        prior=(self.family,self.source,self.core.datum.datum_sha)
        self.long_mixed=long_mixed;self.reshape=long_mixed.reshape
        self.history=long_mixed.history;self.core=self.reshape.core
        self.ctx=c=self.reshape.ctx;self.family=self.reshape.family;self.source=self.reshape.source
        if prior!=(self.family,self.source,self.core.datum.datum_sha):
            raise ValueError("Same core/family/pressure datum required")
        name=PREFIX+"actual_long_reshape_mixed_C4_check.json"
        self.current_long_check=accepted(name,self.family,self.source,
                                       "current_actual_long_reshape_mixed4_available")
        if self.current_long_check["datum_enclosure_sha256"]!=self.core.datum.datum_sha:
            raise ValueError("Current Rsh pressure datum changed")
        self.hashes.update(long_mixed.hashes)
        self.hashes.update(self.current_long_check["input_hashes"]);self.hashes[name]=sha(name)
        self.cache={};self.kernel_cache={};self.proofs=[]
        self.loggap=self.reshape.logref-self.reshape.T
        if endpoints(self.loggap)[0]<=8:
            raise ValueError("Original Rsh-to-reference ordering changed")
        self.logcap=-1000*c.ln(10)-2*self.core.logP;self.tailcap=c.exp(self.logcap)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    @source_precision
    def inputs(self,Z):
        c=self.ctx;Z=c.mpf(Z);key=Z._mpi_
        if key in self.cache:
            return self.cache[key]
        parent=self.reshape.evaluate(Z,phase=1);jet=lambda row:IntervalTaylor(c,row)
        moments={n:jet(row) for n,row in parent["actual_normalized_moment_shape_axial5_coefficients"].items()}
        z=IntervalTaylor.variable(c,Z,5)
        centered=dict(mean_error=moments["mean"]-z*4,
            angular_error=moments["theta"]-c.mpf(5)/8,
            mixed_error=moments["theta_z"]-(z*moments["theta"])*4,
            axial_square=moments["axial"]-(z*moments["mean"])*8+square(z)*16,
            swirl_error=moments["swirl"]-c.mpf(5)/6,pressure_error=moments["pressure"]-5)
        row=centered["axial_square"];lower,upper=endpoints(row[0])
        centered["axial_square"]=IntervalTaylor(c,[c.mpf([max(mp.mpf(0),lower),upper])]+list(row.coefficients[1:]))
        core_packet=current_centered_core(self.history,Z)
        core_E=core_packet["E"]
        incoming=self.history.upstream.packet(Z,1,"macro")
        bridge_delta=IntervalTaylor(c,incoming["actual_delta_V_axial5"])
        switch=self.history.switch.inputs(Z)
        first_delta=IntervalTaylor(c,switch["velocity_increment"])
        E=core_E+bridge_delta+first_delta
        raw_E=E
        ledger=self.reshape.switch.bridge.records["K1_ledger"]
        gate=self.reshape.switch.bridge.records["global_exit_certificate"]
        rho_bridge=read_interval(c,gate["rho_bridge_C2_upper"])
        delta=bridge_delta+first_delta
        delta_C2=sum((c.mpf(max(abs(v) for v in endpoints(delta[k])))*math.factorial(k)
                      for k in range(3)),c.mpf(0))
        if endpoints(delta_C2)[1]>endpoints(rho_bridge)[0]:
            raise ValueError("Current bridge/first-switch C2 enclosure exceeds source theorem")
        if (not gate["actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified"]
                or not gate["short_switch_moments_included_by_exact_velocity_averaging"]
                or "includes both short-switch C2 velocity increments" not in gate["proof"]["axial_budget"]):
            raise ValueError("Same original global exit theorem must include both switches")
        rho=read_interval(c,ledger["rho_core_C2_bound"])+rho_bridge
        coefficients=[]
        for k,q in enumerate(E.coefficients):
            if k<=2:
                bound=symmetric(c,rho)/math.factorial(k)+(self.core.j if k==0 else 0)
                q=intersection(c,q,bound)
            coefficients.append(q)
        E=IntervalTaylor(c,coefficients)
        result=dict(parent=parent,z=z,E=E,rho=rho,source_centered_Rsh=centered,
            original_axis_pressure=jet(parent["pressure_axis_axial5_coefficients"]),
            current_E_component_enclosures=dict(core=list(core_E.coefficients),
                bridge=list(bridge_delta.coefficients),first_switch=list(first_delta.coefficients),
                raw_sum_before_C2_theorem=list(raw_E.coefficients)),
            current_E_source_function=self.long_mixed.shared_axial_source["E_V110_minus_4Z"],
            current_E_increment_C2_proof=dict(ordinary_derivative_sum_upper=delta_C2,
                admitted_bridge_and_both_switches_C2_upper=rho_bridge,
                verified_before_E_intersection=True))
        self.cache[key]=result
        return result

    @source_precision
    def report(self):
        result=super().report()
        result["centered_axial_source_proof"]=dict(
            exact_source=self.long_mixed.shared_axial_source["E_V110_minus_4Z"],
            component_enclosures=self.inputs([-1,1])["current_E_component_enclosures"],
            C2_error_about_j_upper=self.inputs([-1,1])["rho"],
            caps_do_not_define_E=True,independent_baseline_subtraction_not_used=True)
        result.update(current_actual_long_reshape_and_E_used=True,
                      actual_point_moment_history_recovered=False,
                      full_implicit_leading_inputs_recomputed=False)
        return result


class CompliantActualReferenceRestoreMixedC4(MIXED):
    @source_precision
    def __init__(self):
        self.long_mixed=CompliantActualLongReshapeMixedC4()
        self.reference=CompliantActualReferenceRestoreProfiles(self.long_mixed)
        self.ctx=c=self.reference.ctx;self.family=self.reference.family;self.source=self.reference.source
        self.hashes=dict(self.reference.hashes)
        name=PREFIX+"reference_restore_mixed_C4_check.json"
        self.original_check=accepted(name,self.family,self.source,"actual_reference_restore_mixed4_available")
        self.hashes.update(self.original_check["input_hashes"]);self.hashes[name]=sha(name)
        self.invP2=c.exp(-2*self.reference.core.logP)
        self.shared_axial_source=self.long_mixed.shared_axial_source
        self.source_bindings=current_E_source_bindings(self)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    @source_precision
    def report(self):
        result=super().report()
        inp=self.reference.inputs([-1,1])
        result.update(datum_enclosure_sha256=self.reference.core.datum.datum_sha,
            current_centered_E_source_bindings=self.source_bindings,
            current_centered_E_component_enclosures=inp["current_E_component_enclosures"],
            current_centered_E_axial5=list(inp["E"].coefficients),
            current_E_increment_C2_proof=inp["current_E_increment_C2_proof"],
            shared_exact_axial_source=self.shared_axial_source,
            current_actual_Rsh_parent=inp["parent"],
            current_actual_reference_restore_mixed4_available=True,
            current_correlated_E_installed=True,
            current_original_restore_end_exact_4Z=True,
            current_Rsh_source_functional_join_certified=False,
            current_actual_moment_patch_installed=False,
            full_implicit_leading_inputs_recomputed=False,
            actual_point_moment_history_recovered=False,
            global_completed_tensor_admissibility=False,
            current_unpatched_five_defects_at_Rm=self.reference.defects([-1,1],-6))
        return result


@source_precision
def run():
    result=CompliantActualReferenceRestoreMixedC4().report()
    Path(__file__).with_suffix(".json").write_text(json.dumps(encode(pack(result)),indent=2)+"\n",encoding="utf8")
    print("Current correlated E and original reference/restoration mixed4 generated through Rm",flush=True)
    return result


if __name__=="__main__":
    run()

"""Bind actual physical source ports without promoting undeclared derivatives.

The accepted source-germ theorem supplies U through Z2 and U_y through Z1.
The actual affine physical mean recipes are checked, but no modified-m C2
ordinary-Z provenance receipt has been supplied. Those four obligations stay
open; a FunctionJet array or C5 declaration is not a new derivative proof.
"""
from __future__ import annotations

import ast
from collections.abc import Mapping, Sequence
from copy import deepcopy
import hashlib
import json

SCHEMA = "original_Rm_physical_primitive_Z_jets_v1"
GRAPH_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
PREFIX_DIGEST = "c7ab7da8cc5c51cd782a5794e3076a6cf021a45235aae1c75bd246b01491b639"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
RECIPES = {
    PREFIX+"modified_Q4_pressure_C3_functions.py":"7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    PREFIX+"true_Rm_patch_finite_N_C5_functions.py":"c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    PREFIX+"core_inlet_axial6_generic_source.py":"29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890",
}
RECEIPT_NAME = PREFIX+"Rm_outer_mixed_source_joins_check.json"
RECEIPT_SHA = "14d1cc6b641435c60a7e3662c996a260c07b328df5b8f486c4344f3d7eeb4867"
SOURCE_FAMILY = {
    "actual_five_defect_family_sha256":"3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256":"5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256":"dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
CONTEXTS = {"postrestore":dict(coordinate=84528,endpoint=1), "patch":dict(coordinate=185676,endpoint=1)}
PORTS = {
    "postrestore":dict(U=[183268,183269,183270],U_y=[183274,183275,183276],m=[183280,183281,183282],
        m_y=[183292,183293,183294],V=[121571,121572,121352],V_y=[182757,182758,182759],
        modified_m=[121720,121721,121722],negative_m=[183286,183287,183288]),
    "patch":dict(U=[204317,204318,204319],U_y=[204323,204324,204325],m=[204329,204330,204331],
        m_y=[204341,204342,204343],V=[198423,198424,198425],V_y=[203806,203807,203808],
        modified_m=[198591,198592,198593],negative_m=[204335,204336,204337]),
}


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=True,allow_nan=False).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def node(nodes, index):
    require(type(index) is int and 0 <= index < 205147, "source handle outside immutable finite prefix")
    row = nodes[index]
    require(isinstance(row,Mapping), "source node is not a mapping")
    return row


def exact_row(nodes,index,operation,**fields):
    row=node(nodes,index)
    require(row.get("operation")==operation and all(row.get(k)==v for k,v in fields.items()),
        "actual affine source recipe changed at "+str(index))
    return dict(node=index,operation=operation,checked_fields=deepcopy(fields))


def read_record(records,name,expected):
    record=records.get(name)
    require(isinstance(record,Mapping) and record.get("path")=="experiments/root_st073/"+name,
        "required source record absent: "+name)
    text=record.get("source_text")
    require(isinstance(text,str) and record.get("sha256")==expected
        and hashlib.sha256(text.encode("utf-8")).hexdigest()==expected,"source record hash changed: "+name)
    return text


def source_recipes(records):
    require(isinstance(records,Mapping),"source_recipe_records must be a mapping")
    witnesses={}
    for name,expected in RECIPES.items():
        text=read_record(records,name,expected)
        tree=ast.parse(text,filename=name)
        wanted=("physical_pressure_fields",) if "modified_Q4" in name else (("CurrentOriginalTrueRmPatchFiniteNC5Functions",) if "finite_N" in name else ("JetAlgebra",))
        found={n.name:n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in wanted}
        require(set(found)==set(wanted),"source recipe declaration absent: "+name)
        witnesses[name]=dict(sha256=expected,ast_definitions={k:canonical_digest(ast.dump(v,include_attributes=False)) for k,v in found.items()})
    return witnesses


def admitted_velocity_source(records,actual_ports):
    data=json.loads(read_record(records,RECEIPT_NAME,RECEIPT_SHA))
    require(data.get("all_passed") is True and data.get("source_family")==SOURCE_FAMILY,
        "accepted C5 velocity source theorem missing or changed")
    checks=data.get("checks",{})
    for name in ("source_contract","original_radius_and_phase","actual_flat_controls","actual_background_germs",
        "actual_memory_and_leading_inlets","actual_own_mixed_joins"):
        require(checks.get(name,{}).get("passed") is True,"accepted source theorem check missing: "+name)
    memory=checks["actual_memory_and_leading_inlets"]
    finite=memory["actual_finite_N_endpoint_source_germs"]
    calculus=memory["actual_U_source_radial_derivative_calculus"]
    for key in ("root_C5_Z_and_physical_y_germs_from_actual_background_IVP_and_original_generic_operator",
        "original_ordinary_calculus_and_physical_J_inverse_chain_rule_then_supply_U_y0_to4",
        "actual_A_B_inverse_integrals_and_velocity_outputs_not_assigned_common_values"):
        require(finite.get(key) is True,"accepted finite velocity source rule missing: "+key)
    for key in ("actual_U_y_first_row_uses_original_slow_partial_plus_fast_phase_recovery",
        "actual_U_y2_to4_are_derivatives_of_actual_preceding_rows_with_variable_J",
        "accepted_unique_inverse_composition_and_equal_background_germs_then_identify_high_U_germs"):
        require(calculus.get(key) is True,"accepted physical velocity calculus rule missing: "+key)
    bound=[]
    source_rows=calculus["actual_U_y0_to4_C5_source_binding_rows"]
    u_rows=finite["actual_inverse_loop_and_finite_N_source_pairs"]
    for side,label in (("postrestore","postrestore"),("patch","Rm_patch")):
        ports=actual_ports[side]
        for order in (0,1):
            expected=[label,1,order,ports["V_y"][order],ports["U_y"][order]]
            require(source_rows.count(expected)==1,"actual Uy source theorem port absent/duplicated")
            bound.append(expected)
    for order in (0,1,2):
        expected=["postrestore_patch","physical","U",order,actual_ports["postrestore"]["U"][order],actual_ports["patch"]["U"][order]]
        require(u_rows.count(expected)==1,"actual U source theorem port absent/duplicated")
        bound.append(expected)
    return dict(receipt_sha256=RECEIPT_SHA,source_family=deepcopy(SOURCE_FAMILY),actual_source_rows=bound,
        U_through_ordinary_Z2_source_bound=True,Uy_through_ordinary_Z1_source_bound=True,
        modified_m_source_derivative_not_supplied=True,old_theorem_or_checker_replayed=False)


def derive(nodes,actual_source_bindings,source_recipe_records,*,expected_nodes_digest):
    require(expected_nodes_digest==GRAPH_DIGEST and len(nodes)==205300 and canonical_digest(nodes)==GRAPH_DIGEST,
        "actual source graph changed")
    require(canonical_digest(nodes[:205147])==PREFIX_DIGEST,"immutable source prefix changed")
    require(isinstance(actual_source_bindings,Mapping) and set(actual_source_bindings)==set(CONTEXTS),"actual chart sides changed")
    require(node(nodes,20)==dict(operation="bound_variable",name="Z"),"ordinary Z source variable changed")
    exact_row(nodes,33572,"analytic_unary",name="exp",argument=5)
    recipe_witnesses=source_recipes(source_recipe_records)
    affine=[]
    for side,p in PORTS.items():
        binding=actual_source_bindings[side]
        require(isinstance(binding,Mapping),"invalid actual chart binding")
        physical=binding.get("physical",{});radial=binding.get("radial",{})
        for field in ("U","U_y","m","m_y"):
            require(isinstance(physical.get(field),Sequence) and list(physical[field][:3])==p[field],"physical jet ports changed: "+side+"."+field)
        require(list(radial.get("V_y",[])[:3])==p["V_y"],"physical radial Vy source jets changed")
        chart=binding.get("chart",{})
        if chart:
            require(chart.get("coordinate")==CONTEXTS[side]["coordinate"] and chart.get("physical_measure")==187215,
                "actual finite chart coordinate/J changed")
        for order in (0,1,2):
            affine.extend((exact_row(nodes,p["U"][order],"product",arguments=[p["V"][order],33572]),
                exact_row(nodes,p["U_y"][order],"product",arguments=[p["V_y"][order],33572]),
                exact_row(nodes,p["m"][order],"product",arguments=[p["modified_m"][order],33572]),
                exact_row(nodes,p["negative_m"][order],"negative",argument=p["m"][order]),
                exact_row(nodes,p["m_y"][order],"sum",arguments=[p["U"][order],p["negative_m"][order]])))
    velocity=admitted_velocity_source(source_recipe_records,PORTS)
    rows=[];gates=[]
    for side,p in PORTS.items():
        rows.append(dict(side=side,field="U_y",ordinary_order=0,source_handle=p["U_y"][0],derivative_handle=p["U_y"][1],
            source_derivative_proved=True,regularity_proved=True,context=deepcopy(CONTEXTS[side]),
            proof_origin="admitted C5 ordinary-Z/physical-y source germ theorem, actual bound Uy rows",
            source_receipt_sha256=RECEIPT_SHA))
        for order in (0,1):
            dependency=dict(side=side,source_handle=p["modified_m"][order],derivative_handle=p["modified_m"][order+1],
                physical_source_handle=p["m"][order],physical_derivative_handle=p["m"][order+1],
                missing="actual ordinary-D_Z source/regularity receipt for modified m including its history")
            rows.append(dict(side=side,field="physical_mean_y",ordinary_order=order,source_handle=p["m_y"][order],
                derivative_handle=p["m_y"][order+1],source_derivative_proved=False,regularity_proved=False,
                context=deepcopy(CONTEXTS[side]),affine_source_recipe_checked=True,
                velocity_U_source_jet_backed=True,modified_m_dependency=dependency))
            gates.append(dict(gate="modified_m_ordinary_Z_source_receipt_missing",ordinary_order=order,**dependency))
    return dict(schema=SCHEMA,passed=False,source_only=True,nodes_digest=GRAPH_DIGEST,Z_variable=20,
        contexts=deepcopy(CONTEXTS),rows=rows,source_affine_recipes_checked=True,affine_source_witness=affine,
        admitted_velocity_source=velocity,producer_source_recipe_witnesses=recipe_witnesses,
        source_derivative_rows_proved=2,source_derivative_rows_open=4,unresolved_gates=gates,
        derivative_lift_authorized=False,numeric_values_computed=False,differentiated_opaque_alias=False,
        old_operators_replayed=False,full_RQy_germ_proved=False,full_goal_complete=False,
        input_data_digests=dict(actual_source_bindings=canonical_digest(actual_source_bindings),
            source_recipe_records=canonical_digest(source_recipe_records)),
        limits="Affine jet construction is checked; ordinary modified-m derivative/regularity is not inferred from arrays or C5 declarations.")

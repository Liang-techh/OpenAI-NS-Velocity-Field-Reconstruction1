"""Signed leading-viscous source interfaces over the saved 183-term ledger.

Only new source operators are evaluated.  Accepted numerical producers are
not reconstructed, and inherited error rows remain in the active ledger.
"""
from __future__ import annotations

import gzip
import ast
import hashlib
import json
from pathlib import Path

from .original_core_numeric import ROOT, DIRECTORY, PREFIX

PARENT = PREFIX + "compliant_current_original_macro_endpoint_m15_numeric.json"
PARENT_SHA256 = "16eb888da56fad8e5696b5fd8400b827ff7634e2c8f524110267cfcc6f525b80"


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()


class OriginalLeadingViscousM15:
    def __init__(self, *, directory=DIRECTORY):
        self.directory = Path(directory)
        raw = (self.directory/PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest() != PARENT_SHA256:
            raise ValueError("Accepted macro endpoint report required")
        self.parent = json.loads(raw)
        if self.parent["checks"]["passed"] is not True:
            raise ValueError("Passed parent report required")
        self.family = self.parent["source_family"]
        self.input_hashes = dict(self.parent["input_hashes"])
        self.input_hashes[PARENT] = PARENT_SHA256
        self.implementation_hashes = dict(self.parent["implementation_hashes"])
        for base, group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for relative, expected in group.items():
                path = base/relative
                if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                    raise ValueError("Changed inherited source pin: "+relative)
        self.high_archive = self.archive("prescribed_high_radial_source_functions")
        self.nodes = self.high_archive["exact_graph_nodes"]
        self.loop_archive = self.archive("prescribed_loop_C5_functions")
        self.modified_archive = self.archive("modified_Q4_pressure_C3_functions")
        self.prefix_checks = {}
        for name, packet in (("loop",self.loop_archive),("modified",self.modified_archive)):
            nodes = packet["exact_graph_nodes"]
            if self.nodes[:len(nodes)] != nodes:
                raise ValueError("Identical source graph prefix required")
            self.prefix_checks[name] = {"passed":True,"node_count":len(nodes)}

    def archive(self, suffix):
        name = PREFIX+"compliant_current_original_"+suffix+".json.gz"
        if name not in self.input_hashes:
            raise ValueError("Missing inherited archive pin: "+name)
        with gzip.open(self.directory/name,"rt",encoding="utf-8") as stream:
            return json.load(stream)

    def bind_pins(self, packet):
        for row in packet.get("source_pins",{}).values():
            if not isinstance(row,dict) or row.get("passed") is not True:
                raise ValueError("Passed new source pin required")
            relative = row["path"]
            path = ROOT/relative
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != row["sha256"]:
                raise ValueError("Changed new source pin: "+relative)
            if relative.startswith("experiments/root_st073/"):
                target = self.input_hashes; key = path.name
            else:
                target = self.implementation_hashes; key = relative
            if target.get(key,actual) != actual:
                raise ValueError("New source pin conflicts with parent")
            target[key] = actual

    def preservation(self):
        samples = []
        for row in self.parent["samples"]:
            raw = row["raw_m15_C0_C1_terms"]
            reduced = row["correlated_reduced_m15_C0_C1_terms"]
            if len(raw)!=213 or len(reduced)!=183:
                raise ValueError("Inherited raw213/reduced183 packet required")
            if sum(term.get("name")=="reserved_n1_own_swirl_once" for term in reduced)!=1:
                raise ValueError("Reserved contribution must remain exactly once")
            samples.append({"Z":row["Z"],"raw_count":213,"reduced_count":183,
                            "raw_sha256":digest(raw),"reduced_sha256":digest(reduced),
                            "signed_inner_lower_six_sha256":digest(raw[:6]),
                            "independent_target_sha256":digest(row["independent_target"]),
                            "reserved_count":1,"terms_replaced_or_removed":0})
        return {"passed":True,"report":PARENT,"report_sha256":PARENT_SHA256,
                "samples":samples,"ledger_hydrated_by_reference":True,
                "new_source_interfaces_do_not_change_error_budget":True}

    def leading_ftc(self):
        """Bind the leading FTC to independent source rows, with no graph alias."""
        high = self.high_archive["prescribed_high_radial_source_bindings"]
        loop = self.loop_archive["prescribed_loop_C5_bindings"]
        pstar = high["downstream"]["Pstar"]
        parameters = self.parent["source_seam_proofs"]["radius_parameter_guard"]
        if not all(parameters.get(key) is True for key in (
                "passed","ordinary_R_Z_derivatives_zero","Pstar_Z_independent","delta_Z_independent")):
            raise ValueError("Independent axial-constant pressure scales required")
        inverse = self.nodes[69245]
        if (inverse.get("operation")!="positive_quotient" or inverse.get("numerator")!=1
                or inverse.get("denominator")!=pstar):
            raise ValueError("Independent leading normalization 1/Pstar required")
        charts = {}
        for name in ("second_micro","frozen_macro","first_switch"):
            radius = self.parent["native_source_guards"]["chart_guards"][name]["radius_phase_J"]
            if radius.get("d_native_log_R_equals_J") is not True:
                raise ValueError("Actual independent physical derivative R_y=R required")
            bg = high["backgrounds"][name]
            arrays = {"V":bg["fields"]["V"][:3],"m":bg["fields"]["m"][:3],
                      "Vy":bg["radial"]["V_y"][:3]}
            ode_rows = []
            normalization = []
            for order in range(3):
                actual = high["background_y"][name]["m"][order]
                row = self.nodes[actual]
                negatives = [i for i in row.get("arguments",()) if
                             self.nodes[i].get("operation")=="negative"
                             and self.nodes[i].get("argument")==arrays["m"][order]]
                valid = (row.get("operation")=="sum" and len(negatives)==1
                         and sorted(row["arguments"])==sorted([arrays["V"][order],negatives[0]]))
                if not valid:
                    raise ValueError("Independent leading m_y=V-m source row changed")
                ode_rows.append({"ordinary_Z_order":order,"mean_y_node":actual,
                                 "V_node":arrays["V"][order],"m_node":arrays["m"][order],"passed":True})
                for quantity,values in arrays.items():
                    node = values[order]; source = self.nodes[node]
                    zero = source.get("operation")=="exact_rational" and source.get("numerator")==0
                    valid = zero or (source.get("operation")=="product"
                                     and source.get("arguments",[]).count(69245)==1)
                    if not valid:
                        raise ValueError("Independent normalized leading port changed")
                    normalization.append({"quantity":quantity,"ordinary_Z_order":order,
                                          "node":node,"inverse_Pstar_once":True,"passed":True})
            chart = loop["charts"][name]
            charts[name] = {"coordinate":chart["coordinate"],"radius":chart["radius"],
                            "native_J":chart["physical_measure"],"phase":chart["actual_phase"],
                            "independent_normalized_inputs":arrays,"mean_ODE_rows":ode_rows,
                            "normalization_rows":normalization,"passed":True}
        seam = self.parent["source_derivative_theorem"]["seams"]["macro_first"]
        leading = seam["leading_Qy_inputs"]
        if not leading["passed"] or len(leading["rows"])!=8:
            raise ValueError("Passed independent leading seam inputs required")
        for row in leading["rows"]:
            quantity = {"V_y":"Vy"}.get(row["quantity"],row["quantity"])
            order = row["ordinary_Z_order"]
            if (row["left_node"]!=charts["frozen_macro"]["independent_normalized_inputs"][quantity][order]
                    or row["right_node"]!=charts["first_switch"]["independent_normalized_inputs"][quantity][order]
                    or row["actual_function_identity"].get("namespace")!="leading"
                    or row["actual_function_identity"].get("passed") is not True):
                raise ValueError("Foreign modified port in independent leading seam")
        geometry = seam["radius_phase_identity"]
        if not geometry["passed"] or not geometry["common_radius_source_identity_proved"]:
            raise ValueError("Common physical-radius source identity required")
        pressure = self.parent["source_seam_proofs"]["physical_pressure_formula"]
        name = PREFIX+"compliant_current_original_pulse_radial_pressure_C3_functions.py"
        path = self.directory/name; raw = path.read_bytes(); text = raw.decode("utf-8")
        if self.input_hashes.get(name)!=hashlib.sha256(raw).hexdigest():
            raise ValueError("Pinned physical density/FTC source required")
        module = ast.parse(text); functions = {}
        for name,patterns in (("pressure_density",(
                "Py=a.add(a.scale(S,g.mul(g.constant('-1/2'),R)),a.trim(fields['Q_y']),a.trim(fields['Q_yy']))",
                "N1p=a.scale(Py,g.quotient(g.one,R,'positive actual pulse radius'))",
                "m15_y=a.scale(Py,g.neg(R))")),
                ("pressure_FTC",("R_q_y=a.scale(a.add(Qy,Qyy),R)",))):
            function = next(item for item in module.body if isinstance(item,ast.FunctionDef) and item.name==name)
            body = ast.get_source_segment(text,function)
            if not all(pattern in body for pattern in patterns):
                raise ValueError("Actual physical viscous source formula changed")
            functions[name] = {"function_AST_sha256":hashlib.sha256(ast.dump(function,include_attributes=False).encode()).hexdigest(),
                               "source_executed":False,"passed":True}
        endpoints = {}
        for label,chart,query in (("macro0","frozen_macro",0),("macro1","frozen_macro",1),
                                  ("first0","first_switch",0),("first1","first_switch",1)):
            spec = charts[chart]
            endpoints[label] = {"operation":"local_independent_leading_RQy_source_recipe",
                                "chart":chart,"coordinate_node":spec["coordinate"],
                                "native_query_node":query,"native_query_value":str(query),
                                "physical_R_node":spec["radius"],"physical_Pstar_node":pstar,
                                "positive_L_node":25272,"delta_node":25220,"ordinary_Z_node":20,
                                "normalized_inputs":spec["independent_normalized_inputs"],
                                "Qy_formula":"[2Z Vy-(1-delta)Z(V-m)-(1-Z^2)(V_Z-m_Z)]/L",
                                "H_formula":"R*Pstar*Qy_normalized_leading",
                                "ordinary_Z_orders":[0,1],"recipe_is_archived_graph_node":False,
                                "source_value_evaluated":False,"modified_Qy_alias":False}
        return {"passed":True,"independent_leading_charts":charts,"endpoint_recipes":endpoints,
                "physical_density_FTC_source":path.relative_to(ROOT).as_posix(),
                "physical_density_FTC_sha256":hashlib.sha256(raw).hexdigest(),
                "physical_density_FTC_functions":functions,"pressure_formula":pressure,
                "macro1_first0_leading_input_identity_sha256":digest(leading),
                "common_radius_source_identity_sha256":digest(geometry),
                "local_source_identity":"Vmacro+Vfirst=(Hmacro0-Hmacro1)+(Hfirst0-Hfirst1)=Hmacro0-Hfirst1",
                "mean_ODE_partition":"Hlower-Hupper=Elower-Eupper+integral R*C(Uy) dy",
                "endpoint_operator":"E=2Z*R*Uy/L",
                "source_operator":"C(Uy)=[(1-delta)Z*Uy+(1-Z^2)*Uy_Z]/L",
                "physical_R_once":True,"physical_Pstar_once":True,"ordinary_Z1":True,
                "absolute_mean_retained_in_H":True,"first1_H_not_zeroed":True,
                "new_graph_nodes_registered":False,"numeric_endpoint_or_integral_values_evaluated":False,
                "old_error_caps_cancelled_or_removed":False}

    def build(self):
        from .original_macro_leading_source import build_macro_leading_source
        from .original_switch_leading_source import build_switch_leading_source
        arguments = dict(nodes=self.nodes,loop_archive=self.loop_archive,
                         modified_archive=self.modified_archive,high_archive=self.high_archive,
                         parent_report=self.parent,directory=self.directory,root=ROOT)
        macro = build_macro_leading_source(**arguments)
        switch = build_switch_leading_source(**arguments)
        if switch.get("integral_C0_interface_passed") is not True:
            raise ValueError("Advertised phase-one signed source integral C0 interface failed")
        for packet in (macro,switch):
            if (packet.get("source_bindings_passed") is not True
                    or packet.get("endpoint_interface_passed") is not True
                    or not all(row.get("passed") is True for row in packet["family_checks"].values())):
                raise ValueError("Advertised source binding/endpoint arithmetic interface failed")
            self.bind_pins(packet)
        for name in ("original_leading_viscous_m15.py","original_leading_viscous_kernel.py",
                     "original_macro_leading_source.py","original_switch_leading_source.py"):
            relative = "src/openai_ns_reconstruction/"+name
            self.implementation_hashes[relative] = hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        return {"source_family":self.family,"parent_report":PARENT,"parent_report_sha256":PARENT_SHA256,
                "source_prefix_checks":self.prefix_checks,"macro_source":macro,"switch_source":switch,
                "independent_leading_FTC":self.leading_ftc(),
                "inherited_m15_packet":self.preservation(),
                "accepted_interface_scope":"signed endpoint arithmetic, exact-endpoint-mass source integral C0 and independent leading source-functional partition",
                "numeric_integral_scope_complete":False,
                "contract":{"ancestor_constructors_producers_calculations_not_run":True,
                            "point_approximations_not_exact_source_values":True,
                            "positive_remainders_not_signed_centers":True,
                            "raw213_reduced183_and_target_retained":True},
                "limits":{"original_full_signed_source_values":False,"numeric_Rp_closure":False,
                          "global_all_five_common_N_closure":False,"exact_heat_energy_stress_flat":False,
                          "true_temporal_n_recursion_summation":False,"mean_oscillatory_correction":False,
                          "independent_corrected_Cartesian_NS":False,"full_goal_complete":False}}

"""New source-leaf jets and bounded pulse mass over an unchanged m15 ledger.

Accepted reports are hydrated and hash checked.  No accepted constructor,
producer, checker or fixture is invoked by this owner.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_leading_viscous_kernel import (
    decode_pair, encode_pair, endpoint_coefficient, mean_source_coefficient, mean_source_value,
)

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT/"experiments/root_st073"
PREFIX = "lei_ren_part1_paper_"
PARENT = PREFIX+"compliant_current_original_leading_viscous_source_numeric.json"
PARENT_SHA256 = "d1aa3644ea8158fe78dccb210a0a75a4c81ad257b8504758cb5797845b8003d4"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()


class OriginalSourceJetMass:
    def __init__(self, *, directory=DIRECTORY):
        self.directory = Path(directory)
        self.parent_path = self.directory/PARENT
        if sha(self.parent_path)!=PARENT_SHA256:
            raise ValueError("Pinned accepted leading source report required")
        self.parent = json.loads(self.parent_path.read_bytes())
        if self.parent["checks"]["passed"] is not True:
            raise ValueError("Passed accepted leading source report required")
        self.input_hashes = dict(self.parent["input_hashes"])
        self.input_hashes[PARENT] = PARENT_SHA256
        self.implementation_hashes = dict(self.parent["implementation_hashes"])
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:
                    raise ValueError("Changed inherited pin: "+name)

    def archive(self,suffix):
        name = PREFIX+"compliant_current_original_"+suffix+".json.gz"
        if name not in self.input_hashes:
            raise ValueError("Unpinned source archive")
        with gzip.open(self.directory/name,"rt",encoding="utf-8") as stream:
            return json.load(stream)

    def bind_pins(self,packet):
        """Collect nested source pins, including a source AST attachment."""
        count = 0
        def visit(value):
            nonlocal count
            if isinstance(value,dict):
                if "path" in value and "sha256" in value:
                    relative = value["path"]
                    if not isinstance(relative,str) or not relative:
                        raise ValueError("New source pin path required")
                    path = ROOT/relative
                    actual = sha(path)
                    if value.get("passed") is not True or actual!=value["sha256"]:
                        raise ValueError("New source pin failed: "+relative)
                    group,key = ((self.input_hashes,path.name) if relative.startswith("experiments/root_st073/")
                                 else (self.implementation_hashes,relative))
                    if key in group and group[key]!=actual:
                        raise ValueError("New source pin conflicts with accepted source")
                    group[key] = actual; count += 1
                for item in value.values():
                    visit(item)
            elif isinstance(value,(list,tuple)):
                for item in value:
                    visit(item)
        visit(packet.get("source_pins",{}))
        return count

    def interior_C0(self,mass):
        """Evaluate the signed raw46 main component; keep raw47 separately."""
        c = mp.mp.clone(); c.dps = 120
        iv = MPIntervalContext(); iv.dps = 160
        delta = c.exp(-4*(c.exp(40)+11)-30)
        di = iv.exp(-4*(iv.exp(40)+11)-30)
        mass_index = {c.mpf(row["phase"]):row for row in mass["rows"]}
        source_index = {}
        integral_index = {}
        for row in self.parent["switch_source"]["point_remainder_rows"]:
            coordinate = row.get("coordinate",{}).get("value")
            if coordinate is None:
                continue
            key = (row["Z"],c.mpf(coordinate))
            if row.get("name")=="first_switch_viscous_lower":
                source_index[key] = row
            elif row.get("name")=="first_switch_viscous_source_integral":
                integral_index[key] = row
        rows = []; source_comparisons = []; endpoint_comparisons = []
        for key,original in integral_index.items():
            ztext,phase = key
            source = source_index[key]
            inputs = [decode_pair(c,iv,row) for row in source["source_point"]]
            bp = [row[0] for row in inputs]; bi = [row[1] for row in inputs]
            cp = mean_source_value(c.mpf(ztext),delta,bp)
            ci = mean_source_value(iv.mpf(ztext),di,bi)
            m = mass_index[phase]; mpnt,mi = decode_pair(c,iv,m)
            point = cp*mpnt; truth = ci*mi
            if "C0_recomputed_without_mass" in original:
                old = c.mpf(original["C0_recomputed_without_mass"]["approximation"])
                difference = abs(cp-old)/max(c.mpf(1),abs(cp),abs(old))
                comparison = {"Z":ztext,"phase":str(phase),"relative_difference":c.nstr(difference,30),
                              "fixed_point_tolerance":"1e-70","passed":difference<c.mpf("1e-70")}
                if not comparison["passed"]:
                    raise ValueError("Independent saved C(B) C0 comparison failed")
                source_comparisons.append(comparison)
            old = c.mpf(original["independent_saved_row"][0]["approximation"])
            relative = abs(point-old)/max(c.mpf(1),abs(point),abs(old))
            comparison = {"relative_point_difference":c.nstr(relative,35),
                          "point_difference_is_not_functional_equality":True,
                          "source_error_radius_used_as_tolerance":False,
                          "interval_overlap_used_as_equality":False}
            if phase in (0,1):
                comparison.update(scope="exact_mass_endpoint_main_component_only",fixed_point_tolerance="1e-70",
                                  passed=relative<c.mpf("1e-70"))
                if not comparison["passed"]:
                    raise ValueError("Independent saved raw46 endpoint C0 comparison failed")
                endpoint_comparisons.append(comparison)
            else:
                comparison.update(scope="different_numerical_mass_approximations_observed",passed=None)
            rows.append({"Z":ztext,"phase":c.nstr(phase,10),"raw_number":46,
                         "source_binding":"same saved first_switch_lower_source_Uy ordinary B/B_Z",
                         "mass_binding":"new original gp-ratio integral with directed rectangle certification",
                         "mean_source_C0":encode_pair(c,iv,cp,ci),
                         "main_component_C0":encode_pair(c,iv,point,truth),
                         "physical_scale":original["physical_scale"],
                         "independent_saved_C0_comparison":comparison,
                         "C1_computed":False,"raw47_weighted_remainder_retained":True,
                         "main_component_is_not_full_weighted_source_integral":True})
        return {"passed":bool(rows) and all(row["passed"] for row in source_comparisons+endpoint_comparisons),
                "row_count":len(rows),"rows":rows,"independent_C_without_mass_comparisons":source_comparisons,
                "independent_exact_mass_endpoint_comparisons":endpoint_comparisons,
                "C0_main_component_evaluated":True,"C1_evaluated":False,
                "signed_full_primitive_evaluated":False,"old_error_caps_removed":False}

    def lower_operators(self,lower):
        c = mp.mp.clone(); c.dps = 120
        iv = MPIntervalContext(); iv.dps = 180
        delta = c.exp(-4*(c.exp(40)+11)-30)
        di = iv.exp(-4*(iv.exp(40)+11)-30)
        rows = []
        for row in lower["rows"]:
            inputs = [decode_pair(c,iv,value) for value in row["result"]]
            point = [value[0] for value in inputs]; truth = [value[1] for value in inputs]
            z = row["Z"]
            ep = endpoint_coefficient(c.mpf(z),delta,point)
            ei = endpoint_coefficient(iv.mpf(z),di,truth)
            cp = mean_source_coefficient(c.mpf(z),delta,point)
            ci = mean_source_coefficient(iv.mpf(z),di,truth)
            rows.append({"Z":z,"coordinate":"second_micro native2; macro lower endpoint",
                         "signed_endpoint_E_coefficient_C0_C1":[encode_pair(c,iv,p,t) for p,t in zip(ep,ei)],
                         "local_mean_source_C0_C1":[encode_pair(c,iv,p,t) for p,t in zip(cp,ci)],
                         "physical_scale":{"reference_radius":"4*epsilon","epsilon_once_outside_R":True,
                                           "width":"hbB","formula":"4*epsilon^2*hbB times endpoint coefficient",
                                           "actual_radius_discrepancy_and_raw24_retained":True},
                         "macro_J_applied":False,"ordinary_Z1":True,"B_ZZ_input_required_and_present":True,
                         "signed_full_macro_viscous_integral_evaluated":False})
        return {"passed":len(rows)==3,"rows":rows,"raw24_positive_bound_replaced_or_removed":False,
                "source_integral_over_macro_evaluated":False}

    def build(self):
        from .original_interior_mass_source import build_interior_mass_source
        from .original_lower_dy_source import build_lower_dy_source
        mass = build_interior_mass_source(parent_report=self.parent_path,directory=self.directory,root=ROOT)
        if mass.get("passed") is not True:
            raise ValueError("New bounded original mass interface failed")
        self.bind_pins(mass)
        c0 = self.interior_C0(mass)
        high = self.archive("prescribed_high_radial_source_functions")
        lower = build_lower_dy_source(nodes=high["exact_graph_nodes"],high_archive=high,
                                     parent_report=self.parent,directory=self.directory,root=ROOT)
        if lower.get("passed") is not True or len(lower.get("rows",()))!=3:
            raise ValueError("New transported lower source interface failed")
        self.bind_pins(lower)
        lower_operators = self.lower_operators(lower)
        inherited = self.parent["inherited_m15_packet"]
        return {"source_family":self.parent["source_family"],"parent_report":PARENT,"parent_report_sha256":PARENT_SHA256,
                "interior_mass":mass,"signed_source_integral_C0":c0,
                "transported_lower_dy":lower,"lower_local_operators":lower_operators,
                "inherited_m15_packet":inherited,"inherited_m15_packet_digest":digest(inherited),
                "scope":{"new_bounded_source_arithmetic":True,"interior_mass_evaluated":True,
                         "raw46_C0_main_component_evaluated":True,"transported_lower_dy_Z0_Z1_Z2_evaluated":True,
                         "lower_local_endpoint_and_mean_operators_evaluated":True,"ledger_replaced_or_removed":False},
                "contract":{"accepted_ancestor_replayed":False,"positive_caps_selected_as_values":False,
                            "interval_midpoints_claimed_as_actual_values":False,
                            "physical_R_epsilon_and_distinct_widths_retained":True,
                            "raw213_reduced183_and_independent_target_retained":True},
                "limits":dict(self.parent["limits"],source_integral_C1_complete=False,
                              signed_weighted_primitive_remainder=False,full_goal_complete=False)}

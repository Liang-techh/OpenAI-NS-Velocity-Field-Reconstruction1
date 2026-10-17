"""Signed source-integral C1 over hydrated original mass and source jets.

The weighted primitive remainder and the complete m15 ledger are retained.
Spatial source derivatives do not implement temporal coefficient recursion.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import fzero, mpf_abs, mpf_cmp

from .original_leading_viscous_kernel import decode_pair, encode_pair, mean_source_coefficient
from .original_source_jet_kernel import ordinary_product, ordinary_quotient

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT/"experiments/root_st073"
PARENT = "lei_ren_part1_paper_compliant_current_original_source_jet_mass_numeric.json"
PARENT_SHA256 = "c3b2e2a974ebca402e2d2762ab2c508150fdd72f80aa2882152e944bdad84efe"


def correlated_phi_ratio(ctx,iv,phi,continuation_caps,hcap,log_h,*,shared_phi_identity_proved,
                         cap_linear_width_proved,error_threshold="1e-1000"):
    """Enclose Phi0/(Phi0+Delta) using the exact shared-Phi identity.

    Each saved ordinary continuation cap is hcap times a positive bound
    for Delta/h, as certified by its source formula.  First form G from
    Delta/h and the actual denominator envelope.  Combining log(h)+log|G|
    before selecting a representable error cap retains the defining width.
    No width or source value is selected from its enclosure.
    """
    if len(phi)!=3 or len(continuation_caps)!=3:
        raise ValueError("Shared Phi0 and ordinary native2 continuation rows through Z2 required")
    if shared_phi_identity_proved is not True or cap_linear_width_proved is not True:
        raise ValueError("Same-source shared-Phi and linear-width continuation proofs are required")
    def magnitude(value):
        a,b=(mpf_abs(row) for row in value._mpi_)
        return ctx.make_mpf(a if mpf_cmp(a,b)>=0 else b)
    def symmetric(value):return iv.mpf([-1,1])*iv.mpf(magnitude(value))
    # Positivity of hcap is checked by the same directed source quotient.
    cap_jet=[hcap,iv.mpf(0),iv.mpf(0)]
    delta=[symmetric(value) for value in continuation_caps]
    delta_per_h=ordinary_quotient(delta,cap_jet)
    denominator=[a+b for a,b in zip(phi,delta)]
    coefficient=ordinary_quotient(delta_per_h,denominator)
    threshold=iv.mpf(error_threshold)
    log_threshold=iv.log(threshold)
    proof=[];truth=[]
    for order,value in enumerate(coefficient):
        bound=magnitude(value)
        if bound==0:
            error=iv.mpf(0);log_bound=None
        else:
            log_bound=log_h+iv.log(iv.mpf(bound))
            if mpf_cmp(log_bound._mpi_[1],log_threshold._mpi_[0])>=0:
                raise ValueError("The original log-width does not admit the requested ratio-error cap")
            error=iv.mpf([-1,1])*threshold
        truth.append(iv.mpf(1 if order==0 else 0)+error)
        proof.append({"ordinary_Z_order":order,
                      "coefficient_absolute_upper_exact_mpf_tuple":list(bound._mpf_),
                      "combined_original_width_log_bound_exact_mpf_tuples":None if log_bound is None else [list(v) for v in log_bound._mpi_],
                      "representable_absolute_error_cap":error_threshold,
                      "source_width_materialized_or_selected":False,"passed":True})
    return [ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)],truth,{
        "identity":"Phi0/(Phi0+Delta)=1-h*(Delta/h)/(Phi0+Delta)",
        "shared_Phi0_cancellation_is_algebraic":True,
        "continuation_caps_used_only_as_bounds":True,
        "saved_cap_linear_hcap_factor_required":True,
        "original_width_log_combined_before_error_cap":True,
        "ratio_approximation_is_exact_source_value":False,"orders":proof,"passed":True}


def pressure_drive_ordinary(z,delta,p0):
    """Ordinary Z0..2 of the prescribed pressure direction, using p0 Z0..3."""
    if len(p0)!=4 or any(value is None for value in p0):
        raise ValueError("Pressure source ordinary Z2 requires actual p0 through ordinary Z3")
    p,pz,pzz,pzzz=p0
    d=1-z*z
    numerator=[2*(1+delta)*z*p-d*pz,
               2*(1+delta)*p+(4+2*delta)*z*pz-d*pzz,
               (6+4*delta)*pz+(6+2*delta)*z*pzz-d*pzzz]
    denominator=[2*(1-delta*z*z),-4*delta*z,-4*delta]
    return ordinary_quotient(numerator,denominator)


def pressure_source_ordinary(z,delta,p0,quotient,physical_scale):
    """Before-width pressure source -R100*(Pstar^2/epsilon)*q*drive.

    R100, epsilon and Pstar are constant under the ordinary Z derivative.
    The macro angular exp(ell)-1 remainder must be added by the caller.
    """
    if len(quotient)!=3:
        raise ValueError("The original macro quotient needs ordinary Z0..2")
    drive=pressure_drive_ordinary(z,delta,p0)
    return [-physical_scale*value for value in ordinary_product(quotient,drive)]


def macro_nonlinear_enclosure(ctx,iv,source,weight,hcap,*,macro_norm_bound_proved):
    """Retain ordinary order-two coupling in B*(exp(ell)-1)."""
    if len(source)!=3 or macro_norm_bound_proved is not True:
        raise ValueError("Three ordinary source rows and the macro ell-norm source proof are required")
    if any(mpf_cmp(value._mpi_[0],fzero)<=0 for value in (weight,hcap)):
        raise ValueError("Positive source weight and macro cap enclosures are required")
    norm=iv.mpf(0)
    for order,row in enumerate(source):
        a,b=(mpf_abs(v) for v in row._mpi_)
        upper=ctx.make_mpf(a if mpf_cmp(a,b)>=0 else b)
        norm+=iv.mpf(upper)*weight**order/math.factorial(order)
    remainder=norm*6*hcap*iv.exp(6*hcap)
    errors=[remainder*math.factorial(order)/weight**order for order in range(3)]
    truth=[row+iv.mpf([-1,1])*error for row,error in zip(source,errors)]
    return truth,{
        "formula":"norm_order2(B)*6*cap*exp(6*cap); ordinary error k!/weight^k",
        "source_norm_exact_mpf_tuples":[list(v) for v in norm._mpi_],
        "ordinary_errors_exact_mpf_tuples":[[list(v) for v in error._mpi_] for error in errors],
        "macro_ell_source_bound_used":True,"first_switch_ell_not_substituted":True,
        "point_unchanged":True,"passed":True}


def main_coefficient(z,delta,source,mass,*,mass_Z_independent):
    """Ordinary Z0/Z1 of C(B)*M at a fixed original pulse phase."""
    if mass_Z_independent is not True:
        raise ValueError("The original source mass must be proved Z-independent")
    if len(source)!=3 or any(value is None for value in source):
        raise ValueError("Actual ordinary B/B_Z/B_ZZ inputs are required")
    return [value*mass for value in mean_source_coefficient(z,delta,source)]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()


class OriginalSignedC1Source:
    def __init__(self,*,directory=DIRECTORY):
        self.directory = Path(directory)
        self.parent_path = self.directory/PARENT
        if sha(self.parent_path)!=PARENT_SHA256:
            raise ValueError("Pinned accepted source-jet/mass report required")
        self.parent = json.loads(self.parent_path.read_bytes())
        if self.parent["checks"]["passed"] is not True:
            raise ValueError("Passed accepted parent report required")
        self.input_hashes = dict(self.parent["input_hashes"])
        self.input_hashes[PARENT] = PARENT_SHA256
        self.implementation_hashes = dict(self.parent["implementation_hashes"])
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for relative,expected in group.items():
                if sha(base/relative)!=expected:
                    raise ValueError("Changed accepted pin: "+relative)

    def bind_pins(self,packet):
        def visit(value):
            if isinstance(value,dict):
                if "path" in value and "sha256" in value:
                    relative = value["path"]; path = ROOT/relative
                    actual = sha(path)
                    if value.get("passed") is not True or value["sha256"]!=actual:
                        raise ValueError("New source pin failed: "+relative)
                    group,key = ((self.input_hashes,path.name) if relative.startswith("experiments/root_st073/")
                                 else (self.implementation_hashes,relative))
                    if group.get(key,actual)!=actual:
                        raise ValueError("New source pin conflicts with accepted source")
                    group[key] = actual
                for item in value.values():
                    visit(item)
            elif isinstance(value,(list,tuple)):
                for item in value:
                    visit(item)
        visit(packet.get("source_pins",{}))

    def integral_rows(self,source):
        c = mp.mp.clone(); c.dps = 120
        iv = MPIntervalContext(); iv.dps = 180
        delta = c.exp(-4*(c.exp(40)+11)-30)
        di = iv.exp(-4*(iv.exp(40)+11)-30)
        mass = {c.mpf(row["phase"]):row for row in self.parent["interior_mass"]["rows"]}
        source_index = {c.mpf(row["Z"]):row for row in source["rows"]}
        rows = []
        for inherited in self.parent["signed_source_integral_C0"]["rows"]:
            ztext = inherited["Z"]; phase = inherited["phase"]
            source_row = source_index[c.mpf(ztext)]
            inputs = [decode_pair(c,iv,row) for row in source_row["source_rows"]]
            point = [value[0] for value in inputs]; truth = [value[1] for value in inputs]
            mpnt,mi = decode_pair(c,iv,mass[c.mpf(phase)])
            cp = main_coefficient(c.mpf(ztext),delta,point,mpnt,mass_Z_independent=True)
            ci = main_coefficient(iv.mpf(ztext),di,truth,mi,mass_Z_independent=True)
            rows.append({"Z":ztext,"phase":phase,"ordinary_Z_orders":[0,1],"raw_number":46,
                         "source_scope":source["source_scope"],
                         "main_component_C0_C1":[encode_pair(c,iv,p,t) for p,t in zip(cp,ci)],
                         "physical_scale":inherited["physical_scale"],
                         "mass_Z_independent":True,"source_B_ZZ_required_and_present":True,
                         "weighted_raw47_remainder_retained":True,"full_weighted_primitive_evaluated":False,
                         "source_equality_from_interval_overlap":False})
        return {"passed":len(rows)==19,"rows":rows,"row_count":len(rows),
                "ordinary_Z1_retained":True,"physical_R_epsilon_and_hbS_squared_retained":True,
                "full_source_scope":source["full_source_scope"],"full_weighted_primitive_evaluated":False}

    def build(self):
        from .original_modal_source_point import build_modal_source_point
        from .original_factored_source_jet import build_factored_source_jet
        source = build_modal_source_point(parent_report=self.parent_path,directory=self.directory,root=ROOT)
        if source.get("passed") is not True or len(source.get("rows",()))!=3:
            raise ValueError("New same-generator numerical source interface is not ready")
        if source.get("source_family")!=self.parent["source_family"]:
            raise ValueError("Different original source family")
        self.bind_pins(source)
        rejected = build_factored_source_jet(nodes=None,high_archive={"source_family":self.parent["source_family"]},
                                             parent_report=self.parent,directory=self.directory,root=ROOT)
        if rejected.get("passed") is not False or rejected.get("rows"):
            raise ValueError("The rejected R100/macro source alias must stay blocked")
        self.bind_pins(rejected)
        integrals = self.integral_rows(source)
        return {"source_family":self.parent["source_family"],"parent_report":PARENT,"parent_report_sha256":PARENT_SHA256,
                "source_jet":source,"signed_source_integral":integrals,"rejected_draft_source_origin":rejected,
                "inherited_m15_packet":self.parent["inherited_m15_packet"],
                "inherited_m15_packet_digest":self.parent["inherited_m15_packet_digest"],
                "inherited_mass_packet_digest":digest(self.parent["interior_mass"]),
                "scope":{"source_main_component_C0_C1_evaluated":True,"source_scope":source["source_scope"],
                         "full_source_B_ZZ_numeric":source["full_source_scope"],"ledger_replaced_or_removed":False},
                "contract":{"accepted_ancestor_replayed":False,"source_caps_selected_as_values":False,
                            "midpoint_as_function_value":False,"saved_mass_recomputed":False,
                            "P1_m13_used_as_m15":False,"raw213_reduced183_and_target_retained":True},
                "limits":dict(self.parent["limits"],numeric_integral_scope_complete=False,
                              source_primitive_weighted_remainder_complete=False,full_goal_complete=False)}

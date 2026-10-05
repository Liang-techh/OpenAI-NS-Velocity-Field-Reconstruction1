"""Current physical source ownership, scale ledgers and Cartesian/time scope."""
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_downstream_physical_assembly import (
    CurrentDownstreamPhysicalAssembly, CHARTS, OPEN, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_global_physical_assembly import INDICES, COMPONENTS, MICRO, UZ, UT, UR, P, cartesian_templates, physical_operators
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"current_downstream_physical_assembly.json"


def patch_fixed_unit_fixture():
    """Differentiate full radial fields first; freeze the basepoint unit."""
    from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110
        y0=mp.mpf(".31");z0=mp.mpf(".23");rm=mp.mpf(2);tol=mp.mpf("1e-65")
        x0=mp.exp(y0);r0=rm*x0
        field=object.__new__(CompliantGlobalPhysicalAssembly)
        field.ctx=c;field.logP=c.ln(3)
        fields=(lambda y,z:mp.sqrt(rm/2)*mp.exp(y/2),
            lambda y,z:mp.sqrt(rm/2)*mp.exp(y/2)*
                (mp.exp(-mp.mpf(".13")*y)*(1+z+z**4)+mp.mpf(".2")*y**3*z**2))
        count=0;first=None
        for full in fields:
            raw={};expected={}
            for k in range(5):
                for n in range(5-k):
                    key="y"+str(k)+"_Z"+str(n)
                    actual=mp.diff(full,(y0,z0),(k,n))
                    normalized=actual/mp.sqrt(rm/2)
                    raw[key]=c.mpf([normalized-tol,normalized+tol])
                    expected[k,n]=actual/mp.sqrt(r0/2)
            packet={"physical_velocity_pressure_y_Z_mixed4":{
                "Ur_over_sqrt_Rm_over_2":raw,
                "Uz":raw,"Utheta_over_Pstar":raw,"P_over_Pstar2":raw}}
            grids,logs,amplitudes=field.normalized_sources(
                "actual_patch",c.mpf(z0),c.mpf(x0),packet,None,c.ln(c.mpf(r0)))
            if amplitudes[UR]!={5:mp.mpf(".5"),6:mp.mpf("-.5")}:
                raise ValueError("Physical radial frozen unit exponents changed")
            for index,target in expected.items():
                terms=grids[UR][index]
                if len(terms)!=1 or any(terms[0][0]):
                    raise ValueError("Ordinary patch source acquired unexpected factors")
                lo,hi=endpoints(terms[0][1])
                if not lo<=target<=hi:
                    raise ArithmeticError("Patch full-derivative/frozen-current-unit conversion failed")
                restored=terms[0][1]*c.exp((logs[5]-logs[6])/2)
                actual=mp.diff(full,(y0,z0),index)
                if not endpoints(restored)[0]<=actual<=endpoints(restored)[1]:
                    raise ArithmeticError("Patch radial physical unit restoration failed")
                count+=1
            if first is None:
                first=grids[UR][1,0][0][1]
                if not endpoints(first)[0]<=mp.mpf(".5")<=endpoints(first)[1]:
                    raise ArithmeticError("Constant Q full radial derivative must retain one half")
        return dict(independently_differentiated_full_patch_mixed4_rows=count,
            fixed_Rm_to_fixed_current_R_unit_conversion_checked=True,
            physical_radial_amplitude_restoration_checked=True,
            constant_shape_full_radial_first_derivative_nonzero=True,
            finite_fixture_only=True,passed=True)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentDownstreamPhysicalAssembly(require_checked=False)
        expected=encode(pack(field.manifest()))
        if {k:v for k,v in raw.items() if k!="whole_current_downstream_physical_maps"}!=expected:
            raise ValueError("Current physical source/owner/parameter/operator proof changed")
        if (not all(raw["current_provider_graph_identity"].values())
                or not all(raw["unchanged_coordinate_source_methods"].values())
                or not raw["exact_source_and_divergence_identity"]["passed"]
                or raw["current_downstream_cartesian_spatial4_time1_certified"]):
            raise ValueError("Current nested graph/operator/producer acceptance scope differs")
        registry=field.dispatch.manifest()["ordered_current_chart_registry"]
        c=MPIntervalContext();c.dps=240
        indices={"x"+str(i)+"_y"+str(j)+"_z"+str(k) for i,j,k in INDICES}
        native=json.loads((HERE/(PREFIX+"current_matched_source_dispatcher.json")).read_bytes())
        _verify_hashes(native)
        names={UZ:"Uz",UT:"Utheta_over_current_Utheta",UR:"Ur_over_current_sqrt_R_over_2",P:"P_over_Pstar2"}
        availability={}
        for chart in ("switch_first","switch_second"):
            ledger=native["whole_current_chart_evaluations"][chart]["source_packet"]["final_factored_physical_row_ledgers"]
            available={}
            for label,name in names.items():
                keys={}
                for row in ledger:
                    if row["physical_row"].startswith(name+"/"):
                        key=row["physical_row"][len(name)+1:]
                        k,n=(int(part[1:]) for part in key.split("_"))
                        keys[k,n]=bool(row["terms"])
                if len(keys)!=15:raise ValueError("Current uncapped micro source ledger incomplete")
                available[label]=keys
            availability[chart]=available
        def expected_parts(chart,component,index):
            i,j,b=(int(part[1:]) for part in index.split("_"))
            labels=set()
            for label,a,q in cartesian_templates()[component,i,j,b]:
                if chart not in MICRO or any(availability[chart][label][k,n]
                        for k,n in physical_operators()[a,b]):
                    labels.add(label)
            return labels
        patch_fixture=patch_fixed_unit_fixture()
        counts={};terms=0;zeros=0;structural_zeros=0
        def check(row):
            nonlocal terms,zeros
            if bool(row["exact_zero"])!=(row["log_absolute_upper"] is None):
                raise ValueError("Exact zero row conflicts with logarithmic bound")
            if row["exact_zero"]:
                if row["terms"]:raise ValueError("Exact zero row has nonzero source terms")
                zeros+=1;return
            lo,hi=endpoints(read_interval(c,row["log_absolute_upper"]))
            if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                raise ArithmeticError("Nonfinite current Cartesian/time source bound")
            if not row["positive_source_exponentials_not_materialized"]:
                raise ValueError("Huge source exponentials must remain formal")
            for term in row["terms"]:
                if len(term["source_log_exponents"])!=7:
                    raise ValueError("Original source base count differs")
                lo,hi=endpoints(read_interval(c,term["signed_coefficient"]))
                if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                    raise ArithmeticError("Nonfinite retained current signed coefficient")
                terms+=1
        packets=raw["whole_current_downstream_physical_maps"]
        if tuple(packets)!=CHARTS:raise ValueError("All fourteen current owners required in original order")
        for chart,packet in packets.items():
            if (packet["current_source_owner"]!=registry[chart]["provider"]
                    or packet["current_source_acceptance_receipt"]!=registry[chart]["acceptance_receipt"]
                    or packet["datum_enclosure_sha256"]!=field.datum_sha
                    or not packet["current_source_dispatcher_used"]
                    or packet["current_downstream_cartesian_spatial4_time1_certified"]
                    or packet["current_Rp_external_pulse_join_certified"]):
                raise ValueError("Current physical chart owner/acceptance differs")
            spatial=packet["physical_spatial_cartesian_mixed4"];time=packet["first_fixed_x_physical_time_derivative"]
            if set(spatial)!=indices or set(time)!=set(COMPONENTS):
                raise ValueError("All spatial multiindices/fixed-x time components required")
            count=0
            for index,components in spatial.items():
                if set(components)!=set(COMPONENTS):raise ValueError("Cartesian component omitted")
                for label,parts in components.items():
                    expected_labels=expected_parts(chart,label,index)
                    if set(parts)!=expected_labels:
                        raise ValueError("Moving basis/current uncapped source contribution differs")
                    normal={UR,UT} if label in ("ux","uy") else {UZ} if label=="uz" else {P}
                    structural_zeros+=len(normal-expected_labels)
                    for row in parts.values():check(row);count+=1
            for label,parts in time.items():
                if set(parts)!=({UR,UT} if label in ("ux","uy") else {UZ} if label=="uz" else {P}):
                    raise ValueError("Fixed-x time component contribution omitted")
                for row in parts.values():check(row);count+=1
            expected_count=sum(len(expected_parts(chart,label,index)) for index in indices for label in COMPONENTS)+6
            if count!=expected_count:
                raise ValueError("Current source-supported spatial/time contribution count differs")
            if (not packet["moving_cylindrical_basis_differentiated"]
                    or not packet["normalization_not_differentiated_twice"]
                    or packet["microscope_phase_to_logR_applied_before_source_bound"]!=(chart in MICRO)
                    or packet["global_source_row_mode"]!=("uncapped_factored_rows" if chart in MICRO else "provider_prebounded_mixed_rows")
                    or any(packet[k] for k in SCOPES+OPEN)):
                raise ValueError("Original source factoring or unfinished physical/global scope changed")
            counts[chart]=count
        for chart in ("core","bridge_first","pulse_entrance","heat_exterior"):
            try:field.evaluate(chart,0,0)
            except ValueError:pass
            else:raise ValueError("Current downstream adapter silently accepted unadmitted whole-field owner")
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_input_hashes_checked=len(raw["input_hashes"]),
            current_provider_graph_identity=raw["current_provider_graph_identity"],
            current_physical_chart_owners_checked=14,current_physical_source_contributions_checked=counts,
            total_current_spatial_and_time_source_contributions_checked=sum(counts.values()),
            retained_signed_current_source_terms_checked=terms,exact_zero_source_contributions_checked=zeros,
            omitted_exact_zero_spatial_contributions_verified_from_current_micro_ledgers=structural_zeros,
            independent_current_patch_fixed_unit_fixture=patch_fixture,
            unchanged_independent_cartesian_coordinate_fixture=raw["reused_independent_coordinate_fixture"],
            unchanged_independent_microscopic_scale_fixture=raw["reused_independent_micro_scale_fixture"],
            current_downstream_cartesian_spatial4_time1_certified=True,
            current_source_divergence_identity_and_fixed_units_verified=True,
            moving_basis_and_fixed_x_time_source_rows_retained=True,
            current_Rp_external_pulse_join_certified=False,current_core_axis_physical_owner_installed=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current14 physical maps: Cartesian spatial4/fixed-x time1, current owners, retained signed source scales",flush=True)
    return result


if __name__=="__main__":
    run()

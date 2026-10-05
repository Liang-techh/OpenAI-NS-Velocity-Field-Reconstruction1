"""Current pulse physical source factors, coordinates and Cartesian/time rows."""
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_pulse_physical_assembly import (
    CurrentPulsePhysicalAssembly, BASE, PULSE_CHARTS, CHARTS, UNIFORM,
    SCOPES, OPEN, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    INDICES, COMPONENTS, UZ, UT, UR, P, cartesian_templates)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"current_pulse_physical_assembly.json"


def independent_pulse_unit_fixture():
    """Differentiate full finite fields first and freeze each native unit."""
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;tol=mp.mpf("1e-65")
        mu=mp.mpf(".07");lp=mp.log(2);rp=mp.log(3);z0=mp.mpf(".23")
        field=object.__new__(CurrentPulsePhysicalAssembly)
        field.ctx=c;field.params=SimpleNamespace(mu=c.mpf(mu))
        field.logP=c.ln(2);field.logRp=c.ln(3)
        raw_pressure=lambda y,z:1+z**4+(1+z*z)**-2*(1-mp.exp(-(1+2*mu)*(y-rp)))/(2*(1+2*mu))
        full={
            UT:lambda y,z:2*(1+z*z)**-1*mp.exp(-(.5+mu)*(y-rp)),
            UZ:lambda y,z:2*(1+z*z)**-1*mp.exp(-(.5+mu)*(y-rp))*
                (mp.sin(z)+mp.mpf(".2")*(y-rp)**3*z*z),
            UR:lambda y,z:mp.sqrt(mp.exp(y)/2)*2*(1+z*z)**-1*mp.exp(-(.5+mu)*(y-rp))*
                (1+z+mp.mpf(".1")*(y-rp)**2*z**3),
            P:lambda y,z:4*raw_pressure(y,z)}
        count=0;radius_count=0
        for chart,value,t0 in (("pulse_entrance",mp.mpf(".01"),mp.mpf(".01")),
                ("pulse_main",mp.mpf(".8"),mp.mpf(".8")/mu),
                ("pulse_gap_end",mp.mpf(-5),13/mu-5)):
            y0=rp+t0;raw={};targets={}
            units={UT:2*mp.exp(-(.5+mu)*t0),UZ:2*mp.exp(-(.5+mu)*t0),
                UR:2*mp.exp(-(.5+mu)*t0)*mp.sqrt(mp.exp(y0)/2),P:mp.mpf(4)}
            for label,fn in full.items():
                rows={}
                for k in range(5):
                    for n in range(5-k):
                        target=mp.diff(fn,(y0,z0),(k,n));v=target/units[label]
                        rows["y"+str(k)+"_Z"+str(n)]=c.mpf([v-tol,v+tol])
                        targets[label,k,n]=target
                raw[label]=rows
            packet={"physical_mixed_derivatives_total_order_le4":raw}
            logR,_=field.radius(chart,c.mpf(value),packet,None)
            with mp.workdps(160):
                exact_t=value if chart=="pulse_entrance" else value/mu if chart=="pulse_main" else 13/mu+value
                exact_y=mp.log(3)+exact_t
            if not endpoints(logR)[0]<=exact_y<=endpoints(logR)[1]:
                raise ArithmeticError("Independent pulse coordinate radius mismatch")
            radius_count+=1
            grids,logs,amplitudes=field.normalized_sources(chart,c.mpf(z0),c.mpf(value),packet,None,logR)
            for label in full:
                unit=c.exp(sum((logs[i]*power for i,power in amplitudes[label].items()),c.mpf(0)))
                for (k,n),terms in grids[label].items():
                    if len(terms)!=1 or any(terms[0][0]):
                        raise ValueError("Native ordinary source acquired derivative factors")
                    restored=terms[0][1]*unit;target=targets[label,k,n]
                    if not endpoints(restored)[0]<=target<=endpoints(restored)[1]:
                        raise ArithmeticError("Full pulse derivative/frozen-unit conversion failed")
                    count+=1
        # A finite log fixture with tiny mu exposes the lost shared-offset
        # correlation without materializing either gigantic exponential.
        field.params.mu=c.mpf("1e-1000");field.logP=c.mpf(".3");field.logRp=c.mpf(2)
        raw={label:{"y"+str(k)+"_Z"+str(n):c.mpf(1) for k in range(5) for n in range(5-k)}
            for label in (UT,UZ,UR,P)}
        packet={"physical_mixed_derivatives_total_order_le4":raw}
        logR,_=field.radius("pulse_gap",12,packet,None)
        _,oldlogs,_=BASE.normalized_sources(field,"pulse_gap",0,12,packet,None,logR)
        lost=oldlogs[4]+(oldlogs[5]-oldlogs[6])/2
        _,logs,amplitudes=field.normalized_sources("pulse_gap",0,12,packet,None,logR)
        reduced=sum((logs[i]*power for i,power in amplitudes[UR].items()),c.mpf(0))
        with mp.workdps(160):target=mp.mpf(".3")+1-12-mp.log(2)/2
        lo,hi=endpoints(reduced)
        if not lo<=target<=hi or hi-lo>mp.mpf("1e-80") or endpoints(lost)[1]<mp.mpf("1e800"):
            raise ArithmeticError("Giant pulse offset was not cancelled before bounds")
        return dict(independently_differentiated_full_pulse_mixed4_rows=count,
            independent_radial_coordinate_families_checked=radius_count,
            full_radial_prefactors_and_fixed_native_units_checked=True,
            absolute_pressure_Pstar_squared_units_checked=True,
            independent_giant_offset_log_fixture_checked=True,
            tiny_mu_giant_radius_correlation_preserved=True,
            finite_fixture_only=True,passed=True)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentPulsePhysicalAssembly(require_checked=False)
        expected=encode(pack(field.manifest()))
        extra=("whole_current_pulse_physical_maps","whole_current_gap_overlap_physical_map")
        if {k:v for k,v in raw.items() if k not in extra}!=expected:
            raise ValueError("Current pulse physical source, unit graph or operators changed")
        if (not all(raw["current_provider_graph_identity"].values())
                or not all(raw["original_radius_and_cartesian_time_operators_retained"].values())
                or not raw["pulse_fixed_unit_rebase_proof"]["passed"]
                or not raw["current_native_parameter_defining_source_bridge"]["passed"]
                or not raw["current_native_parameter_defining_source_bridge"]["exact_Md40_delta_choice_branch_proved"]
                or raw["pulse_fixed_unit_rebase_proof"]["exact_component_factor_identities_including_radial_orders0_through4"]!=20
                or raw["current_pulse_cartesian_spatial4_time1_certified"]
                or not raw["current_Rp_external_pulse_join_certified"]
                or raw[UNIFORM] or raw["full_pulse_C4_installed"]
                or tuple(raw["current_physical_chart_owners"])!=CHARTS):
            raise ValueError("Current physical graph/limited acceptance scope differs")
        registry=field.source_owners;c=MPIntervalContext();c.dps=240
        indices={"x"+str(i)+"_y"+str(j)+"_z"+str(k) for i,j,k in INDICES}
        packets=raw["whole_current_pulse_physical_maps"]
        if tuple(packets)!=PULSE_CHARTS:raise ValueError("Six new physical pulse maps required")
        cases={**packets,"pulse_gap_overlap":raw["whole_current_gap_overlap_physical_map"]}
        counts={};terms=0;zeros=0
        def check(row,label):
            nonlocal terms,zeros
            if bool(row["exact_zero"])!=(row["log_absolute_upper"] is None):
                raise ValueError("Source exact zero differs from bound")
            if row["exact_zero"]:
                if row["terms"]:raise ValueError("Exact zero has signed source terms")
                zeros+=1;return
            lo,hi=endpoints(read_interval(c,row["log_absolute_upper"]))
            if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                raise ArithmeticError("Nonfinite Cartesian/time pulse log bound")
            if (not row["positive_source_exponentials_not_materialized"]
                    or row["source_row_mode"]!="provider_prebounded_mixed_rows"
                    or endpoints(read_interval(c,row["physical_lambda_exponent"]))[1]>=0):
                raise ValueError("Original pulse source/lambda bound scope changed")
            for term in row["terms"]:
                powers=term["source_log_exponents"]
                if len(powers)!=7:raise ValueError("Full source base count differs")
                decoded=[]
                for value in powers:
                    scalar=mp.make_mpf(tuple(value["exact_mpf_tuple"]))
                    decoded.append((scalar,scalar))
                unit_ok=(decoded[2]==(1,1) and decoded[4]==(1,1) and decoded[5][1]<=0) if label!=P else (
                    decoded[1]==(2,2) and decoded[2]==(0,0) and decoded[4]==(0,0))
                if not unit_ok:
                    raise ValueError("Pulse radius/velocity factor rebase or pressure unit differs")
                lo,hi=endpoints(read_interval(c,term["signed_coefficient"]))
                if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                    raise ArithmeticError("Nonfinite signed physical source coefficient")
                terms+=1
        for chart,packet in cases.items():
            owner="pulse_gap" if chart=="pulse_gap_overlap" else chart
            if (packet["current_source_owner"]!=registry[owner]["provider"]
                    or packet["current_source_acceptance_receipt"]!=registry[owner]["acceptance_receipt"]
                    or packet["datum_enclosure_sha256"]!=field.datum_sha
                    or not packet["current_Rp_external_pulse_join_certified"]
                    or packet["current_pulse_cartesian_spatial4_time1_certified"]
                    or not packet["pulse_correlated_radius_velocity_factors_combined_before_bounds"]
                    or packet[UNIFORM] or packet["full_pulse_C4_installed"]):
                raise ValueError("Current physical pulse source/owner/scope differs")
            if chart=="pulse_gap_overlap" and (not packet["same_current_gap_owner_used"]
                    or packet["original_chart_for_physical_operators"]!="pulse_gap"):
                raise ValueError("Supplemental physical coverage must retain original gap owner")
            spatial=packet["physical_spatial_cartesian_mixed4"]
            time=packet["first_fixed_x_physical_time_derivative"];count=0
            if set(spatial)!=indices or set(time)!=set(COMPONENTS):
                raise ValueError("Physical spatial4/time1 scope omitted")
            for index,components in spatial.items():
                i,j,b=(int(v[1:]) for v in index.split("_"))
                if set(components)!=set(COMPONENTS):raise ValueError("Cartesian component omitted")
                for component,parts in components.items():
                    labels={label for label,a,q in cartesian_templates()[component,i,j,b]}
                    if set(parts)!=labels:raise ValueError("Moving-basis physical contribution omitted")
                    for label,row in parts.items():check(row,label);count+=1
            for component,parts in time.items():
                labels={UR,UT} if component in ("ux","uy") else {UZ} if component=="uz" else {P}
                if set(parts)!=labels:raise ValueError("Fixed-position time contribution omitted")
                for label,row in parts.items():check(row,label);count+=1
            if count!=216:raise ValueError("Pulse spatial/time source contribution count differs")
            if (not packet["moving_cylindrical_basis_differentiated"]
                    or not packet["normalization_not_differentiated_twice"]
                    or any(packet[k] for k in SCOPES+OPEN)):
                raise ValueError("Original operator or unfinished global scope differs")
            counts[chart]=count
        for chart in ("core","bridge_first","flatten","heat_exterior"):
            try:field.evaluate(chart,0,0)
            except ValueError:pass
            else:raise ValueError("Unadmitted physical whole-field owner accepted")
        try:field.evaluate("pulse_gap",0,"12.0001")
        except ValueError:pass
        else:raise ValueError("Private gap overlap widened public route")
        fixture=independent_pulse_unit_fixture()
        retained=raw["retained_fourteen_chart_physical_evidence"]
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_input_hashes_checked=len(raw["input_hashes"]),
            current_physical_chart_owners_checked=20,
            newly_mapped_native_pulse_chart_owners_checked=6,
            retained_fourteen_chart_physical_evidence=retained,
            new_current_pulse_spatial_and_time_source_contributions_checked=counts,
            total_new_source_contributions_including_overlap_checked=sum(counts.values()),
            twenty_owner_regular_source_contributions_checked=retained["total_source_contributions"]+sum(counts[k] for k in PULSE_CHARTS),
            retained_signed_pulse_source_terms_checked=terms,exact_zero_pulse_source_contributions_checked=zeros,
            exact_pulse_fixed_unit_factor_rebase_identities_checked=20,
            independent_pulse_fixed_unit_and_coordinate_fixture=fixture,
            unchanged_independent_cartesian_coordinate_fixture=raw["reused_independent_coordinate_fixture"],
            current_pulse_cartesian_spatial4_time1_certified=True,
            current_twenty_downstream_physical_source_ownership_certified=True,
            current_Rp_external_pulse_join_certified=True,
            **{UNIFORM:False},full_pulse_C4_installed=False,
            current_core_axis_physical_owner_installed=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current20 physical owners: six new pulse spatial4/time1 maps, exact factor rebase, independent unit fixture",flush=True)
    return result


if __name__=="__main__":
    run()

"""Current implicit patch to Rh reference: shared analytic datum and mixed4.

A companion adapter owns only the nine admitted downstream source charts.
The accepted eight-chart dispatcher and historical pre-pulse files are
unchanged. No pressure representative or production point is selected.
"""
import ast
import json
from pathlib import Path

import mpmath as mp
import sympy as s

import lei_ren_part1_paper_compliant_five_moment_repair as repair_module
import lei_ren_part1_paper_compliant_core_physical_field as core_module
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum, CompliantOuterParameters
from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum, BETA2, BETA0
from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 import CompliantPrePulseMixedC4
from lei_ren_part1_paper_compliant_current_matched_source_dispatcher import (
    CurrentMatchedSourceDispatcher, CURRENT_ROUTES, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import (
    class_assignment, keyword_binding, return_call_binding, expression_call_binding)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings, actual_R100_trace_bindings
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

RH_ROUTE=dict(provider=PREFIX+"pre_pulse_mixed_C4.CompliantPrePulseMixedC4",
    method="reference",coverage_coordinate="offset=log(R/Rref)",domain=[-5,0],
    acceptance_receipt=PREFIX+"actual_Rh_source_join_check.json")


def tree(full_stem):
    return ast.parse((HERE/(full_stem+".py")).read_text(encoding="utf8"))


def import_binding(full_stem, module, name, alias=None):
    matches=[v for n in tree(full_stem).body if isinstance(n,ast.ImportFrom)
        and n.module==module for v in n.names if v.name==name and v.asname==alias]
    if len(matches)!=1:raise ValueError("Canonical pressure import changed: "+full_stem)
    return True


def external_assignment(full_stem, class_name, method, expected):
    cls=next(n for n in tree(full_stem).body if isinstance(n,ast.ClassDef) and n.name==class_name)
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==method)
    out={}
    for target,expression in expected.items():
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(v)==target for v in n.targets)]
        want=ast.dump(ast.parse(expression,mode="eval").body)
        if len(values)!=1 or ast.dump(values[0])!=want:
            raise ValueError("Canonical pressure source changed: "+full_stem+"."+target)
        out[target]=True
    return out


def pressure_function_source_bindings():
    imports=dict(
        core_uses_compliant_datum=import_binding(PREFIX+"core_physical_field",PREFIX+"pressure_source","CompliantPressureDatum"),
        repair_alias_is_compliant_datum=import_binding(PREFIX+"five_moment_repair",PREFIX+"pressure_source","CompliantPressureDatum","LogarithmicPressureDatum"),
        pre_pulse_buffer_import=import_binding(PREFIX+"pre_pulse_mixed_C4",PREFIX+"outer_buffer","SharedOuterBuffer"),
        buffer_initial_import=import_binding(PREFIX+"outer_buffer",PREFIX+"outer_initial","SharedOuterInitial"),
        initial_repair_import=import_binding(PREFIX+"outer_initial",PREFIX+"five_moment_repair","SharedFiveMomentRepair"),
        inherited_q_jet_import=import_binding("lei_ren_part1_paper_logarithmic_pressure_datum",
            "lei_ren_part1_paper_candidate_pressure_function","q_jets"))
    constructors=dict(
        core=class_assignment("core_physical_field","CompliantCorePhysicalField","__init__","self.datum","CompliantPressureDatum('40',160)"),
        repair=class_assignment("five_moment_repair","SharedFiveMomentRepair","__init__","self.datum","LogarithmicPressureDatum('40',precision=160)"),
        buffer=class_assignment("outer_buffer","SharedOuterBuffer","__init__","self.initial","SharedOuterInitial()"),
        initial=class_assignment("outer_initial","SharedOuterInitial","__init__","self.repair","SharedFiveMomentRepair()"),
        initial_datum=class_assignment("outer_initial","SharedOuterInitial","__init__","self.datum","self.repair.datum"),
        pre_pulse=class_assignment("pre_pulse_mixed_C4","CompliantPrePulseMixedC4","__init__","self.buffer","SharedOuterBuffer()"),
        pre_pulse_initial=class_assignment("pre_pulse_mixed_C4","CompliantPrePulseMixedC4","__init__","self.initial","self.buffer.initial"),
        pre_pulse_datum=class_assignment("pre_pulse_mixed_C4","CompliantPrePulseMixedC4","__init__","self.datum","self.initial.datum"))
    compliant=external_assignment(PREFIX+"pressure_source","CompliantOuterParameters","__init__",{
        "self.log_epsilon":"c.ln(c.mpf('.001'))+self.log_delta",
        "self.epsilon":"c.exp(self.log_epsilon)"})
    stages=external_assignment(PREFIX+"pressure_source","CompliantPressureDatum","__init__",{
        "self.parameters":"CompliantOuterParameters(Md,precision)",
        "self.stages":"{n:dict(beta=2,mass=box) for n in BETA2}",
        "self.stages['reference_extension']['mass']":"c.mpf('2.5')",
        "self.stages['slope_transition_ref']['mass']":"read_interval(c,refined['stages']['slope_transition_ref']['refined_mass'])",
        "self.stages['axial_turnoff']['mass']":"c.exp(c.mpf('.6'))*(c.exp(-1)-c.exp(-p.yd))/2",
        "self.stages['z_flatten']":"dict(beta='variable [0,2]',mass=box)",
        "self.m2":"sum((self.stages[n]['mass'] for n in BETA2),c.mpf(0))",
        "self.m0":"sum((self.stages[n]['mass'] for n in BETA0),c.mpf(0))",
        "self.rho":"c.mpf('.25')",
        "self.flatten_complex_upper":"self.tail_upper/(1-self.rho**2)**2"})
    cls=next(n for n in tree(PREFIX+"pressure_source").body
        if isinstance(n,ast.ClassDef) and n.name=="CompliantPressureDatum")
    update=ast.dump(ast.parse("self.stages.update({n:dict(beta=0,mass=box) for n in BETA0})",mode="eval").body)
    if sum(isinstance(n,ast.Call) and ast.dump(n)==update for n in ast.walk(cls))!=1:
        raise ValueError("Common beta0 stage sources changed")
    # The three remainder assignments are mutually exclusive branches;
    # bind their exact occurrences rather than a single assignment target.
    fn=next(n for n in ast.walk(tree("lei_ren_part1_paper_logarithmic_pressure_datum"))
        if isinstance(n,ast.FunctionDef) and n.name=="normalized_jets")
    expected=["self.stages['z_flatten']['mass']","c.mpf(0)",
        "c.mpf([endpoints(-bound)[0],endpoints(bound)[1]])"]
    actual=[ast.dump(n.value) for n in ast.walk(fn) if isinstance(n,ast.Assign)
        and any(ast.unparse(v)=="remainder" for v in n.targets)]
    if sorted(actual)!=sorted(ast.dump(ast.parse(v,mode="eval").body) for v in expected):
        raise ValueError("Inherited flatten enclosure branches changed")
    expected_pressure_loop=ast.parse(
        "for n,a in enumerate(q_jets(c,z,order)):\n"
        "    if n==0: remainder=self.stages['z_flatten']['mass']\n"
        "    elif n%2 and endpoints(z)==(mp.mpf(0),mp.mpf(0)): remainder=c.mpf(0)\n"
        "    else:\n"
        "        bound=self.flatten_complex_upper/self.rho**n\n"
        "        remainder=c.mpf([endpoints(-bound)[0],endpoints(bound)[1]])\n"
        "    rows.append(-(self.m2*a+(self.m0 if n==0 else 0)+remainder))").body[0]
    if sum(isinstance(n,ast.For) and ast.dump(n)==ast.dump(expected_pressure_loop) for n in ast.walk(fn))!=1:
        raise ValueError("Normalized P0 prefix equation/parity/flatten source changed")
    expressions=["rows.append(-(self.m2*a+(self.m0 if n==0 else 0)+remainder))"]
    for expr in expressions:
        want=ast.dump(ast.parse(expr,mode="eval").body)
        if sum(isinstance(n,ast.Call) and ast.dump(n)==want for n in ast.walk(fn))!=1:
            raise ValueError("Inherited normalized P0 formula changed")
    qfn=next(n for n in ast.walk(tree("lei_ren_part1_paper_candidate_pressure_function"))
        if isinstance(n,ast.FunctionDef) and n.name=="q_jets")
    expected_loop=ast.parse("for n in range(order):\n    current=coefficients[n]\n    coefficients.append(-(b*(n+2)*current+(n+3)*previous)/(a*(n+1)))\n    previous=current").body[0]
    if sum(isinstance(n,ast.For) and ast.dump(n)==ast.dump(expected_loop) for n in ast.walk(qfn))!=1:
        raise ValueError("q_jets must remain order-independent prefix recurrence")
    qinit={}
    for target,expr in dict(a="1+center**2",b="2*center",coefficients="[1/a**2]",previous="ctx.mpf(0)").items():
        values=[n.value for n in ast.walk(qfn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(v)==ast.dump(ast.parse(expr,mode="eval").body) for v in values)!=1:
            raise ValueError("q_jets initialization changed")
        qinit[target]=True
    truncate=next(n for n in ast.walk(tree("lei_ren_part1_paper_interval_taylor"))
        if isinstance(n,ast.FunctionDef) and n.name=="truncate")
    wanted=ast.dump(ast.parse("IntervalTaylor(self.ctx,self.coefficients[:order+1])",mode="eval").body)
    if sum(isinstance(n,ast.Return) and ast.dump(n.value)==wanted for n in ast.walk(truncate))!=1:
        raise ValueError("Taylor truncation is not first order+1 coefficients")
    bridge=assignment_source_bindings("inner_bridge_profiles","inputs",{
        "datum":"self.core.datum.normalized_jets(endpoints(Z),6)",
        "pressure":"IntervalTaylor(c,[c.mpf(endpoints(value)) for value in datum['normalized_pressure_coefficients']])"})
    trace=actual_R100_trace_bindings()
    pre=assignment_source_bindings("pre_pulse_mixed_C4","packet",{
        "raw":"self.datum.normalized_jets(endpoints(c.mpf(Z)),5)['normalized_pressure_coefficients']",
        "p0":"IntervalTaylor(c,[c.mpf(endpoints(value)) for value in raw])"})
    patch=assignment_source_bindings("actual_moment_patch","evaluate",{
        "p0":"self.reference.inputs(Z)['original_axis_pressure']",
        "P":"p0+square(am)*pressure_moment"})
    return dict(imports=imports,canonical_constructor_paths=constructors,
        epsilon_point001=compliant,common_fourteen_stage_mass_and_flatten_sources=stages,
        same_inherited_normalized_jets_formula_bound=True,same_three_flatten_enclosure_branches_bound=True,
        q_jets_order_independent_prefix_recurrence_bound=True,q_jets_initialization=qinit,
        truncate5_means_coefficients_0_through_5=True,
        current_core_order6_pressure=bridge,current_switch_first_six_projection=trace,
        pre_pulse_direct_order5_pressure=pre,current_patch_retained_axis_pressure=patch)


def pressure_defining_function_proof(left,right):
    if (type(left) is not CompliantPressureDatum or type(right) is not CompliantPressureDatum
            or repair_module.LogarithmicPressureDatum is not CompliantPressureDatum
            or core_module.CompliantPressureDatum is not CompliantPressureDatum
            or CompliantPressureDatum.normalized_jets is not LogarithmicPressureDatum.normalized_jets
            or type(left.parameters) is not CompliantOuterParameters
            or type(right.parameters) is not CompliantOuterParameters):
        raise ValueError("Both Rh owners must use exactly the compliant inherited callable")
    attributes=("definition","source_sha","datum_sha","input_hashes")
    if any(getattr(left,key)!=getattr(right,key) for key in attributes):
        raise ValueError("Same analytic pressure source definition and schedule required")
    if (left.definition["c_epsilon"]!=".001" or left.definition["Md"]!="40"
            or left.ctx.dps!=160 or right.ctx.dps!=160
            or set(left.stages)!=set(BETA2+BETA0+("z_flatten",))):
        raise ValueError("Pinned .001/40/160 fourteen-stage pressure source required")
    def canonical(value):
        return encode(pack(value))
    for key in ("stages","m2","m0","rho","flatten_complex_upper","tail_upper"):
        if canonical(getattr(left,key))!=canonical(getattr(right,key)):
            raise ValueError("Same pressure mass/remainder enclosure construction required: "+key)
    if endpoints(left.stages["z_flatten"]["mass"])[1]<=0:
        raise ValueError("Implicit flatten remainder cannot be zeroed")
    bindings=pressure_function_source_bindings()
    z=s.symbols("Z",real=True);q=(1+z*z)**-2
    qrows=[s.diff(q,z,n)/s.factorial(n) for n in range(7)]
    recurrence=[]
    previous=s.Integer(0)
    for n in range(6):
        actual=-((2*z)*(n+2)*qrows[n]+(n+3)*previous)/((1+z*z)*(n+1))
        if s.simplify(actual-qrows[n+1])!=0:raise ArithmeticError("Pressure q Taylor recurrence failed")
        recurrence.append(n);previous=qrows[n]
    # One unknown analytic flatten function: boxes bound its coefficients,
    # never define six independently selectable source values.
    m2,m0=s.symbols("M2 M0",real=True);f=s.Function("F_flat")(z)
    direct=[-(m2*qrows[n]+(m0 if n==0 else 0)+s.diff(f,z,n)/s.factorial(n)) for n in range(6)]
    extended=[-(m2*qrows[n]+(m0 if n==0 else 0)+s.diff(f,z,n)/s.factorial(n)) for n in range(7)]
    if any(s.simplify(a-b)!=0 for a,b in zip(direct,extended[:6])):
        raise ArithmeticError("First-six defining-function projection failed")
    diagnostics=[]
    for center in (0,".3",[-1,1]):
        full=left.normalized_jets(center,6)["normalized_pressure_coefficients"]
        projected=IntervalTaylor(left.ctx,[left.ctx.mpf(endpoints(v)) for v in full]).truncate(5)
        a=left.normalized_jets(center,5)["normalized_pressure_coefficients"]
        b=right.normalized_jets(center,5)["normalized_pressure_coefficients"]
        for n in range(6):
            if endpoints(projected[n])!=endpoints(a[n]) or endpoints(a[n])!=endpoints(b[n]):
                raise ArithmeticError("Inherited prefix enclosure diagnostic changed")
        diagnostics.append(dict(center=center,first_six_coefficients_checked=6,actual_IntervalTaylor_order6_truncate5_used=True))
    return dict(source_AST_bindings=bindings,compliant_epsilon=".001",
        exact_same_inherited_normalized_jets_callable=True,
        exact_analytic_defining_function="P0/Pstar^2=-(M2*(1+Z^2)^-2+M0+F_flat(Z))",
        normalized_units="P/Pstar^2",same_fourteen_implicit_stage_sources=True,
        common_flatten_function_retained_not_replaced_by_zero=True,
        same_constructor_source_arguments_Md40_precision160=True,
        q_recurrence_symbolic_rows_checked=len(recurrence),
        exact_first_six_Taylor_defining_coefficients_projected=6,
        order5_returns_six_coefficients_order6_returns_seven=True,
        actual_first_six_enclosure_prefix_diagnostics=diagnostics,
        output_enclosures_not_selected_pressure_values=True,
        callable_constructor_and_equations_not_hash_only_proof=True,passed=True)


def Rh_functional_identity():
    y,z=s.symbols("y Z",real=True);Pstar=s.symbols("Pstar",positive=True)
    delta=s.symbols("delta",real=True);P0=s.Function("P0")(z)
    x=s.exp(y);offset=y-6;am=s.exp(-s.Rational(3,5))/(1+z*z)
    u=am*x**s.Rational(1,10);V=4*z;mass=V
    theta=s.Rational(5,8)*x**s.Rational(8,5)
    mixed=4*z*theta;energy=-s.Rational(5,12)*x**s.Rational(6,5)
    pressure=s.Rational(5,2)*x**s.Rational(1,5);G=x*mass
    Q=(2*z*V-(1-delta)*z*mass-(1-z*z)*s.diff(mass,z))/(1-delta*z*z)
    left=dict(Utheta_over_Pstar=u,Uz=V,Ur_over_sqrt_Rm_over_2=s.sqrt(x)*Q,
        P_over_Pstar2=P0+am*am*pressure,
        Mtheta_over_sqrt2_Rm_1p5_Pstar=am*theta,
        Mtheta_z_over_sqrt2_Rm_1p5_Pstar=am*mixed,Mz_over_Rm=G,
        Mztheta_over_Rm_Pstar2=(8*z*G-16*z*z*x)/Pstar**2+am*am*energy,
        Mp_over_Pstar2=am*am*pressure)
    ur=s.exp(offset/10)/(1+z*z);hr=s.Rational(5,8)*ur
    kr=hr*V;er=V*V/Pstar**2-s.Rational(5,12)*ur*ur;pr=s.Rational(5,2)*ur*ur
    Qr=(2*z*V-(1-delta)*z*V-(1-z*z)*s.diff(V,z))/(1-delta*z*z)
    right=dict(Utheta_over_Pstar=ur,Uz=V,Ur_over_sqrt_Rm_over_2=s.sqrt(x)*Qr,
        P_over_Pstar2=P0+pr,
        Mtheta_over_sqrt2_Rm_1p5_Pstar=x**s.Rational(3,2)*hr,
        Mtheta_z_over_sqrt2_Rm_1p5_Pstar=x**s.Rational(3,2)*kr,
        Mz_over_Rm=x*V,Mztheta_over_Rm_Pstar2=x*er,Mp_over_Pstar2=pr)
    rows={}
    for name,value in left.items():
        diff=s.simplify(s.powsimp(value-right[name],force=True))
        if diff!=0:raise ArithmeticError("Current external Rh physical source identity failed: "+name)
        rows[name]=["y"+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)]
    bindings=dict(
        original_reference_amplitude_and_histories=assignment_source_bindings("pre_pulse_mixed_C4","reference",{
            "u":"qi*c.exp(y/10)","V":"z*4","h":"u*c.mpf('.625')",
            "hist":"dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(u)*c.mpf(5)/12,p=square(u)*c.mpf('2.5'))"}),
        reference_radius_keyword=keyword_binding("pre_pulse_mixed_C4","reference","exact_radius_source",
            "'R=Rref*exp(offset), Rh=e^-5 Rref'"),
        original_patch_final_histories=assignment_source_bindings("actual_moment_patch","evaluate",{
            "terminal":"endpoints(x)[0]>=mp.mpf(71)/40",
            "mass":"z*4+defects[0]/x","theta":"defects[2]+x**c.mpf('1.6')*c.mpf('.625')",
            "mixed":"(z*theta)*4+defects[1]","energy":"defects[3]-x**c.mpf('1.2')*c.mpf(5)/12",
            "pressure_moment":"defects[4]+x**c.mpf('.2')*c.mpf('2.5')"}))

    bindings["original_patch_physical_unit_wrapper"]=assignment_source_bindings("actual_patch_mixed_C4","patch_mixed",{
        "U":"[am*row for row in H]",
        "P":"[p0+square(am)*pressure[0]]+[square(am)*row for row in pressure[1:]]",
        "physical_primitives":"dict(Mz_over_Rm=G,Mtheta_over_sqrt2_Rm_1p5_Pstar=[am*row for row in theta],Mtheta_z_over_sqrt2_Rm_1p5_Pstar=[am*row for row in mixed],Mztheta_over_Rm_Pstar2=[((z*G[k])*8-square(z)*(16*xder[k]))*invP2+square(am)*energy[k] for k in range(5)],Mp_over_Pstar2=[square(am)*row for row in pressure])"})
    bindings["original_reference_physical_unit_wrapper"]=assignment_source_bindings("pre_pulse_mixed_C4","physical_mixed",{
        "physical":"dict(Utheta_over_Pstar=U,Uz=V,Ur_over_current_sqrt_R_over_2=[rate_rows(Q,c.mpf('.5'),j) for j in range(5)],P_over_Pstar2=[p0+rows['p'][0]]+rows['p'][1:])",
        "primitives":"dict(Mz_over_current_R=[rate_rows(rows['m'],1,j) for j in range(5)],Mtheta_over_current_sqrt2_R_1p5_Pstar=[rate_rows(rows['h'],c.mpf('1.5'),j) for j in range(5)],Mtheta_z_over_current_sqrt2_R_1p5_Pstar=[rate_rows(rows['k'],c.mpf('1.5'),j) for j in range(5)],Mztheta_over_current_R_Pstar2=[rate_rows(rows['e'],1,j) for j in range(5)],Mp_over_Pstar2=rows['p'])"})
    return dict(exact_physical_function_identities=list(left),
        physical_mixed4_rows_implied_by_exact_function_identities=rows,total_implied_physical_mixed4_rows=135,
        source_AST_bindings=bindings,coordinate_identity="x=R/Rm=exp(y), offset=log(R/Rref)=y-6",
        Rh_identity="x=e, y=1, offset=-5; Rh=e*Rm=e^-5*Rref",
        original_last_support_upper_edge="71/40 < e",
        canonical_extension_equals_source_on_each_side_near_Rh=True,
        current_unique_implicit_map_full_weight_closure_used=True,
        physical_components_and_all_five_primitives_share_P0=True,
        log_radial_operator="D_y=x D_x; higher orders use Stirling combinations",
        right_grid_Ur_fixed_unit_factor="sqrt(e)",right_grid_Mz_fixed_unit_factor="e",
        right_grid_angular_fixed_unit_factor="exp(3/2)",
        right_grid_Mztheta_fixed_unit_factor="e",
        common_units_used_after_physical_differentiation=True,
        normalization_multipliers_frozen_at_basepoint_not_differentiated_in_Z=True,
        no_new_pressure_fit_or_moment_reset=True,interval_overlap_is_not_functional_proof=True,passed=True)


class CurrentRhReferenceDispatcher(CurrentMatchedSourceDispatcher):
    CHARTS=tuple(CURRENT_ROUTES)+("Rh_reference",)

    def __init__(self,require_checked=True):
        super().__init__(require_checked=True)
        self.require_Rh_checked=require_checked;self.rh_reference=None
        self.pressure_proof=None;self.Rh_proof=None;self.Rh_acceptance_loaded=False

    @source_precision
    def _load_Rh(self):
        if self.rh_reference is not None:return
        super()._load()
        self.rh_reference=pre=CompliantPrePulseMixedC4()
        core=self.anchor.patch.core
        if (self.join_proof["reference_side_coordinate"]!="y=log(R/Rm)=offset+6"
                or self.join_proof["patch_side_coordinate"]!="x=R/Rm=exp(y)"
                or not self.join_proof["passed"]):
            raise ValueError("Current accepted Rm radius-coordinate source relation required")
        if (pre.family!=self.family or pre.source!=self.source
                or pre.datum.datum_sha!=self.datum_sha
                or pre.datum is not pre.initial.repair.datum):
            raise ValueError("Current Rh analytic reference owner has another source")
        for path,digest in pre.hashes.items():
            if path in self.hashes and self.hashes[path]!=digest:
                raise ValueError("Current Rh common dependency differs: "+path)
        self.hashes.update(pre.hashes)
        name=PREFIX+"pre_pulse_mixed_C4_check.json"
        check=accepted(name,self.family,self.source,"all_pre_pulse_mixed4_available")
        if (check["datum_enclosure_sha256"]!=self.datum_sha
                or not check["independent_closed_integral_fixture"]["passed"]
                or not check["exact_original_production_and_join_identities"]["original_reference_p"]):
            raise ValueError("Unchanged outer reference equations/physical fixture required")
        self.hashes.update(check["input_hashes"]);self.hashes[name]=sha(name)
        if not (self.receipts["actual_feedback_patch_mixed_C4_check"]["current_Rh_same_unique_implicit_function_and_full_support_retained"]
                and self.receipts["actual_feedback_patch_mixed_C4_check"]["current_five_functional_terminal_identities_connected"]):
            raise ValueError("Current actual full-support functional closure required")
        self.pressure_proof=pressure_defining_function_proof(core.datum,pre.datum)
        self.Rh_proof=Rh_functional_identity()
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        if self.require_Rh_checked:
            name=PREFIX+"actual_Rh_source_join_check.json"
            check=accepted(name,self.family,self.source,"current_Rh_external_neighbor_join_certified")
            if check["datum_enclosure_sha256"]!=self.datum_sha:
                raise ValueError("Current external Rh acceptance has another datum")
            self.hashes.update(check["input_hashes"]);self.hashes[name]=sha(name)
            self.Rh_acceptance_loaded=True

    def provider(self,chart):
        if chart=="Rh_reference":
            self._load_Rh();return self.rh_reference
        return super().provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart!="Rh_reference":return super().evaluate(chart,Z,coordinate)
        provider=self.provider(chart);packet=provider.reference(Z,coordinate)
        grids={key:packet[key] for key in ("physical_velocity_pressure_y_Z_mixed4","physical_five_primitive_y_Z_mixed4")}
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,
            implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            source_provider=RH_ROUTE["provider"],acceptance_receipt=RH_ROUTE["acceptance_receipt"],
            source_coverage_coordinate=RH_ROUTE["coverage_coordinate"],source_coordinate_domain=RH_ROUTE["domain"],
            derivative_coordinate="ordinary y=logR,Z; offset only selects radius",
            original_scale_metadata=packet["exact_formal_prefactors"],physical_mixed_grids=grids,source_packet=packet,
            current_Rh_incoming_unique_implicit_map_owner=PREFIX+"actual_feedback_patch_mixed_C4",
            current_Rh_external_neighbor_join_certified=self.Rh_acceptance_loaded,
            current_Rh_external_neighbor_source_join_proved=True,current_pressure_defining_function_and_projection_retained=True,
            output_kind="directed source-field enclosures; no production point selection",
            **dict.fromkeys(SCOPES,False))

    def manifest(self):
        result=super().manifest();self._load_Rh()
        if result["current_Rh_external_neighbor_join_certified"]:
            raise ValueError("Native eight-chart external Rh flag must remain false")
        result["ordered_current_chart_registry"]["Rh_reference"]=RH_ROUTE
        result.update(current_Rh_pressure_defining_function_proof=self.pressure_proof,
            current_external_Rh_functional_join=self.Rh_proof,
            current_Rh_external_neighbor_join_certified=self.Rh_acceptance_loaded,
            current_Rh_external_neighbor_source_join_proved=True,
            current_Rm_coordinate_proof_consumed=True,
            current_Rh_reference_owner_installed=True,current_downstream_chart_owner_count=9,
            native_eight_chart_Rh_false_flag_preserved=True,
            all_profile_source_charts_callable=False,uniform_physical_units_assembled=False,
            **dict.fromkeys(SCOPES,False),input_hashes=dict(self.hashes))
        return result


@source_precision
def build():
    field=CurrentRhReferenceDispatcher(require_checked=False)
    result=field.manifest()
    result["whole_current_Rh_reference_chart"]=field.evaluate("Rh_reference",[-1,1],[-5,0])
    result["input_hashes"]=dict(field.hashes)
    return encode(pack(result))


def run():
    result=build()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("Current Rh external reference source join generated: shared P0 first-six projection, 9 physical identities, 135 implied rows",flush=True)
    return result


if __name__=="__main__":
    run()

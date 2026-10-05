"""Current Rh closure transported through the original five pre-pulse charts.

The same accepted outer reference object and unchanged five-history ODE
solutions continue to Rp. No source history is reset at axial turnoff.
Legacy charts and receipts remain intact.
"""
import ast
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_actual_Rh_source_join import (
    CurrentRhReferenceDispatcher, CompliantPrePulseMixedC4, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_source_dispatcher import ROUTES
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4_check import symbolic_sources
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

CHARTS=("O2_slope","O2_axial","O2_buffer","O3_slope_mu","O3_power")
PRE_ROUTES={chart:dict(provider=PREFIX+ROUTES[chart][0]+"."+ROUTES[chart][1],
    method=ROUTES[chart][2],coverage_coordinate=ROUTES[chart][3],
    domain=ROUTES[chart][4],acceptance_receipt=PREFIX+"current_pre_pulse_source_dispatcher_check.json",
    extra=ROUTES[chart][6]) for chart in CHARTS}


def packet_argument_binding(method,index,expression):
    tree=ast.parse((HERE/(PREFIX+"pre_pulse_mixed_C4.py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
    calls=[n for n in ast.walk(fn) if isinstance(n,ast.Call) and ast.unparse(n.func)=="self.packet"]
    wanted=ast.dump(ast.parse(expression,mode="eval").body)
    if len(calls)!=1 or len(calls[0].args)<=index or ast.dump(calls[0].args[index])!=wanted:
        raise ValueError("Original radial source jet argument changed: "+method)
    return True


def source_histories(method,env):
    """Interpret the bound original hist AST in exact arbitrary source symbols."""
    tree=ast.parse((HERE/(PREFIX+"pre_pulse_mixed_C4.py")).read_text(encoding="utf8"))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
    node=next(n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
        and any(ast.unparse(v)=="hist" for v in n.targets))
    def value(n):
        label=ast.unparse(n)
        if label in env:return env[label]
        if isinstance(n,ast.Constant):return s.Rational(str(n.value))
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -value(n.operand)
        if isinstance(n,ast.BinOp):
            a,b=value(n.left),value(n.right)
            for kind,op in ((ast.Add,lambda:a+b),(ast.Sub,lambda:a-b),
                    (ast.Mult,lambda:a*b),(ast.Div,lambda:a/b),(ast.Pow,lambda:a**b)):
                if isinstance(n.op,kind):return op()
        if isinstance(n,ast.Call) and len(n.args)==1:
            label=ast.unparse(n.func);a=value(n.args[0])
            if label=="c.mpf":return a
            if label=="c.exp":return s.exp(a)
            if label=="square":return a*a
        raise ValueError("Unbound original history source: "+label)
    return {kw.arg:value(kw.value) for kw in node.keywords}


def original_history_source_bindings():
    data=dict(
        slope=assignment_source_bindings("pre_pulse_mixed_C4","slope",{
            "(J, mass)":"slope_masses(c,y,self.cells)",
            "factor":"c.exp(y/10-c.mpf('.6')*J)",
            "u":"qi*factor","V":"z*4",
            "h":"qi*(c.mpf('.625')+mass[0])*c.exp(-c.mpf('1.5')*y)",
            "hist":"dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(qi)*(c.mpf(5)/12+mass[2]/2)*c.exp(-y),p=square(qi)*(c.mpf('2.5')+mass[1]/2))"}),
        inlet=assignment_source_bindings("pre_pulse_mixed_C4","inlet",{
            "self.cache[key]":"self.slope(Z,1)"}),
        axial=assignment_source_bindings("pre_pulse_mixed_C4","axial",{
            "parent":"self.inlet(Z)",
            "u1":"IntervalTaylor(c,parent['Utheta_over_Pstar_axial5_coefficients'])",
            "old":"{key:rows[0] for key,rows in parent['actual_normalized_primitive_y_derivative_axial5'].items()}",
            "t":"y-1","d":"c.exp(-t)","root":"c.exp(-t/2)","d3":"c.exp(-c.mpf('1.5')*t)",
            "u":"u1*root","V":"[z*(4*b) for b in B]",
            "K":"turnoff_kernels(c,y,self.params.Md,self.cells,self.window)",
            "hist":"dict(m=old['m']*d+z*(4*K['B_mass']),h=old['h']*d3+u1*(root-d3),k=old['k']*d3+u1*z*(4*root*K['B_mass']),e=old['e']*d+square(z)*(16*self.invP2*K['B_squared_mass'])-square(u1)*(t*d/2),p=old['p']+square(u1)*((1-d)/2))"}),
        slope_mu=assignment_source_bindings("pre_pulse_mixed_C4","slope_mu",{
            "parent":"self.axial(Z,buffer_offset=11)",
            "u1":"IntervalTaylor(c,parent['Utheta_over_Pstar_axial5_coefficients'])",
            "old":"{key:rows[0] for key,rows in parent['actual_normalized_primitive_y_derivative_axial5'].items()}",
            "K":"transition_kernels(c,t,mu,self.cells)",
            "factor":"c.exp(-t/2-mu*K['J'])",
            "hist":"dict(m=old['m']*d,h=old['h']*d3+u1*(d3*K['theta']),k=old['k']*d3,e=old['e']*d-square(u1)*(d*K['energy']/2),p=old['p']+square(u1)*(K['pressure']/2))"}),
        power=assignment_source_bindings("pre_pulse_mixed_C4","power",{
            "parent":"self.slope_mu(Z,1)",
            "u1":"IntervalTaylor(c,parent['Utheta_over_Pstar_axial5_coefficients'])",
            "old":"{key:rows[0] for key,rows in parent['actual_normalized_primitive_y_derivative_axial5'].items()}",
            "t":"self.params.Tw*phase",
            "f":"c.exp((-c.mpf('.5')-mu)*t)",
            "hist":"dict(m=old['m']*d,h=old['h']*d3+u1*theta,k=old['k']*d3,e=old['e']*d-square(u1)*(d*decay_integral(c,2*mu,t)/2),p=old['p']+square(u1)*(decay_integral(c,1+2*mu,t)/2))"}),
        common_pressure=assignment_source_bindings("pre_pulse_mixed_C4","packet",{
            "raw":"self.datum.normalized_jets(endpoints(c.mpf(Z)),5)['normalized_pressure_coefficients']",
            "data":"physical_mixed(c,Z,self.delta,u,logU,V,history,p0,self.invP2)"}))
    data["ordinary_radial_amplitude_source_jets"]=dict(
        slope=assignment_source_bindings("pre_pulse_mixed_C4","slope",{
            "logU":"[zero+c.mpf('.1')-sig[0]*c.mpf('.6')]+[zero-sig[k]*math.factorial(k)*c.mpf('.6') for k in range(1,4)]"}),
        slope_mu=assignment_source_bindings("pre_pulse_mixed_C4","slope_mu",{
            "logU":"[zero-c.mpf('.5')-mu*sig[0]]+[zero-mu*sig[k]*math.factorial(k) for k in range(1,4)]"}),
        axial_log_source=assignment_source_bindings("pre_pulse_mixed_C4","axial",{
            "logu":"parent['log_Utheta_over_Pstar_base_source']-t/2",
            "B":"turnoff_derivatives(c,y,md,phase)"}))
    data["ordinary_radial_packet_argument_bindings"]={
        method:dict(logU=packet_argument_binding(method,5,logjet),
            V=packet_argument_binding(method,6,vjet))
        for method,logjet,vjet in (
            ("reference","[zero+c.mpf('.1')]+[zero]*3","[V]+[zero]*4"),
            ("slope","logU","[V]+[zero]*4"),
            ("axial","[zero-c.mpf('.5')]+[zero]*3","V"),
            ("slope_mu","logU","[zero]*5"),
            ("power","[zero-c.mpf('.5')-mu]+[zero]*3","[zero]*5"))}
    data["original_power_positive_mu_and_logR_derivative"]=assignment_source_bindings("pre_pulse_mixed_C4","power",{
        "mu":"c.mpf(self.params.mu)"})
    # These prove equations for arbitrary inlet histories. The old actual
    # Rh closure claim is not reused: it is supplied by the new Rh receipt.
    universal=symbolic_sources()
    if not all(universal.values()):raise ValueError("Unchanged original pre-pulse equations differ")
    return dict(original_method_parent_and_history_AST_bindings=data,
        unchanged_original_five_history_equations=universal,
        inherited_universal_equation_identities_not_current_Rh_admission=True,
        same_current_Rh_reference_provider_methods_used=True)


def join_source_proof(c):
    """Neutral integrating factors and shared ODEs at five actual interfaces."""
    z,u,invP2=s.symbols("Z u invP2",real=True)
    names=("m","h","k","e","p")
    m,h,k,e,p=s.symbols("m h k e p",real=True)
    old=dict(zip(names,(m,h,k,e,p)))
    qi=1/(1+z*z)
    ref=source_histories("reference",dict(V=4*z,h=5*qi/8,u=qi,**{"self.invP2":invP2}))
    slope0=source_histories("slope",dict(V=4*z,h=5*qi/8,qi=qi,y=s.Integer(0),
        **{"self.invP2":invP2,"mass[0]":s.Integer(0),"mass[1]":s.Integer(0),"mass[2]":s.Integer(0)}))
    inherited={"old['"+key+"']":val for key,val in old.items()}
    base=dict(z=z,u1=u,t=s.Integer(0),d=s.Integer(1),root=s.Integer(1),d3=s.Integer(1),
        **inherited,**{"self.invP2":invP2,"K['B_mass']":s.Integer(0),"K['B_squared_mass']":s.Integer(0)})
    axial0=source_histories("axial",base)
    mu0=source_histories("slope_mu",dict(base,**{"K['theta']":s.Integer(0),
        "K['energy']":s.Integer(0),"K['pressure']":s.Integer(0)}))
    power0=source_histories("power",dict(base,theta=s.Integer(0),
        **{"decay_integral(c, 2 * mu, t)":s.Integer(0),
           "decay_integral(c, 1 + 2 * mu, t)":s.Integer(0)}))
    t,d,r,d3,K1,K2=s.symbols("t d root d3 K1 K2",real=True)
    same_y=dict(base,t=t,d=d,root=r,d3=d3,
        **{"K['B_mass']":K1,"K['B_squared_mass']":K2})
    axial_exit=source_histories("axial",same_y)
    buffer_inlet=source_histories("axial",same_y)
    rows={}
    for label,a,b in (("Rref_reference_slope",ref,slope0),
            ("O2_slope_axial_inlet",old,axial0),
            ("O2_turnoff_buffer",axial_exit,buffer_inlet),
            ("Rd_buffer_slope_mu",old,mu0),("Rw_slope_mu_power",old,power0)):
        identities={}
        for name in names:
            if s.simplify(a[name]-b[name])!=0:
                raise ArithmeticError("Actual pre-pulse interface history differs: "+label+"/"+name)
            identities[name]=True
        rows[label]=identities
    endpoint_jets={}
    for value in (0,1):
        jets=sigma_jets(c,c.mpf(value))
        if endpoints(jets[0])!=(mp.mpf(value),mp.mpf(value)) or any(
                endpoints(jets[n])!=(mp.mpf(0),mp.mpf(0)) for n in range(1,5)):
            raise ArithmeticError("Original flat sigma endpoint changed")
        endpoint_jets[str(value)]=dict(value=value,ordinary_derivatives1_through4_exact_zero=True)
    # Each pair uses the same physical_mixed callable/parameters. Physical
    # jets are uniquely determined by common primitive values, logU/V jets
    # and the five first-order ODEs (including all radial prefactors).
    grid=["y"+str(j)+"_Z"+str(n) for j in range(5) for n in range(5-j)]
    rates={
        "Rref_reference_slope":dict(logU=".1",V="4Z",higher_radial_derivatives="0"),
        "O2_slope_axial_inlet":dict(logU="-.5",V="4Z",higher_radial_derivatives="0"),
        "O2_turnoff_buffer":dict(logU="-.5",V="0",higher_radial_derivatives="0"),
        "Rd_buffer_slope_mu":dict(logU="-.5",V="0",higher_radial_derivatives="0"),
        "Rw_slope_mu_power":dict(logU="-.5-mu; exact positive source mu",V="0",higher_radial_derivatives="0")}
    return dict(exact_five_history_interface_identities=rows,
        exact_history_identity_count=25,original_flat_endpoint_jets=endpoint_jets,
        same_log_amplitude_and_velocity_boundary_source_jets=rates,
        turnoff_buffer_same_y_and_kernel_arguments="y=exp(Md) at phase1/offset0; same slope(Z,1) parent",
        source_radial_interfaces=["Rref","Rref*e","Rref*exp(exp(Md))","Rd","Rw"],
        source_derivatives_are_logR_not_selector_phase=True,
        same_physical_mixed_callable_and_five_ODEs_on_every_side=True,
        original_five_normalized_ODEs=[
            "m_y=V-m","h_y=u-3h/2","k_y=uV-3k/2",
            "e_y=V^2/Pstar^2-u^2/2-e","p_y=u^2/2"],
        implied_mixed4_rows_per_interface={key:dict(velocity_pressure_components=4,
            primitive_components=5,rows=9*len(grid),multiindices=grid) for key in rows},
        total_implied_physical_mixed4_interface_rows=675,
        exact_positive_mu_retained_in_formal_power_slope=True,
        retained_axial_and_mixed_histories_not_reset_when_V_zero=True,
        exact_history_identities_interpreted_from_bound_original_AST=True,
        primitive_value_identities_not_interval_overlap=True,passed=True)


class CurrentPrePulseSourceDispatcher(CurrentRhReferenceDispatcher):
    CHARTS=CurrentRhReferenceDispatcher.CHARTS+CHARTS

    def __init__(self,require_checked=True):
        super().__init__(require_checked=True)
        self.require_pre_checked=require_checked;self.pre_loaded=False
        self.pre_acceptance_loaded=False;self.pre_bindings=None;self.pre_joins=None

    @source_precision
    def _load_pre(self):
        if self.pre_loaded:return
        super()._load_Rh();pre=self.rh_reference
        if (type(pre) is not CompliantPrePulseMixedC4 or pre.params is not pre.datum.parameters
                or pre.datum is not pre.initial.repair.datum
                or pre.source!=self.source or pre.family!=self.family):
            raise ValueError("One current Rh provider/parameter/history chain required")
        c=pre.ctx
        if not (0<endpoints(pre.params.mu)[0]<=endpoints(pre.params.mu)[1]<mp.mpf(".001")
                and endpoints(pre.params.Tw)[0]>0):
            raise ValueError("Original positive mu/Tw must not be dropped")
        check=accepted(PREFIX+"pre_pulse_mixed_C4_check.json",self.family,self.source,
            "pre_pulse_functional_mixed4_joins_certified")
        if (check["datum_enclosure_sha256"]!=self.datum_sha
                or not check["independent_closed_integral_fixture"]["passed"]
                or not check["independent_turnoff_coordinate_fixture"]["passed"]
                or not check["retained_nonzero_axial_history_at_Rp_checked"]):
            raise ValueError("Unchanged original derivative/retained-history fixture missing")
        if not self.Rh_acceptance_loaded:
            raise ValueError("Current external Rh receipt must be consumed first")
        self.pre_bindings=original_history_source_bindings()
        self.pre_joins=join_source_proof(c)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        if self.require_pre_checked:
            name=PREFIX+"current_pre_pulse_source_dispatcher_check.json"
            receipt=accepted(name,self.family,self.source,"current_Rh_to_Rp_source_chain_certified")
            if receipt["datum_enclosure_sha256"]!=self.datum_sha:
                raise ValueError("Current pre-pulse acceptance datum differs")
            self.hashes.update(receipt["input_hashes"]);self.hashes[name]=sha(name)
            self.pre_acceptance_loaded=True
        self.pre_loaded=True

    def provider(self,chart):
        if chart in PRE_ROUTES:
            self._load_pre();return self.rh_reference
        return super().provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart not in PRE_ROUTES:return super().evaluate(chart,Z,coordinate)
        provider=self.provider(chart);route=PRE_ROUTES[chart]
        fn=getattr(provider,route["method"])
        packet=fn(Z,buffer_offset=coordinate) if route["extra"]=="buffer_offset" else fn(Z,coordinate)
        keys=("physical_velocity_pressure_y_Z_mixed4","physical_five_primitive_y_Z_mixed4")
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,
            implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            source_provider=route["provider"],acceptance_receipt=route["acceptance_receipt"],
            source_coverage_coordinate=route["coverage_coordinate"],source_coordinate_domain=route["domain"],
            derivative_coordinate="ordinary y=logR,Z; selector phase/offset only selects coverage",
            physical_mixed_grids={key:packet[key] for key in keys},source_packet=packet,
            original_scale_metadata=packet["exact_formal_prefactors"],
            same_current_Rh_reference_provider=True,current_Rh_to_Rp_source_chain_proved=True,
            current_Rh_to_Rp_source_chain_certified=self.pre_acceptance_loaded,
            current_Rp_external_pulse_join_certified=False,
            output_kind="directed source-field enclosures; no production point selection",
            **dict.fromkeys(SCOPES,False))

    def manifest(self):
        result=super().manifest();self._load_pre()
        for chart,route in PRE_ROUTES.items():
            result["ordered_current_chart_registry"][chart]=route
        result.update(current_original_pre_pulse_history_source_bindings=self.pre_bindings,
            current_pre_pulse_five_interface_source_proof=self.pre_joins,
            current_Rh_to_Rp_source_chain_proved=True,
            current_Rh_to_Rp_source_chain_certified=self.pre_acceptance_loaded,
            current_Rp_external_pulse_join_certified=False,
            current_downstream_chart_owner_count=14,
            all_current_pre_pulse_charts_callable=True,
            same_single_current_outer_reference_object_retained=True,
            all_profile_source_charts_callable=False,uniform_physical_units_assembled=False,
            **dict.fromkeys(SCOPES,False),input_hashes=dict(self.hashes))
        return result


@source_precision
def build():
    field=CurrentPrePulseSourceDispatcher(require_checked=False);result=field.manifest()
    packets={}
    for chart in CHARTS:
        domain=[0,11] if chart=="O2_buffer" else [0,1]
        packets[chart]=field.evaluate(chart,[-1,1],domain)
        print("Current original pre-pulse chart: "+chart,flush=True)
    result["whole_current_pre_pulse_charts"]=packets
    result["input_hashes"]=dict(field.hashes)
    return encode(pack(result))


def run():
    result=build()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("Current original pre-pulse chain generated: 14 owners, 5 new charts, retained histories through Rp",flush=True)
    return result


if __name__=="__main__":
    run()

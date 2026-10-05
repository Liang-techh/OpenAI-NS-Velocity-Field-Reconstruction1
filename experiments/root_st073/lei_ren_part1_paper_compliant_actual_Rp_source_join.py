"""Direct current Rp-to-native-pulse source join and current entrance owner.

The native pulse obtains inlet constants from live SharedOuterBuffer
callables. Saved power-inlet samples are not its defining functions.
"""
import ast
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace
from functools import lru_cache

import mpmath as mp
import sympy as s

import lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 as pre_source
import lei_ren_part1_paper_compliant_outer_initial as initial_source
import lei_ren_part1_paper_compliant_outer_buffer as buffer_source
import lei_ren_part1_paper_compliant_axial_pulse_field as axial_source
import lei_ren_part1_paper_compliant_flat_pulse_derivatives as flat_source
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import CompliantPulseMixedC4
from lei_ren_part1_paper_compliant_current_pre_pulse_source_dispatcher import (
    CurrentPrePulseSourceDispatcher, source_histories, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

Z=s.Symbol("Z",real=True)
KEYS=("m","h","k","e","p")
LABELS=dict(m="Mz_over_R",h="Mtheta_over_sqrt2_R_3half_Pstar",
    k="Mtheta_z_over_sqrt2_R_3half_Pstar",e="Mztheta_over_R_Pstar_squared",
    p="Mp_over_Pstar_squared")


def keywords_binding(stem,method,target,expected):
    tree=ast.parse((HERE/(PREFIX+stem+".py")).read_text(encoding="utf8"))
    fn=next(v for v in ast.walk(tree) if isinstance(v,ast.FunctionDef) and v.name==method)
    nodes=[v.value for v in ast.walk(fn) if isinstance(v,ast.Assign)
        and any(ast.unparse(t)==target for t in v.targets)]
    if len(nodes)!=1 or not isinstance(nodes[0],ast.Call) or ast.unparse(nodes[0].func)!="dict":
        raise ValueError("Original source dictionary changed: "+stem+"."+target)
    actual={v.arg:ast.dump(v.value) for v in nodes[0].keywords}
    for key,expression in expected.items():
        if actual.get(key)!=ast.dump(ast.parse(expression,mode="eval").body):
            raise ValueError("Original source dictionary field changed: "+stem+"."+key)
    return dict.fromkeys(expected,True)


def source_recipe_bindings():
    result=dict(
        old_slope=assignment_source_bindings("outer_initial","slope",{
            "(J, integrals)":"transition_integrals(c,y,cells)"}),
        old_axial=assignment_source_bindings("outer_initial","axial",{
            "inlet":"self.slope(Z,1,cells)",
            "K":"turnoff_kernels(c,y,self.params.Md,cells,window)"}),
        old_mu=assignment_source_bindings("outer_buffer","slope_mu",{
            "inlet":"self.initial.axial(Z,buffer_offset=11,cells=cells)",
            "kernels":"transition_kernels(c,t,mu,cells)"}),
        old_power=assignment_source_bindings("outer_buffer","power",{
            "inlet":"self.slope_mu(Z,1,cells)","t":"self.params.Tw*phase",
            "theta_kernel":"(f-d3)/(1-mu)"}),
        current_power=assignment_source_bindings("pre_pulse_mixed_C4","power",{
            "theta":"(f-d3)/(1-mu)"}),
        actual_high_constants=assignment_source_bindings("axial_high_jets","_incoming_constants",{
            "initial":"self.base.pulse.initial","buffer":"self.base.pulse.buffer",
            "p0":"buffer.power(0,1,cells=128)",
            "U":"box(p0['Utheta_over_Pstar'][0])","M":"box(p0['Mz_over_R'][1])",
            "K":"box(p0['Mtheta_z_over_sqrt2_R_3half_Pstar'][1])",
            "EQ":"box(p0['Mztheta_over_R_Pstar_squared'][0])",
            "kernels":"turnoff_kernels(initial.ctx,initial.params.yd,initial.params.Md,cells=128)",
            "td":"box(initial.params.yd)-1",
            "EZ":"16*box(initial.invP2)*(c.exp(-td)+box(kernels['B_squared_mass']))",
            "invP":"c.exp(-box(initial.params.logPstar))"}),
        high_constant_units=keywords_binding("axial_high_jets","_incoming_constants","self.constants",{
            "U":"U","M":"M","K":"K","E_Q":"EQ","E_Z":"EZ",
            "C1":"invP*M/U","C2":"invP*K/(U*U)","C0":"EQ/(U*U)",
            "C_E":"EZ/(U*U)","Tw":"box(initial.params.Tw)"}),
        live_pulse_H_P=assignment_source_bindings("pulse_radial_C4","__init__",{
            "p0":"self.pulse.buffer.power(0,1)",
            "self.inlet_H":"box(p0['Mtheta_over_sqrt2_R_3half_Pstar'][0])",
            "self.inlet_P":"box(p0['Mp_over_Pstar_squared'][0])",
            "self.Xp":"self.inlet_H/box(p0['Utheta_over_Pstar'][0])"}),
        live_pulse_data=assignment_source_bindings("pulse_radial_C4","data",{
            "selected":"self.fifth.select(Z)","k":"self.high.constants",
            "u":"r*k['U']","incoming":"selected['incoming']",
            "result":"(selected,inlet,u,incoming['energy_Taylor'],incoming['moment_Taylor'])"}),
        same_pressure_function=assignment_source_bindings("pulse_radial_C4","pressure_moment",{
            "rows":"self.selection.future.angular.initial.datum.normalized_jets(c.mpf(Z),5)['normalized_pressure_coefficients']"}))
    tree=ast.parse((HERE/(PREFIX+"axial_high_jets.py")).read_text(encoding="utf8"))
    fn=next(v for v in ast.walk(tree) if isinstance(v,ast.FunctionDef) and v.name=="_incoming_constants")
    terms=[v for v in ast.walk(fn) if isinstance(v,ast.AugAssign) and ast.unparse(v.target)=="EZ"]
    if len(terms)!=1 or not isinstance(terms[0].op,ast.Mult) or ast.dump(terms[0].value)!=ast.dump(
            ast.parse("c.exp(-1)*c.exp(-box(initial.params.Tw))",mode="eval").body):
        raise ValueError("Actual positive E_Z source multiplier changed")
    result["positive_E_Z_terminal_multiplier_bound"]=True
    result["native_constructor_path"]=assignment_source_bindings("pulse_radial_C4","__init__",{
        "self.fifth":"CompliantFifthAxialJets()",
        "self.high":"SimpleNamespace(ctx=c,base=self.fifth.fourth.base,constants=self.fifth.fourth.constants,select=self.fifth.select,energy=SimpleNamespace(future=self.fifth.future))",
        "self.selection":"_SelectedSource(self.high)","self.pulse":"self.high.base.pulse",
        "self.mu":"box(self.high.base.mu)","self.delta":"box(self.high.base.future.delta)"})
    result["native_fifth_prefix_source"]=assignment_source_bindings("fifth_axial_jets","select",{
        "old":"self.fourth.select(Z)","incoming":"old['incoming']","k":"self.fourth.constants",
        "moments":"[append_fifth(c,j,c.mpf(0)) for j in incoming['moment_Taylor']]",
        "energy":"append_fifth(c,incoming['energy_Taylor'],6*Z*k['C_E'])"})
    result["actual_entrance_mixed_packet"]=assignment_source_bindings("pulse_mixed_C4","_high_packet",{
        "(selected, _, u, _, _)":"self.data(Z)","coord":"point['coordinate']","kind":"coord['kind']",
        "shape":"gp_jets(c,coord['xi'])","ap":"selected['selected_ap_Taylor']",
        "Brows":"[ap*(shape[k]*math.factorial(k)*self.mu**k) for k in range(5)]",
        "mixed":"transport_mixed(self,Z,point,Brows,u)"})
    result["native_pressure_primitive"]=assignment_source_bindings("pulse_radial_C4","pressure_moment",{
        "p":"IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)",
        "p0":"IntervalTaylor(c,[c.mpf(endpoints(v)) for v in rows])"})
    pressure_tree=ast.parse((HERE/(PREFIX+"pulse_radial_C4.py")).read_text(encoding="utf8"))
    pressure_fn=next(v for v in ast.walk(pressure_tree) if isinstance(v,ast.FunctionDef) and v.name=="pressure_moment")
    returned=next(v.value for v in ast.walk(pressure_fn) if isinstance(v,ast.Return))
    fields={kw.arg:ast.dump(kw.value) for kw in returned.keywords}
    for key,expr in dict(P_over_Pstar_squared="p+p0",P_y_over_Pstar_squared="u*u*(decay/2)",
            Mp_over_Pstar_squared="p",P0_over_Pstar_squared="p0").items():
        if fields.get(key)!=ast.dump(ast.parse(expr,mode="eval").body):
            raise ValueError("Actual native absolute pressure publication changed")
    result["native_absolute_pressure_publication_bound"]=True
    result["scalar_velocity_sources"]=dict(
        slope=assignment_source_bindings("outer_initial","slope",{
            "factor":"c.exp(yy/10-c.mpf('.6')*J)","u":"qi*factor"}),
        axial=assignment_source_bindings("outer_initial","axial",{"u":"u1*root_decay"}),
        slope_mu=assignment_source_bindings("outer_buffer","slope_mu",{
            "f":"c.exp(-t/2-mu*kernels['J'])","u":"u1*f"}),
        power=assignment_source_bindings("outer_buffer","power",{
            "slope":"-c.mpf('.5')-mu","f":"c.exp(slope*t)","u":"u1*f"}))
    result["same_scalar_kernel_callables"]=dict(
        slope=pre_source.transition_integrals is initial_source.transition_integrals,
        turnoff=pre_source.turnoff_kernels is initial_source.turnoff_kernels,
        transition=pre_source.transition_kernels is buffer_source.transition_kernels,
        decay=pre_source.decay_integral is buffer_source.decay_integral)
    if not all(result["same_scalar_kernel_callables"].values()):
        raise ValueError("Different scalar source kernels")
    return result


def terminal_shape_proof():
    """Derive all inlet coefficient functions from actual current hist AST."""
    z=Z;q=1+z*z
    mu,Tw,td,J,I0,I1,I2,KB,KB2,JT,KT,KE,KP,D2,DP=s.symbols(
        "mu Tw td J I0 I1 I2 KB KB2 JT KT KE KP D2 DP",real=True)
    invP2=s.Symbol("invP2",positive=True)
    env=dict(z=z,qi=1/q,y=s.Integer(1),J=J,V=4*z,
        **{"self.invP2":invP2,"mass[0]":I0,"mass[1]":I1,"mass[2]":I2})
    factor=assignment("compliant_pre_pulse_mixed_C4","slope","factor",env)
    env["factor"]=factor
    u=assignment("compliant_pre_pulse_mixed_C4","slope","u",env)
    env["h"]=assignment("compliant_pre_pulse_mixed_C4","slope","h",env)
    history=source_histories("slope",env)
    definitions={}
    def shape(stage,u,history):
        coefficients=dict(U=s.cancel(u*q),M=s.cancel(history["m"]/z),
            H=s.cancel(history["h"]*q),K=s.cancel(history["k"]*q/z),
            P_in=s.cancel(history["p"]*q*q))
        poly=s.Poly(s.cancel(history["e"]*q*q),z)
        coefficients.update(E_Z=poly.coeff_monomial(z**6),E_Q=poly.coeff_monomial(1))
        expected=dict(m=coefficients["M"]*z,h=coefficients["H"]/q,
            k=coefficients["K"]*z/q,e=coefficients["E_Z"]*z*z+coefficients["E_Q"]/q**2,
            p=coefficients["P_in"]/q**2)
        for key,value in coefficients.items():
            if z in value.free_symbols:raise ArithmeticError("Nonconstant inlet source coefficient: "+stage+"/"+key)
        for key in KEYS:
            if s.cancel(history[key]-expected[key])!=0:
                raise ArithmeticError("Actual current q-shape failed: "+stage+"/"+key)
        definitions[stage]={key:str(value) for key,value in coefficients.items()}
        return coefficients
    shape("O2_slope_exit",u,history)
    for stage,method,t in (("Rd","axial",td),("Rw","slope_mu",s.Integer(1)),("Rp","power",Tw)):
        base=dict(z=z,u1=u,t=t,d=s.exp(-t),root=s.exp(-t/2),d3=s.exp(-3*t/2),
            mu=mu,**{"self.invP2":invP2,**{"old['"+key+"']":v for key,v in history.items()},
            "K['B_mass']":KB,"K['B_squared_mass']":KB2,
            "K['theta']":KT,"K['energy']":KE,"K['pressure']":KP,"K['J']":JT,
            "decay_integral(c, 2 * mu, t)":D2,
            "decay_integral(c, 1 + 2 * mu, t)":DP})
        if method=="axial":base["root"]=s.exp(-t/2);new_u=u*base["root"]
        elif method=="slope_mu":
            base["factor"]=assignment("compliant_pre_pulse_mixed_C4",method,"factor",base)
            new_u=assignment("compliant_pre_pulse_mixed_C4",method,"u",base)
        else:
            base["f"]=assignment("compliant_pre_pulse_mixed_C4",method,"f",base)
            base["theta"]=(base["f"]-base["d3"])/(1-mu)
            new_u=assignment("compliant_pre_pulse_mixed_C4",method,"u",base)
        history=source_histories(method,base);u=new_u
        coefficients=shape(stage,u,history)
    expectedEZ=16*invP2*(s.exp(-td)+KB2)*s.exp(-1-Tw)
    if s.simplify(s.expand_power_exp(coefficients["E_Z"]-expectedEZ))!=0:
        raise ArithmeticError("Actual current terminal E_Z differs from native high source")
    qU,qM,qH,qK,qEZ,qEQ,qP=s.symbols("U M H K E_Z E_Q P_in",real=True)
    Pstar=s.Symbol("Pstar",positive=True)
    u=qU/q;m=qM*z;h=qH/q;k=qK*z/q;e=qEZ*z*z+qEQ/q**2;p=qP/q**2
    canon=dict(m1=m/(Pstar*u),m2=k/(Pstar*u*u),X=h/u,E=e/(u*u),Mp=p)
    expected=dict(m1=qM/(Pstar*qU)*z*q,m2=qK/(Pstar*qU*qU)*z*q,
        X=qH/qU,E=qEZ/(qU*qU)*z*z*q*q+qEQ/(qU*qU),Mp=qP/q**2)
    if any(s.cancel(canon[key]-expected[key])!=0 for key in canon):
        raise ArithmeticError("Canonical pulse units differ")
    return dict(actual_current_source_coefficients_by_stage=definitions,
        source_shape_identities=4*7,all_five_histories_retained=True,
        actual_positive_E_Z_recipe_identity=True,canonical_unit_identities=5,
        canonical_source_functions={key:str(value) for key,value in expected.items()},
        same_functions_for_arbitrary_Z_not_saved_sample_interpolation=True,
        unclamped_theta_defines_positive_source_for_0_lt_mu_lt_1_Tw_gt_0=True,passed=True)


@lru_cache(None)
def symbolic_coefficient(expr,n):
    return s.cancel(s.diff(expr,Z,n)/math.factorial(n))


class FunctionJet:
    """Exact ordinary Taylor projection of a symbolic defining Z function."""
    def __init__(self,ctx,coefficients):
        self.ctx=ctx;self.expr=s.cancel(coefficients[0]);self.order=len(coefficients)-1
        for n,value in enumerate(coefficients):
            if s.cancel(symbolic_coefficient(self.expr,n)-value)!=0:
                raise ArithmeticError("Symbolic source projection input is not the defining jet")
    @classmethod
    def function(cls,ctx,expr,order=5):
        obj=object.__new__(cls);obj.ctx=ctx;obj.expr=s.cancel(expr);obj.order=order;return obj
    @classmethod
    def variable(cls,ctx,value,order=5):return cls.function(ctx,value,order)
    @classmethod
    def constant(cls,ctx,value,order=5):return cls.function(ctx,value,order)
    def __getitem__(self,n):
        if n>self.order:raise IndexError("Symbolic projection order exceeded")
        return symbolic_coefficient(self.expr,n)
    @property
    def coefficients(self):return [self[n] for n in range(self.order+1)]
    def truncate(self,n):return self.function(self.ctx,self.expr,min(n,self.order))
    def op(self,other,fn):
        return self.function(self.ctx,fn(self.expr,other.expr if isinstance(other,FunctionJet) else other),
            min(self.order,other.order) if isinstance(other,FunctionJet) else self.order)
    def __add__(self,v):return self.op(v,lambda a,b:a+b)
    __radd__=__add__
    def __sub__(self,v):return self.op(v,lambda a,b:a-b)
    def __rsub__(self,v):return self.op(v,lambda a,b:b-a)
    def __mul__(self,v):return self.op(v,lambda a,b:a*b)
    __rmul__=__mul__
    def __truediv__(self,v):return self.op(v,lambda a,b:a/b)
    def reciprocal(self):return self.function(self.ctx,1/self.expr,self.order)
    def __neg__(self):return self.function(self.ctx,-self.expr,self.order)


def replay(stem,name,env,hashes):
    path=HERE/(PREFIX+stem+".py");hashes[path.name]=sha(path.name)
    fn=next(v for v in ast.walk(ast.parse(path.read_text(encoding="utf8")))
        if isinstance(v,ast.FunctionDef) and v.name==name)
    fn=copy.deepcopy(fn);fn.decorator_list=[]
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        str(path),"exec"),env)
    return env[name]


def mixed_source_join_proof():
    """Replay BOTH physical derivative algorithms in exact source functions."""
    z=Z;q=1+z*z
    U,M,H,K,EZ,EQ,Pin=s.symbols("U M H K EZ EQ Pin",real=True)
    Pstar=s.Symbol("Pstar",positive=True);mu,delta=s.symbols("mu delta",real=True)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)) if isinstance(v,(str,int,float)) else v)
    J=lambda expr,order=5:FunctionJet.function(c,expr,order)
    env=dict(math=math,IntervalTaylor=FunctionJet,
        square=lambda v:v*v,derivative=lambda v:J(s.diff(v.expr,z),v.order-1))
    hashes={}
    for stem,name in (("pre_pulse_mixed_C4","rate_rows"),("pre_pulse_mixed_C4","product_rows"),
            ("long_reshape_mixed_C4","exponential_derivatives"),
            ("pulse_mixed_C4","binomial_product")):
        replay(stem,name,env,hashes)
    leftfn=replay("pre_pulse_mixed_C4","physical_mixed",env,hashes)
    radial=replay("pulse_high_jets","radial",env,hashes)
    rightfn=replay("pulse_mixed_C4","transport_mixed",env,hashes)
    zero=J(0);u=J(U/q);p0=J(s.Function("same_P0")(z))
    histories=dict(m=J(M*z),h=J(H/q),k=J(K*z/q),
        e=J(EZ*z*z+EQ/q**2),p=J(Pin/q**2))
    left=leftfn(c,z,delta,u,[zero-(s.Rational(1,2)+mu)]+[zero]*3,
        [zero]*5,histories,p0,1/Pstar**2)
    field=SimpleNamespace(ctx=c,mu=mu,delta=delta,prate=1+2*mu)
    field.radial=lambda zz,b,m:radial(field,zz,b,m)
    point=dict(Mz_over_R_Utheta=histories["m"]/(u*Pstar),
        Mtheta_z_over_sqrt2_R_3half_Utheta_squared=histories["k"]/(u*u*Pstar),
        Mtheta_over_sqrt2_R_3half_Utheta=histories["h"]/u,
        Mztheta_over_R_Utheta_squared=histories["e"]/(u*u),
        pressure=dict(P_over_Pstar_squared=p0+histories["p"],P_y_over_Pstar_squared=u*u/2))
    right=rightfn(field,z,point,[zero]*5,u)
    mapping=dict(Utheta_over_Pstar=("Utheta_over_Pstar_without_common_theta_radial_factor",1),
        Uz=("Uz_over_Pstar_without_common_theta_radial_factor",Pstar),
        Ur_over_current_sqrt_R_over_2=("Ur_over_sqrt_R_over_2_Pstar_without_common_theta_radial_factor",Pstar),
        P_over_Pstar2=("P_over_Pstar_squared",1))
    count=0
    for label,(other,factor) in mapping.items():
        for index,value in left["physical_velocity_pressure_y_Z_mixed4"][label].items():
            if s.cancel(value-factor*right["physical_mixed_derivatives_total_order_le4"][other][index])!=0:
                raise ArithmeticError("Actual Rp physical source row differs: "+label+"/"+index)
            count+=1
    rows=right["primitive_y_derivative_Taylor"];binomial=env["binomial_product"]
    primitive=dict(
        Mz_over_current_R=(rows["Mz_over_R_Utheta"],s.Rational(1,2)-mu,u*Pstar),
        Mtheta_over_current_sqrt2_R_1p5_Pstar=(rows["Mtheta_over_sqrt2_R_3half_Utheta"],1-mu,u),
        Mtheta_z_over_current_sqrt2_R_1p5_Pstar=(rows["Mtheta_z_over_sqrt2_R_3half_Utheta_squared"],s.Rational(1,2)-2*mu,u*u*Pstar),
        Mztheta_over_current_R_Pstar2=(rows["Mztheta_over_R_Utheta_squared"],-2*mu,u*u))
    for label,(values,rate,factor) in primitive.items():
        for k in range(5):
            transformed=factor*binomial(values,rate,k)
            for n in range(5-k):
                index="y"+str(k)+"_Z"+str(n)
                if s.cancel(left["physical_five_primitive_y_Z_mixed4"][label][index]-transformed[n]*math.factorial(n))!=0:
                    raise ArithmeticError("Actual Rp primitive physical units differ: "+label+"/"+index)
                count+=1
    for k in range(5):
        value=rows["P_over_Pstar_squared"][k]-(p0 if k==0 else zero)
        for n in range(5-k):
            index="y"+str(k)+"_Z"+str(n)
            if s.cancel(left["physical_five_primitive_y_Z_mixed4"]["Mp_over_Pstar2"][index]-value[n]*math.factorial(n))!=0:
                raise ArithmeticError("Actual Rp pressure primitive differs")
            count+=1
    return dict(exact_replayed_velocity_pressure_and_primitive_mixed4_rows=count,
        full_physical_radial_prefactors_and_Pstar_units_retained=True,
        both_actual_source_algorithms_replayed=True,ordinary_axial_projection_through_five=True,
        common_arbitrary_pressure_function_retained=True,passed=True,input_hashes=hashes)


def selected_node(stem,method,target):
    tree=ast.parse((HERE/(PREFIX+stem+".py")).read_text(encoding="utf8"))
    fn=next(v for v in ast.walk(tree) if isinstance(v,ast.FunctionDef) and v.name==method)
    nodes=[v.value for v in ast.walk(fn) if isinstance(v,ast.Assign)
        and any(ast.unparse(t)==target for t in v.targets)]
    if len(nodes)!=1:raise ValueError("Unique actual source projection assignment required")
    return nodes[0]


def native_inlet_projection_proof():
    """Replay actual C4 incoming, fifth projection and native data dictionary."""
    z=Z;q=1+z*z
    U,M,H,K,EZ,EQ,Pin=s.symbols("U M H K EZ EQ Pin",real=True)
    Pstar=s.Symbol("Pstar",positive=True)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)) if isinstance(v,(str,int,float)) else v)
    constants=dict(U=U,M=M,K=K,E_Z=EZ,E_Q=EQ,C1=M/(Pstar*U),
        C2=K/(Pstar*U*U),C0=EQ/(U*U),C_E=EZ/(U*U))
    env=dict(IntervalTaylor=FunctionJet,mp=mp,endpoints=lambda value:(-1,1),
        copy_jet=lambda ctx,jet:jet)
    hashes={}
    incomingfn=replay("axial_high_jets","incoming",env,hashes)
    old=incomingfn(SimpleNamespace(ctx=c,constants=constants),z)
    append=replay("fifth_axial_jets","append_fifth",env,hashes)
    env.update(c=c,incoming=old,Z=z,k=constants,append_fifth=append)
    project=lambda target:eval(compile(ast.Expression(copy.deepcopy(
        selected_node("fifth_axial_jets","select",target))),"<actual projection>","eval"),env)
    moments=project("moments");energy=project("energy")
    env.update(moments=moments,energy=energy)
    tree=ast.parse((HERE/(PREFIX+"fifth_axial_jets.py")).read_text(encoding="utf8"))
    fn=next(v for v in ast.walk(tree) if isinstance(v,ast.FunctionDef) and v.name=="select")
    updates=[kw.value for v in ast.walk(fn) if isinstance(v,ast.Call)
        and ast.unparse(v.func)=="result.update" for kw in v.keywords if kw.arg=="incoming"]
    if len(updates)!=1:raise ValueError("Actual fifth incoming publication changed")
    published=eval(compile(ast.Expression(copy.deepcopy(updates[0])),"<actual incoming publication>","eval"),env)
    selected=dict(incoming=published)
    from lei_ren_part1_paper_candidate_pressure_function import q_jets
    env["q_jets"]=q_jets
    hashes["lei_ren_part1_paper_candidate_pressure_function.py"]=sha("lei_ren_part1_paper_candidate_pressure_function.py")
    pulse=SimpleNamespace(ctx=c,high=SimpleNamespace(constants=constants),
        fifth=SimpleNamespace(select=lambda Z:selected),data_cache={},inlet_H=H,inlet_P=Pin)
    _,inlet,u,energy,moments=replay("pulse_radial_C4","data",env,hashes)(pulse,z)
    expected=dict(Utheta_over_Pstar=U/q,Mz_over_R=M*z,
        Mtheta_over_sqrt2_R_3half_Pstar=H/q,Mtheta_z_over_sqrt2_R_3half_Pstar=K*z/q,
        Mztheta_over_R_Pstar_squared=EZ*z*z+EQ/q**2,Mp_over_Pstar_squared=Pin/q**2)
    count=0
    for key,expr in expected.items():
        for n,value in enumerate(inlet[key]):
            if s.cancel(value-symbolic_coefficient(expr,n))!=0:
                raise ArithmeticError("Actual native inlet projection differs: "+key+"/"+str(n))
            count+=1
    canonical=[M/(Pstar*U)*z*q,K/(Pstar*U**2)*z*q,EQ/U**2+EZ/U**2*z*z*q*q]
    for jet,expr in zip(moments+[energy],canonical):
        for n in range(6):
            if s.cancel(jet[n]-symbolic_coefficient(expr,n))!=0:
                raise ArithmeticError("Actual fifth incoming canonical function differs")
            count+=1
    return dict(actual_native_inlet_axial5_coefficient_identities=count,
        original_high_incoming_fifth_projection_publication_and_native_data_AST_replayed=True,
        native_H_and_P_come_from_live_power_callable_not_saved_samples=True,
        positive_E_Z_kept_separate_from_negative_E_Q=True,passed=True,input_hashes=hashes)


def native_source_graph(pre,pulse):
    initial=pulse.pulse.initial;future=pulse.selection.future.angular.initial
    graph=dict(native_buffer_class=type(pulse.pulse.buffer) is buffer_source.SharedOuterBuffer,
        current_buffer_class=type(pre.buffer) is buffer_source.SharedOuterBuffer,
        native_initial_object=pulse.pulse.initial is pulse.pulse.buffer.initial,
        current_initial_object=pre.initial is pre.buffer.initial,
        live_constants_object=pulse.high.constants is pulse.fifth.fourth.constants,
        live_high_pulse_object=pulse.pulse is pulse.fifth.fourth.base.pulse,
        native_source_family_matches=(initial.family,future.family)==(pre.family,pre.family),
        same_analytic_pressure_definitions=initial.datum.definition==future.datum.definition==pre.datum.definition,
        same_analytic_pressure_source_and_datum=(
            initial.datum.source_sha,initial.datum.datum_sha,future.datum.source_sha,future.datum.datum_sha)==(
            pre.datum.source_sha,pre.datum.datum_sha,pre.datum.source_sha,pre.datum.datum_sha),
        same_pressure_callable=(type(initial.datum).normalized_jets is
            type(future.datum).normalized_jets is type(pre.datum).normalized_jets),
        native_parameter_objects=initial.params is initial.datum.parameters and future.params is future.datum.parameters,
        same_Cstar_defining_records=(
            initial.repair.records["compliant_physical_norm_family"]==
            future.repair.records["compliant_physical_norm_family"]==
            pre.initial.repair.records["compliant_physical_norm_family"]),
        native_future_parameter_object=pulse.high.base.future.params is future.params,
        native_pulse_parameter_object=pulse.pulse.buffer.params is initial.params,
        current_positive_mu_Tw=0<endpoints(pre.params.mu)[0]<1 and endpoints(pre.params.Tw)[0]>0,
        actual_original_gp_callable=flat_source.original_gp is axial_source.gp,
        actual_native_gp_callable=axial_source.CompliantAxialPulseField.main.__globals__["gp"] is axial_source.gp,
        actual_native_gp_jet_callable=CompliantPulseMixedC4._high_packet.__globals__["gp_jets"] is flat_source.gp_jets,
        native_entrance_callable=CompliantPulseMixedC4.entrance is axial_source.CompliantAxialPulseField.entrance)
    if not all(graph.values()):raise ValueError("Actual native pulse source graph differs: "+str(graph))
    return graph


def entrance_boundary_proof(c):
    gp=axial_source.gp(c,c.mpf(0));jets=flat_source.gp_jets(c,c.mpf(0))
    if (any(endpoints(gp[key])!=(mp.mpf(0),mp.mpf(0)) for key in ("value","derivative"))
            or any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in jets)
            or endpoints(axial_source.gp_energy(c,c.mpf(0),None))!=(mp.mpf(0),mp.mpf(0))):
        raise ArithmeticError("Actual pulse entrance shape is not flat")
    bindings=dict(
        current_terminal=assignment_source_bindings("pre_pulse_mixed_C4","power",{"t":"self.params.Tw*phase"}),
        native_incoming_moment=assignment_source_bindings("axial_pulse_field","main",{
            "incoming":"rows[row-1]*self.factor(-lam*t,cut)",
            "X":"1/self.rate+(self.Xp-1/self.rate)*self.factor(-self.rate*t)",
            "emu":"(e0*self.mu+ap*ap*K-decay_integral(c,2,xi)/2)*c.exp(2*xi)"}),
        radius=assignment_source_bindings("pulse_physical_bounds","__init__",{
            "self.logP":"c.exp(c.mpf(parameters['Md']))+11",
            "self.logmu":"c.ln(c.mpf(parameters['c_mu']))-4*self.logP",
            "self.Tw":"-60*self.logmu",
            "self.logRp_parts":"dict(logCstar=10*self.logC,logPstar=10*self.logP,finite_outer_offset=c.ln(110)+self.logP+1+self.Tw)"}))
    tree=ast.parse((HERE/(PREFIX+"axial_pulse_field.py")).read_text(encoding="utf8"))
    fn=next(v for v in ast.walk(tree) if isinstance(v,ast.FunctionDef) and v.name=="entrance")
    result=next(v.value for v in ast.walk(fn) if isinstance(v,ast.Return))
    if ast.dump(result)!=ast.dump(ast.parse("self.main(Z,self.mu*t,entrance_t=t)",mode="eval").body):
        raise ValueError("Native entrance-to-main coordinate changed")
    return dict(actual_gp_value_and_ordinary_jets0_through4_exact_zero=True,
        actual_partial_energy_empty_at_Rp=True,source_calls_and_radius_AST=bindings,
        exact_common_Rp_definition="logRp=ln110+10(logCstar+logPstar)+logPstar+1+Tw",
        original_Rp_relative_pulse_coordinate_zero=True,passed=True)


def actual_parameter_formula_bindings():
    name="lei_ren_part1_paper_logarithmic_outer_parameters.py"
    tree=ast.parse((HERE/name).read_text(encoding="utf8"))
    cls=next(v for v in tree.body if isinstance(v,ast.ClassDef) and v.name=="LogarithmicOuterParameters")
    fn=next(v for v in cls.body if isinstance(v,ast.FunctionDef) and v.name=="__init__")
    wanted=dict(self_logPstar="c.exp(self.md)+11",self_yd="c.exp(self.md)+11",
        self_log_mu="c.ln(c.mpf('.001'))-4*self.logPstar",
        self_mu="c.exp(self.log_mu)",self_Tw="-60*self.log_mu")
    for key,expr in wanted.items():
        target=key.replace("self_","self.",1)
        nodes=[v.value for v in ast.walk(fn) if isinstance(v,ast.Assign)
            and any(ast.unparse(t)==target for t in v.targets)]
        expected=ast.dump(ast.parse(expr,mode="eval").body)
        if sum(ast.dump(v)==expected for v in nodes)!=1:
            raise ValueError("Actual parameter defining formula changed: "+target)
    return dict(same_yd_logPstar_exact_definition=True,positive_mu_and_Tw_source_bound=True,
        original_parameter_assignments_checked=len(wanted),
        input_hashes={name:sha(name)},passed=True)


class CurrentRpPulseSourceDispatcher(CurrentPrePulseSourceDispatcher):
    CHARTS=CurrentPrePulseSourceDispatcher.CHARTS+("pulse_entrance",)

    def __init__(self,require_checked=True):
        super().__init__(require_checked=True)
        self.require_Rp_checked=require_checked;self.Rp_loaded=False;self.Rp_acceptance_loaded=False

    @source_precision
    def _load_Rp(self):
        if self.Rp_loaded:return
        super()._load_pre()
        self.native_pulse=CompliantPulseMixedC4()
        self.Rp_graph=native_source_graph(self.rh_reference,self.native_pulse)
        from lei_ren_part1_paper_compliant_source_dispatcher import ROUTES
        if ROUTES["pulse_entrance"][:3]!=("pulse_mixed_C4","CompliantPulseMixedC4","entrance"):
            raise ValueError("Actual native dispatcher entrance route changed")
        self.Rp_graph["native_dispatcher_route_bound"]=True
        name=PREFIX+"pulse_mixed_C4_check.json"
        receipt=accepted(name,self.family,self.source,"factored_all_chart_mixed_derivative_boxes_available")
        for path,digest in {**self.native_pulse.hashes,**receipt["input_hashes"]}.items():
            if path in self.hashes and self.hashes[path]!=digest:
                raise ValueError("Current Rp native source dependency conflict: "+path)
            self.hashes[path]=digest
        self.hashes[name]=sha(name)
        name=PREFIX+"fifth_axial_jets_check.json"
        fifth=accepted(name,self.family,self.source,"actual_selected_ap_c1_c2_C5_available")
        self.hashes.update(fifth["input_hashes"]);self.hashes[name]=sha(name)
        self.Rp_recipe=source_recipe_bindings();self.Rp_shapes=terminal_shape_proof()
        self.Rp_projection=native_inlet_projection_proof();self.Rp_mixed=mixed_source_join_proof()
        self.Rp_boundary=entrance_boundary_proof(self.rh_reference.ctx)
        self.Rp_parameters=actual_parameter_formula_bindings()
        for proof in (self.Rp_projection,self.Rp_mixed,self.Rp_parameters):self.hashes.update(proof["input_hashes"])
        for stem in ("actual_Rp_source_join","pulse_physical_bounds","axial_high_jets",
                "fifth_axial_jets","pulse_radial_C4","flat_pulse_derivatives","axial_pulse_field"):
            name=PREFIX+stem+".py";self.hashes[name]=sha(name)
        if self.require_Rp_checked:
            name=PREFIX+"actual_Rp_source_join_check.json"
            receipt=accepted(name,self.family,self.source,"current_Rp_external_pulse_join_certified")
            if receipt["datum_enclosure_sha256"]!=self.datum_sha:
                raise ValueError("Current Rp acceptance pressure datum differs")
            self.hashes.update(receipt["input_hashes"]);self.hashes[name]=sha(name)
            self.Rp_acceptance_loaded=True
        self.Rp_loaded=True

    def provider(self,chart):
        if chart=="pulse_entrance":self._load_Rp();return self.native_pulse
        return super().provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart!="pulse_entrance":return super().evaluate(chart,Z,coordinate)
        field=self.provider(chart);t=field.ctx.mpf(coordinate)
        # The .02/mu coordinate box encloses the exact source endpoint.
        # Its interval image has the same tiny outward parameter error.
        limit=endpoints(field.ctx.mpf(".02")/field.mu*field.mu)[1]
        if endpoints(t)[0]<0 or endpoints(field.mu*t)[1]>limit:
            raise ValueError("Current native entrance coverage requires 0<=mu*t<=.02")
        packet=field.entrance(Z,t)
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,
            implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            source_provider=PREFIX+"pulse_mixed_C4.CompliantPulseMixedC4",
            acceptance_receipt=PREFIX+"actual_Rp_source_join_check.json",
            source_coverage_coordinate="t=log(R/Rp); xi=mu*t",source_coordinate_domain="[0,.02/mu]",
            derivative_coordinate="ordinary logR,Z; full radial velocity factors differentiated",
            physical_mixed_grids={"physical_mixed_derivatives_total_order_le4":
                packet["physical_mixed_derivatives_total_order_le4"]},source_packet=packet,
            original_scale_metadata=packet["formal_log_Utheta_over_Pstar"],
            current_Rp_external_pulse_join_proved=True,
            current_Rp_external_pulse_join_certified=self.Rp_acceptance_loaded,
            output_kind="current native pulse source enclosures; no production point selection",
            **dict.fromkeys(SCOPES,False))

    def manifest(self):
        result=super().manifest();self._load_Rp()
        result["ordered_current_chart_registry"]["pulse_entrance"]=dict(
            provider=PREFIX+"pulse_mixed_C4.CompliantPulseMixedC4",method="entrance",
            coverage_coordinate="t=log(R/Rp); xi=mu*t",domain="[0,.02/mu]",
            acceptance_receipt=PREFIX+"actual_Rp_source_join_check.json")
        result.update(current_Rp_native_source_graph=self.Rp_graph,
            current_Rp_scalar_recipe_bindings=self.Rp_recipe,current_Rp_terminal_shape_proof=self.Rp_shapes,
            current_Rp_native_inlet_projection_proof=self.Rp_projection,
            current_Rp_exact_physical_mixed_source_join=self.Rp_mixed,
            current_Rp_entrance_boundary_and_radius_proof=self.Rp_boundary,
            current_Rp_parameter_source_proof=self.Rp_parameters,
            current_Rp_external_pulse_join_proved=True,
            current_Rp_external_pulse_join_certified=self.Rp_acceptance_loaded,
            current_downstream_chart_owner_count=15,current_native_pulse_entrance_owner_installed=True,
            all_current_pulse_charts_installed=False,all_profile_source_charts_callable=False,
            full_current_core_to_heat_physical_assembly=False,
            **dict.fromkeys(SCOPES,False),input_hashes=dict(self.hashes))
        return result


@source_precision
def build():
    field=CurrentRpPulseSourceDispatcher(require_checked=False)
    result=field.manifest();c=field.native_pulse.ctx
    result["current_Rp_left_source"]=field.evaluate("O3_power",[-1,1],1)
    result["current_Rp_native_right_source"]=field.evaluate("pulse_entrance",[-1,1],0)
    upper=endpoints(c.mpf(".02")/field.native_pulse.mu)[1]
    result["whole_current_native_entrance"]=field.evaluate("pulse_entrance",[-1,1],[0,upper])
    result["whole_entrance_original_domain_fully_enclosed"]=True
    result["input_hashes"]=dict(field.hashes)
    return encode(pack(result))


def run():
    result=build()
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("Current Rp source join generated: native live inlet, 135 mixed identities, 15 current owners",flush=True)
    return result


if __name__=="__main__":
    run()

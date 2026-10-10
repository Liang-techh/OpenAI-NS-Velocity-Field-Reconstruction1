"""Exact original leading kernels, common C3 units and outer source interface.

Scalar archive covers remain bounds. The functions below use the original
defining integrals, source parameters and parent histories, never cap values.
"""
import ast
import copy
import gzip
import inspect
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_outer_to_repair_C3_bridge as bridge

quiet,outer,band,target,current,source=bridge.quiet,bridge.outer,bridge.band,bridge.target,bridge.current,bridge.source
providers=band.limit.outer
HERE,PREFIX,sha=bridge.HERE,bridge.PREFIX,bridge.sha
NAME=PREFIX+'current_original_outer_leading_C3_identity.json.gz'
RECEIPT=PREFIX+'current_original_outer_leading_C3_identity_check.json'
GATES=('current_original_outer_leading_C3_source_ODE_and_seam_functions_installed',)
OPEN=bridge.OPEN


def integral(g,body,variable,upper,lower=None):
    lower=g.zero if lower is None else lower
    return current.controls.integral(g,body,variable,lower,upper,
        measure='exact original mathematical kernel, not scalar enclosure output',
        original_kernel_function=True,coordinate_and_kernel_Z_independent=True)


def sigma(g,x):
    if x==g.zero:return g.zero
    if x==g.one:return g.one
    return g.node('exact_original_source_flat_sigma',argument=x.node,
        definition='0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside',
        source_binding=current.ast_binding(providers.slope.original.original.sigma_jets),
        smooth_flat_endpoints=True,Z_independent_argument=True)


def at(g,q,variable,value):return outer.substitute(g,q,variable,value)


def scalar_at(g,q,variable,value):
    return g.node('function_substitution',expression=q.node,variable=variable.node,value=value.node,
        Z_independent_substitution=True)


def primitive_sigma(g,end,name):
    u=g.symbol(name+'_sigma_s');return integral(g,sigma(g,u),name+'_sigma_s',end)


def weighted_sigma_kernel(g,end,name,linear,multiplier):
    u=g.symbol(name+'_mass_s');J=primitive_sigma(g,u,name+'_inner')
    body=g.unary('exp',g.sub(g.mul(g.constant(linear),u),g.mul(multiplier,J)))
    return integral(g,body,name+'_mass_s',end)


def decay_mass(g,k,t):return g.mul(t,g.unary('exprel',g.neg(g.mul(k,t))))


def own_first_y(g,E,V,H):
    a=target.C3Algebra(g);EE=a.mul(E,E)
    D=dict(m=V,h=E,k=a.mul(E,V),e=a.add(a.mul(V,V),a.scale(EE,'-1/2')),p=a.scale(EE,'1/2'))
    return {key:a.add(D[key],a.scale(H[key],-rate)) for key,rate in current.RATES.items()}


def build_source(field):
    g=field.graph;a=target.C3Algebra(g);z=target.C3Function(g.symbol('Z'),g.one,g.zero,g.zero)
    one=a.fixed(1);A=a.div(one,a.add(one,a.mul(z,z)));params=field.parameters
    invS=g.unary('exp',g.neg(params['logP']));B=a.scale(z,g.mul(g.constant(4),invS));zero=a.fixed(0)
    BB,AA=a.mul(B,B),a.mul(A,A)
    charts={};kernel_refs={}
    def packet(chart,x,E,V,H,seeds,kernels=None,physical_coordinate=None):
        return dict(chart=chart,native_coordinate=x,physical_coordinate=x if physical_coordinate is None else physical_coordinate,
            E=E,V=V,histories=H,history_y=own_first_y(g,E,V,H),
            ordinary_Z_orders=[0,1,2,3],seeds=seeds,kernels={} if kernels is None else kernels)
    y=g.symbol('native_Rh_reference');E=a.scale(A,g.unary('exp',g.mul(g.constant('1/10'),y)))
    h=a.scale(E,'5/8');EE=a.mul(E,E)
    H=dict(m=B,h=h,k=a.mul(h,B),e=a.add(BB,a.scale(EE,'-5/12')),p=a.scale(EE,'5/2'))
    charts['Rh_reference']=packet('Rh_reference',y,E,B,H,dict(A=A,B=B))
    x=g.symbol('native_O2_slope');J=primitive_sigma(g,x,'original_O2_slope_J')
    factor=g.unary('exp',g.sub(g.mul(g.constant('1/10'),x),g.mul(g.constant('3/5'),J)))
    Kt=weighted_sigma_kernel(g,x,'original_O2_theta','8/5',g.constant('3/5'))
    Ke=weighted_sigma_kernel(g,x,'original_O2_energy','6/5',g.constant('6/5'))
    Kp=weighted_sigma_kernel(g,x,'original_O2_pressure','1/5',g.constant('6/5'))
    d=g.unary('exp',g.neg(x));d3=g.unary('exp',g.mul(g.constant('-3/2'),x))
    h=a.scale(A,g.mul(g.add(g.constant('5/8'),Kt),d3))
    H=dict(m=B,h=h,k=a.mul(h,B),e=a.add(BB,a.scale(AA,g.mul(g.constant(-1),g.add(g.constant('5/12'),g.mul(g.constant('1/2'),Ke)),d))),
        p=a.scale(AA,g.add(g.constant('5/2'),g.mul(g.constant('1/2'),Kp))))
    charts['O2_slope']=packet('O2_slope',x,a.scale(A,factor),B,H,dict(A=A,B=B),dict(J=J,theta=Kt,energy=Ke,pressure=Kp))
    A1=at(g,charts['O2_slope']['E'],x,g.one)
    H1={key:at(g,q,x,g.one) for key,q in H.items()};A1A1=a.mul(A1,A1)
    u=g.symbol('original_O2_turnoff_integrand_y');yc=g.symbol('original_O2_turnoff_endpoint_y')
    beta=sigma(g,g.sub(g.one,g.quotient(g.unary('log',u),g.constant(40),'original positive Md=40')))
    weight=g.unary('exp',g.sub(u,yc))
    KB=integral(g,g.mul(weight,beta),'original_O2_turnoff_integrand_y',yc,g.one)
    KB2=integral(g,g.mul(weight,beta,beta),'original_O2_turnoff_integrand_y',yc,g.one)
    Y=g.unary('exp',g.constant(40))
    def axial_chart(chart):
        x=g.symbol('native_'+chart);y=g.unary('exp',g.mul(g.constant(40),x)) if chart=='O2_axial' else g.add(Y,x)
        t=g.sub(y,g.one);d=g.unary('exp',g.neg(t));root=g.unary('exp',g.mul(g.constant('-1/2'),t));d3=g.unary('exp',g.mul(g.constant('-3/2'),t))
        if chart=='O2_axial':
            kb,kb2=(scalar_at(g,q,yc,y) for q in (KB,KB2))
            beta_y=sigma(g,g.sub(g.one,g.quotient(g.unary('log',y),g.constant(40),'original positive Md=40')))
            V=a.scale(B,beta_y)
        else:
            suffix=g.unary('exp',g.neg(x))
            kb,kb2=(g.mul(suffix,scalar_at(g,q,yc,Y)) for q in (KB,KB2));V=zero
        H=dict(m=a.add(a.scale(H1['m'],d),a.scale(B,kb)),
            h=a.add(a.scale(H1['h'],d3),a.scale(A1,g.sub(root,d3))),
            k=a.add(a.scale(H1['k'],d3),a.scale(a.mul(A1,B),g.mul(root,kb))),
            e=a.add(a.scale(H1['e'],d),a.scale(BB,kb2),a.scale(A1A1,g.mul(g.constant('-1/2'),t,d))),
            p=a.add(H1['p'],a.scale(A1A1,g.mul(g.constant('1/2'),g.sub(g.one,d)))))
        return packet(chart,x,a.scale(A1,root),V,H,dict(A=A1,B=B,H=H1),
            dict(B_mass=kb,B_squared_mass=kb2,cutoff=V,original_y_endpoint=yc,
                original_B_mass=KB,original_B_squared_mass=KB2),y)
    for chart in ('O2_axial','O2_buffer'):charts[chart]=axial_chart(chart)
    buf=charts['O2_buffer'];Ad=at(g,buf['E'],buf['native_coordinate'],g.constant(11))
    Hd={key:at(g,q,buf['native_coordinate'],g.constant(11)) for key,q in buf['histories'].items()}
    x=g.symbol('native_O3_transition');mu=params['mu'];J=primitive_sigma(g,x,'original_O3_J')
    Kt=weighted_sigma_kernel(g,x,'original_O3_theta',1,mu)
    Ke=weighted_sigma_kernel(g,x,'original_O3_energy',0,g.mul(g.constant(2),mu))
    Kp=weighted_sigma_kernel(g,x,'original_O3_pressure',-1,g.mul(g.constant(2),mu))
    d=g.unary('exp',g.neg(x));d3=g.unary('exp',g.mul(g.constant('-3/2'),x));AdAd=a.mul(Ad,Ad)
    E=a.scale(Ad,g.unary('exp',g.add(g.mul(g.constant('-1/2'),x),g.neg(g.mul(mu,J)))))
    H=dict(m=a.scale(Hd['m'],d),h=a.add(a.scale(Hd['h'],d3),a.scale(Ad,g.mul(d3,Kt))),k=a.scale(Hd['k'],d3),
        e=a.add(a.scale(Hd['e'],d),a.scale(AdAd,g.mul(g.constant('-1/2'),d,Ke))),
        p=a.add(Hd['p'],a.scale(AdAd,g.mul(g.constant('1/2'),Kp))))
    charts['O3_transition']=packet('O3_transition',x,E,zero,H,dict(A=Ad,H=Hd),dict(J=J,theta=Kt,energy=Ke,pressure=Kp))
    Aw=at(g,E,x,g.one);Hw={key:at(g,q,x,g.one) for key,q in H.items()}
    x=g.symbol('native_O3_power');flow=quiet.power_flow(g,Aw,Hw,mu,x)
    charts['O3_power']=packet('O3_power',x,flow['profiles_y0_y1_y2_y3_y4_C3'][0]['E'],zero,
        flow['complete_histories_y0_y1_y2_y3_y4_C3'][0],dict(A=Aw,H=Hw))
    terminal={key:at(g,q,x,g.constant(2)) for key,q in charts['O3_power']['histories'].items()}
    geometry=copy.deepcopy(field.report['actual_geometry_bindings'])
    rh_report=outer.rh.read(outer.rh.NAME)
    geometry['Rh_reference']=rh_report['actual_Rh_C3_continuation_functions']['actual_geometry_source_binding']
    return dict(exact_common_Rref_A_C3=A,exact_common_axial_B_C3=B,original_inverse_Pstar=invS,
        charts=charts,actual_source_Rc_leading_C3=terminal,actual_original_source_geometry=geometry,
        exact_integral_functions_distinct_from_all_archived_caps=True,
        normalized_m_k_and_V_divided_by_Pstar_once=True,common_history_units_no_Rm_multiplier=True)


def source_contracts():
    rh,slope,axial,o3=providers.rh.original,providers.slope.original,providers.axial.original,providers.o3.original
    recovery=providers.o3.reference.long
    return dict(original_rh=rh.source_bindings(),original_slope=slope.source_bindings(),
        original_axial_buffer=axial.source_bindings(),original_o3=o3.source_bindings(),
        pure_frontend_bindings=providers.FRONTEND_BINDINGS,
        actual_background_bindings={key:current.ast_binding(module.background_cell)
            for key,module in (('rh',rh),('slope',slope),('axial_buffer',axial),('o3',o3))},
        unchanged_normalized_recovery=current.ast_binding(recovery.generic.recover_inputs),
        actual_recovery_compiler=current.ast_binding(recovery.compile_recovery),
        actual_recovery_adaptation=recovery.RECOVERY_BINDING,
        original_sigma=current.ast_binding(providers.slope.original.original.sigma_jets),
        original_sigma_left=current.ast_binding(providers.slope.original.original.sigma_jets.__globals__['_sigma_left']),
        original_scalar_sigma=current.ast_binding(providers.axial.original.kernels.stable_sigma),
        original_slope_scalar_sigma=current.ast_binding(providers.slope.original.original.transition_integrals.__globals__['sigma_value_derivative']),
        original_slope_alpha=current.ast_binding(providers.slope.original.original.transition_integrals.__globals__['sigma_value_derivative'].__globals__['alpha_box']),
        original_slope_kernel=current.ast_binding(providers.slope.original.original.slope_masses),
        original_slope_integrals=current.ast_binding(providers.slope.original.original.transition_integrals),
        original_turnoff_kernel=current.ast_binding(providers.axial.original.kernels.turnoff_kernels),
        original_turnoff_derivatives=current.ast_binding(providers.slope.original.original.turnoff_derivatives),
        original_transition_kernel=current.ast_binding(providers.o3.original.kernels.transition_kernels),
        accepted_exact_power_flow=current.ast_binding(quiet.power_flow),
        original_native_geometry=current.ast_binding(current.exact_geometry),
        exact_source_math_builder=current.ast_binding(build_source))


class CurrentOuterLeadingC3Identity:
    def __init__(self,require_checked=True):
        self.bridge=bridge.CurrentOuterToRepairC3Bridge();self.hashes=dict(self.bridge.hashes)
        for module in (bridge,):
            self.hashes[module.NAME]=sha(module.NAME);self.hashes[module.RECEIPT]=sha(module.RECEIPT)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.identity=self.bridge.identity;self.report=self.bridge.outer_report
        self.graph=outer.rh.restore_graph(self.report);self.prefix=copy.deepcopy(self.graph.nodes)
        ref=lambda i:source.FunctionRef(self.graph,i)
        base=outer.rh.read(band.control.NAME)['actual_C3_limit_controls_and_Picard_functions']['exact_C2_limit_functions']['exact_C1_limit_adapter']
        raw=outer.rh.read(current.NAME)['exact_current_integral_and_control_graph']
        self.parameters={key:ref(i) for key,i in raw['parameters'].items()}
        # mu belongs to the shared exact original prefix, not an archive interval.
        self.parameters['mu']=ref(base['parameters']['mu'])
        self.contracts=source_contracts();self.functions=build_source(self);self.acceptance_loaded=False
        if require_checked:
            receipt=outer.rh.read(RECEIPT)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual leading source kernels required')
            for name,value in receipt['input_hashes'].items():
                if sha(name)!=value:raise ValueError('Changed original leading identity input '+name)
            self.acceptance_loaded=True


def run():
    began=time.monotonic();field=CurrentOuterLeadingC3Identity(require_checked=False)
    report=dict(candidate_actual_original_outer_leading_C3_source_kernels_constructed=True,source_family=field.identity,
        exact_graph_nodes=field.graph.nodes,original_outer_graph_prefix_length=len(field.prefix),
        actual_original_leading_C3_source_functions=target.encoded(field.functions),
        original_source_contracts=field.contracts,
        current_outer_leading_endpoint_band_seed_function_identity_installed=False,
        current_outer_complete_history_to_repair_inlet_function_identity_installed=False,
        selected_native_pulse_constructor_consumes_current_C3_frame=False,
        current_numeric_point_field_oracle_installed=False,
        global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(report,separators=(',',':'))+'\n').encode(),mtime=0))
    print('Exact original leading C3 source kernels and all outer charts constructed',flush=True)
    return field


if __name__=='__main__':run()

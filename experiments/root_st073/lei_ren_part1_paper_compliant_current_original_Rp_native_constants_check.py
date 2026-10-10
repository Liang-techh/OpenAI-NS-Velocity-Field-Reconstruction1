"""Bind current Rp constants to actual old buffer and constructor AST leaves."""
import ast
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_native_constants as current
import lei_ren_part1_paper_compliant_outer_initial as initial
import lei_ren_part1_paper_compliant_outer_buffer as buffer
import lei_ren_part1_paper_compliant_axial_high_jets as high
import lei_ren_part1_paper_compliant_pulse_radial_C4 as radial
import lei_ren_part1_paper_logarithmic_outer_parameters as parameters_source

frame=current.frame;native=frame.native;leading_check=frame.bridge.leading_check
zero=leading_check.zero
PACKET=dict(E='Utheta_over_Pstar',m='Mz_over_R',h='Mtheta_over_sqrt2_R_3half_Pstar',
    k='Mtheta_z_over_sqrt2_R_3half_Pstar',e='Mztheta_over_R_Pstar_squared',p='Mp_over_Pstar_squared')


class Projector(native.Projection):
    def expression(self,node,env):
        if isinstance(node,ast.Call) and ast.unparse(node.func)=='get':
            assert len(node.args)==1 and not node.keywords
            return env['get'][node.args[0].value]
        if isinstance(node,ast.Call) and ast.unparse(node.func)=='box':
            assert len(node.args)==1 and not node.keywords
            # Exact mathematical lift of the identified source, not its box.
            return self.expression(node.args[0],env)
        return super().expression(node,env)


def source_frontend(owner,q,kernels):
    """Execute actual old scalar field assignments under true kernel bindings."""
    p=Projector();z=q.z;Q=1+z*z;invP2=s.exp(-2*q.logP)
    env=dict(z=z,qi=1/Q,yy=s.Integer(1),J=kernels['J1'],
        integrals=[kernels[key] for key in ('I_theta','I_pressure','I_energy')],
        **{'self.invP2':invP2})
    for key in ('factor','u','V','m','h','k','e','p'):
        env[key]=p.assignment(initial,'slope',key,env)
    old={'E':env['u'],**{key:env[key] for key in ('m','h','k','e','p')}}
    stages=[]
    def getters(values):return {PACKET[key]:value for key,value in values.items()}
    # y_d is exp(Md)+11. Turnoff is zero throughout the last eleven units.
    td=s.exp(40)+10
    env=dict(z=z,t=td,B=0,By=0,decay=s.exp(-td),root_decay=s.exp(-td/2),
        decay3=s.exp(-s.Rational(3,2)*td),get=getters(old),
        K=dict(B_mass=s.exp(-11)*kernels['K_B_at_cutoff'],
               B_squared_mass=s.exp(-11)*kernels['K_B2_at_cutoff']),
        **{'self.invP2':invP2})
    env['u1']=p.assignment(initial,'axial','u1',env)
    for key in ('u','m','h','k','e','p'):env[key]=p.assignment(initial,'axial',key,env)
    old={'E':env['u'],**{key:env[key] for key in ('m','h','k','e','p')}}
    stages.append('actual SharedOuterInitial.axial at y_d=exp40+11')
    env=dict(t=s.Integer(1),mu=q.mu,get=getters(old),
        kernels={key:kernels[name] for key,name in
            (('J','J_transition'),('theta','K_theta'),('energy','K_energy'),('pressure','K_pressure'))})
    for key in ('f','decay','d3','u1','u','m','h','k','e','p'):
        env[key]=p.assignment(buffer,'slope_mu',key,env)
    old={'E':env['u'],**{key:env[key] for key in ('m','h','k','e','p')}}
    stages.append('actual SharedOuterBuffer.slope_mu at t=1')
    Tw=q.at(owner.frame.functions['Tw']);env=dict(t=Tw,mu=q.mu,get=getters(old))
    for rate in (2*q.mu,1+2*q.mu):
        label=('decay_integral(c, 2 * mu, t)' if rate==2*q.mu
            else 'decay_integral(c, 1 + 2 * mu, t)')
        env[label]=(1-s.exp(-rate*Tw))/rate
    for key in ('slope','f','decay','d3','u1','u'):
        env[key]=p.assignment(buffer,'power',key,env)
    env['theta_kernel']=p.assignment(buffer,'power','theta_kernel',env,index=1,count=3)
    for key in ('m','h','k','e','p'):env[key]=p.assignment(buffer,'power',key,env)
    final={'E':env['u'],**{key:env[key] for key in ('m','h','k','e','p')}}
    stages.append('actual SharedOuterBuffer.power at positive Tw')
    coefficients=dict(U=s.cancel(final['E']*Q),M=s.cancel(final['m']/z),
        K=s.cancel(final['k']*Q/z),H=s.cancel(final['h']*Q),Pin=s.cancel(final['p']*Q*Q))
    poly=s.Poly(s.cancel(final['e']*Q*Q),z)
    coefficients.update(E_Z=poly.coeff_monomial(z**6),E_Q=poly.coeff_monomial(1))
    assert all(z not in value.free_symbols for value in coefficients.values())
    return p,final,coefficients,stages


def constructor(owner,q,coefficients,p):
    """Actual high-constant and radial H/P/Xp defining assignments."""
    c=coefficients;Tw=q.at(owner.frame.functions['Tw']);kernels=current.kernel_functions(owner,q)
    packet={PACKET['E']:[c['U']],PACKET['m']:[0,c['M']],PACKET['k']:[0,c['K']],
        PACKET['e']:[c['E_Q']],PACKET['h']:[c['H']],PACKET['p']:[c['Pin']]}
    env=dict(p0=packet,kernels={'B_squared_mass':s.exp(-11)*kernels['K_B2_at_cutoff']},
        **{'initial.params.yd':s.exp(40)+11,'initial.invP2':s.exp(-2*q.logP),
           'initial.params.logPstar':q.logP,'initial.params.Tw':Tw})
    for key in ('U','M','K','EQ','td','EZ','invP'):
        env[key]=p.assignment(high,'_incoming_constants',key,env)
    tree=ast.parse(Path(high.__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='_incoming_constants')
    multipliers=[n for n in ast.walk(fn) if isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='EZ']
    assert len(multipliers)==1 and isinstance(multipliers[0].op,ast.Mult)
    env['EZ']*=p.expression(multipliers[0].value,env)
    dictionaries=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
        and any(ast.unparse(t)=='self.constants' for t in n.targets)]
    assert len(dictionaries)==1 and isinstance(dictionaries[0],ast.Call)
    assert ast.unparse(dictionaries[0].func)=='dict' and not dictionaries[0].args
    values={kw.arg:kw.value for kw in dictionaries[0].keywords}
    needed=('U','M','K','E_Q','E_Z','C1','C2','C0','C_E')
    result={key:p.expression(values[key],env) for key in needed}
    assert all(key in values for key in ('td','Tw','B_squared_mass','retained_far_tail'))
    for target,key in (('self.inlet_H','H'),('self.inlet_P','Pin'),('self.Xp','Xp')):
        value=p.assignment(radial,'__init__',target,env)
        env[target]=value;result[key]=value
    calls=native.ast_assignments('axial_high_jets','CompliantAxialHighJets','_incoming_constants',dict(
        initial='self.base.pulse.initial',buffer='self.base.pulse.buffer',
        p0='buffer.power(0,1,cells=128)',
        kernels='turnoff_kernels(initial.ctx,initial.params.yd,initial.params.Md,cells=128)'))
    radial_call=native.ast_assignments('pulse_radial_C4','CompliantPulseRadialC4','__init__',
        {'p0':'self.pulse.buffer.power(0,1)'})
    definition=owner.frame.reports[native.NAME]['current_P0_source_binding']['analytic_definition']
    assert definition['Md']=='40','Same admitted Md=40 source family required'
    exact={'Md':s.Integer(definition['Md'])};actual_parameters={}
    exact['self.md']=p.assignment(parameters_source,'__init__','self.md',exact)
    assert exact['self.md']==40
    actual_parameters['self.md']=exact['self.md']
    for key in ('self.logPstar','self.log_mu','self.mu','self.Tw','self.yd'):
        value=p.assignment(parameters_source,'__init__',key,exact)
        exact[key]=value;actual_parameters[key]=value
    geometry=leading_check.Interpreter(owner.frame.bridge.leading)
    graph_parameters=owner.frame.bridge.leading.parameters
    geometry.bindings.pop(graph_parameters['logP'].node)
    geometry.bindings.pop(graph_parameters['mu'].node)
    assert zero(geometry.at(graph_parameters['logP'])-exact['self.logPstar'])
    assert zero(geometry.at(graph_parameters['mu'])-exact['self.mu'])
    assert zero(geometry.at(owner.frame.functions['Tw'])-exact['self.Tw'])
    assert exact['self.yd']==s.exp(40)+11
    parameters=dict(actual_parameter_functions={key:s.sstr(value) for key,value in actual_parameters.items()},
        same_admitted_pressure_source_Md=definition['Md'],
        same_actual_graph_logP_mu_Tw_definitions=True,
        input_hashes={Path(parameters_source.__file__).name:current.sha(Path(parameters_source.__file__).name)})
    # Include the full function binding, so the selected dictionary fields and
    # positive EZ multiplier cannot change independently of this proof.
    return result,dict(actual_high_constant_fields={key:ast.unparse(values[key]) for key in needed},
        actual_positive_EZ_multiplier=ast.unparse(multipliers[0]),
        actual_same_buffer_constructor_calls=calls,actual_radial_H_P_callback=radial_call,
        actual_y_d_and_Tw_parameter_source=parameters,
        source_function_boxes_denote_functions_not_chosen_interval_values=True,
        retained_tail_and_extra_source_fields_are_preserved_not_recomputed=True)


def run(owner=None):
    began=time.monotonic();raw=frame.bridge.outer.rh.read(current.NAME)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    owner=owner if owner is not None else current.CurrentOriginalRpNativeConstants(require_checked=False)
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.OPEN)
    q,constants=owner.current_constants();old,kernels=owner.original_native_functions(q)
    assert raw['actual_current_native_constants']=={key:s.sstr(value) for key,value in constants.items()}
    assert raw['actual_original_native_defining_functions']=={key:s.sstr(value) for key,value in old.items()}
    assert raw['actual_original_kernel_endpoint_functions']=={key:s.sstr(value) for key,value in kernels.items()}
    assert raw['historical_recipe_kernel_aliases']==current.LEGACY_KERNEL_ALIASES
    assert raw['turnoff_cutoff']=='exp(40)' and raw['true_yd']=='exp(40)+11'
    literal=leading_check.kernels(owner.frame.bridge.leading)
    predicates=leading_check.source_predicates();shared=native.kernel_bindings()
    assert literal['all_true_integral_nodes_checked']==16
    assert shared['same_original_continuous_kernel_providers']
    for key in constants:assert zero(constants[key]-old[key]),('original source coefficient',key)
    p,source,coefficients,stages=source_frontend(owner,q,kernels)
    for key in ('U','M','K','E_Q','E_Z','H','Pin'):
        assert zero(coefficients[key]-old[key]),('actual native buffer',key)
    installed,bindings=constructor(owner,q,coefficients,p)
    assert set(installed)==current.KEYS
    for key in installed:assert zero(installed[key]-constants[key]),('actual native constructor',key)
    for name,digest in {**p.hashes,**shared['input_hashes'],
            **bindings['actual_y_d_and_Tw_parameter_source']['input_hashes']}.items():
        assert owner.hashes[name]==digest,'Same checked defining source required: '+name
    pending=('selected_native_pulse_constructor_consumes_current_C3_frame',
        'native_interval_inlet_callback_installed','current_numeric_point_field_oracle_installed')
    assert not any(raw[key] for key in pending)
    result=dict(all_passed=True,source_family=owner.identity,**dict.fromkeys(current.GATES,True),
        true_continuous_kernel_integrand_and_endpoint_proof=literal,
        independent_source_kernel_AST_predicates=predicates,original_shared_kernel_provider_proof=shared,
        current_to_original_native_constant_function_identities=len(constants),
        actual_old_buffer_to_original_native_coefficient_identities=7,
        actual_native_constructor_to_current_constant_function_identities=len(installed),
        actual_old_scalar_field_stages=stages,actual_source_assignment_projection=p.bindings,
        actual_native_constructor_source_binding=bindings,
        yd_is_exp40_plus11_and_exact_K_B2_buffer_suffix_is_exp_minus11=True,
        positive_E_Z_from_its_true_source_not_energy_subtraction_or_padding=True,
        actual_source_functions_identified_not_interval_enclosure_equality=True,
        **dict.fromkeys(pending+current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Current/native cached constant functions and actual buffer/constructor identities passed',flush=True)
    return result


if __name__=='__main__':run()

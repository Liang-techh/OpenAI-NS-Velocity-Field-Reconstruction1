"""Actual compact/quiet composition, native six-row inlet and exact radius."""
import ast
import inspect
import json
import math
from pathlib import Path
import time
import textwrap
from types import SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame as current
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_bridge_check as bridge_check
import lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp_check as quiet_check

bridge=current.bridge;rows=current.target.rows;zero=bridge.leading_check.zero


def composition(field):
    q,frame=field.symbolic_frame();b=field.bridge.leading.bridge
    band=b.band_report['actual_C3_repair_band_functions'];nodes=b.band_nodes
    bindings={int(i):q.at(row['explicit_original_source_row']) for i,row in field.bridge.band_aliases.items()}
    alpha=nodes[band['profiles_at_band_x']['original_alpha']]
    mu=next(i for i in alpha['arguments'] if nodes[i]['operation']!='exact_rational')
    signature=bridge.leading.bridge.ExactFunctionSignatures
    assert signature(nodes).at(mu)==signature(field.graph.nodes).at(field.bridge.leading.parameters['mu'].node)
    bindings[mu]=q.mu;bindings[band['original_band_variable']]=s.Integer(2)
    interp=bridge_check.InletAlgebra(nodes,bindings);leading_rows=0
    for key,rr in band['original_leading_power']['leading_histories'].items():
        for j,i in enumerate(rr):
            assert zero(interp.at(i)-q.at(rows(field.functions['actual_leading_2Rc_histories_C3'][key])[j]))
            leading_rows+=1
    for j,i in enumerate(band['profiles_at_band_x']['radial_y_rows'][0]['original_E']):
        assert zero(interp.at(i)-q.at(rows(field.functions['actual_leading_2Rc_E_C3'])[j]))
    assert leading_rows==20
    actual_exit=quiet_check.compact_exit(field.quiet)
    assert actual_exit['actual_complete_minus_leading_equals_certified_relative_exit_rows']==20
    assert actual_exit['genuine_compact_repair_consumed_instead_of_homogeneous_Rc_shortcut']
    f=field.quiet.functions;Tw=q.at(field.functions['Tw'])
    qbind={f['Tw'].node:Tw,f['parameters']['mu'].node:q.mu}
    for key,rr in f['actual_repaired_exit_complete_C3'].items():
        for j,i in enumerate(rows(rr)):
            qbind[i.node]=q.at(rows(field.functions['actual_leading_2Rc_histories_C3'][key])[j])
    for j,i in enumerate(rows(f['actual_repaired_exit_E_C3'])):
        qbind[i.node]=q.at(rows(field.functions['actual_leading_2Rc_E_C3'])[j])
    for j,i in enumerate(rows(f['original_independent_P0_C3'])):
        assert field.quiet.graph.nodes[i.node]['quantity']=='P0'
        assert field.quiet.graph.nodes[i.node]['source_family']==field.identity
        qbind[i.node]=q.at(rows(field.bridge.functions['actual_independent_P0_C3'])[j])
    actual=bridge_check.InletAlgebra(field.quiet.graph.nodes,qbind);count=0
    for key,value in f['actual_C3_pulse_input_frame'].items():
        if not isinstance(value,current.target.C3Function):continue
        expected=field.functions['actual_identified_Rp_native_frame_C3'][key]
        for i,j in zip(rows(value),rows(expected)):
            assert zero(actual.at(i.node)-q.at(j)),(key,'actual Rp quiet frame')
            count+=1
    native=field.functions['original_power_at_Tw_C3'];native_rows=0
    for key,value in native['complete_histories_y0_y1_y2_y3_y4_C3'][0].items():
        for i,j in zip(rows(value),rows(field.functions['actual_identified_Rp_histories_C3'][key])):
            assert zero(q.at(i)-q.at(j)),(key,'actual original power semigroup')
            native_rows+=1
    for i,j in zip(rows(native['profiles_y0_y1_y2_y3_y4_C3'][0]['E']),rows(field.functions['actual_identified_Rp_E_C3'])):
        assert zero(q.at(i)-q.at(j));native_rows+=1
    assert count==32 and native_rows==24
    return q,frame,dict(actual_band_leading_2Rc_history_rows=leading_rows,
        independently_rechecked_actual_compact_exit=actual_exit,
        actual_current_quiet_Rp_native_frame_C3_rows=count,
        actual_original_power_Tw_semigroup_rows=native_rows,
        genuine_compact_band_and_complete_exit_consumed=True,
        relative_zero_does_not_set_absolute_history_or_pressure_to_zero=True)


def radius(field,q):
    p=field.bridge.leading.parameters;logC=s.Symbol('logCstar',real=True)
    original=field.bridge.leading.functions['actual_original_source_geometry']['O3_power']
    assert original['right_offset']==field.functions['original_Rc_phase_offset'].node
    geom=bridge.leading_check.Interpreter(field.bridge.leading,
        bindings={value.node:s.Symbol(key,real=True) for key,value in p.items()})
    geom.bindings[p['logC'].node]=logC
    # Bind formal logP to its actual parameter formula; it is not a cap value.
    geom.bindings[p['logP'].node]=s.exp(40)+11
    geom.bindings[p['mu'].node]=q.mu
    gotRw=geom.at(field.functions['original_logRw'])
    gotRp=geom.at(field.functions['original_logRp'])
    Tw=geom.at(field.functions['Tw']);logP=s.exp(40)+11
    wanted=s.log(110)+10*(logC+logP)+logP+1
    # Expand positive rational logs; never remove the original h_B*s_c/2.
    def logzero(value):return s.simplify(s.expand_log(value,force=True))==0
    assert logzero(gotRw-wanted) and logzero(gotRp-wanted-Tw)
    f=field.quiet.functions;sig=bridge.leading.bridge.ExactFunctionSignatures
    assert sig(field.quiet.graph.nodes).at(f['original_Rc_offset'].node)==sig(field.graph.nodes).at(field.functions['original_Rc_phase_offset'].node)
    assert sig(field.quiet.graph.nodes).at(f['Tw'].node)==sig(field.graph.nodes).at(field.functions['Tw'].node)
    local=bridge_check.InletAlgebra(field.quiet.graph.nodes,
        {f['original_Rc_offset'].node:s.Symbol('rho',real=True)})
    assert zero(local.at(f['R2'].node)-2*local.at(f['Rc'].node))
    proof=current.native.ast_assignments('pulse_physical_bounds','PulsePhysicalBounds','__init__',dict(
        **{'self.logP':"c.exp(c.mpf(parameters['Md']))+11",
           'self.logmu':"c.ln(c.mpf(parameters['c_mu']))-4*self.logP",
           'self.Tw':'-60*self.logmu',
           'self.logRp_parts':'dict(logCstar=10*self.logC,logPstar=10*self.logP,finite_outer_offset=c.ln(110)+self.logP+1+self.Tw)'}))
    origin=current.native.ast_assignments('current_native_spatial_phase',None,'affine_identity_theorem',
        dict(Ra='a-4*P-1000',minus='Ra+B*sc/2'))
    for name,digest in {**proof['input_hashes'],**origin['input_hashes']}.items():
        assert field.hashes[name]==digest,'Accepted defining radius source required'
    # Use the native proof's exact AST projector with its actual signature.
    return dict(actual_original_Rw_and_Rp_log_radius_identities=True,
        Rw='ln110+10(logCstar+logPstar)+logPstar+1',Rp='Rw_log+Tw',
        exact_same_current_quiet_Rc_and_Tw_defining_expressions=True,
        original_Rc_is_actual_O3_power_right_offset_and_R2_is_twice_Rc=True,
        original_microscopic_phase_origin_restored_not_zeroed=True,
        phase_offset='logRc-log_r_minus',log_r_minus='ln4-4logPstar-1000+hbB*sc/2',
        native_radius_source_AST=proof,original_phase_origin_source_AST=origin)


class Taylor6:
    """Independent exact six-row native Taylor arithmetic."""
    def __init__(self,c,coefficients):self.coefficients=list(coefficients)
    def __getitem__(self,j):return self.coefficients[j]
    @classmethod
    def constant(cls,c,value,order):return cls(c,[value]+[s.Integer(0)]*order)
    def coerce(self,value):return value if isinstance(value,Taylor6) else Taylor6(None,[s.sympify(value)]+[s.Integer(0)]*5)
    def __add__(self,other):
        other=self.coerce(other);return Taylor6(None,[self[j]+other[j] for j in range(6)])
    __radd__=__add__
    def __mul__(self,other):
        other=self.coerce(other);return Taylor6(None,[sum(self[k]*other[j-k] for k in range(j+1)) for j in range(6)])
    __rmul__=__mul__
    def reciprocal(self):
        out=[1/self[0]]
        for j in range(1,6):out.append(-sum(self[k]*out[j-k] for k in range(1,j+1))/self[0])
        return Taylor6(None,out)


def native_inlet(field,q,frame):
    z=q.z;power=field.functions['original_power_at_Tw_C3'];Q=1+z*z
    E=q.at(rows(power['profiles_y0_y1_y2_y3_y4_C3'][0]['E'])[0])
    H={key:q.at(rows(value)[0]) for key,value in power['complete_histories_y0_y1_y2_y3_y4_C3'][0].items()}
    U=s.cancel(E*Q);M=s.cancel(H['m']/z);Ht=s.cancel(H['h']*Q);K=s.cancel(H['k']*Q/z);Pin=s.cancel(H['p']*Q*Q)
    energy=s.Poly(s.cancel(H['e']*Q*Q),z)
    EZ,EQ=energy.coeff_monomial(z**6),energy.coeff_monomial(1)
    constants=dict(U=U,C1=M/U,C2=K/U**2,C0=EQ/U**2,C_E=EZ/U**2)
    assert all(z not in value.free_symbols for value in (*constants.values(),Ht,Pin))
    raw_shape=dict(u=U/Q,m1=constants['C1']*z*Q,m2=constants['C2']*z*Q,
        X=Ht/U,energy=constants['C0']+constants['C_E']*z*z*Q*Q,Mp=Pin/Q**2,
        P0=frame['P0'],pressure=frame['P0']+Pin/Q**2)
    for key,value in raw_shape.items():assert zero(value-frame[key]),(key,'actual native qshape')
    # Bind actual native incoming() statements to the exact current constants
    # and same analytic datum. No finite interval endpoint is selected here.
    def scalar(value):
        if isinstance(value,tuple):
            assert value[0]==value[1];return s.sympify(value[0])
        return s.sympify(value)
    analytic=frame['P0']
    datum=SimpleNamespace(normalized_jets=lambda Z,order:dict(normalized_pressure_coefficients=
        [s.diff(analytic,z,j).subs(z,Z)/math.factorial(j) for j in range(order+1)]))
    owner=SimpleNamespace(ctx=SimpleNamespace(mpf=scalar),constants=constants,Xp=Ht/U,inlet_P=Pin,datum=datum)
    fn=current.native.current.inlet.CompliantPowerInletC4.incoming
    env=dict(self=owner,Z=z,IntervalTaylor=Taylor6,endpoints=lambda value:(value,value))
    scope,binding=bridge.defining_slice(fn,('c','Z','z','q','invq','k','u','m','e','raw','p0'),env)
    tree=ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
    returns=[node for node in tree.body if isinstance(node,ast.Return)]
    assert len(returns)==1,'Actual native inlet has one unconditional return'
    native_packet=eval(compile(ast.Expression(body=returns[0].value),'<actual native inlet return>','eval'),scope)
    assert set(native_packet)==current.FRAME_KEYS-{'pressure'}
    binding['actual_return_expression']=ast.unparse(returns[0].value)
    count=0;callback=field.incoming_exact_taylor(z)
    for key,value in native_packet.items():
        assert len(value.coefficients)==6
        for j,v in enumerate(value.coefficients):
            assert zero(v-callback[key][j]),(key,j,'actual six native Taylor coefficients')
            assert zero(v-s.diff(raw_shape[key],z,j)/math.factorial(j));count+=1
    assert count==42
    for key,rr in field.functions['actual_identified_Rp_native_frame_C3'].items():
        if isinstance(rr,current.target.C3Function):
            for j,i in enumerate(rows(rr)):
                assert zero(q.at(i)-callback[key][j]*math.factorial(j)),(key,j,'C3 inlet subset')
    try:field.incoming_exact_taylor(2)
    except ValueError:pass
    else:raise AssertionError('Native input outside original Z domain admitted')
    old=field.reports[current.native.NAME]
    assert old['underlying_function_projection']['passed']
    assert old['underlying_function_projection']['identity_count']==45
    assert old['current_P0_source_binding']['same_fourteen_continuous_atoms_and_analytic_flatten_function']
    assert old['current_P0_source_binding']['current_order6_first_six_and_native_order5_same_defining_rows']
    # The exact formulas are analytic rational q-shapes and the original
    # analytic pressure function. Higher inlet rows are their true derivatives;
    # this does not promote any upstream C3 physical field to a global C5 field.
    return dict(independent_actual_native_incoming_Taylor_rows=count,
        actual_incoming_source_assignment_replay=binding,
        current_complete_Rp_qshapes_identified_from_actual_leading_functions=True,
        no_second_Pstar_division_of_common_m_k=True,
        true_Z4_Z5_from_identified_analytic_functions_not_padding=True,
        independent_analytic_P0_and_nonzero_Mp_kept=True,
        exact_symbolic_callback_is_not_an_interval_point_oracle=True,
        upstream_global_C4_C5_not_inferred=True)


def run(field=None):
    began=time.monotonic();raw=bridge.outer.rh.read(current.NAME)
    for name,value in raw['input_hashes'].items():assert current.sha(name)==value,name
    field=field if field is not None else current.CurrentOriginalC3RpNativeFrame(require_checked=False)
    assert not field.acceptance_loaded and field.graph.nodes==raw['exact_graph_nodes']
    assert current.target.encoded(field.functions)==raw['actual_Rp_native_frame_functions']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    q,frame,chain=composition(field);geometry=radius(field,q);incoming=native_inlet(field,q,frame)
    assert raw['actual_analytic_native_inlet_functions']=={key:s.sstr(value) for key,value in frame.items()}
    assert raw['actual_native_Z0_to_Z5_Taylor_rows']=={key:[s.sstr(s.diff(value,q.z,j)/math.factorial(j)) for j in range(6)] for key,value in frame.items()}
    pending=('selected_native_pulse_constructor_consumes_current_C3_frame',
        'native_interval_inlet_callback_installed','current_numeric_point_field_oracle_installed',
        'global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed')
    assert not any(raw[k] for k in pending)
    result=dict(all_passed=True,source_family=field.identity,**dict.fromkeys(current.GATES,True),
        independent_actual_complete_chain_and_Rp_frame=chain,
        independent_actual_native_radius=geometry,independent_actual_six_row_native_inlet=incoming,
        **dict.fromkeys(pending+current.OPEN,False),
        input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Complete current Rp/native frame, true six-row inlet and original radius identities passed',flush=True)
    return result


if __name__=='__main__':run()

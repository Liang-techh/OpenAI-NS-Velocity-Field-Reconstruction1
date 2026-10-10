"""Independent analytic derivatives and actual native interval caller evidence."""
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_interval_inlet as current
import lei_ren_part1_paper_compliant_current_original_Rp_native_constants_check as constant_check
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rational_interval(expr,c,bindings):
    if expr.is_Rational:return c.mpf(int(expr.p))/int(expr.q)
    if expr.is_Symbol:return bindings[str(expr)]
    if expr.is_Add:return sum((rational_interval(x,c,bindings) for x in expr.args),c.mpf(0))
    if expr.is_Mul:
        result=c.mpf(1)
        for x in expr.args:result*=rational_interval(x,c,bindings)
        return result
    if expr.is_Pow and expr.exp.is_Integer:
        return rational_interval(expr.base,c,bindings)**int(expr.exp)
    raise ValueError('Only exact rational derivative expressions admitted')


def overlap(a,b,label):
    al,ah=current.endpoints(a);bl,bh=current.endpoints(b)
    if max(al,bl)>min(ah,bh):raise ArithmeticError('Independent derivative outside packet: '+label)


def sample_source_identity(owner):
    """Bind actual Z=1/2 H/Pin extraction to the current whole-Z functions."""
    q,constants=owner.exact.current_constants()
    kernels=current.exact.kernel_functions(owner.exact,q)
    projector,source,coefficients,stages=constant_check.source_frontend(owner.exact,q,kernels)
    z=q.z;Q=1+z*z;half=s.Rational(1,2)
    zero=constant_check.zero
    assert zero(source['h']*Q-constants['H'])
    assert zero(source['p']*Q*Q-constants['Pin'])
    env=dict(q=s.Rational(5,4),**{
        "read_interval(c, sample['Mtheta_over_sqrt2_R_3half_Pstar'][0])":source['h'].subs(z,half),
        "read_interval(c, sample['Mp_over_Pstar_squared'][0])":source['p'].subs(z,half),
        'self.constants':{'U':constants['U']}})
    for target,key in (('self.inlet_H','H'),('self.inlet_P','Pin'),('self.Xp','Xp')):
        value=projector.assignment(current.original,'__init__',target,env)
        assert zero(value-constants[key]),('actual sample extraction',key)
        env[target]=value
    sample_file='lei_ren_part1_paper_compliant_outer_buffer.json'
    buffer=json.loads((current.HERE/sample_file).read_bytes());sample=buffer['samples'][-1]
    assert sample['pulse_inlet'] and current.endpoints(current.original.read_interval(owner.ctx,sample['Z']))==(s.Rational(1,2),s.Rational(1,2))
    assert buffer['actual_five_defect_family_sha256']==owner.family
    assert buffer['implicit_source_sha256']==owner.source and buffer['datum_enclosure_sha256']==owner.datum_sha
    for name,digest in projector.hashes.items():assert owner.hashes[name]==digest,name
    assert owner.hashes[sample_file]==current.sha(sample_file)
    return dict(whole_Z_native_H_and_Pin_function_identities=2,
        actual_power_inlet_sample_extraction_function_identities=3,
        actual_sample_Z='1/2',same_source_family_and_buffer_report=True,
        exact_shapes='Mtheta=H/(1+Z^2); Mp=Pin/(1+Z^2)^2',
        actual_source_assignment_projection=projector.bindings,
        underlying_function_identity_not_inferred_from_sample_intervals=True)


@source_precision
def run(owner=None):
    began=time.monotonic();raw=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    owner=owner if owner is not None else current.CurrentOriginalRpIntervalInlet(require_checked=False)
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.PENDING)
    assert raw['source_family']==owner.family_record and all(owner.ownership().values())
    exact=json.loads((current.HERE/current.exact.RECEIPT).read_bytes())
    assert exact['all_passed'] and all(exact[key] for key in current.exact.GATES)
    assert exact['current_to_original_native_constant_function_identities']==12
    assert exact['actual_native_constructor_to_current_constant_function_identities']==12
    assert exact['actual_source_functions_identified_not_interval_enclosure_equality']
    sample_identity=sample_source_identity(owner)
    c=owner.ctx;z=s.Symbol('z');q=1+z*z
    U,C1,C2,C0,CE,Xp,Pin=s.symbols('U C1 C2 C0 C_E Xp Pin')
    functions=dict(u=U/q,m1=C1*z*q,m2=C2*z*q,energy=C0+CE*z*z*q*q,X=Xp,Mp=Pin/q**2)
    derivatives={key:[s.factor(s.diff(expr,z,n)/s.factorial(n)) for n in range(6)]
                 for key,expr in functions.items()}
    checks={};caller_records=[];row_count=0
    for name,view in raw['source_bound_six_row_views'].items():
        Zbox=current.original.read_interval(c,view['Z'])
        packet=owner.current_power(Zbox,0);incoming=packet['canonical_incoming_Taylor']
        assert set(incoming)==current.FRAME_KEYS
        bindings={**owner.interval_constants,'z':Zbox}
        count=0
        for key,rows in derivatives.items():
            assert incoming[key].order==5
            for n,expr in enumerate(rows):
                independent=rational_interval(expr,c,bindings)
                overlap(incoming[key][n],independent,name+'/'+key+'/'+str(n));count+=1
        p0=owner.datum.normalized_jets(Zbox,5)['normalized_pressure_coefficients']
        for n,value in enumerate(p0):
            overlap(incoming['P0'][n],value,name+'/P0/'+str(n));count+=1
            overlap(incoming['pressure'][n],incoming['P0'][n]+incoming['Mp'][n],name+'/pressure/'+str(n));count+=1
        # The actual inherited power method publishes this callback's objects.
        assert packet['pressure']['P0_over_Pstar_squared'] is incoming['P0']
        assert packet['pressure']['original_analytic_datum_preserved']
        assert owner.calls[-1]['caller']=='CompliantPowerInletC4.power'
        assert owner.calls[-1]['rows']==48
        for key,serialized in view['incoming'].items():
            assert len(serialized['coefficients'])==6
        assert count==48
        checks[name]=dict(independent_true_derivative_rows=count,
            actual_native_power_callback_used=True,independent_P0_preserved=True,
            interval_comparison_is_diagnostic_not_function_identity=True)
        caller_records.append(owner.calls[-1]);row_count+=count
    assert len(checks)==3 and row_count==144 and all(owner.ownership().values())
    for key in ('td','Tw','B_squared_mass','retained_far_tail'):
        assert owner.constants[key] is owner.before.constants[key]
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        exact_current_native_constant_function_identity_receipt=current.exact.RECEIPT,
        original_directed_source_enclosure_receipt=current.REUSED,
        actual_same_buffer_sample_H_Pin_identity=sample_identity,
        actual_native_power_callback_records=caller_records,current_interval_inlet_owner_graph=owner.ownership(),
        independent_analytic_derivative_checks=checks,independent_true_Taylor_rows=row_count,
        interval_objects_retained_without_midpoint_or_cap_selection=True,
        positive_E_Z_and_original_retained_tail_source_preserved=True,
        exact_identity_not_inferred_from_interval_overlap=True,
        scope='Actual local native power inlet caller; selected pulse/flatten injection still open',
        **dict.fromkeys(current.PENDING,False),input_hashes={**owner.hashes,
            current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.encode(current.pack(result)),indent=2)+'\n',
        encoding='utf8',newline='\n')
    print('Current interval inlet: 144 independent Taylor rows and actual native caller passed',flush=True)
    return result


if __name__=='__main__':run()

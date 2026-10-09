"""Same-source R100 switch inlet, with physical factors kept formal.

Reuses accepted complete macro endpoint integrals. Q is recovered through
ordinary order4; physical pressure axis and its radial increment are separate.
This adapter does not integrate either microscopic switch or close R110.
"""
import ast
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_bridge_macro_moments as moments
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields=moments.fields
HERE,PREFIX,sha=moments.HERE,moments.PREFIX,moments.sha
prior,base,ep=moments.prior,moments.base,moments.ep
NAME=PREFIX+'current_original_R100_endpoint.json'
RECEIPT=PREFIX+'current_original_R100_endpoint_check.json'
GATE='original_factored_R100_actual_endpoint_switch_inlet_installed'


def source_bindings():
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    bindings=assignment_source_bindings('actual_bridge_integrals','packet',{
        'Q':"(2*z*V-z*own['actual']['M']*(1-self.core.delta)-d*derivative(own['actual']['M']))/L",
        'pressure':"dress(own['actual']['C'],p['inputs']['F0_squared_ratios'])"})
    bindings.update(assignment_source_bindings('inner_switch_profiles','packet',{
        'Q':"(2*z*v-(z*m)*(1-self.core.delta)-d*derivative(m))/L"}))
    tree=ast.parse((HERE/(PREFIX+'actual_bridge_integrals.py')).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='packet')
    expected=ast.dump(ast.parse("{n:list(j.coefficients) for n,j in own['actual'].items()}",mode='eval').body)
    exports=[kw.value for node in ast.walk(fn) if isinstance(node,ast.Call)
             for kw in node.keywords if kw.arg=='actual_own_six_moments_axial5']
    if len(exports)!=1 or ast.dump(exports[0])!=expected:
        raise ValueError('Original materialized ordinary normalized moment export changed')
    return dict(passed=True,original_Q_and_pressure_assignments=bindings,
        ordinary_normalized_inlet_rows_not_hidden_amplitude_packets=True,
        original_radius='q=1; 2hb+(log(100/Ra)-2hb)=log(100/Ra); R=100',
        original_velocity_prefactors=dict(Ur='sqrt(R/2)*Q',Utheta='sqrt(2R)*F0*phi',Uz='V'),
        physical_moment_units='Mtheta=F0*R^2*H; Mz=R*M; Mtheta_z=F0*R^2*K; Mztheta=R*A-R^2*F0^2*B; Mp=R*F0^2*C')


def restore_row(flow,record):
    if (not isinstance(record,dict) or record.get('point_value_selected') is not False
            or record.get('encloses_original_source_function') is not True):
        raise ValueError('Accepted factored original source enclosure required')
    scale=record['formal_positive_scale']
    powers=tuple(scale['source_exponents'])+(scale['radius_power'],)
    value=prior.ScaledEnclosure(prior.FormalScale(flow.logs,powers,
        fields.previous.read_interval(flow.c,scale['additional_log_interval'])),
        fields.previous.read_interval(flow.c,record['coefficient_interval']),flow.ledger)
    if value.zero!=record['exact_zero']:raise ValueError('Saved exact-zero predicate disagrees')
    return value


def recover_endpoint(flow,Z,delta,p0,amplitude_ratios,amplitude_squared_ratios,actual_fields,actual_histories):
    """Original endpoint algebra; input ratios are ordinary coefficients."""
    c=flow.c;Z=c.mpf(Z);delta=c.mpf(delta)
    if ep(Z)[0]<-1 or ep(Z)[1]>1 or ep(delta)[0]<0 or ep(delta)[1]>=1:
        raise ValueError('Original |Z|<=1 and 0<=delta<1 required')
    if set(actual_fields)!=set(('phi','V')) or set(actual_histories)!=set(moments.RATES):
        raise ValueError('Both actual fields and all six own histories required')
    for row in (*actual_fields.values(),*actual_histories.values()):
        if len(row)!=6 or any(not isinstance(v,prior.ScaledEnclosure) or
                v.scale.bases is not flow.logs or v.ledger is not flow.ledger for v in row):
            raise ValueError('One source basis/ledger and ordinary Z0..5 required')
    for jet in (p0,amplitude_ratios,amplitude_squared_ratios):
        if not isinstance(jet,IntervalTaylor) or jet.ctx is not c or jet.order!=5:
            raise ValueError('Same-context original pressure and amplitude jets0..5 required')
        if any(not mp.isfinite(x) for v in jet.coefficients for x in ep(v)):
            raise ValueError('Finite source coefficient enclosures required')
    if ep(amplitude_ratios[0])!=(1,1) or ep(amplitude_squared_ratios[0])!=(1,1):
        raise ValueError('Amplitude ratios are divided by the base amplitude exactly once')
    z=IntervalTaylor.variable(c,Z,5);d=1-z*z;L=1-z*z*delta
    if ep(L[0])[0]<=0:raise ValueError('Positive original radial recovery denominator required')
    zrows=flow.jet(z);drows=flow.jet(d);invL=flow.jet(L.reciprocal())
    M=actual_histories['M'];V=actual_fields['V']
    # Only M0..5 are known. The padded unused derivative5 must be discarded.
    M_Z=[M[n+1]*(n+1) for n in range(5)]+[flow.scalar(0)]
    numerator=flow.add(flow.scale(flow.multiply(zrows,V),2),
        flow.scale(flow.multiply(zrows,M),-(1-delta)),
        flow.scale(flow.multiply(drows,M_Z),-1))
    Q=flow.multiply(numerator,invL)[:5]
    ratio=flow.jet(amplitude_ratios);ratio2=flow.jet(amplitude_squared_ratios)
    dressed_phi=flow.multiply(actual_fields['phi'],ratio)
    dressed={name:flow.multiply(actual_histories[name],ratio if name in ('H','K') else ratio2)
             for name in ('H','K','B','C')}
    R=flow.scalar(100);F0=flow.factor((0,0,.5,0,0));F02=flow.factor((0,0,1,0,0));P2=flow.factor((0,1,0,0,0))
    pressure_axis=flow.scale(flow.jet(p0),P2)
    pressure_increment=flow.scale(dressed['C'],R*F02)
    physical_moments=dict(Mtheta=flow.scale(dressed['H'],R*R*F0),
        Mz=flow.scale(M,R),Mtheta_z=flow.scale(dressed['K'],R*R*F0),
        Mztheta=flow.add(flow.scale(actual_histories['A'],R),flow.scale(dressed['B'],-R*R*F02)),
        Mp=pressure_increment)
    return dict(radius=R,normalized_actual_fields=actual_fields,
        normalized_actual_own_six_moments=actual_histories,
        F_actual_true_axial5_divided_by_F0=dressed_phi,
        actual_Q_axial4_coefficients=Q,
        pressure_axis_over_Pstar_squared_axial5=flow.jet(p0),
        pressure_increment_true_axial5_divided_by_R_F0_squared=dressed['C'],
        original_F0_derivative_ratios_ordinary=ratio,original_F0_squared_derivative_ratios_ordinary=ratio2,
        physical_velocity_axial_coefficients=dict(Ur=[v*c.sqrt(50) for v in Q],
            Utheta=flow.scale(dressed_phi,F0*c.sqrt(200)),Uz=V),
        physical_pressure_axis_axial5=pressure_axis,
        physical_pressure_radial_increment_axial5=pressure_increment,
        physical_total_pressure_axial5=flow.add(pressure_axis,pressure_increment),
        physical_cumulative_moment_axial5=physical_moments,
        analytic_P0_kept_separate=True,actual_moments_not_comparison=True,
        amplitude_derivatives_dressed_once=True,physical_amplitudes_not_materialized=True,
        Q_ordinary_orders=list(range(5)),other_ordinary_orders=list(range(6)),
        same_original_endpoint_first_switch_schema_ready=True)


class OriginalR100Endpoint:
    mode='genuine_original_R100_factored_actual_fields_moments_Q_pressure'
    def __init__(self,dps=500):
        self.upstream=moments.OriginalBridgeMacroMoments(dps);self.c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes);self.saved=None
        for name in (moments.NAME,moments.RECEIPT):
            row=json.loads((HERE/name).read_bytes())
            if not row.get(moments.GATE) or name==moments.RECEIPT and not row.get('all_passed'):
                raise ValueError('Accepted actual six-history functions required')
            if row['source_family']!=self.family:raise ValueError('Same original moment source family required')
            for path,digest in row['input_hashes'].items():fields.previous.bind(self.hashes,path,digest)
            fields.previous.bind(self.hashes,name,sha(name))
            if name==moments.NAME:self.saved=row
        for stem in ('actual_bridge_integrals','inner_switch_profiles','actual_switch_mixed_C4'):
            name=PREFIX+stem+'.py';fields.previous.bind(self.hashes,name,sha(name))
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.bindings=source_bindings();self.owners={}

    def owner(self,label):
        if label in self.owners:return self.owners[label]
        moment_owner,proof=self.upstream.owner(label);flow=moment_owner.flow;c=self.c
        saved=self.saved['packets'][label][-1]
        if (saved['source_frame']!=label or saved['source_family']!=self.family
                or saved['implicit_source_sha256']!=self.upstream.fields.source
                or saved['datum_enclosure_sha256']!=self.upstream.fields.datum
                or saved['evaluation']['geometry']['fraction']!=[1,1]):
            raise ValueError('Same-source actual q1 endpoint required')
        e=saved['evaluation'];radius=restore_row(flow,e['geometry']['R1'])
        if radius.scale.powers!=(0,0,0,0,0) or ep(radius.scale.offset)!=(0,0) or ep(radius.coefficient)!=(100,100):
            raise ValueError('Exact original R100 endpoint radius required')
        restore=lambda rows:[restore_row(flow,v) for v in rows]
        actual_fields={key:restore(row) for key,row in e['same_original_macro_field_functions'].items()}
        actual_histories={key:restore(row) for key,row in e['actual_six_moment_functions'].items()}
        gradient=proof['original_gradient_coefficients']
        Lambda=fields.previous.read_interval(c,self.upstream.fields.records['core_transfer']['Lambda'])
        def ratios(multiplier):
            out=[c.mpf(1)]
            for n in range(1,6):out.append(sum((-multiplier*Lambda*gradient[j]*out[n-1-j]
                                               for j in range(n)),c.mpf(0))/n)
            return IntervalTaylor(c,out)
        p0=IntervalTaylor(c,proof['original_P0_coefficients'][:6])
        delta=c.exp(fields.previous.read_interval(c,self.upstream.fields.records['pressure_source']
            ['compliant_source']['parameter_bounds']['log_delta']))
        value=recover_endpoint(flow,proof['Z'],delta,p0,ratios(1),ratios(2),actual_fields,actual_histories)
        self.owners[label]=(flow,proof,value);return self.owners[label]

    def evaluate(self,label):
        with mp.workdps(self.c.dps+40):flow,proof,value=self.owner(label)
        return dict(mode=self.mode,source_frame=label,Z=proof['Z'],source_family=self.family,
            implicit_source_sha256=self.upstream.fields.source,datum_enclosure_sha256=self.upstream.fields.datum,
            function_evaluation=fields.serialized(value),source_log_basis=list(flow.logs),
            source_ledger=dict(flow.ledger),accepted_complete_macro_endpoint_reused=True,
            cached_values_remain_source_enclosures_not_points=True,
            actual_R100_R110_switch_function_installed=False,whole_axis_functions_installed=False,
            actual_micro_function_provider_installed=False)


def run():
    began=time.monotonic();owner=OriginalR100Endpoint()
    result=dict(**{GATE:True},source_family=owner.family,source_bindings=owner.bindings,
        packets={label:owner.evaluate(label) for label in ('0','.5')},
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Original factored R100 actual F/V, six histories, Q0..4, physical velocity/moments and separate pressure0..5. Switch functions, micro/whole-Z and R110 remain open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original factored R100 actual switch inlet recovered',flush=True);return result


if __name__=='__main__':run()

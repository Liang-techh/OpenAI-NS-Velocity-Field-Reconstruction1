"""True native/Gamma C5 and leading pulse jets feeding n1 pressure C3.

Existing C1/C2/C3 handles are preserved. The new native rows follow exact
whole-function identities, and Gamma rows use full improper kernels.
Numerical ranges and interval endpoints never define a derivative.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_pulse_native_numeric_identity as native

true=native.true
jets=true.own.jets
source,controls,parent=true.source,true.controls,true.parent
c3=jets.c3
HERE,PREFIX,sha=native.HERE,native.PREFIX,native.sha
NAME=PREFIX+'current_original_n1_pulse_pressure_C3_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_pulse_pressure_C3_functions_check.json'
GATES=('current_original_n1_actual_native_and_Gamma_C5_functions_installed',
    'current_original_n1_actual_leading_pulse_U_mean_C5_Q_C4_functions_installed',
    'current_original_n1_true_pulse_pressure_C3_functions_installed')
OPEN=tuple(dict.fromkeys(native.OPEN+(
    'current_original_n1_complete_signed_defect_C3_functions_installed',
    'n1_repair_velocity_C3_functions_installed','required_C3_C5_derivative_mixed_bounds')))


@dataclass(frozen=True)
class C5Function:
    value:source.FunctionRef
    Z:source.FunctionRef
    ZZ:source.FunctionRef
    ZZZ:source.FunctionRef
    Z4:source.FunctionRef
    Z5:source.FunctionRef


@dataclass(frozen=True)
class C4Function:
    value:source.FunctionRef
    Z:source.FunctionRef
    ZZ:source.FunctionRef
    ZZZ:source.FunctionRef
    Z4:source.FunctionRef


def rows(value):
    if type(value) in (C5Function,C4Function,c3.C3Function):return tuple(vars(value).values())
    if type(value) is source.C1Function:return (value.value,value.Z)
    if type(value) is true.own.c2.C2Function:return (value.value,value.Z,value.ZZ)
    raise TypeError('A genuine C1/C2/C3/C4/C5 source function is required')


def encode(value):
    if type(value) in (C5Function,C4Function):return [q.node for q in rows(value)]
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [encode(v) for v in value]
    return jets.encode(value)


def native_C5_functions(g,old):
    """Extend only scalars identified by exact native whole-function bodies."""
    z=g.symbol('Z');c=g.constant;Q=g.add(g.one,g.mul(z,z))
    U=g.mul(old['u'].value,Q)
    den=g.mul(c(2),g.add(g.one,g.mul(c(12),z,z),g.mul(c(15),z,z,z,z)))
    CE=g.quotient(old['energy'].ZZ,den,'2*(1+12Z^2+15Z^4) >= 2 on real Z')
    symbol=sy.Symbol('Z',real=True)
    u_extra=[g.mul(U,source.expression(g,sy.diff(1/(1+symbol**2),symbol,k),{'Z':z})) for k in (4,5)]
    E4=g.mul(CE,g.add(c(48),g.mul(c(360),z,z)));E5=g.mul(c(720),CE,z)
    result={key:C5Function(*rows(old[key]),*extra) for key,extra in dict(
        u=u_extra,m1=(g.zero,g.zero),m2=(g.zero,g.zero),energy=(E4,E5),X=(g.zero,g.zero)).items()}
    return result,dict(actual_U_from_native_value=U,actual_CE_from_native_ordinary_Z2=CE,
        positive_CE_denominator=den,cubic_moment_Z4_Z5_are_proved_zero_not_padding=True,
        X_higher_rows_require_actual_constant_native_X_identity=True,
        existing_native_C3_handles_preserved=True)


def Gamma_C5_functions(g,owner,old):
    """Fourth/fifth Dhat rows from true M4/M5 expectation kernels."""
    a=owner.pieces['a'];S=owner.schedule['S']
    t,v=owner.pieces['offset_variable'],owner.pieces['Gamma_variable']
    xi=owner.pieces['xi'];z=g.symbol('Z');c=g.constant
    normal=g.unary('gamma',g.add(g.one,a));log=g.unary('log1p',g.mul(xi,v))
    moments=[]
    for k in range(1,6):
        ak=g.add(a,c(k))
        body=g.unary('exp',g.add(g.neg(v),g.mul(ak,g.unary('log',v)),g.neg(g.mul(ak,log))))
        integral=controls.integral(g,body,g.nodes[v.node]['name'],g.zero,owner.pieces['infinity'],
            exact_positive_Gamma_expectation=True,full_improper_integral=True,
            Gamma_normalization=normal.node,Gamma_moment_order=k,zero_endpoint_integrand_continuation=0)
        moments.append(g.quotient(integral,normal,'positive Gamma(1+a)'))
    _,M2,M3,M4,M5=moments
    A,B,C,D=(g.add(a,c(k)) for k in (1,2,3,4))
    T=lambda k:g.unary('exp',g.mul(c(-k),t))
    fourth=g.add(g.mul(c(-48),A,S,T(2),M2),
        g.mul(c(-384),A,B,S,S,T(3),z,z,M3),
        g.mul(c(-256),A,B,C,S,S,S,T(4),z,z,z,z,M4))
    fifth=g.add(g.mul(c(-960),A,B,S,S,T(3),z,M3),
        g.mul(c(-2560),A,B,C,S,S,S,T(4),z,z,z,M4),
        g.mul(c(-1024),A,B,C,D,S,S,S,S,T(5),z,z,z,z,z,M5))
    return C5Function(*rows(old),fourth,fifth),dict(positive_improper_moments=moments,
        actual_xi=xi,normalization=normal,offset=t,actual_nonzero_S=S,
        exact_Z4_Z5_chain_rules_from_same_Dhat=True,
        endpoint_Z_plus_minus_one_not_divided_by_one_minus_Z_squared=True,
        Gamma_moment_domination='M_k <= Gamma(1+a+k)/Gamma(1+a), k=1..5, xi>=0')


class HigherZDerivative(jets.OrdinaryZDerivative):
    """The same actual source derivative, with unambiguous integral order labels."""
    def rule(self,i):
        row=self.g.nodes[i]
        result=super().rule(i)
        if row['operation']=='definite_integral' and result!=self.g.zero:
            attrs=copy.deepcopy(self.g.nodes[result.node]);op=attrs.pop('operation')
            base=row.get('original_slow_Z_derivative_base_order',row.get('ordinary_slow_Z_derivative_order',0))
            count=row.get('ordinary_Z_derivative_order',0)+1
            attrs.update(original_slow_Z_derivative_base_order=base,
                ordinary_slow_Z_derivative_order=base+count,ordinary_Z_derivative_order=count,
                actual_higher_source_body_not_interval_endpoint=True)
            return self.g.node(op,**attrs)
        return result
    def seed_rows(self,value):
        for left,right in zip(rows(value),rows(value)[1:]):self.seed(left,right)
    def extend(self,value,order):
        existing=list(rows(value))
        while len(existing)<=order:existing.append(self.derivative(existing[-1]))
        if order==5:return C5Function(*existing)
        if order==4:return C4Function(*existing)
        if order==3:return c3.C3Function(*existing)
        raise ValueError('Only required genuine C3/C4/C5 lifts supported')


class CurrentOriginalN1PulsePressureC3Functions:
    def __init__(self,native_owner,require_checked=True):
        if type(native_owner) is not native.CurrentOriginalN1PulseNativeNumericIdentity or not native_owner.acceptance_loaded:
            raise TypeError('Checked actual native-to-pulse function identity required')
        native_owner.assert_graph();self.native_owner=native_owner;self.pulse=native_owner.pulse
        self.g=g=self.pulse.g;self.family=copy.deepcopy(self.pulse.family);self.hashes=dict(native_owner.hashes)
        for name in (native.NAME,native.RECEIPT):self._pin(name,sha(name))
        for name,digest in parent.linear.read(native.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        self.old=self.pulse.own.jets
        self.native,self.native_evidence=native_C5_functions(g,self.old.native)
        self.Gamma_Dhat,self.Gamma_evidence=Gamma_C5_functions(g,self.old.complete.Gamma,self.old.Gamma_Dhat)
        self.derivative=d=HigherZDerivative(g)
        d.derivatives=dict(self.old.differentiator.derivatives)
        for jet in self.native.values():d.seed_rows(jet)
        d.seed_rows(self.Gamma_Dhat)
        # These are exact Z-independent defining scalars, after the source
        # whole-function identities checked in the separate proof receipt.
        for key in ('actual_U_from_native_value','actual_CE_from_native_ordinary_Z2'):
            d.seed(self.native_evidence[key],g.zero)
        self.selected_pulse={key:d.extend(value,5) for key,value in self.old.functions['pulse'].items()}
        for jet in self.selected_pulse.values():d.seed_rows(jet)
        self.leading_fields={}
        for chart,fields in self.pulse.own.fields.items():
            self.leading_fields[chart]={key:d.extend(fields[key],5) for key in ('Utheta','U','U_y','m','m_y')}
            self.leading_fields[chart].update({key:d.extend(fields[key],4) for key in ('Q','Q_y')})
        self.integrals={chart:{key:d.extend(value,3) for key,value in functions.items()}
            for chart,functions in self.pulse.bindings['integrals'].items()}
        self.endpoints={chart:[{key:d.extend(value,3) for key,value in end.items()} for end in ends]
            for chart,ends in self.pulse.bindings['per_chart_finite_endpoints'].items()}
        self.global_endpoints=[{key:d.extend(value,3) for key,value in end.items()}
            for end in self.pulse.bindings['global_finite_endpoints']]
        self.pressure={key:d.extend(value,3) for key,value in self.pulse.functions.items()}
        self.pieces={chart:{key:d.extend(value,3) for key,value in functions.items()}
            for chart,functions in self.pulse.pieces.items()}
        self.definition=dict(native_whole_function_C5_rows=True,actual_Gamma_M4_M5_improper_kernels=True,
            leading_main_gap_end_U_mean_C5_and_auxiliary_Q_C4=True,
            actual_signed_n1_pulse_pressure_C3_integral_functions=True,
            existing_native_and_pulse_lower_order_handles_preserved=True,
            physical_Pstar_radial_measures_and_finite_endpoints_retained=True,
            leading_auxiliary_Q_mean_weight='1-delta',
            original_axis_pressure_datum_and_other_full_defect_segments_not_changed=True,
            pressure_derivative_numerical_containment_transferred=False,
            complete_five_defect_C3_and_reserved_coefficient_C3_not_installed=True,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True,
            global_n1_and_required_mixed_norms_complete=False)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.native_owner,self.pulse,self.g);self._family=copy.deepcopy(self.family)
        self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES)
                or any(receipt.get(k) for k in OPEN) or receipt['source_family']!=self.family
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or archive!=sha(NAME)
                or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['C5_pulse_C3_pressure_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked genuine pulse pressure C3 functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed higher pulse source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed higher derivative source '+name)
        self.hashes[name]=digest

    def _snapshot(self):
        return encode(dict(native=self.native,native_evidence=self.native_evidence,Gamma_Dhat=self.Gamma_Dhat,
            Gamma_evidence=self.Gamma_evidence,selected_pulse=self.selected_pulse,
            leading_fields=self.leading_fields,integrals=self.integrals,endpoints=self.endpoints,
            global_endpoints=self.global_endpoints,pressure=self.pressure,pieces=self.pieces,definition=self.definition))

    def assert_graph(self):
        self.native_owner.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or (self.native_owner,self.pulse,self.g)!=self._owners or self.family!=self._family):
            raise ValueError('Actual higher pulse functions or source bindings changed')


def run(native_owner):
    began=time.monotonic();owner=CurrentOriginalN1PulsePressureC3Functions(native_owner,False)
    data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,C5_pulse_C3_pressure_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_NATIVE_GAMMA_C5_AND_PULSE_PRESSURE_C3_DEFINED',len(owner.g.nodes),flush=True)
    return owner,data


if __name__=='__main__':raise SystemExit('run(native_owner) reuses the checked live function source identity')

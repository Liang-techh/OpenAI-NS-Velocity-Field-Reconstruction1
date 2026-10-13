"""Directed cut n1 own fields and five inner moments of the true graph.

The true cutoff graph defines the fields. Accepted numerical regular,
radial and cumulative owners enclose those definitions, after explicit
phase, flat-function, axis, signed-density and history identity guards.
No reserved repair or exterior pressure completion is applied here.
"""
import ast
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge as inner
import lei_ren_part1_paper_compliant_current_original_n1_radial_cutoff_extension as legacy
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

full=inner.full;true=full.cutoff;op=legacy.op;ep=legacy.ep
HERE,PREFIX,sha=inner.HERE,inner.PREFIX,inner.sha
NAME=PREFIX+'current_original_n1_cutoff_numeric_bridge.json.gz'
RECEIPT=PREFIX+'current_original_n1_cutoff_numeric_bridge_check.json'
GATE='current_original_true_cutoff_n1_own_fields_and_five_inner_moment_enclosures_installed'
OPEN=inner.OPEN


def normalize_binary(coefficient,exponent=0):
    """Keep m*2**e exact without allocating 2**an_astronomical_exponent."""
    coefficient=Fraction(coefficient)
    if coefficient==0:return Fraction(0),0
    n,d=coefficient.numerator,coefficient.denominator
    tn=(abs(n)&-abs(n)).bit_length()-1;td=(d&-d).bit_length()-1
    return Fraction(n//(1<<tn),d//(1<<td)),exponent+tn-td


def binary_ref(g,ref,cache=None):
    cache={} if cache is None else cache
    if ref.node in cache:return cache[ref.node]
    n=g.nodes[ref.node];operation=n['operation']
    at=lambda i:binary_ref(g,true.source.FunctionRef(g,i),cache)
    if operation=='exact_rational':value=normalize_binary(Fraction(n['numerator'],n['denominator']))
    elif operation=='product':
        coefficient=Fraction(1);exponent=0
        for i in n['arguments']:
            m,e=at(i);coefficient*=m;exponent+=e
        value=normalize_binary(coefficient,exponent)
    elif operation=='negative':
        m,e=at(n['argument']);value=(-m,e)
    elif operation=='positive_quotient':
        m,e=at(n['numerator']);a,b=at(n['denominator']);value=normalize_binary(m/a,e-b)
    elif operation=='analytic_unary' and n['name']=='exp':
        arg=g.nodes[n['argument']]
        factors=arg['arguments'] if arg['operation']=='product' else [n['argument']]
        exponent=Fraction(1);logs=0
        for i in factors:
            factor=g.nodes[i]
            if factor['operation']=='exact_rational':exponent*=Fraction(factor['numerator'],factor['denominator'])
            elif factor['operation']=='analytic_unary' and factor['name']=='log':
                base=g.nodes[factor['argument']]
                if base['operation']!='exact_rational' or Fraction(base['numerator'],base['denominator'])!=2:
                    raise ValueError('Binary sc exponential requires the exact log(2)')
                logs+=1
            else:raise ValueError('Selected sc exponential is not an exact binary recipe')
        if logs!=1 or exponent.denominator!=1:raise ValueError('Integer binary sc exponent required')
        value=normalize_binary(1,int(exponent))
    else:raise ValueError('Selected sc must have its exact binary rational recipe')
    cache[ref.node]=value;return value


@dataclass(frozen=True,eq=False)
class OriginalTrueCutoffNumericPacket:
    scope:str='cut unrepaired original n1, common inner and Rin'


class CurrentOriginalN1CutoffNumericBridge:
    @source_precision
    def __init__(self,accepted_inner,accepted_cutoff,require_checked=True):
        if type(accepted_inner) is not inner.CurrentOriginalN1InnerNumericBridge or not accepted_inner.acceptance_loaded:
            raise ValueError('Checked true common-inner numerical bridge required')
        if type(accepted_cutoff) is not legacy.CurrentOriginalN1RadialCutoffExtension or not accepted_cutoff.acceptance_loaded:
            raise ValueError('Checked original cutoff own-field numerical owner required')
        self.inner,self.numeric=accepted_inner,accepted_cutoff
        self.cut=accepted_inner.full.cutoff;self.g=self.cut.g;self.family=copy.deepcopy(self.cut.family)
        self.c=accepted_cutoff.c
        self.bindings=full.own.encode(dict(fields=self.cut.fields,cumulative=self.cut.cumulative,
            endpoint=self.cut.endpoint,radial_axis=self.cut.radial.axis,geometry=self.cut.geometry))
        self._bindings=copy.deepcopy(self.bindings)
        self._owners=(self.inner,self.numeric,self.cut,self.g,self.c)
        self.source_identity=self.assert_source_identity()
        self._identity=op.snapshot(self.source_identity)
        self.hashes={}
        for owner in (accepted_inner,accepted_cutoff):
            for name,digest in owner.hashes.items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting cut n1 source pin '+name)
                self.hashes[name]=digest
        for module in (inner,legacy,true,legacy.parent,legacy.parent.parent):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):
                digest=sha(name)
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Changed cut bridge dependency '+name)
                self.hashes[name]=digest
        for name in (Path(__file__).name,Path(__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed true cutoff numerical input '+name)
        self._hashes=dict(self.hashes);self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=full.parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or receipt.get('source_family')!=self.family
                or any(receipt.get(k) for k in OPEN) or pins!=self.hashes or archive!=sha(NAME)
                or receipt.get('source_identity')!=op.snapshot(self.source_identity)):
                raise ValueError('Checked true cutoff own-field numerical bridge required')
            self.acceptance_loaded=True

    @source_precision
    def assert_source_identity(self):
        self.inner.assert_graph();self.numeric.assert_graph()
        if (self.numeric.regular is not self.inner.numeric or self.numeric.physical is not self.inner.numeric.physical
            or self.numeric.inputs is not self.inner.numeric.mean.inputs
            or self.numeric.second.regular is not self.numeric.regular
            or self.numeric.moments.regular is not self.numeric.regular
            or self.numeric.family!=self.family or self.inner.family!=self.family):
            raise ValueError('Cutoff numerical owner escaped the true common-inner source')
        self.inner.assert_source_identity()
        p=self.cut.parameters;sc=binary_ref(self.g,p['sc']);quarter=binary_ref(self.g,p['sc_quarter'])
        if quarter!=normalize_binary(sc[0],sc[1]-2):raise ValueError('True cutoff first phase is not the original sc/4')
        bits=self.inner.selection.sc._mpi_[0];sign,mantissa,exponent,_=bits
        selected=normalize_binary((-1)**sign*mantissa,exponent)
        if selected!=sc or self.inner.selection.phase_in._mpi_!=(self.inner.selection.sc/4)._mpi_:
            raise ValueError('Numerical selected phase and true rational sc differ')
        if getattr(self.numeric.selected_phase,'__func__',None) is not legacy.CurrentOriginalN1RadialCutoffExtension.selected_phase:
            raise ValueError('Legacy cutoff selected-phase source callback changed')
        if legacy.domain.sigma_jets is not true.own.c2.flat_source.sigma_jets:
            raise ValueError('Cutoff derivatives use another flat source')
        argument=self.cut.geometry['argument'].node
        Gamma=true.own.jets.complete.Gamma
        expected_value=Gamma.parent.active.current.ast_binding(Gamma.flat_source.stable_sigma)
        expected_derivative=true.own.c2.current.ast_binding(true.own.c2.flat_source.sigma_jets)
        values=[n for n in self.g.nodes if n['operation']=='exact_original_source_flat_sigma' and n['argument']==argument]
        derivatives={n['name']:n for n in self.g.nodes if n['operation']=='analytic_unary'
            and n.get('argument')==argument and n.get('name') in ('original_flat_sigma_prime','original_flat_sigma_second')}
        if (len(values)!=1 or values[0]['source_binding']!=expected_value
            or not values[0]['exact_function_not_its_interval_cap']
            or set(derivatives)!= {'original_flat_sigma_prime','original_flat_sigma_second'}
            or any(n['derivative_source']!=expected_derivative for n in derivatives.values())
            or derivatives['original_flat_sigma_second']['Taylor_to_ordinary_factor']!=2):
            raise ValueError('True/numeric sigma value and ordinary derivatives differ')
        for quantity in ('F1','Uz1','mean1','K1'):
            if any(ref!=self.g.zero for ref in true.c3.rows(self.cut.radial.axis[quantity][0])):
                raise ValueError('Original even axis value must be exactly zero: '+quantity)
        pressure=self.inner.inner.functions['core']['P1']
        if self.cut.fields['core']['pressure_P1']!=pressure:
            raise ValueError('Core cutoff pressure must keep the original regular pressure function')
        for ref in true.c3.rows(pressure):
            node=self.g.nodes[ref.node]
            if (node['operation']!='current_original_unique_regular_n1_component' or node['quantity']!='P1'
                or node['chart']!='core' or node['zero_axis_regular_power']!=2
                or not node['all_words_retained_in_analytic_limit']):
                raise ValueError('True regular n1 pressure must have exact zero axis data')
        receipt=full.parent.linear.read(legacy.RECEIPT)
        flags=('explicit_hb_inverse_inverse_squared_and_sigma_Taylor_second_coefficient_verified',
            'divergence_recovery_from_own_cutoff_axial_primitive_verified',
            'pressure_recovery_and_nonzero_fixed_leading_forcing_tail_signs_verified',
            'full_first_moment_matches_terminal_own_M',
            'remaining_two_leading_force_tails_and_repair_recursion_global_flags_not_admitted')
        if not receipt['all_passed'] or not all(receipt['checks'].get(k) for k in flags):
            raise ValueError('Actual cutoff own-history/signed-pressure proof required')
        if not self.cut.definition['mean_K_Q_V_pressure_recovered_from_own_history']:
            raise ValueError('True own history changed')
        if not self.inner.inner.definition['one_core_to_collar_path_no_seam_reset']:
            raise ValueError('True first inlet must keep the same regular core solution')
        tree=ast.parse(inner.inspect.getsource(legacy.natural.mean.actual_rows).strip())
        inlet=ast.parse('if ep(phase)[1]==0:\n return inlet').body[0]
        if not any(ast.dump(node)==ast.dump(inlet) for node in ast.walk(tree)):
            raise ValueError('Numerical first-zero leading rows must keep the actual core inlet')
        histories=[node for node in self.g.nodes if node['operation']=='definite_integral'
            and (node.get('actual_cutoff_own_mean') or node.get('actual_cutoff_own_pressure'))]
        if not histories or any(node['lower']!=self.g.zero.node
            or node['upper']!=self.inner.inner.coordinates['first'].node for node in histories):
            raise ValueError('True first-zero own mean/pressure histories must have equal zero endpoints')
        return dict(passed=True,same_checked_true_regular_function_and_live_numeric_source=True,
            selected_phase_exact_binary_rational_is_sc_over_four=True,
            same_original_sigma_value_callable_and_derivative_AST=True,
            ordinary_sigma_second_is_twice_Taylor_coefficient=True,
            original_positive_hb_and_inverse_powers_preserved=True,
            original_even_axis_F_U_mean_K_values_are_exact_zero=True,
            original_core_regular_P1_zero_axis_power_is_two=True,
            first_zero_keeps_same_regular_core_rho4_fields_and_zero_own_history_increments=True,
            signed_base_plus_cutoff_increment_retains_leading_N1p=True,
            own_M_mean_K_Q_V_P_reconstructed_with_full_history=True,
            Q_recovered_analytically_without_axis_radius_division=True,
            defining_fields_and_all_five_inner_moments_are_true_graph_refs=True,
            first_three_moments_have_full_velocity_support=True,
            two_pressure_moments_stop_at_Rin_and_need_exterior_sources=True,
            only_supported_derivative_rows_are_exported=True)

    def assert_graph(self):
        if ((self.inner,self.numeric,self.cut,self.g,self.c)!=self._owners or self.hashes!=self._hashes
            or self.family!=self.inner.family or self.family!=self.numeric.family
            or self.bindings!=self._bindings or op.snapshot(self.source_identity)!=self._identity
            or self.bindings!=full.own.encode(dict(fields=self.cut.fields,cumulative=self.cut.cumulative,
                endpoint=self.cut.endpoint,radial_axis=self.cut.radial.axis,geometry=self.cut.geometry))):
            raise ValueError('True cutoff numeric graph/metadata changed')
        self.inner.assert_graph();self.numeric.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate,chart='core',cells=4,reuse_packet=None,reuse_cumulative=None):
        self.assert_graph();c=self.c;value=c.mpf(coordinate);z=c.mpf(Z);lo,hi=ep(value)
        if not isinstance(coordinate,(str,int,float)) and value._mpi_[0]!=value._mpi_[1]:
            raise ValueError('Cutoff own cumulative endpoint must be one exact radial coordinate')
        if chart not in ('core','first') or lo<0 or hi>(4 if chart=='core' else 1):
            raise ValueError('Original core rho[0,4] or first fraction q[0,1] required')
        if isinstance(cells,bool) or not isinstance(cells,int) or not 1<=cells<=64:raise ValueError('Partition must be in [1,64]')
        packet=self.numeric.evaluate(Z,coordinate,chart,cells) if reuse_packet is None else reuse_packet
        view=self.numeric.report(packet);entry=self.numeric._packets[id(packet)]
        if (view['Z']._mpi_!=z._mpi_ or view['query']['coordinate']._mpi_!=value._mpi_
            or view['query']['chart']!=chart):raise ValueError('Cutoff packet belongs to another actual query')
        chi,pairs,physical,M,own,P,totals=entry[2]
        if view['radial_cutoff']['original_selected_phase']._mpi_!=self.inner.selection.phase_in._mpi_:
            raise ValueError('Actual legacy cutoff phase differs from the true sc/4 recipe')
        if entry[3][4] is not self.numeric.changes(z,coordinate if chart=='first' else '0',cells):
            raise ValueError('Cutoff packet uses another original history partition')
        regular_packet=entry[3][1];rentry=self.numeric.regular._packets[id(regular_packet)]
        W,dx,dr,recovered,recovered_dr,system=rentry[5];rho=system['rho'];algebra=W[0].algebra
        Zjet=algebra.lift(op.Jet.variable(c,z,2));d=1-Zjet*Zjet;L=1-Zjet*Zjet*self.numeric.physical.delta
        U=physical['Uz1']['value'];mean=own['mean1']
        Q=op.truncate((U*(Zjet*2)-mean*(Zjet*(1+self.numeric.physical.delta))-op.dz(mean)*d)/L,1)
        axis=ep(rho)[1]==0
        mean1=(physical['Uz1']['rho1']/2 if axis else (U-op.truncate(mean,1))/rho)
        mean2=(physical['Uz1']['rho2']/3 if axis else (physical['Uz1']['rho1']-mean1*2)/rho)
        radial_own=dict(mean1=[op.truncate(mean,1),mean1,mean2],
            K1=[op.truncate(own['K1'],1),mean1-physical['Uz1']['rho1'],mean2-physical['Uz1']['rho2']])
        mean_packet=self.numeric.regular.natural._packets[id(rentry[2])][2]
        rows=self.numeric.regular.mean._packets[id(mean_packet)][4][1]
        leadingF=algebra.shift(rows['Phi'][0]*legacy.weighted.relative_amplitude(algebra,system['gradient'],1,1),(0,1,0,0))
        Np=algebra.shift(system['g4_over_x']/2,(0,0,1,0))
        pressure_rho=op.truncate(algebra.shift(2*leadingF*physical['F1']['value']+Np,(0,0,0,-1)),1)
        cumulative_packet=None
        if chart=='first' and value._mpi_==c.mpf(1)._mpi_:
            cumulants=dict(totals)
        else:
            cumulative_packet=self.numeric.moments.evaluate(Z,coordinate,chart) if reuse_cumulative is None else reuse_cumulative
            cumulative_view=self.numeric.moments.report(cumulative_packet)
            if (cumulative_view['Z']._mpi_!=z._mpi_ or cumulative_view['query']['chart']!=chart
                or cumulative_view['query']['coordinate']._mpi_!=value._mpi_):
                raise ValueError('Cumulative packet belongs to another original point')
            base=self.numeric.moments._packets[id(cumulative_packet)][3][2]
            changes=entry[3][4][0]
            cumulants={name:base[name]+changes[name] for name in legacy.parent.MOMENTS}
        # These two accepted FTC identities use the very same own fields,
        # avoiding independent boxes for the zero-support mean history.
        cumulants.update(m11=op.truncate(M,1),m13=P)
        result=dict(source_family=self.family,Z=z,chart=chart,coordinate=value,
            original_phase=None if chart=='core' else value*self.inner.selection.phase_in,
            true_function_bindings=self.bindings,
            numerical_cutoff=view,
            physical_own_Q1_Z1=op.export(Q,1),
            physical_own_mean_K_rho0_rho1_rho2_Z1={k:[op.export(v,1) for v in vals] for k,vals in radial_own.items()},
            physical_own_pressure_rho1_Z1=op.export(pressure_rho,1),
            physical_own_M_rho1_Z1=op.export(algebra.shift(U,(0,0,0,-1)),1),
            physical_own_M_rho2_Z1=op.export(algebra.shift(physical['Uz1']['rho1'],(0,0,0,-1)),1),
            physical_five_inner_cumulative_Z1={name:op.export(v,1) for name,v in cumulants.items()},
            five_inner_cumulative_budgets={name:{'Z'+str(k):legacy.parent.moment_budget(v,k) for k in (0,1)} for name,v in cumulants.items()},
            full_velocity_support_names=['m11','m12','m14'],
            exterior_pressure_m13_m15_not_added=True,
            actual_cutoff_own_field_and_cumulative_query_executed=True,
            Q_axis_radius_division_not_used=True,
            even_velocity_rho2_and_supported_own_radial_rows_only=True,
            pressure_rho2_and_full_C3_C5_point_accuracy_not_claimed=True,
            repaired_reserved_fields_and_terminal_closure_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalTrueCutoffNumericPacket()
        state=(Q,radial_own,pressure_rho,cumulants)
        self._packets[id(token)]=(token,result,packet,cumulative_packet,state,op.snapshot((result,state)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalTrueCutoffNumericPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live true cutoff own-field packet required')
        self.numeric.report(entry[2])
        if entry[3] is not None:self.numeric.moments.report(entry[3])
        if op.snapshot((entry[1],entry[4]))!=entry[5]:raise ValueError('True cutoff numeric packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_inner,accepted_cutoff,reuse_packets=None,reuse_cumulative=None):
    began=time.monotonic();owner=CurrentOriginalN1CutoffNumericBridge(accepted_inner,accepted_cutoff,require_checked=False)
    packets={}
    for name,coordinate,chart in (('axis','0','core'),('first_inlet','0','first'),('transition','.625','first'),('zero','.75','first'),('Rin','1','first')):
        print('TRUE_CUT_OWN_CUMULATIVE_START',name,flush=True)
        packets[name]=owner.evaluate('.371',coordinate,chart,cells=4,
            reuse_packet=None if reuse_packets is None else reuse_packets[name],
            reuse_cumulative=None if reuse_cumulative is None else reuse_cumulative.get(name))
        print('TRUE_CUT_OWN_CUMULATIVE_COMPLETE',name,flush=True)
    raw=dict(source_family=owner.family,source_identity=owner.source_identity,true_function_bindings=owner.bindings,
        samples={name:owner.report(packet) for name,packet in packets.items()},input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,raw,packets


def restore_checked_archive():
    receipt=full.parent.linear.read(RECEIPT);raw=full.parent.linear.read(NAME)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
    if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
        or archive!=sha(NAME) or pins!=raw['input_hashes'] or receipt['source_family']!=raw['source_family']
        or receipt['source_identity']!=raw['source_identity']):raise ValueError('Checked true cutoff numeric archive required')
    for name,digest in receipt['input_hashes'].items():
        if sha(name)!=digest:raise ValueError('Changed cutoff numeric archive input '+name)
    return raw,receipt

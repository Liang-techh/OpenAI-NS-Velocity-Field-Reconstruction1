"""Exact native Rp source join for the copied selected numerical pulse.

This joins function bodies and source scalars, rather than overlapping boxes.
It does not yet transfer the numerical pressure error bounds to the true pulse.
"""
import ast
import copy
import gzip
import inspect
import json
from pathlib import Path
import textwrap
import time
from types import SimpleNamespace

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions as true
import lei_ren_part1_paper_compliant_current_original_Rp_selected_pulse as selected
import lei_ren_part1_paper_compliant_current_original_Rp_native_constants_check as native_check
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_coefficients_refined as coefficients
import lei_ren_part1_paper_compliant_axial_high_jets as high
import lei_ren_part1_paper_compliant_fifth_axial_jets as fifth
import lei_ren_part1_paper_compliant_power_inlet_C4 as power
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=true.HERE,true.PREFIX,true.sha
NAME=PREFIX+'current_original_n1_pulse_native_numeric_identity.json.gz'
RECEIPT=PREFIX+'current_original_n1_pulse_native_numeric_identity_check.json'
GATE='current_original_n1_true_pulse_native_numeric_source_identity_installed'
OPEN=tuple(dict.fromkeys(true.OPEN+(
    'current_original_true_pulse_pressure_C1_numeric_bridge_installed',
    'uniform_full_pressure_derivative_accuracy_certified',
    'positive_order_five_moment_repair_closed','temporal_recursion',
    'full_corrected_Navier_Stokes_solution')))


def assignment(fn,target):
    tree=ast.parse(textwrap.dedent(inspect.getsource(fn)))
    matches=[n for n in ast.walk(tree) if isinstance(n,ast.Assign)
        and any(ast.unparse(t)==target for t in n.targets)]
    if len(matches)!=1:raise ValueError('Unique original assignment required '+target)
    return copy.deepcopy(matches[0])


def execute_assignments(fn,targets,env,scope=None):
    """Execute only the named mathematical assignments from actual production."""
    tree=ast.Module(body=[assignment(fn,t) for t in targets],type_ignores=[])
    ast.fix_missing_locations(tree)
    namespace=dict(fn.__globals__);namespace.update(scope or {});namespace.update(env)
    exec(compile(tree,'<actual-native-pulse-source-assignments>','exec'),namespace)
    return namespace


def exact_production_rows():
    """Actual Taylor algebra under an exact scalar context, not a new formula."""
    z,U,C1,C2,C0,CE,F1,F2=sy.symbols('Z U C1 C2 C0 C_E F1 F2')
    c=SimpleNamespace(mpf=sy.sympify)
    k=dict(U=U,C1=C1,C2=C2,C0=C0,C_E=CE)
    owner=SimpleNamespace(ctx=c,constants=k)
    p=execute_assignments(power.CompliantPowerInletC4.incoming,
        ('c','Z','z','q','invq','k','u','m','e'),dict(self=owner,Z=z))
    h=execute_assignments(high.CompliantAxialHighJets.incoming,
        ('c','Z','g','h','k'),dict(self=owner,Z=z))
    # Execute the production return itself: it contains the moment/energy
    # scaling and does not call a pressure or future producer.
    fn=ast.parse(textwrap.dedent(inspect.getsource(high.CompliantAxialHighJets.incoming))).body[0]
    returns=[n for n in ast.walk(fn) if isinstance(n,ast.Return)]
    if len(returns)!=1:raise ValueError('Unique original incoming return required')
    h['incoming']=eval(compile(ast.Expression(returns[0].value),'<actual-incoming-return>','eval'),h)
    old=dict(incoming=h['incoming'])
    c5=execute_assignments(fifth.CompliantFifthAxialJets.select,
        ('incoming','k','moments','energy','rows'),
        dict(self=SimpleNamespace(fourth=owner),old=old,c=c,Z=z,
            b=SimpleNamespace(incoming_factor_caps=[F1,F2])),
        {'copy_jet':lambda ctx,jet:high.IntervalTaylor(ctx,list(jet.coefficients)),
         # append_fifth's actual body is used with the same exact copy lift.
         'append_fifth':exact_append_fifth()})
    functions=dict(u=U/(1+z*z),m1=C1*z*(1+z*z),m2=C2*z*(1+z*z),
        energy=C0+CE*z*z*(1+z*z)**2)
    actual=dict(u=p['u'],m1=p['m'][0],m2=p['m'][1],energy=p['e'])
    selected_rows=dict(m1=c5['moments'][0],m2=c5['moments'][1],energy=c5['energy'])
    for key,expr in functions.items():
        for n in range(6):
            want=sy.diff(expr,z,n)/sy.factorial(n)
            if sy.cancel(actual[key][n]-want)!=0:raise ValueError('Native inlet Taylor body differs '+key)
            if key!='u' and sy.cancel(selected_rows[key][n]-want)!=0:
                raise ValueError('Selected C5 incoming body differs '+key)
    for j,factor in enumerate((F1,F2)):
        for n in range(6):
            if sy.cancel(c5['rows'][j][n]-factor*c5['moments'][j][n])!=0:
                raise ValueError('Selected incoming factor not applied once')
    # Execute the refined consumer's actual assignment as well. It must use
    # raw incoming moments, rather than already factor-normalized rows.
    r=execute_assignments(coefficients.CurrentOriginalRpPulseCoefficientsRefined.evaluate,
        ('incoming','factors','rows'),dict(old={'incoming':{'moment_Taylor':c5['moments']}},
            c=c,amplitude=SimpleNamespace(incoming_factor_caps=[F1,F2]),
            self=SimpleNamespace(selected=SimpleNamespace(ctx=c))),
        {'rebox':lambda ctx,v:v})
    for j,factor in enumerate((F1,F2)):
        for n in range(6):
            if sy.cancel(r['rows'][j][n]-factor*c5['moments'][j][n])!=0:
                raise ValueError('Refined consumer incoming factor differs')
    return dict(passed=True,power_inlet_exact_Taylor_rows=24,selected_incoming_exact_Taylor_rows=18,
        selected_and_refined_once_scaled_row_identities=24,
        retained_energy_polynomial='C0+C_E*Z^2*(1+Z^2)^2',
        energy_Z4_Z6_not_discarded=True,actual_production_assignments_executed=True,
        ordinary_derivatives_divided_by_factorial=True,
        pressure_and_future_functions_not_substituted_by_this_proof=True)


def exact_append_fifth():
    tree=ast.parse(textwrap.dedent(inspect.getsource(fifth.append_fifth)))
    namespace=dict(fifth.append_fifth.__globals__)
    namespace['copy_jet']=lambda c,j:high.IntervalTaylor(c,list(j.coefficients))
    exec(compile(tree,'<actual-fifth-extension>','exec'),namespace)
    return namespace['append_fifth']


def scalar_expression(g,ref,bindings):
    cache=dict(bindings)
    def at(i):
        if i in cache:return cache[i]
        row=g.nodes[i];op=row['operation'];child=lambda j:at(j)
        if op=='exact_rational':v=sy.Rational(row['numerator'],row['denominator'])
        elif op=='sum':v=sum(map(child,row['arguments']),sy.Integer(0))
        elif op=='product':v=sy.prod(map(child,row['arguments']))
        elif op=='negative':v=-child(row['argument'])
        elif op=='positive_quotient':v=child(row['numerator'])/child(row['denominator'])
        elif op=='analytic_unary' and row['name'] in ('exp','log'):
            v=getattr(sy,row['name'])(child(row['argument']))
        else:raise ValueError('Unsupported exact scalar factor operation '+op)
        cache[i]=v;return v
    return at(ref.node)


class CurrentOriginalN1PulseNativeNumericIdentity:
    @source_precision
    def __init__(self,pulse_owner,selected_owner,require_checked=True):
        if type(pulse_owner) is not true.CurrentOriginalN1PulsePressureFunctions or not pulse_owner.acceptance_loaded:
            raise ValueError('Live checked true complete pulse functions required')
        if type(selected_owner) is not selected.CurrentOriginalRpSelectedPulse or not selected_owner.acceptance_loaded:
            raise ValueError('Live checked copied selected numerical pulse required')
        self.pulse,self.selected=pulse_owner,selected_owner
        self.g=pulse_owner.g;self.family=copy.deepcopy(pulse_owner.family)
        if selected_owner.family_record!=self.family:raise ValueError('Current native pulse families differ')
        self.hashes=dict(pulse_owner.hashes)
        for name,digest in selected_owner.hashes.items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting native source '+name)
            self.hashes[name]=digest
        for module in (true,selected,selected.inlet.exact):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        callables=(power.CompliantPowerInletC4.incoming,high.CompliantAxialHighJets.incoming,
            high.CompliantAxialHighJets.select,fifth.CompliantFifthAxialJets.select,fifth.append_fifth,
            coefficients.CurrentOriginalRpPulseCoefficientsRefined.evaluate,
            selected.selected.CurrentSelectedEnergySource.__init__,
            selected.CurrentOriginalRpSelectedPulse.__init__)
        self.callable_bindings=[coefficients.function_binding(fn) for fn in callables]
        for b in self.callable_bindings:self.hashes[b['source_file']]=b['source_sha256']
        for name in (Path(__file__).name,Path(__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed native pulse source '+name)
        self.proof=self.source_identity()
        self._proof=copy.deepcopy(self.proof);self._hashes=dict(self.hashes)
        self._bindings=copy.deepcopy(self.callable_bindings)
        self._constant_state=coefficients.record((selected_owner.constants,selected_owner.inlet.interval_constants,
            selected_owner.amplitude.incoming_factor_log_definitions,selected_owner.amplitude.incoming_factor_caps))
        self._owners=(self.pulse,self.selected,self.g)
        self.acceptance_loaded=False
        if require_checked:
            receipt=true.parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or pins!=self.hashes or archive!=sha(NAME) or receipt.get('source_family')!=self.family
                or receipt.get('source_identity')!=self.proof):
                raise ValueError('Checked true native numeric source join required')
            self.acceptance_loaded=True

    @source_precision
    def source_identity(self):
        self.pulse.assert_graph();self.selected.assert_graph()
        inlet=self.selected.inlet;exact=inlet.exact
        if type(exact) is not selected.inlet.exact.CurrentOriginalRpNativeConstants or not exact.acceptance_loaded:
            raise ValueError('Actual checked native scalar source required')
        receipt=true.parent.linear.read(selected.inlet.exact.RECEIPT)
        if (not receipt['actual_source_functions_identified_not_interval_enclosure_equality']
            or receipt['actual_old_buffer_to_original_native_coefficient_identities']!=7
            or receipt['actual_native_constructor_to_current_constant_function_identities']!=12):
            raise ValueError('Actual native scalar provenance proof required')
        before=len(self.g.nodes)
        native=true.own.jets.native_C3_functions(self.g)
        if true.own.jets.encode(native)!=true.own.jets.encode(self.pulse.own.jets.native) or len(self.g.nodes)!=before:
            raise ValueError('Actual C3 native function graph differs')
        q,constants=exact.current_constants();_,data=exact.frame.symbolic_frame()
        z=q.z;Q=1+z*z
        bodies=dict(u=constants['U']/Q,m1=constants['C1']*z*Q,m2=constants['C2']*z*Q,
            energy=constants['C0']+constants['C_E']*z*z*Q*Q)
        for key,body in bodies.items():
            if not native_check.zero(data[key]-body):raise ValueError('Current native scalar formula differs '+key)
        if any(z in value.free_symbols for value in constants.values()):raise ValueError('Native scalars depend on Z')
        if (type(self.selected.axial4) is not high.CompliantAxialHighJets
            or type(self.selected.fifth) is not fifth.CompliantFifthAxialJets
            or self.selected.axial4.incoming.__func__ is not high.CompliantAxialHighJets.incoming
            or self.selected.axial4.select.__func__ is not high.CompliantAxialHighJets.select
            or self.selected.fifth.select.__func__ is not fifth.CompliantFifthAxialJets.select
            or self.selected.fifth.fourth is not self.selected.axial4
            or self.selected.axial4.constants is not self.selected.constants):
            raise ValueError('Actual selected incoming callback or scalar aliases changed')
        c=self.selected.ctx
        for key in ('U','C1','C2','C0','C_E'):
            expected=c.mpf(selected.inlet.endpoints(inlet.interval_constants[key]))
            if self.selected.constants[key]._mpi_!=expected._mpi_:
                raise ValueError('Directed conversion changed native scalar '+key)
        production=exact_production_rows()
        mu,L=sy.symbols('mu logpref',positive=True)
        W=self.pulse.own.jets.complete.Gamma.pulse.functions
        factors=[scalar_expression(self.g,ref,{W['mu'].node:mu,W['log_common_prefactor'].node:L})
            for ref in W['incoming_factors']]
        env=dict(b=SimpleNamespace(mu=mu,pulse=SimpleNamespace(rows={'common_logpref':L})))
        actual=execute_assignments(selected.selected.CurrentSelectedEnergySource.__init__,
            ('b.incoming_factor_log_definitions',),env)['b'].incoming_factor_log_definitions
        for value,log in zip(factors,actual):
            if sy.simplify(value-sy.exp(log))!=0:raise ValueError('True and numeric incoming factor recipes differ')
        if any(selected.inlet.endpoints(v)[0]<=0 for v in self.selected.amplitude.incoming_factor_margins):
            raise ValueError('Positive actual incoming factor enclosure margin required')
        return dict(passed=True,exact_native_scalar_function_identities=4,
            native_scalar_source_receipt=selected.inlet.exact.RECEIPT,
            actual_buffer_and_constructor_identity_counts=[7,12],
            true_C3_native_function_handles=true.own.jets.encode(native),
            graph_replay_created_no_new_nodes=True,production_Taylor_proof=production,
            same_selected_live_callback_and_shared_constant_objects=True,
            exact_directed_context_conversion_for_native_scalars=True,
            exact_true_and_numeric_factor_recipe_identities=2,
            positive_factor_caps_enclose_functions_not_selected_as_constants=True,
            native_common_moments_already_normalized_by_Pstar=True,
            physical_u_and_mean_apply_Pstar_once=True,
            callable_bindings=self.callable_bindings,
            pressure_error_bounds_and_full_five_moments_not_yet_transferred=True)

    @source_precision
    def assert_graph(self):
        if (self._owners!=(self.pulse,self.selected,self.g) or self.hashes!=self._hashes
            or self.proof!=self._proof or self.callable_bindings!=self._bindings
            or coefficients.record((self.selected.constants,self.selected.inlet.interval_constants,
                self.selected.amplitude.incoming_factor_log_definitions,self.selected.amplitude.incoming_factor_caps))!=self._constant_state):
            raise ValueError('Actual native numeric source join changed')
        if (getattr(self.selected.axial4.incoming,'__func__',None) is not high.CompliantAxialHighJets.incoming
            or getattr(self.selected.axial4.select,'__func__',None) is not high.CompliantAxialHighJets.select
            or getattr(self.selected.fifth.select,'__func__',None) is not fifth.CompliantFifthAxialJets.select):
            raise ValueError('Actual native selected source callbacks changed')
        self.pulse.assert_graph();self.selected.assert_graph()


@source_precision
def run(pulse_owner,selected_owner):
    began=time.monotonic();owner=CurrentOriginalN1PulseNativeNumericIdentity(pulse_owner,selected_owner,False)
    raw=dict(source_family=owner.family,source_identity=owner.proof,input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_PULSE_NATIVE_NUMERIC_SOURCE_JOIN_DEFINED',flush=True)
    return owner,raw

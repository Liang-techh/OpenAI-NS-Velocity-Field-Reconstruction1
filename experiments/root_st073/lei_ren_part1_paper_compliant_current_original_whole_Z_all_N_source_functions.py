"""Original full-phase bridge functions parameterized by epsilon=1/N.

The fixed leading input is unchanged. Uniform N-scaled C1 drivers are
integrated from the defining zero correction inlet through the three
bridge charts. No fixed-N correction history is rescaled or reused.
"""
import ast
from dataclasses import dataclass
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_sharp_functions as current

HERE,PREFIX,sha,bind,ep = current.HERE,current.PREFIX,current.sha,current.bind,current.ep
bridge,backend = current.bridge,current.backend
CELLS,RATES,C0,Z,OPEN = current.CELLS,current.RATES,current.C0,current.Z,current.OPEN
relative = current.previous.relative
original = current.previous.switch.parameters.original
NAME = PREFIX+'current_original_whole_Z_all_N_source_functions.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_all_N_source_functions_check.json'
GATE = 'current_original_whole_Z_bridge_full_phase_all_N_Nscaled_C1_functions_and_R100_transport_installed'


@dataclass(frozen=True)
class Pair:
    value: object
    Z: object


def add(a,b):return Pair(a.value+b.value,a.Z+b.Z)
def scale(a,b):return Pair(a.value*b,a.Z*b)
def multiply(a,b):return Pair(a.value*b.value,a.Z*b.value+a.value*b.Z)
def square(a):return Pair(original.square(a.value),a.value*a.Z*2)
def record(a):return [a.value.record(),a.Z.record()]


def uniform_source_union(values):
    """Hull alternatives in one dominating absolute-log coordinate.

    Full-phase alternatives can have independent wide log offsets. Using
    one of those wide scales as the reference loses its own correlation
    under subtraction. A fixed upper log is only a range coordinate;
    every original signed coefficient and directed tail is retained.
    """
    if not values:raise ValueError('Nonempty same-source function alternatives required')
    first=values[0];c=first.ctx
    for row in values:first.coerce(row)
    nonzero=[row for row in values if not row.zero]
    if not nonzero:return first.scalar(0)
    logs=[]
    for row in nonzero:
        magnitude=max(abs(v) for v in ep(row.coefficient))
        logs.append((row,c.mpf(magnitude),row.scale.evaluate()+c.ln(c.mpf(magnitude))))
    upper=c.mpf(max(ep(log)[1] for _,_,log in logs))
    prior=backend.signed.density.local.prior
    reference=prior.FormalScale(first.scale.bases,offset=upper)
    coefficients=[c.mpf(0)] if len(nonzero)!=len(values) else []
    for row,magnitude,log in logs:
        coefficients.append((row.coefficient/magnitude)*row.bounded_exp(log-upper))
    hull=c.mpf((min(ep(v)[0] for v in coefficients),max(ep(v)[1] for v in coefficients)))
    return prior.ScaledEnclosure(reference,hull,first.ledger)


def epsilon_cover(c,N0,epsilon=None):
    if type(N0) is not int or N0<1:raise ValueError('Positive fixed leading-input prefix N0 required')
    upper=c.mpf(1)/N0
    value=c.mpf((0,ep(upper)[1])) if epsilon is None else c.mpf(epsilon)
    if ep(value)[0]<0 or ep(value)[1]>ep(upper)[1]:
        raise ValueError('epsilon must stay in the admitted [0,1/N0] parameter cover')
    return value


def normalized_drivers(E,V,primitives,N0,epsilon=None):
    """Exact N-scaled density functions, smoothly extended at epsilon=0.

    exprel(x)=integral_0^1 exp(t*x)dt. Its signed argument is bounded
    before materialization; the tiny epsilon factor remains formal.
    U_Z uses the derivative of the defining exponential, not a cap.
    """
    c=E.value.ctx; eps=epsilon_cover(c,N0,epsilon)
    A=Pair(primitives['A'],primitives['A_Z'])
    B=Pair(primitives['B_over_Pstar'],primitives['B_Z_over_Pstar'])
    for q in (E,V,A,B):
        E.value.coerce(q.value); E.value.coerce(q.Z)
    argument=A.value*eps
    finite=bridge.bounded_signed_exponent(argument)
    lo,hi=ep(finite)
    exprel=c.mpf((ep(c.exp(c.mpf(min(lo,0))))[0],ep(c.exp(c.mpf(max(hi,0))))[1]))
    exponential=c.exp(finite)
    U=Pair(E.value*A.value*exprel,
        E.Z*A.value*exprel+E.value*A.Z*exponential)
    theta=multiply(E,U); theta_square=scale(square(U),eps/2)
    drivers=dict(m=B,h=U,
        k=add(add(multiply(V,U),multiply(E,B)),scale(multiply(U,B),eps)),
        e=add(add(scale(multiply(V,B),2),scale(square(B),eps)),scale(add(theta,theta_square),-1)),
        p=add(theta,theta_square))
    return dict(drivers=drivers,epsilon=eps,normalized_swirl_increment=U,
        epsilon_times_A=argument,bounded_epsilon_times_A=finite,
        original_A_and_A_Z_retained=True,epsilon_Z_exact_zero=True,
        defining_U='E*A*exprel(epsilon*A)',
        defining_U_Z='E_Z*A*exprel(epsilon*A)+E*A_Z*exp(epsilon*A)',
        smooth_epsilon_zero_extension=True,all_quadratic_cross_terms_retained=True)


def density_binding():
    """Bind the original finite-N code and exact normalized identities."""
    density=current.backend.current.reference.phase.densities
    code=inspect.getsource(density.density_Z_kernels)
    tree=ast.parse(textwrap.dedent(code)).body[0]
    expected=("A=primitives['A']*factor", "A_Z=primitives['A_Z']*factor",
        "deltaE=E*increment", "deltaE_Z=E_Z*increment+E*exponential*A_Z",
        "deltaV=primitives['B_over_Pstar']*factor",
        "deltaV_Z=primitives['B_Z_over_Pstar']*factor")
    for statement in expected:
        wanted=ast.dump(ast.parse(statement).body[0])
        if sum(ast.dump(node)==wanted for node in ast.walk(tree))!=1:
            raise ValueError('Defining original density assignment changed: '+statement)
    eps,E,V,U,B=sy.symbols('epsilon E V U B')
    namespace=dict(E=E,V=V,deltaE=eps*U,deltaV=eps*B,theta_cross=E*eps*U,
        theta_square=(eps*U)**2,axial_square=(eps*B)**2)
    kernel_tree=ast.parse(textwrap.dedent(inspect.getsource(density.density.signed_density_kernels))).body[0]
    returned=next(node.value for node in ast.walk(kernel_tree) if isinstance(node,ast.Return))
    if not isinstance(returned,ast.Call) or not isinstance(returned.func,ast.Name) or returned.func.id!='dict':
        raise ValueError('Original five signed density dict required')
    expected_rows=dict(m=B,h=U,k=V*U+E*B+eps*U*B,
        e=2*V*B+eps*B**2-E*U-eps*U**2/2,p=E*U+eps*U**2/2)
    # The sole context scalar in the original expression is exactly .5.
    class Scalar(ast.NodeTransformer):
        def visit_Call(self,node):
            if ast.unparse(node.func)=='E.ctx.mpf' and len(node.args)==1:
                if not isinstance(node.args[0],ast.Constant) or node.args[0].value!='.5':
                    raise ValueError('Unexpected original density scalar')
                return ast.copy_location(ast.parse('half',mode='eval').body,node)
            return self.generic_visit(node)
    namespace['half']=sy.Rational(1,2)
    matched=[]
    for row in returned.keywords:
        expression=Scalar().visit(row.value)
        actual=eval(compile(ast.fix_missing_locations(ast.Expression(expression)),'original_density_identity','eval'),
            {'__builtins__':{}},namespace)
        if row.arg not in expected_rows or sy.expand(actual/eps-expected_rows[row.arg])!=0:
            raise ValueError('N-scaled source identity failed: '+str(row.arg))
        matched.append(row.arg)
    if set(matched)!=set(RATES):raise ValueError('All original five density identities required')
    modules=(Path(density.__file__).name,Path(density.density.__file__).name)
    return dict(original_density_module_hashes={name:sha(name) for name in modules},
        original_density_Z_AST_sha256=hashlib.sha256(ast.dump(tree).encode()).hexdigest(),
        original_five_density_normalization_symbolic_identities=matched,
        original_velocity_increments='deltaE=epsilon*U, deltaV=epsilon*B',
        no_saved_fixed_N_total_multiplied_by_N=True)


class WholeZAllNBridgeSourceFunctions:
    def __init__(self,dps=500):
        self.owner=current.WholeZSharpBridgeFunctions(dps)
        self.c,self.identity,self.N0=self.owner.c,self.owner.identity,self.owner.N
        receipt=json.loads((HERE/current.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(current.GATE) \
                or receipt['source_family']!=self.identity or receipt['candidate_N']!=self.N0:
            raise ValueError('Checked same original whole-Z leading input required')
        self.hashes=dict(self.owner.hashes)
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed original input prerequisite: '+name)
            bind(self.hashes,name,digest)
        bind(self.hashes,current.RECEIPT,sha(current.RECEIPT))
        self.binding=density_binding()
        for name,digest in self.binding['original_density_module_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.cache={}

    def primitive_cell(self,ends,chart,left,right=None):
        key=tuple(ends)
        if key not in CELLS or chart not in bridge.CHARTS:
            raise ValueError('Admitted original whole-Z bridge domain required')
        l=self.owner.coordinate(chart,left);r=self.owner.coordinate(chart,left if right is None else right)
        cachekey=(key,chart,l,r)
        if cachekey in self.cache:return self.cache[cachekey]
        packet=self.owner.source_query(key,chart,l,r);op=self.owner.owner(key)
        if packet['source_identity']!=self.identity or packet['candidate_N']!=self.N0:
            raise ValueError('Same fixed leading-input family and prefix N0 required')
        qrows=packet['original_q_C0_Z']
        source=dict(q=qrows[C0],roots=packet['original_roots'],
            original_C_B_E_C0_Z=packet['original_C_B_E_C0_Z'],
            original_eta_log=packet['original_eta_log'],
            original_full_source_quotients=packet['original_full_source_quotients'])
        if packet['exact_common_P0_axial5'] is not op.P0 or not packet['actual_phase_Z_exact_zero']:
            raise ValueError('Same independent P0 and phase-held Z identity required')
        roots=source['roots'];f=op.flow
        relative.same_source(f,[v for row in roots.values() for v in row.values()]+list(qrows.values()))
        _,branches,empty=backend.signed.signed_u_branches(source,self.owner.source.dstar_log)
        alternatives=[]
        def phase_piece(a,b,depth=0):
            phi=self.c.mpf((ep(self.c.mpf(a.numerator)/a.denominator)[0],
                ep(self.c.mpf(b.numerator)/b.denominator)[1]))
            got=[current.sharp_inverse_branch(source,qrows,self.owner.source.dstar_log,phi,branch) for branch in branches]
            if any(row['values'] is None for row in got):
                if depth>=5:raise ArithmeticError('Full original phase cover needs further source/phase refinement')
                middle=(a+b)/2;phase_piece(a,middle,depth+1);phase_piece(middle,b,depth+1);return
            alternatives.extend(dict(exact_phase_left=[a.numerator,a.denominator],
                exact_phase_right=[b.numerator,b.denominator],full_phase_parameter_box=phi,
                branch=branch['name'],primitive_C0_Z_phi=row['values'],original_inverse_proof=row['record'])
                for branch,row in zip(branches,got))
        from fractions import Fraction
        for i in range(4):phase_piece(Fraction(i,4),Fraction(i+1,4))
        E,V=packet['original_generic_source']['common_velocity_E_axial5'],packet['original_generic_source']['common_velocity_V_axial5']
        result=dict(source_identity=self.identity,source_input_fixed_prefix_N0=self.N0,
            exact_Z_cell=list(key),actual_chart=chart,exact_left=list(l),exact_right=list(r),
            exact_common_P0_axial5=op.P0,original_source_packet=packet,
            original_E_V_C0_Z=dict(E=Pair(*E[:2]),V=Pair(*V[:2])),
            full_phase_inverse_function_alternatives=alternatives,empty_signed_branches=empty,
            explicit_full_phase_parameter_cover=[0,1],
            full_cover_contains_original_frac_N_radius_for_every_integer_N=True,
            phase_is_a_range_parameter_not_a_selected_spatial_angle=True,
            phase_Z_exact_zero=True,original_roots_q_velocity_P0_widths_unchanged=True)
        self.cache[cachekey]=result
        return result

    def query(self,ends,chart,left,right=None,epsilon=None):
        source=self.primitive_cell(ends,chart,left,right)
        velocity=source['original_E_V_C0_Z'];eps=epsilon_cover(self.c,self.N0,epsilon)
        rows=[normalized_drivers(velocity['E'],velocity['V'],row['primitive_C0_Z_phi'],self.N0,eps)
            for row in source['full_phase_inverse_function_alternatives']]
        union=uniform_source_union
        drivers={name:Pair(union([row['drivers'][name].value for row in rows]),
            union([row['drivers'][name].Z for row in rows])) for name in RATES}
        return dict(source=source,epsilon=eps,N_scaled_density_C0_Z=drivers,
            alternative_normalized_driver_functions=rows,uniform_for_all_integer_N_ge_N0=True,
            phase_Z_zero_and_genuine_source_Z_product_rules=True,
            epsilon_zero_is_smooth_limit_not_a_selected_finite_N=True)

    def at_N(self,ends,chart,left,right=None,N=None):
        if type(N) is not int or N<self.N0:raise ValueError('Explicit integer N>=fixed input N0 required')
        return self.query(ends,chart,left,right,epsilon=self.c.mpf(1)/N)

    def transport(self,ends):
        op=self.owner.owner(ends);f=op.flow
        inlet=self.query(ends,'first_micro','inlet')
        if not inlet['source']['original_source_packet']['source_owned_collar_q_C0_Z_exact_zero'] \
                or any(not q.value.zero or not q.Z.zero for q in inlet['N_scaled_density_C0_Z'].values()):
            raise ValueError('Defining unchanged-left correction inlet must be exact zero for every N')
        incoming={name:Pair(f.scalar(0),f.scalar(0)) for name in RATES}
        windows=[]
        for chart,partition,_ in self.owner.partitions():
            local={name:Pair(f.scalar(0),f.scalar(0)) for name in RATES};cells=[]
            for left,right in zip(partition,partition[1:]):
                query=self.query(ends,chart,left,right);contributions={};weights={}
                for name,rate in RATES.items():
                    mass,decay,suffix=self.owner.weights(ends,chart,left,right,rate)
                    contribution=scale(query['N_scaled_density_C0_Z'][name],mass*suffix)
                    local[name]=add(local[name],contribution);contributions[name]=contribution
                    weights[name]=dict(original_mass=mass,cell_decay=decay,suffix_decay=suffix)
                cells.append(dict(exact_left=list(self.owner.coordinate(chart,left)),
                    exact_right=list(self.owner.coordinate(chart,right)),
                    full_phase_inverse_function_alternatives=query['source']['full_phase_inverse_function_alternatives'],
                    original_E_V_C0_Z={k:record(v) for k,v in query['source']['original_E_V_C0_Z'].items()},
                    epsilon=query['epsilon'],N_scaled_density_C0_Z={k:record(v) for k,v in query['N_scaled_density_C0_Z'].items()},
                    actual_N_scaled_contribution_C0_Z={k:record(v) for k,v in contributions.items()},
                    original_own_rate_weights=weights,
                    physical_Jacobian_applied_once=True,all_N_phase_cover_not_fixed_N_phase_reuse=True))
            memory={name:self.owner.weights(ends,chart,partition[0],partition[-1],rate)[1] for name,rate in RATES.items()}
            outgoing={name:add(scale(incoming[name],memory[name]),local[name]) for name in RATES}
            windows.append(dict(actual_chart=chart,actual_all_N_source_cells=cells,
                actual_N_scaled_incoming_C0_Z={k:record(v) for k,v in incoming.items()},
                original_incoming_own_rate_memory=memory,
                actual_N_scaled_local_C0_Z={k:record(v) for k,v in local.items()},
                actual_N_scaled_outgoing_C0_Z={k:record(v) for k,v in outgoing.items()}))
            incoming=outgoing
            print('Actual all-N full-phase bridge transport: '+str(ends)+' '+chart,flush=True)
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),fixed_leading_input_N0=self.N0,
            exact_common_P0_axial5=op.P0,epsilon=epsilon_cover(self.c,self.N0),
            actual_N_scaled_zero_inlet_C0_Z={name:[f.scalar(0),f.scalar(0)] for name in RATES},
            actual_all_N_bridge_windows=windows,
            actual_uniform_N_scaled_R100_correction_C0_Z={k:record(v) for k,v in incoming.items()},
            source_owned_zero_inlet_proved_for_every_N=True,
            original_positive_widths_and_rate_zero_pressure_memory_retained=True,
            fixed_N_histories_not_rescaled_or_hydrated=True,
            actual_uniform_N_scaled_Rc_targets_installed=False,
            actual_global_frequency_admitted=False,**dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZAllNBridgeSourceFunctions()
        rows=[owner.transport(ends) for ends in CELLS]
        report=dict(**{GATE:True},source_family=owner.identity,fixed_leading_input_N0=owner.N0,
            exact_Z_partition=CELLS,all_N_source_density_binding=owner.binding,
            actual_four_Z_all_N_bridge_transports=rows,
            actual_full_phase_all_N_bridge_driver_functions_installed=True,
            actual_uniform_N_scaled_Rc_targets_installed=False,
            actual_all_regions_all_N_adapter_installed=False,
            actual_N_squared_averaging_cancellation_proved=False,
            actual_terminal_controls_installed=False,actual_global_frequency_admitted=False,
            physical_original_exterior_five_targets_closed=False,**dict.fromkeys(OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual unchanged whole-Z bridge source, full phase parameter and epsilon=1/N '
                'family of genuine N-scaled C1 density functions. Original positive kernel measures '
                'rebuild zero-inlet normalized corrections through R100 for every N>=N0. '
                'No downstream uniform Rc target, N^-2 cancellation, global frequency, solved '
                'control, heat/stress, temporal recursion or full NS admission.')
        def serialized(value):
            if isinstance(value,Pair):return record(value)
            if isinstance(value,dict):return {k:serialized(v) for k,v in value.items()}
            if isinstance(value,(list,tuple)):return [serialized(v) for v in value]
            return current.serialized(value)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(current.encode(serialized(report)),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()

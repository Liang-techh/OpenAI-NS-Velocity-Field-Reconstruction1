"""Whole original pulse entrance with unchanged five-moment/pressure source.

This companion extends the admitted main source exporter to xi[0,.02].
Only its coordinate guard and domain metadata change. Original nonzero
incoming histories, selected energy cancellation and exact source logs remain.
"""
import ast
import copy
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 as main
from lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 import (
    CompliantPulseMainExitSimilarityC4, PREFIX, source_precision)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import gp_jets
from lei_ren_part1_paper_compliant_axial_pulse_field import gp_energy, gp as original_gp
import lei_ren_part1_paper_compliant_flat_pulse_derivatives as flat_source
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent
DOMAIN=dict(entrance_xi=['0','.02'], entrance_y='[0,.02/mu]', Z=[-1,1],
    ordinary_derivative='d_y=d_logR=mu*d_xi',
    angular_normalization='Mtheta/(sqrt(2)*R^(3/2)*Utheta)')
FALSE_FLAGS=('pulse_entrance_physical_decomposition_constructed','pulse_entrance_cone_certified',
    'pulse_entrance_upstream_completed_physical_interface_certified',
    'whole_outer_cone_certified','completed_full_tensor_cone_certified',
    'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
    'physical_energy_integral_certified','full_background_NS_validation',
    'temporal_recursion','production_exact_point_parameters_selected')


def compiled_entrance_exporter():
    """Replay the admitted entire exporter; change only its domain guard."""
    asts=SourceAST()
    original=asts.method('pulse_main_exit_similarity_C4','main_exit')
    fn=copy.deepcopy(original)
    fn.decorator_list=[]
    expected=ast.parse("endpoints(xi)[0]<endpoints(c.mpf('.02'))[0] or endpoints(xi)[1]>11",mode='eval').body
    matches=[node for node in ast.walk(fn) if isinstance(node,ast.If)
        and ast.dump(node.test)==ast.dump(expected)]
    if len(matches)!=1:raise ValueError('Original coordinate guard changed')
    guard=matches[0]
    guard.test=ast.parse("endpoints(xi)[0]<0 or endpoints(xi)[1]>endpoints(c.mpf('.02'))[1]",mode='eval').body
    # Leave the original exception body and all source arithmetic unchanged.
    fn.name='entrance_source'
    restored=copy.deepcopy(fn)
    restored.name=original.name
    restored.decorator_list=copy.deepcopy(original.decorator_list)
    next(node for node in ast.walk(restored) if isinstance(node,ast.If)
        and ast.dump(node.test)==ast.dump(guard.test)).test=copy.deepcopy(expected)
    if ast.dump(restored)!=ast.dump(original):
        raise ValueError('Entrance exporter changed source arithmetic')
    env=dict(vars(main));env['DOMAIN']=DOMAIN
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<unchanged original pulse source on entrance>', 'exec'),env)
    return env['entrance_source'],dict(
        entire_original_exporter_AST_unchanged_except_coordinate_guard=True,
        same_original_moment_pressure_velocity_and_stress_algorithms=True,
        input_hashes=asts.hashes)


def forward_entrance_energy_rows(c,mu,xi,ap,ein,partial_energy,Bh):
    """Original native forward energy, anchored at the actual nonzero inlet."""
    weighted=(ein*mu+ap*ap*partial_energy-main.decay_integral(c,2,xi)/2)*c.exp(2*xi)
    rows=[ein if xi==0 else weighted/mu]
    square=main.product_rows(Bh,Bh)
    for j in range(4):rows.append(square[j]-(c.mpf('.5') if j==0 else 0)+rows[j]*(2*mu))
    return rows


def entrance_source_proof(records):
    """Source-functional inlet identities, not overlap of interval values."""
    asts=SourceAST();checks={}
    def zero(name,value):
        if s.cancel(s.expand(s.expand_power_exp(value)))!=0:
            raise ArithmeticError('Entrance source identity failed: '+name)
        checks[name]=True
    mu,delta,z,xi=s.symbols('mu delta Z xi',real=True)
    U,Pin,Xp,D2,Rp,Pstar=s.symbols('U Pin Xp exact_end_square Rp Pstar',positive=True)
    C=1/(1+z*z);r=1-mu;p=1+2*mu
    ap=s.Function('same_selected_ap')(z)
    incoming=[s.Function('same_incoming_m'+str(i))(z) for i in (1,2)]
    ein=s.Function('same_incoming_energy')(z)
    future=s.Function('same_complete_future_half')(z)
    loss=s.Function('same_selected_beta_energy')(z)
    terminalP=s.Function('same_absolute_pressure_at_Rv_over_Bv_squared')(z)
    K=( -mu*ein+(1-s.exp(-26))/4+mu*s.exp(-26)*(future-D2*loss) )/ap**2
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp,expm1=lambda v:s.exp(v)-1)
    env=dict(math=math,mp=SimpleNamespace(mpf=lambda v:s.Rational(str(v))),
        gp_jets=lambda ctx,coord:[s.Integer(0)]*5,
        axial_derivative=lambda v:s.diff(v,z))
    for stem,name in (('collar_Gamma_C4','product_rows'),('collar_stress_C3','shifted_rows'),
        ('pulse_end_stress_C3','pulse_coefficients'),
        ('pulse_end_physical_C2','pulse_velocity_rows'),
        ('pulse_main_exit_similarity_C4','split_main_exit_stress')):
        asts.replay(stem,name,env)
    shape_fn=asts.replay('pulse_main_exit_similarity_C4','main_exit_shapes',env)
    shape=shape_fn(c,mu,delta,z,C,Xp,ap,incoming,[s.Integer(0)]*2,K,
        future,loss,terminalP,s.Integer(0))
    full_energy=[a-D2*b for a,b in zip(shape['e0'],shape['J'])]
    expected_e=[ein]
    for j in range(4):expected_e.append(2*mu*expected_e[j]-(s.Rational(1,2) if j==0 else 0))
    for j in range(5):
        zero('actual_selected_energy_inlet_ordinary_row'+str(j),full_energy[j]-expected_e[j])
        zero('actual_flat_entrance_axial_input_row'+str(j),shape['Bh'][j])
        for i,key in enumerate(('mi','ni'),1):
            rate=s.Rational(1,2)-i*mu
            zero('actual_nonzero_incoming_m'+str(i)+'_ordinary_row'+str(j),
                shape[key][j]-incoming[i-1]*(-rate)**j)
    forwardenv=dict(main=SimpleNamespace(decay_integral=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate,
        product_rows=env['product_rows']))
    forward_fn=asts.replay('pulse_entrance_similarity_C4','forward_entrance_energy_rows',forwardenv)
    forward_rows=forward_fn(c,mu,s.Integer(0),ap,ein,s.Integer(0),shape['Bh'])
    for j in range(5):
        zero('actual_anchored_forward_energy_inlet_row'+str(j),forward_rows[j]-expected_e[j])
    asts.expression('axial_pulse_field','main','emu',
        wanted='(e0*self.mu+ap*ap*K-decay_integral(c,2,xi)/2)*c.exp(2*xi)')
    asts.expression('pulse_entrance_similarity_C4','forward_entrance_energy_rows','weighted',
        wanted='(ein*mu+ap*ap*partial_energy-main.decay_integral(c,2,xi)/2)*c.exp(2*xi)')
    asts.method('axial_pulse_field','gp_energy')
    # Both directed shape providers enclose the SAME defining primitive.
    if flat_source.original_gp is not original_gp or main.gp is not original_gp:
        raise ValueError('Entrance scalar and derivative gp sources differ')
    checks['actual_scalar_and_derivative_providers_share_original_gp_object']=True
    asts.expression('flat_pulse_derivatives','_gp_piece','primitive',
        wanted="original_gp(c,xi)['value'] if hi<=10 else None")
    asts.expression('flat_pulse_derivatives','_gp_piece','entrance',wanted='sigma_jets(c,50*xi)')
    asts.expression('axial_pulse_field','gp','derivative',wanted='sigma_enclosure(c,50*xi)*cut-P*slope')
    v,arg=s.symbols('original_entrance_variable sigma_argument',real=True)
    sigma=s.Function('same_original_sigma')
    primitive=s.Integral(sigma(50*v),(v,0,xi))
    sigma_rows=[s.diff(sigma(arg),arg,j).subs(arg,50*xi)/s.factorial(j) for j in range(5)]
    flatfn=asts.method('flat_pulse_derivatives','_gp_piece')
    appends=[node for node in ast.walk(flatfn) if isinstance(node,ast.Call)
        and ast.unparse(node.func)=='coefficients.append']
    if len(appends)!=1:raise ValueError('Original gp derivative recurrence changed')
    asts.bindings['flat_pulse_derivatives._gp_piece.actual_primitive_derivative_append']=True
    coefficients=[primitive]
    for n in range(1,5):
        ns=dict(coefficients=coefficients,entrance=sigma_rows,n=n)
        exec(compile(ast.fix_missing_locations(ast.Module(body=[ast.Expr(value=copy.deepcopy(appends[0]))],
            type_ignores=[])),'<actual original gp derivative recurrence>','exec'),ns)
        zero('actual_entrance_primitive_derivative_Taylor_order'+str(n),
            coefficients[n]-s.diff(primitive,xi,n)/s.factorial(n))
    # Cutoff symmetry gives integral_0^.02 sigma(50v)dv=.01.
    odds=-1/arg**2+1/(1-arg)**2
    zero('actual_sigma_log_odds_reflection_antisymmetry',odds.subs(arg,1-arg)+odds)
    E=s.symbols('positive_sigma_odds',positive=True)
    zero('actual_sigma_reflection_complements',E/(1+E)+1/(1+E)-1)
    asts.expression('axial_pulse_field','point','odds',wanted='-1/y**2+1/(1-y)**2')
    sigma_context=MPIntervalContext();sigma_context.dps=100
    right_sigma=flat_source.sigma_jets(sigma_context,1)
    for j in range(5):
        if endpoints(right_sigma[j])!=((1,1) if j==0 else (0,0)):
            raise ValueError('Original sigma right endpoint not flat')
        checks['actual_xi_point02_gp_derivative'+str(j)+'_from_flat_sigma_and_symmetry']=True
    # Positive forward and future quadratures are partitions of the same
    # exact energy integrand. Derivative plus one anchor proves additivity.
    for target,wanted in (('a','xi*i/cells'),('b','xi*(i+1)/cells')):
        asts.expression('axial_pulse_field','gp_energy',target,wanted=wanted)
    asts.expression('pulse_main_exit_similarity_C4','partial_future_energy','length',wanted='11-xi')
    asts.expression('axial_pulse_field','gp_energy','shape',
        wanted="gp(c,c.mpf([endpoints(a)[0],endpoints(b)[1]]))['value']")
    asts.expression('axial_pulse_field','gp_energy','result',augmented=True,
        wanted='c.exp(-2*a)*decay_integral(c,2,xi/cells)*shape**2')
    asts.expression('pulse_main_exit_similarity_C4','partial_future_energy','value',
        wanted="gp(c,c.mpf([endpoints(a)[0],endpoints(b)[1]]))['value']")
    asts.expression('pulse_main_exit_similarity_C4','partial_future_energy','total',augmented=True,
        wanted='c.exp(-2*a)*decay_integral(c,2,length/cells)*value**2')
    energy_integrand=s.exp(-2*v)*s.Function('same_original_gp')(v)**2
    part=s.Integral(energy_integrand,(v,0,xi));remaining=s.Integral(energy_integrand,(v,xi,11))
    complete=s.Integral(energy_integrand,(v,0,11))
    zero('actual_forward_future_energy_partition_FTC_derivative',s.diff(part+remaining-complete,xi))
    zero('actual_forward_future_energy_partition_initial_anchor',(part+remaining-complete).subs(xi,0).doit())
    I=s.symbols('actual_partial_energy',real=True)
    forward_native=asts.evaluate(asts.expression('axial_pulse_field','main','emu'),
        dict(e0=ein,self=SimpleNamespace(mu=mu),ap=ap,K=I,
            decay_integral=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate,c=c,xi=xi))/mu
    distance=13-xi
    actual_back=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','e0'),
        dict(future=future,K=s.exp(-2*distance),c=c,distance=distance,mu=mu,ap=ap,
            xi=xi,future_energy=K-I))[0]-D2*loss*s.exp(-2*distance)
    zero('actual_selected_forward_and_backward_energy_same_entire_entrance_source',actual_back-forward_native)
    for key in ('ml','nl'):
        for j in range(5):zero('actual_only_local_'+key+'_inlet_row'+str(j)+'_zero',shape[key][j])
    # Bind the original flat source rather than assume an arbitrary input is flat.
    ctx=MPIntervalContext();ctx.dps=100
    flat=gp_jets(ctx,0)
    for j in range(5):
        if endpoints(flat[j])!=(0,0):raise ArithmeticError('Actual gp inlet not flat')
        checks['actual_original_gp_inlet_flat_order'+str(j)]=True
    asts.method('flat_pulse_derivatives','gp_jets')
    for label,key in (('radial','radial'),('axial','axial')):
        for j in range(5):
            zero('actual_local_velocity_'+label+'_inlet_row'+str(j)+'_zero',shape['velocity_local'][key][j])
    # All incoming radial velocity derivatives survive; establish their
    # original source formula before checking equality to upstream power.
    L=1-delta*z*z;d=1-z*z
    radial=(-(1-delta)*z*C*incoming[0]-d*s.diff(C*incoming[0],z))/L
    for j in range(5):
        zero('actual_nonzero_radial_incoming_inlet_row'+str(j),
            shape['velocity_incoming']['radial'][j]-radial*(-s.Rational(1,2))**j)
    # Reconstruct the ACTUAL exporter moment sums at xi=0. This binds the
    # five inlet comparisons to its source expressions and exact prefactors.
    owner=SimpleNamespace(logRp=s.log(Rp),logP=s.log(Pstar),logU=s.log(U),
        Pin=Pin,Xp=Xp,incoming=incoming,J0=loss,mu=mu)
    c.ln=s.log
    exportenv=dict(c=c,self=owner,C=C,r=r,p=p,xi=s.Integer(0),zero=s.Integer(0),
        rows=shape,shifted_rows=env['shifted_rows'],logR=s.log(Rp),
        logB=dict(logPstar=s.log(Pstar),actual_log_inlet_U=s.log(U),inverse_mu=0,finite=0),
        logH=s.Integer(0),logD0=s.log(D2)/2,
        extras=dict(one=0,incoming1=0,incoming2=0,Q=-13*p/mu,end_square=s.log(D2)))
    exportenv['logs']=asts.replay('pulse_main_exit_similarity_C4','logs',exportenv)
    exportenv['moment']=lambda parts,ordinary:s.exp(sum(parts.values()))*ordinary[0]
    for name in ('raw_in1','raw_in2','angular_memory','energy_constant','pinpart','swirllimit'):
        exportenv[name]=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit',name),exportenv)
    exported=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit','moments'),exportenv)
    raw={name:sum(parts.values()) for name,parts in exported.items()}
    expected=dict(theta=s.sqrt(2)*Rp**s.Rational(3,2)*Pstar*U*C*Xp,
        z=Rp*Pstar*U*C*incoming[0],
        theta_z=s.sqrt(2)*Rp**s.Rational(3,2)*(Pstar*U*C)**2*incoming[1],
        z_theta=Rp*(Pstar*U*C)**2*ein,p=Pstar**2*C*C*Pin)
    for name in raw:
        for n in range(6):
            zero('actual_five_raw_inlet_'+name+'_axial'+str(n),s.diff(raw[name]-expected[name],z,n))
    X=1/r+(Xp-1/r)*s.exp(-r*xi/mu)
    zero('actual_signed_angular_inlet_not_reset',X.subs(xi,0)-Xp)
    # Same original absolute datum: terminal pressure is the admitted FTC
    # transport of the inlet primitive plus the canonical analytic datum.
    datum=s.Function('same_analytic_P0_over_Pstar_squared')(z)
    original_terminal=(Pin*C*C+U*U*C*C*(1-s.exp(-13*p/mu))/(2*p)+datum)/(
        U*U*s.exp(-13*p/mu))
    pulse_pressure=-(U*C)**2/(2*p)+(U*U*s.exp(-13*p/mu))*(terminalP+C*C/(2*p))
    for n in range(6):
        zero('actual_absolute_pressure_inlet_same_original_datum_axial'+str(n),
            s.diff(pulse_pressure.subs(terminalP,original_terminal)-(Pin*C*C+datum),z,n))
    zero('actual_raw_Mp_inlet_no_pressure_tail_added',raw['p']/Pstar**2-Pin*C*C)
    for target,wanted in (('self.ap',"decode_jet(self.ctx,fifth['selected_ap_Taylor'])"),
        ('self.incoming',"[decode_jet(self.ctx,v) for v in fifth['incoming']['moment_Taylor']]"),
        ('self.incoming_energy',"decode_jet(self.ctx,fifth['incoming']['energy_Taylor'])")):
        asts.expression('pulse_main_exit_similarity_C4','__init__',target,wanted=wanted)
    gates=records['power_inlet_C4_check']['exact_functional_production_and_join_identities']
    for flag in ('canonical_incoming_m1','canonical_incoming_m2','canonical_incoming_energy',
        'canonical_Xp_constant','canonical_Mp_constant',
        'actual_O2_O3_source_chain_has_exact_canonical_whole_Z_shapes',
        'same_source_primitive_ODEs_identify_y_jets_through4',
        'original_flat_B_entrance_jets_are_exact_zero',
        'paper_3_9_and_physical_prefactors_preserve_join'):
        if not gates[flag]:raise ValueError('Canonical original power inlet missing: '+flag)
        checks['consumed_'+flag]=True
    pressure_gates=records['pulse_end_flatten_join']['source_endpoint_binding']['identities']
    for flag in ('actual_pulse_P0_getter_is_same_canonical_flatten_pressure_over_C0_squared',
        'consumed_same_absolute_pressure_identified_by_original_FTC_and_power_datum'):
        if not pressure_gates[flag]:raise ValueError('Original pressure source not identified: '+flag)
        checks['consumed_'+flag]=True
    for flag in ('exact_uncapped_selected_sources_used','exact_functional_main_gap_and_gap_end_identities_certified'):
        if not records['pulse_interface_certificate'][flag]:raise ValueError('Exact selected source not available')
        checks['consumed_'+flag]=True
    return dict(identities=checks,input_hashes=asts.hashes,
        actual_source_AST_bindings=asts.bindings,
        only_local_input_and_local_forward_moments_zero_at_inlet=True,
        nonzero_incoming_moments_radial_velocity_and_energy_preserved=True,
        original_selected_energy_equation_identifies_backward_and_forward_sources=True,
        inlet_equality_is_functional_not_interval_overlap=True,
        xi_point02_join_uses_identical_original_source_function=True)


class CompliantPulseEntranceSimilarityC4(CompliantPulseMainExitSimilarityC4):
    @source_precision
    def __init__(self,cells=64):
        super().__init__(cells)
        for stem in ('pulse_main_exit_similarity_C4','pulse_main_exit_similarity_C4_check'):
            name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
            if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Entrance and accepted main source family differ')
            if 'all_passed' in record and not record['all_passed']:raise ValueError('Unaccepted original main source')
            for path,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Accepted original main source changed: '+path)
            self.hashes.update(record['input_hashes'])
            self.hashes[name]=hashlib.sha256(raw).hexdigest()
        self.exporter,self.exporter_proof=compiled_entrance_exporter()
        self.hashes.update(self.exporter_proof['input_hashes'])
        self.entranceproof=entrance_source_proof(self.records)
        self.hashes.update(self.entranceproof['input_hashes'])

    @source_precision
    def entrance(self,Z,xi):
        packet=self.exporter(self,Z,xi)
        c=self.ctx;xi=c.mpf(xi);z0=c.mpf(Z)
        z=main.IntervalTaylor.variable(c,z0,5)
        C=main.IntervalTaylor(c,[1+z0**2,2*z0,1,0,0,0]).reciprocal()
        # The full-pulse K argument is unused for xi<=.02; the original
        # positive integral from zero is evaluated directly, without K-total.
        partial=gp_energy(c,xi,None,self.cells)
        Bh=packet['main_exit_source_rows']['Bh']
        energy=forward_entrance_energy_rows(c,self.mu,xi,self.ap,self.incoming_energy,partial,Bh)
        zero=C*0;zeros=[zero]*5
        coefficients=main.pulse_coefficients(self.delta,self.mu,z,C,self.Xp,
            zeros,zeros,zeros,energy,zeros,zeros)['axial']['full_energy_and_pressure']
        packet['original_partial_entrance_energy_bound']=partial
        packet['forward_full_energy_source_rows']=energy
        packet['forward_full_energy_raw_moment']=dict(
            exact_source_log_parts=packet['five_raw_cumulative_moment_log_sectors']['z_theta']['full_unperturbed_energy']['exact_source_log_parts'],
            full_moment_mixed4_coefficient_enclosures=main.ordinary_grid(
                main.shifted_rows([C*C*v for v in energy],-2*self.mu,4),4))
        packet['forward_full_energy_only_axial_stress']=dict(
            exact_source_log_parts=packet['full_meridional_stress_log_sectors']['axial']['full_energy_and_pressure']['exact_source_log_parts'],
            full_stress_mixed3_coefficient_enclosures=main.ordinary_grid(coefficients['full_derivative_rows'],3))
        packet.update(**{flag:False for flag in FALSE_FLAGS})
        packet.update(whole_original_entrance_similarity_source=True,
            incoming_histories_and_signed_memory_not_zeroed=True,
            selected_backward_energy_exactly_same_as_original_forward_inlet_energy=True)
        return packet

    @source_precision
    def report(self):
        c=self.ctx
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,domain=DOMAIN,
            whole_original_entrance=self.entrance([-1,1],[0,'.02']),
            original_inlet=self.entrance([-1,1],0),
            early_ordinary_y_chart=self.entrance([-1,1],self.mu*c.mpf([0,1])),
            common_entrance_main=self.entrance([-1,1],'.02'),
            source_function_proof=self.entranceproof,
            unchanged_original_exporter_proof=self.exporter_proof,
            current_selected_ap_C5_enclosure=self.ap,current_incoming_C5_enclosures=self.incoming,
            current_incoming_energy_C5_enclosure=self.incoming_energy,
            actual_original_whole_entrance_similarity_companion_constructed=True,
            original_inlet_and_entrance_main_similarity_functional_joins_verified=True,
            source_caps_used_as_defining_field_values=False,
            input_hashes=dict(self.hashes,**{Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}),
            **{flag:False for flag in FALSE_FLAGS})


@source_precision
def run():
    result=CompliantPulseEntranceSimilarityC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole original entrance five-moment/pressure/velocity/stress companion generated; physical/cone pending',flush=True)
    return result


if __name__=='__main__':run()

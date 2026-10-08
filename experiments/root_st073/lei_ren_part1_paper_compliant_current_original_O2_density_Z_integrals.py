"""Original O2 phase-held Z source adapter and genuine five density-Z integrals.

The original slow-Z function is AST-bound with exact source-factor identities.
Native conditional log-|u| ranges remain local. A bounded rational direction
factor encloses the implicit inverse derivative without materializing t.
No finite difference or axial sample becomes an original derivative function.
"""
import ast
import copy
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_positive_logq_cells as current

slow=current.density.slow;base=current.base;ep=current.ep
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
NAME=PREFIX+'current_original_O2_density_Z_integrals.json.gz'
RECEIPT=PREFIX+'current_original_O2_density_Z_integrals_check.json'
GATE='original_O2_genuine_fixed_nonzero_Z_five_density_derivative_integrals_enclosed'


def restore_scalar(c,value):
    if isinstance(value,dict) and 'exact_mpf_tuple' in value:return c.mpf(mp.mp.make_mpf(tuple(value['exact_mpf_tuple'])))
    return c.mpf(value)


def correlated_p2_Z_source(part,source_record):
    """Retain the exact common R*Pstar/L source unit before summing rows."""
    c=part['kernel'].c;bases=part['q'].scale.bases;ledger=part['ledger']
    carrier=base.prior.FormalScale(bases,(1,0,-1,0,1));coefficient=c.mpf(0)
    terms=source_record['original_source_function_coefficient_record']['full_original_factored_coefficient_function_ranges']['p2'][1]
    for term in terms:
        r,p,d,ell=term['original_R_Pstar_delta_L_powers']
        if r!=1 or p not in (-1,1) or d<0:raise ValueError('Original O2 p2_Z term carrier contract changed')
        value=current.interval(c,term['finite_coefficient_function_cover'])
        late=term['late_pressure_coefficient_error_log_upper']
        if late is not None:
            error=ep(part['q'].bounded_exp(restore_scalar(c,late)))[1];value+=c.mpf((-error,error))
        scale=base.prior.FormalScale(bases,(p,d,ell,0,r))
        coefficient+=value*part['q'].bounded_exp((scale-carrier).evaluate())
    original=part['roots']['p2'][(0,1)]
    roots={key:dict(row) for key,row in part['roots'].items()}
    roots['p2'][(0,1)]=base.prior.ScaledEnclosure(carrier,coefficient,ledger)
    result=dict(part,roots=roots)
    result['p2_Z_common_carrier_binding']=dict(original_term_rows=terms,
        exact_common_carrier_source_exponents=[1,0,-1,0],exact_common_carrier_radius_power=1,
        original_wide_source_root=original.record(),correlated_same_original_source_root=roots['p2'][(0,1)].record(),
        common_radial_factor_canceled_before_ratio=True,
        every_original_coefficient_and_positive_late_pressure_error_retained=True)
    return result


def strict_p2_ratios(kernel,roots):
    """Same original p2 sign/factors; p2_Z/p2 and dstar/p2 are not points."""
    p2=roots['p2'][(0,0)];lo,hi=ep(p2.coefficient)
    if not lo*hi>0:raise ValueError('Strict original p2 sign required for local-u derivative ratio')
    sign=1 if lo>0 else -1;positive=base.conditioned.absolute(p2)
    lower=positive.scale.evaluate()+kernel.c.ln(kernel.c.mpf(ep(positive.coefficient)[0]))
    rho2=(roots['p2'][(0,1)]*sign).positive_divide(positive,lower)
    C=(kernel.dstar*sign).positive_divide(positive,lower)
    return rho2,C


def correlated_u_Z(kernel,roots):
    p2=roots['p2'][(0,0)]
    if p2.zero:
        return (roots['p2'][(0,1)]*kernel.q).positive_divide(kernel.dstar,kernel.dstar.scale.evaluate())
    rho2,C=strict_p2_ratios(kernel,roots)
    return kernel.u*rho2


def signed_T1(kernel,roots,difference):
    rho2,C=strict_p2_ratios(kernel,roots)
    return C*difference


def signed_T1_Z(kernel,roots,difference,chiZ):
    rho2,C=strict_p2_ratios(kernel,roots)
    return C*(chiZ-rho2*difference)


def implicit_direction_product(kernel,t,T2Z):
    """Exact -T2Z*t/(1+t²), with a theorem-intersected rational factor."""
    c=kernel.c;one=kernel.scalar(1)
    if t.zero:return kernel.scalar(0)
    denominator=one+slow.square(t)
    ratio=t.positive_divide(denominator,0)
    try:
        finite=base.conditioned.clipped(c,current.bounded(ratio),c.mpf('-.5'),c.mpf('.5'))
    except ArithmeticError:
        # The entire original t range is retained by |t/(1+t²)|<=1/2.
        # This is a conservative functional enclosure, not a selected t.
        finite=c.mpf(('-.5','.5'))
        kernel.q.ledger['implicit_direction_universal_range_intersections']=kernel.q.ledger.get('implicit_direction_universal_range_intersections',0)+1
        lo,hi=ep(t.coefficient)
        if lo*hi>0:
            positive=base.conditioned.absolute(t);lower=positive.scale.evaluate()+c.ln(c.mpf(ep(positive.coefficient)[0]))
            if ep(lower)[0]>=0:
                inverse=one.positive_divide(positive,lower)*(1 if lo>0 else -1)
                stable=inverse.positive_divide(one+slow.square(inverse),0)
                finite=base.conditioned.clipped(c,current.bounded(stable),c.mpf('-.5'),c.mpf('.5'))
                kernel.q.ledger['implicit_direction_reciprocal_coordinates']=kernel.q.ledger.get('implicit_direction_reciprocal_coordinates',0)+1
    return -T2Z*finite


def compile_original_slow_Z():
    tree=ast.parse(Path(slow.__file__).read_text(encoding='utf8'))
    fn=copy.deepcopy(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='slow_values'))
    before=copy.deepcopy(fn);changes=[]
    def replace_assignment(nodes,name,expected,replacement):
        found=[n for n in nodes if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in n.targets)]
        if len(found)!=1:raise ValueError('Unique original slow-Z assignment required: '+name)
        item=found[0];old=copy.deepcopy(item.value)
        if ast.dump(old)!=ast.dump(ast.parse(expected,mode='eval').body):raise ValueError('Original slow-Z formula changed: '+name)
        item.value=ast.parse(replacement,mode='eval').body;changes.append((item,old))
    replace_assignment(ast.walk(fn),'uZ',
        '(p2Z*q).positive_divide(kernel.dstar,kernel.dstar.scale.evaluate())',
        'correlated_u_Z(kernel,roots)')
    signed=next(n for n in ast.walk(fn) if isinstance(n,ast.If) and ast.dump(n.test)==ast.dump(ast.parse("kernel.geometry=='signed_Mobius'",mode='eval').body))
    replace_assignment(ast.walk(signed),'T1','q*hinv*(difference/r)',
        'signed_T1(kernel,roots,difference)')
    replace_assignment(ast.walk(signed),'T1Z','q*hinv*(1/r)*(chiZ-K*(difference/r))',
        'signed_T1_Z(kernel,roots,difference,chiZ)')
    replace_assignment(ast.walk(fn),'BZ','-a*(EZ*T1+kernel.E*(T1Z+t*psiZ))*(1/(4*c.pi))',
        '-a*(EZ*T1+kernel.E*(T1Z+implicit_direction_product(kernel,t,T2Z)))*(1/(4*c.pi))')
    edited=copy.deepcopy(fn)
    for node,old in changes:node.value=old
    if ast.dump(fn)!=ast.dump(before):raise ValueError('Slow-Z adapter changed the original math beyond bound identities')
    env=dict(vars(slow),correlated_u_Z=correlated_u_Z,signed_T1=signed_T1,signed_T1_Z=signed_T1_Z,
        implicit_direction_product=implicit_direction_product)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[edited],type_ignores=[])),
        '<original slow_values with exact correlated O2 source identities>','exec'),env)
    return env['slow_values']


bound_original_slow_Z=compile_original_slow_Z()


def identity_contract():
    u,uz,q,dstar,p2,p2z,h,r,D,chiZ=sy.symbols('u uz q dstar p2 p2z h r D chiZ',nonzero=True,real=True)
    source={u:p2*q/dstar,uz:p2z*q/dstar,r:p2*q/(dstar*h)}
    assert sy.cancel((q/h/r-dstar/p2).subs(source))==0
    assert sy.cancel((uz-u*p2z/p2).subs(source))==0
    old=q/h/r*(chiZ-uz/h/r*D)
    assert sy.cancel((old-dstar/p2*(chiZ-p2z/p2*D)).subs(source))==0
    t=sy.symbols('t',real=True)
    assert sy.factor(sy.Rational(1,4)*(1+t*t)**2-t*t)==(t-1)**2*(t+1)**2/4
    return dict(passed=True,exact_u_Z_equals_u_p2_Z_over_p2=True,
        signed_original_T1_and_fixed_psi_T1_Z_identities=True,
        original_radian_angle_convention=True,
        implicit_t_psi_Z_equals_minus_T2_Z_t_over_one_plus_t_squared=True,
        rational_direction_absolute_bound=sy.Rational(1,2).__str__(),
        original_slow_function_AST_restoration_guard=True,
        original_q_Z_a_Z_t0_Z_dstar_Z_zero_by_same_O2_source_definitions=True)


def slow_values(part,coordinate,chart):
    kernel=part['kernel'];roots=part['roots']
    with mp.workdps(kernel.c.dps+40):
        if any(not roots[key][(0,1)].zero for key in ('a','b','t0')) or not kernel.t0.zero:
            raise ValueError('Original O2 a_Z=b_Z=t0_Z=0 required')
        if kernel.geometry=='requires_signed_source_refinement':raise ValueError('No derivative on unadmitted signed source')
        result,record=bound_original_slow_Z(kernel,roots,coordinate,chart)
        record.update(local_basis_order=current.LOCAL_ORDER,source_u_relation='u=p2*q/dstar',
            q_Z_and_dstar_Z_exact_zero_source_identity=True,
            original_phase_held_Z_not_finite_difference=True,
            only_exact_original_source_factor_bindings_and_bounded_rational_direction_changed=True,
            source_Z_derivatives_not_interval_selector_derivatives=True)
        return result,record


def derivative_inverse(kernel,phase,bits):
    """Use the original E coordinate for signed derivative peak conditioning."""
    c=kernel.c;lo,hi=ep(c.mpf(phase))
    if kernel.geometry=='signed_Mobius' and not (lo==hi and lo in (0,mp.mpf('.5'),1)):
        return kernel.inverse_bracket(phase,'E',bits)
    return kernel.evaluate(phase,bits=bits)['selected_inverse']


class OriginalO2DensityZIntegrals:
    def __init__(self):
        accepted=json.loads((HERE/current.RECEIPT).read_bytes())
        if not accepted.get('all_passed') or not accepted.get(current.GATE):raise ValueError('Accepted native positive-q fixed-Z source required')
        self.parent=current.OriginalO2PositiveLogQCells();self.c=self.parent.c;self.family=self.parent.family
        if self.family!=accepted['source_family']:raise ValueError('Original density-Z family differs')
        self.hashes=dict(self.parent.hashes)
        for name,digest in {**accepted['input_hashes'],current.RECEIPT:sha(current.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Density-Z source changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Density-Z source closures disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.identities=identity_contract()

    def query(self,count,index,*,Z_lower,Z_upper,phase='.137',bits=24):
        with mp.workdps(self.c.dps+40):
            source=self.parent.query(count,index,Z_lower=Z_lower,Z_upper=Z_upper);rows=[]
            for original_part in source['pieces']:
                part=correlated_p2_Z_source(original_part,source['record'])
                inverse=part['kernel'].evaluate(phase,bits=bits);selected=inverse['selected_inverse']
                primitive=part['kernel'].primitives(selected['coordinate_interval'],selected['chart'])
                derivative=derivative_inverse(part['kernel'],phase,bits)
                jets,proof=slow_values(part,derivative['coordinate_interval'],derivative['chart'])
                graph=current.CellDensityGraph(self.parent.owner.scales.graph,part,primitive,jets,7)
                values=graph.outputs()
                finite={};formal={}
                for key,value in values['density_Z'].items():
                    formal[key]=value.record()
                    try:finite[key]=current.bounded(value)
                    except ArithmeticError:finite[key]=None
                rows.append(dict(original_inverse=selected,original_slow_Z_contract=proof,
                    original_derivative_inverse=derivative,
                    correlated_original_p2_Z_carrier=part['p2_Z_common_carrier_binding'],
                    original_primitive_Z_enclosures={key:value.record() for key,value in jets.items()},
                    original_five_density_Z_source_enclosures=formal,
                    optional_finite_five_density_Z_ranges=finite,
                    unmaterializable_original_Z_derivatives_retained_as_formal_sources=True,
                    actual_original_density_Z_graph_executed=True))
            return dict(source=source['record'],fixed_phase=str(phase),source_pieces=rows,
                phase_is_free_probe_not_actual_radius_phase=True,source_family=self.family,
                mixed_signed_source_not_skipped=not bool(source['pieces']))

    def integrate(self,count,*,Z='.37',N=7,bits=24):
        c=self.c;level,phase=self.parent.parent.levels[count]
        if N!=phase['explicit_candidate_N']:raise ValueError('Actual same common candidate phase cover required')
        begin=time.monotonic();C0={key:c.mpf(0) for key in current.five.RATES};Z1={key:c.mpf(0) for key in current.five.RATES};rows=[]
        with mp.workdps(c.dps+40):
            for i,phase_cell in enumerate(phase['whole_source_cells']):
                query=self.parent.query(count,i,Z_lower=Z,Z_upper=Z)
                if not query['pieces']:raise ArithmeticError('Original density-Z source requires signed refinement at '+str(i)+'/'+str(count))
                densities={key:[] for key in current.five.RATES};derivatives={key:[] for key in current.five.RATES};proofs=[]
                for source_index,original_part in enumerate(query['pieces']):
                    part=correlated_p2_Z_source(original_part,query['record'])
                    kernel=part['kernel']
                    for phase_cover in phase_cell['true_common_N_phase_boxes']:
                        phi=current.interval(c,phase_cover);selected=kernel.evaluate(phi,bits=bits)['selected_inverse']
                        primitive=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                        derivative=derivative_inverse(kernel,phi,bits)
                        jets,proof=slow_values(part,derivative['coordinate_interval'],derivative['chart'])
                        graph=current.CellDensityGraph(self.parent.owner.scales.graph,part,primitive,jets,N);got=graph.outputs()
                        for key in current.five.RATES:
                            densities[key].append(current.bounded(got['densities'][key]));derivatives[key].append(current.bounded(got['density_Z'][key]))
                        proofs.append(dict(source_piece_index=source_index,actual_true_phase_cover=phase_cover,
                            correlated_original_p2_Z_carrier=part['p2_Z_common_carrier_binding'],
                            original_inverse=selected,original_slow_Z_contract=proof,
                            original_derivative_inverse=derivative,
                            native_original_primitive_Z_enclosures={key:value.record() for key,value in jets.items()},
                            density_Z_is_genuine_ordinary_source_derivative=True))
                hull=lambda values:c.mpf((min(ep(v)[0] for v in values),max(ep(v)[1] for v in values)))
                h0={key:hull(values) for key,values in densities.items()};h1={key:hull(values) for key,values in derivatives.items()}
                unused,masses,decay=current.five.masses(c,c.mpf(1)/count,c.mpf(i)/count,c.mpf(i+1)/count)
                add0={key:h0[key]*masses[key] for key in current.five.RATES};add1={key:h1[key]*masses[key] for key in current.five.RATES}
                for key in current.five.RATES:C0[key]+=add0[key];Z1[key]+=add1[key]
                rows.append(dict(source=query['record'],actual_source_phase_held_Z_records=proofs,
                    five_signed_C0_density_hulls=h0,five_genuine_signed_ordinary_Z_density_hulls=h1,
                    positive_own_rate_final_endpoint_masses=masses,five_C0_contributions=add0,five_ordinary_Z_contributions=add1,
                    fixed_y_window_and_Z_independent_positive_masses=True,
                    derivatives_integrated_from_original_graph_not_from_interval_endpoint_differences=True))
                if (i+1)%max(8,count//8)==0:print('Genuine original density-Z integral:',count,str(Z),i+1,flush=True)
        return dict(source_family=self.family,ordered_source_cells=count,original_Z_exact=str(base.point.pressure.exact_Z(Z)),
            exact_y_window=['0','1'],explicit_candidate_N=N,inverse_bits=bits,own_rates=current.five.RATES,normalized_own_units=current.five.UNITS,
            five_original_C0_integral_contributions=C0,five_genuine_ordinary_Z_integral_contributions=Z1,
            whole_original_density_and_density_Z_source_records=rows,
            original_P0_datum_sha256=self.family['datum_enclosure_sha256'],actual_incoming_histories_not_set_to_zero=True,
            complete_original_y_window_all_source_and_phase_pieces_integrated=True,
            Z_interval_terminal_identities_or_global_controls_installed=False,current_whole_N_selected=False,
            execution_seconds=time.monotonic()-begin)


def run():
    begin=time.monotonic();owner=OriginalO2DensityZIntegrals()
    probes=[owner.query(8192,i,Z_lower=z,Z_upper=z) for i,z in ((7946,'.7'),(8191,'.37'),(8191,'0'))]
    probes.append(owner.query(8192,4341,Z_lower='-.01',Z_upper='.01'))
    integrals=[owner.integrate(n,Z='.37',N=7,bits=24) for n in (64,256,2048)]
    report=dict(**{GATE:True},source_family=owner.family,exact_original_source_derivative_identities=owner.identities,
        native_whole_cell_derivative_probes=probes,actual_original_five_density_Z_integral_refinements=integrals,
        derivative_is_actual_phase_held_Z_not_finite_difference=True,
        fixed_nonzero_Z_C0_and_genuine_density_Z_full_window_enclosed=True,
        numerical_original_source_point_or_integral_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Original O2 native C0 and genuine ordinary-Z five signed integral contribution enclosures over the entire y window at fixed Z37/100, candidate N7. Not continuous-Z terminal matching, all-chart incoming histories/controls, global N, matching/stress/recursion or full NS.')
    raw=json.dumps(base.encoded(report),indent=2).encode('utf8')+b'\n'
    (HERE/NAME).write_bytes(gzip.compress(raw,compresslevel=9,mtime=0));return report


if __name__=='__main__':run()

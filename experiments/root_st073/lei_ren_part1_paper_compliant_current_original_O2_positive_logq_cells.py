"""Native positive O2 log-q cells and nonzero-Z own-rate phase integrals.

The original cutoff complement is evaluated logarithmically. A local fourth
log-|u| slot gives stable conditional coordinates at extraordinary scales.
It is not the standard packet basis; normalized finite density ranges are
exported before any cross-cell summation. No source endpoint becomes a field.
"""
import ast
import copy
import gzip
import hashlib
import json
from pathlib import Path
import time
from types import FunctionType
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_ordered_source_cells as ordered

five=ordered.five;base=ordered.base;density=five.density;ep=ordered.ep
HERE,PREFIX,sha=ordered.HERE,ordered.PREFIX,ordered.sha
NAME=PREFIX+'current_original_O2_positive_logq_cells.json.gz'
RECEIPT=PREFIX+'current_original_O2_positive_logq_cells_check.json'
GATE='original_O2_native_positive_logq_source_cells_and_nonzero_Z_phase_transport_connected'
LOCAL_ORDER=['logPstar','logdelta','logL','logabsu','logR']


def interval(c,value):return ordered.interval(c,value)


def bounded(value):
    """Enclose a bounded value even if only its log lower is microscopic."""
    if value.zero:return value.ctx.mpf(0)
    return value.coefficient*value.bounded_exp(value.scale.evaluate())


class PositiveLogQPhase(base.conditioned.ConditionedPhase):
    """Unchanged source method code, wider bounded-log arithmetic callback."""

    def primitives(self,coordinate,chart='psi'):
        result=base.conditioned.ConditionedPhase.primitives(self,coordinate,chart)
        if self.geometry=='signed_Mobius':
            if not self.t0.zero:raise ValueError('Original O2 t0=0 source required by canceled B identity')
            p2=self.roots['p2'][(0,0)]
            if not ep(p2.coefficient)[0]*ep(p2.coefficient)[1]>0:raise ValueError('Strict original p2 sign required')
            sign=1 if ep(p2.coefficient)[0]>0 else -1
            if sign!=self.sign:raise ValueError('Original r and p2 source signs differ')
            positive=base.conditioned.absolute(p2)
            ratio=self.dstar.positive_divide(positive,positive.scale.evaluate()+self.c.ln(self.c.mpf(ep(positive.coefficient)[0])))
            psi,E=self.angles(self.c.mpf(coordinate),chart)
            # Exact source identity q/(h*r)=dstar/p2. The original
            # small-r expression is retained; no division near r=0 occurs.
            result['B_over_Pstar']=-self.E*self.a*ratio*((E-psi)*sign/2)
        return result


SPECIALIZED_METHODS=[]
for _name,_method in vars(base.conditioned.ConditionedPhase).items():
    if isinstance(_method,FunctionType) and 'bounded_value' in _method.__code__.co_names:
        _globals=dict(_method.__globals__,bounded_value=bounded)
        _copy=FunctionType(_method.__code__,_globals,_method.__name__,_method.__defaults__,_method.__closure__)
        _copy.__kwdefaults__=_method.__kwdefaults__
        setattr(PositiveLogQPhase,_name,_copy);SPECIALIZED_METHODS.append(_name)


def compile_original_u_binding():
    tree=ast.parse(Path(base.conditioned.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='ConditionedPhase')
    fn=copy.deepcopy(next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__'))
    before=copy.deepcopy(fn)
    assignment=next(n for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='self.u' for t in n.targets))
    expected="(self.roots['p2'][ZERO]*self.q).positive_divide(self.dstar,dstar_log)"
    if ast.dump(assignment.value)!=ast.dump(ast.parse(expected,mode='eval').body):raise ValueError('Original u definition changed')
    assignment.value=ast.parse("query['original_u_source']",mode='eval').body
    restored=copy.deepcopy(fn)
    target=next(n for n in ast.walk(restored) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='self.u' for t in n.targets))
    target.value=ast.parse(expected,mode='eval').body
    if ast.dump(restored)!=ast.dump(before):raise ValueError('Original phase math changed beyond equivalent source-u binding')
    env=dict(vars(base.conditioned),bounded_value=bounded)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original ConditionedPhase; equivalent original u source binding>', 'exec'),env)
    return env['__init__']


PositiveLogQPhase.__init__=compile_original_u_binding()
SPECIALIZED_METHODS.remove('__init__')


def small_exp(c,log):
    cutoff=-2*c.dps*c.ln(10);lo,hi=ep(log)
    if hi<ep(cutoff)[0]:return c.mpf((0,ep(c.exp(cutoff))[1]))
    if lo<ep(cutoff)[0]:return c.mpf((0,ep(c.exp(c.mpf(hi)))[1]))
    return c.exp(log)


def log_one_plus_exp(c,nonpositive):
    if ep(nonpositive)[1]>0:raise ValueError('Nonpositive log input required')
    return c.ln(1+small_exp(c,nonpositive))


def log_sum(c,left,right):
    # The upper endpoint only selects arithmetic coordinates, never a source
    # value. Both original log ranges are retained inside the exponential.
    reference=c.mpf(max(ep(left)[1],ep(right)[1]))
    return reference+c.ln(small_exp(c,left-reference)+small_exp(c,right-reference))


def endpoint_logq(c,y,eta_log):
    y=sy.Rational(y);iy=ordered.rational(c,y)
    if y==0:
        return (log_sum(c,c.ln(c.mpf(3)/5),eta_log)-c.ln(c.mpf(4)/5))/2
    if y==1:return (eta_log-c.ln(2))/2
    odds=1/(iy*iy)-1/(1-iy)**2
    if y>=sy.Rational(1,2):
        # log(1-sigma)=odds-log1p(exp(odds)); keep the small
        # complement instead of subtracting a rounded sigma from one.
        odds=base.conditioned.clipped(c,odds,ep(odds)[0],0)
        log_complement=odds-log_one_plus_exp(c,odds)
    else:
        minus=base.conditioned.clipped(c,-odds,ep(-odds)[0],0)
        log_complement=-log_one_plus_exp(c,minus)
    # a=2-1.2*(1-sigma). A finite source cover is sufficient in
    # this denominator; it does not replace the positive complement.
    complement=small_exp(c,log_complement)
    a=base.conditioned.clipped(c,2-c.mpf(6)/5*complement,c.mpf(4)/5,2)
    return (log_sum(c,c.ln(c.mpf(3)/5)+log_complement,eta_log)-c.ln(a))/2


def q_identity():
    sigma,eta=sy.symbols('sigma eta',real=True);a=sy.Rational(4,5)+sy.Rational(6,5)*sigma
    assert sy.cancel((2+2*eta-a)/(2*a)-(sy.Rational(3,5)*(1-sigma)+eta)/a)==0
    s=sy.Symbol('s',real=True);expr=(sy.Rational(3,5)*(1-s)+eta)/(sy.Rational(4,5)+sy.Rational(6,5)*s)
    assert sy.cancel(sy.diff(expr,s)+sy.Rational(6,5)*(1+eta)/(sy.Rational(4,5)+sy.Rational(6,5)*s)**2)==0
    return dict(passed=True,exact_original_q_squared_identity=True,
        original_cutoff_equals_one_by_a_at_most_two_and_eta_positive=True,
        derivative='d(q²)/dsigma=-(6/5)*(1+eta)/a²<0',
        original_endpoints='q²(0)=(3/5+eta)/(4/5); q²(1)=eta/2',
        monotone_q_log_endpoint_hull_covers_entire_source_y_cell=True)


class CellDensityGraph(density.BoundDensityGraph):
    def modulation_exp(self,index):
        """Same original stable graph operation with bounded-log arithmetic."""
        if index in self.expcache:return self.expcache[index]
        x=self.evaluate(index);c=self.c
        finite=base.conditioned.clipped(c,bounded(x),-1,1)
        M=64;mean=c.mpf(1);power=c.mpf(1)
        for k in range(1,M+1):power*=finite;mean+=power/c.factorial(k+1)
        tail=ep(c.exp(1)/c.factorial(M+2))[1]
        mean=base.conditioned.clipped(c,mean+c.mpf((-tail,tail)),c.exp(-1),c.exp(1))
        changed=x*mean;result=(self.scalar(1)+changed,changed)
        self.ledger['bounded_modulation_exp_taylor_tail_operations']=self.ledger.get('bounded_modulation_exp_taylor_tail_operations',0)+1
        self.expcache[index]=result;return result


class OriginalO2PositiveLogQCells:
    def __init__(self):
        accepted=json.loads((HERE/ordered.RECEIPT).read_bytes())
        if not accepted.get('all_passed') or not accepted.get(ordered.GATE):raise ValueError('Accepted original ordered source-cell ranges required')
        self.parent=ordered.OriginalO2OrderedSourceCells();self.owner=self.parent.owner;self.c=self.owner.ctx;self.family=self.parent.family
        if accepted['source_family']!=self.family:raise ValueError('Positive log-q source family differs')
        self.hashes=dict(self.parent.hashes)
        for name,digest in {**accepted['input_hashes'],ordered.RECEIPT:sha(ordered.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Positive log-q source changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Positive log-q hash closures disagree')
            self.hashes[name]=digest
        saved=accepted['compressed_producer_report'];raw=gzip.decompress((HERE/saved['filename']).read_bytes())
        if hashlib.sha256(raw).hexdigest()!=saved['lossless_original_json_sha256']:raise ValueError('Ordered source archive digest differs')
        self.saved=json.loads(raw);self.identity=q_identity();self.logs={}
        c=self.c
        for level,pressure in zip(self.saved['actual_original_ordered_source_levels'],
                self.parent.parent.saved['actual_original_pressure_integral_refinements'],strict=True):
            count=level['ordered_source_cells']
            if count!=pressure['ordered_source_cells'] or level['source_family']!=self.family:raise ValueError('Ordered/phase source levels differ')
            live=dict(level);live['whole_source_cells']=[dict(row,whole_original_radial_profile_covers={
                key:interval(c,value) for key,value in row['whole_original_radial_profile_covers'].items()})
                for row in level['whole_source_cells']]
            self.parent.levels[count]=(live,pressure)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def log_nodes(self,count):
        if count not in self.logs:
            self.logs[count]=[endpoint_logq(self.c,sy.Rational(i,count),self.owner.scales.logs['eta']) for i in range(count+1)]
        return self.logs[count]

    def query(self,count,index,*,Z_lower,Z_upper):
        # Absolute endpoint values and diagnostic records use mpmath's global
        # context too. Keep them at the same guarded native precision as the
        # full integral, including a standalone public query.
        with mp.workdps(self.c.dps+40):
            return self._query(count,index,Z_lower=Z_lower,Z_upper=Z_upper)

    def _query(self,count,index,*,Z_lower,Z_upper):
        original=self.parent.coefficient_cell(count,index,Z_lower=Z_lower,Z_upper=Z_upper);c=self.c
        node=self.log_nodes(count);logq=ordered.hull(c,node[index+1],node[index])
        old_bases=original['q'].scale.bases;old_ledger=original['ledger']
        if ep(old_bases[3])!=(0,0):raise ValueError('Only original zero fourth-slot roots can be locally rebound')
        p2=original['roots']['p2'][(0,0)];pl,ph=ep(p2.coefficient)
        sign=1 if pl>0 else -1 if ph<0 else 0
        abslo=min(abs(pl),abs(ph)) if sign else 0;abshi=max(abs(pl),abs(ph))
        logu=None if not sign else p2.scale.evaluate()+c.ln(c.mpf((abslo,abshi)))+logq-self.owner.scales.logs['d_star']
        def bind(logu_cover):
            bases=old_bases[:3]+(c.mpf(0) if logu_cover is None else logu_cover,)+old_bases[4:];bound_roots={}
            for key,row in original['roots'].items():
                bound_roots[key]={}
                for order,value in row.items():
                    if value.scale.powers[3]!=0:raise ValueError('Original root uses the reserved local source-u slot')
                    bound_roots[key][order]=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,value.scale.powers,value.scale.offset),value.coefficient,old_ledger)
            bound_q=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,offset=logq),1,old_ledger)
            if sign:
                bound_u=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,(0,0,0,1,0)),sign,old_ledger)
            else:
                dstar=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,offset=self.owner.scales.logs['d_star']),1,old_ledger)
                bound_u=(bound_roots['p2'][(0,0)]*bound_q).positive_divide(dstar,self.owner.scales.logs['d_star'])
            bound_kernel=PositiveLogQPhase(dict(q=bound_q,roots=bound_roots,original_u_source=bound_u),self.owner.scales.logs['d_star'])
            return dict(roots=bound_roots,q=bound_q,kernel=bound_kernel,ledger=old_ledger,original_u_source=bound_u)
        full=bind(logu);roots,q,kernel=full['roots'],full['q'],full['kernel']
        pieces=[full] if kernel.geometry!='requires_signed_source_refinement' else [];split=None
        if not pieces and sign:
            margin=c.ln(2)/100
            small_limit=ep(c.ln(c.mpf('.25'))-margin)[0]
            large_limit=ep(c.ln(c.mpf('.125'))+margin)[1]
            ulo,uhi=ep(logu);covers=[]
            if ulo<=small_limit:covers.append(c.mpf((ulo,min(uhi,small_limit))))
            if large_limit<=uhi:covers.append(c.mpf((max(ulo,large_limit),uhi)))
            if large_limit>small_limit:raise ArithmeticError('Original small/large coordinate overlap lost')
            for cover in covers:
                part=bind(cover)
                if part['kernel'].geometry=='requires_signed_source_refinement':raise ArithmeticError('Conditioned original source-u piece needs refinement')
                pieces.append(part)
            split=dict(small_source_logabsu_upper=small_limit,large_source_logabsu_lower=large_limit,
                overlapping_original_small_and_large_u_coordinates=True,piece_logabsu_intervals=covers,
                entire_native_logabsu_source_range_covered=True,
                conditional_original_u_equals_p2_q_over_dstar_on_each_piece=True)
        record=dict(source_family=self.family,ordered_source_cell_level=count,ordered_source_cell_index=index,
            exact_y_cell=original['record']['exact_y_cell'],exact_Z_range=original['record']['exact_Z_range'],
            original_source_function_coefficient_record=original['record'],whole_native_positive_logq_interval=logq,
            original_selected_eta_log=self.owner.scales.logs['eta'],q_identity=self.identity,
            local_basis_contract=dict(order=LOCAL_ORDER,original_root_fourth_powers_all_zero=True,
                original_roots_unchanged_by_rebinding_zero_powers=True,q_fourth_power_exactly_zero=True,
                signed_u_fourth_power_exactly_one=bool(sign),strict_original_p2_sign=sign,
                original_u_relation='u=p2*q/dstar; magnitude retained as a source log range, not a selected coefficient',
                original_phase_math_unchanged_except_equivalent_u_source_binding=True,
                one_shared_basis_and_ledger_for_this_cell=True,standard_packet_basis_compatible=False,
                normalized_finite_density_ranges_exported_before_cross_cell_sum=True),
            native_positive_q=q.record(),conditioned_geometry=kernel.geometry_record(),
            conditioned_source_piece_geometry=[piece['kernel'].geometry_record() for piece in pieces],
            overlapping_source_coordinate_split=split,original_method_code_objects_unchanged=SPECIALIZED_METHODS,
            original_positive_q_not_replaced_by_flat_or_field_point=True,
            inverse_or_global_controls_installed=False)
        return dict(record=record,roots=roots,q=q,kernel=kernel,ledger=old_ledger,pieces=pieces)

    def integrate(self,count,*,Z,N=7,bits=24):
        c=self.c;level,phase=self.parent.levels[count]
        if N!=phase['explicit_candidate_N']:raise ValueError('Same explicit common N as accepted true phase cover required')
        begin=time.monotonic();changes={key:c.mpf(0) for key in five.RATES};records=[]
        with mp.workdps(c.dps+40):
            for i,saved in enumerate(phase['whole_source_cells']):
                query=self.query(count,i,Z_lower=Z,Z_upper=Z)
                if not query['pieces']:
                    raise ArithmeticError('Native source needs refinement at original cell '+str(i)+'/'+str(count)+' Z='+str(Z))
                pieces={key:[] for key in five.RATES};inverse_rows=[]
                for piece_index,part in enumerate(query['pieces']):
                    kernel=part['kernel']
                    for source_phase in saved['true_common_N_phase_boxes']:
                        phi=interval(c,source_phase);inverse=kernel.evaluate(phi,bits=bits)
                        if inverse['status']!='enclosed':raise ArithmeticError('Actual source phase inverse not enclosed')
                        selected=inverse['selected_inverse'];primitives=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                        graph=CellDensityGraph(self.owner.scales.graph,part,primitives,None,N)
                        for key,value in graph.values()['densities'].items():pieces[key].append(bounded(value))
                        inverse_rows.append(dict(source_coordinate_piece_index=piece_index,true_source_phase=source_phase,
                            selected_original_inverse=selected,native_primitives={key:value.record() for key,value in primitives.items()},
                            local_logabsu_basis_order=LOCAL_ORDER))
                hulls={key:c.mpf((min(ep(v)[0] for v in values),max(ep(v)[1] for v in values))) for key,values in pieces.items()}
                left,right=c.mpf(i)/count,c.mpf(i+1)/count;local,masses,decay=five.masses(c,c.mpf(1)/count,left,right)
                contributions={key:hulls[key]*masses[key] for key in five.RATES}
                for key in five.RATES:changes[key]+=contributions[key]
                records.append(dict(source=query['record'],actual_native_inverse_phase_pieces=inverse_rows,
                    five_signed_normalized_density_hulls=hulls,positive_own_rate_final_endpoint_masses=masses,
                    five_signed_own_rate_contributions=contributions,
                    native_B_original_V_and_all_cross_terms_retained=True,
                    finite_enclosures_summed_after_leaving_local_u_factor_basis=True))
                if (i+1)%max(16,count//8)==0:print('Native nonzero-Z O2 integral:',count,str(Z),i+1,flush=True)
        return dict(source_family=self.family,exact_y_window=['0','1'],original_Z_exact=str(base.point.pressure.exact_Z(Z)),
            explicit_candidate_N=N,ordered_source_cells=count,inverse_bits=bits,own_rates=five.RATES,normalized_own_units=five.UNITS,
            five_signed_original_nonzero_Z_own_rate_integral_contributions=changes,whole_source_phase_density_records=records,
            original_P0_datum_sha256=self.family['datum_enclosure_sha256'],incoming_five_histories_and_P0_not_reset=True,
            contribution_only_no_default_actual_incoming_defects=True,
            complete_original_source_window_and_true_phase_unions_integrated=True,
            current_whole_N_selected=False,actual_five_controls_installed=False,Z_functional_terminal_matching_installed=False,
            execution_seconds=time.monotonic()-begin)


def run():
    begin=time.monotonic();owner=OriginalO2PositiveLogQCells();count=8192;queries=[]
    for i,zl,zh in ((7946,'.7','.7'),(8191,'.37','.37'),(8191,'0','0'),(4341,'-.01','.01')):
        queries.append(owner.query(count,i,Z_lower=zl,Z_upper=zh)['record'])
    # An initial full-window native nonmidplane integral, then refinements.
    integrals=[]
    for level in (64,256,2048):integrals.append(owner.integrate(level,Z='.37',N=7,bits=24))
    report=dict(**{GATE:True},source_family=owner.family,native_whole_cell_positive_logq_queries=queries,
        actual_original_nonzero_Z_own_integral_refinements=integrals,
        native_positive_logq_source_formula_and_local_cancellation_basis_installed=True,
        one_fixed_nonzero_Z_complete_original_window_integral_enclosed=True,
        standard_packet_basis_compatible=False,actual_changed_five_moment_integral_evaluated=False,
        numerical_original_source_point_or_integral_oracle_installed=False,actual_five_controls_installed=False,
        current_whole_N_selected=False,**dict.fromkeys(base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Native positive O2 log-q source cells and full-window fixed nonzero-Z candidate-N own-rate contributions using original true phase/inverse and both primitives. Local logabsu basis is not a standard packet basis; finite normalized densities are exported before summation. No all-Z functional matching, all-chart oracle, installed controls, admitted global N, stress or recursion.')
    data=json.dumps(base.encoded(report),indent=2).encode('utf8')+b'\n'
    (HERE/NAME).write_bytes(gzip.compress(data,compresslevel=9,mtime=0))
    return report


if __name__=='__main__':run()

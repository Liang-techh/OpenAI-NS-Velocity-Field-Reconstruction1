"""Original ordinary q_y/q_Z and mixed y2,Z1 derivatives on native boxes.

Flat cutoff branches are lazy and exactly zero. Active derivatives use the
same correlated a/Delta roots, positive gamma and original sigma derivatives.
Unresolved branch/normalization covers request subdivision, never fake jets.
"""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_correlated_shear_q as current

prior=current.prior;native=current.native;packets=current.packets
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep
ORDERS=current.ORDERS;ZERO=(0,0)
NAME=PREFIX+'current_native_q_slow_jets.json'
RECEIPT=PREFIX+'current_native_q_slow_jets_check.json'
GATE='current_original_native_q_ordinary_y2_Z1_slow_jet_functions_executed'


def bounded(value):
    if value.zero:return value.ctx.mpf(0)
    if ep(value.scale.evaluate())[1]<-1000:return value.coefficient*value.bounded_exp(value.scale.evaluate())
    return value.finite_interval()


def sqrt_jet(gamma,log_root_lower):
    """Differentiate R^2=gamma; divide by the same positive original 2R."""
    root=current.nonnegative_sqrt(gamma[ZERO]).positive_intersection(log_root_lower)
    result={ZERO:root};zero=root.scalar(0)
    for j,k in sorted(ORDERS,key=lambda order:(sum(order),order)):
        if (j,k)==ZERO:continue
        correction=[]
        for i in range(j+1):
            for ell in range(k+1):
                if (i,ell) in (ZERO,(j,k)):continue
                a,b=result[(i,ell)],result[(j-i,k-ell)]
                term=current.square(a) if a is b else a*b
                correction.append(term*(math.comb(j,i)*math.comb(k,ell)))
        numerator=gamma[(j,k)]-sum(correction,zero)
        result[(j,k)]=numerator.positive_divide(root*2,log_root_lower+root.ctx.ln(2))
    return result


def compose_sigma(s,ordinary):
    """Ordinary mixed chain rule; sigma_jets Taylor rows are multiplied by n!."""
    scalar=s[ZERO].scalar;d=[scalar(v) for v in ordinary]
    return {ZERO:d[0],(1,0):d[1]*s[(1,0)],
        (2,0):d[2]*current.square(s[(1,0)])+d[1]*s[(2,0)],
        (0,1):d[1]*s[(0,1)],
        (1,1):d[2]*s[(1,0)]*s[(0,1)]+d[1]*s[(1,1)],
        (2,1):d[3]*current.square(s[(1,0)])*s[(0,1)]
            +d[2]*(s[(2,0)]*s[(0,1)]+s[(1,0)]*s[(1,1)]*2)+d[1]*s[(2,1)]}


def mixed_chain_identity():
    y,Z=sy.symbols('y Z');s=sy.Function('s')(y,Z);F=sy.Function('sigma')
    first=sy.Subs(sy.diff(F(sy.Symbol('x')),sy.Symbol('x')),sy.Symbol('x'),s)
    # SymPy supplies an independent multivariate chain expansion.
    expected=sy.diff(F(s),y,2,Z)
    actual=sy.diff(F(s),s,3)*sy.diff(s,y)**2*sy.diff(s,Z)+sy.diff(F(s),s,2)*(sy.diff(s,y,2)*sy.diff(s,Z)+2*sy.diff(s,y)*sy.diff(s,y,Z))+sy.diff(F(s),s)*sy.diff(s,y,2,Z)
    if sy.simplify(expected-actual)!=0:raise ArithmeticError('Original y2 Z1 cutoff chain identity failed')
    return dict(passed=True,original_y2_Z1_cutoff_chain_identity=True,
        original_sigma_Taylor_coefficients_multiplied_by_factorials=True,
        positive_square_root_ordinary_derivative_recurrence='D^alpha(R^2)=2R*D^alpha R+sum_proper binomial(alpha,beta)*D^beta R*D^(alpha-beta)R')


def original_q_jet(roots,eta_log,log_a_lower,loop):
    a=roots['a'];Delta=roots['kappa_minus2'];c=a[ZERO].ctx;scalar=a[ZERO].scalar
    branch=loop['branch'];zero={order:scalar(0) for order in ORDERS}
    if branch=='flat':
        return dict(status='enclosed',rows=zero,branch=branch,cutoff_scope='Delta>=eta: original flat cutoff and every derivative exactly zero',
            active_body_enclosure_evaluated=False)
    if branch!='active':
        return dict(status='requires_source_branch_subdivision',rows=None,branch=branch,
            active_body_enclosure_evaluated=False)
    eta=prior.ScaledEnclosure(prior.FormalScale(a[ZERO].scale.bases,offset=eta_log),1,a[ZERO].ledger)
    if Delta[ZERO].zero or ep(Delta[ZERO].coefficient)[1]<=0:
        cutoff={order:scalar(1 if order==ZERO else 0) for order in ORDERS}
        scope='Delta<=0: original sigma=1, all slow cutoff derivatives exactly zero'
        normalized_argument=None
    else:
        normalized={order:value.positive_divide(eta,eta_log) for order,value in Delta.items()}
        s={order:(scalar(1)-value if order==ZERO else -value) for order,value in normalized.items()}
        try:normalized_argument=bounded(s[ZERO])
        except ArithmeticError:
            return dict(status='requires_source_cutoff_normalization_subdivision',rows=None,branch=branch,
                active_body_enclosure_evaluated=False)
        coefficients=prior.sigma_jets(c,normalized_argument)
        cutoff=compose_sigma(s,[coefficients[n]*math.factorial(n) for n in range(4)])
        scope='original sigma(1-Delta/eta), ordinary chain derivatives through y2 Z1'
    gamma={order:(eta*2 if order==ZERO else scalar(0))-value for order,value in Delta.items()}
    gamma[ZERO]=loop['positive_active_gamma_enclosure']
    ratio=current.quotient_jet(gamma,{order:value*2 for order,value in a.items()},log_a_lower+c.ln(2))
    log_a_upper=a[ZERO].record()['log_absolute_upper']
    # On this same active source box, 2eta-Delta>eta and a<=a_max.
    # Endpoint bounds here are proof coordinates, not selected source values.
    log_root_lower=(eta_log-c.mpf(ep(log_a_upper)[1])-c.ln(2))/2
    ratio[ZERO]=ratio[ZERO].positive_intersection(log_root_lower*2)
    root=sqrt_jet(ratio,log_root_lower)
    rows=current.multiply_jet(cutoff,root)
    return dict(status='enclosed',rows=rows,branch=branch,cutoff_scope=scope,
        normalized_original_cutoff_argument=normalized_argument,active_body_enclosure_evaluated=True,
        positive_active_root_lower_log=log_root_lower,eta_has_no_slow_derivatives=True)


def correlated_excess_rows(roots,log_a_lower):
    """Delta=(a-2)+b^2/a, subtracting the shared constant before tiny addition.

    In original O2 axial/buffer a=2 exactly. Adding a tiny positive quotient
    to2 before subtracting2 can erase its lower bound. This source-equivalent
    order retains the genuine positive b^2/a in its own scale.
    """
    a,b=roots['a'],roots['b']
    quotient=current.quotient_jet(current.multiply_jet(b,b),a,log_a_lower)
    return {order:(value-(2 if order==ZERO else 0))+quotient[order] for order,value in a.items()}


class NativeQSlowJets:
    def __init__(self,owner):
        if type(owner) is not current.NativeCorrelatedShearQ:raise ValueError('Same original correlated shear/q owner required')
        receipt=json.loads((HERE/current.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[current.GATE] or receipt['source_family']!=owner.family:
            raise ValueError('Accepted same-source correlated shear/q required')
        self.owner=owner;self.native=owner.native;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (current.RECEIPT,Path(__file__).name,
            PREFIX+'flat_pulse_derivatives.py',PREFIX+'current_generic_shear_loop.py')})
        self.identity=mixed_chain_identity()

    @native.inlet.source_precision
    def query(self,chart,Z,coordinate):
        source=self.owner.query(chart,Z,coordinate)
        positive=self.owner.owner.decode(self.owner.owner.inventory[chart]['actual_positive_denominator_theorem'])
        eta=packets.interval(self.ctx,self.owner.owner.scales['selected_positive_eta_log'])
        roots=source['roots'];refined=False
        if not chart.startswith('O3_'):
            roots=dict(roots,kappa_minus2=correlated_excess_rows(roots,positive['log_actual_a_positive_lower']))
            loop=current.q_enclosure(roots['a'][ZERO],roots['kappa_minus2'][ZERO],eta,positive['log_actual_a_positive_lower'])
            source_record=dict(source['record'],
                pre_refinement_original_q_branch=source['loop']['branch'],
                pre_refinement_original_q_C0_enclosure=source['q'].record(),
                correlated_shear_and_signed_root_enclosures={key:{'y%d_Z%d'%order:value.record() for order,value in row.items()} for key,row in roots.items()},
                original_q_enclosure={key:value.record() if isinstance(value,prior.ScaledEnclosure) else value for key,value in loop.items()},
                numerical_arithmetic_ledger=dict(roots['a'][ZERO].ledger),
                exact_constant_excess_cancelled_before_tiny_quotient=True,
                original_q_definition_unchanged=True)
            source=dict(source,roots=roots,loop=loop,q=loop['q'],record=source_record);refined=True
        got=original_q_jet(source['roots'],eta,positive['log_actual_a_positive_lower'],source['loop'])
        rows=got['rows'];record={k:v for k,v in got.items() if k!='rows'}
        record.update(source_family=self.family,chart=chart,source_provenance=source['packet'].provenance,
            original_correlated_shear_and_q=source['record'],
            exact_constant_excess_cancelled_before_tiny_quotient=refined,
            current_refined_original_excess_ordinary_rows={'y%d_Z%d'%order:value.record() for order,value in roots['kappa_minus2'].items()},
            current_refined_original_q_C0=source['q'].record(),
            original_q_ordinary_slow_derivative_enclosures=None if rows is None else {'y%d_Z%d'%order:value.record() for order,value in rows.items()},
            original_mixed_derivative_orders=[list(order) for order in ORDERS],
            derivatives_use_original_ordinary_log_radius_rows=True,native_width_and_Pstar_conversion_not_reapplied=True,
            source_point_values_or_zero_derivatives_not_invented=True,
            inverse_phase_or_A_B_slow_jets_installed=False,global_C1_histories_or_Rc_repair_admitted=False,
            **dict.fromkeys(packets.OPEN,False))
        return dict(record=record,rows=rows,source=source)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeQSlowJets(current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))))
        saved=json.loads((HERE/current.NAME).read_bytes())['actual_correlated_shear_and_q_records'];records={}
        for chart,old in saved.items():
            p=old['source_provenance'];Z=packets.interval(owner.ctx,p['Z_box']);v=packets.interval(owner.ctx,p['coordinate_box'])
            records[chart]=owner.query(chart,Z,v)['record']
            print('Original ordinary q slow jets:',chart,records[chart]['status'],flush=True)
        local=owner.query('O2_slope',('.49','.51'),owner.ctx.mpf(('.13369999','.13370001')))
        print('Original q slow jets on actual integral radial/Z cell:',local['record']['status'],flush=True)
        subdivisions={chart:owner.query(chart,('.49','.51'),'.1337')['record'] for chart in ('bridge_first','O2_axial')}
        for chart,record in subdivisions.items():print('Actual original q source Z subdivision:',chart,record['status'],record['branch'],flush=True)
    result=dict(source_family=owner.family,**{GATE:True},native_query_chart_count=len(records),
        actual_native_q_slow_jet_records=records,actual_integral_cell_q_slow_jet_record=local['record'],
        actual_source_Z_subdivision_q_slow_jet_records=subdivisions,
        original_cutoff_mixed_chain_identity=owner.identity,ordinary_derivative_order='y<=2,Z<=1; six rows including y2Z1',
        q_slow_derivative_functions_installed=True,whole_chart_or_implicit_phase_jets_or_C1_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Original six ordinary q slow-jet functions on active/flat accepted native query boxes and the actual O2 integral radial/Z cell; unresolved cutoff/normalization boxes explicitly request subdivision. No whole-chart coverage, phase/A/B slow jets, C1 integral histories, Rc repair/common N/cone/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()

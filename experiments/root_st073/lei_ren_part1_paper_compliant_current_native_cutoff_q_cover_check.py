"""Independent original smooth-cutoff derivatives and source-support checks."""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_cutoff_q_cover as current
import lei_ren_part1_paper_compliant_current_native_phase_first_jets_check as accepted

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
ep=current.ep;require=accepted.require;contains=accepted.contains


def independent_cutoff_checks(c):
    p=accepted.mp.mp.clone();p.dps=110;y,Z=sy.symbols('y Z');comparisons=0;regions={}
    eta=p.mpf('.1');eta_log=c.ln(c.mpf('.1'));allowance=c.mpf(('-1e-70','1e-70'))
    for label,base in (('negative','-.9'),('zero','0'),('transition','.037'),('eta_endpoint','.1'),('flat','.3')):
        bases=tuple(c.mpf(0) for _ in range(5));ledger=accepted.qchecks.new_ledger()
        scalar=lambda value:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),value,ledger)
        aa=sy.Rational(4,5)+sy.Rational(11,100)*y+sy.Rational(7,100)*Z
        dd=sy.Rational(base)+sy.Rational(3,100)*y+sy.Rational(1,25)*Z+sy.Rational(1,20)*y*Z
        roots={name:{order:scalar(str(sy.diff(expr,y,order[0],Z,order[1]).subs({y:0,Z:0})))
            for order in current.ORDERS} for name,expr in (('a',aa),('kappa_minus2',dd))}
        branches,_=current.conditional_cutoff_branches(roots['kappa_minus2'][(0,0)],eta_log)
        def original_q(yy,zz):
            a=p.mpf('.8')+p.mpf('.11')*yy+p.mpf('.07')*zz
            D=p.mpf(base)+p.mpf('.03')*yy+p.mpf('.04')*zz+p.mpf('.05')*yy*zz
            if D>=eta:return p.mpf(0)
            return accepted.original.flat_step(p,1-D/eta)*p.sqrt((2*eta-D)/(2*a))
        references={order:p.diff(original_q,(p.mpf(0),p.mpf(0)),order) for order in current.ORDERS}
        for branch in branches:
            got=current.conditional_q_jet(roots,eta_log,c.ln(c.mpf('.7')),branch)
            for order,value in got['rows'].items():
                require(contains(current.cover.phase.bounded_value(value)+allowance,c.mpf(p.nstr(references[order],115))),
                    'Original cutoff derivative outside '+label+' '+branch['name']+' '+str(order))
                comparisons+=1
            if label=='zero':require(not got['rows'][(0,1)].zero,'Delta=0 must retain genuine sqrt/source derivatives')
            if branch['name']=='flat':require(all(v.zero for v in got['rows'].values()),'Conditional flat jets must be exactly zero')
        regions[label]=dict(branches=[b['name'] for b in branches],comparisons=len(current.ORDERS)*len(branches))
    return dict(passed=True,independent_original_cutoff_ordinary_derivative_comparisons=comparisons,
        regions=regions,original_scalar_function='original flat_step * sqrt((2eta-Delta)/(2a)); lazy zero forDelta>=eta',
        derivative_method='mpmath multivariate diff110 dps; comparison allowance1e-70',
        synthetic_source_fixtures_not_original_field_values=True)


@current.native.inlet.source_precision
def run(original_owner,report=None):
    began=time.monotonic();owner=current.NativeCutoffQCover(original_owner);c=owner.ctx
    report=json.loads((HERE/current.NAME).read_bytes()) if report is None else report
    require(report[current.GATE] and report['source_family']==owner.family,'Same original cutoff producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed original prerequisite: '+name)
    independent=independent_cutoff_checks(c);queries={}
    for chart in ('bridge_first','O2_axial'):
        got=owner.query(chart,(-1,1),'.1337');queries[chart]=got['record']
        require(got['record']['status']=='enclosed' and got['source']['q'] is got['rows'][(0,0)],'Original q/jet C0 consistency required')
        require(all(branch['query']['source']['q'] is branch['query']['rows'][(0,0)] for branch in got['branches']),
            'Conditional q/jet C0 consistency required')
        require(all(v.ctx is c and v.ledger is got['source']['q'].ledger for v in got['rows'].values()),'Cutoff branch union changed source arithmetic')
    sample=owner.query('O2_slope',('.5','.5'),'.1337')['source']['q']
    # Huge/very small logs exercise only range arithmetic. They are not
    # native A values. The native source theorem licenses the bounded factor.
    preserved=[]
    for sign,log in ((1,'-100000'),(-1,'-100000'),(1,'100000'),(-1,'100000')):
        coefficient=sign if log.startswith('-') else c.mpf((0,1) if sign>0 else (-1,0))
        value=current.prior.ScaledEnclosure(current.prior.FormalScale(sample.scale.bases,offset=c.mpf(log)),coefficient,sample.ledger)
        got=owner.supported_expm1(value,160)
        require(got.scale.powers==value.scale.powers and got.scale.offset._mpi_==value.scale.offset._mpi_,
            'Supported expm1 must retain the original formal increment factor')
        signed=(ep(got.coefficient)[0]>=0 and ep(got.coefficient)[1]>0) if sign>0 else (ep(got.coefficient)[1]<=0 and ep(got.coefficient)[0]<0)
        require(got.ledger is value.ledger and signed,
            'Supported expm1 must retain original sign and ledger')
        preserved.append(dict(sign=sign,log_coordinate=log,original_formal_factor_retained=True))
    zero=sample.scalar(0);require(owner.supported_expm1(zero,160) is zero,'Exact zero increment must remain exact')
    impossible=current.prior.ScaledEnclosure(current.prior.FormalScale(sample.scale.bases,offset=c.mpf('100000')),1,sample.ledger)
    try:owner.supported_exponential(impossible,160)
    except ArithmeticError:pass
    else:raise AssertionError('Contradictory source exponent cannot be hidden by a cap')
    for bad in (0,159,True):
        try:owner.supported_exponential(sample,bad)
        except ValueError:pass
        else:raise AssertionError('Invalid original candidate N accepted')
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},
        independent_original_cutoff_derivative_checks=independent,additional_actual_original_full_Z_cutoff_queries=queries,
        original_primitive_identity_checks=current.primitive_identity_checks(),
        supported_expm1_factor_preservation_checks=preserved,range_arithmetic_fixtures_not_native_A_values=True,
        original_q_and_jet_ZERO_same_object=True,original_primitive_A_or_A_Z_not_replaced_by_cap=True,
        actual_density_or_global_histories_controls_or_recursion_admitted=False,
        **dict.fromkeys(current.packets.OPEN,False),execution_seconds=time.monotonic()-began,
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='Original cutoff source jets on two full-Z boxes, independent mixed derivatives through y2/Z1 at cutoff boundaries, exact primitive/support identities and formal increment factor preservation. Branchwise actual density and full-route closure remain separate.')
    (HERE/current.RECEIPT).write_text(json.dumps(current.packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original cutoff q cover focused check PASS:',independent['independent_original_cutoff_ordinary_derivative_comparisons'],flush=True)
    return result

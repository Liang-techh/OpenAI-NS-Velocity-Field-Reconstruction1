"""Original smooth-cutoff q jets on negative, transition and flat subdomains.

Conditional C0 ranges preserve the original a/Delta derivative rows. All
possible cutoff branches are retained; branch hulls do not define a field.
The same source support also gives the exact primitive bound |A|<=5/4.
"""
import json
import math
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_signed_u_phase_cover as cover

slow=cover.first.slow;current=slow.current;prior=slow.prior
native=slow.native;packets=slow.packets;ep=slow.ep;ORDERS=slow.ORDERS;ZERO=slow.ZERO
HERE,PREFIX,sha=slow.HERE,slow.PREFIX,slow.sha
NAME=PREFIX+'current_native_cutoff_q_cover.json'
RECEIPT=PREFIX+'current_native_cutoff_q_cover_check.json'
GATE='current_original_conditional_cutoff_q_slow_jet_cover_executed'


def signed_log_upper(value):
    hi=ep(value.coefficient)[1]
    return None if hi<=0 else ep(value.scale.evaluate()+value.ctx.ln(value.ctx.mpf(hi)))[1]


def signed_log_lower(value):
    lo=ep(value.coefficient)[0]
    return None if lo<=0 else ep(value.scale.evaluate()+value.ctx.ln(value.ctx.mpf(lo)))[0]


def conditional_cutoff_branches(Delta,eta_log):
    c=Delta.ctx;lo,hi=ep(Delta.coefficient);eta_lo,eta_hi=ep(eta_log)
    lower=signed_log_lower(Delta);upper=signed_log_upper(Delta);branches=[];empty=[]
    if lo<=0:
        # a>0 and kappa=a+b^2/a imply Delta>-2 on the true source.
        if Delta.zero:negative=c.mpf(0)
        else:
            try:negative=cover.phase.clipped(c,Delta.coefficient*Delta.bounded_exp(Delta.scale.evaluate()),-2,0)
            except ArithmeticError:negative=c.mpf((-2,0))
        branches.append(dict(name='negative',condition='Delta<=0',Delta=Delta.scalar(negative)))
    else:empty.append(dict(name='negative',proof='original Delta range strictly positive'))
    if hi<0 or (lower is not None and lower>eta_hi):
        empty.append(dict(name='transition',proof='original Delta range disjoint from[0,eta]'))
    else:
        theta=c.mpf((0,1))
        if Delta.zero:theta=c.mpf(0)
        else:
            try:
                ratio=prior.ScaledEnclosure(Delta.scale-prior.FormalScale(Delta.scale.bases,offset=eta_log),Delta.coefficient,Delta.ledger)
                theta=cover.phase.clipped(c,ratio.coefficient*ratio.bounded_exp(ratio.scale.evaluate()),0,1)
            except ArithmeticError:pass
        eta=prior.ScaledEnclosure(prior.FormalScale(Delta.scale.bases,offset=eta_log),1,Delta.ledger)
        branches.append(dict(name='transition',condition='0<=Delta<=eta',Delta=eta*theta,theta=theta))
    if hi<=0 or (upper is not None and upper<eta_lo):
        empty.append(dict(name='flat',proof='original Delta range strictly below positive eta'))
    else:branches.append(dict(name='flat',condition='Delta>=eta',Delta=Delta))
    if not branches:raise ArithmeticError('Original cutoff cover unexpectedly empty')
    return branches,empty


def conditional_q_jet(roots,eta_log,log_a_lower,branch):
    a=roots['a'];c=a[ZERO].ctx;scalar=a[ZERO].scalar
    if branch['name']=='flat':
        rows={order:scalar(0) for order in ORDERS}
        return dict(status='enclosed',rows=rows,cutoff_scope='Conditional Delta>=eta; original q and every derivative exactly zero')
    eta=prior.ScaledEnclosure(prior.FormalScale(a[ZERO].scale.bases,offset=eta_log),1,a[ZERO].ledger)
    Delta=dict(roots['kappa_minus2']);Delta[ZERO]=branch['Delta']
    if branch['name']=='negative':
        cutoff={order:scalar(1 if order==ZERO else 0) for order in ORDERS}
        gamma={order:(eta*2 if order==ZERO else scalar(0))-value for order,value in Delta.items()}
        gamma_lower=eta_log+c.ln(2);scope='Conditional Delta<=0; original sigma=1; original sqrt jets remain nonconstant'
    else:
        # Collect the defining eta/eta cancellation in the C0 geometry.
        # Nonzero theta derivatives remain the original Delta derivatives/eta.
        theta={order:(scalar(branch['theta']) if order==ZERO else value.positive_divide(eta,eta_log))
            for order,value in Delta.items()}
        s={order:(scalar(1)-value if order==ZERO else -value) for order,value in theta.items()}
        coefficients=prior.sigma_jets(c,1-branch['theta'])
        cutoff=slow.compose_sigma(s,[coefficients[n]*math.factorial(n) for n in range(4)])
        gamma={order:(eta*(2-branch['theta']) if order==ZERO else -value) for order,value in Delta.items()}
        gamma_lower=eta_log;scope='Conditional original sigma(1-theta), theta=Delta/eta; original slow derivative rows'
    gamma[ZERO]=gamma[ZERO].positive_intersection(gamma_lower)
    ratio=current.quotient_jet(gamma,{order:value*2 for order,value in a.items()},log_a_lower+c.ln(2))
    log_a_upper=c.mpf(ep(a[ZERO].record()['log_absolute_upper'])[1])
    log_root_lower=(gamma_lower-log_a_upper-c.ln(2))/2
    ratio[ZERO]=ratio[ZERO].positive_intersection(log_root_lower*2)
    root=slow.sqrt_jet(ratio,log_root_lower)
    rows=current.multiply_jet(cutoff,root)
    return dict(status='enclosed',rows=rows,cutoff_scope=scope,positive_conditional_gamma_lower_log=gamma_lower,
        positive_conditional_root_lower_log=log_root_lower,
        original_a_and_Delta_nonzero_derivative_rows_unchanged=True,
        original_sigma_Taylor_rows_multiplied_by_factorials=True,
        geometric_eta_over_eta_cancelled_before_interval_evaluation=branch['name']=='transition')


def primitive_identity_checks():
    a,t,q,h,s,W1,W2,x,nu=sy.symbols('a t q h s W1 W2 x nu',nonzero=True)
    nu0=1+t*t+2*q*q
    F=((1+t*t)*2*sy.pi*x+4*t*q*h*W1+4*q*q*s*W2)/(nu0*2*sy.pi)
    A=a*(4*t*q*h*W1+q*q*(4*s*W2-4*sy.pi*x))/(nu0*4*sy.pi)
    if sy.simplify(A-a*(F-x)/2)!=0:raise ArithmeticError('Original small-chart primitive identity failed')
    r,E=sy.symbols('r E',nonzero=True)
    F=((1+t*t)*x+2*t*q*h/r*(E-x)+q*q/r**2*((2-3*s)*E+s*x+r/sy.pi*sy.sin(2*sy.pi*E)))/nu0
    A=a/2*(2*t*q*h/r*(E-x)+q*q/r**2*((2-3*s)*(E-x)+r/sy.pi*sy.sin(2*sy.pi*E)))/nu0
    if sy.simplify((A-a*(F-x)/2).subs(s,1-r*r))!=0:raise ArithmeticError('Original Mobius-chart primitive identity failed')
    return dict(passed=True,original_small_and_Mobius_primitive_identities=2,
        original_identity='A=a/2*(Phi-psi_fraction); original inverse Phi=phi',
        support_proof='q!=0 => Delta<eta; a>0 => a<=kappa=2+Delta<2+eta<=5/2; phi,psi_fraction in[0,1]',
        conclusion='|A|<=5/4 globally, including exact q=0; N>=160 implies |A/N|<=1/128',
        constant_bound_is_outer_range_not_defining_field_value=True)


class NativeCutoffQCover:
    def __init__(self,owner):
        if type(owner) is not slow.NativeQSlowJets:raise TypeError('Existing accepted original q-jet owner required')
        receipt=json.loads((HERE/slow.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(slow.GATE) or receipt['source_family']!=owner.family:
            raise ValueError('Accepted same-original-family q jets required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service;self.native=owner.native
        self.eta_log=packets.interval(self.ctx,owner.owner.owner.scales['selected_positive_eta_log'])
        if ep(self.eta_log)[1]>ep(-self.ctx.ln(2))[0]:raise ValueError('Original selected eta<=1/2 theorem required')
        checked=json.loads((HERE/cover.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(cover.GATE) or checked['source_family']!=self.family:raise ValueError('Checked original directed source context bridge required')
        self.hashes={**receipt['input_hashes'],**checked['input_hashes'],slow.RECEIPT:sha(slow.RECEIPT),
            cover.RECEIPT:sha(cover.RECEIPT),Path(__file__).name:sha(Path(__file__).name)}
        self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def query(self,chart,Z,coordinate):
        query=cover.bind_source_context(self.owner.query(chart,Z,coordinate),self.ctx)
        source=query['source'];roots=source['roots'];Delta=roots['kappa_minus2'][ZERO]
        proof=self.owner.owner.owner.decode(self.owner.owner.owner.inventory[chart]['actual_positive_denominator_theorem'])
        candidates,empty=conditional_cutoff_branches(Delta,self.eta_log);branches=[]
        for branch in candidates:
            got=conditional_q_jet(roots,self.eta_log,proof['log_actual_a_positive_lower'],branch)
            rows=got['rows'];restricted=dict(roots,kappa_minus2=dict(roots['kappa_minus2']))
            restricted['kappa_minus2'][ZERO]=branch['Delta']
            record={key:value for key,value in got.items() if key!='rows'}
            record.update(conditional_branch=branch['name'],condition=branch['condition'],source_family=self.family,
                original_conditional_q_slow_jet_enclosures={'y%d_Z%d'%order:value.record() for order,value in rows.items()},
                original_unconditioned_source_record_reference='original_unconditioned_q_source in parent cutoff record',
                original_source_packet_provenance=source['packet'].provenance,
                original_a_Delta_nonzero_derivative_rows_identical_to_parent_query=True,
                conditional_C0_range_not_new_source_field=True,branch_boundary_derivatives_not_introduced=True,
                original_q_and_jet_ZERO_same_object=True,**dict.fromkeys(packets.OPEN,False))
            branch_query=dict(query,rows=rows,source=dict(source,roots=restricted,q=rows[ZERO]),record=record)
            branches.append(dict(record=record,query=branch_query))
        union=cover.density.local.same_source_union
        rows={order:union([branch['query']['rows'][order] for branch in branches]) for order in ORDERS}
        record=dict(status='enclosed',chart=chart,source_family=self.family,
            original_unconditioned_q_source=query['record'],original_eta_log=self.eta_log,
            cutoff_branch_cover='Delta<=0 union 0<=Delta<=eta union Delta>=eta covers R',
            conditional_branches=[branch['record'] for branch in branches],branches_proved_empty=empty,
            q_C0_and_first_mixed_ordinary_rows={'y%d_Z%d'%order:value.record() for order,value in rows.items()},
            source_q_and_jet_ZERO_same_object=True,conditional_density_must_precede_cutoff_branch_union=True,
            aggregate_source_roots_are_unconditioned_packet_metadata=True,
            nonlinear_phase_and_density_must_use_conditional_branch_query_objects=True,
            global_histories_controls_or_recursion_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(query,source=dict(source,q=rows[ZERO]),rows=rows,record=record,branches=branches)

    def supported_exponent_argument(self,value,N):
        """Bound only an exponential factor; leave the original A/N formal."""
        N=cover.density.density.candidate_integer(N)
        if N<160:raise ValueError('Original source support uses integer N>=160')
        c=value.ctx;cap=c.mpf(5)/(4*N)
        if value.zero:return c.mpf(0)
        lo,hi=ep(value.coefficient)
        if lo>0 or hi<0:
            lower=ep(value.scale.evaluate()+c.ln(c.mpf(min(abs(lo),abs(hi)))))[0]
            if lower>ep(c.ln(cap))[1]:raise ArithmeticError('Original source exponent range contradicts the primitive support theorem')
        try:finite=value.coefficient*value.bounded_exp(value.scale.evaluate())
        except ArithmeticError:return c.mpf((-ep(cap)[1],ep(cap)[1]))
        return cover.phase.clipped(c,finite,-cap,cap)

    def supported_expm1(self,value,N):
        c=value.ctx;x=self.supported_exponent_argument(value,N)
        if value.zero:return value
        lo,hi=ep(x)
        factor=c.mpf((ep(c.exp(c.mpf(min(lo,0))))[0],ep(c.exp(c.mpf(max(hi,0))))[1]))
        return value*factor

    def supported_exponential(self,value,N):
        return value.ctx.exp(self.supported_exponent_argument(value,N))


@native.inlet.source_precision
def run(original_owner):
    began=time.monotonic();owner=NativeCutoffQCover(original_owner)
    records={chart:owner.query(chart,(-1,1),'.1337')['record'] for chart in ('bridge_first','O2_axial')}
    result=dict(source_family=owner.family,**{GATE:True},actual_original_full_Z_cutoff_queries=records,
        original_primitive_support_bound=primitive_identity_checks(),original_eta_log=owner.eta_log,
        ordinary_mixed_q_orders=[list(order) for order in ORDERS],actual_original_query_count=2,
        source_function_ranges_not_midpoint_values=True,actual_density_or_integral_or_controls_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=owner.hashes,
        scope='Original conditional q ordinary y2/Z1 jets on two formerly unresolved full-Z native boxes, plus exact original primitive/support identity. Branch-local density, full-route integrals, actual controls and recursion remain separate.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result

"""Source-sector component budgets for the genuine canonical n1 series.

Normalize bounds by an original monomial and a directed coefficient norm,
not by a selected field value. Subtract exponent tuples before combining
logs. Huge positive relative tails stay logarithmic; no F0 is resolved.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_local_leading_dispatch as leading
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

selection=leading.selection;coupled=selection.coupled;op=leading.op;domain=leading.domain
HERE,PREFIX,sha=leading.HERE,leading.PREFIX,leading.sha
NAME=PREFIX+'current_original_n1_component_budget.json.gz'
RECEIPT=PREFIX+'current_original_n1_component_budget_check.json'
GATE='current_original_n1_source_sector_component_error_budgets_computed'
OPEN=leading.OPEN;ep=leading.ep


def exponent_difference(left,right):return tuple(a-b for a,b in zip(left,right))


def log_ratio(c,logs,powers,reference,magnitude,reference_norm):
    # Equal huge source logs cancel as integer exponents, before arithmetic.
    diff=exponent_difference(powers,reference)
    return sum((log*p for log,p in zip(logs,diff) if p),c.mpf(0))+c.ln(magnitude/reference_norm)


def exp_cover(c,logarithm,floor_digits=80):
    """Directed norm upper; never exponentiate an enormous source log."""
    floor=c.mpf(10)**(-floor_digits)
    if ep(logarithm)[1]<=ep(c.ln(floor))[0]:return floor,'directed_small_log_cover'
    if ep(logarithm)[1]>1500:return None,'positive_relative_scale_retained_logarithmically'
    return domain.upper(c,c.exp(c.mpf(ep(logarithm)[1]))),'directed_exponential_norm_upper'


def component_budget(value,k,tail_norm,tail_powers=(0,0,0,0)):
    c=value.ctx;logs=value.algebra.logs
    if k>value.order:raise ValueError('Budget requires a genuine available derivative')
    rows={p:jet[k] for p,jet in value.terms.items() if ep(jet[k])!=(0,0)}
    if not rows:
        return dict(prefix_exact_zero=True,relative_width_certified=False,
            total_exact_zero=ep(tail_norm)[1]==0,remaining_tail_norm=tail_norm,
            reference_source_exponents=None,source_and_cell_uncertainty_separated=False)
    # Endpoint mp scalars can retain more bits than the global mp context.
    # Taking their scalar abs first would round before directed conversion.
    norm=lambda v:domain.upper(c,abs(v))
    reference=next(iter(rows));M=norm(rows[reference])
    for p,v in rows.items():
        if ep(log_ratio(c,logs,p,reference,norm(v),M))[0]>0:reference=p;M=norm(v)
    total=c.mpf(0);cells=[];unresolved=False
    for p,v in sorted(rows.items()):
        lo,hi=ep(v);limits=[];rules=[]
        for endpoint in (lo,hi):
            if endpoint==0:limits.append(c.mpf(0));rules.append('exact_zero_coefficient');continue
            logbound=log_ratio(c,logs,p,reference,abs(c.mpf(endpoint)),M)
            upper,rule=exp_cover(c,logbound);rules.append(rule)
            if upper is None:limits.append(None);continue
            # Small-log bounds only supply an upper magnitude, hence the
            # signed interval includes zero on that side.
            if rule=='directed_small_log_cover':
                covered=c.mpf([0,ep(upper)[1]])
                limits.append(-covered if endpoint<0 else covered)
            else:
                low=c.mpf(0) if ep(logbound)[0]<-1500 else domain.lower(c,c.exp(c.mpf(ep(logbound)[0])))
                exact=c.mpf([ep(low)[0],ep(upper)[1]])
                limits.append(-exact if endpoint<0 else exact)
        if any(v is None for v in limits):
            interval=None;unresolved=True
        else:
            interval=c.mpf([ep(limits[0])[0],ep(limits[1])[1]]);total+=interval
        cells.append(dict(source_exponents=list(p),reference_exponent_difference=list(exponent_difference(p,reference)),
            signed_normalized_coefficient_enclosure=interval,endpoint_bound_rules=rules))
    tail_log=None;tail_upper=c.mpf(0);tail_rule='exact_zero_tail'
    if ep(tail_norm)[1]>0:
        tail_log=log_ratio(c,logs,tail_powers,reference,domain.upper(c,tail_norm),M)
        tail_upper,tail_rule=exp_cover(c,tail_log)
    prefix_width=None
    if not unresolved and (ep(total)[0]>0 or ep(total)[1]<0):
        prefix_width=domain.upper(c,(total.b-total.a)/domain.lower(c,abs(total)))
    full=None;relative_width=None;passed=False
    if not unresolved and tail_upper is not None:
        full=total+c.mpf([-1,1])*tail_upper
        lo,hi=ep(full)
        if lo>0 or hi<0:
            lower=domain.lower(c,abs(full));relative_width=domain.upper(c,(full.b-full.a)/lower)
            passed=ep(relative_width)[1]<=c.mpf('.001')
    return dict(prefix_exact_zero=False,reference_source_exponents=list(reference),
        directed_reference_coefficient_norm=M,reference_is_a_bound_normalization_not_a_selected_field_value=True,
        signed_prefix_normalized_enclosure=None if unresolved else total,source_sector_contributions=cells,
        remaining_tail_source_exponents=list(tail_powers),remaining_tail_normalized_log_upper=tail_log,
        remaining_tail_normalized_norm_upper=tail_upper,remaining_tail_bound_rule=tail_rule,
        prefix_width_over_abs_lower=prefix_width,
        full_series_normalized_enclosure=full,total_width_over_abs_lower=relative_width,
        relative_total_width_target='.001',relative_width_certified=passed,
        source_and_cell_uncertainty_separated=False,
        source_amplitude_materialized=False,exact_exponent_subtraction_before_log_arithmetic=True)


@dataclass(frozen=True,eq=False)
class OriginalN1ComponentBudgetPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1ComponentBudget:
    @source_precision
    def __init__(self,local,require_checked=True):
        if type(local) is not leading.CurrentOriginalN1LocalLeadingDispatch or not local.acceptance_loaded:
            raise ValueError('Accepted original local leading dispatcher required')
        local.assert_graph();self.local=local;self.inner=local.inner;self.base=self.inner.base
        self.c=local.c;self.family=copy.deepcopy(local.family);self.hashes=dict(local.hashes);self._packets={}
        bind=domain.parent.parent.original.bind
        for name in (leading.NAME,leading.RECEIPT,Path(__file__).name):bind(self.hashes,name,sha(name))
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Component budget receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Component budget dependency changed '+name)
                bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._local=local;self._inner=self.inner;self._base=self.base;self._c=self.c
        self._family=copy.deepcopy(self.family);self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.local is not self._local or self.inner is not self._inner or self.base is not self._base
            or self.c is not self._c or self.inner is not self.local.inner or self.base is not self.inner.base
            or not self.local.acceptance_loaded or self.family!=self._family or self.family!=self.local.family
            or self.hashes!=self._hashes):raise ValueError('Component budget source graph changed')
        return self.local.assert_graph()

    @source_precision
    def evaluate(self,inner_packet):
        self.assert_graph();c=self.c;view=self.inner.report(inner_packet)
        child=self.inner._packets[id(inner_packet)][2];state=self.base._packets[id(child)][2]
        field=view['canonical_n1_field'];query=view['query'];z=field['Z'];delta=self.base.source.core.delta
        local_packet=self.local.evaluate(z,query['coordinate'],query['chart']);local=self.local.report(local_packet)
        geometry=local['actual_leading_source']['geometry'];rho=geometry.get('rho',geometry.get('rho_enclosure_only'))
        algebra=self.base.inputs.frame(z)[1];radial=coupled.recover_radial(algebra,state,z,rho,delta)
        tail=field['analytic_preconditioned_tail'];E=tail['omitted_sum_bound'];E1=field['derivative_tail']['remaining_first_Z_derivative_norm_upper']
        if query['chart']=='core' and ep(query['coordinate'])[1]==0:E=E1=c.mpf(0)
        Lmin=domain.lower(c,1-delta*z*z)
        if ep(Lmin)[0]<=0:raise ValueError('Real radial tail denominator must be positive')
        Qerror=domain.upper(c,((abs(1-delta)+abs(1+delta))*abs(z)*E+2*abs(1-z*z)*E1)/Lmin)
        budgets={name:{'Z0':component_budget(value,0,E),'Z1':component_budget(value,1,E1)} for name,value in state.items()}
        budgets['mean1']={'Z0':component_budget(radial['mean1'],0,2*E),'Z1':component_budget(radial['mean1'],1,2*E1)}
        budgets['Q1']={'Z0':component_budget(radial['Q1'],0,Qerror)}
        budgets['V1']={'Z0':component_budget(radial['V1'],0,domain.upper(c,rho*Qerror),(0,0,0,-1))}
        count=sum(row['relative_width_certified'] for values in budgets.values() for row in values.values())
        result=dict(source_family=self.family,Z=z,query=query,component_budgets=budgets,
            passing_relative_total_width_rows=count,total_component_derivative_rows=sum(map(len,budgets.values())),
            omitted_W_tail_norm=E,omitted_dZ_W_tail_norm=E1,omitted_Q1_real_query_tail_norm=Qerror,
            radial_tail_scope='real queried Z; W/dZ W tail comes from the common eta/16 analytic tube',
            local_leading_branch_bound=True,combined_source_and_cell_enclosure_reported_separately_from_series_tail=True,
            source_uncertainty_and_integration_error_still_combined=True,sector_preserving_analytic_tail_still_open=True,
            physical_point_values_not_selected=True,whole_Z_or_full_inner_solution_not_promoted=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalN1ComponentBudgetPacket();self._packets[id(packet)]=(packet,result,inner_packet,local_packet,op.snapshot(result))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1ComponentBudgetPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live source component budget packet required')
        self.inner.report(entry[2]);self.local.report(entry[3])
        if op.snapshot(entry[1])!=entry[4]:raise ValueError('Component budget packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(local,inner_packet):
    began=time.monotonic();owner=CurrentOriginalN1ComponentBudget(local,require_checked=False)
    packet=owner.evaluate(inner_packet);view=owner.report(packet)
    raw=dict(source_family=owner.family,samples={'selected_actual_inner_exit':view},input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_COMPONENT_ERROR_BUDGETS_COMPUTED',view['passing_relative_total_width_rows'],view['total_component_derivative_rows'],flush=True)
    return owner,packet

"""Original correlated swirl source integral through the actual first collar.

Core errors, homogeneous core-to-exit transport and the actual normalized
collar source are kept separate. No remaining coupled feedback is solved.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source as core
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

weighted=core.weighted;op=core.op;domain=core.domain;ep=core.ep
HERE,PREFIX,sha=core.HERE,core.PREFIX,core.sha
NAME=PREFIX+'current_original_n1_correlated_inner_source.json.gz'
RECEIPT=PREFIX+'current_original_n1_correlated_inner_source_check.json'
GATE='current_original_n1_correlated_swirl_initial_source_integral_on_selected_inner_interval_enclosed'
OPEN=core.OPEN;FIELDS=core.FIELDS


def vector_norm(values,errors):
    """Per-sector infinity norm for the two-component core Cauchy data."""
    algebra=values[FIELDS[0]].algebra;c=algebra.ctx
    rows={field:weighted.absolute_jet(values[field])+errors[field] for field in FIELDS}
    powers=set().union(*(v.terms for v in rows.values()))
    return op.FactoredJet(algebra,{p:op.Jet(c,[domain.upper(c,c.mpf(max(
        ep(v.terms[p][k])[1] if p in v.terms else 0 for v in rows.values()))) for k in (0,1)]) for p in powers},1)


def homogeneous_transfer_error(c,vector_bound,chi_Z,exit_x):
    """Cauchy-data difference from x=2 to X, including the ordinary Z1 jet.

    A=[[0,1],[-2chi,-3/x]], ||A||infinity<=7/2 for x>=2, |chi|<=1.
    A_Z has norm<=2|chi_Z|. X is independent of Z at fixed original phase.
    """
    if ep(exit_x)[0]<2 or ep(exit_x)[1]>ep(c.sqrt(c.mpf('4.1')))[1]:
        raise ValueError('Original first-collar exit x in[2,sqrt(4.1)] required')
    length=domain.upper(c,c.mpf([0,ep(exit_x-2)[1]]));C=c.mpf(7)/2;D=2*chi_Z
    exponential=domain.upper(c,c.exp(C*length))
    gain=op.Jet(c,[domain.upper(c,C*length*exponential),
                   domain.upper(c,D*length*exponential*(1+C*length))])
    return vector_bound*gain,dict(length_upper=length,matrix_infinity_norm_upper=C,
        matrix_Z_derivative_norm_upper=D,exponential_upper=exponential,
        original_phase_radius_independent_of_Z=True,
        proof='integrate homogeneous 2x2 Cauchy system and its Z derivative; bound exp(C*l)-1 by C*l*exp(C*l)')


@dataclass(frozen=True,eq=False)
class OriginalN1CorrelatedInnerSourcePacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1CorrelatedInnerSource:
    @source_precision
    def __init__(self,accepted_core,require_checked=True):
        if type(accepted_core) is not core.CurrentOriginalN1AnalyticCoreSource or not accepted_core.acceptance_loaded:
            raise ValueError('Accepted original analytic core source owner required')
        accepted_core.assert_graph();self.core=accepted_core;self.weighted=accepted_core.weighted
        self.inputs=accepted_core.inputs;self.local=self.weighted.local;self.c=accepted_core.c
        self.family=copy.deepcopy(accepted_core.family);self.hashes=dict(accepted_core.hashes)
        for name in (core.NAME,core.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._core=accepted_core;self._weighted=self.weighted;self._inputs=self.inputs;self._local=self.local;self._c=self.c
        self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Correlated inner source receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Correlated inner source dependency changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.core is not self._core or self.weighted is not self._weighted or self.inputs is not self._inputs
            or self.local is not self._local or self.c is not self._c or self.weighted is not self.core.weighted
            or self.inputs is not self.core.inputs or self.local is not self.weighted.local
            or not self.core.acceptance_loaded or self.family!=self.core.family or self.hashes!=self._hashes):
            raise ValueError('Correlated inner source graph changed')
        return self.core.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate,collar_cells=4,N=12):
        self.assert_graph();c=self.c;fraction=c.mpf(coordinate)
        if not isinstance(coordinate,(str,int,float)) and ep(fraction)[0]!=ep(fraction)[1]:raise ValueError('One exact outer fraction required')
        if ep(fraction)[0]<0 or ep(fraction)[1]>1:raise ValueError('Selected completed-leading inner fraction in[0,1] required')
        if isinstance(collar_cells,bool) or not isinstance(collar_cells,int) or not 1<=collar_cells<=128:raise ValueError('Bounded collar partition required')
        weighted.kernel_polynomials(N);began=time.monotonic()
        lp=self.local.evaluate(Z,coordinate,'first');lv=self.local.report(lp);z=lv['Z'];query=lv['query']
        phase=query['actual_source_coordinate'];p=self.inputs.evaluate(z,phase,'first');v=self.inputs.report(p)
        if v['geometry']['exact_rho']!='4*exp(hb*s)':raise ValueError('Original fixed-phase rho identity required')
        rho=v['geometry']['rho_enclosure_only'];X=c.sqrt(rho)
        cp=self.core.evaluate(z,'4');cv=self.core.report(cp);centry=self.core._packets[id(cp)]
        values,core_errors=centry[3:5];algebra=self.inputs._packets[id(p)][5]
        if any(val.algebra is not algebra for val in values.values()):raise ValueError('Original amplitude anchors differ')
        axis=self.weighted.physical.axis_inputs(z);chi=axis['chi'].truncate(1)
        chi_Z=domain.upper(c,abs(chi[1]));gradient=op.Jet(c,self.weighted.source.amplitude.gradient_jets(z))
        transfer,transfer_proof=homogeneous_transfer_error(c,vector_norm(values,core_errors),chi_Z,X)
        collar={field:algebra.lift(0,1) for field in FIELDS};tails={field:algebra.lift(0,1) for field in FIELDS}
        kernel_tail=weighted.potential_kernel_tail(c,self.weighted.a,chi_Z,N)
        cells=[];sources=[]
        if ep(phase)[1]>0:
            for j in range(collar_cells):
                s=phase*weighted.budget.coupled.partition(c,j,collar_cells)
                cell=self.inputs.evaluate(z,s,'first');sv=self.inputs.report(cell);entry=self.inputs._packets[id(cell)]
                rad=sv['geometry']['rho_enclosure_only'];x=c.sqrt(rad)
                source=core.normalized_theta_sectors(algebra,entry[2]['Phi'],algebra.lift(op.Jet.variable(c,z,4)),
                                                     rad,self.weighted.physical.delta,algebra.lift(gradient))
                g=sum(source.values(),algebra.lift(0,1))
                ratio=x/X;t=c.mpf([max(0,ep(ratio)[0]),min(1,ep(ratio)[1])])
                kernels=dict(zip(FIELDS,weighted.potential_kernel(chi,t,X,N,1)))
                measure=x*phase/(2*collar_cells)
                for field in FIELDS:
                    actual=algebra.width(g*kernels[field]*measure)
                    omitted=algebra.width(weighted.absolute_jet(g)*kernel_tail[field]*measure)
                    if any(p[0]<0 for val in (actual,omitted) for p in val.terms):raise ValueError('Actual collar source inverse hb did not cancel')
                    collar[field]+=op.truncate(actual,1);tails[field]+=op.truncate(omitted,1)
                cells.append(cell);sources.append(dict(phase_cell=s,rho_enclosure=rad,
                    normalized_theta_sectors={label:op.export(val,1) for label,val in source.items()},
                    actual_normalized_Phi_used_before_F0_dressing=True,dx_measure='hb*x*ds/2',inverse_hb_cancelled_before_enclosure=True))
        sums={field:values[field]+collar[field] for field in FIELDS}
        errors={field:core_errors[field]+transfer+tails[field] for field in FIELDS}
        physical={field:weighted.transform(val,gradient,(0,1,0,1),1) for field,val in sums.items()}
        physical_errors={field:algebra.shift(weighted.absolute_jet(weighted.relative_amplitude(algebra,gradient,1,1))*val,(0,1,0,1)) for field,val in errors.items()}
        diagnostics={field:{'Z'+str(k):core.integral_budget(sums[field],errors[field],k) for k in (0,1)} for field in FIELDS}
        physical_budgets={field:{'Z'+str(k):core.integral_budget(physical[field],physical_errors[field],k) for k in (0,1)} for field in FIELDS}
        result=dict(source_family=self.family,Z=z,query=query,actual_exit_geometry=v['geometry'],
            normalized_correlated_initial_integrals={field:op.export(val,1) for field,val in sums.items()},
            normalized_core_contribution={field:op.export(val,1) for field,val in values.items()},
            normalized_actual_collar_contribution={field:op.export(val,1) for field,val in collar.items()},
            normalized_core_errors={field:op.export(val,1) for field,val in core_errors.items()},
            normalized_homogeneous_core_to_exit_error=op.export(transfer,1),homogeneous_transfer_proof=transfer_proof,
            normalized_collar_kernel_tail={field:op.export(val,1) for field,val in tails.items()},
            normalized_total_errors={field:op.export(val,1) for field,val in errors.items()},
            physical_correlated_initial_integrals={field:op.export(val,1) for field,val in physical.items()},
            physical_total_errors={field:op.export(val,1) for field,val in physical_errors.items()},
            normalized_initial_integral_budgets=diagnostics,physical_initial_integral_budgets=physical_budgets,
            passing_normalized_initial_width_rows=sum(row['relative_initial_integral_width_certified'] for f in diagnostics.values() for row in f.values()),
            actual_collar_source_cells=sources,core_checkpoint_gate=cv[core.GATE],
            ordinary_Z_initial_derivative_order=1,
            original_actual_collar_not_comparison_or_Bessel_model=True,
            original_normalized_Phi_derivatives_extracted_before_amplitude_dressing=True,
            original_width_cancellation_and_amplitude_anchors_retained=True,
            source_quadrature_enclosure_and_kernel_tail_both_included=True,
            core_transport_source_radius_enclosure_not_selected_endpoint=True,
            initial_source_integral_only=True,remaining_normalized_B0_B1_feedback_not_solved=True,
            old_coupled_prefix_not_replaced_by_initial_integral=True,whole_Z_or_full_inner_solution_not_promoted=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1CorrelatedInnerSourcePacket()
        self._packets[id(token)]=(token,result,lp,p,cp,cells,sums,errors,op.snapshot((result,sums,errors)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1CorrelatedInnerSourcePacket or entry is None or entry[0] is not packet:raise ValueError('Live correlated source packet required')
        self.local.report(entry[2]);self.inputs.report(entry[3]);self.core.report(entry[4])
        for cell in entry[5]:self.inputs.report(cell)
        if op.snapshot((entry[1],entry[6],entry[7]))!=entry[8]:raise ValueError('Correlated inner source packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_core):
    owner=CurrentOriginalN1CorrelatedInnerSource(accepted_core,require_checked=False)
    packets={name:owner.evaluate('.371',fraction) for name,fraction in (('actual_first_inlet','0'),('selected_actual_inner_exit','1'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(p) for name,p in packets.items()},input_hashes=owner.hashes,
             **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_CORRELATED_INNER_SOURCE_COMPUTED',len(packets),flush=True)
    return owner,packets

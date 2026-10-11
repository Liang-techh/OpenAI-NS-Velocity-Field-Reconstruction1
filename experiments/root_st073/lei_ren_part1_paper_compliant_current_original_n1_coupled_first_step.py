"""Actual first coupled n1 iterate, with a radius-aware remaining tail.

Wbase=(I+J)Gg; H=(I+J)G(B0-Akin+B1*d_Z).
Compute Wbase and H Wbase from source-owned whole-cell enclosures.
Regular normalized core cells include the axis without division by x.
Original hb, F0, Pstar^2 and Lambda remain formal source sectors.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_initial_integral as initial
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=initial.operator;domain=op.domain
HERE,PREFIX,sha=initial.HERE,initial.PREFIX,initial.sha
NAME=PREFIX+'current_original_n1_coupled_first_step.json.gz'
RECEIPT=PREFIX+'current_original_n1_coupled_first_step_check.json'
GATE='current_original_n1_first_coupled_regular_iteration_computed'
OPEN=initial.OPEN;FIELDS=initial.FIELDS;ep=initial.ep


def partition(c,j,N):return c.mpf([c.mpf(j).a/N,c.mpf(j+1).b/N])


def normalized_kernels(c,t,endpoint):
    """Exact substitution s=x*t also covers endpoint intervals at x=0."""
    f=endpoint*t*(1-t*t)/2
    return dict(F1=f,Uz1=endpoint*initial.log_kernel(c,t,c.mpf(1)),K1=-f,
                P1=c.mpf(1),d_x_F1=t**3,d_x_Uz1=t)


def weighted(algebra,rhs,weights,measure,collar=False,order=1):
    out={}
    for field,row in (('F1','5'),('Uz1','6'),('K1','6'),('P1','4'),('d_x_F1','5'),('d_x_Uz1','6')):
        value=rhs[row]*weights[field]*measure
        if collar:value=algebra.width(value)
        if any(key[0]<0 for key in value.terms):raise ValueError('Uncancelled inverse width in integrated n1 term')
        out[field]=op.truncate(value,order)
    return out


def feedback(algebra,system,base):
    """Off-kinematic rows, separately collecting B0 and genuine B1 d_Z."""
    zero=algebra.lift(0,0);out=[]
    for label in ('B0','B1'):
        rhs={str(row):zero for row in (4,5,6)}
        for row in (4,5,6):
            for col,field in enumerate(FIELDS,1):
                entry=system[label].get(f'{row},{col}')
                if entry is not None:
                    value=base[field] if label=='B0' else op.dz(base[field])
                    rhs[str(row)]=rhs[str(row)]+entry*value
        out.append(rhs)
    return tuple(out)


def perturbation_bound(c,b,epsilon,delta,eta,rho_max,Lmin):
    """Sum of explicit positive off-kinematic row bounds; no C0-1 trick."""
    upper=domain.upper;z=1+eta;d=1+z*z;x=c.sqrt(rho_max);rZ=eta/2
    value=lambda rows,i:rows[i].numeric_upper(c.mpf('1e-50000'))
    f,fr,w,wr,q=(value(b[key],i) for key,i in (('F',0),('F',1),('Uz',0),('Uz',1),('Q',0)))
    Cf=f+rho_max*fr
    Zf=((2+delta)*z*f+d*f/rZ+2*z*rho_max*fr)/Lmin
    Zw=((1+delta)*z*w+d*w/rZ+2*z*rho_max*wr)/Lmin
    H=(1-delta)*z/2+d*w;B=q+(1+2*z*w)/Lmin;E=z*w+c.mpf('.5')
    r4=4*epsilon*x*f
    r5=epsilon*(2*(q+(2-delta)*E/Lmin)+2*Zf+2*(1-delta)*z*Cf/Lmin+2*(1+delta)*z*Cf/Lmin+x*B)
    r6=epsilon*(8*epsilon*z*rho_max*f/Lmin+2*((1-delta)*E/Lmin+Zw+(1-delta)*z*rho_max*wr/Lmin)
        +2*(1+delta)*z*rho_max*wr/Lmin+4*z/Lmin+x*B)
    d5=epsilon*(2*H/Lmin+4*d*Cf/Lmin)
    d6=epsilon*(2*(H+d*rho_max*wr)/Lmin+2*d*rho_max*wr/Lmin+2*d/Lmin)
    return dict(off_kinematic_row_bounds=[upper(c,v) for v in (r4,r5,r6)],
                Cp=upper(c,r4+r5+r6),C1=upper(c,max(ep(d5)[1],ep(d6)[1])))


def preconditioned_tail(c,M0,a,Cp,C1,loss,N=1):
    """Nested Cauchy radii, k radial simplices, and all derivative words.

    For each k choose k gaps loss/k. Every step may differentiate once.
    Akin resummation adds optional radial integrations bounded by a.
    k! >= (k/e)^k yields a geometric majorant for every k>=N+1.
    """
    if isinstance(N,bool) or not isinstance(N,int) or N<0:raise ValueError('Nonnegative retained preconditioned index required')
    if any(ep(v)[0]<0 for v in (M0,a,Cp,C1)) or ep(loss)[0]<=0:raise ValueError('Positive analytic majorants and loss required')
    S=1+a;Mbase=domain.upper(c,S*M0)
    beta=domain.upper(c,c.e*a*S*(Cp/(N+1)+C1/loss))
    if ep(beta)[1]>=1:raise ValueError('Preconditioned geometric tail not admitted')
    return dict(retained_preconditioned_indices=[0,N],initial_base_common_tube_bound=Mbase,
        kinematic_J_norm_upper=a,I_plus_J_norm_upper=S,off_kinematic_B0_norm_upper=Cp,B1_norm_upper=C1,
        total_available_Z_loss=loss,equal_nested_radius_loss_for_length_k='loss/k',
        arbitrary_k_term_formula='Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!',
        factorial_to_geometric_formula='k! >= (k/e)^k',geometric_tail_ratio_upper=beta,
        omitted_sum_bound=domain.upper(c,Mbase*beta**(N+1)/(1-beta)),
        all_k_derivative_words_allowed=True,half_derivative_block_reduction_used=False,
        old_unpreconditioned_tail_reused=False)


@dataclass(frozen=True,eq=False)
class OriginalN1CoupledFirstStepPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1CoupledFirstStep:
    @source_precision
    def __init__(self,base,require_checked=True):
        if type(base) is not initial.CurrentOriginalN1InitialIntegral or not base.acceptance_loaded:
            raise ValueError('Accepted initial n1 integral owner required')
        base.assert_graph();self.base=base;self.inputs=base.inputs;self.source=base.source;self.c=base.c
        self.bounds=self.inputs.bounds;self.family=copy.deepcopy(base.family);self.hashes=dict(base.hashes)
        for name in (initial.NAME,initial.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._base=base;self._inputs=self.inputs;self._bounds=self.bounds;self._source=self.source;self._c=self.c;self._family=copy.deepcopy(self.family)
        self.ranges={};self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[key] for key in OPEN):
                raise ValueError('Coupled n1 receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Coupled n1 input changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.base is not self._base or self.inputs is not self._inputs or self.bounds is not self._bounds
            or self.bounds is not self.inputs.bounds or self.source is not self._source or self.c is not self._c
            or self.c is not self.base.c or self.source is not self.base.source or self.inputs is not self.base.inputs
            or not self.base.acceptance_loaded or self.family!=self._family or self.family!=self.base.family or self.hashes!=self._hashes):
            raise ValueError('Coupled n1 source graph changed')
        return self.base.assert_graph()

    def _system(self,Z,coordinate,chart):
        packet=self.inputs.evaluate(Z,coordinate,chart);view=self.inputs.report(packet)
        return self.inputs._packets[id(packet)][4],view

    @source_precision
    def _initial_range(self,Z,coordinate,chart,core_cells=16,collar_cells=4):
        self.assert_graph();c=self.c;z,algebra,_=self.inputs.frame(Z);coordinate=c.mpf(coordinate)
        key=(z._mpi_,coordinate._mpi_,chart,core_cells,collar_cells)
        if key in self.ranges:
            values,proof=self.ranges[key]
            if op.snapshot(values)!=proof:raise ValueError('Cached n1 initial range changed')
            return values
        if chart=='core':
            if ep(coordinate)[0]<0 or ep(coordinate)[1]>4:raise ValueError('Actual core interval in[0,4] required')
            X=c.sqrt(coordinate);core_end=X;phase=c.mpf(0)
        elif chart=='first':
            if ep(coordinate)[0]<0 or ep(coordinate)[1]>.5:raise ValueError('Actual first-collar interval in[0,1/2] required')
            _,end=self._system(z,coordinate,'first');X=c.sqrt(end['geometry']['rho_enclosure_only']);core_end=c.mpf(2);phase=coordinate
        else:raise ValueError('Actual core/first chart required')
        totals={name:algebra.lift(0,1) for name in FIELDS}
        if ep(core_end)[1]>0:
            for j in range(core_cells):
                t=partition(c,j,core_cells)
                rho=coordinate*t*t if chart=='core' else 4*t*t
                system,_=self._system(z,rho,'core')
                weights=normalized_kernels(c,t,X) if chart=='core' else initial.kernels(c,2*t,X)
                values=weighted(algebra,system['g'],weights,core_end/core_cells)
                totals={name:totals[name]+values[name] for name in FIELDS}
        if ep(phase)[1]>0:
            for j in range(collar_cells):
                s=phase*partition(c,j,collar_cells)
                system,view=self._system(z,s,'first');cell=c.sqrt(view['geometry']['rho_enclosure_only'])
                values=weighted(algebra,system['g'],initial.kernels(c,cell,X),cell*phase/(2*collar_cells),True)
                totals={name:totals[name]+values[name] for name in FIELDS}
        self.ranges[key]=(totals,op.snapshot(totals))
        return totals

    @source_precision
    def tail(self):
        self.assert_graph();c=self.c;packet=self.bounds.domain();view=self.bounds.report(packet)
        eta=view['common_tubes']['outer_Z_radius'];rho_max=view['actual_inner_geometry']['real_rho_upper_bound']
        restore=lambda rows:[domain.WidthBound(c,{int(p):item['coefficient'] for p,item in row.items()}) for row in rows]
        core={key:restore(value) for key,value in view['actual_core_radial_bounds'].items()}
        bridge={key:restore(value) for key,value in view['actual_collar_radial_bounds'].items()}
        args=(view['original_axis_amplitude']['complex_modulus_upper'],self.source.core.delta,eta,rho_max,
              view['denominator_bounds']['L_modulus_lower'])
        cb=domain.hierarchy_bounds(c,core['Phi'],core['Uz'],core['mean'],*args)
        bb=domain.hierarchy_bounds(c,bridge['Phi'],bridge['Uz'],bridge['mean'],*args)
        together={key:[a+b for a,b in zip(cb[key],bb[key])] for key in ('F','Uz','Q')}
        rows=perturbation_bound(c,together,self.source.core.epsilon,self.source.core.delta,eta,rho_max,args[-1])
        if ep(rows['C1'])[1]<ep(view['regular_B1_infinity_norm_upper'])[1]:raise ValueError('Inherited B1 norm not covered')
        tail=preconditioned_tail(c,view['initial_Gg_common_tube_bound'],c.sqrt(rho_max),rows['Cp'],rows['C1'],eta/8,1)
        return dict(**tail,off_kinematic_rows=rows['off_kinematic_row_bounds'],
                    initial_function_Z_radius=eta/4,solution_Z_radius=eta/8,matrix_Z_radius=eta/2,
                    scope='uniform analytic series tail for canonical core/actual-first-collar problem; later completed-leading supports unbound')

    @source_precision
    def evaluate(self,Z,coordinate,chart='core',outer_core_cells=8,outer_collar_cells=2):
        self.assert_graph();c=self.c;z,algebra,_=self.inputs.frame(Z);coordinate=c.mpf(coordinate)
        for count in (outer_core_cells,outer_collar_cells):
            if isinstance(count,bool) or not isinstance(count,int) or not 1<=count<=256:raise ValueError('Bounded positive outer partition required')
        if chart=='core':
            if ep(coordinate)[0]<0 or ep(coordinate)[1]>4:raise ValueError('Actual core endpoint required')
            X=c.sqrt(coordinate);core_end=X;phase=c.mpf(0)
        elif chart=='first':
            if ep(coordinate)[0]<0 or ep(coordinate)[1]>.5:raise ValueError('Actual first-collar endpoint required')
            _,end=self._system(z,coordinate,chart);X=c.sqrt(end['geometry']['rho_enclosure_only']);core_end=c.mpf(2);phase=coordinate
        else:raise ValueError('Only actual inner charts supported')
        began=time.monotonic();base=self._initial_range(z,coordinate,chart)
        corrections=[{name:algebra.lift(0,0) for name in FIELDS} for _ in range(2)];ledger=[]
        if ep(core_end)[1]>0:
            for j in range(outer_core_cells):
                t=partition(c,j,outer_core_cells);rho=coordinate*t*t if chart=='core' else 4*t*t
                cellbase=self._initial_range(z,rho,'core');system,_=self._system(z,rho,'core')
                rhs=feedback(algebra,system,cellbase)
                weights=normalized_kernels(c,t,X) if chart=='core' else initial.kernels(c,2*t,X)
                for index in range(2):
                    values=weighted(algebra,rhs[index],weights,core_end/outer_core_cells,order=0)
                    corrections[index]={name:corrections[index][name]+values[name] for name in FIELDS}
                ledger.append(dict(chart='core',rho_cell=rho,nested_initial_core_cells=16))
        if ep(phase)[1]>0:
            for j in range(outer_collar_cells):
                s=phase*partition(c,j,outer_collar_cells);cellbase=self._initial_range(z,s,'first')
                system,view=self._system(z,s,'first');cell=c.sqrt(view['geometry']['rho_enclosure_only'])
                rhs=feedback(algebra,system,cellbase)
                for index in range(2):
                    values=weighted(algebra,rhs[index],initial.kernels(c,cell,X),cell*phase/(2*outer_collar_cells),True,0)
                    corrections[index]={name:corrections[index][name]+values[name] for name in FIELDS}
                ledger.append(dict(chart='first',phase_cell=s,nested_initial_core_cells=16,nested_initial_collar_cells=4))
        total={name:op.truncate(base[name],0)+corrections[0][name]+corrections[1][name] for name in FIELDS}
        tail=self.tail()
        error=tail['omitted_sum_bound'];error_box=c.mpf([-ep(error)[1],ep(error)[1]])
        if chart=='core' and ep(coordinate)[1]==0:error_box=c.mpf(0)
        series={name:total[name]+algebra.lift(error_box,0) for name in FIELDS}
        result=dict(source_family=self.family,Z=z,endpoint=dict(chart=chart,coordinate=coordinate,x_enclosure_only=X),
            computed_preconditioned_indices=[0,1],computed_equation='W=Wbase+H W; H=(I+J)G(B0-Akin+B1*d_Z)',
            computed_prefix={name:op.export(value,0) for name,value in total.items()},
            canonical_inner_series_value_enclosures={name:op.export(value,0) for name,value in series.items()},
            canonical_inner_series_values_include_remaining_uniform_absolute_tail=True,
            B0_feedback={name:op.export(value,0) for name,value in corrections[0].items()},
            B1_derivative_feedback={name:op.export(value,0) for name,value in corrections[1].items()},
            initial_function={name:op.export(value,1) for name,value in base.items()},
            ordinary_Z_prefix_derivative_order=0,initial_Z_derivative_order=1,
            actual_whole_cell_initial_function_used=True,normalized_axis_cells_without_singular_division=True,
            pressure_2F0F1_feedback_included=True,B1_derivative_terms_computed_not_zeroed=True,
            actual_mean_and_independent_pressure_retained=True,inverse_width_cancelled_before_enclosure=True,
            source_log_basis_names=['hb','F0_anchor','Pstar_squared','Lambda'],source_log_bases=list(algebra.logs),
            analytic_preconditioned_tail=tail,integration_cells=ledger,
            error_budget=dict(computed_prefix='combined source and entire-cell integration enclosure',
                omitted_coupled_series='uniform absolute W norm bound in analytic_preconditioned_tail',
                useful_component_relative_accuracy_not_certified=True),
            completed_leading_support_binding_still_open=True,higher_prefix_derivatives_still_open=True,
            whole_matched_positive_order_field_not_delivered=True,execution_seconds=time.monotonic()-began,**dict.fromkeys(OPEN,False))
        packet=OriginalN1CoupledFirstStepPacket()
        self._packets[id(packet)]=(packet,result,total,corrections,op.snapshot(result),op.snapshot((total,corrections)))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1CoupledFirstStepPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live first-coupled n1 packet from this owner required')
        if op.snapshot(entry[1])!=entry[4] or op.snapshot((entry[2],entry[3]))!=entry[5]:raise ValueError('Coupled n1 packet changed')
        return dict(**copy.deepcopy(entry[1]),**{GATE:self.acceptance_loaded})


@source_precision
def run(base):
    began=time.monotonic();owner=CurrentOriginalN1CoupledFirstStep(base,require_checked=False)
    packets={name:owner.evaluate('.371',value,chart) for name,chart,value in (('core_exit','core','4'),('actual_first_keep','first','.5'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(value) for name,value in packets.items()},input_hashes=owner.hashes,
             **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_COUPLED_FIRST_ITERATION_COMPUTED',len(packets),flush=True)
    return owner,packets

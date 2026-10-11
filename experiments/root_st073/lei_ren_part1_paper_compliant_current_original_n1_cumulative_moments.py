"""Own temporal n1 common-inner cumulative moments, with radial cell ranges.

Two cumulants use exact own-mean/pressure FTC. The remaining three use
signed interval range integrals of the same original function. These
are partial common-inner defects, not a full-support repaired coefficient.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_second_radial_function as parent
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

regular=parent.parent;op=parent.op;natural=parent.natural;mean=parent.mean
core=regular.core;weighted=parent.weighted;domain=parent.domain;ep=parent.ep
HERE,PREFIX,sha=parent.HERE,parent.PREFIX,parent.sha
NAME=PREFIX+'current_original_n1_cumulative_moments.json.gz'
RECEIPT=PREFIX+'current_original_n1_cumulative_moments_check.json'
GATE='current_original_n1_common_inner_five_cumulative_moment_enclosures_installed'
OPEN=parent.OPEN
MOMENTS=('m11','m12','m13','m14','m15')


def core_swirl_range(owner,algebra,Z,rho):
    """The admitted analytic W0 and genuine errors, uniformly on a rho cell."""
    c=owner.c;axis=owner.physical.axis_inputs(Z);chi=axis['chi'].truncate(1)
    gradient=op.Jet(c,owner.regular.weighted.source.amplitude.gradient_jets(Z))
    m=(gradient*axis['d']/axis['L']).truncate(1);alpha=m*m
    model,tail=core.model_chi_derivatives(c,chi,rho,48)
    value=algebra.lift(model['F1']*alpha);error=algebra.lift(tail['F1']*core.norm_jet(alpha))
    full=c.mpf([0,ep(rho)[1]]);raw=owner.physical.normalized_jets(full,Z)
    key=domain.parent.parent.core_module.gridkey
    phi=[algebra.lift(op.Jet(c,[raw['Phi'][key(i,k)]/math.factorial(k) for k in range(4-i)])) for i in range(3)]
    source=core.normalized_theta_sectors(algebra,phi,algebra.lift(op.Jet.variable(c,Z,3)),
        full,owner.physical.delta,algebra.lift(gradient.truncate(3)))
    correction=algebra.lift(op.Jet(c,[raw['nonlinear_correction_bounds'][key(0,k)]/math.factorial(k)
                                    for k in range(2)]))*core.norm_jet(alpha)*2
    norm=core.inverse_norms(c,c.sqrt(c.mpf(ep(rho)[1])),domain.upper(c,abs(chi[1])))['F1']
    for residual in (correction,weighted.absolute_jet(source['Lambda_minus1']),
                     weighted.absolute_jet(source['Lambda_minus2'])):
        error+=op.truncate(residual*norm,1)
    return value,error,gradient


def collar_swirl_range(owner,algebra,Z,phase,rho,cells,dependencies):
    """Variable-endpoint version of the accepted core-plus-collar integral.

    For every S in phase, substitute s=S*t on a fixed t partition. Source,
    kernel and measure boxes therefore contain the whole cumulative path.
    """
    c=owner.c;cp=owner.regular.core.evaluate(Z,'4');owner.regular.core.report(cp)
    dependencies.append((owner.regular.core,cp));entry=owner.regular.core._packets[id(cp)]
    values,errors=entry[3:5]
    if values['F1'].algebra is not algebra:raise ValueError('Same original source frame required')
    axis=owner.physical.axis_inputs(Z);chi=axis['chi'].truncate(1);chiZ=domain.upper(c,abs(chi[1]))
    gradient=op.Jet(c,owner.regular.weighted.source.amplitude.gradient_jets(Z));X=c.sqrt(rho)
    correlated=natural.mean.axial.inner
    transfer,proof=correlated.homogeneous_transfer_error(c,correlated.vector_norm(values,errors),chiZ,X)
    value=values['F1'];error=errors['F1']+transfer
    tail=weighted.potential_kernel_tail(c,owner.regular.weighted.a,chiZ,12)['F1']
    if ep(phase)[1]>0:
        for j in range(cells):
            s=phase*weighted.budget.coupled.partition(c,j,cells)
            packet=owner.inputs.evaluate(Z,s,'first');view=owner.inputs.report(packet)
            dependencies.append((owner.inputs,packet));pe=owner.inputs._packets[id(packet)]
            rad=view['geometry']['rho_enclosure_only'];x=c.sqrt(rad)
            source=core.normalized_theta_sectors(algebra,pe[2]['Phi'],algebra.lift(op.Jet.variable(c,Z,4)),
                rad,owner.physical.delta,algebra.lift(gradient))
            g=sum(source.values(),algebra.lift(0,1));ratio=x/X
            t=c.mpf([max(0,ep(ratio)[0]),min(1,ep(ratio)[1])])
            kernel=weighted.potential_kernel(chi,t,X,12,1)[0]
            measure=x*phase/(2*cells)
            actual=algebra.width(g*kernel*measure)
            omitted=algebra.width(weighted.absolute_jet(g)*tail*measure)
            if any(p[0]<0 for item in (actual,omitted) for p in item.terms):
                raise ValueError('Actual collar source inverse hb did not cancel')
            value+=op.truncate(actual,1);error+=op.truncate(omitted,1)
    return value,error,gradient,proof


def physical_values(algebra,W,gradient):
    return {name:weighted.transform(value,gradient,scale,1)
            for name,value,scale in zip(regular.FIELDS[:4],W,natural.WEIGHTS[:4])}


def rho_integrands(algebra,leadingF,leadingU,physical,rho):
    """dR=epsilon drho; m15 uses exact integration by parts."""
    F,U,P=(physical[k] for k in ('F1','Uz1','P1'))
    return dict(m12=algebra.shift(F*(2*rho),(0,0,0,-2)),
                m14=algebra.shift((leadingU*F+U*leadingF)*(2*rho),(0,0,0,-2)),
                m15_parts=algebra.shift(leadingU*U*2+P,(0,0,0,-1)))


def endpoint_cumulants(algebra,W,rho):
    return dict(m11=algebra.shift((W[1]+W[2])*rho,(0,0,0,-2)),
                m13=algebra.shift(W[3],(0,0,1,-1)),
                m15_boundary=algebra.shift(W[3]*(-rho),(0,0,1,-2)))


def core_model_integrals(owner,algebra,Z,rho,gradient):
    """Integrate the correlated signed core model before enclosing remainders."""
    c=owner.c;axis=owner.physical.axis_inputs(Z);chi=axis['chi'].truncate(1)
    alpha=(gradient*axis['d']/axis['L']).truncate(1)**2
    I=op.Jet.constant(c,0,1);power=op.Jet.constant(c,1,1)
    for n in range(1,49):
        coefficient=c.mpf((-1)**n)/(2**n*math.factorial(n-1)*math.factorial(n+1)*(n+2))
        I+=power*(coefficient*rho**(n+2));power*=chi
    physicalI=weighted.transform(algebra.lift(I*alpha),gradient,natural.WEIGHTS[0],1)
    rows=mean.core_rows(algebra,owner.physical,c.mpf(0),Z)
    z=algebra.lift(axis['z']);u=algebra.lift(axis['u']);d=1-z*z;L=1-z*z*owner.physical.delta
    zero=algebra.lift(0,5);C=op.truncate(-2*op.axial_second([u,zero,zero],z,c.mpf(0),owner.physical.delta,-1-owner.physical.delta),1)
    q=op.truncate((rows['slope']*(z*(3+owner.physical.delta))-op.dz(rows['slope'])*d)/(2*L),1)
    return dict(m12=algebra.shift(physicalI*2,(0,0,0,-2)),
                m14=algebra.shift(physicalI*op.truncate(u,1)*2,(0,0,0,-2)),
                m15_parts=algebra.shift((op.truncate(u,1)*C/4+q)*(rho*rho),(0,0,0,-2)))


def core_residual_integrands(owner,algebra,Z,rho,rows,W,physical,errors,gradient,leadingF):
    """Remainders in separated source atoms; no boxed-model subtraction."""
    c=owner.c;axis=rows['axis'];z=algebra.lift(axis['z']);u=op.truncate(algebra.lift(axis['u']),1)
    d=1-z*z;L=1-z*z*owner.physical.delta
    errorF=weighted.transform(regular.symmetric_error(errors[0]),gradient,natural.WEIGHTS[0],1)
    errorU=weighted.transform(regular.symmetric_error(errors[1]),gradient,natural.WEIGHTS[1],1)
    raw=owner.physical.normalized_jets(rho,Z);key=domain.parent.parent.core_module.gridkey
    nonlinear=[algebra.lift(op.Jet(c,[c.mpf([-1,1])*raw['nonlinear_correction_bounds'][key(i,k)]/math.factorial(k)
                                      for k in range(6-i)])) for i in (0,1)]
    qerror=[(P*(2*z)-(P/2)*(z*(1-owner.physical.delta))-op.dz(P/2)*d)/L for P in nonlinear]
    pressure_error=op.truncate(algebra.shift(qerror[0]+qerror[1]*rho,(0,0,0,-1)),1)
    pressure_error+=weighted.transform(regular.symmetric_error(errors[3]),gradient,natural.WEIGHTS[3],1)
    psi=op.truncate(rows['Psi'][0],1);F,U=physical['F1'],physical['Uz1']
    return dict(m12=algebra.shift(errorF*(2*rho),(0,0,0,-2)),
                m14=algebra.shift((u*errorF+U*leadingF)*(2*rho),(0,0,0,-2))
                    +algebra.shift(psi*F*(2*rho),(0,0,0,-3)),
                m15_parts=algebra.shift(u*errorU*2+pressure_error,(0,0,0,-1))
                          +algebra.shift(psi*U*2,(0,0,0,-2)))


def moment_budget(value,k):
    """Normalize an already-enclosed physical cumulant, preserving its width."""
    result=core.integral_budget(value,value.algebra.lift(0,1),k)
    interval=result.get('initial_source_integral_normalized_enclosure')
    return dict(reference_source_exponents=result.get('reference_source_exponents'),
        signed_physical_partial_moment_normalized_enclosure=interval,
        total_width_over_abs_lower=result.get('total_width_over_abs_lower'),
        relative_width_target='.001',
        relative_partial_moment_width_certified=result.get('relative_initial_integral_width_certified',False),
        exact_zero=result.get('initial_integral_exact_zero',False),
        sign=('negative' if interval is not None and ep(interval)[1]<0 else
              'positive' if interval is not None and ep(interval)[0]>0 else 'zero_or_unresolved'),
        full_function_and_feedback_enclosure_used=True,terminal_full_support_closure_not_implied=True)


@dataclass(frozen=True,eq=False)
class OriginalN1CumulativeMomentPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1CumulativeMoments:
    @source_precision
    def __init__(self,accepted_second,require_checked=True):
        if type(accepted_second) is not parent.CurrentOriginalN1SecondRadialFunction or not accepted_second.acceptance_loaded:
            raise ValueError('Accepted original second-radial owner required')
        accepted_second.assert_graph();self.second=accepted_second;self.regular=accepted_second.regular
        self.inputs=self.regular.mean.inputs;self.c=accepted_second.c;self.physical=accepted_second.physical
        self.family=copy.deepcopy(accepted_second.family);self.hashes=dict(accepted_second.hashes)
        for name in (parent.NAME,parent.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._parents=(self.second,self.regular,self.inputs);self._identities=(self.c,self.physical)
        self._packets={};self._cells={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Own cumulative moment receipt scope differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Own cumulative moment input changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if ((self.second,self.regular,self.inputs)!=self._parents or (self.c,self.physical)!=self._identities or
            self.regular is not self.second.regular or self.inputs is not self.regular.mean.inputs or
            not self.second.acceptance_loaded or self.family!=self.second.family or self.hashes!=self._hashes):
            raise ValueError('Own cumulative moment source graph changed')
        self.second.assert_graph()

    @source_precision
    def range_cell(self,Z,coordinate,chart,source_cells=2):
        """Uniform value/Z1 function range, without any 1/x axis division."""
        self.assert_graph();c=self.c;interval=c.mpf(coordinate);z=c.mpf(Z)
        if isinstance(source_cells,bool) or not isinstance(source_cells,int) or not 1<=source_cells<=64:
            raise ValueError('Source partition in[1,64] required')
        key=(z._mpi_,interval._mpi_,chart,source_cells)
        if key in self._cells:
            data,frozen=self._cells[key]
            if op.snapshot(data[:-1])!=frozen:raise ValueError('Radial cell range changed')
            for source,packet in data[-1]:source.report(packet)
            return data
        packet=self.inputs.evaluate(z,interval,chart);view=self.inputs.report(packet);entry=self.inputs._packets[id(packet)]
        algebra=entry[5];rho=view['geometry'].get('rho',view['geometry'].get('rho_enclosure_only'))
        dependencies=[(self.inputs,packet)];uniform,proof=mean.complex_mean_bounds(algebra,self.physical,
            self.regular.mean.bound,self.regular.weighted.tube)
        rows=(mean.core_rows(algebra,self.physical,rho,z) if chart=='core' else
              mean.actual_rows(algebra,self.physical,entry[2],rho,interval,z,uniform,proof['original_outer_eta']))
        fields=mean.factored_leading_rows(algebra,rows,self.physical)
        if chart=='core':fvalue,ferror,gradient=core_swirl_range(self,algebra,z,rho);transfer=None
        elif chart=='first':fvalue,ferror,gradient,transfer=collar_swirl_range(self,algebra,z,interval,rho,source_cells,dependencies)
        else:raise ValueError('Original core/first source range required')
        seeds=natural.seeds(algebra,rows,fields,z,rho,self.physical.delta);data=self.regular.natural._data
        centers=[fvalue,*[op.truncate(seeds[k],1) for k in ('Uz1','K1','P1')]]
        gap=self.regular.tail['derivative_Cauchy_gap'];shape=domain.upper(c,rho)
        feedback=algebra.lift(op.Jet(c,[domain.upper(c,shape*self.regular.tail['value_remainder_upper']/gap**k) for k in (0,1)]))
        errors=[ferror+feedback]
        for k in ('Uz1','K1','P1'):
            bound=data['initial_errors']['P1' if k=='P1' else 'axial']
            residual=algebra.lift(op.Jet(c,[domain.upper(c,bound/gap**j) for j in (0,1)]))
            errors.append(residual+feedback)
        W=[v+regular.symmetric_error(e) for v,e in zip(centers,errors)]
        physical=physical_values(algebra,W,gradient)
        leadingF=algebra.shift(rows['Phi'][0]*weighted.relative_amplitude(algebra,gradient,1,1),(0,1,0,0))
        integrands=rho_integrands(algebra,op.truncate(leadingF,1),op.truncate(fields['Uz'][0],1),physical,rho)
        residuals=(core_residual_integrands(self,algebra,z,rho,rows,W,physical,errors,gradient,op.truncate(leadingF,1))
                   if chart=='core' else integrands)
        result=(rho,rows,fields,W,physical,integrands,centers,errors,transfer,residuals,dependencies)
        self._cells[key]=(result,op.snapshot(result[:-1]));return result

    @source_precision
    def evaluate(self,Z,coordinate,chart='core',core_cells=4,collar_cells=2,source_cells=2):
        self.assert_graph();began=time.monotonic();c=self.c
        for count in (core_cells,collar_cells,source_cells):
            if isinstance(count,bool) or not isinstance(count,int) or not 1<=count<=64:
                raise ValueError('Radial and source partitions in[1,64] required')
        p=self.regular.evaluate(Z,coordinate,chart);pv=self.regular.report(p);pe=self.regular._packets[id(p)]
        W,_,_,_,_,system=pe[5];algebra=W[0].algebra;z=pv['Z'];rho=system['rho']
        endpoint=endpoint_cumulants(algebra,[op.truncate(v,1) for v in W],rho)
        integrals={k:algebra.lift(0,1) for k in ('m12','m14','m15_parts')};records=[];cells=[]
        upper=c.mpf(coordinate) if chart=='core' else c.mpf(4)
        if ep(upper)[1]>0:
            gradient=system['gradient'];integrals=core_model_integrals(self,algebra,z,upper,gradient)
            for j in range(core_cells):
                interval=upper*weighted.budget.coupled.partition(c,j,core_cells)
                cell=self.range_cell(z,interval,'core',source_cells);cells.append(cell)
                for k in integrals:integrals[k]+=cell[9][k]*(upper/core_cells)
                records.append(dict(chart='core',coordinate=interval,rho=cell[0],
                    measure='drho',integrands={k:op.export(v,1) for k,v in cell[5].items()},
                    separated_residual_integrands={k:op.export(v,1) for k,v in cell[9].items()},
                    normalized_function_range={k:op.export(v,1) for k,v in zip(regular.FIELDS[:4],cell[3])}))
        if chart=='first':
            phase=pv['query']['actual_source_coordinate']
            if ep(phase)[1]>0:
                for j in range(collar_cells):
                    interval=phase*weighted.budget.coupled.partition(c,j,collar_cells)
                    cell=self.range_cell(z,interval,'first',source_cells);cells.append(cell)
                    for k in integrals:integrals[k]+=algebra.width(cell[5][k]*(cell[0]*phase/collar_cells))
                    records.append(dict(chart='first',coordinate=interval,rho=cell[0],
                        measure='drho=hb*rho*ds',actual_fixed_phase_independent_of_Z=True,
                        integrands={k:op.export(v,1) for k,v in cell[5].items()},
                        normalized_function_range={k:op.export(v,1) for k,v in zip(regular.FIELDS[:4],cell[3])}))
        cumulative={k:endpoint[k] for k in ('m11','m13')}
        cumulative.update(m12=integrals['m12'],m14=integrals['m14'],
                          m15=integrals['m15_parts']+endpoint['m15_boundary'])
        budgets={k:{'Z'+str(j):moment_budget(v,j) for j in (0,1)} for k,v in cumulative.items()}
        result=dict(source_family=self.family,Z=z,query=pv['query'],geometry=pv['geometry'],
            physical_partial_cumulative_moments={k:op.export(cumulative[k],1) for k in MOMENTS},
            partial_moment_value_Z1_budgets=budgets,
            passing_partial_moment_value_Z1_rows=sum(v['relative_partial_moment_width_certified'] for rows in budgets.values() for v in rows.values()),
            integrated_parts={k:op.export(v,1) for k,v in integrals.items()},
            exact_endpoint_terms={k:op.export(v,1) for k,v in endpoint.items()},radial_cells=records,
            ordinary_Z_cumulative_derivative_order=1,radial_cell_range_enclosures_not_point_samples=True,
            source_range_variable_endpoint_substitution='s=S*t; t partition fixed for all S in cell',
            original_core_history_retained_at_actual_collar=True,actual_collar_Jacobian_and_hb_retained=True,
            axis_cell_avoids_singular_equation_divisions=True,leading_times_n1_cross_terms_used=True,
            full_pressure_sources_and_remaining_feedback_included=True,
            m11_own_mean_FTC_and_m13_full_pressure_FTC=True,m15_full_pressure_integration_by_parts=True,
            correlated_core_model_integrated_analytically=True,
            core_model_factorial_nonlinear_pressure_inertia_feedback_errors_retained=True,
            no_independent_boxed_model_subtraction_in_core_remainders=True,
            local_partial_cumulative_moments_not_terminal_full_support_defects=True,
            whole_domain_relative_moment_accuracy_not_certified=True,full_support_repair_and_n2_still_open=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1CumulativeMomentPacket();state=(endpoint,integrals,cumulative,records)
        self._packets[id(token)]=(token,result,p,state,cells,op.snapshot((result,state)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1CumulativeMomentPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live own n1 cumulative moment packet required')
        self.regular.report(entry[2])
        if op.snapshot((entry[1],entry[3]))!=entry[5]:raise ValueError('Own cumulative moment packet changed')
        for cell in entry[4]:
            if not any(value[0] is cell and value[1]==op.snapshot(cell[:-1]) for value in self._cells.values()):
                raise ValueError('Own cumulative moment radial range changed')
            for source,token in cell[-1]:source.report(token)
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_second):
    owner=CurrentOriginalN1CumulativeMoments(accepted_second,require_checked=False)
    packets={name:owner.evaluate(z,r,chart) for name,z,r,chart in (
        ('core_axis','.371','0','core'),('core_interior','.371','2','core'),('core_exit','.371','4','core'),
        ('actual_first_inlet','.371','0','first'),('selected_actual_inner_exit','.371','1','first'))}
    raw=dict(source_family=owner.family,samples={k:owner.report(v) for k,v in packets.items()},
             input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_OWN_CUMULATIVE_MOMENTS_COMPUTED',len(packets),flush=True)
    return owner,packets

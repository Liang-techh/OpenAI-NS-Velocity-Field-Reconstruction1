"""Correlated original core swirl source integral, with genuine errors.

The Bessel model solves the dominant forced equation analytically.  The
original nonlinear correction and both lower-Lambda forcing sectors are
bounded separately.  This is an initial source integral, not a coupled
n1 solution; the actual prescribed-shear collar is outside this adapter.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator as weighted
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=weighted.op;domain=weighted.domain;budget=weighted.budget;ep=weighted.ep
HERE,PREFIX,sha=weighted.HERE,weighted.PREFIX,weighted.sha
NAME=PREFIX+'current_original_n1_analytic_core_source.json.gz'
RECEIPT=PREFIX+'current_original_n1_analytic_core_source_check.json'
GATE='current_original_n1_analytic_dominant_core_source_integral_enclosed'
OPEN=weighted.OPEN
FIELDS=('F1','d_x_F1')
ZERO=(0,0,0,0)


def plain_dz(value):
    if value.order<1:raise ValueError('Genuine ordinary source derivative required')
    return op.Jet(value.ctx,[(k+1)*value[k+1] for k in range(value.order)])


def norm_jet(value,order=1):
    if value.order<order:raise ValueError('Missing source derivative for norm')
    return op.Jet(value.ctx,[domain.upper(value.ctx,abs(v)) for v in value.coefficients[:order+1]])


def normalized_theta_sectors(algebra,phi_rows,z,rho,delta,gradient):
    """Exact g5/(F0*Lambda), with F0 derivative cancellations exposed first."""
    f,fr,_=phi_rows;d=1-z*z;L=1-z*z*delta;m=gradient*d/L
    a=-2-delta;b=-3
    cross=m*(f*(z*(a+b))+op.dz(f)*(2*d)-fr*(z*(4*rho)))/L+op.dz(m)*d*f/L
    second=op.axial_second(phi_rows,z,rho,delta,a)
    terms={'dominant':op.truncate(-2*m*m*f,1),
           'Lambda_minus1':algebra.shift(op.truncate(2*cross,1),(0,0,0,-1)),
           'Lambda_minus2':algebra.shift(op.truncate(-2*second,1),(0,0,0,-2))}
    return terms


def model_chi_derivatives(c,chi,rho,N=48):
    """Partial_chi B(chi*rho/2) and its x derivative, plus infinite tails."""
    if isinstance(N,bool) or not isinstance(N,int) or not 4<=N<=128:
        raise ValueError('Bounded factorial-model cutoff in[4,128] required')
    if ep(rho)[0]<0 or ep(rho)[1]>4 or ep(chi[0])[0]<0 or ep(chi[0])[1]>1:
        raise ValueError('Original real core radius and chi in[0,1] required')
    chi=chi.truncate(1);x=c.sqrt(rho)
    f=op.Jet.constant(c,0,1);dx=op.Jet.constant(c,0,1);power=op.Jet.constant(c,1,1)
    for n in range(1,N+1):
        coefficient=c.mpf((-1)**n)/(2**n*math.factorial(n-1)*math.factorial(n+1))
        f+=power*(coefficient*rho**n)
        dx+=power*(coefficient*(2*n)*x*rho**(n-1))
        power*=chi
    if ep(rho)[1]==0:
        return {field:op.Jet.constant(c,0,1) for field in FIELDS},{field:op.Jet.constant(c,0,1) for field in FIELDS}
    n=N+1;r=domain.upper(c,rho);u=domain.upper(c,x);cz=domain.upper(c,abs(chi[1]))
    first=r**n/(2**n*math.factorial(n-1)*math.factorial(n+1))
    first_dx=2*n*u*r**(n-1)/(2**n*math.factorial(n-1)*math.factorial(n+1))
    q=domain.upper(c,r/(2*n*(n+2)))
    ratios=(q,domain.upper(c,q*n/(n-1)),domain.upper(c,q*(n+1)/n),domain.upper(c,q*(n+1)/(n-1)))
    if any(ep(v)[1]>=1 for v in ratios):raise ValueError('Model derivative tail ratio not admitted')
    tails={'F1':op.Jet(c,[domain.upper(c,first/(1-ratios[0])),domain.upper(c,first*(n-1)*cz/(1-ratios[1]))]),
           'd_x_F1':op.Jet(c,[domain.upper(c,first_dx/(1-ratios[2])),domain.upper(c,first_dx*(n-1)*cz/(1-ratios[3]))])}
    return dict(zip(FIELDS,(f,dx))),tails


def inverse_norms(c,x,chi_Z,N=12):
    """Absolute integral/derivative norms of the full Lchi inverse, |chi|<=1."""
    weighted.kernel_polynomials(N)
    out={field:op.Jet.constant(c,0,1) for field in FIELDS}
    for j in range(N+1):
        base=2**j*x**(2*j+2)/(4**(j+1)*math.factorial(j+1)*math.factorial(j+2))
        dxbase=2**j*(2*j+2)*x**(2*j+1)/(4**(j+1)*math.factorial(j+1)*math.factorial(j+2))
        out['F1']+=op.Jet(c,[base,j*chi_Z*base])
        out['d_x_F1']+=op.Jet(c,[dxbase,j*chi_Z*dxbase])
    if ep(x)[1]>0:
        tail=weighted.potential_kernel_tail(c,domain.upper(c,x),chi_Z,N)
        for field in FIELDS:out[field]+=tail[field]*x
    return {field:norm_jet(v) for field,v in out.items()}


def integral_budget(value,error,k):
    """Signed initial source integral enclosure; no remaining-feedback claim."""
    c=value.ctx;prefix=budget.component_budget(value,k,c.mpf(0))
    if prefix['prefix_exact_zero']:
        return dict(initial_integral_exact_zero=not error.terms,relative_initial_integral_width_certified=False,
                    remaining_B0_B1_feedback_included=False)
    ref=prefix['reference_source_exponents'];M=prefix['directed_reference_coefficient_norm']
    total_error=c.mpf(0);sectors=[];unresolved=False
    for powers,jet in sorted(error.terms.items()):
        mag=domain.upper(c,abs(jet[k]))
        if ep(mag)[1]==0:continue
        log=budget.log_ratio(c,value.algebra.logs,powers,ref,mag,M)
        upper,rule=budget.exp_cover(c,log)
        sectors.append(dict(source_exponents=list(powers),normalized_log_upper=log,norm_upper=upper,bound_rule=rule))
        if upper is None:unresolved=True
        else:total_error+=upper
    enclosure=None;width=None;passed=False
    if not unresolved and prefix['signed_prefix_normalized_enclosure'] is not None:
        enclosure=prefix['signed_prefix_normalized_enclosure']+c.mpf([-1,1])*domain.upper(c,total_error)
        if ep(enclosure)[0]>0 or ep(enclosure)[1]<0:
            width=domain.upper(c,(enclosure.b-enclosure.a)/domain.lower(c,abs(enclosure)))
            passed=ep(width)[1]<=ep(c.mpf('.001'))[0]
    return dict(reference_source_exponents=ref,directed_reference_coefficient_norm=M,
        signed_model_integral_normalized_enclosure=prefix['signed_prefix_normalized_enclosure'],
        initial_source_integral_normalized_enclosure=enclosure,total_width_over_abs_lower=width,
        absolute_error_normalized_reference_norm_upper=None if unresolved else domain.upper(c,total_error),
        error_source_sector_contributions=sectors,relative_initial_integral_width_target='.001',
        relative_initial_integral_width_certified=passed,initial_integral_exact_zero=False,
        source_bound_normalization_not_selected_value=True,remaining_B0_B1_feedback_included=False)


@dataclass(frozen=True,eq=False)
class OriginalN1AnalyticCoreSourcePacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1AnalyticCoreSource:
    @source_precision
    def __init__(self,accepted_weighted,require_checked=True):
        if type(accepted_weighted) is not weighted.CurrentOriginalN1WeightedPropagator or not accepted_weighted.acceptance_loaded:
            raise ValueError('Accepted original weighted propagator owner required')
        accepted_weighted.assert_graph();self.weighted=accepted_weighted;self.inputs=accepted_weighted.inputs
        self.physical=accepted_weighted.physical;self.c=accepted_weighted.c
        self.family=copy.deepcopy(accepted_weighted.family);self.hashes=dict(accepted_weighted.hashes)
        for name in (weighted.NAME,weighted.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._weighted=accepted_weighted;self._inputs=self.inputs;self._physical=self.physical;self._c=self.c
        self._parameters=op.snapshot((self.physical.correction,self.physical.h,self.physical.delta))
        self.acceptance_loaded=False;self._packets={}
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Analytic core source receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Analytic core source dependency changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.weighted is not self._weighted or self.inputs is not self._inputs or self.physical is not self._physical
            or self.c is not self._c or self.inputs is not self.weighted.inputs or self.physical is not self.weighted.physical
            or self.c is not self.weighted.c or not self.weighted.acceptance_loaded or self.family!=self.weighted.family
            or self.hashes!=self._hashes or self._parameters!=op.snapshot((self.physical.correction,self.physical.h,self.physical.delta))):
            raise ValueError('Analytic core source graph changed')
        return self.weighted.assert_graph()

    @source_precision
    def evaluate(self,Z,rho,model_degree=48):
        self.assert_graph();c=self.c;r=c.mpf(rho)
        if ep(r)[0]<0 or ep(r)[1]>4:raise ValueError('Original core endpoint rho in[0,4] required; actual collar is separate')
        if not isinstance(rho,(str,int,float)) and ep(r)[0]!=ep(r)[1]:raise ValueError('One exact core outer endpoint required')
        began=time.monotonic();p=self.inputs.evaluate(Z,r,'core');view=self.inputs.report(p)
        z=view['Z'];algebra=self.inputs._packets[id(p)][5]
        axis=self.physical.axis_inputs(z);chi=axis['chi'].truncate(1);gradient=op.Jet(c,self.weighted.source.amplitude.gradient_jets(z))
        d=1-axis['z']**2;L=axis['L'];m=(gradient*d/L).truncate(1);alpha=m*m
        model,model_tail=model_chi_derivatives(c,chi,r,model_degree)
        values={field:algebra.lift(model[field]*alpha) for field in FIELDS}
        errors={field:algebra.lift(model_tail[field]*norm_jet(alpha)) for field in FIELDS}
        pieces={field:{'model_factorial_tail':errors[field]} for field in FIELDS}
        uniform=self.physical.normalized_jets(c.mpf([0,ep(r)[1]]),z)
        key=domain.parent.parent.core_module.gridkey
        rows=[algebra.lift(op.Jet(c,[uniform['Phi'][key(i,k)]/math.factorial(k) for k in range(4-i)])) for i in range(3)]
        zjet=algebra.lift(op.Jet.variable(c,z,3));gjet=algebra.lift(gradient.truncate(3))
        forcing=normalized_theta_sectors(algebra,rows,zjet,c.mpf([0,ep(r)[1]]),self.physical.delta,gjet)
        nonlinear=algebra.lift(op.Jet(c,[uniform['nonlinear_correction_bounds'][key(0,k)]/math.factorial(k) for k in range(2)]))*norm_jet(alpha)*2
        norms=inverse_norms(c,c.sqrt(r),domain.upper(c,abs(chi[1])))
        for field in FIELDS:
            for label,source in (('genuine_nonlinear_core_correction',nonlinear),
                                 ('Lambda_minus1_forcing',weighted.absolute_jet(forcing['Lambda_minus1'])),
                                 ('Lambda_minus2_forcing',weighted.absolute_jet(forcing['Lambda_minus2']))):
                contribution=op.truncate(source*norms[field],1)
                pieces[field][label]=contribution;errors[field]+=contribution
        physical={field:weighted.transform(v,gradient,(0,1,0,1),1) for field,v in values.items()}
        physical_errors={field:algebra.shift(weighted.absolute_jet(weighted.relative_amplitude(algebra,gradient,1,1))*v,(0,1,0,1))
                         for field,v in errors.items()}
        diagnostics={field:{'Z'+str(k):integral_budget(values[field],errors[field],k) for k in (0,1)} for field in FIELDS}
        physical_diagnostics={field:{'Z'+str(k):integral_budget(physical[field],physical_errors[field],k) for k in (0,1)} for field in FIELDS}
        result=dict(source_family=self.family,Z=z,rho=r,chart='original_core_only',
            normalized_source_identity='gstar5=-2*m^2*Phi+2*Lambda^-1*(Tb(m*Phi)+m*Ta(Phi))-2*Lambda^-2*Tb(Ta(Phi))',
            axial_weights={'a':'-2-delta','b':'-3'},m=list(m.coefficients),alpha=list(alpha.coefficients),
            model_identity='Lchi Phi_model=0; Lchi(partial_chi Phi_model)=-2 Phi_model',
            analytic_dominant_integral='m^2*partial_chi B(chi*rho/2)',model_degree=model_degree,
            full_inverse_absolute_norms={field:list(v.coefficients) for field,v in norms.items()},
            raw_normalized_core_forcing_sectors={label:op.export(v,1) for label,v in forcing.items()},
            normalized_analytic_model_integrals={field:op.export(v,1) for field,v in values.items()},
            normalized_error_components={field:{label:op.export(v,1) for label,v in row.items()} for field,row in pieces.items()},
            normalized_total_error_bounds={field:op.export(v,1) for field,v in errors.items()},
            physical_analytic_model_integrals={field:op.export(v,1) for field,v in physical.items()},
            physical_total_error_bounds={field:op.export(v,1) for field,v in physical_errors.items()},
            normalized_initial_integral_budgets=diagnostics,physical_initial_integral_budgets=physical_diagnostics,
            passing_normalized_initial_width_rows=sum(row['relative_initial_integral_width_certified'] for f in diagnostics.values() for row in f.values()),
            initial_source_integral_ordinary_Z_derivative_order=1,
            original_nonlinear_field_enclosed_not_replaced_by_model=True,
            dominant_source_integrated_analytically_without_radial_cell_boxing=True,
            lower_Lambda_sectors_extracted_before_boxing=True,F0_and_Lambda_not_materialized=True,
            full_potential_inverse_included_in_core_error_bounds=True,
            actual_prescribed_shear_collar_not_model_extension=True,actual_first_collar_integral_still_open=True,
            core_initial_source_integral_only=True,remaining_normalized_B0_B1_feedback_not_solved=True,
            old_coupled_prefix_not_replaced_by_initial_integral=True,whole_Z_or_full_inner_solution_not_promoted=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalN1AnalyticCoreSourcePacket()
        self._packets[id(packet)]=(packet,result,p,values,errors,op.snapshot((result,values,errors)))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1AnalyticCoreSourcePacket or entry is None or entry[0] is not packet:
            raise ValueError('Live original analytic core source packet required')
        self.inputs.report(entry[2])
        if op.snapshot((entry[1],entry[3],entry[4]))!=entry[5]:raise ValueError('Analytic core source packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_weighted):
    owner=CurrentOriginalN1AnalyticCoreSource(accepted_weighted,require_checked=False)
    packets={name:owner.evaluate('.371',rho) for name,rho in (('core_axis','0'),('core_interior','2'),('core_exit','4'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(p) for name,p in packets.items()},input_hashes=owner.hashes,
             **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_ANALYTIC_CORE_SOURCE_COMPUTED',len(packets),flush=True)
    return owner,packets

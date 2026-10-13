"""Actual leading pressure integrals from Rin to the modification inlet.

Use the original core continuation as a reference with proved collar errors.
Flat sigma and its first derivative keep two powers of the original hb;
they are magnitude enclosures of the actual functions, not chosen values.
Direct integration of N1p retains hb in dR and avoids endpoint cancellation.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_leading_mean_epsilon_factor as mean
import lei_ren_part1_paper_compliant_current_original_n1_cumulative_moments as moments

op, domain, ep = mean.op, mean.domain, mean.ep
HERE, PREFIX, sha, source_precision = mean.HERE, mean.PREFIX, mean.sha, mean.source_precision
NAME = PREFIX+'current_original_n1_unchanged_pressure_prefix.json.gz'
RECEIPT = PREFIX+'current_original_n1_unchanged_pressure_prefix_check.json'
GATE = 'current_original_n1_unchanged_Rin_to_rminus_pressure_integrals_installed'
OPEN = tuple(dict.fromkeys(mean.OPEN+('complete_signed_pressure_tails_certified',
    'positive_order_five_moment_repair_closed','temporal_recursion',
    'full_corrected_Navier_Stokes_solution')))


def collar_derivative_corrections(corepsi, corephi, chi, chir, ratio, dpsi0, dmean0, rho):
    """Exact difference equations; coefficients remain source-correlated."""
    A=corephi[1]/corephi[0]
    dpsi1=(ratio*chi-1)*corepsi[1]
    dpsi2=ratio*chir*corepsi[1]+ratio*chi*(chi-1)*A*corepsi[1]+(ratio*chi-1)*corepsi[2]
    dmean1=(dpsi0-dmean0)/rho
    dmean2=(dpsi1-2*dmean1)/rho
    return dict(dpsi1=dpsi1,dpsi2=dpsi2,dmean1=dmean1,dmean2=dmean2)


def physical_pressure_range_integrals(algebra, forcing, rho, lower_phase, upper_phase):
    """Integrate a genuine whole-cell range with its positive physical measure."""
    length=upper_phase-lower_phase
    return dict(m13=algebra.width(algebra.shift(forcing*(rho*length),(0,0,0,-1))),
                m15=-algebra.width(algebra.shift(forcing*(rho*rho*length),(0,0,0,-2))))


def factored_sigma_bounds(algebra, phase_upper):
    """sigma^(j)<=Cj exp(4-s^-2)s^-3j <= hb²*1e-60, j=0,1."""
    c=algebra.ctx;upper=c.mpf(ep(phase_upper)[1]);logh=algebra.logs[0]
    if not 0<ep(upper)[0] or ep(upper)[1]>ep(c.mpf('.5'))[0]:
        raise ValueError('Original positive left flat-tail interval required')
    cap=c.mpf('1e-60');proofs=[];values=[]
    for order,constant in ((0,1),(1,4)):
        logarithm=c.ln(constant)+4-1/(upper*upper)-3*order*c.ln(upper)
        relative=logarithm-2*logh
        if ep(relative)[1]>=ep(c.ln(cap))[0]:
            raise ValueError('Original sigma tail does not admit two hb powers')
        values.append(algebra.width(algebra.lift(c.mpf([0,ep(cap)[1]]),4),2))
        proofs.append(dict(ordinary_phase_derivative_order=order,
            original_flat_tail_log_upper=logarithm,log_sigma_over_hb_squared_upper=relative,
            nonzero_coefficient_upper=cap,source_hb_power=2,
            source_value_not_materialized_or_replaced=True,
            tail_bound='Cj*exp(4-s^-2)*s^(-3j); C0=1,C1=4,0<s<=1/2'))
    return values,proofs


def complex_core_slope_bounds(algebra, physical, bound):
    c=algebra.ctx;eta=bound['common_tubes']['outer_Z_radius'];z=1+eta;d=1+z*z
    L=bound['denominator_bounds']['L_modulus_lower'];dt=abs(physical.delta)
    u=4*z+abs(physical.j);H=z*(1+dt)/2+d*u
    p=bound['original_pressure']['normalized_modulus_upper'];pz=2*p/eta
    hydro=domain.upper(c,((1+2*z*u)*u*(1+dt)/2+4*H)/(2*L))
    pressure=domain.upper(c,(d*pz+2*z*p*(1+dt))/(2*L))
    slope=algebra.lift(hydro,0)+algebra.shift(algebra.lift(pressure,0),(0,0,1,0))
    correction=domain.upper(c,4*physical.correction)
    return slope+algebra.lift(2*correction,0)


def refined_original_rows(owner, algebra, Z, rho, phase_upper):
    """Reference core plus actual ODE/own-mean errors on the original collar."""
    c=algebra.ctx;physical=owner.physical;b=owner.mean.bound
    reference=mean.core_rows(algebra,physical,rho,Z)
    uniform,proof=mean.complex_mean_bounds(algebra,physical,b,owner.mean.weighted.tube)
    sigma,sigma_proof=factored_sigma_bounds(algebra,phase_upper)
    sig0,sig1=sigma
    hb=algebra.width(algebra.lift(1,4))
    chi=algebra.lift(1,4)+sig0*(hb-1)
    chir=op.truncate(sig1/rho-algebra.width(sig1/rho,-1),3)
    rmax=domain.upper(c,rho);smax=domain.upper(c,phase_upper)
    A=proof['actual_core_logarithmic_derivative_modulus_upper']
    ratio_coefficient=domain.upper(c,4*rmax*smax*A*c.mpf('1e-60'))
    ratio_log=3*algebra.logs[0]+c.ln(ratio_coefficient)
    if ep(ratio_log)[1]>=ep(c.ln(c.mpf('.5')))[0]:
        raise ValueError('Actual complex ratio needs a smaller Gronwall interval')
    ratio_bound=algebra.width(algebra.lift(ratio_coefficient,0),3)
    gap=proof['point_increment_Taylor_Cauchy_radius']
    ratio_error=mean.cauchy_signed_rows(ratio_bound,4,gap)
    ratio=algebra.lift(1,4)+ratio_error
    C=complex_core_slope_bounds(algebra,physical,b)
    # |chi-1|<=(1+hb)*sigma<=2*hb²*cap. Integrate over
    # drho=hb*rho*ds; the original ratio error contributes hb^4.
    psi_bound=algebra.width(C*(rmax*smax)*(algebra.width(algebra.lift(c.mpf('2e-60'),0),2)+ratio_bound))
    mean_bound=algebra.width(psi_bound*(rmax*smax/4))
    dpsi0=mean.cauchy_signed_rows(psi_bound,5,gap)
    dmean0=mean.cauchy_signed_rows(mean_bound,5,gap)
    corepsi,coremean=reference['Psi'],reference['meanPsi']
    differences=collar_derivative_corrections(corepsi,reference['Phi'],chi,chir,ratio,dpsi0,dmean0,rho)
    dpsi1=op.truncate(differences['dpsi1'],4);dpsi2=op.truncate(differences['dpsi2'],3)
    dmean1=op.truncate(differences['dmean1'],4);dmean2=op.truncate(differences['dmean2'],3)
    A=op.truncate(reference['Phi'][1]/reference['Phi'][0],3)
    phi=[ratio*reference['Phi'][0],ratio*chi*reference['Phi'][1],
         op.truncate(ratio*(chir*reference['Phi'][1]+chi*(chi-1)*A*reference['Phi'][1]+chi*reference['Phi'][2]),3)]
    rows={**reference,
        'Phi':phi,
        'Psi':[corepsi[0]+dpsi0,corepsi[1]+dpsi1,corepsi[2]+dpsi2],
        'meanPsi':[coremean[0]+dmean0,coremean[1]+dmean1,coremean[2]+dmean2]}
    fields=mean.factored_leading_rows(algebra,rows,physical)
    identities=dict(original_sigma_bounds=sigma_proof,
        exact_original_chi='1+(hb-1)*sigma(s)',
        exact_original_ratio_equation='d_rho log(Phi_actual/Phi_core)=(chi-1)*Phi_core_rho/Phi_core',
        exact_original_Psi_equation='Psi_actual_rho=(Phi_actual/Phi_core)*chi*Psi_core_rho',
        exact_original_own_mean_equation='rho*meanPsi_actual=4*meanPsi_core(4)+integral_4^rho Psi_actual drho',
        reference_core_is_not_selected_as_actual_field=True,
        ratio_Gronwall_log_upper=ratio_log,ratio_error_modulus_per_hb_cubed=ratio_coefficient,
        ratio_error_bound=op.export(ratio_bound,0),Psi_error_bound=op.export(psi_bound,0),
        own_mean_error_bound=op.export(mean_bound,0),Cauchy_Taylor_gap=gap,
        actual_core_complex_slope_bound=op.export(C,0),common_Z_tube_preserved=True,
        actual_comparison_errors_retained=True,phase_and_radial_map_Z_independent=True)
    return rows,fields,identities


@dataclass(frozen=True,eq=False)
class OriginalUnchangedPressurePrefixPacket:
    leading_hierarchy_order:int=0


class CurrentOriginalN1UnchangedPressurePrefix:
    @source_precision
    def __init__(self,accepted_mean,require_checked=True):
        if type(accepted_mean) is not mean.CurrentOriginalN1LeadingMeanEpsilonFactor or not accepted_mean.acceptance_loaded:
            raise ValueError('Accepted original leading mean field required')
        accepted_mean.assert_graph();self.mean=accepted_mean;self.inputs=accepted_mean.inputs
        self.physical=accepted_mean.physical;self.local=accepted_mean.local;self.c=accepted_mean.c
        self.family=copy.deepcopy(accepted_mean.family);self.hashes=dict(accepted_mean.hashes)
        for name in (mean.NAME,mean.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.phase_in=self.local.inner.phase_in;self.modification_inlet=2*self.phase_in
        if self.local.inner.sc._mpi_!=(4*self.phase_in)._mpi_:
            raise ValueError('Original common Rin=sc/4 and rminus=sc/2 required')
        self._identities=(self.mean,self.inputs,self.physical,self.local,self.c)
        self._phase=op.snapshot((self.phase_in,self.modification_inlet))
        self._packets={};self.acceptance_loaded=False
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked['all_passed'] or not checked[GATE] or checked['source_family']!=self.family or any(checked[key] for key in OPEN):
                raise ValueError('Unchanged pressure prefix receipt/scope differs')
            for name,digest in checked['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Unchanged pressure prefix input changed '+name)
                self.hashes[name]=digest
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if ((self.mean,self.inputs,self.physical,self.local,self.c)!=self._identities or self.hashes!=self._hashes
            or not self.mean.acceptance_loaded or self.family!=self.mean.family
            or op.snapshot((self.phase_in,self.modification_inlet))!=self._phase):
            raise ValueError('Unchanged original pressure prefix graph changed')
        self.mean.assert_graph()

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();began=time.monotonic();c=self.c
        z,algebra,_=self.inputs.frame(Z);a,b=self.phase_in,self.modification_inlet
        selector=self.local.selector
        # The class module owns the canonical selector, not the live packet.
        from lei_ren_part1_paper_compliant_current_original_n1_local_leading_dispatch import original
        support=original.OriginalMicroFiniteN.modification_support(selector,'first_micro',0,b)
        if not support['increments_exact_zero']:
            raise ValueError('The endpoint must precede actual finite-N modification')
        loghbphase=algebra.logs[0]+c.ln(b)
        cap=c.mpf('1e-80')
        if ep(loghbphase)[1]>=ep(c.ln(cap/2))[0]:
            raise ValueError('Actual rho enclosure needs a sharper radial source bound')
        rho=c.mpf([4,ep(4*(1+cap))[1]])
        if ep(rho)[1]>ep(c.mpf('4.1'))[0]:raise ValueError('Original core continuation bound required')
        if ep(rho)[1]>ep(self.mean.bound['actual_inner_geometry']['real_rho_upper_bound'])[1]:
            raise ValueError('Refined rho must remain inside the accepted core bound')
        rows,fields,proof=refined_original_rows(self,algebra,z,rho,b)
        gradient=op.Jet(c,self.mean.weighted.source.amplitude.gradient_jets(z))
        F=[algebra.shift(phi*mean.weighted.relative_amplitude(algebra,gradient,1,n),(0,1,0,0))
           for phi,n in zip(rows['Phi'],mean.ORDERS)]
        originals=dict(F=F,Uz=fields['Uz'],Q=fields['Q'],mean=fields['mean'])
        forcing=op.truncate(op.known_forcing(algebra,originals,algebra.lift(op.Jet.variable(c,z,5)),rho,self.physical.delta)['N1p'],1)
        values=physical_pressure_range_integrals(algebra,forcing,rho,a,b)
        budgets={};scaled_values={}
        for name,value in values.items():
            budgets[name]={};scaled_values[name]={}
            for k in (0,1):
                row=moments.moment_budget(value,k)
                row.pop('full_function_and_feedback_enclosure_used',None)
                row['actual_leading_source_and_collar_error_enclosure_used']=True
                budgets[name]['Z'+str(k)]=row
                scaled_values[name]['Z'+str(k)]=moments.core.integral_budget(value,algebra.lift(0,1),k)
        result=dict(source_family=self.family,Z=z,query=dict(chart='original_unchanged_first_collar',
            lower_phase=a,upper_phase=b,lower_is_actual_Rin=True,upper_is_actual_rminus=True),
            geometry=dict(actual_rho_recipe='4*exp(hb*s)',rho_enclosure_only=rho,
                actual_hb_phase_log_upper=loghbphase,positive_radial_allowance=cap,
                nonzero_radius_allowance_not_source_value=True,
                physical_measure='dR=Lambda^-1*hb*rho*ds',weighted_physical_measure='R*dR=Lambda^-2*hb*rho^2*ds'),
            original_modification_support=support,actual_source_reference_and_error_proof=proof,
            original_leading_Uz_rho_rows=[op.export(value) for value in fields['Uz']],
            original_leading_Q_rho_rows=[op.export(value) for value in fields['Q']],
            actual_known_pressure_forcing_Z1=op.export(forcing,1),
            physical_unchanged_pressure_integrals_Z1={name:op.export(value,1) for name,value in values.items()},
            signed_source_scaled_pressure_budgets=budgets,
            signed_source_scaled_pressure_values=scaled_values,
            direct_same_source_pressure_integral_not_endpoint_subtraction=True,
            original_hb_measure_retained_once=True,original_m15_minus_R_N1p_sign_retained=True,
            absolute_radius_hb_and_N_not_materialized=True,actual_own_mean_and_pressure_atoms_retained=True,
            core_comparison_has_nonzero_source_bound_errors=True,
            complete_pressure_tails_and_n1_repair_not_claimed=True,whole_Z_precision_not_certified=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalUnchangedPressurePrefixPacket();state=(rows,fields,forcing,values)
        self._packets[id(packet)]=(packet,result,state,op.snapshot((result,state)))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalUnchangedPressurePrefixPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live original unchanged pressure packet required')
        if op.snapshot((entry[1],entry[2]))!=entry[3]:raise ValueError('Unchanged pressure packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_mean):
    owner=CurrentOriginalN1UnchangedPressurePrefix(accepted_mean,False)
    packets={name:owner.evaluate(Z) for name,Z in (('actual_unchanged_tail','.371'),
                 ('actual_unchanged_tail_independent_Z','.211'),('actual_unchanged_tail_second_Z','.639'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(p) for name,p in packets.items()},
             input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_UNCHANGED_PRESSURE_PREFIX_COMPUTED',len(packets),flush=True)
    return owner,packets

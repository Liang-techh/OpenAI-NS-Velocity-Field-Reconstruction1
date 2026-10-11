"""Original axial correction with epsilon and pressure factors exposed.

Uses the core model plus admitted nonlinear correction and the exact actual
first-collar drive identity. It bounds the residual swirl gauge on the
original complex tube; it does not solve remaining n1 feedback.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_correlated_inner_source as inner
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

core=inner.core;weighted=inner.weighted;op=inner.op;domain=inner.domain;ep=inner.ep
HERE,PREFIX,sha=inner.HERE,inner.PREFIX,inner.sha
NAME=PREFIX+'current_original_n1_axial_epsilon_factor.json.gz'
RECEIPT=PREFIX+'current_original_n1_axial_epsilon_factor_check.json'
GATE='current_original_n1_axial_epsilon_factor_and_residual_gauge_bounds_installed'
OPEN=inner.OPEN


def core_psi_rows(algebra,physical,rho,z,order=2):
    """Genuine Psi=slope*rho+nonlinear correction; expose Pstar^2 in slope."""
    c=algebra.ctx;axis=physical.axis_inputs(z)
    zj=axis['z'];u=axis['u'];d=axis['d'];L=axis['L'];H=axis['H'];p=axis['pressure']
    hydro=((1-2*zj*u)*u*((1+physical.delta)/2)+4*H)/(2*L)
    pressure=(d*core.plain_dz(p)-2*zj*p*(1+physical.delta))/(2*L)
    slope=algebra.lift(hydro.truncate(order))+algebra.shift(algebra.lift(pressure.truncate(order)),(0,0,1,0))
    raw=physical.normalized_jets(rho,z);key=domain.parent.parent.core_module.gridkey
    rows=[]
    for i in (0,1):
        errors=op.Jet(c,[c.mpf([-1,1])*raw['nonlinear_correction_bounds'][key(i,k)]/math.factorial(k) for k in range(order+1)])
        rows.append(op.truncate(slope*(rho if i==0 else 1)+algebra.lift(errors),order))
    if ep(c.mpf(rho))[1]==0:rows[0]=algebra.lift(0,order)
    return rows


def complex_psi_bounds(algebra,physical,bound,chi_tube):
    """On eta/2, pressure Cauchy derivative comes from the original eta tube."""
    c=algebra.ctx;eta=bound['common_tubes']['outer_Z_radius'];z=1+eta/2
    d=1+z*z;dt=abs(physical.delta);L=bound['denominator_bounds']['L_modulus_lower']
    u=4*z+abs(physical.j);H=z*(1+dt)/2+d*u
    p=bound['original_pressure']['normalized_modulus_upper'];pz=2*p/eta
    hydro=domain.upper(c,((1+2*z*u)*u*(1+dt)/2+4*H)/(2*L))
    pressure=domain.upper(c,(d*pz+2*z*p*(1+dt))/(2*L))
    slope=algebra.lift(hydro,0)+algebra.shift(algebra.lift(pressure,0),(0,0,1,0))
    radius=bound['actual_inner_geometry']['real_rho_upper_bound']
    correction=domain.upper(c,4*physical.correction)
    initial=slope*radius+algebra.lift(correction,0)
    radial=slope+algebra.lift(2*correction,0)  # rho Cauchy radius 1/2
    floor=bound['denominator_bounds']['Phi_core_modulus_lower']
    phi=bound['actual_collar_radial_bounds']['Phi'][0]['0']['coefficient']
    ratio=domain.upper(c,phi/floor);cap=bound['actual_inner_geometry']['numerical_cap_bounds_only']
    # |chiSpatial|<=1+cap, d_rho=hb*rho*ds, phase length<=1/2.
    actual=initial+algebra.width(radial*(ratio*(1+cap)*radius/2))
    b=chi_tube['H_imaginary_modulus_upper'];S=domain.lower(c,physical.sigma**2-b*b)
    if ep(S)[0]<=0:raise ValueError('Original complex H denominator gap required')
    K=domain.upper(c,1/(2*c.sqrt(S))+b/S)
    gauge=algebra.shift(actual*(2*d*K),(0,0,0,-1))
    gap=eta/4
    derivative=op.FactoredJet(algebra,{p:op.Jet(c,[row[0],domain.upper(c,row[0]/gap)]) for p,row in gauge.terms.items()},1)
    return dict(Psi_actual=actual,residual_gauge=derivative),dict(
        Psi_value_Z_tube_radius=eta/2,gauge_value_Z_tube_radius=eta/2,
        gauge_Z1_Z_tube_radius=eta/4,ordinary_Z1_Cauchy_gap=gap,
        original_outer_eta_preserved=eta,pressure_normalized_modulus_upper=p,
        pressure_Z_modulus_upper_on_eta_over2=pz,actual_Phi_over_core_modulus_upper=ratio,
        H_over_H_squared_plus_sigma_squared_modulus_upper=K,H_denominator_gap_lower=S,
        proof_K='|u+iv|/(u^2+sigma^2-v^2)<=1/(2sqrt(S))+b/S; S=sigma^2-b^2',
        core_nonlinear_correction_joint_value_factor=4,core_nonlinear_rho_derivative_Cauchy_factor=2,
        source_pressure_scale_not_materialized=True,complex_bounds_are_majorants_not_selected_values=True)


def small_sector_bound(value,k):
    c=value.ctx;total=c.mpf(0);terms=[];unresolved=False
    for powers,row in sorted(value.terms.items()):
        mag=domain.upper(c,abs(row[k]))
        if ep(mag)[1]==0:continue
        log=core.budget.log_ratio(c,value.algebra.logs,powers,(0,0,0,0),mag,c.mpf(1))
        cover,rule=core.budget.exp_cover(c,log)
        terms.append(dict(source_exponents=list(powers),log_modulus_upper=log,modulus_cover=cover,bound_rule=rule))
        if cover is None:unresolved=True
        else:total+=cover
    return dict(modulus_upper=None if unresolved else domain.upper(c,total),source_sector_bounds=terms,
                bound_is_absolute_not_relative_to_an_unknown_field=True)


@dataclass(frozen=True,eq=False)
class OriginalN1AxialEpsilonPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1AxialEpsilonFactor:
    @source_precision
    def __init__(self,accepted_inner,require_checked=True):
        if type(accepted_inner) is not inner.CurrentOriginalN1CorrelatedInnerSource or not accepted_inner.acceptance_loaded:
            raise ValueError('Accepted original correlated initial source required')
        accepted_inner.assert_graph();self.inner=accepted_inner;self.weighted=accepted_inner.weighted
        self.inputs=accepted_inner.inputs;self.local=accepted_inner.local;self.physical=self.weighted.physical;self.c=accepted_inner.c
        self.family=copy.deepcopy(accepted_inner.family);self.hashes=dict(accepted_inner.hashes)
        for name in (inner.NAME,inner.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.domain_packet=self.inputs.bounds.domain();self.bound=self.inputs.bounds.report(self.domain_packet)
        self._inner=accepted_inner;self._inputs=self.inputs;self._local=self.local;self._physical=self.physical;self._c=self.c
        self._bound=op.snapshot(self.bound);self.acceptance_loaded=False;self._packets={}
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Axial epsilon receipt scope/family differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Axial epsilon source changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.inner is not self._inner or self.inputs is not self._inputs or self.local is not self._local
            or self.physical is not self._physical or self.c is not self._c or self.inputs is not self.inner.inputs
            or self.local is not self.inner.local or self.physical is not self.inner.weighted.physical
            or self.weighted is not self.inner.weighted or not self.inner.acceptance_loaded
            or self.family!=self.inner.family or self.hashes!=self._hashes or op.snapshot(self.bound)!=self._bound):
            raise ValueError('Axial epsilon graph changed')
        self.inputs.bounds.report(self.domain_packet);return self.inner.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate,chart='first',cells=4):
        self.assert_graph();c=self.c;value=c.mpf(coordinate)
        if not isinstance(coordinate,(str,int,float)) and ep(value)[0]!=ep(value)[1]:raise ValueError('One exact endpoint required')
        if isinstance(cells,bool) or not isinstance(cells,int) or not 1<=cells<=128:raise ValueError('Bounded actual collar partition required')
        if chart not in ('core','first'):raise ValueError('Original core or selected first collar required')
        began=time.monotonic();lp=self.local.evaluate(Z,coordinate,chart);lv=self.local.report(lp)
        z=lv['Z'];q=lv['query'];actual=q['actual_source_coordinate'];p=self.inputs.evaluate(z,actual,chart);pv=self.inputs.report(p)
        algebra=self.inputs._packets[id(p)][5];rho=pv['geometry'].get('rho',pv['geometry'].get('rho_enclosure_only'))
        axis=self.physical.axis_inputs(z);d=algebra.lift(axis['d'].truncate(1));L=algebra.lift(axis['L'].truncate(1))
        gradient=algebra.lift(op.Jet(c,self.weighted.source.amplitude.gradient_jets(z)).truncate(1))
        inlet=core_psi_rows(algebra,self.physical,(rho if chart=='core' else c.mpf(4)),z,2)[0]
        psi=op.truncate(inlet,1);increments=algebra.lift(0,1);ledger=[];packets=[]
        if chart=='first' and ep(actual)[1]>0:
            for j in range(cells):
                phase=actual*weighted.budget.coupled.partition(c,j,cells)
                cell=self.inputs.evaluate(z,phase,'first');sv=self.inputs.report(cell);entry=self.inputs._packets[id(cell)]
                r=sv['geometry']['rho_enclosure_only'];raw=self.physical.normalized_jets(r,z)
                phi_core=op.Jet(c,[raw['Phi']['rho0_Z'+str(k)]/math.factorial(k) for k in range(3)])
                ratio=op.truncate(entry[2]['Phi'][0]/phi_core,1)
                sigma=domain.sigma_jets(c,phase)[0]
                chiSpatial=algebra.lift(1-sigma,1)+algebra.width(algebra.lift(sigma,1))
                psir=op.truncate(core_psi_rows(algebra,self.physical,r,z,2)[1],1)
                integral=algebra.width(ratio*chiSpatial*psir*(r*actual/cells))
                increments+=integral;packets.append(cell)
                ledger.append(dict(phase_cell=phase,rho_enclosure=r,integrand=op.export(ratio*chiSpatial*psir,1),
                                   exact_rho_measure='hb*rho*ds',source_is_actual_Phi_over_core_Phi=True))
            psi+=increments
        delta_Uz=algebra.shift(psi,(0,0,0,-1));delta_H=delta_Uz*d
        gauge=op.truncate(-2*delta_H*gradient/L,1)
        uniform,proof=complex_psi_bounds(algebra,self.physical,self.bound,self.weighted.tube)
        result=dict(source_family=self.family,Z=z,query=q,geometry=pv['geometry'],
            exact_core_identity='Uz_core-uaxis=epsilon*Psi_core',
            exact_collar_identity='Uz_actual-uaxis=epsilon*[Psi_core(4)+integral_4^rho chiSpatial*(Phi_actual/Phi_core)*Psi_core_rho d_rho]',
            exact_gauge_residual='-2*(Hactual-Haxis)*Gprime/L=-2*epsilon*d*Psi_actual*Haxis/(Haxis^2+sigma^2)',
            factored_Psi=op.export(psi,1),factored_Psi_inlet=op.export(inlet,1),factored_actual_collar_Psi_increment=op.export(increments,1),
            factored_Uz_minus_uaxis=op.export(delta_Uz,1),factored_Hactual_minus_Haxis=op.export(delta_H,1),
            factored_remaining_swirl_gauge=op.export(gauge,1),actual_collar_Psi_source_cells=ledger,
            uniform_complex_Psi_actual_bound=op.export(uniform['Psi_actual'],0),
            uniform_complex_residual_gauge_bound=op.export(uniform['residual_gauge'],1),uniform_complex_bound_proof=proof,
            uniform_gauge_absolute_budgets={'Z'+str(k):small_sector_bound(uniform['residual_gauge'],k) for k in (0,1)},
            original_Pstar_squared_and_epsilon_sectors_retained=True,original_complex_tube_not_epsilon_shrunk=True,
            actual_delta_V_not_used_as_total_Uz_minus_axis=True,original_actual_field_not_replaced_by_comparison=True,
            exact_source_width_not_materialized=True,remaining_normalized_B0_B1_feedback_not_solved=True,
            full_weighted_residual_matrix_majorant_still_open=True,whole_Z_or_full_inner_solution_not_promoted=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1AxialEpsilonPacket();state=(result,psi,delta_Uz,delta_H,gauge,uniform)
        self._packets[id(token)]=(token,result,lp,p,packets,state,op.snapshot(state))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1AxialEpsilonPacket or entry is None or entry[0] is not packet:raise ValueError('Live original axial epsilon packet required')
        self.local.report(entry[2]);self.inputs.report(entry[3])
        for cell in entry[4]:self.inputs.report(cell)
        if op.snapshot(entry[5])!=entry[6]:raise ValueError('Axial epsilon packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_inner):
    owner=CurrentOriginalN1AxialEpsilonFactor(accepted_inner,require_checked=False)
    packets={name:owner.evaluate('.371',coordinate,chart) for name,coordinate,chart in
             (('core_axis','0','core'),('core_exit','4','core'),('actual_first_inlet','0','first'),('selected_actual_inner_exit','1','first'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(p) for name,p in packets.items()},input_hashes=owner.hashes,
             **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_AXIAL_EPSILON_FACTOR_COMPUTED',len(packets),flush=True)
    return owner,packets

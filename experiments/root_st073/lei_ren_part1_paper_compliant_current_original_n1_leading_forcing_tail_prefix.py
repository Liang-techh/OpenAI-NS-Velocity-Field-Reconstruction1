"""Original fixed leading-pressure forcing from Rin to r_minus.

Integrate only S and R*S; retain exact viscous endpoint identities. This
extends two own n1 pressure-moment histories, not their complete support.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_radial_cutoff_extension as parent
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=parent.op;mean=parent.parent.mean;weighted=parent.weighted;domain=parent.domain
ep=parent.ep;HERE,PREFIX,sha=parent.HERE,parent.PREFIX,parent.sha
NAME=PREFIX+'current_original_n1_leading_forcing_tail_prefix.json.gz'
RECEIPT=PREFIX+'current_original_n1_leading_forcing_tail_prefix_check.json'
GATE='current_original_n1_fixed_leading_forcing_tail_prefix_to_rminus_installed'
OPEN=parent.OPEN


def forcing_parts(algebra,fields,Z,rho,delta):
    """Physical S, and finite endpoint terms with common Qaxis removed."""
    c=algebra.ctx;z=algebra.lift(op.Jet.variable(c,Z,4));d=1-z*z;L=1-delta*z*z
    q,qr=fields['Q'][:2];u=fields['Uz'][0]
    time_term=(q+op.dz(q)*(z*((1-delta)/2))+qr*rho)/L
    radial=q*(q/2+qr*rho)
    transport=u*(-q*(2*z)+op.dz(q)*d-qr*(2*z*rho))/L
    S=op.truncate(time_term+radial+transport,1)
    qp,qpr=fields['Qpsi'][:2]
    # Qaxis is independent of R, so it cancels exactly in [Q+R Q_R].
    B0=algebra.shift(op.truncate(qp+qpr*rho,1),(0,0,0,-1))
    B1=algebra.shift(op.truncate(qpr*rho*rho,1),(0,0,0,-2))
    return dict(S=S,viscous_endpoint_without_common_Qaxis=B0,
        weighted_viscous_endpoint=B1,
        time=op.truncate(time_term,1),radial=op.truncate(radial,1),transport=op.truncate(transport,1))


def reduce_integrals(S_integral,RS_integral,lower,upper):
    return dict(N1p_integral=-S_integral/2+upper['viscous_endpoint_without_common_Qaxis']-lower['viscous_endpoint_without_common_Qaxis'],
        R_N1p_integral=-RS_integral/2+upper['weighted_viscous_endpoint']-lower['weighted_viscous_endpoint'])


@dataclass(frozen=True,eq=False)
class OriginalN1LeadingForcingTailPrefixPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1LeadingForcingTailPrefix:
    @source_precision
    def __init__(self,accepted_cutoff,require_checked=True):
        if type(accepted_cutoff) is not parent.CurrentOriginalN1RadialCutoffExtension or not accepted_cutoff.acceptance_loaded:
            raise ValueError('Accepted original radial-cutoff owner required')
        accepted_cutoff.assert_graph();self.cutoff=accepted_cutoff;self.inputs=accepted_cutoff.inputs
        self.mean=accepted_cutoff.regular.mean;self.physical=accepted_cutoff.physical;self.c=accepted_cutoff.c
        self.family=copy.deepcopy(accepted_cutoff.family);self.hashes=dict(accepted_cutoff.hashes)
        for name in (parent.NAME,parent.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._identities=(self.cutoff,self.inputs,self.mean,self.physical,self.c)
        self._cells={};self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Leading forcing tail-prefix receipt or scope differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Leading tail-prefix source changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if ((self.cutoff,self.inputs,self.mean,self.physical,self.c)!=self._identities or
            self.inputs is not self.cutoff.inputs or self.mean is not self.cutoff.regular.mean or
            not self.cutoff.acceptance_loaded or self.family!=self.cutoff.family or self.hashes!=self._hashes):
            raise ValueError('Leading forcing tail-prefix graph changed')
        self.cutoff.assert_graph()

    @source_precision
    def source_cell(self,Z,phase):
        self.assert_graph();c=self.c;z=c.mpf(Z);s=c.mpf(phase)
        if ep(s)[0]<0 or ep(s)[1]>ep(c.mpf('.5'))[0]:raise ValueError('Original first phase in[0,1/2] required')
        inlet,_=self.cutoff.selected_phase(Z)
        if ep(s)[1]>ep(2*inlet)[0]:
            raise ValueError('Modified finite-N leading Q/mean source required beyond r_minus=2*Sin')
        key=(z._mpi_,s._mpi_)
        if key in self._cells:
            entry=self._cells[key]
            if op.snapshot(entry[:6])!=entry[7]:raise ValueError('Leading forcing cell changed')
            self.inputs.report(entry[6]);return entry
        packet=self.inputs.evaluate(z,s,'first');view=self.inputs.report(packet);pe=self.inputs._packets[id(packet)]
        algebra=pe[5];rho=view['geometry']['rho_enclosure_only']
        uniform,proof=mean.complex_mean_bounds(algebra,self.physical,self.mean.bound,self.mean.weighted.tube)
        rows=mean.actual_rows(algebra,self.physical,pe[2],rho,s,z,uniform,proof['original_outer_eta'])
        fields=mean.factored_leading_rows(algebra,rows,self.physical)
        parts=forcing_parts(algebra,fields,z,rho,self.physical.delta)
        gradient=op.Jet(c,self.mean.weighted.source.amplitude.gradient_jets(z))
        F=[algebra.shift(phi*weighted.relative_amplitude(algebra,gradient,1,n),(0,1,0,0))
           for phi,n in zip(rows['Phi'],mean.ORDERS)]
        originals=dict(F=F,Uz=fields['Uz'],Q=fields['Q'],mean=fields['mean'])
        parts['N1p_direct']=op.truncate(op.known_forcing(algebra,originals,
            algebra.lift(op.Jet.variable(c,z,5)),rho,self.physical.delta)['N1p'],1)
        entry=(z,s,rho,rows,fields,parts,packet,None);entry=(*entry[:7],op.snapshot(entry[:6]))
        self._cells[key]=entry;return entry

    @source_precision
    def evaluate(self,Z,upper_phase=None,cells=4):
        self.assert_graph();began=time.monotonic();c=self.c
        a,local=self.cutoff.selected_phase(Z);b=2*a if upper_phase is None else c.mpf(upper_phase)
        if isinstance(cells,bool) or not isinstance(cells,int) or not 1<=cells<=64:
            raise ValueError('Leading tail partition in[1,64] required')
        if upper_phase is not None and not isinstance(upper_phase,(str,int,float)) and ep(b)[0]!=ep(b)[1]:raise ValueError('One exact tail endpoint required')
        lower=self.source_cell(Z,a);upper=self.source_cell(Z,b)
        if ep(b)[0]<ep(a)[1]:raise ValueError('First tail endpoint must follow original Rin')
        algebra=lower[5]['S'].algebra;S=algebra.lift(0,1);RS=algebra.lift(0,1)
        records=[];sources=[];length=b-a
        for j in range(cells):
            phase=a+length*weighted.budget.coupled.partition(c,j,cells);cell=self.source_cell(Z,phase)
            sources.append(cell);rho=cell[2];value=cell[5]['S']
            # rho=4*exp(hb*s): dR=Lambda^-1*hb*rho*ds.
            # Keep the original hb once in both physical measures.
            S+=algebra.width(algebra.shift(value*(rho*length/cells),(1,0,0,-1)))
            RS+=algebra.width(algebra.shift(value*(rho*rho*length/cells),(1,0,0,-2)))
            records.append(dict(actual_first_phase_cell=phase,rho=rho,physical_S_Z1=op.export(value,1),
                physical_measure='dR=epsilon*hb*rho*ds',weighted_physical_measure='R*dR=epsilon^2*hb*rho^2*ds'))
        integrals=reduce_integrals(S,RS,lower[5],upper[5]);totals,mp,mchange=self.cutoff.support_moments(Z)
        prefixes=dict(m13=totals['m13']+integrals['N1p_integral'],m15=totals['m15']-integrals['R_N1p_integral'])
        result=dict(source_family=self.family,Z=lower[0],query=dict(chart='original_actual_first_collar_tail',
            lower_actual_phase=a,upper_actual_phase=b,lower_is_original_Rin=True,upper_is_not_full_leading_support=True),
            geometry=dict(lower_rho=lower[2],upper_rho=upper[2],same_original_rho_map='rho=4*exp(hb*s); R=epsilon*rho'),
            physical_fixed_leading_forcing_integrals_Z1={k:op.export(v,1) for k,v in integrals.items()},
            physical_reduced_S_and_RS_integrals_Z1=dict(S=op.export(S,1),RS=op.export(RS,1)),
            physical_updated_pressure_moment_prefixes_Z1={k:op.export(v,1) for k,v in prefixes.items()},
            pressure_moment_prefix_budgets={k:{'Z'+str(j):parent.parent.moment_budget(v,j) for j in (0,1)} for k,v in prefixes.items()},
            unchanged_full_velocity_support_moments={k:op.export(totals[k],1) for k in ('m11','m12','m14')},
            endpoint_identities=dict(unweighted='integral N1p=-integral S/2+[Q+R*Q_R]',
                weighted='integral R*N1p=-integral R*S/2+[R^2*Q_R]',common_Qaxis_cancels_analytically=True),
            original_N1p_known_forcing_bound_to_same_leading_field=True,
            own_Psi_mean_inlet_and_original_pressure_atoms_preserved=True,
            leading_tail_stops_at_same_family_finite_N_modification_origin=True,
            modified_finite_N_Q_mean_and_high_Z_source_required_beyond_rminus=True,
            radial_viscous_second_derivative_integrated_by_exact_FTC=True,
            finite_upper_boundary_terms_not_dropped=True,
            remaining_outer_leading_forcing_tails_uncomputed=True,
            pressure_prefixes_not_full_support_defects_or_repair=True,
            original_pressure_moment_m15_tail_sign_retained=True,
            genuine_whole_cell_ranges_and_actual_phase_Jacobians_used=True,
            whole_Z_relative_width_not_certified=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1LeadingForcingTailPrefixPacket();state=(S,RS,integrals,prefixes,totals,records)
        dependencies=(local,lower,upper,sources,mp,mchange)
        frozen=op.snapshot((result,state,[cell[:6] for cell in (lower,upper,*sources)]))
        self._packets[id(token)]=(token,result,state,dependencies,frozen);return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1LeadingForcingTailPrefixPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live original leading forcing prefix packet required')
        local,lower,upper,sources,mp,mchange=entry[3]
        if op.snapshot((entry[1],entry[2],[cell[:6] for cell in (lower,upper,*sources)]))!=entry[4]:
            raise ValueError('Leading forcing prefix or retained source range changed')
        self.mean.local.report(local)
        for cell in (lower,upper,*sources):
            if op.snapshot(cell[:6])!=cell[7]:raise ValueError('Bound leading forcing cell changed')
            self.inputs.report(cell[6])
        self.cutoff.moments.report(mp)
        if parent.change_snapshot(mchange)!=mchange[4]:raise ValueError('Incoming original cutoff moments changed')
        for cell in mchange[3]:
            for source,token in cell[-1]:source.report(token)
        self.mean.local.report(mchange[5])
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_cutoff):
    owner=CurrentOriginalN1LeadingForcingTailPrefix(accepted_cutoff,require_checked=False)
    a,_=owner.cutoff.selected_phase('.371')
    packets={name:owner.evaluate('.371',a*factor) for name,factor in (('unchanged_tail_midpoint',owner.c.mpf('1.5')),('rminus_endpoint',owner.c.mpf(2)))}
    raw=dict(source_family=owner.family,samples={k:owner.report(v) for k,v in packets.items()},
        input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_LEADING_FORCING_TAIL_PREFIX_COMPUTED',len(packets),flush=True);return owner,packets

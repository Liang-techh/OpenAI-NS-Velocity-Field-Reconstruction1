"""Paper section 15.1 n1 cutoff, own recovery and velocity-support moments.

Only F1/Uz1 are cut directly. Mean/V and pressure are recovered from their
own integrals. Three moments have full velocity support; two still require
actual leading pressure-forcing tails outside the common inner interval.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_cumulative_moments as parent
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=parent.op;regular=parent.regular;natural=parent.natural;weighted=parent.weighted
domain=parent.domain;ep=parent.ep;HERE,PREFIX,sha=parent.HERE,parent.PREFIX,parent.sha
NAME=PREFIX+'current_original_n1_radial_cutoff_extension.json.gz'
RECEIPT=PREFIX+'current_original_n1_radial_cutoff_extension_check.json'
GATE='current_original_n1_radial_cutoff_own_recovery_and_three_velocity_support_moments_installed'
OPEN=parent.OPEN


def cutoff(c,fraction):
    """Smooth, Z-independent cutoff: one to 1/2, zero from 3/4."""
    q=c.mpf(fraction)
    if ep(q)[1]<=ep(c.mpf('.5'))[0]:return c.mpf(1)
    if ep(q)[0]>=ep(c.mpf('.75'))[1]:return c.mpf(0)
    return 1-domain.sigma_jets(c,q*4-2)[0]


def cutoff_radial_jets(algebra,fraction,rho,inlet_phase,order=2):
    """Retain hb^-1/-2 and the original selected phase scale explicitly."""
    c=algebra.ctx;q=c.mpf(fraction);value=cutoff(c,q);zero=algebra.lift(0,order)
    if ep(q)[1]<=ep(c.mpf('.5'))[0] or ep(q)[0]>=ep(c.mpf('.75'))[1]:
        return [algebra.lift(value,order),zero,zero]
    sigma=domain.sigma_jets(c,q*4-2)
    first=algebra.width(algebra.lift(-4*sigma[1]/(inlet_phase*rho),order),-1)
    second=algebra.width(algebra.lift(4*sigma[1]/(inlet_phase*rho*rho),order),-1)
    second-=algebra.width(algebra.lift(32*sigma[2]/(inlet_phase*inlet_phase*rho*rho),order),-2)
    return [algebra.lift(value,order),first,second]


def change_densities(algebra,physical,leadingF,leadingU,rho,change):
    """Physical moment variations, expressed per drho, at fixed leading field."""
    F,U=physical['F1'],physical['Uz1']
    return dict(m11=algebra.shift(U*change,(0,0,0,-1)),
        m12=algebra.shift(F*(change*2*rho),(0,0,0,-2)),
        m13=algebra.shift(leadingF*F*(change*2),(0,0,0,-1)),
        m14=algebra.shift((leadingU*F+U*leadingF)*(change*2*rho),(0,0,0,-2)),
        m15=algebra.shift(leadingU*U*(change*2),(0,0,0,-1))
            -algebra.shift(leadingF*F*(change*2*rho),(0,0,0,-2)))


def recover_velocity(algebra,U,M,Z,rho,delta,axis=False):
    """Paper (15.16); radial primitive is reconstructed, never multiplied by chi."""
    c=algebra.ctx;z=algebra.lift(op.Jet.variable(c,Z,2));d=1-z*z;L=1-z*z*delta
    R=algebra.shift(algebra.lift(rho,2),(0,0,0,-1))
    V=op.truncate((U*(z*2)*R-M*(z*(1+delta))-op.dz(M)*d)/L,1)
    mean=algebra.lift(0,2) if axis else algebra.shift(M/rho,(0,0,0,1))
    return dict(mean1=mean,K1=mean-U,V1=V)


def cut_pair(value,first,second,chi):
    return dict(value=op.truncate(value,1)*op.truncate(chi[0],1),
        rho1=op.truncate(first,1)*op.truncate(chi[0],1)+op.truncate(value,1)*op.truncate(chi[1],1),
        rho2=op.truncate(second,1)*op.truncate(chi[0],1)+op.truncate(first,1)*op.truncate(chi[1],1)*2
             +op.truncate(value,1)*op.truncate(chi[2],1))


def change_snapshot(entry):
    return op.snapshot((entry[:3],[cell[:-1] for cell in entry[3]]))


@dataclass(frozen=True,eq=False)
class OriginalN1RadialCutoffPacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1RadialCutoffExtension:
    @source_precision
    def __init__(self,accepted_moments,require_checked=True):
        if type(accepted_moments) is not parent.CurrentOriginalN1CumulativeMoments or not accepted_moments.acceptance_loaded:
            raise ValueError('Accepted original n1 cumulative moment owner required')
        accepted_moments.assert_graph();self.moments=accepted_moments;self.second=accepted_moments.second
        self.regular=accepted_moments.regular;self.inputs=accepted_moments.inputs
        self.c=accepted_moments.c;self.physical=accepted_moments.physical
        self.family=copy.deepcopy(accepted_moments.family);self.hashes=dict(accepted_moments.hashes)
        for name in (parent.NAME,parent.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._parents=(self.moments,self.second,self.regular,self.inputs);self._identities=(self.c,self.physical)
        self._packets={};self._changes={};self._support={};self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Radial cutoff receipt scope differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Radial cutoff source changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if ((self.moments,self.second,self.regular,self.inputs)!=self._parents or (self.c,self.physical)!=self._identities or
            self.regular is not self.moments.regular or self.second is not self.moments.second or self.inputs is not self.moments.inputs or
            not self.moments.acceptance_loaded or self.family!=self.moments.family or self.hashes!=self._hashes):
            raise ValueError('Radial cutoff source graph changed')
        self.moments.assert_graph()

    def selected_phase(self,Z):
        packet=self.regular.mean.local.evaluate(Z,'1','first');view=self.regular.mean.local.report(packet)
        return view['query']['actual_source_coordinate'],packet

    @source_precision
    def changes(self,Z,fraction,cells=4):
        self.assert_graph();c=self.c;q=c.mpf(fraction)
        if isinstance(cells,bool) or not isinstance(cells,int) or not 1<=cells<=64:
            raise ValueError('Cutoff partition in[1,64] required')
        if ep(q)[0]<0 or ep(q)[1]>1:raise ValueError('Selected fraction in[0,1] required')
        if not isinstance(fraction,(str,int,float)) and ep(q)[0]!=ep(q)[1]:raise ValueError('One exact cutoff cumulative endpoint required')
        key=(c.mpf(Z)._mpi_,q._mpi_,cells)
        if key in self._changes:
            entry=self._changes[key]
            if change_snapshot(entry)!=entry[4]:raise ValueError('Cutoff cumulative changes or cell ranges altered')
            self.regular.mean.local.report(entry[5])
            for cell in entry[3]:
                for source,token in cell[-1]:source.report(token)
            return entry
        phase,lp=self.selected_phase(Z);z,algebra,_=self.inputs.frame(Z)
        sums={k:algebra.lift(0,1) for k in parent.MOMENTS};M2=algebra.lift(0,2);records=[];sources=[]
        if ep(q)[1]>ep(c.mpf('.5'))[0]:
            length=q-c.mpf('.5')
            for j in range(cells):
                native=c.mpf('.5')+length*weighted.budget.coupled.partition(c,j,cells)
                cell=self.moments.range_cell(z,phase*native,'first');sources.append(cell)
                rho,rows,fields,W,physical=cell[:5];change=cutoff(c,native)-1
                gradient=op.Jet(c,self.regular.weighted.source.amplitude.gradient_jets(z))
                leadingF=algebra.shift(rows['Phi'][0]*weighted.relative_amplitude(algebra,gradient,1,1),(0,1,0,0))
                densities=change_densities(algebra,physical,op.truncate(leadingF,1),op.truncate(fields['Uz'][0],1),rho,change)
                measure=rho*phase*length/cells
                for name,value in densities.items():sums[name]+=algebra.width(value*measure)
                seeds=natural.seeds(algebra,rows,fields,z,rho,self.physical.delta)
                gap=self.regular.tail['derivative_Cauchy_gap'];bound=self.regular.natural._data['initial_errors']['axial']
                bound+=parent.domain.upper(c,rho)*self.regular.tail['value_remainder_upper']
                error=algebra.lift(op.Jet(c,[parent.domain.upper(c,bound/gap**k) for k in range(3)]))
                U2=algebra.shift(op.truncate(seeds['Uz1'],2)+regular.symmetric_error(error),(0,0,0,-1))
                M2+=algebra.width(algebra.shift(U2*change,(0,0,0,-1))*measure)
                records.append(dict(selected_fraction_cell=native,actual_phase_cell=phase*native,
                    cutoff_range=cutoff(c,native),rho=rho,measure='dR=epsilon*hb*rho*ds',
                    physical_moment_change_densities_per_drho={k:op.export(v,1) for k,v in densities.items()},
                    own_M_change_Z2=op.export(M2,2)))
        entry=(sums,M2,records,sources,None,lp);entry=(*entry[:4],change_snapshot(entry),lp)
        self._changes[key]=entry;return entry

    @source_precision
    def support_moments(self,Z,cells=4):
        self.assert_graph()
        if isinstance(cells,bool) or not isinstance(cells,int) or not 1<=cells<=64:
            raise ValueError('Cutoff partition in[1,64] required')
        key=(self.c.mpf(Z)._mpi_,cells)
        if key in self._support:
            values,p,change,frozen=self._support[key]
            if op.snapshot(values)!=frozen:raise ValueError('Full velocity-support moments altered')
            self.moments.report(p)
            if self.changes(Z,'1',cells) is not change:raise ValueError('Full support source ledger replaced')
            return values,p,change
        p=self.moments.evaluate(Z,'1','first');self.moments.report(p)
        source=self.moments._packets[id(p)][3][2];change=self.changes(Z,'1',cells)
        values={k:source[k]+change[0][k] for k in parent.MOMENTS}
        self._support[key]=(values,p,change,op.snapshot(values))
        return values,p,change

    @source_precision
    def evaluate(self,Z,coordinate,chart='core',cells=4):
        self.assert_graph();began=time.monotonic();c=self.c
        p=self.second.evaluate(Z,coordinate,chart);sv=self.second.report(p);se=self.second._packets[id(p)]
        rp=se[2];rv=self.regular.report(rp);re=self.regular._packets[id(rp)]
        W,dx,dr,recovery,_,system=re[5];z=rv['Z'];rho=system['rho'];algebra=W[0].algebra
        phase,lp=self.selected_phase(z);fraction=c.mpf(coordinate) if chart=='first' else c.mpf(0)
        chi=(cutoff_radial_jets(algebra,fraction,rho,phase) if chart=='first' else
             [algebra.lift(1,2),algebra.lift(0,2),algebra.lift(0,2)])
        drr=se[3][9]
        pairs={name:cut_pair(W[i],dr[i],drr[i],chi) for name,i in (('F1',0),('Uz1',1))}
        physical={name:{kind:weighted.transform(v,system['gradient'],natural.WEIGHTS[i],1) for kind,v in rows.items()}
                  for name,i in (('F1',0),('Uz1',1)) for rows in (pairs[name],)}
        U2=algebra.shift(op.truncate(W[1],2)*chi[0],(0,0,0,-1))
        own_change=self.changes(z,coordinate if chart=='first' else '0',cells)
        Mbase=algebra.shift(op.truncate(W[1]+W[2],2)*rho,(0,0,0,-2));M=Mbase+own_change[1]
        # M is constant after the velocity cutoff. Use its common terminal
        # primitive rather than subtracting independently boxed late values.
        terminal_packet=None;terminal_change=None
        if chart=='first' and ep(fraction)[0]>=ep(c.mpf('.75'))[1]:
            terminal_packet=self.regular.evaluate(z,'.75','first');self.regular.report(terminal_packet)
            te=self.regular._packets[id(terminal_packet)];tw=te[5][0];trho=te[5][5]['rho']
            terminal_change=self.changes(z,'.75',cells)
            M=algebra.shift(op.truncate(tw[1]+tw[2],2)*trho,(0,0,0,-2))+terminal_change[1]
        own=recover_velocity(algebra,U2,M,z,rho,self.physical.delta,ep(rho)[1]==0)
        P=algebra.shift(op.truncate(W[3],1),(0,0,1,-1))+own_change[0]['m13']
        totals,mp,mchange=self.support_moments(z,cells)
        result=dict(source_family=self.family,Z=z,query=rv['query'],geometry=rv['geometry'],
            radial_cutoff=dict(keep_fraction='1/2',zero_fraction='3/4',value=ep(chi[0].terms.get((0,0,0,0),op.Jet.constant(c,0,2))[0]),
                original_selected_phase=phase,Z_independent=True,
                zero_radius_is_geometric_mean_of_Rkeep_Rin=True,zero_before_arithmetic_mean_by_AM_GM=True),
            normalized_cutoff_rho_jets=[op.export(v,1) for v in chi],
            physical_cut_F1_Uz1_value_rho1_rho2_Z1={k:{j:op.export(v,1) for j,v in rows.items()} for k,rows in physical.items()},
            physical_own_M_Z2=op.export(M,2),physical_own_mean_K_V={k:op.export(v) for k,v in own.items()},
            physical_recovered_pressure_on_common_inner=op.export(P,1),
            physical_cutoff_moment_values_at_Rin={k:op.export(v,1) for k,v in totals.items()},
            full_velocity_support_moments={k:op.export(totals[k],1) for k in ('m11','m12','m14')},
            full_velocity_support_moment_budgets={k:{'Z'+str(j):parent.moment_budget(totals[k],j) for j in (0,1)} for k in ('m11','m12','m14')},
            unresolved_leading_pressure_forcing_tails={
                'm13':'m13_full=m13_cut(Rin)+integral_Rin^infinity N1p0 dR',
                'm15':'m15_full=m15_cut(Rin)-integral_Rin^infinity R*N1p0 dR',
                'definition':'N1p0=-Omega0/(2R), original known_forcing N1p; F1_cut=Uz1_cut=0 outside Rin',
                'both_tails_uncomputed_and_not_set_zero':True},
            full_velocity_support_fixed_before_original_inner_exit=True,
            mean_V_pressure_recovered_not_multiplied_by_cutoff=True,
            original_hb_inverse_and_inverse_squared_cutoff_derivatives_retained=True,
            reconstructed_V_satisfies_own_n1_radial_recovery=True,
            untouched_core_and_keep_branch_preserved=True,
            leading_support_tail_two_moments_repair_n2_and_stress_still_open=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1RadialCutoffPacket();state=(chi,pairs,physical,M,own,P,totals)
        parents=(p,rp,lp,mp,own_change,mchange,terminal_packet,terminal_change)
        self._packets[id(token)]=(token,result,state,parents,op.snapshot((result,state)))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1RadialCutoffPacket or entry is None or entry[0] is not packet:
            raise ValueError('Live own n1 radial cutoff packet required')
        if op.snapshot((entry[1],entry[2]))!=entry[4]:raise ValueError('Own n1 cutoff packet changed')
        p,rp,lp,mp,change,mchange,tp,tchange=entry[3]
        self.second.report(p);self.regular.report(rp);self.regular.mean.local.report(lp);self.moments.report(mp)
        if tp is not None:self.regular.report(tp)
        for data in (change,mchange,tchange):
            if data is None:continue
            if change_snapshot(data)!=data[4]:raise ValueError('Bound cutoff correction ledger or range changed')
            self.regular.mean.local.report(data[5])
            for cell in data[3]:
                for source,token in cell[-1]:source.report(token)
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_moments):
    owner=CurrentOriginalN1RadialCutoffExtension(accepted_moments,require_checked=False)
    packets={name:owner.evaluate('.371',r,chart) for name,r,chart in (
        ('core_axis','0','core'),('core_exit','4','core'),('keep_endpoint','.5','first'),
        ('cutoff_interior','.625','first'),('zero_endpoint','.75','first'),('inner_exit','1','first'))}
    raw=dict(source_family=owner.family,samples={k:owner.report(v) for k,v in packets.items()},
        input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_RADIAL_CUTOFF_EXTENSION_COMPUTED',len(packets),flush=True);return owner,packets

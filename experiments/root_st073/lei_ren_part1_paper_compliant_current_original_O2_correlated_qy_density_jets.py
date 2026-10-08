"""Source-correlated q_y, actual mixed phase and nonlinear slow-y/yZ density.

The original positive q and eta are retained. A source-proved finite outer
range replaces only the excessively wide q_y enclosure. Bounds are not
selected source functions and are never differentiated.
"""
from dataclasses import replace
import gzip
import json
from pathlib import Path
from types import MappingProxyType
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_predicate_five_integrals as integrals

phase=integrals.phase;frames=phase.frames;source=phase.source
base,prior,ep=phase.base,phase.prior,phase.ep
MixedJet=phase.MixedJet;C0,Y,Z,YZ=phase.C0,phase.Y,phase.Z,phase.YZ
HERE,PREFIX,sha=phase.HERE,phase.PREFIX,phase.sha
NAME=PREFIX+'current_original_O2_correlated_qy_density_jets.json.gz'
RECEIPT=PREFIX+'current_original_O2_correlated_qy_density_jets_check.json'
GATE='original_O2_correlated_qy_actual_mixed_phase_and_five_slow_y_yZ_densities_enclosed'
ORDERS=(C0,Y,Z,YZ)


def rational(c,x):
    x=s.Rational(x);return c.mpf(int(x.p))/int(x.q)


def original_log_sigma_complement(c,y):
    """Directed logs of the exact original logistic odds at an interior point."""
    y=s.Rational(y)
    if not 0<y<1:raise ValueError('Interior exact source point required')
    yy=rational(c,y);odds=1/yy**2-1/(1-yy)**2
    if ep(odds)[0]>=0:
        tail=c.log1p(c.exp(-odds));return -odds-tail,-tail
    if ep(odds)[1]<=0:
        tail=c.log1p(c.exp(odds));return -tail,odds-tail
    raise ArithmeticError('Exact source odds sign required')


def correlated_qy_log_cap(c,count,index):
    """Paired flat products at endpoints; directed monotone factors inside.

    q^2 >= (3/5)*(1-sigma)/a, eta<=1, a>=4/5 gives
    |q_y| <= (12/5)*sqrt(5/3)*(5/4)^(3/2)
                *sigma*sqrt(1-sigma)*(y^-3+(1-y)^-3).
    The apparent endpoint singularities are bounded as paired products.
    """
    if type(count) is not int or count<4 or type(index) is not int or not 0<=index<count:
        raise ValueError('Original source count>=4 and cell index required')
    left,right=s.Rational(index,count),s.Rational(index+1,count)
    nn=c.mpf(count);reflection=nn/(nn-1)
    if index==0:
        Hlog=reflection**2-nn**2+c.ln(nn**3+reflection**3)
        recipe='sigma<=exp((1-y)^-2-y^-2); paired y^-3*exp(-1/y^2) increasing on first cell'
    elif index==count-1:
        Hlog=(reflection**2-nn**2)/2+c.ln(nn**3+reflection**3)
        recipe='1-sigma(y)=sigma(1-y); paired (1-y)^-3*exp(-1/(2*(1-y)^2)) increasing in reflected first cell'
    else:
        logsigma,_=original_log_sigma_complement(c,right)
        _,logcomplement=original_log_sigma_complement(c,left)
        S=lambda yy:yy**(-3)+(1-yy)**(-3)
        Smax=max(ep(S(rational(c,left)))[1],ep(S(rational(c,right)))[1])
        Hlog=logsigma+logcomplement/2+c.ln(c.mpf(Smax))
        recipe='sigma increasing; complement decreasing; S convex, maximum at an endpoint; no singular endpoint sampled'
    coefficient=c.mpf(12)/5*c.sqrt(c.mpf(5)/3)*(c.mpf(5)/4)**(c.mpf(3)/2)
    logcap=c.ln(coefficient)+c.mpf(ep(Hlog)[1])
    return logcap,dict(exact_y_cell=[str(left),str(right)],source_product='sigma*sqrt(1-sigma)*(y^-3+(1-y)^-3)',
        source_product_log_upper=ep(Hlog)[1],ordinary_qy_magnitude_log_upper=ep(logcap)[1],recipe=recipe,
        a_lower='4/5',positive_eta_upper='1 (accepted original eta<=1/2 is stronger)',
        source_identity='q^2=(3*(1-sigma)/5+eta)/a; q_y=-(1+eta)*a_y/(2*a^2*q)',
        paired_product_monotonicity_thresholds='y^2<=2/3 on left; (1-y)^2<=1/3 on right; count>=4',
        exact_flat_endpoint_qy_zero_but_whole_adjacent_cell_not_zero=True,
        no_evaluation_of_endpoint_singular_S=True,no_q_floor_or_eta_zero_replacement=True,
        bound_of_actual_derivative_not_derivative_of_cap=True)


class CorrelatedO2PredicateSources(frames.OriginalO2PredicateSources):
    def __init__(self):
        super().__init__();self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def frame(self,count,index,*,branch):
        original=super().frame(count,index,branch=branch)
        if original.record.get('correlated_qy_source_bound_installed'):return original
        a=original.roots['q'].atlas;c=a.ctx
        with mp.workdps(c.dps+40):
            logcap,proof=correlated_qy_log_cap(c,count,index)
            qy=prior.ScaledEnclosure(prior.FormalScale(a.bases,offset=logcap),c.mpf((-1,0)),a.ledger)
            roots=dict(original.roots);oldq=roots['q']
            roots['q']=MixedJet(a,{**oldq.rows,Y:qy})
            q=roots['q'][C0];lower=q.scale.evaluate()+c.ln(c.mpf(ep(q.coefficient)[0]))
            gamma_q=qy.positive_divide(q,lower)
            d=original.kernel.dstar
            u=MixedJet(a,{C0:original.u[C0],Y:original.u[C0]*a.add(original.gamma_g,gamma_q),
                Z:original.u[Z],YZ:a.add(roots['p2'][YZ]*q,roots['p2'][Z]*qy).positive_divide(d,d.scale.evaluate())})
            record=dict(original.record,actual_root_jets={name:row.record() for name,row in roots.items()},
                actual_u_mixed=u.record(),same_source_gamma_q=gamma_q.record(),correlated_qy_source_bound_installed=True,
                original_qy_uncorrelated_native_source=oldq[Y].record(),correlated_qy_source_proof=proof,
                original_positive_q_C0_and_eta_unchanged=True,original_nu_y_exact_source_cancellation_retained=True,
                original_transverse_q_rows_exact_zero_unchanged=True,
                actual_q_y_function_reenclosed_not_replaced_or_selected=True,
                correlated_qy_native_source=qy.record())
            result=replace(original,roots=MappingProxyType(roots),u=u,gamma_q=gamma_q,record=record)
            del self.frames[id(original)];self.frames[id(result)]=result
            self.cache[(count,index,branch)]=result;return result


class CorrelatedO2MixedDensity(phase.OriginalO2PredicateMixed):
    def __init__(self):
        self.source=CorrelatedO2PredicateSources();self.ctx=self.source.ctx;self.family=self.source.family
        self.hashes=dict(self.source.hashes)
        # Bind the accepted phase algorithms and completed C0/Z graph contract;
        # this does not rerun any ancestor suite or inherit a selected source.
        for module in (phase,integrals):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted same-source original mixed phase/density contract required')
            for name,digest in {**receipt['input_hashes'],module.RECEIPT:sha(module.RECEIPT)}.items():
                if sha(name)!=digest:raise ValueError('Original correlated source dependency changed: '+name)
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original correlated source closures disagree')
                self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def density(self,frame,*,N=160):
        integrals.candidate_N(N);self.source.describe(frame);a=frame.roots['q'].atlas
        primitive=self.primitive(frame,a.ctx.mpf((0,1)),phase=a.ctx.mpf((0,1)))
        changed=mixed_changed_densities(a,frame.roots['E'],frame.roots['V'],primitive['A'],primitive['B'],N=N)
        return dict(primitive=primitive,changed=changed)


def mixed_exponential_change(a,A,*,N):
    """Exact mixed chain, C0 expm1 via its directed defining mean."""
    integrals.candidate_N(N);c=a.ctx;x=A*(c.mpf(1)/N)
    actual=integrals.positive.bounded(x[C0])
    if ep(actual)[0]<-1 or ep(actual)[1]>1:
        raise ValueError('Source theorem |A/N|<=1 required before bounded exponential')
    finite=base.conditioned.clipped(c,actual,-1,1)
    mean=c.mpf(1);power=c.mpf(1)
    for k in range(1,65):power*=finite;mean+=power/c.factorial(k+1)
    tail=ep(c.exp(1)/c.factorial(66))[1]
    mean=base.conditioned.clipped(c,mean+c.mpf((-tail,tail)),c.exp(-1),c.exp(1))
    expx=a.scalar(c.exp(finite))
    return MixedJet(a,{C0:x[C0]*mean,Y:x[Y]*expx,Z:x[Z]*expx,
        YZ:a.add(x[YZ],x[Y]*x[Z])*expx})


def mixed_changed_densities(a,E,V,A,B,*,N):
    integrals.candidate_N(N);dE=E*mixed_exponential_change(a,A,N=N);dV=B*(a.ctx.mpf(1)/N)
    return dict(m=dV,h=dE,k=V*dE+E*dV+dE*dV,
        e=V*dV*2+dV*dV-E*dE-dE*dE*(a.ctx.mpf(1)/2),
        p=E*dE+dE*dE*(a.ctx.mpf(1)/2))


def run():
    began=time.monotonic();owner=CorrelatedO2MixedDensity();count=min(owner.source.source.parent.parent.levels)
    records=[];logcaps=[]
    with mp.workdps(owner.ctx.dps+40):
        for index in range(count):
            row={}
            for branch in integrals.BRANCHES:
                frame=owner.source.frame(count,index,branch=branch);result=owner.density(frame,N=160)
                proof=frame.record['correlated_qy_source_proof'];logcaps.append(proof['ordinary_qy_magnitude_log_upper'])
                row[branch]=dict(source=owner.source.describe(frame),actual_whole_phase_mixed=result['primitive']['record'],
                    actual_five_changed_density_C0_y_Z_yZ={name:jet.record() for name,jet in result['changed'].items()},
                    finite_N=160,normalized_B_already_over_Pstar=True,
                    all_original_density_nonlinear_cross_terms_retained=True,
                    slow_y_and_yZ_at_fixed_true_phase_not_full_rapid_spatial_derivative=True,
                    large_q_u_P_C_L_factors_remain_native_not_uniform_terminal_caps=True)
            records.append(row)
            if (index+1)%16==0:print('Source-correlated original O2 mixed densities:',index+1,'/',count,flush=True)
        report=dict(**{GATE:True},source_family=owner.family,explicit_candidate_N=160,
            exact_outer_source_domain=dict(y=['0','1'],Z=['-1','1'],phi=['0','1']),
            whole_original_O2_correlated_qy_mixed_density_cells=records,
            global_qy_source_magnitude_upper=owner.ctx.exp(owner.ctx.mpf(max(logcaps))),
            qy_ordinary_finite_outer_range_not_an_exact_function=True,
            source_ordinary_qy_q_power_zero_instead_of_minus_one=True,
            same_source_gamma_q_q_power_minus_one_instead_of_minus_two=True,
            native_complete_changed_five_slow_y_yZ_density_layer_installed=True,
            source_correlated_qy_propagated_into_actual_mixed_phase_and_density=True,
            source_positive_q_eta_nu_y_pressure_and_histories_unchanged=True,
            sharp_O2_phase_averaging_installed=False,actual_all_route_incoming_histories_installed=False,
            functional_terminal_identity_solved=False,current_whole_N_selected=False,
            all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,
            **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Source-correlated original q_y finite bound, actual varying-q fixed-phase A/B and complete five nonlinear C0/y/Z/yZ density ranges on the full original predicate union. Ordinary slow derivatives remain native; not sharp averaging, useful terminal norm bounds, actual all-route incoming, common global N, recursion or full reconstruction.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Original correlated q_y and actual five mixed densities generated',flush=True);return report


if __name__=='__main__':run()

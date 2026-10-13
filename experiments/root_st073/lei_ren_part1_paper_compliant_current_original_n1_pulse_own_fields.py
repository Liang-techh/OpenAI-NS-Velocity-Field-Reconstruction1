"""Actual main/gap/end pulse velocity, own mean and radial recovery.

Forward main history and backward end history use the selected original
coefficients. Physical Pstar is applied once. Inactive axial gaps keep their
nonzero cumulative mean; the auxiliary Q is not mistaken for Ur/R.
"""
import copy
import gzip
import json
from pathlib import Path
import time
import scipy.integrate  # Initialize the legacy dependency before its deep source chain.

import lei_ren_part1_paper_compliant_current_original_n1_pulse_C3_functions as jets
import lei_ren_part1_paper_compliant_current_original_C2_phase_primitives as c2

parent, source, controls, c3 = jets.parent, jets.source, jets.controls, jets.c3
HERE, PREFIX, sha = jets.HERE, jets.PREFIX, jets.sha
NAME=PREFIX+'current_original_n1_pulse_own_fields.json.gz'
RECEIPT=PREFIX+'current_original_n1_pulse_own_fields_check.json'
GATE='current_original_n1_main_gap_end_own_velocity_mean_Q_functions_installed'
OPEN=jets.OPEN


def encode(value):
    if type(value) is c2.C2Function:return [value.value.node,value.Z.node,value.ZZ.node]
    if isinstance(value,dict):return {key:encode(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [encode(item) for item in value]
    return jets.encode(value)


def radial_Q(g,delta,U,m):
    """Q=[2ZU-(1-delta)Zm-(1-Z^2)mZ]/(1-delta Z^2), as C2."""
    a=c2.C2Algebra(g);z=c2.C2Function(g.symbol('Z'),g.one,g.zero)
    one=a.fixed(1);d=a.add(one,a.neg(a.mul(z,z)))
    de=a.fixed(delta);L=a.add(one,a.neg(a.mul(de,a.mul(z,z))))
    u=c2.C2Function(U.value,U.Z,U.ZZ);mean=c2.C2Function(m.value,m.Z,m.ZZ)
    mz=c2.C2Function(m.Z,m.ZZ,m.ZZZ)
    numerator=a.add(a.scale(a.mul(z,u),2),
        a.neg(a.mul(a.add(one,a.neg(de)),a.mul(z,mean))),a.neg(a.mul(d,mz)))
    return a.div(numerator,L)


def pulse_fields(g,mu,delta,Pstar,native,selected,W,logRp):
    g.ids([mu,delta,Pstar,logRp])
    a=c3.C3Algebra(g);c=g.constant
    bp=g.add(c('1/2'),mu);lams=[g.sub(c('1/2'),g.mul(c(i),mu)) for i in (1,2)]
    Ep=a.scale(native['u'],Pstar)
    gpvar=W['original_gp_variable']
    gp_at=lambda x:g.node('function_substitution',expression=W['original_gp'].node,
        variable=gpvar.node,value=x.node,Z_independent_substitution=True)
    gpnode=g.nodes[W['original_gp'].node]
    primitive=next(source.FunctionRef(g,i) for i in gpnode['arguments']
        if g.nodes[i]['operation']=='definite_integral')
    primitive_at=lambda x:g.node('function_substitution',expression=primitive.node,
        variable=gpvar.node,value=x.node,Z_independent_substitution=True)
    sigma=lambda x:jets.complete.Gamma.sigma(g,x)
    sigmap=lambda x:g.node('analytic_unary',name='original_flat_sigma_prime',argument=x.node,
        definition='sigma(x)*(1-sigma(x))*(2/x^3+2/(1-x)^3) for 0<x<1; zero otherwise',
        derivative_source=c2.current.ast_binding(c2.flat_source.sigma_jets),
        ordinary_derivative='same original sigma Taylor coefficient[1]',
        exact_flat_source_derivative_not_a_bound=True)
    gp_prime=lambda x:g.sub(g.mul(sigma(g.mul(c(50),x)),sigma(g.sub(c(11),x))),
        g.mul(primitive_at(x),sigmap(g.sub(c(11),x))))
    integral=lambda body,var,lo,hi:controls.integral(g,body,g.nodes[var.node]['name'],lo,hi,
        actual_complete_own_pulse_history=True,endpoint_independent_of_Z=True,
        no_finite_history_cutoff=True,numerical_cap_not_defining_integrand=True)
    exponential=lambda x:g.unary('exp',x)
    E=lambda y:a.scale(Ep,exponential(g.neg(g.mul(bp,y))))
    ell,norm=W['ell'],W['normalization']
    beta=lambda local:g.quotient(g.node('compact_raw_beta',argument=g.quotient(local,ell,
        'positive original ell').node,definition='exp(-1/(1-r^2)) on abs(r)<1; zero otherwise'),
        g.mul(ell,norm),'positive original ell and beta normalization')
    beta_prime=lambda local:g.quotient(g.node('analytic_unary',name='compact_raw_beta_derivative',
        argument=g.quotient(local,ell,'positive original ell').node,
        definition='-2*r*exp(-1/(1-r^2))/(1-r^2)^2 on abs(r)<1; zero otherwise',
        derivative_source=c2.current.ast_binding(c2.flat_source.beta_jets),
        ordinary_derivative='same original beta Taylor coefficient[1]',
        exact_original_derivative_not_its_global_cap=True),
        g.mul(ell,ell,norm),'positive ell squared and beta normalization')
    endscale=exponential(W['log_end_scale'])
    centers=(-3,-1);C=(selected['C1'],selected['C2']);out={}
    for chart in ('main','gap','end'):
        query=g.symbol('original_own_pulse_'+('end_s' if chart=='end' else chart+'_xi'))
        y=(g.add(g.quotient(c(13),mu,'positive original mu'),query) if chart=='end'
            else g.quotient(query,mu,'positive original mu'))
        theta=E(y);means=[]
        if chart=='main':
            B=a.scale(selected['ap'],gp_at(query))
            By=a.scale(selected['ap'],g.mul(mu,gp_prime(query)))
            v=g.symbol('original_own_main_history_v')
            for key,lam in zip(('m1','m2'),lams):
                kernel=integral(g.mul(exponential(g.neg(g.mul(lam,g.sub(y,v)))),gp_at(g.mul(mu,v))),
                    v,g.zero,y)
                means.append(a.add(a.scale(native[key],exponential(g.neg(g.mul(lam,y)))),a.scale(selected['ap'],kernel)))
            domain=dict(coordinate='xi=mu*log(R/Rp)',lower=g.zero,upper=c(11),dR_over_dxi='R/mu')
        else:
            v=g.symbol('original_own_end_history_v')
            B=a.fixed(0);By=a.fixed(0)
            if chart=='end':
                for coeff,center in zip(C,centers):
                    local=g.sub(query,c(center))
                    B=a.add(B,a.scale(coeff,beta(local)))
                    By=a.add(By,a.scale(coeff,beta_prime(local)))
                B=a.scale(B,endscale);By=a.scale(By,endscale)
                for lam in lams:
                    own=a.fixed(0)
                    for coeff,center in zip(C,centers):
                        kernel=integral(g.mul(exponential(g.mul(lam,g.sub(v,query))),beta(g.sub(v,c(center)))),
                            v,query,g.zero)
                        own=a.add(own,a.scale(coeff,g.neg(kernel)))
                    means.append(a.scale(own,endscale))
                domain=dict(coordinate='s=log(R/Rv)',lower=c(-4),upper=g.zero,dR_over_ds='R')
            else:
                distance=g.quotient(g.sub(c(13),query),mu,'positive original mu')
                weights=(W['first_row'],[g.add(A,g.mul(mu,D)) for A,D in zip(W['first_row'],W['divided_row'])])
                for lam,ww in zip(lams,weights):
                    own=a.add(*(a.scale(coeff,g.neg(weight)) for coeff,weight in zip(C,ww)))
                    # Combine the exact logarithms before applying exp.
                    factor=exponential(g.add(W['log_end_scale'],g.mul(lam,distance)))
                    means.append(a.scale(own,factor))
                domain=dict(coordinate='xi=mu*log(R/Rp)',lower=c(11),
                    upper=g.sub(c(13),g.mul(c(4),mu)),dR_over_dxi='R/mu')
        U=a.mul(theta,B);Uy=a.mul(theta,a.add(By,a.neg(a.scale(B,bp))))
        m=a.mul(theta,means[0]);my=a.add(U,a.neg(m))
        Q=radial_Q(g,delta,U,m);Qy=radial_Q(g,delta,Uy,my)
        logR=g.add(logRp,y);R=exponential(logR)
        b=c2.C2Algebra(g)
        Ur=b.scale(Q,g.unary('sqrt',g.mul(c('1/2'),R)))
        # The original pressure datum and cumulative swirl history are
        # absolute histories, retained rather than reset at a pulse boundary.
        prate=g.add(g.one,g.mul(c(2),mu))
        mass=g.mul(y,g.unary('exprel',g.neg(g.mul(prate,y))))
        pressure=a.scale(a.add(native['P0'],native['Mp'],a.scale(a.mul(native['u'],native['u']),g.mul(c('1/2'),mass))),
            g.mul(Pstar,Pstar))
        out[chart]=dict(query=query,domain=domain,y=y,log_radius=logR,radius=R,
            Utheta=theta,U=U,U_y=Uy,m=m,m_y=my,Q=Q,Q_y=Qy,Ur=Ur,
            normalized_axial=B,normalized_axial_y=By,normalized_linear_histories=means,
            absolute_swirl_base_pressure=pressure,
            pressure_is_not_n1_pressure_forcing=True,physical_Pstar_applied_exactly_once=True,
            Q_is_auxiliary_Ur_over_sqrt_R_half=True,own_history_not_zeroed_on_gap=True)
    return out


class CurrentOriginalN1PulseOwnFields:
    def __init__(self,jets_owner=None,require_checked=True):
        self.jets=jets_owner or jets.CurrentOriginalN1PulseC3Functions()
        if not self.jets.acceptance_loaded:raise ValueError('Checked actual pulse C3 owner required')
        self.jets.assert_graph();self.g=self.jets.g;self.family=copy.deepcopy(self.jets.family)
        self.hashes=dict(self.jets.hashes)
        for name,digest in parent.linear.read(jets.RECEIPT)['input_hashes'].items():
            if sha(name)!=digest or (name in self.hashes and self.hashes[name]!=digest):
                raise ValueError('Changed actual pulse C3 source '+name)
            self.hashes[name]=digest
        self.hashes[jets.NAME]=sha(jets.NAME);self.hashes[jets.RECEIPT]=sha(jets.RECEIPT)
        complete=self.jets.complete;Gamma=complete.Gamma
        logP=source.FunctionRef(self.g,complete.Gamma.pulse.preRp.operator.background['parameters']['logP'])
        self.fields=pulse_fields(self.g,Gamma.schedule['mu'],Gamma.schedule['delta'],self.g.unary('exp',logP),
            self.jets.native,self.jets.functions['pulse'],Gamma.pulse.functions,Gamma.schedule['logRp'])
        self.hashes[Path(c2.__file__).name]=sha(Path(c2.__file__).name)
        self.hashes[Path(c2.flat_source.__file__).name]=sha(Path(c2.flat_source.__file__).name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(self.g.nodes);self.snapshot=encode(self.fields);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if require_checked:
            receipt,raw=parent.linear.read(RECEIPT),parent.linear.read(NAME)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family
                    or any(receipt.get(key) for key in OPEN) or raw['exact_graph_nodes']!=self.g.nodes
                    or raw['own_main_gap_end_fields']!=self.snapshot):
                raise ValueError('Checked actual own pulse fields required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual own pulse field source '+name)
            self.acceptance_loaded=True

    def assert_graph(self):
        self.jets.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or encode(self.fields)!=self.snapshot or self.hashes!=self._hashes:
            raise ValueError('Actual own pulse field graph changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1PulseOwnFields(require_checked=False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,own_main_gap_end_fields=owner.snapshot,
        actual_selected_ap_end_coefficients_consumed=True,native_nonzero_mean_and_pressure_inlet_retained=True,
        no_extra_Pstar_division=True,missing_Q_third_derivative_not_zero_padded=True,
        **dict.fromkeys((GATE,)+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_ACTUAL_OWN_PULSE_FIELDS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':run()

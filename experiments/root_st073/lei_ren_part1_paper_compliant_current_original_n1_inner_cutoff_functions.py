"""True n1 flat cutoff, own recovery and five signed inner endpoint functions.

Only F1 and Uz1 are cut. Mean/M, K/Q/V and pressure are recovered with
their complete core-to-first histories; leading N1p remains unchanged.
The five inner endpoint functions are not the full repaired five defects.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions as radial

inner=radial.inner;parent,source,controls,c3,own=radial.parent,radial.source,radial.controls,radial.c3,radial.own
HERE,PREFIX,sha=radial.HERE,radial.PREFIX,radial.sha
NAME=PREFIX+'current_original_n1_inner_cutoff_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_inner_cutoff_functions_check.json'
GATE='current_original_n1_true_inner_cutoff_and_five_endpoint_functions_installed'
OPEN=radial.OPEN
MOMENTS=('m11','m12','m13','m14','m15')


def original_cutoff_rows(g,s,rho,parameters):
    """Ordinary rho derivatives of chi=1-sigma(4s/sc_quarter-2)."""
    a=c3.C3Algebra(g);c=g.constant;phase=parameters['sc_quarter'];hb=parameters['hb']
    q=g.quotient(s,phase,'positive original common-inner phase')
    arg=g.sub(g.mul(c(4),q),c(2));sigma=own.jets.complete.Gamma.sigma(g,arg)
    provider=own.c2.current.ast_binding(own.c2.flat_source.sigma_jets)
    prime=g.node('analytic_unary',name='original_flat_sigma_prime',argument=arg.node,
        derivative_source=provider,ordinary_derivative='same original sigma Taylor coefficient[1]',
        exact_flat_source_derivative_not_a_bound=True)
    second=g.node('analytic_unary',name='original_flat_sigma_second',argument=arg.node,
        derivative_source=provider,ordinary_derivative='2*same original sigma Taylor coefficient[2]',
        Taylor_to_ordinary_factor=2,exact_flat_source_derivative_not_a_bound=True)
    div=lambda n,d:g.quotient(n,d,'positive original phase hb and rho factors')
    value=g.sub(g.one,sigma)
    first=g.neg(div(g.mul(c(4),prime),g.mul(phase,hb,rho)))
    second_rho=g.sub(div(g.mul(c(4),prime),g.mul(phase,hb,rho,rho)),
        div(g.mul(c(16),second),g.mul(phase,phase,hb,hb,rho,rho)))
    return [a.fixed(value),a.fixed(first),a.fixed(second_rho)],dict(fraction=q,argument=arg,
        keep_fraction='1/2',zero_fraction='3/4',phase_is_sc_quarter=True,
        ordinary_sigma_second_uses_two_times_Taylor_coefficient=True,
        exact_hb_inverse_and_inverse_squared_retained=True,Z_independent_cutoff=True)


def multiply_cut_rows(g,rows,chi):
    a=c3.C3Algebra(g)
    return [a.mul(rows[0],chi[0]),
        a.add(a.mul(rows[1],chi[0]),a.mul(rows[0],chi[1])),
        a.add(a.mul(rows[2],chi[0]),a.scale(a.mul(rows[1],chi[1]),2),a.mul(rows[0],chi[2]))]


def signed_rho_integrands(g,leading,cut,system,epsilon):
    """Paper (15.1), all five genuine signed densities per drho."""
    a=c3.C3Algebra(g);r=system['rho'];eps2=g.mul(epsilon,epsilon)
    F,U=cut['F1'][0],cut['Uz1'][0]
    f,u=leading['F_rho_rows'][0],leading['U_rho_rows'][0]
    pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])
    return dict(m11=a.scale(U,epsilon),m12=a.scale(F,g.mul(g.constant(2),eps2,r)),
        m13=a.scale(pressure,epsilon),
        m14=a.scale(a.add(a.mul(u,F),a.mul(U,f)),g.mul(g.constant(2),eps2,r)),
        m15=a.add(a.scale(a.mul(u,U),g.mul(g.constant(2),epsilon)),
            a.neg(a.scale(pressure,g.mul(eps2,r)))))


def own_cutoff_fields(g,inner_owner,radial_rows,cutoff_rows):
    """Cut the two velocities, integrate the mean and pressure once."""
    a=c3.C3Algebra(g);c=g.constant;p=inner_owner.parameters;eps,hb=p['epsilon'],p['hb']
    fields={};integrands={};cumulative={};bodies={}
    rcore=inner_owner.coordinates['core'];s=inner_owner.coordinates['first']
    for chart in ('core','first'):
        system=inner_owner.systems[chart];r=system['rho'];state=inner_owner.functions[chart]
        cut={key:multiply_cut_rows(g,radial_rows[chart][key],cutoff_rows[chart]) for key in ('F1','Uz1')}
        if chart=='core':
            mean=radial_rows[chart]['mean1'][0];pressure=state['P1']
        else:
            v=g.symbol('original_n1_cutoff_mean_v');rv=g.mul(c(4),g.unary('exp',g.mul(hb,v)))
            inlet=inner.collar.at(g,radial_rows['core']['mean1'][0],rcore,c(4))
            body=a.scale(inner.collar.at(g,cut['Uz1'][0],s,v),g.mul(hb,rv))
            history=inner.collar.integral(g,body,v,g.zero,s,actual_cutoff_own_mean=True,
                history_constant_after_velocity_zero_by_flat_support=True,core_inlet_not_reset=True)
            mean=a.div(a.add(a.scale(inlet,4),history),a.fixed(r))
            t=g.symbol('original_n1_cutoff_pressure_t');rt=g.mul(c(4),g.unary('exp',g.mul(hb,t)))
            pr=radial.pressure_radial_function(g,inner_owner.leading['first']['F_rho_rows'][0],cut['F1'][0],system['N1p'],eps)
            pbody=a.scale(inner.collar.at(g,pr,s,t),g.mul(hb,rt))
            pressure=a.add(inner.collar.at(g,inner_owner.functions['core']['P1'],rcore,c(4)),
                inner.collar.integral(g,pbody,t,g.zero,s,actual_cutoff_own_pressure=True,
                    leading_N1p_not_cut=True,nonzero_core_pressure_inlet_retained=True))
            bodies[chart]=dict(mean_integrand=body,pressure_integrand=pbody)
        cut_state=dict(Uz1=cut['Uz1'][0],K1=a.add(mean,a.neg(cut['Uz1'][0])))
        if chart=='core':
            # Preserve accepted even axis continuations; the core cutoff is
            # exactly one, so no quotient by rho=0 belongs in this branch.
            cut.update(mean1=radial_rows[chart]['mean1'],K1=radial_rows[chart]['K1'])
        else:cut.update(radial.own_mean_rows(g,cut_state,cut,r))
        recovered=radial.recovered_radial_rows(g,cut,r,eps,p['delta'])
        fields[chart]=dict(rows=cut,radial_primitive_M1=a.scale(mean,g.mul(eps,r)),
            pressure_P1=pressure,**recovered,
            pressure_rho=radial.pressure_radial_function(g,inner_owner.leading[chart]['F_rho_rows'][0],
                cut['F1'][0],system['N1p'],eps))
        integrands[chart]=signed_rho_integrands(g,inner_owner.leading[chart],cut,system,eps)
        if chart=='core':
            t=g.symbol('original_n1_cutoff_core_moment_t')
            functions={key:inner.collar.integral(g,inner.collar.at(g,q,rcore,t),t,g.zero,rcore,
                actual_signed_n1_inner_moment=True,paper_moment=key,axis_integrand_regular=True)
                for key,q in integrands[chart].items() if key not in ('m11','m13')}
        else:
            t=g.symbol('original_n1_cutoff_first_moment_t');rt=g.mul(c(4),g.unary('exp',g.mul(hb,t)))
            functions={key:a.add(inner.collar.at(g,cumulative['core'][key],rcore,c(4)),
                inner.collar.integral(g,a.scale(inner.collar.at(g,q,s,t),g.mul(hb,rt)),t,g.zero,s,
                    actual_signed_n1_inner_moment=True,paper_moment=key,
                    full_core_to_first_history_retained=True,leading_N1p_not_cut=True))
                for key,q in integrands[chart].items() if key not in ('m11','m13')}
        cumulative[chart]={**functions,'m11':fields[chart]['radial_primitive_M1'],'m13':pressure}
    return fields,integrands,cumulative,bodies


class CurrentOriginalN1InnerCutoffFunctions:
    def __init__(self,radial_owner=None,require_checked=True):
        self.radial=radial_owner or radial.CurrentOriginalN1InnerRadialFunctions()
        if not self.radial.acceptance_loaded:raise ValueError('Checked actual n1 radial functions required')
        self.radial.assert_graph();self.inner=self.radial.inner;self.g=g=self.radial.g
        self.family=copy.deepcopy(self.radial.family);self.hashes=dict(self.radial.hashes)
        for name,digest in parent.linear.read(radial.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (radial.NAME,radial.RECEIPT):self._pin(name,sha(name))
        self.parameters=self.inner.parameters;s=self.inner.coordinates['first'];a=c3.C3Algebra(g)
        first,geometry=original_cutoff_rows(g,s,self.inner.systems['first']['rho'],self.parameters)
        self.cutoff_rows=dict(core=[a.fixed(1),a.fixed(0),a.fixed(0)],first=first)
        self.fields,self.integrands,self.cumulative,self.history_bodies=own_cutoff_fields(g,self.inner,self.radial.rows,self.cutoff_rows)
        self.endpoint={key:inner.collar.at(g,self.cumulative['first'][key],s,self.parameters['sc_quarter']) for key in MOMENTS}
        self.geometry=geometry
        self.definition=dict(only_F1_Uz1_cut=True,mean_K_Q_V_pressure_recovered_from_own_history=True,
            n1_temporal_mean_weight='1+delta',leading_pressure_forcing_unchanged=True,
            five_signed_inner_endpoint_functions=True,endpoint_is_Rin=True,
            m13_includes_two_F0_F1_and_N1p=True,m15_includes_negative_R_pressure_derivative=True,
            velocity_moment_histories_constant_after_cutoff_zero=True,
            pressure_and_two_pressure_moments_keep_leading_source_after_cutoff_zero=True,
            exterior_m13_m15_leading_pressure_tails_still_required=True,
            full_five_defects_and_linear_repair_installed=False,original_numeric_n1_solution_oracle_executed=False)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if require_checked:
            raw_data,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or raw_data['exact_graph_nodes']!=g.nodes
                or raw_data['inner_cutoff_bindings']!=self.snapshot):raise ValueError('Checked actual n1 cutoff functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed true cutoff source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed cutoff dependency '+name)
        self.hashes[name]=digest

    def _snapshot(self):return own.encode(dict(cutoff_rows=self.cutoff_rows,fields=self.fields,integrands=self.integrands,
        cumulative=self.cumulative,endpoint=self.endpoint,history_bodies=self.history_bodies,
        geometry=self.geometry,definition=self.definition))

    def assert_graph(self):
        self.radial.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes:
            raise ValueError('True cutoff graph changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1InnerCutoffFunctions(require_checked=False)
    data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,inner_cutoff_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys((GATE,)+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_TRUE_INNER_CUTOFF_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,data


if __name__=='__main__':run()

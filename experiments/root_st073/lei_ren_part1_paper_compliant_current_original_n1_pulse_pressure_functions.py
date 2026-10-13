"""Complete original Rp-to-Rv pulse pressure value/first-Z functions.

Actual own U C3 and auxiliary Q C2 feed the ten signed reduced-S
integrands. Physical radial measures appear once; internal endpoints
telescope by the selected source moment identities. The paper m15 already
includes its minus relative to the literal weighted forcing integral.
"""
import copy
import gzip
import json
from pathlib import Path
import time
import scipy.integrate

import lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields as own

parent,source,controls=own.parent,own.source,own.controls
physical=parent.active.physical
HERE,PREFIX,sha=own.HERE,own.PREFIX,own.sha
NAME=PREFIX+'current_original_n1_pulse_pressure_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_pulse_pressure_functions_check.json'
GATES=('current_original_n1_Rp_to_Rv_pulse_pressure_C1_functions_installed',
       'current_original_n1_full_pulse_C1_function_adapter_installed')
OPEN=tuple(key for key in own.OPEN if key not in GATES)


def pulse_pressure_functions(g,delta,mu,fields,outer_endpoints=None):
    """All three complete charts; C1 derivatives use genuine Q_ZZ."""
    g.ids([delta,mu]);a=parent.C1Algebra(g);products={};integrals={};endpoints={};pieces={}
    for chart in ('main','gap','end'):
        f=fields[chart];x=f['query'];R=f['radius'];domain=f['domain']
        g.ids([x,R,domain['lower'],domain['upper']])
        jacobian=R if chart=='end' else g.quotient(R,mu,'positive original mu')
        products[chart]={key:source.C1Function(p.value,p.Z)
            for key,p in physical.nonradial_integrands(g,f['U'],f['Q']).items()}
        integrals[chart]={}
        for key,p in products[chart].items():
            weight=jacobian if key.startswith('I0_') else g.mul(R,jacobian)
            integrals[chart][key]=source.C1Function(*(controls.integral(g,g.mul(weight,row),
                g.nodes[x.node]['name'],domain['lower'],domain['upper'],
                original_complete_pulse_chart=chart,original_pulse_own_history_retained=True,
                auxiliary_pressure_q_equals_own_Q=True,physical_Pstar_not_reapplied=True,
                measure=('R*ds' if chart=='end' else '(R/mu)*dxi') if key.startswith('I0_')
                    else ('R^2*ds' if chart=='end' else '(R^2/mu)*dxi'),
                endpoint_independent_of_Z=True,ordinary_slow_Z_derivative_order=n,
                no_finite_history_or_chart_cutoff=True,
                signed_integrand_not_its_magnitude_cap=True)
                for n,row in enumerate((p.value,p.Z))))
        ends=[]
        for edge in (domain['lower'],domain['upper']):
            at=lambda p:a.at(source.C1Function(p.value,p.Z),x,edge)
            radius=a.at(a.fixed(R),x,edge).value
            ends.append(parent.endpoint_functions(g,at(f['U']),at(f['Q']),at(f['Q_y']),radius))
        endpoints[chart]=ends
        pieces[chart]=parent.pressure_functions(g,delta,integrals[chart],*ends)
    total={key:a.add(*(integrals[chart][key] for chart in ('main','gap','end')))
        for key in integrals['main']}
    outer=outer_endpoints or (endpoints['main'][0],endpoints['end'][1])
    functions=parent.pressure_functions(g,delta,total,*outer)
    return functions,pieces,dict(products=products,integrals=integrals,total_integrals=total,
        per_chart_finite_endpoints=endpoints,global_finite_endpoints=outer,
        global_functions_use_telescoped_internal_endpoints=True,
        original_m15_negative_weighted_integral_sign_applied_once=True)


def actual_outer_endpoints(owner):
    """Exact flat startup and empty-support end, retaining native Rp mean."""
    g=owner.g;a=parent.C1Algebra(g);c3=own.c3.C3Algebra(g)
    Gamma=owner.jets.complete.Gamma;delta=Gamma.schedule['delta']
    logP=source.FunctionRef(g,Gamma.pulse.preRp.operator.background['parameters']['logP'])
    Ep=c3.scale(owner.jets.native['u'],g.unary('exp',logP))
    mean=c3.mul(Ep,owner.jets.native['m1'])
    Q=own.radial_Q(g,delta,c3.fixed(0),mean)
    q=source.C1Function(Q.value,Q.Z);zero=a.fixed(0)
    Rp=g.unary('exp',Gamma.schedule['logRp'])
    lower=parent.endpoint_functions(g,zero,q,a.neg(q),Rp)
    lower['q_plus_q_y']=zero  # Exact startup U=Uy=UZ=0; the native Q remains nonzero.
    upper={key:zero for key in lower}
    certificate=dict(startup_gp_and_gp_prime_exact_zero=True,
        actual_native_physical_mean_retained=own.encode(mean),
        incoming_Q=own.encode(Q),startup_Q_y_equals_minus_Q=True,
        terminal_beta_centers=[-3,-1],terminal_beta_radius=own.encode(Gamma.pulse.functions['ell']),
        terminal_s=0,terminal_beta_and_beta_prime_exact_zero=True,
        backward_own_histories_have_identical_zero_integration_endpoints=True,
        terminal_U_m_Q_Qy_all_Z_functions_exact_zero=True,
        internal_row_1='A[C]=-mu*(incoming1+ap*Pi1)',
        internal_row_2='(A+mu*D)[C]=-mu*(incoming2+ap*Pi2)',
        exact_second_row='A+mu*D',
        selected_inverse_and_kernel_integrals_from_checked_ancestor=True,
        mean_radial_FTC='m_y=U-m',
        q_plus_q_y_boundary_telescopes=True,R_q_y_boundary_telescopes=True)
    return (lower,upper),certificate


class CurrentOriginalN1PulsePressureFunctions:
    def __init__(self,own_owner=None,require_checked=True):
        self.own=own_owner or own.CurrentOriginalN1PulseOwnFields()
        if not self.own.acceptance_loaded:raise ValueError('Checked actual own pulse fields required')
        self.own.assert_graph();self.g=self.own.g;self.family=copy.deepcopy(self.own.family)
        self.hashes=dict(self.own.hashes)
        for name,digest in parent.linear.read(own.RECEIPT)['input_hashes'].items():
            if sha(name)!=digest or (name in self.hashes and self.hashes[name]!=digest):
                raise ValueError('Changed own pulse field source '+name)
            self.hashes[name]=digest
        self.hashes[own.NAME]=sha(own.NAME);self.hashes[own.RECEIPT]=sha(own.RECEIPT)
        ends,self.endpoint_certificate=actual_outer_endpoints(self.own)
        schedule=self.own.jets.complete.Gamma.schedule
        self.functions,self.pieces,self.bindings=pulse_pressure_functions(self.g,
            schedule['delta'],schedule['mu'],self.own.fields,ends)
        for module in (physical,parent):self.hashes[Path(module.__file__).name]=sha(Path(module.__file__).name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(self.g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            if (not receipt.get('all_passed') or not all(receipt.get(key) for key in GATES)
                    or any(receipt.get(key) for key in OPEN) or receipt['source_family']!=self.family
                    or raw['exact_graph_nodes']!=self.g.nodes or raw['pulse_pressure_bindings']!=self.snapshot):
                raise ValueError('Checked actual complete pulse pressure C1 functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed pulse pressure source '+name)
            self.acceptance_loaded=True

    def _snapshot(self):
        return own.encode(dict(functions=self.functions,pieces=self.pieces,bindings=self.bindings,
            endpoint_certificate=self.endpoint_certificate))

    def assert_graph(self):
        self.own.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes:
            raise ValueError('Actual complete pulse pressure graph changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1PulsePressureFunctions(require_checked=False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,
        pulse_pressure_bindings=owner.snapshot,radial_scope='original Rp<=R<=Rv, all main/gap/end charts',
        true_physical_pressure_value_and_first_Z_functions=True,
        complete_nonzero_own_history_and_finite_mu_retained=True,
        base_swirl_pressure_not_substituted_for_n1_forcing=True,
        numerical_enclosure_centers_not_function_values=True,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_COMPLETE_PULSE_PRESSURE_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':run()

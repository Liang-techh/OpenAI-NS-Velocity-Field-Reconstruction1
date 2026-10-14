"""True transported n1 terminal integrals and their certified zero tails.

The defining terminal functions stay attached to the original signed
integrals. A checked exact change of variables identifies those functions
with the selected physical repair response before any zero is exposed.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_reserved_own_fields as own

source,controls,parent=own.source,own.controls,own.parent
HERE,PREFIX,sha=own.HERE,own.PREFIX,own.sha
NAME=PREFIX+'current_original_n1_terminal_integral_identities.json.gz'
RECEIPT=PREFIX+'current_original_n1_terminal_integral_identities_check.json'
GATES=('current_original_n1_actual_terminal_integrals_reduce_to_reserved_response',
    'n1_repaired_own_terminal_five_function_zeros',
    'current_original_n1_repaired_mass_C1_and_pressure_C1_tail_functions_installed')
OPEN=tuple(k for k in own.OPEN if k not in GATES)
MOMENTS=own.MOMENTS


def mass_tail_fields(g,mass,R,delta):
    """Recovery after the true compact n1 velocities have all vanished."""
    a=parent.C1Algebra(g);mean=a.div(mass,a.fixed(R),'positive physical mass-tail radius')
    z=g.symbol('Z');d=g.sub(g.one,g.mul(z,z));L=g.sub(g.one,g.mul(delta,z,z))
    Q=g.quotient(g.neg(g.add(g.mul(g.add(g.one,delta),z,mean.value),g.mul(d,mean.Z))),L,
        'original positive L=1-delta*Z^2 on common axial domain')
    zero=a.fixed(0)
    return dict(M1=mass,mean1=mean,K1=mean,Uz1=zero,Utheta1=zero,F1=zero,Q1=Q,V1=g.mul(R,Q))


class CurrentOriginalN1TerminalIntegralIdentities:
    def __init__(self,own_owner,require_checked=True):
        if type(own_owner) is not own.CurrentOriginalN1ReservedOwnFields or not own_owner.acceptance_loaded:
            raise TypeError('Checked actual reserved own-field histories required')
        own_owner.assert_graph();self.own=own_owner;self.full=own_owner.full
        self.g=g=own_owner.g;self.family=copy.deepcopy(own_owner.family);self.hashes=dict(own_owner.hashes)
        for name in (own.NAME,own.RECEIPT):self._pin(name,sha(name))
        receipt=parent.linear.read(own.RECEIPT)
        for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
        if (not receipt['checks']['densities']['passed'] or not receipt['checks']['recovery']['passed']
            or not receipt['checks']['source']['full_band_pressure_prefix_identical_to_original_four_C1_rows']):
            raise ValueError('True signed densities/recovery/pressure prefix proof required')
        a=parent.C1Algebra(g);end=g.unary('log',g.constant(2))
        self.patch_endpoint={key:a.at(pair,own_owner.s,end) for key,pair in own_owner.patch_histories.items()}
        self.terminal_sources=dict(own_owner.terminal)
        self.reduced_patch={key:a.add(*(a.mul(B,C) for B,C in zip(self.full.response_matrix[key],
            self.full.coefficients.values()))) for key in MOMENTS}
        self.reduced_terminal={key:a.add(self.full.moments[key],self.reduced_patch[key]) for key in MOMENTS}
        schedule=self.full.Gamma.schedule
        logRv=g.add(schedule['logRp'],g.quotient(g.constant(13),schedule['mu'],'original positive pulse mu'))
        self.tail_origins=dict(mass_log_radius=self.full.linear.functions['log2Rc'],pressure_log_radius=logRv,
            pressure_source_support=copy.deepcopy(self.full.beyond_Rv_support))
        self.mass_coordinate=g.symbol('original_repaired_n1_log_R_over_2Rc')
        self.pressure_coordinate=g.symbol('original_repaired_n1_log_R_over_Rv')
        R=g.unary('exp',g.add(self.tail_origins['mass_log_radius'],self.mass_coordinate))
        self.mass_tail=mass_tail_fields(g,self.terminal_sources['m11'],R,self.full.cutoff.inner.parameters['delta'])
        self.pressure_tail=dict(P1=self.terminal_sources['m13'])
        self.terminal_zero_functions={key:a.fixed(0) for key in MOMENTS}
        self.mass_tail_zero_functions={key:a.fixed(0) if type(pair) is source.C1Function else g.zero
            for key,pair in self.mass_tail.items()}
        self.pressure_tail_zero_functions=dict(P1=a.fixed(0))
        self.definition=dict(actual_terminal_functions_retained_as_signed_integrals=True,
            exact_translated_bump_reduction_required_before_zero_certificate=True,
            five_terminal_value_and_Z1_functions_zero=True,
            velocity_moment_terminal_radius='2Rc (all new supports have already ended)',
            pressure_moment_terminal_radius='Rv (actual quiet and pulse contributions retained)',
            mass_mean_K_C1_and_QV_values_zero_after_2Rc=True,n1_pressure_C1_zero_after_Rv=True,
            original_leading_angular_and_heat_pressure_not_zeroed=True,
            tail_zero_values_follow_actual_terminal_integrals_not_enclosure_overlap=True,
            global_n1_derivatives_joins_and_mixed_bounds_complete=False,
            interior_QV_Z1_not_installed=True,temporal_recursion=False,
            full_corrected_Navier_Stokes_solution=False)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.own,self.full,self.g);self._family=copy.deepcopy(self.family)
        self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES)
                or any(receipt.get(k) for k in OPEN) or receipt['source_family']!=self.family
                or raw['terminal_identity_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]
                or not receipt['checks']['actual_integrals']['passed']
                or not receipt['checks']['transported_partition']['passed']
                or not receipt['checks']['selected_inverse']['passed']
                or not receipt['checks']['tail_recovery']['passed']):
                raise ValueError('Checked true terminal integral/zero identities required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed terminal proof source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed terminal dependency '+name)
        self.hashes[name]=digest

    def _snapshot(self):
        return parent.encode(dict(patch_endpoint=self.patch_endpoint,terminal_sources=self.terminal_sources,
            reduced_patch=self.reduced_patch,reduced_terminal=self.reduced_terminal,tail_origins=self.tail_origins,
            mass_coordinate=self.mass_coordinate,pressure_coordinate=self.pressure_coordinate,
            mass_tail=self.mass_tail,pressure_tail=self.pressure_tail,terminal_zero_functions=self.terminal_zero_functions,
            mass_tail_zero_functions=self.mass_tail_zero_functions,pressure_tail_zero_functions=self.pressure_tail_zero_functions,
            definition=self.definition))

    def assert_graph(self):
        self.own.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot
            or self.hashes!=self._hashes or (self.own,self.full,self.g)!=self._owners
            or self.family!=self._family):raise ValueError('Actual terminal or zero-tail bindings changed')

    def tail_at_ratio(self,kind,ratio):
        """Actual tail function handles, backed by a checked zero identity."""
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Tail zeros require the completed integral proof')
        if kind not in ('mass','pressure'):raise ValueError('Choose the mass or pressure tail')
        if not isinstance(ratio,(int,str,Fraction)) or isinstance(ratio,bool):
            raise TypeError('An exact rational radius ratio is required')
        q=Fraction(ratio)
        if q<1:raise ValueError('Tail radius ratio must be at least one')
        g=self.g;a=parent.C1Algebra(g);x=self.mass_coordinate if kind=='mass' else self.pressure_coordinate
        value=g.unary('log',g.constant(q));functions=self.mass_tail if kind=='mass' else self.pressure_tail
        zeros=self.mass_tail_zero_functions if kind=='mass' else self.pressure_tail_zero_functions
        fields={key:a.at(pair,x,value) if type(pair) is source.C1Function else
            g.node('function_substitution',expression=pair.node,variable=x.node,value=value.node,
                Z_independent_substitution=True) for key,pair in functions.items()}
        return dict(kind=kind,exact_radius_ratio=str(q),origin='2Rc' if kind=='mass' else 'Rv',
            actual_source_functions=fields,certified_equal_zero_functions=zeros,
            five_actual_terminal_integral_identities_consumed=True,interior_QV_Z1_claimed=False,
            leading_angular_and_heat_pressure_unchanged=True)


def run(own_owner):
    began=time.monotonic();owner=CurrentOriginalN1TerminalIntegralIdentities(own_owner,require_checked=False)
    data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,terminal_identity_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_N1_TERMINAL_INTEGRAL_AND_TAIL_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,data


if __name__=='__main__':raise SystemExit('run(own_owner) consumes checked retained integral fields')

"""True analytic core C5 and lower-collar pressure C3 on original geometry.

The core derivatives belong to the already admitted unique analytic nonlinear
solution. No new solution, model, interval endpoint or numerical center is
chosen. Existing source bodies and lower-order handles remain unchanged.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions as lower
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions as higher

source,parent,controls=lower.source,lower.parent,lower.controls
HERE,PREFIX,sha=lower.HERE,lower.PREFIX,lower.sha
NAME=PREFIX+'current_original_n1_analytic_core_lower_pressure_C3_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_analytic_core_lower_pressure_C3_functions_check.json'
GATES=('current_original_n1_actual_analytic_core_required_C5_functions_installed',
    'current_original_n1_actual_lower_U_mean_C5_Q_C4_functions_installed',
    'current_original_n1_true_lower_pressure_C3_functions_installed')
OPEN=tuple(dict.fromkeys(higher.OPEN+('current_original_n1_complete_preRp_pressure_C3_functions_installed',)))


def component_definition(row):
    out=copy.deepcopy(row);out.pop('ordinary_slow_Z_derivative_order')
    return json.dumps(out,sort_keys=True,separators=(',',':'))


class HigherCoreZDerivative(higher.HigherZDerivative):
    """The same exact analytic component, with an admitted-source guard."""
    COMPONENT='current_original_unique_analytic_core_component'
    def __init__(self,lower_owner):
        super().__init__(lower_owner.g);self.core=lower_owner.core
        self.allowed={component_definition(self.g.nodes[jet.value.node])
            for jets in lower_owner.core_radial_rows.values() for jet in jets}
        self._allowed=set(self.allowed)
    def operands(self,i):
        if self.g.nodes[i]['operation']==self.COMPONENT:return []
        return super().operands(i)
    def rule(self,i):
        row=self.g.nodes[i]
        if row['operation']==self.COMPONENT:
            if component_definition(row) not in self.allowed or self.allowed!=self._allowed:
                raise ValueError('Higher derivative must belong to the same admitted analytic core component')
            if (row['source_family']!=self.core.family or not row['unique_analytic_Banach_solution']
                or not row['actual_twenty_term_operator'] or not row['admitted_solution_not_a_model_or_enclosure_center']):
                raise ValueError('Exact original analytic core admission required')
            attrs=copy.deepcopy(row);operation=attrs.pop('operation')
            attrs['ordinary_slow_Z_derivative_order']+=1
            return self.g.node(operation,**attrs)
        return super().rule(i)


class CurrentOriginalN1AnalyticCoreLowerPressureC3Functions:
    def __init__(self,lower_owner,require_checked=True):
        if type(lower_owner) is not lower.CurrentOriginalN1LowerCollarFunctions or not lower_owner.acceptance_loaded:
            raise TypeError('The checked actual lower-collar source owner is required')
        lower_owner.assert_graph();self.lower=lower_owner;self.core=lower_owner.core
        self.g=g=lower_owner.g;self.family=copy.deepcopy(lower_owner.family);self.hashes=dict(lower_owner.hashes)
        for name in (lower.NAME,lower.RECEIPT,Path(higher.__file__).name,higher.NAME,higher.RECEIPT):self._pin(name,sha(name))
        for module in (lower,higher):
            receipt=parent.linear.read(module.RECEIPT)
            if not receipt['all_passed'] or receipt['source_family']!=self.family:
                raise ValueError('Same-family accepted lower source and higher calculus required')
            for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
        self.derivative=d=HigherCoreZDerivative(lower_owner)
        # The actual anchored gauge has a moving axial endpoint. Its already
        # accepted FTC rows seed the derivative; fixed-endpoint integration
        # is never incorrectly applied to that gauge.
        names=('Phi','Psi','mean_Psi','U0','gradient','chi','G','logF0','F0')
        for name in names:d.seed_rows(self.core.functions[name])
        for jets in lower_owner.core_radial_rows.values():
            for jet in jets:d.seed_rows(jet)
        self.core_functions={name:d.extend(self.core.functions[name],5) for name in names}
        for jet in self.core_functions.values():d.seed_rows(jet)
        self.core_radial_rows={name:[d.extend(jet,5) for jet in jets]
            for name,jets in lower_owner.core_radial_rows.items()}
        for jets in self.core_radial_rows.values():
            for jet in jets:d.seed_rows(jet)
        f=lower_owner.fields
        self.fields={key:d.extend(f[key],5) for key in ('U','U_y','m','m_y','Utheta','ratio','logratio')}
        self.fields.update({key:d.extend(f[key],4) for key in ('Q','Q_y','Q_yy')})
        self.radial_fields={key:[d.extend(jet,5) for jet in f[key]]
            for key in ('Phi_rho_rows','Psi_rho_rows','mean_Psi_rho_rows','F_rho_rows','U_rho_rows','m_rho_rows')}
        self.radial_fields['Q_rho_rows']=[d.extend(jet,4) for jet in f['Q_rho_rows']]
        self.N1p=d.extend(lower_owner.pressure_bindings['N1p'],3)
        self.signed_bodies={key:d.extend(pair,3)
            for key,pair in lower_owner.pressure_bindings['signed_physical_integrands'].items()}
        self.pressure={key:d.extend(pair,3) for key,pair in lower_owner.functions.items()}
        self.definition=dict(same_unique_analytic_nonlinear_core_components=True,
            actual_twenty_term_source_and_Green_inverse_not_changed=True,
            actual_source_component_radial_orders=[0,1,2],ordinary_Z_rows_through=5,
            anchored_gauge_higher_rows_from_true_FTC_gradient=True,
            lower_ratio_Psi_and_own_mean_histories_not_reset=True,
            original_binary_sc_and_nonzero_hb_retained=True,
            physical_U_m_Q_no_extra_Pstar=True,leading_Q_weight='1-delta',
            actual_direct_N1p_radial_viscous_terms_retained=True,
            original_m15_minus_applied_once=True,finite_original_endpoints_Z_independent=True,
            source_derivative_functions_not_numeric_enclosure_values=True,
            numerical_core_Z4_Z5_and_lower_pressure_C3_containment=False,
            actual_P0_higher_than_C3_not_installed=True,
            complete_other_pressure_and_five_defect_C3=False,
            global_n1_and_required_mixed_norms_complete=False,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.lower,self.core,self.g);self._family=copy.deepcopy(self.family)
        self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            checker_name=Path(__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES)
                or any(receipt.get(k) for k in OPEN) or receipt['source_family']!=self.family
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or archive!=sha(NAME)
                or checker!=sha(checker_name) or raw['core_lower_C5_C3_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual analytic core/lower pressure C3 functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed higher core/collar source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed higher collar input '+name)
        self.hashes[name]=digest

    def _snapshot(self):
        return higher.encode(dict(core_functions=self.core_functions,core_radial_rows=self.core_radial_rows,
            fields=self.fields,radial_fields=self.radial_fields,N1p=self.N1p,
            signed_bodies=self.signed_bodies,pressure=self.pressure,definition=self.definition))

    def assert_graph(self):
        self.lower.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or (self.lower,self.core,self.g)!=self._owners or self.family!=self._family
            or self.derivative.allowed!=self.derivative._allowed):
            raise ValueError('Actual higher core/lower function source changed')


def run(lower_owner):
    began=time.monotonic();owner=CurrentOriginalN1AnalyticCoreLowerPressureC3Functions(lower_owner,False)
    data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,core_lower_C5_C3_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_ANALYTIC_CORE_C5_AND_LOWER_PRESSURE_C3_DEFINED',len(owner.g.nodes),flush=True)
    return owner,data


if __name__=='__main__':raise SystemExit('run(lower_owner) reuses the checked actual core/collar graph')

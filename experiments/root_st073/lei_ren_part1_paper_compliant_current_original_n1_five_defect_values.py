"""Complete pre-repair n1 moment values from own cutoff plus actual tails.

Velocity moments retain their signed factored physical units. Pressure
moments add the original inner cutoff terms to both external integrals in
one fixed positive unit per moment. No enclosed function is selected as an
affine midpoint; repair and complete pressure derivative accuracy stay open.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_radial_cutoff_extension as cut
import lei_ren_part1_paper_compliant_current_original_n1_external_pressure_values as external

HERE, PREFIX, sha, ep = external.HERE, external.PREFIX, external.sha, external.ep
record, interval, source_precision = external.record, external.interval, external.source_precision
NAME = PREFIX+'current_original_n1_five_defect_values.json.gz'
RECEIPT = PREFIX+'current_original_n1_five_defect_values_check.json'
GATE = 'current_original_n1_full_five_defect_value_enclosures_installed'
OPEN = tuple(dict.fromkeys(cut.OPEN+external.OPEN))
CELLS = (('-1','-.5'),('-.5','0'),('0','.5'),('.5','1'))
VELOCITY = ('m11','m12','m14')
PRESSURE = ('m13','m15')


def pressure_unit_log(c, scales, power):
    """Log(Phi*A*Rv^p), with shared inverse-mu powers combined exactly."""
    if power not in (1,2):raise ValueError('Pressure moment radius power required')
    lp,lu,lc=[interval(c,scales[name]) for name in ('logP','logUp','logRc')]
    lm=-4*lp-c.ln(1000)
    phi=c.mpf('10.99')*c.sqrt(2*c.pi)*(c.exp(5)-c.exp(3))
    two_thirds=c.mpf(2)/3
    return (c.ln(phi)+lp+lu+power*lc-13-2*power
            +(13*power-c.mpf('7.5'))*c.exp(-lm)
            -3*c.exp(-two_thirds*(c.ln(4)+lm))
            -(60*power+c.mpf(4)/3)*lm+two_thirds*c.ln(4)-c.ln(6)/2)


def signed_tiny_term(c, term, unit):
    """Keep the actual sign, log source and a directed nonzero magnitude."""
    coefficient=interval(c,term['ordinary_Z_derivatives'][0])
    source_log=interval(c,term['combined_source_log'])
    magnitude=external.absolute_upper(c,coefficient)
    if ep(magnitude)[1]==0:
        return dict(source_exponents=term['source_exponents'],source_term=term,
            exact_coefficient_zero=True,normalized_value=c.mpf(0),absolute_upper=c.mpf(0))
    logarithm=source_log-unit+c.ln(magnitude)
    if ep(logarithm)[1]>=-1000:
        raise ValueError('Original inner cut contribution not below pressure-unit budget')
    epsilon=c.exp(-1000);lo,hi=ep(coefficient)
    value=c.mpf([0,ep(epsilon)[1]]) if lo>=0 else (
        c.mpf([-ep(epsilon)[1],0]) if hi<=0 else c.mpf([-ep(epsilon)[1],ep(epsilon)[1]]))
    return dict(source_exponents=term['source_exponents'],source_term=term,
        exact_coefficient_zero=False,normalized_log_absolute_upper=logarithm,
        normalized_value=value,absolute_upper=epsilon,
        actual_source_coefficient_sign_retained=True,nonzero_bound_not_selected_zero=True)


def inner_pressure_values(c, view, scales):
    result={}
    for key,power in (('m13',1),('m15',2)):
        actual=view['physical_cutoff_moment_values_at_Rin'][key]
        if actual['ordinary_Z_derivative_order']<1:
            raise ValueError('Original inner cut signed derivative records required')
        unit=pressure_unit_log(c,scales,power)
        terms=[signed_tiny_term(c,term,unit) for term in actual['terms']]
        value=sum((term['normalized_value'] for term in terms),c.mpf(0))
        bound=sum((term['absolute_upper'] for term in terms),c.mpf(0))
        result[key]=dict(unit_log=unit,original_signed_Z1_sectors=actual,
            normalized_signed_value=value,absolute_upper=bound,terms=terms,
            ordinary_inner_Z1_source_preserved=True,full_pressure_derivative_accuracy_not_claimed=True)
    return result


def checked_cut_views(cutoff, views=None):
    """Use genuine retained source packets on each entire original Z cell."""
    cutoff.assert_graph();result=[]
    for ends in CELLS:
        candidates=[entry for entry in cutoff._packets.values()
            if tuple(ep(entry[1]['Z']))==tuple(ep(cutoff.c.mpf(list(ends))))
            and entry[1]['query'].get('chart')=='first']
        # A retained full inner-exit packet is required, not a core/keep query.
        candidates=[entry for entry in candidates
            if entry[1]['radial_cutoff']['value']==(0,0)
            and entry[1]['query'].get('selected_fraction') in ('1',1,None)]
        if views is None or not candidates:
            packet=cutoff.evaluate(list(ends),'1','first',cells=4)
            candidates=[cutoff._packets[id(packet)]]
        entry=candidates[-1];view=cutoff.report(entry[0])
        if views is not None and cut.op.snapshot(views[ends])!=cut.op.snapshot(view):
            raise ValueError('Retained whole-Z cutoff view changed')
        if (view['source_family']!=cutoff.family or
            not view['mean_V_pressure_recovered_not_multiplied_by_cutoff'] or
            not view['full_velocity_support_fixed_before_original_inner_exit'] or
            not view['reconstructed_V_satisfies_own_n1_radial_recovery']):
            raise ValueError('Own original cutoff recovery required')
        result.append(dict(ends=list(ends),view=cut.op.snapshot(view)))
    cutoff.assert_graph()
    return result


def compose_pressure_models(c, parents, source_cells):
    result=[]
    for index, parent in enumerate(parents):
        cell=source_cells[index//128]
        rows={}
        for key in PRESSURE:
            old=parent['fields'][key];extra=interval(c,cell['inner_pressure'][key]['absolute_upper'])
            row=dict(old)
            row['inherited_external_source_error']=interval(c,old['full_nonzero_source_error'])
            row['inner_cut_nonzero_absolute_error']=extra
            row['full_nonzero_source_error']=row['inherited_external_source_error']+extra
            row['absolute_error']=interval(c,old['absolute_error'])+extra
            row['absolute_accuracy_certified']=ep(row['absolute_error'])[1]<c.mpf('.001')
            row['affine_coefficients_describe_external_part_only']=True
            row['inner_cut_sign_enclosure']=interval(c,cell['inner_pressure'][key]['normalized_signed_value'])
            rows[key]=row
        model=dict(parent);model['fields']=rows
        model['inner_cut_source_cell']=cell['ends'];model['all_inner_and_external_pressure_parts_included']=True
        result.append(model)
    if len(result)!=512:raise ValueError('All original external pressure models required')
    return result


class CurrentOriginalN1FiveDefectValues:
    @source_precision
    def __init__(self,cutoff,tails,views=None,require_checked=True):
        if type(cutoff) is not cut.CurrentOriginalN1RadialCutoffExtension or not cutoff.acceptance_loaded:
            raise ValueError('Checked actual n1 cutoff owner required')
        if type(tails) is not external.CurrentOriginalN1ExternalPressureValues or not tails.acceptance_loaded:
            raise ValueError('Checked original complete external pressure required')
        cutoff.assert_graph();tails.assert_graph()
        if cutoff.family!=tails.family:raise ValueError('Common original cut and tail family required')
        self.cutoff,self.tails,self.c=cutoff,tails,tails.c
        self.family=copy.deepcopy(cutoff.family);self.hashes=dict(tails.hashes)
        for name,digest in cutoff.hashes.items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting original source '+name)
            self.hashes[name]=digest
        for module in (cut,external):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.bindings={name:external.upper.pulse.end.coefficients.function_binding(fn) for name,fn in (
            ('fixed_pressure_units',pressure_unit_log),('signed_cut_sector_value',signed_tiny_term),
            ('actual_inner_pressure_values',inner_pressure_values),('retained_cut_packets',checked_cut_views),
            ('complete_pressure_models',compose_pressure_models),('original_cut_own_recovery',cut.recover_velocity),
            ('original_cut_moment_changes',cut.change_densities),('original_cumulative_endpoints',cut.parent.endpoint_cumulants))}
        self.acceptance_loaded=False
        if require_checked:
            raw=json.loads(gzip.decompress((HERE/NAME).read_bytes()))
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Checked complete five-defect values required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed five-defect source '+name)
            self.source_cells=raw['source_cells'];self.models=raw['pressure_models']
            self.hashes.update(receipt['input_hashes']);self.acceptance_loaded=True
        else:
            self.source_cells=checked_cut_views(cutoff,views)
            for cell in self.source_cells:
                cell['inner_pressure']=record(inner_pressure_values(self.c,cell['view'],tails.scale_data))
            self.models=record(compose_pressure_models(self.c,tails.models,self.source_cells))
        self._packets={};self._state=self.state()

    def state(self):return external.snapshot((self.family,self.hashes,self.source_cells,self.models,self.bindings))

    @source_precision
    def assert_graph(self):
        self.cutoff.assert_graph();self.tails.assert_graph()
        if self.state()!=self._state or sha(Path(__file__).name)!=self.hashes[Path(__file__).name]:
            raise ValueError('Complete five-defect source changed')

    def _make_view(self,Z):
        c=self.c;z=c.mpf(Z);tail=self.tails._make_view(Z);active=[];pressure={key:[] for key in PRESSURE}
        for cell in self.source_cells:
            left,right=[c.mpf(s) for s in cell['ends']]
            if ep(z)[0]>ep(right)[1] or ep(z)[1]<ep(left)[0]:continue
            active.append(cell)
            for key in PRESSURE:
                pressure[key].append(interval(c,cell['inner_pressure'][key]['normalized_signed_value']))
        if not active:raise ValueError('Original complete axial domain required')
        inner={key:c.mpf([min(ep(v)[0] for v in values),max(ep(v)[1] for v in values)]) for key,values in pressure.items()}
        full={key:tail['full_pressure_values'][key]+inner[key] for key in PRESSURE}
        source_velocity={key:[dict(exact_Z_cell=cell['ends'],
            signed_physical_sectors=cell['view']['physical_cutoff_moment_values_at_Rin'][key],
            signed_source_budgets=cell['view']['full_velocity_support_moment_budgets'][key]) for cell in active] for key in VELOCITY}
        return dict(Z=z,source_family=self.family,moment_order=['m11','m12','m13','m14','m15'],
            full_velocity_moment_source_values=source_velocity,full_pressure_moment_values=full,
            inner_cut_pressure_values=inner,external_pressure_values=tail['full_pressure_values'],
            original_pressure_units=tail['original_physical_units'],active_source_cells=[cell['ends'] for cell in active],
            full_radial_support_included=True,original_own_M_V_P_and_independent_P0_preserved=True,
            original_m15_external_minus_sign_applied_once=True,pre_repair_defects_not_terminal_zero=True,
            velocity_values_remain_factored_not_zeroed=True,full_pressure_derivative_accuracy_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();view=self._make_view(Z);packet=object()
        self._packets[id(packet)]=(packet,view,external.snapshot(view));return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if entry is None or entry[0] is not packet or entry[2]!=external.snapshot(entry[1]):
            raise ValueError('Live full five-defect packet required')
        return external.upper.pulse.end.copy_view(self.c,entry[1])


@source_precision
def run(cutoff,tails,views=None):
    began=time.monotonic();owner=CurrentOriginalN1FiveDefectValues(cutoff,tails,views,False)
    owner.assert_graph()
    samples={z:owner._make_view(z) for z in ('-1','-.639','0','.211','.371','.639','1')}
    owner.assert_graph()
    raw=dict(source_family=owner.family,source_cells=record(owner.source_cells),pressure_models=record(owner.models),
        samples=record(samples),source_bindings=owner.bindings,input_hashes=owner.hashes,
        full_radial_support_included=True,three_velocity_moments_remain_signed_physical_sectors=True,
        both_pressure_cut_and_external_parts_included=True,inner_cut_first_derivatives_not_full_pressure_accuracy=True,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_FIVE_DEFECT_VALUES_COMPUTED',len(owner.source_cells),len(owner.models),flush=True)
    return owner,raw

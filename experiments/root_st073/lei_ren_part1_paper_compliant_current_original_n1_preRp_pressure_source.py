"""Original r_minus-to-Rp pressure bound and r_minus-to-infinity models.

The 17 active charts use genuine full-Z physical value caps, exact reduced
pressure identities and finite endpoints. Original compact/quiet signed
pieces are reused once. Common Rc powers cancel before scale arithmetic.
No bounded omitted contribution is assigned a signed zero value.
"""
import copy
import gzip
import importlib
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_source as pulse
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields as physical
import lei_ren_part1_paper_compliant_current_original_n1_preRc_pressure_endpoint as endpoint
import lei_ren_part1_paper_compliant_current_original_n1_repair_band_pressure_source as band
import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_source as quiet
import lei_ren_part1_paper_compliant_current_original_n1_radial_free_pressure_integrals as reduced
import lei_ren_part1_paper_compliant_current_original_n1_pressure_forcing_support as support

# The endpoint adapter is stateless. A fresh producer load binds its current
# code even when a long computation keeps the checked parent owners alive.
endpoint=importlib.reload(endpoint)

HERE,PREFIX,sha,ep=pulse.HERE,pulse.PREFIX,pulse.sha,pulse.ep
source_precision=pulse.source_precision
rg=physical.rg
NAME=PREFIX+'current_original_n1_preRp_pressure_source.json.gz'
RECEIPT=PREFIX+'current_original_n1_preRp_pressure_source_check.json'
GATE='current_original_n1_rminus_to_infinity_pressure_values_installed'
OPEN=pulse.OPEN


def record(value):
    if isinstance(value,rg.LogUpper):return record(value.record())
    if isinstance(value,dict):return {key:record(v) for key,v in value.items()}
    if isinstance(value,(list,tuple)):return [record(v) for v in value]
    return pulse.record(value)


def interval(c,value):
    if hasattr(value,'_mpi_'):return c.make_mpf(value._mpi_)
    return c.make_mpf((tuple(value['lower_exact_mpf_tuple']),tuple(value['upper_exact_mpf_tuple'])))


def log_cap(c,value):
    return rg.LogUpper(c) if value.get('exact_zero') else rg.LogUpper(c,interval(c,value['log_absolute_upper']))


def maximum(c,values):
    logs=[ep(v.log)[1] for v in values if v.log is not None]
    return rg.LogUpper(c,c.mpf(max(logs))) if logs else rg.LogUpper(c)


def preRc_coefficients(c,views,endpoints):
    """Absolute reduced-S/viscous caps after factoring Rc and Rc squared."""
    bd=rg.Bounds(c)
    U=maximum(c,[log_cap(c,v['whole_cell_physical_U_Z3_magnitude'][0]) for v in views])
    UZ=maximum(c,[log_cap(c,v['whole_cell_physical_U_Z3_magnitude'][1]) for v in views])
    q=maximum(c,[log_cap(c,v['whole_cell_physical_q_Z2_magnitude'][0]) for v in views])
    qZ=maximum(c,[log_cap(c,v['whole_cell_physical_q_Z2_magnitude'][1]) for v in views])
    common=bd.sum(qZ,bd.scale(U*qZ,2),bd.scale(U*UZ,8))
    S0=bd.sum(bd.scale(q,4),q*q,bd.scale(U*U,20),common)
    S1=bd.sum(bd.scale(q,6),bd.scale(q*q,c.mpf('1.5')),bd.scale(U*U,24),common)
    lower,upper=endpoints['lower']['records'],endpoints['upper']['records']
    visc0=bd.sum(log_cap(c,lower['q_plus_q_y']),log_cap(c,upper['q_plus_q_y']))
    visc1=bd.sum(log_cap(c,lower['q_y']),log_cap(c,upper['q_y']))
    # Rc>=1, so an unweighted finite endpoint is bounded by Rc*cap,
    # and |R*q_y|<=Rc*cap<=Rc^2*cap. Neither endpoint is dropped.
    return dict(fields=dict(U=U,U_Z=UZ,q=q,q_Z=qZ),
        reduced_S=dict(unweighted=S0,weighted=S1),viscous=dict(m13=visc0,m15=visc1),
        m13=bd.sum(bd.scale(S0,c.mpf('.5')),visc0),
        m15=bd.sum(bd.scale(S1,c.mpf('.5')),visc1),
        bounds_scope='every original active R in [r_minus,Rc]',
        physical_measures='integral dR<=Rc; integral R*dR<=Rc^2',
        L_lower='.5',Z_absolute_upper=1,delta_upper='.5',
        value_rows_only=True,no_q_ZZ_needed=True)


def Rc_unit_log(c,owner,power):
    """Log(Rc^p/(Phi*Aphys*Rv^p)) with shared Rc cancelled exactly."""
    if power not in (1,2):raise ValueError('Original pressure radius power required')
    end=owner.end
    logE=pulse.rebox(c,owner.coefficients.selected.amplitude.log_end_scale)
    return (-end.logP-end.logUp-logE+13-c.ln(owner.reference)
        +power*(2+60*c.ln(owner.mu))+(c.mpf('6.5')-13*power)/owner.mu)


def checked_archive(module,gate,family):
    data=(HERE/module.NAME).read_bytes()
    raw=json.loads(gzip.decompress(data) if module.NAME.endswith('.gz') else data)
    receipt=json.loads((HERE/module.RECEIPT).read_bytes())
    if not receipt.get('all_passed') or not receipt.get(gate) or receipt['source_family']!=family or raw['source_family']!=family:
        raise ValueError('Checked same-family radial pressure archive required '+module.NAME)
    hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if sha(name)!=digest:raise ValueError('Changed radial pressure source '+name)
    for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):hashes[name]=sha(name)
    return raw,hashes


def checked_views(fields,views=None):
    """Reuse live checked packets; preserve all four Z cells and 17 charts."""
    fields.assert_graph();out={}
    for ends in pulse.CELLS:
        for chart in fields.charts:
            key=(ends,chart)
            entries=[entry for entry in fields._packets.values()
                if tuple(entry[1]['view']['exact_Z_cell'])==ends and entry[1]['view']['chart']==chart]
            if not entries:
                packet=fields.range_chart(ends,chart)
                entries=[fields._packets[id(packet)]]
            entry=entries[-1];state=entry[1]
            actual=state['view']
            payload=(actual,physical._bound_record(state['bounds']),physical.physical._z1_encode(state['functions']))
            if physical.fingerprint(payload)!=entry[2]:raise ValueError('Changed original chart packet')
            view={k:v for k,v in actual.items() if not k.startswith('signed_')}
            if views is not None and record(views[key])!=record(view):raise ValueError('Cached chart view differs')
            if not (view['source_family']==fields.family and view['whole_cell_range_not_endpoint_value']
                    and view['exact_original_source_cell_union'] and view['physical_common_scale_Pstar_applied_once']):
                raise ValueError('Original full-chart physical source required')
            partition=fields.partitions[chart]
            if record(view['exact_native_cell'])!=record([partition[0],partition[-1]]):
                raise ValueError('Complete original chart union required')
            out[key]=view
    fields.assert_graph()
    return out


def build_data(fields,whole,views):
    c=whole.c
    geometry=endpoint.geometry_proof(fields,whole)
    if ep(whole.end.logRc)[0]<=0:raise ValueError('Actual original Rc>=1 required')
    geometry['original_Rc_ge_one']=True
    geometry['actual_logRc_positive']=whole.end.logRc
    band_raw,band_hashes=checked_archive(band,band.GATE,whole.family)
    quiet_raw,quiet_hashes=checked_archive(quiet,quiet.GATE,whole.family)
    zero_raw,zero_hashes=checked_archive(support,support.GATE,whole.family)
    if set(zero_raw['future_charts'])!=set(support.FUTURE):raise ValueError('Complete original beyond-Rv support required')
    if zero_raw['endpoint']['definition']!='Rv=exp(original_logRp+13/original_mu)':
        raise ValueError('Same original Rv support endpoint required')
    if not zero_raw['endpoint']['angular_and_pressure_not_zeroed']:raise ValueError('Meridional support scope required')
    for view in zero_raw['future_charts'].values():
        if view['exact_zero_unweighted_and_R_weighted_forcing_integrals']!=[0,0]:
            raise ValueError('Original beyond-Rv forcing integral must be exactly zero')
    cells=[]
    for ends in pulse.CELLS:
        endpoints=endpoint.endpoint_caps(fields,ends)
        caps=preRc_coefficients(c,[views[(ends,chart)] for chart in fields.charts],endpoints)
        original_band=next(item['view'] for item in band_raw['samples'].values() if tuple(item['exact_Z_cell'])==ends)
        original_quiet=quiet_raw['actual_quiet_source']
        rows={}
        for key,power in (('m13',1),('m15',2)):
            unit=Rc_unit_log(c,whole,power)
            band_log=interval(c,original_band['signed_pressure_'+key]['log_absolute_coefficient'])
            qkey='log_coefficient' if key=='m13' else 'log_absolute_coefficient'
            quiet_log=interval(c,original_quiet[key][qkey])
            logs=dict(active=caps[key].log+unit,compact_band=band_log+unit,quiet=quiet_log+unit)
            if any(ep(value)[1]>=-1000 for value in logs.values()):
                raise ValueError('Actual pre-Rp pressure cap not below source error budget')
            rows[key]=dict(normalized_log_caps=logs,
                each_nonzero_absolute_upper=c.exp(-1000),
                total_nonzero_absolute_upper=c.exp(-1000)*3,
                no_signed_zero_assigned=True,all_pieces_counted_once=True)
        cells.append(dict(ends=ends,endpoints=record(endpoints),active=record(caps),fields=rows))
    return dict(source_family=whole.family,geometry=record(geometry),cells=record(cells),
        accepted_beyond_Rv_support=dict(endpoint=zero_raw['endpoint'],proof=zero_raw['proof'],
            future_charts=zero_raw['future_charts'],acceptance_receipt=support.RECEIPT),
        actual_chart_views=[dict(ends=ends,chart=chart,view=record(view)) for (ends,chart),view in views.items()],
        scope=dict(lower='original r_minus',upper='original Rp',
            active='17 charts r_minus to Rc',compact='Rc to 2Rc',quiet='2Rc to Rp',
            original_shared_Rc_cancelled_before_boxing=True,
            bounded_preRp_is_nonzero_remainder_not_signed_center=True),
        extra_hashes={**band_hashes,**quiet_hashes,**zero_hashes})


class CurrentOriginalN1PreRpPressureSource:
    @source_precision
    def __init__(self,fields,whole,views=None,require_checked=True):
        if type(fields) is not physical.CurrentOriginalN1AllChartPhysicalPressureFields or not fields.acceptance_loaded:
            raise ValueError('Checked all-chart physical owner required')
        if type(whole) is not pulse.CurrentOriginalN1PulsePressureSource or not whole.acceptance_loaded:
            raise ValueError('Checked entire outer-pulse owner required')
        if fields.family!=whole.family:raise ValueError('Same original family required')
        fields.assert_graph();whole.assert_graph()
        self.fields,self.pulse,self.c,self.family=fields,whole,whole.c,copy.deepcopy(whole.family)
        if require_checked:
            raw=json.loads(gzip.decompress((HERE/NAME).read_bytes()))
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Accepted rminus-to-infinity pressure receipt required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed combined pre-Rp input '+name)
            self.data=raw['source_data']
            self.hashes=dict(receipt['input_hashes'])
        else:
            self.views=checked_views(fields,views)
            self.data=build_data(fields,whole,self.views)
            self.hashes=dict(fields.hashes)
            for inputs in (whole.hashes,self.data['extra_hashes']):
                for name,digest in inputs.items():
                    if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Conflicting source '+name)
                    self.hashes[name]=digest
            for module in (physical,endpoint,pulse,reduced):
                self.hashes[Path(module.__file__).name]=sha(Path(module.__file__).name)
            for module in (physical,pulse):
                self.hashes[module.NAME]=sha(module.NAME);self.hashes[module.RECEIPT]=sha(module.RECEIPT)
            self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.preRp_error={key:self.c.exp(-1000)*3 for key in ('m13','m15')}
        self.acceptance_loaded=require_checked
        self._state=record((self.data,self.family,self.hashes,self.preRp_error))
        self._packets={}

    @source_precision
    def assert_graph(self):
        self.fields.assert_graph();self.pulse.assert_graph()
        if self._state!=record((self.data,self.family,self.hashes,self.preRp_error)):
            raise ValueError('Pre-Rp pressure data/source changed')

    def _make_view(self,Z):
        previous=self.pulse._make_view(Z)
        c=self.c;values={key:previous['full_pressure_values'][key]+c.mpf([-ep(error)[1],ep(error)[1]])
            for key,error in self.preRp_error.items()}
        view=dict(Z=c.mpf(Z),source_family=self.family,radial_scope='r_minus to infinity',
            full_pressure_values=values,preRp_nonzero_error=self.preRp_error,
            original_physical_units=self.pulse.units,
            outer_pulse_view=previous,lower_Rin_to_rminus_not_included=True,
            preRp_bound_is_not_signed_zero=True,pressure_derivative_accuracy_not_claimed=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        return view

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();view=self._make_view(Z)
        packet=object();self._packets[id(packet)]=(packet,view,record(view))
        return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if entry is None or entry[0] is not packet or entry[2]!=record(entry[1]):raise ValueError('Live pre-Rp packet required')
        return pulse.end.copy_view(self.c,entry[1])

    @source_precision
    def affine_models(self):
        self.assert_graph()
        previous=json.loads(gzip.decompress((HERE/pulse.NAME).read_bytes()))
        models=copy.deepcopy(previous['models'])
        for model in models:
            for key,field in model['fields'].items():
                field['preRp_nonzero_source_error']=self.preRp_error[key]
                field['full_nonzero_source_error']=interval(self.c,field['full_nonzero_source_error'])+self.preRp_error[key]
                field['absolute_error']=interval(self.c,field['absolute_error'])+self.preRp_error[key]
                field['absolute_accuracy_certified']=ep(field['absolute_error'])[1]<ep(self.c.mpf('.001'))[0]
            model['radial_scope']='r_minus to infinity'
        self.assert_graph()
        return models


@source_precision
def run(fields,whole,views=None):
    began=time.monotonic();owner=CurrentOriginalN1PreRpPressureSource(fields,whole,views,False)
    models=owner.affine_models()
    samples={z:owner._make_view(z) for z in ('-1','-.639','0','.211','.371','.639','1')}
    owner.assert_graph()
    raw=dict(source_family=owner.family,source_data=owner.data,models=record(models),samples=record(samples),
        original_physical_units=owner.pulse.units,
        source_bindings={name:pulse.end.coefficients.function_binding(fn) for name,fn in (
            ('active_reduced_pressure_caps',preRc_coefficients),('shared_Rc_unit_ratio',Rc_unit_log),
            ('complete_preRp_source_data',build_data),('checked_chart_views',checked_views),
            ('original_endpoint_caps',endpoint.endpoint_caps),('original_geometry',endpoint.geometry_proof),
            ('original_exact_beyond_Rv_support',support.source_proof),
            ('original_reduced_pressure',reduced.reduced_S_integrals))},
        input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False),
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(record(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_RMINUS_TO_INFINITY_PRESSURE_MODELS_COMPUTED',len(models),flush=True)
    return owner,raw

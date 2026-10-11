"""Actual finite background remainder relative to signed stress divergence.

Ratios use complete original physical log scales, including nonzero-Z
lambda. They are finite-sector diagnostics, not all-orders time flatness.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_stress_cone as cone
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

bg=cone.background
HERE,PREFIX,sha=cone.HERE,cone.PREFIX,cone.sha
NAME=PREFIX+'current_original_Rp_remainder_ratios.json.gz'
RECEIPT=PREFIX+'current_original_Rp_remainder_ratios_check.json'
GATE='current_original_Rp_source_correlated_finite_remainder_relative_indicators_installed'
OPEN=cone.OPEN
PAIRS=dict(theta_remainder_over_theta_divergence=('Etheta','divtheta'),
    axial_remainder_over_axial_divergence=('Ez','divz'),
    radial_remainder_over_theta_reference=('Er','divtheta'),
    radial_transport_over_theta_reference=('Er_transport','divtheta'),
    radial_axial_viscosity_over_theta_reference=('Er_axial_viscosity','divtheta'))


@dataclass(frozen=True,eq=False)
class OriginalRemainderRatios:
    chart:str


class CurrentOriginalRpRemainderRatios:
    @source_precision
    def __init__(self,before,require_checked=True):
        if type(before) is not cone.CurrentOriginalRpStressCone or not before.acceptance_loaded:
            raise ValueError('Accepted original point cone owner required')
        self.before,self.background=before,before.before
        self.product,self.arithmetic=self.background.product,self.background.arithmetic
        self.graph,self.ctx,self.family_record=before.graph,before.ctx,before.family_record
        self.pairs=dict(PAIRS)
        self.hashes=dict(before.hashes)
        for name in (cone.NAME,cone.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self._fields={}
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or any(receipt[key] for key in OPEN) or \
                    receipt['source_family']!=self.family_record or receipt['original_ratio_pairs']!={k:list(v) for k,v in PAIRS.items()}:
                raise ValueError('Original finite remainder ratio receipt/scope differs')
            for name in (Path(__file__).name,Path(__file__).stem+'_check.py',NAME):
                if receipt['input_hashes'].get(name)!=sha(name):raise ValueError('Unbound finite remainder source '+name)
            bg.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        checks=dict(accepted_live_original_point_cone=self.before.acceptance_loaded,
            accepted_live_original_full_background=self.background.acceptance_loaded,
            same_actual_graph_and_scale_arithmetic=self.graph is self.arithmetic.graph is self.product.graph,
            original_ratio_pairs_unchanged=PAIRS==self.pairs,
            exact_radial_sector_exponents=tuple(beta for beta,_ in self.background.operators['cylindrical_remainder']['radial'])
                ==(-3,-3+2*bg.DELTA))
        if not all(checks.values()):raise ValueError('Original finite remainder source differs: '+str(checks))
        return checks

    def _issue(self,label,row,delivery,reader):
        coefficient=row['common_scale_coefficient_enclosure'];sign=row['signed_log_value']['sign']
        if coefficient is None or sign not in (-1,1):raise ValueError('Strict signed finite sector required: '+label)
        scale=bg.box.pulse.radius.FunctionRef(self.graph,row['exact_reference_log_scale_function'])
        original=dict(row,component=label)
        record=dict(component=label,unit='acceleration',scale=scale,scale_polynomial=reader.polynomial(scale),
            coefficient_tuple=coefficient._mpi_,sign=sign,scale_dag=self.arithmetic._scale_dag(scale),
            root=object(),multiplier=Fraction(1),delivery=delivery,bindings=dict(reader.bindings),
            aliases=dict(reader.aliases),allowed=reader.allowed_exponentials,
            alias_snapshot=dict(reader.aliases),allowed_snapshot=frozenset(reader.allowed_exponentials),
            binding_snapshot={node:value._mpi_ for node,value in reader.bindings.items()},
            absolute_width_satisfied=row['physical_accuracy'].get('ordinary_numeric_relative_width_satisfied',False),
            original_materialized=row['ordinary_numeric_materialized'],original_row=original)
        return self.arithmetic._issue(record)

    def _ratio_dag(self,ratios):
        pending=[result[key] for result in ratios.values()
            for key in ('exact_scale_log_difference','exact_positive_scale_ratio_function','exact_ratio_identity')]
        snapshots={}
        while pending:
            node=pending.pop()
            if node in snapshots:continue
            data=self.graph.nodes[node];snapshots[node]=copy.deepcopy(data);op=data['operation']
            if op in ('sum','product'):pending.extend(data['arguments'])
            elif op in ('negative','analytic_unary'):pending.append(data['argument'])
            elif op=='positive_quotient':pending.extend((data['numerator'],data['denominator']))
            elif op=='exact_operator_integer_power':pending.append(data['base'])
            elif op=='exact_original_Rp_source_scaled_ratio_cancellation':
                pending.extend(data[key] for key in ('numerator_scale','denominator_scale','exact_log_scale_difference'))
        return snapshots

    @source_precision
    def evaluate(self,packet,relative_width_target='1/1000'):
        self.assert_graph();view=self.background.report(packet)
        field=self.background._fields[id(packet)][1]
        record=self.background.differential._require(field);delivery=record['delivery']
        target=Fraction(relative_width_target)
        if not 0<target<1:raise ValueError('Exact relative-width target required')
        source=self.product.before.source(delivery);reader=self.product.reader(delivery)
        _,_,request=self.product.amplitude.before._validate_delivery(delivery)
        parts=self.background.operators['cylindrical_remainder']['radial']
        radial={name:self.background._enclose([part],source,reader,request.forward_coordinates,target)
            for name,part in zip(('Er_transport','Er_axial_viscosity'),parts)}
        rows=dict(Etheta=view['cylindrical_remainder']['theta'],Ez=view['cylindrical_remainder']['axial'],
            Er=view['cylindrical_remainder']['radial'],divtheta=view['cylindrical_divergence']['theta'],
            divz=view['cylindrical_divergence']['axial'],**radial)
        values={name:self._issue(name,row,delivery,reader) for name,row in rows.items()}
        for name in ('divtheta','divz'):
            if self.arithmetic._require(values[name])['sign']!=1:
                raise ValueError('Actual strictly positive same-unit divergence reference required')
        ratios={name:self.arithmetic.ratio(values[a],values[b]) for name,(a,b) in PAIRS.items()}
        value=OriginalRemainderRatios(record['chart'])
        self._fields[id(value)]=(value,packet,rows,values,ratios,
            self.background._source_dag({'radial_sectors':radial}),self._ratio_dag(ratios))
        return value

    @source_precision
    def report(self,value):
        self.assert_graph();entry=self._fields.get(id(value))
        if type(value) is not OriginalRemainderRatios or entry is None or entry[0] is not value:
            raise ValueError('Live original finite remainder ratios issued by this owner required')
        packet_view=self.background.report(entry[1])
        if packet_view['chart']!=value.chart or any(self.graph.nodes[node]!=data for snapshots in entry[5:]
                for node,data in snapshots.items()):raise ValueError('Original finite remainder operators changed')
        for label,issued in entry[3].items():
            record=self.arithmetic._require(issued)
            if record['coefficient_tuple']!=entry[2][label]['common_scale_coefficient_enclosure']._mpi_ or \
                    record['scale'].node!=entry[2][label]['exact_reference_log_scale_function']:
                raise ValueError('Original signed finite remainder value changed')
        field=self.background._fields[id(entry[1])][1]
        delivery=self.background.differential._require(field)['delivery']
        _,_,request=self.product.amplitude.before._validate_delivery(delivery)
        return dict(chart=value.chart,source_family=copy.deepcopy(self.family_record),
            source_correlated_ratios=copy.deepcopy(entry[4]),radial_sectors={name:copy.deepcopy(entry[2][name])
                for name in ('Er_transport','Er_axial_viscosity')},
            original_ratio_pairs={name:list(pair) for name,pair in PAIRS.items()},
            all_references_strictly_positive_and_same_physical_unit=True,physical_unit='acceleration/acceleration',
            actual_complete_lambda_at_nonzero_Z_retained=request.Z!=0,original_viscosity=1,
            radial_divergence_is_zero_not_used_as_denominator=True,
            radial_theta_reference_is_a_scale_indicator_not_radial_balance=True,
            finite_background_sector_relative_indicator=True,all_orders_time_flatness_certified=False,
            flatness_not_certified=True,physical_absolute_accuracy_not_promoted=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))


@source_precision
def run(before,packets):
    began=time.monotonic();owner=CurrentOriginalRpRemainderRatios(before,require_checked=False)
    values={name:owner.evaluate(packet) for name,packet in packets.items()}
    views={name:owner.report(value) for name,value in values.items()}
    result=dict(source_family=owner.family_record,actual_original_finite_remainder_ratios=views,
        original_ratio_pairs={name:list(pair) for name,pair in PAIRS.items()},input_hashes=owner.hashes,
        source_assertions=owner.assert_graph(),execution_seconds=time.monotonic()-began,
        **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(bg.correlated.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,values,views

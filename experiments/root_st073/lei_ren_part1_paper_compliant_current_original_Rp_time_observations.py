"""Actual multi-time values along one original similarity-coordinate path.

Cross-time scale cancellation retains both issued coefficient enclosures.
This path is neither a measured vortex core nor a recursive correction order.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_remainder_ratios as relative
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

bg=relative.bg
HERE,PREFIX,sha=relative.HERE,relative.PREFIX,relative.sha
NAME=PREFIX+'current_original_Rp_time_observations.json.gz'
RECEIPT=PREFIX+'current_original_Rp_time_observations_check.json'
GATE='current_original_Rp_actual_multitime_similarity_path_observations_installed'
OPEN=relative.OPEN


@dataclass(frozen=True,eq=False)
class OriginalTimeSample:
    field:object
    axial:object
    background:object
    cone:object
    remainder:object


@dataclass(frozen=True,eq=False)
class OriginalTimeObservations:
    chart:str


class CurrentOriginalRpTimeObservations:
    @source_precision
    def __init__(self,before,require_checked=True):
        if type(before) is not relative.CurrentOriginalRpRemainderRatios or not before.acceptance_loaded:
            raise ValueError('Accepted live finite remainder ratio owner required')
        self.before,self.background,self.cone=before,before.background,before.before
        self.axial=self.background.before;self.differential=self.background.differential
        self.product,self.arithmetic=before.product,before.arithmetic
        self.graph,self.ctx,self.family_record=before.graph,before.ctx,before.family_record
        self.hashes=dict(before.hashes)
        for name in (relative.NAME,relative.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self._fields={}
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or any(receipt[key] for key in OPEN) or \
                    receipt['source_family']!=self.family_record:
                raise ValueError('Actual multi-time observation receipt/scope differs')
            for name in (Path(__file__).name,Path(__file__).stem+'_check.py',NAME):
                if receipt['input_hashes'].get(name)!=sha(name):raise ValueError('Unbound multi-time source '+name)
            bg.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        checks=dict(accepted_original_relative_owner=self.before.acceptance_loaded,
            accepted_original_derivative_and_axial_owners=self.differential.acceptance_loaded and self.axial.acceptance_loaded,
            accepted_original_background_and_cone_owners=self.background.acceptance_loaded and self.cone.acceptance_loaded,
            identical_live_graph=self.graph is self.arithmetic.graph is self.product.graph,
            original_guard_unchanged=bg.signed.EXP_LOG_LIMIT==1000)
        if not all(checks.values()):raise ValueError('Original time observation source differs: '+str(checks))
        return checks

    def _validate_sample(self,sample,return_views=False):
        if type(sample) is not OriginalTimeSample:raise ValueError('Actual typed original time sample required')
        entries=(self.axial._issued.get(id(sample.axial)),self.background._fields.get(id(sample.background)),
            self.cone._fields.get(id(sample.cone)),self.before._fields.get(id(sample.remainder)))
        handles=(sample.axial,sample.background,sample.cone,sample.remainder)
        if any(entry is None or entry[0] is not handle for entry,handle in zip(entries,handles)) or \
                any(entry[1] is not sample.field for entry in entries[:3]) or entries[3][1] is not sample.background:
            raise ValueError('All sample packets must belong to the same actual field delivery')
        record=self.differential._require(sample.field)
        self.axial.report(sample.axial)
        # The relative owner validates its actual background packet and field.
        # Retain these returned views within this call instead of validating
        # the same background/cone/relative packet a second time for export.
        cone_view=self.cone.report(sample.cone);relative_view=self.before.report(sample.remainder)
        if not cone_view[relative.cone.POINT_CONE]:raise ValueError('Actual strictly admissible sample cone required')
        _,_,request=self.product.amplitude.before._validate_delivery(record['delivery'])
        if request.chart!='pulse_exit' or request.Z==0:raise ValueError('Supported nonzero-Z original pulse path required')
        return (record,request,dict(cone=cone_view,remainder=relative_view)) if return_views else (record,request)

    def _merged_reader(self,deliveries,ctx=None):
        bindings,aliases={},{};allowed=None
        for delivery in deliveries:
            reader=self.product.reader(delivery)
            if allowed is None:allowed=reader.allowed_exponentials
            elif allowed!=reader.allowed_exponentials:raise ValueError('Cross-time exponential whitelists differ')
            for node,bound in reader.bindings.items():
                if node in bindings and bindings[node]._mpi_!=bound._mpi_:
                    raise ValueError('Cross-time shared source binding differs')
                bindings[node]=bound
            for node,value in reader.aliases.items():
                if node in aliases and aliases[node]!=value:raise ValueError('Cross-time original aliases conflict')
                aliases[node]=value
        return bg.correlated.CorrelatedGraphBounds(bg.correlated.locator.CancelledGraphBounds(
            self.graph,self.ctx if ctx is None else ctx,bindings,allowed),aliases)

    def _profile_row(self,source,reader,request,label,beta):
        return self.background._enclose([(beta,bg.F[label])],source,reader,
            request.forward_coordinates,Fraction(1,1000))

    def _observables(self,sample,record,request):
        g=self.graph;delivery=record['delivery'];coordinates=request.forward_coordinates
        source=self.product.before.source(delivery);local_reader=self.product.reader(delivery)
        one=self.ctx.mpf(1);result={}
        def insert(name,scale,coefficient,unit,exponent=None,row=None):
            lo,hi=bg.ends(coefficient)
            sign=1 if lo>0 else -1 if hi<0 else None
            if sign is None:raise ValueError('Strict nonzero time observable required: '+name)
            result[name]=dict(scale=scale,coefficient=coefficient,sign=sign,unit=unit,
                expected_tau_exponent=exponent,row=row)
        delta=self.product.delta
        q=lambda value:g.constant(Fraction(value))
        insert('tau',coordinates['log_tau'],one,'time',q(1))
        insert('lambda',coordinates['loglambda'],one,'length',q('1/2'))
        insert('radius',coordinates['log_r'],one,'length',q('1/2'))
        insert('abs_axial_coordinate',request.physical_log_abs_z,one,'length',
            g.mul(q('1/2'),g.sub(g.one,delta)))
        insert('axial_to_radial_coordinate_ratio',g.sub(request.physical_log_abs_z,coordinates['log_r']),
            one,'dimensionless',g.mul(q('-1/2'),delta))
        for label,beta,unit in (('Ur',-1,'velocity'),('Utheta',-1-bg.DELTA,'velocity'),
                ('Uz',-1-bg.DELTA,'velocity'),('pressure',-2-2*bg.DELTA,'pressure')):
            row=self._profile_row(source,local_reader,request,label,beta)
            exponent=self.background._exact(bg.s.sympify(beta)/2,coordinates)
            insert(label,bg.box.pulse.radius.FunctionRef(g,row['exact_reference_log_scale_function']),
                row['common_scale_coefficient_enclosure'],unit,exponent,row)
        for label,parts in dict(u=[(-1,bg.CS*bg.F['Ur']),(-1-bg.DELTA,-bg.SN*bg.F['Utheta'])],
                v=[(-1,bg.SN*bg.F['Ur']),(-1-bg.DELTA,bg.CS*bg.F['Utheta'])],
                w=[(-1-bg.DELTA,bg.F['Uz'])]).items():
            row=self.background._enclose(parts,source,local_reader,coordinates,Fraction(1,1000))
            # u/v have two original lambda powers; no single exponent is imposed.
            exponent=None if label in ('u','v') else self.background._exact((-1-bg.DELTA)/2,coordinates)
            insert(label,bg.box.pulse.radius.FunctionRef(g,row['exact_reference_log_scale_function']),
                row['common_scale_coefficient_enclosure'],'velocity',exponent,row)
        axial=self.arithmetic._require(sample.axial)
        insert('omega_z',axial['scale'],self.arithmetic._coefficient(axial),axial['unit'],
            g.mul(q('1/2'),g.sub(q(-2),delta)))
        view=self.background.report(sample.background)
        selected=(('stress_rtheta','cylindrical_stress','r_theta',-2-bg.DELTA,'stress'),
            ('stress_rz','cylindrical_stress','r_z',-2-bg.DELTA,'stress'),
            ('stress_thetatheta','cylindrical_stress','theta_theta',-2,'stress'),
            ('divtheta','cylindrical_divergence','theta',-3-bg.DELTA,'acceleration'),
            ('divz','cylindrical_divergence','axial',-3-bg.DELTA,'acceleration'),
            ('Etheta','cylindrical_remainder','theta',-3+bg.DELTA,'acceleration'),
            ('Ez','cylindrical_remainder','axial',-3+bg.DELTA,'acceleration'))
        for name,section,label,beta,unit in selected:
            row=view[section][label]
            insert(name,bg.box.pulse.radius.FunctionRef(g,row['exact_reference_log_scale_function']),
                row['common_scale_coefficient_enclosure'],unit,
                self.background._exact(bg.s.sympify(beta)/2,coordinates),row)
        ratios=self.before._fields[id(sample.remainder)]
        for name,beta in (('Er_transport',-3),('Er_axial_viscosity',-3+2*bg.DELTA)):
            row=ratios[2][name]
            insert(name,bg.box.pulse.radius.FunctionRef(g,row['exact_reference_log_scale_function']),
                row['common_scale_coefficient_enclosure'],'acceleration',
                self.background._exact(bg.s.sympify(beta)/2,coordinates),row)
        row=ratios[2]['Er']
        insert('Er_total',bg.box.pulse.radius.FunctionRef(g,row['exact_reference_log_scale_function']),
            row['common_scale_coefficient_enclosure'],'acceleration',row=row)
        return result

    def _compare(self,left,right,reader,delta_tau):
        if left['unit']!=right['unit'] or left['sign']!=right['sign']:
            raise ValueError('Same-unit, same-sign cross-time observable required')
        difference=self.graph.sub(right['scale'],left['scale'])
        bound=reader.at(difference)
        coefficient=right['coefficient']/left['coefficient']
        if bg.ends(coefficient)[0]<=0:raise ValueError('Positive magnitude ratio required')
        scale,method=bg.signed.ratio_enclosure(self.ctx,bound)
        magnitude=bound+self.ctx.ln(coefficient)
        fit=magnitude/(self.ctx.mpf(delta_tau.numerator)/delta_tau.denominator)
        exponent=left['expected_tau_exponent'];expected=None;identity=None
        if exponent is not None:
            if reader.polynomial(self.graph.sub(exponent,right['expected_tau_exponent'])):
                raise ValueError('Same original similarity exponent required')
            identity=self.graph.sub(difference,self.graph.mul(self.graph.constant(delta_tau),exponent))
            if reader.polynomial(identity):raise ValueError('Actual scale change differs from original source exponent')
            expected=reader.at(exponent)
        lo,hi=bg.ends(magnitude)
        return dict(sign=right['sign'],physical_unit=left['unit'],
            exact_scale_log_change=difference.node,directed_scale_log_change=bound,
            actual_coefficient_magnitude_ratio=coefficient,coefficient_ratio_not_replaced_by_one=True,
            ordinary_magnitude_ratio=None if scale is None else coefficient*scale,scale_ratio_method=method,
            directed_log_magnitude_ratio=magnitude,fitted_tau_exponent=fit,
            original_expected_tau_exponent=None if exponent is None else exponent.node,
            directed_original_expected_tau_exponent=expected,
            exact_scale_exponent_identity=None if identity is None else identity.node,
            actual_magnitude_increase_certified=lo>0,actual_magnitude_decrease_certified=hi<0,
            absolute_physical_materialization_not_required=True,absolute_width_flags_not_promoted=True)

    @source_precision
    def evaluate(self,samples):
        if type(samples) is not dict or len(samples)<3:raise ValueError('At least three named actual time samples required')
        if any(not isinstance(name,str) or not name for name in samples):raise ValueError('Nonempty sample names required')
        self.assert_graph()
        validated={name:self._validate_sample(sample) for name,sample in samples.items()}
        requests=[request for _,request in validated.values()];first=requests[0]
        if any((r.chart,r.native,r.Z,r.theta)!=(first.chart,first.native,first.Z,first.theta) for r in requests):
            raise ValueError('One fixed original similarity-coordinate path required')
        offsets=[r.finite_log_time_offset for r in requests]
        if any(b>=a for a,b in zip(offsets,offsets[1:])):raise ValueError('Time samples must approach the critical time in order')
        deliveries=[record['delivery'] for record,_ in validated.values()]
        reader=self._merged_reader(deliveries)
        observations={name:self._observables(samples[name],record,request) for name,(record,request) in validated.items()}
        comparisons={};names=list(samples)
        for left,right in zip(names,names[1:]):
            q=validated[right][1].finite_log_time_offset-validated[left][1].finite_log_time_offset
            tau_delta=self.graph.sub(validated[right][1].forward_coordinates['log_tau'],validated[left][1].forward_coordinates['log_tau'])
            if reader.polynomial(self.graph.sub(tau_delta,self.graph.constant(q))):raise ValueError('Actual time difference differs')
            fits={label:self._compare(observations[left][label],observations[right][label],reader,q)
                for label in observations[left]}
            # Sector ratios are indicators at each issued point. Cross-time changes
            # are computed from their two actual source values, not from rounded ratios.
            relative_fits={}
            for label,(numerator,denominator) in relative.PAIRS.items():
                a=observations[left][numerator if numerator!='Er' else 'Er_total']
                b=observations[left][denominator]
                c=observations[right][numerator if numerator!='Er' else 'Er_total']
                d=observations[right][denominator]
                def quotient(x,y):
                    return dict(unit='dimensionless',sign=x['sign']*y['sign'],
                        scale=self.graph.sub(x['scale'],y['scale']),coefficient=x['coefficient']/y['coefficient'],
                        expected_tau_exponent=None if x['expected_tau_exponent'] is None else
                            self.graph.sub(x['expected_tau_exponent'],y['expected_tau_exponent']))
                relative_fits[label]=self._compare(quotient(a,b),quotient(c,d),reader,q)
            comparisons[left+'__to__'+right]=dict(delta_log_tau=str(q),observable_fits=fits,
                finite_remainder_relative_fits=relative_fits)
        scales=[item['scale'] for obs in observations.values() for item in obs.values()]
        scales.extend(bg.box.pulse.radius.FunctionRef(self.graph,fit[key]) for comparison in comparisons.values()
            for fits in (comparison['observable_fits'],comparison['finite_remainder_relative_fits']) for fit in fits.values()
            for key in ('exact_scale_log_change','original_expected_tau_exponent','exact_scale_exponent_identity') if fit[key] is not None)
        dag={node:data for scale in scales for node,data in self.arithmetic._scale_dag(scale).items()}
        value=OriginalTimeObservations(first.chart)
        snapshots={name:{label:(item['coefficient']._mpi_,item['scale'].node,item['sign'],item['unit'])
            for label,item in obs.items()} for name,obs in observations.items()}
        self._fields[id(value)]=(value,dict(samples),validated,observations,comparisons,dag,snapshots,
            copy.deepcopy(bg.correlated.report(comparisons)))
        return value

    @source_precision
    def sample(self,offsets,native=Fraction(21,2),Z=Fraction(371,1000),theta='7/10'):
        offsets=[Fraction(q) for q in offsets]
        if len(offsets)<3 or any(b>=a for a,b in zip(offsets,offsets[1:])):
            raise ValueError('At least three decreasing exact time offsets required')
        source=self.product.amplitude.before
        deliveries=[source.source(source.invert(source.physical_input('pulse_exit',native,Z,str(q),theta))) for q in offsets]
        samples={}
        for q,delivery in zip(offsets,deliveries):
            field=self.differential.evaluate(delivery,'1/1000')
            packet=self.background.evaluate(field,'1/1000')
            samples[str(q)]=OriginalTimeSample(field,self.axial.evaluate(field),packet,
                self.cone.evaluate(field,'1/1000'),self.before.evaluate(packet,'1/1000'))
        return self.evaluate(samples)

    @source_precision
    def report(self,value):
        entry=self._fields.get(id(value))
        if type(value) is not OriginalTimeObservations or entry is None or entry[0] is not value:
            raise ValueError('Live multi-time observations issued by this owner required')
        if any(self.graph.nodes[node]!=data for node,data in entry[5].items()):raise ValueError('Time scale DAG changed')
        if bg.correlated.report(entry[4])!=entry[7]:raise ValueError('Time comparisons changed')
        self.assert_graph()
        points={}
        for name,sample in entry[1].items():
            record,request,packet_views=self._validate_sample(sample,return_views=True)
            if record is not entry[2][name][0] or request is not entry[2][name][1]:raise ValueError('Time delivery changed')
            for label,item in entry[3][name].items():
                if (item['coefficient']._mpi_,item['scale'].node,item['sign'],item['unit'])!=entry[6][name][label]:
                    raise ValueError('Time observable changed')
            relative_view=packet_views['remainder']['source_correlated_ratios']
            log_magnitudes={}
            for label,ratio in relative_view.items():
                bound=ratio['directed_scale_log_difference']+self.ctx.ln(abs(ratio['signed_coefficient_ratio_enclosure']))
                lo,hi=bg.ends(bound)
                log_magnitudes[label]=dict(directed_log_absolute_relative_value=bound,
                    magnitude_strictly_below_one=hi<0,magnitude_strictly_above_one=lo>0,
                    ordinary_absolute_value_not_required=True)
            points[name]=dict(chart=request.chart,native=str(request.native),Z=str(request.Z),theta=str(request.theta),
                finite_log_time_offset=str(request.finite_log_time_offset),
                values={label:dict(exact_log_scale=item['scale'].node,signed_source_coefficient=item['coefficient'],
                    sign=item['sign'],physical_unit=item['unit']) for label,item in entry[3][name].items()},
                strict_point_cone_certified=packet_views['cone'][relative.cone.POINT_CONE],
                finite_remainder_relative_indicators=relative_view,
                finite_remainder_relative_log_magnitudes=log_magnitudes)
        return dict(source_family=copy.deepcopy(self.family_record),actual_time_samples=points,
            time_comparisons=copy.deepcopy(entry[4]),fixed_original_similarity_coordinates=True,
            actual_multitime_comparisons_available=True,
            actual_time_growth_observed=any(comparison['observable_fits'][label]['actual_magnitude_increase_certified']
                for comparison in entry[4].values() for label in ('Ur','Utheta','Uz','omega_z')),
            complete_nonzero_Z_lambda_retained=True,
            independently_measured_vortex_core_width=False,material_trajectory_or_winding=False,
            coefficient_recursion_implemented=False,all_orders_time_flatness_certified=False,
            global_or_regional_cone_certified=False,physical_absolute_accuracy_not_promoted=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))


@source_precision
def run(before,samples):
    began=time.monotonic();owner=CurrentOriginalRpTimeObservations(before,require_checked=False)
    value=owner.evaluate(samples);view=owner.report(value)
    result=dict(source_family=owner.family_record,actual_multitime_original_observations=view,
        input_hashes=owner.hashes,source_assertions=owner.assert_graph(),execution_seconds=time.monotonic()-began,
        **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(bg.correlated.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,value,view

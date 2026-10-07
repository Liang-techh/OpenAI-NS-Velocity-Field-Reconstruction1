"""Whole-source Duhamel defect definitions and N-polynomial log majorants.

The union of checked original covers bounds the intended source function.
Overlapping chart widths are never added. Actual signed cumulative defects
have a source-linked integral contract, not selected values from covers.
This stage supplies quantitative whole-domain bounds and quiet transport,
not a point backend, evaluated integral, new repair or finite N admission.
"""
import ast
import gzip
import json
from pathlib import Path
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_loop_function_sources as source
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as recovery

HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
packets=source.packets;LogUpper=source.current.bounds.LogUpper
NAME=PREFIX+'current_generic_five_defect_bounds.json'
RECEIPT=PREFIX+'current_generic_five_defect_bounds_check.json'
GATE='current_whole_generic_loop_five_defect_Duhamel_log_majorants_certified'
FUNCTION_GATE='current_original_generic_loop_defect_integral_functions_installed'
OPEN=source.OPEN
RATES=dict(recovery.RATES)
ORDERS=('value','Z','y','yZ','yy')


def read_poly(c,row):
    terms={}
    for term in row['terms']:
        power=term['N_power'];cap=term['coefficient']
        if type(power) is not int or power in terms or cap['exact_zero']:
            raise ValueError('Unique nonzero integer frequency terms required')
        terms[power]=LogUpper(c,packets.interval(c,cap['log_absolute_upper']))
    if row['exact_zero']!=(not bool(terms)):raise ValueError('Exact-zero frequency contract differs')
    return source.FrequencyLogBound(c,terms)


def union_poly(c,rows):
    """Per-power max, a cover union rather than a sum of chart integrals."""
    polys={chart:read_poly(c,row) for chart,row in rows.items()}
    result={};witnesses={}
    ep=packets.recovery.endpoints
    for power in sorted({p for poly in polys.values() for p in poly.terms}):
        chart=max((key for key,poly in polys.items() if power in poly.terms),
                  key=lambda key:ep(polys[key].terms[power].log)[1])
        result[power]=LogUpper(c,polys[chart].terms[power].log);witnesses[str(power)]=chart
    return source.FrequencyLogBound(c,result),witnesses


def cumulative_bounds(c,values,Zsources,ysources,log_length_upper):
    """Zero-inlet D, D_Z and ODE derivatives for five actual signed rates."""
    C=lambda value:LogUpper.constant(c,value)
    length=LogUpper(c,log_length_upper)
    result={};quiet={}
    for key,rate0 in RATES.items():
        rate=c.mpf(rate0);mass=length if rate0==0 else C(1/rate)
        D=values[key].scale(mass);DZ=Zsources[key].scale(mass)
        Dy=values[key]+D.scale(C(rate))
        DyZ=Zsources[key]+DZ.scale(C(rate))
        Dyy=ysources[key]+Dy.scale(C(rate))
        result[key]=dict(value=D.record(),Z=DZ.record(),y=Dy.record(),yZ=DyZ.record(),yy=Dyy.record())
        # log(Rc/r_plus)=1, before the new repair beginning at Rc.
        decay=LogUpper(c,-rate);inlet=D.scale(decay);inletZ=DZ.scale(decay)
        quiet[key]=dict(value=inlet.record(),Z=inletZ.record(),
            y=inlet.scale(C(rate)).record(),yZ=inletZ.scale(C(rate)).record(),
            yy=inlet.scale(C(rate*rate)).record(),
            signed_value_transport='D_Rc=exp(-rate)*D_rplus',
            pressure_memory_exactly_preserved=rate0==0)
    return result,quiet


def integral_contracts(records):
    """Source-linked transport contracts; covers are not defining functions."""
    return {
        key:dict(rate=str(rate),inlet_defect_exact_zero=True,
            actual_original_inlet_history_unchanged=True,
            required_source_function='one original Section11 changed signed density, with common_N and global phi; owner-injected seam identity still required',
            source_cover_references={chart:dict(views=source.VIEWS,
                density_node=record['five_signed_increment_rate_roots'][key],
                Z_density_node=record['five_signed_increment_rate_first_derivatives']['Z'][key],
                y_density_node=record['five_signed_increment_rate_first_derivatives']['y'][key])
                for chart,record in records.items()},
            coordinate='y=log(R/r_minus), y_plus=log(r_plus/r_minus)',
            zero_extension_before_y0=True,zero_extension_after_yplus=True,
            inside_function=dict(operation='signed_Duhamel_integral',bound_variable='s',
                lower='0',upper='y',kernel='exp(-rate*(y-s))',
                integrand='delta_density_j(r_minus*exp(s),Z,phi=fractional_part(common_N*s))'),
            after_function='D_j(y,Z)=exp(-rate*(y-y_plus))*D_j(y_plus,Z)',
            ordinary_derivatives=dict(Z='integral_0^y exp(-rate*(y-s))*delta_density_j_Z(s,Z)ds',
                y='delta_density_j(y,Z)-rate*D_j(y,Z)',
                yZ='delta_density_j_Z(y,Z)-rate*D_j_Z(y,Z)',
                yy='delta_density_j_y_total(y,Z)-rate*D_j_y(y,Z)'),
            unchanged_original_axis_pressure_P0=True,
            own_history_function='original_normalized_history_j+D_j',
            source_chart_selector_not_a_fast_phase_or_independent_function=True,
            source_cover_references_do_not_define_point_functions=True,
            actual_source_owner_and_seam_identity_admitted=False,
            numerical_integral_or_point_field_evaluated=False)
        for key,rate in RATES.items()}


def exact_theorem():
    y,Z,t,rate=s.symbols('y Z t rate',real=True)
    f=s.Function('signed_source');checks={}
    D=s.Integral(s.exp(-rate*(y-t))*f(t,Z),(t,0,y))
    identities={
        'zero_inlet':D.subs(y,0),
        'D_y_plus_rate_D_is_source':s.diff(D,y)+rate*D-f(y,Z),
        'D_Z_same_kernel':s.diff(D,Z)-s.Integral(s.exp(-rate*(y-t))*s.diff(f(t,Z),Z),(t,0,y)),
        'D_yZ_is_source_Z_minus_rate_D_Z':s.diff(D,y,Z)-s.diff(f(y,Z),Z)+rate*s.diff(D,Z),
        'D_yy_is_total_source_y_minus_rate_D_y':s.diff(D,y,2)-s.diff(f(y,Z),y)+rate*s.diff(D,y)}
    for label,value in identities.items():
        if s.simplify(value)!=0:raise ArithmeticError('Signed Duhamel identity failed: '+label)
        checks[label]=True
    inlet,zjet,w=s.symbols('signed_inlet signed_inlet_Z positive_quiet_width')
    for key,rate0 in RATES.items():
        rr=s.Rational(str(rate0));quiet=s.exp(-rr*w)*inlet
        if s.diff(quiet,w)+rr*quiet!=0:raise ArithmeticError('Quiet defect memory reset')
        checks['quiet_'+key+'_signed_transport']=True
        if rate0==0 and quiet!=inlet:raise ArithmeticError('Pressure memory must be exact')
    tree=ast.parse((HERE/(PREFIX+'current_generic_shear_moment_recovery.py')).read_text(encoding='utf8'))
    rate_node=next(node.value for node in tree.body if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='RATES' for t in node.targets))
    declared={kw.arg:ast.literal_eval(kw.value) for kw in rate_node.keywords}
    if declared!=RATES:raise ValueError('Original normalized recovery rates changed')
    checks['original_normalized_five_rates_AST_bound']=True
    return dict(passed=True,exact_identities=checks,
        positive_rate_kernel_mass_bound='integral_0^width exp(-rate*t)dt<=1/rate for rate>0',
        pressure_rate0_kernel_mass_bound='width<=logRw+1-logRa on the original modification interval',
        kernel_nonnegative_no_absolute_signed_density_substitution=True,
        Z_independent_original_radial_endpoints_required=True,
        overlapping_source_covers_not_integrated_separately=True)


class CurrentFiveDefectBounds:
    def __init__(self):
        self.ctx=MPIntervalContext();self.ctx.dps=240;self.hashes={}
        self.manifests={};self.receipts={}
        prerequisites={
            'current_generic_loop_function_sources':source.GATE,
            'current_generic_shear_loop_domain':'current_original_generic_loop_two_sided_collars_and_repair_geometry_certified',
            'current_generic_shear_uniform_inputs':'current_original_whole_generic_input_margins_and_log_scales_certified'}
        for stem,gate in prerequisites.items():
            name=PREFIX+stem+'.json';check=PREFIX+stem+'_check.json'
            manifest=json.loads((HERE/name).read_bytes());receipt=json.loads((HERE/check).read_bytes())
            if not receipt['all_passed'] or not receipt[gate] or not manifest[gate] or any(receipt.get(k) for k in OPEN):
                raise ValueError('Checked current source/domain/whole-cover stage required: '+stem)
            for filename,digest in receipt['input_hashes'].items():
                if sha(filename)!=digest or filename in self.hashes and self.hashes[filename]!=digest:
                    raise ValueError('Changed or conflicting whole-source prerequisite: '+filename)
                self.hashes[filename]=digest
            self.hashes[name]=sha(name);self.hashes[check]=sha(check)
            self.manifests[stem]=manifest;self.receipts[stem]=receipt
        self.family=self.manifests['current_generic_loop_function_sources']['source_family']
        if any(row['source_family']!=self.family for row in (*self.manifests.values(),*self.receipts.values())):
            raise ValueError('One actual source/family/P0 across the whole integral required')
        self.records=json.loads(gzip.decompress((HERE/source.VIEWS).read_bytes()))
        uniform=self.manifests['current_generic_shear_uniform_inputs'];whole=uniform['whole_actual_original_generic_input_margin']
        inventory=uniform['current_actual_logarithmic_loop_scales']['whole_source_chart_inventory']
        grouped={chart for row in whole['chart_groups'].values() for chart in row['source_charts']}
        if set(self.records)!=set(inventory) or grouped!=set(inventory) or len(inventory)!=17:
            raise ValueError('Same whole original17-chart source cover required')
        if not whole['original_weak_and_strong_signed_branches_on_entire_modification_domain_certified']:
            raise ValueError('Whole original r_minus..r_plus cover missing')
        self.domain=self.manifests['current_generic_shear_loop_domain']['current_original_generic_loop_domain']
        if self.domain!=uniform['current_chosen_domain']:raise ValueError('Different original microscopic endpoints')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.theorem=exact_theorem()

    def compute(self):
        c=self.ctx;ep=packets.recovery.endpoints;read=lambda row:packets.interval(c,row)
        logRa=read(self.domain['exact_log_Ra_without_lost_offset'])
        logRw=read(self.domain['exact_log_Rw'])
        # The positive left offset is dropped only in this upper estimate.
        # Keep the exact endpoint tree in integral_contracts/domain_binding.
        length_upper=logRw+c.mpf(1)-logRa
        if ep(length_upper)[0]<=0:raise ArithmeticError('Positive whole original log-radius span required')
        log_length_upper=c.ln(c.mpf(ep(length_upper)[1]))
        envelopes={};witnesses={}
        for order in ('value','Z','y'):
            envelopes[order]={};witnesses[order]={}
            for key in RATES:
                rows={chart:record['actual_source_frequency_majorants']['five_increment_rate_log_polynomials'][key]
                    if order=='value' else record['actual_source_frequency_majorants']['five_increment_rate_first_log_polynomials'][order][key]
                    for chart,record in self.records.items()}
                envelopes[order][key],witnesses[order][key]=union_poly(c,rows)
        required_logN=c.mpf(max(ep(read(v['actual_source_frequency_majorants']['required_positive_log_N_lower']))[1] for v in self.records.values()))
        bounds,quiet=cumulative_bounds(c,envelopes['value'],envelopes['Z'],envelopes['y'],log_length_upper)
        return dict(source_family=self.family,**{GATE:True,FUNCTION_GATE:False},**dict.fromkeys(OPEN,False),
            actual_whole_source_chart_count=len(self.records),defect_rates=RATES,
            source_linked_defect_integral_contracts=integral_contracts(self.records),
            source_union_frequency_envelopes={order:{key:poly.record() for key,poly in row.items()} for order,row in envelopes.items()},
            source_union_coefficient_max_witnesses=witnesses,
            actual_whole_defect_frequency_majorants=bounds,
            actual_new_repair_inlet_Rc_pre_repair_defect_majorants=quiet,
            required_positive_log_N_lower=required_logN,
            original_domain_binding=dict(manifest=PREFIX+'current_generic_shear_loop_domain.json',
                original_family=self.family,formal_radii=self.domain['formal_radii'],
                exact_left_offset=self.domain['formal_left_log_radius_offset'],
                original_endpoint_coordinate_trees_preserved=True,
                actual_log_radius_span_upper=length_upper,log_log_radius_span_upper=log_length_upper,
                upper_span_formula='logRw+1-logRa; only positive hb*s_c/2 discarded for an upper bound',
                exact_quiet_log_Rc_over_rplus='1',Z_independent_original_radius_boundaries=True,
                loop_source_exact_zero_after_rplus_by_checked_original_power_branch=True,
                post_Rc_new_repair_not_installed=True),
            exact_signed_Duhamel_and_quiet_transport_theorem=self.theorem,
            whole_source_value_Z_y_yZ_defects_have_only_negative_N_powers=True,
            second_y_defects_retain_fast_source_N_power0=True,
            ordinary_y_yZ_yy_derivatives_scope='open original chart interiors and one-sided traces; seam equality not admitted',
            actual_source_domain_partition_or_overlap_double_count_not_required_for_union_supremum=True,
            unchanged_original_axis_pressure_P0=True,pressure_quiet_memory_preserved=True,
            signed_current_point_loop_or_inverse_jets_installed=False,
            actual_changed_five_moment_transport_integrated=False,
            actual_own_histories_radial_pressure_stress_point_evaluated=False,
            sufficient_mixed4_velocity_or_mixed3_stress_admitted=False,
            current_whole_N_selected=False,source_graph_ancestor_constructors_called=False,
            scope='Rigorous17-chart union N-polynomial log majorants for the source-linked zero-inlet five-defect Duhamel contract, chart-interior ordinary y,Z,yZ,yy and signed quiet carry to the new repair inlet Rc. Actual source-function replay and seam identities are not supplied by covers; no installed/evaluated integral/history, repair, common finite N/global cone or actual recursion.',
            input_hashes=self.hashes)

    def run(self):
        row=self.compute();(HERE/NAME).write_text(json.dumps(packets.encode(row),indent=2)+'\n',encoding='utf8')
        print('Whole original five Duhamel defect log majorants:17-chart union; source-function/seam admission remains open',flush=True)
        return row


def run():return CurrentFiveDefectBounds().run()


if __name__=='__main__':run()

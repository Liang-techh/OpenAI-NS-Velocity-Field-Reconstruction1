"""Uniform absolute pressure-prefix models on the original real Z domain.

The affine coefficients enclose genuine values at their centers. A bound
for the same original function's second derivative supplies a Taylor
remainder, including at zeros. No model coefficient is selected as the
actual field; all original collar errors and physical factors remain.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_unchanged_pressure_prefix as prior

op, domain, ep = prior.op, prior.domain, prior.ep
HERE, PREFIX, sha, source_precision = prior.HERE, prior.PREFIX, prior.sha, prior.source_precision
budget = prior.moments.core.budget
NAME = PREFIX+'current_original_n1_unchanged_axial_pressure_models.json.gz'
RECEIPT = PREFIX+'current_original_n1_unchanged_axial_pressure_models_check.json'
GATE = 'current_original_n1_unchanged_full_Z_absolute_pressure_values_installed'
OPEN = tuple(dict.fromkeys(prior.OPEN+('whole_Z_relative_accuracy_certified',
    'whole_Z_first_derivative_absolute_accuracy_certified')))


def read_interval(c, value):
    return c.make_mpf((tuple(value['lower_exact_mpf_tuple']),tuple(value['upper_exact_mpf_tuple'])))


def norm(c, value):
    return domain.upper(c,abs(value))


def hull(c, values):
    return c.mpf([min(ep(v)[0] for v in values),max(ep(v)[1] for v in values)])


def normalized_row(value, order, reference, normalization):
    """Normalize actual source sectors with one fixed positive bound scale."""
    c=value.ctx
    if ep(normalization)[0]<=0:raise ValueError('Positive bound normalization required')
    total=c.mpf(0)
    for powers,jet in sorted(value.terms.items()):
        endpoints=[]
        for endpoint in ep(jet[order]):
            if endpoint==0:endpoints.append(c.mpf(0));continue
            logarithm=budget.log_ratio(c,value.algebra.logs,powers,reference,
                abs(c.mpf(endpoint)),normalization)
            upper,rule=budget.exp_cover(c,logarithm)
            if upper is None:raise ValueError('Unresolved positive relative source scale')
            lower=c.mpf(0) if rule=='directed_small_log_cover' else domain.lower(c,c.exp(c.mpf(ep(logarithm)[0])))
            covered=c.mpf([ep(lower)[0],ep(upper)[1]])
            endpoints.append(-covered if endpoint<0 else covered)
        total+=c.mpf([ep(endpoints[0])[0],ep(endpoints[1])[1]])
    return total


def affine_error(c, coefficients, halfwidth, curvature_coefficient):
    """Jet row 2 is f''/2; bound both Taylor remainder and coefficient spread."""
    radius=domain.upper(c,norm(c,curvature_coefficient)*halfwidth**2)
    ambiguity=domain.upper(c,coefficients[0].b-coefficients[0].a
        +halfwidth*(coefficients[1].b-coefficients[1].a))
    return radius,ambiguity,domain.upper(c,radius+ambiguity)


@source_precision
def original_values(owner, Z):
    """Same original pressure recipe as the accepted prefix, now through Z2."""
    base=owner.prefix;c=base.c;z,algebra,_=base.inputs.frame(Z)
    a,b=base.phase_in,base.modification_inlet
    from lei_ren_part1_paper_compliant_current_original_n1_local_leading_dispatch import original
    support=original.OriginalMicroFiniteN.modification_support(base.local.selector,'first_micro',0,b)
    if not support['increments_exact_zero']:raise ValueError('Original unchanged support required')
    cap=c.mpf('1e-80');loghbphase=algebra.logs[0]+c.ln(b)
    if ep(loghbphase)[1]>=ep(c.ln(cap/2))[0]:raise ValueError('Original rho allowance not proved')
    rho=c.mpf([4,ep(4*(1+cap))[1]])
    if ep(rho)[1]>ep(base.mean.bound['actual_inner_geometry']['real_rho_upper_bound'])[1]:
        raise ValueError('Original core bound required')
    rows,fields,proof=prior.refined_original_rows(base,algebra,z,rho,b)
    gradient=op.Jet(c,base.mean.weighted.source.amplitude.gradient_jets(z))
    F=[algebra.shift(phi*prior.mean.weighted.relative_amplitude(algebra,gradient,1,n),(0,1,0,0))
       for phi,n in zip(rows['Phi'],prior.mean.ORDERS)]
    forcing=op.known_forcing(algebra,dict(F=F,Uz=fields['Uz'],Q=fields['Q'],mean=fields['mean']),
        algebra.lift(op.Jet.variable(c,z,5)),rho,base.physical.delta)['N1p']
    values=prior.physical_pressure_range_integrals(algebra,forcing,rho,a,b)
    for value in values.values():
        if value.order!=2 or any(p[1]!=0 for p in value.terms):
            raise ValueError('Genuine Z2 rows and no Z-dependent F0 source factor required')
    return values,proof


class CurrentOriginalN1UnchangedAxialPressureModels:
    @source_precision
    def __init__(self,accepted_prefix,require_checked=True):
        if type(accepted_prefix) is not prior.CurrentOriginalN1UnchangedPressurePrefix or not accepted_prefix.acceptance_loaded:
            raise ValueError('Accepted original pressure prefix required')
        accepted_prefix.assert_graph();self.prefix=accepted_prefix;self.c=accepted_prefix.c
        self.family=copy.deepcopy(accepted_prefix.family);self.hashes=dict(accepted_prefix.hashes)
        for name in (prior.NAME,prior.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self.data=None
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Axial pressure model acceptance/scope differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Axial pressure model input changed '+name)
                self.hashes[name]=digest
            self.data=json.loads(gzip.decompress((HERE/NAME).read_bytes()))
            if self.data['source_family']!=self.family:raise ValueError('Axial pressure model source differs')
            self.acceptance_loaded=True
        self._prefix=accepted_prefix;self._c=self.c;self._hashes=copy.deepcopy(self.hashes)
        self._family=copy.deepcopy(self.family);self._data=op.snapshot(self.data)

    def assert_graph(self):
        if self.prefix is not self._prefix or self.c is not self._c or self.family!=self._family or self.hashes!=self._hashes:
            raise ValueError('Axial pressure model owner changed')
        if op.snapshot(self.data)!=self._data:raise ValueError('Axial pressure model data changed')
        self.prefix.assert_graph()

    @source_precision
    def evaluate(self,Z):
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Checked complete axial models required')
        return evaluate_models(self.c,self.data,Z,self.acceptance_loaded)


@source_precision
def evaluate_models(c,data,Z,accepted=False):
    z=c.mpf(Z)
    if ep(z)[0]<-1 or ep(z)[1]>1:raise ValueError('Original real Z domain [-1,1] required')
    values={k:[] for k in ('m13','m15')};derivatives={k:[] for k in values};active=[]
    for index,cell in enumerate(data['cells']):
        left,right=cell['exact_Z_cell'];lo=max(ep(z)[0],c.mpf(left).a);hi=min(ep(z)[1],c.mpf(right).b)
        if lo>hi:continue
        offset=c.mpf([lo,hi])-c.mpf(cell['exact_center'])
        active.append(index)
        for name,model in cell['models'].items():
            a,b=[read_interval(c,v) for v in model['affine_coefficients_Z0_Z1']]
            radius=read_interval(c,model['Taylor_remainder_abs_upper'])
            curvature=read_interval(c,data['normalizations'][name]['full_domain_Z2_Taylor_coefficient'])
            values[name].append(a+b*offset+c.mpf([-1,1])*radius)
            derivatives[name].append(b+2*curvature*offset)
    if not active:raise ValueError('Missing original axial model coverage')
    return dict(Z=z,active_cells=active,normalized_values={n:hull(c,v) for n,v in values.items()},
        normalized_first_derivatives={n:hull(c,v) for n,v in derivatives.items()},
        normalizations=copy.deepcopy(data['normalizations']),
        value_uniform_absolute_error_target=data['normalized_absolute_error_target'],
        value_error_relative_to_global_bound_norm_not_local_nonzero_value=True,
        interval_range_width_includes_actual_axial_variation=True,
        first_derivative_range_enclosure_only=True,
        **{GATE:accepted},**dict.fromkeys(OPEN,False))


@source_precision
def run(accepted_prefix):
    began=time.monotonic();owner=CurrentOriginalN1UnchangedAxialPressureModels(accepted_prefix,False);c=owner.c
    full,proof=original_values(owner,['-1','1']);normalizations={};curvatures={}
    bases=next(iter(full.values())).algebra.logs
    for name,value in full.items():
        b=budget.component_budget(value,0,c.mpf(0));reference=b['reference_source_exponents'];M=b['directed_reference_coefficient_norm']
        if reference[1]!=0 or ep(M)[0]<=0:raise ValueError('Positive Z-independent pressure normalization required')
        curvatures[name]=normalized_row(value,2,reference,M)
        normalizations[name]=dict(reference_source_exponents=reference,directed_reference_coefficient_norm=M,
            full_domain_value_range=normalized_row(value,0,reference,M),
            full_domain_first_derivative_range=normalized_row(value,1,reference,M),
            full_domain_Z2_Taylor_coefficient=curvatures[name],
            physical_scale='hb^a*Pstar_squared^b*Lambda^c times the positive coefficient norm',
            normalization_is_a_bound_not_a_selected_field_value=True,
            source_basis_Z_independent_indices=[0,2,3])
    target=c.mpf('.001');count=4
    while any(ep(norm(c,v)/count**2)[1]>ep(target/4)[0] for v in curvatures.values()):
        count*=2
        if count>4096:raise ValueError('Uniform absolute model needs a different source bound')
    cells=[]
    for index in range(count):
        left=c.mpf(-1)+c.mpf(2)*index/count;right=c.mpf(-1)+c.mpf(2)*(index+1)/count;center=(left+right)/2
        point,_=original_values(owner,center);models={}
        for name,value in point.items():
            if any(value.algebra.logs[j]._mpi_!=bases[j]._mpi_ for j in (0,2,3)):
                raise ValueError('Original physical pressure factors must be constant in Z')
            scale=normalizations[name];ref=scale['reference_source_exponents'];M=scale['directed_reference_coefficient_norm']
            coefficients=[normalized_row(value,k,ref,M) for k in (0,1)]
            radius,ambiguity,total=affine_error(c,coefficients,(right-left)/2,curvatures[name])
            if ep(total)[1]>ep(target)[0]:raise ValueError('Uniform absolute pressure model error exceeds target')
            models[name]=dict(affine_coefficients_Z0_Z1=coefficients,Taylor_remainder_abs_upper=radius,
                coefficient_ambiguity_abs_upper=ambiguity,total_absolute_model_error_upper=total,
                true_center_value_and_derivative_enclosures=True)
        cells.append(dict(exact_Z_cell=[str(ep(left)[0]),str(ep(right)[1])],exact_center=str(ep(center)[0]),models=models))
        if (index+1)%16==0:print('AXIAL_PRESSURE_MODELS',index+1,'of',count,flush=True)
    raw=dict(source_family=owner.family,exact_Z_domain=['-1','1'],cells=cells,normalizations=normalizations,
        normalized_absolute_error_target='.001',original_source_comparison_proof=proof,
        original_source_log_bases=[bases[j] for j in (0,2,3)],
        original_physical_measures_and_m15_sign_retained=True,
        Z2_is_half_the_ordinary_second_derivative=True,continuous_domain_including_zeros=True,
        affine_model_is_an_enclosure_not_the_actual_field=True,input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    snapshot=op.snapshot(raw)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(snapshot,separators=(',',':'))+'\n').encode(),mtime=0))
    owner.assert_graph()
    print('UNCHANGED_FULL_Z_ABSOLUTE_PRESSURE_MODELS_COMPUTED',count,flush=True)
    return owner,snapshot

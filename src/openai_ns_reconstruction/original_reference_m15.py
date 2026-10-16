"""Original reference m15, preserving exact long-end coefficient correlations."""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_long_m15 import OriginalLongM15
from .original_reference_pressure_rows import OriginalReferencePressureRows


class OriginalReferenceM15(OriginalLongM15):
    def __init__(self, *, rows=None, **kwargs):
        if rows is not None and kwargs:
            raise ValueError('Supply original reference pressure source or options')
        self.rows = rows if rows is not None else OriginalReferencePressureRows(**kwargs)
        if not isinstance(self.rows, OriginalReferencePressureRows) or self.rows.chart != 'reference':
            raise ValueError('The original constant-U reference chart is required')
        self.long_engine = OriginalLongM15(rows=self.rows.long_rows)
        for name in ('t','o','ctx','iv','J','audit','directory','family'):
            setattr(self,name,getattr(self.long_engine,name))
        self.hashes,self.code_hashes=dict(self.long_engine.hashes),dict(self.long_engine.code_hashes)
        self.hashes.update(self.rows.hashes);self.code_hashes.update(self.rows.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_long_m15_numeric.json')
        raw=path.read_bytes();data=json.loads(raw)
        if (not data['checks']['passed'] or data['source_family']!=self.family
                or not data['limits']['long_cumulative_signed_m15_with_complete_factored_errors']):
            raise ValueError('Complete checked original finite-long pressure inlet required')
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        self.inlets={Fraction(p['Z']):p for p in data['samples']
                     if p['coordinate']==dict(kind='exact_decimal_long_fraction',value='1')}
        self._packets={};self._leading={};self.chart='reference'
        name='src/openai_ns_reconstruction/original_reference_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            self.inlets[key]=self.long_engine.evaluate(Z,'1')
            self.hashes.update(self.long_engine.hashes);self.code_hashes.update(self.long_engine.code_hashes)
        return self.inlets[key]

    def leading(self,Z,q):
        key=Z,q._mpi_
        if key not in self._leading:
            # The same constant-U source discrepancies remain valid after
            # any further positive mean kernel: each target is inside the
            # existing error envelope. Rescale only the exact integration
            # bound1/2, never an original source parameter or width cap.
            base=self.long_engine.leading(Z,self.iv.mpf(1))
            if base['native_positive_mass_bound']._mpi_!=self.iv.mpf('.5')._mpi_:
                raise ValueError('Checked full long window must admit positive mass1/2')
            gap=self.rows.geometry(q)['local']
            mass=self.iv.mpf(min(_bounds(self.audit,gap)[1],self.audit.mpf('.5')))
            self._leading[key]={kind:{d:row*(2*mass) for d,row in base[kind].items()}
                                for kind in ('inlet_errors','source_errors')}
        return self._leading[key]

    def finite_errors(self,Z,q,packet):
        _,base,geo=self.rows.background(Z,q);owner=self.rows.long_rows
        lowpacket=owner.evaluate(Z,'1');_,lowbase,lowgeo=owner.background(Z,self.iv.mpf(1))
        up=self.correction_densities(Z,packet,base);low=self.correction_densities(Z,lowpacket,lowbase)
        volume=self.correction_densities(Z,packet,self.rows.uniform_background(Z,q))['volume']
        result={}
        for name,power in (('endpoint',2),('viscous',1)):
            result[name]=[self.rows.sum_logs((None if up[name][k] is None else up[name][k]+power*geo['logR'],
                None if low[name][k] is None else low[name][k]+power*lowgeo['logR'])) for k in (0,1)]
        mass=self.iv.mpf(min(_bounds(self.audit,geo['local'])[1],self.audit.mpf('.5')))
        result['volume']=[None if v is None else v+2*geo['logR']+self.iv.log(mass) for v in volume]
        return result,dict(endpoint=self.rows.log_rows(low['endpoint']),viscous=self.rows.log_rows(low['viscous']),
            log_radius=self.o._summary(lowgeo['logR_point'],lowgeo['logR']))

    def at_reference_gap(self,Z,gap):
        return self.evaluate(Z,'g:'+gap)

    def evaluate(self,Z,coordinate):
        packet=self.rows.evaluate(Z,coordinate);q,_,recipe=self.rows._coordinate(coordinate,False)
        key=Z,coordinate
        if key in self._packets:return self._packets[key]
        c,iv,J,o=self.ctx,self.iv,self.J,self.o;incoming=self.inlet(Z)
        if len(incoming['m15_C0_C1_terms'])!=106 or incoming['actual_selected_integer']!=packet['selected_integer']:
            raise ValueError('Same-N complete106-term long moment inlet required')
        terms=copy.deepcopy(incoming['m15_C0_C1_terms']);source={t['name']:t for t in incoming['m15_C0_C1_terms'][88:]}
        geo=self.rows.geometry(q);gp,gap=geo['local_point'],geo['local'];empty=gap._mpi_==iv.mpf(0)._mpi_
        zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1);zero_rows=o._jet_summary(zero,zeroi)
        if not empty:lead=self.leading(Z,q);finite,lower=self.finite_errors(Z,q,packet)
        else:lead=None;finite={name:[None,None] for name in ('endpoint','volume','viscous')};lower=None
        def growth(k):
            return [o._summary(a*k,b*k) for a,b in geo['log_offset_terms']]
        def physical(d):
            return o._summary(2*c.log(110)-d*o.core.logLambda,2*iv.log(110)-d*o.core.ll_iv)
        def width(name):
            atom=o._summary(sum(_bounds(c,o.logh))/2,o.logh)
            atom.update(source_log_atom_only=True,original_native_parameter_name=name,
                        distinct_graph_nodes_with_checked_common_numeric_source=True)
            return atom
        hbB,hbS=width('hbB'),width('hbS')
        def append(name,rows,atoms,**binding):
            terms.append(dict(name='reference_'+name,coefficient_C0_C1=copy.deepcopy(rows),
                log_scale_terms=copy.deepcopy(atoms),actual_function_binding=dict(chart='reference',
                    incoming106_terms_preserved_once=True,original_finite_reference_length=True,**binding)))
        for d in range(3):
            for n in range(d+1):
                k=2-n
                if k:
                    upper=source[f'long_half_S1_epsilon_{d}_mean_{n}_upper']
                    low=source[f'long_half_S1_epsilon_{d}_mean_{n}_lower']
                    group=upper['actual_function_binding']['correlated_coefficient_id']
                    append(f'half_S1_epsilon_{d}_mean_{n}_upper',zero_rows if empty else upper['coefficient_C0_C1'],
                        low['log_scale_terms']+growth(k),correlated_coefficient_id=group,correlated_coefficient_sign=1,
                        same_original_long_coefficient=True,unit_reference_offset_kept_in_separate_atom=True)
                    append(f'half_S1_epsilon_{d}_mean_{n}_lower',zero_rows if empty else low['coefficient_C0_C1'],
                        upper['log_scale_terms'],correlated_coefficient_id=group,correlated_coefficient_sign=-1,
                        same_original_long_coefficient=True,exact_internal_long_endpoint_telescope=True)
                else:
                    src=source[f'long_half_S1_epsilon_{d}_mean_{n}_resonance'];atoms=src['log_scale_terms'][:1]
                    if not empty:
                        atom=o._summary(c.log(gp),iv.log(gap));atom.update(positive_native_log_measure_only=True,
                            original_measure='gap=Lref*q',reference_gap_not_multiplied_into_coefficient=True);atoms=atoms+[atom]
                    append(f'half_S1_epsilon_{d}_mean_{n}_resonance',zero_rows if empty else src['coefficient_C0_C1'],
                           atoms,exact_resonant_kernel_gap=True)
        for d in (1,2):
            for kind,atoms in (('inlet',[hbB]),('source',[hbS,hbS])):
                error=zeroi if empty else lead[kind+'_errors'].get(d,zeroi)
                append(f'{kind}_error_epsilon_{d}',o._jet_summary(zero,error),[physical(d)]+growth(2)+atoms,
                    error_enclosure_only=True,positive_native_measure_already_integrated=True)
        for name in ('endpoint','volume','viscous'):
            bound=self.rows.sum_logs(finite[name]);atoms=[];truth=zeroi
            if bound is not None:
                top=iv.mpf(_bounds(self.audit,bound)[1]);atom=o._summary(_bounds(c,top)[1],top)
                atom.update(positive_error_bound_atom_only=True,source_value_selected=False);atoms.append(atom)
                truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
            append('finite_N_'+name+'_error',o._jet_summary(zero,truth),atoms,error_enclosure_only=True,
                   both_nonflat_endpoints=True,reference_Ldq_measure_included_once=True)
        out=dict(Z=Z,coordinate=recipe,chart='reference',source_family=self.family,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']),m15_C0_C1_terms=terms,
            incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']),full_prefix_cover=copy.deepcopy(packet['source_cover']),
            finite_N_error_logs_in_physical_units={name:self.rows.log_rows(v) for name,v in finite.items()},
            actual_lower_finite_N_endpoint_budget=lower,constant_U_viscous_increment_exact_zero=True,
            contract=dict(complete_long106term_inlet_once=True,original_reference_gap_and_measure=True,
                exact_shared_long_coefficients_and_ids=True,unit_log_offsets_retained=True,
                inherited_mean_atoms_not_rounded=True,huge_exp_gap_not_materialized=True,
                positive_mass_encloses_finite_Ldq=True,both_finite_N_endpoints=True,
                full_prefix_and_all_quadratic_errors=True,actual_constant_U_mean_ODE=True),
            limits=dict(reference_cumulative_signed_m15_with_complete_factored_errors=True,
                restoration_postrestore_m15=False,resolved_finite_N_signed_centers=False,
                useful_total_relative_accuracy=False,independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out

"""Full original O2 finite-N five C0/Z contributions on the predicate union.

Original mixed ranges are reused as function enclosures. The accepted
signed-density graph executes before branch union; the large ordinary-Z
unit remains formal. This is a direct integral enclosure, not terminal
matching, sharp phase averaging or admission of a global frequency.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_predicate_mixed_phase as phase
import lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport as history

frames=phase.frames;source=phase.source;base,prior,ep=phase.base,phase.prior,phase.ep
positive=source.positive;five=positive.five;density=five.density
MixedJet=phase.MixedJet;C0,Y,Z,YZ=phase.C0,phase.Y,phase.Z,phase.YZ
HERE,PREFIX,sha=phase.HERE,phase.PREFIX,phase.sha
NAME=PREFIX+'current_original_O2_predicate_five_integrals.json.gz'
RECEIPT=PREFIX+'current_original_O2_predicate_five_integrals_check.json'
GATE='original_O2_full_predicate_finite_N_five_C0_Z_direct_integrals_enclosed'
BRANCHES=('regular','positive','negative')
KEYS=tuple(five.RATES)


def interval(c,value):return source.ordered.interval(c,value)


def restore_value(a,row):
    scale=row['formal_positive_scale'];powers=tuple(scale['source_exponents'])+(scale['radius_power'],)
    result=prior.ScaledEnclosure(prior.FormalScale(a.bases,powers,interval(a.ctx,scale['additional_log_interval'])),
        interval(a.ctx,row['coefficient_interval']),a.ledger)
    if result.zero!=row['exact_zero'] or row['point_value_selected'] or not row['encloses_original_source_function']:
        raise ValueError('Accepted actual original function range required')
    return result


class CollectedDensityGraph(positive.CellDensityGraph):
    """Same original graph; sums collect powers in the issued O2 atlas."""
    def __init__(self,atlas,*args):
        self.atlas=atlas;super().__init__(*args)
    def evaluate(self,index):
        if index in self.bindings:return self.bindings[index]
        if index in self.cache:return self.cache[index]
        node=self.nodes[index]
        if node['operation']=='sum':
            self.operations['sum']=self.operations.get('sum',0)+1
            result=self.atlas.sum(self.evaluate(i) for i in node['arguments'])
            self.cache[index]=result;return result
        return super().evaluate(index)
    def modulation_exp(self,index):
        if index in self.expcache:return self.expcache[index]
        x=self.evaluate(index);c=self.c
        finite=base.conditioned.clipped(c,positive.bounded(x),-1,1)
        mean=c.mpf(1);power=c.mpf(1)
        for k in range(1,65):power*=finite;mean+=power/c.factorial(k+1)
        tail=ep(c.exp(1)/c.factorial(66))[1]
        mean=base.conditioned.clipped(c,mean+c.mpf((-tail,tail)),c.exp(-1),c.exp(1))
        changed=x*mean;result=(self.atlas.add(self.scalar(1),changed),changed)
        self.expcache[index]=result
        self.ledger['bounded_modulation_exp_taylor_tail_operations']=self.ledger.get('bounded_modulation_exp_taylor_tail_operations',0)+1
        return result


def original_history_rows(a,roots):
    E,V=roots['E'][C0],roots['V'][C0];EZ,VZ=roots['E'][Z],roots['V'][Z]
    values=dict(m=V,h=E,k=E*V,e=a.add(base.current.square(V),-base.current.square(E)*.5),p=base.current.square(E)*.5)
    derivatives=dict(m=VZ,h=EZ,k=a.add(EZ*V,E*VZ),e=a.add(V*VZ*2,-E*EZ),p=E*EZ)
    return values,derivatives


def export_range(a,value,*,derivative_unit):
    """Export a local q/u cover after dividing the exact common Z unit.

    Unit=Lambda0/L(Z)^2 depends only on Z and original constants. A
    normalized derivative cap is not differentiated as a source function.
    No source large factor is replaced by the unit's interval endpoint.
    """
    if value.scale.bases is not a.bases or value.ledger is not a.ledger:
        raise ValueError('Same local source basis and ledger required')
    unit=(11,10,-2,0,0) if derivative_unit else (0,0,0,0,0)
    remaining=tuple(p-u for p,u in zip(value.scale.powers,unit,strict=True))
    if remaining[0]>0 or remaining[1]>0 or remaining[2]<0 or remaining[3]<0 or remaining[4]>0:
        raise ValueError('Collect original C0/Z large unit before local q/u export: '+str(remaining))
    positive_unit=prior.ScaledEnclosure(prior.FormalScale(a.bases,unit),1,a.ledger)
    ratio=value.positive_divide(positive_unit,positive_unit.scale.evaluate())
    result=positive.bounded(ratio)
    if any(not mp.isfinite(v) for v in ep(result)):raise ArithmeticError('Finite normalized export lost')
    return result


def candidate_N(N):
    if type(N) is not int or N<160:raise ValueError('Explicit integer candidate N>=160 required; not global admission')
    return N


def original_full_inlet(a,z):
    """Same source-defined inlet functions, with exact square at the axis.

    The old strict-sign helper used z*z, whose independent interval product
    touches-1 on [-1,1]. The exact pointwise square is nonnegative and Q>=1.
    """
    c=a.ctx;z=c.mpf(z)
    if not -1<=ep(z)[0]<=ep(z)[1]<=1:raise ValueError('Original closed axial inlet domain required')
    zz=z**2;C=1/(1+zz);scalar=a.scalar
    inverse=prior.ScaledEnclosure(prior.FormalScale(a.bases,(-1,0,0,0,0)),1,a.ledger)
    inv2=base.current.square(inverse)
    old=dict(m=inverse*(4*z),h=scalar(c.mpf(5)/8*C),k=inverse*(c.mpf(5)/2*z*C),
        e=a.add(inv2*(16*zz),-scalar(c.mpf(5)/12*C*C)),p=scalar(c.mpf(5)/2*C*C))
    derivative=dict(m=inverse*4,h=scalar(-c.mpf(5)/4*z*C*C),
        k=inverse*(c.mpf(5)/2*(1-zz)*C*C),
        e=a.add(inv2*(32*z),scalar(c.mpf(5)/3*z*C**3)),p=scalar(-10*z*C**3))
    return dict(original_inlet={key:positive.bounded(value) for key,value in old.items()},
        original_inlet_Z={key:positive.bounded(value) for key,value in derivative.items()},
        native_original_inlet_sources={key:value.record() for key,value in old.items()},
        native_original_inlet_Z_sources={key:value.record() for key,value in derivative.items()},
        original_source_defined_inlet_recipe_unchanged=True,exact_pointwise_Z_square_and_Q_at_least_one=True,
        microscopic_positive_inverse_Pstar_retained=True,whole_Z_rational_source_functions_not_axial_samples=True)


class OriginalO2PredicateFiveIntegrals:
    def __init__(self):
        self.source=frames.OriginalO2PredicateSources();self.ctx=self.source.ctx;self.family=self.source.family
        self.hashes=dict(self.source.hashes)
        for module in (phase,history):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted same original predicate phase and history unit contracts required')
            for name,digest in {**receipt['input_hashes'],module.RECEIPT:sha(module.RECEIPT)}.items():
                if sha(name)!=digest:raise ValueError('Original O2 integral dependency changed: '+name)
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('O2 integral source closures disagree')
                self.hashes[name]=digest
            if module is phase:self.phase_receipt=receipt
        raw=gzip.decompress((HERE/phase.NAME).read_bytes())
        if hashlib.sha256(raw).hexdigest()!=self.phase_receipt['compressed_producer_report']['lossless_original_json_sha256']:
            raise ValueError('Accepted original mixed source archive differs')
        saved=json.loads(raw)
        if saved['source_family']!=self.family or not saved[phase.GATE]:raise ValueError('Same original mixed archive required')
        self.saved=saved['whole_original_O2_predicate_mixed_cells']
        self.count=min(self.source.source.parent.parent.levels)
        if len(self.saved)!=self.count:raise ValueError('Full original source partition required')
        self.out=source.O2MixedAtlas(self.source.owner.inputs.frame,lower=-1,upper=1,logq=self.ctx.mpf(0))
        self.Z_unit=prior.ScaledEnclosure(prior.FormalScale(self.out.bases,(11,10,-2,0,0)),1,self.out.ledger)
        self.inlet_contract=history.original_inlet_contract()
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.archive_cache={};self.issued_integrals={}

    def primitives(self,frame):
        self.source.describe(frame);a=frame.roots['q'].atlas;key=id(frame)
        if key in self.archive_cache:return self.archive_cache[key]
        row=self.saved[frame.index][frame.branch];record=row['source'];proof=row['actual_whole_phase_mixed']
        if frame.count!=self.count or record['source_family']!=frame.family or record['branch']!=frame.branch:
            raise ValueError('Issued frame must match original archive source cell and predicate')
        if record['source_level']!=frame.count or record['source_index']!=frame.index or record['exact_outer_Z_bounds']!=frame.record['exact_outer_Z_bounds']:
            raise ValueError('Original archive partition or axial domain changed')
        atlas=record['atlas']
        if atlas['source_basis_order']!=['logPstar','selected_logCstar','logL','original_logq','defined_original_logabsu_or_zero'] or not atlas['legacy_scale_record_radius_power_field_is_logabsu_power_on_this_atlas'] or not atlas['actual_radius_power_is_exact_zero_and_not_in_fifth_slot']:
            raise ValueError('Original archive q/u basis semantics changed')
        if record['exact_y_cell']!=frame.record['exact_y_cell'] or proof['defining_original_source_predicate']!=frame.record['actual_original_source_predicate']:
            raise ValueError('Original source function domain changed')
        for live,saved in zip(a.bases,record['atlas']['defining_basis'],strict=True):
            if ep(live)!=ep(interval(a.ctx,saved)):raise ValueError('Original q/u function bases differ')
        if ep(interval(a.ctx,proof['fixed_true_phase_fraction']))!=(0,1):raise ValueError('Accepted whole true inverse graph cover required')
        jets={name:MixedJet(a,{order:restore_value(a,rows[str(order)]) for order in (C0,Y,Z,YZ)})
            for name,rows in proof['complete_original_A_B_mixed'].items()}
        self.archive_cache[key]=jets;return jets

    def cell(self,index,*,N):
        N=candidate_N(N);c=self.ctx;rows=[];exports=[{key:[] for key in KEYS} for unused in range(4)]
        for branch in BRANCHES:
            frame=self.source.frame(self.count,index,branch=branch);a=frame.roots['q'].atlas;primitives=self.primitives(frame)
            part=dict(kernel=frame.kernel,roots={name:{C0:row[C0],Z:row[Z]} for name,row in frame.roots.items()},ledger=a.ledger)
            graph=CollectedDensityGraph(a,self.source.owner.scales.graph,part,
                dict(A=primitives['A'][C0],B_over_Pstar=primitives['B'][C0]),
                dict(A_Z_slow=primitives['A'][Z],B_Z_slow=primitives['B'][Z]),N)
            got=graph.outputs();old,oldZ=original_history_rows(a,frame.roots)
            groups=(got['densities'],got['density_Z'],old,oldZ)
            normalized=[]
            for j,values in enumerate(groups):
                ranges={key:export_range(a,value,derivative_unit=j==1) for key,value in values.items()}
                for key in KEYS:exports[j][key].append(ranges[key])
                normalized.append(ranges)
            rows.append(dict(branch=branch,actual_original_predicate=frame.record['actual_original_source_predicate'],
                exact_outer_Z_bounds=frame.record['exact_outer_Z_bounds'],original_mixed_archive_source_index=index,
                original_mixed_archive_filename=phase.NAME,original_mixed_receipt_filename=phase.RECEIPT,
                A_B_units='A dimensionless; mixed B is original B_over_Pstar, from B=-a*E*T1/(4*pi) with slope t0=0',
                original_A_B_C0_Z_inputs={name:{str(order):jet[order].record() for order in (C0,Z)} for name,jet in primitives.items()},
                original_graph_density_and_ordinary_Z={label:{key:value.record() for key,value in values.items()}
                    for label,values in zip(('changed_C0','changed_Z','unmodulated_C0','unmodulated_Z'),groups,strict=True)},
                exported_normalized_ranges={label:values for label,values in zip(('changed_C0','changed_Z_over_Lambda0_Lminus2','unmodulated_C0','unmodulated_Z'),normalized,strict=True)},
                changed_Z_large_source_unit_retained_formally=True,source_function_derivatives_not_derivatives_of_exported_caps=True,
                original_density_graph_operations=graph.operations,all_original_nonlinear_cross_terms_executed=True,
                common_S_Pstar_and_original_own_rate_units=five.UNITS))
        hull=lambda values:c.mpf((min(ep(v)[0] for v in values),max(ep(v)[1] for v in values)))
        union=[{key:hull(values) for key,values in group.items()} for group in exports]
        return dict(rows=rows,union=union)

    def integrate(self,*,N=160):
        N=candidate_N(N);c=self.ctx;began=time.monotonic();count=self.count
        if N<count:raise ValueError('This accepted whole-phase archive route requires N/count>=1')
        totals=[{key:c.mpf(0) for key in KEYS} for unused in range(4)];rows=[]
        with mp.workdps(c.dps+40):
            radius=self.source.owner.radius
            endpoints=[radius.evaluate(y=s.Rational(index,count),N=N) for index in range(count+1)]
            origin=endpoints[0]
            for index in range(count):
                left,right=c.mpf(index)/count,c.mpf(index+1)/count
                actual_phase=five.pressure.phase_boxes(c,origin,N,left,right)
                if len(actual_phase)!=1 or ep(actual_phase[0])!=(0,1):raise ValueError('Actual original common-N whole-period source cell required')
                local=self.cell(index,N=N);unused,masses,decay=five.masses(c,c.mpf(1)/count,left,right)
                additions=[{key:values[key]*masses[key] for key in KEYS} for values in local['union']]
                for total,add in zip(totals,additions,strict=True):
                    for key in KEYS:total[key]+=add[key]
                rows.append(dict(exact_y_cell=[str(s.Rational(index,count)),str(s.Rational(index+1,count))],
                    actual_original_common_N_phase_cover=actual_phase,actual_original_phase_endpoint_indices=[index,index+1],
                    predicate_density_records=local['rows'],branch_union_not_sum_or_duplicate_integral=True,
                    normalized_four_density_unions=local['union'],positive_own_rate_final_endpoint_masses=masses,
                    normalized_four_final_endpoint_contributions=additions,
                    whole_source_function_ranges_integrated_not_point_quadrature=True))
                if (index+1)%16==0:print('Full original O2 predicate five integrals:',index+1,'/',count,flush=True)
            changed={key:self.out.scalar(totals[0][key]) for key in KEYS}
            changedZ={key:self.Z_unit*totals[1][key] for key in KEYS}
            old={key:self.out.scalar(totals[2][key]) for key in KEYS};oldZ={key:self.out.scalar(totals[3][key]) for key in KEYS}
            inlet=original_full_inlet(self.out,c.mpf((-1,1)))
            original={key:self.out.scalar(interval(c,inlet['original_inlet'][key])*c.exp(-c.mpf(five.RATES[key]))+totals[2][key]) for key in KEYS}
            originalZ={key:self.out.scalar(interval(c,inlet['original_inlet_Z'][key])*c.exp(-c.mpf(five.RATES[key]))+totals[3][key]) for key in KEYS}
            report=dict(source_family=self.family,exact_y_window=['0','1'],exact_Z_range=['-1','1'],explicit_candidate_N=N,
                ordered_source_cells=count,own_rates=five.RATES,normalized_own_units=five.UNITS,
                source_defined_original_inlet=inlet,source_defined_inlet_contract=self.inlet_contract,
                actual_original_radius_phase_endpoints=endpoints,actual_global_phase_endpoint_indices=[0,count],
                original_pressure_datum_sha256=self.family['datum_enclosure_sha256'],
                changed_C0_integral={key:value.record() for key,value in changed.items()},
                changed_ordinary_Z_integral={key:value.record() for key,value in changedZ.items()},
                original_unmodulated_C0_integral={key:value.record() for key,value in old.items()},
                original_unmodulated_ordinary_Z_integral={key:value.record() for key,value in oldZ.items()},
                source_defined_original_histories_at_y1={key:value.record() for key,value in original.items()},
                source_defined_original_history_Z_at_y1={key:value.record() for key,value in originalZ.items()},
                changed_ordinary_Z_unit=self.Z_unit.record(),output_original_source_atlas=self.out.record(),
                all_four_normalized_integral_enclosures=totals,whole_source_density_phase_mass_records=rows,
                fixed_y_window_and_Z_independent_Duhamel_masses=True,no_extra_R_Jacobian=True,
                large_original_ordinary_Z_factor_not_materialized_or_discarded=True,
                direct_density_integration_no_phase_averaging_or_boundary_term_drop=True,
                original_pressure_separate_from_cumulative_pressure_and_incoming_histories_not_reset=True,
                actual_incoming_must_be_supplied_no_default_zero=True,
                complete_original_y_Z_and_true_phase_union_enclosed=True,
                finite_N_direct_integral_enclosure_not_sharp_terminal_norm=True,
                functional_terminal_identity_solved=False,current_whole_N_selected=False,execution_seconds=time.monotonic()-began)
        result=dict(report=report,changed=changed,changedZ=changedZ,original=original,originalZ=originalZ)
        self.issued_integrals[id(result)]=result;return result

    def apply_incoming(self,integral,*,incoming,incoming_Z,source_family,original_P0_datum_sha256):
        """Explicit native affine C1 transport; caller owns its source covers."""
        if self.issued_integrals.get(id(integral)) is not integral:
            raise ValueError('Issued same-owner original integral required')
        if source_family!=self.family or original_P0_datum_sha256!=self.family['datum_enclosure_sha256']:
            raise ValueError('Same original incoming family and separate P0 datum required')
        if integral['report']['source_family']!=self.family or integral['report']['exact_Z_range']!=['-1','1']:
            raise ValueError('Same original continuous-Z integral required')
        if set(incoming)!=set(KEYS) or set(incoming_Z)!=set(KEYS):raise ValueError('All five explicit actual incoming C0/Z covers required')
        for value in (*incoming.values(),*incoming_Z.values()):
            if type(value) is not prior.ScaledEnclosure or value.scale.bases is not self.out.bases or value.ledger is not self.out.ledger:
                raise ValueError('Same native output atlas incoming function covers required')
        c=self.ctx;out={};outZ={}
        for key,rate in five.RATES.items():
            decay=c.exp(-c.mpf(rate))
            out[key]=self.out.sum((integral['original'][key],integral['changed'][key],incoming[key]*decay))
            outZ[key]=self.out.sum((integral['originalZ'][key],integral['changedZ'][key],incoming_Z[key]*decay))
        return dict(values=out,Z=outZ,record=dict(source_family=self.family,original_P0_datum_sha256=original_P0_datum_sha256,
            explicit_incoming_source_function_covers_are_caller_obligations=True,
            actual_incoming_histories_not_defaulted_to_zero=True,original_P0_not_added_to_cumulative_p=True,
            complete_own_C0={key:value.record() for key,value in out.items()},complete_own_Z={key:value.record() for key,value in outZ.items()},
            exact_pressure_zero_rate_memory_preserved=True,functional_terminal_identity_solved=False))


def run():
    began=time.monotonic();owner=OriginalO2PredicateFiveIntegrals();integral=owner.integrate(N=160)
    report=dict(**{GATE:True},source_family=owner.family,actual_full_original_O2_five_integrals=integral['report'],
        original_predicate_finite_N_C0_Z_density_layer_installed=True,full_original_O2_changed_five_direct_integrals_installed=True,
        explicit_native_incoming_C1_transport_interface_installed=True,
        actual_all_route_incoming_histories_installed=False,sharp_O2_phase_averaging_installed=False,
        functional_terminal_identity_solved=False,all_17_chart_or_24_cell_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Full original O2 candidate-N160 five signed changed and original C0/ordinary-Z direct Duhamel integral enclosures over y[0,1], Z[-1,1], including zero and both signs. Actual source predicates, original nonlinear graph, common radius phase and nonzero original inlet. Z unit Lambda0/L^2 retained formal; explicit incoming API, no default history. Not sharp averaging, all-route incoming data, terminal controls, global N, scale recursion or full reconstruction.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Full original O2 finite-N five C0/Z direct contribution enclosures generated',flush=True);return report


if __name__=='__main__':run()

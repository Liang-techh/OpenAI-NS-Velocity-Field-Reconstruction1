"""Original O2 direct regular carrier calculus and native slow-jet norms.

Actual functions and their fixed-true-phase derivatives are retained.
Regular derivatives use raw original p2 carriers, not dstar*u/q covers.
Signed geometry keeps its actual source predicate and reciprocal factors.
"""
from dataclasses import replace
import gzip
import json
from pathlib import Path
from types import MappingProxyType
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_uniform_mean_bias as current

mixed=current.mixed;integrals=current.integrals;phase=mixed.phase;frames=mixed.frames
source=mixed.source;base,prior,ep=mixed.base,mixed.prior,mixed.ep
MixedJet=mixed.MixedJet;C0,Y,Z,YZ=mixed.ORDERS
HERE,PREFIX,sha=mixed.HERE,mixed.PREFIX,mixed.sha
NAME=PREFIX+'current_original_O2_collected_slow_jets.json.gz'
RECEIPT=PREFIX+'current_original_O2_collected_slow_jets_check.json'
GATE='original_O2_direct_regular_carrier_all_five_native_C0_y_Z_yZ_norms_installed'
UNITS={C0:(0,0,0,0,0),Y:(11,10,-1,0,0),Z:(11,10,-2,0,0),YZ:(22,20,-3,0,0)}


def direct_u_jets(a,p2,q,d,u0):
    """Ordinary source product rules with a same-source predicate C0 cover.

    The C0 cover and derivative covers enclose the original function; they
    are not asserted to form an independently selected compatible field.
    q_Z=q_yZ=0 only on O2 slope.
    """
    if not q[Z].zero or not q[YZ].zero:raise ValueError('Exact original O2 transverse q identities required')
    divide=lambda value:value.positive_divide(d,d.scale.evaluate())
    return MixedJet(a,{C0:u0,Y:divide(a.add(p2[Y]*q[C0],p2[C0]*q[Y])),
        Z:divide(p2[Z]*q[C0]),YZ:divide(a.add(p2[YZ]*q[C0],p2[Z]*q[Y]))})


def source_unit_theorem():
    y,z=s.symbols('y z',real=True);d=s.Symbol('dstar',positive=True)
    p=s.Function('p2')(y,z);q=s.Function('q')(y)
    u=p*q/d
    targets={Y:(s.diff(p,y)*q+p*s.diff(q,y))/d,Z:s.diff(p,z)*q/d,
        YZ:(s.diff(p,y,z)*q+s.diff(p,z)*s.diff(q,y))/d}
    for order,target in targets.items():assert s.simplify(s.diff(u,y,order[0],z,order[1])-target)==0
    return dict(passed=True,exact_direct_original_u_y_Z_yZ_identities=3,
        dominating_normalization_units={str(order):list(powers) for order,powers in UNITS.items()},
        unit_names={str(C0):'1',str(Y):'Lambda0/L',str(Z):'Lambda0/L^2',str(YZ):'Lambda0^2/L^3'},
        Lambda0='Pstar^11*Cstar^10',units_are_Z_only_not_differentiated_source_covers=True,
        raw_p2_C0_y_powers=[11,10,-1,0,0],raw_p2_Z_yZ_powers=[11,10,-2,0,0],
        original_dstar_offset_preserved_not_substituted_as_selected_value=True,
        signed_abs_u_lower='3/16',signed_q_inverse_implication='q^-1 <= (16/3)*Lambda0*max|normalized g|/dstar',
        necessary_signed_function_inequalities_not_claimed_on_whole_outer_rectangle=True,
        negative_u_power_upper_bound='|u|^b <= (3/16)^b for b<0 on signed predicate',
        positive_u_power_upper_bound='|u|^b <= (Lambda0*max|normalized g|*q/dstar)^b for b>0',
        transformed_monomials_are_absolute_bounds_not_redefined_source_functions=True,
        normalization_units_dominate_output_jets_not_exact_raw_derivative_powers=True,
        independent_interval_rows_not_selected_compatible_fields=True,
        no_q_floor_or_logu_as_radius=True,slow_fixed_true_phase_not_total_fast_spatial_y=True)


def regular_predicate_C0(a,raw):
    """Same-source |u|<=1/4 intersected with the direct native cover's sign.

    The huge raw scale is kept for provenance; there is no scalarization of
    Lambda0*q/dstar. This cover contains the intersection, not a field value.
    """
    if raw.scale.bases is not a.bases or raw.ledger is not a.ledger:raise ValueError('Same direct native u source required')
    if raw.zero:return a.scalar(0)
    c=a.ctx;lo,hi=ep(raw.coefficient);quarter=c.mpf(1)/4
    return a.scalar(c.mpf((0 if lo>=0 else -quarter,0 if hi<=0 else quarter)))


class CollectedO2PredicateSources(mixed.CorrelatedO2PredicateSources):
    def frame(self,count,index,*,branch):
        old=super().frame(count,index,branch=branch)
        if old.record.get('collected_original_slow_carrier_installed'):return old
        a=old.roots['q'].atlas;c=a.ctx;roots=dict(old.roots);kernel=old.kernel;u=old.u;raw_u0=None
        with mp.workdps(c.dps+40):
            raw=self.source.source_frame(count,index,Z_lower=a.bounds[0],Z_upper=a.bounds[1])
            p2=MixedJet(a,{order:frames.rebase_original_root(a,value,raw_frame=raw,source_owner=self.source)
                for order,value in raw.roots['p2'].rows.items()})
            for order,value in p2.rows.items():
                expected=(11,10,-1,0,0) if order in (C0,Y) else (11,10,-2,0,0)
                if value.scale.powers!=expected:raise ValueError('Raw original p2 source carrier unit changed')
            if branch=='regular':
                roots['p2']=p2
                logd=a.copy_interval(self.owner.scales.logs['d_star'])
                raw_u0=(p2[C0]*roots['q'][C0]).positive_divide(old.kernel.dstar,logd)
                u0=regular_predicate_C0(a,raw_u0)
                query=dict(q=roots['q'][C0],roots={name:{C0:row[C0],Z:row[Z]} for name,row in roots.items()},
                    original_u_source=u0,regular_predicate=True,signed_predicate=False)
                kernel=frames.O2PredicatePhase(query,logd)
                if kernel.geometry!='small_r_series':raise ValueError('Same original regular source predicate required')
                kernel.nu=roots['nu'][C0];kernel.n0=kernel.normalized(kernel.scalar(1),1,True)
                kernel.nt=kernel.normalized(kernel.scalar(0),1,True)
                kernel.nq=kernel.normalized(base.current.square(roots['q'][C0]),c.mpf('.5'),True)
                kernel.ntq=kernel.normalized(kernel.scalar(0),1/(2*c.sqrt(2)),False)
                u=direct_u_jets(a,p2,roots['q'],kernel.dstar,u0)
            record=dict(old.record,collected_original_slow_carrier_installed=True,
                raw_original_p2_C0_y_Z_yZ=p2.record(),actual_root_jets={name:row.record() for name,row in roots.items()},
                actual_u_mixed=u.record(),regular_raw_p2_replaces_dstar_u_over_q_for_derivatives=branch=='regular',
                regular_kernel_and_u_reissued_on_same_original_basis=branch=='regular',
                regular_direct_native_u_C0_before_predicate=None if raw_u0 is None else raw_u0.record(),
                regular_u_C0_is_cover_of_direct_source_intersected_with_named_predicate=branch=='regular',
                independently_selected_compatible_C0_derivative_fields_claimed=False,
                signed_source_kernel_and_u_unchanged=branch!='regular',
                old_regular_p2_derivative_cover=old.roots['p2'].record() if branch=='regular' else None,
                kernel_geometry=kernel.geometry_record(),source_unit_theorem=source_unit_theorem(),
                original_pressure_q_eta_nu_and_P0_unchanged=True)
            new=replace(old,roots=MappingProxyType(roots),u=u,kernel=kernel,record=record)
            del self.frames[id(old)];self.frames[id(new)]=new;self.cache[(count,index,branch)]=new;return new


def absolute_monomial_bound(owner,frame,value):
    """Same-predicate bound, with native large factors still formal.

    The output is a positive upper envelope for |value|, not a signed
    source range. Its lower endpoint is not a lower bound for |value|.
    Reciprocal q is collected only by the necessary actual signed predicate.
    """
    owner.source.describe(frame);a=frame.roots['q'].atlas;c=a.ctx
    if value.scale.bases is not a.bases or value.ledger is not a.ledger:raise ValueError('Same original native source required')
    if value.zero:return a.scalar(0)
    magnitude=base.conditioned.absolute(value);powers=list(magnitude.scale.powers)
    coefficient=magnitude.coefficient;offset=magnitude.scale.offset
    guard=c.mpf(3)/16;G=c.mpf(-ep(owner.source.g)[0]);logd=frame.kernel.dstar.scale.offset
    upower=powers[4]
    if upower:
        if frame.branch=='regular':raise ValueError('No original logabsu factor on regular atlas')
        if upower<0:coefficient*=guard**upower
        else:
            powers[0]+=11*upower;powers[1]+=10*upower;powers[3]+=upower
            offset-=logd*upower;coefficient*=G**upower
        powers[4]=0
    qpower=powers[3]
    if qpower<0:
        if frame.branch=='regular':raise ValueError('Regular source may not use reciprocal q bound')
        k=-qpower;powers[0]+=11*k;powers[1]+=10*k;powers[3]=0
        offset-=logd*k;coefficient*=(G/guard)**k
    return prior.ScaledEnclosure(prior.FormalScale(a.bases,tuple(powers),offset),coefficient,a.ledger)


def native_norm(owner,frame,value,order):
    if order not in UNITS:raise ValueError('Original slow order required')
    a=frame.roots['q'].atlas;c=a.ctx;bound=absolute_monomial_bound(owner,frame,value)
    unit=prior.ScaledEnclosure(prior.FormalScale(a.bases,UNITS[order]),1,a.ledger)
    remaining=tuple(p-u for p,u in zip(bound.scale.powers,UNITS[order],strict=True))
    if remaining[0]>0 or remaining[1]>0 or remaining[2]<0 or remaining[3]<0 or remaining[4]!=0:
        raise ValueError('Collect actual source slow unit before bounded export: '+str((order,remaining)))
    ratio=bound.positive_divide(unit,unit.scale.evaluate());finite=integrals.positive.bounded(ratio)
    if ep(finite)[0]<0 or any(not mp.isfinite(v) for v in ep(finite)):raise ArithmeticError('Finite nonnegative source norm required')
    return dict(cap=c.mpf((0,ep(finite)[1])),collected_absolute_native_bound=bound,
        exact_positive_normalization_unit=unit)


def leading_densities(a,E,V,A,B):
    """Exact signed first-order original changed-density functions and jets."""
    return dict(m=B,h=E*A,k=E*(V*A+B),e=V*B*2-E*E*A,p=E*E*A)


class OriginalO2CollectedSlowJets(phase.OriginalO2PredicateMixed):
    def __init__(self):
        self.source=CollectedO2PredicateSources();self.ctx=self.source.ctx;self.family=self.source.family
        self.hashes=dict(self.source.hashes)
        for module in (mixed,current):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted same-source original derivative and all-N mean prerequisites required')
            for name,digest in {**receipt['input_hashes'],module.RECEIPT:sha(module.RECEIPT)}.items():
                if sha(name)!=digest:raise ValueError('Original collected derivative dependency changed: '+name)
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original collected derivative closures differ')
                self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.theorem=source_unit_theorem()

    def collect(self,frame):
        self.source.describe(frame);a=frame.roots['q'].atlas;c=a.ctx
        result=self.primitive(frame,c.mpf((0,1)),phase=c.mpf((0,1)))
        # Equivalent original canceled C0 identities enclose the entire free
        # angle interval; derivative rows are still actual fixed-phi rows.
        stronger=frame.kernel.primitives(c.mpf((0,1)),'psi' if frame.branch=='regular' else 'E')
        A=MixedJet(a,{**result['A'].rows,C0:stronger['A']})
        B=MixedJet(a,{**result['B'].rows,C0:stronger['B_over_Pstar']})
        first=leading_densities(a,frame.roots['E'],frame.roots['V'],A,B)
        norm=lambda jets:{name:{str(order):native_norm(self,frame,jet[order],order) for order in mixed.ORDERS} for name,jet in jets.items()}
        primitive_norms=norm(dict(A=A,B=B));density_norms=norm(first)
        def records(norms):return {name:{order:dict(normalized_absolute_cap=row['cap'],
            collected_absolute_native_bound=row['collected_absolute_native_bound'].record(),
            exact_positive_normalization_unit=row['exact_positive_normalization_unit'].record(),
            positive_envelope_lower_endpoint_not_source_magnitude_lower_bound=True,
            bound_only_not_source_function_or_derivative_of_cap=True) for order,row in rows.items()} for name,rows in norms.items()}
        return dict(A=A,B=B,leading=first,primitive_norms=primitive_norms,density_norms=density_norms,
            record=dict(source=self.source.describe(frame),actual_full_phase_mixed_before_C0_strengthening=result['record'],
                actual_same_source_A_B_C0_y_Z_yZ=dict(A=A.record(),B=B.record()),
                native_same_source_signed_first_order_density_C0_y_Z_yZ={name:jet.record() for name,jet in first.items()},
                primitive_native_norms=records(primitive_norms),first_order_density_native_norms=records(density_norms),
                exact_y_source_cell=frame.record['exact_y_cell'],exact_outer_Z_bounds=frame.record['exact_outer_Z_bounds'],
                exact_true_phase_domain=['0','1'],original_phase_reflection_receipt=current.current.RECEIPT,
                exact_leading_density_true_phase_mean_and_slow_derivatives_zero=True,
                stronger_C0_cover_not_differentiated=True,all_signed_density_product_and_cross_terms_retained=True,
                original_radius_common_N_phase_not_yet_spatially_integrated=True,
                same_source_norms_not_point_values=True))


def run():
    began=time.monotonic();owner=OriginalO2CollectedSlowJets();count=64;rows=[]
    maxima={name:{str(order):owner.ctx.mpf(0) for order in mixed.ORDERS} for name in ('A','B',*integrals.KEYS)}
    with mp.workdps(owner.ctx.dps+40):
        for index in range(count):
            cells={}
            for branch in integrals.BRANCHES:
                frame=owner.source.frame(count,index,branch=branch);result=owner.collect(frame);cells[branch]=result['record']
                for norms in (result['primitive_norms'],result['density_norms']):
                    for name,orders in norms.items():
                        for order,row in orders.items():
                            upper=max(ep(maxima[name][order])[1],ep(row['cap'])[1]);maxima[name][order]=owner.ctx.mpf((0,upper))
            rows.append(cells)
            if (index+1)%16==0:print('Collected actual original O2 slow norms:',index+1,'/',count,flush=True)
        report=dict(**{GATE:True},source_family=owner.family,whole_original_O2_collected_slow_source_cells=rows,
            original_source_cells=count,original_predicate_frames=count*3,
            exact_outer_source_domain=dict(y=['0','1'],Z=['-1','1'],phi=['0','1']),
            original_slow_unit_theorem=owner.theorem,whole_source_normalized_slow_norms=maxima,
            same_source_five_first_order_density_functions_and_four_slow_norms_installed=True,
            regular_raw_p2_direct_derivative_calculus_installed=True,
            original_mean_bias_component_unchanged=True,
            actual_zero_mean_phase_primitive_evaluator_installed=False,
            oscillatory_spatial_integral_remainder_enclosed=False,sharp_O2_phase_averaging_installed=False,
            actual_all_route_incoming_histories_installed=False,functional_terminal_identity_solved=False,
            current_whole_N_selected=False,all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,
            **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual original O2 regular raw carrier product calculus and native same-source A/B plus signed first-order five-density C0/y/Z/yZ norms on the full predicate union. Z-only units 1,Lambda0/L,Lambda0/L^2,Lambda0^2/L^3 remain formal. Not a selected field, phase primitive evaluation, sharp spatial averaging, terminal/global N, actual all-route driver, scale recursion or corrected NS.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Original O2 collected native slow jets generated',flush=True);return report


if __name__=='__main__':run()

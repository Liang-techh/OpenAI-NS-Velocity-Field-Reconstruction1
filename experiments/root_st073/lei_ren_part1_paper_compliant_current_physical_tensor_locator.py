"""Directed physical coordinate inversion over the checked 33-region graph.

Source interval uncertainty can give several candidate regions. Their union,
not a midpoint branch choice, encloses the physical tensor. This is still a
source-function enclosure, not a converged velocity/pressure point solver.
"""
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_tensor_registry import (
    CurrentTensorRegistry,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN as REGISTRY_OPEN)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator import (
    V,SYMBOLS,DOMAINS,original_native_radius_recipes,original_radial_cover_theorem,
    implicit_physical_source_theorem,implicit_log_coordinate_map,evaluate_expression,
    constrained_core_program,intersect,nonpositive_exp)

NAME=PREFIX+'current_physical_tensor_locator.json'
RECEIPT=PREFIX+'current_physical_tensor_locator_check.json'
GATES=('current_actual_implicit_physical_coordinates_directed',
       'current_33_region_original_radius_inverse_candidate_locator_available')
# Numeric candidate cover is not yet a source-bound global physical cover
# certificate: exact positivity and all ambiguities must be independently
# audited in that later task. Native old full-view gates remain unchanged.
OPEN=REGISTRY_OPEN+('current_exact_physical_boundary_trace_selection_certified',)


class CurrentPhysicalTensorLocator:
    @source_precision
    def __init__(self,registry=None,require_checked=True):
        self.registry=registry if registry is not None else CurrentTensorRegistry()
        if type(self.registry) is not CurrentTensorRegistry or not self.registry.acceptance_loaded:
            raise ValueError('Checked one-graph full tensor registry required')
        self.ctx=self.registry.ctx;self.physical=self.registry.physical
        self.family=self.registry.family;self.source=self.registry.source;self.datum_sha=self.registry.datum_sha
        self.recipes,self.radius_program_proof=original_native_radius_recipes()
        self.cover_theorem=original_radial_cover_theorem();self.implicit_theorem=implicit_physical_source_theorem()
        self.core_program,self.core_program_proof=constrained_core_program()
        self.values,self.value_bindings=self._source_values()
        asts=SourceAST()
        for stem,method in (('current_complete_physical_assembly','__init__'),('inner_bridge_profiles','__init__'),
            ('inner_switch_profiles','__init__'),('actual_reference_restore_mixed_C4','__init__'),
            ('current_power_angular_source','__init__'),('current_steep_waiting_source','__init__')):
            asts.method(stem,method)
        # The outer parameter source is outside the compliant filename prefix.
        parameter_name='lei_ren_part1_paper_logarithmic_outer_parameters.py'
        self.hashes={**self.registry.hashes,**self.cover_theorem['input_hashes'],**self.implicit_theorem['input_hashes'],
            **self.core_program_proof['input_hashes'],**asts.hashes,parameter_name:sha(parameter_name),
            Path(__file__).name:sha(Path(__file__).name)}
        operator_name=PREFIX+'current_physical_tensor_locator_operator.py';self.hashes[operator_name]=sha(operator_name)
        self.acceptance_loaded=False;self.assert_graph()
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Physical coordinate locator admission exceeds source enclosure scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def _source_values(self):
        c=self.ctx;p=self.physical;a=SYMBOLS
        bridge=p.dispatch.provider('bridge_first');switch=p.dispatch.provider('switch_first')
        reference=p.dispatch.provider('inner_reference');reshape=p.dispatch.provider('reshape')
        outer=p.dispatch.provider('outer_power');steep=p.dispatch.provider('steep_entry')
        if endpoints(bridge.logh)!=endpoints(switch.logh):raise ValueError('Same original bridge/switch positive hb source required')
        values=dict(log_epsilon=self.registry.owners['axis'].log_epsilon,log_Ra=c.ln(bridge.bridge.r),
            h_bridge=bridge.h,h_switch=switch.h,log_Rref=p.logRref,log_Rp=p.logRp,log_Rv=p.logRv,
            mu=p.params.mu,T=reshape.reshape.T,log_gap=reference.reference.loggap,Md=c.mpf(p.params.Md),
            Tw=p.params.Tw,Lrel=outer.Lrel,Ts=steep.Ts,wait=steep.wait,log_C=p.logC,log_P=p.logP)
        for name in ('Md','T','Tw','Ts','wait'):
            if endpoints(values[name])[0]<=0:raise ValueError('Original positive radius slope guard failed: '+name)
        if not 0<endpoints(values['mu'])[0] or endpoints(values['mu'])[1]>=mp.mpf('.25'):
            raise ValueError('Original reciprocal gap monotonicity requires 0<mu<1/4')
        if endpoints(values['log_gap'])[0]<=8 or endpoints(values['Lrel'])[0]<=4:
            raise ValueError('Original reference/outer power radius widths must be positive')
        for name in ('h_bridge','h_switch'):
            if endpoints(values[name])[0]<0 or any(not mp.isfinite(v) for v in endpoints(values[name])):
                raise ValueError('Nonnegative finite enclosure of original positive hb required')
        if endpoints(c.ln(100)-values['log_Ra']-2*values['h_bridge'])[0]<=0 or endpoints(c.ln(c.mpf(110)/100)-2*values['h_switch'])[0]<=0:
            raise ValueError('Original macro bridge/switch power widths must be positive')
        if any(not mp.isfinite(v) for v in endpoints(bridge.logh)):
            raise ValueError('Exact finite log source establishes hb>0 before its cap enclosure')
        bindings=dict(log_epsilon='registry.owners[axis].log_epsilon: same original core epsilon',
            log_Ra='physical.dispatch.provider(bridge_first).bridge.r: exact original Ra=4epsilon',
            h_bridge='bridge.h bounds only; exact source exp(bridge.logh)',
            h_switch='switch.h bounds only; same exact positive source as bridge.logh',
            log_Rref='physical.logRref=log110+10(logC+logP)',log_Rp='physical.logRp=logRref+logP+1+Tw',
            log_Rv='physical.logRv=logRp+13/mu',mu='physical.params.mu',T='physical reshape provider.reshape.T',
            log_gap='physical reference provider.reference.loggap=10(logC+logP)-T',Md='physical.params.Md',
            Tw='physical.params.Tw',Lrel='physical outer power provider.Lrel',Ts='physical steep entry provider.Ts',
            wait='physical steep entry provider.wait: exact positive waiting root enclosure',
            log_C='physical.logC: same original nonlinear source',log_P='physical.logP: original exp(Md)+11')
        bindings['exact_positive_hb_log']=bridge.logh
        return {a[name]:value for name,value in values.items()},bindings

    def assert_graph(self):
        self.registry.assert_graph()
        if self.physical is not self.registry.physical or self.ctx is not self.registry.ctx or not self.registry.acceptance_loaded:
            raise ValueError('Foreign physical locator graph')
        if self.recipes is not original_native_radius_recipes()[0] or self.core_program is not constrained_core_program()[0]:
            raise ValueError('Foreign original radius/Cartesian source program')
        live,bindings=self._source_values()
        if bindings.keys()!=self.value_bindings.keys() or any(endpoints(live[key])!=endpoints(value) for key,value in self.values.items()):
            raise ValueError('Current original source radius parameters changed')

    def _number(self,expr):return evaluate_expression(self.ctx,s.sympify(expr),self.values)

    @source_precision
    def original_log_radius(self,region,coordinate):
        self.assert_graph();coordinate=self._number(coordinate) if isinstance(coordinate,s.Basic) else self.ctx.mpf(coordinate)
        domain=self._number(DOMAINS[region][0]),self._number(DOMAINS[region][1])
        if endpoints(coordinate)[0]<endpoints(domain[0])[0] or endpoints(coordinate)[1]>endpoints(domain[1])[1]:
            raise ValueError('Native radius coordinate outside its original source domain')
        return evaluate_expression(self.ctx,self.recipes[region],{**self.values,V:coordinate})

    def _inverse_candidate(self,region,logR):
        c=self.ctx;expr=self.recipes[region];lo,hi=DOMAINS[region]
        domain=c.mpf([endpoints(self._number(lo))[0],endpoints(self._number(hi))[1]])
        radius_range=c.mpf([endpoints(self.original_log_radius(region,lo))[0],
            endpoints(self.original_log_radius(region,hi))[1]])
        overlap=intersect(c,logR,radius_range)
        if overlap is None:return None
        status='directed_original_radius_inverse'
        if region in ('core_positive_radius','actual_patch'):
            base=s.simplify(expr-s.log(V));logcoord=overlap-self._number(base)
            log_domain=c.ln(domain);logcoord=intersect(c,logcoord,log_domain)
            if logcoord is None:return None
            # Core rho can be arbitrarily small; evaluate only exp(logrho)
            # after shifting by its finite maximum. The source stays exact.
            shift=endpoints(logcoord)[1]
            coordinate=c.exp(c.mpf(shift))*nonpositive_exp(c,logcoord-shift)
        elif region=='O2_axial':
            base=s.simplify(expr-s.exp(SYMBOLS['Md']*V))
            exp_coordinate=intersect(c,overlap-self._number(base),c.mpf([1,endpoints(c.exp(self.values[SYMBOLS['Md']]))[1]]))
            if exp_coordinate is None:return None
            coordinate=c.ln(exp_coordinate)/self.values[SYMBOLS['Md']]
        else:
            base=s.simplify(expr.subs(V,0));width=s.simplify(s.diff(expr,V))
            if V in width.free_symbols:raise ValueError('Unsupported original nonlinear radius inverse')
            width_box=self._number(width)
            if endpoints(width_box)[0]<=0:
                coordinate=domain;status='positive_source_width_below_numeric_resolution'
            else:coordinate=(overlap-self._number(base))/width_box
        coordinate=intersect(c,coordinate,domain)
        if coordinate is None:return None
        cap_participates=bool({SYMBOLS['h_bridge'],SYMBOLS['h_switch']} & expr.free_symbols)
        return dict(region=region,native_coordinate_enclosure=coordinate,source_logR_overlap=overlap,
            source_logR_range=radius_range,exact_original_radius_source=str(expr),
            inverse_status=status,native_coordinate_width=c.mpf(endpoints(coordinate)[1]-endpoints(coordinate)[0]),
            bound_used_as_source=False,numeric_radius_range_and_inverse_are_enclosures_only=True,
            hb_cap_participates_in_numeric_candidate_enclosure=cap_participates,
            exact_hb_source_expression='exp(exact_positive_hb_log)' if cap_participates else None,
            exact_positive_hb_log=self.value_bindings['exact_positive_hb_log'] if cap_participates else None)

    @source_precision
    def locate_log_radius(self,log_r_phys,z_phys,log_tau,theta=0,viscosity=1):
        self.assert_graph();c=self.ctx;lr=c.mpf(log_r_phys);angle=c.mpf(theta)
        if any(mp.isnan(v) or v==mp.inf for v in endpoints(lr)) or any(not mp.isfinite(v) for v in endpoints(angle)):
            raise ValueError('Finite log radius (or -inf at axis) and angle required')
        mapping=implicit_log_coordinate_map(c,z_phys,log_tau,viscosity,self.physical.delta)
        logR=2*lr-c.ln(2)-c.ln(mapping['physical_viscosity'])-2*mapping['actual_log_lambda']
        if endpoints(lr)==(-mp.inf,-mp.inf):
            candidates=[dict(region='core_positive_radius',native_coordinate_enclosure=c.mpf(0),
                inverse_status='exact_nonsingular_axis',source_logR_overlap=logR,native_coordinate_width=c.mpf(0),bound_used_as_source=False)]
        else:
            candidates=[candidate for region in self.recipes if (candidate:=self._inverse_candidate(region,logR)) is not None]
        if not candidates:raise ValueError('No original source region encloses the physical coordinate')
        return dict(**mapping,physical_log_radius=lr,theta=angle,source_logR_enclosure=logR,candidates=candidates,
            candidate_region_count=len(candidates),region_ambiguity_retained=len(candidates)>1,
            defining_field_never_replaced_by_midpoint_or_radius_cap=True,
            source_parameter_and_hb_cap_enclosures_used_only_for_candidate_ranges=True,
            boundary_policy='Union of all directed source candidates; exact trace selection remains open',
            actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            **dict.fromkeys(OPEN,False))

    def _refined_groups(self,view,q):
        c=self.ctx;result={}
        for label,parts in view['canonical_signed_component_groups'].items():
            result[label]=[]
            for row in parts:
                refined=dict(row);factor=row['physical_lambda_exponent']*q
                magnitude=max(abs(v) for v in endpoints(row['signed_coefficient']))
                upper=None if not magnitude else (sum(row['actual_source_log_parts'].values(),c.mpf(0))+
                    row['radial_log_prefactor']+factor+row['physical_viscosity_log_prefactor']+c.ln(c.mpf(magnitude)))
                refined.update(actual_log_lambda=q,actual_lambda_log_prefactor=factor,
                    actual_lambda_log_absolute_upper=upper,
                    original_native_sqrt_tau_log_absolute_upper=row['log_absolute_upper'],
                    actual_log_tau_retained=True,exact_lambda_source='lambda^2*(1-Z^2)=tau')
                result[label].append(refined)
        return result

    @source_precision
    def tensor_from_location(self,location):
        self.assert_graph();c=self.ctx;views=[]
        # Only locations from this live source are admitted; recompute the
        # cheap coordinate map so forged/changed candidate lists cannot route.
        expected=self.locate_log_radius(location['physical_log_radius'],location['physical_z'],location['requested_log_tau'],
            location['theta'],location['physical_viscosity'])
        if encode(pack(expected))!=encode(pack(location)):raise ValueError('Foreign or changed physical locator result')
        for candidate in location['candidates']:
            region=candidate['region'];v=candidate['native_coordinate_enclosure']
            if region=='core_positive_radius':
                root=c.sqrt(2*v);X=root*c.cos(location['theta']);Y=root*c.sin(location['theta'])
                source=self.core_program(self.registry.owners['axis'],location['Z'],X,Y,
                    location['requested_log_tau'],location['physical_viscosity'],rho_source=v)
                view=self.registry._view(region,source,'Cartesian_core')
                view['joint_polar_Cartesian_source_constraint']=self.core_program_proof
            elif region=='heat_exterior':
                # The original admitted Gamma T/E are identically zero on
                # the entire tail. Avoid evaluating exp of an astronomical
                # offset from a wide absolute-radius parameter enclosure.
                view=self.registry.exterior(location['Z'],location['requested_log_tau'],location['theta'],location['physical_viscosity'])
                view['same_source_unbounded_Gamma_zero_tensor_used']=True
            else:view=self.registry.native(region,location['Z'],v,location['requested_log_tau'],location['theta'],location['physical_viscosity'])
            views.append(dict(candidate=candidate,original_full_native_view=view,
                actual_lambda_signed_component_groups=self._refined_groups(view,location['actual_log_lambda'])))
        return dict(location=location,candidate_full_tensor_views=views,
            return_kind='Union of directed complete original source-function tensor enclosures',
            original_source_sectors_and_actual_lambda_factors_retained=True,
            arbitrary_candidate_or_midpoint_not_selected=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def cartesian(self,x_phys,y_phys,z_phys,time,viscosity=1,terminal_time=1,tensors=True):
        c=self.ctx;x=c.mpf(x_phys);y=c.mpf(y_phys);t=c.mpf(time);terminal=c.mpf(terminal_time)
        if any(not mp.isfinite(v) for v in endpoints(x)+endpoints(y)+endpoints(t)+endpoints(terminal)):
            raise ValueError('Finite physical Cartesian coordinates and times required')
        tau=terminal-t
        if endpoints(tau)[0]<=0:raise ValueError('Strictly positive physical time-to-terminal required')
        radius=c.sqrt(x**2+y**2);axis=endpoints(radius)==(mp.mpf(0),mp.mpf(0))
        angle=c.mpf(0) if axis else c.atan2(y,x)
        location=self.locate_log_radius(c.mpf('-inf') if axis else c.ln(radius),z_phys,c.ln(tau),angle,viscosity)
        return self.tensor_from_location(location) if tensors else location

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            original_radius_source_recipes={name:str(expr) for name,expr in self.recipes.items()},
            current_original_radius_source_parameter_bindings=self.value_bindings,
            original_exact_radial_endpoint_theorem=self.cover_theorem,actual_implicit_physical_source_theorem=self.implicit_theorem,
            exact_joint_core_Cartesian_source=self.core_program_proof,
            actual_current_tensor_regions_available=list(self.recipes),
            scope='Directed actual implicit coordinate map and all 33 original radius inverse candidate formulas on the checked full tensor graph. Finite positive times and constant nu>0; candidate union preserves uncertainty. Exact physical seam selection/global cover, resolved point u/v/w/p, cone/lift, corrected NS, energy and completed recursion/flatness remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPhysicalTensorLocator(require_checked=False)
    field.assert_graph();result=field.manifest()
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current physical tensor locator: directed implicit lambda and 33 original radius inverse candidate sources',flush=True)
    return result


if __name__=='__main__':run()

"""Five same-current pulse joins in raw mixed source units.

Formal reciprocal boundaries are exact graph substitutions. Supplemental
native routes cover rounding strips without redefining those boundaries.
This does not install a uniform norm, numerical point or physical time field.
"""
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time
from types import CodeType

import lei_ren_part1_paper_compliant_current_original_Rp_mixed_transport as mixed
import lei_ren_part1_paper_compliant_current_pulse_interfaces as interfaces
from lei_ren_part1_paper_compliant_axial_pulse_field import CompliantAxialPulseField
from lei_ren_part1_paper_compliant_pulse_high_jets import CompliantPulseHighJets
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import CompliantPulseMixedC4
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

pulse=mixed.pulse
HERE,PREFIX,sha=mixed.HERE,mixed.PREFIX,mixed.sha
NAME=PREFIX+'current_original_Rp_pulse_mixed_seams.json.gz'
RECEIPT=PREFIX+'current_original_Rp_pulse_mixed_seams_check.json'
GATES=('current_original_Rp_five_pulse_raw_mixed4_function_seams_installed',
       'current_original_Rp_formal_reciprocal_mixed_source_views_installed',
       'current_original_Rp_pulse_rounded_coverage_supplemental_routes_installed')
OPEN=mixed.OPEN
SEAMS=('entrance_main','main_exit','exit_gap','gap_coordinate','gap_end')
OVERLAPS={'entrance_main_overlap':('pulse_main',Fraction(19999,1000000),Fraction(1,50)),
          'gap_coordinate_overlap':('pulse_gap',Fraction(12),Fraction(120001,10000))}


def semantic_code(method):
    """Include nested code semantics, excluding source filename/line tables."""
    def signature(code):
        def constant(value):
            if isinstance(value,CodeType):return signature(value)
            if isinstance(value,tuple):return tuple(constant(v) for v in value)
            return value
        return (code.co_code,tuple(constant(v) for v in code.co_consts),code.co_names,
            code.co_varnames,code.co_freevars,code.co_cellvars,code.co_argcount,
            code.co_posonlyargcount,code.co_kwonlyargcount,code.co_flags,
            code.co_stacksize,code.co_exceptiontable)
    return signature(getattr(method,'__wrapped__',method).__code__)


def selected_energy_branch_theorem():
    """The actual main energy switch preserves the whole defining function."""
    s=interfaces.s
    xi,mu=s.symbols('xi mu',positive=True);z=s.Symbol('Z',real=True)
    source=lambda name:s.Function(name)(z)
    a,e0,F,end=[source(name) for name in ('selected_ap','incoming_energy','complete_future_half','exact_end_energy')]
    K=s.Symbol('original_total_gp_energy',positive=True)
    tail=s.Function('original_remaining_gp_energy')(xi)
    partial=K-tail;D=13-xi
    target=(1-s.exp(-26))/4-mu*e0+mu*s.exp(-26)*(F-end)
    forward=s.exp(2*xi)*(e0+a*a*partial/mu-(1-s.exp(-2*xi))/(4*mu))
    backward=F*s.exp(-2*D)+(1-s.exp(-2*D))/(4*mu)-a*a*s.exp(2*xi)*tail/mu-end*s.exp(-2*D)
    difference=s.simplify(s.expand_power_exp(s.expand(forward-backward).subs(a*a*K,target)))
    if difference!=0:raise ArithmeticError('Selected forward/backward energy function differs')
    specs=(('main',"K=gp_energy(c,xi,self.selection.K)"),
        ('main',"emu=(e0*self.mu+ap*ap*K-decay_integral(c,2,xi)/2)*c.exp(2*xi)"),
        ('main',"e=(future*c.exp(-2*D)+decay_integral(c,2,D)/(2*self.mu)-ap*ap*(c.exp(2*xi)*gp_future_energy(c,xi)/self.mu)-end_energy*(c.exp(-2*D)*self.E2cap))"),
        ('gp_energy',"result+=c.exp(-2*a)*decay_integral(c,2,xi/cells)*shape**2"),
        ('gp_future_energy',"total+=c.exp(-2*a)*decay_integral(c,2,length/cells)*shape**2"))
    bindings={method+':'+text:interfaces.statement('compliant_axial_pulse_field',method,text) for method,text in specs}
    return dict(passed=True,forward_backward_energy_difference_identically_zero=True,
        source_integral_partition='K_partial(xi)+K_future(xi)=integral_0^11 exp(-2v)*gp(v)^2 dv=K',
        complete_future_half_and_exact_nonzero_end_energy_retained=True,
        equality_holds_as_function_before_y_Z_differentiation=True,
        branch_at_xi10_does_not_create_a_derivative_jump=True,production_bindings=bindings)


class CurrentOriginalRpPulseMixedSeams:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else mixed.CurrentOriginalRpMixedTransport()
        if type(self.before) is not mixed.CurrentOriginalRpMixedTransport or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current fifteen-chart mixed owner required')
        self.owner=self.before.owner;self.selected=self.owner.selected
        self.pulse=self.selected.pulse;self.radius=self.owner.radius
        self.graph=self.before.graph;self.ctx=self.before.ctx;self.family_record=self.before.family_record
        self.hashes=dict(self.before.hashes)
        for name in (mixed.NAME,mixed.RECEIPT,Path(__file__).name,Path(interfaces.__file__).name):
            self.hashes[name]=sha(name)
        # The generic theorem uses arbitrary incoming/selected Z functions;
        # its old numerical physical owner is never constructed or consumed.
        self.production_bindings=interfaces.pulse_source_bindings()
        self.canonical=interfaces.functional_identities()
        self.theorem=interfaces.mixed_source_theorem()
        self.energy_branch=selected_energy_branch_theorem()
        proof_name=PREFIX+'pulse_interface_certificate.json'
        proof=json.loads((HERE/proof_name).read_bytes())
        if (proof['actual_five_defect_family_sha256'],proof['implicit_source_sha256'])!=(
                self.selected.family,self.selected.source):raise ValueError('Pulse theorem family differs')
        if proof['source_bound_functional_pulse_identities']!=self.canonical or not all(self.canonical.values()):
            raise ValueError('Unchanged exact uncapped pulse theorem required')
        flags=('all_passed','exact_uncapped_selected_sources_used','numerical_caps_used_only_as_enclosures',
            'exact_functional_main_gap_and_gap_end_identities_certified',
            'internal_O4_chart_function_equality_not_sample_overlap','functional_axial_derivatives_through5_identified')
        self.uncapped_certificate_flags={key:proof[key] for key in flags}
        if not all(self.uncapped_certificate_flags.values()):raise ValueError('Accepted uncapped pulse source certificate required')
        current_name=PREFIX+'current_selected_energy_source_check.json'
        selected_proof=json.loads((HERE/current_name).read_bytes())
        fields=('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')
        if dict((key,selected_proof[key]) for key in fields)!=self.family_record:
            raise ValueError('Current selected equation certificate family differs')
        self.selection_equation_proof=pulse.radius.post.selected.source_proof(self.selected.seed)
        if self.selection_equation_proof!=selected_proof['current_future_selected_source_proof'] or \
                self.selection_equation_proof!=self.selected.current_selection_source_proof:
            raise ValueError('Actual current positive selection equation source differs')
        required=('passed','checked_defining_parameter_identity_consumed',
            'same_current_positive_quadratic_gives_forward_full_future_half',
            'nonzero_formal_end_energy_and_incoming_terms_retained',
            'original_source_integrals_not_replaced_by_enclosure_values')
        if not all(self.selection_equation_proof[key] for key in required) or not all(selected_proof[key]
                for key in ('all_passed','current_zero_meridional_terminal_histories','current_positive_terminal_full_future_half')):
            raise ValueError('Accepted exact current terminal/uncapped selection equations required')
        pulse.radius.post.selected.inlet.add_hashes(self.hashes,proof['input_hashes'])
        self.hashes[proof_name]=sha(proof_name)
        pulse.radius.post.selected.inlet.add_hashes(self.hashes,selected_proof['input_hashes'])
        self.hashes[current_name]=sha(current_name)
        self.acceptance_loaded=False;self.call_trace=[]
        self.assert_graph();self.coverage()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed'] or not all(record[k] for k in GATES) or any(record[k] for k in OPEN):
                raise ValueError('Current pulse mixed-seam receipt/scope differs')
            if record['source_family']!=self.family_record:raise ValueError('Pulse mixed-seam family differs')
            pulse.radius.post.selected.inlet.add_hashes(self.hashes,record['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        p=self.pulse;s=self.selected;ends=pulse.radius.post.selected.inlet.endpoints
        result=dict(accepted_current_mixed_source=self.before.acceptance_loaded and all(self.before.assert_graph().values()),
            same_current_selected_future_and_independent_P0=all(s.assert_graph().values())
                and p.selection.high is p.high and p.high.select.__self__ is s.fifth
                and p.high.energy.future.__self__ is s.fifth and s.datum is s.flatten.inlet.datum,
            same_source_graph=self.graph is self.owner.graph is self.radius.graph,
            same_six_current_providers=all(self.owner.post.provider(chart) is p for chart in pulse.PULSE),
            original_entrance=p.entrance.__func__ is CompliantAxialPulseField.entrance,
            original_main=p.main.__func__ is CompliantPulseHighJets.main,
            original_gap=p.gap.__func__ is CompliantAxialPulseField.gap,
            original_gap_end=p.gap_from_end.__func__ is CompliantAxialPulseField.gap_from_end,
            original_reduced_gap=p._gap.__func__ is CompliantPulseHighJets._gap,
            original_end=p.end.__func__ is CompliantPulseHighJets.end,
            original_mixed_source_recovery=p._high_packet.__func__ is CompliantPulseMixedC4._high_packet,
            current_rate_parameter_copies=ends(p.mu)==ends(s.future.mu)
                and ends(p.rate)==ends(1-p.mu) and ends(p.prate)==ends(1+2*p.mu),
            current_exact_end_scale_copy=ends(p.logE)==ends(s.amplitude.log_end_scale)
                and ends(s.amplitude.log_end_scale)==ends(p.pulse.logscale),
            positive_caps_only_enclose_exact_end_source=all(ends(getattr(p,name))[0]==0
                and ends(getattr(p,name))[1]>0
                and ends(self.ctx.ln(self.ctx.mpf(ends(getattr(p,name))[1]))-degree*p.logE)[0]>0
                for name,degree in (('Ecap',1),('E2cap',2))),
            current_positive_selection_source=s.current_selection_source_proof['passed'],
            accepted_uncapped_selection_and_terminal_equation_certificate=all(self.uncapped_certificate_flags.values())
                and self.selection_equation_proof==s.current_selection_source_proof,
            arbitrary_source_five_seam_theorem=self.theorem['passed'],
            actual_main_energy_representation_switch=self.energy_branch['passed'])
        if not all(result.values()):raise ValueError('Current pulse mixed seam graph differs: '+str(result))
        return result

    @source_precision
    def coverage(self):
        """Legal numerical strips supplement the exact function endpoints."""
        c=self.ctx;ends=pulse.radius.post.selected.inlet.endpoints
        entrance_upper=ends(c.mpf('1/50')/self.pulse.mu)[0]
        entrance_image=self.pulse.mu*c.mpf(entrance_upper)
        gap_start=-ends(1/self.pulse.mu)[0]
        gap_image=13+self.pulse.mu*c.mpf(gap_start)
        if ends(entrance_image)[0]<ends(c.mpf('19999/1000000'))[1]:
            raise ArithmeticError('Supplemental entrance strip misses legal entrance coverage')
        if ends(c.mpf(gap_start)+1/self.pulse.mu)[0]<0 or ends(gap_image)[1]>ends(c.mpf('120001/10000'))[0]:
            raise ArithmeticError('Supplemental gap strip misses legal end coverage')
        return dict(legal_entrance_upper_bound=entrance_upper,legal_entrance_main_image=entrance_image,
            legal_gap_end_start_bound=gap_start,legal_gap_end_main_image=gap_image,
            supplemental_source_domains={key:[str(a),str(b)] for key,(_,a,b) in OVERLAPS.items()},
            exact_entrance_boundary='1/(50*mu)',exact_gap_end_boundary='-1/mu',
            rounded_bounds_only_choose_legal_coverage=True,
            supplemental_routes_do_not_redefine_function_endpoints=True,
            uniform_derivative_norm_admission=False)

    def geometry(self,chart,coordinate_node,definition):
        actual=self.radius._map(chart,coordinate_node)
        return dict(chart=chart,exact_native_coordinate_function=coordinate_node.node,
            exact_native_coordinate_definition=definition,
            **{key:value.node for key,value in actual.items()},
            absolute_radius_is_a_lazy_function=True,numerical_absolute_radius_point_value_installed=False)

    def from_packet(self,chart,Z,packet,geometry,scale_chart,scale_coordinate,evidence):
        """Use the admitted algebra on an exact source substitution packet."""
        u=packet['inlet_Utheta_over_Pstar_Taylor'];c=self.ctx
        jets=dict(Mz=u*packet['Mz_over_R_Utheta'],
            Mtheta=u*packet['Mtheta_over_sqrt2_R_3half_Utheta']*c.sqrt(2),
            Mtheta_z=u*u*packet['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']*c.sqrt(2),
            Mztheta=u*u*packet['Mztheta_over_R_Utheta_squared'],
            Mp=packet['pressure']['Mp_over_Pstar_squared'],P0=packet['pressure']['P0_over_Pstar_squared'],
            pressure=packet['pressure']['P_over_Pstar_squared'],Utheta=u,Uz=u*packet['Uz_over_Utheta'],
            Ur=u.truncate(4)*packet['Ur_over_sqrt_R_over_2_Utheta']/c.sqrt(2))
        base={key:self.owner.factor(key,scale_chart,scale_coordinate,pulse.POWERS[key],jet) for key,jet in jets.items()}
        velocities,velocity_evidence=self.before.velocity_rows(chart,Z,None,dict(source_packet=packet))
        histories=self.before.histories(velocities);logrows={};nativerows={}
        for key,value in base.items():
            logrows[key]={};nativerows[key]={}
            for k in range(5):
                row=value.coefficients if k==0 else velocities[key][k] if key in velocities else histories[k][key]
                pp=(0,2,2) if k and key in ('Mp','pressure') else value.powers
                for j in range(5-k):
                    coefficient=mixed.IntervalTaylor.constant(c,row[j]*math.factorial(j),0)
                    logrows[key]['y%d_Z%d'%(k,j)]=self.before.factor(key,scale_chart,scale_coordinate,geometry,pp,coefficient,k,j,False)
                    nativerows[key]['n%d_Z%d'%(k,j)]=self.before.factor(key,scale_chart,scale_coordinate,geometry,pp,coefficient,k,j,True)
        self.call_trace.append(dict(chart=chart,source_substitution=evidence,same_selected_provider=True))
        return dict(chart=chart,Z=packet['Z'],geometry=geometry,original_factorized_values=base,
            log_radius_mixed_rows=logrows,native_coordinate_mixed_rows=nativerows,
            velocity_radial_source_rows=velocities,original_forward_source_packet=packet,
            velocity_source_evidence=velocity_evidence,exact_source_substitution=evidence,
            source_endpoint_is_not_a_rounded_coverage_bound=True)

    @source_precision
    def supplemental(self,name,Z,coordinate):
        self.assert_graph()
        if name not in OVERLAPS:raise ValueError('Requested supplemental current source route required')
        chart,lower,upper=OVERLAPS[name];v=pulse.radius.exact_coordinate(coordinate)
        if not lower<=v<=upper:raise ValueError('Supplemental native xi domain required')
        packet=(self.pulse.main(Z,self.ctx.mpf(v.numerator)/v.denominator) if chart=='pulse_main'
            else self.pulse.gap(Z,self.ctx.mpf(v.numerator)/v.denominator))
        geometry=self.geometry(chart,self.graph.constant(v),str(v))
        return self.from_packet(chart,Z,packet,geometry,chart,v,'legal supplemental '+name)

    @source_precision
    def evaluate(self,name,Z):
        self.assert_graph();g=self.graph;mu=self.radius.functions['mu']
        if name not in SEAMS:raise ValueError('One of the five current pulse seams required')
        if name=='entrance_main':
            right=self.before.evaluate('pulse_main',Z,'1/50')
            node=g.quotient(g.constant('1/50'),mu,'accepted current positive mu')
            geometry=self.geometry('pulse_entrance',node,'1/(50*mu)')
            # entrance calls main(Z,mu*t,entrance_t=t). Reduce mu*t=1/50
            # as an exact function before numeric evaluation. The actual
            # main packet is the exact endpoint source observation, not a
            # second calculation with a rounded t and widened mu*t.
            left=self.from_packet('pulse_entrance',Z,right['original_forward_source_packet'],geometry,
                'pulse_main','1/50','proved entrance t=1/(50*mu) -> main xi=1/50')
        elif name=='gap_coordinate':
            left=self.before.evaluate('pulse_gap',Z,12)
            node=g.neg(g.quotient(g.one,mu,'accepted current positive mu'))
            geometry=self.geometry('pulse_gap_end',node,'-1/mu')
            p=self.pulse;c=self.ctx;L=p.pulse.rows['saddle_L'];u0=p.pulse.rows['saddle_u0']
            finite=-3*L+2*c.ln(u0)-c.ln(6)/2-2*c.ln(p.mu)
            packet=p._gap(Z,c.mpf(1),-1/(2*p.mu)+finite,
                dict(kind='gap_end',offset_from_Rv=-1/p.mu),dict(inverse_mu_term=-12*p.prate/p.mu))
            right=self.from_packet('pulse_gap_end',Z,packet,geometry,'pulse_gap',12,
                'proved s=-1/mu -> D=1, xi=12 before numeric arithmetic')
        else:
            l,r,v={'main_exit':('pulse_main','pulse_exit',10),
                'exit_gap':('pulse_exit','pulse_gap',11),'gap_end':('pulse_gap_end','pulse_end',-4)}[name]
            left=self.before.evaluate(l,Z,v);right=self.before.evaluate(r,Z,v)
        return dict(seam=name,Z=self.ctx.mpf(Z),left=left,right=right,
            exact_function_join_from_arbitrary_current_source_theorem=True,
            native_rows_compare_after_exact_Jacobian_conversion=True,
            interval_overlap_is_only_a_consistency_diagnostic=True,
            uniform_global_physical_time_and_point_admission=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


def report(view):
    return {**{k:v for k,v in view.items() if k not in ('left','right')},
        'left':mixed.view_report(view['left']),'right':mixed.view_report(view['right'])}


@source_precision
def run(before=None,observed_owner=None,observed_views=None,observed_strips=None):
    began=time.monotonic();owner=CurrentOriginalRpPulseMixedSeams(before,require_checked=False);views={}
    if observed_owner is not None:
        # A stricter admission guard can refresh provenance without repeating
        # unchanged source observations. Verify every output-producing method,
        # the same live graph and every original defining byte dependency.
        if observed_owner.before is not before or observed_owner.graph is not owner.graph:
            raise ValueError('Same accepted typed source owner required for observation reuse')
        for method in ('geometry','from_packet','supplemental','evaluate'):
            old=getattr(type(observed_owner),method);new=getattr(type(owner),method)
            if semantic_code(old)!=semantic_code(new):
                raise ValueError('Output-producing source method changed: '+method)
        differing={name for name in set(observed_owner.hashes)|set(owner.hashes)
            if observed_owner.hashes.get(name)!=owner.hashes.get(name)}
        if differing-{Path(__file__).name}:raise ValueError('Original source observation dependencies changed: '+str(differing))
        if set(observed_views)!=set(SEAMS) or set(observed_strips)!=set(OVERLAPS):
            raise ValueError('Complete typed source observations required')
        views=observed_views;strips=observed_strips;owner.call_trace=list(observed_owner.call_trace)
    else:
        for name in SEAMS:
            views[name]=owner.evaluate(name,'.521');print('Actual current raw pulse seam:',name,flush=True)
        strips={name:owner.supplemental(name,'.521',(a+b)/2) for name,(_,a,b) in OVERLAPS.items()}
    result=dict(source_family=owner.family_record,candidate_current_pulse_raw_mixed_seams_constructed=True,
        actual_source_graph=owner.assert_graph(),current_production_bindings=owner.production_bindings,
        source_bound_canonical_theorem=owner.canonical,arbitrary_current_source_mixed4_theorem=owner.theorem,
        actual_selected_energy_branch_theorem=owner.energy_branch,
        accepted_uncapped_source_certificate_flags=owner.uncapped_certificate_flags,
        actual_current_positive_selection_terminal_equation_proof=owner.selection_equation_proof,
        actual_five_pulse_mixed_seam_views={k:report(v) for k,v in views.items()},
        actual_supplemental_coverage_views={k:mixed.view_report(v) for k,v in strips.items()},
        rounded_coverage_ledger=owner.coverage(),actual_source_call_trace=owner.call_trace,
        fresh_Z='.521',source_scope='Five raw logR/Z mixed4 function joins, exact formal source endpoint substitutions and legal supplemental coverage; no uniform/global/physical point admission',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(pulse.raw.packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    return owner,views,strips


if __name__=='__main__':run()

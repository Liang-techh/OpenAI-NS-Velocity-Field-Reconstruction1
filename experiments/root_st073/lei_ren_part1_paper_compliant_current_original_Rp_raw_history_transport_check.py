"""Independent units, cumulative transport and current live-call evidence."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_raw_history_transport as current
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as generic
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

PAIRS=(('flatten_exit','power_inlet'),('power_exit','angular_inlet'),
    ('angular_exit','entry_inlet'),('entry_exit','steep_inlet'),
    ('steep_exit','exit_inlet'),('exit_exit','waiting_inlet'),
    ('waiting_exit','collar_inlet'),('collar_exit','exterior_inlet'))


def same(a,b,label):assert s.cancel(s.expand(a-b))==0,label


def interpreter(owner,bindings=None):
    field=SimpleNamespace(graph=owner.graph,parameters=owner.radius.parameters,
        contracts=owner.radius.frame.bridge.leading.contracts)
    return current.radius.RadiusInterpreter(field,False,bindings)


def exact_scale_and_cumulative_proof(owner):
    # Read the complete inlet function, rather than the native numerical U box.
    q0,constants=owner.post.before.inlet.exact.current_constants()
    Tw=s.Symbol('Tw',positive=True)
    q=interpreter(owner,{owner.radius.frame.functions['Tw'].node:Tw})
    same(q.at(owner.U0),constants['U'],'same actual current U0 defining function')
    assert q.z not in q.at(owner.U0).free_symbols
    same(sum(q.at(v) for v in owner.logEv0_parts.values()),
        q.logP+s.log(constants['U'])-s.Rational(13,2)/q.mu-13,'exact positive Ev0 function')
    native_flat=owner.post.before.flatten
    bindings=current.radius.post.downstream.inlet.identity.ast_assignments(
        'current_pulse_flatten_source','CurrentFlattenMixedC4','__init__',{
            'self.U':"self.inlet.constants['U']",
            'self.logEv2_parts':"dict(inlet_log=2*c.ln(self.U),inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))"})
    assert native_flat.U is native_flat.inlet.constants['U']
    assert current.radius.post.selected.inlet.endpoints(native_flat.logEv2_parts['finite_offset'])==(-26,-26)
    # Exact source scales are inspected under algebraic bindings; astronomical
    # absolute logarithms are never added numerically to a finite displacement.
    A,L,U,P,mu=s.symbols('absolute_logRp local_offset logU0 logP mu',real=True)
    q=interpreter(owner,{owner.radius.logRp.node:A,owner.logU0.node:U,
        owner.logP.node:P,owner.mu.node:mu})
    r,e,p=s.symbols('r e p',real=True)
    expected=dict(Mz=(1,1,0),Mtheta=(s.Rational(3,2),1,0),
        Mtheta_z=(s.Rational(3,2),2,0),Mztheta=(1,2,0),
        Mp=(0,0,2),P0=(0,0,2),pressure=(0,0,2),Utheta=(0,1,0))
    identities=0;cancelled=0
    for chart in current.CHARTS:
        # Bind the nonconstant chart offset. Binding a numeric zero offset
        # would incorrectly replace the graph's shared exact-zero node.
        geometry={'offset':owner.radius.maps[chart]['offset'].node}
        q.bindings[geometry['offset']]=L
        for name,powers in expected.items():
            parts=owner.scale(chart,geometry,powers);rp,ep,pp=powers
            got=sum(q.at(v) for v in parts.values())
            wanted=rp*(A+13/mu+L)+ep*(P+U-13/(2*mu)-13)+pp*P
            same(got,wanted,(chart,name));identities+=1
            if name=='Mztheta':
                assert parts['combined_inverse_mu_pulse'].node==owner.graph.zero.node
                cancelled+=1
            if name=='Mtheta':same(q.at(parts['combined_inverse_mu_pulse']),13/mu,'angular scale pulse coefficient')
    # Source equations are read from the actual general radial-history
    # implementation, including the sign and 1/2 in the energy history.
    E,V=s.symbols('E V',real=True)
    densities=generic.history_densities(E,V)
    assert densities==dict(m=V,h=E,k=E*V,e=V**2-E**2/2,p=E**2/2)
    y,z,v=s.symbols('y z v',real=True);J=s.Symbol('native_Jacobian',positive=True)
    R0,Ev0=s.symbols('R0 Ev0',positive=True);theta=s.Function('actual_theta')
    seedA,seedE,seedP=s.symbols('nonzero_angular_memory nonzero_energy_memory nonzero_pressure_memory')
    primal=dict(Mz=s.Integer(0),Mtheta_z=s.Integer(0),
        Mtheta=s.sqrt(2)*R0**s.Rational(3,2)*Ev0*(seedA+s.Integral(s.exp(s.Rational(3,2)*v)*theta(v,z),(v,0,y))),
        Mztheta=R0*Ev0**2*(seedE-s.Integral(s.exp(v)*theta(v,z)**2/2,(v,0,y))),
        Mp=seedP+Ev0**2*s.Integral(theta(v,z)**2/2,(v,0,y)))
    expected_derivatives=dict(Mz=0,Mtheta_z=0,
        Mtheta=s.sqrt(2)*(R0*s.exp(y))**s.Rational(3,2)*Ev0*theta(y,z)*J,
        Mztheta=-(R0*s.exp(y))*Ev0**2*theta(y,z)**2*J/2,
        Mp=Ev0**2*theta(y,z)**2*J/2)
    for name,value in primal.items():same(s.diff(value,y)*J,expected_derivatives[name],name+' cumulative derivative')
    # Pressure has an independent axis datum, whose native derivative is zero.
    P0=s.Function('independent_P0')(z)
    same(s.diff(P0+primal['Mp'],y)*J,expected_derivatives['Mp'],'pressure memory and derivative')
    return dict(exact_current_U0_function_bound=True,exact_Ev0_function_bound=True,
        original_native_flatten_amplitude_assignments=bindings,
        independent_absolute_scale_identity_count=identities,
        nine_energy_scale_inverse_mu_cancellations_are_exact=cancelled==9,
        general_signed_cumulative_density_function_identified=True,
        five_native_derivatives_from_independent_primal_integrals=True,
        retained_nonzero_seed_memory_and_independent_P0=True,passed=True)


def identical(a,b,label):
    assert a.order==b.order==5,label
    for n in range(6):assert a[n]._mpi_==b[n]._mpi_,(label,n)


def overlaps(a,b,label):
    count=0
    for n in range(6):
        al,ah=current.radius.post.selected.inlet.endpoints(a[n]);bl,bh=current.radius.post.selected.inlet.endpoints(b[n])
        assert max(al,bl)<=min(ah,bh),(label,n);count+=1
    return count


@source_precision
def run(owner=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    owner=owner if owner is not None else current.CurrentOriginalRpRawHistoryTransport(require_checked=False)
    assert not owner.acceptance_loaded and raw['candidate_current_raw_history_transport_constructed']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record and all(owner.assert_graph().values())
    # A new candidate owner gives deterministic expression IDs and fresh actual
    # calls; no stored truth flag or box endpoint supplies a field value.
    assert len(owner.prefix)==raw['accepted_radius_graph_prefix_length']
    assert owner.U0.node==raw['exact_current_inlet_U0_node']
    views={};rows=0;zero_rows=0;width_rows=0
    for name,chart,coordinate in current.VIEWS:
        view=owner.evaluate(chart,raw['fresh_Z'],coordinate);views[name]=view
        assert current.packed(current.view_report(view))==raw['actual_twenty_source_views'][name],name
        packet=view['current_source_packet'];h=view['source_normalized_histories'];theta=packet['theta_over_Ev0_Taylor']
        expected=dict(Mz=theta*h['Mz_over_R_Utheta'],
            Mtheta=theta*h['Mtheta_over_sqrt2_R_3half_Utheta']*owner.ctx.sqrt(2),
            Mtheta_z=theta*theta*h['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']*owner.ctx.sqrt(2),
            Mztheta=theta*theta*h['Mztheta_over_R_Utheta_squared'],
            Mp=h['Mp_over_Pstar_squared'],P0=h['P0_over_Pstar_squared'],
            pressure=h['P_over_Pstar_squared'],Utheta=theta)
        for key,value in view['raw_histories'].items():
            identical(value.in_exact_units(current.POWERS[key]),expected[key],(name,key));rows+=6
            assert len(value.log_scale_parts)==6
            assert all(part.graph is owner.graph for _,part in value.log_scale_parts)
            for bound in value.report()['directed_scaled_coefficient_width_bounds']:
                assert current.radius.post.selected.inlet.endpoints(bound)[0]>=0;width_rows+=1
            try:value.in_exact_units((0,0,0))
            except ValueError:pass
            else:raise AssertionError('A physical scale was silently discarded')
        for key in ('Mz','Mtheta_z'):
            for n in range(6):
                assert current.radius.post.selected.inlet.endpoints(view['raw_histories'][key].coefficients[n])==(0,0)
                assert current.radius.post.selected.inlet.endpoints(view['first_native_radial_derivatives'][key].coefficients[n])==(0,0)
                zero_rows+=2
        jac=view['directed_native_Jacobian_bound']
        expected_d=dict(Mtheta=theta*owner.ctx.sqrt(2)*jac,
            Mztheta=-theta*theta*jac/2,Mp=theta*theta*jac/2,pressure=theta*theta*jac/2)
        for key,value in expected_d.items():identical(view['first_native_radial_derivatives'][key].coefficients,value,(name,'derivative',key));rows+=6
        assert view['first_native_radial_derivatives']['Mp'].powers==tuple(map(Fraction,(0,2,0)))
        assert all(current.radius.post.selected.inlet.endpoints(v)==(0,0) for v in view['first_native_radial_derivatives']['P0'].coefficients.coefficients)
        assert view['raw_histories']['P0'].coefficients is h['P0_over_Pstar_squared']
        assert view['raw_histories']['pressure'].coefficients is h['P_over_Pstar_squared']
        assert packet['actual_terminal_zero_linear_and_radial_histories_inherited']
        assert not h['physical_raw_moment_evaluation_installed']
    assert current.packed(owner.graph.nodes[len(owner.prefix):])==raw['appended_exact_function_nodes']
    assert current.packed(owner.post.heat_constants(raw['fresh_Z']))==raw['retained_heat_terminal_constants']
    proof=exact_scale_and_cumulative_proof(owner)
    # Same physical radius/Ev0 scale on both sides of each seam. Interval
    # coefficient overlaps remain diagnostics, not terminal closure proofs.
    seam_rows=0
    for left,right in PAIRS:
        a,b=views[left],views[right]
        # The accepted absolute origin, inlet amplitude and true waiting
        # function are shared on both sides. Keep those unchanged dependencies
        # atomic rather than replaying the earlier source-parameter graph.
        qa=interpreter(owner,{owner.radius.logRp.node:s.Symbol('same_absolute_logRp'),
            owner.logU0.node:s.Symbol('same_logU0'),
            owner.radius.functions['waiting'].node:s.Symbol('same_true_waiting')})
        for key in ('Mz','Mtheta','Mtheta_z','Mztheta','Mp','P0','pressure'):
            av=a['raw_histories'][key];bv=b['raw_histories'][key]
            same(sum(qa.at(v) for _,v in av.log_scale_parts),sum(qa.at(v) for _,v in bv.log_scale_parts),(left,right,key,'scale'))
            seam_rows+=overlaps(av.coefficients,bv.coefficients,(left,right,key))
    rejected=[]
    for chart,coordinate in (('pulse_end',0),('flatten',-1),('waiting',2),('heat_exterior',2)):
        try:owner.evaluate(chart,'.439',coordinate)
        except ValueError:rejected.append([chart,coordinate])
        else:raise AssertionError('Invalid physical history domain accepted')
    try:current.CurrentOriginalRpRawHistoryTransport(owner.post,require_checked=False)
    except ValueError:pass
    else:raise AssertionError('A relative-only postpulse owner bypassed absolute geometry admission')
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        independent_exact_scale_and_cumulative_equation_proof=proof,
        actual_twenty_source_call_coefficient_identities=rows,
        inherited_exact_linear_zero_coefficient_checks=zero_rows,
        directed_scaled_error_width_checks=width_rows,
        eight_seam_scale_identities_and_coefficient_overlap_diagnostics=seam_rows,
        actual_P0_absolute_pressure_and_heat_constants_retained=True,
        wrong_units_relative_only_owner_and_invalid_domains_rejected=rejected,
        no_interval_caps_endpoints_or_midpoints_define_absolute_source_scales=True,
        actual_current_unique_repair_and_provider_graph_preserved=True,
        native_derivative_rows_are_directed_enclosures_not_exact_point_coefficients=True,
        scope='Factorized raw source histories and first native radial transport; unrestricted numerical absolute field, terminal closure and global mixed contracts remain open',
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name),
            Path(generic.__file__).name:current.sha(Path(generic.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_RAW_HISTORY_TRANSPORT 9 charts; raw histories and first native derivatives',flush=True)
    return result


if __name__=='__main__':run()

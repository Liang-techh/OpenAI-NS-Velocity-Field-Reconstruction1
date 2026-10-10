"""Focused current inverse admission, independent roots and log-input checks."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_physical_inverse as current
from lei_ren_part1_paper_compliant_current_original_Rp_physical_source_map_check import interpretation,PhysicalInterpreter
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_check import independent_root
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


class InverseInterpreter(PhysicalInterpreter):
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings:return self.bindings[i]
        n=self.g.nodes[i]
        if n['operation']=='exact_physical_atan2':
            assert n['branch']=='(-pi,pi]; negative x with exact y=0 uses +pi'
            assert n['positive_radius_certificate']=='x^2+y^2>0'
            return s.atan2(self.at(n['y']),self.at(n['x']))
        return super().at(i)


def exact_function_checks(owner,views):
    base,delta=interpretation(owner.before)
    identities=nodes=0
    for name,view in views.items():
        # Different physical inputs can share a root node. Keep each symbolic
        # root binding and its memoized descendants local to this observation.
        reader=InverseInterpreter(type('Field',(),dict(graph=owner.graph,parameters=owner.radius.parameters,
            contracts=owner.radius.frame.bridge.leading.contracts))(),False,dict(base.bindings))
        refs=view['coordinate_functions'];qref=refs['log_lambda'];inputs=view['exact_input_record']
        sign=inputs.get('sign_z') if view['input_kind']=='exact_log_cylindrical' else (1 if Fraction(inputs['z'])>0 else -1 if Fraction(inputs['z'])<0 else 0)
        if sign:
            node=owner.graph.nodes[qref.node]
            assert node['operation']=='exact_original_physical_log_lambda_inverse'
            assert node['delta']==owner.delta_function.node and node['viscosity']==owner.graph.one.node
            assert node['physical_z_sign']==sign and node['log_tau']==refs['log_tau'].node
            assert node['definition']=='unique finite q: 2q=logaddexp(log_tau,2log_abs_z+2delta*q)'
            assert node['positive_root_certificate']=='tau=exp(log_tau)>0; 0<delta<1; F_prime>=2(1-delta)>0'
            q=s.Symbol('unique_actual_log_lambda_'+name,real=True);reader.bindings[qref.node]=q;nodes+=1
            lz=reader.at(node['log_abs_z'])
            assert s.simplify(reader.at(refs['Z'])-sign*s.exp(lz-(1-delta)*q))==0
            assert s.simplify(reader.at(refs['log_abs_Z'])-lz+(1-delta)*q)==0
            identities+=2
        else:
            q=reader.at(refs['log_tau'])/2
            assert s.simplify(reader.at(qref)-q)==0 and reader.at(refs['Z'])==0
            assert refs['log_abs_Z'] is None;identities+=2
        lr=reader.at(refs['log_r']);lt=reader.at(refs['log_tau'])
        assert s.simplify(reader.at(refs['logR'])-(2*lr-s.log(2)-2*q))==0
        assert s.simplify(reader.at(refs['log_one_minus_Z_squared'])-(lt-2*q))==0
        assert s.simplify(reader.at(refs['lambda_value'])-s.exp(q))==0
        if view['input_kind']=='exact_Cartesian':
            x,y,z,t=(s.Rational(inputs[key]) for key in ('x','y','z','t'))
            assert s.simplify(lr-s.log(x*x+y*y)/2)==0 and s.simplify(lt-s.log(1-t))==0
            assert s.simplify(reader.at(refs['theta'])-s.atan2(y,x))==0
            assert all(reader.at(refs[key])==value for key,value in (('x',x),('y',y),('z',z),('requested_t',t)))
            identities+=7
        identities+=4
        assert all(ref is None or ref.graph is owner.graph for ref in refs.values())
        assert not view['numerical_absolute_radius_or_lambda_materialized']
        assert view['numerical_bounds_are_not_defining_coordinate_functions']
    return dict(passed=True,exact_inverse_function_nodes=nodes,exact_coordinate_function_identities=identities,
        exact_current_delta_and_nu1_source_retained=True,zero_axial_inverse_is_exact=True)


def number(q):
    q=Fraction(q);return mp.mpf(q.numerator)/q.denominator


def independent_current_points(owner,views):
    result={};ends=current.ends
    for name,args in current.CASES.items():
        x,y,z,t=map(number,args);tau=1-t;lt=mp.log(tau);lr=mp.log(x*x+y*y)/2
        view=views[name];mapping=view['directed_inverse_mapping'];lo,hi=ends(mapping['actual_log_lambda'])
        rl,rh=ends(view['source_logR_enclosure']);zl,zh=ends(mapping['Z']);al,ah=ends(view['theta_enclosure'])
        count=0
        for delta in ends(owner.delta):
            ql,qh=independent_root(z,lt,1,delta)
            assert lo<=ql<=qh<=hi,(name,'independent current delta corner root')
            for q in (ql,qh):
                assert rl<=2*lr-mp.log(2)-2*q<=rh
                Z=mp.sign(z)*mp.sqrt(max(mp.mpf(0),1-mp.exp(lt-2*q)))
                assert zl<=Z<=zh,(name,'independent Z')
            count+=1
        assert al<=mp.atan2(y,x)<=ah,(name,'directed angle')
        assert mapping['source_delta']._mpi_==owner.delta._mpi_
        assert mapping['physical_viscosity']._mpi_==owner.ctx.mpf(1)._mpi_
        assert mapping['log_root_enclosure_width'].b<mp.mpf('1e-50')
        result[name]=dict(independent_direct_exponential_current_delta_corner_roots=count,
            logR_Z_and_quadrant_angle_inside=True,
            log_root_width=mapping['log_root_enclosure_width'],solver_status=mapping['solver_status'])
    return result


def log_input_checks(owner,views):
    c=owner.ctx;ends=current.ends;fixtures={};count=0
    for label,lz,lt,d,sign in (('positive','-2','-4','.03',1),
            ('negative','2','-7','.8',-1),('tiny_axial','-200','-2','.4',1)):
        out=current.log_coordinate_map(c,lz,lt,d,sign)
        z=sign*mp.exp(mp.mpf(lz));ql,qh=independent_root(z,lt,1,d)
        lo,hi=ends(out['actual_log_lambda']);assert lo<=ql<=qh<=hi
        old=current.original.implicit_log_coordinate_map(c,sign*c.exp(c.mpf(lz)),lt,1,d,'1e-60',512)
        al,ah=ends(old['actual_log_lambda']);assert max(lo,al)<=min(hi,ah)
        fixtures[label]=dict(independent_root_inside=True,original_direct_z_program_overlap=True,
            fixture_only_not_current_profile_coefficients=True);count+=1
    for name in current.LOG_CASES:
        view=views[name];mapping=view['directed_inverse_mapping']
        assert all(mp.isfinite(v) for v in ends(mapping['actual_log_lambda'])+ends(view['source_logR_enclosure']))
        assert mapping['source_delta']._mpi_==owner.delta._mpi_
        assert mapping['exp_truncation_is_enclosure_only_not_source']
        assert not view['numerical_absolute_radius_or_lambda_materialized']
    assert not views['near_time_nonzero_z']['strict_numerical_Z_interior_resolved']
    # Loss of numeric precision at Z=+/-1 is kept visible; the exact
    # finite physical inverse still has strict |Z|<1 by its theorem.
    return dict(passed=True,independent_log_input_fixtures=fixtures,
        original_nonzero_solver_program_preserved_fixtures=count,
        actual_current_astronomical_log_input_cases=4,
        infinity_limit_touching_enclosure_not_admitted_as_point=True,
        no_astronomical_physical_z_lambda_or_radius_exponent_materialized=True)


@source_precision
def run(before=None,observed_owner=None,observed_views=None):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpPhysicalInverse(before,require_checked=False)
    assert not any(candidate[key] for key in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record and candidate['actual_source_graph']==owner.assert_graph()
    assert candidate['original_log_axial_program_binding']==owner.program_binding
    assert candidate['current_forward_unit_viscosity_binding']==owner.forward_unit_binding
    assert candidate['original_implicit_inverse_theorem']==owner.theorem
    assert candidate['independent_unique_axial_inverse_theorem']==owner.independent_theorem
    assert candidate['accepted_original_root_fixture_evidence']==owner.original_root_fixture_evidence
    if observed_owner is None:
        views={name:owner.cartesian(*args) for name,args in current.CASES.items()}
        views.update({name:owner.log_cylindrical(*args) for name,args in current.LOG_CASES.items()})
    else:
        assert type(observed_owner) is type(owner) and observed_owner.before is before
        assert observed_owner.graph is owner.graph and observed_owner.ctx is owner.ctx
        assert observed_owner.hashes==owner.hashes and all(observed_owner.assert_graph().values())
        views=observed_views
    assert set(views)==set(current.CASES)|set(current.LOG_CASES)
    graph=candidate['exact_current_inverse_expression_graph'];assert graph==owner.graph.nodes[:len(graph)]
    for name,view in views.items():
        assert current.physical.mixed.pulse.raw.packed(current.report(view))==candidate['actual_current_inverse_views'][name]
        assert not any(view[key] for key in current.GATES+current.OPEN)
    exact=exact_function_checks(owner,views)
    points=independent_current_points(owner,views);logs=log_input_checks(owner,views)
    rejected=[]
    for name,fn in (
        ('axis_requires_current_core',lambda:owner.cartesian(0,0,'.3','.9')),
        ('terminal_time',lambda:owner.cartesian(1,0,0,1)),
        ('after_terminal',lambda:owner.cartesian(1,0,0,2)),
        ('nonfinite_coordinate',lambda:owner.cartesian(float('inf'),0,0,0)),
        ('interval_is_not_physical_point_input',lambda:owner.cartesian(owner.ctx.mpf(1),0,0,0)),
        ('zero_z_has_no_finite_log',lambda:owner.log_cylindrical(0,0,0,-10)),
        ('nonzero_z_requires_log',lambda:owner.log_cylindrical(0,None,1,-10)),
        ('invalid_axial_sign',lambda:owner.log_cylindrical(0,0,2,-10)),
        ('zero_tolerance',lambda:owner.log_cylindrical(0,0,1,-10,relative_tolerance='0'))):
        try:fn()
        except (ValueError,TypeError):rejected.append(name)
        else:raise AssertionError('Wrong inverse input accepted: '+name)
    bad=copy.copy(owner);bad.delta_function=owner.graph.constant(0)
    try:bad.assert_graph()
    except ValueError:rejected.append('delta_enclosure_or_zero_as_defining_function')
    else:raise AssertionError('Changed exact inverse delta accepted')
    foreign=current.physical.mixed.pulse.radius.FunctionRef(copy.deepcopy(owner.graph),0)
    try:owner.inverse_functions(foreign,None,owner.graph.zero,owner.graph.zero,0)
    except TypeError:rejected.append('foreign_physical_input_graph')
    else:raise AssertionError('Foreign inverse graph accepted')
    sample=views['ordinary']
    try:owner.assemble(sample['input_kind'],sample['exact_input_record'],sample['coordinate_functions'],
        sample['directed_inverse_mapping'],sample['physical_log_radius_enclosure'],sample['theta_enclosure'])
    except ValueError:rejected.append('unwitnessed_or_foreign_solver_output')
    else:raise AssertionError('Unwitnessed coordinate inverse accepted')
    assert current.exact_scalar(0.1)==Fraction.from_float(0.1)!=Fraction('0.1')
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        original_log_axial_program_binding=owner.program_binding,
        current_forward_unit_viscosity_binding=owner.forward_unit_binding,
        exact_actual_inverse_function_transfer=exact,independent_current_Cartesian_point_checks=points,
        explicit_log_axial_inverse_checks=logs,
        binary_float_inputs_retain_exact_binary_value=True,
        original_implicit_inverse_theorem=owner.theorem,independent_unique_axial_inverse_theorem=owner.independent_theorem,
        original_N_definition=owner.original_N_definition,original_delta_enclosure_not_a_point_value=True,
        no_old_locator_or_tensor_registry_owner_constructed=True,
        checked_result_requires_unchanged_live_current_solver_witness=True,
        native_chart_and_point_source_admission_left_open=True,rejected_sources_and_domains=rejected,
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            current.PREFIX+'current_original_Rp_physical_source_map_check.py':current.sha(current.PREFIX+'current_original_Rp_physical_source_map_check.py')},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.physical.mixed.pulse.raw.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_PHYSICAL_INVERSE exact Cartesian/log inputs and directed roots',flush=True)
    return result


if __name__=='__main__':run()

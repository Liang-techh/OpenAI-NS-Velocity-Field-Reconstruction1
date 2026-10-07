"""Independent signed cone, flat-cutoff and current collar source checks."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_inner_exit_strict_collar as source
from lei_ren_part1_paper_compliant_current_O2_reference_slope_relaxed_cone import relaxed_cone_margins


def independent_signed_cone_fixtures(c):
    """Paper D/J directly versus the inherited perturbative dot/cross bound.

    Moderate independent data exercise both axial signs, off-direction
    errors and vanishing stress amplitudes. They do not define source K/hb.
    """
    count=0
    with mp.workdps(150):
        for kv in ('2','1000000','1e10'):
            K=mp.mpf(kv)
            for dv,ev in (('3','4'),('3','-4'),('3','0'),('.4','1'),('.4','-1')):
                d,e=mp.mpf(dv),mp.mpf(ev);H=(d*d+e*e)/d
                if H>K**10 or d<1/(2*K):raise ArithmeticError('Independent fixture outside inherited source hypotheses')
                for ov in ('.01','.1','.25','1e-40'):
                    omega=mp.mpf(ov);chi=1-omega;a=chi*d;b=-chi*e;kappa=a+b*b/a
                    if kappa<=2:raise ArithmeticError('Independent fixture outside strict branch')
                    for xv,yv in (('-.7','-.7'),('-1','0'),('1','0'),('0','1'),('0','-1'),('.7','.7'),('0','0')):
                        et=mp.mpf('.99')*omega/(40*K**6)*mp.mpf(xv)
                        ez=mp.mpf('.99')*omega/(40*K**6)*mp.mpf(yv)
                        theta=omega*d+et;axial=omega*e+ez
                        D=theta-b*axial/a;J=axial+b*theta/a
                        Q=2*D*D-(kappa-2)*J*J
                        dot=(theta*d+axial*e)/d;cross=(d*axial-e*theta)/d
                        if abs((D-dot)/(omega*H))>mp.mpf('1e-110') or abs((J-cross)/(omega*H))>mp.mpf('1e-110'):
                            raise ArithmeticError('Independent paper signed dot/cross units differ')
                        if D/(omega*H)<=mp.mpf('.95') or Q/(omega*H)**2<=mp.mpf('1.8'):
                            raise ArithmeticError('Independent full cone direction/quadratic bound fails')
                        if mp.sqrt(theta**2+axial**2)/(omega*mp.sqrt(d*d+e*e))<=mp.mpf('.95'):
                            raise ArithmeticError('Nonzero stress relative lower lost')
                        values=[mp.nstr(v,140) for v in (a,b,theta,axial)]
                        generic=relaxed_cone_margins(c,*values)
                        if not generic['admitted'] or not generic['strict_source_shear_for_entire_box']:
                            raise ArithmeticError('Independent paper cone adapter rejects strict source fixture')
                        count+=1
    return count


def independent_flat_cutoff_fixtures():
    count=0
    with mp.workdps(150):
        for L in (mp.log(2),mp.mpf(10),mp.mpf(1000)):
            gamma=mp.mpf('.01');X=10*L+mp.log(10/gamma);sc=1/(4*mp.sqrt(X))
            for f in (mp.mpf('.01'),mp.mpf('.5'),mp.mpf(1)):
                x=f*sc;a=mp.exp(-1/x**2);b=mp.exp(-1/(1-x)**2)
                sigma=a/(a+b);odds=1/(1-x)**2-1/x**2
                if abs(sigma/(mp.exp(odds)/(1+mp.exp(odds)))-1)>mp.mpf('1e-110'):
                    raise ArithmeticError('Independent original flat sigma definitions differ')
                if mp.log(sigma)>4-1/sc**2 or mp.log(sigma)+10*L>=mp.log(gamma/10):
                    raise ArithmeticError('Independent whole-subcollar log bound fails')
                count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Current strict collar dependency changed: '+name)
    field=source.CurrentInnerExitStrictCollar(require_checked=False);c=field.ctx
    for key,value in (('exact_current_exit_source_attachment',field.theorem),
        ('cross_graph_amplitude_packet_pressure_width_attachment',field.attachment),
        ('explicit_current_inner_exit_strict_collar',field.proof)):
        if source.encode(value)!=data[key]:raise ValueError('Current exact collar source proof changed: '+key)
    digest=hashlib.sha256(json.dumps(field.attachment,sort_keys=True).encode()).hexdigest()
    if digest!=field.attachment_digest or digest!=data['cross_graph_attachment_sha256']:
        raise ValueError('Cross-graph amplitude/packet/pressure/width digest differs')
    for name,value in field.proof['positive_margins'].items():
        if source.endpoints(value)[0]<=0:raise ArithmeticError('Whole source strict collar inequality missing: '+name)
    if source.endpoints(field.proof['correlated_K10_relative_error_squared_upper'])!=source.endpoints(c.mpf(1)/400):
        raise ArithmeticError('K^10/K^10 error correlation lost')
    specs=(('closed_core_and_collar',(0,1),(-1,1)),('exact_zero_inlet',0,(-1,1)),
        ('strict_left_support_collar',('.5','1'),(-1,1)),('strict_positive_midplane',('.25','.75'),0))
    for name,q,z in specs:
        query=field.query(q,z)
        if source.encode(query)!=data['examples'][name]:raise ValueError('Current strict collar query changed: '+name)
        if query['strict_nonzero_stress_cone_certified_for_entire_query_box']!=name.startswith('strict'):
            raise ArithmeticError('Exact zero stress admitted as nonzero strict collar')
    invalid=((-1,0),(2,0),(('.99','1.01'),0),('.5',2),('.5',('-1.01','0')),(mp.inf,0),('.5',mp.inf))
    for q,z in invalid:
        try:field.query(q,z)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Outside current strict collar domain admitted')
    if not data[source.GATE] or any(data[key] is not False for key in source.OPEN):
        raise ValueError('Scoped current collar promoted unfinished global work')
    cones=independent_signed_cone_fixtures(c);cutoff=independent_flat_cutoff_fixtures()
    result=dict(all_passed=True,actual_five_defect_family_sha256=field.family,
        implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum,
        cross_graph_attachment_sha256=field.attachment_digest,
        exact_source_AST_bindings=len(field.theorem['AST_bindings']),
        exact_original_unit_signed_cone_and_cutoff_identities=len(field.theorem['exact_source_identities']),
        strict_whole_source_subcollar_inequalities=len(field.proof['positive_margins']),
        independent_full_signed_paper_cone_fixtures=cones,independent_flat_cutoff_source_fixtures=cutoff,
        current_source_queries_checked=len(specs),invalid_queries_rejected=len(invalid),
        source_stress_omega_factor_retained_not_materialized=True,
        six_actual_inlet_atoms_and_original_pressure_and_width_attached_by_source_functions=True,
        exact_zero_core_inlet_not_promoted_nonzero=True,
        **{source.GATE:True},**{key:False for key in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(source.encode(result),indent=2)+'\n').encode())
    print('Current source inner strict collar PASS:',cones,'signed cone cases;',cutoff,'independent cutoff cases',flush=True)
    return result


if __name__=='__main__':run()

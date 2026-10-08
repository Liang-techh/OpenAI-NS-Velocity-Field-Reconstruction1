"""Independent original axial derivatives/widths and saved actual source replay."""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_axial_endpoint_q_source as source
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source_check as compare
from lei_ren_part1_paper_compliant_current_transition_complete_prefix_check import exact_replay_equal

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv


def independent_original_axial_checks(family):
    c=MPIntervalContext();c.dps=100;p=mp.mp.clone();p.dps=300
    comparisons=0;widths=0;crossings=0
    M=p.mpf(4);logP=p.exp(M)+11
    for H in (2048,8192):
        eta_log=-p.mpf(H)-2*logP;eta=p.exp(eta_log)
        coordinates=source.native.HalfPstarCoordinates(c,2*(c.exp(4)+11),family)
        coord=source.AxialCoordinates(coordinates,-c.mpf(H)-coordinates.logP_squared,4)
        def sigprime(phase):
            if phase<=0 or phase>=1:return p.mpf(0)
            r=min(phase,1-phase);L=1/(1-r)**2-1/r**2
            tail=p.exp(L);return tail/(1+tail)**2*(2/(1-r)**3+2/r**3)
        def original_b(y,z):
            phase=p.log(y)/M
            if y==1 or y==p.exp(M):return p.mpf(0)
            return -8*p.exp(p.mpf(1)/5)*z*(1+z*z)/M*sigprime(phase)*p.exp(-logP+(y-1)/2-p.log(y))
        def original_rho(y,z):return original_b(y,z)**2/(2*eta)
        def scalar_q_rho_jets(rho):
            if rho>=1:return [p.mpf(0)]*4
            root=p.sqrt(eta*(2-rho)/4)
            chi=p.mpf(1);chi1=chi2=chi3=p.mpf(0)
            if rho:
                C=1/rho**2-1/(1-rho)**2
                H1=2/rho**3+2/(1-rho)**3;H2=6/rho**4-6/(1-rho)**4;H3=24/rho**5+24/(1-rho)**5
                if C>0 and p.log(C)>4*p.dps*p.log(10):
                    assert C-p.log(1+H1+abs(H2)+H3+H1**3+abs(H1*H2))>(p.dps+100)*p.log(10)
                else:
                    tail=p.exp(-C);chi=1/(1+tail);cc=tail/(1+tail);prod=chi*cc
                    chi1=-prod*H1;chi2=prod*(H2+(1-2*chi)*H1**2)
                    chi3=-prod*(H3+3*(1-2*chi)*H1*H2+(1-6*prod)*H1**3)
            ell1=-1/(2*(2-rho));ell2=-1/(2*(2-rho)**2);ell3=-1/(2-rho)**3
            rr1=root*ell1;rr2=root*(ell2+ell1**2);rr3=root*(ell3+3*ell1*ell2+ell1**3)
            return [chi*root,chi1*root+chi*rr1,chi2*root+2*chi1*rr1+chi*rr2,
                chi3*root+3*chi2*rr1+3*chi1*rr2+chi*rr3]
        for Z in (('.36','.38'),('-.38','-.36')):
            for side in ('left','right'):
                print('Independent original axial source:',H,Z,side,flush=True)
                A=p.mpf(0) if side=='left' else p.exp(M)-1-2*M
                Q=p.mpf(H)+3*p.log(H)+A
                def radius(k):return p.sqrt(2/(Q+k))
                def phi(r):return r if side=='left' else 1-r
                zm=(p.mpf(Z[0])+p.mpf(Z[1]))/2
                # Independent root of the original source selects only a
                # finite test cell crossing rho=1, never production data.
                center=p.findroot(lambda k:p.log(original_rho(p.exp(M*phi(radius(k))),zm)),(0,-5))
                kupper=Fraction(int(p.ceil(10*center)),10);klower=kupper-Fraction(1,10)
                geometries=[coord.geometry(side,'xi','0','3/4'),coord.geometry(side,'mixed','3/4','8'),
                    coord.geometry(side,'k','8','6'),coord.geometry(side,'k',kupper,klower),
                    coord.geometry(side,'k',klower-Fraction(1,2),klower-Fraction(3,5))]
                for geometry in geometries:
                    norm=source.normalized_excess(coord,Z,geometry)
                    q=source.original_axial_q(coord,norm['rows'])
                    if geometry['kind']=='xi':
                        ra=p.mpf(0);rb=p.sqrt(2/Q)*p.mpf('.75')
                    elif geometry['kind']=='mixed':ra=3*p.sqrt(2/Q)/4;rb=radius(8)
                    else:
                        ka=p.mpf(str(source.ep(geometry['left'])[0]));kb=p.mpf(str(source.ep(geometry['right'])[0]))
                        ra=radius(ka);rb=radius(kb)
                    compare.enclosed(p,c,geometry['width'],rb-ra);widths+=1
                    pa,pb=(ra,rb) if side=='left' else (1-rb,1-ra)
                    compare.enclosed(p,c,geometry['physical_width'],p.exp(M*pb)-p.exp(M*pa));widths+=1
                    for r in (ra,(ra+rb)/2,rb):
                        y=p.exp(M*phi(r))
                        for z in (zm,):
                            brefs={o:p.diff(original_b,(y,z),o) if r else p.mpf(0) for o in source.ORDERS}
                            b=brefs[source.ZERO];by=brefs[(1,0)];byy=brefs[(2,0)]
                            bz=brefs[(0,1)];byz=brefs[(1,1)];byyz=brefs[(2,1)]
                            # Independent derivatives of the original b
                            # supply the exact rho=b^2/(2eta) product rule.
                            # No second numerical differentiation is needed.
                            refs={source.ZERO:b*b/(2*eta),(1,0):b*by/eta,(2,0):(by*by+b*byy)/eta,
                                (0,1):b*bz/eta,(1,1):(bz*by+b*byz)/eta,
                                (2,1):(2*by*byz+bz*byy+b*byyz)/eta}
                            F=scalar_q_rho_jets(refs[source.ZERO]);ry=refs[(1,0)];ryy=refs[(2,0)];rz=refs[(0,1)]
                            ryz=refs[(1,1)];ryyz=refs[(2,1)]
                            qrefs={source.ZERO:F[0],(1,0):F[1]*ry,(2,0):F[2]*ry**2+F[1]*ryy,
                                (0,1):F[1]*rz,(1,1):F[2]*ry*rz+F[1]*ryz,
                                (2,1):F[3]*ry**2*rz+F[2]*(ryy*rz+2*ry*ryz)+F[1]*ryyz}
                            for o in source.ORDERS:
                                compare.enclosed(p,c,norm['b'][o],brefs[o]);comparisons+=1
                                compare.enclosed(p,c,norm['rows'][o],refs[o]);comparisons+=1
                                compare.enclosed(p,c,q['rows'][o],qrefs[o]);comparisons+=1
                    if q['record']['branch']=='smooth_cutoff_seam':crossings+=1
                # Direct original scalar rho supports the independent
                # middle lower theorem at two boundary and middle samples.
                proof=source.middle_flat_proof(coord,Z)
                star=p.sqrt(2/(H-p.log(H)-200))
                lower=p.mpf(str(ep(proof['original_log_rho_uniform_lower'])[0]))
                for phase in (star,p.mpf('.25'),p.mpf('.5'),p.mpf('.75'),1-star):
                    for z in (p.mpf(Z[0]),p.mpf(Z[1])):
                        assert p.log(original_rho(p.exp(M*phase),z))>=lower;comparisons+=1
    assert crossings
    return dict(passed=True,independent_original_b_rho_q_ordinary_y2_Z1_and_middle_comparisons=comparisons,
        independent_true_phase_and_ordinary_y_width_comparisons=widths,
        independent_original_source_cutoff_crossing_queries=crossings,
        original_scalar_b_differentiated_in_ordinary_y_with_independent_precision=True,
        independent_q_reference_analytic_with_proved_subprecision_cutoff_guard=True,
        finite_fixture_parameters_not_substituted_for_actual_original_source=True)


def run():
    began=time.monotonic();m=json.loads((HERE/source.NAME).read_bytes())
    assert m[source.GATE] and m['original_whole_axial_phase_interval']==['0','1']
    assert not m['full_Z_axis_actual_axial_density_integrals_or_Rc_targets_controls_admitted']
    for key in source.common.current.FLAGS:assert m[key] is False
    for name,digest in m['input_hashes'].items():assert sha(name)==digest,name
    c=MPIntervalContext();c.dps=240;count=0;qrows=0
    with mp.workdps(300):
        independent=independent_original_axial_checks(m['source_family'])
        for archive in m['original_axial_q_source_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            data=json.loads(raw);par=data['original_parameter_sources']
            coordinates=source.native.HalfPstarCoordinates(c,iv(c,data['saved_original_coordinate_theorem']['common_log_bases'][1]),m['source_family'])
            coord=source.AxialCoordinates(coordinates,iv(c,par['original_eta_log']),iv(c,par['original_Md']))
            replay=source.execute_tile(coord,data['exact_Z_range'])
            for key,value in replay.items():exact_replay_equal(json.loads(json.dumps(source.encode(value))),data[key],key)
            assert source.encode(source.source_identity())==data['original_axial_source_identity']
            actual=data['original_endpoint_source_queries']
            for side in ('left','right'):
                selected=[r for r in actual if r['endpoint']==side]
                assert [r['label'] for r in selected]==[r[0] for r in source.PLAN]
                for row,plan in zip(selected,source.PLAN):
                    geom=row['original_true_geometry'];assert geom['exact_declared_typed_endpoints']==list(plan[2:])
                    assert ep(iv(c,geom['positive_true_phase_width']['coefficient_interval']))[0]>0
                    assert ep(iv(c,geom['positive_true_ordinary_y_width']['coefficient_interval']))[0]>0
                    assert geom['microscopic_original_y_increment_retained'] and geom['physical_width_Jacobian_once']
                    assert row['original_q_source']['status']=='enclosed'
                    sign=-1 if Fraction(data['exact_Z_range'][0])>0 else 1
                    coef=ep(iv(c,row['original_b_ordinary_y2_Z1_rows']['y0_Z0']['coefficient_interval']))
                    assert coef[1]<=0 if sign<0 else coef[0]>=0
                    count+=1;qrows+=len(source.ORDERS)
                assert selected[-1]['original_q_source']['branch']=='flat'
                assert selected[0]['original_q_source']['branch']=='active'
                assert selected[-2]['original_q_source']['branch']=='smooth_cutoff_seam'
            assert data['original_middle_flat_source_theorem']['entire_original_middle_Delta_strictly_exceeds_eta']
            assert data['original_whole_phase_partition']['endpoint_q_nonzero_and_interior_flat_distinguished']
    hashes=dict(m['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=m['source_family'],
        independent_original_axial_source_checks=independent,actual_original_axial_endpoint_queries_replayed=count,
        actual_original_q_ordinary_y2_Z1_rows_replayed=qrows,
        exact_original_full_phase_partition_two_endpoints_and_flat_middle_checked=True,
        actual_b_sign_and_nonzero_q_endpoints_and_Z_jets_checked=True,
        source_true_width_and_ordinary_derivative_conversions_once_checked=True,
        accepted_upstream_or_transition_producers_not_rerun=True,
        full_Z_density_integrals_targets_controls_recursion_and_NS_not_claimed=True,
        **dict.fromkeys(source.common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Whole original strict-sign axial b/q endpoint sources and flat middle PASS; density/transport open',flush=True)
    return result


if __name__=='__main__':run()

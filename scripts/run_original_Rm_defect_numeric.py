"""Check actual source Rm signed defects against independent raw primitives."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_Rm_defect_numeric import OriginalRmDefectNumeric,centered_frame
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_Rm_defect_numeric.json'


def independent_primitive_checks():
    """Direct radial integration of independently modified velocity densities.

    R=s^10 removes endpoint fractional powers. Five own moment primitives
    integrate actual E,V, rather than centered defect or sector recurrences.
    """
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    c=mp.mp.clone();c.dps=60;checks=0
    for Z in ('-.37','0','.23'):
        z=J.variable(c,c.mpf(Z),2);am=c.exp(c.mpf('-.6'))/(1+z*z)
        P=c.mpf('3.1');j=c.mpf('.03');eps=c.mpf('.2');psi=1+z*z
        a=(1+z)*c.mpf('.017');b=(1-z+z*z)*c.mpf('-.011')
        expa=a.exp();zero=J.constant(c,0,2)
        def density(s,modified):
            r=s**10;alpha=1-r*r
            V=z*4+(psi*eps+j)*alpha
            E=am*(P*s)
            if modified:V+=b*P;E=E*expa
            return dict(m=V*(10*s**9/P),h=E*(10*s**14/P),
                k=E*V*(10*s**14/(P*P)),e=(V*V-E*E/2)*(10*s**9/(P*P)),
                p=(E*E)*(5/(P*P*s)) if s else zero)
        def integrate(modified):
            return {n:J(c,[c.quad(lambda s:density(s,modified)[n][k],[0,1])
                            for k in range(3)]) for n in ('m','h','k','e','p')}
        raw=integrate(False);changed=integrate(True);inc={n:changed[n]-raw[n] for n in raw}
        # Independent exact integrals of alpha=1-R² and its square.
        mean=c.mpf(2)/3;square=1-c.mpf(2)/3+c.mpf(1)/5
        mixed=c.mpf(5)/8-c.mpf(5)/18
        six=dict(Mpsi=psi*mean,Kpsi_over_Phi=psi*(2*mixed),Apsi=psi*psi*square)
        base=dict(gamma1=J.constant(c,mean,2),gamma2=J.constant(c,square,2),
                  Kalpha_over_Phi=J.constant(c,2*mixed,2),Dpsi=psi*square)
        angular=dict.fromkeys(('H','B','C'),zero)
        result=centered_frame(c,J,Z,j,eps,1/P,am,six,base,angular,
            increments=inc,inverse_exp=((-a).exp()-1,(-a*2).exp()-1))
        modified_am=am*expa
        direct=[changed['m']*P-z*4,
            (changed['k']*P-z*changed['h']*4)/modified_am,
            changed['h']/modified_am-c.mpf(5)/8,
            (changed['e']-z*changed['m']*(8/P)+z*z*(16/(P*P)))/(modified_am*modified_am)+c.mpf(5)/12,
            changed['p']/(modified_am*modified_am)-c.mpf(5)/2]
        for aa,bb in zip(result,direct):
            for k in range(3):
                assert abs(aa[k]-bb[k])<c.mpf('1e-45'),(Z,k,aa[k],bb[k]);checks+=1
        # Constant axial velocity with exact power-law swirl closes all five
        # background identities, not just its two linear moment rows.
        null=centered_frame(c,J,Z,0,0,1/P,am,
            dict(Mpsi=zero,Kpsi_over_Phi=zero,Apsi=zero),
            dict(gamma1=zero,gamma2=zero,Kalpha_over_Phi=zero,Dpsi=zero),angular)
        assert all(v==0 for row in null for v in row.coefficients);checks+=15
    return checks


def run():
    start=time.perf_counter();independent=independent_primitive_checks()
    print('Independent modified radial primitive and five-identity checks pass',flush=True)
    f=OriginalRmDefectNumeric();samples=[];moments=memory=scale=datum=0
    old=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_reference_velocity.json')).read_text(encoding='utf-8'))
    for Z in ('-.3','0','.3','.17'):
        p=f.evaluate(Z);samples.append(p)
        assert not p['limits']['finite_N_active_patch_uniform_inverse']
        assert p['source_contract']['analytic_centering_before_materialization']
        assert p['source_contract']['unknown_signed_errors_not_field_values']
        assert p['source_contract']['P0_added_once_and_excluded_from_defects']
        assert len(p['independent_P0_over_Pstar_squared'])==6;datum+=1
        if Z!='.17':
            previous=next(x for x in old['samples'] if x['Z']==Z and x['chart']=='postrestore' and x['coordinate']['value']=='1')
            assert p['selected_integer']==previous['selected_integer'];memory+=1
            assert p['original_five_correction_history_bounds']==previous['original_five_correction_history_bounds'];memory+=1
            assert p['positive_order_n1_own_history']==previous['positive_order_n1_own_history'];memory+=1
        for group in ('signed_unmodified_defects','actual_finite_N_defect_enclosures','actual_finite_N_canonical_frame_defect_enclosures'):
            for row in p[group]:
                assert len(row)==3
                for coefficient in row:
                    lo,hi=_bounds(f.ctx,f.reference.previous._restore(coefficient))
                    assert f.ctx.isfinite(lo) and f.ctx.isfinite(hi);moments+=1
        # Original finite pressure datum is shared and added exactly once.
        for k,row in enumerate(p['own_pressure_over_Pstar_squared']):
            pp=f.reference.previous._restore(p['independent_P0_over_Pstar_squared'][k])
            inc=f.reference.previous._restore(p['generic_own_moment_enclosures']['p'][k])
            actual=f.reference.previous._restore(row)
            a,b=_bounds(f.ctx,pp+inc);l,h=_bounds(f.ctx,actual)
            assert max(a,l)<=min(b,h);datum+=1
        # Stable centering must keep angular uncertainty far below the
        # actual original angular control unit, exp(-2logP)*1e10*t².
        repair=json.loads((DIRECTORY/(PREFIX+'compliant_five_moment_repair.json')).read_text(encoding='utf-8'))
        t=f.iv.mpf(mp.make_mpf(tuple(repair['axial_box_scale']['lower_exact_mpf_tuple'])))
        unit=f.iv.exp(-2*f.logP)*(10**10*t*t)
        for n in (2,4):
            for row in p['actual_finite_N_defect_enclosures'][n]:
                bound=max(abs(x) for x in _bounds(f.ctx,f.reference.previous._restore(row)))
                assert bound<_bounds(f.ctx,unit)[0];scale+=1
        print('Actual original Rm signed frame: Z='+Z,flush=True)
    # Separate source/scalar normalization check using a core packet, not
    # another export of the new datum routine.
    core=f.reference.leading.core;packet=core.coefficients('.17',degree=4,depth=1)
    scalar=core.evaluate(packet,'0')['profiles']['P0_over_Pstar_squared']
    row=samples[-1]['independent_P0_over_Pstar_squared'][0]
    overlap(row,scalar);datum+=1
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_radial_primitive_and_exact_five_identity_coefficients=independent,
            original_Rm_signed_frame_packets=len(samples),finite_centered_moment_coefficients=moments,
            original_integer_and_incoming_history_checks=memory,pressure_datum_and_once_only_checks=datum,
            angular_defect_below_original_control_unit_checks=scale),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(signed_finite_N_correction_centers_resolved=False,finite_N_active_patch_uniform_inverse=False,
            patch_fields=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for n,d in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==d,n
    for n,d in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==d,n
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    r=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=r['checks'],execution_seconds=r['execution_seconds']),indent=2),flush=True)

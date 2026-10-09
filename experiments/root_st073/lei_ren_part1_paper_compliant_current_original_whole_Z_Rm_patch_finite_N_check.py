"""Check genuine Rm incoming, variable shear, signed drivers and Rh transport."""
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_finite_N as current
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow
from lei_ren_part1_paper_interval_taylor import IntervalTaylor


def encoded(value):
    return current.upstream.upstream.source.bridge._encode(current.serialized(value))


def finite(row):
    return row.ctx.mpf(0) if row.zero else row.coefficient*row.ctx.exp(row.scale.evaluate())


def contains(row,value):
    lo,hi=current.ep(finite(row))
    assert lo<=value<=hi


def independent_fixture():
    """Nonconstant angular shear and axial y derivative in finite source units."""
    c=MPIntervalContext();c.dps=180
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    zz=mp.mpf('.2');xx=mp.mpf(5)/4;tol=mp.mpf('1e-140')
    def jet(fn):
        return f.jet(IntervalTaylor(c,[c.mpf([mp.diff(fn,zz,n)/math.factorial(n)-tol,
            mp.diff(fn,zz,n)/math.factorial(n)+tol]) for n in range(6)]))
    eps=lambda z:(1+z+z**3)/1000
    am=lambda z:2+z/10
    H=lambda z:xx**mp.mpf('.1')+eps(z)*xx**3
    E=lambda z:am(z)*H(z)
    Ey=lambda z:am(z)*(xx**mp.mpf('.1')/10+3*eps(z)*xx**3)
    rv=lambda z:mp.mpf('.5')+z/10+xx**2*(mp.mpf('.03')+z/50)
    rvy=lambda z:2*xx**2*(mp.mpf('.03')+z/50)
    hist=dict(m=lambda z:1+z+z*z,h=lambda z:2-z+z**3/5,
        k=lambda z:3+2*z/5-z*z/5,e=lambda z:7*z*z/5-mp.mpf('.5'),
        p=lambda z:-z/5+z**3+mp.mpf('.3'))
    P0=lambda z:mp.mpf('.9')+3*z/10-z*z/5
    zj=IntervalTaylor.variable(c,c.mpf('.2'),5);zero=[f.scalar(0)]*6
    op=SimpleNamespace(flow=f,c=c,reference=SimpleNamespace(z=zj,delta=c.mpf('.2')),
        zrows=f.jet(zj),P0=jet(P0),Pstar=f.factor((0,.5,0,0,0)),Rm_factor=f.scalar(7),
        controls=[zero,zero,jet(eps),zero,zero])
    packet=dict(geometry=dict(point=True,exact_x=[5,4]),original_P0_normalized_axial5=op.P0,
        actual_H_x_derivative_axial5=[jet(H)],
        actual_gamma_ordinary_x_derivatives=[[c.mpf(xx)**3,3*c.mpf(xx)**2,6*c.mpf(xx),c.mpf(6),c.mpf(0)],
            [c.mpf(0)]*5,[c.mpf(0)]*5],
        raw_current_radius_y_derivative_axial_coefficients=dict(
            histories={key:[jet(fn),zero] for key,fn in hist.items()},
            velocity=dict(theta=[jet(E),jet(Ey)],axial=[jet(rv),jet(rvy)])))
    proxy,raw,recovered,a,shear=current.recover_patch(op,packet)
    roots,qr,proof=current.switch.general_quotients(proxy,recovered,a,c.mpf(-30))
    def expected(z):
        de=mp.mpf('.2');d=1-z*z;L=1-de*z*z
        V=rv(z)/3;Vy=rvy(z)/3
        m=hist['m'](z)/3;mz=mp.diff(hist['m'],z)/3
        h=hist['h'](z);hz=mp.diff(hist['h'],z)
        k=hist['k'](z)/3;kz=mp.diff(hist['k'],z)/3
        e=hist['e'](z);ez=mp.diff(hist['e'],z)
        pressure=P0(z)+hist['p'](z);pz=mp.diff(P0,z)+mp.diff(hist['p'],z)
        transport=m*z*(1-de)+mz*d
        itl=(-E(z)+h*(1-de/2)-hz*z*(1-de)/2)/L
        itq=(k*z*(2*de-1)-kz*d+E(z)*transport)/L
        izl=(-V+(m-mz*z)*(1-de)/2)/L
        izq=(V*transport+2*de*z*e-ez*d+2*(1+de)*z*pressure-pz*d)/L
        C=E(z)-2*Ey(z);B=2*Vy
        aa=C/E(z);b=B/E(z);t=-b/aa;kap=aa+b*b/aa;Delta=kap-2
        q=mp.sqrt((2*mp.exp(-30)-Delta)/(2*aa))
        return dict(E=E(z),C=C,B=B,V=V,pressure=pressure,a=aa,b=b,t0=t,
            kappa=kap,Delta=Delta,q=q,Q=(2*z*V-transport)/L,
            theta_linear=itl,theta_quadratic=itq,axial_linear=izl,axial_quadratic=izq,
            p2=(izl+3*izq)*7*xx/E(z))
    targets={key:recovered['actual_generic_source_numerators'][key] for key in ('E','C','B')}
    targets.update(V=recovered['common_velocity_V_axial5'],pressure=recovered['common_absolute_pressure_axial5'],
        Q=recovered['common_radial_Q_axial4'],**recovered['full_signed_inertial_sectors_axial4'])
    targets.update({key:proof['actual_'+key+'_axial5'] for key in ('a','b','t0','kappa','Delta')})
    targets['p2']=proof['full_signed_p2_axial4'];count=0
    for key,rows in targets.items():
        for n,row in enumerate(rows):
            contains(row,mp.diff(lambda z:expected(z)[key],zz,n)/math.factorial(n));count+=1
    for index,n in ((current.C0,0),(current.Z,1)):
        contains(qr[index],mp.diff(lambda z:expected(z)['q'],zz,n));count+=1
    assert abs(expected(zz)['a']-mp.mpf('.8'))>mp.mpf('.001') and expected(zz)['b']!=0
    assert shear['common_axial_amplitude_canceled_analytically'] and raw['original_P0_normalized_axial5'] is op.P0
    # Independent finite-N candidate formulas retain every quadratic cross term.
    funcs=dict(E=lambda z:2+z/10,V=lambda z:3+z/5,A=lambda z:mp.mpf('.2')+3*z/100,
        B=lambda z:mp.mpf('.4')-z/50)
    density=current.phase.densities.density_Z_kernels(f.scalar(funcs['E'](zz)),f.scalar(mp.mpf('.1')),
        f.scalar(funcs['V'](zz)),f.scalar(mp.mpf('.2')),
        dict(A=f.scalar(funcs['A'](zz)),A_Z=f.scalar(mp.mpf('.03')),
            B_over_Pstar=f.scalar(funcs['B'](zz)),B_Z_over_Pstar=f.scalar(-mp.mpf('.02'))),17)
    def kernels(z):
        e,v=funcs['E'](z),funcs['V'](z);de=e*mp.expm1(funcs['A'](z)/17);dv=funcs['B'](z)/17
        return dict(m=dv,h=de,k=(e+de)*(v+dv)-e*v,
            e=(v+dv)**2-v*v-((e+de)**2-e*e)/2,p=((e+de)**2-e*e)/2)
    for name in current.RATES:
        contains(density['kernels'][name],kernels(zz)[name])
        contains(density['Z_derivatives'][name],mp.diff(lambda z:kernels(z)[name],zz))
    masses=0
    for left,right in zip(current.PARTITION,current.PARTITION[1:]):
        lo=mp.log(mp.mpf(left[0])/left[1]);hi=mp.mpf(1) if right=='Rh' else mp.log(mp.mpf(right[0])/right[1])
        for rate in set(current.RATES.values()):
            width,suffix,mass,decay,tail=current.patch_weights(f,left,right,rate)
            lam=mp.mpf(rate.numerator)/rate.denominator
            contains(mass,mp.quad(lambda s:mp.exp(-lam*s),[0,hi-lo]))
            contains(decay,mp.exp(-lam*(hi-lo)));contains(tail,mp.exp(-lam*(1-hi)));masses+=3
    return dict(passed=True,nonconstant_shear_and_nonzero_b_generic_quotient_comparisons=count,
        independent_candidate_density_C0_Z_comparisons=10,independent_kernel_quadrature_and_memory_comparisons=masses,
        finite_diagnostic_source_units_only=True,actual_source_certificates_are_separate=True)


def run():
    began=time.monotonic();saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and saved['full_same_N_correction_prefix_through_Rh_installed']
    assert saved['actual_patch_finite_N_C0_Z_density_oracle_installed']
    assert not saved['actual_high_order_finite_N_correction_jets_installed']
    assert all(saved[key] is False for key in current.OPEN)
    for name,digest in saved['input_hashes'].items():assert current.sha(name)==digest,'Source changed: '+name
    counts=dict(whole_Z_cells=0,closed_radial_source_cells=0,signed_density_C0_Z_rows=0,
        own_rate_weight_rows=0,genuine_Rm_incoming_rows=0,actual_Rh_correction_rows=0,
        same_live_P0_queries=0,strict_actual_A_N_budgets=0,exact_directed_interval_rows=0,typed_rejections=0)
    with mp.workdps(540):independent=independent_fixture()
    with mp.workdps(540),patch.object(current.long,'OriginalLongReshapeFiniteN',forbidden), \
            patch.object(current.upstream.mixed,'OriginalRmPatchMixed4Cells',forbidden), \
            patch.object(current.prefix.WholeZLongRmFiniteN,'query',forbidden):
        owner=current.WholeZRmPatchFiniteN();c=owner.c
        assert saved['source_family']==owner.identity and saved['candidate_N']==owner.N
        assert owner.N==2**3981 and owner.N.bit_length()==3982
        assert saved['original_source_bindings']==encoded(owner.bindings)
        assert saved['exact_x_partition']==encoded(current.PARTITION)
        for ends,stored in zip(current.upstream.upstream.source.CELLS,saved['source_cells']):
            op=owner.owner(ends);f=op.flow
            live=owner.contribution(ends);same_source(encoded(live),stored,counts)
            assert live['exact_total_log_length']==1
            assert live['leading_controls_not_replaced_by_finite_N_correction']
            assert live['exact_common_P0_axial5'] is op.P0 and op.P0 is op.reference.P0
            inlet=live['actual_Rm_incoming_binding']
            assert inlet['correction_only_not_complete_own_history'] and inlet['no_saved_label_or_N257_density_transplant']
            counts['genuine_Rm_incoming_rows']+=10
            widthsum=c.mpf(0)
            for left,right,cell in zip(current.PARTITION,current.PARTITION[1:],live['actual_source_cells']):
                source=cell['source'];r=source['original_generic_source'];proof=source['original_full_source_quotients']
                assert source['exact_common_P0_axial5'] is op.P0 and r['common_original_P0_axial5'] is op.P0
                assert source['candidate_N']==owner.N and source['actual_A_below_same_candidate_N']
                assert r['actual_source_C_is_same_function_E_minus2Ey'] and r['source_ledger_is_same_object']
                assert source['actual_correlated_shear_source']['common_axial_amplitude_canceled_analytically']
                assert r['full_inertial_linear_quadratic_pressure_meridional_sectors_retained']
                assert r['derivative_orders']==dict(profile_and_histories=5,radial_and_full_inertial=4)
                assert source['actual_phase_Z_exact_zero'] and source['phase_average_cancellation_not_claimed']
                assert not source['actual_high_order_finite_N_correction_jets_installed']
                assert len(f.logs)==5 and r['original_fixed_source_log_bases'] is f.logs
                for rows in (*r['common_own_five_histories_axial5'].values(),
                        *r['actual_generic_source_numerators'].values(),*r['full_signed_inertial_sectors_axial4'].values()):
                    current.parameters.same_source(f,rows)
                current.parameters.same_source(f,[v for group in source['original_roots'].values() for v in group.values()])
                assert proof['one_actual_radius_factor'] is r['original_physical_radius']
                assert proof['actual_positive_E']['strict_positive_complete_source_enclosure']
                assert proof['actual_positive_a']['strict_positive_complete_source_enclosure']
                if source['actual_original_A_log_absolute_upper'] is not None:
                    assert current.ep(source['actual_original_A_log_absolute_upper'])[1]<current.ep(c.ln(owner.N))[0]
                counts['strict_actual_A_N_budgets']+=1;counts['same_live_P0_queries']+=1
                assert len(cell['signed_cell_driver_C0_Z'])==5
                for name,rate in current.RATES.items():
                    weight=cell['own_rate_weights'][name]
                    assert weight['own_rate']==str(rate) and weight['true_radial_measure_applied_once']
                    if rate==0:assert weight['true_cell_decay'].record()==f.scalar(1).record() and weight['true_suffix_decay'].record()==f.scalar(1).record()
                widthsum+=cell['actual_log_width']
                counts['closed_radial_source_cells']+=1;counts['signed_density_C0_Z_rows']+=10;counts['own_rate_weight_rows']+=15
            assert current.ep(widthsum)[0]<=1<=current.ep(widthsum)[1]
            assert live['incoming_own_rate_memory']['p'].record()==f.scalar(1).record()
            for name,rate in current.RATES.items():
                expected=f.scalar(1) if not rate else f.factor((0,0,0,0,0),-c.mpf(rate.numerator)/rate.denominator)
                assert live['incoming_own_rate_memory'][name].record()==expected.record()
            counts['actual_Rh_correction_rows']+=10;counts['whole_Z_cells']+=1
            print('Whole-Z actual same-N Rm..Rh correction audit: '+str(ends),flush=True)
        assert counts['closed_radial_source_cells']==44
        for call in (lambda:owner.owner(('0','0')),lambda:owner.query(('-1','-.5'),(0,1)),
                lambda:owner.query(('-1','-.5'),(2,1),(1,1)),lambda:owner.query(('-1','-.5'),(3,1))):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Invalid actual patch domain accepted')
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
        actual_variable_shear_generic_source_and_same_N_signed_C0_Z_drivers_checked=True,
        genuine_Rm_correction_only_incoming_true_own_rate_memory_and_Rh_exit_checked=True,
        independent_nonconstant_shear_and_density_kernel_fixture=encoded(independent),
        full_closed_real_axial_and_Rm_Rh_radial_cover_checked=True,
        actual_high_order_finite_N_correction_jets_installed=False,
        conservative_full_phase_covers_not_phase_average_cancellation=True,
        replay_counts=counts,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Whole-Z genuine same-N Rm..Rh correction checks passed',flush=True)
    return receipt


if __name__=='__main__':run()

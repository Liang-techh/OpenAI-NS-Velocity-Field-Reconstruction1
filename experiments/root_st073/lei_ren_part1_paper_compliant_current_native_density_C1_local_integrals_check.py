"""Independent full-velocity difference derivatives and native local C1 replay."""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_density_C1_local_integrals as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks
import lei_ren_part1_paper_compliant_current_native_q_slow_jets_check as qchecks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep;packets=current.packets
require=checks.require;contains=checks.contains;prior=current.prior;first=current.first;slow=current.slow


def independent_density_C1_checks(c):
    p=mp.mp.clone();p.dps=c.dps+20;count=0;integration_count=0
    for sign in (-1,1):
        bases=tuple(c.mpf(0) for _ in range(5));ledger=qchecks.new_ledger()
        scalar=lambda value:prior.ScaledEnclosure(prior.FormalScale(bases),value,ledger)
        Z=c.mpf(('-.01','.01'));N=1024
        E=scalar(2+Z*c.mpf('.07')+Z**2*c.mpf('.02'));E_Z=scalar(c.mpf('.07')+Z*c.mpf('.04'))
        V=scalar(-c.mpf('.3')+Z*c.mpf('.13'));V_Z=scalar('.13')
        A=scalar(sign*c.mpf('.03')+Z*c.mpf('.017'));A_Z=scalar('.017')
        B=scalar(-c.mpf('.02')+Z*c.mpf('.03')+Z**2*c.mpf('.011'));B_Z=scalar(c.mpf('.03')+Z*c.mpf('.022'))
        got=current.density_Z_kernels(E,E_Z,V,V_Z,dict(A=A,A_Z=A_Z,B_over_Pstar=B,B_Z_over_Pstar=B_Z),N)
        # Independent baseline: differences of the full modified products and
        # squares, without the implementation's expanded cross-term formulas.
        def reference(z,key):
            E=2+p.mpf('.07')*z+p.mpf('.02')*z*z;V=-p.mpf('.3')+p.mpf('.13')*z
            A=sign*p.mpf('.03')+p.mpf('.017')*z;B=-p.mpf('.02')+p.mpf('.03')*z+p.mpf('.011')*z*z
            EN=E*p.exp(A/N);VN=V+B/N
            return dict(m=VN-V,h=EN-E,k=EN*VN-E*V,
                e=VN*VN-V*V-(EN*EN-E*E)/2,p=(EN*EN-E*E)/2)[key]
        for point in ('-.01','-.005','0','.005','.01'):
            z=p.mpf(point)
            for key in current.RATES:
                for order,rows in ((0,got['kernels']),(1,got['Z_derivatives'])):
                    ref=reference(z,key) if not order else p.diff(lambda x:reference(x,key),z)
                    require(contains(rows[key].finite_interval(),c.mpf(p.nstr(ref,c.dps+15))),
                        'Independent full modified velocity difference not enclosed: '+str(sign)+' '+key+' Z'+str(order))
                    count+=1
    # Unmaterializable nonzero A/B must retain both increments and Z jets.
    bases=tuple(c.mpf(0) for _ in range(5));ledger=qchecks.new_ledger();scalar=lambda value:prior.ScaledEnclosure(prior.FormalScale(bases),value,ledger)
    A=prior.ScaledEnclosure(prior.FormalScale(bases,offset='-1e40'),1,ledger);B=A*3
    got=current.density_Z_kernels(scalar(2),scalar('.4'),scalar('.7'),scalar('-.1'),
        dict(A=A,A_Z=A*3,B_over_Pstar=B,B_Z_over_Pstar=B*4),1024)
    require(all(not value.zero for value in got['kernels'].values()) and all(not value.zero for value in got['Z_derivatives'].values()),
        'Tiny original density values/Z derivatives incorrectly replaced by zero')
    ratio=got['velocities']['deltaE_Z'].positive_divide(A,c.mpf('-1e40')).finite_interval()
    require(contains(ratio,c.mpf('6.4')/1024),'Factored tiny exponential Z chain lost its analytic leading coefficient')
    require(contains(got['velocities']['deltaV_Z'].positive_divide(B,c.mpf('-1e40')+c.ln(3)).finite_interval(),c.mpf(4)/1024),
        'Original normalized B Z factor/N applied incorrectly')
    zero=scalar(0);flat=current.density_Z_kernels(scalar(2),scalar('.4'),scalar('.7'),scalar('-.1'),
        dict(A=zero,A_Z=zero,B_over_Pstar=zero,B_Z_over_Pstar=zero),1024)
    require(all(v.zero for rows in (flat['kernels'],flat['Z_derivatives']) for v in rows.values()),
        'Exact flat primitive jets must give exact flat density changes')
    # Analytic finite integrals with a genuine varying Z source: f(s,Z)=
    # sign*(2+s)*(1+Z+Z²); compare both integral and its Z derivative.
    for rate in (Fraction(0),Fraction(1),Fraction(3,2)):
        w=p.mpf('.00000002');r=p.mpf(rate.numerator)/rate.denominator
        mass=current.local.positive_kernel_mass(c,c.mpf('.00000002'),rate)
        exact_mass=w if not r else -p.expm1(-r*w)/r
        first_moment=w*w/2 if not r else (1-(1+r*w)*p.exp(-r*w))/(r*r)
        linear=(2+w)*exact_mass-first_moment
        zbox=c.mpf(('.49','.51'));radial=c.mpf((2,ep(2+c.mpf('.00000002'))[1]))
        for sign in (-1,1):
            covers=(sign*radial*(1+zbox+zbox**2)*mass,sign*radial*(1+2*zbox)*mass)
            for z in (p.mpf('.49'),p.mpf('.5'),p.mpf('.51')):
                refs=(sign*linear*(1+z+z*z),sign*linear*(1+2*z))
                for covered,ref in zip(covers,refs):
                    require(contains(covered,c.mpf(p.nstr(ref,c.dps+15))),'Independent actual Z-dependent C1 integral not enclosed')
                    integration_count+=1
    return dict(passed=True,independent_full_modified_velocity_difference_C0_Z_comparisons=count,
        independent_analytic_Z_dependent_integral_C0_Z_comparisons=integration_count,
        positive_and_negative_exponent_and_original_nonzero_V_checked=True,
        factored_unmaterializable_density_and_Z_increments_not_zeroed=True,
        tiny_exponential_Z_chain_and_original_B_normalization_checked=True,
        exact_flat_density_and_Z_changes_zero=True,
        derivatives_built_from_independent_full_velocity_products_and_squares=True,
        scalar_fixtures_are_not_native_field_values=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_point_query_count']==5 and saved['native_local_integral_Z_query_count']==2,
        'Five actual density boxes plus two whole C1 integral cell queries required')
    require(not any(saved.get(key) for key in packets.OPEN),'Local C1 cannot complete global stages')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed native density C1 prerequisite: '+name)
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeDensityC1LocalIntegrals(first.NativePhaseFirstJets(slow.NativeQSlowJets(current.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(current.native.NativeGenericSourcePackets(bridge))))))
        independent=independent_density_C1_checks(owner.ctx);density_rows=0;density_Z_rows=0;integral_rows=0;integral_Z_rows=0;regions={}
        def inspect_cells(got):
            nonlocal density_rows,density_Z_rows
            for cell in got['cells']:
                require(cell['values'] is not None and cell['record']['status']=='enclosed','Actual signed density derivative source unresolved')
                values=cell['values'];require(set(values['kernels'])==set(current.RATES) and set(values['Z_derivatives'])==set(current.RATES),'All five original signed C0/Z kernels required')
                require(cell['record']['same_original_factor_basis_and_ledger'] and cell['record']['source_native_width_or_Pstar_conversion_not_reapplied'],
                    'Ordinary native axial derivatives need one normalization/basis')
                require(cell['record']['original_radius_phase_Z_derivative_exactly_zero'],'Actual phase Z derivative cannot be invented')
                base=values['velocities']['original_E']
                for rows in (values['kernels'],values['Z_derivatives'],values['velocities']):
                    for value in rows.values():require(value.scale.bases is base.scale.bases and value.ledger is base.ledger,'Density Z output basis/ledger changed')
                if cell['record']['original_phase_first_jet_source']['geometry']=='flat':
                    require(all(value.zero for rows in (values['kernels'],values['Z_derivatives']) for value in rows.values()),'Native flat density changes/Z rows must be exact zero')
                require(not any(cell['record'].get(key) for key in packets.OPEN),'Density source cannot admit global stages')
                density_rows+=5;density_Z_rows+=5
        for chart,old in saved['native_spatial_signed_density_Z_records'].items():
            geometry=old['actual_original_radius_phase'];prov=old['source_provenance'];exact=geometry['declared_exact_coordinate']
            coordinate={'original_power_offset':'.537'} if chart=='O3_power' else str(exact['numerator'])+'/'+str(exact['denominator'])
            got=owner.spatial_query(chart,packets.interval(owner.ctx,prov['Z_box']),coordinate,saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual native signed density Z replay changed: '+chart);inspect_cells(got)
            if chart=='O2_buffer':require(all(not value.zero for value in got['cells'][0]['values']['kernels'].values()),'Tiny actual buffer signed densities lost')
            regions[chart]=dict(actual_phase_cells=len(got['cells']),signed_density_Z_rows=5*len(got['cells']))
            print('Actual original signed density Z checked:',chart,flush=True)
        signs=dict(m='negative',h='negative',k='negative',e='positive',p='negative')
        for name,old in saved['actual_local_C1_signed_integral_records'].items():
            got=owner.contribution(Z=packets.interval(owner.ctx,old['Z_box']),left='.13369999',right='.13370001',N=saved['candidate_N'])
            require(packets.encode(got['record'])==old,'Actual local C1 integral replay changed: '+name)
            require(Fraction(**old['exact_log_radius_width_fraction'])==Fraction(1,50000000),'Actual Z-independent O2 width changed')
            require(old['Z_independent_endpoints_and_weights_no_boundary_terms'] and old['actual_local_C1_Z_integral_functions_installed'],
                'Local C1 requires actual derivative-under-integral source/geometry theorem')
            require(old['local_contributions_are_not_global_defect_histories'] and old['original_incoming_history_not_assumed_or_reset'],
                'Local C1 cannot reset or invent incoming histories')
            require(not old['global_C1_histories_or_Rc_targets_or_repair_admitted'] and not any(old.get(key) for key in packets.OPEN),
                'Local C1 cannot admit global Rc/repair/cone/recursion')
            require(set(got['Z_derivatives'])==set(current.RATES),'Five actual local integral Z functions required')
            for key,value in got['contributions'].items():
                require(not value.zero and value.record()['sign']==signs[key],'Actual signed C0 integral lost: '+key)
                require(ep(got['record']['original_positive_Duhamel_kernel_mass_covers'][key])[0]>0,'Original actual positive mass required')
            # Inspect source cells via the exact returned C1 record. The local
            # query's C0 and derivative hulls are distinct, never differentiated hull endpoints.
            density_rows+=5;density_Z_rows+=5;integral_rows+=5;integral_Z_rows+=5
            regions[name]=dict(Z_box=old['Z_box'],local_integral_Z_rows=5,
                derivative_enclosure_signs={key:value.record()['sign'] for key,value in got['Z_derivatives'].items()})
            print('Actual local C1 signed integral checked:',name,flush=True)
        for left,right,chart in (('.2','.1','O2_slope'),(owner.ctx.mpf('.1'),'.2','O2_slope'),('.1','.2','reshape')):
            try:owner.contribution(Z=('.5','.5'),left=left,right=right,N=1024,chart=chart)
            except ValueError:continue
            raise ArithmeticError('Invalid/unsupported local C1 geometry accepted')
    require((density_rows,density_Z_rows,integral_rows,integral_Z_rows)==(35,35,10,10),'Declared native signed C1 row scope changed')
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},candidate_N=saved['candidate_N'],
        native_source_query_count=7,original_signed_density_C0_rows_checked=density_rows,original_signed_density_Z_rows_checked=density_Z_rows,
        actual_local_signed_integral_C0_rows_checked=integral_rows,actual_local_signed_integral_Z_rows_checked=integral_Z_rows,
        actual_whole_radial_Z_C1_cell_checked=True,all_ten_original_C0_local_integrals_nonzero=True,
        independent_original_density_and_integral_C1_checks=independent,regions=regions,
        global_C1_histories_or_Rc_or_common_N_or_repair_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual signed density/local C1 PASS:',density_Z_rows,'density Z rows;',integral_Z_rows,'integral Z rows',flush=True)
    return result


if __name__=='__main__':run()

"""Independent physical bump quadrature and one actual control-map replay."""
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rc_control_operator as current
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

encoded = lambda v:current.current.encode(current.current.serialized(v))


def independent_physical_controls():
    c = MPIntervalContext(); c.dps=150
    p = mp.mp.clone(); p.dps=190
    f = MacroFlow(c,-30,c.ln(4),c.ln(9),c.ln(5),c.ln(100)-c.ln(5),'.01')
    comparisons=0
    for mu_text in ('.02','.1'):
        mu=p.mpf(mu_text); N=512; L=p.ln(2); ell=L/40
        centers=[L/5,L/2,4*L/5]
        raw=lambda t:p.exp(-1/(1-t*t)) if abs(t)<1 else p.mpf(0)
        normal=p.quad(raw,[-1,0,1])
        def integrate(i,fun):
            def value(t):
                x=p.exp(centers[i]+ell*t); g=raw(t)/(ell*normal*x)
                return fun(x,g)*ell*x
            return p.quad(value,[-1,0,1])
        B=p.matrix(5,5)
        for j,i in enumerate((0,2)):
            B[0,j]=integrate(i,lambda x,g:g)
            B[1,j]=integrate(i,lambda x,g:p.expm1(-mu*p.ln(x))*g/mu)
        for i in range(3):
            B[2,i+2]=integrate(i,lambda x,g:p.sqrt(x)*g)
            B[3,i+2]=integrate(i,lambda x,g:-x**(-p.mpf('.5')-mu)*g)
            B[4,i+2]=integrate(i,lambda x,g:x**(-p.mpf('1.5')-mu)*g)
        inverse=B**-1
        cross=[integrate(i,lambda x,g:p.sqrt(x)*g*g) for i in (0,2)]
        energy=[integrate(i,lambda x,g:g*g) for i in range(3)]
        pressure=[integrate(i,lambda x,g:g*g/x) for i in range(3)]
        def Q(h):
            a,b,e,g,k=h
            return p.matrix([0,(cross[0]*a*e+cross[1]*b*k)/mu,0,
                energy[0]*a*a+energy[2]*b*b-(energy[0]*e*e+energy[1]*g*g+energy[2]*k*k)/2,
                (pressure[0]*e*e+pressure[1]*g*g+pressure[2]*k*k)/2])
        target=lambda z:p.matrix([p.mpf('.01')*(i+1)+p.mpf('.003')*(i-2)*z+p.mpf('.002')*z*z for i in range(5)])
        def exact(z,j):
            h=p.matrix([0]*5)
            for unused in range(j):h=-inverse*(target(z)+Q(h)/N)
            return h
        z=p.mpf('.37')
        row=lambda v:f.scalar(c.mpf(p.nstr(v,185)))
        d=[current.Pair(row(target(z)[i]),row(p.diff(lambda v:target(v)[i],z))) for i in range(5)]
        W=current.repair.fresh_weights(c,c.mpf(mu_text),cells=512)
        matrix=current.repair.fresh_linear_inverse(c,c.mpf(mu_text),W)
        sequence,residual=current.picard_ranges(c,c.mpf(mu_text),c.ln(c.mpf(mu_text)),N,matrix,d,3)
        for j,controls in enumerate(sequence):
            for i,q in enumerate(controls):
                for order,value in enumerate((q.value,q.Z)):
                    want=p.diff(lambda v:exact(v,j)[i],z,order)
                    lo,hi=current.ep(value.finite_interval(max_log=1000))
                    assert lo <= want <= hi,(mu_text,j,i,order)
                    comparisons+=1
        for i,q in enumerate(residual):
            for order,value in enumerate((q.value,q.Z)):
                fun=lambda v:(target(v)+B*exact(v,3)+Q(exact(v,3))/N)[i]
                want=p.diff(fun,z,order)
                stable=p.diff(lambda v:((Q(exact(v,3))-Q(exact(v,2)))/N)[i],z,order)
                lo,hi=current.ep(value.finite_interval(max_log=1000))
                assert lo <= stable <= hi,(mu_text,'stable_residual',i,order)
                # The physical reference uses a 190-digit numerical matrix
                # inverse. Direct B*h+d cancellation can leave roundoff in
                # rows whose mathematical residual is exactly zero. This
                # tolerance belongs only to the synthetic scalar reference;
                # no actual-source interval is enlarged.
                reference_roundoff=p.mpf('1e-170')*(1+abs(stable))
                assert abs(want-stable) <= reference_roundoff,(mu_text,'direct_residual_roundoff',i,order)
                comparisons+=2
    return dict(passed=True,independent_physical_bump_quadrature_control_Z_residual_comparisons=comparisons,
        physical_x_measures_and_all_quadratic_cross_terms_used=True,
        stable_physical_residual_strictly_enclosed=True,
        synthetic_190_digit_direct_residual_roundoff_tolerance='1e-170*(1+abs(reference))',
        nonzero_polynomial_target_Z=True,synthetic_reference_not_actual_source_claim=True)


def forbidden(*args,**kwargs):raise AssertionError('Unchanged source integrations must not be repeated')


def run():
    began=time.monotonic()
    saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name,digest in saved['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed source prerequisite: '+name)
    with mp.workdps(540):
        fixture=independent_physical_controls()
        owner=current.WholeZRcControlOperator()
        assert owner.graph_digest==saved['source_binding_sha256'] and owner.identity==saved['source_family']
        nodes=saved['exact_function_graph_nodes']
        built=owner.build()
        assert encoded(built['graph'].nodes)==nodes
        roots=lambda pairs,labels:{key:dict(value=q.value.node,Z=q.Z.node) for key,q in zip(labels,pairs)}
        assert [roots(h,current.CONTROLS) for h in built['finite_picard_sequence']] \
            ==saved['exact_finite_control_roots']
        assert roots(built['control_residual'],current.ROWS)==saved['exact_control_residual_roots']
        assert {k:dict(value=v.value.node,Z=v.Z.node) for k,v in built['band_profiles'].items()} \
            ==saved['exact_trial_band_roots']
        sources=[row for row in nodes if row['operation']=='actual_whole_Z_Rc_source_function']
        assert len(sources)==12 and all(row['source_binding_sha256']==owner.graph_digest for row in sources)
        assert all(row['numerical_value_not_installed'] for row in sources)
        counts=0
        for row in saved['actual_four_Z_control_map_cells']:
            source=owner.cells[tuple(row['exact_Z_cell'])]
            assert row['actual_N_scaled_target_C0_Z']==source['actual_fixed_N_N_scaled_target_C0_Z']
            assert row['exact_common_P0_axial5']==source['exact_common_P0_axial5']
            assert row['candidate_N']==owner.N and not row['actual_fixed_N_solution_certified']
            assert len(row['finite_control_sequence_C0_Z'])==4
            for x in ('1','2'):
                band=row['original_supported_trial_band_C0_Z'][x]
                assert all(v['exact_zero'] for key in ('F','G','V') for v in band[key])
                assert band['E']==band['leading_E']
            counts+=5
        with patch.object(current.current.WholeZSharpBridgeFunctions,'retransport',forbidden):
            live=owner.controls(current.CELLS[0],mode='accepted')
            assert encoded([current.pairs_record(h,current.CONTROLS) for h in live['sequence']]) \
                ==saved['actual_four_Z_control_map_cells'][0]['finite_control_sequence_C0_Z']
            assert encoded(current.pairs_record(live['residual'],current.ROWS)) \
                ==saved['actual_four_Z_control_map_cells'][0]['actual_finite_iterate_residual_C0_Z']
        for bad in (('2','3'),):
            try:owner.controls(bad,mode='accepted')
            except ValueError:pass
            else:raise AssertionError('Unadmitted source cell accepted')
        wrong=dict(live,candidate_N=owner.N+1)
        try:owner.band(wrong,'1.5')
        except ValueError:pass
        else:raise AssertionError('Different N band accepted')
        for flag in (*current.OPEN,'actual_terminal_controls_installed','actual_fixed_N_solution_certified',
                'actual_global_frequency_admitted','physical_original_exterior_five_targets_closed'):
            assert saved[flag] is False
        hashes=dict(saved['input_hashes'])
        for name in (current.NAME,Path(__file__).name):current.bind(hashes,name,current.sha(name))
        receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
            independent_physical_control_fixture=fixture,exact_actual_target_pairs_bound=counts,
            exact_source_function_nodes=12,actual_first_Z_three_iterates_and_residual_replayed=True,
            exact_source_bound_graph_and_control_residual_band_roots_rebuilt=True,
            all_four_Z_exact_supported_band_endpoint_C0_Z_joins=True,
            unchanged_source_integrations_not_repeated=True,actual_terminal_controls_installed=False,
            actual_fixed_N_solution_certified=False,actual_global_frequency_admitted=False,
            physical_original_exterior_five_targets_closed=False,**dict.fromkeys(current.OPEN,False),
            input_hashes=hashes,execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.current.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS: actual whole-Z Rc finite control functions',flush=True)
    return receipt


if __name__=='__main__':run()

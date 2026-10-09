"""Independent finite nonlinear bridge integrals and native source contracts."""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_bridge_macro_functions as current
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

cbase=current.base;ep=current.ep


def finite(row):
    c=row.ctx
    return row.coefficient*c.exp(row.scale.evaluate()) if not row.zero else c.mpf(0)


def contain(row,value):
    lo,hi=ep(finite(row));assert lo<=value<=hi,(mp.nstr(lo,16),mp.nstr(value,16),mp.nstr(hi,16))


def independent_checks():
    c=MPIntervalContext();c.dps=100
    z=sy.symbols('z');a=sy.symbols('a1:6')
    coefficients=sy.series(sy.exp(sum(a[n-1]*z**n for n in range(1,6))),z,0,6).removeO()
    bell=[sy.lambdify(a,coefficients.coeff(z,n),'mpmath') for n in range(6)]
    kernel_checks=field_checks=mass_checks=0;records=[]
    with mp.workdps(130):
        for sign,htext in ((1,'.0001'),(-1,'.0003')):
            h=mp.mpf(htext);Ra=mp.mpf('.7');Y=c.ln(100)-c.ln(c.mpf('.7'))
            flow=current.MacroFlow(c,c.ln(c.mpf(htext)),c.mpf('.3'),c.mpf('-.2'),c.ln(c.mpf('.7')),Y,c.mpf('.1'))
            def jet(rows):return IntervalTaylor(c,[c.mpf(str(v)) for v in rows])
            directions=[jet([sign*mp.mpf('.002')/(j+1),'.0004','-.0003','.0002','-.0001','.00005']) for j in range(3)]
            drives={part:[jet([mp.mpf('.02')*(j+1),'.006','-.004','.003','-.002','.001']) for j in range(3)] for part in current.PARTS}
            q=jet(['1.1','.03','-.02','.01','-.005','.002']);phi=jet(['.8','.02','-.01','.005','-.002','.001']);V0=jet(['.3','.04','-.02','.01','-.004','.002'])
            lin=jet([h*mp.mpf('.01'),h*mp.mpf('.004'),h*mp.mpf('-.003'),h*mp.mpf('.002'),h*mp.mpf('-.001'),h*mp.mpf('.0005')])
            vin=jet([h*mp.mpf('.02')]+[h*mp.mpf('.001')]*5)
            flow.set_sources(directions,drives,q,flow.jet(lin),phi,V0,flow.jet(vin))
            def points(jet):return [sum(ep(v))/2 for v in jet.coefficients]
            dr=[points(v) for v in directions];cr={p:[points(v) for v in rows] for p,rows in drives.items()}
            qr,phr,vr,lr,vir=map(points,(q,phi,V0,lin,vin))
            for fraction in ((1,2),(1,1)):
                value=flow.evaluate(fraction);S,R0box,R1box,empty=flow.geometry(fraction)
                qs=mp.mpf(fraction[0])/fraction[1];R0=Ra*mp.exp(2*h);Strue=qs*(mp.log(100/Ra)-2*h)
                def primitive(t,l):return R0*t if l==1 else R0*mp.expm1((1-l)*t)/(1-l)
                def ellrows(t):return [lr[n]-h*sum(dr[l][n]*primitive(t,l) for l in range(3))/2 for n in range(6)]
                def erows(t):
                    e=ellrows(t);return [mp.exp(e[0])*fn(*e[1:]) for fn in bell]
                cut=[0,Strue/3,2*Strue/3,Strue]
                for p in (1,2):
                    for j in range(3):
                        I=mp.quad(lambda t:R0**p*mp.exp((p-j)*t),cut)
                        contain(flow.I(p,j,S,R0box,R1box),I);mass_checks+=1
                        for l in range(3):
                            J=mp.quad(lambda t:R0**p*mp.exp((p-j)*t)*primitive(t,l),cut)
                            contain(flow.J(p,j,l,S,R0box,R1box),J);mass_checks+=1
                        part='swirl' if p==2 else 'hydro'
                        for n in range(6):
                            exact=mp.quad(lambda t:R0**p*mp.exp((p-j)*t)*erows(t)[n],cut)
                            contain(value['kernels'][part][j]['coefficients'][n],exact);kernel_checks+=1
                # Independently integrate the full product, including q and
                # each driving source jet, before comparing its V coefficient.
                dv=vir[:]
                for part,p,scale in (('hydro',1,1),('pressure',1,mp.exp(mp.mpf('.3'))),('swirl',2,mp.exp(mp.mpf('-.2')))):
                    product=[]
                    for n in range(6):
                        def integrand(t):
                            e=erows(t)
                            return sum(R0**p*mp.exp((p-j)*t)*sum(qr[k]*cr[part][j][l]*e[n-k-l]
                                for k in range(n+1) for l in range(n-k+1)) for j in range(3))
                        product.append(-h*scale*mp.quad(integrand,cut))
                        contain(value['parts'][part][n],product[-1]);field_checks+=1
                        dv[n]+=product[-1]
                e=erows(Strue)
                for n in range(6):
                    contain(value['V'][n],vr[n]+dv[n]);field_checks+=1
                    exact=sum(phr[k]*e[n-k] for k in range(n+1))
                    contain(value['phi'][n],exact);field_checks+=1
                records.append(dict(sign=sign,h=htext,fraction=list(fraction),complete_original_nonlinear_integral_checked=True,diagnostic_fixture_only=True))
    return dict(passed=True,independent_complete_nonlinear_kernel_Taylor_comparisons=kernel_checks,
        independent_full_FV_source_product_comparisons=field_checks,
        independent_resonant_and_nonresonant_IJ_mass_quadratures=mass_checks,records=records)


def native_checks(owner):
    c=owner.c;derivative_rows=errors=scale_checks=inlet_checks=0
    for label in ('0','.5'):
        flow,proof=owner.owner(label)
        assert proof['actual_G_not_Gbar_used'] and proof['normalized_swirl_rows_include_true_F0_squared_derivatives_once']
        assert proof['original_micro_inlet_errors_retained'] and proof['no_original_ancestor_constructors_or_producers_executed']
        amplitude=owner.records['anchored_axis_amplitude']['anchored_amplitude_packets'][label]
        native_log=2*current.previous.read_interval(c,amplitude['logF0'])
        assert native_log._mpi_==flow.logs[2]._mpi_;scale_checks+=1
        assert proof['comparison_phi']['directed_nonzero_remainder_norm']['exact_zero'] is False
        assert proof['comparison_V']['directed_nonzero_remainder_norm']['exact_zero'] is False
        # Pressure is the same original14-atom datum, with a separate P0 and
        # every true derivative retained by the original cached recipe.
        original=owner.records['actual_bridge_integrals']['packets'][label]['core']['pressure_axis_axial5']
        for n,v in enumerate(original):
            lo,hi=ep(current.previous.read_interval(c,v));p=proof['original_P0_coefficients'][n]
            assert max(lo,ep(p)[0])<=min(hi,ep(p)[1]);inlet_checks+=1
        for q in ((0,1),(1,2),(1,1)):
            result=flow.evaluate(q)
            assert result['full_original_macro_exponential_integral_enclosed']
            assert result['source_contract']=='q=phi_core_exit/barphi2; ell includes micro angular prefix exactly once'
            assert result['geometry']['fixed_radial_endpoints_Z_independent']
            for key in ('ell','delta_phi','deltaV','phi','V'):
                for row in result[key]:
                    assert row.scale.bases is flow.logs and row.ledger is flow.ledger
                    assert not row.record()['point_value_selected'];derivative_rows+=1
            if q[0]:
                for part in current.PARTS:
                    bound=result['macro_absolute_error_norms'][part]
                    assert not bound.zero and ep(bound.coefficient)[0]>=0;errors+=1
                if q[0]==q[1]:assert ep(finite(result['geometry']['R1']))==(100,100)
            else:
                for n,row in enumerate(result['deltaV']):
                    # Exact zero macro-length leaves the source micro inlet;
                    # do not replace it by an interval selector or reset it.
                    assert row.record()==flow.deltaV_in[n].record()
    return dict(passed=True,genuine_source_FV_and_log_Taylor_records=derivative_rows,
        nonzero_full_nonlinear_macro_error_norms=errors,actual_anchored_scale_identity_checks=scale_checks,
        same_original_pressure_coefficient_overlaps=inlet_checks,micro_inlet_not_reset=True,
        full_bridge_moments_first_switch_and_whole_axis_remain_open=True)


def guards(owner):
    flow,_=owner.owner('0');count=0
    def reject(fn):
        nonlocal count
        try:fn()
        except (ValueError,ArithmeticError):count+=1;return
        raise AssertionError('Invalid original function request accepted')
    for q in ((-1,1),(2,1),(1,0),(True,1),[1,2],(.5,1)):reject(lambda q=q:flow.evaluate(q))
    reject(lambda:owner.owner('anchor'))
    reject(lambda:current.OriginalBridgeMacroFunctions(dps=50))
    c=MPIntervalContext();c.dps=80
    reject(lambda:current.MacroFlow(c,-10,0,0,0,1,.1))
    reject(lambda:flow.set_sources([],{},None,[],None,None,[]))
    assert count==10
    return dict(passed=True,invalid_fraction_source_context_and_precision_requests_rejected=count)


def run():
    began=time.monotonic();owner=current.OriginalBridgeMacroFunctions()
    producer=json.loads((current.HERE/current.NAME).read_bytes())
    assert producer[current.GATE] and producer['source_family']==owner.family
    for name in current.previous.OPEN:assert producer[name] is False
    for path,digest in producer['input_hashes'].items():assert current.sha(path)==digest
    refs=independent_checks();print('Independent complete nonlinear macro F/V integrals PASS',flush=True)
    native=native_checks(owner);bad=guards(owner)
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name,'lei_ren_part1_paper_logarithmic_pressure_datum.py'):
        current.previous.bind(hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_complete_original_flow=refs,genuine_native_source_functions=native,guards=bad,
        **dict.fromkeys(current.previous.OPEN,False),input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Original frozen macro complete nonlinear F/V integral functions at admitted native0,.5 frames, with source micro inlet errors. Full bridge moments, switch, controls/global N and recursion remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(cbase.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original frozen macro source functions and nonzero errors PASS',flush=True);return result


if __name__=='__main__':run()

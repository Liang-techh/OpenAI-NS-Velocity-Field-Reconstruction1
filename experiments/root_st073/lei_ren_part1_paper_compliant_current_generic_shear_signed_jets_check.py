"""Signed quotient DAG identities, actual source leaves and independent fixtures."""
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_shear_signed_jets as source

packets=source.packets


def independent_fixture():
    """Modest signed multimode sources, independent symbolic mixed derivatives."""
    c=MPIntervalContext();c.dps=110;ep=packets.recovery.endpoints
    algebra=packets.FactoredAlgebra(c,(c.mpf(-2),c.ln(9),c.ln(4),c.mpf(-1)),[])
    y,Z=s.symbols('y Z');R0=s.Rational(7,3)
    E1=(1+Z*Z/10)*s.exp(y/10);E2=(s.Rational(1,10)+Z/50)*s.exp(-y/20)
    VV=(s.Rational(3,10)+Z/10)*s.exp(-3*y/20)
    It=(s.Rational(2,5)+3*Z/10+Z*Z/20)*s.exp(y/10)
    Itq=Z*s.exp(-y/5)/20
    Iz=(s.Rational(4,5)-Z/5+y*Z/50)*s.exp(y/20)
    Izq=(s.Rational(1,5)-Z*Z/10)*s.exp(y/5)
    mode1=(0,0,0,.5);mode2=(-1,0,.5,0);v_mode=(0,-.5,0,0);linear_mode=(-1,0,.5,.5)
    def weight(mode):return s.exp(sum(s.Rational(str(p))*q for p,q in zip(mode,(-2,s.log(9),s.log(4),-1))))
    ee=E1*weight(mode1)+E2*weight(mode2);vv=VV*weight(v_mode)
    pp1=R0*s.exp(y)*(It*weight(linear_mode)+3*Itq)/ee
    pp2=R0*s.exp(y)*(Iz*weight(linear_mode)+3*Izq)/ee
    aa=1-2*s.diff(ee,y)/ee;bb=2*s.diff(vv,y)/ee
    refs=dict(E=ee,a=aa,b=bb,p1=pp1,p2=pp2,t0=-bb/aa,kappa_minus2=aa+bb*bb/aa-2)
    count=0
    for yp,zp in ((s.Rational(0),s.Rational(1,5)),(s.Rational(1,50),-s.Rational(3,10)),(-s.Rational(1,50),s.Rational(2,5))):
        def row(expr,order,mode,physical=False):
            yy=s.diff(s.exp(y/2)*expr,y,order)/s.exp(y/2) if physical else s.diff(expr,y,order)
            coeffs=[c.mpf(str(s.N(s.diff(yy,Z,k).subs({y:yp,Z:zp})/math.factorial(k),105))) for k in range(6)]
            return packets.FactoredJet(algebra,{mode:packets.IntervalTaylor(c,coeffs)},5)
        Erows=tuple(row(E1,j,mode1)+row(E2,j,mode2) for j in range(5))
        Vrows=tuple(row(VV,j,v_mode) for j in range(5))
        full={}
        for component,lin,quad in (('theta',It,Itq),('axial',Iz,Izq)):
            full['inertial_'+component+'_linear']=tuple(row(lin,j,linear_mode,True) for j in range(3))
            full['inertial_'+component+'_quadratic']=tuple(row(quad,j,(0,0,0,0),True) for j in range(3))
        packet=SimpleNamespace(chart='fixture',source_family={'fixture':True},algebra=algebra,
            velocity=dict(theta=Erows,axial=Vrows),provenance=dict(original_radius_source={'fixture_R0':str(R0)}),
            recover_original=lambda delta:dict(full_signed_stress_ordinary_y_rows=full))
        result=source.from_packet(packet,0,dict(log_E_positive_lower=c.mpf(-3),log_C_positive_lower=c.mpf(-3),log_actual_a_positive_lower=c.mpf(-3)))
        record=result['jet_expression_dag'];pool=[]
        radius=c.mpf(str(s.N(R0*s.exp(yp),105)))
        def leaf_value(name):
            total=c.mpf(0)
            for term in record['source_derivative_leaves'][name]['original_radius_power_terms']:
                modal=term['coefficient'];value=c.mpf(0)
                for mode,jet in modal.terms.items():
                    value+=jet[0]*c.exp(sum((base*power for base,power in zip(algebra.logs,mode)),c.mpf(0)))
                total+=radius**term['radius_power']*value
            return total
        for node in record['nodes']:
            op=node['operation']
            if op=='constant':value=c.mpf(node['value'])
            elif op=='source_derivative':value=leaf_value(node['name'])
            elif op=='negative':value=-pool[node['argument']]
            elif op=='sum':value=sum((pool[k] for k in node['arguments']),c.mpf(0))
            elif op=='product':
                value=c.mpf(1)
                for k in node['arguments']:value*=pool[k]
            elif op=='positive_function_quotient':
                denominator=pool[node['denominator']]
                if ep(denominator)[0]<=0:raise ArithmeticError('Independent fixture denominator lost positivity')
                value=pool[node['numerator']]/denominator
            else:raise ValueError('Unknown signed expression operation')
            pool.append(value)
        for key,expr in refs.items():
            for j,k in source.ORDERS:
                reference=mp.mpf(str(s.N(s.diff(expr,y,j,Z,k).subs({y:yp,Z:zp}),100)))
                actual=pool[record['roots'][key]['y%d_Z%d'%(j,k)]];lo,hi=ep(actual)
                allowance=mp.mpf('1e-85')*(1+abs(reference))
                if reference<lo-allowance or reference>hi+allowance:
                    raise ArithmeticError('Independent signed mixed jet differs: '+key+str((j,k)))
                count+=1
    if algebra.proofs or algebra.final_rows:raise ArithmeticError('Fixture changed source resolution ledger')
    return dict(comparisons=count,signed_nonzero_axial_and_quadratic_sources=True,
        four_source_log_bases_and_single_R_checked=True,
        production_source_factor_materialization=False,
        tolerance='1e-85 relative/absolute in modest symbolic fixture only')


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in manifest['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed signed source prerequisite: '+name)
    saved=json.loads(gzip.decompress((source.HERE/source.VIEWS).read_bytes()))
    provider=source.CurrentSignedInputJets();count=0;coefficients=0;quotients=0
    if set(saved)!=set(provider.inventory) or len(saved)!=17:raise ArithmeticError('Signed positive source cover incomplete')
    for chart,record in saved.items():
        actual=provider.saved(chart)
        if packets.encode(actual)!=record:raise ValueError('Original signed source-jet record differs: '+chart)
        dag=actual['jet_expression_dag'];nodes=dag['nodes']
        for index,node in enumerate(nodes):
            dependencies=[]
            if node['operation'] in ('sum','product'):dependencies=node['arguments']
            elif node['operation']=='negative':dependencies=[node['argument']]
            elif node['operation']=='positive_function_quotient':
                dependencies=[node['numerator'],node['denominator']];quotients+=1
                if node['positive_certificate'] not in dag['positive_denominator_certificates']:
                    raise ArithmeticError('Uncertified signed denominator')
            elif node['operation']=='source_derivative':
                if node['name'] not in dag['source_derivative_leaves']:raise ArithmeticError('Unbound original source derivative leaf')
            elif node['operation']!='constant':raise ArithmeticError('Unknown signed arithmetic operation')
            if any(type(k) is not int or not 0<=k<index for k in dependencies):raise ArithmeticError('Signed derivative DAG not acyclic')
        for key,leaf in dag['source_derivative_leaves'].items():
            expected_power=1 if key.startswith(('nt_','nz_')) else 0
            for term in leaf['original_radius_power_terms']:
                if term['radius_power']!=expected_power:raise ArithmeticError('Original radius factor applied incorrectly')
                row=term['coefficient']
                if row.order!=0:raise ArithmeticError('Ordinary derivative leaf must be coefficient order0')
                for values in row.terms.values():
                    for q in values.coefficients:
                        if any(not mp.isfinite(v) for v in packets.recovery.endpoints(q)):raise ArithmeticError('Signed coefficient cover nonfinite')
                        coefficients+=1
        for key,roots in dag['roots'].items():
            if set(roots)!={'y%d_Z%d'%k for k in source.ORDERS}:raise ArithmeticError('Signed derivative order incomplete')
            if any(not 0<=k<len(nodes) for k in roots.values()):raise ArithmeticError('Unbound signed input root')
            count+=len(roots)
    bad_scope=json.loads((source.HERE/source.current.RECEIPT).read_bytes());bad_scope[source.OPEN[0]]=True
    guards=(lambda:source.validate_source_receipt(bad_scope,provider.family),
        lambda:provider.saved('core'),lambda:provider.query('bridge_first',0,'.5'),
        lambda:provider.query('O3_slope_mu',0,'.5'),lambda:provider.query('O3_power',0,'.5'),
        lambda:source.CurrentSignedInputJets(O3_owners={'unknown':object()}))
    for call in guards:
        try:call()
        except ValueError:pass
        else:raise ArithmeticError('Invalid cover/live-owner source query admitted')
    if count!=714 or any(manifest[k] for k in source.OPEN):raise ArithmeticError('Signed expressions promoted missing physical stages')
    if packets.encode(provider.theorem)!=manifest['exact_signed_source_jet_theorem']:raise ValueError('Signed source derivative theorem differs')
    with mp.workdps(120):fixture=independent_fixture()
    result=dict(all_passed=True,source_family=provider.family,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
        source_cover_charts=17,signed_input_jet_roots=count,signed_original_leaf_coefficients=coefficients,
        deferred_positive_function_quotient_nodes=quotients,
        exact_signed_mixed_quotient_and_radius_identities=len(provider.theorem['exact_signed_ordinary_quotient_and_radius_identities']),
        independent_signed_source_fixture=fixture,invalid_saved_or_live_owner_guards=len(guards),
        ordinary_axial_factorials_original_y_rows_and_one_R_preserved=True,
        source_cover_enclosures_not_reclassified_as_point_field_values=True,
        actual_radius_width_amplitude_or_denominator_not_evaluated=True,
        successful_arbitrary_coordinate_live_owner_query_tested=False,
        signed_current_point_loop_or_inverse_jets_installed=False,
        source_graph_ancestor_constructors_called=False,
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Signed source derivative DAG PASS:714 roots,17 charts; independent signed multimode fixtures PASS',flush=True)
    return result


if __name__=='__main__':run()

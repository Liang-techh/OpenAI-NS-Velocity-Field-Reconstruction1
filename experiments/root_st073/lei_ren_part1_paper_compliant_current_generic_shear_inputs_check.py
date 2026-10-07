"""Check full physical native stress conversion and source invariant recipes."""
import gzip
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_shear_inputs as source
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows


def overlap(a,b,label):
    al,ah=source.packets.recovery.endpoints(a);bl,bh=source.packets.recovery.endpoints(b)
    if ah<bl or bh<al:raise ArithmeticError('Full original stress/input mismatch: '+label)


def fixture():
    c=MPIntervalContext();c.dps=110;P=source.packets
    logs=(c.mpf('-2'),2*c.ln(3),2*c.ln(2),c.mpf('-1'))
    algebra=P.FactoredAlgebra(c,logs,[]);Z=P.IntervalTaylor.variable(c,'.2',5)
    make=lambda value:algebra.lift(value)
    original={k:algebra.shift(1+Z/10,(j%2,-.5,j%2,0)) for j,k in enumerate(P.recovery.RATES)}
    P0=make(-2+Z/3);R=c.mpf('2.3');S=c.mpf(3);count=0
    def coefficient_value(row):
        out=P.IntervalTaylor.constant(c,0,row.order)
        for powers,jet in row.terms.items():
            out=out+jet*c.exp(sum((log*x for log,x in zip(logs,powers)),c.mpf(0)))
        return out
    def polynomial_value(poly):
        out=P.IntervalTaylor.constant(c,0,4)
        for power,row in poly.terms.items():out=out+coefficient_value(row)*R**power
        return out
    def ratio(quotient):return polynomial_value(quotient.numerator)/polynomial_value(quotient.denominator)
    for sign in (1,-1):
        E=[make(1+Z*Z/10),make('.02')+make(Z/100)]+[make('.03')]*3
        V=[make('.7')+make(Z/4),make(str(sign*.03))+make(Z/100)]+[make('.01')]*3
        family=dict(zip(P.FAMILY_KEYS,('fixture','fixture_source','fixture_datum')))
        packet=P.CurrentSourcePacket('fixture',family,algebra,Z,dict(theta=E,axial=V,radial=[make(0)]*5),
            {k:[v]*5 for k,v in original.items()},[P0+original['p']]+[make(0)]*4,P0,{}, {},
            dict(original_radius_source='artificial finite R',cache_cover=True))
        inputs=source.from_packet(packet,'.001');field=packet.recover_original('.001')
        # Independent original raw-pre program, restoring native units.
        erows=[coefficient_value(v) for v in E];vrows=[coefficient_value(v)*S for v in V]
        history={k:[coefficient_value(v)*(S if k in ('m','k') else 1)
                    for v in field['own_normalized_history_ordinary_y_rows'][k]] for k in P.recovery.RATES}
        pressure=[coefficient_value(v) for v in field['physical_velocity_pressure_ordinary_y_rows']['pressure']]
        raw=raw_pre_stress_rows(c,c.mpf('.001'),Z,erows,vrows,history,pressure)
        F=erows[0]*(S/c.sqrt(2*R));native={}
        for component in ('theta','axial'):
            total=erows[0]*0
            for part in raw[component].values():
                rp,sp,_,_=part['mode']
                total=total+part['full_derivative_rows'][0]*(R**c.mpf(str(rp))*S**sp/c.sqrt(2))
            native[component]=total/F
        q=inputs.quotients();a,b,p1,p2=(ratio(q[k]) for k in ('a','b','p1','p2'))
        pairs=((p1-a,native['theta']),(p2+b,native['axial']),
               (ratio(q['a']),1-2*erows[1]/erows[0]),(ratio(q['b']),2*vrows[1]/(erows[0]*S)))
        for left,right in pairs:
            for k in range(min(left.order,right.order)+1):overlap(left[k],right[k],'native physical units');count+=1
        t0=-b/a;kap=a+b*b/a;H=p1+p2*t0;D=H-kap;J=p2-p1*t0
        pairs=((ratio(q['kappa']),kap),(ratio(q['H0_minus2']),H-2),
               (ratio(q['D']),D),(ratio(q['J']),J),
               (polynomial_value(inputs.quadratic_numerator())/polynomial_value(inputs.denominator)**3,
                2*D*D-(kap-2)*J*J))
        for left,right in pairs:
            for k in range(min(left.order,right.order)+1):overlap(left[k],right[k],'division-free signed invariant');count+=1
        if not P.recovery.endpoints(b.terms[0][0] if hasattr(b,'terms') else b[0])[0]*sign>0:
            raise ArithmeticError('Nonzero signed axial fixture was lost')
        if algebra.proofs or algebra.final_rows:raise ArithmeticError('Production source algebra resolved factors')
    return count


def guards():
    P=source.packets;c=MPIntervalContext();c.dps=90
    algebra=P.FactoredAlgebra(c,(c.mpf(0),)*4,[]);one=source.RadiusPolynomial(algebra,{0:1})
    foreign=source.RadiusPolynomial(P.FactoredAlgebra(c,(c.mpf(0),)*4,[]),{0:1})
    calls=[lambda:source.SourceQuotient(one,source.RadiusPolynomial(algebra,{}),'positive'),
        lambda:one+foreign,lambda:source.RadiusPolynomial(algebra,{-1:1}),lambda:one**-1]
    for call in calls:
        try:call()
        except ValueError:pass
        else:raise ArithmeticError('Invalid source quotient or factor basis admitted')
    return len(calls)


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in manifest['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed generic input dependency: '+name)
    records=json.loads(gzip.decompress((source.HERE/source.VIEWS).read_bytes()))
    owner=source.CurrentGenericShearInputs()
    if list(records)!=list(source.packets.CHARTS):raise ValueError('Actual current source chart inventory differs')
    for chart in records:
        inputs=owner.saved(chart);record=source.packets.encode(inputs.record())
        if records[chart]!=record:raise ValueError('Actual full input record changed: '+chart)
        if any(record[k] for k in source.OPEN) or record['whole_box_quotient_bounds_certified']:
            raise ValueError('Unproved whole-current loop input admission')
        if inputs.packet.algebra.proofs or inputs.packet.algebra.final_rows:raise ValueError('Source factor resolution occurred')
    theorem=source.exact_theorem()
    if theorem!=manifest['exact_full_source_input_theorem']:raise ValueError('Full source input identity changed')
    result=dict(all_passed=True,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
        source_family=owner.service.family,current_original_chart_count=len(records),
        exact_full_input_identities=len(theorem['identities']),
        independent_full_native_stress_and_signed_invariant_fixture_coefficients=fixture(),
        invalid_quotient_and_factor_basis_guards=guards(),
        production_source_factors_resolved=False,ancestor_constructors_called=False,
        whole_box_quotient_bounds_certified=False,
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual full generic inputs and division-free invariants PASS:',len(records),'charts',flush=True)
    return result


if __name__=='__main__':run()

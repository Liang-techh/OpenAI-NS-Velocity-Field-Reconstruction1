"""Full reference axial source covered by overlapping original-u predicates.

The local source variable log|u| is defined by u=p2*q/dstar, never selected.
C0/y/Z/yZ templates keep the original pressure remainder and formal scales.
Only the constant-parameter Rh_reference window is covered.
"""
from dataclasses import dataclass
import ast
import copy
import gzip
import json
from pathlib import Path
from types import FunctionType,MappingProxyType,SimpleNamespace
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_near_midplane_mixed_C1_integrals as regular

near=regular.near;whole=regular.whole;points=regular.points;base=regular.base;prior=regular.prior
mixed=regular.mixed;phase=mixed.phase;ep=regular.ep
HERE,PREFIX,sha=regular.HERE,regular.PREFIX,regular.sha
C0,Y,Z,YZ=mixed.ORDERS;MixedJet=mixed.MixedJet
NAME=PREFIX+'current_original_reference_axial_predicate_C1.json.gz'
RECEIPT=PREFIX+'current_original_reference_axial_predicate_C1_check.json'
GATE='original_full_reference_axial_predicate_source_and_C1_integrals_enclosed'


class AxialSourceAtlas(whole.factors.OriginalSourceFactorAtlas):
    """Exact physical Z=Pstar^kp*Cstar^kc*xi; interval Q/L are C0 hulls."""
    def __init__(self,frame,*,lower,upper,carrier=(0,0),logu=None):
        lo,hi=(points.point.source.exact_rational(v) for v in (lower,upper))
        if not lo<hi or carrier not in ((0,0),(-11,-10)):
            raise ValueError('Strict exact interval and supported physical carrier required')
        if carrier==(0,0) and not -1<=lo<hi<=1:
            raise ValueError('Original physical Z domain is [-1,1]')
        super().__init__(frame,Z=0)
        c=self.ctx
        with mp.workdps(c.dps+40):
            self.Z=None;self.bounds=(lo,hi);self.carrier=carrier
            xi=c.mpf((ep(self.rational(lo))[0],ep(self.rational(hi))[1]))
            z2=prior.ScaledEnclosure(prior.FormalScale(self.bases,
                (2*carrier[0],2*carrier[1],0,0,0)),xi**2,self.ledger)
            Q=1+whole.conditioned.bounded(z2)
            L=1-whole.conditioned.bounded(self.parameter('delta')*z2)
            if ep(Q)[0]<1 or ep(L)[0]<=0:raise ArithmeticError('Original Q/L positivity lost')
            self.bases=self.bases[:2]+(c.ln(L),c.mpf(0) if logu is None else c.mpf(ep(logu)),self.bases[4])
            self.Q=Q;self.L=L
            self.logLambda=11*self.bases[0]+10*self.bases[1]
            if ep(self.bases[0])[0]<=0 or ep(self.bases[1])[0]<=0 or ep(self.bases[2])[1]>0:
                raise ValueError('Original positive P/C and nonpositive logL basis required')
            if ep(self.bases[3])[0]<-2 or ep(self.bases[4])!=(0,0):
                raise ValueError('Bounded lower log|u| and exact zero radius basis required')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def zterm(self,powers,coefficient,*,coordinate,xi,Z_power):
        if type(Z_power) is not int or Z_power<0:raise ValueError('Original nonnegative Z power required')
        value=self.term(powers,coefficient*self.ctx.mpf(xi)**Z_power,coordinate=coordinate)
        shift=prior.FormalScale(self.bases,
            (self.carrier[0]*Z_power,self.carrier[1]*Z_power,0,0,0))
        return regular.finite_offset_anchor(self,prior.ScaledEnclosure(
            value.scale+shift,value.coefficient,self.ledger))

    def add(self,left,right):
        left.scale.pair(right.scale)
        if left.ledger is not self.ledger or right.ledger is not self.ledger:
            raise ValueError('Same original axial atlas ledger required')
        if left.zero:return right
        if right.zero:return left
        # Componentwise dominating MONOMIAL, before evaluating any huge log.
        # logP/logC>0, logL<=0; log|u| has a finite lower even though its
        # upper is astronomical. Larger log|u| exponent is the safe anchor.
        lp,rp=left.scale.powers,right.scale.powers
        if lp[4]!=0 or rp[4]!=0:raise ValueError('Axial atlas radius powers must already be collected')
        powers=(max(lp[0],rp[0]),max(lp[1],rp[1]),min(lp[2],rp[2]),
            max(lp[3],rp[3]),max(lp[4],rp[4]))
        anchor=prior.FormalScale(self.bases,powers)
        lc=left.coefficient*left.bounded_exp((left.scale-anchor).evaluate())
        rc=right.coefficient*right.bounded_exp((right.scale-anchor).evaluate())
        self.ledger['axial_componentwise_formal_scale_sums']=self.ledger.get('axial_componentwise_formal_scale_sums',0)+1
        return prior.ScaledEnclosure(anchor,lc+rc,self.ledger)

    def record(self):
        return dict(source_family=self.family,exact_coordinate_bounds=[str(v) for v in self.bounds],
            physical_Z_carrier=list(self.carrier),physical_Z='Pstar^kp*Cstar^kc*xi',
            Q_definition='1+Z^2',L_definition='1-delta*Z^2',Q=self.Q,L=self.L,
            defining_basis=self.bases,Q_L_hulls_not_exact_arithmetic_correlation=True,
            derivatives_from_original_templates_not_hull_derivatives=True,
            physical_Z_not_native_float_or_zero=True)


def source_values(seed,a,left,right):
    c=a.ctx;y=c.mpf((ep(a.rational(left))[0],ep(a.rational(right))[1]))
    f=c.exp(y/10)
    values=(f,c.mpf(5)/8*f,c.mpf(5)/12*f*f,c.mpf(5)/2*f*f,
        a.copy_interval(seed.owner.inputs.alpha_enclosure),a.Q)
    xi=c.mpf((ep(a.rational(a.bounds[0]))[0],ep(a.rational(a.bounds[1]))[1]))
    return y,xi,values


def compiled_row(seed,a,terms,left,right,*,drop_Z=0):
    y,xi,values=source_values(seed,a,left,right);c=a.ctx;out=a.scalar(0)
    for row in terms:
        if row['k']<drop_Z:raise ValueError('Original p2 odd carrier missing')
        coefficient=c.mpf(row['fn'](*values))
        if row['error_order'] is not None:coefficient*=c.exp(c.mpf(3)/5)*c.mpf((-1,1))
        out=a.add(out,a.zterm(row['powers'],coefficient,coordinate=y,xi=xi,Z_power=row['k']-drop_Z))
    return out


def normalize(a,value,powers):
    if value.zero:return a.ctx.mpf(0)
    shift=prior.FormalScale(a.bases,powers)
    try:return value.coefficient*value.bounded_exp((value.scale-shift).evaluate())
    except ArithmeticError as error:
        raise ArithmeticError('Axial normalization powers '+str((value.scale.powers,powers))+
            ' offset '+str(ep(value.scale.offset))+' log range '+str(ep((value.scale-shift).evaluate()))) from error


def hull(c,values):
    return c.mpf((min(ep(v)[0] for v in values),max(ep(v)[1] for v in values)))


def carrier_certificate(seed,*,radial_cells=16,axial_cells=16):
    """Prove p2/Z has one sign on ALL physical Z, not sampled points."""
    frame=seed.owner.inputs.frame;records=[];gb=[];gyb=[]
    z=seed.owner.z;f,H,D,P=seed.owner.templates['inputs'][1:5]
    alpha=seed.owner.alpha;Q=1+z*z
    radial={H:s.Rational(5,8)*f,D:s.Rational(5,12)*f*f,P:s.Rational(5,2)*f*f}
    p0,p1,p2=seed.owner.templates['pressure_symbols']
    baseline={p0:-alpha/Q**2,p1:4*alpha*z/Q**3,p2:alpha*(4-20*z*z)/Q**4}
    errors=(s.Rational(5,2)/Q**2,10*z/Q**3,22/Q**2)
    compiled={}
    # Collect the complete exact carrier numerator before interval evaluation;
    # expanded independent Z powers needlessly lose its sign near |Z|=1.
    for key,rows in ((C0,seed.owner.templates['rows'][('p2',0)]),(Y,seed.templates[('p2',0)])):
        compiled[key]=[]
        for powers,expr in rows:
            pieces=[(None,expr.xreplace(baseline))]
            pieces.extend((order,s.diff(expr,S)*error) for order,(S,error) in
                enumerate(zip(seed.owner.templates['pressure_symbols'],errors,strict=True)))
            for order,value in pieces:
                if value==0:continue
                value=s.cancel(value.subs(radial,simultaneous=True)/z)
                numerator,denominator=s.fraction(value)
                value=s.horner(numerator,f)/s.factor(denominator)
                if s.simplify(value.subs(z,-z)-value)!=0:
                    raise ValueError('Full original carrier must be even in Z')
                actual=powers if order is None else (powers[0],powers[1]-1,powers[2],powers[3])
                compiled[key].append((actual,s.lambdify((f,z,alpha),value,
                    modules=[{'mpf':seed.ctx.mpf},'mpmath']),order))
    for key,terms in ((C0,seed.owner.terms[('p2',0)]),(Y,seed.terms[('p2',0)])):
        if any(row['k']<1 or row['k']%2!=1 for row in terms):
            raise ValueError('Full original p2 and p2_y must have an odd Z carrier')
    for i in range(radial_cells):
        left=-5+s.Rational(5*i,radial_cells);right=-5+s.Rational(5*(i+1),radial_cells)
        for j in range(axial_cells):
            lo=s.Rational(j,axial_cells);hi=s.Rational(j+1,axial_cells)
            a=AxialSourceAtlas(frame,lower=lo,upper=hi)
            with mp.workdps(a.ctx.dps+40):
                y,zi,values=source_values(seed,a,left,right)
                def carrier(order):
                    terms=[]
                    for powers,fn,error in compiled[order]:
                        coefficient=a.ctx.mpf(fn(values[0],zi,values[-2]))
                        if error is not None:coefficient*=a.ctx.exp(a.ctx.mpf(3)/5)*a.ctx.mpf((-1,1))
                        terms.append(regular.finite_offset_anchor(a,a.term(powers,coefficient,coordinate=y)))
                    return a.sum(terms)
                g=carrier(C0);gy=carrier(Y)
                gn=normalize(a,g,(11,10,0,0,0));gyn=normalize(a,gy,(11,10,0,0,0))
                if ep(gn)[1]>=0:raise ArithmeticError('Full original carrier sign needs refinement: '+str((i,j,ep(gn))))
                gb.append(gn);gyb.append(gyn)
                records.append(dict(exact_y_cell=[str(left),str(right)],exact_abs_Z_cell=[str(lo),str(hi)],
                    full_p2_over_Z_divided_by_Lambda0=gn,full_p2_y_over_Z_divided_by_Lambda0=gyn,
                    all_original_pressure_error_terms_retained=True))
    a=AxialSourceAtlas(frame,lower=-1,upper=1);c=a.ctx
    with mp.workdps(c.dps+40):
        g=hull(c,gb);gy=hull(c,gyb);gamma=gy/g
        amin=c.mpf(4)/5
        ar=a.scalar(amin);eta=a.copy_interval(seed.owner.scales.logs['eta'])
        q=base.current.q_enclosure(ar,ar-2,eta,a.copy_interval(seed.owner.scales.logs['a_min']))['q']
        qf=whole.conditioned.bounded(q)
        if ep(qf)[0]<c.mpf('.5'):raise ArithmeticError('Reference weighted curvature q hypothesis lost')
        dlog=a.copy_interval(seed.owner.scales.logs['d_star']);d=c.exp(dlog)
        abs_g=c.mpf((-ep(g)[1],-ep(g)[0]));cmag=abs_g*qf/d
        # Exact integer coordinate radius; choose it from a certified lower,
        # not from a fitted/sample value. Original Z remains formal.
        needed=ep((c.mpf(1)/4)/c.mpf(ep(cmag)[0]))[1]
        zeta_max=s.Integer(int(mp.ceil(needed))+1)
        if ep(a.rational(zeta_max)*c.mpf(ep(cmag)[0]))[0]<c.mpf('.25'):
            raise ArithmeticError('Conditional central coverage failed')
        return dict(passed=True,source_family=seed.family,
            exact_y_window=['-5','0'],exact_physical_Z_window=['-1','1'],
            sign=-1,Lambda0='Pstar^11*Cstar^10',g_normalized=g,g_y_normalized=gy,
            gamma_y=gamma,q=qf,dstar=d,u_over_Lambda0_abs_Z=cmag,
            exact_central_zeta_max=str(zeta_max),whole_source_partition=records,
            all_negative_Z_covered_by_exact_even_g_and_gy_parity=True,
            source_identities=['p2=Z*g','p2_y=Z*g_y','u=Z*g*q/dstar','u_y/u=g_y/g'],
            q_and_dstar_y_Z_derivatives_exact_zero_on_reference=True,
            separate_pressure_error_boxes_not_exact_joint_correlation=True,
            conditional_small_u_implies_abs_Z_at_most_zeta_max_over_Lambda0=True)


def compile_predicate_phase():
    """The sole new phase guard is a proved closed regular-u predicate."""
    source=whole.conditioned.base.conditioned
    tree=ast.parse(Path(source.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='ConditionedPhase')
    fn=copy.deepcopy(next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__'))
    before=copy.deepcopy(fn)
    assign=next(n for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='self.u' for t in n.targets))
    old_u=copy.deepcopy(assign.value);assign.value=ast.parse("query['original_u_source']",mode='eval').body
    guard=next(n for n in ast.walk(fn) if isinstance(n,ast.If) and ast.unparse(n.test).startswith('logupper <='))
    old_guard=copy.deepcopy(guard.test)
    regular_guard=ast.BoolOp(op=ast.Or(),values=[ast.parse("query.get('regular_predicate',False)",mode='eval').body,guard.test])
    guard.test=ast.BoolOp(op=ast.And(),values=[ast.parse("not query.get('signed_predicate',False)",mode='eval').body,regular_guard])
    restored=copy.deepcopy(fn)
    next(n for n in ast.walk(restored) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='self.u' for t in n.targets)).value=old_u
    next(n for n in ast.walk(restored) if isinstance(n,ast.If) and isinstance(n.test,ast.BoolOp)).test=old_guard
    if ast.dump(restored)!=ast.dump(before):raise ValueError('Original phase changed beyond source binding and closed guard')
    env=dict(vars(source),bounded_value=whole.conditioned.bounded)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<original closed conditional regular guard>','exec'),env)
    return env['__init__']


PREDICATE_PHASE_INIT=compile_predicate_phase()
class PredicatePhase(whole.conditioned.PositiveLogQPhase):
    def __init__(self,query,dstar_log):
        if query.get('regular_predicate',False) and query.get('signed_predicate',False):
            raise ValueError('Choose one original coordinate backend on overlap')
        if query.get('regular_predicate',False):
            box=whole.conditioned.bounded(query['original_u_source']);c=query['q'].ctx
            if max(abs(v) for v in ep(box))>ep(c.mpf(1)/4)[1]:
                raise ValueError('Proved closed regular-u source predicate required')
        if query.get('signed_predicate',False):
            value=query['original_u_source'];c=value.ctx;lo,hi=ep(value.coefficient)
            if not lo*hi>0:raise ValueError('Strict signed source predicate required')
            lower=value.scale.evaluate()+c.ln(c.mpf(min(abs(lo),abs(hi))))
            if ep(lower)[0]<ep(c.ln(c.mpf(3)/16))[0]:
                raise ValueError('Original signed |u|>=3/16 predicate required')
        PREDICATE_PHASE_INIT(self,query,dstar_log)



def compile_correlated_first_rows(function,*,signed):
    """Same mixed formulas; enclose first implicit rows before slow sums.

    At fixed true phase, psi_i=-T2_i/(1+t^2) and
    B_i=-a*(E_i*T1+E*(T1_i-t*T2_i/(1+t^2)))/(4*pi).
    The original mixed implementation already computed that total derivative.
    """
    tree=ast.parse(Path(function.__code__.co_filename).read_text(encoding='utf8'))
    fn=copy.deepcopy(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==function.__name__))
    fn.name='correlated_'+function.__name__
    first=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.For) and ast.unparse(n.target)=='order')
    if ast.unparse(fn.body[first].iter)!='(Y, Z)':
        raise ValueError('Original first slow derivative sites changed')
    del fn.body[first]
    inv='R0' if signed else 'invD'
    for node in fn.body:
        if not isinstance(node,ast.Assign):continue
        targets=[ast.unparse(t) for t in node.targets]
        if 'A' in targets:
            fy="roots['a'][C0]*"+inv+"*T2[Y]*(c.mpf(1)/(4*c.pi))"
            fz="roots['a'][C0]*"+inv+"*T2[Z]*(c.mpf(1)/(4*c.pi))"
            if not signed:fy="a.scalar(0) if symmetry else "+fy;fz="a.scalar(0) if symmetry else "+fz
            node.value=ast.parse("MixedJet(a,{C0:primitive['A'],Y:("+fy+"),Z:("+fz+"),YZ:AYZ})",mode='eval').body
        if 'B' in targets:
            fy="Bmixed[Y]";fz="Bmixed[Z]";fyz="Bmixed[YZ]" if signed else "Bmixed_zero"
            if not signed:fy="a.scalar(0) if symmetry else "+fy;fz="a.scalar(0) if symmetry else "+fz
            node.value=ast.parse("MixedJet(a,{C0:primitive['B_over_Pstar'],Y:("+fy+"),Z:("+fz+"),YZ:"+fyz+"})",mode='eval').body
    env=dict(function.__globals__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original mixed function; correlated first implicit rows>','exec'),env)
    return env[fn.name]


REGULAR_MIXED=compile_correlated_first_rows(regular.regular_fixed_and_implicit_mixed,signed=False)
SIGNED_MIXED=compile_correlated_first_rows(mixed.fixed_and_implicit_mixed,signed=True)


@dataclass(frozen=True)
class AxialPredicateFrame:
    owner:object
    family:object
    left:object
    right:object
    branch:str
    roots:object
    query:object
    record:dict


class PredicateBranch(mixed.OriginalReferenceMixedC1):
    def __init__(self,seed,certificate,branch):
        if branch not in ('regular','positive','negative'):raise ValueError('Original-u predicate required')
        self.seed=seed;self.certificate=certificate;self.branch=branch
        frame=seed.owner.inputs.frame
        ordinary=AxialSourceAtlas(frame,lower=-1,upper=1);c=ordinary.ctx
        with mp.workdps(c.dps+40):
            if branch=='regular':
                radius=s.Rational(certificate['exact_central_zeta_max'])
                a=AxialSourceAtlas(frame,lower=-radius,upper=radius,carrier=(-11,-10))
            else:
                upper=ordinary.logLambda+c.ln(c.mpf(ep(certificate['u_over_Lambda0_abs_Z'])[1]))
                logu=c.mpf((ep(c.ln(c.mpf(3)/16))[0],ep(upper)[1]))
                a=AxialSourceAtlas(frame,lower=-1 if branch=='positive' else 0,
                    upper=0 if branch=='positive' else 1,logu=logu)
        self.atlas=a;self.ctx=a.ctx;self.family=seed.family;self.Z=None
        self.parent=SimpleNamespace(owner=seed.owner)
        self.frames={};self.source_cache={};self.primitive_cache={};self.cache={}
        self.hashes=dict(seed.hashes);self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def source_frame(self,left,right):
        left=points.reference.reference_coordinate(left);right=points.reference.reference_coordinate(right)
        if not -5<=left<right<=0:raise ValueError('Original reference source domain required')
        if (left,right) in self.source_cache:return self.source_cache[(left,right)]
        a=self.atlas;c=a.ctx;seed=self.seed
        with mp.workdps(c.dps+40):
            roots={}
            for name in ('E','V','b','p1','p2'):
                roots[name]=MixedJet(a,{
                    C0:compiled_row(seed,a,seed.owner.terms[(name,0)],left,right),
                    Z:compiled_row(seed,a,seed.owner.terms[(name,1)],left,right),
                    Y:compiled_row(seed,a,seed.terms[(name,0)],left,right),
                    YZ:compiled_row(seed,a,seed.terms[(name,1)],left,right)})
            for name,value in (('a',c.mpf(4)/5),('t0',0)):roots[name]=MixedJet.constant(a,a.scalar(value))
            ar=roots['a'][C0];dlog=a.copy_interval(seed.owner.scales.logs['d_star'])
            d=prior.ScaledEnclosure(prior.FormalScale(a.bases,offset=dlog),1,a.ledger)
            q=base.current.q_enclosure(ar,ar-2,a.copy_interval(seed.owner.scales.logs['eta']),
                a.copy_interval(seed.owner.scales.logs['a_min']))['q']
            original_u=(roots['p2'][C0]*q).positive_divide(d,dlog)
            if self.branch=='regular':
                raw=whole.conditioned.bounded(original_u)
                lo=max(ep(raw)[0],-c.mpf(1)/4);hi=min(ep(raw)[1],c.mpf(1)/4)
                if lo>hi:raise ArithmeticError('Regular source predicate has empty enclosure')
                u=a.scalar(c.mpf((ep(lo)[0] if hasattr(lo,'_mpi_') else lo,
                    ep(hi)[1] if hasattr(hi,'_mpi_') else hi)))
            else:
                u=prior.ScaledEnclosure(prior.FormalScale(a.bases,(0,0,0,1,0)),
                    1 if self.branch=='positive' else -1,a.ledger)
            # EXACT source relations, enclosed on the defining branch predicate.
            # No derivative is taken of its interval endpoints.
            p2=d*u
            p2=p2.positive_divide(q,q.scale.evaluate()+c.ln(c.mpf(ep(q.coefficient)[0])))
            gamma=a.scalar(a.copy_interval(self.certificate['gamma_y']))
            original_p2_C0=roots['p2'][C0];original_p2_Y=roots['p2'][Y]
            roots['p2']=MixedJet(a,{C0:p2,Y:p2*gamma,Z:roots['p2'][Z],YZ:roots['p2'][YZ]})
            query=dict(q=q,roots={name:{C0:row[C0],Z:row[Z]} for name,row in roots.items()},
                original_u_source=u,regular_predicate=self.branch=='regular',signed_predicate=self.branch!='regular')
            kernel=PredicatePhase(query,dlog)
            expected='small_r_series' if self.branch=='regular' else 'signed_Mobius'
            if kernel.geometry!=expected:raise ArithmeticError('Closed conditional branch geometry changed')
            record=dict(source_family=self.family,exact_reference_cell=[str(left),str(right)],
                branch=self.branch,predicate='abs(u)<=1/4' if self.branch=='regular' else
                ('u>=3/16' if self.branch=='positive' else 'u<=-3/16'),
                atlas=a.record(),defined_original_u='p2*q/dstar',actual_original_u=u.record(),
                original_unconditioned_p2_C0=original_p2_C0.record(),
                original_unconditioned_p2_y=original_p2_Y.record(),
                same_source_rebased_p2_C0=p2.record(),source_gamma_y=gamma.record(),
                original_p2_Z_yZ_template_rows_retained=True,
                exact_source_relations=['p2=dstar*u/q','p2_y=p2*(g_y/g)'],
                source_predicate_not_field_or_derivative_of_cap=True,
                conditional_domain_not_claimed_entire_physical_rectangle=True,
                actual_root_jets={name:row.record() for name,row in roots.items()},
                kernel_geometry=kernel.geometry_record())
            result=AxialPredicateFrame(self,self.family,left,right,self.branch,MappingProxyType(roots),
                dict(kernel=kernel),record)
            self.frames[id(result)]=result;self.source_cache[(left,right)]=result;return result

    def primitive(self,frame,coordinate):
        if type(frame) is not AxialPredicateFrame or frame.owner is not self or self.frames.get(id(frame)) is not frame:
            raise ValueError('Issued original conditional source frame required')
        with mp.workdps(self.ctx.dps+40):
            x=self.ctx.mpf(coordinate);key=(id(frame),ep(x))
            if key not in self.primitive_cache:
                function=REGULAR_MIXED if self.branch=='regular' else SIGNED_MIXED
                self.primitive_cache[key]=function(self.atlas,frame.query['kernel'],frame.roots,x)
            return self.primitive_cache[key]


class OriginalReferenceAxialPredicates(mixed.OriginalReferenceMixedC1):
    def __init__(self):
        self.seed=regular.OriginalReferenceNearMidplaneMixed()
        receipt=json.loads((HERE/regular.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(regular.GATE):
            raise ValueError('Accepted original regular mixed source required')
        self.hashes=dict(self.seed.hashes)
        for name,digest in {**receipt['input_hashes'],regular.RECEIPT:sha(regular.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original axial predicate dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Source families differ')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.certificate=carrier_certificate(self.seed)
        self.branches={name:PredicateBranch(self.seed,self.certificate,name) for name in ('regular','positive','negative')}
        self.atlas=AxialSourceAtlas(self.seed.owner.inputs.frame,lower=-1,upper=1)
        self.ctx=self.atlas.ctx;self.family=self.seed.family;self.Z=None
        self.cache={};self.parent=SimpleNamespace(owner=self.seed.owner)

    def cell(self,left,right,*,N):
        left=points.reference.reference_coordinate(left);right=points.reference.reference_coordinate(right)
        key=(left,right,N)
        if key in self.cache:return self.cache[key]
        a=self.atlas;c=a.ctx
        with mp.workdps(c.dps+40):
            pieces={name:owner.cell(left,right,N=N) for name,owner in self.branches.items()}
            caps={jet:{} for jet in (C0,Z)}
            for jet in caps:
                powers=(0,0,0,0,0) if jet==C0 else (11,10,0,0,0)
                for name in points.exact.RATES:
                    caps[jet][name]={}
                    for kind in ('leading','slow_y','nonlinear_remainder','direct_density'):
                        boxes=[a.copy_interval(normalize(owner.atlas,pieces[b]['caps'][jet][name][kind],powers))
                            for b,owner in self.branches.items()]
                        if any(ep(v)[0]<0 for v in boxes):raise ArithmeticError('Positive branch magnitude required')
                        bound=max(ep(v)[1] for v in boxes)
                        caps[jet][name][kind]=prior.ScaledEnclosure(prior.FormalScale(a.bases,powers),c.mpf((0,bound)),a.ledger)
            record=dict(exact_reference_cell=[str(left),str(right)],candidate_N=N,
                original_conditional_source_cells={name:part['record'] for name,part in pieces.items()},
                same_original_function_branch_caps_union_not_sum=True,
                entire_physical_Z_domain_covered_by_three_original_u_predicates=True,
                C0_Z_averaging_caps={str(jet):{name:{kind:v.record() for kind,v in row.items()}
                    for name,row in rows.items()} for jet,rows in caps.items()})
            result=dict(caps=caps,record=record);self.cache[key]=result;return result

    def integrate(self,*,count,N):
        with mp.workdps(self.ctx.dps+40):
            result=super().integrate(count=count,N=N);result.pop('original_Z_exact')
            result.update(exact_physical_Z_window=['-1','1'],
                actual_mixed_source_and_C1_averaging_installed_on_this_fixed_Z_reference_window=False,
                actual_mixed_source_and_C1_averaging_installed_on_whole_reference_Z_window=True,
                no_cross_chart_seam_assumed=True,
                original_u_coordinate_overlap_same_source_identity_required=True,
                branch_predicates=['abs(u)<=1/4','u>=3/16','u<=-3/16'],
                overlapping_predicate_caps_unioned_not_summed=True,
                ordinary_Z_normalization_only='Lambda0=Pstar^11*Cstar^10',
                q_nu_constant_reference_only=True,O2_nonconstant_parameter_extension_installed=False)
            return result


def run():
    begin=time.monotonic();owner=OriginalReferenceAxialPredicates();levels=[]
    print('Full original p2/Z carrier sign proved on256 source rectangles',flush=True)
    for count,N in ((4,160),(16,16384)):
        levels.append(owner.integrate(count=count,N=N))
        print('Original full reference axial predicate C1 integration:',count,'cells N',N,flush=True)
    report=dict(**{GATE:True},source_family=owner.family,carrier_certificate=owner.certificate,
        actual_full_reference_axial_C1_levels=levels,exact_reference_domain=dict(y=['-5','0'],Z=['-1','1']),
        original_four_source_jets_and_full_pressure_errors_retained=True,
        same_original_function_regular_signed_overlap=True,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(near.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Full physical Z[-1,1] on constant-parameter Rh_reference[-5,0], conditional original-u regular/signed charts, source p2/Z sign, correlated p2_y, original Z/yZ roots, mixed primitives and actual-N own-rate C0/Z contribution bounds. Not all-route/O2, terminal closure, global N, stress, recursion or corrected NS.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    return report


if __name__=='__main__':run()

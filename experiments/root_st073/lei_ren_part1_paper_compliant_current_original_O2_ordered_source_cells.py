"""Original ordered O2 masses and whole-cell C0/ordinary-Z source rows.

Accepted J cells feed all three defining mass integrals. Prefix/suffix
integrals bind continuous f/H/D/P and pressure cancellation before interval
arithmetic. Source factor powers and positive late-pressure errors remain.
These are function ranges, not point fields or matched five controls.
"""
import ast
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_five_own_integrals as five

base=five.base;ep=five.ep;HERE,PREFIX,sha=five.HERE,five.PREFIX,five.sha
NAME=PREFIX+'current_original_O2_ordered_source_cells.json.gz'
RECEIPT=PREFIX+'current_original_O2_ordered_source_cells_check.json'
GATE='original_O2_ordered_defining_masses_and_whole_cell_C0_Z_coefficient_ranges_connected'
MASS_SPECS=((sy.Rational(8,5),1),(sy.Rational(1,5),2),(sy.Rational(6,5),2))


def interval(c,value):return five.interval(c,value)


def hull(c,lower,upper):return c.mpf((ep(lower)[0],ep(upper)[1]))


def rational(c,value):return c.mpf(int(value.p))/int(value.q)


def source_alpha_binding():
    name=PREFIX+'current_original_pressure_point_jets.py'
    tree=ast.parse((HERE/name).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='OriginalNormalizedPressurePointJets')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    rows={ast.unparse(t):n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) for t in n.targets}
    expected=dict(mass="c.quad(lambda v:c.exp(v/5-c.mpf('1.2')*self.frame.owner.profiles.J(v)),[0,c.mpf('.5'),1])",
        **{'self.alpha':"c.mpf(5)/2+mass/2+c.exp(-c.mpf(2)/5)/2"})
    for key,value in expected.items():
        if ast.dump(rows[key])!=ast.dump(ast.parse(value,mode='eval').body):raise ValueError('Original pressure alpha source definition changed')
    return dict(passed=True,original_alpha_defining_assignments_bound=True,source=name,sha256=sha(name))


def source_basis_binding():
    name=Path(base.__file__).name;tree=ast.parse((HERE/name).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='OriginalO2ConditionedPrimitives')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='query')
    node=next(n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='bases' for t in n.targets))
    if ast.dump(node)!=ast.dump(ast.parse('(logP,logdelta,c.ln(L),c.mpf(0),logR)',mode='eval').body):
        raise ValueError('Original conditioned consumer semantic scale basis changed')
    return dict(passed=True,source=name,sha256=sha(name),
        source_factor_order=['R','Pstar','delta','L'],consumer_log_basis_order=['logPstar','logdelta','logL','zero','logR'],
        power_mapping='(r,p,d,ell) -> (p,d,ell,0,r)',third_slot_is_logL_not_zero=True)


def pressure_correlation(templates):
    """Substitute exact P0 baseline using P(y)+remaining-pressure-mass(y)."""
    z,f,H,D,P,p0,p0Z,p0ZZ=templates['inputs'];W=sy.Symbol('remaining_pressure_mass',real=True)
    q=1+z*z;alpha=sy.Symbol('same_alpha',real=True)
    pressure=(-alpha/q**2,4*alpha*z/q**3,(4-20*z*z)*alpha/q**4)
    baseline={symbol:value.subs(alpha,P+W) for symbol,value in zip((p0,p0Z,p0ZZ),pressure)}
    rows={};identities={}
    for key,terms in templates['rows'].items():
        rebuilt=[]
        for index,(powers,expr) in enumerate(terms):
            correlated=sy.cancel(expr.subs(baseline))
            # W=alpha-P is an exact function relation at fixed y, not a
            # subtraction of independently selected numeric coefficients.
            assert sy.cancel(correlated.subs(W,alpha-P)-expr.subs(dict(zip((p0,p0Z,p0ZZ),pressure))))==0
            sensitivities=tuple(sy.cancel(sy.diff(expr,v).subs(baseline)) for v in (p0,p0Z,p0ZZ))
            rebuilt.append((powers,correlated,sensitivities))
            identities[str(key)+'_'+str(index)]=True
        rows[key]=tuple(rebuilt)
    return rows,(z,f,H,D,P,W),dict(passed=True,exact_original_coefficient_substitution_identities=identities,
        original_alpha_equals_current_P_plus_remaining_pressure_mass=True,
        independent_large_P_P0_interval_subtraction_avoided=True,
        original_late_pressure_remainder_not_zeroed=True)


class OriginalO2OrderedSourceCells:
    def __init__(self,mass_dps=90):
        if type(mass_dps) is not int or mass_dps<90:raise ValueError('At least 90 directed defining-mass digits required')
        accepted=json.loads((HERE/five.RECEIPT).read_bytes())
        if not accepted.get('all_passed') or not accepted.get(five.GATE):raise ValueError('Accepted original O2 five integral source required')
        self.parent=five.OriginalO2FiveOwnIntegrals();self.owner=self.parent.owner;self.family=self.parent.family;self.c=self.owner.ctx
        self.mass_c=base.MPIntervalContext();self.mass_c.dps=mass_dps
        self.precisions=dict(directed_defining_mass_digits=mass_dps,native_source_factor_digits=self.c.dps,
            separate_directed_contexts_with_ranges_rounded_outward=True)
        if accepted['source_family']!=self.family:raise ValueError('Ordered mass source family differs')
        self.hashes=dict(self.parent.hashes)
        for name,digest in {**accepted['input_hashes'],five.RECEIPT:sha(five.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Ordered mass source changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Ordered source hashes disagree')
            self.hashes[name]=digest
        self.alpha_binding=source_alpha_binding();self.basis_binding=source_basis_binding()
        self.rows,self.symbols,self.correlation=pressure_correlation(self.owner.inputs.templates)
        self.compiled={};c=self.c
        for key,rows in self.rows.items():
            self.compiled[key]=tuple((powers,
                sy.lambdify(self.symbols,expr,modules=[{'mpf':c.mpf},'mpmath']),
                tuple(sy.lambdify(self.symbols,v,modules=[{'mpf':c.mpf},'mpmath']) for v in sensitivities))
                for powers,expr,sensitivities in rows)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.levels={}

    def integrate(self,saved):
        c=self.mass_c;count=saved['ordered_source_cells'];began=time.monotonic()
        if saved['source_family']!=self.family or saved['exact_y_window']!=['0','1']:raise ValueError('Same original source window required')
        with mp.workdps(c.dps+40):
            width=c.mpf(1)/count;contributions=[];weights=[]
            prefix=[[c.mpf(0)]*3];suffix=[[c.mpf(0)]*3 for _ in range(count+1)]
            for i,row in enumerate(saved['whole_source_cells']):
                left,right=c.mpf(i)/count,c.mpf(i+1)/count;J=interval(c,row['J_source_cell'])
                if row['exact_y_cell']!=[str(i)+'/'+str(count),str(i+1)+'/'+str(count)]:raise ValueError('Ordered J source gap')
                masses=[];cell_weights=[]
                for rate,power in MASS_SPECS:
                    r=rational(c,rate)
                    # Here dy>=1/8192 in shipped levels, so a 90-digit
                    # directed subtraction is well conditioned and positive.
                    weight=(c.exp(r*right)-c.exp(r*left))/r
                    if ep(weight)[0]<=0:raise ArithmeticError('Positive finite original exponential cell mass lost')
                    masses.append(weight*c.exp(-c.mpf(3)*power*J/5));cell_weights.append(weight)
                contributions.append(masses);weights.append(cell_weights)
                prefix.append([a+b for a,b in zip(prefix[-1],masses)])
            for i in reversed(range(count)):suffix[i]=[a+b for a,b in zip(suffix[i+1],contributions[i])]
            alpha=c.mpf(5)/2+prefix[-1][1]/2+c.exp(-c.mpf(2)/5)/2
            accepted_alpha=self.owner.inputs.alpha_enclosure
            old=hull(c,accepted_alpha,accepted_alpha)
            assert max(ep(alpha)[0],ep(old)[0])<=min(ep(alpha)[1],ep(old)[1])
            # This is an independent source approximation, never a field
            # cover midpoint or a choice of native pressure parameters.
            approximate=self.owner.inputs.alpha;assert ep(alpha)[0]<=approximate<=ep(alpha)[1]
            records=[];nodes=[]
            for i in range(count+1):
                y=c.mpf(i)/count
                H=(c.mpf(5)/8+prefix[i][0])*c.exp(-3*y/2)
                D=(c.mpf(5)/12+prefix[i][2]/2)*c.exp(-y)
                P=c.mpf(5)/2+prefix[i][1]/2
                W=suffix[i][1]/2+c.exp(-c.mpf(2)/5)/2
                nodes.append(dict(exact_y=str(i)+'/'+str(count),ordered_original_mass_prefixes=prefix[i],
                    ordered_original_mass_suffixes=suffix[i],H=H,D=D,P=P,remaining_pressure_mass=W))
                if i==count:continue
                source=saved['whole_source_cells'][i];ys=hull(c,y,c.mpf(i+1)/count)
                mass_cells=[hull(c,prefix[i][j],prefix[i+1][j]) for j in range(3)]
                profiles=dict(f=interval(c,source['f_source_cell']),
                    H=(c.mpf(5)/8+mass_cells[0])*c.exp(-3*ys/2),
                    D=(c.mpf(5)/12+mass_cells[2]/2)*c.exp(-ys),P=c.mpf(5)/2+mass_cells[1]/2,
                    a=interval(c,source['a_source_cell']),
                    remaining_pressure_mass=hull(c,suffix[i+1][1],suffix[i][1])/2+c.exp(-c.mpf(2)/5)/2)
                assert all(ep(profiles[k])[0]>0 for k in ('f','H','D','P','a','remaining_pressure_mass'))
                records.append(dict(original_source_cache_cell_index=i,exact_y_cell=source['exact_y_cell'],
                    original_J_source_cover=source['J_source_cell'],
                    exact_positive_original_mass_weights=weights[i],directed_original_mass_cell_contributions=contributions[i],
                    whole_original_radial_profile_covers=profiles,whole_cell_mass_covers=mass_cells,
                    source_ranges_not_selected_as_field_points=True,
                    nonmonotone_H_D_f_enclosed_on_entire_y_mass_cell=True))
        result=dict(source_family=self.family,exact_y_window=['0','1'],ordered_source_cells=count,
            accepted_whole_source_cache=dict(filename=self.parent.compact,source_cell_level=count),
            original_mass_specs=[dict(rate=str(rate),power=power) for rate,power in MASS_SPECS],
            ordered_nodes=nodes,whole_source_cells=records,original_alpha_directed_source_enclosure=alpha,
            accepted_independent_alpha_enclosure_overlap=True,original_alpha_binding=self.alpha_binding,
            original_pressure_cancellation_identity=self.correlation,
            original_normalized_midplane_reference_histories_at_y1=dict(m=c.mpf(0),h=nodes[-1]['H'],k=c.mpf(0),e=-nodes[-1]['D'],p=nodes[-1]['P']),
            no_nested_original_quadrature_or_ancestor_producer_reexecuted=True,effective_numerical_precisions=self.precisions,
            source_cell_enclosures_not_original_field_point_values=True,execution_seconds=time.monotonic()-began)
        self.levels[count]=(result,saved);return result

    def coefficient_cell(self,count,index,*,Z_lower,Z_upper):
        if count not in self.levels or type(index) is not int or not 0<=index<count:raise ValueError('Explicit admitted ordered source cell required')
        level,saved=self.levels[count];cell=level['whole_source_cells'][index];c=self.c
        lower=base.point.pressure.exact_Z(Z_lower);upper=base.point.pressure.exact_Z(Z_upper)
        if lower>upper:raise ValueError('Ordered real Z range required')
        with mp.workdps(c.dps+40):
            z=hull(c,rational(c,lower),rational(c,upper));q=1+z**2
            profiles={k:c.mpf(ep(v)) for k,v in cell['whole_original_radial_profile_covers'].items()}
            values=(z,*(profiles[k] for k in ('f','H','D','P','remaining_pressure_mass')))
            rows={};logP=c.exp(40)+11
            tail_logs=[ep(c.mpf(3)/5-logP+c.ln(t)-c.ln(2)-2*c.ln(q))[1] for t in (5,10,44)]
            exact_midplane=lower==0 and upper==0
            if exact_midplane:tail_logs[1]=None
            for key,terms in self.compiled.items():
                row=[]
                for powers,fn,sensitivities in terms:
                    cover=c.mpf(fn(*values));logs=[]
                    for sensitivity,tail in zip(sensitivities,tail_logs):
                        bound=ep(abs(c.mpf(sensitivity(*values))))[1]
                        if bound and tail is not None:logs.append(ep(c.ln(c.mpf(bound))+c.mpf(tail))[1])
                    late=None if not logs else ep(c.mpf(max(logs))+c.ln(len(logs)))[1]
                    if ep(cover)!=(0,0) or late is not None:
                        row.append(dict(original_R_Pstar_delta_L_powers=list(powers),finite_coefficient_function_cover=cover,
                            late_pressure_coefficient_error_log_upper=late,
                            no_field_point_or_coefficient_midpoint_selected=True))
                rows[key]=row
            rows[('a',0)]=[dict(original_R_Pstar_delta_L_powers=[0,0,0,0],
                finite_coefficient_function_cover=profiles['a'],late_pressure_coefficient_error_log_upper=None,
                no_field_point_or_coefficient_midpoint_selected=True)]
            rows[('a',1)]=[]
            ys=hull(c,c.mpf(index)/count,c.mpf(index+1)/count)
            logdelta=-4*logP-30;delta=c.exp(logdelta);L=1-delta*z**2
            assert ep(L)[0]>0
            logC=c.mpf(mp.mp.make_mpf(self.owner.inputs.frame.selected_logCstar_mpf_tuple))
            bases=(logP,logdelta,c.ln(L),c.mpf(0),c.ln(110)+10*(logC+logP)+ys)
            ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
                positive_function_root_intersections=0,directed_independent_log_rescalings=0)
            scalar=lambda value:base.prior.ScaledEnclosure(base.prior.FormalScale(bases),value,ledger)
            roots={}
            for key,row in rows.items():
                result=scalar(0)
                for term in row:
                    r,p,d,ell=term['original_R_Pstar_delta_L_powers'];cover=term['finite_coefficient_function_cover']
                    if term['late_pressure_coefficient_error_log_upper'] is not None:
                        error=ep(scalar(1).bounded_exp(c.mpf(term['late_pressure_coefficient_error_log_upper'])))[1]
                        cover+=c.mpf((-error,error))
                    result+=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,(p,d,ell,0,r)),cover,ledger)
                roots.setdefault(key[0],{})[(0,key[1])]=result
            roots['t0']={(0,0):scalar(0),(0,1):scalar(0)}
            qrecord=saved['whole_source_cells'][index]['q_source_cell'];qscale=qrecord['formal_positive_scale']
            if any(qscale['source_exponents']) or qscale['radius_power']:
                raise ValueError('Z-independent original q range acquired L/radius powers')
            qsource=five.restored_range(qrecord,bases,ledger)
            if qsource.zero:raise ValueError('Original positive eta q source cannot become an exact flat branch')
            kernel=base.conditioned.ConditionedPhase(dict(q=qsource,roots=roots),self.owner.scales.logs['d_star'])
            if exact_midplane:
                assert roots['p2'][(0,0)].zero and roots['V'][(0,0)].zero
                assert not roots['p2'][(0,1)].zero and not roots['V'][(0,1)].zero
        source=dict(source_family=self.family,source_cell_level=count,source_cell_index=index,
            exact_y_cell=cell['exact_y_cell'],exact_Z_range=[str(lower),str(upper)],
            full_original_factored_coefficient_function_ranges={k:[rows[(k,0)],rows[(k,1)]] for k in ('E','V','a','b','p1','p2')},
            native_source_root_enclosures={k:{'y0_Z%d'%order:value.record() for (y,order),value in row.items()} for k,row in roots.items()},
            original_basis_contract=dict(order=['logPstar','logdelta','logL','zero','logR'],semantic_consumer_binding=self.basis_binding,
                original_L='1-original_delta*Z²',original_R='Rref*exp(y)',positive_L_interval=L,
                delta_and_R_Z_independent=True,physical_factors_not_materialized=True),
            pressure_prefix_suffix_cancellation_proof=self.correlation,original_P0_datum_sha256=self.family['datum_enclosure_sha256'],
            ordinary_Z_values_not_Taylor_coefficients=True,exact_Z0_parity_and_nonzero_p2_Z_preserved=exact_midplane,
            whole_source_ranges_not_field_points=True,conditioned_geometry=kernel.geometry_record(),
            inverse_or_modulated_nonmidplane_integral_installed=False,
            numerical_original_source_point_or_integral_oracle_installed=False,
            actual_five_controls_installed=False,current_whole_N_selected=False)
        return dict(record=source,roots=roots,q=qsource,kernel=kernel,ledger=ledger)


def run():
    began=time.monotonic();owner=OriginalO2OrderedSourceCells();levels=[]
    for saved in owner.parent.saved['actual_original_pressure_integral_refinements']:
        levels.append(owner.integrate(saved));print('Original ordered O2 masses:',saved['ordered_source_cells'],flush=True)
    count=levels[-1]['ordered_source_cells'];records=[]
    for y,zl,zh in (('.23','.37','.37'),('.53','-.37','-.37'),('.53','0','0'),('.53','.36','.38'),
                    ('.53','-.38','-.36'),('.53','-.01','.01'),('.53','-1','1'),('.97','.7','.7'),('1','.37','.37')):
        index=min(int(mp.mpf(y)*count),count-1)
        records.append(owner.coefficient_cell(count,index,Z_lower=zl,Z_upper=zh)['record'])
    result=dict(**{GATE:True},source_family=owner.family,
        actual_original_ordered_source_levels=levels,whole_cell_C0_and_ordinary_Z_coefficient_queries=records,
        effective_numerical_precisions=owner.precisions,
        original_ordered_J_and_three_mass_profile_cells_installed=True,
        exact_remaining_pressure_mass_cancellation_installed=True,
        all_physical_source_factors_and_positive_late_pressure_errors_preserved=True,
        numerical_original_source_point_or_integral_oracle_installed=False,
        inverse_or_modulated_nonmidplane_integral_installed=False,
        actual_changed_five_moment_integral_evaluated=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(base.point.source.inertial.profiles.loop.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Original O2 ordered defining masses, continuous y-cell f/H/D/P and pressure suffix, full fixed/nonzero/Z-range C0 and ordinary-Z coefficient covers with exact source factors and original pressure error. Not a point field, nonmidplane inverse/integral, terminal functional matching, all-chart oracle, global N, controls, stress or recursion.')
    content=json.dumps(base.encoded(result),indent=2).encode('utf8')+b'\n'
    (HERE/NAME).write_bytes(gzip.compress(content,compresslevel=9,mtime=0))
    return result


if __name__=='__main__':run()

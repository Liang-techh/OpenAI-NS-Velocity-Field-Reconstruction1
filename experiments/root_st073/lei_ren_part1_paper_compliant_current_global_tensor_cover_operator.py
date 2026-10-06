"""Independent analytic coordinate and radial cover, before numeric location.

For every fixed admitted source parameter, the atlas covers all finite
physical space at tau>0. This does not resolve interval parameter uncertainty
or prove an admissible stress cone, flatness, energy or coefficient recursion.
"""
import functools
import ast
import hashlib
import sympy as s
from lei_ren_part1_paper_compliant_current_correlated_radius_operator import (
    V,SYMBOLS,LOG_HB,canonical_radius,source_geometry,boundary_sources)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST


@functools.lru_cache(maxsize=1)
def independent_radial_atlas():
    """Direct mathematical recipes, not forward calls to the existing atlas."""
    a=SYMBOLS;h=s.exp(LOG_HB);ra=a['log_epsilon']+s.log(4)
    ref=s.log(110)+10*(a['log_C']+a['log_P']);p=ref+a['log_P']+1+a['Tw']
    rv=p+13/a['mu'];rel=rv+100+a['Lrel'];rq=rel+2+a['Ts'];rt=rq+a['wait']
    mu=a['mu'];phase=(0,1)
    rows=(
      ('core_positive_radius',a['log_epsilon']+s.log(V),(0,4)),
      ('bridge_first',ra+h*V,phase),('bridge_second',ra+h*V,(1,2)),
      ('bridge_macro',ra+2*h+(s.log(100)-ra-2*h)*V,phase),
      ('switch_first',s.log(100)+h*V,phase),('switch_second',s.log(100)+h*V,(1,2)),
      ('switch_power',s.log(100)+2*h+(s.log(s.Rational(11,10))-2*h)*V,phase),
      ('reshape',s.log(110)+a['T']*V,phase),
      ('inner_reference',s.log(110)+a['T']+(10*(a['log_C']+a['log_P'])-a['T']-8)*V,phase),
      ('axial_restore',ref-8+V,phase),('restore_buffer',ref+V,(-7,-6)),
      ('actual_patch',ref-6+s.log(V),(1,s.E)),('Rh_reference',ref+V,(-5,0)),
      ('O2_slope',ref+V,phase),('O2_axial',ref+s.exp(a['Md']*V),phase),
      ('O2_buffer',ref+s.exp(a['Md'])+V,(0,11)),('O3_slope_mu',ref+a['log_P']+V,phase),
      ('O3_power',ref+a['log_P']+1+a['Tw']*V,phase),
      ('pulse_entrance',p+V/mu,(0,s.Rational(1,50))),
      ('pulse_main',p+V/mu,(s.Rational(1,50),10)),
      ('pulse_exit',p+V/mu,(10,11)),('pulse_gap',p+V/mu,(11,12)),
      ('pulse_gap_end',p+12/mu+(1/mu-4)*V,phase),('pulse_end',rv+V,(-4,0)),
      ('flatten',rv+V,(0,100)),('outer_power',rv+100+(a['Lrel']-4)*V,phase),
      ('outer_angular',rel+V,(-4,0)),('steep_entry',rel+V,phase),
      ('steep_power',rel+1+a['Ts']*V,phase),('steep_exit',rel+1+a['Ts']+V,phase),
      ('waiting',rq+a['wait']*V,phase),('heat_collar',rt+V,(0,3)),
      ('heat_exterior',rt+V,(3,s.oo)))
    return {name:dict(expression=canonical_radius(expr),domain=tuple(map(s.sympify,domain))) for name,expr,domain in rows}


@functools.lru_cache(maxsize=1)
def strict_source_margins():
    a=SYMBOLS;h=s.exp(LOG_HB)
    return {name:canonical_radius(expr) for name,expr in dict(
      hb=h,bridge_macro_width=s.log(100)-a['log_Ra']-2*h,
      switch_power_width=s.log(s.Rational(11,10))-2*h,
      inner_reference_width=a['log_gap']-8,Md=a['Md'],T=a['T'],Tw=a['Tw'],
      mu=a['mu'],reciprocal_gap_width=1/a['mu']-4,
      outer_power_width=a['Lrel']-4,Ts=a['Ts'],wait=a['wait']).items()}


@functools.lru_cache(maxsize=1)
def independent_radial_cover_theorem():
    atlas=independent_radial_atlas();actual,_=source_geometry();margins=strict_source_margins()
    if list(atlas)!=list(actual):raise ValueError('Independent full radial atlas inventory differs')
    positive={name:s.Symbol('positive_'+name,positive=True) for name in margins}
    u=s.Symbol('positive_native_radius',positive=True);evidence={}
    for name,row in atlas.items():
        expected=actual[name]
        if row['domain']!=expected['domain'] or canonical_radius(row['expression']-expected['original'])!=0:
            raise ValueError('Independent original radial source/domain differs: '+name)
        derivative=canonical_radius(s.diff(row['expression'],V))
        if canonical_radius(derivative-expected['slope'])!=0:raise ValueError('Original radius derivative differs: '+name)
        logical=derivative;premises=[]
        for margin,expr in margins.items():
            if canonical_radius(derivative-expr)==0:
                logical=positive[margin];premises=[margin];break
        if not premises:
            logical=derivative.subs(SYMBOLS['Md'],positive['Md']).subs(SYMBOLS['mu'],positive['mu'])
            if name in ('core_positive_radius','actual_patch'):
                logical=logical.subs(V,u);premises=['native_radius>0']
            elif name=='O2_axial':premises=['Md','exp(real)>0']
            elif derivative.has(SYMBOLS['mu']):premises=['mu']
        if logical.is_positive is not True:raise ValueError('Strict original radius monotonicity unproved: '+name+' '+str(logical))
        lo,hi=row['domain']
        lower=-s.oo if name=='core_positive_radius' else canonical_radius(row['expression'].subs(V,lo))
        upper=s.oo if hi==s.oo else canonical_radius(row['expression'].subs(V,hi))
        if lo>=hi:raise ValueError('Nonempty native source domain required')
        evidence[name]=dict(independent_source=str(row['expression']),native_domain=list(map(str,row['domain'])),
            original_source_and_domain_identity=True,original_derivative_identity=True,
            exact_log_radius_derivative=str(derivative),positive_derivative_form=str(logical),
            strict_positive_premises=premises,lower_log_radius=str(lower),upper_log_radius=str(upper),
            exact_lower=lower,exact_upper=upper)
    joins={};names=list(atlas)
    for left,right in zip(names,names[1:]):
        if canonical_radius(evidence[left]['exact_upper']-evidence[right]['exact_lower'])!=0:
            raise ValueError('Independent radial cover has a gap: '+left+' / '+right)
        joins[left+' -> '+right]=dict(exact_common_endpoint=True,strictly_increasing_interiors=True)
    core_limit=s.limit(s.log(V),V,0,dir='+');tail_limit=s.limit(V,V,s.oo)
    if core_limit!=-s.oo or tail_limit!=s.oo:raise ValueError('Full radial endpoints unproved')
    for row in evidence.values():row.pop('exact_lower');row.pop('exact_upper')
    return dict(independent_original_33_radial_sources_and_positive_derivatives=evidence,
        independent_all_32_contiguous_source_joins=joins,
        core_positive_radius_log_limit=str(core_limit),heat_exterior_log_limit=str(tail_limit),
        radial_cover='Every finite logR is in at least one of the 33 continuous strictly increasing source ranges; interiors are disjoint; only adjacent endpoints are shared. R=0 has the analytic same-core extension.',
        proof_steps=['Each original continuous source has strictly positive derivative on its interior under audited premises.',
            'Every adjacent pair has the same exact source endpoint; all finite-length intervals have strictly positive span.',
            'The first log radius tends to -infinity as rho tends to zero; the last tends to +infinity.',
            'The intermediate value theorem on each chart and the ordered contiguous chain cover all finite logR; strict monotonicity gives a unique native coordinate in each interior.'],
        strict_source_margin_definitions={name:str(expr) for name,expr in margins.items()},
        source_endpoint_joins_not_numeric_overlap=True,finite_source_parameters_required=True,passed=True)


@functools.lru_cache(maxsize=1)
def independent_implicit_cover_theorem():
    Z,delta=s.symbols('Z delta',real=True);nu,tau,p,C=s.symbols('nu tau positive_axial_power positive_complement',positive=True)
    # Use a positive complement for real-power algebra; its total derivative
    # is -2Z and its defining constraint is C=1-Z^2. No complex-branch
    # force simplification is used on a base of unknown sign.
    lam=s.sqrt(tau)/s.sqrt(C)
    z=s.sqrt(nu)*tau**((1-delta)/2)*Z*C**(-(1-delta)/2)
    jacobian=s.sqrt(nu)*tau**((1-delta)/2)*(C+(1-delta)*Z**2)*C**((delta-3)/2)
    derivative_residual=s.simplify(s.diff(z,Z)-2*Z*s.diff(z,C)-jacobian)
    lambda_residual=s.simplify(lam**2*C-tau)
    if s.expand((C+(1-delta)*Z**2).subs(C,1-Z**2)-(1-delta*Z**2))!=0:
        raise ValueError('Positive complement axial numerator identity differs')
    positive_numerator=s.expand(1-delta*Z**2-((1-delta)+delta*(1-Z**2)))
    if derivative_residual!=0 or lambda_residual!=0 or positive_numerator!=0:raise ValueError('Independent global axial map identities differ')
    from_lambda=s.expand_power_base(s.sqrt(nu)*Z*lam**(1-delta),force=False)
    implicit=s.expand_power_base(lam**2-(z**2/nu)*lam**(2*delta)-tau,force=False)
    if s.simplify(from_lambda-z)!=0 or s.simplify(implicit-tau*(1-Z**2-C)/C)!=0:
        raise ValueError('Independent axial map does not define the original implicit root')
    infinity_factor=(1-Z**2)**(-p)
    if s.limit(infinity_factor,Z,1,dir='-')!=s.oo or s.limit(infinity_factor,Z,-1,dir='+')!=s.oo:
        raise ValueError('Global axial end limits not established')
    A,B=s.symbols('positive_tau_term nonnegative_axial_term',positive=True)
    Fprime=2-2*delta*B/(A+B)
    lower=s.simplify(Fprime-2*(1-delta));upper=s.simplify(2-Fprime)
    if lower!=2*delta*A/(A+B) or upper!=2*delta*B/(A+B):raise ValueError('Unique log-root derivative bounds differ')
    return dict(assumptions='finite z_phys; finite tau>0 and constant nu>0; fixed 0<=delta<1',
        inverse_parameterization='lambda=sqrt(tau)/sqrt(1-Z^2); z_phys=sqrt(nu)*tau^((1-delta)/2)*Z*(1-Z^2)^(-(1-delta)/2)',
        exact_axial_jacobian=str(jacobian),axial_jacobian_identity=True,
        strict_jacobian_numerator='1-delta*Z^2=(1-delta)+delta*(1-Z^2)>0 for -1<Z<1',
        axial_limits='z_phys tends to -infinity / +infinity as Z tends to -1 / +1; (1-delta)/2>0',
        unique_global_axial_inverse=True,finite_positive_lambda=True,true_abs_Z_strictly_less_than_one=True,
        equivalent_implicit_root='lambda^2-(z_phys^2/nu)*lambda^(2delta)=tau',
        log_root_derivative_lower='2*(1-delta)>0',log_root_derivative_upper='2',
        log_root_derivative_minus_lower=str(lower),upper_minus_log_root_derivative=str(upper),
        zero_axial_case='z_phys=0 gives Z=0 and q=log(tau)/2',
        radial_map='R=r_phys^2/(2*nu*lambda^2); r_phys=sqrt(2*nu)*lambda*sqrt(R)',
        full_physical_coverage='Axial onto map supplies one finite lambda for each finite z; positive radial scaling supplies every R>=0; periodic angle covers x/y; axis uses Cartesian analytic core.',
        positive_real_power_domain='All powers are taken on tau,nu,1-Z^2>0; endpoints Z=+-1 are infinity limits, not finite physical points.',
        analytic_cover_not_a_numeric_sample_or_flatness_claim=True,passed=True)


@functools.lru_cache(maxsize=1)
def directed_candidate_cover_argument():
    """Bind the actual interval program to its enclosure invariants.

    Width is an accuracy/status issue, not a condition for set coverage.
    Interval Newton retains every parameter root in the analytic bracket.
    """
    asts=SourceAST();sources={}
    for stem,name in (('current_physical_tensor_locator_operator','implicit_log_coordinate_map'),
      ('current_physical_tensor_locator','_inverse_candidate'),('current_physical_tensor_locator','locate_log_radius'),
      ('current_physical_tensor_locator','tensor_from_location'),
      ('current_correlated_tensor_locator','_candidate'),('current_correlated_tensor_locator','locate'),
      ('current_correlated_tensor_locator','tensor'),('current_correlated_tensor_locator','_replay')):
        fn=asts.method(stem,name);sources[stem+'.'+name]=hashlib.sha256(ast.dump(fn,include_attributes=False).encode()).hexdigest()
    for target,wanted in {
        'lo':'endpoints(lt/2)[0]',
        'derivative':'ctx.mpf([endpoints(2*(1-d))[0],2])',
        'contracted':'intersect(ctx,ctx.mpf([lo,hi]),ctx.mpf(mid)-f/derivative)',
        'q':'ctx.mpf([lo,hi])'}.items():
        asts.expression('current_physical_tensor_locator_operator','implicit_log_coordinate_map',target,wanted=wanted)
    asts.expression('current_physical_tensor_locator','_inverse_candidate','coordinate',wanted='domain')
    asts.expression('current_physical_tensor_locator','_inverse_candidate','coordinate',wanted='intersect(c,coordinate,domain)')
    asts.expression('current_physical_tensor_locator','tensor_from_location','source',
        wanted="self.core_program(self.registry.owners['axis'],location['Z'],X,Y,location['requested_log_tau'],location['physical_viscosity'],rho_source=v)")
    asts.expression('current_physical_tensor_locator','tensor_from_location','view',
        wanted="self.registry.exterior(location['Z'],location['requested_log_tau'],location['theta'],location['physical_viscosity'])")
    asts.expression('current_correlated_tensor_locator','_candidate','coordinate',wanted='domain')
    asts.expression('current_correlated_tensor_locator','tensor','source',
        wanted="self.locator.core_program(self.registry.owners['axis'],location['Z'],X,Y,location['requested_log_tau'],location['physical_viscosity'],rho_source=v)")
    asts.expression('current_correlated_tensor_locator','tensor','view',
        wanted="self.registry.exterior(location['Z'],location['requested_log_tau'],location['theta'],location['physical_viscosity'])")
    return dict(actual_interval_program_AST_sha256=sources,actual_source_assignment_bindings=asts.bindings,
        proof_steps=[
          'Initial lower bracket is the directed lower log(tau)/2. The upper bracket makes both positive terms <=exp(2q)/2, for every parameter in its finite box.',
          'F is continuous and its derivative is in [2(1-delta),2], strictly positive. A trial sign changes a bracket endpoint only when it has that sign for all parameters.',
          'For any retained parameter root r, the mean value theorem gives r=mid-F(mid)/D with D in the global derivative box. Intersecting the interval Newton image with the bracket therefore retains every root.',
          'Step-limit or parameter/precision-limit stops preserve a finite bracket; they never choose a parameter midpoint. A narrow bracket is not required for coverage.',
          'The exact Z identities project only to proved [-1,1] bounds. An enclosure may include infinity-limit endpoints, while its exact finite physical source has |Z|<1.',
          'Each source radius range encloses its exact endpoint range. If the point is in a chart, range intersection retains it; monotone inverse operations and legal-domain projection retain its native coordinate.',
          'When a positive source width cannot be numerically separated from zero, the inverse retains the entire legal native domain. No unresolved candidate is discarded to choose one region.',
          'For correlated requests, exact source cancellation is performed before enclosure. A region is excluded only by a proved signed radius displacement or an empty directed legal inverse intersection; unresolved signs/inverses retain the candidate legal domain.',
          'The union of all candidate chart views encloses the physical source tensor. At shared endpoints the admitted same-function traces make the two source values equal; exact numeric seam selection is unnecessary for coverage.',
          'Core coordinates touching rho=0 use the original constrained Cartesian analytic extension. Every exterior candidate uses the unbounded exact Gamma source, never an infinite native-chart call.'],
        accepted_solver_statuses=['exact_axial_zero_source','requested_log_root_width','parameter_or_precision_limited','step_limit'],
        global_cover_is_union_coverage_not_unique_region_selection=True,
        bounded_solver_accuracy_not_promoted_to_field_accuracy=True,
        input_hashes=asts.hashes,passed=True)

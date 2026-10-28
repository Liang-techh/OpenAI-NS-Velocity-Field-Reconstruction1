"""Source-defined pulse endpoint normal forms, before inlet transport.

This module does not construct or evaluate an archived source graph.  It
checks pinned source definitions, then derives both radial moment endpoints
over exact symbols.  Inlet replacement and actual Z jets remain separate gates.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
KERNEL = PREFIX + 'n1_pulse_kernel_functions.py'
OWN = PREFIX + 'n1_pulse_own_fields.py'
RECEIPT = PREFIX + 'n1_pulse_own_fields_check.json'
EXPECTED = {
    KERNEL: '6ea98ad47f0f3f86307067627ac6af9af9201efe320cf0191c06e641d4f7368b',
    OWN: '1400066d2da5373a95f6cd30edad19359256211d0358b454203a92fac53c5ced',
}


def recover(root: Path):
    folder = root / 'experiments/root_st073'
    receipt = json.loads((folder / RECEIPT).read_bytes())
    if receipt.get('all_passed') is not True:
        raise ValueError('own-fields receipt did not pass')
    code = {}
    for name, expected in EXPECTED.items():
        raw = (folder / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected or receipt['input_hashes'].get(name) != expected:
            raise ValueError('source/receipt hash mismatch: ' + name)
        code[name] = raw.decode('utf-8').replace('\r\n', '\n')
    # These literals are anchors to the audited original definitions, not an
    # invocation or replay of the corresponding source constructors.
    anchors = {
        KERNEL: [
            "lam = [g.sub(g.constant('1/2'), g.mul(g.constant(i), mu)) for i in (1, 2)]",
            "g.unary('exprel', g.neg(g.mul(mu, s))))",
            "logE = g.sub(logpref, logmu)",
            "g.nodes[xi.node]['name'], g.zero, g.constant(11)), mu, 'positive original mu'))",
            "gp = g.mul(primitive, sigma(g.sub(g.constant(11), xi)))",
            "definition='0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside'",
        ],
        OWN: [
            "y=(g.add(g.quotient(c(13),mu,'positive original mu'),query) if chart=='end'",
            "else g.quotient(query,mu,'positive original mu'))",
            "kernel=integral(g.mul(exponential(g.neg(g.mul(lam,g.sub(y,v)))),gp_at(g.mul(mu,v))),",
            "weights=(W['first_row'],[g.add(A,g.mul(mu,D)) for A,D in zip(W['first_row'],W['divided_row'])])",
            "factor=exponential(g.add(W['log_end_scale'],g.mul(lam,distance)))",
            "logR=g.add(logRp,y);R=exponential(logR)",
            "m=a.mul(theta,means[0]);my=a.add(U,a.neg(m))",
            "B=a.scale(selected['ap'],gp_at(query))",
            "By=a.scale(selected['ap'],g.mul(mu,gp_prime(query)))",
            "B=a.fixed(0);By=a.fixed(0)",
        ],
    }
    for name, pieces in anchors.items():
        for piece in pieces:
            if piece not in code[name]:
                raise ValueError('source definition anchor missing: ' + piece)

    mu = s.symbols('mu', positive=True)
    lam, alpha, xi, q = s.symbols('lambda alpha xi q', real=True)
    J, M, ap, A0, A1, D0, D1, C0, C1 = s.symbols('J M ap A0 A1 D0 D1 C0 C1')
    E = s.exp(-11 * lam / mu)
    W = s.exp(-13 * lam / mu - alpha)
    H = s.exp(alpha - s.log(mu) + 2 * lam / mu)
    Pi = W * J / mu
    # Main integral: xi=mu*v, so dv=dxi/mu. No quadrature is performed.
    main = E * (M + ap * J / mu)
    gap = -H * (A0*C0 + A1*C1)
    selected_first_row = -mu * (W*M + ap*Pi)
    checks = {
        'source_Pi_factorization': s.simplify(s.exp(lam*(xi-13)/mu-alpha) - W*s.exp(lam*xi/mu)) == 0,
        'source_main_integrand_change_of_variable': s.simplify(s.exp(-lam*(11/mu-xi/mu)) - E*s.exp(lam*xi/mu)) == 0,
        'source_gap_factor': s.simplify(H*mu*W-E) == 0,
        'first_row_endpoint_given_selector': s.simplify(-H*selected_first_row-main) == 0,
        # exprel(q)=(exp(q)-1)/q with its continuous value 1 at q=0.
        'finite_mu_divided_row_identity': s.simplify(1 + (s.exp(q)-1) - s.exp(q)) == 0,
        'lambda2_difference': s.expand((s.Rational(1,2)-2*mu) - ((s.Rational(1,2)-mu)-mu)) == 0,
        'second_row_linear_weights': s.expand((A0+mu*D0)*C0+(A1+mu*D1)*C1-(A0*C0+A1*C1+mu*(D0*C0+D1*C1))) == 0,
        'second_row_endpoint_given_selector': s.simplify(-H*selected_first_row-main) == 0,
    }
    return {
        'schema': 'frame_pulse_two_moment_endpoint_v1',
        'status': 'source_definitions_bound_conditional_two_moment_normal_forms',
        'checks': checks,
        'passed': all(checks.values()),
        'source_hashes': EXPECTED,
        'receipt_source_family': receipt['source_family'],
        'definitions': {
            'lambda_k': '1/2-k*mu, k=1,2',
            'J_k': 'integral_0^11 exp(lambda_k*xi/mu)*gp(xi) dxi; not evaluated',
            'Pi_k': 'W_k*J_k/mu',
            'main_normalized_endpoint': 'exp(-11*lambda_k/mu)*(M_k+ap*J_k/mu)',
            'gap_normalized_endpoint': '-exp(alpha-log(mu)+2*lambda_k/mu)*sum(C_j*row_kj)',
            'row_1j': 'A_j',
            'row_2j': 'A_j+mu*D_j = integral beta_j(s)*exp(lambda_2*s) ds',
            'common_radius': 'exp(logRp+11/mu), from same main/gap formula',
            'common_theta': 'Pstar*native_u*exp(-(1/2+mu)*11/mu), before inlet rewrite',
            'flat_endpoint': 'sigma(0)=sigma_prime(0)=0 implies gp(11)=gp_xi(11)=0; gp_y=mu*gp_xi',
        },
        'premises': ['mu>0', 'the selected row identities hold for the actual functions',
                     'the fixed integral definitions have their stated semantics'],
        'not_proved': ['actual rewritten endpoint graph semantics', 'capture-free inlet substitution',
                       'ordinary Z derivative transport or C2 regularity',
                       'actual frame endpoint identity', 'numerical values or error reduction'],
        'source_constructor_called': False,
        'original_source_evaluated': False,
    }

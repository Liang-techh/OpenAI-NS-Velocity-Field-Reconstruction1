# Original flat shapes and actual O.4 mixed derivatives

The unchanged compliant .001 source now supplies sigma, gp and normalized
end beta derivatives through radial order four, including boxes crossing
their compact support endpoints. Those data feed the selected ap/end
functions and source ODEs, providing all 15 (y,Z) mixed derivatives of each
physical velocity component and original pressure with total order <=4.
Here y=log(R), R is the paper similarity coordinate, and Z is its axial
coordinate. Physical radial/swirl/axial velocity factors are differentiated
before axial derivatives. This is not yet a Cartesian derivative certificate.

The work completes F38-D3c1g1a through g1d. Quantitative flat shape support
majorants are also installed. Full pulse interface/regularity acceptance,
post-pulse mixed C4, whole outer C4, physical energy, stress cone, admissible
stress/flat remainder and genuine temporal recursion remain incomplete.

## Original support derivatives

For the original cutoff

    sigma(x)=exp(-1/x²)/(exp(-1/x²)+exp(-1/(1-x)²)), 0<x<1,

the code retains its exact constant extensions, left logistic expression
and reflection sigma(1-x)=1-sigma(x). With m=min(x,1-x), q=sigma(1-sigma),
q<=exp(4-m^-2), and the first four logarithmic odds derivatives have bounds
|Lj|<=2(j+1)!m^(-j-2). The resulting ordinary derivative envelopes are

    |sigma^(n)| <= Cn exp(4-m^-2)m^(-3n),
    C1=4, C2=28, C3=256, C4=3104.

For a support-tail box 0<m<=W, maximize this expression at
min(W,sqrt(2/(3n))). The order-zero tail bounds sigma near zero and
1-sigma near one. All exact analytic envelopes tend to zero at the endpoint.

For beta(r)=exp(-1/w), w=1-r², |r|<1:

    beta^(n)=beta Rn(r)/w^(2n), R0=1,
    R_(n+1)=w² Rn' +(4n r w-2r)Rn.

The coefficient norm An=||Rn||1 bounds |Rn| on [-1,1]. The envelope
An exp(-1/w)w^(-2n) peaks at w=1/(2n). A tail box uses
min(W,1/(2n)); using a square-root maximizer would be incorrect.
Every derivative has its exact zero extension outside support. The end
shape beta(s/.15)/(.15*normalization) keeps the original positive
normalization and the derivative scale .15^(-(n+1)).

For gp(xi)=P(xi)*sigma(11-xi), P'=sigma(50xi), its primitive is exactly
xi-.01 for xi>=.02 by original cutoff symmetry. Below .02 the original
positive directed quadrature is retained. P^(n)=50^(n-1)sigma^(n-1)(50xi)
and Leibniz/exit signs supply all required derivatives. Broad boxes split at
0, .02, 10 and 11. Main y derivatives multiply the xi derivatives by mu^n.
End derivatives use the exact retained s offsets, never a rounded xi=13.

Exponentials with a huge negative log receive a positive outward numerical
cap. The cap is only an enclosure of the original positive function;
neither the source definition nor its analytic flat limit is changed.
The tail API uses W*exp(-1000) only after proving the full log envelope
is <=-1000+log(W). Its returned positive upper cap therefore also tends
to zero with endpoint distance; it has no fixed numerical floor. Sigma
order zero always means deviation from the respective endpoint constant.

## Source ODE recovery and physical factors

The five existing partial histories are retained as functions of Z:

    (m1)_y=B-(.5-mu)m1,
    (m2)_y=B-(.5-2mu)m2,
    X_y=1-(1-mu)X,
    e_y=B²-.5+2mu e,
    P_y=Utheta²/2.

Repeated differentiation uses the exact binomial derivatives of B².
Original P0+Mp supplies pressure order zero; all higher y derivatives
come from the same swirl balance. Inactive gaps have B and all its
derivatives zero while accumulated moment/energy/angular histories remain.

Paper (3.9) linearly recovers A(B,m1) from the same Mz derivative, with
inlet Utheta_Z/Utheta=-2Z/(1+Z²). Fifth axial input permits A through axial
order four for every needed y order. Utheta has y rate -(.5+mu), whereas
sqrt(R/2)*Utheta has y rate -mu. Thus the actual velocity derivatives use

    Uz_y^k/Utheta = sum_j binom(k,j)(-.5-mu)^(k-j) B_y^j,
    Ur_y^k/(sqrt(R/2)Utheta) = sum_j binom(k,j)(-mu)^(k-j) A_y^j.

The inlet axial factor is multiplied BEFORE taking Z derivatives.
Reported mixed values multiply ordinary Taylor coefficient n by n!.
All 15 multiindices k+n<=4 are stored separately for Uz, Utheta, Ur and P.
Common exact positive radial factors remain described in the segmented
source coordinates; factored values are not full physical norms.

## Evidence and reproduction

Under experiments/root_st073, with prefix lei_ren_part1_paper_:

- compliant_flat_pulse_derivatives.py/.json and _check.py/.json;
- compliant_pulse_mixed_C4.py/.json and _check.py/.json;
- compliant_reconstruction.py --stage flatjets;
- compliant_reconstruction.py --stage mixedjets.

The ordered pipeline contains 59 modules. Earlier source modules and
receipts remain unchanged, and transitive source/family hashes are checked.

Independent symbolic differentiation checks all five beta polynomials,
the coefficient norms, both envelope maximizers, sigma term constants and
all needed analytic flat limits. Independent original-shape point and
support-crossing derivative fixtures check sigma/beta/gp, including a box
spanning every gp chart. End normalization powers and tiny positive caps
are checked separately. The finite-precision shape fixture allowance is
explicitly separate from actual directed source admission.

An independent closed-form fixture solves the moment ODEs analytically,
then differentiates the physical velocities and primitives directly in
both variables. It checks 480 mixed derivatives, including the true
sqrt(R/2) and Utheta product rules. Actual-source checks cover all 15
multiindices in every sampled chart and broad whole-Z entrance/main/exit,
both gap and end boxes covering the entire O.4 interval,
earlier derivative agreement, three chart overlaps, exact terminal Uz/Ur
mixed zeros and retained positive selected future energy.

A read-only GPT-5.6 Luna/max worker reviewed the derivative envelope and
transport formulas. These derivatives belong to the leading profile;
they are not n-dependent coefficient recursion or evidence of blow-up.

## Next executable work

1. Convert the complete factored entrance/main/exit/gap/end box bounds
   into a scale-aware physical bound ledger. Compose exact source identities
   across all interfaces, with quantitative flat velocity compatibility;
   finite sample overlap alone is insufficient for full pulse C4.
2. Extend high mixed derivatives through the corrected post-Rv flatten,
   angular bumps, steep transitions, waiting, epsilon collar and entire
   Gamma exterior while preserving absolute pressure/angular closure.
3. Accept whole-field C4 only with core/inner/outer interfaces. Compute
   physical energy and stress-cone margins from complete physical fields.
4. Construct admissible stress and independently flat remainder; recover
   actual n=1 and n>=2 coefficient equations and separate repairs before
   oscillatory cancellation and full Cartesian residual diagnostics.

No cell Taylor polynomial is used here, and no missing cell remainder is
inferred from jet length. Full pulse/outer C4 and all later-stage acceptance
flags remain false.

# Lei–Ren Part I regular core replacement

This note extracts the actual inner construction in the cached Part I source
`work_paper_cache/lei_ren_part1.txt`. The reference branch is

\[
 U^\theta_{\rm ref}(R,Z)=\frac{P_*}{1+Z^2}
       \left(\frac R{R_{\rm ref}}\right)^{1/10},\qquad
 U^z_{\rm ref}(R,Z)=4Z.
\]

The regular core replaces the nonregular angular exponent `1/10` by
\(U^\theta=\sqrt{2R}\,F\), with \(F\) analytic and positive at the axis,
and initially uses \(U^z=4Z+j+O(R)\). The shift \(j>0\) is structural: it
separates the zero of the axial transport coefficient from the even pressure
symmetry point. The final gluing budget fixes a much smaller `j`; the value
`j=1/20` in the linear discussion is only illustrative.

The cached Part I file contains Sections 8–10 but no Appendix A heading. Its
only appendix references are to Appendix B/C of reference [21], so no local
Appendix A equations are available to quote.

## Actual regular core equations

Section 8.1, equations (8.1)–(8.5), defines

\[
 d=1-Z^2,\quad L=1-\delta Z^2,\quad
 W=1-\frac{(1-\delta)Z M^z+d\,\partial_ZM^z}{R},\quad
 H=\frac{1-\delta}{2}Z+dU^z,
\]

with
\[
 M^z(R,Z)=\int_0^R U^z(\rho,Z)\,d\rho,\qquad
 \partial_RP=F^2,\qquad
 P(R,Z)=P_0(Z)+\int_0^R F(\rho,Z)^2\,d\rho.
\]

The stress-free nonlinear core solves

\[
\begin{aligned}
2L(RF_{RR}+2F_R)&=W(F+RF_R)+\frac\delta2(1-2ZU^z)F+HF_Z,\\
2L(RU^z_{RR}+U^z_R)&=WRU^z_R+\frac{1+\delta}{2}(1-2ZU^z)U^z+HU^z_Z\\
&\qquad+dP_Z-2(1+\delta)ZP-2ZRF^2.
\end{aligned}
\]

Regularity at the axis gives the exact first jets (8.7). For
\(F_0=F(0,Z)\), \(U^z_0=U^z(0,Z)\), and
\(H_0=(1-\delta)Z/2+dU^z_0\),

\[
 F_R(0,Z)=\frac{(1+\delta/2-ZU^z_0-dU^z_{0,Z})F_0+H_0F_{0,Z}}{4L},
\]

\[
 U^z_R(0,Z)=\frac{\frac{1+\delta}{2}(1-2ZU^z_0)U^z_0
+ H_0U^z_{0,Z}+dP_{0,Z}-2(1+\delta)ZP_0}{2L}.
\]

The axis data used by the analytic construction are (8.9)

\[
 U^z_0(Z)=4Z+j,\qquad U^z_{0,Z}=4,\qquad 0<j\leq\frac1{20}.
\]

Let \(Z_0\) be the unique zero of \(H_0\). With
\(\sigma_0=j/500\),

\[
 \chi(Z)=\frac{H_0(Z)^2}{H_0(Z)^2+\sigma_0^2},\qquad
 G(Z)=\int_{Z_0}^Z
 \frac{L(w)H_0(w)}{H_0(w)^2+\sigma_0^2}\,dw,\qquad
 F_0(Z)=C_*^{-1}e^{-\Lambda G(Z)}.
\]

The linear model (8.20) has

\[
 2(RF_{RR}+2F_R)=-(\Lambda\chi+\beta)F,\qquad
 U^z_R=-\frac{g(Z)}{2L},
\]

where \(\beta=(3-\delta/2+jZ)/L\) and `g` is the pressure/axial forcing in
equation (8.16). Its explicit regular solution (8.22) is

\[
 F^{\rm m}=F_0 B(q),\quad q=\frac R2(\Lambda\chi+\beta),\quad
 B(q)=\sum_{n\ge0}\frac{(-q)^n}{n!(n+1)!},
\]

\[
 U^{z,\rm m}=4Z+j-\frac{g(Z)}{2L}R.
\]

The theorem's actual nonlinear solution is written (8.47) as

\[
 F=F_0(Z)\Phi(\Lambda R,Z),\qquad
 U^z=U^z_0(Z)+\Lambda^{-1}\Psi(\Lambda R,Z),\qquad
 \Phi(0,Z)=1,\ \Psi(0,Z)=0,
\]

and is obtained by an analytic contraction of the exact rescaled equations
(8.49)–(8.53). The first angular correction (8.54) is retained because the
leading angular derivative vanishes at \(Z_0\); estimate (8.57) then
preserves \(F_R<0\) there.

The free parameter restrictions are substantive:

\[
 \Lambda\ge\max\{500,j^{-2}\},\qquad
 C_*\ge\Lambda^2e^{\Lambda A_\Omega},\qquad
 R_a=4/\Lambda.
\]

Theorem 8.1 gives a jointly analytic solution on
\([0,4.1/\Lambda]\times[-1,1]\), \(F>0\), \(F_R<0\), and
\(\kappa(R_a,Z)\ge9/4>2\). The core has zero total stress,
\(\mathcal T=0\), and uses the prescribed axis pressure \(P_0\); it does
not determine a new pressure datum by integrating backward from the heat
tail.

## Section 9 matching jets and connection

Set \(r=R_a\), \(f(Z)=F_c(r,Z)\), \(v(Z)=U_c^z(r,Z)\), and let
\(m_j(Z)=M_{j,c}(r,Z)\). The core-to-reference connection has two exact jet
requirements (Section 9.1):

* at \(R_a\), every radial jet of \(F,U^z\) agrees with the core:
  \(\partial_R^kF(R_a,Z)=\partial_R^kF_c(R_a,Z)\) and
  \(\partial_R^kU^z(R_a,Z)=\partial_R^kU_c^z(R_a,Z)\) for every \(k\);
* at \(R_h=e^{-5}R_{\rm ref}\), every jet agrees with
  \(U^\theta_{\rm ref},4Z\).

The five moment targets at \(R_h\) are (9.2)

\[
\begin{aligned}
M_o^\theta&=\frac58R_h\sqrt{2R_h}\,u_h,&M_o^z&=4ZR_h,\\
M_o^{\theta z}&=4ZM_o^\theta,&
M_o^{z\theta}&=16Z^2R_h-\frac5{12}R_hu_h^2,&
M_o^p&=\frac52u_h^2,
\end{aligned}
\qquad u_h=\frac{e^{-1/2}P_*}{1+Z^2}.
\]

The connection must satisfy the five integral identities (9.3), equivalently
the actual cumulative moments at \(R_h\) must equal the outer moments.
Velocity matching alone does not imply this. Equation (9.9) requires

\[
\|v-4Z\|_{C_Z^2}+\|m_z/r-4Z\|_{C_Z^2}\le\epsilon_0,\qquad
H_v\partial_Z\log f\le\frac1{20},
\]

where \(H_v=(1-\delta)Z/2+(1-Z^2)v\). The frozen-profile test (9.14)
requires \(D_f>0\), \(H_f\ge2+4\gamma\) on \([r,110]\), and
\(D_f\ge4\) on \([100,110]\), using the actual frozen moments (9.12) and
the exact inertial formulas (9.13).

The flat inner collar is explicit. For \(y=\log(R/r)\), first define the
comparison profile on \([0,2h_b]\) by (9.23),

\[
 \alpha(y)=1-\sigma((y-h_b)/h_b),\quad
 \partial_y\log\bar F=\alpha\partial_y\log F_c,\quad
 \partial_y\bar V=\alpha\partial_yV_c,
\]

with \((\bar F,\bar V)(0)=(f,v)\), then freeze it for larger \(y\). Its
normalized inertial direction is \(\mathbf q=(\bar D,\bar E)\). The actual
connection on \(r<R\le100\) uses (9.25)–(9.26),

\[
 \chi_b(y)=1-(1-\varepsilon_b)\sigma(y/h_b),\quad
 \partial_y\log F=-\frac{\chi_b}{2}\bar D,\quad
 \partial_yV=-\chi_b\sqrt{R/2}\,F\bar E.
\]

Equivalently,

\[
 F=f\exp\left[-\frac12\int_0^y\chi_b\bar D\,dt\right],\qquad
 V=v-\int_0^y\chi_b\sqrt{re^t/2}\,F\bar E\,dt.
\]

Because \(1-\chi_b\) is flat at zero, all core radial jets match at
\(R_a\). The error estimate (9.27) retains the vanishing factor
\(1-\chi_b\), which is essential because the stress itself vanishes at the
core endpoint. Later steps reshape the swirl by (9.30) to the power profile,
then restore \(V=4Z\) by (9.38). The connection only proves a relaxed cone;
the admissible collar is the initial subinterval.

The parameter ordering and scale restrictions are also part of the input:

\[
 R_{\rm ref}=110(C_*P_*)^{10},\quad C_*\ge e^{4A},\quad
 R_{\rm ref}\ge110e^{T+10}(1+A)^{10},\quad
 R_z=e^{-8}R_{\rm ref}\ge110(1+K)^2/P_*^2.
\]

These are hypotheses under which the estimates in Sections 9–10 apply.

## Section 10 exact five-moment correction

The correction is supported strictly inside the terminal reference interval.
Set (10.2)

\[
 R_m=e^{-6}R_{\rm ref},\qquad R_h=eR_m,\qquad x=R/R_m,\qquad
 A_m^\theta=\frac{e^{-3/5}P_*}{1+Z^2}.
\]

Use the fixed nonnegative bump \(\beta=\beta_{1/40}\), supported in
\([-1/40,1/40]\), with centers

\[
 s_1=5/4,\quad s_2=3/2,\quad s_3=7/4,\qquad
 \gamma_j(x)=\beta(x-s_j),\qquad b_1=\gamma_1,\ b_2=\gamma_3.
\]

For coefficient functions \(h=(c_1,c_2,\xi_1,\xi_2,\xi_3)\),

\[
 f_h=\sum_{j=1}^3\xi_j\gamma_j,\qquad
 g_h=\sum_{i=1}^2c_i b_i,
\]

and on \(R_m<R<2R_m\),

\[
 \widehat u=A_m^\theta(x^{1/10}+f_h),\qquad
 \widehat V=4Z+g_h.
\]

The centered normalized defect is (10.3)

\[
\mathbf d=\left(
\frac{\Delta_z}{R_m},
\frac{\Delta_{\theta z}-4Z\Delta_\theta}{\sqrt2R_m^{3/2}A_m^\theta},
\frac{\Delta_\theta}{\sqrt2R_m^{3/2}A_m^\theta},
\frac{\Delta_{z\theta}-8Z\Delta_z}{R_m(A_m^\theta)^2},
\frac{\Delta_p}{(A_m^\theta)^2}\right),\qquad
\mathfrak e=\|\mathbf d\|_1.
\]

The exact five-component map (10.8), with all integrals over `1 < x < 2`, is

\[
\mathbf F_Z(h)=
\begin{pmatrix}
\int g_h\\
\int x^{3/5}g_h+\int\sqrt{x}f_hg_h\\
\int\sqrt{x}f_h\\
-\int x^{1/10}f_h+\int((A_m^\theta)^{-2}g_h^2-f_h^2/2)\\
\int x^{-9/10}f_h+\int f_h^2/(2x)
\end{pmatrix}
=\mathsf A h+\mathcal Q_Z(h,h).
\]

The two axial columns are separated supports with weights `1` and `x^(3/5)`;
the three angular columns have weights `x^(1/2)`, `-x^(1/10)`, and
`x^(-9/10)`. The resulting fixed matrix is invertible. With constants
`C_A`, `C_Q`, the exact contraction is

\[
 h\longmapsto-\mathsf A^{-1}[\mathbf d+\mathcal Q_Z(h,h)],\qquad
 \|h\|_1\le2C_A\mathfrak e,
\]

provided \(\mathfrak e\le\mathfrak e_*\), where
\(\mathfrak e_*\le(8C_A^2C_Q)^{-1}\). It restores all five moments while
leaving the core and its inner collar unchanged. The paper claims only the
relaxed cone on the bump support; the admissible cone is inherited away from
it.

## Existing extended pressure core versus the source parameters

The current provisional implementation is
`experiments/root_st073/lei_ren_part1_extended_pressure_core.py`, backed by
`src/openai_ns_reconstruction/paper_core_series.py` and
`paper_core_reference.py`. Its saved candidate
`lei_ren_part1_extended_pressure_core.json` uses

~~~text
h = .001, j = .02, sigma = .3, Lambda = 10, C = 2,
pressure_scale = 1, eta_nodes = 257, radial_degree = 4,
R_core = .005, R_anchor = .01.
~~~

It is a finite degree-4 floating-point radial prefix and a shared-pressure
iteration into the current heat collar (`R_b=2048`, inner radius about
`1242.17`). The receipt records a collocation pressure defect `4.31e-10`
and an off-grid holdout defect `7.05e-9`, but explicitly records
`whole_cone_validated=false` and `five_moments_closed=false`. Those numbers
only establish self-consistency of that finite provisional candidate.

It cannot serve as the paper's same-pressure core for the current source
parameters `logPstar=14`, `logRref=10`:

1. The paper requires \(P=P_0+\int_0^RF^2\) with the *completed outer*
   \(P_0=-\int_0^\infty(U_o^\theta)^2/(2R)\,dR\). The source reference
   contribution alone is
   \(M^p(R_{\rm ref})=(5/2)P_*^2/(1+Z^2)^2\), with \(P_*=e^{14}\), hence
   order \(10^{12}\) pressure. The corrected-source receipt gives normalized
   \(P_0/P_*^2\approx-2.78985\) at `Z=.3`, about
   \(-4.0\times10^{12}\) in physical source units. The saved extended core
   has `pressure_scale=1` and axis pressure about `-0.0989` at `Z=0`; it
   never receives the source `P_0`.
2. Its parameters violate the analytic-core hypotheses:
   `Lambda=10 < 500`, `Lambda < j^(-2)=2500`, `C=2 < Lambda^2`, and
   `j=.02` is far larger than the final Section 9 budget
   \(j=\eta_{\rm tol}/8\). The finite `PaperCoreSeries` is explicitly a
   diagnostic prefix, not the jointly analytic nonlinear solution of (8.2).
3. The source-scale compatibility already fails before cone or moment
   questions: Section 9 requires \(R_{\rm ref}=110(C_*P_*)^{10}\).
   With `logPstar=14` and `C_*>=1`, its logarithm is at least
   \(\log(110)+140\approx144.7\), whereas the current candidate has
   `logRref=10`.
4. The radial geometry is different: the theorem's
   \(R_h=e^{-5}R_{\rm ref}\approx148.4\) and
   \(R_m=e^{-6}R_{\rm ref}\approx54.6\) for `logRref=10`, while the current
   heat-collar interface is near `1242`. Its current blend and pressure
   receipt do not implement the Section 9 matching jets or Section 10 bump
   interval.

Thus the existing extended pressure core is useful as a finite numerical
pressure-iteration scaffold, but it is not a same-pressure provisional core
for the `logPstar=14`, `logRref=10` corrected source. A valid source-bound
core would first need the completed source \(P_0(Z)\), the analytic core
parameter regime, and the Section 9 jet/scale compatibility; the five
Section 10 moment defects must then be measured from the actual connection,
not replaced by target values.

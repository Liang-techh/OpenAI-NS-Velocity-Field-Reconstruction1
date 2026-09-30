# Next source connection implementation

Pinned source: arXiv2609.35406v1, Sections 9.4-9.6. The current actual
Step 2 continuation reaches R=100; comparison alone reaches R=110.

Define a=-2 d_y log(F), b=2 d_y V/Utheta, Utheta=sqrt(2R)F.
On 100<=R<=100 exp(hb), keep a=epsilon*barD and set
b=-epsilon*barE*(1-sigma(log(R/100)/hb)). On the next interval up to
100 exp(2hb), put s=log(R/(100 exp(hb)))/hb, b=0 and
a=(1-sigma(s))*epsilon*barD+(4/5)*sigma(s).
Integrate d_y logF=-a/2 and d_y V=b*sqrt(2R)*F/2, together with actual
moments and Z jets. Then a=4/5,b=0 until R=110. Continue moments;
do not replace them by endpoint targets. Check D(110)>3 numerically.

The long reshape has Treshape=400A and Rsh=110 exp(Treshape).
For y=log(R/110),
Utheta=exp(y/10)*u1^(1-sigma(y/Treshape))
        *[1/(Cstar*(1+Z^2))]^sigma(y/Treshape), V=v1.
After Rsh, retain the reference angular power law. Restore V to 4Z
on Rz..eRz, where Rz=exp(-8) Rref, via the source smooth switch.
The actual five repair defects at Rh=exp(-5)Rref are
Delta_j=M_j(Rh)-M_outer,j(Rh). Section 10 must repair these actual
defects on a subinterval of Rm..Rh, Rm=exp(-6)Rref.

## Parameter gate before the long reshape

The source requires Rref=110*(Cstar*Pstar)^10,
Cstar>=exp(4A), and Rref>=110*exp(400A+10)*(1+A)^10,
as well as its earlier complex-domain core constraints. Therefore
100 exp(2hb)<=110<Rsh<Rz<eRz<Rm<2Rm<Rh must be checked for the
same input candidate. The existing logCstar=2 log(Lambda) candidate
was explicitly a necessary lower-bound demonstration; it does not
certify these inequalities. Do not blindly append the long reshape
to its current Rref. Rechoose shared parameters and recompute the
shared pressure/core if needed. The source also chooses
hb=epsilon=cstar*K^(-100); the current independently selected hb and
epsilon are numerical collar inputs, not this certified regime.

The next implementation should preserve explicit failures, sampled
relaxed versus admissible tests and separate finite-precision value
bounds from unproved Z/pressure bounds. Neither radial matching nor
these source profiles alone establish temporal scale recursion.

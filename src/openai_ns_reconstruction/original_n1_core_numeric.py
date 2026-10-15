"""Finite radial prefixes of the original, fully coupled temporal n=1 core.

This computes actual coefficient enclosures and numerical approximations of
finite prefixes. It is not a point oracle for the infinite radial solution:
the radial Taylor remainder is deliberately left open. In particular, the
existing Volterra-word remainder is not reused as a radial Taylor remainder.
"""
import ast
import hashlib
import json

import sympy as sp

from .original_core_numeric import OriginalCoreNumeric, ROOT, PREFIX, _bounds


RING_SOURCE = "lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py"
RAW_SOURCE = "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"


def pure_definitions(path, names, environment):
    """Load named pure definitions, without importing historical owners."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    definitions = [node for node in tree.body
                   if isinstance(node, (ast.ClassDef, ast.FunctionDef)) and node.name in names]
    if {node.name for node in definitions} != set(names):
        raise ValueError("Required pure source definitions changed")
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(path), "exec"), environment)
    return environment


def linear_rhs(v):
    """Exact conjugation of raw rows 4..6 in local F0(Z) units.

    f=epsilon*F1/F0, u=Uz1/epsilon, k=K1/epsilon,
    p=P1/(epsilon*Pstar^2). The F0 derivative enters exactly once.
    Operands can be symbolic scalars or finite Taylor rectangles.
    """
    z,r,de,e,P2,F2 = (v[key] for key in ("z","rho","delta","epsilon","P2","F2"))
    L,d,A,C,g = (v[key] for key in ("L","d","A","C","g"))
    E,H,B,Zf,Zw,Q = (v[key] for key in ("E","H","B","Zf_without_gauge","Zw","Q"))
    f,u,k,p,fz,uz,kz,pz,Ef,Eu = (v[key] for key in
        ("f","u","k","p","fZ","uZ","kZ","pZ","Euler_f","Euler_u"))
    rf = (2*e*(Q+E*(-2+de)/L)-2*H*g/L)*f
    rf += (e**3*(2*Zf-2*(-1+de)*z*C/L)-2*e**2*d*g*A/L)*u
    rf -= 2*e**3*(1+de)*z*C*k/L
    rf += 2*e*B*Ef+2*e*H*fz/L-2*e**3*d*C*(uz+kz)/L
    ru = -8*F2*z*r*A*f/L
    ru += 2*e*(E*(-1+de)/L+Zw-v["Euler_U"]*z*(-1+de)/L)*u
    ru -= 2*e*(1+de)*z*v["Euler_U"]*k/L
    ru -= 4*(e*P2)*z*p/L
    ru += 2*e*B*Eu+2*e*(H-d*v["Euler_U"])*uz/L
    ru -= 2*e*d*v["Euler_U"]*kz/L
    ru += 2*(e*P2)*d*pz/L
    rp = 2*F2/(e*P2)*A*f
    return rf,ru,rp


class OriginalN1CoreNumeric:
    """Source-bound n1 radial coefficients; finite-prefix evaluation only."""
    def __init__(self, *, precision=60, quadrature_nodes=128, directory=None):
        self.leading = OriginalCoreNumeric(precision=precision,
            quadrature_nodes=quadrature_nodes, directory=directory)
        self.ctx,self.iv = self.leading.ctx,self.leading.iv
        self.directory = self.leading.directory
        self.precision = precision
        self.family = self.leading.family
        self.hashes = dict(self.leading.hashes)
        self.code_hashes = {}
        accepted_path = self.directory/(PREFIX+"compliant_current_original_core_numeric_approximation.json")
        accepted = json.loads(accepted_path.read_text(encoding="utf-8"))
        if not accepted["checks"]["passed"] or accepted["source_family"] != self.family:
            raise ValueError("Accepted same-family numerical leading core required")
        self.hashes[accepted_path.name] = hashlib.sha256(accepted_path.read_bytes()).hexdigest()
        for relative,digest in accepted["implementation_hashes"].items():
            path = ROOT/relative
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError("Accepted numerical leading implementation changed")
            self.code_hashes[relative] = digest
        for name in (RING_SOURCE,RAW_SOURCE):
            path = self.directory/name
            self.hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.Rectangle = pure_definitions(self.directory/RING_SOURCE,
            ("TaylorRectangle",), {"sp":sp})["TaylorRectangle"]

    def _rhs(self, fixed, leading, state, n, axial, center, delta, eps, P2, c):
        """Coefficient RHS through rho^n, with two actual extra Z orders."""
        R = self.Rectangle
        def known(rows, radial, order):
            if len(rows) < radial+1 or any(len(row) < order+1 for row in rows[:radial+1]):
                raise ValueError("An unavailable mixed source coefficient may not be padded with zero")
            return R(rows,radial,order,c.mpf)
        def variable(radial,order):
            return R([[center,c.mpf(1)]],radial,order,c.mpf)
        def crop(value,radial,order):
            if radial > value.radial or order > value.axial:
                raise ValueError("Unavailable rectangle requested")
            return R(value.rows,radial,order,c.mpf)
        # Q needs one Z derivative of the same leading mean; Np also consumes
        # Q_(n+1). Do not discard that last radial row before differentiating.
        Aw = known(leading["A"],n+1,axial+2)
        Uw = known(leading["Uz"],n+1,axial+2)
        zw = variable(n+1,axial+2)
        mw = Uw.radial_average()
        cutw = lambda q:q.truncate_axial(axial+1)
        zz,uu,mm = map(cutw,(zw,Uw,mw))
        dd,LL = 1-zz*zz,1-zz*zz*delta
        Qw = (2*zz*uu-zz*mm*(1-delta)-dd*mw.dz())/LL
        A2,U2 = crop(Aw,n,axial+2),crop(Uw,n,axial+2)
        z = variable(n,axial)
        d,L = 1-z*z,1-z*z*delta
        A,U,Q = crop(Aw,n,axial),crop(Uw,n,axial),crop(Qw,n,axial)
        QZ = crop(Qw.dz(),n,axial)
        g2 = -known([fixed["ell_Z_taylor"]],0,axial+2)
        # Promote fixed Z jets by exact radial constancy, not source padding.
        g2 = R(g2.rows,n,axial+2,c.mpf)
        F2 = R(known([fixed["S_Z_taylor"]],0,axial).rows,n,axial,c.mpf)/(eps*eps)
        g = g2.truncate_axial(axial)
        rho = R([[0],[1]],n,axial,c.mpf)
        def T(weight,value):
            order = value.axial-1
            if order < 0:
                raise ValueError("Actual extra Z coefficient required by weighted derivative")
            zv = variable(n,order)
            val = value.truncate_axial(order)
            return (zv*val*weight+(1-zv*zv)*value.dz()-2*zv*val.radial_euler())/(1-zv*zv*delta)
        a = -2-delta
        m2 = (1-variable(n,axial+2)**2)*g2/(1-variable(n,axial+2)**2*delta)
        m = m2.truncate_axial(axial)
        TaA = T(a,A2)
        force_f = -2*m*m*A+(T(-3,(m2*A2).truncate_axial(axial+1))
            +m*TaA.truncate_axial(axial))*(2*eps)-T(-3,TaA)*(2*eps**2)
        force_u = -2*T(-2,T(-1-delta,U2))
        eq = Q.radial_euler()
        S = (Q+z*QZ*((1-delta)/2)+eq)/L+Q*(Q/2+eq)
        S += U*(-2*z*Q+d*QZ-2*z*eq)/L
        pressure_head = R([[(i+1)*(i+2)*v for v in Qw.rows[i+1][:axial+1]]
                           for i in range(n+1)],n,axial,c.mpf)
        Np = -S/(2*P2)+pressure_head/(eps*P2)
        force_u -= 4*z*rho*Np*(eps*P2)/L
        values = {name:known(state[name],n,axial+1) for name in ("f","u","k","p")}
        v = {name:q.truncate_axial(axial) for name,q in values.items()}
        v.update({name+"Z":q.dz() for name,q in values.items()})
        EU,EA = U.radial_euler(),A.radial_euler()
        v.update(z=z,rho=rho,delta=z.constant(delta),epsilon=z.constant(eps),P2=z.constant(P2),F2=F2,L=L,d=d,A=A,C=A+EA,g=g,
            E=z*U-c.mpf(1)/2,H=z*((1-delta)/2)+d*U,B=Q+(1-2*z*U)/L,Q=Q,
            Zf_without_gauge=(z*A*(-2-delta)+d*A2.dz().truncate_axial(axial)-2*z*EA)/L,
            Zw=(z*U*(-1-delta)+d*U2.dz().truncate_axial(axial)-2*z*EU)/L,
            Euler_f=v["f"].radial_euler(),Euler_u=v["u"].radial_euler(),Euler_U=EU)
        rf,ru,rp = linear_rhs(v)
        return rf+force_f,ru+force_u,rp+Np

    def _advance(self, fixed, leading, state, n, degree, depth, center, delta, eps, P2, c):
        axial = degree+depth-n-1
        rf,ru,rp = self._rhs(fixed,leading,state,n,axial,center,delta,eps,P2,c)
        k = n+1
        state["f"].append([v/(4*k*(k+1)) for v in rf.rows[n]])
        state["u"].append([v/(4*k*k) for v in ru.rows[n]])
        state["p"].append([v/k for v in rp.rows[n]])
        state["k"].append([-c.mpf(k)*v/(k+1) for v in state["u"][-1]])

    def coefficients(self,Z,*,degree=12,depth=1):
        if (type(Z) is not str or not -1 < self.ctx.mpf(Z) < 1
                or type(degree) is not int or not 3 <= degree <= 22
                or type(depth) is not int or not 1 <= depth <= 2):
            raise ValueError("Use an interior decimal Z, n1 radial degree3..22, depth1..2")
        packet = dict(Z=Z,degree=degree,depth=depth,source_family=self.family)
        for directed,label in ((False,"point"),(True,"directed")):
            c = self.iv if directed else self.ctx
            # N+1 leading rows with depth+1 contain exactly the mixed orders
            # needed by both the axial second source and pressure recovery.
            lead_degree,lead_depth = degree+1,depth+1
            fixed,leading,delta,eps = self.leading._seed(Z,lead_degree,lead_depth,directed=directed)
            for n in range(lead_degree):
                self.leading.advance(c,fixed,leading,n,c.mpf(Z),delta,eps)
            lp = self.leading.lp_iv if directed else self.leading.logP
            P2 = c.exp(2*lp)
            state = {name:[[c.mpf(0)]*(degree+depth+1)] for name in ("f","u","k","p")}
            for n in range(degree):
                self._advance(fixed,leading,state,n,degree,depth,c.mpf(Z),delta,eps,P2,c)
            packet[label] = state
        return packet

    def evaluate_prefix(self,packet,rho,*,degree=None):
        """Evaluate a polynomial prefix, never label its error a full n1 error."""
        c,iv = self.ctx,self.iv
        N = packet["degree"] if degree is None else degree
        if (packet["source_family"] != self.family or type(N) is not int or not 0 <= N <= packet["degree"]
                or type(rho) is not str or not 0 <= c.mpf(rho) <= 4):
            raise ValueError("Same-family available prefix and rho decimal in [0,4] required")
        result = {}
        for directed,label in ((False,"point"),(True,"directed")):
            ctx = iv if directed else c
            z,r = ctx.mpf(packet["Z"]),ctx.mpf(rho)
            de = self.leading.delta_iv if directed else self.leading.delta
            rows = packet[label]
            fields = {key:sum((rows[key][n][0]*r**n for n in range(1,N+1)),ctx.mpf(0))
                      for key in ("f","u","k","p")}
            mean = sum((rows["u"][n][0]*r**n/(n+1) for n in range(1,N+1)),ctx.mpf(0))
            meanZ = sum((rows["u"][n][1]*r**n/(n+1) for n in range(1,N+1)),ctx.mpf(0))
            fields["mean_u"] = mean
            # n1 has 1+delta here; the leading mean recovery has 1-delta.
            fields["q"] = (2*z*fields["u"]-(1+de)*z*mean-(1-z*z)*meanZ)/(1-de*z*z)
            result[label] = fields
        values = {}
        for field,point in result["point"].items():
            decimal = c.nstr(point,self.precision)
            box = result["directed"][field]
            error = max(abs(v) for v in _bounds(c,box-iv.mpf(decimal)))
            values[field] = dict(approximation=decimal,
                finite_prefix_absolute_error_upper=c.nstr(error*(1+c.mpf("1e-30")),40),
                finite_prefix_error_exact_mpf_tuple=list(error._mpf_),
                finite_prefix_source_enclosure_exact_mpf_tuples=[list(v) for v in box._mpi_],
                error_applies_to_serialized_decimal=True)
        return dict(Z=packet["Z"],rho=rho,radial_degree=N,normalized_prefix=values,
            normalization=dict(f="epsilon*F1/F0(Z)",u="Uz1/epsilon",k="K1/epsilon",
                p="P1/(epsilon*Pstar^2)",q="Q1/epsilon",mean_u="(Uz1+K1)/epsilon"),
            original_coefficient_enclosures=True,finite_radial_prefix_only=True,
            infinite_radial_remainder_bounded=False,full_original_n1_point_error_available=False,
            full_original_n1_point_oracle=False,temporal_n_ge_2_recursion=False)

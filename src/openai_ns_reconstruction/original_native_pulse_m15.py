"""Actual repaired weighted m15 on the original selected Rp-to-Rv pulse.

The inlet is the negative complete pulse integral, from the checked repaired
terminal identity. Main partial integrals are computed separately. End and
gap use finite remaining sources; no m13/P1 inlet is reused as m15.
"""
import copy
from fractions import Fraction
import gzip
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_native_pulse_pressure import OriginalNativePulsePressureVelocity


class OriginalNativePulseM15Numeric:
    def __init__(self, *, pressure_cells=128, **kwargs):
        self.velocity = OriginalNativePulsePressureVelocity(pressure_cells=pressure_cells, **kwargs)
        self.pressure = e = self.velocity.pressure
        self.native = n = e.native
        self.ctx, self.iv, self.J = e.ctx, e.iv, e.J
        self.mu, self.delta, self.family = e.mu, e.delta, e.family
        self.a0, self.a1, self.c = e.a, self.iv.mpf('1.5')-self.mu, 1-2*self.mu
        self.cells, self.directory = pressure_cells, n.directory
        self.hashes, self.code_hashes = dict(e.hashes), dict(e.code_hashes)
        self._k1, self._main_matrix, self._end_matrix, self._packets = {}, {}, {}, {}
        stem = PREFIX+'compliant_current_original_n1_pulse_numeric_bridge'
        p = self.directory/(stem+'_check.json')
        receipt = json.loads(p.read_text(encoding='utf-8'))
        if (not receipt['all_passed'] or not receipt['current_original_true_pulse_pressure_C1_numeric_bridge_installed']
                or receipt['source_family'] != self.family):
            raise ValueError('Checked complete actual pulse C1 source required')
        for name, digest in receipt['input_hashes'].items():
            if name in self.hashes:
                if self.hashes[name] != digest: raise ValueError('Conflicting m15 source: '+name)
            elif hashlib.sha256((self.directory/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Changed m15 source: '+name)
            self.hashes[name] = digest
        self.hashes[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        p = self.directory/(stem+'.json.gz'); raw = p.read_bytes()
        self.hashes[p.name] = hashlib.sha256(raw).hexdigest()
        data = json.loads(gzip.decompress(raw)); identity = data['source_identity']
        if (data['source_family'] != self.family or not identity['exact_main_gap_end_pressure_C1_graph_and_endpoints_replayed']
                or not identity['all_three_charts_and_native_mean_histories_retained']
                or not identity['Pstar_applied_once_and_signed_m15_once']):
            raise ValueError('Same original selected pulse, physical measures and finite histories required')
        self.complete = data['samples']
        self.source_bindings = dict(complete_pulse=stem+'.json.gz',
            actual_terminal='current_original_n1_repaired_own_C3_terminal_functions::terminal[m15]',
            inlet_identity='m15(Rp)=-complete_pulse_increment, since actual repaired m15(Rv)=0',
            signed_weight_already_in_source=True, independent_preRp_numeric_seam=False)
        name = 'src/openai_ns_reconstruction/original_native_pulse_m15.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def restore(self, row):
        return self.native.quiet.numeric._jet([row], True)[0]

    def positive(self, x): return self.pressure.positive(x)

    def box(self, left, right):
        return self.iv.mpf([_bounds(self.ctx, left)[0], _bounds(self.ctx, right)[1]])

    def units(self, packet, kind, extra=()):
        iv, n = self.iv, self.native
        rp = packet['coordinate']['Rp_global_y_log_terms']+[n._row(iv.log(110))]
        r2 = [n._row(self.restore(v)*2) for v in rp]
        lp = iv.mpf(n.amplitude.lp_iv); lm = n.quiet.numeric.outer.logmui
        saddle = iv.exp(-iv.mpf(2)/3*(iv.log(4)+lm))
        finite = iv.mpf(2)/3*iv.log(4)-iv.log(6)/2
        if kind == 'r2_linear': base = r2+[n._row(lp)]
        elif kind == 'r2_quadratic': base = r2+[n._row(2*lp)]
        elif kind == 'r1_linear': base = rp+[n._row(lp)]
        elif kind in ('end_linear','end_linear_gap_query','fixed_inlet'):
            inverse = iv.mpf('18.5')/self.mu if kind != 'end_linear_gap_query' else extra[0]
            if kind == 'end_linear_gap_query': extra = extra[1:]
            base = r2+n._rows(self.J(iv,[lp,inverse,-3*saddle,-iv.mpf(4)/3*lm,finite-13]))
            if kind == 'fixed_inlet':
                reference = iv.mpf('10.99')*iv.sqrt(2*iv.pi)*(iv.exp(5)-iv.exp(3))
                base += [n._row(iv.log(reference)),n._row(iv.log(n.constants['U']))]
        elif kind == 'end_quadratic':
            base = r2+n._rows(self.J(iv,[2*lp,11/self.mu,-6*saddle,-iv.mpf(8)/3*lm,2*finite-26]))
        elif kind == 'end_viscous':
            base = rp+n._rows(self.J(iv,[lp,iv.mpf('5.5')/self.mu,-3*saddle,-iv.mpf(4)/3*lm,finite-13]))
        else: raise ValueError('Unknown weighted m15 unit')
        return copy.deepcopy(base)+[n._row(v) for v in extra]

    def term(self, packet, name, value, kind, extra=()):
        return dict(name=name,coefficient_C0_C1=self.native._rows(value.truncate(1)),
            log_scale_terms=self.units(packet,kind,extra),actual_function_binding=self.source_bindings)

    def inlet(self, Z, packet):
        # A complete integral supplies the boundary value, not a partial cap.
        rows = self.complete['whole_Z']['normalized_true_pulse_C1_enclosures']['m15']
        for sample in self.complete.values():
            row = sample['Z']
            lo,hi = _bounds(self.ctx,self.native._decode(row))
            zl,zh = _bounds(self.ctx,self.iv.mpf(Z))
            if lo == hi and lo <= zl <= zh <= hi:
                rows = sample['normalized_true_pulse_C1_enclosures']['m15']; break
        value = -self.J(self.iv,[self.native._decode(v) for v in rows])
        return [self.term(packet,'actual_Rp_inlet_from_full_repaired_terminal',value,'fixed_inlet')]

    def kernel1(self, xi, t):
        e, iv, c = self.pressure, self.iv, self.ctx
        key = xi._mpi_,t._mpi_
        if key in self._k1: return self._k1[key]
        lo,hi = _bounds(c,t); length = min(lo,c.mpf(1000)); result = iv.mpf(0)
        if length:
            width = iv.mpf(length)/self.cells
            for j in range(self.cells):
                left,right = width*j,width*(j+1)
                g = e.gp(iv,xi-self.mu*self.box(left,right))['value']
                result += iv.exp(-self.a1*left)*self.native.decay_integral(iv,self.a1,width)*g
        if hi > length:
            edge = iv.mpf(length)
            tail = ((iv.mpf(hi)-edge)*iv.exp(-self.a1*edge) if hi <= 1000
                    else iv.exp(-self.a1*edge)/self.a1)
            result += iv.mpf([0,_bounds(c,11*tail)[1]])
        self._k1[key] = result
        return result

    def main_matrix(self, xi, t):
        """Positive local I1 basis matrix for (g,k0), with genuine partial tails."""
        key = xi._mpi_,t._mpi_
        if key in self._main_matrix: return self._main_matrix[key]
        iv,c,e = self.iv,self.ctx,self.pressure
        lo,hi = _bounds(c,t); length = min(lo,c.mpf(1000))
        matrix = [[iv.mpf(0) for _ in range(2)] for _ in range(2)]
        if length:
            width = iv.mpf(length)/self.cells
            for j in range(self.cells):
                left,right = width*j,width*(j+1); v = self.box(left,right)
                x = xi-self.mu*v; g = e.gp(iv,x)['value']
                error = 89*self.mu/(self.a0*self.a0)
                radius = _bounds(c,error)[1]
                k = self.positive(g/self.a0+iv.mpf([-radius,radius]))
                # Before the falling cutoff g is nondecreasing: k0<=g/a0.
                if _bounds(c,x)[1] <= 10:
                    upper = _bounds(c,g/self.a0)[1]
                    kl,kh = _bounds(c,k); k = iv.mpf([kl,min(kh,upper)])
                fields = (g,k)
                weight = iv.exp(-self.c*left)*self.native.decay_integral(iv,self.c,width)
                for a in range(2):
                    for b in range(a,2): matrix[a][b] += weight*fields[a]*fields[b]
        if hi > length:
            edge = iv.mpf(length)
            tail = ((iv.mpf(hi)-edge)*iv.exp(-self.c*edge) if hi <= 1000
                    else iv.exp(-self.c*edge)/self.c)
            bounds = (iv.mpf(11),11/self.a0)
            for a in range(2):
                for b in range(a,2): matrix[a][b] += iv.mpf([0,_bounds(c,tail*bounds[a]*bounds[b])[1]])
        matrix[1][0] = matrix[0][1]
        self._main_matrix[key] = matrix
        return matrix

    def main(self, Z, chart, coordinate, packet):
        iv,e,n = self.iv,self.pressure,self.native
        out = self.inlet(Z,packet)
        if Fraction(coordinate) == 0: return out
        selected,_,u,_,_ = n._data(Z); z,d,L,thetaZ = e.geometry(Z)
        ap = selected['selected_ap_Taylor'].truncate(4)
        bb = z*ap*2/L; cc = -(z*ap*(1-self.delta)+d*e.D(ap,thetaZ))/L
        Db,Dc,Dap = e.D(bb,thetaZ),e.D(cc,thetaZ),e.D(ap,thetaZ)
        t = iv.mpf(coordinate) if chart == 'entrance' else iv.mpf(coordinate)/self.mu
        xi = self.mu*t if chart == 'entrance' else iv.mpf(coordinate)
        k0,kv = e.kernel_pair(xi,t); k1 = self.kernel1(xi,t)
        shape = e.gp(iv,xi); g = shape['value']; A = bb*g+cc*k0; B = ap*g
        Iq = bb*k1+cc*(k0-k1); IqZ = Db*k1+Dc*(k0-k1)
        linear = ((A-Iq)/L+z*IqZ*(1-self.delta)/(2*L))/2
        matrix = self.main_matrix(xi,t); gg,gk,kk = matrix[0][0],matrix[0][1],matrix[1][1]
        QQ = bb*bb*gg+bb*cc*(2*gk)+cc*cc*kk
        UQz = ap*Db*gg+ap*Dc*gk; UU = ap*ap*gg; UUz = ap*Dap*gg
        quadratic_integral = (-QQ/2+d*UQz/L+z*z*UU*(2*(1-self.delta))/(L*L)+z*d*UUz*2/(L*L))/2
        quadratic_upper = A*A/4-z*z*B*B/(L*L)
        c0 = -(1+z*z)*n.constants['C1']
        cross = c0*(A-bb*k0-cc*kv)/2
        viscous = -(bb*(shape['derivative']*self.mu-self.pressure.b*g)+cc*(g-k0))
        out += [self.term(packet,'main_weighted_linear_generated',u*linear,'r2_linear',(t*iv.mpf('1.5'),-xi)),
                self.term(packet,'main_weighted_quadratic_upper',u*u*quadratic_upper,'r2_quadratic',(t,-2*xi)),
                self.term(packet,'main_weighted_quadratic_integral',u*u*quadratic_integral,'r2_quadratic',(t,-2*xi)),
                self.term(packet,'main_weighted_inherited_cross',u*u*cross,'r2_quadratic',(t/2,-xi)),
                self.term(packet,'main_weighted_finite_viscous',u*viscous,'r1_linear',(t/2,-xi))]
        # u^2*c0^2=(U*C1)^2 is exactly Z independent in the canonical inlet.
        coefficient = self.J.constant(iv,-(n.constants['U']*n.constants['C1'])**2/4,1)
        if chart == 'entrance': out.append(self.term(packet,'main_weighted_original_mean_log',coefficient*t,'r2_quadratic'))
        else: out.append(self.term(packet,'main_weighted_original_mean_log',coefficient*xi,'r2_quadratic',(-n.quiet.numeric.outer.logmui,)))
        return out

    def end_matrix(self, s):
        """All four beta/history basis products, including gaps and H_j*H_k."""
        s = Fraction(s)
        if s in self._end_matrix: return self._end_matrix[s]
        iv,e,c = self.iv,self.pressure,self.ctx
        cuts = sorted(set([s,Fraction(0)]+[v for v in map(Fraction,('-3.15','-2.85','-1.15','-.85')) if s<v<0]))
        cells = []
        for left,right in zip(cuts,cuts[1:]):
            width = (right-left)/self.cells
            for j in range(self.cells):
                l,r = left+j*width,left+(j+1)*width
                v = self.box(e.rat(l),e.rat(r)); ds = e.rat(width)
                beta = [e.raw_beta(iv,(v-center)/e.ell)/(e.ell*e.normal) for center in (-3,-1)]
                masses = [ds*iv.exp(self.a0*(v-center))*q for q,center in zip(beta,(-3,-1))]
                cells.append((v,ds,beta,masses))
        matrix = [[iv.mpf(0) for _ in range(4)] for _ in range(4)]
        future = [iv.mpf(0),iv.mpf(0)]
        for v,ds,beta,masses in reversed(cells):
            H = []
            for j,center in enumerate((-3,-1)):
                total = future[j]+masses[j]
                W = iv.mpf([_bounds(c,future[j])[0],_bounds(c,total)[1]])
                H.append(-iv.exp(self.a0*(center-v))*W);future[j] = total
            fields = beta+H; weight = ds*iv.exp(self.c*v)
            for a in range(4):
                for b in range(a,4): matrix[a][b] += weight*fields[a]*fields[b]
        for a in range(4):
            for b in range(a): matrix[a][b] = matrix[b][a]
        self._end_matrix[s] = matrix
        return matrix

    def end(self, Z, s, packet):
        iv,e,n = self.iv,self.pressure,self.native
        S = e.rat(s);selected,_,u,_,_ = n._data(Z)
        z,d,L,thetaZ = e.geometry(Z); C = [v.truncate(4) for v in selected['selected_scaled_end_coefficient_Taylor']]
        bb = [z*v*2/L for v in C];cc = [-(z*v*(1-self.delta)+d*e.D(v,thetaZ))/L for v in C]
        Db,Dc = [e.D(v,thetaZ) for v in bb],[e.D(v,thetaZ) for v in cc]
        data = [e.bump(s,j) for j in (-3,-1)]; zero = self.J.constant(iv,0,3)
        B = sum((v*q['beta'] for v,q in zip(C,data)),zero)
        A = sum((v*q['beta']+w*q['H'] for v,w,q in zip(bb,cc,data)),zero)
        Ay = sum((v*q['by']+w*(q['beta']-q['H']*self.a0) for v,w,q in zip(bb,cc,data)),zero)
        Iq,IqZ = zero,zero
        for j,center in enumerate((-3,-1)):
            W0 = data[j]['W']*iv.exp(self.a0*center)
            W1 = sum((row['ds']*iv.exp(self.a1*row['v'])*row['beta'] for row in data[j]['cells']),iv.mpf(0))
            H = -W1+iv.exp(S)*W0
            Iq += bb[j]*W1+cc[j]*H; IqZ += Db[j]*W1+Dc[j]*H
        matrix = self.end_matrix(s); coeff = bb+cc; Dcoeff = Db+Dc
        QQ = sum((coeff[a]*coeff[b]*matrix[a][b] for a in range(4) for b in range(4)),zero)
        UQz = sum((C[a]*Dcoeff[b]*matrix[a][b] for a in range(2) for b in range(4)),zero)
        UU = sum((C[a]*C[b]*matrix[a][b] for a in range(2) for b in range(2)),zero)
        UUz = sum((C[a]*e.D(C[b],thetaZ)*matrix[a][b] for a in range(2) for b in range(2)),zero)
        linear = (A*iv.exp(self.a1*S)/L+Iq/L-z*IqZ*(1-self.delta)/(2*L))/2
        w = iv.exp(self.c*S)
        quadratic = -(-A*A*(w/2)-QQ/2+d*UQz/L+z*z*B*B*(2*w)/(L*L)
                      +z*z*UU*(2*(1-self.delta))/(L*L)+z*d*UUz*2/(L*L))/2
        viscous = -(Ay-A*e.b)*iv.exp(self.a0*S)
        out = [self.term(packet,'end_weighted_terminal_linear',u*linear,'end_linear'),
               self.term(packet,'end_weighted_terminal_quadratic',u*u*quadratic,'end_quadratic'),
               self.term(packet,'end_weighted_finite_viscous',u*viscous,'end_viscous')]
        return out,A,e.D(A,thetaZ)

    def gap(self, Z, chart, coordinate, packet):
        iv,e,n = self.iv,self.pressure,self.native
        out,A,DA = self.end(Z,'-4',packet)
        if chart == 'gap_end' and Fraction(coordinate) == -4: return out
        _,_,u,_,_ = n._data(Z);z,d,L,thetaZ = e.geometry(Z)
        K,DK = A*iv.exp(-4*self.a0),DA*iv.exp(-4*self.a0)
        linear = u*z*DK*(1-self.delta)/(4*L)
        out.append(self.term(packet,'gap_weighted_linear_at_minus4',-linear*iv.exp(-4),'end_linear'))
        quadratic = u*u*K*K/4
        if chart == 'gap_end':
            S = iv.mpf(coordinate)
            out.append(self.term(packet,'gap_weighted_linear_at_query',linear,'end_linear',(S,)))
            out.append(self.term(packet,'gap_weighted_quadratic_log_history',quadratic*(-S-4),'end_quadratic'))
        else:
            xi = iv.mpf(coordinate); D = 13-xi
            out.append(self.term(packet,'gap_weighted_linear_at_query',linear,'end_linear_gap_query',((xi+iv.mpf('5.5'))/self.mu,)))
            out.append(self.term(packet,'gap_weighted_quadratic_inverse_mu',quadratic*D,'end_quadratic',(-n.quiet.numeric.outer.logmui,)))
            out.append(self.term(packet,'gap_weighted_quadratic_finite',quadratic*(-4),'end_quadratic'))
        return out

    def evaluate(self, Z, chart, coordinate):
        c,iv,n = self.ctx,self.iv,self.native
        if type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1:
            raise ValueError('Interior finite decimal Z required')
        if type(coordinate) is not str or not c.isfinite(c.mpf(coordinate)):
            raise ValueError('Finite decimal pulse coordinate required')
        x = c.mpf(coordinate)
        valid = ((chart=='entrance' and x>=0 and _bounds(c,self.mu*iv.mpf(coordinate))[1]<=11)
            or (chart=='main' and 0<=x<=11)
            or (chart=='gap' and 11<=x<13 and _bounds(c,(13-iv.mpf(coordinate))/self.mu)[0]>=4)
            or (chart=='gap_end' and x<=-4 and _bounds(c,iv.mpf(coordinate)+1/self.mu)[0]>=0)
            or (chart=='end' and -4<=x<=0))
        if not valid: raise ValueError('Original main/gap/end chart domain required')
        key = Z,chart,coordinate
        if key in self._packets: return self._packets[key]
        n._data(Z);inlet = n.inlets[Z]
        packet = dict(coordinate=dict(Rp_global_y_log_terms=copy.deepcopy(inlet['coordinate']['global_y_log_terms'])))
        with mp.workdps(iv.dps+30):
            if chart=='end': terms = self.end(Z,coordinate,packet)[0];mode='backward_actual_terminal'
            elif chart in ('gap','gap_end') or (chart=='main' and x==11):
                terms = self.gap(Z,'gap' if chart=='main' else chart,coordinate,packet);mode='backward_gap_and_actual_end'
            else: terms = self.main(Z,chart,coordinate,packet);mode='actual_terminal_derived_Rp_plus_main_prefix'
        out = dict(Z=Z,native_pulse_chart=chart,native_pulse_coordinate=coordinate,source_family=self.family,
            selected_integer=inlet['selected_integer'],coordinate=packet['coordinate'],m15_C0_C1_terms=terms,
            ordinary_Z_order=1,representation=mode,source_bindings=self.source_bindings,
            contract=dict(actual_repaired_terminal_m15_value_zero_consumed=True,actual_Rp_inlet_not_m13_P1=True,
                complete_pulse_only_used_for_actual_boundary_value=True,main_partial_prefix_integrated_separately=True,
                physical_weight_R_dR=True,weighted_quadratic_rate_one_minus_2mu=True,
                finite_R_Qy_endpoints_retained=True,nonzero_inherited_mean_log_retained=True,
                end_Q_squared_includes_history_products_off_bump_supports=True,gap_anchored_at_actual_minus4_end=True,
                Pstar_Rp_Rv_units_applied_once=True,signed_weight_not_negated_twice=True),
            limits=dict(weighted_m15_C0_C1=True,independent_preRp_numeric_seam=False,signed_point_centers_resolved=False,
                relative_pressure_cancellation_accuracy=False,global_physical_gradients=False,global_five_moment_numeric_join=False,
                common_global_N=False,exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key] = out
        return out

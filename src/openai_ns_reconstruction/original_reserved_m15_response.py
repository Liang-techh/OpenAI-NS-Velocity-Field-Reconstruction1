"""Source-correlated reserved swirl response, before coefficient boxing.

The selected three angular coefficients solve one fixed matrix system.
Contract the partial pressure row with that inverse first. At full support
the row is exactly the third matrix row, so its response is (0,0,1).
This restores the source correlation lost by adding three boxed beta terms.
It does not establish an independent numerical preRp pressure closure.
"""
import copy
from fractions import Fraction
import gzip
import hashlib
import json

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT, DIRECTORY, PREFIX, _bounds
from .original_compact_n1_pressure import weighted_beta_prefix
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets


def inverse3(matrix, determinant):
    """Cofactors of the original common matrix, with its positive determinant."""
    out=[]
    for i in range(3):
        row=[]
        for j in range(3):
            a,b=[k for k in range(3) if k!=j]
            c,d=[k for k in range(3) if k!=i]
            row.append((-1)**(i+j)*(matrix[a][c]*matrix[b][d]-matrix[a][d]*matrix[b][c])/determinant)
        out.append(row)
    return out


class OriginalReservedM15Response:
    def __init__(self, *, directory=DIRECTORY, precision=100, cells=128):
        if type(precision) is not int or precision<80 or type(cells) is not int or not 16<=cells<=1024:
            raise ValueError('Use precision>=80 and16..1024 original bump cells')
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30
        self.J=IntervalTaylor;self.beta_jets=beta_jets;self.cells=cells
        self.hashes={};self.code_hashes={};self._pressure_integrals={}
        self.operator=self.load('n1_reserved_linear_operator')
        self.coefficients=self.load('n1_reserved_coefficient_enclosures')
        self.full=self.load('n1_full_five_numeric_bridge')
        self.terminal=self.load('n1_terminal_integral_identities')
        terminal_receipt=self.receipts['n1_terminal_integral_identities']
        if not all(terminal_receipt['checks'][k]['passed'] for k in
                   ('actual_integrals','transported_partition','selected_inverse')):
            raise ValueError('Checked actual transported and selected inverse identities required')
        if not self.coefficients['source_identity']['passed'] or not self.full['source_identity']['passed']:
            raise ValueError('Original same-source selected action and signed full moments required')
        source=self.coefficients['source_identity']
        if not source['actual_selected_five_coefficient_graph_action_replayed'] or not source['actual_arc_source_identity']['passed']:
            raise ValueError('Actual selected angular inverse and original A_rc source identity required')
        self.arc_identity=copy.deepcopy(source['actual_arc_source_identity'])
        self.matrix_source=self.operator['exact_normalized_basis_and_matrix_functions']
        values=self.operator['actual_normalized_matrix_enclosures'];iv=self.iv
        self.mu=self.decode(values['mu']);self.ell=self.decode(values['ell'])
        self.centers=[self.decode(v) for v in values['centers']]
        self.normal=self.decode(values['base_normalizer'])
        self.matrix=[[self.decode(v) for v in row] for row in values['angular_rows']]
        det=self.decode(values['positive_angular_determinant'])
        if _bounds(self.ctx,det)[0]<=0 or _bounds(self.ctx,self.mu)[0]<=0:
            raise ValueError('Original positive determinant and nonzero mu required')
        self.inverse=inverse3(self.matrix,det)
        self.source_cells=self.coefficients['samples']['whole_Z']['cells']
        for name in ('original_reserved_m15_response.py','original_compact_n1_pressure.py',
                     'original_core_numeric.py'):
            name='src/openai_ns_reconstruction/'+name
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def merge(self, group, base, target):
        for name,digest in group.items():
            if target.get(name,digest)!=digest:
                raise ValueError('Conflicting original response source: '+name)
            if name not in target and hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed original response source: '+name)
            target[name]=digest

    def load(self, suffix):
        if not hasattr(self,'receipts'):self.receipts={}
        stem=PREFIX+'compliant_current_original_'+suffix
        path=self.directory/(stem+'.json.gz');raw=path.read_bytes()
        receipt_path=self.directory/(stem+'_check.json');receipt_raw=receipt_path.read_bytes()
        receipt=json.loads(receipt_raw);data=json.loads(gzip.decompress(raw))
        if not receipt['all_passed'] or data['source_family']!=receipt['source_family']:
            raise ValueError('Checked original response source required: '+suffix)
        if hasattr(self,'family') and data['source_family']!=self.family:
            raise ValueError('One original response family required')
        self.family=data['source_family'];self.merge(receipt['input_hashes'],self.directory,self.hashes)
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        self.hashes[receipt_path.name]=hashlib.sha256(receipt_raw).hexdigest()
        self.receipts[suffix]=receipt
        return data

    def decode(self,row):
        if 'source_enclosure_exact_mpf_tuples' in row:
            values=row['source_enclosure_exact_mpf_tuples']
        else:values=[row[k+'_exact_mpf_tuple'] for k in ('lower','upper')]
        return self.iv.mpf([self.ctx.make_mpf(tuple(v)) for v in values])

    def record(self,value):
        lo,hi=_bounds(self.ctx,value)
        return dict(lower=self.ctx.nstr(lo,self.ctx.dps),upper=self.ctx.nstr(hi,self.ctx.dps),
            source_enclosure_exact_mpf_tuples=[list(v) for v in value._mpi_])

    def intersect(self,a,b):
        lo,hi=_bounds(self.ctx,a);ll,hh=_bounds(self.ctx,b)
        if max(lo,ll)>min(hi,hh):raise ArithmeticError('Same-function response enclosures are disjoint')
        return self.iv.mpf([max(lo,ll),min(hi,hh)])

    def partial_row(self,x):
        """Forward/backward original positive bump integrals, physical ds once."""
        if type(x) is not str or not self.ctx.isfinite(self.ctx.mpf(x)) or not 1<=Fraction(x)<=2:
            raise ValueError('Original compact decimal x in[1,2] required')
        iv,c=self.iv,self.ctx;s=iv.log(iv.mpf(x));rate=iv.mpf('.5')-self.mu
        forward=[];backward=[]
        for j in (1,2,3):
            r=(s-self.centers[j])/self.ell;lo,hi=_bounds(c,r)
            factor=iv.exp(-self.mu*self.centers[j])/self.normal
            def integral(a,b,k):
                lower=weighted_beta_prefix(self,a,k);upper=weighted_beta_prefix(self,b,k)
                return iv.mpf([_bounds(c,lower)[0],_bounds(c,upper)[1]])*factor
            forward.append(iv.mpf(0) if hi<=-1 else integral(max(c.mpf(-1),lo),min(c.mpf(1),hi),rate))
            backward.append(iv.mpf(0) if lo>=1 else integral(max(c.mpf(-1),-hi),min(c.mpf(1),-lo),-rate))
        exact_empty=all(v._mpi_==iv.mpf(0)._mpi_ for v in forward)
        exact_full=all(v._mpi_==iv.mpf(0)._mpi_ for v in backward)
        # The row is the source matrix row k=-mu, not a numerically inferred
        # equality of complete integrals. Its exact inverse action is e_2.
        if exact_empty:response=[iv.mpf(0)]*3
        elif exact_full:response=[iv.mpf(0),iv.mpf(0),iv.mpf(1)]
        else:
            response=[]
            for k in range(3):
                f=sum((forward[j]*self.inverse[j][k] for j in range(3)),iv.mpf(0))
                b=iv.mpf(int(k==2))-sum((backward[j]*self.inverse[j][k] for j in range(3)),iv.mpf(0))
                response.append(self.intersect(f,b))
        return response,dict(exact_empty=exact_empty,exact_full=exact_full,
            partial_integrals=[self.record(v) for v in forward],suffix_integrals=[self.record(v) for v in backward])

    def evaluate(self,Z,x):
        iv,c=self.iv,self.ctx
        if type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1:
            raise ValueError('Original interior decimal Z required')
        response,diagnostic=self.partial_row(x);z=iv.mpf(Z);zl,zh=_bounds(c,z)
        cells=[p for p in self.source_cells if _bounds(c,self.decode(p['clipped_Z']))[0]<=zh
            and _bounds(c,self.decode(p['clipped_Z']))[1]>=zl]
        if not cells:raise ValueError('Original selected coefficient source cover required')
        # A_rc / exp(logSigma) is an original exact rational source, so the
        # same A cancels theta_2=p15/A before taking any independent boxes.
        constant=iv.exp(-iv.mpf(17)/10-iv.mpf(5)*self.mu/2)
        arc=[constant/(1+z*z),-2*z*constant/(1+z*z)**2]
        pressure=self.full['samples']['whole_Z']['normalized_full_pressure_moment_C1_enclosures']['m15']
        inputs=[None,None,[self.decode(v) for v in pressure]]
        for i in (0,1):
            pieces=[]
            for p in cells:
                theta=[self.decode(v) for v in p['normalized_theta_inputs'][i]]
                pieces.append([arc[0]*theta[0],arc[1]*theta[0]+arc[0]*theta[1]])
            inputs[i]=[iv.mpf([min(_bounds(c,v[k])[0] for v in pieces),max(_bounds(c,v[k])[1] for v in pieces)]) for k in (0,1)]
        result=[-sum((response[i]*inputs[i][k] for i in range(3)),iv.mpf(0)) for k in (0,1)]
        return dict(Z=Z,native_compact_x=x,source_family=self.family,
            coefficient_C0_C1=[self.record(v) for v in result],
            angular_inverse_response=[self.record(v) for v in response],
            amplitude_times_theta_C0_C1=[[self.record(v) for v in pair] for pair in inputs],
            original_source_bindings=dict(matrix=self.matrix_source,
                selected_coefficient_functions=self.coefficients['true_function_bindings'],
                full_defect_functions=self.full['true_function_bindings'],
                terminal_identity=self.terminal['terminal_identity_bindings']['reduced_patch']['m15']),
            diagnostics=diagnostic,contract=dict(original_physical_minus_R_theta0_Utheta1_ds=True,
                original_rate_half_minus_mu_and_three_disjoint_supports=True,
                inverse_action_before_independent_beta_enclosures=True,
                exact_original_A_theta2_quotient_cancelled_before_boxing=True,
                ordinary_Z1_product_rule_on_other_inputs=True,
                full_support_response_is_exact_matrix_row_identity=True,
                independent_complete_pulse_not_used_to_define_this_response=True,
                no_source_midpoint_or_huge_exponential_selected=True),
            limits=dict(actual_reserved_m15_response_C0_C1=True,
                signed_source_function_centers_resolved=False,independent_numeric_Rp_seam=False))

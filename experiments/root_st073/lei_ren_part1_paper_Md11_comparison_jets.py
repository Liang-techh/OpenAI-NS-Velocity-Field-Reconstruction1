"""New-source core loader; old comparison equations reused with fresh data.

Production loading requires all124 orders and the new pressure-majorant
admission. It never promotes a partial state or relabels an old core.
"""
import hashlib
import json
import operator
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_Md11_interval_core import STEM,PARAMETERS
from lei_ren_part1_paper_candidate_gauge_core import _unpack_fixed,_unpack_rows
from lei_ren_part1_paper_candidate_shared_inlet import radial_product
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
TAIL='lei_ren_part1_paper_Md11_core_tail_admission.json'


class Md11ComparisonJets(IntervalComparisonJets):
    def __init__(self,steps=4):
        self.steps=operator.index(steps)
        if self.steps<4:raise ValueError('need at least four comparison steps')
        raw=(HERE/(STEM+'_state.json')).read_bytes();self.state_hash=hashlib.sha256(raw).hexdigest()
        state=json.loads(raw);identity=state['identity']
        tail=json.loads((HERE/TAIL).read_bytes())
        for record in (state,tail):
            for name,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('new core/tail dependency changed:'+name)
        if identity['parameters']!=PARAMETERS or identity['new_schedule_sha256']!=tail['new_schedule_sha256']:
            raise ValueError('new pressure/core/tail identity mismatch')
        self.degree=operator.index(state['completed_radial_order'])
        if self.degree!=124 or identity['target_radial_degree']!=124 or not tail['target_met']:
            raise ValueError('production new-source loader requires124 orders and analytic tail gate')
        if not tail['analytic_fixed_point_exists_for_new_fixed_datum']:
            raise ValueError('new analytic source not admitted')
        self.ctx=MPIntervalContext();self.ctx.dps=identity['precision'];self.precision=self.ctx.dps
        self.center_family=identity['center'];self.diagnostic_partial=False
        self.new_schedule_sha256=identity['new_schedule_sha256']
        with mp.workdps(self.precision+40):
            c=self.ctx
            self.tail=dict(target_met=True,rows=[dict(scaled_radial_order=r['scaled_radial_order'],
                axial_order=r['axial_order'],Phi_tail=read_interval(c,r['Phi_tail']),
                Psi_tail=read_interval(c,r['Psi_tail'])) for r in tail['tail_rows']],
                maximum_normalized_mixed_C3_tail=read_interval(c,tail['maximum_normalized_mixed_C3_tail']))
            fixed=_unpack_fixed(c,state['fixed'])
            rows={name:_unpack_rows(c,state['rows'][name]) for name in ('A','Uz','P')}
            if any(len(rs)!=125 or any(len(row)!=128-n for n,row in enumerate(rs)) for rs in rows.values()):
                raise ValueError('radial/axial depth mismatch')
            self.lam=c.mpf(PARAMETERS['Lambda']);self.eps=1/self.lam;self.delta=c.mpf(PARAMETERS['delta'])
            self.z=IntervalTaylor.variable(c,c.mpf(self.center_family),3)
            jet=lambda row:IntervalTaylor(c,row[:4])
            self.phi=[jet(row)/self.lam**n for n,row in enumerate(rows['A'])]
            self.u=[jet(row)/self.lam**n for n,row in enumerate(rows['Uz'])]
            self.ell=jet(fixed['ell_Z_taylor']);self.S=jet(fixed['S_Z_taylor']);self.p0=jet(fixed['P0_Z_taylor'])
            self.F0=fixed['F0_interval'];self.hb=mp.mpf('.005')
            if endpoints(self.F0)[0]<=0:raise ValueError('amplitude positivity lost')
            self.initial_phi=self.radial(self.phi,c.mpf(4))
            if endpoints(self.initial_phi[0])[0]<=0:raise ValueError('fresh inlet positivity unresolved')
            self.phi_squared=radial_product(self.phi,self.phi)
            self.phi_u=radial_product(self.phi,self.u);self.u_squared=radial_product(self.u,self.u)
        self.cache={};self.radial_cache={}

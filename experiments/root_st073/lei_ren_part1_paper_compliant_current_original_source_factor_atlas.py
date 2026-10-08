"""Canonical original factors collected symbolically before logarithmic evaluation.

logdelta=-4logPstar-30 and logR=log110+10logCstar+10logPstar+coordinate.
One fixed-Z, same-family atlas owns context/bases/ledger for cross-cell sums.
No astronomical radius-log subtraction or native radius exponential occurs.
"""
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_two_chart_all_N_coefficients as points

base=points.base;point=points.point;prior=points.prior;ep=points.ep
HERE,PREFIX,sha=points.HERE,points.PREFIX,points.sha


class OriginalSourceFactorAtlas:
    order=('logPstar','selected_logCstar','log_original_L','reserved_zero','reserved_zero')
    def __init__(self,frame,*,Z,dps=260):
        if type(frame) is not point.source.OriginalO2SourceParameterFrame:
            raise TypeError('Original validated parameter frame required')
        if type(dps) is not int or dps<260:raise ValueError('At least260 directed digits required')
        self.frame=frame;self.family=frame.family;self.Z=point.pressure.exact_Z(Z)
        definitions=frame.definitions
        if s.simplify(definitions['log_delta']+4*definitions['logPstar']+30)!=0:
            raise ValueError('Exact original delta relation required')
        if s.expand(definitions['logRref']-s.log(110)-10*(definitions['logCstar']+definitions['logPstar']))!=0:
            raise ValueError('Exact original reference radius relation required')
        self.ctx=c=MPIntervalContext();c.dps=dps
        self.ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
            positive_function_root_intersections=0,directed_independent_log_rescalings=0)
        with mp.workdps(dps+40):
            logP=c.exp(40)+11
            logC=c.mpf(mp.mp.make_mpf(frame.selected_logCstar_mpf_tuple))
            self.log110=c.ln(110);iz=c.mpf(int(self.Z.p))/int(self.Z.q)
            temporary=(logP,logC,c.mpf(0),c.mpf(0),c.mpf(0))
            unit=prior.ScaledEnclosure(prior.FormalScale(temporary),1,self.ledger)
            delta_box=unit.bounded_exp(-4*logP-30)
            L=1-delta_box*iz*iz
            if ep(L)[0]<=0:raise ArithmeticError('Positive original L source lost')
            self.bases=(logP,logC,c.ln(L),c.mpf(0),c.mpf(0))
        self.hashes={**frame.hashes,Path(__file__).name:sha(Path(__file__).name)}
    def scalar(self,value):return prior.ScaledEnclosure(prior.FormalScale(self.bases),value,self.ledger)
    def rational(self,value):
        value=point.source.exact_rational(value)
        return self.ctx.mpf(int(value.p))/int(value.q)
    def copy_interval(self,value):
        if hasattr(value,'_mpi_'):return self.ctx.mpf(ep(value))
        return self.ctx.mpf(value)
    def scale(self,original_powers,coordinate,offset=0):
        if len(original_powers)!=5 or original_powers[3]!=0:
            raise ValueError('Original (logP,logdelta,logL,zero,logR) basis required')
        p,d,ell,unused,r=original_powers
        coordinate=self.copy_interval(coordinate)
        # Exact affine source identities are applied to powers and offsets.
        # No numeric difference of native radius or delta logarithms is used.
        powers=(p-4*d+10*r,10*r,ell,0,0)
        extra=self.copy_interval(offset)-30*d+r*(self.log110+coordinate)
        return prior.FormalScale(self.bases,powers,extra)
    def term(self,source_powers,coefficient,*,coordinate,offset=0):
        r,p,d,ell=source_powers
        return prior.ScaledEnclosure(self.scale((p,d,ell,0,r),coordinate,offset),
            self.copy_interval(coefficient),self.ledger)
    def parameter(self,name):
        if name=='Pstar':return self.term((0,1,0,0),1,coordinate=0)
        if name=='delta':return self.term((0,0,1,0),1,coordinate=0)
        if name=='mu':return self.term((0,-4,0,0),1,coordinate=0,offset=-self.ctx.ln(1000))
        if name=='Rref':return self.term((1,0,0,0),1,coordinate=0)
        raise ValueError('Unsupported original atlas parameter: '+str(name))
    def add(self,left,right):
        left.scale.pair(right.scale)
        if left.ledger is not self.ledger or right.ledger is not self.ledger:
            raise ValueError('Canonical atlas ledger required for addition')
        if left.zero:return right
        if right.zero:return left
        # Compare collected relative powers, never two astronomical absolute
        # logarithms whose common terms can hide the actual source ordering.
        relative=(right.scale-left.scale).evaluate()
        if ep(relative)[0]>0:
            left,right=right,left;relative=(right.scale-left.scale).evaluate()
        if ep(relative)[1]>2*self.ctx.dps*mp.log(10):
            return left+right  # Existing directed independent-range fallback.
        self.ledger['collected_relative_scale_additions']=self.ledger.get('collected_relative_scale_additions',0)+1
        return prior.ScaledEnclosure(left.scale,left.coefficient+right.coefficient*left.bounded_exp(relative),self.ledger)
    def sum(self,values):
        result=self.scalar(0)
        for value in values:result=self.add(result,value)
        return result
    def rebase_piece_value(self,dispatcher,piece,value):
        if type(dispatcher) is not points.TwoChartOriginalPointSources or dispatcher.issued_pieces.get(id(piece)) is not piece:
            raise ValueError('Issued actual point frame required')
        if piece.family!=self.family or piece.Z!=self.Z or piece.graph_sha256!=dispatcher.source_graph_sha256:
            raise ValueError('Same original family, fixed Z and source graph required')
        source_frame=(dispatcher.reference.owner if piece.chart=='Rh_reference' else dispatcher.O2.owner).inputs.frame
        if source_frame.selected_logCstar_mpf_tuple!=self.frame.selected_logCstar_mpf_tuple or any(
                source_frame.definitions[key]!=self.frame.definitions[key]
                for key in ('logPstar','log_delta','logRref','L')):
            raise ValueError('Identical selected original parameter definitions required')
        if type(value) is not prior.ScaledEnclosure or value.ctx is not piece.ctx:
            raise TypeError('Live original directed source value required')
        anchor=piece.values['all_N_original_E_C0'];anchor.scale.pair(value.scale)
        if value.ledger is not anchor.ledger:raise ValueError('Original point ledger required')
        return prior.ScaledEnclosure(self.scale(value.scale.powers,self.rational(piece.coordinate),value.scale.offset),
            self.copy_interval(value.coefficient),self.ledger)
    def record(self):
        return dict(source_family=self.family,original_Z_exact=str(self.Z),canonical_basis_order=self.order,
            defining_basis=self.bases,exact_delta_identity='logdelta=-4logPstar-30',
            exact_radius_identity='logR=log110+10logCstar+10logPstar+coordinate',
            source_power_map='(p,d,ell,0,r)->(p-4d+10r,10r,ell,0,0)',
            offset_map='offset-30d+r*(log110+coordinate)',
            original_L_definition='1-exp(-4logPstar-30)*Z²',
            tiny_positive_delta_retained_as_source_factor_and_directed_L_budget=True,
            same_context_basis_and_ledger_for_cross_cell_sums=True,
            astronomical_log_differences_or_radius_exponentials_used=False,
            ordinary_Z_rows_rebased_as_functions_no_second_L_derivative_applied=True)

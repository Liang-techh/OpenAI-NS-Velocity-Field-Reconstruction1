"""Original two-component stress cone, without interval division by swirl.

H0/H1/H2 are the original strict cone margins times positive source factors.
The tensor-completion diagonal is not a third cone component.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import gzip
import json
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_background_stress as background
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=background.HERE,background.PREFIX,background.sha
NAME=PREFIX+'current_original_Rp_stress_cone.json.gz'
RECEIPT=PREFIX+'current_original_Rp_stress_cone_check.json'
GATE='current_original_Rp_source_bound_signed_stress_cone_point_margins_installed'
POINT_CONE='current_original_Rp_admissible_stress_cone_at_supported_point_certified'
OPEN=background.OPEN


@lru_cache(maxsize=1)
def cone_operators():
    operators,_,_,_=background.source_operators()
    bt,Ttheta=operators['cylindrical_stress']['r_theta'][0]
    bz,Tz=operators['cylindrical_stress']['r_z'][0]
    if s.cancel(bt-bz)!=0 or bt!=-2-background.DELTA:
        raise ValueError('Original two cone components require one positive physical prefactor')
    Ut=background.F['Utheta'];B=2*background.dy(background.F['Uz']);a=2+2*background.MU
    H0=(a-2)*a*Ut**2+B**2
    H1=a*Ut*Ttheta-B*Tz
    Hperp=a*Ut*Tz+B*Ttheta
    H2=2*a*Ut**2*H1**2-H0*Hperp**2
    result=dict(a=[(0,a)],Utheta_positive=[(0,Ut)],H0_v_minus_2=[(0,H0)],
        H1_first_margin=[(0,H1)],H2_quadratic_margin=[(0,H2)])
    definitions=background.operator_definitions({'cone_margins':result})
    return result,definitions


def cone_equivalence():
    a,U,B,Ttheta,Tz=s.symbols('a U B Ttheta Tz',real=True)
    t=-B/(a*U);v=a+B**2/(a*U**2)
    H0=(a-2)*a*U**2+B**2
    H1=a*U*Ttheta-B*Tz;Hperp=a*U*Tz+B*Ttheta
    H2=2*a*U**2*H1**2-H0*Hperp**2
    G1=Ttheta+t*Tz;G2=2*G1**2-(v-2)*(Tz-t*Ttheta)**2
    checks={name:s.cancel(expression)==0 for name,expression in dict(
        H0_equals_a_U_squared_times_v_minus_2=H0-a*U**2*(v-2),
        H1_equals_a_U_times_G1=H1-a*U*G1,
        H2_equals_a_cubed_U_fourth_times_G2=H2-a**3*U**4*G2).items()}
    if not all(checks.values()):raise ArithmeticError('Original stress cone equivalence differs')
    return dict(identities=checks,paper_reference='OpenAI equations (4.20)-(4.23)',
        positive_factors_required=['a','Utheta','lambda**(-2-delta)'],
        source_mapping='F=Utheta/sqrt(2R); a=1-2*Utheta_y/Utheta=2+2mu; b_s=2*Uz_y/Utheta',
        admitted_components=['r_theta','r_z'],completed_diagonal_excluded=True,
        finite_point_margins_do_not_certify_regional_cone=True)


@dataclass(frozen=True,eq=False)
class OriginalStressConePoint:
    chart:str


class CurrentOriginalRpStressCone:
    @source_precision
    def __init__(self,before,require_checked=True):
        if type(before) is not background.CurrentOriginalRpBackgroundStress or not before.acceptance_loaded:
            raise ValueError('Accepted live original background stress owner required')
        self.before,self.product,self.graph,self.ctx=before,before.product,before.graph,before.ctx
        self.differential,self.family_record=before.differential,before.family_record
        self.operators,self.definitions=cone_operators()
        self._definitions=copy.deepcopy(self.definitions)
        self.equivalence=cone_equivalence();self._equivalence=copy.deepcopy(self.equivalence)
        self.hashes=dict(before.hashes)
        for name in (background.NAME,background.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self._fields={};self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or any(receipt[key] for key in OPEN) or \
                    receipt['source_family']!=self.family_record or receipt['original_cone_operator_definitions']!=self.definitions or \
                    receipt['original_cone_equivalence']!=self.equivalence:
                raise ValueError('Original stress cone receipt/source/scope differs')
            for name in (Path(__file__).name,Path(__file__).stem+'_check.py',NAME):
                if receipt['input_hashes'].get(name)!=sha(name):raise ValueError('Unbound original cone source '+name)
            background.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        checks=dict(accepted_original_background=self.before.acceptance_loaded,
            actual_original_two_component_cone_definitions=background.operator_definitions({'cone_margins':self.operators})==self._definitions,
            original_cone_definitions_unchanged=self.definitions==self._definitions,
            original_positive_factor_equivalence_unchanged=self.equivalence==self._equivalence,
            original_mu_source_function=self.before.before.mu is self.product.amplitude.mu,
            same_live_original_graph=self.graph is self.before.graph)
        if not all(checks.values()):raise ValueError('Original stress cone source differs: '+str(checks))
        return checks

    @source_precision
    def evaluate(self,field,relative_width_target='1/1000'):
        self.assert_graph();record=self.differential._require(field)
        if record['chart']!='pulse_exit':raise ValueError('This original cone adapter admits pulse_exit')
        target=Fraction(relative_width_target)
        if not 0<target<1:raise ValueError('Exact relative-width target required')
        source=self.product.before.source(record['delivery']);self.before._validate_source(source)
        Ut=source['log_radius_mixed_rows']['Utheta']['y0_Z0']
        if background.ends(Ut.coefficients[0])[0]<=0:
            raise ValueError('Original strictly positive swirl source required for cone equivalence')
        reader=self.product.reader(record['delivery'])
        if background.ends(reader.at(self.before.before.mu))[0]<=0:
            raise ValueError('Original strictly positive mu source required')
        _,_,request=self.product.amplitude.before._validate_delivery(record['delivery'])
        margins={name:self.before._enclose(parts,source,reader,request.forward_coordinates,target)
            for name,parts in self.operators.items()}
        value=OriginalStressConePoint(record['chart'])
        fingerprint=self.product.amplitude.before._fingerprint(background.box.report(source))
        self._fields[id(value)]=(value,field,source,fingerprint,margins,
            self.before._source_dag({'cone_margins':margins}))
        return value

    @source_precision
    def report(self,value):
        self.assert_graph();entry=self._fields.get(id(value))
        if type(value) is not OriginalStressConePoint or entry is None or entry[0] is not value or value.chart!=entry[1].chart:
            raise ValueError('Live original cone point issued by this owner required')
        record=self.differential._require(entry[1]);source=self.product.before.source(record['delivery'])
        self.before._validate_source(source)
        if source is not entry[2] or self.product.amplitude.before._fingerprint(background.box.report(source))!=entry[3] or \
                any(self.graph.nodes[node]!=data for node,data in entry[5].items()):
            raise ValueError('Original cone source products or defining operators changed')
        margins=copy.deepcopy(entry[4])
        positive=all(row.get('signed_log_value',{}).get('sign')==1 for row in margins.values())
        return dict(chart=value.chart,source_family=copy.deepcopy(self.family_record),cone_margins=margins,
            original_cone_equivalence=copy.deepcopy(self.equivalence),positive_mu_and_swirl_source_required=True,
            signed_margins_certify_original_two_component_cone_at_this_point=positive,
            completed_diagonal_excluded_from_cone=True,original_physical_common_prefactor_retained=True,
            unrestricted_physical_point_API=False,regional_cone_certified=False,
            **{GATE:self.acceptance_loaded,POINT_CONE:self.acceptance_loaded and positive},**dict.fromkeys(OPEN,False))


@source_precision
def run(before,fields):
    began=time.monotonic();owner=CurrentOriginalRpStressCone(before,require_checked=False)
    values={name:owner.evaluate(field) for name,field in fields.items()}
    views={name:owner.report(value) for name,value in values.items()}
    result=dict(source_family=owner.family_record,actual_original_cone_points=views,
        original_cone_operator_definitions=owner.definitions,original_cone_equivalence=owner.equivalence,
        source_assertions=owner.assert_graph(),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **{GATE:False,POINT_CONE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(background.correlated.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,values,views

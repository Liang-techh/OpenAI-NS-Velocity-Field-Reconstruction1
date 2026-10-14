"""Actual reserved n1 own histories, with C1 mean/pressure and C0 Q/V.

The selected coefficients are retained as functions. All five patch
histories are genuine signed integrals, and the incoming pressure includes
the original core datum and fixed leading-order forcing. Higher derivatives
and global terminal closure are deliberately separate requirements.
"""
import copy
import gzip
import json
from fractions import Fraction
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_reserved_coefficient_enclosures as coefficients
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions as full

parent=full.parent
source,controls=full.source,full.controls
HERE,PREFIX,sha=coefficients.HERE,coefficients.PREFIX,coefficients.sha
NAME=PREFIX+'current_original_n1_reserved_own_fields.json.gz'
RECEIPT=PREFIX+'current_original_n1_reserved_own_fields_check.json'
GATE='current_original_n1_reserved_own_C1_histories_and_C0_QV_installed'
OPEN=tuple(dict.fromkeys(coefficients.OPEN+(
    'n1_repaired_own_M_V_P_installed','n1_reserved_own_QV_C1_installed',
    'n1_repaired_own_terminal_five_function_zeros',
    'n1_repair_velocity_C3_functions_installed',
    'required_C3_C5_derivative_mixed_bounds')))
MOMENTS=full.MOMENTS


def integral_pair(g,pair,variable,lower,upper,**metadata):
    """Ordinary Z differentiation commutes with this fixed-coordinate integral."""
    g.ids([pair.value,pair.Z,variable,lower,upper])
    return source.C1Function(*(controls.integral(g,row,g.nodes[variable.node]['name'],lower,upper,
        endpoint_independent_of_Z=True,ordinary_slow_Z_derivative_order=k,
        actual_signed_repaired_n1_history=True,numerical_cap_not_defining_integrand=True,
        **metadata) for k,row in enumerate((pair.value,pair.Z))))


def pressure_prefix_in_band(g,true_owner,s,upper=None):
    """Partial original Rc-to-Rc*exp(s) pressure, preserving actual bodies.

    Only the upper integration limit changes. The native leading functions,
    physical dR/R*dR measures and finite viscous endpoints remain the same.
    This is the fixed N1p contribution, not a forcing rebuilt from n1 patches.
    """
    a=parent.C1Algebra(g)
    binding=true_owner.preRp.raw['defining_bindings']['compact']
    ref=lambda n:source.FunctionRef(g,n)
    restore=lambda ids:source.C1Function(*map(ref,ids))
    x=ref(binding['actual_C3_original_band_variable'])
    upper=upper if upper is not None else g.unary('exp',s)
    integrals={}
    for key,ids in binding['physical_integrals'].items():
        rows=[]
        for order,node in enumerate(ids):
            attrs=copy.deepcopy(g.nodes[node])
            if (attrs.pop('operation')!='definite_integral'
                or attrs['variable']!=g.nodes[x.node]['name']
                or attrs['lower']!=g.one.node
                or g.nodes[attrs['upper']]!=dict(operation='exact_rational',numerator=2,denominator=1)
                or attrs['ordinary_slow_Z_derivative_order']!=order):
                raise ValueError('Original physical band pressure integral changed')
            attrs['upper']=upper.node
            # Identical metadata lets the s=log(2) check recover the old
            # complete integrals exactly, without replacing their bodies.
            rows.append(g.node('definite_integral',**attrs))
        integrals[key]=source.C1Function(*rows)
    raw=true_owner.linear.reports[full.linear.quiet.band.NAME]['actual_C3_repair_band_functions']
    rows=raw['profiles_at_band_x']['radial_y_rows']
    mapping=true_owner.preRp.raw['defining_bindings']['active']['import_mapping']
    scales={node for node in mapping.values() if g.nodes[node]['operation']=='current_original_physical_source_scale'
        and g.nodes[node].get('name')=='Pstar' and g.nodes[node].get('source_family')==true_owner.family}
    if len(scales)!=1:raise ValueError('Unique actual imported physical Pstar source required')
    scale=ref(scales.pop())
    alg=parent.active.c3.C3Algebra(g)
    S=parent.active.c3.C3Function(scale,g.zero,g.zero,g.zero)
    jet=lambda ids:alg.mul(S,parent.active.c3.C3Function(*map(ref,ids)))
    U=jet(rows[0]['V']);Uy=jet(rows[1]['V'])
    mean=jet(raw['complete_histories']['m']);my=jet(raw['complete_history_y']['m'])
    delta=true_owner.cutoff.inner.parameters['delta']
    # This is the leading n=0 source, whose weight is 1-delta. The repaired
    # n1 recovery below has the distinct weight 1+delta.
    q=parent.active.partial.recover_Q_increment(g,U,mean,delta)
    qy=parent.active.partial.recover_Q_increment(g,Uy,my,delta)
    Rc=g.unary('exp',ref(true_owner.linear.background['logRc']))
    at=lambda p:a.at(source.C1Function(p.value,p.Z),x,upper)
    end=parent.endpoint_functions(g,at(U),at(q),at(qy),g.mul(Rc,upper))
    lower={key:restore(ids) for key,ids in binding['finite_endpoints'][0].items()}
    functions=parent.pressure_functions(g,delta,integrals,lower,end)
    return functions,dict(variable=x,upper=upper,integrals=integrals,lower=lower,upper_fields=end,
        original_binding=copy.deepcopy(binding),original_physical_scale=scale,
        fixed_leading_N1p_retained=True,
        original_physical_Pstar_applied_once=True)


def patch_histories(g,true_owner):
    """Integrate the physical Section 15 densities with dR=R ds."""
    a=parent.C1Algebra(g);patch=true_owner.repair_velocity
    s=patch['log_R_over_Rc'];R=patch['physical_radius']
    bg=true_owner.linear.background;W=true_owner.linear.functions
    A=source.C1Function(*(source.FunctionRef(g,i) for i in bg['Astar_ordinary_Z0_Z3'][:2]))
    power=g.unary('exp',g.neg(g.mul(g.add(g.constant('1/2'),W['mu']),
        g.add(g.unary('log',g.constant(2)),source.FunctionRef(g,bg['logRc']),s))))
    theta0=a.scale(power,A)
    sqrt2R=g.unary('sqrt',g.mul(g.constant(2),R))
    U,theta,F=patch['Uz1'],patch['Utheta1'],patch['F1']
    # Local leading Uz0 is exactly zero on all new supports. The fixed
    # N1p and nonlocal leading mean are retained in the unpatched prefix.
    per_dR=dict(m11=U,m12=a.scale(g.mul(g.constant(2),R),F),
        m13=a.scale(g.quotient(g.one,R,'positive reserved physical radius'),a.mul(theta0,theta)),
        m14=a.scale(sqrt2R,a.mul(theta0,U)),m15=a.neg(a.mul(theta0,theta)))
    per_ds={key:a.scale(R,pair) for key,pair in per_dR.items()}
    t=g.symbol('original_n1_reserved_own_history_t')
    bodies={key:a.at(pair,s,t) for key,pair in per_ds.items()}
    histories={key:integral_pair(g,pair,t,g.zero,s,physical_measure='dR=R ds',paper_moment=key,
        original_selected_alpha_beta_functions=True,history_not_reset_between_bumps=True,
        no_leading_N1p_subtraction=True) for key,pair in bodies.items()}
    return histories,dict(per_dR=per_dR,per_ds=per_ds,bodies=bodies,variable=t,
        leading_Utheta0=theta0,physical_F1=F)


def recovery_fields(g,patch,mass,pressure,delta,patch_pressure_y):
    """The actual n1 value recovery needs mean_Z, but not mean_ZZ."""
    a=parent.C1Algebra(g);R=patch['physical_radius']
    mean=a.div(mass,a.fixed(R),'positive nonaxis physical radius R')
    K=a.sub(mean,patch['Uz1'])
    z=g.symbol('Z');d=g.sub(g.one,g.mul(z,z));L=g.sub(g.one,g.mul(delta,z,z))
    Q=g.quotient(g.sub(g.sub(g.mul(g.constant(2),z,patch['Uz1'].value),
        g.mul(g.add(g.one,delta),z,mean.value)),g.mul(d,mean.Z)),L,
        'original positive L=1-delta*Z^2 on the common axial domain')
    return dict(**{key:patch[key] for key in ('Uz1','Utheta1','F1')},
        mean1=mean,K1=K,M1=mass,P1=pressure,Q1=Q,V1=g.mul(R,Q),
        mean1_y=a.sub(patch['Uz1'],mean),M1_y=a.scale(R,patch['Uz1']),
        P1_y_patch=patch_pressure_y)


class CurrentOriginalN1ReservedOwnFields:
    def __init__(self,coefficient_owner,require_checked=True):
        if (type(coefficient_owner) is not coefficients.CurrentOriginalN1ReservedCoefficientEnclosures
            or not coefficient_owner.acceptance_loaded):
            raise TypeError('Checked actual reserved coefficient provider required')
        coefficient_owner.assert_graph()
        self.coefficients=coefficient_owner;self.full=coefficient_owner.full.full
        self.g=g=coefficient_owner.g;self.family=copy.deepcopy(coefficient_owner.family)
        self.hashes=dict(coefficient_owner.hashes)
        for name in (coefficients.NAME,coefficients.RECEIPT):self._pin(name,sha(name))
        self.incoming_span=len(g.nodes)
        a=parent.C1Algebra(g);self.patch=self.full.repair_velocity
        self.s=self.patch['log_R_over_Rc'];self.R=self.patch['physical_radius']
        self.patch_histories,self.patch_bodies=patch_histories(g,self.full)
        self.base_pressure_prefix,self.base_pressure_bodies=pressure_prefix_in_band(g,self.full,self.s)
        self.incoming={key:self.full.moments[key] for key in ('m11','m12','m14')}
        for key in ('m13','m15'):
            self.incoming[key]=a.add(self.full.segments['inner'][key],self.full.lower.functions[key],
                self.full.preRp.pieces['active'][key])
        self.cumulative={key:a.add(self.incoming[key],self.patch_histories[key],
            self.base_pressure_prefix[key]) if key in ('m13','m15')
            else a.add(self.incoming[key],self.patch_histories[key]) for key in MOMENTS}
        delta=self.full.cutoff.inner.parameters['delta']
        self.fields=recovery_fields(g,self.patch,self.cumulative['m11'],self.cumulative['m13'],
            delta,self.patch_bodies['per_ds']['m13'])
        end=g.unary('log',g.constant(2))
        self.band_endpoint={key:a.at(value,self.s,end) for key,value in self.cumulative.items()}
        self.terminal={key:a.add(self.band_endpoint[key],self.full.preRp.pieces['quiet'][key],
            self.full.pulse.functions[key]) if key in ('m13','m15') else self.band_endpoint[key]
            for key in MOMENTS}
        self.definition=dict(actual_selected_coefficient_functions=True,
            reserved_band='Rc <= R <= 2Rc; s=log(R/Rc)',
            C1_mean_M_K_P_and_five_signed_cumulative_functions=True,
            Q_V_genuine_value_functions_using_mean_Z=True,Q_V_Z1_not_installed=True,
            core_first_unchanged_checked_cutoff_fields=True,
            pressure_core_axis_datum_and_lower_active_history_retained_once=True,
            fixed_leading_N1p_retained=True,own_swirl_pressure_patch_integrated=True,
            physical_Utheta1_over_sqrt_2R_once=True,
            leading_pure_power_used_only_under_disjoint_new_patch_supports=True,
            all_bump_histories_accumulated_without_resets=True,
            public_queries_restricted_to_exact_s_over_log2_in_closed_unit_interval=True,
            P1_y_patch_is_not_total_pressure_radial_derivative=True,
            terminal_integral_functions_defined=True,terminal_zeros_proved=False,
            global_own_fields_and_required_higher_derivative_bounds_complete=False)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or raw['own_field_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked true reserved own histories required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed reserved own source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed own-field source '+name)
        self.hashes[name]=digest

    def _snapshot(self):
        return parent.encode(dict(patch_histories=self.patch_histories,patch_bodies=self.patch_bodies,
            base_pressure_prefix=self.base_pressure_prefix,base_pressure_bodies=self.base_pressure_bodies,
            incoming=self.incoming,cumulative=self.cumulative,fields=self.fields,
            band_endpoint=self.band_endpoint,terminal=self.terminal,definition=self.definition))

    def assert_graph(self):
        self.coefficients.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot
            or self.hashes!=self._hashes):raise ValueError('Reserved own-field graph or history changed')


    def at_band_fraction(self,fraction):
        """Source-function query at s=fraction*log(2), confined to this chart.

        Returned handles denote actual functions of Z. They are not interval
        centers or numerical values, and Q/V have no claimed Z1 row yet.
        """
        self.assert_graph()
        if not isinstance(fraction,(int,str,Fraction)) or isinstance(fraction,bool):
            raise TypeError('An exact rational band fraction is required')
        q=Fraction(fraction)
        if not 0<=q<=1:raise ValueError('Own-field chart requires 0 <= s/log(2) <= 1')
        g=self.g;a=parent.C1Algebra(g);value=g.mul(g.constant(q),g.unary('log',g.constant(2)))
        fields={key:a.at(pair,self.s,value) if type(pair) is source.C1Function else
            g.node('function_substitution',expression=pair.node,variable=self.s.node,
                value=value.node,Z_independent_substitution=True)
            for key,pair in self.fields.items()}
        return dict(exact_s_over_log2=str(q),fields=fields,
            C1_mean_M_K_P=True,Q_V_value_only=True,total_pressure_radial_derivative_not_installed=True)


def run(coefficient_owner):
    began=time.monotonic();owner=CurrentOriginalN1ReservedOwnFields(coefficient_owner,require_checked=False)
    data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,own_field_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys((GATE,)+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_RESERVED_OWN_HISTORIES_DEFINED',len(owner.g.nodes),flush=True)
    return owner,data


if __name__=='__main__':
    raise SystemExit('run(coefficient_owner) reuses the checked live source owners')

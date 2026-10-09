"""Exact weighted dx/x references and actual Rm terminal integral replay."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from fractions import Fraction
import lei_ren_part1_paper_compliant_current_original_Rm_terminal_density_integrals as current
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def finite_fixture():
    c=MPIntervalContext();c.dps=200;p=mp.mp.clone();p.dps=240
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    P0=[f.scalar(1),f.scalar('.2')];Rm=f.scalar(7)
    partition=((71,40),(2,1),'Rh');comparisons=0
    def pc(value):return p.e if value=='Rh' else p.mpf(value[0])/value[1]
    def reference(left,right,rate,Z,sign,derivative=False):
        yl,yr=p.log(pc(left)),p.log(pc(right));w=yr-yl;suffix=1-yr;r=p.mpf(rate.numerator)/rate.denominator
        mass=w if not rate else -p.expm1(-r*w)/r
        first_u=w*w/2 if not rate else (1-(1+r*w)*p.exp(-r*w))/(r*r)
        first_s=w*mass-first_u
        constant=3+2*yl if derivative else 2+yl+Z*(3+2*yl)
        slope=2 if derivative else 1+2*Z
        return sign*p.exp(-r*suffix)*(constant*mass+slope*first_s)
    for Z in (p.mpf('.3'),p.mpf('-.4')):
        for sign in (-1,1):
            packets=[]
            for left,right in zip(partition,partition[1:]):
                xl,xr=current.coordinate(c,left),current.coordinate(c,right)
                y=c.mpf([ep(c.ln(xl))[0],ep(c.ln(xr))[1]])
                values={key:f.scalar((2+y+c.mpf(str(Z))*(3+2*y))*sign) for key in current.RATES}
                jets={key:f.scalar((3+2*y)*sign) for key in current.RATES}
                source=dict(source_family='finite diagnostic',actual_closed_radial_source_cell=True,actual_original_spatial_Z_density_interface_installed=True,
                    actual_original_Rm_radius_phase=dict(source_geometry=current.exact_geometry(left,right),candidate_N=257,
                        exact_source_Rm_factor=Rm,
                        actual_Rm_radius_phase_Z_independent=True,source_radial_logarithmic_cell_width=c.ln(xr/xl)),
                    actual_source_bound_phase_density_cells=[dict(source_family='finite diagnostic',source_geometry=current.exact_geometry(left,right),
                        original_common_P0_axial5=P0,exact_same_shared_radius_factor=Rm*c.mpf([ep(xl)[0],ep(xr)[1]]),
                        actual_original_five_signed_density_C0_Z=dict(kernels=values,Z_derivatives=jets))])
                packet=current.integrate_source_cell(f,source,left,right,'Rh',source_family='finite diagnostic',P0=P0,Rm_factor=Rm);packets.append(packet)
                for key,rate in current.RATES.items():
                    for n,row in enumerate(packet['actual_cell_to_target_integral_C0_Z'][key]):
                        expected=reference(left,right,rate,Z,sign,n==1)
                        lo,hi=ep(row.finite_interval());assert lo<=expected<=hi,(Z,sign,key,n);comparisons+=1
            total={key:[sum((packet['actual_cell_to_target_integral_C0_Z'][key][n] for packet in packets),f.scalar(0)) for n in range(2)] for key in current.RATES}
            width=1-c.ln(current.coordinate(c,partition[0]));incoming={key:[f.scalar((i+1)*sign),f.scalar((i+2)*sign)] for i,key in enumerate(current.RATES)}
            carried=current.affine_transport(f,incoming,total,width)
            for i,(key,rate) in enumerate(current.RATES.items()):
                for n,row in enumerate(carried[key]):
                    r=p.mpf(rate.numerator)/rate.denominator
                    expected=sum(reference(a,b,rate,Z,sign,n==1) for a,b in zip(partition,partition[1:]))
                    expected+=(i+1+n)*sign*p.exp(-r*(1-p.log(pc(partition[0]))))
                    lo,hi=ep(row.finite_interval());assert lo<=expected<=hi,(Z,sign,key,n,'incoming');comparisons+=1
            assert ep(carried['p'][0].coefficient)[0]!=0
    tiny=c.mpf(p.make_mpf((0,1,-10000,1)))
    mass=current.local.positive_kernel_mass(c,tiny,Fraction(3,2));assert ep(mass)[0]>0
    assert ep(current.local.positive_kernel_mass(c,c.mpf('.3'),0))==ep(c.mpf('.3'))
    rejected=0
    for partition in (('Rh',),((2,1),(71,40),'Rh'),((1,1),'Rh'),((71,40),(2,1))):
        try:current.validate_partition(c,partition)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid terminal source partition accepted')
    foreign=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    bad={key:[foreign.scalar(1),foreign.scalar(1)] for key in current.RATES}
    try:current.affine_transport(f,bad,total,width)
    except ValueError:rejected+=1
    else:raise AssertionError('Foreign incoming source basis accepted')
    for changed in ({**source,'source_family':'other family'},{**source,'actual_closed_radial_source_cell':False},
        {**source,'actual_source_bound_phase_density_cells':[{**source['actual_source_bound_phase_density_cells'][0],'original_common_P0_axial5':list(P0)}]}):
        try:current.integrate_source_cell(f,changed,left,right,'Rh',source_family='finite diagnostic',P0=P0,Rm_factor=Rm)
        except ValueError:rejected+=1
        else:raise AssertionError('Foreign source family, point cover or P0 object accepted')
    try:current.affine_transport(f,incoming,total,-1)
    except ValueError:rejected+=1
    else:raise AssertionError('Negative transport width accepted')
    return dict(passed=True,independent_signed_dx_over_x_cell_suffix_and_nonzero_incoming_comparisons=comparisons,
        exact_rate_zero_pressure_memory_checked=True,microscopic_positive_mass_not_dropped=True,
        invalid_partition_and_foreign_incoming_rejections=rejected,finite_reference_fixtures_only=True)


def native(owner,report):
    rows=cell_rows=masses=memory=0;signs={}
    for label in ('0','.5'):
        packet=owner.contribution(label)
        assert current.base.encoded(fields.serialized(packet))==report['frames'][label]
        op=owner.upstream.upstream.upstream.upstream.owner(label).op;f=op.flow
        assert packet['exact_common_P0_axial5'] is op.P0
        assert packet['finite_N_incoming_defect_history_is_unsupplied_affine_argument']
        assert packet['complete_prefix_or_corrected_Rh_history_not_claimed']
        assert packet['original_pressure_incoming_decay_exactly_one']
        assert ep(packet['original_incoming_to_Rh_decays']['p'])==(1,1)
        assert len(packet['actual_cells_to_Rh'])==4
        for left,right,cell in zip(current.PARTITION,current.PARTITION[1:],packet['actual_cells_to_Rh']):
            assert cell['source_geometry']==current.exact_geometry(left,right)
            assert cell['actual_phase_source']['exact_source_Rm_factor'] is op.Rm_factor
            assert cell['actual_phase_source']['actual_Rm_radius_phase_Z_independent']
            assert cell['no_extra_R_Pstar_or_N_factor'] and cell['source_cells_and_phase_unions_not_samples']
            assert ep(cell['source_radial_logarithmic_cell_width'])[0]>0
            assert ep(cell['source_logarithmic_suffix_to_target'])[0]>=0
            for key in current.RATES:
                assert ep(cell['original_positive_own_rate_masses'][key])[0]>0;masses+=1
                assert ep(cell['original_own_rate_suffix_decays'][key])[0]>0
                if key=='p':assert ep(cell['original_own_rate_suffix_decays'][key])==(1,1)
                for row in cell['actual_cell_to_target_integral_C0_Z'][key]:
                    assert row.ctx is f.c and row.scale.bases is f.logs and row.ledger is f.ledger;cell_rows+=1
            assert all(cell[key] is False for key in fields.previous.OPEN)
        for key,values in packet['actual_terminal_local_defect_integral_C0_Z'].items():
            assert len(values)==2
            for row in values:
                assert row.ctx is f.c and row.scale.bases is f.logs and row.ledger is f.ledger;rows+=1
            for n in (0,1):
                added=sum((cell['actual_cell_to_target_integral_C0_Z'][key][n] for cell in packet['actual_cells_to_Rh']),f.scalar(0))
                assert current.base.encoded(fields.serialized(added))==current.base.encoded(fields.serialized(values[n]))
        for source in (packet['actual_leading_inlet_memory'],packet['actual_leading_Rh_memory']):
            assert set(source)==set(current.RATES)
            assert all(len(values)==6 and all(v.scale.bases is f.logs and v.ledger is f.ledger for v in values) for values in source.values());memory+=1
        assert all(packet[key] is False for key in fields.previous.OPEN)
        signs[label]={key:[row.record()['sign'] for row in values] for key,values in packet['actual_terminal_local_defect_integral_C0_Z'].items()}
        try:owner.transport_supplied_incoming(label,None)
        except ValueError:pass
        else:raise AssertionError('Unspecified finite-N incoming history reset')
    return dict(passed=True,actual_whole_terminal_local_C0_Z_integral_rows=rows,
        actual_source_cell_to_Rh_C0_Z_rows=cell_rows,positive_original_own_rate_cell_masses=masses,
        actual_leading_inlet_and_Rh_five_history_memories_retained=memory,
        conditional_integral_enclosure_signs=signs,actual_radial_cells=8,
        true_dx_over_x_measure_and_cell_suffixes_checked=True,
        unspecified_incoming_history_rejected_at_both_frames=True,
        complete_prefix_Rc_repair_whole_axis_and_global_N_not_admitted=True)


def run():
    began=time.monotonic()
    with mp.workdps(260):fixture=finite_fixture()
    print('Independent signed dx/x integration, suffixes and nonzero inlet transport PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmTerminalDensityIntegrals(require_checked=False);evidence=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_exact_signed_integral_and_incoming_fixture=fixture,actual_live_source=evidence,
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rm terminal five-density C0/Z Duhamel integral enclosures PASS',flush=True);return result


if __name__=='__main__':run()

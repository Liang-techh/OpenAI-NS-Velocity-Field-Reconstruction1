"""Focused acceptance of the two-source Rsh physical query adapter."""
import gzip
import json
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_Rsh_physical_velocity_interface as source


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed Rsh physical source: '+name)
    field=source.CurrentRshPhysicalVelocityInterface(require_checked=False)
    if data['Rsh_physical_velocity_interface_definition_sha256']!=field.definition:
        raise ValueError('Foreign Rsh source/unit definition')
    if source.numeric.encode(source.numeric.pack(field.theorem))!=data['exact_physical_adapter_theorem']:
        raise ValueError('Rsh physical source linearity theorem differs')
    views=json.loads(gzip.decompress((source.HERE/source.VIEWS_NAME).read_bytes()))
    expected=dict(whole_Z=field.interface(),fresh_time_nu=field.interface(Z=('.1','.4'),log_tau=(-12,-11),viscosity='.8'),
        actual_physical_point=field.physical_interface('.237','-2.337','.337','1.3'))
    if source.numeric.encode(source.numeric.pack(expected))!=views:
        raise ValueError('Rsh saved independent source trace views differ')
    count=0
    for name in ('whole_Z','fresh_time_nu'):
        view=expected[name];common=view['common_physical_velocity_pressure_rows']
        if len(common)!=216 or view['all_global_velocity_interfaces_certified']:
            raise ValueError('Incomplete or overclaimed physical Rsh interface')
        for side in ('left','right'):
            trace=view[side+'_complete_velocity_pressure_source']
            if trace['named_source_side']!=side or trace['primitive_source']!=source.PACKETS[side]+'#/actual_Rsh_exit':
                raise ValueError('Distinct actual source germs lost')
            for key,row in view[side+'_physical_contribution_rows'].items():
                shared=common[key]
                if row['exact_zero']:continue
                if shared['exact_zero'] or source.t.endpoints(shared['log_absolute_upper'])[1]<source.t.endpoints(row['log_absolute_upper'])[1]:
                    raise ArithmeticError('Common Rsh bound fails to cover a proved source side')
                count+=1
    point=expected['actual_physical_point'];c=field.ctx;z=point['physical_z_source_coefficient']
    q=point['actual_log_lambda'];nu=c.mpf('1.3');lam=c.exp(q)
    physical_z=z*c.exp(point['physical_z_log_scale'])
    tau=lam**2-physical_z**2/nu*c.exp(2*field.data['source']['delta']*q)
    requested=c.exp(c.mpf('-2.337'));tl,th=source.t.endpoints(tau);rl,rh=source.t.endpoints(requested)
    if not tl<=rl<=rh<=th or source.t.endpoints(q-c.mpf('-2.337')/2)[0]<=0:
        raise ArithmeticError('Rsh physical query lost the actual implicit lambda equation')
    # Independently keep radial coordinate logarithmic: its exponent is
    # astronomical, while this exact coordinate identity stays computable.
    radial=point['physical_log_r']-(c.ln(nu)+c.ln(2)+field.data['logR']+2*q)/2
    if not source.t.endpoints(radial)[0]<=0<=source.t.endpoints(radial)[1]:
        raise ArithmeticError('Actual Rsh radius map differs')
    invalid=(lambda:field.interface('reshape'),lambda:field.interface(Z=(-2,0)),
        lambda:field.interface(log_tau=mp.inf),lambda:field.interface(viscosity=0),
        lambda:field.interface(theta=mp.nan),lambda:field.physical_interface(1,-2),
        lambda:field.physical_interface(-1,-2),lambda:field.physical_interface('.2',-2,viscosity=-1))
    for request in invalid:
        try:request()
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid Rsh physical request accepted')
    hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)}
    result=dict(all_passed=True,source_family=data['source_family'],
        Rsh_physical_velocity_interface_definition_sha256=field.definition,
        physical_source_map_and_unit_identities=len(field.theorem['identities']),
        common_Rsh_physical_spatial4_time1_contribution_count=216,
        independent_side_nonzero_common_upper_comparisons=count,
        fresh_time_nu_sector_and_actual_lambda_queries_regenerated=True,
        independent_signed_axial_implicit_equation_and_log_radial_map_enclosed=True,
        two_actual_frozen_source_normalizations_retained=True,invalid_requests_rejected=len(invalid),
        current_Rsh_physical_spatial4_time1_interface_available=True,
        all_global_velocity_interfaces_energy_cones_recursion_NS_certified=False,input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Rsh physical interface PASS:',len(field.theorem['identities']),'source/map identities;',count,'two-sided bounds;8 guards',flush=True)
    return result


if __name__=='__main__':run()

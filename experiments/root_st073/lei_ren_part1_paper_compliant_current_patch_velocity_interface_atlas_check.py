"""Exercise named patch exports, new sectors, actual lambda and nu guards."""
import gzip
import json
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_patch_velocity_interface_atlas as source


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest: raise ValueError('Changed named patch interface input: '+name)
    field=source.CurrentPatchVelocityInterfaceAtlas(require_checked=False)
    expected=field.manifest();expected['input_hashes'][source.VIEWS_NAME]=source.sha(source.VIEWS_NAME)
    if source.numeric.encode(source.numeric.pack(expected))!=data:
        raise ValueError('Current named patch interface manifest differs')
    views=json.loads(gzip.decompress((source.HERE/source.VIEWS_NAME).read_bytes()))
    counts={}
    for name in source.SEAMS:
        view=views[name];rows=view['common_physical_velocity_pressure_rows']
        if len(rows)!=216 or view['current_common_physical_contribution_count']!=216:
            raise ValueError('Incomplete common spatial4/time1 interface layout')
        if view['left_defining_source_germ']['named_source_side']!='left' or view['right_defining_source_germ']['named_source_side']!='right':
            raise ValueError('Independent source germs were lost')
        if view['patch_velocity_interface_atlas_definition_sha256']!=field.definition or view['all_global_velocity_interfaces_certified']:
            raise ValueError('Foreign or overclaimed patch provider')
        counts[name]=len(rows)
    fresh=field.interface('patch_support_59',Z=('-0.7','0.8'),log_tau=('-12','-11'),viscosity='.8')
    if source.numeric.encode(source.numeric.pack(fresh))!=views['fresh_time_nu_sector']:
        raise ValueError('Fresh arbitrary time/nu sector differs')
    physical=field.physical_interface('patch_support_71','.237','-2.337','.337','1.3')
    if source.numeric.encode(source.numeric.pack(physical))!=views['actual_physical_point']:
        raise ValueError('Actual lambda query differs')
    # Independent exact-coordinate source: delta is retained, including
    # the signed Z factor in the axial coordinate. A directed enclosure
    # of the implicit equation must include the exact requested tau.
    c=field.ctx;z=c.mpf('.237');lt=c.mpf('-2.337');nu=c.mpf('1.3')
    q=physical['physical_coordinates']['actual_log_lambda']
    coordinates=physical['physical_coordinates']
    lam=c.exp(q)
    actual_z=coordinates['physical_z_source_coefficient']*c.exp(coordinates['physical_z_log_scale'])
    tau_from_physical=lam**2-(actual_z**2/nu)*c.exp(2*field.source['delta']*q)
    tau=c.exp(lt)
    if not source.t.endpoints(tau_from_physical)[0]<=source.t.endpoints(tau)[0]<=source.t.endpoints(tau)[1]<=source.t.endpoints(tau_from_physical)[1]:
        raise ArithmeticError('Independent physical coordinates do not satisfy the actual implicit map')
    if source.t.endpoints(q-lt/2)[0]<=0:
        raise ArithmeticError('Nonzero axial location was replaced by lambda=sqrt(tau)')
    bound=physical['actual_physical_velocity_pressure_bounds']
    for i,j,b in source.t.INDICES:
        index='x%d_y%d_z%d'%(i,j,b)
        for component,rows in bound['actual_lambda_viscosity_spatial4_bounds'][index].items():
            expected_power=c.mpf(str(source.viscosity_power(component,i+j+b)))
            for row in rows.values():
                if row['constant_viscosity_power']._mpi_!=expected_power._mpi_:
                    raise ArithmeticError('Wrong spatial viscosity power')
    invalid=(lambda:field.interface('patch_support_50'),
        lambda:field.interface('patch_support_49',Z=('-1.1',0)),
        lambda:field.interface('patch_support_49',viscosity=0),
        lambda:field.interface('patch_support_49',viscosity=-1),
        lambda:field.interface('patch_support_49',viscosity=mp.inf),
        lambda:field.interface('patch_support_49',log_tau=mp.inf),
        lambda:field.interface('patch_support_49',theta=mp.nan),
        lambda:field.physical_interface('patch_support_49',1,-2),
        lambda:field.physical_interface('patch_support_49',(-1,0),-2))
    for query in invalid:
        try:query()
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Unsupported or nonfinite physical source request accepted')
    hashes=dict(data['input_hashes']);hashes[source.NAME]=source.sha(source.NAME);hashes[Path(__file__).name]=source.sha(Path(__file__).name)
    receipt=dict(all_passed=True,patch_velocity_interface_atlas_definition_sha256=field.definition,
        exact_registry_coordinate_and_viscosity_identities=len(field.theorem['identities']),
        registry_named_common_trace_group_counts=counts,
        new_time_viscosity_sector_and_actual_lambda_query_recomputed=True,
        independent_implicit_physical_coordinate_equation_enclosed=True,
        finite_physical_requests_at_Z_endpoints_rejected=True,invalid_requests_rejected=len(invalid),
        current_patch_physical_velocity_interface_provider_available=True,
        all_global_velocity_interfaces_certified=False,
        energy_common_N_cones_recursion_corrected_NS_certified=False,input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Named patch velocity atlas PASS:6x216 common groups; actual lambda/new time/nu queries;9 guards',flush=True)
    return receipt


if __name__=='__main__':run()

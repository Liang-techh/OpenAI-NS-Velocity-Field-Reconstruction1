"""Inspect/restore seven actual Rm and outer source interfaces."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_Rm_outer_mixed_source_joins as current


class RestoredRmOuterJoinArchive:
    """Archived endpoint source views, without numerical or ancestor owners."""
    def __init__(self,raw):
        self.g=g=current.source.FunctionTransportGraph();g.nodes=raw['exact_graph_nodes']
        g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
        assert len(g.keys)==len(g.nodes)
        assert g.nodes[g.zero.node]==dict(operation='exact_rational',numerator=0,denominator=1)
        assert g.nodes[g.one.node]==dict(operation='exact_rational',numerator=1,denominator=1)
        self.bindings=raw['Rm_outer_mixed_source_join_bindings']
        self.numerical_oracle_available=False
    def endpoint(self,interface,side):
        if interface not in {key for key,*_ in current.JOINS} or side not in ('left','right'):
            raise ValueError('Actual original interface and side required')
        packet=self.bindings['joins'][interface][side];g=self.g
        def ref(node):
            assert 0<=node<len(g.nodes)
            return current.source.FunctionRef(g,node)
        jet=lambda ids:current.outer.FunctionJet(tuple(map(ref,ids)))
        jets=lambda values:list(map(jet,values))
        return dict(native_query=ref(packet['native_query']),radius=jet(packet['radius']),
            physical={key:jets(value) for key,value in packet['physical'].items()},N_y=jets(packet['N_y']),
            own={group:{field:jets(value) for field,value in fields.items()} for group,fields in packet['own'].items()},
            numerical_oracle_available=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    args=parser.parse_args()
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    assert rec['all_passed'] and len(rec['checks'])==6 and all(proof['passed'] for proof in rec['checks'].values())
    assert all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert rec['actual_own_mixed_rows']==1722 and rec['actual_Rm_outer_interfaces']==[key for key,*_ in current.JOINS]
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed source',name)
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    assert len(raw['exact_graph_nodes'])==rec['exact_graph_nodes']
    queries=0
    if args.restore_source:
        archive=RestoredRmOuterJoinArchive(raw)
        for key,*_ in current.JOINS:
            for side in ('left','right'):
                packet=archive.endpoint(key,side)
                assert not packet['numerical_oracle_available']
                assert len(packet['physical']['U'])==len(packet['physical']['Q'])==5
                assert len(packet['N_y'])==3 and len(packet['own']['fields']['P1'])==4
                count=0
                for spec in archive.bindings['contract']:
                    values=packet['own'][spec['group']][spec['field']]
                    assert len(values)==spec['physical_y_order']+1
                    assert all(len(current.outer.rows(jet))==spec['ordinary_Z_order']+1 for jet in values)
                    count+=sum(len(current.outer.rows(jet)) for jet in values)
                assert count==246;queries+=1
        for key,side in (('missing','left'),('postrestore_patch','center')):
            try:archive.endpoint(key,side)
            except ValueError:pass
            else:raise AssertionError('Unsupported endpoint query accepted')
    print(json.dumps(dict(checked=True,source_only_restored=args.restore_source,full_owner_restored=False,
        actual_Rm_outer_interfaces=rec['actual_Rm_outer_interfaces'],actual_own_mixed_rows=1722,
        actual_endpoint_source_queries=queries,exact_graph_nodes=rec['exact_graph_nodes'],
        numerical_oracle_available=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()

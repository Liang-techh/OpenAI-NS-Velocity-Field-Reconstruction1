"""Inspect or restore the accepted two-switch endpoint source archive.

--restore-source hydrates the saved function graph and four endpoint views.
It does not rebuild the ancestors or restore any numerical oracle.
--rebuild-source instead reconstructs the full owner and runs its checker.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins as current


class RestoredSwitchSourceArchive:
    """Checksummed source-function handles at the four accepted endpoints."""
    def __init__(self,raw,receipt):
        self.g=g=current.source.FunctionTransportGraph()
        g.nodes=raw['exact_graph_nodes']
        g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
        assert len(g.keys)==len(g.nodes)
        assert g.nodes[g.zero.node]==dict(operation='exact_rational',numerator=0,denominator=1)
        assert g.nodes[g.one.node]==dict(operation='exact_rational',numerator=1,denominator=1)
        self.bindings=raw['switch_mixed_source_join_bindings']
        self.source_family=raw['source_family'];self.numerical_oracle_available=False
        self.acceptance_loaded=True
        self.joins={}
        for key,*_ in current.JOINS:
            self.joins[key]={}
            for side in ('left','right'):
                packet=self.bindings['joins'][key][side]
                self.joins[key][side]={group:{field:[current.outer.FunctionJet(tuple(
                    current.source.FunctionRef(g,node) for node in ids)) for ids in jets]
                    for field,jets in fields.items()} for group,fields in packet.items()}
        evidence=receipt['checks']['actual_own_mixed_joins']['actual_two_switch_own_mixed_rows']
        assert len(evidence)==492
        for row in evidence:
            for side in ('left','right'):
                jet=self.joins[row['interface']][side][row['group']][row['field']][row['physical_y_order']]
                ref=current.outer.rows(jet)[row['ordinary_Z_order']]
                assert ref.node==row[side+'_node'] and 0<=ref.node<len(g.nodes)

    def at_seam(self,interface,side):
        if interface not in self.joins or side not in ('left','right'):
            raise ValueError('An accepted switch interface and left/right endpoint are required')
        return dict(interface=interface,side=side,source_functions=self.joins[interface][side],
            numerical_oracle_available=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    parser.add_argument('--rebuild-source',action='store_true')
    args=parser.parse_args()
    if args.rebuild_source:
        import lei_ren_part1_paper_compliant_current_original_switch_mixed_source_joins_check as checker
        owner,_=current.run();checker.run(owner)
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    assert rec['all_passed'] and len(rec['checks'])==5 and all(x['passed'] for x in rec['checks'].values())
    assert all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert rec['actual_switch_interfaces']==[key for key,*_ in current.JOINS] and rec['actual_own_mixed_rows']==492
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed source',name)
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    assert len(raw['exact_graph_nodes'])==rec['exact_graph_nodes']
    queries=0
    if args.restore_source:
        archive=RestoredSwitchSourceArchive(raw,rec)
        for key,*_ in current.JOINS:
            for side in ('left','right'):
                packet=archive.at_seam(key,side)
                assert not packet['numerical_oracle_available'];queries+=1
        try:archive.at_seam('unaccepted','left')
        except ValueError:pass
        else:raise AssertionError('Unaccepted endpoint exposed')
    print(json.dumps(dict(checked=True,source_only_restored=args.restore_source,
        ancestor_owners_rebuilt=args.rebuild_source,actual_switch_interfaces=rec['actual_switch_interfaces'],
        actual_own_mixed_pairs=492,actual_endpoint_source_queries=queries,exact_graph_nodes=rec['exact_graph_nodes'],
        numerical_oracle_available=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()

"""Restore checked owners once, evaluate four cutoff cells, compose defects."""
import ast
import inspect
from pathlib import Path
import sys
import textwrap
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values as current
import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values_check as check


def build_checked(cls,cache=None):
    """Follow constructor type guards; never run an ancestor producer/checker."""
    cache={} if cache is None else cache
    if cls in cache:return cache[cls]
    module=sys.modules[cls.__module__]
    signature=inspect.signature(cls)
    required=[p for p in signature.parameters.values()
              if p.default is inspect.Parameter.empty
              and p.kind not in (p.VAR_POSITIONAL,p.VAR_KEYWORD)]
    tree=ast.parse(textwrap.dedent(inspect.getsource(cls.__init__)))
    kwargs={}
    for parameter in required:
        match=None
        for node in ast.walk(tree):
            if (isinstance(node,ast.Compare) and isinstance(node.left,ast.Call)
                and isinstance(node.left.func,ast.Name) and node.left.func.id=='type'
                and len(node.left.args)==1 and isinstance(node.left.args[0],ast.Name)
                and node.left.args[0].id==parameter.name and len(node.comparators)==1):
                target=node.comparators[0];path=[]
                while isinstance(target,ast.Attribute):path.append(target.attr);target=target.value
                if isinstance(target,ast.Name):
                    value=vars(module).get(target.id)
                    for attribute in reversed(path):value=getattr(value,attribute)
                    if inspect.isclass(value):match=value;break
        if match is None:raise ValueError('Missing checked owner type for '+cls.__name__+'.'+parameter.name)
        kwargs[parameter.name]=build_checked(match,cache)
    owner=cls(**kwargs)
    if hasattr(owner,'acceptance_loaded') and not owner.acceptance_loaded:
        raise ValueError('Checked acceptance not loaded '+cls.__name__)
    cache[cls]=owner
    return owner


def main():
    with mp.workdps(540):
        cutoff=build_checked(current.cut.CurrentOriginalN1RadialCutoffExtension)
        tails=current.external.CurrentOriginalN1ExternalPressureValues()
        owner,raw=current.run(cutoff,tails)
        receipt=check.run(owner,raw)
    with mp.workdps(15):
        public=current.CurrentOriginalN1FiveDefectValues(cutoff,tails)
        view=public.report(public.evaluate('.371'))
        assert view[current.GATE] and receipt[current.GATE] and mp.mp.dps==15
    print('ACCEPTED_N1_FIVE_DEFECT_VALUES_PUBLIC_QUERY',flush=True)


if __name__=='__main__':main()

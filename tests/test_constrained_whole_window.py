from types import SimpleNamespace
import json
from openai_ns_reconstruction import constrained_whole_window as module


def test_core_equality_run_preserves_initial_path_in_metadata(monkeypatch,tmp_path):
    def one_evaluation(fun,x0,**kwargs):
        fun(x0)
        return SimpleNamespace(message='controlled single-evaluation test')
    monkeypatch.setattr(module,'minimize',one_evaluation)
    initial='artifacts/constrained/whole_window_equalities/candidate.json'
    module.run(True,initial=initial,output=str(tmp_path),call_budget=2,max_iterations=1)
    report=json.loads((tmp_path/'training.json').read_text())
    assert report['initial_artifact']==initial
    assert report['calls']==1
    assert (tmp_path/'candidate.json').exists()

"""Mutation checks for publication links; no scientific result is changed."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('publication_layout_check', ROOT/'scripts/check_research_layout.py')
LAYOUT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LAYOUT)


def write_json(root, name, value):
    path = root/name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding='utf-8')


@pytest.fixture
def checkout(tmp_path):
    root = tmp_path/'checkout'
    for name in ('docs/CONSTRAINT_SOURCES.md', 'docs/AGENT_TASKS.md',
                 'research_baseline/_spacetime.py',
                 'artifacts/research/ST006/evidence/report.json',
                 'artifacts/research/ST006/candidate.json',
                 'experiments/root_st030/result_summary.json'):
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('{}', encoding='utf-8')
    write_json(root, 'configs/constraints.json', {'source_map': 'docs/CONSTRAINT_SOURCES.md'})
    write_json(root, 'project_status.json', {
        'active_task_file': 'docs/AGENT_TASKS.md',
        'candidate_manifest': 'artifacts/research/ST006/manifest.json',
    })
    write_json(root, 'artifacts/research/ST006/manifest.json', {
        'candidate_file': 'artifacts/research/ST006/candidate.json',
        'source_directory': 'experiments/root_st030',
        'runtime_sha256': {'research_baseline/_spacetime.py': 'placeholder'},
        'evidence_sha256': {'artifacts/research/ST006/evidence/report.json': 'placeholder'},
    })
    write_json(root, 'artifacts/research/experiment_index.json', {
        'comparison_source': '../../experiments/root_st030/result_summary.json',
        'entries': [{
            'id': 'ST006', 'manifest': 'ST006/manifest.json',
            'candidate': 'ST006/candidate.json', 'evidence': 'ST006/evidence/report.json',
        }],
    })
    return root


def test_complete_reference_map(checkout):
    result = LAYOUT.check_metadata_references(checkout)
    assert result['status'] == 'ok'
    assert result['files_checked'] == 10 and result['directories_checked'] == 1
    assert 'not scientific acceptance' in result['scope']


def test_missing_source_document(checkout):
    (checkout/'docs/CONSTRAINT_SOURCES.md').unlink()
    with pytest.raises(ValueError, match=r'constraints.source_map: missing file'):
        LAYOUT.check_metadata_references(checkout)


@pytest.mark.parametrize('value', ['', None, 'https://example.org/source', '../outside', '/tmp/source', r'C:\source'])
def test_invalid_or_external_source_map(checkout, value):
    write_json(checkout, 'configs/constraints.json', {'source_map': value})
    with pytest.raises(ValueError, match='constraints.source_map'):
        LAYOUT.check_metadata_references(checkout)


def test_symlink_escape(checkout):
    path = checkout/'docs/CONSTRAINT_SOURCES.md'
    path.unlink()
    outside = checkout.parent/'outside.md'
    outside.write_text('not in the repository', encoding='utf-8')
    path.symlink_to(outside)
    with pytest.raises(ValueError, match='path exits the repository'):
        LAYOUT.check_metadata_references(checkout)


def test_missing_indexed_evidence(checkout):
    index = json.loads((checkout/'artifacts/research/experiment_index.json').read_text())
    index['entries'][0]['evidence'] = 'ST006/evidence/not_present.json'
    write_json(checkout, 'artifacts/research/experiment_index.json', index)
    with pytest.raises(ValueError, match=r'index.evidence: missing file'):
        LAYOUT.check_metadata_references(checkout)


def test_directory_cannot_replace_source_document(checkout):
    path = checkout/'docs/CONSTRAINT_SOURCES.md'
    path.unlink()
    path.mkdir()
    with pytest.raises(ValueError, match=r'constraints.source_map: missing file'):
        LAYOUT.check_metadata_references(checkout)

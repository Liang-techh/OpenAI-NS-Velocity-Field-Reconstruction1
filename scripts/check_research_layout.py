"""Offline publication inventory and bounded JSON-reference checks.

This checks repository-local metadata paths, not Markdown anchors, external
websites, all research branches, or the mathematical correctness of a source.
"""
from __future__ import annotations

import json
from pathlib import Path, PurePosixPath
import sys
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def check_metadata_references(root: Path = ROOT) -> dict:
    """Reject missing, nonlocal, or escaping paths in publication metadata."""
    root = Path(root).resolve()
    checked: dict[str, str] = {}

    def target(base: Path, value: object, label: str, directory: bool = False) -> Path:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'{label}: expected a nonempty relative path')
        parsed = urlsplit(value)
        if (PurePosixPath(value).is_absolute() or parsed.scheme or parsed.netloc
                or parsed.query or parsed.fragment or '\\' in value):
            raise ValueError(f'{label}: expected a repository-local relative path')
        path = (base / value).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f'{label}: path exits the repository')
        exists = path.is_dir() if directory else path.is_file()
        if not exists:
            kind = 'directory' if directory else 'file'
            raise ValueError(f'{label}: missing {kind}: {path.relative_to(root)}')
        checked[path.relative_to(root).as_posix()] = 'directory' if directory else 'file'
        return path

    def document(path: Path) -> dict:
        value = json.loads(path.read_text(encoding='utf-8'))
        if not isinstance(value, dict):
            raise ValueError(f'{path.name}: expected a JSON object')
        return value

    config = document(target(root, 'configs/constraints.json', 'configuration'))
    target(root, config.get('source_map'), 'constraints.source_map')
    status = document(target(root, 'project_status.json', 'project status'))
    target(root, status.get('active_task_file'), 'project_status.active_task_file')
    manifest_path = target(root, status.get('candidate_manifest'), 'project_status.candidate_manifest')
    manifest = document(manifest_path)
    target(root, manifest.get('candidate_file'), 'manifest.candidate_file')
    target(root, manifest.get('source_directory'), 'manifest.source_directory', directory=True)
    for group in ('runtime_sha256', 'evidence_sha256'):
        entries = manifest.get(group)
        if not isinstance(entries, dict) or not entries:
            raise ValueError(f'manifest.{group}: expected a nonempty path map')
        for name in entries:
            target(root, name, f'manifest.{group}')

    index_path = target(root, 'artifacts/research/experiment_index.json', 'experiment index')
    index = document(index_path)
    target(index_path.parent, index.get('comparison_source'), 'index.comparison_source')
    entries = index.get('entries')
    if not isinstance(entries, list) or not entries:
        raise ValueError('index.entries: expected a nonempty list')
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError('index.entries: each entry must be an object')
        for key in ('manifest', 'candidate', 'evidence', 'original_comparison_evidence'):
            if key in entry:
                target(index_path.parent, entry[key], f'index.{key}')
    return {
        'status': 'ok',
        'files_checked': sum(kind == 'file' for kind in checked.values()),
        'directories_checked': sum(kind == 'directory' for kind in checked.values()),
        'scope': 'Local JSON metadata references only; not scientific acceptance',
    }


def main() -> None:
    from research_baseline import verify_integrity

    verify_integrity()
    required = [
        'README.md', 'docs/README.md', 'docs/REPOSITORY_GUIDE.md',
        'docs/RESEARCH_STATUS.md', 'docs/BRANCH_AND_PR_GUIDE.md',
        'docs/CONSTRAINT_SOURCES.md', 'docs/PUBLICATION_RECEIPT.md',
        'artifacts/research/README.md', 'artifacts/research/experiment_index.json',
        'research_baseline/__init__.py', 'experiments/root_st030/validate.py',
        'tests/test_published_baseline.py', 'tests/test_research_layout.py',
    ]
    missing = [name for name in required if not (ROOT / name).is_file()]
    if missing:
        raise ValueError('Missing publication paths: ' + ', '.join(missing))
    references = check_metadata_references()
    print(json.dumps({
        'layout': 'ok', 'required_paths': len(required),
        'metadata_references': references,
        'top_level_file_counts': {
            path.name: sum(item.is_file() for item in path.rglob('*')
                           if '.git' not in item.parts and '__pycache__' not in item.parts)
            for path in ROOT.iterdir() if path.is_dir() and not path.name.startswith('.')
        },
        'pde_validated': False,
    }, indent=2))


if __name__ == '__main__':
    main()

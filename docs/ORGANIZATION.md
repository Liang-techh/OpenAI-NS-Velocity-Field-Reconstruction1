# ST063 repository organization record

Date: 2026-09-21. User-requested maintenance, issue #1068.

## Scope

Research repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`.
Base main: `502936d119a0b582e6f271813804581abf9cc9e5`.
Base tree: `3316ce2038405449319daec99320604cfdd54ea6`.

The previous homepage emphasized only ST006 and the ST030-era release. The new homepage and current checkpoint separate latest ST063 geometry work, ST061 numerical controls, the main ST054 viewer, and the backward-compatible ST006 API. Directory responsibilities, branch routing and data availability are made explicit.

## Changed and preserved

English public landing pages, ABOUT.txt, a latest-results catalog and byte-preserved study summaries are added or updated. The existing compatibility fields in project_status.json are preserved and augmented with pointers to the new catalog.

No numerical source, candidate array, raw validation report, configuration, test, workflow or integration task is moved, deleted or modified. Legacy filenames, relative imports and open PR dependencies remain intact. No re-fit, new scientific validation, native MATLAB execution or proof is performed for this documentation task.

The source snapshots are existing Git blobs, not rewritten summaries:

| Local snapshot | Original commit / path | Blob identity |
|---|---|---|
| `research_snapshots/ST061.md` | `ad0e6dacf3851a12f4272bb4f6b282cf49506d8e`, `experiments/root_st061/README.md` | `c6e07bf8ce2711033f57d22881a27d1c9b219fb5` |
| `research_snapshots/ST063.md` | `a3d04d3467361bab7c9fc7c8c6aedf31987ee069`, `experiments/root_st063/README.md` | `f699535a5fba0d31930fbdfd4149590d6f9716c2` |

Source snapshots keep historical wording, including their then-current publication boundary. Their numeric records are not silently upgraded when this documentation PR is merged.

## Validation scope

Review the exact PR diff to confirm this scope. Compare the new paired numerical tables, raw identities and archive hashes with the pinned source records and supplied archives. Navigation checks cover the new landing pages, not every historical link across all branches. No repository-wide numerical or native-viewer test pass is inferred from documentation checks or a queued workflow.

## Description metadata

README.md contains the updated public introduction. ABOUT.txt contains the requested short English description. The separate GitHub repository About/description field needs an administration write action, which the discovered connector does not expose. No sidebar metadata success is inferred from committing ABOUT.txt.

## Continued maintenance

For a new study, first add an identity-bound catalog entry and its exact availability, then update current results and visual guidance. Do not change the compatibility loader or overwrite a frozen candidate as part of an index update. Candidate promotion, full asset publication and any future physical restructuring require their own reviewed changes.

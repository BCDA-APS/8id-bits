# Pre-migration state log — 2026-09-14

Restore points captured before migrating five capabilities from `~/ophyd` into `~/bluesky`.
Both pushed to `BCDA-APS/8id-bits`.

## Target: ~/bluesky
- Branch: `main`
- HEAD: `d5fed29` ("Merge pull request #146 ... setup-python-7")
- Working tree: **clean**
- Restore tag: **`checkpoint/bluesky-main-20260914`** -> `d5fed29` (pushed)

Restore with:
```bash
cd ~/bluesky && git checkout main && git reset --hard checkpoint/bluesky-main-20260914
```

## Source: ~/ophyd
- Branch: `dev`
- HEAD: `49d3ec0` ("Name align scans from their arguments, not from sample_info.yaml")
- `dev` was already fully pushed (local == origin/dev), so committed history needed no backup.
- Working tree: **dirty** — 10 modified tracked files, 202 deletions
  (`src/id8_common_dev/**`, `src/legacy/**`, `src/output.png`), 19 untracked paths.
- Restore branch: **`backup/ophyd-worktree-20260914`** -> `8dba7b0` (pushed)
  - Parent is `49d3ec0`; tree is an exact snapshot of the dirty working tree
    including untracked files.
  - Built with a temporary `GIT_INDEX_FILE` so `~/ophyd` HEAD, index and working
    tree were never modified (this checkout is live on pearl over NFS).

Restore/inspect with:
```bash
git show backup/ophyd-worktree-20260914 --stat
git checkout backup/ophyd-worktree-20260914 -- <path>   # recover a single file
```

## Capabilities being migrated (ophyd -> bluesky)
1. Modifying both `sample_info` and `measurement_info` in CSV files
2. Writing metadata after Ctrl+C out of an acquisition
3. Eiger External Series and Rigaku Fast Transfer Mode
4. Start a session without disturbing a running acquisition
   (no Eiger warm-up / throwaway frame)
5. The trio detector acquisition plan

## Host constraints
- `kouga` (164.54.116.55) — GitHub reachable, **no beamline PVs**
- `pearl` (10.*) — beamline PVs, **no internet**; runs Bluesky/Ophyd
- `~/bluesky` and `~/ophyd` on kouga are shared to pearl over NFS: edits here are
  live there immediately. Do not switch branches while a session is running.

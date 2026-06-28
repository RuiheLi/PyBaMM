# -*- coding: utf-8 -*-
"""
Backup P3_R14 and reorganize code by manuscript figure number.

Usage (from P3_R14):
  python reorganize_by_figure.py
  python reorganize_by_figure.py --dry-run
"""
from __future__ import annotations

import argparse
import shutil
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
BACKUP = PARENT / f"P3_R14_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

# (relative source from ROOT, figure id, optional note)
# figure id: main "03".."10", SI "S01".."S16", shared "shared", archive "archive"
MOVES: list[tuple[str, str, str]] = [
    # --- Main text ---
    ("1EC1DMC_Val_GITT_2C.ipynb", "03", "Main Fig. 3 — 2C GITT validation"),
    ("1EC1DMC_Val_GITT_2C_Xitilde_func.ipynb", "03", "Main Fig. 3 — GITT 2C (Xi_tilde variant)"),
    ("Reproduce_sol_seg.ipynb", "04", "Main Fig. 4 — solvent segregation / Hittorf"),
    ("Reproduce_sol_seg_Xitilde_func.ipynb", "04", "Main Fig. 4 — sol. seg. (Xi_tilde)"),
    ("Reproduce_sol_seg_Xitilde_func.executed.ipynb", "04", "Main Fig. 4 — executed copy"),
    ("run_fig5_45C_parity.py", "05", "Main Fig. 5 — 4.5C flux/source (ts_dis parity)"),
    ("1EC1DMC_Rate_Performance.ipynb", "05", "Main Fig. 5–9 hub notebook (see README)"),
    ("Debug_Xitilde.ipynb", "06", "Main Fig. 6 — internal gradients / electrolyte properties"),
    ("1EC1DMC_Rate_Performance_AnSol.ipynb", "07", "Main Fig. 7 — C-rate profiles (from saved sol)"),
    ("run_fig8_discharge_voltage.py", "08", "Main Fig. 8 — discharge voltage curves"),
    ("Rate_Performance_Get_Sol.ipynb", "09", "Main Fig. 9 — run rate sweep / save solutions"),
    ("1EC1DMC_OneCycAge_1C.ipynb", "10", "Main Fig. 10 — ageing @ 1C"),
    ("1EC1DMC_OneCycAge_2C.ipynb", "10", "Main Fig. 10 — ageing @ 2C charge"),
    ("1EC1DMC_OneCycAge_3C_sol_diff.ipynb", "10", "Main Fig. 10 — ageing, sol-diff SEI"),
    ("Reload_Age.ipynb", "10", "Main Fig. 10 — reload ageing results"),
    # --- SI ---
    ("1EC1DMC_Val_GITT_1C.ipynb", "S01", "SI Fig. S1 — 1C GITT (double-solvent)"),
    ("Para/Electrolyte_parameters_Paper.ipynb", "S10", "SI Fig. S10 — electrolyte κ, D, t+"),
    ("Para/Case_2_Para.ipynb", "S12", "SI Fig. S12/S13/S14/S15/S16 — parameter & sensitivity"),
    ("Para/Compare_3Cases.ipynb", "S09", "SI Fig. S9 — Nyman κ,D vs model"),
    ("Debug_Dec_e.ipynb", "S11", "SI Fig. S11 — EC diffusivity sensitivity"),
    # --- Shared tools ---
    ("debug_long_sims.py", "shared", "CLI: reproduce_seg + onecyc (Fig. 4–5)"),
    ("run_all_notebooks.py", "shared", "Batch notebook runner"),
    ("Long_int.py", "shared", "Long integration helper"),
    ("smoke_test_p3_simulation.py", "shared", "Smoke test"),
    ("smoke_short_one_cycle_age.py", "shared", "Short ageing smoke"),
    ("parity_source.py", "shared", "Parity reference runs"),
    ("parity_source_m6.py", "shared", "Parity m6"),
    ("parity_source_m7.py", "shared", "Parity m7"),
    ("parity_source_m8.py", "shared", "Parity m8"),
    ("parity_source_m9.py", "shared", "Parity m9"),
    ("parity_source_m9b.py", "shared", "Parity m9b"),
    ("parity_compare.py", "shared", "Parity compare"),
    ("parity_compare_m6.py", "shared", "Parity compare m6"),
    ("parity_compare_m7.py", "shared", "Parity compare m7"),
    ("parity_compare_m8.py", "shared", "Parity compare m8"),
    ("parity_compare_m9.py", "shared", "Parity compare m9"),
    ("parity_compare_m9b.py", "shared", "Parity compare m9b"),
    ("parity_compare_m10.py", "shared", "Parity compare m10"),
    ("extract_docx_equations.py", "shared", "Docx equation extract"),
    ("reorganize_by_figure.py", "shared", "This reorganisation script"),
    # --- Archive (logs, patches, docs, data) ---
    ("_apply_notebook_fixes.py", "archive", "One-off patch script"),
    ("_inject_nb_progress.py", "archive", "One-off patch script"),
    ("_patch_debug_xitilde_findclose.py", "archive", "One-off patch script"),
    ("_patch_long_notebooks.py", "archive", "One-off patch script"),
    ("hourly_progress.log", "archive", "Old log"),
    ("hourly_progress_tracker.ps1", "archive", "Old tracker"),
    ("run_all_notebooks.log", "archive", "Old log"),
    ("run_all_notebooks_stdout.txt", "archive", "Old log"),
    ("run_all_notebooks_summary.txt", "archive", "Old log"),
    ("run_fixed_three_stdout.txt", "archive", "Old log"),
    ("run_seg_19_resume.log", "archive", "Old log"),
    ("LONG_RUN_RESUME.md", "archive", "Notes"),
    ("ECdrag2_dimensionalized_equations.md", "archive", "Migration doc"),
    ("ECdrag2_migration_status.md", "archive", "Migration doc"),
    ("ECdrag2_symbol_mapping_for_port.md", "archive", "Migration doc"),
    ("DoubleTransport_Dimensionalization.md", "archive", "Migration doc"),
    ("parity_source.npz", "archive", "Parity data"),
    ("parity_source_m6.npz", "archive", "Parity data"),
    ("parity_source_m7.npz", "archive", "Parity data"),
    ("parity_source_m8.npz", "archive", "Parity data"),
    ("parity_source_m9.npz", "archive", "Parity data"),
    ("parity_source_m9b.npz", "archive", "Parity data"),
]

# Double_SimSave subfolder -> figure output folder
SIMSAVE_MOVES: list[tuple[str, str]] = [
    ("Fig5_parity", "05"),
    ("Fig8_discharge_voltage", "08"),
    ("Double_TransRate_Performance", "05"),  # primary; also 6/7/9 outputs
    ("Double_Trans", "04"),
    ("Double_TransOneCycAge_240325", "10"),
    ("Double_TransOneCycAge_240319_SolDiff", "10"),
    ("Double_TransDebug_Xi_tilde", "06"),
    ("Double_TransRate_Performance_Debug_EC", "S11"),
]

README = """# P3_R14 — code organised by manuscript figure

Backup of the pre-reorganisation tree: `{backup}`

## Main text (current numbering after revision)

| Fig. | Folder / prefix | Content |
|------|-------------------|---------|
| 1 | — | Schematic (outside P3_R14; see 审稿意见 `电池模型示意图.pptx`) |
| 2 | — | Conceptual summary (no dedicated simulation notebook) |
| 3 | `fig_03_*` | 2C GITT model validation |
| 4 | `fig_04_*` | Solvent segregation / Hittorf-style discharge |
| 5 | `fig_05_*` | 4.5C EC/Li⁺ flux & source terms |
| 6 | `fig_06_*` | End-of-discharge internal gradients (κ, η, c profiles) |
| 7 | `fig_07_*` | High-Dx case vs C-rate (concentration / overpotential) |
| 8 | `fig_08_*` | Discharge voltage profiles (0.1–2C; 4.5C optional) |
| 9 | `fig_09_*` | Rate performance & end-of-discharge overpotentials |
| 10 | `fig_10_*` | Ageing: c(EC), j_SEI, capacity loss |

### Multi-figure notebook

`figures/main/fig_05_1EC1DMC_Rate_Performance.ipynb` also contains cells for **Fig. 6, 7, 9** (and SI flux panels). See cell comments `# Fig. 4` / `# Fig. 5` / `Fig_S2` in the notebook.

## Supplementary Information

| SI Fig. | Prefix | Primary file |
|---------|--------|--------------|
| S1 | `fig_S01_*` | 1C GITT (double-solvent) |
| S2 | — | Same GITT notebooks, single-solvent branch |
| S3–S7 | `fig_05_*` / `fig_07_*` | Low/single cases in Rate_Performance notebooks |
| S8 | `fig_10_*` | Ageing without solvent consumption (variant in OneCycAge) |
| S9 | `fig_S09_*` | Compare_3Cases / Nyman electrolyte data |
| S10 | `fig_S10_*` | Electrolyte_parameters_Paper |
| S11 | `fig_S11_*` | Debug_Dec_e |
| S12–S16 | `fig_S12_*` | Case_2_Para (conductivity, OCV, solid D, Ve fit) |

## Layout

```
P3_R14/
  figures/main/     fig_XX_<original_name>
  figures/SI/       fig_SXX_<original_name>
  figures/shared/   tools used across figures
  figures/archive/  logs, parity snapshots, migration notes
  outputs/main/     simulation figures & pickles by fig number
  outputs/SI/
```

## Notes

- Run scripts from **P3_R14** root so `Fun_P3` / PyBaMM paths still resolve.
- After moving, open notebooks once and check `BasicPath` / `Save_Fig` paths.
- Regenerate Fig. 8 PNG: `python figures/main/fig_08_run_fig8_discharge_voltage.py --plot-only`
- Regenerate Fig. 5 profile: `python figures/main/fig_05_run_fig5_45C_parity.py --plot-only`
"""


def _dest_path(fig_id: str, src_name: str) -> Path:
    name = Path(src_name).name
    if fig_id == "shared":
        return ROOT / "figures" / "shared" / name
    if fig_id == "archive":
        return ROOT / "figures" / "archive" / name
    if fig_id.startswith("S"):
        return ROOT / "figures" / "SI" / f"fig_{fig_id}_{name}"
    return ROOT / "figures" / "main" / f"fig_{fig_id}_{name}"


def _sim_dest(fig_id: str, folder_name: str) -> Path:
    if fig_id.startswith("S"):
        return ROOT / "outputs" / "SI" / f"fig_{fig_id}_{folder_name}"
    return ROOT / "outputs" / "main" / f"fig_{fig_id}_{folder_name}"


def backup_tree(dry_run: bool) -> Path:
    print(f"Backup -> {BACKUP}")
    if dry_run:
        return BACKUP
    if BACKUP.exists():
        raise SystemExit(f"Backup path already exists: {BACKUP}")
    shutil.copytree(ROOT, BACKUP, ignore=shutil.ignore_patterns("P3_R14_backup_*"))
    print("Backup complete.")
    return BACKUP


def reorganize(dry_run: bool) -> None:
    manifest_lines = ["# Move manifest", ""]
    for src_rel, fig_id, note in MOVES:
        src = ROOT / src_rel
        if not src.exists():
            manifest_lines.append(f"SKIP (missing): {src_rel}  # {note}")
            continue
        dst = _dest_path(fig_id, src_rel)
        manifest_lines.append(f"{src_rel} -> {dst.relative_to(ROOT)}  # {note}")
        if dry_run:
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            raise SystemExit(f"Destination exists: {dst}")
        shutil.move(str(src), str(dst))

    sim_root = ROOT / "Double_SimSave"
    if sim_root.is_dir():
        for folder, fig_id in SIMSAVE_MOVES:
            src = sim_root / folder
            if not src.exists():
                manifest_lines.append(f"SKIP sim: Double_SimSave/{folder}")
                continue
            dst = _sim_dest(fig_id, folder)
            manifest_lines.append(f"Double_SimSave/{folder} -> {dst.relative_to(ROOT)}")
            if dry_run:
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists():
                raise SystemExit(f"Sim dest exists: {dst}")
            shutil.move(str(src), str(dst))
        # move any remaining Double_SimSave items
        remaining = list(sim_root.iterdir()) if sim_root.exists() else []
        if remaining and not dry_run:
            misc = ROOT / "outputs" / "misc"
            misc.mkdir(parents=True, exist_ok=True)
            for item in remaining:
                shutil.move(str(item), str(misc / item.name))
            try:
                sim_root.rmdir()
            except OSError:
                pass

    # Para/ leftovers
    para = ROOT / "Para"
    if para.is_dir() and any(para.iterdir()):
        misc_para = ROOT / "figures" / "archive" / "Para_leftover"
        if not dry_run:
            misc_para.mkdir(parents=True, exist_ok=True)
            for item in para.iterdir():
                shutil.move(str(item), str(misc_para / item.name))
            try:
                para.rmdir()
            except OSError:
                pass
        manifest_lines.append("Para/ -> figures/archive/Para_leftover/")

    readme = README.format(backup=BACKUP)
    manifest_lines.extend(["", readme])
    out_readme = ROOT / "README_FIGURE_INDEX.md"
    out_manifest = ROOT / "figures" / "MOVE_MANIFEST.md"
    if not dry_run:
        out_readme.write_text(readme, encoding="utf-8")
        out_manifest.parent.mkdir(parents=True, exist_ok=True)
        out_manifest.write_text("\n".join(manifest_lines), encoding="utf-8")
        print(f"Wrote {out_readme.name} and {out_manifest}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-backup", action="store_true", help="Only reorganise (not recommended)")
    args = ap.parse_args()
    if args.dry_run:
        print("=== DRY RUN ===")
        reorganize(dry_run=True)
        return
    if not args.skip_backup:
        backup_tree(dry_run=False)
    reorganize(dry_run=False)
    print("Reorganisation complete.")


if __name__ == "__main__":
    main()

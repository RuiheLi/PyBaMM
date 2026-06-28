# ECdrag2 — P3_R14 manuscript workflow

**Primary path for reproducing the paper:** [`wip/Rio_Code/P3_R14/figures/`](figures/)

Notebooks and helper scripts there are organised by main-text and supplementary figure number. Generated plots are written under [`outputs/`](outputs/) (same figure numbering).

## Quick start

1. Install this fork’s PyBaMM (repo root) and dependencies (`requirements.txt`).
2. `cd wip/Rio_Code/P3_R14`
3. Open the notebook or script listed for the figure you need — see [`README_FIGURE_INDEX.md`](README_FIGURE_INDEX.md).
4. Run from **P3_R14 root** so `Fun_P3` and PyBaMM import paths resolve.

Example (Fig. 8 discharge voltage, plot-only):

```bash
python figures/main/fig_08_run_fig8_discharge_voltage.py --plot-only
```

Simulation checkpoints (`*.pkl`) under `outputs/` are not tracked in git (too large); regenerate them by running the corresponding notebook or script.

## Layout

| Path | Role |
|------|------|
| `figures/main/` | Main-text figures (Fig. 3–10) |
| `figures/SI/` | Supplementary figures |
| `figures/shared/` | Cross-figure tools (smoke tests, parity, long-run debug) |
| `figures/archive/` | Migration notes, logs, one-off patches |
| `outputs/main/`, `outputs/SI/` | Simulation outputs aligned with figure numbers |

Full figure-to-file mapping: [`README_FIGURE_INDEX.md`](README_FIGURE_INDEX.md).

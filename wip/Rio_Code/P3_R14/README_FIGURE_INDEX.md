# P3_R14 — code organised by manuscript figure

> **Paper reproduction (ECdrag2):** start in [`figures/`](figures/) — notebooks and scripts are grouped by figure number; see [`README.md`](README.md) and the table below.

Backup of the pre-reorganisation tree: `D:\LRHWork\Model_RH\PyBaMM-ECdrag2\wip\Rio_Code\P3_R14_backup_20260614_180652`

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

"""
Canonical path bootstrap for P3_R14 notebooks and scripts after figure reorganisation.

Usage in a notebook (first path cell):

    from p3_bootstrap import bootstrap, outputs_dir
    paths = bootstrap()          # inserts ECdrag2 pybamm + Fun_P3 on sys.path
    BasicPath = outputs_dir("03", "Validation_240318_Fun_Xi_tidle")
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

# figures/shared/p3_bootstrap.py -> P3_R14
P3_ROOT = Path(__file__).resolve().parents[2]
# -> PyBaMM-ECdrag2 repo root
REPO_ROOT = P3_ROOT.parents[2]
FUN_P3 = REPO_ROOT / "wip" / "Rio_Code" / "Fun_P3"


def bootstrap(*, force_ecdrag_pybamm: bool = True) -> dict[str, Path]:
    """Prepend ECdrag2 fork and Fun_P3; return path dict."""
    repo = str(REPO_ROOT)
    fun = str(FUN_P3)
    if repo not in sys.path:
        sys.path.insert(0, repo)
    if fun not in sys.path:
        sys.path.insert(0, fun)

    if force_ecdrag_pybamm:
        for mod in [m for m in list(sys.modules) if m == "pybamm" or m.startswith("pybamm.")]:
            sys.modules.pop(mod, None)

    return {
        "P3_ROOT": P3_ROOT,
        "REPO_ROOT": REPO_ROOT,
        "FUN_P3": FUN_P3,
    }


def outputs_dir(fig_id: str, subfolder: str = "") -> str:
    """
    fig_id: main figure '03'..'10' or SI 'S01' etc.
    Returns absolute path with trailing subfolder, creates directory.
    """
    if fig_id.upper().startswith("S"):
        base = P3_ROOT / "outputs" / "SI" / f"fig_{fig_id.upper()}"
    else:
        base = P3_ROOT / "outputs" / "main" / f"fig_{fig_id.zfill(2)}"
    path = base / subfolder if subfolder else base
    os.makedirs(path, exist_ok=True)
    return str(path) + ("/" if subfolder and not subfolder.endswith("/") else "")


# Explicit re-exports for static analysis (Pylance / Pyright)
def import_fun_p3():
    bootstrap(force_ecdrag_pybamm=False)
    from Fun_P3 import (  # noqa: WPS433
        Add_var,
        Para_init_Dict,
        Plot_Comp_GITT_Overall,
        recursive_scan,
    )
    return Add_var, Para_init_Dict, Plot_Comp_GITT_Overall, recursive_scan

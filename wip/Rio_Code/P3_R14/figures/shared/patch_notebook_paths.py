# -*- coding: utf-8 -*-
"""Fix REPO_ROOT in shared scripts and patch notebook path cells after reorganisation."""
from __future__ import annotations

import json
import re
from pathlib import Path

P3 = Path(__file__).resolve().parents[2]
SHARED = P3 / "figures" / "shared"

OLD_PATH_CELL = """import sys  
str_path_0 = os.path.abspath(os.path.join(pybamm.__path__[0],'..'))
str_path_1 = os.path.abspath(
   os.path.join(str_path_0,"wip/Rio_Code/Fun_P3"))
sys.path.append(str_path_1) 
from Fun_P3 import *"""

NEW_PATH_CELL = """# P3_R14 path bootstrap (post figure reorganisation)
import sys
sys.path.insert(0, r"{shared}")
from p3_bootstrap import bootstrap, outputs_dir
bootstrap()
from Fun_P3 import *  # noqa: F403
"""

# fig number hints for BasicPath replacement
FIG_OUTPUT_SUB = {
    "fig_03": ("03", "Validation_240318_Fun_Xi_tidle/"),
    "fig_04": ("04", "Reproduce_Sol_Seg/"),
    "fig_05": ("05", "Double_TransRate_Performance/"),
    "fig_06": ("06", "Double_TransDebug_Xi_tilde/"),
    "fig_07": ("07", "Double_TransRate_Performance/"),
    "fig_09": ("09", "Double_TransRate_Performance/"),
    "fig_10": ("10", "Double_TransOneCycAge/"),
    "fig_S01": ("S01", "Validation_240318/"),
}


def _repair_mangled_source(src: str) -> tuple[str, bool]:
    """Split accidental 'F403BasicPath' glue from a bad patch merge."""
    if "F403BasicPath" not in src:
        return src, False
    fixed = src.replace(
        "from Fun_P3 import *  # noqa: F403BasicPath",
        "from Fun_P3 import *  # noqa: F403\nBasicPath",
    )
    return fixed, fixed != src


def _patch_py_repo_depth(path: Path) -> bool:
    """No-op placeholder kept for scripts that now use p3_bootstrap directly."""
    return False


def _patch_notebook(nb_path: Path) -> list[str]:
    changes: list[str] = []
    data = json.loads(nb_path.read_text(encoding="utf-8"))
    name = nb_path.name
    fig_key = name.split("_", 2)[0] + "_" + name.split("_", 2)[1] if name.startswith("fig_") else ""

    for cell in data.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        src = "".join(cell.get("source", []))
        repaired, did_repair = _repair_mangled_source(src)
        if did_repair:
            cell["source"] = [line + "\n" for line in repaired.splitlines()]
            if cell["source"] and not cell["source"][-1].endswith("\n"):
                cell["source"][-1] += "\n"
            changes.append("repair_f403_basicpath")
            src = repaired
        if "from Fun_P3 import" in src and "p3_bootstrap" not in src:
            if OLD_PATH_CELL.replace("\n", "\n") in src or "wip/Rio_Code/Fun_P3" in src:
                new_src = NEW_PATH_CELL.format(shared=str(SHARED).replace("\\", "/"))
                # preserve BasicPath / Target lines below import block
                rest = re.sub(
                    r"import sys.*?from Fun_P3 import \*\n",
                    "",
                    src,
                    count=1,
                    flags=re.DOTALL,
                )
                if "BasicPath" in rest and "Double_SimSave" in rest:
                    hint = FIG_OUTPUT_SUB.get(fig_key[:6] if fig_key.startswith("fig_S") else fig_key[:6])
                    if hint:
                        fig_id, sub = hint
                        rest = re.sub(
                            r'BasicPath = os\.path\.join\([^)]+\)\n',
                            f'BasicPath = outputs_dir("{fig_id}", "{sub}")\n',
                            rest,
                            count=1,
                        )
                merged = new_src.rstrip("\n") + "\n" + rest.lstrip("\n")
                cell["source"] = [line + "\n" for line in merged.splitlines()]
                if cell["source"] and not cell["source"][-1].endswith("\n"):
                    cell["source"][-1] += "\n"
                changes.append("path_cell")
        elif "Double_SimSave" in src:
            new_src = src.replace("Double_SimSave", "outputs/main")
            if new_src != src:
                cell["source"] = [line + "\n" for line in new_src.splitlines()]
                changes.append("double_simsave")

    if changes:
        nb_path.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
    return changes


def main() -> None:
    patched_py = []
    for py in (P3 / "figures").rglob("*.py"):
        if _patch_py_repo_depth(py):
            patched_py.append(str(py.relative_to(P3)))

    nb_changes = {}
    for nb in (P3 / "figures").rglob("*.ipynb"):
        ch = _patch_notebook(nb)
        if ch:
            nb_changes[str(nb.relative_to(P3))] = ch

    print("Patched .py REPO depth:", len(patched_py))
    for p in patched_py:
        print(" ", p)
    print("Patched notebooks:", len(nb_changes))
    for k, v in sorted(nb_changes.items()):
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()

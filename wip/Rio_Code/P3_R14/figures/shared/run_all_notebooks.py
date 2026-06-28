"""Execute P3_R14 figure notebooks in dependency order (env_ecdrag / ECdrag2 pybamm).

Logs to ``figures/shared/run_all_notebooks.log``.

Examples:
  python run_all_notebooks.py --list
  python run_all_notebooks.py --only fig_03_1EC1DMC_Val_GITT_2C.ipynb
  python run_all_notebooks.py --smoke   # import checks only
"""
from __future__ import annotations

import argparse
import copy
import datetime as _dt
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
P3_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
REPO_ROOT = os.path.abspath(os.path.join(P3_ROOT, "..", "..", ".."))
FUN_P3 = os.path.join(REPO_ROOT, "wip", "Rio_Code", "Fun_P3")

# Relative to P3_ROOT — dependency-friendly order
ORDER = [
    "figures/SI/fig_S10_Electrolyte_parameters_Paper.ipynb",
    "figures/SI/fig_S09_Compare_3Cases.ipynb",
    "figures/main/fig_04_Reproduce_sol_seg.ipynb",
    "figures/main/fig_09_Rate_Performance_Get_Sol.ipynb",
    "figures/main/fig_03_1EC1DMC_Val_GITT_2C.ipynb",
    "figures/SI/fig_S01_1EC1DMC_Val_GITT_1C.ipynb",
    "figures/main/fig_03_1EC1DMC_Val_GITT_2C_Xitilde_func.ipynb",
    "figures/main/fig_06_Debug_Xitilde.ipynb",
    "figures/main/fig_07_1EC1DMC_Rate_Performance_AnSol.ipynb",
    "figures/main/fig_05_1EC1DMC_Rate_Performance.ipynb",
    "figures/main/fig_10_1EC1DMC_OneCycAge_1C.ipynb",
    "figures/main/fig_10_1EC1DMC_OneCycAge_2C.ipynb",
    "figures/main/fig_10_1EC1DMC_OneCycAge_3C_sol_diff.ipynb",
    "figures/main/fig_10_Reload_Age.ipynb",
    "figures/SI/fig_S11_Debug_Dec_e.ipynb",
    "figures/SI/fig_S12_Case_2_Para.ipynb",
    "figures/main/fig_04_Reproduce_sol_seg_Xitilde_func.ipynb",
]

SETUP_SRC = (
    "# auto-injected by run_all_notebooks.py\n"
    "import sys, os\n"
    f"_repo = r{REPO_ROOT!r}\n"
    f"_fun = r{FUN_P3!r}\n"
    f"_shared = r{HERE!r}\n"
    f"_p3 = r{P3_ROOT!r}\n"
    "for _p in (_repo, _fun, _shared):\n"
    "    if _p not in sys.path:\n"
    "        sys.path.insert(0, _p)\n"
    "for _mod in [m for m in list(sys.modules) if m == 'pybamm' or m.startswith('pybamm.')]:\n"
    "    sys.modules.pop(_mod, None)\n"
    "import pybamm\n"
    "from p3_bootstrap import bootstrap\n"
    "bootstrap(force_ecdrag_pybamm=False)\n"
    "print('pybamm at', pybamm.__path__[0])\n"
)


def _stamp():
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _resolve_nb(name: str) -> str:
    if os.path.isabs(name):
        return name
    if os.path.sep in name or "/" in name:
        return os.path.join(P3_ROOT, name.replace("/", os.path.sep))
    # bare filename: search figures/
    for root, _, files in os.walk(os.path.join(P3_ROOT, "figures")):
        if name in files:
            return os.path.join(root, name)
    return os.path.join(P3_ROOT, name)


def smoke_test() -> int:
    print("P3_ROOT  ", P3_ROOT)
    print("REPO_ROOT", REPO_ROOT)
    print("FUN_P3   ", FUN_P3)
    sys.path.insert(0, REPO_ROOT)
    sys.path.insert(0, FUN_P3)
    sys.path.insert(0, HERE)
    from p3_bootstrap import bootstrap  # noqa: WPS433

    bootstrap()
    import pybamm  # noqa: WPS433
    from Fun_P3 import recursive_scan, Para_init_Dict  # noqa: WPS433

    assert "PyBaMM-ECdrag2" in pybamm.__path__[0] or "Model_RH" in pybamm.__path__[0]
    print("pybamm   ", pybamm.__path__[0])
    print("recursive_scan OK")
    # quick script smoke
    for script in ("smoke_test_p3_simulation.py",):
        path = os.path.join(HERE, script)
        if os.path.isfile(path):
            print(f"running {script} ...")
            import subprocess

            r = subprocess.run([sys.executable, path], cwd=P3_ROOT, check=False)
            if r.returncode != 0:
                return r.returncode
    print("SMOKE OK")
    return 0


def _execute(nb_path, log, cell_timeout):
    import nbformat
    from nbclient.exceptions import CellExecutionError
    from nbconvert.preprocessors import ExecutePreprocessor

    name = os.path.relpath(nb_path, P3_ROOT)
    log.write(f"[{_stamp()}] START {name}\n")
    log.flush()
    t0 = time.time()

    nb = nbformat.read(nb_path, as_version=4)
    setup_cell = nbformat.v4.new_code_cell(source=SETUP_SRC)
    setup_cell["metadata"]["_run_all_notebooks_injected"] = True
    nb_to_run = copy.deepcopy(nb)
    nb_to_run.cells.insert(0, setup_cell)

    ep = ExecutePreprocessor(
        timeout=cell_timeout,
        kernel_name="python3",
        allow_errors=False,
    )
    status = "ok"
    err_tail = ""
    try:
        ep.preprocess(nb_to_run, {"metadata": {"path": os.path.dirname(nb_path)}})
    except CellExecutionError as exc:
        status = "fail"
        err_tail = "\n".join(str(exc).strip().splitlines()[-25:])
    except Exception as exc:  # noqa: BLE001
        status = "fail"
        err_tail = f"{type(exc).__name__}: {exc}"
    dt = time.time() - t0

    cleaned = copy.deepcopy(nb_to_run)
    cleaned.cells = [
        c
        for c in cleaned.cells
        if not c.get("metadata", {}).get("_run_all_notebooks_injected")
    ]
    nbformat.write(cleaned, nb_path)

    if status == "ok":
        log.write(f"[{_stamp()}] OK in {dt:.1f}s\n\n")
    else:
        log.write(
            f"[{_stamp()}] {status.upper()} in {dt:.1f}s\n  err tail:\n{err_tail}\n\n"
        )
    log.flush()
    return status, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timeout", type=int, default=7200)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.list:
        for rel in ORDER:
            full = _resolve_nb(rel)
            mark = "OK" if os.path.isfile(full) else "MISSING"
            print(f"  [{mark}] {rel}")
        return 0

    if args.smoke:
        return smoke_test()

    log_path = os.path.join(HERE, "run_all_notebooks.log")
    summary_path = os.path.join(HERE, "run_all_notebooks_summary.txt")

    todo = args.only if args.only else ORDER
    results = []
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(f"[{_stamp()}] BATCH START, {len(todo)} notebook(s)\n")
        log.write(f"  P3_ROOT   = {P3_ROOT}\n")
        log.write(f"  REPO_ROOT = {REPO_ROOT}\n\n")
        log.flush()
        for name in todo:
            nb = _resolve_nb(name)
            if not os.path.exists(nb):
                log.write(f"[{_stamp()}] SKIP {name} (missing @ {nb})\n\n")
                results.append((name, "missing", 0.0))
                continue
            status, dt = _execute(nb, log, args.timeout)
            results.append((name, status, dt))
        log.write(f"[{_stamp()}] BATCH END\n")

    with open(summary_path, "w", encoding="utf-8") as summary:
        summary.write(f"# Notebook batch summary - {_stamp()}\n\n")
        summary.write(f"{'status':10s}  {'time[s]':>10s}  notebook\n")
        for name, status, dt in results:
            summary.write(f"{status:10s}  {dt:10.1f}  {name}\n")
        ok = sum(1 for _, s, _ in results if s == "ok")
        summary.write(f"\n{ok}/{len(results)} succeeded\n")

    print(f"log     -> {log_path}")
    print(f"summary -> {summary_path}")
    for name, status, dt in results:
        print(f"  {status:8s}  {dt:8.1f}s  {name}")
    return 0 if all(s in ("ok", "missing") for _, s, _ in results) else 1


if __name__ == "__main__":
    sys.exit(main() or 0)

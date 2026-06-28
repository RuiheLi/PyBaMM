"""
Rerun 4.5C high-Dx discharge (Fig. 5 case) with configurable experiment period.

Default: ts_dis=0.1 s only. Saves .pkl after each solve.

  python run_fig5_45C_parity.py --ts-dis 0.1 --mesh-r-n 80
  python run_fig5_45C_parity.py --plot-only --ts-dis 0.5 0.1
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_SHARED = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "shared"))
if _SHARED not in sys.path:
    sys.path.insert(0, _SHARED)
from p3_bootstrap import bootstrap  # noqa: E402

bootstrap()

import pybamm  # noqa: E402

from Fun_P3 import Add_var, Para_init_Dict, recursive_scan  # noqa: E402

_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = _ROOT / "outputs" / "main" / "fig_05_Fig5_parity"
STATUS_LOG = OUT_DIR / "run_status.log"
RATE = 4.5
CROSS_D = 1e-11


def _log(msg: str) -> None:
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(STATUS_LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def _build_para(mesh_r_n: int) -> dict:
    para_dd = {
        "Mesh list": [[20, 10, 20, mesh_r_n, 20]],
        "Para_Set": ["Li2023_ECdrag"],
        "Contact resistance [Ohm]": [6e-3],
        "Initial Neg SOC": [0.8841301667966484],
        "Initial Pos SOC": [0.23552755074598045],
        "Model option": [
            {
                "SEI": "constant",
                "SEI film resistance": "distributed",
                "SEI porosity change": "true",
                "solvent diffusion": "double spatial consume w refill",
                "electrolyte conductivity": "sol full",
                "contact resistance": "true",
            }
        ],
        "Lithium ion EC cross diffusivity [m2.s-1]": [CROSS_D],
        "EC transference number": ["EC_transference_number_3"],
    }
    stack: list = []
    recursive_scan(stack, para_dd, list(para_dd.keys()), {})
    return stack[0]


def _solve(para_i: dict, ts_dis: float | None, mesh_r_n: int):
    cycle_pack, para_used = Para_init_Dict(para_i)
    mesh_list, model_options = cycle_pack
    model = pybamm.lithium_ion.DFN(options=model_options)
    model = Add_var(para_used, model)
    v_max, v_min = 4.2, 2.5
    if ts_dis is None:
        step = f"Discharge at {RATE} C until {v_min} V"
        label = "default"
    else:
        step = f"Discharge at {RATE} C until {v_min} V ({ts_dis} second period)"
        label = f"period={ts_dis}s"
    exp = pybamm.Experiment([(f"Hold at {v_max} V until C/100", step)])
    sim = pybamm.Simulation(
        model,
        experiment=exp,
        parameter_values=para_used,
        solver=pybamm.CasadiSolver(
            return_solution_if_failed_early=True,
            extra_options_setup={"max_num_steps": 1_000_000},
        ),
        var_pts={
            "x_n": mesh_list[0],
            "x_s": mesh_list[1],
            "x_p": mesh_list[2],
            "r_n": mesh_list[3],
            "r_p": mesh_list[4],
        },
    )
    t0 = time.perf_counter()
    _log(f"fig5 solve {label} BEGIN")
    sol = sim.solve()
    wall = time.perf_counter() - t0
    _log(f"fig5 solve {label} OK wall={wall:.1f}s")
    return sol, label, wall


def _profile_at_end(step, key: str) -> tuple[np.ndarray, np.ndarray]:
    x = step["x [m]"].entries[:, 0] * 1e3
    y = step[key].entries[:, -1]
    return x, y


def _plot_profiles(solutions: list, out_path: str) -> None:
    var_keys = [
        ("Li+ source term [mol.m-3.s-1]", "Li+ source"),
        ("Minus div Li+ flux by migration [mol.m-3.s-1]", "−∇·(Li+ migration)"),
        ("Minus div Li+ flux by solvent [mol.m-3.s-1]", "−∇·(Li+ cross-diffusion)"),
    ]
    fig, axs = plt.subplots(len(var_keys), 1, figsize=(7, 2.2 * len(var_keys)), sharex=True)
    t_end = 0.0
    for sol, label, _ in solutions:
        step = sol.cycles[0].steps[1]
        t_end = float(step["Time [s]"].entries[-1])
        for ax, (key, title) in zip(axs, var_keys):
            x_mm, y = _profile_at_end(step, key)
            ax.plot(x_mm, y, label=label, lw=1.2)
            ax.set_ylabel(title, fontsize=9)
            ax.grid(True, alpha=0.3)
    axs[-1].set_xlabel("Position along cell [mm]")
    axs[0].legend(fontsize=8)
    fig.suptitle(
        f"Fig. 5 @ {RATE}C, high D_e,EC — end of discharge (t={t_end:.0f}s)",
        fontsize=10,
    )
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    _log(f"fig5 plot saved {out_path}")


def _pkl_path(ts: float | None, mesh_r_n: int) -> str:
    tag = "default" if ts is None else f"period{ts}s"
    return os.path.join(OUT_DIR, f"fig5_45C_mesh{mesh_r_n}_{tag}.pkl")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh-r-n", type=int, default=80)
    ap.add_argument("--ts-dis", type=float, nargs="*", default=[0.1])
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    _log(f"fig5 START ts_dis={args.ts_dis} mesh_r_n={args.mesh_r_n} plot_only={args.plot_only}")

    para = _build_para(args.mesh_r_n)
    results = []
    for ts in args.ts_dis:
        pkl = _pkl_path(ts, args.mesh_r_n)
        label = "default" if ts is None else f"period={ts}s"
        if os.path.isfile(pkl) and (args.plot_only or args.skip_existing):
            with open(pkl, "rb") as f:
                sol = pickle.load(f)
            _log(f"fig5 load {pkl}")
            results.append((sol, label, 0.0))
            continue
        if args.plot_only:
            _log(f"fig5 missing pickle {pkl}")
            continue
        sol, label, wall = _solve(para, ts, args.mesh_r_n)
        with open(pkl, "wb") as f:
            pickle.dump(sol, f)
        _log(f"fig5 saved {pkl}")
        results.append((sol, label, wall))

    if results:
        tag = "_".join(f"{t}s" for t in args.ts_dis if t is not None) or "default"
        out = os.path.join(OUT_DIR, f"fig5_45C_mesh{args.mesh_r_n}_{tag}_profile.png")
        _plot_profiles(results, out)
    _log("fig5 DONE")


if __name__ == "__main__":
    main()

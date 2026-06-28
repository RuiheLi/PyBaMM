"""
Generate main-text Fig. 8: discharge voltage vs capacity (single vs dual-solvent).

Default: 0.1C, 1C, 2C only (4.5C optional — very slow).

  python run_fig8_discharge_voltage.py --mesh-r-n 60
  python run_fig8_discharge_voltage.py --rates 0.1 1 2 --skip-existing
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

_ROOT = Path(__file__).resolve().parents[2]  # P3_R14
OUT_DIR = _ROOT / "outputs" / "main" / "fig_08_Fig8_discharge_voltage"
STATUS_LOG = os.path.join(OUT_DIR, "run_status.log")
DEFAULT_RATES = [0.1, 1.0, 2.0]


def _log(msg: str) -> None:
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(STATUS_LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def _para(mesh_r_n: int, cross_d: float) -> dict:
    base = {
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
        "Lithium ion EC cross diffusivity [m2.s-1]": [cross_d],
        "EC transference number": ["EC_transference_number_3"],
    }
    stack: list = []
    recursive_scan(stack, base, list(base.keys()), {})
    return stack[0]


def _pkl_path(rate: float, model: str, mesh_r_n: int) -> str:
    return os.path.join(OUT_DIR, f"fig8_{rate}C_{model}_mesh{mesh_r_n}.pkl")


def _solve(para_i: dict, rate: float) -> pybamm.Solution:
    cycle_pack, para_used = Para_init_Dict(para_i)
    mesh_list, model_options = cycle_pack
    model = pybamm.lithium_ion.DFN(options=model_options)
    model = Add_var(para_used, model)
    v_max, v_min = 4.2, 2.5
    ts = 2.0 if rate <= 1 else 0.1
    exp = pybamm.Experiment(
        [
            (
                f"Hold at {v_max} V until C/100",
                f"Discharge at {rate} C until {v_min} V ({ts} second period)",
            )
        ]
    )
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
    _log(f"fig8 solve {rate}C BEGIN")
    sol = sim.solve()
    _log(f"fig8 solve {rate}C OK wall={time.perf_counter() - t0:.1f}s")
    return sol


def _get_or_solve(
    rate: float, model: str, cross: float, mesh_r_n: int, skip_existing: bool
) -> pybamm.Solution:
    pkl = _pkl_path(rate, model, mesh_r_n)
    if skip_existing and os.path.isfile(pkl):
        with open(pkl, "rb") as f:
            _log(f"fig8 load {rate}C {model} from {pkl}")
            return pickle.load(f)
    sol = _solve(_para(mesh_r_n, cross), rate)
    with open(pkl, "wb") as f:
        pickle.dump(sol, f)
    _log(f"fig8 saved {pkl}")
    return sol


def _plot(rates: list[float], mesh_r_n: int) -> str:
    n = len(rates)
    fig, axs = plt.subplots(1, n, figsize=(3.2 * n, 3.2), squeeze=False)
    axs = axs[0]
    colors = {"single": "#1f77b4", "dual": "#d62728"}
    for ax, rate in zip(axs, rates):
        for label, cross, ls in (("single", 0.0, "--"), ("dual", 1e-11, "-")):
            pkl = _pkl_path(rate, label, mesh_r_n)
            if not os.path.isfile(pkl):
                _log(f"fig8 plot skip missing {pkl}")
                continue
            with open(pkl, "rb") as f:
                sol = pickle.load(f)
            step = sol.cycles[0].steps[1]
            cap = (
                step["Discharge capacity [A.h]"].entries
                - step["Discharge capacity [A.h]"].entries[0]
            )
            vol = step["Terminal voltage [V]"].entries
            ax.plot(cap * 1000, vol, ls=ls, color=colors[label], lw=1.5, label=label)
        ax.set_title(f"{rate}C")
        ax.set_xlabel("Discharge capacity [mAh]")
        ax.set_ylabel("Voltage [V]")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(
        "Fig. 8 — discharge voltage (single dashed, dual solid)",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "Fig8_discharge_voltage_profiles.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    _log(f"fig8 plot saved {out}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh-r-n", type=int, default=60)
    ap.add_argument("--rates", type=float, nargs="+", default=DEFAULT_RATES)
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    _log(f"fig8 START rates={args.rates} mesh_r_n={args.mesh_r_n} plot_only={args.plot_only}")

    if not args.plot_only:
        for rate in args.rates:
            for label, cross in (("single", 0.0), ("dual", 1e-11)):
                _get_or_solve(rate, label, cross, args.mesh_r_n, args.skip_existing)

    _plot(args.rates, args.mesh_r_n)
    _log("fig8 DONE")


if __name__ == "__main__":
    main()

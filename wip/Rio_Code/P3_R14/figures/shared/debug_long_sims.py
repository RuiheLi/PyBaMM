"""
命令行调试长耗时仿真：`Reproduce_sol_seg` 与 `OneCycAge` 类协议。

单个用例：
  python debug_long_sims.py reproduce_seg --rate 1 --mesh-r-n 80

完整扫描（多块率 / 参数组合 + **进度与 ETA预估**）:
  python debug_long_sims.py reproduce_seg_full --mesh-r-n 80
  python debug_long_sims.py reproduce_seg_full --mesh-r-n 80 --rate 1.9

单次 aging：`onecyc`、`onecyc --light`。

三次串联 aging（等同 `1EC1DMC_OneCycAge_2C` 主求解格 + **每场打印耗时与 ETA**）:
  python debug_long_sims.py onecyc_full
  python debug_long_sims.py onecyc_full --light

环境变量：`P3_REPRO_SOL_SEG_RATE`、`P3_ONECYC_CRATE`（CLI 优先级更高）。

重定向日志时若不设 `PYTHONUNBUFFERED`，本脚本仍会在导入 PyBaMM 前尽量把 stdout/stderr
设为行缓冲，便于 Tee / `>>` 下即时看到 `_log()` 与时间戳输出。
"""
from __future__ import annotations

import argparse
import io
import os
import sys
import time
from datetime import datetime


def _ensure_line_buffered_stdio() -> None:
    """Line-buffer stdout/stderr when piped/Tee'd (helps live logs on Windows)."""
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if not callable(reconfigure):
            continue
        try:
            reconfigure(line_buffering=True)
        except (OSError, ValueError, TypeError, AttributeError, io.UnsupportedOperation):
            continue


_ensure_line_buffered_stdio()

_SHARED = os.path.dirname(os.path.abspath(__file__))
if _SHARED not in sys.path:
    sys.path.insert(0, _SHARED)
from p3_bootstrap import bootstrap  # noqa: E402

bootstrap()

import pybamm  # noqa: E402
from Fun_P3 import Add_var, Para_init_Dict, recursive_scan  # noqa: E402


def _log(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def _casadi_solver():
    # Notebook uses 200k; high-C segregation can exhaust IDA steps (mxstep/tout failure).
    return pybamm.CasadiSolver(
        return_solution_if_failed_early=True,
        extra_options_setup={"max_num_steps": 1_000_000},
    )


def _eta_str(step_done_1indexed: int, elapsed_batch_s: float, total_steps: int) -> str:
    if step_done_1indexed <= 0 or total_steps <= 1:
        return "ETA=n/a"
    avg = elapsed_batch_s / step_done_1indexed
    rem = avg * max(0, total_steps - step_done_1indexed)
    return f"ETA_remain~{rem / 60:.1f}min ({rem:.0f}s) avg over {step_done_1indexed}/{total_steps} done"


def _build_para_dd_seg(mesh_r_n: int) -> list:
    para_same = {
        "Mesh list": [[20, 10, 20, mesh_r_n, 20]],
        "Para_Set": ["Li2023_ECdrag"],
        "Contact resistance [Ohm]": [6e-3],
        "Initial Neg SOC": [0.8841301667966484],
        "Initial Pos SOC": [0.23552755074598045],
    }
    para_dd_only = {
        "Model option": [
            {
                "SEI": "constant",
                "SEI film resistance": "distributed",
                "SEI porosity change": "true",
                "solvent diffusion": "double spatial consume w refill",
                "electrolyte conductivity": "sol full",
                "contact resistance": "true",
            },
        ],
        "Lithium ion EC cross diffusivity [m2.s-1]": [1e-11],
        "EC transference number": ["EC_transference_number_3"],
    }
    para_dd = {**para_same, **para_dd_only}
    stack: list = []
    recursive_scan(stack, para_dd, list(para_dd.keys()), {})
    return stack


def _solve_reproduce_seg_one(para_i: dict, rate: float, mesh_r_n: int) -> tuple[float, float]:
    cycle_pack, para_used = Para_init_Dict(para_i)
    mesh_list, model_options = cycle_pack
    model = pybamm.lithium_ion.DFN(options=model_options)
    v_max, v_min = 4.2, 2.5
    ts_dis = 0.5 if rate > 1 else 1
    exp = pybamm.Experiment(
        [
            (
                f"Hold at {v_max} V until C/100",
                f"Discharge at {rate} C until {v_min} V ({ts_dis} second period)",
            )
        ]
        * 1
    )
    model = Add_var(para_used, model)
    var_pts = {
        "x_n": mesh_list[0],
        "x_s": mesh_list[1],
        "x_p": mesh_list[2],
        "r_n": mesh_list[3],
        "r_p": mesh_list[4],
    }
    sim = pybamm.Simulation(
        model,
        experiment=exp,
        parameter_values=para_used,
        solver=_casadi_solver(),
        var_pts=var_pts,
    )
    ws = time.perf_counter()
    sol = sim.solve()
    wall = time.perf_counter() - ws
    t_end = float(sol["Time [s]"].entries[-1])
    return t_end, wall


def run_reproduce_seg(rate: float | None, mesh_r_n: int) -> None:
    rate = rate if rate is not None else float(
        os.environ.get("P3_REPRO_SOL_SEG_RATE", "1.5")
    )
    stack = _build_para_dd_seg(mesh_r_n)
    para_i = stack[0]
    wall0 = time.perf_counter()
    _log(
        f"[reproduce_seg] CasADi integrate BEGIN rate={rate}C mesh_r_n={mesh_r_n} "
        f"(solve phase may stay quiet for a long time; normal)"
    )
    t_end, wall = _solve_reproduce_seg_one(para_i, rate, mesh_r_n)
    _log(
        f"[reproduce_seg] OK rate={rate}C mesh_r_n={mesh_r_n} "
        f"t_phys_end[s]={t_end:.3f} wall_s={wall:.1f} (prep+solve total {time.perf_counter() - wall0:.1f}s)"
    )


def run_reproduce_seg_full(rate: float | None, mesh_r_n: int) -> None:
    rate = rate if rate is not None else float(
        os.environ.get("P3_REPRO_SOL_SEG_RATE", "1.5")
    )
    stack = _build_para_dd_seg(mesh_r_n)
    n = len(stack)
    _log(f"===== reproduce_seg_full START n={n} rate={rate}C mesh_r_n={mesh_r_n} =====")
    t_batch = time.perf_counter()
    for i in range(n):
        de = stack[i]["Lithium ion EC cross diffusivity [m2.s-1]"]
        _log(f"--- [{i + 1}/{n}] cross-D={de}: solve BEGIN ---")
        t0 = time.perf_counter()
        t_end, wall = _solve_reproduce_seg_one(stack[i], rate, mesh_r_n)
        elapsed = time.perf_counter() - t0
        batch_elapsed = time.perf_counter() - t_batch
        _log(
            f"--- [{i + 1}/{n}] OK t_phys[s]={t_end:.3f} wall_this={elapsed:.1f}s | "
            f"{_eta_str(i + 1, batch_elapsed, n)} | batch_elapsed={batch_elapsed / 60:.2f} min ---"
        )
    _log(
        f"===== reproduce_seg_full DONE batch wall={time.perf_counter() - t_batch:.1f}s "
        f"({(time.perf_counter() - t_batch) / 60:.2f} min) ====="
    )


def _build_dd_sd_lists(light: bool) -> tuple[list, list]:
    mesh = [8, 4, 8, 60, 16] if light else [10, 5, 10, 100, 20]

    para_same = {
        "Mesh list": [mesh],
        "Para_Set": ["Li2023_ECdrag"],
        "Contact resistance [Ohm]": [6e-3],
        "Initial Neg SOC": [0.8841301667966484],
        "Initial Pos SOC": [0.23552755074598045],
        "Inner SEI lithium interstitial diffusivity [m2.s-1]": [5e-19],
    }
    para_dd_only = {
        "Model option": [
            {
                "SEI": "interstitial-diffusion limited",
                "SEI film resistance": "distributed",
                "SEI porosity change": "true",
                "solvent diffusion": "double spatial consume w refill",
                "electrolyte conductivity": "sol full",
                "contact resistance": "true",
            },
        ],
        "Lithium ion EC cross diffusivity [m2.s-1]": [0, 1e-11],
    }
    para_sd_only = {
        "Model option": [
            {
                "SEI": "interstitial-diffusion limited",
                "SEI film resistance": "distributed",
                "SEI porosity change": "true",
                "solvent diffusion": "single no consume wo refill",
                "electrolyte conductivity": "full",
                "contact resistance": "true",
            },
        ],
    }
    para_dd = {**para_same, **para_dd_only}
    para_sd = {**para_same, **para_sd_only}
    dd_list: list = []
    sd_list: list = []
    recursive_scan(dd_list, para_dd, list(para_dd.keys()), {})
    recursive_scan(sd_list, para_sd, list(para_sd.keys()), {})
    return dd_list, sd_list


def _solve_onecyc_one(para_pick: dict, crate: float, label: str) -> tuple[float, float]:
    cycle_pack, para_used = Para_init_Dict(para_pick)
    mesh_list, model_options = cycle_pack
    model = pybamm.lithium_ion.DFN(options=model_options)
    model = Add_var(para_used, model)
    v_max, v_min = 4.2, 2.5
    exp = pybamm.Experiment(
        [
            (
                f"Hold at {v_max} V until C/20",
                f"Discharge at 1 C until {v_min} V",
                f"Charge at {crate} C until {v_max} V",
                f"Hold at {v_max} V until C/20",
            )
        ]
    )
    var_pts = {
        "x_n": mesh_list[0],
        "x_s": mesh_list[1],
        "x_p": mesh_list[2],
        "r_n": mesh_list[3],
        "r_p": mesh_list[4],
    }
    sim = pybamm.Simulation(
        model,
        experiment=exp,
        parameter_values=para_used,
        solver=_casadi_solver(),
        var_pts=var_pts,
    )
    _log(f"[{label}] CasADi integrating ...")
    t0 = time.perf_counter()
    sol = sim.solve()
    wall = time.perf_counter() - t0
    t_end = float(sol["Time [s]"].entries[-1])
    _log(f"[{label}] integrate DONE wall={wall:.1f}s t_phys[s]={t_end:.3f}")
    return t_end, wall


def run_onecyc(which: str, light: bool) -> None:
    dd_list, sd_list = _build_dd_sd_lists(light)
    crate = 1 if light else float(os.environ.get("P3_ONECYC_CRATE", "2"))

    if which == "ldx":
        para_pick = dd_list[0]
        tag = "DD_LDx"
    elif which == "hdx":
        para_pick = dd_list[1]
        tag = "DD_HDx"
    else:
        para_pick = sd_list[0]
        tag = "SD"

    wall0 = time.perf_counter()
    t_end, wall = _solve_onecyc_one(para_pick, crate, tag)
    _log(
        f"[onecyc:{which}] OK {tag} light={light} Crate={crate} "
        f"t_phys[s]={t_end:.3f} wall={wall:.1f}s (prep+build+total_wall {time.perf_counter() - wall0:.1f}s)"
    )


def run_onecyc_full(light: bool) -> None:
    dd_list, sd_list = _build_dd_sd_lists(light)
    crate = 1 if light else float(os.environ.get("P3_ONECYC_CRATE", "2"))

    jobs = [
        ("DD_HDx Para_DD[1] Dx=1e-11", dd_list[1]),
        ("DD_LDx Para_DD[0] Dx=0", dd_list[0]),
        ("SD_single_solvent", sd_list[0]),
    ]

    mesh = [8, 4, 8, 60, 16] if light else [10, 5, 10, 100, 20]
    _log(
        f"===== onecyc_full START light={light} Crate={crate} mesh={mesh} "
        "order HDx -> LDx -> SD (same as notebook) ====="
    )
    t_batch = time.perf_counter()
    walls: list[float] = []
    for j, (name, para_pick) in enumerate(jobs):
        _log(f"--- [{j + 1}/3] {name}: FULL experiment solve BEGIN ---")
        t_end, wall = _solve_onecyc_one(para_pick, crate, name)
        walls.append(wall)
        batch_elapsed = time.perf_counter() - t_batch
        avg_done = batch_elapsed / (j + 1)
        rem = avg_done * (3 - j - 1)
        _log(
            f"--- [{j + 1}/3] OK t_phys[s]={t_end:.3f} wall_this={wall:.1f}s | "
            f"ETA_remain_crude~{rem / 60:.1f}min (avg of {j + 1} done) | batch_elapsed {(batch_elapsed) / 60:.2f} min ---"
        )
    tw = sum(walls)
    _log(
        f"===== onecyc_full DONE sum_walls={tw:.1f}s ({tw / 60:.2f} min) "
        f"check={sum(walls):.1f}s ====="
    )


def main() -> None:
    p = argparse.ArgumentParser(description="Debug long P3_R14 simulations (+ progress)")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("reproduce_seg", help="单个 Reproduce_sol_seg 放电")
    r.add_argument("--rate", type=float, default=None)
    r.add_argument("--mesh-r-n", type=int, default=80)

    rf = sub.add_parser(
        "reproduce_seg_full",
        help="扫全部 Para_DD 组合（与原 notebook for 循环一致），带 ETA",
    )
    rf.add_argument("--rate", type=float, default=None)
    rf.add_argument("--mesh-r-n", type=int, default=80)

    o = sub.add_parser("onecyc", help="单场 OneCycAge")
    o.add_argument("--which", choices=("ldx", "hdx", "sd"), default="hdx")
    o.add_argument("--light", action="store_true")

    of = sub.add_parser(
        "onecyc_full",
        help="连续三场 aging（对齐 1EC1DMC_OneCycAge_2C），带 ETA",
    )
    of.add_argument("--light", action="store_true")

    args = p.parse_args()
    if args.cmd == "reproduce_seg":
        run_reproduce_seg(args.rate, args.mesh_r_n)
    elif args.cmd == "reproduce_seg_full":
        run_reproduce_seg_full(args.rate, args.mesh_r_n)
    elif args.cmd == "onecyc":
        run_onecyc(args.which, args.light)
    else:
        run_onecyc_full(args.light)


if __name__ == "__main__":
    main()

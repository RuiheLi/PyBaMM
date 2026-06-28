"""Short smoke for the *one-cycle aging* protocol (same experiment string as the OneCycAge notebooks).

Default mode matches `smoke_test_p3_simulation.py` (proven finish time on typical desktops).

Set ``P3_ONE_CYC_STRESS=1`` to run a coarser-mesh 2C charge leg (same chemistry as `1EC1DMC_OneCycAge_2C.ipynb`,
but still much cheaper than three full solves inside the notebook — may still take a long time).

Set ``P3_ONE_CYC_STRESS=2`` for the aggressive coarse grid + 2C configuration used in the first revision of this
script (very stiff; mainly for overnight / cluster).
"""
from __future__ import annotations

import os
import sys
import time

_SHARED = os.path.dirname(os.path.abspath(__file__))
if _SHARED not in sys.path:
    sys.path.insert(0, _SHARED)
from p3_bootstrap import bootstrap  # noqa: E402

bootstrap()

import pybamm  # noqa: E402
from Fun_P3 import Add_var, Para_init_Dict, recursive_scan  # noqa: E402


def _build_para(mesh: list[int]) -> tuple:
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
        "Lithium ion EC cross diffusivity [m2.s-1]": [1e-11],
    }
    para_dd = {**para_same, **para_dd_only}
    stack: list = []
    recursive_scan(stack, para_dd, list(para_dd.keys()), {})
    return Para_init_Dict(stack[0])


def main() -> None:
    stress = os.environ.get("P3_ONE_CYC_STRESS", "0")
    if stress == "0":
        mesh = [10, 5, 10, 100, 20]
        crate = 1
        tag = "quick-default"
    elif stress == "1":
        mesh = [8, 4, 8, 60, 16]
        crate = 2
        tag = "stress-2c-moderate"
    else:
        mesh = [6, 4, 6, 40, 12]
        crate = 2
        tag = "stress-2c-coarse"

    cycle_pack, para_used = _build_para(mesh)
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
        * 1
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
        solver=pybamm.CasadiSolver(return_solution_if_failed_early=True),
        var_pts=var_pts,
    )
    t0 = time.perf_counter()
    sol = sim.solve()
    dt = time.perf_counter() - t0
    print(
        "smoke_short_one_cycle_age OK",
        tag,
        "crate",
        crate,
        "mesh",
        mesh,
        "t_end[s]",
        float(sol["Time [s]"].entries[-1]),
        "wall_s",
        round(dt, 1),
    )


if __name__ == "__main__":
    main()

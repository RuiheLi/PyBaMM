"""Smoke test: Para_init_Dict + one short simulation (no notebook)."""
import os
import sys

_SHARED = os.path.dirname(os.path.abspath(__file__))
if _SHARED not in sys.path:
    sys.path.insert(0, _SHARED)

from p3_bootstrap import bootstrap  # noqa: E402

bootstrap()

import pybamm  # noqa: E402
from Fun_P3 import Add_var, Para_init_Dict, recursive_scan  # noqa: E402

Para_dict_Same = {
    "Mesh list": [[10, 5, 10, 100, 20]],
    "Para_Set": ["Li2023_ECdrag"],
    "Contact resistance [Ohm]": [6e-3],
    "Initial Neg SOC": [0.8841301667966484],
    "Initial Pos SOC": [0.23552755074598045],
    "Inner SEI lithium interstitial diffusivity [m2.s-1]": [5e-19],
}
Para_dict_DD_ONLY = {
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
Para_dict_DD = {**Para_dict_Same, **Para_dict_DD_ONLY}
Para_DD = []
recursive_scan(Para_DD, Para_dict_DD, list(Para_dict_DD.keys()), {})

CyclePack, para_used = Para_init_Dict(Para_DD[0])
Mesh_list, model_options = CyclePack
model = pybamm.lithium_ion.DFN(options=model_options)
model = Add_var(para_used, model)
Crate = 1
V_max, V_min = 4.2, 2.5
Exp_1 = pybamm.Experiment(
    [
        (
            f"Hold at {V_max} V until C/20",
            f"Discharge at 1 C until {V_min} V",
            f"Charge at {Crate} C until {V_max} V",
            f"Hold at {V_max} V until C/20",
        )
    ]
    * 1
)
var_pts = {
    "x_n": Mesh_list[0],
    "x_s": Mesh_list[1],
    "x_p": Mesh_list[2],
    "r_n": Mesh_list[3],
    "r_p": Mesh_list[4],
}
sim = pybamm.Simulation(
    model,
    experiment=Exp_1,
    parameter_values=para_used,
    solver=pybamm.CasadiSolver(return_solution_if_failed_early=True),
    var_pts=var_pts,
)
sol = sim.solve()
print("OK", sol["Time [s]"].entries[-1])

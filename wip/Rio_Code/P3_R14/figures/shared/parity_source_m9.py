"""
Parity reference for M9: end-to-end P3_R14 protocol on the source repo.

Mirrors ``smoke_test_p3_simulation.py``:
- DFN with Li2023_ECdrag,
- Mesh [10, 5, 10, 100, 20],
- ``SEI = "interstitial-diffusion limited"``,
- ``SEI film resistance = "distributed"``,
- ``SEI porosity change = "true"``,
- ``solvent diffusion = "double spatial consume w refill"``,
- ``electrolyte conductivity = "sol full"``,
- ``contact resistance = "true"``,
- 4-stage experiment: 4.2 V hold -> 1C dis -> 1C ch -> 4.2 V hold.

Saves t / V / c_e_avg / c_EC_avg / Q_sei to ``parity_source_m9.npz``.
"""

import os
import sys
import time

import numpy as np

import pybamm


def main():
    print("source pybamm at", os.path.dirname(pybamm.__file__))

    options = {
        "SEI": "interstitial-diffusion limited",
        "SEI film resistance": "distributed",
        "SEI porosity change": "true",
        "solvent diffusion": "double spatial consume w refill",
        "electrolyte conductivity": "sol full",
        "contact resistance": "true",
    }
    model = pybamm.lithium_ion.DFN(options=options)

    p = pybamm.ParameterValues("Li2023_ECdrag")

    p.update(
        {
            "Contact resistance [Ohm]": 6e-3,
            "Inner SEI lithium interstitial diffusivity [m2.s-1]": 5e-19,
            "Initial concentration in negative electrode [mol.m-3]":
                0.8841301667966484
                * float(p["Maximum concentration in negative electrode [mol.m-3]"]),
            "Initial concentration in positive electrode [mol.m-3]":
                0.23552755074598045
                * float(p["Maximum concentration in positive electrode [mol.m-3]"]),
            "Lithium ion EC cross diffusivity [m2.s-1]": 1e-11,
        },
        check_already_exists=False,
    )

    Crate = 1
    V_max, V_min = 4.2, 2.5
    Exp_1 = pybamm.Experiment(
        [
            (
                f"Hold at {V_max} V until C/20",
                f"Discharge at {Crate} C until {V_min} V",
                f"Charge at {Crate} C until {V_max} V",
                f"Hold at {V_max} V until C/20",
            )
        ]
        * 1
    )

    var_pts = {"x_n": 10, "x_s": 5, "x_p": 10, "r_n": 100, "r_p": 20}

    sim = pybamm.Simulation(
        model,
        experiment=Exp_1,
        parameter_values=p,
        var_pts=var_pts,
        solver=pybamm.CasadiSolver(return_solution_if_failed_early=True),
    )

    t0 = time.time()
    sol = sim.solve()
    t1 = time.time()

    t_phys = sol["Time [s]"].entries
    V = sol["Terminal voltage [V]"].entries
    c_e_avg = sol["X-averaged electrolyte concentration"].entries
    c_e_typ = float(p["Typical electrolyte concentration [mol.m-3]"])
    c_e_avg_phys = c_e_avg * c_e_typ
    c_EC_avg = sol["X-averaged EC concentration"].entries
    c_ec_typ = float(p["Typical EC concentration [mol.m-3]"])
    c_EC_avg_phys = c_EC_avg * c_ec_typ
    Q_sei = sol["Loss of lithium to SEI [mol]"].entries

    print(f"source M9 solve {t1 - t0:.2f}s, t_end={float(t_phys[-1]):.2f}s")
    print(f"  V[0]={float(V[0]):.4f}, V[-1]={float(V[-1]):.4f}, "
          f"V_min={float(V.min()):.4f}, V_max={float(V.max()):.4f}")
    print(f"  c_e_avg [{float(c_e_avg_phys.min()):.2f}, "
          f"{float(c_e_avg_phys.max()):.2f}] mol/m3")
    print(f"  c_EC_avg [{float(c_EC_avg_phys.min()):.2f}, "
          f"{float(c_EC_avg_phys.max()):.2f}] mol/m3")
    print(f"  Q_sei[-1] = {float(Q_sei[-1]):.6e} mol")

    out_path = os.path.join(os.path.dirname(__file__), "parity_source_m9.npz")
    np.savez(
        out_path,
        t=t_phys,
        V=V,
        c_e_avg=c_e_avg_phys,
        c_EC_avg=c_EC_avg_phys,
        Q_sei=Q_sei,
        c_e_typ=c_e_typ,
        c_ec_typ=c_ec_typ,
    )
    print(f"saved {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

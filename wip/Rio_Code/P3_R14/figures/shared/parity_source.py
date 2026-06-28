"""
Parity reference: short 1C discharge on the source ECdrag2 repo with
sol-full conductivity but a *constant* c_EC (= c_ec_init = 6209.49).

This isolates the LJP-extended Ohm's law (the M3 port deliverable). The
matching port-side run uses ``solvent diffusion = "constant"`` so both
solvers see identical dLJP physics with c_EC fixed at its initial value.

Saves ``parity_source.npz`` next to this script.

Run with cwd = D:/LRHWork/Model_RH/PyBaMM-ECdrag2 and
PYTHONPATH = D:/LRHWork/Model_RH/PyBaMM-ECdrag2 so the local pybamm
(the ECdrag2 fork) is picked up.
"""

import os
import sys
import time

import numpy as np

import pybamm


def main():
    print("source pybamm at", os.path.dirname(pybamm.__file__))

    options = {
        "SEI": "none",
        "lithium plating": "none",
        "solvent diffusion": "single no consume wo refill",
        "electrolyte conductivity": "sol full",
    }
    model = pybamm.lithium_ion.DFN(options=options)

    p = pybamm.ParameterValues("Li2023_ECdrag")
    nominal_capacity = float(p["Nominal cell capacity [A.h]"])
    p.update(
        {
            "Contact resistance [Ohm]": 0.0,
            "Initial concentration in negative electrode [mol.m-3]":
                0.8841301667966484
                * float(p["Maximum concentration in negative electrode [mol.m-3]"]),
            "Initial concentration in positive electrode [mol.m-3]":
                0.23552755074598045
                * float(p["Maximum concentration in positive electrode [mol.m-3]"]),
            "Current function [A]": nominal_capacity,
        },
        check_already_exists=False,
    )

    var_pts = {"x_n": 10, "x_s": 5, "x_p": 10, "r_n": 30, "r_p": 20}

    sim = pybamm.Simulation(
        model,
        parameter_values=p,
        var_pts=var_pts,
        solver=pybamm.CasadiSolver(),
    )
    t_grid = np.linspace(0, 60, 31)
    t0 = time.time()
    try:
        sol = sim.solve(t_grid)
    except Exception as e:
        print("source solve EXCEPTION:", repr(e))
        raise
    t1 = time.time()

    V = sol["Terminal voltage [V]"].entries
    c_e_avg = sol["X-averaged electrolyte concentration"].entries
    c_e_typ = float(p["Typical electrolyte concentration [mol.m-3]"])
    c_e_avg_phys = c_e_avg * c_e_typ
    c_EC_avg = sol["X-averaged EC concentration"].entries
    c_ec_typ = float(p["Typical EC concentration [mol.m-3]"])
    c_EC_avg_phys = c_EC_avg * c_ec_typ
    t = sol.t
    t_phys = sol["Time [s]"].entries

    print(f"source solve {t1 - t0:.2f}s")
    print(f"  sol.t shape={t.shape}, min={float(t.min()):.4f}, max={float(t.max()):.4f}")
    print(f"  Time [s] shape={t_phys.shape}, min={float(t_phys.min()):.2f}, max={float(t_phys.max()):.2f}")
    print(f"  V series: V[0]={float(V[0]):.4f}, V[-1]={float(V[-1]):.4f}, n={len(V)}")
    print(f"  c_e_avg_end = {float(c_e_avg_phys[-1]):.2f} mol/m3")
    print(f"  c_EC_avg_end = {float(c_EC_avg_phys[-1]):.2f} mol/m3")
    print(f"  c_e_typ = {c_e_typ}, c_ec_typ = {c_ec_typ}")

    out_path = os.path.join(os.path.dirname(__file__), "parity_source.npz")
    np.savez(
        out_path,
        t=t_phys,
        V=V,
        c_e_avg=c_e_avg_phys,
        c_EC_avg=c_EC_avg_phys,
        c_e_typ=c_e_typ,
        c_ec_typ=c_ec_typ,
        nominal_capacity=nominal_capacity,
    )
    print(f"saved {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

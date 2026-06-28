"""
Parity reference for M8: 0-300 s 1C discharge with sol-full conductivity AND
``solvent diffusion = "double spatial consume w refill"`` AND SEI kinetics
ON ("reaction limited" with the Li2023_ECdrag SEI parameters).

Adds the SEI consumption + refill physics on top of M7's migration setup.

Saves ``parity_source_m8.npz`` next to this script.
"""

import os
import sys
import time

import numpy as np

import pybamm


def main():
    print("source pybamm at", os.path.dirname(pybamm.__file__))

    options = {
        "SEI": "reaction limited",
        "lithium plating": "none",
        "solvent diffusion": "double spatial consume w refill",
        "electrolyte conductivity": "sol full",
    }
    model = pybamm.lithium_ion.DFN(options=options)

    p = pybamm.ParameterValues("Li2023_ECdrag")
    nominal_capacity = float(p["Nominal cell capacity [A.h]"])
    p.update(
        {
            "Contact resistance [Ohm]": 0.0,
            "EC diffusivity in electrolyte [m2.s-1]":
                lambda c_e, c_EC, T: 5e-10 + 0 * c_e + 0 * c_EC + 0 * T,
            "EC Lithium ion cross diffusivity [m2.s-1]":
                lambda c_e, c_EC, T: 5e-12 + 0 * c_e + 0 * c_EC + 0 * T,
            "Lithium ion EC cross diffusivity [m2.s-1]": 0,
            # Xi = 3 (default), c_T closure = Fun_c_T (default).
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
    t_grid = np.linspace(0, 300, 31)
    t0 = time.time()
    sol = sim.solve(t_grid)
    t1 = time.time()

    V = sol["Terminal voltage [V]"].entries
    c_e_avg = sol["X-averaged electrolyte concentration"].entries
    c_e_typ = float(p["Typical electrolyte concentration [mol.m-3]"])
    c_e_avg_phys = c_e_avg * c_e_typ
    c_EC_avg = sol["X-averaged EC concentration"].entries
    c_ec_typ = float(p["Typical EC concentration [mol.m-3]"])
    c_EC_avg_phys = c_EC_avg * c_ec_typ
    Q_sei = sol["Loss of lithium to SEI [mol]"].entries
    t_phys = sol["Time [s]"].entries

    print(f"source M8 solve {t1 - t0:.2f}s")
    print(f"  Time [s] [{float(t_phys.min()):.2f}, {float(t_phys.max()):.2f}]")
    print(f"  V[0]={float(V[0]):.4f}, V[-1]={float(V[-1]):.4f}")
    print(f"  c_e_avg_end = {float(c_e_avg_phys[-1]):.2f} mol/m3")
    print(f"  c_EC_avg_end = {float(c_EC_avg_phys[-1]):.2f} mol/m3")
    print(f"  Q_sei[-1] = {float(Q_sei[-1]):.6e} mol")

    out_path = os.path.join(os.path.dirname(__file__), "parity_source_m8.npz")
    np.savez(
        out_path,
        t=t_phys,
        V=V,
        c_e_avg=c_e_avg_phys,
        c_EC_avg=c_EC_avg_phys,
        Q_sei=Q_sei,
        c_e_typ=c_e_typ,
        c_ec_typ=c_ec_typ,
        nominal_capacity=nominal_capacity,
    )
    print(f"saved {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

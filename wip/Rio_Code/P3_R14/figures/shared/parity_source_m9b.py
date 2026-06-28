"""M9b: same as M9 but with SEI porosity change OFF, to isolate the EC drag
parity from the porosity-change feedback loop."""
import os
import sys
import time

import numpy as np

import pybamm


def main():
    options = {
        "SEI": "interstitial-diffusion limited",
        "SEI film resistance": "distributed",
        "SEI porosity change": "false",
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
    Crate, V_max, V_min = 1, 4.2, 2.5
    Exp_1 = pybamm.Experiment(
        [
            (
                f"Hold at {V_max} V until C/20",
                f"Discharge at {Crate} C until {V_min} V",
                f"Charge at {Crate} C until {V_max} V",
                f"Hold at {V_max} V until C/20",
            )
        ]
    )
    var_pts = {"x_n": 10, "x_s": 5, "x_p": 10, "r_n": 100, "r_p": 20}
    sim = pybamm.Simulation(
        model, experiment=Exp_1, parameter_values=p, var_pts=var_pts,
        solver=pybamm.CasadiSolver(return_solution_if_failed_early=True),
    )
    t0 = time.time()
    sol = sim.solve()
    print(f"source M9b solve {time.time()-t0:.2f}s, t_end={float(sol['Time [s]'].entries[-1]):.2f}s")

    t_phys = sol["Time [s]"].entries
    V = sol["Terminal voltage [V]"].entries
    c_e_typ = float(p["Typical electrolyte concentration [mol.m-3]"])
    c_ec_typ = float(p["Typical EC concentration [mol.m-3]"])
    c_e_avg_phys = sol["X-averaged electrolyte concentration"].entries * c_e_typ
    c_EC_avg_phys = sol["X-averaged EC concentration"].entries * c_ec_typ
    Q_sei = sol["Loss of lithium to SEI [mol]"].entries
    print(f"  V[0]={float(V[0]):.4f}, V_min={float(V.min()):.4f}, V_max={float(V.max()):.4f}")
    print(f"  c_e_avg [{float(c_e_avg_phys.min()):.2f}, {float(c_e_avg_phys.max()):.2f}]")
    print(f"  c_EC_avg [{float(c_EC_avg_phys.min()):.2f}, {float(c_EC_avg_phys.max()):.2f}]")
    print(f"  Q_sei[-1]={float(Q_sei[-1]):.4e}")

    out_path = os.path.join(os.path.dirname(__file__), "parity_source_m9b.npz")
    np.savez(
        out_path, t=t_phys, V=V,
        c_e_avg=c_e_avg_phys, c_EC_avg=c_EC_avg_phys, Q_sei=Q_sei,
        c_e_typ=c_e_typ, c_ec_typ=c_ec_typ,
    )
    print(f"saved {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

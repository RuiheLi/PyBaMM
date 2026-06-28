"""
Parity reference for M7: 0-60 s 1C discharge with sol-full conductivity AND
``solvent diffusion = "double no consume wo refill"`` AND the EC migration
term active (``Xi = 3.0`` constant from ``EC_transference_number_3``), while
SEI is still off so we can run side-by-side with the port without the SEI
hooks.

Differences vs the M6 source script:
- ``Xi`` left at the default (constant 3.0) -> migration term ON.
- Subclass ``DoubleNoSEIWithMigration`` keeps the migration term but stubs
  ``Loss of lithium to SEI [mol]`` so SEI=none works.

Saves ``parity_source_m7.npz`` next to this script.
"""

import os
import sys
import time

import numpy as np

import pybamm
from pybamm.models.submodels.solvent_diffusion.Double_NoConsume_wo_refill import (
    Double_NoConsume_wo_refill,
)


class DoubleNoSEIWithMigration(Double_NoConsume_wo_refill):
    """Same as :class:`Double_NoConsume_wo_refill` but with the SEI hooks
    removed so it can be used with ``SEI = "none"``.

    The EC migration term (``C_RA_typ * c_EC / c_tot * Xi * i_e``) is left
    active so the source-side N_EC is

        N_EC = -tor * D_ec * grad(c_EC)
             - (tau_ec / tau_cross / gamma_e_ec) * tor * D_ec_Li_cross * grad(c_e)
             + C_RA_typ * c_EC / c_tot * Xi * i_e

    which non-dimensionally matches the port's
    ``solvent diffusion = "full with migration"`` mode.
    """

    def get_coupled_variables(self, variables):
        c_EC_dict = {}
        for domain in self.options.whole_cell_domains:
            Domain = domain.capitalize()
            eps_k = variables[f"{Domain} porosity"]
            eps_c_EC_k = variables[f"{Domain} porosity times EC concentration"]
            c_EC_k = eps_c_EC_k / eps_k
            c_EC_dict[domain] = c_EC_k
        variables.update(self._get_standard_EC_concentration_variables(c_EC_dict))

        eps_c_EC = variables["Porosity times EC concentration"]
        c_e = variables["Electrolyte concentration"]
        tor = variables["Electrolyte transport efficiency"]
        i_e = variables["Electrolyte current density"]
        T = variables["Cell temperature"]
        c_EC = variables["EC concentration"]

        param = self.param

        N_EC_diffusion = -tor * param.D_ec(c_e, c_EC, T) * pybamm.grad(c_EC)
        N_cross_diffusion = -(
            param.tau_ec / param.tau_cross / param.gamma_e_ec
            * tor * param.D_ec_Li_cross(c_e, c_EC, T) * pybamm.grad(c_e)
        )
        N_EC_migration = (
            param.C_RA_typ * c_EC / param.c_tot(c_e, c_EC, T)
            * param.Xi(c_e, c_EC, T) * i_e
        )

        N_EC = N_EC_diffusion + N_cross_diffusion + N_EC_migration

        sign_2_n = pybamm.FullBroadcast(
            pybamm.Scalar(0), "negative electrode",
            auxiliary_domains={"secondary": "current collector"})
        sign_2_s = pybamm.FullBroadcast(
            pybamm.Scalar(0), "separator",
            auxiliary_domains={"secondary": "current collector"})
        sign_2_p = pybamm.FullBroadcast(
            pybamm.Scalar(0), "positive electrode",
            auxiliary_domains={"secondary": "current collector"})
        sign_2 = pybamm.concatenation(sign_2_n, sign_2_s, sign_2_p)

        source_terms_ec = sign_2
        source_terms_refill = sign_2

        variables.update(self._get_standard_EC_flux_variables(
            N_EC, N_EC_diffusion, N_EC_migration, N_cross_diffusion,
            source_terms_ec, source_terms_refill,
        ))
        variables.update(
            self._get_total_EC_concentration_electrolyte(eps_c_EC, pybamm.Scalar(0))
        )
        return variables


def main():
    print("source pybamm at", os.path.dirname(pybamm.__file__))

    options = {
        "SEI": "none",
        "lithium plating": "none",
        "solvent diffusion": "double no consume wo refill",
        "electrolyte conductivity": "sol full",
    }
    model = pybamm.lithium_ion.DFN(options=options)
    model.submodels["solvent diffusion"] = DoubleNoSEIWithMigration(
        model.param, model.options
    )

    p = pybamm.ParameterValues("Li2023_ECdrag")
    nominal_capacity = float(p["Nominal cell capacity [A.h]"])
    p.update(
        {
            "Contact resistance [Ohm]": 0.0,
            # Force constant transport functions so source <-> port match.
            "EC diffusivity in electrolyte [m2.s-1]":
                lambda c_e, c_EC, T: 5e-10 + 0 * c_e + 0 * c_EC + 0 * T,
            "EC Lithium ion cross diffusivity [m2.s-1]":
                lambda c_e, c_EC, T: 5e-12 + 0 * c_e + 0 * c_EC + 0 * T,
            "Lithium ion EC cross diffusivity [m2.s-1]": 0,
            # NOTE: Xi default = EC_transference_number_3 -> 3.0, kept ON.
            # NOTE: Total concentration default = Fun_c_T (mass-balance closure).
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
    sol = sim.solve(t_grid)
    t1 = time.time()

    V = sol["Terminal voltage [V]"].entries
    c_e_avg = sol["X-averaged electrolyte concentration"].entries
    c_e_typ = float(p["Typical electrolyte concentration [mol.m-3]"])
    c_e_avg_phys = c_e_avg * c_e_typ
    c_EC_avg = sol["X-averaged EC concentration"].entries
    c_ec_typ = float(p["Typical EC concentration [mol.m-3]"])
    c_EC_avg_phys = c_EC_avg * c_ec_typ
    t_phys = sol["Time [s]"].entries

    print(f"source M7 solve {t1 - t0:.2f}s")
    print(f"  Time [s] [{float(t_phys.min()):.2f}, {float(t_phys.max()):.2f}]")
    print(f"  V[0]={float(V[0]):.4f}, V[-1]={float(V[-1]):.4f}")
    print(f"  c_e_avg_end = {float(c_e_avg_phys[-1]):.2f} mol/m3")
    print(f"  c_EC_avg_end = {float(c_EC_avg_phys[-1]):.2f} mol/m3")

    out_path = os.path.join(os.path.dirname(__file__), "parity_source_m7.npz")
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

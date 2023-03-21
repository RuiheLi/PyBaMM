#
# Class for electrolyte diffusion employing stefan-maxwell
#
import pybamm
import numpy as np
from .base_electrolyte_diffusion import BaseElectrolyteDiffusion


class Full(BaseElectrolyteDiffusion):
    """Class for conservation of mass in the electrolyte employing the
    Stefan-Maxwell constitutive equations. (Full refers to unreduced by
    asymptotic methods)

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    options : dict, optional
        A dictionary of options to be passed to the model.

    **Extends:** :class:`pybamm.electrolyte_diffusion.BaseElectrolyteDiffusion`
    """

    def __init__(self, param, options=None):
        super().__init__(param, options)

    def get_fundamental_variables(self):
        eps_c_e_dict = {}
        for domain in self.options.whole_cell_domains:
            Domain = domain.capitalize()
            eps_c_e_k = pybamm.Variable(
                f"{Domain} porosity times concentration",
                domain=domain,
                auxiliary_domains={"secondary": "current collector"},
                bounds=(0, np.inf),
            )
            eps_c_e_k.print_name = f"eps_c_e_{domain[0]}"
            eps_c_e_dict[domain] = eps_c_e_k

        variables = self._get_standard_porosity_times_concentration_variables(
            eps_c_e_dict
        )

        return variables

    def get_coupled_variables(self, variables):

        c_e_dict = {}
        for domain in self.options.whole_cell_domains:
            Domain = domain.capitalize()
            eps_k = variables[f"{Domain} porosity"]
            eps_c_e_k = variables[f"{Domain} porosity times concentration"]
            c_e_k = eps_c_e_k / eps_k
            c_e_dict[domain] = c_e_k
        variables.update(self._get_standard_concentration_variables(c_e_dict))


        eps_c_e = variables["Porosity times concentration"]
        c_e = variables["Electrolyte concentration"]
        tor = variables["Electrolyte transport efficiency"]
        i_e = variables["Electrolyte current density"]
        v_box = variables["Volume-averaged velocity"]
        T = variables["Cell temperature"]
        c_EC  = variables["EC concentration"]

        param = self.param

        N_e_diffusion =   -   tor * param.D_e(c_e, c_EC,T) * pybamm.grad(c_e)
        N_e_migration = param.C_e * param.t_plus(c_e,c_EC, T) * i_e / param.gamma_e
        N_e_convection = param.C_e * c_e * v_box

        if self.options["solvent diffusion"] in ["single no consume wo refill", 
            "single spatial consume w refill","single spatial consume wo refill"]:
            N_cross_diffusion = pybamm.FullBroadcastToEdges(0,
                [domain for domain in self.options.whole_cell_domains],
                "current collector",)
          
        elif self.options["solvent diffusion"] in ["double spatial consume w refill",
            "double spatial consume wo refill","double no consume wo refill"]:
            N_cross_diffusion =  -   (
                param.e_ratio_Rio * param.gamma_e_ec_Rio * 
                param.tau_diffusion_e / param.tau_cross_Rio * 
                tor * param.D_ec_Li_cross(c_e, c_EC,T) * pybamm.grad(c_EC) )
        
        N_e = N_e_diffusion + N_cross_diffusion + N_e_migration + N_e_convection

        # Mark Ruihe block start
        sign_2_n = pybamm.FullBroadcast(
                pybamm.Scalar(0), "negative electrode", 
                auxiliary_domains={"secondary": "current collector"}
            )
        sign_2_s = pybamm.FullBroadcast(
                pybamm.Scalar(0), "separator", 
                auxiliary_domains={"secondary": "current collector"}
            )
        sign_2_p = pybamm.FullBroadcast(
                pybamm.Scalar(0), "positive electrode", 
                auxiliary_domains={"secondary": "current collector"}
            )
        sign_2 = pybamm.concatenation(sign_2_n, sign_2_s, sign_2_p )

        a_p = variables["Positive electrode surface area to volume ratio"]
        a_n = variables["Negative electrode surface area to volume ratio"]
        zero_s = pybamm.FullBroadcast(0, "separator", "current collector")
        a = pybamm.concatenation(a_n, zero_s, a_p)

        j_inner =  variables["Inner SEI interfacial current density"]
        j_outer =  variables["Outer SEI interfacial current density"]
        j_SEI = j_inner + j_outer
        j_sign_SEI = pybamm.concatenation(j_SEI, sign_2_s, sign_2_p )

        ratio_sei_li = - 1/param.z_sei  ;# change to 1 for now , initially is 0.5
        ratio_ec_li  = 1 ; 

        sum_s_j = variables["Sum of electrolyte reaction source terms"]
        sum_s_j.print_name = "a"
        source_terms = sum_s_j / self.param.gamma_e
        if self.options["solvent diffusion"] in ["single no consume wo refill",
            "double no consume wo refill",   # Ruihe add 230320 
            "double spatial consume wo refill","single spatial consume wo refill" ]:
            source_terms_refill = sign_2
        # Mark Ruihe update 221021 - ignore Li+ (de-)intercalation to 
        # avoid differences between w and wo refill when there are no SEI 
        elif self.options["solvent diffusion"] in ["double spatial consume w refill",
            "single spatial consume w refill"]:
            source_terms_refill = - (
                    a * j_sign_SEI / param.gamma_e * param.Vmolar_Li * param.c_e_init_dimensional  +
                    a * j_sign_SEI / param.gamma_e * param.c_e_init_dimensional * (
                        param.Vmolar_ec*ratio_ec_li +param.Vmolar_CH2OCO2Li2*ratio_sei_li
            ))
        
        variables.update(self._get_standard_flux_variables(
            N_e,N_e_diffusion,N_e_migration,N_cross_diffusion,
            source_terms,source_terms_refill))


        return variables

    def set_rhs(self, variables):

        param = self.param

        eps_c_e = variables["Porosity times concentration"]
        c_e = variables["Electrolyte concentration"]
        c_EC  = variables["EC concentration"]
        N_e   = variables["Li+ flux"]
        div_Vbox = variables["Transverse volume-averaged acceleration"]

        source_terms = variables["Li+ source term"]
        source_terms_refill = variables["Li+ source term refill"]


        self.rhs = {
            eps_c_e: -pybamm.div(N_e) / param.C_e + source_terms 
            - c_e * div_Vbox
            # source term due to replenishment
            + source_terms_refill
        } 


    def set_initial_conditions(self, variables):

        eps_c_e = variables["Porosity times concentration"]

        self.initial_conditions = {
            eps_c_e: self.param.epsilon_init * self.param.c_e_init
        }

    def set_boundary_conditions(self, variables):
        param = self.param
        c_e = variables["Electrolyte concentration"]
        c_EC  = variables["EC concentration"]
        T = variables["Cell temperature"]
        tor = variables["Electrolyte transport efficiency"]
        i_boundary_cc = variables["Current collector current density"]

        def flux_bc(side):
            # returns the flux at a separator/electrode interface
            # assuming v_box = 0 for now
            return (
                pybamm.boundary_value(
                    -(1 - param.t_plus(c_e, T))
                    / (tor * param.gamma_e * param.D_e(c_e, T)),
                    side,
                )
                * i_boundary_cc
                * param.C_e
            )

        if self.options.whole_cell_domains[0] == "negative electrode":
            # left bc at anode/current collector interface
            lbc = pybamm.Scalar(0)
        elif self.options.whole_cell_domains[0] == "separator":
            # left bc at anode/separator interface
            lbc = flux_bc("left")
        if self.options.whole_cell_domains[-1] == "positive electrode":
            # right bc at cathode/current collector interface
            rbc = pybamm.Scalar(0)
        # elif self.options.whole_cell_domains[-1] == "separator":
        #     # right bc at separator/cathode interface
        #     rbc = flux_bc("right")

        self.boundary_conditions = {
            c_e: {"left": (lbc, "Neumann"), "right": (rbc, "Neumann")},
        }
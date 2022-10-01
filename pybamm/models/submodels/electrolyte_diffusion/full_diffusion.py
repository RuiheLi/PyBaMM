#
# Class for electrolyte diffusion employing stefan-maxwell
#
import pybamm

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
        if self.half_cell:
            eps_c_e_n = None
        else:
            eps_c_e_n = pybamm.standard_variables.eps_c_e_n
        eps_c_e_s = pybamm.standard_variables.eps_c_e_s
        eps_c_e_p = pybamm.standard_variables.eps_c_e_p

        variables = self._get_standard_porosity_times_concentration_variables(
            eps_c_e_n, eps_c_e_s, eps_c_e_p
        )

        return variables

    def get_coupled_variables(self, variables):

        if self.half_cell:
            c_e_n = None
        else:
            eps_n = variables["Negative electrode porosity"]
            eps_c_e_n = variables["Negative electrode porosity times concentration"]
            c_e_n = eps_c_e_n / eps_n

        eps_s = variables["Separator porosity"]
        eps_p = variables["Positive electrode porosity"]
        eps_c_e_s = variables["Separator porosity times concentration"]
        eps_c_e_p = variables["Positive electrode porosity times concentration"]
        c_e_s = eps_c_e_s / eps_s
        c_e_p = eps_c_e_p / eps_p

        variables.update(
            self._get_standard_concentration_variables(c_e_n, c_e_s, c_e_p)
        )

        eps_c_e = variables["Porosity times concentration"]
        c_e = variables["Electrolyte concentration"]
        tor = variables["Electrolyte transport efficiency"]
        i_e = variables["Electrolyte current density"]
        v_box = variables["Volume-averaged velocity"]
        T = variables["Cell temperature"]
        c_EC  = variables["EC concentration"]

        param = self.param

        N_e_diffusion = -tor * param.D_e(c_e, c_EC,T) * pybamm.grad(c_e)
        N_cross_diffusion = -(
            param.e_ratio_Rio * param.gamma_e_ec_Rio * 
            param.tau_diffusion_e / param.tau_cross_Rio * 
            tor * param.D_ec_Li_cross * pybamm.grad(c_EC) )
        N_e_migration = param.C_e * param.t_plus(c_e,c_EC, T) * i_e / param.gamma_e
        N_e_convection = param.C_e * c_e * v_box

        if self.options["solvent diffusion"] == "none":
            N_e = N_e_diffusion + N_e_migration + N_e_convection
        elif self.options["solvent diffusion"] == "EC":
            N_e = N_e_diffusion + N_cross_diffusion + N_e_migration

        variables.update(self._get_standard_flux_variables(N_e))
        variables.update(self._get_total_concentration_electrolyte(eps_c_e))

        return variables

    def set_rhs(self, variables):

        param = self.param


        # Mark Ruihe block start
        sign_2_n = pybamm.FullBroadcast(
                pybamm.Scalar(1), "negative electrode", 
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

        eps_c_e = variables["Porosity times concentration"]
        c_e = variables["Electrolyte concentration"]
        c_EC  = variables["EC concentration"]
        N_e = variables["Electrolyte flux"]
        tor = variables["Electrolyte transport efficiency"]
        T = variables["Cell temperature"]
        j_inner =  variables["Inner SEI interfacial current density"]
        j_outer =  variables["Outer SEI interfacial current density"]
        j_SEI = j_inner + j_outer
        j_sign_SEI = pybamm.concatenation(j_SEI, sign_2_s, sign_2_p )
        div_Vbox = variables["Transverse volume-averaged acceleration"]

        sum_s_j = variables["Sum of electrolyte reaction source terms"]
        sum_s_j.print_name = "a"
        source_terms = sum_s_j / self.param.gamma_e

        ratio_sei_li = -0.5 ; # change to 1 for now , initially is 0.5
        ratio_ec_li  = 1 ; 

        if self.options["solvent diffusion"] == "none":
            self.rhs = {
                eps_c_e: -pybamm.div(N_e) / param.C_e + source_terms - c_e * div_Vbox
            }
        elif self.options["solvent diffusion"] == "EC w refill":
            self.rhs = {
            eps_c_e: -pybamm.div(N_e) / param.C_e + source_terms 
            - c_e * div_Vbox
            # source term due to replenishment
            - (
                source_terms * param.Vmolar_Li * param.c_e_init_dimensional  +
                a * j_sign_SEI / param.gamma_e * param.c_e_init_dimensional * (
                    param.Vmolar_ec*ratio_ec_li +param.Vmolar_CH2OCO2Li2*ratio_sei_li
                )
            )
            #-  (    
            #    param.c_e_init_dimensional / param.gamma_e * 
            #    (param.Vmolar_ec - 
            #    ratio_sei_li*param.Vmolar_CH2OCO2Li2  # + param.Vmolar_Li # ignore volume of lithium
            #    ) * a * j_sign_SEI     ) 
            } 
        elif self.options["solvent diffusion"] == "EC wo refill":
            self.rhs = {
            eps_c_e: -pybamm.div(N_e) / param.C_e + source_terms 
            - c_e * div_Vbox
            } 
        
        # Mark Ruihe block start

    def set_initial_conditions(self, variables):

        eps_c_e = variables["Porosity times concentration"]

        self.initial_conditions = {
            eps_c_e: self.param.epsilon_init * self.param.c_e_init
        }

    def set_boundary_conditions(self, variables):
        param = self.param
        c_e = variables["Electrolyte concentration"]
        c_EC  = variables["EC concentration"]

        if self.half_cell:
            # left bc at anode/separator interface
            # assuming v_box = 0 for now
            T = variables["Cell temperature"]
            tor = variables["Electrolyte transport efficiency"]
            i_boundary_cc = variables["Current collector current density"]
            lbc = (
                pybamm.boundary_value(
                    -(1 - param.t_plus(c_e,c_EC, T))
                    / (tor * param.gamma_e * param.D_e(c_e,c_EC, T)),
                    "left",
                )
                * i_boundary_cc
                * param.C_e
            )
        else:
            # left bc at anode/current collector interface
            lbc = pybamm.Scalar(0)

        self.boundary_conditions = {
            c_e: {"left": (lbc, "Neumann"), "right": (pybamm.Scalar(0), "Neumann")},
        }
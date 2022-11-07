#
# Class for electrolyte diffusion employing stefan-maxwell and consider diffusion of both EC and lithium ion 
#
import pybamm
import numpy as np
from .base_solvent_diffusion import BaseSolventDiffusion


class Single_SpatialConsume_w_refill(BaseSolventDiffusion):
    """Class for conservation of solvent in the electrolyte employing the
    Stefan-Maxwell constitutive equations. and 
    consider diffusion of both solvent and lithium ion 

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
        eps_c_EC_dict = {}
        for domain in self.options.whole_cell_domains:
            Domain = domain.capitalize()
            eps_c_EC_k = pybamm.Variable(
                f"{Domain} porosity times EC concentration",
                domain=domain,
                auxiliary_domains={"secondary": "current collector"},
                bounds=(0, np.inf),
            )
            eps_c_EC_k.print_name = f"eps_c_EC_{domain[0]}"
            eps_c_EC_dict[domain] = eps_c_EC_k

        variables = self._get_standard_porosity_times_EC_concentration_variables(
            eps_c_EC_dict
        )

        return variables

    def get_coupled_variables(self, variables):

        c_EC_dict = {}
        for domain in self.options.whole_cell_domains:
            Domain = domain.capitalize()
            eps_k = variables[f"{Domain} porosity"]
            eps_c_EC_k = variables[f"{Domain} porosity times EC concentration"]
            c_EC_k = eps_c_EC_k / eps_k
            c_EC_dict[domain] = c_EC_k
        variables.update(self._get_standard_EC_concentration_variables(c_EC_dict))


        sign_2_s = pybamm.FullBroadcast(
                pybamm.Scalar(0), "separator", 
                auxiliary_domains={"secondary": "current collector"})
        sign_2_p = pybamm.FullBroadcast(
                pybamm.Scalar(0), "positive electrode", 
                auxiliary_domains={"secondary": "current collector"})
        

        eps_c_EC = variables["Porosity times EC concentration"]
        Q_sei = variables["Loss of lithium to SEI [mol]"]

        param = self.param

        N_EC = pybamm.FullBroadcastToEdges(
            0,
            [domain for domain in self.options.whole_cell_domains],
            "current collector",)
        

        a_p = variables["Positive electrode surface area to volume ratio"]
        a_n = variables["Negative electrode surface area to volume ratio"]
        zero_s = pybamm.FullBroadcast(0, "separator", "current collector")
        a = pybamm.concatenation(a_n, zero_s, a_p)

        j_inner =  variables["Inner SEI interfacial current density"]
        j_outer =  variables["Outer SEI interfacial current density"]
        j_SEI = j_inner + j_outer
        j_sign_SEI = pybamm.concatenation(j_SEI, sign_2_s, sign_2_p )

        ## EC:lithium:SEI=2:2:1    ; change to EC:lithium:SEI=2:1:1 for now
        ratio_sei_li = - 1/param.z_sei  ; # change to 1 for now , initially is 0.5
        ratio_ec_li  = 1 ; 

        source_terms_ec =(
            a * j_sign_SEI / param.gamma_e 
            / param.gamma_e_ec_Rio * ratio_ec_li)
        # Mark Ruihe update 221021 - ignore Li+ (de-)intercalation to 
        # avoid differences between w and wo refill when there are no SEI 
        source_terms_refill = - (
            a * j_sign_SEI / param.gamma_e / param.gamma_e_ec_Rio 
            * param.c_ec_0_dim * (
                ratio_ec_li * param.Vmolar_ec 
                +ratio_sei_li * param.Vmolar_CH2OCO2Li2)
            + a * j_sign_SEI / param.gamma_e  / param.gamma_e_ec_Rio 
            * param.c_ec_0_dim * param.Vmolar_Li ) 
        """ source_terms_refill = (
            a * j_sign_SEI / param.gamma_e / param.gamma_e_ec_Rio 
            * 4541 * (
                ratio_ec_li * param.Vmolar_ec 
                +ratio_sei_li * param.Vmolar_CH2OCO2Li2)
            + a * j_sign_SEI / param.gamma_e  / param.gamma_e_ec_Rio 
            * 4541 * param.Vmolar_Li ) """
        
        variables.update(self._get_standard_EC_flux_variables(
            N_EC,N_EC,N_EC,N_EC,
            source_terms_ec,source_terms_refill))
        variables.update(
            self._get_total_EC_concentration_electrolyte(eps_c_EC,Q_sei))

        return variables

    def set_rhs(self, variables):   # Mark Ruihe Li modify  ? can we have more than one variables in this file?

        param = self.param

        eps_c_EC = variables["Porosity times EC concentration"]

        N_EC = variables["EC flux"]
        div_Vbox = variables["Transverse volume-averaged acceleration"]
        source_terms_ec =  variables["EC source term (SEI)"]
        source_terms_refill = variables["EC source term refill"]

        #print('using double spatial consume w refill for EC')
        self.rhs = {
            eps_c_EC: - pybamm.div(N_EC) * param.tau_discharge / param.tau_ec_Rio 
            # source term due to SEI
            + source_terms_ec
            # source term due to replenishment 
            + source_terms_refill
        } 

    def set_initial_conditions(self, variables):

        eps_c_EC= variables["Porosity times EC concentration"]

        self.initial_conditions = {
            eps_c_EC: self.param.epsilon_init * self.param.c_ec_init,   # This needs to be changed to be the same as Andrew's results
        }

    def set_boundary_conditions(self, variables):
        param = self.param
        c_EC = variables["EC concentration"]

        self.boundary_conditions = {
            c_EC: {"left": (pybamm.Scalar(0), "Neumann"), "right": (pybamm.Scalar(0), "Neumann")}, # This needs to be changed to be the same as Andrew's results
        }

#
# Class for electrolyte diffusion employing stefan-maxwell and consider diffusion of both EC and lithium ion 
#
import pybamm
import numpy as np
from .base_solvent_diffusion import BaseSolventDiffusion


class OneSolventDiffusion_w_Refill(BaseSolventDiffusion):
    """Class for conservation of solvent in the electrolyte employing the
    Stefan-Maxwell constitutive equations. and consider diffusion of both solvent and lithium ion 

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
            eps_c_EC_n    = None 
        else:
            eps_c_EC_n    = pybamm.standard_variables.eps_c_EC_n     
        eps_c_EC_s    = pybamm.standard_variables.eps_c_EC_s         
        eps_c_EC_p    = pybamm.standard_variables.eps_c_EC_p         

        variables = self._get_standard_porosity_times_EC_concentration_variables( 
            eps_c_EC_n, eps_c_EC_s, eps_c_EC_p
        )
        return variables

    def get_coupled_variables(self, variables):

        param = self.param

        if self.half_cell:   # Mark Ruihe Li comment
            c_EC_n    = None 
        else:
            eps_n = variables["Negative electrode porosity"]
            eps_c_EC_n= variables["Negative electrode porosity times EC concentration"]
            c_EC_n = eps_c_EC_n / eps_n 

        eps_n = variables["Negative electrode porosity"]
        eps_c_EC_n= variables["Negative electrode porosity times EC concentration"]
        c_EC_n = eps_c_EC_n / eps_n

        eps_s = variables["Separator porosity"]
        eps_p = variables["Positive electrode porosity"]
        eps_c_EC_s = variables["Separator porosity times EC concentration"]
        eps_c_EC_p = variables["Positive electrode porosity times EC concentration"]
        c_EC_s = eps_c_EC_s / eps_s
        c_EC_p = eps_c_EC_p / eps_p
        variables.update(
            self._get_standard_EC_concentration_variables(c_EC_n, c_EC_s, c_EC_p) 
        )

        tor = variables["Electrolyte transport efficiency"]
        i_e = variables["Electrolyte current density"]
        T = variables["Cell temperature"]

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


        
        c_EC  = variables["EC concentration"]
        c_e = variables["Electrolyte concentration"]

        eps_c_EC = variables["Porosity times EC concentration"]

        Q_sei = variables["Loss of lithium to SEI [mol]"]

        N_EC_diffusion = -tor * param.D_ec * pybamm.grad(c_EC)
        N_cross_diffusion = -(
            param.tau_ec_Rio  * param.EC_ratio_Rio / 
            param.tau_cross_Rio /param.gamma_e_ec_Rio* 
            tor * param.D_ec_Li_cross * pybamm.grad(c_e) )
        # Mark Ruihe: most important correction that change again everything... update 221021
        #there should be a minus here!
        N_EC_migration = - param.C_e* param.Xi * i_e / param.gamma_e *(
            param.tau_ec_Rio / param.gamma_e_ec_Rio / param.tau_diffusion_e
        )

        N_EC = N_EC_diffusion + N_cross_diffusion + N_EC_migration

        a_p = variables["Positive electrode surface area to volume ratio"]
        a_n = variables["Negative electrode surface area to volume ratio"]
        zero_s = pybamm.FullBroadcast(0, "separator", "current collector")
        a = pybamm.concatenation(a_n, zero_s, a_p)

        j_inner =  variables["Inner SEI interfacial current density"]
        j_outer =  variables["Outer SEI interfacial current density"]
        j_SEI = j_inner + j_outer
        j_sign_SEI = pybamm.concatenation(j_SEI, sign_2_s, sign_2_p )

        sum_s_j = variables["Sum of electrolyte reaction source terms"]
        sum_s_j.print_name = "a"
        source_terms = sum_s_j / self.param.gamma_e
        ## EC:lithium:SEI=2:2:1    ; change to EC:lithium:SEI=2:1:1 for now

        ratio_sei_li = - 1/param.z_sei  ; # change to 1 for now , initially is 0.5
        ratio_ec_li  = 1 ; 

        source_terms_ec =(
            a * j_sign_SEI / param.gamma_e 
            / param.gamma_e_ec_Rio * ratio_ec_li)
        source_terms_refill = (
            a * j_sign_SEI / param.gamma_e / param.gamma_e_ec_Rio 
            * param.c_ec_0_dim * (
                ratio_ec_li * param.Vmolar_ec 
                +ratio_sei_li * param.Vmolar_CH2OCO2Li2)
            + source_terms / param.gamma_e_ec_Rio 
            * param.c_ec_0_dim * param.Vmolar_Li )
        
        variables.update(self._get_standard_EC_flux_variables(
            N_EC,N_EC_diffusion,N_EC_migration,N_cross_diffusion,
            source_terms_ec,source_terms_refill))
        variables.update(
            self._get_total_EC_concentration_electrolyte(eps_c_EC,Q_sei))

        return variables

    def set_rhs(self, variables):   # Mark Ruihe Li modify  ? can we have more than one variables in this file?

        param = self.param
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

        a_p = variables["Positive electrode surface area to volume ratio"]
        a_n = variables["Negative electrode surface area to volume ratio"]
        zero_s = pybamm.FullBroadcast(0, "separator", "current collector")
        a = pybamm.concatenation(a_n, zero_s, a_p)
        c_e = variables["Electrolyte concentration"]
        c_EC  = variables["EC concentration"]
        eps_c_EC = variables["Porosity times EC concentration"]
        tor = variables["Electrolyte transport efficiency"]
        T = variables["Cell temperature"]
        N_EC = variables["EC flux"]
        div_Vbox = variables["Transverse volume-averaged acceleration"]
        source_terms_ec =  variables["EC source term (SEI)"]
        source_terms_refill = variables["EC source term refill"]

        #print('using EC w refill for EC')
        self.rhs = {
            eps_c_EC: - pybamm.div(N_EC) * param.tau_discharge / param.tau_ec_Rio 
            # source term due to SEI
            + source_terms_ec
            # source term due to replenishment 
            + source_terms_refill

            #(   old ones:
            #    param.EC_ratio_Rio / param.gamma_e_ec_Rio * param.tau_discharge  / 
            #    param.tau_cross_Rio * pybamm.div(tor * param.D_ec_Li_cross * pybamm.grad(c_e))
            #    )
            #+ param.tau_discharge  / param.tau_ec_Rio * pybamm.div(tor * param.D_ec * pybamm.grad(c_EC))
            ## source term due to migration:
            #+    source_terms / param.gamma_e_ec_Rio * param.Xi   
            
            ## source term due to SEI
            #+ a * j_sign_SEI /  param.gamma_e / param.gamma_e_ec_Rio 
            # source term due to replenishment
            #- source_terms / param.gamma_e_ec_Rio * param.Vmolar_Li * param.c_ec_0_dim   # ignore volume of lithium
            #+ (  
            #a * j_sign_SEI /  param.gamma_e / param.gamma_e_ec_Rio *   
            #(  ratio_sei_li*param.Vmolar_CH2OCO2Li2*param.c_ec_0_dim 
            #- param.Vmolar_ec*param.c_ec_0_dim ) ) 
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

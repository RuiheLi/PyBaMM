#
# Submodel for no convection in transverse directions
#
import pybamm
import numpy as np
from .base_solvent_diffusion import BaseSolventDiffusion


class NoSolventDiffusion(BaseSolventDiffusion):
    """
    Submodel for no convection in transverse directions

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    options : dict, optional
        A dictionary of options to be passed to the model.

    **Extends:** :class:`pybamm.convection.through_cell.BaseTransverseModel`
    """

    def __init__(self, param, options=None):
        super().__init__(param, options=options)

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
            self._get_standard_EC_concentration_variables(c_EC_n, c_EC_s, c_EC_p)  # Already Written
        )

        eps_c_EC = variables["Porosity times EC concentration"]

        variables.update(
            self._get_total_EC_concentration_electrolyte(eps_c_EC))
            
        return variables

    def set_rhs(self, variables): 

        eps_c_EC = variables["Porosity times EC concentration"]
        if self.half_cell:
            deps_c_EC_n_dt = None
        else:
            deps_c_EC_n_dt = pybamm.FullBroadcast(
                0, "negative electrode", "current collector")
        deps_c_EC_s_dt = pybamm.FullBroadcast(0, "separator", "current collector")
        deps_c_EC_p_dt = pybamm.FullBroadcast(0, "positive electrode", "current collector")

        deps_c_EC_dt = pybamm.concatenation(
            deps_c_EC_n_dt, deps_c_EC_s_dt, deps_c_EC_p_dt)

        self.rhs = {eps_c_EC: deps_c_EC_dt} 

        
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

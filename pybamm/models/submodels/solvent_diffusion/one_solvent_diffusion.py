#
# Class for electrolyte diffusion employing stefan-maxwell and consider diffusion of both EC and lithium ion 
#
import pybamm
import numpy as np
from .base_solvent_diffusion import BaseSolventDiffusion


class OneSolventDiffusion(BaseSolventDiffusion):
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
            eps_c_EC_n    = None # Mark Ruihe Li add
        else:
            eps_c_EC_n    = pybamm.standard_variables.eps_c_EC_n  # Mark Ruihe Li add   already add to standard_variables     can it not be standard_variables?
        eps_c_EC_s    = pybamm.standard_variables.eps_c_EC_s  # Mark Ruihe Li add       already add to standard_variables     can it not be standard_variables?
        eps_c_EC_p    = pybamm.standard_variables.eps_c_EC_p  # Mark Ruihe Li add       already add to standard_variables     can it not be standard_variables?

        variables = self._get_standard_porosity_times_EC_concentration_variables( 
            eps_c_EC_n, eps_c_EC_s, eps_c_EC_p
        )


        return variables

    def get_coupled_variables(self, variables):

        if self.half_cell:   # Mark Ruihe Li comment
            c_EC_n    = None # Mark Ruihe Li add
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
        c_EC = variables["EC concentration"]

        #N_EC = c_EC* 122333    ###     need to write an expression for EC flux, but give up for now

        #variables.update(self._get_standard_flux_variables(N_EC))
        variables.update(
            self._get_total_EC_concentration_electrolyte(eps_c_EC))

        return variables

    def set_rhs(self, variables):   # Mark Ruihe Li modify  ? can we have more than one variables in this file?

        param = self.param
        eps_c_EC = variables["Porosity times EC concentration"]

        
        self.rhs = {
            eps_c_EC: 1
        } 


    def set_initial_conditions(self, variables):

        eps_c_EC= variables["Porosity times EC concentration"]

        self.initial_conditions = {
            eps_c_EC: 0.1,   # This needs to be changed to be the same as Andrew's results
        }

    def set_boundary_conditions(self, variables):
        param = self.param
        c_EC = variables["EC concentration"]

        self.boundary_conditions = {
            c_EC: {"left": (pybamm.Scalar(0), "Neumann"), "right": (pybamm.Scalar(0), "Neumann")}, # This needs to be changed to be the same as Andrew's results
        }

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
        c_EC_typ_n = pybamm.FullBroadcast(
            1.0, "negative electrode",
            "current collector") # 
        c_EC_typ_s = pybamm.FullBroadcast(
            1.0, "separator",
            "current collector")
        c_EC_typ_p = pybamm.FullBroadcast(
            1.0, "positive electrode",
            "current collector")
        variables=self._get_standard_EC_concentration_variables(
            c_EC_typ_n,c_EC_typ_s,c_EC_typ_p)
        print("have you come here?")

        param = self.param

        param.s.epsilon_init

        variables.update(
            self._get_standard_porosity_times_EC_concentration_variables(
            param.n.epsilon_init,
            param.s.epsilon_init,
            param.p.epsilon_init)
        )

        return variables


    def get_coupled_variables(self, variables):
        
        eps_c_EC = variables["Porosity times EC concentration"]
        #Q_sei = variables["Loss of lithium to SEI [mol]"]

        N_EC_n = pybamm.FullBroadcast(
            0.0, "negative electrode",
            "current collector") 
        N_EC_s = pybamm.FullBroadcast(
            0.0, "separator",
            "current collector")  
        N_EC_p = pybamm.FullBroadcast(
            0.0, "positive electrode",
            "current collector") 
        N_EC = pybamm.concatenation(N_EC_n, N_EC_s, N_EC_p )

        variables.update(self._get_standard_EC_flux_variables(N_EC))

        variables.update(
            self._get_total_EC_concentration_electrolyte(eps_c_EC,0))

        return variables

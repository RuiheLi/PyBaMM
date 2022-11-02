#
# Submodel for no convection in transverse directions
#
import pybamm
import numpy as np
from .base_solvent_diffusion import BaseSolventDiffusion


class NoSolventDiffusion(BaseSolventDiffusion):
    """Class for constant EC concentration of electrolyte

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
        c_EC_dict = {
            domain: pybamm.FullBroadcast(1, domain, "current collector")
            for domain in self.options.whole_cell_domains
        }
        variables=self._get_standard_EC_concentration_variables(c_EC_dict)
        #print("have you come here?")
        N_EC = pybamm.FullBroadcastToEdges(
            0,
            [domain for domain in self.options.whole_cell_domains],
            "current collector",)
        
        variables.update(self._get_standard_EC_flux_variables(
            N_EC,N_EC,N_EC,N_EC,N_EC,N_EC))

        return variables

    def get_coupled_variables(self, variables):
        eps_c_EC_dict = {}
        for domain in self.options.whole_cell_domains:
            Domain = domain.capitalize()
            eps_k = variables[f"{Domain} porosity"]
            c_EC_k = variables[f"{Domain.split()[0]} EC concentration"]
            eps_c_EC_dict[domain] = eps_k * c_EC_k
        variables.update(
            self._get_standard_porosity_times_EC_concentration_variables(eps_c_EC_dict)
        )

        return variables

    def set_boundary_conditions(self, variables):
        """
        We provide boundary conditions even though the concentration is constant
        so that the gradient of the concentration has the correct shape after
        discretisation.
        """

        c_EC = variables["EC concentration"]

        self.boundary_conditions = {
            c_EC: {
                "left": (pybamm.Scalar(0), "Neumann"),
                "right": (pybamm.Scalar(0), "Neumann"),
            }
        }

    def set_events(self, variables):
        # No event since the concentration is constant
        pass

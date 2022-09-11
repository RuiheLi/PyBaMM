#
# Submodel for no convection in transverse directions
#
import pybamm
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

        param = self.param
        c_ec_typ = 4541
        c_ec_typ_n = pybamm.PrimaryBroadcast(
            pybamm.Scalar(c_ec_typ), "negative electrode",
            "current collector")
        c_ec_typ_s = pybamm.PrimaryBroadcast(
            pybamm.Scalar(c_ec_typ), "separator",
            "current collector")
        c_ec_typ_p = pybamm.PrimaryBroadcast(
            pybamm.Scalar(c_ec_typ), "positive electrode",
            "current collector")

        variables=self._get_standard_EC_concentration_variables(
            c_ec_typ_n,c_ec_typ_s,c_ec_typ_p)

        variables.update(
            self._get_standard_porosity_times_EC_concentration_variables(
            c_ec_typ_n * 0.25,
            c_ec_typ_s * 0.47,
            c_ec_typ_p * 0.335)
        )

        eps_c_EC = variables["Porosity times EC concentration"]
        variables.update(
            self._get_total_EC_concentration_electrolyte(eps_c_EC))

        print("have you come here?")
       
        return variables

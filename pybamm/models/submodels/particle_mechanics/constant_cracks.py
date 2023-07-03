#
# Class for no mechanics
#
import pybamm
from .base_mechanics import BaseMechanics


class ConstantCracks(BaseMechanics):
    """
    Class for swelling with constant crack length.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain of the model either 'Negative' or 'Positive'
    options: dict
        A dictionary of options to be passed to the model.
        See :class:`pybamm.BaseBatteryModel`
    phase : str, optional
        Phase of the particle (default is "primary")

    **Extends:** :class:`pybamm.particle_mechanics.BaseMechanics`
    """

    def __init__(self, param, domain, options, phase="primary"):
        super().__init__(param, domain, options, phase)

    def get_fundamental_variables(self):
        domain, Domain = self.domain_Domain

        l_cr_0 = self.param.l_cr_0
        variables = self._get_standard_variables(l_cr_0)
        zero = pybamm.FullBroadcast(
            pybamm.Scalar(0), f"{domain} electrode", "current collector"
        )
        zero_av = pybamm.x_average(zero)
        variables.update(
            {
                f"{Domain} particle cracking rate [m.s-1]": zero,
                f"X-averaged {domain} particle cracking rate [m.s-1]": zero_av,
            }
        )
        return variables

    def get_coupled_variables(self, variables):
        variables.update(self._get_standard_surface_variables(variables))
        variables.update(self._get_mechanical_results(variables))
        return variables

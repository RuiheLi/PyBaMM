#
# Base class for solvent diffusion (consider EC first)
#
import pybamm


class BaseSolventDiffusion(pybamm.BaseSubModel):
    """Base class for conservation of solvent in the electrolyte.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    options : dict, optional
        A dictionary of options to be passed to the model.

    **Extends:** :class:`pybamm.BaseSubModel`
    """

    def __init__(self, param, options=None):
        super().__init__(param, options=options)


    def _get_standard_EC_concentration_variables(self, c_EC_n, c_EC_s, c_EC_p):
        """
        A private function to obtain the standard variables which
        can be derived from the concentration in the EC.

        Parameters
        ----------
        c_EC_n : :class:`pybamm.Symbol`
            The EC concentration in the negative electrode.
        c_EC_s : :class:`pybamm.Symbol`
            The EC concentration in the separator.
        c_EC_p : :class:`pybamm.Symbol`
            The EC concentration in the positive electrode.

        Returns
        -------
        variables : dict
            The variables which can be derived from the concentration in the
            EC.
        """

        c_ec_typ = self.param.c_ec_typ
        c_EC = pybamm.concatenation(c_EC_n, c_EC_s, c_EC_p)

        if self.half_cell:
            # overwrite c_EC_n to be the boundary value of c_EC_s
            c_EC_n = pybamm.boundary_value(c_EC_s, "left")

        c_EC_n_av = pybamm.x_average(c_EC_n)
        c_EC_av = pybamm.x_average(c_EC)
        c_EC_s_av = pybamm.x_average(c_EC_s)
        c_EC_p_av = pybamm.x_average(c_EC_p)

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
 
        variables = {
            "EC concentration": c_EC,
            "EC concentration [mol.m-3]": c_ec_typ * c_EC,
            "EC concentration [Molar]": c_ec_typ * c_EC / 1000,
            "X-averaged EC concentration": c_EC_av,
            "X-averaged EC concentration [mol.m-3]": c_ec_typ * c_EC_av,
            "X-averaged EC concentration [Molar]": c_ec_typ * c_EC_av / 1000,
            "Negative EC concentration": c_EC_n,
            "Negative EC concentration [mol.m-3]": c_ec_typ * c_EC_n,
            "Negative EC concentration [Molar]": c_ec_typ * c_EC_n / 1000,
            "Separator EC concentration": c_EC_s,
            "Separator EC concentration [mol.m-3]": c_ec_typ * c_EC_s,
            "Separator EC concentration [Molar]": c_ec_typ * c_EC_s / 1000,
            "Positive EC concentration": c_EC_p,
            "Positive EC concentration [mol.m-3]": c_ec_typ * c_EC_p,
            "Positive EC concentration [Molar]": c_ec_typ * c_EC_p / 1000,
            "X-averaged negative EC concentration": c_EC_n_av,
            "X-averaged negative EC concentration [mol.m-3]": c_ec_typ
            * c_EC_n_av,
            "X-averaged separator EC concentration": c_EC_s_av,
            "X-averaged separator EC concentration [mol.m-3]": c_ec_typ
            * c_EC_s_av,
            "X-averaged positive EC concentration": c_EC_p_av,
            "X-averaged positive EC concentration [mol.m-3]": c_ec_typ
            * c_EC_p_av,
            
        }

        # Override print_name
        c_EC.print_name = "c_EC"

        return variables
    
    def _get_standard_EC_flux_variables(self, 
        N_EC,N_EC_diffusion,N_EC_migration,N_cross_diffusion):
        """
        A private function to obtain the standard variables which
        can be derived from the mass flux of EC in the electrolyte.

        Parameters
        ----------
        N_EC : :class:`pybamm.Symbol`
            The EC flux in the electrolyte.

        Returns
        -------
        variables : dict
            The variables which can be derived from the flux in the
            electrolyte.
        """

        param = self.param
        flux_EC_scale = param.D_ec_typ * param.c_ec_typ / param.L_x

        variables = {
            "EC flux": N_EC,
            "EC flux [mol.m-2.s-1]": N_EC * flux_EC_scale,
            "EC flux by diffusion": N_EC_diffusion,
            "EC flux by diffusion [mol.m-2.s-1]": N_EC_diffusion * flux_EC_scale,
            "EC flux by migration": N_EC_migration,
            "EC flux by migration [mol.m-2.s-1]": N_EC_migration * flux_EC_scale,
            "EC flux by Li+": N_cross_diffusion,
            "EC flux by Li+ [mol.m-2.s-1]": N_cross_diffusion * flux_EC_scale,
        }

        return variables


    
    def _get_standard_porosity_times_EC_concentration_variables(
        self, eps_c_EC_n, eps_c_EC_s, eps_c_EC_p
    ):
        eps_c_EC = pybamm.concatenation(eps_c_EC_n, eps_c_EC_s, eps_c_EC_p)

        variables = {
            "Porosity times EC concentration": eps_c_EC,
            "Separator porosity times EC concentration": eps_c_EC_s,
            "Positive electrode porosity times EC concentration": eps_c_EC_p,
        }

        if not self.half_cell:
            variables.update(
                {"Negative electrode porosity times EC concentration": eps_c_EC_n}
            )
        return variables

    def _get_total_EC_concentration_electrolyte(self, eps_c_EC,Q_sei):
        """
        A private function to obtain the total EC concentration in the electrolyte.

        Parameters
        ----------
        eps_c_EC : :class:`pybamm.Symbol`
            Porosity times EC concentration

        Returns
        -------
        variables : dict
            The "Total EC in electrolyte [mol]" variable.
        """
        
        c_ec_typ = self.param.c_ec_typ
        L_x = self.param.L_x
        A = self.param.A_cc

        eps_c_EC_av = pybamm.yz_average(pybamm.x_average(eps_c_EC))

        variables = {
            "Total EC in electrolyte": eps_c_EC_av,
            "Total EC in electrolyte [mol]": c_ec_typ * L_x * A * eps_c_EC_av,
            "Total EC in electrolyte and SEI [mol]": c_ec_typ * L_x * A * eps_c_EC_av + Q_sei,
        }

        return variables
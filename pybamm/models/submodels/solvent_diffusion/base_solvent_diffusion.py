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

    def _get_standard_EC_concentration_variables(self, c_EC_dict):
        """
        A private function to obtain the standard variables which
        can be derived from the concentration in the EC.

        Parameters
        ----------
        c_EC_dict : dict of :class:`pybamm.Symbol`
            EC concentrations in the various domains

        Returns
        -------
        variables : dict
            The variables which can be derived from the concentration in the
            EC.
        """

        c_ec_typ = self.param.c_ec_typ
        c_EC = pybamm.concatenation(*c_EC_dict.values())
        # Override print_name
        c_EC.print_name = "c_EC"

        variables = {
            "EC concentration": c_EC,
            "X-averaged EC concentration": pybamm.x_average(c_EC),
        }

        # Case where an electrode is not included (half-cell)
        if "negative electrode" not in self.options.whole_cell_domains:
            c_EC_s = c_EC_dict["separator"]
            c_EC_dict["negative electrode"] = pybamm.boundary_value(c_EC_s, "left")
 
        for domain, c_EC_k in c_EC_dict.items():
            domain = domain.split()[0]
            Domain = domain.capitalize()
            c_EC_k_av = pybamm.x_average(c_EC_k)
            variables.update(
                {
                    f"{Domain} EC concentration": c_EC_k,
                    f"X-averaged {domain} EC concentration": c_EC_k_av,
                }
            )

        # Calculate dimensional variables
        variables_nondim = variables.copy()
        for name, var in variables_nondim.items():
            variables.update(
                {
                    f"{name} [mol.m-3]": c_ec_typ * var,
                    f"{name} [Molar]": c_ec_typ * var / 1000,
                }
            )

        return variables
    
    def _get_standard_porosity_times_EC_concentration_variables(self, eps_c_EC_dict):
        eps_c_EC = pybamm.concatenation(*eps_c_EC_dict.values())
        variables = {"Porosity times EC concentration": eps_c_EC}

        for domain, eps_c_EC_k in eps_c_EC_dict.items():
            Domain = domain.capitalize()
            variables[f"{Domain} porosity times EC concentration"] = eps_c_EC_k

        # Total lithium concentration in electrolyte
        # Mark: Ruihe change, need to add Q_sei
        # variables.update(self._get_total_EC_concentration_electrolyte(eps_c_e))

        return variables


    def _get_standard_EC_flux_variables(self, 
        N_EC,N_EC_diffusion,N_EC_migration,N_cross_diffusion,
        source_terms_ec,source_terms_refill):
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
        source_EC_scale = flux_EC_scale / param.L_x

        variables = {
            "EC flux": N_EC,
            "Minus div EC flux": - pybamm.div(N_EC) * param.tau_discharge / param.tau_ec_Rio ,
            "Minus div EC flux [mol.m-3.s-1]": (
            - pybamm.div(N_EC) * param.tau_discharge 
            / param.tau_ec_Rio *source_EC_scale),
            "EC source term (SEI)": source_terms_ec,
            "EC source term refill": source_terms_refill,
            "EC source term (SEI) [mol.m-3.s-1]": source_terms_ec*source_EC_scale,
            "EC source term refill [mol.m-3.s-1]": source_terms_refill*source_EC_scale,
            "EC flux [mol.m-2.s-1]": N_EC * flux_EC_scale,
            "EC flux by diffusion": N_EC_diffusion,
            "Minus div EC flux by diffusion": (
            - pybamm.div(N_EC_diffusion) * param.tau_discharge 
            / param.tau_ec_Rio)  ,
            "Minus div EC flux by diffusion [mol.m-3.s-1]": (
            - pybamm.div(N_EC_diffusion) * param.tau_discharge 
            / param.tau_ec_Rio) * source_EC_scale ,
            "EC flux by diffusion [mol.m-2.s-1]": N_EC_diffusion * flux_EC_scale,
            "EC flux by migration": N_EC_migration,
            "Minus div EC flux by migration": (
            - pybamm.div(N_EC_migration) * param.tau_discharge 
            / param.tau_ec_Rio)  ,
            "Minus div EC flux by migration [mol.m-3.s-1]": (
            - pybamm.div(N_EC_migration) * param.tau_discharge 
            / param.tau_ec_Rio)  * source_EC_scale ,
            "EC flux by migration [mol.m-2.s-1]": N_EC_migration * flux_EC_scale,
            "EC flux by Li+": N_cross_diffusion,
            "Minus div EC flux by Li+": (
            - pybamm.div(N_cross_diffusion) * param.tau_discharge 
            / param.tau_ec_Rio ) ,
            "Minus div EC flux by Li+ [mol.m-3.s-1]": (
            - pybamm.div(N_cross_diffusion) * param.tau_discharge 
            / param.tau_ec_Rio ) * source_EC_scale,
            "EC flux by Li+ [mol.m-2.s-1]": N_cross_diffusion * flux_EC_scale,
        }

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
            "Total EC in electrolyte and SEI [mol]": (
                c_ec_typ * L_x * A * eps_c_EC_av + Q_sei
            ),
        }

        return variables
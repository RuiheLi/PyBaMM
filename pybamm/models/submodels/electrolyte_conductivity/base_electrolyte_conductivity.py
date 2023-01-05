#
# Base class for electrolyte conductivity
#

import pybamm


class BaseElectrolyteConductivity(pybamm.BaseSubModel):
    """Base class for conservation of charge in the electrolyte.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str, optional
        The domain in which the model holds
    options : dict, optional
        A dictionary of options to be passed to the model.

    **Extends:** :class:`pybamm.BaseSubModel`
    """

    def __init__(self, param, domain=None, options=None):
        super().__init__(param, domain, options=options)

    def _get_standard_potential_variables(self, phi_e_dict):
        """
        A private function to obtain the standard variables which
        can be derived from the potential in the electrolyte.

        Parameters
        ----------
        phi_e_dict : dict of :class:`pybamm.Symbol`
            Dictionary of electrolyte potentials in the relevant domains

        Returns
        -------
        variables : dict
            The variables which can be derived from the potential in the
            electrolyte.
        """

        param = self.param
        pot_scale = param.potential_scale
        U_ref = param.n.U_ref

        phi_e = pybamm.concatenation(*phi_e_dict.values())

        # Case where negative electrode is not included (half-cell)
        if "negative electrode" not in self.options.whole_cell_domains:
            phi_e_s = phi_e_dict["separator"]
            phi_e_dict["negative electrode"] = pybamm.boundary_value(phi_e_s, "left")

        eta_e_av = pybamm.x_average(
            phi_e_dict["positive electrode"]
        ) - pybamm.x_average(phi_e_dict["negative electrode"])
        phi_e_av = pybamm.x_average(phi_e)

        variables = {
            "Electrolyte potential": phi_e,
            "Electrolyte potential [V]": -U_ref + pot_scale * phi_e,
            "X-averaged electrolyte potential": phi_e_av,
            "X-averaged electrolyte potential [V]": -U_ref + pot_scale * phi_e_av,
            "X-averaged electrolyte overpotential": eta_e_av,
            "X-averaged electrolyte overpotential [V]": pot_scale * eta_e_av,
            "Gradient of electrolyte potential": pybamm.grad(phi_e),
        }

        for domain, phi_e_k in phi_e_dict.items():
            name = f"{domain.split()[0]} electrolyte potential"
            Name = name.capitalize()
            phi_e_k_av = pybamm.x_average(phi_e_k)
            variables.update(
                {
                    f"{Name}": phi_e_k,
                    f"{Name} [V]": -U_ref + pot_scale * phi_e_k,
                    f"X-averaged {name}": phi_e_k_av,
                    f"X-averaged {name} [V]": -U_ref + pot_scale * phi_e_k_av,
                }
            )
            if domain in self.options.whole_cell_domains:
                variables[f"Gradient of {name}"] = pybamm.grad(phi_e_k)

        return variables

    def _get_standard_current_variables(self, i_e):
        """
        A private function to obtain the standard variables which
        can be derived from the current in the electrolyte.

        Parameters
        ----------
        i_e : :class:`pybamm.Symbol`
            The current in the electrolyte.

        Returns
        -------
        variables : dict
            The variables which can be derived from the current in the
            electrolyte.
        """

        i_typ = self.param.i_typ
        variables = {
            "Electrolyte current density": i_e,
            "Electrolyte current density [A.m-2]": i_typ * i_e,
        }

        if isinstance(i_e, pybamm.Concatenation):
            if self.options.whole_cell_domains == [
                "negative electrode",
                "separator",
                "positive electrode",
            ]:
                i_e_n, _, i_e_p = i_e.orphans
            elif self.options.whole_cell_domains == ["separator", "positive electrode"]:
                _, i_e_p = i_e.orphans
                i_e_n = None

            if i_e_n is not None:
                variables.update(
                    {
                        "Negative electrolyte current density": i_e_n,
                        "Negative electrolyte current density [A.m-2]": i_e_n * i_typ,
                    }
                )
            if i_e_p is not None:
                variables.update(
                    {
                        "Positive electrolyte current density": i_e_p,
                        "Positive electrolyte current density [A.m-2]": i_e_p * i_typ,
                    }
                )

        return variables

    def _get_split_overpotential(self, eta_c_av, eta_cEC_av, delta_phi_e_av):
        """
        A private function to obtain the standard variables which
        can be derived from the electrode-averaged concentration
        overpotential and Ohmic losses in the electrolyte.

        Parameters
        ----------
        eta_c_av : :class:`pybamm.Symbol`
            The electrode-averaged concentration overpotential - due to c(Li+) only!
        eta_cEC_av : :class:`pybamm.Symbol`
            The electrode-averaged concentration overpotential - due to c(EC) only! - Add by Ruihe Li
        delta_phi_e_av: :class:`pybamm.Symbol`
            The electrode-averaged electrolyte Ohmic losses

        Returns
        -------
        variables : dict
            The variables which can be derived from the electrode-averaged
            concentration overpotential and Ohmic losses in the electrolyte
            electrolyte.
        """

        param = self.param
        pot_scale = param.potential_scale

        variables = {
            "X-averaged concentration overpotential": eta_c_av,
            "X-averaged EC concentration overpotential": eta_cEC_av,
            "X-averaged electrolyte ohmic losses": delta_phi_e_av,
            "X-averaged concentration overpotential [V]": pot_scale * eta_c_av,
            "X-averaged EC concentration overpotential [V]": pot_scale * eta_cEC_av,
            "X-averaged electrolyte ohmic losses [V]": pot_scale * delta_phi_e_av,
        }

        return variables

    def _get_standard_average_surface_potential_difference_variables(
        self, delta_phi_av
    ):
        """
        A private function to obtain the standard variables which
        can be derived from the surface potential difference.

        Parameters
        ----------
        delta_phi_av : :class:`pybamm.Symbol`
            The x-averaged surface potential difference.

        Returns
        -------
        variables : dict
            The variables which can be derived from the surface potential difference.
        """
        domain = self.domain
        ocp_ref = self.domain_param.U_ref

        variables = {
            f"X-averaged {domain} electrode surface potential difference": delta_phi_av,
            f"X-averaged {domain} electrode surface potential difference [V]": ocp_ref
            + delta_phi_av * self.param.potential_scale,
        }

        return variables

    def _get_standard_surface_potential_difference_variables(self, delta_phi):
        """
        A private function to obtain the standard variables which
        can be derived from the surface potential difference.

        Parameters
        ----------
        delta_phi : :class:`pybamm.Symbol`
            The surface potential difference.

        Returns
        -------
        variables : dict
            The variables which can be derived from the surface potential difference.
        """
        domain, Domain = self.domain_Domain

        ocp_ref = self.domain_param.U_ref

        # Broadcast if necessary
        if delta_phi.domain == ["current collector"]:
            delta_phi = pybamm.PrimaryBroadcast(delta_phi, f"{domain} electrode")

        variables = {
            f"{Domain} electrode surface potential difference": delta_phi,
            f"{Domain} electrode surface potential difference [V]": ocp_ref
            + delta_phi * self.param.potential_scale,
        }

        return variables

    def _get_electrolyte_overpotentials(self, variables): # Mark Ruihe change the whole thing
        """
        A private function to obtain the electrolyte overpotential and Ohmic losses.
        Note: requires 'variables' to contain the potential, electrolyte concentration
        and temperature the subdomains: 'negative electrode', 'separator', and
        'positive electrode'.

        Parameters
        ----------
        variables : dict
            The variables that have been set in the rest of the model.

        Returns
        -------
        variables : dict
            The variables including the whole-cell electrolyte potentials
            and currents.
        """
        param = self.param
        phi_e_p = variables["Positive electrolyte potential"]

        c_e_s = variables["Separator electrolyte concentration"]
        c_EC_s = variables["Separator EC concentration"]# Mark: Ruihe add
        c_e_p = variables["Positive electrolyte concentration"]
        c_EC_p = variables["Positive EC concentration"]# Mark: Ruihe add

        T_s = variables["Separator temperature"]
        T_p = variables["Positive electrode temperature"]
        phi_e_n = variables["Negative electrolyte potential"]
        c_e_n = variables["Negative electrolyte concentration"]
        c_EC_n = variables["Negative EC concentration"] # Mark: Ruihe add
        T_n = variables["Negative electrode temperature"]

        if self.options.electrode_types["negative"] == "planar":
            # No concentration overpotential in the counter electrode
            phi_e_n = pybamm.Scalar(0)
            indef_integral_n = pybamm.Scalar(0)
            # concentration overpotential
            indef_integral_s = pybamm.IndefiniteIntegral(
                param.chiT_over_c(c_e_s, c_EC_s, T_s) * pybamm.grad(c_e_s),
                pybamm.standard_spatial_vars.x_s,
            )
            indef_integral_p = pybamm.IndefiniteIntegral(
                param.chiT_over_c(c_e_p,c_EC_p, T_p) * pybamm.grad(c_e_p),
                pybamm.standard_spatial_vars.x_p,
            )
            eta_cEC_av = pybamm.Scalar(0)                               # Mark Ruihe add
        elif self.options["electrolyte conductivity"] == "sol full":    # Mark Ruihe add    
            # concentration overpotential
            indef_integral_n = pybamm.IndefiniteIntegral(
                param.chiT_over_c(c_e_n, c_EC_n,T_n) * pybamm.grad(c_e_n) ,#* param.c_0_back / param.ce_tot
                pybamm.standard_spatial_vars.x_n,
            )
            # concentration overpotential
            indef_integral_s = pybamm.IndefiniteIntegral(
                param.chiT_over_c(c_e_s, c_EC_s, T_s) * pybamm.grad(c_e_s) ,#* param.c_0_back / param.ce_tot
                pybamm.standard_spatial_vars.x_s,
            )
            indef_integral_p = pybamm.IndefiniteIntegral(
                param.chiT_over_c(c_e_p,c_EC_p, T_p) * pybamm.grad(c_e_p) ,#* param.c_0_back / param.ce_tot
                pybamm.standard_spatial_vars.x_p,
            )
            # Mark Ruihe add for c(EC) induced concentration overpotential
            indef_integral_EC_n = pybamm.IndefiniteIntegral((                                 
                2 * param.TDF_EC(c_e_n,c_EC_n, T_n)  # * param.c_0_back / param.ce_tot  
                * param.Xi_0 * param.c_ec_typ / param.c_ec_0_dim
                * pybamm.grad(c_EC_n)   ),pybamm.standard_spatial_vars.x_n,
            )
            indef_integral_EC_s = pybamm.IndefiniteIntegral((                                 
                2 * param.TDF_EC(c_e_s,c_EC_s, T_s)  #  * param.c_0_back / param.ce_tot  
                * param.Xi_0 * param.c_ec_typ / param.c_ec_0_dim
                * pybamm.grad(c_EC_s)   ),pybamm.standard_spatial_vars.x_s,
            )
            indef_integral_EC_p = pybamm.IndefiniteIntegral((                                 
                2 * param.TDF_EC(c_e_p,c_EC_p, T_p)  #  * param.c_0_back / param.ce_tot  
                * param.Xi_0 * param.c_ec_typ / param.c_ec_0_dim
                * pybamm.grad(c_EC_p)   ),pybamm.standard_spatial_vars.x_p,
            )
            # process here:
            integral_EC_n = indef_integral_EC_n
            integral_EC_s = indef_integral_EC_s + pybamm.boundary_value(integral_EC_n, "right")
            integral_EC_p = indef_integral_EC_p + pybamm.boundary_value(integral_EC_s, "right")
            eta_cEC_av = pybamm.x_average(integral_EC_p) - pybamm.x_average(integral_EC_n)
        else:
            # concentration overpotential
            indef_integral_n = pybamm.IndefiniteIntegral(
                param.chiT_over_c(c_e_n, c_EC_n,T_n) * pybamm.grad(c_e_n),
                pybamm.standard_spatial_vars.x_n,
            )
            # concentration overpotential
            indef_integral_s = pybamm.IndefiniteIntegral(
                param.chiT_over_c(c_e_s, c_EC_s, T_s) * pybamm.grad(c_e_s),
                pybamm.standard_spatial_vars.x_s,
            )
            indef_integral_p = pybamm.IndefiniteIntegral(
                param.chiT_over_c(c_e_p,c_EC_p, T_p) * pybamm.grad(c_e_p),
                pybamm.standard_spatial_vars.x_p,
            )
            indef_integral_n_fake = pybamm.IndefiniteIntegral(
                pybamm.grad(pybamm.FullBroadcast(
                pybamm.Scalar(0), "negative electrode", 
                auxiliary_domains={"secondary": "current collector"})),
                pybamm.standard_spatial_vars.x_n,
            )
            indef_integral_p_fake = pybamm.IndefiniteIntegral(
                pybamm.grad(pybamm.FullBroadcast(
                pybamm.Scalar(0), "positive electrode", 
                auxiliary_domains={"secondary": "current collector"})),
                pybamm.standard_spatial_vars.x_p,
            )
            eta_cEC_av = (
                pybamm.x_average(indef_integral_n_fake)
                -
                pybamm.x_average(indef_integral_p_fake) ) # Mark Ruihe add

        
        integral_n = indef_integral_n
        integral_s = indef_integral_s + pybamm.boundary_value(integral_n, "right")
        integral_p = indef_integral_p + pybamm.boundary_value(integral_s, "right")

        eta_c_av = pybamm.x_average(integral_p) - pybamm.x_average(integral_n)

        delta_phi_e_av = (
            pybamm.x_average(phi_e_p) - pybamm.x_average(phi_e_n) - eta_c_av - eta_cEC_av # Mark Ruihe change 
        )
        # print(type(eta_c_av),type(eta_cEC_av),type(delta_phi_e_av))
        # Mark Ruihe change 
        #print(len(eta_c_av),len(eta_cEC_av),len(delta_phi_e_av))
        variables.update(   self._get_split_overpotential( eta_c_av, eta_cEC_av, delta_phi_e_av)   ) # 

        return variables

    def set_boundary_conditions(self, variables):
        phi_e = variables["Electrolyte potential"]

        if self.options.electrode_types["negative"] == "planar":
            phi_e_ref = variables["Lithium metal interface electrolyte potential"]
            lbc = (phi_e_ref, "Dirichlet")
        else:
            lbc = (pybamm.Scalar(0), "Neumann")
        self.boundary_conditions = {
            phi_e: {"left": lbc, "right": (pybamm.Scalar(0), "Neumann")}
        }

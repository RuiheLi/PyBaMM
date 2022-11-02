class LithiumIonParameters(BaseParameters):
    """
    Standard parameters for lithium-ion battery models

    Layout:
        1. Dimensional Parameters
        2. Dimensional Functions
        3. Scalings
        4. Dimensionless Parameters
        5. Dimensionless Functions
        6. Input Current

    Parameters
    ----------

    options : dict, optional
        A dictionary of options to be passed to the parameters. The options that
        can be set are listed below.

            * "particle shape" : str, optional
                Sets the model shape of the electrode particles. This is used to
                calculate the surface area to volume ratio. Can be "spherical"
                (default). TODO: implement "cylindrical" and "platelet".
            * "working electrode": str
                Which electrode(s) intercalates and which is counter. If "both"
                (default), the model is a standard battery. Otherwise can be "negative"
                or "positive" to indicate a half-cell model.

    """

    def __init__(self, options=None):
        self.options = options

        # Get geometric, electrical and thermal parameters
        self.geo = pybamm.GeometricParameters(options)
        self.elec = pybamm.electrical_parameters
        self.therm = pybamm.thermal_parameters

        # Initialize domain parameters
        self.n = DomainLithiumIonParameters("negative", self)
        self.s = DomainLithiumIonParameters("separator", self)
        self.p = DomainLithiumIonParameters("positive", self)
        self.domain_params = {
            "negative": self.n,
            "separator": self.s,
            "positive": self.p,
        }

        # Set parameters and scales
        self._set_dimensional_parameters()
        self._set_scales()
        self._set_dimensionless_parameters()

        # Set input current
        self._set_input_current()

    def _set_dimensional_parameters(self):
        """Defines the dimensional parameters"""
        # Physical constants
        self.R = pybamm.constants.R
        self.F = pybamm.constants.F
        self.k_b = pybamm.constants.k_b
        self.q_e = pybamm.constants.q_e

        self.T_ref = self.therm.T_ref
        self.T_init_dim = self.therm.T_init_dim
        self.T_init = self.therm.T_init

        # Macroscale geometry
        self.L_x = self.geo.L_x
        self.L = self.geo.L
        self.L_y = self.geo.L_y
        self.L_z = self.geo.L_z
        self.r_inner_dimensional = self.geo.r_inner_dimensional
        self.r_outer_dimensional = self.geo.r_outer_dimensional
        self.A_cc = self.geo.A_cc
        self.A_cooling = self.geo.A_cooling
        self.V_cell = self.geo.V_cell

        # Electrical
        self.I_typ = self.elec.I_typ
        self.Q = self.elec.Q
        self.C_rate = self.elec.C_rate
        self.n_electrodes_parallel = self.elec.n_electrodes_parallel
        self.n_cells = self.elec.n_cells
        self.i_typ = self.elec.i_typ
        self.voltage_low_cut_dimensional = self.elec.voltage_low_cut_dimensional
        self.voltage_high_cut_dimensional = self.elec.voltage_high_cut_dimensional

        # Domain parameters
        for domain in self.domain_params.values():
            domain._set_dimensional_parameters()

        # Electrolyte properties
        self.c_e_typ = pybamm.Parameter("Typical electrolyte concentration [mol.m-3]")

        self.epsilon_init = pybamm.concatenation(
            *[
                self.domain_params[domain.split()[0]].epsilon_init
                for domain in self.options.whole_cell_domains
            ]
        )
 
        # Lithium plating parameters
        self.V_bar_plated_Li = pybamm.Parameter(
            "Lithium metal partial molar volume [m3.mol-1]"
        )
        self.c_Li_typ = pybamm.Parameter(
            "Typical plated lithium concentration [mol.m-3]"
        )
        self.c_plated_Li_0_dim = pybamm.Parameter(
            "Initial plated lithium concentration [mol.m-3]"
        )

        # Initial conditions
        # Note: the initial concentration in the electrodes can be set as a function
        # of through-cell position, so is defined later as a function
        self.c_e_init_dimensional = pybamm.Parameter(
            "Initial concentration in electrolyte [mol.m-3]"
        )

        self.alpha_T_cell_dim = pybamm.Parameter(
            "Cell thermal expansion coefficient [m.K-1]"
        )

        # Total lithium
        # Electrolyte
        c_e_av_init = pybamm.xyz_average(self.epsilon_init) * self.c_e_typ
        self.n_Li_e_init = c_e_av_init * self.L_x * self.A_cc

        self.n_Li_particles_init = self.n.n_Li_init + self.p.n_Li_init
        self.n_Li_init = self.n_Li_particles_init + self.n_Li_e_init

        # Reference OCP based on initial concentration
        self.ocv_ref = self.p.U_ref - self.n.U_ref
        self.ocv_init_dim = self.p.prim.U_init_dim - self.n.prim.U_init_dim

        # Mark Ruihe block start
        self.Xi = pybamm.Parameter(
            "EC transference number") 
        """ self.c_ec_0_dim = pybamm.Parameter(
            "EC initial concentration in electrolyte [mol.m-3]")  """
        self.c_ec_typ = pybamm.Parameter("Typical EC concentration [mol.m-3]") 
        self.D_ec_Li_cross_dim = pybamm.Parameter(
            "EC Lithium ion cross diffusivity [m2.s-1]")
        self.D_ec_Li_cross_typ = pybamm.Parameter(
            "Typical EC Lithium ion cross diffusivity [m2.s-1]") 
        self.D_ec_dim = pybamm.Parameter(
            "EC diffusivity in electrolyte [m2.s-1]") #may want to be a function in the future
        self.D_ec_typ = pybamm.Parameter(
            "Typical EC diffusivity in electrolyte [m2.s-1]")
        
        self.Xi = pybamm.Parameter("EC transference number") 
        self.Vmolar_ec = pybamm.Parameter(
            "EC partial molar volume [m3.mol-1]")
        self.Vmolar_Li = pybamm.Parameter(
            "Li partial molar volume [m3.mol-1]")
        self.Vmolar_CH2OCO2Li2 = pybamm.Parameter(
            "CH2OCO2Li2 partial molar volume [m3.mol-1]")
        # Total electrolyte concentration [mol.m-3] - just for double diffusion
        self.ce_tot = self.c_e_typ *2 + self.c_ec_typ
        # Mark Ruihe block end

    def D_e_dimensional(self, c_e, c_EC , T):
        """Dimensional diffusivity in electrolyte"""
        tol = pybamm.settings.tolerances["D_e__c_e"]
        c_e = pybamm.maximum(c_e, tol)
        inputs = {
            "Electrolyte concentration [mol.m-3]": c_e, 
            "EC concentration [mol.m-3]": c_EC, 
            "Temperature [K]": T}
        return pybamm.FunctionParameter("Electrolyte diffusivity [m2.s-1]", inputs)

    def kappa_e_dimensional(self, c_e,c_EC , T):
        """Dimensional electrolyte conductivity"""
        tol = pybamm.settings.tolerances["D_e__c_e"]
        c_e = pybamm.maximum(c_e, tol)
        inputs = {
            "Electrolyte concentration [mol.m-3]": c_e, 
            "EC concentration [mol.m-3]": c_EC, 
            "Temperature [K]": T}
        return pybamm.FunctionParameter("Electrolyte conductivity [S.m-1]", inputs)

    def j0_stripping_dimensional(self, c_e, c_Li, T):
        """Dimensional exchange-current density for stripping [A.m-2]"""
        inputs = {
            "Electrolyte concentration [mol.m-3]": c_e,
            "Plated lithium concentration [mol.m-3]": c_Li,
            "Temperature [K]": T,
        }
        return pybamm.FunctionParameter(
            "Exchange-current density for stripping [A.m-2]", inputs
        )

    def j0_plating_dimensional(self, c_e, c_Li, T):
        """Dimensional exchange-current density for plating [A.m-2]"""
        inputs = {
            "Electrolyte concentration [mol.m-3]": c_e,
            "Plated lithium concentration [mol.m-3]": c_Li,
            "Temperature [K]": T,
        }
        return pybamm.FunctionParameter(
            "Exchange-current density for plating [A.m-2]", inputs
        )
    
    def dead_lithium_decay_rate_dimensional(self, L_sei):
        """Dimensional dead lithium decay rate [s-1]"""
        inputs = {"Total SEI thickness [m]": L_sei}
        return pybamm.FunctionParameter("Dead lithium decay rate [s-1]", inputs)

    def _set_scales(self):
        """Define the scales used in the non-dimensionalisation scheme"""
        # Concentration
        self.electrolyte_concentration_scale = self.c_e_typ

        # Electrical
        # Both potential scales are the same but they have different units
        self.potential_scale = self.R * self.T_ref / self.F  # volts
        self.potential_scale_eV = self.k_b / self.q_e * self.T_ref  # eV
        self.current_scale = self.i_typ
        self.current_scale.print_name = "I_typ"

        # Thermal
        self.Delta_T = self.therm.Delta_T

        # Velocity scale
        self.velocity_scale = pybamm.Scalar(1)

        # Discharge timescale
        if self.options["working electrode"] == "positive":
            self.c_max = self.p.prim.c_max
        else:
            self.c_max = self.n.prim.c_max
        self.tau_discharge = self.F * self.c_max * self.L_x / self.i_typ

        # Electrolyte diffusion timescale
        self.D_e_typ = self.D_e_dimensional(self.c_e_typ,self.c_ec_typ, self.T_ref)
        self.tau_diffusion_e = self.L_x ** 2 / self.D_e_typ

        # Thermal diffusion timescale
        self.tau_th_yz = self.therm.tau_th_yz

        # Choose discharge timescale
        if self.options["timescale"] == "default":
            self.timescale = self.tau_discharge
        else:
            self.timescale = pybamm.Scalar(self.options["timescale"])

        for domain in self.domain_params.values():
            domain._set_scales()

    def _set_dimensionless_parameters(self):
        """Defines the dimensionless parameters"""
        # Timescale ratios
        self.C_e = self.tau_diffusion_e / self.timescale
        self.C_th = self.tau_th_yz / self.timescale

        # Concentration ratios
        self.gamma_e = (self.tau_discharge / self.timescale) * self.c_e_typ / self.c_max

        # Macroscale Geometry
        self.l_x = self.geo.l_x
        self.l_y = self.geo.l_y
        self.l_z = self.geo.l_z
        self.r_inner = self.geo.r_inner
        self.r_outer = self.geo.r_outer
        self.a_cc = self.geo.a_cc
        self.a_cooling = self.geo.a_cooling
        self.v_cell = self.geo.v_cell
        self.l = self.geo.l
        self.delta = self.geo.delta

        for domain in self.domain_params.values():
            domain._set_dimensionless_parameters()

        # Electrolyte Properties
        self.beta_surf = pybamm.Scalar(0)
        
        # Mark Ruihe block start
        self.D_ec_Li_cross = self.D_ec_Li_cross_dim  / self.D_ec_Li_cross_typ 
        self.D_ec = self.D_ec_dim / self.D_ec_typ
        self.gamma_e_ec_Rio = self.c_ec_typ / self.c_e_typ
        self.tau_cross_Rio = self.L_x *  self.L_x / self.D_ec_Li_cross_typ
        self.tau_ec_Rio = self.L_x *  self.L_x / self.D_ec_typ
        self.EC_ratio_Rio = self.c_ec_typ / self.ce_tot
        self.e_ratio_Rio = self.c_e_typ / self.ce_tot  
        # Mark Ruihe block end

        # Electrical
        self.voltage_low_cut = (
            self.voltage_low_cut_dimensional - self.ocv_ref
        ) / self.potential_scale
        self.voltage_high_cut = (
            self.voltage_high_cut_dimensional - self.ocv_ref
        ) / self.potential_scale

        # Thermal
        self.Theta = self.therm.Theta
        self.T_amb_dim = self.therm.T_amb_dim
        self.T_amb = self.therm.T_amb

        self.h_edge = self.therm.h_edge
        self.h_total = self.therm.h_total
        self.rho = self.therm.rho

        self.B = (
            self.i_typ
            * self.R
            * self.T_ref
            * self.tau_th_yz
            / (self.therm.rho_eff_dim_ref * self.F * self.Delta_T * self.L_x)
        )

        # lithium plating parameters
        self.c_plated_Li_0 = self.c_plated_Li_0_dim / self.c_Li_typ

        self.alpha_plating = pybamm.Parameter("Lithium plating transfer coefficient")
        self.alpha_stripping = 1 - self.alpha_plating

        # ratio of lithium plating reaction scaled to intercalation reaction
        self.Gamma_plating = (
            self.n.prim.a_typ * self.n.prim.j_scale * self.timescale
        ) / (self.F * self.c_Li_typ)

        # Initial conditions
        self.c_e_init = self.c_e_init_dimensional / self.c_e_typ
        self.ocv_init = (self.ocv_init_dim - self.ocv_ref) / self.potential_scale

        # Dimensionless mechanical parameters
        self.t0_cr = 3600 / (self.C_rate * self.timescale)

    def chi(self, c_e, c_EC, T):
        """
        Thermodynamic factor:
            (1-2*t_plus) is for Nernst-Planck,
            2*(1-t_plus) for Stefan-Maxwell,
        see Bizeray et al (2016) "Resolving a discrepancy ...".
        """
        return (2 * (1 - self.t_plus(c_e,c_EC, T))) * (self.one_plus_dlnf_dlnc(c_e, T))

    def chiT_over_c(self, c_e,c_EC, T):
        """
        chi * (1 + Theta * T) / c,
        as it appears in the electrolyte potential equation
        """
        tol = pybamm.settings.tolerances["chi__c_e"]
        c_e = pybamm.maximum(c_e, tol)
        return self.chi(c_e, c_EC, T) * (1 + self.Theta * T) / c_e

    def t_plus(self, c_e, c_EC, T):
        """Cation transference number (dimensionless)"""
        inputs = {
            "Electrolyte concentration [mol.m-3]": c_e * self.c_e_typ,
            "EC concentration [mol.m-3]": c_EC * self.c_ec_typ,
            "Temperature [K]": self.Delta_T * T + self.T_ref,
        }
        return pybamm.FunctionParameter("Cation transference number", inputs)

    def one_plus_dlnf_dlnc(self, c_e, T):
        """Thermodynamic factor (dimensionless)"""
        inputs = {
            "Electrolyte concentration [mol.m-3]": c_e * self.c_e_typ,
            "Temperature [K]": self.Delta_T * T + self.T_ref,
        }
        return pybamm.FunctionParameter("1 + dlnf/dlnc", inputs)

    def D_e(self, c_e,c_EC, T):
        """Dimensionless electrolyte diffusivity"""
        c_e_dimensional = c_e * self.c_e_typ
        c_EC_dimensional = c_EC * self.c_ec_typ
        T_dim = self.Delta_T * T + self.T_ref
        return self.D_e_dimensional(
            c_e_dimensional, c_EC_dimensional,T_dim) / self.D_e_typ

    def kappa_e(self, c_e,c_EC, T):
        """Dimensionless electrolyte conductivity"""
        c_e_dimensional = c_e * self.c_e_typ
        c_EC_dimensional = c_EC * self.c_ec_typ
        kappa_scale = self.F**2 * self.D_e_typ * self.c_e_typ / (self.R * self.T_ref)
        T_dim = self.Delta_T * T + self.T_ref
        return self.kappa_e_dimensional(
            c_e_dimensional, c_EC_dimensional,T_dim) / kappa_scale

    def j0_stripping(self, c_e, c_Li, T):
        """Dimensionless exchange-current density for stripping"""
        c_e_dim = c_e * self.c_e_typ
        c_Li_dim = c_Li * self.c_Li_typ
        T_dim = self.Delta_T * T + self.T_ref

        return (
            self.j0_stripping_dimensional(c_e_dim, c_Li_dim, T_dim)
            / self.n.prim.j_scale
        )

    def j0_plating(self, c_e, c_Li, T):
        """Dimensionless reverse plating current"""
        c_e_dim = c_e * self.c_e_typ
        c_Li_dim = c_Li * self.c_Li_typ
        T_dim = self.Delta_T * T + self.T_ref

        return (
            self.j0_plating_dimensional(c_e_dim, c_Li_dim, T_dim) / self.n.prim.j_scale
        )

    def dead_lithium_decay_rate(self, L_sei):
        """Dimensionless exchange-current density for stripping"""
        L_sei_dim = L_sei * self.n.prim.L_sei_0_dim

        return self.dead_lithium_decay_rate_dimensional(L_sei_dim) * self.timescale

    def _set_input_current(self):
        """Set the input current"""

        self.dimensional_current_with_time = pybamm.FunctionParameter(
            "Current function [A]", {"Time [s]": pybamm.t * self.timescale}
        )
        self.dimensional_current_density_with_time = (
            self.dimensional_current_with_time
            / (self.n_electrodes_parallel * self.geo.A_cc)
        )
        self.current_with_time = (
            self.dimensional_current_with_time / self.I_typ * pybamm.sign(self.I_typ)
        )

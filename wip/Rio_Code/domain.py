class DomainLithiumIonParameters(BaseParameters):
    def __init__(self, domain, main_param):
        self.domain = domain
        self.main_param = main_param

        self.geo = getattr(main_param.geo, domain.lower()[0])
        self.therm = getattr(main_param.therm, domain.lower()[0])

        if domain != "separator":
            self.prim = ParticleLithiumIonParameters("primary", self)
            phases_option = int(getattr(main_param.options, domain)["particle phases"])
            if phases_option >= 2:
                self.sec = ParticleLithiumIonParameters("secondary", self)
            else:
                self.sec = NullParameters()
        else:
            self.prim = NullParameters()
            self.sec = NullParameters()
        
        self.phase_params = {"primary": self.prim, "secondary": self.sec}


    def _set_dimensional_parameters(self):
        main = self.main_param
        domain, Domain = self.domain_Domain

        if domain == "separator":
            x = pybamm.standard_spatial_vars.x_s * main.L_x
            self.epsilon_init = pybamm.FunctionParameter(
                "Separator porosity", {"Through-cell distance (x) [m]": x}
            )
            self.epsilon_inactive = 1 - self.epsilon_init
            self.b_e = self.geo.b_e
            self.L = self.geo.L
            return

        x = (
            pybamm.SpatialVariable(
                f"x_{domain[0]}",
                domain=[f"{domain} electrode"],
                auxiliary_domains={"secondary": "current collector"},
                coord_sys="cartesian",
            )
            * main.L_x
        )

        # Macroscale geometry
        self.L_cc = self.geo.L_cc
        self.L = self.geo.L

        for phase in self.phase_params.values():
            phase._set_dimensional_parameters()

        # Tab geometry (for pouch cells)
        self.L_tab = self.geo.L_tab
        self.Centre_y_tab = self.geo.Centre_y_tab
        self.Centre_z_tab = self.geo.Centre_z_tab
        self.A_tab = self.geo.A_tab

        # Particle properties
        self.sigma_cc_dimensional = pybamm.Parameter(
            f"{Domain} current collector conductivity [S.m-1]"
        )
        if main.options.electrode_types[domain] == "porous":
            self.epsilon_init = pybamm.FunctionParameter(
                f"{Domain} electrode porosity", {"Through-cell distance (x) [m]": x}
            )
            epsilon_s_tot = sum(phase.epsilon_s for phase in self.phase_params.values())
            self.epsilon_inactive = 1 - self.epsilon_init - epsilon_s_tot

            self.cap_init = sum(phase.cap_init for phase in self.phase_params.values())
            # Use primary phase to set the reference potential
            self.U_ref = self.prim.U_dimensional(self.prim.c_init_av, main.T_ref)
        else:
            self.U_ref = pybamm.Scalar(0)

        self.n_Li_init = sum(phase.n_Li_init for phase in self.phase_params.values())

        # Tortuosity parameters
        self.b_e = self.geo.b_e
        self.b_s = self.geo.b_s

        self.C_dl_dimensional = pybamm.Parameter(
            f"{Domain} electrode double-layer capacity [F.m-2]"
        )

        # Mechanical parameters
        self.nu = pybamm.Parameter(f"{Domain} electrode Poisson's ratio")
        self.E = pybamm.Parameter(f"{Domain} electrode Young's modulus [Pa]")
        self.c_0_dim = pybamm.Parameter(
            f"{Domain} electrode reference concentration for free of deformation "
            "[mol.m-3]"
        )
        self.Omega = pybamm.Parameter(
            f"{Domain} electrode partial molar volume [m3.mol-1]"
        )
        self.l_cr_0 = pybamm.Parameter(f"{Domain} electrode initial crack length [m]")
        self.w_cr = pybamm.Parameter(f"{Domain} electrode initial crack width [m]")
        self.rho_cr_dim = pybamm.Parameter(
            f"{Domain} electrode number of cracks per unit area [m-2]"
        )
        self.b_cr = pybamm.Parameter(f"{Domain} electrode Paris' law constant b")
        self.m_cr = pybamm.Parameter(f"{Domain} electrode Paris' law constant m")
        self.Eac_cr = pybamm.Parameter(
            f"{Domain} electrode activation energy for cracking rate [kJ.mol-1]"
        )
        # intermediate variables  [K*m^3/mol]
        self.theta_dim = (
            (self.Omega / main.R) * 2 * self.Omega * self.E / 9 / (1 - self.nu)
        )

        # Loss of active material parameters
        self.m_LAM = pybamm.Parameter(
            f"{Domain} electrode LAM constant exponential term"
        )
        self.beta_LAM_dimensional = pybamm.Parameter(
            f"{Domain} electrode LAM constant proportional term [s-1]"
        )
        self.stress_critical_dim = pybamm.Parameter(
            f"{Domain} electrode critical stress [Pa]"
        )
        self.beta_LAM_sei_dimensional = pybamm.Parameter(
            f"{Domain} electrode reaction-driven LAM factor [m3.mol-1]"
        )

        # utilisation parameters
        self.u_init = pybamm.Parameter(
            f"Initial {domain} electrode interface utilisation"
        )
        self.beta_utilisation_dimensional = pybamm.Parameter(
            f"{Domain} electrode current-driven interface utilisation factor [m3.mol-1]"
        )

    def sigma_dimensional(self, T):
        """Dimensional electrical conductivity in electrode"""
        inputs = {"Temperature [K]": T}
        Domain = self.domain.capitalize()
        return pybamm.FunctionParameter(
            f"{Domain} electrode conductivity [S.m-1]", inputs
        )

    def _set_scales(self):
        """Define the scales used in the non-dimensionalisation scheme"""
        for phase in self.phase_params.values():
            phase._set_scales()

        if self.domain == "separator":
            return

    def _set_dimensionless_parameters(self):
        for phase in self.phase_params.values():
            phase._set_dimensionless_parameters()

        main = self.main_param

        if self.domain == "separator":
            self.l = self.geo.l
            self.rho = self.therm.rho
            self.lambda_ = self.therm.lambda_
            return

        # Macroscale Geometry
        self.l_cc = self.geo.l_cc
        self.l = self.geo.l

        # Thermal
        self.rho_cc = self.therm.rho_cc
        self.rho = self.therm.rho
        self.lambda_cc = self.therm.lambda_cc
        self.lambda_ = self.therm.lambda_
        self.h_tab = self.therm.h_tab
        self.h_cc = self.therm.h_cc

        # Tab geometry (for pouch cells)
        self.l_tab = self.geo.l_tab
        self.centre_y_tab = self.geo.centre_y_tab
        self.centre_z_tab = self.geo.centre_z_tab

        # Electrochemical Reactions
        self.C_dl = (
            self.C_dl_dimensional
            * main.potential_scale
            / self.prim.j_scale
            / main.timescale
        )
        # Electrode Properties
        self.sigma_cc = (
            self.sigma_cc_dimensional * main.potential_scale / main.i_typ / main.L_x
        )
        self.sigma_cc_prime = self.sigma_cc * main.delta**2
        self.sigma_cc_dbl_prime = self.sigma_cc_prime * main.delta

        # Electrolyte Properties
        self.beta_surf = pybamm.Scalar(0)

        # Utilisation factors
        self.beta_utilisation = (
            self.beta_utilisation_dimensional
            * self.prim.a_typ
            * self.prim.j_scale
            * main.timescale
        ) / main.F

        if main.options.electrode_types[self.domain] == "planar":
            return

        # Dimensionless mechanical parameters
        self.rho_cr = self.rho_cr_dim * self.l_cr_0 * self.w_cr
        self.theta = self.theta_dim * self.prim.c_max / main.T_ref
        self.c_0 = self.c_0_dim / self.prim.c_max
        self.beta_LAM = self.beta_LAM_dimensional * main.timescale
        # normalised typical time for one cycle
        self.stress_critical = self.stress_critical_dim / self.E
        # Reaction-driven LAM parameters
        self.beta_LAM_sei = (
            self.beta_LAM_sei_dimensional
            * self.prim.a_typ
            * self.prim.j_scale
            * main.timescale
        ) / main.F

    def sigma(self, T):
        """Dimensionless electrode electrical conductivity"""
        main = self.main_param
        T_dim = self.main_param.Delta_T * T + self.main_param.T_ref
        return (
            self.sigma_dimensional(T_dim) * main.potential_scale / main.i_typ / main.L_x
        )

    def sigma_prime(self, T):
        """Rescaled dimensionless electrode electrical conductivity"""
        return self.sigma(T) * self.main_param.delta

    def k_cr(self, T):
        """
        Dimensionless cracking rate for the electrode;
        """
        Domain = self.domain.capitalize()
        T_dim = self.main_param.Delta_T * T + self.main_param.T_ref
        delta_k_cr = self.E**self.m_cr * self.l_cr_0 ** (self.m_cr / 2 - 1)
        return (
            pybamm.FunctionParameter(
                f"{Domain} electrode cracking rate", {"Temperature [K]": T_dim}
            )
            * delta_k_cr
        )


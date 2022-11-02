class ParticleLithiumIonParameters(BaseParameters):
    def __init__(self, phase, domain_param):
        self.domain_param = domain_param
        self.domain = domain_param.domain
        self.main_param = domain_param.main_param
        self.phase = phase
        self.set_phase_name()
        if self.phase == "primary":
            self.geo = domain_param.geo.prim
        elif self.phase == "secondary":
            self.geo = domain_param.geo.sec

    def _set_dimensional_parameters(self):
        main = self.main_param
        domain, Domain = self.domain_Domain
        phase_name = self.phase_name
        pref = self.phase_prefactor

        # Electrochemical reactions
        self.ne = pybamm.Parameter(f"{pref}{Domain} electrode electrons in reaction")

        # Intercalation kinetics
        self.mhc_lambda_dimensional = pybamm.Parameter(
            f"{pref}{Domain} electrode reorganization energy [eV]"
        )
        self.alpha_bv = pybamm.Parameter(
            f"{pref}{Domain} electrode Butler-Volmer transfer coefficient"
        )

        if self.domain == "negative":
            # SEI parameters
            self.V_bar_inner_dimensional = pybamm.Parameter(
                f"{pref}Inner SEI partial molar volume [m3.mol-1]"
            )
            self.V_bar_outer_dimensional = pybamm.Parameter(
                f"{pref}Outer SEI partial molar volume [m3.mol-1]"
            )

            self.m_sei_dimensional = pybamm.Parameter(
                f"{pref}SEI reaction exchange current density [A.m-2]"
            )

            self.R_sei_dimensional = pybamm.Parameter(f"{pref}SEI resistivity [Ohm.m]")
            self.D_sol_dimensional = pybamm.Parameter(
                f"{pref}Outer SEI solvent diffusivity [m2.s-1]"
            )
            self.c_sol_dimensional = pybamm.Parameter(
                f"{pref}Bulk solvent concentration [mol.m-3]"
            )
            self.U_inner_dimensional = pybamm.Parameter(
                f"{pref}Inner SEI open-circuit potential [V]"
            )
            self.U_outer_dimensional = pybamm.Parameter(
                f"{pref}Outer SEI open-circuit potential [V]"
            )
            self.kappa_inner_dimensional = pybamm.Parameter(
                f"{pref}Inner SEI electron conductivity [S.m-1]"
            )
            self.D_li_dimensional = pybamm.Parameter(
                f"{pref}Inner SEI lithium interstitial diffusivity [m2.s-1]"
            )
            self.c_li_0_dimensional = pybamm.Parameter(
                f"{pref}Lithium interstitial reference concentration [mol.m-3]"
            )
            self.L_inner_0_dim = pybamm.Parameter(
                f"{pref}Initial inner SEI thickness [m]"
            )
            self.L_outer_0_dim = pybamm.Parameter(
                f"{pref}Initial outer SEI thickness [m]"
            )
            self.L_sei_0_dim = self.L_inner_0_dim + self.L_outer_0_dim
            self.E_sei_dimensional = pybamm.Parameter(
                f"{pref}SEI growth activation energy [J.mol-1]"
            )

            # EC reaction
            self.c_ec_0_dim = pybamm.Parameter(
                f"{pref}EC initial concentration in electrolyte [mol.m-3]"
            )
            self.D_ec_sei_dim = pybamm.Parameter(
                f"{pref}EC diffusivity in SEI [m2.s-1]")#may want to be a function in the future
            
            self.D_ec_sei_typ = pybamm.Parameter(
                "Typical EC diffusivity in SEI [m2.s-1]") 
            self.D_ec_sei = self.D_ec_sei_dim / self.D_ec_sei_typ



            self.k_sei_dim = pybamm.Parameter(
                f"{pref}SEI kinetic rate constant [m.s-1]"
            )
            self.U_sei_dim = pybamm.Parameter(
                f"{pref}SEI open-circuit potential [V]")

        if main.options.electrode_types[domain] == "planar":
            self.n_Li_init = pybamm.Scalar(0)
            self.U_init_dim = pybamm.Scalar(0)
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
        r = (
            pybamm.SpatialVariable(
                f"r_{domain[0]}",
                domain=[f"{domain} {self.phase_name}particle"],
                auxiliary_domains={
                    "secondary": f"{domain} electrode",
                    "tertiary": "current collector",
                },
                coord_sys="spherical polar",
            )
            * self.geo.R_typ
        )

        # Macroscale geometry
        # Note: the surface area to volume ratio is defined later with the function
        # parameters. The particle size as a function of through-cell position is
        # already defined in geometric_parameters.py
        self.R_dimensional = self.geo.R_dimensional

        # Particle properties
        self.c_max = pybamm.Parameter(
            f"{pref}Maximum concentration in {domain} electrode [mol.m-3]"
        )

        # Particle-size distribution parameters
        self.R_min_dim = self.geo.R_min_dim
        self.R_max_dim = self.geo.R_max_dim
        self.sd_a_dim = self.geo.sd_a_dim
        self.f_a_dist_dimensional = self.geo.f_a_dist_dimensional

        self.epsilon_s = pybamm.FunctionParameter(
            f"{pref}{Domain} electrode active material volume fraction",
            {"Through-cell distance (x) [m]": x},
        )
        self.c_init = (
            pybamm.FunctionParameter(
                f"{pref}Initial concentration in {domain} electrode [mol.m-3]",
                {
                    "Radial distance (r) [m]": r,
                    "Through-cell distance (x) [m]": pybamm.PrimaryBroadcast(
                        x, f"{domain} {phase_name}particle"
                    ),
                },
            )
            / self.c_max
        )
        self.c_init_av = pybamm.xyz_average(pybamm.r_average(self.c_init))
        eps_c_init_av = pybamm.xyz_average(
            self.epsilon_s * pybamm.r_average(self.c_init)
        )
        self.n_Li_init = eps_c_init_av * self.c_max * self.domain_param.L * main.A_cc

        eps_s_av = pybamm.xyz_average(self.epsilon_s)
        self.elec_loading = eps_s_av * self.domain_param.L * self.c_max * main.F / 3600
        self.cap_init = self.elec_loading * main.A_cc

        self.U_init_dim = self.U_dimensional(self.c_init_av, main.T_init_dim)

    def D_dimensional(self, sto, T):
        """Dimensional diffusivity in particle. Note this is defined as a
        function of stochiometry"""
        Domain = self.domain.capitalize()
        inputs = {
            f"{self.phase_prefactor}{Domain} particle stoichiometry": sto,
            "Temperature [K]": T,
        }
        return pybamm.FunctionParameter(
            f"{self.phase_prefactor}{Domain} electrode diffusivity [m2.s-1]",
            inputs,
        )

    def j0_dimensional(self, c_e, c_s_surf, T):
        """Dimensional exchange-current density [A.m-2]"""
        domain, Domain = self.domain_Domain
        inputs = {
            "Electrolyte concentration [mol.m-3]": c_e,
            f"{Domain} particle surface concentration [mol.m-3]": c_s_surf,
            f"{self.phase_prefactor}Maximum {domain} particle "
            "surface concentration [mol.m-3]": self.c_max,
            "Temperature [K]": T,
            f"{self.phase_prefactor}Maximum {domain} particle "
            "surface concentration [mol.m-3]": self.c_max,
        }
        return pybamm.FunctionParameter(
            f"{self.phase_prefactor}{Domain} electrode "
            "exchange-current density [A.m-2]",
            inputs,
        )

    def U_dimensional(self, sto, T, lithiation=None):
        """Dimensional open-circuit potential [V]"""
        # bound stoichiometry between tol and 1-tol. Adding 1/sto + 1/(sto-1) later
        # will ensure that ocp goes to +- infinity if sto goes into that region
        # anyway
        Domain = self.domain.capitalize()
        tol = pybamm.settings.tolerances["U__c_s"]
        sto = pybamm.maximum(pybamm.minimum(sto, 1 - tol), tol)
        if lithiation is None:
            lithiation = ""
        else:
            lithiation = lithiation + " "
        inputs = {f"{self.phase_prefactor}{Domain} particle stoichiometry": sto}
        u_ref = pybamm.FunctionParameter(
            f"{self.phase_prefactor}{Domain} electrode {lithiation}OCP [V]", inputs
        )
        # add a term to ensure that the OCP goes to infinity at 0 and -infinity at 1
        # this will not affect the OCP for most values of sto
        # see #1435
        u_ref = u_ref + 1e-6 * (1 / sto + 1 / (sto - 1))
        dudt_dim_func = self.dUdT_dimensional(sto)
        d = self.domain[0]
        dudt_dim_func.print_name = r"\frac{dU_{" + d + r"}}{dT}"
        return u_ref + (T - self.main_param.T_ref) * dudt_dim_func

    def dUdT_dimensional(self, sto):
        """
        Dimensional entropic change of the open-circuit potential [V.K-1]
        """
        domain, Domain = self.domain_Domain
        inputs = {
            f"{Domain} particle stoichiometry": sto,
            f"{self.phase_prefactor}Maximum {domain} particle "
            "surface concentration [mol.m-3]": self.c_max,
        }
        return pybamm.FunctionParameter(
            f"{self.phase_prefactor}{Domain} electrode OCP entropic change [V.K-1]",
            inputs,
        )

    def _set_scales(self):
        """Define the scales used in the non-dimensionalisation scheme"""
        domain = self.domain
        main = self.main_param

        # Scale for interfacial current density in A/m2
        if main.options.electrode_types[domain] == "planar":
            # planar electrode (boundary condition between negative and separator)
            self.a_typ = 1
            self.j_scale = main.i_typ
            return

        # Microscale
        self.R_typ = self.geo.R_typ

        if main.options["particle shape"] == "spherical":
            self.a_typ = 3 * pybamm.xyz_average(self.epsilon_s) / self.R_typ

        # porous electrode
        self.j_scale = main.i_typ / (self.a_typ * main.L_x)

        # Concentration
        self.particle_concentration_scale = self.c_max

        # Reference exchange-current density
        self.j0_ref_dimensional = (
            self.j0_dimensional(main.c_e_typ, self.c_max / 2, main.T_ref) * 2
        )

        # Reaction timescales
        self.tau_r = main.F * self.c_max / (self.j0_ref_dimensional * self.a_typ)
        # Particle diffusion timescales
        self.D_typ_dim = self.D_dimensional(pybamm.Scalar(1), main.T_ref)
        self.tau_diffusion = self.R_typ**2 / self.D_typ_dim

    def _set_dimensionless_parameters(self):
        main = self.main_param
        domain_param = self.domain_param
        pref = self.phase_prefactor

        # Intercalation kinetics
        self.mhc_lambda = self.mhc_lambda_dimensional / main.potential_scale_eV

        if self.domain == "negative":
            # SEI parameters
            self.inner_sei_proportion = pybamm.Parameter(
                f"{pref}Inner SEI reaction proportion"
            )

            self.z_sei = pybamm.Parameter(
                f"{pref}Ratio of lithium moles to SEI moles")

            self.E_over_RT_sei = self.E_sei_dimensional / main.R / main.T_ref

            self.C_sei_reaction = (self.j_scale / self.m_sei_dimensional) * pybamm.exp(
                -(main.F * domain_param.U_ref / (2 * main.R * main.T_ref))
            )
            # Mark Ruihe block start
            self.C_sei_solvent = (
                self.j_scale
                * self.L_sei_0_dim
                / (self.c_ec_0_dim * main.F * self.D_ec_sei_dim)
            )  
            # Mark Ruihe : change from D_sol_dimensional to D_ec_sei_dim
            #              change from c_sol_dimensional to c_ec_0_dim
            # Mark Ruihe block end

            self.C_sei_electron = (
                self.j_scale
                * main.F
                * self.L_sei_0_dim
                / (self.kappa_inner_dimensional * main.R * main.T_ref)
            )

            self.C_sei_inter = (
                self.j_scale
                * self.L_sei_0_dim
                / (self.D_li_dimensional * self.c_li_0_dimensional * main.F)
            )

            self.U_inner_electron = (
                main.F * self.U_inner_dimensional / main.R / main.T_ref
            )

            self.R_sei = (
                main.F
                * self.j_scale
                * self.R_sei_dimensional
                * self.L_sei_0_dim
                / main.R
                / main.T_ref
            )

            self.v_bar = self.V_bar_outer_dimensional / self.V_bar_inner_dimensional
            self.c_sei_scale = (
                self.L_sei_0_dim * self.a_typ / self.V_bar_inner_dimensional
            )
            self.c_sei_outer_scale = (
                self.L_sei_0_dim * self.a_typ / self.V_bar_outer_dimensional
            )

            self.L_inner_0 = self.L_inner_0_dim / self.L_sei_0_dim
            self.L_outer_0 = self.L_outer_0_dim / self.L_sei_0_dim

            # Dividing by 10000 makes initial condition effectively zero
            # without triggering division by zero errors
            self.L_inner_crack_0 = self.L_inner_0 / 10000
            self.L_outer_crack_0 = self.L_outer_0 / 10000

            # ratio of SEI reaction scale to intercalation reaction
            self.Gamma_SEI = (
                self.V_bar_inner_dimensional * self.j_scale * main.timescale
            ) / (main.F * self.z_sei * self.L_sei_0_dim)

            # EC reaction
            self.C_ec = (
                self.L_sei_0_dim
                * self.j_scale
                / (main.F * self.c_ec_0_dim * self.D_ec_sei_dim)
            )
            self.C_sei_ec = (
                main.F
                * self.k_sei_dim
                * self.c_ec_0_dim
                / self.j_scale
                * (
                    pybamm.exp(
                        -(
                            main.F
                            * (domain_param.U_ref - self.U_sei_dim)
                            / (2 * main.R * main.T_ref)
                        )
                    )
                )
            )
            self.c_sei_init = self.c_ec_0_dim / self.c_sei_outer_scale

        # Initial conditions
        if main.options.electrode_types[self.domain] == "planar":
            self.U_init = pybamm.Scalar(0)
            return
        else:
            self.U_init = (
                self.U_init_dim - self.domain_param.U_ref
            ) / main.potential_scale

        # Timescale ratios
        self.C_diff = self.tau_diffusion / main.timescale
        self.C_r = self.tau_r / main.timescale

        # Microscale geometry
        self.R = self.geo.R
        self.a_R = self.a_typ * self.R_typ

        # Particle-size distribution geometry
        self.R_min = self.geo.R_min
        self.R_max = self.geo.R_max
        self.sd_a = self.geo.sd_a
        self.f_a_dist = self.geo.f_a_dist

        # Concentration ratios
        # In most cases gamma_n will be equal to 1
        self.gamma = (main.tau_discharge / main.timescale) * self.c_max / main.c_max

        # Electrolyte Properties
        self.beta_surf = pybamm.Scalar(0)

    def D(self, c_s, T):
        """Dimensionless particle diffusivity"""
        sto = c_s
        T_dim = self.main_param.Delta_T * T + self.main_param.T_ref
        return self.D_dimensional(sto, T_dim) / self.D_typ_dim

    def j0(self, c_e, c_s_surf, T):
        """Dimensionless exchange-current density"""
        tol = pybamm.settings.tolerances["j0__c_e"]
        c_e = pybamm.maximum(c_e, tol)
        tol = pybamm.settings.tolerances["j0__c_s"]
        c_s_surf = pybamm.maximum(pybamm.minimum(c_s_surf, 1 - tol), tol)
        c_e_dim = c_e * self.main_param.c_e_typ
        c_s_surf_dim = c_s_surf * self.c_max
        T_dim = self.main_param.Delta_T * T + self.main_param.T_ref

        return self.j0_dimensional(c_e_dim, c_s_surf_dim, T_dim) / self.j_scale

    def U(self, c_s, T, lithiation=None):
        """Dimensionless open-circuit potential in the electrode"""
        main = self.main_param
        sto = c_s
        T_dim = self.main_param.Delta_T * T + self.main_param.T_ref
        return (
            self.U_dimensional(sto, T_dim, lithiation) - self.domain_param.U_ref
        ) / main.potential_scale

    def dUdT(self, c_s):
        """Dimensionless entropic change in open-circuit potential"""
        main = self.main_param
        sto = c_s
        return self.dUdT_dimensional(sto) * main.Delta_T / main.potential_scale

    def t_change(self, sto):
        """
        Dimensionless volume change for the electrode;
        sto should be R-averaged
        """
        domain, Domain = self.domain_Domain
        return pybamm.FunctionParameter(
            f"{Domain} electrode volume change",
            {
                "Particle stoichiometry": sto,
                f"{self.phase_prefactor}Maximum {domain} particle "
                "surface concentration [mol.m-3]": self.c_max,
            },
        )

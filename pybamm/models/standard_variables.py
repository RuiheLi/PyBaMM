#
# Standard variables for the models
#
import pybamm


class StandardVariables:
    def __init__(self):
        # Discharge capacity and energy
        self.Q_Ah = pybamm.Variable("Discharge capacity [A.h]")
        self.Q_Wh = pybamm.Variable("Discharge energy [W.h]")
        # Throughput capacity and energy (cumulative)
        self.Qt_Ah = pybamm.Variable("Throughput capacity [A.h]")
        self.Qt_Wh = pybamm.Variable("Throughput energy [W.h]")
        # Electrolyte concentration
        self.c_e_n = pybamm.Variable(
            "Negative electrolyte concentration",
            domain="negative electrode",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )
        self.c_e_s = pybamm.Variable(
            "Separator electrolyte concentration",
            domain="separator",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )
        self.c_e_p = pybamm.Variable(
            "Positive electrolyte concentration",
            domain="positive electrode",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )

        self.c_e_av = pybamm.Variable(
            "X-averaged electrolyte concentration",
            domain="current collector",
            bounds=(0, np.inf),
        )

        # Electrolyte porosity times concentration
        self.eps_c_e_n = pybamm.Variable(
            "Negative electrode porosity times concentration",
            domain="negative electrode",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )
        self.eps_c_e_s = pybamm.Variable(
            "Separator porosity times concentration",
            domain="separator",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )
        self.eps_c_e_p = pybamm.Variable(
            "Positive electrode porosity times concentration",
            domain="positive electrode",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )

        self.eps_c_e_av = pybamm.Variable(
            "X-averaged porosity times concentration",
            domain="current collector",
            bounds=(0, np.inf),
        )


        # Mark Ruihe block start
        # EC concentration
        self.c_EC_n = pybamm.Variable(
            "Negative EC concentration",
            domain="negative electrode",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )
        self.c_EC_s = pybamm.Variable(
            "Separator EC concentration",
            domain="separator",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )
        self.c_EC_p = pybamm.Variable(
            "Positive EC concentration",
            domain="positive electrode",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )

        self.c_EC_av = pybamm.Variable(
            "X-averaged EC concentration",
            domain="current collector",
            bounds=(0, np.inf),
        )

        # Electrolyte porosity times EC concentration
        self.eps_c_EC_n = pybamm.Variable(
            "Negative electrode porosity times EC concentration",
            domain="negative electrode",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )
        self.eps_c_EC_s = pybamm.Variable(
            "Separator porosity times EC concentration",
            domain="separator",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )
        self.eps_c_EC_p = pybamm.Variable(
            "Positive electrode porosity times EC concentration",
            domain="positive electrode",
            auxiliary_domains={"secondary": "current collector"},
            bounds=(0, np.inf),
        )

        self.eps_c_EC_av = pybamm.Variable(
            "X-averaged porosity times EC concentration",
            domain="current collector",
            bounds=(0, np.inf),
        )
        # Mark Ruihe block end

        # Electrolyte potential
        self.phi_e_n = pybamm.Variable(
            "Negative electrolyte potential",
            domain="negative electrode",
            auxiliary_domains={"secondary": "current collector"},
        )
        self.phi_e_s = pybamm.Variable(
            "Separator electrolyte potential",
            domain="separator",
            auxiliary_domains={"secondary": "current collector"},
        )
        self.phi_e_p = pybamm.Variable(
            "Positive electrolyte potential",
            domain="positive electrode",
            auxiliary_domains={"secondary": "current collector"},
        )

        # Electrode potential
        self.phi_s_n = pybamm.Variable(
            "Negative electrode potential",
            domain="negative electrode",
            auxiliary_domains={"secondary": "current collector"},
        )
        self.phi_s_p = pybamm.Variable(
            "Positive electrode potential",
            domain="positive electrode",
            auxiliary_domains={"secondary": "current collector"},
        )

        # Potential difference
        self.delta_phi_n = pybamm.Variable(
            "Negative electrode surface potential difference",
            domain="negative electrode",
            auxiliary_domains={"secondary": "current collector"},
        )
        self.delta_phi_p = pybamm.Variable(
            "Positive electrode surface potential difference",
            domain="positive electrode",
            auxiliary_domains={"secondary": "current collector"},
        )

        self.delta_phi_n_av = pybamm.Variable(
            "X-averaged negative electrode surface potential difference",
            domain="current collector",
        )
        self.delta_phi_p_av = pybamm.Variable(
            "X-averaged positive electrode surface potential difference",
            domain="current collector",
        )

        # current collector variables
        self.phi_s_cn = pybamm.Variable(
            "Negative current collector potential", domain="current collector"
        )
        self.phi_s_cp = pybamm.Variable(
            "Positive current collector potential", domain="current collector"
        )
        self.i_boundary_cc = pybamm.Variable(
            "Current collector current density", domain="current collector"
        )
        self.phi_s_cn_composite = pybamm.Variable(
            "Composite negative current collector potential", domain="current collector"
        )
        self.phi_s_cp_composite = pybamm.Variable(
            "Composite positive current collector potential", domain="current collector"
        )
        self.i_boundary_cc_composite = pybamm.Variable(
            "Composite current collector current density", domain="current collector"
        )

    def __setattr__(self, name, value):
        value.print_name = name
        super().__setattr__(name, value)


standard_variables = StandardVariables()

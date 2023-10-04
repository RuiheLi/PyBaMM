def calculate_initial_SEI_thickness(V_SEI):
    """
    Initial SEI thickness as a function of partial molar volume for the LG M50 cell

    Parameters
    ----------

    V_SEI: SEI partial molar volume [m3.mol-1]

    Returns
    -------

    Initial inner SEI thickness [m]

    Note: Initial inner and outer SEI thicknesses are assumed to be equal 
    """

    delta_Q_SEI = 3.6 * (5000 - 4864.91)
    F = 96485.3
    A = 1.58 * 0.065
    z_SEI = 2
    L_neg = 8.52e-5
    eps_act_neg = 0.75
    R_neg = 5.86e-6
    l_cr_init = 2e-8
    w_cr = 1.5e-8
    rho_cr = 3.18e15
    roughness = 1 + 2 * l_cr_init * w_cr * rho_cr
    L_SEI_init = delta_Q_SEI * R_neg * V_SEI / (
        z_SEI * F * A * L_neg * 3 * eps_act_neg * roughness
    )
    L_inner_init = L_SEI_init / 2

    return L_inner_init

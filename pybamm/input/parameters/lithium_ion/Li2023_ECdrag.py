import pybamm
from pybamm import exp,sqrt,tanh
# Mark Ruihe Li add for ECdrag2 branch and double diffusion model, mostly based on Chen2020

def graphite_LGM50_ocp_Chen2020(sto):
    """
    LG M50 Graphite open-circuit potential as a function of stochiometry, fit taken
    from [1].

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open circuit potential
    """

    u_eq = (
        1.9793 * pybamm.exp(-39.3631 * sto)
        + 0.2482
        - 0.0909 * pybamm.tanh(29.8538 * (sto - 0.1234))
        - 0.04478 * pybamm.tanh(14.9159 * (sto - 0.2769))
        - 0.0205 * pybamm.tanh(30.4444 * (sto - 0.6103))
    )

    return u_eq


def graphite_LGM50_electrolyte_exchange_current_density_Chen2020(
    c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between graphite and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    m_ref = 6.48e-7  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    E_r = 35000
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return (
        m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5
    )


def nmc_LGM50_ocp_Chen2020(sto):
    """
    LG M50 NMC open circuit potential as a function of stochiometry, fit taken
    from [1].

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open circuit potential
    """

    u_eq = (
        -0.8090 * sto
        + 4.4875
        - 0.0428 * pybamm.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * pybamm.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * pybamm.tanh(15.9308 * (sto - 0.3120))
    )

    return u_eq


def nmc_LGM50_electrolyte_exchange_current_density_Chen2020(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between NMC and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    m_ref = 3.42e-6  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    E_r = 17800
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return (
        m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5
    )


def electrolyte_diffusivity_Nyman2008(c_e,c_EC, T):
    """
    Diffusivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
    comes from [1]

    References
    ----------
    .. [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and
    modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte,"
    Electrochim. Acta, vol. 53, no. 22, pp. 6356–6365, 2008.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    D_c_e = 8.794e-11 * (c_e / 1000) ** 2 - 3.972e-10 * (c_e / 1000) + 4.862e-10

    # Nyman et al. (2008) does not provide temperature dependence

    return D_c_e


def electrolyte_conductivity_Nyman2008(c_e, c_EC,T):
    """
    Conductivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
    comes from [1].

    References
    ----------
    .. [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and
    modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte,"
    Electrochim. Acta, vol. 53, no. 22, pp. 6356–6365, 2008.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    sigma_e = (
        0.1297 * (c_e / 1000) ** 3 - 2.51 * (c_e / 1000) ** 1.5 + 3.329 * (c_e / 1000)
    )

    # Nyman et al. (2008) does not provide temperature dependence

    return sigma_e

def electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs):
    c = (
        (c_e >= 0 ) * c_e / 1000 
        +  (c_e < 0 ) * 0  
        )# mol.m-3 -> mol.l
    p1, p2, p3, p4, p5, p6 = coeffs
    A = p1 * (1 + (T - p2))
    B = 1 + p3 * pybamm.sqrt(c) + p4 * (1 + p5 * pybamm.exp(1000 / T)) * c
    C = 1 + c**4 * (p6 * pybamm.exp(1000 / T))
    sigma_e = A * c * B / C  # mS.cm-1
    return sigma_e / 10
def electrolyte_conductivity_EC_DMC_1_1_Landesfeind2019_Con(c_e,c_EC, T):
    coeffs = np.array([7.98e-1, 2.28e2, -1.22, 5.09e-1, -4e-3, 3.79e-3])
    c_e_constant = 4000
    sigma = (
        (c_e <= c_e_constant) * electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)
        +
        (c_e >  c_e_constant) * electrolyte_conductivity_base_Landesfeind2019(c_e_constant,c_EC, T, coeffs)
    )
    return sigma

def electrolyte_conductivity_Wang2021(c_e,c_EC,T): # ACS Energy Letter
    c_T = Fun_c_T(c_e, c_EC,T)
    y_e = c_e / c_T
    sigma =  (
        48.93*y_e**1.5 - 
        284.8*y_e**2.5 + 
        817.7*y_e**4) **2
    return sigma

# Mark Ruihe block start   - add my revised electrolyte properties
def electrolyte_conductivity_Nyman2008Constant(c_e,c_EC, T):
    
    sigma_e = (c_e <= 2000) * (
        0.1297 * (c_e / 1000) ** 3 - 2.51 * (c_e / 1000) ** 1.5 + 3.329 * (c_e / 1000)
    ) + (c_e > 2000) *  (
        0.1297 * (2000 / 1000) ** 3 - 2.51 * (2000 / 1000) ** 1.5 + 3.329 * (2000 / 1000)
    )
    return sigma_e
def electrolyte_diffusivity_Nyman2008Constant(c_e,c_EC, T):
    D_c_e = (
        (c_e <= 2000) * (
            8.794e-11 * (c_e / 1000) ** 2 - 3.972e-10 * (c_e / 1000) + 4.862e-10) 
        + (c_e > 2000) * (
        8.794e-11 * (2000 / 1000) ** 2 - 3.972e-10 * (2000 / 1000) + 4.862e-10)
    )
    return D_c_e

def EC_diffusivity_5E_10(c_e, c_EC , T):
    # Ruihe add: set Heaviside(c_EC) * 5e-10 
    D_ec_dim = (
        (c_EC >= 0 ) * 5e-10 
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim
def EC_diffusivity_3E_10(c_e, c_EC , T):
    # Ruihe add: set Heaviside(c_EC) * 3e-10 
    D_ec_dim = (
        (c_EC >= 0 ) * 3e-10 
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim

# Update 240125 Main change of this round
def Fun_rho(c_e, c_EC,T):
    rho_0 = 1006.1 
    rho_1 = 0.02235185918895445 
    rho_2 = 0.10065156540490541
    return rho_0  + rho_1 * c_EC + rho_2 * c_e
def Fun_c_T(c_e, c_EC,T):
    m_bar_EC = 88.062*1e-3 #   kg/mol
    m_bar_0  = 104.105*1e-3 #   kg/mol
    m_bar_e  = 151.905*1e-3 #   kg/mol
    b = (
        (m_bar_EC-m_bar_0) * c_EC 
        + (m_bar_e-2*m_bar_0)*c_e  ) 
    c_T= (Fun_rho(c_e, c_EC,T) - b) / m_bar_0
    return c_T

def Fun_D_0_EC_e(c_e, c_EC , D_0_e_EC,T):   # Eq. (41)
    # Pre-condition
    D_0_EC_EC  = 5E-10
    D_0_e_e = 3e-10
    c_T = Fun_c_T(c_e, c_EC,T)
    c_0 = c_T - 2*c_e - c_EC 
    V_bar_EC = 6.5312e-05
    V_bar_0  = 1.0347e-04
    V_bar_e  = 5.0943e-05
    C_11 = 1 + c_e/c_0 + (V_bar_EC*c_EC)/(V_bar_0*c_0)
    C_12 = (V_bar_e-2*V_bar_0)*c_EC / (V_bar_0*c_0)
    C_21 = (V_bar_EC- V_bar_0)*c_e  / (V_bar_0*c_0)
    C_22 = 1 + c_EC/c_0 + (V_bar_e*c_e)/(V_bar_0*c_0)
    D_0_EC_e = (
        D_0_EC_EC*C_12/C_11          + 
        2 * c_EC/c_e * (C_22/C_11 - C_12*C_21/(C_11*C_11)  ) * (
            D_0_e_EC -  (D_0_e_EC*C_12 - D_0_e_e*C_11)*C_21 / (C_21*C_12-C_22*C_11)     
        )
    )
    return D_0_EC_e


def Dimensional_EC_Lithium_ion_cross_diffusivity(c_e, c_EC , T):
    D_0_e_EC = pybamm.Parameter("Lithium ion EC cross diffusivity [m2.s-1]") 
    return Fun_D_0_EC_e(c_e, c_EC , D_0_e_EC,T)

def EC_transference_number_3(c_e,c_EC, T):# Mark Ruihe add update 221212
    Xi_0 =   3.0     # pybamm.Parameter("EC transference number zero") 
    Xi = ( 
        (c_EC < 0 ) * Xi_0
        + 
        (c_EC >= 0) * Xi_0
    )
    return Xi

def Fun_y_e_u_e(y_e, y_EC): # Eq. (12) in Jung 2023 paper
    T  = 298.15
    R = 8.31446261815324; F = 96485.3321
    d_deltaU_dye = 2.30205508980545*y_EC*(-2*y_e - y_EC + 1) - 0.0513851582545859*y_EC*(89.6*y_e + 6.81*y_EC - 12.6) + y_EC*(24.5569671298666*y_e - 2.26403007269705 + 0.0776943592809339/y_e)/(1 - 2*y_e) + 2*y_EC*(12.2784835649333*y_e**2 - 2.26403007269705*y_e + 0.0776943592809339*pybamm.log(y_e) + 0.211527003955003)/(1 - 2*y_e)**2 - 2*y_EC*(-1.10889171513396*y_e**0.5 - 10.4851415418483*y_e**3 + 6.0865719952557*y_e**2 + 4.76340417020011*y_e - 10.3386938408227*y_e**1.5 + 6.51820732459422*y_e**2.5 + 64.462681030378*y_e**3.5 - 74.1487833613674*y_e**4.5 + 0.0301630878954419*pybamm.log(y_e) + 0.184138714605309)/(1 - 2*y_e)**2 + (-y_EC/(1 - 2*y_e) + 1)*(-0.554445857566982/y_e**0.5 - 15.508040761234*y_e**0.5 - 31.4554246255448*y_e**2 + 12.1731439905114*y_e + 16.2955183114856*y_e**1.5 + 225.619383606323*y_e**2.5 - 333.669525126154*y_e**3.5 + 4.76340417020011 + 0.0301630878954419/y_e)

    u_e = F/(R*T) * d_deltaU_dye 
    return u_e * y_e
def Fun_t_0plus_Wang2021_yBased(y_e,y_EC):
    return 0.4107-1.487*y_e+2.547*y_e**2
def Fun_X_e_e_yBased(y_e,y_EC):
    u_e_y_e = Fun_y_e_u_e(y_e, y_EC)
    t_0plus = Fun_t_0plus_Wang2021_yBased(y_e,y_EC)
    return u_e_y_e / (1-t_0plus) - 2
def Fun_X_EC_EC_yBased(y_e,y_EC):
    return 0.5* Fun_X_e_e_yBased(y_e,y_EC)
def Fun_u_EC(y_e, y_EC): # Eq. (12) in Jung 2023 paper
    T  = 298.15
    R = 8.31446261815324; F = 96485.3321
    d_deltaU_dyEC = 0.174966463856865*y_EC*(-2*y_e - y_EC + 1) - 0.0256925791272929*y_EC*(89.6*y_e + 6.81*y_EC - 12.6) + 0.0256925791272929*(-2*y_e - y_EC + 1)*(89.6*y_e + 6.81*y_EC - 12.6) + (12.2784835649333*y_e**2 - 2.26403007269705*y_e + 0.0776943592809339*pybamm.log(y_e) + 0.211527003955003)/(1 - 2*y_e) - (-1.10889171513396*y_e**0.5 - 10.4851415418483*y_e**3 + 6.0865719952557*y_e**2 + 4.76340417020011*y_e - 10.3386938408227*y_e**1.5 + 6.51820732459422*y_e**2.5 + 64.462681030378*y_e**3.5 - 74.1487833613674*y_e**4.5 + 0.0301630878954419*pybamm.log(y_e) + 0.184138714605309)/(1 - 2*y_e)
    u_EC = F/(R*T) * d_deltaU_dyEC   # .subs({x:y_e,  y:y_EC })
    return u_EC
def Fun_Xi_tidle(c_e,c_EC,T):
    T = 298.15
    c_T = Fun_c_T(c_e, c_EC,T)
    y_e = c_e / c_T
    y_EC =c_EC/ c_T
    u_EC = Fun_u_EC(y_e, y_EC)
    X_oo = Fun_X_EC_EC_yBased(y_e,y_EC)
    Xi_tilde = - u_EC / (1+X_oo)
    return Xi_tilde

def t_0plus_linkEC(c_e,c_EC, T):# Mark Ruihe add update 221214
    c_EC_0 = pybamm.Parameter("Typical EC concentration [mol.m-3]")
    t_0plus = ( 
        (c_EC < 0 ) * 0 
        + 
        (c_EC <= c_EC_0) * (c_e >= 0) * (0.34 * c_EC / c_EC_0 )
        +
        (c_EC > c_EC_0 ) * 0.34  
    )
    return t_0plus

def electrolyte_conductivity_Ding2001(c_e, c_EC,  T):
    # c_e is lithium ion concentration in electrolyte in mol/m3, need to change to mol/kg
    # also be careful that T here is deg, while others are K
    rho_electrolyte = 1300 # in kg/m3
    c_e_kg = c_e / rho_electrolyte     # in mol/kg 
    M_LiPF6 = 151.905/1000  # kg/mol
    M_EC = 88.062/1000  # kg/mol
    M_EMC = 104.104/1000 # kg/mol
    x_EC = 1 / (1+ ( rho_electrolyte - c_e*M_LiPF6 - c_EC*M_EC  )/M_EMC/c_EC   )
    kai = -3.37115 + 12.5608*c_e_kg - 7.89593*c_e_kg**2 + 3.51922*c_e_kg**3-1.15471*c_e_kg**4 +18.1863*x_EC - 6.22756*c_e_kg*x_EC - 13.6916*c_e_kg**2*x_EC +8.43904*c_e_kg**3*x_EC - 7.83732*x_EC**2 + 19.607*c_e_kg*x_EC**2  - 18.4529*c_e_kg**2*x_EC**2 -30.6369*x_EC**3 + 29.2*c_e_kg*x_EC**3 - 0.0429918*T + 0.180877*c_e_kg*T -0.0836202*c_e_kg**2*T + 0.0230098*c_e_kg**3*T + 0.195946*T*x_EC +0.0676686*c_e_kg*x_EC*T - 0.14134*c_e_kg**2*x_EC*T + 0.147429*x_EC**2*T  +0.173059*c_e_kg*x_EC**2*T - 0.51634*x_EC**3*T - 0.000223097*T**2 +0.000111233*c_e_kg*T**2 + 0.0000495286*c_e_kg**2*T**2  +0.000952777*x_EC*T**2 + 0.00117334 *c_e_kg*x_EC*T**2-0.000619157*x_EC**2*T**2 - 3.46897E-7*T**3 - 2.75041E-6*c_e_kg*T**3 -5.57653E-6*x_EC*T**3 
    if kai < 0:
        kai = 0
    return kai / 10 

def diff_3E_10(c_e, c_EC , T):
    D_Li = (
        (c_EC >= 0 ) * 3e-10
        +  (c_EC < 0 ) * 3e-10 
    )
    return D_Li
def cond_0p5(c_e, c_EC , T):
    cond = (
        (c_EC >= 0 ) * 0.5
        +  (c_EC < 0 ) * 0.5
    )
    return cond

def t_0plus_0p24(c_e, c_EC , T):
    t_0plus = (
        (c_EC >= 0 ) * 0.24
        +  (c_EC < 0 ) * 0.24 
    )
    return t_0plus
def Fun_t_0plus_Wang2021(c_e,c_EC, T):
    c_T = Fun_c_T(c_e, c_EC,T)
    y_e = c_e / c_T
    return 0.4107-1.487*y_e+2.547*y_e**2

#### Use measured LJP
# add Ruihe Li update 230201
def dLJP_dcEC_Nyman_2011(c_e, c_EC , T):
    dLJP_dcEC = -5.394e-6 - 3.616e-2 / c_EC
    return dLJP_dcEC
def dLJP_dce_Nyman_2011(c_e, c_EC , T):
    dLJP_dce = 5.326e-5 + 2.47e-2 / c_e
    return dLJP_dce

# add Ruihe Li update 240125
def dLJP_2_Species_dc_EC_np(c_e,c_EC,T):
    return -7.89754965808084e-9*c_e*c_EC*(1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)**2*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.0256925791272929*c_EC*(-3.53810224682022e-7*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.68911565857653e-8*c_EC/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 6.81/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0256925791272929*c_EC*(7.89754965808084e-9*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 3.94877482904042e-9*c_EC/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 1/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 1.01454209750984e-10*c_EC*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + c_EC*(-1.00338484700878e-11*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 8.94014496325669e-9*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 3.06797530286975e-10*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) - 3.94877482904042e-9*c_EC*(1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2) + (-c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 1)*(1.32989943307757e-15*c_e**3/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**4 - 4.97388303528732e-12*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 - 1.88096104878324e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 2.18938184642623e-9*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**0.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 6.12377610056938e-8*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**1.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 6.43473325345614e-8*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**2.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 8.90920142928264e-7*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**3.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1.31758582203603e-6*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**4.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 1.19107242247655e-10*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2) + (7.89754965808084e-9*c_e*c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)**2*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 3.94877482904042e-9*c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2) - 1/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)))*(-1.16162410618092e-11*c_e**3/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 6.51678874139771e-8*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 4.76340417020011*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 1.10889171513396*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**0.5 - 10.3386938408227*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**1.5 + 6.51820732459422*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**2.5 + 64.462681030378*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**3.5 - 74.1487833613674*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**4.5 + 0.0301630878954419*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.184138714605309) + 0.0256925791272929*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + (1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) 

def dLJP_2_Species_dc_e_np(c_e,c_EC,T):
    return c_EC*(-3.22848499937632e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 2/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)**2*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.0256925791272929*c_EC*(-1.44636127972059e-6*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 1.09929914228764e-7*c_EC/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 89.6/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0256925791272929*c_EC*(3.22848499937632e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1.61424249968816e-8*c_EC/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 4.14740531538772e-10*c_EC*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + c_EC*(-4.1017949457966e-11*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 2.99474190980181e-7*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*(-1.61424249968816e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) - 1.61424249968816e-8*c_EC*(1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2) + (-c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 1)*(5.43657280581555e-15*c_e**3/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**4 - 5.5181747304683e-11*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 5.34428802806659e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 4.76340417020011/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 1.10889171513396*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**0.5*(-8.07121249844079e-9*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 0.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e - 10.3386938408227*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**1.5*(-2.42136374953224e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e + 6.51820732459422*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**2.5*(-4.0356062492204e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 2.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e + 64.462681030378*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**3.5*(-5.64984874890856e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 3.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e - 74.1487833613674*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**4.5*(-7.26409124859671e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 4.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e + 0.0301630878954419*(-1.61424249968816e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e) + (-c_EC*(-3.22848499937632e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 2/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)**2*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 1.61424249968816e-8*c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2))*(-1.16162410618092e-11*c_e**3/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 6.51678874139771e-8*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 4.76340417020011*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 1.10889171513396*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**0.5 - 10.3386938408227*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**1.5 + 6.51820732459422*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**2.5 + 64.462681030378*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**3.5 - 74.1487833613674*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**4.5 + 0.0301630878954419*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.184138714605309) 

def dLJP_1_Specie_dc_e_np(c_e,c_EC,T):
    return -9.01770046675714e-12*c_e**2/(0.000125961489893485*c_e + 1)**3 + 7.37994000756151e-8*c_e/(0.000125961489893485*c_e + 1)**2 - 0.209856986311729/(1.50767557182561*c_e + 11969.3373990775) + 0.0383076354787938*(1.50767557182561*c_e + 11969.3373990775)*(-1.05236811106346e-8*c_e/(0.000125961489893485*c_e + 1)**2 + 1/(1.50767557182561*c_e + 11969.3373990775))/c_e 


import numpy as np

def nmc_LGM50_entropic_change_ORegan2022(sto, c_s_max):
    """
    LG M50 NMC 811 entropic change in open circuit potential (OCP) at a temperature of
    298.15K as a function of the stochiometry. The fit is taken from [1].

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
       Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
       Entropic change [V.K-1]
    """
    a1 = 0.04006
    a2 = -0.06656
    b1 = 0.2828
    b2 = 0.8032
    c1 = 0.0009855
    c2 = 0.02179

    dUdT = (
        a1 * pybamm.exp(-((sto - b1) ** 2) / c1)
        + a2 * pybamm.exp(-((sto - b2) ** 2) / c2)
    ) / 1000
    # fit in mV / K
    return dUdT

def graphite_LGM50_entropic_change_ORegan2022(sto, c_s_max):
    """
    LG M50 Graphite entropic change in open circuit potential (OCP) at a temperature of
    298.15K as a function of the stochiometry. The fit is taken from [1].

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
       Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
       Entropic change [V.K-1]
    """

    a0 = -0.1112
    a1 = -0.09002 * 0  # fixed fit (see discussion O'Regan et al 2021)
    a2 = 0.3561
    b1 = 0.4955
    b2 = 0.08309
    c0 = 0.02914
    c1 = 0.1122
    c2 = 0.004616
    d1 = 63.9

    dUdT = (
        a0 * sto
        + c0
        + a2 * pybamm.exp(-((sto - b2) ** 2) / c2)
        + a1
        * (pybamm.tanh(d1 * (sto - (b1 - c1))) - pybamm.tanh(d1 * (sto - (b1 + c1))))
    ) / 1000  # fit in mV / K

    return dUdT

def separator_LGM50_heat_capacity_ORegan2022(T):
    """
    Wet separator specific heat capacity as a function of the temperature from [1].

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Specific heat capacity
    """

    # value for the dry porous separator (i.e. separator + air, and we neglect the air
    # contribution to density)
    cp_dry = 1.494e-3 * T**3 - 1.444 * T**2 + 475.5 * T - 5.13e4
    rho_dry = 946
    theta_dry = rho_dry * cp_dry

    # value for the bulk electrolyte
    rho_e = 1280
    cp_e = 229
    eps_e = pybamm.Parameter("Separator porosity")
    theta_e = rho_e * cp_e

    # value for the wet separator
    theta_wet = theta_dry + theta_e * eps_e
    rho_wet = rho_dry + rho_e * eps_e
    cp_wet = theta_wet / rho_wet

    return cp_wet
def nmc_LGM50_heat_capacity_ORegan2022(T):

    # value for the dry porous electrode (i.e. electrode + air, and we neglect the air
    # contribution to density)
    cp_dry = -8.414e-4 * T**3 + 0.7892 * T**2 - 241.3 * T + 2.508e4
    rho_dry = 3270
    theta_dry = rho_dry * cp_dry

    # value for the bulk electrolyte
    rho_e = 1280
    cp_e = 229
    eps_e = pybamm.Parameter("Positive electrode porosity")
    theta_e = rho_e * cp_e

    # value for the wet separator
    theta_wet = theta_dry + theta_e * eps_e
    rho_wet = rho_dry + rho_e * eps_e
    cp_wet = theta_wet / rho_wet

    return cp_wet
def graphite_LGM50_heat_capacity_ORegan2022(T):
    # value for the dry porous electrode (i.e. electrode + air, and we neglect the air
    # contribution to density)
    cp_dry = 4.932e-4 * T**3 - 0.491 * T**2 + 169.4 * T - 1.897e4
    rho_dry = 1740
    theta_dry = rho_dry * cp_dry

    # value for the bulk electrolyte
    rho_e = 1280
    cp_e = 229
    eps_e = pybamm.Parameter("Negative electrode porosity")
    theta_e = rho_e * cp_e

    # value for the wet separator
    theta_wet = theta_dry + theta_e * eps_e
    rho_wet = rho_dry + rho_e * eps_e
    cp_wet = theta_wet / rho_wet

    return cp_wet
def aluminium_heat_capacity_CRC(T):
    """
    Aluminium specific heat capacity as a function of the temperature from [1].

    References
    ----------
    .. [1] William M. Haynes (Ed.). "CRC handbook of chemistry and physics". CRC Press
    (2014).

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Specific heat capacity
    """

    cp = 4.503e-6 * T**3 - 6.256e-3 * T**2 + 3.281 * T + 355.7

    return cp


def copper_thermal_conductivity_CRC(T):
    """
    Copper thermal conductivity as a function of the temperature from [1].

    References
    ----------
    .. [1] William M. Haynes (Ed.). "CRC handbook of chemistry and physics". CRC Press
    (2014).

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Thermal conductivity
    """

    lambda_th = -5.409e-7 * T**3 + 7.054e-4 * T**2 - 0.3727 * T + 463.6

    return lambda_th

def copper_thermal_conductivity_CRC(T):
    """
    Copper thermal conductivity as a function of the temperature from [1].

    References
    ----------
    .. [1] William M. Haynes (Ed.). "CRC handbook of chemistry and physics". CRC Press
    (2014).

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Thermal conductivity
    """

    lambda_th = -5.409e-7 * T**3 + 7.054e-4 * T**2 - 0.3727 * T + 463.6

    return lambda_th
def graphite_LGM50_thermal_conductivity_ORegan2022(T):
    """
    Wet negative electrode thermal conductivity as a function of the temperature from
    [1].

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Thermal conductivity
    """

    lambda_wet = -2.61e-4 * T**2 + 0.1726 * T - 24.49

    return lambda_wet
def nmc_LGM50_thermal_conductivity_ORegan2022(T):
    """
    Wet positive electrode thermal conductivity as a function of the temperature from
    [1].

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Thermal conductivity
    """

    lambda_wet = 2.063e-5 * T**2 - 0.01127 * T + 2.331

    return lambda_wet

def copper_heat_capacity_CRC(T):
    """
    Copper specific heat capacity as a function of the temperature from [1].

    References
    ----------
    .. [1] William M. Haynes (Ed.). "CRC handbook of chemistry and physics". CRC Press
    (2014).

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Specific heat capacity
    """

    cp = 1.445e-6 * T**3 - 1.946e-3 * T**2 + 0.9633 * T + 236

    return cp
def graphite_LGM50_diffusivity_ORegan2022(sto, T):
    """
    LG M50 Graphite diffusivity as a function of stochiometry, in this case the
    diffusivity is taken to be a constant. The value is taken from [1].

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
       Electrode stochiometry
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Solid diffusivity
    """

    a0 = 11.17
    a1 = -1.553
    a2 = -6.136
    a3 = -9.725
    a4 = 1.85
    b1 = 0.2031
    b2 = 0.5375
    b3 = 0.9144
    b4 = 0.5953
    c0 = -15.11
    c1 = 0.0006091
    c2 = 0.06438
    c3 = 0.0578
    c4 = 0.001356
    d = 2092

    D_ref = (
        10
        ** (
            a0 * sto
            + c0
            + a1 * pybamm.exp(-((sto - b1) ** 2) / c1)
            + a2 * pybamm.exp(-((sto - b2) ** 2) / c2)
            + a3 * pybamm.exp(-((sto - b3) ** 2) / c3)
            + a4 * pybamm.exp(-((sto - b4) ** 2) / c4)
        )
        * 3.0321  # correcting factor (see O'Regan et al 2021)
    )

    E_D_s = d * pybamm.constants.R
    arrhenius = pybamm.exp(E_D_s / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius
def nmc_LGM50_electronic_conductivity_ORegan2022(T):
    """
    Positive electrode electronic conductivity as a function of the temperature from
    [1].

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Thermal conductivity
    """

    E_r = 3.5e3
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    sigma = 0.8473 * arrhenius

    return sigma
def nmc_LGM50_diffusivity_ORegan2022(sto, T):
    """
    NMC diffusivity as a function of stoichiometry, in this case the
    diffusivity is taken to be a constant. The value is taken from [1].

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
       Electrode stochiometry
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Solid diffusivity
    """

    a1 = -0.9231
    a2 = -0.4066
    a3 = -0.993
    b1 = 0.3216
    b2 = 0.4532
    b3 = 0.8098
    c0 = -13.96
    c1 = 0.002534
    c2 = 0.003926
    c3 = 0.09924
    d = 1449

    D_ref = (
        10
        ** (
            c0
            + a1 * pybamm.exp(-((sto - b1) ** 2) / c1)
            + a2 * pybamm.exp(-((sto - b2) ** 2) / c2)
            + a3 * pybamm.exp(-((sto - b3) ** 2) / c3)
        )
        * 2.7  # correcting factor (see O'Regan et al 2021)
    )

    E_D_s = d * pybamm.constants.R
    arrhenius = pybamm.exp(E_D_s / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius

def graphite_LGM50_electrolyte_exchange_current_density_ORegan2022(
    c_e, c_s_surf, c_s_max, T
):
    """
    Exchange-current density for Butler-Volmer reactions between graphite and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """

    i_ref = 2.668  # (A/m2)
    alpha = 0.792
    E_r = 4e4
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    c_e_ref = pybamm.Parameter("Typical electrolyte concentration [mol.m-3]")

    return (
        i_ref
        * arrhenius
        * (c_e / c_e_ref) ** (1 - alpha)
        * (c_s_surf / c_s_max) ** alpha
        * (1 - c_s_surf / c_s_max) ** (1 - alpha)
    )
def nmc_LGM50_electrolyte_exchange_current_density_ORegan2022(
    c_e, c_s_surf, c_s_max, T
):
    """
    Exchange-current density for Butler-Volmer reactions between NMC and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parameters of a high energy lithium-ion
    cylindrical battery." Electrochimica Acta 425 (2022): 140700

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    i_ref = 5.028  # (A/m2)
    alpha = 0.43
    E_r = 2.401e4
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    c_e_ref = pybamm.Parameter("Typical electrolyte concentration [mol.m-3]")

    return (
        i_ref
        * arrhenius
        * (c_e / c_e_ref) ** (1 - alpha)
        * (c_s_surf / c_s_max) ** alpha
        * (1 - c_s_surf / c_s_max) ** (1 - alpha)
    )

def nmc_LGM50_lithiation_ocp_OKane2023(sto):
    """
    LG M50 NMC lithiation open-circuit potential as a function of stoichiometry.
    Fitted to unpublished measurements by Kieran O'Regan.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry
    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential
    """

    U = (
        -0.7983 * sto
        + 4.513
        - 0.03269 * pybamm.tanh(19.83 * (sto - 0.5424))
        - 18.23 * pybamm.tanh(14.33 * (sto - 0.2771))
        + 18.05 * pybamm.tanh(14.46 * (sto - 0.2776))
    )

    return U

def graphite_LGM50_delithiation_ocp_OKane2023(sto):
    """
    LG M50 Graphite delithiation open-circuit potential as a function of stochiometry.
    Fitted to unpublished measurements taken by Kieran O'Regan.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open circuit potential
    """

    u_eq = (
        1.051 * pybamm.exp(-26.76 * sto)
        + 0.1916
        - 0.05598 * pybamm.tanh(35.62 * (sto - 0.1356))
        - 0.04483 * pybamm.tanh(14.64 * (sto - 0.2861))
        - 0.02097 * pybamm.tanh(26.28 * (sto - 0.6183))
        - 0.02398 * pybamm.tanh(38.1 * (sto - 1))
    )

    return u_eq


# Mark Ruihe block end

# Call dict via a function to avoid errors when editing in place
def get_parameter_values():
    """
    Parameters for an LG M50 cell, from the paper

    > Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. ["Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery
    Models."](https://iopscience.iop.org/article/10.1149/1945-7111/ab9050) Journal of
    the Electrochemical Society 167 (2020): 080534

    and references therein.

    SEI parameters are example parameters for SEI growth from the papers:

    > Ramadass, P., Haran, B., Gomadam, P. M., White, R., & Popov, B. N. (2004).
    Development of first principles capacity fade model for Li-ion cells. Journal of the
     Electrochemical Society, 151(2), A196-A203.
    > Ploehn, H. J., Ramadass, P., & White, R. E. (2004). Solvent diffusion model for
    aging of lithium-ion battery cells. Journal of The Electrochemical Society, 151(3),
    A456-A462.
    > Single, F., Latz, A., & Horstmann, B. (2018). Identifying the mechanism of
    continued growth of the solid–electrolyte interphase. ChemSusChem, 11(12),
    1950-1955.
    > Safari, M., Morcrette, M., Teyssot, A., & Delacour, C. (2009). Multimodal Physics-
    Based Aging Model for Life Prediction of Li-Ion Batteries. Journal of The
    Electrochemical Society, 156(3),
    > Yang, X., Leng, Y., Zhang, G., Ge, S., Wang, C. (2017). Modeling of lithium
    plating induced aging of lithium-ion batteries: Transition from linear to nonlinear
    aging. Journal of Power Sources, 360, 28-40.

    Note: this parameter set does not claim to be representative of the true parameter
    values. Instead these are parameter values that were used to fit SEI models to
    observed experimental data in the referenced papers.
    """

    return {
        # sei
        "Ratio of lithium moles to SEI moles": 2.0,
        "Inner SEI reaction proportion": 0.5,
        "Inner SEI partial molar volume [m3.mol-1]": 9.585e-05,
        "Outer SEI partial molar volume [m3.mol-1]": 9.585e-05,
        "SEI reaction exchange current density [A.m-2]": 1.5e-07,
        "SEI resistivity [Ohm.m]": 200000.0,
        "Outer SEI solvent diffusivity [m2.s-1]": 2.5000000000000002e-22,
        "Bulk solvent concentration [mol.m-3]": 2636.0,
        "Inner SEI open-circuit potential [V]": 0.1,
        "Outer SEI open-circuit potential [V]": 0.8,
        "Inner SEI electron conductivity [S.m-1]": 8.95e-14,
        "Inner SEI lithium interstitial diffusivity [m2.s-1]": 1e-19,
        "Lithium interstitial reference concentration [mol.m-3]": 15.0,
        "Initial inner SEI thickness [m]": 1.2362e-08, # Simon new: 1.2362e-08,
        "Initial outer SEI thickness [m]": 1.2362e-08, # Simon new: 1.2362e-08,
        "EC diffusivity [m2.s-1]": 2e-18,
        "SEI kinetic rate constant [m.s-1]": 1e-12,
        "SEI open-circuit potential [V]": 0.4,
        "SEI growth activation energy [J.mol-1]": 0.0, # Simon: 38000.0,
        "Negative electrode reaction-driven LAM factor [m3.mol-1]": 0.0,
        "Positive electrode reaction-driven LAM factor [m3.mol-1]": 0.0,
        # cell
        "Negative current collector thickness [m]": 1.2e-05,
        "Negative electrode thickness [m]": 8.52e-05,
        "Separator thickness [m]": 1.2e-05,
        "Positive electrode thickness [m]": 7.56e-05,
        "Positive current collector thickness [m]": 1.6e-05,
        "Electrode height [m]": 0.065,
        "Electrode width [m]": 1.58,
        "Cell cooling surface area [m2]": 0.00531,
        "Cell volume [m3]": 2.42e-05,
        "Cell thermal expansion coefficient [m.K-1]": 1.1e-06,
        "Negative current collector conductivity [S.m-1]": 58411000.0,
        "Positive current collector conductivity [S.m-1]": 36914000.0,
        "Negative current collector density [kg.m-3]": 8960.0,
        "Positive current collector density [kg.m-3]": 2700.0,
        "Negative current collector specific heat capacity [J.kg-1.K-1]"
        "": copper_heat_capacity_CRC,
        "Positive current collector specific heat capacity [J.kg-1.K-1]"
        "": aluminium_heat_capacity_CRC,
        "Negative current collector thermal conductivity [W.m-1.K-1]"
        "": copper_thermal_conductivity_CRC,
        "Positive current collector thermal conductivity [W.m-1.K-1]": 237.0,
        "Nominal cell capacity [A.h]": 5.0,
        "Typical current [A]": 5.0,
        "Current function [A]": 5.0,
        "Contact resistance [Ohm]": 0.013, # Mark Ruihe add contact resistance
        # negative electrode
        "Negative electrode conductivity [S.m-1]": 215.0,
        "Maximum concentration in negative electrode [mol.m-3]": 32544.0,# Mark Ruihe change from 33133.0,
        "Negative electrode diffusivity [m2.s-1]"
        "": graphite_LGM50_diffusivity_ORegan2022,
        # Mark Ruihe change "Negative electrode OCP [V]": graphite_LGM50_ocp_Chen2020,
        "Negative electrode OCP [V]": graphite_LGM50_delithiation_ocp_OKane2023,
        "Negative electrode porosity": 0.240507,
        "Negative electrode active material volume fraction": 0.75,
        "Negative particle radius [m]": 5.86e-06,
        "Negative electrode Bruggeman coefficient (electrolyte)": 1.5,
        "Negative electrode Bruggeman coefficient (electrode)": 0.0, # Mark Ruihe, change from 1.5  
        "Negative electrode cation signed stoichiometry": -1.0,
        "Negative electrode electrons in reaction": 1.0,
        "Negative electrode charge transfer coefficient": 0.5,
        "Negative electrode double-layer capacity [F.m-2]": 0.2,
        "Negative electrode exchange-current density [A.m-2]"
        "": graphite_LGM50_electrolyte_exchange_current_density_ORegan2022, # ask Simon change here
        "Negative electrode density [kg.m-3]": 2060.0,
        "Negative electrode specific heat capacity [J.kg-1.K-1]"
        "": graphite_LGM50_heat_capacity_ORegan2022,
        "Negative electrode thermal conductivity [W.m-1.K-1]"
        "": graphite_LGM50_thermal_conductivity_ORegan2022,
        "Negative electrode OCP entropic change [V.K-1]"
        "": graphite_LGM50_entropic_change_ORegan2022,
        # positive electrode
        "Positive electrode conductivity [S.m-1]"
        "": nmc_LGM50_electronic_conductivity_ORegan2022,
        "Maximum concentration in positive electrode [mol.m-3]": 52787.0, # Mark Ruihe change from 63104.0,
        "Positive electrode diffusivity [m2.s-1]": nmc_LGM50_diffusivity_ORegan2022,
        "Positive electrode OCP [V]": nmc_LGM50_lithiation_ocp_OKane2023, #  nmc_LGM50_ocp_Chen2020,
        "Positive electrode porosity": 0.335,
        "Positive electrode active material volume fraction": 0.665,
        "Positive particle radius [m]": 5.22e-06,
        "Positive electrode Bruggeman coefficient (electrolyte)": 1.5,
        "Positive electrode Bruggeman coefficient (electrode)": 0.0,  # Mark Ruihe, change from 1.5  
        "Positive electrode cation signed stoichiometry": -1.0,
        "Positive electrode electrons in reaction": 1.0,
        "Positive electrode charge transfer coefficient": 0.5,
        "Positive electrode double-layer capacity [F.m-2]": 0.2,
        "Positive electrode exchange-current density [A.m-2]"
        "": nmc_LGM50_electrolyte_exchange_current_density_ORegan2022, # ask Simon
        "Positive electrode density [kg.m-3]": 3699.0,
        "Positive electrode specific heat capacity [J.kg-1.K-1]"
        "": nmc_LGM50_heat_capacity_ORegan2022,
        "Positive electrode thermal conductivity [W.m-1.K-1]"
        "": nmc_LGM50_thermal_conductivity_ORegan2022,
        "Positive electrode OCP entropic change [V.K-1]"
        "": nmc_LGM50_entropic_change_ORegan2022,
        # separator
        "Separator porosity": 0.47,
        "Separator Bruggeman coefficient (electrolyte)": 1.5,
        "Separator density [kg.m-3]": 1548.0,
        "Separator specific heat capacity [J.kg-1.K-1]"
        "": separator_LGM50_heat_capacity_ORegan2022,
        "Separator thermal conductivity [W.m-1.K-1]": 0.3344,
        # electrolyte
        "Typical electrolyte concentration [mol.m-3]": 1000.0,
        "Initial concentration in electrolyte [mol.m-3]": 1000.0,
        "Cation transference number": Fun_t_0plus_Wang2021,   # from Andrew 
        "1 + dlnf/dlnc": 1.0,
        "TDF of EC": 5.0, 
        "Measured dLJP_dcEC": dLJP_2_Species_dc_EC_np,
        "Measured dLJP_dce": dLJP_2_Species_dc_e_np,
        "Electrolyte diffusivity [m2.s-1]": diff_3E_10,
        "Electrolyte conductivity [S.m-1]": electrolyte_conductivity_Wang2021,

        # Mark Ruihe block start
        # Typical electrolyte: 1M LiPF6 in EMC:EC=1:1 wt%, 
        #           c_e=1000, c_EMC=5252.58, c_EC=6209.49, c_T=13462.07
        "EC transference number": Fun_Xi_tidle,# Update 240313 - becomes a function
        "EC transference number zero": 0.7  , # from Andrew": 
        "EC initial concentration in electrolyte [mol.m-3]": 6209.49,
        "Typical EC concentration [mol.m-3]": 6209.49, 
        #"Background solvent concentration [mol.m-3]": Fun_c_EMC,  # should from Andrew, add temperoaliy
        "Typical total concentration [mol.m-3]":13462.07,
        "Total concentration [mol.m-3]":Fun_c_T,
        # Update 240125: Distinguish EC,e cross diffusivity and e,EC cross diffusivity
        "EC Lithium ion cross diffusivity [m2.s-1]": Dimensional_EC_Lithium_ion_cross_diffusivity,  # D_(EC,e)^0
        "Lithium ion EC cross diffusivity [m2.s-1]": 1.5e-11,  # Free parameter - D_(e,EC)^0
        "Typical EC Lithium ion cross diffusivity [m2.s-1]": 1.5e-12,
        "EC diffusivity in electrolyte [m2.s-1]": EC_diffusivity_5E_10,     #from Andrew
        "Typical EC diffusivity in electrolyte [m2.s-1]": 5E-10, #from Andrew
        "EC diffusivity in SEI [m2.s-1]": 2E-18               ,  #from Yang2017": 
        "Typical EC diffusivity in SEI [m2.s-1]": 2E-18,

        "EC partial molar volume [m3.mol-1]": 6.5312e-5     ,     # Mark Ruihe Li add
        "Li partial molar volume [m3.mol-1]": 1.3e-5   ,         #Mark Ruihe Li add
        "CH2OCO2Li2 partial molar volume [m3.mol-1]": 9.585e-5 , # Mark Ruihe Li add
        # Mark Ruihe block end

        # experiment
        "Reference temperature [K]": 298.15,
        "Total heat transfer coefficient [W.m-2.K-1]": 10.0,
        "Ambient temperature [K]": 298.15,
        "Number of electrodes connected in parallel to make a cell": 1.0,
        "Number of cells connected in series to make a battery": 1.0,
        "Lower voltage cut-off [V]": 2.5,
        "Upper voltage cut-off [V]": 4.4,
        "Initial concentration in negative electrode [mol.m-3]": 28543.0,
        "Initial concentration in positive electrode [mol.m-3]": 12727.0,
        "Initial temperature [K]": 298.15,
        # citations
        "citations": ["ORegan2022","Chen2020","OKane2022", "OKane2020",],
    }

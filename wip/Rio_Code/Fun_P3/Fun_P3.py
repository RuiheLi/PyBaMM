import pybamm;import pandas as pd   ;import numpy as np;import os;
import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;
import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'
from textwrap import wrap
import openpyxl
import traceback
import random;import time, signal

from pybamm import tanh,exp,sqrt


def EC_diffusivity_7E_10(c_e, c_EC , T):
    D_ec_dim = (
        (c_EC >= 0 ) * 7e-10
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim

def EC_diffusivity_1E_10(c_e, c_EC , T):
    D_ec_dim = (
        (c_EC >= 0 ) * 1e-10
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim
def EC_diffusivity_5E_11(c_e, c_EC , T):
    D_ec_dim = (
        (c_EC >= 0 ) * 5e-11
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim

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


def Diff_Andrew_ACS(c_e,c_EC, T):
    M_EMC = 104.105/1e3 # kg/mol
    M_e   = 151.905/1e3 # kg/mol
    c_e_constant = 3800 # mol/m3
    # Get rho using Table S1:
    c_e_used =(
        (c_e <= c_e_constant) * c_e/ 1e3
        +
        (c_e > c_e_constant)  * c_e_constant/ 1e3
    )
    rho = (
        1007.1+114.2*c_e_used 
        - 8.121*np.power(c_e_used,1.5) 
        - 4.013e-5*np.power(c_e_used,10) ) # kg/m3
    # get y using Eq. (S16)
    y = M_EMC*c_e_used*1e3  / (
        rho + (2*M_EMC-M_e) * c_e_used*1e3 
    )
    #print("diffusivity:",y)
    Diff = (4.998-29.96*y+53.78*np.power(y,2) ) / 1e10
    return Diff
def Density_Andrew_ACS(c_e,c_EC, T):
    M_EMC = 104.105/1e3 # kg/mol
    M_e   = 151.905/1e3 # kg/mol
    c_e_constant = 3800 # mol/m3
    # Get rho using Table S1:
    c_e_used =(
        (c_e <= c_e_constant) * c_e / 1e3
        +
        (c_e > c_e_constant)  * c_e_constant / 1e3
    )
    rho = (
        1007.1+114.2*c_e_used 
        - 8.121*np.power(c_e_used,1.5) 
        - 4.013e-5*np.power(c_e_used,10) ) # kg/m3
    return rho
def Cond_Andrew_ACS(c_e,c_EC, T):
    M_EMC = 104.105/1e3 # kg/mol
    M_e   = 151.905/1e3 # kg/mol
    c_e_constant = 3800 # mol/m3
    # Get rho using Table S1:
    c_e_used =(
        (c_e <= c_e_constant) * c_e / 1e3
        +
        (c_e > c_e_constant)  * c_e_constant / 1e3
    )
    rho = (
        1007.1+114.2*c_e_used 
        - 8.121*np.power(c_e_used,1.5) 
        - 4.013e-5*np.power(c_e_used,10) ) # kg/m3
    # get y using Eq. (S16)
    y = M_EMC*c_e_used *1e3 / (
        rho + (2*M_EMC-M_e) * c_e_used * 1e3
    )
    #print("conductivity:",y)
    cond = np.power( (
        48.93*np.power(y,1.5)
        - 284.8* np.power(y,2.5 )
        + 817.7* np.power(y,4 ) ), 2 )
    return cond
def t_0plus_Andrew_ACS(c_e,c_EC, T):
    M_EMC = 104.105/1e3 # kg/mol
    M_e   = 151.905/1e3 # kg/mol
    c_e_constant = 3800 # mol/m3
    # Get rho using Table S1:
    c_e_used =(
        (c_e <= c_e_constant) * c_e/ 1e3
        +
        (c_e > c_e_constant)  * c_e_constant/ 1e3
    )
    rho = (
        1007.1+114.2*c_e_used 
        - 8.121*np.power(c_e_used,1.5) 
        - 4.013e-5*np.power(c_e_used,10) ) # kg/m3
    # get y using Eq. (S16)
    y = M_EMC*c_e_used*1e3  / (
        rho + (2*M_EMC-M_e) * c_e_used*1e3 
    )
    #print("transference:",y)
    t_0plus = 0.4107 - 1.487*y + 2.547* np.power(y,2)
    return t_0plus

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

def EC_Lithium_ion_cross_diffusivity_1e_11(c_e, c_EC , T):
    D_0_EC_e = (
        (c_EC >= 0 ) * (1e-11 )
        +  (c_EC < 0 ) * 0 
    )
    return D_0_EC_e


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


def EC_transference_number(c_e,c_EC, T):# Mark Ruihe add update 221212
    c_EC_0 = pybamm.Parameter("Typical EC concentration [mol.m-3]")
    Xi_0 =   pybamm.Parameter("EC transference number zero")
    
    Xi = ( 
        (c_EC < 0 ) * 0 
        + 
        (c_EC >= 0) * (Xi_0 * c_EC / c_EC_0 )
        #+
        #(c_EC > c_EC_0 ) * Xi_0
    )
    return Xi

def EC_diffusivity_5E_5(c_e, c_EC , T):
    D_ec_dim = (
        (c_EC >= 0 ) * 5e-5 
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim
def EC_diffusivity_5E_10(c_e, c_EC , T):
    D_ec_dim = (
        (c_EC >= 0 ) * 5e-10
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim

def EC_diffusivity_1E_10(c_e, c_EC , T):
    D_ec_dim = (
        (c_EC >= 0 ) * 1e-10
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim

def EC_diffusivity_3E_10(c_e, c_EC , T):
    D_ec_dim = (
        (c_EC >= 0 ) * 3e-10
        +  (c_EC < 0 ) * 0 
    )
    return D_ec_dim
def Cross_diffusivity_1p5E_12(c_e, c_EC , T):
    D_x_dim_1 = (
        (c_EC >= 0 ) * 1.5e-12 
        +  (c_EC < 0 ) * 0 
    )
    D_x_dim = (
        (c_e >= 0 ) * D_x_dim_1 
        +  (c_e < 0 ) * 0 
    )
    return D_x_dim

def Cross_diffusivity_1p5E_11(c_e, c_EC , T):
    D_x_dim_1 = (
        (c_EC >= 0 ) * 1.5e-11
        +  (c_EC < 0 ) * 0 
    )
    D_x_dim = (
        (c_e >= 0 ) * D_x_dim_1 
        +  (c_e < 0 ) * 0 
    )
    return D_x_dim

def Cross_diffusivity_1p5E_10(c_e, c_EC , T):
    D_x_dim_1 = (
        (c_EC >= 0 ) * 1.5e-10
        +  (c_EC < 0 ) * 0 
    )
    D_x_dim = (
        (c_e >= 0 ) * D_x_dim_1 
        +  (c_e < 0 ) * 0 
    )
    return D_x_dim

def Cross_diffusivity_1p5E_09(c_e, c_EC , T):
    D_x_dim_1 = (
        (c_EC >= 0 ) * 1.5e-9
        +  (c_EC < 0 ) * 0 
    )
    D_x_dim = (
        (c_e >= 0 ) * D_x_dim_1 
        +  (c_e < 0 ) * 0 
    )
    return D_x_dim

def t_0plus_constant(c_e, c_EC , T):
    t_0plus = (
        (c_EC >= 0 ) * 0.24
        +  (c_EC < 0 ) * 0.24 
    )
    return t_0plus

def elely_TDF_15(c_e, c_EC , T):
    Oneplus_dlnfdlnc = (
        (c_EC >= 0 ) * 15
        +  (c_EC < 0 ) * 15
    )
    return Oneplus_dlnfdlnc

def nmc_LGM50_diffusivity_ORegan2022(sto, T):
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

def graphite_LGM50_diffusivity_ORegan2022(sto, T):
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

import numpy as np
# Landesfeind electrolyte properties 
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

def electrolyte_diffusivity_base_Landesfeind2019(c_e,c_EC, T, coeffs):
    c = (
        (c_e >= 0 ) * c_e / 1000 
        +  (c_e < 0 ) * 0  
        )# mol.m-3 -> mol.l
    p1, p2, p3, p4 = coeffs
    A = p1 * pybamm.exp(p2 * c)
    B = pybamm.exp(p3 / T)
    C = pybamm.exp(p4 * c / T)
    D_e = A * B * C * 1e-10  # m2/s
    return D_e

def electrolyte_TDF_base_Landesfeind2019(c_e,c_EC, T, coeffs):
    c = c_e / 1000  # mol.m-3 -> mol.l
    p1, p2, p3, p4, p5, p6, p7, p8, p9 = coeffs
    tdf = (
        p1
        + p2 * c
        + p3 * T
        + p4 * c**2
        + p5 * c * T
        + p6 * T**2
        + p7 * c**3
        + p8 * c**2 * T
        + p9 * c * T**2
    )
    return tdf

def electrolyte_transference_number_base_Landesfeind2019(c_e,c_EC, T, coeffs):
    c = c_e / 1000  # mol.m-3 -> mol.l
    p1, p2, p3, p4, p5, p6, p7, p8, p9 = coeffs
    tplus = (
        p1
        + p2 * c
        + p3 * T
        + p4 * c**2
        + p5 * c * T
        + p6 * T**2
        + p7 * c**3
        + p8 * c**2 * T
        + p9 * c * T**2
    )
    return tplus

def electrolyte_transference_number_EC_EMC_3_7_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array(
        [
            -1.28e1,
            -6.12,
            8.21e-2,
            9.04e-1,
            3.18e-2,
            -1.27e-4,
            1.75e-2,
            -3.12e-3,
            -3.96e-5,
        ]
    )
    return electrolyte_transference_number_base_Landesfeind2019(c_e,c_EC, T, coeffs)
def electrolyte_transference_number_EC_EMC_3_7_Landesfeind2019_Con(c_e,c_EC, T):
    coeffs = np.array(
        [
            -1.28e1,
            -6.12,
            8.21e-2,
            9.04e-1,
            3.18e-2,
            -1.27e-4,
            1.75e-2,
            -3.12e-3,
            -3.96e-5,
        ]
    )
    c_e_constant = 4000
    t0plus = (
        (c_e <= c_e_constant) * electrolyte_transference_number_base_Landesfeind2019(c_e,c_EC, T, coeffs)
        + 
        (c_e > c_e_constant) * electrolyte_transference_number_base_Landesfeind2019(c_e_constant,c_EC, T, coeffs)
    )
    return t0plus

def electrolyte_TDF_EC_EMC_3_7_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array(
        [2.57e1, -4.51e1, -1.77e-1, 1.94, 2.95e-1, 3.08e-4, 2.59e-1, -9.46e-3, -4.54e-4]
    )
    return electrolyte_TDF_base_Landesfeind2019(c_e,c_EC, T, coeffs)


def electrolyte_diffusivity_EC_EMC_3_7_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array([1.01e3, 1.01, -1.56e3, -4.87e2])
    return electrolyte_diffusivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)
def electrolyte_diffusivity_EC_EMC_3_7_Landesfeind2019_Con(c_e,c_EC, T):
    coeffs = np.array([1.01e3, 1.01, -1.56e3, -4.87e2])
    c_e_constant = 4000
    diff_f = (
        (c_e <= c_e_constant) * electrolyte_diffusivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)
        +
        (c_e > c_e_constant) * electrolyte_diffusivity_base_Landesfeind2019(c_e_constant,c_EC, T, coeffs)
    )
    return diff_f


def electrolyte_conductivity_EC_EMC_3_7_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array([5.21e-1, 2.28e2, -1.06, 3.53e-1, -3.59e-3, 1.48e-3])
    return electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)
def electrolyte_conductivity_EC_EMC_3_7_Landesfeind2019_Con(c_e,c_EC, T):
    coeffs = np.array([5.21e-1, 2.28e2, -1.06, 3.53e-1, -3.59e-3, 1.48e-3])
    c_e_constant = 4000
    sigma = (
        (c_e <= c_e_constant) * electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)
        +
        (c_e > c_e_constant) * electrolyte_conductivity_base_Landesfeind2019(c_e_constant,c_EC, T, coeffs)
    )
    return sigma

def electrolyte_diffusivity_EC_DMC_1_1_Landesfeind2019_Con(c_e, c_EC,T):
    coeffs = np.array([1.47e3, 1.33, -1.69e3, -5.63e2])
    c_e_constant = 4000
    diff_f = (
        (c_e <= c_e_constant) * electrolyte_diffusivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)
        +
        (c_e > c_e_constant) * electrolyte_diffusivity_base_Landesfeind2019(c_e_constant,c_EC, T, coeffs)
    )
    return diff_f

def electrolyte_diffusivity_EC_DMC_1_1_Landesfeind2019(c_e, c_EC,T):
    coeffs = np.array([1.47e3, 1.33, -1.69e3, -5.63e2])
    return electrolyte_diffusivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)

def electrolyte_conductivity_EC_DMC_1_1_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array([7.98e-1, 2.28e2, -1.22, 5.09e-1, -4e-3, 3.79e-3])
    return electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)

def electrolyte_conductivity_EC_DMC_1_1_Landesfeind2019_Con(c_e,c_EC, T):
    coeffs = np.array([7.98e-1, 2.28e2, -1.22, 5.09e-1, -4e-3, 3.79e-3])
    c_e_constant = 4000
    sigma = (
        (c_e <= c_e_constant) * electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)
        +
        (c_e >  c_e_constant) * electrolyte_conductivity_base_Landesfeind2019(c_e_constant,c_EC, T, coeffs)
    )
    return sigma

def electrolyte_TDF_EC_DMC_1_1_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array(
        [-5.58, 7.17, 3.80e-2, 1.91, -6.65e-2, -5.08e-5, 1.1e-1, -6.10e-3, 1.51e-4]
    )
    return electrolyte_TDF_base_Landesfeind2019(c_e,c_EC, T, coeffs)

def electrolyte_transference_number_EC_DMC_1_1_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array(
        [
            -7.91,
            2.45e-1,
            5.28e-2,
            6.98e-1,
            -1.08e-2,
            -8.21e-5,
            7.43e-4,
            -2.22e-3,
            3.07e-5,
        ]
    )
    return electrolyte_transference_number_base_Landesfeind2019(c_e,c_EC, T, coeffs)
def electrolyte_transference_number_EC_DMC_1_1_Landesfeind2019_Con(c_e,c_EC, T):
    coeffs = np.array(
        [
            -7.91,
            2.45e-1,
            5.28e-2,
            6.98e-1,
            -1.08e-2,
            -8.21e-5,
            7.43e-4,
            -2.22e-3,
            3.07e-5,
        ]
    )
    c_e_constant = 4000
    t0plus = (
        (c_e <= c_e_constant) * electrolyte_transference_number_base_Landesfeind2019(c_e,c_EC, T, coeffs)
        + 
        (c_e > c_e_constant) * electrolyte_transference_number_base_Landesfeind2019(c_e_constant,c_EC, T, coeffs)
    )
    return t0plus

def electrolyte_diffusivity_EMC_FEC_19_1_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array([5.86e2, 1.33, -1.38e3, -5.82e2])
    return electrolyte_diffusivity_base_Landesfeind2019(c_e, c_EC,T, coeffs)

def electrolyte_TDF_EMC_FEC_19_1_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array(
        [3.22, -1.01e1, -1.58e-2, 6.12, 2.96e-2, 2.42e-5, -2.22e-1, -1.57e-2, 6.30e-6]
    )
    return electrolyte_TDF_base_Landesfeind2019(c_e,c_EC, T, coeffs)

def electrolyte_transference_number_EMC_FEC_19_1_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array(
        [-1.22e1, -3.05, 8.38e-2, 1.78, 1.51e-3, -1.37e-4, -2.45e-2, -5.15e-3, 2.14e-5]
    )
    return electrolyte_transference_number_base_Landesfeind2019(c_e, c_EC,T, coeffs)

def electrolyte_conductivity_EMC_FEC_19_1_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array([2.51e-2, 1.75e2, 1.23, 2.05e-1, -8.81e-2, 2.83e-3])
    return electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)

def electrolyte_conductivity_Valoen2005Constant(c_e,c_EC, T):# Mark Ruihe change
    # T = T + 273.15
    # mol/m3 to molar
    c_e = c_e / 1000
    c_e_constant = 4500/1000
    sigma = (c_e <= c_e_constant ) * (
        (1e-3 / 1e-2) * (
        c_e
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )) + (c_e > c_e_constant ) *  (
        (1e-3 / 1e-2) * (
        c_e_constant
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e_constant * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e_constant ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    ))
    # mS/cm to S/m
    return sigma
def electrolyte_diffusivity_Valoen2005Constant(c_e,c_EC, T): 
    # T = T + 273.15
    # mol/m3 to molar
    c_e = c_e / 1000
    c_e_constant = 4500/1000

    T_g = 229 + 5 * c_e
    T_g_constant = 229 + 5 * c_e_constant
    D_0 = -4.43 - 54 / (T - T_g)
    D_0_constant = -4.43 - 54 / (T - T_g_constant)
    D_1 = -0.22

    D_final = (c_e <= c_e_constant) * (
        (10 ** (D_0 + D_1 * c_e)) * 1e-4
    ) + (c_e > c_e_constant) *  (
        (10 ** (D_0_constant + D_1 * c_e_constant)) * 1e-4
    )

    # cm2/s to m2/s
    # note, in the Valoen paper, ln means log10, so its inverse is 10^x
    return D_final
def electrolyte_conductivity_Nyman2008Exp(c_e,c_EC, T):
    sigma_e = (
        0.1 * 0.06248 * (1+298.15-0.05559) * 
        (c_e/1e3) * (1 - 3.084 *sqrt(c_e/1e3) 
        + 1.33 *(1+ 0.03633 *(exp(1000/298.15))*c_e/1e3)   ) 
        / (1+(c_e/1e3)**4*( 0.00795 *exp(1000/298.15))) 
    )
    return sigma_e
def electrolyte_diffusivity_Nyman2008Exp(c_e,c_EC, T):
    D_c_e = (
        6 * exp( -1 *(c_e/1000)) 
        * exp(-5/298.15) 
        * exp(-95/298.15*(c_e/1000)) * 1e-10 
    )
    return D_c_e

def electrolyte_conductivity_Andrew2022_poly4(x1,y1, T):# x:Li+,y:ec
    p00 =      0.9608 # (0.9188, 1.003)
    p10 =      0.1502 # (0.09883, 0.2015)
    p01 =       0.173 # (0.1481, 0.1979)
    p20 =     -0.3934 # (-0.4666, -0.3202)
    p11 =     -0.1179 # (-0.1627, -0.07308)
    p02 =     -0.1472 # (-0.1767, -0.1176)
    p30 =     0.04244 # (0.01066, 0.07422)
    p21 =     -0.1197 # (-0.1423, -0.09704)
    p12 =   -0.003226 # (-0.02331, 0.01686)
    p40 =     0.01664 # (-0.01452, 0.0478)
    p31 =     0.05145 # (0.02312, 0.07978)
    p22 =     0.05983 # (0.0346, 0.08506)
    x1 = (x1<2300)*x1 + (x1>=2300)*2300
    x= (x1-971.4 )/733.1
    y= (y1-4737) / 3439
    kai  = (
        (x1 > 0 )  * (
        p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y + p30*x*x*x + p21*x*x*y 
        + p12*x*y*y + p40*x*x*x*x + p31*x*x*x*y + p22*x*x*y*y
        )
        + (x1 <= 0 ) * 0 )
    kai_final = (kai>0) * kai + (kai<0) * 0
    return kai_final


def diff_constant(c_e, c_EC , T):
    D_Li = (
        (c_EC >= 0 ) * 3e-10
        +  (c_EC < 0 ) * 3e-10 
    )
    return D_Li
def cond_constant(c_e, c_EC , T):
    cond = (
        (c_EC >= 0 ) * 0.7
        +  (c_EC < 0 ) * 0.7
    )
    return cond

def electrolyte_conductivity_Andrew2022(x,y, T):# x:Li+,y:ec
    p00 =     -0.2524;
    p10 =    0.001402;
    p01 =   0.0001142 ;
    p20 =  -5.379e-07  ;
    p11 =  -1.399e-08 ;
    p02 =  -8.137e-09  ;
    kai  = (
        (x > 0 )  *(
        (p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y)
        + (x <= 0 ) * 0 )
    )
    kai_final = (kai>0) * kai + (kai<0) * 0
    return kai_final

def electrolyte_diffusivity_Valoen2005Constant(c_e,c_EC, T): 
    # T = T + 273.15
    # mol/m3 to molar
    c_e = c_e / 1000
    c_e_constant = 4500/1000

    T_g = 229 + 5 * c_e
    T_g_constant = 229 + 5 * c_e_constant
    D_0 = -4.43 - 54 / (T - T_g)
    D_0_constant = -4.43 - 54 / (T - T_g_constant)
    D_1 = -0.22

    D_final = (c_e <= c_e_constant) * (
        (10 ** (D_0 + D_1 * c_e)) * 1e-4
    ) + (c_e > c_e_constant) *  (
        (10 ** (D_0_constant + D_1 * c_e_constant)) * 1e-4
    )

    # cm2/s to m2/s
    # note, in the Valoen paper, ln means log10, so its inverse is 10^x
    return D_final

def electrolyte_conductivity_Valoen2005Constant(c_e,c_EC, T):# Mark Ruihe change
    # T = T + 273.15
    # mol/m3 to molar
    c_e = c_e / 1000
    c_e_constant = 4500/1000
    sigma = (c_e <= c_e_constant ) * (c_e > 0 ) * (
        (1e-3 / 1e-2) * (
        c_e
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )) + (c_e > c_e_constant ) *  (
        (1e-3 / 1e-2) * (
        c_e_constant
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e_constant * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e_constant ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )) + (c_e <= 0 ) * 0 
    # mS/cm to S/m
    return sigma
    
def electrolyte_conductivity_Valoen2005Constant_wEC_Haya(c_e,c_EC, T):
    # mol/m3 to molar
    c_e = c_e / 1000
    # T = T + 273.15
    sigma = (c_e <= 4.5) * (
        (1e-3 / 1e-2) * (
        c_e
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )) + (c_e > 4.5) *  (
        (1e-3 / 1e-2) * (
        4.5
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + 4.5 * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + 4.5 ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    ))
    a=1.092; b=-6.497e-6; c=-0.7877; d=-0.0004808
    ratio= (
        a*exp(b*c_EC)+c*exp(d*c_EC) )
    return sigma*ratio
def electrolyte_conductivity_Valoen2005Constant_ECtanh500_1(c_e,c_EC, T):# Mark Ruihe change
    # mol/m3 to molar
    c_e = c_e / 1000
    # T = T + 273.15
    sigma = (c_e <= 4.5) * (
        (1e-3 / 1e-2) * (
        c_e
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + c_e * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + c_e ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )) + (c_e > 4.5) *  (
        (1e-3 / 1e-2) * (
        4.5
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T ** 2)
            + 4.5 * (0.668 - 0.0178 * T + 2.80e-5 * T ** 2)
            + 4.5 ** 2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    ))
    coff = 1
    ratio = ( (1-coff)+ coff/2 + coff/2 *  tanh((c_EC-4500*0.5)/500))
    return sigma*ratio

def dLJP_2_Species_dc_EC_np(c_e,c_EC,T):
    return -7.89754965808084e-9*c_e*c_EC*(1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)**2*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.0256925791272929*c_EC*(-3.53810224682022e-7*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.68911565857653e-8*c_EC/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 6.81/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0256925791272929*c_EC*(7.89754965808084e-9*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 3.94877482904042e-9*c_EC/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 1/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 1.01454209750984e-10*c_EC*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + c_EC*(-1.00338484700878e-11*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 8.94014496325669e-9*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 3.06797530286975e-10*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) - 3.94877482904042e-9*c_EC*(1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2) + (-c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 1)*(1.32989943307757e-15*c_e**3/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**4 - 4.97388303528732e-12*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 - 1.88096104878324e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 2.18938184642623e-9*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**0.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 6.12377610056938e-8*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**1.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 6.43473325345614e-8*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**2.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 8.90920142928264e-7*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**3.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1.31758582203603e-6*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**4.5*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 1.19107242247655e-10*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2) + (7.89754965808084e-9*c_e*c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)**2*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 3.94877482904042e-9*c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2) - 1/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)))*(-1.16162410618092e-11*c_e**3/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 6.51678874139771e-8*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 4.76340417020011*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 1.10889171513396*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**0.5 - 10.3386938408227*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**1.5 + 6.51820732459422*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**2.5 + 64.462681030378*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**3.5 - 74.1487833613674*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**4.5 + 0.0301630878954419*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.184138714605309) + 0.0256925791272929*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + (1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) 

def dLJP_2_Species_dc_e_np(c_e,c_EC,T):
    return c_EC*(-3.22848499937632e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 2/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)**2*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.0256925791272929*c_EC*(-1.44636127972059e-6*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 1.09929914228764e-7*c_EC/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 89.6/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0256925791272929*c_EC*(3.22848499937632e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1.61424249968816e-8*c_EC/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 4.14740531538772e-10*c_EC*(-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(89.6*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 6.81*c_EC/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 12.6)/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + c_EC*(-4.1017949457966e-11*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 2.99474190980181e-7*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*(-1.61424249968816e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) - 1.61424249968816e-8*c_EC*(1.31463627670492e-7*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 - 2.26403007269705*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 0.0776943592809339*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.211527003955003)/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2) + (-c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 1)*(5.43657280581555e-15*c_e**3/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**4 - 5.5181747304683e-11*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 5.34428802806659e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 4.76340417020011/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 1.10889171513396*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**0.5*(-8.07121249844079e-9*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 0.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e - 10.3386938408227*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**1.5*(-2.42136374953224e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e + 6.51820732459422*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**2.5*(-4.0356062492204e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 2.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e + 64.462681030378*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**3.5*(-5.64984874890856e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 3.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e - 74.1487833613674*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**4.5*(-7.26409124859671e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 4.5/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e + 0.0301630878954419*(-1.61424249968816e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 1/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)/c_e) + (-c_EC*(-3.22848499937632e-8*c_e/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 2/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)**2*(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 1.61424249968816e-8*c_EC/((-2*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) + 1)*(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2))*(-1.16162410618092e-11*c_e**3/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**3 + 6.51678874139771e-8*c_e**2/(0.000156004935299578*c_e + 3.81620705585473e-5*c_EC + 1)**2 + 4.76340417020011*c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266) - 1.10889171513396*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**0.5 - 10.3386938408227*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**1.5 + 6.51820732459422*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**2.5 + 64.462681030378*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**3.5 - 74.1487833613674*(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266))**4.5 + 0.0301630878954419*pybamm.log(c_e/(1.50767557182561*c_e + 0.368808983131977*c_EC + 9664.28125450266)) + 0.184138714605309) 

def dLJP_1_Specie_dc_e_np(c_e,c_EC,T):
    return -9.01770046675714e-12*c_e**2/(0.000125961489893485*c_e + 1)**3 + 7.37994000756151e-8*c_e/(0.000125961489893485*c_e + 1)**2 - 0.209856986311729/(1.50767557182561*c_e + 11969.3373990775) + 0.0383076354787938*(1.50767557182561*c_e + 11969.3373990775)*(-1.05236811106346e-8*c_e/(0.000125961489893485*c_e + 1)**2 + 1/(1.50767557182561*c_e + 11969.3373990775))/c_e 


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
    kai_final = (kai>0) * kai  + (kai<=0) * 0

    return kai_final / 10


def gr_i_ex_2p68(c_e, c_s_surf, c_s_max, T):
    i_ref = 2.668  # (A/m2)
    alpha = 0.792; E_r = 4e4
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))
    c_e_ref = pybamm.Parameter("Typical electrolyte concentration [mol.m-3]")
    return (
        i_ref
        * arrhenius
        * (c_e / c_e_ref) ** (1 - alpha)
        * (c_s_surf / c_s_max) ** alpha
        * (1 - c_s_surf / c_s_max) ** (1 - alpha)
    )
def gr_i_ex_1p68(c_e, c_s_surf, c_s_max, T):
    i_ref = 1.668  # (A/m2)
    alpha = 0.792; E_r = 4e4
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))
    c_e_ref = pybamm.Parameter("Typical electrolyte concentration [mol.m-3]")
    return (
        i_ref
        * arrhenius
        * (c_e / c_e_ref) ** (1 - alpha)
        * (c_s_surf / c_s_max) ** alpha
        * (1 - c_s_surf / c_s_max) ** (1 - alpha)
    )
def gr_i_ex_3p68(c_e, c_s_surf, c_s_max, T):
    i_ref = 3.668  # (A/m2)
    alpha = 0.792; E_r = 4e4
    arrhenius = pybamm.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))
    c_e_ref = pybamm.Parameter("Typical electrolyte concentration [mol.m-3]")
    return (
        i_ref
        * arrhenius
        * (c_e / c_e_ref) ** (1 - alpha)
        * (c_s_surf / c_s_max) ** alpha
        * (1 - c_s_surf / c_s_max) ** (1 - alpha)
    )

def nc_i_ex_5(c_e, c_s_surf, c_s_max, T):
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
def nc_i_ex_3(c_e, c_s_surf, c_s_max, T):
    i_ref = 3.028  # (A/m2)
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
def nc_i_ex_1(c_e, c_s_surf, c_s_max, T):
    i_ref = 1.028  # (A/m2)
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


# DEFINE my callback:
class RioCallback(pybamm.callbacks.Callback):
    def __init__(self, logfile=None):
        self.logfile = logfile
        self.success  = True
        if logfile is None:
            # Use pybamm's logger, which prints to command line
            self.logger = pybamm.logger
        else:
            # Use a custom logger, this will have its own level so set it to the same
            # level as the pybamm logger (users can override this)
            self.logger = pybamm.get_new_logger(__name__, logfile)
            self.logger.setLevel(pybamm.logger.level)
    
    def on_experiment_error(self, logs):
        self.success  = False
    def on_experiment_infeasible(self, logs):
        self.success  = False

def write_excel_xlsx(path, sheet_name, value):
    import numpy as np
    index = len(value)
    workbook = openpyxl.Workbook()  # 新建工作簿（默认有一个sheet？）
    sheet = workbook.active  # 获得当前活跃的工作页，默认为第一个工作页
    sheet.title = sheet_name  # 给sheet页的title赋值
    for i in range(0, index):
        for j in range(0, len(value[i])):
            sheet.cell(row=i + 1, column=j + 1, value=str(value[i][j]))  # 行，列，值 这里是从1开始计数的
    workbook.save(path)  # 一定要保存
    print("Successfully create a excel file")

# define to kill too long runs - Jianbo Huang kindly writes this
class TimeoutError(Exception):
	pass
def handle_signal(signal_num, frame):
	raise TimeoutError

# valid only for single parameter scan:
def PlotDynamics(Sol,str,Para_scan,BasicPath , Target,Save,fs):
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    # plot overall figures
    label = [str+f"={Para_i}" for Para_i in Para_scan]
    output_variables1 = [
        "Terminal voltage [V]",   
        "Discharge capacity [A.h]",
        "EC concentration [mol.m-3]",
        "Electrolyte concentration [mol.m-3]",
        "Li+ flux [mol.m-2.s-1]",
        "EC flux [mol.m-2.s-1]",
    ]
    quick_plot = pybamm.QuickPlot(
        [sol for sol in Sol], 
        output_variables1,label,
        variable_limits='fixed',time_unit='hours',n_rows=2,
        figsize = (12,8)) #     spatial_unit='mm',
    quick_plot.dynamic_plot();
    if Save == 'True':
        quick_plot.create_gif(
            number_of_images=10, duration=2,
            output_filename=BasicPath + Target+"All cases - overview.gif")

    output_variables2 = [
        "Electrolyte concentration",
        "Minus div Li+ flux",
        "Li+ source term",
        "Minus div Li+ flux by diffusion",
        "Minus div Li+ flux by migration",
        "Minus div Li+ flux by solvent",
        #"Li+ source term refill",
    ]
    quick_plot = pybamm.QuickPlot(
        [sol for sol in Sol], 
        output_variables2,label,
        variable_limits='tight',time_unit='hours',n_rows=2,
        figsize = (12,9)) #     spatial_unit='mm',
    quick_plot.dynamic_plot();
    if Save == 'True':
        quick_plot.create_gif(
            number_of_images=10, duration=2,
            output_filename=BasicPath + Target+"All cases - c(Li) and breakdown.gif")
    else:
        pass

    output_variables3 = [
        "EC concentration",
        "EC source term (SEI)",
        "Minus div EC flux",
        "Minus div EC flux by diffusion",
        "Minus div EC flux by migration",
        "Minus div EC flux by Li+",
        #"Li+ source term refill",
    ]
    quick_plot = pybamm.QuickPlot(
        [sol for sol in Sol], 
        output_variables3,label,
        variable_limits='fixed',time_unit='hours',n_rows=2,
        figsize = (12,9)) #     spatial_unit='mm',
    quick_plot.dynamic_plot();
    if Save == 'True':
        quick_plot.create_gif(
            number_of_images=10, duration=2,
            output_filename=BasicPath + Target+"All cases - c(EC) and breakdown.gif")
    else:
        pass
    return 


# Plot a pair of loc dependent varibles - different cycles
def Plot_Loc_Var( key_all, my_dict,colormap ,fs): # for my_dict only
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    Num_subplot = len(key_all); # must have 2+ keys
    fig, axs = plt.subplots(1,Num_subplot, figsize=(7*Num_subplot,5),tight_layout=True)
    for i in range(0,Num_subplot):
        cmap_i = mpl.cm.get_cmap(colormap, len(my_dict[ key_all[i]] ) ) 
        if 'Negative' in key_all[i] or 'negative' in key_all[i]:
            x_loc = "x_n [m]";
        elif 'Positive' in key_all[i] or 'positive' in key_all[i]:
            x_loc = "x_p [m]";
        elif 'Seperator' in key_all[i] or 'seperator' in key_all[i]:
            x_loc = "x_s [m]";
        else:
            x_loc = "x [m]";
        for j in range(0,len(my_dict[ key_all[i] ])):
            axs[i].plot(
                my_dict[x_loc], 
                my_dict[ key_all[i] ][j],'-',
                color=cmap_i(j),)
            axs[i].set_title(key_all[i] ,   fontdict={'family':'DejaVu Sans','size':fs-1})
            #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
            axs[i].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
            
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].ticklabel_format( axis='x', style='sci',scilimits=[-0.01,0.01], useOffset=None, useLocale=None, useMathText=None)
            #axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)  
    return fig, axs 

# Plot a pair of loc dependent varibles - within one step
from textwrap import wrap
def Plot_Loc_Var_sol( sol,x_loc_all, key_all, cycle, step,colormap,fs): # for initial solution object
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    Num_subplot = len(key_all); # must have 2+ keys
    fig, axs = plt.subplots(1,Num_subplot, figsize=(4*Num_subplot,5),tight_layout=True)
    for i in range(0,Num_subplot):
        x_loc=x_loc_all[i]; key=key_all[i];
        LinesNmum = len(sol.cycles[cycle].steps[step][key].entries[0,:] )
        cmap_i = mpl.cm.get_cmap(colormap, LinesNmum) 
        for j in range(0,LinesNmum):
            axs[i].plot(
                sol.cycles[cycle].steps[step][x_loc].entries[:,0], 
                sol.cycles[cycle].steps[step][key].entries[:,j], '-',
                color=cmap_i(j),)
            axs[i].set_title(
                "\n".join(wrap(key, 30)) # to be adjusted
                ,   fontdict={'family':'DejaVu Sans','size':fs-1})
            #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
            axs[i].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
            
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].ticklabel_format( axis='x', style='sci',scilimits=[-0.01,0.01], useOffset=None, useLocale=None, useMathText=None)
            #axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)  
    return  fig, axs 

def Plot_Loc_Var_sol_6( sol,x_loc_all, key_all, cycle, step,colormap,fs): # for initial solution object
    font = {'family' : 'DejaVu Sans','size'   : fs}
    mpl.rc('font', **font)
    Num_subplot = len(key_all); # must have 2+ keys
    if Num_subplot <= 3:
        fig, axs = plt.subplots(1,Num_subplot, figsize=(2.5*Num_subplot,5),tight_layout=True)
        for i in range(0,Num_subplot):
            x_loc=x_loc_all[i]; key=key_all[i];
            LinesNmum = len(sol.cycles[cycle].steps[step][key].entries[0,:] )
            print(f"Total line number is: {temp_L_Num}")
            if temp_L_Num > 26:
                    xx = np.arange(0,temp_L_Num,int(np.rint(temp_L_Num/50)));
            else:
                xx = np.arange(0,temp_L_Num,1);
            cmap_i = mpl.cm.get_cmap(colormap, LinesNmum) 
            for j in range(0,LinesNmum):
                axs[i].plot(
                    sol.cycles[cycle].steps[step][x_loc].entries[:,0], 
                    sol.cycles[cycle].steps[step][key].entries[:,j], '-',
                    color=cmap_i(j),)
                axs[i].set_title(key ,   fontdict={'family':'DejaVu Sans','size':fs-1})
                #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
                axs[i].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
                
                labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
                
                axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
                axs[i].ticklabel_format( axis='x', style='sci',
                scilimits=[-0.01,0.01], useOffset=None, 
                useLocale=None, useMathText=None)
                #axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)  
    elif Num_subplot >3:
        fig, axs = plt.subplots(2,3, figsize=(2.5*Num_subplot,8.5),tight_layout=True)
        Plot_Count = 0
        for m in range(0,2):
            for n in range(0,3):
                x_loc=x_loc_all[Plot_Count]; 
                key=key_all[Plot_Count]; Plot_Count +=1;
                temp_L_Num= len(sol.cycles[cycle].steps[step][key].entries[0,:] )
                print(f"Total line number is: {temp_L_Num}")
                if temp_L_Num > 26:
                    xx = np.arange(0,temp_L_Num,int(np.rint(temp_L_Num/50)));
                else:
                    xx = np.arange(0,temp_L_Num,1);
                xx = xx.tolist()
                if not xx[-1]==temp_L_Num-1:
                    xx.append(temp_L_Num-1)
                cmap_i = mpl.cm.get_cmap(colormap, len(xx)) 
                for index_j,j in zip(  range(0,len(xx)), xx):
                    len_1 = len(sol.cycles[cycle].steps[step][x_loc].entries[:,0])
                    len_2 = len(sol.cycles[cycle].steps[step][key].entries[:,j])
                    if len_1>len_2:
                        axs[m,n].plot(
                            sol.cycles[cycle].steps[step][x_loc].entries[0:len_2,0], 
                            sol.cycles[cycle].steps[step][key].entries[:,j], '-',
                            color=cmap_i(index_j),)
                    elif len_1<len_2:
                        axs[m,n].plot(
                            sol.cycles[cycle].steps[step][x_loc].entries[:,0], 
                            sol.cycles[cycle].steps[step][key].entries[0:len_1,j], '-',
                            color=cmap_i(index_j),)
                    else:    
                        axs[m,n].plot(
                            sol.cycles[cycle].steps[step][x_loc].entries[:,0], 
                            sol.cycles[cycle].steps[step][key].entries[:,j], '-',
                            color=cmap_i(index_j),)
                    axs[m,n].set_title(key ,   fontdict={'family':'DejaVu Sans','size':fs-1})
                    #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
                    axs[m,n].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
                    
                    labels = axs[m,n].get_xticklabels() + axs[m,n].get_yticklabels(); 
                    [label.set_fontname('DejaVu Sans') for label in labels]
                    
                    axs[m,n].tick_params(
                        labelcolor='k', labelsize=fs, width=1) ;  del labels;
                    axs[m,n].ticklabel_format( 
                        axis='x', style='sci',
                        scilimits=[-0.01,0.01], useOffset=None, 
                        useLocale=None, useMathText=None)
    return  fig, axs 

# Fig. 1
def Plot_Fig_1(Full_cycle,my_dict_AGE,
    BasicPath, Target,   Scan_i,  fs,  dpi):
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    fig, axs = plt.subplots(3,3, figsize=(15,12),tight_layout=True)
    axs[0,0].plot(
        Full_cycle, 
        my_dict_AGE["Discharge capacity [A.h]"] ,
        '-o',label=f"Scan = {Scan_i}")
    axs[0,0].set_ylabel("Capacity [A.h]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[0,0].set_title("Discharge capacity",   fontdict={'family':'DejaVu Sans','size':fs+1})

    axs[0,1].plot(Full_cycle, my_dict_AGE["CDend Loss of capacity to SEI [A.h]"] ,'-o',)
    axs[0,1].set_ylabel("Capacity [A.h]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[0,1].set_title("Loss of capacity to SEI",   fontdict={'family':'DejaVu Sans','size':fs+1})

    axs[0,2].plot(Full_cycle, my_dict_AGE["CDend Local ECM resistance [Ohm]"] ,'-o',)
    axs[0,2].set_ylabel("Resistance [Ohm]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[0,2].set_title("Local ECM resistance",   fontdict={'family':'DejaVu Sans','size':fs+1})

    axs[1,0].plot(Full_cycle, my_dict_AGE["CDcyc Positive SOC range"] ,'-o',label="Pos" )
    axs[1,0].plot(Full_cycle, my_dict_AGE["CDcyc Negative SOC range"] ,'-^',label="Neg" )
    axs[1,1].plot(
        Full_cycle, 
        my_dict_AGE["CC LLI to SEI per Ah"],
        '--^',label="CC" )
    axs[1,1].plot(
        Full_cycle, 
        my_dict_AGE["CV LLI to SEI per Ah"],
        '--^',label="CV" )
    axs[1,1].plot(
        Full_cycle, 
        my_dict_AGE["Cha LLI to SEI per Ah"],
        '-s',label="Charge" )
    axs[1,1].plot(
        Full_cycle, 
        my_dict_AGE["Dis LLI to SEI per Ah"],
        '-o',label="Discharge" )
    axs[1,2].plot(
        Full_cycle, 
        my_dict_AGE["CDend Loss of active material in positive electrode [%]"],
        '-o',label="Pos" )
    axs[1,2].plot(
        Full_cycle, 
        my_dict_AGE["CDend Loss of active material in negative electrode [%]"],
        '-^',label="Neg" )
    axs[1,0].set_ylabel("SOC",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1,0].set_title("SOC range (Dis)",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1,1].set_ylabel("Capacity [A.h]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1,1].set_title("LLI to SEI per Ah",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1,1].set_title("LLI to SEI per Ah",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1,2].set_ylabel("LAM %",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1,2].set_title("LAM",   fontdict={'family':'DejaVu Sans','size':fs+1})

    axs[2,0].plot(
        Full_cycle, 
        my_dict_AGE["CDend Total EC in electrolyte [mol]"] ,
        '-o',label="In elely" )
    axs[2,0].plot(
        Full_cycle, 
        my_dict_AGE["CDend Total EC in electrolyte and SEI [mol]"] ,
        '-^',label="In elely and SEI"  )
    axs[2,1].plot(
        Full_cycle, 
        my_dict_AGE["CDend Total lithium in electrolyte [mol]"],
        '-o', )
    axs[2,2].plot(
        Full_cycle, 
        my_dict_AGE["CDend Total lithium in particles [mol]"],
        '-o',)
    axs[2,0].set_ylabel("Quantity [mol]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[2,0].set_title("Total EC",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[2,1].set_ylabel("Quantity [mol]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[2,1].set_title("Total lithium in electrolyte",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[2,2].set_ylabel("Quantity [mol]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[2,2].set_title("Total lithium in particles",   fontdict={'family':'DejaVu Sans','size':fs+1})


    for i in range(0,3):
        for j in range(0,3):
            axs[i,j].set_xlabel("Cycle numbers",   fontdict={'family':'DejaVu Sans','size':fs})
            
            labels = axs[i,j].get_xticklabels() + axs[i,j].get_yticklabels(); 
            [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i,j].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i,j].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)     
    axs[1,1].ticklabel_format( 
        axis='y', style='sci',scilimits=[-0.01,0.01], 
        useOffset=None, useLocale=None, useMathText=None)
    plt.savefig(
        BasicPath + Target+f"{Scan_i}th Scan/" + 
        "Fig. 1 - Cycle based overall change.png", dpi=dpi)

        
def Plot_Loc_Var_2( Full_cycle, key_all, my_dict,fs): # for my_dict only
    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    Num_subplot = len(key_all); # must have 2+ keys
    fig, axs = plt.subplots(3,3, figsize=(15,12),tight_layout=True)
    count_i = 0; 
    for i in range(0,3):
        for j in range(0,3):
            key_ij = key_all[count_i]; count_i += 1
            if 'Negative' in key_ij or 'negative' in key_ij:
                x_loc = "x_n [m]";
            elif 'Positive' in key_ij or 'positive' in key_ij:
                x_loc = "x_p [m]";
            elif 'Seperator' in key_ij or 'seperator' in key_ij:
                x_loc = "x_s [m]";
            else:
                x_loc = "x [m]";
            X_Len = min(len(my_dict[x_loc]),len(my_dict[ key_ij ][0]))
            #print(x_loc,X_Len)
            axs[i,j].plot(my_dict[x_loc][0:X_Len], my_dict[ key_ij ][0][0:X_Len],'-o',label="1st cycle")
            axs[i,j].plot(my_dict[x_loc][0:X_Len], my_dict[ key_ij ][-1][0:X_Len],'-^',label=f"{Full_cycle[-1]}th cycle")
            axs[i,j].set_title(
                "\n".join(wrap(key_ij, 30)) # to be adjusted
                ,   fontdict={'family':'DejaVu Sans','size':fs-3})
            #axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
            axs[i,j].set_xlabel(x_loc,   fontdict={'family':'DejaVu Sans','size':fs})
            
            labels = axs[i,j].get_xticklabels() + axs[i,j].get_yticklabels(); 
            [label.set_fontname('DejaVu Sans') for label in labels]
            
            axs[i,j].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i,j].ticklabel_format( 
                axis='x', style='sci',scilimits=[-0.01,0.01], 
                useOffset=None, useLocale=None, useMathText=None)
            axs[i,j].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)  
    return fig, axs 


def Plot_Last_Single_Step(
    sol,cycle,step,BasicPath, Target,Scan_i,
    index_cyc,Save,colormap,fs,dpi):

    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    # 1st plot Li+ and flux
    i =Scan_i;
    try:
        fig, axs= Plot_Loc_Var_sol_6(
            sol,
            ["x [m]","x [m]","x [m]","x [m]","x [m]","x [m]",], 
            ["Electrolyte concentration",
            "Minus div Li+ flux",
            "Li+ source term",
            "Minus div Li+ flux by diffusion",
            "Minus div Li+ flux by migration",
            "Minus div Li+ flux by solvent",], 
            cycle,step,colormap,fs)
    except Exception as e:
        traceback.print_exc()
    else:
        pass

    if Save == 'True':
        plt.savefig(
            BasicPath + 
            Target +     f"{i}th Scan/" +
            f"after {index_cyc}th cycle-cycle={cycle}, step={step}, c(Li+) and flux.png", dpi=dpi)
    else:
        pass
    
    # 2nd plot: for EC and flux
    try:
        fig, axs = Plot_Loc_Var_sol_6(
            sol,
            ["x [m]","x [m]","x [m]","x [m]","x [m]","x [m]",], 
            ["EC concentration",
            "c(EC) over c(Li+)",
            "Minus div EC flux",
            "Minus div EC flux by diffusion",
            "Minus div EC flux by migration",
            "Minus div EC flux by Li+",], 
            cycle,step,colormap,fs)
    except Exception as e:
        traceback.print_exc()
    else:
        pass
    if Save == 'True':
        plt.savefig(
            BasicPath + 
            Target+ f"{i}th Scan/" +
            f"after {index_cyc}th cycle-cycle={cycle}, step={step}, c(EC) and flux.png", dpi=dpi)
    else:
        pass
    # 3rd plot: porosity, potential
    try:
        fig, axs = Plot_Loc_Var_sol_6(
            sol,
            ["x_n [m]","x_p [m]","x [m]","x [m]","x [m]","x [m]",], 
            [
                "Negative electrode porosity",
                "Positive electrode potential [V]",
                "Electrolyte current density [A.m-2]",
                "Electrolyte potential [V]",
                "Electrolyte diffusivity [m2.s-1]",
                "Electrolyte conductivity [S.m-1]",
            ], 
            cycle,step,colormap,fs)
    except Exception as e:
        traceback.print_exc()
    else:
        pass
    if Save == 'True':
        plt.savefig(
            BasicPath + 
            Target+f"{i}th Scan/" +
            f"after {index_cyc}th cycle-cycle={cycle}, step={step}, Porosity, elely and potential.png", dpi=dpi)
    else:
        pass

    return


def Plot_Single_Static(Sol,str,cycle, step, 
    Para_scan,BasicPath , Target,Save,colormap,fs,dpi):

    font = {'family' : 'DejaVu Sans','size': fs}
    mpl.rc('font', **font)
    for sol, Para_i in zip(Sol,Para_scan):
        # 1st plot
        fig, axs= Plot_Loc_Var_sol(
            sol,
            ["x [m]","x [m]","x [m]",], 
            ["Electrolyte concentration",
            "Minus div Li+ flux",
            "Li+ source term",], 
            cycle,step,colormap,fs)
        if Save == 'True':
            plt.savefig(
                BasicPath + 
                Target +    str +
                f"={Para_i} - cycle={cycle}, step={step}, c(Li+) and flux.png", dpi=dpi)
        else:
            pass
        
        # 2nd plot
        fig, axs= Plot_Loc_Var_sol(
            sol,
            ["x [m]","x [m]","x [m]",], 
            ["Minus div Li+ flux by diffusion",
            "Minus div Li+ flux by migration",
            "Minus div Li+ flux by solvent",], 
            cycle,step,colormap,fs)
        if Save == 'True':
            plt.savefig(
                BasicPath + 
                Target+str+
                f"={Para_i} - cycle={cycle}, step={step}, Li+ flux breakdown.png", dpi=dpi)
        else:
            pass
        

        fig, axs = Plot_Loc_Var_sol(
            sol,
            ["x [m]","x [m]","x [m]",], 
            ["EC concentration",
            "c(EC) over c(Li+)",
            "Minus div EC flux",], 
            cycle,step,colormap,fs)
        if Save == 'True':
            plt.savefig(
                BasicPath + 
                Target+str+
                f"={Para_i} - cycle={cycle}, step={step}, c(EC) and flux.png", dpi=dpi)
        else:
            pass

        fig, axs = Plot_Loc_Var_sol(
            sol,
            ["x [m]","x [m]","x [m]",], 
            ["Minus div EC flux by diffusion",
            "Minus div EC flux by migration",
            "Minus div EC flux by Li+",], 
            cycle,step,colormap,fs)
        if Save == 'True':
            plt.savefig(
                BasicPath + 
                Target+str+
                f"={Para_i} - cycle={cycle}, step={step}, EC flux breakdown.png", dpi=dpi)
        else:
            pass
    return



'''
有一定概率出错的函数
import random
def may_cause_error():
	if (random.random() > 0.8):
	    return 1 / 0
	else:
	    return 0
'''
def recursive_scan(mylist,kvs, key_list, acc):
    # 递归终止条件
    if len(key_list) == 0:
        mylist.append(acc.copy())   # copy.deepcopy(acc) 如果value是个list，就要deep copy了
        return mylist
    # 继续递归
    k = key_list[0]
    for v in kvs[k]:
        acc[k] = v
        # print(v)
        recursive_scan(mylist,kvs, key_list[1:], acc)


def Run_P3_OneCycle(Rate_Dis,Rate_Cha,model,para,str_model,str_para,var_pts):
    para_used = para.copy()
    if (
        model.options["solvent diffusion"]=="double spatial consume w refill" 
        and 
        model.options["electrolyte conductivity"]=='full'):
        para_used.update({"EC diffusivity in electrolyte [m2.s-1]":EC_diffusivity_5E_5})
        # print('Using EC_diffusivity_5E_5')
    V_max = 4.2;        V_min = 2.5
    if Rate_Dis > 4:
        ts_dis = 0.5
    else: 
        ts_dis = 2
    if Rate_Cha > 4:
        ts_cha = 1
    else:
        ts_cha = 5
    Exp_1  = pybamm.Experiment(
    [ (
        f"Hold at {V_max} V until C/100",
        f"Discharge at {Rate_Dis} C until {V_min} V ({ts_dis} second period)", 
        f"Charge at {Rate_Cha} C until {V_max} V ({ts_cha} second period)", 
        f"Hold at {V_max} V until C/100")    ] * 1 )  

    c_e = model.variables["Electrolyte concentration [mol.m-3]"]
    c_EC= model.variables["EC concentration [mol.m-3]"]
    T = model.variables["Cell temperature [K]"]
    D_e = para_used["Electrolyte diffusivity [m2.s-1]"]
    D_EC= para_used["EC diffusivity in electrolyte [m2.s-1]"]
    sigma_e = para_used["Electrolyte conductivity [S.m-1]"]
    Xi = para_used["EC transference number"]
    model.variables["Electrolyte diffusivity [m2.s-1]"] = D_e(c_e,c_EC, T)
    model.variables["EC diffusivity in electrolyte [m2.s-1]"] = D_EC(c_e,c_EC, T)
    model.variables["Electrolyte conductivity [S.m-1]"] = sigma_e(c_e,c_EC, T)
    model.variables["EC transference number"] = Xi(c_e,c_EC, T)
    model.variables["c(EC) over c(Li+)"] = c_EC / c_e
    t_0plus = para_used["Cation transference number"]
    model.variables["Cation transference number"] = t_0plus(c_e,c_EC, T)
    
    sim    = pybamm.Simulation(
        model, experiment = Exp_1,
        parameter_values = para_used,
        solver = pybamm.CasadiSolver(return_solution_if_failed_early=True),
        var_pts=var_pts,
        )       
    sol    = sim.solve()
    MyDict = {}
    try:
        MyDict ['Discharge Capacity'] = abs(
            sol.cycles[0].steps[1]['Discharge capacity [A.h]'].entries[0] - 
            sol.cycles[0].steps[1]['Discharge capacity [A.h]'].entries[-1] )
        MyDict ['T rise'] = abs(
            sol.cycles[0].steps[1]["Volume-averaged cell temperature [C]"].entries[-1] - 
            sol.cycles[0].steps[1]["Volume-averaged cell temperature [C]"].entries[0] )
    except:
        MyDict ['Discharge Capacity'] = 0
    else:
        pass
    try:
        MyDict ['Charge Capacity'] = abs(
            sol.cycles[0].steps[2]['Discharge capacity [A.h]'].entries[0] - 
            sol.cycles[0].steps[3]['Discharge capacity [A.h]'].entries[-1] )#
    except:
        MyDict ['Charge Capacity']    = 0
    else:
        pass
    MyDict['Solution']  = sol
    MyDict['Rate_Dis']  = Rate_Dis
    MyDict['Rate_Cha']  = Rate_Cha
    print("Finish one cycle")
    return MyDict
    
def Scan_Crate(Rate_Dis_All,Rate_Cha_All,model,para,str_model,str_para,var_pts):   
    Case_Dict = {}
    MyDict_All =[]; Cap_Dis_All = []  ; Cap_Cha_All = []; Trise_All = [];
    if str_model == "Model_DFN": # one species only
        para.update({"Measured dLJP_dce":dLJP_One_Specie_dce_Jung2023})
    else: # use default, which is two speices
        pass
        #para.update({"Measured dLJP_dce":dLJP_One_Specie_dce_Jung2023})
    # scan C-rate
    for Rate_Dis in Rate_Dis_All:
        for Rate_Cha in Rate_Cha_All:
            MyDict_All.append(
                Run_P3_OneCycle(
                    Rate_Dis,Rate_Cha,model,para,
                    str_model,str_para,var_pts))
            Cap_Dis_All.append(MyDict_All[-1]['Discharge Capacity'])
            Cap_Cha_All.append(MyDict_All[-1]['Charge Capacity'])
            Trise_All.append(MyDict_All[-1]['T rise'])
    Case_Dict ['Cap_Dis_All'] = Cap_Dis_All
    Case_Dict ['Cap_Cha_All'] = Cap_Cha_All
    Case_Dict ['str_model'] = str_model
    Case_Dict ['str_para']  = str_para
    Case_Dict ['MyDict_All']  = MyDict_All
    Case_Dict ['Rate_Dis_All']  = Rate_Dis_All
    Case_Dict ['Rate_Cha_All']  = Rate_Cha_All
    Case_Dict ['Trise_All']  = Trise_All
    return Case_Dict

def Para_init_Dict(Para_dict):
    Para_dict_used = Para_dict.copy();
    Para_0=pybamm.ParameterValues(Para_dict_used["Para_Set"]  )
    Para_dict_used.pop("Para_Set")

    if Para_dict_used.__contains__("Mesh list"):
        Mesh_list = Para_dict_used["Mesh list"]  
        Para_dict_used.pop("Mesh list")
    if Para_dict_used.__contains__("Model option"):
        model_options = Para_dict_used["Model option"]  
        Para_dict_used.pop("Model option")
    # about initial SOC
    if Para_dict_used.__contains__("Initial Neg SOC"):
        c_Neg1SOC_in = (
            Para_0["Maximum concentration in negative electrode [mol.m-3]"]
            *Para_dict_used["Initial Neg SOC"]  )
        Para_0.update(
            {"Initial concentration in negative electrode [mol.m-3]":
            c_Neg1SOC_in})
        Para_dict_used.pop("Initial Neg SOC")
    if Para_dict_used.__contains__("Initial Pos SOC"):    
        c_Pos1SOC_in = (
            Para_0["Maximum concentration in positive electrode [mol.m-3]"]
            *Para_dict_used["Initial Pos SOC"]  )
        Para_0.update(
            {"Initial concentration in positive electrode [mol.m-3]":
            c_Pos1SOC_in})
        Para_dict_used.pop("Initial Pos SOC")

    CyclePack = [ Mesh_list,model_options];
    
    for key, value in Para_dict_used.items():
        # risk: will update parameter that doesn't exist, 
        # so need to make sure the name is right 
        if isinstance(value, str):
            Para_0.update({key: eval(value)})
            #Para_dict_used.pop(key)
        else:
            Para_0.update({key: value},check_already_exists=False)
    return CyclePack,Para_0

def Read_ExpCrate(Path_Exp_Crate, ):
    #BasicPath = "D:/OneDrive - Imperial College London/SimDataSave/InputData/" 
    #Target = "Ruihe_newLGM50_Crate/"
    Cell_1_2p5to3C = pd.read_csv(
        Path_Exp_Crate + "dicharge_2.5C_3C_ch_a_1_CA1.txt", #engine='python',
        encoding = "shift-jis",  skiprows = 1,
        sep ='\t', header=None)
    Cell_1_2p5to3C.head()
    newNames = [
        "time/s", "Ns","Ecell/V", "I/mA", 
        "(Q-Qo)/mA.h", "Temperature/ｰC", 
        "Q charge/mA.h","Q discharge/mA.h","R/Ohm"]
    oldNames = np.arange(9).tolist()
    Cell_1_2p5to3C=Cell_1_2p5to3C.rename(columns={i:j for i,j in zip(oldNames,newNames)})
    Cell_1_2p5to3C.head()
    # Get Cell-1 2.5C and 3C:
    font = {'family' : 'DejaVu Sans','size'   : 14};   mpl.rc('font', **font)
    Cell_1_2p5C = Cell_1_2p5to3C[(Cell_1_2p5to3C['Ns']==5)]
    df_dc4 = Cell_1_2p5to3C[(Cell_1_2p5to3C['Ns']==12)]
    Cell_1_3C = Cell_1_2p5to3C[(Cell_1_2p5to3C['Ns']==12)]
    # Read Cell-2 2.5C to 3C:
    Cell_2_2p5to3C = pd.read_csv(
        Path_Exp_Crate + "dicharge_2.5C_3C_ch_a_2_CA2.txt", #engine='python',
        encoding = "shift-jis",  skiprows = 1,
        sep ='\t', header=None)
    Cell_2_2p5to3C.head()
    newNames = [
        "time/s", "Ns","Ecell/V", "I/mA", 
        "(Q-Qo)/mA.h", "Temperature/ｰC", 
        "Q charge/mA.h","Q discharge/mA.h","R/Ohm"]
    oldNames = np.arange(9).tolist()
    Cell_2_2p5to3C=Cell_2_2p5to3C.rename(columns={i:j for i,j in zip(oldNames,newNames)})
    Cell_2_2p5to3C.head()
    Cell_2_2p5C = Cell_2_2p5to3C[(Cell_2_2p5to3C['Ns']==5)]
    df_dc4 = Cell_2_2p5to3C[(Cell_2_2p5to3C['Ns']==12)]
    Cell_2_3C = Cell_2_2p5to3C[(Cell_2_2p5to3C['Ns']==12)]
    # Read Cell-1 3.5C to 4C:
    Cell_1_3p5to4C = pd.read_csv(
        Path_Exp_Crate + "Discharge_test_3p5_to_4C_A1_CA1.txt", #engine='python',
        encoding = "shift-jis",  skiprows = 1,
        sep ='\t', header=None)
    Cell_1_3p5to4C.head()
    newNames = [
        "time/s", "Ns","Ecell/V", "I/mA", 
        "(Q-Qo)/mA.h", "Temperature/ｰC", 
        "Q charge/mA.h","Q discharge/mA.h","R/Ohm"]
    oldNames = np.arange(9).tolist()
    Cell_1_3p5to4C=Cell_1_3p5to4C.rename(columns={i:j for i,j in zip(oldNames,newNames)})
    Cell_1_3p5to4C.head()
    Cell_1_3p5C = Cell_1_3p5to4C[(Cell_1_3p5to4C['Ns']==5)]
    df_dc4 = Cell_1_3p5to4C[(Cell_1_3p5to4C['Ns']==12)]
    Cell_1_4C = Cell_1_3p5to4C[(Cell_1_3p5to4C['Ns']==12)]
    # Read Cell-2 3.5C to 4C:
    Cell_2_3p5to4C = pd.read_csv(
        Path_Exp_Crate + "Discharge_test_3p5_to_4C_A2_CA3.txt", #engine='python',
        encoding = "shift-jis",  skiprows = 1,
        sep ='\t', header=None)
    Cell_2_3p5to4C.head()
    newNames = [
        "time/s", "Ns","Ecell/V", "I/mA", 
        "(Q-Qo)/mA.h", "Temperature/ｰC", 
        "Q charge/mA.h","Q discharge/mA.h","R/Ohm"]
    oldNames = np.arange(9).tolist()
    Cell_2_3p5to4C=Cell_2_3p5to4C.rename(columns={i:j for i,j in zip(oldNames,newNames)})
    Cell_2_3p5to4C.head()
    Cell_2_3p5C = Cell_2_3p5to4C[(Cell_2_3p5to4C['Ns']==5)]
    df_dc4 = Cell_2_3p5to4C[(Cell_2_3p5to4C['Ns']==12)]
    Cell_2_4C = Cell_2_3p5to4C[(Cell_2_3p5to4C['Ns']==12)]
    # Read Cell-2 up to 2C:
    Cell_2_UpTo2C = pd.read_csv(
        Path_Exp_Crate + "dicharge_to_2C_ch_a_2_CA2.txt", #engine='python',
        encoding = "shift-jis",  skiprows = 1,
        sep ='\t', header=None)
    Cell_2_UpTo2C.head()
    newNames = [
        "time/s", "Ns","Ecell/V", "I/mA", 
        "(Q-Qo)/mA.h", "Temperature/ｰC", 
        "Q charge/mA.h","Q discharge/mA.h","R/Ohm"]
    oldNames = np.arange(9).tolist()
    Cell_2_UpTo2C=Cell_2_UpTo2C.rename(columns={i:j for i,j in zip(oldNames,newNames)})
    Cell_2_UpTo2C.head()

    Cell_2_0p5C = Cell_2_UpTo2C[(Cell_2_UpTo2C['Ns']==5)]
    Cell_2_1C = Cell_2_UpTo2C[(Cell_2_UpTo2C['Ns']==12)]
    Cell_2_1p25C = Cell_2_UpTo2C[(Cell_2_UpTo2C['Ns']==19)]
    Cell_2_1p5C = Cell_2_UpTo2C[(Cell_2_UpTo2C['Ns']==26)]
    Cell_2_1p75C = Cell_2_UpTo2C[(Cell_2_UpTo2C['Ns']==33)]
    Cell_2_2C = Cell_2_UpTo2C[(Cell_2_UpTo2C['Ns']==40)]
    # Read Cell-1 up to 2C:

    Cell_1_UpTo2C = pd.read_csv(
        Path_Exp_Crate + "dicharge_to_2C_ch_a_1_CA1.txt", #engine='python',
        encoding = "shift-jis",  skiprows = 1,
        sep ='\t', header=None)
    Cell_1_UpTo2C.head()
    newNames = [
        "time/s", "Ns","Ecell/V", "I/mA", 
        "(Q-Qo)/mA.h", "Temperature/ｰC", 
        "Q charge/mA.h","Q discharge/mA.h","R/Ohm"]
    oldNames = np.arange(9).tolist()
    Cell_1_UpTo2C=Cell_1_UpTo2C.rename(columns={i:j for i,j in zip(oldNames,newNames)})
    Cell_1_UpTo2C.head()

    Cell_1_0p5C = Cell_1_UpTo2C[(Cell_1_UpTo2C['Ns']==5)]
    Cell_1_1C = Cell_1_UpTo2C[(Cell_1_UpTo2C['Ns']==12)]
    Cell_1_1p25C = Cell_1_UpTo2C[(Cell_1_UpTo2C['Ns']==19)]
    Cell_1_1p5C = Cell_1_UpTo2C[(Cell_1_UpTo2C['Ns']==26)]
    Cell_1_1p75C = Cell_1_UpTo2C[(Cell_1_UpTo2C['Ns']==33)]
    Cell_1_2C = Cell_1_UpTo2C[(Cell_1_UpTo2C['Ns']==40)]
    # collect all together
    Cell_2_All = [
        Cell_2_0p5C,Cell_2_1C,Cell_2_1p25C,Cell_2_1p5C,Cell_2_1p75C,
        Cell_2_2C,Cell_2_2p5C,Cell_2_3C,Cell_2_3p5C,Cell_2_4C];
    Cell_1_All = [
        Cell_1_0p5C,Cell_1_1C,Cell_1_1p25C,Cell_1_1p5C,Cell_1_1p75C,
        Cell_1_2C,Cell_1_2p5C,Cell_1_3C,Cell_1_3p5C,Cell_1_4C];
    str_Crate = ["0.5","1","1.25","1.5","1.75","2","2.5","3","3.5","4",]
    Num_Crate = [];
    for str in str_Crate:
        Num_Crate.append(float(str)) # print(Num_Crate)
    return Num_Crate, Cell_1_All,Cell_2_All

def Add_var(para_used,model):
    c_e = model.variables["Electrolyte concentration [mol.m-3]"]
    c_EC= model.variables["EC concentration [mol.m-3]"]
    T = model.variables["Cell temperature [K]"]
    D_e = para_used["Electrolyte diffusivity [m2.s-1]"]
    D_EC= para_used["EC diffusivity in electrolyte [m2.s-1]"]
    D_EC_e_cross = para_used["EC Lithium ion cross diffusivity [m2.s-1]"]
    sigma_e = para_used["Electrolyte conductivity [S.m-1]"]
    dLJP_dcEC = para_used["Measured dLJP_dcEC"] # dLJP_Two_Species_dco_Jung2023(x,y,T): # # ~~~~# x: ce; y: co 
    dLJP_dce  = para_used["Measured dLJP_dce"]
    Xi = para_used["EC transference number"]
    c_T = para_used["Total concentration [mol.m-3]"]
    model.variables["Electrolyte diffusivity [m2.s-1]"] = D_e(c_e,c_EC, T)
    model.variables["EC diffusivity in electrolyte [m2.s-1]"] = D_EC(c_e,c_EC, T)
    model.variables["EC Lithium ion cross diffusivity [m2.s-1]"] = D_EC_e_cross(c_e,c_EC, T)
    model.variables["Electrolyte conductivity [S.m-1]"] = sigma_e(c_e,c_EC, T)
    model.variables["EC transference number"] = Xi(c_e,c_EC, T)
    model.variables["Total concentration [mol.m-3]"] = c_T(c_e,c_EC, T)
    model.variables["c(EC) over c(Li+)"] = c_EC / c_e
    model.variables["dLJP_dcEC"] =  dLJP_dcEC(c_e,c_EC, T)
    model.variables["dLJP_dce"] =  dLJP_dce(c_e,c_EC, T)
    # molar mass in Taeho's paper: 
    M_EMC = 104.105*1e-3; M_EC = 88.062*1e-3; M_e = 151.905*1e-3; #   kg/mol
    c_EMC = model.variables["Total concentration [mol.m-3]"] - 2*c_e- c_EC
    model.variables["c(EMC) [mol.m-3]"] =  c_EMC
    
    model.variables["Total electrolyte concentration [mol.m-3]"] = model.variables["Total concentration [mol.m-3]"] 
    model.variables["y_e"] = c_e / model.variables["Total concentration [mol.m-3]"] 
    model.variables["y_EC"] = c_EC / model.variables["Total concentration [mol.m-3]"] 
    model.variables["EC:EMC wt%"] =  (c_EC*M_EC) / (c_EMC*M_EMC) 
    model.variables["EC:EMC %"] =  c_EC / c_EMC 
    t_0plus = para_used["Cation transference number"]
    model.variables["Cation transference number"] = t_0plus(c_e,c_EC, T)
    return model


def Run_P3_OneCycle_Dict(
        index_i, Para_dd_i, Path_pack,Rate_Dis,str_pickle,
        Return_sol,Save_sol):
    Timer = pybamm.Timer()
    count_i = int(index_i);
    print(f'Start Now! C rate: {Rate_Dis}')  
    if len(Path_pack)>2:
        [BasicPath,Target,Path_Exp_Crate,book_name_xlsx,sheet_name_xlsx,] = Path_pack
    else:
        [BasicPath,Target,] = Path_pack
    ##### Initialise Para_0 and model 
    CyclePack,para_used = Para_init_Dict(Para_dd_i)
    [Mesh_list,model_options] = CyclePack
    model = pybamm.lithium_ion.DFN(options=model_options)
    str_model_options = str(model_options)
    V_max = 4.2;        V_min = 2.5
    if Rate_Dis > 1:
        ts_dis = 0.1
    else: 
        ts_dis = 2
    Exp_1  = pybamm.Experiment(
    [ (
        f"Hold at {V_max} V until C/100",
        f"Discharge at {Rate_Dis} C until {V_min} V", #  ({ts_dis} second period)
        #"Rest for 1 hour",
        #f"Charge at 1 C until {V_max} V (2 second period)", 
        #f"Hold at {V_max} V until C/100"
        )    ] * 1 )  

    model = Add_var(para_used,model)
    var_pts = {
        "x_n": Mesh_list[0],  # negative electrode
        "x_s": Mesh_list[1],  # separator 
        "x_p": Mesh_list[2],  # positive electrode
        "r_n": Mesh_list[3],  # negative particle
        "r_p": Mesh_list[4],  # positive particle
    }
    sim    = pybamm.Simulation(
        model, experiment = Exp_1,
        parameter_values = para_used,
        solver = pybamm.CasadiSolver(return_solution_if_failed_early=True),
        var_pts=var_pts,
        )       
    sol    = sim.solve()
    if Return_sol:
        sol_return = sol
    else:
        sol_return = []
    if Save_sol:
        import pickle
        with open(
            BasicPath + Target
            + f'{str_pickle}_Rate={Rate_Dis}C.pkl', 'wb') as file: # str =  f"D_e_EC_cross={D_e_EC_cross}"
            pickle.dump(sol, file)
    else: 
        pass
    print(f'Finish {str_pickle}_Rate={Rate_Dis}C with Spent {Timer.time()}')
    Timer.reset()
    
    # Get temperature rise and discharge capacity as well
    step_i = sol.cycles[0].steps[1]
    cap = (
        step_i["Discharge capacity [A.h]"].entries[-1]
        -step_i["Discharge capacity [A.h]"].entries[0])
    Trise = (
        step_i["Volume-averaged cell temperature [C]"].entries[-1]
        -step_i["Volume-averaged cell temperature [C]"].entries[0])
    t_dis =( 
        step_i["Time [h]"].entries - step_i["Time [h]"].entries[0] 
        ) 
    vol_dis = step_i["Terminal voltage [V]"].entries 
    
    return sol_return,cap,Trise,t_dis,vol_dis


def Scan_Crate_Paper(
        index_i, Para_dd_i, Path_pack ,  str_model,
        Rate_Dis_All, str_pickle,
        Return_sol,Save_sol):   
    Case_Dict = {}
    print('Start Now! Scan %d.' % index_i)  
    Sol_All =[]; Cap_Dis_All = [];  Trise_All = []; 
    Time_dis_All =[];Vol_dis_All =[];
    try:
        for Rate_Dis in Rate_Dis_All:
            sol,cap,Trise,t_dis,vol_dis = Run_P3_OneCycle_Dict(
                index_i, Para_dd_i, Path_pack,Rate_Dis,
                str_pickle,Return_sol,Save_sol) 
            # (index_i, Para_dd_i, Path_pack,Rate_Dis,str_pickle,Return_sol,Save_sol):
            Sol_All.append(sol)
            Cap_Dis_All.append(cap)
            Trise_All.append(Trise)
            Time_dis_All.append(t_dis)
            Vol_dis_All.append(vol_dis)
        Case_Dict ['Cap_Dis_All'] = Cap_Dis_All
        Case_Dict ['Sol_All']  = Sol_All
        Case_Dict ['Rate_Dis_All']  = Rate_Dis_All
        Case_Dict ['Trise_All']  = Trise_All
        Case_Dict ['Time_dis_All']  = Time_dis_All
        Case_Dict ['Vol_dis_All']  = Vol_dis_All
    except:
        print(f"Something went wrong with {str_model} - Scan={index_i}")
        Case_Dict = "Empty";
    else:
        print(f"Finish {str_model} - Scan={index_i}")
    return Case_Dict


def Scan_Crate_Dict(
        index_i, Para_dd_i, Path_pack ,  str_model,
        Rate_Dis_All, Return_sol,SaveFig):   
    Case_Dict = {}
    print('Start Now! Scan %d.' % index_i)  
    Sol_All =[]; Cap_Dis_All = [];  Trise_All = []; 
    Time_dis_All =[];Vol_dis_All =[];
    try:
        for Rate_Dis in Rate_Dis_All:
            sol,cap,Trise,t_dis,vol_dis = Run_P3_OneCycle_Dict(
                index_i, Para_dd_i, Path_pack,Rate_Dis,
                Return_sol)
            Sol_All.append(sol)
            Cap_Dis_All.append(cap)
            Trise_All.append(Trise)
            Time_dis_All.append(t_dis)
            Vol_dis_All.append(vol_dis)
        Case_Dict ['Cap_Dis_All'] = Cap_Dis_All
        Case_Dict ['Sol_All']  = Sol_All
        Case_Dict ['Rate_Dis_All']  = Rate_Dis_All
        Case_Dict ['Trise_All']  = Trise_All
        Case_Dict ['Time_dis_All']  = Time_dis_All
        Case_Dict ['Vol_dis_All']  = Vol_dis_All
        # Plot capacity comparasion for each scan - TODO: voltage comparasion 
        #     Summarize Exp - plot temperature rise and capacity vs C rate
        if SaveFig == True:    
            for i in range(1):
                Niall_Crate = [0.2, 0.3, 0.4, 0.5, 1, 2, 3]; 
                Niall_Cap = [ 4.815, 4.75, 4.82, 4.82, 4.64, 3.298, 1.983]; 
                Ruihe_Crate = [0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0]
                RL_Cap_Cell_1 = [4.715963693555362, 4.612124591646047, 4.525131462728684, 4.357837218694934, 3.9985270604389145, 3.4233452957546677, 2.5758641034678815, 1.988704897546658, 1.6652388357649741, 1.4100456949750437]
                RL_Cap_Cell_2 = [4.716030959172494, 4.617948432702642, 4.5276396589819425, 4.341444372881664, 3.9091685871110027, 3.286849922924804, 2.4498175045751953, 1.901398288847656, 1.57583310209744, 1.345814775172526]
                # Biologic temperature:
                RL_T_Rise_Cell_1 = [1.8637965639005145, 2.6184089732547733, 2.2785406232761183, 2.7299726754489058, 3.4454458756196686, 3.2409775493440662, 4.36836993533397, 4.043625726318339, 4.889784013534822, 5.519928668485186]
                RL_T_Rise_Cell_2 = [1.6080960965963662, 1.5360009591223402, 2.038724388471504, 2.0249414408349473, 2.984494420154249, 3.5199766991718633, 4.342076997341994, 4.719743013070929, 3.79700883048751, 3.071455647673517]
                # from 1C to 3C only
                Pico_T_Rise_Cell_1 =[
                    26.643-25.982, 26.778-25.974,
                    27.133-25.889,27.481-25.917,27.6-25.9,
                    28.023-26,28.229-25.952]
                Pico_T_Rise_Cell_2 =[
                    26.411-25.706,26.671-25.681,27.025-25.611,
                    27.481-25.644,27.7-25.6,28.361-25.739, 28.832-25.633]
            
            ls = "-"
            # Compare experiment and modelling result: temperature rise and capacity vs C rate
            fig, axs = plt.subplots(1,2, figsize=(9.3,6.4),tight_layout=True)
            # experiment:
            axs[0].plot(Ruihe_Crate, RL_Cap_Cell_1 ,linestyle='none',marker ="o", color="gray",) # label="Cell-1"
            axs[0].plot(Ruihe_Crate, RL_Cap_Cell_2 ,linestyle='none',marker ="s",color="gray",) # label="Cell-2"
            axs[0].plot(Niall_Crate[3:], Niall_Cap[3:] ,linestyle='none',marker ="s",color="gray",) # label="Niall"
            axs[1].plot(Ruihe_Crate, RL_T_Rise_Cell_1 ,linestyle='none',marker ="o",color="gray",)  # label="Bio-Cell-1"
            axs[1].plot(Ruihe_Crate, RL_T_Rise_Cell_2 ,linestyle='none',marker ="s",color="gray",) # label="Bio-Cell-2"
            # simulation - double
            axs[0].plot( Rate_Dis_All,Cap_Dis_All,linestyle=ls,marker ="^",color="b",)
            axs[0].plot( Rate_Dis_All,Cap_Dis_All,linestyle=ls,marker ="^",color="b",)
            axs[1].plot( Rate_Dis_All,Trise_All,linestyle=ls,marker ="^",color="b",)
            axs[1].plot( Rate_Dis_All,Trise_All,linestyle=ls,marker ="^",color="b",)
            axs[0].set_ylabel("Capacity [A.h]",fontsize=fs)
            axs[1].set_ylabel("T rise [degC]",fontsize=fs)
            axs[0].set_xlabel("C rate",fontsize=fs-2)
            axs[1].set_xlabel("C rate",fontsize=fs-2)
            fig.suptitle(f'{str_model} - Scan = {index_i}', fontsize=fs)
            [BasicPath,Target,Path_Exp_Crate,
                book_name_xlsx,sheet_name_xlsx,] = Path_pack
            dpi = 600;
        
            plt.savefig( BasicPath + Target+
                f"{str_model} - Scan={index_i} Cap and Temperature rise.png", dpi=dpi)
            # compare voltage - need to read detailed csv file so must prepare:
            Num_Crate, Cell_1_All,Cell_2_All = Read_ExpCrate(Path_Exp_Crate, )

            cmdd = mpl.cm.get_cmap("cool", len(Rate_Dis_All)) 
            cmgray = mpl.cm.get_cmap("gray", len(Rate_Dis_All)) 
            fig, axs = plt.subplots( figsize=(9.3,6.2),tight_layout=True)
            for j in range(len(Rate_Dis_All)):
                # plot simulation:
                axs.plot(Case_Dict['Time_dis_All'][j],Case_Dict['Vol_dis_All'][j], color=cmdd(j),)
                # plot experiment:
                index = Num_Crate.index(Rate_Dis_All[j])
                axs.plot(
                    (Cell_1_All[index]["time/s"]-Cell_1_All[index]["time/s"].iloc[0])/3600,
                    Cell_1_All[index]["Ecell/V"],linestyle='--',
                    color=cmgray(j),)
            fig.suptitle(f'{str_model} - Scan = {index_i}', fontsize=fs)
        if SaveFig == True:   
            plt.savefig( BasicPath + Target+
                f"{str_model} - Scan={index_i} Vol compare.png", dpi=dpi)    
    except:
        print(f"Something went wrong with {str_model} - Scan={index_i}")
        Case_Dict = "Empty";
    else:
        print(f"Finish {str_model} - Scan={index_i}")
    return Case_Dict

def Plot_quick(Case_Para_All, Crate_i, fs):
    font = {'family' : 'DejaVu Sans','size'   : fs}; mpl.rc('font', **font)
    label = ["Normal DFN","Single transport","Double transport",] 
    output_variables3 = [
        #"X-averaged battery open circuit voltage [V]",
        "Battery voltage [V]",
        "X-averaged battery reaction overpotential [V]",
        "X-averaged battery concentration overpotential [V]",
        "X-averaged EC concentration overpotential [V]", # Mark Ruihe add
        "X-averaged battery electrolyte ohmic losses [V]",
        "X-averaged battery solid phase ohmic losses [V]",
    ]
    quick_plot = pybamm.QuickPlot(
        [
            Case_Para_All[i]['MyDict_All'][Crate_i]['Solution'].cycles[0].steps[1] for i in range(0,3)
        ], 
        output_variables3,label,variable_limits='tight',
        time_unit='hours',
        spatial_unit='mm',     #  (“m”, “mm”, or “um”)
        n_rows=2) #figsize = (18,12),
    quick_plot.dynamic_plot()

    output_variables3 = [
        "Electrolyte potential [V]",
        "EC concentration [mol.m-3]",
        "Electrolyte concentration [mol.m-3]",
        "EC transference number",
        "Electrolyte conductivity [S.m-1]",
        "Electrolyte diffusivity [m2.s-1]",
    ]
    quick_plot = pybamm.QuickPlot(
        [
            Case_Para_All[i]['MyDict_All'][Crate_i]['Solution'].cycles[0].steps[1] for i in range(0,3)
        ], 
        output_variables3,label,variable_limits='fixed',
        time_unit='hours',
        spatial_unit='mm',     #  (“m”, “mm”, or “um”)
        n_rows=2) #figsize = (18,12),
    quick_plot.dynamic_plot()

# this function is to initialize the para with a known dict
def Para_init(Para_dict):
    Para_dict_used = Para_dict.copy();
    Para_0=pybamm.ParameterValues(Para_dict_used["Para_Set"]  )
    Para_dict_used.pop("Para_Set")

    if Para_dict_used.__contains__("Total ageing cycles"):
        Total_Cycles = Para_dict_used["Total ageing cycles"]  
        Para_dict_used.pop("Total ageing cycles")
    else:
        Total_Cycles = False
    if Para_dict_used.__contains__("Mesh list"):
        Mesh_list = Para_dict_used["Mesh list"]  
        Para_dict_used.pop("Mesh list")
    else:
        Mesh_list = False
    if Para_dict_used.__contains__("SaveAsList"):
        SaveAsList = Para_dict_used["SaveAsList"]  
        Para_dict_used.pop("SaveAsList")
    else:
        SaveAsList = False
    if Para_dict_used.__contains__("Ageing cycles between RPT"):
        Cycle_bt_RPT = Para_dict_used["Ageing cycles between RPT"]  
        Para_dict_used.pop("Ageing cycles between RPT")
    if Para_dict_used.__contains__("Update cycles for ageing"):
        Update_Cycles = Para_dict_used["Update cycles for ageing"]  
        Para_dict_used.pop("Update cycles for ageing")
    if Para_dict_used.__contains__("Cycles within RPT"):
        RPT_Cycles = Para_dict_used["Cycles within RPT"]  
        Para_dict_used.pop("Cycles within RPT")
    if Para_dict_used.__contains__("Ageing temperature"):
        Temper_i = Para_dict_used["Ageing temperature"]  + 273.15
        Para_dict_used.pop("Ageing temperature")
    else: 
        Temper_i = False
    if Para_dict_used.__contains__("RPT temperature"):
        Temper_RPT = Para_dict_used["RPT temperature"]   + 273.15
        Para_dict_used.pop("RPT temperature")
    if Para_dict_used.__contains__("Particle mesh points"):
        mesh_par = Para_dict_used["Particle mesh points"]  
        Para_dict_used.pop("Particle mesh points")
    if Para_dict_used.__contains__("Exponential mesh stretch"):
        submesh_strech = Para_dict_used["Exponential mesh stretch"]  
        Para_dict_used.pop("Exponential mesh stretch")
    else:
        submesh_strech = "nan";
    if Para_dict_used.__contains__("Model option"):
        model_options = Para_dict_used["Model option"]  
        Para_dict_used.pop("Model option")
    else:
        model_options = False

    if Para_dict_used.__contains__("Initial Neg SOC"):
        c_Neg1SOC_in = (
            Para_0["Maximum concentration in negative electrode [mol.m-3]"]
            *Para_dict_used["Initial Neg SOC"]  )
        Para_0.update(
            {"Initial concentration in negative electrode [mol.m-3]":
            c_Neg1SOC_in})
        Para_dict_used.pop("Initial Neg SOC")
    if Para_dict_used.__contains__("Initial Pos SOC"):    
        c_Pos1SOC_in = (
            Para_0["Maximum concentration in positive electrode [mol.m-3]"]
            *Para_dict_used["Initial Pos SOC"]  )
        Para_0.update(
            {"Initial concentration in positive electrode [mol.m-3]":
            c_Pos1SOC_in})
        Para_dict_used.pop("Initial Pos SOC")

    CyclePack = [ 
        Total_Cycles,SaveAsList,Mesh_list,Temper_i,model_options];
    
    for key, value in Para_dict_used.items():
        # risk: will update parameter that doesn't exist, 
        # so need to make sure the name is right 
        if isinstance(value, str):
            Para_0.update({key: eval(value)})
            #Para_dict_used.pop(key)
        else:
            Para_0.update({key: value},check_already_exists=False)
    return CyclePack,Para_0

def GetSol_dict (my_dict, keys_all, Sol, 
    cycle_no,step_CD, step_CC , step_RE, step_CV ):

    [keys_loc,keys_tim,keys_cyc]  = keys_all;
    # get time_based variables:
    if len(keys_tim): 
        for key in keys_tim:
            step_no = eval("step_{}".format(key[0:2]))
            if key[3:] == "Time [h]":
                my_dict[key].append  (
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries
                    -
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries[0])
            else:
                #print("Solve up to Step",step_no)
                my_dict[key].append  (
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries)
    # get cycle_step_based variables:
    if len(keys_cyc): 
        for key in keys_cyc:
            if key in ["Discharge capacity [A.h]"]:
                step_no = step_CD;
                my_dict[key].append  (
                    Sol.cycles[cycle_no].steps[step_no][key].entries[-1]
                    - 
                    Sol.cycles[cycle_no].steps[step_no][key].entries[0])
            elif key[0:5] in ["CDend","CCend","CVend","REend",
                            "CDsta","CCsta","CVsta","REsta",]:
                step_no = eval("step_{}".format(key[0:2]))
                if key[2:5] == "sta":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[0])
                elif key[2:5] == "end":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[-1])
                
    # get location_based variables:
    if len(keys_loc): 
        for key in keys_loc:
            if key in ["x_n [m]","x [m]","x_s [m]","x_p [m]"]:
                #print("These variables only add once")
                if not len(my_dict[key]):   # special: add only once
                    my_dict[key] = Sol[key].entries[:,-1]
            elif key[0:5] in ["CDend","CCend","CVend","REend",
                            "CDsta","CCsta","CVsta","REsta",]:      
                #print("These variables add multiple times")
                step_no = eval("step_{}".format(key[0:2]))
                if key[2:5] == "sta":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[:,0])
                elif key[2:5] == "end":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[:,-1])
    return my_dict 

# Run model   # mark:
def Run_P3_model(
    index_i, Para_dict_i,   Path_pack , 
    keys_all_AGE,   Exp_AGE_List, exp_index_pack ):

    ModelTimer = pybamm.Timer()
    Para_dict_old = Para_dict_i.copy();
    count_i = int(index_i);
    print('Start Now! Scan %d.' % count_i)  

    # Un-pack data:
    [cycle_no,step_AGE_CD,
        step_AGE_CC,step_AGE_CV, ] = exp_index_pack
    [exp_AGE,exp_AGE,exp_AGE_2, # seems not used but actually used!
        exp_AGE_CD,exp_AGE_CC,exp_AGE_CV] = Exp_AGE_List
    [BasicPath,Target,
        book_name_xlsx,sheet_name_xlsx,] = Path_pack

    ##### Initialise Para_0 and model 
    CyclePack,Para_0 = Para_init(Para_dict_i)
    [Total_Cycles,SaveAsList,Mesh_list,
        Temper_i,model_options] = CyclePack
    model = pybamm.lithium_ion.DFN(options=model_options)
    str_model_options = str(model_options)
    # add electrolyte properties as variables, 
    # only when they are not constant
    c_e = model.variables["Electrolyte concentration [mol.m-3]"]
    c_EC= model.variables["EC concentration [mol.m-3]"]
    T = model.variables["Cell temperature [K]"]
    D_e = Para_0["Electrolyte diffusivity [m2.s-1]"]
    D_EC= Para_0["EC diffusivity in electrolyte [m2.s-1]"]
    sigma_e = Para_0["Electrolyte conductivity [S.m-1]"]
    dLJP_dcEC = Para_0["Measured dLJP_dcEC"] # dLJP_Two_Species_dco_Jung2023(x,y,T): # # ~~~~# x: ce; y: co 
    dLJP_dce  = Para_0["Measured dLJP_dce"]
    Xi = Para_0["EC transference number"]
    c_T = Para_0["Total concentration [mol.m-3]"]
    model.variables["Electrolyte diffusivity [m2.s-1]"] = D_e(c_e,c_EC, T)
    model.variables["EC diffusivity in electrolyte [m2.s-1]"] = D_EC(c_e,c_EC, T)
    model.variables["Electrolyte conductivity [S.m-1]"] = sigma_e(c_e,c_EC, T)
    model.variables["EC transference number"] = Xi(c_e,c_EC, T)
    model.variables["Total concentration [mol.m-3]"] = c_T(c_e,c_EC, T)
    model.variables["c(EC) over c(Li+)"] = c_EC / c_e
    model.variables["dLJP_dcEC"] =  dLJP_dcEC(c_e,c_EC, T)
    model.variables["dLJP_dce"] =  dLJP_dce(c_e,c_EC, T)
    # molar mass in Taeho's paper: unit: g/mol
    M_EMC = 104.105; M_EC = 88.062; M_e = 151.905;
    c_EMC = model.variables["Total concentration [mol.m-3]"] - c_EC - 2*c_e
    model.variables["c(EMC) [mol.m-3]"] =  c_EMC
    model.variables["EC:EMC wt%"] =  (c_EC*M_EC) / (c_EMC*M_EMC) 
    t_0plus = Para_0["Cation transference number"]
    model.variables["Cation transference number"] = t_0plus(c_e,c_EC, T)
    
    # Define experiment- for ageing only, NOT RPT
    str_exp_All = [];Experiment_All = [];
    for SaveAs_i, exp_i in zip(SaveAsList,Exp_AGE_List):
        str_exp_All.append(str(exp_i))
        Experiment_All.append(
            pybamm.Experiment( exp_i * int(SaveAs_i)  ) )


    #####Important: index canot be pre-determined anymore! ######
    if not Mesh_list:
        var_pts = {
            "x_n": 5,  # negative electrode
            "x_s": 5,  # separator 
            "x_p": 5,  # positive electrode
            "r_n": 60,  # negative particle
            "r_p": 20,  # positive particle
        }
    else:
        var_pts = {
        "x_n": Mesh_list[0],  # negative electrode
        "x_s": Mesh_list[1],  # separator 
        "x_p": Mesh_list[2],  # positive electrode
        "r_n": Mesh_list[3],  # negative particle
        "r_p": Mesh_list[4],  # positive particle
        }    
    # initialize my_dict for outputs
    my_dict_AGE = {}; 
    for keys in keys_all_AGE:
        for key in keys:
            my_dict_AGE[key]=[];	
    ###### Biggest loop:  #####
    Succ_Cyc = []; Sol_All = []; 
    i = 0;step_switch = 0 # initialize the loop
    while (i < Total_Cycles):
        print('try to run %d cycles' % SaveAsList[step_switch])
        try:
            Timelimit = 3600*4;
            # the following turns on for HPC only!
            #signal.signal(signal.SIGALRM, handle_signal)
            #signal.alarm(Timelimit)
            
            if i==0: # 1st time or never succeed, run from scratch:
                rioCall = RioCallback()  # define callback
                sim_0    = pybamm.Simulation(
                    model,
                    experiment = Experiment_All[step_switch],
                    parameter_values = Para_0,
                    solver = pybamm.CasadiSolver(),
                    var_pts=var_pts,
                    #submesh_types=submesh_types
                    ) 
                sol_0    = sim_0.solve(
                    calc_esoh=False,
                    save_at_cycles = SaveAsList[step_switch],
                    callbacks=rioCall)
            else: # succeeded before, 
                model_new = model_old.set_initial_conditions_from(
                    sol_old, inplace=False)
                rioCall = RioCallback()  # define callback
                sim_new    = pybamm.Simulation(
                    model_new,        
                    experiment =  Experiment_All[step_switch],
                    parameter_values = Para_0,
                    solver = pybamm.CasadiSolver(),
                    var_pts=var_pts,
                    #submesh_types=submesh_types
                    ) 
                sol_new    = sim_new.solve(
                    calc_esoh=False,
                    save_at_cycles = SaveAsList[step_switch],
                    callbacks=rioCall)   
            if rioCall.success == False:
                1/0
            #except:
            # the following turns on for HPC only!
            """      except (
            TimeoutError,
            ZeroDivisionError,
            pybamm.expression_tree.exceptions.ModelError,
            pybamm.expression_tree.exceptions.SolverError
            ) as e: """
        except (
            ZeroDivisionError,
            pybamm.expression_tree.exceptions.ModelError,
            pybamm.expression_tree.exceptions.SolverError
            ) as e:
            print('Failed or took too long, shorten cycles')
            print(e)
            step_switch += 1
            if (step_switch >= len(SaveAsList)):
                print('Exit as no options left')
                str_error = traceback.format_exc() # final errors 
                print('Finally finish %d cycles' % i)  
                break
        else:        
            if i == 0: 
                model_old = model; sol_old = sol_0    
            else: 
                model_old = model_new; sol_old = sol_new
                del model_new,sol_new
            ### Should do early post-prosessing and delete to 
            ### save RAM in the near future, not now 
            # because how to post-processing depends on step_switch
            Sol_All.append(sol_old)
            Succ_Cyc.append(SaveAsList[step_switch])
            i += SaveAsList[step_switch]
            print('Succeed! Now it is the %dth cycles' % i)  
            if step_switch <= 2:
                step_switch += 0
            elif step_switch > 2 and step_switch < len(SaveAsList)-1:
                step_switch += 1
                print('Succeed a single step and switch to next step normally')
            else:
                print('Finish last single step and Exit as no options left')
                print('Finally finish %d cycles' % i)  
                break
    ###########################################        
    #    Post-prosessing start from here:     #
    ###########################################
    if not os.path.exists(BasicPath + Target + f"{count_i}th Scan/"):
            os.mkdir(BasicPath + Target+  f"{count_i}th Scan/")
        
    j=0;
    while j <len(Sol_All):
        if len(Sol_All[j].cycles[-1].steps)==1:
            break
        j += 1
    if j < len(Sol_All):
        print("Single step starts from %d" %j)
    elif j==len(Sol_All):
        print("No single step")
    else:
        pass
    
    Sol_all_i = Sol_All
    Succ_Cyc_acc_i = np.cumsum(Succ_Cyc).tolist()
    # Plot for last steps
    if j<len(Sol_all_i):
        print("Not all solution has full cycles")
        for j in range(j,len(Sol_all_i)):
            index_cyc = Succ_Cyc_acc_i[j];
            Plot_Last_Single_Step(
                Sol_all_i[j],0,0,BasicPath, 
                Target,count_i,index_cyc,"True","cool",17,200)

    if j > 2 or j==len(Sol_All):   # normal post-proessing, 
        step_RPT_RE = -1
        ###########################################        
        #    11111111111111111111111111111111     #
        ###########################################        
        my_dict_AGE = {}; 
        for keys in keys_all_AGE:
            for key in keys:
                my_dict_AGE[key]=[];	
        
        Full_cycle = []
        # post-prosessing for full cycle range
        # prosess for the 1st solution: how many cycles do you have?
        if j<len(Sol_all_i):
            print("Not all solution has full cycles")
            for m in range(j,len(Sol_all_i)):
                index_cyc = Succ_Cyc_acc_i[m];
                Plot_Last_Single_Step(
                    Sol_all_i[m],0,0,BasicPath, 
                    Target,count_i,index_cyc,"True","cool",17,200)

        for m in range(0,j):# post-prosess for normal full cycles
            if m == 0 and Succ_Cyc_acc_i[m]>1: # first solution:
                # get two solution
                try:
                    cycle_no=0
                    my_dict_AGE = GetSol_dict (
                        my_dict_AGE,keys_all_AGE, Sol_all_i[m], 
                        cycle_no, step_AGE_CD , step_AGE_CC , 
                        step_RPT_RE, step_AGE_CV   ) 
                    Full_cycle.append(0)
                        
                    cycle_no=Succ_Cyc[m]-1
                    my_dict_AGE = GetSol_dict (
                        my_dict_AGE,keys_all_AGE, Sol_all_i[m], 
                        cycle_no, step_AGE_CD , step_AGE_CC , 
                        step_RPT_RE, step_AGE_CV   ) 
                    Full_cycle.append(Succ_Cyc_acc_i[m])
                except:
                    pass
                else:
                    print("Seems no empty solution")     
            else:
                # get only one solution
                cycle_no=Succ_Cyc[m]-1      
                try:  #   possibly have an empty solution in the middle!
                    my_dict_AGE = GetSol_dict (
                        my_dict_AGE,keys_all_AGE, Sol_all_i[m], 
                        cycle_no, step_AGE_CD , step_AGE_CC , 
                        step_RPT_RE, step_AGE_CV   ) 
                    Full_cycle.append(Succ_Cyc_acc_i[m])
                except:
                    pass
                else:
                    print("Seems no empty solution")        
        # add a bit more to my_dict_AGE
        my_dict_AGE["CDcyc Positive SOC range"]=(
            abs(
            np.array(my_dict_AGE["CDsta Positive electrode SOC"])-
            np.array(my_dict_AGE["CDend Positive electrode SOC"])
        )).tolist()
        my_dict_AGE["CDcyc Negative SOC range"]=(
            abs(
            np.array(my_dict_AGE["CDsta Negative electrode SOC"])-
            np.array(my_dict_AGE["CDend Negative electrode SOC"])
        )
        ).tolist()
        my_dict_AGE["LLI to SEI in one CD step [A.h]"]=(
            np.array(my_dict_AGE["CDend Loss of capacity to SEI [A.h]"])-
            np.array(my_dict_AGE["CDsta Loss of capacity to SEI [A.h]"])
        ).tolist()
        my_dict_AGE["LLI to SEI in one CC step [A.h]"]=(
            np.array(my_dict_AGE["CCend Loss of capacity to SEI [A.h]"])-
            np.array(my_dict_AGE["CCsta Loss of capacity to SEI [A.h]"])
        ).tolist()
        my_dict_AGE["LLI to SEI in one CV step [A.h]"]=(
            np.array(my_dict_AGE["CVend Loss of capacity to SEI [A.h]"])-
            np.array(my_dict_AGE["CVsta Loss of capacity to SEI [A.h]"])
        ).tolist()
        my_dict_AGE["LLI to SEI in one Charge step [A.h]"]=(
            np.array(my_dict_AGE["LLI to SEI in one CV step [A.h]"]) +
            np.array(my_dict_AGE["LLI to SEI in one CC step [A.h]"])
        ).tolist()
        # Add for LLI per Ah:
        my_dict_AGE["CC Charge capacity [A.h]"]=(
            abs(
            np.array(my_dict_AGE["CCsta Discharge capacity [A.h]"])-
            np.array(my_dict_AGE["CCend Discharge capacity [A.h]"])
        )).tolist()
        my_dict_AGE["CV Charge capacity [A.h]"]=(
            abs(
            np.array(my_dict_AGE["CVsta Discharge capacity [A.h]"])-
            np.array(my_dict_AGE["CVend Discharge capacity [A.h]"])
        )).tolist()
        my_dict_AGE["Charge capacity [A.h]"]=(
            np.array(my_dict_AGE["CC Charge capacity [A.h]"]) + 
            np.array(my_dict_AGE["CV Charge capacity [A.h]"])
        ).tolist()
        my_dict_AGE["CC LLI to SEI per Ah"]=(
            np.array(my_dict_AGE["LLI to SEI in one CC step [A.h]"])/
            np.array(my_dict_AGE["CC Charge capacity [A.h]"])
        ).tolist()
        my_dict_AGE["CV LLI to SEI per Ah"]=(
            np.array(my_dict_AGE["LLI to SEI in one CV step [A.h]"])/
            np.array(my_dict_AGE["CV Charge capacity [A.h]"])
        ).tolist()
        my_dict_AGE["Dis LLI to SEI per Ah"]=(
            np.array(my_dict_AGE["LLI to SEI in one CD step [A.h]"])/
            np.array(my_dict_AGE["Discharge capacity [A.h]"])
        ).tolist()
        my_dict_AGE["Cha LLI to SEI per Ah"]=(
            np.array(my_dict_AGE["CC LLI to SEI per Ah"]) +
            np.array(my_dict_AGE["CV LLI to SEI per Ah"])
        ).tolist()
        ###########################################        
        #    11111111111111111111111111111111     #
        ###########################################   
        ###########################################        
        #    22222222222222222222222222222222     #
        ###########################################     
        # Plot fig 1~3:
        fs = 17;  dpi=200;
        key_all_CCend = [
            "CCend Negative electrode porosity",
            "CCend Electrolyte concentration [mol.m-3]",
            "CCend EC concentration [mol.m-3]",
            "CCend Electrolyte potential [V]",
            "CCend Positive electrode potential [V]",
            "CCend Electrolyte current density [A.m-2]",
            "CCend Electrolyte diffusivity [m2.s-1]",
            "CCend Electrolyte conductivity [S.m-1]",
            "CCend Negative electrode SEI interfacial current density [A.m-2]",
        ]
        key_all_CDend = [
            "CDend Negative electrode porosity",
            "CDend Electrolyte concentration [mol.m-3]",
            "CDend EC concentration [mol.m-3]",
            "CDend Electrolyte potential [V]",
            "CDend Positive electrode potential [V]",
            "CDend Electrolyte current density [A.m-2]",
            "CDend Electrolyte diffusivity [m2.s-1]",
            "CDend Electrolyte conductivity [S.m-1]",
            "CDend Negative electrode SEI interfacial current density [A.m-2]",
        ]

        try: 
            Plot_Fig_1(Full_cycle,my_dict_AGE,
                BasicPath, Target,   count_i,  fs,  dpi)
            fig, axs = Plot_Loc_Var_2(
                Full_cycle,key_all_CCend,my_dict_AGE,fs)
            plt.savefig(
                BasicPath + Target+f"{count_i}th Scan/" +
                "Fig. 2 - CCend Loc based overall.png", dpi=dpi)
            fig, axs = Plot_Loc_Var_2(
                Full_cycle,key_all_CDend,my_dict_AGE,fs)
            plt.savefig(
                BasicPath + Target+f"{count_i}th Scan/" +
                "Fig. 3 - CDend Loc based overall.png", dpi=dpi)
        except:
            print(f"Something went wrong during plotting Fig. 1~3 for scan {count_i}")
        else:
            pass
        ###########################################        
        #    22222222222222222222222222222222     #
        ###########################################   
        ###########################################        
        #    33333333333333333333333333333333     #
        ###########################################   
        # write into excel:
        str_exp_AGE_text = str(exp_AGE)
        value_list_temp = list(Para_dict_i.values())
        values = []
        for value_list_temp_i in value_list_temp:
            values.append(str(value_list_temp_i))
        values.insert(0,str(count_i));
        #"Cap Loss","LLI to SEI",
        #"LAM to Neg","LAM to Pos",
        #"Error"])
        try: 
            values.extend([
                str_exp_AGE_text,
                str(my_dict_AGE["Discharge capacity [A.h]"][0] 
                - 
                my_dict_AGE["Discharge capacity [A.h]"][-1]),

                str(my_dict_AGE["CDend Loss of capacity to SEI [A.h]"][-1]),
                str(my_dict_AGE["CDend Negative electrode capacity [A.h]"][0] 
                - 
                my_dict_AGE["CDend Negative electrode capacity [A.h]"][-1]),

                str(my_dict_AGE["CDend Positive electrode capacity [A.h]"][0] 
                - 
                my_dict_AGE["CDend Positive electrode capacity [A.h]"][-1]),
            ])
            values = [values,]
            book_name_xlsx_seperate =   str(count_i)+ '_' + book_name_xlsx;
            sheet_name_xlsx =  str(count_i);
            write_excel_xlsx(
                BasicPath + Target+   book_name_xlsx_seperate, 
                sheet_name_xlsx, values)
            
            mdic_cycles = {
                "Full_cycle": Full_cycle,
            }
            midc_merge = {**my_dict_AGE, **mdic_cycles}
            savemat(BasicPath + Target+f"{count_i}th Scan/" +f"{count_i}th Scan" + '-for_AGE_only.mat',midc_merge)  
        except:
            print(f"Something went wrong during saving .mat for scan {count_i}")
            midc_merge = {}
        else: 
            print(f"Scan No. {index_i} succeed! Saving succeed as well!")

    else:
        my_dict_AGE={};mdic_cycles={};
        midc_merge = {}
        str_exp_AGE_text = str(exp_AGE)
        value_list_temp = list(Para_dict_i.values())
        values = []
        for value_list_temp_i in value_list_temp:
            values.append(str(value_list_temp_i))
        values.insert(0,str(count_i));
        str_error = traceback.format_exc()      
        values.extend([
            str_exp_AGE_text,
            "NaN",
            "NaN",
            "NaN",
            "NaN",
            str_error,
        ])
        values = [values,]
        book_name_xlsx_seperate =   str(count_i)+ '_' + book_name_xlsx;
        sheet_name_xlsx =  str(count_i);
        write_excel_xlsx(
            BasicPath + Target+   book_name_xlsx_seperate, 
            sheet_name_xlsx, values)
        print(f"Scan No. {index_i} degrade too fast!" )

    return Sol_all_i,j,midc_merge



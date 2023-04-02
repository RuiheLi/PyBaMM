import pybamm;import pandas as pd   ;import numpy as np;import os;
import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;
import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'
from textwrap import wrap
import openpyxl
import traceback
import random;import time, signal

from pybamm import tanh,exp,sqrt

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
        (c_e > c_e_constant) * electrolyte_conductivity_base_Landesfeind2019(c_e_constant,c_EC, T, coeffs)
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

# add Ruihe Li update 230315
def dLJP_One_Specie_dce_Jung2023(ce,co,T): # co means c_EC here
    # Eq. (13):
    R = 8.31446261815324;  F = 96485.3321
    c_tot = 1.379*ce+1.113e4

    aln = 1.390; a0 = 1.158; a1 = -8.955; a2 = 164.7
    ddelta_U_dce = R*T/F*(
        aln / ce  + a1/c_tot  + 2*a2*ce/c_tot**2
    )
    return ddelta_U_dce

def dLJP_Two_Species_dco_Jung2023(x,y,T): # # ~~~~# x: ce; y: co 
    # co means c_EC here
    # T = 298.15;     # need to be a variable, unit: K
    dLJP_dco = -6.46958282530615e-13*T*x*y*(4.99846855876465e-6*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 88.12*x/(1.4631*x + 0.3589*y + 9778) + 3.024*pybamm.log(x/(1.4631*x + 0.3589*y + 9778)) + 8.233)/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)**2*(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2*(1.4631*x + 0.3589*y + 9778)) + 8.61733326422705e-5*T*y*(-3.36342232204145e-7*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 2.55635111753373e-8*y/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 6.81/(1.4631*x + 0.3589*y + 9778))*(-2*x/(1.4631*x + 0.3589*y + 9778) - y/(1.4631*x + 0.3589*y + 9778) + 1)/(1.4631*x + 0.3589*y + 9778) + 8.61733326422705e-5*T*y*(7.50763911169966e-9*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 3.75381955584983e-9*y/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 1/(1.4631*x + 0.3589*y + 9778))*(89.6*x/(1.4631*x + 0.3589*y + 9778) + 6.81*y/(1.4631*x + 0.3589*y + 9778) - 12.6)/(1.4631*x + 0.3589*y + 9778) - 3.23479141265307e-13*T*y*(-2*x/(1.4631*x + 0.3589*y + 9778) - y/(1.4631*x + 0.3589*y + 9778) + 1)*(89.6*x/(1.4631*x + 0.3589*y + 9778) + 6.81*y/(1.4631*x + 0.3589*y + 9778) - 12.6)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 8.61733326422705e-5*T*y*(-3.6693605353664e-10*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**3 + 3.30786579261487e-7*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 1.13515503368899e-8*(1.4631*x + 0.3589*y + 9778)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2)/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(1.4631*x + 0.3589*y + 9778)) - 3.23479141265307e-13*T*y*(4.99846855876465e-6*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 88.12*x/(1.4631*x + 0.3589*y + 9778) + 3.024*pybamm.log(x/(1.4631*x + 0.3589*y + 9778)) + 8.233)/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2) + 8.61733326422705e-5*T*(-y/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(1.4631*x + 0.3589*y + 9778)) + 1)*(4.80685670883901e-14*x**3/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**4 - 1.8189401775022e-10*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**3 - 6.95958145654558e-7*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 8.10074260152393e-8*(x/(1.4631*x + 0.3589*y + 9778))**0.5*(1.4631*x + 0.3589*y + 9778)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 2.26580548391096e-6*(x/(1.4631*x + 0.3589*y + 9778))**1.5*(1.4631*x + 0.3589*y + 9778)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 2.38086005329775e-6*(x/(1.4631*x + 0.3589*y + 9778))**2.5*(1.4631*x + 0.3589*y + 9778)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 3.29641664296953e-5*(x/(1.4631*x + 0.3589*y + 9778))**3.5*(1.4631*x + 0.3589*y + 9778)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 4.87508545718217e-5*(x/(1.4631*x + 0.3589*y + 9778))**4.5*(1.4631*x + 0.3589*y + 9778)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 4.4069841585677e-9*(1.4631*x + 0.3589*y + 9778)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2) + 8.61733326422705e-5*T*(7.50763911169966e-9*x*y/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)**2*(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2*(1.4631*x + 0.3589*y + 9778)) + 3.75381955584983e-9*y/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2) - 1/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(1.4631*x + 0.3589*y + 9778)))*(-4.36532412919363e-10*x**3/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**3 + 2.47779284697917e-6*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 185.4*x/(1.4631*x + 0.3589*y + 9778) - 43.16*(x/(1.4631*x + 0.3589*y + 9778))**0.5 - 402.4*(x/(1.4631*x + 0.3589*y + 9778))**1.5 + 253.7*(x/(1.4631*x + 0.3589*y + 9778))**2.5 + 2509*(x/(1.4631*x + 0.3589*y + 9778))**3.5 - 2886*(x/(1.4631*x + 0.3589*y + 9778))**4.5 + 1.174*pybamm.log(x/(1.4631*x + 0.3589*y + 9778)) + 7.167) + 8.61733326422705e-5*T*(-2*x/(1.4631*x + 0.3589*y + 9778) - y/(1.4631*x + 0.3589*y + 9778) + 1)*(89.6*x/(1.4631*x + 0.3589*y + 9778) + 6.81*y/(1.4631*x + 0.3589*y + 9778) - 12.6)/(1.4631*x + 0.3589*y + 9778) + 8.61733326422705e-5*T*(4.99846855876465e-6*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 88.12*x/(1.4631*x + 0.3589*y + 9778) + 3.024*pybamm.log(x/(1.4631*x + 0.3589*y + 9778)) + 8.233)/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(1.4631*x + 0.3589*y + 9778))
    return dLJP_dco  # units: V

def dLJP_Two_Species_dce_Jung2023(x,y,T):
    # x: ce; y: co 
    dLJP_dce = 8.61733326422705e-5*T*y*(-3.06058143893223e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 2/(1.4631*x + 0.3589*y + 9778))*(4.99846855876465e-6*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 88.12*x/(1.4631*x + 0.3589*y + 9778) + 3.024*pybamm.log(x/(1.4631*x + 0.3589*y + 9778)) + 8.233)/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)**2*(1.4631*x + 0.3589*y + 9778)) + 8.61733326422705e-5*T*y*(-1.37114048464164e-6*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 1.04212797995642e-7*y/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 89.6/(1.4631*x + 0.3589*y + 9778))*(-2*x/(1.4631*x + 0.3589*y + 9778) - y/(1.4631*x + 0.3589*y + 9778) + 1)/(1.4631*x + 0.3589*y + 9778) + 8.61733326422705e-5*T*y*(3.06058143893223e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 1.53029071946611e-8*y/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 2/(1.4631*x + 0.3589*y + 9778))*(89.6*x/(1.4631*x + 0.3589*y + 9778) + 6.81*y/(1.4631*x + 0.3589*y + 9778) - 12.6)/(1.4631*x + 0.3589*y + 9778) - 1.31870251207933e-12*T*y*(-2*x/(1.4631*x + 0.3589*y + 9778) - y/(1.4631*x + 0.3589*y + 9778) + 1)*(89.6*x/(1.4631*x + 0.3589*y + 9778) + 6.81*y/(1.4631*x + 0.3589*y + 9778) - 12.6)/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 8.61733326422705e-5*T*y*(-1.49585996079537e-9*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**3 + 1.13454292995228e-5*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 88.12/(1.4631*x + 0.3589*y + 9778) + 3.024*(-1.53029071946611e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 1/(1.4631*x + 0.3589*y + 9778))*(1.4631*x + 0.3589*y + 9778)/x)/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(1.4631*x + 0.3589*y + 9778)) - 1.31870251207933e-12*T*y*(4.99846855876465e-6*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 - 88.12*x/(1.4631*x + 0.3589*y + 9778) + 3.024*pybamm.log(x/(1.4631*x + 0.3589*y + 9778)) + 8.233)/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2) + 8.61733326422705e-5*T*(-y/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(1.4631*x + 0.3589*y + 9778)) + 1)*(1.95957426879419e-13*x**3/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**4 - 2.05111057776714e-9*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**3 + 2.11842670006817e-6*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 185.4/(1.4631*x + 0.3589*y + 9778) - 43.16*(x/(1.4631*x + 0.3589*y + 9778))**0.5*(-7.65145359733057e-9*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 0.5/(1.4631*x + 0.3589*y + 9778))*(1.4631*x + 0.3589*y + 9778)/x - 402.4*(x/(1.4631*x + 0.3589*y + 9778))**1.5*(-2.29543607919917e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 1.5/(1.4631*x + 0.3589*y + 9778))*(1.4631*x + 0.3589*y + 9778)/x + 253.7*(x/(1.4631*x + 0.3589*y + 9778))**2.5*(-3.82572679866528e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 2.5/(1.4631*x + 0.3589*y + 9778))*(1.4631*x + 0.3589*y + 9778)/x + 2509*(x/(1.4631*x + 0.3589*y + 9778))**3.5*(-5.3560175181314e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 3.5/(1.4631*x + 0.3589*y + 9778))*(1.4631*x + 0.3589*y + 9778)/x - 2886*(x/(1.4631*x + 0.3589*y + 9778))**4.5*(-6.88630823759751e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 4.5/(1.4631*x + 0.3589*y + 9778))*(1.4631*x + 0.3589*y + 9778)/x + 1.174*(-1.53029071946611e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 1/(1.4631*x + 0.3589*y + 9778))*(1.4631*x + 0.3589*y + 9778)/x) + 8.61733326422705e-5*T*(-y*(-3.06058143893223e-8*x/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 2/(1.4631*x + 0.3589*y + 9778))/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)**2*(1.4631*x + 0.3589*y + 9778)) + 1.53029071946611e-8*y/((-2*x/(1.4631*x + 0.3589*y + 9778) + 1)*(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2))*(-4.36532412919363e-10*x**3/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**3 + 2.47779284697917e-6*x**2/(0.000149631826549397*x + 3.67048476170996e-5*y + 1)**2 + 185.4*x/(1.4631*x + 0.3589*y + 9778) - 43.16*(x/(1.4631*x + 0.3589*y + 9778))**0.5 - 402.4*(x/(1.4631*x + 0.3589*y + 9778))**1.5 + 253.7*(x/(1.4631*x + 0.3589*y + 9778))**2.5 + 2509*(x/(1.4631*x + 0.3589*y + 9778))**3.5 - 2886*(x/(1.4631*x + 0.3589*y + 9778))**4.5 + 1.174*pybamm.log(x/(1.4631*x + 0.3589*y + 9778)) + 7.167)
    return dLJP_dce # units: V

def dLJP_dcEC_Nyman_2011(c_e, c_EC , T):
    dLJP_dcEC = -5.394e-6 - 3.616e-2 / c_EC
    return dLJP_dcEC
def dLJP_dce_Nyman_2011(c_e, c_EC , T):
    dLJP_dce = 5.326e-5 + 2.47e-2 / c_e
    return dLJP_dce

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
    MyDict_All =[]; Cap_Dis_All = []  ; Cap_Cha_All = []  
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
    Case_Dict ['Cap_Dis_All'] = Cap_Dis_All
    Case_Dict ['Cap_Cha_All'] = Cap_Cha_All
    Case_Dict ['str_model'] = str_model
    Case_Dict ['str_para']  = str_para
    Case_Dict ['MyDict_All']  = MyDict_All
    Case_Dict ['Rate_Dis_All']  = Rate_Dis_All
    Case_Dict ['Rate_Cha_All']  = Rate_Cha_All
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
    if Para_dict_used.__contains__("SaveAsList"):
        SaveAsList = Para_dict_used["SaveAsList"]  
        Para_dict_used.pop("SaveAsList")
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

    """ if Para_dict_used.__contains__("Func Electrolyte conductivity [S.m-1]"):
        Para_0.update({
            "Electrolyte conductivity [S.m-1]": 
            eval(Para_dict_used["Func Electrolyte conductivity [S.m-1]"])})  
        Para_dict_used.pop("Func Electrolyte conductivity [S.m-1]")
    if Para_dict_used.__contains__("Func Electrolyte diffusivity [m2.s-1]"):
        Para_0.update({
            "Electrolyte diffusivity [m2.s-1]": 
            eval(Para_dict_used["Func Electrolyte diffusivity [m2.s-1]"])})
        Para_dict_used.pop("Func Electrolyte diffusivity [m2.s-1]")
    if Para_dict_used.__contains__("Func EC transference number"):
        Para_0.update({
            "EC transference number": 
            eval(Para_dict_used["Func EC transference number"])})
        Para_dict_used.pop("Func EC transference number")
    if Para_dict_used.__contains__("Func Cation transference number"):
        Para_0.update({
            "Cation transference number": 
            eval(Para_dict_used["Func Cation transference number"])})
        Para_dict_used.pop("Func Cation transference number")
    if Para_dict_used.__contains__("Func 1 + dlnf/dlnc"):
        Para_0.update({
            "1 + dlnf/dlnc": 
            eval(Para_dict_used["Func 1 + dlnf/dlnc"])})
        Para_dict_used.pop("Func 1 + dlnf/dlnc") """

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
        Total_Cycles,SaveAsList,Temper_i,model_options];
    
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
    index_xlsx, Para_dict_i,   Path_pack , 
    keys_all_AGE,   Exp_AGE_List, exp_index_pack ):

    ModelTimer = pybamm.Timer()
    Para_dict_old = Para_dict_i.copy();
    count_i = int(index_xlsx);
    print('Start Now! Scan %d.' % count_i)  

    # Un-pack data:
    [cycle_no,step_AGE_CD,
        step_AGE_CC,step_AGE_CV, ] = exp_index_pack
    [exp_AGE,exp_AGE,exp_AGE_2, # may be deleted later
        exp_AGE_CD,exp_AGE_CC,exp_AGE_CV] = Exp_AGE_List
    [BasicPath,Target,
        book_name_xlsx,sheet_name_xlsx,] = Path_pack

    ##### Initialise Para_0 and model 
    CyclePack,Para_0 = Para_init(Para_dict_i)
    [Total_Cycles,SaveAsList,
        Temper_i,model_options] = CyclePack
    model_0 = pybamm.lithium_ion.DFN(options=model_options)
    str_model_options = str(model_options)
    # add electrolyte properties as variables, 
    # only when they are not constant
    c_e = model_0.variables["Electrolyte concentration [mol.m-3]"]
    T = model_0.variables["Cell temperature [K]"]
    c_EC = model_0.variables["EC concentration [mol.m-3]"]
    model_0.variables["c(EC) over c(Li+)"] = c_EC / c_e
    if Para_dict_i.__contains__(
        "Electrolyte conductivity [S.m-1]") and isinstance(
            Para_dict_i["Electrolyte conductivity [S.m-1]"], str):
                model_0.variables["Electrolyte conductivity [S.m-1]"] =(
                    Para_0["Electrolyte conductivity [S.m-1]"](c_e,c_EC, T))
    if Para_dict_i.__contains__(
        "Electrolyte diffusivity [m2.s-1]") and isinstance(
            Para_dict_i["Electrolyte diffusivity [m2.s-1]"], str):
                model_0.variables["Electrolyte diffusivity [m2.s-1]"] =(
                    Para_0['Electrolyte diffusivity [m2.s-1]'](c_e,c_EC, T))
    if Para_dict_i.__contains__(
        "Cation transference number") and isinstance(
            Para_dict_i["Cation transference number"], str):
                model_0.variables["Cation transference number"] =(
                    Para_0['Cation transference number'](c_e,c_EC, T))
    if Para_dict_i.__contains__(
        "EC transference number") and isinstance(
            Para_dict_i["EC transference number"], str):
                model_0.variables["EC transference number"] =(
                    Para_0['EC transference number'](c_e,c_EC, T))
    

    # Define experiment- for ageing only, NOT RPT
    str_exp_All = [];Experiment_All = [];
    for SaveAs_i, exp_i in zip(SaveAsList,Exp_AGE_List):
        str_exp_All.append(str(exp_i))
        Experiment_All.append(
            pybamm.Experiment( exp_i * int(SaveAs_i)  ) )


    #####Important: index canot be pre-determined anymore! ######
    var_pts = {
        "x_n": 10,  # negative electrode
        "x_s": 5,  # separator 
        "x_p": 10,  # positive electrode
        "r_n": 30,  # negative particle
        "r_p": 20,  # positive particle
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
                    model_0,
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
                model_old = model_0; sol_old = sol_0    
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
            print(f"Scan No. {index_xlsx} succeed! Saving succeed as well!")

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
        print(f"Scan No. {index_xlsx} degrade too fast!" )

    return Sol_all_i,j,midc_merge



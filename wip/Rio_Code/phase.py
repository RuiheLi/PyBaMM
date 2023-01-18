# Landesfeind electrolyte properties 
def electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs):
    c = c_e / 1000  # mol.m-3 -> mol.l
    p1, p2, p3, p4, p5, p6 = coeffs
    A = p1 * (1 + (T - p2))
    B = 1 + p3 * pybamm.sqrt(c) + p4 * (1 + p5 * pybamm.exp(1000 / T)) * c
    C = 1 + c**4 * (p6 * pybamm.exp(1000 / T))
    sigma_e = A * c * B / C  # mS.cm-1
    return sigma_e / 10

def electrolyte_diffusivity_base_Landesfeind2019(c_e,c_EC, T, coeffs):
    c = c_e / 1000  # mol.m-3 -> mol.l
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

def electrolyte_TDF_EC_EMC_3_7_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array(
        [2.57e1, -4.51e1, -1.77e-1, 1.94, 2.95e-1, 3.08e-4, 2.59e-1, -9.46e-3, -4.54e-4]
    )
    return electrolyte_TDF_base_Landesfeind2019(c_e,c_EC, T, coeffs)


def electrolyte_diffusivity_EC_EMC_3_7_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array([1.01e3, 1.01, -1.56e3, -4.87e2])
    return electrolyte_diffusivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)


def electrolyte_conductivity_EC_EMC_3_7_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array([5.21e-1, 2.28e2, -1.06, 3.53e-1, -3.59e-3, 1.48e-3])
    return electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)

def electrolyte_diffusivity_EC_DMC_1_1_Landesfeind2019(c_e, c_EC,T):
    coeffs = np.array([1.47e3, 1.33, -1.69e3, -5.63e2])
    return electrolyte_diffusivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)

def electrolyte_conductivity_EC_DMC_1_1_Landesfeind2019(c_e,c_EC, T):
    coeffs = np.array([7.98e-1, 2.28e2, -1.22, 5.09e-1, -4e-3, 3.79e-3])
    return electrolyte_conductivity_base_Landesfeind2019(c_e,c_EC, T, coeffs)
    
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

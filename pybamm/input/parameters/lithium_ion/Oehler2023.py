import pybamm
import numpy as np
import os


# Load data in the appropriate format
path, _ = os.path.split(os.path.abspath(__file__))
current_function_Oehler2023_data = pybamm.parameters.process_1D_data(
    "current_function_Oehler2023.csv", path=path
)


def current_function_Oehler2023(sto):
    name, (x, y) = current_function_Oehler2023_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


path, _ = os.path.split(os.path.abspath(__file__))
graphite_ocp_Oehler2023_data = pybamm.parameters.process_1D_data(
    "graphite_ocp_Oehler2023.csv", path=path
)


path, _ = os.path.split(os.path.abspath(__file__))
graphite_diffusivity_Oehler2023_data = pybamm.parameters.process_1D_data(
    "graphite_diffusivity_Oehler2023.csv", path=path
)


def graphite_diffusivity_Oehler2023(sto, T):
    name, (x, y) = graphite_diffusivity_Oehler2023_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


path, _ = os.path.split(os.path.abspath(__file__))
graphite_ocp_Oehler2023_data = pybamm.parameters.process_1D_data(
    "graphite_ocp_Oehler2023.csv", path=path
)


def graphite_ocp_Oehler2023(sto):
    name, (x, y) = graphite_ocp_Oehler2023_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


path, _ = os.path.split(os.path.abspath(__file__))
nmc_diffusivity_Oehler2023_data = pybamm.parameters.process_1D_data(
    "nmc_diffusivity_Oehler2023.csv", path=path
)


def nmc_diffusivity_Oehler2023(sto, T):
    name, (x, y) = nmc_diffusivity_Oehler2023_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


path, _ = os.path.split(os.path.abspath(__file__))
nmc_ocp_Oehler2023_data = pybamm.parameters.process_1D_data(
    "nmc_ocp_Oehler2023.csv", path=path
)


def nmc_ocp_Oehler2023(sto):
    name, (x, y) = nmc_ocp_Oehler2023_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def graphite_electrolyte_exchange_current_density_Ecker2015(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between graphite and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

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

    k_ref = 1.11 * 1e-10

    # multiply by Faraday's constant to get correct units
    m_ref = (
        pybamm.constants.F * k_ref
    )  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    E_r = 53400

    arrhenius = pybamm.exp(-E_r / (pybamm.constants.R * T)) * pybamm.exp(
        E_r / (pybamm.constants.R * 296.15)
    )

    return (
        m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5
    )


def nco_electrolyte_exchange_current_density_Ecker2015(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between NCO and LiPF6 in
    EC:DMC [1, 2, 3].

    References
    ----------
    .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

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

    k_ref = 3.01e-11

    # multiply by Faraday's constant to get correct units
    m_ref = (
        pybamm.constants.F * k_ref
    )  # (A/m2)(m3/mol)**1.5 - includes ref concentrations

    E_r = 4.36e4
    arrhenius = pybamm.exp(-E_r / (pybamm.constants.R * T)) * pybamm.exp(
        E_r / (pybamm.constants.R * 296.15)
    )

    return (
        m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5
    )


def electrolyte_conductivity_base_Landesfeind2019(c_e, T, coeffs):
    """
    Conductivity of LiPF6 in solvent_X as a function of ion concentration and
    temperature. The data comes from [1].

    References
    ----------
    .. [1] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature
    coeffs: :class:`pybamm.Symbol`
        Fitting parameter coefficients

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte conductivity
    """
    c = c_e / 1000  # mol.m-3 -> mol.l
    p1, p2, p3, p4, p5, p6 = coeffs
    A = p1 * (1 + (T - p2))
    B = 1 + p3 * pybamm.sqrt(c) + p4 * (1 + p5 * pybamm.exp(1000 / T)) * c
    C = 1 + c**4 * (p6 * pybamm.exp(1000 / T))
    sigma_e = A * c * B / C  # mS.cm-1

    return sigma_e / 10


def electrolyte_diffusivity_base_Landesfeind2019(c_e, T, coeffs):
    """
    Diffusivity of LiPF6 in solvent_X as a function of ion concentration and
    temperature. The data comes from [1].

    References
    ----------
    .. [1] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature
    coeffs: :class:`pybamm.Symbol`
        Fitting parameter coefficients

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte diffusivity
    """
    c = c_e / 1000  # mol.m-3 -> mol.l
    p1, p2, p3, p4 = coeffs
    A = p1 * pybamm.exp(p2 * c)
    B = pybamm.exp(p3 / T)
    C = pybamm.exp(p4 * c / T)
    D_e = A * B * C * 1e-10  # m2/s

    return D_e


def electrolyte_TDF_base_Landesfeind2019(c_e, T, coeffs):
    """
    Thermodynamic factor (TDF) of LiPF6 in solvent_X as a function of ion concentration
    and temperature. The data comes from [1].

    References
    ----------
    .. [1] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature
    coeffs: :class:`pybamm.Symbol`
        Fitting parameter coefficients

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte thermodynamic factor
    """
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


def electrolyte_t_plus_base_Landesfeind2019(c_e, T, coeffs):
    """
    Transference number of LiPF6 in solvent_X as a function of ion concentration and
    temperature. The data comes from [1].

    References
    ----------
    .. [1] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature
    coeffs: :class:`pybamm.Symbol`
        Fitting parameter coefficients

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte transference number
    """
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


def electrolyte_transference_number_EC_EMC_3_7_Landesfeind2019(c_e, T):
    """
    Transference number of LiPF6 in EC:EMC (3:7 w:w) as a function of ion
    concentration and temperature. The data comes from [1].

    References
    ----------
    .. [1] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte transference number
    """
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

    c_lim = pybamm.Scalar(4000)  # solubility limit

    t_plus = (
        electrolyte_t_plus_base_Landesfeind2019(c_e, T, coeffs) * (c_e <= c_lim) +
        electrolyte_t_plus_base_Landesfeind2019(c_lim, T, coeffs) * (c_e > c_lim)
    )

    return t_plus


def electrolyte_TDF_EC_EMC_3_7_Landesfeind2019(c_e, T):
    """
    Thermodynamic factor (TDF) of LiPF6 in EC:EMC (3:7 w:w) as a function of ion
    concentration and temperature. The data comes from [1].

    References
    ----------
    .. [1] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte thermodynamic factor
    """
    coeffs = np.array(
        [2.57e1, -4.51e1, -1.77e-1, 1.94, 2.95e-1, 3.08e-4, 2.59e-1, -9.46e-3, -4.54e-4]
    )

    c_lim = pybamm.Scalar(4000)  # solubility limit

    TDF = (
        electrolyte_TDF_base_Landesfeind2019(c_e, T, coeffs) * (c_e <= c_lim) +
        electrolyte_TDF_base_Landesfeind2019(c_lim, T, coeffs) * (c_e > c_lim)
    )

    return TDF


def electrolyte_diffusivity_EC_EMC_3_7_Landesfeind2019(c_e, T):
    """
    Diffusivity of LiPF6 in EC:EMC (3:7 w:w) as a function of ion concentration and
    temperature. The data comes from [1].

    References
    ----------
    .. [1] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte diffusivity
    """
    coeffs = np.array([1.01e3, 1.01, -1.56e3, -4.87e2])

    c_lim = pybamm.Scalar(4000)  # solubility limit

    D_e = (
        electrolyte_diffusivity_base_Landesfeind2019(c_e, T, coeffs) * (c_e <= c_lim) +
        electrolyte_diffusivity_base_Landesfeind2019(c_lim, T, coeffs) * (c_e > c_lim)
    )

    return D_e


def electrolyte_conductivity_EC_EMC_3_7_Landesfeind2019(c_e, T):
    """
    Conductivity of LiPF6 in EC:EMC (3:7 w:w) as a function of ion concentration and
    temperature. The data comes from [1].

    References
    ----------
    .. [1] Landesfeind, J. and Gasteiger, H.A., 2019. Temperature and Concentration
    Dependence of the Ionic Transport Properties of Lithium-Ion Battery Electrolytes.
    Journal of The Electrochemical Society, 166(14), pp.A3079-A3097.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte conductivity
    """
    coeffs = np.array([5.21e-1, 2.28e2, -1.06, 3.53e-1, -3.59e-3, 1.48e-3])

    c_lim = pybamm.Scalar(4000)  # solubility limit

    sigma_e = (
        electrolyte_conductivity_base_Landesfeind2019(c_e, T, coeffs) * (c_e <= c_lim) +
        electrolyte_conductivity_base_Landesfeind2019(c_lim, T, coeffs) * (c_e > c_lim)
    )

    return sigma_e


# Call dict via a function to avoid errors when editing in place
def get_parameter_values():
    """
    Parameters for a Kokam SLPB 75106100 cell, from the papers

        Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of a
        lithium-ion battery I. determination of parameters." Journal of the
        Electrochemical Society 162.9 (2015): A1836-A1848.

        Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of a
        lithium-ion battery II. Model validation." Journal of The Electrochemical
        Society 162.9 (2015): A1849-A1857.

    The tab placement parameters are taken from measurements in

        Hales, Alastair, et al. "The cell cooling coefficient: a standard to define heat
        rejection from lithium-ion batteries." Journal of The Electrochemical Society
        166.12 (2019): A2383.

    The thermal material properties are for a 5 Ah power pouch cell by Kokam. The data
    are extracted from

        Zhao, Y., et al. "Modeling the effects of thermal gradients induced by tab and
        surface cooling on lithium ion cell performance."" Journal of The
        Electrochemical Society, 165.13 (2018): A3169-A3178.

    Graphite negative electrode parameters
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    The fits to data for the electrode and electrolyte properties are those provided
    by Dr. Simon O'Kane in the paper:

        Richardson, Giles, et. al. "Generalised single particle models for high-rate
        operation of graded lithium-ion electrodes: Systematic derivation and
        validation." Electrochemica Acta 339 (2020): 135862

    SEI parameters are example parameters for SEI growth from the papers:


        Ramadass, P., Haran, B., Gomadam, P. M., White, R., & Popov, B. N. (2004).
        Development of first principles capacity fade model for Li-ion cells. Journal of
        the Electrochemical Society, 151(2), A196-A203.

        Ploehn, H. J., Ramadass, P., & White, R. E. (2004). Solvent diffusion model for
        aging of lithium-ion battery cells. Journal of The Electrochemical Society,
        151(3), A456-A462.

        Single, F., Latz, A., & Horstmann, B. (2018). Identifying the mechanism of
        continued growth of the solid-electrolyte interphase. ChemSusChem, 11(12),
        1950-1955.

        Safari, M., Morcrette, M., Teyssot, A., & Delacour, C. (2009). Multimodal
        Physics- Based Aging Model for Life Prediction of Li-Ion Batteries. Journal of
        The Electrochemical Society, 156(3),

        Yang, X., Leng, Y., Zhang, G., Ge, S., Wang, C. (2017). Modeling of lithium
        plating induced aging of lithium-ion batteries: Transition from linear to
        nonlinear aging. Journal of Power Sources, 360, 28-40.

    Note: this parameter set does not claim to be representative of the true parameter
    values. Instead these are parameter values that were used to fit SEI models to
    observed experimental data in the referenced papers.
    """

    return {
        "chemistry": "lithium_ion",
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
        "Inner SEI lithium interstitial diffusivity [m2.s-1]": 1e-20,
        "Lithium interstitial reference concentration [mol.m-3]": 15.0,
        "Initial inner SEI thickness [m]": 2.5e-09,
        "Initial outer SEI thickness [m]": 2.5e-09,
        "EC initial concentration in electrolyte [mol.m-3]": 4541.0,
        "EC diffusivity [m2.s-1]": 2e-18,
        "SEI kinetic rate constant [m.s-1]": 1e-12,
        "SEI open-circuit potential [V]": 0.4,
        "SEI growth activation energy [J.mol-1]": 0.0,
        "Negative electrode reaction-driven LAM factor [m3.mol-1]": 0.0,
        "Positive electrode reaction-driven LAM factor [m3.mol-1]": 0.0,
        # cell
        "Negative current collector thickness [m]": 1.4e-05,
        "Negative electrode thickness [m]": 8.22475e-05,
        "Separator thickness [m]": 3.88e-04,
        "Positive electrode thickness [m]": 6.895e-05,
        "Positive current collector thickness [m]": 1.5e-05,
        "Electrode height [m]": 0.05,
        "Electrode width [m]": 0.05,
        "Negative tab width [m]": 0.007,
        "Negative tab centre y-coordinate [m]": 0.0045,
        "Negative tab centre z-coordinate [m]": 0.101,
        "Positive tab width [m]": 0.0069,
        "Positive tab centre y-coordinate [m]": 0.0309,
        "Positive tab centre z-coordinate [m]": 0.101,
        "Cell cooling surface area [m2]": 0.0172,
        "Cell volume [m3]": 1.52e-06,
        "Negative current collector conductivity [S.m-1]": 58411000.0,
        "Positive current collector conductivity [S.m-1]": 36914000.0,
        "Negative current collector density [kg.m-3]": 8933.0,
        "Positive current collector density [kg.m-3]": 2702.0,
        "Negative current collector specific heat capacity [J.kg-1.K-1]": 385.0,
        "Positive current collector specific heat capacity [J.kg-1.K-1]": 903.0,
        "Negative current collector thermal conductivity [W.m-1.K-1]": 398.0,
        "Positive current collector thermal conductivity [W.m-1.K-1]": 238.0,
        "Nominal cell capacity [A.h]": 0.084289862057800647,
        "Typical current [A]": 0.084289862057800647,
        "Current function [A]": current_function_Oehler2023,
        "Contact resistance [Ohm]": 0,
        # negative electrode
        "Negative electrode conductivity [S.m-1]": 100,
        "Maximum concentration in negative electrode [mol.m-3]": 3.563234689881730e+04,
        "Negative electrode diffusivity [m2.s-1]": graphite_diffusivity_Oehler2023,
        "Negative electrode OCP [V]": graphite_ocp_Oehler2023,
        "Negative electrode porosity": 0.369,
        "Negative electrode active material volume fraction": 0.587387813348071,
        "Negative particle radius [m]": 1e-05,
        "Negative electrode Bruggeman coefficient (electrolyte)": 2.8235890158434283,
        "Negative electrode Bruggeman coefficient (electrode)": 0.0,
        "Negative electrode cation signed stoichiometry": -1.0,
        "Negative electrode electrons in reaction": 1.0,
        "Negative electrode exchange-current density [A.m-2]"
        "": graphite_electrolyte_exchange_current_density_Ecker2015,
        "Negative electrode density [kg.m-3]": 1555.0,
        "Negative electrode specific heat capacity [J.kg-1.K-1]": 1437.0,
        "Negative electrode thermal conductivity [W.m-1.K-1]": 1.58,
        "Negative electrode OCP entropic change [V.K-1]": 0.0,
        # positive electrode
        "Positive electrode conductivity [S.m-1]": 68.1,
        "Maximum concentration in positive electrode [mol.m-3]": 4.625677190445059e+04,
        "Positive electrode diffusivity [m2.s-1]": nmc_diffusivity_Oehler2023,
        "Positive electrode OCP [V]": nmc_ocp_Oehler2023,
        "Positive electrode porosity": 0.325,
        "Positive electrode active material volume fraction": 0.605925614624921,
        "Positive particle radius [m]": 1.4e-05,
        "Positive electrode Bruggeman coefficient (electrolyte)": 1.8938165155165517,
        "Positive electrode Bruggeman coefficient (electrode)": 0.0,
        "Positive electrode exchange-current density [A.m-2]"
        "": nco_electrolyte_exchange_current_density_Ecker2015,
        "Positive electrode cation signed stoichiometry": -1.0,
        "Positive electrode electrons in reaction": 1.0,
        "Positive electrode density [kg.m-3]": 2895.0,
        "Positive electrode specific heat capacity [J.kg-1.K-1]": 1270.0,
        "Positive electrode thermal conductivity [W.m-1.K-1]": 1.04,
        "Positive electrode OCP entropic change [V.K-1]": 0.0,
        # separator
        "Separator porosity": 0.8637,
        "Separator Bruggeman coefficient (electrolyte)": 1.5,
        "Separator density [kg.m-3]": 1017.0,
        "Separator specific heat capacity [J.kg-1.K-1]": 1978.0,
        "Separator thermal conductivity [W.m-1.K-1]": 0.34,
        # electrolyte
        "Typical electrolyte concentration [mol.m-3]": 1000.0,
        "Initial concentration in electrolyte [mol.m-3]": 1000.0,
        "Cation transference number"
        "": electrolyte_transference_number_EC_EMC_3_7_Landesfeind2019,
        "Thermodynamic factor": electrolyte_TDF_EC_EMC_3_7_Landesfeind2019,
        "Electrolyte diffusivity [m2.s-1]"
        "": electrolyte_diffusivity_EC_EMC_3_7_Landesfeind2019,
        "Electrolyte conductivity [S.m-1]"
        "": electrolyte_conductivity_EC_EMC_3_7_Landesfeind2019,
        # experiment
        "Reference temperature [K]": 298.15,
        "Negative current collector surface heat transfer coefficient [W.m-2.K-1]"
        "": 10.0,
        "Positive current collector surface heat transfer coefficient [W.m-2.K-1]"
        "": 10.0,
        "Negative tab heat transfer coefficient [W.m-2.K-1]": 10.0,
        "Positive tab heat transfer coefficient [W.m-2.K-1]": 10.0,
        "Edge heat transfer coefficient [W.m-2.K-1]": 10.0,
        "Total heat transfer coefficient [W.m-2.K-1]": 10.0,
        "Ambient temperature [K]": 298.15,
        "Number of electrodes connected in parallel to make a cell": 1.0,
        "Number of cells connected in series to make a battery": 1.0,
        "Lower voltage cut-off [V]": 2.5,
        "Upper voltage cut-off [V]": 4.2,
        "Initial concentration in negative electrode [mol.m-3]": 340,
        "Initial concentration in positive electrode [mol.m-3]": 37112.565805573424,
        "Initial temperature [K]": 298.15,
        # citations
        "citations": [
            "Ecker2015i",
            "Ecker2015ii",
            "Zhao2018",
            "Hales2019",
            "Richardson2020",
        ],
    }

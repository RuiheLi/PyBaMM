from pybamm import exp
def electrolyte_diffusivity_Nyman2008Exp(c_e,c_EC, T):
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

    D_c_e = (
        6 * exp( -1 *(c_e/1000)) 
        * exp(-5/298.15) 
        * exp(-95/298.15*(c_e/1000)) * 1e-10 
    )

    # Nyman et al. (2008) does not provide temperature dependence

    return D_c_e

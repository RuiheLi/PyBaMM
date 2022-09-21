from pybamm import exp, tanh


def graphite_LGM50_ocp_Chen2020(sto):
    """
    LG M50 Graphite open-circuit potential as a function of stochiometry. The fit is
    taken from [1] with an extra term added by Simon O'Kane to capture behaviour in
    the high stoichiometry range.

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
        1.051 * exp(-26.76 * sto)
        + 0.1916
        - 0.05598 * tanh(35.62 * (sto - 0.1356))
        - 0.04483 * tanh(14.64 * (sto - 0.2861))
        - 0.02097 * tanh(26.28 * (sto - 0.6183))
        - 0.02398 * tanh(38.1 * (sto - 1))
    )

    return u_eq

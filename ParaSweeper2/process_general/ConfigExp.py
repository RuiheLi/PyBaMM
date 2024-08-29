"""A module for the ConfigExp class."""


class ConfigExp:
    def __init__(
        self,
        age_T_degC: int | float,
        rpt_T_degC: int | float,
        tot_age_cyc: int,
        age_cyc: int,
    ) -> None:
        """Create a ConfigExp class.

        Args:
            age_T_degC (int | float): ageing temperature in degC
            rpt_T_degC (int | float): reference performance test temperature in degC
            tot_age_cyc (int): total ageing cycles for the whole ageing experiments
            age_cyc (int): number of ageing cycles between two rpt tests.

        For example, if tot_age_cyc = 1000, age_cyc = 100, 
        assuming rpts are performed at both beginning and end of test, the total rpt will be 11.
        """

""" Main function to run one case


"""


# import general packages:
import os
import pybamm as pb

from .load_sweep_task import load_sweep_task
from .process_general.ConfigExp import ConfigExp
from .process_general.ConfigPara import ConfigPara
from .process_general.ConfigModel import ConfigModel
from .process_general.ConfigSol import ConfigSol


path_sweep_task = "wip/Rio_Code/ParaSweeper"

# from layer-1 to layer-2
case = load_sweep_task(path_sweep_task)
config_exp = ConfigExp(case)
config_para = ConfigPara(case)
config_model = ConfigModel(case)
config_sol = ConfigSol(case)

# from layer-2 to layer-3: start to define pybamm objects
rpt_experiments 














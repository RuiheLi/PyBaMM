# Load modules
import pybamm as pb;import pandas as pd;import numpy as np;
import os, json,openpyxl,traceback,multiprocessing,scipy.optimize,sys,gc
import matplotlib.pyplot as plt;
import pickle,imageio,timeit,random,time, signal
from scipy.io import savemat,loadmat;
from pybamm import constants,exp;import matplotlib as mpl


str_path_0 = os.path.abspath(os.path.join(pb.__path__[0],'..'))
str_path_1 = os.path.abspath(
    os.path.join(str_path_0,"wip/Rio_Code/Fun_Upgrade"))
sys.path.append(str_path_1) 

from wip.Rio_Code.ParaSweeper.TaskResult import *
from wip.Rio_Code.ParaSweeper.plot import *
# from wip.Rio_Code.Fun_Upgrade.Run_model import *
from wip.Rio_Code.ParaSweeper.main import *
from wip.Rio_Code.ParaSweeper.post_process import *
from wip.Rio_Code.ParaSweeper.TaskConfig import *
from wip.Rio_Code.ParaSweeper.get_input_para import *

print("import finish?")
# reduce variables that require define outside
On_HPC = False
if On_HPC:
    case_no = 10
    Path_Input = "InputData/" 
    BasicPath_Save=os.getcwd() 
else:
    str_path_0 = os.path.abspath(os.path.join(pb.__path__[0],'..'))
    str_path_1 = os.path.abspath(
        os.path.join(str_path_0,"wip/Rio_Code/Fun_Upgrade"))
    sys.path.append(str_path_1) 
    
    case_no = 10
    Path_Input = os.path.expanduser(
        "~/EnvPBGEM_NC/SimSave/InputData/") # for Linux
    BasicPath_Save =  os.path.expanduser(
        "~/EnvPBGEM_NC/SimSave/P2_R9_Dim")

# define class for all inputs and customized settings:
# start to do configuration:
path_config = Path_config(
    On_HPC=On_HPC, Path_Input=Path_Input,
    BasicPath_Save=BasicPath_Save,
    purpose_i="Full_Exp1235_NC", case_no=case_no,
    rows_per_file= 1,Scan_end_end= 12,
    )
global_config = Global_config(
    On_HPC=On_HPC, Runshort="GEM-2", 
    Plot_Exp=True, Timeout=True, Return_Sol=True, 
    Check_Short_Time=True, R_from_GITT=True,
    fs=13, dpi=100, Re_No=0, Timelimit=int(3600*48),
    Timeout_text = 'I timed out'   )

exp_config = Exp_config()
para_config = Para_config()
model_config = Model_config()
expData_config = ExpData_config(cap_0 = 4.86491)
config = TaskConfig(
    path_config=path_config, global_config=global_config, 
    exp_config=exp_config,para_config=para_config,
    model_config=model_config,expData_config=expData_config)
# Load input file

# normal version
Para_dict_list = load_combinations_from_csv(config.path_config.Para_file) 
pool_no = len(Para_dict_list) # do parallel computing if needed

# Assuming we only do the first one:
config.Get_Para_dict_i(Para_dict_list[0])

# config.para_config.Para_dict_i['Exp No.'] = 8

config.Get_Exp_config()
config.Initialize_exp_text()
config.Get_Exp_Pack()

# Run the model, assuming only one case:
if config.global_config.Re_No == 0:
    task_result, config = Run_One_Case(config) 
elif config.global_config.Re_No > 0:
    pass

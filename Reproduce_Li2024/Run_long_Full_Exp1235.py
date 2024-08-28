""" 
This script is the baseline to run long simuation for 
the paper Li, Ruihe, Kirkaldy, Niall D., Oehler, Fabian, 
Marinescu, Monica, Offer, Gregory J., & O'Kane, Simon E. J. (2023). 
Lithium-ion battery degradation: using degradation mode analysis 
to validate lifetime prediction modelling. arXiv, 
[arXiv:2311.05482v2]. doi:10.48550/arXiv.2311.05482, 
which is currently under review and will be published soon. 

this specific example run cases defined folders "InputData/Full_Exp1235_NC"
to run other cases, just change this to like 
"Full_Exp23_Paper_11_fine"
"SEI_Dry_Exp23_Paper_11_fine"
"""
# Load modules
import pybamm as pb;import pandas as pd;import numpy as np;
import os, json,openpyxl,traceback,multiprocessing,scipy.optimize,sys
import matplotlib.pyplot as plt;
import pickle,imageio,timeit,random,time, signal
from scipy.io import savemat,loadmat;
from pybamm import constants,exp;import matplotlib as mpl

########################     Global settings!!!
rows_per_file = 1;  

Scan_end_end = 12  
purpose_i = "Full_Exp1235_NC" 


# define options:
On_HPC =  False;        Runshort="GEM-2";    Add_Rest = False
Plot_Exp=True;          Timeout=True;     Return_Sol=True;   
Check_Small_Time=True;  R_from_GITT = True
fs = 13; dpi = 100; Re_No =0
Options = [ 
    On_HPC,Runshort,Add_Rest,
    Plot_Exp,Timeout,Return_Sol,
    Check_Small_Time,R_from_GITT,
    dpi,fs]
Timelimit = int(3600*48) # give 48 hours!


if On_HPC:
    i_bundle = int(os.environ["PBS_ARRAY_INDEX"])
else:
    i_bundle = 10; 
Scan_start = (i_bundle-1)*rows_per_file+1;    
Scan_end   = min(Scan_start + rows_per_file-1, Scan_end_end)    
purpose = f"{purpose_i}_Case_{Scan_start}_{Scan_end}"
Target  = f'/{purpose}/'
# interpetation: Simnon suggested, with cracking activation, heat transfer
para_csv = f"Bundle_{i_bundle}.csv"  # name of the random file to get parameters

# Path setting:
if On_HPC:                          # Run on HPC
    Path_csv = f"InputData/{purpose_i}/" 
    Path_Input = "InputData/" 
    BasicPath=os.getcwd() 
    Para_file = Path_csv +  para_csv
else:
    # Add path to system to ensure Fun_P2 can be used
    import sys  
    str_path_0 = os.path.abspath(os.path.join(pb.__path__[0],'..'))
    str_path_1 = os.path.abspath(
        os.path.join(str_path_0,"Reproduce_Li2024"))
    sys.path.append(str_path_1) 
    Path_Input = os.path.expanduser(
        "~/EnvPBGEM_NC/SimSave/InputData/") # for Linux
    BasicPath =  os.path.expanduser(
        "~/EnvPBGEM_NC/SimSave/P2_R9_Dim")
    Para_file = Path_Input+f'{purpose_i}/'+para_csv
# import all functions 
from Fun_NC import * 


# Load input file
Para_dict_list = load_combinations_from_csv(Para_file)
pool_no = len(Para_dict_list) # do parallel computing if needed

midc_merge_all = [];  Sol_RPT_all = [];  Sol_AGE_all = []
Path_List = [BasicPath, Path_Input,Target,purpose] 
# Run the model
if Re_No == 0:
    midc_merge,Sol_RPT,Sol_AGE,DeBug_Lists = Run_P2_Excel (
        Para_dict_list[0], Path_List, 
        Re_No, Timelimit, Options) 
elif Re_No > 0:
    pass

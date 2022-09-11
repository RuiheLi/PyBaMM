import pybamm as pb;import pandas as pd   ;import numpy as np;import os;import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;import matplotlib as mpl; fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'
for k in range(0,1):
    mpl.rcParams["axes.labelsize"] = fs
    mpl.rcParams["axes.titlesize"] = fs
    mpl.rcParams["xtick.labelsize"] =  fs
    mpl.rcParams["ytick.labelsize"] =  fs
    mpl.rcParams["legend.fontsize"] =  fs
    mpl.rcParams['font.sans-serif'] = ['Times New Roman']
    mpl.rcParams['font.serif'] = ['Times New Roman']
    mpl.rcParams['axes.unicode_minus'] = False
import openpyxl
import traceback
import multiprocessing


from Fun_P2_Union_crack import (
    GetScan,
    Run_model_w_dry_out ,
    Run_model_wo_Dry_out,
    Cal_new_con_Update,
    Run_Model_Base_On_Last_Solution,
    Run_Model_Base_On_Last_Solution_RPT,
    write_excel_xlsx,)


# Run model with dry-out:
for i in range(0,1):  
    Total_Cycles = 4; Cycle_bt_RPT = 2; Update_Cycles = 2;  
    CyclePack = [Total_Cycles,Cycle_bt_RPT,Update_Cycles];
    # Key scan parameters:
    Ratio_excess = [2.00,];
    cs_Neg_Init = [29866,];      # default: 29866
    Diff_SEI = [1.7e-19,1.7e-20,1.7e-22];   
    R_SEI = [2E5,];   
    Bulk_Sol_Con =[ 4541.0,];
    D_Li_inSEI = [  1.0e-20,];    # default: 1e-20 
    c_Li_inte_ref = [15,];    # default: 15
    Diff_EC = [2e-20,];       # default 2e-18
    k_SEI   = [1e-16,];       # default 1e-12
    LAMcr_prop=[2.7778e-06,2.7778e-07,2.7778e-08,];# default 2.7778e-07
    Crack_rate=[3.9e-19,3.9e-20,3.9e-22]     # default 3.9e-20
    Couple_SEI_LiP = [1e-6,1e-7,1e-9,]; # default: 1e-6
    k_LiP = [1E-9,1E-10,1E-11];         # default: 1e-10
    Temper = [25,]; 
    MESH_PAR = [120,];

    (TotalScan, DatePack_scan) = GetScan(Ratio_excess,cs_Neg_Init,Diff_SEI,
        R_SEI,Bulk_Sol_Con,D_Li_inSEI,c_Li_inte_ref,
        Diff_EC,k_SEI,LAMcr_prop,Crack_rate,
        Couple_SEI_LiP,k_LiP,Temper,MESH_PAR);
    #BasicPath = 'D:/OneDrive - Imperial College London/SimDataSave/P2R3'; 
    BasicPath=os.getcwd()
    Target  = '/Test_crack_HPC_1/'
    if not os.path.exists(BasicPath + Target):
        os.mkdir(BasicPath + Target);
    book_name_xlsx = 'Solv_SEI_4cycles.xlsx';

    sheet_name_xlsx = 'Results'
    value3 = [
        ["Index", "Ratio_ex","cs_Neg_Init", "Diff_SEI", "R_SEI", 
        "Bulk_Sol_Con","D_Li_inSEI", "c_Li_inte_ref",
        "Diff_EC","k_SEI","LAMcr_prop","Crack_rate",
        "Couple_SEI_LiP","k_LiP","Temper_i","mesh_par",
        "Cap Loss","LLI to LiP",
        "LLI to SEI","LLI to sei-on-cracks",
        "LAM Neg", "LAM Pos", 
        "Vol_Elely_Tot Final", "Vol_Elely_JR Final","Width Final","Error"],
        ]
    write_excel_xlsx(BasicPath + Target+book_name_xlsx, sheet_name_xlsx, value3)

""" for DatePack_scan_jj in DatePack_scan:
    Run_model_w_dry_out(
            CyclePack , DatePack_scan_jj, BasicPath , Target, book_name_xlsx,
        )  """
if __name__ == "__main__":
    pool = multiprocessing.Pool(4)
    processes = [pool.apply_async(
        Run_model_w_dry_out, 
        args=(
            CyclePack , DatePack_scan_jj, BasicPath , Target, book_name_xlsx,
        )
            ) 
            for DatePack_scan_jj in DatePack_scan]
    result = [p.get() for p in processes]

# For GEM-2 NC paper
import csv, random, os, gc
import pybamm as pb;import pandas as pd   ;import numpy as np;import os;
import matplotlib.pyplot as plt;import os;#import imageio;import timeit
from scipy.io import savemat,loadmat;from pybamm import constants,exp,sqrt;
import io
import matplotlib as mpl; import pickle
from multiprocessing import Queue, Process, set_start_method
from queue import Empty
import openpyxl
import traceback
import random;import time, signal
from wip.Rio_Code.Fun_Upgrade.TaskResult import TaskResult
from wip.Rio_Code.Fun_Upgrade.Plot import *
from wip.Rio_Code.Fun_Upgrade.Run_model import *
from wip.Rio_Code.Fun_Upgrade.Post_process import *

###################################################################
#############    New functions from P3 - 221113        ############
###################################################################



def handle_signal(signal_num, frame):
	raise TimeoutError



# Define a function to calculate concentration change, whether electrolyte being squeezed out or added in
def Cal_new_con_Update(Sol,Para):   # subscript r means the reservoir
    # Note: c_EC_r is the initial EC  concentraiton in the reservoir; 
    #       c_e_r  is the initial Li+ concentraiton in the reservoir;
    #############################################################################################################################  
    ###############################           Step-1 Prepare parameters:        #################################################
    #############################################################################################################################  
    L_p   =Para["Positive electrode thickness [m]"]
    L_n   =Para["Negative electrode thickness [m]"]
    L_s   =Para["Separator thickness [m]"]
    L_y   =Para["Electrode width [m]"]   # Also update to change A_cc and therefore Q
    L_z   =Para["Electrode height [m]"]
    c_EC_r_old=Para["Current solvent concentration in the reservoir [mol.m-3]"]   # add two new parameters for paper 2
    c_e_r_old =Para["Current electrolyte concentration in the reservoir [mol.m-3]"]
    # Initial EC concentration in JR
    c_EC_JR_old =Para["Bulk solvent concentration [mol.m-3]"]  # Be careful when multiplying with pore volume to get amount in mole. Because of electrolyte dry out, that will give more amount than real.   
    # LLI due to electrode,Ratio of EC and lithium is 1:1 -> EC amount consumed is LLINegSEI[-1]
    LLINegSEI = (
        Sol["Loss of lithium to SEI [mol]"].entries[-1] 
        - Sol["Loss of lithium to SEI [mol]"].entries[0] )
    LLINegSEIcr = (
        Sol["Loss of lithium to SEI on cracks [mol]"].entries[-1]
        - 
        Sol["Loss of lithium to SEI on cracks [mol]"].entries[0]
        )
    LLINegDeadLiP = (
        Sol["Loss of lithium to dead lithium plating [mol]"].entries[-1] 
        - Sol["Loss of lithium to dead lithium plating [mol]"].entries[0])
    LLINegLiP = (
        Sol["Loss of lithium to lithium plating [mol]"].entries[-1] 
        - Sol["Loss of lithium to lithium plating [mol]"].entries[0])
    cLi_Xavg  = Sol["X-averaged electrolyte concentration [mol.m-3]"].entries[-1] 
    # Pore volume change with time:
    PoreVolNeg_0 = Sol["X-averaged negative electrode porosity"].entries[0]*L_n*L_y*L_z;
    PoreVolSep_0 = Sol["X-averaged separator porosity"].entries[0]*L_s*L_y*L_z;
    PoreVolPos_0 = Sol["X-averaged positive electrode porosity"].entries[0]*L_p*L_y*L_z;
    PoreVolNeg_1 = Sol["X-averaged negative electrode porosity"].entries[-1]*L_n*L_y*L_z;
    PoreVolSep_1 = Sol["X-averaged separator porosity"].entries[-1]*L_s*L_y*L_z;
    PoreVolPos_1 = Sol["X-averaged positive electrode porosity"].entries[-1]*L_p*L_y*L_z;
    #############################################################################################################################  
    ##### Step-2 Determine How much electrolyte is added (from the reservoir to JR) or squeezed out (from JR to reservoir) ######
    #######################   and finish electrolyte mixing     #################################################################  
    Vol_Elely_Tot_old = Para["Current total electrolyte volume in whole cell [m3]"] 
    Vol_Elely_JR_old  = Para["Current total electrolyte volume in jelly roll [m3]"] 
    if Vol_Elely_Tot_old - Vol_Elely_JR_old < 0:
        print('Model error! Electrolyte in JR is larger than in the cell!')
    Vol_Pore_tot_old  = PoreVolNeg_0 + PoreVolSep_0 + PoreVolPos_0    # pore volume at start time of the run
    Vol_Pore_tot_new  = PoreVolNeg_1 + PoreVolSep_1 + PoreVolPos_1    # pore volume at end   time of the run, intrinsic variable 
    Vol_Pore_decrease = Vol_Elely_JR_old  - Vol_Pore_tot_new # WHY Vol_Elely_JR_old not Vol_Pore_tot_old here? Because even the last state of the last solution (n-1) and the first state of the next solution (n) can be slightly different! 
    # EC:lithium:SEI=2:2:1     for SEI=(CH2OCO2Li)2, but because of too many changes are required, change to 2:1:1 for now
    # update 230703: Simon insist we should do: EC:lithium:SEI  =  2:2:1 and change V_SEI accordingly 
    # Because inner and outer SEI partial molar volume is the same, just set one for whole SEI
    VmolSEI   = Para["Outer SEI partial molar volume [m3.mol-1]"] # 9.8e-5,
    VmolLiP   = Para["Lithium metal partial molar volume [m3.mol-1]"] # 1.3e-05
    VmolEC    = Para["EC partial molar volume [m3.mol-1]"]
    #################   KEY EQUATION FOR DRY-OUT MODEL                   #################
    # UPDATE 230525: assume the formation of dead lithium doesn’t consumed EC
    Vol_EC_consumed  =  ( LLINegSEI + LLINegSEIcr   ) * 1 * VmolEC    # Mark: either with 2 or not, will decide how fast electrolyte being consumed!
    Vol_Elely_need   = Vol_EC_consumed - Vol_Pore_decrease
    Vol_SEILiP_increase = 0.5*(
        (LLINegSEI+LLINegSEIcr) * VmolSEI 
        #+ LLINegLiP * VmolLiP
        )    #  volume increase due to SEI+total LiP 
    #################   KEY EQUATION FOR DRY-OUT MODEL                   #################

    Test_V = Vol_SEILiP_increase - Vol_Pore_decrease  #  This value should always be zero, but now not, which becomes the source of error!
    Test_V2= (Vol_Pore_tot_old - Vol_Elely_JR_old) / Vol_Elely_JR_old * 100; # Porosity errors due to first time step
    
    # Start from here, there are serveral variables to be determined:
    #   1) Vol_Elely_add, or Vol_Elely_squeezed; depends on conditions. easier for 'squeezed' condition
    #   2) Vol_Elely_Tot_new, should always equals to Vol_Elely_Tot_old -  Vol_EC_consumed;
    #   3) Vol_Elely_JR_new, for 'add' condition: see old code; for 'squeezed' condition, equals to pore volume in JR
    #   4) Ratio_Dryout, for 'add' condition: see old code; for 'squeezed' condition, equals to Vol_Elely_Tot_new/Vol_Pore_tot_new 
    #   5) Ratio_CeEC_JR and Ratio_CeLi_JR: for 'add' condition: see old code; for 'squeezed' condition, equals to 1 (unchanged)
    #   6) c_e_r_new and c_EC_r_new: for 'add' condition: equals to old ones (unchanged); for 'squeezed' condition, need to carefully calculate     
    #   7) Width_new: for 'add' condition: see old code; for 'squeezed' condition, equals to L_y (unchanged)
    Vol_Elely_Tot_new = Vol_Elely_Tot_old - Vol_EC_consumed;

    
    if Vol_Elely_need < 0:
        print('Electrolyte is being squeezed out, check plated lithium (active and dead)')
        Vol_Elely_squeezed = - Vol_Elely_need;   # Make Vol_Elely_squeezed>0 for simplicity 
        Vol_Elely_add = 0.0;
        Vol_Elely_JR_new = Vol_Pore_tot_new; 
        Ratio_Dryout = Vol_Elely_Tot_new / Vol_Elely_JR_new
        Ratio_CeEC_JR = 1.0;    # the concentrations in JR don't need to change
        Ratio_CeLi_JR = 1.0; 
        SqueezedLiMol   =  Vol_Elely_squeezed*cLi_Xavg;
        SqueezedECMol   =  Vol_Elely_squeezed*c_EC_JR_old;
        Vol_Elely_reservoir_old = Vol_Elely_Tot_old - Vol_Elely_JR_old; 
        LiMol_reservoir_new = Vol_Elely_reservoir_old*c_e_r_old  + SqueezedLiMol;
        ECMol_reservoir_new = Vol_Elely_reservoir_old*c_EC_r_old + SqueezedECMol;
        c_e_r_new= LiMol_reservoir_new / (Vol_Elely_reservoir_old+Vol_Elely_squeezed);
        c_EC_r_new= ECMol_reservoir_new / (Vol_Elely_reservoir_old+Vol_Elely_squeezed);
        Width_new = L_y; 
    else:   # this means Vol_Elely_need >= 0, therefore the folliwing script should works for Vol_Elely_need=0 as well!
        # Important: Calculate added electrolyte based on excessive electrolyte, can be: 1) added as required; 2) added some, but less than required; 3) added 0 
        Vol_Elely_squeezed = 0;  
        if Vol_Elely_Tot_old > Vol_Elely_JR_old:                             # This means Vol_Elely_JR_old = Vol_Pore_tot_old ()
            if Vol_Elely_Tot_old-Vol_Elely_JR_old >= Vol_Elely_need:         # 1) added as required
                Vol_Elely_add     = Vol_Elely_need;  
                Vol_Elely_JR_new  = Vol_Pore_tot_new;  # also equals to 'Vol_Pore_tot_new', or Vol_Pore_tot_old - Vol_Pore_decrease
                Ratio_Dryout = 1.0;
            else:                                                            # 2) added some, but less than required;                                                         
                Vol_Elely_add     = Vol_Elely_Tot_old - Vol_Elely_JR_old;   
                Vol_Elely_JR_new  = Vol_Elely_Tot_new;                       # This means Vol_Elely_JR_new <= Vol_Pore_tot_new
                Ratio_Dryout = Vol_Elely_JR_new/Vol_Pore_tot_new;
        else:                                                                # 3) added 0 
            Vol_Elely_add = 0;
            Vol_Elely_JR_new  = Vol_Elely_Tot_new; 
            Ratio_Dryout = Vol_Elely_JR_new/Vol_Pore_tot_new;

        # Next: start mix electrolyte based on previous defined equation
        # Lithium amount in liquid phase, at initial time point
        TotLi_Elely_JR_Old = Sol["Total lithium in electrolyte [mol]"].entries[-1]; # remember this is only in JR
        # Added lithium and EC amount in the added lithium electrolyte: - 
        AddLiMol   =  Vol_Elely_add*c_e_r_old
        AddECMol   =  Vol_Elely_add*c_EC_r_old
        # Total amount of Li and EC in current electrolyte:
        # Remember Li has two parts, initial and added; EC has three parts, initial, consumed and added
        TotLi_Elely_JR_New   = TotLi_Elely_JR_Old + AddLiMol
        TotECMol_JR   = Vol_Elely_JR_old*c_EC_JR_old - LLINegSEI + AddECMol  # EC:lithium:SEI=2:2:1     for SEI=(CH2OCO2Li)2
        Ratio_CeEC_JR  = TotECMol_JR    /   Vol_Elely_JR_new   / c_EC_JR_old
        Ratio_CeLi_JR  = TotLi_Elely_JR_New    /   TotLi_Elely_JR_Old   /  Ratio_Dryout # Mark, change on 21-11-19
        c_e_r_new  = c_e_r_old;
        c_EC_r_new = c_EC_r_old;
        Width_new   = Ratio_Dryout * L_y;
        
    # Collect the above parameter in DryOut_List to shorten the code  
    DryOut_List   = [
        Vol_EC_consumed, 
        Vol_Elely_need, 
        Test_V, 
        Test_V2, 
        Vol_Elely_add, 
        Vol_Elely_Tot_new, 
        Vol_Elely_JR_new, 
        Vol_Pore_tot_new, 
        Vol_Pore_decrease, 
        c_e_r_new, c_EC_r_new,
        Ratio_Dryout, Ratio_CeEC_JR, 
        Ratio_CeLi_JR,
        Width_new, 
        ]  # 16 in total
    #print('Loss of lithium to negative SEI', LLINegSEI, 'mol') 
    #print('Loss of lithium to dead lithium plating', LLINegDeadLiP, 'mol') 
    #############################################################################################################################  
    ###################       Step-4 Update parameters here        ##############################################################
    #############################################################################################################################
    Para.update(   
        {'Bulk solvent concentration [mol.m-3]':  
         c_EC_JR_old * Ratio_CeEC_JR  })
    Para.update(
        {'EC initial concentration in electrolyte [mol.m-3]':
         c_EC_JR_old * Ratio_CeEC_JR },) 
    Para.update(   
        {'Ratio of Li-ion concentration change in electrolyte consider solvent consumption':  
        Ratio_CeLi_JR }, check_already_exists=False) 
    Para.update(   
        {'Current total electrolyte volume in whole cell [m3]':  Vol_Elely_Tot_new  }, 
        check_already_exists=False)
    Para.update(   
        {'Current total electrolyte volume in jelly roll [m3]':  Vol_Elely_JR_new  }, 
        check_already_exists=False)
    Para.update(   
        {'Ratio of electrolyte dry out in jelly roll':Ratio_Dryout}, 
        check_already_exists=False)
    Para.update(   {'Electrode width [m]':Width_new})    
    Para.update(   
        {'Current solvent concentration in the reservoir [mol.m-3]':c_EC_r_new}, 
        check_already_exists=False)     
    Para.update(   
        {'Current electrolyte concentration in the reservoir [mol.m-3]':c_e_r_new}, 
        check_already_exists=False)             
    return DryOut_List,Para




# Update 2023-10-04 
def Overwrite_Initial_L_SEI_0_Neg_Porosity(Para_0,config):
    """ 
    This is to overwrite the initial negative electrode porosity 
    and initial SEI thickness (inner, outer) to be consistent 
    with the initial capacity loss 
    """
    cap_loss = (
        Para_0["Nominal cell capacity [A.h]"] 
        - config.expData_config.cap_0)
    delta_Q_SEI = cap_loss * 3600
    V_SEI = Para_0["Outer SEI partial molar volume [m3.mol-1]"] 
    # do this when finish updating
    F = 96485.3
    A = Para_0["Electrode width [m]"] * Para_0["Electrode height [m]"]
    z_SEI = Para_0["Ratio of lithium moles to SEI moles"] 
    L_neg = Para_0["Negative electrode thickness [m]"] 
    eps_act_neg = Para_0["Negative electrode active material volume fraction"]
    R_neg =   Para_0["Negative particle radius [m]"]
    l_cr_init = Para_0["Negative electrode initial crack length [m]"]
    w_cr = Para_0["Negative electrode initial crack width [m]"]
    rho_cr = Para_0["Negative electrode number of cracks per unit area [m-2]"]
    a_neg = (3 * eps_act_neg / R_neg)
    roughness = 1 + 2 * l_cr_init * w_cr * rho_cr
    L_SEI_init = delta_Q_SEI  * V_SEI / (
        z_SEI * F * A * L_neg * a_neg * roughness)

    delta_epi = (L_SEI_init ) * roughness * a_neg
    L_inner_init = L_SEI_init / 2
    epi = 0.25 - delta_epi
    # print(L_inner_init,epi)
    # important: update here!
    Para_0["Negative electrode porosity"] = epi
    Para_0["Initial outer SEI thickness [m]"] = L_inner_init + 2.5e-9
    Para_0["Initial inner SEI thickness [m]"] = L_inner_init + 2.5e-9
    print(
        f"""Has Overwritten Initial outer SEI thickness [m] to be 
        {(L_inner_init+2.5e-9):.2e} and Negative electrode 
        porosity to be {epi:.3f} to account for initial 
        capacity loss of {cap_loss:.3f} Ah""")

    return Para_0

# this function is to initialize the para with a known dict
def Para_init(config):
    Para_dict_used = config.para_config.Para_dict_i.copy()
    Para_0 = pb.ParameterValues(Para_dict_used["Para_Set"]  ) # Note: this is a pybamm.object
    Para_dict_used.pop("Para_Set")
    import ast,json
    if Para_dict_used.__contains__("Scan No"):
        Para_dict_used.pop("Scan No")
    if Para_dict_used.__contains__("Exp No."):
        Para_dict_used.pop("Exp No.")

    if Para_dict_used.__contains__("Total ageing cycles"):
        Para_dict_used.pop("Total ageing cycles")
    if Para_dict_used.__contains__("Ageing cycles between RPT"):
        Para_dict_used.pop("Ageing cycles between RPT")
    if Para_dict_used.__contains__("Update cycles for ageing"):
        Para_dict_used.pop("Update cycles for ageing")
    if Para_dict_used.__contains__("Cycles within RPT"):
        RPT_Cycles = Para_dict_used["Cycles within RPT"]  
        Para_dict_used.pop("Cycles within RPT")

    if Para_dict_used.__contains__("Ageing temperature"):
        Age_T_in_K = Para_dict_used["Ageing temperature"]  + 273.15 # update: change to K
        Para_dict_used.pop("Ageing temperature")
    if Para_dict_used.__contains__("RPT temperature"):
        RPT_T_in_K = Para_dict_used["RPT temperature"] + 273.15 # update: change to K 
        Para_dict_used.pop("RPT temperature")
    if Para_dict_used.__contains__("Mesh list"):
        mesh_list = Para_dict_used["Mesh list"]
        mesh_list = json.loads(mesh_list)  
        Para_dict_used.pop("Mesh list")
    if Para_dict_used.__contains__("Exponential mesh stretch"):
        submesh_strech = Para_dict_used["Exponential mesh stretch"]  
        Para_dict_used.pop("Exponential mesh stretch")
    else:
        submesh_strech = "nan";
    if Para_dict_used.__contains__("Model option"):
        model_options = Para_dict_used["Model option"] 
        # careful: need to change str to dict 
        model_options = ast.literal_eval(model_options)
        Para_dict_used.pop("Model option")
    
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
    # 'Outer SEI partial molar volume [m3.mol-1]'
    if Para_dict_used.__contains__(
        "Outer SEI partial molar volume [m3.mol-1]"): 
        Para_0.update(
            {"Outer SEI partial molar volume [m3.mol-1]":
            Para_dict_used["Outer SEI partial molar volume [m3.mol-1]"]})
        Para_0.update(
            {"Inner SEI partial molar volume [m3.mol-1]":
            Para_dict_used["Outer SEI partial molar volume [m3.mol-1]"]})
        
        Para_dict_used.pop("Outer SEI partial molar volume [m3.mol-1]")

    config.exp_config.RPT_Cycles = RPT_Cycles
    config.exp_config.Age_T_in_K = Age_T_in_K
    config.exp_config.RPT_T_in_K = RPT_T_in_K
    config.model_config.mesh_list = mesh_list
    config.model_config.submesh_strech = submesh_strech
    config.model_config.model_options = model_options
    # Mark Ruihe
    for key, value in Para_dict_used.items():
        # risk: will update parameter that doesn't exist, 
        # so need to make sure the name is right 
        if isinstance(value, str):
            Para_0.update({key: eval(value)})
            #Para_dict_used.pop(key)
        else:
            Para_0.update({key: value},check_already_exists=False)

    Para_0 = Overwrite_Initial_L_SEI_0_Neg_Porosity(Para_0,config)
    return config,Para_0

# Input: Para_0
# Output: mdic_dry, Para_0
def Initialize_mdic_dry(Para_0):
    temp_Int_ElelyExces_Ratio =  Para_0[
        "Initial electrolyte excessive amount ratio"] 
    # used to calculate ce_EC_All
    ce_EC_0 = Para_0[
        'EC initial concentration in electrolyte [mol.m-3]'] 
    Para_0.update(
        {
            'EC initial concentration in electrolyte (Unchanged) [mol.m-3]':
            ce_EC_0}, 
        check_already_exists=False)
    if temp_Int_ElelyExces_Ratio < 1:
        Int_ElelyExces_Ratio = -1
        DryOut = "Off"
    else:
        Int_ElelyExces_Ratio = temp_Int_ElelyExces_Ratio;
        DryOut = "On"
    
    if DryOut == "On":  
        mdic_dry = {
            "CeEC_All": [],
            "c_EC_r_new_All": [],
            "c_e_r_new_All": [],
            "Ratio_CeEC_All":[],
            "Ratio_CeLi_All":[],
            "Ratio_Dryout_All":[],

            "Vol_Elely_Tot_All": [],
            "Vol_Elely_JR_All":[],
            "Vol_Pore_tot_All": [],        
            "Vol_EC_consumed_All":[],
            "Vol_Elely_need_All":[],
            "Width_all":[],
            "Vol_Elely_add_All":[],
            "Vol_Pore_decrease_All":[],
            "Test_V_All":[],
            "Test_V2_All":[],
        }
        T_0                  =  Para_0['Initial temperature [K]']
        Porosity_Neg_0       =  Para_0['Negative electrode porosity']  
        Porosity_Pos_0       =  Para_0['Positive electrode porosity']  
        Porosity_Sep_0       =  Para_0['Separator porosity']  
        cs_Neg_Max           =  Para_0["Maximum concentration in negative electrode [mol.m-3]"];
        L_p                  =  Para_0["Positive electrode thickness [m]"]
        L_n                  =  Para_0["Negative electrode thickness [m]"]
        L_s                  =  Para_0["Separator thickness [m]"]
        L_y                  =  Para_0["Electrode width [m]"]
        Para_0.update({'Initial Electrode width [m]':L_y}, check_already_exists=False)
        L_y_0                =  Para_0["Initial Electrode width [m]"]
        L_z                  =  Para_0["Electrode height [m]"]
        Para_0.update({'Initial Electrode height [m]':L_z}, check_already_exists=False)
        L_z_0                =  Para_0["Initial Electrode height [m]"]
        
        Vol_Elely_Tot        = (
            ( L_n*Porosity_Neg_0 +  L_p*Porosity_Pos_0  +  L_s*Porosity_Sep_0  )  
            * L_y_0 * L_z_0 * Int_ElelyExces_Ratio
        ) # Set initial electrolyte amount [L] 
        Vol_Elely_JR         =  (
            ( L_n*Porosity_Neg_0 +  L_p*Porosity_Pos_0  +  L_s*Porosity_Sep_0  )  
            * L_y_0 * L_z_0 )
        Vol_Pore_tot         =  (
            ( L_n*Porosity_Neg_0 +  L_p*Porosity_Pos_0  +  L_s*Porosity_Sep_0  )  
            * L_y_0 * L_z_0 )
        Ratio_CeEC = 1.0 
        Ratio_CeLi = 1.0  
        Ratio_Dryout = 1.0
        Vol_EC_consumed = 0
        Vol_Elely_need = 0
        Vol_Elely_add = 0
        Vol_Pore_decrease = 0
        Test_V2 = 0 
        Test_V = 0
        print('Initial electrolyte amount is ', Vol_Elely_Tot*1e6, 'mL') 
        Para_0.update(
            {'Current total electrolyte volume in jelly roll [m3]':
            Vol_Elely_JR}, check_already_exists=False)
        Para_0.update(
            {'Current total electrolyte volume in whole cell [m3]':
            Vol_Elely_Tot}, check_already_exists=False)   
        
        mdic_dry["Vol_Elely_Tot_All"].append(Vol_Elely_Tot*1e6);            
        mdic_dry["Vol_Elely_JR_All"].append(Vol_Elely_JR*1e6);     
        mdic_dry["Vol_Pore_tot_All"].append(Vol_Pore_tot*1e6);           
        mdic_dry["Ratio_CeEC_All"].append(Ratio_CeEC);                      
        mdic_dry["Ratio_CeLi_All"].append(Ratio_CeLi);             
        mdic_dry["Ratio_Dryout_All"].append(Ratio_Dryout);
        mdic_dry["Vol_EC_consumed_All"].append(Vol_EC_consumed*1e6);        
        mdic_dry["Vol_Elely_need_All"].append(Vol_Elely_need*1e6);     
        mdic_dry["Width_all"].append(L_y_0);
        mdic_dry["Vol_Elely_add_All"].append(Vol_Elely_add*1e6);            
        mdic_dry["Vol_Pore_decrease_All"].append(Vol_Pore_decrease*1e6);
        mdic_dry["Test_V_All"].append(Test_V*1e6); 
        mdic_dry["Test_V2_All"].append(Test_V2*1e6); 
        mdic_dry["c_e_r_new_All"].append(Para_0[
            "Initial concentration in electrolyte [mol.m-3]"]); 
        mdic_dry["c_EC_r_new_All"].append(Para_0[
            'EC initial concentration in electrolyte [mol.m-3]'])
    else:
        mdic_dry ={}
    
    return Para_0, DryOut, mdic_dry
    

    

# judge ageing shape: - NOT Ready yet!
def check_concave_convex(x_values, y_values):
    # Check if the series has at least three points
    if len(x_values) < 3 or len(y_values) < 3:
        return ["None","None""None"]

    # Calculate the slopes between adjacent points based on "x" and "y" values
    slopes = []
    for i in range(1, len(y_values)):
        slope = (y_values[i] - y_values[i - 1]) / (x_values[i] - x_values[i - 1])
        slopes.append(slope)

    # Determine the indices for the start, middle, and end groups
    start_index = 0
    middle_index = len(y_values) // 2 - 1
    end_index = len(y_values) - 3

    # Check if the start group is concave, convex, or linear
    start_group = slopes[:3]
    start_avg = sum(start_group[:2]) / 2
    start_point = start_group[1]
    # TODO we still need to consider a proper index to judge linear
    if abs(start_avg - start_point) <= 0.005 * start_point:
        start_shape = "Linear"
    elif start_avg <= start_point:
        start_shape = "Sub"
    else:
        start_shape = "Super"

    # Check if the middle group is concave, convex, or linear
    middle_group = slopes[middle_index:middle_index + 3]
    middle_avg = sum(middle_group[:2]) / 2
    middle_point = middle_group[1]
    if abs(middle_avg - middle_point) <= 0.005 * middle_point:
        middle_shape = "Linear"
    elif middle_avg <= middle_point:
        middle_shape = "Sub"
    else:
        middle_shape = "Super"

    # Check if the end group is concave, convex, or linear
    end_group = slopes[end_index:]
    end_avg = sum(end_group[:2]) / 2
    end_point = end_group[1]
    if abs(end_avg - end_point) <= 0.005 * end_point:
        end_shape = "Linear"
    elif end_avg <= end_point:
        end_shape = "Sub"
    else:
        end_shape = "Super"

    # Return the shape results
    shape_results = [start_shape, middle_shape, end_shape]
    return shape_results

# Ensure some variables are integral:
def Get_controller_from_Para_dict_i(Para_dict_i,config):
    Para_dict_i["Scan No"] = int(Para_dict_i["Scan No"])
    Para_dict_i["Exp No."] = int(Para_dict_i["Exp No."])
    config.Scan_No = Para_dict_i["Scan No"]
    config.Exp_No = Para_dict_i["Exp No."]
    config.Age_T = Para_dict_i["Ageing temperature"]  
    return Para_dict_i,config

def Get_PyBaMM_Experiment(config):
    # Define experiment: return three pybamm.experiment objects, 
    #    be careful not to include pybamm objects in config
    # uppack:
    
    update  = config.exp_config.update
    exp_AGE_text = config.exp_config.exp_AGE_text
    exp_breakin_text = config.exp_config.exp_breakin_text
    exp_GITT_text = config.exp_config.exp_GITT_text
    exp_refill = config.exp_config.exp_refill
    exp_adjust_before_age = config.exp_config.exp_adjust_before_age
    exp_RPT_text = config.exp_config.exp_RPT_text

    Experiment_Long   = pb.Experiment( 
        exp_AGE_text * update  )
    if config.global_config.R_from_GITT: 
        Experiment_Breakin= pb.Experiment( 
            exp_breakin_text * 1
            + exp_GITT_text*24  
            + exp_refill * 1    # only do refil if have GITT
            + exp_adjust_before_age*1) 
        Experiment_RPT    = pb.Experiment( 
            exp_RPT_text * 1
            + exp_GITT_text*24  
            + exp_refill * 1    # only do refil if have GITT
            + exp_adjust_before_age*1) 
        
    else:   # then get resistance from C/2
        Experiment_Breakin= pb.Experiment( 
            exp_breakin_text * 1
            + exp_adjust_before_age*1) 
        Experiment_RPT    = pb.Experiment( 
            exp_RPT_text * 1
            + exp_adjust_before_age*1) 
    return Experiment_Long,Experiment_Breakin,Experiment_RPT


def Process_DryOut(Sol_Dry_old, Para_0_Dry_old,task_result,config):
    if config.model_config.DryOut == "On":
        DryOut_List,Paraupdate = Cal_new_con_Update(  
            Sol_Dry_old, Para_0_Dry_old)
        task_result.DryOut_List = DryOut_List
    else:
        Paraupdate = Para_0_Dry_old
    return task_result, Paraupdate


def Run_One_Case(config): 
    ##########################################################
    ##############    Part-1: Initialization    ##############
    ##########################################################
    # start counting time
    ModelTimer = pb.Timer() 
    SmallTimer = pb.Timer()
    # unpack frequently used variable
    Scan_No =  config.Scan_No
    Re_No = config.global_config.Re_No

    print(f'Start Now! Scan {Scan_No} Re {Re_No}')  
    config.Create_folders() # create folders

    log_buffer = io.StringIO()

    # Initialize PyBaMM parameter and further update config
    config,Para_0 = Para_init(config) 

    # Get PyBaMM experiment
    (
        Experiment_Long,Experiment_Breakin,
        Experiment_RPT) = Get_PyBaMM_Experiment(config)

    # initialize my_dict for outputs
    task_result = TaskResult()

    # Gey output keys:
    task_result.Get_Output_Keys(config.model_config.model_options)

    # Initialize my own disctionary for pybamm results:
    task_result.Initialize_my_dict()

    # Initialize dry-out model - whether dry-out or not
    Para_0, DryOut, mdic_dry = Initialize_mdic_dry(Para_0)
    config.model_config.DryOut = DryOut
    task_result.Add_Dryout_Dict(mdic_dry)
    task_result.Keys_error = ["Error tot %","Error SOH %","Error LLI %",
        "Error LAM NE %","Error LAM PE %",
        "Error Res %","Error ageT %","Punish"]

    print(
        f'Scan {Scan_No} Re {Re_No}: Spent {SmallTimer.time()} '
        f'on Initialization, DryOut = {DryOut}')
    SmallTimer.reset()

    ##########################################################
    ##############    Part-2: Run model         ##############
    ##########################################################
    ##########################################################
    # initialize the PyBaMM model
    Model_0 = pb.lithium_ion.DFN(options=config.model_config.model_options)

    # run break-in cycle and get results:
    (
        task_result, config, Para_0, Model_0, SmallTimer, Sol_0
        ) = Run_Breakin_lump_with_try_catch(
        task_result,config,Experiment_Breakin, Para_0, 
        Model_0, SmallTimer)
    
    #############################################################
    #######   2-2: Write a big loop to finish the long experiment    
    if task_result.Flag_Breakin == True: 
        k=0
        # Para_All.append(Para_0);Model_All.append(Model_0);Sol_All_i.append(Sol_0); 
        Para_0_Dry_old = Para_0     
        Model_Dry_old = Model_0 
        Sol_Dry_old = Sol_0   
        del Para_0, Model_0, Sol_0 
        gc.collect() 
        while k < config.exp_config.RPT_num:    
            i=0  
            task_result.my_dict_AGE["avg_Age_T"] = []
            while i < config.exp_config.Runs_bt_RPT:
                # process dry-out model - contain non-dry-out option
                task_result, Paraupdate = Process_DryOut(
                    Sol_Dry_old, Para_0_Dry_old,task_result,config)
                
                # Run aging cycle:
                (
                    task_result,config, Para_0_Dry_old, 
                    Model_Dry_old, SmallTimer, Sol_Dry_old
                    ) = Run_Age_lump_with_try_catch(
                    task_result,config,Experiment_Long, Para_0_Dry_old, 
                    Model_Dry_old, Sol_Dry_old, SmallTimer, k, i)
                i += 1;   ### Finish small loop and add 1 to i 
                if task_result.Flag_AGE == False:
                    break

            ##### Run RPT, and also update dry-out 
            #     parameters (otherwise will have problems)

            # process dry-out model - contain non-dry-out option
            task_result, Paraupdate = Process_DryOut(
                Sol_Dry_old, Para_0_Dry_old, task_result, config)  
            # Run RPT cycle
            (
                task_result, config, Para_0_Dry_old, 
                Model_Dry_old, SmallTimer, Sol_Dry_old
                ) = Run_RPT_lump_with_try_catch(
                task_result, config, Experiment_RPT, Paraupdate, 
                Model_Dry_old, Sol_Dry_old, SmallTimer, k, i)
            
            if task_result.Flag_AGE == False or \
                task_result.Flag_partial_AGE == True:
                break

            # update 240730 - increase save frequency here
            ##########################################################

            # Add model LLI, LAM manually 
            task_result = Get_SOH_LLI_LAM(task_result,config)
            print(f"Scan {Scan_No} Re {Re_No}: Getting extra variables "
                f"within {SmallTimer.time()}")
            SmallTimer.reset()

            # Evaluate errors systematically
            task_result = Evaluate_MPE_with_ExpData(task_result,config)

            #########      3-1: Plot cycle,location, Dryout related 
            # set plotting
            font = {'family' : 'DejaVu Sans','size'   : config.global_config.fs}
            mpl.rc('font', **font)

            Plot_Cyc_RPT_4(task_result, config)
            Plot_DMA_Dec(task_result, config)
            if len(task_result.my_dict_AGE["CDend Porosity"])>1:
                Plot_Loc_AGE_4(task_result, config)
                Plot_HalfCell_V(task_result, config)
            if DryOut == "On":
                Plot_Dryout(task_result, config)
            print(f"Scan {Scan_No} Re {Re_No}: Finish all plots"
                f" within {SmallTimer.time()}")
            SmallTimer.reset()

            log_messages = log_buffer.getvalue()
            log_messages = log_messages.strip()
            Write_Dict_to_Excel(task_result, config, log_messages)
            #########      3-2: Save data as .mat or .json
            Fun_Save_for_Reload(task_result,Model_Dry_old,Para_0_Dry_old,config)
            Fun_Save_for_debug(task_result, config)
            # update 231217: save ageing solution if partially succeed in ageing set
            Fun_Save_Partial_Age(task_result, config)
            # update 231217: save ageing solution if partially succeed in ageing set
            print(
                f"Scan {Scan_No} Re {Re_No}: Save Early for RPT {k} "
                f"within {SmallTimer.time()}")
            SmallTimer.reset()

            k += 1 

    
    ############################################################# 
    #########   An extremely bad case: cannot even finish breakin
    if task_result.Flag_Breakin == False: 
        log_messages = log_buffer.getvalue()
        log_messages = log_messages.strip()

        task_result.mpe_all = [np.nan] * 8
        task_result.Pass_Fail = "Fail at break-in cycle"
        # Update 240430 new script: 
        Write_Dict_to_Excel(task_result, config, log_messages,)
        for key in task_result.Keys_error:
            task_result.my_dict_RPT[key] =  np.nan
    
        with open(
            config.path_config.BasicPath_Save 
            + config.path_config.Target
            +"Mats/" + str(config.Scan_No)+ 
            f'-DeBug_List_Breakin_Re_{Re_No}.pkl', 'wb') as file:
            pickle.dump(task_result.DeBug_List_Breakin, file)

        return task_result, config
    ##########################################################
    ##############   Part-3: Post-prosessing    ##############
    ##########################################################
    else:
        ##########################################################
        #########      3-1: Plot cycle,location, Dryout related 
        # update 23-05-25 there is a bug in Cyc_Update_Index, need to slide a bit:
        Cyc_Update_Index = task_result.my_dict_RPT["Cyc_Update_Index"]
        Cyc_Update_Index.insert(0,0); del Cyc_Update_Index[-1]
        task_result.my_dict_RPT["Cyc_Update_Index"] = Cyc_Update_Index

        # set plotting
        font = {'family' : 'DejaVu Sans','size'   : config.global_config.fs}
        mpl.rc('font', **font)

        Plot_Cyc_RPT_4(task_result, config)
        Plot_DMA_Dec(task_result, config)
        if len(task_result.my_dict_AGE["CDend Porosity"])>1:
            Plot_Loc_AGE_4(task_result, config)
            Plot_HalfCell_V(task_result, config)
        if DryOut == "On":
            Plot_Dryout(task_result, config)

        print(f"Scan {Scan_No} Re {Re_No}: Finish all plots"
            f" within {SmallTimer.time()}")
        SmallTimer.reset()

        ##########################################################
        ##########################################################
        #########      3-3: Save summary to excel 
        log_messages = log_buffer.getvalue()
        log_messages = log_messages.strip()
        Write_Dict_to_Excel(task_result, config, log_messages)

        #########      3-2: Save data as .mat or .json
        Fun_Save_for_Reload(task_result,Model_Dry_old,Para_0_Dry_old,config)
        Fun_Save_for_debug(task_result, config)
        # update 231217: save ageing solution if partially succeed in ageing set
        Fun_Save_Partial_Age(task_result, config)
        # update 231217: save ageing solution if partially succeed in ageing set
        print(f"Scan {Scan_No} Re {Re_No}: Try saving within {SmallTimer.time()}")
        SmallTimer.reset()

        print("Succeed doing something in {}".format(ModelTimer.time()))
        print(f'This is the end of No. {Scan_No} scan, Re {Re_No}')
        return task_result, config




""" Main function to run one case"""

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
from wip.Rio_Code.Fun_Upgrade.TaskResult import *
from wip.Rio_Code.Fun_Upgrade.Plot import *
from wip.Rio_Code.Fun_Upgrade.Run_model import *
from wip.Rio_Code.Fun_Upgrade.Post_process import *
from wip.Rio_Code.Fun_Upgrade.TaskConfig import *
from wip.Rio_Code.Fun_Upgrade.Get_input import *

###################################################################
#############    New functions from P3 - 221113        ############
###################################################################



def handle_signal(signal_num, frame):
	raise TimeoutError






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

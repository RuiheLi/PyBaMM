    ##########################################################
    ##############    Part-0: Log of the scripts    ##########
    ##########################################################
    # add 221205: if Timeout=='True', use Patrick's version, disable pool
    #             else, use pool to accelerate 
    # add Return_Sol, on HPC, always set to False, as it is useless, 
    # add 230221: do sol_new['Throughput capacity [A.h]'].entries += sol_old['Throughput capacity [A.h]'].entries 
    #             and for "Throughput energy [W.h]", when use Model.set_initial_conditions_from
    #             this is done inside the two functions Run_Model_Base_On_Last_Solution(_RPT)

    ##########################################################
    ##############    Part-1: Initialization    ##############
    ##########################################################
    font = {'family' : 'DejaVu Sans','size'   : fs}
    mpl.rc('font', **font)
    ModelTimer = pb.Timer()
    Scan_i = int(index_xlsx)
    print('Start Now! Scan %d.' % Scan_i)  
    Sol_RPT = [];  Sol_AGE = [];
    # pb.set_logging_level('INFO') # show more information!
    # set_start_method('fork') # from Patrick

    # Un-pack data:
    [cycle_no,step_AGE_CD,step_AGE_CC,step_AGE_CV,
        step_RPT_CD,step_RPT_RE , step_RPT_CC ] = exp_index_pack;
    [
        exp_AGE_text,       exp_RPT_0p1C_text,  
        exp_RPT_refill_text,exp_RPT_GITT_text, 
        exp_preAge_text     ] = exp_text_list;
    [BasicPath,Target,book_name_xlsx,sheet_name_xlsx,] = Path_pack
    CyclePack,Para_0 = Para_init(Para_dict_i)
    [Total_Cycles,Cycle_bt_RPT,Update_Cycles,RPT_Cycles,
        Temper_i,Temper_RPT,mesh_list,submesh_strech,model_options] = CyclePack;
    [keys_all_RPT,keys_all_AGE] = keys_all
    str_exp_AGE_text  = str(exp_AGE_text);
    str_exp_RPT_text  = str(exp_text_list[1:]);

    # define experiment
    Experiment_Long   = pb.Experiment( exp_AGE_text * Update_Cycles  )  
    # update 230312 - add GITT and compare R_0
    Experiment_RPT    = pb.Experiment( 
        exp_RPT_0p1C_text*1 
        + exp_RPT_refill_text*1
        + exp_RPT_GITT_text*24 
        + exp_preAge_text     ) 
    Experiment_Breakin= Experiment_RPT

    #####  index definition ######################
    Small_Loop =  int(Cycle_bt_RPT/Update_Cycles);   
    SaveTimes = int(Total_Cycles/Cycle_bt_RPT);   

    # initialize my_dict for outputs
    my_dict_RPT = {}
    for keys in keys_all_RPT:
        for key in keys:
            my_dict_RPT[key]=[];
    my_dict_AGE = {}; 
    for keys in keys_all_AGE:
        for key in keys:
            my_dict_AGE[key]=[];
    my_dict_RPT["Cycle_RPT"] = []; 
    my_dict_RPT["Mean_Res_0p1s"] = []; # one numer for one RPT
    my_dict_RPT["Res_0p1s"] = [];      # one list for one RPT
    my_dict_RPT["SOC_GITT"] = [];      # one list for one RPT
    my_dict_AGE["Cycle_AGE"] = []; 
    Cyc_Update_Index     =[]; 
            
    # update 220924: merge DryOut and Int_ElelyExces_Ratio
    temp_Int_ElelyExces_Ratio =  Para_0["Initial electrolyte excessive amount ratio"] 
    ce_EC_0 = Para_0['EC initial concentration in electrolyte [mol.m-3]'] # used to calculate ce_EC_All
    if temp_Int_ElelyExces_Ratio < 1:
        Int_ElelyExces_Ratio = -1;
        DryOut = "Off";
    else:
        Int_ElelyExces_Ratio = temp_Int_ElelyExces_Ratio;
        DryOut = "On";
    print(f"Scan {Scan_i}: DryOut = {DryOut}")
    if DryOut == "On":  
        mdic_dry,Para_0 = Initialize_mdic_dry(Para_0,Int_ElelyExces_Ratio)
    else:
        mdic_dry ={}

    ##########################################################
    ##############    Part-2: Run model         ##############
    ##########################################################
    ##########################################################
    Timeout_text = 'I timed out'
    ##########    2-1: Define model and run break-in cycle
    try:  
        Timelimit = int(3600*2)
        # the following turns on for HPC only!
        if Timeout == True:
            timeout_RPT = TimeoutFunc(
                Run_Breakin, 
                timeout=Timelimit, 
                timeout_val=Timeout_text)
            Result_list_breakin  = timeout_RPT(
                model_options, Experiment_Breakin, 
                Para_0, mesh_list, submesh_strech)
        else:
            Result_list_breakin  = Run_Breakin(
                model_options, Experiment_Breakin, 
                Para_0, mesh_list, submesh_strech)
        [Model_0,Sol_0,Call_Breakin] = Result_list_breakin
        if Return_Sol == True:
            Sol_RPT.append(Sol_0)
        if Call_Breakin.success == False:
            print("Fail due to Experiment error or infeasible")
            1/0
        if Sol_0 == Timeout_text: # to do: distinguish different failure cases
            print("Fail due to Timeout")
            1/0
        if Sol_0 == "Model error or solver error":
            print("Fail due to Model error or solver error")
            1/0
    except ZeroDivisionError as e:
        str_error_Breakin = str(e)
        print(f"Scan {Scan_i}: Fail break-in cycle, need to exit the whole scan now due to {str_error_Breakin} but do not know how!")
        
        Flag_Breakin = False
    else:
        print(f"Scan {Scan_i}: Finish break-in cycle")
        # post-process for break-in cycle - 0.1C only
        my_dict_RPT = GetSol_dict (my_dict_RPT,keys_all_RPT, Sol_0, 
            0, step_RPT_CD , step_RPT_CC , step_RPT_RE, step_AGE_CV   )
        # update 230312 - Get GITT result -need to make sure index goes like this
        cap_full = 5; Index = np.arange(2,26,1) # index = 2:25
        Mean_Res_0p1s,Res_0p1s,SOC = Get_0p1s_R0(Sol_0,Index,cap_full)
        my_dict_RPT["Mean_Res_0p1s"].append(Mean_Res_0p1s) # one numer for one RPT
        my_dict_RPT["Res_0p1s"].append(Res_0p1s)           # one list for one RPT
        my_dict_RPT["SOC_GITT"].append(SOC)                # one list for one RPT
        del Mean_Res_0p1s,Res_0p1s,SOC 
        cycle_count =0; 
        my_dict_RPT["Cycle_RPT"].append(cycle_count)
        Cyc_Update_Index.append(cycle_count);
        Flag_Breakin = True
        
    Flag_AGE = True; str_error_AGE_final = "Empty";   str_error_RPT = "Empty";
    #############################################################
    #######   2-2: Write a big loop to finish the long experiment    
    if Flag_Breakin == True: 
        k=0
        # Para_All.append(Para_0);Model_All.append(Model_0);Sol_All_i.append(Sol_0); 
        Para_0_Dry_old = Para_0;     Model_Dry_old = Model_0  ; Sol_Dry_old = Sol_0;   del Model_0,Sol_0
        while k < SaveTimes:    
            i=0    
            while i < Small_Loop:
                if DryOut == "On":
                    Data_Pack,Paraupdate   = Cal_new_con_Update (  Sol_Dry_old,   Para_0_Dry_old )
                if DryOut == "Off":
                    Paraupdate = Para_0
                # Run aging cycle:
                try:
                    Timelimit = int(3600*2)
                    if Timeout == True:
                        timeout_AGE = TimeoutFunc(
                            Run_Model_Base_On_Last_Solution, 
                            timeout=Timelimit, 
                            timeout_val=Timeout_text)
                        Result_list_AGE = timeout_AGE( 
                            Model_Dry_old  , Sol_Dry_old , Paraupdate ,Experiment_Long, 
                            Update_Cycles,Temper_i,mesh_list,submesh_strech )
                    else:
                        Result_list_AGE = Run_Model_Base_On_Last_Solution( 
                            Model_Dry_old  , Sol_Dry_old , Paraupdate ,Experiment_Long, 
                            Update_Cycles,Temper_i,mesh_list,submesh_strech )
                    [Model_Dry_i, Sol_Dry_i , Call_Age ] = Result_list_AGE
                    if Return_Sol == True:
                        Sol_AGE.append(Sol_Dry_i)
                    #print(f"Temperature for ageing is now: {Temper_i}")  
                    if Call_Age.success == False:
                        print("Fail due to Experiment error or infeasible")
                        str_error_AGE = "Experiment error or infeasible"
                        1/0
                    if Sol_Dry_i == Timeout_text: # fail due to timeout
                        print("Fail due to Timeout")
                        str_error_AGE = "Timeout"
                        1/0
                    if Sol_Dry_i == "Model error or solver error":
                        print("Fail due to Model error or solver error")
                        str_error_AGE = "Model error or solver error"
                        1/0
                except ZeroDivisionError as e:
                    print(f"Scan {Scan_i}: Fail during No.{Cyc_Update_Index[-1]} ageing cycles due to {str_error_AGE}")
                    Flag_AGE = False
                    str_error_AGE_final = str_error_AGE
                    break
                else:
                    Para_0_Dry_old = Paraupdate;       Model_Dry_old = Model_Dry_i;      Sol_Dry_old = Sol_Dry_i;   
                    del Paraupdate,Model_Dry_i,Sol_Dry_i
                    # post-process for first ageing cycle and every -1 ageing cycle
                    if k==0 and i==0:    
                        my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                            0, step_AGE_CD , step_AGE_CC , step_RPT_RE, step_AGE_CV   )     
                        my_dict_AGE["Cycle_AGE"].append(1)
                    my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                        cycle_no, step_AGE_CD , step_AGE_CC , step_RPT_RE, step_AGE_CV   )    
                    cycle_count +=  Update_Cycles; 
                    my_dict_AGE["Cycle_AGE"].append(cycle_count)           
                    Cyc_Update_Index.append(cycle_count)
                    print(f"Scan {Scan_i}: Finish for No.{Cyc_Update_Index[-1]} ageing cycles")
                    if DryOut == "On":
                        mdic_dry = Update_mdic_dry(Data_Pack,mdic_dry)
                    i += 1;   
            # run RPT, and also update parameters (otherwise will have problems)
            if DryOut == "On":
                Data_Pack , Paraupdate  = Cal_new_con_Update (  Sol_Dry_old,   Para_0_Dry_old   )
            if DryOut == "Off":
                Paraupdate = Para_0     
            try:
                Timelimit = int(3600*2)
                if Timeout == True:
                    timeout_RPT = TimeoutFunc(
                        Run_Model_Base_On_Last_Solution_RPT, 
                        timeout=Timelimit, 
                        timeout_val=Timeout_text)
                    Result_list_RPT = timeout_RPT(
                        Model_Dry_old  , Sol_Dry_old ,   
                        Paraupdate,      Experiment_RPT, RPT_Cycles, 
                        Temper_RPT ,mesh_list ,submesh_strech
                    )
                else:
                    Result_list_RPT = Run_Model_Base_On_Last_Solution_RPT(
                        Model_Dry_old  , Sol_Dry_old ,   
                        Paraupdate,      Experiment_RPT, RPT_Cycles, 
                        Temper_RPT ,mesh_list ,submesh_strech
                    )
                [Model_Dry_i, Sol_Dry_i,Call_RPT]  = Result_list_RPT
                if Return_Sol == True:
                    Sol_RPT.append(Sol_Dry_i)
                #print(f"Temperature for RPT is now: {Temper_RPT}")  
                if Call_RPT.success == False:
                    print("Fail due to Experiment error or infeasible")
                    str_error_RPT = "Experiment error or infeasible"
                    1/0 
                if Sol_Dry_i == Timeout_text:
                    print("Fail due to Timeout")
                    str_error_RPT = "Timeout"
                    1/0
                if Sol_Dry_i == "Model error or solver error":
                    print("Fail due to Model error or solver error")
                    str_error_RPT = "Model error or solver error"
                    1/0
            except ZeroDivisionError as e:
                print(f"Scan {Scan_i}: Fail during No.{Cyc_Update_Index[-1]} RPT cycles, due to {str_error_RPT}")
                break
            else:
                my_dict_RPT = GetSol_dict (my_dict_RPT,keys_all_RPT, Sol_Dry_i, 
                    0, step_RPT_CD , step_RPT_CC , step_RPT_RE, step_AGE_CV   )
                my_dict_RPT["Cycle_RPT"].append(cycle_count)
                Cyc_Update_Index.append(cycle_count)
                # update 230312 - Get GITT result -need to make sure index goes like this
                cap_full = 5; Index = np.arange(2,26,1) # index = 2:25
                Mean_Res_0p1s,Res_0p1s,SOC = Get_0p1s_R0(Sol_Dry_i,Index,cap_full)
                my_dict_RPT["Mean_Res_0p1s"].append(Mean_Res_0p1s) # one numer for one RPT
                my_dict_RPT["Res_0p1s"].append(Res_0p1s)           # one list for one RPT
                my_dict_RPT["SOC_GITT"].append(SOC)                # one list for one RPT
                del Mean_Res_0p1s,Res_0p1s,SOC 
                print(f"Scan {Scan_i}: Finish for No.{Cyc_Update_Index[-1]} RPT cycles")
                if DryOut == "On":
                    mdic_dry = Update_mdic_dry(Data_Pack,mdic_dry)
                Para_0_Dry_old = Paraupdate;    Model_Dry_old = Model_Dry_i  ;     Sol_Dry_old = Sol_Dry_i    ;   
                del Paraupdate,Model_Dry_i,Sol_Dry_i
                if Flag_AGE == False:
                    break
            k += 1 
    ############################################################# 
    #########   An extremely bad case: cannot even finish breakin
    if Flag_Breakin == False: 
        value_list_temp = list(Para_dict_i.values())
        values = []
        for value_list_temp_i in value_list_temp:
            values.append(str(value_list_temp_i))
        values.insert(0,str(Scan_i));
        values.insert(1,DryOut);
        values.extend([
            str_exp_AGE_text,
            str_exp_RPT_text,
            "nan","nan",
            "nan","nan", 
            "nan","nan",
            "nan","nan",
            "nan",str_error_Breakin])
        values = [values,]
        print(str_error_Breakin)
        print("Fail in {}".format(ModelTimer.time())) 
        book_name_xlsx_seperate =   str(Scan_i)+ '_' + book_name_xlsx;
        sheet_name_xlsx =  str(Scan_i);
        write_excel_xlsx(
            BasicPath + Target + book_name_xlsx_seperate, 
            sheet_name_xlsx, values)
        
        midc_merge = {**my_dict_RPT, **my_dict_AGE,**mdic_dry}

        return midc_merge,Sol_RPT,Sol_AGE
    ##########################################################
    ##############   Part-3: Post-prosessing    ##############
    ##########################################################
    # Newly add (220517): save plots, not just a single line in excel file:     
    # Newly add (221114): make plotting as functions
    # Niall_data = loadmat( 'Extracted_all_cell.mat'
    else:
        if not os.path.exists(BasicPath + Target + str(Scan_i)):
            os.mkdir(BasicPath + Target + str(Scan_i) );
        dpi= 100;
        # Update 230221 - Add model LLI, LAM manually 
        my_dict_RPT['Throughput capacity [kA.h]'] = (
            np.array(my_dict_RPT['Throughput capacity [A.h]'])/1e3).tolist()
        my_dict_RPT['CDend SOH [%]'] = ((
            np.array(my_dict_RPT["Discharge capacity [A.h]"])
            /my_dict_RPT["Discharge capacity [A.h]"][0])*100).tolist()
        my_dict_RPT["CDend LAM_ne [%]"] = ((1-
            np.array(my_dict_RPT['CDend Negative electrode capacity [A.h]'])
            /my_dict_RPT['CDend Negative electrode capacity [A.h]'][0])*100).tolist()
        my_dict_RPT["CDend LAM_pe [%]"] = ((1-
            np.array(my_dict_RPT['CDend Positive electrode capacity [A.h]'])
            /my_dict_RPT['CDend Positive electrode capacity [A.h]'][0])*100).tolist()
        my_dict_RPT["CDend LLI [%]"] = ((1-
            np.array(my_dict_RPT["CDend Total lithium capacity in particles [A.h]"])
            /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
        if model_options.__contains__("SEI"):
            my_dict_RPT["CDend LLI SEI [%]"] = ((
                np.array(
                    my_dict_RPT["CDend Loss of capacity to SEI [A.h]"]-
                    my_dict_RPT["CDend Loss of capacity to SEI [A.h]"][0]
                    )
                /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
        if model_options.__contains__("SEI on cracks"):
            my_dict_RPT["CDend LLI SEI on cracks [%]"] = ((
                np.array(
                    my_dict_RPT["CDend Loss of capacity to SEI on cracks [A.h]"]-
                    my_dict_RPT["CDend Loss of capacity to SEI on cracks [A.h]"][0]
                    )
                /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
        if model_options.__contains__("lithium plating"):
            my_dict_RPT["CDend LLI lithium plating [%]"] = ((
                np.array(
                    my_dict_RPT["CDend Loss of capacity to lithium plating [A.h]"]-
                    my_dict_RPT["CDend Loss of capacity to lithium plating [A.h]"][0]
                    )
                /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
  
        ##########################################################
        #########      3-1: Plot cycle,location, Dryout related 
        Plot_Cyc_RPT_4(
            my_dict_RPT,
            Exp_Any_AllData,Temp_Cell_Exp, Plot_Exp ,  # =True or False,
            Scan_i,Temper_i,model_options,BasicPath, Target,fs,dpi)
        if len(my_dict_AGE["CDend Porosity"])>1:
            Plot_Loc_AGE_4(my_dict_AGE,Scan_i,model_options,BasicPath, Target,fs,dpi)
        if DryOut == "On":
            Plot_Dryout(Cyc_Update_Index,mdic_dry,ce_EC_0,Scan_i,BasicPath, Target,fs,dpi)
        ##########################################################
        #########      3-2: Save data as .mat 
        my_dict_RPT["Cyc_Update_Index"] = Cyc_Update_Index
        my_dict_RPT["SaveTimes"]    = SaveTimes
        midc_merge = {**my_dict_RPT, **my_dict_AGE,**mdic_dry}
        savemat(BasicPath + Target+"Mats/" + str(Scan_i)+ '-StructDara_for_Mat.mat',midc_merge)  
        ##########################################################
        #########      3-3: Save summary to excel 
        values=Get_Values_Excel(
            model_options,my_dict_RPT,mdic_dry,
            DryOut,Scan_i,Para_dict_i,str_exp_AGE_text,
            str_exp_RPT_text,
            str_error_AGE_final,
            str_error_RPT)
        values = [values,]
        book_name_xlsx_seperate =   str(Scan_i)+ '_' + book_name_xlsx;
        sheet_name_xlsx =  str(Scan_i);
        write_excel_xlsx(
            BasicPath + Target + book_name_xlsx_seperate, 
            sheet_name_xlsx, values)
        print("Succeed doing something in {}".format(ModelTimer.time()))
        print('This is the end of No.', Scan_i, ' scan')
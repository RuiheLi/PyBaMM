"""Module for Get inputs functions"""

import csv, random, os
import numpy as np
import pandas as pd

# read scan files:
def load_combinations_from_csv(Para_file):
    dataframe = pd.read_csv(Para_file)
    parameter_names = dataframe.columns.tolist()
    combinations = dataframe.values.tolist()
    # Get all para
    Para_dict_list = []
    # get all dictionaries
    for combination in combinations:
        input_dict = {}
        for parameter_name,para_value in zip(parameter_names,combination ):
            input_dict[parameter_name] = para_value
        Para_dict_list.append(input_dict)
    print(f"Total scan case is {len(Para_dict_list)}")

    return Para_dict_list



def recursive_scan(mylist,kvs, key_list, acc):
    # 递归终止条件
    if len(key_list) == 0:
        mylist.append(acc.copy())   # copy.deepcopy(acc) 如果value是个list，就要deep copy了
        return mylist
    # 继续递归
    k = key_list[0]
    for v in kvs[k]:
        acc[k] = v
        # print(v)
        recursive_scan(mylist,kvs, key_list[1:], acc)


def save_rows_to_csv(file_path, rows, header):
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)  # Write parameter names as the header row
        writer.writerows(rows)

      
def generate_combinations(Bounds, Num_tot):
    lower_bounds = []
    upper_bounds = []
    for bound in Bounds:
        lower_bounds.append(bound[0])
        upper_bounds.append(bound[1])
    combinations = []
    for _ in range(Num_tot):
        combination = []
        for lower, upper in zip(lower_bounds, upper_bounds):
            value = random.uniform(lower, upper)
            combination.append(value)
        combinations.append(combination)
    return combinations

def save_combinations_to_csv(combinations, parameter_names, filename):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(parameter_names)  # Write parameter names as the first row
        for combination in combinations:
            writer.writerow(combination)

def Get_Scan_files(
        BasicPath_Save,Target_name,model_options,
        parameter_names,para_short_name,
        Pack,   num,
        rows_per_file,Bundle):
    
    import itertools
    para_dict_Same = {
        "Cycles within RPT":1,
        "RPT temperature":25,
        
        "Para_Set": "OKane2023",
        "Model option":model_options,
        "Current solvent concentration in the reservoir [mol.m-3]":4541.0,
        "Current electrolyte concentration in the reservoir [mol.m-3]":1000,
        "Ratio of Li-ion concentration change in " 
        "electrolyte consider solvent consumption":1.0,
        'EC initial concentration in electrolyte [mol.m-3]':4541.0,
        'Typical EC concentration in electrolyte [mol.m-3]':4541.0, 
        }
    unchange_key2 = list(para_dict_Same.keys())
    unchange_val2 = list(para_dict_Same.values())
    short_pack = [lst for lst in Pack if len(lst) > 1]
    selected_indices = [i for i, lst in enumerate(Pack) if len(lst) > 1]
    shortList_para_short_name = [para_short_name[i] for i in selected_indices]
    shortList_para_short_name.insert(0,"No")
    really_change_val =  [
        list(comb) for comb in itertools.product(*short_pack)]

    change_val =  [
        list(comb) for comb in itertools.product(*Pack)]
    combinations = [[i+1,*elem, *unchange_val2] for i,elem in enumerate(change_val)]
    comb_short   = [[i+1,*elem] for i,elem in enumerate(really_change_val)]
    parameter_names = [*parameter_names,*unchange_key2]
    print("Total cases number is",len(combinations))
    if Bundle:
        # Specify the total number of cases
        total_cases = len(combinations)
        # Specify the number of rows per CSV file, rows_per_file
        # Calculate the number of files needed
        num_files = (total_cases - 1) // rows_per_file + 1
        # Create the target folder
        folder_path = os.path.join(BasicPath_Save, "Get_Random_sets", Target_name)
        os.makedirs(folder_path, exist_ok=True)
        # Write data to each CSV file
        for i in range(num_files):
            file_name = f"Bundle_{i+1}.csv"
            file_path = os.path.join(folder_path, file_name)
            start_row = i * rows_per_file
            end_row = min(start_row + rows_per_file, total_cases)
            rows = combinations[start_row:end_row]
            save_rows_to_csv(file_path, rows, parameter_names)
        filename = BasicPath_Save+f"/Get_Random_sets/{Target_name}/"+f'{Target_name}.csv'
        filename_short = BasicPath_Save+f"/Get_Random_sets/{Target_name}/"+f'{Target_name}_s.csv'
        save_combinations_to_csv(combinations, parameter_names, filename)
        save_combinations_to_csv(comb_short, shortList_para_short_name, filename_short)
        print(f"Combinations saved to '{Target_name}.csv' file.") 
        print(f"CSV files created in folder '{Target_name}'.")
    else:
        filename = BasicPath_Save+"/Get_Random_sets/"+f'{Target_name}.csv'
        save_combinations_to_csv(combinations, parameter_names, filename)
        print(f"Combinations saved to '{Target_name}.csv' file.") 
    return len(combinations)

def get_list_from_tuple(d, num):
    if d[1] > d[0]:
        if d[1] > 100 * d[0]:
            result_list = (np.exp(np.linspace(np.log(d[0]), np.log(d[1]), num=num))).tolist()
        else:
            result_list = (np.linspace(d[0], d[1], num=num)).tolist()
    else:
        result_list = []
    return result_list
def Get_Scan_Orth_Latin(
        BasicPath_Save,Target_name,model_options,
        parameter_names,para_short_name,
        Pack, num,
        rows_per_file,Bundle):
    
    import itertools; from pyDOE import lhs
    para_dict_Same = {
        "Cycles within RPT":1,
        "RPT temperature":25,
        #"Mesh list":[5,5,5,60,20],  
        "Para_Set": "OKane2023",
        "Model option":model_options,
        "Current solvent concentration in the reservoir [mol.m-3]":4541.0,
        "Current electrolyte concentration in the reservoir [mol.m-3]":1000,
        "Ratio of Li-ion concentration change in " 
        "electrolyte consider solvent consumption":1.0,
        'EC initial concentration in electrolyte [mol.m-3]':4541.0,
        'Typical EC concentration in electrolyte [mol.m-3]':4541.0, 
        }
    unchange_key2 = list(para_dict_Same.keys())
    unchange_val2 = list(para_dict_Same.values())
    
    Pack_tuple = []; Pack_tuple_index = []
    Pack_list = [];  Pack_list_index  = []
    for i,item in enumerate(Pack):
        if isinstance(item, tuple):
            Pack_tuple.append(item)
            Pack_tuple_index.append(i)
        elif isinstance(item, list):
            Pack_list.append(item)
            Pack_list_index.append(i)
    com_tuple = []; comb_tu_list =[]
    if len(Pack_tuple) > 1:
        for tuple_i in Pack_tuple:
            com_tuple.append( get_list_from_tuple(tuple_i, num) )
        # apply Latin Hypercube:
        #print(com_tuple)
        samples = lhs(len(com_tuple), samples=num)
        for sample in samples:
            combination = []
            for i, candidate_list in enumerate(com_tuple):
                index = int(sample[i] * num)
                combination.append(candidate_list[index])
            comb_tu_list.append(combination)
    else:
        print("error! Pack_tuple must has 2 elements")
    # apply product sampling:
    comb_li_list = [list(comb) for comb in itertools.product(*Pack_list)]
    #print(comb_tu_list)
    #print(comb_li_list)
    Big_Comb = []
    for comb_tu in comb_tu_list:
        for comb_li in comb_li_list:
            big_comb = [0] * (len(comb_tu)+len(comb_li))
            for comb_tu_i,index in zip(comb_tu,Pack_tuple_index):
                big_comb[index] = comb_tu_i
            for comb_li_i,index in zip(comb_li,Pack_list_index):
                big_comb[index] = comb_li_i
            Big_Comb.append(big_comb)
    #print(Big_Comb)
    Big_Comb
    combinations = [[i+1,*elem, *unchange_val2] for i,elem in enumerate(Big_Comb)]

    parameter_names = [*parameter_names,*unchange_key2]
    print("Total cases number is",len(combinations))
    if Bundle:
        # Specify the total number of cases
        total_cases = len(combinations)
        # Specify the number of rows per CSV file, rows_per_file
        # Calculate the number of files needed
        num_files = (total_cases - 1) // rows_per_file + 1
        # Create the target folder
        folder_path = os.path.join(BasicPath_Save, "Get_Random_sets2", Target_name)
        os.makedirs(folder_path, exist_ok=True)
        # Write data to each CSV file
        for i in range(num_files):
            file_name = f"Bundle_{i+1}.csv"
            file_path = os.path.join(folder_path, file_name)
            start_row = i * rows_per_file
            end_row = min(start_row + rows_per_file, total_cases)
            rows = combinations[start_row:end_row]
            save_rows_to_csv(file_path, rows, parameter_names)
        filename = BasicPath_Save+f"/Get_Random_sets2/{Target_name}/"+f'{Target_name}.csv'
        # filename_short = BasicPath_Save+f"/Get_Random_sets2/{Target_name}/"+f'{Target_name}_s.csv'
        save_combinations_to_csv(combinations, parameter_names, filename)
        # save_combinations_to_csv(comb_short, shortList_para_short_name, filename_short)
        print(f"Combinations saved to '{Target_name}.csv' file.") 
        print(f"CSV files created in folder '{Target_name}'.")
    else:
        filename = BasicPath_Save+"/Get_Random_sets2/"+f'{Target_name}.csv'
        save_combinations_to_csv(combinations, parameter_names, filename)
        print(f"Combinations saved to '{Target_name}.csv' file.") 
    return len(combinations)

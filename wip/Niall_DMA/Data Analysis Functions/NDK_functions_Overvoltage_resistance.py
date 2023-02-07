#!/usr/bin/env python
# coding: utf-8

# In[1]:


#importing Pandas etc
import warnings
warnings.simplefilter('ignore', FutureWarning)

import pandas as pd
import numpy as np


# In[2]:

#Define some variables to be used

time = "Time (s)"
V = "Voltage (V)"
I = "Current (mA)"
Qtot = "Charge Step (mA.h)"
T = "Temperature (degC)"
Qdis = "Charge (mA.h)"
SOC = "SOC (%)"
r0 = "R0 (Ohms)"
NS = "Ns changes"
OCV = "OCV (V)"
Rdyn = "Reff-R0 (Ohms)"
Rtotal = "Reff (Ohms)"
OVol = "Over-Voltage (V)"

#Calculate resistance from instantaneous voltage drop on application of a current. Can be used on any type of discharge test (CC or GITT), with the number of output values (rows in output df) equal to th number of 'pulses' (where CC is a single pulse). Alongside the 'R0' resistance, it also outputs the OCV and SoC values immediately prior to the pulse commencing.
def r0_calc_dis(df):
    #creates empty lists for ocv, r0, and soc, which will later be iteratively populated for each pulse
    ocv_vals = []
    r0_vals = []
    soc_vals = []
    
    pulse_indices = df[(df[NS] == 1) & (df[I] < 0)].index
    numPulses = len(pulse_indices)
    print(numPulses)
    counter = 0
    for row in pulse_indices:
        ocv_vals.append(df[V].iloc[(row - 10):(row - 1)].mean())
        r0_vals.append((ocv_vals[counter] - df[V].iloc[row])/((df[I].iloc[row - 1] - df[I].iloc[row])/1000))
        soc_vals.append(df[SOC].iloc[row - 1].round(2))
        counter += 1
        
    #creates a new df called "r0_output" which contains the values of the SOC, R0, and OCV lists created above. The function then returns this new df as the output.
    data = [[i,j,k] for i,j,k in zip(soc_vals, r0_vals, ocv_vals)]
    r0_output_dis = pd.DataFrame(data, columns = ["SOC (%)", "R0 (Ohms)", "OCV (V)"])
    return r0_output_dis

#Same function as 'r0_calc_dis' above, but modified for charging pulses instead of discharge.
def r0_calc_cha(df):
    #creates empty lists for ocv, r0, and soc, which will later be iteratively populated for each pulse
    ocv = []
    r0 = []
    soc = []
    
    pulse_indices = df[(df[NS] == 1) & (df[I] > 0)].index
    numPulses = len(pulse_indices)
    print(numPulses)
    counter = 0
    for row in pulse_indices:
        ocv.append(df[V].iloc[(row - 10):(row - 1)].mean())
        r0.append((ocv[counter] - df[V].iloc[row])/((df[I].iloc[row - 1] - df[I].iloc[row])/1000))
        soc.append(df[SOC].iloc[row - 1].round(2))
        counter += 1
    
    #creates a new df called "r0_output" which contains the values of the SOC, R0, and OCV lists created above. The function then returns this new df as the output.
    data = [[i,j,k] for i,j,k in zip(soc, r0, ocv)]
    r0_output_cha = pd.DataFrame(data, columns = ["SOC (%)", "R0 (Ohms)", "OCV (V)"])
    return r0_output_cha


#Same function as 'r0_calc_dis' above, but modified to work with data from the Maccor. Requires the number of pulses to be input as an argument.
def r0_calc_dis_Maccor(df, numPulses):
    #creates empty lists for ocv, r0, and soc, which will later be iteratively populated for each pulse
    ocv_vals = []
    r0_vals = []
    soc_vals = []
   
    pulse_indices = []
    for pulse in range(1,(numPulses+1)):
        pulse_indices.append(df[df['cycle'] == (pulse)].index.min())
    
    counter = 0
    for row in pulse_indices:
        ocv_vals.append(df[V].iloc[(row - 10):(row - 1)].mean())
        r0_vals.append((ocv_vals[counter] - df[V].iloc[row])/((df[I].iloc[row - 1] - df[I].iloc[row])/1000))
        soc_vals.append(df[SOC].iloc[row - 1].round(2))
        counter += 1
   
    data = [[i,j,k] for i,j,k in zip(soc_vals, r0_vals, ocv_vals)]
    r0_output_dis = pd.DataFrame(data, columns = ["SOC (%)", "R0 (Ohms)", "OCV (V)"])
    return r0_output_dis



#A function which takes an input dataset to be analysed (any type of charge or discharge, signalled by the 'c_type' argument with values of 'c' or 'd') alongside a reference GITT dataset which contains OCV and R0 vs SOC data. The output is a pandas DF with 10 columns, includig overvoltage, total resistance, dynamic resistance, etc.
def over_voltage_new(df_input, GITT_data, c_type='d'):
    # Create copies of dataframes so originals are not affected by operations in this function
    df_input = df_input.copy()
    GITT_data = GITT_data.copy()
    # Selecting only discharge data or only charge
    if c_type == 'c':
        df_input_discharge = df_input.loc[df_input[I] > 0, :].copy()
    else:
        df_input_discharge = df_input.loc[df_input[I] < 0, :].copy()  # The copy method makes a new dataframe instead of a slice of df_input

    # Create new df with linearly spaced SOC list
    df_ocv_r0 = pd.DataFrame(columns=[SOC, r0, V], data=None)
    df_ocv_r0[SOC] = np.linspace(100,0,10001).round(2)

    # Set indexes as SOC columns for each dataframe using the set_index() method
    df_ocv_r0.set_index(SOC, inplace=True)
    GITT_data.loc[:, SOC] = GITT_data.loc[:, SOC].round(2)  # Rounding to 2 decimal places
    GITT_data.set_index(SOC, inplace=True)
    df_input_discharge.loc[:, SOC] = df_input_discharge.loc[:, SOC].round(2)
    df_input_discharge.set_index(SOC, inplace=True)
    df_input_discharge = df_input_discharge.loc[~df_input_discharge.index.duplicated(), :] # Theres a duplicated SOC row at the end, get rid of it

    # Now pandas will automatically do operations where common indices (SOC points) are found
    df_ocv_r0[r0] = GITT_data[r0]
    df_ocv_r0[V] = GITT_data[OCV]

    # Do interpolations
    df_ocv_r0[r0] = df_ocv_r0[r0].astype("float64")
    df_ocv_r0[r0].interpolate(method="polynomial", order=2, inplace=True)
    df_ocv_r0[V] = df_ocv_r0[V].astype("float64")
    df_ocv_r0[V].interpolate(method="polynomial", order=2, inplace=True)

    # Create the columns to go into the final dataframe
    time_out = df_input_discharge.loc[:, time]  # Series for the voltage
    I_out = df_input_discharge.loc[:, I]
    Q_out = df_input_discharge.loc[:, Qdis]
    T_out = df_input_discharge.loc[:, T]
    V_out = df_input_discharge.loc[:, V]  # Series for the voltage
    V_op_out = df_ocv_r0[V] - df_input_discharge[V] # Series for the voltage but contains lots of empty rows where index doesnt match up
    V_op_out.dropna(inplace=True) # Drops all the rows where the SOC index didnt match up
    V_op_out.rename("Over-Voltage (V)", inplace=True)  # Rename here as easier,get error if try concatenating with columns of same name
    OCV_out = V_out + V_op_out
    OCV_out.rename(OCV, inplace=True)
    r_eff_out = V_op_out/-df_input_discharge[I]*1000  # Same as before, lots of empty rows from df_ocv_r0
    r_eff_out.dropna(inplace=True)
    r_eff_out.rename(Rtotal, inplace=True)
    r_eff_corr_out = r_eff_out - df_ocv_r0[r0]  # Same as before, lots of empty rows from df_ocv_r0
    r_eff_corr_out.dropna(inplace=True)
    r_eff_corr_out.rename(Rdyn, inplace=True)
    r0_out = r_eff_out - r_eff_corr_out
    r0_out.rename(r0, inplace=True)

    # Concatenate together the above series into a final dataframe
    df_out = pd.concat([time_out, I_out, Q_out, OCV_out, V_out, V_op_out, r0_out, r_eff_out, r_eff_corr_out, T_out], axis=1)
    # Finally reset the index so SOC becomes a column that can be indexed normally again
    #df_out.reset_index(inplace = True)
    df_out.sort_values(time, inplace=True)
    df_out.reset_index(inplace=True)
    df_out[time] = df_out[time] - df_out[time].iloc[0] #index the time
    return df_out


#Duplicate of "over_voltage_new"? Should maybe be deleted.
def over_voltage_cha(df_input, GITT_data, c_type='d'):
    # Create copies of dataframes so originals are not affected by operations in this function
    df_input = df_input.copy()
    GITT_data = GITT_data.copy()
    # Selecting only discharge data
    df_input_charge = df_input.loc[df_input[I] > 0, :].copy()  # The copy method makes a new dataframe instead of a slice of df_input

    # Create new df with linearly spaced SOC list
    df_ocv_r0 = pd.DataFrame(columns=[SOC, r0, V], data=None)
    df_ocv_r0[SOC] = np.linspace(100,0,10001).round(2)

    # Set indexes as SOC columns for each dataframe using the set_index() method
    df_ocv_r0.set_index(SOC, inplace=True)
    GITT_data.loc[:, SOC] = GITT_data.loc[:, SOC].round(2)  # Rounding to 2 decimal places
    GITT_data.set_index(SOC, inplace=True)
    df_input_charge.loc[:, SOC] = df_input_charge.loc[:, SOC].round(2)
    df_input_charge.set_index(SOC, inplace=True)
    df_input_charge = df_input_charge.loc[~df_input_charge.index.duplicated(), :] # Theres a duplicated SOC row at the end, get rid of it

    # Now pandas will automatically do operations where common indices (SOC points) are found
    df_ocv_r0[r0] = GITT_data[r0]
    df_ocv_r0[V] = GITT_data[OCV]

    # Do interpolations
    df_ocv_r0[r0] = df_ocv_r0[r0].astype("float64")
    df_ocv_r0[r0].interpolate(method="polynomial", order=2, inplace=True)
    df_ocv_r0[V] = df_ocv_r0[V].astype("float64")
    df_ocv_r0[V].interpolate(method="polynomial", order=2, inplace=True)

    # Create the columns to go into the final dataframe
    time_out = df_input_charge.loc[:, time]  # Series for the voltage
    I_out = df_input_charge.loc[:, I]
    Q_out = df_input_charge.loc[:, Qdis]
    T_out = df_input_charge.loc[:, T]
    V_out = df_input_charge.loc[:, V]  # Series for the voltage
    V_op_out = df_ocv_r0[V] - df_input_charge[V] # Series for the voltage but contains lots of empty rows where index doesnt match up
    V_op_out.dropna(inplace=True) # Drops all the rows where the SOC index didnt match up
    V_op_out.rename("Over-Voltage (V)", inplace=True)  # Rename here as easier,get error if try concatenating with columns of same name
    OCV_out = V_out + V_op_out
    OCV_out.rename(OCV, inplace=True)
    r_eff_out = V_op_out/-df_input_charge[I]*1000  # Same as before, lots of empty rows from df_ocv_r0
    r_eff_out.dropna(inplace=True)
    r_eff_out.rename(Rtotal, inplace=True)
    r_eff_corr_out = r_eff_out - df_ocv_r0[r0]  # Same as before, lots of empty rows from df_ocv_r0
    r_eff_corr_out.dropna(inplace=True)
    r_eff_corr_out.rename(Rdyn, inplace=True)
    r0_out = r_eff_out - r_eff_corr_out
    r0_out.rename(r0, inplace=True)

    # Concatenate together the above series into a final dataframe
    df_out = pd.concat([time_out, I_out, Q_out, OCV_out, V_out, V_op_out, r0_out, r_eff_out, r_eff_corr_out, T_out], axis=1)
    # Finally reset the index so SOC becomes a column that can be indexed normally again
    #df_out.reset_index(inplace = True)
    df_out.sort_values(time, inplace=True)
    df_out.reset_index(inplace=True)
    df_out[time] = df_out[time] - df_out[time].iloc[0] #index the time
    return df_out

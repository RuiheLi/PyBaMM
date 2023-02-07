#!/usr/bin/env python
# coding: utf-8

# In[1]:

#importing Pandas etc
import warnings
warnings.simplefilter('ignore', FutureWarning)

import pandas as pd

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
halfCycle = "half cycle"
OCV = "OCV (V)"


#Defines a function which loads the BT-Lab raw mpt files into Panda dataframes based on the built-in "read_csv" function, using some default arguments for standard properties of the BT-Lab datafiles
def read_mpt(filename, rowsToSkip = 0, colsToKeep = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes"]):
    return pd.read_csv(filename, encoding = "shift-jis", skiprows = rowsToSkip, sep = "\t", usecols = colsToKeep)

#Renames the columns of all three dataframes to better titles
#A function is defined to swap the names of the columns in a df
def rename_cols(df, oldNames = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes"], newNames = [time, V, I, Qtot, T, Qdis, NS]):
    return df.rename(columns={i:j for i,j in zip(oldNames,newNames)})

#Creates a new column called "SOC" based on the charge passed during discharge (so doesn't old true for charging periods, where it will display 100% SOC regardless of the real SOC)
#A function is defined which takes the charge passed during discharge away from the theoretical cell capacity ('cap'=5Ah), and normalises to be in percent
def create_SOC(df, cap=5):
    return (100 + (df[Qtot] - df[Qtot].max())/(cap*10))

#Indexes the time column in each dataframe to start from time = 0 by subtracting the first value in the time column from all subsequent values.
#A function is defined to do this
def index_time(df):
    return df[time] - df[time].iloc[0]

#This function splits up a df into new, smaller df's based on the value of the "half cycle" column.
#It takes a df and a "half cycle" number as arguments, and returns a section of the df where "half cycle" matches the input argument.
#Since "half cycle" doesn't make sense as an input (the discharge data for drive cycle 1 is in half cycle 1, but the discharge data for drive cycle 2 is in half cycle 3, etc.), a 2n-1 rule is applied.
def split_cycles(df, cycleNum):
    return df[((df["half cycle"] == (2*cycleNum)-1) & (df[I] < 0))]

#A function which combines some of the useful functions above to make loading/formatting data into a pandas dataframe easier.
def combo_function_fast(filename, numCycles=1, rowsToSkip=0, colsToKeep = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes"], oldNames = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes"], newNames = [time, V, I, Qtot, T, Qdis, NS]):
    df = read_mpt(filename, rowsToSkip, colsToKeep)
    df = rename_cols(df, oldNames, newNames)
    df[SOC] = create_SOC(df)
    df[time] = index_time(df)
    
    return df

#A function which combines some of the useful functions above to make loading/formatting data into a pandas dataframe easier. This one also does OCV & R0 calculation to output accompanying processed datasets. Requires the "r0_calc_dis' function to be imported previously.
def combo_function_full(filename, numCycles=1, rowsToSkip=0, colsToKeep = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes"], oldNames = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes"], newNames = [time, V, I, Qtot, T, Qdis, NS]):
    df = read_mpt(filename, rowsToSkip, colsToKeep)
    df = rename_cols(df)
    df[SOC] = create_SOC(df)
    df[time] = index_time(df)
    df_r0_dis = r0_calc_dis(df)
    df_r0_cha = r0_calc_cha(df)
    
    return df, df_r0_dis, df_r0_cha

#A function which combines some of the useful functions above to make loading/formatting data into a pandas dataframe easier. This one splits a single datafile into several different pandas DFs based on the 'numFiles' argument.
def combo_function_splitter(filename, numCycles=1, rowsToSkip=0, colsToKeep = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes"], oldNames = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes"], newNames = [time, V, I, Qtot, T, Qdis, NS], numFiles=1, splitPoints=[0]):
    df = read_mpt(filename, rowsToSkip, colsToKeep)
    df = rename_cols(df, oldNames, newNames)
    if numFiles == 1:
        df_list = [df]
    else:
        df_list = []
        df['Step No'] = df[NS].cumsum()
        splitPoints.append(df['Step No'].max())
        for fileNum in range(numFiles):
            df_split = df[(df['Step No'] > splitPoints[fileNum]) & (df['Step No'] <= splitPoints[fileNum+1])]
            df_split.reset_index(inplace=True, drop=True)
            df_list.append(df_split)
    for split_df in df_list:
        split_df[SOC] = create_SOC(split_df)
        split_df[time] = index_time(split_df)
    
    return df_list

#A function which combines some of the useful functions above to make loading/formatting data into a pandas dataframe easier. This one splits a single datafile into several different pandas DFs based on the 'numFiles' argument.
def combo_function_cyc_splitter(filename, numCycles=1, rowsToSkip=0, colsToKeep = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes", 'cycle number'], oldNames = ["time/s", "Ecell/V", "I/mA", "(Q-Qo)/mA.h", "Temperature/ｰC", "Q discharge/mA.h", "Ns changes", 'cycle number'], newNames = [time, V, I, Qtot, T, Qdis, NS, 'cycle number'], numFiles=1, splitPoints=[0]):
    df = read_mpt(filename, rowsToSkip, colsToKeep)
    df = rename_cols(df, oldNames, newNames)
    if numFiles == 1:
        df_list = [df]
    else:
        df_list = []
        splitPoints.append(df['cycle number'].max())
        for fileNum in range(numFiles):
            df_split = df[(df['cycle number'] == splitPoints[fileNum])]
            df_split.reset_index(inplace=True, drop=True)
            df_list.append(df_split)
            
    for split_df in df_list:
        split_df[SOC] = create_SOC(split_df)
        split_df[time] = index_time(split_df)
    
    return df_list

#Simple function which takes a raw Biologic .mpt datafile and returns values of the maximum charge passed during discharge and charge, respectively.
def get_capacity(filename, rowsToSkip=0):
    colsToKeep = ["time/s", "Ecell/V", "I/mA", "Q discharge/mA.h", "Q charge/mA.h"]
    oldNames = colsToKeep
    newNames = [time, V, I, "Q Discharge (mAh)", "Q Charge (mAh)"]
    df = read_mpt(filename, rowsToSkip, colsToKeep)
    df = rename_cols(df, oldNames, newNames)
    
    disCap = df["Q Discharge (mAh)"].max()
    chaCap = df["Q Charge (mAh)"].max()
    
    return disCap, chaCap            

#Simple function which takes a pandas DF containing a "Charge (mA.h)" column and returns a value of the maximum discharge capacity measured.
def discharge_capacity(df):
    return df[Qdis].max()


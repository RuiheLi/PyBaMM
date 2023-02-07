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


#A function is defined to swap the names of the columns in a df
def rename_cols(df, oldNames = ['TestTime(s)', 'Amps', 'Volts', 'Amp-hr', 'Cyc#', 'Step'], newNames = [time, I, V, Qdis, 'cycle', NS]):
    return df.rename(columns={i:j for i,j in zip(oldNames,newNames)})


#A function for parsing Maccor data into a formatted pandas DataFrame
def parser_Maccor(filename, cellCapacity, rowsToSkip=2, colsToKeep = ['TestTime(s)', 'Amps', 'Volts', 'Amp-hr', 'Cyc#', 'Step'], newNames = [time, I, V, Qdis, 'cycle', NS]):
    """A function for loading Maccor data (in csv format) into a pandas dataframe and formatting columns etc into the standard
       format used in the other functions contained in this library.
       
       filename: the name of the csv file to be parsed in (including the directory if not the current folder).
       cellCapacity: the nominal cell capacity in A.h (used for SoC calculation).
       rowsToSkip: the number of header rows in the csv file which aren't variable names/data. Default = 2 rows.
       colsToKeep: the variables/column titles you would like to import from the csv file. If values other than default are used, make sure you update the 'newNames' arguement to match.
       newNames: a list of names which will be used as column titles in place of the 'colsToKeep' names. These should ideally consist of the defined variable names at the top of this module (to ensure compatibility with other modules/functions).
       
       The function returns a pandas DataFrame object containing the parsed data.
       """
    
    df = pd.read_csv(filename, skiprows=rowsToSkip, usecols=colsToKeep)
    if len(colsToKeep) == len(newNames):
        df = rename_cols(df, oldNames=colsToKeep, newNames=newNames)
    else:
        print('Error: the old and new column name lists are unequal lengths; columns have not been renamed')
    # Create an SOC column
    df[Qdis] = df[Qdis]*1000
    df[Qtot] = df[Qdis]*np.sign(df[I])
    raw_dats = df.iloc[1:].copy()
    raw_dats.reset_index(inplace=True, drop=True)
    raw_dats.loc[raw_dats.shape[0]] = [0]*7
    df['dQ'] = -(df[Qtot] - raw_dats[Qtot])
    # I can't remember why I made the next two rows... some kind of noise in the data perhaps?
    df['dQ'][df['dQ'] < -1] = 0
    df['dQ'][df['dQ'] > 1] = 0
    df[Qtot] = np.cumsum(df['dQ'])
    df[SOC] = 100 + (df[Qtot] - df[Qtot].max())/(cellCapacity*10)
    
    return df

#A function for parsing Basytec data into a formatted pandas DataFrame
def parser_Basytec(filename, cellCapacity, rowsToSkip=12, colsToKeep = ['~Time[s]', 'I[A]', 'U[V]', 'Ah[Ah]', 'Ah-Cyc-Discharge', 'Line'], newNames = [time, I, V, Qtot, Qdis, NS]):
    """A function for loading Basytec data (in csv format) into a pandas dataframe and formatting columns etc into the standard
       format used in the other functions contained in this library.
       
       filename: the name of the csv file to be parsed in (including the directory if not the current folder).
       cellCapacity: the nominal cell capacity in A.h (used for SoC calculation).
       rowsToSkip: the number of header rows in the csv file which aren't variable names/data. Default = 2 rows.
       colsToKeep: the variables/column titles you would like to import from the csv file. If values other than default are used, make sure you update the 'newNames' arguement to match.
       newNames: a list of names which will be used as column titles in place of the 'colsToKeep' names. These should ideally consist of the defined variable names at the top of this module (to ensure compatibility with other modules/functions).
       
       The function returns a pandas DataFrame object containing the parsed data.
       """
    
    df = pd.read_csv(filename, encoding='ansi', skiprows=rowsToSkip, usecols=colsToKeep)
    if len(colsToKeep) == len(newNames):
        df = rename_cols(df, oldNames=colsToKeep, newNames=newNames)
    else:
        print('Error: the old and new column name lists are unequal lengths; columns have not been renamed')
    # Create an SOC column and change units of charge and current
    df[Qdis] = df[Qdis]*1000
    df[Qtot] = df[Qtot]*1000
    df[I] = df[I]*1000
    df[SOC] = 100 + (df[Qtot] - df[Qtot].max())/(cellCapacity*10)
    
    return df

#Indexes the time column in each dataframe to start from time = 0 by subtracting the first value in the time column from all subsequent values.
#A function is defined to do this
def index_time(df):
    return df[time] - df[time].iloc[0]


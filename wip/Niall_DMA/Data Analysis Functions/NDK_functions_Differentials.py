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
T = "Temperature (degC)"
SOC = "SOC (%)"
r0 = "R0 (Ohms)"
dQ_dV = "dSOC/dV (%/V)"
OCV = "OCV (V)"
Rdyn = "Reff-R0 (Ohms)"



#A function for calculating dV/dQ data from a pandas DF containing "SOC", "I", and "V" columns. The function calculates dV/dQ for both discharge and charge, outputting two separate DFs. Each DF contains "SOC", "V", and "dV/dQ" columns. The function uses a finite-difference method, with default dQ (or actually dSOC) values of 1% (i.e. 100 intervals).
def dVdQ(df_data, dQ_range=1, Qtotal=100):
    
    df_dis = df_data[(df_data[I]<0)]
    df_cha = df_data[((df_data[I]>0) & (df_data[time]>df_dis[time].iloc[0]))]
    
    Segments = int(Qtotal/dQ_range)
    dVdQ_dis = []
    soc_dis = []
    dVdQ_cha = []
    soc_cha = []
    V_dis = []
    V_cha = []
    
    for i in range(Segments):
        sub_df_dis = df_dis[(df_dis[SOC] < (Qtotal-dQ_range*i)) & (df_dis[SOC] > (Qtotal-dQ_range*(i+1)))]
        #sub_df.set_index((range(len(sub_df))), inplace=True)
        Qcol = sub_df_dis[SOC]
        Vcol = sub_df_dis[V]
        if len(Qcol) < 1:
            continue
        dQ = -(Qcol.iloc[-1] - Qcol.iloc[0])
        dV = Vcol.iloc[-1] - Vcol.iloc[0]
        dVdQ_dis.append(dV/dQ)
        soc_dis.append(Qtotal-dQ_range*i)
        V_dis.append(Vcol.iloc[0])
    
    dVdQ_dis_data = [[k,l,m] for k, l, m in zip(soc_dis, V_dis, dVdQ_dis)]
    dVdQ_dis_data_output = pd.DataFrame(dVdQ_dis_data, columns=[SOC, V, "dV/dSOC (V/%)"])
    
    for i in range(Segments):
        sub_df_cha = df_cha[(df_cha[SOC] > (dQ_range*i)) & (df_cha[SOC] < (dQ_range*(i+1)))]
        #sub_df.set_index((range(len(sub_df))), inplace=True)
        Qcol = sub_df_cha[SOC]
        Vcol = sub_df_cha[V]
        if len(Qcol) < 1:
            continue
        dQ = (Qcol.iloc[-1] - Qcol.iloc[0])
        dV = Vcol.iloc[-1] - Vcol.iloc[0]
        dVdQ_cha.append(dV/dQ)
        soc_cha.append(dQ_range*i)
        V_cha.append(Vcol.iloc[0])
    
    dVdQ_cha_data = [[k,l,m] for k, l, m in zip(soc_cha, V_cha, dVdQ_cha)]
    dVdQ_cha_data_output = pd.DataFrame(dVdQ_cha_data, columns=[SOC, V, "dV/dSOC (V/%)"])
    
    return dVdQ_dis_data_output, dVdQ_cha_data_output


#A function for calculating dQ/dV data from a pandas DF containing "SOC", "I", and "V" columns. The function can calculate dQ/dV for discharge ('d'), charge ('c') or both ('b'), outputting either a single DF or a list of 2 DFs. Each DF contains "SOC", "V", and "dQ/dV" columns. The function uses a finite-difference method, with default dV values of 5 mV with a total range of 1700 mV (i.e. 4.2 to 2.5 V). 
def dQdV(df_data, dV_range=0.005, Vtotal=1.700, I_type = 'd'):
    #The function is more complicated due to the possibility of having discharge, charge, or both. The opertions themselves are very simple though.
    if I_type in ['d','b']:
        df_dis = df_data[(df_data[I]<0)].copy()
        Vmax_dis = df_dis[V].max()
        dQdV_dis = []
        V_dis = []
        soc_dis = []
    if I_type in ['c']:
        df_cha = df_data[df_data[I]>0].copy()
        Vmin_cha = df_cha[V].min()
        dQdV_cha = []
        V_cha = []
        soc_cha = []
    if I_type in ['b']:
        df_cha = df_data[((df_data[I]>0) & (df_data[time]>df_dis[time].iloc[0]))].copy()
        Vmin_cha = df_cha[V].min()
        dQdV_cha = []
        V_cha = []
        soc_cha = []
    
    #Determine the total number of 'segments' (datapoints in the output) based on the optional arguments of the voltage range ('Vtotal') and voltage interval ('dV_range').
    Segments = int(Vtotal/dV_range)
    output_data = []  
    
    #Again this looks more complicated due to the options of difference charge/discharge datasets, but the function just cycles through the number of segments, calculating the 'dQ' and 'dV' for each, alongside the corresponding voltage and SOC values at the beginning of each segment.
    for i in range(Segments):
        
        if I_type in ['d','b']:
            sub_df_dis = df_dis[(df_dis[V] < (Vmax_dis-dV_range*i)) & (df_dis[V] > (Vmax_dis-dV_range*(i+1)))]
            Qcol_dis = sub_df_dis[SOC]
            Vcol_dis = sub_df_dis[V]
        
        if I_type in ['c','b']:
            sub_df_cha = df_cha[((df_cha[V] > (Vmin_cha+(dV_range*i))) & (df_cha[V] < (Vmin_cha+(dV_range*(i+1)))))]
            Qcol_cha = sub_df_cha[SOC]
            Vcol_cha = sub_df_cha[V]
        
        if I_type in ['b']:
            if len(Vcol_dis) < 1:
                if len(Vcol_cha) < 1:
                    continue
                dQ_cha = (Qcol_cha.iloc[-1] - Qcol_cha.iloc[0])
                dV_cha = Vcol_cha.iloc[-1] - Vcol_cha.iloc[0]
                dQdV_cha.append(dQ_cha/dV_cha)
                V_cha.append(Vmin_cha+(dV_range*i))
                soc_cha.append(Qcol_cha.iloc[0])
                continue
        
        if I_type in ['d','b']:
            if len(Vcol_dis) < 1:
                continue
            dQ_dis = -(Qcol_dis.iloc[-1] - Qcol_dis.iloc[0])
            dV_dis = Vcol_dis.iloc[-1] - Vcol_dis.iloc[0]
            dQdV_dis.append(dQ_dis/dV_dis)
            V_dis.append(Vmax_dis-dV_range*i)
            soc_dis.append(Qcol_dis.iloc[0])
        
        if I_type in ['c','b']:
            if len(Vcol_cha) < 1:
                continue
            dQ_cha = (Qcol_cha.iloc[-1] - Qcol_cha.iloc[0])
            dV_cha = Vcol_cha.iloc[-1] - Vcol_cha.iloc[0]
            dQdV_cha.append(dQ_cha/dV_cha)
            V_cha.append(Vmin_cha+(dV_range*i))
            soc_cha.append(Qcol_cha.iloc[0])
    
    if I_type in ['d','b']:
        dQdV_dis_data = [[k,l,m] for k, l, m in zip(soc_dis, V_dis, dQdV_dis)]
        dQdV_dis_data_output = pd.DataFrame(dQdV_dis_data, columns=[SOC, V, "dSOC/dV (%/V)"])
        output_data.append(dQdV_dis_data_output)
    
    if I_type in ['c','b']:
        dQdV_cha_data = [[k,l,m] for k, l, m in zip(soc_cha, V_cha, dQdV_cha)]
        dQdV_cha_data_output = pd.DataFrame(dQdV_cha_data, columns=[SOC, V, "dSOC/dV (%/V)"])
        output_data.append(dQdV_cha_data_output)
    
    if I_type in ['d','c']:
        return output_data[0]
    else:
        return output_data


#Special dQdV function where the voltage range is taken from the data. Only works for specific datasets.    
def dQdV_special(df_data, dV_range=0.005):
    
    df_dis = df_data.copy()
    df_dis[SOC] = df_dis[SOC]*100
    Vmax_dis = df_dis['OCV'].max()
    Vmin_dis = df_dis['OCV'].min()
    dQdV_dis = []
    V_dis = []
    soc_dis = []
    Vtotal = Vmax_dis - Vmin_dis
    
    Segments = int(Vtotal/dV_range)
    
    for i in range(Segments):
        
        sub_df_dis = df_dis[(df_dis['OCV'] < (Vmax_dis-dV_range*i)) & (df_dis['OCV'] > (Vmax_dis-dV_range*(i+1)))]
        Qcol_dis = sub_df_dis[SOC]
        Vcol_dis = sub_df_dis['OCV']
              
        if len(Vcol_dis) < 1:
            continue
        dQ_dis = -(Qcol_dis.iloc[-1] - Qcol_dis.iloc[0])
        dV_dis = Vcol_dis.iloc[-1] - Vcol_dis.iloc[0]
        dQdV_dis.append(dQ_dis/dV_dis)
        V_dis.append(Vmax_dis-dV_range*i)
        soc_dis.append(Qcol_dis.iloc[0])

    dQdV_dis_data = [[k,l,m] for k, l, m in zip(soc_dis, V_dis, dQdV_dis)]
    dQdV_dis_data_output = pd.DataFrame(dQdV_dis_data, columns=[SOC, V, "dSOC/dV (%/V)"])
   
    return dQdV_dis_data_output


#A function for calculating DTV (dT/dV) spectra from two datasets: a temperature profile from e.g. picolog, and a voltage profile from e.g. Biologic. Datasets must be time-indexed in advance. This is done by calculating dV/dt and dT/dt using a finite difference method, similar to the "dQdV" function. Outputs a pandas DF with the DTV values as well as the constituent dV/dt and dT/dt etc.
def DTV(df_T_data, df_V_data, dt_range=30):
    
    t_total = df_T_data[time].max() #Finds how long the T data was recorded for
    #print(t_total)
    
    Segments = int(t_total/dt_range) #How many segments to split the data into, defined by the dt_range argument
    dTdt = [] #dT/dt list
    soc_lst = [] #SOC list
    dVdt = [] #dV/dt list
    V_lst = [] #V list
    time_lst = [] #time list
    DTV_lst = [] #dT/dV list (DTV)
    
    #Loop through the number of segments that the data will be split into
    for i in range(Segments):
        #Create sub-df's from the temperature and voltage data
        sub_df_T = df_T_data[(df_T_data[time] > (dt_range*i)) & (df_T_data[time] < (dt_range*(i+1)))]
        sub_df_V = df_V_data[(df_V_data[time] > (dt_range*i)) & (df_V_data[time] < (dt_range*(i+1)))]
        #sub_df.set_index((range(len(sub_df))), inplace=True)
        Qcol = sub_df_V[SOC] #SOC from the V sub-df
        Vcol = sub_df_V[V] #V from the V sub-df
        Tcol = (sub_df_T[T]+273.15) #T from the T sub-df (converted to Kelvin scale)
        timeTcol = sub_df_T[time] #time from T sub-df
        timeVcol = sub_df_V[time] #time from V sub-df
        #print(len(Qcol))
        if len(Qcol) < 1: #If no data is present in the sub-df, skip to the next segment (next loop iteration)
            continue
        dT = Tcol.iloc[-1] - Tcol.iloc[0] #Calculate dT from start to end of segment
        dV = Vcol.iloc[-1] - Vcol.iloc[0] #Calculate dV from start to end of segment
        dt_T = timeTcol.iloc[-1] - timeTcol.iloc[0] #Calculate dt from start to end of segment (for T data)
        dt_V = timeVcol.iloc[-1] - timeVcol.iloc[0] #Calculate dT from start to end of segment (for V data)
        dTdt.append(dT/dt_T) #Calculate dT/dt and add to list
        dVdt.append(dV/dt_V) #Calculate dV/dt and add to list
        soc_lst.append(Qcol.iloc[0]) #Add the current SOC to a list
        V_lst.append(Vcol.iloc[0]) #Add the current V to a list
        time_lst.append(dt_range*i) #Add the time to a list
        #print(dT)
    
    #Calculate the dT/dV from the dT/dt and dV/dt data
    for n in range(len(time_lst)):
        DTV_lst.append(dTdt[n]/dVdt[n]) 
    
    #Form a new df with the calculated DTV data and use as output
    DTV_data = [[a,b,c,d,e,f] for a,b,c,d,e,f in zip(time_lst, V_lst, soc_lst, dTdt, dVdt, DTV_lst)]
    DTV_data_output = pd.DataFrame(DTV_data, columns=[time, V, SOC, "dT/dt (K/s)", "dV/dt(V/s)", "dT/dV (K/V)"])
    
    return DTV_data_output

#A function for calculating DTV (dT/dV) spectra from a combined pandas DF which contains voltage and temperature data. This function also requires a 'c_rate' argument for correcting for different values. The temperature column must be named 'Temperature av (degC)'. This is done by calculating dV/dt and dT/dt using a finite difference method, similar to the "dQdV" function. Outputs a pandas DF with the DTV values as well as the constituent dV/dt and dT/dt etc.
def DTV_new(df_data, c_rate, dt=20, sample_rate=1):
    
    T_av = 'Temperature av (degC)'
    df = df_data.copy()
    df.reset_index(inplace=True, drop=True)
    
    
    dt_interval = int(dt/(sample_rate*c_rate)) #This modifies the dt_interval based on the 'c_rate' argument
    total_t = df[time].max() - df[time].min()
    segments = int(total_t/dt_interval)
    
    dTdV = []
    dTdt = []
    dVdt = []
    SOC_list = []
    V_list = []
    time_list = []
    
    for i in range(segments):
        dTdt.append((df.loc[((i+1)*dt_interval), T_av] - df.loc[(i*dt_interval), T_av]))
        dVdt.append((df.loc[((i+1)*dt_interval), V] - df.loc[(i*dt_interval), V]))
        SOC_list.append(df.loc[((i)*dt_interval):((i+1)*dt_interval), SOC].mean())
        V_list.append(df.loc[((i)*dt_interval):((i+1)*dt_interval), V].mean())
        time_list.append(df.loc[((i)*dt_interval):((i+1)*dt_interval), time].mean())
        dTdV.append(dTdt[i]/dVdt[i])
    
    DTV_data = [[a,b,c,d,e,f] for a,b,c,d,e,f in zip(time_list, V_list, SOC_list, dTdt, dVdt, dTdV)]
    DTV_data_output = pd.DataFrame(DTV_data, columns=[time, V, SOC, "dT/dt (K/s)", "dV/dt(V/s)", "dT/dV (K/V)"])
    
    return DTV_data_output      
        

#A function for calculating the differential resistance (dR/dQ). This function requires an "Reff-R0 (Ohms)" column in the input DF (i.e. like that output from the 'overvoltage_new' function. The function uses a finite-difference method with default dQ (or dSOC) value of 1%. The output is a DF with 'SOC' and 'dR/dSOC (Ohms/%)' columns.
def dRdQ_calc(df_data, dQ_range = 1, Qtotal = 100):
    
    Segments = int(Qtotal/dQ_range)
    dRdQ = []
    soc = []
    
    #dRdQ = [[r[i]/q[i] for q, r in df_data] for i in  ]
    for i in range(Segments):
        sub_df = df_data[(df_data[SOC] < (Qtotal-dQ_range*i)) & (df_data[SOC] > (Qtotal-dQ_range*(i+1)))]
        #sub_df.set_index((range(len(sub_df))), inplace=True)
        Qcol = sub_df[SOC]
        Rcol = sub_df[Rdyn]
        if len(Qcol) < 1:
            break
        dQ = -(Qcol.iloc[-1] - Qcol.iloc[0])
        dR = Rcol.iloc[-1] - Rcol.iloc[0]
        dRdQ.append(dR/dQ)
        soc.append(Qtotal-dQ_range*i)
    
    data = [[l,m] for l, m in zip(soc, dRdQ)]
    data_output = pd.DataFrame(data, columns=[SOC, "dR/dSOC (Ohms/%)"])
    return data_output        

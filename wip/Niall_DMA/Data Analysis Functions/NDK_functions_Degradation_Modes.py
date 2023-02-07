#!/usr/bin/env python
# coding: utf-8

# In[1]:


#importing Pandas etc
import warnings
warnings.simplefilter('ignore', FutureWarning)

import pandas as pd
import numpy as np
import timeit
from scipy import optimize

# In[2]:

#Define some variables to be used

time = "Time (s)"
V = "Voltage (V)"
I = "Current (mA)"
Qtot = "Charge Step (mA.h)"
T = "Temperature (degC)"
Qdis = "Charge (mA.h)"
r0 = "R0 (Ohms)"
NS = "Ns changes"
SOC = "SOC (%)"
OCV = "OCV (V)"




def calc_full_cell_OCV_standalone(half_data, SOC_points, z_pe_lo, z_ne_lo, z_pe_hi, z_ne_hi):
    
    #Unpack halfdata into anode and cathode
    anode, cathode = half_data

        #Unpack the input data, make copies, and make 'z' the index column ('z' rounded to 4 dp's)
        #z_pe_lo, z_ne_lo, z_pe_hi, z_ne_hi = z_vals
    pe_data=cathode.copy()
    pe_data['z'] = pe_data['z'].round(5)
    pe_data.set_index('z', inplace=True)
    pe_data = pe_data.loc[~pe_data.index.duplicated(), :]
    ne_data=anode.copy()
    ne_data['z'] = ne_data['z'].round(5)
    ne_data.set_index('z', inplace=True)
    ne_data = ne_data.loc[~ne_data.index.duplicated(), :]

        #make linearly spaced values of z for neg and pos electrode, with 10001 points in each
    z_ne = np.linspace(z_ne_lo, z_ne_hi, 100001)
    z_pe = np.linspace(z_pe_lo, z_pe_hi, 100001)

        #make new dataframes for the interpolation
    ne_data_int = pd.DataFrame(data=z_ne, columns=['z']) # Make a new DF which has the linearly spaced 'z' values
        #ne_data_int = ne_data_int[(ne_data_int['z'] >= z_ne_lo) & (ne_data_int['z'] <= z_ne_hi)]
    ne_data_int['z'] = ne_data_int['z'].round(5) # Get rid of rounding errors
    ne_data_int.set_index('z', inplace=True) # Make 'z' the index for matching 'OCV' against
    ne_data_int['OCV'] = ne_data['OCV'] # Get 'OCV' values from the 1/2 cell data by matching 'z' values (the index)
    ne_data_int.interpolate(inplace=True) # Interpolate the 'OCV' for all 'z' values (instead of nly matched ones)
    ne_data_int.reset_index(inplace=True) # Reset index for later calculations
        #print(ne_data_int.head(60))

        #same for pos electrode
    pe_data_int = pd.DataFrame(data=z_pe, columns=['z']) # Exactly the same as above, but for positive electrode
        #pe_data_int = pe_data_int[(pe_data_int['z'] >= z_ne_lo) & (pe_data_int['z'] <= z_ne_hi)]
    pe_data_int['z'] = pe_data_int['z'].round(5)
    pe_data_int.set_index('z', inplace=True)
    pe_data_int['OCV'] = pe_data['OCV']
    pe_data_int.interpolate(inplace=True)
    pe_data_int.reset_index(inplace=True)

        #calculate full cell ocv from 1/2 cells
    cell_ocv = pd.DataFrame(data=None) # New DF for calculating the full cell OCV
    cell_ocv['OCV'] = pe_data_int['OCV'] - ne_data_int['OCV'] # Calculate the potential difference between PE and NE based on indices (NOT 'z')
    cell_ocv[SOC] = np.linspace(0,1, len(cell_ocv['OCV'])) # Create an SOC column, ranging from 0 to 1
    cell_ocv.set_index(SOC, inplace=True) # Make the SOC column the index for next section
        #print(cell_ocv.tail())

        #interpolate newly calculated full cell OCV based on SOC values input into function
        #soc_data = np.unique(SOC_points)
        #SOC_points.drop_duplicates(inplace=True)
    print(len(SOC_points))
    cell_out = pd.DataFrame(data=None, index=SOC_points.round(5)) # Make a new DF, using the SOC of the input full cell data as the index
    cell_out['OCV'] = cell_ocv['OCV'] # Get 'OCV' values from the calculated full cell data based on the SOC (index) of the input data
    #cell_out.drop_duplicates(inplace=True)
    cell_out.interpolate(inplace=True) # Interpolate over all SOC values from input data (index)
    cell_out.reset_index(inplace=True) # Reset the index
        
    #output = np.array(cell_out['OCV'])

    return cell_out, ne_data_int, pe_data_int
    
    
def stoich_OCV_fit(anode_data, cathode_data, full_cell, z_guess=[0.1, 0.02, 0.9, 0.85]):
    
    # Convert V and SOC data from full_cell dataset into numpy arrays (seems to be required for the scipy.optimize.curve_fit function)
    cell_V = np.array(full_cell[V])
    cell_SOC = np.array(full_cell[SOC])
    
    # Make copies of the anode and cathode input datasets, round the 'z' values to 5 dp's, then make that the index (required for index-wise operations in 'calc_full_cell_OCV' function)
    pe_data=cathode_data.copy()
    pe_data['z'] = pe_data['z'].round(5)
    pe_data.set_index('z', inplace=True)
    pe_data = pe_data.loc[~pe_data.index.duplicated(), :]
    ne_data=anode_data.copy()
    ne_data['z'] = ne_data['z'].round(5)
    ne_data.set_index('z', inplace=True)
    ne_data = ne_data.loc[~ne_data.index.duplicated(), :]
    
    # Define the function which calculates a full cell OCV based on the 1/2 cell data (exact copy of 'calc_full_cell_OCV_standalone')
    # This function is subsequently passed into the scipy.optimize.curve_fit function below
    # As well as the arguments listed in the brackets, it also requires 1/2 cell data formatted above
    # 'SOC_points' taken from full_cell dataset, 'z_pe_lo' etc taken from z_guess
    def calc_full_cell_OCV(SOC_points, z_pe_lo, z_ne_lo, z_pe_hi, z_ne_hi):

        #make linearly spaced values of z for neg and pos electrode, with 10001 points in each
        z_ne = np.linspace(z_ne_lo, z_ne_hi, 100001) # It should be noted that 'z' for the PE does not correspond directy to lithiation fraction here (it's 1-lithiation fraction)
        z_pe = np.linspace(z_pe_lo, z_pe_hi, 100001)

        #make new dataframes for the interpolation
        ne_data_int = pd.DataFrame(data=z_ne, columns=['z']) # Make a new DF which has the linearly spaced 'z' values
        ne_data_int['z'] = ne_data_int['z'].round(5) # Get rid of rounding errors
        ne_data_int.set_index('z', inplace=True) # Make 'z' the index for matching 'OCV' against
        ne_data_int['OCV'] = ne_data['OCV'] # Get 'OCV' values from the 1/2 cell data by matching 'z' values (the indices)
        ne_data_int.interpolate(inplace=True) # Interpolate the 'OCV' for all 'z' values (instead of only matched ones)
        ne_data_int.reset_index(inplace=True) # Reset index for later calculations

        #same for pos electrode
        pe_data_int = pd.DataFrame(data=z_pe, columns=['z']) # Exactly the same as above, but for positive electrode
        pe_data_int['z'] = pe_data_int['z'].round(5)
        pe_data_int.set_index('z', inplace=True) # Note: Collapse these 3 lines into 1?
        pe_data_int['OCV'] = pe_data['OCV']
        pe_data_int.interpolate(inplace=True)
        pe_data_int.reset_index(inplace=True)

        #calculate full cell ocv from 1/2 cell datasets
        cell_ocv = pd.DataFrame(data=None) # New DF for calculating the full cell OCV
        cell_ocv['OCV'] = pe_data_int['OCV'] - ne_data_int['OCV'] # Calculate the potential difference between PE and NE based on indices (NOT 'z')
        cell_ocv[SOC] = np.linspace(0,1, len(cell_ocv['OCV'])) # Create an SOC column, ranging from 0 to 1
        cell_ocv.set_index(SOC, inplace=True) # Make the SOC column the index for next section

        # The calculated full cell data is then made to match the input (measured) full cell data in terms of SOC (so that they have the same number of datapoints for the optimisation)
        cell_out = pd.DataFrame(data=None, index=SOC_points.round(5)) # Make a new DF, using the SOC of the input (measured) full cell data as the index
        cell_out['OCV'] = cell_ocv['OCV'] # Get 'OCV' values from the calculated full cell data based on the SOC (index) of the input data
        cell_out.interpolate(inplace=True) # Interpolate over all SOC values from input data (index)
        cell_out.reset_index(inplace=True) # Reset the index (not necessary - delete?)

        # For the scipy.optimize.curve_fit function, the output needs to be in the form of a numpy array for some reason
        output = np.array(cell_out['OCV'])

        return output
    
    # This is the actual optimisation
    # It takes the above function, x- and y-values from the full_cell dataset, and the initial_guess of the z_values provided in z_guess
    # It returns the fitted z_values (lithiation fractions), and a covariance matrix
    z_out, z_cov = optimize.curve_fit(calc_full_cell_OCV, xdata=cell_SOC, ydata=cell_V, p0=z_guess, bounds=([0.,0.,0.,0.], [1.,1.,1.,1.]), diff_step=0.001)
    
    # The fitted z_values are then used for calculating the capacities of the full and 1/2 cells, and the offset
    PE_lo, NE_lo, PE_hi, NE_hi = z_out
    cell_capacity = full_cell[Qdis].max() # Cell capacity read directly from full_cell dataset
    ano_cap = cell_capacity/(NE_hi - NE_lo) # NE capacity is full cell capacity divided by the lithiation range of the NE (upper limit minus lower limit)
    cat_cap = cell_capacity/(PE_hi - PE_lo) # PE capacity is full cell capacity divided by the lithiation range of the PE (upper limit minus lower limit)
    offset = (z_out[0] * cat_cap) - (z_out[1] * ano_cap) # Offset is (lower bound of PE lithitation * PE cap) minus (lower bound of NE lithitation * NE cap)
    stoic_params = [cell_capacity, cat_cap, ano_cap, offset] # Save the above paramters into a list
    
    return z_out, z_cov, stoic_params # Returns the z_values, covariance of the fit, and the stoich parameters


def DM_calc(neg_el, pos_el, BoL_cell, aged_cells):
# Function takes in 1/2 cell datasets along with a BoL and some aged cell datasets
    
    #Format the BoL full cell data
    BoL_cell_data = BoL_cell[BoL_cell[I]<0].loc[:, ['Charge (mA.h)', V]]
    BoL_cell_data.reset_index(inplace=True, drop=True)
    BoL_cell_data[SOC] = 1 - BoL_cell_data['Charge (mA.h)']/BoL_cell_data['Charge (mA.h)'].max()
    
    # Perform the optimisation fit for the BoL data
    z_BoL, z_BoL_cov, BoL_params = stoich_OCV_fit(anode_data=neg_el, cathode_data=pos_el, full_cell=BoL_cell_data)
    
    # Make a list of z_values and cell_parameters and populate with BoL data
    z_list = [z_BoL]
    param_list = [BoL_params]
    
    # Perform the optimisation fit for each of the aged cell datasets (using the z_BoL calculated above as the initial z_guess)
    for aged_data in aged_cells:
        # Format the aged full cell data
        aged_cell_data = aged_data[aged_data[I]<0].loc[:, ['Charge (mA.h)', V]]
        aged_cell_data.reset_index(inplace=True, drop=True)
        aged_cell_data[SOC] = 1 - aged_cell_data['Charge (mA.h)']/aged_cell_data['Charge (mA.h)'].max()
        
        #Perform the optimisation fit for the aged full cell data
        z_EoL, z_EoL_cov, aged_params = stoich_OCV_fit(anode_data=neg_el, cathode_data=pos_el, full_cell=aged_cell_data, z_guess=z_BoL)
        z_list.append(z_EoL) # Add the aged data to the list of z_values
        param_list.append(aged_params) # Add the aged data to the list of cell_parameters
    
    # Make dataframes using the lists complied above
    z_parameter_df = pd.DataFrame(data=z_list, columns=['PE_lo', 'NE_lo', 'PE_hi', 'NE_hi'])
    stoic_parameter_df = pd.DataFrame(data=param_list, columns=['Cell Capacity', 'PE Capacity', 'NE Capacity', 'Offset'])
    
    # Calculate DM parameters from the stoic_parameter dataframe above
    SoH = stoic_parameter_df['Cell Capacity']/stoic_parameter_df['Cell Capacity'][0]
    LAM_pe = 1 - (stoic_parameter_df['PE Capacity']/stoic_parameter_df['PE Capacity'][0])
    LAM_ne = 1 - (stoic_parameter_df['NE Capacity']/stoic_parameter_df['NE Capacity'][0])
    LLI = (stoic_parameter_df['PE Capacity'][0] - stoic_parameter_df['PE Capacity'] - (stoic_parameter_df['Offset'][0]-stoic_parameter_df['Offset']))/stoic_parameter_df['Cell Capacity'][0]
    
    # Compile the DM parameters into a dataframe
    DM_df = pd.DataFrame(data={'SoH':SoH, 'LAM PE':LAM_pe, 'LAM NE':LAM_ne, 'LLI':LLI})
    
    return DM_df, stoic_parameter_df # Output 2 dataframes: one with degradation modes, and another with capacities and offsets  


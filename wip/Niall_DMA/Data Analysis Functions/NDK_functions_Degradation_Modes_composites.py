#!/usr/bin/env python
# coding: utf-8

# Mark Ruieh: In this script, multi or multi_comp means two phase in graphite, compared with single phase
#                             long means do optimization / others for different RPTs/BoL
#                             standalone means only one cell one RPT

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

#The first 3 functions are for electrode-level calculations; all subsequent functions are for full cell calculations.

#Function for calculating the composition of a composite electrode based on the 1/2 cell voltage curves of the composite electrode and the two components. Should be updated to make use of some of the functions further below.
def composite_fit(component1_data, component2_data, composite_el_data, gr_cap_guess=0.85):
    
    #make a copy of the data for the composite electrode
    composite_el_data_edit = composite_el_data.copy()
    
    # Make copies of the two component input datasets (component1 and component2), round the 'z' values to 5 dp's, then make that the index (required for index-wise operations in 'calc_full_cell_OCV' function)
    c1_data=component1_data.copy()
    c1_data['OCV'] = c1_data['OCV'].round(4)
    c1_data.drop_duplicates(subset=['OCV'], inplace=True) # Duplicates cause issues when part of an index
    c1_data.set_index('OCV', inplace=True)
    c2_data=component2_data.copy()
    c2_data['OCV'] = c2_data['OCV'].round(4)
    c2_data.drop_duplicates(subset=['OCV'], inplace=True)
    c2_data.set_index('OCV', inplace=True)
    
    # Find the upper and lower voltage values used in the data fit; this is determined by the dataset with the lowest/highest voltage values, respectively
    if composite_el_data_edit['OCV'].max() < component1_data['OCV'].max() and composite_el_data_edit['OCV'].max() < component2_data['OCV'].max():
        V_upper = composite_el_data_edit['OCV'].max()
    elif component1_data['OCV'].max() < composite_el_data_edit['OCV'].max() and component1_data['OCV'].max() < component2_data['OCV'].max():
        V_upper = component1_data['OCV'].max()
    else:
        V_upper = component2_data['OCV'].max()
    
    if composite_el_data_edit['OCV'].min() > component1_data['OCV'].min() and composite_el_data_edit['OCV'].min() > component2_data['OCV'].min():
        V_lower = composite_el_data_edit['OCV'].min()
    elif component1_data['OCV'].min() > composite_el_data_edit['OCV'].min() and component1_data['OCV'].min() > component2_data['OCV'].min():
        V_lower = component1_data['OCV'].min()
    else:
        V_lower = component2_data['OCV'].min()
    
    # Set the voltage range used in the fitting, based on the upper and lower voltage limits specified above
    V_range = np.linspace(V_lower, V_upper, 10001)
    
    # Convert V and SOC data from full_cell dataset into numpy arrays (seems to be required for the scipy.optimize.curve_fit function)
    composite_el_data_edit = composite_el_data_edit[((composite_el_data_edit['OCV']<V_upper) & (composite_el_data_edit['OCV']>V_lower))]
    el_V = np.array(composite_el_data_edit['OCV'])
    el_cap = np.array(composite_el_data_edit['z'])
    
    # Define the function which calculates a composite V vs Q curve based on the two component curves (exact copy of 'calc_electrode_curve' function)
    # This function is subsequently passed into the scipy.optimize.curve_fit function below
    # As well as the arguments listed in the brackets, it also requires the V_range specified above
    # 'z_points' taken from 'composite_el_data' dataset (el_cap), 'cap_comp1' and 'cap_comp2' taken from 'relative_cap_guess'
    def calc_electrode_cap(z_points, cap_comp1):
        
        cap_comp2 = 1. - cap_comp1

        #make new dataframes for the interpolation (component1 here)
        c1_data_int = pd.DataFrame(data=V_range, columns=['OCV']) # Make a new DF which has the linearly spaced 'OCV' values between the voltage limits specified above
        c1_data_int['OCV'] = c1_data_int['OCV'].round(4) # Get rid of rounding errors
        c1_data_int.drop_duplicates(subset=['OCV'], inplace=True) # Get rid of duplicates
        c1_data_int.set_index('OCV', inplace=True) # Make 'OCV' the index for matching 'z' against
        c1_data_int['z'] = c1_data['z']*cap_comp1 # Multiply the 'z' data of component1 by the guessed (fractional) capacity of component1 at values where the OCVs match (i.e. index-based)
        c1_data_int.interpolate(inplace=True) # Interpolate the 'z' for all 'OCV' values (instead of only matched ones)

        #same for component2
        c2_data_int = pd.DataFrame(data=V_range, columns=['OCV']) # Exactly the same as above, but for component2
        c2_data_int['OCV'] = c2_data_int['OCV'].round(4)
        c2_data_int.drop_duplicates(subset=['OCV'], inplace=True)
        c2_data_int.set_index('OCV', inplace=True)
        c2_data_int['z'] = c2_data['z']*cap_comp2
        c2_data_int.interpolate(inplace=True)

        #calculate the composite electrode 'z' values from the two components
        cell_cap = pd.DataFrame(data=V_range, columns=['OCV']) # New DF for calculating the composite electrode 'z'; uses the same V_range values as the composites
        cell_cap['OCV'] = cell_cap['OCV'].round(4) # Get rid of rounding errors
        cell_cap.drop_duplicates(subset=['OCV'], inplace=True) # Remove any duplicates
        cell_cap.set_index('OCV', inplace=True) # Set the 'OCV' to be the index
        cell_cap['z'] = c1_data_int['z'] + c2_data_int['z'] # Calculate the sum of the capacities (or 'z' values) based on the OCV values (index)
        cell_cap['z'] = cell_cap['z'].round(4) # Get rid of rounding errors
        cell_cap.reset_index(inplace=True) # Makes 'OCV' a normal data column again instead of the index
        cell_cap.drop_duplicates(subset=['z'], inplace=True) # Get rid of any duplicates in the 'z' values
        cell_cap.set_index('z', inplace=True) # Make the 'z' column the index for next section

        # The calculated composite electrode data is then made to match the input (measured) composite electrode data in terms of 'z' (so that they have the same number of datapoints for the optimisation)
        cell_out = pd.DataFrame(data=None, index=z_points.round(4)) # Make a new DF, using the 'z' of the input (measured) composite electrode data as the index
        cell_out['OCV'] = cell_cap['OCV'] # Get 'OCV' values from the calculated composite electrode data based on the 'z' (index) of the input data
        cell_out.interpolate(inplace=True) # Interpolate over all 'z' values from input data (index)
        cell_out.reset_index(inplace=True) # Reset the index (not necessary - delete?)

        # For the scipy.optimize.curve_fit function, the output needs to be in the form of a numpy array for some reason
        output = np.array(cell_out['OCV'])
        
        print(cap_comp1, cap_comp2)
        
        # An additional statement which gives a huge cost (i.e. error in the fit) if the relative capacities of the two components don't add up to 1 (they are meant ot be %s of total)
        #if (cap_comp1 + cap_comp2) != 1.:
        #    output = np.ones((len(el_cap)))*10e5

        return output
    
    # This is the actual optimisation
    # It takes the above function, x- and y-values from the composite_el_data dataset, and the initial_guess of the capacities of the components provided in relative_cap_guess
    # It returns the fitted relative_caps (capacity fractions), and a covariance matrix
    z_out, z_cov = optimize.curve_fit(calc_electrode_cap, xdata=el_cap, ydata=el_V, p0=gr_cap_guess, bounds=(0., 1.), diff_step=0.1)
    
    return z_out, z_cov#, Cap_out # Returns the z_values, covariance of the fit, and the stoich parameters


#Function for calculation the residual between the real and calculated 1/2 cell voltage curves of a composite electrode (in terms of 'z').
def el_fit_error_check(NE_comp1_data, NE_comp2_data, electrode_data, fit_results):
    
    el_data = electrode_data.copy()
    
    el_calc_data = calc_electrode_curve2(NE_comp1_data, NE_comp2_data, el_data, *fit_results)
    
    el_data['OCV'] = el_data['OCV'].round(4)
    el_calc_data['OCV'] = el_calc_data['OCV'].round(4)
    print(el_calc_data.head())
    
    el_calc_data.drop_duplicates(subset=['OCV'], inplace=True)
    el_data.drop_duplicates(subset=['OCV'], inplace=True)
    
    if len(el_calc_data['OCV']) < len(el_data['OCV']):
        short_data = el_calc_data
        diff = pd.DataFrame(data=short_data['OCV'])
        diff.set_index('OCV', inplace=True)
        el_data.set_index('OCV', inplace=True)
        el_calc_data.set_index('OCV', inplace=True)
        
        diff['z error'] = el_data['z'] - el_calc_data['z']
        
        err_array = np.array(diff['z error'])
        rmse_result = np.sqrt(np.square(err_array).mean())
    
    else:
        rmse_result = 10e5
    
    return rmse_result, diff


#Function for calculation the residual between the real and calculated 1/2 cell voltage curves of a composite electrode (in terms of 'z').
def el_fit_error_check_V(NE_comp1_data, NE_comp2_data, electrode_data, fit_results):
    
    el_data = electrode_data.copy()
    
    el_calc_data = calc_electrode_curve2(NE_comp1_data, NE_comp2_data, el_data, *fit_results)
    
    el_data['z'] = el_data['z'].round(3)
    el_calc_data['z'] = el_calc_data['z'].round(3)
    print(el_calc_data.head())
    
    el_calc_data.drop_duplicates(subset=['z'], inplace=True)
    el_data.drop_duplicates(subset=['z'], inplace=True)
    
    if len(el_calc_data['z']) < len(el_data['z']):
        short_data = el_calc_data
        diff = pd.DataFrame(data=short_data['z'])
        diff = diff[diff['z']>0.04]
        diff.set_index('z', inplace=True)
        el_data.set_index('z', inplace=True)
        el_calc_data.set_index('z', inplace=True)
        
        diff['V error'] = el_data['OCV'] - el_calc_data['OCV']
        
        err_array = np.array(diff['V error'])
        rmse_result = np.sqrt(np.square(err_array).mean())
    
    else:
        rmse_result = 10e5
    
    return rmse_result, diff

#The below are for full cell calculations

#Function for formatting data for two components of a composite half cell. Replaced by "format_el_component_data1" function. The difference between the two is whether the voltage limits are determined by the lowest or highest high and low voltage. This version takes the highest high and lowest low, with extrapolation for the limiting dataset.
def format_el_component_data(component1_data, component2_data):
    
    # Make copies of the anode and cathode input datasets, round the 'z' values to 5 dp's, then make that the index (required for index-wise operations in 'calc_full_cell_OCV' function)
    c1_data=component1_data.copy()
    c2_data=component2_data.copy()
    
    
    # Add points to the end of the datasets with higher/lower lower/upper voltage limits, so that the 'calc_electrode_curve function can operate over the full voltage range. Also, use those limits to make a linearly-spaced voltage series between those limits called 'V_ranges'
    if component1_data['OCV'].max() < component2_data['OCV'].max():
        V_upper = component2_data['OCV'].max()
        c1_data.loc[0,'OCV'] = V_upper
    else:
        V_upper = component1_data['OCV'].max()
        c2_data.loc[0,'OCV'] = V_upper
    
    if component1_data['OCV'].min() > component2_data['OCV'].min():
        V_lower = component2_data['OCV'].min()
        c1_data.loc[c1_data.index.max(),'OCV'] = V_lower
    else:
        V_lower = component1_data['OCV'].min()
        c2_data.loc[c2_data.index.max(),'OCV'] = V_lower
    
    V_ranges = np.unique(np.linspace(V_lower, V_upper, 10001).round(decimals=4))
    
    #c1_data['z'] = (c1_data['z'] - c1_data['z'].min())/(c1_data['z'].max() - c1_data['z'].min())
    c1_data['OCV'] = c1_data['OCV'].round(4)
    c1_data.drop_duplicates(subset=['OCV'], inplace=True) # Duplicates cause issues when part of an index
    c1_data.set_index('OCV', inplace=True)
    
    #c2_data['z'] = (c2_data['z'] - c2_data['z'].min())/(c2_data['z'].max() - c2_data['z'].min())
    c2_data['OCV'] = c2_data['OCV'].round(4)
    c2_data.drop_duplicates(subset=['OCV'], inplace=True)
    c2_data.set_index('OCV', inplace=True)

    return c1_data, c2_data, V_ranges


#Function for formatting data for two components of a composite half cell. This replaces the "format_el_component_data" function. The difference between the two is whether the voltage limits are determined by the lowest or highest high and low voltage. This version takes the lowest high and highest low, thereby not requiring extrapolation.
def format_el_component_data1(component1_data, component2_data):
    
    # Make copies of the anode and cathode input datasets, round the 'z' values to 5 dp's, then make that the index (required for index-wise operations in 'calc_full_cell_OCV' function)
    c1_data=component1_data.copy()
    c2_data=component2_data.copy()
    
    
    # Add points to the end of the datasets with higher/lower lower/upper voltage limits, so that the 'calc_electrode_curve function can operate over the full voltage range. Also, use those limits to make a linearly-spaced voltage series between those limits called 'V_ranges'
    if component1_data['OCV'].max() < component2_data['OCV'].max():
        V_upper = component1_data['OCV'].max()
    else:
        V_upper = component2_data['OCV'].max()
    
    if component1_data['OCV'].min() > component2_data['OCV'].min():
        V_lower = component1_data['OCV'].min()
    else:
        V_lower = component2_data['OCV'].min()
    
    c1_data = c1_data[((c1_data['OCV']<=V_upper) & (c1_data['OCV']>=V_lower))]
    c2_data = c2_data[((c2_data['OCV']<=V_upper) & (c2_data['OCV']>=V_lower))]
    
    V_ranges = np.unique(np.linspace(V_lower, V_upper, 10001).round(decimals=4))
    
    c1_data['z'] = (c1_data['z'] - c1_data['z'].min())/(c1_data['z'].max() - c1_data['z'].min())
    c1_data['OCV'] = c1_data['OCV'].round(4)
    c1_data.drop_duplicates(subset=['OCV'], inplace=True) # Duplicates cause issues when part of an index
    c1_data.set_index('OCV', inplace=True)
    
    c2_data['z'] = (c2_data['z'] - c2_data['z'].min())/(c2_data['z'].max() - c2_data['z'].min())
    c2_data['OCV'] = c2_data['OCV'].round(4)
    c2_data.drop_duplicates(subset=['OCV'], inplace=True)
    c2_data.set_index('OCV', inplace=True)

    return c1_data, c2_data, V_ranges



def calc_electrode_curve(component1_data, component2_data, cap_comp1, V_range):
    
    cap_comp2 = 1. - cap_comp1
    
    #make new dataframes for the interpolation (component1 here)
    c1_data_int = pd.DataFrame(data=V_range, columns=['OCV']) # Make a new DF which has the linearly spaced 'OCV' values between the voltage limits specified above
    c1_data_int.set_index('OCV', inplace=True) # Make 'OCV' the index for matching 'z' against
    c1_data_int['z'] = component1_data['z']*cap_comp1 # Multiply the 'z' data of component1 by the guessed (fractional) capacity of component1 at values where the OCVs match (i.e. index-based)
    c1_data_int.interpolate(inplace=True) # Interpolate the 'z' for all 'OCV' values (instead of only matched ones)
    
    #same for component2
    c2_data_int = pd.DataFrame(data=V_range, columns=['OCV']) # Exactly the same as above, but for component2
    c2_data_int.set_index('OCV', inplace=True)
    c2_data_int['z'] = component2_data['z']*cap_comp2
    c2_data_int.interpolate(inplace=True)
    

    #calculate the composite electrode 'z' values from the two components
    electrode_curve = pd.DataFrame(data=V_range, columns=['OCV']) # New DF for calculating the composite electrode 'z'; uses the same V_range values as the composites
    electrode_curve.set_index('OCV', inplace=True) # Set the 'OCV' to be the index
    electrode_curve['z'] = c1_data_int['z'] + c2_data_int['z'] # Calculate the sum of the capacities (or 'z' values) based on the OCV values (index)
    electrode_curve['z'] = electrode_curve['z'].round(5) # Get rid of rounding errors
    electrode_curve.reset_index(inplace=True) # Makes 'OCV' a normal data column again instead of the index
    electrode_curve.drop_duplicates(subset=['z'], inplace=True) # Get rid of any duplicates in the 'z' values
    electrode_curve.set_index('z', inplace=True)

    return electrode_curve




def stoich_OCV_fit_multi_comp(
    anode_comp1_data, anode_comp2_data, cathode_data, full_cell, 
    z_guess=[0.1, 0.002, 0.95, 0.85, 0.84], diff_step_size=0.01):
    
    # Format the negative electrode data (both components)
    ano_c1_data, ano_c2_data, V_range = format_el_component_data1(anode_comp1_data, anode_comp2_data)
    
    # Convert V and SOC data from full_cell dataset into numpy arrays (seems to be required for the scipy.optimize.curve_fit function)
    cell_V = np.array(full_cell[V])
    cell_SOC = np.array(full_cell[SOC])
    
    # Make copies of the positive electrode input dataset, round the 'z' values to 5 dp's, then make that the index (required for index-wise operations in 'calc_full_cell_OCV' function)
    pe_data=cathode_data.copy()
    pe_data['z'] = pe_data['z'].round(5)
    pe_data.set_index('z', inplace=True)
    pe_data = pe_data.loc[~pe_data.index.duplicated(), :]
    
    # Define the function which calculates a full cell OCV based on the 1/2 cell data (exact copy of 'calc_full_cell_OCV_standalone')
    # This function is subsequently passed into the scipy.optimize.curve_fit function below
    # As well as the arguments listed in the brackets, it also requires 1/2 cell data formatted above
    # 'SOC_points' taken from full_cell dataset, 'z_pe_lo' etc taken from z_guess
    def calc_full_cell_OCV_multi(SOC_points, z_pe_lo, z_ne_lo, z_pe_hi, z_ne_hi, comp1_frac):
        
        # Calculate the negative electrode curve based on the fraction of each component (give by 'comp1_frac')
        ne_data = calc_electrode_curve(ano_c1_data, ano_c2_data, comp1_frac, V_range)

        #make linearly spaced values of z for neg and pos electrode, with 10001 points in each
        z_ne = np.unique(np.linspace(z_ne_lo, z_ne_hi, 10001).round(decimals=5)) # It should be noted that 'z' for the PE does not correspond directy to lithiation fraction here (it's 1-lithiation fraction)
        z_pe = np.unique(np.linspace(z_pe_lo, z_pe_hi, 10001).round(decimals=5))

        #make new dataframes for the interpolation
        ne_data_int = pd.DataFrame(data=z_ne, columns=['z']) # Make a new DF which has the linearly spaced 'z' values
        ne_data_int.set_index('z', inplace=True) # Make 'z' the index for matching 'OCV' against
        ne_data_int['OCV'] = ne_data['OCV'] # Get 'OCV' values from the 1/2 cell data by matching 'z' values (the indices)
        ne_data_int.interpolate(inplace=True) # Interpolate the 'OCV' for all 'z' values (instead of only matched ones)
        ne_data_int.reset_index(inplace=True) # Reset index for later calculations

        #same for pos electrode
        pe_data_int = pd.DataFrame(data=z_pe, columns=['z']) # Exactly the same as above, but for positive electrode
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

        # For the scipy.optimize.curve_fit function, the output needs to be in the form of a numpy array for some reason
        output = np.array(cell_out['OCV'])

        return output
     
    # This is the actual optimisation
    # It takes the above function, x- and y-values from the full_cell dataset, 
    # and the initial_guess of the z_values provided in z_guess
    # It returns the fitted z_values (lithiation fractions), and a covariance matrix
    z_out, z_cov = optimize.curve_fit(
        calc_full_cell_OCV_multi, xdata=cell_SOC, ydata=cell_V, 
        p0=z_guess, bounds=([0.,0.,0.,0.,0.5], [1.,1.,1.,1.,1.]), # Mark Ruihe add: constrain Gr_frac to 0.5~1 as we know 
        diff_step=diff_step_size)
    
    # The fitted z_values are then used for calculating the capacities of the full and 1/2 cells, and the offset
    PE_lo, NE_lo, PE_hi, NE_hi, comp1_frac = z_out
    cell_capacity = full_cell[Qdis].max() # Cell capacity read directly from full_cell dataset
    ano_tot_cap = cell_capacity/(NE_hi - NE_lo)
    ano_comp1_cap = ano_tot_cap*comp1_frac # NE capacity is full cell capacity divided by the lithiation range of the NE (upper limit minus lower limit)
    ano_comp2_cap = ano_tot_cap*(1. - comp1_frac) # NE capacity is full cell capacity divided by the lithiation range of the NE (upper limit minus lower limit)
    cat_cap = cell_capacity/(PE_hi - PE_lo) # PE capacity is full cell capacity divided by the lithiation range of the PE (upper limit minus lower limit)
    offset = (z_out[0] * cat_cap) - (z_out[1] * ano_tot_cap) # Offset is (lower bound of PE lithitation * PE cap) minus (lower bound of NE lithitation * NE cap)
    Cap_out = [
        cell_capacity, cat_cap, ano_tot_cap, 
        ano_comp1_cap, ano_comp2_cap, offset] # Save the above paramters into a list
    
    return z_out, z_cov, Cap_out # Returns the z_values, covariance of the fit, and the stoich parameters

# Mark Ruihe: Calculate fitted cell voltage based on the optimised value (one cell, one RPT)
def calc_full_cell_OCV_multi_standalone(
    anode_comp1_data, anode_comp2_data, cathode_data, 
    SOC_points, z_pe_lo, z_ne_lo, z_pe_hi, z_ne_hi, comp1_frac):
    
    ano_c1_data, ano_c2_data, V_range = format_el_component_data1(anode_comp1_data, anode_comp2_data)
    
    # Make copies of the positive electrode input dataset, round the 'z' values to 5 dp's, then make that the index (required for index-wise operations in 'calc_full_cell_OCV' function)
    pe_data=cathode_data.copy()
    pe_data['z'] = pe_data['z'].round(5)
    pe_data.set_index('z', inplace=True)
    pe_data = pe_data.loc[~pe_data.index.duplicated(), :]
    
    # Calculate the negative electrode curve based on the fraction of each component (give by 'comp1_frac')
    ne_data = calc_electrode_curve(ano_c1_data, ano_c2_data, comp1_frac, V_range)

    #make linearly spaced values of z for neg and pos electrode, with 10001 points in each
    z_ne = np.unique(np.linspace(z_ne_lo, z_ne_hi, 10001).round(decimals=5)) # It should be noted that 'z' for the PE does not correspond directy to lithiation fraction here (it's 1-lithiation fraction)
    z_pe = np.unique(np.linspace(z_pe_lo, z_pe_hi, 10001).round(decimals=5))

    #make new dataframes for the interpolation
    ne_data_int = pd.DataFrame(data=z_ne, columns=['z']) # Make a new DF which has the linearly spaced 'z' values
    ne_data_int.set_index('z', inplace=True) # Make 'z' the index for matching 'OCV' against
    ne_data_int['OCV'] = ne_data['OCV'] # Get 'OCV' values from the 1/2 cell data by matching 'z' values (the indices)
    ne_data_int.interpolate(inplace=True) # Interpolate the 'OCV' for all 'z' values (instead of only matched ones)
    ne_data_int.reset_index(inplace=True) # Reset index for later calculations

    #same for pos electrode
    pe_data_int = pd.DataFrame(data=z_pe, columns=['z']) # Exactly the same as above, but for positive electrode
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
    cell_out.reset_index(inplace=True)
    # For the scipy.optimize.curve_fit function, the output needs to be in the form of a numpy array for some reason
    #output = np.array(cell_out['OCV'])

    return cell_out, ne_data_int, ne_data, pe_data_int

# Mark Ruihe: Seems that this function actually doesn't get used!
def DM_calc_multi_comp(neg_el_comp1, neg_el_comp2, pos_el, BoL_cell, aged_cells):
# Function takes in 1/2 cell datasets along with a BoL and some aged cell datasets
    
    #Format the BoL full cell data
    BoL_cell_data = BoL_cell[BoL_cell[I]<0].loc[:, ['Charge (mA.h)', V]]
    BoL_cell_data.reset_index(inplace=True, drop=True)
    BoL_cell_data[SOC] = 1 - BoL_cell_data['Charge (mA.h)']/BoL_cell_data['Charge (mA.h)'].max()
    
    # Perform the optimisation fit for the BoL data
    z_BoL, z_BoL_cov, BoL_params = stoich_OCV_fit_multi_comp(
        anode_comp1_data=neg_el_comp1, anode_comp2_data=neg_el_comp2, 
        cathode_data=pos_el, full_cell=BoL_cell_data)
    
    # Make a list of z_values and cell_parameters and populate with BoL data
    z_list = [z_BoL]
    param_list = [BoL_params]
    cov_matrix = [np.sqrt(np.diag(z_BoL_cov))]
    
    counter_val = 0
    # Perform the optimisation fit for each of the aged cell datasets (using the z_BoL calculated above as the initial z_guess)
    for aged_data in aged_cells:
        # Format the aged full cell data
        aged_cell_data = aged_data[aged_data[I]<0].loc[:, ['Charge (mA.h)', V]]
        aged_cell_data.reset_index(inplace=True, drop=True)
        aged_cell_data[SOC] = 1 - aged_cell_data['Charge (mA.h)']/aged_cell_data['Charge (mA.h)'].max()
        
        #Perform the optimisation fit for the aged full cell data
        z_EoL, z_EoL_cov, aged_params = stoich_OCV_fit_multi_comp(
            anode_comp1_data=neg_el_comp1, anode_comp2_data=neg_el_comp2, 
            cathode_data=pos_el, full_cell=aged_cell_data, 
            z_guess=z_list[counter_val]
            )
        z_list.append(z_EoL) # Add the aged data to the list of z_values
        param_list.append(aged_params) # Add the aged data to the list of cell_parameters
        cov_matrix.append(np.sqrt(np.diag(z_EoL_cov)))
        counter_val += 1
    
    # Make dataframes using the lists complied above
    z_parameter_df = pd.DataFrame(
        data=z_list, columns=['PE_lo', 'NE_lo', 'PE_hi', 'NE_hi', 'Gr_frac'])
    stoic_parameter_df = pd.DataFrame(
        data=param_list, 
        columns=[
        'Cell Capacity', 'PE Capacity', 'NE(tot) Capacity', 
        'NE(Gr) Capacity', 'NE(Si) Capacity', 'Offset'])
    
    # Calculate DM parameters from the stoic_parameter dataframe above
    SoH = stoic_parameter_df['Cell Capacity']/stoic_parameter_df['Cell Capacity'][0]
    LAM_pe = 1 - (stoic_parameter_df['PE Capacity']/stoic_parameter_df['PE Capacity'][0])
    LAM_ne_tot = 1 - (stoic_parameter_df['NE(tot) Capacity']/stoic_parameter_df['NE(tot) Capacity'][0])
    LAM_ne_Gr = 1 - (stoic_parameter_df['NE(Gr) Capacity']/stoic_parameter_df['NE(Gr) Capacity'][0])
    LAM_ne_Si = 1 - (stoic_parameter_df['NE(Si) Capacity']/stoic_parameter_df['NE(Si) Capacity'][0])
    LLI = (
        stoic_parameter_df['PE Capacity'][0] 
        - stoic_parameter_df['PE Capacity'] 
        - (stoic_parameter_df['Offset'][0]-stoic_parameter_df['Offset'])
        )/stoic_parameter_df['Cell Capacity'][0]
    
    # Compile the DM parameters into a dataframe
    DM_df = pd.DataFrame(
        data={
        'SoH':SoH, 'LAM PE':LAM_pe, 'LAM NE_tot':LAM_ne_tot, 
        'LAM NE_Gr':LAM_ne_Gr, 'LAM NE_Si':LAM_ne_Si, 'LLI':LLI})
    
    return DM_df, stoic_parameter_df#, cov_matrix # Output 2 dataframes: one with degradation modes, and another with capacities and offsets


# Mark Ruihe: Calculate error used the fitted result; for one cell one RPT
def DM_error_check(NE_comp1_data, NE_comp2_data, PE_data, cell_data, fit_results):
    
    cell_calc_data, _, _, _ = calc_full_cell_OCV_multi_standalone(
        NE_comp1_data, NE_comp2_data, PE_data, 
        cell_data[SOC], *fit_results)
    
    if len(cell_calc_data['OCV']) == len(cell_data[V]):
        diff = pd.DataFrame(data=cell_data[SOC])
        diff['V error'] = cell_data[V] - cell_calc_data['OCV']
        
        err_array = np.array(diff['V error'])
        rmse_result = np.sqrt(np.square(err_array).mean())
    
    else:
        rmse_result = 10e5
    
    return rmse_result



# Mark Ruihe: Main function for one cell, all RPTs and BoL --> Big changes made!
def DM_calc_multi_comp_long(neg_el_comp1, neg_el_comp2, pos_el, BoL_cell, aged_cells, carry_guess=True):
# Function takes in 1/2 cell datasets along with a BoL and some aged cell datasets
    
    #Format the BoL full cell data
    BoL_cell_data = BoL_cell[BoL_cell[I]<0].loc[:, ['Charge (mA.h)', V]]
    BoL_cell_data.reset_index(inplace=True, drop=True)
    BoL_cell_data[SOC] = (
        1 - 
        BoL_cell_data['Charge (mA.h)']/BoL_cell_data['Charge (mA.h)'].max())
    
    # Perform the optimisation fit for the BoL data
    # Mark Ruihe: Change! Repeat this step for more times (1000?)
    z_BoL, z_BoL_cov, BoL_params = stoich_OCV_fit_multi_comp(
        anode_comp1_data=neg_el_comp1, 
        anode_comp2_data=neg_el_comp2, 
        cathode_data=pos_el, full_cell=BoL_cell_data)
    
    # Calculate error of the fit
    err_BoL = DM_error_check(
        neg_el_comp1, neg_el_comp2, 
        pos_el, BoL_cell_data, z_BoL)
    
    # Make a list of z_values and cell_parameters and populate with BoL data
    # Mark Ruihe: this is for later optimization - choose the best from the list
    z_list = [z_BoL]
    param_list = [BoL_params]
    #cov_matrix = [np.sqrt(np.diag(z_BoL_cov))]
    err_list = [err_BoL]
    
    counter_val = 0
    
    # Perform the optimisation fit for each of the aged cell datasets (
    #           using the z_BoL calculated above as the initial z_guess)
    for aged_data in aged_cells:
        print(counter_val)
        # Format the aged full cell data
        aged_cell_data = aged_data[aged_data[I]<0].loc[:, ['Charge (mA.h)', V]] # choose only the 
        aged_cell_data.reset_index(inplace=True, drop=True)
        aged_cell_data[SOC] = (
            1 - 
            aged_cell_data['Charge (mA.h)']
            /aged_cell_data['Charge (mA.h)'].max()   )
        
        if not carry_guess:
            guess_values = [0.1, 0.002, 0.95, 0.85, 0.84]
        else:
            guess_values = z_list[counter_val]
        
        #Perform the optimisation fit for the aged full cell data
        z_EoL, z_EoL_cov, aged_params = stoich_OCV_fit_multi_comp(
            anode_comp1_data=neg_el_comp1, anode_comp2_data=neg_el_comp2, 
            cathode_data=pos_el, full_cell=aged_cell_data, z_guess=guess_values)
        
        # Calculate error of fit (RMSE of fitted curve minus actual data)
        err_EoL = DM_error_check(neg_el_comp1, neg_el_comp2, pos_el, aged_cell_data, z_EoL)
        
        iter_val=0  
        while (
            aged_params[4] > param_list[counter_val][4] 
            or 
            aged_params[2] > param_list[counter_val][2]) and iter_val < 10:  # repeat for 10 times！
            iter_val += 1
            z_EoL, z_EoL_cov, aged_params = stoich_OCV_fit_multi_comp(
                anode_comp1_data=neg_el_comp1, 
                anode_comp2_data=neg_el_comp2, cathode_data=pos_el, 
                full_cell=aged_cell_data, 
                z_guess=z_list[counter_val], 
                diff_step_size=0.1)
            # Calculate error of fit (RMSE of fitted curve minus actual data)
            err_EoL = DM_error_check(
                neg_el_comp1, neg_el_comp2, 
                pos_el, aged_cell_data, z_EoL)
            
        
        z_list.append(z_EoL) # Add the aged data to the list of z_values
        param_list.append(aged_params) # Add the aged data to the list of cell_parameters
        #cov_matrix.append(np.sqrt(np.diag(z_EoL_cov)))
        err_list.append(err_EoL)
        counter_val += 1
    
    # Make dataframes using the lists complied above
    z_parameter_df = pd.DataFrame(data=z_list, columns=['PE_lo', 'NE_lo', 'PE_hi', 'NE_hi', 'Gr_frac'])
    stoic_parameter_df = pd.DataFrame(
        data=param_list, 
        columns=[
        'Cell Capacity', 'PE Capacity', 'NE(tot) Capacity', 
        'NE(Gr) Capacity', 'NE(Si) Capacity', 'Offset'])
    err_df = pd.DataFrame(data={'RMSE (V)':err_list})   # RMSE for one cell, all RPTs
    
    # Calculate DM parameters from the stoic_parameter dataframe above
    SoH = stoic_parameter_df['Cell Capacity']/stoic_parameter_df['Cell Capacity'][0]
    LAM_pe = 1 - (stoic_parameter_df['PE Capacity']/stoic_parameter_df['PE Capacity'][0])
    LAM_ne_tot = 1 - (stoic_parameter_df['NE(tot) Capacity']/stoic_parameter_df['NE(tot) Capacity'][0])
    LAM_ne_Gr = 1 - (stoic_parameter_df['NE(Gr) Capacity']/stoic_parameter_df['NE(Gr) Capacity'][0])
    LAM_ne_Si = 1 - (stoic_parameter_df['NE(Si) Capacity']/stoic_parameter_df['NE(Si) Capacity'][0])
    LLI = (stoic_parameter_df['PE Capacity'][0] - stoic_parameter_df['PE Capacity'] - (stoic_parameter_df['Offset'][0]-stoic_parameter_df['Offset']))/stoic_parameter_df['Cell Capacity'][0]
    
    # Compile the DM parameters into a dataframe
    DM_df = pd.DataFrame(
        data={
        'SoH':SoH, 'LAM PE':LAM_pe, 'LAM NE_tot':LAM_ne_tot, 
        'LAM NE_Gr':LAM_ne_Gr, 'LAM NE_Si':LAM_ne_Si, 'LLI':LLI})
    #, cov_matrix # Output 2 dataframes: one with degradation modes, and another with capacities and offsets
    return DM_df, stoic_parameter_df, err_df



# In[ ]:





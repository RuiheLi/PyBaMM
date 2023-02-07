# Echem-analysis-functions

Repo of useful functions to analyse battery cycler data

There are various functions I have written for analysing data from battery electrochemistry experiments. They were primarily written for data from the Biologic, but can handle data from other cyclers if it is formated correctly.

There are 5 main modules: 
* Biologic data loading functions for bringing raw data into pandas dataframes and formating etc. 
* Functions for calculating overvoltage and resistance. 
* Functions for calculating 'differentials' (i.e. ICA, DVA, DTV, DRA, etc). 
* Functions for doing OCV-fitting and calculating degradation modes (standard). 
* Functions for doing OCV-fitting and calculating degradation modes (for cells containing composite electrodes).

Each module contains a bunch of different functions. Most are reliant on the pandas, numpy, and scipy libraries.

I have also included an example jupyter notebook for using the degradation mode analysis (DMA) functions, along with some data (both full and half cell) to be used with it. The notebook explains what each function does and how you can use them for your own data. It covers both the 'traditional' OCV-fitting method as well as the modified version for use with cells which have composite electrodes.

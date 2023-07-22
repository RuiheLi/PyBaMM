import csv, random, os
import pybamm as pb;import pandas as pd;import numpy as np;
import os, json,openpyxl,traceback,multiprocessing,scipy.optimize
import matplotlib.pyplot as plt;
import imageio,timeit,random,time, signal
from scipy.io import savemat,loadmat;
from pybamm import constants,exp;import matplotlib as mpl; 
fs=17; 
font = {'family' : 'DejaVu Sans','size'   : fs}
mpl.rc('font', **font)


parameter_names = [
    "Scan No",
    "Exp No.", # 2,3,5
    "Ageing temperature", # 10,25,40
    'Inner SEI lithium interstitial diffusivity [m2.s-1]',#1e-19~5e-18
    'Dead lithium decay constant [s-1]',
    'Lithium plating kinetic rate constant [m.s-1]',
    'Negative electrode LAM constant proportional term [s-1]',
    'Positive electrode LAM constant proportional term [s-1]',
    'Negative electrode cracking rate',
    'Outer SEI partial molar volume [m3.mol-1]',
    "SEI growth activation energy [J.mol-1]", # 1e4~3.8e4
    "Negative cracking growth activation energy [J.mol-1]",#0
    "Negative electrode diffusivity activation energy [J.mol-1]",#1.7e4
    "Positive electrode diffusivity activation energy [J.mol-1]",#1.2e4
    "Contact resistance [Ohm]",#0.010
    'Total heat transfer coefficient [W.m-2.K-1]',#20
    'Initial electrolyte excessive amount ratio', # 1.0 or 0.99
    "Ratio of lithium moles to SEI moles", # 2.0
]


import tkinter as tk
from tkinter import ttk

def on_selection(event):
    selected_parameter = dropdown.get()
    print("Selected Parameter:", selected_parameter)

root = tk.Tk()
root.title("GUI with Dropdown List")

# Create the dropdown list
dropdown_label = ttk.Label(root, text="Select a Parameter:")
dropdown_label.pack(pady=10)

dropdown = ttk.Combobox(root, values=parameter_names)
dropdown.pack()

dropdown.bind("<<ComboboxSelected>>", on_selection)

root.mainloop()

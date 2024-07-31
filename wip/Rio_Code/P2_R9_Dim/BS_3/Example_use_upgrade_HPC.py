# Load modules
import pybamm as pb;import pandas as pd;import numpy as np;
import os, json,openpyxl,traceback,multiprocessing,scipy.optimize,sys,gc
import matplotlib.pyplot as plt;
import pickle,imageio,timeit,random,time, signal
from scipy.io import savemat,loadmat;
from pybamm import constants,exp;import matplotlib as mpl


str_path_0 = os.path.abspath(os.path.join(pb.__path__[0],'..'))
str_path_1 = os.path.abspath(
    os.path.join(str_path_0,"wip/Rio_Code/Fun_Upgrade"))
print(str_path_1)

from wip.Rio_Code.Fun_Upgrade.TaskResult import TaskResult
from wip.Rio_Code.Fun_Upgrade.Plot import *
from wip.Rio_Code.Fun_Upgrade.Run_model import *
from wip.Rio_Code.Fun_Upgrade.Fun_Upgrade import *
from wip.Rio_Code.Fun_Upgrade.Post_process import *

print("import finish?")
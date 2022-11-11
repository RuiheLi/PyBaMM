import pybamm;import pandas as pd   ;import numpy as np;import os;import matplotlib.pyplot as plt;import os;#import imageio
from scipy.io import savemat,loadmat;from pybamm import constants,exp;
import matplotlib as mpl; from textwrap import wrap
fs=17; # or we can set import matplotlib.pyplot as plt then say 'mpl.rc...'
import openpyxl
import traceback
import multiprocessing
font = {'family' : 'DejaVu Sans','size'   : fs}
mpl.rc('font', **font)

import random
def may_cause_error():
	if (random.random() > 0.8):
		return 1 / 0
	else:
		return 0
# 思路：设置若干秒后发出信号量，收到信号量以后抛出异常，在异常捕获语句中保存数据
import time, signal

'''
定义“超时异常”，这个异常被抛出，就表示运行超时/卡住了
'''
class TimeoutError(Exception):
	pass
'''
接收到信号，就抛出超时异常
'''	
def handle_signal(signal_num, frame):
	raise TimeoutError
Index = [1,2,3,4,5]
for i in Index:
    signal.signal(signal.SIGALRM, handle_signal)
    signal.alarm(2)
    try:
        time.sleep(random.uniform(1,3))
        may_cause_error()
    except TimeoutError:
        print(f'Takes too long for Scan {i}, kill and save your data here')
    except ZeroDivisionError:
        print(f'Fail for Scan {i}')
    else:
        print(f'Succeed for Scan {i}!')

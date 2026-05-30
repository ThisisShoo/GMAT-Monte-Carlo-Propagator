"""Experimental script for testing new functions and features."""
from datetime import datetime
import numpy as np
import pandas as pd
import multiprocessing as mp
from scipy import stats
from astropy import units as u
from matplotlib import pyplot as plt

from config import *
import MC_handler as MC
from main import obj_name, start_epoch, step, duration, n

# Define the constants
MC_data_folder = 'MCWorkspace/Data'
obj_name = '2024 PT5'
n = 500

EH = 1.4714e6

# Load the data
i = 0
data_filename = f"{MC_data_folder}/{obj_name}_ephem{i}.csv"

# Isolate the columns for altitude and C3 energy
altitude_col_name = f"Ast{obj_name.replace(' ', '')} Earth Altitude"
eng_col_name = f"Ast{obj_name.replace(' ', '')} Earth C3Energy"

data_file = pd.read_csv(data_filename, header=0)
print('File read')

altitude_col = data_file[altitude_col_name].values / EH
eng_col = data_file[eng_col_name].values

# time_col = data_file[f"Ast{obj_name.replace(' ', '')} UTCGregorian"].values
time_col = []
for t in data_file[f"Ast{obj_name.replace(' ', '')} UTCGregorian"].values:
	date = datetime.strptime(t[:-2], gmat_date_format)
	time_col.append(date)
time_col = np.array(time_col)

print('Plotting data...')
plt.plot(time_col, altitude_col, label='Altitude')
plt.show()

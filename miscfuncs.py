"""Miscellaneous functions for the project."""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from datetime import datetime
from astropy import units as u
from uncertainties import ufloat
import csv

from config import *

def load_report(reprot_path, csv_savepath):
    """Loads the report file and returns a pandas dataframe. Also saves the data as a csv file.
    """
    with open(reprot_path, 'r') as file:
        input = file.readlines()

    data = []
    first_line = True
    for line in input:
        line = line.strip('\n').split(',')
        line = filter(lambda x: x != '', line)
        line = list(line)

        if first_line:
            line = [i.replace('.', ' ') for i in line]
            first_line = False

        data.append(line)

    header = data[0]
    data = data[1:]
    data = np.array(data).T
    
    data = {header[i]: data[i] for i in range(len(header))}

    data = pd.DataFrame(data)

    if csv_savepath is not None:
        data.to_csv(csv_savepath, index=False)

    return data

def load_state_vectors(savepath = None, elements = None, vectors = None):
    """Loads previously saved data from a JPL Horizon query. The data is saved in a csv file."""
    if savepath is not None:
        data = pd.read_csv(savepath)
        obj_name = savepath.split('/')[-1].split('.')[0]
    
        elements = {
            'datetime_str': data['datetime_str'],
            "e": data['e'],
            "a": data['a'],
            "incl": data['incl'],
            "Omega": data['Omega'],
            "w": data['w'],
            "nu": data['nu']
        }

        vectors = {
            'datetime_str': data['datetime_str'],
            'x': data['x'],
            'y': data['y'],
            'z': data['z'],
            'vx': data['vx'],
            'vy': data['vy'],
            'vz': data['vz'],
        }


    query_time = datetime.strptime(elements['datetime_str'][0][5:], horizon_date_format)
    gmat_init_time = query_time.strftime(gmat_date_format)[:-3]

    # Obtain the elements
    ECC = elements['e'][0]
    SMA = (elements['a'][0] * u.au).to(u.km).value
    INC = elements['incl'][0]
    RAAN = elements['Omega'][0]
    AOP = elements['w'][0]
    TA = elements['nu'][0]

    obj_name = 'Ast' + obj_name.replace(' ', '')

    # Create the GMAT object
    obj_state_kepler = {
        'Name': obj_name,
        'Orbit': {
            'Epoch': gmat_init_time,
            'DateFormat': 'UTCGregorian',
            'CoordSys': 'EarthMJ2000Ec',
            'DisplayStateType': 'Keplerian',
            'ECC': ECC,
            'SMA': SMA,
            'INC': INC,
            'RAAN': RAAN,
            'AOP': AOP,
            'TA': TA,
        }
    }

    obj_state_cartesian = {
        'Name': obj_name,
        'Orbit': {
            'Epoch': gmat_init_time,
            'DateFormat': 'UTCGregorian',
            'CoordSys': 'EarthMJ2000Ec',
            'DisplayStateType': 'Cartesian',
            'X': vectors['x'][0],
            'Y': vectors['y'][0],
            'Z': vectors['z'][0],
            'VX': vectors['vx'][0],
            'VY': vectors['vy'][0],
            'VZ': vectors['vz'][0] 
        }
    }

    return obj_name, obj_state_kepler, obj_state_cartesian

def eccentric_to_true(E, e):
    """Both input and output are in degrees"""
    E = E / 180 * np.pi
    return 2 * arctan(((1+e)/(1-e))**0.5 * tan(E/2)) * 180 / np.pi

def mean_to_eccentric(M, e):
    M = M / 180 * np.pi

    if M > np.pi:
        E = M + e/2
    else:
        E = M - e/2

    count = 0
    E_list = []
    end_next = False
    while count < 1e3:
        f_E = E - e * sin(E) - M
        df_dE = (1 - e * cos(E))
        E = E - (f_E / df_dE)
        
        count += 1
        E_list.append(E)

        if (f_E / df_dE) < 1e-8 and end_next == False:
            end_next = True
            continue
        elif end_next:
            break

    return E * 180 / np.pi

def sin(ufloat_i):
    """Find the uncertainty in sin(i) given the uncertainty in i."""
    if type(ufloat_i) in [float, int, np.float64]:
        return np.sin(ufloat_i)
    
    i = ufloat_i.nominal_value
    sigma_i = ufloat_i.std_dev
    
    return ufloat(np.sin(i), np.abs(np.cos(i) * sigma_i))

def cos(ufloat_i):
    """Find the uncertainty in cos(i) given the uncertainty in i."""
    if type(ufloat_i) in [float, int, np.float64]:
        return np.cos(ufloat_i)
    
    i = ufloat_i.nominal_value
    sigma_i = ufloat_i.std_dev
    
    return ufloat(np.cos(i), np.abs(np.sin(i) * sigma_i))

def tan(ufloat_i):
    """Find the uncertainty in tan(i) given the uncertainty in i."""
    if type(ufloat_i) in [float, int, np.float64]:
        return np.tan(ufloat_i)
    
    i = ufloat_i.nominal_value
    sigma_i = ufloat_i.std_dev
    
    return ufloat(np.tan(i), np.abs(1 / np.cos(i)**2 * sigma_i))

def arctan(ufloat_i):
    """Find the uncertainty in arctan(i) given the uncertainty in i."""
    if type(ufloat_i) in [float, int, np.float64]:
        return np.arctan(ufloat_i)
    
    i = ufloat_i.nominal_value
    sigma_i = ufloat_i.std_dev
    
    return ufloat(np.arctan(i), sigma_i / (1 + i**2))


if __name__ == '__main__':
    pass

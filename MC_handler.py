"""Functions for handling the Monte Carlo simulations."""
import numpy as np
from astroquery.jplsbdb import SBDB
from uncertainties import ufloat
from numpy.random import normal
from scipy.interpolate import CubicSpline
import os
from datetime import datetime
from astropy import units as u
from astropy import constants as cons

from miscfuncs import *
import gmat_handler as gh

from matplotlib import pyplot as plt

def import_uncertainties(obj_name:str):
    """Imports the uncertainties for any small object from the JPL SBDB"""
    data = SBDB.query_async(obj_name).json()

    orbit = data['orbit']
    elements_raw = orbit['elements']
    elements_sigma = {}
    for i in elements_raw:
        {elements_sigma.update({('sigma_'+i['name']): eval(i['sigma'])})}
    
    elements_sigma['sigma_a'] = elements_sigma['sigma_a']# * u.AU.to(u.km)

    e = float(elements_raw[0]['value'])
    M = float(elements_raw[6]['value'])

    e_sigma = elements_sigma['sigma_e']
    M_sigma = elements_sigma['sigma_ma']

    E = mean_to_eccentric(ufloat(M, M_sigma), ufloat(e, e_sigma))
    TA = eccentric_to_true(E, ufloat(e, e_sigma))

    TA_sigma = TA.s
    elements_sigma['sigma_ta'] = TA_sigma

    return elements_sigma

def make_random_elements(elements, sigma, n):
    """Makes random Keplarian orbital elements based on the
        nominal values and uncertainties of the input elements.

        Output is an array of random orbital elements, in the following order:
        SMA, ECC, INC, LAN, AOP, TA

    Args:
        elements (dict): Dictionary of nominal orbital elements.
        sigma (dict): Dictionary of uncertainties in the orbital elements.
        n (int): Number of random elements to generate.

    Returns:
        np.ndarray: Array of random orbital elements.    
    """
    SMA_MC = normal(elements['a'][0], sigma['sigma_a'], n)
    ECC_MC = normal(elements['e'][0], sigma['sigma_e'], n)
    INC_MC = normal(elements['incl'][0], sigma['sigma_i'], n)
    LAN_MC = normal(elements['Omega'][0], sigma['sigma_om'], n)
    AOP_MC = normal(elements['w'][0], sigma['sigma_w'], n)
    TA_MC = normal(elements['nu'][0], sigma['sigma_ta'], n)

    # results = np.array([SMA_MC, ECC_MC, INC_MC, LAN_MC, AOP_MC, TA_MC])

    output = []
    for i in range(n):
        MC_elements = {'a': (SMA_MC[i], 'km'), 'e': (ECC_MC[i], ''),
                       'incl': (INC_MC[i], 'deg'), 'Omega': (LAN_MC[i], 'deg'),
                       'w': (AOP_MC[i], 'deg'), 'nu': (TA_MC[i], 'deg')}
        output.append(MC_elements)

    return output#, results

def MC_thread(iter_num, obj_name, sim_args:tuple, path_args:tuple, orb_args:tuple):
    """Runs a single Monte Carlo simulation.
    
    Args:
        iter_num (int): Iteration number.
        obj_name (str): Name of the object being simulated.
        sim_args (tuple): Arguments for the simulation. (start, step, duration)
        path_args (tuple): Arguments for the output path.
            (template path, report path, script savepath)
        orb_args (tuple): Arguments for the orbital elements. (elements, sigma)
    """
    # Unpack the arguments
    start, step, duration = sim_args
    template_path, report_name, script_savepath, csv_folder = path_args
    elements, vectors = orb_args

    alt_column_name = f"Ast{obj_name.replace(' ', '')} Earth Altitude"
    eng_column_name = f"Ast{obj_name.replace(' ', '')} Earth C3Energy"

    func_start = datetime.now()

    # Load the orbital parameters to script
    script_modules = gh.load_sim_params(template_path, obj_name,
                                        start, step, duration,
                                        report_name, elements, vectors)

    # Run GMAT
    gh.execute_script(script_modules, script_savepath)

    simulation_end = datetime.now()
    # print(f"Finished simulation {iter_num} in {simulation_end - func_start}")

    # return 
    # Read the report
    report_path = gh.GmatInstall.replace('\\', '/') + '/output/' + report_name
    
    if csv_folder is None:
        csv_path = None

    csv_path = f'{csv_folder}/{obj_name}_ephem{iter_num}.csv'

    ephemeris = load_report(report_path, csv_path)

    # Check if the altitude column has any negative values
    if alt_column_name in ephemeris.columns:
        alt_column = ephemeris[alt_column_name].values
        alt_column = alt_column.astype(np.float64)

        if min(alt_column) < (1*u.R_earth).to(u.km).value:
            print(f"Clone {iter_num} ran into Earth. Skipping.")
            os.remove(csv_path)

            ephemeris = None

        alt_column = alt_column + (1 * u.R_earth).to(u.km).value

        ephemeris[alt_column_name] = alt_column

    # Compute the absolute velocity from geocentric energy, then check if the 
    # velocity is negative
    elif eng_column_name in ephemeris.columns:
        eng_column = ephemeris[eng_column_name].values
        eng_column = eng_column.astype(np.float64) * u.km**2/u.s**2
        
        grav_pot = -2 * cons.GM_earth / (alt_column * u.km)
        grav_pot = grav_pot.to(u.km**2/u.s**2)

        if sum((eng_column < grav_pot).astype(int)) > 0:
            print(f"Clone {iter_num} entered a forbidden region. Skipping.")
            os.remove(csv_path)

            ephemeris = None
            return ephemeris

    # Convert the angular coordinates to values within the correct range
    coords_list = ['INC', 'RAAN', 'AOP', 'TA']
    for i in coords_list:
        helio_column_name = f"Ast{obj_name.replace(' ', '')} Heliocentric {i}"

        try:
            temp_column = ephemeris[helio_column_name].values
        except KeyError:
            helio_column_name = f"Ast{obj_name.replace(' ', '')} Sun {i}"
            temp_column = ephemeris[helio_column_name].values

        temp_column = temp_column.astype(np.float64) % (360 - 180 * (i == 'INC'))
        ephemeris[helio_column_name] = temp_column 

        geo_column_name = f"Ast{obj_name.replace(' ', '')} EarthMJ2000Eq {i}"
        
        try:
            temp_column = ephemeris[geo_column_name].values
        except KeyError:
            geo_column_name = f"Ast{obj_name.replace(' ', '')} Earth {i}"
            temp_column = ephemeris[geo_column_name].values

        temp_column = temp_column.astype(np.float64) % (360 - 180 * (i == 'INC'))
        ephemeris[geo_column_name] = temp_column

    # Remove old reports to save space
    if os.path.exists(report_path):
        os.remove(report_path)

    print(f"Clone {iter_num} used {datetime.now() - func_start} to run.")

    return ephemeris

def merge_data(data, column, start_epoch, time_column_name, time_axis):
    """Merges the orbital elements from each simulation into a single 2D array.
    For when the data contains uneven lengths, use interpolation.
    
    Args:
        data (list): List of dataframes from each simulation.
        column (str): Column name to merge.
        start_epoch (str): Start epoch of the simulation.
        time_column_name (str): Name of the time column in the dataframe.
        time_axis (np.ndarray): Time axis for interpolation.

    Returns:
        tuple: Column name and the merged data as a 2D numpy array.
    """
    start = datetime.now()
    col_data = []
    count = 0
    start_time = datetime.strptime(start_epoch + " 00:00:00", "%Y-%m-%d %H:%M:%S")
    for result in data:
        time_column = result[time_column_name].values
        # time_column = [(datetime.strptime(i, gmat_date_format) - 
        #                 start_time).total_seconds()/86400
        #                 for i in time_column]
        time_column = pd.to_datetime(time_column, format=gmat_date_format)
        time_column = (time_column - start_time).total_seconds()/86400

        # Convert the time into days since start of simulation
        interp = CubicSpline(time_column, result[column].values)
        col_data.append(interp(time_axis))

        count += 1

    print(f"Done merging {column} in {datetime.now() - start}.")

    return column, np.array(col_data)


if __name__ == '__main__':
    # Test MC_thread
    template_path = f'{script_folder}/template.script'
    MC_data_folder = 'MCWorkspace/Data'
    
    obj_name = '2006 RH120'
    start_epoch = '2006-01-01'

    elements, _ = gh.query_jpl_horizon(obj_name, start_epoch)
    elements_sigma = import_uncertainties(obj_name)
    elements_MC = make_random_elements(elements, elements_sigma, 10)

    sim_args = (start_epoch, 1 * u.day, 5 * u.year)
    path_args = (template_path,
                 f'Report{obj_name}_0.txt',
                 f'MCWorkspace/Scripts/{obj_name}_0.script',
                 MC_data_folder)
    
    MC_thread(0, obj_name, sim_args, path_args, (elements_MC[0], None))


    
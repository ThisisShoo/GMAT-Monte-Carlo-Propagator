"""The main execution script"""
import numpy as np
from astropy import units as u
from datetime import datetime, timedelta
import multiprocessing as mp
import pandas as pd
from scipy import stats

from config import *

obj_name = '2022 NX1'
start_epoch = '2021-07-01'
step = 1 * u.day
duration = 5 * u.year
n = 500

# Define all the constants here
all_columns = [f"Ast{obj_name.replace(' ', '')} {i.replace('.', ' ')}" for i in default_params]
time_column_name = all_columns[0]
alt_column_name = all_columns[2]
    
columns_toread = ["UTCGregorian", "Earth.C3Energy", "Earth.Altitude", "Earth.SMA",
 "Earth.ECC", "EarthMJ2000Eq.INC", "EarthMJ2000Eq.RAAN", "EarthMJ2000Eq.AOP",
 "Earth.TA", 'Sun.SMA', 'Sun.ECC', 'Heliocentric.INC', 'Heliocentric.RAAN',
 'Heliocentric.AOP', 'Sun.TA']
columns_toread = [f"Ast{obj_name.replace(' ', '')} {i.replace('.', ' ')}" for i in columns_toread]
template_path = f'{script_folder}/template.script'
MC_data_folder = 'MCWorkspace/Data'

if __name__ == '__main__':
    import gmat_handler as gh
    import MC_handler as MC

    # Script starts
    start = datetime.now()

    # Get the conditions
    print('Querying starting conditions...')
    elements, _ = gh.query_jpl_horizon(obj_name, start_epoch)
    elements_sigma = MC.import_uncertainties(obj_name)
    elements_MC = MC.make_random_elements(elements, elements_sigma, n)

    # Run the MC simulation
    sim_args = (start_epoch, step, duration)

    print('Running MC simulation...')
    pool = mp.Pool(processes=mp.cpu_count())

    results = []
    for i in range(n):
        path_args = (template_path, f'Report{obj_name}_{i}.txt',
                       f'MCWorkspace/Scripts/{obj_name}_{i}.script',
                       f'MCWorkspace/Data')
        process = pool.apply_async(MC.MC_thread,
                                   args=(i, obj_name,
                                         sim_args, path_args,
                                         (elements_MC[i], None)))
        results.append(process)

    all_results = []
    for i, process in enumerate(results):
        try:
            result = process.get()
            if result is not None:
                all_results.append(result)
        except Exception as e:
            print(f"Error in simulation {i}: {e}")

    # Compute the mode of the lengths of the results
    lengths = [len(result) for result in all_results]
    expected_length = int(stats.mode(lengths)[0])

    # Close the pool
    pool.close()
    pool.join()

    # Script ends
    mark = datetime.now()

    print(f"Simulations ended in {mark - start}")
    
    time_axis = np.linspace(0, duration.to(u.day).value, expected_length)
    # Make a 2D array for each column, with each 2D array's dimensions be (simulation duration, number of simulations)
    # Use multiprocessing
    pool = mp.Pool(mp.cpu_count())

    # processes = []
    data = {}
    for column in columns_toread:
        if column == time_column_name:
            continue
        # try: 
        #     process = pool.apply_async(MC.merge_data,
        #                            args=(all_results, column, start_epoch, time_column_name, time_axis))
        #     processes.append(process)
        # except Exception as e:
        #     print(f"Error in simulation {i}: {e}")

        col_name, col_data = MC.merge_data(all_results, column, start_epoch,
                                            time_column_name, time_axis)
        data[col_name] = col_data
    
    print("Finished collecting the Monte Carlo results")

    start_time = datetime.strptime(start_epoch + " 00:00:00.00", "%Y-%m-%d %H:%M:%S.%f")
    # time_column = [start_time + timedelta(days=i) for i in time_axis]
    time_column = []
    for i in time_axis:
        time = start_time + timedelta(days=i)
        time = time.strftime(gmat_date_format)
        time_column.append(time)

    # Compute the averages and standard deviations
    ephemeris = {time_column_name: time_column}
    for column in columns_toread:
        if column == time_column_name:
            continue
        elif column == alt_column_name:
            # Check if the altitude is always positive. If not, discard the whole MC clone
            pass
        try:
            ephemeris[column] = np.mean(data[column], axis=0)
        except Exception as e:
            print(data)
        ephemeris[column + '_sigma'] = np.std(data[column], axis=0)

    # Save the results to a CSV
    print("Saving results to CSV")
    df = pd.DataFrame(ephemeris)#, orient='index')
    # df = df.transpose()
    df.to_csv(f"{data_folder}/Full ephemeris {obj_name}.csv", index=False)

    print(f"Script took {datetime.now() - start} to run")

    
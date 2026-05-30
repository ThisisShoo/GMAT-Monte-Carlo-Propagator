"""Analysis functions"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from datetime import datetime
from astropy import constants as cons
import astropy.units as u

from config import *
from main import start_epoch, obj_name

# Load the data
def extract_column(data: pd.DataFrame, column_name: str):
	"""Plots the data in the column against time."""
	column = data[column_name]
	column_sigma = data[column_name + "_sigma"]

	return column, column_sigma

if __name__ == "__main__":
	obj_name = '2024 PT5'

	filepath = f'{data_folder}/Full ephemeris {obj_name}.csv'
	columns = [f"{obj_name} {param.replace('.', ' ')}" for param in default_params]
	time_column_name = f"Ast{obj_name.replace(' ', '')} {default_params[0]}"
	
	# Load the data
	data = pd.read_csv(filepath)

	# time_column = [datetime.strptime(t[:-3], gmat_date_format) for t in data[time_column_name].values]
	time_column = []
	for t in data[time_column_name].values:
		date = datetime.strptime(t[:-3], gmat_date_format)
		# time_elapsed = (date - datetime.strptime(start_epoch, "%Y-%m-%d")).total_seconds() * u.s
		# time_elapsed = time_elapsed.to(u.year).value
		time_column.append(date)
	time_column = np.array(time_column)

	# Plot the heliocentric data
	print('Plotting heliocentric data...')
	fig, ((sma, ecc), (inc, raan), (aop, ta)) = plt.subplots(3, 2, figsize=(9, 9))
	
	sun_sma, sun_sma_sigma = extract_column(data, 
								f"Ast{obj_name.replace(' ', '')} Sun SMA")
	sma.plot(time_column, sun_sma, linestyle='-', color = 'black')#, label='Heliocentric SMA')
	sma.fill_between(time_column, 
					 sun_sma - sun_sma_sigma, sun_sma + sun_sma_sigma,
						color='black', alpha=0.2)
	sma.set_ylabel('Semi-Major Axis [km]')
	sma.plot(time_column, np.ones(len(time_column)) * 1.4959787e8,
		  linestyle='--', color = 'black', label='1 AU')
	# 1.4714e6
	sma.plot(time_column, np.ones(len(time_column)) * 1.4959787e8 + 1.4714e6,
		  linestyle='-.', color = 'black', label='+/- 1 EH')
	sma.plot(time_column, np.ones(len(time_column)) * 1.4959787e8 - 1.4714e6,
		  linestyle='-.', color = 'black')
	sma.plot(time_column, np.ones(len(time_column)) * 1.4959787e8 + 3 * 1.4714e6,
		  linestyle=':', color = 'black', label='+/- 3 EH')
	sma.plot(time_column, np.ones(len(time_column)) * 1.4959787e8 - 3 * 1.4714e6,
		  linestyle=':', color = 'black')
	sma.tick_params(axis='x', labelrotation = 80)
	sma.grid()
	sma.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol = 3)

	sun_ecc, sun_ecc_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} Sun ECC")
	ecc.plot(time_column, sun_ecc, linestyle='-', color = 'black', label='Heliocentric ECC')
	ecc.fill_between(time_column, 
					 sun_ecc - sun_ecc_sigma, sun_ecc + sun_ecc_sigma,
						color='black', alpha=0.2)
	ecc.set_ylabel('Eccentricity')
	ecc.tick_params(axis='x', labelrotation = 80)
	ecc.grid()

	helio_inc, helio_inc_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} Heliocentric INC")
	inc.plot(time_column, helio_inc, linestyle='-', color = 'black', label='Heliocentric INC')
	inc.fill_between(time_column, 
					 helio_inc - helio_inc_sigma, helio_inc + helio_inc_sigma,
						color='black', alpha=0.2)
	inc.set_ylabel('Inclination [deg]')
	inc.tick_params(axis='x', labelrotation = 80)
	inc.grid()

	helio_raan, helio_raan_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} Heliocentric RAAN")
	raan.plot(time_column, helio_raan, linestyle='-', color = 'black', label='Heliocentric RAAN')
	raan.fill_between(time_column, 
					 helio_raan - helio_raan_sigma, helio_raan + helio_raan_sigma,
						color='black', alpha=0.2)
	raan.set_ylabel('RAAN [deg]')
	raan.tick_params(axis='x', labelrotation = 80)
	raan.grid()

	helio_aop, helio_aop_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} Heliocentric AOP")
	aop.plot(time_column, helio_aop, linestyle='-', color = 'black', label='Heliocentric AOP')
	aop.fill_between(time_column, 
					 helio_aop - helio_aop_sigma, helio_aop + helio_aop_sigma,
						color='black', alpha=0.2)
	aop.set_ylabel('AOP [deg]')
	aop.tick_params(axis='x', labelrotation = 80)
	aop.grid()

	sun_ta, sun_ta_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} Sun TA")
	ta.plot(time_column, sun_ta, linestyle='-', color = 'black', label='Heliocentric TA')
	ta.fill_between(time_column, 
					 sun_ta - sun_ta_sigma, sun_ta + sun_ta_sigma,
						color='black', alpha=0.2)
	ta.set_ylabel('True Anomaly [deg]')
	ta.tick_params(axis='x', labelrotation = 80)
	ta.grid()
	
	# fig.suptitle(f"{obj_name} Heliocentric Ephemeris")
	fig.supxlabel('Date')
	# fig.supxlabel(f'Years since {start_epoch}')
	fig.tight_layout()
	fig.savefig(f"{data_folder}/{obj_name} heliocentric plots.png")

	# Plot the geocentric data
	print('Plotting geocentric data...')

	fig, ((c3eng, c3eng_zoom), (dist, dist_zoom), (inc, raan)) = plt.subplots(3, 2, figsize=(7, 7))
	
	c3, c3_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} Earth C3Energy")
	c3eng.plot(time_column, c3, linestyle='-', color = 'black',
			label='Geocentric C3Energy')
	c3eng.fill_between(time_column, 
					 c3 - c3_sigma, c3 + c3_sigma,
						color='black', alpha=0.2)
	c3eng.set_ylabel('C3 Energy [km^2/s^2]')
	c3eng.tick_params(axis='x', labelrotation = 80)
	c3eng.grid()

	c3eng_zoom.plot(time_column, c3, linestyle='-', color = 'black',
			label='Geocentric C3Energy')
	c3eng_zoom.fill_between(time_column, 
					 c3 - c3_sigma, c3 + c3_sigma,
						color='black', alpha=0.2)
	c3eng_zoom.grid()
	c3eng_zoom.set_ylim(-0.5, 0.75)
	c3eng_zoom.hlines(0, time_column[0], time_column[-1],
				   color='black', linestyle='--', label='C3 = 0', linewidth=1)
	c3eng_zoom.tick_params(axis='x', labelrotation = 80)

	earth_dist, earth_dist_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} Earth Altitude")
	earth_dist, earth_dist_sigma = np.array([earth_dist, earth_dist_sigma]) / (1.4714e6)
	dist.plot(time_column, earth_dist, linestyle='-', color = 'black', 
		   label='Geocentric Distance')
	dist.fill_between(time_column, 
					 earth_dist - earth_dist_sigma, earth_dist + earth_dist_sigma,
						color='black', alpha=0.2)
	dist.set_ylabel('Distance from Earth [R_Earth Hill]')
	dist.tick_params(axis='x', labelrotation = 80)
	dist.grid()

	dist_zoom.plot(time_column, earth_dist, linestyle='-', color = 'black')#,
			#label='Geocentric Altitude')
	dist_zoom.fill_between(time_column, 
					 earth_dist - earth_dist_sigma, earth_dist + earth_dist_sigma,
						color='black', alpha=0.2)
	dist_zoom.grid()
	dist_zoom.set_ylim(-1, 8)
	dist_zoom.hlines(1, time_column[0], time_column[-1],
				   color='black', linestyle='--', label='1 EH')
	dist_zoom.hlines(3, time_column[0], time_column[-1],
				   color='black', linestyle=':', label='3 EH')
	dist_zoom.legend(bbox_to_anchor=(1, 1), loc='upper right')
	dist_zoom.tick_params(axis='x', labelrotation = 80)
	
	earth_inc, earth_inc_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} EarthMJ2000Eq INC")
	inc.plot(time_column, earth_inc, linestyle='-', color = 'black',
			label='Geocentric INC')
	inc.fill_between(time_column, 
					 earth_inc - earth_inc_sigma, earth_inc + earth_inc_sigma,
						color='black', alpha=0.2)
	inc.set_ylabel('Inclination [deg]')
	inc.tick_params(axis='x', labelrotation = 80)
	inc.grid()

	earth_raan, earth_raan_sigma = extract_column(data,
								f"Ast{obj_name.replace(' ', '')} EarthMJ2000Eq RAAN")
	raan.plot(time_column, earth_raan, linestyle='-', color = 'black', 
		   label='Geocentric RAAN')
	raan.fill_between(time_column, 
					 earth_raan - earth_raan_sigma, earth_raan + earth_raan_sigma,
						color='black', alpha=0.2)
	raan.set_ylabel('RAAN [deg]')
	raan.tick_params(axis='x', labelrotation = 80)
	raan.grid()

	# fig.suptitle(f"{obj_name} Geocentric Ephemeris")
	fig.supxlabel('Date')
	# fig.supxlabel(f'Years since {start_epoch}')
	fig.tight_layout()
	fig.savefig(f"{data_folder}/{obj_name} geocentric plots.png")
	

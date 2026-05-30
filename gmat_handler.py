"""Utility functions for the project."""
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from datetime import timedelta
from astropy import units as u
from astroquery.jplhorizons import Horizons
import pandas as pd
import os

from load_gmat import gmat, GmatInstall
from load_gmat import gmat_py_simple as gpy
from config import *
from miscfuncs import *

ssl_cert_path = 'ssd.jpl.nasa.gov.crt'

def query_jpl_horizon(obj_name, query_time, savepath = None):
    """Queries JPL Horizon for the orbital elements of a target at a given time.
    Returns the osculating Keplarian orbital elements in geocentric ecliptic frame.

    target: str
        The name of the target object. Must be in the JPL Horizon database.
    query_time: str
        The time at which to query the elements. Must be in the format 'YYYY-MM-DD' or 'DD MMM YYYY'.
    savepath: str
        The path to save the output csv file. If None, the output is not saved.
    """
    try: 
        os.environ['REQUESTS_CA_BUNDLE'] = ssl_cert_path
    except:
        pass
    month_dict = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06',
                    'Jul': '07', 'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
    # Convert the query time to the format for datetime
    query_time = query_time.split()
    if len(query_time) == 3:
        query_time = query_time[2] + '-' + month_dict[query_time[1]] + '-' + query_time[0]
    else:
        query_time = query_time[0]

    query_time = datetime.strptime(query_time, '%Y-%m-%d')

    # Query JPL Horizon
    end_date = query_time + timedelta(days=1)
    epochs = {'start': query_time.strftime('%Y-%m-%d'), 'stop': end_date.strftime('%Y-%m-%d'), 'step': '1d'}
    obj = Horizons(id=obj_name, location='500@10', epochs=epochs)

    # Get the osculating elements
    elements = obj.elements(refplane='ecliptic')
    vector = obj.vectors(refplane='ecliptic')

    if savepath is not None:
        # Combine the elements and vectors and save to the same csv file
        output = elements|vector
        output.write(savepath, format='csv', overwrite=True)

    return elements, vector

template_path = 'scripts/Template.script'

class Command():
    """Parses GMAT Command from a script"""
    def __init__(self, command = None):
        """If the line is a GMAT command, it should follow either of the following formats:
        
        Create <object> <name>;
        GMAT <name>.<property> = <value>;

        Then when the iterations reach the "Mission sequence" block, the following format is used:
        BeginMissionSequence;
        Propagate <PropagatorName>(<ObjectName>) {<SpacecraftName>.<ElapsedSecs> = <value>};
        """
        self.command = command
        split_command = command.split(' ')

        self.type = split_command[1]
        self.name = split_command[2].replace(';', '')

    def __str__(self):
        """Prints all members of the class and variables as a string."""
        output = ''
        for key, value in self.__dict__.items():
            output += f'{key}: {value}\n'
        return output
    
    def export(self):
        """Exports the command as a full script."""
        keys_to_skip = ['name', 'header', 'command', 'type', 'params']

        if self.type == 'MissionSequence':
            script_text = self.header + '\n' + 'BeginMissionSequence;\n'
            script_text += f'Propagate {self.propagator}({self.object})' 
            # script_text += ' {'+f"{self.object}.{self.variable} = {self.duration}"+'};'
            script_text += ' {'
            for key, value in self.variable.items():
                script_text += f"{key} = {value}, "

            script_text = script_text[:-2] + '};\n'

            return script_text

        script_text = self.header + '\n' #+ self.command + '\n'
        script_text += f"Create {self.type} {self.name};\n"

        for key, value in self.__dict__.items():
            if key in keys_to_skip:
                continue
            script_text += f'GMAT {self.name}.{key} = {value};\n'

        return script_text


def load_template(template_path):
    """Loads the GMAT template script and returns the template as a dictionary.
    
    template_path: str
        The path to the GMAT template script.

    Returns:
        dict: A dictionary containing the GMAT template script.
    """
    # Read the script
    with open(template_path, 'r') as f:
        input = f.read().split('\n')

    # Load the template script
    header = ''
    header_line_count = 0
    mission_sequence = False
    module_names = []
    for i, line in enumerate(input):
        # Skip empty lines
        if line == '':
            continue

        # Collect the header lines
        if line.startswith('%-') and header_line_count < 3:
            header += line + '\n'
            header_line_count += 1
            continue
        elif mission_sequence is False:
            header_line_count = 0
            header_temp = header
            header = ''

        split_line = line.split(' ')

        if line == 'BeginMissionSequence;' or mission_sequence:
            mission_sequence = True
            mission_header = header_temp
        
        if split_line[0] == 'Create':
            current_module = Command(line)
            globals()[split_line[1]] = current_module
            current_module.header = header_temp
            module_names.append(split_line[1])
            
        elif split_line[0] == 'GMAT' and current_module.name == split_line[1].split('.')[0]:
            variable = '.'.join(split_line[1].split('.')[1:])
            value = ' '.join(split_line[3:])[:-1]

            setattr(current_module, variable, value)

        elif mission_sequence and split_line[0] == 'Propagate':
            sc_name = globals()['Spacecraft'].name
            prop_command = Command(line)
            prop_command.header = mission_header
            prop_command.type = 'MissionSequence'

            prop_command.propagator = globals()['Propagator'].name
            prop_command.object = sc_name
            prop_command.variable = {
                f'{sc_name}.ElapsedSecs': 0,
                f'{sc_name}.Earth.Altitude': 100
            }

    # Clean up some of the variables
    report_params = globals()['ReportFile'].Add[1:-1].split(', ')
    report_params = ['.'.join(param.split('.')[1:]) for param in report_params]
    globals()['ReportFile'].params = report_params

    template = {
        'Spacecraft': globals()['Spacecraft'],
        'ForceModel': globals()['ForceModel'],
        'Propagator': globals()['Propagator'],
        'CoordinateSystem': globals()['CoordinateSystem'],
        'ReportFile': globals()['ReportFile'],
        'MissionSequence': prop_command
    }

    return template

def load_sim_params(template_path, obj_name, start_epoch, step_size, duration,
                    report_filename=None, elements=None, vector=None, min_distance=100):
    """Loads the simulation parameters into the GMAT template script.

    template_path: str
        The path to the GMAT template script.
    obj_name: str
        The name of the object to be simulated.
    start_epoch: str
        The starting epoch of the simulation.
    step_size: astropy.Quantity
        The step size of the simulation.
    duration: astropy.Quantity
        The duration of the simulation.
    elements: dict
        The osculating Keplerian elements of the object.
    vector: dict
        The state vector of the object.

    Returns:
        dict: A dictionary containing the GMAT script with the updated parameters.
    """
    template = load_template(template_path)

    if elements is None and vector is None:
        raise ValueError('Either elements or vector must be provided')
    
    Spacecraft = template['Spacecraft']
    Propagator = template['Propagator']
    MissionSequence = template['MissionSequence']
    ReportFile = template['ReportFile']

    # Convert the time format to GMAT format (UTC Gregorian)
    start_epoch = datetime.strptime(start_epoch, '%Y-%m-%d')
    start_epoch = start_epoch.strftime(gmat_date_format)[:-3]

    # Load the simulation parameters
    Spacecraft.name = 'Ast' + (obj_name.replace(' ', ''))
    Spacecraft.Epoch = start_epoch

    if elements is not None:
        Spacecraft.SMA = str(elements['a'][0] * u.AU.to(u.km))
        Spacecraft.ECC = str(elements['e'][0])
        Spacecraft.INC = str(elements['incl'][0])
        Spacecraft.RAAN = str(elements['Omega'][0])
        Spacecraft.AOP = str(elements['w'][0])
        Spacecraft.TA = str(elements['nu'][0])
    
    elif vector is not None:
        Spacecraft.X = str(vector['x'][0] * u.AU.to(u.km))
        Spacecraft.Y = str(vector['y'][0] * u.AU.to(u.km))
        Spacecraft.Z = str(vector['z'][0] * u.AU.to(u.km))
        Spacecraft.VX = str(vector['vx'][0] * u.AU.to(u.km) / u.day.to(u.s))
        Spacecraft.VY = str(vector['vy'][0] * u.AU.to(u.km) / u.day.to(u.s))
        Spacecraft.VZ = str(vector['vz'][0] * u.AU.to(u.km) / u.day.to(u.s))
    
    Propagator.InitialStepSize = str(step_size.to(u.s).value)

    MissionSequence.object = Spacecraft.name
    MissionSequence.variable = {
        f'{Spacecraft.name}.ElapsedSecs': str(duration.to(u.s).value),
        f'{Spacecraft.name}.Earth.Altitude': str(min_distance)
    }

    param_list = [f"{Spacecraft.name}.{param}" for param in ReportFile.params]
    ReportFile.Add = '{' + ', '.join(param_list) + '}'

    if type(report_filename) is not str:
        report_filename = f"TempReport.txt"

    ReportFile.Filename = report_filename

    # Gather all the modified values
    output = template
    output['Spacecraft'] = Spacecraft
    output['Propagator'] = Propagator
    output['MissionSequence'] = MissionSequence

    return output

def execute_script(script_modules, script_savepath):
    """Assembles the GMAT script from the script modules."""
    time_now = datetime.today().strftime(horizon_date_format.split('.')[0])

    script = f'%General Mission Analysis Tool(GMAT) Script \n%Created: {time_now}\n \n'
    for _, value in script_modules.items():
        script += value.export() + '\n'

    with open(script_savepath, 'w') as f:
        f.write(script)

    gmat_global = gmat.GmatGlobal.Instance()

    gmat_global.SetWriteParameterInfo(True)
    gmat_global.SetWriteFilePathInfo(False)
    gmat_global.SetCommandEchoMode(True)

    # gmat.UseLogFile(f"gmat_log{iter}.txt")
    # gmat.EchoLogFile()
    # print('Echoing GMAT log file to terminal\n')

    gmat.LoadScript(script_savepath)
    gmat.RunScript()

    os.remove(script_savepath)

if __name__ == '__main__':
    test = load_template(template_path)
    print(test['MissionSequence'].export())


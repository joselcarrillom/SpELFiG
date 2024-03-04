#!/usr/bin/env python3

import os, sys, re
import time
import pandas as pd

spectype = sys.argv[1]

if len(sys.argv) < 2:
    print('Please provide the path of the script to import.')
    sys.exit(1)

# &&&& ································ USER ATTENTION HERE ··································· &&&&
# The following function defines the pattern of the ID's of the objects. Modify the pattern
# variable as necessary in accordance of the dataset you are providing to fit.

def find_ID(pattern, file):
    match = re.search(pattern, file)
    if match:
        return match.group(0)

ID_pattern = r'manga-(\d+)-(\d+)' # Change as needed


# Produce the list of files:
basepath = os.getcwd()

if spectype == 'S':
    input_directory = basepath+'/input_data/Spectra'
    exe_path = basepath+'/SpELFiG_exe.py'
elif spectype == 'C':
    input_directory = basepath+'/input_data/Cubes'
    exe_path = basepath+'/SpELFiG_exe_cubes.py'

files_list = []
IDS = []
for root, dirs, files in os.walk(input_directory):
    for file in files:
        file_path = os.path.join(root, file)
        ID_OBJ = find_ID(ID_pattern, file)
        files_list.append(file_path)
        IDS.append(ID_OBJ)


# Path to the existing script (Spelfic_V8_exe) and the configuration script
config_path = basepath+'/SpELFiG_Config_Temp.py'

if spectype == 'S': specdir = 'Spectra'
elif spectype == 'C': specdir = 'Cubes'

output_config_path = basepath+f'/Configs/{specdir}'
timefile_path = output_config_path+'/Execution_Times.txt'
# Specify the path to save modified configuration files:


if not os.path.exists(output_config_path):
    os.makedirs(output_config_path)

# Producing and opening the timefile log:
timefile = open(timefile_path, 'w')

for i, obj in enumerate(IDS):
    input_file = files_list[i]
    # Create a new configuration script for each galaxy
    with open(config_path, 'r') as config_file:
        config_content = config_file.read()

    config_content = config_content.replace("ID_OBJ = '[ID_OBJ]'", f"ID_OBJ = '{obj}'")
    config_content = config_content.replace("INPUT_FILE = '[INPUT_FILE]'", f"INPUT_FILE = '{input_file}'")
    config_content = config_content.replace("spectype = '[SPECTYPE_I]'", f"spectype = '{spectype}'")

    # Save the modified configuration script with a specific name
    output_config_name = f'Spelfic_Config_{obj}.py'
    output_config_file_path = os.path.join(output_config_path, output_config_name)

    with open(output_config_file_path, 'w') as output_config_file:
        output_config_file.write(config_content)

    # Execute the existing script with the modified configuration
    print(f'>>>> ······· >>>> EXECUTING SpELFiG FOR GALAXY: {obj} <<<< ······ <<<<')
    start_time = time.time()

    try:
        os.system(f'python {exe_path} {output_config_file_path}')
    except Exception as e:
        print(f'An execution error ocurred: {e}')
    # Record the end time
    end_time = time.time()
    # Calculate the elapsed time in seconds
    time_seconds = end_time - start_time
    time_minutes = time_seconds / 60  # 60 seconds in a minute
    time_hours = time_minutes / 60  # 60 minutes in an hour
    line_to_write = f"For galaxy {obj}, execution time was: {time_minutes} minutes\n"
    timefile.write(line_to_write)

timefile.close()

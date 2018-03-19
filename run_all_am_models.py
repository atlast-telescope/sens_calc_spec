import os

# delete the existing files
os.system('rm am_runs/*')


# run each am configuration
os.system('/Users/sbryan/Documents/am-9.0/am am_models/ALMA_JJA_50.amc 0 GHz 2550 GHz 10 MHz 30 deg 1.0 > am_runs/JJA_25.dat')
os.system('/Users/sbryan/Documents/am-9.0/am am_models/ALMA_JJA_50.amc 0 GHz 2550 GHz 10 MHz 30 deg 1.02 > am_runs/JJA_25_1p02.dat')
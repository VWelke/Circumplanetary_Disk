# CASA_set_obs_C10.py

##############################################
### SIMOBSERVE() + SIMANALYZE() AUTOMATION ###
### Luke Keyte 05/23                       ###
##############################################

import os

# Define the data dictionary
data = {
    'Sky coordinates': ['14h08m10s'] * 9,
    'Polarisation': ['Dual'] * 9,
    'Bandwidth per Polarization': ['7.5 GHz'] * 6 + ['15 GHz'] * 3,
    'Column Density': ['0.913mm (3rd Octile)'] * 9,
    'Resolution_C10 (arcsecond)': [0.1100, 0.0420, 0.0280, 0.0230, 0.0180, 0.0120, 0.0091, 0.0065, 0.0048],
    'Resolution_C9 (arcsecond)': [0.1400, 0.0570, 0.0380, 0.0310, 0.0250, 0.0170, 0.0120, 0.0088, 0.0066],
    'Resolution_C8 (arcsecond)': [0.240, 0.096, 0.064, 0.052, 0.042, 0.028, 0.021, 0.015, 0.011],
    'Resolution_C7 (arcsecond)': [0.530, 0.210, 0.140, 0.110, 0.092, 0.061, 0.046, 0.033, 0.024],
    'Band': [1, 3, 4, 5, 6, 7, 8, 9, 10],
    "Wavelength (micron)": [7494.81,2997.92,1998.62,1620.50,1304.45,868.96,651.72,461.22,344.59],
    'Wavelength (GHz)': [40, 100, 150, 185, 230, 345, 460, 650, 870],
    'beam_ALMA_C8 (arcsec^2)': [0.065266, 0.010443, 0.004641, 0.003064, 0.001999, 0.000888, 0.000500, 0.000255, 0.000137],
    'beam_ALMA_C7 (arcsec^2)': [0.318285, 0.049969, 0.022209, 0.013710, 0.009590, 0.004216, 0.002398, 0.001234, 0.000653],
    'Noise_C10 (Jy/pixel)': [7.703e-07, 3.689e-06, 6.300e-06, 9.125e-06, 9.573e-06, 1.926e-05, 2.999e-05, 3.040e-05, 2.864e-05],
    'Noise_C9 (Jy/pixel)': [1.248e-06, 6.795e-06, 1.160e-05, 1.658e-05, 1.847e-05, 3.866e-05, 5.214e-05, 5.571e-05, 5.414e-05],
    'Noise_C8 (Jy/pixel)': [3.667e-06, 1.927e-05, 3.291e-05, 4.664e-05, 5.212e-05, 1.049e-04, 1.597e-04, 1.619e-04, 1.504e-04],
    'Noise_C7 (Jy/pixel)': [1.788e-05, 9.223e-05, 1.575e-04, 2.087e-04, 2.501e-04, 4.978e-04, 7.662e-04, 7.835e-04, 7.159e-04],
    'Time_C10 (s)': [265880, 17000, 6800, 339900, 5140, 5100, 25320, 348350, 2274680],
    'Time_C9 (s)': [101300, 5010, 2010, 102960, 1380, 1270, 8380, 103730, 636550],
    'Time_C8 (s)': [11740, 630, 250, 13020, 180, 180, 900, 12290, 82490],
    'Time_C7 (s)': [500, 30, 20, 650, 10, 10, 40, 530, 3650]
}

# Define the base project name and output directory
base_project_name = 'Cycle11_C10'
output_dir = 'observations_C10'

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Loop through each index to create observations for Configuration C10
for i in range(len(data['Band'])):
    project_name = f"{base_project_name}_{data['Time_C10 (s)'][i]}_Band{data['Band'][i]}"
    output_fits = f"{project_name}.fits"
    project_dir = os.path.join(output_dir, project_name)
    
    # Create the project directory if it doesn't exist
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)

    ###############################
    ## OBSERVE EXECUTION BLOCK 1 ##
    ###############################

    default('simobserve')

    # Define the skymodel file based on the band frequency
    skymodel = f'RADMC3d_fits/CPD_PPD_PDS70_{data["Wavelength (micron)"][i]}.fits'  # Model FITS file for each band
    indirection = 'J2000 14h08m10.11s -41d23m52.95s'  # PDS 70 location from NASA exoplanet archive
    incell = '0.02arcsec'  # Model image pixel size (in arcsec)
    inbright = ''
    incenter = f"{data['Wavelength (GHz)'][i]} GHz"  # Central frequency of observation
    mapsize = '2arcsec'  # Total FITS image width
    integration = '30s'  # Reasonable time for each pointing
    obsmode = 'int'
    antennalist = 'alma.cycle11.10.cfg'   # Configuration
    totaltime = f"{data['Time_C10 (s)'][i]}s"  # Total integration time
    refdate = '2025/01/01'
    hourangle = 'transit'  # 04:55
    graphics = 'none'

    simobserve(
        project=project_name,
        skymodel=skymodel,
        indirection=indirection,
        incell=incell,
        inbright=inbright,
        incenter=incenter,
        integration=integration,
        obsmode=obsmode,
        antennalist=antennalist,
        totaltime=totaltime,
        refdate=refdate,
        graphics=graphics,
        hourangle=hourangle,
        overwrite=True
    )

    print(f'SIMOBSERVE COMPLETE for {project_name}')

    ################
    ## SIMANALYZE ##
    ################

    default("simanalyze")

    simanalyze(
        project=project_name,
        vis=f'{project_name}/{project_name}.alma.cycle11.C10.noisy.ms',
        cell='0.02arcsec',
        analyze=True,
        overwrite=True
    )

    print(f'SIMANALYZE COMPLETE for {project_name}')

    #################
    ## EXPORT FITS ##
    #################

    exportfits(
        imagename=f'{project_name}/{project_name}.alma.cycle11.C10.noisy.image',
        fitsimage=os.path.join(project_dir, output_fits)
    )

    print(f'*** FINISHED {project_name} ***')
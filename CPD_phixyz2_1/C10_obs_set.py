# Change output fits


##############################################
### SIMOBSERVE() + SIMANALYZE() AUTOMATION ###
### Luke Keyte 05/23                       ###
##############################################




# Inputs for each observation
frequencies_GHz = [40,100,150,185,230,345,460,650,870]
fig  = [1,3,4,5,6,7,8,9,10]  # Figure number for each frequency


# Loop through each observation setup
for i in range(len(frequencies_GHz)):
    freq       = frequencies_GHz[i]

    project = f'C10_{fig[i]}_36000_pwv1'  # Project and folder name
    
    output_fits = f'{project}.fits'
    skymodel = f'CPD_PPD_PDS70_{fig[i]}.fits'
    incenter = f'{freq} GHz'
    vis = f'{project}/{project}.alma.cycle11.10.noisy.ms'
    imagename = f'{project}/{project}.alma.cycle11.10.noisy.image'

    #####################
    ## SIMOBSERVE BLOCK ##
    #####################

    default('simobserve')

    simobserve(
        project=project,
        skymodel=skymodel,
        indirection='J2000 14h08m10.11s -41d23m52.95s',
        incell='0.002arcsec',
        inbright='',
        incenter=incenter,
        integration='30s',
        obsmode='int',
        antennalist='alma.cycle11.10.cfg',
        totaltime='36000s',
        refdate='2025/01/01',
        thermalnoise='tsys-atm',
        user_pwv=0.472,  # 0.913 for 3rd octile as needed, 0.472 is 1st octile
        graphics='none',
        hourangle='transit',
        overwrite=True
    )

    print(f'SIMOBSERVE COMPLETE for {project}')


    #####################
    ## SIMANALYZE BLOCK ##
    #####################

    default("simanalyze")

    simanalyze(
        project=project,
        vis=vis,
        cell='0.002arcsec',
        analyze=True,
        overwrite=True
    )

    print(f'SIMANALYZE COMPLETE for {project}')


    ##################
    ## EXPORT TO FITS ##
    ##################

    exportfits(imagename=imagename, fitsimage=output_fits)
    print(f'*** EXPORT COMPLETE for {output_fits} ***\n')

print('*** ALL SIMULATIONS COMPLETE ***')
# Change output fits


##############################################
### SIMOBSERVE() + SIMANALYZE() AUTOMATION ###
### Luke Keyte 05/23                       ###
##############################################

# Name of you output FITS file
# Ref the format: output_fits = 'output87_5000s_455_config3_4000.fits'
output_fits = 'C10_6.fits'                                         #######################CHNAGE###########################
###############################
## OBSERVE EXECUTION BLOCK 1 ##
###############################

# 'Time_C10 (s)':  5100, 25320, 348350, 2274680]

# "Wavelength (micron)": [868.96,651.72,461.22,344.59],
# 'Wavelength (GHz)': [345, 460, 650, 870],

##########################################################################################################

default('simobserve')

project = 'C10_6'                              #######################CHANGE####################################
  

skymodel     = 'RADMC3d_fits/CPD_PPD_PDS70_868.96.fits'   # Model FITS file    ##############CHANGE##########################
indirection  = 'J2000 14h08m10.11s -41d23m52.95s'    # PDS 70 location from NASA exoplanet archive                                # RA and DEC coords
incell       = '0.02arcsec'                        # Model image pixel size (in arcsec)
inbright     = ''
incenter     = '345 GHz'    # Central frequency of observation #Band 5 ############################CHANGE#######################################
mapsize      = '2arcsec'          # Total FITS image width
integration  = '30s'  # reasonable time for each pointing
obsmode      = 'int'
#antennalist  = 'alma;0.01arcsec'  # Just set configuration
antennalist  = 'alma.cycle11.10.cfg' # C-10 config
totaltime    = '5100s'          ##################CHANGE###################
refdate      = '2025/01/01'
hourangle    = 'transit'  #04:55
graphics     = 'none'
 
simobserve(
	project=project,
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
	overwrite=True)

print('SIMOBSERVE COMPLETE')


################
## SIMANALYZE ##
################

default("simanalyze")

project = 'C10_6'    #########################CHANGE####################

vis     = 'C10_6/C10_6.alma.cycle11.10.noisy.ms'    #######################CHANGE######################
cell    = '0.02arcsec'  
analyze = True

simanalyze(project=project,
	vis = vis,
	cell=cell,
	analyze=True,
	overwrite=True)

print('SIMANALYZE COMPLETE')



#################
## EXPORT FITS ##
#################

exportfits(imagename='C10_6/C10_6.alma.cycle11.10.noisy.image', fitsimage=output_fits)         ################CHANGE####################

print('*** FINISHED ***')


#exportfits(imagename='Cycle11_9000s_C10_Band5/Cycle11_9000s_C10_Band5.alma.cycle11.10.noisy.image', fitsimage='Cycle11_9000s_C10_Band5_2.fits')
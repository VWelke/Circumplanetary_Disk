# Change output fits


##############################################
### SIMOBSERVE() + SIMANALYZE() AUTOMATION ###
### Luke Keyte 05/23                       ###
##############################################

# Name of you output FITS file
# Ref the format: output_fits = 'output87_5000s_455_config3_4000.fits'
output_fits = 'C10_1_36000.fits'   # (13456789)                                     #######################CHNAGE###########################
###############################
## OBSERVE EXECUTION BLOCK 1 ##
###############################
#timeC10 (first Octile)  [x , x, x, 67968s, ]
#'Time_C10 (s)': [265880, 17000,  6800,   339900,  5140,   5100,  25320,    348350,  2274680]
#  'Wavelength (GHz)': [40, 100, 150, 185, 230, 345, 460, 650, 870],
# 'Wavelength (micron)': [7494.81, 2997.92,1998.62, 1620.50,  1304.45,   868.96, 651.72, 461.22,344.59],
##########################################################################################################

default('simobserve')

project = 'C10_1_36000'                              #######################CHANGE####################################
  

skymodel     = 'CPD_PPD_PDS70_7494.81.fits'   # Model FITS file    ##############CHANGE##########################
indirection  = 'J2000 14h08m10.11s -41d23m52.95s'    # PDS 70 location from NASA exoplanet archive                                # RA and DEC coords
incell       = '0.002arcsec'                        # Model image pixel size (in arcsec)
inbright     = ''
incenter     = '40 GHz'    # Central frequency of observation #Band 5 ############################CHANGE#######################################
mapsize      = '2arcsec'          # Total FITS image width
integration  = '30s'  # reasonable time for each pointing
obsmode      = 'int'
#antennalist  = 'alma;0.01arcsec'  # Just set configuration
antennalist  = 'alma.cycle11.10.cfg' # C-10 config
totaltime    = '36000s'          ##################CHANGE###################
refdate      = '2025/01/01'
thermalnoise='tsys-atm'  # Enables system temperature and atmospheric noise modeling
user_pwv= 0.913 # 0.913 mm for 3rd octile, 0.472mm for 1st octile
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
    thermalnoise=thermalnoise,
    user_pwv=user_pwv,
	graphics=graphics,
	hourangle=hourangle,
	overwrite=True)

print('SIMOBSERVE COMPLETE')


################
## SIMANALYZE ##
################

default("simanalyze")

project = 'C10_1_36000'    #########################CHANGE####################

vis     = 'C10_1_36000/C10_1_36000.alma.cycle11.10.noisy.ms'    #######################CHANGE######################
cell    = '0.002arcsec'  
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

exportfits(imagename='C10_1_36000/C10_1_36000.alma.cycle11.10.noisy.image', fitsimage=output_fits)         ################CHANGE####################

print('*** FINISHED ***')


#exportfits(imagename='Cycle11_9000s_C10_Band5/Cycle11_9000s_C10_Band5.alma.cycle11.10.noisy.image', fitsimage='Cycle11_9000s_C10_Band5_2.fits')
##############################################
### SIMOBSERVE() + SIMANALYZE() AUTOMATION ###
### Luke Keyte 05/23                       ###
##############################################

# Name of you output FITS file
# Ref the format: output_fits = 'output87_5000s_455_config3_4000.fits'
output_fits = 'Cycle11_9000s_C10_Band5.fits'
###############################
## OBSERVE EXECUTION BLOCK 1 ##
###############################

default('simobserve')

project = 'Cycle11_9000s_C10_Band5'


skymodel     = 'RADMC3d_fits/CPD_PPD_PDS70_1621.fits'   # Model FITS file
indirection  = 'J2000 14h08m10.11s -41d23m52.95s'    # PDS 70 location from NASA exoplanet archive                                # RA and DEC coords
incell       = '0.02arcsec'                        # Model image pixel size (in arcsec)
inbright     = ''
incenter     = '185 GHz'    # Central frequency of observation #Band 5
mapsize      = '2arcsec'          # Total FITS image width
integration  = '30s'  # reasonable time for each pointing
obsmode      = 'int'
#antennalist  = 'alma;0.01arcsec'  # Just set configuration
antennalist  = 'alma.cycle11.10.cfg' # C-10 config
totaltime    = '24120s'  #6.7 hours
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

project = 'Cycle11_9000s_C10_Band5'

vis     = 'Cycle11_9000s_C10_Band5/Cycle11_9000s_C10_Band5.alma.cycle11.10.noisy.ms'
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

exportfits(imagename='Cycle11_9000s_C10_Band5/Cycle11_9000s_C10_Band5.alma.cycle11.10.noisy.image', fitsimage=output_fits)

print('*** FINISHED ***')


exportfits(imagename='Cycle11_9000s_C10_Band5/Cycle11_9000s_C10_Band5.alma.cycle11.10.noisy.image', fitsimage='Cycle11_9000s_C10_Band5_2.fits')
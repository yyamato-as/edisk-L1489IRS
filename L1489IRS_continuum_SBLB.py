"""
eDisk data reduction script
This script was written for CASA 6.1.1/6.2
Originally derived from DSHARP reduction scripts


Datasets calibrated (in order of date observed):
SB1: 2013.1.01086.S (20-sep-2015: uid___A002_Xaa5cf7_X5fc1.ms.split.cal) **archival data** 
	 also covered by DDT project - not yet observed (data deliverly will be ~December?)
 

LB1: 2019.1.00261.L (20-Aug-2021: uid___A002_Xef6d27_X1023.ms)
LB2: 2019.1.00261.L (20-Aug-2021: uid___A002_Xef6d27_X168c.ms)
     

reducer: Yoshihide Yamato
"""

### When you restart this script from a certain point, you need to do:
### - import stuffs
### - setting up the data_param dictionary
### CAUTION: When you interrupt this script at a certain point, do not forget to dump the data_param dictionary to a pickle 

### Import statements
#sys.path.append('/home/casa/contrib/AIV/science/analysis_scripts/')
import analysisUtils as au
import analysisUtils as aU
import string
import os
import glob
import numpy as np
import sys
import pickle
execfile('/home/yamato/Project/edisk_analysis/reduction_utils3.py', globals())

#############################################################################################
### have to go to the directry of WD_path to execute the script below by copy and pasting ###
#############################################################################################
WD_path = '/raid/work/yamato/eDisk_data/L1489IRS/'
os.system('cd ' + WD_path)


###############################################################
################ SETUP/METADATA SPECIFICATION #################
################ USERS NEED TO SET STUFF HERE #################
###############################################################

### Use MPI CASA for faster imaging 
### start with
### > ~/Application/CASA/casa-6.2.0-124-eDisk/bin/mpicasa -n 8 ~/Application/CASA/casa-6.2.0-124-eDisk/bin/casa
### I should set the alias to .bashrc...
parallel=True  

### if True, can run script non-interactively if later parameters properly set
skip_plots = True

### start this script from the top or restart from a certain point?
restart = False

### Add field names (corresponding to the field in the MS) here and prefix for 
### filenameing (can be different but try to keep same)
### Only make different if, for example, the field name has a space
# SB_field = 'L1489_IRS'
field   = 'L1489IRS' # different field names between LB and SB since SB is archival data
prefix  = 'L1489IRS' # adopt LB field name

### always include trailing slashes!!
SB_path = WD_path+'ALMA_pipeline_reduced_data/SB/'
LB_path = WD_path+'ALMA_pipeline_reduced_data/LB/'

### scales for multi-scale clean
SB_scales = [0, 5] #[0, 5, 10, 20]
LB_scales = [0, 5, 30]  #[0, 5, 30, 100, 200]

### Add additional dictionary entries if need, i.e., SB2, SB3, LB1, LB2, etc. for each execution
### Note that C18O and 13CO have different spws in the DDT vis LP os the spw ordering
### is different for data that were originally part of the DDT than the LP
### DDT 2019.A.00034.S SB data need 'spws': '25,31,29,27,33,35,37'
### LP  2019.1.00261.L SB data need 'spws': '25,27,29,31,33,35,37'
pl_data_params={'SB1': {'vis': SB_path+'uid___A002_Xaa5cf7_X5fc1.ms.split.cal',
                        'spws': '0,1,2,3,4,5'},
				'LB1': {'vis': LB_path+'uid___A002_Xef6d27_X1023.ms',
                        'spws': '25,27,29,31,33,35,37'},
				'LB2': {'vis': LB_path+'uid___A002_Xef6d27_X168c.ms',
                        'spws': '25,27,29,31,33,35,37'},
               }

### first we need to modify the field name of SB data to be consistent with LB's





### Dictionary defining necessary metadata for each execution
### SiO at 217.10498e9 excluded because of non-detection
### Only bother specifying simple species that are likely present in all datasets
### Hot corino lines (or others) will get taken care of by using the cont.dat

if not restart:
	### may able to make a script to automatically produce these dictionaries?
	data_params = {'SB1': {'vis' : WD_path+prefix+'_SB1.ms',
						   'name' : 'SB1',
						   'field': 'L1489_IRS',
						   'line_spws': np.array([0,1,2,3,4]), # line SPWs, get from listobs
						   'line_freqs': np.array([219.94944200e9,219.56035410e9,220.39868420e9,
												   230.538e9,231.3218e9]), #restfreqs for targeted lines
						   'line_names': ['SO','C18O','13CO','12CO','N2D+'], 
						   'flagrange': np.array([[-5.5,14.5],[-5.5,14.5],[-5.5,14.5],[-5.5,14.5],
												  [-5.5,14.5]]),
						   'orig_spw_map': {0:0, 1:1, 2:2, 3:3, 4:4, 5:5},  # mapping of old spws to new spws (needed for cont.dat to work)
						   'cont_spws':  np.array([0,1,2,3,4,5]),  #spws to use for continuum
						   'cont_avg_width':  np.array([960,960,480,480,480,2]), #n channels to average; approximately aiming for 30 MHz channels
						   'phasecenter': '',
						   'timerange': '2015/09/20/06:50:00~2015/09/20/08:50:00',
# 						   'contdotdat' : SB_path + 'cont.dat' #no cont.dat for SB since it's archival data 
						  },
				   'LB1': {'vis' : WD_path+prefix+'_LB1.ms',
						   'name' : 'LB1',
						   'field': 'L1489IRS',
						   'line_spws': np.array([0,1,2,3,4,6,4,4,4,4,4]), # line SPWs, get from listobs
						   'line_freqs': np.array([218.76006600e9,220.39868420e9,219.94944200e9,219.56035410e9,
												   217.82215e9,230.538e9,217.94005e9,218.16044e9,217.2386e9,
												   218.22219200e9,218.47563200e9]), #restfreqs
						   'line_names': ['H2CO','13CO','SO','C18O','c-C3H2','12CO','c-C3H2','c-C3H2','DCN','H2CO','H2CO'], #restfreqs
						   'flagrange': np.array([[-5.5,14.5],[-5.5,14.5],[-5.5,14.5],[-5.5,14.5],
												  [-5.5,14.5],[-5.5,14.5],[-5.5,14.5],[-5.5,14.5],
												  [-5.5,14.5],[-5.5,14.5],[-5.5,14.5]]),
						   'orig_spw_map': {25:0, 27:1, 29:2, 31:3, 33:4, 35:5, 37:6},  # mapping of old spws to new spws (needed for cont.dat to work)
						   'cont_spws':  np.array([0,1,2,3,4,5,6]),  #spws to use for continuum
						   'cont_avg_width':  np.array([480,480,480,480,60,60,60]), #n channels to average; approximately aiming for 30 MHz channels
						   'phasecenter': '',
						   'timerange': '2021/08/21/09:30:00~2021/08/21/11:00:00',
						   'contdotdat' : LB_path + 'cont.dat'
						  },
				   'LB2': {'vis' : WD_path+prefix+'_LB2.ms',
						   'name' : 'LB2',
						   'field': 'L1489IRS',
						   'line_spws': np.array([0,1,2,3,4,6,4,4,4,4,4]), # line SPWs, get from listobs
						   'line_freqs': np.array([218.76006600e9,220.39868420e9,219.94944200e9,219.56035410e9,
												   217.82215e9,230.538e9,217.94005e9,218.16044e9,217.2386e9,
												   218.22219200e9,218.47563200e9]), #restfreqs
						   'line_names': ['H2CO','13CO','SO','C18O','c-C3H2','12CO','c-C3H2','c-C3H2','DCN','H2CO','H2CO'], #restfreqs
						   'flagrange': np.array([[-5.5,14.5],[-5.5,14.5],[-5.5,14.5],[-5.5,14.5],
												  [-5.5,14.5],[-5.5,14.5],[-5.5,14.5],[-5.5,14.5],
												  [-5.5,14.5],[-5.5,14.5],[-5.5,14.5]]),
						   'orig_spw_map': {25:0, 27:1, 29:2, 31:3, 33:4, 35:5, 37:6},  # mapping of old spws to new spws (needed for cont.dat to work)
						   'cont_spws':  np.array([0,1,2,3,4,5,6]),  #spws to use for continuum
						   'cont_avg_width':  np.array([480,480,480,480,60,60,60]), #n channels to average; approximately aiming for 30 MHz channels
						   'phasecenter': '',
						   'timerange': '2021/08/21/07:40:00~2021/08/21/09:20:00',
						   'contdotdat' : LB_path+'cont.dat'
						  },
				   }

	
else:
	with open(prefix+'.pickle', 'rb') as handle:
		data_params = pickle.load(handle)
		
		
		
### Flag range corresponds to velocity range in each spw that should be flagged. 
### Velocity range should correspond to 
### approximate width of the line contamination

#save data params to a pickle
with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)

###############################################################
#################### DATA PREPARATION #########################
###############################################################

### split out each pipeline-calibrated dataset into an MS only containing the target data 
for i in pl_data_params.keys():
	if os.path.exists(prefix+'_'+i+'.ms'):
		flagmanager(vis=prefix+'_'+i+'.ms', mode="restore", \
                  versionname="starting_flags")
	else:
		if i == 'SB1':
			datacolumn = 'data'
		else:
			datacolumn = 'corrected'
		split(vis=pl_data_params[i]['vis'],outputvis=prefix+'_'+i+'.ms',spw=pl_data_params[i]['spws'],field=data_params[i]['field'],datacolumn=datacolumn)
		
### Backup the the flagging state at start of reduction
for i in data_params.keys():
	if not os.path.exists(data_params[i]['vis'] + ".flagversions/flags.starting_flags"):
		flagmanager(vis=data_params[i]['vis'], mode = 'save', versionname = 'starting_flags', comment = 'Flag　states at start of reduction')
		
### Inspect data in each spw for each dataset
#### OPTIONAL #####
if not skip_plots:
    for i in data_params.keys():
        plotms(vis=data_params[i]['vis'], xaxis='frequency', yaxis='amplitude', 
               field=data_params[i]['field'], ydatacolumn='data', 
               avgtime='1e8', avgscan=True, avgbaseline=True, iteraxis='spw',
               transform=True,freqframe='LSRK')
        input("Press Enter key to advance to next MS/Caltable...")
#### END OPTIONAL ###

### Flag spectral regions around lines and do spectral averaging to make a smaller continuum MS 
for i in data_params.keys():      
    flagchannels_string = get_flagchannels(data_params[i], prefix)
    s=' '  # work around for Python 3 port of following string generating for loops
    print(i) 
    avg_cont(data_params[i], prefix, flagchannels=flagchannels_string,contspws=s.join(str(elem) for elem in data_params[i]['cont_spws'].tolist()).replace(' ',','),width_array=data_params[i]['cont_avg_width'])
    data_params[i]['vis_avg']=prefix+'_'+i+'_initcont.ms'
		

###############################################################
############## INITIAL IMAGING FOR ALIGNMENT ##################
###############################################################


### Image each dataset individually to get source position in each image
### Images are saved in the format prefix+'_name_initcont_exec#.ms'
outertaper='2000klambda' # taper if necessary to align using larger-scale uv data, small-scale may have subtle shifts from phase noise
for i in data_params.keys():
       print('Imaging MS: ',i) 
       if 'LB' in i:
          image_each_obs(data_params[i], prefix, scales=LB_scales,  uvtaper=outertaper,
                   nsigma=5.0, sidelobethreshold=2.5, smoothfactor=1.5,interactive=False,parallel=parallel) 
       else:
          image_each_obs(data_params[i], prefix, scales=SB_scales, 
                   nsigma=5.0, sidelobethreshold=2.5, interactive=False,parallel=parallel)

       #check masks to ensure you are actually masking the image, lower sidelobethreshold if needed
	
""" Fit Gaussians to roughly estimate centers, inclinations, PAs """
""" Loops through each dataset specified """
###default fit region is blank for an obvious single source
fit_region=''

# ###specify manual mask on brightest source if Gaussian fitting fails due to confusion

# mask_ra  =  '19h01m56.419s'.replace('h',':').replace('m',':').replace('s','')
# mask_dec = '-36d57m28.690s'.replace('d','.').replace('m','.').replace('s','')
# mask_pa  = 90.0 	# position angle of mask in degrees
# mask_maj = 0.76	# semimajor axis of mask in arcsec
# mask_min = 0.75 	# semiminor axis of mask in arcsec
# fit_region = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % \
#               (mask_ra, mask_dec, mask_maj, mask_min, mask_pa)

for i in data_params.keys():
       print(i)
       data_params[i]['phasecenter']=fit_gaussian(prefix+'_'+i+'_initcont_exec0.image', region=fit_region,mask=prefix+'_'+i+'_initcont_exec0.mask')

### output on the terminal ###
# SB1
# 04h04m43.070030s +26d18m56.20026s
# #Peak of Gaussian component identified with imfit: J2000 04h04m43.070030s +26d18m56.20026s
# #PA of Gaussian component: 49.10 deg
# #Inclination of Gaussian component: 58.79 deg
# #Pixel coordinates of peak: x = 351.364 y = 446.699
# LB1
# 04h04m43.079755s +26d18m56.11819s
# #Peak of Gaussian component identified with imfit: ICRS 04h04m43.079755s +26d18m56.11819s
# 04h04m43.079755s +26d18m56.11819s
# Separation: radian = 7.84366e-08, degrees = 0.000004 = 4.49409e-06, arcsec = 0.016179 = 0.0161787
# #Peak in J2000 coordinates: 04:04:43.08036, +026:18:56.104205
# #PA of Gaussian component: 29.88 deg
# #Inclination of Gaussian component: 74.34 deg
# #Pixel coordinates of peak: x = 470.055 y = 1439.628
# LB2
# 04h04m43.079780s +26d18m56.11748s
# #Peak of Gaussian component identified with imfit: ICRS 04h04m43.079780s +26d18m56.11748s
# 04h04m43.079780s +26d18m56.11748s
# Separation: radian = 7.82732e-08, degrees = 0.000004 = 4.48473e-06, arcsec = 0.016145 = 0.016145
# #Peak in J2000 coordinates: 04:04:43.08038, +026:18:56.103495
# #PA of Gaussian component: 6.99 deg
# #Inclination of Gaussian component: 77.90 deg
# #Pixel coordinates of peak: x = 469.941 y = 1439.392

### Check phase center fits in viewer, if centers appear too shifted from the Gaussian fit, 
### manually set the phase center dictionary entry by eye

""" The emission centers are slightly misaligned.  So we split out the 
    individual executions, shift the peaks to the phase center, and reassign 
    the phase centers to a common direction. """

### Set common direction for each EB using one as reference (typically best looking LB image)
# choose SB1 for reference since it looks best in imview
for i in data_params.keys():
       #################### MANUALLY SET THIS ######################
       data_params[i]['common_dir']='J2000 04h04m43.070030s +26d18m56.20026s'

### save updated data params to a pickle

with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)
	

###############################################################
#################### SHIFT PHASE CENTERS ######################
###############################################################

for i in data_params.keys():
   print(i)
   data_params[i]['vis_avg_shift']=prefix+'_'+i+'_initcont_shift.ms'
   os.system('rm -rf '+data_params[i]['vis_avg_shift'])
   fixvis(vis=data_params[i]['vis_avg'], outputvis=data_params[i]['vis_avg_shift'], 
       field=data_params[i]['field'], 
       phasecenter='J2000 '+data_params[i]['phasecenter'])
   ### fix planets may throw an error, usually safe to ignore
   fixplanets(vis=data_params[i]['vis_avg_shift'], field=data_params[i]['field'], 
           direction=data_params[i]['common_dir'])
	
# SB1
# LB1
# 2021-11-26 03:56:29	WARN	fixplanets::::casa	The three FIELD table direction reference frame entries for field 0 are not identical in the input data: 0, 21, 21. Will try to continue ...
# LB2
# 2021-11-26 03:57:10	WARN	fixplanets::::casa	The three FIELD table direction reference frame entries for field 0 are not identical in the input data: 0, 21, 21. Will try to continue ...

###############################################################
############### REIMAGING TO CHECK ALIGNMENT ##################
###############################################################
for i in data_params.keys():
       print(i)
       if 'SB' in i:
          scales=SB_scales
       else:
          scales=LB_scales
       for suffix in ['image','mask','mode','psf','pb','residual','sumwt']:
          os.system('rm -rf '+prefix+'_'+i+'_initcont_shift.'+suffix)
       image_each_obs_shift(data_params[i]['vis_avg_shift'], prefix, scales=scales, 
                   nsigma=5.0, sidelobethreshold=2.5, interactive=False,parallel=parallel)

for i in data_params.keys():
      print(i)     
      data_params[i]['phasecenter_new']=fit_gaussian(prefix+'_'+i+'_initcont_shift.image',\
                                                     region=fit_region,mask=prefix+'_'+i+'_initcont_shift.mask')
      print('Phasecenter new: ',data_params[i]['phasecenter_new'])
      print('Phasecenter old: ',data_params[i]['phasecenter'])

### save updated data params to a pickle
with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)


### output on the terminal ###
# SB1
# 04h04m43.070001s +26d18m56.20011s
# #Peak of Gaussian component identified with imfit: J2000 04h04m43.070001s +26d18m56.20011s
# #PA of Gaussian component: 49.05 deg
# #Inclination of Gaussian component: 58.87 deg
# #Pixel coordinates of peak: x = 449.984 y = 449.994
# Phasecenter new:  04h04m43.070001s +26d18m56.20011s
# Phasecenter old:  04h04m43.070030s +26d18m56.20026s
# LB1
# 04h04m43.070171s +26d18m56.19898s
# #Peak of Gaussian component identified with imfit: J2000 04h04m43.070171s +26d18m56.19898s
# #PA of Gaussian component: 28.04 deg
# #Inclination of Gaussian component: 62.99 deg
# #Pixel coordinates of peak: x = 1499.081 y = 1499.562
# Phasecenter new:  04h04m43.070171s +26d18m56.19898s
# Phasecenter old:  04h04m43.08036s +026d18m56.104205s
# LB2
# 04h04m43.070143s +26d18m56.20069s
# #Peak of Gaussian component identified with imfit: J2000 04h04m43.070143s +26d18m56.20069s
# #PA of Gaussian component: 176.29 deg
# #Inclination of Gaussian component: 61.18 deg
# #Pixel coordinates of peak: x = 1499.206 y = 1500.132
# Phasecenter new:  04h04m43.070143s +26d18m56.20069s
# Phasecenter old:  04h04m43.08038s +026d18m56.103495s

###############################################################
############### PLOT UV DATA TO CHECK SCALING #################
###############################################################

### Assign rough emission geometry parameters; keep 0, 0
PA, incl = 0, 0

### Export MS contents into Numpy save files 
export_vislist=[]
for i in data_params.keys():
   export_MS(data_params[i]['vis_avg_shift'])
   export_vislist.append(data_params[i]['vis_avg_shift'].replace('.ms','.vis.npz'))

if not skip_plots:
    ### Plot deprojected visibility profiles for all data together """
    plot_deprojected(export_vislist,
                     fluxscale=[1.0]*len(export_vislist), PA=PA, incl=incl, 
                     show_err=False)

### Now inspect offsets by comparing against a reference 
### Set reference data using the dictionary key.
### Use LB2 as a reference since it has best noise level and looks

#################### MANUALLY SET THIS ######################

refdata='LB2'

reference=prefix+'_'+refdata+'_initcont_shift.vis.npz'
for i in data_params.keys():
   print(i)
   if i != refdata:
      data_params[i]['gencal_scale']=estimate_flux_scale(reference=reference, 
                        comparison=prefix+'_'+i+'_initcont_shift.vis.npz', 
                        incl=incl, PA=PA)
   else:
      data_params[i]['gencal_scale']=1.0
   print(' ')

### output on the terminal ###
# SB1
# #The ratio of the fluxes of L1489IRS_SB1_initcont_shift.vis.npz to L1489IRS_LB2_initcont_shift.vis.npz is 0.99079
# #The scaling factor for gencal is 0.995 for your comparison measurement
# #The error on the weighted mean ratio is 9.843e-03, although it's likely that the weights in the measurement sets are off by some constant factor
 
# LB1
# #The ratio of the fluxes of L1489IRS_LB1_initcont_shift.vis.npz to L1489IRS_LB2_initcont_shift.vis.npz is 0.80145
# #The scaling factor for gencal is 0.895 for your comparison measurement
# #The error on the weighted mean ratio is 4.546e-03, although it's likely that the weights in the measurement sets are off by some constant factor
 
# LB2


###############################################################
############### SCALE DATA RELATIVE TO ONE EB #################
###############################################################

os.system('rm -rf *_rescaled.ms')
for i in data_params.keys():
   rescale_flux(data_params[i]['vis_avg_shift'], [data_params[i]['gencal_scale']])
   rescale_flux(data_params[i]['vis_avg'], [data_params[i]['gencal_scale']])
   data_params[i]['vis_avg_shift_rescaled']=data_params[i]['vis_avg_shift'].replace('.ms','_rescaled.ms')
   data_params[i]['vis_avg_rescaled']=data_params[i]['vis_avg'].replace('.ms','_rescaled.ms')

###############################################################
############## PLOT UV DATA TO CHECK RE-SCALING ###############
###############################################################

if not skip_plots:
    ### Assign rough emission geometry parameters; keep 0, 0
   PA, incl = 0, 0

   ### Check that rescaling did what we expect
   export_vislist_rescaled=[]
   for i in data_params.keys():
      export_MS(data_params[i]['vis_avg_shift_rescaled'])
      export_vislist_rescaled.append(data_params[i]['vis_avg_shift_rescaled'].replace('.ms','.vis.npz'))

   plot_deprojected(export_vislist_rescaled,
                     fluxscale=[1.0]*len(export_vislist_rescaled), PA=PA, incl=incl, 
                     show_err=False)
   ### Make sure differences are no longer significant
   refdata='LB2'
   reference=prefix+'_'+refdata+'_initcont_shift.vis.npz'
   for i in data_params.keys():
      if i != refdata:
         estimate_flux_scale(reference=reference, 
                        comparison=prefix+'_'+i+'_initcont_shift_rescaled.vis.npz', 
                        incl=incl, PA=PA)
		
#The ratio of the fluxes of L1489IRS_SB1_initcont_shift_rescaled.vis.npz to L1489IRS_LB2_initcont_shift.vis.npz is 1.00000
#The scaling factor for gencal is 1.000 for your comparison measurement
#The error on the weighted mean ratio is 9.934e-03, although it's likely that the weights in the measurement sets are off by some constant factor
#The ratio of the fluxes of L1489IRS_LB1_initcont_shift_rescaled.vis.npz to L1489IRS_LB2_initcont_shift.vis.npz is 1.00000
#The scaling factor for gencal is 1.000 for your comparison measurement
#The error on the weighted mean ratio is 5.672e-03, although it's likely that the weights in the measurement sets are off by some constant factor

### Save data params to a pickle
with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)
	

###############################################################
################ SELF-CALIBRATION PREPARATION #################
###############################################################
#selectedVis='vis_avg_rescaled'
selectedVis='vis_avg_shift_rescaled'

### determine best reference antennas based on geometry and flagging
for i in data_params.keys():
   data_params[i]["refant"] = rank_refants(data_params[i][selectedVis])

############### CHECK THESE, SHOULD BE FINE #################
SB_spwmap=[0,0,0,0,0,0] # mapping of spws from the one which has calibration solution to the one which is to be applied the solution
SB_contspws = '' # mean all spws 


### Make a list of EBs to image
vislist=[]
for i in data_params.keys():
      if ('LB' in i): # skip over LB EBs if in SB-only mode
         continue
      vislist.append(data_params[i][selectedVis])


""" Set up a clean mask """

mask_ra  =  data_params[i]['common_dir'].split()[1].replace('h',':').replace('m',':').replace('s','')
mask_dec = data_params[i]['common_dir'].split()[2].replace('d','.').replace('m','.').replace('s','')
mask_pa  = 50.0 	# position angle of mask in degrees
mask_maj = 3.0	# semimajor axis of mask in arcsec
mask_min = 2.0 	# semiminor axis of mask in arcsec

common_mask = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % \
              (mask_ra, mask_dec, mask_maj, mask_min, mask_pa)
""" Define a noise annulus, measure the peak SNR in map """
noise_annulus = "annulus[[%s, %s],['%.2farcsec', '8.0arcsec']]" % \
                (mask_ra, mask_dec, 2.0*mask_maj) 


###############################################################
###################### SELF-CALIBRATION #######################
###############################################################

### Initial dirty map to assess DR
tclean_wrapper(vis=vislist, imagename=prefix+'_dirty', 
               scales=SB_scales, niter=0,parallel=parallel,cellsize='0.025arcsec',imsize=1600, nterms=1)
estimate_SNR(prefix+'_dirty.image.tt0', disk_mask=common_mask, 
             noise_mask=noise_annulus)

#nterms=1
#L1489IRS_dirty.image.tt0
#Beam 0.236 arcsec x 0.157 arcsec (33.22 deg)
#Flux inside disk mask: 16.71 mJy
#Peak intensity of source: 4.54 mJy/beam
#rms: 9.12e-02 mJy/beam
#Peak SNR: 49.79

# nterms=2 (suboptimal SNR, so use nterms = 1)
#L1489IRS_dirty.image.tt0
#Beam 0.236 arcsec x 0.157 arcsec (33.22 deg)
#Flux inside disk mask: 14.19 mJy
#Peak intensity of source: 4.49 mJy/beam
#rms: 1.35e-01 mJy/beam
#Peak SNR: 33.22

### Image produced by iter 0 has not selfcal applied, it's used to set the initial model
### only images >0 have self-calibration applied

### Run self-calibration command set
### 0. Split off corrected data from previous selfcal iteration (except iteration 0)
### 1. Image data to specified nsigma depth, set model column
### 2. Calculate self-cal gain solutions
### 3. Apply self-cal gain solutions to MS

############# USERS MAY NEED TO ADJUST NSIGMA AND SOLINT FOR EACH SELF-CALIBRATION ITERATION ##############
iteration=0
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='p',nsigma=15.0,solint='inf',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel, nterms=1)


### Plot gain corrections, loop through each
if not skip_plots:
   for i in data_params.keys():
     if 'SB' in i:
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True, plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")

		
#L1489IRS_SB-only_p0.image.tt0
#Beam 0.236 arcsec x 0.157 arcsec (33.22 deg)
#Flux inside disk mask: 18.81 mJy
#Peak intensity of source: 4.54 mJy/beam
#rms: 8.36e-02 mJy/beam
#Peak SNR: 54.29

#L1489IRS_SB-only_p0_post.image.tt0
#Beam 0.236 arcsec x 0.157 arcsec (33.22 deg)
#Flux inside disk mask: 20.70 mJy
#Peak intensity of source: 5.59 mJy/beam
#rms: 8.00e-02 mJy/beam
#Peak SNR: 69.83

iteration=1
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='p',nsigma=3.0,solint='30s',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel, nterms=1, finalimageonly=True)

if not skip_plots:
   for i in data_params.keys():
     if 'SB' in i:
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True, plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")
		
#L1489IRS_SB-only_p1.image.tt0
#Beam 0.236 arcsec x 0.157 arcsec (33.22 deg)
#Flux inside disk mask: 26.98 mJy
#Peak intensity of source: 5.64 mJy/beam
#rms: 7.39e-02 mJy/beam
#Peak SNR: 76.31

"""
#L1489IRS_SB-only_p1_post.image.tt0
#Beam 0.236 arcsec x 0.157 arcsec (33.22 deg)
#Flux inside disk mask: 27.66 mJy
#Peak intensity of source: 5.63 mJy/beam
#rms: 7.44e-02 mJy/beam
#Peak SNR: 75.72
"""

### self-cal stopped here since the SNR go down.

"""

### SNR slightly down: see also the threshold SNR to work selfcal well > https://casaguides.nrao.edu/index.php/Self-Calibration_Template
### here, N=35 (number of anntenas), t_int ~ 30 min, t_solint ~ 2 min (scan length)
### -> needed SNR for selfcal ~ 70; but actual SNR is 50...

### Changing self-cal mode here to ap, see use of prevselfcalmode to ensure proper split to see if any improvement
iteration=2
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='ap',prevselfcalmode='p',nsigma=3.0,solint='inf',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel,nterms=1)

if not skip_plots:
   for i in data_params.keys():
     if 'SB' in i:
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_ap'+str(iteration)+'.g'), xaxis='time',
              yaxis='amp',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True,plotrange=[0,0,0,2])
       input("Press Enter key to advance to next MS/Caltable...")
		
#L1489IRS_SB-only_p2.image.tt0
#Beam 0.236 arcsec x 0.157 arcsec (33.22 deg)
#Flux inside disk mask: 27.24 mJy
#Peak intensity of source: 5.58 mJy/beam
#rms: 7.37e-02 mJy/beam
#Peak SNR: 75.66

#L1489IRS_SB-only_p2_post.image.tt0
#Beam 0.236 arcsec x 0.159 arcsec (32.57 deg)
#Flux inside disk mask: 24.43 mJy
#Peak intensity of source: 6.11 mJy/beam
#rms: 7.71e-02 mJy/beam
#Peak SNR: 79.29

iteration=3
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='ap',nsigma=3.0,solint='18s',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel,finalimageonly=True,nterms=1)

if not skip_plots:
   for i in data_params.keys():
     if 'SB' in i:
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_ap'+str(iteration)+'.g'), xaxis='time',
              yaxis='amp',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True,plotrange=[0,0,0,2])
       input("Press Enter key to advance to next MS/Caltable...")
		
#L1489IRS_SB-only_ap3.image.tt0
#Beam 0.236 arcsec x 0.159 arcsec (32.57 deg)
#Flux inside disk mask: 25.08 mJy
#Peak intensity of source: 6.11 mJy/beam
#rms: 7.66e-02 mJy/beam
#Peak SNR: 79.72

"""

for i in data_params.keys():
   if 'SB' in i:
      data_params[i]['selfcal_spwmap_SB-only']=data_params[i]['selfcal_spwmap'].copy()
      data_params[i]['selfcal_tables_SB-only']=data_params[i]['selfcal_tables'].copy()
      data_params[i]['vis_avg_selfcal_SB-only']=(data_params[i]['vis_avg_selfcal']+'.')[:-1]  ## trick to copy the string

### Save data params to a pickle
with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)

	
###############################################################
################### SELF-CALIBRATION SB+LB ####################
###############################################################
LB_spwmap=[0,0,0,0,0,0,0]
LB_contspws = '' 


### Make a list of EBs to image
vislist=[]
for i in data_params.keys():
   vislist.append(data_params[i][selectedVis])

### Initial dirty map to assess DR
tclean_wrapper(vis=vislist, imagename=prefix+'_LB+SB_dirty', 
               scales=SB_scales, niter=0,parallel=parallel,cellsize='0.003arcsec',imsize=7200, nterms=1)
estimate_SNR(prefix+'_LB+SB_dirty.image.tt0', disk_mask=common_mask, 
             noise_mask=noise_annulus)

# nterms=1
#L1489IRS_LB+SB_dirty.image.tt0
#Beam 0.060 arcsec x 0.036 arcsec (20.78 deg)
#Flux inside disk mask: 27.40 mJy
#Peak intensity of source: 3.26 mJy/beam
#rms: 1.82e-02 mJy/beam
#Peak SNR: 178.67

# nterms=2 (suboptimal use nterms=1)
#L1489IRS_LB+SB_dirty.image.tt0
#Beam 0.061 arcsec x 0.036 arcsec (22.27 deg)
#Flux inside disk mask: 39.26 mJy
#Peak intensity of source: 4.13 mJy/beam
#rms: 2.38e-02 mJy/beam
#Peak SNR: 173.39

### Image produced by iter 0 has not selfcal applied, it's used to set the initial model
### only images >0 have self-calibration applied

### Run self-calibration command set
### 0. Split off corrected data from previous selfcal iteration (except iteration 0)
### 1. Image data to specified nsigma depth, set model column
### 2. Calculate self-cal gain solutions
### 3. Apply self-cal gain solutions to MS

############# USERS MAY NEED TO ADJUST NSIGMA AND SOLINT FOR EACH SELF-CALIBRATION ITERATION ##############
############################ CONTINUE SELF-CALIBRATION ITERATIONS UNTIL ###################################
#################### THE S/N BEGINS TO DROP OR SOLINTS ARE AS LOW AS POSSIBLE #############################

#################### LOOK FOR ERRORS IN GAINCAL CLAIMING A FREQUENCY MISMATCH #############################
####### IF FOUND, CHANGE SOLINT, MAYBE TRY TO ALIGN WITH A CERTAIN NUMBER OF SCANS AND TRY AGAIN ##########
########################## IF ALL ELSE FAILS, SIMPLY START WITH solint='inf' ##############################

iteration=0
self_calibrate(prefix,data_params,selectedVis,mode='LB+SB',iteration=iteration,selfcalmode='p',nsigma=20.0,solint='inf',
               noisemasks=[common_mask,noise_annulus],
SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,LB_contspws=LB_contspws,LB_spwmap=LB_spwmap,combine='spw,scan',parallel=parallel,smoothfactor=2.0,imsize=7200,nterms=1)#, noisethreshold=3.0, lownoisethreshold=1.0)


### Plot gain corrections, loop through each
if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_LB+SB_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True, plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")

#L1489IRS_LB+SB_p0.image.tt0
#Beam 0.060 arcsec x 0.036 arcsec (20.78 deg)
#Flux inside disk mask: 30.42 mJy
#Peak intensity of source: 3.27 mJy/beam
#rms: 1.75e-02 mJy/beam
#Peak SNR: 186.49

#L1489IRS_LB+SB_p0_post.image.tt0
#Beam 0.060 arcsec x 0.036 arcsec (20.78 deg)
#Flux inside disk mask: 31.11 mJy
#Peak intensity of source: 3.74 mJy/beam
#rms: 1.73e-02 mJy/beam
#Peak SNR: 216.69



iteration=1
self_calibrate(prefix,data_params,selectedVis,mode='LB+SB',iteration=iteration,selfcalmode='p',nsigma=5.0,solint='inf',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,LB_contspws=LB_contspws,LB_spwmap=LB_spwmap,combine='spw,scan',parallel=parallel,smoothfactor=2.0,imsize=7200,nterms=1)#,noisethreshold=3.0, lownoisethreshold=1.0)

if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_LB+SB_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True, plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")

#L1489IRS_LB+SB_p1.image.tt0
#Beam 0.061 arcsec x 0.036 arcsec (22.27 deg)
#Flux inside disk mask: 34.02 mJy
#Peak intensity of source: 4.73 mJy/beam
#rms: 2.17e-02 mJy/beam
#Peak SNR: 217.67

#L1489IRS_LB+SB_p1_post.image.tt0
#Beam 0.061 arcsec x 0.036 arcsec (22.27 deg)
#Flux inside disk mask: 33.45 mJy
#Peak intensity of source: 5.50 mJy/beam
#rms: 2.14e-02 mJy/beam
#Peak SNR: 256.51


############################# finished up to here in 2021/11/10 ###################################

iteration=2
self_calibrate(prefix,data_params,selectedVis,mode='LB+SB',iteration=iteration,selfcalmode='p',nsigma=5.0,solint='12s',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,LB_contspws=LB_contspws,LB_spwmap=LB_spwmap,combine='spw,scan',parallel=parallel,smoothfactor=2.0,imsize=7200,nterms=1,noisethreshold=3.0, lownoisethreshold=1.0,finalimageonly=True)

#L1489IRS_LB+SB_p2.image.tt0
#Beam 0.061 arcsec x 0.036 arcsec (22.27 deg)
#Flux inside disk mask: 34.01 mJy
#Peak intensity of source: 5.48 mJy/beam
#rms: 2.13e-02 mJy/beam
#Peak SNR: 256.96

"""


iteration=2
self_calibrate(prefix,data_params,selectedVis,mode='LB+SB',iteration=iteration,selfcalmode='p',nsigma=5.0,solint='18s',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,LB_contspws=LB_contspws,LB_spwmap=LB_spwmap,combine='spw,scan',parallel=parallel,smoothfactor=2.0,imsize=7200,nterms=1,noisethreshold=3.0, lownoisethreshold=1.0)

if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_LB+SB_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True,plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")

#L1489IRS_LB+SB_p2.image.tt0
#Beam 0.061 arcsec x 0.036 arcsec (22.27 deg)
#Flux inside disk mask: 34.00 mJy
#Peak intensity of source: 5.49 mJy/beam
#rms: 2.13e-02 mJy/beam
#Peak SNR: 257.45

#L1489IRS_LB+SB_p2_post.image.tt0
#Beam 0.061 arcsec x 0.036 arcsec (22.27 deg)
#Flux inside disk mask: 23.73 mJy
#Peak intensity of source: 3.08 mJy/beam
#rms: 2.16e-02 mJy/beam
#Peak SNR: 142.62 --- something going wrong...


### change to amplitude selfcal
iteration=3
self_calibrate(prefix,data_params,selectedVis,mode='LB+SB',iteration=iteration,selfcalmode='ap',prevselfcalmode='p',nsigma=3.0,solint='inf',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,LB_contspws=LB_contspws,LB_spwmap=LB_spwmap,parallel=parallel,combine='spw',smoothfactor=2.0,imsize=7200, nterms=1,noisethreshold=3.0, lownoisethreshold=1.0)

if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_LB+SB_ap'+str(iteration)+'.g'), xaxis='time',
              yaxis='amp',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True,plotrange=[0,0,0,2])
       input("Press Enter key to advance to next MS/Caltable...")

#L1489IRS_LB+SB_p3.image.tt0
#Beam 0.061 arcsec x 0.036 arcsec (22.27 deg)
#Flux inside disk mask: 22.91 mJy
#Peak intensity of source: 3.06 mJy/beam
#rms: 2.13e-02 mJy/beam
#Peak SNR: 143.46
	

#L1489IRS_LB+SB_p3_post.image.tt0
#Beam 0.061 arcsec x 0.036 arcsec (22.63 deg)
#Flux inside disk mask: 22.54 mJy
#Peak intensity of source: 2.46 mJy/beam
#rms: 2.16e-02 mJy/beam
#Peak SNR: 114.14

# getting worse, stop
"""

###Backup gain table list for LB+SB runs
for i in data_params.keys():
   if 'SB' in i:
      data_params[i]['selfcal_spwmap']=data_params[i]['selfcal_spwmap_SB-only']+data_params[i]['selfcal_spwmap']
      data_params[i]['selfcal_tables']=data_params[i]['selfcal_tables_SB-only']+data_params[i]['selfcal_tables']

#save data params to a pickle
with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)



###############################################################
################# SPLIT OFF FINAL CONT DATA ###################
###############################################################

for i in data_params.keys():
   os.system('rm -rf '+prefix+'_'+i+'_continuum.ms '+prefix+'_'+i+'_continuum.ms.tgz')
   split(vis=data_params[i]['vis_avg_selfcal'], outputvis=prefix+'_'+i+'_continuum.ms',
      datacolumn='data')
   data_params[i]['vis_final']=prefix+'_'+i+'_continuum.ms'
   os.system('tar cvzf '+prefix+'_'+i+'_continuum.ms.tgz '+prefix+'_'+i+'_continuum.ms')

#save data params to a pickle
with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)


###############################################################
################## RUN A FINAL IMAGE SET ######################
###############################################################

### Generate a vislist
vislist=[]
for i in data_params.keys():
   vislist.append(data_params[i]['vis_final'])
scales = SB_scales
imsize=7200
cell='0.003arcsec'
# for robust in [-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0]:
for robust in [0.5]:
    imagename=prefix+'_SBLB_continuum_robust_'+str(robust)
    os.system('rm -rf '+imagename+'*')

    sigma = get_sensitivity(data_params, specmode='mfs',imsize=imsize,robust=robust,cellsize=cell)
	print(sigma)
    # adjust sigma to corrector for the irregular noise in the images if needed
    # correction factor may vary or may not be needed at all depending on source
#     if robust == 2.0 or robust == 1.0:
#       sigma=sigma*1.75
    tclean_wrapper(vis=vislist, imagename=imagename, sidelobethreshold=2.0,
            smoothfactor=2.0, scales=scales, threshold=3.0*sigma, 
            noisethreshold=3.0, lownoisethreshold=1.0, robust=robust, parallel=parallel, 
            cellsize=cell, imsize=imsize, nterms=1, uvtaper='2000klambda')#,phasecenter=data_params['SB1']['common_dir'].replace('J2000','ICRS'))

    imagename=imagename+'.image.tt0'
    exportfits(imagename=imagename, fitsimage=imagename+'.fits',overwrite=True,dropdeg=True)

###############################################################
########################### CLEANUP ###########################
###############################################################

### Remove extra image products
os.system('rm -rf *.residual* *.psf* *.model* *dirty* *.sumwt* *.gridwt* *.workdirectory')

### put selfcalibration intermediate images somewhere safe
os.system('rm -rf initial_images')
os.system('mkdir initial_images')
os.system('mv *initcont*.image *_p*.image* *_ap*.image* initial_images')
os.system('mv *initcont*.mask *_p*.mask *_ap*.mask initial_images')
os.system('rm -rf *_p*.alpha* *_p*.pb.tt0 *_ap*.alpha* *_ap*.pb.tt0')

### Remove intermediate selfcal MSfiles
os.system("rm -rf *p{0..99}.ms")
os.system("rm -rf *p{0..99}.ms.flagversions")
### Remove rescaled selfcal MSfiles
os.system('rm -rf *rescaled.ms')
os.system('rm -rf *rescaled.ms.flagversions')
### Remove rescaled selfcal MSfiles
os.system('rm -rf *initcont*.ms')
os.system('rm -rf *initcont*.ms.flagversions')

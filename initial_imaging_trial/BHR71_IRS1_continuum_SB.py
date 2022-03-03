"""
eDisk data reduction script
This script was written for CASA 6.1.1/6.2
Originally derived from DSHARP reduction scripts


Datasets calibrated (in order of date observed):
SB1: 2019.1.00261.L (04-May-2021; uid___A002_Xebb7f0_X690.ms)
SB2: 2019.1.00261.L (09-May-2021; uid___A002_Xebd1e8_X38be.ms)
 
LB1: 
     

reducer: Yoshihide Yamato
"""

#########################################################################
### have to go to the directry of WD_path to execute the script below ###
#########################################################################

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
skip_plots = False	

### start this script from the top or restart from a certain point?
restart = False

### Add field names (corresponding to the field in the MS) here and prefix for 
### filenameing (can be different but try to keep same)
### Only make different if, for example, the field name has a space
field   = 'BHR71_IRS1'
prefix  = 'BHR71_IRS1' 

### always include trailing slashes!!
WD_path = '/raid/work/yamato/eDisk_data/BHR71_IRS1/'
SB_path = WD_path+'SB/'
LB_path = WD_path+'LB/'

### scales for multi-scale clean
SB_scales = [0, 5] #[0, 5, 10, 20]
LB_scales = [0, 5, 30]  #[0, 5, 30, 100, 200]

### Add additional dictionary entries if need, i.e., SB2, SB3, LB1, LB2, etc. for each execution
### Note that C18O and 13CO have different spws in the DDT vis LP os the spw ordering
### is different for data that were originally part of the DDT than the LP
### DDT 2019.A.00034.S SB data need 'spws': '25,31,29,27,33,35,37'
### LP  2019.1.00261.L SB data need 'spws': '25,27,29,31,33,35,37'
pl_data_params={'SB1': {'vis': SB_path+'uid___A002_Xebb7f0_X690.ms',
                        'spws': '25,27,29,31,33,35,37'},
				'SB2': {'vis': SB_path+'uid___A002_Xebd1e8_X38be.ms',
                        'spws': '25,27,29,31,33,35,37'},
               }

### Dictionary defining necessary metadata for each execution
### SiO at 217.10498e9 excluded because of non-detection
### Only bother specifying simple species that are likely present in all datasets
### Hot corino lines (or others) will get taken care of by using the cont.dat

if not restart:
	### may able to make a script to automatically produce these dictionaries?
	data_params = {'SB1': {'vis' : WD_path+prefix+'_SB1.ms',
						   'name' : 'SB1',
						   'field': field,
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
						   'timerange': '2021/05/04/00:27:00~2021/05/04/01:35:00',
						   'contdotdat' : 'SB/cont.dat'
						  },
				   'SB2': {'vis' : WD_path+prefix+'_SB2.ms',
						   'name' : 'SB2',
						   'field': field,
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
						   'timerange': '2021/05/09/01:28:00~2021/05/09/02:35:00',
						   'contdotdat' : 'SB/cont.dat'
						  },
				   }

	'''    
	No LB yet                 
				   'LB1': {'vis' : LB1_path,
						   'name' : 'LB1',
						   'field' : 'Ced110IRS4',
						   'line_spws': np.array([]), # CO SPWs 
						   'line_freqs': np.array([]),
						   'flagrange': np.array([]), 
						   'cont_spws':  np.array([]), # CO SPWs 
						   'cont_avg_width':  np.array([]), 
						   'phasecenter': ' ',
						   'timerange': '2015/11/01/00:00:00~2015/11/02/00:00:00',
						  }

				   }

	'''
	
elif restart:
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
### do the below at WD_path directory
for i in pl_data_params.keys():
   if os.path.exists(prefix+'_'+i+'.ms'):
      flagmanager(vis=prefix+'_'+i+'.ms', mode="restore", \
                  versionname="starting_flags")
   else:
      split(vis=pl_data_params[i]['vis'],outputvis=prefix+'_'+i+'.ms',spw=pl_data_params[i]['spws'],field=field,datacolumn='corrected')

### BHR71_IRS1_SBn.ms (n=1,2) are produced

### Backup the the flagging state at start of reduction
for i in data_params.keys():
    if not os.path.exists(data_params[i]['vis']+\
            ".flagversions/flags.starting_flags"):
       flagmanager(vis=data_params[i]['vis'], mode = 'save', versionname = 'starting_flags', comment = 'Flag states at start of reduction')

### BHR71_IRS1_SBn.ms (n=1,2) are produced


### Inspect data in each spw for each dataset

### may be able to write a couple of script to plot the visibilities more fastly? (plotms is slow)
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

### BHR71_IRS1_SBn_initcont.ms (n=1,2) are produced
      
###############################################################
############## INITIAL IMAGING FOR ALIGNMENT ##################
###############################################################


### Image each dataset individually to get source position in each image
### Images are saved in the format prefix+'_name_initcont_exec#.ms'
outertaper='2000klambda' # taper if necessary to align using larger-scale uv data, small-scale may have subtle shifts from phase noise
for i in data_params.keys():
       print('Imaging SB: ',i) 
       if 'LB' in i:
          image_each_obs(data_params[i], prefix, scales=LB_scales,  uvtaper=outertaper,
                   nsigma=5.0, sidelobethreshold=2.5, smoothfactor=1.5,interactive=False,parallel=parallel) 
       else:
          image_each_obs(data_params[i], prefix, scales=SB_scales, 
                   nsigma=5.0, sidelobethreshold=2.5, interactive=False,parallel=parallel)

### BHR71_IRS1_SBn_initcont_exec0 (n=1,2) with the extenstions of [gridwt, image, mask, mode, pb, psf, residual, sumwt] are produced

       #check masks to ensure you are actually masking the image by imview, lower sidelobethreshold if needed 

""" Fit Gaussians to roughly estimate centers, inclinations, PAs """
""" Loops through each dataset specified """
###default fit region is blank for an obvious single source -> OK for BHR71 IRS1
fit_region=''
for i in data_params.keys():
       print(i)
       data_params[i]['phasecenter']=fit_gaussian(prefix+'_'+i+'_initcont_exec0.image', region=fit_region,mask=prefix+'_'+i+'_initcont_exec0.mask')

### output on the terminal
# SB1
# 12h01m36.474422s -65d08m49.35978s
# #Peak of Gaussian component identified with imfit: ICRS 12h01m36.474422s -65d08m49.35978s
# 12h01m36.474422s -65d08m49.35978s
# Separation: radian = 1.00182e-07, degrees = 0.000006 = 5.73998e-06, arcsec = 0.020664 = 0.0206639
# #Peak in J2000 coordinates: 12:01:36.47636, -065:08:49.343115
# #PA of Gaussian component: 104.78 deg
# #Inclination of Gaussian component: 33.98 deg
# #Pixel coordinates of peak: x = 520.533 y = 445.308
# SB2
# 12h01m36.478689s -65d08m49.35456s
# #Peak of Gaussian component identified with imfit: ICRS 12h01m36.478689s -65d08m49.35456s
# 12h01m36.478689s -65d08m49.35456s
# Separation: radian = 1.00232e-07, degrees = 0.000006 = 5.74287e-06, arcsec = 0.020674 = 0.0206743
# #Peak in J2000 coordinates: 12:01:36.48063, -065:08:49.337896
# #PA of Gaussian component: 93.69 deg
# #Inclination of Gaussian component: 38.71 deg
# #Pixel coordinates of peak: x = 519.637 y = 445.482


### Check phase center fits in viewer, if centers appear too shifted from the Gaussian fit, 
### manually set the phase center dictionary entry by eye

""" The emission centers are slightly misaligned.  So we split out the 
    individual executions, shift the peaks to the phase center, and reassign 
    the phase centers to a common direction. """

### Set common direction for each EB using one as reference (typically best looking LB image)

for i in data_params.keys():
       #################### MANUALLY SET THIS ######################
       data_params[i]['common_dir']='J2000 12h01m36.474422s -65d08m49.35978s'
### manually set to the gaussian peak of SB1

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

### BHR71_IRS1_SBn_initcont_shift (n=1,2) with the extenstions of [gridwt, image, mask, mode, pb, psf, residual, sumwt] are produced
### output on the terminal
# SB1
# 12h01m36.474408s -65d08m49.35990s
# #Peak of Gaussian component identified with imfit: J2000 12h01m36.474408s -65d08m49.35990s
# #PA of Gaussian component: 104.50 deg
# #Inclination of Gaussian component: 34.08 deg
# #Pixel coordinates of peak: x = 449.999 y = 449.979
# Phasecenter new:  12h01m36.474408s -65d08m49.35990s
# Phasecenter old:  12h01m36.47636s -065d08m49.343115s
# SB2
# 12h01m36.474408s -65d08m49.35972s
# #Peak of Gaussian component identified with imfit: J2000 12h01m36.474408s -65d08m49.35972s
# #PA of Gaussian component: 93.70 deg
# #Inclination of Gaussian component: 37.98 deg
# #Pixel coordinates of peak: x = 449.999 y = 449.985
# Phasecenter new:  12h01m36.474408s -65d08m49.35972s
# Phasecenter old:  12h01m36.48063s -065d08m49.337896s
### new phasecenters are aligned well

	
### save updated data params to a pickle
with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)



###############################################################
############### PLOT UV DATA TO CHECK SCALING #################
###############################################################

### Assign rough emission geometry parameters; keep 0, 0
PA, incl = 95, 35
### from the above estimates

### Export MS contents into Numpy save files 
export_vislist=[]
for i in data_params.keys():
   export_MS(data_params[i]['vis_avg_shift'])
   export_vislist.append(data_params[i]['vis_avg_shift'].replace('.ms','.vis.npz'))
#Measurement set exported to BHR71_IRS1_SB1_initcont_shift.vis.npz
#Measurement set exported to BHR71_IRS1_SB2_initcont_shift.vis.npz

if not skip_plots:
    ### Plot deprojected visibility profiles for all data together """
    plot_deprojected(export_vislist,
                     fluxscale=[1.0]*len(export_vislist), PA=PA, incl=incl, 
                     show_err=False)

### Now inspect offsets by comparing against a reference 
### Set reference data using the dictionary key.
### Using SB1 as reference because it looks the nicest by far

#################### MANUALLY SET THIS ######################
refdata='SB1'

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

#The ratio of the fluxes of BHR71_IRS1_SB2_initcont_shift.vis.npz to BHR71_IRS1_SB1_initcont_shift.vis.npz is 0.97572
#The scaling factor for gencal is 0.988 for your comparison measurement
#The error on the weighted mean ratio is 1.669e-04, although it's likely that the weights in the measurement sets are off by some constant factor


###############################################################
############### SCALE DATA RELATIVE TO ONE EB #################
###############################################################

os.system('rm -rf *_rescaled.ms')
for i in data_params.keys():
   rescale_flux(data_params[i]['vis_avg_shift'], [data_params[i]['gencal_scale']])
   rescale_flux(data_params[i]['vis_avg'], [data_params[i]['gencal_scale']])
   data_params[i]['vis_avg_shift_rescaled']=data_params[i]['vis_avg_shift'].replace('.ms','_rescaled.ms')
   data_params[i]['vis_avg_rescaled']=data_params[i]['vis_avg'].replace('.ms','_rescaled.ms')
	
#Splitting out rescaled values into new MS: BHR71_IRS1_SB1_initcont_shift_rescaled.ms
#Splitting out rescaled values into new MS: BHR71_IRS1_SB1_initcont_rescaled.ms
#Splitting out rescaled values into new MS: BHR71_IRS1_SB2_initcont_shift_rescaled.ms
#Splitting out rescaled values into new MS: BHR71_IRS1_SB2_initcont_rescaled.ms


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
   refdata='SB1'
   reference=prefix+'_'+refdata+'_initcont_shift.vis.npz'
   for i in data_params.keys():
      if i != refdata:
         estimate_flux_scale(reference=reference, 
                        comparison=prefix+'_'+i+'_initcont_shift_rescaled.vis.npz', 
                        incl=incl, PA=PA)

### Save data params to a pickle
with open(prefix+'.pickle', 'wb') as handle:
    pickle.dump(data_params, handle, protocol=pickle.HIGHEST_PROTOCOL)



###############################################################
################ SELF-CALIBRATION PREPARATION #################
###############################################################
selectedVis='vis_avg_rescaled'
#selectedVis='vis_avg_shift_rescaled'

### determine best reference antennas based on geometry and flagging
for i in data_params.keys():
   data_params[i]["refant"] = rank_refants(data_params[i][selectedVis])

'''Find reference antenna, pick 2 near array center'''
'''
if not skip_plots:
   for i in data_params.keys():
      if 'LB' in i:
         continue
      listobs(data_params[i]['vis'])
      plotants(data_params[i]['vis'])
      input("Press Enter key to advance to next MS/Caltable...")
'''

'''antenna name is DV/DA/PMXX'''
'''pad number is @AXXX '''
'''want antenna that is on the same pad if possible, list multiple in case one drops out'''
'''check listobs and fill in the SB_refant field '''
'''with the antenna name (DAXX, DVXX, or PMXX) @ pad number (AXXX)'''
'''so make a comma separated list like: DA43@A035,DV07@A011,...'''


#################### MANUALLY SET THIS ######################
#SB_refant   = 'DA43@A035,DV07@A011,DV05@A042' 

############### CHECK THESE, SHOULD BE FINE #################
SB_spwmap=[0,0,0,0,0,0,0]
SB_contspws = '' 


### Make a list of EBs to image
vislist=[]
for i in data_params.keys():
      if ('LB' in i): # skip over LB EBs if in SB-only mode
         continue
      vislist.append(data_params[i][selectedVis])


""" Set up a clean mask """

mask_ra  =  data_params[i]['common_dir'].split()[1].replace('h',':').replace('m',':').replace('s','')
mask_dec = data_params[i]['common_dir'].split()[2].replace('d','.').replace('m','.').replace('s','')
mask_pa  = 100.0 	# position angle of mask in degrees
mask_maj = 1.2	# semimajor axis of mask in arcsec
mask_min = 1.0 	# semiminor axis of mask in arcsec

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
               scales=SB_scales, niter=0,parallel=parallel,cellsize='0.025arcsec',imsize=1600)
estimate_SNR(prefix+'_dirty.image.tt0', disk_mask=common_mask, 
             noise_mask=noise_annulus)

#BHR71_IRS1_dirty.image.tt0
#Beam 0.263 arcsec x 0.200 arcsec (9.67 deg)
#Flux inside disk mask: 737.77 mJy
#Peak intensity of source: 173.68 mJy/beam
#rms: 9.74e-01 mJy/beam
#Peak SNR: 178.40


### Image produced by iter 0 has not selfcal applied, it's used to set the initial model
### only images >0 have self-calibration applied

### Run self-calibration command set
### 0. Split off corrected data from previous selfcal iteration (except iteration 0)
### 1. Image data to specified nsigma depth, set model column
### 2. Calculate self-cal gain solutions
### 3. Apply self-cal gain solutions to MS
### 4. Check S/N before and after

############# USERS MAY NEED TO ADJUST NSIGMA AND SOLINT FOR EACH SELF-CALIBRATION ITERATION ##############
iteration=0
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='p',nsigma=80.0,solint='inf',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel)


### Plot gain corrections, loop through each
if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True, plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")

### Make note of key metrics of image in each round
#Ced110IRS4_SB-only_p0.image.tt0
#Beam 0.447 arcsec x 0.259 arcsec (11.47 deg)
#Flux inside disk mask: 83.51 mJy
#Peak intensity of source: 30.41 mJy/beam
#rms: 1.69e-01 mJy/beam
#Peak SNR: 180.34
iteration=1
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='p',nsigma=25.0,solint='30s',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel)

if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True, plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")

#Ced110IRS4_SB-only_p1.image.tt0
#Beam 0.447 arcsec x 0.259 arcsec (11.47 deg)
#Flux inside disk mask: 84.55 mJy
#Peak intensity of source: 37.08 mJy/beam
#rms: 4.84e-02 mJy/beam
#Peak SNR: 766.55
iteration=2
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='p',nsigma=5.0,solint='6s',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel)

if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True,plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")

#Ced110IRS4_SB-only_p2.image.tt0
#Beam 0.447 arcsec x 0.259 arcsec (11.47 deg)
#Flux inside disk mask: 84.11 mJy
#Peak intensity of source: 41.78 mJy/beam
#rms: 3.65e-02 mJy/beam
#Peak SNR: 1143.40

iteration=3
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='p',nsigma=3.0,solint='int',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel)

if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_p'+str(iteration)+'.g'),
               xaxis='time', yaxis='phase',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True,plotrange=[0,0,-180,180]) 
       input("Press Enter key to advance to next MS/Caltable...")

#Ced110IRS4_SB-only_p3.image.tt0
#Beam 0.447 arcsec x 0.259 arcsec (11.47 deg)
#Flux inside disk mask: 84.84 mJy
#Peak intensity of source: 43.52 mJy/beam
#rms: 3.55e-02 mJy/beam
#Peak SNR: 1227.21


### Changing self-cal mode here to ap, see use of prevselfcalmode to ensure proper split

iteration=4
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='ap',prevselfcalmode='p',nsigma=3.0,solint='inf',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel)

if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_ap'+str(iteration)+'.g'), xaxis='time',
              yaxis='amp',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True,plotrange=[0,0,0,2])
       input("Press Enter key to advance to next MS/Caltable...")

#Ced110IRS4_SB-only_p4.image.tt0
#Beam 0.447 arcsec x 0.259 arcsec (11.47 deg)
#Flux inside disk mask: 84.85 mJy
#Peak intensity of source: 43.53 mJy/beam
#rms: 3.54e-02 mJy/beam
#Peak SNR: 1229.56

iteration=5
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='ap',nsigma=3.0,solint='18s',
               noisemasks=[common_mask,noise_annulus],
               SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,parallel=parallel)

if not skip_plots:
   for i in data_params.keys():
       plotms(vis=data_params[i][selectedVis].replace('.ms','_SB-only_ap'+str(iteration)+'.g'), xaxis='time',
              yaxis='amp',gridrows=4,gridcols=1,iteraxis='antenna', xselfscale=True,plotrange=[0,0,0,2])
       input("Press Enter key tto advance to next MS/Caltable...")

#Ced110IRS4_SB-only_ap5.image.tt0
#Beam 0.449 arcsec x 0.260 arcsec (11.40 deg)
#Flux inside disk mask: 84.54 mJy
#Peak intensity of source: 43.76 mJy/beam
#rms: 3.35e-02 mJy/beam
#Peak SNR: 1304.93

### Make the final image, will not run another self-calibration
iteration=6
self_calibrate(prefix,data_params,selectedVis,mode='SB-only',iteration=iteration,selfcalmode='ap',nsigma=3.0,solint='18s',
               noisemasks=[common_mask,noise_annulus],SB_contspws=SB_contspws,SB_spwmap=SB_spwmap,
               parallel=parallel,finalimageonly=True)

#Ced110IRS4_SB-only_ap6.image.tt0
#Beam 0.450 arcsec x 0.259 arcsec (11.47 deg)
#Flux inside disk mask: 84.50 mJy
#Peak intensity of source: 43.85 mJy/beam
#rms: 3.35e-02 mJy/beam
#Peak SNR: 1308.68

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

for robust in [2.0,1.0,0.5,0.0,-0.5,-1.0,-2.0]:
    imagename=prefix+'_SB_continuum_robust_'+str(robust)
    os.system('rm -rf '+imagename+'*')

    sigma = get_sensitivity(data_params, specmode='mfs')

    tclean_wrapper(vis=vislist, imagename=imagename, sidelobethreshold=2.0, 
            smoothfactor=1.5, scales=scales, threshold=3.0*sigma, 
            noisethreshold=3.0, robust=robust, parallel=parallel, 
            cellsize='0.025arcsec', imsize=1600,phasecenter=data_params['SB1']['common_dir'].replace('J2000','ICRS'))

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







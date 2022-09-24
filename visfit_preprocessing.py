############ A script for pre-processing of the visibility fit ##############
###                                                                       ###  
### 1. align the phase center of each execution block. This is really     ###
###    important to derive the correct visibility profile especially for  ###
###    the case of long baseline data.                                    ###
###                                                                       ###
### 2. split out the data without any flagged data. Depeding on the data  ###
###    size, one may consider some binning/averaging of the raw data.     ###
###                                                                       ###
#############################################################################

from casatasks import split, fixvis, fixplanets
import casatools
from eDisk_source_dict import source_dict
import os

ms = casatools.ms()

source = "L1489IRS"
datafilepath = "/works/yamato/eDisk/L1489IRS/ShortBaseline/eDisk_calibrated_data/"
# savefilepath = "/works/yamato/eDisk/L1489IRS/visibility_analysis/"
savefilepath = "/works/yamato/eDisk/L1489IRS/ShortBaseline/eDisk_calibrated_data/"
dataindex = ["SB1", "SB2", "LB1", "LB2"]
dataindex = ["SB1", "SB2", "SB3"]
# center_coord = "J2000 " + source_dict[source]["radec"]
center_coord = 'J2000 04h04m43.080s 26d18m56.104s'

### 1. align the phase center to the center derived from Gaussian fit
for i in dataindex:
    filename = "L1489IRS_{:s}_continuum.ms".format(i)
    vis = datafilepath + filename

    print(f"Shifting the phase center of dataset {i}...")
    fixvis(
        vis=vis,
        outputvis=savefilepath + filename.replace(".ms", "_shift.ms"),
        field=source,
        phasecenter=center_coord,
    )
    fixplanets(
        vis=savefilepath + filename.replace(".ms", "_shift.ms"),
        field=source,
        direction=center_coord,
    )
    print("Done.")

### 2. split out the data
# for i in dataindex:
#     vis = savefilepath + "L1489IRS_{:s}_continuum_shift.ms".format(i)

#     ms.open(vis)
#     metadata = ms.metadata()
#     width = [metadata.nchan(j) for j in metadata.spwfordatadesc()]
#     ms.close()

#     print(f"Splitting out the dataset {i}...")
#     outputvis = vis.replace(".ms", ".split.ms")
#     os.system("rm -r " + outputvis)
#     split(
#         vis=vis,
#         outputvis=outputvis,
#         datacolumn="data",
#         keepflags=False,
#         width=width,
#         # timebin="30s",
#         # combine="state,scan"
#     )
#     print("Done.")


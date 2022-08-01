import analysis_utils as au
import casatasks

source = "L1489IRS"
config = "SBLB"
line = "C18O"
robust = 1.0

imagename = au.customimagepath + au.get_image_basename(source, config, line, robust=robust, dv=0.2)

function = "boxcar"
width = 3
outfile = imagename.replace(".fits", ".specsmooth_width" + str(width))

# specsmooth
print("Start smoothing...")
casatasks.specsmooth(imagename, outfile=outfile, function=function, width=width)

print("Exporting to a fits...")
casatasks.exportfits(imagename=outfile, fitsimage=outfile + ".fits", dropdeg=True, overwrite=True)

import os
os.system("rm -r " + outfile)

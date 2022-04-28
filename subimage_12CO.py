from black import out
import casatasks
from analysis_utils import get_corresp_channels, get_spectral_coord
import astropy.io.fits as fits
import subprocess

molecular_species = ["12CO", "13CO", "C18O", "SO"]

for ms in molecular_species:

    imagename = "/raid/work/yamato/eDisk_data/L1489IRS/ALMA_pipeline_reduced_data/try1_continuum_nterms1/L1489IRS_SBLB_{:s}_robust_0.5.image.fits".format(
        ms
    )

    header = fits.getheader(imagename)

    v_range = (-12.0, 28.0)
    velax = get_spectral_coord(header, which="vel")
    chans = "~".join([str(get_corresp_channels(velax, v)) for v in v_range])
    print("Chosen channels: " + chans)

    print("Processing imsubimage...")
    outfile = imagename.replace(".fits", ".sub")
    subprocess.run(["rm", "-r", outfile])
    casatasks.imsubimage(imagename=imagename, outfile=outfile, chans=chans)
    print("Exporting into fits...")
    casatasks.exportfits(
        imagename=outfile,
        fitsimage=outfile + ".fits",
        overwrite=True,
        dropdeg=True,
    )

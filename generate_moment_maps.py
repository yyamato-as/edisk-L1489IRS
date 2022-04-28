from eDisk_source_dict import source_dict
from analysis_utils import  get_corresp_channels, get_spectral_coord
import bettermoments as bm
import astropy.io.fits as fits
import sys

ms = sys.argv[1]
source = "L1489IRS"
source_info = source_dict["L1489IRS"]
clip_th = (-1000, 3)

moment_maps = {}

# for i, ms in enumerate(molecular_species):
imagepath = "/raid/work/yamato/eDisk_data/L1489IRS/ALMA_pipeline_reduced_data/try1_continuum_nterms1/L1489IRS_SBLB_{:s}_robust_0.5.image.fits".format(ms)
if ms == "12CO":
    imagepath = imagepath.replace(".fits", ".sub.fits")


header = fits.getheader(imagepath)
# data = fits.getdata(imagepath).squeeze()
# velax = get_spectral_coord(header, which="vel")

data, velax = bm.load_cube(imagepath)

rms = bm.estimate_RMS(data=data, N=3)
print("rms = {} mJy/beam".format(rms*1e3))

firstchannel, lastchannel = [get_corresp_channels(velax, v*1e3) for v in source_info["emission_extent"][ms]]

data = data[firstchannel:lastchannel, :, :]
velax = velax[firstchannel:lastchannel]

mask = bm.get_threshold_mask(data=data, clip=clip_th, noise_channels=3)

data *= mask

method_func = {"zeroth": bm.collapse_zeroth, "first": bm.collapse_first, "eighth": bm.collapse_eighth}
for method in method_func.keys():
    print("Processing {:s} {:s} moment...".format(ms, method))
    moments = method_func[method](velax=velax, data=data, rms=rms)
    bm.save_to_FITS(moments=moments, method=method, path=imagepath)

    # command = "bettermoments " + imagepath
    # command += " -firstchannel {:s} -lastchannel {:s}".format(str(firstchannel), str(lastchannel))
    # command += " -noisechannels 3"
    # # command += " -clip {:s}".format(str(clip_th))

    # methods = ["zeroth", "first", "eighth"]
    # for m in methods:
    #     command += " -method " + m
    #     print("processing {:s} {:s} moment...".format(ms, m))
    #     subprocess.run(command, shell=True)

#     # plotting
#     for j, mom in enumerate([M0, M8]):
#         norm = ImageNormalize(mom, vmin=0.0)
#         pcolorfast_kw = {'cmap': cpal, 'norm': norm}
#         plot_2D_map(data=mom, X=x, Y=y, ax=axes[j, i], beam=beam, scale=scale, pcolorfast_kw=pcolorfast_kw, contour=False)
#         # norm = ImageNormalize(M8, vmin=0.0)
#         # plot_2D_map(data=M8, X=x, Y=y, ax=axes[1, i], title=ms, beam=beam, scale=scale, pcolorfast_kw=pcolorfast_kw, contour=False)

# for ax in axes.ravel():
#     ax.set(xlim=(5,-5), ylim=(-5, 5))
#     ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
#     ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        

# plt.show()

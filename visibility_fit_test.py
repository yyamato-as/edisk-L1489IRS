import matplotlib.pyplot as plt
import numpy as np
from galario.double import sampleProfile, get_image_size
import astropy.constants as ac
import astropy.units as units
import pickle
import corner
import multiprocessing
from visfit import load_ms
import emcee
from mcmc_tools import EmceeHammer, log_prior
from fileio import import_mms

c = ac.c.to(units.m / units.s).value
pi = np.pi
deg = pi / 180.0  # in rad
arcsec = pi / 180.0 / 3600.0  # in rad

def set_grid(u, v, rmax, verbose=True):
    # gridding parameter
    # print(np.min(np.sqrt(u**2+v**2)))
    nxy, dxy = get_image_size(u, v, verbose=verbose)  # in rad
    # print(nxy, dxy/arcsec)

    # condition for GALARIO interpolation: dxy/2 - Rmin > 5 * dR
    # not to make dR too small, adopt dxy/2/1000 as Rmin
    rmin = dxy / 2.0 * 1e-3  # in rad
    # print(rmin/arcsec)
    # dR = (dxy/2 - Rmin) / 5.0001
    dr = (dxy / 2.0 - rmin) / 5.0001  # in rad
    # dr = 0.01*arcsec
    # print(dr/arcsec)

    r_pad = 2 * dxy  # padding parameter
    rmax = dxy * nxy / np.sqrt(2) + r_pad if rmax is None else rmax # in rad
    # rmax = 10*arcsec
    # print(rmax/arcsec)

    # radial grid on image plane
    r = np.arange(rmin, rmax, dr)  # in rad

    return nxy, dxy, r, rmin, dr

def setup_params(param_dict):

    param_name = [
        key for key in param_dict.keys() if not param_dict[key]["fixed"]
    ]
    fixed_param_name = [
        key for key in param_dict.keys() if param_dict[key]["fixed"]
    ]
    bound = [
        param_dict[key]["bound"]
        for key in param_dict.keys()
        if not param_dict[key]["fixed"]
    ]
    initial_state = [
        param_dict[key]["p0"]
        for key in param_dict.keys()
        if not param_dict[key]["fixed"]
    ]

    return param_name, fixed_param_name, bound, initial_state

def Gaussian1d(r, F, sigma):
    return 10 ** F / (np.sqrt(2 * np.pi) * sigma) * np.exp(-r**2 / (2 * sigma**2))

def GaussianRing1d(r, r0, F, sigma):
    return 10 ** F / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(r - r0)**2 / (2 * sigma**2))



geometrical_param = ["PA", "incl", "dRA", "dDec"]
gp_default = {
    "PA": {"p0": 0.0, "bound": (0.0, 180.0), "fixed": True},
    "incl": {"p0": 0.0, "bound": (0.0, 90.0), "fixed": True},
    "dRA": {"p0": 0.0, "bound": (-5.0, 5.0), "fixed": True},
    "dDec": {"p0": 0.0, "bound": (-5.0, 5.0), "fixed": True},
}
accepted_chain_length = 100


################## Target name ###################
source = "L1489IRS"
dataset = "SBLB"
################## Define the model function and parameters ##################
# model function on the image plane; modify as needed
# def model_func_1d(r, I_g, sigma_g, I_p, sigma_p):
#     return 10**I_g * np.exp(-0.5 * r**2 / sigma_g**2) + 10**I_p * np.exp(
#         -0.5 * r**2 / sigma_p**2
#     )


# def model_func_1d(r, I_p, sigma_p, I_r, r_r, sigma_r):
#     return 10**I_p * np.exp(-0.5 * r**2 / sigma_p**2) + 10**I_r * np.exp(
#         -0.5 * (r - r_r) ** 4 / sigma_r**4
#     )


# model function derived from image fit
model_name = "PointSource_GaussianRing_Gaussian"


def model_func_1d(r, F_c, sigma_c, r0_r, F_r, sigma_r, F_b, sigma_b):
    return (
        Gaussian1d(r, F_c, sigma_c)
        + GaussianRing1d(r, r0_r, F_r, sigma_r)
        + Gaussian1d(r, F_b, sigma_b)
    )


# dictionary of parameters including initial guess, bound, and fixed/free
param_dict = {
    # "I_g": {"p0": 8.5, "bound": (-2.0, 20.0), "fixed": False},
    # "sigma_g": {"p0": 2.0, "bound": (0.1, 10), "fixed": False},
    "F_c": {"p0": 10.59, "bound": (5, 15), "fixed": False},
    "sigma_c": {"p0": 0.01, "bound": (1e-5, 0.2), "fixed": False},
    "r0_r": {"p0": 0.47, "bound": (0.2, 0.7), "fixed": False},
    "F_r": {"p0": 8.31, "bound": (3, 13), "fixed": False},
    "sigma_r": {"p0": 0.20, "bound": (1e-2, 1.5), "fixed": False},
    "F_b": {"p0": 9.14, "bound": (4, 14), "fixed": False},
    "sigma_b": {"p0": 1.59, "bound": (0.3, 5), "fixed": False},
    "PA": {"p0": 68.9, "bound": (0, 180), "fixed": False},
    "incl": {"p0": 72.86, "bound": (0, 90), "fixed": False},
    "dRA": {"p0": 0.0, "bound": (-2, 2), "fixed": False},
    "dDec": {"p0": 0.0, "bound": (-2, 2), "fixed": False},
}
###############################################################################


########## load the visibility and define various functions for MCMC ##########
# load the visibility
print("Loading visibility data...")
# datafilepath = "/raid/work/yamato/edisk_data/edisk_calibrated_data_old/DDT_LP/"
# # uvtabfilename = datafilepath + "L1489IRS_continuum.ms.uvtab"
# # visibility = np.loadtxt(uvtabfilename, unpack=True)
# # u, v, real, imag, weight, freqs = np.ascontiguousarray(visibility)
# path = "/works/yamato/eDisk/L1489IRS/ALMA_pipeline_calibrated_data/"
# vislist = [path + "L1489IRS_{:s}_continuum_shift.bin_30s.split.ms".format(i) for i in ["SB1", "SB2", "SB3", "LB1", "LB2"]]
# u, v, real, imag, weight, freqs = import_mms(vislist=vislist, remove_flag=True)
# filename = "./oldfileiotest_bin30s.npz"
# np.savez(filename, u=u, v=v, real=real, imag=imag, weight=weight, freqs=freqs)

# data = np.load("./visibility/L1489IRS_continuum_shift_oldDDT.split.npz")
data = np.load("./visibility/L1489IRS_continuum_shift.split.bin_30s.npz")
u = data["u"]
v = data["v"]
real = data["V"].real
imag = data["V"].imag
weight = data["weight"]
# freqs = data["freqs"]
print("Done.")
# print(np.all(np.isnan(u)))


# datafilepath = "/raid/work/yamato/edisk_data/edisk_calibrated_data/"
# uvtabfilename = datafilepath + "L1489IRS_continuum.bin_30s.ms.uvtab"
# visibility = np.loadtxt(uvtabfilename, unpack=True)
# u, v, real, imag, weight, freqs = np.ascontiguousarray(visibility)
# print("Done.")

# datafilepath = "/raid/work/yamato/edisk_data/edisk_calibrated_data/"
# msfilename = datafilepath + "L1489IRS_continuum.bin_30s.ms"
# u, v, real, imag, weight, freqs = import_ms(msfilename)
# path = "/works/yamato/eDisk/L1489IRS/ALMA_pipeline_calibrated_data/"
# msfilenames = [path + f"L1489IRS_{i}_continuum_shift.split.ms" for i in ["SB1", "SB2", "SB3", "LB1", "LB2"]]
# vis = load_ms(msfilenames)
# vis = np.load("./visibility/L1489IRS_continuum_shift_bin_10klambda.npz")
# u = vis["u"]
# v = vis["v"]
# # data = vis["data"]
# data = vis["V"]
# weight = vis["weight"]

# get gridding parameters
nxy, dxy, r, rmin, dr = set_grid(u, v, rmax=None)
# nxy = 1024
# dxy = 0.02 * arcsec
# rmin = 0.0005 * arcsec
# dr = 0.001 * arcsec
# rmax = 15 * arcsec
# r = np.arange(rmin, rmax, dr)

# setup the fitting parameters
param_name, fixed_param_name, bound, initial_state = setup_params(param_dict)

# function for visibility sampling
def sample_vis(param_dict):

    # retrieve geometrical params
    PA = param_dict.pop("PA", gp_default["PA"]["p0"])
    incl = param_dict.pop("incl", gp_default["incl"]["p0"])
    dRA = param_dict.pop("dRA", gp_default["dRA"]["p0"])
    dDec = param_dict.pop("dDec", gp_default["dDec"]["p0"])

    # get model array
    model = model_func_1d(r / arcsec, **param_dict)

    # sampling by GALARIO
    V = sampleProfile(
        intensity=model,
        Rmin=rmin,
        dR=dr,
        nxy=nxy,
        dxy=dxy,
        u=u,
        v=v,
        dRA=dRA * arcsec,
        dDec=dDec * arcsec,
        PA=PA * deg,
        inc=incl * deg,
        check=False,
    )

    return V


# likelihood function
def log_likelihood(param):

    param_dict = {name: p for name, p in zip(param_name, param)}

    # update fixed param
    param_dict.update({name: param_dict[name]["p0"] for name in fixed_param_name})

    # sample model visibility
    # print("start")
    model_vis = sample_vis(param_dict)
    # print("stop")

    # compute log likelihood
    rms = np.sqrt(1.0 / weight)
    ll = -0.5 * np.sum(
        ((model_vis.real - real) ** 2 + (model_vis.imag - imag) ** 2) / rms**2
        + np.log(2 * pi * rms**2)
    )

    return ll


# probability function
def log_probability(param):
    lp = log_prior(param, bound)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(param)
    return lp + ll


###############################################################


################### Actual MCMC run ###########################
nwalker = 32
nstep = 5000
progress = True


hammer = EmceeHammer(log_probability=log_probability, initial_state=initial_state, nwalker=nwalker, nstep=nstep)
with multiprocessing.Pool(processes=8) as pool:
# pool = None
    hammer.run(pool=pool, savefilename="L1489IRS_continuum_shift.split.bin_30s.GaussianSingleRing.h5", save=True)
"""

# # with multiprocessing.Pool(processes=8) as pool:
# pool = None
# sampler = run_emcee(
#     log_probability=log_probability,
#     initial_state=initial_state,
#     nwalker=nwalker,
#     nstep=nstep,
#     pool=pool,
#     progress=progress,
#     # nthreads=8
# )

# save the EmsambleSampler object into a pickle
# with open("./L1489IRS_continuum_sampler_{:s}.pkl".format(model_name), "wb") as f:
#     pickle.dump(sampler, f, protocol=pickle.HIGHEST_PROTOCOL)


# load back the sampler
from mcmc_tools import EmceeHammer
nburnin = 4000

hammer = EmceeHammer()
hammer.load_backend("./L1489IRS_continuum_shift.split.bin_30s.GaussianSingleRing.h5")
hammer.plot_corner(nburnin=nburnin, labels=param_name, title_fmt=".5f")
# hammer.plot_walker(nburnin=nburnin, labels=param_name)
plt.show()


MAP_params = hammer.get_MAP_params(nburnin=nburnin)
print(MAP_params)
param_dict = {name: p for name, p in zip(param_name, MAP_params)}

# # update fixed param
# param_dict.update({name: param_dict[name]["p0"] for name in fixed_param_name})

# # best-fit vis
# MAP_vis = sample_vis(param_dict.copy())
# filename = "./L1489IRS_continuum_shift_oldDDT.split.GaussianSingleRing.best.npz"
# np.savez(filename, u=u, v=v, V=MAP_vis, weight=weight)


# from visibility import Visibility

# modelvis = Visibility(filename)
# modelvis.deproject(PA=param_dict.get("PA"), incl=param_dict.get("incl"))
# modeluv, modelV, modelerr = modelvis.bin_1D()

# filename = "./visibility/L1489IRS_continuum_shift_oldDDT.split.npz"
# obsvis = Visibility(filename)
# obsvis.deproject(PA=param_dict.get("PA"), incl=param_dict.get("incl"))
# obsuv, obsV, obserr = obsvis.bin_1D()

# plt.figure()
# plt.errorbar(obsuv, obsV.real, yerr=obserr, fmt="o")
# plt.plot(modeluv, modelV)
# plt.axhline(y=0.0, color="grey", ls="dashed")
# plt.xscale("log")
# plt.xlabel("uv distance [lambda]")
# plt.ylabel("Real [Jy]")
# plt.show()





"""






"""
################################################################


########## plot various figures #############

### TODO: implement the autocorrelation analysis ###
# autocorrelation analysis
# tau = sampler.get_autocorr_time()
# print("Autocorreation time: {:g}".format(tau))

# nburnin = int(3 * tau)
# thin = int(tau / 2)

# load back the EnsambleSampler


with open(
    "./L1489IRS_{:s}_continuum_sampler_{:s}.pkl".format(dataset, model_name), "rb"
) as f:
    sampler = pickle.load(f)

# get the flatted and discarded sample
nburnin = 700
thin = 1

sample = sampler.get_chain()
sample_flat = sampler.get_chain(flat=True)
sample_flat_disc = sampler.get_chain(discard=nburnin, thin=thin, flat=True)


# corner plot
corner_fig = corner.corner(
    sample_flat_disc,
    labels=param_name,
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
)

# chain plot
walker_figs = plot_walker(sample=sample, nburnin=nburnin, labels=param_name)

# overplot MAP model
MAP_param = sample_flat[np.argmax(sampler.get_log_prob(flat=True))]
param_dict = {name: p for name, p in zip(param_name, MAP_param)}

# update fixed param
param_dict.update({name: param_dict[name]["p0"] for name in fixed_param_name})

with open("./L1489IRS_{:s}_continuum_{:s}_MAP_param.npy".format(dataset, model_name), "wb") as f:
    pickle.dump(param_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

# get the visibility for MAP model
MAP_vis = sample_vis(param_dict.copy())

np.save("./L1489IRS_{:s}_continuum_MAP_vis_{:s}_noflag.npy".format(dataset, model_name), MAP_vis)

# real/imag vs. baseline plot for observe values
model_fig = plot_uv(
    u=u,
    v=v,
    real=real,
    imag=imag,
    weight=weight,
    incl=param_dict["incl"],
    PA=param_dict["PA"],
    fmt="o",
    capsize=3,
    markersize=5,
    zorder=-100,
    binsize=50e3,
    uvrange=(0, 1000e3),
)

# deprojection and binning for model visibility
_, _, uvdist_deproj = deproject_vis(u, v, incl=param_dict["incl"], PA=param_dict["PA"])
uvdist_deproj_binned, real_binned, imag_binned = bin_vis(
    uvdist=uvdist_deproj,
    real=MAP_vis.real,
    imag=MAP_vis.imag,
    bins=np.arange(0, 1000e3, 5e3),
    weighted_average=False,
)

# plot model visibility
for ax, V in zip(model_fig.axes, [real_binned, imag_binned]):
    ax.plot(uvdist_deproj_binned / 1e3, V)
    ax.grid()

# ancillary stuffs
model_fig.axes[0].set(ylim=(0.0, 0.05), xlim=(0, 1000))
model_fig.axes[1].set(ylim=(-0.0055, 0.0055), xlim=(0, 1000))

plt.show()


################## save the figures into a directory #####################
import subprocess


# figpath = "./fig_{:s}_{:s}_{:s}/".format(source, dataset, model_name)
# subprocess.run(["mkdir", figpath])
# corner_fig.savefig(figpath + "corner.png", bbox_inches="tight", dpi=300)
# for i, fig in enumerate(walker_figs):
#     fig.savefig(
#         figpath + "walker_{}.png".format(param_name[i]),
#         bbox_inches="tight",
#         dpi=300,
#     )
# model_fig.savefig(figpath + "uvplot.png", bbox_inches="tight", dpi=300)
"""
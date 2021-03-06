import numpy as np
import astropy.constants as ac
import astropy.units as u
import emcee
from galario.double import get_image_size, check_image_size, sampleProfile
from uvplot import UVTable, export_uvtable, COLUMNS_V0
import casatools
import casatasks
import subprocess

c = ac.c.to(u.m / u.s).value
pi = np.pi
deg = pi / 180.0  # in rad
arcsec = pi / 180.0 / 3600.0  # in rad

geometrical_param = ["PA", "incl", "dRA", "dDec"]
gp_default = {
    "PA": {"p0": 90.0, "bound": (0.0, 180.0), "fixed": True},
    "incl": {"p0": 45.0, "bound": (0.0, 90.0), "fixed": True},
    "dRA": {"p0": 0.0, "bound": (-5.0, 5.0), "fixed": True},
    "dDec": {"p0": 0.0, "bound": (-5.0, 5.0), "fixed": True},
}


tb = casatools.table()
ms = casatools.ms()


class VisFit1D:
    def __init__(self, obsmsfile, verbose=False):

        self.obsmsfile = obsmsfile
        self.verbose = verbose

    def load_vis(self, split_kwargs={"datacolumn": "data"}):
        # get number of channels for each spw
        tb.open(self.obsmsfile + "/SPECTRAL_WINDOW")
        self.nchans = tb.getcol("NUM_CHAN")
        tb.close()

        # split out to one channel per an spw  ### TODO: reconsider if this is proper: don't we need to consider the frequency differences between channels? ###
        tempvis = self.obsmsfile + ".split.temp"
        subprocess.run(["rm", "-rf", tempvis])
        casatasks.split(vis=self.obsmsfile, outputvis=tempvis, keepflags=False, width=self.nchans, **split_kwargs)

        tb.open(tempvis)

        # get uvw coordinates
        self.u, self.v, self.w = np.require(tb.getcol("UVW"), requirements="C")

        # get weight
        weight_ori = tb.getcol("WEIGHT")

        datacolumn = split_kwargs.get("datacolumn", "data")
        if datacolumn.upper() in tb.colnames():
            data = tb.getcol(datacolumn.upper())
        else:
            raise KeyError("datacolumn {} not found.".format(datacolumn))

        # average over dual polarization
        V_XX = data[0, 0, :]
        V_YY = data[1, 0, :]
        weight_XX = weight_ori[0, :]
        weight_YY = weight_ori[1, :]

        self.V = np.require((V_XX * weight_XX + V_YY * weight_YY) / (weight_XX + weight_YY), requirements="C")
        self.weight = np.require(weight_XX + weight_YY, requirements="C")

        tb.close()

        # get frequency
        tb.open(tempvis + "/SPECTRAL_WINDOW")
        freqs = tb.getcol("CHAN_FREQ")  # Hz
        tb.close()

        self.obs_nu = freqs.mean()
        self.obs_wle = c / freqs.mean()  # m

        # remove temporary ms file
        subprocess.run(["rm", "-rf", tempvis])

    def export_vis(self, filename):
        np.savetxt(
            filename,
            np.column_stack([self.u, self.v, self.V.real, self.V.imag, self.weight]),
            fmt="%10.6e",
            delimiter="\t",
            header="Extracted from {}.\nwavelength[m] = {}\nColumns:\tu[m]\tv[m]\tRe(V)[Jy]\tIm(V)[Jy]\tweight".format(
                self.obsmsfile, self.obs_wle
            ),
        )

    def fit_vis(
        self, model_func_1d, param_dict, dish_D=12, nwalker=32, nstep=500, pool=None, blobs_dtype=None, progress=True
    ):
        # setup radial grid
        self.set_grid(dish_D=dish_D)

        self.param_dict = param_dict

        # add required geometrical params as needed
        # for gp in geometrical_param:
        #     if not gp in self.param_dict:
        #         self.param_dict.update(gp=gp_default[gp])

        # get lists of free/fixed parameter names
        self.param_name = [key for key in self.param_dict.keys() if not self.param_dict[key]["fixed"]]
        self.fixed_param_name = [key for key in self.param_dict.keys() if self.param_dict[key]["fixed"]]

        def log_likelihood(param):

            param_dict = {name: p for name, p in zip(self.param_name, param)}

            # update fixed param
            param_dict.update({name: self.param_dict[name]["p0"] for name in self.fixed_param_name})

            # sample model visibility
            model_vis = self.sample_vis(model_func_1d, param_dict)

            # compute log likelihood
            rms = np.sqrt(1.0 / self.weight)
            ll = -0.5 * np.sum(
                ((model_vis.real - self.V.real) ** 2 + (model_vis.imag - self.V.imag) ** 2) / rms ** 2
                + np.log(2 * pi * rms ** 2)
            )

            return ll

        def log_probability(param):
            lp = log_prior(param, self.bound)
            if not np.isfinite(lp):
                return -np.inf
            ll = log_likelihood(param)
            return lp + ll

        self.bound = [
            self.param_dict[key]["bound"] for key in self.param_dict.keys() if not self.param_dict[key]["fixed"]
        ]
        self.initial_state = [self.param_dict[key]["p0"] for key in self.param_dict.keys() if not self.param_dict[key]["fixed"]]

        ndim = len(self.initial_state)
        p0 = self.initial_state + 1e-4 * np.random.randn(nwalker, ndim)

        # set smapler
        self.sampler = emcee.EnsembleSampler(
            nwalker, ndim, log_probability, pool=pool, blobs_dtype=blobs_dtype
        )

        # run
        print(
            "starting to run the MCMC sampling with: \n \t initial state:",
            self.initial_state,
            "\n \t number of walkers:",
            nwalker,
            "\n \t number of steps:",
            nstep,
        )
        self.sampler.run_mcmc(p0, int(nstep), progress=progress)

        

    def sample_vis(self, model_func_1d, param_dict):

        # retrieve geometrical params
        PA = param_dict.pop("PA", gp_default["PA"]["p0"])
        incl = param_dict.pop("incl", gp_default["incl"]["p0"])
        dRA = param_dict.pop("dRA", gp_default["dRA"]["p0"])
        dDec = param_dict.pop("dDec", gp_default["dDec"]["p0"])

        # get model array
        model = model_func_1d(self.r, **param_dict)

        # sampling by GALARIO
        vis = sampleProfile(
            intensity=model,
            Rmin=self.rmin,
            dR=self.dr,
            nxy=self.nxy,
            dxy=self.dxy,
            u=self.u,
            v=self.v,
            dRA=dRA * arcsec,
            dDec=dDec * arcsec,
            PA=PA * deg,
            inc=incl * deg,
            check=False,
        )

        return vis

    def set_grid(self, dish_D=12):
        # gridding parameter
        self.nxy, self.dxy = get_image_size(
            self.u, self.v, PB=1.22 * self.obs_wle / dish_D, verbose=self.verbose
        )  # in rad

        # condition for GALARIO interpolation: dxy/2 - Rmin > 5 * dR
        # not to make dR too small, adopt dxy/2/1000 as Rmin
        self.rmin = self.dxy / 2.0 * 1e-3  # in rad

        # dR = (dxy/2 - Rmin) / 5.1
        self.dr = (self.dxy / 2.0 - self.rmin) / 5.1  # in rad

        r_pad = 2 * self.dxy  # padding parameter
        self.rmax = self.dxy * self.nxy / np.sqrt(2) + r_pad  # in rad

        # radial grid in image plane
        self.r = np.arange(self.rmin, self.rmax, self.dr)  # in rad

def condition(p, b):
    if (p >= b[0]) and (p <= b[1]):
        return True
    return False

def log_prior(param, bound):
    for p, b in zip(param, bound):
        if not condition(p, b):
            return -np.inf
    return 0.0

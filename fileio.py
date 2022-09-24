from os import remove
import casatools 
import casatasks
import numpy as np
import shutil
import astropy.constants as ac
import astropy.units as u

c = ac.c.to(u.m/u.s).value

tb = casatools.table()
ms = casatools.ms()

def import_mms(vislist, export_uvtable=False, filename=None, remove_flag=False):
    """Import the visibility from multiple CASA measurement sets as 1D arrays using casatools. 
    Parameters
    ----------
    msfilename : str
        The measurement set filename which you want to import.
    export_uvtable : bool, optional
        If write the output into a text file (similar to UVTable), by default True
    filename : str or None, optional
        The filename of output text file, by default None. Relevant only when export_uvtable=True.
    Returns
    -------
    u, v, real, imag, weight, freqs: six 1D numpy arrays
        visibility
    """

    print("Will import {}...".format(vislist))

    u = []
    v = []
    real = []
    imag = []
    weight = []
    freqs = []

    for vis in vislist:
        _u, _v, _real, _imag, _weight, _freqs = import_ms(vis, export_uvtable=False, remove_flag=remove_flag)
        u.append(_u)
        v.append(_v)
        real.append(_real)
        imag.append(_imag)
        weight.append(_weight)
        freqs.append(_freqs)

    # concatenate 
    u = np.ascontiguousarray(np.concatenate(u))
    v = np.ascontiguousarray(np.concatenate(v))
    real = np.ascontiguousarray(np.concatenate(real))
    imag = np.ascontiguousarray(np.concatenate(imag))
    weight = np.ascontiguousarray(np.concatenate(weight))
    freqs = np.ascontiguousarray(np.concatenate(freqs))

    print("Done.")

    if export_uvtable:
        if filename is None:
            filename = "output.uvtab"

        np.savetxt(
            filename,
            np.column_stack([u, v, real, imag, weight, freqs]),
            fmt="%10.6e",
            delimiter="\t",
            header="u [lambda]\t v [lambda]\t real [Jy] \t imag [Jy]\t weight \t nu [Hz]",
        )

    return u, v, real, imag, weight, freqs



def import_ms(msfilename, export_uvtable=True, filename=None, remove_flag=False):
    """Import the visibility from a CASA measurement set as 1D arrays using casatools. 
    Parameters
    ----------
    msfilename : str
        The measurement set filename which you want to import.
    export_uvtable : bool, optional
        If write the output into a text file (similar to UVTable), by default True
    filename : str or None, optional
        The filename of output text file, by default None. Relevant only when export_uvtable=True.
    Returns
    -------
    u, v, real, imag, weight, freqs: six 1D numpy arrays
        visibility
    """

    ### remove flagged data just in case
    if remove_flag:
        casatasks.split(vis=msfilename, outputvis=msfilename+".split", keepflags=False, datacolumn="data")
        # subprocess.run("rm -r " + outms, shell=True)
        # subprocess.run("mv " + outms+".tmp " + outms, shell=True)

        msfilename += ".split"

    ms.open(msfilename)

    print("Loading {:s}...".format(msfilename))

    spw = [key for key in ms.getspectralwindowinfo()]

    data = {}
    for i in spw:
        ms.selectinit(datadescid=int(i))
        #ms.selectpolarization(corr)
        data[i] = ms.getdata(["u" ,"v", "data", "weight", "axis_info", "flag"])
        ms.reset()

    ms.close()

    # manipulate read visibilities
    u = []
    v = []
    V = []
    weight = []
    freqs = []

    for spw in data.keys():
        print("spw" + spw, data[spw]["data"].shape)
        
        # average over polarization
        if data[spw]["data"].shape[0] == 2:
            _V = np.sum(data[spw]["data"]*data[spw]["weight"][:,None,:], axis=0) / np.sum(data[spw]["weight"], axis=0)
            _weight = np.sum(data[spw]["weight"], axis=0)
            good = np.any(data[spw]["flag"], axis=0) == False
            # V_XX = data[spw]["data"][0,:,:]
            # V_YY = data[spw]["data"][1,:,:]
            # weight_XX = data[spw]["weight"][0,:]
            # weight_YY = data[spw]["weight"][1,:]
            # _weight = weight_XX + weight_YY
            # _V = (V_XX * weight_XX + V_YY * weight_YY) / _weight

        else:
            _weight = data[spw]["weight"].squeeze()
            _V = data[spw]["data"].squeeze()
            good = data[spw]["flag"] == False

        nchan, nuv = _V.shape
        _freqs = data[spw]["axis_info"]["freq_axis"]["chan_freq"]

        _freqs = np.tile(_freqs, nuv)
        _wles = c / _freqs
        _u = np.tile(data[spw]["u"], (nchan, 1)) / _wles # in lmabda
        _v = np.tile(data[spw]["v"], (nchan, 1)) / _wles # in lambda
        _weight = np.tile(_weight, (nchan, 1))

        # remove flagged data
        # _u = _u[good]
        # _v = _v[good]
        # _V = _V[good]
        # _weight = _weight[good]
        # _freqs = _freqs[good]

        # remove the autocorrelation; here all the variables are flattened
        # -> but this removement causes an issue if you want the model visibility to get back the original ms file after fitting... so quit
        # cc = uvdist != 0

        # u = u[cc]
        # v = v[cc]
        # V = V[cc]
        # weight = weight[cc]
        # freqs = freqs[cc]

        # append each component with flatten
        u.append(_u.ravel())
        v.append(_v.ravel())
        V.append(_V.ravel())
        weight.append(_weight.ravel())
        freqs.append(_freqs.ravel())

    # concatenate 
    u = np.ascontiguousarray(np.concatenate(u))
    v = np.ascontiguousarray(np.concatenate(v))
    V = np.concatenate(V)
    real = np.ascontiguousarray(V.real)
    imag = np.ascontiguousarray(V.imag)
    weight = np.ascontiguousarray(np.concatenate(weight))
    freqs = np.ascontiguousarray(np.concatenate(freqs))

    print("Done.")

    if export_uvtable:
        if filename is None:
            filename = "output.uvtab"

        np.savetxt(
            filename,
            np.column_stack([u, v, real, imag, weight, freqs]),
            fmt="%10.6e",
            delimiter="\t",
            header="u [lambda]\t v [lambda]\t real [Jy] \t imag [Jy]\t weight \t nu [Hz]",
        )

    return u, v, real, imag, weight, freqs


def get_number_of_data_point(msfilename, include_pol=False):
    ms.open(msfilename)
    md = ms.metadata()

    ndatapoint = {}
    for i in md.spwfordatadesc():
        ms.selectinit(datadescid=i)
        ndata = md.nchan(spw=i)*ms.nrow(selected=True)
        if include_pol:
            ndata *= md.ncorrforpol(polid=0)
        ndatapoint[i] = ndata
        ms.reset()
    ms.close()

    return ndatapoint


def export_mms(basemslist, outmslist, real, imag, weight):

    for basems, outms in zip(basemslist, outmslist):
        # print(real.size)
        # shutil.copytree(basems, outms)
        # casatasks.split(vis=outms, outputvis=outms+".tmp", keepflags=False, datacolumn="data")
        # subprocess.run("rm -r " + outms, shell=True)
        # subprocess.run("mv " + outms+".tmp " + outms, shell=True)

        # ndatapoint = np.sum([i for i in get_number_of_data_point(outms).values()])
        # print(ndatapoint)
        # good_array = get_good_array(outms)
        # ngood = np.sum(good_array)
        ndatapoint = get_number_of_data_point(basems)
        ndata = np.sum([i for i in ndatapoint.values()])
        

        # split out relevant data
        _real, real = np.split(real, [ndata])
        _imag, imag = np.split(imag, [ndata])
        _weight, weight = np.split(weight, [ndata])
        # print(_real.size, real.size)

        # export
        export_ms(basems, outms, _real, _imag, _weight)

    # check all data are exported
    assert real.size == imag.size == weight.size == 0

    return


def export_ms(basems, outms, real, imag, weight):
    """Export the visibility (1D arrays) to a CASA measurement set file
    Parameters
    ----------
    basems : str
        Measurement set filename from which you get the data by import_ms function. Must have the same number of data points as input.
    outms : str
        The filename to which you want to export.
    real : 1D numpy array
        Real part of the visibility.
    imag : 1D numpy array
        Imaginary part of the visibility.
    weight : 1D numpy array
        Weight.
    """

    shutil.copytree(basems, outms)

    assert real.shape == imag.shape == weight.shape
    assert real.ndim == imag.ndim == weight.ndim == 1

    datasize = len(real)

    ms.open(outms, nomodify=False)

    spw = [key for key in ms.getspectralwindowinfo()]

    # check the number of datapoint consistency
    ndata = {}
    nchan = {}
    spw_array = []
    for i in spw:
        ms.selectinit(datadescid=int(i))
        #ms.selectpolarization(corr)

        rec = ms.getdata(["data"])

        nc = rec["data"].shape[1]
        nd = rec["data"].shape[2]

        nchan[i] = nc
        ndata[i] = nc * nd

        spw_array.append(np.full(ndata[i], int(i)))

        ms.reset()

    spw_array = np.concatenate(spw_array)

    if datasize != np.sum([i for i in ndata.values()]):
        raise ValueError("Data size is not consistent with the base measurement set.")

    # put the data onto each spectral window
    print("Exporting into {:s}...".format(outms))

    V = real + imag*1.0j

    for i in spw:
        print("processing spw="+i)
        ms.selectinit(datadescid=int(i))
        #ms.selectpolarization(corr)

        rec = ms.getdata(["data", "weight"])

        # print(rec["data"].shape, rec["weight"].shape)

        v = V[spw_array == int(i)].reshape(nchan[i], -1)
        w = weight[spw_array == int(i)].reshape(nchan[i], -1)

        if rec["data"].shape[0] == 2 and rec["data"].ndim == 3:
            #rec["data"] = np.tile(v, (2, 1, 1))
            rec["data"][0,:,:] = v 
            rec["data"][1,:,:] = v
            #rec["weight"] = np.tile(np.mean(w, axis=0), (2, 1, 1))
            rec["weight"][0,:] = np.mean(w, axis=0) 
            rec["weight"][1,:] = np.mean(w, axis=0)

        else:
            rec["data"][0,:,:] = v
            rec["weight"][0,:] = np.mean(w, axis=0)

        # print(rec["data"].shape, rec["weight"].shape)

        ms.putdata(rec)
        ms.reset()

    ms.close()

    print("Writing done.")



if __name__ == '__main__':
    import numpy as np
    # from fileio import export_ms

    datafilepath = "/raid/work/yamato/edisk_data/edisk_calibrated_data/"
    msfilename = datafilepath + "L1489IRS_SB1_continuum.bin_30s.ms"
    u, v, real, imag, weight, freqs = import_ms(msfilename, export_uvtable=True, filename=msfilename+".uvtab")

    # MAP_vis = np.load("./L1489IRS_SB1_continuum_PointSource_GaussianRing_Gaussian_MAP_vis.npy")

    export_ms(
        basems="/raid/work/yamato/edisk_data/edisk_calibrated_data/L1489IRS_SB1_continuum.bin_30s.ms",
        outms="/raid/work/yamato/edisk_data/edisk_calibrated_data/L1489IRS_SB1_continuum.bin_30s.model.ms",
        real=real,
        imag=imag,
        weight=np.ones(real.shape),
)
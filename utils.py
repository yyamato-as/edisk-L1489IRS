import subprocess
import os
import glob
import re
import pandas as pd
import numpy as np
import warnings
from casaplotms import plotms

hdr_col = {'scan': (3,), 'field': (4,), 'spw': (11,), 'baseline': (7, 8,), 'antenna': (7,), 'time': (9,), 'corr': (12,),}

# see https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
def atoi(text):
    return int(text) if text.isdigit() else text

# see https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


def pycsplit(inputfile, outname=None, split_str=None, niter=None, verbose=False):

    fn, ext = os.path.splitext(inputfile)

    if outname is None:
        outname = fn

    if split_str is None:
        split_str = '#'

    if niter is None:
        niter = '*'

    cmd = 'csplit'
    cmd += ' -s -f {:s} -b "_%u{:s}"'.format(outname, ext)
    cmd += ' {:s}'.format(inputfile)
    cmd += ' "/{:s}/"'.format(split_str)
    cmd += ' {' + niter + '}'
    
    subprocess.call(cmd, shell=True)
    
    fl = sorted(glob.glob('{:s}_*{:s}'.format(outname, ext)), key=natural_keys) # numerical order

    return fl


import matplotlib.pyplot as plt

class pyplotms:
    
#     plotms_header = 'x y chan scan field ant1 ant2 ant1name ant2name time freq spw corr obs'.split()
    

    def __init__(self,
        vis, maxnplots=100, plotms_kw={"xaxis": "frequency", "yaxis": "amplitude"}, verbose=True
    ):
        self.datafile_prefix = vis 
        self.maxnplots = maxnplots
        self.plotms_kw = plotms_kw
        self.xaxis = self.plotms_kw.get("xaxis")
        self.yaxis = self.plotms_kw.get("yaxis")
        self.misc = self.plotms_kw.get("showatm", False) or self.plotms_kw.get("showtsky", False) or self.plotms_kw.get("showimage", False)
        self.verbose = verbose
       
    def run_plotms(self,):
        # preremove the datafile
        #os.system("rm {:s}.txt".format(self.datafile_prefix))
        subprocess.run("rm {:s}.txt".format(self.datafile_prefix), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        if self.verbose:
            print("Running plotms...")
        # get the data from plotms in CASA
        plotms(
            vis=self.datafile_prefix,
            gridrows=self.maxnplots,  # need to set due to plotms bug... typically 100 is maxmum limit for plotms to work
            gridcols=1,
            xaxis=self.xaxis,
            xdatacolumn=self.plotms_kw.get("xdatacolumn", "data"),
            xframe=self.plotms_kw.get("xframe", "icrs"),
            xinterp=self.plotms_kw.get("xinterp", "cubic spline"),
            yaxis=self.yaxis,
            ydatacolumn=self.plotms_kw.get("ydatacolumn", "data"),
            yframe=self.plotms_kw.get("yframe", "icrs"),
            yinterp=self.plotms_kw.get("yinterp", "cubic spline"),
            selectdata=True,
            field=self.plotms_kw.get("field", ""),
            spw=self.plotms_kw.get("spw", ""),
            timerange=self.plotms_kw.get("timerange", ""),
            uvrange=self.plotms_kw.get("uvrange", ""),
            antenna=self.plotms_kw.get("antenna", ""),
            scan=self.plotms_kw.get("scan", ""),
            correlation=self.plotms_kw.get("correlation", ""),
            array=self.plotms_kw.get("array", ""),
            observation=self.plotms_kw.get("observation", ""),
            intent=self.plotms_kw.get("intent", ""),
            feed=self.plotms_kw.get("feed", ""),
            msselect=self.plotms_kw.get("msselect", ""),
            averagedata=True,
            avgchannel=self.plotms_kw.get("avgchannel", ""),
            avgtime=self.plotms_kw.get("avgtime", ""),
            avgscan=self.plotms_kw.get("avgscan", False),
            avgfield=self.plotms_kw.get("avgfield", False),
            avgbaseline=self.plotms_kw.get("avgbaseline", False),
            avgantenna=self.plotms_kw.get("avgantenna", False),
            avgspw=self.plotms_kw.get("avgspw", False),
            scalar=self.plotms_kw.get("scalar", False),
            transform=True,
            freqframe=self.plotms_kw.get("freqframe", ""),
            restfreq=self.plotms_kw.get("restfreq", ""),
            veldef=self.plotms_kw.get("veldef", "RADIO"),
            #phasecenter=self.plotms_kw.get("phasecenter", ""), # deprecated?
            shift=self.plotms_kw.get("shift", [0.0, 0.0]),
            extendflag=self.plotms_kw.get("extendflag", False),
            extcorr=self.plotms_kw.get("extcorr", False),
            extchannel=self.plotms_kw.get("extchannel", False),
            iteraxis=self.plotms_kw.get("iteraxis", ""),
            plotfile=self.datafile_prefix + ".txt",
            expformat="txt",
            verbose=True,
            exprange="all",  # `all` not work for .txt output... maybe bug? -> set gridrows = maxnplot
            overwrite=False,
            showgui=False,
            showatm=self.plotms_kw.get("showatm", False),
            showtsky=self.plotms_kw.get("showtsky", False),
            showimage=self.plotms_kw.get("showimage", False),
        )
        
        if self.verbose:
            print("Done.")

        # split out each plot data; will be overwritten
        #os.system('rm {:s}'.format(self.datafile_prefix + '_*.txt'))
        subprocess.run('rm {:s}'.format(self.datafile_prefix + '_*.txt'), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        filelist = pycsplit(self.datafile_prefix + ".txt", split_str="# From plot")

        # rename the file
        self.headerfile = filelist[0].replace('_0.txt', '_header.txt')
        #os.system('mv {:s} {:s}'.format(filelist[0], self.headerfile))
        subprocess.run('mv {:s} {:s}'.format(filelist[0], self.headerfile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        self.datafilelist = []
        self.miscfilelist = []
        for i, file in enumerate(filelist[1:]):
            filenames = pycsplit(file, split_str="# x y")
            
            datafile = filenames[1].replace('_{:d}_1.txt'.format(i+1), '_data{:d}.txt'.format(i))
            #os.system('rm {:s}'.format(filenames[0]))
            subprocess.run('rm {:s}'.format(filenames[0]), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #os.system('mv {:s} {:s}'.format(filenames[1], datafile))
            subprocess.run('mv {:s} {:s}'.format(filenames[1], datafile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.datafilelist.append(datafile)
            
            if self.misc:
                assert len(filenames) == 3
                miscfile = filenames[2].replace('_{:d}_2.txt'.format(i+1), '_data{:d}_misc.txt'.format(i))
                #os.system('mv {:s} {:s}'.format(filenames[2], miscfile))
                subprocess.run('mv {:s} {:s}'.format(filenames[2], miscfile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.miscfilelist.append(miscfile)
            
            #os.system('rm {:s}'.format(file))
            subprocess.run('rm {:s}'.format(file), shell=True)
    
    def load_data(self,):
        
        self.iteraxis = self.plotms_kw.get('iteraxis', '')
        self.nplot = len(self.datafilelist)
        
#         dflist = []
#         for datafile in self.datafilelist:
#             df = pd.read_table(datafile, comment='#', sep='\s+', names=plotms_header, index_col=False)
#             dflist.append(df)
        
#         if self.misc:
#             miscdflist = []
#             for miscfile in self.miscfilelist:
#                 df = pd.read_table(datafile, comment='#', sep='\s+', names=plotms_header, index_col=False)
#                 miscdflist.append(df)
        
        self.data = {}
        
        if self.iteraxis == '':
            assert self.nplot == 1
            df = pd.read_table(self.datafilelist[0], comment='#', sep='\s+', usecols=(0,1), names=['x', 'y'])
            self.data[self.iteraxis] = {'data': {self.xaxis: df['x'],
                                  self.yaxis: df['y'], 
                                 },
                        }
            if self.misc:
                df = pd.read_table(self.miscfilelist[0], comment='#', sep='\s+', usecols=(0,1), names=['x', 'y'])
                self.data[self.iteraxis]['misc'] = {'x': df['x'],
                                     'y': df['y'], 
                                    }
        else:
            for i, datafile in enumerate(self.datafilelist):
                names = ['x', 'y', 'ant1name', 'ant2name'] if self.iteraxis == 'baseline' else ['x', 'y', self.iteraxis]
                df = pd.read_table(datafile, comment='#', sep='\s+', usecols=(0,1,*hdr_col[self.iteraxis]), names=names)
                
                if len(df) == 0:
                    #idx = 
                    data = {'data': {self.xaxis: np.array([]),
                                 self.yaxis: np.array([]), 
                                },
                       } 
                
                else:
                    idx = df['ant1name'] + np.array(['-']) + df['ant2name'] if self.iteraxis == 'baseline' else df[self.iteraxis]
                    idx = np.asscalar(np.unique(idx))
                    data = {'data': {self.xaxis: df['x'],
                                     self.yaxis: df['y'], 
                                    },
                           } 
                       
                if self.misc:
                    df = pd.read_table(self.miscfilelist[i], comment='#', sep='\s+', usecols=(0,1), names=['x', 'y'])
                    if len(df) == 0:
                        data['misc'] = {'x': np.array([]),
                                        'y': np.array([]), 
                                        }
                    else:
                        data['misc'] = {'x': df['x'],
                                        'y': df['y'], 
                                        }
                        
                self.data['{:s}_{}'.format(self.iteraxis, idx)] = data
            
                

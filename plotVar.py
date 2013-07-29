import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

import matplotlib as mpl
mpl.rc('font', family='serif')
from matplotlib.backends.backend_pdf import PdfPages

import scipy as sp
import numpy as np
import pylab as pl
import scipy.io as io
import scipy.optimize as opt
import scipy.stats as st
import scipy.special as spec
import os
import string
import sys

pl.rcParams['xtick.labelsize'] = 14
pl.rcParams['ytick.labelsize'] = 14
fontsize1 = 16
fontsize2 = 20

#filename = 'text/APS8inchTests/Apr05/v2/L0/ButterworthFilterMCP800Trigger800/PC2MCP-300_900_300_1100_400@3_15.txt'
filename = '/psec/data1/vostrikov/text/APS8inchTests/Apr05/v2/L0/PC2MCP-300_900_300_1100_400@3_15.txt'

a = np.genfromtxt(filename)

frame = a[:, 0]
A_1, Npulses_1, I_1, m_1, tm_1, ta_1, tb_1, sgn_1, rt_1, ft_1, UV_1, sigmaL_1, sigmaR_1, noiseL_1, noiseR_1 = a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5], a[:, 6], a[:, 7], a[:, 8], a[:, 9], a[:, 10], a[:, 11], a[:, 12], a[:, 13], a[:, 14], a[:, 15]
w10_1, w25_1, w50_1, w75_1, fw10_1, fw25_1, fw50_1, fw75_1, I10_1, I25_1, I50_1, I75_1 = a[:, 16], a[:, 17], a[:, 18], a[:, 19], a[:, 20], a[:, 21], a[:, 22], a[:, 23], a[:, 24], a[:, 25], a[:, 26], a[:, 27]
I2_1 = a[:, 28]
I3_1 = a[:, 29]
t10_1, ft10_1, t25_1, ft25_1, t50_1, ft50_1, t75_1, ft75_1 = a[:, 30], a[:, 31], a[:, 32], a[:, 33], a[:, 34], a[:, 35], a[:, 36], a[:, 37]
Afit_1, Tfit_1, Sfit_1, Xfit_1, ni_1, nf_1 = a[:, 38], a[:, 39], a[:, 40], a[:, 41], a[:, 42], a[:, 43]

pl.figure()
pl.grid(True)
pl.hist2d(ni_1, Tfit_1)
#counts, bins, patches = pl.hist(1e9*t0, bins=99, range=[15.5, 16])
#pl.hist(x, bins=99, color='b', range=[0, 500], lw=3)
#leg = pl.legend(loc='upper right', prop={'size':fontsize1},
#                fancybox=True, shadow=True)
#pl.xlabel(r"Number of electrons, millions", fontsize=fontsize1)
#pl.ylabel(r"Events", fontsize=fontsize1)
#pl.title(r"$\mu_{\mathrm{fit}} = %.1f$ psec, $\sigma_{\mathrm{fit}} = %.1f$ psec" % (p1[1], p1[2]), fontsize=fontsize2)
#pl.xlim(0, 600)
#pl.ylim(0, 400)
#pl.text(x=62, y=560, s=r'%d events' % (len(I)), backgroundcolor='w', fontsize=fontsize2)
#pl.savefig('I50_CH1.eps', pad_inches=0.0)
pl.show()


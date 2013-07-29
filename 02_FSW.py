import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

import matplotlib as mpl
mpl.use('Agg')
mpl.rc('font', family='serif')
#mpl.rc('text', usetex=True)

import scipy as sp
import numpy as np
import pylab as pl
import scipy.io as io
import scipy.optimize as opt
import os
import sys

pl.rcParams['xtick.labelsize'] = 16
pl.rcParams['ytick.labelsize'] = 16
fontsize1 = 18
fontsize2 = 22

data = 'data/'
text = 'text/'
dataset = ''
plot = 'plot/'
triggerID = 1
Nchannels = 4
dt = 100e-12
logBoolean = False

if not os.path.exists(plot + dataset):
	os.makedirs(plot + dataset)

ls = os.listdir(text + dataset)
i = 0
while i < len(ls):
	if os.path.isfile(text + dataset + ls[i]):
		if ls[i].split('.')[-1] == "dat":
			i += 1
		else:
			del ls[i]
	else:
		del ls[i]

for filename in ls:
	print filename
	a = np.genfromtxt(text + dataset + filename)
	channel = a[:, 1]
	Nparameters = int(filename.split('@')[1].split('.')[0].split('_')[1])
	Nfractions = int((len(a[0, :]) - Nparameters - Nchannels - 5) / 5)
	FSW = a[:, Nparameters+Nchannels+3*Nfractions+2] - a[:, Nparameters+Nchannels+3*Nfractions+1]

	for channelID in range(Nchannels):
		if channelID + 1 != triggerID:
			pl.figure()
			pl.grid(True)
			Nbins = int(max(FSW[np.where(channel+1!=triggerID)]) / dt) + 1
			bins = np.linspace(-0.5 * dt * 1e9, (Nbins+1) * dt * 1e9 - 0.5 * dt * 1e9, Nbins+2)
			counts, bins, patches = pl.hist(1e9*FSW[np.where(channel==channelID)],
							bins=bins, log=logBoolean)
			axesRange = list(pl.axis())
			axesRange[0] = 0.0
			axesRange[1] = 6.0
			
			pl.xlabel('Full signal width (nsec)', fontsize=fontsize1)
			pl.ylabel('Events / (%d psec)' % (1e12*dt), fontsize=fontsize1)
			pl.title('Channel #%d: %d events' % (channelID+1, len(FSW[np.where(channel==channelID)])), fontsize=fontsize2)
			
			locs, labels = pl.xticks()
			pl.axis(axesRange)
			pl.savefig(plot + dataset + 'FSW_CH' + str(channelID+1) + '__' + '.'.join(filename.split('.')[:-1]) + '.eps', pad_inches=0.0)
			pl.close()


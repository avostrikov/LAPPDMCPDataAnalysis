import matplotlib as mpl
mpl.use('TkAgg')
mpl.rc('font', family='serif')
mpl.rc('text', usetex=True)

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

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

data = 'data/'
text = 'text/'
plot = 'plot/'
#dataset = 'APS8inchTests/Apr05/v2/L0/'
dataset = 'sim/ButterworthFilterMCP800TriggerInf/'
#dataset = 'APS8inchTests/2012-06-22_timing_optimization/v3/L0/'

map = [0, 1, 2, 3]	# for June data
map = [1, -1, 0]	# for April data
map = [-1, 1]
a0 = [10, 50, 70, 90, 110, 130, 150, 180, 1000]	# for April data
a0 = [10, 50, 90, 130, 180, 1000]		# for June data
a0 = [0, 1000]
a = 300 / np.arange(14.5, 0, -1)
a = (a[1:] + a[:-1]) / 2
a0 = [0]
for i in range(len(a)):
	a0.append(a[i])
a0.append(1000)

Nchannels = len(map)
if 0 in map:
	triggerID = map.index(0)
else:
	triggerID = -1

if not os.path.exists(plot + dataset):
	os.makedirs(plot + dataset)

ls = os.listdir(text + dataset)
i = 0
while i < len(ls):
	if os.path.isfile(text + dataset + ls[i]):
		if ls[i].split('.')[-1] == "txt":
			i += 1
		else:
			del ls[i]
	else:
		del ls[i]

################################################################################
################################################################################
################################################################################
################################################################################

def resolution(h, dt, xlabel, filename, coeff):
	dt /= coeff
	pl.figure()
        pl.grid(True)
        h = h[h==h]
	if len(h) < 10:
		return 0, 0
	try:
		bins = np.arange(min(h)-50*dt, max(h)+50*dt, dt/0.99)
	except:
		if coeff == 1:
                        bins = np.arange(-1000, 1000, dt/0.99)
                elif coeff == 1e3:
                        bins = np.arange(0, 100, dt/0.99)
        if len(bins) > 2000:
                if coeff == 1:
                        bins = np.arange(-1000, 1000, dt/0.99)
                elif coeff == 1e3:
                        bins = np.arange(0, 100, dt/0.99)
			#bins = np.arange(40, 60, dt/0.99)
        counts, bins, patches = pl.hist(h, bins=bins)
	fitfunc = lambda p, x: p[0] * np.exp(-(x-p[1])**2/2/p[2]**2)
        errfunc = lambda p, x, y: fitfunc(p, x) - y
	Range = np.all([h>45, h<55], axis=0)
        p0 = [np.max(counts), np.mean(h[Range]), np.std(h[Range])]
	p0 = [np.max(counts), np.mean(h), np.std(h)]
	nSigma = 2
        x = bins[:-1] + 0.5*np.diff(bins)[0]
        y = counts
        len0 = len(x)
        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
        y = y[np.where(np.abs(x - p1[1]) < nSigma * np.abs(p1[2]))]
        x = x[np.where(np.abs(x - p1[1]) < nSigma * np.abs(p1[2]))]
        len1 = len(x)
        while len1 < len0:
                len0 = len(x)
                try:
                        p1, success = opt.leastsq(errfunc, p1, args=(x, y))
                        y = y[np.where(np.abs(x - p1[1]) < nSigma * np.abs(p1[2]))]
                        x = x[np.where(np.abs(x - p1[1]) < nSigma * np.abs(p1[2]))]
                        len1 = len(x)
                except:
                        break
        tmp = np.linspace(min(bins), max(bins), 100*len(bins))
        pl.plot(tmp, fitfunc(p1, tmp), '-r', lw=3)
	pl.xlim(p1[1]-50*dt, p1[1]+50*dt)
        pl.xlabel(xlabel, fontsize=fontsize1)
        pl.ylabel(r'Events / (%d psec)' % (coeff*np.diff(bins)[0]), fontsize=fontsize1)
        pl.title(r'$\sigma=%.2f$ psec' % (coeff*np.abs(p1[2])), fontsize=fontsize2)
        pl.figtext(0.7, 0.85, r'%d events' % (len(h)), fontsize=fontsize1, backgroundcolor='w')
        pl.savefig(plot + dataset + filename + '.eps', pad_inches=0.0)
        #pl.savefig(plot + dataset + filename + '.png', pad_inches=0.0)
        pl.close()
	return coeff*p1[1], coeff*np.abs(p1[2])

################################################################################
################################################################################
################################################################################
################################################################################

print "Loading data..."
a = []
f = []
maxFrame = 0
for i in range(len(ls)):
	filename = ls[i]
	print filename
	if int(filename.split('@')[1].split('.')[0].split('_')[0]) != Nchannels:
		print "Wrong number of channels."
	a.append(np.genfromtxt(text + dataset + filename))
	maxFrame = max(maxFrame, max(a[i][:, 0]))
	if 'MHz' in filename:
		f.append(int(filename.split('MHz')[0].split('.')[-1]))
	else:
		f.append(0)

Nparameters = int(filename.split('@')[1].split('.')[0].split('_')[1])
Ncolumns = (len(a[0][0, :]) - 1) / Nchannels
Nfractions = (Ncolumns - 2 - Nchannels - Nparameters) / 5

print "Drawing histograms..."
singleChannelData = []
striplineData = []
data = []
for n in np.argsort(f):
        data.append([])
	for channel in range(Nchannels):
		data[-1].append([])
		data[-1][-1].append(np.abs(a[n][:, 1+channel*Ncolumns+Nparameters+5*Nfractions+Nchannels-1]))
		data[-1][-1].append(a[n][:, 1+channel*Ncolumns+13])
		if channel == triggerID:
			data[-1][-1].append(a[n][:, 1+channel*Ncolumns+Nparameters+3*Nfractions+Nchannels+3])
			data[-1][-1].append(a[n][:, 1+channel*Ncolumns+Nparameters+5*Nfractions+Nchannels])
		elif triggerID >= 0:
			data[-1][-1].append(a[n][:, 1+channel*Ncolumns+Nparameters+3*Nfractions+Nchannels+3] - a[n][:, 1+triggerID*Ncolumns+Nparameters+3*Nfractions+Nchannels+3])
                        data[-1][-1].append(a[n][:, 1+channel*Ncolumns+Nparameters+5*Nfractions+Nchannels] - a[n][:, 1+triggerID*Ncolumns+Nparameters+3*Nfractions+Nchannels+3])
		else:
			data[-1][-1].append(a[n][:, 1+channel*Ncolumns+Nparameters+3*Nfractions+Nchannels+3])
                        data[-1][-1].append(a[n][:, 1+channel*Ncolumns+Nparameters+5*Nfractions+Nchannels])

f = np.sort(f)
for n in range(len(f)):
	striplineData.append([])
	singleChannelData.append([])
	if f[n] > 0:
        	s = "__%dMHz" % (f[n])
        else:
	        s = ""
	for channel in range(Nchannels):
		singleChannelData[-1].append([])
		if channel != triggerID:
			if (-map[channel] in map) and (map[channel] > 0):
				striplineData[-1].append([])
			for i in range(len(a0)-1):
				print "Channel#%d: Amplitude from %d to %d mV." % (channel+1, a0[i], a0[i+1])
				A = data[0][channel][0]
				Range = np.all([A >= 1e-3*a0[i], A < 1e-3*a0[i+1]], axis=0)
				x = data[n][channel][-2][Range]
				y = data[n][channel][-1][Range]
				NtoS = data[0][channel][1][Range] / A[Range]
				tCFD, sCFD = resolution(1e9*x, dt=10, xlabel=r'Time, nsec',
								filename='t_CFD50_CH' + str(channel+1) + '__%.2f-%.2f%s' % (a0[i], a0[i+1], s),
								coeff=1e3)
				tFIT, sFIT = resolution(1e9*y, dt=10, xlabel=r'Time, nsec',
								filename='t_FIT_CH' + str(channel+1) + '__%.2f-%.2f%s' % (a0[i], a0[i+1], s),
								coeff=1e3)
				singleChannelData[-1][channel].append([tCFD, sCFD, tFIT, sFIT,
									np.mean(A[Range]), np.std(A[Range]), np.mean(NtoS)])
				if (-map[channel] in map) and (map[channel] > 0):
					channel2 = map.index(-map[channel])
					A = 0.5*(data[0][channel][0] + data[0][channel2][0])
					Range = np.all([A >= 1e-3*a0[i], A < 1e-3*a0[i+1]], axis=0)
					x = data[n][channel][-2][Range] - data[n][channel2][-2][Range]
					y = data[n][channel][-1][Range] - data[n][channel2][-1][Range]
					NtoS = 0.5*(data[0][channel][1][Range] + data[0][channel2][1][Range]) / A[Range]
					dtCFD, dsCFD = resolution(1e12*x, dt=2, xlabel=r'Time difference, psec',
										filename='dt_CFD50_CH' + str(channel+1) + '-CH' + str(channel2+1) + '__%.2f-%.2f%s' % (a0[i], a0[i+1], s),
										coeff=1)
					dtFIT, dsFIT = resolution(1e12*y, dt=2, xlabel=r'Time difference, psec',
										filename='dt_FIT_CH' + str(channel+1) + '-CH' + str(channel2+1) + '__%.2f-%.2f%s' % (a0[i], a0[i+1], s),
										coeff=1)
					x = data[n][channel][-2][Range] + data[n][channel2][-2][Range]
                                        y = data[n][channel][-1][Range] + data[n][channel2][-1][Range]
                                        tCFD, sCFD = resolution(1e9*x, dt=10, xlabel=r'Average time, nsec',
									filename='dt_CFD50_CH' + str(channel+1) + '+CH' + str(channel2+1) + '__%.2f-%.2f%s' % (a0[i], a0[i+1], s),
									coeff=1e3)
                                        tFIT, sFIT = resolution(1e9*y, dt=10, xlabel=r'Average time, nsec',
									filename='dt_FIT_CH' + str(channel+1) + '+CH' + str(channel2+1) + '__%.2f-%.2f%s' % (a0[i], a0[i+1], s),
									coeff=1e3)
					striplineData[-1][-1].append([channel, channel2,
									dtCFD, dsCFD, dtFIT, dsFIT,
									tCFD, sCFD, tFIT, sFIT,
									np.mean(A[Range]), np.std(A[Range]), np.mean(NtoS)])
			singleChannelData[-1][channel] = np.array(singleChannelData[-1][channel])
striplineData = np.array(striplineData)

if len(a0) <= 2:
	sys.exit()

def writeTXT(txt, line):
	for i in range(len(line)):
		txt.write(str(line[i]) + '\t')
	txt.write('\n')

print "Drawing plots..."
c = ['b', 'r', 'g', 'c']
##### Single channel data ####################################################
for channel in range(Nchannels):
	if channel != triggerID:
		##### sigma vs. amplitude (single channel) ###################
		pl.figure()
		txt = open(plot + dataset + 'sigma_vs_amplitude_CH%d.txt' % (channel+1), 'w')
		pl.grid(True)
		for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
			pl.plot(1e3*singleChannelData[n][channel][:, 4], singleChannelData[n][channel][:, 1], '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
			pl.plot(1e3*singleChannelData[n][channel][:, 4], singleChannelData[n][channel][:, 3], '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
			writeTXT(txt, singleChannelData[n][channel][:, 4])
			writeTXT(txt, singleChannelData[n][channel][:, 1])
			writeTXT(txt, singleChannelData[n][channel][:, 4])
			writeTXT(txt, singleChannelData[n][channel][:, 3])
		pl.ylim(0, pl.axis()[3])
		leg = pl.legend(loc='upper right', prop={'size':fontsize1}, fancybox=True, shadow=True)
		pl.xlabel(r'Pulse amplitude, mV', fontsize=fontsize1)
		pl.ylabel(r'Time resolution, psec', fontsize=fontsize1)
		pl.savefig(plot + dataset + 'sigma_vs_amplitude_CH%d.eps' % (channel+1), pad_inches=0.0)
		txt.close()
		pl.close()
		
		##### sigma vs. inverse amplitude (single channel) ###########
		pl.figure()
                pl.grid(True)
		for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
                        fitfunc = lambda p, x: p[0] * x + p[1]
                        errfunc = lambda p, x, y: fitfunc(p, x) - y
			
                        x = 1.0/singleChannelData[n][channel][:, 4]
                        y = singleChannelData[n][channel][:, 1]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.15, 0.28-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n], fontsize=fontsize1, backgroundcolor='w')

                        x = 1.0/singleChannelData[n][channel][:, 4]
                        y = singleChannelData[n][channel][:, 3]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.15, 0.23-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n+1], fontsize=fontsize1, backgroundcolor='w')
                pl.ylim(0, pl.axis()[3])
                leg = pl.legend(loc='lower right', prop={'size':fontsize1}, fancybox=True, shadow=True)
                pl.xlabel(r'Inverse pulse amplitude, V$^{-1}$', fontsize=fontsize1)
                pl.ylabel(r'Time resolution, psec', fontsize=fontsize1)
                pl.savefig(plot + dataset + 'sigma_vs_inverseamplitude_CH%d.eps' % (channel+1), pad_inches=0.0)
                pl.close()

		##### sigma vs. noise-to-signal (single channel) #############
		pl.figure()
		txt = open(plot + dataset + 'sigma_vs_noisetosignal_CH%d.txt' % (channel+1), 'w')
                pl.grid(True)
		for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
                        fitfunc = lambda p, x: p[0] * x + p[1]
                        errfunc = lambda p, x, y: fitfunc(p, x) - y
			
                        x = 100*singleChannelData[n][channel][:, 6]
                        y = singleChannelData[n][channel][:, 1]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
			writeTXT(txt, x)
			writeTXT(txt, y)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.65, 0.28-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n], fontsize=fontsize1, backgroundcolor='w')

                        x = 100*singleChannelData[n][channel][:, 6]
                        y = singleChannelData[n][channel][:, 3]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
			writeTXT(txt, x)
			writeTXT(txt, y)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.65, 0.23-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n+1], fontsize=fontsize1, backgroundcolor='w')
                pl.ylim(0, pl.axis()[3])
                leg = pl.legend(loc='upper left', prop={'size':fontsize1}, fancybox=True, shadow=True)
                pl.xlabel(r'Noise to signal ratio, \%', fontsize=fontsize1)
                pl.ylabel(r'Time resolution, psec', fontsize=fontsize1)
                pl.savefig(plot + dataset + 'sigma_vs_noisetosignal_CH%d.eps' % (channel+1), pad_inches=0.0)
		txt.close()
                pl.close()

if len(striplineData[0]) > 0:
        for stripline in range(len(striplineData[0, :, 0, 0])):
                pl.figure()
		txt = open(plot + dataset + 'sigma_vs_amplitude_CH%d-CH%d.txt' % (striplineData[0, stripline, 0, 0]+1, striplineData[0, stripline, 0, 1]+1), 'w')
                pl.grid(True)
                for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
                        pl.plot(1e3*striplineData[n, stripline, :, 10], striplineData[n, stripline, :, 3], '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
                        pl.plot(1e3*striplineData[n, stripline, :, 10], striplineData[n, stripline, :, 5], '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
			writeTXT(txt, striplineData[n, stripline, :, 10])
			writeTXT(txt, striplineData[n, stripline, :, 3])
			writeTXT(txt, striplineData[n, stripline, :, 10])
			writeTXT(txt, striplineData[n, stripline, :, 5])
                pl.ylim(0, pl.axis()[3])
                leg = pl.legend(loc='upper right', prop={'size':fontsize1}, fancybox=True, shadow=True)
                pl.xlabel(r'Pulse amplitude, mV', fontsize=fontsize1)
                pl.ylabel(r'Differential time resolution, psec', fontsize=fontsize1)
                pl.savefig(plot + dataset + 'sigma_vs_amplitude_CH%d-CH%d.eps' % (striplineData[0, stripline, 0, 0]+1, striplineData[0, stripline, 0, 1]+1), pad_inches=0.0)
		txt.close()
                pl.close()

                pl.figure()
                pl.grid(True)
                for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
                        fitfunc = lambda p, x: p[0] * x + p[1]
                        errfunc = lambda p, x, y: fitfunc(p, x) - y
			
                        x = 1.0/striplineData[n, stripline, :, 10]
                        y = striplineData[n, stripline, :, 3]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.65, 0.28-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n], fontsize=fontsize1, backgroundcolor='w')

                        x = 1.0/striplineData[n, stripline, :, 10]
                        y = striplineData[n, stripline, :, 5]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.65, 0.23-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n+1], fontsize=fontsize1, backgroundcolor='w')
                pl.ylim(0, pl.axis()[3])
                leg = pl.legend(loc='upper left', prop={'size':fontsize1}, fancybox=True, shadow=True)
                pl.xlabel(r'Inverse pulse amplitude, V$^{-1}$', fontsize=fontsize1)
                pl.ylabel(r'Differential time resolution, psec', fontsize=fontsize1)
                pl.savefig(plot + dataset + 'sigma_vs_inverseamplitude_CH%d-CH%d.eps' % (striplineData[0, stripline, 0, 0]+1, striplineData[0, stripline, 0, 1]+1), pad_inches=0.0)
                pl.close()

                pl.figure()
		txt = open(plot + dataset + 'sigma_vs_noisetosignal_CH%d-CH%d.txt' % (striplineData[0, stripline, 0, 0]+1, striplineData[0, stripline, 0, 1]+1), 'w')
                pl.grid(True)
                for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
                        fitfunc = lambda p, x: p[0] * x + p[1]
                        errfunc = lambda p, x, y: fitfunc(p, x) - y
			
                        x = 100*striplineData[n, stripline, :, 12]
                        y = striplineData[n, stripline, :, 3]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
			writeTXT(txt, x)
			writeTXT(txt, y)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.65, 0.28-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n], fontsize=fontsize1, backgroundcolor='w')

                        x = 100*striplineData[n, stripline, :, 12]
                        y = striplineData[n, stripline, :, 5]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
			writeTXT(txt, x)
			writeTXT(txt, y)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.65, 0.23-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n+1], fontsize=fontsize1, backgroundcolor='w')
                pl.ylim(0, pl.axis()[3])
                leg = pl.legend(loc='upper left', prop={'size':fontsize1}, fancybox=True, shadow=True)
                pl.xlabel(r'Noise to signal ratio, \%', fontsize=fontsize1)
                pl.ylabel(r'Differential time resolution, psec', fontsize=fontsize1)
                pl.savefig(plot + dataset + 'sigma_vs_noisetosignal_CH%d-CH%d.eps' % (striplineData[0, stripline, 0, 0]+1, striplineData[0, stripline, 0, 1]+1), pad_inches=0.0)
		txt.close()
                pl.close()

                pl.figure()
                pl.grid(True)
                for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
                        pl.plot(1e3*striplineData[n, stripline, :, 10], striplineData[n, stripline, :, 7], '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
                        pl.plot(1e3*striplineData[n, stripline, :, 10], striplineData[n, stripline, :, 9], '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
                pl.ylim(0, pl.axis()[3])
                leg = pl.legend(loc='upper right', prop={'size':fontsize1}, fancybox=True, shadow=True)
                pl.xlabel(r'Pulse amplitude, mV', fontsize=fontsize1)
                pl.ylabel(r'Average time resolution, psec', fontsize=fontsize1)
                pl.savefig(plot + dataset + 'sigma_vs_amplitude_CH%d+CH%d.eps' % (striplineData[0, stripline, 0, 0]+1, striplineData[0, stripline, 0, 1]+1), pad_inches=0.0)
                pl.close()

                pl.figure()
                pl.grid(True)
                for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
                        fitfunc = lambda p, x: p[0] * x + p[1]
                        errfunc = lambda p, x, y: fitfunc(p, x) - y
			
                        x = 1.0/striplineData[n, stripline, :, 10]
                        y = striplineData[n, stripline, :, 7]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.15, 0.28-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n], fontsize=fontsize1, backgroundcolor='w')

                        x = 1.0/striplineData[n, stripline, :, 10]
                        y = striplineData[n, stripline, :, 9]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.15, 0.23-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n+1], fontsize=fontsize1, backgroundcolor='w')
                pl.ylim(0, pl.axis()[3])
                leg = pl.legend(loc='lower right', prop={'size':fontsize1}, fancybox=True, shadow=True)
                pl.xlabel(r'Inverse pulse amplitude, V$^{-1}$', fontsize=fontsize1)
                pl.ylabel(r'Average time resolution, psec', fontsize=fontsize1)
                pl.savefig(plot + dataset + 'sigma_vs_inverseamplitude_CH%d+CH%d.eps' % (striplineData[0, stripline, 0, 0]+1, striplineData[0, stripline, 0, 1]+1), pad_inches=0.0)
                pl.close()

                pl.figure()
                pl.grid(True)
                for n in range(len(f)):
                        if f[n] > 0:
				s = r" \& Filter"
                        else:
                                s = r""
                        fitfunc = lambda p, x: p[0] * x + p[1]
                        errfunc = lambda p, x, y: fitfunc(p, x) - y
			
                        x = 100*striplineData[n, stripline, :, 12]
                        y = striplineData[n, stripline, :, 7]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n], label='CFD' + s)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.15, 0.28-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n], fontsize=fontsize1, backgroundcolor='w')

                        x = 100*striplineData[n, stripline, :, 12]
                        y = striplineData[n, stripline, :, 9]
                        pl.plot(x, y, '.-', ms=10, lw=3, c=c[2*n+1], label='Fit' + s)
                        p0 = [(y[-1] - y[0]) / (x[-1] - x[0]), 0]
                        p1, success = opt.leastsq(errfunc, p0, args=(x, y))
                        pl.figtext(0.15, 0.23-n*0.1,
				   r'$y=%.1f+%.1fx$' % (p1[1], p1[0]),
				   color=c[2*n+1], fontsize=fontsize1, backgroundcolor='w')
                pl.ylim(0, pl.axis()[3])
                leg = pl.legend(loc='lower right', prop={'size':fontsize1}, fancybox=True, shadow=True)
                pl.xlabel(r'Noise to signal ratio, \%', fontsize=fontsize1)
                pl.ylabel(r'Average time resolution, psec', fontsize=fontsize1)
                pl.savefig(plot + dataset + 'sigma_vs_noisetosignal_CH%d+CH%d.eps' % (striplineData[0, stripline, 0, 0]+1, striplineData[0, stripline, 0, 1]+1), pad_inches=0.0)
                pl.close()

#!/usr/bin/env python -W ignore::FutureWarning
##############################################################################
# External libraries import ##################################################
##############################################################################
#-----------------------------------------------------------------------------
# Import the library for plotting
import matplotlib as mpl
mpl.use('TkAgg')
mpl.rc('font', family='serif')	# set font as in LaTeX
mpl.rc('text', usetex=True)	# allow to use LaTeX commands
#-----------------------------------------------------------------------------
# Import the library to suppress warning messages (not necessary)
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
#-----------------------------------------------------------------------------
# Import data analysis libraries
import scipy as sp
import scipy.io as io
import scipy.optimize as opt
import scipy.interpolate as inter
import scipy.stats as st
import numpy as np
import pylab as pl
#-----------------------------------------------------------------------------
# Import system libraries
import sys	# to use system commands
import os	# to access folders contents etc.
import shutil	# to copy files
#-----------------------------------------------------------------------------
# Import libraries to work with time and dates
import time	# to get current time
import datetime	# to convert time to string
#-----------------------------------------------------------------------------
# Import of self-made library to read data collected from MCP
import mcplib as mcp

##############################################################################
# Constants definition #######################################################
##############################################################################
startTime = time.time()	# time when the script was run
#-----------------------------------------------------------------------------
# Paths to the data.
## All the data are supposed to be splitted to datasets. Dataset is a separate
## folder containing data files and not other folders. It may contain folders,
## but they will be ignored. Files with unknown extentions will be also
## ignored. The program creates a folder for internal data to transfer them
## between the scripts. The structure of this folder is similar to the data
## folder, but with data files created by the program. In the following
## scripts the program creates a folder to store the plots.
data_folder = 'data/'			# path to the MCP data folder
text_folder = 'text/'			# path to the program data folder
dataset = 'APS8inchTests/Apr05/v2/L0'	# dataset name
#-----------------------------------------------------------------------------
# Analysis parameters
f0 = 800e6			# cut-off frequency in Hertz
x = [0.10, 0.25, 0.5, 0.75]	# thresholds for timing in fractions
Nplots = 100			# number of waveforms to display
precision = 1e-12		# spline precision in seconds
#-----------------------------------------------------------------------------
# Plotting parameters
fontsize = {'axis_name': 12, 'axis_tick': 10, 'legend': 10,
	    'subtitle': 14, 'title': 16}		# fontsizes
pl.rcParams['xtick.labelsize'] = fontsize['axis_tick']	# set x ticks fontsize
pl.rcParams['ytick.labelsize'] = fontsize['axis_tick']	# set y ticks fontsize

##############################################################################
# Functions definition #######################################################
##############################################################################
## Prints out process time of the script.
def printTime():
        print '[' + str(datetime.timedelta(seconds=time.time()-startTime))[:9] + ']',
#-----------------------------------------------------------------------------
## Input arguments analysis
def inputAnalysis():
	command = {}
	command['interpolation'] = False
	for i in range(1, len(sys.argv)):
		key = sys.argv[i]
		val = 0
		if '=' in sys.argv[i]:
			key = sys.argv[i][:sys.argv[i].find('=')]
			val = sys.argv[i][sys.argv[i].find('=')+1:].split(',')
		if key == 'window':
			for j in range(len(val)):
				val[j] = 1e-9 * float(val[j])
		elif key == 'triggerID':
			val = int(val[0]) - 1
		elif key == 'interpolation':
			val = True
		else:
			print 'ERROR: Unknown command:', key
			print printTime(), 'Termination.'
			exit()
		command[key] = val
	return command
#-----------------------------------------------------------------------------
# Analysis parameters setting.
## Converts 6-numbers time windows input to 2D list
def timeWindow(time, window, triggerID, Nchannels):
	tiT = window[0]		# ti for trigger
	tfT = window[1]		# tf for trigger
	dt = window[2]		# window width for signal channels
	TW = []
	i = 3
	for channel in range(Nchannels):
		# window time is set a half tick bigger than requested on both sides to include more points
		if channel == triggerID:
			TW.append([np.where(time > tiT - 0.5*np.diff(time)[0])[0][0],
				np.where(time < tfT + 0.5*np.diff(time)[0])[0][-1] + 1])
		else:
			TW.append([np.where(time > window[i] - 0.5*np.diff(time)[0])[0][0],
				np.where(time < window[i] + dt + 0.5*np.diff(time)[0])[0][-1] + 1])
			i += 1
	return TW
## Definition of 
def zeroLevelDefinition(V, TW):
	if np.ndim(V) == 2:
		Nchannels, Npoints = np.shape(V)
		m = np.zeros((Nchannels, Npoints))
		for channel in range(Nchannels):
			m[channel] = np.zeros(Npoints) + np.mean(V[channel, :TW[channel][0]-3*(TW[channel][1]-TW[channel][0])])
		return m
	Nchannels, Nframes, Npoints = np.shape(V)
	m = np.zeros((Nchannels, Nframes, Npoints))
	for frame in range(Nframes):
		m[:, frame, :] = zeroLevelDefinition(V[:, frame, :], TW)
	return m

def parameterDefinition(time, data, command):
	def subplotStructure(n):
		h = int(np.sqrt(float(n)))
		d = n/h + (n % h != 0)
		return 100*h + 10*d + 1
	Nchannels, Nframes, Npoints = np.shape(data)

	if ('triggerID' not in command.keys()) or ('window' not in command.keys()):
		pl.figure()
		for channel in range(Nchannels):
			pl.subplot(subplotStructure(Nchannels) + channel)
			pl.grid()
			pl.title(r'Channel \#' + str(channel+1), fontsize=fontsize['subtitle'])
			for frame in range(Nplots):
				pl.plot(1e9*time, 1e3*data[channel, frame, :], color='b')
			pl.xlabel(r'Time (nsec)', fontsize=fontsize['axis_name'])
			pl.ylabel(r'Voltage (mV)', fontsize=fontsize['axis_name'])
		pl.tight_layout()
		pl.show()
		pl.close()
		
		if 'triggerID' not in command.keys():
			command['triggerID'] = int(raw_input('Trigger ID: ')) - 1
		if 'window' not in command.keys():
			command['window'] = 1e-9 * np.array(raw_input('Time window (nsec): ').split(','), dtype=float)

	command['window'] = timeWindow(time, command['window'], command['triggerID'], len(data))
	Mdata = zeroLevelDefinition(data, command['window'])

	return command, data, Mdata

###############################################################################
###############################################################################
###############################################################################
def GetPulse(v, m, sign, ti, tf):
    if max(sign*(v - m)[ti:tf]) > 0:
        i = np.where(max(sign*(v - m)[ti:tf]) == sign*(v - m)[ti:tf])[0][0] + ti
        a = i
        while ((v[a] - m[a]) * (v[i] - m[a]) > 0) and (a > 0):
            a -= 1
        b = i
        while ((v[b] - m[b]) * (v[i] - m[b]) > 0) and (b < len(v) - 1):
            b += 1
        c = b
        while ((v[c] - m[c]) * (v[i] - m[c]) < 0) and (c < len(v) - 1):
            c += 1
        return a, b, c
    else:
        return -1, 0, 0

def RemovePulse(v, m, a, c):
    vv = np.copy(v)
    vv[a:c+1] = m[a:c+1]
    return vv

def GetPulseFromChannel(V, m, sign, TW):
    v = np.copy(V)
    A = []
    B = []
    C = []
    a, b, c = GetPulse(v, m, sign, TW[0], TW[1])
    while a >= 0:
        A.append(a)
        B.append(b)
        C.append(c)
        v = RemovePulse(v, m, a, b)
        a, b, c = GetPulse(v, m, sign, TW[0], TW[1])
    return A, B, C

def GetPulses(V, m, command):
    triggerID = command['triggerID']
    TW = command['window']
    Nchannels, Npoints = np.shape(m)
    a = []
    b = []
    c = []
    for channel in range(Nchannels):
        a.append([])
        b.append([])
        c.append([])
        a[channel], b[channel], c[channel] = GetPulseFromChannel(V[channel, :], m[channel, :], -1 + 2 * (channel == triggerID), TW[channel])
    return a, b, c

###############################################################################
###############################################################################
###############################################################################
def GetSpline(x, y):
    if len(x) >= 4:
        tck = inter.splrep(x, y)
        xx = np.linspace(min(x), max(x), int((max(x) - min(x)) / precision) + 1)
        yy = inter.splev(xx, tck)
        return xx, yy
    return x, y

def GetCFDTiming(V, time, m, a, b):
	Nchannels, Npoints = np.shape(m)
	T = []
	for channel in range(Nchannels):
		Npulses = len(a[channel])
		T.append(np.zeros((Npulses, 2*len(x))))
		for pulse in range(Npulses):
			v = V[channel, a[channel][pulse]:b[channel][pulse]+1]
			v = v - m[channel, a[channel][pulse]:b[channel][pulse]+1]
			t = time[a[channel][pulse]:b[channel][pulse]+1]
			tt, vv = GetSpline(t, v)
			for fraction in range(len(x)):
				T[channel][pulse, 2*fraction] = tt[np.where(np.abs(vv) >= x[fraction] * max(np.abs(vv)))[0][0]]
				T[channel][pulse, 2*fraction+1] = tt[np.where(np.abs(vv) >= x[fraction] * max(np.abs(vv)))[0][-1]]
	return T

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def chi2(x, y, f, p, sigmaL):
	SSR = np.sum((y - f(p, x))**2 / sigmaL**2)
	if len(x) - len(p) - 1 > 0:
		return SSR / (len(x) - len(p) - 1)
	return SSR

def fit(x, y, fitfunc, p0, sigma):
	errfunc = lambda p, x, y: fitfunc(p, x) - y
	p1, success = opt.leastsq(errfunc, p0, args=(x, y))
	x2 = chi2(x, y, fitfunc, p1, sigma)
	return p1, x2

def getGaussFit(V, time, m, a, b, par):
	Nchannels, Npoints = np.shape(m)
	fitfunc = lambda p, x: p[0] * np.exp(- (x-p[1])**2 / (2*p[2]**2))	# v = A * np.exp(-(t - B)**2/(2*C**2))
	fitfunc2 = lambda p, x: p[0] * np.exp(- (x-p[1])**2 / (2*p[2]**2)) + p[3]
        A, B, C, X2, ni, nf = [], [], [], [], [], []
        for channel in range(Nchannels):
                Npulses = len(a[channel])
                A.append(np.zeros(Npulses))
                B.append(np.zeros(Npulses))
		C.append(np.zeros(Npulses))
                X2.append(np.zeros(Npulses))
		ni.append(np.zeros(Npulses))
		nf.append(np.zeros(Npulses))
                for pulse in range(Npulses):
                        v = V[channel, a[channel][pulse]:b[channel][pulse]+1] - m[channel, a[channel][pulse]:b[channel][pulse]+1]
                        t = time[a[channel][pulse]:b[channel][pulse]+1]
			im = np.where(np.abs(v) == max(np.abs(v)))[0][0]
			tt, vv = t[im-2:im+3], v[im-2:im+3]
			if (len(tt) == 5) and (len(t) > 10):
				pPeak, x2Peak = fit(tt, vv, fitfunc2, [vv[2], tt[2], tt[-1]-tt[0], 0], par[channel][pulse][11])
				tm = pPeak[1]
        	                v = v[np.where(t <= tm + 3*np.diff(t)[0])]
				t = t[np.where(t <= tm + 3*np.diff(t)[0])]
				nf[channel][pulse] = a[channel][pulse] + len(t) - 1
				Am = pPeak[0] + pPeak[3]
				t = t[np.where(np.abs(v) >= 0.15*np.abs(Am))]
				v = v[np.where(np.abs(v) >= 0.15*np.abs(Am))]
				ni[channel][pulse] = nf[channel][pulse] - len(t) + 1
                	        if len(t) >= 5:
                        	        p0 = [Am, tm, pPeak[2]]
                                	p1, X2[channel][pulse] = fit(t, v, fitfunc, p0, par[channel][pulse][11])
	                                A[channel][pulse] = p1[0]
        	                        B[channel][pulse] = p1[1]
					C[channel][pulse] = np.abs(p1[2])
					if (channel == 0) and (np.abs(A[channel][pulse]) > 100e-3) and (np.abs(A[channel][pulse]) < 500e-3) and False:
						pl.rcParams['xtick.labelsize'] = 20
						pl.rcParams['ytick.labelsize'] = 20
						pl.figure()
						pl.grid(True)
						pl.plot(1e9*time[a[channel][pulse]-20:b[channel][pulse]+1], 1e3*V[channel, a[channel][pulse]-20:b[channel][pulse]+1] - m[channel, a[channel][pulse]-20:b[channel][pulse]+1],
							'.b', ms=10)
						tmp = np.arange(time[a[channel][pulse]-20], time[b[channel][pulse]+1], 1e-12)
						pl.plot(1e9*tmp, 1e3*fitfunc(p1, tmp), '-r', lw=3)
						tmp = pl.axis()
						#pl.vlines([1e9*(min(t)-0.5*np.diff(t)[0]), 1e9*(max(t)+0.5*np.diff(t)[0])], ymin=tmp[2], ymax=tmp[3], lw=3)
						#pl.axis(tmp)
						pl.xlabel(r'Time, nsec', fontsize=24)
						pl.ylabel(r'Voltage, mV', fontsize=24)
						pl.tight_layout()
						pl.savefig('pulseFit.eps', pad_inches=0.0)
						pl.show()
						pl.close()
				else:
	                                A[channel][pulse], B[channel][pulse], C[channel][pulse], X2[channel][pulse] = float('NaN'), float('NaN'), float('NaN'), float('NaN')
			else:
				A[channel][pulse], B[channel][pulse], C[channel][pulse], X2[channel][pulse], ni[channel][pulse], nf[channel][pulse] = float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')
        return [A, B, C, X2, ni, nf]

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def integral(t, v):
	I = 0.0
	for i in range(len(t)-1):
		I += 0.5 * (v[i+1] + v[i]) * (t[i+1] - t[i])
	return I

def getPulseParameters(V, time, m, a, b, command):
	triggerID = command['triggerID']
	TW = command['window']
	Nchannels, Npoints = np.shape(m)
	P = []
	sigmaL = np.zeros(Nchannels)
	Nparameters = 15

	UV = 0
	if triggerID >= 0:
		UV = max(np.abs(V[triggerID, TW[triggerID][1]:] - m[triggerID, TW[triggerID][1]:]))

	for channel in range(Nchannels):
		Npulses = len(a[channel])
		P.append([])
		for pulse in range(Npulses):
			v = V[channel, a[channel][pulse]:b[channel][pulse]+1] - m[channel, a[channel][pulse]:b[channel][pulse]+1]
			t = time[a[channel][pulse]:b[channel][pulse]+1]
			if command['interpolation']:
                                t, v = GetSpline(t, v)
			sgn = -1 + 2 * (channel == triggerID)
			A = max(sgn*v)
			itm = np.where(sgn*v == A)[0][0]
                        tm = t[itm]
			I = integral(t, v)
			Ix, w, fw = [], [], []
			for fraction in x:
				if len(np.where(sgn*v >= fraction*A)[0]) == 0:
					print channel, pulse, fraction, t[-1] - t[0]
					pl.plot(t, v)
					pl.show()
				t1 = t[np.where(sgn*v >= fraction*A)[0][0]]
	                        t2 = t[np.where(sgn*v >= fraction*A)[0][-1]]
        	                w.append(tm - t1)
                	        fw.append(t2 - t1)
				Ix.append(integral(t[np.where(sgn*v >= fraction*A)[0][0]:np.where(sgn*v >= fraction*A)[0][-1]], v[np.where(sgn*v >= fraction*A)[0][0]:np.where(sgn*v >= fraction*A)[0][-1]]))
			I2 = []
			for channel2 in range(Nchannels):
				if channel2 != channel:
					dt = TW[channel2][0] - TW[channel][0]
					v2 = V[channel2, a[channel][pulse]+dt:b[channel][pulse]+1+dt] - m[channel2, a[channel][pulse]+dt:b[channel][pulse]+1+dt]
					t2 = time[a[channel][pulse]+dt:b[channel][pulse]+1+dt]
					if command['interpolation']:
						t2, v2 = GetSpline(t2, v2)
					I2.append(integral(t2, v2))

			ta = t[0]
			if len(np.where(sgn*v[:itm]<=0)[0]) > 0:
				ta = t[np.where(sgn*v[:itm]<=0)[0][-1]]
			tb = t[-1]
			if len(np.where(sgn*v[itm:]<=0)[0]) > 0:
				tb = t[np.where(sgn*v[itm:]<=0)[0][0] + itm]
			rt = tm - ta
			ft = tb - tm

			dt = TW[channel][1] - TW[channel][0]
			noiseL = max(np.abs(V[channel, TW[channel][0] - 2*dt:TW[channel][0] - dt] - m[channel, TW[channel][0] - 2*dt:TW[channel][0] - dt]))
			noiseR = max(np.abs(V[channel, TW[channel][1] + dt:TW[channel][1] + 2*dt] - m[channel, TW[channel][1] + dt:TW[channel][1] + 2*dt]))
			sigmaL = np.std(V[channel, TW[channel][0] - 2*dt:TW[channel][0] - dt] - m[channel, TW[channel][0] - 2*dt:TW[channel][0] - dt])
			sigmaR = np.std(V[channel, TW[channel][0] + dt:TW[channel][0] + 2*dt] - m[channel, TW[channel][0] + dt:TW[channel][0] + 2*dt])

			P[channel].append([])
			P[channel][pulse].append(A)
			P[channel][pulse].append(Npulses)
			P[channel][pulse].append(I)
			P[channel][pulse].append(m[channel][0])
			P[channel][pulse].append(tm)
			P[channel][pulse].append(ta)
			P[channel][pulse].append(tb)
			P[channel][pulse].append(sgn)
			P[channel][pulse].append(rt)
			P[channel][pulse].append(ft)
			P[channel][pulse].append(UV)
			P[channel][pulse].append(sigmaL)
                        P[channel][pulse].append(sigmaR)
			P[channel][pulse].append(noiseL)
			P[channel][pulse].append(noiseR)
			P[channel][pulse].append(w)
                        P[channel][pulse].append(fw)
			P[channel][pulse].append(Ix)
			P[channel][pulse].append(I2)
	return P, Nparameters

###############################################################################
###############################################################################
###############################################################################
def WriteFile(filename, data, Nparameters):
	Nchannels = len(data[0][0])
	f = open(filename, 'w')
	for frame in range(len(data)):
		p = data[frame][0]
		tCFD = data[frame][1]
		tFit = data[frame][2]
		for channel in range(Nchannels):
			Npulses, NfractionsX2 = np.shape(tCFD[channel])
			Nfractions = NfractionsX2 / 2
			for pulse in range(Npulses):
				f.write('%7d\t%2d\t' % (frame, channel))
				for parameter in range(Nparameters):
					f.write('%1.8e\t' % (p[channel][pulse][parameter]))
				for parameter in range(Nparameters, Nparameters+3):
					for fraction in range(Nfractions):
						f.write('%1.8e\t' % (p[channel][pulse][parameter][fraction]))
				for channel2 in range(Nchannels-1):
					f.write('%1.8e\t' % (p[channel][pulse][-1][channel2]))
				for fraction in range(NfractionsX2):
					f.write('%1.8e\t' % (tCFD[channel][pulse, fraction]))
				for i in range(len(tFit)):
					f.write('%1.8e\t' % (tFit[i][channel][pulse]))
				f.write('\n')
	f.close()
	shutil.move(filename,
		    '.'.join(filename.split('.')[:-1]) +
		    '@' + str(Nchannels) + '_' + str(Nparameters) + '.' + filename.split('.')[-1])

###############################################################################
###############################################################################
###############################################################################
import scipy.fftpack
def noiseFilter(v, t):
	FFT = np.fft.rfft(v)
	f = np.fft.fftfreq(v.size, np.diff(t)[0])
	f = np.sort(-f[f<=0])
	newFFT = np.copy(FFT)
	#newFFT[np.any([np.abs(f) >= f0], axis=0)] = 0
	G = lambda f, n: np.sqrt(1.0 / (1.0 + (f/f0)**(2*n)))
	newFFT *= G(f, 1)
	#V = scipy.fftpack.irfft(newFFT)
	V = np.fft.irfft(newFFT)
	return V

###############################################################################
###############################################################################
###############################################################################
def CreateDatFile(filename, command):
	global f0
	print 80*'='
	print 'Creating text files:\n\t' + '.'.join(filename.split('.')[:-1]) + '.%dMHz.dat' % (f0/1e6)

	print printTime(), '1. Load data.'
	time, data = mcp.readFileWin(data_folder + filename)
	print printTime(), '   Complete.'

	print printTime(), '2. Reshape data.'
	command, data, M = parameterDefinition(time, data, command)
	Nchannels, Nframes, Npoints = np.shape(data)
	print printTime(), '   Complete.'

	print printTime(), '3. Process data.'
	res = []
	filteredRes = []
	Nframes = 100
	for frame in range(Nframes):
		V = data[:, frame, :]
		m = M[:, frame, :]
		filteredV = np.zeros(np.shape(V))

		for channel in range(Nchannels):
			if (channel == command['triggerID']):
				filteredV[channel, :] = V[channel, :]
				# this is the trigger cut, if I want to filter the trigger
				f0 = 1200e6
				#comment this out if you don't want to filter the trigger
				filteredV[channel, :] = noiseFilter(V[channel, :], time)
				#this is the cut for filtering everything else
				f0 = 800e6
			else:
				filteredV[channel, :] = V[channel, :]
				filteredV[channel, :] = noiseFilter(V[channel, :], time)
		
		a, b, c = GetPulses(V, m, command)
		par, Nparameters = getPulseParameters(V, time, m, a, b, command)
		tCFD = GetCFDTiming(V, time, m, a, b)
		tFit = getGaussFit(V, time, m, a, b, par)
		res.append([par, tCFD, tFit])

		a, b, c = GetPulses(filteredV, m, command)
                par_f, Nparameters = getPulseParameters(filteredV, time, m, a, b, command)
                tCFD_f = GetCFDTiming(filteredV, time, m, a, b)
		tFit_f = getGaussFit(filteredV, time, m, a, b, par_f)
		filteredRes.append([par_f, tCFD_f, tFit_f])
		#res.append([par, tCFD, tFit, par_f, tCFD_f, tFit_f])

		if ((frame + 1) % 500 == 0) or (frame + 1 == Nframes):
			print printTime(), '  ', frame + 1, 'out of', Nframes, 'events processed.'
	print printTime(), '   Complete.'
	
	print printTime(), '4. Write text files.'
	WriteFile(text_folder + '.'.join(filename.split('.')[:-1]) + '.dat', res, Nparameters)
	WriteFile(text_folder + '.'.join(filename.split('.')[:-1]) + '.%dMHz.dat' % (f0/1e6), filteredRes, Nparameters)
	print printTime(), '   Complete.'

	print 'File ' + '.'.join(filename.split('.')[:-1]) + '.%dMHz.dat' % (f0/1e6), 'are created.'
	print 80*'='

def ProcessFolder(dirname, command):
    dirname += '/'
    ls = os.listdir(data_folder + dirname)
    if not os.path.exists(text_folder + dirname):
        os.makedirs(text_folder + dirname)

    for filename in ls:
        if filename[-4:] in ['.txt', '.mat', '.npz']:
            CreateDatFile(dirname + filename, command)
    
    print 'Folder', dirname[:-1], 'is processed.'

command = inputAnalysis()
ProcessFolder(dataset, command)


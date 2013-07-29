#!/usr/bin/env python -W ignore::FutureWarning
import mcplib

import matplotlib as mpl
mpl.use('TkAgg')
mpl.rc('font', family='serif')
mpl.rc('text', usetex=True)

import scipy as sp
import numpy as np
import pylab as pl
import scipy.io as io
import scipy.optimize as opt
import scipy.stats as st
import scipy.special as spec
import os
import string

pl.rcParams['xtick.labelsize'] = 16
pl.rcParams['ytick.labelsize'] = 16
fs1 = 20
fs2 = 24

def drawSignal(t, v, p, x2, ni, nf, xi, xf, filename):
        fitfunc = lambda p, x: p[0] * np.exp(- (x-p[1])**2 / (2*p[2]**2)) + p[3]
        pl.figure()
        pl.grid(True)
        tmp = np.arange(t[0], t[-1], 1e-12)
        pl.plot(1e9*tmp, 1e3*fitfunc(p, tmp), '-g', lw=3)
        tmp = np.arange(t[ni], t[nf+1], 1e-12)
        pl.plot(1e9*tmp, 1e3*fitfunc(p, tmp), '-r', lw=3)
        pl.plot(1e9*t, 1e3*v, '.b', lw=3)
        pl.xlim(xi, xf)
        pl.xlabel(r'Time, nsec', fontsize=fs1)
        pl.ylabel(r'Voltage, mV', fontsize=fs1)
        pl.figtext(x=0.2, y=0.25, s=r"$\chi^2=%.2f$" % (x2), fontsize=fs1, backgroundcolor='w')
        pl.tight_layout()
        pl.savefig(filename, pad_inches=0.0)
        pl.close()

def drawTrigger(t, v, xi, xf, filename):
        pl.figure()
        pl.grid(True)
        pl.plot(1e9*t, 1e3*v, '.b', lw=3)
        pl.xlim(xi, xf)
        pl.xlabel(r'Time, nsec', fontsize=fs1)
        pl.ylabel(r'Voltage, mV', fontsize=fs1)
        pl.tight_layout()
        pl.savefig(filename, pad_inches=0.0)
        pl.close()

filename = 'data/transverse_0mm.mat'
t, v = mcplib.readFileWin(filename)
filename = 'text/transverse_0mm@4_15.txt'
a = np.genfromtxt(filename)

Frame = a[:, 0]
channel = 1
m = a[:, 1+channel*44+3]
A = a[:, 1+channel*44+38]
mu = a[:, 1+channel*44+39]
sigma = a[:, 1+channel*44+40]
X2 = a[:, 1+channel*44+41]
Ni = a[:, 1+channel*44+42]
Nf = a[:, 1+channel*44+43]

#Range = np.all([Ni==Ni, Nf==Nf, mu==mu], axis=0)
#pl.hist(Ni[Range], bins=99)
#pl.show()

"""
for nf in range(444, 452):
        Range = np.all([Ni==433, Nf==nf, mu==mu], axis=0)
        frame = Frame[Range]
        for i in range(10):
                drawSignal(t, v[channel, frame[i], :], [A[frame[i]], mu[frame[i]], sigma[frame[i]], m[frame[i]]],
                           ni=Ni[frame[i]], nf=Nf[frame[i]], xi=10, xf=20,
                           filename="figs/ch%d_Ni=433_Nf=%d_%s.eps" % (channel, Nf[frame[i]], str(int(frame[i])).zfill(4)))
                drawTrigger(t, v[0, frame[i], :], xi=-2, xf=8,
                            filename="figs/ch0_Ni=433_Nf=%d_%s.eps" % (Nf[frame[i]], str(int(frame[i])).zfill(4)))
"""

for frame in range(10):
        drawSignal(t, v[channel, frame, :], [A[frame], mu[frame], sigma[frame], m[frame]], X2[frame],
                        ni=Ni[frame], nf=Nf[frame], xi=10, xf=20,
                        filename="figs/ch%d_%s.eps" % (channel, str(int(frame)).zfill(4)))
        drawTrigger(t, v[0, frame, :], xi=-2, xf=8,
                        filename="figs/ch0_%s.eps" % (str(int(frame)).zfill(4)))

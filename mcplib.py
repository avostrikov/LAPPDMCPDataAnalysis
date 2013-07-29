#!/usr/bin/env python -W ignore::FutureWarning

import scipy as sp
import scipy.io as io
import numpy as np

################################################################################
################################################################################
################################################################################
def readFile(filename):
	if filename[-4:] == '.mat':
		mat = io.matlab.loadmat(filename)

		Nchannels = len(mat['CH'][0])
		if len(np.shape(mat['CH'][0][0].data)) == 3:
			Nframes = len(mat['CH'][0][0].data)
			Npoints = len(mat['CH'][0][0].data[0])
			Nruns = len(mat['CH'][0][0].data[0][0])

			time = np.array(mat['CH'][0][0].time)[:, 0, 0]
                        data = []
                        for channel in range(Nchannels):
                                data.append(mat['CH'][0][channel].data)
                        Nframes *= Nruns
                        data = sp.swapaxes(np.array(data), 2, 3)
                        data = np.reshape(data, (Nchannels, Nframes, Npoints))
                elif len(np.shape(mat['CH'][0][0].data)) == 2:
                        Nframes = len(mat['CH'][0][0].data)
                        Npoints = len(mat['CH'][0][0].data[0])

                        time = np.array(mat['CH'][0][0].time)[:, 0]
                        data = []
                        for channel in range(Nchannels):
                                data.append(mat['CH'][0][channel].data)
                        data = np.array(data)
        elif filename[-4:] == '.txt':
                a = open(filename, 'r').readlines()
                a = ''.join(a).split('\n\n')[:-1]
                b = []
                for s in a:
                        b.append(s.split('\n'))
                a = []
                for i in range(len(b)):
                        a.append([])
                        for j in range(len(b[i])):
                                a[i].append(b[i][j].split('\t')[:-1])
                a = sp.swapaxes(np.array(a, dtype=float), 0, 1)
                time = a[0, 0, :]
                data = a[1:, :, :]
        return time, data

def readFileWin(filename):
	if filename[-4:] == '.mat':
		mat = io.matlab.loadmat(filename)
		
		Nchannels = len(mat['CH'][0])
		
		if len(np.shape(mat['CH'][0][0][0])) == 3:
			Nframes, Npoints, Nruns = np.shape(mat['CH'][0][0][0])
			
			time = mat['CH'][0][0][1][:, 0, 0]
                        data = []
                        for channel in range(Nchannels):
                                data.append(mat['CH'][0][channel][0])
                        Nframes *= Nruns
                        data = sp.swapaxes(np.array(data), 2, 3)
                        data = np.reshape(data, (Nchannels, Nframes, Npoints))
                elif len(np.shape(mat['CH'][0][0][0])) == 2:
			Nframes, Npoints = np.shape(mat['CH'][0][0][0])
                        
                        time = mat['CH'][0][0][1][:, 0]
                        data = []
                        for channel in range(Nchannels):
                                data.append(mat['CH'][0][channel][0])
                        data = np.array(data)
        elif filename[-4:] == '.txt':
                a = open(filename, 'r').readlines()
                a = ''.join(a).split('\n\n')[:-1]
                b = []
                for s in a:
                        b.append(s.split('\n'))
                a = []
                for i in range(len(b)):
                        a.append([])
                        for j in range(len(b[i])):
                                a[i].append(b[i][j].split('\t')[:-1])
                a = sp.swapaxes(np.array(a, dtype=float), 0, 1)
                time = a[0, 0, :]
                data = a[1:, :, :]
	elif filename[-4:] == '.npz':
		npz = np.load(filename)
		time = npz['time']
		data = npz['data']
        return time, data


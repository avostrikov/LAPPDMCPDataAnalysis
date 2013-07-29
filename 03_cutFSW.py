import numpy as np
import os
import sys

data = 'data/'
text = 'text/'
dataset = ''
plot = 'plot/'
triggerID = 1
Nchannels = 4
cutFSW = 2.9e-9
ti = [0, 13, 20, 17]
tf = [2, 16, 23, 20]
#ti = [0, 10.5, 16, 16]
#tf = [1, 12, 17.5, 17.5]
#ti = [49.25, 49.25]
#tf = [50.75, 50.75]

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

for i in range(len(ls)):
	filename = ls[i]
	print filename

	a = np.genfromtxt(text + dataset + filename)
	len0 = len(a)
	# Cut on full signal width
	Nparameters = int(filename.split('@')[1].split('.')[0].split('_')[1])
	Nfractions = int((len(a[0, :]) - Nparameters - Nchannels - 5) / 5)
	Nfractions = 4
	FSW = a[:, Nparameters+Nchannels+3*Nfractions+2] - a[:, Nparameters+Nchannels+3*Nfractions+1]
	a = a[np.any([FSW>=cutFSW, a[:, 1]==triggerID-1], axis=0), :]
	len1 = len(a)
	# Maximum integral only
	j = 0
	while j < len(a):
		if (a[j, 6] >= 1e-9*ti[int(a[j, 1])]) and (a[j, 6] <= 1e-9*tf[int(a[j, 1])]):
			b = np.where(np.all([a[:j, 0]==a[j, 0], a[:j, 1]==a[j, 1]], axis=0))[0]
			if len(b) == 0:
				j += 1
			elif len(b) == 1:
				k = b[0]
				if a[j, 9] * a[j, 4] > a[k, 9] * a[k, 4]:
					a = np.delete(a, k, axis=0)
				else:
					a = np.delete(a, j, axis=0)
			else:
				print 'Error.'
				exit()
		else:
			a = np.delete(a, j, axis=0)
		if ((j+1) % 5000 == 0) or (j+1 == len(a)):
			print '\tMax integral:', (j+1), 'out of', len(a)
	len2 = len(a)
	# Integral analysis
	for frame in np.unique(a[:, 0]):
		if ((frame+1) % 1000 == 0) or (frame == max(np.unique(a[:, 0]))):
			print frame+1, "out of", max(np.unique(a[:, 0]))+1
		b = np.where(a[:, 0]==frame)[0]
		for j in b:
			for channel in range(Nchannels):
				if channel != a[j, 1]:
					if channel in a[b, 1]:
						k = np.where(np.all([a[:, 0]==frame, a[:, 1]==channel], axis=0))[0][0]
						a[j, Nparameters+3*Nfractions+2+channel-(channel>a[j, 1])] = a[k, 4]

	coeff = 2.0 / 50.0 / 1.6e-19
	a[:, 4] *= coeff
	for channel in range(Nchannels-1):
		a[:, Nparameters+3*Nfractions+2+channel] *= coeff

	A = []
	for frame in range(int(max(a[:, 0]))+1):
		if (frame+1) % 1000 == 0:
			print (frame+1), "out of", int(max(a[:, 0]))+1
		b = a[np.all([a[:, 0] == frame], axis=0), :]
		if len(b) > 0:
			c = np.array([[frame]])
			for channel in range(Nchannels):
				d = b[np.all([b[:, 1] == channel], axis=0), 2:]
				if len(d) == 0:
					d = np.zeros((1, len(a[0, :])-2)) * float('NaN')
				c = np.hstack((c, d))
			if A == []:
				A = np.copy(c)
			else:
				A = np.vstack((A, c))
		else:
			c = np.array([[frame]])
			d = np.zeros((1, Nchannels*(len(a[0, :])-2))) * float('NaN')
			c = np.hstack((c, d))
			if A == []:
                                A = np.copy(c)
                        else:
                                A = np.vstack((A, c))
	A = np.array(A)

	f = open(text + dataset + '.'.join(filename.split('.')[:-1]) + '.txt', 'w')
	for j in range(len(A)):
		for k in range(len(A[j])):
			f.write('%1.8e\t' % (A[j, k]))
		f.write('\n')
	f.close()

	print '\tCut efficiency: %d, %d, %d' % (len2, len1, len0)


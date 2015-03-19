import scipy.io
import matplotlib.pyplot as plt
import pylab
from pprint import pprint

mat = scipy.io.loadmat('20111024x3_EEG_ATM.mat')
print (mat.keys())
print (mat['EEGdata'].dtype)
print (mat['EEGdata']['data'][0][0].shape)

with open('t1.txt') as f:
    array = [float(line) for line in f]

print (len(array))

NFFT = 1500      	# the length of the windowing segments
noverlap = 500		# the size of non - overlaping
Fs = 500  			# the sampling frequency
shape = mat['EEGdata']['data'][0][0].shape	# shape of data matrix

# extracting single column (single electrode of data)
for x in range(0, shape[1] - 2): #last two samples are EOG and ECG
	print (x)
	column = mat['EEGdata']['data'][0][0][:, x]
	fig = plt.figure(x)
	ax1 = fig.add_subplot(2,1,1)
	ax1.plot(column)
	ax2 = fig.add_subplot(2,1,2)
	ax2.specgram(column, NFFT=NFFT, Fs=Fs, noverlap=noverlap)

# pylab.show()
import numpy as np
import scipy as sp
import scipy.signal as signal
from scipy.io import loadmat
#from matplotlib.mlab import specgram, psd, window_hanning
from sklearn.decomposition import FastICA
from sklearn import preprocessing
import string
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata


# retrieve data
data = loadmat("20111024x3_EEG_ATM.mat", squeeze_me=False)
data = data['EEGdata']
idx = data['R128'].flatten()[0].flatten()
eeg = data['data'][0][0][(idx[0]+1500):]

no_electrodes = eeg.shape[1] - 2 
print "No of electrodes: " + str(no_electrodes)
print "Shape of eeg array: " + str(eeg.shape)

# change the EEG to absolute values (starting from 0)
min_eeg = np.min(eeg)

if min_eeg < 0:
    eeg = eeg - min_eeg
else:
    eeg = eeg + min_eeg


# generate GFP vector
gfp = np.std(eeg, axis=1)
print "GFP shape: " + str(gfp.shape)
'''
fig = plt.figure("GFP")
ax1 = fig.add_subplot(2,1,1)
ax1.plot(gfp)
'''
# find peaks in GFP
peaks = signal.find_peaks_cwt(gfp, np.arange(20,50))
print "Number of peaks: " + str(len(peaks))

# prepare microstates array
microstates = np.zeros((len(peaks), no_electrodes))
j = 0
for i in peaks:
    microstates[j, :] = eeg[i, 0:no_electrodes]
    j = j + 1

print "Shape of microstates: " + str(microstates.shape)

# perform temporal ICA

# ax2 = fig.add_subplot(2,1,2)
# ax2.plot(peaks, np.arange(0, len(peaks)),'bs')
# plt.show()

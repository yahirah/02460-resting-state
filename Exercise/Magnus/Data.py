__author__ = 'Magnus'

import numpy as np
import scipy as sp
from scipy.io import loadmat
from matplotlib.mlab import specgram

data = loadmat("D:/Skole/10 Semester/02460 Advanced Machine Learning/Project/20111024x3_EEG_ATM.mat", squeeze_me=False)
data = data['EEGdata']

idx = data['R128'].flatten()[0].flatten() # fMRI
eeg = data['data'][0][0] # EEG



fs = 500 # sample frequency
start = idx[0]-1500

s,f,t = specgram(eeg[start:,1],NFFT=2048,Fs=500,noverlap=1000)
# A spectrogram is a visual representation of the spectrum of frequencies
# Returns the tuple (spectrum, freqs, t)
# NFFT the length of the windowing segments
# fs the sampling frequency

# power spectrum
ps = np.abs(s)**2



s.shape

eeg[start:,1].shape
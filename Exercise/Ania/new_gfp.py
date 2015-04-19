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

# change the EEG to absolute values (starting from 0)
def open_file(path):
    # retrieve data
    data = loadmat(path, squeeze_me=False)
    data = data['EEGdata']
    idx = data['R128'].flatten()[0].flatten()
    eeg = data['data'][0][0][(idx[0]+1500):]

    no_electrodes = eeg.shape[1] - 2 
    print "No of electrodes: " + str(no_electrodes)
    print "Shape of eeg array: " + str(eeg.shape)
    return eeg, no_electrodes

def absolute_values(data):
    min_eeg = np.min(data)
    print "Minimal eeg: " + str(min_eeg)
    if min_eeg < 0:
        result = data - min_eeg
    else:
        result = data + min_eeg
    return result

# generate GFP vector
def generate_gfp(data):
    gfp = np.std(eeg, axis=1)
    print "GFP shape: " + str(gfp.shape)
    return gfp

# find peaks in GFP
def find_peaks(data):
    peaks = signal.find_peaks_cwt(gfp, np.arange(20,50))
    print "Number of peaks: " + str(len(peaks))
    return peaks

def generate_microstates(peeks, data, no_electrodes):
    microstates = np.zeros((len(peaks), no_electrodes))
    j = 0
    for i in peaks:
        microstates[j, :] = eeg[i, 0:no_electrodes]
        j = j + 1
    print "Shape of microstates: " + str(microstates.shape)
    return microstates


path = "20111024x3_EEG_ATM.mat"
data, no_electrodes = open_file(path)
eeg = absolute_values(data)
gfp = generate_gfp(eeg)
peaks = find_peaks(gfp)
microstates = generate_microstates(peaks, eeg, no_electrodes)
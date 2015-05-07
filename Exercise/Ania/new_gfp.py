import numpy as np
import scipy as sp
import scipy.signal as signal
from scipy.io import loadmat
#from matplotlib.mlab import specgram, psd, window_hanning
from sklearn.decomposition import FastICA
from sklearn import preprocessing
import string
import pickle
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata

# change the EEG to absolute values (starting from 0)
def open_file(path):
    # retrieve data
    data = loadmat(path, squeeze_me=False)
    data = data['EEGdata']
    idx = data['R128'].flatten()[0].flatten()
    eeg = data['filtered_data'][0][0][(idx[0]+1500):]

    no_electrodes = eeg.shape[1] # - 2 
    print "No of electrodes: " + str(no_electrodes)
    print "Shape of eeg array: " + str(eeg.shape)
    return eeg, no_electrodes



def standarize(data):
    for i in range(1, data.shape[1]):
        row = data[:, i]
        mean = np.mean(row)
        std = np.std(row)
        data[:, i] = (row - mean)/std
        # print "Variance: " + str(np.var(data[:, i])) + ", mean: " + str(np.mean(data[:,i]))
    return data;

def downsample(data):
    return data[0:data.shape[0]:2, :]

# generate GFP vector
def generate_gfp(data):
    gfp = np.std(eeg, axis=1)
    print "GFP shape: " + str(gfp.shape)
    return gfp

# find peaks in GFP
def find_peaks(data):
    peaks = signal.find_peaks_cwt(gfp, np.arange(1,15), noise_perc=50)
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

def concatenate(main, append):
    print "Shape of main:" + str(main.shape)
    print "Shape of append" + str(append.shape)

    result = np.append(main, append, axis=0)
    return result

def save_dict(name, obj):
    with open(name,'w') as f:
        for path in paths:
            line = ' '.join(str(e) for e in obj[path])
            f.write(line+'\n') # python will convert \n to os.linesep
    f.close()
gen_path = "E:\Dropbox\EEG_data\\filtered"
paths = ["20110822x4_EEG_ATM_filt", "20110823x1_EEG_ATM_filt"] #, "20110912x4_EEG_ATM_filt",
''' "20110926x3_EEG_ATM_filt", "20111024x3_EEG_ATM_filt", "20111121x5_EEG_ATM_filt",
         "20111205x6_EEG_ATM_filt", "20120130x5_EEG_ATM_filt", "20120213x5_EEG_ATM_filt",
         "20120312x6_EEG_ATM_filt", "20130819x2_EEG_ATM_filt", "20130826x1_EEG_ATM_filt",
         "20130903x1_EEG_ATM_filt", "20130908x1_EEG_ATM_filt", "20130908x2_EEG_ATM_filt",
         "20130908x3_EEG_ATM_filt", "20131114x1_EEG_ATM_filt", "20131114x2_EEG_ATM_filt"]
'''
result = np.zeros((0,30))
res_eeg = np.zeros((0,30))
res_peaks = { path : np.zeros((0,0)) for path in paths}
for path in paths:
    print "*** Path no " + str(paths.index(path) + 1) + " of " + str(len(paths))
    p = gen_path + "\\" + path;
    eeg, no_electrodes = open_file(p)
    eeg = downsample(eeg)
    eeg = standarize(eeg)
    gfp = generate_gfp(eeg)
    peaks = find_peaks(gfp)
    microstates = eeg[peaks, :]
    result = concatenate(result, microstates)
    res_eeg = concatenate(res_eeg, eeg)
    res_peaks[path] =  peaks;
np.savetxt("microstates.txt", result)
np.savetxt("eeg_standarized.txt", res_eeg)
save_dict("peaks.txt", res_peaks)
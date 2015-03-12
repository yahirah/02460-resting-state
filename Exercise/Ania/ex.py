import scipy.io
import matplotlib.pyplot as mpl
import pylab

mat = scipy.io.loadmat('20111024x3_EEG_ATM.mat')
print mat.keys()
print mat['EEGdata'].dtype
print mat['EEGdata']['R128']
print mat['EEGdata']['Fs'][0][0][0][0]
print mat['EEGdata']['data'][0][0].shape

with open('t1.txt') as f:
    array = [float(line) for line in f]

print len(array)
# extracting single column (single electrode of data)
print mat['EEGdata']['data'][0][0][:, 0]

mpl.plot(mat['EEGdata']['data'][0][0][:, 0])
pylab.show()
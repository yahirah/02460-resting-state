import scipy.io
import matplotlib.pyplot as plt
import pylab
import math
from sklearn import linear_model
import numpy as np

mat = scipy.io.loadmat('20111024x3_EEG_ATM.mat')
print mat['EEGdata'].dtype
fmri = mat['EEGdata']['R128'][0][0]
print len(fmri)
print fmri[0]
# print mat['EEGdata']['data'][0][0].shape


with open('t1.txt') as f:
    array = [float(line) for line in f]

print len(array)
y = array[:200]

NFFT = 1500      	# the length of the windowing segments
noverlap = 1000		# the size of non - overlaping
Fs = 500  			# the sampling frequency
shape = mat['EEGdata']['data'][0][0].shape	# shape of data matrix

# array rotation
test = [[1,2,3],[4,5,6]]
result = zip(*test)[::-1]
print result

# creating features array without first 3 seconds - delay window

features_array = np.zeros((600, 200*(shape[1]-2)))
for x in range(0, 5): #last two samples are EOG and ECG
	print (x)
	column = mat['EEGdata']['data'][0][0][fmri[0][0]:, x]
	fig = plt.figure(x)
	ax1 = fig.add_subplot(3,1,1)
	ax1.plot(column)
	ax2 = fig.add_subplot(3,1,2)
	Pxx, freqs, bins, im = ax2.specgram(column, NFFT=NFFT, Fs=Fs, noverlap=noverlap)
	spec_shape = Pxx.shape
	print spec_shape
	for i in xrange(0, spec_shape[1]):
		feature =Pxx[:200,i]
		# print len(feature)
		# rotated_feature = zip(*feature)[::-1]
		# print len(rotated_feature)
		features_array[i,x*200:(x+1)*200] = feature

	# print ("######")
	# print (Pxx.shape)
	# print (Pxx[0].shape)
	ax3 =fig.add_subplot(3,1,3)
	ax3.plot(y)
	print ("######")

print features_array.shape
# print features_array

clf = linear_model.Ridge (alpha = .5)
print clf.fit (features_array[:200, :], y) 
# print clf.coef_
# print clf.intercept_ 
# print clf.score(features_array[:200, :], y)

y_hat = clf.predict(features_array[:200, :])

error = np.sqrt(pow((y_hat - y),2) / (200))
# print error

pylab.show()
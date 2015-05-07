import numpy as np
import scipy as sp
from scipy.io import loadmat
import scipy.signal as signal
from sklearn import preprocessing
import matplotlib.pyplot as plt
from sklearn.decomposition import FastICA
from matplotlib.mlab import griddata

def open_file(path):
    # retrieve data
    data = loadmat(path, squeeze_me=False)
    data = data['EEGdata']
    idx = data['R128'].flatten()[0].flatten()
    eeg = data['filtered_data'][0][0][(idx[0]+1500):]

    no_electrodes = eeg.shape[1] # - 2 
    print "Filtered data loaded"
    print "No of electrodes: " + str(no_electrodes)
    print "Shape of eeg array: " + str(eeg.shape)
    return eeg, no_electrodes

gen_path = "/home/frans/documents/eegdata/"
paths = ["20110822x4_EEG_ATM_filt"]
path = paths[0]
for path in paths:
    print "*** Path no " + str(paths.index(path) + 1) + " of " + str(len(paths))
    p = gen_path + path;
    eeg, no_electrodes = open_file(p)
    eeg = preprocessing.scale(eeg) # mean 0, std 1
    gfp = np.std(eeg, axis=1)
    peaks = signal.find_peaks_cwt(gfp, np.arange(1,15),noise_perc=50)
    microstates = eeg[peaks,:]

#np.savetxt("microstates.txt", microstates)
#np.savetxt("eeg_standarized.txt", eeg)

#plt.clf()
#plt.plot(gfp[0:1000])
#plt.stem(peaks, np.ones(len(peaks))*2, linefmt='r-', markerfmt='bo', basefmt='r-')
#plt.savefig("gfp.png")


chanlocs = loadmat("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat", squeeze_me=True)
coords = np.zeros((30,2))
channames = []

for channel in range(chanlocs['chanlocs'].shape[0]-2):
    coords[channel,0] = chanlocs['chanlocs'][channel][4]
    coords[channel,1] = chanlocs['chanlocs'][channel][5]
    channames.append(chanlocs['chanlocs'][channel][0])

rot = np.matrix([[0,-1],[1,0]])
rot_coords = rot*coords.T
coords = np.asarray(rot_coords.T) # rotate pi/2 counterclockwise

def topomap(A):
    """Plot spatial maps
       Input: A mixing matrix from ICA
       Output: A saved png image of tiled spatial maps"""
    plt.clf()
    xi = np.linspace(-90, 90, 100)
    yi = np.linspace(-90, 90, 200)
    nICs = A.shape[1]
    for i in range(1,nICs+1):
        ax = plt.subplot(6,5,i)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.xlim(-90, 90)
        plt.ylim(-90, 90)
        # grid the data.
        if i <= nICs:
            zi = griddata(coords[:,0], coords[:,1], A[:,i-1], xi, yi, interp='linear')
            # contour the gridded data, plotting dots at the nonuniform data points.
            #ax.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
            ax.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow, vmax=abs(zi).max(), vmin=-abs(zi).max())
            ax.set_xticklabels([])
            ax.set_yticklabels([])
    #        plt.colorbar()  # draw colorbar
            # plot data points.
            # uncomment next two lines for channellabels
            #for c in range(len(channames)):
            #    ax.scatter(coords[c,0], coords[c,1], marker=r"$ {} $".format(channames[c]) , c='b', s=50, zorder=10, edgecolors='none')
        plt.title('IC%d' %(i),size=8)
        ax.set_aspect('equal', 'datalim')

    #plt.title('Microstates?')
    plt.savefig("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/topomap.png",dpi=600)

def tempomap(intensities, name):
    plt.clf()
    for i in range(intensities.shape[0]):
        ax = plt.subplot(intensities.shape[0],1,i+1)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.plot(range(intensities.shape[1]),intensities[i,:],linewidth=0.2)
        plt.title('IC%d' %(i+1),size=8)
    #plt.title('Temporal components')
    plt.savefig("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/%s.png" %(name),dpi=600)

nICs = 30
ica = FastICA(n_components=nICs)
S = ica.fit_transform(microstates) # temporal
A = ica.mixing_ # spatial 
print "Done with ICA"

# topomap(A)

#The same decomposition matrix was applied to the continuous EEG, which yielded an intensity value for each IC at each time point.

intensities = np.dot(A,eeg.T)


norm_time_course = np.argmax(abs(intensities),0)
binaryseq = np.zeros(eeg.shape)
binaryseq[range(len(norm_time_course)),norm_time_course] = 1

plt.clf()
for i in range(4):
    ax = plt.subplot(4,1,i+1)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.plot(range(5000),binaryseq[0:5000,i],linewidth=0.2)
#    plt.title('IC%d' %(i+1),size=8)
#plt.title('Temporal components')
plt.savefig("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/regressors.png",dpi=600)

hrf = loadmat("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/hrf.mat", squeeze_me=True)
#plt.clf()
#plt.plot(hrf['hrf'])
#plt.savefig("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/hrf.png")


regressors = np.apply_along_axis(lambda x: signal.fftconvolve(x,hrf['hrf'],mode='same'), axis=0, arr=binaryseq)
print "Done with convolution"

plt.clf()
for i in range(4):
    ax = plt.subplot(4,1,i+1)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.plot(range(regressors.shape[0]),regressors[:,i],linewidth=0.2)
#    plt.title('IC%d' %(i+1),size=8)
#plt.title('Temporal components')
plt.savefig("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/regressors.png",dpi=600)

np.savetxt("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/regressors.txt", regressors)
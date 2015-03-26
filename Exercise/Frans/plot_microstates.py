import numpy as np
import scipy as sp
from scipy.io import loadmat
#from matplotlib.mlab import specgram, psd, window_hanning
from sklearn.decomposition import FastICA
from sklearn import preprocessing
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata

data = loadmat("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/20111024x3_EEG_ATM.mat", squeeze_me=False)
data = data['EEGdata']

idx = data['R128'].flatten()[0].flatten()
eeg = data['data'][0][0]

chanlocs = loadmat("/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat", squeeze_me=True)
coords = np.zeros((30,2))

for channel in range(chanlocs['chanlocs'].shape[0]-2):
    coords[channel,0] = chanlocs['chanlocs'][channel][4]
    coords[channel,1] = chanlocs['chanlocs'][channel][5]

rot = np.matrix([[0,-1],[1,0]])
rot_coords = rot*coords.T
coords = np.asarray(rot_coords.T) # rotate pi/2 counterclockwise

fs = 500
start = idx[0]+1500

X = eeg[start:,0:30]
X_scaled = preprocessing.scale(X)

ica = FastICA(n_components=13)
S = ica.fit_transform(X_scaled) # temporal
A = ica.mixing_ # spatial 
print "Done with ICA"

# plot "microstates"
plt.clf()
xi = np.linspace(-90, 90, 100)
yi = np.linspace(-90, 90, 200)
for i in range(1,17):
    ax = plt.subplot(4,4,i)
    # grid the data.
    if i <= 13:
        zi = griddata(coords[:,0], coords[:,1], A[:,i-1], xi, yi, interp='linear')
        # contour the gridded data, plotting dots at the nonuniform data points.
        CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
        CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow, vmax=abs(zi).max(), vmin=-abs(zi).max())
        ax.set_xticklabels([])
        ax.set_yticklabels([])
#        plt.colorbar()  # draw colorbar
        # plot data points.
        plt.scatter(coords[:,0], coords[:,1], marker='o', c='b', s=5, zorder=10)
        plt.xlim(-90, 90)
        plt.ylim(-90, 90)
plt.title('Microstates?')
plt.savefig("spatial_white_13.png")

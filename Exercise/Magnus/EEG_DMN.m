clear all
close all
clc

load('20111024x3_EEG_ATM.mat')
load('chanlocs.mat')

eeg = EEGdata.data(:,1:30)'; % from andreas
EEG_data = EEGdata.data';
% dim(EEG) = (Channels,time)



%% plotting first 9 channels
for i = 1:9
    subplot(3,3,i)
    plot(EEG_data(i,:))
end

%% Lars Kai ICAMS
% demo the data spectrogram
 %load 20111024x3_EEG_O2
%load 20111024x3_EEG_CO2
load t1.txt
a=EEGdata.data;
[S,Aest]=icaMS(a',0,1);


%% One can use EEG lab to remove artefacts
% http://sccn.ucsd.edu/wiki/Chapter_01:_Rejecting_Artifacts

%% Balance variance across subjects before concatenation?? page 4
% Data is standarized before use

%% Reception of GFP EEG data is presumed here. All concatenated
close all
clc
M = dlmread('result_normalized.txt', ' '); % result.txt from python (Anna)
EEG_GFP = M(2:end, :)'; % first one is zero ;)
size(EEG_GFP)
% dim(EEG) = (Channels,time) = X

%% Apply ICA to decompose these microstates received above into several sources
% Single group ICA analysis using the EEGLAB toolbox
% http://sccn.ucsd.edu/eeglab/
% http://sccn.ucsd.edu/wiki/Chapter_01:_Rejecting_Artifacts
% http://sccn.ucsd.edu/wiki/Chapter_09:_Decomposing_Data_Using_ICA
% File > Set path to include the EEGLAB path


% EEG_GFP = EEG_data(:, 30001:31000); % this is not correct, just using first 1000. REPLACE WITH EEGGFP



%% Type ">> eeglab" on the Matlab command line.
eeglab

% file - import - Matlab array
% Matlab variable: EEG_GFP
% 500 hz
% 32 channels

%% chanlocs
EEG.chanlocs = chanlocs(:,1:30); % run this line and the channels will appear in EEG lab
% edit -> channel locations % tehn eeglab finds the information.



%% Single group ICA with 30 IC's
% http://sccn.ucsd.edu/wiki/EEGLAB - I. Single Subject Data Processing

% tools -> run ica -> choose algorithm
% We receive "main microstates" these are seen: plots -> component maps -> in 2d

%% (YUAN) The same decomposition matrix was applied to the continuous
% EEG, which yielded an intensity value for each IC at each time
% point
% Yuan: same decomposition matrix (A - the mixing matrix) to continous EEG
% (without GFP)

% WHat we want from this:
% S =inv(A)X
% Rows of the S matrix are the time course of the component activity


% A = EEG.icaweights % From EEGLAB
% inv(A) = EEG.icawinv % From EEGLAB
S = EEG.icawinv * EEG_GFP;
figure(3)
plot(S)

figure(4)
plot(S(2,:)) % when is component 2 active

%% Find A from EEGLAB
% http://sccn.ucsd.edu/wiki/Chapter_09:_Decomposing_Data_Using_ICA

% In the case of ICA decomposition, the independent component filters are
% chosen to produce the maximally temporally independent signals available
% in the channel data. These are, in effect, information sources in the
% data whose mixtures, via volume conduction, have been recorded at the
% scalp channels. The mixing process (for EEG, by volume conduction) is
% passive, linear, and adds no information to the data. On the contrary,
% it mixes and obscures the functionally distinct and independent
% source contributions.

% We are satisfied that Infomax ICA (runica/binica) gives stable decompositions with up to hundreds of channels

% ICA components of EEG data are maximally temporally independent,
% but spatially unconstrained -- and therefore able to find maps representing
% the projection of a partially synchronized domain / island / patch / region of cortex,
% no matter how much it may overlap the projections of other (relatively independent) EEG sources.

%% Studying and removing ICA components
% Tools > Reject data using ICA > Reject components by map

%  Component property figures can also be accessed directly by selecting Plot > Component properties.
% (There is an equivalent menu item for channels, Plot > Channel properties).
% Artifactual components are also relatively easy to identify by visual inspection of
% component time course (menu Plot > Component activations (scroll)

% Component activations (time courses), select Plot > Component activations (scroll)

%% Subtracting ICA components from data
% Tools > Remove components
% Tools > Reject using ICA > Reject components by map

%% page 557 Elements, 591 machine
% X =AS
% dim(X) = [electrodes , time]
% A = loadings (coordinate system) (mixing to go back to X)
% dim(A) = [electrodes , components]
% S = latent (underlying latent variables/factors) (components)
% dim(S) = [components , time]

% activity of e.g. the 2nd component is found by the 2nd row of S and 2nd
% column of A
% A[:,2] * S[2,:] = XC2

% XC2 =  where is component 2 active (in the brain)
% S[2,:] = when is component 2 active

% Columns of the A matrix which are the scalp projection of the components
% Rows of the S matrix which are the time course of the component activity

% Now, if one want to remove component number 2 from the data (for instance if component number 2 proved to be an artifact),
% one can simply subtract the matrix above (XC2) from the original data X.





%% label microstates based on the maximum absolute intensity value 
% should give series of 0 and 1's in for microstates over time

%% Accomodate different temporal scales by applying
% normalized time course of microstates convolved with double gamma
% hemodynamic response function (HRF)

%% Correlation with DMN

%% Regression


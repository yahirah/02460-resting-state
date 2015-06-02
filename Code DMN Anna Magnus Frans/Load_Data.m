%% load eeg data
clear all
close all
clc
load('chanlocs.mat')
load('hrf.mat')

%% Load EEG
nSubjects = 8;
eeg = zeros(30,300000*3,nSubjects);
files = dir('D:/Skole/10 Semester/02460 Advanced Machine Learning/Project/Data/filteredICA_conc/*.mat');
s = 1;
for k = 1:length(files)
    load(strcat('D:/Skole/10 Semester/02460 Advanced Machine Learning/Project/Data/filteredICA_conc/',files(k).name));
    if size(X,2) == 299999
        X = [X,X(:,end)];
    end
    if size(X,2) == 300001
        X = X(:,1:end-1);
    end
    size(X,2)
    if mod(k,3) == 1
        subject = X;
    elseif mod(k,3) == 2
        subject = [subject,X];
    else
        eeg(:,:,s) = [subject,X];
        s = s+1;
    end
end
disp('done loading eeg')

% All eeg for the 8 subjects
eegfold = reshape(eeg,30,[])'; % All eeg data in two dimensions
% dim(eegfold) = 30 X 300000*8

%% load fmri data
Y = zeros(nSubjects*600,60);
files = dir('D:/Skole/10 Semester/02460 Advanced Machine Learning/Project/Data/All condition 20 (inter) subject analysis/*.txt');
mask = [((11-1)*600+1):(18*600)]; % 11-18
for k = 1:length(files)
    y = csvread(strcat('D:/Skole/10 Semester/02460 Advanced Machine Learning/Project/Data/All condition 20 (inter) subject analysis/',files(k).name));
    tal = regexp(files(k).name,['\d+'],'match');
    y = y(mask);
    str2num(tal{1});
    Y(:,str2num(tal{1})) = y;
end
Y37DF = Y(:,[3 7]); % default modes
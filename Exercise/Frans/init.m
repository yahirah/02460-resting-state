addpath(genpath('/home/frans/MATLAB packages/eeglab13_4_4b/'))
addpath(genpath('/home/frans/MATLAB packages/FastICA_25/'))

% pipeline:
% load data for 5 subjects
% make leave-one-subject-out cross-validation to estimate {1. robustness of
% found microstates, 2. length of microstate dominance intervals}
% for each fold
    % find peaks
    % do ICA
    % order by explained variance
% find similarity between found microstates across folds
load('/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat')
%% load eeg data
eeg = zeros(300000,30,5);
load('/home/frans/Dropbox/filteredICA/20110926x3_EEG_ATM_filtICA.mat');
eeg(1:size(X,2),:,1) = X';
load('/home/frans/Dropbox/filteredICA/20111024x3_EEG_ATM_filtICA.mat');
eeg(1:size(X,2),:,2) = X';
load('/home/frans/Dropbox/filteredICA/20111121x5_EEG_ATM_filtICA.mat');
eeg(1:size(X,2),:,3) = X';
load('/home/frans/Dropbox/filteredICA/20120213x5_EEG_ATM_filtICA.mat');
eeg(1:size(X,2),:,4) = X';
load('/home/frans/Dropbox/filteredICA/20120312x6_EEG_ATM_filtICA.mat');
eeg(1:size(X,2),:,5) = X';
% eeg = zeros(300000,30,18);
% files = dir('/home/frans/documents/AdvancedMachineLearningStuff/filtered/*.mat');
% for k = 1:length(files)
%     load(strcat('/home/frans/documents/AdvancedMachineLearningStuff/filtered/',files(k).name));
%     size(EEGdata.filtered_data,1)
%     eeg(1:size(EEGdata.filtered_data,1),:,k) = EEGdata.filtered_data;
% end
%% load fmri data
Y = zeros(3000,60);
files = dir('/home/frans/documents/AdvancedMachineLearningStuff/fmriics/*.txt');
mask = [((4-1)*600+1):(6*600) ((9-1)*600+1):(10*600)]; % 4,5,6,9,10
for k = 1:length(files)
    y = csvread(strcat('/home/frans/documents/AdvancedMachineLearningStuff/fmriics/',files(k).name));
    tal = regexp(files(k).name,['\d+'],'match');
    %y(((11-1)*600+1):(11*600)) = []; % array is now 600 smaller, so subject 14 is now 13
    %y(((13-1)*600+1):(13*600)) = []; 
    y = y(mask);
    % the last 8 subjects had first 3 seconds removed from eeg
    %for i=10:17
    %    y((i*600+1):(i*600+3)) = [];
    %end
    str2num(tal{1});
    Y(:,str2num(tal{1})) = y;
end
%% do ICA on peak topographies
eegfold = [];
for subject=1:size(eeg,3)
    eegfold = cat(1,eegfold,eeg(:,:,subject));
end

gfp = std(eegfold,0,2);
figure()
plot(gfp)
[~, locs] = findpeaks(gfp,'MinPeakDist',50);
    % findpeaks(gfp(1:1500)),'MinPeakWidth',3,'MaxPeakWidth',60);
microstates = eegfold(locs,:)';
[A, W] = fastica(microstates);%,'approach','symm'); % using cubic polynomial

S = W*eegfold'; % apply seperation matrix to original data
variance_of_component = zeros(30,1);
for i=1:30
    back_proj = (A(:,i)*S(i,:))';
    variance_of_component(i) = 1 - var(eegfold(:) - back_proj(:))/var(eegfold(:));
end
[~, order1] = sort(variance_of_component);

figure()
for i=1:30
    subplot(5,6,order1(i))
    topoplot(A(:,i),chanlocs,'electrodes','off','style','map');%,'plotchans',chanlocs(1:30).labels);
end

Asorted = A(:,order1);
%%
nICs = 30;
As = zeros(30,nICs,5);
vars = zeros(nICs,5);
for fold=1:5
    tmp = 1:5;
    tmp(fold) = [];
    eegfold = [];
    for subject=tmp
        eegfold = cat(1,eegfold,eeg(:,:,subject));
    end
    gfp = std(eegfold,0,2);
    [~, locs] = findpeaks(gfp,'MinPeakDist',50);
    microstates = eegfold(locs,:)';
    [icasig, A, W] = fastica(microstates,'numOfIC',nICs); % using cubic polynomial
    S = W*eegfold'; % apply seperation matrix to original data
    variance_of_component = zeros(nICs,1);
    for i=1:nICs
        back_proj = (A(:,i)*S(i,:))';
        variance_of_component(i) = 1 - var(eegfold(:) - back_proj(:))/var(eegfold(:));
    end
    [~, order] = sort(variance_of_component);
    As(:,:,fold) = A(:,order);
    vars(:,fold) = variance_of_component(order);
end
%%
dims = factor(nICs);
if size(dims,2)==3
    nrow = dims(1)*dims(2);
    ncol = dims(3);
elseif size(dims,2)==2
    nrow = dims(1);
    ncol = dims(2);
else
    nrow = 1;
    ncol = dims(1);
end
for j=1:5
    figure()
    for i=1:nICs
        subplot(nrow,ncol,(nICs+1)-i)
        topoplot(As(:,i,j),chanlocs,'electrodes','off','style','map');%,'plotchans',chanlocs(1:30).labels);
        title(num2str(vars(i,j)))
    end
end
%% try kmeans
[idx C] = kmeans(microstates',30,'MaxIter',1000);
for i=1:30
    subplot(5,6,order1(i))
    topoplot(C(i,:),chanlocs);
end
%% calculate variance explained by each component
variance_of_component = zeros(30,1);
for i=1:30
    back_proj = (A(:,i)*S(i,:))';
    variance_of_component(i) = 1 - var(eeg(:) - back_proj(:))/var(eeg(:));
end
[~, order] = sort(variance_of_component);
%%
idx = 1:2000:95000;
figure()
c = 1
for i=idx
    subplot(5,10,c)
    topoplot(microstates(:,i),chanlocs);
    c = c+1;
end
%% plot temporal components
% should be same kind of plot as yuan appendix C
figure()
subplot(4,1,1)
plot(gfp(1:250))
c = 2;
for i=1:3
    subplot(4,1,c)
    plot(abs(S(i,1:250)))
    c = c+1;
end

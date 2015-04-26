%% load data
load('../20111024x3_EEG_ATM.mat');
dmn = csvread('/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/t1.txt');
chanlocs = load('/home/frans/Dropbox/DTU/AdvancedMachineLearning/Project/chanlocs.mat');
y = dmn(1:200);
eeg = EEGdata.data(:,1:30);
eog = EEGdata.data(:,31);
ecg = EEGdata.data(:,32);
idx = EEGdata.R128; % every 1500 corresponds to fmri
eeg_less = eeg(idx(1)+1500:end,:);
start = idx(1)+1500; % number of sample where t= 3s after first fmri
plot(eeg_less(start:end,1))

[N,M] = size(eeg);

%% remove mean
X = bsxfun(@minus,eeg(start:end,:),mean(eeg(start:end,:),1));
plot(X(:,1))
%X = eeg(start:end,:);

%% band pass 1-40 Hz
spectrogram(X(:,1),1500,1000,2048,500,'yaxis')
[s,f,t,p] = spectrogram(X(:,1),1500,1000,2048,500,'yaxis');
% Cut-off frequency
%d = fdesign.lowpass('Fp,Fst,Ap,Ast',40,45,1,40,128);

%figure()
%spectrogram(xh,1500,1000,2048,500,'yaxis')
%% ICA
[S,Aica]=icaMS(X',20,1);
figure()
nplots = 30;
for j=1:3
    figure()
    for i=1:10
        subplot(10,1,i)
        plot(S(j*i,1000:9000))
    end
end
figure()
subplot(2,1,1)
plot(eog(start+1000:start+9000))
subplot(2,1,2)
plot(ecg(start+1000:start+9000))

%% scalp maps
[X2,Y] = meshgrid(-80:80,-80:80);
x = [chanlocs.chanlocs.X];
y = [chanlocs.chanlocs.Y];
rot = [0,-1;1,0];
coords = rot*[x; y];
x = coords(1,:);
y = coords(2,:);
for i=1:10
    figure()
    z = griddata(x, y, Aica(i,:),X2,Y);
    contourf(X2,Y,z,'o')
end
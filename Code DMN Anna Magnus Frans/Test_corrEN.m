%% Peaks and corr with default mode
%% Here we calculate the correlation of the microstates with the default mode
close all
clc
clearvars -except eegfold eeg chanlocs Y Y37DF hrf

EEG = eegfold; % used the name EEG in code

gfp = std(eegfold,0,2); %GFP
steps = 3; % has been completed for steps = 3
max_corr = [0 0 0];
max_corr_all = zeros(steps,steps);
IC = 30; % 30 independent components

% Counters
j = 0;
h=0
w=0
max_gfp = 2.5;

%for d=ceil(linspace(0, 50, 5)) % Min distance of GFP peaks (0ms,25ms,50ms,75ms,100ms) 1 sample = 2 ms
hi = 0;
for h=linspace(min(gfp), max_gfp, steps) % Min height of GFP peaks
    hi = hi +1
    wi = 0;
    for w=ceil(linspace(0, 50, steps)) % Min width of GFP peaks
        j = j +1
        wi = wi +1
        
        
        [~, locs] = findpeaks(gfp, 'MinPeakHeight',h, 'MinPeakWidth',w); % Find the peaks, for this design
        microstates = EEG(locs,:)'; % Pre-ICA-microstates
%         rng(32135484);
        rng(84894624);
        [icasig, A, W] = fastica(microstates,'numOfIC', IC); % using cubic polynomial as default. IC is the number of independent components
        S = W*EEG'; % apply seperation matrix to original data (get the sources)
        
        % calculate power, to order A matrix
        % power explained by each component c
        A2 = sum(A.^2,1);
        S2 = sum(S.^2,2);
        [J T] = size(EEG');
        p = zeros(IC,1);
        for c=1:IC;
            p(c) = (A2(c) * S2(c)) / (J*T);
        end
        power = p / sum(p);
        [~, order1] = sort(power,'descend');
        A = A(:,order1); % save all mixing matrices for use afterwards
        S = S(order1,:); % order S and save such that we can use it for elastic net later
        
        % Convolution of hdrf
        [maxes, idx] = max(abs(S'),[],2);
        sequence = zeros(size(eegfold));
        
        for s=1:size(eegfold,1)
            sequence(s,idx(s)) = 1;
        end
        
        regressorsdummy = cellfun(@(x)(conv(x,hrf,'same')), num2cell(sequence,1), 'UniformOutput', false);
        regressors = cell2mat(regressorsdummy);
        
        
        % CALCULATE CORRELATION WITH Y THE DEFAULT MODE NETWORK
        Xds = downsample(regressors,1500);
        [cormat, pval] = corr(Xds,Y37DF);
        correlation = cormat;
        
        max_corr_all(hi,wi) = max(max(abs(correlation)))
        
        if max_corr_all(hi,wi) > max_corr(1)
            max_corr(:) = [max_corr_all(hi,wi) hi wi];
        end
        
    end
end
disp('DONE!')
max_corr(1)
max_corr_all % contains all max correlations (absolute value)

h=linspace(min(gfp), max_gfp, steps);
w=ceil(linspace(0, 50, steps));
h(max_corr(2))
w(max_corr(3))

% figure
% findpeaks(gfp,'MinPeakDistance',d(max_corr(2)), 'MinPeakHeight',h(max_corr(3)), 'MinPeakWidth',w(max_corr(4))); % Find the peaks, for this design
%% Output
% DONE!
% 
% ans =
% 
%     0.0682
% 
% 
% max_corr_all =
% 
%     0.0484    0.0457    0.0682
%     0.0472    0.0511    0.0399
%     0.0459    0.0611    0.0564
% 
% 
% ans =
% 
%     0.0226
% 
% 
% ans =
% 
%     50


% max_corr(1)=0.0682
% max_corr(2)=1
% max_corr(3)=3
%% Calculate regressors of optimal design settings
% Same as inner loop above
EEG = eegfold;
gfp = std(eegfold,0,2); %GFP

[~, locs] = findpeaks(gfp, 'MinPeakHeight',h(max_corr(2)), 'MinPeakWidth',w(max_corr(3))); % Find the peaks, for this design
microstates = EEG(locs,:)'; % Pre-ICA-microstates
rng(84894624); % to reproduce the results
[icasig, A, W] = fastica(microstates,'numOfIC', IC); % using cubic polynomial as default. IC is the number of independent components
S = W*EEG'; % apply seperation matrix to original data (get the sources)

% calculate power, to order A matrix
% power explained by each component c
A2 = sum(A.^2,1);
S2 = sum(S.^2,2);
[J T] = size(EEG');
p = zeros(IC,1);
for c=1:IC;
    p(c) = (A2(c) * S2(c)) / (J*T);
end
power = p / sum(p);
[~, order1] = sort(power,'descend');
A = A(:,order1); % save all mixing matrices for use afterwards
S = S(order1,:); % order S and save such that we can use it for elastic net later

% Convolution of hdrf
[maxes, idx] = max(abs(S'),[],2);
sequence = zeros(size(eegfold));

for s=1:size(eegfold,1)
    sequence(s,idx(s)) = 1;
end

regressorsdummy = cellfun(@(x)(conv(x,hrf,'same')), num2cell(sequence,1), 'UniformOutput', false);
regressors = cell2mat(regressorsdummy);

% CALCULATE CORRELATION WITH Y THE DEFAULT MODE NETWORK
Xds = downsample(regressors,1500);
[cormat, pval] = corr(Xds,Y37DF);
correlation = cormat;
max(abs(correlation))

%% correlations between microstate time components and Ys
figure
[corrmatatrix, pval] = corr(Xds,Y);
corrmatatrix(pval>0.05) = 0;
%corrmatatrix(pval>0.001) = 0;
figure
pcolor(corrmatatrix);
colormap('jet');colorbar;
xlabel('fMRI Independent Components','Fontsize',40), ylabel('Microstates','Fontsize',40)
title('Correlation Matrix','Fontsize',40)
set(gca,'fontsize',20)

%% Elastic Net
% See which Y has the best correlation
max(abs(correlation))

% set which one to use
yall = Y37DF(:,2); % y3 or y7


%% Downsample X (S') in order to fit with y
clc

Xall = downsample(regressors,1500);
size(Xall)

% Using all
X = Xall;
y = yall;

max(abs(corr(X,y)))

[n,p]=size(X); % Extract dimensions
h = 100; % For values of lambda
lambdas = [0 logspace(-4,2,h)]; % Define h values of lambda on a log scale (include 0 for entire lasso path until least squares)
figure
plot(lambdas,'.')

%% Elastic Net - start of computations
close all
clc

% K is number for variables included (for each run)
K = 1:p;
CV = 5;
rng(84894624); % reproduce results
Idx = crossvalind('Kfold',n,CV); % kfold CV with n observations and CV folds.

% Explanation:
% For each fold, for each lambda:
% Find entire EN path (one by one LARS algorithm)
% Beta contains entire path for all (30) variables
% we loop this to find the respectgively error estimate (for that fold,
% lambda and number of parameters)

% The idea is that we ask elesticnet to find path with e.g. 30 variables
% available but maybe only 23 are active (not equal 0).

for i = 1:CV % For each  k-fold
    disp([num2str(i) ' of ' num2str(CV)]);
    
    % Training and validation set
    ytr = y(Idx~=i); ytst = y(Idx==i);
    Xtr = X(Idx~=i,:); Xtst = X(Idx==i,:);
    
    % Center y, normalize X
    [ytr,my] = center(ytr); % center response
    ytst = ytst-my; % center test response
    [Xtr,mx,varx] = normalize(Xtr); % normalize
    Xtst=normalizetest(Xtst,mx,varx); % normalize test
    
    % save ytst for that fold
    ytst_fold(1:size(ytst,1),i) = ytst;
    
    for k=1:length(lambdas) % for each lambda
        
        % Same Beta format as with lasso but we have a new beta for each lambda.
        % Order: [variables, # of variables included]
        [Beta info] = elasticnet(Xtr,ytr,lambdas(k),-K(end)); % run until we have -K(end) variables included. size= p x k(end)
        
        
        for j=1:length(K); % for each number for parameters included in the EN path
            beta_EL = Beta(:,K(j)); % extract the betas
            yhatTr = Xtr*beta_EL; % compute the estimated training response
            yhatTst = Xtst*beta_EL; % compute the estimated test response
            
            % Store training and test error
            % Order: [# of variables, lambda, fold]
            Err_trEN(j,k,i) = ((yhatTr-ytr)'*(yhatTr-ytr))/length(ytr); % training error (MSE)
            Err_tstEN(j,k,i) = ((yhatTst-ytst)'*(yhatTst-ytst))/length(ytst); % test error (MSE)
            K_EN(j,k,i)=sum(beta_EL~=0); % how many ACTIVE variables are there at the jth iteration
        end
    end
end
% Average over folds
% order: [# of variables, lambda, fold]
err_trEN = mean(Err_trEN,3);
err_tstEN = mean(Err_tstEN,3);
K_nonzEN = mean(K_EN,3); % K is number of variables included

% Best_no_of_para has to be multiplied by the step for k. i.e. if the step length for k is
% not one multiply with the step length. (only if not: k=1:30)
[Best_no_of_para lambda_no_EN] = find(err_tstEN == min(min(err_tstEN)))
Best_lambda_EN = lambdas(lambda_no_EN)
Err_EN = err_tstEN(Best_no_of_para,lambda_no_EN)

% Make model with best # of parameters
[yc,my] = center(y); % center response (do not have a bias in model)
[Xc,mx,varx] = normalize(X); % normalize
Beta = elasticnet(Xc,yc,Best_lambda_EN,-Best_no_of_para);
beta_EN = Beta(:,end); % get the last entry since we have used no_of_para above
Err_EN_test = (y-X*beta_EN)'*(y-X*beta_EN)/length(y-X*beta_EN)

figure(1)
subplot(1,3,1)
imagesc(K, lambdas, log(err_trEN)), colorbar, hold on
title('log training error')
xlabel('k (parameters)'), ylabel('\lambda')
% NOTE: Check that training error and test error is not too far apart at
%       "optimal model" (avoid overfitting in to large a degree)

subplot(1,3,2)
imagesc(K, lambdas, log(err_tstEN)), colorbar
title('log test error')
xlabel('k (parameters)'), ylabel('\lambda')
% NOTE: find minimum (and with a model as simple as possible) test error

subplot(1,3,3)
imagesc(K, lambdas, K_nonzEN), colorbar
title('mean no. of variables in solution')
xlabel('k (parameters)'), ylabel('\lambda')

figure(2)
plot(err_tstEN) % both lambda and no of parameters
xlabel('k (parameters)','Fontsize',30), ylabel('Validation Error (for each lambda)','Fontsize',30)
title('Number of Parameters','Fontsize',40)
set(gca,'fontsize',20)

figure(3)
plot(y,'r')
hold on
plot(X*beta_EN,'b')
xlabel('Time (sample points)','Fontsize',40), ylabel('Activation','Fontsize',40)
title('Prediction','Fontsize',40)
set(gca,'fontsize',20)

% ytst_fold contains the y (test) for each fold that was used to find the
% error
disp(' ')
disp('########## Under 1? #########')
Err_EN / mean(var(ytst_fold))


% Plots the microstaes that were selected
figure(4)
para = find(Beta(:,end));
for i=1:size(para,1)
    subplot(4,4,i)
    topoplot(A(:,para(i)),chanlocs); % first describes most
    title([num2str(para(i))])
    set(gca,'fontsize',20)
end

figure(5)
plot(err_tstEN / mean(var(ytst_fold))) % both lambda and no of parameters
xlabel('Parameters','Fontsize',30), ylabel('Validation Error (for each lambda)','Fontsize',30)
title('Number of Parameters','Fontsize',40)
set(gca,'fontsize',20)

Best_no_of_para
%% Plots all the considered microstates
figure(6)
for i=1:30
    subplot(5,6,i)
    topoplot(A(:,i),chanlocs); % first describes most
    title([num2str(i)])
    set(gca,'fontsize',20)
end
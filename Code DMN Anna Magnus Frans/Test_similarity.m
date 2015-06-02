%% Peaks and similarity
%% Here we calculate the robustness of the microstates
close all
clc
clearvars -except eegfold eeg chanlocs Y Y37DF hrf

% Initialize
gfp = std(eegfold,0,2); %GFP
steps = 3;
IC = 15; % 15 independent components (cant converge with 30)

% Counters
j = 0;
h=0
w=0


hi = 0;
% for h=linspace(min(gfp), max(gfp)/2, steps) % Min height of GFP peaks
for h=linspace(min(gfp), 2.5, steps) % Min height of GFP peaks
    hi = hi +1
    wi = 0;
    for w=ceil(linspace(0, 50, steps)) % Min width of GFP peaks 1 sample = 2 ms
        j = j +1
        wi = wi +1
        
        
        for k = 1:8; % For each fold
            j
            k
            
            idx = 1:8;
            idx(k) = []; % remove ith subject
            EEG = reshape(eeg(:,:,idx),30,[])'; % concatenate subjects
            
            gfp = std(EEG,0,2); % GFP
            [~, locs] = findpeaks(gfp,'MinPeakHeight',h, 'MinPeakWidth',w); % Find the peaks, for this design
            microstates = EEG(locs,:)'; % Pre-ICA-microstates
            [icasig, A, W] = fastica(microstates,'numOfIC', IC); % using cubic polynomial as default. IC is the number of independent components
            S = W*EEG'; % apply seperation matrix to original data (get the sources)
            
            % calculate power, to order A matrix
            % power explained by each component c
            A2 = sum(A.^2,1);
            S2 = sum(S.^2,2);
            [J T] = size(EEG');
            p = zeros(IC,1);
            for c=1:IC;
                p(c) = (A2(c) * S2(c)) / (IC*T);
            end
            power = p / sum(p);
            [~, order1] = sort(power,'descend');
            A8(:,:,k) = A(:,order1); % save all mixing matrices (ordered)
            S8(:,:,k) = S(order1,:); % order S and save
            
            % You now have the 8 A matrices for each 8-fold in A8 (ordered)
            % These might have a wrong sign X = (-A)(-S) = AS
            % Align the matrices by the very first column of A so they have correct sign
            
            % compare first row of a matrix for correlation. First row
            % explain the most (ordered above). Align by negative sign
            if corr(A8(:,1,1),A8(:,1,k)) < 0
                A8(:,:,k) = -A8(:,:,k);
                S8(:,:,k) = -S8(:,:,k);
            end
        end % end of 8-fold
        
        % Calculate similarity of these 8 folds
        % We calculate the Frobenius Norm of the deviation from the mean of
        % the 8-fold to compares across design settings
        % and take the mean of that
        
        % new 8-fold for robustness
        meanfolds = mean(A8,3); % mean of the 8 folds
        for k = 1:8 % For each fold
            deviation(:,:,k) = meanfolds - A8(:,:,k); % deviation from the mean
            frobnorm(k,hi,wi) = norm(deviation(:,:,k),'fro');
        end
        similarity(hi,wi) = sum(frobnorm(:,hi,wi)); % sum of the frobinious of deviation
    end
end

% Plot the similarity between then different design settings
figure
plot(similarity(:))
xlabel('Run'), ylabel('Sum of Frobenius Norm')
%% print table for report
% Print to latex
% printthis=[[2:11]' similarity]
printthis=[reshape(similarity(:,:),[3,3])]

clc;
d1 = digits(3); % records and sets accuracy
latex(vpa(sym(printthis))) % prints to latex
digits(d1); % restore previous accuracy


%%
close all
figure
subplot(1,2,1)
for i=1:steps
    wi=i;
    hold on;
    plot(similarity(:,wi))
    xlabel('Height Run'), ylabel('Sum of Frobenius Norm')
end

subplot(1,2,2)
for i=1:steps
    hi=i;
    hold on;
    plot(reshape(similarity(hi,:),[1 steps]))
    xlabel('Weight Run'), ylabel('Sum of Frobenius Norm')
end


%% Plot micro
close all
hi=2;
wi=3;
for k = 1:8
    figure
    for i=1:IC
        subplot(3,5,i)
        topoplot(A8(:,i,k),chanlocs); % first describes most
    end
end


%% PROJECT 2 - GRUPPO V
clear; close all; clc

%% 1) Load the data and apply pre-processing steps: filtering (1-30 Hz),
% bad channel replacement, ICA for artifact removal, average reference.

%% PRE
load('resting_pre.mat')
data_pre = EEG_import';

fs = 250; %[Hz]
T = 1/fs; %[s]
Ns = size(data_pre,1);
Nch = size(data_pre,2);

% filtering
% high pass at 1Hz
fc = 1; %Hz
[b,a] = butter(4, fc/(fs/2), 'high') ;
data_pre_filt = filtfilt(b,a,data_pre);

% low pass at 30Hz
fc = 30; %[Hz]
[b,a] = butter(4, fc/(fs/2), 'low') ;
data_pre_filt = filtfilt(b,a,data_pre_filt);

% selecting bad channels considering SD
averageSD = mean(std(data_pre_filt));
badCh = std(data_pre_filt) > 3*averageSD;

if ~any(badCh)
    disp("there aren't bad channels")
else
    chToInterpolate = find(badCh);
    disp(['bad channels to be interpolated: ', num2str(chToInterpolate)])
end

data_pre_filt = data_pre_filt';

save("data_pre_filt.mat", 'data_pre_filt')
% eeglab % interpolation, ICA, average reference

data_pre_pp = ALLEEG.data';
save("data_pre_pp.mat", 'data_pre_pp')

%% POST
load('resting_post.mat')
data_post = EEG_import';
Ns_post = size(data_post,1);

% filtering
% high pass at 1Hz
fc = 1; %Hz
[b,a] = butter(4, fc/(fs/2), 'high') ;
data_post_filt = filtfilt(b,a,data_post);

% low pass at 30Hz
fc = 30; %[Hz]
[b,a] = butter(4, fc/(fs/2), 'low') ;
data_post_filt = filtfilt(b,a,data_post_filt);

% selecting bad channels considering SD
averageSD = mean(std(data_post_filt));
badCh = std(data_post_filt) > 3*averageSD;

if ~any(badCh)
    disp("there aren't bad channels")
else
    chToInterpolate = find(badCh);
    disp(['bad channels to be interpolated: ', num2str(chToInterpolate)])
end

data_post_filt = data_post_filt';

save("data_post_filt.mat", 'data_post_filt')
% eeglab % interpolation, ICA, average reference

data_post_pp = ALLEEG.data';
save("data_post_pp.mat", 'data_post_pp')

%% 2) Compute the average spectrum

%% PRE
% loading
load('data_pre_pp.mat')
fs = 250; %[Hz]
Ns = size(data_pre_pp,1);
Nch = size(data_pre_pp,2);

frequency = 0:0.01:30;

f = frequency ;
window = 1*fs;
noverlap = 0.5*fs;

[pxx,f] = pwelch(data_pre_pp,window,noverlap,f,fs);  % PSD using Welch's overlapped segment averaging estimator.

%% POST
% loading
load('data_post_pp.mat')
Ns_post = size(data_post_pp,1);

frequency = 0:0.01:30;
f_post = frequency ;
[pxx_post,f_post] = pwelch(data_post_pp,window,noverlap,f_post,fs);  % PSD using Welch's overlapped segment averaging estimator.

%% Plot
figure()
subplot(211)
plot(f,pxx)
xlabel('frequency [Hz]')
ylabel('PSD \mu V^2 /Hz')
title('Average spectrum (before)')
subplot(212)
plot(f_post,pxx_post)
xlabel('frequency [Hz]')
ylabel('PSD \mu V^2 /Hz')
title('Average spectrum (after)')

%% 3) Compute the relative power spectra in the canonical EEG frequency bands for each channel

%% PRE
% delta
f1_delta = find(f >= 1, 1, 'first');
f2_delta = find(f >= 4, 1, 'first');

X_delta = f(f1_delta:f2_delta);
Y_delta = pxx(f1_delta:f2_delta, :);

power_delta = trapz(X_delta, Y_delta);

% theta
f1_theta = find(f >= 4, 1, 'first');
f2_theta = find(f >= 8, 1, 'first');

X_theta = f(f1_theta:f2_theta);
Y_theta = pxx(f1_theta:f2_theta, :);

power_theta = trapz(X_theta, Y_theta);

% alpha
f1_alpha = find(f >= 8, 1, 'first');
f2_alpha = find(f >= 13, 1, 'first');

X_alpha = f(f1_alpha:f2_alpha);
Y_alpha = pxx(f1_alpha:f2_alpha, :);

power_alpha = trapz(X_alpha, Y_alpha);

% beta
f1_beta = find(f >= 13, 1, 'first');
f2_beta = find(f >= 30, 1, 'first');

X_beta = f(f1_beta:f2_beta);
Y_beta = pxx(f1_beta:f2_beta, :);

power_beta = trapz(X_beta, Y_beta);

% power all
f1_all = find(f >= 1, 1, 'first');
f2_all = find(f >= 30, 1, 'first');

X_all = f(f1_all:f2_all);
Y_all = pxx(f1_all:f2_all, :);

power_all = trapz(X_all, Y_all);

% relative powers
power_delta_rel = power_delta./power_all;
power_theta_rel = power_theta./power_all;
power_alpha_rel = power_alpha./power_all;
power_beta_rel = power_beta./power_all;

%% POST
% delta
f1_delta = find(f_post >= 1, 1, 'first');
f2_delta = find(f_post >= 4, 1, 'first');

X_delta = f_post(f1_delta:f2_delta);
Y_delta = pxx_post(f1_delta:f2_delta, :);

power_delta = trapz(X_delta, Y_delta);

% theta
f1_theta = find(f_post >= 4, 1, 'first');
f2_theta = find(f_post >= 8, 1, 'first');

X_theta = f_post(f1_theta:f2_theta);
Y_theta = pxx_post(f1_theta:f2_theta, :);

power_theta = trapz(X_theta, Y_theta);

% alpha
f1_alpha = find(f_post >= 8, 1, 'first');
f2_alpha = find(f_post >= 13, 1, 'first');

X_alpha = f_post(f1_alpha:f2_alpha);
Y_alpha = pxx_post(f1_alpha:f2_alpha, :);

power_alpha = trapz(X_alpha, Y_alpha);

% beta
f1_beta = find(f_post >= 13, 1, 'first');
f2_beta = find(f_post >= 30, 1, 'first');

X_beta = f_post(f1_beta:f2_beta);
Y_beta = pxx_post(f1_beta:f2_beta, :);

power_beta = trapz(X_beta, Y_beta);

% power all
f1_all = find(f_post >= 1, 1, 'first');
f2_all = find(f_post >= 30, 1, 'first');

X_all = f_post(f1_all:f2_all);
Y_all = pxx_post(f1_all:f2_all, :);

power_all = trapz(X_all, Y_all);

% relative powers
power_delta_rel_post = power_delta./power_all;
power_theta_rel_post = power_theta./power_all;
power_alpha_rel_post = power_alpha./power_all;
power_beta_rel_post = power_beta./power_all;


%% box-plot relative power spectra
grp = [zeros(1,Nch),ones(1,Nch)];

figure
subplot(221)
boxplot([power_delta_rel, power_delta_rel_post], grp, 'Notch','on', 'Labels', {'Before EAG','After EAG'})
title('Relative power in DELTA band [1-4] Hz')
subplot(222)
boxplot([power_theta_rel, power_theta_rel_post], grp, 'Notch','on', 'Labels', {'Before EAG','After EAG'})
title('Relative power in THETA band [4-8] Hz')
subplot(223)
boxplot([power_alpha_rel, power_alpha_rel_post], grp, 'Notch','on', 'Labels', {'Before EAG','After EAG'})
title('Relative power in ALPHA band [8-13] Hz')
subplot(224)
boxplot([power_beta_rel, power_beta_rel_post], grp, 'Notch','on', 'Labels', {'Before EAG','After EAG'})
title('Relative power in BETA band [13-30] Hz')

%% statistical test relative power spectra
if ~lillietest(power_delta_rel) && ~lillietest(power_delta_rel_post) % if the variables are normally distributed
    [h_delta, p_delta] = ttest2(power_delta_rel, power_delta_rel_post);
    disp('ttest2')
else
    [p_delta, h_delta] = ranksum(power_delta_rel, power_delta_rel_post);
    disp('ranksum')
end

if ~lillietest(power_theta_rel) && ~lillietest(power_theta_rel_post) % if the variables are normally distributed
    [h_theta, p_theta] = ttest2(power_delta_rel, power_delta_rel_post);
    disp('ttest2')
else
    [p_theta, h_theta] = ranksum(power_delta_rel, power_delta_rel_post);
    disp('ranksum')
end

if ~lillietest(power_alpha_rel) && ~lillietest(power_alpha_rel_post) % if the variables are normally distributed
    [h_alpha, p_alpha] = ttest2(power_delta_rel, power_delta_rel_post);
    disp('ttest2')
else
    [p_alpha, h_alpha] = ranksum(power_delta_rel, power_delta_rel_post);
    disp('ranksum')
end

if ~lillietest(power_beta_rel) && ~lillietest(power_beta_rel_post) % if the variables are normally distributed
    [h_beta, p_beta] = ttest2(power_beta_rel, power_beta_rel_post);
    disp('ttest2')
else
    [p_beta, h_beta] = ranksum(power_beta_rel, power_beta_rel_post);
    disp('ranksum')
end


%% 4) Visualize the power spectra values for each channel through the eeglab function topoplot
eeglab

EEG.chanlocs= readlocs('gtech_64.sfp');

% delta
figure()
subplot(121)
topoplot(power_delta_rel, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('Before EAG')
sgtitle('Relative power in DELTA band [1-4] Hz')
colormap('jet')
colorbar
caxis([0 max(power_delta_rel)])
subplot(122)
topoplot(power_delta_rel_post, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('After EAG')
colormap('jet')
colorbar
caxis([0 max(power_delta_rel_post)])

% theta

figure()
subplot(121)
topoplot(power_theta_rel, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('Before EAG')
sgtitle('Relative power in THETA band [4-8] Hz')
colormap('jet')
colorbar
caxis([0 max(power_theta_rel)])
subplot(122)
topoplot(power_theta_rel_post, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('After EAG')
colormap('jet')
colorbar
caxis([0 max(power_theta_rel_post)])

% alpha

figure()
subplot(121)
topoplot(power_alpha_rel, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('Before EAG')
sgtitle('Relative power in ALPHA band [8-13] Hz')
colormap('jet')
colorbar
caxis([0 max(power_alpha_rel)])
subplot(122)
topoplot(power_alpha_rel_post, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('After EAG')
colormap('jet')
colorbar
caxis([0 max(power_alpha_rel_post)])

% beta

figure()
subplot(121)
topoplot(power_beta_rel, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('Before EAG')
sgtitle('Relative power in BETA band [13-30] Hz')
colormap('jet')
colorbar
caxis([0 max(power_beta_rel)])
subplot(122)
topoplot(power_beta_rel_post, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('After EAG')
colormap('jet')
colorbar
caxis([0 max(power_beta_rel_post)])

%% 5) Compute Higuchi fractal measure for all the electrodes and visualize its topographical distribution
klin = 6;
kmax = 35;
Lepoch = 4*fs; % 4s
overlap = Lepoch/2; % 2s

%% PRE
H = zeros(1,Nch);
for i = 1:Nch
    count = 1;

    for k = 1+overlap:overlap:Ns-overlap
        segments(count,:) = data_pre_pp(k-overlap:k+overlap-1,i);
        count = count + 1;
    end

    [Higuchi, output_lnk, output_lnLk] = featuresExtraction2(segments, klin, kmax);
    H_all(i, :) = Higuchi;
    H(i) = median(Higuchi);
    if i == 30
        figure
        plot(output_lnk(1,:), median(output_lnLk), 'r', 'LineWidth', 1)
        xlabel('log(k)')
        ylabel('log(L_k)')
        title('log(L_k) vs log(k) (Cz)')
        hold on
        
        [f,xi] = ksdensity(Higuchi); 
        
    end
end

%% POST
H_post = zeros(1,Nch);
for i = 1:Nch
    count = 1;

    for k = 1+overlap:overlap:Ns_post-overlap
        segments(count,:) = data_post_pp(k-overlap:k+overlap-1,i);
        count = count + 1;
    end

    [Higuchi, output_lnk, output_lnLk] = featuresExtraction2(segments, klin, kmax);
    H_all_post(i, :) = Higuchi;
    H_post(i) = median(Higuchi);
    if i == 30
        plot(output_lnk(1,:), median(output_lnLk), 'b', 'LineWidth', 1)
        xline(log(klin), 'k', 'LineWidth', 1)
        legend('Before EAG','After EAG', 'log(k_{lin})')

        [f_post,xi_post] = ksdensity(Higuchi); 
    end
end

figure
plot(xi,f,  'r', 'LineWidth', 1)
hold on
plot(xi_post,f_post,  'b', 'LineWidth', 1)
xlabel("Higuchi's FD")
ylabel('Density')
title("Higuchi's FD distribution (Cz)")
legend('Before EAG', 'After EAG')

%% topoplot Higuchi's FD median values
%eeglab

EEG.chanlocs= readlocs('gtech_64.sfp');

figure()
subplot(121)
topoplot(H, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('Before EAG')
sgtitle("Higuchi's FD")
colormap('jet')
colorbar
caxis([1 max(H)])
subplot(122)
topoplot(H_post, EEG.chanlocs, 'style','both','electrodes','labelpoint')
title('After EAG')
colormap('jet')
colorbar
caxis([1 max(H_post)])

%% box-plot considering median values of Higuchi's FD
grp = [zeros(1,Nch),ones(1,Nch)];

figure
boxplot([H, H_post], grp, 'Notch','on', 'Labels', {'Before EAG','After EAG'})
title("Median Higuchi's FD")

%% statistical test considering median values of Higuchi's FD
if ~ lillietest(H) && ~lillietest(H_post) % if the variables are normally distributed
    [h, p] = ttest2(H, H_post);
    disp('ttest2')
else
    [p, h] = ranksum(H, H_post);
    disp('ranksum')
end

clear h p

%% box-plot for Higuchi's FD values of all channels

for i = 1:Nch
    if length(EEG.chanlocs(i).labels) == 5
        lab{i} = EEG.chanlocs(i).labels(2:4);
    else
        lab{i} = EEG.chanlocs(i).labels(2:3);
    end
end

figure
h = boxplot(H_all(1:32,:)','Notch','on', 'Labels', lab(1:32), 'Colors', 'g', 'Symbol', 'r.');
set(h,{'linew'},{1})
hold on
h = boxplot(H_all_post(1:32,:)', 'Notch','on', 'Labels',lab(1:32), 'Colors', 'm', 'Symbol', 'r.');
set(h,{'linew'},{1})
title("Higuchi's FD")


figure
h = boxplot(H_all(33:end,:)','Notch','on', 'Labels', lab(33:end), 'Colors', 'g', 'Symbol', 'r.');
set(h,{'linew'},{1})
hold on
h = boxplot(H_all_post(33:end,:)', 'Notch','on', 'Labels',lab(33:end), 'Colors', 'm', 'Symbol', 'r.');
set(h,{'linew'},{1})
title("Higuchi's FD")
clear h
%% statistical test for Higuchi's FD of all channels
for i = 1:Nch
    if ~ lillietest(H_all(i,:)) && ~lillietest(H_all_post(i,:)) % if the variables are normally distributed
        [h(i), p(i)] = ttest2(H_all(i,:), H_all_post(i,:));
        flag(i) = 1;
        disp('ttest2')
    else
        [p(i), h(i)] = ranksum(H_all(i,:), H_all_post(i,:));
        disp('ranksum')
    end
end

figure
stem(find(flag==1), h(find(flag==1)),'r', 'LineWidth', 1);
hold on
stem(find(flag==0), h(find(flag==0)),'b', 'LineWidth', 1);
xticks(1:64)
xticklabels(lab)
legend('ttest2', 'ranksum')
title("Statistical analysis on Higuchi's FD measures")
axis tight


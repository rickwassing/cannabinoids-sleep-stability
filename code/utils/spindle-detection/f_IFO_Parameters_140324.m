function [mean_gauss_PF,mean_gauss_Peak, bandwidth_gauss, mean_gauss_AUC_specified, mean_gauss_AUC_standardrange,power_normalized_mean_FFT] = f_IFO_Parameters_140324(traces, s_Fs, name,b_ToPlot)

if nargin<4
    b_ToPlot=1;
end
%%
% Initialize
%Low Pass Filter and Downsample
[bhi] = fir1(1000, 0.1/(s_Fs/2),'low'); % at -6dB, around 0.09 Hz; 10 to 100 s of Infraslow timescale

[bhi_Antialiasing] = fir1(100, 5/(s_Fs/2),'low'); % Antialiasing filter for the resampling at 10 Hz
s_Fs_down = 10;

% forder = 100;
forder = pop_firwsord('hamming', s_Fs_down, 0.04);
[bhi3] = fir1(forder, 0.12/(s_Fs_down/2),'low'); % Filter at low order for 0.1 Hz at 10 Hz Fs

% To use for the Cheby
% bhi2 = f_getFilterISO(0.1,0.2,1,40,s_Fs);

% Filter for hilbert test
% [bhiSp] = fir1(1000, [12 15]/(s_Fs/2),'bandpass'); % at -6dB, around 0.09 Hz; 10 to 100 s of Infraslow timescale


% FFT by Bout
power_FFT = zeros(length(traces),201);
zerotofive_2000_Hz = linspace(0, 5, 2000);
ptc_shift_plot = [];
length_bouts = [];

% Power Time Course Processing
for Idx_bout = 1:length(traces)

    % Wavelet
    [powertimecourse] = f_MGT(traces{Idx_bout},s_Fs,12,16,.2);
    powertimecourse = mean(abs(powertimecourse))'; % average the power time course across all frequencies of interest
%     powertimecourse = filtfilt(bhi,1,powertimecourse);
    % Cheby II
%         powertimecourse = filtfilt(bhi2.SOSMatrix,1,powertimecourse);

% Hilbert test
    %     v_Signal = traces{Idx_bout};
    %     [bhiSp] = fir1(1000, [12 15]/(s_Fs/2),'bandpass'); 
    %     powertimecourse = filtfilt(bhiSp,1,v_Signal);
    %     powertimecourse =abs(hilbert(powertimecourse));
    %     powertimecourse = filtfilt(bhi,1,powertimecourse);

    % For the filter after (Antialiazing)
    powertimecourse = filtfilt(bhi_Antialiasing,1,powertimecourse);
    powertimecourse_onebout = resample(powertimecourse,s_Fs_down,s_Fs); %downsample
    % For the filter after
    powertimecourse_onebout = filtfilt(bhi3,1,powertimecourse_onebout);

    zerotofive_lengthofbout_Hz = linspace(0,5,(length(powertimecourse_onebout)/2)+1);
    powertimecourse_onebout_shifted = powertimecourse_onebout-mean(powertimecourse_onebout);
    power_bout = (2*abs(fft(powertimecourse_onebout_shifted))/length(powertimecourse_onebout)).^2;
    power_bout = power_bout(1:length(zerotofive_lengthofbout_Hz));
    power_bout = interp1(zerotofive_lengthofbout_Hz,power_bout,zerotofive_2000_Hz);
    power_FFT(Idx_bout,:)=power_bout(1:201);

%     ptc_shift_plot = [ptc_shift_plot powertimecourse_onebout_shifted];
%     length_bouts = [length_bouts length(powertimecourse_onebout)/s_Fs_down]; %put the length in seconds
    clear powertimecourse_onebout powertimecourse_onebout_shifted power_bout zerotofive_lengthofbout_Hz
end

freqz_power_FFT = zerotofive_2000_Hz(1:201);
power_FFT = power_FFT(:,freqz_power_FFT<=0.1);
power_FFT = power_FFT';
power_normalized_bouts_FFT =(power_FFT)./repmat(mean(abs(power_FFT(freqz_power_FFT<=0.1,:))),size(power_FFT,1),1); % RELATIVE POWER
power_normalized_bouts_FFT = power_normalized_bouts_FFT'; % bouts by frequencies
power_normalized_mean_FFT = mean(power_normalized_bouts_FFT,1);

%% Gaussian Fit for the Mean

% Account for ultra low-frequency artifact by creating artificial line between 0 and 0.0075 Hz (assume that at 0 Hz, 0 power)
for_gauss_freqz = freqz_power_FFT(freqz_power_FFT<=0.1);
% for_gauss_freqz = freqz_power_FFT;
for_gauss_mean = power_normalized_mean_FFT;
% for_gauss_mean = for_gauss_mean-mean(for_gauss_mean(for_gauss_freqz>0.05&for_gauss_freqz<0.1));% If we do a stronger filter, do we need this?
makeline = linspace(0,for_gauss_mean(1,4), 4); 
for_gauss_mean(1:4) = makeline;
threshold_foroscillationpeak = std(for_gauss_mean)*1.5;

% Fit a Gaussian 
try
    cfit_Gaussian = fit(for_gauss_freqz', for_gauss_mean', 'gauss2', 'Exclude', isnan(for_gauss_freqz));
catch
    mean_gauss_PF = nan;
    mean_gauss_Peak = nan;
    bandwidth_gauss = nan; 
    mean_gauss_AUC_specified = nan;
    mean_gauss_AUC_standardrange= nan;
    mean_freq= nan;
    disp(['Did not find a periodic component'])    
    if b_ToPlot
%         figure
%         fig = gcf; % Get the current figure handle
        titleText = sprintf('Mean Spectral Plot Properties for %s', name);
%         annotation(fig, 'textbox', [0.1, 0.93, 0.8, 0.05], 'String', titleText, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        plot(for_gauss_freqz, for_gauss_mean, 'k') %raw data
        makeline(5:40) = NaN;
        makeline(1)=0;
        hold on
        plot(for_gauss_freqz, makeline, 'b')
        xline(0.0075, 'k--');
        legend('Raw', 'Artificial Line', 'Inclusion Limit','box','off')
        title('Compare Gaussians to Raw Data, not fond')
        xlim([0, 0.1])
        hold off
        set(gca,'Box','off','TickDir','out')
%         cd('L:\Somnus-Data\User\Maria_Dimitriades\Projects\Infraslow - Schizophrenia\Participants\allcontrols\plots\Updated_MeanFit_16224\')
%         ov = ['IFO_Mean_' name '.png']
%         saveas(1, ov)
%         close all
    end
    return;
end
gauss_mean = cfit_Gaussian(for_gauss_freqz');

coefficients = coeffvalues(cfit_Gaussian);
std_dev = coefficients(3);
mean_freq = coefficients(2);
PeakofGaussian = coefficients(1);
if PeakofGaussian < threshold_foroscillationpeak % THERE IS NOT PERIODIC COMPONENT    
    bandwidth_gauss = nan;
    mean_gauss_PF = nan;
    mean_gauss_Peak = nan;
    mean_gauss_AUC_specified = nan;
    mean_gauss_AUC_standardrange= nan;
    mean_freq= nan;
    disp(['Did not find a periodic component'])    
    % %Plot
    if b_ToPlot
%         figure
        fig = gcf; % Get the current figure handle
        titleText = sprintf('Mean Spectral Plot Properties for %s', name);
%         annotation(fig, 'textbox', [0.1, 0.93, 0.8, 0.05], 'String', titleText, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        plot(for_gauss_freqz, gauss_mean, 'g'); % Gaussian 1
        hold on
        power_normalized_mean_FFT(1:3) = NaN;
        plot(for_gauss_freqz, for_gauss_mean, 'k') %raw data
        makeline(5:40) = NaN;
        makeline(1)=0;
        plot(for_gauss_freqz, makeline, 'b')
        xline(0.0075, 'k--');
        yline(threshold_foroscillationpeak, 'm--');
        legend('Data with Gaussian 2','Raw', 'Artificial Line', 'Inclusion Limit','Threshold','box','off')
        title('Compare Gaussians to Raw Data, not fond')
        xlim([0, 0.1])
        hold off
        set(gca,'Box','off','TickDir','out')
%         cd('L:\Somnus-Data\User\Maria_Dimitriades\Projects\Infraslow - Schizophrenia\Participants\allcontrols\plots\Updated_MeanFit_16224\')
%         ov = ['IFO_Mean_' name '.png']
%         saveas(1, ov)
%         close all
    end
    return;
end
upper_limit = std_dev + mean_freq;
lower_limit = mean_freq - std_dev; % 2 standard deviations are 95%; 1 is 68%; here we chose 1 
bandwidth_gauss = upper_limit - lower_limit;

[pks_mean,locs_mean,~,~] = findpeaks(gauss_mean,for_gauss_freqz,'Annotate','extents');

mean_gauss_PF = [];
mean_gauss_Peak = [];

if sum(pks_mean) >= 1
    mean_gauss_PF=locs_mean(pks_mean==max(pks_mean));
    mean_gauss_Peak=pks_mean(pks_mean==max(pks_mean));
else
    mean_gauss_PF = NaN;
    mean_gauss_Peak=nan;
end


% The lower bound should be above 0.0074 to calculate AUC
if lower_limit < 0.0075
    lb = 0.0075;
elseif lower_limit >= 0.0075
    lb = lower_limit;
end

[~, indexL] = min(abs(for_gauss_freqz - lb));
lb= for_gauss_freqz(indexL);

[~, indexU] = min(abs(for_gauss_freqz - upper_limit));
upper_limit = for_gauss_freqz(indexU);

mean_gauss_AUC_specified = trapz(gauss_mean(for_gauss_freqz>=lb & for_gauss_freqz<=upper_limit));
mean_gauss_AUC_standardrange= trapz(gauss_mean(for_gauss_freqz>=0.0075 & for_gauss_freqz<=0.04));

% %Plot
if b_ToPlot
%     figure
    fig = gcf; % Get the current figure handle
    titleText = sprintf('Mean Spectral Plot Properties for %s', name);
%     annotation(fig, 'textbox', [0.1, 0.93, 0.8, 0.05], 'String', titleText, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    plot(for_gauss_freqz, gauss_mean, 'g'); % Gaussian 2
    hold on
    power_normalized_mean_FFT(1:3) = NaN;
    plot(for_gauss_freqz, for_gauss_mean, 'k') %raw data
    makeline(5:40) = NaN;
%     makeline(5:201) = NaN;
    makeline(1)=0;
    plot(for_gauss_freqz, makeline, 'b')
    xline(0.0075, 'k--');
    plot(locs_mean,pks_mean,'o')
    yline(threshold_foroscillationpeak, 'm--');
    xline(lower_limit, 'r--');
    xline(upper_limit, 'r--');
    legend('Data with Gaussian 2','Raw', 'Artificial Line', 'Inclusion Limit', 'Peak Frequency', 'Threshold','box','off')
    title('Compare Gaussians to Raw Data')
    xlim([0, 0.1])
    hold off
        set(gca,'Box','off','TickDir','out')
%     cd('L:\Somnus-Data\User\Maria_Dimitriades\Projects\Infraslow - Schizophrenia\Participants\allcontrols\plots\Updated_MeanFit_16224\')
%     ov = ['IFO_Mean_' name '.png']
%     saveas(1, ov)
%     close all
end

end


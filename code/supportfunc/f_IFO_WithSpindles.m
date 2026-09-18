function [st_Abs,st_Rel, v_SpPhase] = f_IFO_WithSpindles(m_Traces, s_Fs , v_TimeBouts , v_TimeSpindles , str_name , b_ToPlot)
    % f_IFO is the basic function to do the infraslow anaylsis. It uses the
    % dynamics of sigma from the f_MGT between 12 and 16.
    % Inputs:
    %   m_Traces: Cell array that contians the EEG within N2 with at least 
    %               300 s of continous EEG. Including MA. 
    %
    %   s_Fs: Sampling frequency (in my case 250 Hz). 
    %
    %   v_TimeBouts: Cell array containing the time of the bouts. 
    %
    %   v_TimeSpindles: Vector containing the time at the maximum amplitude of
    %                 the spindles.
    %
    %   str_name: String with name of the file just used to do the title of 
    %               the figures whenever you want to plot.
    %
    %   b_ToPlot: boolean variable (1/0, true/false) if you want to plot.
    %
    % Outputs:
    %   st_Abs: Information of the infraslow metrics for the absolute
    %           spectral power.
    %           
    %   st_Rel: Information of the infraslow metrics for the relative
    %           spectral power.
    %
    %   v_SpPhase: Phase of spindles within the Infraslow Oscillation of
    %           sigma.
    %       
    %   AO: Sept. 2024

    if nargin < 4
        b_ToPlot = 1;
    end

    %% Initialization
    % Filter infraslow

    [bhi] = fir1(100, 0.1/(s_Fs/2),'low'); % at -6dB, around 0.09 Hz; 10 to 100 s of Infraslow timescale

    % Set up the frequency vector
    v_f250 = linspace(0,s_Fs/2,(s_Fs*300/2)+1); % Freuency vector for the welch.
    power_FFT_at250Hz = zeros(length(m_Traces),find(v_f250>=0.1,1));
    powerSp_FFT_at250Hz = zeros(length(m_Traces),find(v_f250>=0.1,1));
    v_SpPhase = [];

    %% Infraslow dynamics per bout
    for Idx_bout = 1:length(m_Traces)
        
        v_CurrTrace = m_Traces{Idx_bout};
        v_CurrTime = v_TimeBouts{Idx_bout};
        v_CurrSp = find(v_TimeSpindles>v_CurrTime(1)&v_TimeSpindles<v_CurrTime(end));
        
        v_CurrSp = v_TimeSpindles(v_CurrSp);
        v_IdxSpn = zeros(size(v_CurrSp));
        
        for idx = 1:length(v_IdxSpn)
            v_IdxSpn(idx) = find(v_CurrTime>=v_CurrSp(idx),1,"first");
        end
        
        %%% Signal with spindles %%% 
        for idxSp = 1 : length(v_CurrSp);
            v_Sp1and0(v_CurrTime==v_CurrSp(idxSp))=1; %#ok<AGROW>
        end
        v_SpSignal = movsum(v_Sp1and0,s_Fs*5);
        v_SpSignal = filtfilt(bhi,1,v_SpSignal);

        v_SpSignal = v_SpSignal - mean(v_SpSignal); 
        [v_powerSpSignal_FFT,v_f250] = pwelch(v_SpSignal,s_Fs*300,s_Fs*150,(s_Fs*300)+1,s_Fs);
        v_powerSpSignal_FFT = v_powerSpSignal_FFT(v_f250<=0.1);
        v_powerSpSignal_FFT(v_f250<=0.0075)=mean(v_powerSpSignal_FFT(v_f250>0.06&v_f250<0.1));
        
        powerSp_FFT_at250Hz(Idx_bout,:)=v_powerSpSignal_FFT';


        % Get the bout's sigma dynamics
        [v_powertimecourse] = f_MGT(v_CurrTrace,s_Fs,12,16,.2); % Wavelet
        
        v_powertimecourse = mean(abs(v_powertimecourse))'; % Average the power time course across all frequencies of interest   
        v_powertimecourse = v_powertimecourse - mean(v_powertimecourse);
        v_CurrPhase = angle(hilbert(v_powertimecourse));

        v_SpPhase=[v_SpPhase;v_CurrPhase(v_IdxSpn)];
    
        % Welch analysis
        [v_power_FFT,v_f250] = pwelch(v_powertimecourse,s_Fs*300,s_Fs*150,(s_Fs*300)+1,s_Fs);
        v_power_FFT=v_power_FFT(v_f250<=0.1);       
        v_power_FFT(v_f250<=0.0075)=mean(v_power_FFT(v_f250>0.06&v_f250<0.1));
        v_f250=v_f250(v_f250<=0.1);
        power_FFT_at250Hz(Idx_bout,:)=v_power_FFT;
        clear v_power_FFT

    end
    


    %% Gaussian fit and parameters
    % Absolute power
    if b_ToPlot
        subplot(1,2,1)
    end
    
    v_BatchOfBouts = mean(power_FFT_at250Hz);
    [mean_gauss_PF,mean_gauss_Peak,bandwidth_gauss,mean_gauss_AUC_specified] = f_GausseanFitting(v_f250,v_BatchOfBouts,str_name,'Abs. Pow',b_ToPlot);
    
    st_Abs.mean_gauss_PF = mean_gauss_PF;
    st_Abs.mean_gauss_Peak = mean_gauss_Peak;
    st_Abs.bandwidth_gauss = bandwidth_gauss;
    st_Abs.mean_gauss_AUC = mean_gauss_AUC_specified;
    st_Abs.v_MeanSpectra = v_BatchOfBouts;
    st_Abs.m_IndSpectra = power_FFT_at250Hz;
    st_Abs.m_IndSpectraSp = powerSp_FFT_at250Hz;

    
    %     % Relative
    if b_ToPlot
        subplot(1,2,1)       
    end

    power_FFT_at250Hz = power_FFT_at250Hz';
    Relpower_FFT_at250Hz =(power_FFT_at250Hz)./repmat(mean(power_FFT_at250Hz(v_f250>.0075,:)),size(power_FFT_at250Hz,1),1); % RELATIVE POWER
    Relpower_FFT_at250Hz = Relpower_FFT_at250Hz';
    v_BatchOfBouts = mean(Relpower_FFT_at250Hz);      
    
    [mean_gauss_PF,mean_gauss_Peak,bandwidth_gauss,mean_gauss_AUC_specified] = f_GausseanFitting(v_f250,v_BatchOfBouts,str_name,'Rel. Pow',b_ToPlot);
    
    st_Rel.mean_gauss_PF = mean_gauss_PF;
    st_Rel.mean_gauss_Peak = mean_gauss_Peak;
    st_Rel.bandwidth_gauss = bandwidth_gauss;
    st_Rel.mean_gauss_AUC = mean_gauss_AUC_specified;
    st_Rel.v_MeanSpectra = v_BatchOfBouts;
    st_Rel.m_IndSpectra = Relpower_FFT_at250Hz;



    powerSp_FFT_at250Hz = powerSp_FFT_at250Hz';

    RelpowerSp_FFT_at250Hz =(powerSp_FFT_at250Hz)./repmat(mean(powerSp_FFT_at250Hz(v_f250>.0075,:)),size(powerSp_FFT_at250Hz,1),1); % RELATIVE POWER
    RelpowerSp_FFT_at250Hz = RelpowerSp_FFT_at250Hz';
    v_BatchOfBouts = mean(RelpowerSp_FFT_at250Hz);   

    st_Rel.m_IndSpectraSp = RelpowerSp_FFT_at250Hz;

end
% In case one wants to use FFT instad of Welch we can use this
function v_power_FFT = f_DoSpectralAnalysisForBout(powertimecourse_onebout,s_Fs,zerotoHalf_ForInterpolation)
    % Create vector of Freq for bout
    zerotoHalfFreq_lengthofbout_Hz = linspace(0,s_Fs/2,(length(powertimecourse_onebout)/2)+1);
    % Simply zero mean
    powertimecourse_onebout_shifted = powertimecourse_onebout-mean(powertimecourse_onebout);
    % SimpleFFT
    power_bout = (2*abs(fft(powertimecourse_onebout_shifted))/length(powertimecourse_onebout)).^2;
    % Interpolation
    power_bout = power_bout(1:length(zerotoHalfFreq_lengthofbout_Hz));
    % Only if you want to eliminate the begining
    power_bout = interp1(zerotoHalfFreq_lengthofbout_Hz,power_bout,zerotoHalf_ForInterpolation);
    v_power_FFT =power_bout(1:length(zerotoHalf_ForInterpolation));
end


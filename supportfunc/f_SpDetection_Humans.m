function [m_Spindles,st_Spindles] = f_SpDetection_Humans(v_Signals, s_Fs, v_Hyp, varargin)
% f_SpindleDetection
% This function detects sleep spindles using the method of envelope and
% the method of LF (Fernandez et al., 2018) with the peaks. With this
% function either you include all the four parameters for or select a
% file to detect spindles.
%
% Inputs:
%   v_Signals:
%       Signal to detect spindles.
%
%   s_Fs:
%       Sampling frequency.
%
%   v_Hyp: Hypnogram vector one number per point of data:
%       W:   0
%       N1: -1
%       N2: -2
%       N3: -3
%       R:  -4
%
%   s_Th_SD: Standard deviation threshold over the mean. (default: 1.5)
%
%   s_NCycles: Minimum number of cycles for spindles (Default: 5).
%
%   s_TimeToMerge: Maximum time to merge spindles (Default 0.5 s)
%
% Outputs:
%   m_Spindles:
%       Matrix (2,n) containing the indices of initial and final
%       positions of the "n" spindles detected.
%	v_Sp_Amp:
%       Spindle amplitudes
%   v_Sp_Speed:
%       Spindle speed in hz
%   v_Sp_NCycles:
%       Number of cycles
%   v_Sp_Time:
%       length of spindles.
%   v_Sp_Pow:
%       Power of spindles
%   v_Sp_InterDistance:
%       Inter Distance time of spindles
%
%   See also f_PhaseAnalysis, f_FFTPowAsLab, f_FFTFromBouts
%
% Original:
% Alejandro Osorio-Forero 2019 from the animals in the Luthi Lab-
% 2023 Accomodated it for the humans in the MS project
%
% Updated 26.11.2023:
% Only use the activity of N2 and apply it to both N2 and N3

s_Th = 1.5; % In SD if b_DefinedThreshold = 0; in absolute number if b_DefinedThreshold=1
s_NCycles = 5; % Number of minimum cycles to detect spindles
s_TimeToMerge = 1/9; % Time in seconds of spindles to merge.
v_Boundries = [10, 16]; % Frequency cutoffs for bandpass filter
b_DefinedThreshold=0;
if varargin{2} == 9999 % Pz
    b_DoPlot = true;
else
    b_DoPlot = false;
end

v_Signals   = v_Signals-nanmean(v_Signals);
if size(v_Signals,2)<size(v_Signals,1)
    v_Signals=v_Signals';
end
v_PosNan = find(isnan(v_Signals));
v_Signals(v_PosNan) = 0;

disp('Filtering...')

% Create filter
forder = round(0.70383+s_Fs*1.6502);
bhi_spfreqs = fir1(forder, v_Boundries/(s_Fs/2));
bhi_hpfreqs = fir1(forder, [10, 30]/(s_Fs/2));

% Arrange trace
v_MovStd = movstd(v_Signals, 3*s_Fs);
v_IsGoodSignal = v_MovStd <= prctile(v_MovStd, 99);
v_Filtered = filtfilt(bhi_spfreqs, 1, v_Signals);
v_HiPassFilt = filtfilt(bhi_hpfreqs, 1, v_Signals);
v_FilteredSquared = v_Filtered.^2;
% Only use the activity of N2 and apply it to both N2 and N3
%     v_InNR = find(v_Hyp==-2|v_Hyp==-3);
v_InNR = find(v_Hyp==-2 & v_IsGoodSignal);

% Fernandez et al., 2018's method
disp('Detecting...')
if b_DefinedThreshold
    s_Threshold = s_Th; % Threshold based only in sleep signal
else
    s_Threshold = mean(v_FilteredSquared(v_InNR))+s_Th*std(v_FilteredSquared(v_InNR)); % Threshold based only in sleep signal
end
st_Spindles.s_Threshold = s_Threshold;

v_InNR = find(v_Hyp==-2|v_Hyp==-3); % Apply it to both N2 and N3
[v_PeakVal,v_PeakLoc] = findpeaks(v_FilteredSquared); % Find the peaks
[v_PeakLoc,v_4Val,~] = intersect(v_PeakLoc,v_InNR);   % Keep only peaks in NR
v_PeakVal=v_PeakVal(v_4Val);
%     v_AllPeaksVal = v_PeakVal;
v_AllPeaksPos = v_PeakLoc;

if isscalar(s_Threshold)
    v_PeakLoc(v_PeakVal<s_Threshold)=[];              % Keep only peaks above the threshold
    v_PeakVal(v_PeakVal<s_Threshold)=[];
else % is adaptive threshold vector
    rm_idx = asrow(v_PeakVal) < asrow(s_Threshold(v_PeakLoc));
    v_PeakLoc(rm_idx)=[];              % Keep only peaks above the threshold
    v_PeakVal(rm_idx)=[];
end
[~,v_DropsLoc] = findpeaks(v_FilteredSquared*-1); % Find the peaks
if v_PeakLoc(1)<v_DropsLoc(1)
    v_DropsLoc=[1 v_DropsLoc];
end
% To get the frequencies
v_SecondsBetweenPeaks = (v_PeakLoc(2:end)-v_PeakLoc(1:end-1))/s_Fs;
v_TooLong = find(v_SecondsBetweenPeaks>1/(2*16));

% First, based on the long distances get the start and end of the
% spindles
try
    v_PosibleEnd = [v_PeakLoc(v_TooLong);v_PeakLoc(end)];
catch
    v_PosibleEnd = [v_PeakLoc(v_TooLong),v_PeakLoc(end)]';
end
try
    v_PosibleStart = [ v_PeakLoc(1);v_PeakLoc(v_TooLong+1)];
catch
    v_PosibleStart = [ v_PeakLoc(1),v_PeakLoc(v_TooLong+1)]';
end
% Get the two cycles before and after the starts for real start and end
v_RealStart = zeros(2,length(v_PosibleStart));
v_RealEnd = v_RealStart;
for idxFind = 1:length(v_PosibleStart)
    if ~isempty(find(v_DropsLoc<v_PosibleStart(idxFind),2,'last'))
        v_RealStart(:,idxFind) = v_DropsLoc(find(v_DropsLoc<v_PosibleStart(idxFind),2,'last'));
    end
    if ~isempty(find(v_DropsLoc>v_PosibleEnd(idxFind),2,'first'))
        v_RealEnd(:,idxFind)   = v_DropsLoc(find(v_DropsLoc>v_PosibleEnd(idxFind),2,'first'));
    end
end
v_RealEnd(1,:)=[];
v_RealStart(2,:)=[];


% Time in sec between spindles.
v_TimeBetweenSpindlesLF = (v_RealStart(2:end)-v_RealEnd(1:end-1))/s_Fs;
v_MergeHere = find(v_TimeBetweenSpindlesLF<s_TimeToMerge); % Time in sec between spindles.
v_RealEnd(v_MergeHere)=[];
v_RealStart(v_MergeHere+1)=[];

v_AllSpLF = zeros(1,length(v_Signals));
for idx = 1:numel(v_RealStart)
    v_AllSpLF(v_RealStart(idx):v_RealEnd(idx))=max(v_Signals(v_RealStart(idx):v_RealEnd(idx)))*.9;
end

% Now only keep those spindles to which the internal frequency is
% between 9 and 16 Hz AND that there are more than 3 cycles
idxLF = 1;
while idxLF <=length(v_RealEnd)
    v_PeaksInSpindle = find(v_AllPeaksPos>v_RealStart(idxLF)&v_AllPeaksPos<v_RealEnd(idxLF));
    v_Speeds = (v_AllPeaksPos(v_PeaksInSpindle(2:end))-v_AllPeaksPos(v_PeaksInSpindle(1:end-1)))/s_Fs;
    if length(v_Speeds)>s_NCycles*2+1 % more than 5 cycles
        v_Speeds = 1./v_Speeds;
        v_Speeds = mean(v_Speeds);
        if (v_RealEnd(idxLF) - v_RealStart(idxLF)) > 3*s_Fs
            % Spindles must be less than 3 seconds
            v_RealStart(idxLF)=[];
            v_RealEnd(idxLF)=[];
        elseif v_Speeds <9*2||v_Speeds>16*2
            v_RealStart(idxLF)=[];
            v_RealEnd(idxLF)=[];
        else
            idxLF = idxLF + 1;
        end
    else
        v_RealStart(idxLF)=[];
        v_RealEnd(idxLF)=[];
    end
end

v_SpLF = zeros(1,length(v_Signals));
for idx = 1:numel(v_RealStart)
    v_SpLF(v_RealStart(idx):v_RealEnd(idx))=max(v_Signals(v_RealStart(idx):v_RealEnd(idx)))*1.2;
end
m_Spindles  = [v_RealStart;v_RealEnd];% In the points depending on the sample frequency 1000 or 200

% Features
disp('Calculating features...')
v_Sp_Amp = nan(1,length(m_Spindles));
v_Sp_Speed = nan(1,length(m_Spindles));
v_Sp_NCycles = nan(1,length(m_Spindles));
v_Sp_Time = nan(1,length(m_Spindles));
v_Sp_Pow =nan(1,length(m_Spindles));

for idxSp =1:size(m_Spindles,2)
    v_SingleSpindle = v_Filtered(m_Spindles(1,idxSp):m_Spindles(2,idxSp));
    v_SingleSpindle = v_SingleSpindle - mean(v_SingleSpindle);

    v_Sp_Amp(idxSp) = max(abs(v_SingleSpindle));

    [~,v_PosPeak]       = findpeaks(v_SingleSpindle);
    v_TempSpeed  = (v_PosPeak(2:end)-v_PosPeak(1:end-1))/s_Fs;

    v_Sp_NCycles(idxSp) = numel(v_PosPeak);
    v_Sp_Speed(idxSp)   = 1/mean(v_TempSpeed);% intraspindle frequency
    v_Sp_Time(idxSp)    = ((m_Spindles(2,idxSp)-m_Spindles(1,idxSp))/s_Fs);% How long is the spindle
    v_Sp_Pow(idxSp)     = sum((v_SingleSpindle).^2);
    v_SpindlesLoc(idxSp)     = round((m_Spindles(2,idxSp)+m_Spindles(1,idxSp))/2);
    %         v_Sp_Speed(idxSp)   = v_Sp_NCycles(idxSp)/v_Sp_Time(idxSp);
end


if b_DoPlot
    % Apply broadband filter
    forder = pop_firwsord('hamming', s_Fs, 2);
    bhi_spfreqs = fir1(forder, [1, 30]/(s_Fs/2));
    v_Times = 1/s_Fs:1/s_Fs:length(v_Signals)/s_Fs;
    v_BroadBandFilt = filtfilt(bhi_spfreqs, 1, v_Signals);
    % Initalize figure
    close all
    Fig = figure();
    Fig.Position = [1 200 1920 200];
    Fig.Color = [.61 .62 .64];
    Ax = axes(Fig);
    Ax.NextPlot = 'add';
    Ax.Color = Fig.Color;
    Ax.Layer = 'top';
    Ax.TickLength = [0, 0];
    Ax.Box = 'on';
    Ax.XGrid = 'on';
    Ax.GridColor = 'w';
    Ax.GridAlpha = 0.4;
    Ax.YTick = [];
    Ax.YLim = [-20, 300];
    Ax.XTick = unique(round(v_Times));
    Ax.XLim = [0, 30];
    % Plot
    plot(v_Times, v_FilteredSquared, '-', 'Color', [0.9 .9 .9])
    plot(v_Times, (v_BroadBandFilt)+200, '-k')
    plot([0, max(v_Times)], [s_Threshold, s_Threshold], ':r')
    plot(v_Times(v_RealStart), v_FilteredSquared(v_RealStart), 'or')
    plot(v_Times(v_RealEnd), v_FilteredSquared(v_RealEnd), 'or')
    plot(v_Times(v_SpindlesLoc), v_FilteredSquared(v_SpindlesLoc), '.r', 'MarkerSize', 6)
    text(v_Times(round(mean(m_Spindles))), ones(1, length(v_SpindlesLoc)), arrayfun(@(str) sprintf('%.1f', str), v_Sp_Speed, 'UniformOutput', false), 'HorizontalAlignment', 'center')
    v_IsNrem = nan(size(v_Signals));
    v_IsNrem(v_InNR) = -15;
    plot(v_Times, v_IsNrem, '-', 'Color', [.3 .4 .5], 'LineWidth', 3);
    % Shift the x-axis limits to view each spindle and export to disk
    Ax.OuterPosition = [0.01, 0.01, 0.98, 0.98];
    for s = 1:size(m_Spindles, 2)
        Ax.XLim = [0, 30] + m_Spindles(1, s)/s_Fs - 15;
        kv = filename2struct(varargin{1}.setname);
        kv.desc = sprintf('spindle%i', s);
        kv.filetype = 'fernandez.png';
        filepath = sprintf('%s/images/spindledet/fernandez', varargin{1}.filepath);
        fullfilepath = sprintf('%s/%s', filepath, struct2filename(kv));
        if exist(filepath, 'dir') == 0
            mkdir(filepath)
        end
        exportgraphics(Fig, fullfilepath, 'Resolution', 144)
    end
end
%% Save
if nargout == 2
    disp('Saving Spindles...')
    %     st_Spindles.m_SpindlesIn200Hz       = m_Spindles;
    st_Spindles.v_Sp_Amp                = v_Sp_Amp;
    st_Spindles.v_Sp_Speed              = v_Sp_Speed;
    st_Spindles.v_Sp_NCycles            = v_Sp_NCycles;
    st_Spindles.v_Sp_Time               = v_Sp_Time;
    st_Spindles.v_Sp_Pow                = v_Sp_Pow;
    %          st_Spindles.v_Sp_InterDistance      = v_AllDistances;
    %          st_Spindles.v_BetweenSpDistance     = v_AllDistances;
end
end


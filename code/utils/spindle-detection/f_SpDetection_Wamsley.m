function [m_Spindles, v_Duration, v_SpindleFreqs] = f_SpDetection_Wamsley(v_Signals, s_Fs, v_Hyp, varargin)

m_Spindles = 0;
v_Duration = 0;
v_SpindleFreqs = 0;

if varargin{2} == 9999 % Pz
    b_DoPlot = true;
else
    b_DoPlot = false;
end

v_Signals = detrend(v_Signals, 0); % Demean
forder = pop_firwsord('hamming', s_Fs, 2);
bhi = fir1(forder, [1, 30]/(s_Fs/2));
v_BroadBandFilt = filtfilt(bhi, 1, v_Signals);
bhi = fir1(forder, [10, 16]/(s_Fs/2));
v_SpindleFilt = filtfilt(bhi, 1, v_Signals);

fc=13.5;
n=7;
ss=n/(2*pi*fc);
tp=2*ss^2;
tStart = -4;
tStop = 4;
timeVector = linspace(tStart,tStop, (tStop-tStart)*s_Fs );
psiWavelet = (pi*tp)^(-0.5).*exp(2*1i*pi*fc.*timeVector).*exp(-timeVector.^2/tp);

input = psiWavelet;
Nfft = 10 * 2^nextpow2(length(input));
psd = 20.*log10(fftshift(abs(fft(input,Nfft))));
freqs = [0:Nfft - 1].*(s_Fs/Nfft);
freqs(freqs >= s_Fs/2) = freqs(freqs >= s_Fs/2) - s_Fs;
freqs = fftshift(freqs);

% computing baseline
t=[1:length(v_BroadBandFilt)]/s_Fs;
tt=t(~isnan(v_BroadBandFilt));
s=v_BroadBandFilt(~isnan(v_BroadBandFilt));
hyps2=v_Hyp(~isnan(v_BroadBandFilt));
fsig = conv(s,psiWavelet,'same');
fsig = abs(real(fsig));
v_MovStd = movstd(fsig, 3*s_Fs);
s_ThBaselineData = prctile(v_MovStd, 99);
ssig=smooth(fsig, 0.1*s_Fs);
baseline=mean(ssig(hyps2==-2 & v_MovStd < s_ThBaselineData), 'omitnan');

% This line used to be
% >> th=2*baseline;
% But the paper described this needed to be 4.5 times the baseline
th = 4.5*baseline;
th_upper = 12*baseline;
cores = ascolumn(ssig) > th & ascolumn(hyps2) <= -2;

% detect consecutive ones between 0.3 and 3 s
coress=num2str(ascolumn(cores));
%pattern=num2str(ones(round(s_Fs*0.3),1));
pattern=num2str(1);
beg=strfind(coress', asrow(pattern)); %start points of each candidate core
b1=[0;diff(beg')];
beg(b1==1)=[];

% store duration of each core
dur=zeros(size(beg));
for jj=1:length(beg)
    k=1;
    while(cores(beg(jj)+k)==1 && length(cores)>beg(jj)+k)
        k=k+1;
    end
    dur(jj)=k;
end
beg(dur>3*s_Fs)=[];
dur(dur>3*s_Fs)=[];

% now extend the spindles: at least 0.5 s above th2=2*baseline
th2=2*baseline;
candidates= ascolumn(ssig) > th2 & ascolumn(hyps2) <= -2;

% detect consecutive ones greater than 0.5 s
candidatess=num2str(candidates)';
pattern=num2str(ones(round(s_Fs*0.5),1));
begc=strfind(candidatess,pattern'); %start points of each candidate core
b1=[0;diff(begc')];
begc(b1==1)=[];

% store duration of each extension
durc=zeros(size(begc));
for jj=1:length(begc)
    k=1;
    while(candidates(begc(jj)+k)==1 && length(candidates)>begc(jj)+k)
        k=k+1;
    end
    durc(jj)=k;
end

% now combine the two classifications
C=zeros(size(cores));
for jj=1:length(beg)
    C(beg(jj):beg(jj)+dur(jj))=1;
end

E=zeros(size(cores));
for jj=1:length(begc)
    E(begc(jj):begc(jj)+durc(jj))=1;
end
% intersect
EC=E&C;
% adjust duration
b1=diff(EC');
begec=find(b1==1);

ind=zeros(size(begec));

for jj=1:length(begec)
    truebeg=find(begc<=begec(jj));
    if ~isempty(truebeg)
        ind(jj)=truebeg(end);
    else ind(jj)=NaN;
    end
end
ind(isnan(ind))=[];
BEG=begc(ind);
DUR=durc(ind);

% now that I found the spindles, next step is merging those too close
% (closer than 0.5 sec)
EPTS=BEG+DUR;
distance=-EPTS(1:end-1)+BEG(2:end);
while (any(distance<(0.5*s_Fs)))
    ind=find(distance<(0.5*s_Fs));
    ind=ind(1);
    if (DUR(ind)+distance(ind)+DUR(ind+1))<3*s_Fs
        BEG(ind+1)=[];
        DUR(ind)=DUR(ind)+distance(ind)+DUR(ind+1);
        DUR(ind+1)=[];
        distance=BEG(2:end)-(BEG(1:end-1)+DUR(1:end-1)); %update distance
    else
        distance(ind)=999; %do not merge
    end
end

% Check again the spindles are not more than 3 seconds
BEG(DUR > 3*s_Fs) = [];
DUR(DUR > 3*s_Fs) = [];

FIN = BEG+DUR;

% Check for outlier spindles (artefactual ones)
idx_rm = [];
for s = 1:length(BEG)
    if any(ssig(BEG(s):FIN(s)) > th_upper)
        idx_rm = [idx_rm, s];
    end
end
BEG(idx_rm) = [];
DUR(idx_rm) = [];
FIN(idx_rm) = [];

m_Spindles = [asrow(BEG); asrow(FIN)];
v_Duration = DUR;

for s = 1:size(m_Spindles, 2)
    v_SingleSpindle = v_SpindleFilt(m_Spindles(1,s):m_Spindles(2,s));
    v_SingleSpindle = v_SingleSpindle - mean(v_SingleSpindle);

    [~,v_PosPeak] = findpeaks(v_SingleSpindle);
    v_TempSpeed = (v_PosPeak(2:end)-v_PosPeak(1:end-1))/s_Fs;

    v_SpindleFreqs(s)   = 1/mean(v_TempSpeed);% intraspindle frequency
end


if b_DoPlot
    % Apply broadband filter
    forder = pop_firwsord('hamming', s_Fs, 2);
    bhi_spfreqs = fir1(forder, [1, 30]/(s_Fs/2));
    v_Times = 1/s_Fs:1/s_Fs:length(v_BroadBandFilt)/s_Fs;
    v_BroadBandFilt = filtfilt(bhi_spfreqs, 1, v_BroadBandFilt);
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
    Ax.YLim = [-30, 1000];
    Ax.XTick = unique(round(v_Times));
    Ax.XLim = [0, 30];
    % Plot
    plot(v_Times, ssig, '-', 'Color', [0.9 .9 .9])
    plot(v_Times, (v_BroadBandFilt.*2)+750, '-k')
    plot([0, max(v_Times)], [th, th], ':r')
    plot([0, max(v_Times)], [th2, th2], ':r')
    plot(v_Times(BEG), ssig(BEG), 'or')
    plot(v_Times(FIN), ssig(FIN), 'or')
    text(v_Times(round(mean(m_Spindles))), ones(1, size(m_Spindles, 2)), arrayfun(@(str) sprintf('%.1f', str), v_SpindleFreqs, 'UniformOutput', false), 'HorizontalAlignment', 'center')
    v_IsNrem = nan(size(v_BroadBandFilt));
    v_IsNrem(v_Hyp <= -2) = -15;
    plot(v_Times, v_IsNrem, '-', 'Color', [.3 .4 .5], 'LineWidth', 3);
    % Shift the x-axis limits to view each spindle and export to disk
    Ax.OuterPosition = [0.01, 0.01, 0.98, 0.98];
    for s = 1:size(m_Spindles, 2)
        Ax.XLim = [0, 30] + m_Spindles(1, s)/s_Fs - 15;
        kv = filename2struct(varargin{1}.setname);
        kv.desc = sprintf('spindle%i', s);
        kv.filetype = 'wamsley.png';
        filepath = sprintf('%s/images/spindledet/wamsley', varargin{1}.filepath);
        fullfilepath = sprintf('%s/%s', filepath, struct2filename(kv));
        if exist(filepath, 'dir') == 0
            mkdir(filepath)
        end
        exportgraphics(Fig, fullfilepath, 'Resolution', 144)
    end
end

end
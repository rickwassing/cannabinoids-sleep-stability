close all

sub = 'r017';
ses = 'etc120';

EEG = LoadDataset(['derivatives/EEG-preproc/sub-', sub, '/ses-', ses, '/sub-', sub, '_ses-', ses, '_task-psg_desc-preproc_eeg.set'], 'all');
SIGMA = LoadDataset(['derivatives/EEG-segmented/sub-', sub, '/ses-', ses, '/sub-', sub ,'_ses-', ses, '_task-psg_desc-sigmanrembout_pow.set'], 'all');

% -------------------------------------------------------------------------
% Extract hypnogram information
Thyp = css_eeglab2hypnogram(EEG);
% -------------------------------------------------------------------------
% Select EEG during N2 and N3
nrembouts = getnrembouts(Thyp, EEG.srate, 300); % '300' is minimum bout duration in seconds
NREM = pop_select(EEG, 'time', nrembouts);
NREM = eeg_checkset(NREM, 'eventconsistency');
NREM.times = NREM.xmin:1/NREM.srate:NREM.xmax; % undo 'pop_select' convertion to milliseconds

%%

chan1 = 78; % good one
chan2 = 72; % bad one

i = 0;

Fig = figure; 
Fig.Position = [1, 660, 1920, 450];

i = i+1;
Ax(i) = axes();
Ax(i).NextPlot = 'add';
plot(Ax(i), NREM.times, NREM.data(chan1,:)); 
plot(Ax(i), NREM.times, NREM.data(chan2,:)+150);
Ax(i).Color = [0.97, 0.98, 0.99];
Ax(i).XLabel.String = 'time (s)';
Ax(i).YLabel.String = 'EEG';
Ax(i).XTick = 0:60:NREM.xmax;
Ax(i).XLim = [0, NREM.xmax];
Ax(i).YLim = [-150, 300];
Ax(i).TickLength = [0, 0];
Ax(i).XGrid = 'on';
Ax(i).YTick = [0, 100];
Ax(i).OuterPosition = [0, 0, 0.9, 0.5];

i = i+1;
Ax(i) = axes();
Ax(i).NextPlot = 'add';
plot(Ax(i), SIGMA.times, SIGMA.data(chan1,:)); 
plot(Ax(i), SIGMA.times, SIGMA.data(chan2,:)+100);
Ax(i).Color = [0.97, 0.98, 0.99];
Ax(i).XLabel.String = 'time (s)';
Ax(i).YLabel.String = 'Sigma power';
Ax(i).XTick = 0:60:SIGMA.xmax;
Ax(i).XLim = [0, SIGMA.xmax];
Ax(i).YLim = [0, 200];
Ax(i).TickLength = [0, 0];
Ax(i).XGrid = 'on';
Ax(i).YTick = [0, 100];
Ax(i).OuterPosition = [0, 0.5, 0.9, 0.5];

i = i+1;
Ax(i) = axes();
Ax(i).NextPlot = 'add';
plot(Ax(i), ISF.freqs, ISF.data(chan1,:))
plot(Ax(i), ISF.freqs, ISF.data(chan2,:))
Ax(i).Color = [0.97, 0.98, 0.99];
Ax(i).XLabel.String = 'freq (Hz)';
Ax(i).YLabel.String = 'ISF power';
Ax(i).XTick = 0:0.01:0.1;
Ax(i).YLim = [0, 8];
Ax(i).TickLength = [0, 0];
Ax(i).XGrid = 'on';
Ax(i).YTick = 0:12;
Ax(i).OuterPosition = [0.85, 0.5, 0.15, .5];

linkaxes(Ax(1:2), 'x')


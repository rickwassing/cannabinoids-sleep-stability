function D = load_isf_topography_data(type)
% -------------------------------------------------------------------------
% Load everything needed to draw the ISF topography figure (Figure 2,
% panels A-H): channel locations, group-level GLM results, colormaps, bout
% descriptives, an example EEG/sigma/HR recording, and every subject's ISF
% and cross-correlation first-level outcomes.
% -------------------------------------------------------------------------
D = struct();
D.chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({D.chanlocs.labels});
incl = hdeeg_scalpchannels(system);
D.chanlocs = D.chanlocs(ismember({D.chanlocs.labels}, incl));
D.chanlocs = channel_clusters(D.chanlocs, 'mff');
% -------------------------------------------------------------------------
% Load group-level results, colormap, and chanlocs
grp = load(sprintf('group-level/a1c_pairttest_cmass/%ssigma/glm.mat', type));
D.GLM = grp.GLM;
grpx = load('group-level/a1c_pairttest_cmass/xcorr/glm.mat');
D.GLMXC = grpx.GLM;
cmap_roma = load('colormap_roma.mat');
D.roma = cmap_roma.roma;
cmap_batlow = load('colormap_batlow.mat');
D.batlow = cmap_batlow.batlow;
D.BOUTS.a = readtable('group-level/nrembout_number.csv');
D.BOUTS.b = readtable('group-level/nrembout_duration.csv');
% -------------------------------------------------------------------------
% Set limits
switch type
    case 'abs'
        D.AmpYLim = [0, 2];
    case 'norm'
        D.AmpYLim = [0, 4];
end
% -------------------------------------------------------------------------
% Load example dataset
D.SIGMA = LoadDataset('derivatives/EEG-segmented/sub-r005/ses-etc120/sub-r005_ses-etc120_task-psg_desc-sigmanrembout_pow.set', 'all');
D.EEG = LoadDataset('derivatives/EEG-processed/sub-r005/ses-etc120/sub-r005_ses-etc120_task-psg_desc-sigma_pow.set', 'all');
D.HR = LoadDataset('derivatives/EEG-preproc/sub-r005/ses-etc120/sub-r005_ses-etc120_task-psg_desc-preprochr_hr.set', 'all');
D.bouts = getnrembouts(css_eeglab2hypnogram(D.EEG), D.EEG.srate, 300);
D.EEG = pop_select(D.EEG, 'time', [0, D.bouts(end, 2)+300]);
D.EEG.times = linspace(D.EEG.xmin, D.EEG.xmax, D.EEG.pnts);
D.HR = pop_select(D.HR, 'time', [0, D.bouts(end, 2)+300]);
D.HR.times = linspace(D.HR.xmin, D.HR.xmax, D.HR.pnts);
% -------------------------------------------------------------------------
% Load data from both groups
Files = dir(sprintf('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*%ssigma*interp_fstlvl.mat', type));
D.ISF = [];
for i = 1:length(Files)
    if i == 1
        D.ISF = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    else
        D.ISF(i) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix'); %#ok<AGROW>
    end
end
% -------------------------------------------------------------------------
% Cross correlation between Sigma power and HR timeseries
Files = dir('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*_desc-nremboutxcorr120s_fstlvl.mat');
D.XC = [];
for i = 1:length(Files)
    tmp = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    if isempty(D.XC)
        D.XC = tmp;
    else
        D.XC(i) = tmp; %#ok<AGROW>
    end
end

end

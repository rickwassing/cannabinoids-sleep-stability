% function css_analyse_1b()
% -------------------------------------------------------------------------
% STEP 1:
% Use 300-second bouts of continuous N2 sleep to describe the sigma/spindle
% ISF between conditions:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SUB: b
% - Confirm that sigma power/spindle density is similar between conditions
force = false;
SigmaFiles = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*_desc-sigmanrembout*.set');
SpdFiles = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*_desc-spindlefdznrembout*.set');
TBL = table();
for f = 1:length(Files)
    hdr = LoadDataset(fullfile(SigmaFiles(f).folder, SigmaFiles(f).name), 'header');
    kv = filename2struct(hdr.setname);
    outfname = sprintf('sub-%s_ses-%s_task-psg_desc-a1b_fstlvl.mat', kv.sub, kv.ses);
    if ~force && exist(sprintf('derivatives/EEG-output-fstlvl/sub-%s/ses-%s/%s', kv.sub, kv.ses, outfname), 'file') == 2
        fprintf('Output exists: skipping %s.\n', outfname)
        continue
    end
    % Load sigma power
    EEG = LoadDataset(fullfile(SigmaFiles(f).folder, SigmaFiles(f).name), 'all');
    Features = struct();
    Features(1).label = 'sigma';
    Features(1).type = 'power';
    Features(1).data = double(mean(EEG.data, 2));
    % Load Spindle events
    EEG = LoadDataset(fullfile(SigmaFiles(f).folder, SpdFiles(f).name), 'all'); % using 'SigmaFiles(f).folder' to make sure its the same subject and session
    SpdDens = zeros(EEG.nbchan, 1);
    RecDurMin = EEG.pnts / (EEG.srate*60);
    for c = 1:EEG.nbchan
        SpdDens(c) = length(find(diff([0, EEG.data(c,:) > 0.5] == 1))) / RecDurMin;
    end
    Features(2).label = 'spindle';
    Features(2).type = 'density';
    Features(2).data = SpdDens;
    % Save output
    css_createfstlvloutput(outfname, Features);
end
% -------------------------------------------------------------------------
% Create supplementary figure showing mean sigma power and spindle density
% in each condition
% -------------------------------------------------------------------------
% Init
close all
Fig = figure;
Fig.Units = 'centimeters';
Fig.Position(3:4) = [9, 9];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CMapRoma = load('colormap_roma.mat');
CMapRoma = CMapRoma.roma;
CMapBatlow = load('colormap_batlow.mat');
CMapBatlow = CMapBatlow.batlow;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ChanLocs = readlocs('group-level/a1b_pairttest_cmass/input/chanlocs.sfp');
DesMat = readmatrix('group-level/a1b_pairttest_cmass/input/desmat.csv');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clear Ax
ai = 0;
% -------------------------------------------------------------------------
% Sigma power
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% p-values
TData = readmatrix('group-level/a1b_pairttest_cmass/output/palm_tstat1_dpv_tstat_m1_d1_c1.csv');
PData = readmatrix('group-level/a1b_pairttest_cmass/output/palm_tstat1_clusterm_tstat_fwep_m1_d1_c1.csv');
fprintf('##############################################################\n')
fprintf('# REPORT\n')
fprintf('No significant differences in sigma power (PALM, cluster-mass FWE-corrected p > %.2f).\n', min(PData))
fprintf('##############################################################\n')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ai = ai+1;
Ax(ai) = axes(Fig);
Ax(ai).OuterPosition = [-0.08 0.5 0.5 0.5];
Ax(ai).Title.String = 'PLACEBO';
Ax(ai).Title.FontSize = 10;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CData = readmatrix('group-level/a1b_pairttest_cmass/input/depvar_1.csv');
CData = mean(CData(DesMat(:, 1) == -1, :), 1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
topoplot(CData, ChanLocs, ...
    'emarker2', {find(PData < 0.05), '.', [1 1 1], 8}, ...
    'contourvals', (TData > -1.9623) + (TData > 1.9623), ...
    'numcontour', length(unique((TData > -1.9623) + (TData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ai = ai+1;
Ax(ai) = axes(Fig);
Ax(ai).OuterPosition = [0.34 0.5 0.5 0.5];
Ax(ai).Title.String = 'ETC120';
Ax(ai).Title.FontSize = 10;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CData = readmatrix('group-level/a1b_pairttest_cmass/input/depvar_1.csv');
CData = mean(CData(DesMat(:, 1) == 1, :), 1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
topoplot(CData, ChanLocs, ...
    'emarker2', {find(PData < 0.05), '.', [1 1 1], 8}, ...
    'contourvals', (TData > -1.9623) + (TData > 1.9623), ...
    'numcontour', length(unique((TData > -1.9623) + (TData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% -------------------------------------------------------------------------
% Spindle density
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% p-values
TData = readmatrix('group-level/a1b_pairttest_cmass/output/palm_tstat1_dpv_tstat_m2_d1_c1.csv');
PData = readmatrix('group-level/a1b_pairttest_cmass/output/palm_tstat1_clusterm_tstat_fwep_m2_d1_c1.csv');
fprintf('##############################################################\n')
fprintf('# REPORT\n')
fprintf('No significant differences in spindle density (PALM, cluster-mass FWE-corrected p > %.2f).\n', min(PData))
fprintf('##############################################################\n')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ai = ai+1;
Ax(ai) = axes(Fig);
Ax(ai).OuterPosition = [-0.08 0 0.5 0.5];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CData = readmatrix('group-level/a1b_pairttest_cmass/input/depvar_2.csv');
CData = mean(CData(DesMat(:, 1) == -1, :), 1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
topoplot(CData, ChanLocs, ...
    'emarker2', {find(PData < 0.05), '.', [1 1 1], 8}, ...
    'contourvals', (TData > -1.9623) + (TData > 1.9623), ...
    'numcontour', length(unique((TData > -1.9623) + (TData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ai = ai+1;
Ax(ai) = axes(Fig);
Ax(ai).OuterPosition = [0.34 0 0.5 0.5];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CData = readmatrix('group-level/a1b_pairttest_cmass/input/depvar_2.csv');
CData = mean(CData(DesMat(:, 1) == 1, :), 1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
topoplot(CData, ChanLocs, ...
    'emarker2', {find(PData < 0.05), '.', [1 1 1], 8}, ...
    'contourvals', (TData > -1.9623) + (TData > 1.9623), ...
    'numcontour', length(unique((TData > -1.9623) + (TData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% -------------------------------------------------------------------------
Ax(1).CLim = [8, 11];
Ax(1).Colormap = CMapBatlow;
Ax(2).CLim = [8, 11];
Ax(2).Colormap = CMapBatlow;
Ax(3).CLim = [8, 13];
Ax(3).Colormap = CMapBatlow;
Ax(4).CLim = [8, 13];
Ax(4).Colormap = CMapBatlow;
% -------------------------------------------------------------------------
clear CBar
CBar(1) = colorbar(Ax(1), 'west');
CBar(1).Position = [0.85, 0.6, 0.03, 0.3];
CBar(1).Label.String = 'SIGMA POWER';
CBar(1).Ticks = Ax(1).CLim;
CBar(2) = colorbar(Ax(3), 'west');
CBar(2).Position = [0.85, 0.1, 0.03, 0.3];
CBar(2).Label.String = 'SPINDLE DENSITY (N/min)';
CBar(2).Ticks = Ax(3).CLim;

% -------------------------------------------------------------------------
% SAVE FIGURE
disp('Saving figure...')
exportgraphics(Fig, 'figures/supp_1.png', 'Resolution', 1200)
disp('Done saving')
% end
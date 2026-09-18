function plot_selected_parietal_channels(chanlocs, idx_chan)
% -------------------------------------------------------------------------
% Create supplementary figure showing which parietal channels were
% selected, coloured by their mean ISF amplitude (feature 2).
% -------------------------------------------------------------------------
type = 'norm';
Files = dir(sprintf('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*%ssigma*interp_fstlvl.mat', type));
ISF = [];
for i = 1:length(Files)
    if i == 1
        ISF = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    else
        ISF(i) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix'); %#ok<SAGROW>
    end
end
YData = arrayfun(@(s) s.features(2).data, ISF, 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Fig = figure();
Fig.Units = 'centimeters';
Fig.Position(3:4) = [8 8];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot
load('colormap_batlow.mat')
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [0, 2], ...
    'numcontour', 0, ...
    'conv', 'on', ...
    'whitebk', 'on', ...
    'emarker2', {idx_chan, '.', 'w', 15, 1});
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exportgraphics(Fig, 'figures/supp_selparietalchans.png', 'Resolution', 600)

end

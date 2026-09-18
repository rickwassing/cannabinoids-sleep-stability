function css_analyse_filter_edge_artefact()
clc
clear
close all
% -------------------------------------------------------------------------
% Determine the impact of the filter edge artefact on the phase estimate of
% the ISF
% -------------------------------------------------------------------------
% Define the set of parietal channels to focus on
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
chanlocs = channel_clusters(chanlocs, 'mff');
idx_chan = [77, 78, 79, 85, 86, 87, 92, 93, 94, 99, 104, 110, 111, 112, 113, 120, 121, 122];
chanlabels = {chanlocs(idx_chan).labels}; %#ok<NASGU>
append_type = 'zeros';
ford_mult = 2;
% -------------------------------------------------------------------------
% Create supplementary figure showing which channels were selected
plot_selected_parietal_channels(chanlocs, idx_chan);
% -------------------------------------------------------------------------
% Load the SIGMA power timeseries
Files = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*_desc-sigmanrembout*.set');
clear SIGMA
for i = 1:length(Files)
    SIGMA(i) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'all'); %#ok<SAGROW> 
    SIGMA(i) = pop_select(SIGMA(i), 'channel', idx_chan); %#ok<SAGROW> % Select parietal channels
end
% -------------------------------------------------------------------------
% Get emperical values for the filter cutoffs
clear filtcfg
filtcfg.cutoff = [0.0138 0.0259];
filtcfg.wintype = 'kaiser';
filtcfg.transbw = 1/50;
filtcfg.rippledev = 0.05;
filtcfg.warg = 6;
filtcfg.adj = [0 0];
filtcfg.order = pop_firwsord(filtcfg.wintype, SIGMA(1).srate, filtcfg.transbw, filtcfg.rippledev); % transition bandwidth of 0.01 Hz
% -------------------------------------------------------------------------
% Plot the filter response
plot_filter_response_supp(SIGMA, filtcfg, append_type);
% -------------------------------------------------------------------------
% Compute the control/test signals and their phase/amplitude difference
[SIGA, SIGB, ANGA, ANGB, DANG, DAMP, sig_control, sig_test] = compute_filter_edge_artefact(SIGMA, filtcfg, append_type, ford_mult);
%%
save('analysis_2a.mat', 'SIGA', 'SIGB', 'ANGA', 'ANGB', 'DANG', 'DAMP', 'filtcfg', 'SIGMA', 'idx_chan', 'sig_test', 'sig_control', 'append_type', 'ford_mult');
%%
N%% -------------------------------------------------------------------------
% Create supplementary figure
plot_filter_edge_artefact_supp(SIGMA, filtcfg, sig_test, sig_control, SIGA, SIGB, ANGA, ANGB, DANG, ford_mult, append_type);

end

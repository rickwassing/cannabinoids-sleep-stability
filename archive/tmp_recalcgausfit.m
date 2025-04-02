clc
sigmaf = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*_ses-*_desc-sigmanrembout_eeg.set'); 

for i = 1:length(sigmaf)
    try
        SIGMA = LoadDataset(fullfile(sigmaf(i).folder, sigmaf(i).name), 'all');
        % Keep only EEG chans
        SIGMA = pop_select(SIGMA, 'channel', find(strcmpi({SIGMA.chanlocs.type}, 'eeg')));
        % Interpolate outliers
        SIGMA = css_sigmainterp(SIGMA, 24);
        plotallsigmachannels(SIGMA)
        close all
        FIT = infraslowmodpowerspect(SIGMA, true);
        kv = filename2struct(SIGMA.setname);
        kv.desc = 'ismfit';
        kv.filetype = 'struct.mat';
        save(fullfile(SIGMA.filepath, struct2filename(kv)), 'FIT', '-v7.3');
    catch ME
        disp(getReport(ME))
        break
    end
end
disp('done')
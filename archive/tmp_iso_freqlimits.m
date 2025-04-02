Files = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*_desc-ismfit_struct.mat');


ISO = struct();
ISO.mu = [];
ISO.ll = [];
ISO.ul = [];
rt = now();
for i = 1:length(Files)
    rt = remainingTime(rt, length(Files));
    load(fullfile(Files(i).folder, Files(i).name));
    ISO.mu = [ISO.mu; FIT(1:end-1).mu];
    ISO.ll = [ISO.ll; FIT(1:end-1).lower];
    ISO.ul = [ISO.ul; FIT(1:end-1).upper];
end

%%
chanlocs = readlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
chanlocs = channel_clusters(chanlocs, 'MFF');

%%
% Frontal = 2
% Central = 3
% Parietal = 4
% Occipital = 5
% Temporal right = 1
% Temporal left = 6

idx_f = [chanlocs.cluster] == 2;
idx_to = [chanlocs.cluster] == 5 | [chanlocs.cluster] == 1 | [chanlocs.cluster] == 6;
idx_cp = [chanlocs.cluster] == 3 | [chanlocs.cluster] == 4;

mu_f = ISO.mu(:, idx_f);
mu_to = ISO.mu(:, idx_to);
mu_cp = ISO.mu(:, idx_cp);

%%
disp('Frontal')
fprintf('mu = %.3f [95%% of estimates between %.3f - %.3f]\n', mean(mu_f(:)), prctile(mu_f(:), 2.5), prctile(mu_f(:), 97.5))

disp('Temporooccipital')
fprintf('mu = %.3f [95%% of estimates between %.3f - %.3f]\n', mean(mu_to(:)), prctile(mu_to(:), 2.5), prctile(mu_to(:), 97.5))

disp('Centroparietal')
fprintf('mu = %.3f [95%% of estimates between %.3f - %.3f]\n', mean(mu_cp(:)), prctile(mu_cp(:), 2.5), prctile(mu_cp(:), 97.5))

%%

close all

figure
histogram(ISO.mu(:))
%%
mean(ISO.mu(:))
%%
prctile(ISO.mu(:), 2.5)
prctile(ISO.mu(:), 97.5)
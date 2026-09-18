Files = dir('derivatives/EEG-preproc/sub-*/ses*/sub-*_desc-ismfit_struct.mat');

BW = struct();
for i = 1:length(Files)
    load(fullfile(Files(i).folder, Files(i).name))
    BW(i).name = Files(i).name;
    BW(i).lower = [FIT.lower];
    BW(i).upper = [FIT.upper];
end

%%

close all

histogram([BW.lower])
prctile([BW.lower], 1)

figure

histogram([BW.upper])
prctile([BW.lower], 99)
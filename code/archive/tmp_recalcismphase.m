clear
clc
eegf = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*_ses-*preproc_eeg.set');
nremf = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*_ses-*_desc-nrem_eeg.set');
sigmaf = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*_ses-*_desc-sigma_eeg.set'); 
fitf = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*_ses-*_desc-ismfit_struct.mat'); 

for i = 9:length(eegf)
    try
        EEG = LoadDataset(fullfile(eegf(i).folder, eegf(i).name), 'header');
        NREM = LoadDataset(fullfile(nremf(i).folder, nremf(i).name), 'header');
        SIGMA = LoadDataset(fullfile(sigmaf(i).folder, sigmaf(i).name), 'all');
        FIT = load(fullfile(fitf(i).folder, fitf(i).name));
        FIT = FIT.FIT;
        
        [EEG, NREM, SIGMA] = extractismphase(EEG, NREM, SIGMA, FIT);
        
        EEG = SaveDataset(EEG, 'header');
        NREM = SaveDataset(NREM, 'header');
        SIGMA = SaveDataset(SIGMA, 'header');
    catch ME
        disp(getReport(ME))
        break
    end
end
disp('done or error')
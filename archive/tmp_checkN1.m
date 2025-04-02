Files = dir('/Volumes/sleep/Sleep/3. ACTIVE STUDIES/CUPID/Arousal paper backup_CANSLEEP/derivatives/EEG-preproc/sub-*/ses-*/sub-*_hypno.csv');

n1boutlength = [];
n1wboutlength = [];
for i = 1:length(Files)
    t = readtable(fullfile(Files(i).folder, Files(i).name));

    tmp = diff([find(diff([1;t.sleepstage;1] == -1) == 1), find(diff([1;t.sleepstage;1] == -1) == -1)]');
    n1boutlength = [n1boutlength, tmp];

    tmp = diff([find(diff([9;t.sleepstage;9] == -1 | [9;t.sleepstage;9] == 1) == 1), find(diff([9;t.sleepstage;9] == -1 | [9;t.sleepstage;9] == 1) == -1)]');
    n1wboutlength = [n1wboutlength, tmp];
end

disp('done')

%%

close all
histogram(n1wboutlength, 'BinMethod', 'integers')

%%
close all
hold on

stairs(t.sleepstage)
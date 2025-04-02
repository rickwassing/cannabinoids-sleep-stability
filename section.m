function [errors] = section(Proc, Files, cfg)
% -------------------------------------------------------------------------
% Init output variable
errors = {};
% -------------------------------------------------------------------------
% Check input
def.do_parallel = false;
def.hosts = {};
def.this_host = '';
def.force = false;
fnames = fieldnames(cfg);
for i = 1:length(fnames)
    def.(fnames{i}) = cfg.(fnames{i});
end
cfg = def;
% -------------------------------------------------------------------------
% Read the phenotype table
PHEN = readtable('Scoring log and notes_CANSLEEP arousal.xlsx');
% -------------------------------------------------------------------------
% Set the output filepaths
for i = 1:length(Files)
    % Get key-values from filename
    kv = filename2struct(Files(i).name);
    % Check if we need to replace the session value
    if ismember(kv.ses, {'1', '2'})
        % replace with condition
        idx_phen = strcmpi(PHEN.folder_name, kv.sub);
        kv.ses = lower(PHEN.condition{idx_phen});
        % remove the last 3 digits from the subject id
        kv.sub = kv.sub(1:end-3);
    end
    % We don't need the 'run' key-value, remove
    if isfield(kv, 'run')
        kv = rmfield(kv, 'run');
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Set the output folder name
    if strcmpi(cfg.derivativeout, 'n/a')
        Files(i).outfolder = 'n/a';
    else
        Files(i).outfolder = sprintf('%s/derivatives/%s/sub-%s/ses-%s', pwd, cfg.derivativeout, kv.sub, kv.ses);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Set the output filenames
    for j = 1:length(Proc)
        if strcmpi(cfg.derivativeout, 'n/a')
            Files(i).outname{j} = 'imsureidontexist.txt';
        else
            % Set the 'desc' and 'filetype' value
            this_kv = kv;
            this_kv.desc = strrep(cfg.kv.desc{j}, '<desc>', kv.desc);
            this_kv.filetype = strrep(cfg.kv.filetype{j}, '<filetype>', kv.filetype);
            if contains(this_kv.filetype, 'fstlvl')
                if ~contains(this_kv.filetype, '.mat')
                    this_kv.filetype = [this_kv.filetype, '.mat'];
                end
            else
                if ~contains(this_kv.filetype, '.set')
                    this_kv.filetype = [this_kv.filetype, '.set'];
                end
            end
            Files(i).outname{j} = struct2filename(this_kv);
        end
    end
end
% -------------------------------------------------------------------------
% Check that output file does not exist yet
for i = 1:length(Files)
    for j = 1:length(Proc)
        if ~cfg.force && exist(fullfile(Files(i).outfolder, Files(i).outname{j}), 'file') == 2
            fprintf('File ''%s'' exists.\n', fullfile(Files(i).outfolder, Files(i).outname{j}))
            Files(i).outname{j} = 'exists';
        end
    end
end
% Remove files of which all output files exist
idx_rm = arrayfun(@(f) all(strcmpi(f.outname, 'exists')), Files);
Files(idx_rm) = [];
% -------------------------------------------------------------------------
% Assign a host to each file
if cfg.do_parallel
    for i = 1:length(Files)
        if all(strcmpi(Files(i).outname, 'exists'))
            Files(i).host = 'none';
        else
            m = mod(i, length(cfg.hosts))+1;
            Files(i).host = cfg.hosts{m};
        end
    end
end
% -------------------------------------------------------------------------
% Print to screen
if cfg.do_parallel
    fprintf('=========================================================================\n')
    fprintf('RUNNING PARALLEL PROCESSES ON %i HOSTS\n', length(cfg.hosts))
else
    fprintf('=========================================================================\n')
    fprintf('RUNNING PROCESSES ON THIS COMPUTER ONLY\n')
end
fprintf('-------------------------------------------------------------------------\n')
fprintf('Running the following processes:\n')
for i = 1:length(Proc)
    fprintf('> %s\n', Proc{i})
end
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
fprintf('On %i files:\n', length(Files))
for i = 1:length(Files)
    if cfg.do_parallel
        hstr = Files(i).host;
    else
        hstr = 'this computer';
    end
    fprintf('> %s in %s on %s.\n', Files(i).name, Files(i).folder, hstr)
end
fprintf('-------------------------------------------------------------------------\n')
% -------------------------------------------------------------------------
% Run across files
rt = now();
for i = 1:length(Files)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fprintf('\n')
    fprintf('#########################################################################\n')
    fprintf('File %i of %i, filename: ''%s''.\n', i, length(Files), Files(i).name)
    rt = remainingTime(rt, length(Files), true);
    fprintf('#########################################################################\n')
    fprintf('\n')
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % If we do parallel processing, then check if this file is assigned to
    % this computer, if not continue to the next file
    if cfg.do_parallel
        if ~strcmpi(cfg.this_host, Files(i).host)
            fprintf('Skipping file ''%s'', assigned to host ''%s''.\n', Files(i).name, Files(i).host)
            continue
        end
    end
    if ~cfg.force && all(strcmpi(Files(i).outname, 'exists'))
        fprintf('Skipping file ''%s'', all output files exist.\n', Files(i).name)
        continue
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Run
    try
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Load the input data
        EEG = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'all');
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Apply functions
        for j = 1:length(Proc)
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Check outout file
            if ~cfg.force && strcmpi(Files(i).outname{j}, 'exists')
                fprintf('Skipping process ''%s'' on file ''%s'', output file exists.\n', Proc{j}, Files(i).name)
                continue
            end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Check what process to run
            switch Proc{j}
                case 'preproc'
                    param = struct();
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_preproc(EEG, param);
                case 'detectspindlesusingfernandez'
                    param = struct();
                    param.method = 'fernandez';
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_extractspindles(EEG, param);
                case 'detectspindlesusingferrarelli'
                    param = struct();
                    param.method = 'ferrarelli';
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_extractspindles(EEG, param);
                case 'detectspindlesusingwamsley'
                    param = struct();
                    param.method = 'wamsley';
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_extractspindles(EEG, param);
                case 'getsigmapowerusingwavelet'
                    param = struct();
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_getsigmapowerusingwavelet(EEG, param);
                case 'getdeltapowerusingwavelet'
                    param = struct();
                    param.MinFreq = 1;
                    param.MaxFreq = 4; 
                    param.Cycles = 3;
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_getspecpowerusingwavelet(EEG, param);
                case 'getthetapowerusingwavelet'
                    param = struct();
                    param.MinFreq = 4;
                    param.MaxFreq = 8; 
                    param.Cycles = 6;
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_getspecpowerusingwavelet(EEG, param);
                case 'inspectspindles'
                    param = struct();
                    css_inspectspindles(EEG, param);
                case 'extractnrembouts'
                    param = struct();
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_extractnrembouts(EEG, param);
                case 'extractnremarousals'
                    param = struct();
                    param.stage = -2;
                    param.cutoff = 100;
                    param.allowalpha = false;
                    param.allowmicro = false;
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_extractarousalbouts(EEG, param);
                case 'extractprerembouts'
                    param = struct();
                    param.cutoff = 300;
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_extractprerembouts(EEG, param);
                case 'infraslowfluctpowerspect_norm'
                    param.doPlot = true;
                    param.doFilter = true;
                    param.normalize = true;
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_infraslowfluctpowerspect(EEG, param);
                case 'infraslowfluctpowerspect_abs'
                    param.doPlot = true;
                    param.doFilter = true;
                    param.normalize = false;
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_infraslowfluctpowerspect(EEG, param);
                case 'crosscorr'
                    param.doPlot = true;
                    param.outfilepath = fullfile(Files(i).outfolder, Files(i).outname{j});
                    css_crosscorr(EEG, param);
                otherwise
                    error('Do not know what process ''%s'' is.', Proc{j})
            end
        end
    catch ME
        % -----------------------------------------------------------------
        % Something has gone wrong
        disp('****************************************************')
        disp(getReport(ME))
        disp('****************************************************')
        err = struct();
        err.it = i;
        err.File = Files(i);
        err.ME = ME;
        errors = [errors; {err}]; %#ok<AGROW>
    end
end
% -------------------------------------------------------------------------
% Save errors
if cfg.do_parallel && ~isempty(errors)
    save(['errors_', cfg.this_host, '_', char(datetime(), 'yyyyMMdd''T''HHmmss'), '.mat'], 'errors', '-mat', '-v7.3');
end
end
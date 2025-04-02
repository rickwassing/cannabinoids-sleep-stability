% BIDS_WEBSITE
% Creates the website for data auditing purposes and simple file inspection.
%
% Usage:
%   >> bids_website(bidsroot)
%
% Inputs:
%   'bidsroot' - [char] full or relative path to the root of a BIDS directory
%
% Outputs:
%   none

% Authors:
%   Rick Wassing, Woolcock Institute of Medical Research, Sydney, Australia
%
% History:
%   Created 2023-07-20, Rick Wassing

% (C) 2023 by Rick Wassing, licensed under
% Attribution-NonCommercial-ShareAlike 4.0 International
% This license requires that reusers give credit to the creator. It allows
% reusers to distribute, remix, adapt, and build upon the material in any
% medium or format, for noncommercial purposes only. If others modify or
% adapt the material, they must license the modified material under
% identical terms.

function bids_website(bidsroot)
knownmodalities = {'beh', 'eeg', 'fnirs'};
knownextensions = {'edf', 'set', 'snirf', 'txt', 'csv', 'tsv'};
% =========================================================================
% INIT
% -------------------------------------------------------------------------
cd(bidsroot)
DatasetInformation = jsondecode(fileread('./rawdata/dataset_description.json'));
Participants = ft_read_tsv('./rawdata/participants.tsv');
[~, idx_sort] = sort(Participants.participant_id);
Participants = Participants(idx_sort, :);
% =========================================================================
% CREATE INDEX PAGE
% -------------------------------------------------------------------------
props = struct();
props.bidsroot = bidsroot;
props.title = 'NeuroVosa';
props.filename = 'index.html';
props.breadcrumb.li.label = 'NeuroVosa';
props.col(1).title = 'Dataset';
props.col(1).component = 'StructList';
props.col(1).props.info = DatasetInformation;
props.col(2).title = 'Raw data';
props.col(2).component = 'RawDataList';
props.col(2).props.participants = Participants;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fprintf('>> ======================================================\n');
fprintf('>> WEB: Creating root index.html\n');
PageTwoCols(props);
% =========================================================================
% FOR EACH SUBJECT CREATE A SUBJECT INDEX PAGE
% -------------------------------------------------------------------------
for i = 1:size(Participants, 1)
    % Read all scan files for this subject
    scansfiles = dir(['./rawdata/', Participants.participant_id{i}, '/**/sub-*_scans.tsv']);
    % Set properties
    props = struct();
    props.bidsroot = strrep(bidsroot, '\', '/');
    props.title = Participants.participant_id{i};
    props.filename = ['./rawdata/', Participants.participant_id{i}, '.html'];
    props.breadcrumb.li(1).href = '../index.html';
    props.breadcrumb.li(1).label = 'NeuroVosa';
    props.breadcrumb.li(2).label = Participants.participant_id{i};
    props.col(1).title = Participants.participant_id{i};
    props.col(1).component = 'StructList';
    props.col(1).props.info = struct('age', Participants.age(i), 'sex', Participants.sex{i});
    props.col(2).title = 'Scans';
    props.col(2).component = 'ScansList';
    for j = 1:length(scansfiles)
        props.col(2).props(j).title = strrep(strrep(scansfiles(j).folder, bidsroot, '.'), '\', '/');
        props.col(2).props(j).scans = ft_read_tsv(fullfile(scansfiles(j).folder, scansfiles(j).name));
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fprintf('>> ======================================================\n');
    fprintf('>> WEB: Creating page for subject ''%s''\n', props.title);
    PageTwoCols(props);
end

% =========================================================================
% FOR EACH SCAN IN EACH SUBJECT CREATE A SCAN QC PAGE
% -------------------------------------------------------------------------
for i = 1:size(Participants, 1)
    % Read all scan files for this subject
    scansfiles = dir(['./rawdata/', Participants.participant_id{i}, '/**/sub-*_scans.tsv']);
    for j = 1:length(scansfiles)
        % Read scans
        scans = ft_read_tsv(fullfile(scansfiles(j).folder, scansfiles(j).name));
        folder = strrep(scansfiles(j).folder, bidsroot, '.');
        for k = 1:size(scans, 1)
            [~, filename] = fileparts(scans.filename{k});
            filename = strrep(filename, '\', '/');
            % Set properties
            props = struct();
            props.bidsroot = bidsroot;
            props.title = strrep(scans.filename{k}, '\', '/');
            props.filename = [folder, filesep, filename, '.html'];
            props.breadcrumb.li(1).href = strjoin([repmat({'../'}, 1, length(strsplit(folder, filesep))-1), {'index.html'}], '');
            props.breadcrumb.li(1).label = 'NeuroVosa';
            props.breadcrumb.li(2).href = strjoin([repmat({'../'}, 1, length(strsplit(folder, filesep))-2), {[Participants.participant_id{i}, '.html']}], '');
            props.breadcrumb.li(2).label = Participants.participant_id{i};
            props.breadcrumb.li(3).label = props.title;
            
            if iscell(scans.acq_time)
                props.subtitle = iso2human(scans.acq_time{k}, 'omitmilliseconds', true);
            elseif isdatetime(scans.acq_time)
                props.subtitle = iso2human(scans.acq_time(k), 'omitmilliseconds', true);
            else
                props.subtitle = 'n/a';
            end

            r = 0;

            methodsfile = dir(sprintf('%s/%s.txt', folder, strrep(filename, '_nirs', '_methods')));
            if length(methodsfile) == 1
                r = r+1;
                props.row(r).component = 'EmbedTxt';
                props.row(r).props.header = 'Methods';
                props.row(r).props.file = fullfile(methodsfile.folder, methodsfile.name);
            end

            citationsfile = dir(sprintf('%s/%s.txt', folder, strrep(filename, '_nirs', '_cites')));
            if length(citationsfile) == 1
                r = r+1;
                props.row(r).component = 'EmbedTxt';
                props.row(r).props.header = 'References (DOI)';
                props.row(r).props.file = fullfile(citationsfile.folder, citationsfile.name);
            end

            r = r+1;
            props.row(r).component = 'ImageList';
            props.row(r).props.title = 'Channel quality';
            props.row(r).props.class = 'col-12';
            props.row(r).props.imagefiles = filterfiles(dir(fullfile(folder, 'qc', strrep(filename, '_nirs', '*.png'))), 'contains', '_desc-qc');

            r = r+1;
            props.row(r).component = 'ImageList';
            props.row(r).props.title = 'Nuisanace regression';
            props.row(r).props.class = 'col-12';
            props.row(r).props.imagefiles = filterfiles(dir(fullfile(folder, 'qc', strrep(filename, '_nirs', '*.png'))), 'contains', '_desc-glmtimeseries');
            
            r = r+1;
            props.row(r).component = 'ImageList';
            props.row(r).props.title = 'Average trials (No nuisance regression)';
            props.row(r).props.class = 'd-inline';
            props.row(r).props.imagefiles = filterfiles(dir(fullfile(folder, 'qc', strrep(filename, '_nirs', '*.png'))), 'contains', {'desc-trialsacrosschans', 'source-dc', '_event'});
            
            r = r+1;
            props.row(r).component = 'ImageList';
            props.row(r).props.title = 'Average trials (After nuisance regression)';
            props.row(r).props.class = 'd-inline';
            props.row(r).props.imagefiles = filterfiles(dir(fullfile(folder, 'qc', strrep(filename, '_nirs', '*.png'))), 'contains', {'desc-trialsacrosschans', 'source-glm', '_event'});
            
            r = r+1;
            props.row(r).component = 'ImageList';
            props.row(r).props.title = 'Individual trials (No nuisance regression)';
            props.row(r).props.class = 'd-inline';
            props.row(r).props.imagefiles = filterfiles(dir(fullfile(folder, 'qc', strrep(filename, '_nirs', '*.png'))), 'contains', {'desc-trialswithinchan', 'source-dc', '_event'});

            r = r+1;
            props.row(r).component = 'ImageList';
            props.row(r).props.title = 'Individual trials (After nuisance regression)';
            props.row(r).props.class = 'd-inline';
            props.row(r).props.imagefiles = filterfiles(dir(fullfile(folder, 'qc', strrep(filename, '_nirs', '*.png'))), 'contains', {'desc-trialswithinchan', 'source-glm', '_event'});
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            fprintf('>> ======================================================\n');
            fprintf('>> WEB: Creating page for subject ''%s'' and scan ''%s''\n', Participants.participant_id{i}, props.title);
            PageOneCol(props);
        end
    end
end
% =========================================================================
% Open file
fprintf('>> ======================================================\n');
fprintf('>> WEB: Done, launching website now...\n');
web(fullfile(bidsroot, 'index.html'), '-browser');

end
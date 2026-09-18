function H = extract_prerem_isf_traces(SigmaFiles, chans, filtcfg, evcfg)
% -------------------------------------------------------------------------
% Load all REM onsets and filter the sigma-power ISF, and detect the
% peaks/troughs at and preceding the REM transition, for every subject and
% bout in 'SigmaFiles'.
% -------------------------------------------------------------------------
% Init
H = struct();
cnt = 0;
% -------------------------------------------------------------------------
% Load and process files
for i = 1:length(SigmaFiles)
    % ---------------------------------------------------------------------
    % Load file and select parietal channels only
    SIGMA = LoadDataset(fullfile(SigmaFiles(i).folder, SigmaFiles(i).name), 'all');
    SIGMA = pop_select(SIGMA, 'channel', chans.idx);
    kv = filename2struct(SIGMA.setname);
    % ---------------------------------------------------------------------
    % Normalization factor
    normfact = median(SIGMA.data, 2, 'omitnan');
    % ---------------------------------------------------------------------
    % Extract bout indexes
    bouts = ([SIGMA.event(strcmpi({SIGMA.event.type}, 'boundary')).latency]);
    bouts = [0.5, bouts, SIGMA.pnts+0.5]; %#ok<AGROW>
    % For each bout
    for b = 1:length(bouts)-1
        % -----------------------------------------------------------------
        % Cut and append the data
        sig = pop_select(SIGMA, 'point', [bouts(b), bouts(b+1)]);
        sig.data = 10.*log10(sig.data./repmat(normfact, 1, sig.pnts));
        fsig = sig;
        fsig.data = detrend(fsig.data', 0, 'omitnan')'; % demean
        fsig = signalappend(fsig, 'zeros', 2, filtcfg);
        % -----------------------------------------------------------------
        % Prepend and append data to cover for filter edge artefact
        fsig = executeappending(fsig);
        % -----------------------------------------------------------------
        % Filter the data
        fsig = pop_firws(fsig, ...
            'fcutoff', filtcfg.cutoff, ...
            'ftype', 'bandpass', ...
            'wtype', filtcfg.wintype, ...
            'warg', 2, ...
            'forder', filtcfg.order, ...
            'plotfresp', false, ...
            'minphase', 0);
        % -----------------------------------------------------------------
        % Apply Hilbert
        fsig.cmx = hilbert(detrend(fsig.data', 0))';
        fsig.ang = angle(fsig.cmx);
        fsig.amp = abs(fsig.cmx);
        % -----------------------------------------------------------------
        % Remove appending
        fsig = executeappending(fsig, 'remove');
        % -----------------------------------------------------------------
        % Adjust time vector where zero indicates REM onset
        if b == 1
            remeplat = sig.event(find(strcmpi({sig.event.type}, 'remeps'), 1, 'first')).origlatency;
        else
            thislat = sig.event(find(strcmpi({sig.event.type}, 'remeps'), 1, 'first')).origlatency;
            remeplat = thislat - (remeplat+remepdur);
        end
        tzero = sig.event(find(strcmpi({sig.event.type}, 'remeps'), 1, 'first')).latency;
        remepdur = sig.event(find(strcmpi({sig.event.type}, 'remeps'), 1, 'first')).duration;
        sig.times = ((0:sig.pnts-1)-tzero)./sig.srate;
        % -----------------------------------------------------------------
        % Store
        cnt = cnt+1;
        H(cnt).sub = kv.sub;
        H(cnt).cond = kv.ses;
        H(cnt).rawsigma = sig.data';
        H(cnt).filtsigma = fsig.data';
        H(cnt).times = sig.times;
        H(cnt).cmx = fsig.cmx';
        H(cnt).ang = fsig.ang';
        H(cnt).amp = fsig.amp';
        H(cnt).remepdur = remepdur;
        H(cnt).remeplat = remeplat;
        % -----------------------------------------------------------------
        % Find peaks and troughs at the REM transition, and the ones before that
        H(cnt).event = prerempeakstroughs(H(cnt), evcfg);
    end
end
disp('done')

end

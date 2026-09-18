function html = RawDataList(props)

html = '';
for i = 1:size(props.participants, 1)
    % Read the scans file for this participant
    scansfiles = dir(['./rawdata/', props.participants.participant_id{i}, '/**/sub-*_scans.tsv']);
    numscans = 0;
    for j = 1:length(scansfiles)
        scans = ft_read_tsv(fullfile(scansfiles(j).folder, scansfiles(j).name));
        numscans = numscans + size(scans, 1);
    end
    html = [html, '\n', '<ul class="bg-white border-bottom shadow list-group list-group-horizontal list-group-flush align-items-center">']; %#ok<*AGROW>
    html = [html, '\n', '   <li class="border border-0 list-group-item">'];
    html = [html, '\n', '       <span class="fw-bold">'];
    html = [html, '\n', props.participants.participant_id{i}];
    html = [html, '\n', '       </span>'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '   <li class="border border-0 list-group-item">'];
    html = [html, '\n', '       <small class="text-body-secondary">Age</small>'];
    html = [html, '\n', '       <p class="m-0 fw-bold">'];
    html = [html, '\n', num2str(props.participants.age(i))];
    html = [html, '\n', '       </p>'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '   <li class="border border-0 list-group-item">'];
    html = [html, '\n', '       <small class="text-body-secondary">Sex</small>'];
    html = [html, '\n', '       <p class="m-0 fw-bold">'];
    html = [html, '\n', props.participants.sex{i}];
    html = [html, '\n', '       </p>'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '   <li class="border border-0 list-group-item d-flex flex-grow-1">'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '   <li class="border border-0 list-group-item">'];
    html = [html, '\n', '       <span class="badge bg-secondary rounded-pill">', num2str(numscans), '</span>'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '   <li class="border border-0 list-group-item">'];
    html = [html, '\n', '       <a href="rawdata/', props.participants.participant_id{i},'.html" target="_self" class="btn btn-primary">Go</a>'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '</ul>'];
end
html = sprintf(html);

end
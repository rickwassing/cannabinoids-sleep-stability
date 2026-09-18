function html = ScansList(props)

html = '';
html = [html, '\n', '<p class="lead mt-3">'];
html = [html, '\n', props.title];
html = [html, '\n', '</p>'];
for i = 1:size(props.scans, 1)
    folder = strrep(props.title, './rawdata', '.');
    props.scans.filename{i} = strrep(props.scans.filename{i}, '\', '/');
    [~, filename] = fileparts(props.scans.filename{i});
    html = [html, '\n', '<ul class="bg-white border-bottom shadow list-group list-group-horizontal list-group-flush align-items-center">']; %#ok<*AGROW>
    html = [html, '\n', '   <li class="border border-0 list-group-item ">'];
    html = [html, '\n', '       <span class="fw-bold">'];
    html = [html, '\n', strrep(props.scans.filename{i}, '\', '/')];
    html = [html, '\n', '       </span>'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '   <li class="border border-0 list-group-item">'];
    html = [html, '\n', '       <small class="text-body-secondary">Acquisition time</small>'];
    html = [html, '\n', '       <p class="m-0 fw-bold">'];
    if iscell(props.scans.acq_time)
        html = [html, '\n', iso2human(props.scans.acq_time{i}, 'omitmilliseconds', true)];
    elseif isdatetime(props.scans.acq_time)
        html = [html, '\n', iso2human(props.scans.acq_time(i), 'omitmilliseconds', true)];
    else
        html = [html, '\n', 'n/a'];
    end
    html = [html, '\n', '       </p>'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '   <li class="border border-0 list-group-item d-flex flex-grow-1">'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '   <li class="border border-0 list-group-item">'];
    html = [html, '\n', '       <a href="', folder, '/', filename, '.html" target="_self" class="btn btn-primary">Go</a>'];
    html = [html, '\n', '   </li>'];
    html = [html, '\n', '</ul>'];

end

end
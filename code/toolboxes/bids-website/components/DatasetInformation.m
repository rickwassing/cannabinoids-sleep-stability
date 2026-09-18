function html = DatasetInformation(props)

keys = fieldnames(props.info);
html = '';
html = [html, '\n', '<div class="list-group mb-3">'];
for i = 1:length(keys)
    html = [html, '\n', '<span class="list-group-item">'];
    html = [html, '\n', '<h6 class="mb-1">'];
    html = [html, '\n', keys{i}];
    html = [html, '\n', '</h6>'];
    html = [html, '\n', '<p class="mb-1">'];
    html = [html, '\n', props.info.(keys{i})];
    html = [html, '\n', '</p>'];
    html = [html, '\n', '</span>'];
end
html = [html, '\n', '</div>'];

end
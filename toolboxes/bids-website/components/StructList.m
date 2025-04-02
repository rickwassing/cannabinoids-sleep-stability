function html = StructList(props)

keys = fieldnames(props.info);
html = '';
html = [html, '\n', '<div class="list-group mb-3">'];
for i = 1:length(keys)
    val = props.info.(keys{i});
    if isstruct(val)
        p = struct();
        p.info = val;
        html = [html, '\n', StructList(p)];
        continue
    end
    if iscell(val)
        val = [val{:}];
    end
    if isnumeric(val)
        val = num2str(val);
    end
    html = [html, '\n', '<span class="list-group-item list-group-item-action">'];
    html = [html, '\n', '<h6 class="mb-1">'];
    html = [html, '\n', keys{i}];
    html = [html, '\n', '</h6>'];
    html = [html, '\n', '<p class="mb-1">'];
    html = [html, '\n', val];
    html = [html, '\n', '</p>'];
    html = [html, '\n', '</span>'];
end
html = [html, '\n', '</div>'];

end
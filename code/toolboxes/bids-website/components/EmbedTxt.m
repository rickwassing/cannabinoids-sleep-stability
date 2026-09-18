function html = EmbedTxt(props)

paragraph = readlines(props.file);
html = '';
html = [html, '\n', '<h2>', props.header, '</h2>'];
for i = 1:size(paragraph, 1)
    html = [html, '\n', '<p>', paragraph{i}, '</p>']; %#ok<AGROW>
end

end
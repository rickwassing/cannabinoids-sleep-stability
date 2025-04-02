function html = Breadcrumb(props)

html = '';
html = [html, '\n', '<nav aria-label="breadcrumb">'];
html = [html, '\n', '   <ol class="breadcrumb">'];
for i = 1:length(props.li)
    if i < length(props.li)
        html = [html, '\n', '       <li class="breadcrumb-item"><a href="', props.li(i).href, '">', props.li(i).label, '</a></li>']; %#ok<AGROW>
    else
        html = [html, '\n', '       <li class="breadcrumb-item active" aria-current="page">', props.li(i).label, '</li>']; %#ok<AGROW>
    end
end
html = [html, '\n', '   </ol>'];
html = [html, '\n', '</nav>'];

html = sprintf(html);

end
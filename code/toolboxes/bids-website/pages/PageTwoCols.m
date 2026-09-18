function html = PageTwoCols(props)
% =========================================================================
% Init
html = '';
% =========================================================================
% Start of the HTML code
p = struct();
p.title = props.title;
html = [html, '\n', Header(p)];
% =========================================================================
% Breadcrumbs
html = [html, '\n', Breadcrumb(props.breadcrumb)];
% =========================================================================
% Breadcrumb
html = [html, '\n', '           <div class="row">'];
html = [html, '\n', '               <div class="col-xs-12 col-sm-12 col-md-4 col-lg-3 col-xl-3">'];
html = [html, '\n', '                   <h1>', props.col(1).title, '</h1>'];
switch props.col(1).component
    case 'StructList'
        html = [html, '\n', StructList(props.col(1).props)];
end
html = [html, '\n', '               </div>'];
html = [html, '\n', '               <div class="mb-4 col-xs-12 col-sm-12 col-md-8 col-lg-9 col-xl-9">'];
html = [html, '\n', '                   <h1>', props.col(2).title, '</h1>'];
switch props.col(2).component
    case 'RawDataList'
        html = [html, '\n', RawDataList(props.col(2).props)];
    case 'ScansList'
        for i = 1:length(props.col(2).props)
            html = [html, '\n', ScansList(props.col(2).props(i))]; %#ok<AGROW>
        end
end
html = [html, '\n', '               </div>'];
html = [html, '\n', '           </div>'];
% =========================================================================
% End of the HTML code
html = [html, '\n', Footer([])];
% =========================================================================
% Write to file
try
html = sprintf(html);
catch ME
    keyboard
end
fid = fopen(props.filename, 'w');
fprintf(fid, '%s', html);
fclose(fid);

end
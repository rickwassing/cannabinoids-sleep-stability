function html = ImageList(props)

id = char(java.util.UUID.randomUUID);
html = '';
html = [html, '\n', '<div class="accordion" id="accordion-', id, '">'];
html = [html, '\n', '   <div class="accordion-item mb-3">'];
html = [html, '\n', '       <h2 class="accordion-header" id="heading-', id, '">'];
html = [html, '\n', '           <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#collapse-', id, '" aria-expanded="true" aria-controls="collapse-', id, '">'];
html = [html, '\n', props.title];
html = [html, '\n', '           </button>'];
html = [html, '\n', '       </h2>'];
html = [html, '\n', '       <div id="collapse-', id, '" class="accordion-collapse collapse" aria-labelledby="heading-', id, '" data-bs-parent="accordion-', id, '">'];
html = [html, '\n', '       <div class="accordion-body">'];
html = [html, '\n', '           <div class="container-fluid g-0 m-0">'];
html = [html, '\n', '               <div class="row">'];
for i = 1:length(props.imagefiles)
    html = [html, '\n', '<div class="', props.class, '">']; %#ok<*AGROW>
    if isfield(props, 'fullfilepath')
        if props.fullfilepath
            html = [html, '\n', '   <img src="', props.imagefiles(i).folder, filesep, props.imagefiles(i).name,'" class="img-fluid py-2 border-bottom">'];
        else
            html = [html, '\n', '   <img src="qc/', props.imagefiles(i).name,'" class="img-fluid py-2 border-bottom">'];
        end
    else
        html = [html, '\n', '   <img src="qc/', props.imagefiles(i).name,'" class="img-fluid py-2 border-bottom">'];
    end
    html = [html, '\n', '</div>'];
end
html = [html, '\n', '               </div>'];
html = [html, '\n', '           </div>'];
html = [html, '\n', '       </div>'];
html = [html, '\n', '   </div>'];
html = [html, '\n', '</div>'];




end
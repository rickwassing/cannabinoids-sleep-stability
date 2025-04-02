
clc

props = [];
props.title = 'Sigma ISF';
props.filename = './inspect-isf.html';

sub = dir('./derivatives/EEG-output-fstlvl/sub-*');
r = 0;
for i = 1:length(sub)
    ses = dir(sprintf('./derivatives/EEG-output-fstlvl/%s/ses-*', sub(i).name));
    for j = 1:length(ses)
        for k = {'a1cabssigma', 'a1cnormsigma'}
            r = r+1;
            p = [];
            p.class = 'bg-white';
            p.fullfilepath = true;
            p.title = sprintf('%s %s (%s)', sub(i).name, ses(j).name, k{:});
            p.imagefiles = dir(sprintf('./derivatives/EEG-output-fstlvl/%s/%s/images/sub*%s*.png', sub(i).name, ses(j).name, k{:}));
            props.row(r).component = 'ImageList';
            props.row(r).props = p;
        end
    end
end

html = PageOneCol(props);
